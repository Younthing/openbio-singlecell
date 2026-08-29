from __future__ import annotations

import time
from typing import Any

from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .artifact_codecs import TABLE_CODEC, read_anndata, write_anndata, write_table
from .artifact_envelope import result_metadata
from .cnv_analysis import (
    analyze_cnv_pca,
    analyze_cnv_score,
    analyze_infer_cnv,
    cnv_pca_code,
    cnv_score_code,
    infer_cnv_code,
)
from .expression_source import ExpressionSourceSpec
from .operations_input import (
    ANNDATA_CODEC,
    ANNDATA_KIND,
    require_artifact_input,
    require_input_names,
    require_parameters,
)
from .staged_state_codec import CNV_STATE_CODEC, read_cnv_state, write_cnv_state
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

CNV_STATE_KIND = "OPENBIO_CNV_STATE"
CNV_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="Full-gene normalized log-expression source for CNV inference",
    default="layer",
    include_raw=False,
    layer_default="log1p_norm",
)


def _write_anndata_output(context: OperationContext, adata: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("adata")
    write_anndata(root, adata)
    return {
        "type": "artifact",
        "name": "adata",
        "kind": ANNDATA_KIND,
        "codec": ANNDATA_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _anndata_root(inputs: dict[str, JSONValue], name: str) -> Any:
    root = require_artifact_input(inputs, name, kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    return read_anndata(root)


def _state_root(inputs: dict[str, JSONValue]) -> Any:
    root = require_artifact_input(
        inputs,
        "cnv_state",
        kind=CNV_STATE_KIND,
        codec=CNV_STATE_CODEC,
    )
    return read_cnv_state(root)


def _summary_result(
    summary: dict[str, Any],
    *,
    title: str,
    operation: str,
    started_at: float,
    cells: int,
    genes: int,
    random_seed: int = 0,
) -> Any:
    return make_summary_result(
        summary=summary,
        title=title,
        operation=operation,
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=random_seed,
    )


@register_operation("openbio.node.infercnv")
def infer_cnv(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Infer CNV"
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(
        parameters,
        {
            "source",
            "reference_key",
            "reference_categories",
            "sample_key",
            "genome_assembly",
            "window_size",
            "step",
            "lfc_clip",
            "dynamic_threshold",
            "exclude_chromosomes",
            "output_key",
            "minimum_reference_cells",
            "chunksize",
            "n_jobs",
            "max_output_gib",
            "overwrite_existing",
        },
        operation=operation,
    )
    adata = _anndata_root(inputs, "adata")
    source_value = parameters["source"]
    if not isinstance(source_value, dict):
        raise ProtocolError("Infer CNV source must be a DynamicCombo JSON object.")
    try:
        expression = CNV_EXPRESSION_SOURCE.resolve(adata, source_value)
    except (TypeError, ValueError) as error:
        raise ProtocolError(str(error)) from error
    analysis_parameters = {
        "source_kind": expression.kind,
        "layer_name": expression.layer_name,
        **{key: value for key, value in parameters.items() if key != "source"},
    }
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    state, summary = analyze_infer_cnv(adata, _owned=True, **analysis_parameters)
    report = _summary_result(
        summary,
        title="Expression-derived CNV inference summary",
        operation="infer_cnv",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    root = context.create_output_directory("cnv_state")
    write_cnv_state(root, state)
    return [
        {
            "type": "artifact",
            "name": "cnv_state",
            "kind": CNV_STATE_KIND,
            "codec": CNV_STATE_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": infer_cnv_code(**analysis_parameters)},
    ]


@register_operation("openbio.node.cnvpca")
def cnv_pca(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "CNV PCA"
    require_input_names(inputs, {"cnv_state"}, operation=operation)
    require_parameters(
        parameters,
        {"n_comps", "output_key", "overwrite_existing", "max_output_gib", "random_seed"},
        operation=operation,
    )
    state = _state_root(inputs)
    state_adata = state._owned_adata()
    cells, genes = int(state_adata.n_obs), int(state_adata.n_vars)
    started_at = time.perf_counter()
    output, summary = analyze_cnv_pca(state, _owned=True, **parameters)
    random_seed = parameters["random_seed"]
    if type(random_seed) is not int:
        raise ProtocolError("CNV PCA random_seed must be an integer.")
    finish_adata(
        output,
        "cnv_pca",
        dict(summary["parameters"]),
        cells,
        genes,
        started_at,
        random_seed=random_seed,
        warnings=[str(value) for value in summary["warnings"]],
    )
    report = _summary_result(
        summary,
        title="CNV PCA summary",
        operation="cnv_pca",
        started_at=started_at,
        cells=cells,
        genes=genes,
        random_seed=random_seed,
    )
    return [
        _write_anndata_output(context, output),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": cnv_pca_code(**parameters)},
    ]


@register_operation("openbio.node.cnvscore")
def cnv_score(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "CNV Score"
    require_input_names(inputs, {"cnv_state", "adata"}, operation=operation)
    require_parameters(
        parameters,
        {"groupby", "output_key", "overwrite_existing"},
        operation=operation,
    )
    state = _state_root(inputs)
    adata = _anndata_root(inputs, "adata")
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    output, table, summary = analyze_cnv_score(state, adata, _owned=True, **parameters)
    finish_adata(
        output,
        "cnv_score",
        dict(summary["parameters"]),
        cells,
        genes,
        started_at,
        warnings=[str(value) for value in summary["warnings"]],
    )
    groupby = parameters["groupby"]
    if not isinstance(groupby, str):
        raise ProtocolError("CNV Score groupby must be a string.")
    table_result = make_table_result(
        title=f"CNV group scores by {groupby}",
        operation="cnv_score",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        table=table,
    )
    report = _summary_result(
        summary,
        title="CNV group score summary",
        operation="cnv_score",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    table_root = context.create_output_directory("table")
    write_table(table_root, table, result_metadata(table_result))
    return [
        _write_anndata_output(context, output),
        {
            "type": "artifact",
            "name": "table",
            "kind": "OPENBIO_SINGLE_CELL_TABLE",
            "codec": TABLE_CODEC,
            "payload": table_root.relative_to(context.output_root).as_posix(),
        },
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": cnv_score_code(**parameters)},
    ]


__all__ = [
    "CNV_EXPRESSION_SOURCE",
    "CNV_STATE_KIND",
    "cnv_pca",
    "cnv_score",
    "infer_cnv",
]
