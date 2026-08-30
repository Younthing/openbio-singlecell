from __future__ import annotations

import time
from typing import Any

from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .cnv_analysis import (
    analyze_cnv_pca,
    analyze_cnv_score,
    analyze_infer_cnv,
    cnv_pca_code,
    cnv_score_code,
    infer_cnv_code,
)
from .expression_source import _CNV_SPEC
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_artifact_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_table_output,
)
from .staged_state_codec import CNV_STATE_CODEC, read_cnv_state, write_cnv_state
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

CNV_STATE_KIND = "OPENBIO_CNV_STATE"
CNV_EXPRESSION_SOURCE = _CNV_SPEC


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
    adata = read_anndata_input(inputs)
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
    return analysis_outputs(
        report,
        infer_cnv_code(**analysis_parameters),
        {
            "type": "artifact",
            "name": "cnv_state",
            "kind": CNV_STATE_KIND,
            "codec": CNV_STATE_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
    )


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
    return analysis_outputs(report, cnv_pca_code(**parameters), write_anndata_output(context, output))


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
    adata = read_anndata_input(inputs)
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
    return analysis_outputs(
        report,
        cnv_score_code(**parameters),
        write_anndata_output(context, output),
        write_table_output(context, table_result, kind="OPENBIO_SINGLE_CELL_TABLE"),
    )


__all__ = [
    "CNV_EXPRESSION_SOURCE",
    "CNV_STATE_KIND",
    "cnv_pca",
    "cnv_score",
    "infer_cnv",
]
