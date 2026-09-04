from __future__ import annotations

import time
from typing import Any

from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata, make_plot_result, make_summary_result, make_table_result
from .artifact_envelope import table_from_metadata
from .cnv_analysis import (
    analyze_cnv_pca,
    analyze_cnv_score,
    analyze_infer_cnv,
    cnv_pca_code,
    cnv_score_code,
    infer_cnv_code,
)
from .cnv_plotting import (
    cnv_heatmap_plot_code,
    cnv_score_plot_code,
    run_cnv_heatmap_plot,
    run_cnv_score_plot,
)
from .expression_source import _CNV_SPEC
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    read_table_input,
    require_artifact_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
    write_table_output,
)
from .staged_state_codec import CNV_STATE_CODEC, read_cnv_state, write_cnv_state
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

CNV_STATE_KIND = "OPENBIO_CNV_STATE"
CNV_EXPRESSION_SOURCE = _CNV_SPEC
CNV_REFERENCE = AnalysisReference(
    citation=(
        "Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic "
        "melanoma by single-cell RNA-seq. Science. 2016;352:189-196."
    ),
    doi="10.1126/science.aad0501",
    url="https://doi.org/10.1126/science.aad0501",
    kind="method",
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)


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
    evidence = summary["key_results"]
    table_parameters = {
        **{key: value for key, value in summary["parameters"].items() if value is not None},
        "input_cnv_state_fingerprint_sha256": evidence["input_cnv_state"]["artifact_fingerprint_sha256"],
        "input_cnv_matrix_fingerprint_sha256": evidence["input_cnv_state"]["cnv_matrix_fingerprint_sha256"],
        "partition_fingerprint_sha256": evidence["partition"]["fingerprint_sha256"],
        "output_fingerprint_sha256": evidence["output"]["output_fingerprint_sha256"],
        "table_fingerprint_sha256": evidence["table_fingerprint_sha256"],
    }
    table_result = make_table_result(
        title=f"CNV group scores by {groupby}",
        operation="cnv_score",
        parameters=table_parameters,
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


def cnv_heatmap_plot_owned(
    state: Any,
    *,
    groupby: str = "",
    view: Any = None,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_cnv_heatmap_plot(state, groupby=groupby, view=view)
    parameters = {"groupby": details["groupby"], "view": details["view_parameters"]}
    code = cnv_heatmap_plot_code(groupby=details["groupby"], view=details["view_parameters"])
    description = (
        f"Rendered {details['row_count']:,} {details['view'].replace('_', ' ')} row(s) across "
        f"{details['window_count']:,} stored genome-ordered inferred-CNV windows."
    )
    warnings = [
        "Expression-derived inferred-CNV signal is indirect and must be interpreted against the declared reference and genome assembly."
    ]
    plotted = make_plot_result(
        title=details["title"],
        operation="cnv_heatmap_plot",
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=details["input_cells"],
        input_genes=details["input_genes"],
        started_at=started_at,
        png=png,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCNVHeatmapPlot",
        title="CNV heatmap plot summary",
        operation="cnv_heatmap_plot",
        methods=(
            "Validated the exact typed Infer CNV producer, immutable artifact envelope, observation and feature axes, "
            "CNV matrix fingerprint, genome assembly, contiguous chromosome-window coordinates, categorical row "
            "partition, and finite displayed values before read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"title", "view_parameters"}},
        parameters=parameters,
        references=(CNV_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("anndata", "numpy", "pandas", "scipy", "infercnvpy", "matplotlib"),
        warnings=warnings,
        limitations=(
            "Expression-derived inferred-CNV signal is indirect and does not replace DNA copy-number measurement.",
            "Color represents the stored centered window signal, not absolute copy number, tumor probability, or a calibrated cutoff.",
            "Grouping and row order are descriptive; no Sample-level Condition test or causal trajectory was performed.",
        ),
        input_cells=details["input_cells"],
        input_genes=details["input_genes"],
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


def cnv_score_plot_owned(
    table: Any,
    *,
    view: Any = None,
) -> tuple[Any, Any, str]:
    from .contracts import TableResult

    if not isinstance(table, TableResult):
        raise TypeError("CNV Score Plot requires the table from CNV Group Score.")
    if not isinstance(table.source, dict) or table.parameters != table.source.get("parameters"):
        raise ValueError("CNV Score Plot requires internally consistent producer provenance.")
    producer = {
        "operation": table.source.get("operation"),
        "parameters": table.parameters,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "random_seed": table.random_seed,
    }
    started_at = time.perf_counter()
    png, details = run_cnv_score_plot(table.table, producer=producer, view=view)
    parameters = {"view": details["view_parameters"]}
    code = cnv_score_plot_code(producer=producer, view=details["view_parameters"])
    description = (
        f"Rendered {details['plotted_groups']:,} upstream-ranked groups using stored descriptive "
        f"{details['view'].replace('_', ' ')} CNV-score evidence."
    )
    warnings = list(table.warnings)
    plotted = make_plot_result(
        title=details["title"],
        operation="cnv_score_plot",
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        random_seed=table.random_seed,
        png=png,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCNVScorePlot",
        title="CNV score plot summary",
        operation="cnv_score_plot",
        methods=(
            "Validated the exact CNV Group Score producer, canonical group identity/rank axis, positive group "
            "counts, ordered absolute-value summaries, upstream state/matrix/partition/output fingerprints, and "
            "complete score-table fingerprint before preserving upstream order in read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"title", "view_parameters"}},
        parameters=parameters,
        references=(CNV_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("numpy", "pandas", "infercnvpy", "matplotlib"),
        warnings=warnings,
        limitations=(
            "The score is a descriptive group mean; cells within a group are not independent score observations.",
            "Scores are not absolute copy number, tumor probabilities, calibrated cutoffs, or hypothesis tests.",
            "The plot performs no Sample-level Condition inference and does not classify cells or groups.",
        ),
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        code=code,
        random_seed=table.random_seed,
    )
    return plotted, report, code


@register_operation("openbio.node.cnvheatmapplot")
def cnv_heatmap_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "CNV Heatmap Plot"
    require_input_names(inputs, {"cnv_state"}, operation=operation)
    require_parameters(parameters, {"groupby", "view"}, operation=operation)
    plotted, report, code = cnv_heatmap_plot_owned(_state_root(inputs), **parameters)
    return analysis_outputs(
        report,
        code,
        write_plot_output(context, plotted, kind="OPENBIO_SINGLE_CELL_PLOT"),
    )


@register_operation("openbio.node.cnvscoreplot")
def cnv_score_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "CNV Score Plot"
    require_input_names(inputs, {"table"}, operation=operation)
    require_parameters(parameters, {"view"}, operation=operation)
    table, metadata = read_table_input(inputs, "table", kind="OPENBIO_SINGLE_CELL_TABLE")
    plotted, report, code = cnv_score_plot_owned(table_from_metadata(metadata, table), **parameters)
    return analysis_outputs(
        report,
        code,
        write_plot_output(context, plotted, kind="OPENBIO_SINGLE_CELL_PLOT"),
    )


__all__ = [
    "CNV_EXPRESSION_SOURCE",
    "CNV_STATE_KIND",
    "cnv_heatmap_plot",
    "cnv_heatmap_plot_owned",
    "cnv_pca",
    "cnv_score",
    "cnv_score_plot",
    "cnv_score_plot_owned",
    "infer_cnv",
]
