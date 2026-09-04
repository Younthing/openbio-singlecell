from __future__ import annotations

import time
from typing import Any

from . import PLUGIN_VERSION
from .abundance_artifact_codecs import (
    COMPOSITION_MODEL_CODEC,
    COMPOSITION_MODEL_KIND,
    MILO_RESULT_CODEC,
    MILO_RESULT_KIND,
    read_composition_model_result,
    read_milo_result,
    write_composition_model_result,
    write_milo_result,
)
from .abundance_plotting import (
    milo_differential_abundance_plot_code,
    run_milo_differential_abundance_plot,
    run_sccoda_differential_composition_plot,
    run_tasccoda_differential_composition_plot,
    sccoda_differential_composition_plot_code,
    tasccoda_differential_composition_plot_code,
)
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import make_plot_result, make_summary_result, make_table_result
from .artifact_envelope import table_from_metadata
from .composition_modeling import (
    run_sccoda_differential_composition,
    run_tasccoda_differential_composition,
    sccoda_differential_composition_code,
    tasccoda_differential_composition_code,
)
from .composition_plot import run_sample_composition_plot, sample_composition_plot_code
from .composition_result import build_composition_model_result
from .composition_summary import run_sample_composition_summary, sample_composition_summary_code
from .milo_analysis import milo_differential_abundance_code, run_milo_differential_abundance
from .milo_result import build_milo_result
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    read_table_input,
    require_artifact_input,
    require_input_names,
    require_parameters,
    write_plot_output,
    write_table_output,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)
PANDAS_REFERENCE = AnalysisReference(
    citation="The pandas development team. pandas software documentation.",
    url="https://pandas.pydata.org/docs/",
    kind="software_documentation",
)


def _results(
    table: Any,
    summary: dict[str, Any],
    *,
    table_title: str,
    summary_title: str,
    operation: str,
    input_cells: int,
    input_genes: int,
    started_at: float,
    random_seed: int | None = None,
) -> tuple[Any, Any]:
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    common = {
        "operation": operation,
        "parameters": parameters,
        "description": summary["results"],
        "warnings": warnings,
        "input_cells": input_cells,
        "input_genes": input_genes,
        "started_at": started_at,
    }
    if random_seed is not None:
        common["random_seed"] = random_seed
    return (
        make_table_result(table=table, title=table_title, **common),
        make_summary_result(summary=summary, title=summary_title, **common),
    )


def _records(context: OperationContext, result: Any, report: Any, code: str) -> list[JSONValue]:
    return analysis_outputs(report, code, write_table_output(context, result, kind=TABLE_KIND))


def _milo_result_output(context: OperationContext, artifact: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("result")
    write_milo_result(root, artifact)
    return {
        "type": "artifact",
        "name": "result",
        "kind": MILO_RESULT_KIND,
        "codec": MILO_RESULT_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _composition_result_output(context: OperationContext, artifact: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("result")
    write_composition_model_result(root, artifact)
    return {
        "type": "artifact",
        "name": "result",
        "kind": COMPOSITION_MODEL_KIND,
        "codec": COMPOSITION_MODEL_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def sample_composition_owned(
    adata: Any,
    *,
    sample_key: str = "sample",
    condition_key: str = "condition",
    annotation_key: str = "cell_type",
    annotation_status: str = "unknown",
    max_output_rows: int = 2_000_000,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    table, summary = run_sample_composition_summary(
        adata,
        sample_key=sample_key,
        condition_key=condition_key,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        max_output_rows=max_output_rows,
        openbio_version=PLUGIN_VERSION,
    )
    result, report = _results(
        table,
        summary,
        table_title=f"Sample composition by {annotation_key}",
        summary_title="Sample composition summary",
        operation="sample_composition_summary",
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    code = sample_composition_summary_code(
        sample_key=sample_key,
        condition_key=condition_key,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        max_output_rows=max_output_rows,
        openbio_version=PLUGIN_VERSION,
    )
    return result, report, code


def sample_composition_plot_owned(
    table: Any, *, value: str = "proportion", input_genes: int = 0
) -> tuple[Any, Any, str]:
    if isinstance(input_genes, bool) or not isinstance(input_genes, int) or input_genes < 0:
        raise TypeError("Sample Composition Plot input_genes provenance must be a nonnegative integer.")
    started_at = time.perf_counter()
    png, details = run_sample_composition_plot(table, value=value)
    parameters = {"value": details["value"]}
    code = sample_composition_plot_code(value=details["value"])
    warnings = list(details["annotation_readability_warnings"])
    common = {
        "parameters": parameters,
        "description": (
            f"Rendered {details['plotted_samples']:,} Sample bars grouped by Condition with "
            f"{details['plotted_annotations']:,} stacked annotation categories."
        ),
        "warnings": warnings,
        "input_cells": details["input_cells"],
        "input_genes": input_genes,
        "started_at": started_at,
    }
    plotted = make_plot_result(
        png=png,
        title="Sample composition proportions" if details["value"] == "proportion" else "Sample composition counts",
        operation="sample_composition_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellSampleCompositionPlot",
        title="Sample composition plot summary",
        operation="sample_composition_plot",
        methods=(
            "Validated and rendered the canonical complete Sample-by-annotation table as one stacked bar per "
            "Sample. Samples were ordered into Condition groups in first-observed order; no Condition aggregation "
            "or inferential calculation was performed."
        ),
        results=common["description"],
        key_results=details,
        parameters=parameters,
        references=[PANDAS_REFERENCE, MATPLOTLIB_REFERENCE],
        software_packages=["pandas", "numpy", "matplotlib"],
        warnings=warnings,
        limitations=[
            "Displayed proportions describe captured-cell composition and do not estimate absolute tissue abundance.",
            "The plot is descriptive; Samples remain the independent biological units and no Condition contrast is tested.",
        ],
        input_cells=details["input_cells"],
        input_genes=input_genes,
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


def milo_owned(adata: Any, **parameters: Any) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    table, summary, evidence = run_milo_differential_abundance(
        adata, **parameters, openbio_version=PLUGIN_VERSION
    )
    artifact = build_milo_result(table=table, **evidence)
    result, report = _results(
        table,
        summary,
        table_title="Milo differential abundance",
        summary_title="Milo differential abundance summary",
        operation="milo_differential_abundance",
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=int(parameters["random_seed"]),
    )
    code = milo_differential_abundance_code(**parameters, openbio_version=PLUGIN_VERSION)
    return artifact, result, report, code


def milo_differential_abundance_plot_owned(
    result: Any,
    *,
    view: Any = None,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_milo_differential_abundance_plot(result, view=view)
    parameters = {"view": details["view_parameters"]}
    code = milo_differential_abundance_plot_code(view=details["view_parameters"])
    description = (
        f"Rendered Milo {details['view'].replace('_', ' ')} for {details['neighborhoods']:,} retained "
        f"neighborhoods in the {details['comparison_condition']} versus {details['reference_condition']} "
        "Sample-level differential-abundance contrast."
    )
    warnings = []
    if details["annotation_status"] == "provisional":
        warnings.append("Neighborhood annotations are provisional and require expert review.")
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": details["observations"],
        "input_genes": 0,
        "started_at": started_at,
        "random_seed": result.provenance["random_seed"],
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="milo_differential_abundance_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMiloDifferentialAbundancePlot",
        title="Milo differential-abundance plot summary",
        operation="milo_differential_abundance_plot",
        methods=(
            "Validated the complete typed Milo result, current-content fingerprints, observation and neighborhood "
            "axes, membership matrix, overlap graph, representative coordinates, and canonical result table before "
            "read-only rendering. Differential-abundance evidence remains neighborhood-level with Samples as the "
            "replicate unit."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=[MATPLOTLIB_REFERENCE],
        software_packages=["numpy", "pandas", "scipy", "matplotlib"],
        warnings=warnings,
        limitations=(
            "Spatial FDR is Milo neighborhood-overlap evidence and is not a cell-level p-value.",
            "Overlapping neighborhoods are not independent, and selected neighborhoods do not define discrete cell populations.",
            "Majority annotations inherit the declared provisional or curated annotation status.",
        ),
        input_cells=details["observations"],
        input_genes=0,
        started_at=started_at,
        code=code,
        random_seed=result.provenance["random_seed"],
    )
    return plotted, report, code


def sccoda_owned(adata: Any, **parameters: Any) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    table, summary, evidence = run_sccoda_differential_composition(
        adata, **parameters, openbio_version=PLUGIN_VERSION
    )
    artifact = build_composition_model_result(table=table, **evidence)
    result, report = _results(
        table,
        summary,
        table_title="scCODA differential composition",
        summary_title="scCODA differential composition summary",
        operation="sccoda_differential_composition",
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=int(parameters["random_seed"]),
    )
    code = sccoda_differential_composition_code(**parameters, openbio_version=PLUGIN_VERSION)
    return artifact, result, report, code


def _composition_differential_plot_owned(
    result: Any,
    *,
    expected_method: str,
    view: Any = None,
) -> tuple[Any, Any, str]:
    display_method = "scCODA" if expected_method == "sccoda" else "tascCODA"
    node_id = (
        "OpenBioSingleCellSccodaDifferentialCompositionPlot"
        if expected_method == "sccoda"
        else "OpenBioSingleCellTasccodaDifferentialCompositionPlot"
    )
    runner = (
        run_sccoda_differential_composition_plot
        if expected_method == "sccoda"
        else run_tasccoda_differential_composition_plot
    )
    code_builder = (
        sccoda_differential_composition_plot_code
        if expected_method == "sccoda"
        else tasccoda_differential_composition_plot_code
    )
    started_at = time.perf_counter()
    png, details = runner(result, view=view)
    parameters = {"view": details["view_parameters"]}
    code = code_builder(view=details["view_parameters"])
    description = (
        f"Rendered {display_method} {details['view'].replace('_', ' ')} from {details['posterior_draws']:,} "
        f"retained posterior draws for the {details['comparison_condition']} versus "
        f"{details['reference_condition']} relative compositional contrast; effects are defined against "
        f"reference cell type {details['reference_cell_type']!r}."
    )
    warnings = [
        "Only one posterior chain was retained; between-chain mixing and R-hat are unavailable.",
        "Pertpy 1.3.0 did not retain divergence counts for this audited composition runner.",
    ]
    if details["annotation_status"] == "provisional":
        warnings.append("Cell-type annotations are provisional and formal biological interpretation is invalid.")
    random_seed = result.model_metadata["model"].get("random_seed", 0)
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": 0,
        "input_genes": 0,
        "started_at": started_at,
        "random_seed": random_seed,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation=f"{expected_method}_differential_composition_plot",
        **common,
    )
    method_details = (
        "scCODA posterior expected-FDR selection over cell-type coefficients"
        if expected_method == "sccoda"
        else "tascCODA tree-adaptive selection over direct hierarchy-node effects with propagated leaf effects"
    )
    report, code = make_analysis_report(
        node_id=node_id,
        title=f"{display_method} differential-composition plot summary",
        operation=f"{expected_method}_differential_composition_plot",
        methods=(
            f"Validated the complete typed {display_method} result, exact producer/method identity, table, "
            "posterior and sampler current-content fingerprints, cell-type and hierarchy axes, and retained "
            f"model metadata before read-only rendering. The scientific contract is {method_details}; biological "
            "Samples, not cells, are the replicate unit."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=[MATPLOTLIB_REFERENCE],
        software_packages=["numpy", "pandas", "matplotlib"],
        warnings=warnings,
        limitations=(
            "Composition effects are relative to the declared cell-type reference and are not absolute abundance effects.",
            "The one-chain fit cannot establish between-chain mixing; R-hat and divergence counts are unavailable.",
            *(
                (
                    "Hierarchy-node credibility applies to direct modeled nodes; propagated leaf effects are not independent discoveries.",
                )
                if expected_method == "tasccoda"
                else ()
            ),
            "Observed Condition associations are not causal without an appropriate experimental or causal design.",
        ),
        input_cells=0,
        input_genes=0,
        started_at=started_at,
        code=code,
        random_seed=random_seed,
    )
    return plotted, report, code


def sccoda_differential_composition_plot_owned(
    result: Any,
    *,
    view: Any = None,
) -> tuple[Any, Any, str]:
    return _composition_differential_plot_owned(result, expected_method="sccoda", view=view)


def tasccoda_differential_composition_plot_owned(
    result: Any,
    *,
    view: Any = None,
) -> tuple[Any, Any, str]:
    return _composition_differential_plot_owned(result, expected_method="tasccoda", view=view)


def tasccoda_owned(adata: Any, **parameters: Any) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    table, summary, evidence = run_tasccoda_differential_composition(
        adata, **parameters, openbio_version=PLUGIN_VERSION
    )
    artifact = build_composition_model_result(table=table, **evidence)
    result, report = _results(
        table,
        summary,
        table_title="tascCODA differential composition",
        summary_title="tascCODA differential composition summary",
        operation="tasccoda_differential_composition",
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=int(parameters["random_seed"]),
    )
    code = tasccoda_differential_composition_code(**parameters, openbio_version=PLUGIN_VERSION)
    return artifact, result, report, code


@register_operation("openbio.node.samplecompositionsummary")
def sample_composition_summary(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Sample Composition Summary")
    require_parameters(
        parameters,
        {"sample_key", "condition_key", "annotation_key", "annotation_status", "max_output_rows"},
        operation="Sample Composition Summary",
    )
    return _records(context, *sample_composition_owned(read_anndata_input(inputs), **parameters))


@register_operation("openbio.node.samplecompositionplot")
def sample_composition_plot(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"table"}, operation="Sample Composition Plot")
    require_parameters(parameters, {"value"}, operation="Sample Composition Plot")
    table, metadata = read_table_input(inputs, "table", kind=TABLE_KIND)
    table_result = table_from_metadata(metadata, table)
    if (
        not isinstance(table_result.source, dict)
        or table_result.source.get("operation") != "sample_composition_summary"
    ):
        raise ValueError("Sample Composition Plot requires a Sample Composition Summary table artifact.")
    plotted, report, code = sample_composition_plot_owned(
        table_result.table, input_genes=table_result.input_genes, **parameters
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.milodifferentialabundance")
def milo_differential_abundance(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Milo Differential Abundance")
    require_parameters(
        parameters,
        {
            "sample_key",
            "condition_key",
            "reference_condition",
            "comparison_condition",
            "technical_batch_key",
            "categorical_covariate_keys_json",
            "continuous_covariate_keys_json",
            "annotation_key",
            "annotation_status",
            "representation_key",
            "n_neighbors",
            "neighborhood_proportion",
            "mixed_annotation_threshold",
            "spatial_fdr_threshold",
            "min_abs_log2_fold_change",
            "random_seed",
        },
        operation="Milo Differential Abundance",
    )
    artifact, result, report, code = milo_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(
        report,
        code,
        _milo_result_output(context, artifact),
        write_table_output(context, result, kind=TABLE_KIND),
    )


@register_operation("openbio.node.milodifferentialabundanceplot")
def milo_differential_abundance_plot(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"result"}, operation="Milo Differential Abundance Plot")
    require_parameters(parameters, {"view"}, operation="Milo Differential Abundance Plot")
    root = require_artifact_input(
        inputs,
        "result",
        kind=MILO_RESULT_KIND,
        codec=MILO_RESULT_CODEC,
    )
    plotted, report, code = milo_differential_abundance_plot_owned(read_milo_result(root), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.sccodadifferentialcomposition")
def sccoda_differential_composition(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="scCODA Differential Composition")
    require_parameters(
        parameters,
        {
            "sample_key",
            "annotation_key",
            "annotation_status",
            "condition_key",
            "reference_condition",
            "comparison_condition",
            "adjustment_covariate_keys_json",
            "reference_cell_type",
            "estimated_fdr",
            "num_samples",
            "num_warmup",
            "random_seed",
        },
        operation="scCODA Differential Composition",
    )
    artifact, result, report, code = sccoda_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(
        report,
        code,
        _composition_result_output(context, artifact),
        write_table_output(context, result, kind=TABLE_KIND),
    )


@register_operation("openbio.node.sccodadifferentialcompositionplot")
def sccoda_differential_composition_plot(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"result"}, operation="scCODA Differential Composition Plot")
    require_parameters(parameters, {"view"}, operation="scCODA Differential Composition Plot")
    root = require_artifact_input(
        inputs,
        "result",
        kind=COMPOSITION_MODEL_KIND,
        codec=COMPOSITION_MODEL_CODEC,
    )
    plotted, report, code = sccoda_differential_composition_plot_owned(
        read_composition_model_result(root),
        **parameters,
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.tasccodadifferentialcomposition")
def tasccoda_differential_composition(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="tascCODA Differential Composition")
    require_parameters(
        parameters,
        {
            "sample_key",
            "annotation_key",
            "annotation_status",
            "hierarchy_keys_json",
            "condition_key",
            "reference_condition",
            "comparison_condition",
            "adjustment_covariate_keys_json",
            "reference_cell_type",
            "aggregation_bias",
            "num_samples",
            "num_warmup",
            "random_seed",
        },
        operation="tascCODA Differential Composition",
    )
    artifact, result, report, code = tasccoda_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(
        report,
        code,
        _composition_result_output(context, artifact),
        write_table_output(context, result, kind=TABLE_KIND),
    )


@register_operation("openbio.node.tasccodadifferentialcompositionplot")
def tasccoda_differential_composition_plot(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"result"}, operation="tascCODA Differential Composition Plot")
    require_parameters(parameters, {"view"}, operation="tascCODA Differential Composition Plot")
    root = require_artifact_input(
        inputs,
        "result",
        kind=COMPOSITION_MODEL_KIND,
        codec=COMPOSITION_MODEL_CODEC,
    )
    plotted, report, code = tasccoda_differential_composition_plot_owned(
        read_composition_model_result(root),
        **parameters,
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


__all__ = [
    "milo_differential_abundance",
    "milo_differential_abundance_plot",
    "milo_differential_abundance_plot_owned",
    "milo_owned",
    "sample_composition_owned",
    "sample_composition_plot",
    "sample_composition_plot_owned",
    "sample_composition_summary",
    "sccoda_differential_composition",
    "sccoda_differential_composition_plot",
    "sccoda_differential_composition_plot_owned",
    "sccoda_owned",
    "tasccoda_differential_composition",
    "tasccoda_differential_composition_plot",
    "tasccoda_differential_composition_plot_owned",
    "tasccoda_owned",
]
