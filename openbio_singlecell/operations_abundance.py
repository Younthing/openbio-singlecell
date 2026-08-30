from __future__ import annotations

import time
from typing import Any

from . import PLUGIN_VERSION
from .analysis_utils import make_summary_result, make_table_result
from .composition_modeling import (
    run_sccoda_differential_composition,
    run_tasccoda_differential_composition,
    sccoda_differential_composition_code,
    tasccoda_differential_composition_code,
)
from .composition_summary import run_sample_composition_summary, sample_composition_summary_code
from .milo_analysis import milo_differential_abundance_code, run_milo_differential_abundance
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_table_output,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"


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


def milo_owned(adata: Any, **parameters: Any) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    table, summary = run_milo_differential_abundance(adata, **parameters, openbio_version=PLUGIN_VERSION)
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
    return result, report, code


def sccoda_owned(adata: Any, **parameters: Any) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    table, summary = run_sccoda_differential_composition(adata, **parameters, openbio_version=PLUGIN_VERSION)
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
    return result, report, code


def tasccoda_owned(adata: Any, **parameters: Any) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    table, summary = run_tasccoda_differential_composition(adata, **parameters, openbio_version=PLUGIN_VERSION)
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
    return result, report, code


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
    return _records(context, *milo_owned(read_anndata_input(inputs), **parameters))


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
    return _records(context, *sccoda_owned(read_anndata_input(inputs), **parameters))


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
    return _records(context, *tasccoda_owned(read_anndata_input(inputs), **parameters))


__all__ = [
    "milo_differential_abundance",
    "milo_owned",
    "sample_composition_owned",
    "sample_composition_summary",
    "sccoda_differential_composition",
    "sccoda_owned",
    "tasccoda_differential_composition",
    "tasccoda_owned",
]
