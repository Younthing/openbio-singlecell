from __future__ import annotations

import copy
import time
from typing import Any

from .analysis_utils import make_plot_result, make_summary_result, make_table_result
from .artifact_codecs import PLOT_CODEC, TABLE_CODEC, read_anndata, write_plot, write_table
from .artifact_envelope import result_metadata
from .augur import (
    AUGUR_ARTIFACT_TYPE,
    AUGUR_VIEWS,
    AugurResult,
    augur_code,
    augur_results_code,
    run_augur_artifact,
    select_augur_view,
)
from .augur_codec import AUGUR_CODEC, read_augur, write_augur
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .operations_input import (
    ANNDATA_CODEC,
    ANNDATA_KIND,
    require_artifact_input,
    require_input_names,
    require_parameters,
)
from .population_correlation import analyze_population_centroid_correlation, population_correlation_code
from .worker_protocol import JSONValue, OperationContext, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
AUGUR_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="Augur count source",
    default="X",
    include_raw=True,
    layer_default="counts",
)


def _anndata_input(inputs: dict[str, JSONValue]) -> Any:
    root = require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    return read_anndata(root)


def _table_output(context: OperationContext, result: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("table")
    write_table(root, result.table, result_metadata(result))
    return {
        "type": "artifact",
        "name": "table",
        "kind": TABLE_KIND,
        "codec": TABLE_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _plot_output(context: OperationContext, result: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("plot")
    write_plot(root, result.png, result_metadata(result))
    return {
        "type": "artifact",
        "name": "plot",
        "kind": PLOT_KIND,
        "codec": PLOT_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def augur_artifact_owned(
    adata: Any,
    *,
    sample_key: str = "sample",
    population_key: str = "cell_type",
    condition_key: str = "condition",
    control: str = "",
    treatment: str = "",
    classifier: str = "random_forest_classifier",
    source: DynamicExpressionSource | None = None,
    annotation_status: str = "unknown",
    technical_batch_key: str = "",
    n_subsamples: int = 50,
    subsample_size: int = 20,
    folds: int = 3,
    n_threads: int = 1,
    random_seed: int = 123,
    max_result_rows: int = 10_000_000,
    max_result_mib: float = 1024.0,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], Any, str]:
    started_at = time.perf_counter()
    expression = AUGUR_EXPRESSION_SOURCE.resolve(adata, source)
    tables, summary, metadata = run_augur_artifact(
        adata,
        sample_key=sample_key,
        population_key=population_key,
        condition_key=condition_key,
        control=control,
        treatment=treatment,
        classifier=classifier,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        annotation_status=annotation_status,
        technical_batch_key=technical_batch_key,
        n_subsamples=n_subsamples,
        subsample_size=subsample_size,
        folds=folds,
        n_threads=n_threads,
        random_seed=random_seed,
        max_result_rows=max_result_rows,
        max_result_mib=max_result_mib,
    )
    parameters = dict(summary["parameters"])
    report = make_summary_result(
        summary=summary,
        title="Augur exploratory population priorities",
        operation="augur",
        parameters=parameters,
        description=str(summary["results"]),
        warnings=[str(warning) for warning in summary["warnings"]],
        input_cells=int(adata.n_obs),
        input_genes=int(summary["key_results"]["selected_source_features"]),
        started_at=started_at,
        random_seed=random_seed,
    )
    return tables, summary, metadata, report, augur_code(**parameters)


def augur_legacy_owned(adata: Any, **parameters: Any) -> tuple[AugurResult, Any, str]:
    tables, summary, metadata, report, code = augur_artifact_owned(adata, **parameters)
    return AugurResult(tables=tables, summary=summary, metadata=metadata), report, code


def _augur_view_results(table: Any, summary: dict[str, Any], *, view: str, started_at: float) -> tuple[Any, Any, str]:
    selected = summary["selected_view"]
    key_results = summary["key_results"]
    parameters = {
        "view": view,
        "artifact_fingerprint_sha256": selected["artifact_fingerprint_sha256"],
        "table_fingerprint_sha256": selected["table_fingerprint_sha256"],
    }
    common = {
        "operation": "augur_results",
        "parameters": parameters,
        "warnings": [str(warning) for warning in summary["warnings"]],
        "input_cells": int(key_results["input_cells"]),
        "input_genes": int(key_results["selected_source_features"]),
        "started_at": started_at,
    }
    table_result = make_table_result(
        table=table,
        title=f"Augur {view.replace('_', ' ')}",
        description=f"Selected the canonical {view!r} view from a validated Augur result artifact.",
        **common,
    )
    report = make_summary_result(
        summary=summary,
        title=f"Augur {view.replace('_', ' ')} summary",
        description=f"Selected the canonical {view!r} view without rerunning or reinterpreting Augur.",
        **common,
    )
    return table_result, report, augur_results_code(view=view)


def augur_results_artifact_owned(
    tables: dict[str, Any],
    parent_summary: dict[str, Any],
    metadata: dict[str, Any],
    *,
    view: str = "priorities",
) -> tuple[Any, Any, str]:
    if view not in AUGUR_VIEWS:
        raise ValueError(f"Unsupported Augur result view: {view!r}; expected one of {list(AUGUR_VIEWS)}.")
    started_at = time.perf_counter()
    summary = copy.deepcopy(parent_summary)
    table = tables[view]
    summary["selected_view"] = {
        "node_id": "OpenBioSingleCellAugurResults",
        "view": view,
        "rows": int(table.shape[0]),
        "columns": int(table.shape[1]),
        "column_names": list(table.columns),
        "table_fingerprint_sha256": metadata["table_fingerprints_sha256"][view],
        "artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
        "operation": "pure_artifact_view",
    }
    return _augur_view_results(table, summary, view=view, started_at=started_at)


def augur_results_legacy_owned(result: AugurResult, *, view: str = "priorities") -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    table, summary = select_augur_view(result, view=view)
    return _augur_view_results(table, summary, view=view, started_at=started_at)


def population_correlation_owned(
    adata: Any,
    *,
    population_key: str = "cell_type",
    representation_key: str = "X_pca",
    correlation_method: str = "pearson",
    linkage_method: str = "complete",
    n_dimensions: int = 0,
    annotation_status: str = "unknown",
    color_map: str = "RdYlBu",
    show_numbers: bool = False,
    max_groups: int = 200,
    max_output_rows: int = 100000,
) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    table, png, summary = analyze_population_centroid_correlation(
        adata,
        population_key=population_key,
        representation_key=representation_key,
        correlation_method=correlation_method,
        linkage_method=linkage_method,
        n_dimensions=n_dimensions,
        annotation_status=annotation_status,
        color_map=color_map,
        show_numbers=show_numbers,
        max_groups=max_groups,
        max_output_rows=max_output_rows,
    )
    parameters = dict(summary["parameters"])
    description = str(summary["results"])
    warnings = [str(warning) for warning in summary["warnings"]]
    input_cells, input_genes = int(adata.n_obs), int(adata.n_vars)
    common = {
        "operation": "population_centroid_correlation",
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": input_cells,
        "input_genes": input_genes,
        "started_at": started_at,
    }
    table_result = make_table_result(
        table=table,
        title=f"Population centroid correlations by {population_key}",
        **common,
    )
    plot_result = make_plot_result(
        png=png,
        title=f"Population centroid correlation by {population_key}",
        **common,
    )
    report = make_summary_result(
        summary=summary,
        title="Population centroid correlation summary",
        **common,
    )
    code = population_correlation_code(
        population_key=population_key,
        representation_key=representation_key,
        correlation_method=correlation_method,
        linkage_method=linkage_method,
        n_dimensions=n_dimensions,
        annotation_status=annotation_status,
        color_map=color_map,
        show_numbers=show_numbers,
        max_groups=max_groups,
        max_output_rows=max_output_rows,
    )
    return table_result, plot_result, report, code


@register_operation("openbio.node.augur")
def augur(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Augur Cell Prioritization")
    require_parameters(
        parameters,
        {
            "sample_key",
            "population_key",
            "condition_key",
            "control",
            "treatment",
            "classifier",
            "source",
            "annotation_status",
            "technical_batch_key",
            "n_subsamples",
            "subsample_size",
            "folds",
            "n_threads",
            "random_seed",
            "max_result_rows",
            "max_result_mib",
        },
        operation="Augur Cell Prioritization",
    )
    tables, summary, metadata, report, code = augur_artifact_owned(_anndata_input(inputs), **parameters)
    root = context.create_output_directory("result")
    write_augur(root, tables, summary, metadata)
    return [
        {
            "type": "artifact",
            "name": "result",
            "kind": AUGUR_ARTIFACT_TYPE,
            "codec": AUGUR_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.augurresults")
def augur_results(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"result"}, operation="Augur Results")
    require_parameters(parameters, {"view"}, operation="Augur Results")
    root = require_artifact_input(inputs, "result", kind=AUGUR_ARTIFACT_TYPE, codec=AUGUR_CODEC)
    result, report, code = augur_results_artifact_owned(*read_augur(root), view=parameters["view"])
    return [
        _table_output(context, result),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.celltypecorrelation")
def population_centroid_correlation(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Population Centroid Correlation")
    require_parameters(
        parameters,
        {
            "population_key",
            "representation_key",
            "correlation_method",
            "linkage_method",
            "n_dimensions",
            "annotation_status",
            "color_map",
            "show_numbers",
            "max_groups",
            "max_output_rows",
        },
        operation="Population Centroid Correlation",
    )
    table, plot, report, code = population_correlation_owned(_anndata_input(inputs), **parameters)
    return [
        _table_output(context, table),
        _plot_output(context, plot),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


__all__ = [
    "AUGUR_EXPRESSION_SOURCE",
    "augur",
    "augur_artifact_owned",
    "augur_legacy_owned",
    "augur_results",
    "augur_results_artifact_owned",
    "augur_results_legacy_owned",
    "population_centroid_correlation",
    "population_correlation_owned",
]
