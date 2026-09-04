from __future__ import annotations

import importlib
import time
from collections.abc import Mapping
from typing import Any

from . import PLUGIN_VERSION, dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import make_plot_result, make_summary_result, make_table_result
from .artifact_envelope import table_from_metadata
from .differential_plotting import (
    pseudobulk_condition_contrast_plot_code,
    pseudobulk_qc_plot_code,
    run_pseudobulk_condition_contrast_plot,
    run_pseudobulk_qc_plot,
    run_scvi_population_de_evidence_plot,
    scvi_population_de_evidence_plot_code,
)
from .expression_source import _PSEUDOBULK_SPEC
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
from .pseudobulk import build_pseudobulk_summary, pseudobulk_code, run_pseudobulk
from .pseudobulk_artifact_codec import (
    PSEUDOBULK_CODEC,
    PSEUDOBULK_KIND,
    read_pseudobulk,
    write_pseudobulk,
)
from .sample_design import (
    build_pseudobulk_engine_summary,
    pseudobulk_engine_code,
    run_pseudobulk_edger,
    run_pseudobulk_pydeseq2,
)
from .scvi_de_evidence import (
    analyze_scvi_model_de_evidence,
    resolve_scvi_de_mode,
    resolve_scvi_population_scope,
    scvi_de_code,
)
from .scvi_model import SCVIModel
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
SCVI_MODEL_KIND = "OPENBIO_SCVI_MODEL"
SCVI_MODEL_CODEC = "scvi-native-directory"
PSEUDOBULK_SOURCE = _PSEUDOBULK_SPEC
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)


def _strict_source(value: JSONValue, adata: Any) -> Any:
    if not isinstance(value, dict):
        raise ProtocolError("Pseudobulk raw-count source must be a JSON object.")
    kind = value.get("source")
    expected = {"source", PSEUDOBULK_SOURCE.layer_input_id} if kind == "layer" else {"source"}
    if set(value) != expected:
        raise ProtocolError(
            f"Pseudobulk source fields must be exactly {sorted(expected)!r}; got {sorted(value)!r}."
        )
    try:
        return PSEUDOBULK_SOURCE.resolve(adata, value)
    except (TypeError, ValueError) as error:
        raise ProtocolError(str(error)) from error


def _strict_dynamic_combo(
    value: JSONValue,
    *,
    discriminator: str,
    variants: Mapping[str, set[str]],
    description: str,
) -> dict[str, JSONValue]:
    if not isinstance(value, dict):
        raise ProtocolError(f"{description} must be a DynamicCombo JSON object.")
    variant = value.get(discriminator)
    if not isinstance(variant, str) or variant not in variants:
        raise ProtocolError(f"{description} has an unsupported {discriminator}: {variant!r}.")
    expected = variants[variant]
    if set(value) != expected:
        raise ProtocolError(
            f"{description} fields for {variant!r} must be exactly {sorted(expected)!r}; "
            f"got {sorted(value)!r}."
        )
    return value


def run_pseudobulk_owned(
    adata: Any,
    *,
    sample_key: str = "sample",
    population_key: str = "cell_type",
    condition_key: str = "condition",
    technical_batch_key: str = "",
    categorical_covariate_keys: str = "",
    continuous_covariate_keys: str = "",
    source: Mapping[str, object] | None = None,
    inference_mode: str = "formal",
    min_cells: int = 10,
    min_counts: int = 1000,
) -> tuple[Any, Any, str]:
    expression = PSEUDOBULK_SOURCE.resolve(adata, source)
    dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells = int(adata.n_obs)
    artifact, diagnostics = run_pseudobulk(
        adata,
        sample_key=sample_key,
        population_key=population_key,
        condition_key=condition_key,
        technical_batch_key=technical_batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        inference_mode=inference_mode,
        min_cells=min_cells,
        min_counts=min_counts,
        openbio_version=PLUGIN_VERSION,
        decoupler_module=None,
    )
    code = pseudobulk_code(
        sample_key=sample_key,
        population_key=population_key,
        condition_key=condition_key,
        technical_batch_key=technical_batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        inference_mode=inference_mode,
        min_cells=min_cells,
        min_counts=min_counts,
        openbio_version=PLUGIN_VERSION,
    )
    summary = build_pseudobulk_summary(diagnostics)
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    report = make_summary_result(
        summary=summary,
        title="Sample-by-population pseudobulk counts",
        operation="pseudobulk",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=cells,
        input_genes=int(diagnostics["input_genes"]),
        started_at=started_at,
    )
    return artifact, report, code


def pseudobulk_qc_plot_owned(artifact: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_pseudobulk_qc_plot(artifact, view=view)
    parameters = {"view": details["view_parameters"]}
    code = pseudobulk_qc_plot_code(view=details["view_parameters"])
    description = (
        f"Rendered descriptive {details['view'].replace('_', ' ')} evidence for "
        f"{details['plotted_profiles']:,} retained Sample-by-population profiles across "
        f"{details['plotted_samples']:,} independent Sample declarations."
    )
    warnings: list[str] = []
    input_cells = int(artifact.metadata["input_dimensions"]["cells"])
    input_genes = int(details["plotted_genes"])
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": input_cells,
        "input_genes": input_genes,
        "started_at": started_at,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="pseudobulk_qc_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellPseudobulkQCPlot",
        title="Pseudobulk QC plot summary",
        operation="pseudobulk_qc_plot",
        methods=(
            "Validated the complete immutable Sample-by-population pseudobulk artifact and its current-content "
            "fingerprint before read-only Matplotlib rendering of retained profile cell counts, raw-count library "
            "sizes, and Sample-by-population coverage already stored by the producer."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=[MATPLOTLIB_REFERENCE],
        software_packages=["numpy", "pandas", "scipy", "matplotlib"],
        warnings=warnings,
        limitations=(
            "These are descriptive quality diagnostics; they do not test a Condition contrast.",
            "Sample identity, count-source identity, and population annotation remain analyst declarations.",
            "Cell-count and library-size differences can reflect recovery, Technical batch, composition, or biology.",
        ),
        input_cells=input_cells,
        input_genes=input_genes,
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


def _run_pseudobulk_engine_owned(
    artifact: Any,
    *,
    engine: str,
    population: str = "",
    reference_condition: str = "",
    comparison_condition: str = "",
    categorical_covariate_keys: str = "",
    continuous_covariate_keys: str = "",
    fdr_threshold: float = 0.05,
    min_abs_log2_fold_change: float = 0.0,
    min_count: int = 10,
    min_total_count: int = 15,
    large_n: int = 10,
    min_prop: float = 0.7,
    n_cpus: int = 1,
) -> tuple[Any, Any, str]:
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    runtime_parameters: dict[str, Any] = {
        "population": population,
        "reference_condition": reference_condition,
        "comparison_condition": comparison_condition,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "fdr_threshold": fdr_threshold,
        "min_abs_log2_fold_change": min_abs_log2_fold_change,
        "min_count": min_count,
        "min_total_count": min_total_count,
        "large_n": large_n,
        "min_prop": min_prop,
        "openbio_version": PLUGIN_VERSION,
    }
    if engine == "pydeseq2":
        runtime_parameters["n_cpus"] = n_cpus
        table, diagnostics = run_pseudobulk_pydeseq2(artifact, **runtime_parameters, pertpy_module=None)
        operation = "pseudobulk_deseq2"
        title = f"PyDESeq2: {comparison_condition.strip()} vs {reference_condition.strip()}"
        report_title = "Pseudobulk PyDESeq2 Condition contrast"
    elif engine == "edger":
        table, diagnostics = run_pseudobulk_edger(artifact, **runtime_parameters, pertpy_module=None)
        operation = "pseudobulk_edger"
        title = f"edgeR QL: {comparison_condition.strip()} vs {reference_condition.strip()}"
        report_title = "Pseudobulk edgeR QL Condition contrast"
    else:
        raise ValueError(f"Unsupported pseudobulk engine: {engine!r}")
    summary = build_pseudobulk_engine_summary(diagnostics)
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    metadata = artifact.metadata
    input_cells = int(metadata["input_dimensions"]["cells"])
    input_genes = int(diagnostics["expression_filter"]["genes_before"])
    result = make_table_result(
        title=title,
        operation=operation,
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=input_cells,
        input_genes=input_genes,
        started_at=started_at,
        table=science.pd.DataFrame(table),
    )
    code = pseudobulk_engine_code(engine=engine, **runtime_parameters)
    report = make_summary_result(
        summary=summary,
        title=report_title,
        operation=operation,
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=input_cells,
        input_genes=input_genes,
        started_at=started_at,
    )
    return result, report, code


def run_pseudobulk_edger_owned(artifact: Any, **parameters: Any) -> tuple[Any, Any, str]:
    return _run_pseudobulk_engine_owned(artifact, engine="edger", **parameters)


def run_pseudobulk_deseq2_owned(artifact: Any, **parameters: Any) -> tuple[Any, Any, str]:
    return _run_pseudobulk_engine_owned(artifact, engine="pydeseq2", **parameters)


def pseudobulk_condition_contrast_plot_owned(
    table: Any,
    *,
    view: Any = None,
) -> tuple[Any, Any, str]:
    from .contracts import TableResult, _metadata_value

    if not isinstance(table, TableResult):
        raise TypeError(
            "Pseudobulk Condition Contrast Plot requires the table result from Pseudobulk edgeR or PyDESeq2."
        )
    if not isinstance(table.source, dict) or _metadata_value(table.parameters) != table.source.get("parameters"):
        raise ValueError("Pseudobulk Condition Contrast Plot requires internally consistent producer provenance.")
    producer = {
        "operation": table.source.get("operation"),
        "parameters": table.parameters,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "random_seed": table.random_seed,
    }
    started_at = time.perf_counter()
    png, details = run_pseudobulk_condition_contrast_plot(table.table, producer=producer, view=view)
    parameters = {"view": details["view_parameters"]}
    code = pseudobulk_condition_contrast_plot_code(producer=producer, view=details["view_parameters"])
    description = (
        f"Rendered {details['engine']} {details['view'].replace('_', ' ')} evidence for "
        f"{details['tested_genes']:,} genes in the Sample-level "
        f"{details['comparison_condition']} versus {details['reference_condition']} Condition contrast "
        f"within population {details['population']!r}."
    )
    warnings = list(table.warnings)
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "started_at": started_at,
        "random_seed": table.random_seed,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="pseudobulk_condition_contrast_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellPseudobulkConditionContrastPlot",
        title="Pseudobulk Condition contrast plot summary",
        operation="pseudobulk_condition_contrast_plot",
        methods=(
            f"Validated the exact {details['engine']} producer identity, canonical full tested-gene table, "
            "producer-generated current-content fingerprint, Condition direction, thresholds, and adjusted-p-value "
            "universe before read-only rendering. The independent Sample is the replicate unit."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=[MATPLOTLIB_REFERENCE],
        software_packages=["numpy", "pandas", "matplotlib"],
        warnings=warnings,
        limitations=(
            "The plot visualizes the upstream frequentist test; it does not refit the model or change its tested universe.",
            "Reported log2 fold changes are unshrunk and can be unstable for low-count genes or few Samples.",
            "Biological Sample identity and count-source identity remain analyst declarations.",
        ),
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        code=code,
        random_seed=table.random_seed,
    )
    return plotted, report, code


def run_scvi_differential_owned(
    adata: Any,
    model: SCVIModel,
    *,
    groupby: str = "group",
    group1: str = "",
    group2: str = "",
    population_scope: Mapping[str, object] | None = None,
    mode: Mapping[str, object] | None = None,
    batch_handling: str = "shared_technical_batches",
    n_samples_overall: int = 5000,
    random_seed: int = 0,
) -> tuple[Any, Any, str]:
    subset_column, subset_value = resolve_scvi_population_scope(population_scope)
    mode_name, delta, fdr_target = resolve_scvi_de_mode(mode)
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    table, summary = analyze_scvi_model_de_evidence(
        adata,
        model,
        groupby=groupby,
        group1=group1,
        group2=group2,
        subset_column=subset_column,
        subset_value=subset_value,
        mode=mode_name,
        delta=delta,
        fdr_target=fdr_target,
        batch_handling=batch_handling,
        n_samples_overall=n_samples_overall,
        random_seed=random_seed,
    )
    parameters = dict(summary["parameters"])
    description = summary["results"]
    warnings = [str(warning) for warning in summary["warnings"]]
    input_cells = int(summary["comparison"]["selected_scope_cells"])
    result = make_table_result(
        title=f"scVI model evidence: {group1} vs {group2}",
        operation="scvi_model_de_evidence",
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=input_cells,
        input_genes=len(model.var_names),
        started_at=started_at,
        table=science.pd.DataFrame(table),
    )
    code = scvi_de_code(
        groupby=groupby,
        group1=group1,
        group2=group2,
        subset_column=subset_column,
        subset_value=subset_value,
        mode=mode_name,
        delta=delta,
        fdr_target=fdr_target,
        batch_handling=batch_handling,
        n_samples_overall=n_samples_overall,
        random_seed=random_seed,
    )
    report = make_summary_result(
        summary=summary,
        title="scVI model DE evidence summary",
        operation="scvi_model_de_evidence",
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=input_cells,
        input_genes=len(model.var_names),
        started_at=started_at,
        random_seed=random_seed,
    )
    return result, report, code


def scvi_population_de_evidence_plot_owned(
    table: Any,
    *,
    view: Any = None,
) -> tuple[Any, Any, str]:
    from .contracts import TableResult, _metadata_value

    if not isinstance(table, TableResult):
        raise TypeError("scVI Population DE Evidence Plot requires the table from scVI Model DE Evidence.")
    if not isinstance(table.source, dict) or _metadata_value(table.parameters) != table.source.get("parameters"):
        raise ValueError("scVI Population DE Evidence Plot requires internally consistent producer provenance.")
    producer = {
        "operation": table.source.get("operation"),
        "parameters": table.parameters,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "random_seed": table.random_seed,
    }
    started_at = time.perf_counter()
    png, details = run_scvi_population_de_evidence_plot(table.table, producer=producer, view=view)
    parameters = {"view": details["view_parameters"]}
    code = scvi_population_de_evidence_plot_code(producer=producer, view=details["view_parameters"])
    description = (
        f"Rendered scVI {details['view'].replace('_', ' ')} for {details['tested_fitted_features']:,} "
        f"fitted features comparing {details['group1']!r} with {details['group2']!r}; this is "
        "cell/model-conditional posterior evidence, not Sample-level Condition inference."
    )
    warnings = list(table.warnings)
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "started_at": started_at,
        "random_seed": table.random_seed,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="scvi_population_de_evidence_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellSCVIPopulationDEEvidencePlot",
        title="scVI population DE evidence plot summary",
        operation="scvi_population_de_evidence_plot",
        methods=(
            "Validated the exact scVI Model DE Evidence producer, canonical mode-specific fitted-feature table, "
            "producer-generated current-content fingerprint, population direction, posterior probabilities, "
            "decoded-expression quantities, and retained upstream ranking before read-only rendering. The "
            "statistical unit is the selected cells and fitted-model posterior draws, not independent Samples."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=[MATPLOTLIB_REFERENCE],
        software_packages=["numpy", "pandas", "matplotlib"],
        warnings=warnings,
        limitations=(
            "No between-Sample variance is estimated; use pseudobulk edgeR or PyDESeq2 for Condition inference.",
            "Posterior expected-FDR tags are not p-values or Benjamini–Hochberg adjusted p-values.",
            "Evidence is conditional on the fitted model, fitted feature view, selected cells, and nuisance declarations.",
        ),
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        code=code,
        random_seed=table.random_seed,
    )
    return plotted, report, code


def _table_record(context: OperationContext, result: Any) -> dict[str, JSONValue]:
    table = result.table
    science = dependencies.require_scientific_dependencies()
    nullable_float_columns = [
        name
        for name in table.columns
        if science.pd.api.types.is_float_dtype(table[name].dtype) and bool(table[name].isna().any())
    ]
    if nullable_float_columns:
        # Missing statistical estimates are scientific nulls, not non-finite JSON numbers. Only the
        # affected columns are materialized into pandas' portable nullable dtype for JSONL encoding.
        table = table.copy(deep=False)
        for name in nullable_float_columns:
            table[name] = table[name].astype(science.pd.Float64Dtype())
    return write_table_output(context, result, kind=TABLE_KIND, table=table)


@register_operation("openbio.node.pseudobulk")
def pseudobulk(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    expected = {
        "sample_key",
        "population_key",
        "condition_key",
        "technical_batch_key",
        "categorical_covariate_keys",
        "continuous_covariate_keys",
        "source",
        "inference_mode",
        "min_cells",
        "min_counts",
    }
    require_input_names(inputs, {"adata"}, operation="Pseudobulk")
    require_parameters(parameters, expected, operation="Pseudobulk")
    adata = read_anndata_input(inputs)
    expression = _strict_source(parameters["source"], adata)
    owned_parameters = dict(parameters)
    owned_parameters["source"] = expression.parameters()
    artifact, report, code = run_pseudobulk_owned(adata, **owned_parameters)
    root = context.create_output_directory("pseudobulk")
    write_pseudobulk(root, artifact)
    return analysis_outputs(
        report,
        code,
        {
            "type": "artifact",
            "name": "pseudobulk",
            "kind": PSEUDOBULK_KIND,
            "codec": PSEUDOBULK_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
    )


@register_operation("openbio.node.pseudobulkqcplot")
def pseudobulk_qc_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"pseudobulk"}, operation="Pseudobulk QC Plot")
    require_parameters(parameters, {"view"}, operation="Pseudobulk QC Plot")
    root = require_artifact_input(
        inputs,
        "pseudobulk",
        kind=PSEUDOBULK_KIND,
        codec=PSEUDOBULK_CODEC,
    )
    plotted, report, code = pseudobulk_qc_plot_owned(read_pseudobulk(root), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


_ENGINE_PARAMETERS = {
    "population",
    "reference_condition",
    "comparison_condition",
    "categorical_covariate_keys",
    "continuous_covariate_keys",
    "fdr_threshold",
    "min_abs_log2_fold_change",
    "min_count",
    "min_total_count",
    "large_n",
    "min_prop",
}


def _pseudobulk_engine_input(
    inputs: dict[str, JSONValue], parameters: dict[str, JSONValue], *, operation: str, deseq2: bool
) -> Any:
    require_input_names(inputs, {"pseudobulk"}, operation=operation)
    require_parameters(parameters, _ENGINE_PARAMETERS | ({"n_cpus"} if deseq2 else set()), operation=operation)
    root = require_artifact_input(
        inputs,
        "pseudobulk",
        kind=PSEUDOBULK_KIND,
        codec=PSEUDOBULK_CODEC,
    )
    return read_pseudobulk(root)


@register_operation("openbio.node.pseudobulkedger")
def pseudobulk_edger(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    artifact = _pseudobulk_engine_input(inputs, parameters, operation="Pseudobulk edgeR", deseq2=False)
    result, report, code = run_pseudobulk_edger_owned(artifact, **parameters)
    return analysis_outputs(report, code, _table_record(context, result))


@register_operation("openbio.node.pseudobulkdeseq2")
def pseudobulk_deseq2(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    artifact = _pseudobulk_engine_input(inputs, parameters, operation="Pseudobulk DESeq2", deseq2=True)
    result, report, code = run_pseudobulk_deseq2_owned(artifact, **parameters)
    return analysis_outputs(report, code, _table_record(context, result))


@register_operation("openbio.node.pseudobulkconditioncontrastplot")
def pseudobulk_condition_contrast_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"table"}, operation="Pseudobulk Condition Contrast Plot")
    require_parameters(parameters, {"view"}, operation="Pseudobulk Condition Contrast Plot")
    table, metadata = read_table_input(inputs, "table", kind=TABLE_KIND)
    plotted, report, code = pseudobulk_condition_contrast_plot_owned(
        table_from_metadata(metadata, table),
        **parameters,
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


def _scvi_training_parameters(adata: Any) -> dict[str, Any]:
    metadata = adata.uns.get("openbio_singlecell")
    history = metadata.get("analysis_history") if isinstance(metadata, Mapping) else None
    entries = list(history.values()) if isinstance(history, Mapping) else []
    for entry in reversed(entries):
        if isinstance(entry, Mapping) and entry.get("operation") == "scvi_integration":
            parameters = entry.get("parameters")
            if isinstance(parameters, Mapping):
                return dict(parameters)
    raise ValueError("Native scVI artifact is missing its OpenBio training provenance.")


def _load_scvi_model(root: Any, adata: Any) -> SCVIModel:
    scvi = importlib.import_module("scvi")
    model_namespace = getattr(scvi, "model", None)
    model_class = getattr(model_namespace, "SCVI", None)
    load = getattr(model_class, "load", None)
    if not callable(load):
        raise RuntimeError("Installed scvi-tools does not expose the required SCVI.load API.")
    model = load(str(root), adata=adata)
    if getattr(model, "adata", None) is not adata:
        raise RuntimeError("Loaded scVI model did not attach to the worker-owned AnnData input.")
    return SCVIModel(model, adata, _scvi_training_parameters(adata))


@register_operation("openbio.node.scvidifferentialexpression")
def scvi_differential_expression(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    expected = {
        "groupby",
        "group1",
        "group2",
        "population_scope",
        "mode",
        "batch_handling",
        "n_samples_overall",
        "random_seed",
    }
    require_input_names(inputs, {"adata", "model"}, operation="scVI Model DE Evidence")
    require_parameters(parameters, expected, operation="scVI Model DE Evidence")
    owned_parameters = dict(parameters)
    owned_parameters["population_scope"] = _strict_dynamic_combo(
        parameters["population_scope"],
        discriminator="population_scope",
        variants={
            "all": {"population_scope"},
            "obs_value": {"population_scope", "subset_column", "subset_value"},
        },
        description="scVI population_scope",
    )
    owned_parameters["mode"] = _strict_dynamic_combo(
        parameters["mode"],
        discriminator="mode",
        variants={
            "change": {"mode", "delta", "fdr_target"},
            "vanilla": {"mode"},
        },
        description="scVI mode",
    )
    adata = read_anndata_input(inputs)
    model_root = require_artifact_input(inputs, "model", kind=SCVI_MODEL_KIND, codec=SCVI_MODEL_CODEC)
    model = _load_scvi_model(model_root, adata)
    result, report, code = run_scvi_differential_owned(adata, model, **owned_parameters)
    return analysis_outputs(report, code, _table_record(context, result))


@register_operation("openbio.node.scvipopulationdeevidenceplot")
def scvi_population_de_evidence_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"table"}, operation="scVI Population DE Evidence Plot")
    require_parameters(parameters, {"view"}, operation="scVI Population DE Evidence Plot")
    table, metadata = read_table_input(inputs, "table", kind=TABLE_KIND)
    plotted, report, code = scvi_population_de_evidence_plot_owned(
        table_from_metadata(metadata, table),
        **parameters,
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


__all__ = [
    "PSEUDOBULK_CODEC",
    "PSEUDOBULK_KIND",
    "pseudobulk",
    "pseudobulk_condition_contrast_plot",
    "pseudobulk_condition_contrast_plot_owned",
    "pseudobulk_deseq2",
    "pseudobulk_edger",
    "pseudobulk_qc_plot",
    "pseudobulk_qc_plot_owned",
    "run_pseudobulk_deseq2_owned",
    "run_pseudobulk_edger_owned",
    "run_pseudobulk_owned",
    "run_scvi_differential_owned",
    "scvi_differential_expression",
    "scvi_population_de_evidence_plot",
    "scvi_population_de_evidence_plot_owned",
]
