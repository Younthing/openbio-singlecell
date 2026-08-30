from __future__ import annotations

import importlib
import time
from collections.abc import Mapping
from typing import Any

from . import PLUGIN_VERSION, dependencies
from .analysis_utils import make_summary_result, make_table_result
from .artifact_codecs import TABLE_CODEC, read_anndata, write_table
from .artifact_envelope import result_metadata
from .expression_source import ExpressionSourceSpec
from .operations_input import (
    ANNDATA_CODEC,
    ANNDATA_KIND,
    require_artifact_input,
    require_input_names,
    require_parameters,
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

PSEUDOBULK_SOURCE = ExpressionSourceSpec(
    description="Pseudobulk raw-count source",
    default="layer",
    include_raw=True,
    layer_default="counts",
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


def _table_record(context: OperationContext, result: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("table")
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
    write_table(root, table, result_metadata(result))
    return {
        "type": "artifact",
        "name": "table",
        "kind": TABLE_KIND,
        "codec": TABLE_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _small_records(report: Any, code: str) -> list[JSONValue]:
    return [
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


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
    adata = read_anndata(require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC))
    expression = _strict_source(parameters["source"], adata)
    owned_parameters = dict(parameters)
    owned_parameters["source"] = expression.parameters()
    artifact, report, code = run_pseudobulk_owned(adata, **owned_parameters)
    root = context.create_output_directory("pseudobulk")
    write_pseudobulk(root, artifact)
    return [
        {
            "type": "artifact",
            "name": "pseudobulk",
            "kind": PSEUDOBULK_KIND,
            "codec": PSEUDOBULK_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
        *_small_records(report, code),
    ]


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
    return [_table_record(context, result), *_small_records(report, code)]


@register_operation("openbio.node.pseudobulkdeseq2")
def pseudobulk_deseq2(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    artifact = _pseudobulk_engine_input(inputs, parameters, operation="Pseudobulk DESeq2", deseq2=True)
    result, report, code = run_pseudobulk_deseq2_owned(artifact, **parameters)
    return [_table_record(context, result), *_small_records(report, code)]


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
    adata = read_anndata(require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC))
    model_root = require_artifact_input(inputs, "model", kind=SCVI_MODEL_KIND, codec=SCVI_MODEL_CODEC)
    model = _load_scvi_model(model_root, adata)
    result, report, code = run_scvi_differential_owned(adata, model, **owned_parameters)
    return [_table_record(context, result), *_small_records(report, code)]


__all__ = [
    "PSEUDOBULK_CODEC",
    "PSEUDOBULK_KIND",
    "pseudobulk",
    "pseudobulk_deseq2",
    "pseudobulk_edger",
    "run_pseudobulk_deseq2_owned",
    "run_pseudobulk_edger_owned",
    "run_pseudobulk_owned",
    "run_scvi_differential_owned",
    "scvi_differential_expression",
]
