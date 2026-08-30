from __future__ import annotations

import copy
import inspect
import json
import math
from collections.abc import Mapping, Sequence
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .analysis_reporting import _package_version, collect_software_versions

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


SCVI_DE_SUMMARY_SCHEMA = "openbio-singlecell/scvi-model-de-evidence/v1"
SCVI_DE_NODE_ID = "OpenBioSingleCellSCVIDifferentialExpression"

SCVI_DE_REFERENCES = tuple(
    MappingProxyType(reference)
    for reference in (
        {
            "citation": (
                "Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell "
                "transcriptomics. Nature Methods. 2018;15:1053-1058."
            ),
            "doi": "10.1038/s41592-018-0229-2",
            "url": "https://doi.org/10.1038/s41592-018-0229-2",
            "kind": "method",
        },
        {
            "citation": (
                "Boyeau P, Lopez R, Regier J, Gayoso A, Jordan MI, Yosef N. Deep generative models for detecting "
                "differential expression in single cells. bioRxiv. 2019."
            ),
            "doi": "10.1101/794289",
            "url": "https://doi.org/10.1101/794289",
            "kind": "method",
        },
        {
            "citation": (
                "Gayoso A, Lopez R, Xing G, et al. scvi-tools: a library for deep probabilistic analysis of "
                "single-cell omics data. Nature Biotechnology. 2022;40:163-166."
            ),
            "doi": "10.1038/s41587-021-01206-w",
            "url": "https://doi.org/10.1038/s41587-021-01206-w",
            "kind": "software",
        },
        {
            "citation": "scvi-tools 1.5.0 differential expression public API and user guide.",
            "doi": None,
            "url": "https://docs.scvi-tools.org/en/1.5.0/api/reference/scvi.model.SCVI.html#scvi.model.SCVI.differential_expression",
            "kind": "software_documentation",
        },
        {
            "citation": (
                "Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential "
                "expression. Nature Communications. 2021;12:5692."
            ),
            "doi": "10.1038/s41467-021-25960-2",
            "url": "https://doi.org/10.1038/s41467-021-25960-2",
            "kind": "practice",
        },
        {
            "citation": (
                "Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias in "
                "single-cell studies. Nature Communications. 2021;12:738."
            ),
            "doi": "10.1038/s41467-021-21038-1",
            "url": "https://doi.org/10.1038/s41467-021-21038-1",
            "kind": "practice",
        },
    )
)

_COMMON_NUMERIC_COLUMNS = (
    "scale1",
    "scale2",
    "raw_mean1",
    "raw_mean2",
    "non_zeros_proportion1",
    "non_zeros_proportion2",
    "raw_normalized_mean1",
    "raw_normalized_mean2",
)
_CHANGE_COLUMNS = (
    "proba_de",
    "proba_not_de",
    "bayes_factor",
    "scale1",
    "scale2",
    "pseudocounts",
    "delta",
    "lfc_mean",
    "lfc_median",
    "lfc_std",
    "lfc_min",
    "lfc_max",
    "raw_mean1",
    "raw_mean2",
    "non_zeros_proportion1",
    "non_zeros_proportion2",
    "raw_normalized_mean1",
    "raw_normalized_mean2",
    "is_de_fdr",
)
_VANILLA_COLUMNS = (
    "proba_m1",
    "proba_m2",
    "bayes_factor",
    *_COMMON_NUMERIC_COLUMNS,
)


def resolve_scvi_population_scope(value: Mapping[str, object] | None) -> tuple[str | None, str | None]:
    if value is None:
        return None, None
    if not isinstance(value, Mapping):
        raise TypeError("scVI population_scope must be a DynamicCombo value.")
    variant = value.get("population_scope")
    if variant == "all":
        unknown = sorted(str(key) for key in value if key != "population_scope")
        if unknown:
            raise ValueError(f"scVI all-population scope received inactive or unknown parameters: {unknown!r}.")
        return None, None
    if variant != "obs_value":
        raise ValueError(f"Unsupported scVI population scope: {variant!r}.")
    unknown = sorted(str(key) for key in value if key not in {"population_scope", "subset_column", "subset_value"})
    if unknown:
        raise ValueError(f"scVI obs_value scope received unknown parameters: {unknown!r}.")
    subset_column = value.get("subset_column")
    subset_value = value.get("subset_value")
    if not isinstance(subset_column, str) or not subset_column or subset_column != subset_column.strip():
        raise ValueError("scVI obs_value subset_column must be a canonical nonempty string.")
    if not isinstance(subset_value, str) or not subset_value or subset_value != subset_value.strip():
        raise ValueError("scVI obs_value subset_value must be a canonical nonempty string.")
    return subset_column, subset_value


def resolve_scvi_de_mode(value: Mapping[str, object] | None) -> tuple[str, float | None, float | None]:
    if value is None:
        return "change", 0.25, 0.05
    if not isinstance(value, Mapping):
        raise TypeError("scVI mode must be a DynamicCombo value.")
    variant = value.get("mode")
    if variant == "vanilla":
        unknown = sorted(str(key) for key in value if key != "mode")
        if unknown:
            raise ValueError(f"scVI vanilla mode received inactive or unknown parameters: {unknown!r}.")
        return "vanilla", None, None
    if variant != "change":
        raise ValueError(f"Unsupported scVI model-DE mode: {variant!r}.")
    unknown = sorted(str(key) for key in value if key not in {"mode", "delta", "fdr_target"})
    if unknown:
        raise ValueError(f"scVI change mode received unknown parameters: {unknown!r}.")
    delta = value.get("delta", 0.25)
    fdr_target = value.get("fdr_target", 0.05)
    if isinstance(delta, bool) or not isinstance(delta, (int, float)) or not math.isfinite(float(delta)):
        raise TypeError("scVI change-mode delta must be a finite positive number.")
    if float(delta) <= 0:
        raise ValueError("scVI change-mode delta must be positive.")
    if isinstance(fdr_target, bool) or not isinstance(fdr_target, (int, float)) or not math.isfinite(float(fdr_target)):
        raise TypeError("scVI change-mode fdr_target must be a finite number between zero and one.")
    if not 0 < float(fdr_target) < 1:
        raise ValueError("scVI change-mode fdr_target must be strictly between zero and one.")
    return "change", float(delta), float(fdr_target)


def _finite_numeric_column(table: DataFrame, column: str) -> Any:
    import numpy as np
    import pandas as pd

    series = table[column]
    if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
        raise TypeError(f"scVI backend column {column!r} must have a non-Boolean numeric dtype.")
    values = np.asarray(series, dtype=float)
    if values.ndim != 1 or values.size != len(table) or not bool(np.isfinite(values).all()):
        raise RuntimeError(f"scVI backend column {column!r} must contain only finite numeric values.")
    return values


def _validate_exact_text_column(table: DataFrame, column: str, expected: str) -> None:
    values = table[column]
    if bool(values.isna().any()) or not all(isinstance(value, str) for value in values):
        raise RuntimeError(f"scVI backend column {column!r} must contain exact non-missing strings.")
    unexpected = sorted({value for value in values if value != expected})
    if unexpected:
        raise RuntimeError(
            f"scVI backend returned unexpected {column!r} labels {unexpected!r}; expected only {expected!r}."
        )


def canonicalize_scvi_de_table(
    raw_table: Any,
    *,
    fitted_features: Sequence[str],
    group1: str,
    group2: str,
    mode: str,
    delta: float | None,
    fdr_target: float | None,
) -> DataFrame:
    """Validate the released scvi-tools output and return the stable OpenBio schema."""
    import numpy as np
    import pandas as pd

    if not isinstance(raw_table, pd.DataFrame):
        raise TypeError("scVI differential_expression must return a pandas DataFrame.")
    table = raw_table.copy(deep=True)
    if table.empty:
        raise RuntimeError("scVI differential_expression returned an empty table.")
    if "gene" in table.columns:
        raise RuntimeError("scVI backend output contains a colliding 'gene' column.")
    genes = list(table.index)
    if not all(isinstance(gene, str) and gene and gene == gene.strip() for gene in genes):
        raise RuntimeError("scVI backend feature index must contain canonical nonempty string identifiers.")
    if len(genes) != len(set(genes)):
        raise RuntimeError("scVI backend feature index contains duplicates.")
    fitted = list(fitted_features)
    if len(genes) != len(fitted) or set(genes) != set(fitted):
        missing = [gene for gene in fitted if gene not in set(genes)][:5]
        unknown = [gene for gene in genes if gene not in set(fitted)][:5]
        raise RuntimeError(
            "scVI backend output must contain exactly one row per fitted feature; "
            f"missing={missing!r}, unknown={unknown!r}."
        )

    dynamic_fdr_columns = [str(column) for column in table.columns if str(column).startswith("is_de_fdr_")]
    if mode == "change":
        if len(dynamic_fdr_columns) != 1:
            raise RuntimeError("scVI change-mode output must contain exactly one dynamic is_de_fdr_<target> column.")
        if "is_de_fdr" in table.columns:
            raise RuntimeError("scVI backend output contains a colliding canonical is_de_fdr column.")
        try:
            returned_target = float(dynamic_fdr_columns[0].removeprefix("is_de_fdr_"))
        except ValueError as error:
            raise RuntimeError("scVI backend returned an unparsable posterior-FDR column name.") from error
        if fdr_target is None or not math.isclose(returned_target, fdr_target, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(
                f"scVI backend posterior-FDR target {returned_target!r} does not match requested {fdr_target!r}."
            )
        table = table.rename(columns={dynamic_fdr_columns[0]: "is_de_fdr"})
        scientific_columns = _CHANGE_COLUMNS
    else:
        if dynamic_fdr_columns or "is_de_fdr" in table.columns:
            raise RuntimeError("scVI vanilla output must not contain a posterior-FDR tag.")
        scientific_columns = _VANILLA_COLUMNS

    expected_columns = {"comparison", "group1", "group2", *scientific_columns}
    missing_columns = sorted(expected_columns - set(table.columns))
    unexpected_columns = sorted(set(table.columns) - expected_columns)
    if missing_columns or unexpected_columns:
        raise RuntimeError(
            "scVI backend output does not match the pinned stable schema; "
            f"missing={missing_columns!r}, unexpected={unexpected_columns!r}."
        )
    comparison = f"{group1} vs {group2}"
    _validate_exact_text_column(table, "comparison", comparison)
    _validate_exact_text_column(table, "group1", group1)
    _validate_exact_text_column(table, "group2", group2)

    numeric_columns = [column for column in scientific_columns if column != "is_de_fdr"]
    numeric = {column: _finite_numeric_column(table, column) for column in numeric_columns}
    for column in ("scale1", "scale2", "raw_mean1", "raw_mean2", "raw_normalized_mean1", "raw_normalized_mean2"):
        if bool((numeric[column] < 0).any()):
            raise RuntimeError(f"scVI backend column {column!r} cannot contain negative values.")
    for column in ("non_zeros_proportion1", "non_zeros_proportion2"):
        if bool(((numeric[column] < 0) | (numeric[column] > 1)).any()):
            raise RuntimeError(f"scVI backend probability/proportion column {column!r} must lie in [0, 1].")

    if mode == "change":
        for column in ("proba_de", "proba_not_de"):
            if bool(((numeric[column] < 0) | (numeric[column] > 1)).any()):
                raise RuntimeError(f"scVI backend probability column {column!r} must lie in [0, 1].")
        if not bool(np.allclose(numeric["proba_de"] + numeric["proba_not_de"], 1.0, atol=1e-6, rtol=0.0)):
            raise RuntimeError("scVI proba_de and proba_not_de must be complementary within 1e-6.")
        if delta is None or not bool(np.allclose(numeric["delta"], delta, atol=1e-12, rtol=0.0)):
            raise RuntimeError("scVI backend delta values do not match the requested practical-change threshold.")
        if bool((numeric["pseudocounts"] <= 0).any()) or not bool(
            np.allclose(numeric["pseudocounts"], numeric["pseudocounts"][0], atol=0.0, rtol=0.0)
        ):
            raise RuntimeError("scVI realized pseudocounts must be one constant positive value.")
        if bool((numeric["lfc_std"] < 0).any()):
            raise RuntimeError("scVI lfc_std cannot be negative.")
        if bool((numeric["lfc_min"] > numeric["lfc_max"]).any()):
            raise RuntimeError("scVI lfc_min cannot exceed lfc_max.")
        for column in ("lfc_mean", "lfc_median"):
            if bool(((numeric[column] < numeric["lfc_min"]) | (numeric[column] > numeric["lfc_max"])).any()):
                raise RuntimeError(f"scVI {column} must lie within the reported LFC sample range.")
        if not pd.api.types.is_bool_dtype(table["is_de_fdr"].dtype) or bool(table["is_de_fdr"].isna().any()):
            raise TypeError("scVI canonical is_de_fdr must contain non-missing Boolean values.")
        sort_column = "proba_de"
    else:
        for column in ("proba_m1", "proba_m2"):
            if bool(((numeric[column] < 0) | (numeric[column] > 1)).any()):
                raise RuntimeError(f"scVI backend probability column {column!r} must lie in [0, 1].")
        if not bool(np.allclose(numeric["proba_m1"] + numeric["proba_m2"], 1.0, atol=1e-6, rtol=0.0)):
            raise RuntimeError("scVI proba_m1 and proba_m2 must be complementary within 1e-6.")
        sort_column = "bayes_factor"

    feature_order = {gene: index for index, gene in enumerate(fitted)}
    table.insert(0, "gene", genes)
    table["__feature_order"] = [feature_order[gene] for gene in genes]
    table = table.sort_values(
        [sort_column, "__feature_order"], ascending=[False, True], kind="mergesort", ignore_index=True
    ).drop(columns="__feature_order")
    ordered = ["gene", "comparison", "group1", "group2", *scientific_columns]
    return table.loc[:, ordered].copy()


def _top_effects(table: DataFrame, *, mode: str, limit: int = 10) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in table.head(limit).to_dict(orient="records"):
        if mode == "change":
            rows.append(
                {
                    "gene": record["gene"],
                    "lfc_mean": float(record["lfc_mean"]),
                    "proba_de": float(record["proba_de"]),
                    "bayes_factor": float(record["bayes_factor"]),
                    "posterior_fdr_tagged": bool(record["is_de_fdr"]),
                    "raw_mean_group1": float(record["raw_mean1"]),
                    "raw_mean_group2": float(record["raw_mean2"]),
                }
            )
        else:
            rows.append(
                {
                    "gene": record["gene"],
                    "bayes_factor": float(record["bayes_factor"]),
                    "probability_group1_higher": float(record["proba_m1"]),
                    "raw_mean_group1": float(record["raw_mean1"]),
                    "raw_mean_group2": float(record["raw_mean2"]),
                }
            )
    return rows


def _standalone_accelerator_runtime(model_evidence: Mapping[str, Any]) -> dict[str, Any]:
    diagnostics = model_evidence.get("training_diagnostics")
    if not isinstance(diagnostics, Mapping):
        diagnostics = {}
    runtime: dict[str, Any] = {
        "training_device": str(diagnostics.get("device", "unknown")),
        "cuda_available": None,
        "torch_cuda_compiled_version": None,
        "cudnn_version": None,
        "mps_available": None,
    }
    try:
        import torch
    except ImportError:
        return runtime

    cuda = getattr(torch, "cuda", None)
    is_cuda_available = getattr(cuda, "is_available", None)
    if callable(is_cuda_available):
        runtime["cuda_available"] = bool(is_cuda_available())
    torch_version = getattr(torch, "version", None)
    compiled_cuda = getattr(torch_version, "cuda", None)
    runtime["torch_cuda_compiled_version"] = None if compiled_cuda is None else str(compiled_cuda)
    backends = getattr(torch, "backends", None)
    cudnn = getattr(backends, "cudnn", None)
    cudnn_version = getattr(cudnn, "version", None)
    if callable(cudnn_version):
        value = cudnn_version()
        runtime["cudnn_version"] = None if value is None else str(value)
    mps = getattr(backends, "mps", None)
    is_mps_available = getattr(mps, "is_available", None)
    if callable(is_mps_available):
        runtime["mps_available"] = bool(is_mps_available())
    return runtime


def _analyze_scvi_model_de_evidence_core(
    adata: AnnData,
    model: Any,
    *,
    groupby: str,
    group1: str,
    group2: str,
    subset_column: str | None,
    subset_value: str | None,
    mode: str,
    delta: float | None,
    fdr_target: float | None,
    batch_handling: str,
    n_samples_overall: int,
    random_seed: int,
    openbio_version: str,
) -> tuple[DataFrame, dict[str, Any]]:
    """Run and report one exploratory scVI model-based population contrast."""
    missing_interface = [
        name for name in ("differential_expression_evidence", "var_names", "evidence") if not hasattr(model, name)
    ]
    if missing_interface or not callable(getattr(model, "differential_expression_evidence", None)):
        raise TypeError(
            "scVI Model DE Evidence requires a fitted-model artifact exposing callable "
            "differential_expression_evidence plus var_names and evidence; "
            f"missing={missing_interface!r}."
        )
    run = model.differential_expression_evidence(
        adata,
        groupby=groupby,
        group1=group1,
        group2=group2,
        subset_column=subset_column,
        subset_value=subset_value,
        mode=mode,
        delta=delta,
        fdr_target=fdr_target,
        batch_handling=batch_handling,
        n_samples_overall=n_samples_overall,
        random_seed=random_seed,
    )
    table = canonicalize_scvi_de_table(
        run.table,
        fitted_features=model.var_names,
        group1=group1,
        group2=group2,
        mode=mode,
        delta=delta,
        fdr_target=fdr_target,
    )
    if not isinstance(run.comparison_evidence, Mapping) or not isinstance(run.backend_parameters, Mapping):
        raise TypeError("scVI fitted-model evidence must expose mapping comparison and backend parameters.")
    comparison = copy.deepcopy(dict(run.comparison_evidence))
    model_evidence_value = model.evidence
    if not isinstance(model_evidence_value, Mapping):
        raise TypeError("scVI fitted-model evidence property must be a mapping.")
    model_evidence = copy.deepcopy(dict(model_evidence_value))
    accelerator_runtime = _standalone_accelerator_runtime(model_evidence)
    model_evidence["accelerator_runtime"] = accelerator_runtime
    tagged = int(table["is_de_fdr"].sum()) if mode == "change" else None
    if mode == "change":
        tagged_table = table.loc[table["is_de_fdr"]]
        tagged_positive = tagged_table.loc[tagged_table["lfc_mean"] > 0]
        tagged_negative = tagged_table.loc[tagged_table["lfc_mean"] < 0]
        tagged_zero = tagged_table.loc[tagged_table["lfc_mean"] == 0]
    else:
        tagged_positive = table.iloc[0:0]
        tagged_negative = table.iloc[0:0]
        tagged_zero = table.iloc[0:0]
    results_text = (
        f"The fitted scVI model evaluated {len(table):,} fitted features for {group1!r} versus {group2!r} "
        f"using {comparison['group1_cells']:,} and {comparison['group2_cells']:,} cells, respectively."
    )
    if mode == "change":
        results_text += (
            f" {tagged:,} features were tagged at posterior expected FDR {float(fdr_target):g} "
            f"({len(tagged_positive):,} positive and {len(tagged_negative):,} negative decoded log2-fold-change "
            "directions); this is "
            "model-conditional exploratory evidence, not replicate-aware significance."
        )
    else:
        results_text += " Vanilla mode reports directional posterior evidence and does not provide an FDR tag."

    warnings = [
        "Cells/posterior draws are not independent biological Samples; this output is not formal Condition inference.",
        "Decoded-expression evidence is conditional on the fitted model, fitted feature view, selected cells, and nuisance declarations.",
        "Training completion and finite loss do not prove adequate fit, Technical-batch removal, or biological conservation.",
        "The supplied model is process-local; generated code requires the same live fitted SCVIModel and does not retrain it.",
        "A fixed seed does not guarantee bitwise equality across devices, hardware, or software versions.",
    ]
    if batch_handling == "observed_technical_batches":
        warnings.append(
            "Observed Technical-batch composition remains part of the estimand because decoder batch standardization was disabled."
        )
    if min(comparison["group1_cells"], comparison["group2_cells"]) < 30:
        warnings.append("At least one compared population has fewer than 30 cells; posterior evidence may be unstable.")
    small_batch_rows = [
        row
        for row in comparison["group_by_technical_batch_counts"]
        if (0 < row["group1_cells"] < 5) or (0 < row["group2_cells"] < 5)
    ]
    if small_batch_rows:
        warnings.append("Some supported group-by-Technical-batch cells contain fewer than five cells.")
    if n_samples_overall < 5000:
        warnings.append("Posterior sampling used fewer than the released scvi-tools default of 5,000 draws.")
    if model_evidence["fitted_state_fingerprint_sha256"] is None:
        warnings.append(
            "The backend did not expose a hashable fitted state_dict; exact fitted weights lack an independent hash."
        )
    if len(tagged_zero):
        warnings.append(
            "Some posterior-FDR-tagged features had exactly zero mean decoded log2 fold change and are reported "
            "separately from positive and negative directions."
        )

    parameters = {
        "groupby": groupby,
        "group1": group1,
        "group2": group2,
        "population_scope": comparison["population_scope"],
        "mode": mode,
        "delta": delta,
        "fdr_target": fdr_target,
        "batch_handling": batch_handling,
        "n_samples_overall": n_samples_overall,
        "random_seed": random_seed,
        "fixed_backend_parameters": copy.deepcopy(dict(run.backend_parameters)),
    }
    method = {
        "name": "scVI model-based differential-expression evidence",
        "mode": mode,
        "hypothesis": (
            f"abs(decoded log2 fold change) >= {float(delta):g} with test_mode='two'"
            if mode == "change"
            else "directional decoded expression in group1 greater than group2"
        ),
        "statistical_unit": "selected cells and model posterior mixture; not independent biological Samples",
        "expression_source": "model-owned registered count view and decoded normalized expression; downstream adata.X is ignored",
        "feature_scope": "exact ordered fitted-model feature view",
        "batch_estimand": (
            "both populations decoded on the ordered intersection of actually observed primary Technical batches"
            if batch_handling == "shared_technical_batches"
            else "each cell decoded in its observed primary Technical batch"
        ),
        "posterior_sampling": {"draws": n_samples_overall, "weights": "uniform", "random_seed": random_seed},
        "test_mode": "two" if mode == "change" else None,
        "pseudocount_policy": "released data-driven estimate (pseudocounts=None)" if mode == "change" else None,
        "realized_pseudocount": float(table["pseudocounts"].iloc[0]) if mode == "change" else None,
        "posterior_fdr": (
            {
                "target": fdr_target,
                "meaning": "posterior expected FDR under the fitted model; not a p-value or Benjamini-Hochberg FDR",
            }
            if mode == "change"
            else None
        ),
        "all_stats": True,
        "filter_outlier_cells": False,
        "use_permutation": False,
        "sorting": "descending proba_de" if mode == "change" else "descending bayes_factor",
    }
    limitations = [
        "No Sample-level between-replicate variance is estimated; use pseudobulk edgeR/DESeq2 for Condition inference.",
        "Posterior-FDR tags are not p-values and do not establish frequentist error control across biological Samples.",
        "Results cover the fitted model feature view, often HVGs, and must not be described as full-gene evidence.",
        "Shared-batch decoding standardizes only the primary registered Technical batch, not every nuisance covariate.",
        "The same cells may have contributed to model fitting and clustering, so evidence is exploratory and model-dependent.",
        "Equivalent code requires the same live fitted weights, registered count view, software/device environment, "
        "and observation metadata; it does not reconstruct or retrain the model.",
    ]
    methods_text = (
        "Exploratory scVI model-based differential-expression evidence compared the selected populations over the "
        "fitted model's exact feature view, using its model-owned registered count source and decoded normalized "
        "expression. Posterior draws used uniform weights; the declared Technical batch handling defines the decoder "
        "estimand. This is cell/model-conditional evidence rather than replicate-aware Condition inference."
    )
    summary = {
        "schema_version": SCVI_DE_SUMMARY_SCHEMA,
        "status": "exploratory_model_evidence",
        "node_id": SCVI_DE_NODE_ID,
        "methods": methods_text,
        "results": results_text,
        "method_details": method,
        "comparison": comparison,
        "model_evidence": model_evidence,
        "key_results": {
            "tested_fitted_features": len(table),
            "posterior_fdr_tagged_features": tagged,
            "posterior_fdr_tagged_positive_features": len(tagged_positive) if mode == "change" else None,
            "posterior_fdr_tagged_negative_features": len(tagged_negative) if mode == "change" else None,
            "posterior_fdr_tagged_zero_direction_features": len(tagged_zero) if mode == "change" else None,
            "leading_tagged_positive_effects": _top_effects(tagged_positive, mode=mode),
            "leading_tagged_negative_effects": _top_effects(tagged_negative, mode=mode),
            "top_model_evidence": (_top_effects(table, mode=mode) if mode == "vanilla" or tagged == 0 else []),
            "top_ranked_effects": _top_effects(table, mode=mode),
        },
        "parameters": parameters,
        "references": [dict(reference) for reference in SCVI_DE_REFERENCES],
        "software_versions": collect_software_versions(
            ("anndata", "numpy", "pandas", "scvi-tools", "torch", "lightning"),
            openbio_version=openbio_version,
        ),
        "runtime_environment": {"accelerator": copy.deepcopy(accelerator_runtime)},
        "warnings": warnings,
        "limitations": limitations,
    }
    json.dumps(summary, allow_nan=False, ensure_ascii=False)
    return table, summary


def analyze_scvi_model_de_evidence(
    adata: AnnData,
    model: Any,
    *,
    groupby: str,
    group1: str,
    group2: str,
    subset_column: str | None,
    subset_value: str | None,
    mode: str,
    delta: float | None,
    fdr_target: float | None,
    batch_handling: str,
    n_samples_overall: int,
    random_seed: int,
) -> tuple[DataFrame, dict[str, Any]]:
    return _analyze_scvi_model_de_evidence_core(
        adata,
        model,
        groupby=groupby,
        group1=group1,
        group2=group2,
        subset_column=subset_column,
        subset_value=subset_value,
        mode=mode,
        delta=delta,
        fdr_target=fdr_target,
        batch_handling=batch_handling,
        n_samples_overall=n_samples_overall,
        random_seed=random_seed,
        openbio_version=PLUGIN_VERSION,
    )


def scvi_de_code(
    *,
    groupby: str,
    group1: str,
    group2: str,
    subset_column: str | None,
    subset_value: str | None,
    mode: str,
    delta: float | None,
    fdr_target: float | None,
    batch_handling: str,
    n_samples_overall: int,
    random_seed: int,
) -> str:
    """Equivalent source: it intentionally consumes the same fitted process-local model."""
    embedded_functions = (
        _finite_numeric_column,
        _validate_exact_text_column,
        canonicalize_scvi_de_table,
        _top_effects,
        _package_version,
        collect_software_versions,
        _standalone_accelerator_runtime,
        _analyze_scvi_model_de_evidence_core,
    )
    source_sections = [
        "from __future__ import annotations",
        (
            "# Numerical equivalence requires the same live fitted weights, registered count view, "
            "software/device environment, and observation metadata."
        ),
        "import copy",
        "import json",
        "import math",
        "from collections.abc import Mapping, Sequence",
        f"SCVI_DE_SUMMARY_SCHEMA = {SCVI_DE_SUMMARY_SCHEMA!r}",
        f"SCVI_DE_NODE_ID = {SCVI_DE_NODE_ID!r}",
        f"SCVI_DE_REFERENCES = {tuple(dict(reference) for reference in SCVI_DE_REFERENCES)!r}",
        f"_COMMON_NUMERIC_COLUMNS = {_COMMON_NUMERIC_COLUMNS!r}",
        f"_CHANGE_COLUMNS = {_CHANGE_COLUMNS!r}",
        f"_VANILLA_COLUMNS = {_VANILLA_COLUMNS!r}",
        *(inspect.getsource(function).strip() for function in embedded_functions),
        f'''def scvi_model_de_evidence(adata, fitted_model):
    """Return validated evidence using the same live fitted-model artifact, without retraining."""
    return _analyze_scvi_model_de_evidence_core(
        adata,
        fitted_model,
        groupby={groupby!r},
        group1={group1!r},
        group2={group2!r},
        subset_column={subset_column!r},
        subset_value={subset_value!r},
        mode={mode!r},
        delta={delta!r},
        fdr_target={fdr_target!r},
        batch_handling={batch_handling!r},
        n_samples_overall={n_samples_overall!r},
        random_seed={random_seed!r},
        openbio_version={PLUGIN_VERSION!r},
    )''',
    ]
    return "\n\n".join(source_sections) + "\n"


__all__ = [
    "SCVI_DE_NODE_ID",
    "SCVI_DE_REFERENCES",
    "SCVI_DE_SUMMARY_SCHEMA",
    "analyze_scvi_model_de_evidence",
    "canonicalize_scvi_de_table",
    "resolve_scvi_de_mode",
    "resolve_scvi_population_scope",
    "scvi_de_code",
]
