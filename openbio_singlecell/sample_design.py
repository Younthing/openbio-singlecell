from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .analysis_reporting import _package_version, _plain_json, collect_software_versions, summarize_numeric
from .pseudobulk import (
    _standalone_artifact_fingerprint,
    _standalone_validate_pseudobulk_artifact,
)

EDGER_COLUMNS = [
    "gene",
    "log2_fold_change",
    "log_counts_per_million",
    "quasi_likelihood_f",
    "p_value",
    "p_adjusted",
]
PYDESEQ2_COLUMNS = [
    "gene",
    "base_mean",
    "log2_fold_change",
    "log2_fold_change_standard_error",
    "wald_statistic",
    "p_value",
    "p_adjusted",
]


def _standalone_bh_adjust(p_values):
    """Return exact Benjamini-Hochberg adjusted p-values in original order."""
    import numpy as np

    values = np.asarray(p_values, dtype=float).ravel()
    if values.size == 0:
        return values.copy()
    order = np.argsort(values, kind="mergesort")
    ranked = values[order] * float(values.size) / np.arange(1, values.size + 1, dtype=float)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 0.0, 1.0)
    return adjusted


def _standalone_parse_key_list(value, *, description):
    if not isinstance(value, str):
        raise TypeError(f"{description} must be a comma-separated string.")
    keys = [item.strip() for item in value.split(",") if item.strip()]
    if len(keys) != len(set(keys)):
        raise ValueError(f"{description} contains duplicate keys.")
    return keys


def _standalone_backend_matrix_fingerprint(matrix):
    """Hash the exact validated count content without depending on sparse storage layout."""
    import hashlib

    import numpy as np
    from scipy import sparse

    shape = tuple(int(value) for value in matrix.shape)
    digest = hashlib.sha256()
    digest.update(repr(shape).encode("ascii"))
    rows = matrix.tocsr() if sparse.issparse(matrix) else np.asarray(matrix)
    dtype = rows.dtype
    if np.issubdtype(dtype, np.integer) and not np.issubdtype(dtype, np.bool_):
        canonical_kind = b"signed-int64"
        canonical_dtype = "<i8"
    elif np.issubdtype(dtype, np.floating):
        canonical_kind = b"float64"
        canonical_dtype = "<f8"
    else:
        canonical_kind = f"unsupported:{dtype.str}".encode("ascii", errors="backslashreplace")
        canonical_dtype = dtype
    digest.update(canonical_kind)
    for row_index in range(shape[0]):
        row = rows.getrow(row_index).toarray().ravel() if sparse.issparse(rows) else rows[row_index]
        canonical = np.asarray(row, dtype=canonical_dtype)
        digest.update(int(canonical.size).to_bytes(8, "little", signed=False))
        digest.update(canonical.tobytes(order="C"))
    return digest.hexdigest()


def _standalone_backend_input_snapshot(adata, design, contrast):
    """Snapshot scientific backend inputs before an external Adapter can mutate them."""
    import numpy as np

    return {
        "obs_names": tuple(adata.obs_names),
        "var_names": tuple(adata.var_names),
        "obs": adata.obs.copy(deep=True),
        "var": adata.var.copy(deep=True),
        "count_fingerprint": _standalone_backend_matrix_fingerprint(adata.X),
        "design_columns": tuple(design.columns),
        "design_index": tuple(design.index),
        "design_values": np.asarray(design, dtype=float).copy(),
        "contrast": np.asarray(contrast, dtype=float).copy(),
    }


def _standalone_validate_backend_input_snapshot(snapshot, adata, design, contrast, *, operation):
    """Require an external Adapter to preserve every pre-existing scientific input."""
    import numpy as np
    import pandas as pd

    if tuple(adata.obs_names) != snapshot["obs_names"] or tuple(adata.var_names) != snapshot["var_names"]:
        raise RuntimeError(f"{operation} backend mutated the modeled Sample or tested-gene axis.")
    if _standalone_backend_matrix_fingerprint(adata.X) != snapshot["count_fingerprint"]:
        raise RuntimeError(f"{operation} backend mutated the modeled count values.")
    obs_columns = list(snapshot["obs"].columns)
    var_columns = list(snapshot["var"].columns)
    if any(column not in adata.obs for column in obs_columns) or any(column not in adata.var for column in var_columns):
        raise RuntimeError(f"{operation} backend removed pre-existing Sample or feature metadata.")
    try:
        pd.testing.assert_frame_equal(adata.obs.loc[:, obs_columns], snapshot["obs"], check_exact=True)
        pd.testing.assert_frame_equal(adata.var.loc[:, var_columns], snapshot["var"], check_exact=True)
    except AssertionError as exc:
        raise RuntimeError(f"{operation} backend mutated pre-existing Sample or feature metadata.") from exc
    if (
        tuple(design.columns) != snapshot["design_columns"]
        or tuple(design.index) != snapshot["design_index"]
        or not np.array_equal(np.asarray(design, dtype=float), snapshot["design_values"])
    ):
        raise RuntimeError(f"{operation} backend mutated the encoded design matrix.")
    if not np.array_equal(np.asarray(contrast, dtype=float), snapshot["contrast"]):
        raise RuntimeError(f"{operation} backend mutated the requested contrast vector.")


def _standalone_validate_pydeseq2_dds_snapshot(snapshot, dds, *, operation):
    """Validate the fitted DeseqDataSet that owns the actual PyDESeq2 model inputs."""
    import numpy as np
    import pandas as pd

    if dds is None:
        raise RuntimeError(f"{operation} backend did not expose its fitted DeseqDataSet.")
    required_attributes = ("X", "obs", "var", "obs_names", "var_names", "obsm")
    missing_attributes = [name for name in required_attributes if not hasattr(dds, name)]
    if missing_attributes:
        raise RuntimeError(
            f"{operation} fitted DeseqDataSet is missing required scientific state: {missing_attributes}."
        )
    if tuple(dds.obs_names) != snapshot["obs_names"] or tuple(dds.var_names) != snapshot["var_names"]:
        raise RuntimeError(f"{operation} fitted DeseqDataSet mutated the modeled Sample or tested-gene axis.")
    if _standalone_backend_matrix_fingerprint(dds.X) != snapshot["count_fingerprint"]:
        raise RuntimeError(f"{operation} fitted DeseqDataSet mutated the modeled count values.")
    if not isinstance(dds.obs, pd.DataFrame) or not isinstance(dds.var, pd.DataFrame):
        raise RuntimeError(f"{operation} fitted DeseqDataSet did not expose pandas Sample/feature metadata.")
    obs_columns = list(snapshot["obs"].columns)
    var_columns = list(snapshot["var"].columns)
    if any(column not in dds.obs for column in obs_columns) or any(column not in dds.var for column in var_columns):
        raise RuntimeError(f"{operation} fitted DeseqDataSet removed pre-existing Sample or feature metadata.")
    try:
        pd.testing.assert_frame_equal(dds.obs.loc[:, obs_columns], snapshot["obs"], check_exact=True)
        pd.testing.assert_frame_equal(dds.var.loc[:, var_columns], snapshot["var"], check_exact=True)
    except AssertionError as exc:
        raise RuntimeError(f"{operation} fitted DeseqDataSet mutated pre-existing Sample or feature metadata.") from exc
    try:
        dds_design = dds.obsm["design_matrix"]
    except (KeyError, TypeError) as exc:
        raise RuntimeError(f"{operation} fitted DeseqDataSet omitted its encoded design matrix.") from exc
    if not isinstance(dds_design, pd.DataFrame):
        raise RuntimeError(f"{operation} fitted DeseqDataSet design matrix did not preserve Sample/column identity.")
    if (
        tuple(dds_design.columns) != snapshot["design_columns"]
        or tuple(dds_design.index) != snapshot["design_index"]
        or not np.array_equal(np.asarray(dds_design, dtype=float), snapshot["design_values"])
    ):
        raise RuntimeError(f"{operation} fitted DeseqDataSet mutated the encoded design matrix.")
    required_controls = (
        "fit_type",
        "size_factors_fit_type",
        "control_genes",
        "min_mu",
        "min_disp",
        "max_disp",
        "refit_cooks",
        "min_replicates",
        "beta_tol",
        "quiet",
        "low_memory",
    )
    missing_controls = [name for name in required_controls if not hasattr(dds, name)]
    if missing_controls:
        raise RuntimeError(f"{operation} fitted DeseqDataSet omitted fixed controls: {missing_controls}.")
    controls_preserved = (
        dds.fit_type == "parametric"
        and dds.size_factors_fit_type == "ratio"
        and dds.control_genes is None
        and float(dds.min_mu) == 0.5
        and float(dds.min_disp) == 1e-8
        and float(dds.max_disp) == max(10.0, float(len(snapshot["obs_names"])))
        and bool(dds.refit_cooks) is True
        and int(dds.min_replicates) == 7
        and float(dds.beta_tol) == 1e-8
        and bool(dds.quiet) is False
        and bool(dds.low_memory) is False
    )
    if not controls_preserved:
        raise RuntimeError(f"{operation} fitted DeseqDataSet did not preserve the fixed construction controls.")


def _standalone_normalize_backend_warnings(caught):
    """Aggregate repeated backend warnings into deterministic strict-JSON records."""
    records = []
    positions = {}
    for warning in caught:
        category = getattr(getattr(warning, "category", None), "__name__", "Warning")
        message = " ".join(str(getattr(warning, "message", warning)).split())
        key = (category, message)
        if key in positions:
            records[positions[key]]["count"] += 1
        else:
            positions[key] = len(records)
            records.append({"category": category, "message": message, "count": 1})
    return records


def _standalone_prepare_pseudobulk_design(
    artifact,
    *,
    population,
    reference_condition,
    comparison_condition,
    categorical_covariate_keys="",
    continuous_covariate_keys="",
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
):
    """Purely select, encode, validate, and weak-expression-filter one Sample-level contrast."""
    import math

    import numpy as np
    import pandas as pd
    from scipy import sparse

    operation = "Pseudobulk design preflight"

    def required_label(value, description):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{operation} {description} cannot be empty.")
        return value.strip()

    adata, metadata = _standalone_validate_pseudobulk_artifact(artifact, copy_result=False)
    roles = metadata["role_keys"]
    sample_key = roles["sample"]
    population_key = roles["population"]
    condition_key = roles["condition"]
    technical_batch_key = roles["technical_batch"]
    population = required_label(population, "population")
    reference_condition = required_label(reference_condition, "reference Condition")
    comparison_condition = required_label(comparison_condition, "comparison Condition")
    if reference_condition == comparison_condition:
        raise ValueError(f"{operation} reference and comparison Conditions must be different.")
    categorical_keys = _standalone_parse_key_list(
        categorical_covariate_keys,
        description=f"{operation} categorical covariate keys",
    )
    continuous_keys = _standalone_parse_key_list(
        continuous_covariate_keys,
        description=f"{operation} continuous covariate keys",
    )
    declared = metadata["declared_covariates"]
    undeclared_categorical = sorted(set(categorical_keys) - set(declared["categorical"]))
    undeclared_continuous = sorted(set(continuous_keys) - set(declared["continuous"]))
    if undeclared_categorical or undeclared_continuous:
        raise ValueError(
            f"{operation} covariates were not validated by the aggregation artifact; "
            f"categorical={undeclared_categorical}, continuous={undeclared_continuous}."
        )
    available_populations = list(dict.fromkeys(adata.obs[population_key].tolist()))
    if population not in available_populations:
        raise ValueError(f"{operation} population {population!r} is absent; available={available_populations[:30]}.")
    population_mask = adata.obs[population_key].to_numpy(dtype=object) == population
    population_profile_count = int(population_mask.sum())
    available_conditions = list(dict.fromkeys(adata.obs.loc[population_mask, condition_key].tolist()))
    missing_conditions = [
        condition for condition in (reference_condition, comparison_condition) if condition not in available_conditions
    ]
    if missing_conditions:
        raise ValueError(
            f"{operation} requested Conditions are absent after population selection: {missing_conditions}; "
            f"available={available_conditions}."
        )
    condition_mask = adata.obs[condition_key].isin([reference_condition, comparison_condition]).to_numpy()
    # This single materialization is the scientific model workspace: one population and the two
    # requested Conditions. It replaces the former population copy followed by another Condition copy.
    model_data = adata[population_mask & condition_mask].copy()
    if bool(model_data.obs.duplicated(subset=[sample_key]).any()):
        duplicated = model_data.obs.loc[model_data.obs.duplicated(subset=[sample_key], keep=False), sample_key]
        raise ValueError(
            f"{operation} selected population contains duplicate independent Sample rows: "
            f"{sorted(set(duplicated.tolist()))[:20]}."
        )
    condition_counts = {
        reference_condition: int((model_data.obs[condition_key] == reference_condition).sum()),
        comparison_condition: int((model_data.obs[condition_key] == comparison_condition).sum()),
    }
    counts = model_data.X.tocsr() if sparse.issparse(model_data.X) else np.asarray(model_data.X)
    library_sizes = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    if not bool(np.isfinite(library_sizes).all()) or bool((library_sizes <= 0).any()):
        raise ValueError(f"{operation} selected Sample libraries must have positive finite totals.")

    nuisance_categorical = []
    if technical_batch_key:
        nuisance_categorical.append(technical_batch_key)
    nuisance_categorical.extend(categorical_keys)
    nuisance_categorical = list(dict.fromkeys(nuisance_categorical))
    design_columns = ["intercept"]
    design_vectors = [np.ones(model_data.n_obs, dtype=float)]
    categorical_references = {}
    categorical_levels = {}
    omitted_constant_terms = []
    categorical_condition_overlap = []
    condition_values = model_data.obs[condition_key].to_numpy(dtype=object)
    for key in nuisance_categorical:
        values = model_data.obs[key].tolist()
        if any(not isinstance(value, str) or not value for value in values):
            raise ValueError(f"{operation} categorical term {key!r} contains invalid labels.")
        levels = sorted(set(values))
        categorical_levels[key] = levels
        if len(levels) < 2:
            omitted_constant_terms.append(key)
            continue
        categorical_references[key] = levels[0]
        array = np.asarray(values, dtype=object)
        missing_cells = [
            {"level": level, "condition": condition}
            for level in levels
            for condition in (reference_condition, comparison_condition)
            if not bool(((array == level) & (condition_values == condition)).any())
        ]
        if missing_cells:
            categorical_condition_overlap.append({"key": key, "missing_cells": missing_cells})
        for level in levels[1:]:
            design_columns.append(f"{key}[{level}]")
            design_vectors.append((array == level).astype(float))
    continuous_summaries = {}
    continuous_condition_overlap = []
    for key in continuous_keys:
        try:
            values = pd.to_numeric(model_data.obs[key], errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{operation} continuous term {key!r} must be numeric.") from exc
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} continuous term {key!r} must be finite.")
        if bool(np.allclose(values, values[0], rtol=0.0, atol=0.0)):
            omitted_constant_terms.append(key)
            continuous_summaries[key] = summarize_numeric(values)
            continue
        design_columns.append(key)
        design_vectors.append(values)
        continuous_summaries[key] = summarize_numeric(values)
        reference_values = values[condition_values == reference_condition]
        comparison_values = values[condition_values == comparison_condition]
        if float(reference_values.max()) < float(comparison_values.min()) or float(comparison_values.max()) < float(
            reference_values.min()
        ):
            continuous_condition_overlap.append(
                {
                    "key": key,
                    "reference_range": [float(reference_values.min()), float(reference_values.max())],
                    "comparison_range": [float(comparison_values.min()), float(comparison_values.max())],
                }
            )
    condition_column = f"{condition_key}[{comparison_condition} vs {reference_condition}]"
    design_columns.append(condition_column)
    design_vectors.append((condition_values == comparison_condition).astype(float))
    design_array = np.column_stack(design_vectors).astype(float, copy=False)
    if not bool(np.isfinite(design_array).all()):
        raise ValueError(f"{operation} encoded design contains non-finite values.")
    rank = int(np.linalg.matrix_rank(design_array))
    if rank != design_array.shape[1]:
        raise ValueError(
            f"{operation} design is rank deficient (rank={rank}, columns={design_array.shape[1]}); "
            f"Condition is confounded or nuisance terms are aliased: {design_columns}."
        )
    residual_df = int(design_array.shape[0] - rank)
    if residual_df <= 0:
        raise ValueError(
            f"{operation} design leaves no positive residual degrees of freedom "
            f"(Samples={design_array.shape[0]}, rank={rank})."
        )
    contrast = np.zeros(design_array.shape[1], dtype=float)
    contrast[-1] = 1.0
    if not bool(np.isfinite(contrast).all()) or bool(np.allclose(contrast, 0.0)):
        raise RuntimeError(f"{operation} constructed an invalid comparison-minus-reference contrast.")
    contrast_projection = design_array @ contrast
    if bool(np.allclose(contrast_projection, 0.0)):
        raise ValueError(f"{operation} comparison-minus-reference contrast is not estimable.")
    design = pd.DataFrame(design_array, index=model_data.obs_names.copy(), columns=design_columns)
    design_condition_number = float(np.linalg.cond(design_array))

    for name, value, minimum in (
        ("min_count", min_count, 0),
        ("min_total_count", min_total_count, 0),
        ("large_n", large_n, 0),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{operation} {name} must be an integer >= {minimum}.")
    if (
        isinstance(min_prop, bool)
        or not isinstance(min_prop, (int, float))
        or not math.isfinite(float(min_prop))
        or not 0.0 <= float(min_prop) <= 1.0
    ):
        raise ValueError(f"{operation} min_prop must be finite in [0, 1].")
    min_prop = float(min_prop)
    median_library_size = float(np.median(library_sizes))
    if median_library_size <= 0:
        raise ValueError(f"{operation} median library size must be positive.")
    cpm_cutoff = float(min_count / median_library_size * 1_000_000.0)
    sample_size = np.zeros(model_data.n_vars, dtype=int)
    total_count = np.asarray(counts.sum(axis=0, dtype=np.float64)).ravel()
    for row_index, library_size in enumerate(library_sizes):
        row = counts.getrow(row_index).toarray().ravel() if sparse.issparse(counts) else np.asarray(counts[row_index])
        sample_size += (row / library_size * 1_000_000.0 >= cpm_cutoff).astype(int)
    min_sample_size = float(min(condition_counts.values()))
    if min_sample_size > large_n:
        min_sample_size = float(large_n + (min_sample_size - large_n) * min_prop)
    keep = (sample_size >= (min_sample_size - 1e-14)) & (total_count >= (min_total_count - 1e-14))
    retained_indices = np.flatnonzero(keep)
    if retained_indices.size == 0:
        raise ValueError(f"{operation} weak-expression filtering removed every gene.")
    genes_before = int(model_data.n_vars)
    removed_genes_preview = model_data.var_names[~keep].tolist()[:100]
    # model_data is already a private scientific subset, so feature filtering can happen in place.
    model_data._inplace_subset_var(retained_indices)
    filtered = model_data
    if filtered.n_vars != retained_indices.size or not filtered.var_names.is_unique:
        raise RuntimeError(f"{operation} failed to retain a unique deterministic tested-gene universe.")
    filter_diagnostics = {
        "policy": "decoupler_2_2_filter_by_expr_equivalent",
        "group": condition_key,
        "min_count": min_count,
        "min_total_count": min_total_count,
        "large_n": large_n,
        "min_prop": min_prop,
        "median_library_size": median_library_size,
        "cpm_cutoff": cpm_cutoff,
        "effective_min_sample_size": min_sample_size,
        "genes_before": genes_before,
        "genes_retained": int(filtered.n_vars),
        "genes_removed": int(genes_before - filtered.n_vars),
        "retained_genes": filtered.var_names.tolist(),
        "removed_genes_preview": removed_genes_preview,
        "removed_genes_preview_truncated": int((~keep).sum()) > 100,
    }
    diagnostics = {
        "artifact_fingerprint_sha256": metadata["artifact_fingerprint"],
        "aggregation_backend": metadata["backend"],
        "count_source": metadata["count_source"],
        "annotation_status": metadata["annotation_status"],
        "inference_status": metadata["inference_status"],
        "declarations": metadata["declarations"],
        "role_aliases": metadata["role_aliases"],
        "population": population,
        "available_populations": available_populations,
        "reference_condition": reference_condition,
        "comparison_condition": comparison_condition,
        "effect_direction": "positive log2 fold change means higher in comparison versus reference",
        "condition_counts": condition_counts,
        "modeled_samples": int(model_data.n_obs),
        "modeled_sample_ids": model_data.obs[sample_key].tolist(),
        "excluded_other_condition_profiles": int(population_profile_count - model_data.n_obs),
        "design": {
            "requested_formula_like": "~ " + " + ".join([*nuisance_categorical, *continuous_keys, condition_key]),
            "formula_like": "~ "
            + " + ".join(
                [
                    *[key for key in nuisance_categorical if key not in omitted_constant_terms],
                    *[key for key in continuous_keys if key not in omitted_constant_terms],
                    condition_key,
                ]
            ),
            "columns": design_columns,
            "rank": rank,
            "residual_df": residual_df,
            "condition_number": design_condition_number,
            "categorical_references": categorical_references,
            "categorical_levels": categorical_levels,
            "continuous_summaries": continuous_summaries,
            "omitted_constant_terms": list(dict.fromkeys(omitted_constant_terms)),
            "categorical_condition_overlap": categorical_condition_overlap,
            "continuous_condition_overlap": continuous_condition_overlap,
            "contrast_vector": contrast.tolist(),
            "contrast_column": condition_column,
            "row_order": model_data.obs[sample_key].tolist(),
        },
        "expression_filter": filter_diagnostics,
        "library_size_distribution_before_filter": summarize_numeric(library_sizes),
    }
    formal_invalid_reasons = list(metadata["formal_interpretation_invalid_reasons"])
    if any(count < 3 for count in condition_counts.values()):
        formal_invalid_reasons.append("condition_replication_below_three")
    if categorical_condition_overlap:
        formal_invalid_reasons.append("incomplete_categorical_condition_overlap")
    if continuous_condition_overlap:
        formal_invalid_reasons.append("nonoverlapping_continuous_covariate_ranges")
    formal_invalid_reasons = list(dict.fromkeys(formal_invalid_reasons))
    diagnostics["formal_interpretation_invalid"] = bool(formal_invalid_reasons)
    diagnostics["formal_interpretation_invalid_reasons"] = formal_invalid_reasons
    return {
        "adata": filtered,
        "design": design,
        "contrast": contrast,
        "metadata": metadata,
        "diagnostics": diagnostics,
        "categorical_covariates": categorical_keys,
        "continuous_covariates": continuous_keys,
    }


def _standalone_backend_preflight(module, *, engine):
    """Fail closed on the audited Pertpy 1.3 engine interfaces and runtime dependencies."""
    import importlib
    import inspect as runtime_inspect
    from collections.abc import Mapping
    from importlib import metadata as importlib_metadata

    if module is None:
        try:
            module = importlib.import_module("pertpy")
        except (ImportError, OSError) as exc:
            raise RuntimeError(f"{engine} requires Pertpy 1.3; install the compatible Pertpy extra.") from exc
    pertpy_version_value = getattr(module, "__version__", None)
    if pertpy_version_value is None:
        try:
            pertpy_version_value = importlib_metadata.version("pertpy")
        except importlib_metadata.PackageNotFoundError:
            pertpy_version_value = "unknown"
    pertpy_version = str(pertpy_version_value)
    version_parts = pertpy_version.split(".")
    try:
        version_pair = tuple(int(part) for part in version_parts[:2])
    except ValueError as exc:
        raise RuntimeError(f"{engine} could not parse Pertpy version {pertpy_version!r}.") from exc
    if version_pair != (1, 3):
        raise RuntimeError(f"{engine} requires audited Pertpy >=1.3,<1.4; detected {pertpy_version!r}.")
    class_name = "EdgeR" if engine == "edgeR" else "PyDESeq2"
    backend_class = getattr(getattr(module, "tl", None), class_name, None)
    if backend_class is None:
        raise RuntimeError(f"{engine} requires pertpy.tl.{class_name}.")
    for method_name in ("fit", "test_contrasts"):
        if not callable(getattr(backend_class, method_name, None)):
            raise RuntimeError(f"{engine} Pertpy adapter is missing {method_name}().")
    try:
        init_parameters = runtime_inspect.signature(backend_class).parameters
        fit_parameters = runtime_inspect.signature(backend_class.fit).parameters
        test_parameters = runtime_inspect.signature(backend_class.test_contrasts).parameters
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{engine} could not inspect the Pertpy adapter interface.") from exc
    if "design" not in init_parameters or not any(
        parameter.kind == runtime_inspect.Parameter.VAR_KEYWORD for parameter in fit_parameters.values()
    ):
        raise RuntimeError(f"{engine} Pertpy adapter has an incompatible construction/fit interface.")
    if ("contrasts" not in test_parameters and "contrast" not in test_parameters) or not any(
        parameter.kind == runtime_inspect.Parameter.VAR_KEYWORD for parameter in test_parameters.values()
    ):
        raise RuntimeError(f"{engine} Pertpy adapter has an incompatible test_contrasts interface.")

    override = getattr(module, "_openbio_backend_versions", None)
    if override is not None:
        if not isinstance(override, Mapping):
            raise RuntimeError(f"{engine} fake/backend version disclosure must be a mapping.")
        versions = {str(key): str(value) for key, value in override.items()}
        versions.setdefault("pertpy", pertpy_version)
        if engine == "edgeR":
            required_r_versions = {"rpy2", "R", "edgeR", "BiocParallel", "RhpcBLASctl"}
            missing_r_versions = sorted(required_r_versions - set(versions))
            if missing_r_versions:
                raise RuntimeError(f"edgeR backend version disclosure is incomplete; missing={missing_r_versions}.")
        else:
            pydeseq_parts = versions.get("pydeseq2", "").split(".")
            try:
                pydeseq_pair = tuple(int(part) for part in pydeseq_parts[:2])
            except ValueError as exc:
                raise RuntimeError(f"PyDESeq2 could not parse version {versions.get('pydeseq2')!r}.") from exc
            if pydeseq_pair != (0, 5):
                raise RuntimeError(
                    f"PyDESeq2 requires audited pydeseq2 >=0.5,<0.6; detected {versions.get('pydeseq2')!r}."
                )
        return module, backend_class, versions
    versions = {"pertpy": pertpy_version}
    if engine == "edgeR":
        try:
            versions["rpy2"] = importlib_metadata.version("rpy2")
            from rpy2 import robjects as ro
            from rpy2.robjects.packages import importr

            versions["R"] = str(ro.r("R.version.string")[0])
            for package in ("edgeR", "BiocParallel", "RhpcBLASctl"):
                importr(package)
                versions[package] = str(ro.r(f"as.character(packageVersion('{package}'))")[0])
        except Exception as exc:
            raise RuntimeError(
                "edgeR requires rpy2, a working R runtime, and Bioconductor edgeR (plus BiocParallel and "
                f"RhpcBLASctl); dependency preflight failed ({type(exc).__name__}: {exc})."
            ) from exc
    else:
        try:
            versions["pydeseq2"] = importlib_metadata.version("pydeseq2")
        except importlib_metadata.PackageNotFoundError as exc:
            raise RuntimeError("PyDESeq2 requires pydeseq2 >=0.5,<0.6.") from exc
        pydeseq_parts = versions["pydeseq2"].split(".")
        try:
            pydeseq_pair = tuple(int(part) for part in pydeseq_parts[:2])
        except ValueError as exc:
            raise RuntimeError(f"PyDESeq2 could not parse version {versions['pydeseq2']!r}.") from exc
        if pydeseq_pair != (0, 5):
            raise RuntimeError(f"PyDESeq2 requires audited pydeseq2 >=0.5,<0.6; detected {versions['pydeseq2']!r}.")
    return module, backend_class, versions


def _standalone_engine_software_versions(*, openbio_version, backend_versions, decoupler_version, engine):
    backend_versions = {str(key): str(value) for key, value in backend_versions.items()}
    packages = ["anndata", "formulaic-contrasts", "numpy", "pandas", "scipy"]
    if engine == "pydeseq2" and "pydeseq2" not in backend_versions:
        packages.append("pydeseq2")
    discovered = collect_software_versions(packages, openbio_version=str(openbio_version))
    return {
        "python": discovered.pop("python"),
        "openbio-singlecell": discovered.pop("openbio-singlecell"),
        **backend_versions,
        "decoupler": str(decoupler_version),
        **{package: version for package, version in discovered.items() if package not in backend_versions},
    }


def _standalone_engine_summary(diagnostics):
    """Build the complete strict EdgeR or PyDESeq2 scientific report."""
    import json
    from collections.abc import Mapping

    if not isinstance(diagnostics, Mapping):
        raise TypeError("Pseudobulk contrast summary diagnostics must be a mapping.")
    engine = diagnostics.get("engine")
    if engine not in {"edger", "pydeseq2"}:
        raise ValueError(f"Unsupported pseudobulk contrast summary engine: {engine!r}.")
    for field in ("parameters", "warnings", "software_versions"):
        if field not in diagnostics:
            raise ValueError(f"Pseudobulk contrast summary diagnostics are missing {field!r}.")
    key_results = {
        str(key): _plain_json(value)
        for key, value in diagnostics.items()
        if key not in {"parameters", "warnings", "software_versions", "engine"}
    }
    parameters = _plain_json(diagnostics["parameters"])
    warnings = [str(value) for value in diagnostics["warnings"]]
    direction = (
        f"{parameters['comparison_condition']} versus {parameters['reference_condition']} within population "
        f"{parameters['population']!r}"
    )
    condition_counts = key_results["condition_counts"]
    if engine == "edger":
        methods = (
            "Selected one population from a fingerprint-validated Sample-by-population artifact whose count source "
            "and biological Sample identity are analyst declarations, built "
            "and rank-checked an additive Sample-level design, and applied the documented Decoupler 2.2 "
            "filter_by_expr-equivalent prefilter on the modeled libraries. Pertpy 1.3 executed edgeR DGEList, TMM "
            "calcNormFactors, estimateDisp, robust glmQLFit, glmQLFTest for the explicit comparison-minus-reference "
            "vector, and topTags with Benjamini-Hochberg FDR over every retained gene."
        )
        method_name = "edgeR quasi-likelihood"
        node_id = "OpenBioSingleCellPseudobulkEdgeR"
        references = [
            {
                "citation": "Robinson MD, McCarthy DJ, Smyth GK. edgeR. Bioinformatics. 2010;26:139-140.",
                "url": "https://doi.org/10.1093/bioinformatics/btp616",
                "kind": "method",
                "doi": "10.1093/bioinformatics/btp616",
            },
            {
                "citation": "Robinson MD, Oshlack A. A scaling normalization method for RNA-seq data. Genome Biology. 2010;11:R25.",
                "url": "https://doi.org/10.1186/gb-2010-11-3-r25",
                "kind": "method",
                "doi": "10.1186/gb-2010-11-3-r25",
            },
            {
                "citation": "Lun ATL, Chen Y, Smyth GK. edgeR quasi-likelihood methods. Methods Mol Biol. 2016;1418:391-416.",
                "url": "https://doi.org/10.1007/978-1-4939-3578-9_19",
                "kind": "method",
                "doi": "10.1007/978-1-4939-3578-9_19",
            },
            {
                "citation": "Pertpy developers. Pertpy 1.3 EdgeR adapter source and differential-expression tutorial.",
                "url": "https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html",
                "kind": "software_documentation",
                "doi": None,
            },
        ]
    else:
        realized_dispersion_fit = parameters["realized_policy"]["dispersion_fit_type"]
        realized_size_factor_fit = parameters["realized_policy"]["size_factors_fit_type"]
        methods = (
            "Selected one population from a fingerprint-validated Sample-by-population artifact whose count source "
            "and biological Sample identity are analyst declarations, built "
            "and rank-checked an additive Sample-level design, and applied the documented Decoupler 2.2 "
            "filter_by_expr-equivalent prefilter on the modeled libraries. Pertpy 1.3 executed PyDESeq2 with a "
            f"requested parametric dispersion trend and realized {realized_dispersion_fit} trend, requested ratio "
            f"size factors and realized {realized_size_factor_fit} size-factor fitting, "
            "refit_cooks=True, a negative-binomial Wald test for "
            "the explicit comparison-minus-reference vector, Cook's filtering, independent filtering at alpha, and "
            "unshrunk maximum-likelihood log2 fold changes."
        )
        method_name = "PyDESeq2 Wald"
        node_id = "OpenBioSingleCellPseudobulkDESeq2"
        references = [
            {
                "citation": "Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15:550.",
                "url": "https://doi.org/10.1186/s13059-014-0550-8",
                "kind": "method",
                "doi": "10.1186/s13059-014-0550-8",
            },
            {
                "citation": "Muzellec B, et al. PyDESeq2. Bioinformatics. 2023;39:btad547.",
                "url": "https://doi.org/10.1093/bioinformatics/btad547",
                "kind": "software",
                "doi": "10.1093/bioinformatics/btad547",
            },
            {
                "citation": "Pertpy developers. Pertpy 1.3 PyDESeq2 adapter source and differential-expression tutorial.",
                "url": "https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html",
                "kind": "software_documentation",
                "doi": None,
            },
        ]
    references.extend(
        [
            {
                "citation": "Benjamini Y, Hochberg Y. Controlling the false discovery rate. JRSS B. 1995;57:289-300.",
                "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
                "kind": "method",
                "doi": "10.1111/j.2517-6161.1995.tb02031.x",
            },
            {
                "citation": "Squair JW, et al. Confronting false discoveries in single-cell differential expression. Nature Communications. 2021;12:5692.",
                "url": "https://doi.org/10.1038/s41467-021-25960-2",
                "kind": "practice",
                "doi": "10.1038/s41467-021-25960-2",
            },
            {
                "citation": "Heumos L, et al. Pertpy. Nature Methods. 2025.",
                "url": "https://doi.org/10.1038/s41592-025-02909-7",
                "kind": "software",
                "doi": "10.1038/s41592-025-02909-7",
            },
            {
                "citation": "Decoupler developers. decoupler.pp.filter_by_expr official documentation.",
                "url": "https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.pp.filter_by_expr.html",
                "kind": "software_documentation",
                "doi": None,
            },
            {
                "citation": "Virshup I, et al. anndata: Annotated data. Journal of Open Source Software. 2021;6:4371.",
                "url": "https://doi.org/10.21105/joss.04371",
                "kind": "software",
                "doi": "10.21105/joss.04371",
            },
        ]
    )
    results = (
        f"Within {key_results['annotation_status'].capitalize()} population {parameters['population']!r}, "
        f"{method_name} compared {parameters['comparison_condition']!r} "
        f"(n={int(condition_counts[parameters['comparison_condition']])}) with "
        f"{parameters['reference_condition']!r} (n={int(condition_counts[parameters['reference_condition']])}) "
        f"across independent Samples. {int(key_results['genes_tested']):,} genes were tested; "
        f"{int(key_results['significant_total']):,} met FDR <= {float(parameters['fdr_threshold']):g} and "
        f"absolute log2 fold change >= {float(parameters['min_abs_log2_fold_change']):g} "
        f"({int(key_results['significant_up']):,} higher and {int(key_results['significant_down']):,} lower in "
        f"the comparison Condition)."
    )
    if key_results["formal_interpretation_invalid"]:
        results += (
            " The calculation completed, but formal interpretation is flagged invalid for: "
            + ", ".join(key_results["formal_interpretation_invalid_reasons"])
            + "."
        )
    limitations = [
        "The contrast is associative and does not establish a causal Condition effect.",
        "The bounded additive design does not support paired/repeated Samples, random effects, interactions, splines, or nested technical replicates.",
        "Weak-expression filtering and low replicate counts can materially affect the tested universe, dispersion estimates, power, and FDR.",
        "Population annotation, profile QC, and missing Sample-by-population combinations can alter which biological Samples are modeled.",
        "Reported log2 fold changes are unshrunk and may be unstable for low-count genes or small Sample sizes.",
    ]
    summary = {
        "schema_version": 1,
        "node_id": node_id,
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": limitations,
        "references": references,
        "software_versions": _plain_json(diagnostics["software_versions"]),
    }
    if direction not in summary["key_results"]["effect_direction"]:
        summary["key_results"]["comparison_statement"] = direction
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def _standalone_run_edger(
    artifact,
    *,
    population,
    reference_condition,
    comparison_condition,
    categorical_covariate_keys="",
    continuous_covariate_keys="",
    fdr_threshold=0.05,
    min_abs_log2_fold_change=0.0,
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
    openbio_version="not-installed",
    pertpy_module=None,
):
    """Run one audited Pertpy/edgeR quasi-likelihood contrast."""
    import json
    import math

    import numpy as np
    import pandas as pd

    operation = "Pseudobulk edgeR"
    for name, value in (
        ("fdr_threshold", fdr_threshold),
        ("min_abs_log2_fold_change", min_abs_log2_fold_change),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError(f"{operation} {name} must be finite.")
    fdr_threshold = float(fdr_threshold)
    min_abs_log2_fold_change = float(min_abs_log2_fold_change)
    if not 0.0 <= fdr_threshold <= 1.0 or min_abs_log2_fold_change < 0.0:
        raise ValueError(f"{operation} FDR must lie in [0, 1] and the magnitude threshold must be non-negative.")
    prepared = _standalone_prepare_pseudobulk_design(
        artifact,
        population=population,
        reference_condition=reference_condition,
        comparison_condition=comparison_condition,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        min_count=min_count,
        min_total_count=min_total_count,
        large_n=large_n,
        min_prop=min_prop,
    )
    pertpy_module, backend_class, backend_versions = _standalone_backend_preflight(
        pertpy_module,
        engine="edgeR",
    )
    model_data = prepared["adata"]
    # The prepared subset is an engine-specific worker-owned workspace; the backend may consume it directly.
    backend_data = model_data
    backend_design = prepared["design"].copy()
    backend_contrast = np.asarray(prepared["contrast"], dtype=float).copy()
    backend_snapshot = _standalone_backend_input_snapshot(backend_data, backend_design, backend_contrast)
    try:
        model = backend_class(backend_data, design=backend_design)
        model.fit(robust=True)
        backend_result = model.test_contrasts(backend_contrast)
    except Exception as exc:
        raise RuntimeError(f"{operation} backend failed ({type(exc).__name__}: {exc}).") from exc
    _standalone_validate_backend_input_snapshot(
        backend_snapshot,
        backend_data,
        backend_design,
        backend_contrast,
        operation=operation,
    )
    model_adata = getattr(model, "adata", backend_data)
    model_design = getattr(model, "design", backend_design)
    if model_adata is not backend_data or model_design is not backend_design:
        _standalone_validate_backend_input_snapshot(
            backend_snapshot,
            model_adata,
            model_design,
            backend_contrast,
            operation=operation,
        )
    if not isinstance(backend_result, pd.DataFrame):
        raise RuntimeError(f"{operation} backend returned a non-DataFrame result.")
    required = {"variable", "log_fc", "logCPM", "F", "p_value", "adj_p_value"}
    optional = {"contrast"}
    missing = sorted(required - set(backend_result.columns))
    unknown = sorted(set(backend_result.columns) - required - optional)
    if missing or unknown:
        raise RuntimeError(f"{operation} backend result schema mismatch; missing={missing}, unknown={unknown}.")
    genes = backend_result["variable"].tolist()
    if any(not isinstance(gene, str) or not gene for gene in genes) or len(genes) != len(set(genes)):
        raise RuntimeError(f"{operation} backend genes must be unique nonblank strings.")
    expected_genes = list(backend_snapshot["var_names"])
    if set(genes) != set(expected_genes) or len(genes) != len(expected_genes):
        raise RuntimeError(f"{operation} backend result does not contain the exact retained tested-gene universe.")
    numeric_columns = ["log_fc", "logCPM", "F", "p_value", "adj_p_value"]
    numeric = {}
    for column in numeric_columns:
        if pd.api.types.is_bool_dtype(backend_result[column].dtype):
            raise RuntimeError(f"{operation} backend column {column!r} must be numeric and non-boolean.")
        try:
            values = pd.to_numeric(backend_result[column], errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{operation} backend column {column!r} must be numeric.") from exc
        if not bool(np.isfinite(values).all()):
            raise RuntimeError(f"{operation} backend column {column!r} contains non-finite values.")
        numeric[column] = values
    if bool((numeric["F"] < 0).any()):
        raise RuntimeError(f"{operation} quasi-likelihood F statistics must be non-negative.")
    for column in ("p_value", "adj_p_value"):
        if bool(((numeric[column] < 0.0) | (numeric[column] > 1.0)).any()):
            raise RuntimeError(f"{operation} backend {column} values must lie in [0, 1].")
    expected_adjusted = _standalone_bh_adjust(numeric["p_value"])
    if not bool(np.allclose(numeric["adj_p_value"], expected_adjusted, rtol=1e-8, atol=1e-12)):
        raise RuntimeError(
            f"{operation} backend adjusted p-values disagree with Benjamini-Hochberg correction over the "
            "complete retained tested-gene universe."
        )
    table = pd.DataFrame(
        {
            "gene": genes,
            "log2_fold_change": numeric["log_fc"],
            "log_counts_per_million": numeric["logCPM"],
            "quasi_likelihood_f": numeric["F"],
            "p_value": numeric["p_value"],
            "p_adjusted": numeric["adj_p_value"],
        }
    ).sort_values(["p_value", "gene"], kind="mergesort", ignore_index=True)
    calls = (table["p_adjusted"] <= fdr_threshold) & (table["log2_fold_change"].abs() >= min_abs_log2_fold_change)
    significant_up = int((calls & (table["log2_fold_change"] > 0)).sum())
    significant_down = int((calls & (table["log2_fold_change"] < 0)).sum())
    normalization_factors = getattr(model, "normalization_factors", None)
    normalization_summary = {
        "available_from_pertpy_adapter": normalization_factors is not None,
        "distribution": (
            summarize_numeric(normalization_factors) if normalization_factors is not None else None
        ),
    }
    parameters = {
        "engine": "edgeR",
        "population": prepared["diagnostics"]["population"],
        "condition_key": prepared["metadata"]["role_keys"]["condition"],
        "reference_condition": prepared["diagnostics"]["reference_condition"],
        "comparison_condition": prepared["diagnostics"]["comparison_condition"],
        "technical_batch_key": prepared["metadata"]["role_keys"]["technical_batch"],
        "categorical_covariate_keys": prepared["categorical_covariates"],
        "continuous_covariate_keys": prepared["continuous_covariates"],
        "fdr_threshold": fdr_threshold,
        "min_abs_log2_fold_change": min_abs_log2_fold_change,
        "expression_filter": {
            key: prepared["diagnostics"]["expression_filter"][key]
            for key in ("min_count", "min_total_count", "large_n", "min_prop")
        },
        "design_columns": prepared["diagnostics"]["design"]["columns"],
        "contrast_vector": prepared["diagnostics"]["design"]["contrast_vector"],
        "fixed_policy": {
            "normalization": "TMM via Pertpy edgeR calcNormFactors",
            "dispersion": "Pertpy edgeR estimateDisp default controls",
            "ql_fit": "glmQLFit robust=True",
            "test": "glmQLFTest",
            "correction": "Benjamini-Hochberg over every retained tested gene",
            "lfc_shrinkage": False,
            "result_truncation": False,
        },
    }
    warnings = [
        "Weak-expression filtering is design- and library-size-dependent; changing its thresholds changes the tested FDR universe.",
        "Reported log2 fold changes are unshrunk and may be unstable for low-count genes or small Sample sizes.",
        "Pertpy 1.3 does not expose edgeR normalization factors in its public result table; the report retains input library-size diagnostics and records that limitation.",
        "The count source and biological Sample identity are caller declarations; structural checks cannot prove assay provenance or biological independence.",
    ]
    condition_counts = prepared["diagnostics"]["condition_counts"]
    if any(count == 1 for count in condition_counts.values()):
        warnings.append(
            "At least one Condition has only one independent Sample declaration; edgeR can estimate this full-rank design from the remaining residual variation, but formal interpretation is invalid and highly fragile."
        )
    elif any(count < 3 for count in condition_counts.values()):
        warnings.append(
            "At least one Condition has fewer than three independent Sample declarations; calculation remains estimable, but formal interpretation is invalid and fragile."
        )
    elif any(count == 3 for count in condition_counts.values()):
        warnings.append(
            "At least one Condition has exactly three independent Sample declarations; this is executable and does not invalidate formal interpretation, but dispersion and FDR estimates may be sensitive."
        )
    if prepared["diagnostics"]["design"]["omitted_constant_terms"]:
        warnings.append(
            "Constant nuisance terms were omitted as redundant with the intercept: "
            f"{prepared['diagnostics']['design']['omitted_constant_terms']}."
        )
    if prepared["diagnostics"]["design"]["categorical_condition_overlap"]:
        warnings.append(
            "Categorical nuisance levels have incomplete Condition overlap; the design is full rank but relies on partial adjustment support."
        )
    if prepared["diagnostics"]["design"]["continuous_condition_overlap"]:
        warnings.append(
            "Continuous nuisance ranges do not overlap across Conditions; the design is full rank but adjustment relies on extrapolation."
        )
    if prepared["diagnostics"]["role_aliases"]:
        warnings.append(f"Metadata keys were declared for multiple roles: {prepared['diagnostics']['role_aliases']}.")
    if prepared["diagnostics"]["formal_interpretation_invalid"]:
        warnings.append(
            "Formal interpretation is invalid for the completed calculation: "
            + ", ".join(prepared["diagnostics"]["formal_interpretation_invalid_reasons"])
            + "."
        )
    software_versions = _standalone_engine_software_versions(
        openbio_version=openbio_version,
        backend_versions=backend_versions,
        decoupler_version=prepared["metadata"]["backend"]["version"],
        engine="edger",
    )
    diagnostics = {
        "engine": "edger",
        **prepared["diagnostics"],
        "genes_tested": int(len(table)),
        "null_p_value_count": 0,
        "null_p_adjusted_count": 0,
        "significant_total": int(calls.sum()),
        "significant_up": significant_up,
        "significant_down": significant_down,
        "top_genes": _plain_json(table.head(20).to_dict("records")),
        "backend_input_postconditions": {
            "sample_axis_preserved": True,
            "tested_gene_axis_preserved": True,
            "count_values_preserved": True,
            "design_and_contrast_preserved": True,
        },
        "normalization_factors": normalization_summary,
        "backend": {
            "adapter": "pertpy.tl.EdgeR",
            "fit_kwargs": {"robust": True},
            "test_policy": "glmQLFTest + topTags(n=Inf, adjust.method='BH', sort.by='PValue')",
        },
        "parameters": parameters,
        "warnings": warnings,
        "software_versions": software_versions,
    }
    json.dumps(_plain_json(diagnostics), ensure_ascii=False, allow_nan=False)
    _standalone_validate_pseudobulk_artifact(artifact, copy_result=False)
    return table, diagnostics


def _standalone_run_pydeseq2(
    artifact,
    *,
    population,
    reference_condition,
    comparison_condition,
    categorical_covariate_keys="",
    continuous_covariate_keys="",
    fdr_threshold=0.05,
    min_abs_log2_fold_change=0.0,
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
    n_cpus=1,
    openbio_version="not-installed",
    pertpy_module=None,
):
    """Run one audited Pertpy/PyDESeq2 Wald contrast with null preservation."""
    import json
    import math
    import os
    import warnings as runtime_warnings

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from scipy import sparse

    operation = "Pseudobulk PyDESeq2"
    for name, value in (
        ("fdr_threshold", fdr_threshold),
        ("min_abs_log2_fold_change", min_abs_log2_fold_change),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError(f"{operation} {name} must be finite.")
    fdr_threshold = float(fdr_threshold)
    min_abs_log2_fold_change = float(min_abs_log2_fold_change)
    if not 0.0 <= fdr_threshold <= 1.0 or min_abs_log2_fold_change < 0.0:
        raise ValueError(f"{operation} FDR must lie in [0, 1] and the magnitude threshold must be non-negative.")
    try:
        available_cpus = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        available_cpus = os.cpu_count() or 1
    if isinstance(n_cpus, bool) or not isinstance(n_cpus, int) or not 1 <= n_cpus <= available_cpus:
        raise ValueError(f"{operation} n_cpus must be between 1 and the detected {available_cpus} CPUs.")
    prepared = _standalone_prepare_pseudobulk_design(
        artifact,
        population=population,
        reference_condition=reference_condition,
        comparison_condition=comparison_condition,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        min_count=min_count,
        min_total_count=min_total_count,
        large_n=large_n,
        min_prop=min_prop,
    )
    model_data = prepared["adata"]
    dense_matrix_bytes = int(model_data.n_obs) * int(model_data.n_vars) * 8
    estimated_fit_bytes = dense_matrix_bytes * 12
    memory_budget_bytes = 2 * 1024**3
    if estimated_fit_bytes > memory_budget_bytes:
        raise MemoryError(
            f"{operation} estimated dense fit working memory {estimated_fit_bytes / 1024**3:.3f} GiB exceeds "
            "the fixed 2 GiB safety budget; narrow the population or strengthen the documented prefilter."
        )
    pertpy_module, backend_class, backend_versions = _standalone_backend_preflight(
        pertpy_module,
        engine="PyDESeq2",
    )
    # PyDESeq2 requires its own minimal AnnData schema, but the worker-owned count matrix need not be duplicated.
    backend_data = AnnData(
        X=model_data.X,
        obs=pd.DataFrame(
            {"__openbio_sample_identity__": model_data.obs_names.to_numpy(dtype=object)},
            index=model_data.obs_names.copy(),
        ),
        var=pd.DataFrame(
            {"__openbio_feature_identity__": model_data.var_names.to_numpy(dtype=object)},
            index=model_data.var_names.copy(),
        ),
    )
    backend_design = prepared["design"].copy()
    backend_contrast = np.asarray(prepared["contrast"], dtype=float).copy()
    backend_snapshot = _standalone_backend_input_snapshot(backend_data, backend_design, backend_contrast)
    try:
        with runtime_warnings.catch_warnings(record=True) as caught_warnings:
            runtime_warnings.simplefilter("always")
            model = backend_class(backend_data, design=backend_design)
            model.fit(
                fit_type="parametric",
                size_factors_fit_type="ratio",
                control_genes=None,
                min_mu=0.5,
                min_disp=1e-8,
                max_disp=10.0,
                min_replicates=7,
                beta_tol=1e-8,
                n_cpus=n_cpus,
                quiet=False,
                low_memory=False,
            )
            backend_result = model.test_contrasts(
                backend_contrast,
                alpha=fdr_threshold,
                lfc_shrink=None,
                cooks_filter=True,
                independent_filter=True,
                prior_LFC_var=None,
                lfc_null=0.0,
                alt_hypothesis=None,
                quiet=False,
                n_cpus=n_cpus,
            )
    except Exception as exc:
        raise RuntimeError(f"{operation} backend failed ({type(exc).__name__}: {exc}).") from exc
    backend_warning_records = _standalone_normalize_backend_warnings(caught_warnings)
    _standalone_validate_backend_input_snapshot(
        backend_snapshot,
        backend_data,
        backend_design,
        backend_contrast,
        operation=operation,
    )
    model_adata = getattr(model, "adata", backend_data)
    model_design = getattr(model, "design", backend_design)
    if model_adata is not backend_data or model_design is not backend_design:
        _standalone_validate_backend_input_snapshot(
            backend_snapshot,
            model_adata,
            model_design,
            backend_contrast,
            operation=operation,
        )
    dds = getattr(model, "dds", None)
    _standalone_validate_pydeseq2_dds_snapshot(backend_snapshot, dds, operation=operation)
    realized_max_disp = float(dds.max_disp)
    if sparse.issparse(dds.X):
        count_storage = dds.X.tocsc(copy=True)
        count_storage.eliminate_zeros()
        every_retained_gene_contains_zero = bool((count_storage.getnnz(axis=0) < dds.X.shape[0]).all())
    else:
        every_retained_gene_contains_zero = bool((np.asarray(dds.X) == 0).any(axis=0).all())
    realized_size_factors_fit_type = "iterative" if every_retained_gene_contains_zero else "ratio"
    if not isinstance(backend_result, pd.DataFrame):
        raise RuntimeError(f"{operation} backend returned a non-DataFrame result.")
    required = {"variable", "baseMean", "log_fc", "lfcSE", "stat", "p_value", "adj_p_value"}
    optional = {"contrast"}
    missing = sorted(required - set(backend_result.columns))
    unknown = sorted(set(backend_result.columns) - required - optional)
    if missing or unknown:
        raise RuntimeError(f"{operation} backend result schema mismatch; missing={missing}, unknown={unknown}.")
    genes = backend_result["variable"].tolist()
    if any(not isinstance(gene, str) or not gene for gene in genes) or len(genes) != len(set(genes)):
        raise RuntimeError(f"{operation} backend genes must be unique nonblank strings.")
    expected_genes = list(backend_snapshot["var_names"])
    if set(genes) != set(expected_genes) or len(genes) != len(expected_genes):
        raise RuntimeError(f"{operation} backend result does not contain the exact retained tested-gene universe.")
    numeric = {}
    for column in ("baseMean", "log_fc", "lfcSE", "stat", "p_value", "adj_p_value"):
        if pd.api.types.is_bool_dtype(backend_result[column].dtype):
            raise RuntimeError(f"{operation} backend column {column!r} must be numeric and non-boolean.")
        try:
            numeric[column] = pd.to_numeric(backend_result[column], errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{operation} backend column {column!r} must be numeric or null.") from exc
    if not bool(np.isfinite(numeric["baseMean"]).all()) or bool((numeric["baseMean"] < 0).any()):
        raise RuntimeError(f"{operation} base means must be finite and non-negative.")
    finite_lfc_se = numeric["lfcSE"][np.isfinite(numeric["lfcSE"])]
    if bool((finite_lfc_se < 0).any()):
        raise RuntimeError(f"{operation} finite log2 fold-change standard errors must be non-negative.")
    for column in ("log_fc", "lfcSE", "stat"):
        values = numeric[column]
        if bool(np.isinf(values).any()):
            raise RuntimeError(f"{operation} backend column {column!r} contains infinite values.")
    for column in ("p_value", "adj_p_value"):
        values = numeric[column]
        finite = values[np.isfinite(values)]
        if bool(np.isinf(values).any()) or bool(((finite < 0.0) | (finite > 1.0)).any()):
            raise RuntimeError(f"{operation} finite backend {column} values must lie in [0, 1].")
    adjusted_mask = np.isfinite(numeric["adj_p_value"])
    if bool((adjusted_mask & ~np.isfinite(numeric["p_value"])).any()):
        raise RuntimeError(f"{operation} backend returned an adjusted p-value where the raw p-value is null.")
    expected_adjusted = _standalone_bh_adjust(numeric["p_value"][adjusted_mask])
    if not bool(np.allclose(numeric["adj_p_value"][adjusted_mask], expected_adjusted, rtol=1e-8, atol=1e-12)):
        raise RuntimeError(
            f"{operation} backend adjusted p-values disagree with Benjamini-Hochberg correction over the "
            "independent-filter-retained universe."
        )
    table = pd.DataFrame(
        {
            "gene": genes,
            "base_mean": numeric["baseMean"],
            "log2_fold_change": numeric["log_fc"],
            "log2_fold_change_standard_error": numeric["lfcSE"],
            "wald_statistic": numeric["stat"],
            "p_value": numeric["p_value"],
            "p_adjusted": numeric["adj_p_value"],
        }
    ).sort_values(["p_value", "gene"], kind="mergesort", na_position="last", ignore_index=True)
    valid_calls = table["p_adjusted"].notna() & table["log2_fold_change"].notna()
    calls = (
        valid_calls
        & (table["p_adjusted"] <= fdr_threshold)
        & (table["log2_fold_change"].abs() >= min_abs_log2_fold_change)
    )
    p_null_count = int(table["p_value"].isna().sum())
    padj_null_count = int(table["p_adjusted"].isna().sum())
    independent_filter_null_count = int((table["p_value"].notna() & table["p_adjusted"].isna()).sum())
    significant_up = int((calls & (table["log2_fold_change"] > 0)).sum())
    significant_down = int((calls & (table["log2_fold_change"] < 0)).sum())

    size_factors = None
    dispersions = None
    dds_uns = getattr(dds, "uns", {}) if dds is not None else {}
    realized_dispersion_fit = dds_uns.get("disp_function_type") if hasattr(dds_uns, "get") else None
    if realized_dispersion_fit not in {"parametric", "mean"}:
        raise RuntimeError(
            f"{operation} backend did not expose a supported realized dispersion trend; "
            f"observed={realized_dispersion_fit!r}."
        )
    if "size_factors" in getattr(dds, "obsm", {}):
        size_factors = np.asarray(dds.obsm["size_factors"]).ravel()
    elif "size_factors" in getattr(dds, "obs", {}):
        size_factors = np.asarray(dds.obs["size_factors"]).ravel()
    if "dispersions" in getattr(dds, "var", {}):
        dispersions = np.asarray(dds.var["dispersions"]).ravel()
    parameters = {
        "engine": "PyDESeq2",
        "population": prepared["diagnostics"]["population"],
        "condition_key": prepared["metadata"]["role_keys"]["condition"],
        "reference_condition": prepared["diagnostics"]["reference_condition"],
        "comparison_condition": prepared["diagnostics"]["comparison_condition"],
        "technical_batch_key": prepared["metadata"]["role_keys"]["technical_batch"],
        "categorical_covariate_keys": prepared["categorical_covariates"],
        "continuous_covariate_keys": prepared["continuous_covariates"],
        "fdr_threshold": fdr_threshold,
        "alpha": fdr_threshold,
        "min_abs_log2_fold_change": min_abs_log2_fold_change,
        "expression_filter": {
            key: prepared["diagnostics"]["expression_filter"][key]
            for key in ("min_count", "min_total_count", "large_n", "min_prop")
        },
        "design_columns": prepared["diagnostics"]["design"]["columns"],
        "contrast_vector": prepared["diagnostics"]["design"]["contrast_vector"],
        "n_cpus": n_cpus,
        "fixed_policy": {
            "fit_type": "parametric",
            "size_factors_fit_type": "ratio",
            "control_genes": None,
            "min_mu": 0.5,
            "min_disp": 1e-8,
            "max_disp": 10.0,
            "refit_cooks": True,
            "refit_cooks_owner": "pertpy_1_3_adapter",
            "min_replicates": 7,
            "cooks_refit_min_replicates": 7,
            "beta_tol": 1e-8,
            "cooks_filter": True,
            "independent_filter": True,
            "prior_LFC_var": None,
            "lfc_null": 0.0,
            "alt_hypothesis": None,
            "lfc_shrinkage": None,
            "quiet": False,
            "low_memory": False,
            "result_truncation": False,
        },
        "requested_policy": {
            "dispersion_fit_type": "parametric",
            "size_factors_fit_type": "ratio",
        },
        "realized_policy": {
            "dispersion_fit_type": realized_dispersion_fit,
            "size_factors_fit_type": realized_size_factors_fit_type,
            "max_disp": realized_max_disp,
        },
    }
    warnings = [
        "Weak-expression prefiltering and PyDESeq2 independent filtering are distinct stages; changing either changes power or adjusted-p-value availability.",
        "Reported log2 fold changes are unshrunk maximum-likelihood estimates and may be unstable for low-count genes or small Sample sizes.",
        "Pertpy densifies sparse count matrices for PyDESeq2; a conservative fixed working-memory preflight was applied.",
        "The count source and biological Sample identity are caller declarations; structural checks cannot prove assay provenance or biological independence.",
    ]
    for record in backend_warning_records:
        warnings.append(f"PyDESeq2 backend {record['category']} ({record['count']} occurrence(s)): {record['message']}")
    if p_null_count:
        warnings.append(
            f"PyDESeq2 returned {p_null_count:,} null p-values, consistent with Cook's filtering or model-level unavailable tests; rows were retained."
        )
    if independent_filter_null_count:
        warnings.append(
            f"PyDESeq2 independent filtering left {independent_filter_null_count:,} genes with a p-value but null adjusted p-value; rows were retained."
        )
    condition_counts = prepared["diagnostics"]["condition_counts"]
    if any(count == 1 for count in condition_counts.values()):
        warnings.append(
            "At least one Condition has only one independent Sample declaration; PyDESeq2 can execute this full-rank positive-residual design, but formal interpretation is invalid and highly fragile."
        )
    elif any(count < 3 for count in condition_counts.values()):
        warnings.append(
            "At least one Condition has fewer than three independent Sample declarations; calculation remains estimable, but formal interpretation is invalid and fragile."
        )
    elif any(count == 3 for count in condition_counts.values()):
        warnings.append(
            "At least one Condition has exactly three independent Sample declarations; this is executable and does not invalidate formal interpretation, but dispersion and FDR estimates may be sensitive."
        )
    if any(count < 7 for count in condition_counts.values()):
        warnings.append(
            "At least one Condition/design cohort has fewer than seven Sample declarations; the fixed PyDESeq2 min_replicates=7 policy limits Cook outlier replacement/refitting but is not a fit minimum."
        )
    if prepared["diagnostics"]["design"]["omitted_constant_terms"]:
        warnings.append(
            "Constant nuisance terms were omitted as redundant with the intercept: "
            f"{prepared['diagnostics']['design']['omitted_constant_terms']}."
        )
    if prepared["diagnostics"]["design"]["categorical_condition_overlap"]:
        warnings.append(
            "Categorical nuisance levels have incomplete Condition overlap; the design is full rank but relies on partial adjustment support."
        )
    if prepared["diagnostics"]["design"]["continuous_condition_overlap"]:
        warnings.append(
            "Continuous nuisance ranges do not overlap across Conditions; the design is full rank but adjustment relies on extrapolation."
        )
    if prepared["diagnostics"]["role_aliases"]:
        warnings.append(f"Metadata keys were declared for multiple roles: {prepared['diagnostics']['role_aliases']}.")
    if prepared["diagnostics"]["formal_interpretation_invalid"]:
        warnings.append(
            "Formal interpretation is invalid for the completed calculation: "
            + ", ".join(prepared["diagnostics"]["formal_interpretation_invalid_reasons"])
            + "."
        )
    software_versions = _standalone_engine_software_versions(
        openbio_version=openbio_version,
        backend_versions=backend_versions,
        decoupler_version=prepared["metadata"]["backend"]["version"],
        engine="pydeseq2",
    )
    diagnostics = {
        "engine": "pydeseq2",
        **prepared["diagnostics"],
        "genes_tested": int(len(table)),
        "null_p_value_count": p_null_count,
        "null_p_adjusted_count": padj_null_count,
        "independent_filter_null_count": independent_filter_null_count,
        "significant_total": int(calls.sum()),
        "significant_up": significant_up,
        "significant_down": significant_down,
        "top_genes": _plain_json(table.loc[table["p_value"].notna()].head(20).to_dict("records")),
        "size_factor_distribution": (summarize_numeric(size_factors) if size_factors is not None else None),
        "dispersion_distribution": (summarize_numeric(dispersions) if dispersions is not None else None),
        "requested_dispersion_fit_type": "parametric",
        "realized_dispersion_fit_type": realized_dispersion_fit,
        "requested_size_factors_fit_type": "ratio",
        "realized_size_factors_fit_type": realized_size_factors_fit_type,
        "every_retained_gene_contains_zero": every_retained_gene_contains_zero,
        "backend_warning_records": backend_warning_records,
        "backend_input_postconditions": {
            "sample_axis_preserved": True,
            "tested_gene_axis_preserved": True,
            "count_values_preserved": True,
            "design_and_contrast_preserved": True,
            "fitted_dds_sample_axis_preserved": True,
            "fitted_dds_tested_gene_axis_preserved": True,
            "fitted_dds_count_values_preserved": True,
            "fitted_dds_preexisting_obs_var_preserved": True,
            "fitted_dds_design_preserved": True,
            "fitted_dds_fixed_controls_preserved": True,
        },
        "memory": {
            "dense_matrix_bytes": dense_matrix_bytes,
            "estimated_fit_bytes": estimated_fit_bytes,
            "budget_bytes": memory_budget_bytes,
        },
        "backend": {
            "adapter": "pertpy.tl.PyDESeq2",
            "fit_kwargs": {
                "fit_type": "parametric",
                "size_factors_fit_type": "ratio",
                "control_genes": None,
                "min_mu": 0.5,
                "min_disp": 1e-8,
                "max_disp": 10.0,
                "min_replicates": 7,
                "beta_tol": 1e-8,
                "n_cpus": n_cpus,
                "quiet": False,
                "low_memory": False,
            },
            "test_kwargs": {
                "alpha": fdr_threshold,
                "lfc_shrink": None,
                "cooks_filter": True,
                "independent_filter": True,
                "prior_LFC_var": None,
                "lfc_null": 0.0,
                "alt_hypothesis": None,
                "quiet": False,
                "n_cpus": n_cpus,
            },
        },
        "parameters": parameters,
        "warnings": warnings,
        "software_versions": software_versions,
    }
    json.dumps(_plain_json(diagnostics), ensure_ascii=False, allow_nan=False)
    _standalone_validate_pseudobulk_artifact(artifact, copy_result=False)
    return table, diagnostics


def prepare_pseudobulk_design(artifact: Any, **kwargs: Any) -> dict[str, Any]:
    return _standalone_prepare_pseudobulk_design(artifact, **kwargs)


def run_pseudobulk_edger(artifact: Any, **kwargs: Any):
    return _standalone_run_edger(artifact, **kwargs)


def run_pseudobulk_pydeseq2(artifact: Any, **kwargs: Any):
    return _standalone_run_pydeseq2(artifact, **kwargs)


def build_pseudobulk_engine_summary(diagnostics: dict[str, Any]) -> dict[str, Any]:
    return _standalone_engine_summary(diagnostics)


def pseudobulk_engine_code(*, engine: str, **parameters: Any) -> str:
    if engine not in {"edger", "pydeseq2"}:
        raise ValueError(f"Unsupported pseudobulk engine code target: {engine!r}.")
    runner = _standalone_run_edger if engine == "edger" else _standalone_run_pydeseq2
    function_name = "run_pseudobulk_edger" if engine == "edger" else "run_pseudobulk_pydeseq2"
    source_functions = [
        _plain_json,
        _standalone_artifact_fingerprint,
        _standalone_validate_pseudobulk_artifact,
        summarize_numeric,
        _standalone_bh_adjust,
        _standalone_parse_key_list,
        _standalone_backend_matrix_fingerprint,
        _standalone_backend_input_snapshot,
        _standalone_validate_backend_input_snapshot,
    ]
    if engine == "pydeseq2":
        source_functions.append(_standalone_validate_pydeseq2_dds_snapshot)
    source_functions.extend(
        [
            _standalone_normalize_backend_warnings,
            _standalone_prepare_pseudobulk_design,
            _standalone_backend_preflight,
            _package_version,
            collect_software_versions,
            _standalone_engine_software_versions,
            _standalone_engine_summary,
            runner,
        ]
    )
    implementations = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip() for function in source_functions
    )
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in parameters.items())
    return f"""from __future__ import annotations

{implementations}


def {function_name}(pseudobulk):
    table, diagnostics = {runner.__name__}(
        pseudobulk,
        {rendered}
    )
    return table, _standalone_engine_summary(diagnostics)
"""


__all__ = [
    "EDGER_COLUMNS",
    "PYDESEQ2_COLUMNS",
    "build_pseudobulk_engine_summary",
    "prepare_pseudobulk_design",
    "pseudobulk_engine_code",
    "run_pseudobulk_edger",
    "run_pseudobulk_pydeseq2",
]
