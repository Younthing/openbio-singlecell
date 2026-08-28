from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .score_artifact import (
    SCORE_ARTIFACT_SCHEMA_VERSION,
    SCORE_ARTIFACTS_KEY,
    _canonical_json_sha256,
    _canonical_string_axis,
    score_frame_fingerprint,
    validate_score_artifact,
)

PATHWAY_SCORE_CONTRAST_COLUMNS = [
    "population",
    "pathway",
    "condition_a",
    "condition_b",
    "n_samples_a",
    "n_samples_b",
    "n_cells_a",
    "n_cells_b",
    "mean_a",
    "mean_b",
    "sd_a",
    "sd_b",
    "mean_difference",
    "ci_low",
    "ci_high",
    "t_statistic",
    "degrees_of_freedom",
    "p_value",
    "p_adjusted",
    "rank",
]


def _standalone_pathway_score_contrast(
    adata,
    *,
    sample_key="sample",
    condition_key="condition",
    annotation_key="cell_type",
    population="",
    condition_a="",
    condition_b="",
    score_key="aucell_scores",
    min_cells_per_sample_population=10,
    min_samples_per_condition=3,
    technical_batch_key="",
    confidence_level=0.95,
    annotation_status="unknown",
    openbio_version="unknown",
):
    """Compare one canonical score family using Sample means and a two-sided Welch test."""
    import hashlib
    import importlib.metadata
    import json
    import math
    import platform
    import warnings as runtime_warnings

    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats

    operation = "Pathway score Sample-level Welch contrast"
    if not hasattr(adata, "obs") or not hasattr(adata, "obs_names") or not hasattr(adata, "obsm"):
        raise TypeError(f"{operation} requires an AnnData-like input with obs, obs_names, and obsm.")
    if int(getattr(adata, "n_obs", len(adata.obs))) <= 0:
        raise ValueError(f"{operation} requires at least one observation.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError(f"{operation} annotation_status must be one of 'unknown', 'provisional', or 'curated'.")
    for value, description in (
        (population, "population"),
        (condition_a, "condition_a"),
        (condition_b, "condition_b"),
        (score_key, "score_key"),
    ):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{operation} {description} must be a nonblank, whitespace-canonical string.")
    if condition_a == condition_b:
        raise ValueError(f"{operation} requires two different Conditions.")
    if (
        isinstance(min_cells_per_sample_population, bool)
        or not isinstance(min_cells_per_sample_population, int)
        or min_cells_per_sample_population < 1
    ):
        raise TypeError(f"{operation} min_cells_per_sample_population must be a positive integer.")
    if isinstance(min_samples_per_condition, bool) or not isinstance(min_samples_per_condition, int):
        raise TypeError(f"{operation} min_samples_per_condition must be an integer.")
    if min_samples_per_condition < 2:
        raise ValueError(f"{operation} requires min_samples_per_condition >= 2.")
    if isinstance(confidence_level, bool) or not isinstance(confidence_level, (int, float)):
        raise TypeError(f"{operation} confidence_level must be numeric.")
    confidence_level = float(confidence_level)
    if not math.isfinite(confidence_level) or not 0.0 < confidence_level < 1.0:
        raise ValueError(f"{operation} confidence_level must be finite and strictly between 0 and 1.")

    role_keys = {
        "sample_key": sample_key,
        "condition_key": condition_key,
        "annotation_key": annotation_key,
    }
    if technical_batch_key is None:
        technical_batch_key = ""
    if not isinstance(technical_batch_key, str) or technical_batch_key != technical_batch_key.strip():
        raise ValueError(f"{operation} technical_batch_key must be empty or whitespace-canonical.")
    if technical_batch_key:
        role_keys["technical_batch_key"] = technical_batch_key
    for role, key in role_keys.items():
        if not isinstance(key, str) or not key or key != key.strip():
            raise ValueError(f"{operation} {role} must be a nonblank, whitespace-canonical column name.")
    if len(set(role_keys.values())) != len(role_keys):
        raise ValueError(f"{operation} Sample, Condition, annotation, and Technical batch roles must be distinct.")
    missing_columns = [key for key in role_keys.values() if key not in adata.obs]
    if missing_columns:
        raise ValueError(f"{operation} observation columns not found: {missing_columns}.")

    def normalize_labels(values, role):
        normalized = []
        display_identities = {}
        for position, original in enumerate(values):
            if not pd.api.types.is_scalar(original):
                raise TypeError(f"{operation} {role} value at position {position} must be scalar.")
            try:
                missing = bool(pd.isna(original))
            except (TypeError, ValueError):
                missing = False
            if missing:
                raise ValueError(f"{operation} {role} contains a missing value at position {position}.")
            value = original.item() if isinstance(original, np.generic) else original
            if isinstance(value, str):
                if not value or value != value.strip():
                    raise ValueError(
                        f"{operation} {role} contains a blank or whitespace-padded label at position {position}."
                    )
                label = value
            elif isinstance(value, bool):
                raise TypeError(f"{operation} {role} boolean labels are not supported.")
            elif isinstance(value, int):
                label = str(value)
            elif isinstance(value, float):
                if not math.isfinite(value):
                    raise ValueError(f"{operation} {role} contains a non-finite numeric label.")
                label = str(value)
            else:
                raise TypeError(
                    f"{operation} {role} labels must be strings or finite integer/float scalars; "
                    f"observed {type(value).__name__}."
                )
            identity = (type(value).__module__, type(value).__qualname__, repr(value))
            prior = display_identities.get(label)
            if prior is not None and prior != identity:
                raise ValueError(
                    f"{operation} {role} values collide after display normalization at {label!r}: "
                    f"{prior[2]} versus {identity[2]}."
                )
            display_identities[label] = identity
            normalized.append(label)
        return normalized

    sample_values = normalize_labels(adata.obs[sample_key].tolist(), "Sample")
    condition_values = normalize_labels(adata.obs[condition_key].tolist(), "Condition")
    annotation_values = normalize_labels(adata.obs[annotation_key].tolist(), "annotation")
    technical_batch_values = (
        normalize_labels(adata.obs[technical_batch_key].tolist(), "Technical batch") if technical_batch_key else None
    )

    sample_order = list(dict.fromkeys(sample_values))
    sample_condition = {}
    sample_batches = {}
    for position, (sample, condition) in enumerate(zip(sample_values, condition_values, strict=True)):
        prior_condition = sample_condition.get(sample)
        if prior_condition is not None and prior_condition != condition:
            raise ValueError(
                f"{operation} requires one Condition per Sample; Sample {sample!r} maps to "
                f"both {prior_condition!r} and {condition!r}."
            )
        sample_condition[sample] = condition
        if technical_batch_values is not None:
            batch = technical_batch_values[position]
            batches = sample_batches.setdefault(sample, [])
            if batch not in batches:
                batches.append(batch)

    observed_conditions = list(dict.fromkeys(condition_values))
    missing_conditions = [value for value in (condition_a, condition_b) if value not in observed_conditions]
    if missing_conditions:
        raise ValueError(f"{operation} selected Conditions are not observed: {missing_conditions}.")
    observed_populations = list(dict.fromkeys(annotation_values))
    if population not in observed_populations:
        raise ValueError(f"{operation} population {population!r} is not observed in {annotation_key!r}.")

    scores, score_artifact = validate_score_artifact(
        adata,
        score_key=score_key,
        required_storage="obsm",
        allowed_methods=("AUCell", "GSVA"),
        np=np,
        pd=pd,
    )
    if not scores.index.equals(adata.obs_names):
        raise ValueError(f"{operation} score index must equal the AnnData observation axis in exact order.")
    pathway_names = list(scores.columns)
    score_values = scores.to_numpy(dtype=float, copy=True)
    if not bool(np.isfinite(score_values).all()):
        raise ValueError(f"{operation} score family contains non-finite values.")

    selected_conditions = {condition_a, condition_b}
    eligible_samples = [sample for sample in sample_order if sample_condition[sample] in selected_conditions]
    cell_positions_by_sample = {sample: [] for sample in eligible_samples}
    for position, (sample, condition, annotation) in enumerate(
        zip(sample_values, condition_values, annotation_values, strict=True)
    ):
        if condition in selected_conditions and annotation == population:
            cell_positions_by_sample[sample].append(position)
    selected_cells_before_sample_filter = int(sum(len(values) for values in cell_positions_by_sample.values()))
    if selected_cells_before_sample_filter == 0:
        raise ValueError(
            f"{operation} found no cells for population {population!r} in Conditions "
            f"{condition_a!r} and {condition_b!r}."
        )

    retained_samples = []
    excluded_samples = []
    for sample in eligible_samples:
        cell_count = len(cell_positions_by_sample[sample])
        if cell_count < min_cells_per_sample_population:
            excluded_samples.append(
                {
                    "sample": sample,
                    "condition": sample_condition[sample],
                    "cells": cell_count,
                    "reason": "below_min_cells_per_sample_population",
                }
            )
        else:
            retained_samples.append(sample)
    retained_by_condition = {
        condition_a: [sample for sample in retained_samples if sample_condition[sample] == condition_a],
        condition_b: [sample for sample in retained_samples if sample_condition[sample] == condition_b],
    }
    insufficient = {
        condition: len(samples)
        for condition, samples in retained_by_condition.items()
        if len(samples) < min_samples_per_condition
    }
    if insufficient:
        raise ValueError(
            f"{operation} requires at least {min_samples_per_condition} retained Samples per Condition after "
            f"the population-cell threshold; observed={insufficient}."
        )

    aggregated_rows = []
    for sample in retained_samples:
        positions = cell_positions_by_sample[sample]
        means = score_values[np.asarray(positions, dtype=int), :].mean(axis=0)
        if not bool(np.isfinite(means).all()):
            raise RuntimeError(f"{operation} produced non-finite Sample means for Sample {sample!r}.")
        aggregated_rows.append(means)
    aggregate = pd.DataFrame(aggregated_rows, index=retained_samples, columns=pathway_names, dtype=float)
    aggregate.index.name = "sample"
    if aggregate.shape != (len(retained_samples), len(pathway_names)):
        raise RuntimeError(f"{operation} Sample aggregation shape invariant failed.")

    technical_batch_audit = {
        "technical_batch_key": None,
        "status": "not_audited",
        "formal_interpretation_invalid": None,
        "formal_interpretation_invalid_reasons": ["technical_batch_not_provided"],
        "samples_spanning_batches": [],
        "shared_batches": [],
        "batches_by_condition": {},
        "sample_counts_by_condition_batch": {},
    }
    if technical_batch_key:
        batches_by_condition = {
            condition: list(
                dict.fromkeys(
                    batch
                    for sample in samples
                    for batch in sample_batches[sample]
                )
            )
            for condition, samples in retained_by_condition.items()
        }
        shared_batches = sorted(set(batches_by_condition[condition_a]).intersection(batches_by_condition[condition_b]))
        samples_spanning_batches = [
            {"sample": sample, "technical_batches": sample_batches[sample]}
            for sample in retained_samples
            if len(sample_batches[sample]) > 1
        ]
        formal_interpretation_invalid_reasons = []
        if samples_spanning_batches:
            formal_interpretation_invalid_reasons.append("sample_spans_multiple_technical_batches")
        if not shared_batches:
            formal_interpretation_invalid_reasons.append("perfect_condition_technical_batch_confounding")
        formal_interpretation_invalid = bool(formal_interpretation_invalid_reasons)
        sample_counts_by_condition_batch = {}
        for condition, samples in retained_by_condition.items():
            sample_counts_by_condition_batch[condition] = {
                batch: int(sum(batch in sample_batches[sample] for sample in samples))
                for batch in batches_by_condition[condition]
            }
        technical_batch_audit = {
            "technical_batch_key": technical_batch_key,
            "status": (
                "formal_interpretation_invalid"
                if formal_interpretation_invalid
                else "overlap_present_no_adjustment"
            ),
            "formal_interpretation_invalid": formal_interpretation_invalid,
            "formal_interpretation_invalid_reasons": formal_interpretation_invalid_reasons,
            "samples_spanning_batches": samples_spanning_batches,
            "shared_batches": shared_batches,
            "batches_by_condition": batches_by_condition,
            "sample_counts_by_condition_batch": sample_counts_by_condition_batch,
        }

    rows = []
    backend_runtime_warning_pathways = []
    constant_arm_pathways = []
    samples_a = retained_by_condition[condition_a]
    samples_b = retained_by_condition[condition_b]
    n_cells_a = int(sum(len(cell_positions_by_sample[sample]) for sample in samples_a))
    n_cells_b = int(sum(len(cell_positions_by_sample[sample]) for sample in samples_b))
    for pathway in pathway_names:
        values_a = aggregate.loc[samples_a, pathway].to_numpy(dtype=float, copy=True)
        values_b = aggregate.loc[samples_b, pathway].to_numpy(dtype=float, copy=True)
        variance_a = float(np.var(values_a, ddof=1))
        variance_b = float(np.var(values_b, ddof=1))
        if variance_a == 0.0 or variance_b == 0.0:
            constant_arm_pathways.append(pathway)
        if variance_a == 0.0 and variance_b == 0.0:
            raise ValueError(
                f"{operation} pathway {pathway!r} has zero Sample-level variance in both Conditions; "
                "Welch statistics are not finite and the complete family cannot be tested."
            )
        with runtime_warnings.catch_warnings(record=True) as caught:
            runtime_warnings.simplefilter("always", RuntimeWarning)
            test = stats.ttest_ind(
                values_a,
                values_b,
                equal_var=False,
                nan_policy="raise",
                alternative="two-sided",
            )
            interval = test.confidence_interval(confidence_level=confidence_level)
        if caught:
            backend_runtime_warning_pathways.append(pathway)
        numeric = np.asarray(
            [
                values_a.mean(),
                values_b.mean(),
                values_a.std(ddof=1),
                values_b.std(ddof=1),
                values_a.mean() - values_b.mean(),
                interval.low,
                interval.high,
                test.statistic,
                test.df,
                test.pvalue,
            ],
            dtype=float,
        )
        if not bool(np.isfinite(numeric).all()):
            raise ValueError(
                f"{operation} pathway {pathway!r} yielded a non-finite Welch result; "
                "the complete pathway family was not reduced or partially reported."
            )
        rows.append(
            {
                "population": population,
                "pathway": pathway,
                "condition_a": condition_a,
                "condition_b": condition_b,
                "n_samples_a": len(samples_a),
                "n_samples_b": len(samples_b),
                "n_cells_a": n_cells_a,
                "n_cells_b": n_cells_b,
                "mean_a": float(numeric[0]),
                "mean_b": float(numeric[1]),
                "sd_a": float(numeric[2]),
                "sd_b": float(numeric[3]),
                "mean_difference": float(numeric[4]),
                "ci_low": float(numeric[5]),
                "ci_high": float(numeric[6]),
                "t_statistic": float(numeric[7]),
                "degrees_of_freedom": float(numeric[8]),
                "p_value": float(numeric[9]),
            }
        )
    table = pd.DataFrame.from_records(rows)
    p_values = table["p_value"].to_numpy(dtype=float)
    p_adjusted = np.asarray(stats.false_discovery_control(p_values, method="bh"), dtype=float)
    if p_adjusted.shape != p_values.shape or not bool(np.isfinite(p_adjusted).all()):
        raise RuntimeError(f"{operation} complete-family Benjamini-Hochberg adjustment failed.")
    table["p_adjusted"] = p_adjusted
    table["rank"] = table["p_adjusted"].rank(method="min", ascending=True).astype("int64")
    table = table.sort_values(["p_adjusted", "p_value", "pathway"], kind="mergesort").reset_index(drop=True)
    table = table[PATHWAY_SCORE_CONTRAST_COLUMNS]
    numeric_columns = [
        column
        for column in PATHWAY_SCORE_CONTRAST_COLUMNS
        if column not in {"population", "pathway", "condition_a", "condition_b"}
    ]
    if not bool(np.isfinite(table[numeric_columns].to_numpy(dtype=float)).all()):
        raise RuntimeError(f"{operation} canonical result table contains non-finite values.")
    if len(table) != len(pathway_names) or set(table["pathway"]) != set(pathway_names):
        raise RuntimeError(f"{operation} did not preserve the complete pathway family.")

    def aggregate_fingerprint():
        digest = hashlib.sha256()
        digest.update(b"openbio-singlecell/pathway-score-sample-means/v1\0")
        digest.update(
            json.dumps(retained_samples, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")
        )
        digest.update(
            json.dumps(pathway_names, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")
        )
        counts = [len(cell_positions_by_sample[sample]) for sample in retained_samples]
        digest.update(json.dumps(counts, separators=(",", ":"), allow_nan=False).encode("utf-8"))
        digest.update(np.ascontiguousarray(aggregate.to_numpy(dtype=float), dtype="<f8").tobytes(order="C"))
        return digest.hexdigest()

    retained_sample_set = set(retained_samples)
    selected_observation_ids = [
        str(adata.obs_names[position])
        for position, (sample, condition, annotation) in enumerate(
            zip(sample_values, condition_values, annotation_values, strict=True)
        )
        if sample in retained_sample_set and condition in selected_conditions and annotation == population
    ]
    selected_observation_fingerprint = _canonical_json_sha256(selected_observation_ids)
    smallest = table.sort_values(["p_adjusted", "p_value", "pathway"], kind="mergesort").head(5)
    strongest = (
        table.assign(_absolute=table["mean_difference"].abs())
        .sort_values(["_absolute", "pathway"], ascending=[False, True], kind="mergesort")
        .head(5)
    )
    top_adjusted = [
        {
            "pathway": str(row.pathway),
            "mean_difference": float(row.mean_difference),
            "ci_low": float(row.ci_low),
            "ci_high": float(row.ci_high),
            "p_adjusted": float(row.p_adjusted),
        }
        for row in smallest.itertuples(index=False)
    ]
    strongest_effects = [
        {
            "pathway": str(row.pathway),
            "mean_difference": float(row.mean_difference),
            "ci_low": float(row.ci_low),
            "ci_high": float(row.ci_high),
            "p_adjusted": float(row.p_adjusted),
        }
        for row in strongest.itertuples(index=False)
    ]

    warnings = []
    if annotation_status == "unknown":
        warnings.append("Population annotation status is unknown; interpretation requires annotation review.")
    elif annotation_status == "provisional":
        warnings.append("Population labels are Provisional annotation and are not Curated annotation.")
    else:
        warnings.append("Curated annotation status is a caller assertion unless upstream provenance verifies it.")
    if excluded_samples:
        warnings.append(
            f"Excluded {len(excluded_samples)} Samples with fewer than "
            f"{min_cells_per_sample_population} selected-population cells."
        )
    if min(len(samples_a), len(samples_b)) == 2:
        warnings.append(
            "At least one Condition arm has only two retained biological Samples. Welch statistics are defined, "
            "but uncertainty and variance estimates are fragile; at least three Samples per arm is recommended."
        )
    if technical_batch_key:
        if technical_batch_audit["formal_interpretation_invalid"]:
            warnings.append(
                "The Technical batch audit found Sample/batch multiplicity or perfect Condition/batch confounding. "
                "Arithmetic results are reported, but formal Condition-effect interpretation is invalid without a "
                "design-aware model or corrected experimental design."
            )
        else:
            warnings.append(
                "Technical batch overlap was audited but not adjusted by this two-sample test; residual imbalance can confound effects."
            )
    else:
        warnings.append(
            "Technical batch was not audited or adjusted; use a Sample-level regression when adjustment is needed."
        )
    if constant_arm_pathways:
        warnings.append(
            f"{len(constant_arm_pathways)} pathways had zero Sample-level variance in one arm but finite Welch results."
        )
    if backend_runtime_warning_pathways:
        warnings.append(
            f"SciPy emitted numerically expected RuntimeWarning diagnostics for {len(backend_runtime_warning_pathways)} "
            "constant/near-constant pathway tests; every released result was independently required to be finite."
        )
    warnings.append(
        f"Upstream {score_artifact['method']} scores are method- and resource-specific and are not interchangeable with other score definitions."
    )

    min_p_adjusted = float(table["p_adjusted"].min())
    significant_005 = int((table["p_adjusted"] <= 0.05).sum())
    methods = (
        f"Within population {population!r}, finite cell-level {score_artifact['method']} scores were averaged "
        "arithmetically once per biological Sample after requiring the declared minimum cell support. "
        f"Sample means from Conditions {condition_a!r} and {condition_b!r} were compared with two-sided "
        f"Welch independent-samples t-tests (difference {condition_a} minus {condition_b}); "
        "Benjamini-Hochberg adjustment was applied once across the complete upstream pathway family."
    )
    results = (
        f"The contrast retained {len(samples_a):,} and {len(samples_b):,} biological Samples "
        f"({n_cells_a:,} and {n_cells_b:,} cells) for {condition_a!r} and {condition_b!r}, respectively, "
        f"and tested all {len(pathway_names):,} pathways. The minimum adjusted p-value was "
        f"{min_p_adjusted:.3g}; {significant_005:,} pathways had BH-adjusted p <= 0.05. "
        "Effect estimates are differences of Sample means with the declared Welch confidence intervals."
    )
    if technical_batch_audit["formal_interpretation_invalid"]:
        results += (
            " The Technical batch audit marks formal_interpretation_invalid=true; effect estimates and "
            "p-values must not be reported as an unconfounded Condition effect."
        )
    parameters = {
        "sample_key": sample_key,
        "condition_key": condition_key,
        "annotation_key": annotation_key,
        "population": population,
        "condition_a": condition_a,
        "condition_b": condition_b,
        "score_key": score_key,
        "min_cells_per_sample_population": int(min_cells_per_sample_population),
        "min_samples_per_condition": int(min_samples_per_condition),
        "technical_batch_key": technical_batch_key or None,
        "confidence_level": confidence_level,
        "annotation_status": annotation_status,
        "aggregation": "arithmetic mean per Sample within selected population",
        "test": "two-sided Welch independent-samples t-test",
        "difference_direction": f"{condition_a} - {condition_b}",
        "multiple_testing": "Benjamini-Hochberg across complete pathway family",
    }
    excluded_preview_limit = 100
    retained_sample_support = []
    for sample in retained_samples:
        support = {
            "sample": sample,
            "condition": sample_condition[sample],
            "cells": len(cell_positions_by_sample[sample]),
        }
        if technical_batch_key:
            support["technical_batches"] = sample_batches[sample]
        retained_sample_support.append(support)
    key_results = {
        "inference_unit": "Sample",
        "population": population,
        "annotation_status": annotation_status,
        "upstream_score_method": score_artifact["method"],
        "score_artifact_fingerprint_sha256": score_artifact["artifact_fingerprint_sha256"],
        "score_content_fingerprint_sha256": score_artifact["score_content_sha256"],
        "score_resource": score_artifact["resource"],
        "score_expression": score_artifact["expression"],
        "pathway_count": len(pathway_names),
        "pathway_order": pathway_names,
        "eligible_samples_before_cell_filter": len(eligible_samples),
        "retained_samples": len(retained_samples),
        "retained_sample_support_preview": retained_sample_support[:excluded_preview_limit],
        "retained_sample_support_preview_truncated": len(retained_sample_support) > excluded_preview_limit,
        "samples_per_condition": {condition_a: len(samples_a), condition_b: len(samples_b)},
        "cells_per_condition": {condition_a: n_cells_a, condition_b: n_cells_b},
        "selected_cells_before_sample_filter": selected_cells_before_sample_filter,
        "excluded_samples_count": len(excluded_samples),
        "excluded_samples_preview": excluded_samples[:excluded_preview_limit],
        "excluded_samples_preview_truncated": len(excluded_samples) > excluded_preview_limit,
        "sample_aggregation_fingerprint_sha256": aggregate_fingerprint(),
        "selected_observation_count": len(selected_observation_ids),
        "selected_observation_axis_sha256": selected_observation_fingerprint,
        "technical_batch_audit": technical_batch_audit,
        "formal_interpretation_invalid": technical_batch_audit["formal_interpretation_invalid"],
        "constant_one_arm_pathways_count": len(constant_arm_pathways),
        "constant_one_arm_pathways_preview": constant_arm_pathways[:100],
        "minimum_p_adjusted": min_p_adjusted,
        "significant_bh_0_05": significant_005,
        "top_adjusted_results": top_adjusted,
        "strongest_effects": strongest_effects,
    }
    references = [dict(reference) for reference in score_artifact["references"]]
    references.extend(
        [
            {
                "citation": (
                    "Welch BL. The generalization of Student's problem when several different population "
                    "variances are involved. Biometrika. 1947;34:28-35."
                ),
                "url": "https://doi.org/10.1093/biomet/34.1-2.28",
                "kind": "method",
                "doi": "10.1093/biomet/34.1-2.28",
            },
            {
                "citation": (
                    "Benjamini Y, Hochberg Y. Controlling the false discovery rate. Journal of the Royal "
                    "Statistical Society Series B. 1995;57:289-300."
                ),
                "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
                "kind": "method",
                "doi": "10.1111/j.2517-6161.1995.tb02031.x",
            },
            {
                "citation": (
                    "Zimmerman KD, et al. A practical solution to pseudoreplication bias in single-cell studies. "
                    "Nature Communications. 2021;12:738."
                ),
                "url": "https://doi.org/10.1038/s41467-021-21038-1",
                "kind": "practice",
                "doi": "10.1038/s41467-021-21038-1",
            },
            {
                "citation": (
                    "Squair JW, et al. Confronting false discoveries in single-cell differential expression. "
                    "Nature Communications. 2021;12:5692."
                ),
                "url": "https://doi.org/10.1038/s41467-021-25960-2",
                "kind": "practice",
                "doi": "10.1038/s41467-021-25960-2",
            },
        ]
    )

    def package_version(distribution):
        try:
            return importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            return "not-installed"

    software_versions = {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
        "anndata": package_version("anndata"),
        "pandas": package_version("pandas"),
        "numpy": package_version("numpy"),
        "scipy": str(scipy.__version__),
    }
    if score_artifact["method"] in {"AUCell", "GSVA"}:
        software_versions["decoupler"] = package_version("decoupler")
    else:
        software_versions["scanpy"] = package_version("scanpy")
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellPathwayScoreTTest",
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": [
            "A two-sample Welch test cannot adjust Technical batch or other Sample-level covariates.",
            "Sample means discard within-Sample distributional structure and can be imprecise with few Samples.",
            "Cell capture, quality control, annotation error, and cells per Sample can influence the aggregated scores.",
            "A non-significant result does not establish equivalence, and a positive score is not causal pathway activation.",
            "The tested family and biological interpretation are inseparable from the upstream scoring method, resource, expression state, and observation cohort.",
        ],
        "references": references,
        "software_versions": software_versions,
    }
    json.dumps(summary, allow_nan=False)
    return table, summary


def run_pathway_score_contrast(
    adata: Any,
    *,
    sample_key: str = "sample",
    condition_key: str = "condition",
    annotation_key: str = "cell_type",
    population: str,
    condition_a: str,
    condition_b: str,
    score_key: str = "aucell_scores",
    min_cells_per_sample_population: int = 10,
    min_samples_per_condition: int = 3,
    technical_batch_key: str = "",
    confidence_level: float = 0.95,
    annotation_status: str = "unknown",
    openbio_version: str = "unknown",
) -> tuple[Any, dict[str, Any]]:
    return _standalone_pathway_score_contrast(
        adata,
        sample_key=sample_key,
        condition_key=condition_key,
        annotation_key=annotation_key,
        population=population,
        condition_a=condition_a,
        condition_b=condition_b,
        score_key=score_key,
        min_cells_per_sample_population=min_cells_per_sample_population,
        min_samples_per_condition=min_samples_per_condition,
        technical_batch_key=technical_batch_key,
        confidence_level=confidence_level,
        annotation_status=annotation_status,
        openbio_version=openbio_version,
    )


def pathway_score_contrast_code(
    *,
    sample_key: str,
    condition_key: str,
    annotation_key: str,
    population: str,
    condition_a: str,
    condition_b: str,
    score_key: str,
    min_cells_per_sample_population: int,
    min_samples_per_condition: int,
    technical_batch_key: str,
    confidence_level: float,
    annotation_status: str,
    openbio_version: str,
) -> str:
    helper_names = (
        _canonical_json_sha256,
        _canonical_string_axis,
        score_frame_fingerprint,
        validate_score_artifact,
    )
    helpers = "\n\n".join(textwrap.dedent(inspect.getsource(helper)).strip() for helper in helper_names)
    implementation = textwrap.dedent(inspect.getsource(_standalone_pathway_score_contrast)).strip()
    return f"""from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any

SCORE_ARTIFACTS_KEY = {SCORE_ARTIFACTS_KEY!r}
SCORE_ARTIFACT_SCHEMA_VERSION = {SCORE_ARTIFACT_SCHEMA_VERSION!r}
PATHWAY_SCORE_CONTRAST_COLUMNS = {PATHWAY_SCORE_CONTRAST_COLUMNS!r}

{helpers}

{implementation}


def contrast_pathway_scores(adata):
    return _standalone_pathway_score_contrast(
        adata,
        sample_key={sample_key!r},
        condition_key={condition_key!r},
        annotation_key={annotation_key!r},
        population={population!r},
        condition_a={condition_a!r},
        condition_b={condition_b!r},
        score_key={score_key!r},
        min_cells_per_sample_population={min_cells_per_sample_population!r},
        min_samples_per_condition={min_samples_per_condition!r},
        technical_batch_key={technical_batch_key!r},
        confidence_level={confidence_level!r},
        annotation_status={annotation_status!r},
        openbio_version={openbio_version!r},
    )
"""


__all__ = [
    "PATHWAY_SCORE_CONTRAST_COLUMNS",
    "pathway_score_contrast_code",
    "run_pathway_score_contrast",
]
