from __future__ import annotations

import inspect
import textwrap
from typing import Any


def _standalone_sample_composition_summary(
    adata,
    *,
    sample_key="sample",
    condition_key="condition",
    annotation_key="cell_type",
    annotation_status="unknown",
    max_output_rows=2_000_000,
    openbio_version="unknown",
    _return_diagnostics=False,
):
    """Build one complete Sample-by-annotation descriptive composition grid."""
    import importlib.metadata
    import json
    import math
    import platform

    import numpy as np
    import pandas as pd

    operation = "Sample composition summary"
    if not hasattr(adata, "obs") or not hasattr(adata, "obs_names"):
        raise TypeError(f"{operation} requires an AnnData-like input with obs and obs_names.")
    if int(getattr(adata, "n_obs", len(adata.obs))) <= 0:
        raise ValueError(f"{operation} requires at least one observation.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    if isinstance(max_output_rows, bool) or not isinstance(max_output_rows, int) or max_output_rows <= 0:
        raise TypeError(f"{operation} max_output_rows must be a positive integer.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError(
            f"{operation} annotation_status must be one of 'unknown', 'provisional', or 'curated'."
        )

    role_keys = {
        "sample_key": sample_key,
        "condition_key": condition_key,
        "annotation_key": annotation_key,
    }
    for role, key in role_keys.items():
        if not isinstance(key, str) or not key.strip() or key != key.strip():
            raise ValueError(f"{operation} {role} must be a canonical nonblank column name.")
    if len(set(role_keys.values())) != len(role_keys):
        raise ValueError(f"{operation} Sample, Condition, and annotation columns must be distinct.")
    missing_columns = [key for key in role_keys.values() if key not in adata.obs]
    if missing_columns:
        raise ValueError(f"{operation} observation columns not found: {missing_columns}.")

    def normalize_labels(values, role):
        normalized = []
        identities = {}
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
            prior = identities.get(label)
            if prior is not None and prior != identity:
                raise ValueError(
                    f"{operation} {role} values collide after display normalization at {label!r}: "
                    f"{prior[2]} versus {identity[2]}."
                )
            identities[label] = identity
            normalized.append(label)
        return normalized

    sample_values = normalize_labels(adata.obs[sample_key].tolist(), "Sample")
    condition_values = normalize_labels(adata.obs[condition_key].tolist(), "Condition")

    sample_order = list(dict.fromkeys(sample_values))
    sample_condition = {}
    for sample, condition in zip(sample_values, condition_values, strict=True):
        prior = sample_condition.get(sample)
        if prior is not None and prior != condition:
            raise ValueError(
                f"{operation} requires one Condition per Sample; Sample {sample!r} maps to "
                f"both {prior!r} and {condition!r}."
            )
        sample_condition[sample] = condition
    condition_order = list(dict.fromkeys(sample_condition[sample] for sample in sample_order))

    annotation_values = normalize_labels(adata.obs[annotation_key].tolist(), "annotation")
    observed_annotations = list(dict.fromkeys(annotation_values))
    annotation_series = adata.obs[annotation_key]
    unused_declared_annotations = []
    if isinstance(annotation_series.dtype, pd.CategoricalDtype):
        declared_annotations = normalize_labels(annotation_series.cat.categories.tolist(), "declared annotation category")
        annotation_order = [label for label in declared_annotations if label in set(observed_annotations)]
        unused_declared_annotations = [label for label in declared_annotations if label not in set(observed_annotations)]
    else:
        annotation_order = observed_annotations
    if not annotation_order:
        raise ValueError(f"{operation} found no observed annotation labels.")

    requested_rows = len(sample_order) * len(annotation_order)
    if requested_rows > max_output_rows:
        raise ValueError(
            f"{operation} complete grid requires {requested_rows:,} rows "
            f"({len(sample_order):,} Samples x {len(annotation_order):,} annotations), exceeding "
            f"max_output_rows={max_output_rows:,}."
        )

    counts = {}
    sample_totals = {sample: 0 for sample in sample_order}
    for sample, annotation in zip(sample_values, annotation_values, strict=True):
        key = (sample, annotation)
        counts[key] = counts.get(key, 0) + 1
        sample_totals[sample] += 1
    rows = []
    for sample in sample_order:
        total = sample_totals[sample]
        if total <= 0:
            raise RuntimeError(f"{operation} produced a nonpositive denominator for Sample {sample!r}.")
        for annotation in annotation_order:
            count = int(counts.get((sample, annotation), 0))
            rows.append(
                {
                    "sample": sample,
                    "condition": sample_condition[sample],
                    "annotation": annotation,
                    "cell_count": count,
                    "sample_total_cells": int(total),
                    "proportion": float(count / total),
                }
            )
    table = pd.DataFrame.from_records(
        rows,
        columns=["sample", "condition", "annotation", "cell_count", "sample_total_cells", "proportion"],
    )
    table["sample"] = table["sample"].astype("string")
    table["condition"] = table["condition"].astype("string")
    table["annotation"] = table["annotation"].astype("string")
    table["cell_count"] = table["cell_count"].astype("int64")
    table["sample_total_cells"] = table["sample_total_cells"].astype("int64")
    table["proportion"] = table["proportion"].astype("float64")

    if len(table) != requested_rows or bool(table.duplicated(["sample", "annotation"]).any()):
        raise RuntimeError(f"{operation} complete-grid key invariant failed.")
    if int(table["cell_count"].sum()) != int(adata.n_obs):
        raise RuntimeError(f"{operation} cell-count conservation invariant failed.")
    if bool((table["cell_count"] < 0).any()) or bool((table["sample_total_cells"] <= 0).any()):
        raise RuntimeError(f"{operation} count/denominator range invariant failed.")
    if not bool(np.isfinite(table["proportion"].to_numpy(dtype=float)).all()):
        raise RuntimeError(f"{operation} produced non-finite proportions.")
    if bool(((table["proportion"] < 0.0) | (table["proportion"] > 1.0)).any()):
        raise RuntimeError(f"{operation} proportion range invariant failed.")
    summed_counts = table.groupby("sample", sort=False, observed=True)["cell_count"].sum()
    denominators = table.groupby("sample", sort=False, observed=True)["sample_total_cells"].agg(["min", "max"])
    if any(int(summed_counts.loc[sample]) != sample_totals[sample] for sample in sample_order):
        raise RuntimeError(f"{operation} within-Sample count conservation invariant failed.")
    if bool((denominators["min"] != denominators["max"]).any()):
        raise RuntimeError(f"{operation} denominator consistency invariant failed.")
    proportion_sums = table.groupby("sample", sort=False, observed=True)["proportion"].sum()
    maximum_proportion_sum_error = float(np.max(np.abs(proportion_sums.to_numpy(dtype=float) - 1.0)))
    if maximum_proportion_sum_error > 1e-12:
        raise RuntimeError(f"{operation} within-Sample closure invariant failed.")
    if bool(((table["cell_count"] == 0) & (table["proportion"] != 0.0)).any()):
        raise RuntimeError(f"{operation} zero-count proportion invariant failed.")

    sample_total_array = np.asarray([sample_totals[sample] for sample in sample_order], dtype=float)
    total_quantiles = np.quantile(sample_total_array, [0.0, 0.25, 0.5, 0.75, 1.0])
    samples_per_condition = {
        condition: int(sum(sample_condition[sample] == condition for sample in sample_order))
        for condition in condition_order
    }
    descriptive_limit = 500
    condition_descriptives_count = len(condition_order) * len(annotation_order)
    condition_descriptives = []
    for condition in condition_order:
        if len(condition_descriptives) >= descriptive_limit:
            break
        condition_samples = [sample for sample in sample_order if sample_condition[sample] == condition]
        for annotation in annotation_order:
            if len(condition_descriptives) >= descriptive_limit:
                break
            values = np.asarray(
                [counts.get((sample, annotation), 0) / sample_totals[sample] for sample in condition_samples],
                dtype=float,
            )
            quantiles = np.quantile(values, [0.0, 0.25, 0.5, 0.75, 1.0])
            condition_descriptives.append(
                {
                    "condition": condition,
                    "annotation": annotation,
                    "n_samples": int(values.size),
                    "mean": float(values.mean()),
                    "min": float(quantiles[0]),
                    "q1": float(quantiles[1]),
                    "median": float(quantiles[2]),
                    "q3": float(quantiles[3]),
                    "max": float(quantiles[4]),
                    "zero_samples": int(np.count_nonzero(values == 0.0)),
                }
            )
    mean_proportions = (
        table.groupby("annotation", sort=False, observed=True)["proportion"]
        .mean()
        .sort_values(ascending=False, kind="mergesort")
    )
    leading_annotations = [
        {"annotation": str(annotation), "mean_sample_proportion": float(value)}
        for annotation, value in mean_proportions.iloc[:5].items()
    ]
    zero_count_rows = int((table["cell_count"] == 0).sum())

    warnings = []
    if annotation_status == "unknown":
        warnings.append("Annotation status is unknown; descriptive population labels are not established as curated.")
    elif annotation_status == "provisional":
        warnings.append("Population labels are Provisional annotation and require review before formal interpretation.")
    else:
        warnings.append("Curated annotation status is a caller assertion unless upstream provenance independently verifies it.")
    if len(condition_order) == 1:
        warnings.append("Only one Condition is present; this remains a descriptive composition summary.")
    if any(value < 2 for value in samples_per_condition.values()):
        warnings.append("At least one Condition contains fewer than two Samples; no inferential claim is supported.")

    leading_text = ", ".join(
        f"{record['annotation']} ({record['mean_sample_proportion']:.3f})" for record in leading_annotations[:3]
    ) or "none"
    methods = (
        "All input cells were counted by biological Sample across the complete globally observed annotation universe. "
        "Missing Sample-annotation combinations were retained as zero counts, and proportions used all input cells "
        "within each Sample as the denominator; no smoothing, pseudocount, filtering, or hypothesis test was applied."
    )
    results = (
        f"The complete grid contains {len(sample_order):,} Samples from {len(condition_order):,} Conditions and "
        f"{len(annotation_order):,} annotation categories ({requested_rows:,} rows). Sample depth ranged from "
        f"{int(total_quantiles[0]):,} to {int(total_quantiles[4]):,} cells; {zero_count_rows:,} rows "
        f"({zero_count_rows / requested_rows:.1%}) had zero captured cells. The leading annotations by mean "
        f"Sample proportion were {leading_text}. This is descriptive composition only; no Condition contrast was performed."
    )
    parameters = {
        "sample_key": sample_key,
        "condition_key": condition_key,
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "max_output_rows": int(max_output_rows),
        "denominator": "all input cells within each Sample",
        "missing_metadata_policy": "fail",
        "unused_declared_annotation_policy": "exclude_and_disclose",
    }
    key_results = {
        "input_cells": int(adata.n_obs),
        "samples": len(sample_order),
        "conditions": len(condition_order),
        "annotations": len(annotation_order),
        "sample_order_preview": sample_order[:100],
        "sample_order_preview_truncated": len(sample_order) > 100,
        "condition_order": condition_order,
        "annotation_order_preview": annotation_order[:100],
        "annotation_order_preview_truncated": len(annotation_order) > 100,
        "unused_declared_annotations": unused_declared_annotations[:100],
        "unused_declared_annotations_count": len(unused_declared_annotations),
        "complete_grid_rows": requested_rows,
        "zero_count_rows": zero_count_rows,
        "zero_count_fraction": float(zero_count_rows / requested_rows),
        "samples_per_condition": samples_per_condition,
        "cells_per_sample": {
            "min": int(total_quantiles[0]),
            "q1": float(total_quantiles[1]),
            "median": float(total_quantiles[2]),
            "q3": float(total_quantiles[3]),
            "max": int(total_quantiles[4]),
        },
        "denominator": "all input cells within each Sample",
        "maximum_proportion_sum_error": maximum_proportion_sum_error,
        "annotation_status": annotation_status,
        "leading_annotations": leading_annotations,
        "condition_descriptives": condition_descriptives,
        "condition_descriptives_count": condition_descriptives_count,
        "condition_descriptives_truncated": condition_descriptives_count > descriptive_limit,
    }
    references = [
        {
            "citation": "Virshup I et al. anndata: Annotated data. JOSS. 2021;6:4371.",
            "url": "https://doi.org/10.21105/joss.04371",
            "kind": "software",
            "doi": "10.21105/joss.04371",
        },
        {
            "citation": "The pandas development team. pandas software documentation.",
            "url": "https://pandas.pydata.org/docs/",
            "kind": "software_documentation",
            "doi": None,
        },
        {
            "citation": "Aitchison J. The Statistical Analysis of Compositional Data. JRSS B. 1982;44:139-160.",
            "url": "https://doi.org/10.1111/j.2517-6161.1982.tb01195.x",
            "kind": "method",
            "doi": "10.1111/j.2517-6161.1982.tb01195.x",
        },
        {
            "citation": "Büttner M et al. scCODA is a Bayesian model for compositional single-cell data analysis. 2021.",
            "url": "https://doi.org/10.1038/s41467-021-27150-6",
            "kind": "practice",
            "doi": "10.1038/s41467-021-27150-6",
        },
    ]

    def package_version(distribution):
        try:
            return importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            return "not-installed"

    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellSampleCompositionSummary",
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": [
            "Proportions are closed relative captured-cell composition, not absolute tissue abundance.",
            "Dissociation, capture, quality control, and cells recovered per Sample can change observed composition.",
            "A zero captured count does not prove biological absence.",
            "No statistical Condition contrast was performed; cells are not treated as independent replicates.",
            "Unknown or Provisional annotation does not become Curated annotation through aggregation.",
        ],
        "references": references,
        "software_versions": {
            "python": platform.python_version(),
            "openbio-singlecell": str(openbio_version),
            "anndata": package_version("anndata"),
            "pandas": package_version("pandas"),
            "numpy": package_version("numpy"),
        },
    }
    json.dumps(summary, allow_nan=False)
    if _return_diagnostics:
        return table, summary
    return table


def run_sample_composition_summary(
    adata: Any,
    *,
    sample_key: str = "sample",
    condition_key: str = "condition",
    annotation_key: str = "cell_type",
    annotation_status: str = "unknown",
    max_output_rows: int = 2_000_000,
    openbio_version: str = "unknown",
) -> tuple[Any, dict[str, Any]]:
    return _standalone_sample_composition_summary(
        adata,
        sample_key=sample_key,
        condition_key=condition_key,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        max_output_rows=max_output_rows,
        openbio_version=openbio_version,
        _return_diagnostics=True,
    )


def sample_composition_summary_code(
    *,
    sample_key: str,
    condition_key: str,
    annotation_key: str,
    annotation_status: str,
    max_output_rows: int,
    openbio_version: str,
) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_sample_composition_summary)).strip()
    return f"""from __future__ import annotations

{implementation}


def summarize_sample_composition(adata):
    return _standalone_sample_composition_summary(
        adata,
        sample_key={sample_key!r},
        condition_key={condition_key!r},
        annotation_key={annotation_key!r},
        annotation_status={annotation_status!r},
        max_output_rows={max_output_rows!r},
        openbio_version={openbio_version!r},
        _return_diagnostics=True,
    )
"""


__all__ = [
    "run_sample_composition_summary",
    "sample_composition_summary_code",
]
