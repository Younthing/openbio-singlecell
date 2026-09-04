from __future__ import annotations

import hashlib
import inspect
import json
import textwrap
from typing import Any

PCA_METADATA_COLUMNS = [
    "component",
    "metadata",
    "metadata_type",
    "n_cells_total",
    "n_cells_missing",
    "n_samples_total",
    "n_samples_analyzed",
    "n_samples_missing",
    "n_levels",
    "statistic_name",
    "statistic",
    "estimate_name",
    "estimate",
    "effect_size_name",
    "effect_size",
    "p_value",
    "p_adjusted",
    "significant",
]


def pca_metadata_table_fingerprint(table: Any) -> str:
    import math

    import pandas as pd

    if not isinstance(table, pd.DataFrame) or list(table.columns) != PCA_METADATA_COLUMNS:
        raise ValueError("PCA metadata table fingerprint requires the exact canonical columns.")
    numeric_columns = {
        "n_cells_total",
        "n_cells_missing",
        "n_samples_total",
        "n_samples_analyzed",
        "n_samples_missing",
        "n_levels",
        "statistic",
        "estimate",
        "effect_size",
        "p_value",
        "p_adjusted",
    }
    rows = []
    for row in table.itertuples(index=False, name=None):
        normalized = []
        for column, value in zip(table.columns, row, strict=True):
            if value is None or bool(pd.isna(value)):
                normalized.append(None)
            elif column == "significant":
                normalized.append(bool(value))
            elif column in numeric_columns:
                number = float(value)
                if not math.isfinite(number):
                    raise ValueError("PCA metadata table fingerprint requires finite numeric values.")
                normalized.append(number.hex())
            else:
                normalized.append(str(value))
        rows.append(normalized)
    payload = json.dumps(
        {"schema": "openbio-singlecell/pca-metadata-associations/v1", "rows": rows},
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _standalone_pca_metadata_associations_plot(
    table,
    *,
    producer,
    view=None,
    _return_details=False,
):
    import io
    import re
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    operation = "PCA Metadata Associations Plot"
    if view is None:
        mode = {"view": "association_heatmap"}
    elif not isinstance(view, Mapping):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant == "association_heatmap":
        if set(mode) != {"view"}:
            raise ValueError(f"{operation} heatmap received inactive or unknown parameters.")
    elif variant == "effect_sizes":
        if set(mode) != {"view", "max_associations"}:
            raise ValueError(f"{operation} effect-size view received inactive or unknown parameters.")
        maximum = mode["max_associations"]
        if isinstance(maximum, (bool, np.bool_)) or not isinstance(maximum, (int, np.integer)):
            raise TypeError(f"{operation} max_associations must be an integer.")
        maximum = int(maximum)
        if not 1 <= maximum <= 100:
            raise ValueError(f"{operation} max_associations must be between 1 and 100.")
        mode["max_associations"] = maximum
    else:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")

    if not isinstance(table, pd.DataFrame) or list(table.columns) != PCA_METADATA_COLUMNS:
        raise ValueError(f"{operation} requires the exact canonical PCA metadata-association columns.")
    if table.empty:
        raise ValueError(f"{operation} requires at least one valid stored association.")
    if len(table) > 10_000:
        raise ValueError(f"{operation} supports at most 10,000 stored associations.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError(f"{operation} requires the canonical zero-based row index.")
    if not isinstance(producer, Mapping) or set(producer) != {
        "operation",
        "parameters",
        "input_cells",
        "input_genes",
    }:
        raise ValueError(f"{operation} requires exact producer provenance.")
    if producer["operation"] != "pca_metadata_associations":
        raise ValueError(f"{operation} requires the PCA Metadata Associations producer.")
    parameters = producer["parameters"]
    required_parameters = {
        "use_rep",
        "sample_key",
        "categorical_obs_keys",
        "continuous_obs_keys",
        "alpha",
        "sample_aggregation",
        "continuous_test",
        "categorical_test",
        "categorical_effect_size",
        "p_adjust_method",
        "p_adjust_scope",
        "table_content_fingerprint_sha256",
    }
    if not isinstance(parameters, Mapping) or set(parameters) != required_parameters:
        raise ValueError(f"{operation} producer parameter schema is invalid.")
    if (
        parameters["sample_aggregation"] != "unweighted_mean"
        or parameters["continuous_test"] != "scipy_spearmanr_two_sided"
        or parameters["categorical_test"] != "one_way_anova"
        or parameters["categorical_effect_size"] != "eta_squared"
        or parameters["p_adjust_method"] != "benjamini_hochberg"
        or parameters["p_adjust_scope"] != "all_valid_pc_metadata_pairs"
    ):
        raise ValueError(f"{operation} producer scientific policy is invalid.")
    alpha = parameters["alpha"]
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not np.isfinite(alpha) or not 0 <= alpha <= 1:
        raise ValueError(f"{operation} producer alpha is invalid.")
    for name in ("input_cells", "input_genes"):
        value = producer[name]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{operation} producer {name} must be a nonnegative integer.")

    string_columns = [
        "component",
        "metadata",
        "metadata_type",
        "statistic_name",
        "effect_size_name",
    ]
    for column in string_columns:
        values = table[column].tolist()
        if any(not isinstance(value, str) or not value or value != value.strip() for value in values):
            raise ValueError(f"{operation} {column!r} must contain canonical strings.")
    component_numbers = []
    for component in table["component"]:
        match = re.fullmatch(r"PC([1-9][0-9]*)", component)
        if match is None:
            raise ValueError(f"{operation} components must use canonical one-based PC labels.")
        component_numbers.append(int(match.group(1)))
    if bool(table.duplicated(["component", "metadata"]).any()):
        raise ValueError(f"{operation} contains duplicate PC-by-metadata hypotheses.")
    requested_categorical = parameters["categorical_obs_keys"]
    requested_continuous = parameters["continuous_obs_keys"]
    if not isinstance(requested_categorical, list) or not isinstance(requested_continuous, list):
        raise ValueError(f"{operation} producer metadata axes must be lists.")
    if (
        any(not isinstance(value, str) or not value for value in [*requested_categorical, *requested_continuous])
        or len(set([*requested_categorical, *requested_continuous]))
        != len([*requested_categorical, *requested_continuous])
    ):
        raise ValueError(f"{operation} producer metadata axes are invalid.")
    type_by_metadata = {
        metadata: metadata_type
        for metadata, metadata_type in table[["metadata", "metadata_type"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    }
    if table[["metadata", "metadata_type"]].drop_duplicates()["metadata"].duplicated().any():
        raise ValueError(f"{operation} assigns more than one type to a metadata field.")
    if any(type_ not in {"categorical", "continuous"} for type_ in type_by_metadata.values()):
        raise ValueError(f"{operation} metadata_type must be categorical or continuous.")
    if any(type_by_metadata.get(key) not in {None, "categorical"} for key in requested_categorical) or any(
        type_by_metadata.get(key) not in {None, "continuous"} for key in requested_continuous
    ):
        raise ValueError(f"{operation} metadata types disagree with producer provenance.")
    if any(key not in [*requested_categorical, *requested_continuous] for key in type_by_metadata):
        raise ValueError(f"{operation} contains metadata outside producer provenance.")

    integer_columns = [
        "n_cells_total",
        "n_cells_missing",
        "n_samples_total",
        "n_samples_analyzed",
        "n_samples_missing",
    ]
    integers = {}
    for column in integer_columns:
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_integer_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} must contain canonical integers.")
        integers[column] = series.to_numpy(dtype=np.int64)
    if (
        bool((integers["n_cells_total"] < 1).any())
        or bool((integers["n_cells_missing"] < 0).any())
        or bool((integers["n_samples_total"] < 3).any())
        or bool((integers["n_samples_analyzed"] < 2).any())
        or bool((integers["n_samples_missing"] < 0).any())
        or not bool(
            np.equal(
                integers["n_samples_analyzed"] + integers["n_samples_missing"],
                integers["n_samples_total"],
            ).all()
        )
        or bool((integers["n_cells_missing"] > integers["n_cells_total"]).any())
    ):
        raise ValueError(f"{operation} Sample/cell count columns are inconsistent.")
    if not bool(np.equal(integers["n_cells_total"], producer["input_cells"]).all()):
        raise ValueError(f"{operation} cell counts disagree with producer provenance.")

    for column in ("statistic", "effect_size", "p_value", "p_adjusted"):
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} must contain real numeric values.")
        values = series.to_numpy(dtype=float)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} {column!r} must contain finite values.")
    effect = table["effect_size"].to_numpy(dtype=float)
    raw_p = table["p_value"].to_numpy(dtype=float)
    adjusted_p = table["p_adjusted"].to_numpy(dtype=float)
    if bool(((effect < 0.0) | (effect > 1.0)).any()):
        raise ValueError(f"{operation} effect size must lie in [0, 1].")
    if bool(((raw_p < 0.0) | (raw_p > 1.0)).any()) or bool(
        ((adjusted_p < 0.0) | (adjusted_p > 1.0)).any()
    ):
        raise ValueError(f"{operation} p-values must lie in [0, 1].")
    if bool((adjusted_p + 1e-12 < raw_p).any()):
        raise ValueError(f"{operation} adjusted p-value cannot be smaller than its raw p-value.")
    if not pd.api.types.is_bool_dtype(table["significant"].dtype) or not bool(
        np.array_equal(table["significant"].to_numpy(dtype=bool), adjusted_p <= float(alpha))
    ):
        raise ValueError(f"{operation} significant flags disagree with the producer alpha.")
    for row in table.itertuples(index=False):
        if row.metadata_type == "continuous":
            if (
                row.statistic_name != "spearman_rho"
                or row.estimate_name != "spearman_rho"
                or row.effect_size_name != "rho_squared"
                or row.estimate is None
                or not np.isclose(float(row.estimate), float(row.statistic), rtol=0.0, atol=1e-12)
                or not np.isclose(float(row.effect_size), float(row.statistic) ** 2, rtol=1e-10, atol=1e-12)
            ):
                raise ValueError(f"{operation} continuous association semantics are invalid.")
            if row.n_levels is not None:
                raise ValueError(f"{operation} continuous associations cannot carry n_levels.")
        else:
            if (
                row.statistic_name != "anova_f"
                or row.estimate_name is not None
                or row.estimate is not None
                or row.effect_size_name != "eta_squared"
                or float(row.statistic) < 0.0
                or isinstance(row.n_levels, (bool, np.bool_))
                or not isinstance(row.n_levels, (int, float, np.integer, np.floating))
                or not float(row.n_levels).is_integer()
                or int(row.n_levels) < 2
            ):
                raise ValueError(f"{operation} categorical association semantics are invalid.")
    expected_order = table.sort_values(
        ["p_adjusted", "effect_size", "component", "metadata"],
        ascending=[True, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    if not table.equals(expected_order):
        raise ValueError(f"{operation} rows do not follow the canonical producer evidence order.")
    content_fingerprint = pca_metadata_table_fingerprint(table)
    if parameters["table_content_fingerprint_sha256"] != content_fingerprint:
        raise ValueError(f"{operation} table content fingerprint does not match producer evidence.")

    components = [f"PC{number}" for number in sorted(set(component_numbers))]
    metadata_axis = [
        key for key in [*requested_categorical, *requested_continuous] if key in type_by_metadata
    ]
    effect_size_semantics = {"categorical": "eta_squared", "continuous": "rho_squared"}
    displayed_rows = []
    if variant == "association_heatmap":
        if len(components) > 100 or len(metadata_axis) > 100:
            raise ValueError(f"{operation} heatmap supports at most 100 components and 100 metadata fields.")
        values = np.full((len(metadata_axis), len(components)), np.nan, dtype=float)
        significance = np.zeros(values.shape, dtype=bool)
        component_position = {value: index for index, value in enumerate(components)}
        metadata_position = {value: index for index, value in enumerate(metadata_axis)}
        for row in table.itertuples(index=False):
            position = (metadata_position[row.metadata], component_position[row.component])
            values[position] = float(row.effect_size)
            significance[position] = bool(row.significant)
        figure = Figure(
            figsize=(
                max(7.0, min(20.0, 3.0 + 0.55 * len(components))),
                max(4.5, min(16.0, 2.5 + 0.42 * len(metadata_axis))),
            ),
            constrained_layout=True,
        )
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        masked = np.ma.masked_invalid(values)
        palette = __import__("matplotlib").colormaps["viridis"].with_extremes(bad="#eeeeee")
        image = axis.imshow(masked, aspect="auto", cmap=palette, vmin=0.0, vmax=1.0)
        axis.set_xticks(np.arange(len(components)), components)
        axis.set_yticks(np.arange(len(metadata_axis)), metadata_axis)
        axis.set_xlabel("PCA component")
        axis.set_ylabel("Sample-level metadata")
        axis.set_title("Stored PCA metadata association effect sizes")
        for row, column in zip(*np.where(significance), strict=True):
            axis.text(column, row, "*", ha="center", va="center", color="white", fontsize=10)
        colorbar = figure.colorbar(image, ax=axis)
        colorbar.set_label("Effect size (eta² or rho²)")
        significance_annotation = f"* BH-adjusted p ≤ α (α={float(alpha):.6g})"
        figure.legend(
            [Line2D([], [], marker="$*$", linestyle="none", color="#333333", markersize=10)],
            [significance_annotation.removeprefix("* ")],
            loc="outside lower center",
            frameon=False,
        )
        missing_pairs = int(np.isnan(values).sum())
        title = "PCA metadata association heatmap"
        displayed_rows = table.to_dict(orient="records")
    else:
        selected = table.iloc[: mode["max_associations"]]
        displayed_rows = selected.to_dict(orient="records")
        labels = [f"{row.metadata} · {row.component}" for row in selected.itertuples(index=False)]
        colors = ["#4C78A8" if row.metadata_type == "continuous" else "#F58518" for row in selected.itertuples(index=False)]
        figure = Figure(
            figsize=(8.5, max(4.5, min(16.0, 2.5 + 0.34 * len(selected)))),
            constrained_layout=True,
        )
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        positions = np.arange(len(selected))
        axis.barh(positions, selected["effect_size"].to_numpy(dtype=float), color=colors)
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.set_xlim(0.0, 1.0)
        axis.set_xlabel("Effect size (eta² or rho²)")
        axis.set_title("Canonical PCA metadata evidence order")
        axis.grid(axis="x", alpha=0.2)
        missing_pairs = 0
        significance_annotation = None
        title = "PCA metadata association effect sizes"

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": variant,
        "view_parameters": mode,
        "components": components,
        "metadata": metadata_axis,
        "metadata_types": type_by_metadata,
        "effect_size_semantics": effect_size_semantics,
        "plotted_associations": len(displayed_rows),
        "available_associations": len(table),
        "significant_associations": int(table["significant"].sum()),
        "significance_annotation": significance_annotation,
        "missing_component_metadata_pairs": missing_pairs,
        "displayed_rows": displayed_rows,
        "canonical_order": "p_adjusted ascending, effect_size descending, component then metadata",
        "sample_aggregation": parameters["sample_aggregation"],
        "sample_key": parameters["sample_key"],
        "alpha": float(alpha),
        "input_cells": producer["input_cells"],
        "input_genes": producer["input_genes"],
        "table_content_fingerprint_sha256": content_fingerprint,
        "title": title,
    }
    return (png, details) if _return_details else png


def pca_metadata_associations_plot_code(*, producer: dict[str, Any], view: dict[str, Any]) -> str:
    fingerprint = textwrap.dedent(inspect.getsource(pca_metadata_table_fingerprint)).strip()
    implementation = textwrap.dedent(inspect.getsource(_standalone_pca_metadata_associations_plot)).strip()
    return f"""from __future__ import annotations

import hashlib
import json

PCA_METADATA_COLUMNS = {PCA_METADATA_COLUMNS!r}

{fingerprint}


{implementation}


def plot_pca_metadata_associations(table):
    return _standalone_pca_metadata_associations_plot(
        table,
        producer={producer!r},
        view={view!r},
    )
"""


__all__ = [
    "PCA_METADATA_COLUMNS",
    "_standalone_pca_metadata_associations_plot",
    "pca_metadata_associations_plot_code",
    "pca_metadata_table_fingerprint",
]
