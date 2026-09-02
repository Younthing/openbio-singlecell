from __future__ import annotations

import inspect
import textwrap
from typing import Any


def _standalone_sample_composition_plot(
    table,
    *,
    value="proportion",
    _return_diagnostics=False,
):
    """Render Sample-level composition without aggregating Samples by Condition."""
    import io

    import matplotlib
    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Sample composition plot"
    expected_columns = [
        "sample",
        "condition",
        "annotation",
        "cell_count",
        "sample_total_cells",
        "proportion",
    ]
    if not isinstance(table, pd.DataFrame):
        raise TypeError(f"{operation} requires a pandas DataFrame.")
    if table.columns.tolist() != expected_columns:
        raise ValueError(f"{operation} requires canonical columns {expected_columns!r} in that order.")
    if table.empty:
        raise ValueError(f"{operation} requires at least one canonical composition row.")
    if value not in {"proportion", "cell_count"}:
        raise ValueError(f"{operation} value must be 'proportion' or 'cell_count'.")

    for column in ("sample", "condition", "annotation"):
        labels = table[column].tolist()
        if any(not isinstance(label, str) for label in labels):
            raise TypeError(f"{operation} {column!r} must contain strings only.")
        if any(not label.strip() or label != label.strip() for label in labels):
            raise ValueError(f"{operation} {column!r} values must be canonical nonblank strings.")
    numeric = {}
    for column in ("cell_count", "sample_total_cells"):
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} values must be nonnegative integers.")
        values = series.to_numpy(dtype=float, na_value=np.nan)
        if not bool(np.isfinite(values).all()) or not bool(np.equal(values, np.floor(values)).all()):
            raise ValueError(f"{operation} {column!r} values must be finite nonnegative integers.")
        numeric[column] = values
    if bool((numeric["cell_count"] < 0).any()):
        raise ValueError(f"{operation} 'cell_count' values must be nonnegative integers.")
    if bool((numeric["sample_total_cells"] <= 0).any()):
        raise ValueError(f"{operation} requires one positive denominator per Sample.")
    proportion_series = table["proportion"]
    if pd.api.types.is_bool_dtype(proportion_series.dtype) or not pd.api.types.is_numeric_dtype(
        proportion_series.dtype
    ):
        raise TypeError(f"{operation} 'proportion' values must be finite numbers.")
    proportions = proportion_series.to_numpy(dtype=float, na_value=np.nan)
    if not bool(np.isfinite(proportions).all()):
        raise ValueError(f"{operation} 'proportion' values must be finite.")

    sample_encounter_order = list(dict.fromkeys(table["sample"].tolist()))
    annotation_order = list(dict.fromkeys(table["annotation"].tolist()))
    # ponytail: fixed static-PNG ceiling; replace with pagination or interactive rendering if larger plots are needed.
    max_samples = 100
    max_annotations = 50
    max_segments = 2_000
    segment_count = len(sample_encounter_order) * len(annotation_order)
    if len(sample_encounter_order) > max_samples:
        raise ValueError(
            f"{operation} supports at most {max_samples:,} Samples in one static PNG; "
            f"found {len(sample_encounter_order):,}."
        )
    if len(annotation_order) > max_annotations:
        raise ValueError(
            f"{operation} supports at most {max_annotations:,} annotation categories in one static PNG; "
            f"found {len(annotation_order):,}."
        )
    if segment_count > max_segments:
        raise ValueError(
            f"{operation} supports at most {max_segments:,} stacked bar segments in one static PNG; "
            f"found {segment_count:,}."
        )
    sample_condition = {}
    for sample in sample_encounter_order:
        conditions = list(dict.fromkeys(table.loc[table["sample"] == sample, "condition"].tolist()))
        if len(conditions) != 1:
            raise ValueError(f"{operation} requires one Condition per Sample; Sample {sample!r} is ambiguous.")
        sample_condition[sample] = conditions[0]
    condition_order = list(dict.fromkeys(sample_condition[sample] for sample in sample_encounter_order))
    sample_order = [
        sample
        for condition in condition_order
        for sample in sample_encounter_order
        if sample_condition[sample] == condition
    ]
    expected_grid = [(sample, annotation) for sample in sample_encounter_order for annotation in annotation_order]
    observed_grid = list(table[["sample", "annotation"]].itertuples(index=False, name=None))
    if observed_grid != expected_grid:
        raise ValueError(
            f"{operation} requires the canonical complete Sample-by-annotation grid in nested Sample order."
        )
    for sample in sample_encounter_order:
        mask = table["sample"] == sample
        denominators = numeric["sample_total_cells"][mask.to_numpy()]
        if not bool(np.equal(denominators, denominators[0]).all()):
            raise ValueError(f"{operation} requires one positive denominator per Sample.")
        counts = numeric["cell_count"][mask.to_numpy()]
        if float(counts.sum()) != float(denominators[0]):
            raise ValueError(f"{operation} cell counts must sum to the Sample denominator.")
    expected_proportions = numeric["cell_count"] / numeric["sample_total_cells"]
    if not bool(np.allclose(proportions, expected_proportions, rtol=0.0, atol=1e-12)):
        raise ValueError(f"{operation} proportion must equal cell_count / sample_total_cells for every row.")
    values_by_key = {
        (sample, annotation): float(amount)
        for sample, annotation, amount in table[["sample", "annotation", value]].itertuples(index=False, name=None)
    }

    figure_width = max(6.4, min(20.0, 0.65 * len(sample_order) + 2.5))
    figure = Figure(figsize=(figure_width, 6.0))
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    positions = np.arange(len(sample_order), dtype=float)
    bottom = np.zeros(len(sample_order), dtype=float)
    palette = matplotlib.colormaps["tab20" if len(annotation_order) <= 20 else "turbo"].resampled(
        max(1, len(annotation_order))
    )
    annotation_colors = {
        annotation: matplotlib.colors.to_hex(palette(index)) for index, annotation in enumerate(annotation_order)
    }
    for annotation in annotation_order:
        heights = np.asarray([values_by_key[(sample, annotation)] for sample in sample_order], dtype=float)
        axis.bar(
            positions,
            heights,
            width=0.8,
            bottom=bottom,
            color=annotation_colors[annotation],
            label=annotation,
        )
        bottom += heights

    group_start = 0
    for condition in condition_order:
        group_size = sum(sample_condition[sample] == condition for sample in sample_order)
        group_end = group_start + group_size
        axis.text(
            (group_start + group_end - 1) / 2,
            -0.2,
            condition,
            ha="center",
            va="top",
            transform=axis.get_xaxis_transform(),
        )
        if group_end < len(sample_order):
            axis.axvline(group_end - 0.5, color="#666666", linewidth=0.8)
        group_start = group_end

    axis.set_xticks(positions, sample_order, rotation=45, ha="right")
    y_axis_label = "Proportion" if value == "proportion" else "Cell count"
    axis.set_ylabel(y_axis_label)
    if value == "proportion":
        axis.set_ylim(0.0, 1.0)
    axis.set_title("Sample composition by Condition")
    axis.legend(title="Annotation", bbox_to_anchor=(1.02, 1.0), loc="upper left")
    axis.grid(axis="y", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()

    sample_totals = table.drop_duplicates("sample")["sample_total_cells"]
    annotation_readability_warnings = []
    if len(annotation_order) > 20:
        annotation_readability_warnings.append(
            f"The static PNG contains {len(annotation_order):,} annotation categories; "
            "colors and legend entries may be difficult to distinguish."
        )
    details = {
        "value": value,
        "plotted_samples": len(sample_order),
        "plotted_annotations": len(annotation_order),
        "plotted_segments": segment_count,
        "sample_order": sample_order,
        "condition_order": condition_order,
        "annotation_order": annotation_order,
        "annotation_colors": annotation_colors,
        "annotation_readability_warnings": annotation_readability_warnings,
        "stack_totals": [float(amount) for amount in bottom],
        "y_axis_label": y_axis_label,
        "input_cells": int(sample_totals.sum()),
        "condition_aggregation_performed": False,
    }
    if _return_diagnostics:
        return png, details
    return png


def run_sample_composition_plot(table: Any, *, value: str = "proportion") -> tuple[bytes, dict[str, Any]]:
    return _standalone_sample_composition_plot(table, value=value, _return_diagnostics=True)


def sample_composition_plot_code(*, value: str) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_sample_composition_plot)).strip()
    return f"""from __future__ import annotations

{implementation}


def plot_sample_composition(table):
    return _standalone_sample_composition_plot(table, value={value!r})
"""


__all__ = ["run_sample_composition_plot", "sample_composition_plot_code"]
