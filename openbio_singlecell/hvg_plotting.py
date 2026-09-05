from __future__ import annotations

import inspect
from textwrap import dedent


def render_hvg_selection_plot(adata):
    import io
    from collections.abc import Mapping

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from pandas.api.types import is_bool_dtype

    flavor_columns = {
        "seurat": ("dispersions", "dispersions_norm"),
        "cell_ranger": ("dispersions", "dispersions_norm"),
        "seurat_v3": ("variances", "variances_norm"),
        "seurat_v3_paper": ("variances", "variances_norm"),
        "pearson_residuals": ("variances", "residual_variances"),
    }
    hvg_metadata = adata.uns.get("hvg")
    flavor = hvg_metadata.get("flavor") if isinstance(hvg_metadata, Mapping) else None
    if flavor not in flavor_columns:
        raise ValueError(
            "HVG Selection Plot requires stored flavor metadata in adata.uns['hvg']['flavor'] "
            f"for one of {list(flavor_columns)!r}."
        )
    openbio_metadata = adata.uns.get("openbio_singlecell")
    history = openbio_metadata.get("analysis_history") if isinstance(openbio_metadata, Mapping) else None
    hvg_history = (
        [
            entry
            for entry in history.values()
            if isinstance(entry, Mapping) and entry.get("operation") == "highly_variable_genes"
        ]
        if isinstance(history, Mapping)
        else []
    )
    recorded_input_genes = int(hvg_history[-1]["input_genes"]) if hvg_history else None
    feature_axis_validation = "verified" if recorded_input_genes is not None else "unverified"
    base_variability_column, selection_variability_column = flavor_columns[flavor]
    origin_columns = ("highly_variable_algorithm", "highly_variable_forced")
    has_selection_origin = any(column in adata.var for column in origin_columns)
    required_columns = (
        "means",
        base_variability_column,
        selection_variability_column,
        "highly_variable",
        *(origin_columns if has_selection_origin else ()),
    )
    missing_columns = [column for column in required_columns if column not in adata.var]
    if missing_columns:
        raise ValueError(
            f"HVG Selection Plot requires complete stored evidence; missing adata.var columns {missing_columns!r}."
        )
    metric_values = {}
    for column in ("means", base_variability_column, selection_variability_column):
        try:
            values = np.asarray(adata.var[column], dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"HVG Selection Plot requires finite numeric evidence in adata.var[{column!r}]."
            ) from error
        if bool(np.isinf(values).any()) or (
            (column == "means" or flavor not in {"seurat", "cell_ranger"}) and bool(np.isnan(values).any())
        ):
            raise ValueError(f"HVG Selection Plot requires finite numeric evidence in adata.var[{column!r}].")
        metric_values[column] = values
    means = metric_values["means"]
    selection_values = {}
    for column in ("highly_variable", *(origin_columns if has_selection_origin else ())):
        series = adata.var[column]
        if not is_bool_dtype(series.dtype) or bool(series.isna().any()):
            raise ValueError(f"HVG Selection Plot requires complete boolean evidence in adata.var[{column!r}].")
        selection_values[column] = np.asarray(series, dtype=bool)
    final = selection_values["highly_variable"]
    unselected = ~final
    if has_selection_origin:
        algorithm = selection_values["highly_variable_algorithm"]
        forced_mask = selection_values["highly_variable_forced"]
        if not np.array_equal(final, algorithm | forced_mask):
            raise ValueError(
                "HVG Selection Plot selection masks are inconsistent: highly_variable must equal "
                "highly_variable_algorithm OR highly_variable_forced."
            )
        forced = forced_mask & ~algorithm
        if not bool(np.isfinite(metric_values[selection_variability_column][algorithm]).all()):
            raise ValueError(
                "HVG Selection Plot algorithm-selected features require finite selection evidence in "
                f"adata.var[{selection_variability_column!r}]."
            )
        categories = (
            (unselected, "Unselected", "#BDBDBD"),
            (algorithm, "Algorithm-selected", "#0072B2"),
            (forced, "Forced", "#D55E00"),
        )
        selection_counts = {
            "algorithm_selected": int(algorithm.sum()),
            "forced": int(forced.sum()),
            "unselected": int(unselected.sum()),
        }
    else:
        categories = ((unselected, "Unselected", "#BDBDBD"), (final, "Selected", "#0072B2"))
        selection_counts = {
            "selected": int(final.sum()),
            "algorithm_selected": None,
            "forced": None,
            "unselected": int(unselected.sum()),
        }
    figure = Figure(figsize=(10.0, 4.2))
    FigureCanvasAgg(figure)
    axes = figure.subplots(1, 2)
    plot_warnings = (
        []
        if recorded_input_genes is not None
        else [
            "HVG Selection Plot feature-axis completeness is unverified because no OpenBio "
            "highly_variable_genes history is available."
        ]
    )
    if recorded_input_genes is not None and recorded_input_genes != int(adata.n_vars):
        feature_axis_validation = "different_from_recorded"
        plot_warnings.append(
            f"HVG Selection Plot history recorded {recorded_input_genes:,} input genes, while the current "
            f"AnnData has {int(adata.n_vars):,}; only the current stored features are displayed."
        )
    if not has_selection_origin:
        plot_warnings.append(
            "HVG selection origin is unavailable; displayed the stored highly_variable mask without "
            "classifying features as algorithm-selected or forced."
        )
    plotted_points_by_panel = {}
    omitted_features_by_panel = {}
    for axis, column in zip(axes, (selection_variability_column, base_variability_column), strict=True):
        values = metric_values[column]
        finite = np.isfinite(values)
        if not bool(finite.any()):
            raise ValueError(f"HVG Selection Plot has no finite points for adata.var[{column!r}].")
        omitted = int((~finite).sum())
        plotted_points_by_panel[column] = int(finite.sum())
        omitted_features_by_panel[column] = omitted
        if omitted:
            plot_warnings.append(
                f"HVG Selection Plot {column!r} panel omitted {omitted:,} feature(s) with NaN variability."
            )
        for mask, label, color in categories:
            plotted = mask & finite
            axis.scatter(means[plotted], values[plotted], s=12, alpha=0.8, color=color, linewidths=0, label=label)
        axis.set_xlabel("Mean expression")
        axis.set_ylabel(column.replace("_", " ").capitalize())
        axis.grid(alpha=0.15, linewidth=0.6)
    axes[0].legend(frameon=False)
    figure.suptitle(f"HVG selection ({flavor})")
    figure.tight_layout()
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    details = {
        "flavor": flavor,
        "mean_column": "means",
        "variability_columns": [selection_variability_column, base_variability_column],
        "selection_counts": selection_counts,
        "input_features": int(adata.n_vars),
        "plotted_points_by_panel": plotted_points_by_panel,
        "omitted_features_by_panel": omitted_features_by_panel,
        "feature_axis_validation": feature_axis_validation,
    }
    details["warnings"] = plot_warnings
    return png, details


def hvg_selection_plot_code() -> str:
    renderer = inspect.getsource(render_hvg_selection_plot)
    return (
        renderer
        + "\n\n"
        + dedent(
            """
        def plot_hvg_selection(adata):
            return render_hvg_selection_plot(adata)[0]
        """
        )
    )


__all__ = ["hvg_selection_plot_code", "render_hvg_selection_plot"]
