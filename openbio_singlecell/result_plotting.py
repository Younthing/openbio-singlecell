from __future__ import annotations

import inspect
import threading
from collections.abc import Mapping
from typing import Any

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
MAX_MARKER_PLOT_GENES = 50
MAX_MARKER_PLOT_STAT_ROWS = 10_000
MAX_MARKER_PLOT_BYTES = 1024**3
MARKER_PLOT_RNG_LOCK = threading.RLock()


def _plot_embedding_impl(
    adata,
    *,
    embedding_key="X_umap",
    x_dimension=1,
    y_dimension=2,
    color="leiden",
    color_source="obs",
    source_kind="X",
    layer_name=None,
    color_mode="auto",
    point_size=10.0,
    continuous_color_map="viridis",
    categorical_palette="tab20",
    sort_order=True,
    missing_color="lightgray",
    legend_policy="automatic",
    resolved_category_colors=None,
    resolved_normalization=None,
):
    """Validate and render two dimensions of one stored embedding without mutating AnnData."""
    import io
    import math
    from numbers import Real

    import matplotlib
    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import ListedColormap, to_hex, to_rgba
    from matplotlib.figure import Figure
    from scipy import sparse

    png_signature = b"\x89PNG\r\n\x1a\n"
    legend_limit = 20
    qualitative_palettes = {
        "Accent",
        "Dark2",
        "Paired",
        "Pastel1",
        "Pastel2",
        "Set1",
        "Set2",
        "Set3",
        "tab10",
        "tab20",
        "tab20b",
        "tab20c",
    }

    def scalar_missing(value, label):
        missing = pd.isna(value)
        if not isinstance(missing, (bool, np.bool_)):
            raise ValueError(f"{label} values must be scalar.")
        return bool(missing)

    def strict_string(value, label, *, allow_empty=False):
        if not isinstance(value, str):
            raise TypeError(f"{label} must be a string.")
        result = value.strip()
        if not result and not allow_empty:
            raise ValueError(f"{label} cannot be empty.")
        return result

    def categorical_plan(series):
        missing_mask = np.zeros(len(series), dtype=bool)
        rendered = np.empty(len(series), dtype=object)
        identities_by_label = {}
        values_by_identity = {}
        category_identities = []

        if isinstance(series.dtype, pd.CategoricalDtype):
            raw_categories = list(series.cat.categories)
            category_identities = [
                (type(value).__module__, type(value).__qualname__, repr(value)) for value in raw_categories
            ]
            values_by_identity.update(zip(category_identities, raw_categories, strict=True))
            codes = series.cat.codes.to_numpy(copy=True)
            missing_mask = codes < 0
            for row, code in enumerate(codes):
                if code >= 0:
                    identity = category_identities[int(code)]
                    rendered[row] = str(values_by_identity[identity])
                else:
                    rendered[row] = None
            ordered = bool(series.cat.ordered)
        else:
            codes = np.full(len(series), -1, dtype=int)
            identity_to_code = {}
            for row, value in enumerate(series.tolist()):
                if scalar_missing(value, "Embedding categorical color"):
                    missing_mask[row] = True
                    rendered[row] = None
                    continue
                identity = (type(value).__module__, type(value).__qualname__, repr(value))
                if identity not in identity_to_code:
                    identity_to_code[identity] = len(category_identities)
                    category_identities.append(identity)
                    values_by_identity[identity] = value
                codes[row] = identity_to_code[identity]
                rendered[row] = str(value)
            ordered = False

        labels = []
        for identity in category_identities:
            label = str(values_by_identity[identity])
            identities_by_label.setdefault(label, set()).add(identity)
            labels.append(label)
        collisions = sorted(label for label, identities in identities_by_label.items() if len(identities) > 1)
        if collisions:
            raise ValueError(
                "Embedding categorical color contains distinct values that collapse after string conversion: "
                f"{collisions!r}."
            )
        if len(set(labels)) != len(labels):
            raise ValueError("Embedding categorical color labels must be unique after rendering.")
        if bool(missing_mask.any()) and "(missing)" in labels:
            raise ValueError("Embedding categorical color label '(missing)' collides with the missing-value label.")
        counts = {label: int(np.count_nonzero(rendered == label)) for label in labels}
        return rendered, missing_mask, labels, counts, ordered

    if not hasattr(adata, "n_obs") or not hasattr(adata, "obs") or not hasattr(adata, "obsm"):
        raise TypeError("Embedding Plot requires an in-memory AnnData object.")
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError("Embedding Plot does not support backed AnnData; load it into memory first.")
    n_obs = int(adata.n_obs)
    if n_obs < 1:
        raise ValueError("Embedding Plot requires at least one observation.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError("Embedding Plot requires unique observation identifiers.")

    embedding_key = strict_string(embedding_key, "Embedding embedding_key")
    color = strict_string(color, "Embedding color", allow_empty=True)
    color_mode = strict_string(color_mode, "Embedding color_mode")
    if color_mode not in {"auto", "categorical", "continuous"}:
        raise ValueError(f"Unsupported Embedding color_mode: {color_mode!r}.")
    continuous_color_map = strict_string(continuous_color_map, "Embedding continuous_color_map")
    categorical_palette = strict_string(categorical_palette, "Embedding categorical_palette")
    legend_policy = strict_string(legend_policy, "Embedding legend_policy")
    if legend_policy not in {"automatic", "show", "hide"}:
        raise ValueError(f"Unsupported Embedding legend_policy: {legend_policy!r}.")
    if not isinstance(sort_order, (bool, np.bool_)):
        raise TypeError("Embedding sort_order must be a boolean.")
    sort_order = bool(sort_order)
    if isinstance(point_size, (bool, np.bool_)) or not isinstance(point_size, Real):
        raise TypeError("Embedding point_size must be a finite positive number.")
    point_size = float(point_size)
    if not math.isfinite(point_size) or point_size <= 0:
        raise ValueError("Embedding point_size must be a finite positive number.")
    try:
        missing_rgba = to_rgba(missing_color)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid Embedding missing_color: {missing_color!r}.") from error
    missing_hex = to_hex(missing_rgba, keep_alpha=True)

    if embedding_key not in adata.obsm:
        raise ValueError(f"Embedding coordinates not found in obsm[{embedding_key!r}].")
    stored = adata.obsm[embedding_key]
    if isinstance(stored, pd.DataFrame):
        if not stored.index.equals(adata.obs_names):
            raise ValueError("Embedding coordinate DataFrame index is not exactly aligned to adata.obs_names.")
        coordinate_source = "dataframe_index_verified"
        stored = stored.to_numpy(copy=False)
    else:
        coordinate_source = "anndata_obsm_row_contract"
        if not sparse.issparse(stored):
            stored = np.asarray(stored)
    stored_coordinates = stored
    if stored_coordinates.ndim != 2 or stored_coordinates.shape[0] != n_obs or stored_coordinates.shape[1] < 2:
        raise ValueError(
            f"Embedding coordinates must have shape exactly ({n_obs}, n_dimensions) with n_dimensions >= 2; "
            f"received {stored_coordinates.shape!r}."
        )
    if not np.issubdtype(stored_coordinates.dtype, np.number) or np.issubdtype(
        stored_coordinates.dtype, np.bool_
    ):
        raise TypeError("Embedding coordinates must contain real numeric values.")
    if np.issubdtype(stored_coordinates.dtype, np.complexfloating):
        raise TypeError("Embedding coordinates must contain real numeric values, not complex values.")
    available_dimensions = int(stored_coordinates.shape[1])
    for value, label in ((x_dimension, "x_dimension"), (y_dimension, "y_dimension")):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"Embedding {label} must be an integer.")
        if value < 1 or value > available_dimensions:
            raise ValueError(
                f"Embedding {label} must be between 1 and {available_dimensions}; received {value!r}."
            )
    x_dimension, y_dimension = int(x_dimension), int(y_dimension)
    if x_dimension == y_dimension:
        raise ValueError("Embedding x_dimension and y_dimension must be different.")
    coordinates = stored_coordinates[:, [x_dimension - 1, y_dimension - 1]]
    if sparse.issparse(coordinates):
        coordinates = coordinates.toarray()
    coordinates = np.asarray(coordinates, dtype=float)
    if not bool(np.isfinite(coordinates).all()):
        raise ValueError("The selected Embedding coordinate dimensions contain non-finite values.")

    if color_source not in {"none", "obs", "gene"}:
        raise ValueError(f"Unsupported Embedding color source: {color_source!r}.")
    if color_source == "gene":
        color = strict_string(color, "Embedding gene")
    elif color_source == "none":
        color = ""

    warnings = []
    actual_mode = "none"
    color_details = {
        "target_kind": color_source,
        "target": color or None,
        "requested_mode": color_mode,
        "effective_mode": "none",
        "column": color if color_source == "obs" and color else None,
        "gene": color if color_source == "gene" and color else None,
        "missing_count": 0,
        "fixed_color": "#246bfe" if not color else None,
    }
    category_payload = None
    continuous_payload = None
    if color:
        if color_source == "obs":
            if color not in adata.obs:
                raise ValueError(f"Embedding color column not found in obs: {color!r}.")
            values = adata.obs[color]
            if len(values) != n_obs or not values.index.equals(adata.obs_names):
                raise ValueError("Embedding color column is not exactly aligned to adata.obs_names.")
        else:
            if source_kind == "raw":
                if adata.raw is None:
                    raise ValueError("Embedding gene-color source 'raw' was selected, but adata.raw is unavailable.")
                matrix = adata.raw.X
                var_names = adata.raw.var_names
                layer_name = None
            elif source_kind == "layer":
                layer_name = strict_string(layer_name, "Embedding gene-color layer")
                if layer_name not in adata.layers:
                    raise ValueError(f"Embedding gene-color layer not found: {layer_name!r}.")
                matrix = adata.layers[layer_name]
                var_names = adata.var_names
            elif source_kind == "X":
                if layer_name is not None:
                    raise ValueError("Embedding gene-color X source cannot carry a layer name.")
                matrix = adata.X
                var_names = adata.var_names
            else:
                raise ValueError(f"Unsupported Embedding gene-color source: {source_kind!r}.")
            if not bool(var_names.is_unique):
                raise ValueError("Embedding gene-color source requires unique feature identifiers.")
            expected_shape = (n_obs, len(var_names))
            if tuple(matrix.shape) != expected_shape:
                raise ValueError(
                    f"Embedding gene-color source must have shape {expected_shape!r}; "
                    f"received {tuple(matrix.shape)!r}."
                )
            if color not in var_names:
                raise ValueError(f"Embedding gene not found in the selected expression source: {color!r}.")
            index = var_names.get_loc(color)
            if not isinstance(index, (int, np.integer)):
                raise ValueError("Embedding gene lookup was ambiguous.")
            selected = matrix[:, index]
            numeric = selected.toarray().ravel() if sparse.issparse(selected) else np.asarray(selected).ravel()
            if not np.issubdtype(numeric.dtype, np.number) or np.issubdtype(numeric.dtype, np.bool_):
                raise TypeError("Embedding gene-color source must contain real numeric values.")
            if np.iscomplexobj(numeric):
                raise TypeError("Embedding gene-color source must contain real numeric values, not complex values.")
            if not bool(np.isfinite(numeric).all()):
                raise ValueError("Embedding selected gene contains non-finite expression values.")
            values = pd.Series(numeric, index=adata.obs_names, name=color)
            color_details["expression_source"] = {
                "source": source_kind,
                **({"layer_name": layer_name} if source_kind == "layer" else {}),
            }
        dtype = values.dtype
        if color_source == "gene":
            actual_mode = "continuous"
        elif color_mode == "auto":
            if (
                isinstance(dtype, pd.CategoricalDtype)
                or pd.api.types.is_bool_dtype(dtype)
                or pd.api.types.is_string_dtype(dtype)
                or pd.api.types.is_object_dtype(dtype)
            ):
                actual_mode = "categorical"
            elif pd.api.types.is_numeric_dtype(dtype) and not pd.api.types.is_complex_dtype(dtype):
                actual_mode = "continuous"
            else:
                raise TypeError(
                    f"Embedding color column {color!r} has unsupported dtype {dtype!r}; choose an explicit valid mode."
                )
        else:
            actual_mode = color_mode

        if actual_mode == "categorical":
            rendered, missing_mask, categories, counts, declared_ordered = categorical_plan(values)
            if categorical_palette not in matplotlib.colormaps:
                raise ValueError(f"Unknown Matplotlib categorical palette: {categorical_palette!r}.")
            palette = matplotlib.colormaps[categorical_palette]
            if resolved_category_colors is None:
                if isinstance(palette, ListedColormap) and len(categories) <= len(palette.colors):
                    available = list(palette.colors)
                    colors = [to_hex(value, keep_alpha=True) for value in available[: len(categories)]]
                else:
                    positions = np.linspace(0.0, 1.0, max(len(categories), 1))
                    colors = [to_hex(palette(float(position)), keep_alpha=True) for position in positions]
                category_colors = dict(zip(categories, colors, strict=True))
            else:
                if not isinstance(resolved_category_colors, dict):
                    raise TypeError("Resolved Embedding category colors must be a mapping.")
                if set(resolved_category_colors) != set(categories) or len(resolved_category_colors) != len(categories):
                    raise ValueError("Resolved Embedding category colors do not exactly match the observed categories.")
                category_colors = {}
                for category in categories:
                    try:
                        category_colors[category] = to_hex(
                            to_rgba(resolved_category_colors[category]), keep_alpha=True
                        )
                    except (TypeError, ValueError) as error:
                        raise ValueError(f"Invalid resolved color for Embedding category {category!r}.") from error
            if len(set(category_colors.values())) != len(category_colors):
                warnings.append("The selected categorical palette assigns duplicate rendered colors to different labels.")
            if bool(missing_mask.any()) and missing_hex in set(category_colors.values()):
                warnings.append("The selected missing_color duplicates a categorical palette color.")
            if categorical_palette not in qualitative_palettes:
                warnings.append(
                    f"Palette {categorical_palette!r} is not a recognized qualitative Matplotlib palette; "
                    "category colors should not be interpreted as magnitudes."
                )
            if bool(missing_mask.any()):
                warnings.append(
                    f"Rendered {int(missing_mask.sum()):,} observations with missing categorical color values "
                    f"using the explicit special color {missing_hex}."
                )
            legend_entries = len(categories) + int(bool(missing_mask.any()))
            legend_shown = (
                legend_policy == "show" or (legend_policy == "automatic" and legend_entries <= legend_limit)
            )
            if legend_policy == "automatic" and not legend_shown:
                warnings.append(
                    f"Legend omitted because {legend_entries} entries exceed the automatic limit of {legend_limit}; "
                    "the category-color mapping remains available in summary."
                )
            category_payload = (rendered, missing_mask, categories, category_colors, legend_shown)
            color_details.update(
                {
                    "effective_mode": "categorical",
                    "dtype": str(dtype),
                    "category_order": categories,
                    "declared_categorical_ordered": declared_ordered,
                    "category_counts": counts,
                    "category_colors": category_colors,
                    "palette": categorical_palette,
                    "missing_count": int(missing_mask.sum()),
                    "missing_color": missing_hex,
                    "legend_policy": legend_policy,
                    "legend_entry_limit": legend_limit,
                    "legend_shown": legend_shown,
                    "draw_order": "declared category order; missing values last",
                }
            )
        else:
            if not pd.api.types.is_numeric_dtype(dtype) or pd.api.types.is_bool_dtype(dtype):
                raise TypeError(f"Embedding continuous color target {color!r} must have a non-boolean numeric dtype.")
            if pd.api.types.is_complex_dtype(dtype):
                raise TypeError(f"Embedding continuous color target {color!r} must contain real values.")
            try:
                numeric = values.to_numpy(dtype=float, na_value=np.nan, copy=True)
            except (TypeError, ValueError) as error:
                raise TypeError(f"Embedding continuous color target {color!r} must contain real numeric values.") from error
            finite_mask = np.isfinite(numeric)
            missing_mask = pd.isna(numeric)
            infinity_mask = ~finite_mask & ~missing_mask
            if bool(infinity_mask.any()):
                warnings.append(
                    f"Rendered {int(infinity_mask.sum()):,} observations with infinite continuous color values "
                    f"using the explicit special color {missing_hex}; stored values were left unchanged."
                )
            if not bool(finite_mask.any()):
                raise ValueError(f"Embedding continuous color target {color!r} contains no finite values.")
            if bool(missing_mask.any()):
                warnings.append(
                    f"Rendered {int(missing_mask.sum()):,} observations with missing continuous color values "
                    f"using the explicit special color {missing_hex}."
                )
            if continuous_color_map not in matplotlib.colormaps:
                raise ValueError(f"Unknown Matplotlib continuous colormap: {continuous_color_map!r}.")
            finite_values = numeric[finite_mask]
            observed_min = float(finite_values.min())
            observed_max = float(finite_values.max())
            if resolved_normalization is None:
                if observed_min == observed_max:
                    padding = max(abs(observed_min) * 1e-6, 1e-12)
                    vmin, vmax = observed_min - padding, observed_max + padding
                    warnings.append(
                        f"Continuous color column {color!r} is constant at {observed_min:.6g}; a symmetric "
                        "display-only normalization interval was used."
                    )
                else:
                    vmin, vmax = observed_min, observed_max
            else:
                if not isinstance(resolved_normalization, (tuple, list)) or len(resolved_normalization) != 2:
                    raise TypeError("Resolved Embedding normalization must contain exactly vmin and vmax.")
                vmin, vmax = (float(value) for value in resolved_normalization)
                if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin >= vmax:
                    raise ValueError("Resolved Embedding normalization must be finite with vmin < vmax.")
            finite_indices = np.flatnonzero(finite_mask)
            if sort_order:
                finite_indices = finite_indices[np.argsort(numeric[finite_indices], kind="stable")]
            quantiles = np.quantile(finite_values, [0.0, 0.25, 0.5, 0.75, 1.0])
            continuous_payload = (numeric, ~finite_mask, finite_indices, float(vmin), float(vmax))
            color_details.update(
                {
                    "effective_mode": "continuous",
                    "dtype": str(dtype),
                    "missing_count": int(missing_mask.sum()),
                    "infinite_count": int(infinity_mask.sum()),
                    "masked_nonfinite_count": int((~finite_mask).sum()),
                    "missing_color": missing_hex,
                    "finite_count": int(finite_mask.sum()),
                    "observed_range": [observed_min, observed_max],
                    "observed_quantiles": {
                        "min": float(quantiles[0]),
                        "q1": float(quantiles[1]),
                        "median": float(quantiles[2]),
                        "q3": float(quantiles[3]),
                        "max": float(quantiles[4]),
                    },
                    "normalization": {"kind": "linear", "vmin": float(vmin), "vmax": float(vmax)},
                    "colormap": continuous_color_map,
                    "sort_order": sort_order,
                    "draw_order": (
                        "non-finite color values first, then finite values in stable ascending order"
                        if sort_order
                        else "non-finite color values first, then finite values in observation order"
                    ),
                }
            )

    figure = Figure(figsize=(8.0, 7.0), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    try:
        if actual_mode == "none":
            axis.scatter(
                coordinates[:, 0],
                coordinates[:, 1],
                s=point_size,
                alpha=0.8,
                linewidths=0,
                color="#246bfe",
            )
        elif actual_mode == "categorical":
            rendered, missing_mask, categories, category_colors, legend_shown = category_payload
            for category in categories:
                mask = rendered == category
                if bool(mask.any()):
                    axis.scatter(
                        coordinates[mask, 0],
                        coordinates[mask, 1],
                        s=point_size,
                        alpha=0.8,
                        linewidths=0,
                        color=category_colors[category],
                        label=category,
                    )
            if bool(missing_mask.any()):
                axis.scatter(
                    coordinates[missing_mask, 0],
                    coordinates[missing_mask, 1],
                    s=point_size,
                    alpha=0.8,
                    linewidths=0,
                    color=missing_hex,
                    label="(missing)",
                )
            if legend_shown:
                axis.legend(title=color, bbox_to_anchor=(1.02, 1.0), loc="upper left", markerscale=1.5)
        else:
            numeric, missing_mask, finite_indices, vmin, vmax = continuous_payload
            if bool(missing_mask.any()):
                axis.scatter(
                    coordinates[missing_mask, 0],
                    coordinates[missing_mask, 1],
                    s=point_size,
                    alpha=0.8,
                    linewidths=0,
                    color=missing_hex,
                )
            points = axis.scatter(
                coordinates[finite_indices, 0],
                coordinates[finite_indices, 1],
                c=numeric[finite_indices],
                cmap=continuous_color_map,
                vmin=vmin,
                vmax=vmax,
                s=point_size,
                alpha=0.8,
                linewidths=0,
            )
            figure.colorbar(points, ax=axis, label=color)
        embedding_title = "UMAP" if embedding_key == "X_umap" else embedding_key
        x_label = f"{embedding_title} {x_dimension}"
        y_label = f"{embedding_title} {y_dimension}"
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.set_title(f"{embedding_title} colored by {color}" if color else embedding_title)
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
        png = buffer.getvalue()
    finally:
        figure.clear()
    if not png.startswith(png_signature):
        raise RuntimeError("Embedding renderer did not produce a valid PNG payload.")

    details = {
        "embedding_key": embedding_key,
        "n_observations": n_obs,
        "plotted_observations": n_obs,
        "coordinate_alignment": coordinate_source,
        "available_coordinate_dimensions": available_dimensions,
        "used_coordinate_dimensions": [x_dimension, y_dimension],
        "coordinate_ranges": {
            x_label: [float(coordinates[:, 0].min()), float(coordinates[:, 0].max())],
            y_label: [float(coordinates[:, 1].min()), float(coordinates[:, 1].max())],
        },
        "color": color_details,
        "rendering": {
            "figure_size_inches": [8.0, 7.0],
            "dpi": 120,
            "bbox_inches": "tight",
            "point_size_points_squared": point_size,
            "alpha": 0.8,
            "linewidth_points": 0.0,
            "uncolored_point_color": "#246bfe",
        },
        "warnings": warnings,
    }
    return png, details


def embedding_plot_code(*, parameters: Mapping[str, Any], details: Mapping[str, Any]) -> str:
    color_details = details["color"]
    category_colors = color_details.get("category_colors")
    normalization = color_details.get("normalization")
    resolved_normalization = None
    if isinstance(normalization, Mapping):
        resolved_normalization = (normalization["vmin"], normalization["vmax"])
    implementation = inspect.getsource(_plot_embedding_impl)
    arguments = ",\n        ".join(
        [
            f"embedding_key={parameters['embedding_key']!r}",
            f"x_dimension={parameters['x_dimension']!r}",
            f"y_dimension={parameters['y_dimension']!r}",
            f"color={color_details['target'] or ''!r}",
            f"color_source={color_details['target_kind']!r}",
            f"source_kind={color_details.get('expression_source', {}).get('source', 'X')!r}",
            f"layer_name={color_details.get('expression_source', {}).get('layer_name')!r}",
            f"color_mode={color_details['requested_mode']!r}",
            f"point_size={parameters['point_size']!r}",
            f"continuous_color_map={parameters['continuous_color_map']!r}",
            f"categorical_palette={parameters['categorical_palette']!r}",
            f"sort_order={parameters['sort_order']!r}",
            f"missing_color={parameters['missing_color']!r}",
            f"legend_policy={parameters['legend_policy']!r}",
            f"resolved_category_colors={category_colors!r}",
            f"resolved_normalization={resolved_normalization!r}",
        ]
    )
    return (
        f"{implementation}\n\n"
        "def plot_embedding(adata):\n"
        "    png, _details = _plot_embedding_impl(\n"
        "        adata,\n"
        f"        {arguments},\n"
        "    )\n"
        "    return png\n"
    )


def embedding_provenance(adata, embedding_key: str) -> tuple[dict[str, Any], list[str]]:
    draw_graph_prefix = "X_draw_graph_"
    draw_graph_suffix = embedding_key.removeprefix(draw_graph_prefix) if embedding_key.startswith(draw_graph_prefix) else None
    metadata_key = (
        "draw_graph"
        if draw_graph_suffix
        else {"X_umap": "umap", "X_pca": "pca"}.get(embedding_key, embedding_key)
    )
    provenance: dict[str, Any] = {
        "embedding_key": embedding_key,
        "metadata_key": metadata_key,
        "scanpy_metadata_verified": False,
        "openbio_history_verified": False,
    }
    warnings: list[str] = []
    metadata = adata.uns.get(metadata_key)
    metadata_params = metadata.get("params") if isinstance(metadata, Mapping) else None
    if not isinstance(metadata_params, Mapping):
        metadata_params = None

    openbio = adata.uns.get("openbio_singlecell")
    history = openbio.get("analysis_history") if isinstance(openbio, Mapping) else None
    latest_draw_graph_parameters = None
    if isinstance(history, Mapping):
        for entry in reversed(list(history.values())):
            if not isinstance(entry, Mapping):
                continue
            parameters = entry.get("parameters")
            if not isinstance(parameters, Mapping):
                continue
            operation = entry.get("operation")
            if operation == "force_directed_graph" and latest_draw_graph_parameters is None:
                latest_draw_graph_parameters = parameters
            matches = operation == "umap" and parameters.get("key_added") == embedding_key
            if operation == "force_directed_graph" and draw_graph_suffix:
                layout = parameters.get("layout")
                matches = (
                    parameters.get("key_suffix") == draw_graph_suffix
                    and isinstance(layout, str)
                )
            if not matches:
                continue
            provenance["openbio_history_verified"] = True
            provenance["openbio_operation"] = operation
            provenance["openbio_parameters"] = {
                str(key): value
                for key, value in parameters.items()
                if value is None or isinstance(value, (str, bool, int, float))
            }
            if isinstance(entry.get("random_seed"), int):
                provenance["embedding_random_seed"] = int(entry["random_seed"])
            break
    if metadata_params is not None:
        if draw_graph_suffix is None:
            metadata_matches = True
        elif latest_draw_graph_parameters is None:
            draw_graph_keys = {key for key in adata.obsm if key.startswith(draw_graph_prefix)}
            metadata_matches = draw_graph_keys == {embedding_key} and isinstance(metadata_params.get("layout"), str)
        else:
            metadata_matches = (
                provenance.get("openbio_operation") == "force_directed_graph"
                and latest_draw_graph_parameters.get("key_suffix") == draw_graph_suffix
                and metadata_params.get("layout") == provenance["openbio_parameters"].get("layout")
                and metadata_params.get("layout") == latest_draw_graph_parameters.get("layout")
            )
        if metadata_matches:
            provenance["scanpy_metadata_verified"] = True
            provenance["scanpy_parameters"] = {
                str(key): value
                for key, value in metadata_params.items()
                if value is None or isinstance(value, (str, bool, int, float))
            }
    if not provenance["scanpy_metadata_verified"] and not provenance["openbio_history_verified"]:
        warnings.append(
            f"No matching stored parameters or OpenBio history were found for {embedding_key!r}; "
            "the coordinates were validated but their embedding computation provenance is unverified."
        )
    return provenance, warnings


def _marker_expression_plot_impl(
    adata,
    *,
    genes,
    groupby="leiden",
    plot=None,
    source_kind="layer",
    layer_name="log1p_norm",
    group_order="observed",
    random_seed=0,
    expected_group_order=None,
    _science=None,
    _rng_lock=None,
):
    """Validate, summarize, and render a marker panel from one explicit expression representation."""
    import contextlib
    import io
    import math
    from collections.abc import Mapping
    from numbers import Real
    from types import SimpleNamespace

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    if _science is None:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from scipy import sparse

        _science = SimpleNamespace(ad=ad, np=np, pd=pd, sc=sc, sparse=sparse)
    ad = _science.ad
    np = _science.np
    pd = _science.pd
    sc = _science.sc
    sparse = _science.sparse

    png_signature = b"\x89PNG\r\n\x1a\n"
    max_stat_rows = 10_000
    max_plot_bytes = 1024**3

    def strict_string(value, label):
        if not isinstance(value, str):
            raise TypeError(f"{label} must be a string.")
        result = value.strip()
        if not result:
            raise ValueError(f"{label} cannot be empty.")
        return result

    def scalar_missing(value, label):
        missing = pd.isna(value)
        if not isinstance(missing, (bool, np.bool_)):
            raise ValueError(f"{label} values must be scalar.")
        return bool(missing)

    def parse_plot(value):
        if value is None:
            value = {
                "plot": "dotplot",
                "standard_scale": "var",
                "expression_cutoff": 0.0,
                "mean_only_expressed": False,
            }
        if not isinstance(value, Mapping):
            raise TypeError("Marker Expression plot must be a DynamicCombo value.")
        variant = value.get("plot")
        if not isinstance(variant, str):
            raise TypeError("Marker Expression plot selection must be a string.")
        variant = variant.strip()
        schemas = {
            "dotplot": {"standard_scale", "expression_cutoff", "mean_only_expressed"},
            "matrixplot": {"standard_scale"},
            "tracksplot": set(),
            "violin": {"density_norm", "show_cells", "y_scale"},
        }
        if variant not in schemas:
            raise ValueError(f"Unsupported Marker Expression plot: {variant!r}.")
        allowed = {"plot", *schemas[variant]}
        unknown = sorted(str(key) for key in value if key not in allowed)
        if unknown:
            raise ValueError(f"Marker Expression {variant} received inactive or unknown parameters: {unknown!r}.")
        effective = {"plot": variant}
        if variant in {"dotplot", "matrixplot"}:
            standard_scale = value.get("standard_scale", "var")
            if not isinstance(standard_scale, str) or standard_scale.strip() not in {"none", "var", "group"}:
                raise ValueError(f"Unsupported Marker Expression standard_scale: {standard_scale!r}.")
            effective["standard_scale"] = standard_scale.strip()
        if variant == "dotplot":
            cutoff = value.get("expression_cutoff", 0.0)
            if isinstance(cutoff, (bool, np.bool_)) or not isinstance(cutoff, Real):
                raise TypeError("Marker Expression expression_cutoff must be a finite number.")
            cutoff = float(cutoff)
            if not math.isfinite(cutoff):
                raise ValueError("Marker Expression expression_cutoff must be finite.")
            mean_only = value.get("mean_only_expressed", False)
            if not isinstance(mean_only, (bool, np.bool_)):
                raise TypeError("Marker Expression mean_only_expressed must be a boolean.")
            effective.update(expression_cutoff=cutoff, mean_only_expressed=bool(mean_only))
        elif variant == "violin":
            density_norm = value.get("density_norm", "width")
            if not isinstance(density_norm, str) or density_norm.strip() not in {"width", "area", "count"}:
                raise ValueError(f"Unsupported Marker Expression density_norm: {density_norm!r}.")
            show_cells = value.get("show_cells", False)
            if not isinstance(show_cells, (bool, np.bool_)):
                raise TypeError("Marker Expression show_cells must be a boolean.")
            y_scale = value.get("y_scale", "linear")
            if not isinstance(y_scale, str) or y_scale.strip() not in {"linear", "log"}:
                raise ValueError(f"Unsupported Marker Expression y_scale: {y_scale!r}.")
            effective.update(
                density_norm=density_norm.strip(), show_cells=bool(show_cells), y_scale=y_scale.strip()
            )
        return effective

    def render_groups(series):
        if isinstance(series.dtype, pd.CategoricalDtype):
            raw_categories = list(series.cat.categories)
            codes = series.cat.codes.to_numpy(copy=True)
            if bool((codes < 0).any()):
                raise ValueError("Marker Expression groupby contains missing values.")
            unused = [raw_categories[index] for index in range(len(raw_categories)) if not bool((codes == index).any())]
            if unused:
                raise ValueError(f"Marker Expression groupby contains unused declared categories: {unused!r}.")
            values = raw_categories
            row_identities = [
                (type(raw_categories[int(code)]).__module__, type(raw_categories[int(code)]).__qualname__, repr(raw_categories[int(code)]))
                for code in codes
            ]
        else:
            dtype = series.dtype
            if pd.api.types.is_numeric_dtype(dtype) and not pd.api.types.is_bool_dtype(dtype):
                raise TypeError("Marker Expression groupby must be categorical; numeric grouping requires explicit upstream binning.")
            if not (
                pd.api.types.is_bool_dtype(dtype)
                or pd.api.types.is_string_dtype(dtype)
                or pd.api.types.is_object_dtype(dtype)
            ):
                raise TypeError(f"Marker Expression groupby has unsupported dtype {dtype!r}.")
            identities = []
            value_by_identity = {}
            row_identities = []
            for value in series.tolist():
                if scalar_missing(value, "Marker Expression groupby"):
                    raise ValueError("Marker Expression groupby contains missing values.")
                identity = (type(value).__module__, type(value).__qualname__, repr(value))
                if identity not in value_by_identity:
                    value_by_identity[identity] = value
                    identities.append(identity)
                row_identities.append(identity)
            values = [value_by_identity[identity] for identity in identities]

        identity_to_label = {}
        identities_by_label = {}
        for value in values:
            identity = (type(value).__module__, type(value).__qualname__, repr(value))
            label = str(value).strip()
            if not label:
                raise ValueError("Marker Expression groupby contains blank labels.")
            identity_to_label[identity] = label
            identities_by_label.setdefault(label, set()).add(identity)
        collisions = sorted(label for label, identities in identities_by_label.items() if len(identities) > 1)
        if collisions:
            raise ValueError(
                "Marker Expression groupby contains distinct values that collapse after string conversion: "
                f"{collisions!r}."
            )
        labels = [identity_to_label[(type(value).__module__, type(value).__qualname__, repr(value))] for value in values]
        rendered = np.asarray([identity_to_label[identity] for identity in row_identities], dtype=object)
        counts = {label: int(np.count_nonzero(rendered == label)) for label in labels}
        return rendered, labels, counts

    def figure_from_result(value):
        visited = set()

        def visit(item):
            if item is None or id(item) in visited:
                return None
            visited.add(id(item))
            candidate = getattr(item, "fig", None)
            if candidate is not None and hasattr(candidate, "savefig"):
                return candidate
            candidate = getattr(item, "figure", None)
            if candidate is not None and hasattr(candidate, "savefig"):
                return candidate
            if isinstance(item, Mapping):
                for nested in item.values():
                    found = visit(nested)
                    if found is not None:
                        return found
            elif isinstance(item, (list, tuple)):
                for nested in item:
                    found = visit(nested)
                    if found is not None:
                        return found
            elif isinstance(item, np.ndarray):
                for nested in item.ravel().tolist():
                    found = visit(nested)
                    if found is not None:
                        return found
            return None

        result = visit(value)
        if result is None:
            raise RuntimeError("Scanpy marker-expression backend did not return a Matplotlib figure.")
        return result

    @contextlib.contextmanager
    def scanpy_v1_context():
        settings = getattr(sc, "settings", None)
        if settings is None:
            raise RuntimeError("Scanpy marker-expression backend does not expose settings.")
        override = getattr(settings, "override", None)
        preset = getattr(sc, "Preset", None)
        if callable(override) and preset is not None and hasattr(preset, "ScanpyV1"):
            with override(preset=preset.ScanpyV1, autoshow=False):
                yield
            return
        if not hasattr(settings, "autoshow"):
            raise RuntimeError("Scanpy marker-expression backend cannot locally control autoshow.")
        original_autoshow = settings.autoshow
        settings.autoshow = False
        try:
            yield
        finally:
            settings.autoshow = original_autoshow

    @contextlib.contextmanager
    def isolated_numpy_random(seed):
        state = np.random.get_state()
        np.random.seed(seed)
        try:
            yield
        finally:
            np.random.set_state(state)

    if not hasattr(adata, "n_obs") or not hasattr(adata, "obs") or not hasattr(adata, "var_names"):
        raise TypeError("Marker Expression Plot requires an in-memory AnnData object.")
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError("Marker Expression Plot does not support backed AnnData; load it into memory first.")
    n_obs = int(adata.n_obs)
    if n_obs < 1:
        raise ValueError("Marker Expression Plot requires at least one observation.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError("Marker Expression Plot requires unique observation identifiers.")
    if not isinstance(genes, str):
        raise TypeError("Marker Expression genes must be a comma-separated string.")
    gene_names = list(dict.fromkeys(item.strip() for item in genes.split(",") if item.strip()))
    if not gene_names:
        raise ValueError("Marker Expression Plot requires at least one comma-separated gene.")
    groupby = strict_string(groupby, "Marker Expression groupby")
    group_order = strict_string(group_order, "Marker Expression group_order")
    if group_order not in {"observed", "dendrogram"}:
        raise ValueError(f"Unsupported Marker Expression group_order: {group_order!r}.")
    if isinstance(random_seed, (bool, np.bool_)) or not isinstance(random_seed, (int, np.integer)):
        raise TypeError("Marker Expression random_seed must be an integer.")
    random_seed = int(random_seed)
    if random_seed < 0 or random_seed > 2**31 - 1:
        raise ValueError("Marker Expression random_seed must be between 0 and 2^31 - 1.")
    mode = parse_plot(plot)

    if source_kind not in {"X", "raw", "layer"}:
        raise ValueError(f"Unsupported Marker Expression source: {source_kind!r}.")
    if source_kind == "raw":
        if adata.raw is None:
            raise ValueError("Marker Expression source 'raw' was selected, but adata.raw is unavailable.")
        matrix = adata.raw.X
        var_names = adata.raw.var_names
        source_layer_name = None
    elif source_kind == "layer":
        layer_name = strict_string(layer_name, "Marker Expression layer")
        if layer_name not in adata.layers:
            raise ValueError(f"Marker Expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
        var_names = adata.var_names
        source_layer_name = layer_name
    else:
        if layer_name is not None:
            raise ValueError("Marker Expression X source cannot carry a layer name.")
        matrix = adata.X
        var_names = adata.var_names
        source_layer_name = None
    if not bool(var_names.is_unique):
        raise ValueError("Marker Expression selected source requires unique feature identifiers.")
    expected_shape = (n_obs, len(var_names))
    if tuple(matrix.shape) != expected_shape:
        raise ValueError(
            f"Marker Expression selected source must have shape {expected_shape!r}; received {tuple(matrix.shape)!r}."
        )
    missing_genes = [gene for gene in gene_names if gene not in var_names]
    if missing_genes:
        raise ValueError(f"Marker Expression genes not found in the selected source: {missing_genes!r}.")
    indices = [var_names.get_loc(gene) for gene in gene_names]
    if any(not isinstance(index, (int, np.integer)) for index in indices):
        raise ValueError("Marker Expression selected source feature lookup was ambiguous.")
    selected = matrix[:, indices]
    values = selected.data if sparse.issparse(selected) else np.asarray(selected).ravel()
    values = np.asarray(values)
    if not np.issubdtype(values.dtype, np.number) or np.issubdtype(values.dtype, np.bool_):
        raise TypeError("Marker Expression selected source must contain real numeric values.")
    if np.iscomplexobj(values):
        raise TypeError("Marker Expression selected source must contain real numeric values, not complex values.")
    if values.size and not bool(np.isfinite(values).all()):
        raise ValueError("Marker Expression selected genes contain non-finite expression values.")
    if groupby not in adata.obs:
        raise ValueError(f"Marker Expression groupby column not found in obs: {groupby!r}.")
    groups, observed_order, group_counts = render_groups(adata.obs[groupby])
    statistic_rows = len(observed_order) * len(gene_names)
    if statistic_rows > max_stat_rows:
        raise ValueError(
            "Marker Expression group-by-gene summary would exceed "
            f"the {max_stat_rows:,}-row reporting limit."
        )
    width = float(min(20.0, max(7.0, 3.0 + 0.6 * len(gene_names))))
    height = float(min(16.0, max(5.0, 2.5 + 0.45 * len(observed_order))))
    selected_storage_bytes = (
        int(selected.data.nbytes + selected.indices.nbytes + selected.indptr.nbytes)
        if sparse.issparse(selected)
        else int(np.asarray(selected).nbytes)
    )
    panel_bytes = int(n_obs * len(gene_names) * 8)
    work_matrix_bytes = panel_bytes
    dendrogram_matrix_bytes = panel_bytes if group_order == "dendrogram" else 0
    group_aggregation_bytes = int(len(observed_order) * len(gene_names) * 8)
    correlation_bytes = int(len(observed_order) ** 2 * 8) if group_order == "dendrogram" else 0
    figure_rgba_bytes = int(math.ceil(width * 120) * math.ceil(height * 120) * 4)
    statistics_bytes = int(statistic_rows * 1024)
    group_label_bytes = int(n_obs * 16)
    estimated_peak_bytes = sum(
        (
            selected_storage_bytes,
            panel_bytes,
            work_matrix_bytes,
            dendrogram_matrix_bytes,
            group_aggregation_bytes,
            correlation_bytes,
            figure_rgba_bytes * 3,
            statistics_bytes,
            group_label_bytes,
        )
    )
    resource_estimate = {
        "limit_bytes": max_plot_bytes,
        "estimated_peak_bytes": estimated_peak_bytes,
        "selected_source_slice_bytes": selected_storage_bytes,
        "dense_panel_bytes": panel_bytes,
        "plotting_matrix_copy_bytes": work_matrix_bytes,
        "dendrogram_matrix_copy_bytes": dendrogram_matrix_bytes,
        "group_aggregation_bytes": group_aggregation_bytes,
        "correlation_matrix_bytes": correlation_bytes,
        "figure_and_png_working_bytes": figure_rgba_bytes * 3,
        "summary_rows_estimated_bytes": statistics_bytes,
        "group_labels_estimated_bytes": group_label_bytes,
    }
    if estimated_peak_bytes > max_plot_bytes:
        raise ValueError(
            f"Marker Expression estimated peak working set {estimated_peak_bytes:,} bytes exceeds the "
            f"{max_plot_bytes:,}-byte limit."
        )
    panel = selected.toarray() if sparse.issparse(selected) else np.asarray(selected).copy()
    panel = np.asarray(panel, dtype=float)
    if panel.shape != (n_obs, len(gene_names)) or not bool(np.isfinite(panel).all()):
        raise RuntimeError("Marker Expression selected panel could not be materialized as a finite matrix.")
    warnings = []
    if mode["plot"] == "violin" and mode["y_scale"] == "log" and bool((panel <= 0).any()):
        warnings.append(
            "The selected panel contains nonpositive expression values that cannot be displayed on the log axis; "
            "the original values were passed to Scanpy unchanged."
        )
    dendrogram_details = {
        "enabled": group_order == "dendrogram",
        "source_matches_panel": True,
        "correlation": None,
        "linkage": None,
        "optimal_ordering": None,
    }
    resolved_order = list(observed_order)
    if group_order == "dendrogram":
        if len(observed_order) < 2:
            raise ValueError("Marker Expression dendrogram group order requires at least two observed groups.")
        if len(gene_names) < 2:
            raise ValueError("Marker Expression dendrogram group order requires at least two selected genes.")
        group_means = np.vstack([panel[groups == label].mean(axis=0) for label in observed_order])
        if any(bool(np.allclose(row, row[0], rtol=0.0, atol=0.0)) for row in group_means):
            raise ValueError(
                "Marker Expression dendrogram is undefined because at least one group has a constant selected-gene profile."
            )
        correlations = np.corrcoef(group_means)
        if correlations.shape != (len(observed_order), len(observed_order)) or not bool(
            np.isfinite(correlations).all()
        ):
            raise ValueError("Marker Expression dendrogram produced non-finite group correlations.")
        private_groupby = "_openbio_marker_group"
        minimal = ad.AnnData(
            X=panel.copy(),
            obs=pd.DataFrame(
                {
                    private_groupby: pd.Categorical(
                        groups.copy(), categories=observed_order, ordered=True
                    )
                },
                index=adata.obs_names.copy(),
            ),
            var=pd.DataFrame(index=pd.Index(gene_names, dtype=object)),
        )
        dendrogram_result = sc.tl.dendrogram(
            minimal,
            private_groupby,
            var_names=list(gene_names),
            use_raw=False,
            cor_method="pearson",
            linkage_method="complete",
            optimal_ordering=False,
            key_added="_openbio_marker_panel_dendrogram",
            inplace=False,
        )
        if not isinstance(dendrogram_result, Mapping):
            raise RuntimeError("Scanpy dendrogram backend did not return a result mapping.")
        ordered = dendrogram_result.get("categories_ordered")
        if not isinstance(ordered, (list, tuple)):
            raise RuntimeError("Scanpy dendrogram backend did not return categories_ordered.")
        resolved_order = [str(value) for value in ordered]
        if len(resolved_order) != len(observed_order) or set(resolved_order) != set(observed_order):
            raise RuntimeError("Scanpy dendrogram category order does not match the observed groups.")
        if expected_group_order is not None and list(expected_group_order) != resolved_order:
            raise RuntimeError("Recomputed Marker Expression dendrogram order differs from the expected order.")
        dendrogram_details.update(correlation="pearson", linkage="complete", optimal_ordering=False)
        warnings.append(
            "Group order is descriptive and was derived from complete-linkage clustering of Pearson correlations "
            "between group means for only the selected marker panel."
        )
    elif expected_group_order is not None and list(expected_group_order) != resolved_order:
        raise RuntimeError("Observed Marker Expression group order differs from the expected order.")

    if len(resolved_order) == 1:
        warnings.append("Only one group is present; cross-group visual comparison is unavailable.")
    small_groups = {label: count for label, count in group_counts.items() if count < 10}
    if small_groups:
        warnings.append(f"Groups with fewer than 10 cells may yield unstable visual summaries: {small_groups!r}.")
    if mode.get("standard_scale", "none") != "none":
        warnings.append(
            f"standard_scale={mode['standard_scale']!r} rescales descriptive means and removes absolute magnitude "
            "comparability along the scaled dimension."
        )

    plotted_statistics = []
    empty_expressed = []
    for label in resolved_order:
        group_panel = panel[groups == label]
        for column, gene in enumerate(gene_names):
            gene_values = group_panel[:, column]
            row = {"group": label, "gene": gene, "n_cells": int(gene_values.size)}
            if mode["plot"] == "dotplot":
                expressed = gene_values > mode["expression_cutoff"]
                row["fraction_above_cutoff"] = float(expressed.mean())
                row["expression_cutoff"] = mode["expression_cutoff"]
                if mode["mean_only_expressed"]:
                    row["mean_expression"] = float(gene_values[expressed].mean()) if bool(expressed.any()) else None
                    row["mean_basis"] = "cells_above_cutoff"
                    if not bool(expressed.any()):
                        empty_expressed.append({"group": label, "gene": gene})
                else:
                    row["mean_expression"] = float(gene_values.mean())
                    row["mean_basis"] = "all_cells"
            elif mode["plot"] == "matrixplot":
                row["mean_expression"] = float(gene_values.mean())
            elif mode["plot"] == "tracksplot":
                row.update(
                    mean_expression=float(gene_values.mean()),
                    minimum=float(gene_values.min()),
                    maximum=float(gene_values.max()),
                )
            else:
                quantiles = np.quantile(gene_values, [0.0, 0.25, 0.5, 0.75, 1.0])
                row.update(
                    minimum=float(quantiles[0]),
                    q1=float(quantiles[1]),
                    median=float(quantiles[2]),
                    mean=float(gene_values.mean()),
                    q3=float(quantiles[3]),
                    maximum=float(quantiles[4]),
                )
            plotted_statistics.append(row)
    if empty_expressed:
        warnings.append(
            "Some group-gene combinations had no cells above expression_cutoff; their expressed-cell mean is "
            "reported as null instead of silently replacing Scanpy's undefined value."
        )

    work = ad.AnnData(
        X=panel.copy(),
        obs=pd.DataFrame(
            {groupby: pd.Categorical(groups.copy(), categories=resolved_order, ordered=True)},
            index=adata.obs_names.copy(),
        ),
        var=pd.DataFrame(index=pd.Index(gene_names, dtype=object)),
    )
    existing_figures = set(plt.get_fignums())
    figure = None
    png = None
    common = {
        "use_raw": False,
        "layer": None,
        "log": False,
        "show": False,
        "save": False,
    }
    rng_isolation = contextlib.ExitStack()
    try:
        if _rng_lock is not None:
            rng_isolation.enter_context(_rng_lock)
        rng_isolation.enter_context(isolated_numpy_random(random_seed))
        with scanpy_v1_context():
            if mode["plot"] == "dotplot":
                rendered = sc.pl.dotplot(
                    work,
                    var_names=list(gene_names),
                    groupby=groupby,
                    standard_scale=None if mode["standard_scale"] == "none" else mode["standard_scale"],
                    expression_cutoff=mode["expression_cutoff"],
                    mean_only_expressed=mode["mean_only_expressed"],
                    categories_order=list(resolved_order),
                    dendrogram=False,
                    return_fig=True,
                    **common,
                )
                make_figure = getattr(rendered, "make_figure", None)
                if not callable(make_figure):
                    raise RuntimeError("Scanpy dotplot backend did not return the legacy Matplotlib plot contract.")
                make_figure()
                figure = figure_from_result(rendered)
            elif mode["plot"] == "matrixplot":
                rendered = sc.pl.matrixplot(
                    work,
                    var_names=list(gene_names),
                    groupby=groupby,
                    standard_scale=None if mode["standard_scale"] == "none" else mode["standard_scale"],
                    categories_order=list(resolved_order),
                    dendrogram=False,
                    return_fig=True,
                    **common,
                )
                make_figure = getattr(rendered, "make_figure", None)
                if not callable(make_figure):
                    raise RuntimeError("Scanpy matrixplot backend did not return the legacy Matplotlib plot contract.")
                make_figure()
                figure = figure_from_result(rendered)
            elif mode["plot"] == "tracksplot":
                rendered = sc.pl.tracksplot(
                    work,
                    var_names=list(gene_names),
                    groupby=groupby,
                    dendrogram=False,
                    **common,
                )
                figure = figure_from_result(rendered)
            else:
                rendered = sc.pl.violin(
                    work,
                    keys=list(gene_names),
                    groupby=groupby,
                    density_norm=mode["density_norm"],
                    stripplot=mode["show_cells"],
                    jitter=0.4 if mode["show_cells"] else False,
                    size=1,
                    order=list(resolved_order),
                    multi_panel=len(gene_names) > 1,
                    use_raw=False,
                    layer=None,
                    log=mode["y_scale"] == "log",
                    show=False,
                    save=False,
                )
                figure = figure_from_result(rendered)
        figure_number = getattr(figure, "number", None)
        if figure_number in existing_figures:
            raise RuntimeError("Scanpy marker-expression backend reused a pre-existing caller figure.")
        figure.set_size_inches(width, height, forward=True)
        original_canvas = getattr(figure, "canvas", None)
        FigureCanvasAgg(figure)
        buffer = io.BytesIO()
        try:
            figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
            png = buffer.getvalue()
        finally:
            if original_canvas is not None:
                figure.set_canvas(original_canvas)
    finally:
        rng_isolation.close()
        for number in set(plt.get_fignums()) - existing_figures:
            plt.close(number)
    if not isinstance(png, bytes) or not png.startswith(png_signature):
        raise RuntimeError("Marker Expression renderer did not produce a valid PNG payload.")

    details = {
        "plot": mode,
        "genes": gene_names,
        "groupby": groupby,
        "group_order_policy": group_order,
        "observed_group_order": observed_order,
        "resolved_group_order": resolved_order,
        "group_cell_counts": group_counts,
        "missing_group_cells": 0,
        "source_kind": source_kind,
        "layer_name": source_layer_name,
        "source_gene_count": int(len(var_names)),
        "selected_panel_shape": [n_obs, len(gene_names)],
        "selected_panel_range": [float(panel.min()), float(panel.max())],
        "resource_estimate": resource_estimate,
        "dendrogram": dendrogram_details,
        "plotted_statistics": plotted_statistics,
        "rendering": {
            "figure_size_inches": [width, height],
            "dpi": 120,
            "bbox_inches": "tight",
            "scanpy_preset": "ScanpyV1",
            "autoshow": False,
            "scanpy_save": False,
            "random_seed": random_seed,
            "numpy_global_rng_state_restored": True,
        },
        "warnings": warnings,
    }
    return png, details


def marker_expression_plot_code(
    *,
    genes: str,
    groupby: str,
    plot: Mapping[str, Any],
    source_kind: str,
    layer_name: str | None,
    group_order: str,
    random_seed: int,
    resolved_group_order: list[str],
) -> str:
    implementation = inspect.getsource(_marker_expression_plot_impl)
    return (
        "import threading\n\n"
        "_marker_plot_rng_lock = threading.RLock()\n\n"
        f"{implementation}\n\n"
        "def plot_marker_expression(adata):\n"
        "    png, _details = _marker_expression_plot_impl(\n"
        "        adata,\n"
        f"        genes={genes!r},\n"
        f"        groupby={groupby!r},\n"
        f"        plot={dict(plot)!r},\n"
        f"        source_kind={source_kind!r},\n"
        f"        layer_name={layer_name!r},\n"
        f"        group_order={group_order!r},\n"
        f"        random_seed={random_seed!r},\n"
        f"        expected_group_order={list(resolved_group_order)!r},\n"
        "        _rng_lock=_marker_plot_rng_lock,\n"
        "    )\n"
        "    return png\n"
    )


__all__ = [
    "MAX_MARKER_PLOT_BYTES",
    "MARKER_PLOT_RNG_LOCK",
    "MAX_MARKER_PLOT_GENES",
    "MAX_MARKER_PLOT_STAT_ROWS",
    "PNG_SIGNATURE",
    "_marker_expression_plot_impl",
    "_plot_embedding_impl",
    "marker_expression_plot_code",
    "embedding_plot_code",
    "embedding_provenance",
]
