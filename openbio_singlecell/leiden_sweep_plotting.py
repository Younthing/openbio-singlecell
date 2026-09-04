from __future__ import annotations

import inspect
import textwrap
from typing import Any


def leiden_resolution_metrics_fingerprint(table: Any) -> str:
    import hashlib
    import json
    import math

    import pandas as pd

    if not isinstance(table, pd.DataFrame):
        raise TypeError("Leiden resolution metrics fingerprint requires a pandas DataFrame.")

    def plain(value):
        if value is None or value is pd.NA:
            return ["none", None]
        if hasattr(value, "item"):
            value = value.item()
        if isinstance(value, bool):
            return ["bool", value]
        if isinstance(value, str):
            return ["str", value]
        if isinstance(value, int):
            return ["int", value]
        if isinstance(value, float):
            if math.isnan(value):
                return ["none", None]
            if not math.isfinite(value):
                raise ValueError("Leiden resolution metrics fingerprint rejects infinity.")
            return ["float", value.hex()]
        raise TypeError(
            f"Leiden resolution metrics fingerprint does not support {type(value).__name__} values."
        )

    payload = {
        "schema": "openbio-singlecell/leiden-resolution-metrics/v1",
        "columns": table.columns.tolist(),
        "rows": [[plain(value) for value in row] for row in table.itertuples(index=False, name=None)],
    }
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _standalone_leiden_resolution_sweep_plot(
    table,
    *,
    producer,
    view=None,
    _return_details=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Leiden Resolution Sweep Plot"
    if view is None:
        mode = {"view": "quality_curves"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    if set(mode) != {"view"} or mode.get("view") not in {"quality_curves", "cluster_sizes"}:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if not isinstance(producer, Mapping) or set(producer) != {
        "operation",
        "parameters",
        "input_cells",
        "input_genes",
    }:
        raise ValueError(f"{operation} requires exact Leiden Resolution Sweep producer provenance.")
    if producer["operation"] != "leiden_resolution_sweep":
        raise ValueError(f"{operation} requires the Leiden Resolution Sweep producer.")
    parameters = producer["parameters"]
    required_parameters = {
        "resolutions",
        "keys",
        "stability_repeats",
        "total_runs",
        "fixed_policy",
        "table_content_fingerprint_sha256",
    }
    if not isinstance(parameters, Mapping) or not required_parameters.issubset(parameters):
        raise ValueError(f"{operation} producer parameters are incomplete.")
    input_cells, input_genes = producer["input_cells"], producer["input_genes"]
    for name, value, minimum in (("input_cells", input_cells, 1), ("input_genes", input_genes, 0)):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or int(value) < minimum
        ):
            raise ValueError(f"{operation} {name} provenance is invalid.")
    fixed_policy = parameters["fixed_policy"]
    if not isinstance(fixed_policy, Mapping) or dict(fixed_policy) != {
        "flavor": "igraph",
        "directed": False,
        "use_weights": True,
        "objective_function": "modularity",
    }:
        raise ValueError(f"{operation} requires the canonical Leiden sweep policy provenance.")

    columns = [
        "resolution",
        "key",
        "n_clusters",
        "min_cluster_size",
        "median_cluster_size",
        "max_cluster_size",
        "smallest_cluster_fraction",
        "largest_cluster_fraction",
        "singleton_clusters",
        "modularity",
        "stability_repeats",
        "stability_mean_ari",
        "stability_min_ari",
        "stability_max_ari",
        "adjacent_previous_resolution",
        "adjacent_resolution_ari",
    ]
    if not isinstance(table, pd.DataFrame):
        raise TypeError(f"{operation} requires a pandas DataFrame.")
    if table.columns.tolist() != columns:
        raise ValueError(f"{operation} requires canonical columns {columns!r} in that order.")
    if table.empty or len(table) > 100:
        raise ValueError(f"{operation} requires between 1 and 100 stored resolution rows.")

    def numeric(column):
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} must be real numeric evidence.")
        return series.to_numpy(dtype=float, na_value=np.nan)

    resolutions = numeric("resolution")
    if not bool(np.isfinite(resolutions).all()) or bool((resolutions < 0).any()) or not bool(
        np.all(np.diff(resolutions) > 0)
    ):
        raise ValueError(f"{operation} resolution axis must be finite, nonnegative, unique, and increasing.")
    declared_resolutions = parameters["resolutions"]
    try:
        declared_resolutions = np.asarray(declared_resolutions, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{operation} producer resolution axis is invalid.") from error
    if not np.array_equal(resolutions, declared_resolutions):
        raise ValueError(f"{operation} resolution axis does not match producer provenance.")
    keys = table["key"].tolist()
    if not all(isinstance(key, str) and key for key in keys) or len(set(keys)) != len(keys):
        raise ValueError(f"{operation} stored partition keys must be unique nonblank strings.")
    if keys != list(parameters["keys"]):
        raise ValueError(f"{operation} partition-key axis does not match producer provenance.")

    integer_columns = (
        "n_clusters",
        "min_cluster_size",
        "max_cluster_size",
        "singleton_clusters",
        "stability_repeats",
    )
    integers = {}
    for column in integer_columns:
        values = numeric(column)
        if not bool(np.isfinite(values).all()) or not bool(np.equal(values, np.rint(values)).all()):
            raise ValueError(f"{operation} {column!r} must contain finite integers.")
        integers[column] = values.astype(np.int64)
    n_clusters = integers["n_clusters"]
    minimum_size = integers["min_cluster_size"]
    maximum_size = integers["max_cluster_size"]
    singleton = integers["singleton_clusters"]
    repeats = integers["stability_repeats"]
    median_size = numeric("median_cluster_size")
    if (
        bool((n_clusters < 1).any())
        or bool((minimum_size < 1).any())
        or bool((maximum_size < minimum_size).any())
        or bool((maximum_size > int(input_cells)).any())
        or not bool(np.isfinite(median_size).all())
        or bool(((median_size < minimum_size) | (median_size > maximum_size)).any())
        or bool(((singleton < 0) | (singleton > n_clusters)).any())
        or bool(((singleton > 0) & (minimum_size != 1)).any())
    ):
        raise ValueError(f"{operation} stored cluster-size diagnostics are inconsistent.")
    expected_repeats = parameters["stability_repeats"]
    expected_total_runs = parameters["total_runs"]
    if (
        isinstance(expected_repeats, (bool, np.bool_))
        or not isinstance(expected_repeats, (int, np.integer))
        or int(expected_repeats) < 1
        or not bool(np.equal(repeats, int(expected_repeats)).all())
        or isinstance(expected_total_runs, (bool, np.bool_))
        or not isinstance(expected_total_runs, (int, np.integer))
        or int(expected_total_runs) != len(table) * int(expected_repeats)
    ):
        raise ValueError(f"{operation} repeat-start provenance is inconsistent.")

    smallest_fraction = numeric("smallest_cluster_fraction")
    largest_fraction = numeric("largest_cluster_fraction")
    if (
        not bool(np.isfinite(smallest_fraction).all())
        or not bool(np.isfinite(largest_fraction).all())
        or bool(((smallest_fraction <= 0) | (largest_fraction > 1)).any())
        or bool((smallest_fraction > largest_fraction).any())
        or not bool(np.allclose(smallest_fraction, minimum_size / int(input_cells), rtol=0, atol=1e-12))
        or not bool(np.allclose(largest_fraction, maximum_size / int(input_cells), rtol=0, atol=1e-12))
    ):
        raise ValueError(f"{operation} cluster fractions do not match stored sizes and cell count.")
    modularity = numeric("modularity")
    if not bool(np.isfinite(modularity).all()) or bool(((modularity < -1) | (modularity > 1)).any()):
        raise ValueError(f"{operation} modularity must contain finite values in [-1, 1].")
    stability_mean = numeric("stability_mean_ari")
    stability_min = numeric("stability_min_ari")
    stability_max = numeric("stability_max_ari")
    stability_stack = np.vstack([stability_min, stability_mean, stability_max])
    if int(expected_repeats) == 1:
        if not bool(np.isnan(stability_stack).all()):
            raise ValueError(f"{operation} one-start rows must not claim repeat-start ARI evidence.")
    elif (
        not bool(np.isfinite(stability_stack).all())
        or bool(((stability_stack < -1) | (stability_stack > 1)).any())
        or bool((stability_min > stability_mean).any())
        or bool((stability_mean > stability_max).any())
    ):
        raise ValueError(f"{operation} repeat-start ARI diagnostics are inconsistent.")
    adjacent_previous = numeric("adjacent_previous_resolution")
    adjacent_ari = numeric("adjacent_resolution_ari")
    if not (np.isnan(adjacent_previous[0]) and np.isnan(adjacent_ari[0])):
        raise ValueError(f"{operation} first row cannot have adjacent-resolution evidence.")
    if len(table) > 1 and (
        not bool(np.isfinite(adjacent_previous[1:]).all())
        or not bool(np.array_equal(adjacent_previous[1:], resolutions[:-1]))
        or not bool(np.isfinite(adjacent_ari[1:]).all())
        or bool(((adjacent_ari[1:] < -1) | (adjacent_ari[1:] > 1)).any())
    ):
        raise ValueError(f"{operation} adjacent-resolution ARI diagnostics are inconsistent.")
    content_fingerprint = leiden_resolution_metrics_fingerprint(table)
    if parameters["table_content_fingerprint_sha256"] != content_fingerprint:
        raise ValueError(f"{operation} table current-content fingerprint does not match producer provenance.")

    if mode["view"] == "quality_curves":
        figure = Figure(figsize=(11.0, 8.0), constrained_layout=True)
        FigureCanvasAgg(figure)
        axes = figure.subplots(2, 2, sharex=True)
        axes[0, 0].plot(resolutions, n_clusters, marker="o", color="#4C78A8")
        axes[0, 0].set_title("Communities")
        axes[0, 0].set_ylabel("Clusters")
        axes[0, 1].plot(resolutions, modularity, marker="o", color="#F58518")
        axes[0, 1].set_title("Modularity")
        axes[0, 1].set_ylabel("Modularity")
        if int(expected_repeats) > 1:
            axes[1, 0].fill_between(resolutions, stability_min, stability_max, color="#54A24B", alpha=0.2)
            axes[1, 0].plot(resolutions, stability_mean, marker="o", color="#54A24B")
        else:
            axes[1, 0].text(0.5, 0.5, "Not assessed (one start)", ha="center", va="center", transform=axes[1, 0].transAxes)
        axes[1, 0].set_title("Repeat-start ARI")
        axes[1, 0].set_ylabel("Adjusted Rand index")
        axes[1, 1].plot(resolutions[1:], adjacent_ari[1:], marker="o", color="#E45756")
        axes[1, 1].set_title("Adjacent-resolution ARI")
        axes[1, 1].set_ylabel("Adjusted Rand index")
        for axis in axes.ravel():
            axis.set_xlabel("Resolution")
            axis.set_xticks(resolutions)
            axis.grid(alpha=0.2)
        title = "Leiden resolution quality diagnostics"
    else:
        figure = Figure(figsize=(13.0, 4.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axes = figure.subplots(1, 3)
        axes[0].fill_between(resolutions, minimum_size, maximum_size, color="#4C78A8", alpha=0.2)
        axes[0].plot(resolutions, median_size, marker="o", color="#4C78A8", label="Median")
        axes[0].set_title("Cluster size range")
        axes[0].set_ylabel("Cells")
        axes[0].legend(frameon=False)
        axes[1].fill_between(resolutions, smallest_fraction, largest_fraction, color="#F58518", alpha=0.2)
        axes[1].set_title("Cluster fraction range")
        axes[1].set_ylabel("Fraction of cells")
        axes[2].bar(resolutions.astype(str), singleton, color="#E45756")
        axes[2].set_title("Singleton communities")
        axes[2].set_ylabel("Singleton clusters")
        for axis in axes:
            axis.set_xlabel("Resolution")
            axis.grid(axis="y", alpha=0.2)
        title = "Leiden resolution cluster-size diagnostics"
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")
    details = {
        "view": mode["view"],
        "view_parameters": mode,
        "title": title,
        "resolutions": resolutions.tolist(),
        "partition_keys": keys,
        "cluster_counts": n_clusters.astype(int).tolist(),
        "cluster_size_ranges": [
            [int(low), float(middle), int(high)]
            for low, middle, high in zip(minimum_size, median_size, maximum_size, strict=True)
        ],
        "cluster_fraction_ranges": [
            [float(low), float(high)]
            for low, high in zip(smallest_fraction, largest_fraction, strict=True)
        ],
        "singleton_clusters": singleton.astype(int).tolist(),
        "modularity": modularity.tolist(),
        "stability_repeats": int(expected_repeats),
        "stability_mean_ari": None if int(expected_repeats) == 1 else stability_mean.tolist(),
        "stability_min_ari": None if int(expected_repeats) == 1 else stability_min.tolist(),
        "stability_max_ari": None if int(expected_repeats) == 1 else stability_max.tolist(),
        "adjacent_resolution_ari": [None, *adjacent_ari[1:].tolist()],
        "table_content_fingerprint_sha256": content_fingerprint,
        "automatic_resolution_selection": False,
    }
    return (png, details) if _return_details else png


def run_leiden_resolution_sweep_plot(
    table: Any,
    *,
    producer: dict[str, Any],
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_leiden_resolution_sweep_plot(
        table,
        producer=producer,
        view=view,
        _return_details=True,
    )


def leiden_resolution_sweep_plot_code(*, producer: dict[str, Any], view: dict[str, Any]) -> str:
    fingerprint = textwrap.dedent(inspect.getsource(leiden_resolution_metrics_fingerprint)).strip()
    implementation = textwrap.dedent(inspect.getsource(_standalone_leiden_resolution_sweep_plot)).strip()
    return f"""from __future__ import annotations

{fingerprint}


{implementation}


def plot_leiden_resolution_sweep(resolution_metrics):
    return _standalone_leiden_resolution_sweep_plot(
        resolution_metrics,
        producer={producer!r},
        view={view!r},
    )
"""


__all__ = [
    "leiden_resolution_sweep_plot_code",
    "leiden_resolution_metrics_fingerprint",
    "run_leiden_resolution_sweep_plot",
]
