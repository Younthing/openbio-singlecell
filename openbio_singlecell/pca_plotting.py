from __future__ import annotations

import inspect
from collections.abc import Mapping
from io import BytesIO
from typing import Any

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.ticker import PercentFormatter


def _pca_variance_plot_impl(adata: Any, *, n_pcs: int = 0) -> tuple[bytes, dict[str, Any]]:
    if isinstance(n_pcs, (bool, np.bool_)) or not isinstance(n_pcs, (int, np.integer)):
        raise TypeError("PCA Variance Plot n_pcs must be an integer.")
    n_pcs = int(n_pcs)
    if n_pcs < 0:
        raise ValueError("PCA Variance Plot n_pcs must be non-negative.")

    metadata = adata.uns.get("pca")
    if not isinstance(metadata, Mapping) or "variance_ratio" not in metadata:
        raise ValueError("PCA Variance Plot requires stored adata.uns['pca']['variance_ratio'] values.")
    stored = np.asarray(metadata["variance_ratio"])
    if stored.ndim != 1 or stored.size == 0:
        raise ValueError("PCA Variance Plot variance_ratio must be a non-empty one-dimensional array.")
    if not np.issubdtype(stored.dtype, np.number) or np.issubdtype(stored.dtype, np.bool_) or np.iscomplexobj(stored):
        raise TypeError("PCA Variance Plot variance_ratio must contain real numeric values.")
    ratios = np.asarray(stored, dtype=float)
    if not bool(np.isfinite(ratios).all()):
        raise ValueError("PCA Variance Plot variance_ratio contains non-finite values.")
    if bool((ratios < 0).any()):
        raise ValueError("PCA Variance Plot variance_ratio contains negative values.")
    total_ratio = float(ratios.sum(dtype=float))
    if total_ratio > 1.0 and not bool(np.isclose(total_ratio, 1.0, rtol=1e-6, atol=1e-8)):
        raise ValueError("PCA Variance Plot cumulative variance_ratio cannot exceed 1.")

    available = int(ratios.size)
    if n_pcs > available:
        raise ValueError(f"PCA Variance Plot requested {n_pcs} PCs, but only {available} are stored.")
    plotted = available if n_pcs == 0 else n_pcs
    selected = ratios[:plotted]
    cumulative = np.cumsum(selected, dtype=float)
    components = np.arange(1, plotted + 1)

    width = float(min(14.0, max(7.0, 6.0 + 0.12 * plotted)))
    figure = Figure(figsize=(width, 5.0))
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    axis.bar(components, selected, width=0.8, color="#4c78a8", label="Explained variance")
    axis.set_xlabel("Principal component")
    axis.set_ylabel("Explained variance ratio")
    axis.set_xlim(0.25, plotted + 0.75)
    axis.set_ylim(0.0, min(1.0, max(0.01, float(selected.max()) * 1.12)))
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    tick_step = max(1, int(np.ceil(plotted / 12)))
    ticks = components[::tick_step]
    if int(ticks[-1]) != plotted:
        ticks = np.append(ticks, plotted)
    axis.set_xticks(ticks)
    axis.grid(axis="y", alpha=0.25, linewidth=0.8)

    cumulative_axis = axis.twinx()
    cumulative_axis.plot(
        components,
        cumulative,
        color="#f58518",
        marker="o",
        markersize=3.5,
        linewidth=2.0,
        label="Cumulative variance",
    )
    cumulative_axis.set_ylabel("Cumulative explained variance ratio")
    cumulative_axis.set_ylim(0.0, 1.0)
    cumulative_axis.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    handles, labels = axis.get_legend_handles_labels()
    cumulative_handles, cumulative_labels = cumulative_axis.get_legend_handles_labels()
    axis.legend(handles + cumulative_handles, labels + cumulative_labels, loc="best", frameon=False)
    figure.tight_layout()

    buffer = BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError("PCA Variance Plot renderer did not produce a valid PNG payload.")

    warnings = []
    if plotted < available:
        warnings.append(f"Displayed the first {plotted} of {available} stored principal components.")
    return png, {
        "variance_ratio_source": "adata.uns['pca']['variance_ratio']",
        "requested_n_pcs": n_pcs,
        "available_components": available,
        "plotted_components": plotted,
        "component_numbers": components.tolist(),
        "variance_ratio": selected.tolist(),
        "cumulative_variance_ratio": cumulative.tolist(),
        "available_cumulative_variance_ratio": total_ratio,
        "rendering": {
            "figure_size_inches": [width, 5.0],
            "dpi": 120,
            "bar_series": "explained_variance_ratio",
            "line_series": "cumulative_explained_variance_ratio",
        },
        "warnings": warnings,
    }


def pca_variance_plot_code(*, n_pcs: int) -> str:
    implementation = inspect.getsource(_pca_variance_plot_impl)
    return (
        "from collections.abc import Mapping\n"
        "from io import BytesIO\n"
        "from typing import Any\n\n"
        "import numpy as np\n"
        "from matplotlib.backends.backend_agg import FigureCanvasAgg\n"
        "from matplotlib.figure import Figure\n"
        "from matplotlib.ticker import PercentFormatter\n\n\n"
        f"{implementation}\n\n"
        "def plot_pca_variance(adata):\n"
        f"    png, _details = _pca_variance_plot_impl(adata, n_pcs={n_pcs!r})\n"
        "    return png\n"
    )


__all__ = ["_pca_variance_plot_impl", "pca_variance_plot_code"]
