from __future__ import annotations

import inspect
import textwrap


def _standalone_pca_loadings_plot(adata, *, component=1, n_genes=20, _return_details=False):
    import io

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "PCA Loadings Plot"
    if not hasattr(adata, "varm") or not hasattr(adata, "var_names"):
        raise TypeError(f"{operation} requires an AnnData-like input.")
    if "PCs" not in adata.varm:
        raise ValueError(f"{operation} requires stored adata.varm['PCs'] loadings.")
    if not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique feature identifiers.")
    loadings = np.asarray(adata.varm["PCs"])
    if loadings.ndim != 2 or loadings.shape[0] != int(adata.n_vars) or loadings.shape[1] < 1:
        raise ValueError(f"{operation} loadings must have shape (n_vars, n_components).")
    if (
        not np.issubdtype(loadings.dtype, np.number)
        or np.issubdtype(loadings.dtype, np.bool_)
        or np.iscomplexobj(loadings)
    ):
        raise TypeError(f"{operation} loadings must contain real numeric values.")
    if not bool(np.isfinite(loadings).all()):
        raise ValueError(f"{operation} loadings contain non-finite values.")
    if isinstance(component, (bool, np.bool_)) or not isinstance(component, (int, np.integer)):
        raise TypeError(f"{operation} component must be an integer.")
    component = int(component)
    available_components = int(loadings.shape[1])
    if component < 1 or component > available_components:
        raise ValueError(
            f"{operation} requested component {component}, but only {available_components} components are stored."
        )
    if isinstance(n_genes, (bool, np.bool_)) or not isinstance(n_genes, (int, np.integer)):
        raise TypeError(f"{operation} n_genes must be an integer.")
    n_genes = int(n_genes)
    if n_genes < 1:
        raise ValueError(f"{operation} n_genes must be positive.")
    if n_genes > 50:
        raise ValueError(f"{operation} supports at most 50 genes in one readable static plot.")

    values = np.asarray(loadings[:, component - 1], dtype=float)
    plotted_count = min(n_genes, int(adata.n_vars))
    order = np.argsort(-np.abs(values), kind="stable")[:plotted_count]
    genes = [str(adata.var_names[index]) for index in order]
    if len(set(genes)) != len(genes):
        raise ValueError(f"{operation} feature identifiers collide after display conversion.")
    selected = values[order]

    height = max(4.0, min(14.0, 1.8 + 0.32 * plotted_count))
    figure = Figure(figsize=(8.0, height), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    positions = np.arange(plotted_count)
    colors = np.where(selected >= 0, "#3274a1", "#d65f5f")
    axis.barh(positions, selected, color=colors)
    axis.axvline(0.0, color="#444444", linewidth=0.8)
    axis.set_yticks(positions, genes)
    axis.invert_yaxis()
    axis.set_xlabel("PCA loading")
    axis.set_title(f"PC{component} top absolute loadings")
    axis.grid(axis="x", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")
    details = {
        "component": component,
        "available_components": available_components,
        "requested_genes": n_genes,
        "plotted_genes": plotted_count,
        "genes": genes,
        "loadings": selected.tolist(),
        "ordering": "descending absolute loading; stable feature order for ties",
        "positive_loading_count": int((selected > 0).sum()),
        "negative_loading_count": int((selected < 0).sum()),
    }
    return (png, details) if _return_details else png


def _standalone_neighbor_graph_diagnostics_plot(
    adata,
    *,
    neighbors_key="neighbors",
    max_working_memory_gib=4.0,
    _return_details=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components

    operation = "Neighbor Graph Diagnostics Plot"
    if not hasattr(adata, "uns") or not hasattr(adata, "obsp") or not hasattr(adata, "obs_names"):
        raise TypeError(f"{operation} requires an AnnData-like input.")
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires in-memory AnnData.")
    if int(adata.n_obs) < 2 or not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires at least two uniquely identified observations.")
    if not isinstance(neighbors_key, str) or not neighbors_key.strip():
        raise ValueError(f"{operation} neighbors_key must be a nonblank string.")
    neighbors_key = neighbors_key.strip()
    if isinstance(max_working_memory_gib, (bool, np.bool_)) or not isinstance(
        max_working_memory_gib, (int, float, np.integer, np.floating)
    ):
        raise TypeError(f"{operation} max_working_memory_gib must be numeric.")
    max_working_memory_gib = float(max_working_memory_gib)
    if not np.isfinite(max_working_memory_gib) or max_working_memory_gib <= 0:
        raise ValueError(f"{operation} max_working_memory_gib must be finite and positive.")
    metadata = adata.uns.get(neighbors_key)
    if not isinstance(metadata, Mapping):
        raise ValueError(f"{operation} neighbor metadata not found in uns: {neighbors_key!r}.")

    def matrix_from_pointer(pointer_name, label):
        pointer = metadata.get(pointer_name)
        if not isinstance(pointer, str) or not pointer.strip():
            raise ValueError(f"{operation} uns[{neighbors_key!r}] has no valid {pointer_name}.")
        pointer = pointer.strip()
        if pointer not in adata.obsp:
            raise ValueError(f"{operation} {label} matrix not found in obsp: {pointer!r}.")
        matrix = adata.obsp[pointer]
        if not sparse.issparse(matrix) or tuple(matrix.shape) != (int(adata.n_obs), int(adata.n_obs)):
            raise ValueError(f"{operation} {label} matrix must be sparse with shape (n_obs, n_obs).")
        return pointer, matrix

    connectivity_key, raw_connectivities = matrix_from_pointer("connectivities_key", "connectivity")
    distance_key, raw_distances = matrix_from_pointer("distances_key", "distance")

    def sparse_storage_bytes(matrix):
        total = int(np.asarray(matrix.data).nbytes)
        for name in ("indices", "indptr", "row", "col"):
            values = getattr(matrix, name, None)
            if values is not None:
                total += int(np.asarray(values).nbytes)
        return total

    input_sparse_bytes = sparse_storage_bytes(raw_connectivities) + sparse_storage_bytes(raw_distances)
    estimated_working_bytes = input_sparse_bytes * 6 + int(adata.n_obs) * 64
    working_limit_bytes = int(max_working_memory_gib * 1024**3)
    if estimated_working_bytes > working_limit_bytes:
        raise MemoryError(
            f"{operation} estimated sparse working memory {estimated_working_bytes / 1024**3:.3f} GiB "
            f"exceeds max_working_memory_gib={max_working_memory_gib:.3f}."
        )

    def canonical_matrix(matrix, label):
        matrix = matrix.tocsr(copy=True)
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        values = np.asarray(matrix.data)
        if (
            not np.issubdtype(values.dtype, np.number)
            or np.issubdtype(values.dtype, np.bool_)
            or np.iscomplexobj(values)
        ):
            raise TypeError(f"{operation} {label} matrix must contain real numeric values.")
        if values.size and (not bool(np.isfinite(values).all()) or bool((values < 0).any())):
            raise ValueError(f"{operation} {label} matrix must contain finite nonnegative values.")
        if bool(np.any(matrix.diagonal() != 0)):
            raise ValueError(f"{operation} {label} matrix must have a zero diagonal.")
        return matrix

    connectivities = canonical_matrix(raw_connectivities, "connectivity")
    distances = canonical_matrix(raw_distances, "distance")
    delta = connectivities - connectivities.T
    delta.eliminate_zeros()
    if delta.nnz and float(np.max(np.abs(delta.data))) > 1e-8:
        raise ValueError(f"{operation} connectivity matrix must be symmetric.")
    connectivities = ((connectivities + connectivities.T) * 0.5).tocsr()
    connectivities.eliminate_zeros()
    if connectivities.nnz == 0:
        raise ValueError(f"{operation} connectivity graph has no positive edges.")

    degree = np.asarray((connectivities > 0).sum(axis=1), dtype=int).ravel()
    weighted_degree = np.asarray(connectivities.sum(axis=1), dtype=float).ravel()
    component_count, labels = connected_components(connectivities, directed=False)
    component_sizes = np.sort(np.bincount(labels, minlength=component_count))[::-1]
    distance_values = np.asarray(distances.data, dtype=float)
    if distance_values.size == 0:
        raise ValueError(f"{operation} distance matrix has no positive stored distances.")

    figure = Figure(figsize=(12.0, 4.0), constrained_layout=True)
    FigureCanvasAgg(figure)
    axes = figure.subplots(1, 3)
    axes[0].hist(degree, bins=min(40, max(1, int(degree.max()) + 1)), color="#3274a1")
    axes[0].set_title("Graph degree")
    axes[0].set_xlabel("Positive neighbors")
    axes[0].set_ylabel("Cells")
    axes[1].hist(distance_values, bins=40, color="#55a868")
    axes[1].set_title("Stored neighbor distances")
    axes[1].set_xlabel("Distance")
    axes[1].set_ylabel("Directed edges")
    shown_components = component_sizes[:100]
    axes[2].bar(np.arange(1, len(shown_components) + 1), shown_components, color="#c44e52")
    axes[2].set_title("Connected component sizes")
    axes[2].set_xlabel("Component rank")
    axes[2].set_ylabel("Cells")
    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")

    preview_limit = 100
    details = {
        "neighbors_key": neighbors_key,
        "connectivities_key": connectivity_key,
        "distances_key": distance_key,
        "observations": int(adata.n_obs),
        "positive_directed_connectivity_entries": int(connectivities.nnz),
        "positive_directed_distance_entries": int(distance_values.size),
        "connected_components": int(component_count),
        "component_sizes": component_sizes[:preview_limit].astype(int).tolist(),
        "component_sizes_truncated": len(component_sizes) > preview_limit,
        "degree_values": degree[:preview_limit].astype(int).tolist(),
        "degree_values_count": int(degree.size),
        "degree_values_truncated": degree.size > preview_limit,
        "weighted_degree_range": [float(weighted_degree.min()), float(weighted_degree.max())],
        "positive_distance_values": distance_values[:preview_limit].tolist(),
        "positive_distance_values_count": int(distance_values.size),
        "positive_distance_values_truncated": distance_values.size > preview_limit,
        "neighbor_parameters": {
            str(key): value
            for key, value in (metadata.get("params") or {}).items()
            if value is None or isinstance(value, (str, bool, int, float))
        }
        if isinstance(metadata.get("params"), Mapping)
        else {},
        "memory_preflight": {
            "input_sparse_bytes": input_sparse_bytes,
            "estimated_working_bytes": estimated_working_bytes,
            "max_working_memory_gib": max_working_memory_gib,
        },
        "fixed_render_pixels": 12 * 120 * 4 * 120,
    }
    return (png, details) if _return_details else png


def pca_loadings_plot_code(*, component: int, n_genes: int) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_pca_loadings_plot)).strip()
    return f"""{implementation}


def plot_pca_loadings(adata):
    return _standalone_pca_loadings_plot(adata, component={component!r}, n_genes={n_genes!r})
"""


def neighbor_graph_diagnostics_plot_code(*, neighbors_key: str, max_working_memory_gib: float) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_neighbor_graph_diagnostics_plot)).strip()
    return f"""{implementation}


def plot_neighbor_graph_diagnostics(adata):
    return _standalone_neighbor_graph_diagnostics_plot(
        adata,
        neighbors_key={neighbors_key!r},
        max_working_memory_gib={max_working_memory_gib!r},
    )
"""


__all__ = [
    "_standalone_neighbor_graph_diagnostics_plot",
    "_standalone_pca_loadings_plot",
    "neighbor_graph_diagnostics_plot_code",
    "pca_loadings_plot_code",
]
