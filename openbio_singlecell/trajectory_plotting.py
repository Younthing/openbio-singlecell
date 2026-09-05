from __future__ import annotations

import inspect
import textwrap
from collections.abc import Mapping
from typing import Any

from .trajectory_analysis import (
    _TA_DIFFMAP_PROVENANCE_KEY,
    _TA_DIFFMAP_SCHEMA,
    _TA_DPT_PROVENANCE_KEY,
    _TA_DPT_SCHEMA,
    _TA_PAGA_PROVENANCE_KEY,
    _TA_PAGA_SCHEMA,
    _ta_clean_text,
    _ta_diffmap_output_fingerprint,
    _ta_graph,
    _ta_hash_array,
    _ta_hash_sparse,
    _ta_hash_text_records,
    _ta_json_value,
    _ta_label_record,
    _ta_obs_names,
    _ta_paga_matrix,
    _ta_paga_output_fingerprint,
    _ta_partition,
    _ta_sparse_graph_matrix,
    _ta_transition_matrix,
    _ta_validate_dense_float,
    _ta_validate_diffmap_bundle,
    _ta_validate_diffmap_eigenpairs,
)


def _standalone_validate_paga_plot_input(adata):
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.sparse import csgraph

    operation = "PAGA Plot"
    provenance = adata.uns.get(_TA_PAGA_PROVENANCE_KEY)
    expected_fields = {
        "schema",
        "neighbors_key",
        "graph_fingerprint_sha256",
        "observation_fingerprint_sha256",
        "groupby",
        "categories",
        "group_counts",
        "partition_fingerprint_sha256",
        "model",
        "use_rna_velocity",
        "output_fingerprint_sha256",
        "scanpy_version",
    }
    if not isinstance(provenance, dict) or set(provenance) != expected_fields:
        raise ValueError(f"{operation} requires exact OpenBio PAGA provenance.")
    if (
        provenance["schema"] != _TA_PAGA_SCHEMA
        or provenance["model"] != "v1.2"
        or provenance["use_rna_velocity"] is not False
    ):
        raise ValueError(f"{operation} PAGA producer identity or fixed model policy is invalid.")
    graph = _ta_graph(
        adata,
        neighbors_key=provenance["neighbors_key"],
        operation=operation,
        numpy=np,
        scipy_sparse=sparse,
        scipy_csgraph=csgraph,
    )
    if (
        provenance["graph_fingerprint_sha256"] != graph["graph_fingerprint_sha256"]
        or provenance["observation_fingerprint_sha256"] != graph["observation_fingerprint_sha256"]
    ):
        raise ValueError(f"{operation} named graph or observation axis differs from the producer evidence.")
    partition = _ta_partition(
        adata,
        groupby=provenance["groupby"],
        operation=operation,
        pandas=pd,
        numpy=np,
    )
    stored_categories = provenance["categories"]
    if hasattr(stored_categories, "tolist"):
        stored_categories = stored_categories.tolist()
    if not isinstance(stored_categories, list) or not all(isinstance(record, Mapping) for record in stored_categories):
        raise ValueError(f"{operation} stored category axis is malformed.")
    stored_categories = [_ta_json_value(dict(record)) for record in stored_categories]
    stored_counts = provenance["group_counts"]
    if hasattr(stored_counts, "tolist"):
        stored_counts = stored_counts.tolist()
    if not isinstance(stored_counts, list):
        raise ValueError(f"{operation} stored group counts are malformed.")
    stored_counts = [int(value) for value in stored_counts]
    if (
        stored_categories != partition["categories"]
        or stored_counts != partition["counts"]
        or provenance["partition_fingerprint_sha256"] != partition["partition_fingerprint_sha256"]
    ):
        raise ValueError(f"{operation} categorical partition differs from the producer evidence.")
    bundle = adata.uns.get("paga")
    if not isinstance(bundle, dict) or bundle.get("groups") != partition["groupby"]:
        raise ValueError(f"{operation} canonical uns['paga'] bundle is missing or uses another group axis.")
    n_groups = len(partition["categories"])
    connectivities = _ta_paga_matrix(
        bundle.get("connectivities"),
        n_groups=n_groups,
        name="connectivities",
        symmetric=True,
        numpy=np,
        scipy_sparse=sparse,
    )
    tree = _ta_paga_matrix(
        bundle.get("connectivities_tree"),
        n_groups=n_groups,
        name="connectivities_tree",
        symmetric=False,
        numpy=np,
        scipy_sparse=sparse,
    )
    sizes = np.asarray(adata.uns.get(f"{partition['groupby']}_sizes"))
    if sizes.shape != (n_groups,) or not np.issubdtype(sizes.dtype, np.integer):
        raise ValueError(f"{operation} group-size sidecar is missing or invalid.")
    if [int(value) for value in sizes] != partition["counts"]:
        raise ValueError(f"{operation} group-size sidecar differs from the categorical partition.")
    reciprocal = tree.multiply(tree.T)
    reciprocal.eliminate_zeros()
    if reciprocal.nnz:
        raise ValueError(f"{operation} spanning forest stores reciprocal directed edges.")
    tree_coo = tree.tocoo()
    if tree_coo.nnz:
        parent_weights = np.asarray(connectivities[tree_coo.row, tree_coo.col]).ravel()
        if bool((parent_weights <= 0).any()) or not bool(np.allclose(parent_weights, tree_coo.data, rtol=1e-8, atol=1e-12)):
            raise ValueError(f"{operation} spanning forest is not an exact subgraph of connectivities.")
    graph_components, _ = csgraph.connected_components(connectivities, directed=False)
    forest_components, _ = csgraph.connected_components((tree + tree.T).tocsr(), directed=False)
    if forest_components != graph_components or tree.nnz != n_groups - graph_components:
        raise ValueError(f"{operation} stored tree is not a spanning forest of the abstraction.")
    fingerprint = _ta_paga_output_fingerprint(connectivities, tree, sizes, partition, numpy=np)
    if fingerprint != provenance["output_fingerprint_sha256"]:
        raise ValueError(f"{operation} output fingerprint differs from the stored PAGA bundle.")
    return graph, partition, connectivities, tree, sizes, dict(provenance)


def _standalone_paga_plot(adata, *, min_connectivity=0.0, _return_details=False):
    import io
    import math

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    operation = "PAGA Plot"
    if isinstance(min_connectivity, bool) or not isinstance(min_connectivity, (int, float)):
        raise TypeError(f"{operation} min_connectivity must be numeric.")
    min_connectivity = float(min_connectivity)
    if not math.isfinite(min_connectivity) or not 0.0 <= min_connectivity <= 1.0:
        raise ValueError(f"{operation} min_connectivity must be finite and between zero and one.")
    graph, partition, connectivities, tree, sizes, provenance = _standalone_validate_paga_plot_input(adata)
    groups = [record["display"] for record in partition["categories"]]
    n_groups = len(groups)
    angles = np.linspace(0.0, 2.0 * np.pi, n_groups, endpoint=False) + np.pi / 2.0
    coordinates = np.column_stack([np.cos(angles), np.sin(angles)])
    forest_pairs = {
        tuple(sorted((int(left), int(right))))
        for left, right in zip(tree.tocoo().row, tree.tocoo().col, strict=True)
    }
    available = []
    matrix = connectivities.tocoo()
    for left, right, weight in zip(matrix.row, matrix.col, matrix.data, strict=True):
        left, right = int(left), int(right)
        if left >= right:
            continue
        available.append(
            {
                "group_a": groups[left],
                "group_b": groups[right],
                "connectivity": float(weight),
                "in_spanning_forest": (left, right) in forest_pairs,
                "left_index": left,
                "right_index": right,
            }
        )
    available.sort(key=lambda record: (-record["connectivity"], record["group_a"], record["group_b"]))
    displayed = [record for record in available if record["connectivity"] >= min_connectivity]
    figure = Figure(figsize=(7.0, 7.0), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    for record in displayed:
        left, right = record["left_index"], record["right_index"]
        axis.plot(
            coordinates[[left, right], 0],
            coordinates[[left, right], 1],
            color="#4C78A8" if record["in_spanning_forest"] else "#9E9E9E",
            linewidth=0.7 + 4.0 * record["connectivity"],
            linestyle="-" if record["in_spanning_forest"] else "--",
            alpha=0.85,
            zorder=1,
        )
    marker_sizes = 300.0 + 900.0 * np.sqrt(sizes / sizes.max())
    axis.scatter(coordinates[:, 0], coordinates[:, 1], s=marker_sizes, color="#F2CF5B", edgecolor="#333333", zorder=2)
    for index, group in enumerate(groups):
        axis.text(coordinates[index, 0], coordinates[index, 1], group, ha="center", va="center", zorder=3)
    axis.set_title("PAGA connectivity abstraction (circular display layout)")
    axis.set_aspect("equal")
    axis.axis("off")
    legend_handles = []
    if any(record["in_spanning_forest"] for record in displayed):
        legend_handles.append(Line2D([0], [0], color="#4C78A8", linewidth=2, label="Spanning-forest edge"))
    if any(not record["in_spanning_forest"] for record in displayed):
        legend_handles.append(
            Line2D([0], [0], color="#9E9E9E", linewidth=2, linestyle="--", label="Other retained edge")
        )
    if legend_handles:
        axis.legend(handles=legend_handles, loc="upper right")
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    edge_details = [
        {key: value for key, value in record.items() if key not in {"left_index", "right_index"}}
        for record in displayed[:100]
    ]
    details = {
        "groupby": partition["groupby"],
        "groups": groups,
        "group_sizes": [int(value) for value in sizes],
        "min_connectivity": min_connectivity,
        "available_edges": len(available),
        "plotted_edges": len(displayed),
        "plotted_spanning_forest_edges": sum(record["in_spanning_forest"] for record in displayed),
        "plotted_edge_details": edge_details,
        "plotted_edge_details_truncated": len(displayed) > len(edge_details),
        "layout": "fixed_circular_in_partition_order",
        "edge_style_semantics": {
            "solid_blue": "stored spanning-forest edge",
            "dashed_gray": "other retained stored connectivity edge",
        },
        "connected_components": len(graph["component_sizes"]),
        "graph_fingerprint_sha256": provenance["graph_fingerprint_sha256"],
        "partition_fingerprint_sha256": provenance["partition_fingerprint_sha256"],
        "output_fingerprint_sha256": provenance["output_fingerprint_sha256"],
        "title": "PAGA connectivity abstraction",
    }
    if _return_details:
        return png, details
    return png


def run_paga_plot(adata: Any, *, min_connectivity: float = 0.0) -> tuple[bytes, dict[str, Any]]:
    return _standalone_paga_plot(adata, min_connectivity=min_connectivity, _return_details=True)


def _trajectory_plot_source(*functions: object) -> str:
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence

_TA_DIFFMAP_SCHEMA = {_TA_DIFFMAP_SCHEMA!r}
_TA_PAGA_SCHEMA = {_TA_PAGA_SCHEMA!r}
_TA_DPT_SCHEMA = {_TA_DPT_SCHEMA!r}
_TA_DIFFMAP_PROVENANCE_KEY = {_TA_DIFFMAP_PROVENANCE_KEY!r}
_TA_PAGA_PROVENANCE_KEY = {_TA_PAGA_PROVENANCE_KEY!r}
_TA_DPT_PROVENANCE_KEY = {_TA_DPT_PROVENANCE_KEY!r}

{implementation}
"""


_GRAPH_HELPERS = (
    _ta_clean_text,
    _ta_json_value,
    _ta_obs_names,
    _ta_hash_text_records,
    _ta_hash_sparse,
    _ta_hash_array,
    _ta_sparse_graph_matrix,
    _ta_graph,
)


def paga_plot_code(*, min_connectivity: float) -> str:
    source = _trajectory_plot_source(
        *_GRAPH_HELPERS,
        _ta_label_record,
        _ta_partition,
        _ta_paga_matrix,
        _ta_paga_output_fingerprint,
        _standalone_validate_paga_plot_input,
        _standalone_paga_plot,
    )
    wrapper = f"""

def plot_paga(adata):
    return _standalone_paga_plot(adata, min_connectivity={min_connectivity!r})
"""
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, "<plot_paga>", "exec")
    return code


def _standalone_validate_diffusion_plot_input(adata):
    import numpy as np
    from scipy import sparse
    from scipy.sparse import csgraph

    operation = "Diffusion Spectrum Plot"
    provenance = adata.uns.get(_TA_DIFFMAP_PROVENANCE_KEY)
    expected_fields = {
        "schema",
        "neighbors_key",
        "connectivities_key",
        "distances_key",
        "graph_fingerprint_sha256",
        "observation_fingerprint_sha256",
        "n_comps",
        "random_seed",
        "output_fingerprint_sha256",
        "scanpy_version",
    }
    if provenance is not None:
        if not isinstance(provenance, dict) or set(provenance) != expected_fields:
            raise ValueError(f"{operation} requires exact OpenBio Diffusion Map provenance when present.")
        if provenance["schema"] != _TA_DIFFMAP_SCHEMA:
            raise ValueError(f"{operation} producer schema is unsupported.")
    graph = _ta_graph(
        adata,
        neighbors_key=provenance["neighbors_key"] if provenance is not None else "neighbors",
        operation=operation,
        numpy=np,
        scipy_sparse=sparse,
        scipy_csgraph=csgraph,
    )
    if provenance is not None and (
        provenance["connectivities_key"] != graph["connectivities_key"]
        or provenance["distances_key"] != graph["distances_key"]
    ):
        raise ValueError(f"{operation} graph matrix pointers differ from the producer evidence.")
    coordinates, eigenvalues, validated_provenance, validation = _ta_validate_diffmap_bundle(
        adata,
        graph=graph,
        numpy=np,
    )
    return graph, coordinates, eigenvalues, validated_provenance, validation


def _standalone_diffusion_spectrum_plot(adata, *, _return_details=False):
    import io

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    operation = "Diffusion Spectrum Plot"
    graph, coordinates, eigenvalues, provenance, validation = _standalone_validate_diffusion_plot_input(adata)
    indices = np.arange(len(eigenvalues), dtype=int)
    figure = Figure(figsize=(8.0, 5.0), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    colors = ["#F58518", *(["#4C78A8"] * max(0, len(eigenvalues) - 1))]
    axis.bar(indices, eigenvalues, color=colors)
    axis.plot(indices, eigenvalues, color="#333333", linewidth=0.8, alpha=0.6)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xlabel("Stored diffusion component index")
    axis.set_ylabel("Transition-operator eigenvalue")
    axis.set_title("Graph-bound Diffusion Map spectrum")
    axis.grid(axis="y", alpha=0.2)
    axis.legend(
        handles=[
            Patch(color="#F58518", label="Stationary component (index 0)"),
            Patch(color="#4C78A8", label="Informative diffusion components"),
        ]
    )
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "component_indices": [int(value) for value in indices],
        "eigenvalues": [float(value) for value in eigenvalues],
        "plotted_components": len(eigenvalues),
        "informative_components": len(eigenvalues) - 1,
        "stationary_component_index": 0,
        "stationary_component_is_not_informative": True,
        "color_semantics": {
            "orange": "stationary component index 0",
            "blue": "informative diffusion components",
        },
        "coordinate_shape": [int(value) for value in coordinates.shape],
        "neighbors_key": graph["neighbors_key"],
        "graph_fingerprint_sha256": provenance["graph_fingerprint_sha256"],
        "observation_fingerprint_sha256": provenance["observation_fingerprint_sha256"],
        "output_fingerprint_sha256": provenance["output_fingerprint_sha256"],
        "scientific_validation": validation,
        "source_provenance_available": provenance.get("schema") == _TA_DIFFMAP_SCHEMA,
        "title": "Graph-bound Diffusion Map spectrum",
    }
    if _return_details:
        return png, details
    return png


def run_diffusion_spectrum_plot(adata: Any) -> tuple[bytes, dict[str, Any]]:
    return _standalone_diffusion_spectrum_plot(adata, _return_details=True)


def diffusion_spectrum_plot_code() -> str:
    source = _trajectory_plot_source(
        *_GRAPH_HELPERS,
        _ta_validate_dense_float,
        _ta_diffmap_output_fingerprint,
        _ta_transition_matrix,
        _ta_validate_diffmap_eigenpairs,
        _ta_validate_diffmap_bundle,
        _standalone_validate_diffusion_plot_input,
        _standalone_diffusion_spectrum_plot,
    )
    wrapper = """

def plot_diffusion_spectrum(adata):
    return _standalone_diffusion_spectrum_plot(adata)
"""
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, "<plot_diffusion_spectrum>", "exec")
    return code


def _standalone_validate_dpt_plot_input(adata):
    import numpy as np

    operation = "DPT Gene Trend Plot"
    graph, coordinates, _eigenvalues, diffusion_provenance, _validation = _standalone_validate_diffusion_plot_input(
        adata
    )
    provenance = adata.uns.get(_TA_DPT_PROVENANCE_KEY)
    expected_fields = {
        "schema",
        "neighbors_key",
        "graph_fingerprint_sha256",
        "observation_fingerprint_sha256",
        "diffusion_output_fingerprint_sha256",
        "n_dcs",
        "n_branchings",
        "root",
        "root_index",
        "output_fingerprint_sha256",
        "scanpy_version",
    }
    if not isinstance(provenance, dict) or set(provenance) != expected_fields:
        raise ValueError(f"{operation} requires exact OpenBio DPT provenance.")
    if provenance["schema"] != _TA_DPT_SCHEMA or provenance["n_branchings"] != 0:
        raise ValueError(f"{operation} producer identity or fixed non-branching policy is invalid.")
    if (
        provenance["neighbors_key"] != graph["neighbors_key"]
        or provenance["graph_fingerprint_sha256"] != graph["graph_fingerprint_sha256"]
        or provenance["observation_fingerprint_sha256"] != graph["observation_fingerprint_sha256"]
        or provenance["diffusion_output_fingerprint_sha256"] != diffusion_provenance["output_fingerprint_sha256"]
    ):
        raise ValueError(f"{operation} graph, observation, or diffusion provenance is stale.")
    n_dcs = provenance["n_dcs"]
    if isinstance(n_dcs, bool) or not isinstance(n_dcs, int) or not 1 <= n_dcs <= coordinates.shape[1]:
        raise ValueError(f"{operation} n_dcs does not match the graph-bound diffusion component axis.")
    root_index = provenance["root_index"]
    if isinstance(root_index, bool) or not isinstance(root_index, int) or not 0 <= root_index < int(adata.n_obs):
        raise ValueError(f"{operation} root index is invalid for the observation axis.")
    root = provenance["root"]
    if not isinstance(root, dict) or set(root) != {"mode", "root_cell_id", "rationale", "population"}:
        raise ValueError(f"{operation} root evidence has an invalid exact schema.")
    if (
        root["mode"] not in {"cell_id", "group_medoid"}
        or root["root_cell_id"] != str(adata.obs_names[root_index])
        or not isinstance(root["rationale"], str)
        or not root["rationale"].strip()
        or (root["mode"] == "cell_id" and root["population"] is not None)
        or (root["mode"] == "group_medoid" and not isinstance(root["population"], dict))
    ):
        raise ValueError(f"{operation} root evidence does not align to the selected cell.")
    if adata.uns.get("iroot") != root_index:
        raise ValueError(f"{operation} canonical iroot differs from the producer root index.")
    if "dpt_pseudotime" not in adata.obs:
        raise ValueError(f"{operation} canonical pseudotime column is missing.")
    pseudotime = adata.obs["dpt_pseudotime"].to_numpy()
    reachable = graph["component_labels"] == graph["component_labels"][root_index]
    if (
        pseudotime.shape != (int(adata.n_obs),)
        or not np.issubdtype(pseudotime.dtype, np.floating)
        or not bool(np.isfinite(pseudotime[reachable]).all())
        or not bool(np.isposinf(pseudotime[~reachable]).all())
        or bool((pseudotime[reachable] < -1e-7).any())
        or bool((pseudotime[reachable] > 1.0 + 1e-7).any())
        or not np.isclose(pseudotime[root_index], 0.0, rtol=0.0, atol=1e-7)
    ):
        raise ValueError(f"{operation} pseudotime axis is invalid or misoriented at the stored root.")
    fingerprint = _ta_hash_array(pseudotime, numpy=np)
    if fingerprint != provenance["output_fingerprint_sha256"]:
        raise ValueError(f"{operation} pseudotime fingerprint differs from the producer evidence.")
    return pseudotime.astype(float, copy=False), dict(provenance), graph, diffusion_provenance


def _trajectory_expression_fingerprint(values, *, obs_names, genes, source_kind, layer_name):
    import hashlib
    import json

    import numpy as np

    digest = hashlib.sha256()
    header = json.dumps(
        {
            "obs_names": list(obs_names),
            "genes": list(genes),
            "source_kind": source_kind,
            "layer_name": layer_name,
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    digest.update(header)
    digest.update(np.ascontiguousarray(values, dtype="<f8").tobytes(order="C"))
    return digest.hexdigest()


def _standalone_dpt_gene_trend_plot(
    adata,
    *,
    genes,
    source_kind,
    layer_name,
    n_bins=20,
    _return_details=False,
):
    import io

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from scipy import sparse

    operation = "DPT Gene Trend Plot"
    pseudotime, provenance, graph, diffusion_provenance = _standalone_validate_dpt_plot_input(adata)
    if int(adata.n_obs) > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 cells in one static binned plot.")
    if not isinstance(genes, str):
        raise TypeError(f"{operation} genes must be a comma-separated string.")
    selected_genes = [value.strip() for value in genes.split(",") if value.strip()]
    if not selected_genes or len(selected_genes) > 12 or len(set(selected_genes)) != len(selected_genes):
        raise ValueError(f"{operation} requires between 1 and 12 unique gene names.")
    if source_kind not in {"X", "raw", "layer"}:
        raise ValueError(f"{operation} expression source must be X, raw, or one named layer.")
    if source_kind == "raw":
        if adata.raw is None:
            raise ValueError(f"{operation} raw expression was selected but adata.raw is unavailable.")
        matrix = adata.raw.X
        var_names = adata.raw.var_names
    elif source_kind == "layer":
        if not isinstance(layer_name, str) or not layer_name or layer_name != layer_name.strip():
            raise ValueError(f"{operation} layer source requires one canonical layer name.")
        if layer_name not in adata.layers:
            raise ValueError(f"{operation} expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
        var_names = adata.var_names
    else:
        if layer_name is not None:
            raise ValueError(f"{operation} X source cannot carry a layer name.")
        matrix = adata.X
        var_names = adata.var_names
    if not bool(var_names.is_unique):
        raise ValueError(f"{operation} expression feature identifiers must be unique.")
    feature_names = list(var_names)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in feature_names):
        raise ValueError(f"{operation} expression feature identifiers must be canonical strings.")
    missing = [gene for gene in selected_genes if gene not in var_names]
    if missing:
        raise ValueError(f"{operation} genes were not found exactly in the selected source: {missing!r}.")
    if tuple(matrix.shape) != (int(adata.n_obs), len(feature_names)):
        raise ValueError(f"{operation} expression source does not align to the observation and feature axes.")
    positions = [int(var_names.get_loc(gene)) for gene in selected_genes]
    selected = matrix[:, positions]
    expression = np.asarray(selected.toarray() if sparse.issparse(selected) else selected, dtype=float)
    if expression.shape != (int(adata.n_obs), len(selected_genes)) or not bool(np.isfinite(expression).all()):
        raise ValueError(f"{operation} selected expression values must be finite and axis-aligned.")
    if isinstance(n_bins, bool) or not isinstance(n_bins, int) or not 1 <= n_bins <= 100:
        raise ValueError(f"{operation} n_bins must be an integer between 1 and 100.")
    reachable = np.isfinite(pseudotime)
    codes = np.full(len(pseudotime), -1, dtype=int)
    codes[reachable] = np.minimum(
        np.floor(np.clip(pseudotime[reachable], 0.0, 1.0) * n_bins).astype(int), n_bins - 1
    )
    present = [index for index in range(n_bins) if bool(np.any(codes == index))]
    plot_warnings = []
    if not bool(reachable.all()):
        plot_warnings.append(
            f"Excluded {int((~reachable).sum()):,} unreachable cells with positive-infinite pseudotime from display bins."
        )
    if len(present) == 1:
        plot_warnings.append("One display bin contains all reachable cells; no across-bin trend is shown.")
    centers = []
    counts = []
    means = []
    q1 = []
    q3 = []
    for index in present:
        mask = codes == index
        centers.append(float(pseudotime[mask].mean()))
        counts.append(int(mask.sum()))
        means.append(expression[mask].mean(axis=0))
        q1.append(np.quantile(expression[mask], 0.25, axis=0))
        q3.append(np.quantile(expression[mask], 0.75, axis=0))
    means = np.asarray(means, dtype=float)
    q1 = np.asarray(q1, dtype=float)
    q3 = np.asarray(q3, dtype=float)
    columns = min(3, len(selected_genes))
    rows = int(np.ceil(len(selected_genes) / columns))
    figure = Figure(figsize=(5.0 * columns, 3.8 * rows), constrained_layout=True)
    FigureCanvasAgg(figure)
    axes = np.atleast_1d(figure.subplots(rows, columns, squeeze=False)).reshape(-1)
    for gene_index, gene in enumerate(selected_genes):
        axis = axes[gene_index]
        axis.fill_between(centers, q1[:, gene_index], q3[:, gene_index], color="#4C78A8", alpha=0.2, label="bin IQR")
        axis.plot(centers, means[:, gene_index], color="#4C78A8", marker="o", label="bin mean")
        axis.set_xlabel("Root-dependent diffusion pseudotime")
        axis.set_ylabel("Selected expression")
        axis.set_title(gene)
        axis.grid(alpha=0.2)
        axis.legend()
    for axis in axes[len(selected_genes) :]:
        axis.set_visible(False)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    source = {"source": source_kind}
    if source_kind == "layer":
        source["layer_name"] = layer_name
    details = {
        "genes": selected_genes,
        "source": source,
        "source_gene_count": len(feature_names),
        "requested_bins": n_bins,
        "nonempty_bins": len(present),
        "bin_indices": present,
        "bin_centers": centers,
        "bin_cell_counts": counts,
        "bin_means": means.tolist(),
        "bin_q1": q1.tolist(),
        "bin_q3": q3.tolist(),
        "plotted_gene_bin_marks": int(len(present) * len(selected_genes)),
        "pseudotime_min": float(pseudotime[reachable].min()),
        "pseudotime_max": float(pseudotime[reachable].max()),
        "reachable_cells": int(reachable.sum()),
        "unreachable_cells": int((~reachable).sum()),
        "warnings": plot_warnings,
        "root_cell_id": provenance["root"]["root_cell_id"],
        "root_mode": provenance["root"]["mode"],
        "neighbors_key": graph["neighbors_key"],
        "graph_fingerprint_sha256": provenance["graph_fingerprint_sha256"],
        "diffusion_output_fingerprint_sha256": diffusion_provenance["output_fingerprint_sha256"],
        "pseudotime_fingerprint_sha256": provenance["output_fingerprint_sha256"],
        "expression_fingerprint_sha256": _trajectory_expression_fingerprint(
            expression,
            obs_names=adata.obs_names,
            genes=selected_genes,
            source_kind=source_kind,
            layer_name=layer_name,
        ),
        "title": "Expression summaries along diffusion pseudotime",
    }
    if _return_details:
        return png, details
    return png


def run_dpt_gene_trend_plot(
    adata: Any,
    *,
    genes: str,
    source_kind: str,
    layer_name: str | None,
    n_bins: int,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_dpt_gene_trend_plot(
        adata,
        genes=genes,
        source_kind=source_kind,
        layer_name=layer_name,
        n_bins=n_bins,
        _return_details=True,
    )


def dpt_gene_trend_plot_code(
    *,
    genes: str,
    source_kind: str,
    layer_name: str | None,
    n_bins: int,
) -> str:
    source = _trajectory_plot_source(
        *_GRAPH_HELPERS,
        _ta_validate_dense_float,
        _ta_diffmap_output_fingerprint,
        _ta_transition_matrix,
        _ta_validate_diffmap_eigenpairs,
        _ta_validate_diffmap_bundle,
        _standalone_validate_diffusion_plot_input,
        _standalone_validate_dpt_plot_input,
        _trajectory_expression_fingerprint,
        _standalone_dpt_gene_trend_plot,
    )
    wrapper = f"""

def plot_dpt_gene_trends(adata):
    return _standalone_dpt_gene_trend_plot(
        adata,
        genes={genes!r},
        source_kind={source_kind!r},
        layer_name={layer_name!r},
        n_bins={n_bins!r},
    )
"""
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, "<plot_dpt_gene_trends>", "exec")
    return code


__all__ = [
    "diffusion_spectrum_plot_code",
    "dpt_gene_trend_plot_code",
    "paga_plot_code",
    "run_diffusion_spectrum_plot",
    "run_dpt_gene_trend_plot",
    "run_paga_plot",
]
