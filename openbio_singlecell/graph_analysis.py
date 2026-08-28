from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from . import dependencies

if TYPE_CHECKING:
    from anndata import AnnData


@dataclass(frozen=True, slots=True)
class NamedGraph:
    neighbors_key: str
    distances_key: str | None
    connectivities_key: str
    distances: Any | None
    connectivities: Any
    n_obs: int
    obs_names: tuple[str, ...]
    positive_edges: int
    density: float
    component_labels: tuple[int, ...]
    component_sizes: tuple[int, ...]
    degree: tuple[int, ...]
    weighted_degree: tuple[float, ...]
    isolated_observations: tuple[str, ...]
    neighbor_parameters: dict[str, Any]
    undirected: bool
    symmetric: bool
    zero_diagonal: bool
    self_loop_count: int


@dataclass(frozen=True, slots=True)
class LeidenPartitionResult:
    membership: tuple[str, ...]
    modularity: float | None
    cluster_sizes: tuple[tuple[str, int], ...]
    singleton_count: int


@dataclass(frozen=True, slots=True)
class LeidenStabilityResult:
    starts: int
    seeds: tuple[int, ...]
    pairwise_ari: tuple[float, ...]
    mean_ari: float | None
    min_ari: float | None
    max_ari: float | None


def _clean_key(value: Any, *, name: str, operation: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{operation} {name} must be a string.")
    value = value.strip()
    if not value:
        raise ValueError(f"{operation} {name} cannot be empty.")
    return value


def resolve_named_graph(
    adata: AnnData,
    neighbors_key: str,
    *,
    operation: str,
    require_undirected: bool = True,
    require_zero_diagonal: bool = True,
    require_distances: bool = False,
) -> NamedGraph:
    science = dependencies.require_scientific_dependencies()
    if bool(adata.isbacked):
        raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
    if int(adata.n_obs) < 2:
        raise ValueError(f"{operation} requires at least two observations.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    neighbors_key = _clean_key(neighbors_key, name="neighbors_key", operation=operation)
    metadata = adata.uns.get(neighbors_key)
    if not isinstance(metadata, Mapping):
        raise ValueError(f"{operation} neighbor metadata not found in uns: {neighbors_key!r}")
    connectivities_key = metadata.get("connectivities_key")
    if not isinstance(connectivities_key, str) or not connectivities_key.strip():
        raise ValueError(f"{operation} uns[{neighbors_key!r}] has no valid connectivities_key.")
    connectivities_key = connectivities_key.strip()
    if connectivities_key not in adata.obsp:
        raise ValueError(f"{operation} connectivity matrix not found in obsp: {connectivities_key!r}")
    matrix = _validated_sparse_graph_matrix(
        adata.obsp[connectivities_key],
        n_obs=int(adata.n_obs),
        operation=operation,
        matrix_name="connectivity",
        require_zero_diagonal=require_zero_diagonal,
    )
    diagonal = science.np.asarray(matrix.diagonal(), dtype=float)
    zero_diagonal = bool(science.np.allclose(diagonal, 0.0, rtol=0.0, atol=1e-12))
    self_loop_count = int(science.np.count_nonzero(diagonal))
    delta = matrix - matrix.T
    delta.eliminate_zeros()
    symmetric = not delta.nnz or float(science.np.max(science.np.abs(delta.data))) <= 1e-8
    if require_undirected:
        if not symmetric:
            raise ValueError(f"{operation} requires a symmetric undirected connectivity matrix.")
        matrix = ((matrix + matrix.T) * 0.5).tocsr()
    if require_zero_diagonal:
        matrix.setdiag(0)
        matrix.eliminate_zeros()
    off_diagonal = matrix.copy()
    off_diagonal.setdiag(0)
    off_diagonal.eliminate_zeros()
    positive_edges = int(off_diagonal.nnz // 2 if require_undirected else off_diagonal.nnz)
    if positive_edges == 0:
        raise ValueError(f"{operation} connectivity graph contains no positive off-diagonal edges.")
    from scipy.sparse.csgraph import connected_components

    n_components, labels = connected_components(matrix, directed=not require_undirected, connection="weak")
    sizes = science.np.bincount(labels, minlength=n_components)
    degree = science.np.asarray((off_diagonal > 0).sum(axis=1), dtype=int).ravel()
    weighted_degree = science.np.asarray(off_diagonal.sum(axis=1), dtype=float).ravel()
    isolated_indices = science.np.flatnonzero(degree == 0)
    params = metadata.get("params")
    params = dict(params) if isinstance(params, Mapping) else {}
    distances_key: str | None = None
    distances = None
    if require_distances:
        raw_distances_key = metadata.get("distances_key")
        if not isinstance(raw_distances_key, str) or not raw_distances_key.strip():
            raise ValueError(f"{operation} uns[{neighbors_key!r}] has no valid distances_key.")
        distances_key = raw_distances_key.strip()
        if distances_key not in adata.obsp:
            raise ValueError(f"{operation} distance matrix not found in obsp: {distances_key!r}")
        distances = _validated_sparse_graph_matrix(
            adata.obsp[distances_key],
            n_obs=int(adata.n_obs),
            operation=operation,
            matrix_name="distance",
            require_zero_diagonal=True,
        )
    possible_edges = int(adata.n_obs) * (int(adata.n_obs) - 1) / (2 if require_undirected else 1)
    return NamedGraph(
        neighbors_key=neighbors_key,
        distances_key=distances_key,
        connectivities_key=connectivities_key,
        distances=distances,
        connectivities=matrix,
        n_obs=int(adata.n_obs),
        obs_names=tuple(str(value) for value in adata.obs_names),
        positive_edges=positive_edges,
        density=float(positive_edges / possible_edges) if possible_edges else 0.0,
        component_labels=tuple(int(value) for value in labels),
        component_sizes=tuple(int(value) for value in sorted(sizes, reverse=True)),
        degree=tuple(int(value) for value in degree),
        weighted_degree=tuple(float(value) for value in weighted_degree),
        isolated_observations=tuple(str(adata.obs_names[index]) for index in isolated_indices[:20]),
        neighbor_parameters=params,
        undirected=require_undirected,
        symmetric=symmetric,
        zero_diagonal=zero_diagonal,
        self_loop_count=self_loop_count,
    )


def _validated_sparse_graph_matrix(
    matrix: Any,
    *,
    n_obs: int,
    operation: str,
    matrix_name: str,
    require_zero_diagonal: bool = True,
) -> Any:
    science = dependencies.require_scientific_dependencies()
    if not science.sparse.issparse(matrix):
        raise TypeError(f"{operation} {matrix_name} matrix must be a SciPy sparse matrix.")
    if tuple(matrix.shape) != (n_obs, n_obs):
        raise ValueError(f"{operation} {matrix_name} matrix must have shape (n_obs, n_obs).")
    matrix = matrix.tocsr(copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    values = science.np.asarray(matrix.data)
    if values.size and (
        not bool(science.np.issubdtype(values.dtype, science.np.number)) or bool(science.np.iscomplexobj(values))
    ):
        raise TypeError(f"{operation} {matrix_name} matrix must contain real numeric weights.")
    if values.size and not bool(science.np.isfinite(values).all()):
        raise ValueError(f"{operation} {matrix_name} matrix contains non-finite weights.")
    if values.size and bool((values < 0).any()):
        raise ValueError(f"{operation} {matrix_name} matrix contains negative weights.")
    diagonal = science.np.asarray(matrix.diagonal(), dtype=float)
    if require_zero_diagonal and diagonal.size and not bool(
        science.np.allclose(diagonal, 0.0, rtol=0.0, atol=1e-12)
    ):
        raise ValueError(f"{operation} {matrix_name} matrix must have a zero diagonal.")
    return matrix


def graph_diagnostics(graph: NamedGraph) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()

    def summary(values: tuple[int | float, ...]) -> dict[str, int | float | None]:
        array = science.np.asarray(values, dtype=float)
        quantiles = science.np.quantile(array, [0.0, 0.25, 0.5, 0.75, 1.0])
        return {
            "n": int(array.size),
            "min": float(quantiles[0]),
            "q1": float(quantiles[1]),
            "median": float(quantiles[2]),
            "mean": float(array.mean()),
            "q3": float(quantiles[3]),
            "max": float(quantiles[4]),
        }

    isolated_count = sum(value == 0 for value in graph.degree)
    connectivity_values = tuple(float(value) for value in graph.connectivities.data)
    distance_values = tuple(float(value) for value in graph.distances.data) if graph.distances is not None else None
    return {
        "neighbors_key": graph.neighbors_key,
        "distances_key": graph.distances_key,
        "connectivities_key": graph.connectivities_key,
        "actual_keys": {
            "uns": graph.neighbors_key,
            "distances_obsp": graph.distances_key,
            "connectivities_obsp": graph.connectivities_key,
        },
        "shape": [graph.n_obs, graph.n_obs],
        "storage": "sparse_csr",
        "connectivities": {
            "shape": [graph.n_obs, graph.n_obs],
            "dtype": str(graph.connectivities.dtype),
            "nnz": int(graph.connectivities.nnz),
            "weights": summary(connectivity_values) if connectivity_values else None,
            "symmetric": graph.symmetric,
            "zero_diagonal": graph.zero_diagonal,
            "self_loop_count": graph.self_loop_count,
            "undirected_contract": graph.undirected,
        },
        "distances": (
            {
                "shape": [graph.n_obs, graph.n_obs],
                "dtype": str(graph.distances.dtype),
                "nnz": int(graph.distances.nnz),
                "weights": summary(distance_values) if distance_values else None,
                "zero_diagonal": True,
            }
            if graph.distances is not None
            else None
        ),
        "positive_undirected_edges": graph.positive_edges if graph.undirected else None,
        "positive_directed_edges": graph.positive_edges if not graph.undirected else None,
        "density": graph.density,
        "connected_components": len(graph.component_sizes),
        "component_sizes": list(graph.component_sizes),
        "isolated_observation_count": isolated_count,
        "isolated_observation_examples": list(graph.isolated_observations),
        "degree": summary(graph.degree),
        "weighted_degree": summary(graph.weighted_degree),
        "neighbor_parameters": graph.neighbor_parameters,
    }


def validate_graph_result_key(
    adata: AnnData,
    key_added: str,
    *,
    operation: str,
    obs: bool,
    uns: bool,
    overwrite_existing: bool,
) -> str:
    key_added = _clean_key(key_added, name="key_added", operation=operation)
    if not isinstance(overwrite_existing, bool):
        raise TypeError(f"{operation} overwrite_existing must be a boolean.")
    collisions = []
    if obs and key_added in adata.obs:
        collisions.append(f"obs[{key_added!r}]")
    if uns and key_added in adata.uns:
        collisions.append(f"uns[{key_added!r}]")
    if collisions and not overwrite_existing:
        raise ValueError(f"{operation} result already exists: {', '.join(collisions)}")
    return key_added


def validate_leiden_settings(
    resolution: Any,
    n_iterations: Any,
    stability_repeats: Any,
) -> tuple[float, int, int]:
    science = dependencies.require_scientific_dependencies()
    if isinstance(resolution, bool):
        raise TypeError("Leiden resolution must be a number.")
    try:
        resolution = float(resolution)
    except (TypeError, ValueError) as error:
        raise TypeError("Leiden resolution must be a number.") from error
    if not bool(science.np.isfinite(resolution)) or resolution < 0:
        raise ValueError("Leiden resolution must be finite and non-negative.")
    if isinstance(n_iterations, bool) or not isinstance(n_iterations, int):
        raise TypeError("Leiden n_iterations must be an integer.")
    if isinstance(stability_repeats, bool) or not isinstance(stability_repeats, int):
        raise TypeError("Leiden stability_repeats must be an integer.")
    if stability_repeats < 1:
        raise ValueError("Leiden stability_repeats must be positive.")
    return resolution, n_iterations, stability_repeats


def validate_random_seed(random_seed: Any, *, stability_repeats: int = 1) -> int:
    if isinstance(random_seed, bool) or not isinstance(random_seed, int):
        raise TypeError("random_seed must be an integer.")
    if isinstance(stability_repeats, bool) or not isinstance(stability_repeats, int):
        raise TypeError("Stability repeats must be an integer.")
    if stability_repeats < 1:
        raise ValueError("Stability repeats must include at least one start.")
    if random_seed < 0 or random_seed > 2**31 - stability_repeats:
        raise ValueError("random_seed range must remain between 0 and 2**31 - 1 for every start.")
    return random_seed


def run_leiden_partition(
    adata: AnnData,
    graph: NamedGraph,
    *,
    resolution: float,
    key_added: str,
    random_seed: int,
    n_iterations: int,
) -> LeidenPartitionResult:
    science = dependencies.require_scientific_dependencies()
    validate_random_seed(random_seed)
    if tuple(str(value) for value in adata.obs_names) != graph.obs_names:
        raise ValueError("Leiden AnnData observation order does not match the validated graph.")
    science.sc.tl.leiden(
        adata,
        adjacency=graph.connectivities,
        resolution=resolution,
        key_added=key_added,
        random_state=random_seed,
        flavor="igraph",
        n_iterations=n_iterations,
        directed=False,
        use_weights=True,
        objective_function="modularity",
    )
    labels = adata.obs[key_added]
    if not isinstance(labels.dtype, science.pd.CategoricalDtype):
        raise RuntimeError("Leiden backend did not return categorical memberships.")
    if len(labels) != graph.n_obs or not labels.index.equals(adata.obs_names):
        raise RuntimeError("Leiden memberships do not align to the graph observations.")
    if bool(labels.isna().any()):
        raise RuntimeError("Leiden backend returned missing memberships.")
    observed_values = tuple(str(value) for value in labels.astype(object))
    observed = tuple(dict.fromkeys(observed_values))
    mapping = {label: str(index) for index, label in enumerate(observed)}
    membership = tuple(mapping[label] for label in observed_values)
    categories = tuple(str(index) for index in range(len(observed)))
    labels = science.pd.Series(
        science.pd.Categorical(membership, categories=categories),
        index=adata.obs_names,
        name=key_added,
    )
    adata.obs[key_added] = labels
    counts = labels.value_counts(sort=False)
    cluster_sizes = tuple((str(label), int(count)) for label, count in counts.items())
    if sum(count for _, count in cluster_sizes) != graph.n_obs:
        raise RuntimeError("Leiden cluster sizes do not cover every observation exactly once.")
    if any(count <= 0 for _, count in cluster_sizes) or tuple(label for label, _ in cluster_sizes) != categories:
        raise RuntimeError("Leiden membership categories and cluster sizes are inconsistent.")
    result_metadata = adata.uns.get(key_added)
    backend_value = result_metadata.get("modularity") if isinstance(result_metadata, Mapping) else None
    if isinstance(backend_value, bool):
        raise RuntimeError("Leiden backend modularity is missing or invalid.")
    try:
        backend_modularity = float(backend_value)
    except (TypeError, ValueError) as error:
        raise RuntimeError("Leiden backend modularity is missing or invalid.") from error
    if not bool(science.np.isfinite(backend_modularity)):
        raise RuntimeError("Leiden backend modularity is non-finite.")
    import igraph

    upper = science.sparse.triu(graph.connectivities, k=1).tocoo()
    igraph_graph = igraph.Graph(
        n=graph.n_obs,
        edges=list(zip(upper.row.tolist(), upper.col.tolist(), strict=True)),
        directed=False,
    )
    igraph_graph.es["weight"] = [float(value) for value in upper.data]
    label_ids = {label: index for index, label in enumerate(dict.fromkeys(membership))}
    modularity = float(
        igraph_graph.modularity(
            [label_ids[label] for label in membership],
            weights=igraph_graph.es["weight"],
            resolution=resolution,
            directed=False,
        )
    )
    if not bool(science.np.isfinite(modularity)):
        raise RuntimeError("Leiden modularity is non-finite.")
    if not bool(science.np.isclose(backend_modularity, modularity, rtol=1e-9, atol=1e-12)):
        raise RuntimeError("Leiden backend modularity is inconsistent with the validated graph partition.")
    return LeidenPartitionResult(
        membership=membership,
        modularity=modularity,
        cluster_sizes=cluster_sizes,
        singleton_count=sum(count == 1 for _, count in cluster_sizes),
    )


def assess_leiden_stability(
    adata: AnnData,
    graph: NamedGraph,
    *,
    base_membership: tuple[str, ...],
    resolution: float,
    random_seed: int,
    n_iterations: int,
    repeats: int,
) -> LeidenStabilityResult:
    if len(base_membership) != graph.n_obs:
        raise ValueError("Leiden base membership does not align to the graph observations.")
    validate_random_seed(random_seed, stability_repeats=repeats)
    memberships = [base_membership]
    seeds = [int(random_seed)]
    for offset in range(1, repeats):
        seed = int(random_seed) + offset
        scratch = adata.copy()
        key = f"__openbio_leiden_stability_{offset}"
        collision = 0
        while key in scratch.obs or key in scratch.uns:
            collision += 1
            key = f"__openbio_leiden_stability_{offset}_{collision}"
        result = run_leiden_partition(
            scratch,
            graph,
            resolution=resolution,
            key_added=key,
            random_seed=seed,
            n_iterations=n_iterations,
        )
        memberships.append(result.membership)
        seeds.append(seed)
    pairwise = []
    if len(memberships) > 1:
        from sklearn.metrics import adjusted_rand_score

        for left in range(len(memberships)):
            for right in range(left + 1, len(memberships)):
                pairwise.append(float(adjusted_rand_score(memberships[left], memberships[right])))
    return LeidenStabilityResult(
        starts=repeats,
        seeds=tuple(seeds),
        pairwise_ari=tuple(pairwise),
        mean_ari=float(sum(pairwise) / len(pairwise)) if pairwise else None,
        min_ari=min(pairwise) if pairwise else None,
        max_ari=max(pairwise) if pairwise else None,
    )


__all__ = [
    "LeidenPartitionResult",
    "LeidenStabilityResult",
    "NamedGraph",
    "assess_leiden_stability",
    "graph_diagnostics",
    "resolve_named_graph",
    "run_leiden_partition",
    "validate_graph_result_key",
    "validate_leiden_settings",
    "validate_random_seed",
]
