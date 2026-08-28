from __future__ import annotations

import copy
import math

import pytest

from openbio_singlecell.graph_analysis import (
    resolve_named_graph,
    run_leiden_partition,
    validate_random_seed,
)
from openbio_singlecell.nodes_embedding import (
    OpenBioSingleCellLeiden,
    OpenBioSingleCellNeighbors,
    OpenBioSingleCellUMAP,
    _leiden_code,
)


def _graph_adata(science, matrix):
    adata = science.ad.AnnData(science.np.zeros((matrix.shape[0], 1)))
    adata.obs_names = [f"cell_{index}" for index in range(matrix.shape[0])]
    adata.obsp["connectivities"] = matrix
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": 2},
    }
    return adata


def test_named_graph_canonicalizes_values_inside_tolerance(science):
    matrix = science.sparse.csr_matrix(
        science.np.asarray(
            [
                [1e-13, 1.0, 0.0],
                [1.0 + 1e-9, 0.0, 0.5],
                [0.0, 0.5, 0.0],
            ]
        )
    )
    adata = _graph_adata(science, matrix)

    graph = resolve_named_graph(adata, "neighbors", operation="test")

    assert science.np.array_equal(graph.connectivities.diagonal(), science.np.zeros(3))
    assert (graph.connectivities != graph.connectivities.T).nnz == 0
    assert graph.obs_names == tuple(adata.obs_names)
    assert adata.obsp["connectivities"][0, 0] == pytest.approx(1e-13)


def test_leiden_normalizes_unused_backend_category(monkeypatch, science):
    matrix = science.sparse.csr_matrix(science.np.asarray([[0.0, 1.0], [1.0, 0.0]]))
    adata = _graph_adata(science, matrix)
    graph = resolve_named_graph(adata, "neighbors", operation="test")

    def malformed_backend(target, **kwargs):
        key = kwargs["key_added"]
        target.obs[key] = science.pd.Categorical(["0", "0"], categories=["0", "1"])
        target.uns[key] = {"params": {}, "modularity": 0.0}

    monkeypatch.setattr(science.sc.tl, "leiden", malformed_backend)
    result = run_leiden_partition(
        adata,
        graph,
        resolution=1.0,
        key_added="leiden",
        random_seed=0,
        n_iterations=-1,
    )

    assert result.membership == ("0", "0")
    assert list(adata.obs["leiden"].cat.categories) == ["0"]


@pytest.mark.parametrize("seed", [True, -1, 2**31 - 1])
def test_stability_seed_range_is_strict(seed):
    error = TypeError if seed is True else ValueError
    with pytest.raises(error):
        validate_random_seed(seed, stability_repeats=2)


def test_leiden_defaults_disclose_fast_iterations_and_stability_assessment():
    inputs = {item.id: item for item in OpenBioSingleCellLeiden.GET_SCHEMA().inputs}
    assert inputs["n_iterations"].default == 2
    assert inputs["stability_repeats"].default == 5
    assert inputs["n_iterations"].advanced is True
    assert inputs["stability_repeats"].advanced is True


def test_generated_leiden_matches_runtime_diagnostics(science):
    rng = science.np.random.default_rng(7)
    adata = science.ad.AnnData(science.np.zeros((36, 1)))
    adata.obs_names = [f"cell_{index}" for index in range(36)]
    adata.obsm["X_pca"] = rng.normal(size=(36, 5))
    neighbors = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_dimensions=5,
        n_neighbors=6,
        random_seed=4,
    ).result[0]

    node_output = OpenBioSingleCellLeiden.execute(
        neighbors,
        resolution=0.4,
        stability_repeats=3,
        random_seed=9,
    )
    runtime, report, code = node_output.result
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated = namespace["run_leiden"](neighbors)

    assert list(generated.obs["leiden"].astype(str)) == list(runtime.obs["leiden"].astype(str))
    expected = runtime.uns["leiden"]["openbio_diagnostics"]
    actual = generated.uns["leiden"]["openbio_diagnostics"]
    assert actual["cluster_sizes"] == expected["cluster_sizes"]
    assert actual["singleton_count"] == expected["singleton_count"]
    assert actual["stability"] == expected["stability"]
    sizes = science.np.asarray(list(expected["cluster_sizes"].values()), dtype=int)
    expected_fractions = {label: size / neighbors.n_obs for label, size in expected["cluster_sizes"].items()}
    assert actual["cluster_fractions"] == expected_fractions
    assert actual["min_cluster_size"] == int(sizes.min())
    assert actual["median_cluster_size"] == pytest.approx(float(science.np.median(sizes)))
    assert actual["max_cluster_size"] == int(sizes.max())
    assert actual["largest_cluster_fraction"] == pytest.approx(float(sizes.max() / neighbors.n_obs))
    assert math.isclose(actual["modularity"], expected["modularity"], rel_tol=1e-12, abs_tol=1e-12)
    summary_results = report.summary["key_results"]
    assert summary_results["cluster_sizes"] == expected["cluster_sizes"]
    assert summary_results["cluster_fractions"] == expected_fractions
    assert summary_results["min_cluster_size"] == int(sizes.min())
    assert summary_results["median_cluster_size"] == pytest.approx(float(science.np.median(sizes)))
    assert summary_results["max_cluster_size"] == int(sizes.max())
    assert summary_results["largest_cluster_fraction"] == pytest.approx(float(sizes.max() / neighbors.n_obs))
    assert any(reference["doi"] == "10.1007/BF01908075" for reference in report.summary["references"])
    assert "igraph" in report.summary["software_versions"]
    assert "pandas" in report.summary["software_versions"]
    assert "leidenalg" not in report.summary["software_versions"]
    assert resolve_named_graph(runtime, "neighbors", operation="downstream Leiden result")
    umap = OpenBioSingleCellUMAP.execute(runtime, random_seed=3).result[0]
    assert umap.obsm["X_umap"].shape == (runtime.n_obs, 2)


def test_leiden_open_parameters_and_zero_variable_graph_match_generated_code(science):
    rng = science.np.random.default_rng(71)
    adata = science.ad.AnnData(science.np.empty((30, 0)))
    adata.obs_names = [f"cell_{index}" for index in range(30)]
    adata.obsm["X_latent"] = rng.normal(size=(30, 5))
    neighbors = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_latent",
        n_neighbors=6,
        random_seed=4,
    ).result[0]

    runtime, report, code = OpenBioSingleCellLeiden.execute(
        neighbors,
        resolution=0,
        n_iterations=-7,
        stability_repeats=1,
        random_seed=9,
    ).result

    assert runtime.n_vars == 0
    assert report.summary["parameters"]["resolution"] == 0
    assert report.summary["parameters"]["n_iterations"] == -7
    assert any("resolution=0" in warning for warning in report.summary["warnings"])
    assert any("not assessed" in warning for warning in report.summary["warnings"])
    namespace: dict[str, object] = {}
    exec(code, namespace)
    with pytest.warns(UserWarning):
        generated = namespace["run_leiden"](neighbors)
    science.pd.testing.assert_series_equal(generated.obs["leiden"], runtime.obs["leiden"])
    assert generated.uns["leiden"]["openbio_diagnostics"] == runtime.uns["leiden"]["openbio_diagnostics"]


def test_leiden_accepts_finite_modularity_outside_unit_interval(monkeypatch, science):
    matrix = science.sparse.csr_matrix(science.np.asarray([[0.0, 1.0], [1.0, 0.0]]))
    adata = _graph_adata(science, matrix)
    graph = resolve_named_graph(adata, "neighbors", operation="test")
    resolution = 10.0
    backend_modularity = -5.0

    def backend(target, **kwargs):
        key = kwargs["key_added"]
        target.obs[key] = science.pd.Categorical(["left", "right"])
        target.uns[key] = {"params": {}, "modularity": backend_modularity}

    monkeypatch.setattr(science.sc.tl, "leiden", backend)
    result = run_leiden_partition(
        adata,
        graph,
        resolution=resolution,
        key_added="leiden",
        random_seed=0,
        n_iterations=2,
    )

    assert result.modularity == pytest.approx(backend_modularity)
    assert result.modularity < -1


@pytest.mark.parametrize("overwrite_existing", [False, True])
def test_runtime_and_generated_leiden_reserve_resolved_graph_metadata_key(
    science,
    overwrite_existing,
):
    rng = science.np.random.default_rng(19)
    adata = science.ad.AnnData(science.np.zeros((24, 1)))
    adata.obs_names = [f"cell_{index}" for index in range(24)]
    adata.obsm["X_pca"] = rng.normal(size=(24, 4))
    neighbors = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_dimensions=4,
        n_neighbors=5,
        random_seed=2,
    ).result[0]
    metadata_snapshot = copy.deepcopy(neighbors.uns["neighbors"])
    connectivities_snapshot = neighbors.obsp["connectivities"].copy()
    distances_snapshot = neighbors.obsp["distances"].copy()
    parameters = {
        "resolution": 0.5,
        "key_added": "neighbors",
        "neighbors_key": "neighbors",
        "n_iterations": 2,
        "stability_repeats": 1,
        "overwrite_existing": overwrite_existing,
        "random_seed": 0,
    }
    namespace: dict[str, object] = {}
    exec(_leiden_code(parameters), namespace)

    runners = (
        lambda value: OpenBioSingleCellLeiden.execute(
            value,
            resolution=0.5,
            key_added="neighbors",
            neighbors_key="neighbors",
            stability_repeats=1,
            overwrite_existing=overwrite_existing,
        ),
        namespace["run_leiden"],
    )
    for runner in runners:
        with pytest.raises(ValueError, match="cannot equal.*neighbors_key|graph metadata"):
            runner(neighbors)
        assert neighbors.uns["neighbors"] == metadata_snapshot
        connectivity_delta = neighbors.obsp["connectivities"] - connectivities_snapshot
        connectivity_delta.eliminate_zeros()
        assert connectivity_delta.nnz == 0
        distance_delta = neighbors.obsp["distances"] - distances_snapshot
        distance_delta.eliminate_zeros()
        assert distance_delta.nnz == 0
        resolved = resolve_named_graph(neighbors, "neighbors", operation="post-collision downstream")
        assert resolved.neighbors_key == "neighbors"


@pytest.mark.parametrize(
    ("malformation", "match"),
    [
        ("noncategorical", "categorical memberships"),
        ("nonfinite_modularity", "non-finite"),
        ("inconsistent_modularity", "inconsistent"),
    ],
)
def test_runtime_and_generated_leiden_reject_malformed_backend(science, monkeypatch, malformation, match):
    rng = science.np.random.default_rng(31)
    adata = science.ad.AnnData(science.np.zeros((24, 1)))
    adata.obsm["X_pca"] = rng.normal(size=(24, 4))
    neighbors = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_dimensions=4,
        n_neighbors=5,
    ).result[0]
    _, _, code = OpenBioSingleCellLeiden.execute(neighbors, resolution=0.5, stability_repeats=1).result
    namespace: dict[str, object] = {}
    exec(code, namespace)

    def malformed_backend(target, *, key_added, **kwargs):
        del kwargs
        if malformation == "noncategorical":
            target.obs[key_added] = ["0"] * target.n_obs
            modularity = 0.0
        else:
            target.obs[key_added] = science.pd.Categorical(["0"] * target.n_obs, categories=["0"])
            modularity = science.np.nan if malformation == "nonfinite_modularity" else 0.87654321
        target.uns[key_added] = {"params": {}, "modularity": modularity}

    monkeypatch.setattr(science.sc.tl, "leiden", malformed_backend)
    runners = [
        lambda value: OpenBioSingleCellLeiden.execute(value, resolution=0.5, stability_repeats=1),
        namespace["run_leiden"],
    ]
    for runner in runners:
        with pytest.raises(RuntimeError, match=match):
            runner(neighbors)


def test_neighbors_summary_discloses_complete_graph_bundle(science):
    rng = science.np.random.default_rng(41)
    adata = science.ad.AnnData(science.np.zeros((30, 1)))
    adata.obsm["X_pca"] = rng.normal(size=(30, 5))

    output, report, _ = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_dimensions=5,
        n_neighbors=6,
    ).result
    diagnostics = report.summary["key_results"]["graph"]

    assert diagnostics["actual_keys"] == {
        "uns": "neighbors",
        "distances_obsp": "distances",
        "connectivities_obsp": "connectivities",
    }
    assert diagnostics["distances"]["shape"] == [output.n_obs, output.n_obs]
    assert diagnostics["distances"]["nnz"] == output.obsp["distances"].nnz
    assert diagnostics["connectivities"]["nnz"] == output.obsp["connectivities"].nnz
    assert diagnostics["connectivities"]["symmetric"] is True


@pytest.mark.parametrize(
    ("method", "expected_dois"),
    [
        ("umap", set()),
        ("gauss", {"10.1073/pnas.0500334102", "10.1038/nmeth.3971"}),
        ("jaccard", {"10.1016/j.cell.2015.05.047"}),
    ],
)
def test_neighbors_references_match_connectivity_kernel(science, method, expected_dois):
    rng = science.np.random.default_rng(43)
    adata = science.ad.AnnData(science.np.zeros((30, 1)))
    adata.obsm["X_pca"] = rng.normal(size=(30, 5))

    _, report, _ = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_dimensions=5,
        n_neighbors=6,
        method=method,
    ).result
    actual_dois = {reference["doi"] for reference in report.summary["references"] if reference["doi"]}
    assert expected_dois <= actual_dois
    if method != "umap":
        assert "10.48550/arXiv.1802.03426" not in actual_dois


@pytest.mark.parametrize(
    ("malformation", "match"),
    [
        ("dense_distances", "distance matrix.*sparse"),
        ("nonfinite_distances", "distance matrix.*non-finite"),
        ("negative_connectivities", "connectivity matrix.*negative"),
        ("pointer_mismatch", "pointer"),
    ],
)
def test_neighbors_runtime_and_generated_code_reject_malformed_graph_bundle(science, monkeypatch, malformation, match):
    rng = science.np.random.default_rng(47)
    adata = science.ad.AnnData(science.np.zeros((24, 1)))
    adata.obsm["X_pca"] = rng.normal(size=(24, 4))
    _, _, code = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_dimensions=4,
        n_neighbors=5,
    ).result
    namespace = {}
    exec(code, namespace)

    def backend(target, *, key_added, **kwargs):
        del kwargs
        prefix = key_added or "neighbors"
        distance_key = "distances" if key_added is None else f"{key_added}_distances"
        connectivity_key = "connectivities" if key_added is None else f"{key_added}_connectivities"
        rows = science.np.arange(target.n_obs)
        columns = (rows + 1) % target.n_obs
        distances = science.sparse.csr_matrix(
            (science.np.ones(target.n_obs), (rows, columns)),
            shape=(target.n_obs, target.n_obs),
        )
        connectivities = distances.maximum(distances.T).tocsr()
        if malformation == "dense_distances":
            distances = distances.toarray()
        elif malformation == "nonfinite_distances":
            distances.data[0] = science.np.nan
        elif malformation == "negative_connectivities":
            connectivities.data[0] = -1.0
        elif malformation == "pointer_mismatch":
            distance_key = "unexpected_distances"
            connectivity_key = "unexpected_connectivities"
        target.uns[prefix] = {
            "distances_key": distance_key,
            "connectivities_key": connectivity_key,
            "params": {},
        }
        target.obsp[distance_key] = distances
        target.obsp[connectivity_key] = connectivities

    monkeypatch.setattr(science.sc.pp, "neighbors", backend)
    runners = (
        lambda value: OpenBioSingleCellNeighbors.execute(
            value,
            use_rep="X_pca",
            n_dimensions=4,
            n_neighbors=5,
        ),
        namespace["compute_neighbors"],
    )
    for runner in runners:
        with pytest.raises((TypeError, ValueError, RuntimeError), match=match):
            runner(adata)
