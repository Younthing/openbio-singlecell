from __future__ import annotations

import random

import pytest

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.nodes_embedding import (
    OpenBioSingleCellForceDirectedGraph,
    OpenBioSingleCellNeighbors,
    OpenBioSingleCellPCA,
    OpenBioSingleCellTSNE,
    OpenBioSingleCellUMAP,
)


def _embedding_input(science, *, n_obs=24, n_vars=6):
    rng = science.np.random.default_rng(101)
    adata = science.ad.AnnData(rng.normal(size=(n_obs, n_vars)))
    adata.obs_names = [f"cell_{index}" for index in range(n_obs)]
    adata.var_names = [f"gene_{index}" for index in range(n_vars)]
    adata.obsm["X_pca"] = rng.normal(size=(n_obs, 4))
    return adata


def _neighbor_input(science):
    adata = _embedding_input(science)
    return OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_dimensions=4,
        n_neighbors=5,
        random_seed=3,
    ).result[0]


def _malformed_coordinates(science, n_obs, kind):
    if kind == "shape":
        return science.np.zeros((n_obs, 3), dtype=float)
    if kind == "nonnumeric":
        return science.np.full((n_obs, 2), "invalid", dtype=object)
    values = science.np.zeros((n_obs, 2), dtype=float)
    values[0, 0] = science.np.nan
    return values


def _assert_numpy_rng_state_equal(science, expected, actual):
    assert expected[0] == actual[0]
    assert science.np.array_equal(expected[1], actual[1])
    assert expected[2:] == actual[2:]


@pytest.mark.parametrize(
    ("node", "kwargs", "generated_name", "result_key"),
    [
        (
            OpenBioSingleCellNeighbors,
            {"n_neighbors": 5, "metric": "euclidean"},
            "compute_neighbors",
            "connectivities",
        ),
        (
            OpenBioSingleCellTSNE,
            {"perplexity": 5.0, "metric": "euclidean"},
            "run_tsne",
            "X_tsne",
        ),
    ],
)
def test_x_dimension_slicing_preserves_user_temp_key(science, node, kwargs, generated_name, result_key):
    rng = science.np.random.default_rng(3)
    adata = science.ad.AnnData(rng.normal(size=(24, 6)))
    sentinel = rng.normal(size=(24, 2))
    prefix = (
        "__openbio_neighbors_representation" if node is OpenBioSingleCellNeighbors else "__openbio_tsne_representation"
    )
    adata.obsm[prefix] = sentinel.copy()

    node_output = node.execute(adata, use_rep="X", n_dimensions=3, random_seed=2, **kwargs)
    runtime, code = node_output.result[0], node_output.result[2]
    assert science.np.array_equal(runtime.obsm[prefix], sentinel)
    assert result_key in (runtime.obsp if node is OpenBioSingleCellNeighbors else runtime.obsm)

    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated = namespace[generated_name](adata)
    assert science.np.array_equal(generated.obsm[prefix], sentinel)
    assert result_key in (generated.obsp if node is OpenBioSingleCellNeighbors else generated.obsm)
    if node is OpenBioSingleCellNeighbors:
        difference = generated.obsp[result_key] - runtime.obsp[result_key]
        assert difference.nnz == 0 or science.np.allclose(difference.data, 0.0)
    else:
        assert science.np.allclose(generated.obsm[result_key], runtime.obsm[result_key])


def test_pca_umap_and_force_graph_generated_code_matches(science):
    rng = science.np.random.default_rng(11)
    adata = science.ad.AnnData(rng.normal(size=(30, 8)))
    adata.layers["log1p_norm"] = adata.X.copy()

    pca_output = OpenBioSingleCellPCA.execute(adata, n_comps=4, use_hvg=False)
    pca, pca_code = pca_output.result[0], pca_output.result[2]
    namespace: dict[str, object] = {}
    exec(pca_code, namespace)
    with pytest.warns(UserWarning, match="provenance is unknown"):
        generated_pca = namespace["run_pca"](adata)
    assert science.np.allclose(generated_pca.obsm["X_pca"], pca.obsm["X_pca"])
    assert science.np.allclose(generated_pca.varm["PCs"], pca.varm["PCs"])

    neighbors = OpenBioSingleCellNeighbors.execute(
        pca,
        use_rep="X_pca",
        n_dimensions=4,
        n_neighbors=5,
    ).result[0]
    umap_output = OpenBioSingleCellUMAP.execute(neighbors)
    umap, umap_code = umap_output.result[0], umap_output.result[2]
    namespace = {}
    exec(umap_code, namespace)
    generated_umap = namespace["run_umap"](neighbors)
    assert science.np.allclose(generated_umap.obsm["X_umap"], umap.obsm["X_umap"])

    graph_output = OpenBioSingleCellForceDirectedGraph.execute(neighbors, layout="fr")
    graph, graph_code = graph_output.result[0], graph_output.result[2]
    namespace = {}
    exec(graph_code, namespace)
    generated_graph = namespace["run_force_directed_graph"](neighbors)
    assert science.np.allclose(generated_graph.obsm["X_draw_graph_fr"], graph.obsm["X_draw_graph_fr"])


def test_pca_schema_has_only_computational_bounds():
    inputs = {item.id: item for item in OpenBioSingleCellPCA.GET_SCHEMA().inputs}

    assert inputs["n_comps"].min == 1
    assert inputs["n_comps"].max is None
    assert inputs["max_output_gib"].min is None
    assert inputs["max_output_gib"].max is None


def test_pca_sparse_preflight_does_not_densify_expression(science, monkeypatch):
    rng = science.np.random.default_rng(23)
    matrix = science.sparse.csr_matrix(rng.normal(size=(30, 8)))
    adata = science.ad.AnnData(matrix)
    adata.layers["log1p_norm"] = matrix.copy()

    def reject_toarray(*args, **kwargs):
        raise AssertionError("PCA preflight must not densify sparse expression")

    monkeypatch.setattr(science.sparse.csr_matrix, "toarray", reject_toarray)
    output, _, code = OpenBioSingleCellPCA.execute(adata, n_comps=4, use_hvg=False).result
    assert output.obsm["X_pca"].shape == (30, 4)

    namespace: dict[str, object] = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="provenance is unknown"):
        generated = namespace["run_pca"](adata)
    assert generated.obsm["X_pca"].shape == (30, 4)


def test_pca_generated_code_recognizes_scale_to_layer_provenance(science):
    rng = science.np.random.default_rng(29)
    adata = science.ad.AnnData(rng.normal(size=(30, 8)))
    adata.layers["log1p_norm"] = adata.X.copy()
    _, _, code = OpenBioSingleCellPCA.execute(adata, n_comps=4, use_hvg=False).result
    scaled = adata.copy()
    metadata = ensure_metadata(scaled)
    metadata["analysis_history"] = {
        "scaled": {
            "operation": "scale_to_layer",
            "parameters": {"output_layer": "log1p_norm"},
        }
    }
    scaled.uns["openbio_singlecell"] = metadata

    _, report, _ = OpenBioSingleCellPCA.execute(scaled, n_comps=4, use_hvg=False).result
    assert not any("provenance is unknown" in item for item in report.summary["warnings"])
    assert report.summary["parameters"]["expression_state"] == "scaled"
    assert report.summary["parameters"]["expression_state_evidence"] == "OpenBio Scale history"

    namespace: dict[str, object] = {}
    exec(code, namespace)
    assert namespace["_expression_state"](scaled) == "scaled"
    generated = namespace["run_pca"](scaled)
    assert generated.obsm["X_pca"].shape == (30, 4)


def test_pca_runtime_and_generated_code_disclose_snapshot_expression_x_counts(science):
    rng = science.np.random.default_rng(37)
    valid = science.ad.AnnData(rng.normal(size=(30, 8)))
    valid.uns["log1p"] = {"base": None}
    _, _, code = OpenBioSingleCellPCA.execute(
        valid,
        n_comps=3,
        use_hvg=False,
        source={"source": "X"},
    ).result
    namespace = {}
    exec(code, namespace)

    counts = science.ad.AnnData(rng.poisson(2.0, size=(30, 8)))
    metadata = ensure_metadata(counts)
    metadata["analysis_history"] = {
        "snapshot": {
            "operation": "snapshot_expression",
            "parameters": {"source": "X", "layer_name": "counts"},
        }
    }
    counts.uns["openbio_singlecell"] = metadata

    assert namespace["_expression_state"](counts) == "counts"
    output, report, _ = OpenBioSingleCellPCA.execute(
        counts,
        n_comps=3,
        use_hvg=False,
        source={"source": "X"},
    ).result
    assert output.obsm["X_pca"].shape == (30, 3)
    assert any("explicitly selected 'counts' expression" in warning for warning in report.summary["warnings"])
    with pytest.warns(UserWarning, match="explicitly selected 'counts' expression"):
        generated = namespace["run_pca"](counts)
    assert science.np.allclose(generated.obsm["X_pca"], output.obsm["X_pca"])


def test_pca_output_budget_accounts_for_full_axis_loadings_with_hvg(science):
    rng = science.np.random.default_rng(39)

    def make_adata(n_vars):
        adata = science.ad.AnnData(rng.normal(size=(30, n_vars)))
        adata.layers["log1p_norm"] = adata.X.copy()
        adata.var["highly_variable"] = [True, True, True, *([False] * (n_vars - 3))]
        metadata = ensure_metadata(adata)
        metadata["analysis_history"] = {
            "logged": {
                "operation": "normalize_to_layer",
                "parameters": {"output_layer": "log1p_norm", "transform": "log1p"},
            }
        }
        adata.uns["openbio_singlecell"] = metadata
        return adata

    small = make_adata(8)
    _, report, code = OpenBioSingleCellPCA.execute(
        small,
        n_comps=2,
        use_hvg=True,
        max_output_gib=1e-6,
    ).result
    assert report.summary["key_results"]["estimated_dense_output_bytes"] == 30 * 2 * 4 + 8 * 2 * 8 + 2 * 16
    namespace = {}
    exec(code, namespace)
    large = make_adata(1000)

    runners = (
        lambda value: OpenBioSingleCellPCA.execute(
            value,
            n_comps=2,
            use_hvg=True,
            max_output_gib=1e-6,
        ),
        namespace["run_pca"],
    )
    for runner in runners:
        with pytest.raises(ValueError, match="dense output estimate"):
            runner(large)


def test_tsne_schema_uses_audited_scanpy_learning_rate():
    inputs = {item.id: item for item in OpenBioSingleCellTSNE.GET_SCHEMA().inputs}
    assert inputs["learning_rate"].default == 1000.0


def test_umap_runs_and_warns_when_min_dist_exceeds_spread(science):
    adata = _neighbor_input(science)
    output, report, code = OpenBioSingleCellUMAP.execute(adata, min_dist=1.1, spread=1.0).result

    assert output.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert any("min_dist exceeds spread" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="min_dist exceeds spread"):
        generated = namespace["run_umap"](adata)
    assert science.np.allclose(generated.obsm["X_umap"], output.obsm["X_umap"])


def test_neighbors_discloses_scanpy_effective_k_for_oversized_request(science):
    adata = _embedding_input(science)
    output, report, code = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_pca",
        n_neighbors=100,
        random_seed=3,
    ).result
    summary = report.summary
    effective = output.uns["neighbors"]["params"]["n_neighbors"]

    assert summary["parameters"]["n_neighbors"] == 100
    assert summary["parameters"]["effective_n_neighbors"] == effective
    assert summary["key_results"]["requested_n_neighbors"] == 100
    assert summary["key_results"]["effective_n_neighbors"] == effective
    assert effective != 100
    assert any("resolved to effective" in warning for warning in summary["warnings"])

    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="resolved to effective"):
        generated = namespace["compute_neighbors"](adata)
    assert generated.uns["neighbors"]["params"]["n_neighbors"] == effective


def test_embedding_and_graph_nodes_support_zero_variable_adata_with_obsm(science):
    rng = science.np.random.default_rng(61)
    adata = science.ad.AnnData(science.np.empty((24, 0)))
    adata.obs_names = [f"cell_{index}" for index in range(24)]
    adata.obsm["X_latent"] = rng.normal(size=(24, 4))

    neighbors, _, neighbors_code = OpenBioSingleCellNeighbors.execute(
        adata,
        use_rep="X_latent",
        n_neighbors=5,
        random_seed=2,
    ).result
    tsne, _, tsne_code = OpenBioSingleCellTSNE.execute(
        adata,
        use_rep="X_latent",
        perplexity=5,
        random_seed=2,
    ).result
    umap, _, umap_code = OpenBioSingleCellUMAP.execute(neighbors, random_seed=2).result
    force, _, force_code = OpenBioSingleCellForceDirectedGraph.execute(neighbors, random_seed=2).result

    assert neighbors.n_vars == tsne.n_vars == umap.n_vars == force.n_vars == 0
    namespace = {}
    for code in (neighbors_code, tsne_code, umap_code, force_code):
        exec(code, namespace)
    generated_neighbors = namespace["compute_neighbors"](adata)
    generated_tsne = namespace["run_tsne"](adata)
    generated_umap = namespace["run_umap"](neighbors)
    generated_force = namespace["run_force_directed_graph"](neighbors)
    assert generated_neighbors.n_vars == generated_tsne.n_vars == generated_umap.n_vars == generated_force.n_vars == 0
    assert science.np.allclose(generated_tsne.obsm["X_tsne"], tsne.obsm["X_tsne"])
    assert science.np.allclose(generated_umap.obsm["X_umap"], umap.obsm["X_umap"])
    assert science.np.allclose(generated_force.obsm["X_draw_graph_fr"], force.obsm["X_draw_graph_fr"])


def test_umap_accepts_directed_self_loop_graph_but_force_rejects_it(science):
    adata = _neighbor_input(science)
    graph = adata.obsp["connectivities"].tolil()
    graph[0, 0] = 0.25
    graph[0, 1] = 0.0
    graph[1, 0] = 0.75
    adata.obsp["connectivities"] = graph.tocsr()

    output, report, code = OpenBioSingleCellUMAP.execute(adata, random_seed=4).result
    diagnostics = report.summary["key_results"]["graph"]["connectivities"]
    assert output.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert diagnostics["symmetric"] is False
    assert diagnostics["self_loop_count"] == 1
    assert any("asymmetric" in warning for warning in report.summary["warnings"])
    assert any("self-loop" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(code, namespace)
    assert namespace["run_umap"](adata).obsm["X_umap"].shape == (adata.n_obs, 2)

    with pytest.raises(ValueError, match="zero diagonal|symmetric|undirected"):
        OpenBioSingleCellForceDirectedGraph.execute(adata)


def test_tsne_subunit_perplexity_runs_with_warning(science):
    adata = _embedding_input(science)
    output, report, code = OpenBioSingleCellTSNE.execute(adata, perplexity=0.5, random_seed=7).result
    assert output.obsm["X_tsne"].shape == (adata.n_obs, 2)
    assert any("below 1" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="below 1"):
        generated = namespace["run_tsne"](adata)
    assert science.np.allclose(generated.obsm["X_tsne"], output.obsm["X_tsne"])


@pytest.mark.parametrize(
    ("metric", "malformation", "match"),
    [
        ("correlation", "constant", "constant observation"),
        ("cosine", "zero", "zero-norm observation"),
    ],
)
@pytest.mark.parametrize("sparse_input", [False, True])
def test_neighbors_metric_degeneracy_is_rejected_before_backend(
    science,
    monkeypatch,
    metric,
    malformation,
    match,
    sparse_input,
):
    rng = science.np.random.default_rng(53)
    valid = science.ad.AnnData(science.np.zeros((24, 1)))
    representation = rng.normal(size=(24, 4))
    valid.obsm["X_latent"] = representation
    _, _, code = OpenBioSingleCellNeighbors.execute(
        valid,
        use_rep="X_latent",
        n_dimensions=4,
        n_neighbors=5,
        metric=metric,
    ).result
    namespace = {}
    exec(code, namespace)
    invalid = valid.copy()
    malformed = representation.copy()
    malformed[0] = 2.0 if malformation == "constant" else 0.0
    invalid.obsm["X_latent"] = science.sparse.csr_matrix(malformed) if sparse_input else malformed
    backend_calls = 0

    def backend(*args, **kwargs):
        nonlocal backend_calls
        del args, kwargs
        backend_calls += 1
        raise AssertionError("Neighbors backend must not run for a degenerate metric row.")

    monkeypatch.setattr(science.sc.pp, "neighbors", backend)
    runners = (
        lambda value: OpenBioSingleCellNeighbors.execute(
            value,
            use_rep="X_latent",
            n_dimensions=4,
            n_neighbors=5,
            metric=metric,
        ),
        namespace["compute_neighbors"],
    )
    for runner in runners:
        with pytest.raises(ValueError, match=match):
            runner(invalid)
    assert backend_calls == 0


def test_umap_discloses_and_executes_fixed_hidden_optimizer_policy(science, monkeypatch):
    adata = _neighbor_input(science)
    captured = {}

    def backend(target, **kwargs):
        captured.update(kwargs)
        target.obsm["X_umap"] = science.np.zeros((target.n_obs, 2), dtype=float)
        target.uns["umap"] = {"params": {"a": 1.0, "b": 1.0, "random_state": kwargs["random_state"]}}

    monkeypatch.setattr(science.sc.tl, "umap", backend)
    _, report, _ = OpenBioSingleCellUMAP.execute(adata).result

    assert captured["maxiter"] is None
    assert captured["alpha"] == 1.0
    assert captured["gamma"] == 1.0
    assert captured["negative_sample_rate"] == 5
    assert captured["a"] is None
    assert captured["b"] is None
    assert report.summary["parameters"]["init_pos"] == "spectral"
    assert report.summary["parameters"]["negative_sample_rate"] == 5


@pytest.mark.parametrize("kind", ["shape", "nonnumeric", "nonfinite"])
def test_umap_runtime_and_generated_code_reject_malformed_backend(science, monkeypatch, kind):
    adata = _neighbor_input(science)
    _, _, code = OpenBioSingleCellUMAP.execute(adata).result
    namespace = {}
    exec(code, namespace)

    def backend(target, *, key_added, **kwargs):
        del kwargs
        coordinate_key = "X_umap" if key_added is None else key_added
        parameter_key = "umap" if key_added is None else key_added
        target.obsm[coordinate_key] = _malformed_coordinates(science, target.n_obs, kind)
        target.uns[parameter_key] = {"params": {}}

    monkeypatch.setattr(science.sc.tl, "umap", backend)
    for runner in (OpenBioSingleCellUMAP.execute, namespace["run_umap"]):
        with pytest.raises((TypeError, RuntimeError), match="coordinate"):
            runner(adata)


@pytest.mark.parametrize("kind", ["shape", "nonnumeric", "nonfinite"])
def test_tsne_runtime_and_generated_code_reject_malformed_backend(science, monkeypatch, kind):
    adata = _embedding_input(science)
    _, _, code = OpenBioSingleCellTSNE.execute(adata, perplexity=5.0).result
    namespace = {}
    exec(code, namespace)

    def backend(target, *, key_added, **kwargs):
        del kwargs
        target.obsm[key_added] = _malformed_coordinates(science, target.n_obs, kind)
        target.uns[key_added] = {"params": {}}

    monkeypatch.setattr(science.sc.tl, "tsne", backend)
    runners = (
        lambda value: OpenBioSingleCellTSNE.execute(value, perplexity=5.0),
        namespace["run_tsne"],
    )
    for runner in runners:
        with pytest.raises((TypeError, RuntimeError), match="coordinate"):
            runner(adata)


@pytest.mark.parametrize(
    ("malformation", "match"),
    [
        ("score_shape", "score"),
        ("nonfinite_loadings", "loading"),
        ("integer_variance", "variance"),
        ("excluded_loading", "excluded"),
    ],
)
def test_pca_runtime_and_generated_code_reject_malformed_backend(science, monkeypatch, malformation, match):
    adata = _embedding_input(science, n_obs=30, n_vars=8)
    del adata.obsm["X_pca"]
    adata.layers["log1p_norm"] = adata.X.copy()
    adata.var["highly_variable"] = [True, True, True, True, False, False, False, False]
    metadata = ensure_metadata(adata)
    metadata["analysis_history"] = {
        "scaled": {
            "operation": "scale_to_layer",
            "parameters": {"output_layer": "log1p_norm"},
        }
    }
    adata.uns["openbio_singlecell"] = metadata
    _, _, code = OpenBioSingleCellPCA.execute(adata, n_comps=2, use_hvg=True).result
    namespace = {}
    exec(code, namespace)

    def backend(target, *, n_comps, **kwargs):
        del kwargs
        scores = science.np.zeros((target.n_obs, n_comps), dtype=science.np.float32)
        loadings = science.np.zeros((target.n_vars, n_comps), dtype=float)
        variance = science.np.asarray([2.0, 1.0], dtype=float)
        if malformation == "score_shape":
            scores = science.np.zeros((target.n_obs, n_comps + 1), dtype=science.np.float32)
        elif malformation == "nonfinite_loadings":
            loadings[0, 0] = science.np.inf
        elif malformation == "integer_variance":
            variance = science.np.asarray([2, 1], dtype=int)
        elif malformation == "excluded_loading":
            loadings[4, 0] = 1.0
        target.obsm["X_pca"] = scores
        target.varm["PCs"] = loadings
        target.uns["pca"] = {
            "params": {},
            "variance": variance,
            "variance_ratio": science.np.asarray([0.5, 0.25], dtype=float),
        }

    monkeypatch.setattr(science.sc.pp, "pca", backend)
    runners = (
        lambda value: OpenBioSingleCellPCA.execute(value, n_comps=2, use_hvg=True),
        namespace["run_pca"],
    )
    for runner in runners:
        with pytest.raises(RuntimeError, match=match):
            runner(adata)


def test_pca_requires_unique_variable_names_in_runtime_and_generated_code(science):
    valid = _embedding_input(science, n_obs=30, n_vars=8)
    del valid.obsm["X_pca"]
    valid.layers["log1p_norm"] = valid.X.copy()
    _, _, code = OpenBioSingleCellPCA.execute(valid, n_comps=3, use_hvg=False).result
    namespace = {}
    exec(code, namespace)
    duplicate = valid.copy()
    duplicate.var_names = ["duplicate", "duplicate", *[f"gene_{index}" for index in range(2, 8)]]

    runners = (
        lambda value: OpenBioSingleCellPCA.execute(value, n_comps=3, use_hvg=False),
        namespace["run_pca"],
    )
    for runner in runners:
        with pytest.raises(ValueError, match="unique.*variable"):
            runner(duplicate)


def test_pca_reports_partial_constant_variables(science):
    rng = science.np.random.default_rng(17)
    matrix = rng.normal(size=(30, 8))
    matrix[:, 3] = 7.0
    adata = science.ad.AnnData(matrix)
    adata.var_names = [f"gene_{index}" for index in range(8)]
    adata.layers["log1p_norm"] = matrix.copy()

    _, report, _ = OpenBioSingleCellPCA.execute(adata, n_comps=3, use_hvg=False).result
    results = report.summary["key_results"]
    assert results["constant_selected_variables"] == 1
    assert results["constant_selected_variable_examples"] == ["gene_3"]
    assert any("zero variance" in warning for warning in report.summary["warnings"])


@pytest.mark.parametrize("kind", ["shape", "nonnumeric", "nonfinite"])
def test_force_graph_runtime_and_generated_code_reject_malformed_backend(science, monkeypatch, kind):
    adata = _neighbor_input(science)
    _, _, code = OpenBioSingleCellForceDirectedGraph.execute(adata, layout="fr").result
    namespace = {}
    exec(code, namespace)

    def backend(target, *, layout, key_added_ext, **kwargs):
        del kwargs
        target.obsm[f"X_draw_graph_{key_added_ext or layout}"] = _malformed_coordinates(science, target.n_obs, kind)
        target.uns["draw_graph"] = {"params": {"layout": layout}}

    monkeypatch.setattr(science.sc.tl, "draw_graph", backend)
    runners = (
        lambda value: OpenBioSingleCellForceDirectedGraph.execute(value, layout="fr"),
        namespace["run_force_directed_graph"],
    )
    for runner in runners:
        with pytest.raises((TypeError, RuntimeError), match="coordinate"):
            runner(adata)


def test_force_graph_generated_code_validates_existing_initialization(science):
    adata = _neighbor_input(science)
    adata.obsm["initial"] = science.np.zeros((adata.n_obs, 2), dtype=float)
    _, _, code = OpenBioSingleCellForceDirectedGraph.execute(
        adata,
        layout="fr",
        init_mode="existing",
        init_key="initial",
    ).result
    namespace = {}
    exec(code, namespace)
    invalid = adata.copy()
    invalid.obsm["initial"] = science.np.full((adata.n_obs, 2), science.np.nan)

    runners = (
        lambda value: OpenBioSingleCellForceDirectedGraph.execute(
            value,
            layout="fr",
            init_mode="existing",
            init_key="initial",
        ),
        namespace["run_force_directed_graph"],
    )
    for runner in runners:
        with pytest.raises(ValueError, match="finite"):
            runner(invalid)


def test_force_graph_generated_code_reproduces_paga_initialization_validation(science):
    adata = _neighbor_input(science)
    adata.obs["clusters"] = science.pd.Categorical(
        ["a" if index < adata.n_obs // 2 else "b" for index in range(adata.n_obs)],
        categories=["a", "b"],
    )
    adata.uns["paga"] = {
        "groups": "clusters",
        "pos": science.np.asarray([[0.0, 0.0], [1.0, 1.0]]),
        "connectivities": science.sparse.csr_matrix([[0.0, 1.0], [1.0, 0.0]]),
    }
    _, _, code = OpenBioSingleCellForceDirectedGraph.execute(
        adata,
        layout="fr",
        init_mode="paga",
    ).result
    namespace = {}
    exec(code, namespace)
    invalid = adata.copy()
    invalid.uns["paga"]["pos"] = science.np.asarray([["bad", "0"], ["1", "1"]], dtype=object)

    runners = (
        lambda value: OpenBioSingleCellForceDirectedGraph.execute(
            value,
            layout="fr",
            init_mode="paga",
        ),
        namespace["run_force_directed_graph"],
    )
    for runner in runners:
        with pytest.raises(TypeError, match="real numeric"):
            runner(invalid)


def test_force_graph_runtime_and_generated_code_restore_global_rng(science):
    adata = _neighbor_input(science)
    _, _, code = OpenBioSingleCellForceDirectedGraph.execute(adata, layout="fr").result
    namespace = {}
    exec(code, namespace)
    runners = (
        lambda: OpenBioSingleCellForceDirectedGraph.execute(adata, layout="fr"),
        lambda: namespace["run_force_directed_graph"](adata),
    )
    for runner in runners:
        science.np.random.seed(31415)
        random.seed(92653)
        numpy_before = science.np.random.get_state()
        python_before = random.getstate()
        runner()
        _assert_numpy_rng_state_equal(science, numpy_before, science.np.random.get_state())
        assert random.getstate() == python_before


def test_force_graph_layout_interface_and_kamada_kawai_reference(science):
    inputs = {item.id: item for item in OpenBioSingleCellForceDirectedGraph.GET_SCHEMA().inputs}
    assert set(inputs["layout"].options) == {"fr", "kk", "fa"}
    _, report, _ = OpenBioSingleCellForceDirectedGraph.execute(_neighbor_input(science), layout="kk").result
    assert any(reference["doi"] == "10.1016/0020-0190(89)90102-6" for reference in report.summary["references"])
