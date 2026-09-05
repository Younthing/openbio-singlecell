from __future__ import annotations

import json
import random

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.nodes_embedding import (
    OpenBioSingleCellForceDirectedGraph as ForceDirectedGraphNode,
)
from openbio_singlecell.nodes_embedding import (
    OpenBioSingleCellPCA as PCANode,
)
from openbio_singlecell.nodes_embedding import (
    OpenBioSingleCellTSNE as TSNENode,
)
from openbio_singlecell.operations_embedding import (
    force_directed_graph,
    neighbors,
    pca,
    tsne,
    umap,
)
from openbio_singlecell.operations_input import ANNDATA_CODEC, ANNDATA_KIND
from openbio_singlecell.worker_protocol import OperationContext, ProtocolError
from tests.artifact_operation_harness import run_anndata_operation

_DEFAULT_PARAMETERS = {
    pca: {
        "n_comps": 50,
        "use_hvg": True,
        "source": None,
        "overwrite_existing": False,
        "max_output_gib": 2.0,
        "random_seed": 0,
    },
    neighbors: {
        "use_rep": "X_pca",
        "n_dimensions": 0,
        "n_neighbors": 15,
        "metric": "euclidean",
        "method": "umap",
        "key_added": "neighbors",
        "overwrite_existing": False,
        "random_seed": 0,
    },
    umap: {
        "neighbors_key": "neighbors",
        "min_dist": 0.5,
        "spread": 1.0,
        "key_added": "X_umap",
        "overwrite_existing": False,
        "random_seed": 0,
    },
    tsne: {
        "use_rep": "X_pca",
        "n_dimensions": 0,
        "perplexity": 30.0,
        "metric": "euclidean",
        "early_exaggeration": 12.0,
        "learning_rate": 1000.0,
        "key_added": "X_tsne",
        "overwrite_existing": False,
        "random_seed": 0,
    },
    force_directed_graph: {
        "layout": "fr",
        "init_mode": "random",
        "init_key": "X_draw_graph_fr",
        "neighbors_key": "neighbors",
        "key_suffix": "",
        "overwrite_existing": False,
        "random_seed": 0,
    },
}


def _run(operation, adata, **parameters):
    return run_anndata_operation(operation, adata, _DEFAULT_PARAMETERS[operation] | parameters).result


def _artifact_descriptor(root):
    return {
        "type": "artifact",
        "kind": ANNDATA_KIND,
        "codec": ANNDATA_CODEC,
        "path": str(root.resolve()),
    }


def test_pca_operation_publishes_ordered_artifacts_without_changing_input_file(tmp_path, science):
    rng = science.np.random.default_rng(501)
    adata = science.ad.AnnData(rng.normal(size=(12, 5)))
    adata.obs_names = [f"cell_{index}" for index in range(adata.n_obs)]
    adata.var_names = [f"gene_{index}" for index in range(adata.n_vars)]
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_payload = input_root / ANNDATA_PAYLOAD
    before = input_payload.read_bytes()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext(staging, "00000000-0000-4000-8000-000000000001")
    parameters = {
        "n_comps": 2,
        "use_hvg": False,
        "source": {"source": "X"},
        "overwrite_existing": False,
        "max_output_gib": 1.0,
        "random_seed": 7,
    }

    records = pca(context, {"adata": _artifact_descriptor(input_root)}, parameters)

    assert [(record["type"], record["name"]) for record in records] == [
        ("artifact", "adata"),
        ("summary", "summary"),
        ("string", "code"),
    ]
    assert records[0] == {
        "type": "artifact",
        "name": "adata",
        "kind": ANNDATA_KIND,
        "codec": ANNDATA_CODEC,
        "payload": "outputs/adata",
    }
    assert input_payload.read_bytes() == before
    assert "X_pca" not in read_anndata(input_root).obsm
    written = read_anndata(staging / records[0]["payload"])
    assert written.obsm["X_pca"].shape == (12, 2)
    json.dumps(records, allow_nan=False)

    with pytest.raises(ProtocolError, match="parameters"):
        pca(context, {"adata": _artifact_descriptor(input_root)}, parameters | {"unexpected": True})


def _embedding_input(science, *, n_obs=24, n_vars=6):
    rng = science.np.random.default_rng(101)
    adata = science.ad.AnnData(rng.normal(size=(n_obs, n_vars)))
    adata.obs_names = [f"cell_{index}" for index in range(n_obs)]
    adata.var_names = [f"gene_{index}" for index in range(n_vars)]
    adata.obsm["X_pca"] = rng.normal(size=(n_obs, 4))
    return adata


def _neighbor_input(science):
    adata = _embedding_input(science)
    return _run(
        neighbors,
        adata,
        use_rep="X_pca",
        n_dimensions=4,
        n_neighbors=5,
        random_seed=3,
    )[0]


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
            neighbors,
            {"n_neighbors": 5, "metric": "euclidean"},
            "compute_neighbors",
            "connectivities",
        ),
        (
            tsne,
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
        "__openbio_neighbors_representation" if node is neighbors else "__openbio_tsne_representation"
    )
    adata.obsm[prefix] = sentinel.copy()

    node_output = _run(node, adata, use_rep="X", n_dimensions=3, random_seed=2, **kwargs)
    runtime, code = node_output[0], node_output[2]
    assert science.np.array_equal(runtime.obsm[prefix], sentinel)
    assert result_key in (runtime.obsp if node is neighbors else runtime.obsm)

    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated = namespace[generated_name](adata)
    assert science.np.array_equal(generated.obsm[prefix], sentinel)
    assert result_key in (generated.obsp if node is neighbors else generated.obsm)
    if node is neighbors:
        difference = generated.obsp[result_key] - runtime.obsp[result_key]
        assert difference.nnz == 0 or science.np.allclose(difference.data, 0.0)
    else:
        assert science.np.allclose(generated.obsm[result_key], runtime.obsm[result_key])


def test_pca_umap_and_force_graph_generated_code_matches(science):
    rng = science.np.random.default_rng(11)
    adata = science.ad.AnnData(rng.normal(size=(30, 8)))
    adata.layers["log1p_norm"] = adata.X.copy()

    pca_output = _run(pca, adata, n_comps=4, use_hvg=False)
    pca_result, pca_code = pca_output[0], pca_output[2]
    namespace: dict[str, object] = {}
    exec(pca_code, namespace)
    generated_pca = namespace["run_pca"](adata)
    assert science.np.allclose(generated_pca.obsm["X_pca"], pca_result.obsm["X_pca"])
    assert science.np.allclose(generated_pca.varm["PCs"], pca_result.varm["PCs"])

    neighbor_output = _run(
        neighbors,
        pca_result,
        use_rep="X_pca",
        n_dimensions=4,
        n_neighbors=5,
    )[0]
    umap_output = _run(umap, neighbor_output)
    umap_result, umap_code = umap_output[0], umap_output[2]
    namespace = {}
    exec(umap_code, namespace)
    generated_umap = namespace["run_umap"](neighbor_output)
    assert science.np.allclose(generated_umap.obsm["X_umap"], umap_result.obsm["X_umap"])

    graph_output = _run(force_directed_graph, neighbor_output, layout="fr")
    graph_result, graph_code = graph_output[0], graph_output[2]
    namespace = {}
    exec(graph_code, namespace)
    generated_graph = namespace["run_force_directed_graph"](neighbor_output)
    assert science.np.allclose(generated_graph.obsm["X_draw_graph_fr"], graph_result.obsm["X_draw_graph_fr"])


def test_pca_schema_has_only_computational_bounds():
    inputs = {item.id: item for item in PCANode.define_schema().inputs}

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
    output, _, code = _run(pca, adata, n_comps=4, use_hvg=False)
    assert output.obsm["X_pca"].shape == (30, 4)

    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated = namespace["run_pca"](adata)
    assert generated.obsm["X_pca"].shape == (30, 4)


def test_pca_does_not_infer_expression_state_from_history(science):
    rng = science.np.random.default_rng(29)
    adata = science.ad.AnnData(rng.normal(size=(30, 8)))
    adata.layers["log1p_norm"] = adata.X.copy()
    _, _, code = _run(pca, adata, n_comps=4, use_hvg=False)
    scaled = adata.copy()
    metadata = ensure_metadata(scaled)
    metadata["analysis_history"] = {
        "scaled": {
            "operation": "scale_to_layer",
            "parameters": {"output_layer": "log1p_norm"},
        }
    }
    scaled.uns["openbio_singlecell"] = metadata

    _, report, code = _run(pca, scaled, n_comps=4, use_hvg=False)
    assert "expression_state" not in report.summary["parameters"]
    assert "expression_state_evidence" not in report.summary["parameters"]

    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated = namespace["run_pca"](scaled)
    assert generated.obsm["X_pca"].shape == (30, 4)


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
    _, report, code = _run(
        pca,
        small,
        n_comps=2,
        use_hvg=True,
        max_output_gib=1e-6,
    )
    assert report.summary["key_results"]["estimated_dense_output_bytes"] == 30 * 2 * 4 + 8 * 2 * 8 + 2 * 16
    namespace = {}
    exec(code, namespace)
    large = make_adata(1000)

    runners = (
        lambda value: _run(
            pca,
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
    inputs = {item.id: item for item in TSNENode.define_schema().inputs}
    assert inputs["learning_rate"].default == 1000.0


def test_umap_runs_and_warns_when_min_dist_exceeds_spread(science):
    adata = _neighbor_input(science)
    output, report, code = _run(umap, adata, min_dist=1.1, spread=1.0)

    assert output.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert any("min_dist exceeds spread" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="min_dist exceeds spread"):
        generated = namespace["run_umap"](adata)
    assert science.np.allclose(generated.obsm["X_umap"], output.obsm["X_umap"])


def test_neighbors_discloses_scanpy_effective_k_for_oversized_request(science):
    adata = _embedding_input(science)
    output, report, code = _run(
        neighbors,
        adata,
        use_rep="X_pca",
        n_neighbors=100,
        random_seed=3,
    )
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

    neighbor_output, _, neighbors_code = _run(
        neighbors,
        adata,
        use_rep="X_latent",
        n_neighbors=5,
        random_seed=2,
    )
    tsne_output, _, tsne_code = _run(
        tsne,
        adata,
        use_rep="X_latent",
        perplexity=5,
        random_seed=2,
    )
    umap_output, _, umap_code = _run(umap, neighbor_output, random_seed=2)
    force_output, _, force_code = _run(force_directed_graph, neighbor_output, random_seed=2)

    assert neighbor_output.n_vars == tsne_output.n_vars == umap_output.n_vars == force_output.n_vars == 0
    namespace = {}
    for code in (neighbors_code, tsne_code, umap_code, force_code):
        exec(code, namespace)
    generated_neighbors = namespace["compute_neighbors"](adata)
    generated_tsne = namespace["run_tsne"](adata)
    generated_umap = namespace["run_umap"](neighbor_output)
    generated_force = namespace["run_force_directed_graph"](neighbor_output)
    assert generated_neighbors.n_vars == generated_tsne.n_vars == generated_umap.n_vars == generated_force.n_vars == 0
    assert science.np.allclose(generated_tsne.obsm["X_tsne"], tsne_output.obsm["X_tsne"])
    assert science.np.allclose(generated_umap.obsm["X_umap"], umap_output.obsm["X_umap"])
    assert science.np.allclose(generated_force.obsm["X_draw_graph_fr"], force_output.obsm["X_draw_graph_fr"])


def test_umap_accepts_directed_self_loop_graph_but_force_rejects_it(science):
    adata = _neighbor_input(science)
    graph = adata.obsp["connectivities"].tolil()
    graph[0, 0] = 0.25
    graph[0, 1] = 0.0
    graph[1, 0] = 0.75
    adata.obsp["connectivities"] = graph.tocsr()

    output, report, code = _run(umap, adata, random_seed=4)
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
        _run(force_directed_graph, adata)


def test_tsne_subunit_perplexity_runs_with_warning(science):
    adata = _embedding_input(science)
    output, report, code = _run(tsne, adata, perplexity=0.5, random_seed=7)
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
    _, _, code = _run(
        neighbors,
        valid,
        use_rep="X_latent",
        n_dimensions=4,
        n_neighbors=5,
        metric=metric,
    )
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
        lambda value: _run(
            neighbors,
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


@pytest.mark.parametrize(
    ("operation", "parameters", "function_name"),
    [
        (neighbors, {"n_neighbors": 5}, "compute_neighbors"),
        (tsne, {"perplexity": 5.0}, "run_tsne"),
    ],
)
def test_cosine_zero_rows_follow_backend_convention_with_warning(science, operation, parameters, function_name):
    adata = _embedding_input(science)
    adata.obsm["X_pca"][0] = 0.0
    with pytest.warns(UserWarning, match="zero-norm"):
        output, report, code = _run(operation, adata, metric="cosine", random_seed=3, **parameters)

    assert any("zero-norm" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="zero-norm"):
        generated = namespace[function_name](adata)
    if operation is neighbors:
        science.np.testing.assert_allclose(generated.obsp["connectivities"].toarray(), output.obsp["connectivities"].toarray())
    else:
        science.np.testing.assert_allclose(generated.obsm["X_tsne"], output.obsm["X_tsne"])
    science.np.testing.assert_array_equal(adata.obsm["X_pca"][0], 0.0)


def test_umap_discloses_and_executes_fixed_hidden_optimizer_policy(science, monkeypatch):
    adata = _neighbor_input(science)
    captured = {}

    def backend(target, **kwargs):
        captured.update(kwargs)
        target.obsm["X_umap"] = science.np.zeros((target.n_obs, 2), dtype=float)
        target.uns["umap"] = {"params": {"a": 1.0, "b": 1.0, "random_state": kwargs["random_state"]}}

    monkeypatch.setattr(science.sc.tl, "umap", backend)
    _, report, _ = _run(umap, adata)

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
    _, _, code = _run(umap, adata)
    namespace = {}
    exec(code, namespace)

    def backend(target, *, key_added, **kwargs):
        del kwargs
        coordinate_key = "X_umap" if key_added is None else key_added
        parameter_key = "umap" if key_added is None else key_added
        target.obsm[coordinate_key] = _malformed_coordinates(science, target.n_obs, kind)
        target.uns[parameter_key] = {"params": {}}

    monkeypatch.setattr(science.sc.tl, "umap", backend)
    for runner in (lambda value: _run(umap, value), namespace["run_umap"]):
        with pytest.raises((TypeError, RuntimeError), match="coordinate"):
            runner(adata)


@pytest.mark.parametrize("kind", ["shape", "nonnumeric", "nonfinite"])
def test_tsne_runtime_and_generated_code_reject_malformed_backend(science, monkeypatch, kind):
    adata = _embedding_input(science)
    _, _, code = _run(tsne, adata, perplexity=5.0)
    namespace = {}
    exec(code, namespace)

    def backend(target, *, key_added, **kwargs):
        del kwargs
        target.obsm[key_added] = _malformed_coordinates(science, target.n_obs, kind)
        target.uns[key_added] = {"params": {}}

    monkeypatch.setattr(science.sc.tl, "tsne", backend)
    runners = (
        lambda value: _run(tsne, value, perplexity=5.0),
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
    _, _, code = _run(pca, adata, n_comps=2, use_hvg=True)
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
        lambda value: _run(pca, value, n_comps=2, use_hvg=True),
        namespace["run_pca"],
    )
    for runner in runners:
        with pytest.raises(RuntimeError, match=match):
            runner(adata)


def test_pca_requires_unique_variable_names_in_runtime_and_generated_code(science):
    valid = _embedding_input(science, n_obs=30, n_vars=8)
    del valid.obsm["X_pca"]
    valid.layers["log1p_norm"] = valid.X.copy()
    _, _, code = _run(pca, valid, n_comps=3, use_hvg=False)
    namespace = {}
    exec(code, namespace)
    duplicate = valid.copy()
    duplicate.var_names = ["duplicate", "duplicate", *[f"gene_{index}" for index in range(2, 8)]]

    runners = (
        lambda value: _run(pca, value, n_comps=3, use_hvg=False),
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

    _, report, _ = _run(pca, adata, n_comps=3, use_hvg=False)
    results = report.summary["key_results"]
    assert results["constant_selected_variables"] == 1
    assert results["constant_selected_variable_examples"] == ["gene_3"]
    assert any("zero variance" in warning for warning in report.summary["warnings"])


@pytest.mark.parametrize("kind", ["shape", "nonnumeric", "nonfinite"])
def test_force_graph_runtime_and_generated_code_reject_malformed_backend(science, monkeypatch, kind):
    adata = _neighbor_input(science)
    _, _, code = _run(force_directed_graph, adata, layout="fr")
    namespace = {}
    exec(code, namespace)

    def backend(target, *, layout, key_added_ext, **kwargs):
        del kwargs
        target.obsm[f"X_draw_graph_{key_added_ext or layout}"] = _malformed_coordinates(science, target.n_obs, kind)
        target.uns["draw_graph"] = {"params": {"layout": layout}}

    monkeypatch.setattr(science.sc.tl, "draw_graph", backend)
    runners = (
        lambda value: _run(force_directed_graph, value, layout="fr"),
        namespace["run_force_directed_graph"],
    )
    for runner in runners:
        with pytest.raises((TypeError, RuntimeError), match="coordinate"):
            runner(adata)


def test_force_graph_generated_code_validates_existing_initialization(science):
    adata = _neighbor_input(science)
    adata.obsm["initial"] = science.np.zeros((adata.n_obs, 2), dtype=float)
    _, _, code = _run(
        force_directed_graph,
        adata,
        layout="fr",
        init_mode="existing",
        init_key="initial",
    )
    namespace = {}
    exec(code, namespace)
    invalid = adata.copy()
    invalid.obsm["initial"] = science.np.full((adata.n_obs, 2), science.np.nan)

    runners = (
        lambda value: _run(
            force_directed_graph,
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
    _, _, code = _run(
        force_directed_graph,
        adata,
        layout="fr",
        init_mode="paga",
    )
    namespace = {}
    exec(code, namespace)
    invalid = adata.copy()
    invalid.uns["paga"]["pos"] = science.np.asarray([["bad", "0"], ["1", "1"]], dtype=object)

    runners = (
        lambda value: _run(
            force_directed_graph,
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
    _, _, code = _run(force_directed_graph, adata, layout="fr")
    namespace = {}
    exec(code, namespace)
    runners = (
        lambda: _run(force_directed_graph, adata, layout="fr"),
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
    inputs = {item.id: item for item in ForceDirectedGraphNode.define_schema().inputs}
    assert set(inputs["layout"].options) == {"fr", "kk", "fa"}
    _, report, _ = _run(force_directed_graph, _neighbor_input(science), layout="kk")
    assert any(reference["doi"] == "10.1016/0020-0190(89)90102-6" for reference in report.summary["references"])
