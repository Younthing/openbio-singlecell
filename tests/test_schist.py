from __future__ import annotations

import importlib
import inspect
import json
import os
import sys
import types

import pytest

from openbio_singlecell.node_types import AnnDataType, SummaryResultType
from openbio_singlecell.nodes_abundance import ABUNDANCE_NODE_CLASSES
from openbio_singlecell.nodes_schist import (
    SCHIST_NODE_CLASSES,
    OpenBioSingleCellSchistHierarchyPlot,
    OpenBioSingleCellSchistNestedModel,
)
from openbio_singlecell.operations_schist import schist_nested_model_owned
from openbio_singlecell.schist_analysis import run_schist_nested_model, schist_nested_model_code


def _graph_adata(science):
    obs = science.pd.DataFrame(
        {"sample": ["s1", "s1", "s2", "s2", "s3", "s3"]},
        index=[f"cell_{index}" for index in range(6)],
    )
    adata = science.ad.AnnData(
        science.np.arange(18, dtype=float).reshape(6, 3),
        obs=obs,
        var=science.pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    rows = science.np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0])
    columns = science.np.array([1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 0, 5])
    weights = science.np.array([0.8, 0.8, 0.6, 0.6, 0.9, 0.9, 0.7, 0.7, 0.5, 0.5, 0.4, 0.4])
    adata.obsp["custom_connectivities"] = science.sparse.csc_matrix(
        (weights, (rows, columns)),
        shape=(6, 6),
    )
    adata.uns["custom_neighbors"] = {
        "connectivities_key": "custom_connectivities",
        "params": {"n_neighbors": 2, "use_rep": "X_pca"},
    }
    return adata


def _install_fake_schist(science, monkeypatch, *, malformed=None):
    calls = []

    def fit_model(
        adata,
        nested=True,
        assortative=False,
        collect_marginals=True,
        n_samples=100,
        key_added=None,
        adjacency=None,
        neighbors_key="neighbors",
        constraint_key=None,
        deg_corr=True,
        directed=False,
        use_weights=False,
        bisection=True,
        simple_init=False,
        n_jobs=-1,
        n_iter=10,
        beta=1.0,
        save_model=None,
        copy=False,
        random_seed=None,
    ):
        calls.append(
            {
                "nested": nested,
                "assortative": assortative,
                "collect_marginals": collect_marginals,
                "n_samples": n_samples,
                "key_added": key_added,
                "adjacency": adjacency.copy(),
                "neighbors_key": neighbors_key,
                "constraint_key": constraint_key,
                "deg_corr": deg_corr,
                "directed": directed,
                "use_weights": use_weights,
                "bisection": bisection,
                "simple_init": simple_init,
                "n_jobs": n_jobs,
                "n_iter": n_iter,
                "beta": beta,
                "save_model": save_model,
                "copy": copy,
                "random_seed": random_seed,
            }
        )
        if malformed == "raise":
            raise RuntimeError("fake backend failed")
        memberships = [
            ["0", "0", "1", "1", "2", "2"],
            ["0", "0", "0", "0", "1", "1"],
            ["0", "0", "0", "0", "0", "0"],
        ]
        if malformed == "split_parent":
            memberships[1] = ["0", "1", "0", "1", "1", "1"]
        for level, values in enumerate(memberships):
            adata.obs[f"{key_added}_level_{level}"] = science.pd.Categorical(values)
        if malformed == "gap":
            del adata.obs[f"{key_added}_level_1"]
        elif malformed == "noncategorical":
            adata.obs[f"{key_added}_level_1"] = memberships[1]
        elif malformed == "unused_category":
            adata.obs[f"{key_added}_level_1"] = science.pd.Categorical(memberships[1], categories=["0", "1", "2"])
        elif malformed == "missing_root":
            del adata.obs[f"{key_added}_level_2"]
        marginals = [
            science.np.array(
                [
                    [0.8, 0.1, 0.0],
                    [0.7, 0.2, 0.0],
                    [0.1, 0.8, 0.0],
                    [0.0, 0.9, 0.0],
                    [0.0, 0.1, 0.8],
                    [0.0, 0.0, 0.9],
                ]
            ),
            science.np.array(
                [[0.9, 0.0], [0.9, 0.0], [0.8, 0.1], [0.9, 0.0], [0.1, 0.8], [0.0, 0.9]]
            ),
            science.np.full((6, 1), 0.9),
        ]
        for level, matrix in enumerate(marginals):
            adata.obsm[f"CM_{key_added}_level_{level}"] = matrix
        if malformed == "missing_marginal":
            del adata.obsm[f"CM_{key_added}_level_1"]
        elif malformed == "marginal_negative":
            adata.obsm[f"CM_{key_added}_level_0"][0, 0] = -0.1
        elif malformed == "marginal_over_one":
            adata.obsm[f"CM_{key_added}_level_0"][0, 0] = 1.1
        elif malformed == "marginal_zero_row":
            adata.obsm[f"CM_{key_added}_level_0"][0, :] = 0.0
        elif malformed == "marginal_shape":
            adata.obsm[f"CM_{key_added}_level_1"] = science.np.ones((6, 3)) / 3
        if malformed == "missing_root":
            del adata.obsm[f"CM_{key_added}_level_2"]
        adata.uns.setdefault("schist", {})[key_added] = {
            "stats": {
                "entropy": 12.5,
                "modularity": science.np.array([0.42, 0.21, 0.0]),
                "level_entropy": science.np.array([5.0, 2.0, 0.0]),
            },
            "blocks": {
                "0": science.np.array([0, 0, 1, 1, 2, 2], dtype=int),
                "1": science.np.array([0, 0, 1], dtype=int),
                "2": science.np.array([0, 0], dtype=int),
            },
            "params": {
                "nested": nested,
                "assortative": assortative,
                "neighbors_key": neighbors_key,
                "use_weights": use_weights,
                "key_added": key_added,
                "n_samples": n_samples,
                "collect_marginals": collect_marginals,
                "random_seed": random_seed,
                "deg_corr": deg_corr,
                "directed": directed,
                "n_iter": n_iter,
                "beta": beta,
            },
        }
        bundle = adata.uns["schist"][key_added]
        if malformed == "stats_nonfinite":
            bundle["stats"]["entropy"] = science.np.inf
        elif malformed == "stats_length":
            bundle["stats"]["modularity"] = science.np.array([0.42, 0.21])
        elif malformed == "params_mismatch":
            bundle["params"]["directed"] = True
        elif malformed == "blocks_missing":
            del bundle["blocks"]["1"]
        elif malformed == "blocks_duplicate_key_type":
            bundle["blocks"][1] = bundle["blocks"]["1"]
        elif malformed == "blocks_mismatch":
            bundle["blocks"]["1"] = science.np.array([0, 1, 1], dtype=int)
        elif malformed == "new_obsp":
            adata.obsp["schist_unexpected"] = science.sparse.eye(adata.n_obs, format="csr")
        elif malformed == "change_x":
            adata.X[0, 0] = -999.0
        elif malformed == "adjacency_mutated":
            adjacency.data[0] = adjacency.data[0] / 2
        if malformed == "return_copy":
            return adata

    assert list(inspect.signature(fit_model).parameters) == [
        "adata",
        "nested",
        "assortative",
        "collect_marginals",
        "n_samples",
        "key_added",
        "adjacency",
        "neighbors_key",
        "constraint_key",
        "deg_corr",
        "directed",
        "use_weights",
        "bisection",
        "simple_init",
        "n_jobs",
        "n_iter",
        "beta",
        "save_model",
        "copy",
        "random_seed",
    ]
    schist = types.ModuleType("schist")
    schist.__version__ = "0.10.0"
    schist.__author__ = "Davide Cittaro, Leonardo Morelli"
    schist.__file__ = "/official/dawe/schist/__init__.py"
    schist.inference = types.SimpleNamespace(fit_model=fit_model)
    monkeypatch.setitem(sys.modules, "schist", schist)

    graph_tool = types.ModuleType("graph_tool")
    graph_tool_all = types.ModuleType("graph_tool.all")
    graph_tool_all.__version__ = "3.0.1"
    graph_tool_all.__file__ = "/conda/graph_tool/all.py"
    graph_tool_all.openmp_enabled = lambda: True
    graph_tool_all.openmp_get_num_threads = lambda: 8
    graph_tool.all = graph_tool_all
    monkeypatch.setitem(sys.modules, "graph_tool", graph_tool)
    monkeypatch.setitem(sys.modules, "graph_tool.all", graph_tool_all)
    return calls


def _malform_named_graph(adata, science, case):
    key = "custom_connectivities"
    if case == "missing_metadata":
        del adata.uns["custom_neighbors"]
    elif case == "missing_pointer":
        del adata.uns["custom_neighbors"]["connectivities_key"]
    elif case == "missing_obsp":
        adata.uns["custom_neighbors"]["connectivities_key"] = "absent"
    elif case == "dense":
        adata.obsp[key] = adata.obsp[key].toarray()
    elif case == "wrong_shape":
        adata.obsp._data[key] = science.sparse.csr_matrix((5, 5))
    elif case == "negative":
        matrix = adata.obsp[key].copy()
        matrix.data[0] = -1.0
        adata.obsp[key] = matrix
    elif case == "complex":
        adata.obsp[key] = adata.obsp[key].astype(complex)
    elif case == "nonfinite":
        matrix = adata.obsp[key].copy()
        matrix.data[0] = science.np.inf
        adata.obsp[key] = matrix
    elif case == "diagonal":
        matrix = adata.obsp[key].copy()
        matrix.setdiag(0.25)
        adata.obsp[key] = matrix
    elif case == "asymmetric":
        matrix = adata.obsp[key].tolil()
        matrix[0, 1] = 0.2
        adata.obsp[key] = matrix.tocsc()
    elif case == "no_edges":
        adata.obsp[key] = science.sparse.csr_matrix((6, 6), dtype=float)
    elif case == "isolate":
        matrix = adata.obsp[key].tolil()
        matrix[5, :] = 0
        matrix[:, 5] = 0
        adata.obsp[key] = matrix.tocsr()
    elif case == "duplicate_obs":
        adata.obs_names = ["duplicate", "duplicate", "c2", "c3", "c4", "c5"]
    else:
        raise AssertionError(case)
    return adata


def test_schist_schema_is_one_atomic_nested_hierarchy_analysis():
    schema = OpenBioSingleCellSchistNestedModel.define_schema()

    assert schema.node_id == "OpenBioSingleCellSchistNestedModel"
    assert schema.display_name == "Schist Nested-SBM Hierarchy"
    assert schema.category == "openbio/single-cell/clustering"
    assert [item.id for item in schema.inputs] == [
        "adata",
        "random_seed",
        "neighbors_key",
        "key_added",
        "posterior_samples",
        "degree_correction",
        "overwrite_existing",
        "max_working_memory_gib",
    ]
    inputs = {item.id: item for item in schema.inputs}
    assert list(inspect.signature(schist_nested_model_owned).parameters) == [
        item.id for item in schema.inputs
    ]
    assert all(getattr(item, "advanced", False) for item in schema.inputs[1:])
    assert inputs["neighbors_key"].default == "neighbors"
    assert inputs["key_added"].default == "nsbm"
    assert inputs["posterior_samples"].default == 100
    assert inputs["degree_correction"].default is True
    assert inputs["overwrite_existing"].default is False
    assert inputs["max_working_memory_gib"].default == 8.0
    assert inputs["random_seed"].default == 123
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("adata", AnnDataType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_schist_registry_ownership_is_clustering_not_abundance():
    assert SCHIST_NODE_CLASSES == [OpenBioSingleCellSchistNestedModel, OpenBioSingleCellSchistHierarchyPlot]
    assert all(
        node.define_schema().node_id != "OpenBioSingleCellSchistNestedModel" for node in ABUNDANCE_NODE_CLASSES
    )


def test_schist_runs_exact_reviewed_model_and_generated_code_is_equivalent(science, monkeypatch):
    adata = _graph_adata(science)
    original = adata.copy()
    calls = _install_fake_schist(science, monkeypatch)

    output, report, code = schist_nested_model_owned(
        adata.copy(),
        random_seed=19,
        neighbors_key="custom_neighbors",
        key_added="hierarchy",
        posterior_samples=120,
        degree_correction=False,
        overwrite_existing=False,
        max_working_memory_gib=1.0,
    )

    assert len(calls) == 1
    call = calls[0]
    assert {key: value for key, value in call.items() if key != "adjacency"} == {
        "nested": True,
        "assortative": False,
        "collect_marginals": True,
        "n_samples": 120,
        "key_added": "hierarchy",
        "neighbors_key": "custom_neighbors",
        "constraint_key": None,
        "deg_corr": False,
        "directed": False,
        "use_weights": False,
        "bisection": True,
        "simple_init": False,
        "n_jobs": 1,
        "n_iter": 10,
        "beta": 1.0,
        "save_model": None,
        "copy": False,
        "random_seed": 19,
    }
    assert science.sparse.isspmatrix_csr(call["adjacency"])
    science.np.testing.assert_allclose(call["adjacency"].toarray(), call["adjacency"].toarray().T)
    assert call["adjacency"].dtype == science.np.dtype("float64")
    assert report.summary["analysis_status"] == "exploratory_graph_clustering"
    assert report.summary["key_results"]["cluster_counts_finest_to_root"] == [3, 2, 1]
    assert report.summary["key_results"]["root_level"] == 2
    assert report.summary["completion"]["backend_returned_and_postconditions_passed"] is True
    assert report.summary["completion"]["convergence_diagnostic_available"] is False
    assert report.summary["graph"]["neighbors_key"] == "custom_neighbors"
    assert len(report.summary["graph"]["fingerprint_sha256"]) == 64
    assert report.summary["parameters"]["connectivity_weights_used_by_model"] is False
    assert report.summary["software_versions"]["schist"] == "0.10.0"
    hierarchy_fingerprint = output.uns["schist"]["hierarchy"]["openbio_evidence_sha256"]
    assert len(hierarchy_fingerprint) == 64
    assert report.summary["key_results"]["hierarchy_evidence_sha256"] == hierarchy_fingerprint
    assert report.summary["software_versions"]["graph-tool"] == "3.0.1"
    assert any("differential abundance" in item.lower() for item in report.summary["limitations"])
    assert any("curated annotation" in item.lower() for item in report.summary["limitations"])
    assert "selected_level" not in json.dumps(report.summary)
    json.dumps(report.summary, allow_nan=False)
    assert "from openbio_singlecell" not in code
    assert 'backend["fit_model"](' in code
    assert "save_model=None" in code
    compile(code, "<schist-code>", "exec")

    science.np.testing.assert_array_equal(adata.X, original.X)
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    assert adata.uns.keys() == original.uns.keys()
    assert set(output.obs) == {"sample", "hierarchy_level_0", "hierarchy_level_1", "hierarchy_level_2"}
    history = output.uns["openbio_singlecell"]["analysis_history"]
    latest_history = history[sorted(history)[-1]]
    assert latest_history["operation"] == "schist_nested_sbm_hierarchy"
    assert latest_history["random_seed"] == 19

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["run_schist_nsbm"](adata)
    assert len(calls) == 2
    assert {key: value for key, value in calls[1].items() if key != "adjacency"} == {
        key: value for key, value in calls[0].items() if key != "adjacency"
    }
    science.np.testing.assert_array_equal(calls[1]["adjacency"].toarray(), calls[0]["adjacency"].toarray())
    science.pd.testing.assert_frame_equal(reproduced.obs, output.obs)
    for key in ("CM_hierarchy_level_0", "CM_hierarchy_level_1", "CM_hierarchy_level_2"):
        science.np.testing.assert_array_equal(reproduced.obsm[key], output.obsm[key])
    assert reproduced.uns["schist"]["hierarchy"].keys() == output.uns["schist"]["hierarchy"].keys()
    assert "openbio_singlecell" not in reproduced.uns
    assert report.summary["execution"]["runtime_duration_seconds"] >= 0
    assert reproduced_summary["execution"]["runtime_duration_seconds"] >= 0
    normalized_runtime = json.loads(json.dumps(report.summary))
    normalized_generated = json.loads(json.dumps(reproduced_summary))
    normalized_runtime["execution"]["runtime_duration_seconds"] = "measured"
    normalized_generated["execution"]["runtime_duration_seconds"] = "measured"
    assert normalized_generated == normalized_runtime


@pytest.mark.parametrize(
    ("overrides", "error_type", "message"),
    [
        ({"random_seed": 0}, ValueError, "zero leaves graph-tool unseeded"),
        ({"random_seed": True}, TypeError, "random_seed"),
        ({"posterior_samples": 99}, ValueError, "at least 100"),
        ({"posterior_samples": True}, TypeError, "posterior_samples"),
        ({"degree_correction": 1}, TypeError, "degree_correction"),
        ({"overwrite_existing": 1}, TypeError, "overwrite_existing"),
        ({"neighbors_key": " "}, ValueError, "neighbors_key"),
        ({"key_added": "unsafe-key"}, ValueError, "key_added"),
        ({"max_working_memory_gib": 0}, ValueError, "max_working_memory_gib"),
        ({"max_working_memory_gib": float("nan")}, ValueError, "max_working_memory_gib"),
    ],
)
def test_invalid_schist_settings_fail_before_copy_or_backend_import(
    science,
    monkeypatch,
    overrides,
    error_type,
    message,
):
    adata = _graph_adata(science)
    copied = []
    original_copy = type(adata).copy

    def forbidden_copy(self, *args, **kwargs):
        copied.append(True)
        return original_copy(self, *args, **kwargs)

    imported = []
    original_import_module = importlib.import_module

    def tracked_import(name, package=None):
        if name in {"schist", "graph_tool.all"}:
            imported.append(name)
        return original_import_module(name, package)

    monkeypatch.setattr(type(adata), "copy", forbidden_copy)
    monkeypatch.setattr(importlib, "import_module", tracked_import)
    parameters = {
        "random_seed": 19,
        "neighbors_key": "custom_neighbors",
        "key_added": "hierarchy",
        "posterior_samples": 120,
        "degree_correction": False,
        "overwrite_existing": False,
        "max_working_memory_gib": 1.0,
    }
    parameters.update(overrides)

    with pytest.raises(error_type, match=message):
        schist_nested_model_owned(adata, **parameters)

    assert copied == []
    assert imported == []


@pytest.mark.parametrize(
    ("case", "error_type", "message"),
    [
        ("missing_metadata", ValueError, "metadata not found"),
        ("missing_pointer", ValueError, "connectivities_key"),
        ("missing_obsp", ValueError, "not found in obsp"),
        ("dense", TypeError, "sparse"),
        ("wrong_shape", ValueError, "shape"),
        ("negative", ValueError, "negative"),
        ("complex", TypeError, "real numeric"),
        ("nonfinite", ValueError, "non-finite"),
        ("diagonal", ValueError, "zero diagonal"),
        ("asymmetric", ValueError, "symmetric"),
        ("no_edges", ValueError, "no positive"),
        ("isolate", ValueError, "isolated"),
        ("duplicate_obs", ValueError, "unique observation"),
    ],
)
def test_malformed_named_graph_fails_before_output_copy_or_backend(
    science,
    monkeypatch,
    case,
    error_type,
    message,
):
    adata = _malform_named_graph(_graph_adata(science), science, case)
    copied = []
    original_copy = type(adata).copy

    def forbidden_copy(self, *args, **kwargs):
        copied.append(True)
        return original_copy(self, *args, **kwargs)

    imported = []
    original_import_module = importlib.import_module

    def tracked_import(name, package=None):
        if name in {"schist", "graph_tool.all"}:
            imported.append(name)
        return original_import_module(name, package)

    monkeypatch.setattr(type(adata), "copy", forbidden_copy)
    monkeypatch.setattr(importlib, "import_module", tracked_import)

    with pytest.raises(error_type, match=message):
        schist_nested_model_owned(
            adata,
            random_seed=19,
            neighbors_key="custom_neighbors",
            key_added="hierarchy",
            posterior_samples=120,
            max_working_memory_gib=1.0,
        )

    assert copied == []
    assert imported == []


def test_equivalent_sparse_graph_formats_have_one_canonical_fingerprint(science, monkeypatch):
    base = _graph_adata(science).obsp["custom_connectivities"].tocsr()
    formats = [
        science.sparse.csr_matrix(base),
        science.sparse.csc_matrix(base),
        science.sparse.csr_array(base),
    ]
    fingerprints = []
    passed_matrices = []
    for matrix in formats:
        adata = _graph_adata(science)
        adata.obsp["custom_connectivities"] = matrix
        calls = _install_fake_schist(science, monkeypatch)
        _, report, _ = schist_nested_model_owned(
            adata,
            random_seed=19,
            neighbors_key="custom_neighbors",
            key_added="hierarchy",
            posterior_samples=120,
            max_working_memory_gib=1.0,
        )
        fingerprints.append(report.summary["graph"]["fingerprint_sha256"])
        passed_matrices.append(calls[0]["adjacency"])

    assert len(set(fingerprints)) == 1
    assert all(science.sparse.isspmatrix_csr(matrix) for matrix in passed_matrices)
    for matrix in passed_matrices[1:]:
        science.np.testing.assert_array_equal(matrix.indptr, passed_matrices[0].indptr)
        science.np.testing.assert_array_equal(matrix.indices, passed_matrices[0].indices)
        science.np.testing.assert_array_equal(matrix.data, passed_matrices[0].data)


def test_disconnected_named_graph_is_preserved_and_disclosed(science, monkeypatch):
    adata = _graph_adata(science)
    rows = [0, 1, 1, 2, 2, 0, 3, 4, 4, 5, 5, 3]
    columns = [1, 0, 2, 1, 0, 2, 4, 3, 5, 4, 3, 5]
    adata.obsp["custom_connectivities"] = science.sparse.csr_matrix(
        (science.np.ones(len(rows)), (rows, columns)),
        shape=(6, 6),
    )
    _install_fake_schist(science, monkeypatch)

    _, report, _ = schist_nested_model_owned(
        adata,
        random_seed=19,
        neighbors_key="custom_neighbors",
        key_added="hierarchy",
        posterior_samples=120,
        max_working_memory_gib=1.0,
    )

    assert report.summary["graph"]["connected_components"] == 2
    assert report.summary["graph"]["component_sizes"] == [3, 3]
    assert any("disconnected" in warning.lower() for warning in report.summary["warnings"])


def test_default_scanpy_named_graph_is_explicitly_resolved(science, monkeypatch):
    adata = _graph_adata(science)
    adata.obsp["connectivities"] = adata.obsp.pop("custom_connectivities")
    adata.uns["neighbors"] = adata.uns.pop("custom_neighbors")
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"
    calls = _install_fake_schist(science, monkeypatch)

    _, report, _ = schist_nested_model_owned(
        adata,
        random_seed=19,
        posterior_samples=120,
        max_working_memory_gib=1.0,
    )

    assert calls[0]["neighbors_key"] == "neighbors"
    assert report.summary["graph"]["connectivities_key"] == "connectivities"


def test_success_preserves_all_unowned_anndata_state(science, monkeypatch):
    adata = _graph_adata(science)
    adata.obs["quality"] = science.pd.Categorical(["good", "good", "ok", "ok", "good", "ok"])
    adata.layers["counts"] = science.np.arange(18, dtype=int).reshape(6, 3)
    adata.obsm["X_pca"] = science.np.arange(12, dtype=float).reshape(6, 2)
    adata.varm["loadings"] = science.np.arange(6, dtype=float).reshape(3, 2)
    adata.obsp["other_graph"] = science.sparse.eye(6, format="csr")
    adata.varp["gene_graph"] = science.sparse.eye(3, format="csr")
    adata.uns["provenance"] = {"method": "upstream", "values": science.np.array([1, 2, 3])}
    adata.raw = adata.copy()
    original = adata.copy()
    _install_fake_schist(science, monkeypatch)

    output, _, _ = schist_nested_model_owned(
        adata,
        random_seed=19,
        neighbors_key="custom_neighbors",
        key_added="hierarchy",
        posterior_samples=120,
        max_working_memory_gib=1.0,
    )

    science.np.testing.assert_array_equal(output.X, original.X)
    science.pd.testing.assert_frame_equal(output.obs[["sample", "quality"]], original.obs)
    science.pd.testing.assert_frame_equal(output.var, original.var)
    science.np.testing.assert_array_equal(output.layers["counts"], original.layers["counts"])
    science.np.testing.assert_array_equal(output.obsm["X_pca"], original.obsm["X_pca"])
    science.np.testing.assert_array_equal(output.varm["loadings"], original.varm["loadings"])
    science.np.testing.assert_array_equal(output.obsp["other_graph"].toarray(), original.obsp["other_graph"].toarray())
    science.np.testing.assert_array_equal(output.varp["gene_graph"].toarray(), original.varp["gene_graph"].toarray())
    science.np.testing.assert_array_equal(output.raw.X, original.raw.X)
    science.pd.testing.assert_frame_equal(output.raw.var, original.raw.var)
    assert output.uns["provenance"]["method"] == "upstream"
    science.np.testing.assert_array_equal(output.uns["provenance"]["values"], science.np.array([1, 2, 3]))


def test_too_few_or_backed_observations_fail_before_backend(science, monkeypatch, tmp_path):
    too_small = _graph_adata(science)[:1].copy()
    with pytest.raises(ValueError, match="at least two"):
        schist_nested_model_owned(too_small)

    path = tmp_path / "backed.h5ad"
    _graph_adata(science).write_h5ad(path)
    backed = science.ad.read_h5ad(path, backed="r")
    try:
        with pytest.raises(ValueError, match="in-memory"):
            schist_nested_model_owned(backed)
    finally:
        backed.file.close()


@pytest.mark.parametrize(
    "case",
    [
        "missing_schist",
        "schist_oserror",
        "wrong_pypi_identity",
        "old_version",
        "future_version",
        "missing_inference",
        "missing_fit_model",
        "signature_drift",
        "default_drift",
        "missing_graph_tool",
        "graph_tool_oserror",
    ],
)
def test_schist_dependency_gate_fails_closed_before_output_copy(science, monkeypatch, case):
    adata = _graph_adata(science)
    _install_fake_schist(science, monkeypatch)
    schist = sys.modules["schist"]
    if case == "wrong_pypi_identity":
        schist.__author__ = "Unrelated Knowledge Graph Authors"
    elif case == "old_version":
        schist.__version__ = "0.9.4"
    elif case == "future_version":
        schist.__version__ = "0.10.1"
    elif case == "missing_inference":
        del schist.inference
    elif case == "missing_fit_model":
        schist.inference = types.SimpleNamespace()
    elif case == "signature_drift":
        schist.inference.fit_model = lambda adata, **kwargs: None
    elif case == "default_drift":
        fit_model = schist.inference.fit_model
        signature = inspect.signature(fit_model)
        parameters = [
            parameter.replace(default=1) if parameter.name == "n_jobs" else parameter
            for parameter in signature.parameters.values()
        ]
        fit_model.__signature__ = signature.replace(parameters=parameters)

    original_import_module = importlib.import_module

    def gated_import(name, package=None):
        if name == "schist" and case == "missing_schist":
            raise ImportError("not installed")
        if name == "schist" and case == "schist_oserror":
            raise OSError("binary load failed")
        if name == "graph_tool.all" and case == "missing_graph_tool":
            raise ImportError("graph-tool not installed")
        if name == "graph_tool.all" and case == "graph_tool_oserror":
            raise OSError("graph-tool binary load failed")
        return original_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", gated_import)
    copied = []
    original_copy = type(adata).copy

    def forbidden_copy(self, *args, **kwargs):
        copied.append(True)
        return original_copy(self, *args, **kwargs)

    monkeypatch.setattr(type(adata), "copy", forbidden_copy)

    with pytest.raises(RuntimeError) as captured:
        schist_nested_model_owned(
            adata,
            random_seed=19,
            neighbors_key="custom_neighbors",
            key_added="hierarchy",
            posterior_samples=120,
            max_working_memory_gib=1.0,
        )

    message = str(captured.value)
    assert "https://github.com/dawe/schist" in message
    assert "WSL" in message
    assert "pip install schist" not in message
    assert copied == []


@pytest.mark.parametrize("collision", ["obs", "obsm", "uns"])
def test_output_family_collision_fails_before_copy_or_backend(science, monkeypatch, collision):
    adata = _graph_adata(science)
    if collision == "obs":
        adata.obs["hierarchy_level_9"] = science.pd.Categorical(["0"] * adata.n_obs)
    elif collision == "obsm":
        adata.obsm["CM_hierarchy_level_9"] = science.np.ones((adata.n_obs, 1))
    else:
        adata.uns["schist"] = {"hierarchy": {"stale": True}}
    copied = []
    original_copy = type(adata).copy

    def forbidden_copy(self, *args, **kwargs):
        copied.append(True)
        return original_copy(self, *args, **kwargs)

    imported = []
    original_import_module = importlib.import_module

    def tracked_import(name, package=None):
        if name in {"schist", "graph_tool.all"}:
            imported.append(name)
        return original_import_module(name, package)

    monkeypatch.setattr(type(adata), "copy", forbidden_copy)
    monkeypatch.setattr(importlib, "import_module", tracked_import)

    with pytest.raises(ValueError, match="complete owned family"):
        schist_nested_model_owned(
            adata,
            random_seed=19,
            neighbors_key="custom_neighbors",
            key_added="hierarchy",
            posterior_samples=120,
            max_working_memory_gib=1.0,
        )

    assert copied == []
    assert imported == []


def test_overwrite_replaces_complete_owned_family_on_private_copy(science, monkeypatch):
    adata = _graph_adata(science)
    adata.obs["hierarchy_level_9"] = science.pd.Categorical(["0"] * adata.n_obs)
    adata.obsm["CM_hierarchy_level_9"] = science.np.ones((adata.n_obs, 1))
    adata.uns["schist"] = {"hierarchy": {"stale": True}, "other_model": {"keep": True}}
    original = adata.copy()
    _install_fake_schist(science, monkeypatch)

    output, report, _ = schist_nested_model_owned(
        adata.copy(),
        random_seed=19,
        neighbors_key="custom_neighbors",
        key_added="hierarchy",
        posterior_samples=120,
        overwrite_existing=True,
        max_working_memory_gib=1.0,
    )

    assert "hierarchy_level_9" not in output.obs
    assert "CM_hierarchy_level_9" not in output.obsm
    assert output.uns["schist"]["other_model"] == {"keep": True}
    assert output.uns["schist"]["hierarchy"]["stats"]["entropy"] == 12.5
    assert "hierarchy_level_9" in adata.obs
    assert "CM_hierarchy_level_9" in adata.obsm
    assert adata.uns["schist"]["hierarchy"] == original.uns["schist"]["hierarchy"]
    assert any("pre-existing" in warning for warning in report.summary["warnings"])


@pytest.mark.parametrize(
    ("container", "key"),
    [
        ("obs", "hierarchy_level_bad"),
        ("obsm", "CM_hierarchy_level_bad"),
    ],
)
def test_noncanonical_vendor_prefix_collision_is_never_overwritten(science, container, key):
    adata = _graph_adata(science)
    if container == "obs":
        adata.obs[key] = "reserved"
    else:
        adata.obsm[key] = science.np.ones((adata.n_obs, 1))

    with pytest.raises(ValueError, match="noncanonical prefix"):
        schist_nested_model_owned(
            adata,
            random_seed=19,
            neighbors_key="custom_neighbors",
            key_added="hierarchy",
            posterior_samples=120,
            overwrite_existing=True,
            max_working_memory_gib=1.0,
        )


@pytest.mark.parametrize(
    ("malformed", "message"),
    [
        ("raise", "fake backend failed"),
        ("gap", "consecutive"),
        ("noncategorical", "categorical"),
        ("unused_category", "unused categories"),
        ("missing_root", "one root block"),
        ("split_parent", "nested parent mapping"),
        ("missing_marginal", "one consecutive marginal"),
        ("marginal_negative", "probabilities in"),
        ("marginal_over_one", "probabilities in"),
        ("marginal_zero_row", "positive represented mass"),
        ("marginal_shape", "marginal matrix shape"),
        ("stats_nonfinite", "total entropy"),
        ("stats_length", "modularity"),
        ("params_mismatch", "parameter metadata"),
        ("blocks_missing", "block arrays"),
        ("blocks_duplicate_key_type", "block arrays"),
        ("blocks_mismatch", "block arrays are inconsistent"),
        ("new_obsp", "outside its owned"),
        ("change_x", "outside its owned"),
        ("adjacency_mutated", "mutated the validated adjacency"),
        ("return_copy", "unexpectedly returned"),
    ],
)
def test_runtime_and_generated_code_reject_the_same_malformed_backend_state(
    science,
    monkeypatch,
    malformed,
    message,
):
    adata = _graph_adata(science)
    original = adata.copy()
    _install_fake_schist(science, monkeypatch, malformed=malformed)
    parameters = {
        "random_seed": 19,
        "neighbors_key": "custom_neighbors",
        "key_added": "hierarchy",
        "posterior_samples": 120,
        "degree_correction": False,
        "overwrite_existing": False,
        "max_working_memory_gib": 1.0,
        "openbio_version": "0.2.0",
    }
    code = schist_nested_model_code(**parameters)
    namespace = {}
    exec(code, namespace)

    with pytest.raises(RuntimeError, match=message) as runtime_error:
        run_schist_nested_model(adata.copy(), **parameters)
    with pytest.raises(type(runtime_error.value), match=message) as generated_error:
        namespace["run_schist_nsbm"](adata.copy())

    assert str(generated_error.value) == str(runtime_error.value)
    science.np.testing.assert_array_equal(adata.X, original.X)
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    science.np.testing.assert_array_equal(
        adata.obsp["custom_connectivities"].toarray(),
        original.obsp["custom_connectivities"].toarray(),
    )


def test_full_memory_preflight_fails_before_canonical_csr_copy_backend_or_anndata_copy(science, monkeypatch):
    adata = _graph_adata(science)
    copied = []
    imported = []
    original_copy = type(adata).copy
    original_import_module = importlib.import_module

    def forbidden_adata_copy(self, *args, **kwargs):
        copied.append(True)
        return original_copy(self, *args, **kwargs)

    def forbidden_csr_matrix(*args, **kwargs):
        raise AssertionError("canonical CSR was materialized before the full memory guard")

    def tracked_import(name, package=None):
        if name in {"schist", "graph_tool.all"}:
            imported.append(name)
        return original_import_module(name, package)

    monkeypatch.setattr(type(adata), "copy", forbidden_adata_copy)
    monkeypatch.setattr(science.sparse, "csr_matrix", forbidden_csr_matrix)
    monkeypatch.setattr(importlib, "import_module", tracked_import)

    with pytest.raises(MemoryError, match="components="):
        run_schist_nested_model(
            adata,
            random_seed=19,
            neighbors_key="custom_neighbors",
            key_added="hierarchy",
            posterior_samples=120,
            max_working_memory_gib=0.001,
        )

    assert copied == []
    assert imported == []


def test_real_schist_010_backend_smoke_in_declared_optional_environment(science):
    if os.environ.get("OPENBIO_RUN_SCHIST_010_SMOKE") != "1":
        pytest.skip(
            "Set OPENBIO_RUN_SCHIST_010_SMOKE=1 only in the supported environment containing official "
            "dawe/schist==0.10.0 and graph-tool."
        )
    adata = _graph_adata(science)

    output, summary = run_schist_nested_model(
        adata,
        random_seed=19,
        neighbors_key="custom_neighbors",
        key_added="hierarchy",
        posterior_samples=100,
        max_working_memory_gib=8.0,
        openbio_version="0.2.0",
    )

    assert summary["software_versions"]["schist"] == "0.10.0"
    assert summary["completion"]["backend_returned_and_postconditions_passed"] is True
    assert summary["key_results"]["hierarchy_levels"]
    assert output.obs[summary["key_results"]["output_storage"]["obs_membership_keys"][-1]].nunique() == 1
