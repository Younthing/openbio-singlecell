from __future__ import annotations

import json

import pytest

import openbio_singlecell.trajectory_analysis as trajectory_analysis
from openbio_singlecell.nodes_trajectory import (
    OpenBioSingleCellDiffusionMap,
    OpenBioSingleCellDPT,
    OpenBioSingleCellPAGA,
)
from openbio_singlecell.trajectory_analysis import (
    _ta_select_root,
    analyze_diffusion_map,
    analyze_dpt,
    analyze_paga,
)


def _trajectory_input(science, *, n_obs=36, neighbors_key="neighbors"):
    rng = science.np.random.default_rng(73)
    adata = science.ad.AnnData(rng.normal(size=(n_obs, 7)))
    adata.obs_names = [f"cell_{index:03d}" for index in range(n_obs)]
    adata.var_names = [f"gene_{index}" for index in range(7)]
    thirds = n_obs // 3
    labels = ["early"] * thirds + ["middle"] * thirds + ["late"] * (n_obs - 2 * thirds)
    adata.obs["state"] = science.pd.Categorical(
        labels,
        categories=["early", "unused", "middle", "late"],
        ordered=True,
    )
    science.sc.pp.neighbors(
        adata,
        n_neighbors=6,
        use_rep="X",
        metric="euclidean",
        method="umap",
        random_state=11,
        key_added=None if neighbors_key == "neighbors" else neighbors_key,
    )
    return adata


def _assert_rng_state_equal(science, expected, actual):
    assert expected[0] == actual[0]
    science.np.testing.assert_array_equal(expected[1], actual[1])
    assert expected[2:] == actual[2:]


def _assert_sparse_equal(science, left, right):
    difference = left.tocsr() - right.tocsr()
    difference.eliminate_zeros()
    assert difference.nnz == 0 or science.np.allclose(difference.data, 0.0, rtol=0.0, atol=1e-12)


def _tamper_named_connectivities_with_zero_average(adata, *, neighbors_key="neighbors"):
    key = adata.uns[neighbors_key]["connectivities_key"]
    matrix = adata.obsp[key].astype("float64").tolil(copy=True)
    coo = matrix.tocoo()
    row, column = next(
        (int(left), int(right)) for left, right in zip(coo.row, coo.col, strict=True) if int(left) < int(right)
    )
    matrix[row, column] = float(matrix[row, column]) + 4e-9
    matrix[column, row] = float(matrix[column, row]) - 4e-9
    adata.obsp[key] = matrix.tocsr()


def test_diffusion_map_real_backend_report_state_and_standalone_code(science):
    adata = _trajectory_input(science, neighbors_key="custom_graph")
    connectivity_before = adata.obsp["custom_graph_connectivities"].copy()
    distance_before = adata.obsp["custom_graph_distances"].copy()
    rng_before = science.np.random.get_state()
    print_before = science.np.get_printoptions()

    output, report, code = OpenBioSingleCellDiffusionMap.execute(
        adata,
        neighbors_key="custom_graph",
        n_comps=7,
        random_seed=17,
    ).result

    assert "X_diffmap" not in adata.obsm
    assert "diffmap_evals" not in adata.uns
    _assert_sparse_equal(science, adata.obsp["custom_graph_connectivities"], connectivity_before)
    _assert_sparse_equal(science, adata.obsp["custom_graph_distances"], distance_before)
    _assert_rng_state_equal(science, rng_before, science.np.random.get_state())
    assert science.np.get_printoptions() == print_before
    assert output.obsm["X_diffmap"].shape == (adata.n_obs, 7)
    assert output.uns["diffmap_evals"].shape == (7,)
    assert report.summary["key_results"]["stationary_component_index"] == 0
    assert len(report.summary["key_results"]["eigenvalues"]) == 7
    assert report.summary["software_versions"]["scanpy"] == "1.12.3"
    assert len(report.summary["references"]) >= 3
    json.dumps(report.summary, allow_nan=False)
    assert "openbio_singlecell" not in code
    compile(code, "<diffusion-map-code>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated, generated_summary = namespace["run_diffusion_map"](adata)
    science.np.testing.assert_allclose(generated.obsm["X_diffmap"], output.obsm["X_diffmap"], rtol=0, atol=0)
    science.np.testing.assert_allclose(generated.uns["diffmap_evals"], output.uns["diffmap_evals"], rtol=0, atol=0)
    assert (
        generated_summary["key_results"]["outputs"]["output_fingerprint_sha256"]
        == report.summary["key_results"]["outputs"]["output_fingerprint_sha256"]
    )


def test_diffusion_map_treats_stored_zero_diagonal_as_absent_in_runtime_and_code(science):
    baseline = _trajectory_input(science)
    with_zero_diagonal = baseline.copy()
    storage_before = {}
    for field in ("connectivities_key", "distances_key"):
        key = with_zero_diagonal.uns["neighbors"][field]
        coo = with_zero_diagonal.obsp[key].tocoo()
        matrix = science.sparse.csr_matrix(
            (
                science.np.concatenate([coo.data, science.np.asarray([0.0])]),
                (
                    science.np.concatenate([coo.row, science.np.asarray([0])]),
                    science.np.concatenate([coo.col, science.np.asarray([0])]),
                ),
            ),
            shape=coo.shape,
        )
        diagonal = matrix.tocoo()
        assert bool(((diagonal.row == 0) & (diagonal.col == 0) & (diagonal.data == 0.0)).any())
        with_zero_diagonal.obsp[key] = matrix
        storage_before[key] = (matrix.data.copy(), matrix.indices.copy(), matrix.indptr.copy())

    _, baseline_summary = analyze_diffusion_map(baseline, n_comps=5, random_seed=19)
    _, report, code = OpenBioSingleCellDiffusionMap.execute(
        with_zero_diagonal,
        n_comps=5,
        random_seed=19,
    ).result

    baseline_fingerprint = baseline_summary["key_results"]["graph"]["graph_fingerprint_sha256"]
    assert report.summary["key_results"]["graph"]["graph_fingerprint_sha256"] == baseline_fingerprint
    for key, (data, indices, indptr) in storage_before.items():
        science.np.testing.assert_array_equal(with_zero_diagonal.obsp[key].data, data)
        science.np.testing.assert_array_equal(with_zero_diagonal.obsp[key].indices, indices)
        science.np.testing.assert_array_equal(with_zero_diagonal.obsp[key].indptr, indptr)

    namespace: dict[str, object] = {}
    exec(code, namespace)
    _, generated_summary = namespace["run_diffusion_map"](with_zero_diagonal)
    assert generated_summary["key_results"]["graph"]["graph_fingerprint_sha256"] == baseline_fingerprint


def test_diffusion_map_preserves_off_diagonal_zero_distances_and_rejects_self_loops(science):
    baseline = _trajectory_input(science)
    with_zero_distance = baseline.copy()
    distance_key = with_zero_distance.uns["neighbors"]["distances_key"]
    distance = with_zero_distance.obsp[distance_key].tocsr()
    row, column = next(
        (row, column)
        for row in range(distance.shape[0])
        for column in range(distance.shape[1])
        if row != column and distance[row, column] == 0
    )
    coo = distance.tocoo()
    with_zero_distance.obsp[distance_key] = science.sparse.csr_matrix(
        (
            science.np.concatenate([coo.data, science.np.asarray([0.0])]),
            (
                science.np.concatenate([coo.row, science.np.asarray([row])]),
                science.np.concatenate([coo.col, science.np.asarray([column])]),
            ),
        ),
        shape=coo.shape,
    )

    _, baseline_summary = analyze_diffusion_map(baseline, n_comps=5, random_seed=29)
    _, zero_summary = analyze_diffusion_map(with_zero_distance, n_comps=5, random_seed=29)
    baseline_graph = baseline_summary["key_results"]["graph"]
    zero_graph = zero_summary["key_results"]["graph"]
    assert zero_graph["stored_directed_distance_edges"] == baseline_graph["stored_directed_distance_edges"] + 1
    assert zero_graph["graph_fingerprint_sha256"] != baseline_graph["graph_fingerprint_sha256"]

    self_loop = baseline.copy()
    connectivity_key = self_loop.uns["neighbors"]["connectivities_key"]
    connectivity = self_loop.obsp[connectivity_key].tolil(copy=True)
    connectivity[0, 0] = 0.25
    self_loop.obsp[connectivity_key] = connectivity.tocsr()
    with pytest.raises(ValueError, match="zero diagonal"):
        analyze_diffusion_map(self_loop, n_comps=5)


@pytest.mark.parametrize("failure", ["missing_pointer", "dense", "asymmetric", "nonfinite", "duplicate_axis"])
def test_trajectory_named_graph_preflight_is_fail_closed(science, failure):
    adata = _trajectory_input(science)
    if failure == "missing_pointer":
        del adata.uns["neighbors"]["distances_key"]
        expected = "distances_key"
    elif failure == "dense":
        key = adata.uns["neighbors"]["connectivities_key"]
        adata.obsp[key] = adata.obsp[key].toarray()
        expected = "sparse"
    elif failure == "asymmetric":
        key = adata.uns["neighbors"]["connectivities_key"]
        matrix = adata.obsp[key].tolil(copy=True)
        matrix[0, 1] = float(matrix[0, 1]) + 0.5
        adata.obsp[key] = matrix.tocsr()
        expected = "symmetric"
    elif failure == "nonfinite":
        key = adata.uns["neighbors"]["distances_key"]
        matrix = adata.obsp[key].tocsr(copy=True)
        matrix.data[0] = science.np.nan
        adata.obsp[key] = matrix
        expected = "non-finite"
    else:
        adata.obs_names = ["duplicate", "duplicate", *list(adata.obs_names[2:])]
        expected = "unique"
    with pytest.raises((TypeError, ValueError), match=expected):
        analyze_diffusion_map(adata, n_comps=5)


def test_diffusion_map_fingerprints_and_consumes_the_same_canonical_graph(science):
    base = _trajectory_input(science)
    key = base.uns["neighbors"]["connectivities_key"]
    coo = base.obsp[key].tocoo()
    edge = next((int(left), int(right)) for left, right in zip(coo.row, coo.col, strict=True) if int(left) < int(right))
    left = base.copy()
    right = base.copy()
    for target, sign in ((left, 1.0), (right, -1.0)):
        matrix = target.obsp[key].astype("float64").tolil(copy=True)
        row, column = edge
        matrix[row, column] = float(matrix[row, column]) + sign * 4e-9
        matrix[column, row] = float(matrix[column, row]) - sign * 4e-9
        target.obsp[key] = matrix.tocsr()

    left_output, left_summary = analyze_diffusion_map(left, n_comps=5, random_seed=23)
    right_output, right_summary = analyze_diffusion_map(right, n_comps=5, random_seed=23)

    _assert_sparse_equal(science, left_output.obsp[key], right_output.obsp[key])
    assert (left_output.obsp[key] - left_output.obsp[key].T).nnz == 0
    assert (
        left_summary["key_results"]["graph"]["graph_fingerprint_sha256"]
        == right_summary["key_results"]["graph"]["graph_fingerprint_sha256"]
    )
    assert (
        left_summary["key_results"]["outputs"]["output_fingerprint_sha256"]
        == right_summary["key_results"]["outputs"]["output_fingerprint_sha256"]
    )


def test_diffusion_map_bounds_collisions_overwrite_and_backend_postcondition(science, monkeypatch):
    adata = _trajectory_input(science)
    with pytest.raises(ValueError, match="at least 3|between 3"):
        analyze_diffusion_map(adata, n_comps=2)
    with pytest.raises(ValueError, match="n_obs - 1"):
        analyze_diffusion_map(adata, n_comps=adata.n_obs)
    first, _ = analyze_diffusion_map(adata, n_comps=5)
    with pytest.raises(ValueError, match="already exists"):
        analyze_diffusion_map(first, n_comps=5)
    replaced, _ = analyze_diffusion_map(first, n_comps=6, overwrite_existing=True)
    assert replaced.obsm["X_diffmap"].shape == (adata.n_obs, 6)

    real_backend = science.sc.tl.diffmap

    def corrupt_backend(adata, n_comps=15, *, neighbors_key=None, random_state=0, copy=False):
        result = real_backend(
            adata,
            n_comps=n_comps,
            neighbors_key=neighbors_key,
            random_state=random_state,
            copy=copy,
        )
        adata.uns["diffmap_evals"][[0, 1]] = adata.uns["diffmap_evals"][[1, 0]]
        return result

    monkeypatch.setattr(science.sc.tl, "diffmap", corrupt_backend)
    with pytest.raises(RuntimeError, match="non-increasing"):
        analyze_diffusion_map(adata, n_comps=5)


def test_diffusion_map_rejects_structurally_valid_fake_eigenpairs_and_stale_dpt(science, monkeypatch):
    adata = _trajectory_input(science)

    def fake_backend(adata, n_comps=15, *, neighbors_key=None, random_state=0, copy=False):
        del neighbors_key, random_state, copy
        rng = science.np.random.default_rng(101)
        basis, _ = science.np.linalg.qr(rng.normal(size=(adata.n_obs, n_comps)))
        adata.obsm["X_diffmap"] = basis.astype("float32")
        adata.uns["diffmap_evals"] = science.np.linspace(0.9, 0.1, n_comps, dtype="float32")

    monkeypatch.setattr(science.sc.tl, "diffmap", fake_backend)
    with pytest.raises(RuntimeError, match="eigenpair|eigenbasis"):
        analyze_diffusion_map(adata, n_comps=5)
    assert "X_diffmap" not in adata.obsm
    assert "diffmap_evals" not in adata.uns

    monkeypatch.undo()
    diffmapped = _diffmapped(science)
    dpt_output, _ = analyze_dpt(diffmapped, root_cell_id="cell_000", n_dcs=5)
    with pytest.raises(ValueError, match="downstream DPT"):
        analyze_diffusion_map(dpt_output, n_comps=6, overwrite_existing=True)


def test_diffusion_map_component_summary_budget_is_explicit(science, monkeypatch):
    monkeypatch.setattr(trajectory_analysis, "_TA_MAX_DIFFMAP_COMPONENT_SUMMARIES", 2)
    _output, summary = analyze_diffusion_map(_trajectory_input(science), n_comps=6)
    results = summary["key_results"]
    assert results["informative_component_summary_count"] == 5
    assert results["informative_component_summaries_truncated"] is True
    assert results["informative_component_summaries_limit"] == 2
    assert len(results["informative_component_summaries"]) == 2
    assert any("truncated" in warning for warning in summary["warnings"])


def test_all_trajectory_backends_reject_zero_average_graph_tampering(science, monkeypatch):
    diffmap_input = _trajectory_input(science)
    real_diffmap = science.sc.tl.diffmap

    def tampering_diffmap(adata, n_comps=15, *, neighbors_key=None, random_state=0, copy=False):
        result = real_diffmap(
            adata,
            n_comps=n_comps,
            neighbors_key=neighbors_key,
            random_state=random_state,
            copy=copy,
        )
        _tamper_named_connectivities_with_zero_average(adata, neighbors_key=neighbors_key or "neighbors")
        return result

    with monkeypatch.context() as scoped:
        scoped.setattr(science.sc.tl, "diffmap", tampering_diffmap)
        with pytest.raises(RuntimeError, match="exact canonical named graph"):
            analyze_diffusion_map(diffmap_input, n_comps=5)

    paga_input = _trajectory_input(science)
    real_paga = science.sc.tl.paga

    def tampering_paga(
        adata,
        groups=None,
        *,
        use_rna_velocity=False,
        model="v1.2",
        neighbors_key=None,
        copy=False,
    ):
        result = real_paga(
            adata,
            groups=groups,
            use_rna_velocity=use_rna_velocity,
            model=model,
            neighbors_key=neighbors_key,
            copy=copy,
        )
        _tamper_named_connectivities_with_zero_average(adata, neighbors_key=neighbors_key or "neighbors")
        return result

    with monkeypatch.context() as scoped:
        scoped.setattr(science.sc.tl, "paga", tampering_paga)
        with pytest.raises(RuntimeError, match="exact canonical named graph"):
            analyze_paga(paga_input, groupby="state")

    dpt_input = _diffmapped(science)
    real_dpt = science.sc.tl.dpt

    def tampering_dpt(
        adata,
        n_dcs=10,
        *,
        n_branchings=0,
        min_group_size=0.01,
        allow_kendall_tau_shift=True,
        neighbors_key=None,
        copy=False,
    ):
        result = real_dpt(
            adata,
            n_dcs=n_dcs,
            n_branchings=n_branchings,
            min_group_size=min_group_size,
            allow_kendall_tau_shift=allow_kendall_tau_shift,
            neighbors_key=neighbors_key,
            copy=copy,
        )
        _tamper_named_connectivities_with_zero_average(adata, neighbors_key=neighbors_key or "neighbors")
        return result

    with monkeypatch.context() as scoped:
        scoped.setattr(science.sc.tl, "dpt", tampering_dpt)
        with pytest.raises(RuntimeError, match="exact canonical named graph"):
            analyze_dpt(dpt_input, root_cell_id="cell_000", n_dcs=5)


def test_paga_real_backend_unused_categories_report_and_standalone_code(science):
    adata = _trajectory_input(science, neighbors_key="custom_graph")
    input_categories = list(adata.obs["state"].cat.categories)
    output, report, code = OpenBioSingleCellPAGA.execute(
        adata,
        groupby="state",
        neighbors_key="custom_graph",
    ).result

    assert list(adata.obs["state"].cat.categories) == input_categories
    assert list(output.obs["state"].cat.categories) == ["early", "middle", "late"]
    assert output.uns["paga"]["groups"] == "state"
    assert output.uns["paga"]["connectivities"].shape == (3, 3)
    assert output.uns["state_sizes"].tolist() == [12, 12, 12]
    tree = output.uns["paga"]["connectivities_tree"]
    reciprocal = tree.multiply(tree.T)
    reciprocal.eliminate_zeros()
    assert reciprocal.nnz == 0
    abstract = report.summary["key_results"]["abstract_graph"]
    assert abstract["edge_table_truncated"] is False
    assert len(abstract["all_edges"]) == abstract["nonzero_undirected_edges"]
    assert report.summary["parameters"]["model"] == "v1.2"
    assert report.summary["parameters"]["use_rna_velocity"] is False
    json.dumps(report.summary, allow_nan=False)
    assert "openbio_singlecell" not in code
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated, generated_summary = namespace["run_paga"](adata)
    _assert_sparse_equal(
        science,
        generated.uns["paga"]["connectivities"],
        output.uns["paga"]["connectivities"],
    )
    _assert_sparse_equal(
        science,
        generated.uns["paga"]["connectivities_tree"],
        output.uns["paga"]["connectivities_tree"],
    )
    assert generated.uns["state_sizes"].tolist() == output.uns["state_sizes"].tolist()
    assert generated_summary["key_results"]["abstract_graph"] == abstract


def test_paga_partition_and_collision_guards(science):
    adata = _trajectory_input(science)
    adata.obs["plain"] = adata.obs["state"].astype(str)
    with pytest.raises(TypeError, match="categorical"):
        analyze_paga(adata, groupby="plain")

    one_group = adata.copy()
    one_group.obs["one"] = science.pd.Categorical(["same"] * one_group.n_obs)
    with pytest.raises(ValueError, match="at least two"):
        analyze_paga(one_group, groupby="one")

    missing = adata.copy()
    missing.obs["state"] = missing.obs["state"].astype(object)
    missing.obs.iloc[0, missing.obs.columns.get_loc("state")] = None
    missing.obs["state"] = science.pd.Categorical(missing.obs["state"])
    with pytest.raises(ValueError, match="missing"):
        analyze_paga(missing, groupby="state")

    first, _ = analyze_paga(adata, groupby="state")
    with pytest.raises(ValueError, match="already exists"):
        analyze_paga(first, groupby="state")
    replaced, _ = analyze_paga(first, groupby="state", overwrite_existing=True)
    assert replaced.uns["paga"]["groups"] == "state"


def test_paga_rejects_backend_tree_or_size_tampering(science, monkeypatch):
    adata = _trajectory_input(science)
    real_backend = science.sc.tl.paga

    def corrupt_backend(
        adata,
        groups=None,
        *,
        use_rna_velocity=False,
        model="v1.2",
        neighbors_key=None,
        copy=False,
    ):
        result = real_backend(
            adata,
            groups=groups,
            use_rna_velocity=use_rna_velocity,
            model=model,
            neighbors_key=neighbors_key,
            copy=copy,
        )
        adata.uns[f"{groups}_sizes"][0] += 1
        return result

    monkeypatch.setattr(science.sc.tl, "paga", corrupt_backend)
    with pytest.raises(RuntimeError, match="group-size"):
        analyze_paga(adata, groupby="state")


def test_paga_recomputes_model_v1_2_and_isolates_backend_side_effects(science, monkeypatch):
    adata = _trajectory_input(science)
    adata.obs["unrelated_text"] = [f"text-{index % 2}" for index in range(adata.n_obs)]
    adata.var["unrelated_var_text"] = [f"gene-text-{index % 2}" for index in range(adata.n_vars)]
    obs_dtype = str(adata.obs["unrelated_text"].dtype)
    var_dtype = str(adata.var["unrelated_var_text"].dtype)
    real_backend = science.sc.tl.paga
    clean, _ = analyze_paga(adata, groupby="state")
    assert str(clean.obs["unrelated_text"].dtype) == obs_dtype
    assert str(clean.var["unrelated_var_text"].dtype) == var_dtype
    assert str(adata.obs["unrelated_text"].dtype) == obs_dtype
    assert str(adata.var["unrelated_var_text"].dtype) == var_dtype

    def fake_backend(
        adata,
        groups=None,
        *,
        use_rna_velocity=False,
        model="v1.2",
        neighbors_key=None,
        copy=False,
    ):
        del use_rna_velocity, model, neighbors_key, copy
        categories = adata.obs[groups].cat.categories
        shape = (len(categories), len(categories))
        adata.uns["paga"] = {
            "groups": groups,
            "connectivities": science.sparse.csr_matrix(shape, dtype=float),
            "connectivities_tree": science.sparse.csr_matrix(shape, dtype=float),
        }
        adata.uns[f"{groups}_sizes"] = science.np.asarray(adata.obs[groups].value_counts(sort=False), dtype=int)

    assert real_backend is not fake_backend
    monkeypatch.setattr(science.sc.tl, "paga", fake_backend)
    with pytest.raises(RuntimeError, match="independent model-v1.2"):
        analyze_paga(adata, groupby="state")
    assert "paga" not in adata.uns
    assert "state_sizes" not in adata.uns


def test_paga_overwrite_clears_the_previous_group_size_sidecar(science):
    adata = _trajectory_input(science)
    first, _ = analyze_paga(adata, groupby="state")
    first.obs["phase"] = science.pd.Categorical(
        ["one"] * 18 + ["two"] * 18,
        categories=["one", "two"],
    )
    second, _ = analyze_paga(first, groupby="phase", overwrite_existing=True)
    assert "state_sizes" not in second.uns
    assert second.uns["phase_sizes"].tolist() == [18, 18]

    malformed = _trajectory_input(science)
    malformed.uns["paga"] = {"connectivities": science.sparse.eye(2, format="csr")}
    snapshot = malformed.copy()
    with pytest.raises(ValueError, match="previous canonical group-size sidecar"):
        analyze_paga(malformed, groupby="state", overwrite_existing=True)
    assert malformed.uns["paga"].keys() == snapshot.uns["paga"].keys()


def test_paga_summary_edge_budget_is_explicit_and_strict_json(science, monkeypatch):
    monkeypatch.setattr(trajectory_analysis, "_TA_MAX_PAGA_SUMMARY_EDGES", 1)
    _output, summary = analyze_paga(_trajectory_input(science), groupby="state")
    abstract = summary["key_results"]["abstract_graph"]
    assert abstract["nonzero_undirected_edges"] > 1
    assert abstract["edge_table_truncated"] is True
    assert abstract["edge_table_limit"] == 1
    assert len(abstract["all_edges"]) == 1
    assert any("truncated" in warning for warning in summary["warnings"])
    json.dumps(summary, allow_nan=False)


def _diffmapped(science):
    adata = _trajectory_input(science)
    output, _ = analyze_diffusion_map(adata, n_comps=7, random_seed=19)
    return output


def test_dpt_requires_graph_bound_diffusion_provenance_and_detects_staleness(science):
    bare = _trajectory_input(science)
    science.sc.tl.diffmap(bare, n_comps=7, neighbors_key="neighbors", random_state=19)
    with pytest.raises(ValueError, match="OpenBio Diffusion Map provenance"):
        analyze_dpt(bare, root_cell_id="cell_000", n_dcs=5)

    stale = _diffmapped(science)
    connectivity_key = stale.uns["neighbors"]["connectivities_key"]
    connectivity = stale.obsp[connectivity_key].tolil(copy=True)
    current = float(connectivity[0, 1])
    connectivity[0, 1] = current + 0.01
    connectivity[1, 0] = current + 0.01
    stale.obsp[connectivity_key] = connectivity.tocsr()
    with pytest.raises(ValueError, match="provenance mismatch"):
        analyze_dpt(stale, root_cell_id="cell_000", n_dcs=5)


def test_dpt_exact_cell_root_report_and_standalone_code(science):
    adata = _diffmapped(science)
    output, report, code = OpenBioSingleCellDPT.execute(
        adata,
        neighbors_key="neighbors",
        root_mode={"root_mode": "cell_id", "root_cell_id": "cell_004"},
        n_dcs=5,
    ).result

    assert "dpt_pseudotime" not in adata.obs
    assert output.uns["iroot"] == 4
    assert output.obs["dpt_pseudotime"].iloc[4] == pytest.approx(0.0, abs=1e-7)
    assert output.obs["dpt_pseudotime"].between(0, 1).all()
    assert report.summary["key_results"]["root"]["root_cell_id"] == "cell_004"
    assert report.summary["parameters"]["scanpy_fixed_arguments"]["n_branchings"] == 0
    assert len(report.summary["references"]) >= 2
    json.dumps(report.summary, allow_nan=False)
    assert "openbio_singlecell" not in code
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated, generated_summary = namespace["run_dpt"](adata)
    science.np.testing.assert_allclose(generated.obs["dpt_pseudotime"], output.obs["dpt_pseudotime"], rtol=0, atol=0)
    assert generated_summary["key_results"]["root"] == report.summary["key_results"]["root"]


def test_dpt_group_medoid_uses_exact_typed_label_and_disclosed_rule(science):
    adata = _diffmapped(science)
    adata.obs["numeric_stage"] = science.pd.Categorical(
        [1] * 12 + [2] * 12 + [3] * 12,
        categories=[1, 2, 3],
    )
    output, summary = analyze_dpt(
        adata,
        root_mode="group_medoid",
        root_column="numeric_stage",
        root_value=1,
        n_dcs=5,
    )
    root = summary["key_results"]["root"]
    assert root["population"]["value"] == {"type": "integer", "value": 1, "display": "1"}
    assert root["population"]["cell_count"] == 12
    assert output.obs["numeric_stage"].iloc[output.uns["iroot"]] == 1
    with pytest.raises(ValueError, match="not found exactly"):
        analyze_dpt(
            adata,
            root_mode="group_medoid",
            root_column="numeric_stage",
            root_value="1",
            n_dcs=5,
        )

    node_output, node_report, _code = OpenBioSingleCellDPT.execute(
        adata,
        root_mode={
            "root_mode": "group_medoid",
            "root_column": "numeric_stage",
            "root_value_type": "integer",
            "root_value": "1",
        },
        n_dcs=5,
    ).result
    assert node_output.obs["numeric_stage"].iloc[node_output.uns["iroot"]] == 1
    assert node_report.summary["parameters"]["root_value"]["type"] == "integer"
    adata.obs["numeric_float_stage"] = science.pd.Categorical([1.5] * 12 + [2.5] * 24)
    float_output, float_report, _float_code = OpenBioSingleCellDPT.execute(
        adata,
        root_mode={
            "root_mode": "group_medoid",
            "root_column": "numeric_float_stage",
            "root_value_type": "number",
            "root_value": "1.5",
        },
        n_dcs=5,
    ).result
    assert float_output.obs["numeric_float_stage"].iloc[float_output.uns["iroot"]] == 1.5
    assert float_report.summary["parameters"]["root_value"]["type"] == "number"
    with pytest.raises(ValueError, match="canonical base-10"):
        OpenBioSingleCellDPT.execute(
            adata,
            root_mode={
                "root_mode": "group_medoid",
                "root_column": "numeric_stage",
                "root_value_type": "integer",
                "root_value": "01",
            },
            n_dcs=5,
        )
    with pytest.raises(ValueError, match="finite numeric"):
        OpenBioSingleCellDPT.execute(
            adata,
            root_mode={
                "root_mode": "group_medoid",
                "root_column": "numeric_stage",
                "root_value_type": "number",
                "root_value": "nan",
            },
            n_dcs=5,
        )


def test_dpt_group_medoid_is_cell_id_stable_under_input_reordering(science):
    adata = science.ad.AnnData(science.np.zeros((4, 1)))
    adata.obs_names = ["cell_d", "cell_b", "cell_a", "cell_c"]
    adata.obs["root_group"] = science.pd.Categorical(["root"] * 4)
    coordinates = science.np.asarray(
        [
            [0.0, 1e16],
            [0.0, -1e16],
            [0.0, 1.0],
            [0.0, -1.0],
        ]
    )
    first_index, first = _ta_select_root(
        adata,
        root_mode="group_medoid",
        root_cell_id="",
        root_column="root_group",
        root_value="root",
        coordinates=coordinates,
        n_dcs=2,
        numpy=science.np,
    )
    order = [2, 0, 3, 1]
    reordered = adata[order].copy()
    second_index, second = _ta_select_root(
        reordered,
        root_mode="group_medoid",
        root_cell_id="",
        root_column="root_group",
        root_value="root",
        coordinates=coordinates[order],
        n_dcs=2,
        numpy=science.np,
    )
    assert adata.obs_names[first_index] == reordered.obs_names[second_index]
    assert first["root_cell_id"] == second["root_cell_id"]
    assert first["population"]["tie_tolerance"] > 0


def test_dpt_bounds_disconnected_graph_collisions_and_backend_postcondition(science, monkeypatch):
    adata = _diffmapped(science)
    with pytest.raises(ValueError, match="exceeds"):
        analyze_dpt(adata, root_cell_id="cell_000", n_dcs=8)
    first, _ = analyze_dpt(adata, root_cell_id="cell_000", n_dcs=5)
    with pytest.raises(ValueError, match="already exists"):
        analyze_dpt(first, root_cell_id="cell_000", n_dcs=5)
    replaced, _ = analyze_dpt(first, root_cell_id="cell_001", n_dcs=5, overwrite_existing=True)
    assert replaced.uns["iroot"] == 1

    real_backend = science.sc.tl.dpt

    def corrupt_backend(
        adata,
        n_dcs=10,
        *,
        n_branchings=0,
        min_group_size=0.01,
        allow_kendall_tau_shift=True,
        neighbors_key=None,
        copy=False,
    ):
        result = real_backend(
            adata,
            n_dcs=n_dcs,
            n_branchings=n_branchings,
            min_group_size=min_group_size,
            allow_kendall_tau_shift=allow_kendall_tau_shift,
            neighbors_key=neighbors_key,
            copy=copy,
        )
        adata.obs.iloc[int(adata.uns["iroot"]), adata.obs.columns.get_loc("dpt_pseudotime")] = 0.5
        return result

    monkeypatch.setattr(science.sc.tl, "dpt", corrupt_backend)
    with pytest.raises(RuntimeError, match="zero pseudotime"):
        analyze_dpt(adata, root_cell_id="cell_000", n_dcs=5)


def test_dpt_rejects_structurally_valid_fake_pseudotime(science, monkeypatch):
    adata = _diffmapped(science)

    def fake_backend(
        adata,
        n_dcs=10,
        *,
        n_branchings=0,
        min_group_size=0.01,
        allow_kendall_tau_shift=True,
        neighbors_key=None,
        copy=False,
    ):
        del n_dcs, n_branchings, min_group_size, allow_kendall_tau_shift, neighbors_key, copy
        values = science.np.linspace(0.0, 1.0, adata.n_obs)
        values[int(adata.uns["iroot"])] = 0.0
        adata.obs["dpt_pseudotime"] = values

    monkeypatch.setattr(science.sc.tl, "dpt", fake_backend)
    with pytest.raises(RuntimeError, match="independently recomputed diffusion distance"):
        analyze_dpt(adata, root_cell_id="cell_000", n_dcs=5)
    assert "dpt_pseudotime" not in adata.obs
    assert "openbio_dpt" not in adata.uns


def test_paga_and_dpt_preserve_numpy_process_state(science):
    adata = _trajectory_input(science)
    rng_before = science.np.random.get_state()
    print_before = science.np.get_printoptions()
    analyze_paga(adata, groupby="state")
    _assert_rng_state_equal(science, rng_before, science.np.random.get_state())
    assert science.np.get_printoptions() == print_before

    mapped = _diffmapped(science)
    rng_before = science.np.random.get_state()
    print_before = science.np.get_printoptions()
    analyze_dpt(mapped, root_cell_id="cell_000", n_dcs=5)
    _assert_rng_state_equal(science, rng_before, science.np.random.get_state())
    assert science.np.get_printoptions() == print_before


def test_dpt_rejects_disconnected_named_graph_even_with_valid_diffmap(science):
    n_obs = 12
    adata = science.ad.AnnData(science.np.random.default_rng(3).normal(size=(n_obs, 4)))
    adata.obs_names = [f"d{index}" for index in range(n_obs)]
    rows = []
    columns = []
    for start in (0, 6):
        for index in range(start, start + 5):
            rows.extend([index, index + 1])
            columns.extend([index + 1, index])
    weights = science.np.ones(len(rows), dtype=float)
    matrix = science.sparse.csr_matrix((weights, (rows, columns)), shape=(n_obs, n_obs))
    adata.obsp["block_connectivities"] = matrix.copy()
    adata.obsp["block_distances"] = matrix.copy()
    adata.uns["block"] = {
        "connectivities_key": "block_connectivities",
        "distances_key": "block_distances",
        "params": {"n_neighbors": 2, "method": "manual_test"},
    }
    diffmapped, summary = analyze_diffusion_map(adata, neighbors_key="block", n_comps=5)
    assert summary["key_results"]["graph"]["connected_components"] == 2
    with pytest.raises(ValueError, match="one connected graph component"):
        analyze_dpt(diffmapped, neighbors_key="block", root_cell_id="d0", n_dcs=4)


def test_trajectory_schemas_have_atomic_three_output_contracts():
    schemas = [
        OpenBioSingleCellDiffusionMap.GET_SCHEMA(),
        OpenBioSingleCellPAGA.GET_SCHEMA(),
        OpenBioSingleCellDPT.GET_SCHEMA(),
    ]
    for schema in schemas:
        assert [output.display_name for output in schema.outputs] == ["adata", "summary", "code"]
    assert {value.id for value in schemas[0].inputs} >= {
        "adata",
        "neighbors_key",
        "n_comps",
        "overwrite_existing",
        "random_seed",
    }
    assert {value.id for value in schemas[1].inputs} >= {
        "adata",
        "groupby",
        "neighbors_key",
        "overwrite_existing",
    }
    assert {value.id for value in schemas[2].inputs} >= {
        "adata",
        "neighbors_key",
        "root_mode",
        "n_dcs",
        "overwrite_existing",
    }
    root_mode = next(value for value in schemas[2].inputs if value.id == "root_mode")
    group_medoid = next(option for option in root_mode.options if option.key == "group_medoid")
    assert [value.id for value in group_medoid.inputs] == ["root_column", "root_value_type", "root_value"]


def test_generated_code_is_a_node_specific_dependency_closure(science):
    adata = _trajectory_input(science)
    _diffmap_output, _diffmap_report, diffmap_code = OpenBioSingleCellDiffusionMap.execute(adata, n_comps=5).result
    _paga_output, _paga_report, paga_code = OpenBioSingleCellPAGA.execute(adata, groupby="state").result
    mapped = _diffmapped(science)
    _dpt_output, _dpt_report, dpt_code = OpenBioSingleCellDPT.execute(
        mapped,
        root_mode={"root_mode": "cell_id", "root_cell_id": "cell_000"},
        n_dcs=5,
    ).result
    assert "def _ta_run_paga" not in diffmap_code and "def _ta_run_dpt" not in diffmap_code
    assert "def _ta_run_diffusion_map" not in paga_code and "def _ta_run_dpt" not in paga_code
    assert "def _ta_run_diffusion_map" not in dpt_code and "def _ta_run_paga" not in dpt_code
    assert max(len(diffmap_code), len(paga_code), len(dpt_code)) < 45_000
