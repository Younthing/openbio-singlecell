from __future__ import annotations

import copy
import inspect
import itertools
import json

import igraph
import pytest
from sklearn.metrics import adjusted_rand_score

from openbio_singlecell.graph_analysis import resolve_named_graph
from openbio_singlecell.nodes_integration import (
    OpenBioSingleCellLeidenResolutionSweep,
    _leiden_resolution_sweep_code,
    _parse_resolutions,
)

METRIC_COLUMNS = [
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


def _adata_with_named_graph(science):
    n_obs = 12
    values = science.np.arange(n_obs * 5, dtype=float).reshape(n_obs, 5)
    adata = science.ad.AnnData(
        X=values,
        obs=science.pd.DataFrame(
            {"sample": science.pd.Categorical(["s1"] * 6 + ["s2"] * 6)},
            index=[f"cell_{index}" for index in range(n_obs)],
        ),
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(5)]),
    )
    graph = science.sparse.lil_matrix((n_obs, n_obs), dtype=float)
    for start in (0, 4, 8):
        for left in range(start, start + 4):
            for right in range(left + 1, start + 4):
                graph[left, right] = 1.0
                graph[right, left] = 1.0
    for left, right in ((3, 4), (7, 8)):
        graph[left, right] = 0.05
        graph[right, left] = 0.05
    adata.obsp["custom_connectivities"] = graph.tocsr()
    adata.uns["custom_neighbors"] = {
        "connectivities_key": "custom_connectivities",
        "params": {
            "n_neighbors": 4,
            "metric": "cosine",
            "method": "umap",
            "use_rep": "X_test",
            "random_state": 17,
        },
    }
    return adata


def _execute(adata, **overrides):
    arguments = {
        "resolutions": "0.25,0.5,1.0,2.0",
        "key_prefix": "leiden",
        "neighbors_key": "custom_neighbors",
        "n_iterations": 2,
        "stability_repeats": 5,
        "random_seed": 0,
        "overwrite_existing": False,
    }
    arguments.update(overrides)
    return OpenBioSingleCellLeidenResolutionSweep.execute(adata, **arguments).result


def _metric_table(table_result):
    assert table_result.kind == "table"
    return table_result.table


def _membership_keys(table):
    return [str(value) for value in table["key"].tolist()]


def _pairwise_ari(science, memberships):
    values = [adjusted_rand_score(left, right) for left, right in itertools.combinations(memberships, 2)]
    return {
        "mean": float(science.np.mean(values)),
        "min": float(science.np.min(values)),
        "max": float(science.np.max(values)),
    }


def _independent_modularity(science, connectivities, labels, *, resolution):
    graph = igraph.Graph.Weighted_Adjacency(
        science.np.asarray(connectivities.toarray(), dtype=float).tolist(),
        mode="undirected",
        attr="weight",
        loops=False,
    )
    membership = tuple(str(value) for value in labels)
    label_ids = {label: index for index, label in enumerate(dict.fromkeys(membership))}
    return float(
        graph.modularity(
            [label_ids[label] for label in membership],
            weights="weight",
            resolution=resolution,
            directed=False,
        )
    )


def test_schema_exposes_structured_sweep_contract():
    schema = OpenBioSingleCellLeidenResolutionSweep.GET_SCHEMA()
    assert [item.id for item in schema.inputs] == [
        "adata",
        "resolutions",
        "key_prefix",
        "neighbors_key",
        "n_iterations",
        "stability_repeats",
        "random_seed",
        "overwrite_existing",
    ]
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("adata", "OPENBIO_ANNDATA"),
        ("resolution_metrics", "OPENBIO_SINGLE_CELL_TABLE"),
        ("summary", "OPENBIO_SINGLE_CELL_SUMMARY"),
        ("code", "STRING"),
    ]
    execute_parameters = list(inspect.signature(OpenBioSingleCellLeidenResolutionSweep.execute).parameters)
    assert execute_parameters[-2:] == ["random_seed", "overwrite_existing"]


def test_valid_custom_graph_runs_exactly_resolution_count_times_total_starts(science, monkeypatch):
    adata = _adata_with_named_graph(science)
    original_leiden = science.sc.tl.leiden
    calls = []

    def counted_leiden(*args, **kwargs):
        calls.append(dict(kwargs))
        return original_leiden(*args, **kwargs)

    monkeypatch.setattr(science.sc.tl, "leiden", counted_leiden)
    output, table_result, _, _ = _execute(
        adata,
        resolutions="0.4,0.8,1.2",
        stability_repeats=3,
        random_seed=11,
    )
    table = _metric_table(table_result)

    assert output is not adata
    assert table["resolution"].tolist() == [0.4, 0.8, 1.2]
    assert len(calls) == 3 * 3
    assert sorted(call["resolution"] for call in calls) == [0.4] * 3 + [0.8] * 3 + [1.2] * 3
    assert all(call["flavor"] == "igraph" for call in calls)
    assert all(call["directed"] is False for call in calls)
    assert all(call["use_weights"] is True for call in calls)
    assert all(call["n_iterations"] == 2 for call in calls)
    expected_adjacency = adata.obsp["custom_connectivities"]
    for call in calls:
        assert "neighbors_key" not in call
        adjacency = call["adjacency"]
        assert science.sparse.isspmatrix_csr(adjacency)
        assert adjacency.has_canonical_format
        difference = adjacency - expected_adjacency
        difference.eliminate_zeros()
        assert difference.nnz == 0


@pytest.mark.parametrize(
    ("resolutions", "message"),
    [
        ("", "resolution"),
        ("not-a-number", "resolution"),
        ("0.5,0.5", "duplicate"),
        ("1,1.0", "duplicate"),
        ("nan", "finite"),
        ("inf", "finite"),
        ("-0.5", "non-negative"),
    ],
)
def test_resolution_parser_is_strict_and_bounded(science, monkeypatch, resolutions, message):
    backend_calls = []

    def unexpected_backend_call(*args, **kwargs):
        backend_calls.append((args, kwargs))
        raise AssertionError("Leiden backend ran before resolution validation completed")

    monkeypatch.setattr(science.sc.tl, "leiden", unexpected_backend_call)
    with pytest.raises((TypeError, ValueError), match=message):
        _execute(_adata_with_named_graph(science), resolutions=resolutions)
    assert backend_calls == []


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"key_prefix": ""}, "prefix"),
        ({"random_seed": -1}, "random_seed"),
        ({"random_seed": True}, "random_seed"),
        ({"random_seed": 1.5}, "random_seed"),
        ({"random_seed": 2**31 - 2, "stability_repeats": 3}, "seed"),
    ],
)
def test_sweep_settings_and_total_run_count_are_bounded(science, monkeypatch, overrides, message):
    backend_calls = []

    def unexpected_backend_call(*args, **kwargs):
        backend_calls.append((args, kwargs))
        raise AssertionError("Leiden backend ran before sweep settings validation completed")

    monkeypatch.setattr(science.sc.tl, "leiden", unexpected_backend_call)
    with pytest.raises((TypeError, ValueError), match=message):
        _execute(_adata_with_named_graph(science), **overrides)
    assert backend_calls == []


def test_open_expert_boundaries_and_single_repeat_nulls(science):
    adata = _adata_with_named_graph(science)
    output, table_result, report, code = _execute(
        adata,
        resolutions="0",
        key_prefix="expert result",
        n_iterations=0,
        stability_repeats=1,
    )
    row = _metric_table(table_result).iloc[0]

    assert row["key"] == "expert result_0"
    assert row[["stability_mean_ari", "stability_min_ari", "stability_max_ari"]].isna().all()
    summary_row = report.summary["key_results"]["resolution_metrics"][0]
    assert summary_row["stability_mean_ari"] is None
    assert summary_row["stability_min_ari"] is None
    assert summary_row["stability_max_ari"] is None
    assert output.uns["expert result_0"]["openbio_diagnostics"]["stability"]["pairwise_ari"] == []
    assert any("not assessed" in warning for warning in report.summary["warnings"])
    json.dumps(report.summary, allow_nan=False)

    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning):
        generated_output, generated_table = namespace["leiden_resolution_sweep"](adata)
    generated_row = generated_table.iloc[0]
    assert generated_row[["stability_mean_ari", "stability_min_ari", "stability_max_ari"]].isna().all()
    science.pd.testing.assert_series_equal(generated_output.obs["expert result_0"], output.obs["expert result_0"])


def test_resolution_parser_has_no_arbitrary_count_ceiling():
    values = _parse_resolutions(",".join(str(index / 10) for index in range(257)))
    assert len(values) == 257
    assert values[0] == 0.0


def _malformed_graph(adata, science, kind):
    key = "custom_connectivities"
    graph = adata.obsp[key].copy()
    if kind == "shape":
        adata.obsp._data[key] = science.sparse.csr_matrix((adata.n_obs - 1, adata.n_obs - 1))
        return
    if kind == "nonfinite":
        graph.data[0] = science.np.nan
    elif kind == "negative":
        graph.data[0] = -1.0
    elif kind == "diagonal":
        graph = graph.tolil()
        graph[0, 0] = 1.0
        graph = graph.tocsr()
    elif kind == "asymmetric":
        graph = graph.tolil()
        graph[0, 1] = 0.0
        graph = graph.tocsr()
        graph.eliminate_zeros()
    elif kind == "empty":
        graph = science.sparse.csr_matrix(graph.shape, dtype=float)
    else:  # pragma: no cover - test helper guard
        raise AssertionError(kind)
    adata.obsp[key] = graph


@pytest.mark.parametrize(
    ("kind", "message"),
    [
        ("shape", "shape"),
        ("nonfinite", "finite"),
        ("negative", "negative|non-negative"),
        ("diagonal", "diagonal|self"),
        ("asymmetric", "symmetric|undirected"),
        ("empty", "edge|empty"),
    ],
)
def test_graph_validation_rejects_malformed_connectivities(science, kind, message):
    adata = _adata_with_named_graph(science)
    _malformed_graph(adata, science, kind)
    with pytest.raises(ValueError, match=message):
        _execute(adata)


@pytest.mark.parametrize("missing", ["neighbors", "pointer", "obsp"])
def test_named_graph_contract_requires_metadata_pointer_and_matrix(science, missing):
    adata = _adata_with_named_graph(science)
    if missing == "neighbors":
        del adata.uns["custom_neighbors"]
    elif missing == "pointer":
        del adata.uns["custom_neighbors"]["connectivities_key"]
    else:
        del adata.obsp["custom_connectivities"]
    with pytest.raises(ValueError, match="neighbor|connectivit|graph"):
        _execute(adata)


@pytest.mark.parametrize("container_name", ["obs", "uns"])
def test_all_output_collisions_are_preflighted_before_any_backend_run(science, monkeypatch, container_name):
    adata = _adata_with_named_graph(science)
    if container_name == "obs":
        adata.obs["leiden_1"] = science.pd.Categorical(["existing"] * adata.n_obs)
    else:
        adata.uns["leiden_1"] = {"existing": True}

    def unexpected_backend_call(*args, **kwargs):
        del args, kwargs
        raise AssertionError("Leiden backend ran before collision preflight completed")

    monkeypatch.setattr(science.sc.tl, "leiden", unexpected_backend_call)
    with pytest.raises(ValueError, match="exist|collision|overwrite"):
        _execute(adata, resolutions="0.5,1.0")
    assert "leiden_0_5" not in adata.obs
    assert "leiden_0_5" not in adata.uns


@pytest.mark.parametrize("overwrite_existing", [False, True])
def test_runtime_and_generated_sweep_reserve_resolved_graph_metadata_key(science, overwrite_existing):
    adata = _adata_with_named_graph(science)
    adata.uns["leiden_1"] = adata.uns.pop("custom_neighbors")
    metadata_snapshot = copy.deepcopy(adata.uns["leiden_1"])
    graph_snapshot = adata.obsp["custom_connectivities"].copy()
    parameters = {
        "resolutions": [1.0],
        "keys": ["leiden_1"],
        "neighbors_key": "leiden_1",
        "n_iterations": 2,
        "stability_repeats": 2,
        "random_seed": 0,
        "overwrite_existing": overwrite_existing,
    }
    namespace: dict[str, object] = {}
    exec(_leiden_resolution_sweep_code(parameters), namespace)

    runners = (
        lambda value: _execute(
            value,
            resolutions="1.0",
            key_prefix="leiden",
            neighbors_key="leiden_1",
            stability_repeats=2,
            overwrite_existing=overwrite_existing,
        ),
        namespace["leiden_resolution_sweep"],
    )
    for runner in runners:
        with pytest.raises(ValueError, match="cannot equal.*neighbors_key|graph metadata"):
            runner(adata)
        assert adata.uns["leiden_1"] == metadata_snapshot
        delta = adata.obsp["custom_connectivities"] - graph_snapshot
        delta.eliminate_zeros()
        assert delta.nnz == 0
        graph = resolve_named_graph(adata, "leiden_1", operation="post-collision downstream")
        assert graph.connectivities_key == "custom_connectivities"


def test_memberships_metrics_summary_and_input_immutability(science):
    adata = _adata_with_named_graph(science)
    original_x = adata.X.copy()
    original_obs = adata.obs.copy(deep=True)
    original_uns = copy.deepcopy(adata.uns)
    original_graph = adata.obsp["custom_connectivities"].copy()

    output, table_result, report, code = _execute(
        adata,
        resolutions="0.4,0.8,1.2",
        stability_repeats=3,
        random_seed=11,
    )
    table = _metric_table(table_result)
    summary = report.summary

    assert list(table.columns) == METRIC_COLUMNS
    assert table["resolution"].tolist() == [0.4, 0.8, 1.2]
    assert table["stability_repeats"].tolist() == [3, 3, 3]
    for row in table.to_dict(orient="records"):
        key = row["key"]
        assert key in output.obs
        assert str(output.obs[key].dtype) == "category"
        assert not output.obs[key].isna().any()
        sizes = output.obs[key].value_counts(sort=False)
        assert int(sizes.sum()) == adata.n_obs
        assert row["n_clusters"] == len(sizes)
        assert row["min_cluster_size"] == int(sizes.min())
        assert row["median_cluster_size"] == pytest.approx(float(science.np.median(sizes)))
        assert row["max_cluster_size"] == int(sizes.max())
        assert row["smallest_cluster_fraction"] == pytest.approx(float(sizes.min() / adata.n_obs))
        assert row["largest_cluster_fraction"] == pytest.approx(float(sizes.max() / adata.n_obs))
        assert row["singleton_clusters"] == int((sizes == 1).sum())
        assert row["modularity"] == pytest.approx(float(output.uns[key]["modularity"]))
        independent_modularity = _independent_modularity(
            science,
            original_graph,
            output.obs[key].to_numpy(),
            resolution=row["resolution"],
        )
        assert science.np.isfinite(row["modularity"])
        assert -1.0 <= row["modularity"] <= 1.0
        assert row["modularity"] == pytest.approx(independent_modularity)
        assert -1.0 <= row["stability_min_ari"] <= row["stability_mean_ari"] <= 1.0
        assert row["stability_mean_ari"] <= row["stability_max_ari"] <= 1.0

    keys = _membership_keys(table)
    assert science.pd.isna(table.iloc[0]["adjacent_previous_resolution"])
    assert science.pd.isna(table.iloc[0]["adjacent_resolution_ari"])
    for index in range(1, len(table)):
        expected = adjusted_rand_score(
            output.obs[keys[index - 1]],
            output.obs[keys[index]],
        )
        assert table.iloc[index]["adjacent_previous_resolution"] == table.iloc[index - 1]["resolution"]
        assert table.iloc[index]["adjacent_resolution_ari"] == pytest.approx(expected)

    expected_cluster_sizes = {
        key: {str(label): int(count) for label, count in output.obs[key].value_counts(sort=False).items()}
        for key in keys
    }
    assert summary["key_results"]["cluster_sizes_by_resolution"] == expected_cluster_sizes
    assert summary["parameters"]["fixed_policy"] == {
        "flavor": "igraph",
        "directed": False,
        "use_weights": True,
        "objective_function": "modularity",
    }
    assert summary["key_results"]["graph"]["connectivities_key"] == "custom_connectivities"
    assert summary["key_results"]["graph"]["isolated_cells"] == 0
    assert isinstance(code, str) and code.strip()
    assert not any(str(key).startswith("__openbio_leiden_stability_") for key in output.obs)
    assert not any(str(key).startswith("__openbio_leiden_stability_") for key in output.uns)

    science.np.testing.assert_array_equal(adata.X, original_x)
    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    assert adata.uns == original_uns
    difference = adata.obsp["custom_connectivities"] - original_graph
    assert difference.nnz == 0
    assert resolve_named_graph(output, "custom_neighbors", operation="downstream sweep result")


def test_within_resolution_stability_matches_disclosed_seed_schedule(science):
    adata = _adata_with_named_graph(science)
    output, table_result, report, _ = _execute(
        adata,
        resolutions="0.6,1.0",
        n_iterations=2,
        stability_repeats=3,
        random_seed=23,
    )
    table = _metric_table(table_result)
    seeds = report.summary["key_results"]["stability_seeds"]
    assert seeds == [23, 24, 25]

    for row in table.to_dict(orient="records"):
        memberships = []
        for seed in seeds:
            temporary = adata.copy()
            science.sc.tl.leiden(
                temporary,
                resolution=row["resolution"],
                key_added="_manual_leiden",
                adjacency=adata.obsp["custom_connectivities"],
                random_state=seed,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                use_weights=True,
                objective_function="modularity",
            )
            memberships.append(temporary.obs["_manual_leiden"].to_numpy(copy=True))
        direct = _pairwise_ari(science, memberships)
        assert row["stability_mean_ari"] == pytest.approx(direct["mean"])
        assert row["stability_min_ari"] == pytest.approx(direct["min"])
        assert row["stability_max_ari"] == pytest.approx(direct["max"])
        assert row["key"] in output.obs


@pytest.mark.parametrize(
    ("malformation", "message"),
    [
        ("noncategorical", "categorical|membership"),
        ("nonfinite_modularity", "modularity"),
        ("out_of_range_modularity", "modularity"),
    ],
)
def test_runtime_and_generated_code_reject_malformed_backend_results(
    science,
    monkeypatch,
    malformation,
    message,
):
    _, _, _, code = _execute(
        _adata_with_named_graph(science),
        resolutions="0.5,1.0",
        stability_repeats=2,
    )
    namespace: dict[str, object] = {}
    exec(compile(code, "<leiden-resolution-sweep-code>", "exec"), namespace)

    def malformed_leiden(adata, *, key_added, **kwargs):
        del kwargs
        if malformation == "noncategorical":
            adata.obs[key_added] = ["0"] * adata.n_obs
        elif malformation == "unused_category":
            adata.obs[key_added] = science.pd.Categorical(
                ["0"] * (adata.n_obs - 1) + ["1"],
                categories=["0", "1", "2"],
            )
        elif malformation == "noncanonical_labels":
            adata.obs[key_added] = science.pd.Categorical(
                ["0", "2"] * (adata.n_obs // 2),
                categories=["0", "2"],
            )
        else:
            adata.obs[key_added] = science.pd.Categorical(
                ["0", "1"] * (adata.n_obs // 2),
                categories=["0", "1"],
            )
        modularity = 0.25
        if malformation == "nonfinite_modularity":
            modularity = science.np.nan
        elif malformation == "out_of_range_modularity":
            modularity = 2.0
        adata.uns[key_added] = {"params": {}, "modularity": modularity}

    monkeypatch.setattr(science.sc.tl, "leiden", malformed_leiden)
    runners = [
        lambda adata: _execute(adata, resolutions="0.5,1.0", stability_repeats=2),
        namespace["leiden_resolution_sweep"],
    ]
    for runner in runners:
        with pytest.raises((RuntimeError, ValueError), match=message):
            runner(_adata_with_named_graph(science))


@pytest.mark.parametrize("label_kind", ["unused_category", "noncanonical_labels"])
def test_runtime_and_generated_code_normalize_complete_backend_labels(science, monkeypatch, label_kind):
    _, _, _, code = _execute(
        _adata_with_named_graph(science),
        resolutions="0.5,1.0",
        stability_repeats=2,
    )
    namespace: dict[str, object] = {}
    exec(compile(code, "<leiden-resolution-sweep-code>", "exec"), namespace)

    def backend(adata, *, key_added, resolution, **kwargs):
        del kwargs
        raw = ["A"] * (adata.n_obs - 1) + ["B"]
        categories = ["A", "B", "unused"] if label_kind == "unused_category" else ["A", "B"]
        if label_kind == "noncanonical_labels":
            raw = ["A", "B"] * (adata.n_obs // 2)
        adata.obs[key_added] = science.pd.Categorical(raw, categories=categories)
        adjacency = adata.obsp["custom_connectivities"]
        upper = science.sparse.triu(adjacency, k=1).tocoo()
        graph = igraph.Graph(
            n=adata.n_obs,
            edges=list(zip(upper.row.tolist(), upper.col.tolist(), strict=True)),
            directed=False,
        )
        graph.es["weight"] = [float(value) for value in upper.data]
        mapping = {label: index for index, label in enumerate(dict.fromkeys(raw))}
        modularity = graph.modularity(
            [mapping[label] for label in raw],
            weights=graph.es["weight"],
            resolution=resolution,
            directed=False,
        )
        adata.uns[key_added] = {"params": {}, "modularity": modularity}

    monkeypatch.setattr(science.sc.tl, "leiden", backend)
    runtime = _execute(
        _adata_with_named_graph(science), resolutions="0.5,1.0", stability_repeats=2
    )[0]
    generated, _ = namespace["leiden_resolution_sweep"](_adata_with_named_graph(science))
    for key in ("leiden_0_5", "leiden_1"):
        assert list(runtime.obs[key].cat.categories) == ["0", "1"]
        science.pd.testing.assert_series_equal(generated.obs[key], runtime.obs[key])


def test_summary_is_strict_json_and_makes_no_best_or_truth_claim(science):
    _, table_result, report, _ = _execute(
        _adata_with_named_graph(science),
        resolutions="0.5,1.0",
        stability_repeats=2,
    )
    table = _metric_table(table_result)
    summary = report.summary

    json.dumps(summary, allow_nan=False)
    assert summary["node_id"] == "OpenBioSingleCellLeidenResolutionSweep"
    expected_metrics = table.astype(object).where(table.notna(), None).to_dict(orient="records")
    assert summary["key_results"]["resolution_metrics"] == expected_metrics
    assert summary["references"]
    assert "scanpy" in summary["software_versions"]
    assert "igraph" in summary["software_versions"]
    serialized = json.dumps(summary, sort_keys=True).lower()
    assert "best_resolution" not in serialized
    assert "optimal_resolution" not in serialized
    results_text = str(summary["results"]).lower()
    assert "true number of cell types" not in results_text
    assert "ground truth" not in results_text
    assert "condition difference" not in results_text


def test_generated_code_reproduces_memberships_and_resolution_metrics(science):
    adata = _adata_with_named_graph(science)
    runtime_output, runtime_table_result, _, code = _execute(
        adata,
        resolutions="0.4,0.8,1.2",
        stability_repeats=3,
        random_seed=31,
    )
    runtime_table = _metric_table(runtime_table_result)

    namespace: dict[str, object] = {}
    compile(code, "<leiden-resolution-sweep-code>", "exec")
    exec(code, namespace)
    assert "leiden_resolution_sweep" in namespace
    generated_output, generated_table = namespace["leiden_resolution_sweep"](adata)

    science.pd.testing.assert_frame_equal(generated_table, runtime_table)
    for key in _membership_keys(runtime_table):
        science.pd.testing.assert_series_equal(generated_output.obs[key], runtime_output.obs[key])
    assert not any(str(key).startswith("__openbio_leiden_stability_") for key in runtime_output.obs)
    assert not any(str(key).startswith("__openbio_leiden_stability_") for key in runtime_output.uns)
    assert not any(str(key).startswith("__openbio_leiden_stability_") for key in generated_output.obs)
    assert not any(str(key).startswith("__openbio_leiden_stability_") for key in generated_output.uns)
    assert "leiden_0_4" not in adata.obs
