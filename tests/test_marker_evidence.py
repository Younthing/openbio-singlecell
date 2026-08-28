from __future__ import annotations

import copy
import json
from dataclasses import replace

import pytest
from statsmodels.stats.multitest import multipletests

import openbio_singlecell.marker_evidence as marker_evidence
from openbio_singlecell.contracts import ensure_metadata, record_history
from openbio_singlecell.marker_evidence import (
    MARKER_COLUMNS,
    MARKER_UNIVERSE_COLUMNS,
    _validate_canonical_marker_table,
)
from openbio_singlecell.nodes_results import (
    OpenBioSingleCellFilterMarkerGenes,
    OpenBioSingleCellMarkerGenes,
)


def _output_values(node_output):
    return node_output.result


def _matrix_values(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


def _logged_adata(science, *, sparse_matrix: bool = True, gene_prefix: str = "G"):
    first = science.np.array(
        [
            [12, 2, 1, 1, 2, 1],
            [11, 1, 2, 1, 1, 2],
            [10, 2, 1, 2, 1, 1],
            [9, 1, 2, 1, 2, 1],
            [11, 2, 1, 1, 1, 2],
            [10, 1, 2, 2, 1, 1],
        ],
        dtype=float,
    )
    second = science.np.array(
        [
            [1, 12, 1, 2, 1, 1],
            [2, 11, 2, 1, 1, 2],
            [1, 10, 1, 1, 2, 1],
            [2, 9, 2, 1, 1, 1],
            [1, 11, 1, 2, 2, 1],
            [2, 10, 2, 1, 1, 2],
        ],
        dtype=float,
    )
    logged = science.np.log1p(science.np.vstack([first, second]))
    matrix = science.sparse.csr_matrix(logged) if sparse_matrix else logged
    obs = science.pd.DataFrame(
        {"cluster": science.pd.Categorical(["A"] * 6 + ["B"] * 6, categories=["A", "B"])},
        index=[f"cell_{index}" for index in range(12)],
    )
    var = science.pd.DataFrame(index=[f"{gene_prefix}{index}" for index in range(6)])
    adata = science.ad.AnnData(matrix.copy(), obs=obs, var=var)
    adata.layers["log1p_norm"] = matrix.copy()
    ensure_metadata(adata, display_name="marker fixture", source={"kind": "test"})
    record_history(
        adata,
        "normalize_to_layer",
        {"source": "X", "output_layer": "log1p_norm", "transform": "log1p"},
        int(adata.n_obs),
        int(adata.n_vars),
    )
    return adata


def _rank(adata, *, n_genes: int = 0):
    return _output_values(
        OpenBioSingleCellMarkerGenes.execute(
            adata,
            groupby="cluster",
            method="wilcoxon",
            source={"source": "layer", "layer_name": "log1p_norm"},
            n_genes=n_genes,
            tie_correct=True,
            max_output_rows=10_000,
        )
    )


def _run_code(code: str, function_name: str, *arguments):
    namespace: dict[str, object] = {}
    exec(compile(code, "<openbio-marker-code>", "exec"), namespace)
    return namespace[function_name](*arguments)


def test_marker_schemas_expose_atomic_evidence_universe_summary_and_code():
    marker_schema = OpenBioSingleCellMarkerGenes.GET_SCHEMA()
    filter_schema = OpenBioSingleCellFilterMarkerGenes.GET_SCHEMA()

    assert [output.display_name for output in marker_schema.outputs] == ["table", "universe", "summary", "code"]
    assert [output.display_name for output in filter_schema.outputs] == ["table", "universe", "summary", "code"]
    assert [input_.id for input_ in filter_schema.inputs][:2] == ["table", "universe"]
    assert "max_working_memory_gib" in [input_.id for input_ in marker_schema.inputs]


@pytest.mark.parametrize("sparse_matrix", [False, True])
def test_marker_outputs_are_canonical_bound_and_strict_json(science, sparse_matrix):
    adata = _logged_adata(science, sparse_matrix=sparse_matrix)
    table, universe, summary, code = _rank(adata)

    assert list(table.table.columns) == MARKER_COLUMNS
    assert list(universe.table.columns) == MARKER_UNIVERSE_COLUMNS
    assert universe.table["gene"].tolist() == adata.var_names.tolist()
    assert universe.table["universe_rank"].tolist() == list(range(1, adata.n_vars + 1))
    assert table.parameters["analysis_fingerprint"] == universe.parameters["analysis_fingerprint"]
    assert table.parameters["ranking_fingerprint"] == universe.parameters["ranking_fingerprint"]
    assert table.parameters["universe_fingerprint"] == universe.parameters["universe_fingerprint"]
    assert table.parameters["content_fingerprint"] != universe.parameters["content_fingerprint"]
    assert len(table.table) == adata.n_vars * 2
    assert table.table.groupby("group", sort=False)["rank"].apply(list).tolist() == [
        list(range(1, adata.n_vars + 1)),
        list(range(1, adata.n_vars + 1)),
    ]
    json.dumps(summary.summary, allow_nan=False)
    assert "Cluster marker evidence" in summary.summary["results"]
    assert "normal-approximation" in summary.summary["results"]
    assert "tie correction True" in summary.summary["results"]
    assert "Condition contrast" in " ".join(summary.summary["limitations"])
    assert {"scanpy", "anndata", "numpy", "pandas", "scipy", "statsmodels"}.issubset(
        summary.summary["software_versions"]
    )
    assert any("Statsmodels" in reference["citation"] for reference in summary.summary["references"])
    assert summary.summary["parameters"]["max_working_memory_gib"] == 4.0
    assert summary.summary["key_results"]["total_output_rows"] == len(table.table) + len(universe.table)
    assert "working_memory_preflight" in summary.summary["key_results"]
    assert "rank_marker_genes" in code


def test_marker_dense_and_sparse_results_agree(science):
    dense_table = _rank(_logged_adata(science, sparse_matrix=False))[0].table
    sparse_table = _rank(_logged_adata(science, sparse_matrix=True))[0].table

    science.pd.testing.assert_frame_equal(dense_table, sparse_table, check_exact=False, rtol=1e-6, atol=1e-8)


def test_marker_recomputes_bh_over_complete_family_including_constants(science):
    adata = _logged_adata(science, sparse_matrix=False)
    logged = science.np.asarray(adata.layers["log1p_norm"]).copy()
    logged[:, -1] = 0.0
    adata.layers["log1p_norm"] = logged

    table = _rank(adata)[0].table

    for _, group_table in table.groupby("group", sort=False):
        expected = multipletests(group_table["p_value"], method="fdr_bh")[1]
        science.np.testing.assert_allclose(group_table["p_adjusted"], expected, rtol=1e-12, atol=1e-15)
        assert (group_table["p_adjusted"] >= group_table["p_value"] - 1e-12).all()


def test_marker_is_read_only_and_generated_code_is_equivalent(science):
    adata = _logged_adata(science)
    original_x = _matrix_values(adata.X, science).copy()
    original_layer = _matrix_values(adata.layers["log1p_norm"], science).copy()
    original_obs = adata.obs.copy(deep=True)
    original_uns = copy.deepcopy(adata.uns)

    table, universe, _, code = _rank(adata, n_genes=4)
    generated_table, generated_universe = _run_code(code, "rank_marker_genes", adata)

    science.pd.testing.assert_frame_equal(generated_table, table.table)
    science.pd.testing.assert_frame_equal(generated_universe, universe.table)
    science.np.testing.assert_array_equal(_matrix_values(adata.X, science), original_x)
    science.np.testing.assert_array_equal(_matrix_values(adata.layers["log1p_norm"], science), original_layer)
    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    assert adata.uns == original_uns


def test_marker_calls_scanpy_with_fixed_scientific_policy(science, monkeypatch):
    adata = _logged_adata(science)
    original = science.sc.tl.rank_genes_groups
    received = {}

    def recording_backend(*args, **kwargs):
        received.update(kwargs)
        return original(*args, **kwargs)

    monkeypatch.setattr(science.sc.tl, "rank_genes_groups", recording_backend)
    _rank(adata, n_genes=3)

    assert received["groups"] == "all"
    assert received["reference"] == "rest"
    assert received["corr_method"] == "benjamini-hochberg"
    assert received["rankby_abs"] is False
    assert received["pts"] is True
    assert received["tie_correct"] is True
    assert received["use_raw"] is False
    assert received["layer"] is None


def test_marker_backend_receives_only_variable_expression_and_budget_preflights_first(science, monkeypatch):
    adata = _logged_adata(science, sparse_matrix=False)
    logged = science.np.asarray(adata.layers["log1p_norm"]).copy()
    logged[:, -2:] = 0.0
    adata.layers["log1p_norm"] = logged
    original = science.sc.tl.rank_genes_groups
    observed_n_vars = []

    def recording_backend(work, *args, **kwargs):
        observed_n_vars.append(int(work.n_vars))
        return original(work, *args, **kwargs)

    monkeypatch.setattr(science.sc.tl, "rank_genes_groups", recording_backend)
    _rank(adata)
    assert observed_n_vars == [adata.n_vars - 2]

    observed_n_vars.clear()
    with pytest.raises(ValueError, match="exceeding max_working_memory_gib"):
        OpenBioSingleCellMarkerGenes.execute(
            adata,
            groupby="cluster",
            method="wilcoxon",
            source={"source": "layer", "layer_name": "log1p_norm"},
            n_genes=0,
            tie_correct=True,
            max_output_rows=10_000,
            max_working_memory_gib=1e-12,
        )
    assert observed_n_vars == []


def test_marker_memory_guard_precedes_identifier_and_group_materialization(science, monkeypatch):
    adata = _logged_adata(science)

    def forbidden_identifier_materialization(*_args, **_kwargs):
        raise AssertionError("identifier materialization occurred before the memory guard")

    monkeypatch.setattr(marker_evidence, "_validate_identifier_values", forbidden_identifier_materialization)
    with pytest.raises(ValueError, match="identifier/group validation buffers.*exceeding max_working_memory_gib"):
        OpenBioSingleCellMarkerGenes.execute(
            adata,
            groupby="cluster",
            method="wilcoxon",
            source={"source": "layer", "layer_name": "log1p_norm"},
            n_genes=0,
            tie_correct=True,
            max_output_rows=10_000,
            max_working_memory_gib=1e-12,
        )


@pytest.mark.parametrize("method", ["t-test", "t-test_overestim_var"])
def test_marker_supported_t_test_paths_return_complete_finite_evidence(science, method):
    adata = _logged_adata(science)
    table, universe, summary, code = _output_values(
        OpenBioSingleCellMarkerGenes.execute(
            adata,
            "cluster",
            method,
            {"source": "layer", "layer_name": "log1p_norm"},
            0,
            False,
            10_000,
        )
    )

    assert len(table.table) == len(universe.table) * 2
    assert science.np.isfinite(table.table[MARKER_COLUMNS[2:]].to_numpy(dtype=float)).all()
    assert summary.summary["parameters"]["marker_method"] == method
    generated_table, generated_universe = _run_code(code, "rank_marker_genes", adata)
    science.pd.testing.assert_frame_equal(generated_table, table.table)
    science.pd.testing.assert_frame_equal(generated_universe, universe.table)


@pytest.mark.parametrize("method", ["wilcoxon", "t-test", "t-test_overestim_var"])
@pytest.mark.parametrize("sparse_matrix", [False, True])
def test_marker_retains_zero_and_positive_constant_genes_as_neutral_hypotheses(
    science,
    method,
    sparse_matrix,
):
    adata = _logged_adata(science, sparse_matrix=False)
    logged = science.np.asarray(adata.layers["log1p_norm"]).copy()
    logged[:, -2] = science.np.log1p(3.0)
    logged[:, -1] = 0.0
    adata.layers["log1p_norm"] = science.sparse.csr_matrix(logged) if sparse_matrix else logged

    table, universe, summary, _ = _output_values(
        OpenBioSingleCellMarkerGenes.execute(
            adata,
            "cluster",
            method,
            {"source": "layer", "layer_name": "log1p_norm"},
            0,
            method == "wilcoxon",
            10_000,
            4.0,
        )
    )

    constants = table.table.loc[table.table["gene"].isin(["G4", "G5"])]
    assert len(constants) == 4
    assert (constants[["score", "log2_fold_change_approx"]] == 0.0).all().all()
    assert (constants[["p_value", "p_adjusted"]] == 1.0).all().all()
    assert set(constants.loc[constants["gene"] == "G4", "fraction_in_group"]) == {1.0}
    assert set(constants.loc[constants["gene"] == "G4", "fraction_reference"]) == {1.0}
    assert set(constants.loc[constants["gene"] == "G5", "fraction_in_group"]) == {0.0}
    assert set(constants.loc[constants["gene"] == "G5", "fraction_reference"]) == {0.0}
    assert universe.table["gene"].tolist() == adata.var_names.tolist()
    assert summary.summary["key_results"]["constant_gene_count"] == 2


@pytest.mark.parametrize("sparse_matrix", [False, True])
def test_marker_tiny_nonzero_expression_is_not_misclassified_as_constant(science, sparse_matrix):
    adata = _logged_adata(science, sparse_matrix=False)
    logged = science.np.asarray(adata.layers["log1p_norm"]).copy()
    logged[:, -1] = 0.0
    logged[0, -1] = 1e-8
    adata.layers["log1p_norm"] = science.sparse.csr_matrix(logged) if sparse_matrix else logged

    table, _, summary, code = _rank(adata)
    tiny = table.table.loc[table.table["gene"] == "G5"].set_index("group")

    assert tiny.loc["A", "fraction_in_group"] == pytest.approx(1 / 6)
    assert tiny.loc["A", "fraction_reference"] == pytest.approx(0.0)
    assert tiny.loc["B", "fraction_in_group"] == pytest.approx(0.0)
    assert tiny.loc["B", "fraction_reference"] == pytest.approx(1 / 6)
    assert not (tiny[["score", "p_value"]] == science.np.asarray([0.0, 1.0])).all(axis=1).any()
    assert summary.summary["key_results"]["constant_gene_count"] == 0

    generated_table, _ = _run_code(code, "rank_marker_genes", adata)
    science.pd.testing.assert_frame_equal(generated_table, table.table)


def test_marker_all_constant_logged_expression_returns_complete_neutral_family_without_backend(science, monkeypatch):
    adata = _logged_adata(science, sparse_matrix=False)
    adata.layers["log1p_norm"] = science.np.zeros(adata.shape, dtype=float)

    def forbidden_backend(*args, **kwargs):
        raise AssertionError("backend must not run when every gene is globally constant")

    monkeypatch.setattr(science.sc.tl, "rank_genes_groups", forbidden_backend)
    table, universe, _, _ = _rank(adata)

    assert len(table.table) == adata.n_vars * 2
    assert len(universe.table) == adata.n_vars
    assert (
        (table.table[["score", "log2_fold_change_approx", "fraction_in_group", "fraction_reference"]] == 0).all().all()
    )
    assert (table.table[["p_value", "p_adjusted"]] == 1).all().all()


def test_marker_rejects_malformed_backend_results(science, monkeypatch):
    adata = _logged_adata(science)
    original = science.sc.get.rank_genes_groups_df

    def malformed(*args, **kwargs):
        frame = original(*args, **kwargs).copy()
        frame.loc[0, "pvals_adj"] = science.np.nan
        return frame

    monkeypatch.setattr(science.sc.get, "rank_genes_groups_df", malformed)
    with pytest.raises(RuntimeError, match="backend column 'pvals_adj'.*finite"):
        _rank(adata, n_genes=3)


def test_marker_rejects_inconsistent_backend_metadata(science, monkeypatch):
    adata = _logged_adata(science)
    original = science.sc.tl.rank_genes_groups

    def malformed(work, *args, **kwargs):
        original(work, *args, **kwargs)
        work.uns[kwargs["key_added"]]["params"]["reference"] = "A"

    monkeypatch.setattr(science.sc.tl, "rank_genes_groups", malformed)
    with pytest.raises(RuntimeError, match="inconsistent 'reference'"):
        _rank(adata, n_genes=3)


def test_marker_rejects_invalid_grouping(science):
    adata = _logged_adata(science)
    adata.obs = adata.obs.assign(cluster=["A"] * adata.n_obs)

    with pytest.raises(TypeError, match="must be categorical"):
        _rank(adata, n_genes=3)


def test_marker_external_count_like_and_noninteger_sources_are_disclosed(science):
    template = _logged_adata(science, sparse_matrix=False)
    counts = science.ad.AnnData(
        X=template.X.copy(),
        obs=template.obs.copy(),
        var=template.var.copy(),
    )
    counts.layers["log1p_norm"] = science.np.rint(science.np.expm1(template.layers["log1p_norm"]))
    _count_table, _count_universe, count_summary, _count_code = _rank(counts, n_genes=3)
    assert any("count-like" in warning for warning in count_summary.summary["warnings"])

    unknown = science.ad.AnnData(
        X=template.X.copy(),
        obs=template.obs.copy(),
        var=template.var.copy(),
    )
    unknown.layers["log1p_norm"] = template.layers["log1p_norm"].copy()
    table, universe, summary, code = _rank(unknown, n_genes=3)
    assert summary.summary["parameters"]["expression_state"] == "unknown"
    assert "explicitly selected expression state" in " ".join(summary.summary["warnings"])
    assert "nonstandard abundance" in summary.summary["methods"]
    with pytest.warns(UserWarning, match="expert choice was honored"):
        generated_table, generated_universe = _run_code(code, "rank_marker_genes", unknown)
    science.pd.testing.assert_frame_equal(generated_table, table.table)
    science.pd.testing.assert_frame_equal(generated_universe, universe.table)


def test_marker_rejects_logreg_irrelevant_tie_policy_and_output_budget(science):
    adata = _logged_adata(science)
    common = {
        "adata": adata,
        "groupby": "cluster",
        "source": {"source": "layer", "layer_name": "log1p_norm"},
        "n_genes": 3,
        "max_output_rows": 10_000,
    }
    with pytest.raises(ValueError, match="Unsupported Marker Genes method"):
        OpenBioSingleCellMarkerGenes.execute(method="logreg", tie_correct=False, **common)
    with pytest.raises(ValueError, match="only applicable to the Wilcoxon"):
        OpenBioSingleCellMarkerGenes.execute(method="t-test", tie_correct=True, **common)
    with pytest.raises(ValueError, match="exceeding max_output_rows"):
        OpenBioSingleCellMarkerGenes.execute(
            adata,
            "cluster",
            "wilcoxon",
            {"source": "layer", "layer_name": "log1p_norm"},
            0,
            True,
            5,
        )
    with pytest.raises(ValueError, match=r"marker and universe outputs \(6 \+ 6\)"):
        OpenBioSingleCellMarkerGenes.execute(
            adata,
            "cluster",
            "wilcoxon",
            {"source": "layer", "layer_name": "log1p_norm"},
            3,
            True,
            11,
        )


def test_marker_normalizes_non_string_categories_and_rejects_serialization_collisions(science):
    adata = _logged_adata(science)
    adata.obs = adata.obs.assign(cluster=science.pd.Categorical([1] * 6 + [2] * 6, categories=[1, 2], ordered=True))
    original = adata.obs["cluster"].copy()

    table = _rank(adata, n_genes=3)[0]

    assert list(dict.fromkeys(table.table["group"])) == ["1", "2"]
    science.pd.testing.assert_series_equal(adata.obs["cluster"], original)

    adata.obs = adata.obs.assign(cluster=science.pd.Categorical([1] * 6 + ["1"] * 6, categories=[1, "1"]))
    with pytest.raises(ValueError, match="collide after string serialization"):
        _rank(adata, n_genes=3)


def test_marker_rejects_surrounding_whitespace_in_canonical_identifiers(science):
    adata = _logged_adata(science)
    adata.var_names = [" G0", *adata.var_names[1:].tolist()]
    with pytest.raises(ValueError, match="surrounding whitespace"):
        _rank(adata, n_genes=3)

    adata = _logged_adata(science)
    adata.obs["cluster"] = science.pd.Categorical([" A"] * 6 + ["B"] * 6, categories=[" A", "B"])
    with pytest.raises(ValueError, match="surrounding whitespace"):
        _rank(adata, n_genes=3)


def test_filter_inclusive_thresholds_preserve_order_rank_universe_and_provenance(science):
    table, universe, _, _ = _rank(_logged_adata(science))
    candidate = table.table.iloc[0]
    original_table = table.table.copy(deep=True)
    original_universe = universe.table.copy(deep=True)
    original_universe_source = copy.deepcopy(universe.source)

    filtered, returned_universe, summary, code = _output_values(
        OpenBioSingleCellFilterMarkerGenes.execute(
            table,
            universe,
            min_log2_fold_change=float(candidate["log2_fold_change_approx"]),
            min_fraction_in_group=float(candidate["fraction_in_group"]),
            max_fraction_reference=float(candidate["fraction_reference"]),
            max_p_adjusted=float(candidate["p_adjusted"]),
        )
    )

    assert ((filtered.table["group"] == candidate["group"]) & (filtered.table["gene"] == candidate["gene"])).any()
    positions = [
        original_table.index[(original_table["group"] == row.group) & (original_table["gene"] == row.gene)][0]
        for row in filtered.table.itertuples()
    ]
    assert positions == sorted(positions)
    assert returned_universe is universe
    science.pd.testing.assert_frame_equal(table.table, original_table)
    science.pd.testing.assert_frame_equal(universe.table, original_universe)
    assert universe.source == original_universe_source
    json.dumps(summary.summary, allow_nan=False)
    generated_table, generated_universe = _run_code(code, "filter_marker_genes", table.table, universe.table)
    science.pd.testing.assert_frame_equal(generated_table, filtered.table)
    science.pd.testing.assert_frame_equal(generated_universe, universe.table)

    with pytest.raises(ValueError, match="requires table operation.*marker_genes"):
        OpenBioSingleCellFilterMarkerGenes.execute(filtered, universe, -100.0, 0.0, 1.0, 1.0)


def test_filter_rejects_same_size_different_universe_and_tampered_content(science):
    table_a, universe_a, _, _ = _rank(_logged_adata(science, gene_prefix="A"))
    _, universe_b, _, _ = _rank(_logged_adata(science, gene_prefix="B"))

    with pytest.raises(ValueError, match="analysis fingerprints do not match"):
        OpenBioSingleCellFilterMarkerGenes.execute(table_a, universe_b)

    changed_frame = universe_a.table.copy()
    changed_frame.loc[0, "gene"] = "DIFFERENT"
    tampered_universe = replace(universe_a, table=changed_frame)
    with pytest.raises(ValueError, match="SHA-256 identity|current-content fingerprint"):
        OpenBioSingleCellFilterMarkerGenes.execute(table_a, tampered_universe)

    changed_table = table_a.table.copy()
    changed_table.loc[0, "gene"] = "OUTSIDE_UNIVERSE"
    tampered_table = replace(table_a, table=changed_table)
    with pytest.raises(ValueError, match="current-content fingerprint"):
        OpenBioSingleCellFilterMarkerGenes.execute(tampered_table, universe_a)


def test_filter_rejects_same_axes_and_parameters_from_different_ranking_content(science):
    adata_a = _logged_adata(science, sparse_matrix=False)
    adata_b = _logged_adata(science, sparse_matrix=False)
    changed = science.np.asarray(adata_b.layers["log1p_norm"]).copy()
    changed[:6, 2] += 1.0
    adata_b.layers["log1p_norm"] = changed
    table_a, _, _, _ = _rank(adata_a)
    table_b, universe_b, _, _ = _rank(adata_b)

    assert table_a.parameters["analysis_fingerprint"] == table_b.parameters["analysis_fingerprint"]
    assert table_a.parameters["universe_fingerprint"] == table_b.parameters["universe_fingerprint"]
    assert table_a.parameters["ranking_fingerprint"] != table_b.parameters["ranking_fingerprint"]
    with pytest.raises(ValueError, match="ranking fingerprints do not match"):
        OpenBioSingleCellFilterMarkerGenes.execute(table_a, universe_b)


def test_filter_recomputes_numeric_content_fingerprint_and_checks_bh_invariant(science):
    table, universe, _, _ = _rank(_logged_adata(science))
    changed = table.table.copy()
    changed.loc[0, "score"] += 0.25
    with pytest.raises(ValueError, match="current-content fingerprint"):
        OpenBioSingleCellFilterMarkerGenes.execute(replace(table, table=changed), universe)

    invalid = table.table.copy()
    invalid.loc[0, "p_value"] = 0.8
    invalid.loc[0, "p_adjusted"] = 0.7
    with pytest.raises(ValueError, match="cannot be smaller than raw p-values"):
        OpenBioSingleCellFilterMarkerGenes.execute(replace(table, table=invalid), universe)

    tolerated = table.table.copy()
    tolerated.loc[0, "p_value"] = 0.5
    tolerated.loc[0, "p_adjusted"] = 0.5 - 5e-13
    _validate_canonical_marker_table(tolerated, np=science.np, pd=science.pd)


def test_filter_rejects_wrong_producer_and_incomplete_direct_marker_grid(science):
    table, universe, _, _ = _rank(_logged_adata(science))
    parameters = copy.deepcopy(table.parameters)
    source = copy.deepcopy(table.source)
    parameters["producer_node_id"] = "OpenBioSingleCellUMAPPlot"
    source["parameters"]["producer_node_id"] = "OpenBioSingleCellUMAPPlot"
    wrong_producer = replace(table, parameters=parameters, source=source)

    with pytest.raises(ValueError, match="approved marker-evidence producer"):
        OpenBioSingleCellFilterMarkerGenes.execute(wrong_producer, universe)

    incomplete = replace(table, table=table.table.iloc[1:].reset_index(drop=True))
    with pytest.raises(ValueError, match="current-content fingerprint"):
        OpenBioSingleCellFilterMarkerGenes.execute(incomplete, universe)


def test_filter_valid_empty_output_is_reported_and_cannot_be_filtered_again(science):
    table, universe, _, _ = _rank(_logged_adata(science))

    empty, returned_universe, summary, code = _output_values(
        OpenBioSingleCellFilterMarkerGenes.execute(table, universe, 1_000_000.0, 1.0, 0.0, 0.0)
    )

    assert empty.table.empty
    assert list(empty.table.columns) == MARKER_COLUMNS
    assert returned_universe is universe
    assert "No marker rows passed" in " ".join(summary.summary["warnings"])
    generated_table, generated_universe = _run_code(code, "filter_marker_genes", table.table, universe.table)
    science.pd.testing.assert_frame_equal(generated_table, empty.table)
    science.pd.testing.assert_frame_equal(generated_universe, universe.table)
    with pytest.raises(ValueError, match="requires table operation.*marker_genes"):
        OpenBioSingleCellFilterMarkerGenes.execute(empty, universe, 0.0, 0.0, 1.0, 1.0)


@pytest.mark.parametrize(
    "mutator,error",
    [
        (lambda frame: frame.assign(p_adjusted="bad"), "numeric dtype"),
        (lambda frame: frame.assign(p_adjusted=float("nan")), "finite"),
        (lambda frame: science_duplicate(frame), "duplicate"),
        (lambda frame: frame.assign(rank=0), "positive integers"),
    ],
)
def test_filter_rejects_malformed_marker_tables(science, mutator, error):
    table, universe, _, _ = _rank(_logged_adata(science))
    malformed = replace(table, table=mutator(table.table.copy()))

    with pytest.raises((TypeError, ValueError), match=error):
        OpenBioSingleCellFilterMarkerGenes.execute(malformed, universe)


def science_duplicate(frame):
    duplicate = frame.copy()
    duplicate.iloc[1] = duplicate.iloc[0]
    return duplicate


@pytest.mark.parametrize(
    "threshold,value,error",
    [
        ("min_log2_fold_change", True, "finite number"),
        ("min_fraction_in_group", -0.1, "between 0 and 1"),
        ("max_fraction_reference", float("inf"), "finite"),
        ("max_p_adjusted", 1.1, "between 0 and 1"),
    ],
)
def test_filter_rejects_invalid_thresholds_before_table_values(science, threshold, value, error):
    table, universe, _, _ = _rank(_logged_adata(science))
    malformed = replace(table, table=table.table.assign(p_adjusted="also bad"))
    values = {
        "min_log2_fold_change": 1.0,
        "min_fraction_in_group": 0.25,
        "max_fraction_reference": 0.5,
        "max_p_adjusted": 0.05,
    }
    values[threshold] = value

    with pytest.raises((TypeError, ValueError), match=error):
        OpenBioSingleCellFilterMarkerGenes.execute(malformed, universe, **values)


def test_filter_generated_code_matches_runtime_failures(science):
    table, universe, _, _ = _rank(_logged_adata(science))
    filtered, _, _, code = _output_values(
        OpenBioSingleCellFilterMarkerGenes.execute(table, universe, 0.0, 0.0, 1.0, 1.0)
    )
    generated_table, generated_universe = _run_code(code, "filter_marker_genes", table.table, universe.table)
    science.pd.testing.assert_frame_equal(generated_table, filtered.table)
    science.pd.testing.assert_frame_equal(generated_universe, universe.table)

    malformed = table.table.assign(p_adjusted="not numeric")
    with pytest.raises(TypeError, match="numeric dtype"):
        _run_code(code, "filter_marker_genes", malformed, universe.table)

    changed = table.table.copy()
    changed.loc[0, "score"] += 0.5
    with pytest.raises(ValueError, match="embedded in this code"):
        _run_code(code, "filter_marker_genes", changed, universe.table)
