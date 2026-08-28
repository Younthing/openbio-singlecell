from __future__ import annotations

import json

import pytest
from scipy import stats

from openbio_singlecell.nodes_results import (
    PCA_METADATA_COLUMNS,
    OpenBioSingleCellPCAMetadataAssociations,
)


def _adata(science):
    samples = [f"sample_{index}" for index in range(1, 9)]
    condition_by_sample = ["control"] * 4 + ["treated"] * 4
    age_by_sample = [20.0, 28.0, 35.0, 43.0, 51.0, 60.0, 68.0, 77.0]
    pc1_by_sample = [0.0, 1.5, 0.5, 2.5, 5.0, 8.0, 6.5, 10.0]
    pc2_by_sample = [8.0, 6.0, 7.0, 5.0, 4.0, 2.0, 3.0, 1.0]
    cell_offsets = [-0.2, 0.0, 0.2]

    sample_values = []
    condition_values = []
    age_values = []
    scores = []
    for sample, condition, age, pc1, pc2 in zip(
        samples,
        condition_by_sample,
        age_by_sample,
        pc1_by_sample,
        pc2_by_sample,
        strict=True,
    ):
        for offset in cell_offsets:
            sample_values.append(sample)
            condition_values.append(condition)
            age_values.append(age + offset)
            scores.append([pc1 + offset, pc2 - offset])

    obs_names = [f"cell_{index}" for index in range(len(scores))]
    adata = science.ad.AnnData(
        X=science.np.arange(len(scores) * 3, dtype=float).reshape(len(scores), 3),
        obs=science.pd.DataFrame(
            {
                "sample": sample_values,
                "condition": condition_values,
                "age": age_values,
            },
            index=obs_names,
        ),
        var=science.pd.DataFrame(index=["gene_1", "gene_2", "gene_3"]),
    )
    adata.obsm["X_pca"] = science.np.asarray(scores, dtype=float)
    return adata


def _execute(adata, **overrides):
    arguments = {
        "use_rep": "X_pca",
        "sample_key": "sample",
        "categorical_obs_keys": "condition",
        "continuous_obs_keys": "age",
        "alpha": 0.05,
    }
    arguments.update(overrides)
    return OpenBioSingleCellPCAMetadataAssociations.execute(adata, **arguments).result


def _sample_means(adata, science):
    frame = science.pd.DataFrame(
        adata.obsm["X_pca"],
        index=adata.obs_names,
        columns=["PC1", "PC2"],
    )
    frame["sample"] = adata.obs["sample"].to_numpy()
    return frame.groupby("sample", sort=False)[["PC1", "PC2"]].mean()


def test_schema_exposes_atomic_sample_level_contract():
    schema = OpenBioSingleCellPCAMetadataAssociations.GET_SCHEMA()

    assert [item.id for item in schema.inputs] == [
        "adata",
        "use_rep",
        "sample_key",
        "categorical_obs_keys",
        "continuous_obs_keys",
        "alpha",
    ]
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("table", "OPENBIO_SINGLE_CELL_TABLE"),
        ("summary", "OPENBIO_SINGLE_CELL_SUMMARY"),
        ("code", "STRING"),
    ]


def test_results_match_scipy_and_generated_code(science):
    adata = _adata(science)
    original_obs = adata.obs.copy(deep=True)
    original_scores = adata.obsm["X_pca"].copy()

    table_result, report, code = _execute(adata)
    table = table_result.table

    assert table_result.kind == "table"
    assert table.columns.tolist() == PCA_METADATA_COLUMNS
    assert set(table["metadata_type"]) == {"categorical", "continuous"}
    assert len(table) == 4
    assert table[["statistic", "effect_size", "p_value", "p_adjusted"]].notna().all().all()

    sample_scores = _sample_means(adata, science)
    sample_obs = adata.obs.groupby("sample", sort=False).agg(condition=("condition", "first"), age=("age", "mean"))

    continuous = table[(table["component"] == "PC2") & (table["metadata"] == "age")].iloc[0]
    expected_spearman = stats.spearmanr(sample_obs["age"], sample_scores["PC2"], alternative="two-sided")
    assert continuous["statistic_name"] == "spearman_rho"
    assert continuous["estimate_name"] == "spearman_rho"
    assert continuous["effect_size_name"] == "rho_squared"
    assert continuous["statistic"] == pytest.approx(expected_spearman.statistic)
    assert continuous["estimate"] == pytest.approx(expected_spearman.statistic)
    assert continuous["effect_size"] == pytest.approx(expected_spearman.statistic**2)
    assert continuous["p_value"] == pytest.approx(expected_spearman.pvalue)

    categorical = table[(table["component"] == "PC1") & (table["metadata"] == "condition")].iloc[0]
    groups = [
        sample_scores.loc[sample_obs["condition"] == level, "PC1"].to_numpy()
        for level in sample_obs["condition"].drop_duplicates()
    ]
    expected_anova = stats.f_oneway(*groups)
    all_values = sample_scores["PC1"].to_numpy()
    grand_mean = all_values.mean()
    ss_total = science.np.square(all_values - grand_mean).sum()
    ss_between = sum(group.size * (group.mean() - grand_mean) ** 2 for group in groups)
    assert categorical["statistic_name"] == "anova_f"
    assert categorical["effect_size_name"] == "eta_squared"
    assert categorical["statistic"] == pytest.approx(expected_anova.statistic)
    assert categorical["effect_size"] == pytest.approx(ss_between / ss_total)
    assert categorical["p_value"] == pytest.approx(expected_anova.pvalue)

    expected_adjusted = stats.false_discovery_control(table["p_value"].to_numpy(), method="bh")
    science.np.testing.assert_allclose(table["p_adjusted"], expected_adjusted)
    assert table["significant"].tolist() == (expected_adjusted <= 0.05).tolist()

    assert report.kind == "summary"
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellPCAMetadataAssociations"
    assert summary["parameters"]["sample_aggregation"] == "unweighted_mean"
    assert summary["parameters"]["continuous_test"] == "scipy_spearmanr_two_sided"
    assert summary["parameters"]["p_adjust_scope"] == "all_valid_pc_metadata_pairs"
    assert "spearman_permutation_seed" not in summary["parameters"]
    assert "spearman_p_value_policy" not in summary["parameters"]
    assert summary["key_results"]["samples"] == 8
    assert summary["key_results"]["requested_hypotheses"] == 4
    assert summary["key_results"]["valid_hypotheses"] == 4
    assert summary["key_results"]["skipped_hypothesis_count"] == 0
    assert summary["key_results"]["metadata_diagnostics"]["age"]["test"] == "scipy.stats.spearmanr"
    assert summary["key_results"]["metadata_diagnostics"]["age"]["alternative"] == "two-sided"
    assert "scipy's two-sided spearman" in summary["methods"].lower()
    assert "exploratory" in summary["results"].lower()
    assert "not causal" in summary["results"].lower()
    assert any("Condition" in limitation for limitation in summary["limitations"])
    assert {reference["kind"] for reference in summary["references"]} >= {"method", "software", "practice"}
    assert any(reference["doi"] == "10.2307/1412159" for reference in summary["references"])
    assert not any(reference["doi"] == "10.2202/1544-6115.1585" for reference in summary["references"])
    assert {"anndata", "numpy", "pandas", "scipy"} <= set(summary["software_versions"])
    json.dumps(summary, allow_nan=False)
    assert "permutation" not in code.lower()

    namespace = {}
    compile(code, "<OpenBioSingleCellPCAMetadataAssociations-code>", "exec")
    exec(code, namespace)
    generated = namespace["pca_metadata_associations"](adata)
    science.pd.testing.assert_frame_equal(generated, table)

    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    science.np.testing.assert_array_equal(adata.obsm["X_pca"], original_scores)


def test_samples_are_equal_weight_regardless_of_cell_multiplicity(science):
    balanced = _adata(science)
    expected = _execute(balanced)[0].table

    sample_values = balanced.obs["sample"].to_numpy()
    indices = []
    for index, sample in enumerate(sample_values):
        repeat = 9 if sample == "sample_1" else (4 if sample == "sample_8" else 1)
        indices.extend([index] * repeat)
    expanded_names = [f"expanded_{index}" for index in range(len(indices))]
    expanded_obs = balanced.obs.iloc[indices].copy()
    expanded_obs.index = expanded_names
    unbalanced = science.ad.AnnData(
        X=science.np.asarray(balanced.X)[indices],
        obs=expanded_obs,
        var=balanced.var.copy(),
    )
    unbalanced.obsm["X_pca"] = balanced.obsm["X_pca"][indices]
    observed = _execute(unbalanced)[0].table

    invariant_columns = [
        "component",
        "metadata",
        "metadata_type",
        "n_samples_total",
        "n_samples_analyzed",
        "n_samples_missing",
        "n_levels",
        "statistic_name",
        "statistic",
        "estimate_name",
        "estimate",
        "effect_size_name",
        "effect_size",
        "p_value",
        "p_adjusted",
        "significant",
    ]
    science.pd.testing.assert_frame_equal(observed[invariant_columns], expected[invariant_columns])
    assert observed["n_cells_total"].unique().tolist() == [unbalanced.n_obs]


def test_three_sample_spearman_uses_scipy_two_sided_p_value(science):
    adata = science.ad.AnnData(
        X=science.np.ones((3, 1), dtype=float),
        obs=science.pd.DataFrame(
            {"sample": ["sample_1", "sample_2", "sample_3"], "age": [1.0, 2.0, 3.0]},
            index=["cell_1", "cell_2", "cell_3"],
        ),
        var=science.pd.DataFrame(index=["gene_1"]),
    )
    adata.obsm["X_pca"] = science.np.asarray([[1.0], [2.0], [3.0]])

    table, report, code = _execute(adata, categorical_obs_keys="", continuous_obs_keys="age")
    row = table.table.iloc[0]
    expected = stats.spearmanr([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], alternative="two-sided")

    assert row["statistic"] == pytest.approx(expected.statistic)
    assert row["p_value"] == pytest.approx(expected.pvalue)
    assert row["p_adjusted"] == pytest.approx(expected.pvalue)
    assert bool(row["significant"]) == bool(expected.pvalue <= 0.05)
    diagnostic = report.summary["key_results"]["metadata_diagnostics"]["age"]
    assert diagnostic["test"] == "scipy.stats.spearmanr"
    assert diagnostic["alternative"] == "two-sided"

    namespace = {}
    exec(code, namespace)
    science.pd.testing.assert_frame_equal(namespace["pca_metadata_associations"](adata), table.table)


def test_nonmonotonic_spearman_matches_scipy_and_code(science):
    sample_names = [f"sample_{index}" for index in range(10)]
    adata = science.ad.AnnData(
        X=science.np.ones((10, 1), dtype=float),
        obs=science.pd.DataFrame(
            {"sample": sample_names, "age": science.np.arange(10, dtype=float)},
            index=[f"cell_{index}" for index in range(10)],
        ),
        var=science.pd.DataFrame(index=["gene_1"]),
    )
    adata.obsm["X_pca"] = science.np.asarray(
        [[0.0], [2.0], [1.0], [4.0], [3.0], [6.0], [5.0], [8.0], [7.0], [9.0]],
        dtype=float,
    )

    first, first_report, code = _execute(adata, categorical_obs_keys="", continuous_obs_keys="age")
    row = first.table.iloc[0]
    expected = stats.spearmanr(
        science.np.arange(10, dtype=float),
        adata.obsm["X_pca"][:, 0],
        alternative="two-sided",
    )
    assert row["statistic"] == pytest.approx(expected.statistic)
    assert row["p_value"] == pytest.approx(expected.pvalue)
    diagnostic = first_report.summary["key_results"]["metadata_diagnostics"]["age"]
    assert diagnostic == {
        "metadata_type": "continuous",
        "n_cells_total": 10,
        "n_cells_missing": 0,
        "n_samples_total": 10,
        "n_samples_analyzed": 10,
        "n_samples_missing": 0,
        "test": "scipy.stats.spearmanr",
        "alternative": "two-sided",
    }
    assert any("asymptotic" in limitation.lower() for limitation in first_report.summary["limitations"])

    namespace = {}
    exec(code, namespace)
    science.pd.testing.assert_frame_equal(namespace["pca_metadata_associations"](adata), first.table)


def test_missing_values_are_disclosed_at_cell_and_sample_levels(science):
    adata = _adata(science)
    sample_one = adata.obs["sample"] == "sample_1"
    sample_two_first = adata.obs.index[adata.obs["sample"] == "sample_2"][0]
    sample_eight = adata.obs["sample"] == "sample_8"
    adata.obs.loc[sample_one, "age"] = science.np.nan
    adata.obs.loc[sample_two_first, "age"] = science.np.nan
    adata.obs.loc[sample_eight, "condition"] = None

    table, report, _ = _execute(adata)
    continuous = table.table[table.table["metadata"] == "age"]
    categorical = table.table[table.table["metadata"] == "condition"]

    assert continuous["n_cells_missing"].unique().tolist() == [4]
    assert continuous["n_samples_analyzed"].unique().tolist() == [7]
    assert continuous["n_samples_missing"].unique().tolist() == [1]
    assert categorical["n_cells_missing"].unique().tolist() == [3]
    assert categorical["n_samples_analyzed"].unique().tolist() == [7]
    assert categorical["n_samples_missing"].unique().tolist() == [1]
    assert "analyzed Samples=7" in report.summary["results"]
    assert report.summary["key_results"]["metadata_diagnostics"]["condition"]["group_sizes"] == {
        "control": 4,
        "treated": 3,
    }


def test_sparse_representation_and_degenerate_components_are_handled(science):
    adata = _adata(science)
    scores = adata.obsm["X_pca"]
    adata.obsm["X_pca"] = science.sparse.csr_matrix(science.np.column_stack([scores, science.np.ones(adata.n_obs)]))

    table, report, code = _execute(adata)

    assert set(table.table["component"]) == {"PC1", "PC2"}
    skipped = report.summary["key_results"]["skipped_hypotheses"]
    assert {(item["component"], item["metadata"]) for item in skipped} == {
        ("PC3", "age"),
        ("PC3", "condition"),
    }
    assert any("Skipped 2" in warning for warning in report.summary["warnings"])
    assert report.summary["key_results"]["requested_hypotheses"] == 6
    assert report.summary["key_results"]["valid_hypotheses"] == 4
    assert report.summary["key_results"]["skipped_hypothesis_count"] == 2
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="Skipped 2"):
        generated = namespace["pca_metadata_associations"](adata)
    science.pd.testing.assert_frame_equal(generated, table.table)


@pytest.mark.parametrize("alpha", [0.0, 1.0])
def test_alpha_inclusive_boundary_runs_with_warning(science, alpha):
    table, report, code = _execute(_adata(science), alpha=alpha)
    assert any(f"alpha={alpha:g}" in warning for warning in report.summary["warnings"])
    assert table.table["significant"].tolist() == (table.table["p_adjusted"] <= alpha).tolist()
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match=f"alpha={alpha:g}"):
        generated = namespace["pca_metadata_associations"](_adata(science))
    science.pd.testing.assert_frame_equal(generated, table.table)


@pytest.mark.parametrize("alpha", [float("nan"), float("inf"), -0.1, 1.1])
def test_alpha_must_be_finite_and_inside_inclusive_unit_interval(science, alpha):
    with pytest.raises(ValueError, match="between zero and one, inclusive"):
        _execute(_adata(science), alpha=alpha)


def test_rejects_invalid_metadata_and_representation(science):
    adata = _adata(science)

    with pytest.raises(ValueError, match="at least one categorical or continuous"):
        _execute(adata, categorical_obs_keys="", continuous_obs_keys="")
    with pytest.raises(ValueError, match="both categorical and continuous"):
        _execute(adata, categorical_obs_keys="condition", continuous_obs_keys="condition")
    with pytest.raises(ValueError, match="columns not found"):
        _execute(adata, continuous_obs_keys="missing")
    with pytest.raises(ValueError, match="representation not found"):
        _execute(adata, use_rep="missing")

    invalid = adata.copy()
    invalid.obsm["X_pca"][0, 0] = science.np.inf
    with pytest.raises(ValueError, match="non-finite"):
        _execute(invalid)

    invalid = adata.copy()
    invalid.obs["text"] = ["x"] * invalid.n_obs
    with pytest.raises(TypeError, match="numeric dtype"):
        _execute(invalid, continuous_obs_keys="text")

    invalid = adata.copy()
    invalid.obs["constant"] = 1.0
    table, report, _ = _execute(invalid, continuous_obs_keys="constant")
    assert set(table.table["metadata"]) == {"condition"}
    assert any("constant" in warning.lower() for warning in report.summary["warnings"])


def test_rejects_sample_violations_and_insufficient_replication(science):
    adata = _adata(science)

    invalid = adata.copy()
    invalid.obs.iloc[0, invalid.obs.columns.get_loc("sample")] = None
    with pytest.raises(ValueError, match="missing or blank"):
        _execute(invalid)

    invalid = adata.copy()
    invalid.obs.iloc[0, invalid.obs.columns.get_loc("condition")] = "treated"
    with pytest.raises(ValueError, match="constant within each Sample"):
        _execute(invalid)

    invalid = adata.copy()
    invalid.obs.loc[invalid.obs["sample"] != "sample_1", "condition"] = "treated"
    table, report, _ = _execute(invalid)
    assert not table.table.empty
    assert any("singleton Sample group" in warning for warning in report.summary["warnings"])

    invalid = adata.copy()
    invalid.obs_names = ["duplicate"] * invalid.n_obs
    with pytest.raises(ValueError, match="unique observation names"):
        _execute(invalid)

    two_samples = adata[adata.obs["sample"].isin(["sample_1", "sample_2"])].copy()
    with pytest.raises(ValueError, match="at least three independent Samples"):
        _execute(two_samples)


def test_boolean_continuous_metadata_runs_with_warning(science):
    adata = _adata(science)
    adata.obs["flag"] = adata.obs["condition"].eq("treated")

    table, report, code = _execute(
        adata,
        categorical_obs_keys="",
        continuous_obs_keys="flag",
    )

    assert set(table.table["metadata"]) == {"flag"}
    assert any("boolean" in warning.lower() and "0/1" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="boolean"):
        generated = namespace["pca_metadata_associations"](adata)
    science.pd.testing.assert_frame_equal(generated, table.table)


def test_invalid_metadata_variables_are_skipped_independently(science):
    adata = _adata(science)
    adata.obs["mostly_missing"] = science.np.nan
    for sample, value in (("sample_1", 1.0), ("sample_2", 2.0)):
        adata.obs.loc[adata.obs["sample"] == sample, "mostly_missing"] = value
    adata.obs["single_group"] = "one"

    table, report, _ = _execute(
        adata,
        categorical_obs_keys="condition,single_group",
        continuous_obs_keys="age,mostly_missing",
    )

    assert set(table.table["metadata"]) == {"age", "condition"}
    skipped = report.summary["key_results"]["skipped_hypotheses"]
    assert {(item["component"], item["metadata"]) for item in skipped} >= {
        ("PC1", "mostly_missing"),
        ("PC2", "mostly_missing"),
        ("PC1", "single_group"),
        ("PC2", "single_group"),
    }
    assert any("mostly_missing" in warning and "skipped" in warning for warning in report.summary["warnings"])
    assert any("single_group" in warning and "skipped" in warning for warning in report.summary["warnings"])


def test_all_degenerate_components_fail_instead_of_emitting_nonfinite_results(science):
    adata = _adata(science)
    adata.obsm["X_pca"] = science.np.ones((adata.n_obs, 2))

    with pytest.raises(ValueError, match="no valid hypotheses"):
        _execute(adata)
