from __future__ import annotations

import json
import random
import sys
import types
from importlib import metadata as importlib_metadata

import pytest

from openbio_singlecell.node_types import SummaryResultType, TableResultType
from openbio_singlecell.nodes_abundance import OpenBioSingleCellMiloDifferentialAbundance


class _FakeMuData(dict):
    @property
    def mod(self):
        return self


def _milo_adata(science, *, samples_per_condition=3, include_other=True):
    sample_rows = []
    conditions = ["control", "treated"]
    sample_number = 0
    for condition in conditions:
        for replicate in range(samples_per_condition):
            sample_number += 1
            sample = f"sample_{sample_number}"
            batch = "batch_a" if replicate % 2 == 0 else "batch_b"
            for cell in range(4):
                sample_rows.append((sample, condition, batch, "T" if cell < 3 else "B"))
    if include_other:
        for cell in range(4):
            sample_rows.append(("sample_other", "other", "batch_a", "T" if cell < 2 else "B"))
    obs = science.pd.DataFrame(
        sample_rows,
        columns=["sample", "condition", "batch", "cell_type"],
        index=[f"cell_{index:03d}" for index in range(len(sample_rows))],
    )
    matrix = science.np.arange(len(obs) * 3, dtype=float).reshape(len(obs), 3)
    adata = science.ad.AnnData(matrix, obs=obs, var=science.pd.DataFrame(index=["G1", "G2", "G3"]))
    adata.obsm["X_pca"] = science.np.column_stack(
        [
            science.np.linspace(-1.0, 1.0, len(obs)),
            science.np.cos(science.np.linspace(0.0, 3.0, len(obs))),
        ]
    )
    return adata


def _install_fake_r_versions(monkeypatch):
    rpy2 = types.ModuleType("rpy2")
    rpy2.__version__ = "3.6.4"
    robjects = types.ModuleType("rpy2.robjects")

    def evaluate(expression):
        versions = {
            "getRversion": "4.6.0",
            'packageVersion("edgeR")': "4.10.3",
            'packageVersion("limma")': "3.68.0",
            'packageVersion("statmod")': "1.5.1",
        }
        for token, version in versions.items():
            if token in expression:
                return [version]
        raise AssertionError(f"unexpected R version expression: {expression}")

    robjects.r = evaluate
    monkeypatch.setitem(sys.modules, "rpy2", rpy2)
    monkeypatch.setitem(sys.modules, "rpy2.robjects", robjects)


def _install_fake_milo_backend(science, monkeypatch, *, malformed=None):
    calls = []
    pertpy = types.ModuleType("pertpy")
    pertpy.__version__ = "1.3.0"

    class Milo:
        def load(self, input, feature_key="rna"):
            calls.append(("load", input.obs_names.tolist(), feature_key))
            return _FakeMuData({feature_key: input.copy()})

        def make_nhoods(
            self,
            data,
            *,
            neighbors_key=None,
            feature_key="rna",
            prop=0.1,
            seed=0,
            copy=False,
        ):
            calls.append(("make_nhoods", neighbors_key, feature_key, prop, seed, copy))
            rna = data[feature_key]
            random.seed(seed)
            if malformed == "make_raises":
                science.np.random.random()
                raise RuntimeError("synthetic make_nhoods failure")
            n_obs = rna.n_obs
            index_positions = [0, 1, 2, 3]
            graph = science.sparse.csr_matrix(rna.obsp[f"{neighbors_key}_connectivities"])
            membership = graph[:, index_positions].copy()
            membership.data = science.np.ones_like(membership.data)
            if malformed == "membership_nonbinary":
                membership[0, 0] = 2
            elif malformed == "membership_graph_mismatch":
                membership[0, 0] = 1 - membership[0, 0]
            rna.obsm["nhoods"] = science.sparse.csr_matrix(membership)
            rna.obs["nhood_ixs_random"] = 0
            rna.obs["nhood_ixs_refined"] = 0
            rna.obs["nhood_kth_distance"] = 0.0
            expected_random_count = int(round(n_obs * prop))
            actual_random_count = expected_random_count - 1 if malformed == "random_count" else expected_random_count
            random_positions = list(range(actual_random_count))
            for position in random_positions:
                rna.obs.iloc[position, rna.obs.columns.get_loc("nhood_ixs_random")] = 1
            for offset, position in enumerate(index_positions):
                rna.obs.iloc[position, rna.obs.columns.get_loc("nhood_ixs_refined")] = 1
                rna.obs.iloc[position, rna.obs.columns.get_loc("nhood_kth_distance")] = 0.5 + 0.1 * offset
            rna.uns["nhood_neighbors_key"] = neighbors_key
            if malformed == "kth_zero":
                rna.obs.iloc[0, rna.obs.columns.get_loc("nhood_kth_distance")] = 0.0
            elif malformed == "kth_graph_mismatch":
                rna.obs.iloc[0, rna.obs.columns.get_loc("nhood_kth_distance")] = 0.55
            return None

        def count_nhoods(self, data, sample_col, feature_key="rna"):
            calls.append(("count_nhoods", sample_col, feature_key))
            rna = data[feature_key]
            membership = rna.obsm["nhoods"]
            dummies = science.pd.get_dummies(rna.obs[sample_col])
            counts = science.sparse.csr_matrix(membership).T.dot(science.sparse.csr_matrix(dummies.values)).T
            if malformed == "count_mismatch":
                counts = counts.toarray()
                counts[0, 0] += 1
            elif malformed == "count_noninteger":
                counts = counts.toarray().astype(float)
                counts[0, 0] += 0.5
            elif malformed == "count_negative":
                counts = counts.toarray()
                counts[0, 0] = -1
            milo = science.ad.AnnData(
                counts,
                obs=science.pd.DataFrame(index=dummies.columns.astype(str)),
                var=science.pd.DataFrame(index=[str(index) for index in range(membership.shape[1])]),
            )
            refined = rna.obs["nhood_ixs_refined"].to_numpy() == 1
            milo.var["index_cell"] = rna.obs_names[refined].tolist()
            milo.var["kth_distance"] = rna.obs.loc[refined, "nhood_kth_distance"].to_numpy()
            milo.uns["sample_col"] = sample_col
            if malformed == "count_axis":
                milo = milo[:-1].copy()
            elif malformed == "count_names":
                milo.obs_names = ["unknown", *milo.obs_names.tolist()[1:]]
            elif malformed == "index_cell_duplicate":
                milo.var.loc[milo.var.index[1], "index_cell"] = milo.var.iloc[0]["index_cell"]
            data["milo"] = milo
            if malformed == "count_corrupts_membership":
                data[feature_key].obsm["nhoods"][0, 0] = 1 - data[feature_key].obsm["nhoods"][0, 0]
            return data

        def da_nhoods(
            self,
            mdata,
            *,
            design,
            model_contrasts=None,
            subset_samples=None,
            add_intercept=True,
            feature_key="rna",
            reml=True,
            max_iter=50,
            tol=1e-5,
            solver="pydeseq2",
        ):
            calls.append(
                (
                    "da_nhoods",
                    design,
                    model_contrasts,
                    subset_samples,
                    add_intercept,
                    feature_key,
                    reml,
                    max_iter,
                    tol,
                    solver,
                )
            )
            if malformed == "da_raises":
                random.random()
                science.np.random.random()
                raise RuntimeError("synthetic da_nhoods failure")
            result = mdata["milo"].var
            result["logFC"] = [-1.2, 0.8, 0.2, 1.5]
            result["logCPM"] = [4.0, 5.0, 3.0, 6.0]
            result["F"] = [8.0, 7.0, 0.2, 9.0]
            result["PValue"] = [0.005, 0.01, 0.8, 0.02]
            result["FDR"] = [0.02, 0.02, 0.8, 0.02666666666666667]
            result["SpatialFDR"] = [
                0.01586309523809524,
                0.01730519480519481,
                0.8,
                0.02581113801452785,
            ]
            if malformed == "result_nonfinite":
                result.loc[result.index[0], "logFC"] = science.np.nan
            elif malformed == "result_range":
                result.loc[result.index[0], "SpatialFDR"] = 1.1
            elif malformed == "result_negative_f":
                result.loc[result.index[0], "F"] = -0.1
            elif malformed == "result_boolean":
                result["PValue"] = [True, False, True, False]
            elif malformed == "result_missing_column":
                result.drop(columns="logCPM", inplace=True)
            elif malformed == "fdr_mismatch":
                result.loc[result.index[0], "FDR"] = 0.9
            elif malformed == "spatial_fdr_mismatch":
                result.loc[result.index[0], "SpatialFDR"] = 0.9
            elif malformed == "da_corrupts_counts":
                mdata["milo"].X[0, 0] += 1
            if malformed == "da_returns_value":
                return mdata

        def annotate_nhoods(self, mdata, anno_col, feature_key="rna"):
            calls.append(("annotate_nhoods", anno_col, feature_key))
            rna = mdata[feature_key]
            milo = mdata["milo"]
            labels = science.pd.get_dummies(rna.obs[anno_col])
            fractions = science.sparse.csr_matrix(rna.obsm["nhoods"]).T.dot(labels.to_numpy(dtype=float))
            fractions = science.np.asarray(fractions)
            fractions = fractions / fractions.sum(axis=1, keepdims=True)
            milo.varm["frac_annotation"] = fractions
            milo.uns["annotation_labels"] = labels.columns.astype(str).tolist()
            milo.uns["annotation_obs"] = anno_col
            milo.var["nhood_annotation"] = labels.columns[fractions.argmax(axis=1)].astype(str)
            milo.var["nhood_annotation_frac"] = fractions.max(axis=1)
            if malformed == "annotation_mismatch":
                milo.var.loc[milo.var.index[0], "nhood_annotation"] = "wrong"
            elif malformed == "annotation_fraction_range":
                milo.var.loc[milo.var.index[0], "nhood_annotation_frac"] = 1.1
            elif malformed == "annotation_fraction_missing":
                del milo.varm["frac_annotation"]
            if malformed == "annotate_returns_value":
                return mdata

    pertpy.tl = types.SimpleNamespace(Milo=Milo)
    pertpy.calls = calls
    monkeypatch.setitem(sys.modules, "pertpy", pertpy)
    _install_fake_r_versions(monkeypatch)

    def neighbors(
        adata,
        n_neighbors=15,
        n_pcs=None,
        *,
        distances=None,
        use_rep=None,
        knn=True,
        method="umap",
        transformer=None,
        metric=None,
        metric_kwds=None,
        random_state=0,
        key_added=None,
        copy=False,
    ):
        calls.append(
            (
                "neighbors",
                adata.obs_names.tolist(),
                n_neighbors,
                n_pcs,
                distances,
                use_rep,
                knn,
                method,
                transformer,
                metric,
                metric_kwds,
                random_state,
                key_added,
                copy,
            )
        )
        if malformed == "neighbors_raises":
            random.random()
            science.np.random.random()
            raise RuntimeError("synthetic neighbors failure")
        n_obs = adata.n_obs
        rows = []
        cols = []
        values = []
        for row in range(n_obs):
            for column in (row % 4, (row + 1) % 4):
                rows.append(row)
                cols.append(column)
                values.append(0.5 + 0.1 * row if row < 4 else 1.0)
        distance_matrix = science.sparse.csr_matrix((values, (rows, cols)), shape=(n_obs, n_obs))
        connectivity = distance_matrix.copy()
        connectivity.data = science.np.ones_like(connectivity.data)
        if malformed == "membership_duplicate":
            connectivity_dense = connectivity.toarray()
            connectivity_dense[:, 1] = connectivity_dense[:, 0]
            connectivity = science.sparse.csr_matrix(connectivity_dense)
        adata.obsp[f"{key_added}_distances"] = distance_matrix
        adata.obsp[f"{key_added}_connectivities"] = connectivity
        adata.uns[key_added] = {
            "connectivities_key": f"{key_added}_connectivities",
            "distances_key": f"{key_added}_distances",
            "params": {
                "n_neighbors": n_neighbors,
                "method": method,
                "random_state": random_state,
                "metric": metric,
                "use_rep": use_rep,
            },
        }
        if malformed == "graph_metadata":
            adata.uns[key_added]["params"]["use_rep"] = "wrong_representation"
        elif malformed == "graph_method_metadata":
            adata.uns[key_added]["params"]["method"] = "gauss"
        elif malformed == "graph_seed_metadata":
            adata.uns[key_added]["params"]["random_state"] = random_state + 1
        elif malformed == "graph_metric_metadata":
            adata.uns[key_added]["params"]["metric"] = "cosine"
        elif malformed == "graph_keys":
            adata.uns[key_added]["distances_key"] = "wrong_distances"
        return adata.copy() if copy else None

    monkeypatch.setattr(science.sc.pp, "neighbors", neighbors)
    return pertpy


def _execute_milo(adata, **overrides):
    parameters = {
        "sample_key": "sample",
        "condition_key": "condition",
        "reference_condition": "control",
        "comparison_condition": "treated",
        "technical_batch_key": "batch",
        "categorical_covariate_keys_json": "[]",
        "continuous_covariate_keys_json": "[]",
        "annotation_key": "cell_type",
        "annotation_status": "provisional",
        "representation_key": "X_pca",
        "n_neighbors": 3,
        "neighborhood_proportion": 0.2,
        "mixed_annotation_threshold": 0.75,
        "spatial_fdr_threshold": 0.1,
        "min_abs_log2_fold_change": 0.0,
        "random_seed": 19,
    }
    parameters.update(overrides)
    return OpenBioSingleCellMiloDifferentialAbundance.execute(adata, **parameters).result


def test_milo_schema_exposes_one_pairwise_sample_contrast_and_report_outputs():
    schema = OpenBioSingleCellMiloDifferentialAbundance.GET_SCHEMA()

    assert [item.id for item in schema.inputs] == [
        "adata",
        "sample_key",
        "condition_key",
        "reference_condition",
        "comparison_condition",
        "technical_batch_key",
        "categorical_covariate_keys_json",
        "continuous_covariate_keys_json",
        "annotation_key",
        "annotation_status",
        "representation_key",
        "n_neighbors",
        "neighborhood_proportion",
        "mixed_annotation_threshold",
        "spatial_fdr_threshold",
        "min_abs_log2_fold_change",
        "random_seed",
    ]
    inputs = {item.id: item for item in schema.inputs}
    assert inputs["annotation_status"].default == "provisional"
    assert inputs["categorical_covariate_keys_json"].default == "[]"
    assert inputs["continuous_covariate_keys_json"].default == "[]"
    assert inputs["spatial_fdr_threshold"].default == 0.1
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_milo_runs_one_edge_r_sample_contrast_and_code_is_equivalent(science, monkeypatch):
    adata = _milo_adata(science)
    original = adata.copy()
    fake = _install_fake_milo_backend(science, monkeypatch)
    python_rng_before = random.getstate()
    numpy_rng_before = science.np.random.get_state()

    table_result, report, code = _execute_milo(adata)

    assert table_result.table.columns.tolist() == [
        "neighborhood_id",
        "index_cell",
        "neighborhood_size",
        "kth_distance",
        "reference_condition",
        "comparison_condition",
        "log2_fold_change",
        "log_counts_per_million",
        "quasi_likelihood_f",
        "p_value",
        "p_adjusted_bh",
        "spatial_fdr",
        "majority_annotation",
        "majority_annotation_fraction",
        "reported_annotation",
        "is_mixed",
    ]
    assert table_result.table["log2_fold_change"].tolist() == [-1.2, 0.8, 0.2, 1.5]
    assert table_result.table["majority_annotation"].tolist() == ["B", "T", "T", "B"]
    assert table_result.table["majority_annotation_fraction"].tolist() == [0.5, 1.0, 1.0, 0.5]
    assert table_result.table["reported_annotation"].tolist() == ["Mixed", "T", "T", "Mixed"]
    assert table_result.table["is_mixed"].tolist() == [True, False, False, True]
    assert report.summary["analysis_status"] == "sample_level_condition_inference"
    assert isinstance(report.summary["methods"], str)
    assert "Pairwise Milo" in report.summary["methods"]
    assert report.summary["comparison"]["positive_log2_fold_change"] == "treated enriched versus control"
    assert report.summary["design_evidence"]["samples_by_condition"] == {"control": 3, "treated": 3}
    assert report.summary["key_results"]["tested_neighborhoods"] == 4
    assert report.summary["key_results"]["comparison_enriched_calls"] == 2
    assert report.summary["key_results"]["reference_enriched_calls"] == 1
    assert "analysis_summary" not in report.summary
    assert report.summary["parameters"]["solver"] == "edger"
    assert report.summary["parameters"]["neighbor_transformer"] == "pynndescent"
    assert report.summary["software_versions"]["edgeR"] == "4.10.3"
    assert any("provisional" in warning.lower() for warning in report.summary["warnings"])
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<milo-code>", "exec")

    assert fake.calls[0][0] == "load"
    assert len(fake.calls[0][1]) == 24
    assert all("other" not in cell for cell in fake.calls[0][1])
    da_call = next(call for call in fake.calls if call[0] == "da_nhoods")
    assert da_call[1] == "~ openbio_condition + openbio_batch"
    assert da_call[2] == "openbio_conditioncomparison-openbio_conditionreference"
    assert da_call[3:] == (None, False, "rna", True, 50, 1e-5, "edger")
    make_call = next(call for call in fake.calls if call[0] == "make_nhoods")
    assert make_call[-2:] == (19, False)

    science.np.testing.assert_array_equal(adata.X, original.X)
    science.np.testing.assert_array_equal(adata.obsm["X_pca"], original.obsm["X_pca"])
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    assert random.getstate() == python_rng_before
    current_numpy_state = science.np.random.get_state()
    assert current_numpy_state[0] == numpy_rng_before[0]
    science.np.testing.assert_array_equal(current_numpy_state[1], numpy_rng_before[1])
    assert current_numpy_state[2:] == numpy_rng_before[2:]

    namespace = {}
    exec(code, namespace)
    reproduced_table, reproduced_summary = namespace["run_milo_differential_abundance"](adata)
    science.pd.testing.assert_frame_equal(reproduced_table, table_result.table)
    assert reproduced_summary == report.summary


@pytest.mark.parametrize(
    ("parameter", "value", "message"),
    [
        ("categorical_covariate_keys_json", "{}", "JSON array"),
        ("categorical_covariate_keys_json", '["batch", "batch"]', "duplicate"),
        ("categorical_covariate_keys_json", '[" batch"]', "surrounding whitespace"),
        ("categorical_covariate_keys_json", "[null]", "must be a string"),
        ("continuous_covariate_keys_json", "[NaN]", "non-standard JSON constant"),
        ("continuous_covariate_keys_json", '["condition"]', "distinct columns"),
    ],
)
def test_milo_rejects_ambiguous_nuisance_json_before_backend(science, parameter, value, message):
    adata = _milo_adata(science)

    with pytest.raises((TypeError, ValueError), match=message):
        _execute_milo(adata, **{parameter: value})


def test_milo_rejects_pseudoreplication_and_confounded_design_before_backend(science):
    too_few = _milo_adata(science, samples_per_condition=2, include_other=False)
    with pytest.raises(ValueError, match="at least 3 independent Samples"):
        _execute_milo(too_few)

    ambiguous = _milo_adata(science, include_other=False)
    ambiguous.obs.loc[ambiguous.obs.index[-1], "sample"] = "sample_1"
    with pytest.raises(ValueError, match="maps to multiple Conditions"):
        _execute_milo(ambiguous)

    confounded = _milo_adata(science, include_other=False)
    confounded.obs["batch"] = confounded.obs["condition"].map({"control": "batch_a", "treated": "batch_b"})
    with pytest.raises(ValueError, match="rank deficient"):
        _execute_milo(confounded)


def test_milo_encodes_valid_sample_level_nuisance_covariates_deterministically(science, monkeypatch):
    adata = _milo_adata(science)
    ages = {
        "sample_1": 30.0,
        "sample_2": 50.0,
        "sample_3": 40.0,
        "sample_4": 35.0,
        "sample_5": 55.0,
        "sample_6": 45.0,
        "sample_other": 60.0,
    }
    adata.obs["age"] = adata.obs["sample"].map(ages)
    adata.obs["batch"] = science.pd.Categorical(
        adata.obs["batch"],
        categories=["batch_b", "batch_a", "unused_batch"],
        ordered=True,
    )
    fake = _install_fake_milo_backend(science, monkeypatch)

    _, report, _ = _execute_milo(
        adata,
        technical_batch_key="",
        categorical_covariate_keys_json='["batch"]',
        continuous_covariate_keys_json='["age"]',
    )

    da_call = next(call for call in fake.calls if call[0] == "da_nhoods")
    assert da_call[1] == "~ openbio_condition + openbio_cat_0 + openbio_cont_0"
    assert da_call[2] == "openbio_conditioncomparison-openbio_conditionreference"
    assert report.summary["design_evidence"]["columns"] == [
        "openbio_conditionreference",
        "openbio_conditioncomparison",
        "openbio_cat_0level_1",
        "openbio_cont_0",
    ]
    assert report.summary["design_evidence"]["categorical_covariates"] == [
        {
            "source_key": "batch",
            "reference_level": "batch_b",
            "levels": ["batch_b", "batch_a"],
        }
    ]
    assert report.summary["design_evidence"]["continuous_covariates"] == [
        {"source_key": "age", "min": 30.0, "max": 55.0}
    ]
    assert report.summary["design_evidence"]["rank"] == 4
    assert report.summary["design_evidence"]["residual_degrees_of_freedom"] == 2


def test_milo_uses_collision_free_private_columns_and_graph_without_mutating_input(science, monkeypatch):
    adata = _milo_adata(science)
    adata.obs["openbio_sample"] = "preserve_sample"
    adata.obs["openbio_condition"] = "preserve_condition"
    adata.obs["openbio_batch"] = "preserve_batch"
    adata.uns["_openbio_milo_graph"] = {"sentinel": "preserve"}
    adata.obsp["_openbio_milo_graph_distances"] = science.sparse.eye(adata.n_obs, format="csr")
    adata.obsp["_openbio_milo_graph_connectivities"] = science.sparse.eye(adata.n_obs, format="csr")
    original = adata.copy()
    fake = _install_fake_milo_backend(science, monkeypatch)

    _execute_milo(adata)

    neighbor_call = next(call for call in fake.calls if call[0] == "neighbors")
    assert neighbor_call[-2] == "_openbio_milo_graph_1"
    da_call = next(call for call in fake.calls if call[0] == "da_nhoods")
    assert da_call[1] == "~ openbio_condition_1 + openbio_batch_1"
    assert da_call[2] == "openbio_condition_1comparison-openbio_condition_1reference"
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    assert adata.uns == original.uns
    for key in original.obsp:
        assert (adata.obsp[key] != original.obsp[key]).nnz == 0


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("same_conditions", "must differ"),
        ("absent_condition", "requested Conditions are absent"),
        ("role_collision", "distinct columns"),
        ("numeric_sample", "must contain strings"),
        ("missing_annotation", "contains missing values"),
        ("numeric_annotation", "must contain strings"),
        ("whitespace_obs_id", "surrounding whitespace"),
        ("whitespace_condition", "empty or non-canonical"),
    ],
)
def test_milo_rejects_ambiguous_roles_and_noncanonical_identity(science, case, message):
    adata = _milo_adata(science)
    overrides = {}
    if case == "same_conditions":
        overrides["comparison_condition"] = "control"
    elif case == "absent_condition":
        overrides["comparison_condition"] = "absent"
    elif case == "role_collision":
        overrides["sample_key"] = "condition"
    elif case == "numeric_sample":
        adata.obs["sample"] = range(adata.n_obs)
    elif case == "missing_annotation":
        adata.obs.loc[adata.obs.index[0], "cell_type"] = None
    elif case == "numeric_annotation":
        adata.obs["cell_type"] = range(adata.n_obs)
    elif case == "whitespace_obs_id":
        names = adata.obs_names.tolist()
        names[0] = f" {names[0]}"
        adata.obs_names = names
    elif case == "whitespace_condition":
        adata.obs.loc[adata.obs.index[0], "condition"] = " control"

    with pytest.raises((TypeError, ValueError), match=message):
        _execute_milo(adata, **overrides)


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("missing", "representation not found"),
        ("nonnumeric", "must be numeric"),
        ("nonfinite", "non-finite"),
        ("one_dimension", "at least two dimensions"),
        ("degenerate", "degenerate"),
        ("neighbors_too_large", "n_neighbors < selected cells"),
        ("zero_proportion", "greater than zero"),
        ("bad_mixed_threshold", "at most 1.0"),
        ("negative_lfc_threshold", "at least 0.0"),
        ("boolean_seed", "must be an integer"),
    ],
)
def test_milo_rejects_invalid_representation_and_numeric_controls(science, case, message):
    adata = _milo_adata(science)
    overrides = {}
    if case == "missing":
        overrides["representation_key"] = "missing"
    elif case == "nonnumeric":
        adata.obsm["X_pca"] = adata.obsm["X_pca"].astype(str)
    elif case == "nonfinite":
        adata.obsm["X_pca"][0, 0] = science.np.nan
    elif case == "one_dimension":
        adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :1]
    elif case == "degenerate":
        adata.obsm["X_pca"] = science.np.ones((adata.n_obs, 2))
    elif case == "neighbors_too_large":
        overrides["n_neighbors"] = 24
    elif case == "zero_proportion":
        overrides["neighborhood_proportion"] = 0.0
    elif case == "bad_mixed_threshold":
        overrides["mixed_annotation_threshold"] = 1.1
    elif case == "negative_lfc_threshold":
        overrides["min_abs_log2_fold_change"] = -0.1
    elif case == "boolean_seed":
        overrides["random_seed"] = True

    with pytest.raises((TypeError, ValueError), match=message):
        _execute_milo(adata, **overrides)


def test_milo_reporting_thresholds_do_not_filter_or_refit_the_neighborhood_universe(science, monkeypatch):
    adata = _milo_adata(science)
    fake = _install_fake_milo_backend(science, monkeypatch)

    baseline_table, _, _ = _execute_milo(adata)
    strict_table, strict_report, _ = _execute_milo(
        adata,
        spatial_fdr_threshold=0.0,
        min_abs_log2_fold_change=2.0,
    )

    science.pd.testing.assert_frame_equal(strict_table.table, baseline_table.table)
    assert strict_report.summary["key_results"]["tested_neighborhoods"] == 4
    assert strict_report.summary["key_results"]["comparison_enriched_calls"] == 0
    assert strict_report.summary["key_results"]["reference_enriched_calls"] == 0
    assert strict_report.summary["key_results"]["comparison_enriched"] == []
    assert strict_report.summary["key_results"]["reference_enriched"] == []
    assert len([call for call in fake.calls if call[0] == "da_nhoods"]) == 2


def test_milo_discloses_extreme_cell_count_imbalance_without_treating_cells_as_replicates(science, monkeypatch):
    adata = _milo_adata(science, include_other=False)
    keep = (adata.obs["condition"] == "treated") | (
        (adata.obs["condition"] == "control") & ~adata.obs["sample"].duplicated()
    )
    adata = adata[keep].copy()
    _install_fake_milo_backend(science, monkeypatch)

    _, report, _ = _execute_milo(adata)

    assert report.summary["design_evidence"]["samples_by_condition"] == {"control": 3, "treated": 3}
    assert report.summary["input_evidence"]["condition_cell_counts"] == {"control": 3, "treated": 12}
    assert report.summary["input_evidence"]["condition_cell_imbalance_ratio"] == 4.0
    assert any("4-fold imbalanced" in warning for warning in report.summary["warnings"])


@pytest.mark.parametrize(
    ("malformed", "message"),
    [
        ("graph_metadata", "metadata does not match"),
        ("graph_method_metadata", "metadata does not match"),
        ("graph_seed_metadata", "metadata does not match"),
        ("graph_metric_metadata", "metadata does not match"),
        ("neighbors_raises", "synthetic neighbors failure"),
        ("make_raises", "synthetic make_nhoods failure"),
        ("membership_nonbinary", "must be binary"),
        ("membership_graph_mismatch", "does not match the private graph"),
        ("membership_duplicate", "unique neighborhood memberships"),
        ("kth_zero", "finite positive kth-neighbor"),
        ("kth_graph_mismatch", "do not match the private graph"),
        ("count_axis", "rows do not match"),
        ("count_names", "rows do not match"),
        ("count_mismatch", "do not equal"),
        ("count_noninteger", "must be integers"),
        ("count_negative", "contains negative values"),
        ("index_cell_duplicate", "changed neighborhood index-cell identity"),
        ("count_corrupts_membership", "changed neighborhood membership"),
        ("da_corrupts_counts", "changed Sample-by-neighborhood counts"),
        ("da_raises", "synthetic da_nhoods failure"),
        ("da_returns_value", "da_nhoods violated in-place semantics"),
        ("annotate_returns_value", "annotate_nhoods violated in-place semantics"),
        ("result_missing_column", "missing columns"),
        ("result_boolean", "must be numeric"),
        ("result_negative_f", "cannot be negative"),
        ("result_nonfinite", "complete and finite"),
        ("result_range", "must lie in"),
        ("annotation_mismatch", "majority annotation does not match"),
        ("annotation_fraction_range", "must be complete and lie in"),
        ("annotation_fraction_missing", "full annotation fraction matrix"),
    ],
)
def test_milo_fails_closed_on_malformed_external_stage_outputs(science, monkeypatch, malformed, message):
    adata = _milo_adata(science)
    original = adata.copy()
    _install_fake_milo_backend(science, monkeypatch, malformed=malformed)
    python_rng_before = random.getstate()
    numpy_rng_before = science.np.random.get_state()

    with pytest.raises((TypeError, ValueError, RuntimeError), match=message):
        _execute_milo(adata)

    science.np.testing.assert_array_equal(adata.X, original.X)
    science.np.testing.assert_array_equal(adata.obsm["X_pca"], original.obsm["X_pca"])
    science.pd.testing.assert_frame_equal(adata.obs, original.obs)
    assert random.getstate() == python_rng_before
    current_numpy_state = science.np.random.get_state()
    assert current_numpy_state[0] == numpy_rng_before[0]
    science.np.testing.assert_array_equal(current_numpy_state[1], numpy_rng_before[1])
    assert current_numpy_state[2:] == numpy_rng_before[2:]


@pytest.mark.parametrize("malformed", ["fdr_mismatch", "spatial_fdr_mismatch"])
def test_milo_independently_verifies_bh_and_spatial_fdr(science, monkeypatch, malformed):
    adata = _milo_adata(science)
    _install_fake_milo_backend(science, monkeypatch, malformed=malformed)

    with pytest.raises(RuntimeError, match="independent .*FDR verification"):
        _execute_milo(adata)


@pytest.mark.parametrize(
    ("malformed", "message"),
    [
        ("graph_keys", "bundle key metadata"),
        ("random_count", "candidate index count"),
    ],
)
def test_milo_verifies_graph_bundle_and_random_sampling_contract(science, monkeypatch, malformed, message):
    adata = _milo_adata(science)
    _install_fake_milo_backend(science, monkeypatch, malformed=malformed)

    with pytest.raises(RuntimeError, match=message):
        _execute_milo(adata)


def test_milo_fails_closed_on_pertpy_signature_or_version_drift(science, monkeypatch):
    adata = _milo_adata(science)
    fake = _install_fake_milo_backend(science, monkeypatch)
    compatible_milo = fake.tl.Milo

    class MissingSeedMilo(compatible_milo):
        def make_nhoods(
            self,
            data,
            *,
            neighbors_key=None,
            feature_key="rna",
            prop=0.1,
            copy=False,
        ):
            return super().make_nhoods(
                data,
                neighbors_key=neighbors_key,
                feature_key=feature_key,
                prop=prop,
                seed=0,
                copy=copy,
            )

    fake.tl.Milo = MissingSeedMilo
    with pytest.raises(RuntimeError, match="Milo.make_nhoods.*missing parameters.*seed"):
        _execute_milo(adata)

    fake.tl.Milo = compatible_milo
    real_distribution_version = importlib_metadata.version

    def incompatible_pertpy_version(distribution):
        if distribution == "pertpy":
            return "1.4.0"
        return real_distribution_version(distribution)

    monkeypatch.setattr(importlib_metadata, "version", incompatible_pertpy_version)
    with pytest.raises(RuntimeError, match="supports the audited Pertpy"):
        _execute_milo(adata)


def test_milo_fails_closed_on_scanpy_signature_drift(science, monkeypatch):
    adata = _milo_adata(science)
    _install_fake_milo_backend(science, monkeypatch)

    def missing_transformer_neighbors(
        adata,
        n_neighbors=15,
        n_pcs=None,
        *,
        distances=None,
        use_rep=None,
        knn=True,
        method="umap",
        metric=None,
        metric_kwds=None,
        random_state=0,
        key_added=None,
        copy=False,
    ):
        raise AssertionError("incompatible Scanpy function must not be called")

    monkeypatch.setattr(science.sc.pp, "neighbors", missing_transformer_neighbors)
    with pytest.raises(RuntimeError, match="scanpy.pp.neighbors.*missing parameters.*transformer"):
        _execute_milo(adata)


def test_milo_requires_verified_r_and_bioconductor_versions_for_the_report(science, monkeypatch):
    adata = _milo_adata(science)
    fake = _install_fake_milo_backend(science, monkeypatch)
    broken_robjects = types.ModuleType("rpy2.robjects")

    def unavailable_r(_expression):
        raise RuntimeError("R runtime unavailable")

    broken_robjects.r = unavailable_r
    monkeypatch.setitem(sys.modules, "rpy2.robjects", broken_robjects)

    with pytest.raises(RuntimeError, match="could not verify R, edgeR, limma and statmod versions"):
        _execute_milo(adata)
    assert fake.calls == []
