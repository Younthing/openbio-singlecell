from __future__ import annotations

import os
from importlib import metadata as importlib_metadata
from importlib.util import find_spec

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, concat

from openbio_singlecell import composition_modeling
from openbio_singlecell.composition_modeling import (
    SCCODA_TABLE_COLUMNS,
    TASCCODA_TABLE_COLUMNS,
    run_sccoda_differential_composition,
    run_tasccoda_differential_composition,
    sccoda_differential_composition_code,
    tasccoda_differential_composition_code,
)
from openbio_singlecell.node_types import SummaryResultType, TableResultType
from openbio_singlecell.nodes_abundance import (
    OpenBioSingleCellSccodaDifferentialComposition,
    OpenBioSingleCellTasccodaDifferentialComposition,
)


def test_composition_node_schemas_are_atomic_and_method_specific():
    sccoda = OpenBioSingleCellSccodaDifferentialComposition.define_schema()
    tasccoda = OpenBioSingleCellTasccodaDifferentialComposition.define_schema()

    assert [item.id for item in sccoda.inputs] == [
        "adata",
        "sample_key",
        "annotation_key",
        "annotation_status",
        "condition_key",
        "reference_condition",
        "comparison_condition",
        "adjustment_covariate_keys_json",
        "reference_cell_type",
        "estimated_fdr",
        "num_samples",
        "num_warmup",
        "random_seed",
    ]
    assert [item.id for item in tasccoda.inputs] == [
        "adata",
        "sample_key",
        "annotation_key",
        "annotation_status",
        "hierarchy_keys_json",
        "condition_key",
        "reference_condition",
        "comparison_condition",
        "adjustment_covariate_keys_json",
        "reference_cell_type",
        "aggregation_bias",
        "num_samples",
        "num_warmup",
        "random_seed",
    ]
    assert [item.default for item in sccoda.inputs[1:]] == [
        "sample",
        "cell_type",
        "provisional",
        "condition",
        "",
        "",
        "[]",
        "automatic",
        0.05,
        10000,
        1000,
        0,
    ]
    assert [item.default for item in tasccoda.inputs[1:]] == [
        "sample",
        "cell_type",
        "provisional",
        "[]",
        "condition",
        "",
        "",
        "[]",
        "automatic",
        0.0,
        10000,
        1000,
        0,
    ]
    assert [item.display_name for item in sccoda.outputs] == ["table", "summary", "code"]
    assert [item.display_name for item in tasccoda.outputs] == ["table", "summary", "code"]
    assert [item.io_type for item in sccoda.outputs] == [
        TableResultType.io_type,
        SummaryResultType.io_type,
        "STRING",
    ]
    assert [item.io_type for item in tasccoda.outputs] == [
        TableResultType.io_type,
        SummaryResultType.io_type,
        "STRING",
    ]
    assert sccoda.inputs[3].options == ["provisional", "curated"]
    assert tasccoda.inputs[3].options == ["provisional", "curated"]
    assert all(item.advanced for item in sccoda.inputs[-3:])
    assert all(item.advanced for item in tasccoda.inputs[-4:])


def _composition_adata() -> AnnData:
    rows: list[dict[str, object]] = []
    for sample, condition, batch in [
        ("A1", "A", "x"),
        ("A2", "A", "y"),
        ("B1", "B", "x"),
        ("B2", "B", "y"),
    ]:
        for cell_type, count in [("T", 4), ("B", 2 if condition == "A" else 5)]:
            rows.extend(
                {
                    "sample": sample,
                    "condition": condition,
                    "cell_type": cell_type,
                    "batch": batch,
                }
                for _ in range(count)
            )
    obs = pd.DataFrame(rows, index=[f"cell_{index}" for index in range(len(rows))])
    obs["batch"] = pd.Categorical(obs["batch"], categories=["x", "y"], ordered=True)
    return AnnData(np.zeros((len(obs), 2)), obs=obs, var=pd.DataFrame(index=["g1", "g2"]))


def test_sccoda_rejects_extra_condition_before_backend_work():
    adata = _composition_adata()
    extra = adata[:1].copy()
    extra.obs["sample"] = "C1"
    extra.obs["condition"] = "C"
    combined = concat([adata, extra], index_unique="-")

    class BackendMustNotRun:
        def fit_sccoda(self, context):  # pragma: no cover - a call is the failure
            raise AssertionError("backend must not run")

    with pytest.raises(ValueError, match="exactly the declared two Condition levels"):
        run_sccoda_differential_composition(
            combined,
            sample_key="sample",
            annotation_key="cell_type",
            annotation_status="curated",
            condition_key="condition",
            reference_condition="A",
            comparison_condition="B",
            adjustment_covariate_keys_json="[]",
            reference_cell_type="T",
            estimated_fdr=0.05,
            num_samples=100,
            num_warmup=50,
            random_seed=0,
            _backend=BackendMustNotRun(),
        )


class _FlatFakeBackend:
    def fit_sccoda(self, context):
        assert context["design"]["formula"] == "1 + __openbio_condition_comparison"
        assert context["design"]["focal_column"] == "__openbio_condition_comparison"
        assert context["count_matrix"].shape == (4, 2)
        assert context["resolved_reference"] == "T"
        draws = context["num_samples"]
        effects = np.zeros((draws, 2), dtype=float)
        effects[: int(draws * 0.96), 0] = 0.4
        intercepts = np.tile(np.array([0.1, -0.1]), (draws, 1))
        return {
            "cell_types": list(context["cell_types"]),
            "resolved_reference": context["resolved_reference"],
            "posterior_intercept_samples": intercepts,
            "posterior_focal_effect_samples": effects,
            "backend_count_matrix": context["count_matrix"].copy(),
            "backend_sample_ids": list(context["sample_internal_order"]),
            "upstream_contract": {
                "model_type": "classic",
                "select_type": "spikeslab",
                "covariate_names": list(context["design"]["model_columns"]),
                "reference_index": 1,
                "automatic_reference_absence_threshold": 0.05,
                "pseudocount": 0.5,
                "algorithm": "NUTS",
                "num_samples": context["num_samples"],
                "num_warmup": context["num_warmup"],
                "chain_count": 1,
            },
            "diagnostics": {
                "acceptance_rate": 0.8,
                "ess_bulk": [600.0, 550.0],
                "ess_tail": [500.0, 480.0],
                "mcse_mean": [0.01, 0.02],
                "potential_energy": [10.0, 10.5, 9.8],
                "num_steps": [7.0, 8.0, 7.0],
                "step_size": [0.1, 0.1, 0.1],
                "rhat_available": False,
                "divergences_available": False,
            },
            "backend_runtime": {
                "jax_enable_x64": True,
                "jax_devices": [{"platform": "cpu", "device_kind": "fake-cpu"}],
            },
            "software_versions": {
                "openbio-singlecell": context["openbio_version"],
                "pertpy": "1.3.0",
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "anndata": "0.13.2",
                "jax": "0.7.2",
                "jaxlib": "0.7.2",
                "numpyro": "0.19.0",
                "arviz": "0.22.0",
                "mudata": "0.3.2",
                "patsy": "1.0.1",
            },
        }


def _flat_parameters() -> dict[str, object]:
    return {
        "sample_key": "sample",
        "annotation_key": "cell_type",
        "annotation_status": "curated",
        "condition_key": "condition",
        "reference_condition": "A",
        "comparison_condition": "B",
        "adjustment_covariate_keys_json": "[]",
        "reference_cell_type": "T",
        "estimated_fdr": 0.05,
        "num_samples": 100,
        "num_warmup": 50,
        "random_seed": 7,
        "openbio_version": "test",
    }


def _run_flat(adata: AnnData, backend: object, **overrides: object):
    parameters = _flat_parameters()
    parameters.update(overrides)
    return run_sccoda_differential_composition(adata, _backend=backend, **parameters)


def test_sccoda_fake_backend_has_focal_only_expected_fdr_table_and_code_parity():
    adata = _composition_adata()
    obs_before = adata.obs.copy(deep=True)
    uns_before = dict(adata.uns)

    table, summary = _run_flat(adata, _FlatFakeBackend())

    assert list(table.columns) == SCCODA_TABLE_COLUMNS
    assert table["cell_type"].tolist() == ["B", "T"]
    assert table["credible_effect"].tolist() == [True, False]
    assert table["inclusion_probability"].tolist() == [0.96, 0.0]
    assert table["inclusion_probability_threshold"].tolist() == [0.96, 0.96]
    assert table["realized_expected_fdr"].tolist() == pytest.approx([0.04, 0.04])
    assert table.loc[0, "model_coefficient"] == pytest.approx(0.4)
    assert table.loc[1, "reference_constraint"]
    assert summary["status"] == "limited_single_chain"
    assert summary["question"]["sample_is_experimental_unit"] is True
    assert summary["selection"]["strict_expected_fdr_inequality"] is True
    assert summary["diagnostics"]["rhat_available"] is False
    assert summary["diagnostics"]["divergences_available"] is False
    assert "converged" not in str(summary).lower()
    assert len(summary["references"]) >= 4
    assert summary["software_versions"]["pertpy"] == "1.3.0"
    pd.testing.assert_frame_equal(adata.obs, obs_before)
    assert adata.uns == uns_before

    code = sccoda_differential_composition_code(
        sample_key="sample",
        annotation_key="cell_type",
        annotation_status="curated",
        condition_key="condition",
        reference_condition="A",
        comparison_condition="B",
        adjustment_covariate_keys_json="[]",
        reference_cell_type="T",
        estimated_fdr=0.05,
        num_samples=100,
        num_warmup=50,
        random_seed=7,
        openbio_version="test",
    )
    compile(code, "<sccoda-code>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated_table, generated_summary = namespace["run_sccoda_differential_composition"](
        adata,
        _backend=_FlatFakeBackend(),
    )
    pd.testing.assert_frame_equal(generated_table, table)
    assert generated_summary == summary


def _hierarchy_adata() -> AnnData:
    rows: list[dict[str, object]] = []
    for sample, condition in [("A1", "A"), ("A2", "A"), ("B1", "B"), ("B2", "B")]:
        for cell_type, lineage, count in [
            ("T", "lymphoid", 4),
            ("B", "lymphoid", 2 if condition == "A" else 5),
            ("M", "myeloid", 3),
            ("N", "myeloid", 2),
        ]:
            rows.extend(
                {
                    "sample": sample,
                    "condition": condition,
                    "cell_type": cell_type,
                    "root": "immune",
                    "lineage": lineage,
                }
                for _ in range(count)
            )
    obs = pd.DataFrame(rows, index=[f"tree_cell_{index}" for index in range(len(rows))])
    return AnnData(np.zeros((len(obs), 2)), obs=obs, var=pd.DataFrame(index=["g1", "g2"]))


class _TreeFakeBackend:
    def fit_tasccoda(self, context):
        hierarchy = context["hierarchy"]
        assert hierarchy["levels_passed_to_pertpy"] == ["root", "lineage", "cell_type"]
        assert hierarchy["root"] == "root_immune"
        assert hierarchy["reference_path"] == ["root_immune", "lineage_lymphoid", "T"]
        assert context["aggregation_bias"] == -0.5
        draws = context["num_samples"]
        node_names = hierarchy["node_names"]
        node_samples = np.zeros((draws, len(node_names)), dtype=float)
        node_samples[:, node_names.index("B")] = 0.3
        node_samples[:, node_names.index("lineage_myeloid")] = 0.2
        node_samples[:, node_names.index("M")] = 0.1
        intercepts = np.tile(np.array([0.1, -0.2, 0.0, -0.1]), (draws, 1))
        reference_indices = [node_names.index(name) for name in hierarchy["reference_nodes"]]
        nonreference_indices = [index for index in range(len(node_names)) if index not in reference_indices]
        descendant_counts = np.asarray(hierarchy["ancestor_matrix"].sum(axis=0), dtype=float)
        phi = context["aggregation_bias"]
        scale = 1.0 / (1.0 + np.exp(-phi * (descendant_counts / len(context["cell_types"]) - 0.5)))
        lambda_1_scaled = 2.0 * 5.0 * scale
        return {
            "cell_types": list(context["cell_types"]),
            "resolved_reference": context["resolved_reference"],
            "posterior_intercept_samples": intercepts,
            "posterior_focal_node_effect_samples": node_samples,
            "posterior_theta_samples": np.full(draws, 0.5),
            "backend_count_matrix": context["count_matrix"].copy(),
            "backend_sample_ids": list(context["sample_internal_order"]),
            "upstream_contract": {
                "model_type": "tree_agg",
                "select_type": "sslasso",
                "covariate_names": list(context["design"]["model_columns"]),
                "reference_index": reference_indices,
                "reference_nodes": list(hierarchy["reference_nodes"]),
                "node_names": list(node_names),
                "ancestor_matrix": hierarchy["ancestor_matrix"].copy(),
                "automatic_reference_absence_threshold": 0.05,
                "pseudocount": 0.5,
                "algorithm": "NUTS",
                "num_samples": context["num_samples"],
                "num_warmup": context["num_warmup"],
                "chain_count": 1,
                "lambda_0": 50.0,
                "lambda_1": 5.0,
                "theta": 0.5,
                "phi": phi,
                "node_leaves_nonreference": descendant_counts[nonreference_indices],
                "lambda_1_scaled_nonreference": lambda_1_scaled[nonreference_indices],
            },
            "diagnostics": {
                "acceptance_rate": 0.81,
                "ess_bulk": [500.0, 510.0],
                "ess_tail": [450.0, 440.0],
                "mcse_mean": [0.01, 0.02],
                "potential_energy": [9.0, 9.5, 9.2],
                "num_steps": [7.0, 8.0, 7.0],
                "step_size": [0.1, 0.1, 0.1],
                "rhat_available": False,
                "divergences_available": False,
            },
            "backend_runtime": {
                "jax_enable_x64": True,
                "jax_devices": [{"platform": "cpu", "device_kind": "fake-cpu"}],
            },
            "software_versions": {
                "openbio-singlecell": context["openbio_version"],
                "pertpy": "1.3.0",
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "anndata": "0.13.2",
                "jax": "0.7.2",
                "jaxlib": "0.7.2",
                "numpyro": "0.19.0",
                "arviz": "0.22.0",
                "mudata": "0.3.2",
                "patsy": "1.0.1",
                "ete4": "4.4.0",
                "toytree": "3.0.10",
            },
        }


def _tree_parameters() -> dict[str, object]:
    return {
        "sample_key": "sample",
        "annotation_key": "cell_type",
        "annotation_status": "curated",
        "hierarchy_keys_json": '["root", "lineage"]',
        "condition_key": "condition",
        "reference_condition": "A",
        "comparison_condition": "B",
        "adjustment_covariate_keys_json": "[]",
        "reference_cell_type": "T",
        "aggregation_bias": -0.5,
        "num_samples": 100,
        "num_warmup": 50,
        "random_seed": 11,
        "openbio_version": "test",
    }


def _run_tree(adata: AnnData, backend: object, **overrides: object):
    parameters = _tree_parameters()
    parameters.update(overrides)
    return run_tasccoda_differential_composition(adata, _backend=backend, **parameters)


def test_tasccoda_fake_backend_separates_direct_nodes_and_derived_leaves_with_code_parity():
    adata = _hierarchy_adata()
    obs_before = adata.obs.copy(deep=True)

    table, summary = _run_tree(adata, _TreeFakeBackend())

    assert list(table.columns) == TASCCODA_TABLE_COLUMNS
    direct = table.loc[table["effect_scope"] == "hierarchy_node"].set_index("effect_name")
    leaves = table.loc[table["effect_scope"] == "derived_leaf"].set_index("effect_name")
    assert set(direct.index) == {"lineage_lymphoid", "lineage_myeloid", "B", "M", "N", "T"}
    assert direct.loc["B", "credible_effect"]
    assert direct.loc["lineage_myeloid", "credible_effect"]
    assert direct.loc["T", "reference_constraint"]
    assert leaves.loc["B", "model_effect"] == pytest.approx(0.3)
    assert leaves.loc["M", "model_effect"] == pytest.approx(0.3)
    assert leaves.loc["N", "model_effect"] == pytest.approx(0.2)
    assert leaves.loc["T", "model_effect"] == pytest.approx(0.0)
    assert leaves["credible_effect"].isna().all()
    assert summary["selection"]["lambda_0"] == 50.0
    assert summary["selection"]["lambda_1"] == 5.0
    assert summary["selection"]["theta"] == 0.5
    assert summary["hierarchy"]["reference_path"] == ["root_immune", "lineage_lymphoid", "T"]
    assert summary["hierarchy"]["derived_leaf_effects_are_propagated"] is True
    assert "estimated_fdr" not in summary["parameters"]
    assert "estimated_fdr" not in str(summary)
    pd.testing.assert_frame_equal(adata.obs, obs_before)

    code = tasccoda_differential_composition_code(
        sample_key="sample",
        annotation_key="cell_type",
        annotation_status="curated",
        hierarchy_keys_json='["root", "lineage"]',
        condition_key="condition",
        reference_condition="A",
        comparison_condition="B",
        adjustment_covariate_keys_json="[]",
        reference_cell_type="T",
        aggregation_bias=-0.5,
        num_samples=100,
        num_warmup=50,
        random_seed=11,
        openbio_version="test",
    )
    assert "estimated_fdr" not in code
    assert "est_fdr" not in code
    assert "set_fdr" not in code
    compile(code, "<tasccoda-code>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated_table, generated_summary = namespace["run_tasccoda_differential_composition"](
        adata,
        _backend=_TreeFakeBackend(),
    )
    pd.testing.assert_frame_equal(generated_table, table)
    assert generated_summary == summary


@pytest.mark.parametrize(
    ("value", "error", "match"),
    [
        ('{"batch": true}', TypeError, "decode to a JSON array"),
        ('["batch", "batch"]', ValueError, "duplicate column"),
        ('["condition"]', ValueError, "metadata roles must use distinct"),
        ("[NaN]", ValueError, "non-standard JSON constant"),
    ],
)
def test_adjustment_list_is_strict_json_and_role_disjoint(value, error, match):
    with pytest.raises(error, match=match):
        _run_flat(
            _composition_adata(),
            _FlatFakeBackend(),
            adjustment_covariate_keys_json=value,
        )


def test_sample_metadata_replication_and_rank_fail_before_backend_and_preserve_input():
    class BackendMustNotRun:
        def fit_sccoda(self, context):  # pragma: no cover - a call is the failure
            raise AssertionError("backend must not run")

    conflicting = _composition_adata()
    conflicting.obs.iloc[0, conflicting.obs.columns.get_loc("condition")] = "B"
    before = conflicting.obs.copy(deep=True)
    with pytest.raises(ValueError, match="unique within each biological Sample"):
        _run_flat(conflicting, BackendMustNotRun())
    pd.testing.assert_frame_equal(conflicting.obs, before)

    insufficient = _composition_adata()
    insufficient = insufficient[insufficient.obs["sample"] != "B2"].copy()
    with pytest.raises(ValueError, match="at least two biological Samples per Condition"):
        _run_flat(insufficient, BackendMustNotRun())

    confounded = _composition_adata()
    confounded.obs["condition_copy"] = (confounded.obs["condition"] == "B").astype(float)
    with pytest.raises(ValueError, match="rank deficient or perfectly confounded"):
        _run_flat(
            confounded,
            BackendMustNotRun(),
            adjustment_covariate_keys_json='["condition_copy"]',
        )


def test_paired_subject_like_adjustment_is_rejected_without_random_effect_semantics():
    rows = []
    for pair in range(3):
        for condition in ["A", "B"]:
            sample = f"{condition}{pair}"
            for cell_type in ["B", "T"]:
                rows.extend(
                    {
                        "sample": sample,
                        "condition": condition,
                        "cell_type": cell_type,
                        "subject": f"P{pair}",
                    }
                    for _ in range(2)
                )
    obs = pd.DataFrame(rows, index=[f"paired_{index}" for index in range(len(rows))])
    obs["subject"] = pd.Categorical(obs["subject"], categories=["P0", "P1", "P2"], ordered=True)
    adata = AnnData(np.zeros((len(obs), 1)), obs=obs, var=pd.DataFrame(index=["g1"]))

    with pytest.raises(ValueError, match="paired/repeated-Sample pattern"):
        _run_flat(
            adata,
            _FlatFakeBackend(),
            adjustment_covariate_keys_json='["subject"]',
        )


def test_expected_fdr_uses_strict_inequality_and_official_fallback():
    class BoundaryBackend(_FlatFakeBackend):
        def fit_sccoda(self, context):
            result = super().fit_sccoda(context)
            effects = np.zeros_like(result["posterior_focal_effect_samples"])
            effects[:95, 0] = 0.4
            result["posterior_focal_effect_samples"] = effects
            return result

    table, summary = _run_flat(_composition_adata(), BoundaryBackend())

    assert not table["credible_effect"].any()
    assert table["inclusion_probability_threshold"].tolist() == [1.0, 1.0]
    assert summary["selection"]["fallback_threshold_used"] is True
    assert summary["selection"]["candidate_expected_fdr_before_floor"] == 0.0
    assert summary["selection"]["realized_expected_fdr_after_floor"] == 0.0


@pytest.mark.parametrize("corruption", ["count_matrix", "posterior", "diagnostics"])
def test_backend_corruption_fails_atomically_and_generated_function_matches(corruption):
    class CorruptBackend(_FlatFakeBackend):
        def fit_sccoda(self, context):
            result = super().fit_sccoda(context)
            if corruption == "count_matrix":
                result["backend_count_matrix"][0, 0] += 1
            elif corruption == "posterior":
                result["posterior_focal_effect_samples"][0, 0] = np.nan
            else:
                result["diagnostics"]["rhat_available"] = True
            return result

    adata = _composition_adata()
    before = adata.obs.copy(deep=True)
    with pytest.raises(RuntimeError):
        _run_flat(adata, CorruptBackend())
    pd.testing.assert_frame_equal(adata.obs, before)

    code = sccoda_differential_composition_code(**_flat_parameters())
    namespace: dict[str, object] = {}
    exec(code, namespace)
    with pytest.raises(RuntimeError):
        namespace["run_sccoda_differential_composition"](adata, _backend=CorruptBackend())


@pytest.mark.parametrize("malformation", ["multiple_roots", "multiple_leaf_paths", "unused_category"])
def test_hierarchy_contract_rejects_malformed_declared_trees_before_backend(malformation):
    adata = _hierarchy_adata()
    if malformation == "multiple_roots":
        adata.obs.loc[adata.obs["cell_type"].isin(["M", "N"]), "root"] = "stromal"
        match = "one root"
    elif malformation == "multiple_leaf_paths":
        first_b = adata.obs.index[adata.obs["cell_type"] == "B"][0]
        adata.obs.loc[first_b, "lineage"] = "myeloid"
        match = "more than one hierarchy path"
    else:
        adata.obs["lineage"] = pd.Categorical(
            adata.obs["lineage"],
            categories=["lymphoid", "myeloid", "unused"],
            ordered=True,
        )
        match = "unused categories"
    before = adata.obs.copy(deep=True)

    with pytest.raises(ValueError, match=match):
        _run_tree(adata, _TreeFakeBackend())
    pd.testing.assert_frame_equal(adata.obs, before)


def test_tasccoda_rejects_backend_hierarchy_or_fixed_prior_drift():
    class CorruptTreeBackend(_TreeFakeBackend):
        def fit_tasccoda(self, context):
            result = super().fit_tasccoda(context)
            result["upstream_contract"]["lambda_1"] = 3.5
            return result

    with pytest.raises(RuntimeError, match="prior 'lambda_1'"):
        _run_tree(_hierarchy_adata(), CorruptTreeBackend())


def test_numeric_sample_identifiers_are_collision_safe_for_pertpy_loading():
    adata = _composition_adata()
    mapping = {"A1": 1, "A2": 2, "B1": 3, "B2": 4}
    adata.obs["sample"] = adata.obs["sample"].map(mapping)

    table, summary = _run_flat(adata, _FlatFakeBackend())

    assert len(table) == 2
    assert [row["sample"] for row in summary["input"]["complete_sample_cell_type_counts"][::2]] == [1, 2, 3, 4]


def test_public_nodes_wrap_table_summary_and_specialized_code_atomically(monkeypatch):
    monkeypatch.setattr(
        composition_modeling,
        "_PinnedPertpyBackend",
        lambda method, openbio_version: _FlatFakeBackend() if method == "sccoda" else _TreeFakeBackend(),
    )
    flat_parameters = _flat_parameters()
    flat_parameters.pop("openbio_version")
    flat_table, flat_report, flat_code = OpenBioSingleCellSccodaDifferentialComposition.execute(
        _composition_adata(),
        **flat_parameters,
    ).result
    tree_parameters = _tree_parameters()
    tree_parameters.pop("openbio_version")
    tree_table, tree_report, tree_code = OpenBioSingleCellTasccodaDifferentialComposition.execute(
        _hierarchy_adata(),
        **tree_parameters,
    ).result

    assert list(flat_table.table.columns) == SCCODA_TABLE_COLUMNS
    assert flat_report.summary["node_id"] == "OpenBioSingleCellSccodaDifferentialComposition"
    assert "def run_sccoda_differential_composition(adata" in flat_code
    assert list(tree_table.table.columns) == TASCCODA_TABLE_COLUMNS
    assert tree_report.summary["node_id"] == "OpenBioSingleCellTasccodaDifferentialComposition"
    assert "def run_tasccoda_differential_composition(adata" in tree_code


def test_real_adapter_rejects_any_nonexact_pertpy_version_before_backend_work(monkeypatch):
    real_version = importlib_metadata.version

    def version(distribution):
        return "1.2.9" if distribution == "pertpy" else real_version(distribution)

    monkeypatch.setattr(importlib_metadata, "version", version)
    adata = _composition_adata()
    before = adata.obs.copy(deep=True)

    with pytest.raises(RuntimeError, match="requires exact Pertpy 1.3.0.*1.2.9"):
        run_sccoda_differential_composition(adata, **_flat_parameters())
    pd.testing.assert_frame_equal(adata.obs, before)


_REAL_SCCODA_SMOKE_AVAILABLE = (
    os.environ.get("OPENBIO_RUN_PERTPY_COMPOSITION_SMOKE") == "1"
    and all(find_spec(module) is not None for module in ["pertpy", "jax", "jaxlib", "numpyro", "arviz", "mudata"])
    and importlib_metadata.version("pertpy") == "1.3.0"
)


@pytest.mark.skipif(
    not _REAL_SCCODA_SMOKE_AVAILABLE,
    reason="set OPENBIO_RUN_PERTPY_COMPOSITION_SMOKE=1 in an isolated Pertpy 1.3.0 CPU environment",
)
def test_optional_real_pinned_sccoda_cpu_smoke():
    table, summary = run_sccoda_differential_composition(
        _composition_adata(),
        **{
            **_flat_parameters(),
            "num_samples": 40,
            "num_warmup": 40,
            "random_seed": 0,
        },
    )

    assert list(table.columns) == SCCODA_TABLE_COLUMNS
    assert summary["software_versions"]["pertpy"] == "1.3.0"
    assert summary["diagnostics"]["chain_count"] == 1
