from __future__ import annotations

import builtins
import copy
import json
import random
import sys
import types

import pytest

from openbio_singlecell.operations_differential import run_scvi_differential_owned
from openbio_singlecell.scvi_de_evidence import (
    SCVI_DE_REFERENCES,
    SCVI_DE_SUMMARY_SCHEMA,
    analyze_scvi_model_de_evidence,
    canonicalize_scvi_de_table,
    resolve_scvi_de_mode,
    resolve_scvi_population_scope,
    scvi_de_code,
)
from openbio_singlecell.scvi_model import SCVI_MODEL_ARTIFACT_SCHEMA, SCVIModel


def current_model_parameters(**overrides):
    parameters = {
        "source": "X",
        "count_source_state": "counts",
        "count_source_state_evidence": "test fixture raw counts",
    }
    parameters.update(overrides)
    return parameters


class FullFakeSCVI:
    def __init__(self, registered_adata, science):
        self.adata = registered_adata
        self.is_trained = True
        self.science = science
        self.calls = []
        self.deregister_calls = []
        self.backend_error = None
        self.unrelated_manager = object()
        self.manager_store = {
            id(registered_adata): registered_adata,
            id(self.unrelated_manager): self.unrelated_manager,
        }
        self.class_manager_store = dict(self.manager_store)

    def differential_expression(
        self,
        adata=None,
        groupby=None,
        group1=None,
        group2=None,
        mode="vanilla",
        delta=0.25,
        batch_size=None,
        all_stats=True,
        batch_correction=False,
        batchid1=None,
        batchid2=None,
        fdr_target=0.05,
        silent=False,
        weights="uniform",
        filter_outlier_cells=False,
        importance_weighting_kwargs=None,
        dataloader=None,
        **kwargs,
    ):
        call = {
            "adata": adata,
            "groupby": groupby,
            "group1": group1,
            "group2": group2,
            "mode": mode,
            "delta": delta,
            "batch_size": batch_size,
            "all_stats": all_stats,
            "batch_correction": batch_correction,
            "batchid1": batchid1,
            "batchid2": batchid2,
            "fdr_target": fdr_target,
            "silent": silent,
            "weights": weights,
            "filter_outlier_cells": filter_outlier_cells,
            "importance_weighting_kwargs": importance_weighting_kwargs,
            "dataloader": dataloader,
            **kwargs,
        }
        self.calls.append(call)
        self.manager_store[id(adata)] = adata
        adata.obs["backend_mutation"] = True
        if self.backend_error is not None:
            raise self.backend_error
        group1 = group1[0]
        common = {
            "scale1": [1.0, 2.0, 3.0],
            "scale2": [2.0, 1.0, 2.5],
            "raw_mean1": [2.0, 3.0, 4.0],
            "raw_mean2": [3.0, 2.0, 3.5],
            "non_zeros_proportion1": [0.5, 0.75, 1.0],
            "non_zeros_proportion2": [0.75, 0.5, 1.0],
            "raw_normalized_mean1": [100.0, 200.0, 300.0],
            "raw_normalized_mean2": [200.0, 100.0, 250.0],
            "comparison": [f"{group1} vs {group2}"] * 3,
            "group1": [group1] * 3,
            "group2": [group2] * 3,
        }
        if mode == "vanilla":
            return self.science.pd.DataFrame(
                {
                    "proba_m1": [0.2, 0.95, 0.6],
                    "proba_m2": [0.8, 0.05, 0.4],
                    "bayes_factor": [-1.3, 2.9, 0.4],
                    **common,
                },
                index=["g2", "g1", "g3"],
            )
        target = fdr_target
        return self.science.pd.DataFrame(
            {
                "proba_de": [0.8, 0.95, 0.6],
                "proba_not_de": [0.2, 0.05, 0.4],
                "bayes_factor": [1.3, 2.9, 0.4],
                "pseudocounts": [0.5, 0.5, 0.5],
                "delta": [delta] * 3,
                "lfc_mean": [-1.0, 1.2, 0.2],
                "lfc_median": [-0.9, 1.1, 0.1],
                "lfc_std": [0.2, 0.3, 0.1],
                "lfc_min": [-1.4, 0.5, -0.1],
                "lfc_max": [-0.5, 1.8, 0.4],
                f"is_de_fdr_{target}": [False, True, False],
                **common,
            },
            index=["g2", "g1", "g3"],
        )

    def deregister_manager(self, analysis_adata):
        self.deregister_calls.append(analysis_adata)
        if id(analysis_adata) not in self.class_manager_store:
            raise ValueError("Analysis manager is not registered in the class store.")
        self.manager_store.pop(id(analysis_adata), None)
        self.class_manager_store.pop(id(analysis_adata), None)

    def get_anndata_manager(self, analysis_adata, required=False):
        manager = self.manager_store.get(id(analysis_adata))
        if manager is None and required:
            raise ValueError("Analysis manager is not registered for this model instance.")
        return manager

    def register_manager(self, manager):
        self.class_manager_store[id(manager)] = manager


@pytest.fixture
def scvi_evidence_fixture(science, monkeypatch):
    class FakeDifferentialComputation:
        def get_bayes_factors(
            self,
            *,
            n_samples_overall=5000,
            use_permutation=False,
            pseudocounts=None,
            change_fn=None,
            m1_domain_fn=None,
            test_mode="three",
        ):
            raise AssertionError("Only the installed signature is inspected.")

    class FakeInferenceMode:
        def __init__(self, enabled):
            self.enabled = enabled

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

    torch_state = {"value": b"initial"}
    fake_torch = types.ModuleType("torch")
    fake_torch.random = types.SimpleNamespace(
        get_rng_state=lambda: torch_state["value"],
        set_rng_state=lambda value: torch_state.__setitem__("value", value),
    )
    fake_torch.cuda = types.SimpleNamespace(is_available=lambda: False)
    fake_torch.version = types.SimpleNamespace(cuda=None)
    fake_torch.backends = types.SimpleNamespace(
        cudnn=types.SimpleNamespace(version=lambda: None),
        mps=types.SimpleNamespace(is_available=lambda: False),
    )
    fake_torch.inference_mode = FakeInferenceMode
    fake_scvi = types.ModuleType("scvi")
    fake_scvi.__path__ = []
    fake_scvi.__version__ = "1.5.0.post1"
    fake_scvi.settings = types.SimpleNamespace(seed=None)
    fake_model = types.ModuleType("scvi.model")
    fake_model.__path__ = []
    fake_base = types.ModuleType("scvi.model.base")
    fake_base.DifferentialComputation = FakeDifferentialComputation
    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    monkeypatch.setitem(sys.modules, "scvi.model", fake_model)
    monkeypatch.setitem(sys.modules, "scvi.model.base", fake_base)

    obs = science.pd.DataFrame(
        {
            "batch": science.pd.Categorical(
                ["b1", "b1", "b2", "b3", "b1", "b2", "b2", "b4"],
                categories=["b1", "b2", "b3", "b4", "unused"],
            ),
            "group": science.pd.Categorical(
                ["A", "A", "A", "A", "B", "B", "B", "B"],
                categories=["A", "B", "unused-group"],
            ),
            "population": science.pd.Categorical(["keep"] * 6 + ["drop"] * 2),
        },
        index=[f"c{index}" for index in range(8)],
    )
    var = science.pd.DataFrame(index=["g1", "g2", "g3"])
    registered = science.ad.AnnData(
        science.np.arange(1, 25, dtype=float).reshape(8, 3),
        obs=obs,
        var=var,
    )
    raw_model = FullFakeSCVI(registered, science)
    model = SCVIModel(
        raw_model,
        registered,
        {
            "source": "X",
            "technical_batch_key": "batch",
            "categorical_covariates": [],
            "continuous_covariates": [],
            "size_factor_key": "",
            "random_seed": 17,
            "count_source_state": "counts",
            "count_source_state_evidence": "test fixture raw counts",
        },
        diagnostics={"actual_epochs": 3, "metrics": {"elbo_train": {"last": 5.0}}},
    )
    return registered, raw_model, model


def _run(adata, model, **overrides):
    parameters = {
        "groupby": "group",
        "group1": "A",
        "group2": "B",
        "subset_column": None,
        "subset_value": None,
        "mode": "change",
        "delta": 0.25,
        "fdr_target": 0.05,
        "batch_handling": "shared_technical_batches",
        "n_samples_overall": 5000,
        "random_seed": 23,
    }
    parameters.update(overrides)
    return analyze_scvi_model_de_evidence(adata, model, **parameters)


def test_scvi_dynamic_resolvers_reject_inactive_parameters():
    assert resolve_scvi_population_scope(None) == (None, None)
    assert resolve_scvi_population_scope({"population_scope": "obs_value", "subset_column": "x", "subset_value": "y"}) == (
        "x",
        "y",
    )
    assert resolve_scvi_de_mode(None) == ("change", 0.25, 0.05)
    assert resolve_scvi_de_mode({"mode": "vanilla"}) == ("vanilla", None, None)
    with pytest.raises(ValueError, match="inactive"):
        resolve_scvi_population_scope({"population_scope": "all", "subset_column": "x"})
    with pytest.raises(ValueError, match="inactive"):
        resolve_scvi_de_mode({"mode": "vanilla", "delta": 0.25})


def test_scvi_model_de_uses_only_shared_supported_batches_and_reports(scvi_evidence_fixture, science):
    adata, raw_model, model = scvi_evidence_fixture
    obs_before = adata.obs.copy(deep=True)
    python_state = random.getstate()
    numpy_state = science.np.random.get_state()

    table, summary = _run(adata, model)

    assert len(raw_model.calls) == 1
    call = raw_model.calls[0]
    assert call["group1"] == ["A"]
    assert call["group2"] == "B"
    assert call["batch_correction"] is True
    assert call["batchid1"] == ["b1", "b2"]
    assert call["batchid2"] == ["b1", "b2"]
    assert call["test_mode"] == "two"
    assert call["all_stats"] is True
    assert call["weights"] == "uniform"
    assert call["filter_outlier_cells"] is False
    assert call["use_permutation"] is False
    assert call["pseudocounts"] is None
    assert len(raw_model.deregister_calls) == 1
    assert raw_model.deregister_calls[0] is call["adata"]
    assert raw_model.manager_store == {
        id(adata): adata,
        id(raw_model.unrelated_manager): raw_model.unrelated_manager,
    }
    assert raw_model.class_manager_store == raw_model.manager_store
    assert "backend_mutation" not in adata.obs
    science.pd.testing.assert_frame_equal(adata.obs, obs_before)
    assert random.getstate() == python_state
    current_numpy_state = science.np.random.get_state()
    assert current_numpy_state[0] == numpy_state[0]
    science.np.testing.assert_array_equal(current_numpy_state[1], numpy_state[1])
    assert current_numpy_state[2:] == numpy_state[2:]

    assert table["gene"].tolist() == ["g1", "g2", "g3"]
    assert table.columns[:4].tolist() == ["gene", "comparison", "group1", "group2"]
    assert "is_de_fdr_0.05" not in table.columns
    assert table["is_de_fdr"].tolist() == [True, False, False]
    assert summary["schema_version"] == SCVI_DE_SUMMARY_SCHEMA
    assert summary["status"] == "exploratory_model_evidence"
    assert summary["comparison"]["shared_technical_batches_used"] == ["b1", "b2"]
    assert summary["comparison"]["excluded_nonshared_technical_batches"] == ["b3", "b4"]
    assert summary["comparison"]["unused_declared_group_categories"] == ["unused-group"]
    assert isinstance(summary["methods"], str)
    assert "scVI" in summary["methods"]
    assert isinstance(summary["results"], str)
    assert summary["key_results"]["posterior_fdr_tagged_features"] == 1
    assert summary["key_results"]["posterior_fdr_tagged_positive_features"] == 1
    assert summary["key_results"]["posterior_fdr_tagged_negative_features"] == 0
    assert summary["key_results"]["leading_tagged_positive_effects"] == [
        {
            "gene": "g1",
            "lfc_mean": 1.2,
            "proba_de": 0.95,
            "bayes_factor": 2.9,
            "posterior_fdr_tagged": True,
            "raw_mean_group1": 3.0,
            "raw_mean_group2": 2.0,
        }
    ]
    assert summary["key_results"]["leading_tagged_negative_effects"] == []
    assert summary["key_results"]["top_model_evidence"] == []
    assert summary["key_results"]["top_ranked_effects"][0]["bayes_factor"] == 2.9
    assert summary["model_evidence"]["artifact_schema_version"] == SCVI_MODEL_ARTIFACT_SCHEMA
    assert summary["model_evidence"]["registered_count_state"] == "counts"
    assert summary["model_evidence"]["registered_count_state_evidence"] == "test fixture raw counts"
    assert summary["model_evidence"]["accelerator_runtime"] == {
        "training_device": "unknown",
        "cuda_available": False,
        "torch_cuda_compiled_version": None,
        "cudnn_version": None,
        "mps_available": False,
    }
    assert summary["runtime_environment"]["accelerator"] == summary["model_evidence"]["accelerator_runtime"]
    assert summary["method_details"]["test_mode"] == "two"
    assert "method" not in summary
    assert "analysis_summary" not in summary
    assert any(reference["doi"] == "10.1038/s41587-021-01206-w" for reference in summary["references"])
    assert "not independent biological Samples" in summary["warnings"][0]
    json.dumps(summary, allow_nan=False)


def test_scvi_model_de_generated_source_is_exact(scvi_evidence_fixture, science):
    adata, _, model = scvi_evidence_fixture
    expected_table, expected_summary = _run(adata, model)
    source = scvi_de_code(
        groupby="group",
        group1="A",
        group2="B",
        subset_column=None,
        subset_value=None,
        mode="change",
        delta=0.25,
        fdr_target=0.05,
        batch_handling="shared_technical_batches",
        n_samples_overall=5000,
        random_seed=23,
    )
    assert "from openbio_singlecell" not in source
    assert "import openbio_singlecell" not in source
    original_import = builtins.__import__

    def audited_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.startswith("openbio_singlecell"):
            raise AssertionError(f"Generated code imported internal OpenBio implementation: {name}")
        return original_import(name, globals, locals, fromlist, level)

    generated_builtins = dict(vars(builtins))
    generated_builtins["__import__"] = audited_import
    namespace = {"__builtins__": generated_builtins}
    exec(source, namespace)
    actual_table, actual_summary = namespace["scvi_model_de_evidence"](adata, model)
    science.pd.testing.assert_frame_equal(actual_table, expected_table)
    assert actual_summary == expected_summary


def test_scvi_model_de_node_returns_table_summary_and_equivalent_code(scvi_evidence_fixture, science):
    adata, _, model = scvi_evidence_fixture
    result, report, code = run_scvi_differential_owned(
        adata,
        model,
        groupby="group",
        group1="A",
        group2="B",
        population_scope={"population_scope": "all"},
        mode={"mode": "change", "delta": 0.25, "fdr_target": 0.05},
        batch_handling="shared_technical_batches",
        n_samples_overall=5000,
        random_seed=23,
    )
    assert result.source["operation"] == "scvi_model_de_evidence"
    assert report.summary["status"] == "exploratory_model_evidence"
    assert result.table.shape == (3, 23)
    namespace = {}
    exec(code, namespace)
    reproduced_table, reproduced_summary = namespace["scvi_model_de_evidence"](adata, model)
    science.pd.testing.assert_frame_equal(reproduced_table, result.table)
    assert reproduced_summary == report.summary


def test_scvi_model_de_observed_batches_omits_batch_ids(scvi_evidence_fixture):
    adata, raw_model, model = scvi_evidence_fixture
    _, summary = _run(adata, model, batch_handling="observed_technical_batches")
    call = raw_model.calls[-1]
    assert call["batch_correction"] is False
    assert call["batchid1"] is None
    assert call["batchid2"] is None
    assert summary["comparison"]["shared_technical_batches_used"] == ["b1", "b2"]
    assert any("Observed Technical-batch composition" in warning for warning in summary["warnings"])


def test_scvi_model_de_rejects_complete_batch_confounding_before_backend(scvi_evidence_fixture, science):
    adata, _, _ = scvi_evidence_fixture
    adata = adata.copy()
    adata.obs["batch"] = science.pd.Categorical(["b1"] * 4 + ["b2"] * 4)
    raw_model = FullFakeSCVI(adata, science)
    model = SCVIModel(
        raw_model,
        adata,
        current_model_parameters(
            technical_batch_key="batch",
            categorical_covariates=[],
            continuous_covariates=[],
            size_factor_key="",
        ),
    )
    with pytest.raises(ValueError, match="share no observed Technical-batch"):
        _run(adata, model)
    assert raw_model.calls == []
    assert raw_model.deregister_calls == []


def test_scvi_model_de_uses_model_owned_batch_metadata(scvi_evidence_fixture):
    adata, raw_model, model = scvi_evidence_fixture
    downstream = adata.copy()
    downstream.obs["batch"] = ["tampered"] * downstream.n_obs
    _run(downstream, model)
    assert raw_model.calls[-1]["batchid1"] == ["b1", "b2"]


def test_scvi_model_detects_registered_count_tampering(scvi_evidence_fixture):
    registered, raw_model, model = scvi_evidence_fixture
    registered.X[0, 0] += 1
    with pytest.raises(RuntimeError, match="count matrix changed"):
        _run(registered, model)
    assert raw_model.calls == []


def test_scvi_model_fingerprints_the_declared_count_layer(scvi_evidence_fixture, science):
    registered, _, _ = scvi_evidence_fixture
    registered = registered.copy()
    registered.layers["counts"] = registered.X.copy()
    raw_model = FullFakeSCVI(registered, science)
    model = SCVIModel(
        raw_model,
        registered,
        current_model_parameters(
            source="layer",
            counts_layer="counts",
            technical_batch_key="batch",
            categorical_covariates=[],
            continuous_covariates=[],
            size_factor_key="",
        ),
    )
    registered.X[0, 0] += 100
    assert model.evidence["registered_count_source"] == "layers['counts']"
    registered.layers["counts"][0, 0] += 1
    with pytest.raises(RuntimeError, match="count matrix changed"):
        _ = model.evidence


def test_scvi_model_detects_registered_setup_metadata_tampering(scvi_evidence_fixture):
    registered, _, model = scvi_evidence_fixture
    registered.obs.loc["c0", "batch"] = "b2"
    with pytest.raises(RuntimeError, match="setup metadata changed"):
        _ = model.evidence


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda table: table.drop(index="g3"), "exactly one row per fitted feature"),
        (lambda table: table.assign(proba_de=[1.1, 0.95, 0.6]), "must lie in"),
        (lambda table: table.assign(proba_not_de=[0.3, 0.05, 0.4]), "complementary"),
        (lambda table: table.assign(delta=[0.5, 0.5, 0.5]), "delta values"),
        (lambda table: table.assign(extra=[1.0, 2.0, 3.0]), "unexpected"),
        (lambda table: table.rename(columns={"is_de_fdr_0.05": "is_de_fdr_0.1"}), "does not match"),
    ],
)
def test_scvi_canonicalizer_rejects_malformed_backend(scvi_evidence_fixture, mutation, message):
    adata, raw_model, model = scvi_evidence_fixture
    raw = raw_model.differential_expression(
        adata=adata,
        group1=["A"],
        group2="B",
        mode="change",
        delta=0.25,
        fdr_target=0.05,
    )
    with pytest.raises(RuntimeError, match=message):
        canonicalize_scvi_de_table(
            mutation(copy.deepcopy(raw)),
            fitted_features=model.var_names,
            group1="A",
            group2="B",
            mode="change",
            delta=0.25,
            fdr_target=0.05,
        )


def test_scvi_model_de_rejects_groupby_declared_as_nuisance(scvi_evidence_fixture):
    adata, raw_model, model = scvi_evidence_fixture
    with pytest.raises(ValueError, match="declared as a Technical/nuisance"):
        _run(adata, model, groupby="batch", group1="b1", group2="b2")
    assert raw_model.calls == []


def test_scvi_model_de_rejects_inactive_mode_parameters(scvi_evidence_fixture):
    adata, raw_model, model = scvi_evidence_fixture
    with pytest.raises(ValueError, match="vanilla mode does not accept"):
        _run(adata, model, mode="vanilla")
    assert raw_model.calls == []


def test_scvi_model_de_vanilla_discloses_directional_model_evidence(scvi_evidence_fixture):
    adata, raw_model, model = scvi_evidence_fixture

    table, summary = _run(adata, model, mode="vanilla", delta=None, fdr_target=None)

    assert table["gene"].tolist() == ["g1", "g3", "g2"]
    assert "is_de_fdr" not in table
    assert raw_model.calls[-1]["mode"] == "vanilla"
    assert summary["key_results"]["posterior_fdr_tagged_features"] is None
    assert summary["key_results"]["posterior_fdr_tagged_positive_features"] is None
    assert summary["key_results"]["leading_tagged_positive_effects"] == []
    assert summary["key_results"]["top_model_evidence"][0] == {
        "gene": "g1",
        "bayes_factor": 2.9,
        "probability_group1_higher": 0.95,
        "raw_mean_group1": 3.0,
        "raw_mean_group2": 2.0,
    }


def test_scvi_model_de_summary_references_are_not_shared(scvi_evidence_fixture):
    adata, _, model = scvi_evidence_fixture
    _, first = _run(adata, model)
    original_citation = first["references"][0]["citation"]
    first["references"][0]["citation"] = "caller mutation"

    _, second = _run(adata, model)

    assert second["references"][0]["citation"] == original_citation
    with pytest.raises(TypeError):
        SCVI_DE_REFERENCES[0]["citation"] = "global mutation"


def test_scvi_model_de_warns_when_registered_count_state_is_unknown(scvi_evidence_fixture):
    adata, raw_model, _ = scvi_evidence_fixture
    model = SCVIModel(
        raw_model,
        adata,
        current_model_parameters(
            technical_batch_key="batch",
            categorical_covariates=[],
            continuous_covariates=[],
            size_factor_key="",
            count_source_state="unknown",
            count_source_state_evidence=None,
        ),
    )

    _, summary = _run(adata, model)

    assert summary["model_evidence"]["registered_count_state"] == "unknown"
    assert any("not verified as raw counts" in warning for warning in summary["warnings"])


def test_scvi_model_de_rejects_unaudited_scvi_version_before_backend(scvi_evidence_fixture, monkeypatch):
    adata, raw_model, model = scvi_evidence_fixture
    monkeypatch.setattr(sys.modules["scvi"], "__version__", "1.6.0")

    with pytest.raises(RuntimeError, match=r">=1\.5,<1\.6"):
        _run(adata, model)

    assert raw_model.calls == []


def test_scvi_model_de_rejects_kwargs_only_backend_signature(scvi_evidence_fixture, science):
    class KwargsOnlySCVI(FullFakeSCVI):
        def differential_expression(self, **kwargs):
            return super().differential_expression(**kwargs)

    adata, _, _ = scvi_evidence_fixture
    raw_model = KwargsOnlySCVI(adata, science)
    model = SCVIModel(
        raw_model,
        adata,
        current_model_parameters(
            technical_batch_key="batch",
            categorical_covariates=[],
            continuous_covariates=[],
            size_factor_key="",
        ),
    )

    with pytest.raises(RuntimeError, match="does not support required arguments"):
        _run(adata, model)

    assert raw_model.calls == []


def test_scvi_model_de_backend_failure_cleans_only_analysis_manager(scvi_evidence_fixture):
    adata, raw_model, model = scvi_evidence_fixture
    raw_model.backend_error = ValueError("backend evidence failure")

    with pytest.raises(ValueError, match="backend evidence failure"):
        _run(adata, model)

    assert len(raw_model.deregister_calls) == 1
    assert raw_model.deregister_calls[0] is raw_model.calls[0]["adata"]
    assert raw_model.manager_store == {
        id(adata): adata,
        id(raw_model.unrelated_manager): raw_model.unrelated_manager,
    }
    assert raw_model.class_manager_store == raw_model.manager_store
