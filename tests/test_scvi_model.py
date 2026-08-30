from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

from openbio_singlecell.node_types import AnnDataType, SCVIModelType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_differential import (
    OpenBioSingleCellSCVIDifferentialExpression as _SCVIDifferentialExpressionSchema,
)
from openbio_singlecell.nodes_integration import OpenBioSingleCellSCVIIntegration
from openbio_singlecell.operations_integration import scvi_integration
from openbio_singlecell.scvi_model import SCVIModel
from tests.artifact_operation_harness import run_anndata_operation


def output_values(node_output):
    return node_output.result


def _run_scvi(adata, **overrides):
    parameters = {
        "source": {"source": "layer", "counts_layer": "counts"},
        "technical_batch_key": "batch",
        "categorical_covariates": "",
        "continuous_covariates": "",
        "n_latent": 10,
        "gene_likelihood": "zinb",
        "n_layers": 1,
        "dispersion": "gene",
        "dropout_rate": 0.1,
        "epochs": {"epochs": "automatic"},
        "early_stopping": False,
        "train_size": 0.9,
        "batch_size": 128,
        "size_factor_key": "",
        "accelerator": "auto",
        "output_key": "X_scVI",
        "overwrite_existing": False,
        "random_seed": 0,
    }
    parameters.update(overrides)
    return run_anndata_operation(scvi_integration, adata, parameters, artifact_output="model")


def current_training_parameters(**overrides):
    parameters = {
        "source": "X",
        "count_source_state": "counts",
        "count_source_state_evidence": "test fixture raw counts",
    }
    parameters.update(overrides)
    return parameters


@pytest.fixture
def adata(science):
    obs = science.pd.DataFrame(
        {"batch": science.pd.Categorical(["one", "one", "one", "two", "two", "two"])},
        index=[f"cell_{index}" for index in range(6)],
    )
    var = science.pd.DataFrame(index=[f"gene_{index}" for index in range(4)])
    value = science.ad.AnnData(science.np.arange(24, dtype=float).reshape(6, 4), obs=obs, var=var)
    value.layers["counts"] = value.X.copy()
    return value


class FakeTrainedSCVI:
    def __init__(self, registered_adata, science):
        self.adata = registered_adata
        self.is_trained = True
        self.science = science
        self.differential_expression_calls = []
        self.deregister_manager_calls = []
        self.differential_expression_error = None
        self.deregister_manager_error = None
        self.train_calls = 0

    def train(self, **kwargs):
        self.train_calls += 1
        raise AssertionError("Differential expression must not retrain the supplied scVI model.")

    def differential_expression(self, **kwargs):
        analysis_adata = kwargs["adata"]
        self.differential_expression_calls.append(kwargs)
        analysis_adata.obs["_fake_internal_mutation"] = True
        if self.differential_expression_error is not None:
            raise self.differential_expression_error
        return self.science.pd.DataFrame(
            {"lfc_mean": [1.5, -0.5], "proba_de": [0.95, 0.8]},
            index=["gene_0", "gene_1"],
        )

    def deregister_manager(self, analysis_adata=None):
        self.deregister_manager_calls.append(analysis_adata)
        if self.deregister_manager_error is not None:
            raise self.deregister_manager_error


def test_scvi_schemas_use_the_concrete_model_wire():
    integration = OpenBioSingleCellSCVIIntegration.define_schema()
    integration_inputs = {input_.id: input_ for input_ in integration.inputs}
    assert [(output.display_name, output.io_type) for output in integration.outputs] == [
        ("adata", AnnDataType.io_type),
        ("model", SCVIModelType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert integration_inputs["source"].get_io_type() == "COMFY_DYNAMICCOMBO_V3"
    assert integration_inputs["source"].advanced is not True
    assert [(option.key, [item.id for item in option.inputs]) for option in integration_inputs["source"].options] == [
        ("layer", ["counts_layer"]),
        ("X", []),
    ]
    assert integration_inputs["technical_batch_key"].default == "batch"
    assert integration_inputs["size_factor_key"].advanced is True
    assert integration_inputs["n_layers"].default == 1
    assert integration_inputs["dropout_rate"].default == 0.1
    assert integration_inputs["dispersion"].default == "gene"
    assert "compute_mde" not in integration_inputs
    assert "store_latent_distribution" not in integration_inputs

    differential = _SCVIDifferentialExpressionSchema.define_schema()
    differential_inputs = {input_.id: input_ for input_ in differential.inputs}
    assert [input_.id for input_ in differential.inputs] == [
        "adata",
        "model",
        "groupby",
        "group1",
        "group2",
        "population_scope",
        "mode",
        "batch_handling",
        "n_samples_overall",
        "random_seed",
    ]
    assert differential.inputs[0].io_type == AnnDataType.io_type
    assert differential.inputs[1].io_type == SCVIModelType.io_type
    assert [(output.display_name, output.io_type) for output in differential.outputs] == [
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert differential_inputs["population_scope"].get_io_type() == "COMFY_DYNAMICCOMBO_V3"
    assert [(option.key, [item.id for item in option.inputs]) for option in differential_inputs["population_scope"].options] == [
        ("all", []),
        ("obs_value", ["subset_column", "subset_value"]),
    ]
    assert differential_inputs["mode"].get_io_type() == "COMFY_DYNAMICCOMBO_V3"
    assert [(option.key, [item.id for item in option.inputs]) for option in differential_inputs["mode"].options] == [
        ("change", ["delta", "fdr_target"]),
        ("vanilla", []),
    ]
    assert differential_inputs["batch_handling"].default == "shared_technical_batches"
    assert differential_inputs["n_samples_overall"].advanced is True
    assert differential_inputs["random_seed"].advanced is True


def test_scvi_integration_publishes_the_trained_native_model(adata, science, monkeypatch):
    class FakeSCVI:
        setup_calls = []
        instances = []

        @classmethod
        def setup_anndata(cls, registered_adata, **kwargs):
            cls.setup_calls.append((registered_adata, kwargs))

        def __init__(self, registered_adata, **kwargs):
            self.adata = registered_adata
            self.is_trained = False
            self.init_parameters = kwargs
            self.train_parameters = None
            self.n_latent = kwargs["n_latent"]
            self.latent_calls = []
            self.history = {
                "elbo_train": science.np.asarray([12.0, 8.0, 6.0]),
                "elbo_validation": science.np.asarray([13.0, 9.0, 7.0]),
            }
            self.train_indices = list(range(5))
            self.validation_indices = [5]
            self.test_indices = []
            self.device = "cpu"
            type(self).instances.append(self)

        def train(self, **kwargs):
            self.train_parameters = kwargs
            self.is_trained = True

        def get_latent_representation(self, *, give_mean=True, return_dist=False):
            self.latent_calls.append((give_mean, return_dist))
            return science.np.arange(self.adata.n_obs * self.n_latent, dtype=float).reshape(
                self.adata.n_obs, self.n_latent
            )

        def differential_expression(self, **kwargs):
            raise AssertionError("Integration must not run differential expression.")

        def deregister_manager(self, analysis_adata):
            raise AssertionError("Integration must not create temporary AnnData managers.")

        def save(self, path, *, overwrite, save_anndata):
            destination = Path(path)
            destination.mkdir()
            (destination / "model.pt").write_bytes(b"native-scvi-model")
            assert overwrite is False
            assert save_anndata is True

    fake_scvi = types.ModuleType("scvi")
    fake_scvi.settings = types.SimpleNamespace(seed=None)
    fake_scvi.model = types.SimpleNamespace(SCVI=FakeSCVI)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    adata.obs["total_counts"] = science.np.asarray(adata.X.sum(axis=1)).ravel()
    output, model_artifact, report, code = output_values(
        _run_scvi(
            adata,
            source={"source": "X"},
            technical_batch_key="batch",
            size_factor_key="total_counts",
            n_latent=2,
            epochs={"epochs": "fixed", "max_epochs": 3},
            random_seed=17,
        )
    )

    raw_model = FakeSCVI.instances[0]
    assert model_artifact == {
        "kind": "OPENBIO_SCVI_MODEL",
        "codec": "scvi-native-directory",
        "members": ("model.pt",),
    }
    assert raw_model.adata is not output
    assert list(raw_model.adata.obs_names) == list(output.obs_names)
    assert raw_model.is_trained is True
    assert raw_model.train_parameters == {
        "max_epochs": 3,
        "accelerator": "auto",
        "devices": 1,
        "train_size": 0.9,
        "validation_size": None,
        "batch_size": 128,
        "early_stopping": False,
        "early_stopping_monitor": "elbo_validation",
        "early_stopping_patience": 10,
        "early_stopping_min_delta": 0.0,
        "check_val_every_n_epoch": 1,
    }
    setup_adata, setup_kwargs = FakeSCVI.setup_calls[0]
    assert setup_adata is raw_model.adata
    assert setup_adata is not output
    assert setup_kwargs == {
        "layer": None,
        "batch_key": "batch",
        "size_factor_key": "total_counts",
        "categorical_covariate_keys": None,
        "continuous_covariate_keys": None,
    }
    assert raw_model.latent_calls == [(True, False)]
    assert output.obsm["X_scVI"].shape == (adata.n_obs, 2)
    assert "X_scVI" not in adata.obsm
    assert report.summary["key_results"]["training"]["actual_epochs"] == 3
    assert report.summary["parameters"]["technical_batch_key"] == "batch"
    assert report.summary["parameters"]["max_epochs"] == 3
    assert report.summary["parameters"]["n_latent"] == 2
    assert report.summary["parameters"]["random_seed"] == 17
    assert report.summary["parameters"]["source"] == "X"
    assert "counts_layer" not in report.summary["parameters"]
    assert report.summary["parameters"]["size_factor_key"] == "total_counts"
    compile(code, "<scvi-code>", "exec")
    assert fake_scvi.settings.seed is None


def test_scvi_integration_restores_gradients_inside_comfyui_inference_mode(adata, science, monkeypatch):
    import torch

    inference_modes = []

    class GradientTrainingSCVI:
        instances = []

        @classmethod
        def setup_anndata(cls, registered_adata, **kwargs):
            inference_modes.append(("setup", torch.is_inference_mode_enabled()))

        def __init__(self, registered_adata, **kwargs):
            inference_modes.append(("init", torch.is_inference_mode_enabled()))
            self.adata = registered_adata
            self.is_trained = False
            self.n_latent = kwargs["n_latent"]
            self.weight = torch.nn.Parameter(torch.ones((1, registered_adata.n_vars)))
            self.history = {"elbo_train": science.np.asarray([2.0, 1.0])}
            self.device = "cpu"
            type(self).instances.append(self)

        def train(self, **kwargs):
            inference_modes.append(("train", torch.is_inference_mode_enabled()))
            inputs = torch.ones((1, self.adata.n_vars))
            loss = (inputs @ self.weight.T).square().sum()
            loss.backward()
            self.is_trained = True

        def get_latent_representation(self, *, give_mean=True, return_dist=False):
            inference_modes.append(("latent", torch.is_inference_mode_enabled()))
            return science.np.arange(self.adata.n_obs * self.n_latent, dtype=float).reshape(
                self.adata.n_obs, self.n_latent
            )

        def differential_expression(self, **kwargs):
            raise AssertionError("Integration must not run differential expression.")

        def deregister_manager(self, analysis_adata):
            raise AssertionError("Integration must not create temporary AnnData managers.")

        def save(self, path, *, overwrite, save_anndata):
            destination = Path(path)
            destination.mkdir()
            (destination / "model.pt").write_bytes(b"native-scvi-model")

    fake_scvi = types.ModuleType("scvi")
    fake_scvi.settings = types.SimpleNamespace(seed=None)
    fake_scvi.model = types.SimpleNamespace(SCVI=GradientTrainingSCVI)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    with torch.inference_mode():
        output, model_artifact, _, _ = output_values(
            _run_scvi(
                adata,
                source={"source": "X"},
                n_latent=2,
                epochs={"epochs": "fixed", "max_epochs": 1},
            )
        )
        assert torch.is_inference_mode_enabled() is True

    assert output.obsm["X_scVI"].shape == (adata.n_obs, 2)
    assert model_artifact["codec"] == "scvi-native-directory"
    assert GradientTrainingSCVI.instances[0].weight.grad is not None
    assert inference_modes == [
        ("setup", False),
        ("init", False),
        ("train", False),
        ("latent", False),
    ]


def test_scvi_model_rejects_untrained_or_detached_models(adata, science):
    untrained = FakeTrainedSCVI(adata, science)
    untrained.is_trained = False
    with pytest.raises(ValueError, match="requires a trained"):
        SCVIModel(untrained, adata, current_training_parameters())

    detached = FakeTrainedSCVI(adata.copy(), science)
    with pytest.raises(ValueError, match="attached to the registered"):
        SCVIModel(detached, adata, current_training_parameters())


def test_scvi_model_snapshots_training_parameters(adata, science):
    parameters = current_training_parameters(n_latent=2, categorical_covariates=["batch"])
    model = SCVIModel(FakeTrainedSCVI(adata, science), adata, parameters)
    parameters["categorical_covariates"].append("condition")
    visible_parameters = model.training_parameters
    visible_parameters["categorical_covariates"].append("donor")

    assert model.training_parameters["categorical_covariates"] == ["batch"]


@pytest.mark.parametrize("state", ["normalized", "logged", "scaled", "transformed", "pearson_residuals", "derived"])
def test_scvi_model_preserves_non_count_expression_state_as_evidence(adata, science, state):
    model = SCVIModel(
        FakeTrainedSCVI(adata, science),
        adata,
        current_training_parameters(
            count_source_state=state,
            count_source_state_evidence="explicit test evidence",
        ),
    )

    assert model.training_parameters["count_source_state"] == state
    assert model.evidence["registered_count_state"] == state
    assert model.evidence["registered_count_state_evidence"] == "explicit test evidence"


def test_scvi_model_snapshots_training_diagnostics(adata, science):
    diagnostics = {"actual_epochs": 2, "metrics": {"elbo_train": {"last": 1.0}}}
    model = SCVIModel(
        FakeTrainedSCVI(adata, science),
        adata,
        current_training_parameters(),
        diagnostics=diagnostics,
    )
    diagnostics["metrics"]["elbo_train"]["last"] = 999.0
    visible = model.diagnostics
    visible["metrics"]["elbo_train"]["last"] = -1.0

    assert model.diagnostics["metrics"]["elbo_train"]["last"] == 1.0
