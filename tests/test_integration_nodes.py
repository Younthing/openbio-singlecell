from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import os
import random
import subprocess
import sys
import types
from pathlib import Path

import pytest
import torch

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.node_types import AnnDataType, SCVIModelType, SummaryResultType
from openbio_singlecell.nodes_integration import (
    OpenBioSingleCellHarmonyIntegration as HarmonyIntegrationNode,
)
from openbio_singlecell.nodes_integration import (
    OpenBioSingleCellSCVIIntegration as SCVIIntegrationNode,
)
from openbio_singlecell.operations_input import ANNDATA_CODEC, ANNDATA_KIND
from openbio_singlecell.operations_integration import (
    harmony_integration,
    scvi_integration,
)
from openbio_singlecell.worker_protocol import OperationContext
from tests.artifact_operation_harness import run_anndata_operation


def _run_harmony(adata, **overrides):
    parameters = {
        "technical_batch_keys": "batch",
        "basis": "X_pca",
        "theta": {"theta": "automatic"},
        "ridge_penalty": -1.0,
        "sigma": 0.1,
        "n_clusters": 0,
        "tau": 0.0,
        "adjusted_basis": "X_pca_harmony",
        "overwrite_existing": False,
        "max_iter_harmony": 10,
        "max_iter_kmeans": 4,
        "random_seed": 0,
    }
    parameters.update(overrides)
    return run_anndata_operation(harmony_integration, adata, parameters).result


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
    return run_anndata_operation(scvi_integration, adata, parameters, artifact_output="model").result


def output_values(node_output):
    return node_output


@pytest.fixture(autouse=True)
def _reviewed_harmonypy_distribution(monkeypatch):
    original_version = importlib.metadata.version

    def version(name):
        return "2.0.0" if name == "harmonypy" else original_version(name)

    monkeypatch.setattr(importlib.metadata, "version", version)


@pytest.fixture
def integration_adata(science):
    counts = science.np.arange(1, 60 * 8 + 1, dtype=float).reshape(60, 8) % 17 + 1
    obs = science.pd.DataFrame(
        {
            "batch": science.pd.Categorical(["lane_a"] * 30 + ["lane_b"] * 30),
            "protocol": science.pd.Categorical(["v1", "v2"] * 30),
            "percent_mito": science.np.linspace(1.0, 9.0, 60),
        },
        index=[f"cell_{index}" for index in range(60)],
    )
    var = science.pd.DataFrame(index=[f"gene_{index}" for index in range(8)])
    adata = science.ad.AnnData(counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.obsm["X_pca"] = science.np.arange(60 * 5, dtype=float).reshape(60, 5) / 10
    metadata = ensure_metadata(adata)
    metadata["analysis_history"] = {
        "snapshot": {
            "operation": "snapshot_expression",
            "parameters": {"source": "X", "layer_name": "counts"},
        }
    }
    adata.uns["openbio_singlecell"] = metadata
    return adata


def _fake_harmonypy(science):
    module = types.ModuleType("harmonypy")
    module.__version__ = "2.0.0"
    module.calls = []
    module.mutate_input = False
    theta_sentinel = object()

    def run_harmony(
        data_mat,
        meta_data,
        vars_use,
        theta=theta_sentinel,
        lamb=None,
        sigma=0.1,
        nclust=None,
        tau=0,
        block_size=0.05,
        max_iter_harmony=10,
        max_iter_kmeans=4,
        epsilon_cluster=1e-3,
        epsilon_harmony=1e-2,
        verbose=True,
        random_state=0,
        ncores=0,
    ):
        module.calls.append(
            {
                "shape": tuple(data_mat.shape),
                "obs_names": tuple(meta_data.index),
                "vars_use": list(vars_use),
                "theta": None if theta is theta_sentinel else theta,
                "theta_provided": theta is not theta_sentinel,
                "lamb": lamb,
                "sigma": science.np.asarray(sigma, dtype=float).tolist(),
                "nclust": nclust,
                "tau": tau,
                "block_size": block_size,
                "max_iter_harmony": max_iter_harmony,
                "max_iter_kmeans": max_iter_kmeans,
                "epsilon_cluster": epsilon_cluster,
                "epsilon_harmony": epsilon_harmony,
                "verbose": verbose,
                "random_state": random_state,
                "ncores": ncores,
            }
        )
        if module.mutate_input:
            data_mat[0, 0] = -12345.0
        corrected = science.np.asarray(data_mat, dtype=float) + 0.25
        return types.SimpleNamespace(
            Z_corr=corrected,
            objective_harmony=[10.0, 7.0, 6.0],
            objective_kmeans=[10.0, 8.0, 7.0, 6.5, 6.0],
            kmeans_rounds=[2, 2],
        )

    module.run_harmony = run_harmony
    return module


def test_harmony_schema_and_fake_backend_contract(integration_adata, science, monkeypatch):
    schema = HarmonyIntegrationNode.define_schema()
    inputs = {item.id: item for item in schema.inputs}
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("adata", AnnDataType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert inputs["technical_batch_keys"].default == "batch"
    assert inputs["max_iter_kmeans"].default == 4
    assert inputs["random_seed"].advanced is True
    assert inputs["sigma"].min is None
    assert inputs["sigma"].max is None
    assert inputs["n_clusters"].max is None
    assert inputs["max_iter_harmony"].max is None
    assert inputs["max_iter_kmeans"].max is None
    assert [(option.key, [item.id for item in option.inputs]) for option in inputs["theta"].options] == [
        ("automatic", []),
        ("custom", ["theta_value"]),
    ]

    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)
    original_basis = integration_adata.obsm["X_pca"].copy()
    output, report, code = output_values(
        _run_harmony(
            integration_adata,
            technical_batch_keys="batch, protocol",
            n_clusters=3,
            random_seed=19,
        )
    )

    assert len(fake.calls) == 1
    assert fake.calls[0] == {
        "shape": (60, 5),
        "obs_names": tuple(integration_adata.obs_names),
        "vars_use": ["batch", "protocol"],
        "theta": None,
        "theta_provided": False,
        "lamb": None,
        "sigma": [0.1, 0.1, 0.1],
        "nclust": 3,
        "tau": 0.0,
        "block_size": 0.05,
        "max_iter_harmony": 10,
        "max_iter_kmeans": 4,
        "epsilon_cluster": 0.001,
        "epsilon_harmony": 0.01,
        "verbose": False,
        "random_state": 19,
        "ncores": 1,
    }
    science.np.testing.assert_array_equal(output.obsm["X_pca"], original_basis)
    science.np.testing.assert_allclose(output.obsm["X_pca_harmony"], original_basis + 0.25)
    assert "X_pca_harmony" not in integration_adata.obsm
    assert report.summary["key_results"]["backend_orientation"] == "n_obs_by_components_preserved"
    assert report.summary["parameters"]["harmonypy_version"] == "2.0.0"
    assert report.summary["key_results"]["objective_final"] == 6.0
    assert report.summary["parameters"]["technical_batch_keys"] == ["batch", "protocol"]
    assert report.summary["parameters"]["theta_mode"] == "automatic"
    assert report.summary["parameters"]["theta"] is None
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<harmony-code>", "exec")

    namespace = {}
    exec(code, namespace)
    reproduced = namespace["run_harmony_integration"](integration_adata.copy())
    science.np.testing.assert_allclose(reproduced.obsm["X_pca_harmony"], output.obsm["X_pca_harmony"])
    science.np.testing.assert_array_equal(reproduced.obsm["X_pca"], original_basis)
    metadata_collision = integration_adata.copy()
    metadata_collision.uns["harmony"] = {"X_pca_harmony": {"prior": True}}
    with pytest.raises(ValueError, match="metadata already exists"):
        namespace["run_harmony_integration"](metadata_collision)


def test_harmony_custom_theta_is_explicit_in_runtime_summary_and_code(integration_adata, science, monkeypatch):
    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)

    output, report, code = output_values(
        _run_harmony(
            integration_adata,
            theta={"theta": "custom", "theta_value": 3.5},
            n_clusters=3,
        )
    )

    assert fake.calls[-1]["theta_provided"] is True
    assert fake.calls[-1]["theta"] == 3.5
    assert report.summary["parameters"]["theta_mode"] == "custom"
    assert report.summary["parameters"]["theta"] == 3.5
    namespace = {}
    exec(code, namespace)
    reproduced = namespace["run_harmony_integration"](integration_adata.copy())
    assert fake.calls[-1]["theta_provided"] is True
    assert fake.calls[-1]["theta"] == 3.5
    science.np.testing.assert_allclose(reproduced.obsm["X_pca_harmony"], output.obsm["X_pca_harmony"])


def test_harmony_square_backend_result_preserves_reviewed_200_orientation(science, monkeypatch):
    raw_z_corr = science.np.asarray(
        [
            [1.0, 2.0, 3.0, 4.0],
            [10.0, 11.0, 12.0, 13.0],
            [20.0, 21.0, 22.0, 23.0],
            [30.0, 31.0, 32.0, 33.0],
        ]
    )
    obs = science.pd.DataFrame(
        {"batch": science.pd.Categorical(["a", "a", "b", "b"])},
        index=[f"cell_{index}" for index in range(4)],
    )
    adata = science.ad.AnnData(
        science.np.ones((4, 3)),
        obs=obs,
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(3)]),
    )
    adata.obsm["X_pca"] = science.np.arange(16, dtype=float).reshape(4, 4)
    fake = types.ModuleType("harmonypy")
    fake.__version__ = "2.0.0"

    def run_harmony(data_mat, meta_data, vars_use, **kwargs):
        del data_mat, meta_data, vars_use, kwargs
        return types.SimpleNamespace(
            Z_corr=raw_z_corr.copy(),
            objective_harmony=[3.0, 2.0],
            objective_kmeans=[3.0, 2.0],
            kmeans_rounds=[1],
        )

    fake.run_harmony = run_harmony
    monkeypatch.setitem(sys.modules, "harmonypy", fake)

    runtime, report, code = output_values(_run_harmony(adata, n_clusters=2))
    namespace: dict[str, object] = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="Small Harmony Technical-batch levels"):
        generated = namespace["run_harmony_integration"](adata.copy())

    expected = raw_z_corr
    science.np.testing.assert_array_equal(runtime.obsm["X_pca_harmony"], expected)
    science.np.testing.assert_array_equal(generated.obsm["X_pca_harmony"], expected)
    assert report.summary["key_results"]["backend_orientation"] == "n_obs_by_components_preserved"
    assert report.summary["key_results"]["backend_z_corr_shape"] == [4, 4]
    assert "X_pca_harmony" not in adata.obsm


def test_harmony_runtime_and_generated_code_reject_legacy_transposed_backend_shape(
    integration_adata,
    science,
    monkeypatch,
):
    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)
    _, _, code = output_values(_run_harmony(integration_adata, n_clusters=3))
    namespace: dict[str, object] = {}
    exec(code, namespace)

    def wrong_orientation(data_mat, meta_data, vars_use, **kwargs):
        del meta_data, vars_use, kwargs
        return types.SimpleNamespace(
            Z_corr=science.np.asarray(data_mat, dtype=float).T.copy(),
            objective_harmony=[2.0, 1.0],
        )

    fake.run_harmony = wrong_orientation
    runners = (
        lambda value: _run_harmony(value, n_clusters=3),
        namespace["run_harmony_integration"],
    )
    for runner in runners:
        with pytest.raises(RuntimeError, match=r"Z_corr must preserve \(n_obs, n_components\)"):
            runner(integration_adata)
        assert "X_pca_harmony" not in integration_adata.obsm


def test_harmony_accepts_capability_checked_2_0_patch_and_rejects_identity_drift(
    integration_adata, science, monkeypatch
):
    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)

    fake.__version__ = "2.0.1"
    monkeypatch.setattr(importlib.metadata, "version", lambda name: "2.0.1" if name == "harmonypy" else "test")
    _, report, _ = output_values(_run_harmony(integration_adata, n_clusters=3))
    assert report.summary["parameters"]["harmonypy_version"] == "2.0.1"

    fake.__version__ = "0.2.0"
    monkeypatch.setattr(importlib.metadata, "version", lambda name: "0.2.0" if name == "harmonypy" else "test")
    with pytest.raises(RuntimeError, match=r"harmonypy 2\.0\.x"):
        _run_harmony(integration_adata, n_clusters=3)

    fake.__version__ = "2.0.0"
    monkeypatch.setattr(importlib.metadata, "version", lambda name: "0.2.0" if name == "harmonypy" else "test")
    with pytest.raises(RuntimeError, match="distribution reports '0.2.0'"):
        _run_harmony(integration_adata, n_clusters=3)


@pytest.mark.parametrize(
    ("mutator", "match"),
    [
        (lambda adata, np: adata.obsm.__setitem__("X_pca", np.full((adata.n_obs, 3), np.nan)), "non-finite"),
        (lambda adata, np: adata.obsm.__setitem__("X_pca", np.ones((adata.n_obs, 0))), "at least one component"),
        (lambda adata, np: adata.obsm.__setitem__("X_pca_harmony", np.ones((adata.n_obs, 3))), "already exists"),
        (lambda adata, np: adata.uns.__setitem__("harmony", "invalid"), "must be a mapping"),
        (
            lambda adata, np: adata.uns.__setitem__("harmony", {"X_pca_harmony": {"prior": True}}),
            "metadata already exists",
        ),
    ],
)
def test_harmony_rejects_invalid_state_before_backend(integration_adata, science, mutator, match):
    mutator(integration_adata, science.np)
    with pytest.raises(ValueError, match=match):
        _run_harmony(integration_adata, n_clusters=3)


@pytest.mark.parametrize(
    ("labels", "summary_text", "generated_text"),
    [
        (["only"] * 60, "only one observed level", "only one observed level"),
        (["singleton", *(["other"] * 59)], "singleton levels", "contains singleton levels"),
    ],
)
def test_harmony_advisory_small_strata_run_and_disclose(
    integration_adata, science, monkeypatch, labels, summary_text, generated_text
):
    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)
    integration_adata.obs = integration_adata.obs.assign(batch=labels)

    output, report, code = output_values(_run_harmony(integration_adata, n_clusters=3))
    assert output.obsm["X_pca_harmony"].shape == integration_adata.obsm["X_pca"].shape
    assert any(summary_text in warning for warning in report.summary["warnings"])

    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning) as caught:
        reproduced = namespace["run_harmony_integration"](integration_adata.copy())
    assert any(generated_text in str(item.message) for item in caught)
    science.np.testing.assert_allclose(reproduced.obsm["X_pca_harmony"], output.obsm["X_pca_harmony"])


def test_harmony_one_component_and_zero_variance_are_advisory_with_generated_parity(
    integration_adata, science, monkeypatch
):
    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)
    integration_adata.obsm["X_pca"] = science.np.ones((integration_adata.n_obs, 1), dtype=float)
    original = integration_adata.obsm["X_pca"].copy()

    output, report, code = output_values(_run_harmony(integration_adata, n_clusters=1))

    assert output.obsm["X_pca_harmony"].shape == (integration_adata.n_obs, 1)
    assert report.summary["key_results"]["constant_input_components"] == 1
    assert report.summary["key_results"]["constant_input_component_indices"] == [0]
    assert any("one-dimensional basis" in warning for warning in report.summary["warnings"])
    assert any("zero-variance component" in warning for warning in report.summary["warnings"])
    assert fake.calls[-1]["sigma"] == [0.1]
    science.np.testing.assert_array_equal(integration_adata.obsm["X_pca"], original)

    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning) as caught:
        generated = namespace["run_harmony_integration"](integration_adata.copy())
    messages = [str(item.message) for item in caught]
    assert any("one-dimensional basis" in message for message in messages)
    assert any("zero-variance component" in message for message in messages)
    science.np.testing.assert_allclose(generated.obsm["X_pca_harmony"], output.obsm["X_pca_harmony"])
    science.np.testing.assert_array_equal(generated.obsm["X_pca"], original)


@pytest.mark.parametrize("axis_case", ["no_features", "duplicate_feature_ids"])
def test_harmony_ignores_irrelevant_variable_axis(integration_adata, science, monkeypatch, axis_case):
    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)
    if axis_case == "no_features":
        value = integration_adata[:, []].copy()
    else:
        value = integration_adata.copy()
        value.var_names = ["duplicate", "duplicate", *[f"gene_{index}" for index in range(2, value.n_vars)]]

    output, report, code = output_values(_run_harmony(value, n_clusters=1))

    assert output.obsm["X_pca_harmony"].shape == value.obsm["X_pca"].shape
    assert report.summary["key_results"]["features"] == value.n_vars
    assert fake.calls[-1]["sigma"] == [0.1]
    namespace = {}
    exec(code, namespace)
    generated = namespace["run_harmony_integration"](value.copy())
    science.np.testing.assert_allclose(generated.obsm["X_pca_harmony"], output.obsm["X_pca_harmony"])


def test_harmony_cluster_safety_and_distinct_destination_fail_before_backend(integration_adata, science, monkeypatch):
    fake = _fake_harmonypy(science)
    monkeypatch.setitem(sys.modules, "harmonypy", fake)

    for clusters in (1, integration_adata.n_obs):
        _, report, _ = output_values(
            _run_harmony(integration_adata, n_clusters=clusters)
        )
        assert report.summary["parameters"]["n_clusters"] == clusters
        assert fake.calls[-1]["sigma"] == [0.1] * clusters

    calls = len(fake.calls)
    with pytest.raises(ValueError, match="cannot exceed the number of cells"):
        _run_harmony(integration_adata, n_clusters=integration_adata.n_obs + 1)
    with pytest.raises(ValueError, match="must differ from basis"):
        _run_harmony(
            integration_adata,
            adjusted_basis="X_pca",
            overwrite_existing=True,
            n_clusters=1,
        )
    assert len(fake.calls) == calls


def test_harmony_backend_receives_private_basis_copy(integration_adata, science, monkeypatch):
    fake = _fake_harmonypy(science)
    fake.mutate_input = True
    monkeypatch.setitem(sys.modules, "harmonypy", fake)
    original = integration_adata.obsm["X_pca"].copy()

    output, _, code = output_values(_run_harmony(integration_adata, n_clusters=2))

    science.np.testing.assert_array_equal(integration_adata.obsm["X_pca"], original)
    science.np.testing.assert_array_equal(output.obsm["X_pca"], original)
    namespace = {}
    exec(code, namespace)
    generated = namespace["run_harmony_integration"](integration_adata.copy())
    science.np.testing.assert_array_equal(integration_adata.obsm["X_pca"], original)
    science.np.testing.assert_array_equal(generated.obsm["X_pca"], original)


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"theta": 2.0}, TypeError, "DynamicCombo"),
        ({"theta": {"theta": "other"}}, ValueError, "theta mode"),
        ({"theta": {"theta": "custom", "theta_value": -1.0}}, ValueError, "custom theta"),
        ({"ridge_penalty": "automatic"}, ValueError, "ridge_penalty"),
        ({"ridge_penalty": True}, ValueError, "ridge_penalty"),
        ({"sigma": "0.1"}, ValueError, "sigma"),
        ({"sigma": True}, ValueError, "sigma"),
        ({"random_seed": 2**31}, ValueError, "random_seed"),
        ({"overwrite_existing": 1}, TypeError, "overwrite_existing"),
    ],
)
def test_harmony_rejects_invalid_parameters_before_backend(integration_adata, kwargs, error, match):
    with pytest.raises(error, match=match):
        _run_harmony(integration_adata, n_clusters=3, **kwargs)


def _fake_scvi(science):
    class FakeSCVI:
        setup_calls = []
        instances = []
        report_trained = True
        constant_latent = False

        @classmethod
        def setup_anndata(cls, adata, **kwargs):
            cls.setup_calls.append((adata, kwargs))

        def __init__(self, adata, **kwargs):
            self.adata = adata
            self.parameters = kwargs
            self.is_trained = False
            epoch = science.pd.Index([0, 1], name="epoch")
            self.history = {
                key: science.pd.DataFrame({key: values}, index=epoch.copy())
                for key, values in {
                    "elbo_train": [20.0, 12.0],
                    "elbo_validation": [22.0, 13.0],
                    "reconstruction_loss_train": [15.0, 9.0],
                    "reconstruction_loss_validation": [17.0, 10.0],
                    "kl_local_train": [5.0, 3.0],
                    "kl_local_validation": [5.0, 3.0],
                }.items()
            }
            self.train_indices = list(range(54))
            self.validation_indices = list(range(54, 60))
            self.test_indices = []
            self.device = "cpu"
            self.train_kwargs = None
            self.save_calls = []
            self.deregister_manager_calls = []
            type(self).instances.append(self)

        def train(self, **kwargs):
            random.random()
            science.np.random.random()
            torch.rand(1)
            self.train_kwargs = kwargs
            self.is_trained = type(self).report_trained

        def get_latent_representation(self):
            n_latent = self.parameters["n_latent"]
            if type(self).constant_latent:
                return science.np.zeros((self.adata.n_obs, n_latent), dtype=float)
            return science.np.arange(self.adata.n_obs * n_latent, dtype=float).reshape(self.adata.n_obs, n_latent)

        def differential_expression(self, **kwargs):
            return science.pd.DataFrame({"lfc_mean": [1.0]}, index=[self.adata.var_names[0]])

        def deregister_manager(self, adata=None):
            self.deregister_manager_calls.append(adata)

        def save(self, path, *, overwrite, save_anndata):
            destination = Path(path)
            destination.mkdir()
            (destination / "model.pt").write_bytes(b"native-scvi-model")
            self.save_calls.append({"path": str(destination), "overwrite": overwrite, "save_anndata": save_anndata})

    module = types.ModuleType("scvi")
    module.settings = types.SimpleNamespace(seed=123)
    module.model = types.SimpleNamespace(SCVI=FakeSCVI)
    return module, FakeSCVI


def test_scvi_operation_saves_official_worker_bound_native_directory(
    tmp_path,
    integration_adata,
    science,
    monkeypatch,
):
    fake_scvi, fake_class = _fake_scvi(science)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, integration_adata)
    input_payload = input_root / ANNDATA_PAYLOAD
    before = input_payload.read_bytes()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext(staging, "00000000-0000-4000-8000-000000000002")
    parameters = {
        "source": {"source": "layer", "counts_layer": "counts"},
        "technical_batch_key": "batch",
        "categorical_covariates": "",
        "continuous_covariates": "",
        "n_latent": 3,
        "gene_likelihood": "zinb",
        "n_layers": 1,
        "dispersion": "gene",
        "dropout_rate": 0.1,
        "epochs": {"epochs": "fixed", "max_epochs": 2},
        "early_stopping": False,
        "train_size": 0.9,
        "batch_size": 128,
        "size_factor_key": "",
        "accelerator": "cpu",
        "output_key": "X_scVI",
        "overwrite_existing": False,
        "random_seed": 5,
    }

    records = scvi_integration(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": ANNDATA_KIND,
                "codec": ANNDATA_CODEC,
                "path": str(input_root.resolve()),
            }
        },
        parameters,
    )

    assert [(record["type"], record["name"]) for record in records] == [
        ("artifact", "adata"),
        ("artifact", "model"),
        ("summary", "summary"),
        ("string", "code"),
    ]
    assert records[1] == {
        "type": "artifact",
        "name": "model",
        "kind": "OPENBIO_SCVI_MODEL",
        "codec": "scvi-native-directory",
        "payload": "outputs/model/native",
    }
    model_root = staging / records[1]["payload"]
    assert (model_root / "model.pt").read_bytes() == b"native-scvi-model"
    assert fake_class.instances[-1].save_calls == [{"path": str(model_root), "overwrite": False, "save_anndata": True}]
    assert input_payload.read_bytes() == before
    assert "X_scVI" not in read_anndata(input_root).obsm
    assert read_anndata(staging / records[0]["payload"]).obsm["X_scVI"].shape == (60, 3)
    json.dumps(records, allow_nan=False)


def test_scvi_fake_backend_reporting_code_and_private_model_state(integration_adata, science, monkeypatch):
    fake_scvi, fake_class = _fake_scvi(science)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    python_state = random.getstate()
    numpy_state = science.np.random.get_state()
    torch_state = torch.random.get_rng_state()
    output, model, report, code = output_values(
        _run_scvi(
            integration_adata,
            source={"source": "layer", "counts_layer": "counts"},
            technical_batch_key="batch",
            categorical_covariates="protocol",
            continuous_covariates="percent_mito",
            n_latent=3,
            epochs={"epochs": "fixed", "max_epochs": 12},
            early_stopping=True,
            accelerator="cpu",
            random_seed=7,
        )
    )

    assert random.getstate() == python_state
    restored_numpy_state = science.np.random.get_state()
    assert restored_numpy_state[0] == numpy_state[0]
    science.np.testing.assert_array_equal(restored_numpy_state[1], numpy_state[1])
    assert restored_numpy_state[2:] == numpy_state[2:]
    assert torch.equal(torch.random.get_rng_state(), torch_state)
    raw_model = fake_class.instances[0]
    assert raw_model.parameters == {
        "n_layers": 1,
        "n_latent": 3,
        "gene_likelihood": "zinb",
        "dispersion": "gene",
        "dropout_rate": 0.1,
    }
    assert raw_model.train_kwargs["max_epochs"] == 12
    assert raw_model.train_kwargs["accelerator"] == "cpu"
    assert raw_model.train_kwargs["devices"] == 1
    assert raw_model.train_kwargs["early_stopping"] is True
    assert raw_model.train_kwargs["early_stopping_monitor"] == "elbo_validation"
    assert raw_model.train_kwargs["early_stopping_patience"] == 10
    assert raw_model.train_kwargs["early_stopping_min_delta"] == 0.0
    assert raw_model.train_kwargs["check_val_every_n_epoch"] == 1
    assert fake_scvi.settings.seed == 123
    assert model == {
        "kind": "OPENBIO_SCVI_MODEL",
        "codec": "scvi-native-directory",
        "members": ("model.pt", "openbio-training-diagnostics.json"),
    }
    assert raw_model.adata is not output
    assert set(raw_model.adata.obs["batch"].astype(str)) == {"lane_a", "lane_b"}
    assert output.obsm["X_scVI"].shape == (60, 3)
    assert "X_scVI" not in integration_adata.obsm
    assert report.summary["parameters"]["epochs_mode"] == "fixed"
    assert report.summary["parameters"]["devices"] == 1
    assert report.summary["parameters"]["accelerator"] == "cpu"
    assert report.summary["key_results"]["training"]["device"] == "cpu"
    assert report.summary["key_results"]["training"]["actual_epochs"] == 2
    inferred_state_fields = {"count_source_state", "count_source_state_evidence"}
    assert inferred_state_fields.isdisjoint(report.summary["parameters"])
    assert inferred_state_fields.isdisjoint(report.summary["key_results"])
    assert report.summary["key_results"]["model_session_only"] is True
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<scvi-code>", "exec")

    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_model = namespace["run_scvi_integration"](integration_adata.copy())
    science.np.testing.assert_array_equal(reproduced.obsm["X_scVI"], output.obsm["X_scVI"])
    assert reproduced_model.train_kwargs == raw_model.train_kwargs
    collision = integration_adata.copy()
    collision.obsm["X_scVI"] = science.np.ones((collision.n_obs, 3))
    with pytest.raises(ValueError, match="already exists"):
        namespace["run_scvi_integration"](collision)
    zero_cell = integration_adata.copy()
    zero_cell.layers["counts"] = science.np.vstack(
        [science.np.zeros((1, zero_cell.n_vars)), science.np.asarray(zero_cell.layers["counts"])[1:]]
    )
    with pytest.raises(ValueError, match="cells with zero total counts"):
        namespace["run_scvi_integration"](zero_cell)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
@pytest.mark.filterwarnings("ignore::FutureWarning")
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_scvi_real_backend_open_boundaries_smoke():
    if importlib.util.find_spec("scvi") is None:
        pytest.skip("Real scVI smoke requires the optional scvi-tools backend.")

    plugin_root = Path(__file__).resolve().parents[1]
    comfy_root = next(Path(entry).resolve() for entry in sys.path if entry and (Path(entry) / "main.py").is_file())
    environment = os.environ.copy()
    python_paths = [str(plugin_root), str(comfy_root)]
    if environment.get("PYTHONPATH"):
        python_paths.append(environment["PYTHONPATH"])
    environment["PYTHONPATH"] = os.pathsep.join(python_paths)
    environment["OPENBIO_COMFYUI_ROOT"] = str(comfy_root)

    completed = subprocess.run(
        [sys.executable, str(Path(__file__).with_name("scvi_real_smoke.py"))],
        cwd=plugin_root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    diagnostic = completed.stdout + completed.stderr
    assert completed.returncode == 0, diagnostic
    prefix = "OPENBIO_SCVI_REAL_SMOKE="
    payload_lines = [line for line in completed.stdout.splitlines() if line.startswith(prefix)]
    assert len(payload_lines) == 1, diagnostic
    payload = json.loads(payload_lines[0][len(prefix) :])
    assert payload == {
        "backend_distribution": "scvi-tools",
        "backend_version": importlib.metadata.version("scvi-tools"),
        "code_compiles": True,
        "input_counts_unchanged": True,
        "input_has_no_output_embedding": True,
        "latent_shape": [32, 8],
        "model_class_endswith_scvi": True,
        "no_validation_holdout_warning": True,
        "overcomplete_warning": True,
        "status": "pass",
        "training_diagnostic_metrics": [
            "elbo_train",
            "reconstruction_loss_train",
            "kl_local_train",
            "kl_global_train",
        ],
        "training_plot_code_parity": True,
        "training_plot_png": True,
        "zero_total_genes": 1,
    }


def test_scvi_restores_mps_rng_state_in_runtime_and_generated_code(integration_adata, science, monkeypatch):
    fake_scvi, fake_class = _fake_scvi(science)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    state = {"value": "initial"}
    calls = []

    def get_rng_state():
        calls.append(("get", state["value"]))
        return state["value"]

    def set_rng_state(value):
        calls.append(("set", value))
        state["value"] = value

    original_train = fake_class.train

    def train_and_mutate_mps(self, **kwargs):
        original_train(self, **kwargs)
        state["value"] = "mutated"

    monkeypatch.setattr(torch.backends.mps, "is_available", lambda: True)
    monkeypatch.setattr(torch.mps, "get_rng_state", get_rng_state)
    monkeypatch.setattr(torch.mps, "set_rng_state", set_rng_state)
    monkeypatch.setattr(fake_class, "train", train_and_mutate_mps)

    _, _, _, code = output_values(
        _run_scvi(
            integration_adata,
            source={"source": "layer", "counts_layer": "counts"},
            n_latent=3,
            accelerator="mps",
        )
    )
    assert state["value"] == "initial"
    assert calls == [("get", "initial"), ("set", "initial")]

    calls.clear()
    namespace = {}
    exec(code, namespace)
    namespace["run_scvi_integration"](integration_adata)
    assert state["value"] == "initial"
    assert calls == [("get", "initial"), ("set", "initial")]


def test_scvi_runtime_and_generated_code_reject_untrained_backend(integration_adata, science, monkeypatch):
    fake_scvi, fake_class = _fake_scvi(science)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    _, _, _, code = output_values(
        _run_scvi(
            integration_adata,
            source={"source": "layer", "counts_layer": "counts"},
            n_latent=3,
            accelerator="cpu",
        )
    )
    fake_class.report_trained = False

    with pytest.raises(RuntimeError, match="did not report a trained model"):
        _run_scvi(
            integration_adata,
            source={"source": "layer", "counts_layer": "counts"},
            n_latent=3,
            accelerator="cpu",
        )

    namespace = {}
    exec(code, namespace)
    with pytest.raises(RuntimeError, match="did not report a trained model"):
        namespace["run_scvi_integration"](integration_adata)


@pytest.mark.parametrize(
    ("mutator", "match"),
    [
        (
            lambda adata, np: adata.layers.__setitem__(
                "counts",
                np.where(
                    (np.arange(adata.n_obs)[:, None] == 0) & (np.arange(adata.n_vars)[None, :] == 0),
                    -1.0,
                    np.asarray(adata.layers["counts"]),
                ),
            ),
            "negative",
        ),
        (
            lambda adata, np: adata.layers.__setitem__(
                "counts",
                np.where(
                    (np.arange(adata.n_obs)[:, None] == 0) & (np.arange(adata.n_vars)[None, :] == 0),
                    np.nan,
                    np.asarray(adata.layers["counts"]),
                ),
            ),
            "non-finite",
        ),
        (
            lambda adata, np: adata.layers.__setitem__(
                "counts", np.vstack([np.zeros((1, adata.n_vars)), np.asarray(adata.layers["counts"])[1:]])
            ),
            "zero total counts",
        ),
        (lambda adata, np: adata.obsm.__setitem__("X_scVI", np.ones((adata.n_obs, 2))), "already exists"),
    ],
)
def test_scvi_rejects_invalid_scientific_state_before_backend(integration_adata, science, mutator, match):
    mutator(integration_adata, science.np)
    kwargs = {"continuous_covariates": "percent_mito"} if "nonconstant" in match else {}
    with pytest.raises(ValueError, match=match):
        _run_scvi(integration_adata, n_latent=3, **kwargs)


@pytest.mark.parametrize(
    ("mutator", "kwargs", "summary_text", "generated_text"),
    [
        (
            lambda adata: setattr(adata, "obs", adata.obs.assign(batch=["one"] * adata.n_obs)),
            {},
            "only one observed level",
            "only one observed level",
        ),
        (
            lambda adata: setattr(adata, "obs", adata.obs.assign(percent_mito=[1.0] * adata.n_obs)),
            {"continuous_covariates": "percent_mito"},
            "Constant scVI continuous nuisance",
            "is constant",
        ),
    ],
)
def test_scvi_advisory_noninformative_covariates_run_and_disclose(
    integration_adata, science, monkeypatch, mutator, kwargs, summary_text, generated_text
):
    fake_scvi, _ = _fake_scvi(science)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    mutator(integration_adata)

    output, _, report, code = output_values(
        _run_scvi(integration_adata, n_latent=3, **kwargs)
    )
    assert output.obsm["X_scVI"].shape == (integration_adata.n_obs, 3)
    assert any(summary_text in warning for warning in report.summary["warnings"])

    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match=generated_text):
        generated, _ = namespace["run_scvi_integration"](integration_adata)
    science.np.testing.assert_allclose(generated.obsm["X_scVI"], output.obsm["X_scVI"])


@pytest.mark.parametrize(
    ("train_size", "train_warning"),
    [
        (0.2, "fewer than half of cells"),
        (1.0, "leaves no validation holdout"),
    ],
)
def test_scvi_open_expert_architecture_and_count_advisories_have_generated_parity(
    integration_adata,
    science,
    monkeypatch,
    train_size,
    train_warning,
):
    fake_scvi, fake_class = _fake_scvi(science)
    fake_class.constant_latent = True
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    selected = science.np.asarray(integration_adata.layers["counts"], dtype=float).copy() + 0.5
    selected[:, 0] = 0.0
    integration_adata.layers["counts"] = selected
    original = selected.copy()

    output, _, report, code = output_values(
        _run_scvi(
            integration_adata,
            source={"source": "layer", "counts_layer": "counts"},
            n_latent=10,
            n_layers=21,
            train_size=train_size,
            accelerator="cpu",
        )
    )

    expected_warning_texts = (
        "non-integer values",
        "all-zero gene",
        "overcomplete architecture",
        "unusually deep architecture",
        train_warning,
        "constant latent dimension",
    )
    for text in expected_warning_texts:
        assert any(text in warning for warning in report.summary["warnings"])
    results = report.summary["key_results"]
    assert results["count_source_integer_like"] is False
    assert results["zero_total_genes"] == 1
    assert results["constant_latent_dimensions"] == 10
    assert results["constant_latent_dimension_indices"] == list(range(10))
    expected_validation_interval = None if train_size == 1 else 1
    assert fake_class.instances[-1].train_kwargs["check_val_every_n_epoch"] == expected_validation_interval
    assert report.summary["parameters"]["check_val_every_n_epoch"] == expected_validation_interval
    assert "verified raw UMI" not in report.summary["methods"]
    json.dumps(report.summary, allow_nan=False)
    science.np.testing.assert_array_equal(integration_adata.layers["counts"], original)

    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning) as caught:
        generated, _ = namespace["run_scvi_integration"](integration_adata)
    generated_warnings = [str(item.message) for item in caught]
    for text in expected_warning_texts:
        assert any(text in warning for warning in generated_warnings)
    science.np.testing.assert_array_equal(generated.obsm["X_scVI"], output.obsm["X_scVI"])
    science.np.testing.assert_array_equal(integration_adata.layers["counts"], original)


def test_scvi_dropout_one_is_advisory_with_generated_parity(integration_adata, science, monkeypatch):
    fake_scvi, _ = _fake_scvi(science)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    output, _, report, code = output_values(
        _run_scvi(
            integration_adata,
            dropout_rate=1.0,
            n_latent=3,
            accelerator="cpu",
        )
    )

    assert output.obsm["X_scVI"].shape == (integration_adata.n_obs, 3)
    assert any("dropout_rate=1" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(code, namespace)
    with pytest.warns(UserWarning, match="dropout_rate=1"):
        generated, _ = namespace["run_scvi_integration"](integration_adata)
    science.np.testing.assert_array_equal(generated.obsm["X_scVI"], output.obsm["X_scVI"])


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"dropout_rate": "0.1"}, ValueError, "dropout_rate"),
        ({"dropout_rate": True}, ValueError, "dropout_rate"),
        ({"dropout_rate": -0.01}, ValueError, "dropout_rate"),
        ({"dropout_rate": 1.01}, ValueError, "dropout_rate"),
        ({"early_stopping": 1}, TypeError, "early_stopping"),
        ({"overwrite_existing": 1}, TypeError, "overwrite_existing"),
        ({"train_size": 0.0}, ValueError, "train_size"),
        ({"train_size": 1.01}, ValueError, "train_size"),
        ({"train_size": 1.0, "early_stopping": True}, ValueError, "validation holdout"),
    ],
)
def test_scvi_rejects_invalid_parameters_before_backend(integration_adata, kwargs, error, match):
    with pytest.raises(error, match=match):
        _run_scvi(integration_adata, n_latent=3, **kwargs)


def test_scvi_schema_is_narrow_and_uses_concrete_outputs():
    schema = SCVIIntegrationNode.define_schema()
    inputs = {item.id: item for item in schema.inputs}
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("adata", AnnDataType.io_type),
        ("model", SCVIModelType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert inputs["n_layers"].default == 1
    assert inputs["dispersion"].default == "gene"
    assert inputs["dropout_rate"].default == 0.1
    assert inputs["dropout_rate"].max == 1.0
    assert inputs["early_stopping"].default is False
    assert inputs["n_latent"].max is None
    assert inputs["n_layers"].max is None
    assert inputs["batch_size"].max is None
    assert inputs["train_size"].min == 0.0
    assert inputs["train_size"].max == 1.0
    assert "compute_mde" not in inputs
    assert "store_latent_distribution" not in inputs
    assert [(option.key, [item.id for item in option.inputs]) for option in inputs["epochs"].options] == [
        ("automatic", []),
        ("fixed", ["max_epochs"]),
    ]
    fixed_epochs = next(option for option in inputs["epochs"].options if option.key == "fixed")
    assert fixed_epochs.inputs[0].max is None
