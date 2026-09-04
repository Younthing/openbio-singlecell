from __future__ import annotations

import hashlib
import importlib.metadata
import sys
import types
import uuid
from pathlib import Path

import pytest

from openbio_singlecell.artifact_codecs import read_anndata, read_plot, write_anndata
from openbio_singlecell.node_types import AnnDataType, PlotResultType, SummaryResultType
from openbio_singlecell.nodes_integration import (
    OpenBioSingleCellHarmonyConvergencePlot,
    OpenBioSingleCellSCVITrainingPlot,
)
from openbio_singlecell.operations_integration import (
    harmony_convergence_plot,
    harmony_convergence_plot_owned,
    harmony_integration,
    scvi_integration,
    scvi_training_plot,
    scvi_training_plot_owned,
)
from openbio_singlecell.scvi_model import (
    build_scvi_training_diagnostics,
    read_scvi_training_diagnostics,
    write_scvi_training_diagnostics,
)
from openbio_singlecell.worker_protocol import OperationContext
from tests.artifact_operation_harness import run_anndata_operation


def _harmony_input(science):
    adata = science.ad.AnnData(science.np.ones((6, 2), dtype=float))
    adata.obs_names = [f"cell-{index}" for index in range(6)]
    adata.var_names = ["gene-a", "gene-b"]
    adata.obs["batch"] = science.pd.Categorical(["a", "a", "a", "b", "b", "b"])
    adata.obsm["X_pca"] = science.np.arange(18, dtype=float).reshape(6, 3)
    return adata


def _harmony_parameters():
    return {
        "technical_batch_keys": "batch",
        "basis": "X_pca",
        "theta": {"theta": "automatic"},
        "ridge_penalty": -1.0,
        "sigma": 0.1,
        "n_clusters": 2,
        "tau": 0.0,
        "adjusted_basis": "X_pca_harmony",
        "overwrite_existing": False,
        "max_iter_harmony": 10,
        "max_iter_kmeans": 4,
        "random_seed": 3,
    }


def _install_fake_harmony(monkeypatch, science):
    module = types.ModuleType("harmonypy")
    module.__version__ = "2.0.0"

    def run_harmony(
        data_mat,
        meta_data,
        vars_use,
        *,
        lamb,
        sigma,
        nclust,
        tau,
        block_size,
        max_iter_harmony,
        max_iter_kmeans,
        epsilon_cluster,
        epsilon_harmony,
        verbose,
        random_state,
        ncores,
    ):
        del (
            meta_data,
            vars_use,
            lamb,
            sigma,
            nclust,
            tau,
            block_size,
            max_iter_harmony,
            max_iter_kmeans,
            epsilon_cluster,
            epsilon_harmony,
            verbose,
            random_state,
            ncores,
        )
        return types.SimpleNamespace(
            Z_corr=science.np.asarray(data_mat, dtype=float) + 0.5,
            objective_harmony=[12.0, 8.5, 7.25, 7.0],
            objective_kmeans=[12.0, 10.0, 8.5, 7.25, 7.1, 7.0],
            kmeans_rounds=[2, 1, 2],
        )

    module.run_harmony = run_harmony
    monkeypatch.setitem(sys.modules, "harmonypy", module)
    original_version = importlib.metadata.version
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: "2.0.0" if name == "harmonypy" else original_version(name),
    )


def _scvi_parameters():
    return {
        "source": {"source": "X"},
        "technical_batch_key": "batch",
        "categorical_covariates": "",
        "continuous_covariates": "",
        "n_latent": 2,
        "gene_likelihood": "zinb",
        "n_layers": 1,
        "dispersion": "gene",
        "dropout_rate": 0.1,
        "epochs": {"epochs": "fixed", "max_epochs": 3},
        "early_stopping": False,
        "train_size": 0.9,
        "batch_size": 128,
        "size_factor_key": "",
        "accelerator": "cpu",
        "output_key": "X_scVI",
        "overwrite_existing": False,
        "random_seed": 4,
    }


def _install_fake_scvi(monkeypatch, science):
    class FakeSCVI:
        @classmethod
        def setup_anndata(cls, registered_adata, **kwargs):
            del registered_adata, kwargs

        def __init__(self, registered_adata, **kwargs):
            self.adata = registered_adata
            self.n_latent = kwargs["n_latent"]
            self.is_trained = False
            epoch = science.pd.Index([0, 1, 2], name="epoch")
            values = {
                "elbo_train": [12.0, 8.0, 6.0],
                "elbo_validation": [13.0, 9.0, 7.0],
                "reconstruction_loss_train": [10.0, 7.0, 5.5],
                "reconstruction_loss_validation": [11.0, 8.0, 6.5],
                "kl_local_train": [2.0, 1.0, 0.5],
                "kl_local_validation": [2.2, 1.2, 0.6],
                "kl_global_train": [0.3, 0.2, 0.1],
                "kl_global_validation": [0.35, 0.25, 0.15],
            }
            self.history = {
                key: science.pd.DataFrame({key: metric_values}, index=epoch.copy())
                for key, metric_values in values.items()
            }
            self.train_indices = science.np.arange(5)
            self.validation_indices = science.np.asarray([5])
            self.test_indices = science.np.asarray([], dtype=int)
            self.device = "cpu"

        def train(self, **kwargs):
            del kwargs
            self.is_trained = True

        def get_latent_representation(self):
            return science.np.arange(self.adata.n_obs * self.n_latent, dtype=float).reshape(
                self.adata.n_obs, self.n_latent
            )

        def differential_expression(self, **kwargs):
            del kwargs
            raise AssertionError("integration must not run differential expression")

        def deregister_manager(self, adata=None):
            del adata

        def save(self, path, *, overwrite, save_anndata):
            assert overwrite is False
            assert save_anndata is True
            destination = Path(path)
            destination.mkdir()
            (destination / "model.pt").write_bytes(b"native-model")

    module = types.ModuleType("scvi")
    module.settings = types.SimpleNamespace(seed=None)
    module.model = types.SimpleNamespace(SCVI=FakeSCVI)
    monkeypatch.setitem(sys.modules, "scvi", module)
    return module


def test_harmony_retains_complete_objective_evidence_tied_to_adjusted_basis(science, monkeypatch):
    _install_fake_harmony(monkeypatch, science)
    adata = _harmony_input(science)

    output, report, code = run_anndata_operation(
        harmony_integration,
        adata,
        _harmony_parameters(),
    ).result

    evidence = output.uns["harmony"]["X_pca_harmony"]
    assert evidence["schema_version"] == "openbio-singlecell/harmony-diagnostics/v1"
    assert evidence["payload"]["adjusted_basis"] == "X_pca_harmony"
    assert evidence["payload"]["basis"] == "X_pca"
    assert evidence["payload"]["harmony_objective"] == [
        {"iteration": 0, "value": 12.0},
        {"iteration": 1, "value": 8.5},
        {"iteration": 2, "value": 7.25},
        {"iteration": 3, "value": 7.0},
    ]
    assert evidence["payload"]["kmeans_objective"] == [
        {"update": index, "value": value}
        for index, value in enumerate([12.0, 10.0, 8.5, 7.25, 7.1, 7.0])
    ]
    assert evidence["payload"]["outer_rounds"] == [
        {"round": 1, "kmeans_updates": 2},
        {"round": 2, "kmeans_updates": 1},
        {"round": 3, "kmeans_updates": 2},
    ]
    assert len(evidence["payload"]["input_basis_fingerprint_sha256"]) == 64
    assert len(evidence["payload"]["adjusted_basis_fingerprint_sha256"]) == 64
    assert len(evidence["payload"]["observation_axis_fingerprint_sha256"]) == 64
    assert len(evidence["payload_sha256"]) == 64
    assert report.summary["key_results"]["harmony_diagnostic_evidence"]["payload_sha256"] == evidence[
        "payload_sha256"
    ]

    namespace: dict[str, object] = {}
    exec(code, namespace)
    reproduced = namespace["run_harmony_integration"](adata.copy())
    reproduced_evidence = reproduced.uns["harmony"]["X_pca_harmony"]
    assert reproduced_evidence == evidence


def test_harmony_convergence_plot_schema_is_a_read_only_companion():
    schema = OpenBioSingleCellHarmonyConvergencePlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellHarmonyConvergencePlot"
    assert schema.display_name == "Harmony Convergence Plot"
    assert schema.category == "openbio/single-cell/batch-integration"
    assert [item.id for item in schema.inputs] == ["adata", "adjusted_basis", "max_image_pixels"]
    assert schema.inputs[0].io_type == AnnDataType.io_type
    assert schema.inputs[1].default == "X_pca_harmony"
    assert schema.inputs[2].default == 40_000_000
    assert schema.inputs[2].advanced is True
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("plot", PlotResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_harmony_convergence_plot_reads_both_stored_objective_sequences_without_reintegration(
    tmp_path, science, monkeypatch
):
    import matplotlib

    _install_fake_harmony(monkeypatch, science)
    output, _, _ = run_anndata_operation(
        harmony_integration,
        _harmony_input(science),
        _harmony_parameters(),
    ).result
    original_basis = output.obsm["X_pca"].copy()
    original_adjusted = output.obsm["X_pca_harmony"].copy()
    original_evidence_sha = output.uns["harmony"]["X_pca_harmony"]["payload_sha256"]
    sys.modules["harmonypy"].run_harmony = lambda *_args, **_kwargs: (_ for _ in ()).throw(
        AssertionError("plot must not rerun Harmony")
    )
    numpy_state = science.np.random.get_state()
    rc_state = {key: matplotlib.rcParams[key] for key in ("figure.dpi", "font.size", "axes.grid")}

    plot, report, code = harmony_convergence_plot_owned(output)

    assert plot.png.startswith(b"\x89PNG\r\n\x1a\n")
    details = report.summary["key_results"]
    assert details["harmony_objective_points"] == 4
    assert details["kmeans_objective_points"] == 6
    assert details["kmeans_rounds"] == [2, 1, 2]
    assert details["outer_round_boundaries"] == [0, 2, 3, 5]
    assert details["boundary_line_legend"] == "Outer-round boundary"
    assert details["payload_sha256"] == original_evidence_sha
    assert report.summary["parameters"] == {
        "adjusted_basis": "X_pca_harmony",
        "max_image_pixels": 40_000_000,
    }
    namespace: dict[str, object] = {}
    exec(code, namespace)
    assert namespace["plot_harmony_convergence"](output) == plot.png
    science.np.testing.assert_array_equal(output.obsm["X_pca"], original_basis)
    science.np.testing.assert_array_equal(output.obsm["X_pca_harmony"], original_adjusted)
    assert output.uns["harmony"]["X_pca_harmony"]["payload_sha256"] == original_evidence_sha
    restored_numpy_state = science.np.random.get_state()
    assert restored_numpy_state[0] == numpy_state[0]
    science.np.testing.assert_array_equal(restored_numpy_state[1], numpy_state[1])
    assert restored_numpy_state[2:] == numpy_state[2:]
    assert {key: matplotlib.rcParams[key] for key in rc_state} == rc_state

    input_root = tmp_path / "harmony"
    input_root.mkdir()
    write_anndata(input_root, output)
    before = hashlib.sha256((input_root / "data.h5ad").read_bytes()).hexdigest()
    staging = tmp_path / "plot.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))
    records = harmony_convergence_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"adjusted_basis": "X_pca_harmony", "max_image_pixels": 40_000_000},
    )
    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert read_plot(staging / records[0]["payload"])[0] == plot.png
    assert hashlib.sha256((input_root / "data.h5ad").read_bytes()).hexdigest() == before


def test_harmony_convergence_plot_rejects_payload_and_coordinate_tampering(science, monkeypatch):
    _install_fake_harmony(monkeypatch, science)
    output, _, _ = run_anndata_operation(
        harmony_integration,
        _harmony_input(science),
        _harmony_parameters(),
    ).result

    changed_payload = output.copy()
    changed_payload.uns["harmony"]["X_pca_harmony"]["payload"]["harmony_objective"][1]["value"] = 9.0
    with pytest.raises(ValueError, match="payload fingerprint"):
        harmony_convergence_plot_owned(changed_payload)

    changed_coordinates = output.copy()
    changed_coordinates.obsm["X_pca_harmony"][0, 0] += 1.0
    with pytest.raises(ValueError, match="adjusted basis.*fingerprint"):
        harmony_convergence_plot_owned(changed_coordinates)

    changed_batches = output.copy()
    changed_batches.obs.loc[changed_batches.obs_names[0], "batch"] = "b"
    with pytest.raises(ValueError, match="Technical-batch labels.*fingerprint"):
        harmony_convergence_plot_owned(changed_batches)


def test_harmony_convergence_plot_enforces_the_pixel_budget(science, monkeypatch):
    _install_fake_harmony(monkeypatch, science)
    output, _, _ = run_anndata_operation(
        harmony_integration,
        _harmony_input(science),
        _harmony_parameters(),
    ).result

    with pytest.raises(ValueError, match="exceeding max_image_pixels=1"):
        harmony_convergence_plot_owned(output, max_image_pixels=1)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda result, np: result.objective_harmony.__setitem__(1, np.nan), "only finite"),
        (lambda result, _np: result.objective_kmeans.pop(), "every clustering update"),
        (lambda result, _np: result.objective_kmeans.__setitem__(2, 8.0), "outer-round boundaries"),
        (lambda result, _np: result.kmeans_rounds.__setitem__(0, 5), "clustering-iteration bound"),
    ],
)
def test_harmony_rejects_incomplete_or_inconsistent_objective_contract(
    science, monkeypatch, mutate, message
):
    _install_fake_harmony(monkeypatch, science)
    original = sys.modules["harmonypy"].run_harmony

    def invalid_result(*args, **kwargs):
        result = original(*args, **kwargs)
        mutate(result, science.np)
        return result

    sys.modules["harmonypy"].run_harmony = invalid_result
    with pytest.raises(RuntimeError, match=message):
        run_anndata_operation(harmony_integration, _harmony_input(science), _harmony_parameters())


def test_scvi_retains_complete_native_epoch_identity_in_training_diagnostics(science, monkeypatch):
    _install_fake_scvi(monkeypatch, science)

    _output, _model, report, _code = run_anndata_operation(
        scvi_integration,
        _harmony_input(science),
        _scvi_parameters(),
        artifact_output="model",
    ).result

    diagnostics = report.summary["key_results"]["training"]
    assert diagnostics["schema_version"] == "openbio-singlecell/scvi-training-diagnostics/v1"
    assert diagnostics["actual_epochs"] == 3
    assert diagnostics["metrics"]["elbo_train"]["records"] == [
        {"epoch": 0, "value": 12.0},
        {"epoch": 1, "value": 8.0},
        {"epoch": 2, "value": 6.0},
    ]
    assert diagnostics["metrics"]["elbo_validation"]["records"][-1] == {"epoch": 2, "value": 7.0}
    assert list(diagnostics["metrics"]) == [
        "elbo_train",
        "elbo_validation",
        "reconstruction_loss_train",
        "reconstruction_loss_validation",
        "kl_local_train",
        "kl_local_validation",
        "kl_global_train",
        "kl_global_validation",
    ]
    assert diagnostics["metrics"]["kl_global_train"]["records"][-1] == {"epoch": 2, "value": 0.1}
    assert diagnostics["training_observations"] == 6
    assert diagnostics["fitted_features"] == 2


def test_scvi_native_artifact_round_trips_diagnostics_and_binds_them_to_model_file(
    tmp_path, science, monkeypatch
):
    _install_fake_scvi(monkeypatch, science)
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, _harmony_input(science))
    before = hashlib.sha256((input_root / "data.h5ad").read_bytes()).hexdigest()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = scvi_integration(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        _scvi_parameters(),
    )

    model_record = next(record for record in records if record["name"] == "model")
    model_root = staging / model_record["payload"]
    payload = read_scvi_training_diagnostics(model_root)
    summary = next(record["value"]["summary"] for record in records if record["name"] == "summary")
    assert payload["training_diagnostics"] == summary["key_results"]["training"]
    assert payload["training_parameters"] == summary["parameters"]
    assert payload["software_versions"] == summary["software_versions"]
    assert payload["model_identity"]["path"] == "model.pt"
    assert payload["model_identity"]["size"] == len(b"native-model")
    assert len(payload["model_identity"]["sha256"]) == 64
    assert (model_root / "openbio-training-diagnostics.json").is_file()
    assert read_anndata(input_root).obsm.keys() == {"X_pca"}
    assert hashlib.sha256((input_root / "data.h5ad").read_bytes()).hexdigest() == before

    model_hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in model_root.iterdir()
        if path.is_file()
    }
    plot_staging = tmp_path / "plot.partial"
    plot_staging.mkdir()
    plot_context = OperationContext.from_request_path(plot_staging / "request.json", str(uuid.uuid4()))
    plot_records = scvi_training_plot(
        plot_context,
        {
            "model": {
                "type": "artifact",
                "kind": "OPENBIO_SCVI_MODEL",
                "codec": "scvi-native-directory",
                "path": str(model_root.resolve()),
            }
        },
        {"max_image_pixels": 40_000_000},
    )
    assert [record["name"] for record in plot_records] == ["plot", "summary", "code"]
    runtime_png = read_plot(plot_staging / plot_records[0]["payload"])[0]
    assert runtime_png.startswith(b"\x89PNG\r\n\x1a\n")
    namespace: dict[str, object] = {}
    generated_code = next(record["value"] for record in plot_records if record["name"] == "code")
    exec(generated_code, namespace)
    assert namespace["plot_scvi_training"](model_root) == runtime_png
    assert model_hashes == {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in model_root.iterdir()
        if path.is_file()
    }

    model_path = model_root / "model.pt"
    model_path.write_bytes(b"native-modem")
    with pytest.raises(ValueError, match="model identity"):
        read_scvi_training_diagnostics(model_root)
    model_path.write_bytes(b"native-model")
    sidecar = model_root / "openbio-training-diagnostics.json"
    sidecar.write_text(sidecar.read_text(encoding="utf-8").replace('"value":12.0', '"value":12.1', 1), encoding="utf-8")
    with pytest.raises(ValueError, match="payload fingerprint"):
        read_scvi_training_diagnostics(model_root)


@pytest.mark.parametrize("corruption", ["missing_metric", "step_index", "epoch_gap", "nonfinite"])
def test_scvi_rejects_incomplete_or_noncanonical_native_epoch_history(science, monkeypatch, corruption):
    module = _install_fake_scvi(monkeypatch, science)
    model = module.model.SCVI(_harmony_input(science), n_latent=2)
    if corruption == "missing_metric":
        del model.history["kl_local_train"]
    elif corruption == "step_index":
        model.history["elbo_train"].index.name = "step"
    elif corruption == "epoch_gap":
        model.history["elbo_train"].index = science.pd.Index([0, 2, 3], name="epoch")
    else:
        model.history["elbo_train"].iloc[1, 0] = science.np.nan

    with pytest.raises(RuntimeError, match="training history|training diagnostics"):
        build_scvi_training_diagnostics(model)


def test_scvi_training_plot_schema_consumes_the_native_model_artifact():
    schema = OpenBioSingleCellSCVITrainingPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellSCVITrainingPlot"
    assert schema.display_name == "scVI Training Plot"
    assert schema.category == "openbio/single-cell/batch-integration"
    assert [item.id for item in schema.inputs] == ["model", "max_image_pixels"]
    assert schema.inputs[0].io_type == "OPENBIO_SCVI_MODEL"
    assert schema.inputs[1].default == 40_000_000
    assert schema.inputs[1].advanced is True
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("plot", PlotResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_scvi_training_plot_renders_all_retained_metric_pairs_without_training(tmp_path, science):
    import matplotlib

    epoch = [0, 1, 2]
    values = {
        "elbo_train": [12.0, 8.0, 6.0],
        "elbo_validation": [13.0, 9.0, 7.0],
        "reconstruction_loss_train": [10.0, 7.0, 5.5],
        "reconstruction_loss_validation": [11.0, 8.0, 6.5],
        "kl_local_train": [2.0, 1.0, 0.5],
        "kl_local_validation": [2.2, 1.2, 0.6],
        "kl_global_train": [0.3, 0.2, 0.1],
        "kl_global_validation": [0.35, 0.25, 0.15],
    }
    diagnostics = {
        "schema_version": "openbio-singlecell/scvi-training-diagnostics/v1",
        "actual_epochs": 3,
        "metrics": {
            metric: {
                "records": [
                    {"epoch": epoch_value, "value": metric_value}
                    for epoch_value, metric_value in zip(epoch, metric_values, strict=True)
                ]
            }
            for metric, metric_values in values.items()
        },
        "train_cells": 5,
        "validation_cells": 1,
        "test_cells": 0,
        "training_observations": 6,
        "fitted_features": 2,
        "device": "cpu",
    }
    model_root = tmp_path / "native"
    model_root.mkdir()
    (model_root / "model.pt").write_bytes(b"native-model")
    native_payload = write_scvi_training_diagnostics(
        model_root,
        diagnostics,
        training_parameters={"technical_batch_key": "batch"},
        software_versions={"scvi-tools": "1.5.0.post1"},
    )
    numpy_state = science.np.random.get_state()
    rc_state = {key: matplotlib.rcParams[key] for key in ("figure.dpi", "font.size", "axes.grid")}

    plot, report, code = scvi_training_plot_owned(native_payload)

    assert plot.png.startswith(b"\x89PNG\r\n\x1a\n")
    details = report.summary["key_results"]
    assert details["actual_epochs"] == 3
    assert details["plotted_panels"] == ["elbo", "reconstruction_loss", "kl_local", "kl_global"]
    assert details["plotted_metrics"] == list(values)
    assert details["metric_summaries"]["elbo_train"] == {
        "points": 3,
        "first_epoch": 0,
        "last_epoch": 2,
        "first": 12.0,
        "last": 6.0,
        "best": 6.0,
    }
    assert report.summary["parameters"] == {"max_image_pixels": 40_000_000}
    assert "train(" not in code
    assert "SCVI.load" not in code
    namespace: dict[str, object] = {}
    exec(code, namespace)
    assert namespace["plot_scvi_training"](model_root) == plot.png
    restored_numpy_state = science.np.random.get_state()
    assert restored_numpy_state[0] == numpy_state[0]
    science.np.testing.assert_array_equal(restored_numpy_state[1], numpy_state[1])
    assert restored_numpy_state[2:] == numpy_state[2:]
    assert {key: matplotlib.rcParams[key] for key in rc_state} == rc_state
    assert diagnostics["metrics"]["elbo_train"]["records"][0] == {"epoch": 0, "value": 12.0}

    with pytest.raises(ValueError, match="exceeding max_image_pixels=1"):
        scvi_training_plot_owned(native_payload, max_image_pixels=1)
