"""Run the optional real scVI smoke in an isolated native-runtime process."""

from __future__ import annotations

import json
import tempfile
import warnings
from importlib import import_module, metadata
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.artifact_envelope import summary_from_metadata
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.operations_input import ANNDATA_CODEC, ANNDATA_KIND
from openbio_singlecell.operations_integration import scvi_integration, scvi_training_plot_owned
from openbio_singlecell.scvi_model import read_scvi_training_diagnostics
from openbio_singlecell.worker_protocol import OperationContext


def _require(condition: object, message: str) -> None:
    if not bool(condition):
        raise AssertionError(message)


def main() -> None:
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    rng = np.random.default_rng(71)
    counts = rng.poisson(2.0, size=(32, 4)).astype(float)
    counts[:, 0] = 0.0
    counts[:, 1] += 1.0
    obs = pd.DataFrame(
        {"batch": pd.Categorical(["lane_a"] * 16 + ["lane_b"] * 16)},
        index=[f"real_cell_{index}" for index in range(32)],
    )
    var = pd.DataFrame(index=[f"real_gene_{index}" for index in range(4)])
    adata = ad.AnnData(counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    openbio_metadata = ensure_metadata(adata)
    openbio_metadata["analysis_history"] = {
        "snapshot": {
            "operation": "snapshot_expression",
            "parameters": {"source": "X", "layer_name": "counts"},
        }
    }
    adata.uns["openbio_singlecell"] = openbio_metadata

    with tempfile.TemporaryDirectory(prefix="openbio-real-scvi-") as directory:
        root = Path(directory)
        input_root = root / "input"
        input_root.mkdir()
        write_anndata(input_root, adata)
        input_payload = input_root / ANNDATA_PAYLOAD
        before = input_payload.read_bytes()
        staging = root / "run.partial"
        staging.mkdir()
        records = scvi_integration(
            OperationContext(staging, "00000000-0000-4000-8000-000000000003"),
            {
                "adata": {
                    "type": "artifact",
                    "path": str(input_root.resolve()),
                    "kind": ANNDATA_KIND,
                    "codec": ANNDATA_CODEC,
                }
            },
            {
                "source": {"source": "layer", "counts_layer": "counts"},
                "technical_batch_key": "batch",
                "categorical_covariates": "",
                "continuous_covariates": "",
                "n_latent": 8,
                "gene_likelihood": "zinb",
                "n_layers": 1,
                "dispersion": "gene",
                "dropout_rate": 0.1,
                "epochs": {"epochs": "fixed", "max_epochs": 1},
                "early_stopping": False,
                "train_size": 1.0,
                "batch_size": 16,
                "size_factor_key": "",
                "accelerator": "cpu",
                "output_key": "X_scVI",
                "overwrite_existing": False,
                "random_seed": 11,
            },
        )
        output = read_anndata(staging / records[0]["payload"])
        model_root = staging / records[1]["payload"]
        report = summary_from_metadata(records[2]["value"])
        code = records[3]["value"]
        model = import_module("scvi.model").SCVI.load(str(model_root), adata=output)
        native_diagnostics = read_scvi_training_diagnostics(model_root)
        training_plot, training_report, training_code = scvi_training_plot_owned(native_diagnostics)
        training_namespace = {}
        exec(training_code, training_namespace)
        generated_training_png = training_namespace["plot_scvi_training"](model_root)

        _require(input_payload.read_bytes() == before, "Real scVI operation modified its input artifact.")
        _require(model_root.is_dir(), "Real scVI operation did not publish its native model directory.")
        json.dumps(records, allow_nan=False)

    latent = output.obsm["X_scVI"]
    _require(latent.shape == (32, 8), "Real scVI smoke returned an unexpected latent shape.")
    _require(np.isfinite(latent).all(), "Real scVI smoke returned non-finite latent values.")
    _require(type(model).__name__.endswith("SCVI"), "Real scVI smoke returned the wrong model class.")
    _require(report.summary["key_results"]["zero_total_genes"] == 1, "Zero-total gene accounting drifted.")
    _require(report.summary["parameters"]["n_latent"] == 8, "n_latent reporting drifted.")
    _require(report.summary["parameters"]["train_size"] == 1.0, "train_size reporting drifted.")
    _require(
        "kl_global_train" in native_diagnostics["training_diagnostics"]["metrics"],
        "Real scVI training diagnostics dropped the native global-KL sequence.",
    )
    _require(
        report.summary["parameters"]["check_val_every_n_epoch"] is None,
        "No-holdout training cadence reporting drifted.",
    )
    overcomplete_warning = any(
        "overcomplete architecture" in warning_text for warning_text in report.summary["warnings"]
    )
    no_validation_holdout_warning = any(
        "leaves no validation holdout" in warning_text for warning_text in report.summary["warnings"]
    )
    _require(overcomplete_warning, "Real scVI smoke omitted the overcomplete-architecture warning.")
    _require(no_validation_holdout_warning, "Real scVI smoke omitted the no-holdout warning.")
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<scvi-real-code>", "exec")
    np.testing.assert_array_equal(adata.layers["counts"], counts)
    _require("X_scVI" not in adata.obsm, "Real scVI smoke mutated the caller-owned AnnData.")

    payload = {
        "backend_distribution": "scvi-tools",
        "backend_version": metadata.version("scvi-tools"),
        "code_compiles": True,
        "input_counts_unchanged": True,
        "input_has_no_output_embedding": True,
        "latent_shape": list(latent.shape),
        "model_class_endswith_scvi": True,
        "no_validation_holdout_warning": no_validation_holdout_warning,
        "overcomplete_warning": overcomplete_warning,
        "status": "pass",
        "training_diagnostic_metrics": training_report.summary["key_results"]["plotted_metrics"],
        "training_plot_code_parity": generated_training_png == training_plot.png,
        "training_plot_png": training_plot.png.startswith(b"\x89PNG\r\n\x1a\n"),
        "zero_total_genes": report.summary["key_results"]["zero_total_genes"],
    }
    print("OPENBIO_SCVI_REAL_SMOKE=" + json.dumps(payload, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
