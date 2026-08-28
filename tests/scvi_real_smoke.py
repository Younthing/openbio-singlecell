"""Run the optional real scVI smoke in an isolated native-runtime process."""

from __future__ import annotations

import json
import warnings
from importlib import metadata

import anndata as ad
import numpy as np
import pandas as pd

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.nodes_integration import OpenBioSingleCellSCVIIntegration


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

    output, model, report, code = OpenBioSingleCellSCVIIntegration.execute(
        adata,
        source={"source": "layer", "counts_layer": "counts"},
        n_latent=8,
        epochs={"epochs": "fixed", "max_epochs": 1},
        train_size=1.0,
        batch_size=16,
        accelerator="cpu",
        random_seed=11,
    ).result

    latent = output.obsm["X_scVI"]
    _require(latent.shape == (32, 8), "Real scVI smoke returned an unexpected latent shape.")
    _require(np.isfinite(latent).all(), "Real scVI smoke returned non-finite latent values.")
    _require(model.evidence["model_class"].endswith("SCVI"), "Real scVI smoke returned the wrong model class.")
    _require(report.summary["key_results"]["zero_total_genes"] == 1, "Zero-total gene accounting drifted.")
    _require(report.summary["parameters"]["n_latent"] == 8, "n_latent reporting drifted.")
    _require(report.summary["parameters"]["train_size"] == 1.0, "train_size reporting drifted.")
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
        "zero_total_genes": report.summary["key_results"]["zero_total_genes"],
    }
    print("OPENBIO_SCVI_REAL_SMOKE=" + json.dumps(payload, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
