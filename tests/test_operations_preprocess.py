from __future__ import annotations

import hashlib
import json
import uuid
from pathlib import Path

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.operations_preprocess import (
    highly_variable_genes,
    normalize_to_layer,
    normalize_total,
    pearson_residuals_to_layer,
    scale,
)
from openbio_singlecell.worker_protocol import OperationContext


def _artifact_descriptor(root: Path) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": "OPENBIO_ANNDATA",
        "codec": "anndata-h5ad-v1",
        "path": str(root.resolve()),
    }


def _run_operation(tmp_path: Path, name: str, operation, adata, parameters):
    input_root = tmp_path / f"{name}-input"
    input_root.mkdir()
    ensure_metadata(adata, source={"kind": "test"})
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / f"{name}.partial"
    staging.mkdir()
    records = operation(
        OperationContext(staging, str(uuid.uuid4())),
        {"adata": _artifact_descriptor(input_root)},
        parameters,
    )
    assert hashlib.sha256(input_path.read_bytes()).digest() == before
    assert [record["name"] for record in records] == ["adata", "summary", "code"]
    assert records[0]["kind"] == "OPENBIO_ANNDATA"
    assert records[0]["codec"] == "anndata-h5ad-v1"
    assert records[1]["type"] == "summary"
    assert records[2]["type"] == "string"
    compile(records[2]["value"], f"<{name}-code>", "exec")
    json.dumps(records, allow_nan=False)
    return read_anndata(staging / records[0]["payload"]), records


def _adata(science):
    counts = science.sparse.csr_matrix(
        science.np.asarray(
            [
                [8.0, 2.0, 1.0, 4.0],
                [3.0, 9.0, 2.0, 1.0],
                [4.0, 1.0, 10.0, 3.0],
                [2.0, 5.0, 3.0, 8.0],
            ]
        )
    )
    adata = science.ad.AnnData(X=counts.copy())
    adata.obs_names = [f"cell-{index}" for index in range(4)]
    adata.var_names = [f"gene-{index}" for index in range(4)]
    adata.layers["counts"] = counts.copy()
    adata.layers["log1p_norm"] = counts.log1p()
    adata.raw = adata
    return adata


def test_normalization_operations_publish_new_artifacts_without_mutating_input(tmp_path, science):
    normalized, records = _run_operation(
        tmp_path,
        "normalize-total",
        normalize_total,
        _adata(science),
        {"target_sum": 100.0, "source": {"source": "X"}},
    )
    science.np.testing.assert_allclose(science.np.asarray(normalized.X.sum(axis=1)).ravel(), 100.0)
    assert records[1]["value"]["summary"]["node_id"] == "OpenBioSingleCellNormalizeTotal"

    layered, _ = _run_operation(
        tmp_path,
        "normalize-layer",
        normalize_to_layer,
        _adata(science),
        {
            "source": {"source": "layer", "source_layer": "counts"},
            "target_sum": 100.0,
            "transform": "log1p",
            "output_layer": "normalized",
            "overwrite_existing": False,
        },
    )
    science.np.testing.assert_array_equal(layered.X.toarray(), _adata(science).X.toarray())
    assert "normalized" in layered.layers


def test_feature_modeling_operations_use_only_matrix_workspaces(tmp_path, science):
    pearson, _ = _run_operation(
        tmp_path,
        "pearson",
        pearson_residuals_to_layer,
        _adata(science),
        {
            "source": {"source": "layer", "source_layer": "counts"},
            "theta": 100.0,
            "clipping_mode": "sqrt_n_obs",
            "custom_clip": 10.0,
            "output_layer": "residuals",
            "overwrite_existing": False,
            "max_dense_gib": 2.0,
        },
    )
    assert pearson.layers["residuals"].shape == pearson.shape

    scaled, _ = _run_operation(
        tmp_path,
        "scale",
        scale,
        _adata(science),
        {
            "source": {"source": "layer", "layer_name": "log1p_norm"},
            "zero_center": True,
            "clipping_mode": "custom",
            "custom_max_value": 10.0,
            "output_layer": "scaled",
            "overwrite_existing": False,
            "max_dense_gib": 2.0,
        },
    )
    assert scaled.layers["scaled"].shape == scaled.shape


def test_hvg_subsetting_materializes_only_the_scientific_subset(tmp_path, science):
    selected, records = _run_operation(
        tmp_path,
        "hvg",
        highly_variable_genes,
        _adata(science),
        {
            "n_top_genes": 2,
            "flavor": "seurat",
            "source": {"source": "layer", "layer_name": "log1p_norm"},
            "batch_key": "",
            "always_keep_genes": "gene-3",
            "subset": True,
            "theta": 100.0,
            "clipping_mode": "sqrt_n_obs",
            "custom_clip": 10.0,
            "chunksize": 1000,
            "span": 0.3,
            "n_bins": 20,
            "overwrite_existing": False,
        },
    )
    assert selected.n_vars >= 2
    assert "gene-3" in selected.var_names
    assert records[1]["value"]["summary"]["key_results"]["subset_applied"] is True


def test_preprocess_node_modules_are_schema_only_and_log1p_example_is_in_place():
    package = Path(__file__).parents[1] / "openbio_singlecell"
    node_source = (package / "nodes_preprocess.py").read_text(encoding="utf-8")
    assert "def execute" not in node_source
    assert "dependencies" not in node_source

    source = (package / "operations_preprocess.py").read_text(
        encoding="utf-8"
    )
    assert "output = adata.copy()" not in source
    assert "sc.pp.log1p(adata, base=None, copy=False)" in source


def test_preprocess_operation_rejects_nonfinite_before_creating_output(tmp_path, science):
    adata = _adata(science)
    adata.layers["counts"][0, 0] = science.np.nan
    input_root = tmp_path / "invalid-input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    staging = tmp_path / "invalid.partial"
    staging.mkdir()
    with pytest.raises(ValueError, match="non-finite"):
        normalize_to_layer(
            OperationContext(staging, str(uuid.uuid4())),
            {"adata": _artifact_descriptor(input_root)},
            {
                "source": {"source": "layer", "source_layer": "counts"},
                "target_sum": 100.0,
                "transform": "none",
                "output_layer": "normalized",
                "overwrite_existing": False,
            },
        )
    assert not (staging / "outputs").exists()
