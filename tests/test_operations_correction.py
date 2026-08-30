from __future__ import annotations

import hashlib
import json
import uuid
from pathlib import Path

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.operations_correction import filter_doublets, mark_mad_outliers, scrublet
from openbio_singlecell.worker_protocol import OperationContext


def _descriptor(root: Path) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": "OPENBIO_ANNDATA",
        "codec": "anndata-h5ad-v1",
        "path": str(root.resolve()),
    }


def _run(tmp_path: Path, name: str, operation, adata, parameters):
    input_root = tmp_path / f"{name}-input"
    input_root.mkdir()
    ensure_metadata(adata, source={"kind": "test"})
    write_anndata(input_root, adata)
    payload = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(payload.read_bytes()).digest()
    staging = tmp_path / f"{name}.partial"
    staging.mkdir()
    records = operation(
        OperationContext(staging, str(uuid.uuid4())),
        {"adata": _descriptor(input_root)},
        parameters,
    )
    assert hashlib.sha256(payload.read_bytes()).digest() == before
    assert [record["name"] for record in records] == ["adata", "summary", "code"]
    json.dumps(records, allow_nan=False)
    compile(records[2]["value"], f"<{name}-code>", "exec")
    return read_anndata(staging / records[0]["payload"]), records


def test_mad_marking_annotates_worker_private_anndata_without_full_copy(tmp_path, science):
    adata = science.ad.AnnData(X=science.np.ones((6, 3)))
    adata.obs["sample"] = ["a"] * 3 + ["b"] * 3
    adata.obs["total_counts"] = [10.0, 11.0, 100.0, 20.0, 21.0, 200.0]
    marked, records = _run(
        tmp_path,
        "mad",
        mark_mad_outliers,
        adata,
        {
            "metrics": "total_counts",
            "batch_key": "sample",
            "nmads": 2.0,
            "direction": "upper",
            "output_column": "outlier",
            "scale_mad": False,
            "minimum_group_size": 3,
        },
    )
    assert marked.obs["outlier"].tolist() == [False, False, True, False, False, True]
    assert records[1]["value"]["summary"]["node_id"] == "OpenBioSingleCellMarkMADOutliers"


def test_scrublet_uses_algorithm_workspace_and_writes_scores_back(tmp_path, science):
    rng = science.np.random.default_rng(17)
    matrix = science.sparse.csr_matrix(rng.poisson(2.0, size=(80, 30)).astype(float))
    adata = science.ad.AnnData(X=matrix)
    adata.obs_names = [f"cell-{index}" for index in range(80)]
    adata.var_names = [f"gene-{index}" for index in range(30)]
    scored, records = _run(
        tmp_path,
        "scrublet",
        scrublet,
        adata,
        {
            "batch_key": "",
            "random_seed": 123,
            "expected_doublet_rate": 0.05,
            "threshold_mode": "manual",
            "threshold": 0.25,
            "sim_doublet_ratio": 2.0,
            "n_prin_comps": 5,
            "n_neighbors": 5,
            "source": {"source": "X"},
        },
    )
    assert scored.shape == adata.shape
    assert scored.obs["doublet_score"].notna().all()
    assert scored.obs["predicted_doublet"].dtype == bool
    assert "scrublet" in scored.uns
    assert records[1]["value"]["summary"]["node_id"] == "OpenBioSingleCellScrublet"


def test_filter_doublets_materializes_the_scientific_subset(tmp_path, science):
    adata = science.ad.AnnData(X=science.np.arange(12, dtype=float).reshape(4, 3))
    adata.obs_names = ["a", "b", "c", "d"]
    adata.obs["predicted_doublet"] = [False, True, False, True]
    filtered, _ = _run(
        tmp_path,
        "filter",
        filter_doublets,
        adata,
        {"prediction_column": "predicted_doublet"},
    )
    assert filtered.obs_names.tolist() == ["a", "c"]
    assert filtered.n_vars == adata.n_vars


def test_correction_nodes_are_schema_only_and_operations_have_no_defensive_anndata_copy():
    package = Path(__file__).parents[1] / "openbio_singlecell"
    node_source = (package / "nodes_correction.py").read_text(encoding="utf-8")
    operation_source = (package / "operations_correction.py").read_text(encoding="utf-8")
    assert "def execute" not in node_source
    assert "dependencies" not in node_source
    assert "adata.copy()" not in operation_source
    assert "adata[keep].copy()" in operation_source
    assert "matrix.copy()" in operation_source
