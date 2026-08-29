from __future__ import annotations

import ast
import asyncio
import hashlib
import inspect
import json
import sys
from pathlib import Path

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.operations_preprocess import log1p
from openbio_singlecell.worker_client import WorkerOperationError, probe_python, run_worker

ANNDATA_KIND = "OPENBIO_ANNDATA"
ANNDATA_CODEC = "anndata-h5ad-v1"


@pytest.fixture(scope="module")
def worker_identity():
    return asyncio.run(probe_python(sys.executable))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _file_descriptor(path: Path, portable_path: str) -> dict[str, object]:
    stat = path.stat()
    return {
        "type": "file",
        "path": str(path.resolve()),
        "provenance": {"path": portable_path, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns},
    }


def _artifact_descriptor(root: Path) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": ANNDATA_KIND,
        "codec": ANNDATA_CODEC,
        "path": str(root.resolve()),
    }


def _assert_uncompressed_h5ad(path: Path) -> None:
    h5py = pytest.importorskip("h5py")
    compressions = []
    with h5py.File(path, "r") as handle:
        handle.visititems(
            lambda _name, value: compressions.append(value.compression)
            if isinstance(value, h5py.Dataset)
            else None
        )
    assert compressions
    assert set(compressions) == {None}


def test_load_then_log1p_runs_in_real_one_shot_workers_without_mutating_inputs(
    tmp_path: Path,
    science,
    worker_identity,
):
    source = tmp_path / "sample.h5ad"
    matrix = science.sparse.csr_matrix(
        science.np.asarray(
            [
                [1.0, 2.0, 0.0],
                [4.0, 0.0, 2.0],
            ]
        )
    )
    with pytest.warns(UserWarning, match="Variable names are not unique"):
        original = science.ad.AnnData(
            X=matrix,
            obs=science.pd.DataFrame({"batch": ["a", "b"]}, index=["cell_a", "cell_b"]),
            var=science.pd.DataFrame(index=["gene", "gene", "other"]),
        )
    original.layers["counts"] = matrix.copy()
    original.raw = original
    ensure_metadata(original, display_name="embedded", source={"kind": "external", "path": "prior/source.h5ad"})
    original.write_h5ad(source, compression="gzip")
    source_before = _sha256(source)

    load_staging = tmp_path / "load.partial"
    load_staging.mkdir()
    load_records = asyncio.run(
        run_worker(
            worker_identity,
            load_staging / "request.json",
            "openbio.node.loadh5ad",
            inputs={"path": _file_descriptor(source, "sample.h5ad")},
            parameters={"make_var_names_unique": True},
        )
    )

    assert load_records == [
        {
            "type": "artifact",
            "name": "adata",
            "kind": ANNDATA_KIND,
            "codec": ANNDATA_CODEC,
            "payload": "outputs/adata",
        }
    ]
    json.dumps(load_records, allow_nan=False)
    assert _sha256(source) == source_before
    loaded_root = load_staging / load_records[0]["payload"]
    loaded_path = loaded_root / ANNDATA_PAYLOAD
    assert loaded_path.is_file()
    assert loaded_path.resolve() != source.resolve()
    _assert_uncompressed_h5ad(loaded_path)
    loaded = read_anndata(loaded_root)
    assert loaded.shape == original.shape
    assert loaded.var_names.is_unique
    assert loaded.uns["openbio_singlecell"]["display_name"] == "sample"
    loaded_source = loaded.uns["openbio_singlecell"]["source"]
    assert loaded_source["kind"] == "h5ad"
    assert loaded_source["path"] == "sample.h5ad"
    assert loaded_source["embedded_source"] == {"kind": "external", "path": "prior/source.h5ad"}
    assert loaded_source["axis_names"]["var_names_repaired"] == 1

    loaded_before = _sha256(loaded_path)
    input_values = loaded.X.toarray()
    input_counts = loaded.layers["counts"].toarray()
    input_raw = loaded.raw.X.toarray()
    log_staging = tmp_path / "log.partial"
    log_staging.mkdir()
    log_records = asyncio.run(
        run_worker(
            worker_identity,
            log_staging / "request.json",
            "openbio.node.log1p",
            inputs={"adata": _artifact_descriptor(loaded_root)},
            parameters={},
        )
    )

    assert [record["name"] for record in log_records] == ["adata", "summary", "code"]
    assert log_records[0] == {
        "type": "artifact",
        "name": "adata",
        "kind": ANNDATA_KIND,
        "codec": ANNDATA_CODEC,
        "payload": "outputs/adata",
    }
    assert log_records[1]["type"] == "summary"
    summary_metadata = log_records[1]["value"]
    assert summary_metadata["kind"] == "summary"
    assert summary_metadata["summary"]["node_id"] == "OpenBioSingleCellLog1p"
    assert summary_metadata["summary"]["key_results"]["base"] == "natural"
    assert "private AnnData" in summary_metadata["summary"]["methods"]
    assert log_records[2]["type"] == "string"
    assert log_records[2]["value"].endswith("\n")
    compile(log_records[2]["value"], "<OpenBioSingleCellLog1p-code>", "exec")
    json.dumps(log_records, allow_nan=False)

    assert _sha256(loaded_path) == loaded_before
    logged_root = log_staging / log_records[0]["payload"]
    logged_path = logged_root / ANNDATA_PAYLOAD
    assert logged_path.is_file()
    assert logged_path.resolve() != loaded_path.resolve()
    _assert_uncompressed_h5ad(logged_path)
    logged = read_anndata(logged_root)
    science.np.testing.assert_allclose(logged.X.toarray(), science.np.log1p(input_values))
    science.np.testing.assert_array_equal(logged.layers["counts"].toarray(), input_counts)
    science.np.testing.assert_array_equal(logged.raw.X.toarray(), input_raw)
    history = logged.uns["openbio_singlecell"]["analysis_history"]
    assert history[sorted(history)[-1]]["operation"] == "log1p"
    assert "log1p" in logged.uns

    assert "adata.copy(" not in inspect.getsource(log1p)


def test_log1p_nonfinite_failure_does_not_publish_an_output(
    tmp_path: Path,
    science,
    worker_identity,
):
    artifact_root = tmp_path / "nonfinite-input"
    artifact_root.mkdir()
    adata = science.ad.AnnData(X=science.np.asarray([[1.0, science.np.nan]]))
    ensure_metadata(adata, source={"kind": "test"})
    write_anndata(artifact_root, adata)
    input_path = artifact_root / ANNDATA_PAYLOAD
    input_before = _sha256(input_path)
    staging = tmp_path / "failed.partial"
    staging.mkdir()

    with pytest.raises(WorkerOperationError, match="non-finite"):
        asyncio.run(
            run_worker(
                worker_identity,
                staging / "request.json",
                "openbio.node.log1p",
                inputs={"adata": _artifact_descriptor(artifact_root)},
                parameters={},
            )
        )

    assert _sha256(input_path) == input_before
    assert not (staging / "outputs").exists()
    assert not (staging / "request.json").exists()
    assert not (staging / "response.json").exists()


def test_science_operations_reject_mixed_or_extended_input_descriptors(
    tmp_path: Path,
    science,
    worker_identity,
):
    source = tmp_path / "sample.h5ad"
    science.ad.AnnData(X=science.np.eye(2)).write_h5ad(source)
    descriptor = _file_descriptor(source, "sample.h5ad")
    descriptor["kind"] = ANNDATA_KIND
    staging = tmp_path / "invalid.partial"
    staging.mkdir()

    with pytest.raises(WorkerOperationError, match="descriptor fields"):
        asyncio.run(
            run_worker(
                worker_identity,
                staging / "request.json",
                "openbio.node.loadh5ad",
                inputs={"path": descriptor},
                parameters={"make_var_names_unique": True},
            )
        )

    assert not (staging / "outputs").exists()


def test_worker_science_modules_do_not_import_comfy_or_node_modules():
    package = Path(__file__).parents[1] / "openbio_singlecell"
    paths = [package / "worker_operations.py", *sorted(package.glob("operations_*.py"))]
    for path in paths:
        imports = {
            name
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for name in ([node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            + [alias.name for alias in node.names]
        }
        forbidden = {
            name
            for name in imports
            if name in {"comfy", "comfy_api", "folder_paths"}
            or name.startswith(("comfy.", "comfy_api.", "folder_paths.", "nodes_", "openbio_singlecell.nodes_"))
        }
        assert forbidden == set(), f"{path.name} imports {sorted(forbidden)}"
