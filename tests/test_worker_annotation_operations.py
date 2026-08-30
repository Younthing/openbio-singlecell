from __future__ import annotations

import hashlib
import uuid
from pathlib import Path

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.operations_annotation import map_cluster_annotations
from openbio_singlecell.operations_schist import schist_nested_model
from openbio_singlecell.operations_trajectory import diffusion_map
from openbio_singlecell.worker_protocol import OperationContext, ProtocolError, WorkerResponse


def _artifact_descriptor(root: Path, *, kind: str = "OPENBIO_ANNDATA", codec: str = "anndata-h5ad-v1"):
    return {"type": "artifact", "kind": kind, "codec": codec, "path": str(root.resolve())}


def test_map_cluster_annotations_publishes_a_new_anndata_without_rewriting_input(tmp_path: Path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    adata = science.ad.AnnData(
        X=science.np.eye(3),
        obs=science.pd.DataFrame(
            {"leiden": science.pd.Categorical(["0", "1", "0"], categories=["0", "1"])},
            index=["a", "b", "c"],
        ),
    )
    ensure_metadata(adata, source={"kind": "test"})
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).hexdigest()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = map_cluster_annotations(
        context,
        {"adata": _artifact_descriptor(input_root)},
        {
            "groupby": "leiden",
            "mapping_json": '{"0":"T cell","1":"B cell"}',
            "output_column": "cell_type",
            "unmapped_policy": "error",
            "annotation_status": "curated",
            "overwrite_existing": False,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["adata", "summary", "code"]
    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == before
    output = read_anndata(staging / records[0]["payload"])
    assert output.obs["cell_type"].astype(str).tolist() == ["T cell", "B cell", "T cell"]
    assert records[1]["value"]["summary"]["node_id"] == "OpenBioSingleCellMapClusterAnnotations"


def test_schist_operation_rejects_extended_parameter_objects(tmp_path: Path, science):
    input_root = tmp_path / "schist-input"
    input_root.mkdir()
    write_anndata(input_root, science.ad.AnnData(X=science.np.eye(2)))
    staging = tmp_path / "schist.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    with pytest.raises(ProtocolError, match="parameters must be exactly"):
        schist_nested_model(
            context,
            {"adata": _artifact_descriptor(input_root)},
            {
                "random_seed": 123,
                "neighbors_key": "neighbors",
                "key_added": "nsbm",
                "posterior_samples": 100,
                "degree_correction": True,
                "overwrite_existing": False,
                "max_working_memory_gib": 8.0,
                "unexpected": True,
            },
        )


def test_trajectory_operation_rejects_extended_parameter_objects(tmp_path: Path, science):
    input_root = tmp_path / "trajectory-input"
    input_root.mkdir()
    write_anndata(input_root, science.ad.AnnData(X=science.np.eye(2)))
    staging = tmp_path / "trajectory.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    with pytest.raises(ProtocolError, match="parameters must be exactly"):
        diffusion_map(
            context,
            {"adata": _artifact_descriptor(input_root)},
            {
                "neighbors_key": "neighbors",
                "n_comps": 3,
                "overwrite_existing": False,
                "random_seed": 0,
                "unexpected": True,
            },
        )


def test_annotation_schist_and_trajectory_operations_do_not_import_comfy_modules():
    package = Path(__file__).parents[1] / "openbio_singlecell"
    for filename in ("operations_annotation.py", "operations_schist.py", "operations_trajectory.py"):
        source = (package / filename).read_text(encoding="utf-8")
        assert "comfy_api" not in source
        assert "folder_paths" not in source
        assert "nodes_" not in source
