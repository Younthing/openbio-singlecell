from __future__ import annotations

import hashlib
import uuid
from pathlib import Path

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_table, write_anndata
from openbio_singlecell.operations_abundance import sample_composition_summary
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse


def _artifact_descriptor(root: Path) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": "OPENBIO_ANNDATA",
        "codec": "anndata-h5ad-v1",
        "path": str(root.resolve()),
    }


def test_sample_composition_operation_writes_portable_table_without_rewriting_input(tmp_path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    adata = science.ad.AnnData(
        X=science.np.arange(8, dtype=float).reshape(4, 2),
        obs=science.pd.DataFrame(
            {
                "sample": ["s1", "s1", "s2", "s2"],
                "condition": ["control", "control", "treated", "treated"],
                "cell_type": science.pd.Categorical(["A", "B", "A", "A"], categories=["A", "B"]),
            },
            index=["c1", "c2", "c3", "c4"],
        ),
        var=science.pd.DataFrame(index=["g1", "g2"]),
    )
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).hexdigest()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = sample_composition_summary(
        context,
        {"adata": _artifact_descriptor(input_root)},
        {
            "sample_key": "sample",
            "condition_key": "condition",
            "annotation_key": "cell_type",
            "annotation_status": "curated",
            "max_output_rows": 100,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["table", "summary", "code"]
    assert records[0]["codec"] == "table-jsonl-v1"
    table, metadata = read_table(staging / records[0]["payload"])
    assert table["cell_count"].tolist() == [1, 1, 2, 0]
    assert metadata["kind"] == "table"
    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == before


def test_abundance_and_population_operations_do_not_import_comfy_modules():
    package = Path(__file__).parents[1] / "openbio_singlecell"
    for filename in ("operations_abundance.py", "operations_population.py"):
        source = (package / filename).read_text(encoding="utf-8")
        assert "comfy_api" not in source
        assert "folder_paths" not in source
        assert "nodes_" not in source
