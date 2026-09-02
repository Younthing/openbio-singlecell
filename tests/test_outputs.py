from __future__ import annotations

import asyncio
import csv
import sys
from io import BytesIO
from pathlib import Path

import pytest
from PIL import Image

from openbio_singlecell.artifact_codecs import write_anndata, write_plot, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.artifact_runtime import ArtifactTicket
from openbio_singlecell.artifact_service import current_artifact_runtime, initialize_artifact_service
from openbio_singlecell.contracts import PlotResult, SummaryResult, TableResult
from openbio_singlecell.node_types import PlotResultType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_output import (
    OpenBioSingleCellExportCSV,
    OpenBioSingleCellPreviewResult,
    OpenBioSingleCellSaveH5AD,
    OpenBioSingleCellSavePNG,
    export_csv,
    preview_result,
    save_h5ad,
    save_png,
)


@pytest.fixture(scope="module", autouse=True)
def _artifact_service(tmp_path_factory):
    asyncio.run(initialize_artifact_service(tmp_path_factory.mktemp("output-runtime"), sys.executable))


def _fields():
    return dict(
        title="result",
        parameters={},
        description="description",
        warnings=[],
        input_cells=2,
        input_genes=2,
        random_seed=0,
        elapsed_seconds=0.1,
        source={},
    )


def _summary():
    return SummaryResult(
        summary={
            "schema_version": 1,
            "node_id": "OpenBioSingleCellTest",
            "methods": "Test method.",
            "results": "Test result.",
            "key_results": {},
            "parameters": {},
            "warnings": [],
            "limitations": [],
            "references": [{"citation": "Test.", "url": "https://example.org", "kind": "method"}],
            "software_versions": {"python": "test", "openbio-singlecell": "test"},
        },
        **_fields(),
    )


def _png(color: str = "red") -> bytes:
    buffer = BytesIO()
    Image.new("RGB", (2, 2), color=color).save(buffer, format="PNG")
    return buffer.getvalue()


def _publish(runtime, name: str, kind: str, codec: str, writer) -> ArtifactTicket:
    lease = runtime.begin_run()
    root = lease.staging_path / name
    root.mkdir()
    writer(root)
    return lease.publish({name: {"kind": kind, "codec": codec, "payload": name}})[name]


def _anndata_ticket(runtime, science) -> ArtifactTicket:
    adata = science.ad.AnnData(science.sparse.csr_matrix(science.np.eye(2)))
    adata.obs_names = ["cell-a", "cell-b"]
    adata.var_names = ["gene-a", "gene-b"]
    return _publish(runtime, "adata", "OPENBIO_ANNDATA", "anndata-h5ad-v1", lambda root: write_anndata(root, adata))


def _table_ticket(runtime, science, rows: int = 2) -> ArtifactTicket:
    frame = science.pd.DataFrame(
        {
            "gene": [f"G{index}" for index in range(rows)],
            "effect": science.np.arange(rows, dtype=float),
        },
        index=science.pd.Index([f"cell-{index}" for index in range(rows)], name="cell_id"),
    )
    result = TableResult(table=frame, **_fields())
    return _publish(
        runtime,
        "table",
        "OPENBIO_SINGLE_CELL_TABLE",
        "table-jsonl-v1",
        lambda root: write_table(root, frame, result_metadata(result)),
    )


def _plot_ticket(runtime, color: str = "red") -> ArtifactTicket:
    png = _png(color)
    result = PlotResult(png=png, **_fields())
    return _publish(
        runtime,
        "plot",
        "OPENBIO_SINGLE_CELL_PLOT",
        "plot-png-v1",
        lambda root: write_plot(root, png, result_metadata(result)),
    )


def test_output_schemas_accept_tickets_and_drop_h5ad_recompression():
    h5ad = OpenBioSingleCellSaveH5AD.GET_SCHEMA()
    preview = OpenBioSingleCellPreviewResult.GET_SCHEMA()
    csv_schema = OpenBioSingleCellExportCSV.GET_SCHEMA()
    png = OpenBioSingleCellSavePNG.GET_SCHEMA()

    assert [item.id for item in h5ad.inputs] == ["adata", "filename_prefix", "overwrite"]
    assert preview.inputs[0].get_io_type() == ",".join(
        [SummaryResultType.io_type, TableResultType.io_type, PlotResultType.io_type]
    )
    assert csv_schema.inputs[0].io_type == TableResultType.io_type
    assert png.inputs[0].io_type == PlotResultType.io_type
    assert all(schema.description for schema in (preview, h5ad, csv_schema, png))


def test_summary_preview_stays_in_memory_and_creates_no_file(comfy_directories):
    _, output_dir, temp_dir = comfy_directories
    ui = preview_result(_summary(), "node")

    assert ui["openbio_singlecell"][0]["kind"] == "summary"
    assert not list(output_dir.rglob("*"))
    assert not list(temp_dir.rglob("*"))


def test_table_preview_streams_only_bounded_rows(comfy_directories, science):
    _, output_dir, temp_dir = comfy_directories
    ticket = _table_ticket(current_artifact_runtime(), science, rows=105)

    payload = preview_result(ticket, "table-node")["openbio_singlecell"][0]

    assert payload["kind"] == "table"
    assert payload["total_rows"] == 105
    assert len(payload["rows"]) == 100
    assert payload["columns"] == ["gene", "effect"]
    assert payload["index_name"] == "cell_id"
    assert payload["index_values"] == [f"cell-{index}" for index in range(100)]
    assert any("Preview limited" in warning for warning in payload["warnings"])
    assert not list(output_dir.rglob("*"))
    assert not list(temp_dir.rglob("*"))


def test_plot_preview_copies_ticket_file_and_overwrites_per_node(comfy_directories):
    _, _, temp_dir = comfy_directories
    runtime = current_artifact_runtime()
    first = preview_result(_plot_ticket(runtime, "red"), "same-node")
    second_ticket = _plot_ticket(runtime, "blue")
    second = preview_result(second_ticket, "same-node")

    files = list((temp_dir / "openbio-singlecell").iterdir())
    assert len(files) == 1
    assert files[0].read_bytes() == (runtime.resolve(second_ticket) / "plot.png").read_bytes()
    assert first["images"][0]["filename"] == second["images"][0]["filename"]
    assert first["images"][0]["type"] == "temp"


def test_file_outputs_stream_from_tickets_without_changing_sources(comfy_directories, science):
    runtime = current_artifact_runtime()
    adata = _anndata_ticket(runtime, science)
    table = _table_ticket(runtime, science)
    plot = _plot_ticket(runtime, "orange")
    source_h5ad = runtime.resolve(adata) / "data.h5ad"
    source_png = runtime.resolve(plot) / "plot.png"
    before_h5ad = source_h5ad.read_bytes()
    before_png = source_png.read_bytes()

    h5ad_first = save_h5ad(adata, "adata")
    h5ad_second = save_h5ad(adata, "adata")
    csv_target = export_csv(table, "markers", overwrite=True)
    png_target = save_png(plot, "umap", overwrite=True)

    assert h5ad_first.filename == "adata_00001.h5ad"
    assert h5ad_second.filename == "adata_00002.h5ad"
    assert Path(h5ad_first.path).read_bytes() == before_h5ad
    assert Path(png_target.path).read_bytes() == before_png
    with Path(csv_target.path).open(encoding="utf-8", newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["cell_id", "gene", "effect"],
            ["cell-0", "G0", "0.0"],
            ["cell-1", "G1", "1.0"],
        ]
    assert source_h5ad.read_bytes() == before_h5ad
    assert source_png.read_bytes() == before_png


def test_csv_export_rejects_index_header_that_collides_with_a_data_column(
    comfy_directories,
    science,
):
    frame = science.pd.DataFrame(
        {"cell_id": ["value"]},
        index=science.pd.Index(["cell-0"], name="cell_id"),
    )
    result = TableResult(table=frame, **_fields())
    ticket = _publish(
        current_artifact_runtime(),
        "table",
        "OPENBIO_SINGLE_CELL_TABLE",
        "table-jsonl-v1",
        lambda root: write_table(root, frame, result_metadata(result)),
    )

    with pytest.raises(ValueError, match="collide"):
        export_csv(ticket, "duplicate-header")


def test_output_adapters_reject_wrong_ticket_kinds(comfy_directories, science):
    runtime = current_artifact_runtime()
    table = _table_ticket(runtime, science)
    plot = _plot_ticket(runtime)

    for operation, value, expected in (
        (save_h5ad, table, "OPENBIO_ANNDATA"),
        (export_csv, plot, "OPENBIO_SINGLE_CELL_TABLE"),
        (save_png, table, "OPENBIO_SINGLE_CELL_PLOT"),
    ):
        try:
            operation(value, "wrong")
        except (TypeError, ValueError) as error:
            assert expected in str(error)
        else:
            raise AssertionError("wrong artifact kind was accepted")


def test_manifest_mismatch_preserves_existing_output(comfy_directories, science):
    runtime = current_artifact_runtime()
    ticket = _anndata_ticket(runtime, science)
    target = save_h5ad(ticket, "stable", overwrite=True)
    previous = Path(target.path).read_bytes()
    (runtime.resolve(ticket) / "data.h5ad").write_bytes(b"changed")

    try:
        save_h5ad(ticket, "stable", overwrite=True)
    except ValueError as error:
        assert "payload size" in str(error)
    else:
        raise AssertionError("mutated artifact was accepted")

    assert Path(target.path).read_bytes() == previous
