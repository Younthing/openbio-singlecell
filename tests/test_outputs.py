from __future__ import annotations

from io import BytesIO
from pathlib import Path

import h5py
import pytest
from comfy_api.latest import io
from PIL import Image

from openbio_singlecell.contracts import PlotResult, SummaryResult, TableResult, ensure_metadata, record_history
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

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def result_fields():
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


def make_summary(summary=None):
    return SummaryResult(summary={"description": "summary"} if summary is None else summary, **result_fields())


def make_table(table):
    return TableResult(table=table, **result_fields())


def png_bytes(color="red"):
    buffer = BytesIO()
    Image.new("RGB", (2, 2), color=color).save(buffer, format="PNG")
    return buffer.getvalue()


def make_plot(png=None):
    return PlotResult(png=png_bytes() if png is None else png, **result_fields())


def test_save_h5ad_schema_uses_anndata_contract():
    schema = OpenBioSingleCellSaveH5AD.GET_SCHEMA()

    assert schema.inputs[0].id == "adata"
    assert schema.inputs[1].id == "filename_prefix"
    assert schema.inputs[1].default == "adata"
    assert schema.inputs[3].id == "compression"
    assert schema.inputs[3].default == "gzip"
    assert schema.inputs[3].options == ["gzip", "lzf", "none"]


def test_result_output_schemas_use_concrete_contracts():
    preview_schema = OpenBioSingleCellPreviewResult.GET_SCHEMA()
    csv_schema = OpenBioSingleCellExportCSV.GET_SCHEMA()
    png_schema = OpenBioSingleCellSavePNG.GET_SCHEMA()

    assert preview_schema.inputs[0].id == "result"
    assert preview_schema.inputs[0].get_io_type() == ",".join(
        [SummaryResultType.io_type, TableResultType.io_type, PlotResultType.io_type]
    )
    assert csv_schema.inputs[0].id == "table"
    assert csv_schema.inputs[0].io_type == TableResultType.io_type
    assert png_schema.inputs[0].id == "plot"
    assert png_schema.inputs[0].io_type == PlotResultType.io_type


def test_preview_does_not_persist_tables(comfy_directories):
    _, output_dir, temp_dir = comfy_directories
    ui = preview_result(make_summary(), "node")
    assert "openbio_singlecell" in ui
    assert not list(output_dir.rglob("*"))
    assert not list(temp_dir.rglob("*"))


def test_plot_preview_overwrites_per_node(comfy_directories):
    _, _, temp_dir = comfy_directories
    first_png = png_bytes("red")
    second_png = png_bytes("blue")
    first = preview_result(make_plot(first_png), "same-node")
    second = preview_result(make_plot(second_png), "same-node")
    files = list((temp_dir / "openbio-singlecell").iterdir())
    assert len(files) == 1
    assert files[0].read_bytes() == second_png
    assert first["images"][0]["filename"] == second["images"][0]["filename"]
    assert first["images"][0]["type"] == "temp"


def test_plot_preview_uses_v3_hidden_unique_id(comfy_directories):
    _, _, temp_dir = comfy_directories
    preview_node = OpenBioSingleCellPreviewResult.PREPARE_CLASS_CLONE(
        {"hidden_inputs": {io.Hidden.unique_id: "workflow-node-42"}}
    )
    png = png_bytes("green")
    output = preview_node.execute(make_plot(png))
    expected = preview_result(make_plot(png), "workflow-node-42")

    assert output.ui["images"][0]["filename"] == expected["images"][0]["filename"]
    assert len(list((temp_dir / "openbio-singlecell").iterdir())) == 1


def test_plot_preview_isolated_by_workflow_and_overwrites_within_workflow(comfy_directories):
    _, _, temp_dir = comfy_directories

    def execute(workflow_id, png, workflow_state=None):
        preview_node = OpenBioSingleCellPreviewResult.PREPARE_CLASS_CLONE(
            {
                "hidden_inputs": {
                    io.Hidden.unique_id: "shared-node",
                    io.Hidden.extra_pnginfo: {
                        "workflow": {"id": workflow_id, "state": workflow_state},
                    },
                }
            }
        )
        return preview_node.execute(make_plot(png)).ui["images"][0]["filename"]

    first_a = execute("workflow-a", png_bytes("red"), {"value": 1})
    first_b = execute("workflow-b", png_bytes("green"))
    second_a_png = png_bytes("blue")
    second_a = execute("workflow-a", second_a_png, {"value": 2})

    assert first_a != first_b
    assert first_a == second_a
    files = list((temp_dir / "openbio-singlecell").iterdir())
    assert len(files) == 2
    assert next(path for path in files if path.name == second_a).read_bytes() == second_a_png


def test_plot_preview_without_workflow_id_has_stable_fallback(comfy_directories):
    _, _, temp_dir = comfy_directories
    filenames = []
    for extra_pnginfo in (None, {"workflow": {}}, {"workflow": {"id": ""}}):
        preview_node = OpenBioSingleCellPreviewResult.PREPARE_CLASS_CLONE(
            {
                "hidden_inputs": {
                    io.Hidden.unique_id: "fallback-node",
                    io.Hidden.extra_pnginfo: extra_pnginfo,
                }
            }
        )
        output = preview_node.execute(make_plot(png_bytes("purple")))
        filenames.append(output.ui["images"][0]["filename"])

    assert len(set(filenames)) == 1
    assert len(list((temp_dir / "openbio-singlecell").iterdir())) == 1


def test_csv_and_png_increment_and_overwrite(comfy_directories, science):
    table_result = make_table(science.pd.DataFrame({"gene": ["A", "B"]}))
    plot_result = make_plot(png_bytes("orange"))

    csv_first = export_csv(table_result, "markers")
    csv_second = export_csv(table_result, "markers")
    csv_overwrite = export_csv(table_result, "markers", overwrite=True)
    png_first = save_png(plot_result, "umap")
    png_second = save_png(plot_result, "umap")
    png_overwrite = save_png(plot_result, "umap", overwrite=True)

    assert csv_first.filename == "markers_00001.csv"
    assert csv_second.filename == "markers_00002.csv"
    assert csv_overwrite.filename == "markers.csv"
    assert png_first.filename == "umap_00001.png"
    assert png_second.filename == "umap_00002.png"
    assert png_overwrite.filename == "umap.png"


def test_csv_preserves_named_and_unnamed_row_identity(comfy_directories, science):
    _, _, _ = comfy_directories
    named = science.pd.DataFrame(
        {"gene": ["A", "B"], "effect": [1.5, science.np.nan]},
        index=science.pd.Index(["cell-1", "cell-2"], name="cell_id"),
    )
    named_target = export_csv(make_table(named), "named")
    unnamed_target = export_csv(make_table(named.rename_axis(None)), "unnamed")

    assert Path(named_target.path).read_text(encoding="utf-8").splitlines() == [
        "cell_id,gene,effect",
        "cell-1,A,1.5",
        "cell-2,B,NA",
    ]
    assert Path(unnamed_target.path).read_text(encoding="utf-8").splitlines()[0] == (
        "__openbio_row_id__,gene,effect"
    )


def test_csv_preserves_multiindex_and_rejects_ambiguous_headers(comfy_directories, science):
    _, _, _ = comfy_directories
    index = science.pd.MultiIndex.from_tuples(
        [("sample-1", "A"), ("sample-2", "B")],
        names=["sample_id", None],
    )
    frame = science.pd.DataFrame({"score": [1.0, 2.0]}, index=index)
    target = export_csv(make_table(frame), "multi-index")
    assert Path(target.path).read_text(encoding="utf-8").splitlines()[0] == (
        "sample_id,__openbio_row_id_1__,score"
    )

    with pytest.raises(ValueError, match="collide"):
        export_csv(make_table(frame.rename_axis(["score", "gene"])), "collision")
    with pytest.raises(ValueError, match="single-level column axis"):
        export_csv(
            make_table(
                science.pd.DataFrame(
                    [[1, 2]],
                    columns=science.pd.MultiIndex.from_tuples([("a", "x"), ("b", "y")]),
                )
            ),
            "multi-columns",
        )


def test_png_validation_preserves_previous_file_and_preview(comfy_directories):
    _, output_dir, temp_dir = comfy_directories
    valid = png_bytes("green")
    durable = save_png(make_plot(valid), "plot", overwrite=True)
    preview_result(make_plot(valid), "preview-node")
    corrupt = make_plot(PNG_SIGNATURE + b"truncated")

    with pytest.raises(ValueError, match="invalid, truncated, or unsafe"):
        save_png(corrupt, "plot", overwrite=True)
    with pytest.raises(ValueError, match="invalid, truncated, or unsafe"):
        preview_result(corrupt, "preview-node")

    assert Path(durable.path).read_bytes() == valid
    preview_files = list((temp_dir / "openbio-singlecell").iterdir())
    assert len(preview_files) == 1
    assert preview_files[0].read_bytes() == valid
    assert not list((output_dir / "openbio-singlecell").glob(".*"))
    assert not list((temp_dir / "openbio-singlecell").glob(".*"))


def test_output_adapters_reject_wrong_result_connections(science):
    summary = make_summary()
    table = make_table(science.pd.DataFrame({"gene": ["A"]}))

    with pytest.raises(TypeError, match="OPENBIO_SINGLE_CELL_TABLE"):
        export_csv(summary, "summary")
    with pytest.raises(TypeError, match="OPENBIO_SINGLE_CELL_PLOT"):
        save_png(table, "table")
    with pytest.raises(TypeError, match="summary, table, or plot"):
        preview_result(object())


def test_h5ad_increment_and_overwrite(comfy_directories, science):
    adata = science.ad.AnnData(science.np.eye(2))
    metadata = ensure_metadata(adata, "roundtrip", {"kind": "test"})
    metadata["analysis_history"] = {
        "000001": {"operation": "first"},
        "000004": {"operation": "fourth"},
    }
    metadata["warnings"] = ["roundtrip warning"]
    record_history(adata, "normalize_total", {"target_sum": 10000.0}, 2, 2)
    first = save_h5ad(adata, "adata")
    second = save_h5ad(adata, "adata")
    overwrite = save_h5ad(adata, "adata", overwrite=True)
    assert first.filename == "adata_00001.h5ad"
    assert second.filename == "adata_00002.h5ad"
    assert overwrite.filename == "adata.h5ad"
    loaded = science.ad.read_h5ad(overwrite.path)
    history = loaded.uns["openbio_singlecell"]["analysis_history"]
    assert isinstance(history, dict)
    assert history["000001"]["operation"] == "first"
    assert history["000004"]["operation"] == "fourth"
    assert history["000005"]["operation"] == "normalize_total"
    assert ensure_metadata(loaded)["warnings"] == ["roundtrip warning"]
    with h5py.File(overwrite.path, "r") as handle:
        assert handle["X"].compression == "gzip"


def test_atomic_writers_preserve_existing_destinations_on_failure(comfy_directories, science, monkeypatch):
    _, output_dir, _ = comfy_directories
    table = make_table(science.pd.DataFrame({"gene": ["A"]}))
    adata = science.ad.AnnData(science.np.eye(2))
    csv_target = export_csv(table, "table", overwrite=True)
    h5ad_target = save_h5ad(adata, "adata", overwrite=True)
    original_csv = Path(csv_target.path).read_bytes()
    original_h5ad = Path(h5ad_target.path).read_bytes()

    def broken_to_csv(self, path, *args, **kwargs):
        Path(path).write_text("partial", encoding="utf-8")
        raise RuntimeError("injected CSV writer failure")

    monkeypatch.setattr(science.pd.DataFrame, "to_csv", broken_to_csv)
    with pytest.raises(RuntimeError, match="injected CSV"):
        export_csv(table, "table", overwrite=True)
    assert Path(csv_target.path).read_bytes() == original_csv
    monkeypatch.undo()

    def broken_write_h5ad(self, path, *args, **kwargs):
        Path(path).write_bytes(b"partial")
        raise RuntimeError("injected H5AD writer failure")

    monkeypatch.setattr(science.ad.AnnData, "write_h5ad", broken_write_h5ad)
    with pytest.raises(RuntimeError, match="injected H5AD"):
        save_h5ad(adata, "adata", overwrite=True)
    assert Path(h5ad_target.path).read_bytes() == original_h5ad
    assert not list((output_dir / "openbio-singlecell").glob(".*"))


def test_h5ad_rejects_invalid_compression_without_creating_output(comfy_directories, science):
    _, output_dir, _ = comfy_directories
    with pytest.raises(ValueError, match="compression must be one of"):
        save_h5ad(science.ad.AnnData(science.np.eye(1)), "adata", compression="zstd")
    assert not list(output_dir.rglob("*.h5ad"))


@pytest.mark.parametrize(("compression", "expected"), [("gzip", "gzip"), ("lzf", "lzf"), ("none", None)])
def test_h5ad_uses_only_declared_portable_compressions(comfy_directories, science, compression, expected):
    _, _, _ = comfy_directories
    target = save_h5ad(science.ad.AnnData(science.np.eye(2)), f"adata-{compression}", compression=compression)
    with h5py.File(target.path, "r") as handle:
        assert handle["X"].compression == expected


def test_h5ad_readback_failure_preserves_existing_destination(comfy_directories, science, monkeypatch):
    _, output_dir, _ = comfy_directories
    adata = science.ad.AnnData(science.np.eye(2))
    target = save_h5ad(adata, "verified", overwrite=True)
    previous = Path(target.path).read_bytes()

    def broken_reader(*args, **kwargs):
        raise RuntimeError("injected H5AD verification failure")

    monkeypatch.setattr(science.ad, "read_h5ad", broken_reader)
    with pytest.raises(RuntimeError, match="injected H5AD verification"):
        save_h5ad(adata, "verified", overwrite=True)

    assert Path(target.path).read_bytes() == previous
    assert not list((output_dir / "openbio-singlecell").glob(".*"))
