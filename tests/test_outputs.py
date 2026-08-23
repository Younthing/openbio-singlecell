from __future__ import annotations

import pytest
from comfy_api.latest import io

from openbio_singlecell import dependencies
from openbio_singlecell.contracts import SingleCellResult, ensure_metadata, record_history
from openbio_singlecell.nodes_output import (
    OpenBioSingleCellPreviewResult,
    OpenBioSingleCellSaveH5AD,
    export_csv,
    preview_result,
    save_h5ad,
    save_png,
)


def make_result(kind, *, table=None, png=None):
    return SingleCellResult(
        kind=kind,
        title="result",
        parameters={},
        description="description",
        warnings=[],
        input_cells=2,
        input_genes=2,
        random_seed=0,
        elapsed_seconds=0.1,
        source={},
        table=table,
        png=png,
    )


def test_save_h5ad_schema_uses_anndata_contract():
    schema = OpenBioSingleCellSaveH5AD.GET_SCHEMA()

    assert schema.inputs[0].id == "adata"
    assert schema.inputs[1].id == "filename_prefix"
    assert schema.inputs[1].default == "adata"


def test_preview_does_not_persist_tables(comfy_directories):
    _, output_dir, temp_dir = comfy_directories
    ui = preview_result(make_result("summary"), "node")
    assert "openbio_singlecell" in ui
    assert not list(output_dir.rglob("*"))
    assert not list(temp_dir.rglob("*"))


def test_plot_preview_overwrites_per_node(comfy_directories):
    _, _, temp_dir = comfy_directories
    first = preview_result(make_result("plot", png=b"first"), "same-node")
    second = preview_result(make_result("plot", png=b"second"), "same-node")
    files = list((temp_dir / "openbio-singlecell").iterdir())
    assert len(files) == 1
    assert files[0].read_bytes() == b"second"
    assert first["images"][0]["filename"] == second["images"][0]["filename"]
    assert first["images"][0]["type"] == "temp"


def test_plot_preview_uses_v3_hidden_unique_id(comfy_directories):
    _, _, temp_dir = comfy_directories
    preview_node = OpenBioSingleCellPreviewResult.PREPARE_CLASS_CLONE(
        {"hidden_inputs": {io.Hidden.unique_id: "workflow-node-42"}}
    )
    output = preview_node.execute(make_result("plot", png=b"preview"))
    expected = preview_result(make_result("plot", png=b"preview"), "workflow-node-42")

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
        return preview_node.execute(make_result("plot", png=png)).ui["images"][0]["filename"]

    first_a = execute("workflow-a", b"first-a", {"value": 1})
    first_b = execute("workflow-b", b"first-b")
    second_a = execute("workflow-a", b"second-a", {"value": 2})

    assert first_a != first_b
    assert first_a == second_a
    files = list((temp_dir / "openbio-singlecell").iterdir())
    assert len(files) == 2
    assert next(path for path in files if path.name == second_a).read_bytes() == b"second-a"


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
        output = preview_node.execute(make_result("plot", png=b"fallback"))
        filenames.append(output.ui["images"][0]["filename"])

    assert len(set(filenames)) == 1
    assert len(list((temp_dir / "openbio-singlecell").iterdir())) == 1


def test_csv_and_png_increment_and_overwrite(comfy_directories):
    if not dependencies.AVAILABLE:
        pytest.skip("Scientific dependencies are unavailable.")
    table_result = make_result("table", table=dependencies.pd.DataFrame({"gene": ["A", "B"]}))
    plot_result = make_result("plot", png=b"png")

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


def test_h5ad_increment_and_overwrite(comfy_directories):
    if not dependencies.AVAILABLE:
        pytest.skip("Scientific dependencies are unavailable.")
    adata = dependencies.ad.AnnData(dependencies.np.eye(2))
    metadata = ensure_metadata(adata, "roundtrip", {"kind": "test"})
    metadata["analysis_history"] = {
        "000001": {"operation": "legacy_one"},
        "000004": {"operation": "legacy_four"},
    }
    metadata["warnings"] = ["roundtrip warning"]
    record_history(adata, "normalize_total", {"target_sum": 10000.0}, 2, 2)
    first = save_h5ad(adata, "adata")
    second = save_h5ad(adata, "adata")
    overwrite = save_h5ad(adata, "adata", overwrite=True)
    assert first.filename == "adata_00001.h5ad"
    assert second.filename == "adata_00002.h5ad"
    assert overwrite.filename == "adata.h5ad"
    loaded = dependencies.ad.read_h5ad(overwrite.path)
    history = loaded.uns["openbio_singlecell"]["analysis_history"]
    assert isinstance(history, dict)
    assert history["000001"]["operation"] == "legacy_one"
    assert history["000004"]["operation"] == "legacy_four"
    assert history["000005"]["operation"] == "normalize_total"
    assert ensure_metadata(loaded)["warnings"] == ["roundtrip warning"]


@pytest.mark.parametrize("legacy_type", ["list", "tuple", "ndarray"])
def test_legacy_history_sequences_are_normalized(legacy_type):
    if not dependencies.AVAILABLE:
        pytest.skip("Scientific dependencies are unavailable.")
    entry = {"operation": "legacy"}
    values = [entry]
    if legacy_type == "tuple":
        values = tuple(values)
    elif legacy_type == "ndarray":
        values = dependencies.np.asarray(values, dtype=object)

    adata = dependencies.ad.AnnData(dependencies.np.eye(1))
    adata.uns["openbio_singlecell"] = {"analysis_history": values}
    history = ensure_metadata(adata)["analysis_history"]

    assert isinstance(history, dict)
    assert history == {"000000": entry}
