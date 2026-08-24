from pathlib import Path

WEB_ROOT = Path(__file__).resolve().parents[1] / "web"


def test_preview_uses_extension_lifecycle_without_prototype_patching():
    source = (WEB_ROOT / "openbio_singlecell.js").read_text(encoding="utf-8")

    assert "nodeCreated(node)" in source
    assert "onNodeOutputsUpdated(nodeOutputs)" in source
    assert "beforeRegisterNodeDef" not in source
    assert ".prototype" not in source


def test_preview_renderer_and_styles_are_plugin_owned_assets():
    entrypoint = (WEB_ROOT / "openbio_singlecell.js").read_text(encoding="utf-8")
    renderer = (WEB_ROOT / "single_cell_result_renderer.mjs").read_text(encoding="utf-8")
    stylesheet = (WEB_ROOT / "openbio_singlecell.css").read_text(encoding="utf-8")

    assert "single_cell_result_renderer.mjs" in entrypoint
    assert "openbio_singlecell.css" in entrypoint
    assert "normalizeSingleCellPayload" in renderer
    assert ".openbio-sc-preview" in stylesheet
