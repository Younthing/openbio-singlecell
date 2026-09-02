from pathlib import Path

WEB_ROOT = Path(__file__).resolve().parents[1] / "web"


def test_preview_uses_extension_lifecycle_without_prototype_patching():
    source = (WEB_ROOT / "openbio_singlecell.js").read_text(encoding="utf-8")

    assert "getCustomWidgets()" in source
    assert "coreStudyParameterWidgets" in source
    assert "nodeCreated(node)" in source
    assert "onNodeOutputsUpdated(nodeOutputs)" in source
    assert "workflow_migrations" not in source
    assert "migrateExpressionStateWorkflow" not in source
    assert "beforeConfigureGraph" not in source
    assert "beforeRegisterNodeDef" not in source
    assert ".prototype" not in source
    assert not (WEB_ROOT / "workflow_migrations.mjs").exists()


def test_preview_renderer_and_styles_are_plugin_owned_assets():
    entrypoint = (WEB_ROOT / "openbio_singlecell.js").read_text(encoding="utf-8")
    renderer = (WEB_ROOT / "single_cell_result_renderer.mjs").read_text(encoding="utf-8")
    stylesheet = (WEB_ROOT / "openbio_singlecell.css").read_text(encoding="utf-8")

    assert "single_cell_result_renderer.mjs" in entrypoint
    assert "openbio_singlecell.css" in entrypoint
    assert "normalizeSingleCellPayload" in renderer
    assert ".openbio-sc-preview" in stylesheet
    assert ".openbio-sc-preview__details" in stylesheet
    assert ".openbio-sc-preview__index-cell" in stylesheet
    assert ".openbio-sc-output" in stylesheet
    assert ".openbio-sc-output__action:focus-visible" in stylesheet


def test_core_study_widget_is_a_small_plugin_owned_module():
    entrypoint = (WEB_ROOT / "openbio_singlecell.js").read_text(encoding="utf-8")
    study_widgets = (WEB_ROOT / "study_input_widgets.mjs").read_text(encoding="utf-8")

    assert "study_input_widgets.mjs" in entrypoint
    assert "OPENBIO_CORE_STUDY_PARAMETERS_WIDGET" in study_widgets
    assert 'CORE_STUDY_PARAMETERS_WIDGET = "study_parameters_json"' in study_widgets
    assert "OPENBIO_SAMPLE_SHEET_WIDGET" not in study_widgets
    assert "sample_sheet" not in study_widgets.lower()
    assert not (WEB_ROOT / "editable_table.mjs").exists()
    assert "OPENBIO_STUDY_DESIGN_WIDGET" not in study_widgets
    assert "study_design_json" not in study_widgets
    assert ".prototype" not in study_widgets


def test_input_file_upload_widget_is_plugin_owned_and_uses_comfy_uploads():
    entrypoint = (WEB_ROOT / "openbio_singlecell.js").read_text(encoding="utf-8")
    upload_widget = (WEB_ROOT / "input_file_upload_widget.mjs").read_text(encoding="utf-8")
    stylesheet = (WEB_ROOT / "openbio_singlecell.css").read_text(encoding="utf-8")

    assert "input_file_upload_widget.mjs" in entrypoint
    assert "/scripts/api.js" in entrypoint
    assert "OPENBIO_INPUT_FILE_UPLOAD_WIDGET" in upload_widget
    assert 'fetchApi("/upload/image"' in upload_widget
    assert 'body.append("type", "input")' in upload_widget
    assert "node.onDragOver" in upload_widget
    assert "node.onDragDrop" in upload_widget
    assert ".openbio-input-file__drop-zone" in stylesheet
