from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from comfy_api.latest import io

from openbio_singlecell.extension import NODE_CLASSES
from openbio_singlecell.nodes_preprocess import OpenBioSingleCellNormalizeToLayer
from scripts.generate_example_workflows import _node as generated_node
from scripts.generate_example_workflows import main as generate_example_workflows
from tests.workflow_helpers import selected_widget_names, workflow_execute_kwargs

PLUGIN_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIRECTORY = PLUGIN_ROOT / "example_workflows"
TEMPLATE_NAMES = (
    "Quality Control and Clean Counts",
    "Cell Clustering and Marker Discovery",
    "Sample Composition Comparison",
    "scVI Batch Integration and Contrast",
    "Single-Cell Best Practice",
)
EXPECTED_EXAMPLES = {f"{name}.json" for name in TEMPLATE_NAMES}
EXPECTED_COVERS = {f"{name}.jpg" for name in TEMPLATE_NAMES}
WIRE_TYPES = {
    "STRING",
    "OPENBIO_ANNDATA",
    "OPENBIO_SINGLE_CELL_TABLE",
    "OPENBIO_SINGLE_CELL_PLOT",
    "OPENBIO_SINGLE_CELL_SUMMARY",
    "OPENBIO_AUGUR_RESULT",
    "OPENBIO_DGIDB_RESOURCE",
    "OPENBIO_LIANA_RESULT",
    "OPENBIO_TF_ACTIVITY",
    "OPENBIO_SINGLE_CELL_PSEUDOBULK",
    "OPENBIO_CNMF_RUN",
    "OPENBIO_SCVI_MODEL",
    "OPENBIO_SCENIC_RESULT",
    "OPENBIO_SCENIC_BINARY",
    "OPENBIO_CNV_STATE",
    "OPENBIO_CASSIOPEIA_CHARACTERS",
    "OPENBIO_CASSIOPEIA_TREE",
    "OPENBIO_VELOCITY_STATE",
}


def _input_wire_type(item) -> str:
    return item.get_io_type()


def _is_widget_input(item) -> bool:
    return isinstance(item, io.WidgetInput) or item.get_io_type() == "COMFY_DYNAMICCOMBO_V3"


def _is_wire_input(item, wired_input_names: set[str]) -> bool:
    return not _is_widget_input(item) or item.id in wired_input_names


def _registered_schemas():
    schemas = [node.GET_SCHEMA() for node in NODE_CLASSES]
    return {schema.node_id: schema for schema in schemas}


def _load_examples() -> dict[str, dict]:
    return {
        path.name: json.loads(path.read_text(encoding="utf-8")) for path in sorted(EXAMPLE_DIRECTORY.glob("*.json"))
    }


def _node(workflow: dict, node_type: str) -> dict:
    matches = [node for node in workflow["nodes"] if node["type"] == node_type]
    assert len(matches) == 1, f"Expected exactly one {node_type}, found {len(matches)}"
    return matches[0]


def _nodes(workflow: dict, node_type: str) -> list[dict]:
    return [node for node in workflow["nodes"] if node["type"] == node_type]


def _node_with_widget(workflow: dict, node_type: str, widget: str, value: object) -> dict:
    matches = [
        node for node in _nodes(workflow, node_type) if node.get("widgets_values_named", {}).get(widget) == value
    ]
    assert len(matches) == 1, f"Expected one {node_type} with {widget}={value!r}, found {len(matches)}"
    return matches[0]


def _widgets(workflow: dict, node_type: str) -> dict[str, object]:
    node = _node(workflow, node_type)
    return {name.partition(".")[2] or name: value for name, value in node.get("widgets_values_named", {}).items()}


def _assert_link(
    workflow: dict,
    origin_type: str,
    origin_output: str,
    target_type: str,
    target_input: str,
    wire_type: str,
) -> None:
    origin = _node(workflow, origin_type)
    target = _node(workflow, target_type)
    origin_slot = next(index for index, item in enumerate(origin["outputs"]) if item["name"] == origin_output)
    target_slot = next(index for index, item in enumerate(target["inputs"]) if item["name"] == target_input)

    assert any(
        (
            link["origin_id"],
            link["origin_slot"],
            link["target_id"],
            link["target_slot"],
            link["type"],
        )
        == (origin["id"], origin_slot, target["id"], target_slot, wire_type)
        for link in workflow["links"]
    )


def _assert_node_link(
    workflow: dict,
    origin: dict,
    origin_output: str,
    target: dict,
    target_input: str,
    wire_type: str,
) -> None:
    origin_slot = next(index for index, item in enumerate(origin["outputs"]) if item["name"] == origin_output)
    target_slot = next(index for index, item in enumerate(target["inputs"]) if item["name"] == target_input)

    assert any(
        (
            link["origin_id"],
            link["origin_slot"],
            link["target_id"],
            link["target_slot"],
            link["type"],
        )
        == (origin["id"], origin_slot, target["id"], target_slot, wire_type)
        for link in workflow["links"]
    )


def test_example_workflows_are_the_only_packaged_workflow_source():
    examples = _load_examples()

    assert set(examples) == EXPECTED_EXAMPLES
    assert not list((PLUGIN_ROOT / "web" / "workflows").glob("*.json"))


def test_example_workflows_match_the_schema_driven_generator():
    assert generate_example_workflows(["--check"]) == 0


def test_workflow_generator_cli_resolves_comfyui_before_importing_schemas():
    comfy_root = next(Path(entry) for entry in sys.path if (Path(entry) / "main.py").is_file())
    environment = os.environ.copy()
    environment["OPENBIO_COMFYUI_ROOT"] = str(comfy_root)

    completed = subprocess.run(
        [sys.executable, str(PLUGIN_ROOT / "scripts" / "generate_example_workflows.py"), "--check"],
        cwd=PLUGIN_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Verified 5 generated workflows." in completed.stdout


def test_packaged_workflows_exclude_removed_node_ids():
    removed_node_ids = {
        "OpenBioSingleCellUseExpressionLayer",
        "OpenBioSingleCellNormalizeGeneNames",
        "OpenBioSingleCellMarkerORAAnnotation",
        "OpenBioSingleCellDifferentialCompositionTest",
        "OpenBioSingleCellDecouplerPseudobulkContrast",
        "OpenBioSingleCellRunPySCENIC",
        "OpenBioSingleCellLianaResults",
        "OpenBioSingleCellCNVStructure",
    }
    for workflow in _load_examples().values():
        assert removed_node_ids.isdisjoint(node["type"] for node in workflow["nodes"])


def test_workflow_generator_rejects_inactive_or_unnamespaced_dynamic_overrides():
    for overrides in [
        {"source": "X", "source.source_layer": "counts"},
        {"source.source_layer": "counts"},
    ]:
        with pytest.raises(ValueError, match="inactive for the selected dynamic option"):
            generated_node(
                "normalize",
                "OpenBioSingleCellNormalizeToLayer",
                (0, 0),
                overrides,
            )

    with pytest.raises(ValueError, match="has no widget inputs"):
        generated_node(
            "normalize",
            "OpenBioSingleCellNormalizeToLayer",
            (0, 0),
            {"source": "layer", "source_layer": "counts"},
        )

    node = generated_node(
        "normalize",
        "OpenBioSingleCellNormalizeToLayer",
        (0, 0),
        {"source": "layer", "source.source_layer": "counts"},
    )
    assert node["widgets"][:2] == [("source", "layer"), ("source.source_layer", "counts")]
    kwargs = workflow_execute_kwargs(
        OpenBioSingleCellNormalizeToLayer,
        {"widgets_values_named": dict(node["widgets"])},
    )
    assert kwargs["source"] == {"source": "layer", "source_layer": "counts"}

    x_node = generated_node(
        "normalize",
        "OpenBioSingleCellNormalizeToLayer",
        (0, 0),
        {"source": "X"},
    )
    x_kwargs = workflow_execute_kwargs(
        OpenBioSingleCellNormalizeToLayer,
        {"widgets_values_named": dict(x_node["widgets"])},
    )
    assert x_kwargs["source"] == {"source": "X"}

    invalid_x_values = dict(x_node["widgets"])
    invalid_x_values["source.source_layer"] = "counts"
    with pytest.raises(ValueError, match="unknown or inactive widget values"):
        workflow_execute_kwargs(
            OpenBioSingleCellNormalizeToLayer,
            {"widgets_values_named": invalid_x_values},
        )

    unknown_values = dict(x_node["widgets"])
    unknown_values["unknown"] = "value"
    with pytest.raises(ValueError, match="unknown or inactive widget values"):
        workflow_execute_kwargs(
            OpenBioSingleCellNormalizeToLayer,
            {"widgets_values_named": unknown_values},
        )

    missing_values = dict(x_node["widgets"])
    del missing_values["target_sum"]
    with pytest.raises(ValueError, match="missing active widget values"):
        workflow_execute_kwargs(
            OpenBioSingleCellNormalizeToLayer,
            {"widgets_values_named": missing_values},
        )


def test_example_workflow_covers_match_the_exact_template_set():
    from PIL import Image

    json_paths = list(EXAMPLE_DIRECTORY.glob("*.json"))
    cover_paths = list(EXAMPLE_DIRECTORY.glob("*.jpg"))

    assert {path.name for path in json_paths} == EXPECTED_EXAMPLES
    assert {path.name for path in cover_paths} == EXPECTED_COVERS
    assert {path.stem for path in json_paths} == {path.stem for path in cover_paths}

    for cover_path in cover_paths:
        assert cover_path.stat().st_size <= 1_000_000
        with Image.open(cover_path) as cover:
            assert cover.format == "JPEG"
            assert cover.mode == "RGB"
            assert cover.size == (768, 768)
            cover.verify()


@pytest.mark.parametrize("filename", sorted(EXPECTED_EXAMPLES))
def test_example_workflow_matches_registered_node_schemas(filename):
    workflow = _load_examples()[filename]
    schemas = _registered_schemas()
    nodes = {node["id"]: node for node in workflow["nodes"]}
    links = {link["id"]: link for link in workflow["links"]}

    assert workflow["version"] == 1
    assert workflow["state"] == {
        "lastGroupId": max(group["id"] for group in workflow["groups"]),
        "lastNodeId": max(nodes),
        "lastLinkId": max(links),
        "lastRerouteId": 0,
    }
    assert "last_node_id" not in workflow
    assert "last_link_id" not in workflow
    assert len(nodes) == len(workflow["nodes"])
    assert len(links) == len(workflow["links"])

    for link in links.values():
        assert set(link) == {"id", "origin_id", "origin_slot", "target_id", "target_slot", "type"}

    for node in nodes.values():
        node_id = node["type"]
        assert node_id.startswith("OpenBioSingleCell")
        assert node_id in schemas
        assert node["properties"]["Node name for S&R"] == node_id

        schema = schemas[node_id]
        assert schema.category.startswith("openbio/single-cell/")
        actual_inputs = [(item["name"], item["type"]) for item in node.get("inputs", [])]
        wired_input_names = {name for name, _ in actual_inputs}
        expected_inputs = [
            (item.id, _input_wire_type(item)) for item in schema.inputs if _is_wire_input(item, wired_input_names)
        ]
        assert actual_inputs == expected_inputs

        expected_outputs = [(item.display_name, item.io_type) for item in schema.outputs]
        actual_outputs = [(item["name"], item["type"]) for item in node.get("outputs", [])]
        assert actual_outputs == expected_outputs

        named_values = node.get("widgets_values_named", {})
        widget_inputs = [item for item in schema.inputs if _is_widget_input(item)]
        expected_widget_names = selected_widget_names(widget_inputs, named_values)
        assert list(named_values) == expected_widget_names
        assert dict(zip(expected_widget_names, node.get("widgets_values", []), strict=True)) == named_values
        for input_ in node.get("inputs", []):
            schema_input = next(item for item in schema.inputs if item.id == input_["name"])
            if _is_widget_input(schema_input):
                assert input_["widget"] == {"name": input_["name"]}
            else:
                assert "widget" not in input_

    for link in links.values():
        link_id = link["id"]
        origin = nodes[link["origin_id"]]
        target = nodes[link["target_id"]]
        origin_output = origin["outputs"][link["origin_slot"]]
        target_input = target["inputs"][link["target_slot"]]

        assert link["type"] in WIRE_TYPES
        assert origin_output["type"] == link["type"]
        assert link["type"] in target_input["type"].split(",")
        assert link_id in origin_output["links"]
        assert target_input["link"] == link_id


def test_examples_use_only_explicit_data_artifact_and_domain_object_contracts():
    serialized = "\n".join(json.dumps(value, sort_keys=True) for value in _load_examples().values())

    assert "OPENBIO_ANNDATA" in serialized
    assert "OPENBIO_SAMPLE_SHEET" not in serialized
    assert '"type": "STRING"' in serialized
    assert "OPENBIO_SINGLE_CELL_STUDY_DESIGN" not in serialized
    assert "OPENBIO_SINGLE_CELL_TABLE" in serialized
    assert "OPENBIO_SINGLE_CELL_PLOT" in serialized
    assert "OPENBIO_SINGLE_CELL_SUMMARY" in serialized
    assert "OPENBIO_SINGLE_CELL_RESULT" not in serialized
    assert "OPENBIO_SC_DATASET" not in serialized
    assert "OPENBIO_SC_RESULT" not in serialized
    assert '"name": "dataset"' not in serialized
    assert "OpenBioDatasetSummary" not in serialized


def test_quality_control_template_uses_reviewed_production_thresholds():
    examples = _load_examples()
    workflow = examples["Quality Control and Clean Counts.json"]

    assert _widgets(workflow, "OpenBioSingleCellFilterCells") == {
        "source": "X",
        "min_genes": 90,
        "max_genes": 0,
        "min_counts": 0,
        "max_counts": 0,
        "max_pct_mito": 20.0,
        "mito_column": "pct_counts_mt",
        "enable_min_counts": False,
        "enable_max_counts": False,
        "enable_max_pct_mito": True,
    }
    assert _widgets(workflow, "OpenBioSingleCellFilterGenes") == {
        "source": "X",
        "min_cells": 3,
        "max_cells": 0,
        "min_counts": 0,
        "max_counts": 0,
        "enable_min_counts": False,
        "enable_max_counts": False,
    }
    nodes = {node["id"]: node for node in workflow["nodes"]}
    filter_cells = _node(workflow, "OpenBioSingleCellFilterCells")
    filter_genes = _node(workflow, "OpenBioSingleCellFilterGenes")
    cell_filter_targets = {
        nodes[link["target_id"]]["type"] for link in workflow["links"] if link["origin_id"] == filter_cells["id"]
    }
    gene_filter_targets = {
        nodes[link["target_id"]]["type"] for link in workflow["links"] if link["origin_id"] == filter_genes["id"]
    }
    assert "OpenBioSingleCellQCPlots" in cell_filter_targets
    assert "OpenBioSingleCellQCPlots" not in gene_filter_targets


def test_clustering_template_preserves_counts_and_uses_a_layer_aware_pipeline():
    workflow = _load_examples()["Cell Clustering and Marker Discovery.json"]
    node_types = {node["type"] for node in workflow["nodes"]}

    assert "OpenBioSingleCellSnapshotExpression" not in node_types
    assert "OpenBioSingleCellNormalizeToLayer" in node_types
    assert "OpenBioSingleCellScale" not in node_types
    assert _widgets(workflow, "OpenBioSingleCellNormalizeToLayer") == {
        "source": "X",
        "target_sum": 10000.0,
        "transform": "log1p",
        "output_layer": "log1p_norm",
        "overwrite_existing": False,
    }
    assert _widgets(workflow, "OpenBioSingleCellHighlyVariableGenes")["source"] == "layer"
    assert _widgets(workflow, "OpenBioSingleCellHighlyVariableGenes")["layer_name"] == "log1p_norm"
    assert _widgets(workflow, "OpenBioSingleCellPCA")["layer_name"] == "log1p_norm"
    assert _widgets(workflow, "OpenBioSingleCellPCA")["source"] == "layer"
    assert _widgets(workflow, "OpenBioSingleCellMarkerGenes")["source"] == "layer"
    assert _widgets(workflow, "OpenBioSingleCellMarkerGenes")["layer_name"] == "log1p_norm"
    assert _widgets(workflow, "OpenBioSingleCellUMAPPlot") == {
        "embedding_key": "X_umap",
        "color": "leiden",
        "color_mode": "auto",
        "point_size": 10.0,
        "continuous_color_map": "viridis",
        "categorical_palette": "tab20",
        "sort_order": True,
        "missing_color": "lightgray",
        "legend_policy": "automatic",
    }
    _assert_link(
        workflow,
        "OpenBioSingleCellLoadH5AD",
        "adata",
        "OpenBioSingleCellNormalizeToLayer",
        "adata",
        "OPENBIO_ANNDATA",
    )
    _assert_link(
        workflow,
        "OpenBioSingleCellNormalizeToLayer",
        "adata",
        "OpenBioSingleCellHighlyVariableGenes",
        "adata",
        "OPENBIO_ANNDATA",
    )


def test_composition_template_fans_out_the_demo_study_parameters_as_strings():
    workflow = _load_examples()["Sample Composition Comparison.json"]

    node_types = {node["type"] for node in workflow["nodes"]}
    assert "OpenBioSingleCellSampleSheet" not in node_types
    assert "OpenBioSingleCellApplySampleMetadata" not in node_types
    study_parameters = json.loads(_widgets(workflow, "OpenBioSingleCellCoreStudyParameters")["study_parameters_json"])
    assert study_parameters == {
        "schema_version": 1,
        "sample_column": "sample",
        "condition_column": "condition",
        "batch_column": "batch",
        "annotation_column": "cell_type",
        "reference": "control",
        "comparison": "treated",
    }
    assert _widgets(workflow, "OpenBioSingleCellSampleCompositionSummary") == {
        "sample_key": "sample",
        "condition_key": "condition",
        "annotation_key": "cell_type",
        "annotation_status": "unknown",
        "max_output_rows": 2_000_000,
    }
    _assert_link(
        workflow,
        "OpenBioSingleCellLoadH5AD",
        "adata",
        "OpenBioSingleCellSampleCompositionSummary",
        "adata",
        "OPENBIO_ANNDATA",
    )
    for source_output, target_input in [
        ("sample_column", "sample_key"),
        ("condition_column", "condition_key"),
        ("annotation_column", "annotation_key"),
    ]:
        _assert_link(
            workflow,
            "OpenBioSingleCellCoreStudyParameters",
            source_output,
            "OpenBioSingleCellSampleCompositionSummary",
            target_input,
            "STRING",
        )

    _assert_link(
        workflow,
        "OpenBioSingleCellLoadH5AD",
        "adata",
        "OpenBioSingleCellAnnDataSummary",
        "adata",
        "OPENBIO_ANNDATA",
    )

    core_parameters = _node(workflow, "OpenBioSingleCellCoreStudyParameters")
    assert core_parameters["inputs"] == []
    for output_name in ("batch_column", "reference", "comparison"):
        output = next(item for item in core_parameters["outputs"] if item["name"] == output_name)
        assert output["links"] == []

    composition = _node(workflow, "OpenBioSingleCellSampleCompositionSummary")
    connected_names = {"sample_key", "condition_key", "annotation_key"}
    linked_widgets = {item["name"]: item["widget"] for item in composition["inputs"] if item["name"] in connected_names}
    assert linked_widgets == {name: {"name": name} for name in connected_names}
    summary_link = next(
        link
        for link in workflow["links"]
        if link["origin_id"] == composition["id"] and composition["outputs"][link["origin_slot"]]["name"] == "summary"
    )
    summary_preview = next(node for node in workflow["nodes"] if node["id"] == summary_link["target_id"])
    _assert_node_link(
        workflow,
        composition,
        "summary",
        summary_preview,
        "result",
        "OPENBIO_SINGLE_CELL_SUMMARY",
    )
    assert summary_preview["type"] == "OpenBioSingleCellPreviewResult"


def test_scvi_template_uses_counts_auto_epochs_and_its_concrete_model_consumer():
    workflow = _load_examples()["scVI Batch Integration and Contrast.json"]

    assert "OpenBioSingleCellSnapshotExpression" not in {node["type"] for node in workflow["nodes"]}
    assert _widgets(workflow, "OpenBioSingleCellNormalizeToLayer")["source"] == "X"
    integration = _widgets(workflow, "OpenBioSingleCellSCVIIntegration")
    assert integration["source"] == "X"
    assert "counts_layer" not in integration
    assert integration["technical_batch_key"] == "batch"
    assert integration["epochs"] == "automatic"
    assert integration["output_key"] == "X_scVI"
    assert _widgets(workflow, "OpenBioSingleCellSCVIDifferentialExpression") == {
        "groupby": "cell_type",
        "group1": "T cell",
        "group2": "B cell",
        "population_scope": "all",
        "mode": "change",
        "delta": 0.25,
        "fdr_target": 0.05,
        "batch_handling": "shared_technical_batches",
        "n_samples_overall": 5000,
        "random_seed": 0,
    }
    assert _widgets(workflow, "OpenBioSingleCellNeighbors")["use_rep"] == "X_scVI"
    _assert_link(
        workflow,
        "OpenBioSingleCellFilterGenes",
        "adata",
        "OpenBioSingleCellNormalizeToLayer",
        "adata",
        "OPENBIO_ANNDATA",
    )
    _assert_link(
        workflow,
        "OpenBioSingleCellNormalizeToLayer",
        "adata",
        "OpenBioSingleCellSCVIIntegration",
        "adata",
        "OPENBIO_ANNDATA",
    )
    _assert_link(
        workflow,
        "OpenBioSingleCellSCVIIntegration",
        "model",
        "OpenBioSingleCellSCVIDifferentialExpression",
        "model",
        "OPENBIO_SCVI_MODEL",
    )
    _assert_link(
        workflow,
        "OpenBioSingleCellSCVIIntegration",
        "adata",
        "OpenBioSingleCellSCVIDifferentialExpression",
        "adata",
        "OPENBIO_ANNDATA",
    )


def test_best_practice_template_uses_existing_nodes_for_the_reviewed_two_branch_contract():
    workflow = _load_examples()["Single-Cell Best Practice.json"]
    node_types = {node["type"] for node in workflow["nodes"]}

    assert _widgets(workflow, "OpenBioSingleCellLoadH5AD")["path"] == ("openbio-singlecell/anndata_qc.h5ad")
    study = json.loads(_widgets(workflow, "OpenBioSingleCellCoreStudyParameters")["study_parameters_json"])
    assert study == {
        "schema_version": 1,
        "sample_column": "sample",
        "condition_column": "group",
        "batch_column": "batch",
        "annotation_column": "celltypist_cell_type",
        "reference": "Normal",
        "comparison": "nonDM_ED",
    }

    hard_filter = _node_with_widget(workflow, "OpenBioSingleCellFilterCells", "min_genes", 200)
    mitochondrial_cap = _node_with_widget(workflow, "OpenBioSingleCellFilterCells", "max_pct_mito", 20.0)
    assert hard_filter["widgets_values_named"]["max_pct_mito"] == 0.0
    assert hard_filter["widgets_values_named"]["enable_max_pct_mito"] is False
    assert mitochondrial_cap["widgets_values_named"]["min_genes"] == 0
    assert mitochondrial_cap["widgets_values_named"]["enable_max_pct_mito"] is True
    assert _widgets(workflow, "OpenBioSingleCellFilterGenes")["min_cells"] == 3

    general_mad = _node_with_widget(
        workflow,
        "OpenBioSingleCellMarkMADOutliers",
        "output_column",
        "qc_mad_general",
    )
    mitochondrial_mad = _node_with_widget(
        workflow,
        "OpenBioSingleCellMarkMADOutliers",
        "output_column",
        "qc_mad_mito",
    )
    assert general_mad["widgets_values_named"] == {
        "metrics": "total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
        "batch_key": "sample",
        "nmads": 5.0,
        "direction": "both",
        "output_column": "qc_mad_general",
        "scale_mad": False,
        "minimum_group_size": 3,
    }
    assert mitochondrial_mad["widgets_values_named"] == {
        "metrics": "pct_counts_mt",
        "batch_key": "sample",
        "nmads": 3.0,
        "direction": "upper",
        "output_column": "qc_mad_mito",
        "scale_mad": False,
        "minimum_group_size": 3,
    }
    assert len(_nodes(workflow, "OpenBioSingleCellScrublet")) == 1
    assert len(_nodes(workflow, "OpenBioSingleCellFilterDoublets")) == 1

    snapshot = _node(workflow, "OpenBioSingleCellSnapshotExpression")
    assert snapshot["widgets_values_named"] == {
        "source": "X",
        "overwrite_existing": False,
    }
    assert _widgets(workflow, "OpenBioSingleCellNormalizeToLayer") == {
        "source": "layer",
        "source_layer": "counts",
        "target_sum": 10000.0,
        "transform": "log1p",
        "output_layer": "log1p_norm",
        "overwrite_existing": False,
    }

    hvg_nodes = _nodes(workflow, "OpenBioSingleCellHighlyVariableGenes")
    assert len(hvg_nodes) == 2
    assert {node["widgets_values_named"]["subset"] for node in hvg_nodes} == {False, True}
    for node in hvg_nodes:
        assert node["widgets_values_named"]["n_top_genes"] == 5000
        assert node["widgets_values_named"]["flavor"] == "seurat"
        assert node["widgets_values_named"]["source"] == "layer"
        assert node["widgets_values_named"]["source.layer_name"] == "log1p_norm"

    scvi = _node(workflow, "OpenBioSingleCellSCVIIntegration")
    assert scvi["widgets_values_named"]["source"] == "layer"
    assert scvi["widgets_values_named"]["source.counts_layer"] == "counts"
    assert scvi["widgets_values_named"]["output_key"] == "X_scVI"
    assert _widgets(workflow, "OpenBioSingleCellNeighbors")["use_rep"] == "X_scVI"
    assert "OpenBioSingleCellSCVIDifferentialExpression" not in node_types

    hvg_only = next(node for node in hvg_nodes if node["widgets_values_named"]["subset"] is True)
    hvg_full = next(node for node in hvg_nodes if node["widgets_values_named"]["subset"] is False)
    merge = _node(workflow, "OpenBioSingleCellMergeObservationAnnotations")
    assert merge["widgets_values_named"]["conflict_policy"] == "error"
    leiden = _node(workflow, "OpenBioSingleCellLeiden")
    _assert_node_link(workflow, hvg_only, "adata", scvi, "adata", "OPENBIO_ANNDATA")
    _assert_node_link(workflow, hvg_full, "adata", merge, "adata", "OPENBIO_ANNDATA")
    _assert_node_link(workflow, leiden, "adata", merge, "subset_adata", "OPENBIO_ANNDATA")

    marker_nodes = _nodes(workflow, "OpenBioSingleCellMarkerGenes")
    assert len(marker_nodes) == 2
    assert all(node["widgets_values_named"]["groupby"] == "leiden_scvi" for node in marker_nodes)
    assert all(node["widgets_values_named"]["source.layer_name"] == "log1p_norm" for node in marker_nodes)
    marker_plot = _node(workflow, "OpenBioSingleCellMarkerExpressionPlot")
    assert marker_plot["widgets_values_named"] == {
        "genes": "PECAM1,VWF,KDR,COL1A1,DCN,RGS5,CSPG4,ACTA2,TAGLN,SOX10,S100B,CD68,LYZ,CD3D,CD3E,NKG7,MS4A1,CD79A",
        "groupby": "leiden_scvi",
        "plot": "dotplot",
        "plot.standard_scale": "var",
        "plot.expression_cutoff": 0.0,
        "plot.mean_only_expressed": False,
        "source": "layer",
        "source.layer_name": "log1p_norm",
        "group_order": "observed",
        "random_seed": 0,
    }

    celltypist = _node(workflow, "OpenBioSingleCellCellTypistAnnotation")
    assert celltypist["widgets_values_named"] == {
        "source": "raw",
        "expression_state": "verified_counts",
        "model": "Adult_Human_Vascular.pkl",
        "majority_voting": False,
        "over_clustering_key": "",
        "min_prop": 0.0,
        "label_column": "celltypist_cell_type",
        "confidence_column": "celltypist_confidence",
        "probability_key": "celltypist_probabilities",
        "store_decision_matrix": False,
        "decision_key": "celltypist_decision_scores",
        "metadata_key": "celltypist",
        "overwrite_existing": False,
    }
    pseudobulk = _node(workflow, "OpenBioSingleCellPseudobulk")
    deseq2 = _node(workflow, "OpenBioSingleCellPseudobulkDESeq2")
    assert pseudobulk["widgets_values_named"] == {
        "sample_key": "sample",
        "population_key": "celltypist_cell_type",
        "condition_key": "group",
        "technical_batch_key": "batch",
        "categorical_covariate_keys": "",
        "continuous_covariate_keys": "",
        "source": "layer",
        "source.layer_name": "counts",
        "inference_mode": "exploratory",
        "min_cells": 10,
        "min_counts": 1000,
    }
    assert all(
        node["widgets_values_named"]["missing_policy"] == "exclude"
        for node in _nodes(workflow, "OpenBioSingleCellSubsetObservations")
    )
    assert deseq2["widgets_values_named"] == {
        "population": "smc_pc_intermediate",
        "reference_condition": "Normal",
        "comparison_condition": "nonDM_ED",
        "categorical_covariate_keys": "",
        "continuous_covariate_keys": "",
        "fdr_threshold": 0.05,
        "min_abs_log2_fold_change": 0,
        "min_count": 10,
        "min_total_count": 15,
        "large_n": 10,
        "min_prop": 0.7,
        "n_cpus": 1,
    }
    _assert_node_link(
        workflow,
        pseudobulk,
        "pseudobulk",
        deseq2,
        "pseudobulk",
        "OPENBIO_SINGLE_CELL_PSEUDOBULK",
    )

    save_nodes = _nodes(workflow, "OpenBioSingleCellSaveH5AD")
    assert {node["widgets_values_named"]["filename_prefix"] for node in save_nodes} == {
        "best_practice_full_gene_results",
        "best_practice_hvg_scvi_results",
    }


def test_best_practice_template_only_links_the_two_condition_values_from_study_parameters():
    workflow = _load_examples()["Single-Cell Best Practice.json"]
    nodes_by_id = {node["id"]: node for node in workflow["nodes"]}
    study_parameters = _node(workflow, "OpenBioSingleCellCoreStudyParameters")

    outgoing_links = {
        (
            study_parameters["outputs"][link["origin_slot"]]["name"],
            nodes_by_id[link["target_id"]]["type"],
            nodes_by_id[link["target_id"]]["inputs"][link["target_slot"]]["name"],
        )
        for link in workflow["links"]
        if link["origin_id"] == study_parameters["id"]
    }

    assert outgoing_links == {
        ("reference", "OpenBioSingleCellPseudobulkDESeq2", "reference_condition"),
        ("comparison", "OpenBioSingleCellPseudobulkDESeq2", "comparison_condition"),
    }

    assert _widgets(workflow, "OpenBioSingleCellScrublet")["batch_key"] == "sample"
    assert {
        node["widgets_values_named"]["batch_key"] for node in _nodes(workflow, "OpenBioSingleCellMarkMADOutliers")
    } == {"sample"}
    assert {
        node["widgets_values_named"]["batch_key"] for node in _nodes(workflow, "OpenBioSingleCellHighlyVariableGenes")
    } == {"sample"}
    assert _widgets(workflow, "OpenBioSingleCellSCVIIntegration")["technical_batch_key"] == "batch"
    assert _widgets(workflow, "OpenBioSingleCellSampleCompositionSummary") == {
        "sample_key": "sample",
        "condition_key": "group",
        "annotation_key": "celltypist_cell_type",
        "annotation_status": "provisional",
        "max_output_rows": 2_000_000,
    }
    pseudobulk = _widgets(workflow, "OpenBioSingleCellPseudobulk")
    assert pseudobulk["sample_key"] == "sample"
    assert pseudobulk["population_key"] == "celltypist_cell_type"
    assert pseudobulk["condition_key"] == "group"
    assert pseudobulk["technical_batch_key"] == "batch"
    assert (
        _node_with_widget(
            workflow,
            "OpenBioSingleCellSubsetObservations",
            "values",
            "Normal,nonDM_ED",
        )["widgets_values_named"]["column"]
        == "group"
    )
    deseq2 = _widgets(workflow, "OpenBioSingleCellPseudobulkDESeq2")
    assert deseq2["population"] == "smc_pc_intermediate"
    assert deseq2["reference_condition"] == "Normal"
    assert deseq2["comparison_condition"] == "nonDM_ED"
