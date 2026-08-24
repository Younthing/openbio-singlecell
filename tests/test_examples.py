from __future__ import annotations

import json
from pathlib import Path

import pytest

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
)
EXPECTED_EXAMPLES = {f"{name}.json" for name in TEMPLATE_NAMES}
EXPECTED_COVERS = {f"{name}.jpg" for name in TEMPLATE_NAMES}
WIRE_TYPES = {
    "OPENBIO_ANNDATA",
    "OPENBIO_SINGLE_CELL_TABLE",
    "OPENBIO_SINGLE_CELL_PLOT",
    "OPENBIO_SINGLE_CELL_SUMMARY",
    "OPENBIO_SCVI_MODEL",
    "OPENBIO_SCENIC_NETWORK",
    "OPENBIO_CASSIOPEIA_TREE",
}


def _input_wire_type(item) -> str:
    return item.get_io_type()


def _is_wire_input(item) -> bool:
    return bool(set(_input_wire_type(item).split(",")) & WIRE_TYPES)


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


def _widgets(workflow: dict, node_type: str) -> dict[str, object]:
    node = _node(workflow, node_type)
    return {
        name.partition(".")[2] or name: value
        for name, value in node.get("widgets_values_named", {}).items()
    }


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


def test_example_workflows_are_the_only_packaged_workflow_source():
    examples = _load_examples()

    assert set(examples) == EXPECTED_EXAMPLES
    assert not list((PLUGIN_ROOT / "web" / "workflows").glob("*.json"))


def test_example_workflows_match_the_schema_driven_generator():
    assert generate_example_workflows(["--check"]) == 0


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
        expected_inputs = [(item.id, _input_wire_type(item)) for item in schema.inputs if _is_wire_input(item)]
        actual_inputs = [(item["name"], item["type"]) for item in node.get("inputs", [])]
        assert actual_inputs == expected_inputs

        expected_outputs = [(item.display_name, item.io_type) for item in schema.outputs]
        actual_outputs = [(item["name"], item["type"]) for item in node.get("outputs", [])]
        assert actual_outputs == expected_outputs

        named_values = node.get("widgets_values_named", {})
        widget_inputs = [item for item in schema.inputs if not _is_wire_input(item)]
        expected_widget_names = selected_widget_names(widget_inputs, named_values)
        assert list(named_values) == expected_widget_names
        assert dict(zip(expected_widget_names, node.get("widgets_values", []), strict=True)) == named_values

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
        "min_genes": 90,
        "max_genes": 0,
        "min_counts": 0,
        "max_counts": 0,
        "max_pct_mito": 20.0,
        "mito_column": "pct_counts_mt",
    }
    assert _widgets(workflow, "OpenBioSingleCellFilterGenes") == {
        "min_cells": 3,
        "max_cells": 0,
        "min_counts": 0,
        "max_counts": 0,
    }
    nodes = {node["id"]: node for node in workflow["nodes"]}
    filter_cells = _node(workflow, "OpenBioSingleCellFilterCells")
    filter_genes = _node(workflow, "OpenBioSingleCellFilterGenes")
    cell_filter_targets = {
        nodes[link["target_id"]]["type"]
        for link in workflow["links"]
        if link["origin_id"] == filter_cells["id"]
    }
    gene_filter_targets = {
        nodes[link["target_id"]]["type"]
        for link in workflow["links"]
        if link["origin_id"] == filter_genes["id"]
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
    }
    assert _widgets(workflow, "OpenBioSingleCellHighlyVariableGenes")["source"] == "layer"
    assert _widgets(workflow, "OpenBioSingleCellHighlyVariableGenes")["layer_name"] == "log1p_norm"
    assert _widgets(workflow, "OpenBioSingleCellPCA")["layer_name"] == "log1p_norm"
    assert _widgets(workflow, "OpenBioSingleCellPCA")["source"] == "layer"
    assert _widgets(workflow, "OpenBioSingleCellMarkerGenes")["source"] == "layer"
    assert _widgets(workflow, "OpenBioSingleCellMarkerGenes")["layer_name"] == "log1p_norm"
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


def test_composition_template_names_required_metadata_and_demo_groups():
    workflow = _load_examples()["Sample Composition Comparison.json"]

    assert _widgets(workflow, "OpenBioSingleCellSampleCompositionSummary") == {
        "sample_key": "sample",
        "group_key": "batch",
        "annotation_key": "cell_type",
    }
    assert _widgets(workflow, "OpenBioSingleCellDifferentialCompositionTest") == {
        "sample_key": "sample",
        "group_key": "batch",
        "annotation_key": "cell_type",
        "control_group": "batch_1",
        "comparison_groups": "batch_2",
        "pseudocount": 0.001,
    }


def test_scvi_template_uses_counts_auto_epochs_and_its_concrete_model_consumer():
    workflow = _load_examples()["scVI Batch Integration and Contrast.json"]

    assert "OpenBioSingleCellSnapshotExpression" not in {node["type"] for node in workflow["nodes"]}
    assert _widgets(workflow, "OpenBioSingleCellNormalizeToLayer")["source"] == "X"
    integration = _widgets(workflow, "OpenBioSingleCellSCVIIntegration")
    assert integration["source"] == "X"
    assert "counts_layer" not in integration
    assert integration["batch_key"] == "batch"
    assert integration["max_epochs"] == 0
    assert integration["output_key"] == "X_scVI"
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
