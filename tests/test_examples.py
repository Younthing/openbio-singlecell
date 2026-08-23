from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbio_singlecell.extension import NODE_CLASSES

PLUGIN_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIRECTORY = PLUGIN_ROOT / "example_workflows"
EXPECTED_EXAMPLES = {
    "openbio_singlecell_basic_qc.json",
    "openbio_singlecell_full_analysis.json",
}
WIRE_TYPES = {"OPENBIO_ANNDATA", "OPENBIO_SINGLE_CELL_RESULT"}


def _registered_schemas():
    schemas = [node.GET_SCHEMA() for node in NODE_CLASSES]
    return {schema.node_id: schema for schema in schemas}


def _load_examples() -> dict[str, dict]:
    return {
        path.name: json.loads(path.read_text(encoding="utf-8")) for path in sorted(EXAMPLE_DIRECTORY.glob("*.json"))
    }


def test_example_workflows_are_the_only_packaged_workflow_source():
    examples = _load_examples()

    assert set(examples) == EXPECTED_EXAMPLES
    assert not list((PLUGIN_ROOT / "web" / "workflows").glob("*.json"))


@pytest.mark.parametrize("filename", sorted(EXPECTED_EXAMPLES))
def test_example_workflow_matches_registered_node_schemas(filename):
    workflow = _load_examples()[filename]
    schemas = _registered_schemas()
    nodes = {node["id"]: node for node in workflow["nodes"]}
    links = {link[0]: link for link in workflow["links"]}

    assert workflow["last_node_id"] == max(nodes)
    assert workflow["last_link_id"] == max(links)
    assert len(nodes) == len(workflow["nodes"])
    assert len(links) == len(workflow["links"])

    for node in nodes.values():
        node_id = node["type"]
        assert node_id.startswith("OpenBioSingleCell")
        assert node_id in schemas
        assert node["properties"]["Node name for S&R"] == node_id

        schema = schemas[node_id]
        assert schema.category.startswith("openbio/single-cell/")
        expected_inputs = [(item.id, item.io_type) for item in schema.inputs if item.io_type in WIRE_TYPES]
        actual_inputs = [(item["name"], item["type"]) for item in node.get("inputs", [])]
        assert actual_inputs == expected_inputs

        expected_outputs = [(item.display_name, item.io_type) for item in schema.outputs]
        actual_outputs = [(item["name"], item["type"]) for item in node.get("outputs", [])]
        assert actual_outputs == expected_outputs

        widget_inputs = [item for item in schema.inputs if item.io_type not in WIRE_TYPES]
        assert len(node.get("widgets_values", [])) == len(widget_inputs)

        if node_id == "OpenBioSingleCellSaveH5AD":
            assert node["widgets_values"][0] == "adata"

    for link_id, origin_id, origin_slot, target_id, target_slot, wire_type in links.values():
        origin = nodes[origin_id]
        target = nodes[target_id]
        origin_output = origin["outputs"][origin_slot]
        target_input = target["inputs"][target_slot]

        assert wire_type in WIRE_TYPES
        assert origin_output["type"] == wire_type
        assert target_input["type"] == wire_type
        assert link_id in origin_output["links"]
        assert target_input["link"] == link_id


def test_examples_use_only_the_explicit_anndata_and_result_contracts():
    serialized = "\n".join(json.dumps(value, sort_keys=True) for value in _load_examples().values())

    assert "OPENBIO_ANNDATA" in serialized
    assert "OPENBIO_SINGLE_CELL_RESULT" in serialized
    assert "OPENBIO_SC_DATASET" not in serialized
    assert "OPENBIO_SC_RESULT" not in serialized
    assert '"name": "dataset"' not in serialized
    assert "OpenBioDatasetSummary" not in serialized


def test_examples_use_effective_demo_feature_selection_defaults():
    examples = _load_examples()

    for workflow in examples.values():
        nodes = {node["type"]: node for node in workflow["nodes"]}
        assert nodes["OpenBioSingleCellFilterCells"]["widgets_values"] == [90, 0, 0, 220]
        assert nodes["OpenBioSingleCellFilterGenes"]["widgets_values"] == [100, 0, 0, 0]

    full_nodes = {node["type"]: node for node in examples["openbio_singlecell_full_analysis.json"]["nodes"]}
    assert full_nodes["OpenBioSingleCellHighlyVariableGenes"]["widgets_values"] == [
        200,
        "seurat",
        "X",
        "log1p_norm",
        "",
        False,
    ]
