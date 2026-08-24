from __future__ import annotations

import json
import tomllib
from pathlib import Path

from openbio_singlecell import PLUGIN_VERSION

PLUGIN_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_EXAMPLES = [
    "example_workflows/Quality Control and Clean Counts.json",
    "example_workflows/Cell Clustering and Marker Discovery.json",
    "example_workflows/Sample Composition Comparison.json",
    "example_workflows/scVI Batch Integration and Contrast.json",
]
EXPECTED_COVERS = [
    "example_workflows/Quality Control and Clean Counts.jpg",
    "example_workflows/Cell Clustering and Marker Discovery.jpg",
    "example_workflows/Sample Composition Comparison.jpg",
    "example_workflows/scVI Batch Integration and Contrast.jpg",
]


def test_release_manifest_records_the_public_contract_and_artifacts():
    manifest = json.loads((PLUGIN_ROOT / "release_manifest.json").read_text(encoding="utf-8"))
    custom_node = manifest["release_artifacts"]["custom_node"]

    assert manifest["schema_version"] == 1
    assert manifest["product"] == "openbio-singlecell"
    assert manifest["version"] == PLUGIN_VERSION
    assert custom_node["example_workflows"] == EXPECTED_EXAMPLES
    assert custom_node["workflow_covers"] == EXPECTED_COVERS
    assert custom_node["contract"] == {
        "node_id_prefix": "OpenBioSingleCell",
        "category_prefix": "openbio/single-cell/",
        "anndata_wire_type": "OPENBIO_ANNDATA",
        "artifact_wire_types": {
            "table": "OPENBIO_SINGLE_CELL_TABLE",
            "plot": "OPENBIO_SINGLE_CELL_PLOT",
            "summary": "OPENBIO_SINGLE_CELL_SUMMARY",
        },
        "analysis_object_wire_types": {
            "scvi_model": "OPENBIO_SCVI_MODEL",
            "scenic_network": "OPENBIO_SCENIC_NETWORK",
            "cassiopeia_tree": "OPENBIO_CASSIOPEIA_TREE",
        },
    }
    assert all((PLUGIN_ROOT / relative).is_file() for relative in custom_node["example_workflows"])
    assert all((PLUGIN_ROOT / relative).is_file() for relative in custom_node["workflow_covers"])
    assert [Path(relative).stem for relative in custom_node["example_workflows"]] == [
        Path(relative).stem for relative in custom_node["workflow_covers"]
    ]
    assert manifest["release_artifacts"]["openbio_frontend_dist"]["build_command"] == "corepack pnpm build:openbio"


def test_python_package_metadata_includes_release_readme():
    metadata = tomllib.loads((PLUGIN_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]

    assert metadata["name"] == "openbio-singlecell"
    assert metadata["version"] == PLUGIN_VERSION
    assert metadata["readme"] == "README.md"
    assert metadata["license"]["text"] == "GPL-3.0-or-later"


def test_release_guidance_uses_the_generic_openbio_frontend_mode():
    paths = [
        PLUGIN_ROOT / "README.md",
        PLUGIN_ROOT / "release_manifest.json",
        PLUGIN_ROOT / "scripts" / "start.ps1",
        PLUGIN_ROOT / "scripts" / "start.sh",
    ]
    release_text = "\n".join(path.read_text(encoding="utf-8") for path in paths)

    assert "build:openbio" in release_text
    assert "build:openbio-singlecell" not in release_text
    assert "openbio/single-cell/" in release_text
