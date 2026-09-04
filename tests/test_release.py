from __future__ import annotations

import json
import tomllib
from pathlib import Path

from openbio_singlecell import PLUGIN_VERSION

PLUGIN_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_EXAMPLES = [
    "example_workflows/Quality Control and Clean Counts.json",
    "example_workflows/Cell Clustering and Marker Discovery.json",
    "example_workflows/Subpopulation Reclustering and Annotation.json",
    "example_workflows/Sample Composition Comparison.json",
    "example_workflows/scVI Batch Integration and Contrast.json",
    "example_workflows/Single-Cell Best Practice.json",
]
EXPECTED_COVERS = [
    "example_workflows/Quality Control and Clean Counts.jpg",
    "example_workflows/Cell Clustering and Marker Discovery.jpg",
    "example_workflows/Subpopulation Reclustering and Annotation.jpg",
    "example_workflows/Sample Composition Comparison.jpg",
    "example_workflows/scVI Batch Integration and Contrast.jpg",
    "example_workflows/Single-Cell Best Practice.jpg",
]
EXPECTED_OPTIONAL_DEPENDENCIES = {
    "celltypist": ["celltypist>=1.7,<2"],
    "decoupler": ["decoupler==2.2.0"],
    "harmony": ["harmonypy==2.0.0"],
    "scvi": ["scvi-tools>=1.5,<1.6"],
    "pertpy": ["pertpy==1.3.0"],
    "composition": ["pertpy[tcoda]==1.3.0"],
    "milo": ["pertpy[milo-edger]==1.3.0"],
    "pseudobulk": [
        "decoupler==2.2.0",
        "pertpy[de]==1.3.0",
        "pydeseq2>=0.5,<0.6",
    ],
    "liana": ["liana==1.9.0", "pandas<3"],
    "velocity": ["scvelo==0.3.4"],
    "cnv": ["infercnvpy==0.6.1"],
    "cnmf": ["cnmf==1.7.1"],
    "lineage": ["cassiopeia-mt==2.1.3; platform_system == 'Linux'"],
    "forceatlas2": ["fa2-modified==0.4"],
}
EXPECTED_ARTIFACT_TICKET_WIRE_TYPES = {
    "anndata": "OPENBIO_ANNDATA",
    "table": "OPENBIO_SINGLE_CELL_TABLE",
    "plot": "OPENBIO_SINGLE_CELL_PLOT",
    "augur_result": "OPENBIO_AUGUR_RESULT",
    "milo_result": "OPENBIO_MILO_RESULT",
    "composition_model_result": "OPENBIO_COMPOSITION_MODEL_RESULT",
    "dgidb_resource": "OPENBIO_DGIDB_RESOURCE",
    "liana_result": "OPENBIO_LIANA_RESULT",
    "tf_activity": "OPENBIO_TF_ACTIVITY",
    "pseudobulk": "OPENBIO_SINGLE_CELL_PSEUDOBULK",
    "cnmf_run": "OPENBIO_CNMF_RUN",
    "scvi_model": "OPENBIO_SCVI_MODEL",
    "scenic_result": "OPENBIO_SCENIC_RESULT",
    "cnv_state": "OPENBIO_CNV_STATE",
    "cassiopeia_characters": "OPENBIO_CASSIOPEIA_CHARACTERS",
    "cassiopeia_tree": "OPENBIO_CASSIOPEIA_TREE",
    "velocity_state": "OPENBIO_VELOCITY_STATE",
}


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
        "worker_wire_type": "OPENBIO_WORKER",
        "artifact_ticket_wire_types": EXPECTED_ARTIFACT_TICKET_WIRE_TYPES,
        "reporting_outputs": {
            "summary": "OPENBIO_SINGLE_CELL_SUMMARY",
            "code": "STRING",
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


def test_python_dependency_declarations_match_release_manifest():
    metadata = tomllib.loads((PLUGIN_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    requirements = [
        line.strip()
        for line in (PLUGIN_ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    manifest = json.loads((PLUGIN_ROOT / "release_manifest.json").read_text(encoding="utf-8"))

    assert metadata["dependencies"] == requirements
    assert manifest["python_dependencies"] == {
        "scanpy": requirements[0],
        "anndata": requirements[1],
        "pillow": requirements[2],
    }
    runtime_extras = {
        name: requirements
        for name, requirements in metadata["optional-dependencies"].items()
        if name != "dev"
    }
    assert "ora" not in metadata["optional-dependencies"]
    assert "ora" not in manifest["optional_python_dependencies"]
    assert runtime_extras == EXPECTED_OPTIONAL_DEPENDENCIES
    assert manifest["optional_python_dependencies"] == EXPECTED_OPTIONAL_DEPENDENCIES
    assert manifest["external_backend_contracts"] == {
        "schist": {
            "execution_environment": "comfyui-python",
            "requirements": [
                "conda-forge::schist=0.10.0",
                "conda-forge::graph-tool",
            ],
            "forbidden_pypi_distribution": "schist",
        },
        "pyscenic": {
            "execution_environment": "external-python-3.10-or-aertslab/pyscenic:0.12.1",
            "requirements": [
                "pyscenic==0.12.1",
                "ctxcore==0.2.0",
                "arboreto==0.1.6",
                "loompy==3.0.8",
                "numpy>=1.21,<1.24",
                "pandas>=1.3.5,<2",
            ],
            "integration": "openbio-singlecell/pyscenic-external-run/v1 bundle import",
        },
    }
    assert manifest["excluded_optional_python_dependencies"] == {
        "omicverse": {
            "audited_version": "2.3.1",
            "reason": (
                "omicverse==2.3.1 requires anndata<0.12.0, which conflicts with the core "
                "anndata>=0.13.2,<0.14 contract"
            ),
        }
    }


def test_optional_backend_guidance_is_explicit_and_does_not_expand_core_requirements():
    readme = (PLUGIN_ROOT / "README.md").read_text(encoding="utf-8")
    notices = (PLUGIN_ROOT / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    requirements = (PLUGIN_ROOT / "requirements.txt").read_text(encoding="utf-8")

    assert "conda install -c conda-forge \"schist=0.10.0\" graph-tool" in readme
    assert "package named `schist` on PyPI is an unrelated" in readme
    assert "pyscenic==0.12.1" in readme
    assert "ctxcore==0.2.0" in readme
    assert "numpy>=1.21,<1.24" in readme
    assert "openbio-singlecell/pyscenic-external-run/v1" in readme
    assert "OmicVerse is not declared as an extra" in readme
    assert "DGIdb loader accepts only a user-provided, versioned local CSV/TSV" in readme
    assert "allow_network_access" in readme
    assert "harmonypy==2.0.0" in notices
    assert "cassiopeia-mt==2.1.3; platform_system == 'Linux'" in notices
    assert "scvelo==0.3.4" in notices
    assert "pertpy==1.3.0" in notices
    assert "liana==1.9.0" in notices
    assert "decoupler==2.2.0" in notices
    assert "omicverse" not in requirements.casefold()
    assert "pyscenic" not in requirements.casefold()
    assert "schist" not in requirements.casefold()


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


def test_start_scripts_only_gate_on_the_frontend_entrypoint():
    for relative in ("scripts/start.ps1", "scripts/start.sh"):
        source = (PLUGIN_ROOT / relative).read_text(encoding="utf-8")

        assert "index.html" in source
        assert "THIRD_PARTY_NOTICES.md" not in source
        assert "LICENSE" not in source
