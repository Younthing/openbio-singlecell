from __future__ import annotations

import asyncio
import gc
import hashlib
import json
import sys
import uuid
from pathlib import Path

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from openbio_singlecell.artifact_runtime import ArtifactTicket
from openbio_singlecell.artifact_service import (
    current_artifact_runtime,
    execute_artifact_node,
    initialize_artifact_service,
)
from openbio_singlecell.dgidb_artifact_codec import read_dgidb_resource
from openbio_singlecell.dgidb_resource import validate_dgidb_resource
from openbio_singlecell.nodes_enrichment import OpenBioSingleCellDGIdbAnnotation
from openbio_singlecell.operations_enrichment import (
    gene_panel_scores_operation,
    load_dgidb_resource_operation,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse, registered_operation_ids


def _metadata_json() -> str:
    return json.dumps(
        {
            "name": "DGIdb deterministic test snapshot",
            "version": "5.0-test",
            "release_date": "2026-08-28",
            "download_url": "https://dgidb.org/downloads",
            "organism": "Homo sapiens",
            "gene_identifier_namespace": "HGNC symbol",
            "scope": "research-only deterministic test resource",
            "license": "DGIdb aggregate; reviewed source-specific terms",
            "source_license_review": "reviewed: deterministic fixture sources are redistributable",
            "citation": "Deterministic DGIdb test snapshot.",
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )


def _write_resource(path: Path) -> None:
    path.write_text(
        "drug_claim_name\tgene_claim_name\tinteraction_claim_source\tinteraction_types\n"
        "DrugA\tG1\tSource1\tinhibitor\n"
        "DrugA\tG2\tSource2\tbinder\n"
        "DrugB\tG3\tSource1\tactivator\n",
        encoding="utf-8-sig",
    )


def _gene_set_metadata_json() -> str:
    return json.dumps(
        {
            "name": "reviewed-test-sets",
            "version": "2026.08",
            "date": "2026-08-01",
            "organism": "human",
            "identifier_namespace": "HGNC symbol",
            "scope": "unit-test scientific contract",
            "license": "CC0-1.0",
            "citation": "OpenBio test fixture, 2026.",
        },
        separators=(",", ":"),
    )


def test_load_dgidb_operation_publishes_jsonl_artifact_without_rewriting_input(tmp_path: Path):
    resource_path = tmp_path / "dgidb.tsv"
    _write_resource(resource_path)
    before = hashlib.sha256(resource_path.read_bytes()).hexdigest()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = load_dgidb_resource_operation(
        context,
        {
            "resource_file": {
                "type": "file",
                "path": str(resource_path.resolve()),
                "provenance": {"path": resource_path.name},
            }
        },
        {
            "resource_metadata_json": _metadata_json(),
            "drug_column": "drug_claim_name",
            "gene_column": "gene_claim_name",
            "source_column": "interaction_claim_source",
            "evidence_column": "interaction_types",
            "max_file_bytes": 1_000_000,
            "max_rows": 1_000,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert hashlib.sha256(resource_path.read_bytes()).hexdigest() == before
    assert records[0]["kind"] == "OPENBIO_DGIDB_RESOURCE"
    assert records[0]["codec"] == "dgidb-table-jsonl-v1"
    resource = read_dgidb_resource(staging / records[0]["payload"])
    owned_table = object.__getattribute__(resource, "_table")
    assert validate_dgidb_resource(resource, copy_payload=False)[0] is owned_table
    assert resource.table().loc[:, ["drug", "gene"]].values.tolist() == [
        ["DrugA", "G1"],
        ["DrugA", "G2"],
        ["DrugB", "G3"],
    ]
    assert records[1]["value"]["summary"]["node_id"] == "OpenBioSingleCellDGIdbAnnotation"


def test_gene_panel_operation_publishes_new_anndata_without_rewriting_input(tmp_path: Path, science):
    resource_path = tmp_path / "sets.csv"
    resource_path.write_text(
        "geneset,genesymbol\npanel_x,G0\npanel_x,G1\npanel_x,G2\n",
        encoding="utf-8",
    )
    counts = science.np.arange(1, 41, dtype=float).reshape(8, 5)
    adata = science.ad.AnnData(counts)
    adata.obs_names = [f"cell_{index}" for index in range(8)]
    adata.var_names = [f"G{index}" for index in range(5)]
    science.sc.pp.normalize_total(adata, target_sum=10_000.0)
    science.sc.pp.log1p(adata)
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).hexdigest()
    staging = tmp_path / "score.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = gene_panel_scores_operation(
        context,
        {
            "adata": {
                "type": "artifact",
                "path": str(input_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            },
            "gene_sets_file": {
                "type": "file",
                "path": str(resource_path.resolve()),
                "provenance": {"path": resource_path.name},
            },
        },
        {
            "resource_metadata_json": _gene_set_metadata_json(),
            "panel": "panel_x",
            "source": {"source": "X"},
            "source_column": "geneset",
            "target_column": "genesymbol",
            "output_key": "panel_score",
            "ctrl_size": 1,
            "n_bins": 2,
            "random_seed": 0,
            "overwrite_existing": False,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == before
    output = read_anndata(staging / records[0]["payload"])
    assert "panel_score" in output.obs
    assert records[0]["codec"] == "anndata-h5ad-v1"


def test_all_enrichment_nodes_have_allowlisted_operations_and_schema_only_adapters():
    expected = {
        "openbio.node.aucellscores",
        "openbio.node.dgidbannotation",
        "openbio.node.druggsea",
        "openbio.node.drughypergeometric",
        "openbio.node.drugscores",
        "openbio.node.genepanelscores",
        "openbio.node.genesetoverrepresentation",
        "openbio.node.gsvascores",
        "openbio.node.pathwayscorettest",
        "openbio.node.rankedgsea",
    }
    assert expected.issubset(set(registered_operation_ids()))
    package = Path(__file__).parents[1] / "openbio_singlecell"
    operation_source = (package / "operations_enrichment.py").read_text(encoding="utf-8")
    node_source = (package / "nodes_enrichment.py").read_text(encoding="utf-8")
    assert "comfy_api" not in operation_source
    assert "folder_paths" not in operation_source
    assert "nodes_" not in operation_source
    assert "def execute" not in node_source


def test_dgidb_file_descriptor_runs_in_a_real_one_shot_worker(comfy_directories):
    input_dir, _, temp_dir = comfy_directories
    resource_path = input_dir / "worker-dgidb.tsv"
    _write_resource(resource_path)
    before = resource_path.read_bytes()
    gc.collect()

    async def execute():
        await initialize_artifact_service(temp_dir, sys.executable)
        return await execute_artifact_node(
            OpenBioSingleCellDGIdbAnnotation,
            None,
            {
                "resource_file": resource_path.name,
                "resource_metadata_json": _metadata_json(),
                "drug_column": "drug_claim_name",
                "gene_column": "gene_claim_name",
                "source_column": "interaction_claim_source",
                "evidence_column": "interaction_types",
                "max_file_bytes": 1_000_000,
                "max_rows": 1_000,
            },
        )

    ticket, report, code = asyncio.run(execute()).args

    assert isinstance(ticket, ArtifactTicket)
    assert ticket.kind == "OPENBIO_DGIDB_RESOURCE"
    assert ticket.codec == "dgidb-table-jsonl-v1"
    assert report.summary["node_id"] == "OpenBioSingleCellDGIdbAnnotation"
    assert "load_dgidb_resource" in code
    assert resource_path.read_bytes() == before
    read_dgidb_resource(current_artifact_runtime().resolve(ticket))
