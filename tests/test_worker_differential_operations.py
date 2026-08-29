from __future__ import annotations

import uuid
from pathlib import Path

import pytest

from openbio_singlecell.nodes_differential import DIFFERENTIAL_NODE_CLASSES
from openbio_singlecell.operations_differential import (
    PSEUDOBULK_CODEC,
    PSEUDOBULK_KIND,
    pseudobulk,
    pseudobulk_deseq2,
    pseudobulk_edger,
    scvi_differential_expression,
)
from openbio_singlecell.worker_protocol import OperationContext, ProtocolError, registered_operation_ids


def test_differential_nodes_are_schema_only():
    source = (Path(__file__).parents[1] / "openbio_singlecell" / "nodes_differential.py").read_text(
        encoding="utf-8"
    )
    assert "def execute" not in source
    assert "dependencies" not in source
    assert [node.define_schema().node_id for node in DIFFERENTIAL_NODE_CLASSES] == [
        "OpenBioSingleCellPseudobulk",
        "OpenBioSingleCellPseudobulkEdgeR",
        "OpenBioSingleCellPseudobulkDESeq2",
        "OpenBioSingleCellSCVIDifferentialExpression",
    ]


def test_differential_operations_publish_the_expected_contract_ids():
    assert PSEUDOBULK_KIND == "OPENBIO_SINGLE_CELL_PSEUDOBULK"
    assert PSEUDOBULK_CODEC == "pseudobulk-h5ad-v1"
    assert pseudobulk.__name__ == "pseudobulk"
    assert pseudobulk_edger.__name__ == "pseudobulk_edger"
    assert pseudobulk_deseq2.__name__ == "pseudobulk_deseq2"
    assert scvi_differential_expression.__name__ == "scvi_differential_expression"
    assert {
        "openbio.node.pseudobulk",
        "openbio.node.pseudobulkedger",
        "openbio.node.pseudobulkdeseq2",
        "openbio.node.scvidifferentialexpression",
    }.issubset(registered_operation_ids())


def test_scvi_worker_operation_rejects_inactive_dynamic_combo_fields_before_loading_artifacts(tmp_path):
    staging = tmp_path / "run.partial"
    staging.mkdir()
    adata_root = tmp_path / "adata"
    model_root = tmp_path / "model"
    adata_root.mkdir()
    model_root.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))
    inputs = {
        "adata": {
            "type": "artifact",
            "path": str(adata_root.resolve()),
            "kind": "OPENBIO_ANNDATA",
            "codec": "anndata-h5ad-v1",
        },
        "model": {
            "type": "artifact",
            "path": str(model_root.resolve()),
            "kind": "OPENBIO_SCVI_MODEL",
            "codec": "scvi-native-directory",
        },
    }
    parameters = {
        "groupby": "group",
        "group1": "A",
        "group2": "B",
        "population_scope": {"population_scope": "all", "subset_column": "cell_type"},
        "mode": {"mode": "vanilla"},
        "batch_handling": "observed_technical_batches",
        "n_samples_overall": 100,
        "random_seed": 0,
    }

    with pytest.raises(ProtocolError, match="fields for 'all' must be exactly"):
        scvi_differential_expression(context, inputs, parameters)
