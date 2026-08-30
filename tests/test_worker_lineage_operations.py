from __future__ import annotations

import uuid
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

from openbio_singlecell.artifact_codecs import read_anndata, write_anndata
from openbio_singlecell.cassiopeia_codec import decode_characters
from openbio_singlecell.cassiopeia_tree import _ARTIFACT_TOKEN, CassiopeiaCharacters
from openbio_singlecell.operations_lineage import cassiopeia_lineage_qc, cassiopeia_plasticity
from openbio_singlecell.worker_protocol import OperationContext


def test_lineage_qc_operation_encodes_owned_characters_and_table(tmp_path, monkeypatch):
    source = tmp_path / "alleles.tsv"
    source.write_text("cellBC\tr1\nc1\tA\n", encoding="utf-8")
    original = source.read_bytes()
    artifact = CassiopeiaCharacters(
        {"T1": pd.DataFrame([[1]], index=["c1"], columns=["r1"])},
        {"T1": {0: {1: 1.0}}},
        {"T1": {0: {1: "A"}}},
        {"source": "test"},
        _token=_ARTIFACT_TOKEN,
    )
    monkeypatch.setattr(
        "openbio_singlecell.operations_lineage.lineage_qc_owned",
        lambda _path, **_parameters: (
            artifact,
            SimpleNamespace(table=pd.DataFrame({"status": ["pass"]})),
            object(),
            "code",
        ),
    )
    monkeypatch.setattr(
        "openbio_singlecell.operations_lineage.result_metadata",
        lambda _value: {"kind": "test"},
    )
    staging = tmp_path / "staging"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))
    parameters = {
        "first_column_as_index": True,
        "lineage_column": "Tumor",
        "cell_barcode_column": "cellBC",
        "integration_barcode_column": "intBC",
        "cut_site_columns": "r1",
        "prior_grouping_columns": "Tumor,intBC",
        "missing_data_allele": "",
        "allele_representation_threshold": 0.98,
        "minimum_cells": 1,
        "maximum_missing_fraction": 0.8,
        "maximum_uncut_fraction": 0.8,
        "minimum_unique_fraction": 0.05,
        "minimum_informative_character_fraction": 0.2,
        "max_file_mib": 512,
        "max_matrix_gib": 2.0,
    }

    records = cassiopeia_lineage_qc(
        context,
        {
            "allele_table_file": {
                "type": "file",
                "path": str(source.resolve()),
                "provenance": {"path": "alleles.tsv"},
            }
        },
        parameters,
    )

    assert [record["name"] for record in records] == ["characters", "table", "summary", "code"]
    assert records[0]["codec"] == "openbio-cassiopeia-characters"
    assert records[1]["codec"] == "table-jsonl-v1"
    assert decode_characters(staging / records[0]["payload"])["lineage_ids"] == ("T1",)
    assert source.read_bytes() == original


def test_plasticity_operation_publishes_a_new_anndata_without_rewriting_input(tmp_path, monkeypatch):
    value = ad.AnnData(np.zeros((2, 1), dtype=float))
    value.obs_names = ["c1", "c2"]
    value.var_names = ["g1"]
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, value)
    original = (input_root / "data.h5ad").read_bytes()
    tree_root = tmp_path / "tree"
    tree_root.mkdir()
    captured = {}

    def fake_owned(adata, _tree, **_parameters):
        captured["adata"] = adata
        adata.obs["sc_effective_plasticity"] = [0.0, 0.5]
        return (
            adata,
            SimpleNamespace(table=pd.DataFrame({"cell_id": ["c1", "c2"]})),
            object(),
            "code",
        )

    monkeypatch.setattr("openbio_singlecell.operations_lineage.read_tree", lambda _root: object())
    monkeypatch.setattr("openbio_singlecell.operations_lineage.plasticity_owned", fake_owned)
    monkeypatch.setattr("openbio_singlecell.operations_lineage.result_metadata", lambda _value: {"kind": "test"})
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cassiopeia_plasticity(
        context,
        {
            "adata": {
                "type": "artifact",
                "path": str(input_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            },
            "tree": {
                "type": "artifact",
                "path": str(tree_root.resolve()),
                "kind": "OPENBIO_CASSIOPEIA_TREE",
                "codec": "openbio-cassiopeia-tree",
            },
        },
        {
            "annotation_key": "cell_type",
            "annotation_status": "unknown",
            "analysis_mode": "exploratory",
            "minimum_state_fraction": 0.025,
            "output_key": "sc_effective_plasticity",
            "overwrite_existing": False,
            "max_working_gib": 4.0,
        },
    )

    output = read_anndata(staging / records[0]["payload"])
    assert [record["name"] for record in records] == ["adata", "table", "summary", "code"]
    assert captured["adata"].obs["sc_effective_plasticity"].tolist() == [0.0, 0.5]
    assert output.obs["sc_effective_plasticity"].tolist() == [0.0, 0.5]
    assert (input_root / "data.h5ad").read_bytes() == original
