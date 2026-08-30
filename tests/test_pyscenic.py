from __future__ import annotations

import copy
import hashlib
import json
import socket
import subprocess
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import requests
import scipy

from openbio_singlecell import PLUGIN_VERSION
from openbio_singlecell.node_types import (
    AnnDataType,
    SCENICResultArtifactType,
    SummaryResultType,
)
from openbio_singlecell.nodes_regulatory import (
    OpenBioSingleCellImportPySCENICResults,
    OpenBioSingleCellSCENICActivityBinarization,
    OpenBioSingleCellSCENICRegulonSpecificity,
    OpenBioSingleCellSCENICTFModules,
)
from openbio_singlecell.pyscenic_import import (
    PYSCENIC_MANIFEST_SCHEMA,
    _axis_sha256,
    _matrix_sha256,
    import_pyscenic_bundle,
    pyscenic_import_code,
)
from openbio_singlecell.scenic_artifact import SCENICResultArtifact
from openbio_singlecell.scenic_binarization import (
    HDT_DIP_SIMULATIONS,
    _derive_threshold_exact,
    binarize_scenic_activity,
    scenic_binarization_code,
)
from openbio_singlecell.scenic_membership import scenic_membership_code, scenic_regulon_membership
from openbio_singlecell.scenic_rss import compute_scenic_rss, scenic_rss_code


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_regulons(path: Path, *, context: str | None = None) -> None:
    columns = pd.MultiIndex.from_product(
        [
            ["Enrichment"],
            [
                "AUC",
                "NES",
                "MotifSimilarityQvalue",
                "OrthologousIdentity",
                "Annotation",
                "Context",
                "TargetGenes",
                "RankAtMax",
            ],
        ]
    )
    index = pd.MultiIndex.from_tuples(
        [("TF1", "M1"), ("TF1", "M2"), ("TF2", "M3")], names=["TF", "MotifID"]
    )
    frame = pd.DataFrame(
        [
            [
                0.20,
                3.2,
                0.0,
                1.0,
                "direct",
                context or "frozenset({'activating', 'weight>75.0%'})",
                "[('G1', 0.8), ('G2', 0.4)]",
                100,
            ],
            [
                0.18,
                3.0,
                0.0,
                1.0,
                "orthology",
                "frozenset({'activating', 'weight>75.0%'})",
                "[('G1', 0.9)]",
                90,
            ],
            [
                0.22,
                3.5,
                0.0,
                1.0,
                "direct",
                "frozenset({'repressing', 'weight>75.0%'})",
                "[('G3', 0.7), ('TF1', 0.2)]",
                110,
            ],
        ],
        index=index,
        columns=columns,
    )
    frame.to_csv(path)


def _refresh_manifest_hashes(manifest_path: Path) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    base = manifest_path.parent
    manifest["expression"]["sha256"] = _sha256(base / manifest["expression"]["file"])
    for declaration in manifest["artifacts"].values():
        declaration["sha256"] = _sha256(base / declaration["file"])
    for declaration in manifest["resources"]:
        declaration["sha256"] = _sha256(base / declaration["file"])
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, allow_nan=False, indent=2), encoding="utf-8"
    )


@pytest.fixture
def pyscenic_bundle(tmp_path, science):
    bundle = tmp_path / "pyscenic-bundle"
    bundle.mkdir()
    obs_names = [f"cell_{index}" for index in range(8)]
    gene_names = ["TF1", "G1", "G2", "TF2", "G3"]
    counts = science.np.asarray(
        [
            [5, 1, 0, 2, 4],
            [4, 2, 0, 1, 5],
            [6, 1, 1, 1, 4],
            [5, 0, 2, 2, 3],
            [2, 4, 5, 5, 1],
            [1, 5, 4, 6, 0],
            [2, 6, 5, 5, 1],
            [1, 4, 6, 7, 0],
        ],
        dtype=float,
    )
    adata = science.ad.AnnData(
        counts,
        obs=science.pd.DataFrame(
            {
                "cell_type": science.pd.Categorical(
                    ["A"] * 4 + ["B"] * 4, categories=["A", "B", "unused"]
                )
            },
            index=obs_names,
        ),
        var=science.pd.DataFrame(index=gene_names),
    )

    expression = bundle / "expression.csv"
    science.pd.DataFrame(counts, index=obs_names, columns=gene_names).to_csv(expression)
    adjacency = bundle / "adjacency.csv"
    science.pd.DataFrame(
        {
            "TF": ["TF1", "TF1", "TF2"],
            "target": ["G1", "G2", "G3"],
            "importance": [0.9, 0.4, 0.7],
        }
    ).to_csv(adjacency, index=False)
    regulons = bundle / "regulons.csv"
    _write_regulons(regulons)
    aucell = bundle / "aucell.csv"
    science.pd.DataFrame(
        {
            "TF1(+)": [0.02, 0.04, 0.03, 0.05, 0.70, 0.80, 0.75, 0.90],
            "TF2(-)": [0.75, 0.80, 0.70, 0.85, 0.08, 0.05, 0.10, 0.04],
        },
        index=obs_names,
    ).to_csv(aucell)
    tf_list = bundle / "tfs.txt"
    tf_list.write_text("TF1\nTF2\n", encoding="utf-8")
    ranking = bundle / "ranking.feather"
    ranking.write_bytes(b"ARROW1" + b"audited-feather-v2-fixture" + b"ARROW1")
    motifs = bundle / "motifs.tbl"
    motifs.write_text(
        "#motif_id\tgene_name\tmotif_similarity_qvalue\torthologous_identity\tdescription\n"
        "M1\tTF1\t0.0\t1.0\tdirect\n"
        "M2\tTF1\t0.0\t1.0\torthology\n"
        "M3\tTF2\t0.0\t1.0\tdirect\n",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": PYSCENIC_MANIFEST_SCHEMA,
        "run_id": "external-run-001",
        "organism": "Homo sapiens",
        "genome_build": "GRCh38",
        "gene_namespace": "HGNC-symbol",
        "expression_state": "post_qc_raw_counts",
        "container": {
            "image": "aertslab/pyscenic:0.12.1",
            "digest": f"sha256:{'1' * 64}",
        },
        "software_versions": {
            "pyscenic": "0.12.1",
            "ctxcore": "0.2.0",
            "arboreto": "0.1.6",
            "python": "3.10.13",
        },
        "commands": {
            "grn": [
                "pyscenic",
                "grn",
                "expression.csv",
                "tfs.txt",
                "-o",
                "adjacency.csv",
                "--method",
                "grnboost2",
                "--num_workers",
                "1",
                "--seed",
                "7",
            ],
            "ctx": [
                "pyscenic",
                "ctx",
                "adjacency.csv",
                "ranking.feather",
                "--annotations_fname",
                "motifs.tbl",
                "--expression_mtx_fname",
                "expression.csv",
                "-o",
                "regulons.csv",
                "--num_workers",
                "1",
                "--mode",
                "custom_multiprocessing",
                "--mask_dropouts",
            ],
            "aucell": [
                "pyscenic",
                "aucell",
                "expression.csv",
                "regulons.csv",
                "-o",
                "aucell.csv",
                "--auc_threshold",
                "0.05",
                "--num_workers",
                "1",
                "--seed",
                "7",
            ],
        },
        "expression": {
            "file": "expression.csv",
            "sha256": _sha256(expression),
            "cells": len(obs_names),
            "genes": len(gene_names),
            "observation_ids_sha256": _axis_sha256(obs_names),
            "gene_ids_sha256": _axis_sha256(gene_names),
            "matrix_sha256": _matrix_sha256(
                counts,
                rows=len(obs_names),
                columns=len(gene_names),
                numpy=science.np,
                sparse=science.sparse,
            ),
        },
        "resources": [
            {
                "role": "tf_list",
                "file": "tfs.txt",
                "sha256": _sha256(tf_list),
                "release": "SCENIC protocol fixture",
                "organism": "Homo sapiens",
                "genome_build": "GRCh38",
                "gene_namespace": "HGNC-symbol",
                "license": "fixture-only",
                "citation": "Van de Sande et al. 2020",
            },
            {
                "role": "ranking_database",
                "file": "ranking.feather",
                "sha256": _sha256(ranking),
                "release": "cisTarget Feather-v2 fixture",
                "organism": "Homo sapiens",
                "genome_build": "GRCh38",
                "gene_namespace": "HGNC-symbol",
                "license": "fixture-only",
                "citation": "Van de Sande et al. 2020",
            },
            {
                "role": "motif_annotations",
                "file": "motifs.tbl",
                "sha256": _sha256(motifs),
                "release": "motif2TF fixture",
                "organism": "Homo sapiens",
                "genome_build": "GRCh38",
                "gene_namespace": "HGNC-symbol",
                "license": "fixture-only",
                "citation": "Van de Sande et al. 2020",
            },
        ],
        "artifacts": {
            "adjacency": {"file": "adjacency.csv", "sha256": _sha256(adjacency)},
            "regulons": {"file": "regulons.csv", "sha256": _sha256(regulons)},
            "aucell": {"file": "aucell.csv", "sha256": _sha256(aucell)},
        },
    }
    manifest_path = bundle / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, allow_nan=False, indent=2), encoding="utf-8"
    )
    return adata, manifest_path


def test_pyscenic_node_schemas_are_atomic_and_typed():
    importer = OpenBioSingleCellImportPySCENICResults.define_schema()
    rss = OpenBioSingleCellSCENICRegulonSpecificity.define_schema()
    binary = OpenBioSingleCellSCENICActivityBinarization.define_schema()
    membership = OpenBioSingleCellSCENICTFModules.define_schema()

    assert [item.id for item in importer.inputs] == [
        "adata",
        "run_manifest_json",
        "overwrite",
        "max_file_bytes",
        "max_adjacency_edges",
        "max_regulon_edges",
        "max_dense_bytes",
    ]
    assert [(item.display_name, item.io_type) for item in importer.outputs] == [
        ("adata", AnnDataType.io_type),
        ("scenic_result", SCENICResultArtifactType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert rss.inputs[1].get_io_type() == SCENICResultArtifactType.io_type
    assert binary.inputs[0].get_io_type() == SCENICResultArtifactType.io_type
    assert [(item.display_name, item.io_type) for item in binary.outputs] == [
        ("thresholds", "OPENBIO_SINGLE_CELL_TABLE"),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert all(item.io_type != "OPENBIO_SCENIC_BINARY" for item in binary.outputs)
    assert membership.inputs[0].get_io_type() == SCENICResultArtifactType.io_type
    assert all(item.io_type != "OPENBIO_SCENIC_NETWORK" for item in membership.outputs)


def test_import_validates_bundle_owns_outputs_and_emits_strict_report(pyscenic_bundle, science):
    adata, manifest_path = pyscenic_bundle
    original = adata.copy()
    output, artifact, summary = import_pyscenic_bundle(adata, manifest_path)

    science.np.testing.assert_array_equal(adata.X, original.X)
    assert output is adata
    assert output.obs_names.tolist() == adata.obs_names.tolist()
    assert output.var_names.tolist() == adata.var_names.tolist()
    assert output.obsm["scenic_auc"].columns.tolist() == ["TF1(+)", "TF2(-)"]
    assert artifact.regulon_names == ("TF1(+)", "TF2(-)")
    assert artifact.membership.to_dict(orient="records") == [
        {
            "regulon": "TF1(+)",
            "transcription_factor": "TF1",
            "regulation": "activating",
            "context": "weight>75.0%",
            "target": "G1",
            "target_weight": 0.9,
            "motif_evidence_count": 2,
            "motif_ids": ["M1", "M2"],
            "target_rank": 1,
        },
        {
            "regulon": "TF1(+)",
            "transcription_factor": "TF1",
            "regulation": "activating",
            "context": "weight>75.0%",
            "target": "G2",
            "target_weight": 0.4,
            "motif_evidence_count": 1,
            "motif_ids": ["M1"],
            "target_rank": 2,
        },
        {
            "regulon": "TF2(-)",
            "transcription_factor": "TF2",
            "regulation": "repressing",
            "context": "weight>75.0%",
            "target": "G3",
            "target_weight": 0.7,
            "motif_evidence_count": 1,
            "motif_ids": ["M3"],
            "target_rank": 1,
        },
        {
            "regulon": "TF2(-)",
            "transcription_factor": "TF2",
            "regulation": "repressing",
            "context": "weight>75.0%",
            "target": "TF1",
            "target_weight": 0.2,
            "motif_evidence_count": 1,
            "motif_ids": ["M3"],
            "target_rank": 2,
        },
    ]
    activity_copy = artifact.activity
    membership_copy = artifact.membership
    provenance_copy = artifact.provenance
    activity_copy.iloc[0, 0] = 1.0
    membership_copy.loc[0, "target"] = "changed"
    provenance_copy["run_id"] = "changed"
    assert artifact.activity.iloc[0, 0] == 0.02
    assert artifact.membership.loc[0, "target"] == "G1"
    assert artifact.provenance["run_id"] == "external-run-001"
    assert summary["key_results"]["validation_accounting"]["adjacency"]["edges"] == 3
    assert summary["software_versions"]["numpy"] != "not-installed"
    assert {reference["doi"] for reference in summary["references"]} >= {
        "10.1038/nmeth.4463",
        "10.1038/s41596-020-0336-2",
    }
    json.dumps(summary, allow_nan=False)


def test_scenic_file_codec_round_trip_uses_h5ad_jsonl_and_json(pyscenic_bundle, science, tmp_path):
    from openbio_singlecell.regulatory_artifact_codecs import read_scenic, write_scenic
    from openbio_singlecell.scenic_artifact import validate_scenic_result_artifact

    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    root = tmp_path / "scenic"
    root.mkdir()

    write_scenic(root, artifact)
    portable = read_scenic(root)
    activity, membership, provenance, metadata = validate_scenic_result_artifact(
        portable,
        exact_type=False,
        numpy=science.np,
        pandas=science.pd,
        copy_result=False,
    )

    science.pd.testing.assert_frame_equal(activity, artifact.activity)
    science.pd.testing.assert_frame_equal(membership, artifact.membership)
    assert provenance == artifact.provenance
    assert metadata["artifact_fingerprint_sha256"] == artifact.artifact_fingerprint_sha256
    assert (root / "activity" / "data.h5ad").is_file()
    assert (root / "membership" / "data.jsonl").is_file()
    assert (root / "result.json").is_file()


def test_pyscenic_worker_operation_writes_new_artifacts_without_rewriting_input(
    pyscenic_bundle,
    tmp_path,
):
    import uuid

    from openbio_singlecell.artifact_codecs import read_anndata, write_anndata
    from openbio_singlecell.operations_regulatory import import_pyscenic_results
    from openbio_singlecell.regulatory_artifact_codecs import read_scenic
    from openbio_singlecell.worker_protocol import OperationContext

    adata, manifest_path = pyscenic_bundle
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_bytes = (input_root / "data.h5ad").read_bytes()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = import_pyscenic_results(
        context,
        {
            "adata": {
                "type": "artifact",
                "path": str(input_root.resolve()),
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
            },
            "run_manifest_json": {
                "type": "file",
                "path": str(manifest_path.resolve()),
                "provenance": {},
            },
        },
        {
            "overwrite": False,
            "max_file_bytes": 2_147_483_647,
            "max_adjacency_edges": 100,
            "max_regulon_edges": 100,
            "max_dense_bytes": 1_073_741_824,
        },
    )

    assert [record["name"] for record in records] == ["adata", "scenic_result", "summary", "code"]
    assert (input_root / "data.h5ad").read_bytes() == input_bytes
    assert "scenic_auc" in read_anndata(staging / records[0]["payload"]).obsm
    portable = read_scenic(staging / records[1]["payload"])
    assert portable["regulons"] == ["TF1(+)", "TF2(-)"]


def test_import_generated_code_is_full_runtime_equivalent(pyscenic_bundle, science):
    adata, manifest_path = pyscenic_bundle
    generated_input = adata.copy()
    runtime_output, runtime_artifact, runtime_summary = import_pyscenic_bundle(
        adata, manifest_path, max_adjacency_edges=100, max_regulon_edges=100
    )
    code = pyscenic_import_code(
        manifest_path=str(manifest_path),
        overwrite=False,
        max_file_bytes=2_147_483_647,
        max_adjacency_edges=100,
        max_regulon_edges=100,
        max_dense_bytes=1_073_741_824,
    )
    compile(code, "<generated-pyscenic-import>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    generated_output, generated_artifact, generated_summary = namespace[
        "import_validated_pyscenic_bundle"
    ](generated_input)

    science.pd.testing.assert_frame_equal(
        runtime_output.obsm["scenic_auc"], generated_output.obsm["scenic_auc"]
    )
    assert runtime_artifact.to_portable() == generated_artifact.to_portable()
    assert runtime_summary == generated_summary


@pytest.mark.parametrize(
    ("member", "message"),
    [
        ("adjacency.csv", "SHA-256 mismatch"),
        ("ranking.feather", "SHA-256 mismatch"),
        ("motifs.tbl", "SHA-256 mismatch"),
    ],
)
def test_import_rejects_any_member_hash_mismatch(pyscenic_bundle, member, message):
    adata, manifest_path = pyscenic_bundle
    path = manifest_path.parent / member
    path.write_bytes(path.read_bytes() + b"tamper")
    with pytest.raises(ValueError, match=message):
        import_pyscenic_bundle(adata, manifest_path)


def test_import_rejects_malicious_literals_duplicate_edges_and_unbound_commands(
    pyscenic_bundle, science
):
    adata, manifest_path = pyscenic_bundle
    _write_regulons(
        manifest_path.parent / "regulons.csv",
        context="frozenset(__import__('os').system('echo unsafe'))",
    )
    _refresh_manifest_hashes(manifest_path)
    with pytest.raises(ValueError, match="safe set literal"):
        import_pyscenic_bundle(adata, manifest_path)

    adata, manifest_path = pyscenic_bundle
    adjacency = science.pd.read_csv(manifest_path.parent / "adjacency.csv")
    science.pd.concat([adjacency, adjacency.iloc[[0]]]).to_csv(
        manifest_path.parent / "adjacency.csv", index=False
    )
    _refresh_manifest_hashes(manifest_path)
    with pytest.raises(ValueError, match="duplicate TF-target"):
        import_pyscenic_bundle(adata, manifest_path)

    adata, manifest_path = pyscenic_bundle
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["commands"]["grn"][2] = "other-expression.csv"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="positional inputs"):
        import_pyscenic_bundle(adata, manifest_path)


def test_import_fails_closed_for_axes_identity_and_overwrite_but_ignores_unrelated_files(
    pyscenic_bundle, science
):
    adata, manifest_path = pyscenic_bundle
    reordered = adata[adata.obs_names[::-1], :].copy()
    with pytest.raises(ValueError, match="observation identity fingerprint"):
        import_pyscenic_bundle(reordered, manifest_path)

    changed = adata.copy()
    changed.X[0, 0] += 1
    with pytest.raises(ValueError, match="expression-matrix fingerprint"):
        import_pyscenic_bundle(changed, manifest_path)

    (manifest_path.parent / "undeclared.txt").write_text("extra", encoding="utf-8")
    _extra_output, _extra_artifact, extra_summary = import_pyscenic_bundle(adata, manifest_path)
    assert any("unrelated files" in warning for warning in extra_summary["warnings"])
    (manifest_path.parent / "undeclared.txt").unlink()

    existing = adata.copy()
    existing.obsm["scenic_auc"] = science.pd.DataFrame(
        science.np.ones((adata.n_obs, 1)), index=adata.obs_names, columns=["existing"]
    )
    before = existing.obsm["scenic_auc"].copy(deep=True)
    with pytest.raises(ValueError, match="already exists"):
        import_pyscenic_bundle(existing, manifest_path)
    science.pd.testing.assert_frame_equal(existing.obsm["scenic_auc"], before)


def test_import_accepts_explicit_numeric_state_optional_container_and_cli_defaults(
    pyscenic_bundle, science
):
    adata, manifest_path = pyscenic_bundle
    transformed = science.np.log1p(science.np.asarray(adata.X, dtype=float)) - 1.0
    adata.X = transformed.copy()
    science.pd.DataFrame(
        transformed,
        index=adata.obs_names,
        columns=adata.var_names,
    ).to_csv(manifest_path.parent / "expression.csv", index_label="expert_cell_id")
    aucell = science.pd.read_csv(manifest_path.parent / "aucell.csv", index_col=0)
    aucell.to_csv(manifest_path.parent / "aucell.csv", index_label="expert_activity_cell")
    regulon_path = manifest_path.parent / "regulons.csv"
    regulon_path.write_text(
        regulon_path.read_text(encoding="utf-8").replace("('G2', 0.4)", "('G2', 0.0)"),
        encoding="utf-8",
    )
    adjacency = science.pd.read_csv(manifest_path.parent / "adjacency.csv")
    adjacency.loc[0, "importance"] = 0.0
    adjacency.to_csv(manifest_path.parent / "adjacency.csv", index=False)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["expression_state"] = "expert_scaled_finite_expression"
    manifest["container"] = None
    manifest["software_versions"].update({"numpy": "1.23.5", "loompy": "3.0.8"})
    manifest["expression"]["matrix_sha256"] = _matrix_sha256(
        transformed,
        rows=adata.n_obs,
        columns=adata.n_vars,
        numpy=science.np,
        sparse=science.sparse,
    )

    def omit_options(command, *flags):
        for flag in flags:
            position = command.index(flag)
            del command[position : position + 2]

    omit_options(manifest["commands"]["grn"], "--method", "--num_workers", "--seed")
    omit_options(manifest["commands"]["ctx"], "--num_workers", "--mode")
    omit_options(manifest["commands"]["aucell"], "--auc_threshold", "--num_workers", "--seed")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _refresh_manifest_hashes(manifest_path)
    generated_input = adata.copy()

    output, artifact, summary = import_pyscenic_bundle(adata, manifest_path)
    science.np.testing.assert_array_equal(output.X, transformed)
    assert artifact.provenance["container"] is None
    assert artifact.provenance["expression_state"] == "expert_scaled_finite_expression"
    assert artifact.provenance["declared_software_versions"]["loompy"] == "3.0.8"
    assert summary["parameters"]["commands"]["grn"] == {
        "num_workers": None,
        "seed": None,
        "method": "grnboost2",
    }
    assert summary["parameters"]["commands"]["ctx"]["mode"] == "custom_multiprocessing"
    assert summary["parameters"]["commands"]["aucell"]["auc_threshold"] == 0.05
    assert summary["key_results"]["validation_accounting"]["adjacency"][
        "zero_importance_edges"
    ] == 1
    assert summary["key_results"]["validation_accounting"]["regulons"][
        "zero_weight_edges"
    ] == 1
    assert any("explicit expert choice" in warning for warning in summary["warnings"])
    assert any("No immutable container" in warning for warning in summary["warnings"])
    assert any("omitted --seed" in warning for warning in summary["warnings"])

    code = pyscenic_import_code(
        manifest_path=str(manifest_path),
        overwrite=False,
        max_file_bytes=2_147_483_647,
        max_adjacency_edges=10_000_000,
        max_regulon_edges=2_000_000,
        max_dense_bytes=1_073_741_824,
    )
    namespace: dict[str, object] = {}
    exec(compile(code, "<generated-expert-pyscenic-import>", "exec"), namespace)
    generated_output, generated_artifact, generated_summary = namespace[
        "import_validated_pyscenic_bundle"
    ](generated_input)
    science.np.testing.assert_array_equal(output.X, generated_output.X)
    assert artifact.to_portable() == generated_artifact.to_portable()
    assert summary == generated_summary


def test_import_rejects_path_traversal_and_symlink_members(pyscenic_bundle, tmp_path):
    adata, manifest_path = pyscenic_bundle
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["resources"][0]["file"] = "../tfs.txt"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="normalized relative path|traversal"):
        import_pyscenic_bundle(adata, manifest_path)

    manifest["resources"][0]["file"] = "tfs.txt"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    tf_path = manifest_path.parent / "tfs.txt"
    outside = tmp_path / "outside-tfs.txt"
    outside.write_bytes(tf_path.read_bytes())
    tf_path.unlink()
    try:
        tf_path.symlink_to(outside)
    except OSError:
        pytest.skip("This Windows environment does not permit symlink creation.")
    with pytest.raises(ValueError, match="symlink"):
        import_pyscenic_bundle(adata, manifest_path)


def test_import_never_uses_network_shell_or_temporary_execution(pyscenic_bundle, monkeypatch):
    adata, manifest_path = pyscenic_bundle

    def forbidden(*args, **kwargs):
        del args, kwargs
        raise AssertionError("safe pySCENIC import attempted an external side effect")

    monkeypatch.setattr(subprocess, "run", forbidden)
    monkeypatch.setattr(tempfile, "TemporaryDirectory", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setattr(urllib.request, "urlopen", forbidden)
    monkeypatch.setattr(requests.sessions.Session, "request", forbidden)
    output, artifact, summary = import_pyscenic_bundle(adata, manifest_path)
    assert output.obsm["scenic_auc"].shape == (8, 2)
    assert artifact.regulon_names == ("TF1(+)", "TF2(-)")
    assert summary["status"] == "validated_external_pyscenic_regulatory_evidence"


def test_rss_complete_grid_formula_input_immutability_and_generated_parity(
    pyscenic_bundle, science
):
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    artifact_before = artifact.to_portable()
    table, summary = compute_scenic_rss(
        adata,
        artifact,
        annotation_key="cell_type",
        annotation_status="provisional",
        max_output_rows=20,
    )
    assert len(table) == 4
    assert set(table["group"]) == {"A", "B"}
    assert table.groupby("group")["regulon"].nunique().to_dict() == {"A": 2, "B": 2}
    assert summary["key_results"]["unused_categorical_levels"] == ["unused"]
    assert artifact.to_portable() == artifact_before
    values = artifact.activity["TF1(+)"].to_numpy()
    indicator = science.np.asarray([1.0] * 4 + [0.0] * 4)
    expected = 1.0 - scipy.spatial.distance.jensenshannon(
        values / values.sum(), indicator / indicator.sum()
    )
    observed = table.loc[(table["group"] == "A") & (table["regulon"] == "TF1(+)"), "rss"].item()
    assert observed == pytest.approx(expected, abs=1e-15)

    code = scenic_rss_code(
        parameters={
            "annotation_key": "cell_type",
            "annotation_status": "provisional",
            "max_output_rows": 20,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    namespace: dict[str, object] = {}
    exec(compile(code, "<generated-scenic-rss>", "exec"), namespace)
    generated_table, generated_summary = namespace["scenic_regulon_specificity"](
        adata, artifact.to_portable()
    )
    science.pd.testing.assert_frame_equal(table, generated_table)
    assert summary == generated_summary


def test_rss_warns_for_one_group_and_rejects_bad_annotations_and_row_budget(pyscenic_bundle):
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    with pytest.raises(ValueError, match="complete grid"):
        compute_scenic_rss(adata, artifact, max_output_rows=3)
    one_group = adata.copy()
    one_group.obs["cell_type"] = "A"
    one_group_table, one_group_summary = compute_scenic_rss(one_group, artifact)
    assert one_group_table["group"].unique().tolist() == ["A"]
    assert any("no between-group specificity" in warning for warning in one_group_summary["warnings"])
    missing = adata.copy()
    missing.obs.loc[missing.obs_names[0], "cell_type"] = None
    with pytest.raises(ValueError, match="missing labels"):
        compute_scenic_rss(missing, artifact)


def test_rss_matches_real_pyscenic_0121_backend(pyscenic_bundle, monkeypatch):
    pyscenic = pytest.importorskip("pyscenic")
    assert pyscenic.__version__ == "0.12.1"
    backend = pytest.importorskip("pyscenic.rss")
    monkeypatch.setitem(np.__dict__, "float", float)
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    table, _summary = compute_scenic_rss(adata, artifact, max_output_rows=20)
    expected = backend.regulon_specificity_scores(artifact.activity, adata.obs["cell_type"])
    for group in expected.index:
        for regulon in expected.columns:
            observed = table.loc[
                (table["group"] == group) & (table["regulon"] == regulon), "rss"
            ].item()
            assert observed == expected.loc[group, regulon]


def test_membership_is_bounded_schema_stable_and_generated_equivalent(pyscenic_bundle, science):
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    table, summary = scenic_regulon_membership(
        artifact, transcription_factor="TF1", max_output_rows=10
    )
    assert table["transcription_factor"].tolist() == ["TF1", "TF1"]
    assert summary["key_results"]["complete_rows"] == 4
    empty, empty_summary = scenic_regulon_membership(
        artifact, transcription_factor="missing", max_output_rows=10
    )
    assert empty.empty and empty.columns.tolist() == table.columns.tolist()
    assert any("schema-stable empty" in value for value in empty_summary["warnings"])
    with pytest.raises(ValueError, match="exceeding max_output_rows"):
        scenic_regulon_membership(artifact, max_output_rows=3)

    code = scenic_membership_code(
        parameters={
            "transcription_factor": "TF1",
            "max_output_rows": 10,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    namespace: dict[str, object] = {}
    exec(compile(code, "<generated-scenic-membership>", "exec"), namespace)
    generated_table, generated_summary = namespace["scenic_final_regulon_membership"](
        artifact.to_portable()
    )
    science.pd.testing.assert_frame_equal(table, generated_table)
    assert summary == generated_summary


def test_binarization_overrides_strict_comparison_rng_and_generated_parity(
    pyscenic_bundle, science
):
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    overrides = science.pd.DataFrame(
        {"regulon": ["TF1(+)", "TF2(-)"], "threshold": [0.70, 0.70]}
    )
    overrides_before = overrides.copy(deep=True)
    state_before = copy.deepcopy(science.np.random.get_state())
    binary, thresholds, summary = binarize_scenic_activity(
        artifact,
        random_seed=11,
        threshold_overrides=overrides,
        max_dense_bytes=10_000,
    )
    science.pd.testing.assert_frame_equal(overrides, overrides_before)
    state_after = science.np.random.get_state()
    assert state_before[0] == state_after[0]
    science.np.testing.assert_array_equal(state_before[1], state_after[1])
    assert state_before[2:] == state_after[2:]
    assert bool(binary.binary.loc["cell_4", "TF1(+)"]) is False
    assert bool(binary.binary.loc["cell_5", "TF1(+)"]) is True
    assert thresholds["threshold_source"].tolist() == ["manual_override", "manual_override"]
    assert summary["key_results"]["estimated_peak_dense_bytes"] == adata.n_obs * 2 * 96
    assert summary["parameters"]["dip_test_simulations"] == HDT_DIP_SIMULATIONS

    code = scenic_binarization_code(
        parameters={
            "random_seed": 11,
            "max_dense_bytes": 10_000,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    assert "GPL-3.0-or-later" in code
    assert "ce41b61b6570490949bd12b514e9f6de46d19c1f" in code
    assert "a0e3d448a4b266f54ec63a5b3d5be351fbd1db1c" in code
    namespace: dict[str, object] = {}
    exec(compile(code, "<generated-scenic-binary>", "exec"), namespace)
    generated_binary, generated_thresholds, generated_summary = namespace[
        "binarize_scenic_regulon_activity"
    ](artifact.to_portable(), overrides)
    assert binary.to_portable() == generated_binary.to_portable()
    science.pd.testing.assert_frame_equal(thresholds, generated_thresholds)
    assert summary == generated_summary

    derived_binary, derived_thresholds, derived_summary = binarize_scenic_activity(
        artifact,
        random_seed=11,
        max_dense_bytes=1_000_000,
    )
    assert derived_summary["key_results"]["estimated_peak_dense_bytes"] == (
        adata.n_obs * 2 * 96 + adata.n_obs * HDT_DIP_SIMULATIONS * 16
    )
    derived_code = scenic_binarization_code(
        parameters={
            "random_seed": 11,
            "max_dense_bytes": 1_000_000,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    derived_namespace: dict[str, object] = {}
    exec(compile(derived_code, "<generated-scenic-derived-binary>", "exec"), derived_namespace)
    generated_derived_binary, generated_derived_thresholds, generated_derived_summary = (
        derived_namespace["binarize_scenic_regulon_activity"](artifact.to_portable())
    )
    assert derived_binary.to_portable() == generated_derived_binary.to_portable()
    science.pd.testing.assert_frame_equal(derived_thresholds, generated_derived_thresholds)
    assert derived_summary == generated_derived_summary


def test_binarization_rejects_bad_seed_memory_and_adversarial_backend(pyscenic_bundle):
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    with pytest.raises(ValueError, match="random_seed"):
        binarize_scenic_activity(artifact, random_seed=0)
    with pytest.raises(ValueError, match="dense working bytes"):
        binarize_scenic_activity(artifact, max_dense_bytes=1)
    partial_override = pd.DataFrame(
        {"regulon": ["TF1(+)"], "threshold": [0.5]}
    )
    with pytest.raises(ValueError, match="dense working bytes"):
        binarize_scenic_activity(
            artifact,
            threshold_overrides=partial_override,
            max_dense_bytes=10_000,
        )

    expert_overrides = pd.DataFrame(
        {"regulon": ["TF1(+)", "TF2(-)"], "threshold": [-0.1, 1.1]}
    )
    expert_binary, _expert_thresholds, expert_summary = binarize_scenic_activity(
        artifact,
        threshold_overrides=expert_overrides,
        max_dense_bytes=10_000,
    )
    assert bool(expert_binary.binary["TF1(+)"].all()) is True
    assert bool(expert_binary.binary["TF2(-)"].any()) is False
    assert expert_summary["key_results"]["out_of_range_override_regulons"] == [
        "TF1(+)",
        "TF2(-)",
    ]
    assert any("explicit expert choices" in warning for warning in expert_summary["warnings"])

    zero_activity = artifact.activity
    zero_activity["TF1(+)"] = 0.0
    zero_artifact = SCENICResultArtifact(zero_activity, artifact.membership, artifact.provenance)
    zero_override = pd.DataFrame({"regulon": ["TF2(-)"], "threshold": [0.5]})
    zero_binary, zero_thresholds, zero_summary = binarize_scenic_activity(
        zero_artifact,
        threshold_overrides=zero_override,
        max_dense_bytes=1_000_000,
    )
    zero_row = zero_thresholds.loc[zero_thresholds["regulon"] == "TF1(+)"]
    assert zero_row["threshold"].item() == 0.0
    assert zero_row["active_cells"].item() == 0
    assert bool(zero_binary.binary["TF1(+)"].any()) is False
    assert "TF1(+)" in zero_summary["key_results"]["low_information_regulons"]

    def mutating_backend(frame, regulon, seed):
        del regulon, seed
        frame.iloc[0, 0] = 1.0
        return 0.5

    with pytest.raises(RuntimeError, match="modified its private activity"):
        binarize_scenic_activity(artifact, _threshold_backend=mutating_backend)

    def nonfinite_backend(frame, regulon, seed):
        del frame, regulon, seed
        return float("nan")

    with pytest.raises(RuntimeError, match="non-finite"):
        binarize_scenic_activity(artifact, _threshold_backend=nonfinite_backend)


def test_worker_owned_scenic_consumers_reuse_decoded_activity(
    pyscenic_bundle,
    monkeypatch,
):
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    portable = artifact.to_portable()
    portable["activity"] = pd.DataFrame(
        portable["activity"],
        index=portable["observations"],
        columns=portable["regulons"],
    )
    activity = portable["activity"]
    copied_shapes = []
    copy_true_shapes = []
    original_copy = pd.DataFrame.copy
    original_to_numpy = pd.DataFrame.to_numpy

    def tracked_copy(frame, *args, **kwargs):
        if frame is activity:
            copied_shapes.append(frame.shape)
        return original_copy(frame, *args, **kwargs)

    def tracked_to_numpy(frame, *args, **kwargs):
        if frame is activity and kwargs.get("copy") is True:
            copy_true_shapes.append(frame.shape)
        return original_to_numpy(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "copy", tracked_copy)
    monkeypatch.setattr(pd.DataFrame, "to_numpy", tracked_to_numpy)
    compute_scenic_rss(
        adata,
        portable,
        max_output_rows=20,
        _portable_artifact=True,
        _worker_owned=True,
    )
    overrides = pd.DataFrame(
        {"regulon": portable["regulons"], "threshold": [0.5] * len(portable["regulons"])}
    )
    binary, _thresholds, binary_summary = binarize_scenic_activity(
        portable,
        threshold_overrides=overrides,
        max_dense_bytes=10_000,
        _portable_artifact=True,
        _worker_owned=True,
    )

    assert binary is None
    assert len(binary_summary["key_results"]["binary_artifact_fingerprint_sha256"]) == 64
    assert copied_shapes == []
    assert copy_true_shapes == []


def test_binarization_exact_transcription_matches_real_pyscenic_0121(monkeypatch):
    pyscenic = pytest.importorskip("pyscenic")
    assert pyscenic.__version__ == "0.12.1"
    backend = pytest.importorskip("pyscenic.binarization")
    scipy = pytest.importorskip("scipy")
    sklearn = pytest.importorskip("sklearn")
    monkeypatch.setattr(np, "msort", np.sort, raising=False)
    monkeypatch.setattr(np, "float", float, raising=False)
    values = np.concatenate([np.linspace(0.01, 0.20, 40), np.linspace(0.70, 0.95, 40)])
    frame = pd.DataFrame({"TF1(+)": values})
    expected = backend.derive_threshold(frame, "TF1(+)", seed=17, method="hdt")
    observed = _derive_threshold_exact(
        frame,
        "TF1(+)",
        17,
        numpy=np,
        scipy_stats=scipy.stats,
        minimize_scalar=scipy.optimize.minimize_scalar,
        mixture=sklearn.mixture,
    )
    assert observed == expected


def test_manifest_commands_parse_with_real_pyscenic_0121_cli(pyscenic_bundle, monkeypatch):
    pyscenic = pytest.importorskip("pyscenic")
    assert pyscenic.__version__ == "0.12.1"
    # Tagged 0.12.1 imports its CLI through transform.py, which still names this removed alias.
    # The test-only alias lets us inspect the real parser; production deliberately never imports this module.
    monkeypatch.setitem(np.__dict__, "object", object)
    cli = pytest.importorskip("pyscenic.cli.pyscenic")
    _adata, manifest_path = pyscenic_bundle
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    monkeypatch.chdir(manifest_path.parent)
    parser = cli.create_argument_parser()
    parsed = {}
    for stage, command in manifest["commands"].items():
        arguments = parser.parse_args(command[1:])
        parsed[stage] = arguments
        for value in vars(arguments).values():
            values = value if isinstance(value, (list, tuple)) else [value]
            for item in values:
                if hasattr(item, "close") and not getattr(item, "closed", True):
                    item.close()
    assert parsed["grn"].method == "grnboost2"
    assert parsed["grn"].seed == 7
    assert parsed["ctx"].mode == "custom_multiprocessing"
    assert parsed["ctx"].mask_dropouts is True
    assert parsed["aucell"].seed == 7
    assert parsed["aucell"].auc_threshold == pytest.approx(0.05)


def test_typed_artifact_rejects_forged_provenance_and_payload(pyscenic_bundle):
    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    provenance = artifact.provenance
    provenance["declared_software_versions"]["pyscenic"] = "0.12.0"
    with pytest.raises(ValueError, match="exactly 0.12.1"):
        SCENICResultArtifact(artifact.activity, artifact.membership, provenance)
    portable = artifact.to_portable()
    portable["activity"][0][0] = 0.99
    with pytest.raises(ValueError, match="fingerprint"):
        compute_scenic_rss(adata, portable, _portable_artifact=True)


def test_portable_scenic_validation_does_not_request_full_activity_copy(
    pyscenic_bundle,
    monkeypatch,
):
    from openbio_singlecell.scenic_artifact import validate_scenic_result_artifact

    adata, manifest_path = pyscenic_bundle
    _output, artifact, _summary = import_pyscenic_bundle(adata, manifest_path)
    portable = artifact.to_portable()
    portable["activity"] = pd.DataFrame(
        portable["activity"],
        index=portable["observations"],
        columns=portable["regulons"],
    )
    copy_requests = []
    original_to_numpy = pd.DataFrame.to_numpy

    def tracked_to_numpy(frame, *args, **kwargs):
        if frame is portable["activity"] and kwargs.get("copy") is True:
            copy_requests.append(frame.shape)
        return original_to_numpy(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_numpy", tracked_to_numpy)
    validate_scenic_result_artifact(
        portable,
        exact_type=False,
        numpy=np,
        pandas=pd,
        copy_result=False,
    )

    assert copy_requests == []

