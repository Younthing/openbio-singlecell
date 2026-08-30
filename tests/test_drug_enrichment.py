from __future__ import annotations

import copy
import hashlib
import json
import os
import socket
import sys
import time
import types
from dataclasses import replace
from pathlib import Path

import pytest
from scipy import stats

from openbio_singlecell import PLUGIN_VERSION
from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.dgidb_resource import (
    DGIDB_RESOURCE_COLUMNS,
    DGIdbResource,
    dgidb_resource_code,
    load_dgidb_resource,
    validate_dgidb_resource,
)
from openbio_singlecell.drug_enrichment import (
    DRUG_GSEA_COLUMNS,
    DRUG_ORA_COLUMNS,
    drug_gsea_code,
    drug_ora_code,
    run_drug_gsea,
    run_drug_ora,
)
from openbio_singlecell.drug_score import drug_score_code, run_drug_score
from openbio_singlecell.enrichment_artifacts import (
    GENERIC_RANKED_ARTIFACT_ROLE,
    GENERIC_SELECTED_ARTIFACT_ROLE,
    GENERIC_UNIVERSE_ARTIFACT_ROLE,
    enrichment_universe_content_fingerprint,
    enrichment_universe_identity_fingerprint,
    generic_enrichment_content_fingerprint,
    generic_ranking_fingerprint,
    validate_enrichment_artifact_pair,
)
from openbio_singlecell.files import resolve_input_path
from openbio_singlecell.nodes_enrichment import (
    OpenBioSingleCellDGIdbAnnotation,
)
from openbio_singlecell.operations_enrichment import (
    drug_gsea_owned,
    drug_hypergeometric_owned,
    drug_scores_owned,
    load_dgidb_resource_owned,
)


def _metadata_json(**updates):
    metadata = {
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
    }
    metadata.update(updates)
    return json.dumps(metadata, ensure_ascii=False, separators=(",", ":"))


def _write_resource(path):
    path.write_text(
        "drug_claim_name\tgene_claim_name\tinteraction_claim_source\tinteraction_types\n"
        'DrugA\tG1\tSource1\t"inhibitor,curated"\n'
        "DrugA\tG2\tSource1\tinhibitor\n"
        "DrugA\tG3\tSource2\tbinder\n"
        "DrugA\tG1\tSource2\tclinical\n"
        "DrugB\tG6\tSource1\tbinder\n"
        "DrugB\tG7\tSource1\tbinder\n"
        "DrugB\tG8\tSource1\tbinder\n"
        "DrugMix\tG2\tSource3\tunknown\n"
        "DrugMix\tG7\tSource3\tunknown\n"
        "DrugZero\tG4\t\t\n"
        "DrugZero\tG5\t\t\n",
        encoding="utf-8-sig",
    )
    return path


def _load_resource(path):
    return load_dgidb_resource(
        str(path),
        requested_path=path.name,
        resource_metadata_json=_metadata_json(),
        max_file_bytes=1_000_000,
        max_rows=1_000,
        openbio_version=PLUGIN_VERSION,
    )


def _generic_artifacts(science, *, purpose, inference_unit="Sample", gene_column="gene"):
    genes = [f"G{index}" for index in range(1, 9)]
    universe_frame = science.pd.DataFrame(
        {"gene": genes, "universe_rank": list(range(1, len(genes) + 1))}
    )
    ranked_frame = science.pd.DataFrame(
        [
            {"comparison": comparison, "gene": gene, "score": score, "rank": rank}
            for comparison, scores in (
                ("treated_vs_control", [8.0, 6.0, 5.0, 2.0, 0.0, -1.0, -3.0, -7.0]),
                ("other_vs_control", [7.0, 5.0, 4.0, 1.0, 0.0, -2.0, -4.0, -8.0]),
            )
            for rank, (gene, score) in enumerate(zip(genes, scores, strict=True), start=1)
        ],
        columns=["comparison", "gene", "score", "rank"],
    )
    selected_frame = science.pd.DataFrame(
        [
            {"comparison": "treated_vs_control", "gene": "G1", "selection_score": 8.0},
            {"comparison": "treated_vs_control", "gene": "G2", "selection_score": 6.0},
            {"comparison": "other_vs_control", "gene": "G1", "selection_score": 7.0},
        ],
        columns=["comparison", "gene", "selection_score"],
    )
    if gene_column != "gene":
        ranked_frame = ranked_frame.rename(columns={"gene": gene_column})
        selected_frame = selected_frame.rename(columns={"gene": gene_column})
    ranking_fingerprint = generic_ranking_fingerprint(
        ranked_frame,
        comparison_column="comparison",
        gene_column=gene_column,
        score_column="score",
    )
    universe_fingerprint = enrichment_universe_identity_fingerprint(genes)
    analysis_fingerprint = "a" * 64
    scope = {
        "evidence_scope": "condition_contrast",
        "inference_unit": inference_unit,
        "direction": "treated versus control; positive values favor treated",
        "replicate_aware": True,
        "technical_batch_handling": "technical_batch included as a nuisance term",
        "condition_inference_level": "sample_level" if inference_unit == "Sample" else "cell_level",
        "organism": "Homo sapiens",
        "identifier_namespace": "HGNC symbol",
        "ranking_truncated": False,
    }
    universe_parameters = {
        "enrichment_evidence_schema_version": 1,
        "enrichment_approved": True,
        "producer_node_id": "OpenBioSingleCellTestUniverse",
        "producer_operation": "test_enrichment_universe",
        "artifact_role": GENERIC_UNIVERSE_ARTIFACT_ROLE,
        "analysis_fingerprint": analysis_fingerprint,
        "ranking_fingerprint": ranking_fingerprint,
        "universe_fingerprint": universe_fingerprint,
        "content_fingerprint": enrichment_universe_content_fingerprint(universe_frame),
        **scope,
    }
    if purpose == "ranked":
        frame = ranked_frame
        role = GENERIC_RANKED_ARTIFACT_ROLE
        operation = "test_ranked_evidence"
        extra = {"gene_column": gene_column, "comparison_column": "comparison", "score_column": "score"}
    else:
        frame = selected_frame
        role = GENERIC_SELECTED_ARTIFACT_ROLE
        operation = "test_selected_evidence"
        extra = {
            "gene_column": gene_column,
            "comparison_column": "comparison",
            "selection_applied": True,
            "upstream_content_fingerprint": generic_enrichment_content_fingerprint(
                ranked_frame, artifact_role=GENERIC_RANKED_ARTIFACT_ROLE
            ),
            "upstream_producer_node_id": "OpenBioSingleCellTestRankedEvidence",
        }
    table_parameters = {
        "enrichment_evidence_schema_version": 1,
        "enrichment_approved": True,
        "producer_node_id": (
            "OpenBioSingleCellTestRankedEvidence"
            if purpose == "ranked"
            else "OpenBioSingleCellTestSelectedEvidence"
        ),
        "producer_operation": operation,
        "artifact_role": role,
        "analysis_fingerprint": analysis_fingerprint,
        "ranking_fingerprint": ranking_fingerprint,
        "universe_fingerprint": universe_fingerprint,
        "content_fingerprint": generic_enrichment_content_fingerprint(frame, artifact_role=role),
        **scope,
        **extra,
    }
    started_at = time.perf_counter()
    table = make_table_result(
        table=frame,
        title="evidence",
        operation=operation,
        parameters=table_parameters,
        description="fixture",
        warnings=[],
        input_cells=48,
        input_genes=len(genes),
        started_at=started_at,
    )
    universe = make_table_result(
        table=universe_frame,
        title="universe",
        operation="test_enrichment_universe",
        parameters=universe_parameters,
        description="fixture",
        warnings=[],
        input_cells=48,
        input_genes=len(genes),
        started_at=started_at,
    )
    return table, universe


def _provenance(science, table, universe, *, purpose, gene_column="gene"):
    return validate_enrichment_artifact_pair(
        table,
        universe,
        purpose=purpose,
        selector="treated_vs_control",
        gene_column=gene_column,
        score_column="score" if purpose == "ranked" else None,
        np=science.np,
        pd=science.pd,
    )


def _with_custom_producer_and_missing_scientific_declarations(result, *, producer_node_id):
    parameters = copy.deepcopy(result.parameters)
    for field in (
        "enrichment_approved",
        "producer_operation",
        "evidence_scope",
        "inference_unit",
        "direction",
        "replicate_aware",
        "technical_batch_handling",
        "condition_inference_level",
        "organism",
        "identifier_namespace",
        "ranking_truncated",
        "selection_applied",
        "upstream_content_fingerprint",
        "upstream_producer_node_id",
    ):
        parameters.pop(field, None)
    if producer_node_id is None:
        parameters.pop("producer_node_id", None)
    else:
        parameters["producer_node_id"] = producer_node_id
    source = copy.deepcopy(result.source)
    source["parameters"] = copy.deepcopy(parameters)
    return replace(result, parameters=parameters, source=source)


def _enrichment_parameters(provenance, *, purpose):
    parameters = {
        "artifact_family": provenance["artifact_family"],
        "comparison": provenance["comparison"],
        "comparison_column": provenance["comparison_column"],
        "gene_column": provenance["gene_column"],
        "evidence_scope": provenance["evidence_scope"],
        "inference_unit": provenance["inference_unit"],
        "direction": provenance["direction"],
        "replicate_aware": provenance["replicate_aware"],
        "technical_batch_handling": provenance["technical_batch_handling"],
        "upstream_parameters": provenance["upstream_parameters"],
        "provenance_warnings": provenance["provenance_warnings"],
        "provenance_declarations": provenance["provenance_declarations"],
        "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
        "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
        "expected_universe_fingerprint": provenance["universe_fingerprint"],
        "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
        "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
        "min_targets": 2,
        "max_output_rows": 100,
        "openbio_version": PLUGIN_VERSION,
    }
    if purpose == "ranked":
        parameters.update(
            {"score_column": provenance["score_column"], "max_targets": 4, "n_permutations": 1000, "random_seed": 17}
        )
    else:
        parameters.update({"min_overlap": 1, "max_p_adjusted": 1.0})
    return parameters


def _fake_decoupler(science):
    module = types.ModuleType("decoupler")
    module.__version__ = "2.2.0"
    module.gsea_calls = []
    module.ora_calls = []

    def gsea(*, data, net, tmin, raw, empty, verbose, times, seed):
        module.gsea_calls.append(
            {"tmin": tmin, "raw": raw, "empty": empty, "verbose": verbose, "times": times, "seed": seed}
        )
        sources = list(dict.fromkeys(net["source"].tolist()))
        values = data.iloc[0].to_dict()
        nes = [
            sum(float(values[gene]) for gene in net.loc[net["source"] == source, "target"]) /
            int((net["source"] == source).sum())
            for source in sources
        ]
        adjusted = [0.01 + 0.1 * index for index in range(len(sources))]
        return (
            science.pd.DataFrame([nes], index=data.index, columns=sources),
            science.pd.DataFrame([adjusted], index=data.index, columns=sources),
        )

    def query_set(*, features, net, alternative, n_bg, ha_corr, tmin, verbose):
        module.ora_calls.append(
            {"alternative": alternative, "n_bg": n_bg, "ha_corr": ha_corr, "tmin": tmin, "verbose": verbose}
        )
        selected = set(features)
        rows = []
        p_values = []
        for source in dict.fromkeys(net["source"].tolist()):
            targets = set(net.loc[net["source"] == source, "target"])
            a = len(selected & targets)
            b = len(targets - selected)
            c = len(selected - targets)
            d = n_bg - a - b - c
            log_odds = float(science.np.log(((a + ha_corr) * (d + ha_corr)) / ((b + ha_corr) * (c + ha_corr))))
            p_value = float(stats.fisher_exact([[a, b], [c, d]], alternative=alternative).pvalue)
            rows.append({"source": source, "stat": log_odds, "pval": p_value})
            p_values.append(p_value)
        adjusted = stats.false_discovery_control(science.np.asarray(p_values), method="bh")
        for row, value in zip(rows, adjusted, strict=True):
            row["padj"] = float(value)
        return science.pd.DataFrame(rows)

    module.mt = types.SimpleNamespace(gsea=gsea, query_set=query_set)
    return module


def _fake_pertpy(science):
    module = types.ModuleType("pertpy")
    module.__version__ = "1.3.0"
    module.calls = []

    class Enrichment:
        def score(self, adata, *, layer, targets, nested, method, key_added):
            module.calls.append(
                {"layer": layer, "targets": copy.deepcopy(targets), "nested": nested, "method": method, "key_added": key_added}
            )
            matrix = adata.layers[layer] if layer is not None else adata.X
            drug, requested = next(iter(targets.items()))
            matched = [gene for gene in adata.var_names if gene in set(requested)]
            indexes = [adata.var_names.get_loc(gene) for gene in matched]
            scores = (
                science.np.asarray(matrix[:, indexes].mean(axis=1)).ravel()
                if science.sparse.issparse(matrix)
                else science.np.asarray(matrix)[:, indexes].mean(axis=1)
            )
            adata.uns[f"{key_added}_score"] = scores[:, None]
            adata.uns[f"{key_added}_variables"] = science.pd.Index([drug])
            genes = science.pd.DataFrame(columns=["genes"]).astype(object)
            all_genes = science.pd.DataFrame(columns=["all_genes"]).astype(object)
            genes.loc[drug, "genes"] = "|".join(matched)
            all_genes.loc[drug, "all_genes"] = "|".join(requested)
            adata.uns[f"{key_added}_genes"] = {"var": genes}
            adata.uns[f"{key_added}_all_genes"] = {"var": all_genes}

    module.tl = types.SimpleNamespace(Enrichment=Enrichment)
    return module


def _adata(science, *, sparse=False):
    values = science.np.asarray(
        [
            [0.1, 0.4, 0.7, 0.2, 0.3, 0.8, 0.9, 0.6],
            [0.3, 0.6, 0.9, 0.4, 0.5, 1.0, 1.1, 0.8],
            [0.5, 0.8, 1.1, 0.6, 0.7, 1.2, 1.3, 1.0],
        ],
        dtype=float,
    )
    matrix = science.sparse.csr_matrix(values) if sparse else values
    adata = science.ad.AnnData(matrix)
    adata.obs_names = ["Cell1", "Cell2", "Cell3"]
    adata.var_names = [f"G{index}" for index in range(1, 9)]
    adata.uns["log1p"] = {"base": None}
    return adata


def test_dgidb_loader_is_strict_tamper_evident_and_generated_equivalent(tmp_path):
    path = _write_resource(tmp_path / "dgidb.tsv")
    resource, summary = _load_resource(path)
    table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    assert list(table.columns) == list(DGIDB_RESOURCE_COLUMNS)
    assert table[["drug", "gene"]].values.tolist()[:4] == [
        ["DrugA", "G1"],
        ["DrugA", "G2"],
        ["DrugA", "G3"],
        ["DrugB", "G6"],
    ]
    first = table.iloc[0]
    assert json.loads(first["sources"]) == ["Source1", "Source2"]
    assert json.loads(first["evidence"]) == ["inhibitor,curated", "clinical"]
    assert accounting["raw_pairs"] == 11
    assert accounting["canonical_pairs"] == 10
    assert accounting["duplicate_pairs_removed"] == 1
    assert artifact["raw_file_sha256"] == summary["key_results"]["resource"]["raw_file_sha256"]
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    source = dgidb_resource_code(
        requested_path=path.name,
        resource_metadata_json=_metadata_json(),
        drug_column="drug_claim_name",
        gene_column="gene_claim_name",
        source_column="interaction_claim_source",
        evidence_column="interaction_types",
        max_file_bytes=1_000_000,
        max_rows=1_000,
        expected_sha256=artifact["raw_file_sha256"],
        openbio_version=PLUGIN_VERSION,
    )
    namespace = {}
    exec(compile(source, "<dgidb-resource-code>", "exec"), namespace)
    generated_table, generated_summary = namespace["load_dgidb_resource"](str(path))
    assert table.equals(generated_table)
    assert summary == generated_summary

    returned = resource.table()
    returned.loc[0, "gene"] = "TAMPERED_COPY"
    assert resource.table().loc[0, "gene"] == "G1"
    assert metadata["gene_identifier_namespace"] == "HGNC symbol"


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"organism": "Mus musculus"}, "Homo sapiens"),
        ({"gene_identifier_namespace": "Ensembl"}, "HGNC symbol"),
        ({"release_date": "2026/08/28"}, "ISO 8601"),
        ({"download_url": "latest"}, "absolute HTTP"),
    ],
)
def test_dgidb_loader_rejects_invalid_scientific_metadata(tmp_path, updates, message):
    path = _write_resource(tmp_path / "dgidb.tsv")
    with pytest.raises(ValueError, match=message):
        load_dgidb_resource(
            str(path),
            requested_path=path.name,
            resource_metadata_json=_metadata_json(**updates),
        )


@pytest.mark.parametrize(
    ("updates", "warning_fragment"),
    [
        ({"source_license_review": "pending"}, "does not declare a completed review"),
        ({"license": "unknown"}, "placeholder or explicitly unreviewed"),
    ],
)
def test_dgidb_loader_allows_caller_declared_license_review_with_warning_and_generated_parity(
    tmp_path, updates, warning_fragment
):
    path = _write_resource(tmp_path / "dgidb.tsv")
    metadata_json = _metadata_json(**updates)
    expected_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
    resource, summary = load_dgidb_resource(
        str(path),
        requested_path=path.name,
        resource_metadata_json=metadata_json,
        max_file_bytes=256 * 1024 * 1024,
        max_rows=5_000_000,
        expected_sha256=expected_sha256,
        openbio_version=PLUGIN_VERSION,
    )
    table, _metadata, _accounting, artifact = validate_dgidb_resource(resource)
    assert any(warning_fragment in warning for warning in summary["warnings"])
    assert any("caller attestations" in warning for warning in summary["warnings"])
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    source = dgidb_resource_code(
        requested_path=path.name,
        resource_metadata_json=metadata_json,
        drug_column="drug_claim_name",
        gene_column="gene_claim_name",
        source_column="interaction_claim_source",
        evidence_column="interaction_types",
        max_file_bytes=256 * 1024 * 1024,
        max_rows=5_000_000,
        expected_sha256=expected_sha256,
        openbio_version=PLUGIN_VERSION,
    )
    namespace = {}
    exec(compile(source, "<dgidb-license-attestation>", "exec"), namespace)
    generated_table, generated_summary = namespace["load_dgidb_resource"](str(path))
    assert generated_table.equals(table)
    assert generated_summary == summary


def test_dgidb_loader_adopts_newly_parsed_table_without_defensive_copy(
    tmp_path,
    monkeypatch,
):
    import openbio_singlecell.dgidb_resource as dgidb

    path = _write_resource(tmp_path / "dgidb.tsv")
    copy_flags = []
    original_validate = dgidb._standalone_validate_dgidb_table

    def tracked_validate(table, copy_table=True):
        copy_flags.append(copy_table)
        return original_validate(table, copy_table=copy_table)

    monkeypatch.setattr(dgidb, "_standalone_validate_dgidb_table", tracked_validate)
    _load_resource(path)

    assert copy_flags
    assert not any(copy_flags)


def test_dgidb_loader_rejects_width_collision_guards_and_pinned_hash(tmp_path):
    malformed = tmp_path / "malformed.tsv"
    malformed.write_text(
        "drug_claim_name\tgene_claim_name\tinteraction_claim_source\tinteraction_types\nDrugA\tG1\tSource1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="fields; expected exactly"):
        _load_resource(malformed)
    collision = tmp_path / "collision.tsv"
    collision.write_text(
        "drug_claim_name\tgene_claim_name\tinteraction_claim_source\tinteraction_types\n"
        "DrugA\tG1\tS\tE\n DrugA \tG2\tS\tE\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="normalization collision"):
        _load_resource(collision)
    path = _write_resource(tmp_path / "dgidb.tsv")
    with pytest.raises(ValueError, match="max_file_bytes"):
        load_dgidb_resource(
            str(path), requested_path=path.name, resource_metadata_json=_metadata_json(), max_file_bytes=1
        )
    with pytest.raises(ValueError, match="max_rows"):
        load_dgidb_resource(
            str(path), requested_path=path.name, resource_metadata_json=_metadata_json(), max_rows=1
        )
    with pytest.raises(ValueError, match="pinned expected bytes"):
        load_dgidb_resource(
            str(path),
            requested_path=path.name,
            resource_metadata_json=_metadata_json(),
            expected_sha256="0" * 64,
        )


def test_dgidb_artifact_detects_private_table_metadata_and_accounting_tamper(tmp_path):
    path = _write_resource(tmp_path / "dgidb.tsv")
    resource, _ = _load_resource(path)
    object.__getattribute__(resource, "_table").loc[0, "gene"] = "OTHER"
    with pytest.raises(ValueError, match="tampered|duplicate|accounting"):
        validate_dgidb_resource(resource)
    resource, _ = _load_resource(path)
    object.__getattribute__(resource, "_metadata")["version"] = "changed"
    with pytest.raises(ValueError, match="tampered"):
        validate_dgidb_resource(resource)
    resource, _ = _load_resource(path)
    object.__getattribute__(resource, "_accounting")["drug_count"] = 99
    with pytest.raises(ValueError, match="drug-count|tampered"):
        validate_dgidb_resource(resource)
    resource, _ = _load_resource(path)
    object.__getattribute__(resource, "_accounting")["source_value_count"] += 1
    with pytest.raises(ValueError, match="source-value accounting"):
        validate_dgidb_resource(resource)
    resource, _ = _load_resource(path)
    object.__getattribute__(resource, "_accounting")["rows_without_evidence"] = 99
    with pytest.raises(ValueError, match="rows_without_evidence"):
        validate_dgidb_resource(resource)
    with pytest.raises(AttributeError, match="immutable"):
        resource.new_field = "x"


def test_dgidb_loader_performs_no_network_access(tmp_path, monkeypatch):
    path = _write_resource(tmp_path / "dgidb.tsv")

    def forbidden(*args, **kwargs):
        raise AssertionError("network access is forbidden")

    monkeypatch.setattr(socket, "create_connection", forbidden)
    resource, _ = _load_resource(path)
    assert isinstance(resource, DGIdbResource)


@pytest.mark.parametrize("as_sparse", [False, True])
def test_drug_score_exact_pertpy_call_mean_and_generated_parity(science, tmp_path, as_sparse):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    adata = _adata(science, sparse=as_sparse)
    original = adata.copy()
    backend = _fake_pertpy(science)
    parameters = {
        "drug": "DrugA",
        "source_kind": "X",
        "layer_name": None,
        "output_key": "drug_target_score",
        "min_matched_targets": 2,
        "overwrite_existing": False,
        "openbio_version": PLUGIN_VERSION,
    }
    output, summary = run_drug_score(adata, resource, **parameters, pertpy_module=backend)
    expected = science.np.asarray(original.X[:, :3].mean(axis=1)).ravel()
    science.np.testing.assert_allclose(output.obs["drug_target_score"], expected)
    assert "_openbio_dgidb_single_drug_score" not in output.uns
    assert backend.calls[0]["targets"] == {"DrugA": ["G1", "G2", "G3"]}
    assert backend.calls[0]["nested"] is False
    assert backend.calls[0]["method"] == "mean"
    science.np.testing.assert_allclose(adata.X.toarray() if as_sparse else adata.X, original.X.toarray() if as_sparse else original.X)
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    source = drug_score_code(
        resource_metadata=metadata,
        resource_accounting=accounting,
        resource_artifact_metadata=artifact,
        **parameters,
    )
    namespace = {}
    exec(compile(source, "<drug-score-code>", "exec"), namespace)
    generated, generated_summary = namespace["run_drug_score"](original.copy(), table.copy(), backend)
    science.np.testing.assert_allclose(generated.obs["drug_target_score"], output.obs["drug_target_score"])
    assert generated_summary == summary


def test_drug_score_guards_source_collision_target_resource_and_backend(science, tmp_path):
    path = _write_resource(tmp_path / "dgidb.tsv")
    resource, _ = _load_resource(path)
    adata = _adata(science)
    with pytest.raises(ValueError, match="one of 'X', 'layer', or 'raw'"):
        run_drug_score(adata, resource, drug="DrugA", source_kind="obsm", output_key="score")
    with pytest.raises(ValueError, match="absent"):
        run_drug_score(adata, resource, drug="Missing", source_kind="X", output_key="score")
    with pytest.raises(ValueError, match="below min_matched_targets"):
        run_drug_score(
            adata[:, ["G1"]].copy(),
            resource,
            drug="DrugA",
            source_kind="X",
            output_key="score",
            min_matched_targets=2,
        )
    collision = adata.copy()
    collision.obs["score"] = 0.0
    with pytest.raises(ValueError, match="already exists"):
        run_drug_score(collision, resource, drug="DrugA", source_kind="X", output_key="score")
    drifted = _fake_pertpy(science)

    def opaque_score(self, adata, **kwargs):
        del self, adata, kwargs

    drifted.tl.Enrichment.score = opaque_score
    with pytest.raises(RuntimeError, match="missing required parameters"):
        run_drug_score(
            adata,
            resource,
            drug="DrugA",
            source_kind="X",
            output_key="score",
            pertpy_module=drifted,
        )
    wrong = _fake_pertpy(science)
    valid = wrong.tl.Enrichment.score

    def altered(self, adata, *, layer, targets, nested, method, key_added):
        valid(
            self,
            adata,
            layer=layer,
            targets=targets,
            nested=nested,
            method=method,
            key_added=key_added,
        )
        adata.uns[f"{key_added}_score"] += 1.0

    wrong.tl.Enrichment.score = altered
    with pytest.raises(RuntimeError, match="independently computed target mean"):
        run_drug_score(
            adata,
            resource,
            drug="DrugA",
            source_kind="X",
            output_key="score",
            pertpy_module=wrong,
        )


@pytest.mark.parametrize("source_kind", ["X", "raw"])
def test_drug_score_allows_count_like_and_all_zero_cells_with_cautious_generated_parity(
    science, tmp_path, source_kind
):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    adata = _adata(science)
    counts = science.np.arange(24, dtype=float).reshape(3, 8)
    counts[0, :] = 0.0
    if source_kind == "raw":
        raw = science.ad.AnnData(
            counts,
            obs=science.pd.DataFrame(index=adata.obs_names.copy()),
            var=science.pd.DataFrame(index=adata.var_names.copy()),
        )
        adata.raw = raw
    else:
        adata.X = counts.copy()
        adata.uns.pop("log1p")
    generated_input = adata.copy()
    backend = _fake_pertpy(science)
    parameters = {
        "drug": "DrugA",
        "source_kind": source_kind,
        "layer_name": None,
        "output_key": "score",
        "min_matched_targets": 2,
        "overwrite_existing": False,
        "openbio_version": PLUGIN_VERSION,
    }
    output, summary = run_drug_score(adata, resource, **parameters, pertpy_module=backend)
    expected = counts[:, :3].mean(axis=1)
    science.np.testing.assert_allclose(output.obs["score"], expected)
    assert summary["key_results"]["all_zero_observation_count"] == 1
    assert any("count-like" in warning for warning in summary["warnings"])
    assert any("all zero" in warning for warning in summary["warnings"])

    table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    source = drug_score_code(
        resource_metadata=metadata,
        resource_accounting=accounting,
        resource_artifact_metadata=artifact,
        **parameters,
    )
    namespace = {}
    exec(compile(source, f"<drug-score-{source_kind}-counts>", "exec"), namespace)
    generated, generated_summary = namespace["run_drug_score"](generated_input, table.copy(), backend)
    science.np.testing.assert_allclose(generated.obs["score"], output.obs["score"])
    assert generated_summary == summary


def test_real_pertpy_13_single_drug_smoke(science, tmp_path):
    pytest.importorskip("pertpy")
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    output, summary = run_drug_score(
        _adata(science),
        resource,
        drug="DrugA",
        source_kind="X",
        output_key="score",
        min_matched_targets=3,
        openbio_version=PLUGIN_VERSION,
    )
    assert output.obs["score"].notna().all()
    assert summary["software_versions"]["pertpy"].startswith("1.3.")


def test_drug_ora_complete_family_contingencies_and_generated_parity(science, tmp_path):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    table, universe = _generic_artifacts(science, purpose="selected")
    parameters = _enrichment_parameters(_provenance(science, table, universe, purpose="selected"), purpose="selected")
    backend = _fake_decoupler(science)
    evidence, summary = run_drug_ora(table, universe, resource, **parameters, decoupler_module=backend)
    assert list(evidence.columns) == DRUG_ORA_COLUMNS
    assert evidence["drug"].tolist() == ["DrugA", "DrugMix", "DrugZero", "DrugB"]
    assert len(evidence) == 4
    assert evidence.loc[evidence["drug"] == "DrugZero", "overlap_count"].item() == 0
    for row in evidence.itertuples(index=False):
        assert row.a == row.overlap_count
        assert row.a + row.b == row.set_size_in_universe
        assert row.a + row.c == row.selected_count
        assert row.a + row.b + row.c + row.d == row.universe_count
    assert backend.ora_calls == [{"alternative": "greater", "n_bg": 8, "ha_corr": 0.5, "tmin": 2, "verbose": False}]
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    resource_table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    source = drug_ora_code(
        resource_metadata=metadata,
        resource_accounting=accounting,
        resource_artifact_metadata=artifact,
        **parameters,
    )
    namespace = {}
    exec(compile(source, "<drug-ora-code>", "exec"), namespace)
    generated, generated_summary = namespace["run_drug_hypergeometric"](
        table.table.copy(), universe.table.copy(), resource_table.copy(), backend
    )
    science.pd.testing.assert_frame_equal(evidence, generated)
    assert summary == generated_summary


def test_drug_gsea_complete_family_adjusted_p_targets_and_generated_parity(science, tmp_path):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    table, universe = _generic_artifacts(science, purpose="ranked")
    parameters = _enrichment_parameters(_provenance(science, table, universe, purpose="ranked"), purpose="ranked")
    backend = _fake_decoupler(science)
    evidence, summary = run_drug_gsea(table, universe, resource, **parameters, decoupler_module=backend)
    assert list(evidence.columns) == DRUG_GSEA_COLUMNS
    assert evidence["drug"].tolist() == ["DrugA", "DrugB", "DrugMix", "DrugZero"]
    assert "p_value" not in evidence
    assert all(json.loads(value) for value in evidence["matched_targets"])
    assert backend.gsea_calls == [
        {"tmin": 2, "raw": False, "empty": False, "verbose": False, "times": 1000, "seed": 17}
    ]
    assert any("not a connectivity-map" in limitation for limitation in summary["limitations"])
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    resource_table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    source = drug_gsea_code(
        resource_metadata=metadata,
        resource_accounting=accounting,
        resource_artifact_metadata=artifact,
        **parameters,
    )
    namespace = {}
    exec(compile(source, "<drug-gsea-code>", "exec"), namespace)
    generated, generated_summary = namespace["run_drug_gsea"](
        table.table.copy(), universe.table.copy(), resource_table.copy(), backend
    )
    science.pd.testing.assert_frame_equal(evidence, generated)
    assert summary == generated_summary


def test_drug_gsea_allows_low_permutations_and_seed_zero_with_cautious_generated_parity(science, tmp_path):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    table, universe = _generic_artifacts(science, purpose="ranked")
    parameters = _enrichment_parameters(_provenance(science, table, universe, purpose="ranked"), purpose="ranked")
    parameters.update({"n_permutations": 10, "random_seed": 0})
    evidence, summary = run_drug_gsea(
        table,
        universe,
        resource,
        **parameters,
        decoupler_module=_fake_decoupler(science),
    )
    assert summary["key_results"]["permutation_null_status"] == "not_shuffled_seed_zero"
    assert summary["key_results"]["formal_permutation_inference_valid"] is False
    assert any("Only 10 permutations" in warning for warning in summary["warnings"])
    assert any("disables permutation shuffling" in warning for warning in summary["warnings"])
    assert "computational flags" in summary["results"]

    resource_table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    source = drug_gsea_code(
        resource_metadata=metadata,
        resource_accounting=accounting,
        resource_artifact_metadata=artifact,
        **parameters,
    )
    namespace = {}
    exec(compile(source, "<drug-gsea-low-permutations>", "exec"), namespace)
    generated, generated_summary = namespace["run_drug_gsea"](
        table.table.copy(), universe.table.copy(), resource_table.copy(), _fake_decoupler(science)
    )
    science.pd.testing.assert_frame_equal(evidence, generated)
    assert generated_summary == summary


@pytest.mark.parametrize(
    ("purpose", "runner", "code_builder", "generated_name"),
    [
        ("selected", run_drug_ora, drug_ora_code, "run_drug_hypergeometric"),
        ("ranked", run_drug_gsea, drug_gsea_code, "run_drug_gsea"),
    ],
)
def test_drug_enrichment_keeps_custom_evidence_gene_column_separate_from_canonical_universe(
    science, tmp_path, purpose, runner, code_builder, generated_name
):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    table, universe = _generic_artifacts(science, purpose=purpose, gene_column="symbol")
    provenance = _provenance(science, table, universe, purpose=purpose, gene_column="symbol")
    parameters = _enrichment_parameters(provenance, purpose=purpose)
    backend = _fake_decoupler(science)

    evidence, summary = runner(table, universe, resource, **parameters, decoupler_module=backend)
    resource_table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    source = code_builder(
        resource_metadata=metadata,
        resource_accounting=accounting,
        resource_artifact_metadata=artifact,
        **parameters,
    )
    namespace = {}
    exec(compile(source, "<drug-custom-gene-column-code>", "exec"), namespace)
    generated, generated_summary = namespace[generated_name](
        table.table.copy(), universe.table.copy(), resource_table.copy(), backend
    )

    assert parameters["gene_column"] == "symbol"
    assert universe.table.columns.tolist() == ["gene", "universe_rank"]
    science.pd.testing.assert_frame_equal(evidence, generated)
    assert summary == generated_summary


@pytest.mark.parametrize(
    ("purpose", "runner", "code_builder", "generated_name"),
    [
        ("selected", run_drug_ora, drug_ora_code, "run_drug_hypergeometric"),
        ("ranked", run_drug_gsea, drug_gsea_code, "run_drug_gsea"),
    ],
)
def test_drug_enrichment_custom_producer_and_missing_scientific_metadata_are_cautious(
    science, tmp_path, purpose, runner, code_builder, generated_name
):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    table, universe = _generic_artifacts(science, purpose=purpose)
    table = _with_custom_producer_and_missing_scientific_declarations(
        table,
        producer_node_id="ExternalExpertEvidence",
    )
    universe = _with_custom_producer_and_missing_scientific_declarations(
        universe,
        producer_node_id=None,
    )
    provenance = _provenance(science, table, universe, purpose=purpose)
    parameters = _enrichment_parameters(provenance, purpose=purpose)
    backend = _fake_decoupler(science)
    evidence, summary = runner(
        table,
        universe,
        resource,
        **parameters,
        decoupler_module=backend,
    )
    assert summary["key_results"]["provenance"]["scientific_metadata_status"] == (
        "caller_declared_not_programmatically_verified"
    )
    assert any("caller declarations" in warning for warning in summary["warnings"])
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    resource_table, metadata, accounting, artifact = validate_dgidb_resource(resource)
    source = code_builder(
        resource_metadata=metadata,
        resource_accounting=accounting,
        resource_artifact_metadata=artifact,
        **parameters,
    )
    namespace = {}
    exec(compile(source, f"<drug-{purpose}-custom-producer-code>", "exec"), namespace)
    generated, generated_summary = namespace[generated_name](
        table.table.copy(),
        universe.table.copy(),
        resource_table.copy(),
        _fake_decoupler(science),
    )
    science.pd.testing.assert_frame_equal(evidence, generated)
    assert summary == generated_summary


def test_real_decoupler_22_drug_ora_and_gsea_smoke(science, tmp_path):
    decoupler = pytest.importorskip("decoupler")
    if not getattr(decoupler, "__version__", "").startswith("2.2."):
        pytest.skip("audited decoupler 2.2 is unavailable")
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    selected, selected_universe = _generic_artifacts(science, purpose="selected")
    ora_parameters = _enrichment_parameters(
        _provenance(science, selected, selected_universe, purpose="selected"), purpose="selected"
    )
    ora, _ = run_drug_ora(selected, selected_universe, resource, **ora_parameters)
    assert ora["p_adjusted"].between(0.0, 1.0).all()
    ranked, ranked_universe = _generic_artifacts(science, purpose="ranked")
    gsea_parameters = _enrichment_parameters(
        _provenance(science, ranked, ranked_universe, purpose="ranked"), purpose="ranked"
    )
    gsea, _ = run_drug_gsea(ranked, ranked_universe, resource, **gsea_parameters)
    assert gsea["p_adjusted"].between(0.0, 1.0).all()


def test_drug_enrichment_cell_level_is_cautious_but_cross_pair_and_resource_tamper_fail(
    science, tmp_path
):
    resource, _ = _load_resource(_write_resource(tmp_path / "dgidb.tsv"))
    table, universe = _generic_artifacts(science, purpose="selected", inference_unit="Cell")
    provenance = _provenance(science, table, universe, purpose="selected")
    parameters = _enrichment_parameters(provenance, purpose="selected")
    _, cell_level_summary = run_drug_ora(
        table,
        universe,
        resource,
        **parameters,
        decoupler_module=_fake_decoupler(science),
    )
    assert cell_level_summary["key_results"]["provenance"][
        "formal_condition_interpretation_validated"
    ] is False
    assert any(
        "formal Condition interpretation is not validated" in warning
        for warning in cell_level_summary["warnings"]
    )
    ranked, ranked_universe = _generic_artifacts(science, purpose="ranked")
    parameters = _enrichment_parameters(
        _provenance(science, ranked, ranked_universe, purpose="ranked"), purpose="ranked"
    )
    foreign_parameters = copy.deepcopy(ranked_universe.parameters)
    foreign_parameters["ranking_fingerprint"] = "f" * 64
    foreign_source = copy.deepcopy(ranked_universe.source)
    foreign_source["parameters"] = copy.deepcopy(foreign_parameters)
    foreign = replace(ranked_universe, parameters=foreign_parameters, source=foreign_source)
    with pytest.raises(ValueError, match="ranking_fingerprint changed"):
        run_drug_gsea(ranked, foreign, resource, **parameters, decoupler_module=_fake_decoupler(science))
    object.__getattribute__(resource, "_table").loc[0, "gene"] = "OTHER"
    with pytest.raises(ValueError, match="tampered|accounting"):
        run_drug_gsea(
            ranked,
            ranked_universe,
            resource,
            **parameters,
            decoupler_module=_fake_decoupler(science),
        )


def test_drug_nodes_expose_atomic_typed_contracts_and_hash_resource(
    science, comfy_directories, monkeypatch
):
    input_dir, _, _ = comfy_directories
    path = _write_resource(input_dir / "dgidb.tsv")
    resource_output = load_dgidb_resource_owned(
        str(resolve_input_path(path.name, extensions=(".csv", ".tsv"))),
        requested_path=path.name,
        resource_metadata_json=_metadata_json(),
        max_file_bytes=1_000_000,
        max_rows=1_000,
    )
    resource = resource_output[0]
    assert isinstance(resource, DGIdbResource)
    assert resource_output[1].kind == "summary"
    compile(resource_output[2], "<resource-node-code>", "exec")
    before = OpenBioSingleCellDGIdbAnnotation.fingerprint_inputs(path.name)
    with pytest.raises(ValueError, match="max_file_bytes"):
        OpenBioSingleCellDGIdbAnnotation.fingerprint_inputs(path.name, max_file_bytes=1)
    raw = path.read_bytes()
    stat = path.stat()
    replacement = raw.replace(b"DrugA", b"DrugZ", 1)
    assert len(replacement) == len(raw)
    path.write_bytes(replacement)
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    after = OpenBioSingleCellDGIdbAnnotation.fingerprint_inputs(path.name)
    assert before != after
    _, _, _, original_artifact = validate_dgidb_resource(resource)
    changed_resource, _ = _load_resource(path)
    _, _, _, changed_artifact = validate_dgidb_resource(changed_resource)
    assert changed_artifact["raw_file_sha256"] != original_artifact["raw_file_sha256"]
    assert (
        changed_artifact["canonical_content_fingerprint_sha256"]
        != original_artifact["canonical_content_fingerprint_sha256"]
    )

    score_backend = _fake_pertpy(science)
    monkeypatch.setitem(sys.modules, "pertpy", score_backend)
    score_output = drug_scores_owned(
        _adata(science),
        resource,
        drug="DrugA",
        source={"source": "X"},
        output_key="score",
    )
    assert len(score_output) == 3 and score_output[1].kind == "summary"
    compile(score_output[2], "<score-node-code>", "exec")

    decoupler = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", decoupler)
    selected, selected_universe = _generic_artifacts(science, purpose="selected")
    ora_output = drug_hypergeometric_owned(
        selected,
        selected_universe,
        resource,
        comparison="treated_vs_control",
        min_targets=2,
        min_overlap=1,
        max_p_adjusted=1.0,
    )
    assert len(ora_output) == 3 and ora_output[0].kind == "table" and ora_output[1].kind == "summary"
    compile(ora_output[2], "<ora-node-code>", "exec")
    ranked, ranked_universe = _generic_artifacts(science, purpose="ranked")
    gsea_output = drug_gsea_owned(
        ranked,
        ranked_universe,
        resource,
        comparison="treated_vs_control",
        min_targets=2,
        max_targets=4,
        n_permutations=1000,
        random_seed=17,
    )
    assert len(gsea_output) == 3 and gsea_output[0].kind == "table" and gsea_output[1].kind == "summary"
    compile(gsea_output[2], "<gsea-node-code>", "exec")


def test_drug_nodes_have_no_legacy_remote_or_internal_ranking_paths():
    source = Path("openbio_singlecell/nodes_enrichment.py").read_text(encoding="utf-8")
    drug_source = source[source.index("class OpenBioSingleCellDGIdbAnnotation") : source.index("ENRICHMENT_NODE_CLASSES")]
    assert "rank_genes_groups" not in drug_source
    assert ".dgidb.dictionary" not in drug_source
    assert ".annotate(" not in drug_source
    assert "groupby" not in drug_source
