from __future__ import annotations

import copy
import hashlib
import json
import os
import sys
import time
import types
from dataclasses import replace

import pytest
from scipy import stats

from openbio_singlecell import PLUGIN_VERSION
from openbio_singlecell.analysis_utils import make_table_result
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
from openbio_singlecell.enrichment_resource import load_gene_set_resource
from openbio_singlecell.generic_ora import (
    GENERIC_ORA_COLUMNS,
    build_generic_ora_summary,
    generic_ora_code,
    run_generic_ora_evidence,
)
from openbio_singlecell.marker_evidence import (
    MARKER_COLUMNS,
    marker_table_content_fingerprint,
    marker_universe_content_fingerprint,
)
from openbio_singlecell.nodes_enrichment import (
    OpenBioSingleCellGeneSetOverrepresentation,
    OpenBioSingleCellRankedGSEA,
)
from openbio_singlecell.ranked_enrichment import (
    RANKED_GSEA_COLUMNS,
    build_ranked_gsea_summary,
    ranked_gsea_code,
    run_ranked_gsea_evidence,
)


def _output_values(node_output):
    return node_output.result


def _metadata_json() -> str:
    return json.dumps(
        {
            "name": "Test pathways",
            "version": "2026.08",
            "date": "2026-08-28",
            "organism": "Homo sapiens",
            "identifier_namespace": "HGNC symbol",
            "scope": "test resource",
            "license": "CC-BY-4.0",
            "citation": "Deterministic test resource.",
        },
        separators=(",", ":"),
    )


def _write_resource(path):
    path.write_text(
        "source,target\n"
        "UP,G1\nUP,G2\nUP,G3\nUP,G1\n"
        "DOWN,G6\nDOWN,G7\nDOWN,G8\n"
        "MIX,G2\nMIX,G7\n"
        "ZERO,G4\nZERO,G5\n"
        "BIG,G1\nBIG,G2\nBIG,G3\nBIG,G4\nBIG,G5\n"
        "ALL,G1\nALL,G2\nALL,G3\nALL,G4\nALL,G5\nALL,G6\nALL,G7\nALL,G8\n",
        encoding="utf-8",
    )
    return path


def _generic_artifacts(science, *, purpose="ranked", inference_unit="Sample"):
    genes = [f"G{index}" for index in range(1, 9)]
    universe_frame = science.pd.DataFrame(
        {"gene": genes, "universe_rank": list(range(1, len(genes) + 1))}
    )
    ranked_frame = science.pd.DataFrame(
        [
            {"comparison": comparison, "gene": gene, "score": score, "rank": rank}
            for comparison, ordered_scores in (
                ("treated_vs_control", [8.0, 6.0, 6.0, 2.0, 0.0, -1.0, -3.0, -7.0]),
                ("other_vs_control", [7.0, 5.0, 4.0, 1.0, 0.0, -2.0, -4.0, -8.0]),
            )
            for rank, (gene, score) in enumerate(zip(genes, ordered_scores, strict=True), start=1)
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
    ranking_fingerprint = generic_ranking_fingerprint(
        ranked_frame,
        comparison_column="comparison",
        gene_column="gene",
        score_column="score",
    )
    universe_identity = enrichment_universe_identity_fingerprint(genes)
    analysis_fingerprint = "a" * 64
    scope_fields = {
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
        "universe_fingerprint": universe_identity,
        "content_fingerprint": enrichment_universe_content_fingerprint(universe_frame),
        **scope_fields,
    }
    if purpose == "ranked":
        frame = ranked_frame
        role = GENERIC_RANKED_ARTIFACT_ROLE
        operation = "test_ranked_evidence"
        extra = {
            "gene_column": "gene",
            "comparison_column": "comparison",
            "score_column": "score",
        }
    else:
        frame = selected_frame
        role = GENERIC_SELECTED_ARTIFACT_ROLE
        operation = "test_selected_evidence"
        extra = {
            "gene_column": "gene",
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
        "universe_fingerprint": universe_identity,
        "content_fingerprint": generic_enrichment_content_fingerprint(frame, artifact_role=role),
        **scope_fields,
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


def _fake_decoupler(science):
    module = types.ModuleType("decoupler")
    module.__version__ = "2.2.0"
    module.gsea_calls = []
    module.ora_calls = []

    def gsea(*, data, net, tmin, raw, empty, verbose, times, seed):
        module.gsea_calls.append(
            {
                "data": data.copy(),
                "net": net.copy(),
                "tmin": tmin,
                "raw": raw,
                "empty": empty,
                "verbose": verbose,
                "times": times,
                "seed": seed,
            }
        )
        sources = list(dict.fromkeys(net["source"].tolist()))
        values = {gene: float(data.iloc[0][gene]) for gene in data.columns}
        nes = []
        for source in sources:
            targets = net.loc[net["source"] == source, "target"].tolist()
            nes.append(sum(values[target] for target in targets) / len(targets))
        adjusted = [0.01 + 0.1 * index for index in range(len(sources))]
        return (
            science.pd.DataFrame([nes], index=data.index.copy(), columns=sources),
            science.pd.DataFrame([adjusted], index=data.index.copy(), columns=sources),
        )

    def query_set(*, features, net, alternative, n_bg, ha_corr, tmin, verbose):
        module.ora_calls.append(
            {
                "features": list(features),
                "net": net.copy(),
                "alternative": alternative,
                "n_bg": n_bg,
                "ha_corr": ha_corr,
                "tmin": tmin,
                "verbose": verbose,
            }
        )
        selected = set(features)
        rows = []
        p_values = []
        sources = list(dict.fromkeys(net["source"].tolist()))
        source_rows = []
        for source in sources:
            targets = set(net.loc[net["source"] == source, "target"].tolist())
            a = len(selected & targets)
            b = len(targets - selected)
            c = len(selected - targets)
            d = n_bg - a - b - c
            log_odds = float(science.np.log(((a + ha_corr) * (d + ha_corr)) / ((b + ha_corr) * (c + ha_corr))))
            p_value = float(stats.fisher_exact([[a, b], [c, d]], alternative=alternative).pvalue)
            source_rows.append((source, log_odds, p_value))
            p_values.append(p_value)
        adjusted = stats.false_discovery_control(science.np.asarray(p_values), method="bh")
        for (source, log_odds, p_value), padj in zip(source_rows, adjusted, strict=True):
            rows.append({"source": source, "stat": log_odds, "pval": p_value, "padj": float(padj)})
        return science.pd.DataFrame(rows)

    module.mt = types.SimpleNamespace(gsea=gsea, query_set=query_set)
    return module


def _marker_artifacts(science, *, filtered):
    genes = [f"G{index}" for index in range(1, 9)]
    rows = []
    for group, scores in (
        ("A", [8.0, 6.0, 5.0, 2.0, 0.0, -1.0, -3.0, -7.0]),
        ("B", [7.0, 5.0, 4.0, 1.0, 0.0, -2.0, -4.0, -8.0]),
    ):
        for rank, (gene, score) in enumerate(zip(genes, scores, strict=True), start=1):
            rows.append(
                {
                    "group": group,
                    "gene": gene,
                    "rank": rank,
                    "score": score,
                    "log2_fold_change_approx": 2.0 if rank <= 2 else -0.5,
                    "p_value": 0.001 * rank,
                    "p_adjusted": min(1.0, 0.01 * rank),
                    "fraction_in_group": 0.8 if rank <= 2 else 0.1,
                    "fraction_reference": 0.1 if rank <= 2 else 0.8,
                }
            )
    complete = science.pd.DataFrame(rows, columns=MARKER_COLUMNS)
    frame = complete.loc[
        (complete["log2_fold_change_approx"] >= 0.5)
        & (complete["fraction_in_group"] >= 0.2)
        & (complete["fraction_reference"] <= 0.5)
        & (complete["p_adjusted"] <= 0.05)
    ].reset_index(drop=True) if filtered else complete
    universe_frame = science.pd.DataFrame(
        {"gene": genes, "universe_rank": list(range(1, 9))},
        columns=["gene", "universe_rank"],
    )
    identity = enrichment_universe_identity_fingerprint(genes)
    shared = {
        "marker_evidence_schema_version": 2,
        "analysis_fingerprint": "d" * 64,
        "ranking_fingerprint": "e" * 64,
        "universe_fingerprint": identity,
        "marker_groupby": "cluster",
        "marker_method": "wilcoxon",
        "marker_source": "X",
        "group_labels": ["A", "B"],
        "group_sizes": {"A": 24, "B": 24},
        "source_gene_count": len(genes),
        "requested_n_genes": 0,
        "actual_n_genes_per_group": len(genes),
        "ranking_truncated": False,
        "tie_correct": True,
        "max_output_rows": 1000,
        "max_working_memory_gib": 1.0,
        "constant_gene_count": 0,
        "variable_gene_count": len(genes),
        "bh_hypotheses_per_group": len(genes),
        "marker_output_rows": 2 * len(genes),
        "universe_output_rows": len(genes),
        "total_output_rows": 3 * len(genes),
        "groups": "all",
        "reference": "rest",
        "rankby_abs": False,
        "pts": True,
        "corr_method": "benjamini-hochberg",
        "expression_state": "logged",
        "expression_interpretation": "verified logarithmized abundance",
    }
    if filtered:
        operation = "filter_marker_genes"
        parameters = {
            **shared,
            "producer_node_id": "OpenBioSingleCellFilterMarkerGenes",
            "upstream_producer_node_id": "OpenBioSingleCellMarkerGenes",
            "artifact_role": "marker_table",
            "content_fingerprint": marker_table_content_fingerprint(frame),
            "upstream_content_fingerprint": marker_table_content_fingerprint(complete),
            "min_log2_fold_change": 0.5,
            "min_fraction_in_group": 0.2,
            "max_fraction_reference": 0.5,
            "max_p_adjusted": 0.05,
            "comparison_semantics": "inclusive",
        }
    else:
        operation = "marker_genes"
        parameters = {
            **shared,
            "producer_node_id": "OpenBioSingleCellMarkerGenes",
            "artifact_role": "marker_table",
            "content_fingerprint": marker_table_content_fingerprint(frame),
        }
    universe_parameters = {
        **shared,
        "producer_node_id": "OpenBioSingleCellMarkerGenes",
        "artifact_role": "tested_gene_universe",
        "content_fingerprint": marker_universe_content_fingerprint(universe_frame),
    }
    started = time.perf_counter()
    table = make_table_result(
        table=frame,
        title="markers",
        operation=operation,
        parameters=parameters,
        description="fixture",
        warnings=[],
        input_cells=48,
        input_genes=8,
        started_at=started,
    )
    universe = make_table_result(
        table=universe_frame,
        title="universe",
        operation="marker_genes",
        parameters=universe_parameters,
        description="fixture",
        warnings=[],
        input_cells=48,
        input_genes=8,
        started_at=started,
    )
    return table, universe


def _provenance(science, table, universe, *, purpose, selector="treated_vs_control"):
    return validate_enrichment_artifact_pair(
        table,
        universe,
        purpose=purpose,
        selector=selector,
        gene_column="gene",
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


def _base_parameters(provenance, resource_path, *, purpose):
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
        "resource_path": str(resource_path),
        "requested_resource_path": str(resource_path),
        "resource_metadata_json": _metadata_json(),
        "source_column": "source",
        "target_column": "target",
        "min_targets": 2,
        "max_output_rows": 100,
        "openbio_version": PLUGIN_VERSION,
    }
    if purpose == "ranked":
        parameters.update(
            {
                "score_column": provenance["score_column"],
                "max_targets": 4,
                "n_permutations": 1000,
                "random_seed": 17,
            }
        )
    else:
        parameters.update({"min_overlap": 2, "max_p_adjusted": 0.2})
    return parameters


def test_resource_loader_is_strict_complete_and_ordered(science, tmp_path):
    path = _write_resource(tmp_path / "sets.csv")
    resource = load_gene_set_resource(
        str(path),
        _metadata_json(),
        "source",
        "target",
        [f"G{index}" for index in range(1, 9)],
        2,
    )
    assert resource["source_order"] == ("UP", "DOWN", "MIX", "ZERO", "BIG", "ALL")
    assert resource["accounting"]["duplicate_pairs_removed"] == 1
    assert resource["accounting"]["source_count_after_min_targets"] == 6
    assert list(resource["network"].columns) == ["source", "target"]
    json.dumps(resource["metadata"], allow_nan=False)


@pytest.mark.parametrize("extension", ["tsv", "gmt"])
def test_resource_loader_supports_strict_tsv_and_gmt(science, tmp_path, extension):
    path = tmp_path / f"sets.{extension}"
    if extension == "tsv":
        path.write_text("source\ttarget\nA\tG1\nA\tG2\n", encoding="utf-8-sig")
    else:
        path.write_text("A\tdescription\tG1\tG2\n", encoding="utf-8-sig")
    resource = load_gene_set_resource(
        str(path), _metadata_json(), "source", "target", ["G1", "G2", "G3"], 2
    )
    assert resource["retained_sources"] == ("A",)
    assert resource["targets_in_universe"]["A"] == ("G1", "G2")


def test_resource_loader_rejects_malformed_rows_metadata_and_changed_hash(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("source,target\nA,G1,EXTRA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="logical data row"):
        load_gene_set_resource(str(path), _metadata_json(), "source", "target", ["G1"], 1)
    good = _write_resource(tmp_path / "good.csv")
    original_sha = hashlib.sha256(good.read_bytes()).hexdigest()
    stat = good.stat()
    data = good.read_bytes().replace(b"UP,G1", b"ZZ,G1", 1)
    assert len(data) == stat.st_size
    good.write_bytes(data)
    os.utime(good, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    with pytest.raises(ValueError, match="fingerprint changed"):
        load_gene_set_resource(
            str(good),
            _metadata_json(),
            "source",
            "target",
            [f"G{index}" for index in range(1, 9)],
            1,
            expected_sha256=original_sha,
        )
    duplicate_metadata = _metadata_json().replace('"version":"2026.08"', '"version":"2026.08","version":"x"')
    with pytest.raises(ValueError, match="duplicate key"):
        load_gene_set_resource(str(good), duplicate_metadata, "source", "target", ["G1"], 1)


def test_ranked_gsea_runtime_generated_parity_and_complete_family(science, tmp_path):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="ranked")
    provenance = _provenance(science, table, universe, purpose="ranked")
    parameters = _base_parameters(provenance, resource_path, purpose="ranked")
    backend = _fake_decoupler(science)
    before_table = table.table.copy(deep=True)
    before_universe = universe.table.copy(deep=True)
    evidence, diagnostics = run_ranked_gsea_evidence(
        table, universe, **parameters, decoupler_module=backend
    )
    summary = build_ranked_gsea_summary(diagnostics)
    assert list(evidence.columns) == RANKED_GSEA_COLUMNS
    assert set(evidence["gene_set"]) == {"UP", "DOWN", "MIX", "ZERO"}
    assert "p_value" not in evidence and "p_adjusted" in evidence
    assert diagnostics["resource"]["sources_removed_by_max_targets_count"] == 1
    assert diagnostics["resource"]["whole_universe_sources_removed_count"] == 1
    assert backend.gsea_calls[0]["times"] == 1000
    assert backend.gsea_calls[0]["seed"] == 17
    assert backend.gsea_calls[0]["raw"] is False
    assert summary["node_id"] == "OpenBioSingleCellRankedGSEA"
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    science.pd.testing.assert_frame_equal(table.table, before_table)
    science.pd.testing.assert_frame_equal(universe.table, before_universe)

    generated = ranked_gsea_code(
        **parameters,
        expected_resource_sha256=diagnostics["resource"]["sha256"],
    )
    compile(generated, "<ranked-gsea-code>", "exec")
    namespace = {}
    exec(generated, namespace)
    generated_evidence, generated_summary = namespace["run_ranked_gsea"](
        table.table.copy(), universe.table.copy(), _fake_decoupler(science)
    )
    science.pd.testing.assert_frame_equal(evidence, generated_evidence)
    assert summary == generated_summary


def test_generic_ora_runtime_generated_parity_and_contingencies(science, tmp_path):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="selected")
    provenance = _provenance(science, table, universe, purpose="selected")
    parameters = _base_parameters(provenance, resource_path, purpose="selected")
    backend = _fake_decoupler(science)
    evidence, diagnostics = run_generic_ora_evidence(
        table, universe, **parameters, decoupler_module=backend
    )
    summary = build_generic_ora_summary(diagnostics)
    assert list(evidence.columns) == GENERIC_ORA_COLUMNS
    assert set(evidence["gene_set"]) == {"UP", "DOWN", "MIX", "ZERO", "BIG", "ALL"}
    assert bool((evidence["overlap_count"] == 0).any())
    assert bool(((evidence["a"] + evidence["b"] + evidence["c"] + evidence["d"]) == 8).all())
    assert bool((evidence["p_adjusted"] >= evidence["p_value"] - 1e-15).all())
    assert backend.ora_calls[0]["alternative"] == "greater"
    assert backend.ora_calls[0]["n_bg"] == 8
    assert backend.ora_calls[0]["ha_corr"] == 0.5
    assert summary["node_id"] == "OpenBioSingleCellGeneSetOverrepresentation"
    json.dumps(summary, ensure_ascii=False, allow_nan=False)

    generated = generic_ora_code(
        **parameters,
        expected_resource_sha256=diagnostics["resource"]["sha256"],
    )
    compile(generated, "<generic-ora-code>", "exec")
    namespace = {}
    exec(generated, namespace)
    generated_evidence, generated_summary = namespace["run_gene_set_overrepresentation"](
        table.table.copy(), universe.table.copy(), _fake_decoupler(science)
    )
    science.pd.testing.assert_frame_equal(evidence, generated_evidence)
    assert summary == generated_summary


def test_marker_artifact_contracts_are_supported_only_at_correct_stage(science, tmp_path):
    resource_path = _write_resource(tmp_path / "sets.csv")
    ranked, ranked_universe = _marker_artifacts(science, filtered=False)
    ranked_provenance = _provenance(
        science, ranked, ranked_universe, purpose="ranked", selector="A"
    )
    assert ranked_provenance["artifact_family"] == "marker_v2"
    assert ranked_provenance["evidence_scope"] == "cluster_marker_evidence"
    ranked_parameters = _base_parameters(ranked_provenance, resource_path, purpose="ranked")
    ranked_evidence, ranked_diagnostics = run_ranked_gsea_evidence(
        ranked,
        ranked_universe,
        **ranked_parameters,
        decoupler_module=_fake_decoupler(science),
    )
    assert not ranked_evidence.empty
    assert any("Cluster marker evidence" in warning for warning in ranked_diagnostics["warnings"])

    selected, selected_universe = _marker_artifacts(science, filtered=True)
    selected_provenance = _provenance(
        science, selected, selected_universe, purpose="selected", selector="A"
    )
    selected_parameters = _base_parameters(selected_provenance, resource_path, purpose="selected")
    selected_evidence, selected_diagnostics = run_generic_ora_evidence(
        selected,
        selected_universe,
        **selected_parameters,
        decoupler_module=_fake_decoupler(science),
    )
    assert not selected_evidence.empty
    assert selected_diagnostics["selection"]["evidence_scope"] == "cluster_marker_evidence"

    with pytest.raises(ValueError, match="table operation"):
        _provenance(science, selected, selected_universe, purpose="ranked", selector="A")
    with pytest.raises(ValueError, match="table operation"):
        _provenance(science, ranked, ranked_universe, purpose="selected", selector="A")


def test_real_decoupler_22_official_gsea_and_query_set_contracts(science, tmp_path):
    decoupler = pytest.importorskip("decoupler")
    if not str(decoupler.__version__).startswith("2.2."):
        pytest.skip("This contract test targets the audited decoupler 2.2 family.")
    resource_path = _write_resource(tmp_path / "sets.csv")

    ranked, ranked_universe = _generic_artifacts(science, purpose="ranked")
    ranked_provenance = _provenance(science, ranked, ranked_universe, purpose="ranked")
    ranked_parameters = _base_parameters(ranked_provenance, resource_path, purpose="ranked")
    ranked_evidence, ranked_diagnostics = run_ranked_gsea_evidence(
        ranked,
        ranked_universe,
        **ranked_parameters,
        decoupler_module=decoupler,
    )
    assert list(ranked_evidence.columns) == RANKED_GSEA_COLUMNS
    assert bool(ranked_evidence["p_adjusted"].between(0.0, 1.0).all())
    assert ranked_diagnostics["backend"]["api"] == "decoupler.mt.gsea"
    ranked_source = ranked_gsea_code(
        **ranked_parameters,
        expected_resource_sha256=ranked_diagnostics["resource"]["sha256"],
    )
    ranked_namespace = {}
    exec(compile(ranked_source, "<real-ranked-gsea-code>", "exec"), ranked_namespace)
    generated_ranked, generated_ranked_summary = ranked_namespace["run_ranked_gsea"](
        ranked.table.copy(), ranked_universe.table.copy(), decoupler
    )
    science.pd.testing.assert_frame_equal(ranked_evidence, generated_ranked)
    assert build_ranked_gsea_summary(ranked_diagnostics) == generated_ranked_summary

    selected, selected_universe = _generic_artifacts(science, purpose="selected")
    selected_provenance = _provenance(science, selected, selected_universe, purpose="selected")
    selected_parameters = _base_parameters(selected_provenance, resource_path, purpose="selected")
    ora_evidence, ora_diagnostics = run_generic_ora_evidence(
        selected,
        selected_universe,
        **selected_parameters,
        decoupler_module=decoupler,
    )
    assert list(ora_evidence.columns) == GENERIC_ORA_COLUMNS
    assert bool(ora_evidence["p_adjusted"].between(0.0, 1.0).all())
    assert ora_diagnostics["backend"]["api"] == "decoupler.mt.query_set"
    ora_source = generic_ora_code(
        **selected_parameters,
        expected_resource_sha256=ora_diagnostics["resource"]["sha256"],
    )
    ora_namespace = {}
    exec(compile(ora_source, "<real-generic-ora-code>", "exec"), ora_namespace)
    generated_ora, generated_ora_summary = ora_namespace["run_gene_set_overrepresentation"](
        selected.table.copy(), selected_universe.table.copy(), decoupler
    )
    science.pd.testing.assert_frame_equal(ora_evidence, generated_ora)
    assert build_generic_ora_summary(ora_diagnostics) == generated_ora_summary


def test_ranked_gsea_discloses_backend_tie_policy_and_result_ties_share_rank(science, tmp_path):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="ranked")
    provenance = _provenance(science, table, universe, purpose="ranked")
    parameters = _base_parameters(provenance, resource_path, purpose="ranked")
    backend = _fake_decoupler(science)

    def tied_gsea(*, data, net, tmin, raw, empty, verbose, times, seed):
        del tmin, raw, empty, verbose, times, seed
        sources = list(dict.fromkeys(net["source"].tolist()))
        assert sources == ["UP", "DOWN", "MIX", "ZERO"]
        return (
            science.pd.DataFrame([[2.0, -2.0, 1.0, 1.0]], index=data.index, columns=sources),
            science.pd.DataFrame([[0.05, 0.05, 0.1, 0.1]], index=data.index, columns=sources),
        )

    backend.mt.gsea = tied_gsea
    evidence, diagnostics = run_ranked_gsea_evidence(
        table, universe, **parameters, decoupler_module=backend
    )
    assert evidence["gene_set"].tolist() == ["UP", "DOWN", "MIX", "ZERO"]
    assert evidence["rank"].tolist() == [1, 2, 3, 3]
    assert diagnostics["ranking"]["exact_tie_block_count"] == 1
    assert diagnostics["ranking"]["producer_row_order_fingerprinted"] is True
    assert diagnostics["ranking"]["backend_exact_tie_policy"] == "decoupler_2.2_deterministic_seed_0"
    assert diagnostics["ranking"]["backend_tie_break_seed"] == 0
    assert any("internal seed-0 tie break" in warning for warning in diagnostics["warnings"])
    summary = build_ranked_gsea_summary(diagnostics)
    assert any("internal seed-0 tie breaker" in limitation for limitation in summary["limitations"])


def test_ranked_gsea_rejects_unverifiable_kwargs_passthrough(science, tmp_path):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="ranked")
    provenance = _provenance(science, table, universe, purpose="ranked")
    parameters = _base_parameters(provenance, resource_path, purpose="ranked")
    backend = _fake_decoupler(science)
    valid_gsea = backend.mt.gsea

    def opaque_gsea(*, data, net, tmin, raw, empty, verbose, **kwargs):
        return valid_gsea(
            data=data,
            net=net,
            tmin=tmin,
            raw=raw,
            empty=empty,
            verbose=verbose,
            **kwargs,
        )

    backend.mt.gsea = opaque_gsea
    with pytest.raises(RuntimeError, match="cannot verify.*seed.*times"):
        run_ranked_gsea_evidence(table, universe, **parameters, decoupler_module=backend)


def test_generic_ora_complete_ties_share_rank_and_are_flagged(science, tmp_path):
    resource_path = tmp_path / "ties.csv"
    resource_path.write_text(
        "source,target\nA,G1\nA,G2\nB,G1\nB,G2\nC,G6\nC,G7\n",
        encoding="utf-8",
    )
    table, universe = _generic_artifacts(science, purpose="selected")
    provenance = _provenance(science, table, universe, purpose="selected")
    parameters = _base_parameters(provenance, resource_path, purpose="selected")
    parameters["max_p_adjusted"] = 1.0
    evidence, _ = run_generic_ora_evidence(
        table,
        universe,
        **parameters,
        decoupler_module=_fake_decoupler(science),
    )
    tied = evidence.loc[evidence["gene_set"].isin(["A", "B"])]
    assert tied["rank"].nunique() == 1
    assert set(tied["candidate_status"]) == {"tied_best_candidate"}


@pytest.mark.parametrize("version", ["1.9.2", "2.1.9", "2.3.0", "3.0.0"])
def test_enrichment_requires_audited_decoupler_22_family(science, tmp_path, version):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="ranked")
    provenance = _provenance(science, table, universe, purpose="ranked")
    parameters = _base_parameters(provenance, resource_path, purpose="ranked")
    backend = _fake_decoupler(science)
    backend.__version__ = version
    with pytest.raises(RuntimeError, match="requires decoupler 2.2.x"):
        run_ranked_gsea_evidence(table, universe, **parameters, decoupler_module=backend)


def test_enrichment_cell_level_condition_declaration_is_cautious_not_rejected(
    science, tmp_path
):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="ranked", inference_unit="cell")
    provenance = _provenance(science, table, universe, purpose="ranked")
    assert provenance["inference_unit"] == "cell"
    assert provenance["provenance_declarations"]["formal_condition_interpretation_validated"] is False
    assert any("formal Condition interpretation is not validated" in warning for warning in provenance["provenance_warnings"])
    parameters = _base_parameters(provenance, resource_path, purpose="ranked")
    _, diagnostics = run_ranked_gsea_evidence(
        table,
        universe,
        **parameters,
        decoupler_module=_fake_decoupler(science),
    )
    summary = build_ranked_gsea_summary(diagnostics)
    assert summary["key_results"]["provenance"]["formal_condition_interpretation_validated"] is False
    assert any("formal Condition interpretation is not validated" in warning for warning in summary["warnings"])
    json.dumps(summary, ensure_ascii=False, allow_nan=False)


@pytest.mark.parametrize("purpose", ["ranked", "selected"])
def test_custom_producer_and_missing_scientific_metadata_run_with_cautious_generated_parity(
    science, tmp_path, purpose
):
    resource_path = _write_resource(tmp_path / f"{purpose}.csv")
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
    assert provenance["evidence_scope"] == "unspecified"
    assert provenance["provenance_declarations"]["scientific_metadata_status"] == (
        "caller_declared_not_programmatically_verified"
    )
    assert any("caller declarations" in warning for warning in provenance["provenance_warnings"])
    assert any("no canonical producer" in warning for warning in provenance["provenance_warnings"])
    parameters = _base_parameters(provenance, resource_path, purpose=purpose)
    backend = _fake_decoupler(science)
    if purpose == "ranked":
        evidence, diagnostics = run_ranked_gsea_evidence(
            table,
            universe,
            **parameters,
            decoupler_module=backend,
        )
        summary = build_ranked_gsea_summary(diagnostics)
        generated = ranked_gsea_code(
            **parameters,
            expected_resource_sha256=diagnostics["resource"]["sha256"],
        )
        entrypoint = "run_ranked_gsea"
    else:
        evidence, diagnostics = run_generic_ora_evidence(
            table,
            universe,
            **parameters,
            decoupler_module=backend,
        )
        summary = build_generic_ora_summary(diagnostics)
        generated = generic_ora_code(
            **parameters,
            expected_resource_sha256=diagnostics["resource"]["sha256"],
        )
        entrypoint = "run_gene_set_overrepresentation"
    assert summary["key_results"]["provenance"]["scientific_metadata_status"] == (
        "caller_declared_not_programmatically_verified"
    )
    assert any("caller declarations" in warning for warning in summary["warnings"])
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    namespace = {}
    exec(compile(generated, f"<{purpose}-custom-evidence-code>", "exec"), namespace)
    generated_evidence, generated_summary = namespace[entrypoint](
        table.table.copy(),
        universe.table.copy(),
        _fake_decoupler(science),
    )
    science.pd.testing.assert_frame_equal(evidence, generated_evidence)
    assert summary == generated_summary


def test_enrichment_rejects_numeric_and_universe_tampering(science):
    table, universe = _generic_artifacts(science, purpose="ranked")
    tampered_table = replace(table, table=table.table.copy(deep=True))
    tampered_table.table.loc[0, "score"] = 999.0
    with pytest.raises(ValueError, match="current-content fingerprint"):
        _provenance(science, tampered_table, universe, purpose="ranked")
    tampered_universe = replace(universe, table=universe.table.copy(deep=True))
    tampered_universe.table.loc[0, "gene"] = "OTHER"
    with pytest.raises(ValueError, match="tested-universe identity"):
        _provenance(science, table, tampered_universe, purpose="ranked")


def test_enrichment_rejects_cross_paired_artifacts(science):
    table, universe = _generic_artifacts(science, purpose="ranked")
    foreign_parameters = copy.deepcopy(universe.parameters)
    foreign_parameters["ranking_fingerprint"] = "f" * 64
    foreign_source = copy.deepcopy(universe.source)
    foreign_source["parameters"] = copy.deepcopy(foreign_parameters)
    foreign_universe = replace(universe, parameters=foreign_parameters, source=foreign_source)
    with pytest.raises(ValueError, match="disagree on ranking_fingerprint"):
        _provenance(science, table, foreign_universe, purpose="ranked")


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"random_seed": -1}, "random_seed"),
        ({"n_permutations": 1}, "n_permutations"),
        ({"max_output_rows": 1}, "max_output_rows"),
    ],
)
def test_ranked_gsea_guards(science, tmp_path, updates, message):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="ranked")
    provenance = _provenance(science, table, universe, purpose="ranked")
    parameters = _base_parameters(provenance, resource_path, purpose="ranked")
    parameters.update(updates)
    with pytest.raises(ValueError, match=message):
        run_ranked_gsea_evidence(
            table,
            universe,
            **parameters,
            decoupler_module=_fake_decoupler(science),
        )


def test_ranked_gsea_allows_low_permutations_and_seed_zero_with_cautious_generated_parity(science, tmp_path):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="ranked")
    provenance = _provenance(science, table, universe, purpose="ranked")
    parameters = _base_parameters(provenance, resource_path, purpose="ranked")
    parameters.update({"n_permutations": 10, "random_seed": 0})
    evidence, diagnostics = run_ranked_gsea_evidence(
        table,
        universe,
        **parameters,
        decoupler_module=_fake_decoupler(science),
    )
    summary = build_ranked_gsea_summary(diagnostics)
    assert diagnostics["permutation_null_status"] == "not_shuffled_seed_zero"
    assert diagnostics["formal_permutation_inference_valid"] is False
    assert any("Only 10 permutations" in warning for warning in summary["warnings"])
    assert any("disables permutation shuffling" in warning for warning in summary["warnings"])
    assert "computational flags" in summary["results"]

    generated = ranked_gsea_code(
        **parameters,
        expected_resource_sha256=diagnostics["resource"]["sha256"],
    )
    namespace = {}
    exec(compile(generated, "<ranked-gsea-low-permutations>", "exec"), namespace)
    generated_evidence, generated_summary = namespace["run_ranked_gsea"](
        table.table.copy(), universe.table.copy(), _fake_decoupler(science)
    )
    science.pd.testing.assert_frame_equal(evidence, generated_evidence)
    assert generated_summary == summary


def test_ora_rejects_backend_disagreement(science, tmp_path):
    resource_path = _write_resource(tmp_path / "sets.csv")
    table, universe = _generic_artifacts(science, purpose="selected")
    provenance = _provenance(science, table, universe, purpose="selected")
    parameters = _base_parameters(provenance, resource_path, purpose="selected")
    backend = _fake_decoupler(science)
    valid_query = backend.mt.query_set

    def wrong_query(**kwargs):
        result = valid_query(**kwargs)
        result.loc[0, "stat"] += 1.0
        return result

    backend.mt.query_set = wrong_query
    with pytest.raises(RuntimeError, match="log odds disagree"):
        run_generic_ora_evidence(table, universe, **parameters, decoupler_module=backend)


def test_nodes_expose_table_summary_code_and_hash_resource(
    science, comfy_directories, monkeypatch
):
    input_dir, _, _ = comfy_directories
    resource_path = _write_resource(input_dir / "sets.csv")
    fake = _fake_decoupler(science)
    monkeypatch.setitem(sys.modules, "decoupler", fake)
    ranked_table, ranked_universe = _generic_artifacts(science, purpose="ranked")
    ranked_output = _output_values(
        OpenBioSingleCellRankedGSEA.execute(
            ranked_table,
            ranked_universe,
            gene_sets_file=resource_path.name,
            resource_metadata_json=_metadata_json(),
            comparison="treated_vs_control",
            min_targets=2,
            max_targets=4,
            n_permutations=1000,
            random_seed=5,
        )
    )
    assert len(ranked_output) == 3
    assert ranked_output[0].kind == "table"
    assert ranked_output[1].kind == "summary"
    compile(ranked_output[2], "<node-ranked-code>", "exec")

    selected_table, selected_universe = _generic_artifacts(science, purpose="selected")
    ora_output = _output_values(
        OpenBioSingleCellGeneSetOverrepresentation.execute(
            selected_table,
            selected_universe,
            gene_sets_file=resource_path.name,
            resource_metadata_json=_metadata_json(),
            comparison="treated_vs_control",
            min_targets=2,
            min_overlap=2,
        )
    )
    assert len(ora_output) == 3
    assert ora_output[0].kind == "table"
    assert ora_output[1].kind == "summary"
    compile(ora_output[2], "<node-ora-code>", "exec")
    before = OpenBioSingleCellRankedGSEA.fingerprint_inputs(resource_path.name)
    payload = resource_path.read_bytes().replace(b"UP,G1", b"ZZ,G1", 1)
    stat = resource_path.stat()
    resource_path.write_bytes(payload)
    os.utime(resource_path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    after = OpenBioSingleCellRankedGSEA.fingerprint_inputs(resource_path.name)
    assert before != after


def test_enrichment_node_schemas_match_atomic_contract():
    ranked = OpenBioSingleCellRankedGSEA.define_schema()
    ora = OpenBioSingleCellGeneSetOverrepresentation.define_schema()
    assert [output.display_name for output in ranked.outputs] == ["table", "summary", "code"]
    assert [output.display_name for output in ora.outputs] == ["table", "summary", "code"]
    ranked_inputs = [value.id for value in ranked.inputs]
    ora_inputs = [value.id for value in ora.inputs]
    assert ranked_inputs[:5] == [
        "table",
        "universe",
        "gene_sets_file",
        "resource_metadata_json",
        "comparison",
    ]
    assert ora_inputs[:5] == [
        "table",
        "universe",
        "gene_sets_file",
        "resource_metadata_json",
        "comparison",
    ]
    assert "adata" not in ora_inputs
    assert "selection_column" not in ora_inputs
