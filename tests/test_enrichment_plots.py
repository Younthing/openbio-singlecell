from __future__ import annotations

import hashlib
import json
import time
import uuid

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, read_plot, write_anndata, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.contracts import PlotResult, ensure_metadata
from openbio_singlecell.dgidb_resource import _standalone_dgidb_json_sha256
from openbio_singlecell.enrichment_plotting import table_evidence_fingerprint
from openbio_singlecell.nodes_enrichment import (
    OpenBioSingleCellDrugGSEAPlot,
    OpenBioSingleCellDrugHypergeometricPlot,
    OpenBioSingleCellORAEvidencePlot,
    OpenBioSingleCellPathwayScoreContrastPlot,
    OpenBioSingleCellRankedGSEAPlot,
    OpenBioSingleCellScoreActivityPlot,
)
from openbio_singlecell.operations_enrichment import (
    drug_gsea_plot_owned,
    drug_hypergeometric_plot_owned,
    gene_set_ora_plot_owned,
    pathway_score_contrast_plot_owned,
    ranked_gsea_plot_owned,
    score_activity_plot_owned,
)
from openbio_singlecell.pathway_score_contrast import PATHWAY_SCORE_CONTRAST_COLUMNS
from openbio_singlecell.score_artifact import build_score_artifact, store_score_artifact
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _table_result(science, contract: str):
    if contract == "pathway_score_contrast":
        frame = science.pd.DataFrame(
            [
                ["T", "P1", "treated", "control", 4, 4, 40, 38, 1.5, 0.5, 0.3, 0.2, 1.0, 0.4, 1.6, 3.1, 5.2, 0.01, 0.03, 1],
                ["T", "P2", "treated", "control", 4, 4, 40, 38, 0.2, 0.7, 0.4, 0.3, -0.5, -1.0, 0.0, -2.2, 5.8, 0.04, 0.06, 2],
            ],
            columns=PATHWAY_SCORE_CONTRAST_COLUMNS,
        )
        operation = "pathway_score_sample_welch"
    elif contract == "ranked_gsea":
        frame = science.pd.DataFrame(
            [
                ["treated_vs_control", "UP", 1.8, 0.01, 1, True, 10, 8],
                ["treated_vs_control", "DOWN", -1.2, 0.08, 2, False, 9, 7],
            ],
            columns=[
                "comparison", "gene_set", "normalized_enrichment_score", "p_adjusted", "rank",
                "significant", "set_size_before_universe", "set_size_in_universe",
            ],
        )
        operation = "ranked_gsea"
    elif contract == "gene_set_ora":
        frame = science.pd.DataFrame(
            [
                ["selected", "treated_vs_control", "UP", 5, 20, 8, 7, 4, '["G1","G2","G3","G4"]', 4, 3, 1, 12, 2.3, 0.002, 0.01, 1, True, True, True, True, "candidate"],
                ["selected", "treated_vs_control", "DOWN", 5, 20, 7, 6, 1, '["G5"]', 1, 5, 4, 10, -0.4, 0.4, 0.4, 2, False, False, False, False, "nonpositive"],
            ],
            columns=[
                "scope", "comparison", "gene_set", "selected_count", "universe_count",
                "set_size_before_universe", "set_size_in_universe", "overlap_count", "overlap_genes",
                "a", "b", "c", "d", "log_odds_ratio", "p_value", "p_adjusted", "rank",
                "positive_enrichment", "passes_min_overlap", "significant", "candidate", "candidate_status",
            ],
        )
        operation = "gene_set_overrepresentation"
    elif contract == "drug_ora":
        base = _table_result(science, "gene_set_ora").table.rename(columns={"gene_set": "drug"})
        base.insert(9, "resource_sources", ['["DrugBank"]', '["DrugBank"]'])
        frame = base
        operation = "drug_hypergeometric"
    else:
        frame = _table_result(science, "ranked_gsea").table.rename(columns={"gene_set": "drug"})
        frame["matched_targets"] = ['["G1","G2"]', '["G3"]']
        frame["resource_sources"] = ['["DrugBank"]', '["DrugBank"]']
        operation = "drug_gsea"
    parameters = {
        "plot_evidence_schema": "openbio-singlecell/plot-evidence/v1",
        "plot_evidence_fingerprint_sha256": table_evidence_fingerprint(frame, contract=contract),
    }
    return make_table_result(
        table=frame,
        title=contract,
        operation=operation,
        parameters=parameters,
        description="fixture",
        warnings=[],
        input_cells=78,
        input_genes=200,
        started_at=time.perf_counter(),
    )


def _score_adata(science):
    names = [f"cell_{index}" for index in range(6)]
    adata = science.ad.AnnData(
        science.np.ones((6, 3)),
        obs=science.pd.DataFrame(
            {"cell_type": science.pd.Categorical(["T", "T", "T", "B", "B", "B"], categories=["T", "B"])},
            index=names,
        ),
        var=science.pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    scores = science.pd.DataFrame(
        {"P1": [0.9, 0.8, 0.7, 0.2, 0.1, 0.0], "P2": [0.1, 0.2, 0.3, 0.7, 0.8, 0.9]},
        index=names,
    )
    adata.obsm["aucell_scores"] = scores
    adata.obsm["X_umap"] = science.np.asarray(
        [[0.0, 0.0], [0.2, 0.1], [0.4, 0.2], [1.0, 1.0], [1.2, 1.1], [1.4, 1.2]]
    )
    artifact = build_score_artifact(
        frame=scores,
        feature_names=adata.var_names,
        producer_node="OpenBioSingleCellAUCellScores",
        method="AUCell",
        storage="obsm",
        score_key="aucell_scores",
        resource={"sha256": "a" * 64},
        expression={"source": "X", "feature_axis_sha256": "b" * 64},
        parameters={"output_key": "aucell_scores"},
        references=[{"citation": "fixture", "url": "https://example.test", "kind": "method"}],
        np=science.np,
        pd=science.pd,
    )
    store_score_artifact(adata, artifact)
    return adata


def test_enrichment_plot_node_schemas_are_domain_specific_and_closed():
    score = OpenBioSingleCellScoreActivityPlot.define_schema()
    assert score.category == "openbio/single-cell/enrichment"
    assert [item.id for item in score.inputs] == ["adata", "score_key", "view"]
    assert [option.key for option in score.inputs[2].options] == [
        "distributions",
        "grouped_heatmap",
        "embedding",
    ]
    assert [item.display_name for item in score.outputs] == ["plot", "summary", "code"]

    cases = [
        (OpenBioSingleCellPathwayScoreContrastPlot, "Pathway Score Contrast Plot", ["forest", "effect_significance"]),
        (OpenBioSingleCellRankedGSEAPlot, "Ranked GSEA Plot", ["nes_lollipop", "significance_dot"]),
        (OpenBioSingleCellORAEvidencePlot, "ORA Evidence Plot", ["odds_ratio_dot", "overlap_bar"]),
        (OpenBioSingleCellDrugHypergeometricPlot, "Drug Hypergeometric Plot", ["odds_ratio_dot", "overlap_bar"]),
        (OpenBioSingleCellDrugGSEAPlot, "Drug GSEA Plot", ["nes_lollipop", "significance_dot"]),
    ]
    for node, display_name, views in cases:
        schema = node.define_schema()
        assert schema.display_name == display_name
        assert schema.category == "openbio/single-cell/enrichment"
        assert [item.id for item in schema.inputs] == ["table", "view"]
        assert [option.key for option in schema.inputs[1].options] == views
        assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_score_activity_plot_uses_exact_stored_scores_and_group_order(science):
    adata = _score_adata(science)
    snapshot = adata.copy()

    plotted, report, code = score_activity_plot_owned(
        adata,
        score_key="aucell_scores",
        view={"view": "grouped_heatmap", "groupby": "cell_type", "max_scores": 10},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["producer_node"] == "OpenBioSingleCellAUCellScores"
    assert details["method"] == "AUCell"
    assert details["score_names"] == ["P1", "P2"]
    assert details["groups"] == ["T", "B"]
    science.np.testing.assert_allclose(details["group_means"], [[0.8, 0.2], [0.1, 0.8]])
    assert details["plotted_cells"] == 6
    json.dumps(report.summary, allow_nan=False)
    namespace: dict[str, object] = {}
    exec(compile(code, "<score-activity-plot>", "exec"), namespace)
    assert namespace["plot_score_activity"](adata) == plotted.png
    science.pd.testing.assert_frame_equal(adata.obsm["aucell_scores"], snapshot.obsm["aucell_scores"])


@pytest.mark.parametrize(
    ("method", "producer", "storage"),
    [
        ("AUCell", "OpenBioSingleCellAUCellScores", "obsm"),
        ("GSVA", "OpenBioSingleCellGSVAScores", "obsm"),
        ("Scanpy score_genes", "OpenBioSingleCellGenePanelScores", "obs"),
    ],
)
def test_score_activity_plot_accepts_each_exact_score_producer(science, method, producer, storage):
    adata = _score_adata(science)
    original = adata.uns["openbio_singlecell_score_artifacts"].pop("aucell_scores")
    frame = adata.obsm.pop("aucell_scores")
    score_key = "stored_score"
    if storage == "obs":
        frame = frame[["P1"]].rename(columns={"P1": score_key})
        adata.obs[score_key] = frame[score_key]
    else:
        adata.obsm[score_key] = frame
    artifact = build_score_artifact(
        frame=frame,
        feature_names=adata.var_names,
        producer_node=producer,
        method=method,
        storage=storage,
        score_key=score_key,
        resource=original["resource"],
        expression=original["expression"],
        parameters={"output_key": score_key},
        references=original["references"],
        np=science.np,
        pd=science.pd,
    )
    store_score_artifact(adata, artifact)

    _plotted, report, _code = score_activity_plot_owned(adata, score_key=score_key)

    assert report.summary["key_results"]["method"] == method
    assert report.summary["key_results"]["producer_node"] == producer


def test_score_activity_plot_validates_stored_drug_score_history_and_fingerprint(science, tmp_path):
    adata = _score_adata(science)
    del adata.obsm["aucell_scores"]
    del adata.uns["openbio_singlecell_score_artifacts"]
    score_key = "drug_target_score"
    values = [0.1, 0.2, 0.3, 0.8, 0.9, 1.0]
    adata.obs[score_key] = values
    fingerprint = _standalone_dgidb_json_sha256(
        {
            "schema": "openbio-drug-score/v1",
            "observations": adata.obs_names.tolist(),
            "scores": values,
        }
    )
    ensure_metadata(adata)
    adata.uns["openbio_singlecell"]["analysis_history"]["000000"] = {
        "operation": "dgidb_drug_score",
        "parameters": {"output_key": score_key, "score_fingerprint_sha256": fingerprint},
    }
    root = tmp_path / "drug-score"
    root.mkdir()
    write_anndata(root, adata)
    adata = read_anndata(root)

    _plotted, report, _code = score_activity_plot_owned(adata, score_key=score_key)

    details = report.summary["key_results"]
    assert details["producer_node"] == "OpenBioSingleCellDrugScores"
    assert details["method"] == "DGIdb mean target-expression score"
    assert details["evidence_fingerprint_sha256"] == fingerprint


def test_score_activity_embedding_uses_stored_coordinates_and_first_score_by_default(science):
    adata = _score_adata(science)

    plotted, report, _code = score_activity_plot_owned(
        adata,
        score_key="aucell_scores",
        view={"view": "embedding", "embedding_key": "X_umap", "score_name": ""},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["parameters"]["view"]["score_name"] == "P1"
    assert report.summary["key_results"]["score_names"] == ["P1"]


@pytest.mark.parametrize(
    ("owned", "contract", "view", "function_name", "mark_key", "expected"),
    [
        (pathway_score_contrast_plot_owned, "pathway_score_contrast", {"view": "forest", "max_pathways": 20}, "plot_pathway_score_contrast", "mean_differences", [1.0, -0.5]),
        (ranked_gsea_plot_owned, "ranked_gsea", {"view": "nes_lollipop", "max_gene_sets": 20}, "plot_ranked_gsea", "normalized_enrichment_scores", [1.8, -1.2]),
        (gene_set_ora_plot_owned, "gene_set_ora", {"view": "odds_ratio_dot", "max_gene_sets": 20}, "plot_ora_evidence", "log_odds_ratios", [2.3, -0.4]),
        (drug_hypergeometric_plot_owned, "drug_ora", {"view": "overlap_bar", "max_drugs": 20}, "plot_drug_hypergeometric", "overlap_counts", [4, 1]),
        (drug_gsea_plot_owned, "drug_gsea", {"view": "significance_dot", "max_drugs": 20}, "plot_drug_gsea", "normalized_enrichment_scores", [1.8, -1.2]),
    ],
)
def test_enrichment_table_plots_report_exact_visual_marks_and_code_parity(
    science, owned, contract, view, function_name, mark_key, expected
):
    result = _table_result(science, contract)
    snapshot = result.table.copy(deep=True)

    plotted, report, code = owned(result, view=view)

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details[mark_key] == expected
    assert details["plotted_items"] == 2
    assert details["result_fingerprint_sha256"] == result.parameters["plot_evidence_fingerprint_sha256"]
    namespace: dict[str, object] = {}
    exec(compile(code, f"<{function_name}>", "exec"), namespace)
    assert namespace[function_name](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, snapshot)


@pytest.mark.parametrize(
    ("case", "message"),
    [("producer", "producer"), ("schema", "canonical columns"), ("fingerprint", "fingerprint")],
)
def test_enrichment_table_plot_rejects_wrong_or_tampered_contract(science, case, message):
    result = _table_result(science, "ranked_gsea")
    if case == "producer":
        result.source["operation"] = "other"
    elif case == "schema":
        result.table.drop(columns="p_adjusted", inplace=True)
    else:
        result.table.loc[0, "normalized_enrichment_score"] = 9.0

    with pytest.raises((TypeError, ValueError), match=message):
        ranked_gsea_plot_owned(result)


@pytest.mark.parametrize(
    ("contract", "operation_name", "view"),
    [
        ("pathway_score_contrast", "pathway_score_contrast_plot_operation", {"view": "forest", "max_pathways": 20}),
        ("ranked_gsea", "ranked_gsea_plot_operation", {"view": "nes_lollipop", "max_gene_sets": 20}),
        ("gene_set_ora", "gene_set_ora_plot_operation", {"view": "overlap_bar", "max_gene_sets": 20}),
        ("drug_ora", "drug_hypergeometric_plot_operation", {"view": "odds_ratio_dot", "max_drugs": 20}),
        ("drug_gsea", "drug_gsea_plot_operation", {"view": "significance_dot", "max_drugs": 20}),
    ],
)
def test_enrichment_plot_workers_round_trip_each_table_adapter(
    tmp_path, science, contract, operation_name, view
):
    import openbio_singlecell.operations_enrichment as operations

    result = _table_result(science, contract)
    input_root = tmp_path / contract
    input_root.mkdir()
    write_table(input_root, result.table, result_metadata(result))
    before = {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}
    staging = tmp_path / f"{contract}-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = getattr(operations, operation_name)(
        context,
        {"table": {"type": "artifact", "kind": "OPENBIO_SINGLE_CELL_TABLE", "codec": "table-jsonl-v1", "path": str(input_root.resolve())}},
        {"view": view},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["source"]["operation"].endswith("_plot")
    assert before == {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}


def test_score_activity_plot_worker_reads_anndata_without_rewriting_it(tmp_path, science):
    from openbio_singlecell.operations_enrichment import score_activity_plot_operation

    adata = _score_adata(science)
    input_root = tmp_path / "adata"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "score-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = score_activity_plot_operation(
        context,
        {"adata": {"type": "artifact", "kind": "OPENBIO_ANNDATA", "codec": "anndata-h5ad-v1", "path": str(input_root.resolve())}},
        {"score_key": "aucell_scores", "view": {"view": "distributions", "max_scores": 20}},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    namespace: dict[str, object] = {}
    exec(compile(records[2]["value"], "<worker-score-plot-code>", "exec"), namespace)
    assert namespace["plot_score_activity"](read_anndata(input_root)) == png
    assert hashlib.sha256(input_path.read_bytes()).digest() == before
