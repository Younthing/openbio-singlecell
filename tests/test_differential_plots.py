from __future__ import annotations

import copy
import hashlib
import json
import time
import uuid

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import read_plot, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.differential_evidence import differential_table_content_fingerprint
from openbio_singlecell.nodes_differential import (
    OpenBioSingleCellPseudobulkConditionContrastPlot,
    OpenBioSingleCellPseudobulkQCPlot,
    OpenBioSingleCellSCVIPopulationDEEvidencePlot,
)
from openbio_singlecell.operations_differential import (
    pseudobulk_condition_contrast_plot,
    pseudobulk_condition_contrast_plot_owned,
    pseudobulk_qc_plot,
    pseudobulk_qc_plot_owned,
    run_pseudobulk_owned,
    scvi_population_de_evidence_plot,
    scvi_population_de_evidence_plot_owned,
)
from openbio_singlecell.pseudobulk_artifact_codec import PSEUDOBULK_CODEC, PSEUDOBULK_KIND, write_pseudobulk
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse, registered_operation_ids

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _pseudobulk_fixture(science):
    samples = ["control_1", "control_2", "treated_1", "treated_2"]
    rows = []
    counts = []
    for sample_index, sample in enumerate(samples):
        condition = "control" if sample.startswith("control") else "treated"
        for population_index, population in enumerate(("T", "B")):
            rows.append({"sample": sample, "cell_type": population, "condition": condition})
            counts.append(
                [
                    2 + sample_index + population_index,
                    5 + 2 * sample_index,
                    3 + population_index,
                    8 - population_index,
                ]
            )
    adata = science.ad.AnnData(
        science.np.asarray(counts, dtype=science.np.int64),
        obs=science.pd.DataFrame(rows, index=[f"cell_{index}" for index in range(len(rows))]),
        var=science.pd.DataFrame(index=["g1", "g2", "g3", "g4"]),
    )
    adata.uns["openbio_singlecell"] = {
        "schema_version": 1,
        "version": "0.2.0",
        "display_name": "fixture",
        "source": {},
        "random_seed": 0,
        "warnings": [],
        "analysis_history": {},
        "annotations": {
            "cell_type": {
                "schema_version": 1,
                "operation": "map_cluster_annotations",
                "annotation_status": "curated",
                "curation_assertion": "caller_declared",
                "output_column": "cell_type",
            }
        },
    }
    artifact, _, _ = run_pseudobulk_owned(
        adata,
        sample_key="sample",
        population_key="cell_type",
        condition_key="condition",
        source={"source": "X"},
        min_cells=1,
        min_counts=1,
    )
    return artifact


def test_differential_plot_nodes_expose_separate_scientific_contracts():
    pseudobulk = OpenBioSingleCellPseudobulkQCPlot.define_schema()
    contrast = OpenBioSingleCellPseudobulkConditionContrastPlot.define_schema()
    scvi = OpenBioSingleCellSCVIPopulationDEEvidencePlot.define_schema()

    assert pseudobulk.node_id == "OpenBioSingleCellPseudobulkQCPlot"
    assert pseudobulk.display_name == "Pseudobulk QC Plot"
    assert pseudobulk.category == "openbio/single-cell/differential-expression"
    assert [item.id for item in pseudobulk.inputs] == ["pseudobulk", "view"]
    assert [option.key for option in pseudobulk.inputs[1].options] == ["profile_qc"]

    assert contrast.node_id == "OpenBioSingleCellPseudobulkConditionContrastPlot"
    assert contrast.display_name == "Pseudobulk Condition Contrast Plot"
    assert [item.id for item in contrast.inputs] == ["table", "view"]
    assert [option.key for option in contrast.inputs[1].options] == ["volcano", "ma", "top_genes"]
    assert [item.id for item in contrast.inputs[1].options[2].inputs] == ["n_genes"]

    assert scvi.node_id == "OpenBioSingleCellSCVIPopulationDEEvidencePlot"
    assert scvi.display_name == "scVI Population DE Evidence Plot"
    assert [item.id for item in scvi.inputs] == ["table", "view"]
    assert [option.key for option in scvi.inputs[1].options] == [
        "probability_lfc",
        "bayes_factor",
        "top_genes",
    ]
    assert [item.display_name for item in pseudobulk.outputs] == ["plot", "summary", "code"]
    assert [item.display_name for item in contrast.outputs] == ["plot", "summary", "code"]
    assert [item.display_name for item in scvi.outputs] == ["plot", "summary", "code"]


def test_pseudobulk_qc_plot_uses_validated_profiles_without_mutating_artifact(science):
    artifact = _pseudobulk_fixture(science)
    before = artifact.to_adata()
    metadata_before = copy.deepcopy(artifact.metadata)

    plotted, report, code = pseudobulk_qc_plot_owned(artifact, view={"view": "profile_qc"})

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellPseudobulkQCPlot"
    assert summary["parameters"] == {"view": {"view": "profile_qc"}}
    assert summary["key_results"]["plotted_profiles"] == 8
    assert summary["key_results"]["samples"] == ["control_1", "control_2", "treated_1", "treated_2"]
    assert summary["key_results"]["populations"] == ["T", "B"]
    assert "Sample-level Condition inference" not in summary["results"]
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<pseudobulk-qc-plot-code>", "exec"), namespace)
    assert namespace["plot_pseudobulk_qc"](artifact) == plotted.png
    science.np.testing.assert_array_equal(artifact.to_adata().X, before.X)
    science.pd.testing.assert_frame_equal(artifact.to_adata().obs, before.obs)
    assert artifact.metadata == metadata_before


def _edger_result(science):
    table = science.pd.DataFrame(
        {
            "gene": ["g_up", "g_down", "g_flat"],
            "log2_fold_change": [2.0, -1.5, 0.2],
            "log_counts_per_million": [6.0, 5.0, 7.0],
            "quasi_likelihood_f": [12.0, 8.0, 0.2],
            "p_value": [0.001, 0.02, 0.4],
            "p_adjusted": [0.003, 0.03, 0.4],
        }
    )
    parameters = {
        "engine": "edgeR",
        "population": "T",
        "condition_key": "condition",
        "reference_condition": "control",
        "comparison_condition": "treated",
        "fdr_threshold": 0.05,
        "min_abs_log2_fold_change": 1.0,
        "fixed_policy": {"lfc_shrinkage": None},
        "table_content_fingerprint_sha256": differential_table_content_fingerprint(
            table,
            contract="pseudobulk_edger",
        ),
    }
    return make_table_result(
        table=table,
        title="edgeR QL: treated vs control",
        operation="pseudobulk_edger",
        parameters=parameters,
        description="Sample-level Condition contrast.",
        warnings=[],
        input_cells=800,
        input_genes=3,
        started_at=time.perf_counter(),
    )


def test_pseudobulk_condition_contrast_plot_preserves_frequentist_sample_semantics(science):
    result = _edger_result(science)
    before = result.table.copy(deep=True)

    plotted, report, code = pseudobulk_condition_contrast_plot_owned(
        result,
        view={"view": "volcano"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellPseudobulkConditionContrastPlot"
    assert summary["parameters"] == {"view": {"view": "volcano"}}
    assert summary["key_results"]["engine"] == "edgeR"
    assert summary["key_results"]["reference_condition"] == "control"
    assert summary["key_results"]["comparison_condition"] == "treated"
    assert summary["key_results"]["significant_up"] == 1
    assert summary["key_results"]["significant_down"] == 1
    assert summary["key_results"]["color_legend"] == {
        "Significant: comparison higher": "#D64A4A",
        "Significant: comparison lower": "#3D6FB6",
        "Not called significant": "#B8B8B8",
    }
    assert "independent Sample" in summary["methods"]
    assert "current-content fingerprint" in summary["methods"]
    assert "posterior" not in summary["results"].lower()
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<pseudobulk-condition-plot-code>", "exec"), namespace)
    assert namespace["plot_pseudobulk_condition_contrast"](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, before)


def _scvi_change_result(science):
    table = science.pd.DataFrame(
        {
            "gene": ["g_up", "g_down", "g_flat"],
            "comparison": ["activated vs resting"] * 3,
            "group1": ["activated"] * 3,
            "group2": ["resting"] * 3,
            "proba_de": [0.98, 0.93, 0.55],
            "proba_not_de": [0.02, 0.07, 0.45],
            "bayes_factor": [3.89, 2.59, 0.20],
            "scale1": [0.4, 0.1, 0.2],
            "scale2": [0.1, 0.4, 0.2],
            "pseudocounts": [0.1, 0.1, 0.1],
            "delta": [0.25, 0.25, 0.25],
            "lfc_mean": [1.8, -1.3, 0.05],
            "lfc_median": [1.7, -1.2, 0.04],
            "lfc_std": [0.3, 0.4, 0.2],
            "lfc_min": [0.8, -2.0, -0.4],
            "lfc_max": [2.4, -0.5, 0.5],
            "raw_mean1": [6.0, 2.0, 3.0],
            "raw_mean2": [2.0, 5.0, 3.0],
            "non_zeros_proportion1": [0.9, 0.5, 0.7],
            "non_zeros_proportion2": [0.5, 0.9, 0.7],
            "raw_normalized_mean1": [2.0, 0.5, 1.0],
            "raw_normalized_mean2": [0.5, 2.0, 1.0],
            "is_de_fdr": [True, True, False],
        }
    )
    parameters = {
        "groupby": "state",
        "group1": "activated",
        "group2": "resting",
        "population_scope": {"population_scope": "all"},
        "mode": "change",
        "delta": 0.25,
        "fdr_target": 0.05,
        "batch_handling": "shared_technical_batches",
        "n_samples_overall": 5000,
        "random_seed": 23,
        "table_content_fingerprint_sha256": differential_table_content_fingerprint(
            table,
            contract="scvi_model_de_evidence/change",
        ),
    }
    return make_table_result(
        table=table,
        title="scVI model evidence: activated vs resting",
        operation="scvi_model_de_evidence",
        parameters=parameters,
        description="Cell-level exploratory model evidence.",
        warnings=["Cells are not independent biological Samples."],
        input_cells=500,
        input_genes=3,
        started_at=time.perf_counter(),
        random_seed=23,
    )


def test_scvi_population_plot_keeps_posterior_evidence_distinct_from_condition_inference(science):
    result = _scvi_change_result(science)
    before = result.table.copy(deep=True)

    plotted, report, code = scvi_population_de_evidence_plot_owned(
        result,
        view={"view": "probability_lfc"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellSCVIPopulationDEEvidencePlot"
    assert summary["parameters"] == {"view": {"view": "probability_lfc"}}
    assert summary["key_results"]["mode"] == "change"
    assert summary["key_results"]["posterior_fdr_tagged"] == 2
    assert summary["key_results"]["statistical_unit"] == "cells and fitted-model posterior draws"
    assert "current-content fingerprint" in summary["methods"]
    assert "not Sample-level Condition inference" in summary["results"]
    assert "p-value" not in summary["results"]
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<scvi-population-plot-code>", "exec"), namespace)
    assert namespace["plot_scvi_population_de_evidence"](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, before)


def test_pseudobulk_qc_plot_does_not_offer_or_run_new_projection_analyses(science):
    artifact = _pseudobulk_fixture(science)

    plotted, _, code = pseudobulk_qc_plot_owned(artifact, view={"view": "profile_qc"})

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert "log1p" not in code
    assert "linalg.svd" not in code
    assert "corrcoef" not in code
    for unsupported in ("sample_pca", "sample_correlation"):
        with pytest.raises(ValueError, match="closed supported"):
            pseudobulk_qc_plot_owned(artifact, view={"view": unsupported})


def _deseq2_result(science):
    table = science.pd.DataFrame(
        {
            "gene": ["g_up", "g_down", "g_filtered", "g_flat"],
            "base_mean": [50.0, 30.0, 0.0, 80.0],
            "log2_fold_change": [2.0, -1.5, science.np.nan, 0.2],
            "log2_fold_change_standard_error": [0.3, 0.4, science.np.nan, 0.5],
            "wald_statistic": [6.0, -4.0, science.np.nan, 0.4],
            "p_value": [0.001, 0.02, science.np.nan, 0.4],
            "p_adjusted": [0.002, 0.02, science.np.nan, science.np.nan],
        }
    )
    parameters = {
        "engine": "PyDESeq2",
        "population": "T",
        "condition_key": "condition",
        "reference_condition": "control",
        "comparison_condition": "treated",
        "fdr_threshold": 0.05,
        "min_abs_log2_fold_change": 1.0,
        "table_content_fingerprint_sha256": differential_table_content_fingerprint(
            table,
            contract="pseudobulk_deseq2",
        ),
    }
    return make_table_result(
        table=table,
        title="PyDESeq2: treated vs control",
        operation="pseudobulk_deseq2",
        parameters=parameters,
        description="Sample-level Condition contrast.",
        warnings=[],
        input_cells=800,
        input_genes=4,
        started_at=time.perf_counter(),
    )


@pytest.mark.parametrize(
    ("build", "view"),
    [
        (_edger_result, {"view": "ma"}),
        (_edger_result, {"view": "top_genes", "n_genes": 2}),
        (_deseq2_result, {"view": "volcano"}),
        (_deseq2_result, {"view": "ma"}),
    ],
)
def test_pseudobulk_contrast_views_support_both_frequentist_engines(science, build, view):
    result = build(science)

    plotted, report, code = pseudobulk_condition_contrast_plot_owned(result, view=view)

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["key_results"]["view"] == view["view"]
    namespace: dict[str, object] = {}
    exec(compile(code, "<pseudobulk-contrast-view-code>", "exec"), namespace)
    assert namespace["plot_pseudobulk_condition_contrast"](result.table) == plotted.png


def _scvi_vanilla_result(science):
    table = science.pd.DataFrame(
        {
            "gene": ["g_up", "g_down", "g_flat"],
            "comparison": ["activated vs resting"] * 3,
            "group1": ["activated"] * 3,
            "group2": ["resting"] * 3,
            "proba_m1": [0.95, 0.1, 0.5],
            "proba_m2": [0.05, 0.9, 0.5],
            "bayes_factor": [2.9, 2.2, 0.0],
            "scale1": [0.4, 0.1, 0.2],
            "scale2": [0.1, 0.4, 0.2],
            "raw_mean1": [6.0, 2.0, 3.0],
            "raw_mean2": [2.0, 5.0, 3.0],
            "non_zeros_proportion1": [0.9, 0.5, 0.7],
            "non_zeros_proportion2": [0.5, 0.9, 0.7],
            "raw_normalized_mean1": [2.0, 0.5, 1.0],
            "raw_normalized_mean2": [0.5, 2.0, 1.0],
        }
    )
    parameters = {
        "groupby": "state",
        "group1": "activated",
        "group2": "resting",
        "population_scope": {"population_scope": "all"},
        "mode": "vanilla",
        "delta": None,
        "fdr_target": None,
        "batch_handling": "observed_technical_batches",
        "n_samples_overall": 5000,
        "random_seed": 23,
        "table_content_fingerprint_sha256": differential_table_content_fingerprint(
            table,
            contract="scvi_model_de_evidence/vanilla",
        ),
    }
    return make_table_result(
        table=table,
        title="scVI model evidence: activated vs resting",
        operation="scvi_model_de_evidence",
        parameters=parameters,
        description="Cell-level exploratory model evidence.",
        warnings=[],
        input_cells=500,
        input_genes=3,
        started_at=time.perf_counter(),
        random_seed=23,
    )


@pytest.mark.parametrize(
    ("build", "view"),
    [
        (_scvi_change_result, {"view": "bayes_factor"}),
        (_scvi_change_result, {"view": "top_genes", "n_genes": 2}),
        (_scvi_vanilla_result, {"view": "probability_lfc"}),
        (_scvi_vanilla_result, {"view": "top_genes", "n_genes": 2}),
    ],
)
def test_scvi_population_evidence_views_cover_change_and_vanilla_modes(science, build, view):
    result = build(science)

    plotted, report, code = scvi_population_de_evidence_plot_owned(result, view=view)

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["key_results"]["view"] == view["view"]
    namespace: dict[str, object] = {}
    exec(compile(code, "<scvi-population-view-code>", "exec"), namespace)
    assert namespace["plot_scvi_population_de_evidence"](result.table) == plotted.png


def test_differential_table_plots_reject_wrong_producers_and_tampered_evidence(science):
    wrong = _edger_result(science)
    wrong.source["operation"] = "marker_genes"
    with pytest.raises(ValueError, match="edgeR or PyDESeq2"):
        pseudobulk_condition_contrast_plot_owned(wrong)

    malformed = _edger_result(science)
    malformed.table.loc[0, "p_adjusted"] = 0.9
    with pytest.raises(ValueError, match="Benjamini"):
        pseudobulk_condition_contrast_plot_owned(malformed)

    inconsistent = _scvi_change_result(science)
    inconsistent.table.loc[0, "proba_not_de"] = 0.5
    with pytest.raises(ValueError, match="complementary"):
        scvi_population_de_evidence_plot_owned(inconsistent)

    with pytest.raises(ValueError, match="closed supported"):
        pseudobulk_qc_plot_owned(_pseudobulk_fixture(science), view={"view": "profile_qc", "extra": 1})


def test_pseudobulk_edger_plot_rejects_valid_range_content_tampering(science):
    result = _edger_result(science)
    _, _, code = pseudobulk_condition_contrast_plot_owned(result)
    result.table.loc[0, "log_counts_per_million"] = 6.25

    with pytest.raises(ValueError, match="current-content fingerprint"):
        pseudobulk_condition_contrast_plot_owned(result)
    namespace: dict[str, object] = {}
    exec(compile(code, "<pseudobulk-edger-tamper-code>", "exec"), namespace)
    with pytest.raises(ValueError, match="current-content fingerprint"):
        namespace["plot_pseudobulk_condition_contrast"](result.table)


def test_pseudobulk_pydeseq2_plot_rejects_valid_range_content_tampering(science):
    result = _deseq2_result(science)
    _, _, code = pseudobulk_condition_contrast_plot_owned(result)
    result.table.loc[0, "base_mean"] = 50.25

    with pytest.raises(ValueError, match="current-content fingerprint"):
        pseudobulk_condition_contrast_plot_owned(result)
    namespace: dict[str, object] = {}
    exec(compile(code, "<pseudobulk-pydeseq2-tamper-code>", "exec"), namespace)
    with pytest.raises(ValueError, match="current-content fingerprint"):
        namespace["plot_pseudobulk_condition_contrast"](result.table)


def test_scvi_population_plot_rejects_valid_range_content_tampering(science):
    result = _scvi_change_result(science)
    _, _, code = scvi_population_de_evidence_plot_owned(result)
    result.table.loc[0, "scale1"] = 0.45

    with pytest.raises(ValueError, match="current-content fingerprint"):
        scvi_population_de_evidence_plot_owned(result)
    namespace: dict[str, object] = {}
    exec(compile(code, "<scvi-de-tamper-code>", "exec"), namespace)
    with pytest.raises(ValueError, match="current-content fingerprint"):
        namespace["plot_scvi_population_de_evidence"](result.table)


def _tree_hash(root):
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).digest()
        for path in root.rglob("*")
        if path.is_file()
    }


def test_pseudobulk_qc_plot_worker_reads_typed_artifact_without_rewriting_it(tmp_path, science):
    artifact = _pseudobulk_fixture(science)
    input_root = tmp_path / "pseudobulk"
    input_root.mkdir()
    write_pseudobulk(input_root, artifact)
    before = _tree_hash(input_root)
    staging = tmp_path / "plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = pseudobulk_qc_plot(
        context,
        {
            "pseudobulk": {
                "type": "artifact",
                "kind": PSEUDOBULK_KIND,
                "codec": PSEUDOBULK_CODEC,
                "path": str(input_root.resolve()),
            }
        },
        {"view": {"view": "profile_qc"}},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"view": {"view": "profile_qc"}}
    assert _tree_hash(input_root) == before


@pytest.mark.parametrize(
    ("operation", "build", "view"),
    [
        (pseudobulk_condition_contrast_plot, _edger_result, {"view": "volcano"}),
        (pseudobulk_condition_contrast_plot, _deseq2_result, {"view": "ma"}),
        (scvi_population_de_evidence_plot, _scvi_change_result, {"view": "probability_lfc"}),
        (scvi_population_de_evidence_plot, _scvi_vanilla_result, {"view": "bayes_factor"}),
    ],
)
def test_differential_table_plot_workers_round_trip_png_without_rewriting_table(
    tmp_path,
    science,
    operation,
    build,
    view,
):
    result = build(science)
    input_root = tmp_path / operation.__name__
    input_root.mkdir()
    encoded = result.table.copy()
    for column in encoded.columns:
        if science.pd.api.types.is_float_dtype(encoded[column].dtype) and bool(encoded[column].isna().any()):
            encoded[column] = encoded[column].astype(science.pd.Float64Dtype())
    write_table(input_root, encoded, result_metadata(result))
    before = _tree_hash(input_root)
    staging = tmp_path / f"{operation.__name__}.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = operation(
        context,
        {
            "table": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"view": view},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"view": view}
    assert _tree_hash(input_root) == before


def test_differential_plot_worker_operation_ids_are_registered():
    assert {
        "openbio.node.pseudobulkqcplot",
        "openbio.node.pseudobulkconditioncontrastplot",
        "openbio.node.scvipopulationdeevidenceplot",
    }.issubset(registered_operation_ids())


def test_differential_plots_restore_matplotlib_and_numpy_global_state(science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    artifact = _pseudobulk_fixture(science)
    original_linewidth = matplotlib.rcParams["lines.linewidth"]
    try:
        matplotlib.rcParams["lines.linewidth"] = 3.25
        expected_linewidth = matplotlib.rcParams["lines.linewidth"]
        expected_figures = pyplot.get_fignums()
        expected_numpy = science.np.random.get_state()

        pseudobulk_qc_plot_owned(artifact)
        pseudobulk_condition_contrast_plot_owned(_edger_result(science))
        scvi_population_de_evidence_plot_owned(_scvi_change_result(science))

        assert matplotlib.rcParams["lines.linewidth"] == expected_linewidth
        assert pyplot.get_fignums() == expected_figures
        observed_numpy = science.np.random.get_state()
        assert observed_numpy[0] == expected_numpy[0]
        science.np.testing.assert_array_equal(observed_numpy[1], expected_numpy[1])
        assert observed_numpy[2:] == expected_numpy[2:]
    finally:
        matplotlib.rcParams["lines.linewidth"] = original_linewidth
