from __future__ import annotations

import hashlib
import json
import time
import uuid

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.contracts import PlotResult, ensure_metadata
from openbio_singlecell.nodes_factorization import (
    OpenBioSingleCellCNMFProgramsPlot,
    OpenBioSingleCellCNMFRankPlot,
)
from openbio_singlecell.operations_factorization import (
    cnmf_programs_plot,
    cnmf_programs_plot_owned,
    cnmf_rank_plot,
    cnmf_rank_plot_owned,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _rank_table_result(science):
    table = science.pd.DataFrame(
        {
            "k": [2, 3, 4],
            "stability": [0.81, 0.76, 0.68],
            "prediction_error": [12.0, 9.0, 8.5],
            "statistics_density_threshold": [2.0, 2.0, 2.0],
            "expected_restarts": [4, 4, 4],
            "completed_restarts": [4, 4, 4],
            "combined_components": [8, 12, 16],
        }
    )
    parameters = {
        "source": "layer",
        "layer_name": "counts",
        "candidate_ks": [2, 3, 4],
        "components_min": 2,
        "components_max": 4,
        "n_iter": 4,
        "total_restarts": 12,
        "source_features": 100,
        "input_fingerprint": "sha256:" + "a" * 64,
        "k_metrics_fingerprint_sha256": _rank_metrics_fingerprint(table),
        "fixed_policy": {"statistics_density_threshold": 2.0},
    }
    return make_table_result(
        table=table,
        title="cNMF rank-survey metrics",
        operation="cnmf_rank_survey",
        parameters=parameters,
        description="Complete cNMF K diagnostics; no K was selected.",
        warnings=[],
        input_cells=20,
        input_genes=100,
        started_at=time.perf_counter(),
        random_seed=17,
    )


def _rank_metrics_fingerprint(table):
    rows = [
        [
            int(row.k),
            float(row.stability),
            float(row.prediction_error),
            float(row.statistics_density_threshold),
            int(row.expected_restarts),
            int(row.completed_restarts),
            int(row.combined_components),
        ]
        for row in table.itertuples(index=False)
    ]
    payload = json.dumps(rows, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return f"sha256:{hashlib.sha256(payload).hexdigest()}"


def _axis_fingerprint(values):
    digest = hashlib.sha256()
    for value in values:
        encoded = str(value).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
    return f"sha256:{digest.hexdigest()}"


def _matrix_fingerprint(matrix, rows, columns, science):
    values = science.np.asarray(matrix)
    digest = hashlib.sha256()
    digest.update(str(tuple(int(value) for value in values.shape)).encode("ascii"))
    digest.update(str(values.dtype).encode("ascii"))
    for axis in (rows, columns):
        for value in axis:
            encoded = str(value).encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "little"))
            digest.update(encoded)
    digest.update(science.np.ascontiguousarray(values).tobytes())
    return f"sha256:{digest.hexdigest()}"


def _program_adata(science):
    obs_names = [f"cell_{index}" for index in range(6)]
    var_names = [f"gene_{index}" for index in range(5)]
    programs = ["cNMF_1", "cNMF_2"]
    usage = science.pd.DataFrame(
        [[0.9, 0.1], [0.8, 0.2], [0.7, 0.3], [0.3, 0.7], [0.2, 0.8], [0.1, 0.9]],
        index=obs_names,
        columns=programs,
    )
    scores = science.pd.DataFrame(
        [[5.0, 1.0], [4.0, 2.0], [3.0, 5.0], [2.0, 4.0], [1.0, 3.0]],
        index=var_names,
        columns=programs,
    )
    top_genes = {"cNMF_1": ["gene_0", "gene_1", "gene_2"], "cNMF_2": ["gene_2", "gene_3", "gene_4"]}
    adata = science.ad.AnnData(
        science.np.ones((len(obs_names), len(var_names))),
        obs=science.pd.DataFrame(index=obs_names),
        var=science.pd.DataFrame(index=var_names),
    )
    adata.obsm["X_cnmf_usage"] = usage
    adata.varm["cnmf_gep_scores"] = scores
    adata.varm["cnmf_gep_tpm"] = scores.copy()
    adata.uns["cnmf"] = {
        "schema_version": 2,
        "adapter_version": 4,
        "selected_k": 2,
        "program_names": programs,
        "top_genes": top_genes,
        "storage_keys": {
            "usage": "obsm:X_cnmf_usage",
            "gep_scores": "varm:cnmf_gep_scores",
            "gep_tpm": "varm:cnmf_gep_tpm",
        },
        "observation_axis_fingerprint_sha256": _axis_fingerprint(obs_names),
        "feature_axis_fingerprint_sha256": _axis_fingerprint(var_names),
        "usage_fingerprint_sha256": _matrix_fingerprint(usage.to_numpy(), obs_names, programs, science),
        "gep_scores_fingerprint_sha256": _matrix_fingerprint(
            scores.to_numpy(), var_names, programs, science
        ),
        "top_genes_fingerprint_sha256": "sha256:"
        + hashlib.sha256(
            json.dumps(top_genes, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
        ).hexdigest(),
    }
    ensure_metadata(adata, display_name="cNMF programs", source={"kind": "test"})
    adata.uns["openbio_singlecell"]["analysis_history"]["000000"] = {
        "operation": "cnmf_consensus_programs",
        "parameters": {"selected_k": 2},
        "input_cells": 6,
        "input_genes": 5,
        "random_seed": 17,
        "elapsed_seconds": 0.1,
        "timestamp": "2026-01-01T00:00:00+00:00",
    }
    return adata


def test_cnmf_rank_plot_schema_consumes_exact_rank_table_contract():
    schema = OpenBioSingleCellCNMFRankPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellCNMFRankPlot"
    assert schema.display_name == "cNMF Rank Plot"
    assert schema.category == "openbio/single-cell/factorization"
    assert [item.id for item in schema.inputs] == ["k_metrics", "view"]
    assert [option.key for option in schema.inputs[1].options] == [
        "stability_prediction_error",
        "restart_completion",
    ]
    assert [option.inputs for option in schema.inputs[1].options] == [[], []]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_cnmf_programs_plot_schema_has_closed_usage_and_top_gene_views():
    schema = OpenBioSingleCellCNMFProgramsPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellCNMFProgramsPlot"
    assert schema.display_name == "cNMF Programs Plot"
    assert schema.category == "openbio/single-cell/factorization"
    assert [item.id for item in schema.inputs] == ["adata", "view"]
    assert [option.key for option in schema.inputs[1].options] == [
        "usage_distributions",
        "top_genes",
    ]
    assert [[item.id for item in option.inputs] for option in schema.inputs[1].options] == [
        [],
        ["genes_per_program"],
    ]
    assert schema.inputs[1].options[1].inputs[0].id == "genes_per_program"
    assert schema.inputs[1].options[1].inputs[0].default == 10
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_cnmf_rank_plot_renders_complete_stability_and_error_evidence_without_selecting_k(science):
    result = _rank_table_result(science)
    snapshot = result.table.copy(deep=True)

    plotted, report, code = cnmf_rank_plot_owned(
        result,
        view={"view": "stability_prediction_error"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellCNMFRankPlot"
    assert summary["parameters"] == {"view": {"view": "stability_prediction_error"}}
    assert summary["key_results"]["candidate_ks"] == [2, 3, 4]
    assert summary["key_results"]["plotted_ranks"] == 3
    assert summary["key_results"]["stability"] == [0.81, 0.76, 0.68]
    assert summary["key_results"]["prediction_error"] == [12.0, 9.0, 8.5]
    assert "selected_k" not in summary["key_results"]
    assert "best" not in summary["results"].lower()
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<cnmf-rank-plot-code>", "exec"), namespace)
    assert namespace["plot_cnmf_rank"](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, snapshot)


def test_cnmf_program_usage_plot_validates_stored_program_evidence_and_code_matches(science):
    adata = _program_adata(science)
    snapshot = adata.copy()

    plotted, report, code = cnmf_programs_plot_owned(
        adata,
        view={"view": "usage_distributions"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellCNMFProgramsPlot"
    assert summary["parameters"] == {"view": {"view": "usage_distributions"}}
    assert summary["key_results"]["programs"] == ["cNMF_1", "cNMF_2"]
    assert summary["key_results"]["plotted_programs"] == 2
    assert summary["key_results"]["plotted_cells"] == 6
    assert summary["key_results"]["usage_means"] == pytest.approx([0.5, 0.5])
    assert "cluster" not in summary["results"].lower()
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<cnmf-programs-plot-code>", "exec"), namespace)
    assert namespace["plot_cnmf_programs"](adata) == plotted.png
    science.np.testing.assert_array_equal(adata.X, snapshot.X)
    science.pd.testing.assert_frame_equal(adata.obsm["X_cnmf_usage"], snapshot.obsm["X_cnmf_usage"])
    science.pd.testing.assert_frame_equal(adata.varm["cnmf_gep_scores"], snapshot.varm["cnmf_gep_scores"])


def test_cnmf_program_top_gene_view_uses_stored_order_and_scores(science):
    adata = _program_adata(science)

    plotted, report, code = cnmf_programs_plot_owned(
        adata,
        view={"view": "top_genes", "genes_per_program": 2},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["view"] == "top_genes"
    assert details["plotted_top_genes"] == 4
    assert details["top_gene_evidence"] == [
        {"program": "cNMF_1", "genes": ["gene_0", "gene_1"], "gep_scores": [5.0, 4.0]},
        {"program": "cNMF_2", "genes": ["gene_2", "gene_3"], "gep_scores": [5.0, 4.0]},
    ]
    namespace: dict[str, object] = {}
    exec(compile(code, "<cnmf-top-genes-plot-code>", "exec"), namespace)
    assert namespace["plot_cnmf_programs"](adata) == plotted.png


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("producer", "Rank Survey"),
        ("schema", "canonical columns"),
        ("k_axis", "K values"),
        ("restart", "complete restart"),
        ("density", "fixed producer policy"),
        ("metric_content", "metrics fingerprint"),
        ("fingerprint", "input fingerprint"),
    ],
)
def test_cnmf_rank_plot_rejects_wrong_or_tampered_table_contract(science, case, message):
    result = _rank_table_result(science)
    if case == "producer":
        result.source["operation"] = "other"
    elif case == "schema":
        result.table.drop(columns="prediction_error", inplace=True)
    elif case == "k_axis":
        result.table.loc[1, "k"] = 4
    elif case == "restart":
        result.table.loc[0, "completed_restarts"] = 3
    elif case == "density":
        result.table["statistics_density_threshold"] = 1.5
    elif case == "metric_content":
        result.table.loc[0, "stability"] = 0.79
    else:
        result.parameters["input_fingerprint"] = "not-a-fingerprint"
        result.source["parameters"]["input_fingerprint"] = "not-a-fingerprint"

    with pytest.raises((TypeError, ValueError), match=message):
        cnmf_rank_plot_owned(result)


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("observation_axis", "observation axis fingerprint"),
        ("feature_axis", "feature axis fingerprint"),
        ("usage", "usage fingerprint"),
        ("scores", "GEP score fingerprint"),
        ("top_genes", "top-gene fingerprint"),
        ("history", "provenance"),
        ("schema", "complete stored cNMF"),
    ],
)
def test_cnmf_programs_plot_rejects_tampered_evidence_before_rendering(science, case, message):
    adata = _program_adata(science)
    if case == "observation_axis":
        adata = adata[[1, 0, 2, 3, 4, 5]].copy()
    elif case == "feature_axis":
        adata = adata[:, [1, 0, 2, 3, 4]].copy()
    elif case == "usage":
        adata.obsm["X_cnmf_usage"].iloc[0] = [0.8, 0.2]
    elif case == "scores":
        adata.varm["cnmf_gep_scores"].iloc[0, 0] += 1.0
    elif case == "top_genes":
        adata.uns["cnmf"]["top_genes"]["cNMF_1"] = ["gene_1", "gene_0", "gene_2"]
    elif case == "history":
        adata.uns["openbio_singlecell"]["analysis_history"] = {}
    else:
        del adata.uns["cnmf"]["usage_fingerprint_sha256"]

    with pytest.raises(ValueError, match=message):
        cnmf_programs_plot_owned(adata)


def test_factorization_plots_restore_matplotlib_state(science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    original = matplotlib.rcParams["lines.linewidth"]
    try:
        matplotlib.rcParams["lines.linewidth"] = 3.25
        rc_value = matplotlib.rcParams["lines.linewidth"]
        figure_numbers = pyplot.get_fignums()

        cnmf_rank_plot_owned(_rank_table_result(science), view={"view": "restart_completion"})
        cnmf_programs_plot_owned(
            _program_adata(science), view={"view": "top_genes", "genes_per_program": 2}
        )

        assert matplotlib.rcParams["lines.linewidth"] == rc_value
        assert pyplot.get_fignums() == figure_numbers
    finally:
        matplotlib.rcParams["lines.linewidth"] = original


def test_cnmf_rank_plot_worker_round_trips_png_without_rewriting_table(tmp_path, science):
    result = _rank_table_result(science)
    input_root = tmp_path / "rank"
    input_root.mkdir()
    write_table(input_root, result.table, result_metadata(result))
    before = {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}
    staging = tmp_path / "rank-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cnmf_rank_plot(
        context,
        {
            "k_metrics": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"view": {"view": "stability_prediction_error"}},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"view": {"view": "stability_prediction_error"}}
    assert before == {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}


@pytest.mark.parametrize(
    "view",
    [
        {"view": "usage_distributions"},
        {"view": "top_genes", "genes_per_program": 2},
    ],
)
def test_cnmf_programs_plot_worker_round_trips_png_without_rewriting_anndata(tmp_path, science, view):
    adata = _program_adata(science)
    input_root = tmp_path / "programs"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "programs-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cnmf_programs_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
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
    assert hashlib.sha256(input_path.read_bytes()).digest() == before
