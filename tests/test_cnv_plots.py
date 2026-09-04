from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import read_plot, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.cnv_analysis import (
    CNV_STATE_SCHEMA,
    CNV_STATE_TYPE,
    CNVState,
    _cnv_artifact_fingerprints,
    _cnv_artifact_hash,
    _cnv_score_table_fingerprint,
)
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.node_types import CNVStateType, PlotResultType, TableResultType
from openbio_singlecell.nodes_cnv import OpenBioSingleCellCNVHeatmapPlot, OpenBioSingleCellCNVScorePlot
from openbio_singlecell.operations_cnv import (
    cnv_heatmap_plot,
    cnv_heatmap_plot_owned,
    cnv_score_plot,
    cnv_score_plot_owned,
)
from openbio_singlecell.staged_state_codec import CNV_STATE_CODEC, write_cnv_state
from openbio_singlecell.worker_protocol import OperationContext

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _cnv_state(science):
    np = science.np
    obs_names = [f"cell_{index}" for index in range(6)]
    var_names = [f"gene_{index}" for index in range(4)]
    adata = science.ad.AnnData(
        np.arange(24, dtype=float).reshape(6, 4) / 10,
        obs=science.pd.DataFrame(
            {
                "cell_type": science.pd.Categorical(
                    ["reference"] * 3 + ["query"] * 3,
                    categories=["reference", "query"],
                    ordered=True,
                ),
                "sample": ["S1", "S1", "S2", "S1", "S2", "S2"],
            },
            index=obs_names,
        ),
        var=science.pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1", "chr2", "chr2"],
                "start": science.np.asarray([0, 100, 0, 100], dtype="int64"),
                "end": science.np.asarray([99, 199, 99, 199], dtype="int64"),
            },
            index=var_names,
        ),
    )
    matrix = np.asarray(
        [
            [-0.1, 0.0, 0.1, 0.0],
            [-0.2, 0.0, 0.2, 0.1],
            [-0.1, 0.1, 0.1, 0.0],
            [0.5, 0.6, -0.4, -0.5],
            [0.4, 0.5, -0.3, -0.4],
            [0.6, 0.7, -0.5, -0.6],
        ],
        dtype=float,
    )
    adata.obsm["X_cnv"] = matrix
    adata.uns["cnv"] = {"chr_pos": {"chr1": 0, "chr2": 2}}
    adata.uns["openbio_cnv_state"] = {
        "schema": CNV_STATE_SCHEMA,
        "output_key": "cnv",
        "genome_assembly": "GRCh38",
        "reference_key": "cell_type",
        "reference_categories": ["reference"],
        "sample_key": "sample",
        "windows": [
            {"chromosome": "chr1", "genes": 2, "first_window_index": 0, "window_count": 2},
            {"chromosome": "chr2", "genes": 2, "first_window_index": 2, "window_count": 2},
        ],
        "infercnvpy_version": "0.6.1",
    }
    metadata = {
        "schema": CNV_STATE_SCHEMA,
        "artifact_type": CNV_STATE_TYPE,
        "producer_node_id": "OpenBioSingleCellInferCNV",
        "stage": "inferred",
        "source_kind": "X",
        "layer_name": None,
        "output_key": "cnv",
        "reference_key": "cell_type",
        "sample_key": "sample",
        "genome_assembly": "GRCh38",
        "parameters": {"window_size": 2, "step": 1},
        "fingerprints": {},
        "artifact_fingerprint_sha256": "",
    }
    metadata["fingerprints"] = _cnv_artifact_fingerprints(
        adata,
        metadata,
        numpy=science.np,
        scipy_sparse=science.sparse,
    )
    metadata["artifact_fingerprint_sha256"] = _cnv_artifact_hash(metadata, metadata["fingerprints"])
    return CNVState(adata, metadata)


def _oversized_sparse_cnv_state(science):
    cells, windows = 10_000, 501
    adata = science.ad.AnnData(
        science.sparse.csr_matrix((cells, 1), dtype=float),
        obs=science.pd.DataFrame(
            {
                "cell_type": science.pd.Categorical(["reference"] * cells),
                "sample": ["S1"] * cells,
            },
            index=[f"cell_{index}" for index in range(cells)],
        ),
        var=science.pd.DataFrame(
            {"chromosome": ["chr1"], "start": [0], "end": [99]},
            index=["gene_0"],
        ),
    )
    adata.obsm["X_cnv"] = science.sparse.csr_matrix((cells, windows), dtype=float)
    adata.uns["cnv"] = {"chr_pos": {"chr1": 0}}
    adata.uns["openbio_cnv_state"] = {
        "schema": CNV_STATE_SCHEMA,
        "output_key": "cnv",
        "genome_assembly": "GRCh38",
        "reference_key": "cell_type",
        "reference_categories": ["reference"],
        "sample_key": "sample",
        "windows": [
            {"chromosome": "chr1", "genes": 1, "first_window_index": 0, "window_count": windows}
        ],
        "infercnvpy_version": "0.6.1",
    }
    metadata = {
        "schema": CNV_STATE_SCHEMA,
        "artifact_type": CNV_STATE_TYPE,
        "producer_node_id": "OpenBioSingleCellInferCNV",
        "stage": "inferred",
        "source_kind": "X",
        "layer_name": None,
        "output_key": "cnv",
        "reference_key": "cell_type",
        "sample_key": "sample",
        "genome_assembly": "GRCh38",
        "parameters": {"window_size": 1, "step": 1},
        "fingerprints": {},
        "artifact_fingerprint_sha256": "",
    }
    metadata["fingerprints"] = _cnv_artifact_fingerprints(
        adata,
        metadata,
        numpy=science.np,
        scipy_sparse=science.sparse,
    )
    metadata["artifact_fingerprint_sha256"] = _cnv_artifact_hash(metadata, metadata["fingerprints"])
    return CNVState(adata, metadata)


def _score_result(science):
    table = science.pd.DataFrame(
        {
            "rank": [1, 2],
            "group_type": ["string", "string"],
            "group_value": ["query", "reference"],
            "group_display": ["query", "reference"],
            "cell_count": [3, 3],
            "cnv_score": [0.5, 0.1],
            "absolute_min": [0.3, 0.0],
            "absolute_q1": [0.4, 0.05],
            "absolute_median": [0.5, 0.1],
            "absolute_mean": [0.5, 0.1],
            "absolute_q3": [0.6, 0.15],
            "absolute_max": [0.7, 0.2],
        }
    )
    parameters = {
        "groupby": "cell_type",
        "use_rep": "cnv",
        "output_key": "cnv_score",
        "overwrite_existing": False,
        "inplace": True,
        "independent_validation_rtol": 1e-6,
        "independent_validation_atol": 1e-8,
        "input_cnv_state_fingerprint_sha256": "a" * 64,
        "input_cnv_matrix_fingerprint_sha256": "b" * 64,
        "partition_fingerprint_sha256": "c" * 64,
        "output_fingerprint_sha256": "d" * 64,
        "table_fingerprint_sha256": _cnv_score_table_fingerprint(table),
    }
    return make_table_result(
        table=table,
        title="CNV group scores by cell_type",
        operation="cnv_score",
        parameters=parameters,
        description="Descriptive group scores.",
        warnings=[],
        input_cells=6,
        input_genes=4,
        started_at=time.perf_counter(),
    )


def _context(root: Path) -> OperationContext:
    root.mkdir()
    return OperationContext.from_request_path(root / "request.json", str(uuid.uuid4()))


def _descriptor(root: Path, *, kind: str, codec: str) -> dict[str, object]:
    return {"type": "artifact", "path": str(root.resolve()), "kind": kind, "codec": codec}


def _snapshot(root: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in root.iterdir() if path.is_file()}


def test_cnv_plot_schemas_are_typed_domain_companions():
    heatmap = OpenBioSingleCellCNVHeatmapPlot.define_schema()
    score = OpenBioSingleCellCNVScorePlot.define_schema()

    assert heatmap.node_id == "OpenBioSingleCellCNVHeatmapPlot"
    assert heatmap.display_name == "CNV Heatmap Plot"
    assert heatmap.category == "openbio/single-cell/copy-number"
    assert [item.id for item in heatmap.inputs] == ["cnv_state", "groupby", "view"]
    assert heatmap.inputs[0].io_type == CNVStateType.io_type
    assert [option.key for option in heatmap.inputs[2].options] == ["group_mean", "cells"]
    assert [[item.id for item in option.inputs] for option in heatmap.inputs[2].options] == [
        ["max_groups"],
        ["max_cells"],
    ]
    assert [item.io_type for item in heatmap.outputs] == [
        PlotResultType.io_type,
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]

    assert score.node_id == "OpenBioSingleCellCNVScorePlot"
    assert score.display_name == "CNV Score Plot"
    assert score.category == "openbio/single-cell/copy-number"
    assert [item.id for item in score.inputs] == ["table", "view"]
    assert score.inputs[0].io_type == TableResultType.io_type
    assert [option.key for option in score.inputs[1].options] == ["score_bar", "absolute_interval"]
    assert [[item.id for item in option.inputs] for option in score.inputs[1].options] == [
        ["max_groups"],
        ["max_groups"],
    ]
    assert [item.io_type for item in score.outputs] == [
        PlotResultType.io_type,
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]


@pytest.mark.parametrize(
    "view,expected_rows",
    [
        ({"view": "group_mean", "max_groups": 50}, 2),
        ({"view": "cells", "max_cells": 100}, 6),
    ],
)
def test_cnv_heatmap_plot_uses_verified_genome_windows_without_mutation(science, view, expected_rows):
    state = _cnv_state(science)
    before = state.to_adata()

    plotted, report, code = cnv_heatmap_plot_owned(state, groupby="", view=view)

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["groupby"] == "cell_type"
    assert details["row_count"] == expected_rows
    assert details["window_count"] == 4
    assert details["chromosomes"] == ["chr1", "chr2"]
    assert details["chromosome_start_indices"] == [0, 2]
    assert details["genome_assembly"] == "GRCh38"
    assert details["cnv_matrix_fingerprint_sha256"] == state.metadata["fingerprints"]["cnv_matrix_sha256"]
    assert "tumor" in " ".join(report.summary["limitations"]).lower()
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<cnv-heatmap-plot-code>", "exec"), namespace)
    assert namespace["plot_cnv_heatmap"](state) == plotted.png
    after = state.to_adata()
    science.np.testing.assert_array_equal(after.obsm["X_cnv"], before.obsm["X_cnv"])
    science.pd.testing.assert_frame_equal(after.obs, before.obs)
    assert after.uns["cnv"] == before.uns["cnv"]


@pytest.mark.parametrize("view", ["score_bar", "absolute_interval"])
def test_cnv_score_plot_preserves_descriptive_group_axis_and_intervals(science, view):
    result = _score_result(science)
    before = result.table.copy(deep=True)

    plotted, report, code = cnv_score_plot_owned(
        result,
        view={"view": view, "max_groups": 30},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["groups"] == ["query", "reference"]
    assert details["ranks"] == [1, 2]
    assert details["cnv_scores"] == pytest.approx([0.5, 0.1])
    assert details["absolute_q1"] == pytest.approx([0.4, 0.05])
    assert details["absolute_q3"] == pytest.approx([0.6, 0.15])
    assert details["table_fingerprint_sha256"] == result.parameters["table_fingerprint_sha256"]
    assert "tumor" in " ".join(report.summary["limitations"]).lower()
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<cnv-score-plot-code>", "exec"), namespace)
    assert namespace["plot_cnv_score"](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, before)


def test_cnv_plot_workers_round_trip_typed_inputs_without_rewriting_them(science, tmp_path):
    state_root = tmp_path / "state"
    state_root.mkdir()
    write_cnv_state(state_root, _cnv_state(science))
    state_before = _snapshot(state_root)
    heatmap_context = _context(tmp_path / "heatmap")
    heatmap_records = cnv_heatmap_plot(
        heatmap_context,
        {"cnv_state": _descriptor(state_root, kind="OPENBIO_CNV_STATE", codec=CNV_STATE_CODEC)},
        {"groupby": "", "view": {"view": "group_mean", "max_groups": 50}},
    )

    assert [record["name"] for record in heatmap_records] == ["plot", "summary", "code"]
    assert read_plot(heatmap_context.output_root / heatmap_records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert _snapshot(state_root) == state_before

    score_result = _score_result(science)
    table_root = tmp_path / "table"
    table_root.mkdir()
    write_table(table_root, score_result.table, result_metadata(score_result))
    table_before = _snapshot(table_root)
    score_context = _context(tmp_path / "score")
    score_records = cnv_score_plot(
        score_context,
        {"table": _descriptor(table_root, kind="OPENBIO_SINGLE_CELL_TABLE", codec="table-jsonl-v1")},
        {"view": {"view": "absolute_interval", "max_groups": 30}},
    )

    assert [record["name"] for record in score_records] == ["plot", "summary", "code"]
    assert read_plot(score_context.output_root / score_records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert _snapshot(table_root) == table_before


def test_cnv_plots_reject_fingerprint_tampering_and_wrong_producer(science):
    state = _cnv_state(science)
    state._adata.obsm["X_cnv"][0, 0] += 0.25
    with pytest.raises(ValueError, match="fingerprint"):
        cnv_heatmap_plot_owned(state, view={"view": "group_mean", "max_groups": 50})

    result = _score_result(science)
    result.table.loc[0, "cnv_score"] = 0.9
    with pytest.raises(ValueError, match="fingerprint"):
        cnv_score_plot_owned(result, view={"view": "score_bar", "max_groups": 30})

    wrong = _score_result(science)
    wrong.source["operation"] = "pseudobulk_edger"
    with pytest.raises(ValueError, match="producer|Group Score"):
        cnv_score_plot_owned(wrong, view={"view": "score_bar", "max_groups": 30})


def test_cnv_heatmap_rejects_semantically_invalid_rehashed_genome_windows(science):
    state = _cnv_state(science)
    adata = state.to_adata()
    metadata = state.metadata
    adata.uns["openbio_cnv_state"]["windows"][1]["first_window_index"] = 1
    metadata["fingerprints"] = _cnv_artifact_fingerprints(
        adata,
        metadata,
        numpy=science.np,
        scipy_sparse=science.sparse,
    )
    metadata["artifact_fingerprint_sha256"] = _cnv_artifact_hash(metadata, metadata["fingerprints"])
    forged = CNVState(adata, metadata)

    with pytest.raises(ValueError, match="contiguous"):
        cnv_heatmap_plot_owned(forged, view={"view": "group_mean", "max_groups": 50})


def test_cnv_heatmap_preflights_displayed_values_before_sparse_densification(science, monkeypatch):
    state = _oversized_sparse_cnv_state(science)

    def reject_toarray(*_args, **_kwargs):
        raise AssertionError("sparse matrix was densified before display-size preflight")

    monkeypatch.setattr(science.sparse.csr_matrix, "toarray", reject_toarray)
    with pytest.raises(ValueError, match="5,000,000 displayed"):
        cnv_heatmap_plot_owned(state, view={"view": "cells", "max_cells": 10_000})
