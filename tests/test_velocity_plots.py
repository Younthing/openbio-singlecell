from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import read_plot, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.node_types import PlotResultType, TableResultType, VelocityStateType
from openbio_singlecell.nodes_velocity import (
    OpenBioSingleCellVelocityDynamicsPlot,
    OpenBioSingleCellVelocityGeneRankingPlot,
)
from openbio_singlecell.operations_velocity import (
    velocity_dynamics_plot,
    velocity_dynamics_plot_owned,
    velocity_gene_ranking_plot,
    velocity_gene_ranking_plot_owned,
)
from openbio_singlecell.staged_state_codec import VELOCITY_STATE_CODEC, write_velocity_state
from openbio_singlecell.velocity_analysis import _build_velocity_state, run_velocity_ranking
from openbio_singlecell.velocity_portable import _vs_series_fingerprint, _vs_summary, _vs_write_state
from openbio_singlecell.worker_protocol import OperationContext

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _recovered_state(science):
    np = science.np
    n_obs, n_vars = 8, 3
    obs_names = [f"cell_{index}" for index in range(n_obs)]
    var_names = ["G1", "G2", "G3"]
    adata = science.ad.AnnData(
        np.ones((n_obs, n_vars), dtype=np.float32),
        obs=science.pd.DataFrame(index=obs_names),
        var=science.pd.DataFrame(index=var_names),
    )
    adata.layers["spliced"] = np.arange(n_obs * n_vars, dtype=np.float32).reshape(n_obs, n_vars) / 5 + 1
    adata.layers["unspliced"] = np.arange(n_obs * n_vars, dtype=np.float32).reshape(n_obs, n_vars) / 10 + 0.5
    adata.layers["Ms"] = adata.layers["spliced"].copy()
    adata.layers["Mu"] = adata.layers["unspliced"].copy()
    for key in ("fit_t", "fit_tau", "fit_tau_"):
        adata.layers[key] = np.tile(np.linspace(0.0, 1.0, n_obs)[:, None], (1, n_vars))
    adata.var["openbio_dynamics_selected"] = [True, True, True]
    adata.var["openbio_dynamics_fit_success"] = [True, True, False]
    adata.var["openbio_dynamics_velocity_gene"] = [True, False, True]
    adata.var["fit_likelihood"] = [0.8, 1.2, np.nan]
    adata.var["fit_r2"] = [0.4, 0.7, np.nan]
    for index, key in enumerate(
        (
            "fit_alpha",
            "fit_beta",
            "fit_gamma",
            "fit_t_",
            "fit_scaling",
            "fit_std_u",
            "fit_std_s",
            "fit_u0",
            "fit_s0",
            "fit_pval_steady",
            "fit_steady_u",
            "fit_steady_s",
            "fit_variance",
            "fit_alignment_scaling",
        ),
        start=1,
    ):
        adata.var[key] = [index / 10, index / 5, np.nan]
    adata.varm["loss"] = np.asarray([[0.8, 0.4, 0.2], [0.7, 0.3, 0.1], [np.nan, np.nan, np.nan]])
    adata.uns["recover_dynamics"] = {
        "fit_connected_states": True,
        "fit_basal_transcription": False,
        "use_raw": False,
    }
    adjacency = science.sparse.diags([np.ones(n_obs - 1), np.ones(n_obs - 1)], [-1, 1], shape=(n_obs, n_obs))
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"method": "fixture"},
    }
    adata.obsp["connectivities"] = adjacency.tocsr()
    adata.obsp["distances"] = adjacency.tocsr()
    parameters = {"gene_selection": "all", "max_iter": 3}
    portable = _vs_write_state(
        adata,
        stage="dynamics_recovered",
        producer_node_id="OpenBioSingleCellRecoverDynamics",
        scvelo_version="0.3.4",
        parameters=parameters,
        parent=None,
        neighbors_key="neighbors",
        vkey=None,
        velocity_mode=None,
        has_dynamics=True,
    )
    summary = _vs_summary(
        node_id="OpenBioSingleCellRecoverDynamics",
        status="exploratory_dynamical_parameters_recovered",
        methods="Recovered fixture dynamics.",
        results="Two fits succeeded.",
        key_results={"state_fingerprint_sha256": portable["state_fingerprint_sha256"]},
        parameters=parameters,
        references=[{"citation": "fixture", "url": "https://example.test", "kind": "method"}],
        software_versions={"scvelo": "0.3.4"},
        report_warnings=[],
        limitations=["Model-fit evidence only."],
    )
    return _build_velocity_state(adata, summary)


def _ranking_result(science):
    state = _recovered_state(science)
    table, summary = run_velocity_ranking(state, top_n=10, include_failed=True, max_output_rows=100)
    table = table.convert_dtypes()
    evidence = summary["key_results"]
    parameters = {
        **summary["parameters"],
        "dynamics_fit_fingerprint_sha256": evidence["dynamics_fit_fingerprint_sha256"],
        "upstream_state_fingerprint_sha256": evidence["upstream_state_fingerprint_sha256"],
        "table_fingerprint_sha256": _vs_series_fingerprint(
            table,
            list(table.columns),
            label="dynamics-ranking-plot",
        ),
    }
    return make_table_result(
        table=table,
        title="Recovered dynamics fit ranking",
        operation="rank_recovered_dynamics",
        parameters=parameters,
        description=summary["results"],
        warnings=[],
        input_cells=8,
        input_genes=3,
        started_at=time.perf_counter(),
    )


def _context(root: Path) -> OperationContext:
    root.mkdir()
    return OperationContext.from_request_path(root / "request.json", str(uuid.uuid4()))


def _descriptor(root: Path, *, kind: str, codec: str) -> dict[str, object]:
    return {"type": "artifact", "path": str(root.resolve()), "kind": kind, "codec": codec}


def _snapshot(root: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in root.iterdir() if path.is_file()}


def test_velocity_plot_schemas_are_typed_read_only_companions():
    dynamics = OpenBioSingleCellVelocityDynamicsPlot.define_schema()
    ranking = OpenBioSingleCellVelocityGeneRankingPlot.define_schema()

    assert dynamics.node_id == "OpenBioSingleCellVelocityDynamicsPlot"
    assert dynamics.display_name == "Velocity Dynamics Plot"
    assert dynamics.category == "openbio/single-cell/velocity"
    assert [item.id for item in dynamics.inputs] == ["velocity_state", "gene", "view"]
    assert dynamics.inputs[0].io_type == VelocityStateType.io_type
    assert [option.key for option in dynamics.inputs[2].options] == ["phase_portrait", "fit_loss"]
    assert [item.io_type for item in dynamics.outputs] == [PlotResultType.io_type, "OPENBIO_SINGLE_CELL_SUMMARY", "STRING"]

    assert ranking.node_id == "OpenBioSingleCellVelocityGeneRankingPlot"
    assert ranking.display_name == "Velocity Gene Ranking Plot"
    assert ranking.category == "openbio/single-cell/velocity"
    assert [item.id for item in ranking.inputs] == ["table", "view"]
    assert ranking.inputs[0].io_type == TableResultType.io_type
    assert [option.key for option in ranking.inputs[1].options] == ["fit_quality", "kinetic_parameters"]
    assert [[item.id for item in option.inputs] for option in ranking.inputs[1].options] == [
        ["max_genes"],
        ["max_genes"],
    ]
    assert [item.io_type for item in ranking.outputs] == [PlotResultType.io_type, "OPENBIO_SINGLE_CELL_SUMMARY", "STRING"]


@pytest.mark.parametrize("view", ["phase_portrait", "fit_loss"])
def test_velocity_dynamics_plot_uses_verified_stored_fit_without_mutation(science, view):
    state = _recovered_state(science)
    before = state.portable_adata()

    plotted, report, code = velocity_dynamics_plot_owned(
        state,
        gene="",
        view={"view": view},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["gene"] == "G2"
    assert details["gene_selection"] == "highest_stored_fit_likelihood"
    assert details["fit_likelihood"] == pytest.approx(1.2)
    assert details["dynamics_fit_fingerprint_sha256"] == state.portable_adata().uns[
        "openbio_velocity_state"
    ]["dynamics_fit_fingerprint_sha256"]
    assert "fate" in " ".join(report.summary["limitations"]).lower()
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<velocity-dynamics-plot-code>", "exec"), namespace)
    assert namespace["plot_velocity_dynamics"](state.portable_adata()) == plotted.png
    after = state.portable_adata()
    science.np.testing.assert_array_equal(after.layers["Ms"], before.layers["Ms"])
    science.np.testing.assert_array_equal(after.layers["fit_t"], before.layers["fit_t"])
    science.pd.testing.assert_frame_equal(after.var, before.var)


@pytest.mark.parametrize("view", ["fit_quality", "kinetic_parameters"])
def test_velocity_gene_ranking_plot_preserves_upstream_rank_and_kinetic_evidence(science, view):
    result = _ranking_result(science)
    before = result.table.copy(deep=True)

    plotted, report, code = velocity_gene_ranking_plot_owned(
        result,
        view={"view": view, "max_genes": 30},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["genes"] == ["G2", "G1"]
    assert details["ranks"] == [1, 2]
    assert details["fit_likelihood"] == pytest.approx([1.2, 0.8])
    assert details["table_fingerprint_sha256"] == result.parameters["table_fingerprint_sha256"]
    assert details["dynamics_fit_fingerprint_sha256"] == result.parameters[
        "dynamics_fit_fingerprint_sha256"
    ]
    assert "p-value" in " ".join(report.summary["limitations"])
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<velocity-gene-ranking-plot-code>", "exec"), namespace)
    assert namespace["plot_velocity_gene_ranking"](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, before)


def test_velocity_plot_workers_round_trip_typed_inputs_without_rewriting_them(science, tmp_path):
    state_root = tmp_path / "state"
    state_root.mkdir()
    write_velocity_state(state_root, _recovered_state(science))
    state_before = _snapshot(state_root)
    dynamics_context = _context(tmp_path / "dynamics")

    dynamics_records = velocity_dynamics_plot(
        dynamics_context,
        {"velocity_state": _descriptor(state_root, kind="OPENBIO_VELOCITY_STATE", codec=VELOCITY_STATE_CODEC)},
        {"gene": "G1", "view": {"view": "phase_portrait"}},
    )

    assert [record["name"] for record in dynamics_records] == ["plot", "summary", "code"]
    assert read_plot(dynamics_context.output_root / dynamics_records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert _snapshot(state_root) == state_before

    ranking_result = _ranking_result(science)
    table_root = tmp_path / "table"
    table_root.mkdir()
    write_table(table_root, ranking_result.table, result_metadata(ranking_result))
    table_before = _snapshot(table_root)
    ranking_context = _context(tmp_path / "ranking")
    ranking_records = velocity_gene_ranking_plot(
        ranking_context,
        {"table": _descriptor(table_root, kind="OPENBIO_SINGLE_CELL_TABLE", codec="table-jsonl-v1")},
        {"view": {"view": "fit_quality", "max_genes": 20}},
    )

    assert [record["name"] for record in ranking_records] == ["plot", "summary", "code"]
    assert read_plot(ranking_context.output_root / ranking_records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert _snapshot(table_root) == table_before


def test_velocity_plots_reject_tampered_state_table_and_wrong_producer(science):
    state = _recovered_state(science)
    state._adata.layers["fit_t"][0, 0] += 0.25
    with pytest.raises(ValueError, match="fingerprint|current-content"):
        velocity_dynamics_plot_owned(state, gene="G1", view={"view": "phase_portrait"})

    result = _ranking_result(science)
    result.table.loc[0, "fit_likelihood"] = 9.0
    with pytest.raises(ValueError, match="fingerprint"):
        velocity_gene_ranking_plot_owned(result, view={"view": "fit_quality", "max_genes": 20})

    wrong = _ranking_result(science)
    wrong.source["operation"] = "marker_genes"
    with pytest.raises(ValueError, match="producer|Ranking producer"):
        velocity_gene_ranking_plot_owned(wrong, view={"view": "fit_quality", "max_genes": 20})
