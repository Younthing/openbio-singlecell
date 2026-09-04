from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata
from openbio_singlecell.cell_cycle import _cell_cycle_axis_fingerprint, _cell_cycle_score_fingerprint
from openbio_singlecell.nodes_trajectory import OpenBioSingleCellCellCycleScorePlot
from openbio_singlecell.operations_trajectory import cell_cycle_score_plot, cell_cycle_score_plot_owned
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _scored_adata(science):
    observation_ids = [f"cell_{index}" for index in range(9)]
    s_score = [-1.0, -0.5, -0.1, 0.8, 0.6, 0.4, 0.1, 0.2, 0.2]
    g2m_score = [-1.0, -0.2, -0.3, 0.2, 0.1, 0.2, 0.8, 0.7, 0.5]
    phase = ["G1", "G1", "G1", "S", "S", "S", "G2M", "G2M", "G2M"]
    adata = science.ad.AnnData(
        science.np.ones((9, 2)),
        obs=science.pd.DataFrame(
            {
                "cell_cycle_s_score": s_score,
                "cell_cycle_g2m_score": g2m_score,
                "cell_cycle_phase": science.pd.Categorical(
                    phase,
                    categories=["G1", "S", "G2M"],
                ),
            },
            index=observation_ids,
        ),
        var=science.pd.DataFrame(index=["G0", "G1"]),
    )
    adata.uns["openbio_singlecell"] = {
        "schema_version": 1,
        "analysis_history": {
            "000000": {
                "operation": "cell_cycle_score",
                "parameters": {
                    "source_kind": "layer",
                    "layer_name": "log1p_norm",
                    "gene_set_source": "regev_human_97",
                    "organism": "human",
                    "output_prefix": "cell_cycle",
                    "observation_axis_fingerprint_sha256": _cell_cycle_axis_fingerprint(observation_ids),
                    "score_bundle_fingerprint_sha256": _cell_cycle_score_fingerprint(
                        observation_ids,
                        s_score,
                        g2m_score,
                        phase,
                    ),
                },
            }
        },
    }
    return adata


def _replace_s_scores(adata):
    adata.obs.loc[:, "cell_cycle_s_score"] = 0.0


def test_cell_cycle_score_plot_schema_is_a_read_only_annotation_companion():
    schema = OpenBioSingleCellCellCycleScorePlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellCellCycleScorePlot"
    assert schema.display_name == "Cell Cycle Score Plot"
    assert schema.category == "openbio/single-cell/annotation"
    assert [item.id for item in schema.inputs] == ["adata", "output_prefix"]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_cell_cycle_score_plot_renders_stored_score_plane_and_phase_counts(science):
    adata = _scored_adata(science)
    before = adata.copy()

    plotted, report, code = cell_cycle_score_plot_owned(adata, output_prefix="cell_cycle")

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["phase_order"] == ["G1", "S", "G2M"]
    assert details["phase_counts"] == {"G1": 3, "S": 3, "G2M": 3}
    assert details["plotted_cells"] == 9
    assert details["score_columns"] == {
        "s_score": "cell_cycle_s_score",
        "g2m_score": "cell_cycle_g2m_score",
        "phase": "cell_cycle_phase",
    }
    assert "Condition inference" in " ".join(report.summary["limitations"])
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<cell-cycle-score-plot-code>", "exec"), namespace)
    assert namespace["plot_cell_cycle_scores"](adata) == plotted.png
    science.pd.testing.assert_frame_equal(adata.obs, before.obs)
    assert adata.uns == before.uns


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (_replace_s_scores, "score bundle fingerprint"),
        (lambda adata: setattr(adata, "obs_names", list(reversed(adata.obs_names))), "axis fingerprint"),
        (
            lambda adata: adata.uns["openbio_singlecell"]["analysis_history"]["000000"].__setitem__(
                "operation", "dpt"
            ),
            "producer history",
        ),
    ],
)
def test_cell_cycle_score_plot_rejects_tampered_or_wrong_evidence(science, mutation, message):
    adata = _scored_adata(science)
    mutation(adata)

    with pytest.raises((TypeError, ValueError), match=message):
        cell_cycle_score_plot_owned(adata)


def test_cell_cycle_score_plot_rejects_phase_inconsistent_with_scores(science):
    adata = _scored_adata(science)
    changed = [
        "S" if phase == "G1" else "G1" if phase == "S" else phase
        for phase in adata.obs["cell_cycle_phase"].astype(str)
    ]
    adata.obs["cell_cycle_phase"] = science.pd.Categorical(
        changed,
        categories=["G1", "S", "G2M"],
    )
    history = adata.uns["openbio_singlecell"]["analysis_history"]["000000"]
    history["parameters"]["score_bundle_fingerprint_sha256"] = _cell_cycle_score_fingerprint(
        adata.obs_names,
        adata.obs["cell_cycle_s_score"],
        adata.obs["cell_cycle_g2m_score"],
        adata.obs["cell_cycle_phase"].astype(str),
    )

    with pytest.raises(ValueError, match="phase labels disagree"):
        cell_cycle_score_plot_owned(adata)


def test_cell_cycle_score_plot_matches_scanpy_nonnegative_tie_as_s_phase(science):
    adata = _scored_adata(science)
    adata.obs.loc["cell_8", "cell_cycle_s_score"] = 0.5
    adata.obs.loc["cell_8", "cell_cycle_g2m_score"] = 0.5
    phases = adata.obs["cell_cycle_phase"].astype(str).tolist()
    phases[-1] = "S"
    adata.obs["cell_cycle_phase"] = science.pd.Categorical(phases, categories=["G1", "S", "G2M"])
    history = adata.uns["openbio_singlecell"]["analysis_history"]["000000"]
    history["parameters"]["score_bundle_fingerprint_sha256"] = _cell_cycle_score_fingerprint(
        adata.obs_names,
        adata.obs["cell_cycle_s_score"],
        adata.obs["cell_cycle_g2m_score"],
        phases,
    )

    _, report, _ = cell_cycle_score_plot_owned(adata)

    assert report.summary["key_results"]["phase_counts"] == {"G1": 3, "S": 4, "G2M": 2}
    assert "ties resolve to S" in report.summary["key_results"]["phase_rule"]


def test_cell_cycle_score_plot_worker_preserves_the_input_h5ad(tmp_path, science):
    adata = _scored_adata(science)
    input_root = tmp_path / "cell-cycle"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "cell-cycle-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = cell_cycle_score_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"output_prefix": "cell_cycle"},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"output_prefix": "cell_cycle"}
    assert hashlib.sha256(input_path.read_bytes()).digest() == before


def test_cell_cycle_score_plot_restores_matplotlib_global_state(science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    matplotlib.rcParams["lines.linewidth"] = 3.25
    before_rc = matplotlib.rcParams["lines.linewidth"]
    before_figures = pyplot.get_fignums()

    cell_cycle_score_plot_owned(_scored_adata(science))

    assert matplotlib.rcParams["lines.linewidth"] == before_rc
    assert pyplot.get_fignums() == before_figures
