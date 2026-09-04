from __future__ import annotations

import hashlib
import json
import time
import uuid

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import read_plot, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.nodes_annotation import OpenBioSingleCellMarkerORAEvidencePlot
from openbio_singlecell.operations_annotation import marker_ora_evidence_plot, marker_ora_evidence_plot_owned
from openbio_singlecell.ora_evidence import ORA_EVIDENCE_COLUMNS, marker_ora_table_fingerprint
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _log_odds(science, a, b, c, d):
    return float(science.np.log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))))


def _marker_ora_result(science):
    rows = [
        (
            "A",
            "Type1",
            3,
            10,
            5,
            4,
            2,
            '["G1","G2"]',
            2,
            2,
            1,
            5,
            _log_odds(science, 2, 2, 1, 5),
            0.01,
            0.02,
            0.04,
            1,
            True,
            True,
            True,
            "eligible",
        ),
        (
            "A",
            "Type2",
            3,
            10,
            4,
            3,
            1,
            '["G3"]',
            1,
            2,
            2,
            5,
            _log_odds(science, 1, 2, 2, 5),
            0.30,
            0.30,
            0.40,
            2,
            True,
            False,
            False,
            "insufficient_overlap",
        ),
        (
            "B",
            "Type2",
            2,
            10,
            4,
            3,
            2,
            '["G4","G5"]',
            2,
            1,
            0,
            7,
            _log_odds(science, 2, 1, 0, 7),
            0.02,
            0.04,
            0.08,
            1,
            True,
            True,
            True,
            "eligible",
        ),
        (
            "B",
            "Type1",
            2,
            10,
            5,
            4,
            0,
            "[]",
            0,
            4,
            2,
            4,
            _log_odds(science, 0, 4, 2, 4),
            0.50,
            0.50,
            0.50,
            2,
            False,
            False,
            False,
            "nonpositive",
        ),
    ]
    table = science.pd.DataFrame(rows, columns=ORA_EVIDENCE_COLUMNS)
    parameters = {
        "producer_node_id": "OpenBioSingleCellMarkerORAEvidence",
        "artifact_role": "marker_ora_evidence_table",
        "analysis_fingerprint": "a" * 64,
        "ranking_fingerprint": "b" * 64,
        "universe_fingerprint": "c" * 64,
        "marker_content_fingerprint": "d" * 64,
        "universe_content_fingerprint": "e" * 64,
        "upstream_marker_content_fingerprint": "f" * 64,
        "resource_csv": "marker_sets.csv",
        "resource_sha256": "1" * 64,
        "resource_metadata": {
            "name": "Marker sets",
            "version": "1.0",
            "organism": "human",
            "identifier_namespace": "gene_symbol",
            "citation": "Test marker sets",
            "license": "CC0",
            "retrieved_at": "2026-01-01",
            "url": "https://example.org/marker-sets",
        },
        "source_column": "source",
        "target_column": "target",
        "min_targets": 2,
        "min_overlap": 2,
        "max_p_adjusted": 0.05,
        "alternative": "greater",
        "n_bg": 10,
        "ha_corr": 0.5,
        "within_group_correction": "benjamini-hochberg",
        "global_correction": "benjamini-hochberg",
    }
    parameters["table_content_fingerprint_sha256"] = marker_ora_table_fingerprint(table)
    return make_table_result(
        table=table,
        title="Marker ORA evidence",
        operation="marker_ora_evidence",
        parameters=parameters,
        description="Exploratory Marker ORA evidence.",
        warnings=[],
        input_cells=40,
        input_genes=10,
        started_at=time.perf_counter(),
    )


def test_marker_ora_evidence_plot_schema_has_two_closed_annotation_views():
    schema = OpenBioSingleCellMarkerORAEvidencePlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellMarkerORAEvidencePlot"
    assert schema.display_name == "Marker ORA Evidence Plot"
    assert schema.category == "openbio/single-cell/annotation"
    assert [item.id for item in schema.inputs] == ["table", "view"]
    assert [option.key for option in schema.inputs[1].options] == ["enrichment_dot", "overlap_bar"]
    assert [[item.id for item in option.inputs] for option in schema.inputs[1].options] == [
        ["max_terms_per_group"],
        ["max_terms_per_group"],
    ]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


@pytest.mark.parametrize("view", ["enrichment_dot", "overlap_bar"])
def test_marker_ora_plot_reads_canonical_evidence_without_annotation_or_reanalysis(science, view):
    result = _marker_ora_result(science)
    before = result.table.copy(deep=True)

    plotted, report, code = marker_ora_evidence_plot_owned(
        result,
        view={"view": view, "max_terms_per_group": 2},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["view"] == view
    assert details["groups"] == ["A", "B"]
    assert details["plotted_rows"] == 4
    assert details["plotted_sources"] == ["Type1", "Type2"]
    assert details["eligible_rows"] == 2
    assert details["displayed_rows"][0]["overlap_genes"] == ["G1", "G2"]
    text = " ".join([report.summary["results"], *report.summary["limitations"]])
    assert "Cluster marker evidence" in text
    assert "Curated annotation" in text
    assert "Condition contrast" in text
    assert "resource" not in report.summary["methods"].lower().replace("stored resource", "")
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<marker-ora-plot-code>", "exec"), namespace)
    assert namespace["plot_marker_ora_evidence"](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, before)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda result: result.source.__setitem__("operation", "gene_set_overrepresentation"), "producer"),
        (lambda result: result.table.loc.__setitem__((0, "a"), 1), "contingency"),
        (lambda result: result.table.loc.__setitem__((0, "log_odds_ratio"), 99.0), "log odds"),
        (lambda result: result.table.loc.__setitem__((0, "candidate_status"), "curated"), "candidate status"),
    ],
)
def test_marker_ora_plot_rejects_wrong_or_tampered_evidence(science, mutation, message):
    result = _marker_ora_result(science)
    mutation(result)

    with pytest.raises((TypeError, ValueError), match=message):
        marker_ora_evidence_plot_owned(result)


def test_marker_ora_plot_worker_preserves_table_artifact(tmp_path, science):
    result = _marker_ora_result(science)
    input_root = tmp_path / "marker-ora"
    input_root.mkdir()
    write_table(input_root, result.table, result_metadata(result))
    before = {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}
    staging = tmp_path / "marker-ora-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = marker_ora_evidence_plot(
        context,
        {
            "table": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"view": {"view": "enrichment_dot", "max_terms_per_group": 2}},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"]["view"]["view"] == "enrichment_dot"
    assert before == {
        path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()
    }


def test_marker_ora_plot_rejects_inactive_view_fields(science):
    with pytest.raises(ValueError, match="inactive or unknown"):
        marker_ora_evidence_plot_owned(
            _marker_ora_result(science),
            view={"view": "enrichment_dot", "max_terms_per_group": 2, "other": True},
        )


def test_marker_ora_plot_restores_matplotlib_global_state(science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    matplotlib.rcParams["lines.linewidth"] = 3.25
    before_rc = matplotlib.rcParams["lines.linewidth"]
    before_figures = pyplot.get_fignums()

    marker_ora_evidence_plot_owned(_marker_ora_result(science))

    assert matplotlib.rcParams["lines.linewidth"] == before_rc
    assert pyplot.get_fignums() == before_figures
