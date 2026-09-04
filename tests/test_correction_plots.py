from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.nodes_correction import (
    OpenBioSingleCellMADOutlierPlot,
    OpenBioSingleCellScrubletDiagnosticsPlot,
)
from openbio_singlecell.operations_correction import (
    mad_outlier_plot,
    mad_outlier_plot_owned,
    mark_mad_outliers_owned,
    scrublet_diagnostics_plot,
    scrublet_diagnostics_plot_owned,
    scrublet_owned,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _mad_adata(science):
    return science.ad.AnnData(
        science.np.ones((6, 2), dtype=float),
        obs=science.pd.DataFrame(
            {
                "sample": ["s1"] * 3 + ["s2"] * 3,
                "total_counts": [1.0, 2.0, 9.0, 10.0, 11.0, 18.0],
            },
            index=[f"cell_{index}" for index in range(6)],
        ),
    )


def _scrublet_adata(science):
    return science.ad.AnnData(
        science.np.asarray(
            [[2, 0, 1, 0], [3, 1, 0, 1], [2, 1, 1, 0], [0, 3, 1, 1], [1, 2, 0, 2], [0, 2, 2, 1]],
            dtype=float,
        ),
        obs=science.pd.DataFrame(
            {"sample": ["s1"] * 3 + ["s2"] * 3},
            index=[f"cell_{index}" for index in range(6)],
        ),
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(4)]),
    )


def _install_scrublet_result(monkeypatch, science):
    def fake_scrublet(work, **parameters):
        scores = science.np.asarray([0.1, 0.4, 0.2, 0.3, 0.8, 0.1], dtype=float)
        work.obs["doublet_score"] = scores
        work.obs["predicted_doublet"] = scores > 0.25
        work.uns["scrublet"] = {
            "batches": {
                "s1": {
                    "threshold": 0.25,
                    "doublet_scores_sim": science.np.asarray([0.2, 0.3, 0.4, 0.5]),
                    "parameters": {},
                },
                "s2": {
                    "threshold": 0.25,
                    "doublet_scores_sim": science.np.asarray([0.25, 0.35, 0.45, 0.55]),
                    "parameters": {},
                },
            },
            "batched_by": parameters["batch_key"],
        }

    monkeypatch.setattr(science.sc.pp, "scrublet", fake_scrublet)


def test_mad_outlier_plot_schema_is_domain_specific_and_closed():
    schema = OpenBioSingleCellMADOutlierPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellMADOutlierPlot"
    assert schema.display_name == "MAD Outlier Plot"
    assert schema.category == "openbio/single-cell/correction"
    assert [item.id for item in schema.inputs] == ["adata", "view"]
    assert [option.key for option in schema.inputs[1].options] == [
        "metric_distributions",
        "sample_marked_fraction",
    ]
    assert [option.inputs for option in schema.inputs[1].options] == [[], []]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_scrublet_diagnostics_plot_schema_is_domain_specific_and_closed():
    schema = OpenBioSingleCellScrubletDiagnosticsPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellScrubletDiagnosticsPlot"
    assert schema.display_name == "Scrublet Diagnostics Plot"
    assert schema.category == "openbio/single-cell/correction"
    assert [item.id for item in schema.inputs] == ["adata", "view"]
    assert [option.key for option in schema.inputs[1].options] == [
        "score_distributions",
        "sample_predicted_fraction",
    ]
    assert [option.inputs for option in schema.inputs[1].options] == [[], []]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_mad_outlier_producer_retains_portable_threshold_and_axis_evidence(tmp_path, science):
    adata = _mad_adata(science)

    output, _, _ = mark_mad_outliers_owned(
        adata,
        metrics="total_counts",
        batch_key="sample",
        nmads=3.0,
        direction="upper",
    )
    evidence = output.uns["openbio_mad_outliers"]

    assert set(evidence) == {
        "schema_version",
        "producer_operation",
        "observation_axis_fingerprint_sha256",
        "metrics",
        "batch_key",
        "output_column",
        "nmads",
        "direction",
        "scale_mad",
        "minimum_group_size",
        "thresholds",
        "group_rates",
    }
    assert evidence["producer_operation"] == "mark_mad_outliers"
    assert evidence["observation_axis_fingerprint_sha256"].startswith("sha256:")
    assert evidence["thresholds"].columns.tolist() == [
        "metric",
        "group",
        "status",
        "n",
        "finite",
        "missing",
        "median",
        "mad",
        "lower_threshold",
        "upper_threshold",
        "flagged",
    ]
    assert evidence["thresholds"][["metric", "group", "flagged"]].to_dict("records") == [
        {"metric": "total_counts", "group": "s1", "flagged": 1},
        {"metric": "total_counts", "group": "s2", "flagged": 1},
    ]
    assert evidence["group_rates"].to_dict("records") == [
        {"group": "s1", "cells": 3, "marked": 1, "marked_fraction": 1 / 3},
        {"group": "s2", "cells": 3, "marked": 1, "marked_fraction": 1 / 3},
    ]

    path = tmp_path / "mad.h5ad"
    output.write_h5ad(path)
    restored = science.ad.read_h5ad(path)
    restored_evidence = restored.uns["openbio_mad_outliers"]
    assert restored_evidence["observation_axis_fingerprint_sha256"] == evidence[
        "observation_axis_fingerprint_sha256"
    ]
    science.pd.testing.assert_frame_equal(restored_evidence["thresholds"], evidence["thresholds"])
    science.pd.testing.assert_frame_equal(restored_evidence["group_rates"], evidence["group_rates"])


def test_mad_metric_distribution_plot_reads_stored_evidence_and_code_matches(science):
    output, _, _ = mark_mad_outliers_owned(
        _mad_adata(science),
        metrics="total_counts",
        batch_key="sample",
        nmads=3.0,
        direction="upper",
    )
    snapshot = output.copy()

    plotted, report, code = mad_outlier_plot_owned(
        output,
        view={"view": "metric_distributions"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellMADOutlierPlot"
    assert summary["parameters"] == {"view": {"view": "metric_distributions"}}
    assert summary["key_results"]["view"] == "metric_distributions"
    assert summary["key_results"]["metrics"] == ["total_counts"]
    assert summary["key_results"]["groups"] == ["s1", "s2"]
    assert summary["key_results"]["plotted_observations"] == 6
    assert summary["key_results"]["plotted_thresholds"] == 2
    assert summary["key_results"]["threshold_line_legend"] == "Stored MAD threshold"
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<mad-outlier-plot-code>", "exec"), namespace)
    assert namespace["plot_mad_outliers"](output) == plotted.png

    science.np.testing.assert_array_equal(output.X, snapshot.X)
    science.pd.testing.assert_frame_equal(output.obs, snapshot.obs)
    science.pd.testing.assert_frame_equal(
        output.uns["openbio_mad_outliers"]["thresholds"],
        snapshot.uns["openbio_mad_outliers"]["thresholds"],
    )


def test_scrublet_producer_retains_portable_simulated_score_and_axis_evidence(
    tmp_path, monkeypatch, science
):
    _install_scrublet_result(monkeypatch, science)

    output, _, _ = scrublet_owned(
        _scrublet_adata(science),
        batch_key="sample",
        threshold_mode="manual",
        threshold=0.25,
        n_prin_comps=2,
    )
    evidence = output.uns["openbio_scrublet_diagnostics"]

    assert set(evidence) == {
        "schema_version",
        "producer_operation",
        "observation_axis_fingerprint_sha256",
        "batch_key",
        "score_column",
        "prediction_column",
        "group_statistics",
        "score_evidence_fingerprint_sha256",
    }
    assert evidence["producer_operation"] == "scrublet"
    assert evidence["observation_axis_fingerprint_sha256"].startswith("sha256:")
    assert evidence["score_evidence_fingerprint_sha256"].startswith("sha256:")
    assert evidence["group_statistics"].columns.tolist() == [
        "group",
        "cells",
        "predicted_doublets",
        "predicted_fraction",
        "threshold",
        "simulated_count",
    ]
    assert evidence["group_statistics"].to_dict("records") == [
        {
            "group": "s1",
            "cells": 3,
            "predicted_doublets": 1,
            "predicted_fraction": 1 / 3,
            "threshold": 0.25,
            "simulated_count": 4,
        },
        {
            "group": "s2",
            "cells": 3,
            "predicted_doublets": 2,
            "predicted_fraction": 2 / 3,
            "threshold": 0.25,
            "simulated_count": 4,
        },
    ]

    path = tmp_path / "scrublet.h5ad"
    output.write_h5ad(path)
    restored = science.ad.read_h5ad(path)
    science.np.testing.assert_array_equal(
        restored.uns["scrublet"]["batches"]["s1"]["doublet_scores_sim"],
        [0.2, 0.3, 0.4, 0.5],
    )
    science.pd.testing.assert_frame_equal(
        restored.uns["openbio_scrublet_diagnostics"]["group_statistics"],
        evidence["group_statistics"],
    )


def test_scrublet_score_distribution_plot_uses_observed_simulated_threshold_evidence_and_code_matches(
    monkeypatch, science
):
    _install_scrublet_result(monkeypatch, science)
    output, _, _ = scrublet_owned(
        _scrublet_adata(science),
        batch_key="sample",
        threshold_mode="manual",
        threshold=0.25,
        n_prin_comps=2,
    )
    snapshot = output.copy()

    plotted, report, code = scrublet_diagnostics_plot_owned(
        output,
        view={"view": "score_distributions"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellScrubletDiagnosticsPlot"
    assert summary["parameters"] == {"view": {"view": "score_distributions"}}
    assert summary["key_results"]["view"] == "score_distributions"
    assert summary["key_results"]["groups"] == ["s1", "s2"]
    assert summary["key_results"]["plotted_observed_scores"] == 6
    assert summary["key_results"]["plotted_simulated_scores"] == 8
    assert summary["key_results"]["plotted_thresholds"] == 2
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<scrublet-diagnostics-plot-code>", "exec"), namespace)
    assert namespace["plot_scrublet_diagnostics"](output) == plotted.png

    science.np.testing.assert_array_equal(output.X, snapshot.X)
    science.pd.testing.assert_frame_equal(output.obs, snapshot.obs)
    science.pd.testing.assert_frame_equal(
        output.uns["openbio_scrublet_diagnostics"]["group_statistics"],
        snapshot.uns["openbio_scrublet_diagnostics"]["group_statistics"],
    )


def test_correction_plot_alternate_views_restore_matplotlib_state(monkeypatch, science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    mad, _, _ = mark_mad_outliers_owned(
        _mad_adata(science), metrics="total_counts", batch_key="sample", nmads=3.0, direction="upper"
    )
    _install_scrublet_result(monkeypatch, science)
    scrublet_output, _, _ = scrublet_owned(
        _scrublet_adata(science), batch_key="sample", threshold_mode="manual", threshold=0.25, n_prin_comps=2
    )
    original = matplotlib.rcParams["lines.linewidth"]
    try:
        matplotlib.rcParams["lines.linewidth"] = 3.25
        rc_value = matplotlib.rcParams["lines.linewidth"]
        figure_numbers = pyplot.get_fignums()

        mad_plot, mad_report, _ = mad_outlier_plot_owned(
            mad, view={"view": "sample_marked_fraction"}
        )
        scrublet_plot, scrublet_report, _ = scrublet_diagnostics_plot_owned(
            scrublet_output, view={"view": "sample_predicted_fraction"}
        )

        assert mad_plot.png.startswith(PNG_SIGNATURE)
        assert mad_report.summary["key_results"]["group_rates"][0]["marked_fraction"] == 1 / 3
        assert scrublet_plot.png.startswith(PNG_SIGNATURE)
        assert scrublet_report.summary["key_results"]["group_statistics"][1]["predicted_fraction"] == 2 / 3
        assert matplotlib.rcParams["lines.linewidth"] == rc_value
        assert pyplot.get_fignums() == figure_numbers
    finally:
        matplotlib.rcParams["lines.linewidth"] = original


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("axis", "axis fingerprint"),
        ("schema", "exact stored Mark MAD Outliers"),
        ("history", "provenance"),
        ("threshold", "stored thresholds"),
    ],
)
def test_mad_outlier_plot_rejects_tampered_evidence_before_rendering(science, case, message):
    output, _, _ = mark_mad_outliers_owned(
        _mad_adata(science), metrics="total_counts", batch_key="sample", nmads=3.0, direction="upper"
    )
    if case == "axis":
        output = output[[1, 0, 2, 3, 4, 5]].copy()
    elif case == "schema":
        output.uns["openbio_mad_outliers"]["unexpected"] = True
    elif case == "history":
        output.uns["openbio_singlecell"]["analysis_history"] = {}
    else:
        output.uns["openbio_mad_outliers"]["thresholds"].loc[0, "upper_threshold"] += 1.0

    with pytest.raises(ValueError, match=message):
        mad_outlier_plot_owned(output)


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("axis", "axis fingerprint"),
        ("schema", "exact stored Scrublet"),
        ("history", "provenance"),
        ("simulated", "simulated scores"),
        ("simulated_content", "score evidence fingerprint"),
        ("observed_content", "score evidence fingerprint"),
    ],
)
def test_scrublet_plot_rejects_tampered_evidence_before_rendering(
    monkeypatch, science, case, message
):
    _install_scrublet_result(monkeypatch, science)
    output, _, _ = scrublet_owned(
        _scrublet_adata(science), batch_key="sample", threshold_mode="manual", threshold=0.25, n_prin_comps=2
    )
    if case == "axis":
        output = output[[1, 0, 2, 3, 4, 5]].copy()
    elif case == "schema":
        output.uns["openbio_scrublet_diagnostics"]["unexpected"] = True
    elif case == "history":
        output.uns["openbio_singlecell"]["analysis_history"] = {}
    else:
        if case == "simulated":
            output.uns["scrublet"]["batches"]["s1"]["doublet_scores_sim"][0] = science.np.nan
        elif case == "simulated_content":
            output.uns["scrublet"]["batches"]["s1"]["doublet_scores_sim"][0] = 0.21
        else:
            output.obs.loc["cell_0", "doublet_score"] = 0.11

    with pytest.raises(ValueError, match=message):
        scrublet_diagnostics_plot_owned(output)


@pytest.mark.parametrize(
    ("operation", "parameters"),
    [
        (mad_outlier_plot, {"view": {"view": "metric_distributions"}}),
        (scrublet_diagnostics_plot, {"view": {"view": "score_distributions"}}),
    ],
)
def test_correction_plot_workers_round_trip_png_without_rewriting_input(
    tmp_path, monkeypatch, science, operation, parameters
):
    if operation is mad_outlier_plot:
        output, _, _ = mark_mad_outliers_owned(
            _mad_adata(science), metrics="total_counts", batch_key="sample", nmads=3.0, direction="upper"
        )
    else:
        _install_scrublet_result(monkeypatch, science)
        output, _, _ = scrublet_owned(
            _scrublet_adata(science),
            batch_key="sample",
            threshold_mode="manual",
            threshold=0.25,
            n_prin_comps=2,
        )
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, output)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = operation(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        parameters,
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == parameters
    assert hashlib.sha256(input_path.read_bytes()).digest() == before
