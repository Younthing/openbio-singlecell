from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import read_plot, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.contracts import TableResult
from openbio_singlecell.leiden_sweep_plotting import leiden_resolution_metrics_fingerprint
from openbio_singlecell.nodes_integration import OpenBioSingleCellLeidenResolutionSweepPlot
from openbio_singlecell.operations_integration import (
    leiden_resolution_sweep_plot,
    leiden_resolution_sweep_plot_owned,
)
from openbio_singlecell.worker_protocol import OperationContext

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
METRIC_COLUMNS = [
    "resolution",
    "key",
    "n_clusters",
    "min_cluster_size",
    "median_cluster_size",
    "max_cluster_size",
    "smallest_cluster_fraction",
    "largest_cluster_fraction",
    "singleton_clusters",
    "modularity",
    "stability_repeats",
    "stability_mean_ari",
    "stability_min_ari",
    "stability_max_ari",
    "adjacent_previous_resolution",
    "adjacent_resolution_ari",
]


def _resolution_metrics(science):
    resolutions = [0.25, 0.5, 1.0]
    keys = ["leiden_r0p25", "leiden_r0p5", "leiden_r1"]
    rows = [
        [0.25, keys[0], 2, 30, 50.0, 70, 0.3, 0.7, 0, 0.31, 5, 0.98, 0.96, 1.0, None, None],
        [0.5, keys[1], 3, 20, 30.0, 50, 0.2, 0.5, 0, 0.44, 5, 0.95, 0.90, 0.99, 0.25, 0.81],
        [1.0, keys[2], 5, 5, 15.0, 35, 0.05, 0.35, 0, 0.51, 5, 0.91, 0.84, 0.96, 0.5, 0.72],
    ]
    table = science.pd.DataFrame(rows, columns=METRIC_COLUMNS)
    for column in ("n_clusters", "min_cluster_size", "max_cluster_size", "singleton_clusters", "stability_repeats"):
        table[column] = table[column].astype("int64")
    for column in (
        "stability_mean_ari",
        "stability_min_ari",
        "stability_max_ari",
        "adjacent_previous_resolution",
        "adjacent_resolution_ari",
    ):
        table[column] = science.pd.array(table[column], dtype="Float64")
    parameters = {
        "resolutions": resolutions,
        "keys": keys,
        "key_prefix": "leiden",
        "neighbors_key": "neighbors",
        "n_iterations": 2,
        "stability_repeats": 5,
        "random_seed": 0,
        "overwrite_existing": False,
        "total_runs": 15,
        "fixed_policy": {
            "flavor": "igraph",
            "directed": False,
            "use_weights": True,
            "objective_function": "modularity",
        },
    }
    parameters["table_content_fingerprint_sha256"] = leiden_resolution_metrics_fingerprint(table)
    source = {
        "operation": "leiden_resolution_sweep",
        "parameters": parameters,
        "input_cells": 100,
        "input_genes": 20,
        "random_seed": 0,
        "elapsed_seconds": 1.0,
        "timestamp": "2026-09-04T00:00:00+00:00",
    }
    return TableResult(
        table=table,
        title="Leiden resolution metrics",
        parameters=parameters,
        description="stored sweep",
        warnings=[],
        input_cells=100,
        input_genes=20,
        random_seed=0,
        elapsed_seconds=1.0,
        source=source,
    )


def test_leiden_resolution_sweep_plot_schema_and_quality_curves(science):
    schema = OpenBioSingleCellLeidenResolutionSweepPlot.define_schema()
    assert schema.category == "openbio/single-cell/clustering"
    assert [item.id for item in schema.inputs] == ["resolution_metrics", "view"]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]

    result = _resolution_metrics(science)
    original = result.table.copy(deep=True)
    plotted, report, code = leiden_resolution_sweep_plot_owned(
        result,
        view={"view": "quality_curves"},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["resolutions"] == [0.25, 0.5, 1.0]
    assert details["cluster_counts"] == [2, 3, 5]
    assert details["modularity"] == [0.31, 0.44, 0.51]
    assert details["view"] == "quality_curves"
    json.dumps(report.summary, allow_nan=False)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_leiden_resolution_sweep"](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, original)


def test_leiden_resolution_sweep_plot_cluster_sizes_view(science):
    plotted, report, _ = leiden_resolution_sweep_plot_owned(
        _resolution_metrics(science),
        view={"view": "cluster_sizes"},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["key_results"]["cluster_size_ranges"] == [[30, 50.0, 70], [20, 30.0, 50], [5, 15.0, 35]]


def test_leiden_resolution_sweep_plot_rejects_provenance_and_axis_tampering(science):
    result = _resolution_metrics(science)
    bad_source = {**result.source, "operation": "other"}
    wrong_producer = TableResult(
        table=result.table,
        title=result.title,
        parameters=result.parameters,
        description=result.description,
        warnings=result.warnings,
        input_cells=result.input_cells,
        input_genes=result.input_genes,
        random_seed=result.random_seed,
        elapsed_seconds=result.elapsed_seconds,
        source=bad_source,
    )
    with pytest.raises(ValueError, match="Leiden Resolution Sweep producer"):
        leiden_resolution_sweep_plot_owned(wrong_producer)

    tampered = _resolution_metrics(science)
    tampered.table.loc[1, "resolution"] = 0.75
    with pytest.raises(ValueError, match="resolution axis"):
        leiden_resolution_sweep_plot_owned(tampered)


def test_leiden_resolution_sweep_plot_rejects_valid_range_content_tampering(science):
    tampered = _resolution_metrics(science)
    _plot, _report, code = leiden_resolution_sweep_plot_owned(tampered)
    namespace: dict[str, object] = {}
    exec(code, namespace)
    tampered.table.loc[1, "modularity"] = 0.45

    with pytest.raises(ValueError, match="current-content fingerprint"):
        leiden_resolution_sweep_plot_owned(tampered)
    with pytest.raises(ValueError, match="current-content fingerprint"):
        namespace["plot_leiden_resolution_sweep"](tampered.table)


def test_leiden_resolution_sweep_plot_worker_round_trip(tmp_path, science):
    result = _resolution_metrics(science)
    input_root = tmp_path / "metrics"
    input_root.mkdir()
    write_table(input_root, result.table, result_metadata(result))
    input_files = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_root.iterdir()}
    staging = tmp_path / "plot.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = leiden_resolution_sweep_plot(
        context,
        {
            "resolution_metrics": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"view": {"view": "quality_curves"}},
    )

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert read_plot(staging / records[0]["payload"])[0].startswith(PNG_SIGNATURE)
    assert {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_root.iterdir()} == input_files
