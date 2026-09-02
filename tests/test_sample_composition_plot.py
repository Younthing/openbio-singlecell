from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import TABLE_CODEC, read_plot, write_table
from openbio_singlecell.node_types import PlotResultType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_abundance import OpenBioSingleCellSampleCompositionPlot
from openbio_singlecell.operations_abundance import sample_composition_plot, sample_composition_plot_owned
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _composition_table(science):
    table = science.pd.DataFrame(
        {
            "sample": ["s1", "s1", "s2", "s2", "s3", "s3"],
            "condition": ["control", "control", "treated", "treated", "control", "control"],
            "annotation": ["A", "B", "A", "B", "A", "B"],
            "cell_count": [3, 1, 1, 3, 0, 2],
            "sample_total_cells": [4, 4, 4, 4, 2, 2],
            "proportion": [0.75, 0.25, 0.25, 0.75, 0.0, 1.0],
        }
    )
    table["sample"] = table["sample"].astype("string")
    table["condition"] = table["condition"].astype("string")
    table["annotation"] = table["annotation"].astype("string")
    table["cell_count"] = table["cell_count"].astype("int64")
    table["sample_total_cells"] = table["sample_total_cells"].astype("int64")
    table["proportion"] = table["proportion"].astype("float64")
    return table


def _scaled_composition_table(science, *, samples, annotations):
    rows = [
        {
            "sample": f"s{sample}",
            "condition": f"condition-{sample % 2}",
            "annotation": f"annotation-{annotation}",
            "cell_count": 1,
            "sample_total_cells": annotations,
            "proportion": 1.0 / annotations,
        }
        for sample in range(samples)
        for annotation in range(annotations)
    ]
    table = science.pd.DataFrame.from_records(
        rows,
        columns=["sample", "condition", "annotation", "cell_count", "sample_total_cells", "proportion"],
    )
    table[["sample", "condition", "annotation"]] = table[["sample", "condition", "annotation"]].astype("string")
    table[["cell_count", "sample_total_cells"]] = table[["cell_count", "sample_total_cells"]].astype("int64")
    table["proportion"] = table["proportion"].astype("float64")
    return table


def test_sample_composition_plot_schema():
    schema = OpenBioSingleCellSampleCompositionPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellSampleCompositionPlot"
    assert schema.display_name == "Sample Composition Plot"
    assert schema.category == "openbio/single-cell/visualization"
    assert [item.id for item in schema.inputs] == ["table", "value"]
    assert schema.inputs[0].io_type == TableResultType.io_type
    assert schema.inputs[1].options == ["proportion", "cell_count"]
    assert schema.inputs[1].default == "proportion"
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("plot", PlotResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_sample_composition_plot_keeps_samples_separate_and_groups_them_by_condition(science):
    table = _composition_table(science)
    original = table.copy(deep=True)

    plot, report, code = sample_composition_plot_owned(table)

    assert plot.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    json.dumps(summary, allow_nan=False)
    assert summary["node_id"] == "OpenBioSingleCellSampleCompositionPlot"
    assert summary["parameters"] == {"value": "proportion"}
    assert summary["key_results"]["sample_order"] == ["s1", "s3", "s2"]
    assert summary["key_results"]["condition_order"] == ["control", "treated"]
    assert summary["key_results"]["annotation_order"] == ["A", "B"]
    assert not {
        "sample_order_preview",
        "sample_order_count",
        "sample_order_preview_truncated",
        "annotation_order_preview",
        "annotation_order_count",
        "annotation_order_preview_truncated",
    } & summary["key_results"].keys()
    assert summary["key_results"]["plotted_samples"] == 3
    assert summary["key_results"]["condition_aggregation_performed"] is False
    assert "no Condition aggregation" in summary["methods"]
    compile(code, "<sample-composition-plot-code>", "exec")
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_sample_composition"](table) == plot.png
    science.pd.testing.assert_frame_equal(table, original)


def test_sample_composition_plot_can_render_cell_counts_per_sample(science):
    table = _composition_table(science)
    proportion_plot, _, _ = sample_composition_plot_owned(table)

    count_plot, report, code = sample_composition_plot_owned(table, value="cell_count")

    assert count_plot.png.startswith(PNG_SIGNATURE)
    assert count_plot.png != proportion_plot.png
    assert report.summary["parameters"] == {"value": "cell_count"}
    assert report.summary["key_results"]["sample_order"] == ["s1", "s3", "s2"]
    assert report.summary["key_results"]["stack_totals"] == [4.0, 2.0, 4.0]
    assert report.summary["key_results"]["y_axis_label"] == "Cell count"
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_sample_composition"](table) == count_plot.png


@pytest.mark.parametrize(
    ("samples", "annotations", "message"),
    [
        (101, 1, r"at most 100 Samples.*found 101"),
        (1, 51, r"at most 50 annotation categories.*found 51"),
        (41, 49, r"at most 2,000 stacked bar segments.*found 2,009"),
    ],
)
def test_sample_composition_plot_rejects_unrenderable_scale_before_materializing_grid(
    science, monkeypatch, samples, annotations, message
):
    table = _scaled_composition_table(science, samples=samples, annotations=annotations)

    def fail_if_grid_is_materialized(*args, **kwargs):
        raise AssertionError("composition grid was materialized before the scale guard")

    monkeypatch.setattr(science.pd.DataFrame, "itertuples", fail_if_grid_is_materialized)

    with pytest.raises(ValueError, match=message):
        sample_composition_plot_owned(table)


def test_sample_composition_plot_uses_unique_colors_and_warns_at_high_annotation_cardinality(science):
    table = _scaled_composition_table(science, samples=1, annotations=50)

    plot, report, code = sample_composition_plot_owned(table)

    key_results = report.summary["key_results"]
    colors = key_results["annotation_colors"]
    assert list(colors) == [f"annotation-{index}" for index in range(50)]
    assert len(set(colors.values())) == 50
    warning = "The static PNG contains 50 annotation categories; colors and legend entries may be difficult to distinguish."
    assert key_results["annotation_readability_warnings"] == [warning]
    assert report.summary["warnings"] == [warning]
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_sample_composition"](table) == plot.png


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda table: table.drop(index=table.index[-1]).reset_index(drop=True),
            "complete Sample-by-annotation grid",
        ),
        (
            lambda table: _with_condition_conflict(table),
            "one Condition per Sample",
        ),
    ],
)
def test_sample_composition_plot_rejects_noncanonical_grid_structure(science, mutate, message):
    table = mutate(_composition_table(science))

    with pytest.raises(ValueError, match=message):
        sample_composition_plot_owned(table)


def _with_condition_conflict(table):
    conflicting = table.copy(deep=True)
    conflicting.loc[1, "condition"] = "treated"
    return conflicting


@pytest.mark.parametrize(
    ("column", "row", "value", "message"),
    [
        ("annotation", 0, " A", "canonical nonblank strings"),
        ("cell_count", 0, -1, "nonnegative integers"),
        ("cell_count", 0, 1.5, "nonnegative integers"),
        ("sample_total_cells", 1, 5, "one positive denominator per Sample"),
        ("proportion", 0, float("nan"), "finite"),
        ("proportion", 0, 0.5, "cell_count / sample_total_cells"),
    ],
)
def test_sample_composition_plot_rejects_noncanonical_row_values(science, column, row, value, message):
    table = _composition_table(science)
    if column == "cell_count" and isinstance(value, float):
        table[column] = table[column].astype(float)
    table.loc[row, column] = value

    with pytest.raises((TypeError, ValueError), match=message):
        sample_composition_plot_owned(table)


def test_sample_composition_plot_rejects_non_string_identifiers(science):
    table = _composition_table(science)
    table["sample"] = table["sample"].astype(object)
    table.loc[0, "sample"] = 1

    with pytest.raises(TypeError, match="must contain strings"):
        sample_composition_plot_owned(table)


def test_sample_composition_plot_worker_round_trips_png_without_rewriting_table(tmp_path, science):
    input_root = tmp_path / "table"
    input_root.mkdir()
    metadata = {
        "kind": "table",
        "title": "Sample composition by cell_type",
        "parameters": {},
        "description": "Canonical Sample composition.",
        "warnings": [],
        "input_cells": 10,
        "input_genes": 2,
        "random_seed": 0,
        "elapsed_seconds": 0.1,
        "source": {"operation": "sample_composition_summary"},
    }
    write_table(input_root, _composition_table(science), metadata)
    before = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_root.iterdir()}
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = sample_composition_plot(
        context,
        {
            "table": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": TABLE_CODEC,
                "path": str(input_root.resolve()),
            }
        },
        {"value": "cell_count"},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, output_metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert output_metadata["kind"] == "plot"
    assert output_metadata["parameters"] == {"value": "cell_count"}
    assert (output_metadata["input_cells"], output_metadata["input_genes"]) == (10, 2)
    assert before == {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_root.iterdir()}


def test_sample_composition_plot_worker_rejects_an_unrelated_table_artifact(tmp_path, science):
    input_root = tmp_path / "table"
    input_root.mkdir()
    write_table(
        input_root,
        _composition_table(science),
        {
            "kind": "table",
            "title": "Unrelated table",
            "parameters": {},
            "description": "Not a Sample Composition Summary output.",
            "warnings": [],
            "input_cells": 10,
            "input_genes": 2,
            "random_seed": 0,
            "elapsed_seconds": 0.1,
            "source": {"operation": "other_operation"},
        },
    )
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    with pytest.raises(ValueError, match="Sample Composition Summary table artifact"):
        sample_composition_plot(
            context,
            {
                "table": {
                    "type": "artifact",
                    "kind": "OPENBIO_SINGLE_CELL_TABLE",
                    "codec": TABLE_CODEC,
                    "path": str(input_root.resolve()),
                }
            },
            {"value": "proportion"},
        )
