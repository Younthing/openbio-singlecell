from __future__ import annotations

import inspect
import time
from dataclasses import FrozenInstanceError

import pytest

from openbio_singlecell.analysis_utils import make_plot_result, make_summary_result, make_table_result
from openbio_singlecell.contracts import PlotResult, SummaryResult, TableResult
from openbio_singlecell.payload import MAX_COLUMNS, MAX_ROWS, result_to_payload

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def canonical_summary(**extra):
    return {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellTest",
        "methods": "Test method.",
        "results": "Test result.",
        "key_results": {},
        "parameters": {},
        "warnings": [],
        "limitations": [],
        "references": [{"citation": "Test citation.", "url": "https://example.org", "kind": "method"}],
        "software_versions": {"python": "test", "openbio-singlecell": "test"},
        **extra,
    }


def result_fields():
    return dict(
        title="result",
        parameters={},
        description="description",
        warnings=[],
        input_cells=1,
        input_genes=1,
        random_seed=0,
        elapsed_seconds=0.25,
        source={},
    )


def factory_fields():
    return dict(
        title="result",
        operation="test_result",
        parameters={},
        description="description",
        warnings=[],
        input_cells=1,
        input_genes=1,
        started_at=time.perf_counter(),
    )


def test_table_payload_is_bounded_and_normalizes_nonfinite_values(science):
    rows = [[row * 1000 + column for column in range(70)] for row in range(105)]
    rows[0][0] = float("nan")
    rows[0][1] = float("inf")
    table = science.pd.DataFrame(rows, columns=[f"column_{index}" for index in range(70)])

    payload = result_to_payload(TableResult(table=table, **result_fields()))

    assert len(payload["columns"]) == MAX_COLUMNS
    assert len(payload["rows"]) == MAX_ROWS
    assert payload["total_rows"] == 105
    assert payload["rows"][0][:2] == [None, None]
    assert any("Export CSV" in warning for warning in payload["warnings"])


def test_summary_collections_and_strings_are_bounded():
    summary = canonical_summary(
        values=list(range(100)),
        long="x" * 5000,
        **{f"key_{index}": index for index in range(100)},
    )
    payload = result_to_payload(SummaryResult(summary=summary, **result_fields()))

    assert len(payload["summary"]["values"]) == 65
    assert len(payload["summary"]["long"]) <= 2000
    assert payload["summary"]["…"] == "truncated"


def test_specific_factories_construct_the_concrete_result_types(science):
    summary = make_summary_result(summary=canonical_summary(cells=1), **factory_fields())
    table = make_table_result(table=science.pd.DataFrame({"value": [1]}), **factory_fields())
    plot = make_plot_result(png=PNG_SIGNATURE + b"plot", **factory_fields())

    assert isinstance(summary, SummaryResult) and summary.kind == "summary"
    assert isinstance(table, TableResult) and table.kind == "table"
    assert isinstance(plot, PlotResult) and plot.kind == "plot"


def test_specific_factory_interfaces_cannot_accept_the_wrong_payload(science):
    signatures = {
        "summary": inspect.signature(make_summary_result).parameters,
        "table": inspect.signature(make_table_result).parameters,
        "plot": inspect.signature(make_plot_result).parameters,
    }

    assert (
        "summary" in signatures["summary"]
        and "table" not in signatures["summary"]
        and "png" not in signatures["summary"]
    )
    assert "table" in signatures["table"] and "summary" not in signatures["table"] and "png" not in signatures["table"]
    assert "png" in signatures["plot"] and "summary" not in signatures["plot"] and "table" not in signatures["plot"]
    assert all("kind" not in parameters for parameters in signatures.values())

    table = science.pd.DataFrame({"value": [1]})
    with pytest.raises(TypeError):
        make_summary_result(table=table, **factory_fields())
    with pytest.raises(TypeError):
        make_table_result(summary={"cells": 1}, **factory_fields())
    with pytest.raises(TypeError):
        make_plot_result(table=table, **factory_fields())


def test_concrete_results_expose_only_their_payload(science):
    summary = SummaryResult(summary=canonical_summary(cells=1), **result_fields())
    table = TableResult(table=science.pd.DataFrame({"value": [1]}), **result_fields())
    plot = PlotResult(png=PNG_SIGNATURE + b"plot", **result_fields())

    assert not hasattr(summary, "table") and not hasattr(summary, "png")
    assert not hasattr(table, "summary") and not hasattr(table, "png")
    assert not hasattr(plot, "summary") and not hasattr(plot, "table")
    with pytest.raises(TypeError, match="unexpected keyword argument 'table'"):
        SummaryResult(summary=canonical_summary(cells=1), table=table.table, **result_fields())


def test_concrete_results_reject_invalid_content():
    with pytest.raises(TypeError, match="must be a mapping"):
        SummaryResult(summary=None, **result_fields())
    with pytest.raises(TypeError, match="pandas DataFrame"):
        TableResult(table=object(), **result_fields())
    with pytest.raises(ValueError, match="standard PNG signature"):
        PlotResult(png=b"png", **result_fields())
    with pytest.raises(TypeError, match="must be bytes"):
        PlotResult(png="png", **result_fields())


def test_result_kind_and_fields_are_read_only(science):
    summary = SummaryResult(summary=canonical_summary(cells=1), **result_fields())
    table = TableResult(table=science.pd.DataFrame({"value": [1]}), **result_fields())
    plot = PlotResult(png=PNG_SIGNATURE + b"plot", **result_fields())

    for result in (summary, table, plot):
        with pytest.raises(FrozenInstanceError):
            result.kind = "changed"
        with pytest.raises(FrozenInstanceError):
            result.title = "changed"
