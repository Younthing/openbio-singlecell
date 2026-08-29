from __future__ import annotations

import json
from io import BytesIO

import pytest
from PIL import Image

from openbio_singlecell.artifact_envelope import (
    plot_from_metadata,
    result_metadata,
    summary_from_metadata,
    table_from_metadata,
    write_strict_json,
)
from openbio_singlecell.contracts import PlotResult, SummaryResult, TableResult


def _base_fields():
    return {
        "title": "result",
        "parameters": {"alpha": 0.25},
        "description": "description",
        "warnings": ["warning"],
        "input_cells": 3,
        "input_genes": 2,
        "random_seed": 7,
        "elapsed_seconds": 0.5,
        "source": {"operation": "test"},
    }


def _summary_payload():
    return {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellTest",
        "methods": "Method.",
        "results": "Result.",
        "key_results": {},
        "parameters": {},
        "warnings": [],
        "limitations": [],
        "references": [{"citation": "Citation.", "url": "https://example.org", "kind": "method"}],
        "software_versions": {"python": "test", "openbio-singlecell": "test"},
    }


def _png_bytes():
    output = BytesIO()
    Image.new("RGB", (2, 1), "blue").save(output, format="PNG")
    return output.getvalue()


def test_summary_metadata_round_trips_as_a_small_direct_value():
    original = SummaryResult(summary=_summary_payload(), **_base_fields())

    restored = summary_from_metadata(result_metadata(original))

    assert restored == original


def test_table_and_plot_metadata_exclude_large_payloads(science):
    table = TableResult(table=science.pd.DataFrame({"value": [1, 2]}), **_base_fields())
    plot = PlotResult(png=_png_bytes(), **_base_fields())

    table_metadata = result_metadata(table)
    plot_metadata = result_metadata(plot)

    assert "table" not in table_metadata
    assert "png" not in plot_metadata
    restored_table = table_from_metadata(table_metadata, table.table)
    science.pd.testing.assert_frame_equal(restored_table.table, table.table)
    assert result_metadata(restored_table) == table_metadata
    assert plot_from_metadata(plot_metadata, plot.png) == plot


def test_strict_json_rejects_non_finite_values(tmp_path):
    target = tmp_path / "result.json"

    with pytest.raises(ValueError):
        write_strict_json(target, {"bad": float("nan")})

    assert not target.exists()


def test_result_metadata_is_plain_strict_json(tmp_path):
    target = tmp_path / "result.json"
    original = SummaryResult(summary=_summary_payload(), **_base_fields())

    write_strict_json(target, result_metadata(original))

    assert json.loads(target.read_text(encoding="utf-8"))["kind"] == "summary"
