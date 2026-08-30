from __future__ import annotations

import importlib.metadata
import json
import time
from copy import deepcopy

import pytest

from openbio_singlecell import PLUGIN_VERSION
from openbio_singlecell.analysis_reporting import (
    AnalysisReference,
    collect_software_versions,
    make_analysis_report,
    summarize_numeric,
)
from openbio_singlecell.contracts import SummaryResult
from openbio_singlecell.payload import result_to_payload

REFERENCE = AnalysisReference(
    citation="Example method. Journal. 2026.",
    doi="10.0000/example",
    url="https://doi.org/10.0000/example",
    kind="method",
)


def test_collect_software_versions_uses_the_plugin_contract_version(monkeypatch):
    monkeypatch.setattr(importlib.metadata, "version", lambda _package: "distribution-version")

    assert collect_software_versions(())["openbio-singlecell"] == PLUGIN_VERSION


def _summary_result_fields() -> dict[str, object]:
    return {
        "title": "Example summary",
        "parameters": {},
        "description": "Example result.",
        "warnings": [],
        "input_cells": 1,
        "input_genes": 1,
        "random_seed": 0,
        "elapsed_seconds": 0.0,
        "source": {},
    }


def _canonical_summary() -> dict[str, object]:
    return {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellExample",
        "methods": "The example method was applied.",
        "results": "One result was produced.",
        "key_results": {},
        "parameters": {},
        "warnings": [],
        "limitations": [],
        "references": [
            {
                "citation": "Example method. Journal. 2026.",
                "url": "https://doi.org/10.0000/example",
                "kind": "method",
            }
        ],
        "software_versions": {"python": "3.13.0", "openbio-singlecell": "0.2.0"},
    }


def test_summary_result_rejects_missing_canonical_report_field():
    summary = _canonical_summary()
    del summary["methods"]

    with pytest.raises(ValueError, match="missing required fields.*methods"):
        SummaryResult(summary=summary, **_summary_result_fields())


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("schema_version", True, "schema_version must be a positive integer or non-empty string"),
        ("node_id", " ", "node_id must be a non-empty string"),
        ("methods", ["method"], "methods must be a non-empty string"),
        ("results", "", "results must be a non-empty string"),
        ("key_results", [], "key_results must be a mapping"),
        ("parameters", [], "parameters must be a mapping"),
        ("warnings", "warning", "warnings must be a list of strings"),
        ("limitations", [1], "limitations must be a list of strings"),
        ("software_versions", [], "software_versions must be a mapping"),
    ],
)
def test_summary_result_rejects_invalid_canonical_field_types(field, value, message):
    summary = deepcopy(_canonical_summary())
    summary[field] = value

    with pytest.raises((TypeError, ValueError), match=message):
        SummaryResult(summary=summary, **_summary_result_fields())


@pytest.mark.parametrize(
    ("references", "message"),
    [
        ([], "requires at least one reference"),
        (["citation"], "reference 0 must be a mapping"),
        ([{"citation": "Citation", "url": "https://example.org"}], "reference 0 is missing.*kind"),
        (
            [{"citation": " ", "url": "https://example.org", "kind": "method"}],
            "reference 0 citation must be a non-empty string",
        ),
    ],
)
def test_summary_result_rejects_invalid_references(references, message):
    summary = deepcopy(_canonical_summary())
    summary["references"] = references

    with pytest.raises((TypeError, ValueError), match=message):
        SummaryResult(summary=summary, **_summary_result_fields())


@pytest.mark.parametrize(
    ("software_versions", "message"),
    [
        ({"openbio-singlecell": "0.2.0"}, "missing required software versions.*python"),
        ({"python": "3.13.0"}, "missing required software versions.*openbio-singlecell"),
        (
            {"python": "", "openbio-singlecell": "0.2.0"},
            "software version python must be a non-empty string",
        ),
        (
            {"python": "3.13.0", "openbio-singlecell": 2},
            "software version openbio-singlecell must be a non-empty string",
        ),
    ],
)
def test_summary_result_requires_runtime_and_openbio_versions(software_versions, message):
    summary = deepcopy(_canonical_summary())
    summary["software_versions"] = software_versions

    with pytest.raises((TypeError, ValueError), match=message):
        SummaryResult(summary=summary, **_summary_result_fields())


def test_summary_result_requires_strict_json_values():
    summary = deepcopy(_canonical_summary())
    summary["key_results"] = {"not_finite": float("nan")}

    with pytest.raises(ValueError, match="strict JSON-compatible"):
        SummaryResult(summary=summary, **_summary_result_fields())


def test_analysis_report_is_strict_json_and_records_versions():
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellExample",
        title="Example report",
        operation="example_report",
        methods="The example method was applied.",
        results="One result was produced.",
        key_results={"finite": 1.0, "not_finite": float("nan")},
        parameters={"threshold": 2},
        references=(REFERENCE,),
        software_packages=("anndata",),
        warnings=("Example warning.",),
        limitations=("Example limitation.",),
        input_cells=3,
        input_genes=4,
        started_at=time.perf_counter(),
        code="def example(adata):\n    return adata",
    )

    assert report.summary["key_results"] == {"finite": 1.0, "not_finite": None}
    assert report.summary["software_versions"]["openbio-singlecell"] == "0.2.0"
    assert report.summary["software_versions"]["anndata"] != "not-installed"
    assert json.loads(json.dumps(report.summary, allow_nan=False)) == report.summary
    json.dumps(result_to_payload(report), allow_nan=False)
    assert code == "def example(adata):\n    return adata\n"


def test_analysis_report_rejects_missing_references_and_code():
    fields = {
        "node_id": "OpenBioSingleCellExample",
        "title": "Example report",
        "operation": "example_report",
        "methods": "Methods.",
        "results": "Results.",
        "key_results": {},
        "parameters": {},
        "software_packages": (),
        "warnings": (),
        "limitations": (),
        "input_cells": 1,
        "input_genes": 1,
        "started_at": time.perf_counter(),
    }
    with pytest.raises(ValueError, match="requires at least one"):
        make_analysis_report(references=(), code="pass", **fields)
    with pytest.raises(ValueError, match="code cannot be empty"):
        make_analysis_report(references=(REFERENCE,), code=" ", **fields)


def test_numeric_summary_reports_missing_values(science):
    summary = summarize_numeric(science.np.asarray([1.0, science.np.nan, 3.0]))

    assert summary == {
        "n": 3,
        "missing": 1,
        "min": 1.0,
        "q1": 1.5,
        "median": 2.0,
        "mean": 2.0,
        "q3": 2.5,
        "max": 3.0,
    }


def test_analysis_report_normalizes_numpy_and_pandas_collections(science):
    report, _ = make_analysis_report(
        node_id="OpenBioSingleCellExample",
        title="Collection report",
        operation="collection_report",
        methods="Methods.",
        results="Results.",
        key_results={
            "array": science.np.asarray([1.0, science.np.nan]),
            "series": science.pd.Series([2, 3]),
        },
        parameters={},
        references=(REFERENCE,),
        software_packages=(),
        warnings=(),
        limitations=(),
        input_cells=2,
        input_genes=1,
        started_at=time.perf_counter(),
        code="def example(adata):\n    return adata",
    )

    assert report.summary["key_results"] == {"array": [1.0, None], "series": [2, 3]}
    json.dumps(report.summary, allow_nan=False)
