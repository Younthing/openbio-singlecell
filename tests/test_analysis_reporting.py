from __future__ import annotations

import json
import time

import pytest

from openbio_singlecell.analysis_reporting import (
    AnalysisReference,
    make_analysis_report,
    summarize_numeric,
)
from openbio_singlecell.payload import result_to_payload

REFERENCE = AnalysisReference(
    citation="Example method. Journal. 2026.",
    doi="10.0000/example",
    url="https://doi.org/10.0000/example",
    kind="method",
)


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
