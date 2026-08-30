from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .contracts import PlotResult, SingleCellResult, SummaryResult, TableResult

_BASE_FIELDS = (
    "title",
    "parameters",
    "description",
    "warnings",
    "input_cells",
    "input_genes",
    "random_seed",
    "elapsed_seconds",
    "source",
)


def result_metadata(result: SingleCellResult) -> dict[str, Any]:
    if not isinstance(result, (SummaryResult, TableResult, PlotResult)):
        raise TypeError("Expected an OpenBio single-cell result value.")
    metadata = {"kind": result.kind}
    metadata.update({name: getattr(result, name) for name in _BASE_FIELDS})
    if isinstance(result, SummaryResult):
        metadata["summary"] = result.summary
    return metadata


def _base_metadata(metadata: Mapping[str, Any], *, kind: str) -> dict[str, Any]:
    if not isinstance(metadata, Mapping) or metadata.get("kind") != kind:
        raise ValueError(f"Expected {kind!r} result metadata.")
    missing = [name for name in _BASE_FIELDS if name not in metadata]
    if missing:
        raise ValueError(f"Result metadata is missing fields: {missing!r}.")
    return {name: metadata[name] for name in _BASE_FIELDS}


def summary_from_metadata(metadata: Mapping[str, Any]) -> SummaryResult:
    fields = _base_metadata(metadata, kind="summary")
    if "summary" not in metadata:
        raise ValueError("Summary result metadata is missing its summary payload.")
    return SummaryResult(summary=metadata["summary"], **fields)


def table_from_metadata(metadata: Mapping[str, Any], table: Any) -> TableResult:
    return TableResult(table=table, **_base_metadata(metadata, kind="table"))


def plot_from_metadata(metadata: Mapping[str, Any], png: bytes) -> PlotResult:
    return PlotResult(png=png, **_base_metadata(metadata, kind="plot"))


def write_strict_json(path: str | Path, value: Any) -> None:
    target = Path(path)
    payload = json.dumps(value, allow_nan=False, ensure_ascii=False, separators=(",", ":")) + "\n"
    with target.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(payload)


__all__ = [
    "plot_from_metadata",
    "result_metadata",
    "summary_from_metadata",
    "table_from_metadata",
    "write_strict_json",
]
