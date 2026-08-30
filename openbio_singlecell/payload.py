from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from datetime import date, datetime
from typing import Any, Literal, TypedDict, cast

from . import SCHEMA_VERSION
from .contracts import PlotResult, SingleCellResult, SummaryResult, TableResult

MAX_ROWS = 100
MAX_COLUMNS = 64
MAX_COLLECTION_ITEMS = 64
MAX_WARNINGS = 32
MAX_STRING_LENGTH = 2000
MAX_DEPTH = 4

type JsonScalar = None | bool | int | float | str
type JsonValue = JsonScalar | list[JsonValue] | dict[str, JsonValue]


class _ResultPayloadFields(TypedDict):
    schema_version: int
    title: str
    columns: list[str]
    rows: list[list[JsonScalar]]
    total_rows: int
    summary: JsonValue
    warnings: list[str]
    elapsed_seconds: float | None


class SummaryResultPayload(_ResultPayloadFields):
    kind: Literal["summary"]


class TableResultPayload(_ResultPayloadFields):
    kind: Literal["table"]


class PlotResultPayload(_ResultPayloadFields):
    kind: Literal["plot"]


type SingleCellResultPayload = SummaryResultPayload | TableResultPayload | PlotResultPayload


def _bounded_string(value: Any) -> str:
    text = str(value)
    if len(text) <= MAX_STRING_LENGTH:
        return text
    return text[: MAX_STRING_LENGTH - 1] + "…"


def _json_scalar(value: Any) -> JsonScalar:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, str):
        return _bounded_string(value)
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    if type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if hasattr(value, "item"):
        try:
            return _json_scalar(value.item())
        except (TypeError, ValueError):
            pass
    try:
        unequal = value != value
        if isinstance(unequal, bool) and unequal:
            return None
    except (TypeError, ValueError):
        pass
    return _bounded_string(value)


def _bounded_value(value: Any, depth: int = 0) -> JsonValue:
    if depth >= MAX_DEPTH:
        return _json_scalar(value)
    if isinstance(value, dict):
        result = {}
        for index, (key, item) in enumerate(value.items()):
            if index >= MAX_COLLECTION_ITEMS:
                result["…"] = "truncated"
                break
            result[_bounded_string(key)] = _bounded_value(item, depth + 1)
        return result
    if isinstance(value, (list, tuple)):
        result = [_bounded_value(item, depth + 1) for item in value[:MAX_COLLECTION_ITEMS]]
        if len(value) > MAX_COLLECTION_ITEMS:
            result.append("…")
        return result
    return _json_scalar(value)


def result_to_payload(result: SingleCellResult) -> SingleCellResultPayload:
    if not isinstance(result, (SummaryResult, TableResult, PlotResult)):
        raise TypeError("Expected an OpenBio single-cell result value.")

    columns: list[str] = []
    rows: list[list[JsonScalar]] = []
    total_rows = 0
    warnings = [_bounded_string(warning) for warning in result.warnings[:MAX_WARNINGS]]

    if isinstance(result, TableResult):
        table = result.table
        total_rows = int(len(table))
        columns = [_bounded_string(column) for column in list(table.columns)[:MAX_COLUMNS]]
        preview = table.iloc[:MAX_ROWS, :MAX_COLUMNS]
        rows = [[_json_scalar(value) for value in row] for row in preview.itertuples(index=False, name=None)]
        if len(table.columns) > MAX_COLUMNS or total_rows > MAX_ROWS:
            warnings.append(
                f"Preview limited to {MAX_ROWS} rows and {MAX_COLUMNS} columns; Export CSV for the complete table."
            )

    summary = result.summary if isinstance(result, SummaryResult) else {"description": result.description}

    elapsed = float(result.elapsed_seconds)
    return cast(
        SingleCellResultPayload,
        {
            "schema_version": SCHEMA_VERSION,
            "kind": result.kind,
            "title": _bounded_string(result.title),
            "columns": columns,
            "rows": rows,
            "total_rows": total_rows,
            "summary": _bounded_value(summary),
            "warnings": warnings[:MAX_WARNINGS],
            "elapsed_seconds": elapsed if math.isfinite(elapsed) else None,
        },
    )


def artifact_metadata_to_payload(
    metadata: Mapping[str, Any],
    *,
    columns: Sequence[Any] = (),
    rows: Sequence[Sequence[Any]] = (),
    total_rows: int = 0,
) -> TableResultPayload | PlotResultPayload:
    kind = metadata.get("kind")
    if kind not in {"table", "plot"}:
        raise ValueError("Artifact preview metadata must describe a table or plot result.")
    if not isinstance(metadata.get("title"), str) or not isinstance(metadata.get("description"), str):
        raise ValueError("Artifact preview metadata is missing its title or description.")
    warnings_value = metadata.get("warnings")
    if not isinstance(warnings_value, list):
        raise ValueError("Artifact preview metadata warnings must be a list.")
    bounded_columns = [_bounded_string(value) for value in list(columns)[:MAX_COLUMNS]]
    bounded_rows = [
        [_json_scalar(value) for value in list(row)[:MAX_COLUMNS]]
        for row in list(rows)[:MAX_ROWS]
    ]
    warnings = [_bounded_string(value) for value in warnings_value[:MAX_WARNINGS]]
    if kind == "table" and (len(columns) > MAX_COLUMNS or total_rows > MAX_ROWS):
        warnings.append(
            f"Preview limited to {MAX_ROWS} rows and {MAX_COLUMNS} columns; Export CSV for the complete table."
        )
    elapsed_value = metadata.get("elapsed_seconds")
    elapsed = float(elapsed_value) if isinstance(elapsed_value, (int, float)) else None
    return cast(
        TableResultPayload | PlotResultPayload,
        {
            "schema_version": SCHEMA_VERSION,
            "kind": kind,
            "title": _bounded_string(metadata["title"]),
            "columns": bounded_columns,
            "rows": bounded_rows,
            "total_rows": int(total_rows),
            "summary": _bounded_value({"description": metadata["description"]}),
            "warnings": warnings[:MAX_WARNINGS],
            "elapsed_seconds": elapsed if elapsed is not None and math.isfinite(elapsed) else None,
        },
    )
