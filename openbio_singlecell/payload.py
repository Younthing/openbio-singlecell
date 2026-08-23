from __future__ import annotations

import math
from datetime import date, datetime
from typing import Any, Literal, TypedDict, cast

from . import SCHEMA_VERSION
from .contracts import SingleCellResult

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
    columns: list[str] = []
    rows: list[list[JsonScalar]] = []
    total_rows = 0
    warnings = [_bounded_string(warning) for warning in result.warnings[:MAX_WARNINGS]]

    if result.kind == "table" and result.table is not None:
        table = result.table
        total_rows = int(len(table))
        columns = [_bounded_string(column) for column in list(table.columns)[:MAX_COLUMNS]]
        preview = table.iloc[:MAX_ROWS, :MAX_COLUMNS]
        rows = [[_json_scalar(value) for value in row] for row in preview.itertuples(index=False, name=None)]
        if len(table.columns) > MAX_COLUMNS or total_rows > MAX_ROWS:
            warnings.append(
                f"Preview limited to {MAX_ROWS} rows and {MAX_COLUMNS} columns; Export CSV for the complete table."
            )

    summary = result.summary
    if summary is None:
        summary = {"description": result.description}

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
