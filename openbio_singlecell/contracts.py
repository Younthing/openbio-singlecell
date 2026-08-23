from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from datetime import UTC, datetime
from typing import TYPE_CHECKING, Any, Literal

from . import PLUGIN_VERSION, SCHEMA_VERSION

if TYPE_CHECKING:
    from anndata import AnnData


METADATA_KEY = "openbio_singlecell"
ResultKind = Literal["summary", "table", "plot"]
_MISSING = object()


@dataclass(slots=True)
class SingleCellResult:
    kind: ResultKind
    title: str
    parameters: dict[str, Any]
    description: str
    warnings: list[str]
    input_cells: int
    input_genes: int
    random_seed: int
    elapsed_seconds: float
    source: dict[str, Any]
    table: Any = None
    png: bytes | None = None
    summary: Any = None


def _metadata_value(value: Any) -> Any:
    if value is None:
        return _MISSING
    if isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    if hasattr(value, "item"):
        try:
            return _metadata_value(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            converted = _metadata_value(item)
            if converted is not _MISSING:
                result[str(key)] = converted
        return result
    if isinstance(value, (list, tuple)):
        result = []
        for item in value:
            converted = _metadata_value(item)
            if converted is not _MISSING:
                result.append(converted)
        return result
    return str(value)


def _normalize_history(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        converted = _metadata_value(value)
        return converted if isinstance(converted, dict) else {}

    if hasattr(value, "tolist"):
        value = value.tolist()
    if not isinstance(value, (list, tuple)):
        return {}

    history = {}
    for index, item in enumerate(value):
        converted = _metadata_value(item)
        if converted is not _MISSING:
            history[f"{index:06d}"] = converted
    return history


def _normalize_warnings(value: Any) -> list[str]:
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, str):
        return [value]
    if not isinstance(value, (list, tuple)):
        return []
    return [str(warning) for warning in value]


def _next_history_key(history: dict[str, Any]) -> str:
    indices = []
    for key in history:
        try:
            indices.append(int(key))
        except (TypeError, ValueError):
            continue

    index = max(indices, default=-1) + 1
    key = f"{index:06d}"
    while key in history:
        index += 1
        key = f"{index:06d}"
    return key


def make_analysis_source(
    operation: str,
    parameters: dict[str, Any],
    input_cells: int,
    input_genes: int,
    random_seed: int = 0,
    elapsed_seconds: float = 0.0,
) -> dict[str, Any]:
    return {
        "operation": operation,
        "parameters": _metadata_value(parameters),
        "input_cells": int(input_cells),
        "input_genes": int(input_genes),
        "random_seed": int(random_seed),
        "elapsed_seconds": float(elapsed_seconds),
        "timestamp": datetime.now(UTC).isoformat(),
    }


def ensure_metadata(
    adata: AnnData,
    display_name: str | None = None,
    source: dict[str, Any] | None = None,
) -> dict[str, Any]:
    existing = adata.uns.get(METADATA_KEY)
    metadata = copy.deepcopy(existing) if isinstance(existing, dict) else {}
    metadata.setdefault("schema_version", SCHEMA_VERSION)
    metadata.setdefault("version", PLUGIN_VERSION)
    metadata.setdefault("display_name", display_name or "AnnData")
    metadata.setdefault("source", _metadata_value(source or {}))
    try:
        metadata["random_seed"] = int(metadata.get("random_seed", 0))
    except (TypeError, ValueError, OverflowError):
        metadata["random_seed"] = 0
    metadata["warnings"] = _normalize_warnings(metadata.get("warnings", []))
    metadata["analysis_history"] = _normalize_history(metadata.get("analysis_history", {}))
    adata.uns[METADATA_KEY] = metadata
    return metadata


def record_history(
    adata: AnnData,
    operation: str,
    parameters: dict[str, Any],
    input_cells: int,
    input_genes: int,
    random_seed: int = 0,
    elapsed_seconds: float = 0.0,
    warnings: list[str] | None = None,
) -> dict[str, Any]:
    metadata = ensure_metadata(adata)
    entry = make_analysis_source(
        operation,
        parameters,
        input_cells,
        input_genes,
        random_seed,
        elapsed_seconds,
    )
    history = dict(metadata["analysis_history"])
    history[_next_history_key(history)] = entry
    metadata["analysis_history"] = history
    metadata["random_seed"] = int(random_seed)
    if warnings:
        metadata["warnings"] = list(metadata.get("warnings", [])) + [str(w) for w in warnings]
    adata.uns[METADATA_KEY] = metadata
    return entry
