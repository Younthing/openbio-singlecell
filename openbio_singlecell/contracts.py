from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import TYPE_CHECKING, Any, Literal

from . import PLUGIN_VERSION, SCHEMA_VERSION

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


METADATA_KEY = "openbio_singlecell"
_MISSING = object()
_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


@dataclass(frozen=True, slots=True)
class _SingleCellResultBase:
    title: str
    parameters: dict[str, Any]
    description: str
    warnings: list[str]
    input_cells: int
    input_genes: int
    random_seed: int
    elapsed_seconds: float
    source: dict[str, Any]


@dataclass(frozen=True, slots=True)
class SummaryResult(_SingleCellResultBase):
    summary: Any
    kind: Literal["summary"] = field(init=False, default="summary")

    def __post_init__(self) -> None:
        if self.summary is None:
            raise ValueError("SummaryResult requires summary data.")


@dataclass(frozen=True, slots=True)
class TableResult(_SingleCellResultBase):
    table: DataFrame
    kind: Literal["table"] = field(init=False, default="table")

    def __post_init__(self) -> None:
        from pandas import DataFrame as PandasDataFrame

        if not isinstance(self.table, PandasDataFrame):
            raise TypeError("TableResult table data must be a pandas DataFrame.")


@dataclass(frozen=True, slots=True)
class PlotResult(_SingleCellResultBase):
    png: bytes
    kind: Literal["plot"] = field(init=False, default="plot")

    def __post_init__(self) -> None:
        if not isinstance(self.png, bytes):
            raise TypeError("PlotResult PNG data must be bytes.")
        if not self.png.startswith(_PNG_SIGNATURE):
            raise ValueError("PlotResult PNG data must start with the standard PNG signature.")


type SingleCellResult = SummaryResult | TableResult | PlotResult


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
    if not isinstance(value, dict):
        raise ValueError("OpenBio metadata analysis_history must be a mapping.")
    converted = _metadata_value(value)
    if not isinstance(converted, dict):
        raise ValueError("OpenBio metadata analysis_history must be a mapping.")
    return converted


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
    if METADATA_KEY in adata.uns:
        existing = adata.uns[METADATA_KEY]
        if not isinstance(existing, dict):
            raise ValueError("OpenBio metadata must be a mapping.")
        if existing.get("schema_version") != SCHEMA_VERSION:
            raise ValueError(
                "OpenBio metadata has unsupported schema_version "
                f"{existing.get('schema_version')!r}; expected {SCHEMA_VERSION}."
            )
        metadata = copy.deepcopy(existing)
    else:
        metadata = {}
    metadata["schema_version"] = SCHEMA_VERSION
    metadata["version"] = PLUGIN_VERSION
    if display_name is None:
        metadata.setdefault("display_name", "AnnData")
    else:
        metadata["display_name"] = display_name
    if source is None:
        metadata.setdefault("source", {})
    else:
        metadata["source"] = _metadata_value(source)
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
