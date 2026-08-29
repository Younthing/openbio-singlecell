from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .artifact_codecs import read_table, write_table
from .artifact_envelope import write_strict_json
from .augur import AUGUR_VIEWS, validate_augur_portable
from .worker_protocol import read_json

AUGUR_CODEC = "augur-jsonl-v1"
AUGUR_SCHEMA = "schema.json"
AUGUR_RESULT = "result.json"
_FORMAT = "openbio-augur-jsonl"


def _root(path: str | Path) -> Path:
    root = Path(path).resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError(f"Augur artifact root is not a directory: {root}")
    return root


def write_augur(
    root: str | Path,
    tables: Mapping[str, Any],
    summary: Mapping[str, Any],
    metadata: Mapping[str, Any],
) -> None:
    artifact_root = _root(root)
    owned_tables, owned_summary, owned_metadata = validate_augur_portable(tables, summary, metadata)
    table_root = artifact_root / "tables"
    table_root.mkdir(exist_ok=False)
    for view in AUGUR_VIEWS:
        view_root = table_root / view
        view_root.mkdir(exist_ok=False)
        write_table(view_root, owned_tables[view], {"view": view})
    write_strict_json(
        artifact_root / AUGUR_SCHEMA,
        {"version": 1, "format": _FORMAT, "views": list(AUGUR_VIEWS)},
    )
    write_strict_json(
        artifact_root / AUGUR_RESULT,
        {"summary": owned_summary, "metadata": owned_metadata},
    )


def read_augur(root: str | Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    artifact_root = _root(root)
    schema = read_json(artifact_root / AUGUR_SCHEMA)
    expected_schema = {"version": 1, "format": _FORMAT, "views": list(AUGUR_VIEWS)}
    if schema != expected_schema:
        raise ValueError("Augur artifact codec schema is invalid.")
    result = read_json(artifact_root / AUGUR_RESULT)
    if not isinstance(result, dict) or set(result) != {"summary", "metadata"}:
        raise ValueError("Augur artifact result JSON is invalid.")
    tables = {}
    for view in AUGUR_VIEWS:
        table, view_metadata = read_table(artifact_root / "tables" / view)
        if view_metadata != {"view": view}:
            raise ValueError(f"Augur artifact table metadata is invalid for {view!r}.")
        tables[view] = table
    summary = result["summary"]
    metadata = result["metadata"]
    if not isinstance(summary, dict) or not isinstance(metadata, dict):
        raise ValueError("Augur artifact summary and metadata must be JSON objects.")
    return validate_augur_portable(tables, summary, metadata)


__all__ = ["AUGUR_CODEC", "read_augur", "write_augur"]
