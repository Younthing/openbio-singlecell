from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .artifact_codecs import read_table, write_table
from .dgidb_resource import DGIdbResource, validate_dgidb_resource

DGIDB_CODEC = "dgidb-table-jsonl-v1"


def write_dgidb_resource(root: str | Path, resource: DGIdbResource) -> list[dict[str, str | int]]:
    if type(resource) is not DGIdbResource:
        raise TypeError("Expected an exact DGIdbResource value.")
    # The worker owns the freshly loaded resource; encode that state directly instead of calling copy-returning getters.
    table = object.__getattribute__(resource, "_table")
    metadata = object.__getattribute__(resource, "_metadata")
    accounting = object.__getattribute__(resource, "_accounting")
    artifact = object.__getattribute__(resource, "_artifact_metadata")
    return write_table(
        root,
        table,
        {
            "kind": "dgidb-resource",
            "metadata": metadata,
            "accounting": accounting,
            "artifact_metadata": artifact,
        },
    )


def read_dgidb_resource(root: str | Path) -> DGIdbResource:
    table, envelope = read_table(root)
    expected = {"kind", "metadata", "accounting", "artifact_metadata"}
    if not isinstance(envelope, Mapping) or set(envelope) != expected or envelope["kind"] != "dgidb-resource":
        raise ValueError("DGIdb artifact metadata schema is invalid.")
    resource = DGIdbResource._from_owned(
        table=table,
        metadata=_mapping(envelope["metadata"], "metadata"),
        accounting=_mapping(envelope["accounting"], "accounting"),
        artifact_metadata=_mapping(envelope["artifact_metadata"], "artifact metadata"),
    )
    validate_dgidb_resource(resource, copy_payload=False)
    return resource


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"DGIdb {name} must be a JSON object.")
    return value


__all__ = ["DGIDB_CODEC", "read_dgidb_resource", "write_dgidb_resource"]
