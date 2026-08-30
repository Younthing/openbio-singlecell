from __future__ import annotations

import json
import math
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .artifact_codecs import read_anndata, write_anndata

VELOCITY_STATE_CODEC = "velocity-state-h5ad-json-v1"
CNV_STATE_CODEC = "cnv-state-h5ad-json-v1"
STATE_PAYLOAD = "state.json"


def _root(directory: str | os.PathLike[str]) -> Path:
    root = Path(directory).resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError(f"Staged-state artifact root is not a directory: {root}")
    return root


def _strict_json(value: Any, *, location: str) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{location} contains a non-finite JSON number.")
        return value
    if isinstance(value, list):
        return [_strict_json(item, location=f"{location}[]") for item in value]
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise TypeError(f"{location} JSON object keys must be strings.")
        return {key: _strict_json(item, location=f"{location}.{key}") for key, item in value.items()}
    raise TypeError(f"{location} contains a non-JSON value of type {type(value).__name__}.")


def _state_path(root: str | os.PathLike[str]) -> Path:
    return _root(root) / STATE_PAYLOAD


def _write_state(root: str | os.PathLike[str], value: Mapping[str, Any]) -> dict[str, str | int]:
    path = _state_path(root)
    if os.path.lexists(path):
        raise FileExistsError(f"Artifact payload already exists: {STATE_PAYLOAD}")
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(
            _strict_json(dict(value), location=STATE_PAYLOAD),
            handle,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        handle.write("\n")
    return {"path": STATE_PAYLOAD, "size": path.stat().st_size}


def _reject_constant(value: str) -> None:
    raise ValueError(f"Staged-state JSON contains a non-finite number: {value}")


def _read_state(root: str | os.PathLike[str]) -> dict[str, Any]:
    path = _state_path(root)
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=_reject_constant)
    value = _strict_json(value, location=STATE_PAYLOAD)
    if not isinstance(value, dict):
        raise TypeError("Staged-state metadata must be a JSON object.")
    return value


def write_velocity_state(root: str | os.PathLike[str], state: Any) -> list[dict[str, str | int]]:
    from .velocity_analysis import VelocityState, validate_velocity_state
    from .velocity_portable import VELOCITY_STATE_KEY

    if type(state) is not VelocityState:
        raise TypeError("Expected an exact worker-private VelocityState.")
    adata, summary, metadata, _portable = validate_velocity_state(state, _owned=True)
    portable_state = adata.uns.pop(VELOCITY_STATE_KEY)
    try:
        descriptors = write_anndata(root, adata)
    finally:
        adata.uns[VELOCITY_STATE_KEY] = portable_state
    descriptors.append(
        _write_state(
            root,
            {
                "version": 1,
                "kind": "OPENBIO_VELOCITY_STATE",
                "codec": VELOCITY_STATE_CODEC,
                "metadata": metadata,
                "summary": summary,
                "portable_uns": {VELOCITY_STATE_KEY: portable_state},
            },
        )
    )
    return descriptors


def read_velocity_state(root: str | os.PathLike[str]) -> Any:
    from .velocity_analysis import VelocityState, validate_velocity_state
    from .velocity_portable import VELOCITY_STATE_KEY

    record = _read_state(root)
    if set(record) != {
        "version",
        "kind",
        "codec",
        "metadata",
        "summary",
        "portable_uns",
    }:
        raise ValueError("Velocity state JSON has an invalid exact schema.")
    if (
        record["version"] != 1
        or record["kind"] != "OPENBIO_VELOCITY_STATE"
        or record["codec"] != VELOCITY_STATE_CODEC
        or not isinstance(record["metadata"], dict)
        or not isinstance(record["summary"], dict)
        or not isinstance(record["portable_uns"], dict)
        or set(record["portable_uns"]) != {VELOCITY_STATE_KEY}
        or not isinstance(record["portable_uns"][VELOCITY_STATE_KEY], dict)
    ):
        raise ValueError("Velocity state JSON identity is invalid.")
    adata = read_anndata(root)
    if VELOCITY_STATE_KEY in adata.uns:
        raise ValueError("Velocity H5AD unexpectedly duplicates sidecar-owned state provenance.")
    adata.uns[VELOCITY_STATE_KEY] = record["portable_uns"][VELOCITY_STATE_KEY]
    state = VelocityState._from_owned(
        adata=adata,
        metadata=record["metadata"],
        summary=record["summary"],
    )
    validate_velocity_state(state, _owned=True)
    return state


def write_cnv_state(root: str | os.PathLike[str], state: Any) -> list[dict[str, str | int]]:
    from scipy import sparse

    from .cnv_analysis import _CNV_PROVENANCE_KEY, CNVState, _cnv_validate_state

    if type(state) is not CNVState:
        raise TypeError("Expected an exact worker-private CNVState.")
    import numpy as np

    adata, metadata = _cnv_validate_state(state, numpy=np, scipy_sparse=sparse, _owned=True)
    # AnnData/HDF5 cannot encode the ordered list-of-records used by the CNV
    # provenance. It belongs in the strict JSON sidecar and is restored before
    # validating any worker-private consumer state.
    provenance = adata.uns.pop(_CNV_PROVENANCE_KEY)
    try:
        descriptors = write_anndata(root, adata)
    finally:
        adata.uns[_CNV_PROVENANCE_KEY] = provenance
    descriptors.append(
        _write_state(
            root,
            {
                "version": 1,
                "kind": "OPENBIO_CNV_STATE",
                "codec": CNV_STATE_CODEC,
                "metadata": metadata,
                "portable_uns": {_CNV_PROVENANCE_KEY: provenance},
            },
        )
    )
    return descriptors


def read_cnv_state(root: str | os.PathLike[str]) -> Any:
    import numpy as np
    from scipy import sparse

    from .cnv_analysis import _CNV_PROVENANCE_KEY, CNVState, _cnv_validate_state

    record = _read_state(root)
    if set(record) != {"version", "kind", "codec", "metadata", "portable_uns"}:
        raise ValueError("CNV state JSON has an invalid exact schema.")
    if (
        record["version"] != 1
        or record["kind"] != "OPENBIO_CNV_STATE"
        or record["codec"] != CNV_STATE_CODEC
        or not isinstance(record["metadata"], dict)
        or not isinstance(record["portable_uns"], dict)
        or set(record["portable_uns"]) != {_CNV_PROVENANCE_KEY}
        or not isinstance(record["portable_uns"][_CNV_PROVENANCE_KEY], dict)
    ):
        raise ValueError("CNV state JSON identity is invalid.")
    adata = read_anndata(root)
    if _CNV_PROVENANCE_KEY in adata.uns:
        raise ValueError("CNV H5AD unexpectedly duplicates sidecar-owned state provenance.")
    adata.uns[_CNV_PROVENANCE_KEY] = record["portable_uns"][_CNV_PROVENANCE_KEY]
    state = CNVState(adata, record["metadata"], _owned=True)
    _cnv_validate_state(state, numpy=np, scipy_sparse=sparse, _owned=True)
    return state


__all__ = [
    "CNV_STATE_CODEC",
    "STATE_PAYLOAD",
    "VELOCITY_STATE_CODEC",
    "read_cnv_state",
    "read_velocity_state",
    "write_cnv_state",
    "write_velocity_state",
]
