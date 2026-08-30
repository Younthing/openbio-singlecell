from __future__ import annotations

from pathlib import Path
from typing import Any

from .artifact_codecs import ANNDATA_PAYLOAD, read_anndata, write_anndata
from .pseudobulk import (
    _standalone_plain_json,
    _standalone_validate_pseudobulk_artifact,
    _StandalonePseudobulkArtifact,
)
from .worker_protocol import read_json, write_json

PSEUDOBULK_CODEC = "pseudobulk-h5ad-v1"
PSEUDOBULK_KIND = "OPENBIO_SINGLE_CELL_PSEUDOBULK"
PSEUDOBULK_METADATA = "metadata.json"


def _root(root: str | Path) -> Path:
    path = Path(root).resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(f"Pseudobulk artifact root is not a directory: {path}")
    return path


def write_pseudobulk(root: str | Path, artifact: Any) -> list[dict[str, str | int]]:
    path = _root(root)
    adata, metadata = _standalone_validate_pseudobulk_artifact(artifact, copy_result=False)
    write_anndata(path, adata)
    metadata_path = path / PSEUDOBULK_METADATA
    if metadata_path.exists():
        raise FileExistsError(f"Pseudobulk artifact metadata already exists: {metadata_path}")
    write_json(metadata_path, _standalone_plain_json(metadata))
    return [
        {"path": ANNDATA_PAYLOAD, "size": (path / ANNDATA_PAYLOAD).stat().st_size},
        {"path": PSEUDOBULK_METADATA, "size": metadata_path.stat().st_size},
    ]


def read_pseudobulk(root: str | Path) -> Any:
    path = _root(root)
    metadata = read_json(path / PSEUDOBULK_METADATA)
    if not isinstance(metadata, dict):
        raise TypeError("Pseudobulk artifact metadata must be a JSON object.")
    artifact = _StandalonePseudobulkArtifact._from_owned(read_anndata(path), metadata)
    _standalone_validate_pseudobulk_artifact(artifact, copy_result=False)
    return artifact


__all__ = [
    "PSEUDOBULK_CODEC",
    "PSEUDOBULK_KIND",
    "PSEUDOBULK_METADATA",
    "read_pseudobulk",
    "write_pseudobulk",
]
