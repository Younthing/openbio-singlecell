from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import unicodedata
import uuid
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

from .artifact_runtime import (
    NATIVE_AFFINITY_CODECS,
    ArtifactRuntime,
    ArtifactTicket,
    _is_link_or_reparse,
)

PERSISTED_MANIFEST_VERSION = 1
_INVALID_COMPONENT = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
_WINDOWS_RESERVED = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{index}" for index in range(1, 10)),
    *(f"LPT{index}" for index in range(1, 10)),
}


@dataclass(frozen=True, slots=True)
class PersistedArtifact:
    path: Path


def _safe_component(value: object, description: str, *, normalize: bool = False) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{description} must be a string.")
    candidate = unicodedata.normalize("NFKC", value).strip()
    if not normalize and candidate != value:
        raise ValueError(f"{description} is not a safe path component.")
    if (
        not candidate
        or candidate in {".", ".."}
        or candidate.endswith((" ", "."))
        or _INVALID_COMPONENT.search(candidate)
        or candidate.split(".", 1)[0].upper() in _WINDOWS_RESERVED
    ):
        raise ValueError(f"{description} is not a safe path component.")
    return candidate


def _safe_relative(value: str) -> Path:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError("Persisted artifact file path is invalid.")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or path.as_posix() != value
        or not path.parts
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError("Persisted artifact file path is invalid.")
    for part in path.parts:
        _safe_component(part, "Persisted artifact file name")
    return Path(*path.parts)


def _hash_file(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
            size += len(chunk)
    return size, digest.hexdigest()


def _copy_file(source: Path, destination: Path) -> dict[str, str | int]:
    digest = hashlib.sha256()
    size = 0
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as source_handle, destination.open("xb") as destination_handle:
        while chunk := source_handle.read(1024 * 1024):
            destination_handle.write(chunk)
            digest.update(chunk)
            size += len(chunk)
        destination_handle.flush()
        os.fsync(destination_handle.fileno())
    _fsync_directory(destination.parent)
    return {"size": size, "sha256": digest.hexdigest()}


def _source_files(source: Path, payload_type: str) -> list[tuple[str, Path]]:
    if _is_link_or_reparse(source):
        raise ValueError("Portable artifact payload cannot be a link or reparse point.")
    resolved_source = source.resolve(strict=True)
    if payload_type == "file":
        _safe_relative(source.name)
        return [(source.name, source)]
    files = []
    for candidate in source.rglob("*"):
        if _is_link_or_reparse(candidate):
            raise ValueError("Portable artifacts cannot contain links or reparse points.")
        try:
            candidate.resolve(strict=True).relative_to(resolved_source)
        except ValueError as error:
            raise ValueError("Portable artifact payload escapes its root.") from error
        if candidate.is_file():
            relative = candidate.relative_to(source).as_posix()
            _safe_relative(relative)
            files.append((relative, candidate))
    files.sort(key=lambda item: item[0])
    return files


def _write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, allow_nan=False, sort_keys=True, separators=(",", ":"))
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _fsync_directory(path: Path) -> None:
    if os.name == "nt":
        return
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _load_manifest(path: Path) -> dict[str, Any]:
    def reject_constant(value: str) -> None:
        raise ValueError(f"Persisted manifest contains an invalid JSON number: {value}")

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Persisted manifest contains a duplicate key: {key}")
            result[key] = value
        return result

    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=reject_constant,
            object_pairs_hook=unique_object,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("Persisted artifact manifest is unreadable.") from error
    if not isinstance(value, dict) or set(value) != {"version", "artifact", "payload_type", "files"}:
        raise ValueError("Persisted artifact manifest schema is invalid.")
    artifact = value["artifact"]
    if not isinstance(artifact, dict) or set(artifact) != {"session", "run", "output", "kind", "codec"}:
        raise ValueError("Persisted artifact manifest schema is invalid.")
    if (
        type(value["version"]) is not int
        or value["version"] != PERSISTED_MANIFEST_VERSION
        or value["payload_type"] not in {"file", "directory"}
    ):
        raise ValueError("Persisted artifact manifest schema is invalid.")
    if not all(isinstance(item, str) and item for item in artifact.values()):
        raise ValueError("Persisted artifact manifest schema is invalid.")
    if not isinstance(value["files"], list):
        raise ValueError("Persisted artifact manifest schema is invalid.")
    seen = set()
    for item in value["files"]:
        if (
            not isinstance(item, dict)
            or set(item) != {"path", "size", "sha256"}
            or not isinstance(item["path"], str)
            or not isinstance(item["size"], int)
            or isinstance(item["size"], bool)
            or item["size"] < 0
            or not isinstance(item["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", item["sha256"]) is None
        ):
            raise ValueError("Persisted artifact manifest schema is invalid.")
        _safe_relative(item["path"])
        if item["path"] in seen:
            raise ValueError("Persisted artifact manifest schema is invalid.")
        seen.add(item["path"])
    return value


def _validate_persisted(root: Path, expected: dict[str, Any]) -> dict[str, Any]:
    if _is_link_or_reparse(root) or not root.is_dir():
        raise ValueError("Persisted artifact target is not a regular directory.")
    resolved_root = root.resolve(strict=True)
    manifest = _load_manifest(root / "manifest.json")
    if manifest != expected:
        raise ValueError("Persisted artifact manifest does not match its source artifact.")
    expected_paths = {item["path"] for item in manifest["files"]}
    observed_paths = set()
    for candidate in root.rglob("*"):
        if _is_link_or_reparse(candidate):
            raise ValueError("Persisted artifacts cannot contain links or reparse points.")
        try:
            candidate.resolve(strict=True).relative_to(resolved_root)
        except ValueError as error:
            raise ValueError("Persisted artifact file escapes its root.") from error
        if candidate.is_file() and candidate != root / "manifest.json":
            observed_paths.add(candidate.relative_to(root).as_posix())
    if observed_paths != expected_paths:
        raise ValueError("Persisted artifact files do not match their manifest.")
    for item in manifest["files"]:
        size, digest = _hash_file(root / _safe_relative(item["path"]))
        if size != item["size"] or digest != item["sha256"]:
            raise ValueError("Persisted artifact file hash does not match its manifest.")
    return manifest


def _manifest(ticket: ArtifactTicket, payload_type: object, files: list[dict[str, str | int]]) -> dict[str, Any]:
    return {
        "version": PERSISTED_MANIFEST_VERSION,
        "artifact": {
            "session": ticket.session,
            "run": ticket.run,
            "output": ticket.output,
            "kind": ticket.kind,
            "codec": ticket.codec,
        },
        "payload_type": payload_type,
        "files": files,
    }


def persist_artifact(
    runtime: ArtifactRuntime,
    ticket: ArtifactTicket,
    output_root: str | os.PathLike[str],
    name: str,
) -> PersistedArtifact:
    source = runtime.resolve(ticket)
    record = runtime.record(ticket)
    if ticket.codec in NATIVE_AFFINITY_CODECS:
        raise ValueError("scVI and cNMF session-only native artifacts cannot be persisted.")
    safe_name = _safe_component(name, "Artifact name", normalize=True)
    run = _safe_component(ticket.run, "Artifact run")
    output = _safe_component(ticket.output, "Artifact output")
    root = Path(output_root).resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError(f"Output root is not a directory: {root}")
    artifacts_root = root / "openbio-singlecell" / "artifacts"
    artifacts_root.mkdir(parents=True, exist_ok=True)
    artifacts_root = artifacts_root.resolve(strict=True)
    try:
        artifacts_root.relative_to(root)
    except ValueError as error:
        raise ValueError("Persisted artifact directory escapes the output root.") from error
    target = artifacts_root / f"{safe_name}__{run}-{output}"
    source_files = _source_files(source, str(record["payload_type"]))
    if os.path.lexists(target):
        files = []
        for relative, source_file in source_files:
            size, digest = _hash_file(source_file)
            files.append({"path": relative, "size": size, "sha256": digest})
        expected = _manifest(ticket, record["payload_type"], files)
        try:
            _validate_persisted(target, expected)
        except (OSError, ValueError) as error:
            raise FileExistsError(f"Persisted artifact conflicts with existing target: {target.name}") from error
        return PersistedArtifact(target)

    while True:
        staging = artifacts_root / f".{target.name}.{uuid.uuid4()}.partial"
        try:
            staging.mkdir(exist_ok=False)
        except FileExistsError:
            continue
        break
    try:
        files = []
        for relative, source_file in source_files:
            if relative == "manifest.json":
                raise ValueError("Portable artifact payload reserves manifest.json for persistence metadata.")
            copied = _copy_file(source_file, staging / _safe_relative(relative))
            files.append({"path": relative, **copied})
        manifest = _manifest(ticket, record["payload_type"], files)
        _write_manifest(staging / "manifest.json", manifest)
        _validate_persisted(staging, manifest)
        _fsync_directory(staging)
        try:
            os.rename(staging, target)
        except OSError:
            if not os.path.lexists(target):
                raise
            try:
                _validate_persisted(target, manifest)
            except (OSError, ValueError) as error:
                raise FileExistsError(
                    f"Persisted artifact conflicts with existing target: {target.name}"
                ) from error
        _fsync_directory(artifacts_root)
        return PersistedArtifact(target)
    finally:
        try:
            shutil.rmtree(staging)
        except FileNotFoundError:
            pass
