from __future__ import annotations

import ast
import copy
import csv
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
import platform
import re
import struct
import threading
from collections.abc import Sequence
from pathlib import Path, PurePosixPath
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .scenic_artifact import (
    SCENIC_MEMBERSHIP_COLUMNS,
    SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
    SCENIC_RESULT_ARTIFACT_TYPE,
    SCENIC_RESULT_PRODUCER_NODE_ID,
    SCENIC_RESULT_PRODUCER_SCHEMA,
    SCENICResultArtifact,
    _activity_bytes,
    _canonical_axis,
    _canonical_json,
    _canonical_membership,
    _result_fingerprint,
    _strict_provenance,
)

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


PYSCENIC_VERSION = "0.12.1"
PYSCENIC_MANIFEST_SCHEMA = "openbio-singlecell/pyscenic-external-run/v1"
PYSCENIC_IMPORT_SUMMARY_SCHEMA = "openbio-singlecell/pyscenic-import-summary/v1"
SCENIC_ACTIVITY_KEY = "scenic_auc"
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_CONTAINER_DIGEST_PATTERN = re.compile(r"^sha256:[0-9a-f]{64}$")
_REGULON_PATTERN = re.compile(r"^(.+)\(([+-])\)$")
_CSV_LIMIT_LOCK = threading.Lock()
_MOTIF_COLUMNS = (
    "AUC",
    "NES",
    "MotifSimilarityQvalue",
    "OrthologousIdentity",
    "Annotation",
    "Context",
    "TargetGenes",
    "RankAtMax",
)
PYSCENIC_REFERENCES = [
    {
        "citation": (
            "Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network "
            "analysis. Nature Protocols. 2020;15:2247-2276."
        ),
        "doi": "10.1038/s41596-020-0336-2",
        "url": "https://doi.org/10.1038/s41596-020-0336-2",
        "kind": "method",
    },
    {
        "citation": (
            "Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. "
            "Nature Methods. 2017;14:1083-1086."
        ),
        "doi": "10.1038/nmeth.4463",
        "url": "https://doi.org/10.1038/nmeth.4463",
        "kind": "method",
    },
    {
        "citation": (
            "Moerman T, et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory "
            "networks. Bioinformatics. 2019;35:2159-2161."
        ),
        "doi": "10.1093/bioinformatics/bty916",
        "url": "https://doi.org/10.1093/bioinformatics/bty916",
        "kind": "method",
    },
    {
        "citation": "pySCENIC 0.12.1 tagged CLI, motif-table loader, regulon conversion, and AUCell implementation.",
        "doi": None,
        "url": "https://github.com/aertslab/pySCENIC/tree/0.12.1",
        "kind": "software_documentation",
    },
    {
        "citation": "Virshup I, et al. anndata: Annotated data. Journal of Open Source Software. 2024;9:4371.",
        "doi": "10.21105/joss.04371",
        "url": "https://doi.org/10.21105/joss.04371",
        "kind": "software",
    },
]


def _strict_json_loads(text: str, *, description: str) -> Any:
    def reject_constant(value: str) -> None:
        raise ValueError(f"{description} contains non-finite JSON constant {value!r}.")

    def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{description} contains duplicate JSON key {key!r}.")
            result[key] = value
        return result

    try:
        return json.loads(
            text,
            parse_constant=reject_constant,
            object_pairs_hook=reject_duplicate_keys,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"{description} is not valid strict JSON ({exc}).") from exc


def _require_fields(value: Any, expected: set[str], *, description: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError(f"{description} must be a JSON object.")
    observed = set(value)
    if observed != expected:
        raise ValueError(
            f"{description} fields are invalid; missing={sorted(expected - observed)!r}, "
            f"unknown={sorted(observed - expected)!r}."
        )
    return value


def _canonical_text(value: Any, *, description: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip() or "\x00" in value:
        raise ValueError(f"{description} must be a canonical nonblank string.")
    return value


def _sha256_text(value: Any, *, description: str) -> str:
    value = _canonical_text(value, description=description)
    if _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{description} must be a lowercase SHA-256 digest.")
    return value


def _positive_int(value: Any, *, description: str, maximum: int = 2**63 - 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise ValueError(f"{description} must be an integer in [1, {maximum}].")
    return value


def _sha256_file(path: Path, *, max_file_bytes: int) -> tuple[str, int]:
    if path.is_symlink():
        raise ValueError(f"pySCENIC bundle members cannot be symlinks: {path}.")
    try:
        stat_before = path.stat()
    except OSError as exc:
        raise ValueError(f"Cannot stat pySCENIC bundle member {path}: {exc}.") from exc
    if not path.is_file():
        raise ValueError(f"pySCENIC bundle member is not a regular file: {path}.")
    if stat_before.st_size < 1:
        raise ValueError(f"pySCENIC bundle member is empty: {path}.")
    if stat_before.st_size > max_file_bytes:
        raise ValueError(
            f"pySCENIC bundle member {path.name!r} is {stat_before.st_size:,} bytes, exceeding "
            f"max_file_bytes={max_file_bytes:,}."
        )
    digest = hashlib.sha256()
    observed_size = 0
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                observed_size += len(chunk)
                if observed_size > max_file_bytes:
                    raise ValueError(
                        f"pySCENIC bundle member {path.name!r} grew beyond max_file_bytes during hashing."
                    )
                digest.update(chunk)
    except OSError as exc:
        raise ValueError(f"Cannot hash pySCENIC bundle member {path}: {exc}.") from exc
    stat_after = path.stat()
    identity_before = (stat_before.st_dev, stat_before.st_ino, stat_before.st_size, stat_before.st_mtime_ns)
    identity_after = (stat_after.st_dev, stat_after.st_ino, stat_after.st_size, stat_after.st_mtime_ns)
    if identity_before != identity_after or observed_size != stat_after.st_size:
        raise ValueError(f"pySCENIC bundle member changed during hashing: {path}.")
    return digest.hexdigest(), observed_size


def _resolve_bundle_member(base: Path, relative_name: Any, *, description: str) -> Path:
    name = _canonical_text(relative_name, description=description)
    if "\\" in name or ":" in name:
        raise ValueError(f"{description} must use a normalized relative POSIX path.")
    pure = PurePosixPath(name)
    if pure.is_absolute() or any(part in {"", ".", ".."} for part in pure.parts):
        raise ValueError(f"{description} must be a normalized relative path without traversal.")
    candidate = base.joinpath(*pure.parts)
    current = base
    for part in pure.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"{description} cannot resolve through a symlink.")
    try:
        resolved_base = base.resolve(strict=True)
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(resolved_base)
    except (OSError, ValueError) as exc:
        raise ValueError(f"{description} escapes or is missing from the manifest bundle.") from exc
    if candidate.is_symlink() or resolved.is_symlink():
        raise ValueError(f"{description} cannot resolve through a symlink.")
    return resolved


def _command_flag(command: list[str], flag: str, *, stage: str) -> str:
    positions = [index for index, value in enumerate(command) if value == flag]
    if len(positions) != 1 or positions[0] + 1 >= len(command):
        raise ValueError(f"pySCENIC {stage} command must contain exactly one {flag} value.")
    return command[positions[0] + 1]


def _optional_command_flag(command: list[str], flag: str, *, stage: str) -> str | None:
    positions = [index for index, value in enumerate(command) if value == flag]
    if len(positions) > 1 or (positions and positions[0] + 1 >= len(command)):
        raise ValueError(f"pySCENIC {stage} command contains an invalid repeated or valueless {flag} option.")
    return command[positions[0] + 1] if positions else None


def _optional_command_flag_alias(
    command: list[str], flags: tuple[str, ...], *, stage: str
) -> str | None:
    positions = [index for index, value in enumerate(command) if value in flags]
    if len(positions) > 1 or (positions and positions[0] + 1 >= len(command)):
        raise ValueError(
            f"pySCENIC {stage} command contains invalid repeated or valueless options from {list(flags)!r}."
        )
    return command[positions[0] + 1] if positions else None


def _command_flag_alias(command: list[str], flags: tuple[str, ...], *, stage: str) -> str:
    positions = [index for index, value in enumerate(command) if value in flags]
    if len(positions) != 1 or positions[0] + 1 >= len(command):
        raise ValueError(
            f"pySCENIC {stage} command must contain exactly one of {list(flags)!r} with a value."
        )
    return command[positions[0] + 1]


def _parse_command(command: Any, *, stage: str) -> tuple[list[str], dict[str, Any]]:
    if not isinstance(command, list) or len(command) < 3:
        raise TypeError(f"pySCENIC {stage} command must be a JSON argument array.")
    values = []
    for position, value in enumerate(command):
        values.append(_canonical_text(value, description=f"pySCENIC {stage} argument {position}"))
    if values[:2] != ["pyscenic", stage]:
        raise ValueError(f"pySCENIC {stage} command must begin with ['pyscenic', {stage!r}].")
    workers_text = _optional_command_flag(values, "--num_workers", stage=stage)
    workers: int | None = None
    if workers_text is not None:
        try:
            workers = int(workers_text)
        except ValueError as exc:
            raise ValueError(f"pySCENIC {stage} --num_workers must be an integer.") from exc
        _positive_int(workers, description=f"pySCENIC {stage} --num_workers", maximum=65536)
    parsed: dict[str, Any] = {"num_workers": workers}
    if stage in {"grn", "aucell"}:
        seed_text = _optional_command_flag(values, "--seed", stage=stage)
        seed: int | None = None
        if seed_text is not None:
            try:
                seed = int(seed_text)
            except ValueError as exc:
                raise ValueError(f"pySCENIC {stage} --seed must be an integer.") from exc
            if not 0 <= seed <= 2**31 - 1:
                raise ValueError(f"pySCENIC {stage} --seed must lie in [0, 2^31-1].")
        parsed["seed"] = seed
    if stage == "grn":
        method = _optional_command_flag_alias(values, ("-m", "--method"), stage=stage) or "grnboost2"
        _command_flag_alias(values, ("-o", "--output"), stage=stage)
        if method not in {"grnboost2", "genie3"}:
            raise ValueError("pySCENIC grn --method must be grnboost2 or genie3.")
        parsed["method"] = method
    elif stage == "ctx":
        if "--no_pruning" in values or "-n" in values:
            raise ValueError("pySCENIC ctx must perform motif pruning for a final SCENIC artifact.")
        mode = _optional_command_flag(values, "--mode", stage=stage) or "custom_multiprocessing"
        if mode not in {"custom_multiprocessing", "dask_multiprocessing", "dask_cluster"}:
            raise ValueError("pySCENIC ctx --mode is invalid for 0.12.1.")
        _command_flag(values, "--annotations_fname", stage=stage)
        _command_flag(values, "--expression_mtx_fname", stage=stage)
        _command_flag_alias(values, ("-o", "--output"), stage=stage)
        parsed["mode"] = mode
        parsed["mask_dropouts"] = "--mask_dropouts" in values
    else:
        threshold_text = _optional_command_flag(values, "--auc_threshold", stage=stage)
        try:
            threshold = 0.05 if threshold_text is None else float(threshold_text)
        except ValueError as exc:
            raise ValueError("pySCENIC aucell --auc_threshold must be numeric.") from exc
        if not math.isfinite(threshold) or not 0.0 < threshold <= 1.0:
            raise ValueError("pySCENIC aucell --auc_threshold must lie in (0, 1].")
        _command_flag_alias(values, ("-o", "--output"), stage=stage)
        parsed["auc_threshold"] = threshold
    return values, parsed


def _validate_manifest(
    manifest: Any,
    *,
    base: Path,
) -> tuple[dict[str, Any], dict[str, Path], dict[str, Any]]:
    manifest = _require_fields(
        manifest,
        {
            "schema_version",
            "run_id",
            "organism",
            "genome_build",
            "gene_namespace",
            "expression_state",
            "container",
            "software_versions",
            "commands",
            "expression",
            "resources",
            "artifacts",
        },
        description="pySCENIC manifest",
    )
    if manifest["schema_version"] != PYSCENIC_MANIFEST_SCHEMA:
        raise ValueError(f"pySCENIC manifest schema must be exactly {PYSCENIC_MANIFEST_SCHEMA!r}.")
    for field_name in ("run_id", "organism", "genome_build", "gene_namespace"):
        _canonical_text(manifest[field_name], description=f"pySCENIC manifest {field_name}")
    _canonical_text(manifest["expression_state"], description="pySCENIC manifest expression_state")
    if manifest["container"] is not None:
        container = _require_fields(
            manifest["container"], {"image", "digest"}, description="pySCENIC container"
        )
        _canonical_text(container["image"], description="pySCENIC container image")
        digest = _canonical_text(container["digest"], description="pySCENIC container digest")
        if _CONTAINER_DIGEST_PATTERN.fullmatch(digest) is None:
            raise ValueError("pySCENIC container digest must be an immutable sha256:<64 hex> digest.")
    versions = manifest["software_versions"]
    if not isinstance(versions, dict):
        raise TypeError("pySCENIC declared software versions must be a JSON object.")
    required_versions = {"pyscenic", "ctxcore", "arboreto", "python"}
    if not required_versions.issubset(versions):
        raise ValueError(
            f"pySCENIC declared software versions are missing {sorted(required_versions - set(versions))!r}."
        )
    for name, version in versions.items():
        _canonical_text(name, description="pySCENIC declared software package name")
        _canonical_text(version, description=f"pySCENIC declared {name} version")
    if versions["pyscenic"] != PYSCENIC_VERSION:
        raise ValueError(f"Imported results require pySCENIC exactly {PYSCENIC_VERSION}.")
    commands = _require_fields(
        manifest["commands"], {"grn", "ctx", "aucell"}, description="pySCENIC commands"
    )
    normalized_commands: dict[str, list[str]] = {}
    command_parameters: dict[str, dict[str, Any]] = {}
    for stage in ("grn", "ctx", "aucell"):
        normalized_commands[stage], command_parameters[stage] = _parse_command(
            commands[stage], stage=stage
        )
    expression = _require_fields(
        manifest["expression"],
        {
            "file",
            "sha256",
            "cells",
            "genes",
            "observation_ids_sha256",
            "gene_ids_sha256",
            "matrix_sha256",
        },
        description="pySCENIC expression declaration",
    )
    cells = _positive_int(expression["cells"], description="pySCENIC expression cells")
    genes = _positive_int(expression["genes"], description="pySCENIC expression genes")
    for field_name in ("sha256", "observation_ids_sha256", "gene_ids_sha256", "matrix_sha256"):
        _sha256_text(expression[field_name], description=f"pySCENIC expression {field_name}")
    artifacts = _require_fields(
        manifest["artifacts"],
        {"adjacency", "regulons", "aucell"},
        description="pySCENIC artifacts",
    )
    declarations: list[tuple[str, dict[str, Any]]] = [("expression", expression)]
    for role in ("adjacency", "regulons", "aucell"):
        declaration = _require_fields(
            artifacts[role], {"file", "sha256"}, description=f"pySCENIC {role} artifact"
        )
        _sha256_text(declaration["sha256"], description=f"pySCENIC {role} SHA-256")
        declarations.append((role, declaration))
    resources = manifest["resources"]
    if not isinstance(resources, list) or not resources:
        raise ValueError("pySCENIC manifest resources must be a non-empty list.")
    resource_roles: list[str] = []
    normalized_resources: list[dict[str, Any]] = []
    for position, raw_resource in enumerate(resources):
        resource = _require_fields(
            raw_resource,
            {
                "role",
                "file",
                "sha256",
                "release",
                "organism",
                "genome_build",
                "gene_namespace",
                "license",
                "citation",
            },
            description=f"pySCENIC resource {position}",
        )
        role = _canonical_text(resource["role"], description=f"pySCENIC resource {position} role")
        if role not in {"tf_list", "ranking_database", "motif_annotations"}:
            raise ValueError(f"pySCENIC resource {position} has unsupported role {role!r}.")
        resource_roles.append(role)
        for field_name in ("release", "organism", "genome_build", "gene_namespace", "license", "citation"):
            _canonical_text(resource[field_name], description=f"pySCENIC resource {position} {field_name}")
        if any(resource[name] != manifest[name] for name in ("organism", "genome_build", "gene_namespace")):
            raise ValueError(f"pySCENIC resource {position} organism/build/namespace is inconsistent.")
        _sha256_text(resource["sha256"], description=f"pySCENIC resource {position} SHA-256")
        normalized_resources.append(copy.deepcopy(resource))
        declarations.append((f"resource:{position}:{role}", resource))
    if resource_roles.count("tf_list") != 1 or resource_roles.count("motif_annotations") != 1:
        raise ValueError("pySCENIC resources require exactly one TF list and one motif annotation table.")
    if resource_roles.count("ranking_database") < 1:
        raise ValueError("pySCENIC resources require at least one ranking database.")
    paths: dict[str, Path] = {}
    names: set[str] = set()
    for label, declaration in declarations:
        name = _canonical_text(declaration["file"], description=f"pySCENIC {label} file")
        if name in names:
            raise ValueError(f"pySCENIC bundle file is declared more than once: {name!r}.")
        names.add(name)
        paths[label] = _resolve_bundle_member(base, name, description=f"pySCENIC {label} file")
    if paths["expression"].suffix.lower() != ".csv":
        raise ValueError("pySCENIC expression bundle member must be an official cell-by-gene CSV.")
    if any(paths[role].suffix.lower() != ".csv" for role in ("adjacency", "regulons", "aucell")):
        raise ValueError("pySCENIC adjacency, regulon, and AUCell bundle members must be CSV files.")
    resource_files = {
        resource["role"]: [] for resource in normalized_resources
    }
    for position, resource in enumerate(normalized_resources):
        resource_files[resource["role"]].append(resource["file"])
        expected_suffix = {
            "tf_list": ".txt",
            "ranking_database": ".feather",
            "motif_annotations": ".tbl",
        }[resource["role"]]
        if paths[f"resource:{position}:{resource['role']}"] .suffix.lower() != expected_suffix:
            raise ValueError(
                f"pySCENIC {resource['role']} resource must use the official {expected_suffix} extension."
            )
    tf_file = resource_files["tf_list"][0]
    ranking_files = resource_files["ranking_database"]
    motif_file = resource_files["motif_annotations"][0]
    grn_command = normalized_commands["grn"]
    ctx_command = normalized_commands["ctx"]
    aucell_command = normalized_commands["aucell"]
    if grn_command[2:4] != [expression["file"], tf_file]:
        raise ValueError("pySCENIC grn positional inputs do not match the manifest expression and TF list.")
    if _command_flag_alias(grn_command, ("-o", "--output"), stage="grn") != artifacts["adjacency"]["file"]:
        raise ValueError("pySCENIC grn output does not match the declared adjacency artifact.")
    if ctx_command[2 : 3 + len(ranking_files)] != [artifacts["adjacency"]["file"], *ranking_files]:
        raise ValueError("pySCENIC ctx positional inputs do not match the adjacency and ranking resources.")
    if _command_flag(ctx_command, "--annotations_fname", stage="ctx") != motif_file:
        raise ValueError("pySCENIC ctx motif annotations do not match the declared resource.")
    if _command_flag(ctx_command, "--expression_mtx_fname", stage="ctx") != expression["file"]:
        raise ValueError("pySCENIC ctx expression input does not match the declared expression export.")
    if _command_flag_alias(ctx_command, ("-o", "--output"), stage="ctx") != artifacts["regulons"]["file"]:
        raise ValueError("pySCENIC ctx output does not match the declared regulon artifact.")
    if aucell_command[2:4] != [expression["file"], artifacts["regulons"]["file"]]:
        raise ValueError("pySCENIC aucell positional inputs do not match the expression and regulon artifacts.")
    if _command_flag_alias(aucell_command, ("-o", "--output"), stage="aucell") != artifacts["aucell"]["file"]:
        raise ValueError("pySCENIC aucell output does not match the declared AUCell artifact.")
    normalized = copy.deepcopy(manifest)
    normalized["commands"] = normalized_commands
    normalized["resources"] = normalized_resources
    return normalized, paths, {"cells": cells, "genes": genes, "commands": command_parameters}


def _load_manifest_and_verify(
    manifest_path: str | os.PathLike[str],
    *,
    max_file_bytes: int,
) -> tuple[dict[str, Any], dict[str, Path], dict[str, Any], dict[str, dict[str, Any]], str]:
    manifest_file = Path(manifest_path)
    if manifest_file.is_symlink():
        raise ValueError("pySCENIC manifest cannot be a symlink.")
    manifest_file = manifest_file.resolve(strict=True)
    manifest_hash, manifest_size = _sha256_file(manifest_file, max_file_bytes=max_file_bytes)
    try:
        text = manifest_file.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"Cannot read pySCENIC manifest as UTF-8: {exc}.") from exc
    post_manifest_hash, post_manifest_size = _sha256_file(
        manifest_file, max_file_bytes=max_file_bytes
    )
    if post_manifest_hash != manifest_hash or post_manifest_size != manifest_size:
        raise ValueError("pySCENIC manifest changed between hashing and parsing.")
    manifest = _strict_json_loads(text, description="pySCENIC manifest")
    normalized, paths, details = _validate_manifest(manifest, base=manifest_file.parent)
    declared_by_label: dict[str, str] = {
        "expression": normalized["expression"]["sha256"],
        **{
            role: normalized["artifacts"][role]["sha256"]
            for role in ("adjacency", "regulons", "aucell")
        },
    }
    for position, resource in enumerate(normalized["resources"]):
        declared_by_label[f"resource:{position}:{resource['role']}"] = resource["sha256"]
    accounting: dict[str, dict[str, Any]] = {}
    for label, path in paths.items():
        observed_hash, size = _sha256_file(path, max_file_bytes=max_file_bytes)
        if observed_hash != declared_by_label[label]:
            raise ValueError(
                f"pySCENIC bundle SHA-256 mismatch for {label}: declared={declared_by_label[label]}, "
                f"observed={observed_hash}."
            )
        accounting[label] = {
            "file": str(path.relative_to(manifest_file.parent.resolve())),
            "sha256": observed_hash,
            "bytes": size,
        }
    accounting["manifest"] = {
        "file": manifest_file.name,
        "sha256": manifest_hash,
        "bytes": manifest_size,
    }
    return normalized, paths, details, accounting, manifest_hash


def pyscenic_bundle_cache_fingerprint(
    manifest_path: str | os.PathLike[str],
    *,
    max_file_bytes: int = 2_147_483_647,
) -> tuple[str, tuple[tuple[str, str, int], ...]]:
    _positive_int(max_file_bytes, description="pySCENIC max_file_bytes")
    _, _, _, accounting, manifest_hash = _load_manifest_and_verify(
        manifest_path, max_file_bytes=max_file_bytes
    )
    members = tuple(
        (label, item["sha256"], item["bytes"])
        for label, item in sorted(accounting.items())
        if label != "manifest"
    )
    return manifest_hash, members


def _axis_sha256(values: Sequence[str]) -> str:
    return hashlib.sha256(_canonical_json(list(values)).encode("utf-8")).hexdigest()


def _new_matrix_digest(rows: int, columns: int) -> Any:
    digest = hashlib.sha256()
    digest.update(b"openbio-sparse-row-matrix-v1\0")
    digest.update(struct.pack("<QQ", rows, columns))
    return digest


def _update_matrix_digest(digest: Any, indices: Any, values: Any, *, numpy: Any) -> None:
    indices = numpy.asarray(indices, dtype="<u8")
    values = numpy.asarray(values, dtype="<f8")
    digest.update(struct.pack("<Q", int(indices.size)))
    digest.update(indices.tobytes(order="C"))
    digest.update(values.tobytes(order="C"))


def _matrix_sha256(matrix: Any, *, rows: int, columns: int, numpy: Any, sparse: Any) -> str:
    if getattr(matrix, "shape", None) != (rows, columns):
        raise ValueError("AnnData expression matrix shape does not match its axes.")
    digest = _new_matrix_digest(rows, columns)
    if sparse.issparse(matrix):
        work = matrix.tocsr(copy=True)
        work.sum_duplicates()
        work.sort_indices()
        work.eliminate_zeros()
        data = numpy.asarray(work.data, dtype=float)
        if not bool(numpy.isfinite(data).all()):
            raise ValueError("pySCENIC expression values must be finite.")
        for row in range(rows):
            start, end = int(work.indptr[row]), int(work.indptr[row + 1])
            _update_matrix_digest(
                digest,
                work.indices[start:end],
                data[start:end],
                numpy=numpy,
            )
    else:
        for row in range(rows):
            values = numpy.asarray(matrix[row, :], dtype=float).reshape(-1)
            if values.size != columns:
                raise ValueError("AnnData dense expression row width changed during fingerprinting.")
            if not bool(numpy.isfinite(values).all()):
                raise ValueError("pySCENIC expression values must be finite.")
            indices = numpy.flatnonzero(values != 0.0)
            _update_matrix_digest(digest, indices, values[indices], numpy=numpy)
    return digest.hexdigest()


def _csv_rows(path: Path, *, max_field_bytes: int) -> Any:
    previous_limit = csv.field_size_limit()
    next_limit = min(max(max_field_bytes, previous_limit), 2**31 - 1)
    with _CSV_LIMIT_LOCK:
        csv.field_size_limit(next_limit)
        try:
            try:
                with path.open("r", encoding="utf-8-sig", newline="") as stream:
                    yield from csv.reader(stream, strict=True)
            except (OSError, UnicodeError, csv.Error) as exc:
                raise ValueError(f"Cannot parse strict CSV bundle member {path.name!r}: {exc}.") from exc
        finally:
            csv.field_size_limit(previous_limit)


def _tsv_rows(path: Path, *, max_field_bytes: int) -> Any:
    previous_limit = csv.field_size_limit()
    next_limit = min(max(max_field_bytes, previous_limit), 2**31 - 1)
    with _CSV_LIMIT_LOCK:
        csv.field_size_limit(next_limit)
        try:
            try:
                with path.open("r", encoding="utf-8-sig", newline="") as stream:
                    yield from csv.reader(stream, delimiter="\t", strict=True)
            except (OSError, UnicodeError, csv.Error) as exc:
                raise ValueError(f"Cannot parse strict TSV bundle member {path.name!r}: {exc}.") from exc
        finally:
            csv.field_size_limit(previous_limit)


def _read_adjacency_csv(
    path: Path,
    *,
    genes: tuple[str, ...],
    max_adjacency_edges: int,
    max_file_bytes: int,
) -> dict[str, int]:
    rows = _csv_rows(path, max_field_bytes=max_file_bytes)
    try:
        header = next(rows)
    except StopIteration as exc:
        raise ValueError("pySCENIC adjacency CSV is empty.") from exc
    if header != ["TF", "target", "importance"]:
        raise ValueError("pySCENIC adjacency CSV columns must be exactly ['TF', 'target', 'importance'].")
    gene_set = set(genes)
    seen_edges: set[tuple[str, str]] = set()
    transcription_factors: set[str] = set()
    targets: set[str] = set()
    zero_importance_edges = 0
    for line_number, row in enumerate(rows, start=2):
        if not row or all(value == "" for value in row):
            raise ValueError(f"pySCENIC adjacency CSV contains a blank row at line {line_number}.")
        if len(row) != len(header):
            raise ValueError(f"pySCENIC adjacency CSV line {line_number} has an inconsistent width.")
        tf = _canonical_text(row[0], description="pySCENIC adjacency TF")
        target = _canonical_text(row[1], description="pySCENIC adjacency target")
        if tf not in gene_set or target not in gene_set:
            raise ValueError("pySCENIC adjacency references a TF or target outside the bound expression axis.")
        edge = (tf, target)
        if edge in seen_edges:
            raise ValueError(f"pySCENIC adjacency contains duplicate TF-target edge {edge!r}.")
        seen_edges.add(edge)
        if len(seen_edges) > max_adjacency_edges:
            raise ValueError("pySCENIC adjacency edges exceed max_adjacency_edges.")
        token = row[2]
        if not token or token != token.strip():
            raise ValueError("pySCENIC adjacency importance must be canonical and nonblank.")
        try:
            importance = float(token)
        except ValueError as exc:
            raise ValueError("pySCENIC adjacency importance must be numeric.") from exc
        if not math.isfinite(importance) or importance < 0.0:
            raise ValueError("pySCENIC adjacency importance must be finite and nonnegative.")
        if importance == 0.0:
            zero_importance_edges += 1
        transcription_factors.add(tf)
        targets.add(target)
    if not seen_edges:
        raise ValueError("pySCENIC adjacency CSV contains no TF-target edges.")
    return {
        "edges": len(seen_edges),
        "transcription_factors": len(transcription_factors),
        "targets": len(targets),
        "zero_importance_edges": zero_importance_edges,
    }


def _read_tf_list(path: Path, *, genes: tuple[str, ...]) -> dict[str, int]:
    try:
        stream = path.open("r", encoding="utf-8-sig", newline="")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"Cannot parse pySCENIC TF list as UTF-8: {exc}.") from exc
    names: list[str] = []
    try:
        with stream:
            for line_number, raw_line in enumerate(stream, start=1):
                line = raw_line.rstrip("\r\n")
                if not line or line.lstrip().startswith("#"):
                    continue
                names.append(
                    _canonical_text(line, description=f"pySCENIC TF-list line {line_number}")
                )
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"Cannot parse pySCENIC TF list as UTF-8: {exc}.") from exc
    if not names:
        raise ValueError("pySCENIC TF list contains no identifiers.")
    if len(names) != len(set(names)):
        raise ValueError("pySCENIC TF list contains duplicate identifiers.")
    gene_set = set(genes)
    overlap = sum(name in gene_set for name in names)
    if overlap < 1:
        raise ValueError("pySCENIC TF list has no identifier on the bound expression axis.")
    return {"identifiers": len(names), "identifiers_in_expression": overlap}


def _read_motif_annotations(path: Path, *, max_file_bytes: int) -> dict[str, int]:
    rows = _tsv_rows(path, max_field_bytes=max_file_bytes)
    try:
        header = next(rows)
    except StopIteration as exc:
        raise ValueError("pySCENIC motif-annotation table is empty.") from exc
    required = (
        "#motif_id",
        "gene_name",
        "motif_similarity_qvalue",
        "orthologous_identity",
        "description",
    )
    if any(header.count(name) != 1 for name in required):
        raise ValueError(f"pySCENIC motif annotations require unique columns {list(required)!r}.")
    positions = {name: header.index(name) for name in required}
    records = 0
    for line_number, row in enumerate(rows, start=2):
        if not row or all(value == "" for value in row):
            raise ValueError(f"pySCENIC motif annotations contain a blank row at line {line_number}.")
        if len(row) != len(header):
            raise ValueError(f"pySCENIC motif-annotation line {line_number} has an inconsistent width.")
        motif = _canonical_text(row[positions["#motif_id"]], description="pySCENIC annotation motif ID")
        gene = _canonical_text(row[positions["gene_name"]], description="pySCENIC annotation gene")
        _canonical_text(row[positions["description"]], description="pySCENIC annotation description")
        numeric: dict[str, float] = {}
        for column in ("motif_similarity_qvalue", "orthologous_identity"):
            token = row[positions[column]]
            if not token or token != token.strip():
                raise ValueError(f"pySCENIC annotation {column} must be canonical and nonblank.")
            try:
                numeric[column] = float(token)
            except ValueError as exc:
                raise ValueError(f"pySCENIC annotation {column} must be numeric.") from exc
            if not math.isfinite(numeric[column]):
                raise ValueError(f"pySCENIC annotation {column} must be finite.")
        if not 0.0 <= numeric["motif_similarity_qvalue"] <= 1.0:
            raise ValueError("pySCENIC annotation motif_similarity_qvalue must lie in [0, 1].")
        if numeric["orthologous_identity"] < 0.0:
            raise ValueError("pySCENIC annotation orthologous_identity cannot be negative.")
        del motif, gene
        records += 1
    if records < 1:
        raise ValueError("pySCENIC motif-annotation table contains no records.")
    return {"records": records}


def _read_ranking_database(path: Path) -> dict[str, str]:
    try:
        with path.open("rb") as stream:
            prefix = stream.read(6)
            stream.seek(-6, os.SEEK_END)
            suffix = stream.read(6)
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot inspect pySCENIC ranking database framing: {exc}.") from exc
    if prefix != b"ARROW1" or suffix != b"ARROW1":
        raise ValueError("pySCENIC ranking databases must use Feather-v2/Arrow IPC framing.")
    return {"format": "feather-v2"}


def _read_expression_csv(
    path: Path,
    *,
    expected_observations: tuple[str, ...],
    expected_genes: tuple[str, ...],
    numpy: Any,
    max_file_bytes: int,
) -> str:
    rows = _csv_rows(path, max_field_bytes=max_file_bytes)
    try:
        header = next(rows)
    except StopIteration as exc:
        raise ValueError("pySCENIC expression CSV is empty.") from exc
    if len(header) != len(expected_genes) + 1 or header[1:] != list(expected_genes):
        raise ValueError("pySCENIC expression CSV gene axis/order differs from AnnData.")
    if header[0] and header[0] != header[0].strip():
        raise ValueError("pySCENIC expression CSV index header must be canonical when present.")
    digest = _new_matrix_digest(len(expected_observations), len(expected_genes))
    seen_rows = 0
    for row_number, row in enumerate(rows, start=2):
        if not row or all(value == "" for value in row):
            raise ValueError(f"pySCENIC expression CSV contains a blank row at line {row_number}.")
        if len(row) != len(header):
            raise ValueError(f"pySCENIC expression CSV line {row_number} has an inconsistent width.")
        if seen_rows >= len(expected_observations):
            raise ValueError("pySCENIC expression CSV contains extra cells.")
        if row[0] != expected_observations[seen_rows]:
            raise ValueError("pySCENIC expression CSV cell axis/order differs from AnnData.")
        numeric: list[float] = []
        for column, token in enumerate(row[1:], start=1):
            if not token or token != token.strip():
                raise ValueError(
                    f"pySCENIC expression CSV has a blank/noncanonical value at line {row_number}, column {column + 1}."
                )
            try:
                value = float(token)
            except ValueError as exc:
                raise ValueError("pySCENIC expression CSV contains a nonnumeric value.") from exc
            if not math.isfinite(value):
                raise ValueError("pySCENIC expression CSV values must be finite.")
            numeric.append(value)
        array = numpy.asarray(numeric, dtype=float)
        indices = numpy.flatnonzero(array != 0.0)
        _update_matrix_digest(digest, indices, array[indices], numpy=numpy)
        seen_rows += 1
    if seen_rows != len(expected_observations):
        raise ValueError("pySCENIC expression CSV is missing cells.")
    return digest.hexdigest()


def _read_auc_csv(
    path: Path,
    *,
    expected_observations: tuple[str, ...],
    max_dense_bytes: int,
    max_file_bytes: int,
    numpy: Any,
    pandas: Any,
) -> DataFrame:
    rows = _csv_rows(path, max_field_bytes=max_file_bytes)
    try:
        header = next(rows)
    except StopIteration as exc:
        raise ValueError("pySCENIC AUCell CSV is empty.") from exc
    if len(header) < 2:
        raise ValueError("pySCENIC AUCell CSV must have a cell index and named regulon columns.")
    if header[0] and header[0] != header[0].strip():
        raise ValueError("pySCENIC AUCell CSV index header must be canonical when present.")
    regulons = _canonical_axis(header[1:], name="pySCENIC AUCell regulon axis")
    for regulon in regulons:
        match = _REGULON_PATTERN.fullmatch(regulon)
        if match is None or not match.group(1) or match.group(1) != match.group(1).strip():
            raise ValueError(f"pySCENIC AUCell regulon name is not canonical TF(+/-): {regulon!r}.")
    estimated_peak_dense_bytes = len(expected_observations) * len(regulons) * 96
    if estimated_peak_dense_bytes > max_dense_bytes:
        raise ValueError(
            f"pySCENIC AUCell import requires an estimated {estimated_peak_dense_bytes:,} dense working bytes, exceeding "
            f"max_dense_bytes={max_dense_bytes:,}."
        )
    values: list[list[float]] = []
    for row_number, row in enumerate(rows, start=2):
        if not row or all(value == "" for value in row):
            raise ValueError(f"pySCENIC AUCell CSV contains a blank row at line {row_number}.")
        if len(row) != len(header):
            raise ValueError(f"pySCENIC AUCell CSV line {row_number} has an inconsistent width.")
        position = len(values)
        if position >= len(expected_observations):
            raise ValueError("pySCENIC AUCell CSV contains extra cells.")
        if row[0] != expected_observations[position]:
            raise ValueError("pySCENIC AUCell CSV cell axis/order differs from AnnData.")
        parsed: list[float] = []
        for token in row[1:]:
            if not token or token != token.strip():
                raise ValueError("pySCENIC AUCell CSV contains a blank/noncanonical score.")
            try:
                value = float(token)
            except ValueError as exc:
                raise ValueError("pySCENIC AUCell CSV contains a nonnumeric score.") from exc
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError("pySCENIC AUCell scores must be finite and lie in [0, 1].")
            parsed.append(value)
        values.append(parsed)
    if len(values) != len(expected_observations):
        raise ValueError("pySCENIC AUCell CSV is missing cells.")
    return pandas.DataFrame(values, index=expected_observations, columns=regulons, dtype=float)


def _literal_depth(value: Any, *, limit: int, depth: int = 0) -> int:
    if depth > limit:
        raise ValueError("pySCENIC motif literal nesting exceeds the safety limit.")
    if isinstance(value, (list, tuple, set, frozenset, dict)):
        items = value.items() if isinstance(value, dict) else value
        for item in items:
            if isinstance(value, dict):
                key, nested = item
                _literal_depth(key, limit=limit, depth=depth + 1)
                _literal_depth(nested, limit=limit, depth=depth + 1)
            else:
                _literal_depth(item, limit=limit, depth=depth + 1)
    return depth


def _parse_context(text: str) -> tuple[str, ...]:
    if len(text.encode("utf-8")) > 65536:
        raise ValueError("pySCENIC Context literal exceeds 64 KiB.")
    prefix = "frozenset("
    if not text.startswith(prefix) or not text.endswith(")"):
        raise ValueError("pySCENIC Context must use the official frozenset(<set-literal>) grammar.")
    inner = text[len(prefix) : -1]
    if not inner:
        raise ValueError("pySCENIC Context cannot be empty.")
    try:
        parsed = ast.literal_eval(inner)
    except (SyntaxError, ValueError) as exc:
        raise ValueError("pySCENIC Context contains an invalid safe set literal.") from exc
    _literal_depth(parsed, limit=3)
    if not isinstance(parsed, set) or not parsed:
        raise ValueError("pySCENIC Context inner value must be a non-empty set literal.")
    result = tuple(sorted(_canonical_text(value, description="pySCENIC Context item") for value in parsed))
    signs = set(result).intersection({"activating", "repressing"})
    if len(signs) != 1:
        raise ValueError("pySCENIC Context must contain exactly one activating/repressing marker.")
    return result


def _parse_targets(text: str, *, remaining_edges: int) -> tuple[tuple[str, float], ...]:
    if len(text.encode("utf-8")) > 64 * 1024 * 1024:
        raise ValueError("pySCENIC TargetGenes literal exceeds 64 MiB.")
    try:
        parsed = ast.literal_eval(text)
    except (SyntaxError, ValueError) as exc:
        raise ValueError("pySCENIC TargetGenes contains an invalid safe literal.") from exc
    _literal_depth(parsed, limit=4)
    if not isinstance(parsed, list) or not parsed:
        raise ValueError("pySCENIC TargetGenes must be a non-empty list of (gene, weight) pairs.")
    if len(parsed) > remaining_edges:
        raise ValueError("pySCENIC motif target occurrences exceed max_regulon_edges.")
    result: list[tuple[str, float]] = []
    seen: set[str] = set()
    for item in parsed:
        if not isinstance(item, tuple) or len(item) != 2:
            raise ValueError("pySCENIC TargetGenes entries must be exact (gene, weight) tuples.")
        target = _canonical_text(item[0], description="pySCENIC target gene")
        weight = item[1]
        if isinstance(weight, bool) or not isinstance(weight, (int, float)):
            raise TypeError("pySCENIC target weight must be numeric.")
        weight = float(weight)
        if not math.isfinite(weight) or weight < 0.0:
            raise ValueError("pySCENIC target weights must be finite and nonnegative.")
        if target in seen:
            raise ValueError("pySCENIC TargetGenes contains a duplicate target within one motif row.")
        seen.add(target)
        result.append((target, weight))
    return tuple(result)


def _optional_float(token: str, *, description: str) -> float | None:
    if token == "":
        return None
    if token != token.strip():
        raise ValueError(f"{description} must be canonical.")
    try:
        value = float(token)
    except ValueError as exc:
        raise ValueError(f"{description} must be numeric or blank.") from exc
    if not math.isfinite(value):
        raise ValueError(f"{description} must be finite when present.")
    return value


def _read_regulon_csv(
    path: Path,
    *,
    genes: tuple[str, ...],
    regulon_order: tuple[str, ...],
    max_regulon_edges: int,
    max_file_bytes: int,
    pandas: Any,
) -> DataFrame:
    iterator = _csv_rows(path, max_field_bytes=max_file_bytes)
    try:
        top_header = next(iterator)
        lower_header = next(iterator)
        index_header = next(iterator)
    except StopIteration as exc:
        raise ValueError("pySCENIC regulon CSV is missing its three-row multi-index header.") from exc
    width = 2 + len(_MOTIF_COLUMNS)
    if (
        len(top_header) != width
        or top_header[:2] != ["", ""]
        or top_header[2:] != ["Enrichment"] * len(_MOTIF_COLUMNS)
        or lower_header != ["", "", *_MOTIF_COLUMNS]
        or index_header != ["TF", "MotifID", *([""] * len(_MOTIF_COLUMNS))]
    ):
        raise ValueError("pySCENIC regulon CSV multi-index header/schema is not the official 0.12.1 shape.")
    parsed_rows: list[dict[str, Any]] = []
    exact_rows: set[Any] = set()
    target_occurrences = 0
    for line_number, row in enumerate(iterator, start=4):
        if not row or all(value == "" for value in row):
            raise ValueError(f"pySCENIC regulon CSV contains a blank row at line {line_number}.")
        if len(row) != width:
            raise ValueError(f"pySCENIC regulon CSV line {line_number} has an inconsistent width.")
        tf = _canonical_text(row[0], description="pySCENIC motif TF")
        motif = _canonical_text(row[1], description="pySCENIC motif ID")
        auc_value = _optional_float(row[2], description="pySCENIC motif AUC")
        nes = _optional_float(row[3], description="pySCENIC motif NES")
        if auc_value is None or nes is None:
            raise ValueError("pySCENIC motif AUC and NES cannot be blank.")
        similarity = _optional_float(row[4], description="pySCENIC motif similarity q-value")
        orthology = _optional_float(row[5], description="pySCENIC orthologous identity")
        annotation = _canonical_text(row[6], description="pySCENIC motif annotation")
        context = _parse_context(row[7])
        targets = _parse_targets(row[8], remaining_edges=max_regulon_edges - target_occurrences)
        target_occurrences += len(targets)
        rank_at_max = _optional_float(row[9], description="pySCENIC RankAtMax")
        if rank_at_max is None or not rank_at_max.is_integer() or rank_at_max < 0:
            raise ValueError("pySCENIC RankAtMax must be a nonnegative integer.")
        canonical = (
            tf,
            motif,
            auc_value,
            nes,
            similarity,
            orthology,
            annotation,
            context,
            targets,
            int(rank_at_max),
        )
        if canonical in exact_rows:
            raise ValueError("pySCENIC regulon CSV contains an exact duplicate motif record.")
        exact_rows.add(canonical)
        regulation = "repressing" if "repressing" in context else "activating"
        regulon = f"{tf}{'(-)' if regulation == 'repressing' else '(+)'}"
        parsed_rows.append(
            {
                "regulon": regulon,
                "transcription_factor": tf,
                "regulation": regulation,
                "context_items": context,
                "motif_id": motif,
                "targets": targets,
            }
        )
    if not parsed_rows:
        raise ValueError("pySCENIC regulon CSV contains no motif-pruned records.")
    observed_regulons = {row["regulon"] for row in parsed_rows}
    if observed_regulons != set(regulon_order):
        raise ValueError("pySCENIC motif-pruned regulon family does not equal AUCell columns.")
    gene_set = set(genes)
    by_regulon: dict[str, dict[str, Any]] = {
        regulon: {"contexts": set(), "edges": {}} for regulon in regulon_order
    }
    for row in parsed_rows:
        group = by_regulon[row["regulon"]]
        group["contexts"].update(
            item
            for item in row["context_items"]
            if item not in {"activating", "repressing"} and not item.endswith(".png")
        )
        for target, weight in row["targets"]:
            if target not in gene_set:
                raise ValueError(
                    f"pySCENIC regulon {row['regulon']!r} references target outside the expression axis: {target!r}."
                )
            edge = group["edges"].setdefault(target, {"weight": weight, "motifs": set()})
            edge["weight"] = max(edge["weight"], weight)
            edge["motifs"].add(row["motif_id"])
    output_rows: list[dict[str, Any]] = []
    for regulon in regulon_order:
        match = _REGULON_PATTERN.fullmatch(regulon)
        assert match is not None
        tf, sign = match.groups()
        regulation = "activating" if sign == "+" else "repressing"
        group = by_regulon[regulon]
        contexts = sorted(group["contexts"])
        context = ";".join(contexts) if contexts else "motif_pruned"
        sorted_edges = sorted(
            group["edges"].items(), key=lambda item: (-item[1]["weight"], item[0])
        )
        previous_weight: float | None = None
        rank = 0
        for position, (target, evidence) in enumerate(sorted_edges, start=1):
            if evidence["weight"] != previous_weight:
                rank = position
                previous_weight = evidence["weight"]
            motifs = sorted(evidence["motifs"])
            output_rows.append(
                {
                    "regulon": regulon,
                    "transcription_factor": tf,
                    "regulation": regulation,
                    "context": context,
                    "target": target,
                    "target_weight": float(evidence["weight"]),
                    "motif_evidence_count": len(motifs),
                    "motif_ids": motifs,
                    "target_rank": rank,
                }
            )
    if len(output_rows) > max_regulon_edges:
        raise ValueError("Consolidated pySCENIC regulon edges exceed max_regulon_edges.")
    return pandas.DataFrame(output_rows, columns=SCENIC_MEMBERSHIP_COLUMNS)


def _software_versions(packages: Sequence[str], *, openbio_version: str) -> dict[str, str]:
    result = {"python": platform.python_version(), "openbio-singlecell": openbio_version}
    for package in packages:
        if package in result:
            continue
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = "not-installed"
    return result


def import_pyscenic_bundle(
    adata: AnnData,
    manifest_path: str | os.PathLike[str],
    *,
    overwrite: bool = False,
    max_file_bytes: int = 2_147_483_647,
    max_adjacency_edges: int = 10_000_000,
    max_regulon_edges: int = 2_000_000,
    max_dense_bytes: int = 1_073_741_824,
    openbio_version: str = PLUGIN_VERSION,
    _portable_artifact: bool = False,
) -> tuple[AnnData, Any, dict[str, Any]]:
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if not isinstance(overwrite, bool):
        raise TypeError("pySCENIC overwrite must be boolean.")
    _positive_int(max_file_bytes, description="pySCENIC max_file_bytes")
    _positive_int(max_adjacency_edges, description="pySCENIC max_adjacency_edges")
    _positive_int(max_regulon_edges, description="pySCENIC max_regulon_edges")
    _positive_int(max_dense_bytes, description="pySCENIC max_dense_bytes")
    if adata.n_obs < 1 or adata.n_vars < 1:
        raise ValueError("pySCENIC import requires a non-empty AnnData object.")
    output = adata
    if SCENIC_ACTIVITY_KEY in output.obsm and not overwrite:
        raise ValueError(
            f"AnnData obsm[{SCENIC_ACTIVITY_KEY!r}] already exists; enable overwrite only after reviewing provenance."
        )
    observations = _canonical_axis(output.obs_names.tolist(), name="AnnData observation axis")
    genes = _canonical_axis(output.var_names.tolist(), name="AnnData gene axis")
    manifest, paths, details, accounting, manifest_hash = _load_manifest_and_verify(
        manifest_path, max_file_bytes=max_file_bytes
    )
    if details["cells"] != len(observations) or details["genes"] != len(genes):
        raise ValueError("pySCENIC manifest expression dimensions differ from AnnData.")
    observation_hash = _axis_sha256(observations)
    gene_hash = _axis_sha256(genes)
    matrix_hash = _matrix_sha256(
        output.X,
        rows=len(observations),
        columns=len(genes),
        numpy=np,
        sparse=sparse,
    )
    expression = manifest["expression"]
    if observation_hash != expression["observation_ids_sha256"]:
        raise ValueError("pySCENIC manifest observation identity fingerprint differs from AnnData.")
    if gene_hash != expression["gene_ids_sha256"]:
        raise ValueError("pySCENIC manifest gene identity fingerprint differs from AnnData.")
    if matrix_hash != expression["matrix_sha256"]:
        raise ValueError("pySCENIC manifest expression-matrix fingerprint differs from AnnData.")
    csv_matrix_hash = _read_expression_csv(
        paths["expression"],
        expected_observations=observations,
        expected_genes=genes,
        numpy=np,
        max_file_bytes=max_file_bytes,
    )
    if csv_matrix_hash != matrix_hash:
        raise ValueError("pySCENIC expression CSV values differ from the bound AnnData expression matrix.")
    validation_accounting: dict[str, Any] = {
        "adjacency": _read_adjacency_csv(
            paths["adjacency"],
            genes=genes,
            max_adjacency_edges=max_adjacency_edges,
            max_file_bytes=max_file_bytes,
        )
    }
    ranking_validation: list[dict[str, str]] = []
    for position, resource in enumerate(manifest["resources"]):
        label = f"resource:{position}:{resource['role']}"
        if resource["role"] == "tf_list":
            validation_accounting["tf_list"] = _read_tf_list(paths[label], genes=genes)
        elif resource["role"] == "motif_annotations":
            validation_accounting["motif_annotations"] = _read_motif_annotations(
                paths[label], max_file_bytes=max_file_bytes
            )
        else:
            ranking_validation.append(_read_ranking_database(paths[label]))
    validation_accounting["ranking_databases"] = ranking_validation
    activity = _read_auc_csv(
        paths["aucell"],
        expected_observations=observations,
        max_dense_bytes=max_dense_bytes,
        max_file_bytes=max_file_bytes,
        numpy=np,
        pandas=pd,
    )
    membership = _read_regulon_csv(
        paths["regulons"],
        genes=genes,
        regulon_order=tuple(activity.columns.tolist()),
        max_regulon_edges=max_regulon_edges,
        max_file_bytes=max_file_bytes,
        pandas=pd,
    )
    validation_accounting["regulons"] = {
        "consolidated_edges": len(membership),
        "zero_weight_edges": int((membership["target_weight"] == 0.0).sum()),
    }
    for parsed_label, parsed_path in paths.items():
        post_hash, post_size = _sha256_file(parsed_path, max_file_bytes=max_file_bytes)
        if post_hash != accounting[parsed_label]["sha256"] or post_size != accounting[parsed_label]["bytes"]:
            raise ValueError(f"pySCENIC {parsed_label} changed between verification and parsing.")
    manifest_file = Path(manifest_path)
    if manifest_file.is_symlink():
        raise ValueError("pySCENIC manifest changed to a symlink during parsing.")
    post_manifest_hash, post_manifest_size = _sha256_file(
        manifest_file.resolve(strict=True), max_file_bytes=max_file_bytes
    )
    if (
        post_manifest_hash != accounting["manifest"]["sha256"]
        or post_manifest_size != accounting["manifest"]["bytes"]
    ):
        raise ValueError("pySCENIC manifest changed between verification and completed parsing.")
    provenance = {
        "artifact_type": SCENIC_RESULT_ARTIFACT_TYPE,
        "artifact_schema_version": SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
        "producer_node_id": SCENIC_RESULT_PRODUCER_NODE_ID,
        "producer_schema": SCENIC_RESULT_PRODUCER_SCHEMA,
        "manifest_schema": PYSCENIC_MANIFEST_SCHEMA,
        "manifest_sha256": manifest_hash,
        "run_id": manifest["run_id"],
        "organism": manifest["organism"],
        "genome_build": manifest["genome_build"],
        "gene_namespace": manifest["gene_namespace"],
        "expression_state": manifest["expression_state"],
        "input_dimensions": {"cells": len(observations), "genes": len(genes)},
        "expression_fingerprints": {
            "observation_ids_sha256": observation_hash,
            "gene_ids_sha256": gene_hash,
            "matrix_sha256": matrix_hash,
        },
        "container": copy.deepcopy(manifest["container"]),
        "declared_software_versions": copy.deepcopy(manifest["software_versions"]),
        "commands": copy.deepcopy(manifest["commands"]),
        "resources": copy.deepcopy(manifest["resources"]),
        "bundle_accounting": copy.deepcopy(accounting),
    }
    artifact = (
        SCENICResultArtifact._portable_from_owned(activity, membership, provenance)
        if _portable_artifact
        else SCENICResultArtifact(activity, membership, provenance)
    )
    artifact_fingerprint = (
        artifact["artifact_fingerprint_sha256"]
        if isinstance(artifact, dict)
        else artifact.artifact_fingerprint_sha256
    )
    output.obsm[SCENIC_ACTIVITY_KEY] = activity
    activity_values = activity.to_numpy(dtype=float, copy=False)
    leading_regulons = []
    for regulon in activity.columns[:10]:
        values = activity[regulon].to_numpy(dtype=float)
        leading_regulons.append(
            {
                "regulon": regulon,
                "mean_auc": float(values.mean()),
                "max_auc": float(values.max()),
                "targets": int((membership["regulon"] == regulon).sum()),
            }
        )
    parameters = {
        "manifest_path": str(Path(manifest_path).resolve()),
        "overwrite": overwrite,
        "activity_key": SCENIC_ACTIVITY_KEY,
        "max_file_bytes": max_file_bytes,
        "max_adjacency_edges": max_adjacency_edges,
        "max_regulon_edges": max_regulon_edges,
        "max_dense_bytes": max_dense_bytes,
        "commands": copy.deepcopy(details["commands"]),
    }
    warnings = [
        "External versions are declared-and-hash-bound provenance; the importer does not claim they are locally loaded.",
        "Structural and content binding cannot prove that the declared external computation was honestly executed.",
        "Only manifest-declared files consumed by this importer are hash-bound; unrelated files beside the manifest are ignored.",
    ]
    if manifest["expression_state"] != "post_qc_raw_counts":
        warnings.append(
            f"Expression state {manifest['expression_state']!r} is an explicit expert choice rather than the "
            "recommended post_qc_raw_counts state; transformed or negative values can change GRN and AUCell semantics."
        )
    if manifest["container"] is None:
        warnings.append(
            "No immutable container identity was declared; package versions and command arguments remain recorded, "
            "but the external execution environment is less reproducible."
        )
    for stage in ("grn", "ctx", "aucell"):
        command = manifest["commands"][stage]
        if "--num_workers" not in command:
            warnings.append(
                f"The {stage} command omitted --num_workers and used the pySCENIC 0.12.1 host-dependent default."
            )
    for stage in ("grn", "aucell"):
        if "--seed" not in manifest["commands"][stage]:
            warnings.append(
                f"The {stage} command omitted --seed; its stochastic step is not exactly reproducible from the manifest."
            )
    if not any(flag in manifest["commands"]["grn"] for flag in ("-m", "--method")):
        warnings.append("The grn command omitted --method and therefore used the 0.12.1 default grnboost2.")
    if "--mode" not in manifest["commands"]["ctx"]:
        warnings.append(
            "The ctx command omitted --mode and therefore used the 0.12.1 default custom_multiprocessing."
        )
    if "--auc_threshold" not in manifest["commands"]["aucell"]:
        warnings.append("The aucell command omitted --auc_threshold and therefore used the 0.12.1 default 0.05.")
    if validation_accounting["adjacency"]["zero_importance_edges"]:
        warnings.append(
            f"The GRN adjacency contains {validation_accounting['adjacency']['zero_importance_edges']:,} "
            "zero-importance edge(s); they were retained as backend-valid external output."
        )
    if validation_accounting["regulons"]["zero_weight_edges"]:
        warnings.append(
            f"The final regulon family contains {validation_accounting['regulons']['zero_weight_edges']:,} "
            "zero-weight target edge(s); they were retained as explicit backend output."
        )
    summary = {
        "schema_version": PYSCENIC_IMPORT_SUMMARY_SCHEMA,
        "node_id": SCENIC_RESULT_PRODUCER_NODE_ID,
        "status": "validated_external_pyscenic_regulatory_evidence",
        "methods": (
            "Imported an externally executed three-stage pySCENIC 0.12.1 GRN/ctx/AUCell bundle after strict "
            "manifest, streamed SHA-256, safe motif-literal, exact axis, numeric-expression fingerprint, and complete "
            "regulon-family validation. No pySCENIC code, shell command, or network access ran in this node."
        ),
        "results": (
            f"Validated {len(observations):,} cells, {len(genes):,} genes, {activity.shape[1]:,} motif-pruned "
            f"regulons, and {len(membership):,} consolidated regulon-target edges."
        ),
        "key_results": {
            "scientific_label": "Exploratory cell-level SCENIC regulatory-network evidence",
            "cells": len(observations),
            "genes": len(genes),
            "regulons": int(activity.shape[1]),
            "regulon_edges": len(membership),
            "transcription_factors": int(membership["transcription_factor"].nunique()),
            "motif_identifiers": int(
                len({motif for motifs in membership["motif_ids"] for motif in motifs})
            ),
            "activity_min": float(activity_values.min()),
            "activity_max": float(activity_values.max()),
            "leading_regulons": leading_regulons,
            "artifact_fingerprint_sha256": artifact_fingerprint,
            "manifest_sha256": manifest_hash,
            "bundle_accounting": copy.deepcopy(accounting),
            "validation_accounting": copy.deepcopy(validation_accounting),
            "resource_accounting": copy.deepcopy(manifest["resources"]),
            "external_provenance": {
                "container": copy.deepcopy(manifest["container"]),
                "software_versions": copy.deepcopy(manifest["software_versions"]),
                "run_id": manifest["run_id"],
                "organism": manifest["organism"],
                "genome_build": manifest["genome_build"],
                "gene_namespace": manifest["gene_namespace"],
                "expression_state": manifest["expression_state"],
            },
        },
        "parameters": parameters,
        "references": copy.deepcopy(PYSCENIC_REFERENCES),
        "software_versions": {
            **_software_versions(
                ["anndata", "numpy", "pandas", "scipy"], openbio_version=openbio_version
            ),
            **{
                f"external-declared:{name}": version
                for name, version in manifest["software_versions"].items()
            },
        },
        "warnings": warnings,
        "limitations": [
            "SCENIC combines expression-derived co-expression with motif priors and does not prove direct TF binding or causality.",
            "Cell-level AUCell activity is exploratory evidence, not a biological-Sample-level Condition contrast.",
            "Results inherit organism, build, namespace, resource-release, and input-cell-composition limitations from the external run.",
        ],
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return output, artifact, summary


def pyscenic_import_code(
    *,
    manifest_path: str,
    overwrite: bool,
    max_file_bytes: int,
    max_adjacency_edges: int,
    max_regulon_edges: int,
    max_dense_bytes: int,
    function_name: str = "import_validated_pyscenic_bundle",
) -> str:
    if not isinstance(function_name, str) or not function_name.isidentifier():
        raise ValueError("Generated pySCENIC importer function_name must be a Python identifier.")
    helpers = (
        _canonical_json,
        _canonical_axis,
        _activity_bytes,
        _canonical_membership,
        _strict_provenance,
        _result_fingerprint,
        SCENICResultArtifact,
        _strict_json_loads,
        _require_fields,
        _canonical_text,
        _sha256_text,
        _positive_int,
        _sha256_file,
        _resolve_bundle_member,
        _command_flag,
        _optional_command_flag,
        _optional_command_flag_alias,
        _command_flag_alias,
        _parse_command,
        _validate_manifest,
        _load_manifest_and_verify,
        _axis_sha256,
        _new_matrix_digest,
        _update_matrix_digest,
        _matrix_sha256,
        _csv_rows,
        _tsv_rows,
        _read_adjacency_csv,
        _read_tf_list,
        _read_motif_annotations,
        _read_ranking_database,
        _read_expression_csv,
        _read_auc_csv,
        _literal_depth,
        _parse_context,
        _parse_targets,
        _optional_float,
        _read_regulon_csv,
        _software_versions,
        import_pyscenic_bundle,
    )
    helper_source = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    return f'''from __future__ import annotations

import ast
import copy
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import struct
import threading
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

PLUGIN_VERSION = {PLUGIN_VERSION!r}
PYSCENIC_VERSION = {PYSCENIC_VERSION!r}
PYSCENIC_MANIFEST_SCHEMA = {PYSCENIC_MANIFEST_SCHEMA!r}
PYSCENIC_IMPORT_SUMMARY_SCHEMA = {PYSCENIC_IMPORT_SUMMARY_SCHEMA!r}
SCENIC_ACTIVITY_KEY = {SCENIC_ACTIVITY_KEY!r}
SCENIC_RESULT_ARTIFACT_TYPE = {SCENIC_RESULT_ARTIFACT_TYPE!r}
SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION = {SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION!r}
SCENIC_RESULT_PRODUCER_NODE_ID = {SCENIC_RESULT_PRODUCER_NODE_ID!r}
SCENIC_RESULT_PRODUCER_SCHEMA = {SCENIC_RESULT_PRODUCER_SCHEMA!r}
SCENIC_MEMBERSHIP_COLUMNS = {SCENIC_MEMBERSHIP_COLUMNS!r}
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{{64}}$")
_CONTAINER_DIGEST_PATTERN = re.compile(r"^sha256:[0-9a-f]{{64}}$")
_REGULON_PATTERN = re.compile(r"^(.+)\\(([+-])\\)$")
_CSV_LIMIT_LOCK = threading.Lock()
_MOTIF_COLUMNS = {_MOTIF_COLUMNS!r}
PYSCENIC_REFERENCES = {PYSCENIC_REFERENCES!r}

{helper_source}


def {function_name}(adata):
    """Import one content-bound external pySCENIC bundle without executing pySCENIC."""
    return import_pyscenic_bundle(
        adata,
        {str(Path(manifest_path).resolve())!r},
        overwrite={overwrite!r},
        max_file_bytes={max_file_bytes!r},
        max_adjacency_edges={max_adjacency_edges!r},
        max_regulon_edges={max_regulon_edges!r},
        max_dense_bytes={max_dense_bytes!r},
        openbio_version={PLUGIN_VERSION!r},
    )
'''


__all__ = [
    "PYSCENIC_IMPORT_SUMMARY_SCHEMA",
    "PYSCENIC_MANIFEST_SCHEMA",
    "PYSCENIC_VERSION",
    "SCENIC_ACTIVITY_KEY",
    "import_pyscenic_bundle",
    "pyscenic_bundle_cache_fingerprint",
    "pyscenic_import_code",
]
