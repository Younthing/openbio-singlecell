from __future__ import annotations

import shutil
import stat
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, fields, replace
from pathlib import Path
from typing import Any

from .artifact_codecs import read_anndata, write_anndata
from .artifact_envelope import write_strict_json
from .cnmf_standalone import (
    CNMFFixedPolicy,
    CNMFKMetric,
    CNMFResourceEstimate,
    CNMFRun,
    CNMFRunMetadata,
    _cleanup_preserving_primary,
    _is_reparse,
    _load_backend,
    _require_directory_identity,
    _validate_backend_layout,
    close_run_preserving_primary,
)
from .worker_protocol import read_json

CNMF_NATIVE_CODEC = "cnmf-native-directory"
CODEC_VERSION = 1
CONTROL_DIRECTORY = ".openbio-cnmf"
RUN_RECORD = "run.json"
BASE_DIRECTORY = "base"


def _empty_root(directory: str | Path) -> Path:
    root = Path(directory).absolute()
    _require_directory_identity(root, label="native artifact destination")
    if any(root.iterdir()):
        raise ValueError("cNMF artifact destination must be an empty directory.")
    return root


def _validate_tree(root: Path) -> None:
    _require_directory_identity(root, label="native artifact root")
    for child in root.iterdir():
        value = child.lstat()
        if _is_reparse(value):
            raise RuntimeError(f"cNMF native artifact path traverses a link or reparse point: {child}")
        if stat.S_ISDIR(value.st_mode):
            _validate_tree(child)
        elif not stat.S_ISREG(value.st_mode):
            raise RuntimeError(f"cNMF native artifact member must be a regular file or directory: {child}")


def _copy_tree(source: Path, target: Path, *, omit_control: bool = False) -> None:
    _require_directory_identity(source, label="native artifact source")
    for child in source.iterdir():
        if omit_control and child.name == CONTROL_DIRECTORY:
            continue
        value = child.lstat()
        if _is_reparse(value):
            raise RuntimeError(f"cNMF native artifact path traverses a link or reparse point: {child}")
        destination = target / child.name
        if stat.S_ISDIR(value.st_mode):
            destination.mkdir()
            _copy_tree(child, destination)
        elif stat.S_ISREG(value.st_mode):
            shutil.copy2(child, destination)
        else:
            raise RuntimeError(f"cNMF native artifact member must be a regular file or directory: {child}")


def write_cnmf_run(directory: str | Path, run: CNMFRun) -> None:
    if type(run) is not CNMFRun or run.closed:
        raise TypeError("Expected a live CNMFRun.")
    target = _empty_root(directory)
    source = Path(run.private_root).absolute()
    _copy_tree(source, target)
    control = target / CONTROL_DIRECTORY
    control.mkdir()
    base = control / BASE_DIRECTORY
    base.mkdir()
    write_anndata(base, object.__getattribute__(run, "_base_adata"))
    write_strict_json(
        control / RUN_RECORD,
        {
            "codec": CNMF_NATIVE_CODEC,
            "version": CODEC_VERSION,
            "metadata": asdict(run.metadata),
            "metrics": [asdict(metric) for metric in run.metrics],
        },
    )


def _exact_record(value: Any, record_type: type[Any], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {field.name for field in fields(record_type)}:
        raise ValueError(f"cNMF {label} record schema is invalid.")
    return value


def _metadata(value: Any) -> CNMFRunMetadata:
    data = _exact_record(value, CNMFRunMetadata, "metadata")
    return CNMFRunMetadata(
        **{
            **data,
            "input_advisories": tuple(data["input_advisories"]),
            "candidate_ks": tuple(data["candidate_ks"]),
            "realized_highvar_genes": tuple(data["realized_highvar_genes"]),
            "restart_seeds": tuple(tuple(item) for item in data["restart_seeds"]),
            "requested_resource": CNMFResourceEstimate(
                **_exact_record(data["requested_resource"], CNMFResourceEstimate, "requested resource")
            ),
            "realized_resource": CNMFResourceEstimate(
                **_exact_record(data["realized_resource"], CNMFResourceEstimate, "realized resource")
            ),
            "fixed_policy": CNMFFixedPolicy(**_exact_record(data["fixed_policy"], CNMFFixedPolicy, "fixed policy")),
            "artifact_hashes": tuple(tuple(item) for item in data["artifact_hashes"]),
            "known_upstream_warnings": tuple(tuple(item) for item in data["known_upstream_warnings"]),
        }
    )


def _records(root: Path) -> tuple[CNMFRunMetadata, tuple[CNMFKMetric, ...]]:
    value = read_json(root / CONTROL_DIRECTORY / RUN_RECORD)
    if (
        not isinstance(value, dict)
        or set(value) != {"codec", "version", "metadata", "metrics"}
        or value["codec"] != CNMF_NATIVE_CODEC
        or value["version"] != CODEC_VERSION
        or not isinstance(value["metrics"], list)
    ):
        raise ValueError("cNMF native artifact control record is invalid.")
    metrics = tuple(CNMFKMetric(**_exact_record(item, CNMFKMetric, "metric")) for item in value["metrics"])
    return _metadata(value["metadata"]), metrics


@contextmanager
def checkout_cnmf_run(
    directory: str | Path,
    *,
    checkout_parent: str | Path,
) -> Iterator[CNMFRun]:
    source = Path(directory).absolute()
    _validate_tree(source)
    metadata, metrics = _records(source)
    parent = Path(checkout_parent).absolute()
    _require_directory_identity(parent, label="checkout parent")
    owner = tempfile.TemporaryDirectory(prefix="openbio-cnmf-checkout-", dir=parent)
    checkout = Path(owner.name).resolve(strict=True)
    try:
        _copy_tree(source, checkout, omit_control=True)
        base_adata = read_anndata(source / CONTROL_DIRECTORY / BASE_DIRECTORY)
        api = _load_backend()
        backend = api.cnmf_class(output_dir=str(checkout), name="run")
        metadata = replace(
            metadata,
            backend_paths_fingerprint=_validate_backend_layout(backend, checkout, "run"),
        )
        run = CNMFRun(
            metadata=metadata,
            metrics=metrics,
            backend=backend,
            base_adata=base_adata,
            temporary_directory=owner,
        )
    except BaseException as primary:
        _cleanup_preserving_primary(owner.cleanup, primary, context="cNMF native checkout construction")
        raise
    try:
        yield run
    except BaseException as primary:
        close_run_preserving_primary(run, primary, context="cNMF native checkout")
        raise
    else:
        run.close()


__all__ = ["CNMF_NATIVE_CODEC", "checkout_cnmf_run", "write_cnmf_run"]
