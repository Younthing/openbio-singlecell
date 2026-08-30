from __future__ import annotations

import gc
import hashlib
import json
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from openbio_singlecell import artifact_persist as artifact_persist_module
from openbio_singlecell.artifact_persist import persist_artifact
from openbio_singlecell.artifact_runtime import ArtifactRuntime


def _runtime(tmp_path: Path, session: str = "session", run: str = "run") -> ArtifactRuntime:
    identifiers = iter((session, run))
    return ArtifactRuntime(tmp_path / session, uuid_factory=lambda: next(identifiers))


def test_persist_directory_streams_complete_copy_with_hash_manifest_and_keeps_source(tmp_path):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    source = lease.staging_path / "portable"
    (source / "nested").mkdir(parents=True)
    (source / "data.jsonl").write_bytes(b'{"value":1}\n')
    (source / "nested" / "schema.json").write_bytes(b'{"version":1}\n')
    ticket = lease.publish(
        {
            "result": {
                "kind": "OPENBIO_TABLE",
                "codec": "table-jsonl",
                "payload": "portable",
            }
        }
    )["result"]
    published_source = runtime.resolve(ticket)
    output_root = tmp_path / "output"
    output_root.mkdir()

    persisted = persist_artifact(runtime, ticket, output_root, "report")

    expected_path = output_root / "openbio-singlecell" / "artifacts" / "report__run-result"
    assert persisted.path == expected_path
    assert (expected_path / "data.jsonl").read_bytes() == b'{"value":1}\n'
    assert (expected_path / "nested" / "schema.json").read_bytes() == b'{"version":1}\n'
    assert published_source.is_dir()
    assert (published_source / "data.jsonl").is_file()
    manifest = json.loads((expected_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest == {
        "version": 1,
        "artifact": {
            "session": "session",
            "run": "run",
            "output": "result",
            "kind": "OPENBIO_TABLE",
            "codec": "table-jsonl",
        },
        "payload_type": "directory",
        "files": [
            {
                "path": "data.jsonl",
                "size": 12,
                "sha256": hashlib.sha256(b'{"value":1}\n').hexdigest(),
            },
            {
                "path": "nested/schema.json",
                "size": 14,
                "sha256": hashlib.sha256(b'{"version":1}\n').hexdigest(),
            },
        ],
    }
    assert not list((output_root / "openbio-singlecell" / "artifacts").glob("*.partial"))

    ready_path = lease.ready_path
    del ticket, lease
    gc.collect()
    assert not ready_path.exists()
    assert (expected_path / "data.jsonl").read_bytes() == b'{"value":1}\n'


@pytest.mark.parametrize("codec", ["scvi-native-directory", "cnmf-native-directory"])
def test_persist_rejects_session_native_artifacts_after_ticket_validation(tmp_path, codec):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    source = lease.staging_path / "native"
    source.mkdir()
    (source / "weights.bin").write_bytes(b"native")
    ticket = lease.publish(
        {
            "model": {
                "kind": "OPENBIO_NATIVE_MODEL",
                "codec": codec,
                "payload": "native",
                "worker_affinity": "python:probe",
            }
        }
    )["model"]
    published_source = runtime.resolve(ticket)
    output_root = tmp_path / "output"
    output_root.mkdir()

    with pytest.raises(ValueError, match="session-only native"):
        persist_artifact(runtime, ticket, output_root, "model")

    assert (published_source / "weights.bin").read_bytes() == b"native"
    assert not (output_root / "openbio-singlecell" / "artifacts" / "model__run-model").exists()


def test_persist_retry_of_the_same_unchanged_artifact_is_idempotent(tmp_path):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    (lease.staging_path / "data.h5ad").write_bytes(b"portable-h5ad")
    ticket = lease.publish(
        {
            "adata": {
                "kind": "OPENBIO_ANNDATA",
                "codec": "h5ad",
                "payload": "data.h5ad",
            }
        }
    )["adata"]
    output_root = tmp_path / "output"
    output_root.mkdir()

    first = persist_artifact(runtime, ticket, output_root, "analysis")
    manifest_before = (first.path / "manifest.json").read_bytes()
    second = persist_artifact(runtime, ticket, output_root, "analysis")

    assert second == first
    assert (second.path / "data.h5ad").read_bytes() == b"portable-h5ad"
    assert (second.path / "manifest.json").read_bytes() == manifest_before


def test_persist_fails_without_changing_an_existing_target_with_conflicting_identity(tmp_path):
    first_runtime = _runtime(tmp_path, session="session-a", run="shared-run")
    first_lease = first_runtime.begin_run()
    (first_lease.staging_path / "data.bin").write_bytes(b"first")
    first_ticket = first_lease.publish(
        {"result": {"kind": "OPENBIO_TABLE", "codec": "table-jsonl", "payload": "data.bin"}}
    )["result"]
    second_runtime = _runtime(tmp_path, session="session-b", run="shared-run")
    second_lease = second_runtime.begin_run()
    (second_lease.staging_path / "data.bin").write_bytes(b"second")
    second_ticket = second_lease.publish(
        {"result": {"kind": "OPENBIO_TABLE", "codec": "table-jsonl", "payload": "data.bin"}}
    )["result"]
    output_root = tmp_path / "output"
    output_root.mkdir()

    first = persist_artifact(first_runtime, first_ticket, output_root, "shared")
    original_manifest = (first.path / "manifest.json").read_bytes()

    with pytest.raises(FileExistsError, match="conflicts"):
        persist_artifact(second_runtime, second_ticket, output_root, "shared")

    assert (first.path / "data.bin").read_bytes() == b"first"
    assert (first.path / "manifest.json").read_bytes() == original_manifest
    assert not list(first.path.parent.glob("*.partial"))


def test_concurrent_persist_of_the_same_artifact_commits_once_and_is_idempotent(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    (lease.staging_path / "data.bin").write_bytes(b"shared")
    ticket = lease.publish(
        {"result": {"kind": "OPENBIO_TABLE", "codec": "table-jsonl", "payload": "data.bin"}}
    )["result"]
    output_root = tmp_path / "output"
    output_root.mkdir()
    barrier = threading.Barrier(2)
    original_rename = artifact_persist_module.os.rename

    def simultaneous_commit(source, target):
        if Path(source).name.endswith(".partial"):
            barrier.wait(timeout=5)
        return original_rename(source, target)

    monkeypatch.setattr(artifact_persist_module.os, "rename", simultaneous_commit)

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(
            executor.map(
                lambda _: persist_artifact(runtime, ticket, output_root, "shared"),
                range(2),
            )
        )

    assert results[0] == results[1]
    assert (results[0].path / "data.bin").read_bytes() == b"shared"
    assert not list(results[0].path.parent.glob("*.partial"))


def test_copy_failure_cleans_staging_without_publishing_a_visible_result(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    (lease.staging_path / "data.bin").write_bytes(b"source")
    ticket = lease.publish(
        {"result": {"kind": "OPENBIO_TABLE", "codec": "table-jsonl", "payload": "data.bin"}}
    )["result"]
    source = runtime.resolve(ticket)
    output_root = tmp_path / "output"
    output_root.mkdir()
    original_open = Path.open

    def fail_source_read(path, *args, **kwargs):
        if path == source and args and args[0] == "rb":
            raise OSError("simulated source read failure")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", fail_source_read)

    with pytest.raises(OSError, match="simulated"):
        persist_artifact(runtime, ticket, output_root, "failed")

    artifacts_root = output_root / "openbio-singlecell" / "artifacts"
    assert not (artifacts_root / "failed__run-result").exists()
    assert not list(artifacts_root.glob("*.partial"))


def test_staging_cleanup_failure_is_reported_instead_of_silently_leaking(tmp_path, monkeypatch):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    (lease.staging_path / "data.bin").write_bytes(b"source")
    ticket = lease.publish(
        {"result": {"kind": "OPENBIO_TABLE", "codec": "table-jsonl", "payload": "data.bin"}}
    )["result"]
    output_root = tmp_path / "output"
    output_root.mkdir()

    monkeypatch.setattr(
        artifact_persist_module,
        "_copy_file",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("copy failed")),
    )
    original_rmtree = artifact_persist_module.shutil.rmtree

    def locked_rmtree(path, *args, **kwargs):
        if Path(path).name.endswith(".partial"):
            if kwargs.get("ignore_errors"):
                return None
            raise OSError("staging locked")
        return original_rmtree(path, *args, **kwargs)

    monkeypatch.setattr(artifact_persist_module.shutil, "rmtree", locked_rmtree)

    with pytest.raises(OSError, match="staging locked"):
        persist_artifact(runtime, ticket, output_root, "locked")


@pytest.mark.parametrize("name", ["../escape", r"..\escape", "bad:name", "CON", "trailing."])
def test_persist_rejects_names_that_are_unsafe_on_windows_or_posix(tmp_path, name):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    (lease.staging_path / "data.bin").write_bytes(b"source")
    ticket = lease.publish(
        {"result": {"kind": "OPENBIO_TABLE", "codec": "table-jsonl", "payload": "data.bin"}}
    )["result"]
    output_root = tmp_path / "output"
    output_root.mkdir()

    with pytest.raises(ValueError, match="safe path component"):
        persist_artifact(runtime, ticket, output_root, name)

    assert not (tmp_path / "escape").exists()


def test_persist_keeps_nested_payload_files_named_manifest_json(tmp_path):
    runtime = _runtime(tmp_path)
    lease = runtime.begin_run()
    source = lease.staging_path / "portable" / "nested"
    source.mkdir(parents=True)
    (source / "manifest.json").write_bytes(b'{"payload":true}\n')
    ticket = lease.publish(
        {"result": {"kind": "OPENBIO_TABLE", "codec": "table-jsonl", "payload": "portable"}}
    )["result"]
    output_root = tmp_path / "output"
    output_root.mkdir()

    persisted = persist_artifact(runtime, ticket, output_root, "nested")

    assert (persisted.path / "nested" / "manifest.json").read_bytes() == b'{"payload":true}\n'
