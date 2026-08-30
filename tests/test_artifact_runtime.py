from __future__ import annotations

import copy
import gc
import json
import os
import pickle
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from openbio_singlecell import artifact_runtime as artifact_runtime_module
from openbio_singlecell.artifact_runtime import NATIVE_AFFINITY_CODECS, ArtifactRuntime, ArtifactTicket


def test_published_artifact_ticket_resolves_without_copying_or_pickling(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "data.bin").write_bytes(b"abc")

    ticket = run.publish(
        {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"}}
    )["result"]

    assert runtime.resolve(ticket).read_bytes() == b"abc"
    assert copy.copy(ticket) is ticket
    assert copy.deepcopy(ticket) is ticket
    with pytest.raises(TypeError, match="cannot be pickled"):
        pickle.dumps(ticket)


def test_begin_run_retries_uuid_collision_with_exclusive_directory(tmp_path: Path):
    identifiers = iter(("session", "taken", "fresh"))
    runtime = ArtifactRuntime(tmp_path, uuid_factory=lambda: next(identifiers))
    (runtime.session_path / "taken.partial").mkdir()

    run = runtime.begin_run()

    assert run.run_id == "fresh"
    assert run.staging_path.is_dir()
    assert not (runtime.session_path / "fresh.ready").exists()


def test_concurrent_sessions_never_share_a_directory(tmp_path: Path):
    identifiers = iter(("collision", "collision", "second"))
    lock = threading.Lock()

    def next_identifier():
        with lock:
            return next(identifiers)

    with ThreadPoolExecutor(max_workers=2) as executor:
        runtimes = list(executor.map(lambda _: ArtifactRuntime(tmp_path, uuid_factory=next_identifier), range(2)))

    assert len({runtime.session_path for runtime in runtimes}) == 2
    assert all(runtime.session_path.is_dir() for runtime in runtimes)


def test_concurrent_runs_never_share_a_staging_directory(tmp_path: Path):
    identifiers = iter(("session", "collision", "collision", "second"))
    lock = threading.Lock()

    def next_identifier():
        with lock:
            return next(identifiers)

    runtime = ArtifactRuntime(tmp_path, uuid_factory=next_identifier)
    with ThreadPoolExecutor(max_workers=2) as executor:
        runs = list(executor.map(lambda _: runtime.begin_run(), range(2)))

    assert len({run.staging_path for run in runs}) == 2
    assert all(run.staging_path.is_dir() for run in runs)


def test_concurrent_run_id_stays_reserved_after_partial_is_published(
    tmp_path: Path,
    monkeypatch,
):
    identifiers = iter(("session", "collision", "collision", "fresh"))
    identifier_lock = threading.Lock()

    def next_identifier():
        with identifier_lock:
            return next(identifiers)

    runtime = ArtifactRuntime(tmp_path, uuid_factory=next_identifier)
    collision_barrier = threading.Barrier(2)
    published = threading.Event()
    publisher_lock = threading.Lock()
    publisher_thread: int | None = None
    real_mkdir = Path.mkdir

    def controlled_mkdir(path: Path, *args, **kwargs):
        nonlocal publisher_thread
        if path.name != "collision.partial":
            return real_mkdir(path, *args, **kwargs)
        collision_barrier.wait(timeout=5)
        with publisher_lock:
            if publisher_thread is None:
                publisher_thread = threading.get_ident()
            is_publisher = publisher_thread == threading.get_ident()
        if not is_publisher:
            assert published.wait(timeout=5)
        return real_mkdir(path, *args, **kwargs)

    monkeypatch.setattr(Path, "mkdir", controlled_mkdir)

    def begin_and_publish():
        run = runtime.begin_run()
        if threading.get_ident() != publisher_thread:
            return run, None
        (run.staging_path / "data.bin").write_bytes(b"published")
        ticket = run.publish(
            {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"}}
        )["result"]
        published.set()
        return run, ticket

    with ThreadPoolExecutor(max_workers=2) as executor:
        runs_and_tickets = list(executor.map(lambda _: begin_and_publish(), range(2)))

    assert {run.run_id for run, _ in runs_and_tickets} == {"collision", "fresh"}
    assert not (runtime.session_path / "collision.partial").exists()
    published_ticket = next(ticket for _, ticket in runs_and_tickets if ticket is not None)
    assert runtime.resolve(published_ticket).read_bytes() == b"published"
    next(run for run, ticket in runs_and_tickets if ticket is None).abort()


def test_abort_discards_partial_run_without_publishing_it(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "unfinished.bin").write_bytes(b"unfinished")

    run.abort()

    assert not run.staging_path.exists()
    assert not run.ready_path.exists()
    run.abort()
    with pytest.raises(RuntimeError, match="aborted"):
        run.publish(
            {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "unfinished.bin"}}
        )


def test_one_run_publishes_file_and_directory_outputs_with_size_manifest(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "data.bin").write_bytes(b"abc")
    model = run.staging_path / "model"
    (model / "nested").mkdir(parents=True)
    (model / "weights.bin").write_bytes(b"12345")
    (model / "nested" / "metadata.json").write_bytes(b"{}")

    tickets = run.publish(
        {
            "data": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"},
            "model": {"kind": "OPENBIO_MODEL", "codec": "native-directory", "payload": "model"},
        }
    )

    assert runtime.resolve(tickets["data"]).name == "data.bin"
    assert runtime.resolve(tickets["model"]).name == "model"
    manifest = json.loads((run.ready_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["outputs"]["model"]["files"] == [
        {"path": "model/nested/metadata.json", "size": 2},
        {"path": "model/weights.bin", "size": 5},
    ]


def test_resolve_rejects_manifest_fields_outside_the_strict_schema(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "data.bin").write_bytes(b"abc")
    ticket = run.publish(
        {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"}}
    )["result"]
    manifest_path = run.ready_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"]["result"]["unexpected"] = True
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="manifest schema"):
        runtime.resolve(ticket)


def test_shared_run_is_removed_after_its_last_ticket_and_lease_are_released(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "a.bin").write_bytes(b"a")
    (run.staging_path / "b.bin").write_bytes(b"b")
    tickets = run.publish(
        {
            "a": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "a.bin"},
            "b": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "b.bin"},
        }
    )
    ready_path = run.ready_path

    del run
    first = tickets.pop("a")
    del first
    gc.collect()
    assert ready_path.is_dir()

    last = tickets.pop("b")
    del last
    gc.collect()
    assert not ready_path.exists()


def test_startup_scavenges_only_sessions_with_confirmed_dead_or_reused_pid(tmp_path: Path):
    root = tmp_path / "openbio-singlecell" / "artifact-runtime"
    root.mkdir(parents=True)
    for session, pid, created in (("dead", 11, 1.0), ("alive", 22, 2.0), ("reused", 33, 3.0)):
        session_path = root / f"{session}.session"
        (session_path / "old.ready").mkdir(parents=True)
        (session_path / "session.json").write_text(
            json.dumps({"version": 1, "pid": pid, "process_create_time": created}),
            encoding="utf-8",
        )

    observed_create_times = {11: None, 22: 2.0, 33: 30.0}
    runtime = ArtifactRuntime(
        tmp_path,
        uuid_factory=lambda: "current",
        process_identity=(44, 4.0),
        process_probe=observed_create_times.get,
    )

    assert not (root / "dead.session").exists()
    assert (root / "alive.session").is_dir()
    assert not (root / "reused.session").exists()
    assert json.loads((runtime.session_path / "session.json").read_text(encoding="utf-8")) == {
        "version": 1,
        "pid": 44,
        "process_create_time": 4.0,
    }


def test_failed_orphan_deletion_is_retried_on_the_next_startup(tmp_path: Path, monkeypatch):
    root = tmp_path / "openbio-singlecell" / "artifact-runtime"
    victim = root / "dead.session"
    victim.mkdir(parents=True)
    metadata = {"version": 1, "pid": 11, "process_create_time": 1.0}
    (victim / "session.json").write_text(json.dumps(metadata), encoding="utf-8")
    (victim / "locked.bin").write_bytes(b"locked")
    real_rmtree = artifact_runtime_module.shutil.rmtree
    failed = False

    def fail_once(path):
        nonlocal failed
        if Path(path) == victim and not failed:
            failed = True
            (victim / "session.json").unlink()
            raise PermissionError("still open")
        real_rmtree(path)

    monkeypatch.setattr(artifact_runtime_module.shutil, "rmtree", fail_once)
    process_times = {11: None, 44: 4.0, 55: 5.0}
    ArtifactRuntime(
        tmp_path,
        uuid_factory=lambda: "first",
        process_identity=(44, 4.0),
        process_probe=process_times.get,
    )
    assert json.loads((victim / "session.json").read_text(encoding="utf-8")) == metadata

    monkeypatch.setattr(artifact_runtime_module.shutil, "rmtree", real_rmtree)
    ArtifactRuntime(
        tmp_path,
        uuid_factory=lambda: "second",
        process_identity=(55, 5.0),
        process_probe=process_times.get,
    )
    assert not victim.exists()


def test_empty_run_cannot_publish_a_visible_manifest(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()

    with pytest.raises(ValueError, match="at least one output"):
        run.publish({})

    assert run.staging_path.is_dir()
    assert not run.ready_path.exists()


def test_resolve_rejects_ticket_without_its_process_local_lease(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    forged = ArtifactTicket(
        session=runtime.session_id,
        run="run",
        output="result",
        kind="OPENBIO_TEST",
        codec="binary",
    )

    with pytest.raises(ValueError, match="live run lease"):
        runtime.resolve(forged)


def test_resolve_rejects_changed_payload_size(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    payload = run.staging_path / "data.bin"
    payload.write_bytes(b"abc")
    ticket = run.publish(
        {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"}}
    )["result"]
    runtime.resolve(ticket).write_bytes(b"changed")

    with pytest.raises(ValueError, match="size"):
        runtime.resolve(ticket)


def test_resolve_rejects_manifest_path_escape(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "data.bin").write_bytes(b"abc")
    ticket = run.publish(
        {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"}}
    )["result"]
    manifest_path = run.ready_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"]["result"]["payload"] = "../outside.bin"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="stay inside"):
        runtime.resolve(ticket)


def test_resolve_rejects_replaced_ready_directory_symlink(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "data.bin").write_bytes(b"abc")
    ticket = run.publish(
        {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"}}
    )["result"]
    moved = run.ready_path.with_name(f"{run.run_id}.moved")
    run.ready_path.rename(moved)
    try:
        os.symlink(moved, run.ready_path, target_is_directory=True)
    except OSError as error:
        moved.rename(run.ready_path)
        pytest.skip(f"Directory symlinks are unavailable: {error}")

    with pytest.raises(ValueError, match="ready directory"):
        runtime.resolve(ticket)


def test_publish_rejects_windows_reparse_payload(tmp_path: Path, monkeypatch):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    payload = run.staging_path / "payload"
    payload.mkdir()
    (payload / "data.bin").write_bytes(b"outside")
    original_lstat = Path.lstat

    def reparse_lstat(path, *args, **kwargs):
        value = original_lstat(path, *args, **kwargs)
        if path == payload:
            return value.__class__(
                (*value[:],),
                {"st_file_attributes": 0x400},
            )
        return value

    monkeypatch.setattr(Path, "lstat", reparse_lstat)

    with pytest.raises(ValueError, match="link or reparse"):
        run.publish(
            {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "payload"}}
        )


def test_native_directory_manifest_can_record_worker_affinity(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    model = run.staging_path / "model"
    model.mkdir()
    (model / "weights.bin").write_bytes(b"weights")

    ticket = run.publish(
        {
            "model": {
                "kind": "OPENBIO_SCVI_MODEL",
                "codec": "scvi-native-directory",
                "payload": "model",
                "worker_affinity": "python:probe-fingerprint",
            }
        }
    )["model"]

    manifest = json.loads((run.ready_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["outputs"][ticket.output]["worker_affinity"] == "python:probe-fingerprint"
    assert runtime.resolve(ticket) == run.ready_path / "model"


def test_portable_codec_rejects_worker_affinity(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "data.h5ad").write_bytes(b"h5ad")

    with pytest.raises(ValueError, match="only for native scVI and cNMF"):
        run.publish(
            {
                "data": {
                    "kind": "OPENBIO_ANNDATA",
                    "codec": "h5ad",
                    "payload": "data.h5ad",
                    "worker_affinity": "python:probe-fingerprint",
                }
            }
        )

    assert not run.ready_path.exists()


@pytest.mark.parametrize("codec", sorted(NATIVE_AFFINITY_CODECS))
def test_native_codec_requires_worker_affinity_when_published(tmp_path: Path, codec: str):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    model = run.staging_path / "model"
    model.mkdir()
    (model / "weights.bin").write_bytes(b"weights")

    with pytest.raises(ValueError, match="require.*Worker affinity"):
        run.publish(
            {
                "model": {
                    "kind": "OPENBIO_NATIVE_MODEL",
                    "codec": codec,
                    "payload": "model",
                }
            }
        )


@pytest.mark.parametrize("codec", sorted(NATIVE_AFFINITY_CODECS))
def test_native_codec_manifest_without_worker_affinity_is_rejected(tmp_path: Path, codec: str):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    model = run.staging_path / "model"
    model.mkdir()
    (model / "weights.bin").write_bytes(b"weights")
    ticket = run.publish(
        {
            "model": {
                "kind": "OPENBIO_NATIVE_MODEL",
                "codec": codec,
                "payload": "model",
                "worker_affinity": "python:probe-fingerprint",
            }
        }
    )["model"]
    manifest_path = run.ready_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    del manifest["outputs"]["model"]["worker_affinity"]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="worker affinity"):
        runtime.resolve(ticket)


def test_record_returns_a_defensive_copy_after_full_ticket_validation(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    model = run.staging_path / "model"
    model.mkdir()
    (model / "weights.bin").write_bytes(b"weights")
    ticket = run.publish(
        {
            "model": {
                "kind": "OPENBIO_SCVI_MODEL",
                "codec": "scvi-native-directory",
                "payload": "model",
                "worker_affinity": "python:probe-fingerprint",
            }
        }
    )["model"]

    record = runtime.record(ticket)
    assert record["worker_affinity"] == "python:probe-fingerprint"
    record["files"][0]["size"] = 999

    assert runtime.record(ticket)["files"] == [{"path": "model/weights.bin", "size": 7}]


def test_close_waits_for_live_ticket_and_lease_before_removing_session(tmp_path: Path):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    (run.staging_path / "data.bin").write_bytes(b"abc")
    ticket = run.publish(
        {"result": {"kind": "OPENBIO_TEST", "codec": "binary", "payload": "data.bin"}}
    )["result"]
    ready_path = run.ready_path
    session_path = runtime.session_path

    runtime.close()
    assert session_path.is_dir()
    assert ready_path.is_dir()
    assert runtime.resolve(ticket).read_bytes() == b"abc"

    del run
    runtime.close()
    assert session_path.is_dir()
    assert ready_path.is_dir()

    del ticket
    gc.collect()
    assert not ready_path.exists()
    runtime.close()
    assert not session_path.exists()
    runtime.close()


def test_close_is_best_effort_and_preserves_session_identity_for_retry(tmp_path: Path, monkeypatch):
    runtime = ArtifactRuntime(tmp_path)
    run = runtime.begin_run()
    session_path = runtime.session_path
    expected_metadata = json.loads((session_path / "session.json").read_text(encoding="utf-8"))

    runtime.close()
    assert run.staging_path.is_dir()

    run.abort()
    del run
    gc.collect()
    real_rmtree = artifact_runtime_module.shutil.rmtree

    def fail_after_metadata_removal(path):
        if Path(path) == session_path:
            (session_path / "session.json").unlink()
            raise PermissionError("still open")
        real_rmtree(path)

    monkeypatch.setattr(artifact_runtime_module.shutil, "rmtree", fail_after_metadata_removal)
    assert runtime.close() is False
    assert json.loads((session_path / "session.json").read_text(encoding="utf-8")) == expected_metadata

    monkeypatch.setattr(artifact_runtime_module.shutil, "rmtree", real_rmtree)
    assert runtime.close() is True
    assert not session_path.exists()
