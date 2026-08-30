from __future__ import annotations

import asyncio
import os
import sys
import textwrap
import uuid
from pathlib import Path

import pytest

import openbio_singlecell.worker_client as worker_client_module
from openbio_singlecell.worker_client import (
    OUTPUT_TAIL_BYTES,
    WORKER_ENTRY,
    WorkerOperationError,
    WorkerProcessError,
    _invoke_worker,
    probe_python,
    run_worker,
)
from openbio_singlecell.worker_protocol import (
    PROBE_OPERATION,
    PROTOCOL_VERSION,
    OperationContext,
    ProtocolError,
    WorkerRequest,
    WorkerResponse,
    execute_request,
    read_json,
    register_operation,
    registered_operation_ids,
    response_path_for,
    write_json,
)


def _write_script(path: Path, source: str) -> Path:
    path.write_text(textwrap.dedent(source), encoding="utf-8")
    return path


def test_protocol_is_strict_json_and_registry_is_an_allowlist(tmp_path: Path):
    request_id = str(uuid.uuid4())
    request = WorkerRequest.create(
        "test.add",
        inputs={"source": {"path": str(tmp_path / "input.bin"), "codec": "bytes"}},
        parameters={"left": 2, "right": 3},
        request_id=request_id,
    )
    request_path = tmp_path / "request.json"
    write_json(request_path, request.to_json())
    assert WorkerRequest.from_json(read_json(request_path)) == request
    ordered_path = tmp_path / "ordered.json"
    write_json(ordered_path, {"z": 1, "a": 2})
    ordered_text = ordered_path.read_text(encoding="utf-8")
    assert ordered_text.index('"z"') < ordered_text.index('"a"')

    with pytest.raises(ProtocolError, match="finite"):
        write_json(tmp_path / "nan.json", {"value": float("nan")})
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"version":1,"version":1}', encoding="utf-8")
    with pytest.raises(ProtocolError, match="Duplicate JSON key"):
        read_json(duplicate)
    with pytest.raises(ProtocolError, match="unexpected fields"):
        WorkerRequest.from_json({**request.to_json(), "shell": "powershell"})

    register_operation(
        "test.add",
        lambda context, _inputs, parameters: {
            "value": parameters["left"] + parameters["right"],
            "output_root": str(context.output_root),
        },
    )
    context = OperationContext.from_request_path(request_path, request_id)
    assert execute_request(request, context) == WorkerResponse.success(
        request_id,
        {"value": 5, "output_root": str(tmp_path)},
    )
    register_operation("test.rows", lambda _context, _inputs, _parameters: [{"row": 1}, {"row": 2}])
    assert {PROBE_OPERATION, "test.add", "test.rows"}.issubset(registered_operation_ids())
    rows_request = WorkerRequest.create("test.rows", inputs={}, parameters={})
    rows_context = OperationContext.from_request_path(request_path, rows_request.request_id)
    assert execute_request(rows_request, rows_context).outputs == [
        {"row": 1},
        {"row": 2},
    ]
    empty_rows = WorkerResponse.success(str(uuid.uuid4()), [])
    assert WorkerResponse.from_json(empty_rows.to_json()).outputs == []
    with pytest.raises(ProtocolError, match="missing outputs"):
        WorkerResponse(request_id=str(uuid.uuid4()), ok=True).to_json()
    with pytest.raises(ProtocolError, match="not allowlisted"):
        unknown = WorkerRequest.create("test.unknown", inputs={}, parameters={})
        execute_request(unknown, OperationContext.from_request_path(request_path, unknown.request_id))


def test_operation_context_confines_output_directories_to_the_run_staging_root(tmp_path: Path):
    request_id = str(uuid.uuid4())
    request_path = tmp_path / "request.json"
    context = OperationContext.from_request_path(request_path, request_id)

    output = context.create_output_directory("adata")

    assert output == tmp_path / "outputs" / "adata"
    with pytest.raises(ProtocolError, match="output name"):
        context.create_output_directory("../escape")


def test_probe_uses_isolated_absolute_entry_and_returns_executable_identity():
    identity = asyncio.run(probe_python(sys.executable))

    executable = Path(identity.executable)
    stat = executable.stat()
    assert executable.is_absolute()
    assert WORKER_ENTRY.is_absolute()
    assert identity.protocol_version == PROTOCOL_VERSION
    assert identity.python_version[:2] >= (3, 12)
    assert identity.size == stat.st_size
    assert identity.mtime_ns == stat.st_mtime_ns
    assert {"anndata", "numpy", "pandas", "scanpy", "scipy"} == dict(identity.dependencies).keys()
    assert identity.fingerprint[0] == os.path.normcase(os.path.realpath(sys.executable))


def test_entry_returns_a_protocol_error_for_an_unallowlisted_operation(tmp_path: Path):
    identity = asyncio.run(probe_python(sys.executable))
    request_path = tmp_path / "request.json"

    with pytest.raises(WorkerOperationError, match="not allowlisted"):
        asyncio.run(
            run_worker(
                identity,
                request_path,
                "test.not_registered",
                inputs={},
                parameters={},
            )
        )
    assert not request_path.exists()
    assert not response_path_for(request_path).exists()


def test_nonzero_worker_keeps_only_bounded_output_tails(tmp_path: Path):
    entry = _write_script(
        tmp_path / "no_response.py",
        """
        import sys
        sys.stdout.write("o" * 200_000)
        sys.stdout.flush()
        sys.stderr.write("e" * 200_000 + " final-error")
        sys.stderr.flush()
        raise SystemExit(7)
        """,
    )
    request = WorkerRequest.create(PROBE_OPERATION, inputs={}, parameters={})

    with pytest.raises(WorkerProcessError) as caught:
        asyncio.run(_invoke_worker(sys.executable, entry, tmp_path / "request.json", request))

    error = caught.value
    assert error.returncode == 7
    assert len(error.stdout_tail.encode()) <= OUTPUT_TAIL_BYTES
    assert len(error.stderr_tail.encode()) <= OUTPUT_TAIL_BYTES
    assert error.stderr_tail.endswith(" final-error")


@pytest.mark.skipif(os.name == "nt" and not Path(os.environ.get("SystemRoot", "C:\\Windows"), "System32", "taskkill.exe").is_file(), reason="taskkill is unavailable")
def test_cancellation_terminates_the_worker_process_tree(tmp_path: Path):
    child_marker = tmp_path / "child-finished"
    started_marker = tmp_path / "started"
    entry = _write_script(
        tmp_path / "tree.py",
        """
        import argparse
        import json
        import subprocess
        import sys
        import time
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument("--request", required=True)
        request = json.loads(Path(parser.parse_args().request).read_text(encoding="utf-8"))
        marker = request["parameters"]["marker"]
        started = request["parameters"]["started"]
        child = "import time; from pathlib import Path; time.sleep(2); Path(r'%s').write_text('alive')" % marker
        subprocess.Popen([sys.executable, "-I", "-c", child])
        Path(started).write_text("started", encoding="utf-8")
        time.sleep(60)
        """,
    )
    request = WorkerRequest.create(
        "test.tree",
        inputs={},
        parameters={"marker": str(child_marker), "started": str(started_marker)},
    )

    async def cancel_after_start() -> None:
        task = asyncio.create_task(_invoke_worker(sys.executable, entry, tmp_path / "request.json", request))
        for _ in range(100):
            if started_marker.exists():
                break
            await asyncio.sleep(0.02)
        assert started_marker.exists()
        task.cancel()
        with pytest.raises(asyncio.CancelledError):
            await task
        await asyncio.sleep(2.2)

    asyncio.run(cancel_after_start())
    assert not child_marker.exists()


def test_posix_cancellation_terminates_process_group_after_parent_exit(monkeypatch):
    class ExitedWorker:
        pid = 4242
        returncode = 0

        async def wait(self) -> int:
            return self.returncode

    signals = []
    monkeypatch.setattr(worker_client_module.os, "name", "posix")
    monkeypatch.setattr(worker_client_module.signal, "SIGKILL", 9, raising=False)
    monkeypatch.setattr(worker_client_module, "_POSIX_TERMINATION_GRACE_SECONDS", 0.0, raising=False)
    monkeypatch.setattr(
        worker_client_module.os,
        "killpg",
        lambda pid, sent_signal: signals.append((pid, sent_signal)),
        raising=False,
    )

    asyncio.run(worker_client_module._terminate_process_tree(ExitedWorker()))

    assert signals == [
        (4242, worker_client_module.signal.SIGTERM),
        (4242, 0),
        (4242, worker_client_module.signal.SIGKILL),
    ]


@pytest.mark.skipif(os.name != "nt", reason="Windows orphaned-process tree regression")
def test_cancellation_terminates_children_after_worker_parent_has_exited(tmp_path: Path):
    import psutil

    child_marker = tmp_path / "orphan-child-finished"
    started_marker = tmp_path / "orphan-started"
    entry = _write_script(
        tmp_path / "orphan_tree.py",
        """
        import argparse
        import json
        import os
        import subprocess
        import sys
        import time
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument("--request", required=True)
        request = json.loads(Path(parser.parse_args().request).read_text(encoding="utf-8"))
        marker = request["parameters"]["marker"]
        started = request["parameters"]["started"]
        child = "import time; from pathlib import Path; time.sleep(2); Path(r'%s').write_text('alive')" % marker
        time.sleep(0.25)
        subprocess.Popen([sys.executable, "-I", "-c", child])
        Path(started).write_text(str(os.getpid()), encoding="utf-8")
        """,
    )
    request = WorkerRequest.create(
        "test.orphan-tree",
        inputs={},
        parameters={"marker": str(child_marker), "started": str(started_marker)},
    )

    async def cancel_after_parent_exit() -> None:
        task = asyncio.create_task(_invoke_worker(sys.executable, entry, tmp_path / "request.json", request))
        for _ in range(200):
            if started_marker.exists():
                break
            await asyncio.sleep(0.01)
        assert started_marker.exists()
        parent_pid = int(started_marker.read_text(encoding="utf-8"))
        for _ in range(200):
            if not psutil.pid_exists(parent_pid):
                break
            await asyncio.sleep(0.01)
        assert not psutil.pid_exists(parent_pid)
        assert not task.done()
        task.cancel()
        with pytest.raises(asyncio.CancelledError):
            await task
        await asyncio.sleep(2.2)

    asyncio.run(cancel_after_parent_exit())
    assert not child_marker.exists()


def test_response_must_exist_and_match_the_request(tmp_path: Path):
    silent_entry = _write_script(tmp_path / "silent.py", "raise SystemExit(0)\n")
    silent_request = WorkerRequest.create(PROBE_OPERATION, inputs={}, parameters={})
    with pytest.raises(WorkerProcessError, match="without a response"):
        asyncio.run(_invoke_worker(sys.executable, silent_entry, tmp_path / "silent-request.json", silent_request))

    entry = _write_script(
        tmp_path / "wrong_response.py",
        """
        import argparse
        import json
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument("--request", required=True)
        request_path = Path(parser.parse_args().request)
        response_path = request_path.with_name("response.json")
        response_path.write_text(json.dumps({
            "version": 1,
            "request_id": "00000000-0000-0000-0000-000000000000",
            "ok": True,
            "outputs": {},
        }), encoding="utf-8")
        """,
    )
    request = WorkerRequest.create(PROBE_OPERATION, inputs={}, parameters={})

    with pytest.raises(ProtocolError, match="request ID"):
        asyncio.run(_invoke_worker(sys.executable, entry, tmp_path / "request.json", request))


def test_response_path_is_fixed_next_to_request(tmp_path: Path):
    request_path = tmp_path / "control" / "request.json"
    assert response_path_for(request_path) == request_path.with_name("response.json")
