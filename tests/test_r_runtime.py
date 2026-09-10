from __future__ import annotations

import asyncio
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import openbio_singlecell.r_runtime as r_runtime
from openbio_singlecell.worker_client import _invoke_worker
from openbio_singlecell.worker_protocol import WorkerRequest


def test_probe_records_loaded_package_identity_and_changes_with_library_update(tmp_path, monkeypatch):
    package = tmp_path / "monocle"
    package.mkdir()
    description = package / "DESCRIPTION"
    description.write_text("Version: 2.40.0\n", encoding="utf-8")
    inherited_path = os.environ["PATH"]
    observed = []

    def run(argv, **kwargs):
        request = json.loads(Path(argv[-1]).read_text(encoding="utf-8"))
        observed.append((argv, kwargs, request))
        Path(request["output_dir"], "summary.json").write_text(
            json.dumps({"r_version": "R version 4.4.3", "packages": {
                "monocle": {"version": "2.40.0", "path": str(package)},
            }}), encoding="utf-8",
        )
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", run)
    options = dict(r_home="R home", library_paths="first\nsecond", path_prefix="bin1\nbin2", packages=("monocle",))
    first = r_runtime.probe_r_runtime(sys.executable, **options)
    assert first["executable"] == os.path.normcase(os.path.realpath(sys.executable))
    assert first["r_version"] == "R version 4.4.3"
    assert first["packages"]["monocle"]["description"]["path"] == str(description.resolve())
    assert len(first["packages"]["monocle"]["description"]["sha256"]) == 64
    argv, kwargs, request = observed[-1]
    assert argv[1] == "--vanilla"
    assert kwargs["env"]["R_HOME"] == "R home"
    assert kwargs["env"]["R_LIBS_USER"] == os.pathsep.join(("first", "second"))
    assert kwargs["env"]["PATH"] == os.pathsep.join(("bin1", "bin2", inherited_path))
    assert request["packages"] == ["monocle"]
    assert os.environ["PATH"] == inherited_path
    before = r_runtime.runtime_fingerprint(sys.executable, **options)
    stat = description.stat()
    description.write_text("Version: 2.40.1\n", encoding="utf-8")
    os.utime(description, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert r_runtime.runtime_fingerprint(sys.executable, **options) != before


def _runtime():
    executable = os.path.normcase(os.path.realpath(sys.executable))
    stat = Path(executable).stat()
    return {
        "executable": executable, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
        "r_home": "", "library_paths": "", "path_prefix": "", "r_version": "4.4.3", "packages": {},
    }


def test_driver_preserves_outputs_and_warnings_and_inherits_the_worker_process_tree(tmp_path, monkeypatch, capsys):
    script = tmp_path / "driver with spaces.R"
    script.write_text("# packaged driver", encoding="utf-8")
    observed_requests = []

    def run(argv, **kwargs):
        assert argv[:3] == [_runtime()["executable"], "--vanilla", str(script)]
        assert kwargs.get("shell", False) is False
        assert kwargs.get("start_new_session", False) is False
        assert not kwargs.get("creationflags")
        assert kwargs["timeout"] is None
        request = json.loads(Path(argv[-1]).read_text(encoding="utf-8"))
        observed_requests.append(Path(argv[-1]))
        assert request["formula"] == "~sm.ns(Pseudotime, df=3)"
        Path(request["output_dir"], "summary.json").write_text('{"cells":7}', encoding="utf-8")
        (tmp_path / "trajectory.rds").write_bytes(b"native model")
        kwargs["stdout"].write(b"native R output\n")
        kwargs["stderr"].write(b"Warning: native R warning\n")
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", run)
    request = {"formula": "~sm.ns(Pseudotime, df=3)", "output_dir": str(tmp_path)}
    assert r_runtime.run_r_script(_runtime(), script, request, tmp_path) == {"cells": 7}
    assert r_runtime.run_r_script(_runtime(), script, request, tmp_path) == {"cells": 7}
    assert observed_requests[0] != observed_requests[1]
    assert not any(path.exists() for path in observed_requests)
    assert (tmp_path / "trajectory.rds").read_bytes() == b"native model"
    captured = capsys.readouterr()
    assert "native R output" in captured.out
    assert "Warning: native R warning" in captured.err


def test_driver_errors_include_only_bounded_r_diagnostics_and_clean_private_files(tmp_path, monkeypatch):
    def run(argv, **kwargs):
        kwargs["stdout"].write(b"o" * 200_000)
        kwargs["stderr"].write(b"e" * 200_000 + b"\nError: supplied formula has no estimable terms")
        return subprocess.CompletedProcess(argv, 7)

    monkeypatch.setattr(subprocess, "run", run)
    with pytest.raises(r_runtime.RRuntimeError, match="formula has no estimable terms") as caught:
        r_runtime.run_r_script(_runtime(), tmp_path / "driver.R", {"output_dir": str(tmp_path)}, tmp_path)
    assert caught.value.returncode == 7
    assert len(caught.value.stdout_tail.encode()) <= r_runtime.OUTPUT_TAIL_BYTES
    assert len(caught.value.stderr_tail.encode()) <= r_runtime.OUTPUT_TAIL_BYTES
    assert list(tmp_path.iterdir()) == []


def test_runtime_validation_accepts_native_versions_but_rejects_a_replaced_package(tmp_path, monkeypatch):
    package = tmp_path / "monocle"
    package.mkdir()
    description = package / "DESCRIPTION"
    description.write_text("Version: 2.99.1\n", encoding="utf-8")
    stat = description.stat()
    runtime = {**_runtime(), "packages": {"monocle": {
        "version": "2.99.1", "path": str(package), "description": {
            "path": str(description), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "sha256": hashlib.sha256(description.read_bytes()).hexdigest(),
        },
    }}}
    assert r_runtime.validate_r_runtime(runtime) == runtime
    description.write_text("Version: 2.99.2\n", encoding="utf-8")
    os.utime(description, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    with pytest.raises(ValueError, match="monocle.*changed"):
        r_runtime.validate_r_runtime(runtime)
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: pytest.fail("stale runtime was executed"))
    with pytest.raises(ValueError, match="monocle.*changed"):
        r_runtime.run_r_script(runtime, tmp_path / "driver.R", {"output_dir": str(tmp_path)}, tmp_path)
    with pytest.raises(ValueError, match="library_paths"):
        r_runtime.validate_r_runtime({**_runtime(), "library_paths": 7})


def test_probe_timeout_retains_the_r_diagnostic(tmp_path, monkeypatch):
    def run(argv, **kwargs):
        assert kwargs["timeout"] == 60
        kwargs["stderr"].write(b"Loading required package: monocle\n")
        raise subprocess.TimeoutExpired(argv, kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", run)
    with pytest.raises(r_runtime.RRuntimeError, match="timed out.*60") as caught:
        r_runtime.probe_r_runtime(sys.executable, packages=("monocle",))
    assert "Loading required package: monocle" in str(caught.value)


@pytest.mark.skipif(not os.environ.get("OPENBIO_TEST_RSCRIPT"), reason="Set OPENBIO_TEST_RSCRIPT for native R checks")
def test_native_r_driver_probe_warning_error_and_worker_cancellation(tmp_path):
    runtime = r_runtime.probe_r_runtime(
        os.environ["OPENBIO_TEST_RSCRIPT"], r_home=os.environ.get("OPENBIO_TEST_R_HOME", ""),
        library_paths=os.environ.get("OPENBIO_TEST_R_LIBRARIES", ""),
        path_prefix=os.environ.get("OPENBIO_TEST_R_PATH", ""), packages=("Matrix", "jsonlite"),
    )
    assert "Matrix" in runtime["packages"]
    assert r_runtime.validate_r_runtime(runtime) == runtime
    script = tmp_path / "native.R"
    script.write_text(
        'req <- jsonlite::fromJSON(commandArgs(TRUE)[[1]])\n'
        'warning("native advisory")\n'
        'jsonlite::write_json(list(total=sum(c(1, 2.5, -0.5))), '
        'file.path(req$output_dir, "summary.json"), auto_unbox=TRUE)\n', encoding="utf-8",
    )
    assert r_runtime.run_r_script(runtime, script, {"output_dir": str(tmp_path)}, tmp_path) == {"total": 3}
    script.write_text('stop("native Monocle diagnostic")\n', encoding="utf-8")
    with pytest.raises(r_runtime.RRuntimeError, match="native Monocle diagnostic"):
        r_runtime.run_r_script(runtime, script, {"output_dir": str(tmp_path)}, tmp_path)
    assert not list(tmp_path.glob(".r-*"))

    script.write_text(
        'req <- jsonlite::fromJSON(commandArgs(TRUE)[[1]])\n'
        'writeLines(as.character(Sys.getpid()), file.path(req$output_dir, "started"))\n'
        'Sys.sleep(2)\n'
        'writeLines("alive", file.path(req$output_dir, "finished"))\n', encoding="utf-8",
    )
    entry = tmp_path / "entry.py"
    entry.write_text(
        'import json, sys\nfrom pathlib import Path\n'
        f'sys.path.insert(0, {str(Path(__file__).resolve().parents[1])!r})\n'
        'from openbio_singlecell.r_runtime import run_r_script\n'
        'p = json.loads(Path(sys.argv[-1]).read_text(encoding="utf-8"))["parameters"]\n'
        'run_r_script(p["runtime"], Path(p["script"]), {"output_dir":p["directory"]}, Path(p["directory"]))\n',
        encoding="utf-8",
    )

    async def cancel():
        request = WorkerRequest.create("test.r_cancel", inputs={}, parameters={
            "runtime": runtime, "script": str(script), "directory": str(tmp_path),
        })
        task = asyncio.create_task(_invoke_worker(sys.executable, entry, tmp_path / "worker-request.json", request))
        try:
            async with asyncio.timeout(15):
                while not (tmp_path / "started").exists():
                    if task.done():
                        await task
                    await asyncio.sleep(0.05)
        finally:
            task.cancel()
            with pytest.raises(asyncio.CancelledError):
                await task
        await asyncio.sleep(2.1)
        assert not (tmp_path / "finished").exists(), "R survived Worker cancellation"

    asyncio.run(cancel())
