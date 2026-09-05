"""Run real R child-process success, failure and cancellation assertions."""
from __future__ import annotations

import argparse
import asyncio
import gc
import json
from pathlib import Path
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import psutil

from openbio_singlecell.artifact_runtime import ArtifactRuntime
from openbio_singlecell.worker_client import _invoke_worker
from openbio_singlecell.worker_protocol import WorkerRequest

DELAY_SECONDS = 4
ENTRY = Path(__file__).with_name("lifecycle_entry.py")


def process_alive(pid):
    try:
        process = psutil.Process(pid)
        return process.is_running() and process.status() != psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        return False


def request_for(args, mode, evidence):
    evidence.mkdir()
    return WorkerRequest.create("prototype.r.lifecycle", inputs={}, parameters={
        "rscript": str(Path(args.rscript).resolve(strict=True)), "r_bin": args.r_bin,
        "r_home": args.r_home,
        "mode": mode, "evidence": str(evidence), "delay_seconds": DELAY_SECONDS,
    })


async def completed_run(runtime, args, mode, evidence):
    lease = runtime.begin_run()
    partial, ready = lease.staging_path, lease.ready_path
    response, stdout, stderr = await _invoke_worker(
        sys.executable, ENTRY, partial / "request.json", request_for(args, mode, evidence),
    )
    (evidence / "worker-stderr.txt").write_text(stderr, encoding="utf-8")
    (evidence / "response.json").write_text(json.dumps(response.to_json(), indent=2), encoding="utf-8")
    r_pid = int((evidence / "r.pid").read_text())
    assert not process_alive(r_pid), "R must finish before the worker response is accepted"
    if mode == "fail":
        assert not response.ok
        assert response.error["type"] == "RProcessError"
        assert "Intentional R lifecycle failure" in response.error["message"]
        assert "Intentional R lifecycle failure" in stderr
        assert not ready.exists()
        lease.abort()
        assert not partial.exists() and not ready.exists()
        return {"r_pid": r_pid, "error": response.error, "published": False, "partial_removed": True}
    assert response.ok
    record = response.outputs[0]
    tickets = lease.publish({record["name"]: {key: record[key] for key in ("kind", "codec", "payload")}})
    assert not partial.exists() and ready.is_dir()
    content = runtime.resolve(tickets["result"]).read_text().strip()
    assert content == "complete R output"
    del tickets, lease
    gc.collect()
    assert not ready.exists(), "The final ticket/lease reference must release the published run"
    return {"r_pid": r_pid, "published_content": content, "published": True, "released_run_removed": True}


async def cancelled_run(runtime, args, evidence):
    lease = runtime.begin_run()
    partial, ready = lease.staging_path, lease.ready_path
    task = asyncio.create_task(_invoke_worker(
        sys.executable, ENTRY, partial / "request.json", request_for(args, "sleep", evidence),
    ))
    try:
        deadline = time.monotonic() + 30
        while not (evidence / "started.txt").exists():
            if task.done():
                response, _, stderr = task.result()
                raise AssertionError(f"R exited before the start marker: {response.to_json()} {stderr}")
            assert time.monotonic() < deadline, "R did not start within 30 seconds"
            await asyncio.sleep(0.02)
        r_pid = int((evidence / "r.pid").read_text())
        assert process_alive(r_pid), "Cancellation must target a live R process"
        started = time.monotonic()
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass
        else:
            raise AssertionError("The worker invocation did not propagate cancellation")
        latency = time.monotonic() - started
        assert not process_alive(r_pid), "The R child survived worker cancellation"
    finally:
        if not task.done():
            task.cancel()
            await asyncio.gather(task, return_exceptions=True)
        lease.abort()
    assert not partial.exists() and not ready.exists()
    await asyncio.sleep(DELAY_SECONDS + 0.25)
    assert not (evidence / "late.txt").exists(), "R wrote its delayed marker after cancellation"
    return {
        "r_pid": r_pid, "r_alive_before_cancel": True, "r_alive_after_cancel": False,
        "cancellation_seconds": round(latency, 3), "late_marker_absent": True,
        "late_marker_observation_seconds": DELAY_SECONDS + 0.25, "partial_removed": True,
        "published": False,
    }


async def run(args):
    report_path = Path(args.output).resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    evidence = Path(tempfile.mkdtemp(prefix="lifecycle-evidence-", dir=report_path.parent))
    with tempfile.TemporaryDirectory(prefix="r-lifecycle-runtime-") as temp_root:
        runtime = ArtifactRuntime(temp_root)
        report = {
            "python": sys.executable, "rscript": str(Path(args.rscript).resolve()),
            "evidence_directory": str(evidence),
            "success": await completed_run(runtime, args, "success", evidence / "success"),
            "failure": await completed_run(runtime, args, "fail", evidence / "fail"),
            "cancellation": await cancelled_run(runtime, args, evidence / "sleep"),
        }
        assert runtime.close(), "The experiment left live run leases"
    report["r_version"] = (evidence / "success" / "r-version.txt").read_text().strip()
    report["all_assertions_passed"] = True
    text = json.dumps(report, indent=2)
    report_path.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--r-bin", default="", help="Optional R DLL directory prepended to the child PATH")
    parser.add_argument("--r-home", default="", help="Optional R_HOME set only for the R child")
    parser.add_argument("--output", default=str(Path(__file__).with_name("lifecycle-report.json")))
    asyncio.run(run(parser.parse_args()))
