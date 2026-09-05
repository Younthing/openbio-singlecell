"""One-shot entry for the disposable R process-lifecycle experiment."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from openbio_singlecell.worker_protocol import (
    OperationContext, WorkerRequest, WorkerResponse, execute_request,
    read_json, register_operation, response_path_for, write_json,
)


class RProcessError(RuntimeError):
    pass


@register_operation("prototype.r.lifecycle")
def lifecycle(context, inputs, parameters):
    output = context.create_output_directory("lifecycle")
    evidence = Path(parameters["evidence"])
    environment = os.environ.copy()
    r_paths = [str(Path(parameters["rscript"]).parent)]
    if parameters["r_bin"]:
        r_paths.insert(0, parameters["r_bin"])
    environment["PATH"] = os.pathsep.join([*r_paths, environment["PATH"]])
    if parameters["r_home"]:
        environment["R_HOME"] = parameters["r_home"]
    command = [
        parameters["rscript"], "--vanilla", str(Path(__file__).with_name("lifecycle.R")),
        parameters["mode"], str(output), str(evidence), str(parameters["delay_seconds"]),
    ]
    # Keep the child in its parent's process group / Windows job for cancellation.
    stderr_file = evidence / "r-stderr.txt"
    with stderr_file.open("wb") as stderr:
        with subprocess.Popen(command, env=environment, stderr=stderr) as child:
            (evidence / "launcher.pid").write_text(str(child.pid), encoding="utf-8")
            returncode = child.wait()
    with stderr_file.open("rb") as stderr:
        stderr.seek(max(0, stderr_file.stat().st_size - 64 * 1024))
        error_tail = stderr.read().decode("utf-8", errors="replace")
    if error_tail:
        print(error_tail, file=sys.stderr, flush=True)
    if returncode:
        raise RProcessError(f"Rscript exited with code {returncode}: {error_tail.strip()}")
    return [{
        "type": "artifact", "name": "result", "kind": "OPENBIO_R_LIFECYCLE",
        "codec": "utf8-v1", "payload": (output / "result.txt").relative_to(context.output_root).as_posix(),
    }]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    request_path = Path(parser.parse_args().request).resolve(strict=True)
    request = WorkerRequest.from_json(read_json(request_path))
    try:
        response = execute_request(request, OperationContext.from_request_path(request_path, request.request_id))
    except Exception as error:
        traceback.print_exc(file=sys.stderr)
        response = WorkerResponse.failure(request.request_id, type(error).__name__, str(error))
    write_json(response_path_for(request_path), response.to_json())


if __name__ == "__main__":
    main()
