from __future__ import annotations

import argparse
import importlib
import sys
import traceback
from pathlib import Path

# Isolated mode intentionally omits the script directory. Add only this plugin's
# package root so an alternate trusted Python can execute the checked-in worker.
PLUGIN_ROOT = Path(__file__).resolve().parents[1]
if str(PLUGIN_ROOT) not in sys.path:
    sys.path.insert(0, str(PLUGIN_ROOT))

from openbio_singlecell.worker_protocol import (  # noqa: E402
    PROBE_OPERATION,
    OperationContext,
    WorkerRequest,
    WorkerResponse,
    execute_request,
    read_json,
    response_path_for,
    write_json,
)


def _load_scientific_operations() -> None:
    module_name = "openbio_singlecell.worker_operations"
    try:
        importlib.import_module(module_name)
    except ModuleNotFoundError as error:
        if error.name != module_name:
            raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="openbio-worker")
    parser.add_argument("--request", required=True)
    arguments = parser.parse_args(argv)
    request_path = Path(arguments.request).resolve(strict=True)

    try:
        request = WorkerRequest.from_json(read_json(request_path))
    except Exception as error:
        print(f"Invalid worker request: {error}", file=sys.stderr, flush=True)
        return 2

    try:
        if request.operation != PROBE_OPERATION:
            _load_scientific_operations()
        context = OperationContext.from_request_path(request_path, request.request_id)
        response = execute_request(request, context)
    except Exception as error:
        traceback.print_exc(file=sys.stderr)
        response = WorkerResponse.failure(request.request_id, type(error).__name__, str(error))

    try:
        write_json(response_path_for(request_path), response.to_json())
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
