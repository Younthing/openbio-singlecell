from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Any

from comfy_api.latest import io

from .node_types import WorkerType
from .worker_client import WorkerIdentity, WorkerRuntimeError, canonical_python_executable, probe_python

CATEGORY = "openbio/single-cell/runtime"


@dataclass(frozen=True, slots=True)
class WorkerProfile:
    identity: WorkerIdentity

    @property
    def executable(self) -> str:
        return self.identity.executable


def worker_executable_fingerprint(python: str) -> tuple[str, int, int]:
    canonical, stat = canonical_python_executable(python)
    return canonical, int(stat.st_size), int(stat.st_mtime_ns)


class OpenBioSingleCellPythonWorker(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPythonWorker",
            display_name="Python Worker",
            category=CATEGORY,
            description="Select and validate one trusted local Python 3.12 executable for scientific nodes.",
            inputs=[io.String.Input("python", default=sys.executable)],
            outputs=[WorkerType.Output(display_name="worker")],
        )

    @classmethod
    def validate_inputs(cls, python: str) -> bool | str:
        try:
            canonical_python_executable(python)
        except (TypeError, ValueError, OSError, WorkerRuntimeError) as error:
            return str(error)
        return True

    @classmethod
    def fingerprint_inputs(cls, python: str) -> Any:
        return worker_executable_fingerprint(python)

    @classmethod
    async def execute(cls, python: str) -> io.NodeOutput:
        return io.NodeOutput(WorkerProfile(await probe_python(python)))


WORKER_NODE_CLASSES = [OpenBioSingleCellPythonWorker]


__all__ = [
    "OpenBioSingleCellPythonWorker",
    "WORKER_NODE_CLASSES",
    "WorkerProfile",
    "worker_executable_fingerprint",
]
