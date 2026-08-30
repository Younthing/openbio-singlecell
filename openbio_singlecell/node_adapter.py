from __future__ import annotations

import inspect
import sys
from typing import Any

from comfy_api.latest import io

from .node_types import WorkerType
from .nodes_worker import WorkerProfile, worker_executable_fingerprint

_SCHEMA_CACHE_FIELDS = (
    "SCHEMA",
    "_DESCRIPTION",
    "_CATEGORY",
    "_EXPERIMENTAL",
    "_DEPRECATED",
    "_DEV_ONLY",
    "_API_NODE",
    "_OUTPUT_NODE",
    "_HAS_INTERMEDIATE_OUTPUT",
    "_INPUT_IS_LIST",
    "_OUTPUT_IS_LIST",
    "_RETURN_TYPES",
    "_RETURN_NAMES",
    "_OUTPUT_TOOLTIPS",
    "_NOT_IDEMPOTENT",
    "_ACCEPT_ALL_INPUTS",
)


def _worker_fingerprint(worker: WorkerProfile | None) -> tuple[object, ...]:
    if worker is None:
        return worker_executable_fingerprint(sys.executable)
    if not isinstance(worker, WorkerProfile):
        raise TypeError("worker must be an OPENBIO_WORKER value.")
    return worker.identity.fingerprint


def _override(node: type[io.ComfyNode], name: str):
    descriptor = node.__dict__.get(name)
    return descriptor.__func__ if isinstance(descriptor, classmethod) else None


def adapt_scientific_node(node: type[io.ComfyNode]) -> type[io.ComfyNode]:
    original_fingerprint = _override(node, "fingerprint_inputs")
    original_lazy_status = _override(node, "check_lazy_status")
    original_validation = _override(node, "validate_inputs")

    @classmethod
    def define_schema(cls) -> io.Schema:
        schema = node.define_schema()
        schema.inputs.append(WorkerType.Input("worker", optional=True))
        return schema

    @classmethod
    async def execute(cls, **kwargs: Any) -> io.NodeOutput:
        worker = kwargs.pop("worker", None)
        from .artifact_service import execute_artifact_node

        return await execute_artifact_node(node, worker, kwargs)

    @classmethod
    def fingerprint_inputs(cls, **kwargs: Any) -> Any:
        worker = kwargs.pop("worker", None)
        own = original_fingerprint(cls, **kwargs) if original_fingerprint is not None else None
        return own, _worker_fingerprint(worker)

    @classmethod
    def validate_inputs(cls, **kwargs: Any) -> bool | str:
        worker = kwargs.pop("worker", None)
        if worker is not None and not isinstance(worker, WorkerProfile):
            return "worker must be connected to a Python Worker node."
        return original_validation(cls, **kwargs) if original_validation is not None else True

    @classmethod
    async def check_lazy_status(cls, **kwargs: Any) -> list[str]:
        kwargs.pop("worker", None)
        result = original_lazy_status(cls, **kwargs)
        return await result if inspect.isawaitable(result) else result

    namespace = {
        "OPENBIO_WORKER_ADAPTED": True,
        "ORIGINAL_NODE_CLASS": node,
        "define_schema": define_schema,
        "execute": execute,
        "fingerprint_inputs": fingerprint_inputs,
        "validate_inputs": validate_inputs,
        "__module__": node.__module__,
        "__doc__": node.__doc__,
    }
    if original_lazy_status is not None:
        namespace["check_lazy_status"] = check_lazy_status
    namespace.update({field: None for field in _SCHEMA_CACHE_FIELDS})
    return type(node.__name__, (node,), namespace)


__all__ = ["adapt_scientific_node"]
