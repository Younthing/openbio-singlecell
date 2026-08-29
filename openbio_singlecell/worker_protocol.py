from __future__ import annotations

import importlib
import json
import math
import os
import re
import sys
import uuid
from collections.abc import Callable
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path

PROTOCOL_VERSION = 1
PROBE_OPERATION = "openbio.probe"
CORE_DEPENDENCIES = ("anndata", "numpy", "pandas", "scanpy", "scipy")

type JSONScalar = None | bool | int | float | str
type JSONValue = JSONScalar | list[JSONValue] | dict[str, JSONValue]
type JSONContainer = list[JSONValue] | dict[str, JSONValue]
type Operation = Callable[[OperationContext, dict[str, JSONValue], dict[str, JSONValue]], JSONContainer]

_OPERATION_ID = re.compile(r"[a-z][a-z0-9_.-]{0,127}\Z")
_OPERATIONS: dict[str, Operation] = {}


class ProtocolError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class OperationContext:
    output_root: Path
    request_id: str

    @classmethod
    def from_request_path(cls, request_path: str | os.PathLike[str], request_id: str) -> OperationContext:
        _require_request_id(request_id)
        root = Path(request_path).resolve().parent
        if not root.is_dir():
            raise ProtocolError("Worker request parent must be an existing staging directory.")
        return cls(output_root=root, request_id=request_id)

    def create_output_directory(self, output_name: str) -> Path:
        if not isinstance(output_name, str) or re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", output_name) is None:
            raise ProtocolError("Worker output name is invalid.")
        outputs = self.output_root / "outputs"
        outputs.mkdir(exist_ok=True)
        target = outputs / output_name
        target.mkdir(exist_ok=False)
        return target


def _validate_json(value: object, *, location: str = "JSON value") -> None:
    if value is None or isinstance(value, (bool, str)):
        return
    if isinstance(value, int):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ProtocolError(f"{location} must contain only finite numbers.")
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _validate_json(item, location=f"{location}[{index}]")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ProtocolError(f"{location} object keys must be strings.")
            _validate_json(item, location=f"{location}.{key}")
        return
    raise ProtocolError(f"{location} contains a non-JSON value of type {type(value).__name__}.")


def _reject_constant(value: str) -> None:
    raise ProtocolError(f"JSON must contain only finite numbers, not {value}.")


def _object_without_duplicates(pairs: list[tuple[str, JSONValue]]) -> dict[str, JSONValue]:
    result: dict[str, JSONValue] = {}
    for key, value in pairs:
        if key in result:
            raise ProtocolError(f"Duplicate JSON key: {key!r}.")
        result[key] = value
    return result


def read_json(path: str | os.PathLike[str]) -> JSONValue:
    try:
        value = json.loads(
            Path(path).read_text(encoding="utf-8"),
            parse_constant=_reject_constant,
            object_pairs_hook=_object_without_duplicates,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ProtocolError(f"Could not read strict JSON from {Path(path)}: {error}") from error
    _validate_json(value)
    return value


def write_json(path: str | os.PathLike[str], value: JSONValue) -> None:
    _validate_json(value)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = destination.with_name(f".{destination.name}.{uuid.uuid4()}.tmp")
    try:
        with staging.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(value, handle, allow_nan=False, ensure_ascii=False, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(staging, destination)
    finally:
        try:
            staging.unlink()
        except FileNotFoundError:
            pass


def _require_object(value: object, description: str) -> dict[str, JSONValue]:
    if not isinstance(value, dict):
        raise ProtocolError(f"{description} must be a JSON object.")
    _validate_json(value, location=description)
    return value


def _require_container(value: object, description: str) -> JSONContainer:
    if not isinstance(value, (dict, list)):
        raise ProtocolError(f"{description} must be a JSON object or array.")
    _validate_json(value, location=description)
    return value


def _require_exact_fields(value: dict[str, JSONValue], expected: set[str], description: str) -> None:
    actual = set(value)
    if actual != expected:
        unexpected = sorted(actual - expected)
        missing = sorted(expected - actual)
        details = []
        if unexpected:
            details.append(f"unexpected fields {unexpected}")
        if missing:
            details.append(f"missing fields {missing}")
        raise ProtocolError(f"{description} has {' and '.join(details)}.")


def _require_version(value: object) -> None:
    if type(value) is not int or value != PROTOCOL_VERSION:
        raise ProtocolError(f"Unsupported worker protocol version: {value!r}.")


def _require_request_id(value: object) -> str:
    if not isinstance(value, str):
        raise ProtocolError("Worker request ID must be a UUID string.")
    try:
        parsed = uuid.UUID(value)
    except ValueError as error:
        raise ProtocolError("Worker request ID must be a UUID string.") from error
    if str(parsed) != value:
        raise ProtocolError("Worker request ID must use canonical lowercase UUID form.")
    return value


def _require_operation_id(value: object) -> str:
    if not isinstance(value, str) or _OPERATION_ID.fullmatch(value) is None:
        raise ProtocolError("Worker operation ID is invalid.")
    return value


@dataclass(frozen=True, slots=True)
class WorkerRequest:
    request_id: str
    operation: str
    inputs: dict[str, JSONValue]
    parameters: dict[str, JSONValue]
    version: int = PROTOCOL_VERSION

    @classmethod
    def create(
        cls,
        operation: str,
        *,
        inputs: dict[str, JSONValue],
        parameters: dict[str, JSONValue],
        request_id: str | None = None,
    ) -> WorkerRequest:
        return cls.from_json(
            {
                "version": PROTOCOL_VERSION,
                "request_id": request_id or str(uuid.uuid4()),
                "operation": operation,
                "inputs": inputs,
                "parameters": parameters,
            }
        )

    @classmethod
    def from_json(cls, value: object) -> WorkerRequest:
        data = _require_object(value, "Worker request")
        _require_exact_fields(data, {"version", "request_id", "operation", "inputs", "parameters"}, "Worker request")
        _require_version(data["version"])
        request_id = _require_request_id(data["request_id"])
        operation = _require_operation_id(data["operation"])
        inputs = _require_object(data["inputs"], "Worker request inputs")
        parameters = _require_object(data["parameters"], "Worker request parameters")
        return cls(request_id=request_id, operation=operation, inputs=inputs, parameters=parameters)

    def to_json(self) -> dict[str, JSONValue]:
        return {
            "version": self.version,
            "request_id": self.request_id,
            "operation": self.operation,
            "inputs": self.inputs,
            "parameters": self.parameters,
        }


@dataclass(frozen=True, slots=True)
class WorkerResponse:
    request_id: str
    ok: bool
    outputs: JSONContainer | None = None
    error: dict[str, JSONValue] | None = None
    version: int = PROTOCOL_VERSION

    @classmethod
    def success(cls, request_id: str, outputs: JSONContainer) -> WorkerResponse:
        _require_request_id(request_id)
        return cls(request_id=request_id, ok=True, outputs=_require_container(outputs, "Worker outputs"))

    @classmethod
    def failure(cls, request_id: str, error_type: str, message: str) -> WorkerResponse:
        _require_request_id(request_id)
        if not isinstance(error_type, str) or not error_type:
            raise ProtocolError("Worker error type must be a non-empty string.")
        if not isinstance(message, str):
            raise ProtocolError("Worker error message must be a string.")
        return cls(request_id=request_id, ok=False, error={"type": error_type, "message": message})

    @classmethod
    def from_json(cls, value: object) -> WorkerResponse:
        data = _require_object(value, "Worker response")
        common = {"version", "request_id", "ok"}
        if type(data.get("ok")) is not bool:
            raise ProtocolError("Worker response ok must be a boolean.")
        expected = common | ({"outputs"} if data["ok"] else {"error"})
        _require_exact_fields(data, expected, "Worker response")
        _require_version(data["version"])
        request_id = _require_request_id(data["request_id"])
        if data["ok"]:
            return cls.success(request_id, _require_container(data["outputs"], "Worker outputs"))
        error = _require_object(data["error"], "Worker error")
        _require_exact_fields(error, {"type", "message"}, "Worker error")
        error_type = error["type"]
        message = error["message"]
        if not isinstance(error_type, str) or not error_type or not isinstance(message, str):
            raise ProtocolError("Worker error must contain a non-empty string type and string message.")
        return cls.failure(request_id, error_type, message)

    def to_json(self) -> dict[str, JSONValue]:
        value: dict[str, JSONValue] = {
            "version": self.version,
            "request_id": self.request_id,
            "ok": self.ok,
        }
        if self.ok:
            if self.outputs is None:
                raise ProtocolError("Successful worker response is missing outputs.")
            value["outputs"] = self.outputs
        else:
            if self.error is None:
                raise ProtocolError("Failed worker response is missing an error.")
            value["error"] = self.error
        return value


def register_operation(operation_id: str, operation: Operation | None = None):
    operation_id = _require_operation_id(operation_id)

    def register(candidate: Operation) -> Operation:
        if not callable(candidate):
            raise TypeError("Worker operation must be callable.")
        if operation_id in _OPERATIONS:
            raise ProtocolError(f"Worker operation is already registered: {operation_id}")
        _OPERATIONS[operation_id] = candidate
        return candidate

    return register(operation) if operation is not None else register


def registered_operation_ids() -> tuple[str, ...]:
    return tuple(sorted(_OPERATIONS))


def execute_request(request: WorkerRequest, context: OperationContext) -> WorkerResponse:
    if not isinstance(request, WorkerRequest):
        raise TypeError("Expected a WorkerRequest.")
    operation = _OPERATIONS.get(request.operation)
    if operation is None:
        raise ProtocolError(f"Worker operation is not allowlisted: {request.operation}")
    if not isinstance(context, OperationContext) or context.request_id != request.request_id:
        raise ProtocolError("Worker operation context does not match its request.")
    outputs = operation(context, request.inputs, request.parameters)
    if not isinstance(outputs, (dict, list)):
        raise ProtocolError(f"Worker operation {request.operation} did not return a JSON object or array.")
    return WorkerResponse.success(request.request_id, outputs)


def response_path_for(request_path: str | os.PathLike[str]) -> Path:
    return Path(request_path).with_name("response.json")


@register_operation(PROBE_OPERATION)
def _probe(
    _context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> dict[str, JSONValue]:
    if inputs or parameters:
        raise ProtocolError("Worker probe does not accept inputs or parameters.")
    versions: dict[str, JSONValue] = {}
    for dependency in CORE_DEPENDENCIES:
        importlib.import_module(dependency)
        versions[dependency] = metadata.version(dependency)
    return {
        "protocol_version": PROTOCOL_VERSION,
        "python_version": [sys.version_info.major, sys.version_info.minor, sys.version_info.micro],
        "dependencies": versions,
        "isolated": bool(sys.flags.isolated),
    }


__all__ = [
    "CORE_DEPENDENCIES",
    "JSONContainer",
    "JSONValue",
    "OperationContext",
    "PROBE_OPERATION",
    "PROTOCOL_VERSION",
    "ProtocolError",
    "WorkerRequest",
    "WorkerResponse",
    "execute_request",
    "read_json",
    "register_operation",
    "registered_operation_ids",
    "response_path_for",
    "write_json",
]
