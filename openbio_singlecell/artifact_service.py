from __future__ import annotations

import asyncio
import json
from pathlib import Path
from typing import Any

from comfy_api.latest import io

from .artifact_envelope import summary_from_metadata
from .artifact_runtime import NATIVE_AFFINITY_CODECS, ArtifactRuntime, ArtifactTicket
from .files import (
    input_file_provenance,
    resolve_input_path,
    tenx_mtx_provenance,
    tenx_study_provenance,
)
from .worker_client import WorkerIdentity, probe_python, run_worker

ANNDATA_CODEC = "anndata-h5ad-v1"
_COMFY_INTERRUPT_POLL_SECONDS = 0.25

_runtime: ArtifactRuntime | None = None
_default_worker: WorkerIdentity | None = None

_FILE_INPUT_RULES = {
    ("OpenBioSingleCellLoadH5AD", "path"): (".h5ad",),
    ("OpenBioSingleCellLoad10xH5", "path"): (".h5", ".hdf5"),
    ("OpenBioSingleCellMapGeneIdsFromGTF", "gtf_path"): (".gtf", ".gtf.gz"),
    ("OpenBioSingleCellMarkerORAEvidence", "resource_csv"): (".csv",),
    ("OpenBioSingleCellLianaCommunication", "resource_csv"): (".csv",),
    ("OpenBioSingleCellCollecTRIULM", "network_csv"): (".csv",),
    ("OpenBioSingleCellImportPySCENICResults", "run_manifest_json"): (".json",),
    ("OpenBioSingleCellAUCellScores", "gene_sets_file"): (".csv", ".tsv", ".gmt"),
    ("OpenBioSingleCellGSVAScores", "gene_sets_file"): (".csv", ".tsv", ".gmt"),
    ("OpenBioSingleCellGenePanelScores", "gene_sets_file"): (".csv", ".tsv", ".gmt"),
    ("OpenBioSingleCellRankedGSEA", "gene_sets_file"): (".csv", ".tsv", ".gmt"),
    ("OpenBioSingleCellGeneSetOverrepresentation", "gene_sets_file"): (".csv", ".tsv", ".gmt"),
    ("OpenBioSingleCellDGIdbAnnotation", "resource_file"): (".csv", ".tsv"),
    ("OpenBioSingleCellCassiopeiaLineageQC", "allele_table_file"): (".csv", ".tsv"),
}
_DIRECTORY_INPUT_RULES = {
    ("OpenBioSingleCellLoad10xMTX", "directory"): tenx_mtx_provenance,
    ("OpenBioSingleCellLoad10xStudy", "directory"): tenx_study_provenance,
}


def operation_id_for_node(node_id: str) -> str:
    prefix = "OpenBioSingleCell"
    if not isinstance(node_id, str) or not node_id.startswith(prefix) or len(node_id) == len(prefix):
        raise ValueError(f"Invalid OpenBio scientific node ID: {node_id!r}")
    return f"openbio.node.{node_id.removeprefix(prefix).lower()}"


async def initialize_artifact_service(temp_root: str | Path, python: str) -> None:
    global _default_worker, _runtime

    identity = await probe_python(python)
    previous = _runtime
    if previous is not None and not previous.close():
        raise RuntimeError("Cannot replace the Artifact Runtime while workflow tickets are still live.")
    _runtime = ArtifactRuntime(temp_root)
    _default_worker = identity


def current_artifact_runtime() -> ArtifactRuntime:
    if _runtime is None:
        raise RuntimeError("Artifact Runtime is not initialized; OpenBio extension on_load did not complete.")
    return _runtime


def _current_default_worker() -> WorkerIdentity:
    if _default_worker is None:
        raise RuntimeError("Default Python Worker is not initialized; OpenBio extension on_load did not complete.")
    return _default_worker


def _identity_affinity(identity: WorkerIdentity) -> str:
    return json.dumps(
        {
            "dependencies": [list(item) for item in identity.dependencies],
            "executable": identity.executable,
            "mtime_ns": identity.mtime_ns,
            "protocol_version": identity.protocol_version,
            "python_version": list(identity.python_version),
            "size": identity.size,
        },
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def _identity_from_affinity(value: object) -> WorkerIdentity:
    try:
        data = json.loads(str(value))
    except (TypeError, ValueError) as error:
        raise ValueError("Native artifact Worker affinity is unreadable.") from error
    expected = {"dependencies", "executable", "mtime_ns", "protocol_version", "python_version", "size"}
    if not isinstance(data, dict) or set(data) != expected:
        raise ValueError("Native artifact Worker affinity schema is invalid.")
    dependencies = data["dependencies"]
    python_version = data["python_version"]
    if (
        not isinstance(data["executable"], str)
        or not data["executable"]
        or type(data["size"]) is not int
        or type(data["mtime_ns"]) is not int
        or type(data["protocol_version"]) is not int
        or not isinstance(python_version, list)
        or len(python_version) != 3
        or any(type(part) is not int for part in python_version)
        or not isinstance(dependencies, list)
        or any(
            not isinstance(item, list)
            or len(item) != 2
            or not all(isinstance(part, str) and part for part in item)
            for item in dependencies
        )
    ):
        raise ValueError("Native artifact Worker affinity schema is invalid.")
    return WorkerIdentity(
        executable=data["executable"],
        size=data["size"],
        mtime_ns=data["mtime_ns"],
        protocol_version=data["protocol_version"],
        python_version=tuple(python_version),
        dependencies=tuple(tuple(item) for item in dependencies),
    )


def _select_worker(
    worker: WorkerIdentity | None,
    ticket_inputs: dict[str, ArtifactTicket],
    runtime: ArtifactRuntime,
) -> WorkerIdentity:
    if worker is not None and not isinstance(worker, WorkerIdentity):
        raise TypeError("worker must be an OPENBIO_WORKER value.")
    affinities = {
        str(record["worker_affinity"])
        for ticket in ticket_inputs.values()
        if "worker_affinity" in (record := runtime.record(ticket))
    }
    if len(affinities) > 1:
        raise ValueError("Native artifact inputs require different Python Workers.")
    if not affinities:
        return worker or _current_default_worker()

    affinity = affinities.pop()
    if worker is not None:
        if _identity_affinity(worker) != affinity:
            raise ValueError("Selected Python Worker does not match the native artifact producer.")
        return worker
    return _identity_from_affinity(affinity)


def _file_descriptor(node_id: str, name: str, value: object) -> dict[str, object] | None:
    extensions = _FILE_INPUT_RULES.get((node_id, name))
    if extensions is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{node_id} {name} must be a string path.")
    if not value.strip():
        return None
    return {
        "type": "file",
        "path": resolve_input_path(value, extensions=extensions),
        "provenance": input_file_provenance(value, extensions),
    }


def _directory_descriptor(node_id: str, name: str, value: object) -> dict[str, object] | None:
    provenance = _DIRECTORY_INPUT_RULES.get((node_id, name))
    if provenance is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{node_id} {name} must be a string path.")
    return {
        "type": "directory",
        "path": resolve_input_path(value, kind="directory"),
        "provenance": provenance(value),
    }


def _prepare_request(
    node_id: str,
    kwargs: dict[str, Any],
    runtime: ArtifactRuntime,
) -> tuple[dict[str, object], dict[str, object], dict[str, ArtifactTicket]]:
    inputs: dict[str, object] = {}
    parameters: dict[str, object] = {}
    tickets: dict[str, ArtifactTicket] = {}
    for name, value in kwargs.items():
        if isinstance(value, ArtifactTicket):
            path = runtime.resolve(value)
            inputs[name] = {
                "type": "artifact",
                "path": str(path),
                "kind": value.kind,
                "codec": value.codec,
            }
            tickets[name] = value
            continue
        descriptor = _file_descriptor(node_id, name, value)
        if descriptor is None:
            descriptor = _directory_descriptor(node_id, name, value)
        if descriptor is not None:
            inputs[name] = descriptor
        else:
            parameters[name] = value
    return inputs, parameters, tickets


def _record_fields(record_type: str) -> set[str]:
    if record_type == "artifact":
        return {"type", "name", "kind", "codec", "payload"}
    if record_type in {"summary", "string"}:
        return {"type", "name", "value"}
    if record_type == "input_ref":
        return {"type", "name", "input"}
    raise ValueError(f"Worker returned an unsupported output record type: {record_type!r}")


def _validate_response_records(
    node: type[io.ComfyNode],
    records: object,
    ticket_inputs: dict[str, ArtifactTicket],
) -> list[dict[str, object]]:
    schema = node.define_schema()
    expected = [(output.display_name, output.io_type) for output in schema.outputs]
    if not isinstance(records, list) or len(records) != len(expected):
        raise ValueError(f"Worker returned the wrong number of outputs for {schema.node_id}.")
    validated: list[dict[str, object]] = []
    for position, (record, (expected_name, expected_kind)) in enumerate(zip(records, expected, strict=True)):
        if not isinstance(record, dict) or not isinstance(record.get("type"), str):
            raise ValueError(f"Worker output {position} is not a typed record.")
        record_type = record["type"]
        if set(record) != _record_fields(record_type):
            raise ValueError(f"Worker output {position} has an invalid {record_type!r} record schema.")
        if record["name"] != expected_name:
            raise ValueError(
                f"Worker output {position} is named {record['name']!r}; expected {expected_name!r}."
            )
        if record_type == "artifact":
            for field in ("kind", "codec", "payload"):
                if not isinstance(record[field], str) or not record[field]:
                    raise ValueError(
                        f"Worker artifact {expected_name!r} {field} must be a non-empty string."
                    )
            if record["kind"] != expected_kind:
                raise ValueError(f"Worker artifact {expected_name!r} does not match wire type {expected_kind}.")
        if record_type == "summary" and expected_kind != "OPENBIO_SINGLE_CELL_SUMMARY":
            raise ValueError(f"Worker summary {expected_name!r} does not match wire type {expected_kind}.")
        if record_type == "string":
            if expected_kind != "STRING":
                raise ValueError(f"Worker string {expected_name!r} does not match wire type {expected_kind}.")
            if not isinstance(record["value"], str):
                raise ValueError(f"Worker string {expected_name!r} must contain a string value.")
        if record_type == "input_ref":
            source = record["input"]
            if not isinstance(source, str) or not source:
                raise ValueError(f"Worker input reference {expected_name!r} must name an input ticket.")
            try:
                ticket = ticket_inputs[source]
            except KeyError as error:
                raise ValueError(f"Worker referenced unknown input ticket {source!r}.") from error
            if ticket.kind != expected_kind:
                raise ValueError(
                    f"Worker input reference {expected_name!r} does not match wire type {expected_kind}."
                )
        validated.append(record)
    return validated


async def _run_worker_with_comfy_interrupt(
    identity: WorkerIdentity,
    request_path: Path,
    operation: str,
    *,
    inputs: dict[str, object],
    parameters: dict[str, object],
):
    import comfy.model_management

    worker_task = asyncio.create_task(
        run_worker(
            identity,
            request_path,
            operation,
            inputs=inputs,
            parameters=parameters,
        )
    )
    try:
        while True:
            if comfy.model_management.processing_interrupted():
                worker_task.cancel()
                await asyncio.shield(asyncio.gather(worker_task, return_exceptions=True))
                raise comfy.model_management.InterruptProcessingException()
            done, _pending = await asyncio.wait(
                {worker_task},
                timeout=_COMFY_INTERRUPT_POLL_SECONDS,
            )
            if comfy.model_management.processing_interrupted():
                worker_task.cancel()
                await asyncio.shield(asyncio.gather(worker_task, return_exceptions=True))
                raise comfy.model_management.InterruptProcessingException()
            if done:
                return worker_task.result()
    except asyncio.CancelledError:
        worker_task.cancel()
        await asyncio.shield(asyncio.gather(worker_task, return_exceptions=True))
        raise


async def execute_artifact_node(
    node: type[io.ComfyNode],
    worker: WorkerIdentity | None,
    kwargs: dict[str, Any],
) -> io.NodeOutput:
    if not isinstance(node, type) or not issubclass(node, io.ComfyNode):
        raise TypeError("Expected an OpenBio Comfy node class.")
    node_id = node.define_schema().node_id
    prepare_arguments = getattr(node, "prepare_worker_arguments", None)
    if prepare_arguments is not None:
        kwargs = prepare_arguments(dict(kwargs))
        if not isinstance(kwargs, dict) or not all(isinstance(name, str) for name in kwargs):
            raise TypeError(f"{node_id} worker argument adapter must return a string-keyed dictionary.")
    runtime = current_artifact_runtime()
    inputs, parameters, ticket_inputs = _prepare_request(node_id, kwargs, runtime)
    selected_worker = _select_worker(worker, ticket_inputs, runtime)
    lease = runtime.begin_run()
    published = False
    try:
        response = await _run_worker_with_comfy_interrupt(
            selected_worker,
            lease.staging_path / "request.json",
            operation_id_for_node(node_id),
            inputs=inputs,
            parameters=parameters,
        )
        records = _validate_response_records(node, response, ticket_inputs)
        artifact_specs: dict[str, dict[str, str]] = {}
        for record in records:
            if record["type"] != "artifact":
                continue
            spec = {
                "kind": record["kind"],
                "codec": record["codec"],
                "payload": record["payload"],
            }
            if spec["codec"] in NATIVE_AFFINITY_CODECS:
                spec["worker_affinity"] = _identity_affinity(selected_worker)
            artifact_specs[str(record["name"])] = spec
        tickets = lease.publish(artifact_specs) if artifact_specs else {}
        published = bool(artifact_specs)
        values = []
        for record in records:
            record_type = record["type"]
            if record_type == "artifact":
                values.append(tickets[str(record["name"])])
            elif record_type == "summary":
                values.append(summary_from_metadata(record["value"]))
            elif record_type == "string":
                values.append(record["value"])
            else:
                values.append(ticket_inputs[record["input"]])
        if not published:
            lease.abort()
        return io.NodeOutput(*values)
    except BaseException:
        if not published:
            lease.abort()
        raise


__all__ = [
    "ANNDATA_CODEC",
    "current_artifact_runtime",
    "execute_artifact_node",
    "initialize_artifact_service",
    "operation_id_for_node",
]
