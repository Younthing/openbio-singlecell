from __future__ import annotations

import copy
import json
import os
import shutil
import stat
import uuid
import weakref
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

MANIFEST_VERSION = 1
SESSION_VERSION = 1
NATIVE_AFFINITY_CODECS = frozenset({"scvi-native-directory", "cnmf-native-directory"})


def _is_link_or_reparse(path: Path) -> bool:
    try:
        value = path.lstat()
    except FileNotFoundError:
        return False
    attributes = int(getattr(value, "st_file_attributes", 0))
    reparse_flag = int(getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400))
    return stat.S_ISLNK(value.st_mode) or bool(attributes & reparse_flag)


def _remove_run_paths(staging_path: Path, ready_path: Path) -> None:
    for path in (staging_path, ready_path):
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass
        except OSError:
            # Windows can retain file handles briefly; the session scavenger retries later.
            pass


def _default_process_probe(pid: int) -> float | None:
    import psutil

    try:
        return float(psutil.Process(pid).create_time())
    except psutil.NoSuchProcess:
        return None


def _session_identity(path: Path) -> tuple[int, float] | None:
    try:
        payload = json.loads(
            (path / "session.json").read_text(encoding="utf-8"),
            parse_constant=_reject_json_constant,
        )
    except (OSError, UnicodeError, ValueError):
        return None
    if not isinstance(payload, dict) or set(payload) != {"version", "pid", "process_create_time"}:
        return None
    pid = payload["pid"]
    created = payload["process_create_time"]
    if (
        payload["version"] != SESSION_VERSION
        or not isinstance(pid, int)
        or isinstance(pid, bool)
        or pid <= 0
        or not isinstance(created, (int, float))
        or isinstance(created, bool)
        or created <= 0
    ):
        return None
    return pid, float(created)


def _scavenge_dead_sessions(root: Path, process_probe: Callable[[int], float | None]) -> None:
    for session_path in root.glob("*.session"):
        if not session_path.is_dir() or _is_link_or_reparse(session_path):
            continue
        identity = _session_identity(session_path)
        if identity is None:
            continue
        pid, expected_create_time = identity
        try:
            observed_create_time = process_probe(pid)
        except Exception:
            continue
        if observed_create_time is not None and float(observed_create_time) == expected_create_time:
            continue
        try:
            shutil.rmtree(session_path)
        except OSError:
            # Preserve enough identity to retry after Windows releases open handles.
            try:
                (session_path / "session.json").write_text(
                    json.dumps(
                        {
                            "version": SESSION_VERSION,
                            "pid": pid,
                            "process_create_time": expected_create_time,
                        },
                        allow_nan=False,
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                    encoding="utf-8",
                )
            except OSError:
                pass


def _safe_relative_path(value: object) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError("Artifact payload path must be a non-empty POSIX relative path.")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError("Artifact payload path must stay inside its run directory.")
    return path.as_posix()


def _reject_json_constant(value: str):
    raise ValueError(f"Artifact manifest contains invalid JSON constant {value!r}.")


def _payload_inventory(staging_path: Path, payload: object) -> tuple[str, str, list[dict[str, object]]]:
    relative = _safe_relative_path(payload)
    if _is_link_or_reparse(staging_path):
        raise ValueError("Artifact run directory cannot be a link or reparse point.")
    resolved_root = staging_path.resolve(strict=True)
    path = staging_path.joinpath(*PurePosixPath(relative).parts)
    current = staging_path
    for part in PurePosixPath(relative).parts:
        current /= part
        if _is_link_or_reparse(current):
            raise ValueError(f"Artifact payload cannot traverse a link or reparse point: {relative}")
    try:
        path.resolve(strict=True).relative_to(resolved_root)
    except ValueError as error:
        raise ValueError(f"Artifact payload escapes its run directory: {relative}") from error
    if path.is_file():
        return relative, "file", [{"path": relative, "size": path.stat().st_size}]
    if not path.is_dir():
        raise FileNotFoundError(f"Artifact payload does not exist: {relative}")

    files = []
    for candidate in path.rglob("*"):
        if _is_link_or_reparse(candidate):
            raise ValueError(f"Artifact payload cannot contain a link or reparse point: {relative}")
        try:
            candidate.resolve(strict=True).relative_to(resolved_root)
        except ValueError as error:
            raise ValueError(f"Artifact payload escapes its run directory: {relative}") from error
        if candidate.is_file():
            files.append(
                {
                    "path": candidate.relative_to(staging_path).as_posix(),
                    "size": candidate.stat().st_size,
                }
            )
    files.sort(key=lambda item: str(item["path"]))
    return relative, "directory", files


def _load_manifest(path: Path) -> dict[str, object]:
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_json_constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("Artifact manifest is unreadable.") from error
    if not isinstance(manifest, dict) or set(manifest) != {"version", "session", "run", "outputs"}:
        raise ValueError("Artifact manifest schema is invalid.")
    if manifest["version"] != MANIFEST_VERSION:
        raise ValueError("Artifact manifest schema version is unsupported.")
    if not isinstance(manifest["session"], str) or not isinstance(manifest["run"], str):
        raise ValueError("Artifact manifest schema is invalid.")
    outputs = manifest["outputs"]
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError("Artifact manifest schema is invalid.")
    for output, record in outputs.items():
        if not isinstance(output, str) or not output or not isinstance(record, dict):
            raise ValueError("Artifact manifest schema is invalid.")
        record_keys = frozenset(record)
        required_keys = {"kind", "codec", "payload", "payload_type", "files"}
        if record_keys not in {frozenset(required_keys), frozenset(required_keys | {"worker_affinity"})}:
            raise ValueError("Artifact manifest schema is invalid.")
        if not all(isinstance(record[key], str) and record[key] for key in ("kind", "codec", "payload")):
            raise ValueError("Artifact manifest schema is invalid.")
        if record["codec"] in NATIVE_AFFINITY_CODECS and "worker_affinity" not in record:
            raise ValueError("Native artifact manifest requires worker affinity.")
        if "worker_affinity" in record and (
            record["codec"] not in NATIVE_AFFINITY_CODECS
            or not isinstance(record["worker_affinity"], str)
            or not record["worker_affinity"]
        ):
            raise ValueError("Artifact manifest worker affinity is invalid.")
        _safe_relative_path(record["payload"])
        if record["payload_type"] not in {"file", "directory"} or not isinstance(record["files"], list):
            raise ValueError("Artifact manifest schema is invalid.")
        seen: set[str] = set()
        for item in record["files"]:
            if not isinstance(item, dict) or set(item) != {"path", "size"}:
                raise ValueError("Artifact manifest schema is invalid.")
            relative = _safe_relative_path(item["path"])
            size = item["size"]
            if relative in seen or not isinstance(size, int) or isinstance(size, bool) or size < 0:
                raise ValueError("Artifact manifest schema is invalid.")
            seen.add(relative)
    return manifest


@dataclass(frozen=True, slots=True)
class ArtifactTicket:
    session: str
    run: str
    output: str
    kind: str
    codec: str
    _lease: RunLease | None = field(default=None, init=False, repr=False, compare=False, hash=False)
    _token: object | None = field(default=None, init=False, repr=False, compare=False, hash=False)

    def __copy__(self) -> ArtifactTicket:
        return self

    def __deepcopy__(self, memo: dict[int, Any]) -> ArtifactTicket:
        del memo
        return self

    def __reduce_ex__(self, protocol: int):
        del protocol
        raise TypeError("ArtifactTicket cannot be pickled.")


class RunLease:
    def __init__(self, runtime: ArtifactRuntime, run_id: str, staging_path: Path) -> None:
        self._runtime = runtime
        self.run_id = run_id
        self.staging_path = staging_path
        self.ready_path = staging_path.with_name(f"{run_id}.ready")
        self._published = False
        self._aborted = False
        self._token = object()
        self._finalizer = weakref.finalize(self, _remove_run_paths, staging_path, self.ready_path)

    def publish(self, outputs: dict[str, dict[str, str]]) -> dict[str, ArtifactTicket]:
        if self._aborted:
            raise RuntimeError("Artifact run was aborted.")
        if self._published:
            raise RuntimeError("Artifact run was already published.")
        if not isinstance(outputs, dict) or not outputs:
            raise ValueError("An artifact run must publish at least one output.")
        records: dict[str, dict[str, object]] = {}
        for output, spec in outputs.items():
            if not isinstance(output, str) or not output:
                raise ValueError("Artifact output names must be non-empty strings.")
            required_keys = {"kind", "codec", "payload"}
            spec_keys = frozenset(spec) if isinstance(spec, dict) else frozenset()
            if spec_keys not in {frozenset(required_keys), frozenset(required_keys | {"worker_affinity"})}:
                raise ValueError("Artifact output specifications require kind, codec, and payload.")
            if not all(isinstance(spec[key], str) and spec[key] for key in ("kind", "codec", "payload")):
                raise ValueError("Artifact output specification values must be non-empty strings.")
            affinity = spec.get("worker_affinity")
            if spec["codec"] in NATIVE_AFFINITY_CODECS and affinity is None:
                raise ValueError("Native artifact codecs require Worker affinity.")
            if affinity is not None and (
                spec["codec"] not in NATIVE_AFFINITY_CODECS
                or not isinstance(affinity, str)
                or not affinity
            ):
                raise ValueError("Worker affinity is allowed only for native scVI and cNMF directory codecs.")
            payload, payload_type, files = _payload_inventory(self.staging_path, spec["payload"])
            records[output] = {
                "kind": spec["kind"],
                "codec": spec["codec"],
                "payload": payload,
                "payload_type": payload_type,
                "files": files,
            }
            if affinity is not None:
                records[output]["worker_affinity"] = affinity
        manifest = {
            "version": MANIFEST_VERSION,
            "session": self._runtime.session_id,
            "run": self.run_id,
            "outputs": records,
        }
        (self.staging_path / "manifest.json").write_text(
            json.dumps(manifest, allow_nan=False, sort_keys=True, separators=(",", ":")),
            encoding="utf-8",
        )
        os.rename(self.staging_path, self.ready_path)
        self._published = True

        tickets = {}
        for output, record in records.items():
            ticket = ArtifactTicket(
                session=self._runtime.session_id,
                run=self.run_id,
                output=output,
                kind=str(record["kind"]),
                codec=str(record["codec"]),
            )
            object.__setattr__(ticket, "_lease", self)
            object.__setattr__(ticket, "_token", self._token)
            tickets[output] = ticket
        return tickets

    def abort(self) -> None:
        if self._published:
            raise RuntimeError("A published artifact run cannot be aborted.")
        self._finalizer()
        self._aborted = True
        self._runtime._leases.discard(self)


class ArtifactRuntime:
    def __init__(
        self,
        temp_root: str | Path,
        *,
        uuid_factory: Callable[[], object] = uuid.uuid4,
        process_identity: tuple[int, float] | None = None,
        process_probe: Callable[[int], float | None] = _default_process_probe,
    ) -> None:
        self._uuid_factory = uuid_factory
        self._leases: weakref.WeakSet[RunLease] = weakref.WeakSet()
        self._closed = False
        self.root = Path(temp_root) / "openbio-singlecell" / "artifact-runtime"
        self.root.mkdir(parents=True, exist_ok=True)
        _scavenge_dead_sessions(self.root, process_probe)
        if process_identity is None:
            pid = os.getpid()
            process_create_time = process_probe(pid)
            if process_create_time is None:
                raise RuntimeError("The current process identity could not be determined.")
            process_identity = pid, process_create_time
        while True:
            self.session_id = str(self._uuid_factory())
            self.session_path = self.root / f"{self.session_id}.session"
            try:
                self.session_path.mkdir(exist_ok=False)
            except FileExistsError:
                continue
            break
        pid, process_create_time = process_identity
        self._session_metadata = {
            "version": SESSION_VERSION,
            "pid": pid,
            "process_create_time": process_create_time,
        }
        (self.session_path / "session.json").write_text(
            json.dumps(
                self._session_metadata,
                allow_nan=False,
                sort_keys=True,
                separators=(",", ":"),
            ),
            encoding="utf-8",
        )

    def begin_run(self) -> RunLease:
        if self._closed:
            raise RuntimeError("Artifact runtime is closed.")
        while True:
            run_id = str(self._uuid_factory())
            staging_path = self.session_path / f"{run_id}.partial"
            ready_path = staging_path.with_name(f"{run_id}.ready")
            if ready_path.exists():
                continue
            try:
                staging_path.mkdir(exist_ok=False)
            except FileExistsError:
                continue
            if ready_path.exists():
                staging_path.rmdir()
                continue
            lease = RunLease(self, run_id, staging_path)
            self._leases.add(lease)
            return lease

    def close(self) -> bool:
        if self._leases:
            return False
        self._closed = True
        try:
            shutil.rmtree(self.session_path)
        except FileNotFoundError:
            return True
        except OSError:
            try:
                (self.session_path / "session.json").write_text(
                    json.dumps(
                        self._session_metadata,
                        allow_nan=False,
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                    encoding="utf-8",
                )
            except OSError:
                pass
            return False
        return True

    def _validated_record(self, ticket: ArtifactTicket) -> tuple[RunLease, dict[str, object], str]:
        if not isinstance(ticket, ArtifactTicket) or ticket.session != self.session_id:
            raise ValueError("Artifact ticket does not belong to this runtime session.")
        lease = ticket._lease
        if (
            lease is None
            or ticket._token is not lease._token
            or lease._runtime is not self
            or not lease._published
        ):
            raise ValueError("Artifact ticket has no live run lease.")
        if _is_link_or_reparse(lease.ready_path) or not lease.ready_path.is_dir():
            raise ValueError("Artifact ready directory is missing or invalid.")
        manifest_path = lease.ready_path / "manifest.json"
        manifest = _load_manifest(manifest_path)
        if manifest["session"] != ticket.session or manifest["run"] != ticket.run:
            raise ValueError("Artifact manifest identity does not match its ticket.")
        try:
            record = manifest["outputs"][ticket.output]
        except KeyError as error:
            raise ValueError("Artifact manifest does not contain the ticket output.") from error
        if record["kind"] != ticket.kind or record["codec"] != ticket.codec:
            raise ValueError("Artifact manifest type does not match its ticket.")
        payload, payload_type, files = _payload_inventory(lease.ready_path, record["payload"])
        if payload_type != record["payload_type"] or files != record["files"]:
            raise ValueError("Artifact payload size no longer matches its manifest.")
        return lease, record, payload

    def record(self, ticket: ArtifactTicket) -> dict[str, object]:
        _, record, _ = self._validated_record(ticket)
        return copy.deepcopy(record)

    def resolve(self, ticket: ArtifactTicket) -> Path:
        lease, _, payload = self._validated_record(ticket)
        path = lease.ready_path.joinpath(*PurePosixPath(payload).parts)
        return path
