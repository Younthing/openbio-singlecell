from __future__ import annotations

import asyncio
import os
import signal
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .worker_protocol import (
    CORE_DEPENDENCIES,
    PROBE_OPERATION,
    PROTOCOL_VERSION,
    JSONContainer,
    JSONValue,
    ProtocolError,
    WorkerRequest,
    WorkerResponse,
    read_json,
    response_path_for,
    write_json,
)

OUTPUT_TAIL_BYTES = 64 * 1024
PIPE_DRAIN_TIMEOUT_SECONDS = 5.0
_POSIX_TERMINATION_GRACE_SECONDS = 5.0
WORKER_ENTRY = Path(__file__).with_name("worker_entry.py").resolve()


class WorkerRuntimeError(RuntimeError):
    pass


class WorkerProcessError(WorkerRuntimeError):
    def __init__(self, message: str, *, returncode: int | None, stdout_tail: str, stderr_tail: str) -> None:
        details = [message]
        if stderr_tail:
            details.append(f"stderr tail:\n{stderr_tail}")
        elif stdout_tail:
            details.append(f"stdout tail:\n{stdout_tail}")
        super().__init__("\n".join(details))
        self.returncode = returncode
        self.stdout_tail = stdout_tail
        self.stderr_tail = stderr_tail


class WorkerOperationError(WorkerRuntimeError):
    def __init__(self, operation: str, error_type: str, message: str, *, stdout_tail: str, stderr_tail: str) -> None:
        super().__init__(f"Worker operation {operation!r} failed with {error_type}: {message}")
        self.operation = operation
        self.error_type = error_type
        self.stdout_tail = stdout_tail
        self.stderr_tail = stderr_tail


class WorkerProbeError(WorkerRuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class WorkerIdentity:
    executable: str
    size: int
    mtime_ns: int
    protocol_version: int
    python_version: tuple[int, int, int]
    dependencies: tuple[tuple[str, str], ...]

    @property
    def fingerprint(self) -> tuple[object, ...]:
        return (
            self.executable,
            self.size,
            self.mtime_ns,
            self.protocol_version,
            self.python_version,
            self.dependencies,
        )


def canonical_python_executable(executable: str | os.PathLike[str]) -> tuple[str, os.stat_result]:
    raw_path = os.fspath(executable)
    if not isinstance(raw_path, str):
        raise WorkerProbeError("Python executable path must be a string.")
    if not raw_path.strip():
        raise WorkerProbeError("Python executable path cannot be empty.")
    canonical = os.path.normcase(os.path.realpath(raw_path))
    if not os.path.isabs(raw_path):
        raise WorkerProbeError("Python executable path must be absolute.")
    if not os.path.isfile(canonical):
        raise WorkerProbeError(f"Python executable is not a file: {raw_path}")
    if os.name != "nt" and not os.access(canonical, os.X_OK):
        raise WorkerProbeError(f"Python executable is not executable: {raw_path}")
    return canonical, os.stat(canonical)


async def _drain_tail(stream: asyncio.StreamReader, limit: int) -> bytes:
    tail = bytearray()
    while chunk := await stream.read(64 * 1024):
        tail.extend(chunk)
        if len(tail) > limit:
            del tail[: len(tail) - limit]
    return bytes(tail)


async def _collect_drain_tails(
    stdout_task: asyncio.Task[bytes],
    stderr_task: asyncio.Task[bytes],
) -> tuple[bytes, bytes]:
    tasks = (stdout_task, stderr_task)
    done, pending = await asyncio.wait(tasks, timeout=PIPE_DRAIN_TIMEOUT_SECONDS)
    for task in pending:
        task.cancel()
    if pending:
        await asyncio.gather(*pending, return_exceptions=True)
    return tuple(task.result() if task in done else b"" for task in tasks)


def _decode_tail(value: bytes) -> str:
    return value.decode("utf-8", errors="replace")


def _create_windows_kill_job(pid: int) -> int:
    import ctypes
    from ctypes import wintypes

    class BasicLimitInformation(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_longlong),
            ("PerJobUserTimeLimit", ctypes.c_longlong),
            ("LimitFlags", wintypes.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", wintypes.DWORD),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", wintypes.DWORD),
            ("SchedulingClass", wintypes.DWORD),
        ]

    class IOCounters(ctypes.Structure):
        _fields_ = [
            ("ReadOperationCount", ctypes.c_ulonglong),
            ("WriteOperationCount", ctypes.c_ulonglong),
            ("OtherOperationCount", ctypes.c_ulonglong),
            ("ReadTransferCount", ctypes.c_ulonglong),
            ("WriteTransferCount", ctypes.c_ulonglong),
            ("OtherTransferCount", ctypes.c_ulonglong),
        ]

    class ExtendedLimitInformation(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", BasicLimitInformation),
            ("IoInfo", IOCounters),
            ("ProcessMemoryLimit", ctypes.c_size_t),
            ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t),
            ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CreateJobObjectW.argtypes = (ctypes.c_void_p, wintypes.LPCWSTR)
    kernel32.CreateJobObjectW.restype = wintypes.HANDLE
    kernel32.SetInformationJobObject.argtypes = (
        wintypes.HANDLE,
        ctypes.c_int,
        ctypes.c_void_p,
        wintypes.DWORD,
    )
    kernel32.SetInformationJobObject.restype = wintypes.BOOL
    kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.AssignProcessToJobObject.argtypes = (wintypes.HANDLE, wintypes.HANDLE)
    kernel32.AssignProcessToJobObject.restype = wintypes.BOOL
    kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
    kernel32.CloseHandle.restype = wintypes.BOOL

    job = kernel32.CreateJobObjectW(None, None)
    if not job:
        raise WorkerRuntimeError("Windows could not create a worker process job.") from ctypes.WinError()
    information = ExtendedLimitInformation()
    information.BasicLimitInformation.LimitFlags = 0x00002000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
    try:
        if not kernel32.SetInformationJobObject(
            job,
            9,  # JobObjectExtendedLimitInformation
            ctypes.byref(information),
            ctypes.sizeof(information),
        ):
            raise WorkerRuntimeError("Windows could not configure the worker process job.") from ctypes.WinError()
        process_handle = kernel32.OpenProcess(0x00000101, False, pid)
        if not process_handle:
            raise WorkerRuntimeError("Windows could not open the worker process.") from ctypes.WinError()
        try:
            if not kernel32.AssignProcessToJobObject(job, process_handle):
                raise WorkerRuntimeError("Windows could not assign the worker process job.") from ctypes.WinError()
        finally:
            kernel32.CloseHandle(process_handle)
    except BaseException:
        kernel32.CloseHandle(job)
        raise
    return int(job)


def _terminate_windows_job(job: int) -> bool:
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.TerminateJobObject.argtypes = (wintypes.HANDLE, wintypes.UINT)
    kernel32.TerminateJobObject.restype = wintypes.BOOL
    return bool(kernel32.TerminateJobObject(job, 1))


def _close_windows_job(job: int) -> None:
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
    kernel32.CloseHandle.restype = wintypes.BOOL
    kernel32.CloseHandle(job)


async def _terminate_process_tree(process: asyncio.subprocess.Process, windows_job: int | None = None) -> None:
    if os.name == "nt":
        if windows_job is not None and _terminate_windows_job(windows_job):
            if process.returncode is None:
                await process.wait()
            return
        try:
            taskkill = await asyncio.create_subprocess_exec(
                "taskkill",
                "/PID",
                str(process.pid),
                "/T",
                "/F",
                stdout=asyncio.subprocess.DEVNULL,
                stderr=asyncio.subprocess.DEVNULL,
            )
            await taskkill.wait()
        except OSError:
            if process.returncode is None:
                process.kill()
        if process.returncode is None:
            try:
                await asyncio.wait_for(process.wait(), timeout=5)
            except TimeoutError:
                process.kill()
                await process.wait()
        return

    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        await process.wait()
        return

    deadline = asyncio.get_running_loop().time() + _POSIX_TERMINATION_GRACE_SECONDS
    while True:
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            break
        if asyncio.get_running_loop().time() >= deadline:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            break
        await asyncio.sleep(0.05)
    await process.wait()


async def _invoke_worker(
    executable: str | os.PathLike[str],
    entry_path: str | os.PathLike[str],
    request_path: str | os.PathLike[str],
    request: WorkerRequest,
) -> tuple[WorkerResponse, str, str]:
    canonical_executable, _ = canonical_python_executable(executable)
    entry = Path(entry_path).resolve(strict=True)
    if not entry.is_file():
        raise WorkerRuntimeError(f"Worker entry is not a file: {entry}")
    request_file = Path(request_path).resolve()
    response_file = response_path_for(request_file)
    try:
        response_file.unlink()
    except FileNotFoundError:
        pass
    write_json(request_file, request.to_json())

    process_options: dict[str, object]
    if os.name == "nt":
        process_options = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
    else:
        process_options = {"start_new_session": True}
    process = await asyncio.create_subprocess_exec(
        canonical_executable,
        "-I",
        str(entry),
        "--request",
        str(request_file),
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        **process_options,
    )
    windows_job = None
    if os.name == "nt":
        try:
            windows_job = _create_windows_kill_job(process.pid)
        except BaseException:
            await _terminate_process_tree(process)
            raise
    assert process.stdout is not None
    assert process.stderr is not None
    stdout_task = asyncio.create_task(_drain_tail(process.stdout, OUTPUT_TAIL_BYTES))
    stderr_task = asyncio.create_task(_drain_tail(process.stderr, OUTPUT_TAIL_BYTES))
    try:
        try:
            returncode = await process.wait()
            stdout_bytes, stderr_bytes = await _collect_drain_tails(stdout_task, stderr_task)
        except asyncio.CancelledError:
            await asyncio.shield(_terminate_process_tree(process, windows_job))
            await asyncio.shield(_collect_drain_tails(stdout_task, stderr_task))
            raise
    finally:
        if windows_job is not None:
            _close_windows_job(windows_job)
    stdout_tail = _decode_tail(stdout_bytes)
    stderr_tail = _decode_tail(stderr_bytes)
    if returncode != 0:
        raise WorkerProcessError(
            f"Worker process exited with code {returncode}.",
            returncode=returncode,
            stdout_tail=stdout_tail,
            stderr_tail=stderr_tail,
        )
    if not response_file.is_file():
        raise WorkerProcessError(
            "Worker process exited without a response.",
            returncode=returncode,
            stdout_tail=stdout_tail,
            stderr_tail=stderr_tail,
        )
    response = WorkerResponse.from_json(read_json(response_file))
    if response.request_id != request.request_id:
        raise ProtocolError(
            f"Worker response request ID {response.request_id!r} does not match {request.request_id!r}."
        )
    return response, stdout_tail, stderr_tail


async def probe_python(executable: str | os.PathLike[str]) -> WorkerIdentity:
    canonical, before = canonical_python_executable(executable)
    request = WorkerRequest.create(PROBE_OPERATION, inputs={}, parameters={})
    try:
        with tempfile.TemporaryDirectory(prefix="openbio-worker-probe-") as control_directory:
            response, stdout_tail, stderr_tail = await _invoke_worker(
                canonical,
                WORKER_ENTRY,
                Path(control_directory) / "request.json",
                request,
            )
    except WorkerRuntimeError as error:
        raise WorkerProbeError(f"Python worker probe failed for {canonical}: {error}") from error
    if not response.ok:
        assert response.error is not None
        raise WorkerProbeError(
            f"Python worker probe failed with {response.error['type']}: {response.error['message']}\n{stderr_tail or stdout_tail}"
        )
    assert response.outputs is not None
    outputs = response.outputs
    if not isinstance(outputs, dict):
        raise WorkerProbeError("Python worker probe returned an invalid response object.")
    protocol_version = outputs.get("protocol_version")
    python_version = outputs.get("python_version")
    dependencies = outputs.get("dependencies")
    if outputs.get("isolated") is not True:
        raise WorkerProbeError("Python worker probe did not run in isolated mode.")
    if protocol_version != PROTOCOL_VERSION:
        raise WorkerProbeError(f"Python worker protocol {protocol_version!r} does not match {PROTOCOL_VERSION}.")
    if (
        not isinstance(python_version, list)
        or len(python_version) != 3
        or any(type(part) is not int for part in python_version)
    ):
        raise WorkerProbeError("Python worker probe returned an invalid Python version.")
    version_tuple = tuple(python_version)
    if version_tuple < (3, 12, 0):
        raise WorkerProbeError(f"Python 3.12 or newer is required; worker reported {version_tuple}.")
    if not isinstance(dependencies, dict) or set(dependencies) != set(CORE_DEPENDENCIES):
        raise WorkerProbeError("Python worker probe returned an invalid core dependency set.")
    if any(not isinstance(version, str) or not version for version in dependencies.values()):
        raise WorkerProbeError("Python worker probe returned an invalid dependency version.")
    after = os.stat(canonical)
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise WorkerProbeError("Python executable changed while it was being probed.")
    return WorkerIdentity(
        executable=canonical,
        size=after.st_size,
        mtime_ns=after.st_mtime_ns,
        protocol_version=protocol_version,
        python_version=version_tuple,
        dependencies=tuple(sorted((name, version) for name, version in dependencies.items())),
    )


def _validate_identity(identity: WorkerIdentity) -> None:
    if not isinstance(identity, WorkerIdentity):
        raise TypeError("Expected a probed WorkerIdentity.")
    canonical, stat = canonical_python_executable(identity.executable)
    if canonical != identity.executable or (stat.st_size, stat.st_mtime_ns) != (identity.size, identity.mtime_ns):
        raise WorkerRuntimeError("Python worker executable changed after it was probed.")
    if identity.protocol_version != PROTOCOL_VERSION:
        raise WorkerRuntimeError("Python worker identity has a stale protocol version.")


async def run_worker(
    identity: WorkerIdentity,
    request_path: str | os.PathLike[str],
    operation: str,
    *,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> JSONContainer:
    _validate_identity(identity)
    request = WorkerRequest.create(operation, inputs=inputs, parameters=parameters)
    request_file = Path(request_path).resolve()
    try:
        response, stdout_tail, stderr_tail = await _invoke_worker(
            identity.executable,
            WORKER_ENTRY,
            request_file,
            request,
        )
        if not response.ok:
            assert response.error is not None
            raise WorkerOperationError(
                operation,
                str(response.error["type"]),
                str(response.error["message"]),
                stdout_tail=stdout_tail,
                stderr_tail=stderr_tail,
            )
        assert response.outputs is not None
        return response.outputs
    finally:
        for control_file in (request_file, response_path_for(request_file)):
            try:
                control_file.unlink()
            except FileNotFoundError:
                pass


__all__ = [
    "OUTPUT_TAIL_BYTES",
    "WORKER_ENTRY",
    "WorkerIdentity",
    "WorkerOperationError",
    "WorkerProbeError",
    "WorkerProcessError",
    "WorkerRuntimeError",
    "canonical_python_executable",
    "probe_python",
    "run_worker",
]
