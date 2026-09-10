"""Run the packaged R drivers inside the existing one-shot Worker's process tree."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from .worker_protocol import read_json, write_json

OUTPUT_TAIL_BYTES = 64 * 1024
PROBE_SCRIPT = Path(__file__).with_name("r") / "probe_runtime.R"


class RRuntimeError(RuntimeError):
    def __init__(self, message: str, *, returncode: int | None = None, stdout_tail: str = "", stderr_tail: str = ""):
        diagnostic = stderr_tail or stdout_tail
        super().__init__(f"{message}\n{diagnostic}" if diagnostic else message)
        self.returncode = returncode
        self.stdout_tail = stdout_tail
        self.stderr_tail = stderr_tail


def _child_environment(runtime: dict) -> dict[str, str]:
    environment = os.environ.copy()
    for field, variable in (("r_home", "R_HOME"), ("library_paths", "R_LIBS_USER")):
        if runtime[field]:
            environment[variable] = os.pathsep.join(runtime[field].splitlines())
    if runtime["path_prefix"]:
        environment["PATH"] = os.pathsep.join((*runtime["path_prefix"].splitlines(), environment.get("PATH", "")))
    return environment


def _description_identity(path: Path) -> dict:
    path = path.resolve(strict=True)
    stat = path.stat()
    return {
        "path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def _invoke(runtime: dict, script: Path, request: dict, workdir: Path, *, timeout: float | None = None) -> dict:
    with tempfile.TemporaryDirectory(prefix=".r-", dir=workdir) as private:
        private = Path(private)
        request_path = private / "request.json"
        write_json(request_path, request)
        with (private / "stdout.log").open("w+b") as stdout, (private / "stderr.log").open("w+b") as stderr:
            execution_error = None
            try:
                process = subprocess.run(
                    [runtime["executable"], "--vanilla", str(script.resolve()), str(request_path)],
                    cwd=workdir, env=_child_environment(runtime), stdout=stdout, stderr=stderr, timeout=timeout,
                )
            except (subprocess.TimeoutExpired, OSError) as error:
                execution_error = error
            tails = []
            for stream in (stdout, stderr):
                stream.seek(0, os.SEEK_END)
                stream.seek(max(0, stream.tell() - OUTPUT_TAIL_BYTES))
                tails.append(stream.read().decode("utf-8", errors="replace"))
        stdout_tail, stderr_tail = tails
        if execution_error is not None:
            message = (
                f"R runtime probe timed out after {timeout} seconds."
                if isinstance(execution_error, subprocess.TimeoutExpired)
                else f"Could not start Rscript: {execution_error}"
            )
            raise RRuntimeError(message, stdout_tail=stdout_tail, stderr_tail=stderr_tail) from execution_error
        if process.returncode:
            raise RRuntimeError(
                f"Rscript exited with status {process.returncode}.", returncode=process.returncode,
                stdout_tail=stdout_tail, stderr_tail=stderr_tail,
            )
        sys.stdout.write(stdout_tail)
        sys.stderr.write(stderr_tail)
        summary = read_json(Path(request["output_dir"]) / "summary.json")
        if not isinstance(summary, dict):
            raise RRuntimeError("R driver summary must be a JSON object.")
        return summary


def run_r_script(runtime: dict, script: Path, request: dict, workdir: Path) -> dict:
    return _invoke(validate_r_runtime(runtime), script, request, workdir)


def validate_r_runtime(value: dict) -> dict:
    if not isinstance(value, dict):
        raise ValueError("R runtime must be a JSON object returned by the R Runtime node.")
    for key in ("executable", "r_home", "library_paths", "path_prefix", "r_version"):
        if not isinstance(value.get(key), str):
            raise ValueError(f"R runtime {key} must be a string.")
    executable = Path(value["executable"])
    if not executable.is_absolute() or not executable.is_file():
        raise ValueError("R runtime executable must be an existing absolute file path.")
    stat = executable.stat()
    if (value.get("size"), value.get("mtime_ns")) != (stat.st_size, stat.st_mtime_ns):
        raise ValueError("Rscript executable changed; rerun the R Runtime node.")
    if not isinstance(value.get("packages"), dict):
        raise ValueError("R runtime packages must be a JSON object.")
    for name, package in value["packages"].items():
        if (
            not isinstance(name, str) or not isinstance(package, dict)
            or not isinstance(package.get("version"), str) or not isinstance(package.get("path"), str)
            or not isinstance(package.get("description"), dict)
        ):
            raise ValueError("R runtime package identity is malformed.")
        try:
            current = _description_identity(Path(package["path"]) / "DESCRIPTION")
        except OSError as error:
            raise ValueError(f"R package {name} changed; rerun the R Runtime node: {error}") from error
        if current != package["description"]:
            raise ValueError(f"R package {name} changed; rerun the R Runtime node.")
    return value


def probe_r_runtime(
    rscript: str, *, r_home: str = "", library_paths: str = "", path_prefix: str = "", packages: tuple[str, ...] = (),
) -> dict:
    runtime = {"r_home": r_home, "library_paths": library_paths, "path_prefix": path_prefix}
    executable = shutil.which(rscript, path=_child_environment(runtime).get("PATH"))
    if executable is None:
        raise ValueError(f"Rscript executable was not found: {rscript}")
    executable = os.path.normcase(os.path.realpath(executable))
    stat = Path(executable).stat()
    runtime.update(executable=executable, size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    with tempfile.TemporaryDirectory(prefix="openbio-r-probe-") as directory:
        result = _invoke(
            runtime, PROBE_SCRIPT, {"packages": list(packages), "output_dir": directory}, Path(directory), timeout=60,
        )
    for package in result["packages"].values():
        package["description"] = _description_identity(Path(package["path"]) / "DESCRIPTION")
    return {**runtime, **result}


def runtime_fingerprint(
    rscript: str, *, r_home: str = "", library_paths: str = "", path_prefix: str = "", packages: tuple[str, ...] = (),
) -> str:
    return json.dumps(probe_r_runtime(
        rscript, r_home=r_home, library_paths=library_paths, path_prefix=path_prefix, packages=packages,
    ), sort_keys=True, separators=(",", ":"))
