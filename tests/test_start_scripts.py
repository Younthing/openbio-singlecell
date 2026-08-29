from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _fake_checkout(tmp_path: Path) -> tuple[Path, Path]:
    comfy = tmp_path / "ComfyUI"
    frontend = tmp_path / "ComfyUI_frontend"
    comfy.mkdir()
    (frontend / "dist").mkdir(parents=True)
    (frontend / "dist" / "index.html").write_text("ok", encoding="utf-8")
    (comfy / "main.py").write_text(
        "import json, sys\nprint('OPENBIO_ARGS=' + json.dumps(sys.argv[1:]))\n",
        encoding="utf-8",
    )
    return comfy, frontend


def _reported_args(completed: subprocess.CompletedProcess[str]) -> list[str]:
    assert completed.returncode == 0, completed.stdout + completed.stderr
    line = next(line for line in completed.stdout.splitlines() if line.startswith("OPENBIO_ARGS="))
    return json.loads(line.removeprefix("OPENBIO_ARGS="))


@pytest.mark.parametrize(
    "explicit",
    [[], ["--cache-classic"], ["--cache-none"], ["--cache-lru", "3"], ["--cache-ram=2"], ["--high-ram"]],
)
def test_powershell_launcher_defaults_cache_only_when_unspecified(tmp_path, explicit):
    shell = shutil.which("pwsh") or shutil.which("powershell")
    if shell is None:
        pytest.skip("PowerShell is unavailable")
    comfy, frontend = _fake_checkout(tmp_path)

    completed = subprocess.run(
        [
            shell,
            "-NoProfile",
            "-File",
            str(ROOT / "scripts" / "start.ps1"),
            "-Python",
            sys.executable,
            "-ComfyRoot",
            str(comfy),
            "-FrontendRoot",
            str(frontend),
            *explicit,
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    args = _reported_args(completed)
    if explicit:
        assert "--cache-classic" in args if explicit == ["--cache-classic"] else "--cache-classic" not in args
        assert args[-len(explicit) :] == explicit
    else:
        assert args[-1] == "--cache-classic"


def _posix_path(path: Path) -> str:
    if os.name != "nt":
        return str(path)
    return f"/{path.drive[0].lower()}{path.as_posix()[2:]}"


def _find_sh() -> str | None:
    found = shutil.which("sh")
    if found or os.name != "nt":
        return found
    git_sh = Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "Git" / "bin" / "sh.exe"
    return str(git_sh) if git_sh.is_file() else None


@pytest.mark.parametrize(
    "explicit",
    [[], ["--cache-classic"], ["--cache-none"], ["--cache-lru=3"], ["--cache-ram=2"], ["--high-ram"]],
)
def test_posix_launcher_defaults_cache_only_when_unspecified(tmp_path, explicit):
    shell = _find_sh()
    if shell is None:
        pytest.skip("POSIX shell is unavailable")
    comfy, frontend = _fake_checkout(tmp_path)
    environment = os.environ.copy()
    environment.update(
        OPENBIO_PYTHON=_posix_path(Path(sys.executable)),
        OPENBIO_FRONTEND_ROOT=_posix_path(frontend),
    )

    completed = subprocess.run(
        [shell, _posix_path(ROOT / "scripts" / "start.sh"), "--comfy-root", _posix_path(comfy), *explicit],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    args = _reported_args(completed)
    if explicit:
        assert "--cache-classic" in args if explicit == ["--cache-classic"] else "--cache-classic" not in args
        assert args[-len(explicit) :] == explicit
    else:
        assert args[-1] == "--cache-classic"
