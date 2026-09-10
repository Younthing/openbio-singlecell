from __future__ import annotations

import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


def test_built_wheel_contains_the_native_r_drivers(tmp_path):
    root = Path(__file__).resolve().parents[1]
    source = tmp_path / "source"
    source.mkdir()
    for name in ("pyproject.toml", "README.md", "LICENSE"):
        shutil.copyfile(root / name, source / name)
    shutil.copytree(root / "openbio_singlecell", source / "openbio_singlecell",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    wheel_dir = tmp_path / "wheel"
    wheel_dir.mkdir()
    subprocess.run(
        [sys.executable, "-c", "import setuptools.build_meta as b; b.build_wheel(" + repr(str(wheel_dir)) + ")"],
        cwd=source, check=True, capture_output=True, text=True,
    )
    wheels = list(wheel_dir.glob("*.whl"))
    assert len(wheels) == 1
    with zipfile.ZipFile(wheels[0]) as wheel:
        for name in ("probe_runtime.R", "monocle2.R", "monocle2_compat.R"):
            packaged = f"openbio_singlecell/r/{name}"
            assert wheel.read(packaged) == (root / packaged).read_bytes()
