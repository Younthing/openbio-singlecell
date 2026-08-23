from __future__ import annotations

import asyncio
import os
import subprocess
import sys
import textwrap
from pathlib import Path

from openbio_singlecell import PLUGIN_VERSION, SCHEMA_VERSION
from openbio_singlecell.extension import NODE_CLASSES, OpenBioSingleCellExtension, comfy_entrypoint


def test_package_metadata_is_versioned():
    assert PLUGIN_VERSION == "0.1.0"
    assert SCHEMA_VERSION == 1


def test_empty_extension_loads():
    extension = asyncio.run(comfy_entrypoint())

    assert isinstance(extension, OpenBioSingleCellExtension)
    assert NODE_CLASSES == []
    assert asyncio.run(extension.get_node_list()) == []


def test_extension_loads_without_scientific_dependencies():
    comfy_root = next(Path(entry) for entry in sys.path if (Path(entry) / "main.py").is_file())
    script = textwrap.dedent(
        """
        import asyncio
        import importlib.abc
        import sys

        class BlockAnndata(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname == "anndata" or fullname.startswith("anndata."):
                    raise ModuleNotFoundError("anndata blocked for dependency test")
                return None

        sys.meta_path.insert(0, BlockAnndata())

        from openbio_singlecell import dependencies
        from openbio_singlecell.extension import NODE_CLASSES, comfy_entrypoint

        assert not dependencies.AVAILABLE
        assert NODE_CLASSES == []

        try:
            dependencies.require_scientific_dependencies()
        except RuntimeError as error:
            message = str(error)
            assert "scientific dependencies are unavailable" in message
            assert "pip install -r" in message
        else:
            raise AssertionError("missing dependencies did not produce an installation error")

        extension = asyncio.run(comfy_entrypoint())
        asyncio.run(extension.on_load())
        assert asyncio.run(extension.get_node_list()) == []
        """
    )
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(Path(__file__).resolve().parents[1]), str(comfy_root), environment.get("PYTHONPATH", "")]
    )

    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=Path(__file__).resolve().parents[1],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
