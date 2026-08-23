from __future__ import annotations

import asyncio
import os
import subprocess
import sys
import textwrap
from pathlib import Path

from openbio_singlecell import PLUGIN_VERSION, SCHEMA_VERSION
from openbio_singlecell.extension import NODE_CLASSES, OpenBioSingleCellExtension, comfy_entrypoint
from openbio_singlecell.node_types import AnnDataType, SingleCellResultType

EXPECTED_NODE_IDS = {
    "OpenBioSingleCellLoadH5AD",
    "OpenBioSingleCellLoad10xMTX",
    "OpenBioSingleCellLoad10xH5",
    "OpenBioSingleCellAnnDataSummary",
}


def test_package_metadata_is_versioned():
    assert PLUGIN_VERSION == "0.1.0"
    assert SCHEMA_VERSION == 1
    assert AnnDataType.io_type == "OPENBIO_ANNDATA"
    assert SingleCellResultType.io_type == "OPENBIO_SINGLE_CELL_RESULT"


def test_input_extension_loads():
    extension = asyncio.run(comfy_entrypoint())

    assert isinstance(extension, OpenBioSingleCellExtension)
    assert {node.GET_SCHEMA().node_id for node in NODE_CLASSES} == EXPECTED_NODE_IDS
    assert asyncio.run(extension.get_node_list()) == NODE_CLASSES


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
        assert len(NODE_CLASSES) == 4

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
        assert len(asyncio.run(extension.get_node_list())) == 4
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
