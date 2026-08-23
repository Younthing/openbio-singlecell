from __future__ import annotations

import asyncio

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
