from __future__ import annotations

import logging

from comfy_api.latest import ComfyExtension, io

from . import dependencies
from .nodes_input import INPUT_NODE_CLASSES

NODE_CLASSES: list[type[io.ComfyNode]] = [*INPUT_NODE_CLASSES]


class OpenBioSingleCellExtension(ComfyExtension):
    async def on_load(self) -> None:
        if not dependencies.AVAILABLE:
            logging.warning(
                "openbio-singlecell loaded without scientific dependencies. Run: %s",
                dependencies.INSTALL_COMMAND,
            )

    async def get_node_list(self) -> list[type[io.ComfyNode]]:
        return NODE_CLASSES


async def comfy_entrypoint() -> OpenBioSingleCellExtension:
    return OpenBioSingleCellExtension()


__all__ = ["NODE_CLASSES", "OpenBioSingleCellExtension", "comfy_entrypoint"]
