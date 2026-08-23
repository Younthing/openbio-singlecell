from __future__ import annotations

import logging

from comfy_api.latest import ComfyExtension, io

from . import dependencies
from .nodes_embedding import EMBEDDING_NODE_CLASSES
from .nodes_input import INPUT_NODE_CLASSES
from .nodes_preprocess import PREPROCESS_NODE_CLASSES
from .nodes_qc import QC_NODE_CLASSES
from .nodes_results import RESULT_NODE_CLASSES

NODE_CLASSES: list[type[io.ComfyNode]] = [
    *INPUT_NODE_CLASSES,
    *QC_NODE_CLASSES,
    *PREPROCESS_NODE_CLASSES,
    *EMBEDDING_NODE_CLASSES,
    *RESULT_NODE_CLASSES,
]


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
