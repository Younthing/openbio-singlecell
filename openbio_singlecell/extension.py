from __future__ import annotations

from comfy_api.latest import ComfyExtension, io

from .nodes_data import DATA_NODE_CLASSES
from .nodes_embedding import EMBEDDING_NODE_CLASSES
from .nodes_input import INPUT_NODE_CLASSES
from .nodes_integration import INTEGRATION_NODE_CLASSES
from .nodes_output import OUTPUT_NODE_CLASSES
from .nodes_preprocess import PREPROCESS_NODE_CLASSES
from .nodes_qc import QC_NODE_CLASSES
from .nodes_results import RESULT_NODE_CLASSES

NODE_CLASSES: list[type[io.ComfyNode]] = [
    *INPUT_NODE_CLASSES,
    *DATA_NODE_CLASSES,
    *QC_NODE_CLASSES,
    *PREPROCESS_NODE_CLASSES,
    *INTEGRATION_NODE_CLASSES,
    *EMBEDDING_NODE_CLASSES,
    *RESULT_NODE_CLASSES,
    *OUTPUT_NODE_CLASSES,
]


class OpenBioSingleCellExtension(ComfyExtension):
    async def get_node_list(self) -> list[type[io.ComfyNode]]:
        return NODE_CLASSES


async def comfy_entrypoint() -> OpenBioSingleCellExtension:
    return OpenBioSingleCellExtension()


__all__ = ["NODE_CLASSES", "OpenBioSingleCellExtension", "comfy_entrypoint"]
