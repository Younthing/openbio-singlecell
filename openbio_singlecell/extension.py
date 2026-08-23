from __future__ import annotations

from comfy_api.latest import ComfyExtension, io

from .nodes_annotation import ANNOTATION_NODE_CLASSES
from .nodes_cnv import CNV_NODE_CLASSES
from .nodes_communication import COMMUNICATION_NODE_CLASSES
from .nodes_correction import CORRECTION_NODE_CLASSES
from .nodes_data import DATA_NODE_CLASSES
from .nodes_differential import DIFFERENTIAL_NODE_CLASSES
from .nodes_embedding import EMBEDDING_NODE_CLASSES
from .nodes_enrichment import ENRICHMENT_NODE_CLASSES
from .nodes_input import INPUT_NODE_CLASSES
from .nodes_integration import INTEGRATION_NODE_CLASSES
from .nodes_output import OUTPUT_NODE_CLASSES
from .nodes_population import POPULATION_NODE_CLASSES
from .nodes_preprocess import PREPROCESS_NODE_CLASSES
from .nodes_qc import QC_NODE_CLASSES
from .nodes_results import RESULT_NODE_CLASSES
from .nodes_trajectory import TRAJECTORY_NODE_CLASSES

NODE_CLASSES: list[type[io.ComfyNode]] = [
    *INPUT_NODE_CLASSES,
    *DATA_NODE_CLASSES,
    *QC_NODE_CLASSES,
    *CORRECTION_NODE_CLASSES,
    *PREPROCESS_NODE_CLASSES,
    *INTEGRATION_NODE_CLASSES,
    *EMBEDDING_NODE_CLASSES,
    *ANNOTATION_NODE_CLASSES,
    *COMMUNICATION_NODE_CLASSES,
    *TRAJECTORY_NODE_CLASSES,
    *CNV_NODE_CLASSES,
    *DIFFERENTIAL_NODE_CLASSES,
    *ENRICHMENT_NODE_CLASSES,
    *POPULATION_NODE_CLASSES,
    *RESULT_NODE_CLASSES,
    *OUTPUT_NODE_CLASSES,
]


class OpenBioSingleCellExtension(ComfyExtension):
    async def get_node_list(self) -> list[type[io.ComfyNode]]:
        return NODE_CLASSES


async def comfy_entrypoint() -> OpenBioSingleCellExtension:
    return OpenBioSingleCellExtension()


__all__ = ["NODE_CLASSES", "OpenBioSingleCellExtension", "comfy_entrypoint"]
