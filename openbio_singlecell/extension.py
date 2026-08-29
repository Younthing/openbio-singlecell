from __future__ import annotations

import sys

import folder_paths
from comfy.cli_args import args as comfy_args
from comfy_api.latest import ComfyExtension, io

from .artifact_service import initialize_artifact_service
from .cache_policy import require_file_artifact_cache
from .node_adapter import adapt_scientific_node
from .nodes_abundance import ABUNDANCE_NODE_CLASSES
from .nodes_annotation import ANNOTATION_NODE_CLASSES
from .nodes_cnv import CNV_NODE_CLASSES
from .nodes_communication import COMMUNICATION_NODE_CLASSES
from .nodes_correction import CORRECTION_NODE_CLASSES
from .nodes_data import DATA_NODE_CLASSES
from .nodes_differential import DIFFERENTIAL_NODE_CLASSES
from .nodes_embedding import EMBEDDING_NODE_CLASSES
from .nodes_enrichment import ENRICHMENT_NODE_CLASSES
from .nodes_factorization import FACTORIZATION_NODE_CLASSES
from .nodes_input import INPUT_NODE_CLASSES
from .nodes_integration import INTEGRATION_NODE_CLASSES
from .nodes_lineage import LINEAGE_NODE_CLASSES
from .nodes_output import OUTPUT_NODE_CLASSES
from .nodes_population import POPULATION_NODE_CLASSES
from .nodes_preprocess import PREPROCESS_NODE_CLASSES
from .nodes_qc import QC_NODE_CLASSES
from .nodes_regulatory import REGULATORY_NODE_CLASSES
from .nodes_results import RESULT_NODE_CLASSES
from .nodes_schist import SCHIST_NODE_CLASSES
from .nodes_study import STUDY_NODE_CLASSES
from .nodes_trajectory import TRAJECTORY_NODE_CLASSES
from .nodes_velocity import VELOCITY_NODE_CLASSES
from .nodes_worker import WORKER_NODE_CLASSES

RAW_NODE_CLASSES: list[type[io.ComfyNode]] = [
    *WORKER_NODE_CLASSES,
    *INPUT_NODE_CLASSES,
    *STUDY_NODE_CLASSES,
    *DATA_NODE_CLASSES,
    *QC_NODE_CLASSES,
    *CORRECTION_NODE_CLASSES,
    *PREPROCESS_NODE_CLASSES,
    *FACTORIZATION_NODE_CLASSES,
    *INTEGRATION_NODE_CLASSES,
    *EMBEDDING_NODE_CLASSES,
    *SCHIST_NODE_CLASSES,
    *ANNOTATION_NODE_CLASSES,
    *COMMUNICATION_NODE_CLASSES,
    *TRAJECTORY_NODE_CLASSES,
    *LINEAGE_NODE_CLASSES,
    *VELOCITY_NODE_CLASSES,
    *CNV_NODE_CLASSES,
    *DIFFERENTIAL_NODE_CLASSES,
    *ABUNDANCE_NODE_CLASSES,
    *ENRICHMENT_NODE_CLASSES,
    *REGULATORY_NODE_CLASSES,
    *POPULATION_NODE_CLASSES,
    *RESULT_NODE_CLASSES,
    *OUTPUT_NODE_CLASSES,
]

_MAIN_PROCESS_CLASSES = {*WORKER_NODE_CLASSES, *STUDY_NODE_CLASSES, *OUTPUT_NODE_CLASSES}
NODE_CLASSES = [
    node if node in _MAIN_PROCESS_CLASSES else adapt_scientific_node(node)
    for node in RAW_NODE_CLASSES
]


class OpenBioSingleCellExtension(ComfyExtension):
    async def on_load(self) -> None:
        require_file_artifact_cache(comfy_args)
        await initialize_artifact_service(folder_paths.get_temp_directory(), sys.executable)

    async def get_node_list(self) -> list[type[io.ComfyNode]]:
        return NODE_CLASSES


async def comfy_entrypoint() -> OpenBioSingleCellExtension:
    return OpenBioSingleCellExtension()


__all__ = ["NODE_CLASSES", "RAW_NODE_CLASSES", "OpenBioSingleCellExtension", "comfy_entrypoint"]
