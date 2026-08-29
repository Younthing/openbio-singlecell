from __future__ import annotations

import asyncio
import inspect
import os
import subprocess
import sys
import textwrap
from collections import Counter
from pathlib import Path

from comfy_api.latest import io

from openbio_singlecell import PLUGIN_VERSION, SCHEMA_VERSION, dependencies
from openbio_singlecell.extension import (
    NODE_CLASSES,
    RAW_NODE_CLASSES,
    OpenBioSingleCellExtension,
    comfy_entrypoint,
)
from openbio_singlecell.node_types import (
    AnnDataType,
    AugurResultType,
    CassiopeiaCharactersType,
    CassiopeiaTreeType,
    CNMFRunType,
    CNVStateType,
    DGIdbResourceType,
    LianaResultType,
    PlotResultType,
    PseudobulkType,
    SCENICResultArtifactType,
    SCVIModelType,
    SummaryResultType,
    TableResultType,
    TFActivityArtifactType,
    VelocityStateType,
    WorkerType,
)

EXPECTED_NODE_IDS = {
    "OpenBioSingleCellPythonWorker",
    "OpenBioSingleCellLoadH5AD",
    "OpenBioSingleCellLoad10xMTX",
    "OpenBioSingleCellLoad10xStudy",
    "OpenBioSingleCellLoad10xH5",
    "OpenBioSingleCellCoreStudyParameters",
    "OpenBioSingleCellAnnDataSummary",
    "OpenBioSingleCellSnapshotExpression",
    "OpenBioSingleCellSubsetObservations",
    "OpenBioSingleCellMergeObservationAnnotations",
    "OpenBioSingleCellMapGeneIdsFromGTF",
    "OpenBioSingleCellCalculateQC",
    "OpenBioSingleCellFilterCells",
    "OpenBioSingleCellFilterGenes",
    "OpenBioSingleCellQCPlots",
    "OpenBioSingleCellMarkMADOutliers",
    "OpenBioSingleCellScrublet",
    "OpenBioSingleCellFilterDoublets",
    "OpenBioSingleCellNormalizeTotal",
    "OpenBioSingleCellLog1p",
    "OpenBioSingleCellNormalizeToLayer",
    "OpenBioSingleCellPearsonResidualsToLayer",
    "OpenBioSingleCellHighlyVariableGenes",
    "OpenBioSingleCellScale",
    "OpenBioSingleCellCNMFRankSurvey",
    "OpenBioSingleCellCNMF",
    "OpenBioSingleCellHarmonyIntegration",
    "OpenBioSingleCellSCVIIntegration",
    "OpenBioSingleCellLeidenResolutionSweep",
    "OpenBioSingleCellPCA",
    "OpenBioSingleCellNeighbors",
    "OpenBioSingleCellUMAP",
    "OpenBioSingleCellTSNE",
    "OpenBioSingleCellForceDirectedGraph",
    "OpenBioSingleCellLeiden",
    "OpenBioSingleCellCellTypistAnnotation",
    "OpenBioSingleCellMarkerORAEvidence",
    "OpenBioSingleCellMapClusterAnnotations",
    "OpenBioSingleCellLianaCommunication",
    "OpenBioSingleCellLianaDotPlot",
    "OpenBioSingleCellCellCycleScore",
    "OpenBioSingleCellDiffusionMap",
    "OpenBioSingleCellPAGA",
    "OpenBioSingleCellDPT",
    "OpenBioSingleCellCassiopeiaLineageQC",
    "OpenBioSingleCellReconstructCassiopeiaTree",
    "OpenBioSingleCellCassiopeiaExpansionTest",
    "OpenBioSingleCellCassiopeiaPlasticity",
    "OpenBioSingleCellVelocityFilterAndNormalize",
    "OpenBioSingleCellVelocityMoments",
    "OpenBioSingleCellEstimateVelocity",
    "OpenBioSingleCellVelocityGraph",
    "OpenBioSingleCellRecoverDynamics",
    "OpenBioSingleCellVelocityGeneRanking",
    "OpenBioSingleCellVelocityStreamPlot",
    "OpenBioSingleCellInferCNV",
    "OpenBioSingleCellCNVPCA",
    "OpenBioSingleCellCNVScore",
    "OpenBioSingleCellPseudobulk",
    "OpenBioSingleCellPseudobulkEdgeR",
    "OpenBioSingleCellPseudobulkDESeq2",
    "OpenBioSingleCellSCVIDifferentialExpression",
    "OpenBioSingleCellSampleCompositionSummary",
    "OpenBioSingleCellSchistNestedModel",
    "OpenBioSingleCellMiloDifferentialAbundance",
    "OpenBioSingleCellSccodaDifferentialComposition",
    "OpenBioSingleCellTasccodaDifferentialComposition",
    "OpenBioSingleCellAUCellScores",
    "OpenBioSingleCellGSVAScores",
    "OpenBioSingleCellGenePanelScores",
    "OpenBioSingleCellPathwayScoreTTest",
    "OpenBioSingleCellRankedGSEA",
    "OpenBioSingleCellGeneSetOverrepresentation",
    "OpenBioSingleCellDGIdbAnnotation",
    "OpenBioSingleCellDrugScores",
    "OpenBioSingleCellDrugHypergeometric",
    "OpenBioSingleCellDrugGSEA",
    "OpenBioSingleCellCollecTRIULM",
    "OpenBioSingleCellRankTFActivities",
    "OpenBioSingleCellImportPySCENICResults",
    "OpenBioSingleCellSCENICRegulonSpecificity",
    "OpenBioSingleCellSCENICActivityBinarization",
    "OpenBioSingleCellSCENICTFModules",
    "OpenBioSingleCellAugur",
    "OpenBioSingleCellAugurResults",
    "OpenBioSingleCellCellTypeCorrelation",
    "OpenBioSingleCellMarkerGenes",
    "OpenBioSingleCellUMAPPlot",
    "OpenBioSingleCellFilterMarkerGenes",
    "OpenBioSingleCellMarkerExpressionPlot",
    "OpenBioSingleCellPCAMetadataAssociations",
    "OpenBioSingleCellPreviewResult",
    "OpenBioSingleCellSaveH5AD",
    "OpenBioSingleCellExportCSV",
    "OpenBioSingleCellSavePNG",
    "OpenBioSingleCellPersistArtifact",
}

REMOVED_NODE_IDS = {
    "OpenBioSingleCellUseExpressionLayer",
    "OpenBioSingleCellNormalizeGeneNames",
    "OpenBioSingleCellMarkerORAAnnotation",
    "OpenBioSingleCellDecouplerPseudobulkContrast",
    "OpenBioSingleCellDifferentialCompositionTest",
    "OpenBioSingleCellRunPySCENIC",
    "OpenBioSingleCellLianaResults",
    "OpenBioSingleCellCNVStructure",
}


def test_package_metadata_is_versioned():
    assert PLUGIN_VERSION == "0.2.0"
    assert SCHEMA_VERSION == 1
    assert AnnDataType.io_type == "OPENBIO_ANNDATA"
    assert AugurResultType.io_type == "OPENBIO_AUGUR_RESULT"
    assert DGIdbResourceType.io_type == "OPENBIO_DGIDB_RESOURCE"
    assert LianaResultType.io_type == "OPENBIO_LIANA_RESULT"
    assert TableResultType.io_type == "OPENBIO_SINGLE_CELL_TABLE"
    assert TFActivityArtifactType.io_type == "OPENBIO_TF_ACTIVITY"
    assert PlotResultType.io_type == "OPENBIO_SINGLE_CELL_PLOT"
    assert SummaryResultType.io_type == "OPENBIO_SINGLE_CELL_SUMMARY"
    assert SCVIModelType.io_type == "OPENBIO_SCVI_MODEL"
    assert SCENICResultArtifactType.io_type == "OPENBIO_SCENIC_RESULT"
    assert WorkerType.io_type == "OPENBIO_WORKER"
    assert CassiopeiaCharactersType.io_type == "OPENBIO_CASSIOPEIA_CHARACTERS"
    assert CassiopeiaTreeType.io_type == "OPENBIO_CASSIOPEIA_TREE"
    assert CNMFRunType.io_type == "OPENBIO_CNMF_RUN"
    assert CNVStateType.io_type == "OPENBIO_CNV_STATE"
    assert PseudobulkType.io_type == "OPENBIO_SINGLE_CELL_PSEUDOBULK"
    assert VelocityStateType.io_type == "OPENBIO_VELOCITY_STATE"


def test_scientific_dependency_api_is_minimal():
    assert dependencies.__all__ == [
        "ScientificDependencies",
        "require_scientific_dependencies",
    ]
    assert not hasattr(dependencies, "INSTALL_COMMAND")
    assert not hasattr(dependencies, "REQUIREMENTS_PATH")


def test_input_extension_loads():
    extension = asyncio.run(comfy_entrypoint())

    assert isinstance(extension, OpenBioSingleCellExtension)
    node_ids = [node.GET_SCHEMA().node_id for node in NODE_CLASSES]
    schemas = [node.GET_SCHEMA() for node in NODE_CLASSES]
    assert len(node_ids) == 96
    assert len(node_ids) == len(set(node_ids))
    assert set(node_ids) == EXPECTED_NODE_IDS
    assert REMOVED_NODE_IDS.isdisjoint(node_ids)
    assert all(schema.is_deprecated is False for schema in schemas)
    assert Counter(schema.category for schema in schemas) == {
        "openbio/single-cell/annotation": 4,
        "openbio/single-cell/batch-integration": 2,
        "openbio/single-cell/cell-communication": 2,
        "openbio/single-cell/cell-prioritization": 2,
        "openbio/single-cell/clustering": 3,
        "openbio/single-cell/copy-number": 3,
        "openbio/single-cell/correction": 3,
        "openbio/single-cell/data": 4,
        "openbio/single-cell/diagnostics": 3,
        "openbio/single-cell/differential-abundance": 4,
        "openbio/single-cell/differential-expression": 4,
        "openbio/single-cell/dimension-reduction": 5,
        "openbio/single-cell/enrichment": 10,
        "openbio/single-cell/factorization": 2,
        "openbio/single-cell/input": 4,
        "openbio/single-cell/lineage": 4,
        "openbio/single-cell/marker-evidence": 2,
        "openbio/single-cell/output": 5,
        "openbio/single-cell/preprocessing": 6,
        "openbio/single-cell/qc": 4,
        "openbio/single-cell/regulatory": 6,
        "openbio/single-cell/runtime": 1,
        "openbio/single-cell/study": 1,
        "openbio/single-cell/trajectory": 3,
        "openbio/single-cell/velocity": 7,
        "openbio/single-cell/visualization": 2,
    }
    assert asyncio.run(extension.get_node_list()) == NODE_CLASSES


def test_node_categories_preserve_domain_distinctions():
    categories = {node.GET_SCHEMA().node_id: node.GET_SCHEMA().category for node in NODE_CLASSES}

    assert categories["OpenBioSingleCellMarkerGenes"] == "openbio/single-cell/marker-evidence"
    assert categories["OpenBioSingleCellFilterMarkerGenes"] == "openbio/single-cell/marker-evidence"
    assert categories["OpenBioSingleCellLeiden"] == "openbio/single-cell/clustering"
    assert categories["OpenBioSingleCellAnnDataSummary"] == "openbio/single-cell/diagnostics"
    assert categories["OpenBioSingleCellCellCycleScore"] == "openbio/single-cell/annotation"
    assert categories["OpenBioSingleCellCellTypeCorrelation"] == "openbio/single-cell/diagnostics"


def test_scientific_nodes_are_async_and_expose_one_optional_worker_socket():
    scientific_ids = {
        node.define_schema().node_id
        for node in RAW_NODE_CLASSES
        if node.define_schema().category
        not in {"openbio/single-cell/runtime", "openbio/single-cell/study", "openbio/single-cell/output"}
    }
    registered = {node.GET_SCHEMA().node_id: node for node in NODE_CLASSES}

    for node_id in scientific_ids:
        node = registered[node_id]
        schema = node.GET_SCHEMA()
        worker_inputs = [item for item in schema.inputs if item.id == "worker"]
        assert getattr(node, "OPENBIO_WORKER_ADAPTED", False) is True
        assert inspect.iscoroutinefunction(node.execute)
        assert len(worker_inputs) == 1
        assert schema.inputs[-1] is worker_inputs[0]
        assert worker_inputs[0].get_io_type() == WorkerType.io_type
        assert worker_inputs[0].optional is True

    worker_node = registered["OpenBioSingleCellPythonWorker"]
    worker_schema = worker_node.GET_SCHEMA()
    assert worker_schema.display_name == "Python Worker"
    assert [item.id for item in worker_schema.inputs] == ["python"]
    assert [(item.display_name, item.io_type) for item in worker_schema.outputs] == [
        ("worker", WorkerType.io_type)
    ]
    assert inspect.iscoroutinefunction(worker_node.execute)


def test_worker_registry_exactly_matches_adapted_scientific_nodes():
    from openbio_singlecell import worker_operations as _worker_operations  # noqa: F401
    from openbio_singlecell.artifact_service import operation_id_for_node
    from openbio_singlecell.worker_protocol import registered_operation_ids

    expected = {
        operation_id_for_node(node.GET_SCHEMA().node_id)
        for node in NODE_CLASSES
        if getattr(node, "OPENBIO_WORKER_ADAPTED", False)
    }
    actual = {
        operation_id
        for operation_id in registered_operation_ids()
        if operation_id.startswith("openbio.node.")
    }

    assert actual == expected


def test_node_outputs_use_only_their_concrete_public_contracts():
    expected_types = {
        "adata": AnnDataType.io_type,
        "result": AugurResultType.io_type,
        "table": TableResultType.io_type,
        "plot": PlotResultType.io_type,
        "summary": SummaryResultType.io_type,
        "model": SCVIModelType.io_type,
        "characters": CassiopeiaCharactersType.io_type,
        "tree": CassiopeiaTreeType.io_type,
        "run": CNMFRunType.io_type,
        "pseudobulk": PseudobulkType.io_type,
        "k_metrics": TableResultType.io_type,
        "resource": DGIdbResourceType.io_type,
        "cnv_state": CNVStateType.io_type,
        "activities": TFActivityArtifactType.io_type,
        "scenic_result": SCENICResultArtifactType.io_type,
        "velocity_state": VelocityStateType.io_type,
        "worker": WorkerType.io_type,
    }
    expected_multi_outputs = {
        "OpenBioSingleCellCalculateQC": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellFilterCells": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellFilterGenes": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellQCPlots": [
            ("plot", PlotResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMarkMADOutliers": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellScrublet": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellFilterDoublets": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellNormalizeTotal": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellLog1p": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellNormalizeToLayer": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPearsonResidualsToLayer": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellHighlyVariableGenes": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellScale": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCNMFRankSurvey": [
            ("run", CNMFRunType.io_type),
            ("k_metrics", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCNMF": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellAnnDataSummary": [
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSnapshotExpression": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSubsetObservations": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMergeObservationAnnotations": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMapGeneIdsFromGTF": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPCA": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellHarmonyIntegration": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSCVIIntegration": [
            ("adata", AnnDataType.io_type),
            ("model", SCVIModelType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellNeighbors": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellUMAP": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellTSNE": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellForceDirectedGraph": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellLeiden": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSchistNestedModel": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCellCycleScore": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellDiffusionMap": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPAGA": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellDPT": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCassiopeiaLineageQC": [
            ("characters", CassiopeiaCharactersType.io_type),
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellReconstructCassiopeiaTree": [
            ("tree", CassiopeiaTreeType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCassiopeiaExpansionTest": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCassiopeiaPlasticity": [
            ("adata", AnnDataType.io_type),
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellVelocityFilterAndNormalize": [
            ("velocity_state", VelocityStateType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellVelocityMoments": [
            ("velocity_state", VelocityStateType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellEstimateVelocity": [
            ("velocity_state", VelocityStateType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellVelocityGraph": [
            ("velocity_state", VelocityStateType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellRecoverDynamics": [
            ("velocity_state", VelocityStateType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellVelocityGeneRanking": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellVelocityStreamPlot": [
            ("plot", PlotResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPseudobulk": [
            ("pseudobulk", PseudobulkType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPseudobulkEdgeR": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPseudobulkDESeq2": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSCVIDifferentialExpression": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSampleCompositionSummary": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMiloDifferentialAbundance": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSccodaDifferentialComposition": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellTasccodaDifferentialComposition": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellLeidenResolutionSweep": [
            ("adata", AnnDataType.io_type),
            ("resolution_metrics", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPCAMetadataAssociations": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMarkerGenes": [
            ("table", TableResultType.io_type),
            ("universe", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellFilterMarkerGenes": [
            ("table", TableResultType.io_type),
            ("universe", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellUMAPPlot": [
            ("plot", PlotResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMarkerExpressionPlot": [
            ("plot", PlotResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCellTypeCorrelation": [
            ("table", TableResultType.io_type),
            ("plot", PlotResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellAugur": [
            ("result", AugurResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellAugurResults": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCellTypistAnnotation": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMapClusterAnnotations": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellMarkerORAEvidence": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellLianaCommunication": [
            ("result", LianaResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellLianaDotPlot": [
            ("plot", PlotResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellAUCellScores": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellGSVAScores": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellGenePanelScores": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellPathwayScoreTTest": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellRankedGSEA": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellGeneSetOverrepresentation": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellDGIdbAnnotation": [
            ("resource", DGIdbResourceType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellDrugScores": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellDrugHypergeometric": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellDrugGSEA": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellInferCNV": [
            ("cnv_state", CNVStateType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCNVPCA": [
            ("adata", AnnDataType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCNVScore": [
            ("adata", AnnDataType.io_type),
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCollecTRIULM": [
            ("adata", AnnDataType.io_type),
            ("activities", TFActivityArtifactType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellRankTFActivities": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellImportPySCENICResults": [
            ("adata", AnnDataType.io_type),
            ("scenic_result", SCENICResultArtifactType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSCENICRegulonSpecificity": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSCENICActivityBinarization": [
            ("thresholds", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellSCENICTFModules": [
            ("table", TableResultType.io_type),
            ("summary", SummaryResultType.io_type),
            ("code", "STRING"),
        ],
        "OpenBioSingleCellCoreStudyParameters": [
            ("sample_column", "STRING"),
            ("condition_column", "STRING"),
            ("batch_column", "STRING"),
            ("annotation_column", "STRING"),
            ("reference", "STRING"),
            ("comparison", "STRING"),
        ],
    }
    for node in NODE_CLASSES:
        schema = node.GET_SCHEMA()
        actual_outputs = [(output.display_name, output.io_type) for output in schema.outputs]
        expected_outputs = expected_multi_outputs.get(schema.node_id)
        if expected_outputs is not None:
            assert actual_outputs == expected_outputs
        elif actual_outputs:
            assert len(actual_outputs) == 1
            output_name, output_type = actual_outputs[0]
            assert output_type == expected_types[output_name]
        assert schema.is_output_node is (not actual_outputs)


def test_registered_ports_reject_generic_and_legacy_analysis_contracts():
    allowed_openbio_types = {
        AnnDataType.io_type,
        AugurResultType.io_type,
        CassiopeiaCharactersType.io_type,
        CassiopeiaTreeType.io_type,
        CNMFRunType.io_type,
        CNVStateType.io_type,
        DGIdbResourceType.io_type,
        LianaResultType.io_type,
        PlotResultType.io_type,
        PseudobulkType.io_type,
        SCENICResultArtifactType.io_type,
        SCVIModelType.io_type,
        SummaryResultType.io_type,
        TableResultType.io_type,
        TFActivityArtifactType.io_type,
        VelocityStateType.io_type,
        WorkerType.io_type,
    }
    allowed_output_types = allowed_openbio_types | {"STRING"}
    forbidden_types = {
        "*",
        "MODEL",
        "OBJECT",
        "OPENBIO_DATASET",
        "OPENBIO_SC_DATASET",
        "OPENBIO_SC_RESULT",
        "OPENBIO_SAMPLE_SHEET",
        "OPENBIO_SINGLE_CELL_STUDY_DESIGN",
        "OPENBIO_SINGLE_CELL_RESULT",
        "OPENBIO_SCENIC_BINARY",
    }
    object_consumers = {}

    for node in NODE_CLASSES:
        schema = node.GET_SCHEMA()
        for input_ in schema.inputs:
            port_types = set(input_.get_io_type().split(","))
            assert not port_types & forbidden_types
            assert {item for item in port_types if item.startswith("OPENBIO_")} <= allowed_openbio_types
            assert input_.id.lower() not in {"dataset", "dataset_name"}
            if port_types & {
                AugurResultType.io_type,
                SCVIModelType.io_type,
                CassiopeiaCharactersType.io_type,
                CassiopeiaTreeType.io_type,
                CNMFRunType.io_type,
                CNVStateType.io_type,
                DGIdbResourceType.io_type,
                LianaResultType.io_type,
                PseudobulkType.io_type,
                SCENICResultArtifactType.io_type,
                TFActivityArtifactType.io_type,
                VelocityStateType.io_type,
            }:
                object_consumers[(schema.node_id, input_.id)] = port_types

        for output in schema.outputs:
            assert output.io_type not in forbidden_types
            assert output.io_type in allowed_output_types
            assert output.display_name.lower() not in {"dataset", "dataset_name"}

    assert object_consumers == {
        ("OpenBioSingleCellPersistArtifact", "artifact"): {
            AnnDataType.io_type,
            AugurResultType.io_type,
            CassiopeiaCharactersType.io_type,
            CassiopeiaTreeType.io_type,
            CNVStateType.io_type,
            DGIdbResourceType.io_type,
            LianaResultType.io_type,
            PlotResultType.io_type,
            PseudobulkType.io_type,
            SCENICResultArtifactType.io_type,
            TableResultType.io_type,
            TFActivityArtifactType.io_type,
            VelocityStateType.io_type,
        },
        ("OpenBioSingleCellAugurResults", "result"): {AugurResultType.io_type},
        ("OpenBioSingleCellDrugGSEA", "resource"): {DGIdbResourceType.io_type},
        ("OpenBioSingleCellDrugHypergeometric", "resource"): {DGIdbResourceType.io_type},
        ("OpenBioSingleCellDrugScores", "resource"): {DGIdbResourceType.io_type},
        ("OpenBioSingleCellLianaDotPlot", "result"): {LianaResultType.io_type},
        ("OpenBioSingleCellCNVPCA", "cnv_state"): {CNVStateType.io_type},
        ("OpenBioSingleCellCNVScore", "cnv_state"): {CNVStateType.io_type},
        ("OpenBioSingleCellCNMF", "run"): {CNMFRunType.io_type},
        ("OpenBioSingleCellPseudobulkEdgeR", "pseudobulk"): {PseudobulkType.io_type},
        ("OpenBioSingleCellPseudobulkDESeq2", "pseudobulk"): {PseudobulkType.io_type},
        ("OpenBioSingleCellSCVIDifferentialExpression", "model"): {SCVIModelType.io_type},
        ("OpenBioSingleCellRankTFActivities", "activities"): {TFActivityArtifactType.io_type},
        ("OpenBioSingleCellSCENICRegulonSpecificity", "scenic_result"): {
            SCENICResultArtifactType.io_type
        },
        ("OpenBioSingleCellSCENICActivityBinarization", "scenic_result"): {
            SCENICResultArtifactType.io_type
        },
        ("OpenBioSingleCellSCENICTFModules", "scenic_result"): {
            SCENICResultArtifactType.io_type
        },
        ("OpenBioSingleCellReconstructCassiopeiaTree", "characters"): {
            CassiopeiaCharactersType.io_type
        },
        ("OpenBioSingleCellCassiopeiaExpansionTest", "tree"): {CassiopeiaTreeType.io_type},
        ("OpenBioSingleCellCassiopeiaPlasticity", "tree"): {CassiopeiaTreeType.io_type},
        ("OpenBioSingleCellVelocityMoments", "velocity_state"): {VelocityStateType.io_type},
        ("OpenBioSingleCellEstimateVelocity", "velocity_state"): {VelocityStateType.io_type},
        ("OpenBioSingleCellVelocityGraph", "velocity_state"): {VelocityStateType.io_type},
        ("OpenBioSingleCellRecoverDynamics", "velocity_state"): {VelocityStateType.io_type},
        ("OpenBioSingleCellVelocityGeneRanking", "velocity_state"): {VelocityStateType.io_type},
        ("OpenBioSingleCellVelocityStreamPlot", "velocity_state"): {VelocityStateType.io_type},
    }


def test_parameter_consumers_keep_native_string_widget_inputs():
    expected_inputs = {
        "OpenBioSingleCellSampleCompositionSummary": {"sample_key", "condition_key", "annotation_key"},
        "OpenBioSingleCellPseudobulkEdgeR": {
            "population",
            "reference_condition",
            "comparison_condition",
            "categorical_covariate_keys",
            "continuous_covariate_keys",
        },
        "OpenBioSingleCellPseudobulkDESeq2": {
            "population",
            "reference_condition",
            "comparison_condition",
            "categorical_covariate_keys",
            "continuous_covariate_keys",
        },
        "OpenBioSingleCellAugur": {
            "sample_key",
            "population_key",
            "condition_key",
            "control",
            "treatment",
            "technical_batch_key",
        },
        "OpenBioSingleCellRankTFActivities": {"annotation_key", "reference"},
    }
    schemas = {node.GET_SCHEMA().node_id: node.GET_SCHEMA() for node in NODE_CLASSES}

    for node_id, input_ids in expected_inputs.items():
        inputs = {input_.id: input_ for input_ in schemas[node_id].inputs}
        for input_id in input_ids:
            input_ = inputs[input_id]
            assert isinstance(input_, io.WidgetInput)
            assert input_.get_io_type() == "STRING"
            assert input_.socketless is not True
            assert input_.force_input is not True


def test_extension_loads_without_scientific_dependencies():
    comfy_root = next(Path(entry) for entry in sys.path if (Path(entry) / "main.py").is_file())
    script = textwrap.dedent(
        """
        import asyncio
        import importlib.abc
        import sys

        blocked_roots = {
            "anndata",
            "cassiopeia",
            "decoupler",
            "loompy",
            "matplotlib",
            "omicverse",
            "pandas",
            "pertpy",
            "pyscenic",
            "scanpy",
            "schist",
            "scipy",
            "scvi",
            "graph_tool",
        }
        attempted = []

        class BlockScience(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname.partition(".")[0] in blocked_roots:
                    attempted.append(fullname)
                    raise ModuleNotFoundError(f"{fullname} blocked for dependency test")
                return None

        sys.meta_path.insert(0, BlockScience())

        from comfy.cli_args import args
        from openbio_singlecell import dependencies
        from openbio_singlecell.extension import NODE_CLASSES, comfy_entrypoint

        args.cache_classic = True
        args.cache_none = False
        args.cache_lru = 0
        args.cache_ram = []
        assert NODE_CLASSES
        node_count = len(NODE_CLASSES)
        extension = asyncio.run(comfy_entrypoint())
        asyncio.run(extension.on_load())
        assert len(asyncio.run(extension.get_node_list())) == node_count
        assert attempted == []

        try:
            dependencies.require_scientific_dependencies()
        except RuntimeError as error:
            message = str(error)
            assert "scientific dependencies are unavailable" in message
            assert "pip install -r" in message
        else:
            raise AssertionError("missing dependencies did not produce an installation error")
        assert attempted
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
