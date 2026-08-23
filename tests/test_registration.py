from __future__ import annotations

import asyncio
import os
import subprocess
import sys
import textwrap
from pathlib import Path

from openbio_singlecell import PLUGIN_VERSION, SCHEMA_VERSION, dependencies
from openbio_singlecell.extension import NODE_CLASSES, OpenBioSingleCellExtension, comfy_entrypoint
from openbio_singlecell.node_types import AnnDataType, SingleCellResultType

EXPECTED_NODE_IDS = {
    "OpenBioSingleCellLoadH5AD",
    "OpenBioSingleCellLoad10xMTX",
    "OpenBioSingleCellLoad10xStudy",
    "OpenBioSingleCellLoad10xH5",
    "OpenBioSingleCellAnnDataSummary",
    "OpenBioSingleCellUseExpressionLayer",
    "OpenBioSingleCellNormalizeGeneNames",
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
    "OpenBioSingleCellMarkerORAAnnotation",
    "OpenBioSingleCellMapClusterAnnotations",
    "OpenBioSingleCellLianaCommunication",
    "OpenBioSingleCellLianaResults",
    "OpenBioSingleCellLianaDotPlot",
    "OpenBioSingleCellCellCycleScore",
    "OpenBioSingleCellDiffusionMap",
    "OpenBioSingleCellPAGA",
    "OpenBioSingleCellDPT",
    "OpenBioSingleCellCassiopeiaLineageQC",
    "OpenBioSingleCellCassiopeiaExpansionTest",
    "OpenBioSingleCellVelocityFilterAndNormalize",
    "OpenBioSingleCellVelocityMoments",
    "OpenBioSingleCellEstimateVelocity",
    "OpenBioSingleCellVelocityGraph",
    "OpenBioSingleCellRecoverDynamics",
    "OpenBioSingleCellVelocityGeneRanking",
    "OpenBioSingleCellVelocityStreamPlot",
    "OpenBioSingleCellInferCNV",
    "OpenBioSingleCellCNVStructure",
    "OpenBioSingleCellPseudobulk",
    "OpenBioSingleCellPseudobulkEdgeR",
    "OpenBioSingleCellDecouplerPseudobulkContrast",
    "OpenBioSingleCellPseudobulkDESeq2",
    "OpenBioSingleCellSCVIDifferentialExpression",
    "OpenBioSingleCellSampleCompositionSummary",
    "OpenBioSingleCellDifferentialCompositionTest",
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
    "OpenBioSingleCellRunPySCENIC",
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
}


def test_package_metadata_is_versioned():
    assert PLUGIN_VERSION == "0.1.0"
    assert SCHEMA_VERSION == 1
    assert AnnDataType.io_type == "OPENBIO_ANNDATA"
    assert SingleCellResultType.io_type == "OPENBIO_SINGLE_CELL_RESULT"


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
    assert {node.GET_SCHEMA().node_id for node in NODE_CLASSES} == EXPECTED_NODE_IDS
    assert asyncio.run(extension.get_node_list()) == NODE_CLASSES


def test_node_outputs_use_only_the_public_analysis_contracts():
    expected_types = {
        "adata": AnnDataType.io_type,
        "result": SingleCellResultType.io_type,
    }
    for node in NODE_CLASSES:
        schema = node.GET_SCHEMA()
        output_names = [output.display_name for output in schema.outputs]
        assert output_names in (["adata"], ["result"], [])
        assert schema.is_output_node is (not output_names)
        for output in schema.outputs:
            assert output.io_type == expected_types[output.display_name]


def test_extension_loads_without_scientific_dependencies():
    comfy_root = next(Path(entry) for entry in sys.path if (Path(entry) / "main.py").is_file())
    script = textwrap.dedent(
        """
        import asyncio
        import importlib.abc
        import sys

        blocked_roots = {"anndata", "matplotlib", "pandas", "scanpy", "scipy"}
        attempted = []

        class BlockScience(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname.partition(".")[0] in blocked_roots:
                    attempted.append(fullname)
                    raise ModuleNotFoundError(f"{fullname} blocked for dependency test")
                return None

        sys.meta_path.insert(0, BlockScience())

        from openbio_singlecell import dependencies
        from openbio_singlecell.extension import NODE_CLASSES, comfy_entrypoint

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
