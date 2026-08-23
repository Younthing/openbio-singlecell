from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


INTEGRATION_CATEGORY = "openbio/single-cell/batch-integration"
CLUSTERING_CATEGORY = "openbio/single-cell/clustering"


def _parse_resolutions(value: str) -> list[float]:
    resolutions = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        resolution = float(item)
        if resolution <= 0:
            raise ValueError("Leiden resolutions must be greater than zero.")
        if resolution not in resolutions:
            resolutions.append(resolution)
    if not resolutions:
        raise ValueError("At least one Leiden resolution is required.")
    return resolutions


def _resolution_key(prefix: str, resolution: float) -> str:
    suffix = f"{resolution:g}".replace("-", "m").replace(".", "_")
    return f"{prefix}_{suffix}"


class OpenBioSingleCellHarmonyIntegration(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHarmonyIntegration",
            display_name="Harmony Integration",
            category=INTEGRATION_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("batch_key", default="batch"),
                io.String.Input("basis", default="X_pca"),
                io.String.Input("adjusted_basis", default="X_pca_harmony", advanced=True),
                io.Int.Input("max_iter_harmony", default=10, min=1, max=1000, advanced=True),
                io.Int.Input("max_iter_kmeans", default=100, min=1, max=1000, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        batch_key: str = "batch",
        basis: str = "X_pca",
        adjusted_basis: str = "X_pca_harmony",
        max_iter_harmony: int = 10,
        max_iter_kmeans: int = 100,
    ) -> io.NodeOutput:
        if batch_key not in adata.obs:
            raise ValueError(f"Harmony batch column not found in obs: {batch_key!r}")
        if basis not in adata.obsm:
            raise ValueError(f"Harmony basis not found in obsm: {basis!r}")
        try:
            import scanpy.external as sce
        except ImportError as error:
            raise RuntimeError("Harmony Integration requires scanpy external dependencies and harmonypy.") from error

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        sce.pp.harmony_integrate(
            output,
            key=batch_key,
            basis=basis,
            adjusted_basis=adjusted_basis,
            max_iter_harmony=max_iter_harmony,
            max_iter_kmeans=max_iter_kmeans,
        )
        parameters = {
            "batch_key": batch_key,
            "basis": basis,
            "adjusted_basis": adjusted_basis,
            "max_iter_harmony": max_iter_harmony,
            "max_iter_kmeans": max_iter_kmeans,
        }
        finish_adata(output, "harmony_integration", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellLeidenResolutionSweep(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLeidenResolutionSweep",
            display_name="Leiden Resolution Sweep",
            category=CLUSTERING_CATEGORY,
            description="Run Leiden clustering for a comma-separated list of resolutions.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("resolutions", default="0.25,0.5,1.0,2.0"),
                io.String.Input("key_prefix", default="leiden", advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        resolutions: str = "0.25,0.5,1.0,2.0",
        key_prefix: str = "leiden",
        neighbors_key: str = "neighbors",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if not key_prefix.strip():
            raise ValueError("Leiden key_prefix cannot be empty.")
        values = _parse_resolutions(resolutions)
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        keys = []
        for resolution in values:
            key = _resolution_key(key_prefix.strip(), resolution)
            science.sc.tl.leiden(
                output,
                resolution=resolution,
                key_added=key,
                neighbors_key=neighbors_key,
                random_state=random_seed,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            keys.append(key)

        parameters = {
            "resolutions": values,
            "keys": keys,
            "neighbors_key": neighbors_key,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "leiden_resolution_sweep",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


INTEGRATION_NODE_CLASSES = [
    OpenBioSingleCellHarmonyIntegration,
    OpenBioSingleCellLeidenResolutionSweep,
]
