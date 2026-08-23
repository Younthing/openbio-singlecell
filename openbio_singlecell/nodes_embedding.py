from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/dimension-reduction"


class OpenBioSingleCellPCA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPCA",
            display_name="PCA",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_comps", default=50, min=1, max=4096),
                io.Boolean.Input("use_hvg", default=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        n_comps: int = 50,
        use_hvg: bool = True,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        if use_hvg and "highly_variable" not in adata.var:
            raise ValueError(
                "PCA with use_hvg enabled requires var['highly_variable']; run Highly Variable Genes first."
            )
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        dependencies.sc.pp.pca(
            output,
            n_comps=n_comps,
            mask_var="highly_variable" if use_hvg else None,
            svd_solver="arpack",
            random_state=random_seed,
        )
        parameters = {"n_comps": n_comps, "use_hvg": use_hvg, "random_seed": random_seed}
        finish_adata(
            output,
            "pca",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellNeighbors(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellNeighbors",
            display_name="Neighbors",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_neighbors", default=15, min=2, max=4096),
                io.Int.Input("n_pcs", default=50, min=1, max=4096),
                io.Combo.Input(
                    "metric",
                    options=["cosine", "euclidean", "correlation", "manhattan"],
                    default="cosine",
                ),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        n_neighbors: int = 15,
        n_pcs: int = 50,
        metric: str = "cosine",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        dependencies.sc.pp.neighbors(
            output,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            metric=metric,
            random_state=random_seed,
        )
        parameters = {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "metric": metric,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "neighbors",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellUMAP(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellUMAP",
            display_name="UMAP",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("min_dist", default=0.5, min=0.0, max=1.0, step=0.05),
                io.Float.Input("spread", default=1.0, min=0.000001, step=0.1),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        min_dist: float = 0.5,
        spread: float = 1.0,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        dependencies.sc.tl.umap(
            output,
            min_dist=min_dist,
            spread=spread,
            random_state=random_seed,
        )
        parameters = {"min_dist": min_dist, "spread": spread, "random_seed": random_seed}
        finish_adata(
            output,
            "umap",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellLeiden(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLeiden",
            display_name="Leiden Clustering",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("resolution", default=1.0, min=0.000001, step=0.1),
                io.String.Input("key_added", default="leiden"),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        resolution: float = 1.0,
        key_added: str = "leiden",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        if not key_added:
            raise ValueError("Leiden key_added cannot be empty.")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        dependencies.sc.tl.leiden(
            output,
            resolution=resolution,
            key_added=key_added,
            random_state=random_seed,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        parameters = {
            "resolution": resolution,
            "key_added": key_added,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "leiden",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


EMBEDDING_NODE_CLASSES = [
    OpenBioSingleCellPCA,
    OpenBioSingleCellNeighbors,
    OpenBioSingleCellUMAP,
    OpenBioSingleCellLeiden,
]
