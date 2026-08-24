from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/dimension-reduction"


class OpenBioSingleCellPCA(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(description="PCA source", layer_default="log1p_norm")

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
                cls.EXPRESSION_SOURCE.input(),
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
        source: DynamicExpressionSource | None = None,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if use_hvg and "highly_variable" not in adata.var:
            raise ValueError(
                "PCA with use_hvg enabled requires var['highly_variable']; run Highly Variable Genes first."
            )
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.pp.pca(
            output,
            n_comps=n_comps,
            layer=expression.scanpy_layer,
            mask_var="highly_variable" if use_hvg else None,
            svd_solver="arpack",
            random_state=random_seed,
        )
        parameters = {
            "n_comps": n_comps,
            "use_hvg": use_hvg,
            **expression.parameters(),
            "random_seed": random_seed,
        }
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
                io.String.Input("use_rep", default=""),
                io.String.Input("key_added", default="", advanced=True),
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
        use_rep: str = "",
        key_added: str = "",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.pp.neighbors(
            output,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            use_rep=use_rep or None,
            metric=metric,
            random_state=random_seed,
            key_added=key_added or None,
        )
        parameters = {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "metric": metric,
            "use_rep": use_rep,
            "key_added": key_added,
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
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
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
        neighbors_key: str = "neighbors",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.tl.umap(
            output,
            min_dist=min_dist,
            spread=spread,
            neighbors_key=neighbors_key,
            random_state=random_seed,
        )
        parameters = {
            "min_dist": min_dist,
            "spread": spread,
            "neighbors_key": neighbors_key,
            "random_seed": random_seed,
        }
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


class OpenBioSingleCellTSNE(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellTSNE",
            display_name="t-SNE",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("use_rep", default="X_pca"),
                io.Int.Input("n_pcs", default=50, min=1, max=4096),
                io.Float.Input("perplexity", default=30.0, min=1.0, step=1.0),
                io.String.Input("key_added", default="X_tsne", advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        use_rep: str = "X_pca",
        n_pcs: int = 50,
        perplexity: float = 30.0,
        key_added: str = "X_tsne",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if use_rep and use_rep not in adata.obsm and use_rep != "X":
            raise ValueError(f"t-SNE representation not found in obsm: {use_rep!r}")
        key_added = key_added.strip()
        if not key_added:
            raise ValueError("t-SNE key_added cannot be empty.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.tl.tsne(
            output,
            n_pcs=n_pcs,
            use_rep=use_rep or None,
            perplexity=perplexity,
            random_state=random_seed,
            key_added=key_added,
        )
        parameters = {
            "use_rep": use_rep,
            "n_pcs": n_pcs,
            "perplexity": perplexity,
            "key_added": key_added,
            "random_seed": random_seed,
        }
        finish_adata(output, "tsne", parameters, cells, genes, started_at, random_seed=random_seed)
        return io.NodeOutput(output)


class OpenBioSingleCellForceDirectedGraph(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellForceDirectedGraph",
            display_name="Force-Directed Graph",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("layout", options=["fa", "fr", "kk", "drl", "lgl", "rt"], default="fa"),
                io.String.Input("init_pos", default="paga"),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.String.Input("key_suffix", default="", advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        layout: str = "fa",
        init_pos: str = "paga",
        neighbors_key: str = "neighbors",
        key_suffix: str = "",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.tl.draw_graph(
            output,
            layout=layout,
            init_pos=init_pos or None,
            neighbors_key=neighbors_key or None,
            key_added_ext=key_suffix or None,
            random_state=random_seed,
        )
        parameters = {
            "layout": layout,
            "init_pos": init_pos,
            "neighbors_key": neighbors_key,
            "key_suffix": key_suffix,
            "random_seed": random_seed,
        }
        finish_adata(
            output,
            "force_directed_graph",
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
                io.String.Input("key_added", default="leiden", advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
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
        neighbors_key: str = "neighbors",
        random_seed: int = 0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if not key_added:
            raise ValueError("Leiden key_added cannot be empty.")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.tl.leiden(
            output,
            resolution=resolution,
            key_added=key_added,
            neighbors_key=neighbors_key,
            random_state=random_seed,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        parameters = {
            "resolution": resolution,
            "key_added": key_added,
            "neighbors_key": neighbors_key,
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
    OpenBioSingleCellTSNE,
    OpenBioSingleCellForceDirectedGraph,
    OpenBioSingleCellLeiden,
]
