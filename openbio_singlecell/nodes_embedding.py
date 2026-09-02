from __future__ import annotations

from comfy_api.latest import io

from .expression_source import _PCA_SPEC
from .node_types import AnnDataType, PlotResultType, analysis_outputs

CATEGORY = "openbio/single-cell/dimension-reduction"
CLUSTERING_CATEGORY = "openbio/single-cell/clustering"
MAX_RANDOM_SEED = 2**31 - 1
SUPPORTED_GRAPH_LAYOUTS = ("fr", "kk", "fa")


class OpenBioSingleCellPCA(io.ComfyNode):
    EXPRESSION_SOURCE = _PCA_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPCA",
            display_name="PCA",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_comps", default=50, min=1),
                io.Boolean.Input("use_hvg", default=True),
                cls.EXPRESSION_SOURCE.input(),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_output_gib", default=2.0, step=0.25, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellPCAVariancePlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPCAVariancePlot",
            display_name="PCA Variance Plot",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_pcs", default=0, min=0),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellNeighbors(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellNeighbors",
            display_name="Neighbors",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("use_rep", default="X_pca"),
                io.Int.Input("n_dimensions", default=0, min=0),
                io.Int.Input("n_neighbors", default=15, min=2),
                io.Combo.Input(
                    "metric", options=["euclidean", "cosine", "correlation", "manhattan"], default="euclidean"
                ),
                io.Combo.Input("method", options=["umap", "gauss", "jaccard"], default="umap", advanced=True),
                io.String.Input("key_added", default="neighbors", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellUMAP(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellUMAP",
            display_name="UMAP",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.Float.Input("min_dist", default=0.5, min=0.0, step=0.05),
                io.Float.Input("spread", default=1.0, step=0.1),
                io.String.Input("key_added", default="X_umap", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


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
                io.Int.Input("n_dimensions", default=0, min=0),
                io.Float.Input("perplexity", default=30.0, step=1.0),
                io.Combo.Input(
                    "metric", options=["euclidean", "cosine", "correlation", "manhattan"], default="euclidean"
                ),
                io.Float.Input("early_exaggeration", default=12.0, step=1.0, advanced=True),
                io.Float.Input("learning_rate", default=1000.0, step=50.0, advanced=True),
                io.String.Input("key_added", default="X_tsne", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellForceDirectedGraph(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellForceDirectedGraph",
            display_name="Force-Directed Graph",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("layout", options=list(SUPPORTED_GRAPH_LAYOUTS), default="fr"),
                io.Combo.Input("init_mode", options=["random", "existing", "paga"], default="random"),
                io.String.Input("init_key", default="X_draw_graph_fr", advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.String.Input("key_suffix", default="", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellLeiden(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLeiden",
            display_name="Leiden Clustering",
            category=CLUSTERING_CATEGORY,
            description=(
                "Partition one named neighbor graph with Leiden and report modularity, cluster sizes, and multi-start "
                "stability."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("resolution", default=1.0, min=0.0, step=0.1),
                io.String.Input("key_added", default="leiden", advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.Int.Input("n_iterations", default=2, advanced=True),
                io.Int.Input("stability_repeats", default=5, min=1, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


EMBEDDING_NODE_CLASSES = [
    OpenBioSingleCellPCA,
    OpenBioSingleCellPCAVariancePlot,
    OpenBioSingleCellNeighbors,
    OpenBioSingleCellUMAP,
    OpenBioSingleCellTSNE,
    OpenBioSingleCellForceDirectedGraph,
    OpenBioSingleCellLeiden,
]

__all__ = [
    "EMBEDDING_NODE_CLASSES",
    "OpenBioSingleCellForceDirectedGraph",
    "OpenBioSingleCellLeiden",
    "OpenBioSingleCellNeighbors",
    "OpenBioSingleCellPCA",
    "OpenBioSingleCellPCAVariancePlot",
    "OpenBioSingleCellTSNE",
    "OpenBioSingleCellUMAP",
]
