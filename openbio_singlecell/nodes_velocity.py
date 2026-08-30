from __future__ import annotations

from comfy_api.latest import io

from .node_types import (
    AnnDataType,
    PlotResultType,
    TableResultType,
    VelocityStateType,
    analysis_outputs,
)

CATEGORY = "openbio/single-cell/velocity"


class OpenBioSingleCellVelocityFilterAndNormalize(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityFilterAndNormalize",
            display_name="Prepare Velocity Abundances",
            category=CATEGORY,
            description=(
                "Validate raw spliced/unspliced counts, apply explicit shared-gene filters, and normalize only "
                "canonical velocity layers into a fingerprinted staged artifact."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("spliced_layer", default="spliced"),
                io.String.Input("unspliced_layer", default="unspliced"),
                io.Int.Input("min_shared_counts", default=20, min=0, max=2**31 - 1),
                io.Int.Input("min_shared_cells", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Combo.Input(
                    "normalization_target",
                    options=["median_library"],
                    default="median_library",
                    advanced=True,
                ),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(VelocityStateType.Output(display_name="velocity_state")),
        )


class OpenBioSingleCellVelocityMoments(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityMoments",
            display_name="Velocity Moments",
            category=CATEGORY,
            description="Compute Ms/Mu from one exact named graph without hidden PCA or neighbor construction.",
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.Combo.Input(
                    "mode",
                    options=["connectivities", "distances"],
                    default="connectivities",
                    advanced=True,
                ),
                io.Float.Input("max_dense_gib", default=2.0, min=1e-12, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(VelocityStateType.Output(display_name="velocity_state")),
        )


class OpenBioSingleCellEstimateVelocity(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellEstimateVelocity",
            display_name="Estimate RNA Velocity",
            category=CATEGORY,
            description=(
                "Fit one explicit deterministic, stochastic, or dynamical velocity model on a validated stage."
            ),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.Combo.Input(
                    "mode",
                    options=["stochastic", "deterministic", "dynamical"],
                    default="stochastic",
                ),
                io.String.Input("vkey", default="velocity", advanced=True),
                io.Float.Input("min_r2", default=0.01, min=0.0, advanced=True),
                io.Float.Input("min_likelihood", default=0.001, min=0.0, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(VelocityStateType.Output(display_name="velocity_state")),
        )


class OpenBioSingleCellVelocityGraph(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityGraph",
            display_name="Velocity Directed-Correlation Graph",
            category=CATEGORY,
            description=(
                "Convert one verified velocity estimate into positive/negative directed cosine-alignment evidence."
            ),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.String.Input("vkey", default="velocity", advanced=True),
                io.Combo.Input("xkey", options=["Ms"], default="Ms", advanced=True),
                io.Combo.Input(
                    "mode_neighbors",
                    options=["distances", "connectivities"],
                    default="distances",
                    advanced=True,
                ),
                io.Int.Input("n_jobs", default=1, min=1, max=1024, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(VelocityStateType.Output(display_name="velocity_state")),
        )


class OpenBioSingleCellRecoverDynamics(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRecoverDynamics",
            display_name="Recover RNA Velocity Dynamics",
            category=CATEGORY,
            description=("Recover reusable dynamical splicing parameters for one explicit, deterministic gene set."),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.Combo.Input(
                    "gene_selection",
                    options=["velocity_genes", "all"],
                    default="velocity_genes",
                ),
                io.Int.Input("n_top_genes", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("max_iter", default=10, min=1, max=10_000, advanced=True),
                io.Int.Input("n_jobs", default=1, min=1, max=1024, advanced=True),
                io.Float.Input("max_dense_gib", default=2.0, min=1e-12, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(VelocityStateType.Output(display_name="velocity_state")),
        )


class OpenBioSingleCellVelocityGeneRanking(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityGeneRanking",
            display_name="Recovered Dynamics Fit Ranking",
            category=CATEGORY,
            description=(
                "Pure bounded view of verified recovered fits, ranked by fit likelihood without inferential claims."
            ),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.Int.Input("top_n", default=50, min=1, max=2**31 - 1),
                io.Boolean.Input("include_failed", default=False, advanced=True),
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellVelocityStreamPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityStreamPlot",
            display_name="Velocity Stream Plot",
            category=CATEGORY,
            description="Read-only visualization adapter over one verified velocity-graph artifact.",
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.String.Input("basis", default="umap"),
                io.String.Input("color_key", default="leiden"),
                io.Float.Input("density", default=2.0, min=1e-12, advanced=True),
                io.Float.Input("smooth", default=0.5, min=1e-12, advanced=True),
                io.Float.Input("min_mass", default=1.0, min=0.0, advanced=True),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


VELOCITY_NODE_CLASSES = [
    OpenBioSingleCellVelocityFilterAndNormalize,
    OpenBioSingleCellVelocityMoments,
    OpenBioSingleCellEstimateVelocity,
    OpenBioSingleCellVelocityGraph,
    OpenBioSingleCellRecoverDynamics,
    OpenBioSingleCellVelocityGeneRanking,
    OpenBioSingleCellVelocityStreamPlot,
]


__all__ = ["VELOCITY_NODE_CLASSES"]
