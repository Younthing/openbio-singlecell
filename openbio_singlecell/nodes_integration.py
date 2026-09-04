from __future__ import annotations

from comfy_api.latest import io

from .expression_source import _SCVI_SPEC
from .node_types import AnnDataType, PlotResultType, SCVIModelType, TableResultType, analysis_outputs

INTEGRATION_CATEGORY = "openbio/single-cell/batch-integration"
CLUSTERING_CATEGORY = "openbio/single-cell/clustering"


class OpenBioSingleCellHarmonyIntegration(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHarmonyIntegration",
            display_name="Harmony Integration",
            category=INTEGRATION_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("technical_batch_keys", default="batch"),
                io.String.Input("basis", default="X_pca"),
                io.DynamicCombo.Input(
                    "theta",
                    options=[
                        io.DynamicCombo.Option("automatic", []),
                        io.DynamicCombo.Option(
                            "custom",
                            [io.Float.Input("theta_value", default=2.0, min=0.0)],
                        ),
                    ],
                ),
                io.Float.Input("ridge_penalty", default=-1.0, min=-1.0, advanced=True),
                io.Float.Input("sigma", default=0.1, advanced=True),
                io.Int.Input("n_clusters", default=0, min=0, advanced=True),
                io.Float.Input("tau", default=0.0, min=0.0, advanced=True),
                io.String.Input("adjusted_basis", default="X_pca_harmony", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("max_iter_harmony", default=10, min=1, advanced=True),
                io.Int.Input("max_iter_kmeans", default=4, min=1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellSCVIIntegration(io.ComfyNode):
    EXPRESSION_SOURCE = _SCVI_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCVIIntegration",
            display_name="scVI Integration",
            category=INTEGRATION_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("technical_batch_key", default="batch"),
                io.String.Input("categorical_covariates", default=""),
                io.String.Input("continuous_covariates", default=""),
                io.Int.Input("n_latent", default=10, min=1),
                io.Combo.Input("gene_likelihood", options=["zinb", "nb", "poisson"], default="zinb"),
                io.Int.Input("n_layers", default=1, min=1, advanced=True),
                io.Combo.Input("dispersion", options=["gene", "gene-batch"], default="gene", advanced=True),
                io.Float.Input("dropout_rate", default=0.1, min=0.0, max=1.0, step=0.05, advanced=True),
                io.DynamicCombo.Input(
                    "epochs",
                    options=[
                        io.DynamicCombo.Option("automatic", []),
                        io.DynamicCombo.Option("fixed", [io.Int.Input("max_epochs", default=400, min=1)]),
                    ],
                ),
                io.Boolean.Input("early_stopping", default=False, advanced=True),
                io.Float.Input("train_size", default=0.9, min=0.0, max=1.0, step=0.05, advanced=True),
                io.Int.Input("batch_size", default=128, min=1, advanced=True),
                io.String.Input("size_factor_key", default="", advanced=True),
                io.Combo.Input("accelerator", options=["auto", "cpu", "gpu", "mps"], default="auto", advanced=True),
                io.String.Input("output_key", default="X_scVI", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(
                AnnDataType.Output(display_name="adata"), SCVIModelType.Output(display_name="model")
            ),
        )


class OpenBioSingleCellHarmonyConvergencePlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHarmonyConvergencePlot",
            display_name="Harmony Convergence Plot",
            category=INTEGRATION_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("adjusted_basis", default="X_pca_harmony"),
                io.Int.Input("max_image_pixels", default=40_000_000, min=1, advanced=True),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellSCVITrainingPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCVITrainingPlot",
            display_name="scVI Training Plot",
            category=INTEGRATION_CATEGORY,
            inputs=[
                SCVIModelType.Input("model"),
                io.Int.Input("max_image_pixels", default=40_000_000, min=1, advanced=True),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


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
                io.Int.Input("n_iterations", default=2, advanced=True),
                io.Int.Input("stability_repeats", default=5, min=1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(
                AnnDataType.Output(display_name="adata"),
                TableResultType.Output(display_name="resolution_metrics"),
            ),
        )


class OpenBioSingleCellLeidenResolutionSweepPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLeidenResolutionSweepPlot",
            display_name="Leiden Resolution Sweep Plot",
            category=CLUSTERING_CATEGORY,
            description="Read-only quality and cluster-size diagnostics across stored Leiden resolutions.",
            inputs=[
                TableResultType.Input("resolution_metrics"),
                io.DynamicCombo.Input(
                    "view",
                    options=[
                        io.DynamicCombo.Option("quality_curves", []),
                        io.DynamicCombo.Option("cluster_sizes", []),
                    ],
                ),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


INTEGRATION_NODE_CLASSES = [
    OpenBioSingleCellHarmonyIntegration,
    OpenBioSingleCellHarmonyConvergencePlot,
    OpenBioSingleCellSCVIIntegration,
    OpenBioSingleCellSCVITrainingPlot,
    OpenBioSingleCellLeidenResolutionSweep,
    OpenBioSingleCellLeidenResolutionSweepPlot,
]

__all__ = [
    "INTEGRATION_NODE_CLASSES",
    "OpenBioSingleCellHarmonyConvergencePlot",
    "OpenBioSingleCellHarmonyIntegration",
    "OpenBioSingleCellLeidenResolutionSweep",
    "OpenBioSingleCellLeidenResolutionSweepPlot",
    "OpenBioSingleCellSCVIIntegration",
    "OpenBioSingleCellSCVITrainingPlot",
]
