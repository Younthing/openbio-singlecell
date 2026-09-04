from __future__ import annotations

from comfy_api.latest import io

from .expression_source import _CNV_SPEC
from .node_types import AnnDataType, CNVStateType, PlotResultType, TableResultType, analysis_outputs

CATEGORY = "openbio/single-cell/copy-number"


class OpenBioSingleCellInferCNV(io.ComfyNode):
    EXPRESSION_SOURCE = _CNV_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellInferCNV",
            display_name="Infer CNV",
            category=CATEGORY,
            description=(
                "Infer a fingerprinted expression-derived CNV state from an explicit normal reference, "
                "full-gene transformed expression, and named genome assembly."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("reference_key", default="cell_type"),
                io.String.Input("reference_categories", default=""),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("genome_assembly", default=""),
                io.Int.Input("window_size", default=100, min=2, max=10000),
                io.Int.Input("step", default=10, min=1, max=10000, advanced=True),
                io.Float.Input("lfc_clip", default=3.0, min=1e-12, advanced=True),
                io.Float.Input("dynamic_threshold", default=1.5, min=0.0, advanced=True),
                io.String.Input("exclude_chromosomes", default="chrX,chrY,chrM", advanced=True),
                io.String.Input("output_key", default="cnv", advanced=True),
                io.Int.Input("minimum_reference_cells", default=20, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("chunksize", default=5000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("n_jobs", default=1, min=1, max=256, advanced=True),
                io.Float.Input("max_output_gib", default=4.0, min=0.01, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(CNVStateType.Output(display_name="cnv_state")),
        )


class OpenBioSingleCellCNVPCA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNVPCA",
            display_name="CNV PCA",
            category=CATEGORY,
            description="Reduce one validated inferred-CNV window matrix; graph construction remains separate.",
            inputs=[
                CNVStateType.Input("cnv_state"),
                io.Int.Input("n_comps", default=30, min=1, max=4096),
                io.String.Input("output_key", default="X_cnv_pca", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_output_gib", default=2.0, min=0.01, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellCNVScore(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNVScore",
            display_name="CNV Group Score",
            category=CATEGORY,
            description=(
                "Assign and independently verify descriptive mean absolute inferred-CNV scores for an explicit "
                "categorical partition; no tumor cutoff is inferred."
            ),
            inputs=[
                CNVStateType.Input("cnv_state"),
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cnv_leiden"),
                io.String.Input("output_key", default="cnv_score", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(
                AnnDataType.Output(display_name="adata"), TableResultType.Output(display_name="table")
            ),
        )


class OpenBioSingleCellCNVHeatmapPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNVHeatmapPlot",
            display_name="CNV Heatmap Plot",
            category=CATEGORY,
            description=(
                "Read-only heatmap of the verified genome-ordered inferred-CNV window matrix."
            ),
            inputs=[
                CNVStateType.Input("cnv_state"),
                io.String.Input("groupby", default=""),
                io.DynamicCombo.Input(
                    "view",
                    options=[
                        io.DynamicCombo.Option(
                            "group_mean",
                            [io.Int.Input("max_groups", default=50, min=2, max=256)],
                        ),
                        io.DynamicCombo.Option(
                            "cells",
                            [io.Int.Input("max_cells", default=2000, min=2, max=10000)],
                        ),
                    ],
                ),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellCNVScorePlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNVScorePlot",
            display_name="CNV Score Plot",
            category=CATEGORY,
            description="Read-only descriptive group-score view; no tumor cutoff or hypothesis test is inferred.",
            inputs=[
                TableResultType.Input("table"),
                io.DynamicCombo.Input(
                    "view",
                    options=[
                        io.DynamicCombo.Option(
                            "score_bar",
                            [io.Int.Input("max_groups", default=30, min=1, max=100)],
                        ),
                        io.DynamicCombo.Option(
                            "absolute_interval",
                            [io.Int.Input("max_groups", default=30, min=1, max=100)],
                        ),
                    ],
                ),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


CNV_NODE_CLASSES = [
    OpenBioSingleCellInferCNV,
    OpenBioSingleCellCNVPCA,
    OpenBioSingleCellCNVScore,
    OpenBioSingleCellCNVHeatmapPlot,
    OpenBioSingleCellCNVScorePlot,
]


__all__ = [
    "CNV_NODE_CLASSES",
    "OpenBioSingleCellCNVHeatmapPlot",
    "OpenBioSingleCellCNVPCA",
    "OpenBioSingleCellCNVScore",
    "OpenBioSingleCellCNVScorePlot",
    "OpenBioSingleCellInferCNV",
]
