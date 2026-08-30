from __future__ import annotations

from comfy_api.latest import io

from .expression_source import _SCRUBLET_SPEC
from .node_types import AnnDataType, analysis_outputs

CATEGORY = "openbio/single-cell/correction"
MAX_RANDOM_SEED = 2**31 - 1


class OpenBioSingleCellMarkMADOutliers(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkMADOutliers",
            display_name="Mark MAD Outliers",
            category=CATEGORY,
            description="Mark a union of metric outliers, optionally within each Sample/capture.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input(
                    "metrics",
                    default="total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
                ),
                io.String.Input("batch_key", default="sample"),
                io.Float.Input("nmads", default=5.0, min=0.1, max=100.0, step=0.1),
                io.Combo.Input("direction", options=["both", "upper", "lower"], default="both"),
                io.String.Input("output_column", default="outlier", advanced=True),
                io.Boolean.Input("scale_mad", default=False, advanced=True),
                io.Int.Input("minimum_group_size", default=3, min=1, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellScrublet(io.ComfyNode):
    EXPRESSION_SOURCE = _SCRUBLET_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellScrublet",
            display_name="Run Scrublet",
            category=CATEGORY,
            description="Score and predict capture-level doublets from unnormalized counts without filtering cells.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("batch_key", default="sample"),
                io.Int.Input("random_seed", default=123, min=0, max=MAX_RANDOM_SEED, advanced=True),
                io.Float.Input("expected_doublet_rate", default=0.05, min=0.0001, max=0.9999, step=0.01),
                io.Combo.Input("threshold_mode", options=["automatic", "manual"], default="automatic"),
                io.Float.Input("threshold", default=0.25, step=0.01),
                io.Float.Input("sim_doublet_ratio", default=2.0, min=0.1, max=100.0, step=0.1, advanced=True),
                io.Int.Input("n_prin_comps", default=30, min=1, max=10_000, advanced=True),
                io.Int.Input("n_neighbors", default=0, min=0, max=1_000_000, advanced=True),
                cls.EXPRESSION_SOURCE.input(),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellFilterDoublets(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterDoublets",
            display_name="Filter Predicted Doublets",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("prediction_column", default="predicted_doublet", advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


CORRECTION_NODE_CLASSES = [
    OpenBioSingleCellMarkMADOutliers,
    OpenBioSingleCellScrublet,
    OpenBioSingleCellFilterDoublets,
]

__all__ = ["CORRECTION_NODE_CLASSES"]
