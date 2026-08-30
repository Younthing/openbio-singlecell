from __future__ import annotations

from comfy_api.latest import io

from .node_types import AnnDataType, TableResultType, analysis_outputs

ABUNDANCE_CATEGORY = "openbio/single-cell/differential-abundance"


class OpenBioSingleCellSampleCompositionSummary(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSampleCompositionSummary",
            display_name="Sample Composition Summary",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                ),
                io.Int.Input(
                    "max_output_rows",
                    default=2_000_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellMiloDifferentialAbundance(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMiloDifferentialAbundance",
            display_name="Milo Differential Abundance",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("technical_batch_key", default="", advanced=True),
                io.String.Input("categorical_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("continuous_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
                io.String.Input("representation_key", default="X_pca"),
                io.Int.Input("n_neighbors", default=30, min=2, max=2**31 - 1),
                io.Float.Input(
                    "neighborhood_proportion",
                    default=0.1,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Float.Input(
                    "mixed_annotation_threshold",
                    default=0.6,
                    min=0.0,
                    max=1.0,
                    step=0.05,
                    advanced=True,
                ),
                io.Float.Input(
                    "spatial_fdr_threshold",
                    default=0.1,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Float.Input(
                    "min_abs_log2_fold_change",
                    default=0.0,
                    min=0.0,
                    step=0.1,
                    advanced=True,
                ),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellSccodaDifferentialComposition(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSccodaDifferentialComposition",
            display_name="scCODA Differential Composition",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("adjustment_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("reference_cell_type", default="automatic"),
                io.Float.Input("estimated_fdr", default=0.05, min=0.001, max=0.999, step=0.01),
                io.Int.Input("num_samples", default=10000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("num_warmup", default=1000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellTasccodaDifferentialComposition(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellTasccodaDifferentialComposition",
            display_name="tascCODA Differential Composition",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
                io.String.Input("hierarchy_keys_json", default="[]"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("adjustment_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("reference_cell_type", default="automatic"),
                io.Float.Input("aggregation_bias", default=0.0, step=0.1, advanced=True),
                io.Int.Input("num_samples", default=10000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("num_warmup", default=1000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


ABUNDANCE_NODE_CLASSES = [
    OpenBioSingleCellSampleCompositionSummary,
    OpenBioSingleCellMiloDifferentialAbundance,
    OpenBioSingleCellSccodaDifferentialComposition,
    OpenBioSingleCellTasccodaDifferentialComposition,
]


__all__ = [
    "ABUNDANCE_NODE_CLASSES",
    "OpenBioSingleCellMiloDifferentialAbundance",
    "OpenBioSingleCellSampleCompositionSummary",
    "OpenBioSingleCellSccodaDifferentialComposition",
    "OpenBioSingleCellTasccodaDifferentialComposition",
]
