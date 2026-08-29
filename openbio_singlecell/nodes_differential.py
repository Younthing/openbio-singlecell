from __future__ import annotations

from comfy_api.latest import io

from .expression_source import ExpressionSourceSpec
from .node_types import AnnDataType, PseudobulkType, SCVIModelType, SummaryResultType, TableResultType

CATEGORY = "openbio/single-cell/differential-expression"


class OpenBioSingleCellPseudobulk(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Pseudobulk raw-count source",
        default="layer",
        include_raw=True,
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulk",
            display_name="Pseudobulk",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("population_key", default="cell_type"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("technical_batch_key", default="", advanced=True),
                io.String.Input("categorical_covariate_keys", default="", advanced=True),
                io.String.Input("continuous_covariate_keys", default="", advanced=True),
                cls.EXPRESSION_SOURCE.input(),
                io.Combo.Input("inference_mode", options=["formal", "exploratory"], default="formal"),
                io.Int.Input("min_cells", default=10, min=1, max=2**31 - 1),
                io.Int.Input("min_counts", default=1000, min=0, max=2**31 - 1),
            ],
            outputs=[
                PseudobulkType.Output(display_name="pseudobulk"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellPseudobulkEdgeR(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulkEdgeR",
            display_name="Pseudobulk edgeR",
            category=CATEGORY,
            inputs=[
                PseudobulkType.Input("pseudobulk"),
                io.String.Input("population", default=""),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("categorical_covariate_keys", default="", advanced=True),
                io.String.Input("continuous_covariate_keys", default="", advanced=True),
                io.Float.Input("fdr_threshold", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Float.Input("min_abs_log2_fold_change", default=0.0, min=0.0, step=0.1),
                io.Int.Input("min_count", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("min_total_count", default=15, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("large_n", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Float.Input("min_prop", default=0.7, min=0.0, max=1.0, step=0.05, advanced=True),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellPseudobulkDESeq2(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulkDESeq2",
            display_name="Pseudobulk PyDESeq2 Contrast",
            category=CATEGORY,
            inputs=[
                PseudobulkType.Input("pseudobulk"),
                io.String.Input("population", default=""),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("categorical_covariate_keys", default="", advanced=True),
                io.String.Input("continuous_covariate_keys", default="", advanced=True),
                io.Float.Input("fdr_threshold", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Float.Input("min_abs_log2_fold_change", default=0.0, min=0.0, step=0.1),
                io.Int.Input("min_count", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("min_total_count", default=15, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("large_n", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Float.Input("min_prop", default=0.7, min=0.0, max=1.0, step=0.05, advanced=True),
                io.Int.Input("n_cpus", default=1, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellSCVIDifferentialExpression(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCVIDifferentialExpression",
            display_name="scVI Model DE Evidence",
            category=CATEGORY,
            description=(
                "Compare two cell populations with one supplied fitted scVI model. The result is exploratory "
                "model evidence, not Sample-level Condition inference."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                SCVIModelType.Input("model"),
                io.String.Input("groupby", default="group"),
                io.String.Input("group1", default=""),
                io.String.Input("group2", default=""),
                io.DynamicCombo.Input(
                    "population_scope",
                    options=[
                        io.DynamicCombo.Option("all", []),
                        io.DynamicCombo.Option(
                            "obs_value",
                            [
                                io.String.Input("subset_column", default=""),
                                io.String.Input("subset_value", default=""),
                            ],
                        ),
                    ],
                ),
                io.DynamicCombo.Input(
                    "mode",
                    options=[
                        io.DynamicCombo.Option(
                            "change",
                            [
                                io.Float.Input("delta", default=0.25, min=1e-9, step=0.05),
                                io.Float.Input(
                                    "fdr_target",
                                    default=0.05,
                                    min=1e-9,
                                    max=0.999999999,
                                    step=0.01,
                                    advanced=True,
                                ),
                            ],
                        ),
                        io.DynamicCombo.Option("vanilla", []),
                    ],
                ),
                io.Combo.Input(
                    "batch_handling",
                    options=["shared_technical_batches", "observed_technical_batches"],
                    default="shared_technical_batches",
                ),
                io.Int.Input("n_samples_overall", default=5000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


DIFFERENTIAL_NODE_CLASSES = [
    OpenBioSingleCellPseudobulk,
    OpenBioSingleCellPseudobulkEdgeR,
    OpenBioSingleCellPseudobulkDESeq2,
    OpenBioSingleCellSCVIDifferentialExpression,
]


__all__ = [
    "DIFFERENTIAL_NODE_CLASSES",
    "OpenBioSingleCellPseudobulk",
    "OpenBioSingleCellPseudobulkEdgeR",
    "OpenBioSingleCellPseudobulkDESeq2",
    "OpenBioSingleCellSCVIDifferentialExpression",
]
