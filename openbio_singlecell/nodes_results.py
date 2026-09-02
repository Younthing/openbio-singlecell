from __future__ import annotations

from comfy_api.latest import io

from .expression_source import _EMBEDDING_PLOT_SPEC, _MARKER_PLOT_SPEC, _MARKER_SPEC
from .marker_evidence import MARKER_COLUMNS, MARKER_METHODS
from .node_types import AnnDataType, PlotResultType, TableResultType, analysis_outputs

MARKER_CATEGORY = "openbio/single-cell/marker-evidence"
PLOT_CATEGORY = "openbio/single-cell/visualization"
DIAGNOSTIC_CATEGORY = "openbio/single-cell/diagnostics"


class OpenBioSingleCellMarkerGenes(io.ComfyNode):
    EXPRESSION_SOURCE = _MARKER_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerGenes",
            display_name="Marker Genes",
            category=MARKER_CATEGORY,
            description=(
                "Compute exploratory all-groups-versus-rest Cluster marker evidence; this is not a Sample-level "
                "Condition contrast."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.Combo.Input("method", options=list(MARKER_METHODS), default="wilcoxon"),
                cls.EXPRESSION_SOURCE.input(),
                io.Int.Input("n_genes", default=100, min=0, max=2**31 - 1),
                io.Boolean.Input("tie_correct", default=True, advanced=True),
                io.Int.Input(
                    "max_output_rows",
                    default=1_000_000,
                    min=1,
                    max=100_000_000,
                    advanced=True,
                ),
                io.Float.Input(
                    "max_working_memory_gib",
                    default=4.0,
                    min=0.000001,
                    step=0.25,
                    advanced=True,
                ),
            ],
            outputs=analysis_outputs(
                TableResultType.Output(display_name="table"), TableResultType.Output(display_name="universe")
            ),
        )


class OpenBioSingleCellUMAPPlot(io.ComfyNode):
    EXPRESSION_SOURCE = _EMBEDDING_PLOT_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellUMAPPlot",
            display_name="Embedding Plot",
            category=PLOT_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("embedding_key", default="X_umap"),
                io.Int.Input("x_dimension", default=1, min=1, max=4096),
                io.Int.Input("y_dimension", default=2, min=1, max=4096),
                io.DynamicCombo.Input(
                    "color",
                    options=[
                        io.DynamicCombo.Option(
                            "obs",
                            [
                                io.String.Input("obs_key", default="leiden"),
                                io.Combo.Input(
                                    "color_mode",
                                    options=["auto", "categorical", "continuous"],
                                    default="auto",
                                    advanced=True,
                                ),
                            ],
                        ),
                        io.DynamicCombo.Option(
                            "gene",
                            [io.String.Input("gene", default=""), cls.EXPRESSION_SOURCE.input()],
                        ),
                        io.DynamicCombo.Option("none", []),
                    ],
                ),
                io.Float.Input("point_size", default=10.0, min=0.1, max=1000.0, step=1.0),
                io.String.Input("continuous_color_map", default="viridis", advanced=True),
                io.String.Input("categorical_palette", default="tab20", advanced=True),
                io.Boolean.Input("sort_order", default=True, advanced=True),
                io.String.Input("missing_color", default="lightgray", advanced=True),
                io.Combo.Input(
                    "legend_policy",
                    options=["automatic", "show", "hide"],
                    default="automatic",
                    advanced=True,
                ),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellFilterMarkerGenes(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterMarkerGenes",
            display_name="Filter Marker Genes",
            category=MARKER_CATEGORY,
            description=(
                "Filter a validated Cluster marker evidence table and its tested universe by inclusive thresholds "
                "without recomputing the marker tests."
            ),
            inputs=[
                TableResultType.Input("table"),
                TableResultType.Input("universe"),
                io.Float.Input("min_log2_fold_change", default=1.0, step=0.1),
                io.Float.Input("min_fraction_in_group", default=0.25, min=0.0, max=1.0, step=0.05),
                io.Float.Input("max_fraction_reference", default=0.5, min=0.0, max=1.0, step=0.05),
                io.Float.Input("max_p_adjusted", default=0.05, min=0.0, max=1.0, step=0.01),
            ],
            outputs=analysis_outputs(
                TableResultType.Output(display_name="table"), TableResultType.Output(display_name="universe")
            ),
        )


class OpenBioSingleCellMarkerExpressionPlot(io.ComfyNode):
    EXPRESSION_SOURCE = _MARKER_PLOT_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerExpressionPlot",
            display_name="Grouped Gene Expression Plot",
            category=PLOT_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("genes", default=""),
                io.String.Input("groupby", default="leiden"),
                io.DynamicCombo.Input(
                    "plot",
                    options=[
                        io.DynamicCombo.Option(
                            "dotplot",
                            [
                                io.Combo.Input(
                                    "standard_scale",
                                    options=["none", "var", "group"],
                                    default="var",
                                    advanced=True,
                                ),
                                io.Float.Input("expression_cutoff", default=0.0, step=0.1, advanced=True),
                                io.Boolean.Input("mean_only_expressed", default=False, advanced=True),
                            ],
                        ),
                        io.DynamicCombo.Option(
                            "matrixplot",
                            [
                                io.Combo.Input(
                                    "standard_scale",
                                    options=["none", "var", "group"],
                                    default="var",
                                    advanced=True,
                                )
                            ],
                        ),
                        io.DynamicCombo.Option("tracksplot", []),
                        io.DynamicCombo.Option(
                            "violin",
                            [
                                io.Combo.Input(
                                    "density_norm",
                                    options=["width", "area", "count"],
                                    default="width",
                                    advanced=True,
                                ),
                                io.Boolean.Input("show_cells", default=False, advanced=True),
                                io.Combo.Input(
                                    "y_scale",
                                    options=["linear", "log"],
                                    default="linear",
                                    advanced=True,
                                ),
                            ],
                        ),
                    ],
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.Combo.Input(
                    "group_order",
                    options=["observed", "dendrogram"],
                    default="observed",
                    advanced=True,
                ),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellPCAMetadataAssociations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPCAMetadataAssociations",
            display_name="PCA Metadata Associations",
            category=DIAGNOSTIC_CATEGORY,
            description="Screen PCA covariation with metadata using Samples as independent replicates.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("use_rep", default="X_pca"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("categorical_obs_keys", default=""),
                io.String.Input("continuous_obs_keys", default=""),
                io.Float.Input("alpha", default=0.05, min=0.0, max=1.0, step=0.01, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


RESULT_NODE_CLASSES = [
    OpenBioSingleCellMarkerGenes,
    OpenBioSingleCellUMAPPlot,
    OpenBioSingleCellFilterMarkerGenes,
    OpenBioSingleCellMarkerExpressionPlot,
    OpenBioSingleCellPCAMetadataAssociations,
]

__all__ = [
    "MARKER_COLUMNS",
    "RESULT_NODE_CLASSES",
    "OpenBioSingleCellFilterMarkerGenes",
    "OpenBioSingleCellMarkerExpressionPlot",
    "OpenBioSingleCellMarkerGenes",
    "OpenBioSingleCellPCAMetadataAssociations",
    "OpenBioSingleCellUMAPPlot",
]
