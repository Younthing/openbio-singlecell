from __future__ import annotations

from comfy_api.latest import io

from .expression_source import ExpressionSourceSpec
from .node_types import AnnDataType, PlotResultType, SummaryResultType

CATEGORY = "openbio/single-cell/qc"
MAX_THRESHOLD = 2**31 - 1
CALCULATE_QC_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="QC count expression source",
    default="X",
    include_raw=True,
    layer_input_id="source_layer",
    layer_default="counts",
)
FILTER_CELLS_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="Cell-filter count expression source",
    default="X",
    include_raw=True,
    layer_input_id="source_layer",
    layer_default="counts",
)
FILTER_GENES_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="Gene-filter count expression source",
    default="X",
    include_raw=True,
    layer_input_id="source_layer",
    layer_default="counts",
)
QC_PLOTS_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="QC plot expression source",
    default="X",
    include_raw=True,
    layer_input_id="source_layer",
    layer_default="counts",
)


class OpenBioSingleCellCalculateQC(io.ComfyNode):
    EXPRESSION_SOURCE = CALCULATE_QC_EXPRESSION_SOURCE

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCalculateQC",
            display_name="Calculate QC Metrics",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Boolean.Input("include_ribosomal", default=True),
                io.Boolean.Input("include_hemoglobin", default=True),
                io.String.Input("mitochondrial_prefix", default="MT-", advanced=True),
                io.String.Input("ribosomal_prefixes", default="RPS,RPL", advanced=True),
                io.String.Input("hemoglobin_pattern", default=r"^HB[^(P)]", advanced=True),
                io.String.Input("percent_top", default="20,50,100,200,500", advanced=True),
                io.Boolean.Input("log1p", default=True, advanced=True),
                cls.EXPRESSION_SOURCE.input(),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellFilterCells(io.ComfyNode):
    EXPRESSION_SOURCE = FILTER_CELLS_EXPRESSION_SOURCE

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterCells",
            display_name="Filter Cells",
            category=CATEGORY,
            description=(
                "Filter cells; zero disables detected-gene bounds, while expression and mitochondrial bounds "
                "use explicit enable controls."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("min_genes", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("max_genes", default=0, min=0, max=MAX_THRESHOLD),
                io.Float.Input("min_counts", default=0.0, step=1.0),
                io.Float.Input("max_counts", default=0.0, step=1.0),
                io.Float.Input("max_pct_mito", default=0.0, min=0.0, max=100.0, step=1.0),
                io.String.Input("mito_column", default="pct_counts_mt", advanced=True),
                cls.EXPRESSION_SOURCE.input(),
                io.Boolean.Input("enable_min_counts", default=False),
                io.Boolean.Input("enable_max_counts", default=False),
                io.Boolean.Input("enable_max_pct_mito", default=False),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellFilterGenes(io.ComfyNode):
    EXPRESSION_SOURCE = FILTER_GENES_EXPRESSION_SOURCE

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterGenes",
            display_name="Filter Genes",
            category=CATEGORY,
            description=(
                "Filter genes; zero disables detected-cell bounds, while total-expression bounds use explicit "
                "enable controls."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("min_cells", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("max_cells", default=0, min=0, max=MAX_THRESHOLD),
                io.Float.Input("min_counts", default=0.0, step=1.0),
                io.Float.Input("max_counts", default=0.0, step=1.0),
                cls.EXPRESSION_SOURCE.input(),
                io.Boolean.Input("enable_min_counts", default=False),
                io.Boolean.Input("enable_max_counts", default=False),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellQCPlots(io.ComfyNode):
    EXPRESSION_SOURCE = QC_PLOTS_EXPRESSION_SOURCE

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellQCPlots",
            display_name="QC Plots",
            category=CATEGORY,
            description="Visualize a coherent family of cell-level QC metrics derived from one explicit source.",
            inputs=[AnnDataType.Input("adata"), cls.EXPRESSION_SOURCE.input()],
            outputs=[
                PlotResultType.Output(display_name="plot"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


QC_NODE_CLASSES = [
    OpenBioSingleCellCalculateQC,
    OpenBioSingleCellFilterCells,
    OpenBioSingleCellFilterGenes,
    OpenBioSingleCellQCPlots,
]
