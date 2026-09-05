from __future__ import annotations

from comfy_api.latest import io

from .expression_source import _CELL_CYCLE_SPEC, _MARKER_PLOT_SPEC
from .node_types import AnnDataType, PlotResultType, analysis_outputs

CATEGORY = "openbio/single-cell/trajectory"
ANNOTATION_CATEGORY = "openbio/single-cell/annotation"
CELL_CYCLE_EXPRESSION_SOURCE = _CELL_CYCLE_SPEC


class OpenBioSingleCellCellCycleScore(io.ComfyNode):
    EXPRESSION_SOURCE = CELL_CYCLE_EXPRESSION_SOURCE

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellCycleScore",
            display_name="Cell Cycle Score",
            category=ANNOTATION_CATEGORY,
            description=(
                "Score one versioned human or custom cell-cycle program on an explicit Raw, X, or layer source; "
                "the report discloses departures from recommended full-gene log-normalized input."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.DynamicCombo.Input(
                    "gene_set_source",
                    options=[
                        io.DynamicCombo.Option("regev_human_97", []),
                        io.DynamicCombo.Option(
                            "custom",
                            [
                                io.String.Input("s_genes", default=""),
                                io.String.Input("g2m_genes", default=""),
                            ],
                        ),
                    ],
                ),
                io.Combo.Input("organism", options=["human", "mouse", "other"], default="human"),
                io.String.Input("output_prefix", default="cell_cycle", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellCellCycleScorePlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellCycleScorePlot",
            display_name="Cell Cycle Score Plot",
            category=ANNOTATION_CATEGORY,
            description=(
                "Render the stored S/G2M score plane and phase counts from Cell Cycle Score without rescoring genes."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("output_prefix", default="cell_cycle"),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )

class OpenBioSingleCellDiffusionMap(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDiffusionMap",
            display_name="Diffusion Map",
            category=CATEGORY,
            description=(
                "Compute and fingerprint one diffusion basis from an exact named neighbor graph; no root or "
                "pseudotime is selected."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.Int.Input("n_comps", default=15, min=3, max=4096),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )

class OpenBioSingleCellPAGA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPAGA",
            display_name="PAGA",
            category=CATEGORY,
            description=("Coarse-grain one named undirected neighbor graph over one explicit categorical partition."),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )

class OpenBioSingleCellDPT(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDPT",
            display_name="Diffusion Pseudotime",
            category=CATEGORY,
            description=(
                "Orient an already graph-bound Diffusion Map from one exact cell or deterministic nominated-"
                "population medoid."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.DynamicCombo.Input(
                    "root_mode",
                    options=[
                        io.DynamicCombo.Option(
                            "cell_id",
                            [io.String.Input("root_cell_id", default="")],
                        ),
                        io.DynamicCombo.Option(
                            "group_medoid",
                            [
                                io.String.Input("root_column", default="leiden"),
                                io.Combo.Input(
                                    "root_value_type",
                                    options=["string", "integer", "number"],
                                    default="string",
                                ),
                                io.String.Input("root_value", default=""),
                            ],
                        ),
                    ],
                ),
                io.Int.Input("n_dcs", default=10, min=1, max=4096),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellPAGAPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPAGAPlot",
            display_name="PAGA Plot",
            category=CATEGORY,
            description="Read-only circular layout of the verified stored PAGA abstraction.",
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("min_connectivity", default=0.0, min=0.0, max=1.0, step=0.05),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellDiffusionSpectrumPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDiffusionSpectrumPlot",
            display_name="Diffusion Spectrum Plot",
            category=CATEGORY,
            description="Read-only component spectrum from the graph-bound Diffusion Map bundle.",
            inputs=[AnnDataType.Input("adata")],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellDPTGeneTrendPlot(io.ComfyNode):
    EXPRESSION_SOURCE = _MARKER_PLOT_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDPTGeneTrendPlot",
            display_name="DPT Gene Trend Plot",
            category=CATEGORY,
            description=(
                "Read-only binned expression summaries along verified root-dependent diffusion pseudotime."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("genes", default=""),
                cls.EXPRESSION_SOURCE.input(),
                io.Int.Input("n_bins", default=20, min=1, max=100),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


TRAJECTORY_NODE_CLASSES = [
    OpenBioSingleCellCellCycleScore,
    OpenBioSingleCellCellCycleScorePlot,
    OpenBioSingleCellDiffusionMap,
    OpenBioSingleCellPAGA,
    OpenBioSingleCellDPT,
    OpenBioSingleCellPAGAPlot,
    OpenBioSingleCellDiffusionSpectrumPlot,
    OpenBioSingleCellDPTGeneTrendPlot,
]
