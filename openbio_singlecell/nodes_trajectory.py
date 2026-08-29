from __future__ import annotations

from comfy_api.latest import io

from .expression_source import ExpressionSourceSpec
from .node_types import AnnDataType, SummaryResultType

CATEGORY = "openbio/single-cell/trajectory"
ANNOTATION_CATEGORY = "openbio/single-cell/annotation"
CELL_CYCLE_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="Cell-cycle expression source; full-gene log-normalized values are recommended",
    default="layer",
    include_raw=True,
    layer_default="log1p_norm",
)


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
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
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
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
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
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
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
                io.Int.Input("n_dcs", default=10, min=2, max=4096),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

TRAJECTORY_NODE_CLASSES = [
    OpenBioSingleCellCellCycleScore,
    OpenBioSingleCellDiffusionMap,
    OpenBioSingleCellPAGA,
    OpenBioSingleCellDPT,
]
