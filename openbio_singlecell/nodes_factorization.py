from __future__ import annotations

from comfy_api.latest import io

from .cnmf_standalone import CNMF_LOCAL_NEIGHBORHOOD_SIZE
from .expression_source import ExpressionSourceSpec
from .node_types import AnnDataType, CNMFRunType, SummaryResultType, TableResultType

CATEGORY = "openbio/single-cell/factorization"


class OpenBioSingleCellCNMFRankSurvey(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="cNMF count source",
        default="layer",
        include_raw=True,
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNMFRankSurvey",
            display_name="cNMF Rank Survey",
            category=CATEGORY,
            description="Run the official cNMF rank family and publish a Worker-bound native artifact.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Int.Input("components_min", default=3, min=2),
                io.Int.Input("components_max", default=19, min=2),
                io.Int.Input("n_iter", default=100, min=2),
                io.Int.Input("num_highvar_genes", default=2000, min=1),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                CNMFRunType.Output(display_name="run"),
                TableResultType.Output(display_name="k_metrics"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellCNMF(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNMF",
            display_name="cNMF Consensus Programs",
            category=CATEGORY,
            description="Calculate consensus programs in a private writable checkout of a cNMF native artifact.",
            inputs=[
                CNMFRunType.Input("run"),
                io.Int.Input("selected_k", default=7, min=2),
                io.Float.Input("density_threshold", default=2.0, min=0.0, max=2.0, step=0.05),
                io.Float.Input(
                    "local_neighborhood_size",
                    default=CNMF_LOCAL_NEIGHBORHOOD_SIZE,
                    min=0.001,
                    max=1.0,
                    step=0.05,
                    advanced=True,
                ),
                io.Int.Input("n_top_genes", default=100, min=1, max=2**31 - 1, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


FACTORIZATION_NODE_CLASSES = [OpenBioSingleCellCNMFRankSurvey, OpenBioSingleCellCNMF]


__all__ = [
    "FACTORIZATION_NODE_CLASSES",
    "OpenBioSingleCellCNMF",
    "OpenBioSingleCellCNMFRankSurvey",
]
