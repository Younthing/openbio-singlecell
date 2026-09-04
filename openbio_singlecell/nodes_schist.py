from __future__ import annotations

from comfy_api.latest import io

from .node_types import AnnDataType, PlotResultType, analysis_outputs

SCHIST_CATEGORY = "openbio/single-cell/clustering"


class OpenBioSingleCellSchistNestedModel(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSchistNestedModel",
            display_name="Schist Nested-SBM Hierarchy",
            category=SCHIST_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("random_seed", default=123, min=1, max=2**31 - 1, advanced=True),
                io.String.Input("neighbors_key", default="neighbors", advanced=True),
                io.String.Input("key_added", default="nsbm", advanced=True),
                io.Int.Input("posterior_samples", default=100, min=100, max=2**31 - 1, advanced=True),
                io.Boolean.Input("degree_correction", default=True, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_working_memory_gib", default=8.0, min=0.001, step=0.5, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellSchistHierarchyPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSchistHierarchyPlot",
            display_name="Schist Hierarchy Plot",
            category=SCHIST_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("key_added", default="nsbm"),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


SCHIST_NODE_CLASSES = [OpenBioSingleCellSchistNestedModel, OpenBioSingleCellSchistHierarchyPlot]


__all__ = [
    "OpenBioSingleCellSchistNestedModel",
    "OpenBioSingleCellSchistHierarchyPlot",
    "SCHIST_NODE_CLASSES",
]
