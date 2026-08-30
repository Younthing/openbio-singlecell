from __future__ import annotations

from comfy_api.latest import io

from .expression_source import (
    _HVG_SPEC,
    _NORMALIZE_LAYER_SPEC,
    _NORMALIZE_TOTAL_SPEC,
    _PEARSON_RESIDUAL_SPEC,
    _SCALE_SPEC,
)
from .node_types import AnnDataType, analysis_outputs

CATEGORY = "openbio/single-cell/preprocessing"


class OpenBioSingleCellNormalizeTotal(io.ComfyNode):
    EXPRESSION_SOURCE = _NORMALIZE_TOTAL_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellNormalizeTotal",
            display_name="Normalize Total",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("target_sum", default=10000.0, step=1000.0),
                cls.EXPRESSION_SOURCE.input(),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellLog1p(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLog1p",
            display_name="Log1p",
            category=CATEGORY,
            inputs=[AnnDataType.Input("adata")],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellNormalizeToLayer(io.ComfyNode):
    EXPRESSION_SOURCE = _NORMALIZE_LAYER_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellNormalizeToLayer",
            display_name="Normalize to Layer",
            category=CATEGORY,
            description="Normalize an expression matrix into a named layer without replacing adata.X.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Float.Input("target_sum", default=10000.0, step=1000.0),
                io.Combo.Input("transform", options=["none", "log1p", "sqrt"], default="log1p"),
                io.String.Input("output_layer", default="log1p_norm", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellPearsonResidualsToLayer(io.ComfyNode):
    EXPRESSION_SOURCE = _PEARSON_RESIDUAL_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPearsonResidualsToLayer",
            display_name="Pearson Residuals to Layer",
            category=CATEGORY,
            description="Compute analytic Pearson residuals and store them in an AnnData layer.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Float.Input("theta", default=100.0, step=10.0, advanced=True),
                io.Combo.Input(
                    "clipping_mode",
                    options=["sqrt_n_obs", "custom", "none"],
                    default="sqrt_n_obs",
                ),
                io.Float.Input("custom_clip", default=10.0, min=0.0, step=1.0, advanced=True),
                io.String.Input("output_layer", default="analytic_pearson_residuals", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_dense_gib", default=2.0, step=0.25, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellHighlyVariableGenes(io.ComfyNode):
    EXPRESSION_SOURCE = _HVG_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHighlyVariableGenes",
            display_name="Highly Variable Genes",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_top_genes", default=2000, min=1),
                io.Combo.Input(
                    "flavor",
                    options=["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper", "pearson_residuals"],
                    default="seurat",
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("batch_key", default=""),
                io.String.Input("always_keep_genes", default="", advanced=True),
                io.Boolean.Input("subset", default=False, advanced=True),
                io.Float.Input("theta", default=100.0, step=10.0, advanced=True),
                io.Combo.Input(
                    "clipping_mode",
                    options=["sqrt_n_obs", "custom", "none"],
                    default="sqrt_n_obs",
                    advanced=True,
                ),
                io.Float.Input("custom_clip", default=10.0, min=0.0, step=1.0, advanced=True),
                io.Int.Input("chunksize", default=1000, min=1, advanced=True),
                io.Float.Input("span", default=0.3, max=1.0, step=0.05, advanced=True),
                io.Int.Input("n_bins", default=20, min=1, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


class OpenBioSingleCellScale(io.ComfyNode):
    EXPRESSION_SOURCE = _SCALE_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellScale",
            display_name="Scale",
            category=CATEGORY,
            description="Scale one expression source into a named layer without replacing X or Raw.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Boolean.Input("zero_center", default=True),
                io.Combo.Input("clipping_mode", options=["custom", "none"], default="custom"),
                io.Float.Input("custom_max_value", default=10.0, min=0.0, step=1.0, advanced=True),
                io.String.Input("output_layer", default="scaled", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_dense_gib", default=2.0, step=0.25, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


PREPROCESS_NODE_CLASSES = [
    OpenBioSingleCellNormalizeTotal,
    OpenBioSingleCellLog1p,
    OpenBioSingleCellNormalizeToLayer,
    OpenBioSingleCellPearsonResidualsToLayer,
    OpenBioSingleCellHighlyVariableGenes,
    OpenBioSingleCellScale,
]
