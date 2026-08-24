from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/preprocessing"


class OpenBioSingleCellNormalizeTotal(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellNormalizeTotal",
            display_name="Normalize Total",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("target_sum", default=10000.0, min=0.000001, step=1000.0),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, target_sum: float = 10000.0) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        if "counts" not in output.layers:
            output.layers["counts"] = output.X.copy()
        science.sc.pp.normalize_total(output, target_sum=target_sum, inplace=True)
        parameters = {"target_sum": target_sum}
        finish_adata(output, "normalize_total", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellLog1p(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLog1p",
            display_name="Log1p",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Boolean.Input("set_raw", default=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, set_raw: bool = True) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.pp.log1p(output)
        if set_raw:
            output.raw = output.copy()
        parameters = {"set_raw": set_raw}
        finish_adata(output, "log1p", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellNormalizeToLayer(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Normalize source",
        layer_input_id="source_layer",
        layer_default="counts",
    )

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
                io.Float.Input("target_sum", default=10000.0, min=0.0, step=1000.0),
                io.Combo.Input("transform", options=["none", "log1p", "sqrt"], default="log1p"),
                io.String.Input("output_layer", default="log1p_norm", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        target_sum: float = 10000.0,
        transform: str = "log1p",
        output_layer: str = "log1p_norm",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if not output_layer.strip():
            raise ValueError("Normalize output_layer cannot be empty.")
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        matrix = expression.matrix(output)
        work = science.ad.AnnData(X=matrix.copy(), obs=output.obs.copy(), var=output.var.copy())
        science.sc.pp.normalize_total(work, target_sum=target_sum or None, inplace=True)
        if transform == "log1p":
            science.sc.pp.log1p(work)
        elif transform == "sqrt":
            work.X = work.X.sqrt() if science.sparse.issparse(work.X) else science.np.sqrt(work.X)
        output.layers[output_layer.strip()] = work.X.copy()
        parameters = {
            **expression.parameters(),
            "target_sum": target_sum,
            "transform": transform,
            "output_layer": output_layer.strip(),
        }
        finish_adata(output, "normalize_to_layer", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellPearsonResidualsToLayer(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Pearson residual source",
        layer_input_id="source_layer",
        layer_default="counts",
    )

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
                io.Float.Input("theta", default=100.0, min=0.000001, step=10.0),
                io.Float.Input("clip", default=0.0, min=0.0, step=1.0, advanced=True),
                io.String.Input("output_layer", default="analytic_pearson_residuals", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        theta: float = 100.0,
        clip: float = 0.0,
        output_layer: str = "analytic_pearson_residuals",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        output_layer = output_layer.strip()
        if not output_layer:
            raise ValueError("Pearson residual output_layer cannot be empty.")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        normalized = science.sc.experimental.pp.normalize_pearson_residuals(
            output,
            theta=theta,
            clip=clip or None,
            layer=expression.scanpy_layer,
            inplace=False,
        )
        output.layers[output_layer] = normalized["X"]
        parameters = {
            **expression.parameters(),
            "theta": theta,
            "clip": clip,
            "output_layer": output_layer,
        }
        finish_adata(output, "pearson_residuals_to_layer", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellHighlyVariableGenes(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Highly Variable Genes source",
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHighlyVariableGenes",
            display_name="Highly Variable Genes",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_top_genes", default=2000, min=1, max=2**31 - 1),
                io.Combo.Input(
                    "flavor",
                    options=["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper", "pearson_residuals"],
                    default="seurat",
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("batch_key", default=""),
                io.String.Input("always_keep_genes", default=""),
                io.Boolean.Input("subset", default=False),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        n_top_genes: int = 2000,
        flavor: str = "seurat",
        source: DynamicExpressionSource | None = None,
        batch_key: str = "",
        always_keep_genes: str = "",
        subset: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        expression = cls.EXPRESSION_SOURCE.resolve(output, source)
        if batch_key and batch_key not in output.obs:
            raise ValueError(f"Highly Variable Genes batch column not found in obs: {batch_key!r}")
        requested_keep = list(dict.fromkeys(gene.strip() for gene in always_keep_genes.split(",") if gene.strip()))
        hvg = (
            science.sc.experimental.pp.highly_variable_genes
            if flavor == "pearson_residuals"
            else science.sc.pp.highly_variable_genes
        )
        hvg(
            output,
            layer=expression.scanpy_layer,
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key or None,
            subset=False,
            inplace=True,
        )
        present_keep = [gene for gene in requested_keep if gene in output.var_names]
        if present_keep:
            output.var.loc[present_keep, "highly_variable"] = True
        if subset:
            output = output[:, output.var["highly_variable"].astype(bool)].copy()
        missing_keep = [gene for gene in requested_keep if gene not in output.var_names]
        parameters = {
            "n_top_genes": n_top_genes,
            "flavor": flavor,
            **expression.parameters(),
            "batch_key": batch_key,
            "always_keep_genes": requested_keep,
            "subset": subset,
        }
        warnings = [f"Requested keep genes were not found: {missing_keep}"] if missing_keep else []
        finish_adata(output, "highly_variable_genes", parameters, cells, genes, started_at, warnings=warnings)
        return io.NodeOutput(output)


class OpenBioSingleCellScale(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellScale",
            display_name="Scale",
            category=CATEGORY,
            description="Scale genes using Scanpy; this can densify sparse matrices.",
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("max_value", default=10.0, min=0.0, step=1.0),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, max_value: float = 10.0) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.pp.scale(output, max_value=max_value if max_value > 0 else None)
        parameters = {"max_value": max_value}
        finish_adata(output, "scale", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


PREPROCESS_NODE_CLASSES = [
    OpenBioSingleCellNormalizeTotal,
    OpenBioSingleCellLog1p,
    OpenBioSingleCellNormalizeToLayer,
    OpenBioSingleCellPearsonResidualsToLayer,
    OpenBioSingleCellHighlyVariableGenes,
    OpenBioSingleCellScale,
]
