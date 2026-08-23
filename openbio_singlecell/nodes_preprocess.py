from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata
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
                io.Boolean.Input("set_raw", default=True, advanced=True),
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


class OpenBioSingleCellHighlyVariableGenes(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHighlyVariableGenes",
            display_name="Highly Variable Genes",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_top_genes", default=2000, min=1, max=2**31 - 1),
                io.Combo.Input("flavor", options=["seurat", "cell_ranger"], default="seurat"),
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
        subset: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        science.sc.pp.highly_variable_genes(
            output,
            n_top_genes=n_top_genes,
            flavor=flavor,
            subset=subset,
            inplace=True,
        )
        parameters = {"n_top_genes": n_top_genes, "flavor": flavor, "subset": subset}
        finish_adata(output, "highly_variable_genes", parameters, cells, genes, started_at)
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
    OpenBioSingleCellHighlyVariableGenes,
    OpenBioSingleCellScale,
]
