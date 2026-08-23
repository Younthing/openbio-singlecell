from __future__ import annotations

import importlib
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from .analysis_utils import finish_adata
from .node_types import AnnDataType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/copy-number"


def _require_infercnvpy() -> Any:
    try:
        return importlib.import_module("infercnvpy")
    except (ImportError, OSError) as error:
        raise RuntimeError("Copy-number analysis requires the infercnvpy package.") from error


def _comma_separated_values(value: str) -> list[str]:
    values = list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))
    if not values:
        raise ValueError("At least one reference category is required.")
    return values


class OpenBioSingleCellInferCNV(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellInferCNV",
            display_name="Infer CNV",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("reference_key", default="cell_type"),
                io.String.Input("reference_categories", default=""),
                io.Int.Input("window_size", default=250, min=1, max=2**31 - 1),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        reference_key: str = "cell_type",
        reference_categories: str = "",
        window_size: int = 250,
    ) -> io.NodeOutput:
        if reference_key not in adata.obs:
            raise ValueError(f"InferCNV reference column not found in obs: {reference_key!r}")
        coordinate_columns = ["chromosome", "start", "end"]
        missing_coordinates = [column for column in coordinate_columns if column not in adata.var]
        if missing_coordinates:
            raise ValueError(f"InferCNV gene coordinate columns not found in var: {missing_coordinates}")

        categories = _comma_separated_values(reference_categories)
        infercnvpy = _require_infercnvpy()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        infercnvpy.tl.infercnv(
            output,
            reference_key=reference_key,
            reference_cat=categories,
            window_size=window_size,
        )
        parameters = {
            "reference_key": reference_key,
            "reference_categories": categories,
            "window_size": window_size,
        }
        finish_adata(output, "infercnv", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellCNVStructure(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNVStructure",
            display_name="CNV Structure",
            category=CATEGORY,
            description="Compute the infercnvpy embedding, clusters, and per-cell CNV score.",
            inputs=[AnnDataType.Input("adata")],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData) -> io.NodeOutput:
        infercnvpy = _require_infercnvpy()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        infercnvpy.tl.pca(output)
        infercnvpy.pp.neighbors(output)
        infercnvpy.tl.leiden(output)
        infercnvpy.tl.umap(output)
        infercnvpy.tl.cnv_score(output)
        finish_adata(output, "cnv_structure", {}, cells, genes, started_at)
        return io.NodeOutput(output)


CNV_NODE_CLASSES = [
    OpenBioSingleCellInferCNV,
    OpenBioSingleCellCNVStructure,
]


__all__ = [
    "CNV_NODE_CLASSES",
    "OpenBioSingleCellCNVStructure",
    "OpenBioSingleCellInferCNV",
]
