from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_result
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/differential-expression"


def _require_pertpy() -> Any:
    try:
        import pertpy
    except (ImportError, OSError) as error:
        raise RuntimeError("Pseudobulk analysis requires the pertpy package.") from error
    return pertpy


class OpenBioSingleCellPseudobulk(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulk",
            display_name="Pseudobulk",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("groupby", default="cell_type"),
                io.String.Input("counts_layer", default="counts"),
                io.Combo.Input("mode", options=["sum", "mean"], default="sum"),
                io.Int.Input("min_cells", default=10, min=1, max=2**31 - 1),
                io.Int.Input("min_counts", default=1000, min=0, max=2**31 - 1),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        groupby: str = "cell_type",
        counts_layer: str = "counts",
        mode: str = "sum",
        min_cells: int = 10,
        min_counts: int = 1000,
    ) -> io.NodeOutput:
        if sample_key not in adata.obs:
            raise ValueError(f"Pseudobulk sample column not found in obs: {sample_key!r}")
        if groupby not in adata.obs:
            raise ValueError(f"Pseudobulk group column not found in obs: {groupby!r}")
        if counts_layer not in adata.layers:
            raise ValueError(f"Pseudobulk counts layer not found: {counts_layer!r}")

        pertpy = _require_pertpy()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = pertpy.tl.PseudobulkSpace().compute(
            adata,
            target_col=sample_key,
            groups_col=groupby,
            layer_key=counts_layer,
            mode=mode,
            min_cells=min_cells,
            min_counts=min_counts,
        )
        output.layers["counts"] = output.X.copy()
        parameters = {
            "sample_key": sample_key,
            "groupby": groupby,
            "counts_layer": counts_layer,
            "mode": mode,
            "min_cells": min_cells,
            "min_counts": min_counts,
        }
        finish_adata(output, "pseudobulk", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellPseudobulkEdgeR(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulkEdgeR",
            display_name="Pseudobulk edgeR",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("design", default="~group"),
                io.String.Input("contrast_column", default="group"),
                io.String.Input("baseline", default="control"),
                io.String.Input("comparison", default="treatment"),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        design: str = "~group",
        contrast_column: str = "group",
        baseline: str = "control",
        comparison: str = "treatment",
    ) -> io.NodeOutput:
        if contrast_column not in adata.obs:
            raise ValueError(f"edgeR contrast column not found in obs: {contrast_column!r}")

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        model = pertpy.tl.EdgeR(adata, design=design)
        model.fit()
        contrast = model.contrast(
            column=contrast_column,
            baseline=baseline,
            group_to_compare=comparison,
        )
        table = model.test_contrasts(contrast).reset_index()
        if "index" in table.columns and "gene" not in table.columns:
            table = table.rename(columns={"index": "gene"})

        parameters = {
            "design": design,
            "contrast_column": contrast_column,
            "baseline": baseline,
            "comparison": comparison,
        }
        result = make_result(
            kind="table",
            title=f"edgeR: {comparison} vs {baseline}",
            operation="pseudobulk_edger",
            parameters=parameters,
            description="Pseudobulk differential expression from pertpy's edgeR wrapper.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=science.pd.DataFrame(table),
        )
        return io.NodeOutput(result)


DIFFERENTIAL_NODE_CLASSES = [
    OpenBioSingleCellPseudobulk,
    OpenBioSingleCellPseudobulkEdgeR,
]


__all__ = [
    "DIFFERENTIAL_NODE_CLASSES",
    "OpenBioSingleCellPseudobulk",
    "OpenBioSingleCellPseudobulkEdgeR",
]
