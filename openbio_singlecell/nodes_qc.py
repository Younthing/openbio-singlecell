from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import figure_to_png, finish_adata, make_result, matrix_totals_and_nonzero
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/qc"
MAX_THRESHOLD = 2**31 - 1


class OpenBioSingleCellCalculateQC(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCalculateQC",
            display_name="Calculate QC Metrics",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("mitochondrial_prefix", default="MT-", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, mitochondrial_prefix: str = "MT-") -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        output.var["mt"] = output.var_names.astype(str).str.startswith(mitochondrial_prefix)
        warnings = []
        if not bool(output.var["mt"].any()):
            warnings.append(f"No genes matched mitochondrial prefix {mitochondrial_prefix!r}.")
        dependencies.sc.pp.calculate_qc_metrics(
            output,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )
        parameters = {"mitochondrial_prefix": mitochondrial_prefix}
        finish_adata(output, "calculate_qc", parameters, cells, genes, started_at, warnings=warnings)
        return io.NodeOutput(output)


class OpenBioSingleCellFilterCells(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterCells",
            display_name="Filter Cells",
            category=CATEGORY,
            description="Filter cells; a threshold of zero disables that bound.",
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("min_genes", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("max_genes", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("min_counts", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("max_counts", default=0, min=0, max=MAX_THRESHOLD),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        min_genes: int = 0,
        max_genes: int = 0,
        min_counts: int = 0,
        max_counts: int = 0,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        counts, detected = matrix_totals_and_nonzero(adata.X, axis=1)
        mask = dependencies.np.ones(cells, dtype=bool)
        if min_genes > 0:
            mask &= detected >= min_genes
        if max_genes > 0:
            mask &= detected <= max_genes
        if min_counts > 0:
            mask &= counts >= min_counts
        if max_counts > 0:
            mask &= counts <= max_counts
        output = adata[mask].copy()
        parameters = {
            "min_genes": min_genes,
            "max_genes": max_genes,
            "min_counts": min_counts,
            "max_counts": max_counts,
        }
        finish_adata(output, "filter_cells", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellFilterGenes(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterGenes",
            display_name="Filter Genes",
            category=CATEGORY,
            description="Filter genes; a threshold of zero disables that bound.",
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("min_cells", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("max_cells", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("min_counts", default=0, min=0, max=MAX_THRESHOLD),
                io.Int.Input("max_counts", default=0, min=0, max=MAX_THRESHOLD),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        min_cells: int = 0,
        max_cells: int = 0,
        min_counts: int = 0,
        max_counts: int = 0,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        counts, detected = matrix_totals_and_nonzero(adata.X, axis=0)
        mask = dependencies.np.ones(genes, dtype=bool)
        if min_cells > 0:
            mask &= detected >= min_cells
        if max_cells > 0:
            mask &= detected <= max_cells
        if min_counts > 0:
            mask &= counts >= min_counts
        if max_counts > 0:
            mask &= counts <= max_counts
        output = adata[:, mask].copy()
        parameters = {
            "min_cells": min_cells,
            "max_cells": max_cells,
            "min_counts": min_counts,
            "max_counts": max_counts,
        }
        finish_adata(output, "filter_genes", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellQCPlots(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellQCPlots",
            display_name="QC Plots",
            category=CATEGORY,
            inputs=[AnnDataType.Input("adata")],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(cls, adata: AnnData) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        if cells == 0 or genes == 0:
            raise ValueError("QC plots require at least one cell and one gene.")

        total_counts, detected = matrix_totals_and_nonzero(adata.X, axis=1)
        warnings = []
        if "pct_counts_mt" in adata.obs:
            pct_mt = dependencies.np.asarray(adata.obs["pct_counts_mt"], dtype=float)
        elif "mt" in adata.var:
            mt_mask = dependencies.np.asarray(adata.var["mt"], dtype=bool)
            mt_counts, _ = matrix_totals_and_nonzero(adata.X[:, mt_mask], axis=1)
            pct_mt = dependencies.np.divide(
                mt_counts * 100.0,
                total_counts,
                out=dependencies.np.zeros_like(total_counts, dtype=float),
                where=total_counts != 0,
            )
        else:
            pct_mt = dependencies.np.zeros(cells, dtype=float)
            warnings.append(
                "Mitochondrial annotations are unavailable; percent mitochondrial counts are shown as zero."
            )

        figure = dependencies.Figure(figsize=(10, 7), constrained_layout=True)
        axes = figure.subplots(2, 2)
        axes[0, 0].hist(total_counts, bins=40, color="#246bfe")
        axes[0, 0].set_title("Total counts per cell")
        axes[0, 1].hist(detected, bins=40, color="#17a673")
        axes[0, 1].set_title("Genes detected per cell")
        axes[1, 0].scatter(total_counts, detected, s=8, alpha=0.6, color="#6f4bf2")
        axes[1, 0].set_xlabel("Total counts")
        axes[1, 0].set_ylabel("Genes detected")
        axes[1, 1].scatter(total_counts, pct_mt, s=8, alpha=0.6, color="#e45d3a")
        axes[1, 1].set_xlabel("Total counts")
        axes[1, 1].set_ylabel("Mitochondrial counts (%)")
        png = figure_to_png(figure)
        result = make_result(
            kind="plot",
            title="Quality control plots",
            operation="qc_plots",
            parameters={},
            description="Cell count, feature, and mitochondrial quality metrics.",
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            png=png,
        )
        return io.NodeOutput(result)


QC_NODE_CLASSES = [
    OpenBioSingleCellCalculateQC,
    OpenBioSingleCellFilterCells,
    OpenBioSingleCellFilterGenes,
    OpenBioSingleCellQCPlots,
]
