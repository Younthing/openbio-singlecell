from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .cnv_analysis import (
    CNVState,
    analyze_cnv_pca,
    analyze_cnv_score,
    analyze_infer_cnv,
    cnv_pca_code,
    cnv_score_code,
    infer_cnv_code,
)
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, CNVStateType, SummaryResultType, TableResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/copy-number"


class OpenBioSingleCellInferCNV(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Full-gene normalized log-expression source for CNV inference",
        default="layer",
        include_raw=False,
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellInferCNV",
            display_name="Infer CNV",
            category=CATEGORY,
            description=(
                "Infer a fingerprinted expression-derived CNV state from an explicit normal reference, "
                "full-gene transformed expression, and named genome assembly."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("reference_key", default="cell_type"),
                io.String.Input("reference_categories", default=""),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("genome_assembly", default=""),
                io.Int.Input("window_size", default=100, min=2, max=10000),
                io.Int.Input("step", default=10, min=1, max=10000, advanced=True),
                io.Float.Input("lfc_clip", default=3.0, min=1e-12, advanced=True),
                io.Float.Input("dynamic_threshold", default=1.5, min=0.0, advanced=True),
                io.String.Input("exclude_chromosomes", default="chrX,chrY,chrM", advanced=True),
                io.String.Input("output_key", default="cnv", advanced=True),
                io.Int.Input("minimum_reference_cells", default=20, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("chunksize", default=5000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("n_jobs", default=1, min=1, max=256, advanced=True),
                io.Float.Input("max_output_gib", default=4.0, min=0.01, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                CNVStateType.Output(display_name="cnv_state"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        reference_key: str = "cell_type",
        reference_categories: str = "",
        sample_key: str = "sample",
        genome_assembly: str = "",
        window_size: int = 100,
        step: int = 10,
        lfc_clip: float = 3.0,
        dynamic_threshold: float = 1.5,
        exclude_chromosomes: str = "chrX,chrY,chrM",
        output_key: str = "cnv",
        minimum_reference_cells: int = 20,
        chunksize: int = 5000,
        n_jobs: int = 1,
        max_output_gib: float = 4.0,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        parameters = {
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "reference_key": reference_key,
            "reference_categories": reference_categories,
            "sample_key": sample_key,
            "genome_assembly": genome_assembly,
            "window_size": window_size,
            "step": step,
            "lfc_clip": lfc_clip,
            "dynamic_threshold": dynamic_threshold,
            "exclude_chromosomes": exclude_chromosomes,
            "output_key": output_key,
            "minimum_reference_cells": minimum_reference_cells,
            "chunksize": chunksize,
            "n_jobs": n_jobs,
            "max_output_gib": max_output_gib,
            "overwrite_existing": overwrite_existing,
        }
        cnv_state, summary = analyze_infer_cnv(adata, **parameters)
        report = make_summary_result(
            summary=summary,
            title="Expression-derived CNV inference summary",
            operation="infer_cnv",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
        )
        return io.NodeOutput(cnv_state, report, infer_cnv_code(**parameters))


class OpenBioSingleCellCNVPCA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNVPCA",
            display_name="CNV PCA",
            category=CATEGORY,
            description="Reduce one validated inferred-CNV window matrix; graph construction remains separate.",
            inputs=[
                CNVStateType.Input("cnv_state"),
                io.Int.Input("n_comps", default=30, min=1, max=4096),
                io.String.Input("output_key", default="X_cnv_pca", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_output_gib", default=2.0, min=0.01, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        cnv_state: CNVState,
        n_comps: int = 30,
        output_key: str = "X_cnv_pca",
        overwrite_existing: bool = False,
        max_output_gib: float = 2.0,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        input_adata = cnv_state.to_adata()
        cells, genes = int(input_adata.n_obs), int(input_adata.n_vars)
        parameters = {
            "n_comps": n_comps,
            "output_key": output_key,
            "overwrite_existing": overwrite_existing,
            "max_output_gib": max_output_gib,
            "random_seed": random_seed,
        }
        output, summary = analyze_cnv_pca(cnv_state, **parameters)
        finish_adata(
            output,
            "cnv_pca",
            dict(summary["parameters"]),
            cells,
            genes,
            started_at,
            random_seed=random_seed,
            warnings=[str(value) for value in summary["warnings"]],
        )
        report = make_summary_result(
            summary=summary,
            title="CNV PCA summary",
            operation="cnv_pca",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output, report, cnv_pca_code(**parameters))


class OpenBioSingleCellCNVScore(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNVScore",
            display_name="CNV Group Score",
            category=CATEGORY,
            description=(
                "Assign and independently verify descriptive mean absolute inferred-CNV scores for an explicit "
                "categorical partition; no tumor cutoff is inferred."
            ),
            inputs=[
                CNVStateType.Input("cnv_state"),
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cnv_leiden"),
                io.String.Input("output_key", default="cnv_score", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        cnv_state: CNVState,
        adata: AnnData,
        groupby: str = "cnv_leiden",
        output_key: str = "cnv_score",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        parameters = {
            "groupby": groupby,
            "output_key": output_key,
            "overwrite_existing": overwrite_existing,
        }
        output, table, summary = analyze_cnv_score(cnv_state, adata, **parameters)
        finish_adata(
            output,
            "cnv_score",
            dict(summary["parameters"]),
            cells,
            genes,
            started_at,
            warnings=[str(value) for value in summary["warnings"]],
        )
        table_result = make_table_result(
            title=f"CNV group scores by {groupby}",
            operation="cnv_score",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title="CNV group score summary",
            operation="cnv_score",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
        )
        return io.NodeOutput(output, table_result, report, cnv_score_code(**parameters))


CNV_NODE_CLASSES = [
    OpenBioSingleCellInferCNV,
    OpenBioSingleCellCNVPCA,
    OpenBioSingleCellCNVScore,
]


__all__ = [
    "CNV_NODE_CLASSES",
    "OpenBioSingleCellCNVPCA",
    "OpenBioSingleCellCNVScore",
    "OpenBioSingleCellInferCNV",
]
