from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import PLUGIN_VERSION
from .analysis_utils import make_summary_result, make_table_result
from .composition_modeling import (
    run_sccoda_differential_composition,
    run_tasccoda_differential_composition,
    sccoda_differential_composition_code,
    tasccoda_differential_composition_code,
)
from .composition_summary import run_sample_composition_summary, sample_composition_summary_code
from .milo_analysis import milo_differential_abundance_code, run_milo_differential_abundance
from .node_types import AnnDataType, SummaryResultType, TableResultType

if TYPE_CHECKING:
    from anndata import AnnData


ABUNDANCE_CATEGORY = "openbio/single-cell/differential-abundance"


class OpenBioSingleCellSampleCompositionSummary(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSampleCompositionSummary",
            display_name="Sample Composition Summary",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                ),
                io.Int.Input(
                    "max_output_rows",
                    default=2_000_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        condition_key: str = "condition",
        annotation_key: str = "cell_type",
        annotation_status: str = "unknown",
        max_output_rows: int = 2_000_000,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, summary = run_sample_composition_summary(
            adata,
            sample_key=sample_key,
            condition_key=condition_key,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            max_output_rows=max_output_rows,
            openbio_version=PLUGIN_VERSION,
        )
        code = sample_composition_summary_code(
            sample_key=sample_key,
            condition_key=condition_key,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            max_output_rows=max_output_rows,
            openbio_version=PLUGIN_VERSION,
        )
        result = make_table_result(
            title=f"Sample composition by {annotation_key}",
            operation="sample_composition_summary",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title="Sample composition summary",
            operation="sample_composition_summary",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellMiloDifferentialAbundance(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMiloDifferentialAbundance",
            display_name="Milo Differential Abundance",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("technical_batch_key", default="", advanced=True),
                io.String.Input("categorical_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("continuous_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
                io.String.Input("representation_key", default="X_pca"),
                io.Int.Input("n_neighbors", default=30, min=2, max=2**31 - 1),
                io.Float.Input(
                    "neighborhood_proportion",
                    default=0.1,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Float.Input(
                    "mixed_annotation_threshold",
                    default=0.6,
                    min=0.0,
                    max=1.0,
                    step=0.05,
                    advanced=True,
                ),
                io.Float.Input(
                    "spatial_fdr_threshold",
                    default=0.1,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Float.Input(
                    "min_abs_log2_fold_change",
                    default=0.0,
                    min=0.0,
                    step=0.1,
                    advanced=True,
                ),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        condition_key: str = "condition",
        reference_condition: str = "",
        comparison_condition: str = "",
        technical_batch_key: str = "",
        categorical_covariate_keys_json: str = "[]",
        continuous_covariate_keys_json: str = "[]",
        annotation_key: str = "cell_type",
        annotation_status: str = "provisional",
        representation_key: str = "X_pca",
        n_neighbors: int = 30,
        neighborhood_proportion: float = 0.1,
        mixed_annotation_threshold: float = 0.6,
        spatial_fdr_threshold: float = 0.1,
        min_abs_log2_fold_change: float = 0.0,
        random_seed: int = 123,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        table, summary = run_milo_differential_abundance(
            adata,
            sample_key=sample_key,
            condition_key=condition_key,
            reference_condition=reference_condition,
            comparison_condition=comparison_condition,
            technical_batch_key=technical_batch_key,
            categorical_covariate_keys_json=categorical_covariate_keys_json,
            continuous_covariate_keys_json=continuous_covariate_keys_json,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            representation_key=representation_key,
            n_neighbors=n_neighbors,
            neighborhood_proportion=neighborhood_proportion,
            mixed_annotation_threshold=mixed_annotation_threshold,
            spatial_fdr_threshold=spatial_fdr_threshold,
            min_abs_log2_fold_change=min_abs_log2_fold_change,
            random_seed=random_seed,
            openbio_version=PLUGIN_VERSION,
        )
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        table_result = make_table_result(
            title="Milo differential abundance",
            operation="milo_differential_abundance",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title="Milo differential abundance summary",
            operation="milo_differential_abundance",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        code = milo_differential_abundance_code(
            sample_key=sample_key,
            condition_key=condition_key,
            reference_condition=reference_condition,
            comparison_condition=comparison_condition,
            technical_batch_key=technical_batch_key,
            categorical_covariate_keys_json=categorical_covariate_keys_json,
            continuous_covariate_keys_json=continuous_covariate_keys_json,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            representation_key=representation_key,
            n_neighbors=n_neighbors,
            neighborhood_proportion=neighborhood_proportion,
            mixed_annotation_threshold=mixed_annotation_threshold,
            spatial_fdr_threshold=spatial_fdr_threshold,
            min_abs_log2_fold_change=min_abs_log2_fold_change,
            random_seed=random_seed,
            openbio_version=PLUGIN_VERSION,
        )
        return io.NodeOutput(table_result, report, code)


class OpenBioSingleCellSccodaDifferentialComposition(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSccodaDifferentialComposition",
            display_name="scCODA Differential Composition",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("adjustment_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("reference_cell_type", default="automatic"),
                io.Float.Input("estimated_fdr", default=0.05, min=0.001, max=0.999, step=0.01),
                io.Int.Input("num_samples", default=10000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("num_warmup", default=1000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        annotation_key: str = "cell_type",
        annotation_status: str = "provisional",
        condition_key: str = "condition",
        reference_condition: str = "",
        comparison_condition: str = "",
        adjustment_covariate_keys_json: str = "[]",
        reference_cell_type: str = "automatic",
        estimated_fdr: float = 0.05,
        num_samples: int = 10000,
        num_warmup: int = 1000,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        table, summary = run_sccoda_differential_composition(
            adata,
            sample_key=sample_key,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            condition_key=condition_key,
            reference_condition=reference_condition,
            comparison_condition=comparison_condition,
            adjustment_covariate_keys_json=adjustment_covariate_keys_json,
            reference_cell_type=reference_cell_type,
            estimated_fdr=estimated_fdr,
            num_samples=num_samples,
            num_warmup=num_warmup,
            random_seed=random_seed,
            openbio_version=PLUGIN_VERSION,
        )
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        description = summary["results"]["writing_summary"]
        table_result = make_table_result(
            title="scCODA differential composition",
            operation="sccoda_differential_composition",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title="scCODA differential composition summary",
            operation="sccoda_differential_composition",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        code = sccoda_differential_composition_code(
            sample_key=sample_key,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            condition_key=condition_key,
            reference_condition=reference_condition,
            comparison_condition=comparison_condition,
            adjustment_covariate_keys_json=adjustment_covariate_keys_json,
            reference_cell_type=reference_cell_type,
            estimated_fdr=estimated_fdr,
            num_samples=num_samples,
            num_warmup=num_warmup,
            random_seed=random_seed,
            openbio_version=PLUGIN_VERSION,
        )
        return io.NodeOutput(table_result, report, code)


class OpenBioSingleCellTasccodaDifferentialComposition(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellTasccodaDifferentialComposition",
            display_name="tascCODA Differential Composition",
            category=ABUNDANCE_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
                io.String.Input("hierarchy_keys_json", default="[]"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("adjustment_covariate_keys_json", default="[]", advanced=True),
                io.String.Input("reference_cell_type", default="automatic"),
                io.Float.Input("aggregation_bias", default=0.0, step=0.1, advanced=True),
                io.Int.Input("num_samples", default=10000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("num_warmup", default=1000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        annotation_key: str = "cell_type",
        annotation_status: str = "provisional",
        hierarchy_keys_json: str = "[]",
        condition_key: str = "condition",
        reference_condition: str = "",
        comparison_condition: str = "",
        adjustment_covariate_keys_json: str = "[]",
        reference_cell_type: str = "automatic",
        aggregation_bias: float = 0.0,
        num_samples: int = 10000,
        num_warmup: int = 1000,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        table, summary = run_tasccoda_differential_composition(
            adata,
            sample_key=sample_key,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            hierarchy_keys_json=hierarchy_keys_json,
            condition_key=condition_key,
            reference_condition=reference_condition,
            comparison_condition=comparison_condition,
            adjustment_covariate_keys_json=adjustment_covariate_keys_json,
            reference_cell_type=reference_cell_type,
            aggregation_bias=aggregation_bias,
            num_samples=num_samples,
            num_warmup=num_warmup,
            random_seed=random_seed,
            openbio_version=PLUGIN_VERSION,
        )
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        description = summary["results"]["writing_summary"]
        table_result = make_table_result(
            title="tascCODA differential composition",
            operation="tasccoda_differential_composition",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title="tascCODA differential composition summary",
            operation="tasccoda_differential_composition",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        code = tasccoda_differential_composition_code(
            sample_key=sample_key,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            hierarchy_keys_json=hierarchy_keys_json,
            condition_key=condition_key,
            reference_condition=reference_condition,
            comparison_condition=comparison_condition,
            adjustment_covariate_keys_json=adjustment_covariate_keys_json,
            reference_cell_type=reference_cell_type,
            aggregation_bias=aggregation_bias,
            num_samples=num_samples,
            num_warmup=num_warmup,
            random_seed=random_seed,
            openbio_version=PLUGIN_VERSION,
        )
        return io.NodeOutput(table_result, report, code)


ABUNDANCE_NODE_CLASSES = [
    OpenBioSingleCellSampleCompositionSummary,
    OpenBioSingleCellMiloDifferentialAbundance,
    OpenBioSingleCellSccodaDifferentialComposition,
    OpenBioSingleCellTasccodaDifferentialComposition,
]


__all__ = [
    "ABUNDANCE_NODE_CLASSES",
    "OpenBioSingleCellMiloDifferentialAbundance",
    "OpenBioSingleCellSampleCompositionSummary",
    "OpenBioSingleCellSccodaDifferentialComposition",
    "OpenBioSingleCellTasccodaDifferentialComposition",
]
