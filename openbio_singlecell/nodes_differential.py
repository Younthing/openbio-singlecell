from __future__ import annotations

import time
from collections.abc import Mapping
from typing import TYPE_CHECKING

from comfy_api.latest import io

from . import PLUGIN_VERSION, dependencies
from .analysis_utils import make_summary_result, make_table_result
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, PseudobulkType, SCVIModelType, SummaryResultType, TableResultType
from .pseudobulk import build_pseudobulk_summary, pseudobulk_code, run_pseudobulk
from .sample_design import (
    build_pseudobulk_engine_summary,
    pseudobulk_engine_code,
    run_pseudobulk_edger,
    run_pseudobulk_pydeseq2,
)
from .scvi_de_evidence import (
    analyze_scvi_model_de_evidence,
    resolve_scvi_de_mode,
    resolve_scvi_population_scope,
    scvi_de_code,
)
from .scvi_model import SCVIModel

if TYPE_CHECKING:
    from anndata import AnnData

    from .pseudobulk import PseudobulkArtifact


CATEGORY = "openbio/single-cell/differential-expression"


class OpenBioSingleCellPseudobulk(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Pseudobulk raw-count source",
        default="layer",
        include_raw=True,
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulk",
            display_name="Pseudobulk",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("population_key", default="cell_type"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("technical_batch_key", default="", advanced=True),
                io.String.Input("categorical_covariate_keys", default="", advanced=True),
                io.String.Input("continuous_covariate_keys", default="", advanced=True),
                cls.EXPRESSION_SOURCE.input(),
                io.Combo.Input("inference_mode", options=["formal", "exploratory"], default="formal"),
                io.Int.Input("min_cells", default=10, min=1, max=2**31 - 1),
                io.Int.Input("min_counts", default=1000, min=0, max=2**31 - 1),
            ],
            outputs=[
                PseudobulkType.Output(display_name="pseudobulk"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        population_key: str = "cell_type",
        condition_key: str = "condition",
        technical_batch_key: str = "",
        categorical_covariate_keys: str = "",
        continuous_covariate_keys: str = "",
        source: DynamicExpressionSource | None = None,
        inference_mode: str = "formal",
        min_cells: int = 10,
        min_counts: int = 1000,
    ) -> io.NodeOutput:
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells = int(adata.n_obs)
        artifact, diagnostics = run_pseudobulk(
            adata,
            sample_key=sample_key,
            population_key=population_key,
            condition_key=condition_key,
            technical_batch_key=technical_batch_key,
            categorical_covariate_keys=categorical_covariate_keys,
            continuous_covariate_keys=continuous_covariate_keys,
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            inference_mode=inference_mode,
            min_cells=min_cells,
            min_counts=min_counts,
            openbio_version=PLUGIN_VERSION,
            decoupler_module=None,
        )
        code_parameters = {
            "sample_key": sample_key,
            "population_key": population_key,
            "condition_key": condition_key,
            "technical_batch_key": technical_batch_key,
            "categorical_covariate_keys": categorical_covariate_keys,
            "continuous_covariate_keys": continuous_covariate_keys,
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "inference_mode": inference_mode,
            "min_cells": min_cells,
            "min_counts": min_counts,
            "openbio_version": PLUGIN_VERSION,
        }
        code = pseudobulk_code(**code_parameters)
        summary = build_pseudobulk_summary(diagnostics)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        report = make_summary_result(
            summary=summary,
            title="Sample-by-population pseudobulk counts",
            operation="pseudobulk",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=cells,
            input_genes=int(diagnostics["input_genes"]),
            started_at=started_at,
        )
        return io.NodeOutput(artifact, report, code)


class OpenBioSingleCellPseudobulkEdgeR(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulkEdgeR",
            display_name="Pseudobulk edgeR",
            category=CATEGORY,
            inputs=[
                PseudobulkType.Input("pseudobulk"),
                io.String.Input("population", default=""),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("categorical_covariate_keys", default="", advanced=True),
                io.String.Input("continuous_covariate_keys", default="", advanced=True),
                io.Float.Input("fdr_threshold", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Float.Input("min_abs_log2_fold_change", default=0.0, min=0.0, step=0.1),
                io.Int.Input("min_count", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("min_total_count", default=15, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("large_n", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Float.Input("min_prop", default=0.7, min=0.0, max=1.0, step=0.05, advanced=True),
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
        pseudobulk: PseudobulkArtifact,
        population: str = "",
        reference_condition: str = "",
        comparison_condition: str = "",
        categorical_covariate_keys: str = "",
        continuous_covariate_keys: str = "",
        fdr_threshold: float = 0.05,
        min_abs_log2_fold_change: float = 0.0,
        min_count: int = 10,
        min_total_count: int = 15,
        large_n: int = 10,
        min_prop: float = 0.7,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        runtime_parameters = {
            "population": population,
            "reference_condition": reference_condition,
            "comparison_condition": comparison_condition,
            "categorical_covariate_keys": categorical_covariate_keys,
            "continuous_covariate_keys": continuous_covariate_keys,
            "fdr_threshold": fdr_threshold,
            "min_abs_log2_fold_change": min_abs_log2_fold_change,
            "min_count": min_count,
            "min_total_count": min_total_count,
            "large_n": large_n,
            "min_prop": min_prop,
            "openbio_version": PLUGIN_VERSION,
        }
        table, diagnostics = run_pseudobulk_edger(
            pseudobulk,
            **runtime_parameters,
            pertpy_module=None,
        )
        summary = build_pseudobulk_engine_summary(diagnostics)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        metadata = pseudobulk.metadata
        result = make_table_result(
            title=f"edgeR QL: {comparison_condition.strip()} vs {reference_condition.strip()}",
            operation="pseudobulk_edger",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=int(metadata["input_dimensions"]["cells"]),
            input_genes=int(diagnostics["expression_filter"]["genes_before"]),
            started_at=started_at,
            table=science.pd.DataFrame(table),
        )
        code = pseudobulk_engine_code(engine="edger", **runtime_parameters)
        report = make_summary_result(
            summary=summary,
            title="Pseudobulk edgeR QL Condition contrast",
            operation="pseudobulk_edger",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=int(metadata["input_dimensions"]["cells"]),
            input_genes=int(diagnostics["expression_filter"]["genes_before"]),
            started_at=started_at,
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellPseudobulkDESeq2(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulkDESeq2",
            display_name="Pseudobulk PyDESeq2 Contrast",
            category=CATEGORY,
            inputs=[
                PseudobulkType.Input("pseudobulk"),
                io.String.Input("population", default=""),
                io.String.Input("reference_condition", default=""),
                io.String.Input("comparison_condition", default=""),
                io.String.Input("categorical_covariate_keys", default="", advanced=True),
                io.String.Input("continuous_covariate_keys", default="", advanced=True),
                io.Float.Input("fdr_threshold", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Float.Input("min_abs_log2_fold_change", default=0.0, min=0.0, step=0.1),
                io.Int.Input("min_count", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("min_total_count", default=15, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("large_n", default=10, min=0, max=2**31 - 1, advanced=True),
                io.Float.Input("min_prop", default=0.7, min=0.0, max=1.0, step=0.05, advanced=True),
                io.Int.Input("n_cpus", default=1, min=1, max=2**31 - 1, advanced=True),
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
        pseudobulk: PseudobulkArtifact,
        population: str = "",
        reference_condition: str = "",
        comparison_condition: str = "",
        categorical_covariate_keys: str = "",
        continuous_covariate_keys: str = "",
        fdr_threshold: float = 0.05,
        min_abs_log2_fold_change: float = 0.0,
        min_count: int = 10,
        min_total_count: int = 15,
        large_n: int = 10,
        min_prop: float = 0.7,
        n_cpus: int = 1,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        runtime_parameters = {
            "population": population,
            "reference_condition": reference_condition,
            "comparison_condition": comparison_condition,
            "categorical_covariate_keys": categorical_covariate_keys,
            "continuous_covariate_keys": continuous_covariate_keys,
            "fdr_threshold": fdr_threshold,
            "min_abs_log2_fold_change": min_abs_log2_fold_change,
            "min_count": min_count,
            "min_total_count": min_total_count,
            "large_n": large_n,
            "min_prop": min_prop,
            "n_cpus": n_cpus,
            "openbio_version": PLUGIN_VERSION,
        }
        table, diagnostics = run_pseudobulk_pydeseq2(
            pseudobulk,
            **runtime_parameters,
            pertpy_module=None,
        )
        summary = build_pseudobulk_engine_summary(diagnostics)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        metadata = pseudobulk.metadata
        result = make_table_result(
            title=f"PyDESeq2: {comparison_condition.strip()} vs {reference_condition.strip()}",
            operation="pseudobulk_deseq2",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=int(metadata["input_dimensions"]["cells"]),
            input_genes=int(diagnostics["expression_filter"]["genes_before"]),
            started_at=started_at,
            table=science.pd.DataFrame(table),
        )
        code = pseudobulk_engine_code(engine="pydeseq2", **runtime_parameters)
        report = make_summary_result(
            summary=summary,
            title="Pseudobulk PyDESeq2 Condition contrast",
            operation="pseudobulk_deseq2",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=int(metadata["input_dimensions"]["cells"]),
            input_genes=int(diagnostics["expression_filter"]["genes_before"]),
            started_at=started_at,
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellSCVIDifferentialExpression(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCVIDifferentialExpression",
            display_name="scVI Model DE Evidence",
            category=CATEGORY,
            description=(
                "Compare two cell populations with one supplied fitted scVI model. The result is exploratory "
                "model evidence, not Sample-level Condition inference."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                SCVIModelType.Input("model"),
                io.String.Input("groupby", default="group"),
                io.String.Input("group1", default=""),
                io.String.Input("group2", default=""),
                io.DynamicCombo.Input(
                    "population_scope",
                    options=[
                        io.DynamicCombo.Option("all", []),
                        io.DynamicCombo.Option(
                            "obs_value",
                            [
                                io.String.Input("subset_column", default=""),
                                io.String.Input("subset_value", default=""),
                            ],
                        ),
                    ],
                ),
                io.DynamicCombo.Input(
                    "mode",
                    options=[
                        io.DynamicCombo.Option(
                            "change",
                            [
                                io.Float.Input("delta", default=0.25, min=1e-9, step=0.05),
                                io.Float.Input(
                                    "fdr_target",
                                    default=0.05,
                                    min=1e-9,
                                    max=0.999999999,
                                    step=0.01,
                                    advanced=True,
                                ),
                            ],
                        ),
                        io.DynamicCombo.Option("vanilla", []),
                    ],
                ),
                io.Combo.Input(
                    "batch_handling",
                    options=["shared_technical_batches", "observed_technical_batches"],
                    default="shared_technical_batches",
                ),
                io.Int.Input("n_samples_overall", default=5000, min=1, max=2**31 - 1, advanced=True),
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
        model: SCVIModel,
        groupby: str = "group",
        group1: str = "",
        group2: str = "",
        population_scope: Mapping[str, object] | None = None,
        mode: Mapping[str, object] | None = None,
        batch_handling: str = "shared_technical_batches",
        n_samples_overall: int = 5000,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        subset_column, subset_value = resolve_scvi_population_scope(population_scope)
        mode_name, delta, fdr_target = resolve_scvi_de_mode(mode)
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        table, summary = analyze_scvi_model_de_evidence(
            adata,
            model,
            groupby=groupby,
            group1=group1,
            group2=group2,
            subset_column=subset_column,
            subset_value=subset_value,
            mode=mode_name,
            delta=delta,
            fdr_target=fdr_target,
            batch_handling=batch_handling,
            n_samples_overall=n_samples_overall,
            random_seed=random_seed,
        )
        parameters = dict(summary["parameters"])
        description = summary["results"]
        warnings = [str(warning) for warning in summary["warnings"]]
        input_cells = int(summary["comparison"]["selected_scope_cells"])
        result = make_table_result(
            title=f"scVI model evidence: {group1} vs {group2}",
            operation="scvi_model_de_evidence",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=len(model.var_names),
            started_at=started_at,
            table=science.pd.DataFrame(table),
        )
        code = scvi_de_code(
            groupby=groupby,
            group1=group1,
            group2=group2,
            subset_column=subset_column,
            subset_value=subset_value,
            mode=mode_name,
            delta=delta,
            fdr_target=fdr_target,
            batch_handling=batch_handling,
            n_samples_overall=n_samples_overall,
            random_seed=random_seed,
        )
        report = make_summary_result(
            summary=summary,
            title="scVI model DE evidence summary",
            operation="scvi_model_de_evidence",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=len(model.var_names),
            started_at=started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(result, report, code)


DIFFERENTIAL_NODE_CLASSES = [
    OpenBioSingleCellPseudobulk,
    OpenBioSingleCellPseudobulkEdgeR,
    OpenBioSingleCellPseudobulkDESeq2,
    OpenBioSingleCellSCVIDifferentialExpression,
]


__all__ = [
    "DIFFERENTIAL_NODE_CLASSES",
    "OpenBioSingleCellPseudobulk",
    "OpenBioSingleCellPseudobulkEdgeR",
    "OpenBioSingleCellPseudobulkDESeq2",
    "OpenBioSingleCellSCVIDifferentialExpression",
]
