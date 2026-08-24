from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_table_result
from .node_types import AnnDataType, SCVIModelType, TableResultType
from .scvi_model import SCVIModel

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/differential-expression"


def _require_pertpy() -> Any:
    try:
        import pertpy
    except (ImportError, OSError) as error:
        raise RuntimeError("Pseudobulk analysis requires the pertpy package.") from error
    return pertpy


def _require_decoupler() -> Any:
    try:
        import decoupler
    except (ImportError, OSError) as error:
        raise RuntimeError("Decoupler pseudobulk contrast requires the decoupler package.") from error
    return decoupler


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
                io.Combo.Input("source", options=["X", "layer"], default="layer"),
                io.String.Input("layer_name", default="counts"),
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
        source: str = "layer",
        layer_name: str = "counts",
        mode: str = "sum",
        min_cells: int = 10,
        min_counts: int = 1000,
    ) -> io.NodeOutput:
        if sample_key not in adata.obs:
            raise ValueError(f"Pseudobulk sample column not found in obs: {sample_key!r}")
        if groupby not in adata.obs:
            raise ValueError(f"Pseudobulk group column not found in obs: {groupby!r}")
        if source not in {"X", "layer"}:
            raise ValueError(f"Unsupported pseudobulk expression source: {source!r}")
        layer_name = layer_name.strip()
        if source == "layer" and layer_name not in adata.layers:
            raise ValueError(f"Pseudobulk expression layer not found: {layer_name!r}")

        pertpy = _require_pertpy()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = pertpy.tl.PseudobulkSpace().compute(
            adata,
            target_col=sample_key,
            groups_col=groupby,
            layer_key=layer_name if source == "layer" else None,
            mode=mode,
            min_cells=min_cells,
            min_counts=min_counts,
        )
        if mode == "sum":
            output.layers["counts"] = output.X.copy()
        parameters = {
            "sample_key": sample_key,
            "groupby": groupby,
            "source": source,
            "layer_name": layer_name,
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
                io.String.Input("baseline", default=""),
                io.String.Input("comparison", default=""),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        design: str = "~group",
        contrast_column: str = "group",
        baseline: str = "",
        comparison: str = "",
    ) -> io.NodeOutput:
        if contrast_column not in adata.obs:
            raise ValueError(f"edgeR contrast column not found in obs: {contrast_column!r}")
        if not baseline or not comparison:
            raise ValueError("edgeR baseline and comparison values are required.")

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
        result = make_table_result(
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


class OpenBioSingleCellDecouplerPseudobulkContrast(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDecouplerPseudobulkContrast",
            display_name="Decoupler Pseudobulk Contrast",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("groupby", default="cell_type"),
                io.String.Input("condition_key", default="group"),
                io.String.Input("condition", default=""),
                io.String.Input("reference", default=""),
                io.Combo.Input("source", options=["X", "layer"], default="layer"),
                io.String.Input("layer_name", default="counts"),
                io.Combo.Input("method", options=["t-test", "wilcoxon"], default="t-test"),
                io.Float.Input("min_prop", default=0.1, min=0.0, max=1.0, step=0.05),
                io.Int.Input("min_samples", default=3, min=1, max=2**31 - 1),
                io.Float.Input("target_sum", default=10000.0, min=0.000001, step=1000.0, advanced=True),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        groupby: str = "cell_type",
        condition_key: str = "group",
        condition: str = "",
        reference: str = "",
        source: str = "layer",
        layer_name: str = "counts",
        method: str = "t-test",
        min_prop: float = 0.1,
        min_samples: int = 3,
        target_sum: float = 10000.0,
    ) -> io.NodeOutput:
        required_obs = [sample_key, groupby, condition_key]
        missing_obs = [key for key in required_obs if key not in adata.obs]
        if missing_obs:
            raise ValueError(f"Decoupler pseudobulk observation columns not found: {missing_obs}")
        if source == "layer" and layer_name not in adata.layers:
            raise ValueError(f"Decoupler pseudobulk layer not found: {layer_name!r}")
        if not condition or not reference:
            raise ValueError("Decoupler pseudobulk condition and reference values are required.")

        decoupler = _require_decoupler()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        pdata = decoupler.get_pseudobulk(
            adata,
            sample_col=sample_key,
            groups_col=groupby,
            min_prop=min_prop,
            min_smpls=min_samples,
            layer=layer_name if source == "layer" else None,
        )
        science.sc.pp.normalize_total(pdata, target_sum=target_sum)
        science.sc.pp.log1p(pdata)
        log_fold_changes, p_values = decoupler.get_contrast(
            pdata,
            group_col=groupby,
            condition_col=condition_key,
            condition=condition,
            reference=reference,
            method=method,
        )
        table = science.pd.DataFrame(decoupler.format_contrast_results(log_fold_changes, p_values)).reset_index()

        parameters = {
            "sample_key": sample_key,
            "groupby": groupby,
            "condition_key": condition_key,
            "condition": condition,
            "reference": reference,
            "source": source,
            "layer_name": layer_name,
            "method": method,
            "min_prop": min_prop,
            "min_samples": min_samples,
            "target_sum": target_sum,
        }
        result = make_table_result(
            title=f"Pseudobulk contrast: {condition} vs {reference}",
            operation="decoupler_pseudobulk_contrast",
            parameters=parameters,
            description="Decoupler pseudobulk contrast formatted as one differential-expression table.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellPseudobulkDESeq2(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPseudobulkDESeq2",
            display_name="Pseudobulk DESeq2",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("design", default="~group"),
                io.String.Input("contrast_column", default="group"),
                io.String.Input("baseline", default=""),
                io.String.Input("comparison", default=""),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        design: str = "~group",
        contrast_column: str = "group",
        baseline: str = "",
        comparison: str = "",
    ) -> io.NodeOutput:
        if contrast_column not in adata.obs:
            raise ValueError(f"DESeq2 contrast column not found in obs: {contrast_column!r}")
        if not baseline or not comparison:
            raise ValueError("DESeq2 baseline and comparison values are required.")

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        model = pertpy.tl.PyDESeq2(adata=adata, design=design)
        model.fit()
        table = model.test_contrasts([contrast_column, baseline, comparison]).reset_index()
        if "index" in table.columns and "gene" not in table.columns:
            table = table.rename(columns={"index": "gene"})

        parameters = {
            "design": design,
            "contrast_column": contrast_column,
            "baseline": baseline,
            "comparison": comparison,
        }
        result = make_table_result(
            title=f"DESeq2: {comparison} vs {baseline}",
            operation="pseudobulk_deseq2",
            parameters=parameters,
            description="Pseudobulk differential expression from pertpy's PyDESeq2 wrapper.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=science.pd.DataFrame(table),
        )
        return io.NodeOutput(result)


class OpenBioSingleCellSCVIDifferentialExpression(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCVIDifferentialExpression",
            display_name="scVI Differential Expression",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                SCVIModelType.Input("model"),
                io.String.Input("groupby", default="group"),
                io.String.Input("group1", default=""),
                io.String.Input("group2", default=""),
                io.String.Input("subset_column", default=""),
                io.String.Input("subset_value", default=""),
                io.Combo.Input("mode", options=["vanilla", "change"], default="vanilla"),
                io.Float.Input(
                    "delta",
                    default=0.25,
                    min=0.0,
                    step=0.05,
                    advanced=True,
                    tooltip="Only used when mode is 'change'.",
                ),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        model: SCVIModel,
        groupby: str = "group",
        group1: str = "",
        group2: str = "",
        subset_column: str = "",
        subset_value: str = "",
        mode: str = "vanilla",
        delta: float = 0.25,
    ) -> io.NodeOutput:
        if not isinstance(model, SCVIModel):
            raise TypeError("scVI Differential Expression requires an SCVIModel input.")
        groupby = groupby.strip()
        group1 = group1.strip()
        group2 = group2.strip()
        subset_column = subset_column.strip()
        subset_value = subset_value.strip()
        mode = mode.strip()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        table = model.differential_expression(
            adata,
            groupby=groupby,
            group1=group1,
            group2=group2,
            subset_column=subset_column,
            subset_value=subset_value,
            mode=mode,
            delta=delta,
        ).reset_index()
        input_cells = int(adata.n_obs)
        if subset_column:
            input_cells = int((adata.obs[subset_column].astype(str) == subset_value).sum())
        if "index" in table.columns and "gene" not in table.columns:
            table = table.rename(columns={"index": "gene"})

        parameters = {
            "groupby": groupby,
            "group1": group1,
            "group2": group2,
            "subset_column": subset_column,
            "subset_value": subset_value,
            "mode": mode,
            "delta": delta,
            "model_class": model.model_class,
            "model_training_parameters": dict(model.training_parameters),
        }
        result = make_table_result(
            title=f"scVI DE: {group1} vs {group2}",
            operation="scvi_differential_expression",
            parameters=parameters,
            description="Differential expression evaluated with the supplied trained scVI model.",
            warnings=[],
            input_cells=input_cells,
            input_genes=len(model.var_names),
            started_at=started_at,
            table=science.pd.DataFrame(table),
        )
        return io.NodeOutput(result)


DIFFERENTIAL_NODE_CLASSES = [
    OpenBioSingleCellPseudobulk,
    OpenBioSingleCellPseudobulkEdgeR,
    OpenBioSingleCellDecouplerPseudobulkContrast,
    OpenBioSingleCellPseudobulkDESeq2,
    OpenBioSingleCellSCVIDifferentialExpression,
]


__all__ = [
    "DIFFERENTIAL_NODE_CLASSES",
    "OpenBioSingleCellDecouplerPseudobulkContrast",
    "OpenBioSingleCellPseudobulk",
    "OpenBioSingleCellPseudobulkEdgeR",
    "OpenBioSingleCellPseudobulkDESeq2",
    "OpenBioSingleCellSCVIDifferentialExpression",
]
