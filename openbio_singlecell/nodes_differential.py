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
                io.String.Input("baseline", default=""),
                io.String.Input("comparison", default=""),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
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
            outputs=[SingleCellResultType.Output(display_name="result")],
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
        result = make_result(
            kind="table",
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
                io.Combo.Input("source", options=["X", "layer"], default="layer"),
                io.String.Input("layer_name", default="counts"),
                io.String.Input("groupby", default="group"),
                io.String.Input("group1", default=""),
                io.String.Input("group2", default=""),
                io.String.Input("subset_column", default=""),
                io.String.Input("subset_value", default=""),
                io.Combo.Input("gene_likelihood", options=["zinb", "nb", "poisson"], default="nb"),
                io.Int.Input("n_latent", default=10, min=1, max=4096, advanced=True),
                io.Float.Input("delta", default=0.25, min=0.0, step=0.05, advanced=True),
                io.Int.Input("max_epochs", default=0, min=0, max=100000, advanced=True),
                io.Boolean.Input("early_stopping", default=True, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: str = "layer",
        layer_name: str = "counts",
        groupby: str = "group",
        group1: str = "",
        group2: str = "",
        subset_column: str = "",
        subset_value: str = "",
        gene_likelihood: str = "nb",
        n_latent: int = 10,
        delta: float = 0.25,
        max_epochs: int = 0,
        early_stopping: bool = True,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        if source == "layer" and layer_name not in adata.layers:
            raise ValueError(f"scVI differential-expression layer not found: {layer_name!r}")
        if groupby not in adata.obs:
            raise ValueError(f"scVI differential-expression group column not found in obs: {groupby!r}")
        if not group1 or not group2:
            raise ValueError("scVI differential-expression group1 and group2 values are required.")
        if bool(subset_column) != bool(subset_value):
            raise ValueError("scVI subset_column and subset_value must be provided together.")
        if subset_column and subset_column not in adata.obs:
            raise ValueError(f"scVI subset column not found in obs: {subset_column!r}")

        try:
            import scvi
        except (ImportError, OSError) as error:
            raise RuntimeError("scVI Differential Expression requires the scvi-tools package.") from error

        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        scvi.settings.seed = random_seed
        work = adata.copy()
        scvi.model.SCVI.setup_anndata(work, layer=layer_name if source == "layer" else None)
        model = scvi.model.SCVI(work, n_latent=n_latent, gene_likelihood=gene_likelihood)
        train_kwargs = {"early_stopping": early_stopping}
        if max_epochs:
            train_kwargs["max_epochs"] = max_epochs
        model.train(**train_kwargs)

        analysis_adata = work
        if subset_column:
            mask = work.obs[subset_column].astype(str) == subset_value
            analysis_adata = work[mask].copy()
            if analysis_adata.n_obs == 0:
                raise ValueError(f"scVI subset contains no observations: {subset_column}={subset_value!r}")
        table = model.differential_expression(
            analysis_adata,
            groupby=groupby,
            group1=group1,
            group2=group2,
            delta=delta,
        ).reset_index()
        if "index" in table.columns and "gene" not in table.columns:
            table = table.rename(columns={"index": "gene"})

        parameters = {
            "source": source,
            "layer_name": layer_name,
            "groupby": groupby,
            "group1": group1,
            "group2": group2,
            "subset_column": subset_column,
            "subset_value": subset_value,
            "gene_likelihood": gene_likelihood,
            "n_latent": n_latent,
            "delta": delta,
            "max_epochs": max_epochs,
            "early_stopping": early_stopping,
            "random_seed": random_seed,
        }
        result = make_result(
            kind="table",
            title=f"scVI DE: {group1} vs {group2}",
            operation="scvi_differential_expression",
            parameters=parameters,
            description="scVI differential expression fitted and evaluated in one node.",
            warnings=[],
            input_cells=int(analysis_adata.n_obs),
            input_genes=int(analysis_adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
            table=science.pd.DataFrame(table),
        )
        return io.NodeOutput(result)


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
