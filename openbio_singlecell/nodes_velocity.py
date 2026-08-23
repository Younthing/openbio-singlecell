from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import figure_to_png, finish_adata, make_plot_result, make_table_result
from .node_types import AnnDataType, PlotResultType, TableResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/velocity"


def _require_scvelo() -> Any:
    try:
        import scvelo as scv
    except (ImportError, OSError) as error:
        raise RuntimeError(
            "RNA velocity nodes require the optional dependency scvelo. Install it with: python -m pip install scvelo"
        ) from error
    return scv


def _use_velocity_layers(
    adata: AnnData,
    *,
    spliced_layer: str,
    unspliced_layer: str,
) -> None:
    if spliced_layer not in adata.layers:
        raise ValueError(f"Spliced layer not found: {spliced_layer!r}")
    if unspliced_layer not in adata.layers:
        raise ValueError(f"Unspliced layer not found: {unspliced_layer!r}")
    if spliced_layer == unspliced_layer:
        raise ValueError("Spliced and unspliced layers must be different.")

    if spliced_layer != "spliced":
        adata.layers["spliced"] = adata.layers[spliced_layer].copy()
    if unspliced_layer != "unspliced":
        adata.layers["unspliced"] = adata.layers[unspliced_layer].copy()


class OpenBioSingleCellVelocityFilterAndNormalize(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityFilterAndNormalize",
            display_name="Velocity Filter and Normalize",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("spliced_layer", default="spliced"),
                io.String.Input("unspliced_layer", default="unspliced"),
                io.Int.Input("min_shared_counts", default=20, min=0, max=2**31 - 1),
                io.Int.Input("n_top_genes", default=2000, min=1, max=2**31 - 1),
                io.String.Input("subset_highly_variable", default="seurat"),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        spliced_layer: str = "spliced",
        unspliced_layer: str = "unspliced",
        min_shared_counts: int = 20,
        n_top_genes: int = 2000,
        subset_highly_variable: str = "seurat",
    ) -> io.NodeOutput:
        scv = _require_scvelo()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        _use_velocity_layers(
            output,
            spliced_layer=spliced_layer,
            unspliced_layer=unspliced_layer,
        )
        scv.pp.filter_and_normalize(
            output,
            min_shared_counts=min_shared_counts,
            subset_highly_variable=subset_highly_variable,
            n_top_genes=n_top_genes,
        )
        parameters = {
            "spliced_layer": spliced_layer,
            "unspliced_layer": unspliced_layer,
            "min_shared_counts": min_shared_counts,
            "n_top_genes": n_top_genes,
            "subset_highly_variable": subset_highly_variable,
        }
        finish_adata(
            output,
            "velocity_filter_and_normalize",
            parameters,
            cells,
            genes,
            started_at,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellVelocityMoments(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityMoments",
            display_name="Velocity Moments",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_pcs", default=0, min=0, max=4096),
                io.Int.Input("n_neighbors", default=0, min=0, max=4096),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        n_pcs: int = 0,
        n_neighbors: int = 0,
    ) -> io.NodeOutput:
        scv = _require_scvelo()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        scv.pp.moments(
            output,
            n_pcs=n_pcs or None,
            n_neighbors=n_neighbors or None,
        )
        parameters = {"n_pcs": n_pcs or None, "n_neighbors": n_neighbors or None}
        finish_adata(output, "velocity_moments", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellEstimateVelocity(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellEstimateVelocity",
            display_name="Estimate RNA Velocity",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input(
                    "mode",
                    options=["deterministic", "stochastic", "dynamical"],
                    default="deterministic",
                ),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        mode: str = "deterministic",
    ) -> io.NodeOutput:
        scv = _require_scvelo()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        scv.tl.velocity(output, mode=mode)
        parameters = {"mode": mode}
        finish_adata(output, "estimate_velocity", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellVelocityGraph(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityGraph",
            display_name="Velocity Graph",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_jobs", default=8, min=1, max=1024, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, n_jobs: int = 8) -> io.NodeOutput:
        scv = _require_scvelo()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        scv.tl.velocity_graph(output, n_jobs=n_jobs)
        parameters = {"n_jobs": n_jobs}
        finish_adata(output, "velocity_graph", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellRecoverDynamics(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRecoverDynamics",
            display_name="Recover Velocity Dynamics",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_jobs", default=8, min=1, max=1024, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(cls, adata: AnnData, n_jobs: int = 8) -> io.NodeOutput:
        scv = _require_scvelo()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        scv.tl.recover_dynamics(output, n_jobs=n_jobs)
        parameters = {"n_jobs": n_jobs}
        finish_adata(output, "recover_velocity_dynamics", parameters, cells, genes, started_at)
        return io.NodeOutput(output)


class OpenBioSingleCellVelocityGeneRanking(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityGeneRanking",
            display_name="RNA Velocity Gene Ranking",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("top_n", default=5, min=1, max=2**31 - 1),
                io.String.Input("likelihood_column", default="fit_likelihood", advanced=True),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        top_n: int = 5,
        likelihood_column: str = "fit_likelihood",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if likelihood_column not in adata.var:
            raise ValueError(
                f"Velocity likelihood column not found in var: {likelihood_column!r}; run Recover Dynamics first."
            )

        started_at = time.perf_counter()
        ranked = adata.var[likelihood_column].sort_values(ascending=False).head(top_n)
        table = science.pd.DataFrame(
            {
                "gene": ranked.index.astype(str),
                likelihood_column: ranked.to_numpy(),
            }
        )
        parameters = {
            "top_n": top_n,
            "likelihood_column": likelihood_column,
        }
        result = make_table_result(
            title="RNA velocity gene ranking",
            operation="velocity_gene_ranking",
            parameters=parameters,
            description="Genes ranked by recovered-dynamics fit likelihood.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellVelocityStreamPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityStreamPlot",
            display_name="Velocity Stream Plot",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("basis", default="umap"),
                io.String.Input("groupby", default="leiden"),
            ],
            outputs=[PlotResultType.Output(display_name="plot")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        basis: str = "umap",
        groupby: str = "leiden",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        scv = _require_scvelo()
        started_at = time.perf_counter()
        figure = science.Figure(figsize=(8, 6), constrained_layout=True)
        axis = figure.subplots()
        scv.pl.velocity_embedding_stream(
            adata,
            basis=basis,
            color=groupby or None,
            show=False,
            ax=axis,
        )
        png = figure_to_png(figure)
        parameters = {"basis": basis, "groupby": groupby}
        result = make_plot_result(
            title=f"RNA velocity on {basis}",
            operation="velocity_stream_plot",
            parameters=parameters,
            description="Velocity stream embedding from the computed velocity graph.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            png=png,
        )
        return io.NodeOutput(result)


VELOCITY_NODE_CLASSES = [
    OpenBioSingleCellVelocityFilterAndNormalize,
    OpenBioSingleCellVelocityMoments,
    OpenBioSingleCellEstimateVelocity,
    OpenBioSingleCellVelocityGraph,
    OpenBioSingleCellRecoverDynamics,
    OpenBioSingleCellVelocityGeneRanking,
    OpenBioSingleCellVelocityStreamPlot,
]


__all__ = ["VELOCITY_NODE_CLASSES"]
