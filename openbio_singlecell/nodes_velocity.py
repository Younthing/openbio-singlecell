from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from .analysis_utils import make_plot_result, make_summary_result, make_table_result
from .node_types import (
    AnnDataType,
    PlotResultType,
    SummaryResultType,
    TableResultType,
    VelocityStateType,
)
from .velocity_analysis import (
    VelocityState,
    run_velocity_estimate,
    run_velocity_graph,
    run_velocity_moments,
    run_velocity_prepare,
    run_velocity_ranking,
    run_velocity_recover,
    run_velocity_stream,
    validate_velocity_state,
    velocity_code,
)

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/velocity"


def _summary_result(
    summary: dict,
    *,
    title: str,
    operation: str,
    started_at: float,
    input_cells: int | None = None,
    input_genes: int | None = None,
):
    key_results = summary["key_results"]
    cells = int(key_results.get("input_cells", input_cells))
    genes = int(key_results.get("input_genes", input_genes))
    return make_summary_result(
        summary=summary,
        title=title,
        operation=operation,
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )


class OpenBioSingleCellVelocityFilterAndNormalize(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityFilterAndNormalize",
            display_name="Prepare Velocity Abundances",
            category=CATEGORY,
            description=(
                "Validate raw spliced/unspliced counts, apply explicit shared-gene filters, and normalize only "
                "canonical velocity layers into a fingerprinted staged artifact."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("spliced_layer", default="spliced"),
                io.String.Input("unspliced_layer", default="unspliced"),
                io.Int.Input("min_shared_counts", default=20, min=0, max=2**31 - 1),
                io.Int.Input(
                    "min_shared_cells", default=0, min=0, max=2**31 - 1, advanced=True
                ),
                io.Combo.Input(
                    "normalization_target",
                    options=["median_library"],
                    default="median_library",
                    advanced=True,
                ),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                VelocityStateType.Output(display_name="velocity_state"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        spliced_layer: str = "spliced",
        unspliced_layer: str = "unspliced",
        min_shared_counts: int = 20,
        min_shared_cells: int = 0,
        normalization_target: str = "median_library",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        parameters = {
            "spliced_layer": spliced_layer,
            "unspliced_layer": unspliced_layer,
            "min_shared_counts": min_shared_counts,
            "min_shared_cells": min_shared_cells,
            "normalization_target": normalization_target,
            "overwrite_existing": overwrite_existing,
        }
        result = run_velocity_prepare(adata, **parameters)
        summary = result.summary
        report = _summary_result(
            summary,
            title="Velocity abundance preparation summary",
            operation="prepare_velocity_abundances",
            started_at=started_at,
        )
        return io.NodeOutput(result, report, velocity_code("prepare", **parameters))


class OpenBioSingleCellVelocityMoments(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityMoments",
            display_name="Velocity Moments",
            category=CATEGORY,
            description="Compute Ms/Mu from one exact named graph without hidden PCA or neighbor construction.",
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.String.Input("neighbors_key", default="neighbors"),
                io.Combo.Input(
                    "mode",
                    options=["connectivities", "distances"],
                    default="connectivities",
                    advanced=True,
                ),
                io.Float.Input("max_dense_gib", default=2.0, min=1e-12, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                VelocityStateType.Output(display_name="velocity_state"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        velocity_state: VelocityState,
        neighbors_key: str = "neighbors",
        mode: str = "connectivities",
        max_dense_gib: float = 2.0,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        parameters = {
            "neighbors_key": neighbors_key,
            "mode": mode,
            "max_dense_gib": max_dense_gib,
            "overwrite_existing": overwrite_existing,
        }
        result = run_velocity_moments(velocity_state, **parameters)
        summary = result.summary
        report = _summary_result(
            summary,
            title="Velocity moments summary",
            operation="velocity_moments",
            started_at=started_at,
        )
        return io.NodeOutput(result, report, velocity_code("moments", **parameters))


class OpenBioSingleCellEstimateVelocity(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellEstimateVelocity",
            display_name="Estimate RNA Velocity",
            category=CATEGORY,
            description=(
                "Fit one explicit deterministic, stochastic, or dynamical velocity model on a validated stage."
            ),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.Combo.Input(
                    "mode",
                    options=["stochastic", "deterministic", "dynamical"],
                    default="stochastic",
                ),
                io.String.Input("vkey", default="velocity", advanced=True),
                io.Float.Input("min_r2", default=0.01, min=0.0, advanced=True),
                io.Float.Input("min_likelihood", default=0.001, min=0.0, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                VelocityStateType.Output(display_name="velocity_state"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        velocity_state: VelocityState,
        mode: str = "stochastic",
        vkey: str = "velocity",
        min_r2: float = 0.01,
        min_likelihood: float = 0.001,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        parameters = {
            "mode": mode,
            "vkey": vkey,
            "min_r2": min_r2,
            "min_likelihood": min_likelihood,
            "overwrite_existing": overwrite_existing,
        }
        result = run_velocity_estimate(velocity_state, **parameters)
        summary = result.summary
        report = _summary_result(
            summary,
            title="RNA velocity estimate summary",
            operation="estimate_rna_velocity",
            started_at=started_at,
        )
        return io.NodeOutput(result, report, velocity_code("estimate", **parameters))


class OpenBioSingleCellVelocityGraph(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityGraph",
            display_name="Velocity Directed-Correlation Graph",
            category=CATEGORY,
            description=(
                "Convert one verified velocity estimate into positive/negative directed cosine-alignment evidence."
            ),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.String.Input("vkey", default="velocity", advanced=True),
                io.Combo.Input("xkey", options=["Ms"], default="Ms", advanced=True),
                io.Combo.Input(
                    "mode_neighbors",
                    options=["distances", "connectivities"],
                    default="distances",
                    advanced=True,
                ),
                io.Int.Input("n_jobs", default=1, min=1, max=1024, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                VelocityStateType.Output(display_name="velocity_state"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        velocity_state: VelocityState,
        vkey: str = "velocity",
        xkey: str = "Ms",
        mode_neighbors: str = "distances",
        n_jobs: int = 1,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        parameters = {
            "vkey": vkey,
            "xkey": xkey,
            "mode_neighbors": mode_neighbors,
            "n_jobs": n_jobs,
            "overwrite_existing": overwrite_existing,
        }
        result = run_velocity_graph(velocity_state, **parameters)
        summary = result.summary
        report = _summary_result(
            summary,
            title="Velocity directed-correlation graph summary",
            operation="build_velocity_graph",
            started_at=started_at,
        )
        return io.NodeOutput(result, report, velocity_code("graph", **parameters))


class OpenBioSingleCellRecoverDynamics(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRecoverDynamics",
            display_name="Recover RNA Velocity Dynamics",
            category=CATEGORY,
            description=(
                "Recover reusable dynamical splicing parameters for one explicit, deterministic gene set."
            ),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.Combo.Input(
                    "gene_selection",
                    options=["velocity_genes", "all"],
                    default="velocity_genes",
                ),
                io.Int.Input("n_top_genes", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("max_iter", default=10, min=1, max=10_000, advanced=True),
                io.Int.Input("n_jobs", default=1, min=1, max=1024, advanced=True),
                io.Float.Input("max_dense_gib", default=2.0, min=1e-12, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                VelocityStateType.Output(display_name="velocity_state"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        velocity_state: VelocityState,
        gene_selection: str = "velocity_genes",
        n_top_genes: int = 0,
        max_iter: int = 10,
        n_jobs: int = 1,
        max_dense_gib: float = 2.0,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        parameters = {
            "gene_selection": gene_selection,
            "n_top_genes": n_top_genes,
            "max_iter": max_iter,
            "n_jobs": n_jobs,
            "max_dense_gib": max_dense_gib,
            "overwrite_existing": overwrite_existing,
        }
        result = run_velocity_recover(velocity_state, **parameters)
        summary = result.summary
        report = _summary_result(
            summary,
            title="Recovered RNA velocity dynamics summary",
            operation="recover_velocity_dynamics",
            started_at=started_at,
        )
        return io.NodeOutput(result, report, velocity_code("recover", **parameters))


class OpenBioSingleCellVelocityGeneRanking(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityGeneRanking",
            display_name="Recovered Dynamics Fit Ranking",
            category=CATEGORY,
            description=(
                "Pure bounded view of verified recovered fits, ranked by fit likelihood without inferential claims."
            ),
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.Int.Input("top_n", default=50, min=1, max=2**31 - 1),
                io.Boolean.Input("include_failed", default=False, advanced=True),
                io.Int.Input(
                    "max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True
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
        velocity_state: VelocityState,
        top_n: int = 50,
        include_failed: bool = False,
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        adata, _, _, _ = validate_velocity_state(velocity_state)
        started_at = time.perf_counter()
        parameters = {
            "top_n": top_n,
            "include_failed": include_failed,
            "max_output_rows": max_output_rows,
        }
        table, summary = run_velocity_ranking(velocity_state, **parameters)
        result = make_table_result(
            title="Recovered dynamics fit ranking",
            operation="rank_recovered_dynamics",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        report = _summary_result(
            summary,
            title="Recovered dynamics fit ranking summary",
            operation="rank_recovered_dynamics",
            started_at=started_at,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
        )
        return io.NodeOutput(result, report, velocity_code("ranking", **parameters))


class OpenBioSingleCellVelocityStreamPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellVelocityStreamPlot",
            display_name="Velocity Stream Plot",
            category=CATEGORY,
            description="Read-only visualization adapter over one verified velocity-graph artifact.",
            inputs=[
                VelocityStateType.Input("velocity_state"),
                io.String.Input("basis", default="umap"),
                io.String.Input("color_key", default="leiden"),
                io.Float.Input("density", default=2.0, min=1e-12, advanced=True),
                io.Float.Input("smooth", default=0.5, min=1e-12, advanced=True),
                io.Float.Input("min_mass", default=1.0, min=0.0, advanced=True),
            ],
            outputs=[
                PlotResultType.Output(display_name="plot"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        velocity_state: VelocityState,
        basis: str = "umap",
        color_key: str = "leiden",
        density: float = 2.0,
        smooth: float = 0.5,
        min_mass: float = 1.0,
    ) -> io.NodeOutput:
        adata, _, _, _ = validate_velocity_state(
            velocity_state, allowed_stages=("velocity_graph",)
        )
        started_at = time.perf_counter()
        parameters = {
            "basis": basis,
            "color_key": color_key,
            "density": density,
            "smooth": smooth,
            "min_mass": min_mass,
        }
        png, summary = run_velocity_stream(velocity_state, **parameters)
        result = make_plot_result(
            title=f"RNA velocity stream on {basis}",
            operation="render_velocity_stream",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            png=png,
        )
        report = _summary_result(
            summary,
            title="Velocity stream visualization summary",
            operation="render_velocity_stream",
            started_at=started_at,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
        )
        return io.NodeOutput(result, report, velocity_code("stream", **parameters))


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
