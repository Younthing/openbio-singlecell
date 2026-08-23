from __future__ import annotations

import importlib
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import figure_to_png, finish_adata, make_result
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/cell-communication"
METHODS = ["cellphonedb", "rank_aggregate"]
PLOT_PROFILES = {
    "cellphonedb": {
        "colour": "lr_means",
        "size": "cellphone_pvals",
        "filter": "cellphone_pvals",
        "orderby": "lr_means",
        "orderby_ascending": False,
        "inverse_colour": False,
    },
    "rank_aggregate": {
        "colour": "magnitude_rank",
        "size": "specificity_rank",
        "filter": "specificity_rank",
        "orderby": "magnitude_rank",
        "orderby_ascending": True,
        "inverse_colour": True,
    },
}


def _require_liana() -> Any:
    try:
        return importlib.import_module("liana")
    except (ImportError, OSError) as exc:
        raise RuntimeError(f"The 'liana' package is required for cell communication nodes ({exc}).") from exc


def _required_name(value: str, description: str) -> str:
    name = value.strip()
    if not name:
        raise ValueError(f"{description} cannot be empty.")
    return name


def _liana_table(adata: AnnData, result_key: str) -> Any:
    science = dependencies.require_scientific_dependencies()
    result_key = _required_name(result_key, "LIANA result key")
    if result_key not in adata.uns:
        raise ValueError(f"LIANA results not found in adata.uns[{result_key!r}].")
    table = adata.uns[result_key]
    if not isinstance(table, science.pd.DataFrame):
        raise ValueError(f"adata.uns[{result_key!r}] is not a LIANA result table.")
    return table


def _labels(value: str) -> list[str] | None:
    labels = [label.strip() for label in value.split(",") if label.strip()]
    return labels or None


class OpenBioSingleCellLianaCommunication(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLianaCommunication",
            display_name="LIANA Communication",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cell_type"),
                io.Combo.Input("method", options=METHODS, default="rank_aggregate"),
                io.String.Input("resource_name", default="consensus"),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="X"),
                io.String.Input("layer_name", default="log1p_norm"),
                io.String.Input("result_key", default="liana_res", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "cell_type",
        method: str = "rank_aggregate",
        resource_name: str = "consensus",
        source: str = "X",
        layer_name: str = "log1p_norm",
        result_key: str = "liana_res",
    ) -> io.NodeOutput:
        groupby = _required_name(groupby, "LIANA groupby column")
        resource_name = _required_name(resource_name, "LIANA resource name")
        result_key = _required_name(result_key, "LIANA result key")
        if groupby not in adata.obs:
            raise ValueError(f"LIANA groupby column not found in obs: {groupby!r}")
        if source == "raw" and adata.raw is None:
            raise ValueError("LIANA expression source 'raw' was selected, but adata.raw is unavailable.")
        if source == "layer" and layer_name not in adata.layers:
            raise ValueError(f"LIANA expression layer not found: {layer_name!r}")

        liana = _require_liana()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        if source == "raw":
            work = adata.raw.to_adata()
            work.obs = adata.obs.copy()
        else:
            work = output
            if source == "layer":
                work.X = work.layers[layer_name].copy()

        run_method = getattr(liana.method, method)
        run_method(
            work,
            groupby=groupby,
            resource_name=resource_name,
            return_all_lrs=True,
            use_raw=False,
            verbose=False,
        )
        output.uns[result_key] = work.uns["liana_res"].copy(deep=True)
        parameters = {
            "groupby": groupby,
            "method": method,
            "resource_name": resource_name,
            "source": source,
            "layer_name": layer_name,
            "result_key": result_key,
            "return_all_lrs": True,
        }
        finish_adata(
            output,
            "liana_communication",
            parameters,
            cells,
            genes,
            started_at,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellLianaResults(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLianaResults",
            display_name="LIANA Results",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("result_key", default="liana_res", advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(cls, adata: AnnData, result_key: str = "liana_res") -> io.NodeOutput:
        started_at = time.perf_counter()
        table = _liana_table(adata, result_key).copy(deep=True).reset_index(drop=True)
        result = make_result(
            kind="table",
            title="LIANA communication results",
            operation="liana_results",
            parameters={"result_key": result_key},
            description="Ligand-receptor interactions produced by LIANA.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellLianaDotPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLianaDotPlot",
            display_name="LIANA Dot Plot",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("method", options=METHODS, default="rank_aggregate"),
                io.String.Input("source_labels", default=""),
                io.String.Input("target_labels", default=""),
                io.Float.Input("significance_threshold", default=0.01, min=0.0, max=1.0, step=0.01),
                io.Int.Input("top_n", default=20, min=1, max=2**31 - 1),
                io.String.Input("result_key", default="liana_res", advanced=True),
                io.Float.Input("figure_width", default=12.0, min=1.0, max=100.0, step=1.0, advanced=True),
                io.Float.Input("figure_height", default=8.0, min=1.0, max=100.0, step=1.0, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        method: str = "rank_aggregate",
        source_labels: str = "",
        target_labels: str = "",
        significance_threshold: float = 0.01,
        top_n: int = 20,
        result_key: str = "liana_res",
        figure_width: float = 12.0,
        figure_height: float = 8.0,
    ) -> io.NodeOutput:
        table = _liana_table(adata, result_key)
        liana = _require_liana()
        profile = PLOT_PROFILES[method]
        filter_column = profile["filter"]
        started_at = time.perf_counter()
        work = adata.copy()
        work.uns["liana_res"] = table
        plot = liana.pl.dotplot(
            adata=work,
            colour=profile["colour"],
            size=profile["size"],
            inverse_colour=profile["inverse_colour"],
            inverse_size=True,
            source_labels=_labels(source_labels),
            target_labels=_labels(target_labels),
            filterby=filter_column,
            filter_lambda=lambda value: value <= significance_threshold,
            orderby=profile["orderby"],
            orderby_ascending=profile["orderby_ascending"],
            top_n=top_n,
            figure_size=(figure_width, figure_height),
            size_range=(1, 5),
        )
        figure = plot if hasattr(plot, "savefig") else plot.draw()
        png = figure_to_png(figure)
        parameters = {
            "method": method,
            "source_labels": _labels(source_labels) or [],
            "target_labels": _labels(target_labels) or [],
            "significance_threshold": significance_threshold,
            "top_n": top_n,
            "result_key": result_key,
            "figure_width": figure_width,
            "figure_height": figure_height,
        }
        result = make_result(
            kind="plot",
            title=f"LIANA {method.replace('_', ' ')} dot plot",
            operation="liana_dot_plot",
            parameters=parameters,
            description="Top ligand-receptor interactions between source and target cell groups.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            png=png,
        )
        return io.NodeOutput(result)


COMMUNICATION_NODE_CLASSES = [
    OpenBioSingleCellLianaCommunication,
    OpenBioSingleCellLianaResults,
    OpenBioSingleCellLianaDotPlot,
]


__all__ = [
    "COMMUNICATION_NODE_CLASSES",
    "OpenBioSingleCellLianaCommunication",
    "OpenBioSingleCellLianaDotPlot",
    "OpenBioSingleCellLianaResults",
]
