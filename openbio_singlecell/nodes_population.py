from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import figure_to_png, make_result
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


PRIORITY_CATEGORY = "openbio/single-cell/cell-prioritization"
VISUALIZATION_CATEGORY = "openbio/single-cell/visualization"


def _require_pertpy() -> Any:
    try:
        import pertpy
    except (ImportError, OSError) as error:
        raise RuntimeError("Augur analysis requires the pertpy package.") from error
    return pertpy


class OpenBioSingleCellAugur(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAugur",
            display_name="Augur Cell Prioritization",
            category=PRIORITY_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("cell_type_key", default="cell_type"),
                io.String.Input("condition_key", default="group"),
                io.String.Input("control", default="control"),
                io.String.Input("treatment", default="treatment"),
                io.String.Input("model", default="random_forest_classifier"),
                io.Combo.Input(
                    "result_table",
                    options=["summary_metrics", "full_results", "feature_importances"],
                    default="summary_metrics",
                ),
                io.Int.Input("subsample_size", default=50, min=2, max=2**31 - 1),
                io.Boolean.Input("select_variance_features", default=False),
                io.Float.Input("span", default=0.75, min=0.0, max=1.0, step=0.05, advanced=True),
                io.Int.Input("n_threads", default=1, min=1, max=1024, advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        cell_type_key: str = "cell_type",
        condition_key: str = "group",
        control: str = "control",
        treatment: str = "treatment",
        model: str = "random_forest_classifier",
        result_table: str = "summary_metrics",
        subsample_size: int = 50,
        select_variance_features: bool = False,
        span: float = 0.75,
        n_threads: int = 1,
        random_seed: int = 123,
    ) -> io.NodeOutput:
        if cell_type_key not in adata.obs:
            raise ValueError(f"Augur cell type column not found in obs: {cell_type_key!r}")
        if condition_key not in adata.obs:
            raise ValueError(f"Augur condition column not found in obs: {condition_key!r}")

        pertpy = _require_pertpy()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        augur = pertpy.tl.Augur(model)
        loaded = augur.load(
            adata,
            label_col=condition_key,
            cell_type_col=cell_type_key,
            condition_label=control,
            treatment_label=treatment,
        )
        _, results = augur.predict(
            loaded,
            subsample_size=subsample_size,
            n_threads=n_threads,
            select_variance_features=select_variance_features,
            span=span,
            key_added="augurpy_results",
            random_state=random_seed,
        )
        table = science.pd.DataFrame(results[result_table]).reset_index()
        parameters = {
            "cell_type_key": cell_type_key,
            "condition_key": condition_key,
            "control": control,
            "treatment": treatment,
            "model": model,
            "result_table": result_table,
            "subsample_size": subsample_size,
            "select_variance_features": select_variance_features,
            "span": span,
            "n_threads": n_threads,
            "random_seed": random_seed,
        }
        result = make_result(
            kind="table",
            title=f"Augur {result_table.replace('_', ' ')}",
            operation="augur",
            parameters=parameters,
            description=f"Cell-type prioritization for {treatment} versus {control}.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellCellTypeCorrelation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellTypeCorrelation",
            display_name="Cell Type Correlation",
            category=VISUALIZATION_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="cell_type"),
                io.String.Input("use_rep", default="X_pca"),
                io.Combo.Input("cor_method", options=["pearson", "spearman", "kendall"], default="pearson"),
                io.String.Input("color_map", default="RdYlBu", advanced=True),
                io.Boolean.Input("show_numbers", default=False, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "cell_type",
        use_rep: str = "X_pca",
        cor_method: str = "pearson",
        color_map: str = "RdYlBu",
        show_numbers: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        if groupby not in adata.obs:
            raise ValueError(f"Correlation group column not found in obs: {groupby!r}")
        if use_rep and use_rep not in adata.obsm:
            raise ValueError(f"Correlation representation not found in obsm: {use_rep!r}")

        started_at = time.perf_counter()
        work = adata.copy()
        science.sc.tl.dendrogram(
            work,
            groupby=groupby,
            use_rep=use_rep or None,
            use_raw=False,
            cor_method=cor_method,
        )
        figure = science.Figure(figsize=(7, 6), constrained_layout=True)
        axis = figure.subplots()
        science.sc.pl.correlation_matrix(
            work,
            groupby,
            dendrogram=False,
            show_correlation_numbers=show_numbers,
            cmap=color_map,
            show=False,
            ax=axis,
        )
        png = figure_to_png(figure)
        parameters = {
            "groupby": groupby,
            "use_rep": use_rep,
            "cor_method": cor_method,
            "color_map": color_map,
            "show_numbers": show_numbers,
        }
        result = make_result(
            kind="plot",
            title=f"Cell type correlation by {groupby}",
            operation="cell_type_correlation",
            parameters=parameters,
            description="Correlation matrix between cell groups.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            png=png,
        )
        return io.NodeOutput(result)


POPULATION_NODE_CLASSES = [
    OpenBioSingleCellAugur,
    OpenBioSingleCellCellTypeCorrelation,
]


__all__ = [
    "OpenBioSingleCellAugur",
    "OpenBioSingleCellCellTypeCorrelation",
    "POPULATION_NODE_CLASSES",
]
