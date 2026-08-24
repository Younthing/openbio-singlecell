from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import figure_to_png, finish_adata, make_plot_result, make_table_result
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, PlotResultType, TableResultType

if TYPE_CHECKING:
    from anndata import AnnData


PRIORITY_CATEGORY = "openbio/single-cell/cell-prioritization"
VISUALIZATION_CATEGORY = "openbio/single-cell/visualization"
AUGUR_RESULT_TABLES = ["summary_metrics", "full_results", "feature_importances"]


def _require_pertpy() -> Any:
    try:
        import pertpy
    except (ImportError, OSError) as error:
        raise RuntimeError("Augur analysis requires the pertpy package.") from error
    return pertpy


class OpenBioSingleCellAugur(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Augur source",
        include_raw=True,
        layer_default="log1p_norm",
    )

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
                io.String.Input("control", default=""),
                io.String.Input("treatment", default=""),
                io.Combo.Input(
                    "model",
                    options=[
                        "random_forest_classifier",
                        "logistic_regression_classifier",
                        "random_forest_regressor",
                    ],
                    default="random_forest_classifier",
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.Int.Input("subsample_size", default=50, min=2, max=2**31 - 1),
                io.Boolean.Input("select_variance_features", default=False),
                io.Float.Input("span", default=0.75, min=0.0, max=1.0, step=0.05, advanced=True),
                io.Int.Input("n_threads", default=1, min=1, max=1024, advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
                io.String.Input("result_key", default="augurpy_results", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        cell_type_key: str = "cell_type",
        condition_key: str = "group",
        control: str = "",
        treatment: str = "",
        model: str = "random_forest_classifier",
        source: DynamicExpressionSource | None = None,
        subsample_size: int = 50,
        select_variance_features: bool = False,
        span: float = 0.75,
        n_threads: int = 1,
        random_seed: int = 123,
        result_key: str = "augurpy_results",
    ) -> io.NodeOutput:
        if cell_type_key not in adata.obs:
            raise ValueError(f"Augur cell type column not found in obs: {cell_type_key!r}")
        if condition_key not in adata.obs:
            raise ValueError(f"Augur condition column not found in obs: {condition_key!r}")
        control = control.strip()
        treatment = treatment.strip()
        if not control or not treatment:
            raise ValueError("Augur requires both control and treatment labels.")
        if control == treatment:
            raise ValueError("Augur control and treatment labels must be different.")
        result_key = result_key.strip()
        if not result_key:
            raise ValueError("Augur result key cannot be empty.")

        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        if expression.kind == "raw":
            analysis_adata = adata.raw.to_adata()
            analysis_adata.obs = adata.obs.copy()
        elif expression.kind == "layer":
            analysis_adata = adata.copy()
            analysis_adata.X = analysis_adata.layers[expression.layer_name].copy()
        else:
            analysis_adata = adata

        pertpy = _require_pertpy()
        started_at = time.perf_counter()
        augur = pertpy.tl.Augur(model)
        loaded = augur.load(
            analysis_adata,
            label_col=condition_key,
            cell_type_col=cell_type_key,
            condition_label=control,
            treatment_label=treatment,
        )
        output, _ = augur.predict(
            loaded,
            subsample_size=subsample_size,
            n_threads=n_threads,
            select_variance_features=select_variance_features,
            span=span,
            key_added=result_key,
            random_state=random_seed,
        )
        parameters = {
            "cell_type_key": cell_type_key,
            "condition_key": condition_key,
            "control": control,
            "treatment": treatment,
            "model": model,
            **expression.parameters(),
            "subsample_size": subsample_size,
            "select_variance_features": select_variance_features,
            "span": span,
            "n_threads": n_threads,
            "random_seed": random_seed,
            "result_key": result_key,
        }
        finish_adata(
            output,
            "augur",
            parameters,
            int(adata.n_obs),
            int(adata.n_vars),
            started_at,
            random_seed=random_seed,
        )
        return io.NodeOutput(output)


class OpenBioSingleCellAugurResults(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAugurResults",
            display_name="Augur Results",
            category=PRIORITY_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Combo.Input("result_table", options=AUGUR_RESULT_TABLES, default="summary_metrics"),
                io.String.Input("result_key", default="augurpy_results", advanced=True),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        result_table: str = "summary_metrics",
        result_key: str = "augurpy_results",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        result_key = result_key.strip()
        stored = adata.uns.get(result_key)
        if not isinstance(stored, dict):
            raise ValueError(f"Augur results not found in uns[{result_key!r}].")
        if result_table not in stored:
            raise ValueError(f"Augur result table {result_table!r} not found in uns[{result_key!r}].")
        started_at = time.perf_counter()
        table = science.pd.DataFrame(stored[result_table]).reset_index()
        result = make_table_result(
            title=f"Augur {result_table.replace('_', ' ')}",
            operation="augur_results",
            parameters={"result_table": result_table, "result_key": result_key},
            description="Stored cell-type prioritization results from Augur.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
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
            outputs=[PlotResultType.Output(display_name="plot")],
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
        result = make_plot_result(
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
    OpenBioSingleCellAugurResults,
    OpenBioSingleCellCellTypeCorrelation,
]


__all__ = [
    "OpenBioSingleCellAugur",
    "OpenBioSingleCellAugurResults",
    "OpenBioSingleCellCellTypeCorrelation",
    "POPULATION_NODE_CLASSES",
]
