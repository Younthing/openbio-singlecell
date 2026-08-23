from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import figure_to_png, make_result
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


DIFFERENTIAL_CATEGORY = "openbio/single-cell/differential-expression"
PLOT_CATEGORY = "openbio/single-cell/visualization"
MARKER_COLUMNS = [
    "group",
    "gene",
    "score",
    "logFC",
    "p",
    "p_adj",
    "rank",
    "pct_in_group",
    "pct_rest",
]


def _first_column(frame: Any, names: tuple[str, ...]) -> Any:
    for name in names:
        if name in frame.columns:
            return frame[name]
    return dependencies.pd.Series(dependencies.np.nan, index=frame.index)


class OpenBioSingleCellMarkerGenes(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerGenes",
            display_name="Marker Genes",
            category=DIFFERENTIAL_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.Combo.Input(
                    "method",
                    options=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"],
                    default="wilcoxon",
                ),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="X"),
                io.String.Input("layer_name", default="", advanced=True),
                io.Int.Input("n_genes", default=100, min=1, max=2**31 - 1),
                io.Boolean.Input("pts", default=True, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        groupby: str = "leiden",
        method: str = "wilcoxon",
        source: str = "X",
        layer_name: str = "",
        n_genes: int = 100,
        pts: bool = True,
        random_seed: int = 0,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        if groupby not in adata.obs:
            raise ValueError(f"Marker groupby column not found in obs: {groupby!r}")
        if source == "raw" and adata.raw is None:
            raise ValueError("Marker source 'raw' was selected, but adata.raw is unavailable.")
        if source == "layer" and layer_name not in adata.layers:
            raise ValueError(f"Marker layer not found: {layer_name!r}")

        started_at = time.perf_counter()
        cells = int(adata.n_obs)
        genes = int(adata.raw.n_vars) if source == "raw" else int(adata.n_vars)
        work = adata.copy()
        if not isinstance(work.obs[groupby].dtype, dependencies.pd.CategoricalDtype):
            work.obs[groupby] = dependencies.pd.Categorical(work.obs[groupby].astype(str))

        key = "_openbio_marker_genes"
        use_raw = source == "raw"
        layer = layer_name if source == "layer" else None
        method_options = {"random_state": random_seed} if method == "logreg" else {}
        dependencies.sc.tl.rank_genes_groups(
            work,
            groupby=groupby,
            method=method,
            n_genes=n_genes,
            use_raw=use_raw,
            layer=layer,
            pts=pts,
            key_added=key,
            **method_options,
        )
        ranked = dependencies.sc.get.rank_genes_groups_df(work, group=None, key=key)
        table = dependencies.pd.DataFrame(index=ranked.index)
        table["group"] = _first_column(ranked, ("group",)).astype(str)
        table["gene"] = _first_column(ranked, ("names", "gene")).astype(str)
        table["score"] = _first_column(ranked, ("scores", "score"))
        table["logFC"] = _first_column(ranked, ("logfoldchanges", "logFC"))
        table["p"] = _first_column(ranked, ("pvals", "p"))
        table["p_adj"] = _first_column(ranked, ("pvals_adj", "p_adj"))
        table["rank"] = table.groupby("group", sort=False).cumcount() + 1
        table["pct_in_group"] = _first_column(ranked, ("pct_nz_group", "pct_in_group"))
        table["pct_rest"] = _first_column(ranked, ("pct_nz_reference", "pct_rest"))
        table = table[MARKER_COLUMNS].reset_index(drop=True)

        parameters = {
            "groupby": groupby,
            "method": method,
            "source": source,
            "layer_name": layer_name,
            "n_genes": n_genes,
            "pts": pts,
            "random_seed": random_seed,
        }
        result = make_result(
            kind="table",
            title=f"Marker genes by {groupby}",
            operation="marker_genes",
            parameters=parameters,
            description="Ranked marker genes with a stable output column contract.",
            warnings=[],
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellUMAPPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellUMAPPlot",
            display_name="UMAP Plot",
            category=PLOT_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("color", default="leiden"),
                io.Float.Input("point_size", default=10.0, min=0.1, max=1000.0, step=1.0),
                io.String.Input("color_map", default="viridis", advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        color: str = "leiden",
        point_size: float = 10.0,
        color_map: str = "viridis",
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        if "X_umap" not in adata.obsm:
            raise ValueError("UMAP coordinates not found in obsm['X_umap']; run UMAP first.")
        coordinates = dependencies.np.asarray(adata.obsm["X_umap"])
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError("UMAP coordinates must have at least two columns.")
        if color and color not in adata.obs:
            raise ValueError(f"UMAP color column not found in obs: {color!r}")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        warnings = []
        figure = dependencies.Figure(figsize=(8, 7), constrained_layout=True)
        axis = figure.subplots()

        if not color:
            axis.scatter(coordinates[:, 0], coordinates[:, 1], s=point_size, alpha=0.8, color="#246bfe")
        else:
            values = adata.obs[color]
            if dependencies.pd.api.types.is_numeric_dtype(values.dtype):
                numeric = dependencies.np.asarray(values, dtype=float)
                points = axis.scatter(
                    coordinates[:, 0],
                    coordinates[:, 1],
                    c=numeric,
                    cmap=color_map,
                    s=point_size,
                    alpha=0.8,
                )
                figure.colorbar(points, ax=axis, label=color)
            else:
                categorical = dependencies.pd.Categorical(values.astype(str))
                categories = list(categorical.categories)
                for index, category in enumerate(categories):
                    mask = dependencies.np.asarray(categorical.codes == index)
                    axis.scatter(
                        coordinates[mask, 0],
                        coordinates[mask, 1],
                        s=point_size,
                        alpha=0.8,
                        label=str(category),
                    )
                if len(categories) <= 20:
                    axis.legend(title=color, bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=1.5)
                else:
                    warnings.append(f"Legend omitted because {color!r} contains more than 20 categories.")

        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        axis.set_title(f"UMAP colored by {color}" if color else "UMAP")
        png = figure_to_png(figure)
        parameters = {"color": color, "point_size": point_size, "color_map": color_map}
        result = make_result(
            kind="plot",
            title=axis.get_title(),
            operation="umap_plot",
            parameters=parameters,
            description="UMAP embedding rendered from arrays without modifying AnnData.",
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            png=png,
        )
        return io.NodeOutput(result)


RESULT_NODE_CLASSES = [OpenBioSingleCellMarkerGenes, OpenBioSingleCellUMAPPlot]
