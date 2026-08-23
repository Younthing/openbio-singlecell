from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import figure_to_png, make_result
from .contracts import SingleCellResult
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData


DIFFERENTIAL_CATEGORY = "openbio/single-cell/differential-expression"
PLOT_CATEGORY = "openbio/single-cell/visualization"
DIAGNOSTIC_CATEGORY = "openbio/single-cell/diagnostics"
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


def _first_column(frame: Any, names: tuple[str, ...], science: dependencies.ScientificDependencies) -> Any:
    for name in names:
        if name in frame.columns:
            return frame[name]
    return science.pd.Series(science.np.nan, index=frame.index)


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
                io.String.Input("layer_name", default=""),
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
        science = dependencies.require_scientific_dependencies()
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
        if not isinstance(work.obs[groupby].dtype, science.pd.CategoricalDtype):
            work.obs[groupby] = science.pd.Categorical(work.obs[groupby].astype(str))

        key = "_openbio_marker_genes"
        use_raw = source == "raw"
        layer = layer_name if source == "layer" else None
        method_options = {"random_state": random_seed} if method == "logreg" else {}
        science.sc.tl.rank_genes_groups(
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
        ranked = science.sc.get.rank_genes_groups_df(work, group=None, key=key)
        table = science.pd.DataFrame(index=ranked.index)
        table["group"] = _first_column(ranked, ("group",), science).astype(str)
        table["gene"] = _first_column(ranked, ("names", "gene"), science).astype(str)
        table["score"] = _first_column(ranked, ("scores", "score"), science)
        table["logFC"] = _first_column(ranked, ("logfoldchanges", "logFC"), science)
        table["p"] = _first_column(ranked, ("pvals", "p"), science)
        table["p_adj"] = _first_column(ranked, ("pvals_adj", "p_adj"), science)
        table["rank"] = table.groupby("group", sort=False).cumcount() + 1
        table["pct_in_group"] = _first_column(ranked, ("pct_nz_group", "pct_in_group"), science)
        table["pct_rest"] = _first_column(ranked, ("pct_nz_reference", "pct_rest"), science)
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
        science = dependencies.require_scientific_dependencies()
        if "X_umap" not in adata.obsm:
            raise ValueError("UMAP coordinates not found in obsm['X_umap']; run UMAP first.")
        coordinates = science.np.asarray(adata.obsm["X_umap"])
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError("UMAP coordinates must have at least two columns.")
        if color and color not in adata.obs:
            raise ValueError(f"UMAP color column not found in obs: {color!r}")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        warnings = []
        figure = science.Figure(figsize=(8, 7), constrained_layout=True)
        axis = figure.subplots()

        if not color:
            axis.scatter(coordinates[:, 0], coordinates[:, 1], s=point_size, alpha=0.8, color="#246bfe")
        else:
            values = adata.obs[color]
            if science.pd.api.types.is_numeric_dtype(values.dtype):
                numeric = science.np.asarray(values, dtype=float)
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
                categorical = science.pd.Categorical(values.astype(str))
                categories = list(categorical.categories)
                for index, category in enumerate(categories):
                    mask = science.np.asarray(categorical.codes == index)
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


class OpenBioSingleCellFilterMarkerGenes(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterMarkerGenes",
            display_name="Filter Marker Genes",
            category=DIFFERENTIAL_CATEGORY,
            inputs=[
                SingleCellResultType.Input("result"),
                io.Float.Input("min_logfc", default=1.0, step=0.1),
                io.Float.Input("min_pct_in_group", default=0.25, min=0.0, max=1.0, step=0.05),
                io.Float.Input("max_pct_rest", default=0.5, min=0.0, max=1.0, step=0.05),
                io.Float.Input("max_p_adj", default=0.05, min=0.0, max=1.0, step=0.01),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        result: SingleCellResult,
        min_logfc: float = 1.0,
        min_pct_in_group: float = 0.25,
        max_pct_rest: float = 0.5,
        max_p_adj: float = 0.05,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        required = {"logFC", "pct_in_group", "pct_rest", "p_adj"}
        if result.kind != "table" or not isinstance(result.table, science.pd.DataFrame):
            raise ValueError("Filter Marker Genes requires a table result from Marker Genes.")
        missing = sorted(required.difference(result.table.columns))
        if missing:
            raise ValueError(f"Marker result is missing columns: {missing}")

        started_at = time.perf_counter()
        table = result.table.copy()
        table = table[
            (science.pd.to_numeric(table["logFC"], errors="coerce") >= min_logfc)
            & (science.pd.to_numeric(table["pct_in_group"], errors="coerce") >= min_pct_in_group)
            & (science.pd.to_numeric(table["pct_rest"], errors="coerce") <= max_pct_rest)
            & (science.pd.to_numeric(table["p_adj"], errors="coerce") <= max_p_adj)
        ].reset_index(drop=True)
        parameters = {
            "min_logfc": min_logfc,
            "min_pct_in_group": min_pct_in_group,
            "max_pct_rest": max_pct_rest,
            "max_p_adj": max_p_adj,
        }
        filtered = make_result(
            kind="table",
            title=f"Filtered {result.title}",
            operation="filter_marker_genes",
            parameters=parameters,
            description="Marker genes filtered by fold change, prevalence, and adjusted p-value.",
            warnings=[],
            input_cells=result.input_cells,
            input_genes=result.input_genes,
            started_at=started_at,
            random_seed=result.random_seed,
            table=table,
        )
        return io.NodeOutput(filtered)


class OpenBioSingleCellMarkerExpressionPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerExpressionPlot",
            display_name="Marker Expression Plot",
            category=PLOT_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("genes", default=""),
                io.String.Input("groupby", default="leiden"),
                io.Combo.Input(
                    "plot_type", options=["dotplot", "matrixplot", "tracksplot", "violin"], default="dotplot"
                ),
                io.Combo.Input("source", options=["X", "raw", "layer"], default="raw"),
                io.String.Input("layer_name", default="log1p_norm"),
                io.Combo.Input("standard_scale", options=["none", "var", "group"], default="var", advanced=True),
                io.Boolean.Input("dendrogram", default=False, advanced=True),
                io.Boolean.Input("log", default=False, advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        genes: str = "",
        groupby: str = "leiden",
        plot_type: str = "dotplot",
        source: str = "raw",
        layer_name: str = "log1p_norm",
        standard_scale: str = "var",
        dendrogram: bool = False,
        log: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        gene_names = list(dict.fromkeys(gene.strip() for gene in genes.split(",") if gene.strip()))
        if not gene_names:
            raise ValueError("Marker Expression Plot requires at least one comma-separated gene.")
        if groupby not in adata.obs:
            raise ValueError(f"Marker plot group column not found in obs: {groupby!r}")
        expression = adata.raw.var_names if source == "raw" and adata.raw is not None else adata.var_names
        if source == "raw" and adata.raw is None:
            raise ValueError("Marker plot source 'raw' was selected, but adata.raw is unavailable.")
        if source == "layer" and layer_name not in adata.layers:
            raise ValueError(f"Marker plot layer not found: {layer_name!r}")
        missing = [gene for gene in gene_names if gene not in expression]
        if missing:
            raise ValueError(f"Marker plot genes not found in the selected expression source: {missing}")

        started_at = time.perf_counter()
        plot_kwargs = {
            "use_raw": source == "raw",
            "layer": layer_name if source == "layer" else None,
            "log": log,
            "show": False,
        }
        if plot_type == "dotplot":
            plot = science.sc.pl.dotplot(
                adata,
                gene_names,
                groupby,
                standard_scale=None if standard_scale == "none" else standard_scale,
                dendrogram=dendrogram,
                return_fig=True,
                **plot_kwargs,
            )
            plot.make_figure()
            figure = plot.fig
        elif plot_type == "matrixplot":
            plot = science.sc.pl.matrixplot(
                adata,
                gene_names,
                groupby,
                standard_scale=None if standard_scale == "none" else standard_scale,
                dendrogram=dendrogram,
                return_fig=True,
                **plot_kwargs,
            )
            plot.make_figure()
            figure = plot.fig
        elif plot_type == "tracksplot":
            axes = science.sc.pl.tracksplot(
                adata,
                gene_names,
                groupby,
                dendrogram=dendrogram,
                **plot_kwargs,
            )
            figure = next(iter(axes.values())).figure
        else:
            axis = science.sc.pl.violin(
                adata,
                gene_names,
                groupby=groupby,
                multi_panel=len(gene_names) > 1,
                stripplot=False,
                **plot_kwargs,
            )
            figure = axis.fig if hasattr(axis, "fig") else axis.figure

        parameters = {
            "genes": gene_names,
            "groupby": groupby,
            "plot_type": plot_type,
            "source": source,
            "layer_name": layer_name,
            "standard_scale": standard_scale,
            "dendrogram": dendrogram,
            "log": log,
        }
        plotted = make_result(
            kind="plot",
            title=f"{plot_type}: {', '.join(gene_names)} by {groupby}",
            operation="marker_expression_plot",
            parameters=parameters,
            description="Marker expression rendered with Scanpy without modifying AnnData.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            png=figure_to_png(figure),
        )
        return io.NodeOutput(plotted)


class OpenBioSingleCellPCAMetadataAssociations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPCAMetadataAssociations",
            display_name="PCA Metadata Associations",
            category=DIAGNOSTIC_CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("use_rep", default="X_pca"),
                io.String.Input("obs_keys", default=""),
                io.Float.Input("alpha", default=0.05, min=0.0, max=1.0, step=0.01),
                io.String.Input("p_adjust_method", default="fdr_bh", advanced=True),
            ],
            outputs=[SingleCellResultType.Output(display_name="result")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        use_rep: str = "X_pca",
        obs_keys: str = "",
        alpha: float = 0.05,
        p_adjust_method: str = "fdr_bh",
    ) -> io.NodeOutput:
        use_rep = use_rep.strip()
        if not use_rep:
            raise ValueError("PCA metadata association representation cannot be empty.")
        if use_rep not in adata.obsm:
            raise ValueError(f"PCA metadata association representation not found in obsm: {use_rep!r}")
        selected_obs = list(dict.fromkeys(key.strip() for key in obs_keys.split(",") if key.strip()))
        if not selected_obs:
            raise ValueError("PCA metadata associations require at least one comma-separated obs column.")
        missing_obs = [key for key in selected_obs if key not in adata.obs]
        if missing_obs:
            raise ValueError(f"PCA metadata association columns not found in obs: {missing_obs}")
        try:
            import decoupler
        except (ImportError, OSError) as error:
            raise RuntimeError("PCA Metadata Associations requires the decoupler package.") from error

        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        table = decoupler.get_metadata_associations(
            adata,
            obs_keys=selected_obs,
            obsm_key=use_rep,
            inplace=False,
            alpha=alpha,
            method=p_adjust_method,
            verbose=False,
        )
        table = science.pd.DataFrame(table).reset_index()
        parameters = {
            "use_rep": use_rep,
            "obs_keys": selected_obs,
            "alpha": alpha,
            "p_adjust_method": p_adjust_method,
        }
        result = make_result(
            kind="table",
            title=f"Metadata associations with {use_rep}",
            operation="pca_metadata_associations",
            parameters=parameters,
            description="ANOVA associations between an embedding and selected observation metadata.",
            warnings=[],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


RESULT_NODE_CLASSES = [
    OpenBioSingleCellMarkerGenes,
    OpenBioSingleCellUMAPPlot,
    OpenBioSingleCellFilterMarkerGenes,
    OpenBioSingleCellMarkerExpressionPlot,
    OpenBioSingleCellPCAMetadataAssociations,
]
