from __future__ import annotations

import time
from typing import TYPE_CHECKING

from comfy_api.latest import io

from .analysis_utils import make_plot_result, make_summary_result, make_table_result
from .augur import (
    AUGUR_CLASSIFIERS,
    AUGUR_VIEWS,
    AugurResult,
    augur_code,
    augur_results_code,
    run_augur_analysis,
    select_augur_view,
)
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, AugurResultType, PlotResultType, SummaryResultType, TableResultType
from .population_correlation import analyze_population_centroid_correlation, population_correlation_code

if TYPE_CHECKING:
    from anndata import AnnData


PRIORITY_CATEGORY = "openbio/single-cell/cell-prioritization"
VISUALIZATION_CATEGORY = "openbio/single-cell/visualization"


class OpenBioSingleCellAugur(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Augur count source",
        default="X",
        include_raw=True,
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAugur",
            display_name="Augur Cell Prioritization",
            category=PRIORITY_CATEGORY,
            description=(
                "Prioritize populations with audited two-Condition Pertpy Augur classifier cross-validation."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("population_key", default="cell_type"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("control", default=""),
                io.String.Input("treatment", default=""),
                io.Combo.Input(
                    "classifier",
                    options=list(AUGUR_CLASSIFIERS),
                    default="random_forest_classifier",
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                ),
                io.String.Input("technical_batch_key", default="", advanced=True),
                io.Int.Input("n_subsamples", default=50, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("subsample_size", default=20, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("folds", default=3, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("n_threads", default=1, min=1, max=1024, advanced=True),
                io.Int.Input("random_seed", default=123, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input(
                    "max_result_rows",
                    default=10_000_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
                io.Float.Input(
                    "max_result_mib",
                    default=1024.0,
                    min=1.0,
                    max=1_048_576.0,
                    step=1.0,
                    advanced=True,
                ),
            ],
            outputs=[
                AugurResultType.Output(display_name="result"),
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
        control: str = "",
        treatment: str = "",
        classifier: str = "random_forest_classifier",
        source: DynamicExpressionSource | None = None,
        annotation_status: str = "unknown",
        technical_batch_key: str = "",
        n_subsamples: int = 50,
        subsample_size: int = 20,
        folds: int = 3,
        n_threads: int = 1,
        random_seed: int = 123,
        max_result_rows: int = 10_000_000,
        max_result_mib: float = 1024.0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        result = run_augur_analysis(
            adata,
            sample_key=sample_key,
            population_key=population_key,
            condition_key=condition_key,
            control=control,
            treatment=treatment,
            classifier=classifier,
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            annotation_status=annotation_status,
            technical_batch_key=technical_batch_key,
            n_subsamples=n_subsamples,
            subsample_size=subsample_size,
            folds=folds,
            n_threads=n_threads,
            random_seed=random_seed,
            max_result_rows=max_result_rows,
            max_result_mib=max_result_mib,
        )
        summary = result.summary
        parameters = dict(summary["parameters"])
        report = make_summary_result(
            summary=summary,
            title="Augur exploratory population priorities",
            operation="augur",
            parameters=parameters,
            description=str(summary["results"]),
            warnings=[str(warning) for warning in summary["warnings"]],
            input_cells=int(adata.n_obs),
            input_genes=int(summary["key_results"]["selected_source_features"]),
            started_at=started_at,
            random_seed=random_seed,
        )
        code = augur_code(**parameters)
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellAugurResults(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAugurResults",
            display_name="Augur Results",
            category=PRIORITY_CATEGORY,
            inputs=[
                AugurResultType.Input("result"),
                io.Combo.Input("view", options=list(AUGUR_VIEWS), default="priorities"),
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
        result: AugurResult,
        view: str = "priorities",
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, summary = select_augur_view(result, view=view)
        selected = summary["selected_view"]
        key_results = summary["key_results"]
        parameters = {
            "view": view,
            "artifact_fingerprint_sha256": selected["artifact_fingerprint_sha256"],
            "table_fingerprint_sha256": selected["table_fingerprint_sha256"],
        }
        table_result = make_table_result(
            title=f"Augur {view.replace('_', ' ')}",
            operation="augur_results",
            parameters=parameters,
            description=f"Selected the canonical {view!r} view from a validated Augur result artifact.",
            warnings=[str(warning) for warning in summary["warnings"]],
            input_cells=int(key_results["input_cells"]),
            input_genes=int(key_results["selected_source_features"]),
            started_at=started_at,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title=f"Augur {view.replace('_', ' ')} summary",
            operation="augur_results",
            parameters=parameters,
            description=f"Selected the canonical {view!r} view without rerunning or reinterpreting Augur.",
            warnings=[str(warning) for warning in summary["warnings"]],
            input_cells=int(key_results["input_cells"]),
            input_genes=int(key_results["selected_source_features"]),
            started_at=started_at,
        )
        return io.NodeOutput(table_result, report, augur_results_code(view=view))


class OpenBioSingleCellCellTypeCorrelation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellTypeCorrelation",
            display_name="Population Centroid Correlation",
            category=VISUALIZATION_CATEGORY,
            description=(
                "Compute descriptive correlations between population centroids in one explicit representation, "
                "with a canonical pair table and dendrogram-ordered matrix plot."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("population_key", default="cell_type"),
                io.String.Input("representation_key", default="X_pca"),
                io.Combo.Input(
                    "correlation_method", options=["pearson", "spearman", "kendall"], default="pearson"
                ),
                io.Combo.Input(
                    "linkage_method",
                    options=["complete", "average", "single", "weighted"],
                    default="complete",
                    advanced=True,
                ),
                io.Int.Input("n_dimensions", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Combo.Input(
                    "annotation_status", options=["unknown", "provisional", "curated"], default="unknown"
                ),
                io.String.Input("color_map", default="RdYlBu", advanced=True),
                io.Boolean.Input("show_numbers", default=False, advanced=True),
                io.Int.Input("max_groups", default=200, min=2, max=1000, advanced=True),
                io.Int.Input("max_output_rows", default=100000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                PlotResultType.Output(display_name="plot"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        population_key: str = "cell_type",
        representation_key: str = "X_pca",
        correlation_method: str = "pearson",
        linkage_method: str = "complete",
        n_dimensions: int = 0,
        annotation_status: str = "unknown",
        color_map: str = "RdYlBu",
        show_numbers: bool = False,
        max_groups: int = 200,
        max_output_rows: int = 100000,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, png, summary = analyze_population_centroid_correlation(
            adata,
            population_key=population_key,
            representation_key=representation_key,
            correlation_method=correlation_method,
            linkage_method=linkage_method,
            n_dimensions=n_dimensions,
            annotation_status=annotation_status,
            color_map=color_map,
            show_numbers=show_numbers,
            max_groups=max_groups,
            max_output_rows=max_output_rows,
        )
        parameters = dict(summary["parameters"])
        description = str(summary["results"])
        warnings = [str(warning) for warning in summary["warnings"]]
        input_cells, input_genes = int(adata.n_obs), int(adata.n_vars)
        table_result = make_table_result(
            title=f"Population centroid correlations by {population_key}",
            operation="population_centroid_correlation",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=input_genes,
            started_at=started_at,
            table=table,
        )
        plot_result = make_plot_result(
            title=f"Population centroid correlation by {population_key}",
            operation="population_centroid_correlation",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=input_genes,
            started_at=started_at,
            png=png,
        )
        code = population_correlation_code(
            population_key=population_key,
            representation_key=representation_key,
            correlation_method=correlation_method,
            linkage_method=linkage_method,
            n_dimensions=n_dimensions,
            annotation_status=annotation_status,
            color_map=color_map,
            show_numbers=show_numbers,
            max_groups=max_groups,
            max_output_rows=max_output_rows,
        )
        report = make_summary_result(
            summary=summary,
            title="Population centroid correlation summary",
            operation="population_centroid_correlation",
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=input_genes,
            started_at=started_at,
        )
        return io.NodeOutput(table_result, plot_result, report, code)


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
