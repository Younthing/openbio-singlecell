from __future__ import annotations

from comfy_api.latest import io

from .augur import AUGUR_CLASSIFIERS, AUGUR_VIEWS
from .expression_source import _AUGUR_SPEC
from .node_types import AnnDataType, AugurResultType, PlotResultType, TableResultType, analysis_outputs

PRIORITY_CATEGORY = "openbio/single-cell/cell-prioritization"
DIAGNOSTIC_CATEGORY = "openbio/single-cell/diagnostics"
AUGUR_EXPRESSION_SOURCE = _AUGUR_SPEC


class OpenBioSingleCellAugur(io.ComfyNode):
    EXPRESSION_SOURCE = AUGUR_EXPRESSION_SOURCE

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAugur",
            display_name="Augur Cell Prioritization",
            category=PRIORITY_CATEGORY,
            description=("Prioritize populations with audited two-Condition Pertpy Augur classifier cross-validation."),
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
            outputs=analysis_outputs(AugurResultType.Output(display_name="result")),
        )


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
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellCellTypeCorrelation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellTypeCorrelation",
            display_name="Population Centroid Correlation",
            category=DIAGNOSTIC_CATEGORY,
            description=(
                "Compute descriptive correlations between population centroids in one explicit representation, "
                "with a canonical pair table and dendrogram-ordered matrix plot."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("population_key", default="cell_type"),
                io.String.Input("representation_key", default="X_pca"),
                io.Combo.Input("correlation_method", options=["pearson", "spearman", "kendall"], default="pearson"),
                io.Combo.Input(
                    "linkage_method",
                    options=["complete", "average", "single", "weighted"],
                    default="complete",
                    advanced=True,
                ),
                io.Int.Input("n_dimensions", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Combo.Input("annotation_status", options=["unknown", "provisional", "curated"], default="unknown"),
                io.String.Input("color_map", default="RdYlBu", advanced=True),
                io.Boolean.Input("show_numbers", default=False, advanced=True),
                io.Int.Input("max_groups", default=200, min=2, max=1000, advanced=True),
                io.Int.Input("max_output_rows", default=100000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(
                TableResultType.Output(display_name="table"), PlotResultType.Output(display_name="plot")
            ),
        )


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
