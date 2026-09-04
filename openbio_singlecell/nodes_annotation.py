from __future__ import annotations

import hashlib
from typing import Any

from comfy_api.latest import io

from .annotation_core import celltypist_model_fingerprint
from .expression_source import _CELLTYPIST_SPEC
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType, PlotResultType, TableResultType, analysis_outputs

CATEGORY = "openbio/single-cell/annotation"


class OpenBioSingleCellCellTypistAnnotation(io.ComfyNode):
    EXPRESSION_SOURCE = _CELLTYPIST_SPEC

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellTypistAnnotation",
            display_name="CellTypist Annotation",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Combo.Input(
                    "expression_state",
                    options=["verified_cp10k_log1p", "verified_counts"],
                    default="verified_cp10k_log1p",
                ),
                io.String.Input("model", default=""),
                io.Boolean.Input("majority_voting", default=False),
                io.String.Input("over_clustering_key", default="", advanced=True),
                io.Float.Input("min_prop", default=0.0, min=0.0, max=1.0, step=0.05, advanced=True),
                io.String.Input("label_column", default="celltypist_cell_type", advanced=True),
                io.String.Input("confidence_column", default="celltypist_confidence", advanced=True),
                io.String.Input("probability_key", default="celltypist_probabilities", advanced=True),
                io.Boolean.Input("store_decision_matrix", default=False, advanced=True),
                io.String.Input("decision_key", default="celltypist_decision_scores", advanced=True),
                io.String.Input("metadata_key", default="celltypist", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )

    @classmethod
    def fingerprint_inputs(cls, model: str = "", **kwargs: Any) -> Any:
        return celltypist_model_fingerprint(model)


class OpenBioSingleCellMarkerORAEvidence(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerORAEvidence",
            display_name="Marker ORA Evidence",
            category=CATEGORY,
            inputs=[
                TableResultType.Input("table"),
                TableResultType.Input("universe"),
                io.String.Input("resource_csv", default=""),
                io.String.Input("resource_metadata_json", default="{}"),
                io.String.Input("source_column", default="source", advanced=True),
                io.String.Input("target_column", default="target", advanced=True),
                io.Int.Input("min_targets", default=3, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("min_overlap", default=2, min=1, max=2**31 - 1, advanced=True),
                io.Float.Input("max_p_adjusted", default=0.05, min=0.0, max=1.0, step=0.01),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )

    @classmethod
    def validate_inputs(cls, resource_csv: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(resource_csv, extensions=(".csv",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, resource_csv: str, **kwargs: Any) -> Any:
        path = resolve_input_path(resource_csv, extensions=(".csv",))
        identity = input_file_fingerprint(resource_csv, (".csv",))
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return ("openbio-marker-ora-resource-v1", *identity, digest.hexdigest())


class OpenBioSingleCellMarkerORAEvidencePlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerORAEvidencePlot",
            display_name="Marker ORA Evidence Plot",
            category=CATEGORY,
            description=(
                "Render stored Marker ORA enrichment or overlap evidence without loading a resource, retesting, "
                "or assigning annotation labels."
            ),
            inputs=[
                TableResultType.Input("table"),
                io.DynamicCombo.Input(
                    "view",
                    options=[
                        io.DynamicCombo.Option(
                            "enrichment_dot",
                            [io.Int.Input("max_terms_per_group", default=10, min=1, max=30)],
                        ),
                        io.DynamicCombo.Option(
                            "overlap_bar",
                            [io.Int.Input("max_terms_per_group", default=10, min=1, max=30)],
                        ),
                    ],
                ),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellCellTypistDiagnosticsPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellTypistDiagnosticsPlot",
            display_name="CellTypist Diagnostics Plot",
            category=CATEGORY,
            description=(
                "Render confidence or class-score diagnostics from stored CellTypist Provisional annotation "
                "without rerunning the classifier."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("metadata_key", default="celltypist"),
                io.DynamicCombo.Input(
                    "view",
                    options=[
                        io.DynamicCombo.Option(
                            "confidence_distributions",
                            [io.Int.Input("max_labels", default=30, min=1, max=100)],
                        ),
                        io.DynamicCombo.Option(
                            "probability_heatmap",
                            [
                                io.String.Input("groupby", default="leiden"),
                                io.Int.Input("max_classes", default=30, min=1, max=100),
                            ],
                        ),
                    ],
                ),
            ],
            outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
        )


class OpenBioSingleCellMapClusterAnnotations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMapClusterAnnotations",
            display_name="Map Cluster Annotations",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.String.Input("mapping_json", default="{}"),
                io.String.Input("output_column", default="cell_type", advanced=True),
                io.Combo.Input(
                    "unmapped_policy",
                    options=["error", "preserve_cluster_label", "set_missing"],
                    default="error",
                ),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


ANNOTATION_NODE_CLASSES = [
    OpenBioSingleCellCellTypistAnnotation,
    OpenBioSingleCellCellTypistDiagnosticsPlot,
    OpenBioSingleCellMarkerORAEvidence,
    OpenBioSingleCellMarkerORAEvidencePlot,
    OpenBioSingleCellMapClusterAnnotations,
]


__all__ = [
    "ANNOTATION_NODE_CLASSES",
    "OpenBioSingleCellCellTypistAnnotation",
    "OpenBioSingleCellCellTypistDiagnosticsPlot",
    "OpenBioSingleCellMapClusterAnnotations",
    "OpenBioSingleCellMarkerORAEvidence",
    "OpenBioSingleCellMarkerORAEvidencePlot",
]
