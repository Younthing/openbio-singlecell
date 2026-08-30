from __future__ import annotations

from typing import Any

from comfy_api.latest import io

from .expression_source import ExpressionSourceSpec
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType, SummaryResultType

CATEGORY = "openbio/single-cell/data"
GTF_EXTENSIONS = (".gtf", ".gtf.gz")


class OpenBioSingleCellSnapshotExpression(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Raw snapshot source",
        default="X",
        include_raw=False,
        layer_input_id="source_layer",
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSnapshotExpression",
            display_name="Snapshot Expression",
            category=CATEGORY,
            description=(
                "Snapshot an explicit expression source into conventional counts and Raw storage with advisories."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellSubsetObservations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSubsetObservations",
            display_name="Subset Observations",
            category=CATEGORY,
            description="Subset observations by string values in an obs column.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("column", default="sample"),
                io.String.Input("values", default=""),
                io.Boolean.Input("invert", default=False),
                io.Combo.Input(
                    "missing_policy",
                    options=["exclude", "include", "error"],
                    default="exclude",
                    advanced=True,
                ),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellMergeObservationAnnotations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMergeObservationAnnotations",
            display_name="Merge Observation Annotations",
            category=CATEGORY,
            description="Merge an observation annotation from a subset AnnData by obs_names.",
            inputs=[
                AnnDataType.Input("adata"),
                AnnDataType.Input("subset_adata"),
                io.String.Input("source_column", default="cell_type"),
                io.String.Input("target_column", default="cell_type", advanced=True),
                io.Combo.Input(
                    "conflict_policy",
                    options=["error", "keep_target", "overwrite"],
                    default="error",
                    advanced=True,
                ),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellMapGeneIdsFromGTF(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMapGeneIdsFromGTF",
            display_name="Annotate Gene Identity from GTF",
            category=CATEGORY,
            description=(
                "Promote declared stable gene IDs and annotate GTF gene names without filtering expression features."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gtf_path", default="annotations/genes.gtf.gz"),
                io.Combo.Input(
                    "id_source",
                    options=["var_column", "var_names"],
                    default="var_column",
                ),
                io.String.Input("id_column", default="gene_ids", advanced=True),
                io.String.Input("gene_name_column", default="gene_symbols", advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, gtf_path: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gtf_path, extensions=GTF_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gtf_path: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(gtf_path, GTF_EXTENSIONS)


DATA_NODE_CLASSES = [
    OpenBioSingleCellSnapshotExpression,
    OpenBioSingleCellSubsetObservations,
    OpenBioSingleCellMergeObservationAnnotations,
    OpenBioSingleCellMapGeneIdsFromGTF,
]


__all__ = [
    "DATA_NODE_CLASSES",
    "OpenBioSingleCellMapGeneIdsFromGTF",
    "OpenBioSingleCellMergeObservationAnnotations",
    "OpenBioSingleCellSnapshotExpression",
    "OpenBioSingleCellSubsetObservations",
]
