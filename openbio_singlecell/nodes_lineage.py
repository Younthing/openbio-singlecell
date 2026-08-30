from __future__ import annotations

import hashlib
from typing import Any

from comfy_api.latest import io

from .cassiopeia_tree import ALLELE_TABLE_EXTENSIONS
from .files import resolve_input_path
from .node_types import (
    AnnDataType,
    CassiopeiaCharactersType,
    CassiopeiaTreeType,
    SummaryResultType,
    TableResultType,
)

CATEGORY = "openbio/single-cell/lineage"
MAX_INTEGER = 2**31 - 1


def _bounded_file_fingerprint(relative_path: str, max_file_mib: int) -> tuple[str, int]:
    path = resolve_input_path(relative_path, extensions=ALLELE_TABLE_EXTENSIONS)
    limit = int(max_file_mib) * 1024 * 1024
    if int(max_file_mib) <= 0:
        raise ValueError("max_file_mib must be positive.")
    with open(path, "rb") as handle:
        content = handle.read(limit + 1)
    if len(content) > limit:
        raise ValueError(f"Allele table exceeds max_file_mib={int(max_file_mib)}.")
    return hashlib.sha256(content).hexdigest(), len(content)


class OpenBioSingleCellCassiopeiaLineageQC(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCassiopeiaLineageQC",
            display_name="Prepare Cassiopeia Characters + QC",
            category=CATEGORY,
            description=(
                "Conflict-audit one bounded allele table, compute empirical priors, encode exact cut sites, and emit "
                "a fingerprinted character artifact with disclosed QC pass/warning status."
            ),
            inputs=[
                io.String.Input("allele_table_file", default="openbio-singlecell/allele_table.tsv"),
                io.Boolean.Input("first_column_as_index", default=True, advanced=True),
                io.String.Input("lineage_column", default="Tumor", advanced=True),
                io.String.Input("cell_barcode_column", default="cellBC", advanced=True),
                io.String.Input("integration_barcode_column", default="intBC", advanced=True),
                io.String.Input("cut_site_columns", default="r1,r2,r3"),
                io.String.Input("prior_grouping_columns", default="Tumor,intBC"),
                io.String.Input("missing_data_allele", default="", advanced=True),
                io.Float.Input("allele_representation_threshold", default=0.98, min=1e-12, max=1.0, step=0.01),
                io.Int.Input("minimum_cells", default=2, min=1, max=MAX_INTEGER, advanced=True),
                io.Float.Input("maximum_missing_fraction", default=0.8, min=0.0, max=1.0, step=0.01),
                io.Float.Input("maximum_uncut_fraction", default=0.8, min=0.0, max=1.0, step=0.01),
                io.Float.Input("minimum_unique_fraction", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Float.Input("minimum_informative_character_fraction", default=0.2, min=0.0, max=1.0, step=0.01),
                io.Int.Input("max_file_mib", default=512, min=1, max=MAX_INTEGER, advanced=True),
                io.Float.Input("max_matrix_gib", default=2.0, min=1e-12, advanced=True),
            ],
            outputs=[
                CassiopeiaCharactersType.Output(display_name="characters"),
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, allele_table_file: str, max_file_mib: int = 512, **kwargs: Any) -> bool | str:
        del kwargs
        try:
            _bounded_file_fingerprint(allele_table_file, max_file_mib)
        except (ValueError, FileNotFoundError, OSError) as error:
            return str(error)
        return True

    @classmethod
    def fingerprint_inputs(cls, allele_table_file: str, max_file_mib: int = 512, **kwargs: Any) -> tuple[str, int]:
        del kwargs
        return _bounded_file_fingerprint(allele_table_file, max_file_mib)


class OpenBioSingleCellReconstructCassiopeiaTree(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellReconstructCassiopeiaTree",
            display_name="Reconstruct Cassiopeia VanillaGreedy Tree",
            category=CATEGORY,
            description="Solve one computable lineage character payload with fixed Cassiopeia VanillaGreedy.",
            inputs=[
                CassiopeiaCharactersType.Input("characters"),
                io.String.Input("lineage_id", default=""),
                io.Combo.Input(
                    "prior_transformation",
                    options=["negative_log", "inverse", "square_root_inverse"],
                    default="negative_log",
                ),
                io.Boolean.Input("collapse_mutationless_edges", default=False, advanced=True),
            ],
            outputs=[
                CassiopeiaTreeType.Output(display_name="tree"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellCassiopeiaExpansionTest(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCassiopeiaExpansionTest",
            display_name="Cassiopeia Clade Expansion Test",
            category=CATEGORY,
            description="Test one complete within-tree clade family and apply Benjamini-Hochberg FDR correction.",
            inputs=[
                CassiopeiaTreeType.Input("tree"),
                io.Int.Input("minimum_clade_size", default=10, min=1, max=MAX_INTEGER),
                io.Int.Input("minimum_depth", default=1, min=0, max=MAX_INTEGER, advanced=True),
                io.Float.Input("fdr_threshold", default=0.05, min=0.0, max=1.0, step=0.01),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


class OpenBioSingleCellCassiopeiaPlasticity(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCassiopeiaPlasticity",
            display_name="Cassiopeia EffectivePlasticity",
            category=CATEGORY,
            description="Compute descriptive per-cell EffectivePlasticity for one tree and categorical annotation.",
            inputs=[
                AnnDataType.Input("adata"),
                CassiopeiaTreeType.Input("tree"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                ),
                io.Combo.Input("analysis_mode", options=["exploratory", "report_grade"], default="exploratory"),
                io.Float.Input("minimum_state_fraction", default=0.025, min=0.0, max=1.0, step=0.005),
                io.String.Input("output_key", default="sc_effective_plasticity", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_working_gib", default=4.0, min=1e-12, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


LINEAGE_NODE_CLASSES = [
    OpenBioSingleCellCassiopeiaLineageQC,
    OpenBioSingleCellReconstructCassiopeiaTree,
    OpenBioSingleCellCassiopeiaExpansionTest,
    OpenBioSingleCellCassiopeiaPlasticity,
]

__all__ = ["LINEAGE_NODE_CLASSES"]
