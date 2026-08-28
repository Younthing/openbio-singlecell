from __future__ import annotations

import hashlib
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import PLUGIN_VERSION
from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .cassiopeia_tree import (
    ALLELE_TABLE_EXTENSIONS,
    CassiopeiaCharacters,
    CassiopeiaTree,
    add_cassiopeia_plasticity,
    compute_cassiopeia_expansions,
    prepare_cassiopeia_characters,
    reconstruct_cassiopeia_tree,
)
from .files import resolve_input_path
from .node_types import (
    AnnDataType,
    CassiopeiaCharactersType,
    CassiopeiaTreeType,
    SummaryResultType,
    TableResultType,
)

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/lineage"
MAX_INTEGER = 2**31 - 1


def _summary_result(
    summary: dict[str, Any],
    *,
    title: str,
    operation: str,
    started_at: float,
    input_cells: int,
    input_genes: int = 0,
) -> Any:
    return make_summary_result(
        summary=summary,
        title=title,
        operation=operation,
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=int(input_cells),
        input_genes=int(input_genes),
        started_at=started_at,
    )


def _table_result(
    table: Any,
    summary: dict[str, Any],
    *,
    title: str,
    operation: str,
    started_at: float,
    input_cells: int,
    input_genes: int = 0,
) -> Any:
    return make_table_result(
        table=table,
        title=title,
        operation=operation,
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=int(input_cells),
        input_genes=int(input_genes),
        started_at=started_at,
    )


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
                io.Float.Input(
                    "allele_representation_threshold", default=0.98, min=1e-12, max=1.0, step=0.01
                ),
                io.Int.Input("minimum_cells", default=2, min=1, max=MAX_INTEGER, advanced=True),
                io.Float.Input("maximum_missing_fraction", default=0.8, min=0.0, max=1.0, step=0.01),
                io.Float.Input("maximum_uncut_fraction", default=0.8, min=0.0, max=1.0, step=0.01),
                io.Float.Input("minimum_unique_fraction", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Float.Input(
                    "minimum_informative_character_fraction", default=0.2, min=0.0, max=1.0, step=0.01
                ),
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
    def fingerprint_inputs(
        cls, allele_table_file: str, max_file_mib: int = 512, **kwargs: Any
    ) -> tuple[str, int]:
        del kwargs
        return _bounded_file_fingerprint(allele_table_file, max_file_mib)

    @classmethod
    def execute(
        cls,
        allele_table_file: str,
        first_column_as_index: bool = True,
        lineage_column: str = "Tumor",
        cell_barcode_column: str = "cellBC",
        integration_barcode_column: str = "intBC",
        cut_site_columns: str = "r1,r2,r3",
        prior_grouping_columns: str = "Tumor,intBC",
        missing_data_allele: str = "",
        allele_representation_threshold: float = 0.98,
        minimum_cells: int = 2,
        maximum_missing_fraction: float = 0.8,
        maximum_uncut_fraction: float = 0.8,
        minimum_unique_fraction: float = 0.05,
        minimum_informative_character_fraction: float = 0.2,
        max_file_mib: int = 512,
        max_matrix_gib: float = 2.0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        path = resolve_input_path(allele_table_file, extensions=ALLELE_TABLE_EXTENSIONS)
        artifact, table, summary, code = prepare_cassiopeia_characters(
            path,
            first_column_as_index=first_column_as_index,
            lineage_column=lineage_column,
            cell_barcode_column=cell_barcode_column,
            integration_barcode_column=integration_barcode_column,
            cut_site_columns=cut_site_columns,
            prior_grouping_columns=prior_grouping_columns,
            missing_data_allele=missing_data_allele,
            allele_representation_threshold=allele_representation_threshold,
            minimum_cells=minimum_cells,
            maximum_missing_fraction=maximum_missing_fraction,
            maximum_uncut_fraction=maximum_uncut_fraction,
            minimum_unique_fraction=minimum_unique_fraction,
            minimum_informative_character_fraction=minimum_informative_character_fraction,
            max_file_mib=max_file_mib,
            max_matrix_gib=max_matrix_gib,
            openbio_version=PLUGIN_VERSION,
        )
        input_cells = int(summary["key_results"]["globally_unique_cells"])
        table_result = _table_result(
            table,
            summary,
            title="Cassiopeia character preparation QC",
            operation="prepare_cassiopeia_characters",
            started_at=started_at,
            input_cells=input_cells,
        )
        report = _summary_result(
            summary,
            title="Cassiopeia character preparation summary",
            operation="prepare_cassiopeia_characters",
            started_at=started_at,
            input_cells=input_cells,
        )
        return io.NodeOutput(artifact, table_result, report, code)


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

    @classmethod
    def execute(
        cls,
        characters: CassiopeiaCharacters,
        lineage_id: str,
        prior_transformation: str = "negative_log",
        collapse_mutationless_edges: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        artifact, summary, code = reconstruct_cassiopeia_tree(
            characters,
            lineage_id,
            prior_transformation=prior_transformation,
            collapse_mutationless_edges=collapse_mutationless_edges,
            openbio_version=PLUGIN_VERSION,
        )
        report = _summary_result(
            summary,
            title=f"Cassiopeia VanillaGreedy tree: {artifact.lineage_id}",
            operation="reconstruct_cassiopeia_vanilla_greedy",
            started_at=started_at,
            input_cells=artifact.input_cells,
        )
        return io.NodeOutput(artifact, report, code)


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

    @classmethod
    def execute(
        cls,
        tree: CassiopeiaTree,
        minimum_clade_size: int = 10,
        minimum_depth: int = 1,
        fdr_threshold: float = 0.05,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, summary, code = compute_cassiopeia_expansions(
            tree,
            minimum_clade_size=minimum_clade_size,
            minimum_depth=minimum_depth,
            fdr_threshold=fdr_threshold,
            openbio_version=PLUGIN_VERSION,
        )
        table_result = _table_result(
            table,
            summary,
            title=f"Cassiopeia clade expansion evidence: {tree.lineage_id}",
            operation="cassiopeia_clade_expansion_test",
            started_at=started_at,
            input_cells=tree.input_cells,
        )
        report = _summary_result(
            summary,
            title=f"Cassiopeia clade expansion summary: {tree.lineage_id}",
            operation="cassiopeia_clade_expansion_test",
            started_at=started_at,
            input_cells=tree.input_cells,
        )
        return io.NodeOutput(table_result, report, code)


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
                io.Combo.Input(
                    "analysis_mode", options=["exploratory", "report_grade"], default="exploratory"
                ),
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

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        tree: CassiopeiaTree,
        annotation_key: str = "cell_type",
        annotation_status: str = "unknown",
        analysis_mode: str = "exploratory",
        minimum_state_fraction: float = 0.025,
        output_key: str = "sc_effective_plasticity",
        overwrite_existing: bool = False,
        max_working_gib: float = 4.0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        output, table, summary, code = add_cassiopeia_plasticity(
            adata,
            tree,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            analysis_mode=analysis_mode,
            minimum_state_fraction=minimum_state_fraction,
            output_key=output_key,
            overwrite_existing=overwrite_existing,
            max_working_gib=max_working_gib,
            openbio_version=PLUGIN_VERSION,
        )
        finish_adata(
            output,
            "cassiopeia_effective_plasticity",
            summary["parameters"],
            int(adata.n_obs),
            int(adata.n_vars),
            started_at,
        )
        table_result = _table_result(
            table,
            summary,
            title=f"Cassiopeia EffectivePlasticity cells: {tree.lineage_id}",
            operation="cassiopeia_effective_plasticity",
            started_at=started_at,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
        )
        report = _summary_result(
            summary,
            title=f"Cassiopeia EffectivePlasticity summary: {tree.lineage_id}",
            operation="cassiopeia_effective_plasticity",
            started_at=started_at,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
        )
        return io.NodeOutput(output, table_result, report, code)


LINEAGE_NODE_CLASSES = [
    OpenBioSingleCellCassiopeiaLineageQC,
    OpenBioSingleCellReconstructCassiopeiaTree,
    OpenBioSingleCellCassiopeiaExpansionTest,
    OpenBioSingleCellCassiopeiaPlasticity,
]


__all__ = ["LINEAGE_NODE_CLASSES"]
