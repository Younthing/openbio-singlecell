from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_table_result
from .cassiopeia_tree import (
    ALLELE_TABLE_EXTENSIONS,
    CassiopeiaTree,
    _normalize_columns,
    _read_allele_table,
    _require_cassiopeia,
    add_cassiopeia_plasticity,
    compute_cassiopeia_expansions,
    reconstruct_cassiopeia_tree,
)
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType, CassiopeiaTreeType, TableResultType

if TYPE_CHECKING:
    from anndata import AnnData

CATEGORY = "openbio/single-cell/lineage"
MAX_INTEGER = 2**31 - 1


def _percent_uncut(values: Any, science: dependencies.ScientificDependencies) -> float:
    values = science.np.asarray(values)
    observed = values != -1
    return float(science.np.count_nonzero(values == 0) / max(1, int(science.np.count_nonzero(observed))))


def _percent_indels(character_matrix: Any, science: dependencies.ScientificDependencies) -> float:
    values = science.np.asarray(character_matrix).ravel()
    observed = values != -1
    observed_count = int(science.np.count_nonzero(observed))
    if observed_count == 0:
        return float("nan")
    return float(1.0 - (science.np.count_nonzero(values[observed] == 0) / observed_count))


class OpenBioSingleCellCassiopeiaLineageQC(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCassiopeiaLineageQC",
            display_name="Cassiopeia Lineage QC",
            category=CATEGORY,
            description="Summarize lineage-tracing quality for each tumor in a Cassiopeia allele table.",
            inputs=[
                io.String.Input(
                    "allele_table_file",
                    default="openbio-singlecell/allele_table.tsv",
                ),
                io.Boolean.Input("first_column_as_index", default=True, advanced=True),
                io.String.Input("tumor_column", default="Tumor", advanced=True),
                io.String.Input("cell_barcode_column", default="cellBC", advanced=True),
                io.String.Input("integration_barcode_column", default="intBC", advanced=True),
                io.Int.Input("cut_sites_per_intbc", default=3, min=1, max=MAX_INTEGER, advanced=True),
                io.Float.Input(
                    "minimum_intbc_fraction",
                    default=0.2,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Int.Input(
                    "minimum_cells_for_summary",
                    default=2,
                    min=1,
                    max=MAX_INTEGER,
                    advanced=True,
                ),
                io.Float.Input(
                    "maximum_uncut_fraction",
                    default=0.8,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Float.Input(
                    "allele_representation_threshold",
                    default=0.98,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Int.Input(
                    "lineage_size_threshold",
                    default=100,
                    min=1,
                    max=MAX_INTEGER,
                    advanced=True,
                ),
                io.Float.Input(
                    "percent_unique_threshold",
                    default=0.05,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                    advanced=True,
                ),
                io.Float.Input(
                    "percent_unsaturated_threshold",
                    default=0.2,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                    advanced=True,
                ),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def validate_inputs(cls, allele_table_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(allele_table_file, extensions=ALLELE_TABLE_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as error:
            return str(error)
        return True

    @classmethod
    def fingerprint_inputs(cls, allele_table_file: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(allele_table_file, ALLELE_TABLE_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        allele_table_file: str,
        first_column_as_index: bool = True,
        tumor_column: str = "Tumor",
        cell_barcode_column: str = "cellBC",
        integration_barcode_column: str = "intBC",
        cut_sites_per_intbc: int = 3,
        minimum_intbc_fraction: float = 0.2,
        minimum_cells_for_summary: int = 2,
        maximum_uncut_fraction: float = 0.8,
        allele_representation_threshold: float = 0.98,
        lineage_size_threshold: int = 100,
        percent_unique_threshold: float = 0.05,
        percent_unsaturated_threshold: float = 0.2,
    ) -> io.NodeOutput:
        cassiopeia, lineage_utils = _require_cassiopeia()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        allele_table = _read_allele_table(allele_table_file, first_column_as_index)
        allele_table = _normalize_columns(
            allele_table,
            {
                "Tumor": tumor_column,
                "cellBC": cell_barcode_column,
                "intBC": integration_barcode_column,
            },
        )
        rows = []
        warnings = []
        skipped_small_tumors = 0
        for tumor, tumor_allele_table in allele_table.groupby("Tumor", observed=True, sort=True):
            if int(tumor_allele_table["cellBC"].nunique()) < minimum_cells_for_summary:
                skipped_small_tumors += 1
                continue

            tumor_allele_table = tumor_allele_table.copy()
            tumor_allele_table["lineageGrp"] = tumor_allele_table["Tumor"]
            lineage_group = lineage_utils.filter_intbcs_final_lineages(
                tumor_allele_table,
                min_intbc_thresh=minimum_intbc_fraction,
            )[0]
            number_of_cut_sites = int(lineage_group["intBC"].nunique()) * cut_sites_per_intbc
            if number_of_cut_sites == 0:
                warnings.append(f"Skipped tumor {tumor!r} because no integration barcodes passed filtering.")
                continue

            character_matrix, _, _ = cassiopeia.pp.convert_alleletable_to_character_matrix(
                lineage_group,
                allele_rep_thresh=allele_representation_threshold,
            )
            if character_matrix.shape[1] == 0:
                character_matrix, _, _ = cassiopeia.pp.convert_alleletable_to_character_matrix(
                    lineage_group,
                    allele_rep_thresh=1.0,
                )
            if character_matrix.shape[1] == 0:
                warnings.append(f"Skipped tumor {tumor!r} because no lineage characters remained.")
                continue

            percent_uncut = character_matrix.apply(
                lambda row: _percent_uncut(row.to_numpy(), science),
                axis=1,
            )
            filtered = character_matrix[percent_uncut < maximum_uncut_fraction]
            if filtered.empty:
                percent_unique = float("nan")
                cut_rate = float("nan")
                warnings.append(f"Tumor {tumor!r} has no cells below the maximum uncut fraction.")
            else:
                percent_unique = float(filtered.drop_duplicates().shape[0] / filtered.shape[0])
                cut_rate = _percent_indels(filtered, science)

            saturated_targets = number_of_cut_sites - int(character_matrix.shape[1])
            percent_unsaturated = float(1.0 - (saturated_targets / number_of_cut_sites))
            rows.append(
                {
                    "Tumor": str(tumor),
                    "PercentUnique": percent_unique,
                    "CutRate": cut_rate,
                    "NumSaturatedTargets": saturated_targets,
                    "PercentUnsaturatedTargets": percent_unsaturated,
                    "NumCells": int(filtered.shape[0]),
                }
            )

        if not rows:
            raise ValueError("No tumors produced Cassiopeia lineage QC statistics.")
        table = science.pd.DataFrame.from_records(rows)
        table["PoorQC"] = (
            table["PercentUnique"].isna()
            | (table["PercentUnique"] <= percent_unique_threshold)
            | (table["PercentUnsaturatedTargets"] <= percent_unsaturated_threshold)
        )
        table["SmallLineage"] = table["NumCells"] < lineage_size_threshold
        table["PassesQC"] = ~(table["PoorQC"] | table["SmallLineage"])
        if skipped_small_tumors:
            warnings.append(f"Skipped {skipped_small_tumors} tumors with fewer than {minimum_cells_for_summary} cells.")

        parameters = {
            "allele_table_file": allele_table_file,
            "first_column_as_index": first_column_as_index,
            "tumor_column": tumor_column,
            "cell_barcode_column": cell_barcode_column,
            "integration_barcode_column": integration_barcode_column,
            "cut_sites_per_intbc": cut_sites_per_intbc,
            "minimum_intbc_fraction": minimum_intbc_fraction,
            "minimum_cells_for_summary": minimum_cells_for_summary,
            "maximum_uncut_fraction": maximum_uncut_fraction,
            "allele_representation_threshold": allele_representation_threshold,
            "lineage_size_threshold": lineage_size_threshold,
            "percent_unique_threshold": percent_unique_threshold,
            "percent_unsaturated_threshold": percent_unsaturated_threshold,
        }
        result = make_table_result(
            title="Cassiopeia lineage quality",
            operation="cassiopeia_lineage_qc",
            parameters=parameters,
            description="Tumor-level lineage quality statistics from the Cassiopeia allele-table workflow.",
            warnings=warnings,
            input_cells=int(allele_table["cellBC"].nunique()),
            input_genes=0,
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellReconstructCassiopeiaTree(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellReconstructCassiopeiaTree",
            display_name="Reconstruct Cassiopeia Tree",
            category=CATEGORY,
            description="Reconstruct and solve one reusable Cassiopeia tumor lineage tree.",
            inputs=[
                io.String.Input(
                    "allele_table_file",
                    default="openbio-singlecell/allele_table.tsv",
                ),
                io.String.Input("tumor", default=""),
                io.Boolean.Input("first_column_as_index", default=True, advanced=True),
                io.String.Input("tumor_column", default="Tumor", advanced=True),
                io.String.Input("cell_barcode_column", default="cellBC", advanced=True),
                io.String.Input("integration_barcode_column", default="intBC", advanced=True),
                io.String.Input("mutation_family_column", default="MetFamily", advanced=True),
                io.Float.Input(
                    "allele_representation_threshold",
                    default=0.9,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
            ],
            outputs=[CassiopeiaTreeType.Output(display_name="tree")],
        )

    @classmethod
    def validate_inputs(cls, allele_table_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(allele_table_file, extensions=ALLELE_TABLE_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as error:
            return str(error)
        return True

    @classmethod
    def fingerprint_inputs(cls, allele_table_file: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(allele_table_file, ALLELE_TABLE_EXTENSIONS)

    @classmethod
    def execute(
        cls,
        allele_table_file: str,
        tumor: str = "",
        first_column_as_index: bool = True,
        tumor_column: str = "Tumor",
        cell_barcode_column: str = "cellBC",
        integration_barcode_column: str = "intBC",
        mutation_family_column: str = "MetFamily",
        allele_representation_threshold: float = 0.9,
    ) -> io.NodeOutput:
        tree = reconstruct_cassiopeia_tree(
            allele_table_file=allele_table_file,
            tumor=tumor,
            first_column_as_index=first_column_as_index,
            tumor_column=tumor_column,
            cell_barcode_column=cell_barcode_column,
            integration_barcode_column=integration_barcode_column,
            mutation_family_column=mutation_family_column,
            allele_representation_threshold=allele_representation_threshold,
        )
        return io.NodeOutput(tree)


class OpenBioSingleCellCassiopeiaExpansionTest(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCassiopeiaExpansionTest",
            display_name="Cassiopeia Expansion Test",
            category=CATEGORY,
            description="Test the clades of a reconstructed Cassiopeia tree for expansion.",
            inputs=[
                CassiopeiaTreeType.Input("tree"),
                io.Float.Input(
                    "minimum_clade_fraction",
                    default=0.15,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Int.Input("minimum_depth", default=1, min=0, max=MAX_INTEGER, advanced=True),
                io.Float.Input(
                    "expansion_pvalue_threshold",
                    default=0.01,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
            ],
            outputs=[TableResultType.Output(display_name="table")],
        )

    @classmethod
    def execute(
        cls,
        tree: CassiopeiaTree,
        minimum_clade_fraction: float = 0.15,
        minimum_depth: int = 1,
        expansion_pvalue_threshold: float = 0.01,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, effective_min_clade_size = compute_cassiopeia_expansions(
            tree,
            minimum_clade_fraction=minimum_clade_fraction,
            minimum_depth=minimum_depth,
            expansion_pvalue_threshold=expansion_pvalue_threshold,
        )
        warnings = []
        if not bool(table["is_expansion"].any()):
            warnings.append("No lineage nodes passed the requested expansion p-value threshold.")

        parameters = {
            "tree_provenance": dict(tree.provenance),
            "tumor": tree.tumor,
            "minimum_clade_fraction": minimum_clade_fraction,
            "effective_min_clade_size": effective_min_clade_size,
            "minimum_depth": minimum_depth,
            "expansion_pvalue_threshold": expansion_pvalue_threshold,
            "tree_cells": tree.input_cells,
            "tree_characters": tree.character_count,
        }
        result = make_table_result(
            title=f"Cassiopeia expansions: {tree.tumor}",
            operation="cassiopeia_expansion_test",
            parameters=parameters,
            description="Node-level expansion probabilities from a VanillaGreedy Cassiopeia lineage reconstruction.",
            warnings=warnings,
            input_cells=tree.input_cells,
            input_genes=0,
            started_at=started_at,
            table=table,
        )
        return io.NodeOutput(result)


class OpenBioSingleCellCassiopeiaPlasticity(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCassiopeiaPlasticity",
            display_name="Cassiopeia Plasticity",
            category=CATEGORY,
            description="Add single-cell effective plasticity from a reconstructed Cassiopeia tree to AnnData.",
            inputs=[
                AnnDataType.Input("adata"),
                CassiopeiaTreeType.Input("tree"),
                io.String.Input("annotation_key", default="cell_type"),
                io.String.Input("output_key", default="scPlasticity", advanced=True),
                io.String.Input("summary_key", default="cassiopeia_plasticity", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        tree: CassiopeiaTree,
        annotation_key: str = "cell_type",
        output_key: str = "scPlasticity",
        summary_key: str = "cassiopeia_plasticity",
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        output = add_cassiopeia_plasticity(
            tree,
            adata,
            annotation_key=annotation_key,
            output_key=output_key,
            summary_key=summary_key,
        )
        parameters = {
            "tree_provenance": dict(tree.provenance),
            "tumor": tree.tumor,
            "annotation_key": annotation_key.strip(),
            "output_key": output_key.strip(),
            "summary_key": summary_key.strip(),
            "tree_cells": tree.input_cells,
            "tree_characters": tree.character_count,
        }
        finish_adata(
            output,
            "cassiopeia_plasticity",
            parameters,
            int(adata.n_obs),
            int(adata.n_vars),
            started_at,
        )
        return io.NodeOutput(output)


LINEAGE_NODE_CLASSES = [
    OpenBioSingleCellCassiopeiaLineageQC,
    OpenBioSingleCellReconstructCassiopeiaTree,
    OpenBioSingleCellCassiopeiaExpansionTest,
    OpenBioSingleCellCassiopeiaPlasticity,
]


__all__ = ["LINEAGE_NODE_CLASSES"]
