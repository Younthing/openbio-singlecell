from __future__ import annotations

import importlib
import os
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_utils import finish_adata, make_result
from .files import input_file_fingerprint, resolve_input_path
from .node_types import AnnDataType, SingleCellResultType

if TYPE_CHECKING:
    from anndata import AnnData

CATEGORY = "openbio/single-cell/lineage"
ALLELE_TABLE_EXTENSIONS = (".txt", ".tsv", ".csv")
MAX_INTEGER = 2**31 - 1


def _require_cassiopeia() -> tuple[Any, Any]:
    try:
        cassiopeia = importlib.import_module("cassiopeia")
        lineage_utils = importlib.import_module("cassiopeia.preprocess.lineage_utils")
    except (ImportError, OSError) as error:
        raise RuntimeError("Cassiopeia lineage nodes require the cassiopeia-lineage package.") from error
    return cassiopeia, lineage_utils


def _read_allele_table(path: str, first_column_as_index: bool) -> Any:
    science = dependencies.require_scientific_dependencies()
    resolved = resolve_input_path(path, extensions=ALLELE_TABLE_EXTENSIONS)
    separator = "," if os.path.splitext(resolved)[1].lower() == ".csv" else "\t"
    table = science.pd.read_csv(resolved, sep=separator, index_col=0 if first_column_as_index else None)
    if table.empty:
        raise ValueError("The allele table is empty.")
    return table


def _normalize_cassiopeia_columns(table: Any, columns: dict[str, str]) -> Any:
    normalized = {canonical: source.strip() for canonical, source in columns.items()}
    if any(not source for source in normalized.values()):
        raise ValueError("Allele-table column names cannot be empty.")
    if len(set(normalized.values())) != len(normalized):
        raise ValueError("Allele-table column names must be different.")
    missing = [source for source in normalized.values() if source not in table.columns]
    if missing:
        raise ValueError(f"Allele table is missing columns: {missing}")

    source_columns = set(normalized.values())
    collisions = [
        canonical
        for canonical, source in normalized.items()
        if source != canonical and canonical in table.columns and canonical not in source_columns
    ]
    if collisions:
        raise ValueError(f"Allele-table column mapping conflicts with existing columns: {collisions}")
    return table.rename(columns={source: canonical for canonical, source in normalized.items()})


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
            outputs=[SingleCellResultType.Output(display_name="result")],
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
        allele_table = _normalize_cassiopeia_columns(
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
        result = make_result(
            kind="table",
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


class OpenBioSingleCellCassiopeiaExpansionTest(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCassiopeiaExpansionTest",
            display_name="Cassiopeia Expansion Test",
            category=CATEGORY,
            description="Reconstruct one tumor lineage and test its clades for expansion.",
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
            outputs=[SingleCellResultType.Output(display_name="result")],
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
        minimum_clade_fraction: float = 0.15,
        minimum_depth: int = 1,
        expansion_pvalue_threshold: float = 0.01,
    ) -> io.NodeOutput:
        tumor = tumor.strip()
        if not tumor:
            raise ValueError("Tumor cannot be empty.")

        cassiopeia, _ = _require_cassiopeia()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        allele_table = _read_allele_table(allele_table_file, first_column_as_index)
        allele_table = _normalize_cassiopeia_columns(
            allele_table,
            {
                "Tumor": tumor_column,
                "cellBC": cell_barcode_column,
                "intBC": integration_barcode_column,
                "MetFamily": mutation_family_column,
            },
        )
        tumor_allele_table = allele_table[allele_table["Tumor"].astype(str) == tumor].copy()
        if tumor_allele_table.empty:
            raise ValueError(f"Tumor not found in allele table: {tumor!r}")
        input_cells = int(tumor_allele_table["cellBC"].nunique())
        if input_cells < 2:
            raise ValueError("Cassiopeia expansion testing requires at least two cells in the selected tumor.")

        indel_priors = cassiopeia.pp.compute_empirical_indel_priors(
            allele_table,
            grouping_variables=["intBC", "MetFamily"],
        )
        character_matrix, priors, _ = cassiopeia.pp.convert_alleletable_to_character_matrix(
            tumor_allele_table,
            allele_rep_thresh=allele_representation_threshold,
            mutation_priors=indel_priors,
        )
        if character_matrix.shape[0] < 2 or character_matrix.shape[1] == 0:
            raise ValueError("The selected tumor did not produce enough cells and characters to reconstruct a tree.")

        tree = cassiopeia.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)
        solver = cassiopeia.solver.VanillaGreedySolver()
        solver.solve(tree)
        effective_min_clade_size = minimum_clade_fraction * tree.n_cell
        cassiopeia.tl.compute_expansion_pvalues(
            tree,
            min_clade_size=effective_min_clade_size,
            min_depth=minimum_depth,
        )

        rows = []
        for node in tree.depth_first_traverse_nodes():
            expansion_pvalue = float(tree.get_attribute(node, "expansion_pvalue"))
            depth = 0
            parent = None
            if node != tree.root:
                parent = tree.parent(node)
                current = parent
                depth = 1
                while current != tree.root:
                    current = tree.parent(current)
                    depth += 1
            rows.append(
                {
                    "tumor": tumor,
                    "node": str(node),
                    "parent": None if parent is None else str(parent),
                    "depth": depth,
                    "leaf_count": len(tree.leaves_in_subtree(node)),
                    "expansion_pvalue": expansion_pvalue,
                    "is_expansion": expansion_pvalue < expansion_pvalue_threshold,
                }
            )
        table = science.pd.DataFrame.from_records(rows)
        warnings = []
        if not bool(table["is_expansion"].any()):
            warnings.append("No lineage nodes passed the requested expansion p-value threshold.")

        parameters = {
            "allele_table_file": allele_table_file,
            "tumor": tumor,
            "first_column_as_index": first_column_as_index,
            "tumor_column": tumor_column,
            "cell_barcode_column": cell_barcode_column,
            "integration_barcode_column": integration_barcode_column,
            "mutation_family_column": mutation_family_column,
            "allele_representation_threshold": allele_representation_threshold,
            "minimum_clade_fraction": minimum_clade_fraction,
            "effective_min_clade_size": effective_min_clade_size,
            "minimum_depth": minimum_depth,
            "expansion_pvalue_threshold": expansion_pvalue_threshold,
            "tree_cells": int(tree.n_cell),
            "tree_characters": int(character_matrix.shape[1]),
        }
        result = make_result(
            kind="table",
            title=f"Cassiopeia expansions: {tumor}",
            operation="cassiopeia_expansion_test",
            parameters=parameters,
            description="Node-level expansion probabilities from a VanillaGreedy Cassiopeia lineage reconstruction.",
            warnings=warnings,
            input_cells=input_cells,
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
            description="Reconstruct one tumor lineage and add its single-cell effective plasticity to AnnData.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input(
                    "allele_table_file",
                    default="openbio-singlecell/allele_table.tsv",
                ),
                io.String.Input("tumor", default=""),
                io.String.Input("annotation_key", default="cell_type"),
                io.Float.Input(
                    "allele_representation_threshold",
                    default=0.9,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Boolean.Input("first_column_as_index", default=True, advanced=True),
                io.String.Input("tumor_column", default="Tumor", advanced=True),
                io.String.Input("cell_barcode_column", default="cellBC", advanced=True),
                io.String.Input("integration_barcode_column", default="intBC", advanced=True),
                io.String.Input("mutation_family_column", default="MetFamily", advanced=True),
                io.String.Input("output_key", default="scPlasticity", advanced=True),
                io.String.Input("summary_key", default="cassiopeia_plasticity", advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
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
        adata: AnnData,
        allele_table_file: str,
        tumor: str = "",
        annotation_key: str = "cell_type",
        allele_representation_threshold: float = 0.9,
        first_column_as_index: bool = True,
        tumor_column: str = "Tumor",
        cell_barcode_column: str = "cellBC",
        integration_barcode_column: str = "intBC",
        mutation_family_column: str = "MetFamily",
        output_key: str = "scPlasticity",
        summary_key: str = "cassiopeia_plasticity",
    ) -> io.NodeOutput:
        tumor = tumor.strip()
        annotation_key = annotation_key.strip()
        output_key = output_key.strip()
        summary_key = summary_key.strip()
        if not tumor:
            raise ValueError("Tumor cannot be empty.")
        if not annotation_key or annotation_key not in adata.obs:
            raise ValueError(f"Plasticity annotation column not found in obs: {annotation_key!r}")
        if not output_key:
            raise ValueError("Plasticity output_key cannot be empty.")
        if output_key == annotation_key:
            raise ValueError("Plasticity output_key must differ from annotation_key.")
        if not summary_key:
            raise ValueError("Plasticity summary_key cannot be empty.")

        cassiopeia, _ = _require_cassiopeia()
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        allele_table = _read_allele_table(allele_table_file, first_column_as_index)
        allele_table = _normalize_cassiopeia_columns(
            allele_table,
            {
                "Tumor": tumor_column,
                "cellBC": cell_barcode_column,
                "intBC": integration_barcode_column,
                "MetFamily": mutation_family_column,
            },
        )
        tumor_allele_table = allele_table[allele_table["Tumor"].astype(str) == tumor].copy()
        if tumor_allele_table.empty:
            raise ValueError(f"Tumor not found in allele table: {tumor!r}")

        indel_priors = cassiopeia.pp.compute_empirical_indel_priors(
            allele_table,
            grouping_variables=["intBC", "MetFamily"],
        )
        character_matrix, priors, _ = cassiopeia.pp.convert_alleletable_to_character_matrix(
            tumor_allele_table,
            allele_rep_thresh=allele_representation_threshold,
            mutation_priors=indel_priors,
        )
        if character_matrix.shape[0] < 2 or character_matrix.shape[1] == 0:
            raise ValueError("The selected tumor did not produce enough cells and characters to reconstruct a tree.")

        tree = cassiopeia.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)
        cassiopeia.solver.VanillaGreedySolver().solve(tree)
        output = adata.copy()
        missing_leaves = [leaf for leaf in tree.leaves if leaf not in output.obs_names]
        if missing_leaves:
            raise ValueError(
                f"Expression AnnData is missing {len(missing_leaves)} lineage cells; "
                f"first missing cell: {missing_leaves[0]!r}."
            )

        tree.cell_meta = science.pd.DataFrame(output.obs.loc[tree.leaves, annotation_key].astype(str))
        parsimony = cassiopeia.tl.score_small_parsimony(tree, meta_item=annotation_key)
        effective_plasticity_score = float(parsimony / len(tree.nodes))

        for node in tree.depth_first_traverse_nodes():
            node_parsimony = cassiopeia.tl.score_small_parsimony(
                tree,
                meta_item=annotation_key,
                root=node,
            )
            leaf_count = len(tree.leaves_in_subtree(node))
            tree.set_attribute(node, "effective_plasticity", node_parsimony / leaf_count)

        tree.cell_meta[output_key] = 0.0
        for leaf in tree.leaves:
            ancestor_plasticities = []
            parent = tree.parent(leaf)
            while True:
                ancestor_plasticities.append(tree.get_attribute(parent, "effective_plasticity"))
                if parent == tree.root:
                    break
                parent = tree.parent(parent)
            tree.cell_meta.loc[leaf, output_key] = science.np.mean(ancestor_plasticities)

        output.obs[output_key] = science.np.nan
        output.obs.loc[tree.leaves, output_key] = tree.cell_meta[output_key]
        output.uns[summary_key] = {
            "tumor": tumor,
            "annotation_key": annotation_key,
            "parsimony": int(parsimony),
            "effective_plasticity_score": effective_plasticity_score,
            "tree_nodes": int(len(tree.nodes)),
            "tree_leaves": int(len(tree.leaves)),
            "tree_characters": int(character_matrix.shape[1]),
        }
        parameters = {
            "allele_table_file": allele_table_file,
            "tumor": tumor,
            "annotation_key": annotation_key,
            "allele_representation_threshold": allele_representation_threshold,
            "first_column_as_index": first_column_as_index,
            "tumor_column": tumor_column,
            "cell_barcode_column": cell_barcode_column,
            "integration_barcode_column": integration_barcode_column,
            "mutation_family_column": mutation_family_column,
            "output_key": output_key,
            "summary_key": summary_key,
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
    OpenBioSingleCellCassiopeiaExpansionTest,
    OpenBioSingleCellCassiopeiaPlasticity,
]


__all__ = ["LINEAGE_NODE_CLASSES"]
