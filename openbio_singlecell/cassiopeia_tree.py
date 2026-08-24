from __future__ import annotations

import copy
import importlib
import os
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from . import dependencies
from .files import resolve_input_path

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


ALLELE_TABLE_EXTENSIONS = (".txt", ".tsv", ".csv")


def _require_cassiopeia() -> tuple[Any, Any]:
    try:
        cassiopeia = importlib.import_module("cassiopeia")
        lineage_utils = importlib.import_module("cassiopeia.preprocess.lineage_utils")
    except (ImportError, OSError) as error:
        raise RuntimeError("Cassiopeia lineage nodes require the cassiopeia-lineage package.") from error
    return cassiopeia, lineage_utils


def _read_allele_table(path: str, first_column_as_index: bool) -> DataFrame:
    science = dependencies.require_scientific_dependencies()
    resolved = resolve_input_path(path, extensions=ALLELE_TABLE_EXTENSIONS)
    separator = "," if os.path.splitext(resolved)[1].lower() == ".csv" else "\t"
    table = science.pd.read_csv(resolved, sep=separator, index_col=0 if first_column_as_index else None)
    if table.empty:
        raise ValueError("The allele table is empty.")
    return table


def _normalize_columns(table: DataFrame, columns: Mapping[str, str]) -> DataFrame:
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


@dataclass(frozen=True, slots=True, eq=False)
class CassiopeiaTree:
    tumor: str
    input_cells: int
    character_count: int
    provenance: Mapping[str, Any]
    _solved_tree: Any = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        tumor = self.tumor.strip()
        if not tumor:
            raise ValueError("Cassiopeia tree tumor cannot be empty.")
        if self.input_cells < 2:
            raise ValueError("Cassiopeia tree requires at least two input cells.")
        if self.character_count < 1:
            raise ValueError("Cassiopeia tree requires at least one character.")
        if not isinstance(self.provenance, Mapping):
            raise TypeError("Cassiopeia tree provenance must be a mapping.")
        if not self.provenance:
            raise ValueError("Cassiopeia tree provenance cannot be empty.")

        required_attributes = ("n_cell", "root", "leaves", "nodes")
        required_methods = (
            "copy",
            "depth_first_traverse_nodes",
            "get_attribute",
            "leaves_in_subtree",
            "parent",
            "set_attribute",
        )
        missing = [name for name in required_attributes if not hasattr(self._solved_tree, name)]
        missing.extend(name for name in required_methods if not callable(getattr(self._solved_tree, name, None)))
        if missing:
            raise TypeError(f"Cassiopeia tree is missing required solved-tree capabilities: {sorted(missing)}")
        if int(self._solved_tree.n_cell) < 2:
            raise ValueError("Cassiopeia solved tree requires at least two cells.")
        if self._solved_tree.root is None or len(self._solved_tree.nodes) < 2 or len(self._solved_tree.leaves) < 2:
            raise ValueError("Cassiopeia solver did not produce a usable rooted tree.")

        object.__setattr__(self, "tumor", tumor)
        object.__setattr__(self, "provenance", MappingProxyType(copy.deepcopy(dict(self.provenance))))


def reconstruct_cassiopeia_tree(
    *,
    allele_table_file: str,
    tumor: str,
    first_column_as_index: bool,
    tumor_column: str,
    cell_barcode_column: str,
    integration_barcode_column: str,
    mutation_family_column: str,
    allele_representation_threshold: float,
) -> CassiopeiaTree:
    tumor = tumor.strip()
    if not tumor:
        raise ValueError("Tumor cannot be empty.")
    if not 0.0 <= allele_representation_threshold <= 1.0:
        raise ValueError("allele_representation_threshold must be between 0 and 1.")

    started_at = time.perf_counter()
    cassiopeia, _ = _require_cassiopeia()
    allele_table = _read_allele_table(allele_table_file, first_column_as_index)
    allele_table = _normalize_columns(
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
        raise ValueError("Cassiopeia tree reconstruction requires at least two cells in the selected tumor.")

    indel_priors = cassiopeia.pp.compute_empirical_indel_priors(
        allele_table,
        grouping_variables=["intBC", "MetFamily"],
    )
    character_matrix, priors, _ = cassiopeia.pp.convert_alleletable_to_character_matrix(
        tumor_allele_table,
        allele_rep_thresh=allele_representation_threshold,
        mutation_priors=indel_priors,
    )
    character_count = int(character_matrix.shape[1])
    if character_matrix.shape[0] < 2 or character_count == 0:
        raise ValueError("The selected tumor did not produce enough cells and characters to reconstruct a tree.")

    solved_tree = cassiopeia.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)
    cassiopeia.solver.VanillaGreedySolver().solve(solved_tree)
    parameters = {
        "allele_table_file": allele_table_file,
        "tumor": tumor,
        "first_column_as_index": first_column_as_index,
        "tumor_column": tumor_column,
        "cell_barcode_column": cell_barcode_column,
        "integration_barcode_column": integration_barcode_column,
        "mutation_family_column": mutation_family_column,
        "allele_representation_threshold": allele_representation_threshold,
    }
    provenance = {
        "operation": "reconstruct_cassiopeia_tree",
        "parameters": parameters,
        "input_cells": input_cells,
        "character_count": character_count,
        "elapsed_seconds": time.perf_counter() - started_at,
        "timestamp": datetime.now(UTC).isoformat(),
    }
    return CassiopeiaTree(
        tumor=tumor,
        input_cells=input_cells,
        character_count=character_count,
        provenance=provenance,
        _solved_tree=solved_tree,
    )


def _working_tree(tree: CassiopeiaTree) -> Any:
    if not isinstance(tree, CassiopeiaTree):
        raise TypeError("Expected an OPENBIO_CASSIOPEIA_TREE value.")
    return tree._solved_tree.copy()


def compute_cassiopeia_expansions(
    tree: CassiopeiaTree,
    *,
    minimum_clade_fraction: float,
    minimum_depth: int,
    expansion_pvalue_threshold: float,
) -> tuple[DataFrame, float]:
    if not 0.0 <= minimum_clade_fraction <= 1.0:
        raise ValueError("minimum_clade_fraction must be between 0 and 1.")
    if minimum_depth < 0:
        raise ValueError("minimum_depth cannot be negative.")
    if not 0.0 <= expansion_pvalue_threshold <= 1.0:
        raise ValueError("expansion_pvalue_threshold must be between 0 and 1.")

    cassiopeia, _ = _require_cassiopeia()
    science = dependencies.require_scientific_dependencies()
    solved_tree = _working_tree(tree)
    effective_min_clade_size = minimum_clade_fraction * solved_tree.n_cell
    cassiopeia.tl.compute_expansion_pvalues(
        solved_tree,
        min_clade_size=effective_min_clade_size,
        min_depth=minimum_depth,
    )

    rows = []
    for node in solved_tree.depth_first_traverse_nodes():
        expansion_pvalue = float(solved_tree.get_attribute(node, "expansion_pvalue"))
        depth = 0
        parent = None
        if node != solved_tree.root:
            parent = solved_tree.parent(node)
            current = parent
            depth = 1
            while current != solved_tree.root:
                current = solved_tree.parent(current)
                depth += 1
        rows.append(
            {
                "tumor": tree.tumor,
                "node": str(node),
                "parent": None if parent is None else str(parent),
                "depth": depth,
                "leaf_count": len(solved_tree.leaves_in_subtree(node)),
                "expansion_pvalue": expansion_pvalue,
                "is_expansion": expansion_pvalue < expansion_pvalue_threshold,
            }
        )
    return science.pd.DataFrame.from_records(rows), float(effective_min_clade_size)


def add_cassiopeia_plasticity(
    tree: CassiopeiaTree,
    adata: AnnData,
    *,
    annotation_key: str,
    output_key: str,
    summary_key: str,
) -> AnnData:
    annotation_key = annotation_key.strip()
    output_key = output_key.strip()
    summary_key = summary_key.strip()
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
    solved_tree = _working_tree(tree)
    output = adata.copy()
    missing_leaves = [leaf for leaf in solved_tree.leaves if leaf not in output.obs_names]
    if missing_leaves:
        raise ValueError(
            f"Expression AnnData is missing {len(missing_leaves)} lineage cells; "
            f"first missing cell: {missing_leaves[0]!r}."
        )

    solved_tree.cell_meta = science.pd.DataFrame(output.obs.loc[solved_tree.leaves, annotation_key].astype(str))
    parsimony = cassiopeia.tl.score_small_parsimony(solved_tree, meta_item=annotation_key)
    effective_plasticity_score = float(parsimony / len(solved_tree.nodes))

    for node in solved_tree.depth_first_traverse_nodes():
        node_parsimony = cassiopeia.tl.score_small_parsimony(
            solved_tree,
            meta_item=annotation_key,
            root=node,
        )
        leaf_count = len(solved_tree.leaves_in_subtree(node))
        solved_tree.set_attribute(node, "effective_plasticity", node_parsimony / leaf_count)

    solved_tree.cell_meta[output_key] = 0.0
    for leaf in solved_tree.leaves:
        ancestor_plasticities = []
        parent = solved_tree.parent(leaf)
        while True:
            ancestor_plasticities.append(solved_tree.get_attribute(parent, "effective_plasticity"))
            if parent == solved_tree.root:
                break
            parent = solved_tree.parent(parent)
        solved_tree.cell_meta.loc[leaf, output_key] = science.np.mean(ancestor_plasticities)

    output.obs[output_key] = science.np.nan
    output.obs.loc[solved_tree.leaves, output_key] = solved_tree.cell_meta[output_key]
    output.uns[summary_key] = {
        "tumor": tree.tumor,
        "annotation_key": annotation_key,
        "parsimony": int(parsimony),
        "effective_plasticity_score": effective_plasticity_score,
        "tree_nodes": int(len(solved_tree.nodes)),
        "tree_leaves": int(len(solved_tree.leaves)),
        "tree_characters": tree.character_count,
    }
    return output


__all__ = [
    "ALLELE_TABLE_EXTENSIONS",
    "CassiopeiaTree",
    "add_cassiopeia_plasticity",
    "compute_cassiopeia_expansions",
    "reconstruct_cassiopeia_tree",
]
