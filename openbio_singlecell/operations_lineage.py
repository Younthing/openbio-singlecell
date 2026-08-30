from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from . import PLUGIN_VERSION
from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .cassiopeia_codec import (
    CHARACTERS_CODEC,
    TREE_CODEC,
    read_characters,
    read_tree,
    write_characters,
    write_tree,
)
from .cassiopeia_tree import (
    CassiopeiaCharacters,
    CassiopeiaTree,
    add_cassiopeia_plasticity,
    compute_cassiopeia_expansions,
    prepare_cassiopeia_characters,
    reconstruct_cassiopeia_tree,
)
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_artifact_input,
    require_file_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_table_output,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

CHARACTERS_KIND = "OPENBIO_CASSIOPEIA_CHARACTERS"
TREE_KIND = "OPENBIO_CASSIOPEIA_TREE"
TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"


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


def _artifact_record(
    context: OperationContext,
    *,
    name: str,
    kind: str,
    codec: str,
    writer: Any,
    value: Any,
) -> dict[str, JSONValue]:
    root = context.create_output_directory(name)
    writer(root, value)
    return {
        "type": "artifact",
        "name": name,
        "kind": kind,
        "codec": codec,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def lineage_qc_owned(path: str | Path, **parameters: Any) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    artifact, table, summary, code = prepare_cassiopeia_characters(path, **parameters, openbio_version=PLUGIN_VERSION)
    input_cells = int(summary["key_results"]["globally_unique_cells"])
    return (
        artifact,
        _table_result(
            table,
            summary,
            title="Cassiopeia character preparation QC",
            operation="prepare_cassiopeia_characters",
            started_at=started_at,
            input_cells=input_cells,
        ),
        _summary_result(
            summary,
            title="Cassiopeia character preparation summary",
            operation="prepare_cassiopeia_characters",
            started_at=started_at,
            input_cells=input_cells,
        ),
        code,
    )


def reconstruct_owned(characters: CassiopeiaCharacters, **parameters: Any) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    artifact, summary, code = reconstruct_cassiopeia_tree(
        characters,
        **parameters,
        openbio_version=PLUGIN_VERSION,
        _worker_owned=True,
    )
    report = _summary_result(
        summary,
        title=f"Cassiopeia VanillaGreedy tree: {artifact.lineage_id}",
        operation="reconstruct_cassiopeia_vanilla_greedy",
        started_at=started_at,
        input_cells=artifact.input_cells,
    )
    return artifact, report, code


def expansion_owned(tree: CassiopeiaTree, **parameters: Any) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    table, summary, code = compute_cassiopeia_expansions(
        tree,
        **parameters,
        openbio_version=PLUGIN_VERSION,
        _worker_owned=True,
    )
    return (
        _table_result(
            table,
            summary,
            title=f"Cassiopeia clade expansion evidence: {tree.lineage_id}",
            operation="cassiopeia_clade_expansion_test",
            started_at=started_at,
            input_cells=tree.input_cells,
        ),
        _summary_result(
            summary,
            title=f"Cassiopeia clade expansion summary: {tree.lineage_id}",
            operation="cassiopeia_clade_expansion_test",
            started_at=started_at,
            input_cells=tree.input_cells,
        ),
        code,
    )


def plasticity_owned(adata: Any, tree: CassiopeiaTree, **parameters: Any) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    output, table, summary, code = add_cassiopeia_plasticity(
        adata,
        tree,
        **parameters,
        openbio_version=PLUGIN_VERSION,
        _worker_owned=True,
    )
    finish_adata(
        output,
        "cassiopeia_effective_plasticity",
        summary["parameters"],
        int(adata.n_obs),
        int(adata.n_vars),
        started_at,
    )
    return (
        output,
        _table_result(
            table,
            summary,
            title=f"Cassiopeia EffectivePlasticity cells: {tree.lineage_id}",
            operation="cassiopeia_effective_plasticity",
            started_at=started_at,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
        ),
        _summary_result(
            summary,
            title=f"Cassiopeia EffectivePlasticity summary: {tree.lineage_id}",
            operation="cassiopeia_effective_plasticity",
            started_at=started_at,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
        ),
        code,
    )


@register_operation("openbio.node.cassiopeialineageqc")
def cassiopeia_lineage_qc(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"allele_table_file"}, operation="Cassiopeia Lineage QC")
    require_parameters(
        parameters,
        {
            "first_column_as_index",
            "lineage_column",
            "cell_barcode_column",
            "integration_barcode_column",
            "cut_site_columns",
            "prior_grouping_columns",
            "missing_data_allele",
            "allele_representation_threshold",
            "minimum_cells",
            "maximum_missing_fraction",
            "maximum_uncut_fraction",
            "minimum_unique_fraction",
            "minimum_informative_character_fraction",
            "max_file_mib",
            "max_matrix_gib",
        },
        operation="Cassiopeia Lineage QC",
    )
    path, _provenance = require_file_input(inputs, "allele_table_file")
    characters, table, summary, code = lineage_qc_owned(path, **parameters)
    return analysis_outputs(
        summary,
        code,
        _artifact_record(
            context,
            name="characters",
            kind=CHARACTERS_KIND,
            codec=CHARACTERS_CODEC,
            writer=write_characters,
            value=characters,
        ),
        write_table_output(context, table, kind=TABLE_KIND),
    )


@register_operation("openbio.node.reconstructcassiopeiatree")
def reconstruct_cassiopeia(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"characters"}, operation="Reconstruct Cassiopeia Tree")
    require_parameters(
        parameters,
        {"lineage_id", "prior_transformation", "collapse_mutationless_edges"},
        operation="Reconstruct Cassiopeia Tree",
    )
    root = require_artifact_input(inputs, "characters", kind=CHARACTERS_KIND, codec=CHARACTERS_CODEC)
    tree, summary, code = reconstruct_owned(read_characters(root), **parameters)
    return analysis_outputs(
        summary,
        code,
        _artifact_record(
            context,
            name="tree",
            kind=TREE_KIND,
            codec=TREE_CODEC,
            writer=write_tree,
            value=tree,
        ),
    )


@register_operation("openbio.node.cassiopeiaexpansiontest")
def cassiopeia_expansion(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"tree"}, operation="Cassiopeia Expansion Test")
    require_parameters(
        parameters,
        {"minimum_clade_size", "minimum_depth", "fdr_threshold"},
        operation="Cassiopeia Expansion Test",
    )
    root = require_artifact_input(inputs, "tree", kind=TREE_KIND, codec=TREE_CODEC)
    table, summary, code = expansion_owned(read_tree(root), **parameters)
    return analysis_outputs(
        summary,
        code,
        write_table_output(context, table, kind=TABLE_KIND),
    )


@register_operation("openbio.node.cassiopeiaplasticity")
def cassiopeia_plasticity(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "tree"}, operation="Cassiopeia Plasticity")
    require_parameters(
        parameters,
        {
            "annotation_key",
            "annotation_status",
            "analysis_mode",
            "minimum_state_fraction",
            "output_key",
            "overwrite_existing",
            "max_working_gib",
        },
        operation="Cassiopeia Plasticity",
    )
    tree_root = require_artifact_input(inputs, "tree", kind=TREE_KIND, codec=TREE_CODEC)
    adata, table, summary, code = plasticity_owned(read_anndata_input(inputs), read_tree(tree_root), **parameters)
    return analysis_outputs(
        summary,
        code,
        write_anndata_output(context, adata),
        write_table_output(context, table, kind=TABLE_KIND),
    )


__all__ = [
    "cassiopeia_expansion",
    "cassiopeia_lineage_qc",
    "cassiopeia_plasticity",
    "expansion_owned",
    "lineage_qc_owned",
    "plasticity_owned",
    "reconstruct_cassiopeia",
    "reconstruct_owned",
]
