from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from . import PLUGIN_VERSION
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata, make_plot_result, make_summary_result, make_table_result
from .artifact_envelope import table_from_metadata
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
from .lineage_plotting import (
    cassiopeia_expansion_plot_code,
    cassiopeia_lineage_qc_plot_code,
    cassiopeia_plasticity_plot_code,
    cassiopeia_tree_plot_code,
    lineage_table_fingerprint,
    run_cassiopeia_expansion_plot,
    run_cassiopeia_lineage_qc_plot,
    run_cassiopeia_plasticity_plot,
    run_cassiopeia_tree_plot,
)
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    read_table_input,
    require_artifact_input,
    require_file_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
    write_table_output,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

CHARACTERS_KIND = "OPENBIO_CASSIOPEIA_CHARACTERS"
TREE_KIND = "OPENBIO_CASSIOPEIA_TREE"
TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
CASSIOPEIA_REFERENCE = AnalysisReference(
    citation=(
        "Jones MG et al. Inference of single-cell phylogenies from lineage tracing data using Cassiopeia. "
        "Genome Biology. 2020."
    ),
    doi="10.1186/s13059-020-02000-8",
    url="https://doi.org/10.1186/s13059-020-02000-8",
    kind="method",
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)
YANG_REFERENCE = AnalysisReference(
    citation=(
        "Yang D, Jones MG et al. Lineage tracing reveals the phylodynamics, plasticity, and paths of tumor "
        "evolution. Cell. 2022."
    ),
    doi="10.1016/j.cell.2022.04.015",
    url="https://doi.org/10.1016/j.cell.2022.04.015",
    kind="method",
)


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
    parameters: dict[str, Any] | None = None,
) -> Any:
    return make_table_result(
        table=table,
        title=title,
        operation=operation,
        parameters=summary["parameters"] if parameters is None else parameters,
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


def _table_input(inputs: dict[str, JSONValue], name: str = "table") -> Any:
    table, metadata = read_table_input(inputs, name, kind=TABLE_KIND)
    return table_from_metadata(metadata, table)


def lineage_qc_owned(path: str | Path, **parameters: Any) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    artifact, table, summary, code = prepare_cassiopeia_characters(
        str(path), **parameters, openbio_version=PLUGIN_VERSION
    )
    input_cells = int(summary["key_results"]["globally_unique_cells"])
    table_parameters = {
        **summary["parameters"],
        "characters_fingerprint": artifact.fingerprint,
        "qc_table_sha256": lineage_table_fingerprint(table),
    }
    return (
        artifact,
        _table_result(
            table,
            summary,
            title="Cassiopeia character preparation QC",
            operation="prepare_cassiopeia_characters",
            started_at=started_at,
            input_cells=input_cells,
            parameters=table_parameters,
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
    table_parameters = {
        **summary["parameters"],
        "lineage_id": tree.lineage_id,
        "tree_fingerprint": tree.fingerprint,
        "topology_sha256": tree.topology_fingerprint,
        "table_sha256": lineage_table_fingerprint(table),
    }
    return (
        _table_result(
            table,
            summary,
            title=f"Cassiopeia clade expansion evidence: {tree.lineage_id}",
            operation="cassiopeia_clade_expansion_test",
            started_at=started_at,
            input_cells=tree.input_cells,
            parameters=table_parameters,
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
    table_parameters = {
        **summary["parameters"],
        "lineage_id": tree.lineage_id,
        "tree_fingerprint": tree.fingerprint,
        "topology_sha256": tree.topology_fingerprint,
        "table_sha256": lineage_table_fingerprint(table),
    }
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
            parameters=table_parameters,
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
        write_table_output(
            context,
            table,
            kind=TABLE_KIND,
            json_list_columns=("qc_warnings", "unavailability_reasons"),
        ),
    )


def lineage_qc_plot_owned(characters: Any, table: Any) -> tuple[Any, Any, str]:
    from .contracts import TableResult

    if type(characters) is not CassiopeiaCharacters:
        raise TypeError("Cassiopeia Lineage QC Plot requires an audited character artifact.")
    if not isinstance(table, TableResult):
        raise TypeError("Cassiopeia Lineage QC Plot requires the canonical QC TableResult.")
    source_parameters = table.source.get("parameters") if isinstance(table.source, dict) else None
    if (
        not isinstance(source_parameters, dict)
        or any(
            table.parameters.get(key) != source_parameters.get(key)
            for key in ("characters_fingerprint", "qc_table_sha256")
        )
        or any(getattr(table, key) != table.source.get(key) for key in ("input_cells", "input_genes", "random_seed"))
    ):
        raise ValueError("Cassiopeia Lineage QC Plot requires internally consistent table provenance.")
    producer = {
        "operation": table.source.get("operation"),
        "parameters": table.parameters,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "random_seed": table.random_seed,
    }
    started_at = time.perf_counter()
    png, details = run_cassiopeia_lineage_qc_plot(characters, table.table, producer=producer)
    parameters: dict[str, Any] = {}
    warnings = list(table.warnings)
    code = cassiopeia_lineage_qc_plot_code(producer=producer)
    description = (
        f"Rendered stored lineage QC evidence for {details['lineages_plotted']} lineage(s); "
        f"{details['lineages_available']} retained a computable character payload."
    )
    common = {
        "operation": "cassiopeia_lineage_qc_plot",
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "started_at": started_at,
        "random_seed": table.random_seed,
    }
    plotted = make_plot_result(png=png, title="Cassiopeia lineage QC", **common)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCassiopeiaLineageQCPlot",
        title="Cassiopeia lineage QC plot summary",
        operation="cassiopeia_lineage_qc_plot",
        methods=(
            "Validated the exact Cassiopeia character artifact, canonical QC table producer, table content "
            "fingerprint, lineage identities, retained character axes, and fraction domains before read-only rendering."
        ),
        results=description,
        key_results=details,
        parameters=parameters,
        references=(CASSIOPEIA_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("cassiopeia-mt", "matplotlib", "numpy", "pandas"),
        warnings=warnings,
        limitations=(
            "QC thresholds are protocol-dependent screening rules and do not establish biological validity.",
            "Lineage cells are observational descendants, not independent biological Samples.",
        ),
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        code=code,
        random_seed=table.random_seed,
    )
    return plotted, report, code


@register_operation("openbio.node.cassiopeialineageqcplot")
def cassiopeia_lineage_qc_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"characters", "table"}, operation="Cassiopeia Lineage QC Plot")
    require_parameters(parameters, set(), operation="Cassiopeia Lineage QC Plot")
    characters_root = require_artifact_input(
        inputs,
        "characters",
        kind=CHARACTERS_KIND,
        codec=CHARACTERS_CODEC,
    )
    plotted, report, code = lineage_qc_plot_owned(
        read_characters(characters_root),
        _table_input(inputs),
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


def cassiopeia_tree_plot_owned(
    tree: Any,
    adata: Any,
    *,
    annotation_key: str = "cell_type",
) -> tuple[Any, Any, str]:
    if type(tree) is not CassiopeiaTree:
        raise TypeError("Cassiopeia Tree Plot requires an audited tree artifact.")
    started_at = time.perf_counter()
    png, details = run_cassiopeia_tree_plot(tree, adata, annotation_key=annotation_key)
    parameters = {"annotation_key": details["annotation_key"]}
    code = cassiopeia_tree_plot_code(**parameters)
    description = (
        f"Rendered the stored {details['node_count']}-node Cassiopeia topology for lineage "
        f"{details['lineage_id']!r}, with {details['leaf_count']} leaves colored by {annotation_key!r}."
    )
    warnings = []
    if not details["leaf_labels_rendered"]:
        warnings.append("Leaf text labels were hidden above 80 leaves; every leaf tip remains plotted.")
    common = {
        "operation": "cassiopeia_tree_plot",
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "started_at": started_at,
    }
    plotted = make_plot_result(
        png=png,
        title=f"Cassiopeia tree: {details['lineage_id']}",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCassiopeiaTreePlot",
        title="Cassiopeia tree plot summary",
        operation="cassiopeia_tree_plot",
        methods=(
            "Validated the typed tree artifact, rooted arborescence, topology fingerprint, exact leaf identities, "
            "and categorical AnnData annotations, then rendered the stored topology without reconstruction."
        ),
        results=description,
        key_results=details,
        parameters=parameters,
        references=(CASSIOPEIA_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("anndata", "cassiopeia-mt", "matplotlib", "pandas"),
        warnings=warnings,
        limitations=(
            "The topology is a heuristic reconstruction, not a confidence interval or proof of ancestry.",
            "Leaf colors are descriptive annotations and do not establish lineage-state transitions.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


@register_operation("openbio.node.cassiopeiatreeplot")
def cassiopeia_tree_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "tree"}, operation="Cassiopeia Tree Plot")
    require_parameters(parameters, {"annotation_key"}, operation="Cassiopeia Tree Plot")
    tree_root = require_artifact_input(inputs, "tree", kind=TREE_KIND, codec=TREE_CODEC)
    plotted, report, code = cassiopeia_tree_plot_owned(
        read_tree(tree_root),
        read_anndata_input(inputs),
        **parameters,
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


def cassiopeia_expansion_plot_owned(tree: Any, table: Any) -> tuple[Any, Any, str]:
    from .contracts import TableResult

    if type(tree) is not CassiopeiaTree:
        raise TypeError("Cassiopeia Expansion Plot requires an audited tree artifact.")
    if not isinstance(table, TableResult):
        raise TypeError("Cassiopeia Expansion Plot requires the canonical expansion TableResult.")
    source_parameters = table.source.get("parameters") if isinstance(table.source, dict) else None
    binding_keys = ("lineage_id", "tree_fingerprint", "topology_sha256", "table_sha256")
    if (
        not isinstance(source_parameters, dict)
        or any(table.parameters.get(key) != source_parameters.get(key) for key in binding_keys)
        or any(getattr(table, key) != table.source.get(key) for key in ("input_cells", "input_genes", "random_seed"))
    ):
        raise ValueError("Cassiopeia Expansion Plot requires internally consistent table provenance.")
    producer = {
        "operation": table.source.get("operation"),
        "parameters": table.parameters,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "random_seed": table.random_seed,
    }
    started_at = time.perf_counter()
    png, details = run_cassiopeia_expansion_plot(tree, table.table, producer=producer)
    parameters: dict[str, Any] = {}
    warnings = list(table.warnings)
    code = cassiopeia_expansion_plot_code(producer=producer)
    description = (
        f"Rendered all {details['topology_nodes_plotted']} stored topology nodes and "
        f"{details['eligible_clades']} eligible expansion tests for lineage {details['lineage_id']!r}; "
        f"{details['significant_clades']} met the stored FDR threshold."
    )
    common = {
        "operation": "cassiopeia_expansion_plot",
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "started_at": started_at,
        "random_seed": table.random_seed,
    }
    plotted = make_plot_result(
        png=png,
        title=f"Cassiopeia expansion: {details['lineage_id']}",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCassiopeiaExpansionPlot",
        title="Cassiopeia expansion plot summary",
        operation="cassiopeia_expansion_plot",
        methods=(
            "Validated the exact tree and topology fingerprints, canonical full expansion table, node-parent "
            "relationships, descendant-leaf counts, eligible hypothesis family, and FDR calls before read-only rendering."
        ),
        results=description,
        key_results=details,
        parameters=parameters,
        references=(YANG_REFERENCE, CASSIOPEIA_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("cassiopeia-mt", "matplotlib", "numpy", "pandas"),
        warnings=warnings,
        limitations=(
            "Clades are nested within one inferred tree and are not independent biological Samples.",
            "Expansion significance is conditional on the reconstructed topology and stored hypothesis family.",
        ),
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        code=code,
        random_seed=table.random_seed,
    )
    return plotted, report, code


@register_operation("openbio.node.cassiopeiaexpansionplot")
def cassiopeia_expansion_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"table", "tree"}, operation="Cassiopeia Expansion Plot")
    require_parameters(parameters, set(), operation="Cassiopeia Expansion Plot")
    tree_root = require_artifact_input(inputs, "tree", kind=TREE_KIND, codec=TREE_CODEC)
    plotted, report, code = cassiopeia_expansion_plot_owned(
        read_tree(tree_root),
        _table_input(inputs),
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


def cassiopeia_plasticity_plot_owned(
    adata: Any,
    tree: Any,
    table: Any,
    *,
    output_key: str = "sc_effective_plasticity",
) -> tuple[Any, Any, str]:
    from .contracts import TableResult

    if type(tree) is not CassiopeiaTree:
        raise TypeError("Cassiopeia Plasticity Plot requires an audited tree artifact.")
    if not isinstance(table, TableResult):
        raise TypeError("Cassiopeia Plasticity Plot requires the canonical plasticity TableResult.")
    source_parameters = table.source.get("parameters") if isinstance(table.source, dict) else None
    binding_keys = ("lineage_id", "tree_fingerprint", "topology_sha256", "table_sha256", "output_key")
    if (
        not isinstance(source_parameters, dict)
        or any(table.parameters.get(key) != source_parameters.get(key) for key in binding_keys)
        or any(getattr(table, key) != table.source.get(key) for key in ("input_cells", "input_genes", "random_seed"))
    ):
        raise ValueError("Cassiopeia Plasticity Plot requires internally consistent table provenance.")
    producer = {
        "operation": table.source.get("operation"),
        "parameters": table.parameters,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "random_seed": table.random_seed,
    }
    started_at = time.perf_counter()
    png, details = run_cassiopeia_plasticity_plot(
        adata,
        tree,
        table.table,
        producer=producer,
        output_key=output_key,
    )
    parameters = {"output_key": details["output_key"]}
    warnings = list(table.warnings)
    if details["tree_leaves_without_scores"]:
        warnings.append(
            f"{details['tree_leaves_without_scores']} tree leaf/leaves lacked an included plasticity score and "
            "are shown in gray."
        )
    code = cassiopeia_plasticity_plot_code(producer=producer, output_key=output_key)
    description = (
        f"Rendered {details['included_cells']} stored single-cell EffectivePlasticity scores across "
        f"{len(details['annotation_categories'])} retained annotation state(s) for lineage {details['lineage_id']!r}."
    )
    common = {
        "operation": "cassiopeia_plasticity_plot",
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": table.input_cells,
        "input_genes": table.input_genes,
        "started_at": started_at,
        "random_seed": table.random_seed,
    }
    plotted = make_plot_result(
        png=png,
        title=f"Cassiopeia EffectivePlasticity: {details['lineage_id']}",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCassiopeiaPlasticityPlot",
        title="Cassiopeia plasticity plot summary",
        operation="cassiopeia_plasticity_plot",
        methods=(
            "Validated the exact tree and topology fingerprints, canonical full cell table, AnnData observation "
            "axis, annotation and status columns, stored EffectivePlasticity provenance, and per-cell score equality "
            "before read-only rendering."
        ),
        results=description,
        key_results=details,
        parameters=parameters,
        references=(YANG_REFERENCE, CASSIOPEIA_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("anndata", "cassiopeia-mt", "matplotlib", "numpy", "pandas"),
        warnings=warnings,
        limitations=(
            "EffectivePlasticity is descriptive and conditional on one inferred tree and supplied annotation.",
            "Cells and nested subtrees are not independent biological Samples; no Condition inference is shown.",
        ),
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        code=code,
        random_seed=table.random_seed,
    )
    return plotted, report, code


@register_operation("openbio.node.cassiopeiaplasticityplot")
def cassiopeia_plasticity_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "table", "tree"}, operation="Cassiopeia Plasticity Plot")
    require_parameters(parameters, {"output_key"}, operation="Cassiopeia Plasticity Plot")
    tree_root = require_artifact_input(inputs, "tree", kind=TREE_KIND, codec=TREE_CODEC)
    plotted, report, code = cassiopeia_plasticity_plot_owned(
        read_anndata_input(inputs),
        read_tree(tree_root),
        _table_input(inputs),
        **parameters,
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


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
    "cassiopeia_expansion_plot",
    "cassiopeia_expansion_plot_owned",
    "cassiopeia_lineage_qc",
    "cassiopeia_lineage_qc_plot",
    "cassiopeia_plasticity",
    "cassiopeia_plasticity_plot",
    "cassiopeia_plasticity_plot_owned",
    "cassiopeia_tree_plot",
    "cassiopeia_tree_plot_owned",
    "expansion_owned",
    "lineage_qc_owned",
    "lineage_qc_plot_owned",
    "plasticity_owned",
    "reconstruct_cassiopeia",
    "reconstruct_owned",
]
