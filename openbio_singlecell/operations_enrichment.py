from __future__ import annotations

import time
from typing import Any

from . import PLUGIN_VERSION, dependencies
from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .artifact_envelope import table_from_metadata
from .dgidb_artifact_codec import DGIDB_CODEC, read_dgidb_resource, write_dgidb_resource
from .dgidb_resource import dgidb_resource_code, load_dgidb_resource, validate_dgidb_resource
from .drug_enrichment import drug_gsea_code, drug_ora_code, run_drug_gsea, run_drug_ora
from .drug_score import drug_score_code, run_drug_score
from .enrichment_artifacts import validate_enrichment_artifact_pair
from .expression_source import _AUCELL_SPEC, _DRUG_SCORE_SPEC, _GENE_PANEL_SPEC, _GSVA_SPEC, DynamicExpressionSource
from .gene_set_scoring import (
    gene_set_scoring_code,
    run_aucell_scores,
    run_gene_panel_score,
    run_gsva_scores,
)
from .generic_ora import build_generic_ora_summary, generic_ora_code, run_generic_ora_evidence
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    read_table_input,
    require_artifact_input,
    require_file_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_table_output,
)
from .pathway_score_contrast import pathway_score_contrast_code, run_pathway_score_contrast
from .ranked_enrichment import build_ranked_gsea_summary, ranked_gsea_code, run_ranked_gsea_evidence
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

DGIDB_KIND = "OPENBIO_DGIDB_RESOURCE"
TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
TABLE_CODEC = "table-jsonl-v1"
AUCELL_EXPRESSION_SOURCE = _AUCELL_SPEC
GSVA_EXPRESSION_SOURCE = _GSVA_SPEC
GENE_PANEL_EXPRESSION_SOURCE = _GENE_PANEL_SPEC
DRUG_EXPRESSION_SOURCE = _DRUG_SCORE_SPEC


def _requested_path(provenance: dict[str, JSONValue], *, resource: str) -> str:
    requested = provenance.get("path")
    if not isinstance(requested, str):
        raise ProtocolError(f"{resource} provenance path must be a string.")
    return requested


def _table_input(inputs: dict[str, JSONValue], name: str) -> Any:
    table, metadata = read_table_input(inputs, name, kind=TABLE_KIND)
    return table_from_metadata(metadata, table)


def _dgidb_input(inputs: dict[str, JSONValue]) -> Any:
    root = require_artifact_input(inputs, "resource", kind=DGIDB_KIND, codec=DGIDB_CODEC)
    return read_dgidb_resource(root)


def _anndata_records(context: OperationContext, output: Any, report: Any, code: str) -> list[JSONValue]:
    return analysis_outputs(report, code, write_anndata_output(context, output))


def _table_records(context: OperationContext, result: Any, report: Any, code: str) -> list[JSONValue]:
    return analysis_outputs(report, code, write_table_output(context, result, kind=TABLE_KIND))


def load_dgidb_resource_owned(
    resource_path: str,
    *,
    requested_path: str,
    resource_metadata_json: str = "{}",
    drug_column: str = "drug_claim_name",
    gene_column: str = "gene_claim_name",
    source_column: str = "interaction_claim_source",
    evidence_column: str = "interaction_types",
    max_file_bytes: int = 536_870_912,
    max_rows: int = 2_000_000,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    resource, summary = load_dgidb_resource(
        resource_path,
        requested_path=requested_path,
        resource_metadata_json=resource_metadata_json,
        drug_column=drug_column,
        gene_column=gene_column,
        source_column=source_column,
        evidence_column=evidence_column,
        max_file_bytes=max_file_bytes,
        max_rows=max_rows,
        openbio_version=PLUGIN_VERSION,
    )
    _, _, accounting, artifact = validate_dgidb_resource(resource, copy_payload=False)
    report = make_summary_result(
        summary=summary,
        title=f"DGIdb resource {summary['key_results']['resource']['metadata']['version']}",
        operation="load_dgidb_resource",
        parameters=dict(summary["parameters"]),
        description=summary["results"],
        warnings=list(summary["warnings"]),
        input_cells=0,
        input_genes=int(accounting["gene_count"]),
        started_at=started_at,
    )
    code = dgidb_resource_code(
        requested_path=requested_path,
        resource_metadata_json=resource_metadata_json,
        drug_column=drug_column,
        gene_column=gene_column,
        source_column=source_column,
        evidence_column=evidence_column,
        max_file_bytes=max_file_bytes,
        max_rows=max_rows,
        expected_sha256=artifact["raw_file_sha256"],
        openbio_version=PLUGIN_VERSION,
    )
    return resource, report, code


@register_operation("openbio.node.dgidbannotation")
def load_dgidb_resource_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"resource_file"}, operation="Load DGIdb Resource")
    require_parameters(
        parameters,
        {
            "resource_metadata_json",
            "drug_column",
            "gene_column",
            "source_column",
            "evidence_column",
            "max_file_bytes",
            "max_rows",
        },
        operation="Load DGIdb Resource",
    )
    path, provenance = require_file_input(inputs, "resource_file")
    requested_path = _requested_path(provenance, resource="DGIdb resource")
    resource, report, code = load_dgidb_resource_owned(
        str(path),
        requested_path=requested_path,
        **parameters,
    )
    root = context.create_output_directory("resource")
    write_dgidb_resource(root, resource)
    return analysis_outputs(
        report,
        code,
        {
            "type": "artifact",
            "name": "resource",
            "kind": DGIDB_KIND,
            "codec": DGIDB_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
    )


def aucell_scores_owned(
    adata: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str = "{}",
    source: DynamicExpressionSource | None = None,
    source_column: str = "geneset",
    target_column: str = "genesymbol",
    min_targets: int = 5,
    n_top_features: int = 0,
    batch_size: int = 250_000,
    max_output_rows: int = 2_000_000,
    max_working_memory_gib: float = 4.0,
    output_key: str = "aucell_scores",
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    expression = AUCELL_EXPRESSION_SOURCE.resolve(adata, source)
    output, summary = run_aucell_scores(
        adata,
        gene_sets_path=gene_sets_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        source_column=source_column,
        target_column=target_column,
        min_targets=min_targets,
        n_top_features=n_top_features,
        batch_size=batch_size,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        output_key=output_key,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
    )
    report = make_summary_result(
        summary=summary,
        title="AUCell observation scores",
        operation="aucell_scores",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    code = gene_set_scoring_code(
        function_name="score_aucell",
        parameters={
            "method": "aucell",
            "gene_sets_path": gene_sets_path,
            "requested_resource_path": requested_resource_path,
            "resource_metadata_json": resource_metadata_json,
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "source_column": source_column,
            "target_column": target_column,
            "min_targets": min_targets,
            "output_key": output_key,
            "overwrite_existing": overwrite_existing,
            "expected_resource_sha256": summary["key_results"]["resource"]["sha256"],
            "n_top_features": n_top_features,
            "batch_size": batch_size,
            "max_output_rows": max_output_rows,
            "max_working_memory_gib": max_working_memory_gib,
            "openbio_version": PLUGIN_VERSION,
        },
    )
    return output, report, code


@register_operation("openbio.node.aucellscores")
def aucell_scores_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "gene_sets_file"}, operation="AUCell Scores")
    require_parameters(
        parameters,
        {
            "resource_metadata_json",
            "source",
            "source_column",
            "target_column",
            "min_targets",
            "n_top_features",
            "batch_size",
            "max_output_rows",
            "max_working_memory_gib",
            "output_key",
            "overwrite_existing",
        },
        operation="AUCell Scores",
    )
    path, provenance = require_file_input(inputs, "gene_sets_file")
    return _anndata_records(
        context,
        *aucell_scores_owned(
            read_anndata_input(inputs),
            gene_sets_path=str(path),
            requested_resource_path=_requested_path(provenance, resource="Gene-set resource"),
            **parameters,
        ),
    )


def gsva_scores_owned(
    adata: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str = "{}",
    source: DynamicExpressionSource | None = None,
    source_column: str = "geneset",
    target_column: str = "genesymbol",
    min_targets: int = 10,
    kernel: str = "gaussian_normalized",
    maxdiff: bool = True,
    absrnk: bool = False,
    tau: float = 1.0,
    max_working_memory_gib: float = 4.0,
    batch_size: int = 250_000,
    max_output_rows: int = 2_000_000,
    output_key: str = "gsva_scores",
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    expression = GSVA_EXPRESSION_SOURCE.resolve(adata, source)
    output, summary = run_gsva_scores(
        adata,
        gene_sets_path=gene_sets_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        source_column=source_column,
        target_column=target_column,
        min_targets=min_targets,
        kernel=kernel,
        maxdiff=maxdiff,
        absrnk=absrnk,
        tau=tau,
        batch_size=batch_size,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        output_key=output_key,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
    )
    report = make_summary_result(
        summary=summary,
        title="GSVA observation scores",
        operation="gsva_scores",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    code = gene_set_scoring_code(
        function_name="score_gsva",
        parameters={
            "method": "gsva",
            "gene_sets_path": gene_sets_path,
            "requested_resource_path": requested_resource_path,
            "resource_metadata_json": resource_metadata_json,
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "source_column": source_column,
            "target_column": target_column,
            "min_targets": min_targets,
            "output_key": output_key,
            "overwrite_existing": overwrite_existing,
            "expected_resource_sha256": summary["key_results"]["resource"]["sha256"],
            "kernel": kernel,
            "maxdiff": maxdiff,
            "absrnk": absrnk,
            "tau": tau,
            "batch_size": batch_size,
            "max_output_rows": max_output_rows,
            "max_working_memory_gib": max_working_memory_gib,
            "openbio_version": PLUGIN_VERSION,
        },
    )
    return output, report, code


@register_operation("openbio.node.gsvascores")
def gsva_scores_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "gene_sets_file"}, operation="GSVA Scores")
    require_parameters(
        parameters,
        {
            "resource_metadata_json",
            "source",
            "source_column",
            "target_column",
            "min_targets",
            "kernel",
            "maxdiff",
            "absrnk",
            "tau",
            "max_working_memory_gib",
            "batch_size",
            "max_output_rows",
            "output_key",
            "overwrite_existing",
        },
        operation="GSVA Scores",
    )
    path, provenance = require_file_input(inputs, "gene_sets_file")
    return _anndata_records(
        context,
        *gsva_scores_owned(
            read_anndata_input(inputs),
            gene_sets_path=str(path),
            requested_resource_path=_requested_path(provenance, resource="Gene-set resource"),
            **parameters,
        ),
    )


def gene_panel_scores_owned(
    adata: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str = "{}",
    panel: str = "",
    source: DynamicExpressionSource | None = None,
    source_column: str = "geneset",
    target_column: str = "genesymbol",
    output_key: str = "panel_score",
    ctrl_size: int = 0,
    n_bins: int = 25,
    random_seed: int = 0,
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    expression = GENE_PANEL_EXPRESSION_SOURCE.resolve(adata, source)
    output, summary = run_gene_panel_score(
        adata,
        gene_sets_path=gene_sets_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        panel=panel,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        source_column=source_column,
        target_column=target_column,
        output_key=output_key,
        ctrl_size=ctrl_size,
        n_bins=n_bins,
        random_seed=random_seed,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
    )
    report = make_summary_result(
        summary=summary,
        title=f"Gene panel score: {panel}",
        operation="gene_panel_score",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    code = gene_set_scoring_code(
        function_name="score_gene_panel",
        parameters={
            "method": "panel",
            "gene_sets_path": gene_sets_path,
            "requested_resource_path": requested_resource_path,
            "resource_metadata_json": resource_metadata_json,
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "source_column": source_column,
            "target_column": target_column,
            "min_targets": 1,
            "output_key": output_key,
            "overwrite_existing": overwrite_existing,
            "expected_resource_sha256": summary["key_results"]["resource"]["sha256"],
            "panel": panel,
            "ctrl_size": ctrl_size,
            "n_bins": n_bins,
            "random_seed": random_seed,
            "openbio_version": PLUGIN_VERSION,
        },
    )
    return output, report, code


@register_operation("openbio.node.genepanelscores")
def gene_panel_scores_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "gene_sets_file"}, operation="Gene Panel Score")
    require_parameters(
        parameters,
        {
            "resource_metadata_json",
            "panel",
            "source",
            "source_column",
            "target_column",
            "output_key",
            "ctrl_size",
            "n_bins",
            "random_seed",
            "overwrite_existing",
        },
        operation="Gene Panel Score",
    )
    path, provenance = require_file_input(inputs, "gene_sets_file")
    output, report, code = gene_panel_scores_owned(
        read_anndata_input(inputs),
        gene_sets_path=str(path),
        requested_resource_path=_requested_path(provenance, resource="Gene-set resource"),
        **parameters,
    )
    return _anndata_records(context, output, report, code)


def pathway_score_ttest_owned(
    adata: Any,
    *,
    sample_key: str = "sample",
    condition_key: str = "condition",
    annotation_key: str = "cell_type",
    population: str = "",
    condition_a: str = "",
    condition_b: str = "",
    score_key: str = "aucell_scores",
    min_cells_per_sample_population: int = 10,
    min_samples_per_condition: int = 3,
    technical_batch_key: str = "",
    confidence_level: float = 0.95,
    annotation_status: str = "unknown",
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    table, summary = run_pathway_score_contrast(
        adata,
        sample_key=sample_key,
        condition_key=condition_key,
        annotation_key=annotation_key,
        population=population,
        condition_a=condition_a,
        condition_b=condition_b,
        score_key=score_key,
        min_cells_per_sample_population=min_cells_per_sample_population,
        min_samples_per_condition=min_samples_per_condition,
        technical_batch_key=technical_batch_key,
        confidence_level=confidence_level,
        annotation_status=annotation_status,
        openbio_version=PLUGIN_VERSION,
    )
    result = make_table_result(
        title=f"{condition_a} vs {condition_b} pathway scores in {population}",
        operation="pathway_score_sample_welch",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        table=table,
    )
    report = make_summary_result(
        summary=summary,
        title=f"Pathway score contrast: {condition_a} vs {condition_b}",
        operation="pathway_score_sample_welch",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    code = pathway_score_contrast_code(
        sample_key=sample_key,
        condition_key=condition_key,
        annotation_key=annotation_key,
        population=population,
        condition_a=condition_a,
        condition_b=condition_b,
        score_key=score_key,
        min_cells_per_sample_population=min_cells_per_sample_population,
        min_samples_per_condition=min_samples_per_condition,
        technical_batch_key=technical_batch_key,
        confidence_level=confidence_level,
        annotation_status=annotation_status,
        openbio_version=PLUGIN_VERSION,
    )
    return result, report, code


@register_operation("openbio.node.pathwayscorettest")
def pathway_score_ttest_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Pathway Score T-Test")
    require_parameters(
        parameters,
        {
            "sample_key",
            "condition_key",
            "annotation_key",
            "population",
            "condition_a",
            "condition_b",
            "score_key",
            "min_cells_per_sample_population",
            "min_samples_per_condition",
            "technical_batch_key",
            "confidence_level",
            "annotation_status",
        },
        operation="Pathway Score T-Test",
    )
    return _table_records(context, *pathway_score_ttest_owned(read_anndata_input(inputs), **parameters))


def ranked_gsea_owned(
    table: Any,
    universe: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str = "{}",
    comparison: str = "",
    gene_column: str = "gene",
    score_column: str = "score",
    source_column: str = "source",
    target_column: str = "target",
    min_targets: int = 15,
    max_targets: int = 500,
    n_permutations: int = 1000,
    random_seed: int = 123,
    max_output_rows: int = 100_000,
) -> tuple[Any, Any, str]:
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    provenance = validate_enrichment_artifact_pair(
        table,
        universe,
        purpose="ranked",
        selector=comparison,
        gene_column=gene_column,
        score_column=score_column,
        np=science.np,
        pd=science.pd,
    )
    runtime_parameters = {
        "artifact_family": provenance["artifact_family"],
        "comparison": provenance["comparison"],
        "comparison_column": provenance["comparison_column"],
        "gene_column": provenance["gene_column"],
        "score_column": provenance["score_column"],
        "evidence_scope": provenance["evidence_scope"],
        "inference_unit": provenance["inference_unit"],
        "direction": provenance["direction"],
        "replicate_aware": provenance["replicate_aware"],
        "technical_batch_handling": provenance["technical_batch_handling"],
        "upstream_parameters": provenance["upstream_parameters"],
        "provenance_warnings": provenance["provenance_warnings"],
        "provenance_declarations": provenance["provenance_declarations"],
        "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
        "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
        "expected_universe_fingerprint": provenance["universe_fingerprint"],
        "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
        "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
        "resource_path": gene_sets_path,
        "requested_resource_path": requested_resource_path,
        "resource_metadata_json": resource_metadata_json,
        "source_column": source_column,
        "target_column": target_column,
        "min_targets": min_targets,
        "max_targets": max_targets,
        "n_permutations": n_permutations,
        "random_seed": random_seed,
        "max_output_rows": max_output_rows,
        "openbio_version": PLUGIN_VERSION,
    }
    evidence, diagnostics = run_ranked_gsea_evidence(
        table,
        universe,
        **runtime_parameters,
        decoupler_module=None,
    )
    summary = build_ranked_gsea_summary(diagnostics)
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    result = make_table_result(
        title=f"Ranked GSEA: {provenance['comparison']}",
        operation="ranked_gsea",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        random_seed=random_seed,
        table=evidence,
    )
    report = make_summary_result(
        summary=summary,
        title=f"Ranked GSEA summary: {provenance['comparison']}",
        operation="ranked_gsea",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    code = ranked_gsea_code(
        **runtime_parameters,
        expected_resource_sha256=diagnostics["resource"]["sha256"],
    )
    return result, report, code


@register_operation("openbio.node.rankedgsea")
def ranked_gsea_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"table", "universe", "gene_sets_file"}, operation="Ranked GSEA")
    require_parameters(
        parameters,
        {
            "resource_metadata_json",
            "comparison",
            "gene_column",
            "score_column",
            "source_column",
            "target_column",
            "min_targets",
            "max_targets",
            "n_permutations",
            "random_seed",
            "max_output_rows",
        },
        operation="Ranked GSEA",
    )
    path, provenance = require_file_input(inputs, "gene_sets_file")
    return _table_records(
        context,
        *ranked_gsea_owned(
            _table_input(inputs, "table"),
            _table_input(inputs, "universe"),
            gene_sets_path=str(path),
            requested_resource_path=_requested_path(provenance, resource="Gene-set resource"),
            **parameters,
        ),
    )


def gene_set_overrepresentation_owned(
    table: Any,
    universe: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str = "{}",
    comparison: str = "",
    gene_column: str = "gene",
    source_column: str = "source",
    target_column: str = "target",
    min_targets: int = 3,
    min_overlap: int = 2,
    max_p_adjusted: float = 0.05,
    max_output_rows: int = 100_000,
) -> tuple[Any, Any, str]:
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    provenance = validate_enrichment_artifact_pair(
        table,
        universe,
        purpose="selected",
        selector=comparison,
        gene_column=gene_column,
        score_column=None,
        np=science.np,
        pd=science.pd,
    )
    runtime_parameters = {
        "artifact_family": provenance["artifact_family"],
        "comparison": provenance["comparison"],
        "comparison_column": provenance["comparison_column"],
        "gene_column": provenance["gene_column"],
        "evidence_scope": provenance["evidence_scope"],
        "inference_unit": provenance["inference_unit"],
        "direction": provenance["direction"],
        "replicate_aware": provenance["replicate_aware"],
        "technical_batch_handling": provenance["technical_batch_handling"],
        "upstream_parameters": provenance["upstream_parameters"],
        "provenance_warnings": provenance["provenance_warnings"],
        "provenance_declarations": provenance["provenance_declarations"],
        "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
        "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
        "expected_universe_fingerprint": provenance["universe_fingerprint"],
        "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
        "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
        "resource_path": gene_sets_path,
        "requested_resource_path": requested_resource_path,
        "resource_metadata_json": resource_metadata_json,
        "source_column": source_column,
        "target_column": target_column,
        "min_targets": min_targets,
        "min_overlap": min_overlap,
        "max_p_adjusted": max_p_adjusted,
        "max_output_rows": max_output_rows,
        "openbio_version": PLUGIN_VERSION,
    }
    evidence, diagnostics = run_generic_ora_evidence(
        table,
        universe,
        **runtime_parameters,
        decoupler_module=None,
    )
    summary = build_generic_ora_summary(diagnostics)
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    result = make_table_result(
        title=f"Gene-set overrepresentation: {provenance['comparison']}",
        operation="gene_set_overrepresentation",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        table=evidence,
    )
    report = make_summary_result(
        summary=summary,
        title=f"Gene-set overrepresentation summary: {provenance['comparison']}",
        operation="gene_set_overrepresentation",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
    )
    code = generic_ora_code(
        **runtime_parameters,
        expected_resource_sha256=diagnostics["resource"]["sha256"],
    )
    return result, report, code


@register_operation("openbio.node.genesetoverrepresentation")
def gene_set_overrepresentation_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(
        inputs,
        {"table", "universe", "gene_sets_file"},
        operation="Gene Set Overrepresentation",
    )
    require_parameters(
        parameters,
        {
            "resource_metadata_json",
            "comparison",
            "gene_column",
            "source_column",
            "target_column",
            "min_targets",
            "min_overlap",
            "max_p_adjusted",
            "max_output_rows",
        },
        operation="Gene Set Overrepresentation",
    )
    path, provenance = require_file_input(inputs, "gene_sets_file")
    return _table_records(
        context,
        *gene_set_overrepresentation_owned(
            _table_input(inputs, "table"),
            _table_input(inputs, "universe"),
            gene_sets_path=str(path),
            requested_resource_path=_requested_path(provenance, resource="Gene-set resource"),
            **parameters,
        ),
    )


def drug_scores_owned(
    adata: Any,
    resource: Any,
    *,
    drug: str = "",
    source: DynamicExpressionSource | None = None,
    output_key: str = "drug_target_score",
    min_matched_targets: int = 1,
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    expression = DRUG_EXPRESSION_SOURCE.resolve(adata, source)
    runtime_parameters = {
        "drug": drug,
        "source_kind": expression.kind,
        "layer_name": expression.layer_name,
        "output_key": output_key,
        "min_matched_targets": min_matched_targets,
        "overwrite_existing": overwrite_existing,
        "openbio_version": PLUGIN_VERSION,
    }
    output, summary = run_drug_score(
        adata,
        resource,
        copy_resource=False,
        **runtime_parameters,
        pertpy_module=None,
    )
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    finish_adata(output, "dgidb_drug_score", parameters, cells, genes, started_at, warnings=warnings)
    report = make_summary_result(
        summary=summary,
        title=f"DGIdb target-expression score: {drug}",
        operation="dgidb_drug_score",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    _, resource_metadata, resource_accounting, resource_artifact = validate_dgidb_resource(
        resource, copy_payload=False
    )
    code = drug_score_code(
        resource_metadata=resource_metadata,
        resource_accounting=resource_accounting,
        resource_artifact_metadata=resource_artifact,
        **runtime_parameters,
    )
    return output, report, code


@register_operation("openbio.node.drugscores")
def drug_scores_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "resource"}, operation="DGIdb Single-Drug Target Score")
    require_parameters(
        parameters,
        {"drug", "source", "output_key", "min_matched_targets", "overwrite_existing"},
        operation="DGIdb Single-Drug Target Score",
    )
    return _anndata_records(
        context,
        *drug_scores_owned(read_anndata_input(inputs), _dgidb_input(inputs), **parameters),
    )


def drug_hypergeometric_owned(
    table: Any,
    universe: Any,
    resource: Any,
    *,
    comparison: str = "",
    gene_column: str = "gene",
    min_targets: int = 3,
    min_overlap: int = 2,
    max_p_adjusted: float = 0.05,
    max_output_rows: int = 100_000,
) -> tuple[Any, Any, str]:
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    provenance = validate_enrichment_artifact_pair(
        table,
        universe,
        purpose="selected",
        selector=comparison,
        gene_column=gene_column,
        score_column=None,
        np=science.np,
        pd=science.pd,
    )
    runtime_parameters = {
        "artifact_family": provenance["artifact_family"],
        "comparison": provenance["comparison"],
        "comparison_column": provenance["comparison_column"],
        "gene_column": provenance["gene_column"],
        "evidence_scope": provenance["evidence_scope"],
        "inference_unit": provenance["inference_unit"],
        "direction": provenance["direction"],
        "replicate_aware": provenance["replicate_aware"],
        "technical_batch_handling": provenance["technical_batch_handling"],
        "upstream_parameters": provenance["upstream_parameters"],
        "provenance_warnings": provenance["provenance_warnings"],
        "provenance_declarations": provenance["provenance_declarations"],
        "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
        "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
        "expected_universe_fingerprint": provenance["universe_fingerprint"],
        "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
        "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
        "min_targets": min_targets,
        "min_overlap": min_overlap,
        "max_p_adjusted": max_p_adjusted,
        "max_output_rows": max_output_rows,
        "openbio_version": PLUGIN_VERSION,
    }
    evidence, summary = run_drug_ora(
        table, universe, resource, copy_resource=False, **runtime_parameters
    )
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    result = make_table_result(
        title=f"Drug target overrepresentation: {provenance['comparison']}",
        operation="drug_hypergeometric",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        table=evidence,
    )
    report = make_summary_result(
        summary=summary,
        title=f"Drug target overrepresentation summary: {provenance['comparison']}",
        operation="drug_hypergeometric",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
    )
    _, resource_metadata, resource_accounting, resource_artifact = validate_dgidb_resource(
        resource, copy_payload=False
    )
    code = drug_ora_code(
        resource_metadata=resource_metadata,
        resource_accounting=resource_accounting,
        resource_artifact_metadata=resource_artifact,
        **runtime_parameters,
    )
    return result, report, code


@register_operation("openbio.node.drughypergeometric")
def drug_hypergeometric_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(
        inputs,
        {"table", "universe", "resource"},
        operation="Drug Hypergeometric Enrichment",
    )
    require_parameters(
        parameters,
        {"comparison", "gene_column", "min_targets", "min_overlap", "max_p_adjusted", "max_output_rows"},
        operation="Drug Hypergeometric Enrichment",
    )
    return _table_records(
        context,
        *drug_hypergeometric_owned(
            _table_input(inputs, "table"),
            _table_input(inputs, "universe"),
            _dgidb_input(inputs),
            **parameters,
        ),
    )


def drug_gsea_owned(
    table: Any,
    universe: Any,
    resource: Any,
    *,
    comparison: str = "",
    gene_column: str = "gene",
    score_column: str = "score",
    min_targets: int = 15,
    max_targets: int = 500,
    n_permutations: int = 1000,
    random_seed: int = 123,
    max_output_rows: int = 100_000,
) -> tuple[Any, Any, str]:
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    provenance = validate_enrichment_artifact_pair(
        table,
        universe,
        purpose="ranked",
        selector=comparison,
        gene_column=gene_column,
        score_column=score_column,
        np=science.np,
        pd=science.pd,
    )
    runtime_parameters = {
        "artifact_family": provenance["artifact_family"],
        "comparison": provenance["comparison"],
        "comparison_column": provenance["comparison_column"],
        "gene_column": provenance["gene_column"],
        "score_column": provenance["score_column"],
        "evidence_scope": provenance["evidence_scope"],
        "inference_unit": provenance["inference_unit"],
        "direction": provenance["direction"],
        "replicate_aware": provenance["replicate_aware"],
        "technical_batch_handling": provenance["technical_batch_handling"],
        "upstream_parameters": provenance["upstream_parameters"],
        "provenance_warnings": provenance["provenance_warnings"],
        "provenance_declarations": provenance["provenance_declarations"],
        "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
        "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
        "expected_universe_fingerprint": provenance["universe_fingerprint"],
        "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
        "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
        "min_targets": min_targets,
        "max_targets": max_targets,
        "n_permutations": n_permutations,
        "random_seed": random_seed,
        "max_output_rows": max_output_rows,
        "openbio_version": PLUGIN_VERSION,
    }
    evidence, summary = run_drug_gsea(
        table, universe, resource, copy_resource=False, **runtime_parameters
    )
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    result = make_table_result(
        title=f"Drug GSEA: {provenance['comparison']}",
        operation="drug_gsea",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        random_seed=random_seed,
        table=evidence,
    )
    report = make_summary_result(
        summary=summary,
        title=f"Drug GSEA summary: {provenance['comparison']}",
        operation="drug_gsea",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    _, resource_metadata, resource_accounting, resource_artifact = validate_dgidb_resource(
        resource, copy_payload=False
    )
    code = drug_gsea_code(
        resource_metadata=resource_metadata,
        resource_accounting=resource_accounting,
        resource_artifact_metadata=resource_artifact,
        **runtime_parameters,
    )
    return result, report, code


@register_operation("openbio.node.druggsea")
def drug_gsea_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"table", "universe", "resource"}, operation="Drug GSEA")
    require_parameters(
        parameters,
        {
            "comparison",
            "gene_column",
            "score_column",
            "min_targets",
            "max_targets",
            "n_permutations",
            "random_seed",
            "max_output_rows",
        },
        operation="Drug GSEA",
    )
    return _table_records(
        context,
        *drug_gsea_owned(
            _table_input(inputs, "table"),
            _table_input(inputs, "universe"),
            _dgidb_input(inputs),
            **parameters,
        ),
    )


__all__ = [
    "AUCELL_EXPRESSION_SOURCE",
    "DGIDB_KIND",
    "DRUG_EXPRESSION_SOURCE",
    "GENE_PANEL_EXPRESSION_SOURCE",
    "GSVA_EXPRESSION_SOURCE",
    "TABLE_CODEC",
    "TABLE_KIND",
    "aucell_scores_operation",
    "aucell_scores_owned",
    "drug_gsea_operation",
    "drug_gsea_owned",
    "drug_hypergeometric_operation",
    "drug_hypergeometric_owned",
    "drug_scores_operation",
    "drug_scores_owned",
    "gene_set_overrepresentation_operation",
    "gene_set_overrepresentation_owned",
    "gene_panel_scores_operation",
    "gene_panel_scores_owned",
    "gsva_scores_operation",
    "gsva_scores_owned",
    "load_dgidb_resource_operation",
    "load_dgidb_resource_owned",
    "pathway_score_ttest_operation",
    "pathway_score_ttest_owned",
    "ranked_gsea_operation",
    "ranked_gsea_owned",
]
