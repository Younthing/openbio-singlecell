from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from . import PLUGIN_VERSION, dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .annotation_core import (
    build_celltypist_summary,
    celltypist_annotation_code,
    map_cluster_annotations_code,
    run_celltypist_annotation,
    run_map_cluster_annotations,
)
from .artifact_codecs import read_anndata, read_table, write_table
from .artifact_envelope import result_metadata, table_from_metadata
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .marker_evidence import validate_marker_artifact_pair
from .operations_input import (
    ANNDATA_CODEC,
    ANNDATA_KIND,
    require_artifact_input,
    require_file_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
)
from .ora_evidence import build_marker_ora_summary, marker_ora_evidence_code, run_marker_ora_evidence
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
TABLE_CODEC = "table-jsonl-v1"

CELLTYPIST_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="CellTypist expression source",
    default="X",
    include_raw=True,
    layer_input_id="layer_name",
    layer_default="counts",
)

ANNOTATION_PRACTICE_REFERENCE = AnalysisReference(
    citation=(
        "Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. "
        "Molecular Systems Biology. 2019;15:e8746."
    ),
    doi="10.15252/msb.20188746",
    url="https://doi.org/10.15252/msb.20188746",
    kind="practice",
)
PANDAS_MAPPING_REFERENCE = AnalysisReference(
    citation="pandas developers. pandas.Series.map and categorical data documentation.",
    url="https://pandas.pydata.org/docs/reference/api/pandas.Series.map.html",
    kind="software_documentation",
)
PANDAS_REFERENCE = AnalysisReference(
    citation="McKinney W. Data Structures for Statistical Computing in Python. SciPy 2010.",
    doi="10.25080/Majora-92bf1922-00a",
    url="https://doi.org/10.25080/Majora-92bf1922-00a",
    kind="software",
)
ANNDATA_REFERENCE = AnalysisReference(
    citation="Virshup I, et al. anndata: Annotated data. Journal of Open Source Software. 2024;9:4371.",
    doi="10.21105/joss.04371",
    url="https://doi.org/10.21105/joss.04371",
    kind="software",
)
JSON_REFERENCE = AnalysisReference(
    citation="Bray T. The JavaScript Object Notation (JSON) Data Interchange Format. RFC 8259.",
    doi="10.17487/RFC8259",
    url="https://www.rfc-editor.org/info/rfc8259/",
    kind="software_documentation",
)


def _anndata_input(inputs: dict[str, JSONValue], name: str = "adata") -> Any:
    root = require_artifact_input(inputs, name, kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    return read_anndata(root)


def _table_input(inputs: dict[str, JSONValue], name: str) -> Any:
    root = require_artifact_input(inputs, name, kind=TABLE_KIND, codec=TABLE_CODEC)
    table, metadata = read_table(root)
    return table_from_metadata(metadata, table)


def _table_output(context: OperationContext, result: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("table")
    write_table(root, result.table, result_metadata(result))
    return {
        "type": "artifact",
        "name": "table",
        "kind": TABLE_KIND,
        "codec": TABLE_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _result_records(context: OperationContext, output: Any, report: Any, code: str) -> list[JSONValue]:
    return [
        write_anndata_output(context, output),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


def celltypist_owned(
    adata: Any,
    *,
    source: DynamicExpressionSource | None = None,
    expression_state: str = "verified_cp10k_log1p",
    model: str = "",
    majority_voting: bool = False,
    over_clustering_key: str = "",
    min_prop: float = 0.0,
    label_column: str = "celltypist_cell_type",
    confidence_column: str = "celltypist_confidence",
    probability_key: str = "celltypist_probabilities",
    store_decision_matrix: bool = False,
    decision_key: str = "celltypist_decision_scores",
    metadata_key: str = "celltypist",
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    dependencies.require_scientific_dependencies()
    expression = CELLTYPIST_EXPRESSION_SOURCE.resolve(adata, source)
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, diagnostics = run_celltypist_annotation(
        adata,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        expression_state=expression_state,
        model_identifier=model,
        majority_voting=majority_voting,
        over_clustering_key=over_clustering_key,
        min_prop=min_prop,
        label_column=label_column,
        confidence_column=confidence_column,
        probability_key=probability_key,
        store_decision_matrix=store_decision_matrix,
        decision_key=decision_key,
        metadata_key=metadata_key,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
        celltypist_module=None,
    )
    parameters = dict(diagnostics["parameters"])
    warnings = list(diagnostics["warnings"])
    finish_adata(output, "celltypist_annotation", parameters, cells, genes, started_at, warnings=warnings)
    code = celltypist_annotation_code(
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        expression_state=expression_state,
        resolved_model_path=diagnostics["model"]["resolved_path"],
        requested_model_identifier=diagnostics["model"]["requested_identifier"],
        expected_model_sha256=diagnostics["model"]["sha256"],
        majority_voting=majority_voting,
        over_clustering_key=over_clustering_key.strip(),
        min_prop=float(min_prop),
        label_column=label_column.strip(),
        confidence_column=confidence_column.strip(),
        probability_key=probability_key.strip(),
        store_decision_matrix=store_decision_matrix,
        decision_key=decision_key.strip(),
        metadata_key=metadata_key.strip(),
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
    )
    summary = build_celltypist_summary(diagnostics)
    report = make_summary_result(
        summary=summary,
        title="CellTypist provisional annotation summary",
        operation="celltypist_annotation",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    return output, report, code


def marker_ora_owned(
    table: Any,
    universe: Any,
    *,
    resource_path: str | Path,
    requested_resource_path: str,
    resource_metadata_json: str = "{}",
    source_column: str = "source",
    target_column: str = "target",
    min_targets: int = 3,
    min_overlap: int = 2,
    max_p_adjusted: float = 0.05,
) -> tuple[Any, Any, str]:
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    provenance = validate_marker_artifact_pair(
        table,
        universe,
        np=science.np,
        pd=science.pd,
        allowed_table_operations=("filter_marker_genes",),
    )
    if provenance["table_operation"] != "filter_marker_genes":
        raise RuntimeError("Marker ORA Evidence accepted a non-filtered marker artifact unexpectedly.")
    evidence, diagnostics = run_marker_ora_evidence(
        table,
        universe,
        expected_analysis_fingerprint=provenance["analysis_fingerprint"],
        expected_ranking_fingerprint=provenance["ranking_fingerprint"],
        expected_universe_fingerprint=provenance["universe_fingerprint"],
        expected_marker_content_fingerprint=provenance["table_content_fingerprint"],
        expected_universe_content_fingerprint=provenance["universe_content_fingerprint"],
        group_labels=provenance["group_labels"],
        upstream_parameters=provenance["upstream_parameters"],
        resource_path=resource_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        source_column=source_column,
        target_column=target_column,
        min_targets=min_targets,
        min_overlap=min_overlap,
        max_p_adjusted=max_p_adjusted,
        openbio_version=PLUGIN_VERSION,
        decoupler_module=None,
    )
    summary = build_marker_ora_summary(diagnostics)
    parameters = dict(summary["parameters"])
    warnings = list(summary["warnings"])
    result = make_table_result(
        table=evidence,
        title="Marker ORA evidence",
        operation="marker_ora_evidence",
        parameters=parameters,
        description=(
            "Explicit selected-set over-representation evidence for Cluster marker review; "
            "no AnnData labels were assigned."
        ),
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
    )
    code = marker_ora_evidence_code(
        expected_analysis_fingerprint=provenance["analysis_fingerprint"],
        expected_ranking_fingerprint=provenance["ranking_fingerprint"],
        expected_universe_fingerprint=provenance["universe_fingerprint"],
        expected_marker_content_fingerprint=provenance["table_content_fingerprint"],
        expected_universe_content_fingerprint=provenance["universe_content_fingerprint"],
        group_labels=provenance["group_labels"],
        upstream_parameters=provenance["upstream_parameters"],
        resource_path=resource_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        expected_resource_sha256=diagnostics["resource"]["sha256"],
        source_column=source_column.strip(),
        target_column=target_column.strip(),
        min_targets=min_targets,
        min_overlap=min_overlap,
        max_p_adjusted=float(max_p_adjusted),
        openbio_version=PLUGIN_VERSION,
    )
    report = make_summary_result(
        summary=summary,
        title="Marker ORA evidence summary",
        operation="marker_ora_evidence",
        parameters=parameters,
        description=summary["results"],
        warnings=warnings,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
    )
    return result, report, code


def map_cluster_annotations_owned(
    adata: Any,
    *,
    groupby: str = "leiden",
    mapping_json: str = "{}",
    output_column: str = "cell_type",
    unmapped_policy: str = "error",
    annotation_status: str = "provisional",
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, diagnostics = run_map_cluster_annotations(
        adata,
        groupby=groupby,
        mapping_json=mapping_json,
        output_column=output_column,
        unmapped_policy=unmapped_policy,
        annotation_status=annotation_status,
        overwrite_existing=overwrite_existing,
    )
    normalized_mapping = {**diagnostics["effective_mapping"], **diagnostics["unused_declared_level_mapping"]}
    parameters = {
        "groupby": groupby.strip(),
        "mapping": normalized_mapping,
        "output_column": output_column.strip(),
        "unmapped_policy": unmapped_policy,
        "annotation_status": annotation_status,
        "overwrite_existing": overwrite_existing,
    }
    warnings = []
    if diagnostics["missing_source_cells"]:
        warnings.append(
            f"{diagnostics['missing_source_cells']:,} cells had missing source labels; "
            "true missing values were preserved."
        )
    if unmapped_policy == "preserve_cluster_label" and diagnostics["unmapped_cells"]:
        warnings.append("The output mixes mapped biological labels with preserved raw cluster identifiers.")
    if annotation_status == "curated":
        warnings.append(
            "Curated status is a caller declaration; the software did not verify expert review, marker evidence, "
            "or ontology validity."
        )
    finish_adata(output, "map_cluster_annotations", parameters, cells, genes, started_at, warnings=warnings)
    code = map_cluster_annotations_code(
        groupby=groupby.strip(),
        mapping_json=mapping_json,
        output_column=output_column.strip(),
        unmapped_policy=unmapped_policy,
        annotation_status=annotation_status,
        overwrite_existing=overwrite_existing,
    )
    output_counts_text = ", ".join(
        f"{category!r}={count:,}" for category, count in diagnostics["output_counts"].items()
    ) or "none"
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMapClusterAnnotations",
        title="Cluster annotation mapping summary",
        operation="map_cluster_annotations",
        methods=(
            f"A strict JSON map was applied to observed {groupby.strip()!r} levels. Levels absent from the map "
            f"were handled by policy {unmapped_policy!r}; output categories followed the declared categorical "
            "source order or lexical observed-level order."
        ),
        results=(
            f"Mapped {len(diagnostics['effective_mapping']):,} of {diagnostics['observed_level_count']:,} observed "
            f"clusters, covering {diagnostics['mapped_cells']:,} of {cells:,} cells. Policy {unmapped_policy!r} "
            f"preserved {diagnostics['preserved_cells']:,} cells and set {diagnostics['set_missing_cells']:,} cells "
            f"missing; {diagnostics['missing_source_cells']:,} cells were already missing in the source and "
            f"{diagnostics['output_missing_cells']:,} cells are missing in the output. The output has "
            f"{len(diagnostics['output_categories']):,} categories with counts {output_counts_text}, and the map "
            f"contains {len(diagnostics['many_to_one_merges']):,} many-to-one target merges. The output column "
            f"existed before execution: {'yes' if diagnostics['output_existed'] else 'no'}; overwrite_existing was "
            f"{'true' if overwrite_existing else 'false'}. The result is recorded as "
            f"{'caller-declared Curated annotation' if annotation_status == 'curated' else 'Provisional annotation'}."
        ),
        key_results={
            **diagnostics,
            "annotation_provenance": output.uns["openbio_singlecell"]["annotations"][output_column.strip()],
            "mapping_unknown_keys": [],
            "output_dtype": "category",
            "output_ordered": False,
        },
        parameters=parameters,
        references=(
            ANNOTATION_PRACTICE_REFERENCE,
            PANDAS_MAPPING_REFERENCE,
            PANDAS_REFERENCE,
            ANNDATA_REFERENCE,
            JSON_REFERENCE,
        ),
        software_packages=("pandas", "anndata"),
        warnings=warnings,
        limitations=(
            "Every cell in one source cluster receives the same label, so within-cluster heterogeneity is hidden.",
            "Cluster resolution and preprocessing affect the source partition and therefore the applied annotation.",
            "The mapping operation does not inspect Cluster marker evidence, reference expression, or ontology identifiers.",
            "Provisional annotation requires review before use as Curated annotation in Sample-level Condition contrast.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=code,
    )
    return output, report, code


@register_operation("openbio.node.celltypistannotation")
def celltypist_annotation(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="CellTypist Annotation")
    expected = {
        "source",
        "expression_state",
        "model",
        "majority_voting",
        "over_clustering_key",
        "min_prop",
        "label_column",
        "confidence_column",
        "probability_key",
        "store_decision_matrix",
        "decision_key",
        "metadata_key",
        "overwrite_existing",
    }
    require_parameters(parameters, expected, operation="CellTypist Annotation")
    return _result_records(context, *celltypist_owned(_anndata_input(inputs), **parameters))


@register_operation("openbio.node.markeroraevidence")
def marker_ora_evidence(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"table", "universe", "resource_csv"}, operation="Marker ORA Evidence")
    expected = {
        "resource_metadata_json",
        "source_column",
        "target_column",
        "min_targets",
        "min_overlap",
        "max_p_adjusted",
    }
    require_parameters(parameters, expected, operation="Marker ORA Evidence")
    resource_path, provenance = require_file_input(inputs, "resource_csv")
    requested = provenance.get("path")
    if not isinstance(requested, str):
        raise ProtocolError("Marker ORA resource provenance path must be a string.")
    result, report, code = marker_ora_owned(
        _table_input(inputs, "table"),
        _table_input(inputs, "universe"),
        resource_path=resource_path,
        requested_resource_path=requested,
        **parameters,
    )
    return [
        _table_output(context, result),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.mapclusterannotations")
def map_cluster_annotations(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Map Cluster Annotations")
    expected = {
        "groupby",
        "mapping_json",
        "output_column",
        "unmapped_policy",
        "annotation_status",
        "overwrite_existing",
    }
    require_parameters(parameters, expected, operation="Map Cluster Annotations")
    return _result_records(context, *map_cluster_annotations_owned(_anndata_input(inputs), **parameters))


__all__ = [
    "CELLTYPIST_EXPRESSION_SOURCE",
    "celltypist_annotation",
    "celltypist_owned",
    "map_cluster_annotations",
    "map_cluster_annotations_owned",
    "marker_ora_evidence",
    "marker_ora_owned",
]
