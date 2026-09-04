from __future__ import annotations

import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from . import PLUGIN_VERSION, dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata, make_plot_result, make_summary_result, make_table_result
from .annotation_core import (
    build_celltypist_summary,
    celltypist_annotation_code,
    map_cluster_annotations_code,
    run_celltypist_annotation,
    run_map_cluster_annotations,
)
from .annotation_plotting import _standalone_celltypist_diagnostics_plot, celltypist_diagnostics_plot_code
from .artifact_envelope import table_from_metadata
from .contracts import TableResult
from .expression_source import _CELLTYPIST_SPEC, DynamicExpressionSource
from .marker_evidence import validate_marker_artifact_pair
from .marker_ora_plotting import _standalone_marker_ora_evidence_plot, marker_ora_evidence_plot_code
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    read_table_input,
    require_file_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
    write_table_output,
)
from .ora_evidence import build_marker_ora_summary, marker_ora_evidence_code, run_marker_ora_evidence
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
CELLTYPIST_EXPRESSION_SOURCE = _CELLTYPIST_SPEC

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
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)
CELLTYPIST_REFERENCE = AnalysisReference(
    citation=(
        "Domínguez Conde C, Xu C, Jarvis LB, et al. Cross-tissue immune cell analysis reveals "
        "tissue-specific features in humans. Science. 2022;376:eabl5197."
    ),
    doi="10.1126/science.abl5197",
    url="https://doi.org/10.1126/science.abl5197",
    kind="method",
)
FISHER_REFERENCE = AnalysisReference(
    citation=(
        "Fisher RA. On the interpretation of chi-square from contingency tables, and the calculation of P. "
        "Journal of the Royal Statistical Society. 1922;85:87-94."
    ),
    doi="10.2307/2340521",
    url="https://doi.org/10.2307/2340521",
    kind="method",
)
JSON_REFERENCE = AnalysisReference(
    citation="Bray T. The JavaScript Object Notation (JSON) Data Interchange Format. RFC 8259.",
    doi="10.17487/RFC8259",
    url="https://www.rfc-editor.org/info/rfc8259/",
    kind="software_documentation",
)


def _table_input(inputs: dict[str, JSONValue], name: str) -> Any:
    table, metadata = read_table_input(inputs, name, kind=TABLE_KIND)
    return table_from_metadata(metadata, table)


def _result_records(context: OperationContext, output: Any, report: Any, code: str) -> list[JSONValue]:
    return analysis_outputs(report, code, write_anndata_output(context, output))


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


def celltypist_diagnostics_plot_owned(
    adata: Any,
    *,
    metadata_key: str = "celltypist",
    view: dict[str, object] | None = None,
) -> tuple[Any, Any, str]:
    dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    png, details = _standalone_celltypist_diagnostics_plot(
        adata,
        metadata_key=metadata_key,
        view=view,
        _return_details=True,
    )
    parameters = {
        "metadata_key": details["metadata_key"],
        "view": details["view_parameters"],
    }
    warnings = list(details["warnings"])
    code = celltypist_diagnostics_plot_code(
        metadata_key=details["metadata_key"],
        view=details["view_parameters"],
    )
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="celltypist_diagnostics_plot",
        parameters=parameters,
        description=(
            f"Read-only {details['view'].replace('_', ' ')} over stored CellTypist class scores for "
            f"{details['plotted_cells']:,} cells."
        ),
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCellTypistDiagnosticsPlot",
        title=details["title"],
        operation="celltypist_diagnostics_plot",
        methods=(
            "Validated schema-version-1 CellTypist producer provenance, the observation and model-class axes, "
            "stored label/confidence columns, probability content fingerprints, and OpenBio annotation history. "
            "Rendered only the selected diagnostic view from stored sigmoid scores without loading a model or "
            "rerunning classification."
        ),
        results=(
            f"Displayed {details['view'].replace('_', ' ')} for {details['plotted_cells']:,} cells from "
            f"{details['model_class_count']:,} model classes. These remain CellTypist Provisional annotation "
            "diagnostics and are not a Sample-level Condition contrast or Curated annotation."
        ),
        key_results=details,
        parameters=parameters,
        references=[CELLTYPIST_REFERENCE, MATPLOTLIB_REFERENCE, ANNDATA_REFERENCE],
        software_packages=["anndata", "matplotlib", "numpy", "pandas"],
        warnings=warnings,
        limitations=[
            "CellTypist outputs are Provisional annotation and require marker, tissue, and study-context review "
            "before Curated annotation.",
            "Independent sigmoid class scores are not calibrated multiclass probabilities or biological certainty.",
            "Cell-level score distributions are descriptive and are not a replicate-aware Condition contrast.",
            "The probability heatmap averages within the explicit grouping column and does not infer lineage, "
            "identity, or statistical significance.",
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


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


def marker_ora_evidence_plot_owned(
    table: Any,
    *,
    view: dict[str, object] | None = None,
) -> tuple[Any, Any, str]:
    if not isinstance(table, TableResult):
        raise TypeError("Marker ORA Evidence Plot requires a TableResult input.")
    if not isinstance(table.source, Mapping) or table.source.get("operation") != "marker_ora_evidence":
        raise ValueError("Marker ORA Evidence Plot requires its exact Marker ORA Evidence producer.")
    source_parameters = table.source.get("parameters")
    if not isinstance(source_parameters, Mapping) or dict(source_parameters) != table.parameters:
        raise ValueError("Marker ORA Evidence Plot producer parameters are inconsistent.")
    producer = {
        "operation": table.source["operation"],
        "parameters": dict(table.parameters),
        "input_cells": int(table.input_cells),
        "input_genes": int(table.input_genes),
    }
    started_at = time.perf_counter()
    png, details = _standalone_marker_ora_evidence_plot(
        table.table,
        producer=producer,
        view=view,
        _return_details=True,
    )
    parameters = {"view": details["view_parameters"]}
    warnings = []
    if details["plotted_rows"] < details["available_rows"]:
        warnings.append(
            f"Displayed {details['plotted_rows']:,} of {details['available_rows']:,} stored rows by taking "
            "the first max_terms_per_group rows in upstream rank order."
        )
    code = marker_ora_evidence_plot_code(producer=producer, view=details["view_parameters"])
    title = "Marker ORA evidence"
    plotted = make_plot_result(
        png=png,
        title=title,
        operation="marker_ora_evidence_plot",
        parameters=parameters,
        description=(
            f"Read-only {details['view'].replace('_', ' ')} for {details['plotted_rows']:,} stored "
            "group-source evidence rows."
        ),
        warnings=warnings,
        input_cells=int(table.input_cells),
        input_genes=int(table.input_genes),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMarkerORAEvidencePlot",
        title=title,
        operation="marker_ora_evidence_plot",
        methods=(
            "Validated the exact canonical Marker ORA table, producer operation and policy, upstream marker "
            "fingerprints, contingency identities, Haldane–Anscombe log odds, overlap-gene JSON, descriptive "
            "gates, within-group order, and current-content fingerprint. Rendered only the stored evidence; "
            "no file was opened and no hypothesis was retested."
        ),
        results=(
            f"Displayed {details['plotted_rows']:,} stored rows across {len(details['groups']):,} groups; "
            f"{details['eligible_rows']:,} displayed rows met the producer's descriptive gates. This remains "
            "exploratory Cluster marker evidence and does not assign a Curated annotation or establish a "
            "Condition contrast."
        ),
        key_results=details,
        parameters=parameters,
        references=[FISHER_REFERENCE, MATPLOTLIB_REFERENCE, ANNOTATION_PRACTICE_REFERENCE],
        software_packages=["matplotlib", "numpy", "pandas"],
        warnings=warnings,
        limitations=[
            "Marker ORA is exploratory Cluster marker evidence, not a Curated annotation or cell-type truth.",
            "Marker selection and over-representation evidence reuse the same cells and are not independent "
            "confirmatory tests or a Sample-level Condition contrast.",
            "Odds ratios, overlaps, and adjusted p-values depend on the upstream marker filter, tested-gene "
            "universe, set collection, and descriptive thresholds.",
            "Truncation for display preserves upstream rank order and does not change the stored evidence.",
        ],
        input_cells=int(table.input_cells),
        input_genes=int(table.input_genes),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


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
    return _result_records(context, *celltypist_owned(read_anndata_input(inputs), **parameters))


@register_operation("openbio.node.celltypistdiagnosticsplot")
def celltypist_diagnostics_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="CellTypist Diagnostics Plot")
    require_parameters(parameters, {"metadata_key", "view"}, operation="CellTypist Diagnostics Plot")
    plotted, report, code = celltypist_diagnostics_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


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
    return analysis_outputs(report, code, write_table_output(context, result, kind=TABLE_KIND))


@register_operation("openbio.node.markeroraevidenceplot")
def marker_ora_evidence_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"table"}, operation="Marker ORA Evidence Plot")
    require_parameters(parameters, {"view"}, operation="Marker ORA Evidence Plot")
    plotted, report, code = marker_ora_evidence_plot_owned(_table_input(inputs, "table"), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


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
    return _result_records(context, *map_cluster_annotations_owned(read_anndata_input(inputs), **parameters))


__all__ = [
    "CELLTYPIST_EXPRESSION_SOURCE",
    "celltypist_annotation",
    "celltypist_diagnostics_plot",
    "celltypist_diagnostics_plot_owned",
    "celltypist_owned",
    "map_cluster_annotations",
    "map_cluster_annotations_owned",
    "marker_ora_evidence",
    "marker_ora_evidence_plot",
    "marker_ora_evidence_plot_owned",
    "marker_ora_owned",
]
