from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from . import PLUGIN_VERSION
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import make_plot_result, make_summary_result, make_table_result
from .artifact_envelope import table_from_metadata
from .collectri_ulm import collectri_ulm_code, run_collectri_ulm
from .expression_source import _COLLECTRI_SPEC, DynamicExpressionSource
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
from .pyscenic_import import import_pyscenic_bundle, pyscenic_import_code
from .regulatory_artifact_codecs import (
    SCENIC_CODEC,
    TF_ACTIVITY_CODEC,
    read_scenic,
    read_tf_activity,
    write_scenic,
    write_tf_activity,
)
from .regulatory_plotting import (
    PLOT_EVIDENCE_SCHEMA,
    regulatory_table_plot_code,
    run_regulatory_table_plot,
    run_scenic_activity_plot,
    run_tf_activity_plot,
    scenic_activity_plot_code,
    table_evidence_fingerprint,
    tf_activity_plot_code,
)
from .scenic_binarization import binarize_scenic_activity, scenic_binarization_code
from .scenic_membership import scenic_membership_code, scenic_regulon_membership
from .scenic_rss import compute_scenic_rss, scenic_rss_code
from .tf_activity_ranking import rank_tf_activities, rank_tf_activities_code
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
TF_ACTIVITY_KIND = "OPENBIO_TF_ACTIVITY"
SCENIC_KIND = "OPENBIO_SCENIC_RESULT"
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
COLLECTRI_EXPRESSION_SOURCE = _COLLECTRI_SPEC

MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D graphics environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)


def _plot_evidence_parameters(parameters: dict[str, Any], table: Any, *, contract: str) -> dict[str, Any]:
    return {
        **parameters,
        "plot_evidence_schema": PLOT_EVIDENCE_SCHEMA,
        "plot_evidence_fingerprint_sha256": table_evidence_fingerprint(table, contract=contract),
    }


def _tf_activity_input(inputs: dict[str, JSONValue]) -> dict[str, Any]:
    root = require_artifact_input(
        inputs,
        "activities",
        kind=TF_ACTIVITY_KIND,
        codec=TF_ACTIVITY_CODEC,
    )
    return read_tf_activity(root)


def _scenic_input(inputs: dict[str, JSONValue]) -> dict[str, Any]:
    root = require_artifact_input(
        inputs,
        "scenic_result",
        kind=SCENIC_KIND,
        codec=SCENIC_CODEC,
    )
    return read_scenic(root)


def _table_input(inputs: dict[str, JSONValue], name: str) -> Any:
    table, metadata = read_table_input(inputs, name, kind=TABLE_KIND)
    return table_from_metadata(metadata, table)


def _tf_activity_record(context: OperationContext, artifact: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("activities")
    write_tf_activity(root, artifact)
    return {
        "type": "artifact",
        "name": "activities",
        "kind": TF_ACTIVITY_KIND,
        "codec": TF_ACTIVITY_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _scenic_record(context: OperationContext, artifact: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("scenic_result")
    write_scenic(root, artifact)
    return {
        "type": "artifact",
        "name": "scenic_result",
        "kind": SCENIC_KIND,
        "codec": SCENIC_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def collectri_ulm_owned(
    adata: Any,
    *,
    resource_mode: str = "local_network",
    network_path: str | Path | None = None,
    resource_metadata_json: str = "{}",
    source: DynamicExpressionSource | None = None,
    affiliation_license: str = "academic",
    allow_network_access: bool = False,
    complex_policy: str = "retain",
    min_targets: int = 5,
    batch_size: int = 250_000,
    max_output_rows: int = 2_000_000,
    max_working_memory_gib: float = 4.0,
    overwrite_existing: bool = False,
    decoupler_module: Any | None = None,
) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    expression = COLLECTRI_EXPRESSION_SOURCE.resolve(adata, source)
    resolved_path = None if network_path is None else str(Path(network_path).resolve())
    output, activities, summary = run_collectri_ulm(
        adata,
        resource_mode=resource_mode,
        network_path=resolved_path,
        resource_metadata_json=resource_metadata_json,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        affiliation_license=affiliation_license,
        allow_network_access=allow_network_access,
        complex_policy=complex_policy,
        min_targets=min_targets,
        batch_size=batch_size,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        overwrite_existing=overwrite_existing,
        decoupler_module=decoupler_module,
        openbio_version=PLUGIN_VERSION,
    )
    report = make_summary_result(
        summary=summary,
        title="CollecTRI ULM activity summary",
        operation="collectri_ulm",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    code = collectri_ulm_code(
        parameters={
            "resource_mode": resource_mode,
            "network_path": resolved_path,
            "resource_metadata_json": resource_metadata_json,
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "affiliation_license": affiliation_license,
            "allow_network_access": allow_network_access,
            "complex_policy": complex_policy,
            "min_targets": min_targets,
            "batch_size": batch_size,
            "max_output_rows": max_output_rows,
            "max_working_memory_gib": max_working_memory_gib,
            "overwrite_existing": overwrite_existing,
            "expected_resource_sha256": summary["key_results"]["resource"]["file_sha256"],
            "expected_canonical_network_sha256": summary["key_results"]["resource"][
                "canonical_network_sha256"
            ],
            "openbio_version": PLUGIN_VERSION,
        },
        resolved_resource_provenance=summary["key_results"]["resource"],
    )
    return output, activities, report, code


def rank_tf_activities_owned(
    adata: Any,
    activities: Any,
    *,
    annotation_key: str = "cell_type",
    annotation_status: str = "unknown",
    reference: str = "rest",
    method: str = "t-test_overestim_var",
    report_p_adjusted: float = 0.05,
    max_output_rows: int = 100_000,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    table, summary = rank_tf_activities(
        adata,
        activities,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        reference=reference,
        method=method,
        report_p_adjusted=report_p_adjusted,
        max_output_rows=max_output_rows,
        openbio_version=PLUGIN_VERSION,
        _portable_artifact=True,
        _worker_owned=True,
    )
    common = {
        "parameters": summary["parameters"],
        "description": summary["results"],
        "warnings": summary["warnings"],
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "started_at": started_at,
    }
    result = make_table_result(
        table=table,
        title=f"TF activities by {annotation_key}",
        operation="rank_tf_activities",
        parameters=_plot_evidence_parameters(summary["parameters"], table, contract="tf_ranking"),
        **{key: value for key, value in common.items() if key != "parameters"},
    )
    report = make_summary_result(
        summary=summary,
        title=f"TF activity ranking by {annotation_key}",
        operation="rank_tf_activities",
        **common,
    )
    code = rank_tf_activities_code(
        parameters={
            "annotation_key": annotation_key,
            "annotation_status": annotation_status,
            "reference": reference,
            "method": method,
            "report_p_adjusted": report_p_adjusted,
            "max_output_rows": max_output_rows,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    return result, report, code


def import_pyscenic_owned(
    adata: Any,
    manifest_path: str | Path,
    *,
    overwrite: bool = False,
    max_file_bytes: int = 2_147_483_647,
    max_adjacency_edges: int = 10_000_000,
    max_regulon_edges: int = 2_000_000,
    max_dense_bytes: int = 1_073_741_824,
) -> tuple[Any, Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    resolved_manifest = Path(manifest_path).resolve()
    output, scenic_result, summary = import_pyscenic_bundle(
        adata,
        resolved_manifest,
        overwrite=overwrite,
        max_file_bytes=max_file_bytes,
        max_adjacency_edges=max_adjacency_edges,
        max_regulon_edges=max_regulon_edges,
        max_dense_bytes=max_dense_bytes,
        openbio_version=PLUGIN_VERSION,
        _portable_artifact=True,
    )
    report = make_summary_result(
        summary=summary,
        title="Imported pySCENIC regulatory evidence",
        operation="import_pyscenic_results",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    code = pyscenic_import_code(
        manifest_path=str(resolved_manifest),
        overwrite=overwrite,
        max_file_bytes=max_file_bytes,
        max_adjacency_edges=max_adjacency_edges,
        max_regulon_edges=max_regulon_edges,
        max_dense_bytes=max_dense_bytes,
    )
    return output, scenic_result, report, code


def scenic_rss_owned(
    adata: Any,
    scenic_result: Any,
    *,
    annotation_key: str = "cell_type",
    annotation_status: str = "unknown",
    max_output_rows: int = 100_000,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    table, summary = compute_scenic_rss(
        adata,
        scenic_result,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        max_output_rows=max_output_rows,
        openbio_version=PLUGIN_VERSION,
        _portable_artifact=True,
        _worker_owned=True,
    )
    common = {
        "parameters": summary["parameters"],
        "description": summary["results"],
        "warnings": summary["warnings"],
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "started_at": started_at,
    }
    result = make_table_result(
        table=table,
        title=f"SCENIC regulon specificity by {annotation_key}",
        operation="scenic_regulon_specificity",
        parameters=_plot_evidence_parameters(summary["parameters"], table, contract="scenic_rss"),
        **{key: value for key, value in common.items() if key != "parameters"},
    )
    report = make_summary_result(
        summary=summary,
        title=f"SCENIC regulon specificity by {annotation_key}",
        operation="scenic_regulon_specificity",
        **common,
    )
    code = scenic_rss_code(
        parameters={
            "annotation_key": annotation_key,
            "annotation_status": annotation_status,
            "max_output_rows": max_output_rows,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    return result, report, code


def scenic_binarization_owned(
    scenic_result: Any,
    *,
    random_seed: int = 1,
    threshold_overrides: Any | None = None,
    max_dense_bytes: int = 1_073_741_824,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    overrides = (
        threshold_overrides.table
        if hasattr(threshold_overrides, "table")
        else threshold_overrides
    )
    binary, threshold_table, summary = binarize_scenic_activity(
        scenic_result,
        random_seed=random_seed,
        threshold_overrides=overrides,
        max_dense_bytes=max_dense_bytes,
        openbio_version=PLUGIN_VERSION,
        _portable_artifact=True,
        _worker_owned=True,
    )
    del binary
    provenance = scenic_result["provenance"]
    dimensions = provenance["input_dimensions"]
    common = {
        "parameters": summary["parameters"],
        "description": summary["results"],
        "warnings": summary["warnings"],
        "input_cells": int(dimensions["cells"]),
        "input_genes": int(dimensions["genes"]),
        "started_at": started_at,
        "random_seed": random_seed,
    }
    result = make_table_result(
        table=threshold_table,
        title="SCENIC activity thresholds",
        operation="scenic_activity_binarization",
        parameters=_plot_evidence_parameters(
            summary["parameters"], threshold_table, contract="scenic_binarization"
        ),
        **{key: value for key, value in common.items() if key != "parameters"},
    )
    report = make_summary_result(
        summary=summary,
        title="SCENIC activity binarization",
        operation="scenic_activity_binarization",
        **common,
    )
    code = scenic_binarization_code(
        parameters={
            "random_seed": random_seed,
            "max_dense_bytes": max_dense_bytes,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    return result, report, code


def scenic_membership_owned(
    scenic_result: Any,
    *,
    transcription_factor: str = "",
    max_output_rows: int = 100_000,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    table, summary = scenic_regulon_membership(
        scenic_result,
        transcription_factor=transcription_factor,
        max_output_rows=max_output_rows,
        openbio_version=PLUGIN_VERSION,
        _portable_artifact=True,
    )
    dimensions = scenic_result["provenance"]["input_dimensions"]
    common = {
        "parameters": summary["parameters"],
        "description": summary["results"],
        "warnings": summary["warnings"],
        "input_cells": int(dimensions["cells"]),
        "input_genes": int(dimensions["genes"]),
        "started_at": started_at,
    }
    title = f"SCENIC final regulons: {transcription_factor or 'all TFs'}"
    result = make_table_result(
        table=table,
        title=title,
        operation="scenic_final_regulon_membership",
        parameters=_plot_evidence_parameters(summary["parameters"], table, contract="scenic_membership"),
        **{key: value for key, value in common.items() if key != "parameters"},
    )
    report = make_summary_result(
        summary=summary,
        title=title,
        operation="scenic_final_regulon_membership",
        **common,
    )
    code = scenic_membership_code(
        parameters={
            "transcription_factor": transcription_factor,
            "max_output_rows": max_output_rows,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    return result, report, code


@register_operation("openbio.node.collectriulm")
def collectri_ulm(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_parameters(
        parameters,
        {
            "resource_mode",
            "resource_metadata_json",
            "source",
            "affiliation_license",
            "allow_network_access",
            "complex_policy",
            "min_targets",
            "batch_size",
            "max_output_rows",
            "max_working_memory_gib",
            "overwrite_existing",
        },
        operation="CollecTRI ULM",
    )
    resource_mode = parameters["resource_mode"]
    if resource_mode == "local_network":
        require_input_names(inputs, {"adata", "network_csv"}, operation="CollecTRI ULM")
        network_path, _provenance = require_file_input(inputs, "network_csv")
    elif resource_mode == "official_collectri":
        require_input_names(inputs, {"adata"}, operation="CollecTRI ULM")
        network_path = None
    else:
        raise ProtocolError("CollecTRI resource_mode must select local_network or official_collectri.")
    output, activities, report, code = collectri_ulm_owned(
        read_anndata_input(inputs),
        network_path=network_path,
        **parameters,
    )
    return analysis_outputs(
        report,
        code,
        write_anndata_output(context, output),
        _tf_activity_record(context, activities),
    )


@register_operation("openbio.node.ranktfactivities")
def rank_tf_activity(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "activities"}, operation="Rank TF Activities")
    require_parameters(
        parameters,
        {
            "annotation_key",
            "annotation_status",
            "reference",
            "method",
            "report_p_adjusted",
            "max_output_rows",
        },
        operation="Rank TF Activities",
    )
    result, report, code = rank_tf_activities_owned(
        read_anndata_input(inputs),
        _tf_activity_input(inputs),
        **parameters,
    )
    return analysis_outputs(report, code, write_table_output(context, result, kind=TABLE_KIND))


@register_operation("openbio.node.importpyscenicresults")
def import_pyscenic_results(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "run_manifest_json"}, operation="Import pySCENIC Results")
    require_parameters(
        parameters,
        {
            "overwrite",
            "max_file_bytes",
            "max_adjacency_edges",
            "max_regulon_edges",
            "max_dense_bytes",
        },
        operation="Import pySCENIC Results",
    )
    manifest_path, _provenance = require_file_input(inputs, "run_manifest_json")
    output, scenic_result, report, code = import_pyscenic_owned(
        read_anndata_input(inputs),
        manifest_path,
        **parameters,
    )
    return analysis_outputs(
        report,
        code,
        write_anndata_output(context, output),
        _scenic_record(context, scenic_result),
    )


@register_operation("openbio.node.scenicregulonspecificity")
def scenic_regulon_specificity(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata", "scenic_result"}, operation="SCENIC Regulon Specificity")
    require_parameters(
        parameters,
        {"annotation_key", "annotation_status", "max_output_rows"},
        operation="SCENIC Regulon Specificity",
    )
    result, report, code = scenic_rss_owned(
        read_anndata_input(inputs),
        _scenic_input(inputs),
        **parameters,
    )
    return analysis_outputs(report, code, write_table_output(context, result, kind=TABLE_KIND))


@register_operation("openbio.node.scenicactivitybinarization")
def scenic_activity_binarization(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    if set(inputs) not in ({"scenic_result"}, {"scenic_result", "threshold_overrides"}):
        raise ProtocolError(
            "SCENIC Activity Binarization inputs must contain scenic_result and optional threshold_overrides."
        )
    require_parameters(
        parameters,
        {"random_seed", "max_dense_bytes"},
        operation="SCENIC Activity Binarization",
    )
    overrides = (
        _table_input(inputs, "threshold_overrides")
        if "threshold_overrides" in inputs
        else None
    )
    result, report, code = scenic_binarization_owned(
        _scenic_input(inputs),
        threshold_overrides=overrides,
        **parameters,
    )
    return analysis_outputs(
        report,
        code,
        write_table_output(context, result, name="thresholds", kind=TABLE_KIND),
    )


@register_operation("openbio.node.scenictfmodules")
def scenic_tf_modules(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"scenic_result"}, operation="SCENIC Final Regulon Membership")
    require_parameters(
        parameters,
        {"transcription_factor", "max_output_rows"},
        operation="SCENIC Final Regulon Membership",
    )
    result, report, code = scenic_membership_owned(_scenic_input(inputs), **parameters)
    return analysis_outputs(
        report,
        code,
        write_table_output(context, result, kind=TABLE_KIND, json_list_columns=("motif_ids",)),
    )


def _regulatory_activity_plot_owned(
    adata: Any,
    artifact: Any,
    *,
    node_id: str,
    operation: str,
    label: str,
    runner: Any,
    code_builder: Any,
    view: Any,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = runner(adata, artifact, view=view)
    parameters = {"view": details["view_parameters"]}
    code = code_builder(view=details["view_parameters"])
    description = (
        f"Rendered {details['plotted_activities']:,} of {details['available_activities']:,} stored {label} "
        f"activities across {details['plotted_cells']:,} cells using the {details['view'].replace('_', ' ')} view."
    )
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation=operation,
        parameters=parameters,
        description=description,
        warnings=[],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id=node_id,
        title=f"{label} activity plot summary",
        operation=operation,
        methods=(
            "Validated the exact typed producer artifact, current-content fingerprint, activity axes, and exact "
            "AnnData observation alignment before read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"title", "view_parameters"}},
        parameters=parameters,
        references=[MATPLOTLIB_REFERENCE],
        software_packages=["matplotlib", "numpy", "pandas"],
        warnings=[],
        limitations=(
            "Cell-level activity distributions and arithmetic group means are descriptive and are not Sample-level Condition inference.",
            "Activity scores inherit upstream expression, resource, cell-cohort, and Technical-batch limitations and do not establish TF binding or causality.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


def tf_activity_plot_owned(adata: Any, activities: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    return _regulatory_activity_plot_owned(
        adata,
        activities,
        node_id="OpenBioSingleCellTFActivityPlot",
        operation="tf_activity_plot",
        label="CollecTRI ULM",
        runner=run_tf_activity_plot,
        code_builder=tf_activity_plot_code,
        view=view,
    )


def scenic_activity_plot_owned(adata: Any, scenic_result: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    return _regulatory_activity_plot_owned(
        adata,
        scenic_result,
        node_id="OpenBioSingleCellSCENICActivityPlot",
        operation="scenic_activity_plot",
        label="SCENIC AUCell",
        runner=run_scenic_activity_plot,
        code_builder=scenic_activity_plot_code,
        view=view,
    )


def _regulatory_table_plot_owned(
    result: Any,
    *,
    contract: str,
    node_id: str,
    operation: str,
    function_name: str,
    label: str,
    view: Any,
) -> tuple[Any, Any, str]:
    from .contracts import TableResult

    if not isinstance(result, TableResult):
        raise TypeError(f"{label} Plot requires its canonical TableResult input.")
    source_parameters = result.source.get("parameters") if isinstance(result.source, dict) else None
    if not isinstance(source_parameters, dict) or any(
        result.parameters.get(key) != source_parameters.get(key)
        for key in ("plot_evidence_schema", "plot_evidence_fingerprint_sha256")
    ):
        raise ValueError(f"{label} producer parameters are internally inconsistent.")
    producer = {
        "operation": result.source.get("operation"),
        "parameters": result.parameters,
        "input_cells": result.input_cells,
        "input_genes": result.input_genes,
    }
    started_at = time.perf_counter()
    png, details = run_regulatory_table_plot(
        result.table,
        producer=producer,
        contract=contract,
        view=view,
    )
    parameters = {"view": details["view_parameters"]}
    code = regulatory_table_plot_code(
        function_name=function_name,
        producer=producer,
        contract=contract,
        view=details["view_parameters"],
    )
    description = (
        f"Rendered {details['plotted_marks']:,} visual marks from {details['total_rows']:,} stored {label} "
        "evidence rows using the "
        f"{details['view'].replace('_', ' ')} view."
    )
    warnings = list(result.warnings)
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation=operation,
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=result.input_cells,
        input_genes=result.input_genes,
        started_at=started_at,
        random_seed=result.random_seed,
    )
    report, code = make_analysis_report(
        node_id=node_id,
        title=f"{label} plot summary",
        operation=operation,
        methods=(
            "Validated exact producer identity, canonical schema, scientific table invariants, and the "
            "current-content fingerprint before read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"title", "view_parameters"}},
        parameters=parameters,
        references=[MATPLOTLIB_REFERENCE],
        software_packages=["matplotlib", "numpy", "pandas"],
        warnings=warnings,
        limitations=(
            "The plot presents upstream evidence and does not rerun activity inference, statistical testing, RSS, or binarization.",
            "Regulon and TF evidence does not establish direct binding, causality, cell identity, or a replicate-aware Condition effect.",
        ),
        input_cells=result.input_cells,
        input_genes=result.input_genes,
        started_at=started_at,
        code=code,
        random_seed=result.random_seed,
    )
    return plotted, report, code


def tf_activity_ranking_plot_owned(result: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    return _regulatory_table_plot_owned(
        result,
        contract="tf_ranking",
        node_id="OpenBioSingleCellTFActivityRankingPlot",
        operation="tf_activity_ranking_plot",
        function_name="plot_tf_activity_ranking",
        label="TF activity ranking",
        view=view,
    )


def scenic_rss_plot_owned(result: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    return _regulatory_table_plot_owned(
        result,
        contract="scenic_rss",
        node_id="OpenBioSingleCellSCENICRegulonSpecificityPlot",
        operation="scenic_regulon_specificity_plot",
        function_name="plot_scenic_regulon_specificity",
        label="SCENIC regulon specificity",
        view=view,
    )


def scenic_binarization_plot_owned(result: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    return _regulatory_table_plot_owned(
        result,
        contract="scenic_binarization",
        node_id="OpenBioSingleCellSCENICBinarizationPlot",
        operation="scenic_binarization_plot",
        function_name="plot_scenic_binarization",
        label="SCENIC binarization",
        view=view,
    )


def scenic_membership_plot_owned(result: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    return _regulatory_table_plot_owned(
        result,
        contract="scenic_membership",
        node_id="OpenBioSingleCellSCENICRegulonMembershipPlot",
        operation="scenic_regulon_membership_plot",
        function_name="plot_scenic_regulon_membership",
        label="SCENIC regulon membership",
        view=view,
    )


@register_operation("openbio.node.tfactivityplot")
def tf_activity_plot_operation(context, inputs, parameters):
    require_input_names(inputs, {"adata", "activities"}, operation="TF Activity Plot")
    require_parameters(parameters, {"view"}, operation="TF Activity Plot")
    plotted, report, code = tf_activity_plot_owned(
        read_anndata_input(inputs), _tf_activity_input(inputs), **parameters
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.scenicactivityplot")
def scenic_activity_plot_operation(context, inputs, parameters):
    require_input_names(inputs, {"adata", "scenic_result"}, operation="SCENIC Activity Plot")
    require_parameters(parameters, {"view"}, operation="SCENIC Activity Plot")
    plotted, report, code = scenic_activity_plot_owned(
        read_anndata_input(inputs), _scenic_input(inputs), **parameters
    )
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


def _table_plot_operation(context, inputs, parameters, *, label: str, owned: Any):
    require_input_names(inputs, {"table"}, operation=label)
    require_parameters(parameters, {"view"}, operation=label)
    plotted, report, code = owned(_table_input(inputs, "table"), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.tfactivityrankingplot")
def tf_activity_ranking_plot_operation(context, inputs, parameters):
    return _table_plot_operation(
        context, inputs, parameters, label="TF Activity Ranking Plot", owned=tf_activity_ranking_plot_owned
    )


@register_operation("openbio.node.scenicregulonspecificityplot")
def scenic_rss_plot_operation(context, inputs, parameters):
    return _table_plot_operation(
        context, inputs, parameters, label="SCENIC Regulon Specificity Plot", owned=scenic_rss_plot_owned
    )


@register_operation("openbio.node.scenicbinarizationplot")
def scenic_binarization_plot_operation(context, inputs, parameters):
    return _table_plot_operation(
        context, inputs, parameters, label="SCENIC Binarization Plot", owned=scenic_binarization_plot_owned
    )


@register_operation("openbio.node.scenicregulonmembershipplot")
def scenic_membership_plot_operation(context, inputs, parameters):
    return _table_plot_operation(
        context, inputs, parameters, label="SCENIC Regulon Membership Plot", owned=scenic_membership_plot_owned
    )


__all__ = [
    "collectri_ulm",
    "collectri_ulm_owned",
    "import_pyscenic_owned",
    "import_pyscenic_results",
    "rank_tf_activities_owned",
    "rank_tf_activity",
    "scenic_activity_plot_operation",
    "scenic_activity_plot_owned",
    "scenic_activity_binarization",
    "scenic_binarization_plot_operation",
    "scenic_binarization_plot_owned",
    "scenic_binarization_owned",
    "scenic_membership_owned",
    "scenic_membership_plot_operation",
    "scenic_membership_plot_owned",
    "scenic_regulon_specificity",
    "scenic_rss_plot_operation",
    "scenic_rss_plot_owned",
    "scenic_rss_owned",
    "scenic_tf_modules",
    "tf_activity_plot_operation",
    "tf_activity_plot_owned",
    "tf_activity_ranking_plot_operation",
    "tf_activity_ranking_plot_owned",
]
