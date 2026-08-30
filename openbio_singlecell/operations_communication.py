from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from . import PLUGIN_VERSION
from .analysis_utils import make_plot_result, make_summary_result
from .expression_source import _LIANA_SPEC, DynamicExpressionSource
from .liana_artifact_codec import LIANA_CODEC, read_liana_result, write_liana_result
from .liana_communication import liana_communication_code, run_liana_communication
from .liana_plot import liana_dot_plot_code, render_liana_dot_plot
from .liana_result import validate_liana_result
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_artifact_input,
    require_file_input,
    require_input_names,
    require_parameters,
    write_plot_output,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

LIANA_KIND = "OPENBIO_LIANA_RESULT"
LIANA_EXPRESSION_SOURCE = _LIANA_SPEC


def _label_list(value: str, *, label: str) -> list[str]:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a comma-delimited string.")
    if not value.strip():
        return []
    labels = [part.strip() for part in value.split(",")]
    if any(not item for item in labels):
        raise ValueError(f"{label} contains an empty comma-delimited item.")
    if len(labels) != len(set(labels)):
        raise ValueError(f"{label} contains duplicate labels.")
    return labels


def liana_communication_owned(
    adata: Any,
    *,
    sample_key: str = "sample",
    condition_key: str = "condition",
    identity_key: str = "cell_type",
    annotation_status: str = "unknown",
    organism: str = "Homo sapiens",
    method: str = "rank_aggregate",
    resource_mode: str = "bundled_human",
    resource_name: str = "consensus",
    resource_path: str | Path | None = None,
    resource_metadata_json: str = "{}",
    source: DynamicExpressionSource | None = None,
    expression_proportion: float = 0.1,
    min_cells_per_identity_sample: int = 5,
    permutations: int = 1000,
    random_seed: int = 1337,
    jobs: int = 1,
    max_output_rows: int = 2_000_000,
    max_working_memory_gib: float = 4.0,
    liana_module: Any | None = None,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    expression = LIANA_EXPRESSION_SOURCE.resolve(adata, source)
    resolved_path = None if resource_path is None else str(Path(resource_path).resolve())
    result, summary = run_liana_communication(
        adata,
        copy_table=False,
        sample_key=sample_key,
        condition_key=condition_key,
        identity_key=identity_key,
        annotation_status=annotation_status,
        organism=organism,
        method=method,
        resource_mode=resource_mode,
        resource_name=resource_name,
        resource_path=resolved_path,
        resource_metadata_json=resource_metadata_json,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        expression_proportion=expression_proportion,
        min_cells_per_identity_sample=min_cells_per_identity_sample,
        permutations=permutations,
        random_seed=random_seed,
        jobs=jobs,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        openbio_version=PLUGIN_VERSION,
        liana_module=liana_module,
    )
    report = make_summary_result(
        summary=summary,
        title="Sample-resolved LIANA communication summary",
        operation="liana_communication",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    accounting = summary["key_results"]["resource"]["accounting"]
    code = liana_communication_code(
        parameters={
            "sample_key": sample_key,
            "condition_key": condition_key,
            "identity_key": identity_key,
            "annotation_status": annotation_status,
            "organism": organism,
            "method": method,
            "resource_mode": resource_mode,
            "resource_name": resource_name,
            "resource_path": resolved_path,
            "resource_metadata_json": resource_metadata_json,
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "expression_proportion": expression_proportion,
            "min_cells_per_identity_sample": min_cells_per_identity_sample,
            "permutations": permutations,
            "random_seed": random_seed,
            "jobs": jobs,
            "max_output_rows": max_output_rows,
            "max_working_memory_gib": max_working_memory_gib,
            "openbio_version": PLUGIN_VERSION,
            "expected_file_sha256": accounting["raw_file_sha256"],
            "expected_canonical_resource_sha256": accounting["canonical_resource_sha256"],
            "expected_effective_resource_sha256": accounting["effective_resource_sha256"],
        }
    )
    return result, report, code


def liana_dot_plot_owned(
    result: Any,
    *,
    source_labels: str = "",
    target_labels: str = "",
    selection_method: str = "rank_aggregate",
    selection_threshold: float = 0.05,
    top_n: int = 20,
    figure_width: float = 12.0,
    figure_height: float = 8.0,
    max_plot_rows: int = 100_000,
    max_image_pixels: int = 40_000_000,
    liana_module: Any | None = None,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    _table, provenance, _metadata = validate_liana_result(
        result,
        exact_type=False,
        copy_result=False,
    )
    if selection_method != provenance["method"]:
        raise ValueError(
            f"LIANA plot selection branch {selection_method!r} does not match the typed result method "
            f"{provenance['method']!r}."
        )
    source_label_list = _label_list(source_labels, label="LIANA source_labels")
    target_label_list = _label_list(target_labels, label="LIANA target_labels")
    png, summary = render_liana_dot_plot(
        result,
        exact_type=False,
        source_labels=source_label_list,
        target_labels=target_label_list,
        selection_threshold=selection_threshold,
        top_n=top_n,
        figure_width=figure_width,
        figure_height=figure_height,
        max_plot_rows=max_plot_rows,
        max_image_pixels=max_image_pixels,
        openbio_version=PLUGIN_VERSION,
        liana_module=liana_module,
    )
    cells = int(provenance["expression"]["cells"])
    genes = int(provenance["expression"]["genes"])
    common = {
        "parameters": summary["parameters"],
        "description": summary["results"],
        "warnings": summary["warnings"],
        "input_cells": cells,
        "input_genes": genes,
        "started_at": started_at,
    }
    plotted = make_plot_result(
        title=f"LIANA {provenance['method'].replace('_', ' ')} by-Sample dot plot",
        operation="liana_dot_plot",
        png=png,
        **common,
    )
    report = make_summary_result(
        summary=summary,
        title="LIANA dot-plot selection summary",
        operation="liana_dot_plot",
        **common,
    )
    code = liana_dot_plot_code(
        parameters={
            "source_labels": source_label_list,
            "target_labels": target_label_list,
            "selection_threshold": selection_threshold,
            "top_n": top_n,
            "figure_width": figure_width,
            "figure_height": figure_height,
            "max_plot_rows": max_plot_rows,
            "max_image_pixels": max_image_pixels,
            "openbio_version": PLUGIN_VERSION,
        }
    )
    return plotted, report, code


@register_operation("openbio.node.lianacommunication")
def liana_communication(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_parameters(
        parameters,
        {
            "sample_key",
            "condition_key",
            "identity_key",
            "annotation_status",
            "organism",
            "method",
            "resource_mode",
            "resource_name",
            "resource_metadata_json",
            "source",
            "expression_proportion",
            "min_cells_per_identity_sample",
            "permutations",
            "random_seed",
            "jobs",
            "max_output_rows",
            "max_working_memory_gib",
        },
        operation="LIANA Communication",
    )
    resource_mode = parameters["resource_mode"]
    expected_inputs = {"adata", "resource_csv"} if resource_mode == "local_resource" else {"adata"}
    require_input_names(inputs, expected_inputs, operation="LIANA Communication")
    resource_path = require_file_input(inputs, "resource_csv")[0] if "resource_csv" in inputs else None
    result, report, code = liana_communication_owned(
        read_anndata_input(inputs),
        resource_path=resource_path,
        **parameters,
    )
    root = context.create_output_directory("result")
    write_liana_result(root, result)
    return analysis_outputs(
        report,
        code,
        {
            "type": "artifact",
            "name": "result",
            "kind": LIANA_KIND,
            "codec": LIANA_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
    )


@register_operation("openbio.node.lianadotplot")
def liana_dot_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"result"}, operation="LIANA Dot Plot")
    require_parameters(
        parameters,
        {
            "source_labels",
            "target_labels",
            "selection_method",
            "selection_threshold",
            "top_n",
            "figure_width",
            "figure_height",
            "max_plot_rows",
            "max_image_pixels",
        },
        operation="LIANA Dot Plot",
    )
    result_root = require_artifact_input(
        inputs,
        "result",
        kind=LIANA_KIND,
        codec=LIANA_CODEC,
    )
    plotted, report, code = liana_dot_plot_owned(read_liana_result(result_root), **parameters)
    return analysis_outputs(
        report,
        code,
        write_plot_output(context, plotted, kind="OPENBIO_SINGLE_CELL_PLOT"),
    )


__all__ = [
    "liana_communication",
    "liana_communication_owned",
    "liana_dot_plot",
    "liana_dot_plot_owned",
]
