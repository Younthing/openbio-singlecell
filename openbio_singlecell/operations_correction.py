from __future__ import annotations

import hashlib
import time
from textwrap import dedent
from typing import Any, Literal

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, make_plot_result, validate_count_expression
from .correction_plotting import (
    mad_outlier_plot_code,
    run_mad_outlier_plot,
    run_scrublet_diagnostics_plot,
    scrublet_diagnostics_plot_code,
)
from .expression_source import _SCRUBLET_SPEC, ExpressionSource
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
)
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

CORRECTION_SOFTWARE_PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy")
SCRUBLET_SOURCE = _SCRUBLET_SPEC
SCANPY_REFERENCE = AnalysisReference(
    citation=(
        "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. "
        "Genome Biology. 2018;19:15."
    ),
    doi="10.1186/s13059-017-1382-0",
    url="https://doi.org/10.1186/s13059-017-1382-0",
    kind="software",
)
SCRUBLET_REFERENCE = AnalysisReference(
    citation=(
        "Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in "
        "Single-Cell Transcriptomic Data. Cell Systems. 2019;8(4):281-291.e9."
    ),
    doi="10.1016/j.cels.2018.11.005",
    url="https://doi.org/10.1016/j.cels.2018.11.005",
    kind="method",
)
LEYS_MAD_REFERENCE = AnalysisReference(
    citation=(
        "Leys C, Ley C, Klein O, Bernard P, Licata L. Detecting outliers: Do not use standard deviation "
        "around the mean, use absolute deviation around the median. Journal of Experimental Social Psychology. "
        "2013;49(4):764-766."
    ),
    doi="10.1016/j.jesp.2013.03.013",
    url="https://doi.org/10.1016/j.jesp.2013.03.013",
    kind="method",
)
SCIPY_MAD_REFERENCE = AnalysisReference(
    citation="SciPy documentation: scipy.stats.median_abs_deviation.",
    url="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html",
    kind="software_documentation",
)
ANNDATA_REFERENCE = AnalysisReference(
    citation="AnnData documentation: annotated data matrices and observation subsetting.",
    url="https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html",
    kind="software_documentation",
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)

def _input_anndata(
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
    *,
    operation: str,
    expected_parameters: set[str],
) -> Any:
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(parameters, expected_parameters, operation=operation)
    return read_anndata_input(inputs)


def _records(context: OperationContext, adata: Any, report: Any, code: str) -> list[JSONValue]:
    return analysis_outputs(report, code, write_anndata_output(context, adata))


def _string(value: JSONValue, *, name: str, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise ProtocolError(f"{name} must be a string.")
    result = value.strip()
    if not allow_empty and not result:
        raise ProtocolError(f"{name} cannot be empty.")
    return result


def _number(value: JSONValue, *, name: str, positive: bool = False) -> float:
    science = dependencies.require_scientific_dependencies()
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProtocolError(f"{name} must be a number.")
    number = float(value)
    if not bool(science.np.isfinite(number)) or (positive and number <= 0):
        suffix = " and greater than zero" if positive else ""
        raise ProtocolError(f"{name} must be finite{suffix}.")
    return number


def _integer(value: JSONValue, *, name: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ProtocolError(f"{name} must be an integer no smaller than {minimum}.")
    return value


def _boolean(value: JSONValue, *, name: str) -> bool:
    if not isinstance(value, bool):
        raise ProtocolError(f"{name} must be a boolean.")
    return value


def _expression_source_label(expression: ExpressionSource) -> str:
    if expression.kind == "layer":
        return f"AnnData layer {expression.layer_name!r}"
    if expression.kind == "raw":
        return "AnnData raw expression"
    return "AnnData.X"


def _expression_var(adata: Any, expression: ExpressionSource) -> Any:
    return adata.raw.var if expression.kind == "raw" else adata.var


def _group_positions(
    adata: Any,
    group_key: str,
    *,
    report_context: str,
) -> list[tuple[str, Any]]:
    science = dependencies.require_scientific_dependencies()
    if not group_key:
        return [("all_cells", science.np.arange(adata.n_obs, dtype=int))]
    _string_group_labels(adata.obs[group_key], group_key=group_key, report_context=report_context)
    grouped = adata.obs.groupby(group_key, observed=True, sort=False).indices
    return [(str(label), science.np.asarray(positions, dtype=int)) for label, positions in grouped.items()]


def _string_group_labels(values: Any, *, group_key: str, report_context: str) -> list[str]:
    science = dependencies.require_scientific_dependencies()
    rendered = [str(label) for label in science.pd.unique(values)]
    if len(set(rendered)) != len(rendered):
        raise ValueError(
            f"Sample/group labels in {group_key!r} collide when represented in {report_context}; "
            "use unambiguous labels with one consistent dtype."
        )
    return rendered


def _comma_separated_metrics(value: str) -> list[str]:
    metrics = list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))
    if not metrics:
        raise ValueError("At least one QC metric is required.")
    return metrics


def _observation_axis_fingerprint(names: Any) -> str:
    digest = hashlib.sha256()
    for value in names:
        encoded = str(value).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
    return f"sha256:{digest.hexdigest()}"


def _mad_outlier_statistics(
    values: Any,
    *,
    nmads: float,
    direction: Literal["both", "upper", "lower"],
    scale_mad: bool,
) -> tuple[Any, dict[str, int | float | None]]:
    science = dependencies.require_scientific_dependencies()
    from scipy.stats import median_abs_deviation

    array = science.np.asarray(values, dtype=float)
    finite = science.np.isfinite(array)
    finite_values = array[finite]
    median = float(science.np.median(finite_values))
    mad = float(median_abs_deviation(finite_values, scale="normal" if scale_mad else 1.0))
    lower = -science.np.inf if direction == "upper" else median - nmads * mad
    upper = science.np.inf if direction == "lower" else median + nmads * mad
    flags = finite & ((array < lower) | (array > upper))
    return flags, {
        "n": int(array.size),
        "finite": int(finite.sum()),
        "missing": int(array.size - finite.sum()),
        "median": median,
        "mad": mad,
        "lower_threshold": None if not science.np.isfinite(lower) else float(lower),
        "upper_threshold": None if not science.np.isfinite(upper) else float(upper),
        "flagged": int(flags.sum()),
    }


def _mad_outliers_code(parameters: dict[str, Any]) -> str:
    return dedent(
        f"""
        import numpy as np
        import pandas as pd
        from scipy.stats import median_abs_deviation


        def mark_mad_outliers(adata):
            metrics = {parameters["metrics"]!r}
            group_key = {parameters["batch_key"]!r}
            if group_key:
                if group_key not in adata.obs:
                    raise ValueError(f"QC Sample/group column not found in obs: {{group_key!r}}")
                if adata.obs[group_key].isna().any():
                    raise ValueError(f"QC Sample/group column contains missing labels: {{group_key!r}}")
                rendered = [str(label) for label in pd.unique(adata.obs[group_key])]
                if len(set(rendered)) != len(rendered):
                    raise ValueError(
                        f"Sample/group labels in {{group_key!r}} collide when represented in grouped MAD reports."
                    )
            groups = (
                adata.obs.groupby(group_key, observed=True, sort=False).indices.values()
                if group_key else [np.arange(adata.n_obs)]
            )
            flags = np.zeros(adata.n_obs, dtype=bool)
            for metric in metrics:
                values = np.asarray(adata.obs[metric], dtype=float)
                for positions in groups:
                    positions = np.asarray(positions, dtype=int)
                    selected = values[positions]
                    finite = np.isfinite(selected)
                    if int(finite.sum()) < {parameters["minimum_group_size"]!r}:
                        continue
                    median = float(np.median(selected[finite]))
                    mad = float(median_abs_deviation(
                        selected[finite], scale={"'normal'" if parameters["scale_mad"] else "1.0"}
                    ))
                    lower = -np.inf if {parameters["direction"]!r} == "upper" else median - {parameters["nmads"]!r} * mad
                    upper = np.inf if {parameters["direction"]!r} == "lower" else median + {parameters["nmads"]!r} * mad
                    flags[positions] |= finite & ((selected < lower) | (selected > upper))
            adata.obs[{parameters["output_column"]!r}] = flags
            return adata
        """
    ).strip()


def mark_mad_outliers_owned(
    adata: Any,
    *,
    metrics: str = "total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
    batch_key: str = "sample",
    nmads: float = 5.0,
    direction: str = "both",
    output_column: str = "outlier",
    scale_mad: bool = False,
    minimum_group_size: int = 3,
) -> tuple[Any, Any, str]:
    parameters: dict[str, JSONValue] = {
        "metrics": metrics,
        "batch_key": batch_key,
        "nmads": nmads,
        "direction": direction,
        "output_column": output_column,
        "scale_mad": scale_mad,
        "minimum_group_size": minimum_group_size,
    }
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    if cells == 0 or genes == 0:
        raise ValueError("MAD outlier marking requires at least one cell and one gene.")
    selected_metrics = _comma_separated_metrics(_string(parameters["metrics"], name="MAD metrics"))
    batch_key = _string(parameters["batch_key"], name="MAD batch_key", allow_empty=True)
    nmads = _number(parameters["nmads"], name="MAD multiplier")
    if nmads < 0:
        raise ProtocolError("MAD multiplier must be non-negative.")
    direction = _string(parameters["direction"], name="MAD direction")
    if direction not in {"both", "upper", "lower"}:
        raise ProtocolError(f"Unsupported MAD direction: {direction!r}")
    output_column = _string(parameters["output_column"], name="MAD output_column")
    scale_mad = _boolean(parameters["scale_mad"], name="MAD scale_mad")
    minimum_group_size = _integer(parameters["minimum_group_size"], name="MAD minimum_group_size", minimum=1)
    missing_metrics = [metric for metric in selected_metrics if metric not in adata.obs]
    if missing_metrics:
        raise ValueError(f"QC metrics not found in obs: {missing_metrics}")
    if batch_key and batch_key not in adata.obs:
        raise ValueError(f"QC Sample/group column not found in obs: {batch_key!r}")
    if batch_key and bool(adata.obs[batch_key].isna().any()):
        raise ValueError(f"QC Sample/group column contains missing labels: {batch_key!r}")
    science = dependencies.require_scientific_dependencies()
    for metric in selected_metrics:
        if not science.pd.api.types.is_numeric_dtype(adata.obs[metric].dtype):
            raise ValueError(f"QC metric must be numeric: {metric!r}")
    started_at = time.perf_counter()
    warnings: list[str] = []
    if output_column in adata.obs:
        warnings.append(f"Existing observation column {output_column!r} was overwritten.")
    if not batch_key:
        warnings.append("MAD thresholds were estimated globally because no Sample/group key was supplied.")
    groups = _group_positions(adata, batch_key, report_context="grouped MAD reports")
    flags = science.np.zeros(cells, dtype=bool)
    metric_statistics = []
    for metric in selected_metrics:
        values = science.np.asarray(adata.obs[metric], dtype=float)
        metric_flags = science.np.zeros(cells, dtype=bool)
        group_statistics = []
        for group_label, positions in groups:
            group_values = values[positions]
            finite_count = int(science.np.isfinite(group_values).sum())
            if finite_count < minimum_group_size:
                group_statistics.append(
                    {
                        "group": group_label,
                        "n": int(group_values.size),
                        "finite": finite_count,
                        "missing": int(group_values.size - finite_count),
                        "status": "skipped_insufficient_finite_values",
                        "median": None,
                        "mad": None,
                        "lower_threshold": None,
                        "upper_threshold": None,
                        "flagged": 0,
                    }
                )
                warnings.append(
                    f"Metric {metric!r}, group {group_label!r} had {finite_count} finite observations and was "
                    f"skipped because minimum_group_size={minimum_group_size}."
                )
                continue
            group_flags, statistics = _mad_outlier_statistics(
                group_values,
                nmads=nmads,
                direction=direction,  # type: ignore[arg-type]
                scale_mad=scale_mad,
            )
            metric_flags[positions] = group_flags
            status = "zero_mad" if statistics["mad"] == 0.0 else "evaluated"
            if status == "zero_mad":
                warnings.append(
                    f"Metric {metric!r}, group {group_label!r} had zero MAD; the effective threshold equals "
                    "the median and strict inequalities leave tied values unflagged."
                )
            group_statistics.append({"group": group_label, "status": status, **statistics})
        flags |= metric_flags
        metric_statistics.append(
            {
                "metric": metric,
                "flagged": int(metric_flags.sum()),
                "missing": int((~science.np.isfinite(values)).sum()),
                "groups": group_statistics,
            }
        )
    adata.obs[output_column] = flags
    report_parameters = {
        "metrics": selected_metrics,
        "batch_key": batch_key,
        "nmads": nmads,
        "direction": direction,
        "output_column": output_column,
        "scale_mad": scale_mad,
        "mad_scale": "normal" if scale_mad else "raw",
        "minimum_group_size": minimum_group_size,
    }
    finish_adata(adata, "mark_mad_outliers", report_parameters, cells, genes, started_at, warnings=warnings)
    group_union = [
        {
            "group": group_label,
            "cells": int(len(positions)),
            "flagged": int(flags[positions].sum()),
            "flagged_percent": float(flags[positions].mean() * 100.0),
        }
        for group_label, positions in groups
    ]
    threshold_rows = [
        {"metric": metric["metric"], **group}
        for metric in metric_statistics
        for group in metric["groups"]
    ]
    threshold_frame = science.pd.DataFrame.from_records(
        threshold_rows,
        columns=[
            "metric",
            "group",
            "status",
            "n",
            "finite",
            "missing",
            "median",
            "mad",
            "lower_threshold",
            "upper_threshold",
            "flagged",
        ],
    )
    for column in ("median", "mad", "lower_threshold", "upper_threshold"):
        threshold_frame[column] = science.pd.to_numeric(threshold_frame[column], errors="coerce")
    adata.uns["openbio_mad_outliers"] = {
        "schema_version": 1,
        "producer_operation": "mark_mad_outliers",
        "observation_axis_fingerprint_sha256": _observation_axis_fingerprint(adata.obs_names),
        "metrics": selected_metrics,
        "batch_key": batch_key,
        "output_column": output_column,
        "nmads": nmads,
        "direction": direction,
        "scale_mad": scale_mad,
        "minimum_group_size": minimum_group_size,
        "thresholds": threshold_frame,
        "group_rates": science.pd.DataFrame.from_records(
            (
                {
                    "group": group["group"],
                    "cells": group["cells"],
                    "marked": group["flagged"],
                    "marked_fraction": group["flagged"] / group["cells"],
                }
                for group in group_union
            ),
            columns=["group", "cells", "marked", "marked_fraction"],
        ),
    }
    flagged_cells = int(flags.sum())
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMarkMADOutliers",
        title="MAD outlier marking summary",
        operation="mark_mad_outliers_report",
        methods=(
            f"MAD outliers were marked for {', '.join(selected_metrics)} using a {direction}-direction rule at "
            f"{nmads:g} MADs; metric flags were combined by union."
        ),
        results=(
            f"{flagged_cells:,} of {cells:,} cells ({flagged_cells * 100.0 / cells:.1f}%) were marked in "
            f"obs[{output_column!r}]. No cells were removed."
        ),
        key_results={
            "cells": cells,
            "genes": genes,
            "flagged_cells": flagged_cells,
            "unflagged_cells": cells - flagged_cells,
            "flagged_percent": float(flagged_cells * 100.0 / cells),
            "combination_rule": "union",
            "output_column": output_column,
            "metric_statistics": metric_statistics,
            "group_union": group_union,
        },
        parameters=report_parameters,
        references=(SCIPY_MAD_REFERENCE, LEYS_MAD_REFERENCE),
        software_packages=("scipy", "numpy", "pandas", "anndata"),
        warnings=warnings,
        limitations=(
            "MAD cutoffs are dataset- and Sample-dependent QC flags, not universal low-quality-cell definitions.",
            "Metrics with different biological directions should be evaluated in separate node invocations.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_mad_outliers_code(report_parameters),
    )
    return adata, report, code


def mad_outlier_plot_owned(adata: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_mad_outlier_plot(adata, view=view)
    parameters = {"view": details["view_parameters"]}
    code = mad_outlier_plot_code(view=details["view_parameters"])
    warnings = []
    if details["missing_metric_values"]:
        warnings.append(
            f"Excluded {details['missing_metric_values']:,} non-finite metric value(s) from histogram rendering."
        )
    description = (
        f"Rendered {details['view'].replace('_', ' ')} for {details['plotted_observations']:,} cells, "
        f"{len(details['metrics']):,} metric(s), and {len(details['groups']):,} Sample/group(s)."
    )
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "started_at": started_at,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="mad_outlier_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMADOutlierPlot",
        title="MAD outlier plot summary",
        operation="mad_outlier_plot",
        methods=(
            "Validated the stored Mark MAD Outliers producer provenance, observation axis, threshold grid, "
            "metric values, and union flag before read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=(LEYS_MAD_REFERENCE, SCIPY_MAD_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("numpy", "pandas", "scipy", "matplotlib", "anndata"),
        warnings=warnings,
        limitations=(
            "MAD flags are dataset- and Sample-dependent QC evidence, not universal low-quality-cell definitions.",
            "The plot visualizes stored thresholds and descriptive marked fractions; it does not choose thresholds or remove cells.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


@register_operation("openbio.node.markmadoutliers")
def mark_mad_outliers(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    expected = {"metrics", "batch_key", "nmads", "direction", "output_column", "scale_mad", "minimum_group_size"}
    adata = _input_anndata(
        inputs,
        parameters,
        operation="MAD outlier marking",
        expected_parameters=expected,
    )
    return _records(context, *mark_mad_outliers_owned(adata, **parameters))


@register_operation("openbio.node.madoutlierplot")
def mad_outlier_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="MAD Outlier Plot")
    require_parameters(parameters, {"view"}, operation="MAD Outlier Plot")
    plotted, report, code = mad_outlier_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind="OPENBIO_SINGLE_CELL_PLOT"))


def _scrublet_group_payloads(scrublet_uns: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(scrublet_uns, dict):
        raise RuntimeError("Scanpy Scrublet did not produce an uns['scrublet'] result mapping.")
    batches = scrublet_uns.get("batches")
    if batches is None:
        return {"all_cells": scrublet_uns}
    if not isinstance(batches, dict):
        raise RuntimeError("Scanpy Scrublet produced an invalid batched result mapping.")
    return {str(label): payload for label, payload in batches.items() if isinstance(payload, dict)}


def _scrublet_thresholds(scrublet_uns: Any, expected_groups: list[str]) -> dict[str, float]:
    science = dependencies.require_scientific_dependencies()
    payloads = _scrublet_group_payloads(scrublet_uns)
    if set(payloads) != set(expected_groups):
        missing = sorted(set(expected_groups) - set(payloads))
        raise RuntimeError(f"Scanpy Scrublet did not return results for every Sample/group: {missing}")
    thresholds: dict[str, float] = {}
    missing_thresholds = []
    for label, payload in payloads.items():
        try:
            numeric = float(payload.get("threshold"))
        except (TypeError, ValueError):
            missing_thresholds.append(label)
            continue
        if not bool(science.np.isfinite(numeric)):
            missing_thresholds.append(label)
            continue
        thresholds[label] = numeric
    if missing_thresholds:
        raise RuntimeError(
            "Scrublet did not determine a usable threshold for Sample/group(s) "
            f"{missing_thresholds}. Inspect the simulated score distribution and rerun with a manual threshold."
        )
    return thresholds


def _scrublet_simulated_scores(scrublet_uns: Any, expected_groups: list[str]) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    payloads = _scrublet_group_payloads(scrublet_uns)
    if set(payloads) != set(expected_groups):
        missing = sorted(set(expected_groups) - set(payloads))
        raise RuntimeError(f"Scanpy Scrublet did not return results for every Sample/group: {missing}")
    simulated = {}
    for label in expected_groups:
        values = science.np.asarray(payloads[label].get("doublet_scores_sim"))
        if values.ndim != 1 or values.size < 1:
            raise RuntimeError(
                f"Scanpy Scrublet did not retain a one-dimensional simulated score sequence for Sample/group {label!r}."
            )
        if (
            not bool(science.np.issubdtype(values.dtype, science.np.number))
            or bool(science.np.iscomplexobj(values))
            or not bool(science.np.isfinite(values).all())
        ):
            raise RuntimeError(
                f"Scanpy Scrublet produced invalid simulated doublet scores for Sample/group {label!r}."
            )
        simulated[label] = values.astype(float, copy=False)
    return simulated


def _scrublet_score_evidence_fingerprint(
    scores: Any,
    predictions: Any,
    thresholds: dict[str, float],
    simulated_scores: dict[str, Any],
    groups: list[str],
) -> str:
    science = dependencies.require_scientific_dependencies()
    digest = hashlib.sha256(b"scrublet-score-evidence-v1")
    for values, dtype in ((scores, science.np.float64), (predictions, science.np.uint8)):
        array = science.np.asarray(values, dtype=dtype)
        digest.update(str(tuple(array.shape)).encode("ascii"))
        digest.update(science.np.ascontiguousarray(array).tobytes())
    for group in groups:
        encoded = group.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
        digest.update(science.np.asarray([thresholds[group]], dtype=science.np.float64).tobytes())
        simulated = science.np.asarray(simulated_scores[group], dtype=science.np.float64)
        digest.update(str(tuple(simulated.shape)).encode("ascii"))
        digest.update(science.np.ascontiguousarray(simulated).tobytes())
    return f"sha256:{digest.hexdigest()}"


def _scrublet_code(expression: ExpressionSource, parameters: dict[str, Any]) -> str:
    if expression.kind == "layer":
        matrix = f"adata.layers[{expression.layer_name!r}]"
    elif expression.kind == "raw":
        matrix = "adata.raw.X"
    else:
        matrix = "adata.X"
    var = "adata.raw.var.copy()" if expression.kind == "raw" else "adata.var.copy()"
    layer = expression.layer_name if expression.kind == "layer" else None
    use_raw = expression.kind == "raw"
    return dedent(
        f"""
        import warnings

        import numpy as np
        import pandas as pd
        import scanpy as sc
        from anndata import AnnData
        from scipy import sparse


        def _validated_scrublet_thresholds(scrublet_uns, expected_groups):
            if not isinstance(scrublet_uns, dict):
                raise RuntimeError("Scanpy Scrublet did not produce an uns['scrublet'] result mapping.")
            batches = scrublet_uns.get("batches")
            if batches is None:
                payloads = {{"all_cells": scrublet_uns}}
            elif isinstance(batches, dict):
                payloads = {{str(label): value for label, value in batches.items() if isinstance(value, dict)}}
            else:
                raise RuntimeError("Scanpy Scrublet produced an invalid batched result mapping.")
            if set(payloads) != set(expected_groups):
                missing = sorted(set(expected_groups) - set(payloads))
                raise RuntimeError(f"Scanpy Scrublet did not return results for every Sample/group: {{missing}}")
            missing_thresholds = []
            for label, payload in payloads.items():
                try:
                    threshold = float(payload.get("threshold"))
                except (TypeError, ValueError):
                    missing_thresholds.append(label)
                    continue
                if not np.isfinite(threshold):
                    missing_thresholds.append(label)
            if missing_thresholds:
                raise RuntimeError(
                    "Scrublet did not determine a usable threshold for Sample/group(s) "
                    f"{{missing_thresholds}}. Inspect the simulated score distribution and rerun with a manual threshold."
                )


        def run_scrublet(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("Scrublet requires at least one cell and one gene.")
            layer = {layer!r}
            if layer is not None and layer not in adata.layers:
                raise ValueError(f"Expression layer not found: {{layer!r}}")
            if {use_raw!r} and adata.raw is None:
                raise ValueError("Raw expression is unavailable.")
            group_key = {parameters["batch_key"]!r}
            if group_key:
                if group_key not in adata.obs:
                    raise ValueError(f"Scrublet Sample/group column not found in obs: {{group_key!r}}")
                if adata.obs[group_key].isna().any():
                    raise ValueError(f"Scrublet Sample/group column contains missing labels: {{group_key!r}}")
                rendered_labels = [str(label) for label in pd.unique(adata.obs[group_key])]
                if len(set(rendered_labels)) != len(rendered_labels):
                    raise ValueError(
                        f"Sample/group labels in {{group_key!r}} collide when represented in Scrublet provenance; "
                        "use unambiguous labels with one consistent dtype."
                    )
                group_sizes = adata.obs.groupby(group_key, observed=True, sort=False).size()
                minimum_cells = {parameters["n_prin_comps"]!r} + 1
                too_small = {{
                    str(label): int(size) for label, size in group_sizes.items() if int(size) < minimum_cells
                }}
                if too_small:
                    raise ValueError(
                        f"Scrublet n_prin_comps={parameters['n_prin_comps']!r} requires at least "
                        f"{{minimum_cells}} cells in every Sample/group: {{too_small}}"
                    )
                expected_groups = rendered_labels
            else:
                minimum_cells = {parameters["n_prin_comps"]!r} + 1
                if adata.n_obs < minimum_cells:
                    raise ValueError(
                        f"Scrublet n_prin_comps={parameters['n_prin_comps']!r} requires at least "
                        f"{{minimum_cells}} cells."
                    )
                expected_groups = ["all_cells"]
            matrix = {matrix}
            if matrix.shape[1] <= {parameters["n_prin_comps"]!r}:
                raise ValueError(
                    f"Scrublet n_prin_comps={parameters['n_prin_comps']!r} must be smaller than the selected "
                    f"source feature count ({{matrix.shape[1]}})."
                )
            values = np.asarray(matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel())
            if not np.isfinite(values).all():
                raise ValueError("Selected Scrublet source contains non-finite expression values.")
            if values.size == 0 or not (values > 0).any():
                raise ValueError("Selected Scrublet source contains no positive expression values.")
            if (values < 0).any():
                warnings.warn(
                    "Selected Scrublet source contains negative expression values; the selected expression was retained, "
                    "but it is not a conventional unnormalized count matrix.",
                    UserWarning,
                    stacklevel=2,
                )
            work_obs = adata.obs.copy()
            if not adata.obs_names.is_unique:
                work_obs.index = [f"__openbio_scrublet_row_{{index}}" for index in range(adata.n_obs)]
            work = AnnData(X={matrix}.copy(), obs=work_obs, var={var})
            sc.pp.scrublet(
                work,
                batch_key=group_key or None,
                sim_doublet_ratio={parameters["sim_doublet_ratio"]!r},
                expected_doublet_rate={parameters["expected_doublet_rate"]!r},
                stdev_doublet_rate=0.02,
                synthetic_doublet_umi_subsampling=1.0,
                knn_dist_metric="euclidean",
                normalize_variance=True,
                log_transform=False,
                mean_center=True,
                n_prin_comps={parameters["n_prin_comps"]!r},
                use_approx_neighbors=None,
                get_doublet_neighbor_parents=False,
                n_neighbors={parameters["n_neighbors"]!r},
                threshold={parameters["threshold"]!r},
                verbose=False,
                copy=False,
                random_state={parameters["random_seed"]!r},
            )
            _validated_scrublet_thresholds(work.uns.get("scrublet"), expected_groups)
            if "doublet_score" not in work.obs or "predicted_doublet" not in work.obs:
                raise RuntimeError("Scanpy Scrublet did not produce its documented observation columns.")
            scores = np.asarray(work.obs["doublet_score"], dtype=float)
            predictions = work.obs["predicted_doublet"]
            if not np.isfinite(scores).all():
                raise RuntimeError("Scanpy Scrublet produced non-finite observed doublet scores.")
            if not pd.api.types.is_bool_dtype(predictions.dtype) or predictions.isna().any():
                raise RuntimeError("Scanpy Scrublet predicted_doublet output is not a complete boolean column.")
            adata.obs["doublet_score"] = scores
            adata.obs["predicted_doublet"] = predictions.to_numpy(dtype=bool)
            adata.uns["scrublet"] = work.uns["scrublet"]
            return adata
        """
    ).strip()


def scrublet_owned(
    adata: Any,
    *,
    batch_key: str = "sample",
    random_seed: int = 123,
    expected_doublet_rate: float = 0.05,
    threshold_mode: str = "automatic",
    threshold: float = 0.25,
    sim_doublet_ratio: float = 2.0,
    n_prin_comps: int = 30,
    n_neighbors: int = 0,
    source: dict[str, JSONValue] | None = None,
) -> tuple[Any, Any, str]:
    parameters: dict[str, JSONValue] = {
        "batch_key": batch_key,
        "random_seed": random_seed,
        "expected_doublet_rate": expected_doublet_rate,
        "threshold_mode": threshold_mode,
        "threshold": threshold,
        "sim_doublet_ratio": sim_doublet_ratio,
        "n_prin_comps": n_prin_comps,
        "n_neighbors": n_neighbors,
        "source": {"source": "X"} if source is None else source,
    }
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    if cells == 0 or genes == 0:
        raise ValueError("Scrublet requires at least one cell and one gene.")
    batch_key = _string(parameters["batch_key"], name="Scrublet batch_key", allow_empty=True)
    random_seed = _integer(parameters["random_seed"], name="Scrublet random_seed", minimum=0)
    expected_doublet_rate = _number(parameters["expected_doublet_rate"], name="Scrublet expected_doublet_rate")
    if not 0.0 < expected_doublet_rate < 1.0:
        raise ValueError("Scrublet expected doublet rate must be between zero and one.")
    threshold_mode = _string(parameters["threshold_mode"], name="Scrublet threshold_mode")
    if threshold_mode not in {"automatic", "manual"}:
        raise ProtocolError(f"Unsupported Scrublet threshold mode: {threshold_mode!r}")
    manual_threshold = (
        _number(parameters["threshold"], name="Scrublet manual threshold")
        if threshold_mode == "manual"
        else None
    )
    sim_doublet_ratio = _number(
        parameters["sim_doublet_ratio"],
        name="Scrublet simulated-doublet ratio",
        positive=True,
    )
    n_prin_comps = _integer(parameters["n_prin_comps"], name="Scrublet n_prin_comps", minimum=1)
    n_neighbors = _integer(parameters["n_neighbors"], name="Scrublet n_neighbors", minimum=0)
    source_value = parameters["source"]
    if not isinstance(source_value, dict):
        raise ProtocolError("Scrublet source must be a JSON object.")
    source_kind = source_value.get("source")
    expected_source_fields = (
        {"source", SCRUBLET_SOURCE.layer_input_id} if source_kind == "layer" else {"source"}
    )
    if set(source_value) != expected_source_fields:
        raise ProtocolError(
            "Scrublet source fields must be exactly "
            f"{sorted(expected_source_fields)!r}; got {sorted(source_value)!r}."
        )
    try:
        expression = SCRUBLET_SOURCE.resolve(adata, source_value)
    except (TypeError, ValueError) as error:
        raise ProtocolError(str(error)) from error
    warnings: list[str] = []
    if manual_threshold is not None and not 0.0 <= manual_threshold <= 1.0:
        warnings.append(
            f"Manual Scrublet threshold {manual_threshold:g} lies outside the conventional [0, 1] score "
            "interval; the explicit extreme threshold was retained."
        )
    duplicate_cell_ids = not adata.obs_names.is_unique
    if duplicate_cell_ids:
        warnings.append(
            "Observation identifiers are not unique. Scrublet used temporary positional identifiers; the "
            "original row order and duplicate Cell IDs were preserved in the output."
        )
    overwritten_fields = [
        field
        for field, present in (
            ("obs['doublet_score']", "doublet_score" in adata.obs),
            ("obs['predicted_doublet']", "predicted_doublet" in adata.obs),
            ("uns['scrublet']", "scrublet" in adata.uns),
            ("uns['openbio_scrublet_diagnostics']", "openbio_scrublet_diagnostics" in adata.uns),
        )
        if present
    ]
    if overwritten_fields:
        warnings.append("Existing Scrublet result fields were overwritten: " + ", ".join(overwritten_fields) + ".")
    if batch_key:
        if batch_key not in adata.obs:
            raise ValueError(f"Scrublet Sample/group column not found in obs: {batch_key!r}")
        if bool(adata.obs[batch_key].isna().any()):
            raise ValueError(f"Scrublet Sample/group column contains missing labels: {batch_key!r}")
        expected_groups = _string_group_labels(
            adata.obs[batch_key],
            group_key=batch_key,
            report_context="Scrublet provenance",
        )
        group_sizes = adata.obs.groupby(batch_key, observed=True, sort=False).size()
        too_small = {
            str(label): int(size) for label, size in group_sizes.items() if int(size) < n_prin_comps + 1
        }
        if too_small:
            raise ValueError(
                f"Scrublet n_prin_comps={n_prin_comps} requires at least {n_prin_comps + 1} cells in every "
                f"Sample/group: {too_small}"
            )
    else:
        if cells < n_prin_comps + 1:
            raise ValueError(f"Scrublet n_prin_comps={n_prin_comps} requires at least {n_prin_comps + 1} cells.")
        expected_groups = ["all_cells"]
        warnings.append(
            "Scrublet was run globally because no Sample/group key was supplied; merged captures can distort "
            "simulated doublet composition."
        )
    source_label = _expression_source_label(expression)
    matrix = expression.matrix(adata)
    if matrix.shape[1] <= n_prin_comps:
        raise ValueError(
            f"Scrublet n_prin_comps={n_prin_comps} must be smaller than the selected source feature count "
            f"({matrix.shape[1]})."
        )
    warnings.extend(
        validate_count_expression(matrix, source_label=source_label, require_positive=True, require_nonnegative=False)
    )
    work_obs = adata.obs.copy()
    if duplicate_cell_ids:
        work_obs.index = [f"__openbio_scrublet_row_{index}" for index in range(cells)]
    # Scrublet mutates and filters internally; this dedicated matrix/axis workspace is algorithmically required.
    work = science.ad.AnnData(X=matrix.copy(), obs=work_obs, var=_expression_var(adata, expression).copy())
    try:
        science.sc.pp.scrublet(
            work,
            batch_key=batch_key or None,
            sim_doublet_ratio=sim_doublet_ratio,
            expected_doublet_rate=expected_doublet_rate,
            stdev_doublet_rate=0.02,
            synthetic_doublet_umi_subsampling=1.0,
            knn_dist_metric="euclidean",
            normalize_variance=True,
            log_transform=False,
            mean_center=True,
            n_prin_comps=n_prin_comps,
            use_approx_neighbors=None,
            get_doublet_neighbor_parents=False,
            n_neighbors=n_neighbors or None,
            threshold=manual_threshold,
            verbose=False,
            copy=False,
            random_state=random_seed,
        )
    except ValueError as error:
        if "n_components" in str(error) or "n_prin_comps" in str(error):
            raise ValueError(
                "Scrublet principal-component dimension is too large after internal filtering; reduce "
                "n_prin_comps or provide a larger Sample/capture."
            ) from error
        raise
    scrublet_uns = work.uns.get("scrublet")
    thresholds = _scrublet_thresholds(scrublet_uns, expected_groups)
    simulated_scores = _scrublet_simulated_scores(scrublet_uns, expected_groups)
    if "doublet_score" not in work.obs or "predicted_doublet" not in work.obs:
        raise RuntimeError("Scanpy Scrublet did not produce its documented observation columns.")
    scores = science.np.asarray(work.obs["doublet_score"], dtype=float)
    predictions = work.obs["predicted_doublet"]
    if not bool(science.np.isfinite(scores).all()):
        raise RuntimeError("Scanpy Scrublet produced non-finite observed doublet scores.")
    if not science.pd.api.types.is_bool_dtype(predictions.dtype) or bool(predictions.isna().any()):
        raise RuntimeError("Scanpy Scrublet predicted_doublet output is not a complete boolean column.")
    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = predictions.to_numpy(dtype=bool)
    adata.uns["scrublet"] = work.uns["scrublet"]
    report_parameters = {
        **expression.parameters(),
        "batch_key": batch_key,
        "random_seed": random_seed,
        "expected_doublet_rate": expected_doublet_rate,
        "threshold_mode": threshold_mode,
        "threshold": manual_threshold,
        "sim_doublet_ratio": sim_doublet_ratio,
        "n_prin_comps": n_prin_comps,
        "n_neighbors": n_neighbors or None,
        "stdev_doublet_rate": 0.02,
        "synthetic_doublet_umi_subsampling": 1.0,
        "knn_dist_metric": "euclidean",
        "normalize_variance": True,
        "log_transform": False,
        "mean_center": True,
        "use_approx_neighbors": None,
        "get_doublet_neighbor_parents": False,
        "verbose": False,
        "copy": False,
        "random_api": "random_state",
    }
    predicted_count = int(predictions.sum())
    if predicted_count == 0:
        warnings.append("Scrublet predicted no doublets; verify the score distribution before filtering.")
    finish_adata(
        adata,
        "scrublet",
        report_parameters,
        cells,
        genes,
        started_at,
        random_seed=random_seed,
        warnings=warnings,
    )
    groups = _group_positions(adata, batch_key, report_context="Scrublet reports")
    group_payloads = _scrublet_group_payloads(scrublet_uns)
    per_group = []
    for group_label, positions in groups:
        group_predictions = science.np.asarray(predictions.iloc[positions], dtype=bool)
        group_scores = scores[positions]
        per_group.append(
            {
                "group": group_label,
                "cells": int(len(positions)),
                "predicted_doublets": int(group_predictions.sum()),
                "predicted_percent": float(group_predictions.mean() * 100.0),
                "doublet_score": summarize_numeric(group_scores),
                "threshold": thresholds[group_label],
                "effective_parameters": group_payloads[group_label].get("parameters", {}),
            }
        )
    adata.uns["openbio_scrublet_diagnostics"] = {
        "schema_version": 1,
        "producer_operation": "scrublet",
        "observation_axis_fingerprint_sha256": _observation_axis_fingerprint(adata.obs_names),
        "batch_key": batch_key,
        "score_column": "doublet_score",
        "prediction_column": "predicted_doublet",
        "score_evidence_fingerprint_sha256": _scrublet_score_evidence_fingerprint(
            scores,
            predictions.to_numpy(dtype=bool),
            thresholds,
            simulated_scores,
            [group["group"] for group in per_group],
        ),
        "group_statistics": science.pd.DataFrame.from_records(
            (
                {
                    "group": group["group"],
                    "cells": group["cells"],
                    "predicted_doublets": group["predicted_doublets"],
                    "predicted_fraction": group["predicted_doublets"] / group["cells"],
                    "threshold": group["threshold"],
                    "simulated_count": int(simulated_scores[group["group"]].size),
                }
                for group in per_group
            ),
            columns=[
                "group",
                "cells",
                "predicted_doublets",
                "predicted_fraction",
                "threshold",
                "simulated_count",
            ],
        ),
    }
    predicted_percent = float(predicted_count * 100.0 / cells)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellScrublet",
        title="Scrublet doublet prediction summary",
        operation="scrublet_report",
        methods=(
            f"Scanpy Scrublet scored {cells:,} transcriptomes from {source_label} and simulated doublets "
            f"{'within ' + batch_key if batch_key else 'globally'}."
        ),
        results=(
            f"Scrublet predicted {predicted_count:,} of {cells:,} cells ({predicted_percent:.1f}%) as doublets. "
            "Cells were annotated but not removed."
        ),
        key_results={
            "input_cells": cells,
            "selected_source_genes": int(work.n_vars),
            "output_genes": genes,
            "duplicate_observation_identifiers": duplicate_cell_ids,
            "source": source_label,
            "groups": len(groups),
            "predicted_doublets": predicted_count,
            "predicted_percent": predicted_percent,
            "doublet_score": summarize_numeric(scores),
            "effective_thresholds": thresholds,
            "per_group": per_group,
        },
        parameters=report_parameters,
        references=(SCRUBLET_REFERENCE, SCANPY_REFERENCE),
        software_packages=(*CORRECTION_SOFTWARE_PACKAGES, "scikit-image"),
        warnings=warnings,
        limitations=(
            "Scrublet predictions are model-based QC evidence, not ground-truth doublet labels.",
            "Homotypic doublets can be difficult to detect.",
        ),
        input_cells=cells,
        input_genes=int(work.n_vars),
        started_at=started_at,
        code=_scrublet_code(expression, report_parameters),
    )
    return adata, report, code


def scrublet_diagnostics_plot_owned(adata: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_scrublet_diagnostics_plot(adata, view=view)
    parameters = {"view": details["view_parameters"]}
    code = scrublet_diagnostics_plot_code(view=details["view_parameters"])
    description = (
        f"Rendered {details['view'].replace('_', ' ')} for {details['plotted_observed_scores']:,} observed "
        f"scores, {details['plotted_simulated_scores']:,} simulated scores, and "
        f"{len(details['groups']):,} Sample/group(s)."
    )
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": [],
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "started_at": started_at,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="scrublet_diagnostics_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellScrubletDiagnosticsPlot",
        title="Scrublet diagnostics plot summary",
        operation="scrublet_diagnostics_plot",
        methods=(
            "Validated the Run Scrublet producer provenance, observation axis, observed scores, complete "
            "simulated-score sequences, effective thresholds, predictions, and Sample/group statistics before "
            "read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=(SCRUBLET_REFERENCE, SCANPY_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("scanpy", "numpy", "pandas", "matplotlib", "anndata"),
        warnings=[],
        limitations=(
            "Scrublet predictions are model-based QC evidence, not ground-truth doublet labels.",
            "Observed-versus-simulated separation and predicted fractions do not establish a universal filtering threshold.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


@register_operation("openbio.node.scrublet")
def scrublet(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    expected = {
        "batch_key",
        "random_seed",
        "expected_doublet_rate",
        "threshold_mode",
        "threshold",
        "sim_doublet_ratio",
        "n_prin_comps",
        "n_neighbors",
        "source",
    }
    adata = _input_anndata(inputs, parameters, operation="Scrublet", expected_parameters=expected)
    return _records(context, *scrublet_owned(adata, **parameters))


@register_operation("openbio.node.scrubletdiagnosticsplot")
def scrublet_diagnostics_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Scrublet Diagnostics Plot")
    require_parameters(parameters, {"view"}, operation="Scrublet Diagnostics Plot")
    plotted, report, code = scrublet_diagnostics_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind="OPENBIO_SINGLE_CELL_PLOT"))


def _filter_doublets_code(prediction_column: str) -> str:
    return dedent(
        f"""
        import pandas as pd


        def filter_predicted_doublets(adata):
            column = {prediction_column!r}
            predictions = adata.obs[column]
            if not pd.api.types.is_bool_dtype(predictions.dtype) or predictions.isna().any():
                raise ValueError(f"Doublet prediction column must be complete boolean data: {{column!r}}")
            return adata[~predictions.to_numpy(dtype=bool)].copy()
        """
    ).strip()


def filter_doublets_owned(
    adata: Any,
    *,
    prediction_column: str = "predicted_doublet",
) -> tuple[Any, Any, str]:
    parameters: dict[str, JSONValue] = {"prediction_column": prediction_column}
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    if cells == 0 or genes == 0:
        raise ValueError("Doublet filtering requires at least one cell and one gene.")
    prediction_column = _string(parameters["prediction_column"], name="Doublet prediction_column")
    if prediction_column not in adata.obs:
        raise ValueError(f"Doublet prediction column not found in obs: {prediction_column!r}")
    science = dependencies.require_scientific_dependencies()
    predictions = adata.obs[prediction_column]
    if not science.pd.api.types.is_bool_dtype(predictions.dtype):
        raise ValueError(f"Doublet prediction column must have boolean dtype: {prediction_column!r}")
    if bool(predictions.isna().any()):
        raise ValueError(f"Doublet prediction column contains missing values: {prediction_column!r}")
    predicted_values = predictions.to_numpy(dtype=bool)
    keep = ~predicted_values
    retained_cells = int(keep.sum())
    # AnnData slicing returns a view; filtering must materialize the independent result artifact.
    adata = adata[keep].copy()
    report_parameters = {"prediction_column": prediction_column}
    predicted_count = int(predicted_values.sum())
    warnings: list[str] = []
    if retained_cells == 0:
        warnings.append(
            "Every input cell was marked for removal. The empty result is preserved as the explicit filtering "
            "outcome and will not be usable by most downstream analyses."
        )
    elif predicted_count == 0:
        warnings.append("No cells were marked as doublets; the output retains every input cell.")
    finish_adata(adata, "filter_doublets", report_parameters, cells, genes, started_at, warnings=warnings)
    removed_percent = float(predicted_count * 100.0 / cells)
    retained_percent = float(retained_cells * 100.0 / cells)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellFilterDoublets",
        title="Predicted-doublet filtering summary",
        operation="filter_doublets_report",
        methods=(
            f"Cells with boolean True in obs[{prediction_column!r}] were removed by inverse-mask subsetting. "
            "The prediction column was consumed without coercion."
        ),
        results=(
            f"Removed {predicted_count:,} of {cells:,} cells ({removed_percent:.1f}%) and retained "
            f"{retained_cells:,} cells ({retained_percent:.1f}%); all {genes:,} genes were preserved."
        ),
        key_results={
            "input_cells": cells,
            "predicted_doublets_removed": predicted_count,
            "removed_percent": removed_percent,
            "retained_singlet_candidates": retained_cells,
            "retained_percent": retained_percent,
            "genes": genes,
            "prediction_column": prediction_column,
            "prediction_dtype": str(predictions.dtype),
        },
        parameters=report_parameters,
        references=(ANNDATA_REFERENCE,),
        software_packages=("anndata", "numpy", "pandas"),
        warnings=warnings,
        limitations=(
            "Filtering inherits the assumptions, threshold choice, and errors of the upstream doublet caller.",
            "Predicted singlets can still include homotypic or otherwise undetected doublets.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_filter_doublets_code(prediction_column),
    )
    return adata, report, code


@register_operation("openbio.node.filterdoublets")
def filter_doublets(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    adata = _input_anndata(
        inputs,
        parameters,
        operation="Doublet filtering",
        expected_parameters={"prediction_column"},
    )
    return _records(context, *filter_doublets_owned(adata, **parameters))


__all__ = [
    "filter_doublets",
    "filter_doublets_owned",
    "mad_outlier_plot",
    "mad_outlier_plot_owned",
    "mark_mad_outliers",
    "mark_mad_outliers_owned",
    "scrublet",
    "scrublet_diagnostics_plot",
    "scrublet_diagnostics_plot_owned",
    "scrublet_owned",
]
