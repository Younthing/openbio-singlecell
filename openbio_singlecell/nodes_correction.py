from __future__ import annotations

import copy
import time
import warnings as python_warnings
from textwrap import dedent
from typing import TYPE_CHECKING, Literal

from comfy_api.latest import io

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, validate_count_expression
from .expression_source import DynamicExpressionSource, ExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, SummaryResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/correction"
MAX_RANDOM_SEED = 2**31 - 1
CORRECTION_SOFTWARE_PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy")
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


def _expression_source_label(expression: ExpressionSource) -> str:
    if expression.kind == "layer":
        return f"AnnData layer {expression.layer_name!r}"
    if expression.kind == "raw":
        return "AnnData raw expression"
    return "AnnData.X"


def _expression_matrix_code(expression: ExpressionSource, object_name: str) -> str:
    if expression.kind == "layer":
        return f"{object_name}.layers[{expression.layer_name!r}]"
    if expression.kind == "raw":
        return f"{object_name}.raw.X"
    return f"{object_name}.X"


def _expression_var(adata: AnnData, expression: ExpressionSource) -> object:
    return adata.raw.var if expression.kind == "raw" else adata.var


def _copy_preserving_duplicate_observation_ids(adata: AnnData) -> AnnData:
    with python_warnings.catch_warnings():
        python_warnings.filterwarnings(
            "ignore",
            message="Observation names are not unique.*",
            category=UserWarning,
        )
        return adata.copy()


def _comma_separated_metrics(value: str) -> list[str]:
    metrics = list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))
    if not metrics:
        raise ValueError("At least one QC metric is required.")
    return metrics


def _mad_outlier_statistics(
    values: object,
    *,
    nmads: float,
    direction: Literal["both", "upper", "lower"],
    scale_mad: bool,
) -> tuple[object, dict[str, int | float | None]]:
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


def _group_positions(
    adata: AnnData,
    group_key: str,
    *,
    report_context: str = "grouped reports",
) -> list[tuple[str, object]]:
    science = dependencies.require_scientific_dependencies()
    if not group_key:
        return [("all_cells", science.np.arange(adata.n_obs, dtype=int))]
    _string_group_labels(adata.obs[group_key], group_key=group_key, report_context=report_context)
    grouped = adata.obs.groupby(group_key, observed=True, sort=False).indices
    return [(str(label), science.np.asarray(positions, dtype=int)) for label, positions in grouped.items()]


def _string_group_labels(
    values: object,
    *,
    group_key: str,
    report_context: str,
) -> list[str]:
    science = dependencies.require_scientific_dependencies()
    original = list(science.pd.unique(values))
    rendered = [str(label) for label in original]
    if len(set(rendered)) != len(rendered):
        raise ValueError(
            f"Sample/group labels in {group_key!r} collide when represented in {report_context}; "
            "use unambiguous labels with one consistent dtype."
        )
    return rendered


def _mad_outliers_code(parameters: dict[str, object]) -> str:
    return dedent(
        f"""
        import numpy as np
        import pandas as pd
        from scipy.stats import median_abs_deviation


        def mark_mad_outliers(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("MAD outlier marking requires at least one cell and one gene.")
            output = adata.copy()
            metrics = list({parameters["metrics"]!r})
            group_key = {parameters["batch_key"]!r}
            output_column = {parameters["output_column"]!r}
            missing_metrics = [metric for metric in metrics if metric not in output.obs]
            if missing_metrics:
                raise ValueError(f"QC metrics not found in obs: {{missing_metrics}}")
            for metric in metrics:
                if not pd.api.types.is_numeric_dtype(output.obs[metric].dtype):
                    raise ValueError(f"QC metric must be numeric: {{metric!r}}")
            if group_key:
                if group_key not in output.obs:
                    raise ValueError(f"QC Sample/group column not found in obs: {{group_key!r}}")
                if output.obs[group_key].isna().any():
                    raise ValueError(f"QC Sample/group column contains missing labels: {{group_key!r}}")
                original_labels = list(pd.unique(output.obs[group_key]))
                rendered_labels = [str(label) for label in original_labels]
                if len(set(rendered_labels)) != len(rendered_labels):
                    raise ValueError(
                        f"Sample/group labels in {{group_key!r}} collide when represented in reports; "
                        "use unambiguous labels with one consistent dtype."
                    )
                groups = [
                    np.asarray(positions, dtype=int)
                    for positions in output.obs.groupby(group_key, observed=True, sort=False).indices.values()
                ]
            else:
                groups = [np.arange(output.n_obs, dtype=int)]

            flags = np.zeros(output.n_obs, dtype=bool)
            for metric in metrics:
                values = np.asarray(output.obs[metric], dtype=float)
                metric_flags = np.zeros(output.n_obs, dtype=bool)
                for positions in groups:
                    group_values = values[positions]
                    finite = np.isfinite(group_values)
                    if int(finite.sum()) < {parameters["minimum_group_size"]!r}:
                        continue
                    finite_values = group_values[finite]
                    median = float(np.median(finite_values))
                    mad = float(
                        median_abs_deviation(
                            finite_values,
                            scale={"'normal'" if parameters["scale_mad"] else "1.0"},
                        )
                    )
                    lower = -np.inf if {parameters["direction"]!r} == "upper" else median - {
                        parameters["nmads"]!r
                    } * mad
                    upper = np.inf if {parameters["direction"]!r} == "lower" else median + {
                        parameters["nmads"]!r
                    } * mad
                    metric_flags[positions] = finite & ((group_values < lower) | (group_values > upper))
                flags |= metric_flags
            output.obs[output_column] = flags
            return output
        """
    )


class OpenBioSingleCellMarkMADOutliers(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkMADOutliers",
            display_name="Mark MAD Outliers",
            category=CATEGORY,
            description="Mark a union of metric outliers, optionally within each Sample/capture.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input(
                    "metrics",
                    default="total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
                ),
                io.String.Input("batch_key", default="sample"),
                io.Float.Input("nmads", default=5.0, min=0.1, max=100.0, step=0.1),
                io.Combo.Input("direction", options=["both", "upper", "lower"], default="both"),
                io.String.Input("output_column", default="outlier", advanced=True),
                io.Boolean.Input("scale_mad", default=False, advanced=True),
                io.Int.Input("minimum_group_size", default=3, min=1, max=MAX_RANDOM_SEED, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        metrics: str = "total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
        batch_key: str = "sample",
        nmads: float = 5.0,
        direction: Literal["both", "upper", "lower"] = "both",
        output_column: str = "outlier",
        scale_mad: bool = False,
        minimum_group_size: int = 3,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        if cells == 0 or genes == 0:
            raise ValueError("MAD outlier marking requires at least one cell and one gene.")
        if not bool(science.np.isfinite(nmads)) or nmads <= 0:
            raise ValueError("MAD multiplier must be greater than zero.")
        if direction not in {"both", "upper", "lower"}:
            raise ValueError(f"Unsupported MAD direction: {direction!r}")
        if minimum_group_size < 1:
            raise ValueError("MAD minimum group size must be at least one finite observation.")
        batch_key = batch_key.strip()
        output_column = output_column.strip()
        if not output_column:
            raise ValueError("MAD output column cannot be empty.")

        output = adata.copy()
        selected_metrics = _comma_separated_metrics(metrics)
        missing_metrics = [metric for metric in selected_metrics if metric not in output.obs]
        if missing_metrics:
            raise ValueError(f"QC metrics not found in obs: {missing_metrics}")
        if batch_key and batch_key not in output.obs:
            raise ValueError(f"QC Sample/group column not found in obs: {batch_key!r}")
        if batch_key and bool(output.obs[batch_key].isna().any()):
            raise ValueError(f"QC Sample/group column contains missing labels: {batch_key!r}")
        for metric in selected_metrics:
            if not science.pd.api.types.is_numeric_dtype(output.obs[metric].dtype):
                raise ValueError(f"QC metric must be numeric: {metric!r}")

        warnings = []
        if output_column in output.obs:
            warnings.append(f"Existing observation column {output_column!r} was overwritten.")
        if not batch_key:
            warnings.append("MAD thresholds were estimated globally because no Sample/group key was supplied.")
        groups = _group_positions(output, batch_key, report_context="grouped MAD reports")
        flags = science.np.zeros(cells, dtype=bool)
        metric_statistics = []
        for metric in selected_metrics:
            values = science.np.asarray(output.obs[metric], dtype=float)
            metric_flags = science.np.zeros(cells, dtype=bool)
            group_statistics = []
            for group_label, positions in groups:
                group_values = values[positions]
                finite_count = int(science.np.isfinite(group_values).sum())
                if finite_count < minimum_group_size:
                    missing = int(group_values.size - finite_count)
                    group_statistics.append(
                        {
                            "group": group_label,
                            "n": int(group_values.size),
                            "finite": finite_count,
                            "missing": missing,
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
                        direction=direction,
                        scale_mad=scale_mad,
                )
                metric_flags[positions] = group_flags
                status = "zero_mad" if statistics["mad"] == 0.0 else "evaluated"
                group_statistics.append({"group": group_label, "status": status, **statistics})
                if status == "zero_mad":
                    warnings.append(
                        f"Metric {metric!r}, group {group_label!r} had zero MAD; the effective threshold equals "
                        "the group median."
                    )
            flags |= metric_flags
            metric_statistics.append(
                {
                    "metric": metric,
                    "flagged": int(metric_flags.sum()),
                    "missing": int((~science.np.isfinite(values)).sum()),
                    "groups": group_statistics,
                }
            )
        output.obs[output_column] = flags

        parameters = {
            "metrics": selected_metrics,
            "batch_key": batch_key,
            "nmads": nmads,
            "direction": direction,
            "output_column": output_column,
            "scale_mad": scale_mad,
            "mad_scale": "normal" if scale_mad else "raw",
            "minimum_group_size": minimum_group_size,
        }
        finish_adata(output, "mark_mad_outliers", parameters, cells, genes, started_at, warnings=warnings)

        group_union = [
            {
                "group": group_label,
                "cells": int(len(positions)),
                "flagged": int(flags[positions].sum()),
                "flagged_percent": float(flags[positions].mean() * 100.0),
            }
            for group_label, positions in groups
        ]
        flagged_cells = int(flags.sum())
        methods = (
            f"MAD outliers were marked for {', '.join(selected_metrics)} using a {direction}-direction rule at "
            f"{nmads:g} MADs. Thresholds were estimated "
            f"{'within ' + batch_key if batch_key else 'globally'} with "
            f"{'normal-consistent' if scale_mad else 'raw'} MAD scaling; metric flags were combined by union."
        )
        results = (
            f"{flagged_cells:,} of {cells:,} cells ({flagged_cells * 100.0 / cells:.1f}%) were marked in "
            f"obs[{output_column!r}]. No cells were removed."
        )
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellMarkMADOutliers",
            title="MAD outlier marking summary",
            operation="mark_mad_outliers_report",
            methods=methods,
            results=results,
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
            parameters=parameters,
            references=(SCIPY_MAD_REFERENCE, LEYS_MAD_REFERENCE),
            software_packages=("scipy", "numpy", "pandas", "anndata"),
            warnings=warnings,
            limitations=(
                "MAD cutoffs are dataset- and Sample-dependent QC flags, not universal definitions of low-quality cells.",
                "Metrics with different biological directions should be evaluated in separate node invocations.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            code=_mad_outliers_code(parameters),
        )
        return io.NodeOutput(output, report, code)


def _scrublet_group_payloads(scrublet_uns: object) -> dict[str, dict[str, object]]:
    if not isinstance(scrublet_uns, dict):
        raise RuntimeError("Scanpy Scrublet did not produce an uns['scrublet'] result mapping.")
    batches = scrublet_uns.get("batches")
    if batches is None:
        return {"all_cells": scrublet_uns}
    if not isinstance(batches, dict):
        raise RuntimeError("Scanpy Scrublet produced an invalid batched result mapping.")
    return {str(label): payload for label, payload in batches.items() if isinstance(payload, dict)}


def _scrublet_thresholds(scrublet_uns: object, expected_groups: list[str]) -> dict[str, float]:
    science = dependencies.require_scientific_dependencies()
    payloads = _scrublet_group_payloads(scrublet_uns)
    if set(payloads) != set(expected_groups):
        missing = sorted(set(expected_groups) - set(payloads))
        raise RuntimeError(f"Scanpy Scrublet did not return results for every Sample/group: {missing}")
    thresholds = {}
    missing_thresholds = []
    for label, payload in payloads.items():
        threshold = payload.get("threshold")
        try:
            numeric = float(threshold)
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


def _scrublet_code(expression: ExpressionSource, parameters: dict[str, object]) -> str:
    matrix_code = _expression_matrix_code(expression, "adata")
    var_code = "adata.raw.var.copy()" if expression.kind == "raw" else "adata.var.copy()"
    layer_argument = repr(expression.layer_name) if expression.kind == "layer" else "None"
    use_raw = expression.kind == "raw"
    threshold = parameters["threshold"] if parameters["threshold_mode"] == "manual" else None
    n_neighbors = parameters["n_neighbors"] or None
    return dedent(
        f"""
        import copy
        import warnings

        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
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
            thresholds = {{}}
            missing_thresholds = []
            for label, payload in payloads.items():
                try:
                    value = float(payload.get("threshold"))
                except (TypeError, ValueError):
                    missing_thresholds.append(label)
                    continue
                if not np.isfinite(value):
                    missing_thresholds.append(label)
                    continue
                thresholds[label] = value
            if missing_thresholds:
                raise RuntimeError(
                    "Scrublet did not determine a usable threshold for Sample/group(s) "
                    f"{{missing_thresholds}}. Inspect the simulated score distribution and rerun with a manual threshold."
                )
            return thresholds


        def run_scrublet(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("Scrublet requires at least one cell and one gene.")
            layer = {layer_argument}
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
                original_labels = list(pd.unique(adata.obs[group_key]))
                rendered_labels = [str(label) for label in original_labels]
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

            matrix = {matrix_code}
            if matrix.shape[1] <= {parameters["n_prin_comps"]!r}:
                raise ValueError(
                    f"Scrublet n_prin_comps={parameters['n_prin_comps']!r} must be smaller than the selected "
                    f"source feature count ({{matrix.shape[1]}})."
                )
            values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
            values = np.asarray(values)
            if not np.isfinite(values).all():
                raise ValueError("Selected Scrublet source contains non-finite expression values.")
            if (values < 0).any():
                raise ValueError("Selected Scrublet source contains negative expression values and cannot be used as counts.")
            if values.size == 0 or not (values > 0).any():
                raise ValueError("Selected Scrublet source contains no positive expression values.")
            work_obs = adata.obs.copy()
            if not adata.obs_names.is_unique:
                work_obs.index = [f"__openbio_scrublet_row_{{index}}" for index in range(adata.n_obs)]
            work = ad.AnnData(
                X=matrix.copy(),
                obs=work_obs,
                var={var_code},
            )
            try:
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
                    n_neighbors={n_neighbors!r},
                    threshold={threshold!r},
                    verbose=False,
                    copy=False,
                    random_state={parameters["random_seed"]!r},
                )
            except ValueError as error:
                if "n_components" in str(error) or "n_prin_comps" in str(error):
                    raise ValueError(
                        "Scrublet principal-component dimension is too large after its internal cell, gene, "
                        "and variable-gene filtering; reduce n_prin_comps or provide a larger Sample/capture."
                    ) from error
                raise
            _validated_scrublet_thresholds(work.uns.get("scrublet"), expected_groups)
            if "doublet_score" not in work.obs or "predicted_doublet" not in work.obs:
                raise RuntimeError("Scanpy Scrublet did not produce its documented observation columns.")
            scores = np.asarray(work.obs["doublet_score"], dtype=float)
            predictions = work.obs["predicted_doublet"]
            if not np.isfinite(scores).all():
                raise RuntimeError("Scanpy Scrublet produced non-finite observed doublet scores.")
            if not pd.api.types.is_bool_dtype(predictions.dtype) or predictions.isna().any():
                raise RuntimeError("Scanpy Scrublet predicted_doublet output is not a complete boolean column.")

            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="Observation names are not unique.*",
                    category=UserWarning,
                )
                output = adata.copy()
            output.obs["doublet_score"] = scores
            output.obs["predicted_doublet"] = predictions.to_numpy(dtype=bool)
            output.uns["scrublet"] = copy.deepcopy(work.uns["scrublet"])
            return output
        """
    )


class OpenBioSingleCellScrublet(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Scrublet count expression source",
        default="X",
        include_raw=True,
        layer_input_id="source_layer",
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellScrublet",
            display_name="Run Scrublet",
            category=CATEGORY,
            description="Score and predict capture-level doublets from unnormalized counts without filtering cells.",
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("batch_key", default="sample"),
                io.Int.Input("random_seed", default=123, min=0, max=MAX_RANDOM_SEED, advanced=True),
                io.Float.Input("expected_doublet_rate", default=0.05, min=0.0001, max=0.9999, step=0.01),
                io.Combo.Input("threshold_mode", options=["automatic", "manual"], default="automatic"),
                io.Float.Input("threshold", default=0.25, step=0.01),
                io.Float.Input("sim_doublet_ratio", default=2.0, min=0.1, max=100.0, step=0.1, advanced=True),
                io.Int.Input("n_prin_comps", default=30, min=1, max=10_000, advanced=True),
                io.Int.Input("n_neighbors", default=0, min=0, max=1_000_000, advanced=True),
                cls.EXPRESSION_SOURCE.input(),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        batch_key: str = "sample",
        random_seed: int = 123,
        expected_doublet_rate: float = 0.05,
        threshold_mode: Literal["automatic", "manual"] = "automatic",
        threshold: float = 0.25,
        sim_doublet_ratio: float = 2.0,
        n_prin_comps: int = 30,
        n_neighbors: int = 0,
        source: DynamicExpressionSource | None = None,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        if cells == 0 or genes == 0:
            raise ValueError("Scrublet requires at least one cell and one gene.")
        if not 0.0 < expected_doublet_rate < 1.0:
            raise ValueError("Scrublet expected doublet rate must be between zero and one.")
        if threshold_mode not in {"automatic", "manual"}:
            raise ValueError(f"Unsupported Scrublet threshold mode: {threshold_mode!r}")
        manual_threshold = None
        if threshold_mode == "manual":
            try:
                manual_threshold = float(threshold)
            except (TypeError, ValueError) as error:
                raise ValueError("Scrublet manual threshold must be finite.") from error
            if not bool(science.np.isfinite(manual_threshold)):
                raise ValueError("Scrublet manual threshold must be finite.")
        try:
            sim_doublet_ratio = float(sim_doublet_ratio)
        except (TypeError, ValueError) as error:
            raise ValueError("Scrublet simulated-doublet ratio must be finite and greater than zero.") from error
        if not bool(science.np.isfinite(sim_doublet_ratio)) or sim_doublet_ratio <= 0:
            raise ValueError("Scrublet simulated-doublet ratio must be finite and greater than zero.")
        if n_prin_comps < 1:
            raise ValueError("Scrublet principal components must be at least one.")
        if n_neighbors < 0:
            raise ValueError("Scrublet neighbors cannot be negative; zero selects Scanpy's automatic value.")

        batch_key = batch_key.strip()
        warnings = []
        if manual_threshold is not None and not 0.0 <= manual_threshold <= 1.0:
            warnings.append(
                f"Manual Scrublet threshold {manual_threshold:g} lies outside the conventional [0, 1] score "
                "interval; the explicit extreme threshold was retained."
            )
        duplicate_cell_ids = not adata.obs_names.is_unique
        if duplicate_cell_ids:
            warnings.append(
                "Observation identifiers are not unique. Scrublet was run with temporary positional identifiers; "
                "the original row order and duplicate Cell IDs were preserved in the output."
            )
        overwritten_scrublet_fields = [
            field
            for field, present in (
                ("obs['doublet_score']", "doublet_score" in adata.obs),
                ("obs['predicted_doublet']", "predicted_doublet" in adata.obs),
                ("uns['scrublet']", "scrublet" in adata.uns),
            )
            if present
        ]
        if overwritten_scrublet_fields:
            warnings.append(
                "Existing Scrublet result fields were overwritten: " + ", ".join(overwritten_scrublet_fields) + "."
            )
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
            minimum_cells = n_prin_comps + 1
            too_small = {str(label): int(size) for label, size in group_sizes.items() if int(size) < minimum_cells}
            if too_small:
                raise ValueError(
                    f"Scrublet n_prin_comps={n_prin_comps} requires at least {minimum_cells} cells in every "
                    f"Sample/group: {too_small}"
                )
        else:
            minimum_cells = n_prin_comps + 1
            if cells < minimum_cells:
                raise ValueError(
                    f"Scrublet n_prin_comps={n_prin_comps} requires at least {minimum_cells} cells."
                )
            expected_groups = ["all_cells"]
            warnings.append(
                "Scrublet was run globally because no Sample/group key was supplied; merged captures can distort "
                "simulated doublet composition."
            )

        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        source_label = _expression_source_label(expression)
        matrix = expression.matrix(adata)
        if matrix.shape[1] <= n_prin_comps:
            raise ValueError(
                f"Scrublet n_prin_comps={n_prin_comps} must be smaller than the selected source feature count "
                f"({matrix.shape[1]})."
            )
        warnings.extend(
            validate_count_expression(
                matrix,
                source_label=source_label,
                require_positive=True,
            )
        )
        work_obs = adata.obs.copy()
        if duplicate_cell_ids:
            work_obs.index = [f"__openbio_scrublet_row_{index}" for index in range(cells)]
        work = science.ad.AnnData(
            X=matrix.copy(),
            obs=work_obs,
            var=_expression_var(adata, expression).copy(),
        )
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
                    "Scrublet principal-component dimension is too large after its internal cell, gene, and "
                    "variable-gene filtering; reduce n_prin_comps or provide a larger Sample/capture."
                ) from error
            raise
        scrublet_uns = work.uns.get("scrublet")
        thresholds = _scrublet_thresholds(scrublet_uns, expected_groups)
        if "doublet_score" not in work.obs or "predicted_doublet" not in work.obs:
            raise RuntimeError("Scanpy Scrublet did not produce its documented observation columns.")
        scores = science.np.asarray(work.obs["doublet_score"], dtype=float)
        predictions = work.obs["predicted_doublet"]
        if not bool(science.np.isfinite(scores).all()):
            raise RuntimeError("Scanpy Scrublet produced non-finite observed doublet scores.")
        if not science.pd.api.types.is_bool_dtype(predictions.dtype) or bool(predictions.isna().any()):
            raise RuntimeError("Scanpy Scrublet predicted_doublet output is not a complete boolean column.")

        output = _copy_preserving_duplicate_observation_ids(adata)
        output.obs["doublet_score"] = scores
        output.obs["predicted_doublet"] = predictions.to_numpy(dtype=bool)
        output.uns["scrublet"] = copy.deepcopy(work.uns["scrublet"])
        parameters = {
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
            warnings.append(
                "Scrublet predicted no doublets; verify the score distribution and effective threshold before filtering."
            )
        finish_adata(
            output,
            "scrublet",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
            warnings=warnings,
        )

        groups = _group_positions(output, batch_key, report_context="Scrublet reports")
        group_payloads = _scrublet_group_payloads(scrublet_uns)
        per_group = []
        for group_label, positions in groups:
            group_predictions = science.np.asarray(predictions.iloc[positions], dtype=bool)
            group_scores = scores[positions]
            effective_parameters = group_payloads[group_label].get("parameters", {})
            per_group.append(
                {
                    "group": group_label,
                    "cells": int(len(positions)),
                    "predicted_doublets": int(group_predictions.sum()),
                    "predicted_percent": float(group_predictions.mean() * 100.0),
                    "doublet_score": summarize_numeric(group_scores),
                    "threshold": thresholds[group_label],
                    "effective_parameters": effective_parameters,
                }
            )
        predicted_percent = float(predicted_count * 100.0 / cells)
        methods = (
            f"Scanpy Scrublet scored {cells:,} observed transcriptomes from {source_label} and simulated doublets "
            f"{'within ' + batch_key if batch_key else 'globally'}. Predictions used "
            f"{'automatically selected thresholds' if threshold_mode == 'automatic' else f'a manual threshold of {manual_threshold:g}'}."
        )
        results = (
            f"Scrublet predicted {predicted_count:,} of {cells:,} cells ({predicted_percent:.1f}%) as doublets. "
            "Cells were annotated but not removed."
        )
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellScrublet",
            title="Scrublet doublet prediction summary",
            operation="scrublet_report",
            methods=methods,
            results=results,
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
            parameters=parameters,
            references=(SCRUBLET_REFERENCE, SCANPY_REFERENCE),
            software_packages=(*CORRECTION_SOFTWARE_PACKAGES, "scikit-image"),
            warnings=warnings,
            limitations=(
                "Scrublet predictions are model-based QC evidence, not ground-truth doublet labels; inspect score "
                "separation and embedding localization before filtering.",
                "Homotypic doublets can be difficult to detect, and expected rates can differ between capture libraries.",
                "Run captures separately when their cell-state compositions or loading rates differ materially.",
            ),
            input_cells=cells,
            input_genes=int(work.n_vars),
            started_at=started_at,
            code=_scrublet_code(expression, parameters),
        )
        return io.NodeOutput(output, report, code)


def _filter_doublets_code(prediction_column: str) -> str:
    return dedent(
        f"""
        import numpy as np
        import pandas as pd


        def filter_predicted_doublets(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("Doublet filtering requires at least one cell and one gene.")
            column = {prediction_column!r}
            if column not in adata.obs:
                raise ValueError(f"Doublet prediction column not found in obs: {{column!r}}")
            predictions = adata.obs[column]
            if not pd.api.types.is_bool_dtype(predictions.dtype):
                raise ValueError(f"Doublet prediction column must have boolean dtype: {{column!r}}")
            if predictions.isna().any():
                raise ValueError(f"Doublet prediction column contains missing values: {{column!r}}")
            keep = ~predictions.to_numpy(dtype=bool)
            return adata[keep].copy()
        """
    )


class OpenBioSingleCellFilterDoublets(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellFilterDoublets",
            display_name="Filter Predicted Doublets",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("prediction_column", default="predicted_doublet", advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        prediction_column: str = "predicted_doublet",
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        if cells == 0 or genes == 0:
            raise ValueError("Doublet filtering requires at least one cell and one gene.")
        prediction_column = prediction_column.strip()
        if not prediction_column:
            raise ValueError("Doublet prediction column cannot be empty.")
        if prediction_column not in adata.obs:
            raise ValueError(f"Doublet prediction column not found in obs: {prediction_column!r}")
        predictions = adata.obs[prediction_column]
        if not science.pd.api.types.is_bool_dtype(predictions.dtype):
            raise ValueError(f"Doublet prediction column must have boolean dtype: {prediction_column!r}")
        if bool(predictions.isna().any()):
            raise ValueError(f"Doublet prediction column contains missing values: {prediction_column!r}")
        predicted_values = predictions.to_numpy(dtype=bool)
        keep = ~predicted_values
        retained_cells = int(keep.sum())
        output = adata[keep].copy()
        parameters = {"prediction_column": prediction_column}
        predicted_count = int(predicted_values.sum())
        warnings = []
        if retained_cells == 0:
            warnings.append(
                "Every input cell was marked for removal. The empty result is preserved as the explicit filtering "
                "outcome and will not be usable by most downstream analyses."
            )
        elif predicted_count == 0:
            warnings.append("No cells were marked as doublets; the output retains every input cell.")
        finish_adata(output, "filter_doublets", parameters, cells, genes, started_at, warnings=warnings)

        removed_percent = float(predicted_count * 100.0 / cells)
        retained_percent = float(retained_cells * 100.0 / cells)
        methods = (
            f"Cells with boolean True in obs[{prediction_column!r}] were removed by inverse-mask subsetting. "
            "The prediction column was consumed without coercing strings, numbers, or missing values."
        )
        results = (
            f"Removed {predicted_count:,} of {cells:,} cells ({removed_percent:.1f}%) marked as predicted doublets "
            f"and retained {retained_cells:,} cells ({retained_percent:.1f}%); all {genes:,} genes were preserved."
        )
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellFilterDoublets",
            title="Predicted-doublet filtering summary",
            operation="filter_doublets_report",
            methods=methods,
            results=results,
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
            parameters=parameters,
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
        return io.NodeOutput(output, report, code)


CORRECTION_NODE_CLASSES = [
    OpenBioSingleCellMarkMADOutliers,
    OpenBioSingleCellScrublet,
    OpenBioSingleCellFilterDoublets,
]


__all__ = ["CORRECTION_NODE_CLASSES"]
