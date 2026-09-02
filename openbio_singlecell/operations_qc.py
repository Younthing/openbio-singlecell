from __future__ import annotations

import inspect
import math
import time
from collections.abc import Mapping
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import (
    figure_to_png,
    finish_adata,
    make_plot_result,
    matrix_totals_and_nonzero,
    validate_count_expression,
)
from .expression_source import (
    _CALCULATE_QC_SPEC,
    _FILTER_CELLS_SPEC,
    _FILTER_GENES_SPEC,
    _QC_PLOTS_SPEC,
    DynamicExpressionSource,
    ExpressionSource,
)
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

if TYPE_CHECKING:
    from anndata import AnnData


PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
CALCULATE_QC_EXPRESSION_SOURCE = _CALCULATE_QC_SPEC
FILTER_CELLS_EXPRESSION_SOURCE = _FILTER_CELLS_SPEC
FILTER_GENES_EXPRESSION_SOURCE = _FILTER_GENES_SPEC
QC_PLOTS_EXPRESSION_SOURCE = _QC_PLOTS_SPEC
# ponytail: static PNGs are readable up to 40 groups; add pagination or an interactive renderer for larger sets.
_QC_PLOT_MAX_GROUPS = 40
QC_SOFTWARE_PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy", "matplotlib")
SCANPY_REFERENCE = AnalysisReference(
    citation=(
        "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. "
        "Genome Biology. 2018;19:15."
    ),
    doi="10.1186/s13059-017-1382-0",
    url="https://doi.org/10.1186/s13059-017-1382-0",
    kind="software",
)
SCATER_REFERENCE = AnalysisReference(
    citation=(
        "McCarthy DJ, Campbell KR, Lun ATL, Wills QF. Scater: pre-processing, quality control, normalization "
        "and visualization of single-cell RNA-seq data in R. Bioinformatics. 2017;33(8):1179-1186."
    ),
    doi="10.1093/bioinformatics/btw777",
    url="https://doi.org/10.1093/bioinformatics/btw777",
    kind="method",
)
MITO_QC_REFERENCE = AnalysisReference(
    citation=(
        "Osorio D, Cai JJ. Systematic determination of the mitochondrial proportion in human and mice tissues "
        "for single-cell RNA-sequencing data quality control. Bioinformatics. 2021;37(7):963-967."
    ),
    doi="10.1093/bioinformatics/btaa751",
    url="https://doi.org/10.1093/bioinformatics/btaa751",
    kind="practice",
)
QC_LIMITATIONS = (
    "QC distributions and filtering thresholds are assay-, tissue-, and Sample-dependent; review them by Sample "
    "rather than treating the configured values as universal cutoffs.",
)


def _required_bool(parameters: dict[str, JSONValue], name: str) -> bool:
    value = parameters[name]
    if type(value) is not bool:
        raise TypeError(f"{name} must be a boolean.")
    return value


def _required_int(parameters: dict[str, JSONValue], name: str) -> int:
    value = parameters[name]
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer.")
    return value


def _required_number(parameters: dict[str, JSONValue], name: str) -> float:
    value = parameters[name]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a number.")
    return float(value)


def _required_string(parameters: dict[str, JSONValue], name: str) -> str:
    value = parameters[name]
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string.")
    return value


def _source_parameter(parameters: dict[str, JSONValue]) -> DynamicExpressionSource | None:
    value = parameters["source"]
    if value is not None and not isinstance(value, dict):
        raise TypeError("source must be a DynamicCombo object.")
    return value


def _expression_source_label(expression: ExpressionSource) -> str:
    if expression.kind == "layer":
        return f"AnnData layer {expression.layer_name!r}"
    if expression.kind == "raw":
        return "AnnData raw expression"
    return "AnnData.X"


def _expression_matrix_code(expression: ExpressionSource, object_name: str = "adata") -> str:
    if expression.kind == "layer":
        return f"{object_name}.layers[{expression.layer_name!r}]"
    if expression.kind == "raw":
        return f"{object_name}.raw.X"
    return f"{object_name}.X"


def _expression_var(adata: AnnData, expression: ExpressionSource) -> object:
    return adata.raw.var if expression.kind == "raw" else adata.var


def _validate_qc_log1p_domain(matrix: object, qc_masks: list[object]) -> None:
    science = dependencies.require_scientific_dependencies()
    row_totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
    column_totals, _ = matrix_totals_and_nonzero(matrix, axis=0)
    column_means = science.np.asarray(column_totals, dtype=float) / int(matrix.shape[0])
    if bool((science.np.asarray(row_totals, dtype=float) <= -1.0).any()) or bool((column_means <= -1.0).any()):
        raise ValueError("QC log1p annotations are undefined because a selected expression sum or mean is <= -1.")
    for mask in qc_masks:
        subset_totals, _ = matrix_totals_and_nonzero(matrix[:, mask], axis=1)
        if bool((science.np.asarray(subset_totals, dtype=float) <= -1.0).any()):
            raise ValueError("QC log1p annotations are undefined because a QC-feature expression sum is <= -1.")


def _map_raw_variable_metrics(
    work: AnnData,
    output: AnnData,
    metric_columns: list[str],
    warnings: list[str],
) -> int:
    if not work.var_names.is_unique:
        warnings.append(
            "Raw feature identifiers are not unique; observation QC metrics were calculated, but variable-level "
            "QC metrics could not be mapped safely to current var."
        )
        return 0
    mapped = output.var_names.isin(work.var_names)
    aligned = work.var.loc[:, metric_columns].reindex(output.var_names)
    for column in metric_columns:
        output.var[column] = aligned[column].to_numpy()
    mapped_count = int(mapped.sum())
    if mapped_count != output.n_vars:
        warnings.append(
            f"Raw variable-level QC metrics mapped to {mapped_count} of {output.n_vars} current features; "
            "unmatched current features carry missing variable metrics."
        )
    return mapped_count


def _percentage(numerator: int, denominator: int) -> float:
    return float(numerator * 100.0 / denominator) if denominator else 0.0


def _display_number(value: int | float | None) -> str:
    return "unavailable" if value is None else f"{value:,.3g}"


def _validate_bounds(minimum: int, maximum: int, *, label: str) -> None:
    if minimum < 0 or maximum < 0:
        raise ValueError(f"{label} bounds cannot be negative.")
    if minimum > 0 and maximum > 0 and minimum > maximum:
        raise ValueError(f"{label} minimum cannot exceed its maximum.")


def _resolve_numeric_threshold(value: float, enabled: bool, *, label: str) -> float | None:
    if not enabled:
        return None
    if not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number when enabled.")
    return value


def _validate_numeric_interval(minimum: float | None, maximum: float | None, *, label: str) -> None:
    if minimum is not None and maximum is not None and minimum > maximum:
        raise ValueError(f"{label} minimum cannot exceed its maximum.")


def _comma_separated_prefixes(value: str, name: str) -> tuple[str, ...]:
    prefixes = tuple(dict.fromkeys(item.strip().upper() for item in value.split(",") if item.strip()))
    if not prefixes:
        raise ValueError(f"{name} must contain at least one comma-separated prefix.")
    return prefixes


def _percent_top_values(value: str) -> tuple[int, ...]:
    values = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            position = int(item)
        except ValueError as error:
            raise ValueError(f"percent_top values must be positive integers; received {item!r}.") from error
        if position <= 0:
            raise ValueError("percent_top values must be greater than zero.")
        if position not in values:
            values.append(position)
    return tuple(values)


def _calculate_qc_code(expression: ExpressionSource, parameters: dict[str, object]) -> str:
    layer_argument = repr(expression.layer_name) if expression.kind == "layer" else "None"
    use_raw = expression.kind == "raw"
    return dedent(
        f"""
        import numpy as np
        import scanpy as sc
        from scipy import sparse


        def calculate_qc_metrics(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("QC metric calculation requires at least one cell and one gene.")
            layer = {layer_argument}
            if layer is not None and layer not in adata.layers:
                raise ValueError(f"Expression layer not found: {{layer!r}}")
            if {use_raw!r} and adata.raw is None:
                raise ValueError("Raw expression is unavailable.")
            work = adata.raw.to_adata() if {use_raw!r} else adata
            matrix = work.X if layer is None else work.layers[layer]
            values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
            values = np.asarray(values)
            if values.size and not np.isfinite(values).all():
                raise ValueError("Selected expression source contains non-finite expression values.")
            gene_names = work.var_names.astype(str)
            upper_gene_names = gene_names.str.upper()
            work.var["mt"] = upper_gene_names.str.startswith({str(parameters["mitochondrial_prefix"]).upper()!r})

            qc_vars = ["mt"]
            ribosomal_prefixes = tuple({parameters["ribosomal_prefixes"]!r})
            if {parameters["include_ribosomal"]!r}:
                work.var["ribo"] = upper_gene_names.str.startswith(ribosomal_prefixes)
                qc_vars.append("ribo")
            if {parameters["include_hemoglobin"]!r}:
                work.var["hb"] = gene_names.str.contains(
                    {parameters["hemoglobin_pattern"]!r}, case=False, regex=True, na=False
                )
                qc_vars.append("hb")

            current_names = adata.var_names.astype(str)
            current_upper = current_names.str.upper()
            adata.var["mt"] = current_upper.str.startswith({str(parameters["mitochondrial_prefix"]).upper()!r})
            if {parameters["include_ribosomal"]!r}:
                adata.var["ribo"] = current_upper.str.startswith(ribosomal_prefixes)
            if {parameters["include_hemoglobin"]!r}:
                adata.var["hb"] = current_names.str.contains(
                    {parameters["hemoglobin_pattern"]!r}, case=False, regex=True, na=False
                )

            if {parameters["log1p"]!r}:
                row_totals = np.asarray(matrix.sum(axis=1)).ravel()
                column_means = np.asarray(matrix.sum(axis=0)).ravel() / matrix.shape[0]
                if (row_totals <= -1).any() or (column_means <= -1).any():
                    raise ValueError(
                        "QC log1p annotations are undefined because a selected expression sum or mean is <= -1."
                    )
                for qc_var in qc_vars:
                    subset_totals = np.asarray(matrix[:, np.asarray(work.var[qc_var], dtype=bool)].sum(axis=1)).ravel()
                    if (subset_totals <= -1).any():
                        raise ValueError(
                            "QC log1p annotations are undefined because a QC-feature expression sum is <= -1."
                        )

            requested_percent_top = list({parameters["percent_top"]!r})
            supported_percent_top = [rank for rank in requested_percent_top if rank <= work.n_vars]
            before_metric_columns = set(work.var.columns)
            sc.pp.calculate_qc_metrics(
                work,
                qc_vars=qc_vars,
                percent_top=supported_percent_top or None,
                layer=None if {use_raw!r} else {layer_argument},
                use_raw=False,
                log1p={parameters["log1p"]!r},
                inplace=True,
            )
            for qc_var in qc_vars:
                if not bool(work.var[qc_var].any()):
                    for prefix in ("total_counts_", "log1p_total_counts_", "pct_counts_"):
                        column = f"{{prefix}}{{qc_var}}"
                        if column in work.obs:
                            del work.obs[column]
            if {use_raw!r}:
                adata.obs = work.obs.copy()
                metric_columns = [column for column in work.var.columns if column not in before_metric_columns]
                if work.var_names.is_unique:
                    aligned = work.var.loc[:, metric_columns].reindex(adata.var_names)
                    for column in metric_columns:
                        adata.var[column] = aligned[column].to_numpy()
            if supported_percent_top != requested_percent_top:
                total_counts = np.asarray(adata.obs["total_counts"], dtype=float)
                all_features_percentage = np.full(total_counts.shape, np.nan, dtype=float)
                all_features_percentage[total_counts != 0] = 100.0
                for rank in requested_percent_top:
                    if rank > work.n_vars:
                        adata.obs[f"pct_counts_in_top_{{rank}}_genes"] = all_features_percentage
            return adata
        """
    )


def calculate_qc_owned(
    adata: Any,
    *,
    include_ribosomal: bool = True,
    include_hemoglobin: bool = True,
    mitochondrial_prefix: str = "MT-",
    ribosomal_prefixes: str = "RPS,RPL",
    hemoglobin_pattern: str = r"^HB[^(P)]",
    percent_top: str = "20,50,100,200,500",
    log1p: bool = True,
    source: DynamicExpressionSource | None = None,
) -> tuple[Any, Any, str]:
    parameters: dict[str, JSONValue] = {
        "include_ribosomal": include_ribosomal,
        "include_hemoglobin": include_hemoglobin,
        "mitochondrial_prefix": mitochondrial_prefix,
        "ribosomal_prefixes": ribosomal_prefixes,
        "hemoglobin_pattern": hemoglobin_pattern,
        "percent_top": percent_top,
        "log1p": log1p,
        "source": source,
    }
    include_ribosomal = _required_bool(parameters, "include_ribosomal")
    include_hemoglobin = _required_bool(parameters, "include_hemoglobin")
    mitochondrial_prefix = _required_string(parameters, "mitochondrial_prefix")
    ribosomal_prefixes = _required_string(parameters, "ribosomal_prefixes")
    hemoglobin_pattern = _required_string(parameters, "hemoglobin_pattern")
    percent_top = _required_string(parameters, "percent_top")
    log1p = _required_bool(parameters, "log1p")
    if not mitochondrial_prefix.strip():
        raise ValueError("Mitochondrial gene prefix cannot be empty.")
    if include_hemoglobin and not hemoglobin_pattern.strip():
        raise ValueError("Hemoglobin gene pattern cannot be empty when hemoglobin QC is enabled.")

    cells, genes = int(adata.n_obs), int(adata.n_vars)
    if cells == 0 or genes == 0:
        raise ValueError("QC metric calculation requires at least one cell and one gene.")
    expression = CALCULATE_QC_EXPRESSION_SOURCE.resolve(adata, _source_parameter(parameters))
    source_label = _expression_source_label(expression)
    matrix = expression.matrix(adata)
    warnings = validate_count_expression(matrix, source_label=source_label, require_nonnegative=False)

    started_at = time.perf_counter()
    # Raw can have a broader feature axis, so it needs a temporary AnnData materialization for Scanpy metrics.
    work = adata.raw.to_adata() if expression.kind == "raw" else adata
    gene_names = work.var_names.astype(str)
    upper_gene_names = gene_names.str.upper()
    work.var["mt"] = upper_gene_names.str.startswith(mitochondrial_prefix.strip().upper())

    qc_vars = ["mt"]
    parsed_ribosomal_prefixes: tuple[str, ...] = ()
    if include_ribosomal:
        parsed_ribosomal_prefixes = _comma_separated_prefixes(ribosomal_prefixes, "ribosomal_prefixes")
        work.var["ribo"] = upper_gene_names.str.startswith(parsed_ribosomal_prefixes)
        qc_vars.append("ribo")
    if include_hemoglobin:
        work.var["hb"] = gene_names.str.contains(hemoglobin_pattern, case=False, regex=True, na=False)
        qc_vars.append("hb")

    current_gene_names = adata.var_names.astype(str)
    current_upper_gene_names = current_gene_names.str.upper()
    adata.var["mt"] = current_upper_gene_names.str.startswith(mitochondrial_prefix.strip().upper())
    if include_ribosomal:
        adata.var["ribo"] = current_upper_gene_names.str.startswith(parsed_ribosomal_prefixes)
    if include_hemoglobin:
        adata.var["hb"] = current_gene_names.str.contains(hemoglobin_pattern, case=False, regex=True, na=False)

    requested_percent_top = _percent_top_values(percent_top)
    source_genes = int(matrix.shape[1])
    supported_percent_top = tuple(position for position in requested_percent_top if position <= source_genes)
    for qc_var in qc_vars:
        if not bool(work.var[qc_var].any()):
            warning = (
                f"No genes matched mitochondrial prefix {mitochondrial_prefix!r}."
                if qc_var == "mt"
                else f"No genes matched the {qc_var!r} QC annotation."
            )
            warnings.append(warning)
    if log1p:
        science = dependencies.require_scientific_dependencies()
        _validate_qc_log1p_domain(
            matrix,
            [science.np.asarray(work.var[qc_var], dtype=bool) for qc_var in qc_vars],
        )
    before_metric_columns = set(work.var.columns)
    science = dependencies.require_scientific_dependencies()
    science.sc.pp.calculate_qc_metrics(
        work,
        qc_vars=qc_vars,
        percent_top=supported_percent_top or None,
        layer=None if expression.kind == "raw" else expression.scanpy_layer,
        use_raw=False,
        log1p=log1p,
        inplace=True,
    )
    for qc_var in qc_vars:
        if not bool(work.var[qc_var].any()):
            for prefix in ("total_counts_", "log1p_total_counts_", "pct_counts_"):
                column = f"{prefix}{qc_var}"
                if column in work.obs:
                    del work.obs[column]

    mapped_variable_metrics = genes
    if expression.kind == "raw":
        adata.obs = work.obs.copy()
        metric_columns = [column for column in work.var.columns if column not in before_metric_columns]
        mapped_variable_metrics = _map_raw_variable_metrics(work, adata, metric_columns, warnings)

    undefined_all_feature_percent_cells = 0
    if requested_percent_top and supported_percent_top != requested_percent_top:
        total_counts = science.np.asarray(adata.obs["total_counts"], dtype=float)
        all_features_percentage = science.np.full(total_counts.shape, science.np.nan, dtype=float)
        nonzero_totals = total_counts != 0
        all_features_percentage[nonzero_totals] = 100.0
        undefined_all_feature_percent_cells = int((~nonzero_totals).sum())
        for position in requested_percent_top:
            if position > source_genes:
                adata.obs[f"pct_counts_in_top_{position}_genes"] = all_features_percentage
        warnings.append(
            "Some percent_top positions exceed the number of genes; those percentages include all available genes."
        )
        if undefined_all_feature_percent_cells:
            warnings.append(
                f"All-feature percentages were undefined for {undefined_all_feature_percent_cells:,} cell(s) "
                "with zero selected-source expression sum and were reported as missing rather than 0%."
            )

    report_parameters = {
        **expression.parameters(),
        "include_ribosomal": include_ribosomal,
        "include_hemoglobin": include_hemoglobin,
        "mitochondrial_prefix": mitochondrial_prefix.strip(),
        "ribosomal_prefixes": list(parsed_ribosomal_prefixes),
        "hemoglobin_pattern": hemoglobin_pattern,
        "percent_top": list(requested_percent_top),
        "log1p": log1p,
    }
    finish_adata(adata, "calculate_qc", report_parameters, cells, genes, started_at, warnings=warnings)

    qc_distributions = {
        "total_counts": summarize_numeric(adata.obs["total_counts"]),
        "detected_genes": summarize_numeric(adata.obs["n_genes_by_counts"]),
    }
    if "pct_counts_mt" in adata.obs:
        qc_distributions["mitochondrial_percent"] = summarize_numeric(adata.obs["pct_counts_mt"])
    matched_genes = {qc_var: int(work.var[qc_var].sum()) for qc_var in qc_vars}
    total_median = qc_distributions["total_counts"]["median"]
    detected_median = qc_distributions["detected_genes"]["median"]
    mitochondrial_median = qc_distributions.get("mitochondrial_percent", {}).get("median")
    mitochondrial_result = (
        f" Median mitochondrial expression was {_display_number(mitochondrial_median)}%."
        if mitochondrial_median is not None
        else ""
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCalculateQC",
        title="QC metric calculation summary",
        operation="calculate_qc_report",
        methods=(
            f"Quality-control metrics were calculated for {cells:,} cells and {source_genes:,} selected-source "
            f"features from {source_label} with scanpy.pp.calculate_qc_metrics on the worker's private AnnData. "
            "Configured mitochondrial, ribosomal, and hemoglobin gene annotations were supplied as QC feature "
            "sets when enabled."
        ),
        results=(
            f"QC metrics were added without filtering cells or genes. Median total expression was "
            f"{_display_number(total_median)} and the median number of detected genes was "
            f"{_display_number(detected_median)}.{mitochondrial_result}"
        ),
        key_results={
            "input_cells": cells,
            "input_genes": genes,
            "source_features": source_genes,
            "current_features_with_mapped_variable_metrics": mapped_variable_metrics,
            "matched_qc_genes": matched_genes,
            "requested_percent_top": list(requested_percent_top),
            "effective_percent_top": list(supported_percent_top),
            "undefined_all_feature_percent_cells": undefined_all_feature_percent_cells,
            "distributions": qc_distributions,
        },
        parameters=report_parameters,
        references=(SCANPY_REFERENCE, SCATER_REFERENCE, MITO_QC_REFERENCE),
        software_packages=QC_SOFTWARE_PACKAGES,
        warnings=warnings,
        limitations=QC_LIMITATIONS,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_calculate_qc_code(expression, report_parameters),
    )
    return adata, report, code


@register_operation("openbio.node.calculateqc")
def calculate_qc(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    operation = "Calculate QC"
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(
        parameters,
        {
            "include_ribosomal",
            "include_hemoglobin",
            "mitochondrial_prefix",
            "ribosomal_prefixes",
            "hemoglobin_pattern",
            "percent_top",
            "log1p",
            "source",
        },
        operation=operation,
    )
    output, report, code = calculate_qc_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_anndata_output(context, output))


def _filter_cells_code(expression: ExpressionSource, parameters: dict[str, object]) -> str:
    matrix_code = _expression_matrix_code(expression)
    layer_argument = repr(expression.layer_name) if expression.kind == "layer" else "None"
    use_raw = expression.kind == "raw"
    return dedent(
        f"""
        import numpy as np
        from scipy import sparse


        def filter_cells(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("Cell filtering requires at least one cell and one gene.")
            layer = {layer_argument}
            if layer is not None and layer not in adata.layers:
                raise ValueError(f"Expression layer not found: {{layer!r}}")
            if {use_raw!r} and adata.raw is None:
                raise ValueError("Raw expression is unavailable.")
            matrix = {matrix_code}
            values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
            values = np.asarray(values)
            if values.size and not np.isfinite(values).all():
                raise ValueError("Selected expression source contains non-finite expression values.")
            if sparse.issparse(matrix):
                total_expression = np.asarray(matrix.sum(axis=1)).ravel()
                detected_genes = np.asarray(matrix.count_nonzero(axis=1)).ravel()
            else:
                matrix = np.asarray(matrix)
                total_expression = matrix.sum(axis=1)
                detected_genes = np.count_nonzero(matrix, axis=1)

            keep = np.ones(adata.n_obs, dtype=bool)
            if {parameters["min_genes"]!r} > 0:
                keep &= detected_genes >= {parameters["min_genes"]!r}
            if {parameters["max_genes"]!r} > 0:
                keep &= detected_genes <= {parameters["max_genes"]!r}
            if {parameters["enable_min_counts"]!r}:
                keep &= total_expression >= {parameters["min_counts"]!r}
            if {parameters["enable_max_counts"]!r}:
                keep &= total_expression <= {parameters["max_counts"]!r}
            if {parameters["enable_max_pct_mito"]!r}:
                column = {parameters["mito_column"]!r}
                if column not in adata.obs:
                    raise ValueError(f"Cell mitochondrial percentage column not found in obs: {{column!r}}")
                mitochondrial_percent = np.asarray(adata.obs[column], dtype=float)
                if not np.isfinite(mitochondrial_percent).all():
                    raise ValueError(
                        f"Cell mitochondrial percentage column contains missing or non-finite values: {{column!r}}"
                    )
                keep &= mitochondrial_percent <= {parameters["max_pct_mito"]!r}
            return adata[keep].copy()
        """
    )


def filter_cells_owned(
    adata: Any,
    *,
    min_genes: int = 0,
    max_genes: int = 0,
    min_counts: float = 0.0,
    max_counts: float = 0.0,
    max_pct_mito: float = 0.0,
    mito_column: str = "pct_counts_mt",
    source: DynamicExpressionSource | None = None,
    enable_min_counts: bool = False,
    enable_max_counts: bool = False,
    enable_max_pct_mito: bool = False,
) -> tuple[Any, Any, str]:
    parameters: dict[str, JSONValue] = {
        "min_genes": min_genes,
        "max_genes": max_genes,
        "min_counts": min_counts,
        "max_counts": max_counts,
        "max_pct_mito": max_pct_mito,
        "mito_column": mito_column,
        "source": source,
        "enable_min_counts": enable_min_counts,
        "enable_max_counts": enable_max_counts,
        "enable_max_pct_mito": enable_max_pct_mito,
    }
    min_genes = _required_int(parameters, "min_genes")
    max_genes = _required_int(parameters, "max_genes")
    min_counts = _required_number(parameters, "min_counts")
    max_counts = _required_number(parameters, "max_counts")
    max_pct_mito = _required_number(parameters, "max_pct_mito")
    mito_column = _required_string(parameters, "mito_column").strip()
    enable_min_counts = _required_bool(parameters, "enable_min_counts")
    enable_max_counts = _required_bool(parameters, "enable_max_counts")
    enable_max_pct_mito = _required_bool(parameters, "enable_max_pct_mito")

    _validate_bounds(min_genes, max_genes, label="Detected-gene")
    minimum_expression = _resolve_numeric_threshold(
        min_counts,
        enable_min_counts,
        label="Minimum total-expression threshold",
    )
    maximum_expression = _resolve_numeric_threshold(
        max_counts,
        enable_max_counts,
        label="Maximum total-expression threshold",
    )
    maximum_mitochondrial_percent = _resolve_numeric_threshold(
        max_pct_mito,
        enable_max_pct_mito,
        label="Maximum mitochondrial percentage threshold",
    )
    _validate_numeric_interval(minimum_expression, maximum_expression, label="Total-expression")
    if maximum_mitochondrial_percent is not None and not 0.0 <= maximum_mitochondrial_percent <= 100.0:
        raise ValueError("Cell mitochondrial percentage threshold must be between 0 and 100 when enabled.")

    science = dependencies.require_scientific_dependencies()
    mito_values = None
    if maximum_mitochondrial_percent is not None:
        if mito_column not in adata.obs:
            raise ValueError(f"Cell mitochondrial percentage column not found in obs: {mito_column!r}")
        try:
            mito_values = science.np.asarray(adata.obs[mito_column], dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"Cell mitochondrial percentage column must contain numeric values: {mito_column!r}"
            ) from error
        if not bool(science.np.isfinite(mito_values).all()):
            raise ValueError(
                f"Cell mitochondrial percentage column contains missing or non-finite values: {mito_column!r}"
            )

    cells, genes = int(adata.n_obs), int(adata.n_vars)
    if cells == 0 or genes == 0:
        raise ValueError("Cell filtering requires at least one cell and one gene.")
    expression = FILTER_CELLS_EXPRESSION_SOURCE.resolve(adata, _source_parameter(parameters))
    source_label = _expression_source_label(expression)
    matrix = expression.matrix(adata)
    warnings = validate_count_expression(matrix, source_label=source_label, require_nonnegative=False)
    started_at = time.perf_counter()
    counts, detected = matrix_totals_and_nonzero(matrix, axis=1)
    mask = science.np.ones(cells, dtype=bool)
    if min_genes > 0:
        mask &= detected >= min_genes
    if max_genes > 0:
        mask &= detected <= max_genes
    if minimum_expression is not None:
        mask &= counts >= minimum_expression
    if maximum_expression is not None:
        mask &= counts <= maximum_expression
    if maximum_mitochondrial_percent is not None:
        mask &= mito_values <= maximum_mitochondrial_percent
    retained_cells = int(mask.sum())
    # AnnData slicing returns a view; materialization is scientifically required before publishing an independent file.
    output = adata[mask].copy()
    report_parameters = {
        **expression.parameters(),
        "min_genes": min_genes,
        "max_genes": max_genes,
        "min_counts": minimum_expression,
        "max_counts": maximum_expression,
        "max_pct_mito": maximum_mitochondrial_percent,
        "mito_column": mito_column,
        "enable_min_counts": minimum_expression is not None,
        "enable_max_counts": maximum_expression is not None,
        "enable_max_pct_mito": maximum_mitochondrial_percent is not None,
    }
    active_filters = [
        name
        for name, enabled in (
            ("min_genes", min_genes > 0),
            ("max_genes", max_genes > 0),
            ("min_counts", minimum_expression is not None),
            ("max_counts", maximum_expression is not None),
            ("max_pct_mito", maximum_mitochondrial_percent is not None),
        )
        if enabled
    ]
    if retained_cells == 0:
        warnings.append(
            "The configured criteria removed every cell. The empty result is preserved as the explicit "
            "filtering outcome and will not be usable by most downstream analyses."
        )
    elif not active_filters:
        warnings.append("No cell-filter criterion was enabled; the output contains every input cell.")
    elif retained_cells == cells:
        warnings.append("All cells passed the configured filtering criteria.")
    finish_adata(output, "filter_cells", report_parameters, cells, genes, started_at, warnings=warnings)

    removed_cells = cells - retained_cells
    retention_percent = _percentage(retained_cells, cells)
    zero_total_expression_cells = int((counts == 0).sum())
    mitochondrial_metric = {
        "column": mito_column,
        "available": mito_column in adata.obs,
        "used_for_filtering": maximum_mitochondrial_percent is not None,
        "source": f"obs[{mito_column!r}]" if maximum_mitochondrial_percent is not None else None,
        "distribution": summarize_numeric(mito_values) if mito_values is not None else None,
    }
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellFilterCells",
        title="Cell filtering summary",
        operation="filter_cells_report",
        methods=(
            f"Cells were filtered by retaining the intersection of enabled criteria. Total expression and detected "
            f"genes were calculated from {source_label}; when enabled, mitochondrial percentage was read from "
            f"obs[{mito_column!r}]. The selected AnnData view was materialized once for the output artifact."
        ),
        results=(
            f"Of {cells:,} input cells, {retained_cells:,} ({retention_percent:.1f}%) were retained and "
            f"{removed_cells:,} were removed; all {genes:,} genes were preserved. "
            f"The input contained {zero_total_expression_cells:,} cells with zero total expression in {source_label}."
        ),
        key_results={
            "input_cells": cells,
            "retained_cells": retained_cells,
            "removed_cells": removed_cells,
            "retention_percent": retention_percent,
            "genes": genes,
            "source_features": int(matrix.shape[1]),
            "active_filters": active_filters,
            "criterion_sources": {
                "total_expression": source_label,
                "detected_genes": source_label,
                "mitochondrial_percent": mitochondrial_metric["source"],
            },
            "zero_total_expression_cells": zero_total_expression_cells,
            "mitochondrial_metric": mitochondrial_metric,
            "input_distributions": {
                "total_expression": summarize_numeric(counts),
                "detected_genes": summarize_numeric(detected),
            },
        },
        parameters=report_parameters,
        references=(SCANPY_REFERENCE, SCATER_REFERENCE, MITO_QC_REFERENCE),
        software_packages=QC_SOFTWARE_PACKAGES,
        warnings=warnings,
        limitations=QC_LIMITATIONS,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_filter_cells_code(expression, report_parameters),
    )
    return output, report, code


@register_operation("openbio.node.filtercells")
def filter_cells(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    operation = "Filter Cells"
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(
        parameters,
        {
            "min_genes",
            "max_genes",
            "min_counts",
            "max_counts",
            "max_pct_mito",
            "mito_column",
            "source",
            "enable_min_counts",
            "enable_max_counts",
            "enable_max_pct_mito",
        },
        operation=operation,
    )
    output, report, code = filter_cells_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_anndata_output(context, output))


def _gene_filter_statistics(
    adata: AnnData,
    expression: ExpressionSource,
) -> tuple[object, object, object, list[str]]:
    science = dependencies.require_scientific_dependencies()
    matrix = expression.matrix(adata)
    source_counts, source_detected = matrix_totals_and_nonzero(matrix, axis=0)
    if expression.kind != "raw" or adata.raw.var_names.equals(adata.var_names):
        return (
            science.np.asarray(source_counts, dtype=float),
            science.np.asarray(source_detected, dtype=float),
            science.np.ones(adata.n_vars, dtype=bool),
            [],
        )
    if not adata.raw.var_names.is_unique:
        raise ValueError("Gene filtering cannot map a broader Raw feature axis with duplicate Raw feature identifiers.")
    source_positions = {str(name): index for index, name in enumerate(adata.raw.var_names)}
    counts = science.np.full(adata.n_vars, science.np.nan, dtype=float)
    detected = science.np.full(adata.n_vars, science.np.nan, dtype=float)
    evaluated = science.np.zeros(adata.n_vars, dtype=bool)
    for current_index, name in enumerate(adata.var_names):
        source_index = source_positions.get(str(name))
        if source_index is None:
            continue
        counts[current_index] = source_counts[source_index]
        detected[current_index] = source_detected[source_index]
        evaluated[current_index] = True
    unevaluated = int((~evaluated).sum())
    warnings = []
    if unevaluated:
        warnings.append(
            f"The selected Raw source did not contain {unevaluated} current feature identifier(s); those features "
            "were retained because no selected-source filtering evidence was available."
        )
    return counts, detected, evaluated, warnings


def _filter_genes_code(expression: ExpressionSource, parameters: dict[str, object]) -> str:
    matrix_code = _expression_matrix_code(expression)
    layer_argument = repr(expression.layer_name) if expression.kind == "layer" else "None"
    use_raw = expression.kind == "raw"
    return dedent(
        f"""
        import numpy as np
        from scipy import sparse


        def filter_genes(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("Gene filtering requires at least one cell and one gene.")
            layer = {layer_argument}
            if layer is not None and layer not in adata.layers:
                raise ValueError(f"Expression layer not found: {{layer!r}}")
            if {use_raw!r} and adata.raw is None:
                raise ValueError("Raw expression is unavailable.")
            matrix = {matrix_code}
            values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
            values = np.asarray(values)
            if values.size and not np.isfinite(values).all():
                raise ValueError("Selected expression source contains non-finite expression values.")
            if sparse.issparse(matrix):
                source_total_expression = np.asarray(matrix.sum(axis=0)).ravel()
                source_detected_cells = np.asarray(matrix.count_nonzero(axis=0)).ravel()
            else:
                matrix = np.asarray(matrix)
                source_total_expression = matrix.sum(axis=0)
                source_detected_cells = np.count_nonzero(matrix, axis=0)

            if {use_raw!r} and not adata.raw.var_names.equals(adata.var_names):
                if not adata.raw.var_names.is_unique:
                    raise ValueError(
                        "Gene filtering cannot map a broader Raw feature axis with duplicate Raw feature identifiers."
                    )
                positions = {{str(name): index for index, name in enumerate(adata.raw.var_names)}}
                total_expression = np.full(adata.n_vars, np.nan, dtype=float)
                detected_cells = np.full(adata.n_vars, np.nan, dtype=float)
                evaluated = np.zeros(adata.n_vars, dtype=bool)
                for current_index, name in enumerate(adata.var_names):
                    source_index = positions.get(str(name))
                    if source_index is not None:
                        total_expression[current_index] = source_total_expression[source_index]
                        detected_cells[current_index] = source_detected_cells[source_index]
                        evaluated[current_index] = True
            else:
                total_expression = np.asarray(source_total_expression, dtype=float)
                detected_cells = np.asarray(source_detected_cells, dtype=float)
                evaluated = np.ones(adata.n_vars, dtype=bool)

            keep = np.ones(adata.n_vars, dtype=bool)
            if {parameters["min_cells"]!r} > 0:
                keep[evaluated] &= detected_cells[evaluated] >= {parameters["min_cells"]!r}
            if {parameters["max_cells"]!r} > 0:
                keep[evaluated] &= detected_cells[evaluated] <= {parameters["max_cells"]!r}
            if {parameters["enable_min_counts"]!r}:
                keep[evaluated] &= total_expression[evaluated] >= {parameters["min_counts"]!r}
            if {parameters["enable_max_counts"]!r}:
                keep[evaluated] &= total_expression[evaluated] <= {parameters["max_counts"]!r}
            return adata[:, keep].copy()
        """
    )


def filter_genes_owned(
    adata: Any,
    *,
    min_cells: int = 0,
    max_cells: int = 0,
    min_counts: float = 0.0,
    max_counts: float = 0.0,
    source: DynamicExpressionSource | None = None,
    enable_min_counts: bool = False,
    enable_max_counts: bool = False,
) -> tuple[Any, Any, str]:
    parameters: dict[str, JSONValue] = {
        "min_cells": min_cells,
        "max_cells": max_cells,
        "min_counts": min_counts,
        "max_counts": max_counts,
        "source": source,
        "enable_min_counts": enable_min_counts,
        "enable_max_counts": enable_max_counts,
    }
    min_cells = _required_int(parameters, "min_cells")
    max_cells = _required_int(parameters, "max_cells")
    min_counts = _required_number(parameters, "min_counts")
    max_counts = _required_number(parameters, "max_counts")
    enable_min_counts = _required_bool(parameters, "enable_min_counts")
    enable_max_counts = _required_bool(parameters, "enable_max_counts")

    _validate_bounds(min_cells, max_cells, label="Detected-cell")
    minimum_expression = _resolve_numeric_threshold(
        min_counts,
        enable_min_counts,
        label="Minimum total-expression threshold",
    )
    maximum_expression = _resolve_numeric_threshold(
        max_counts,
        enable_max_counts,
        label="Maximum total-expression threshold",
    )
    _validate_numeric_interval(minimum_expression, maximum_expression, label="Total-expression")
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    if cells == 0 or genes == 0:
        raise ValueError("Gene filtering requires at least one cell and one gene.")
    expression = FILTER_GENES_EXPRESSION_SOURCE.resolve(adata, _source_parameter(parameters))
    source_label = _expression_source_label(expression)
    matrix = expression.matrix(adata)
    warnings = validate_count_expression(matrix, source_label=source_label, require_nonnegative=False)
    started_at = time.perf_counter()
    counts, detected, evaluated, alignment_warnings = _gene_filter_statistics(adata, expression)
    warnings.extend(alignment_warnings)
    science = dependencies.require_scientific_dependencies()
    mask = science.np.ones(genes, dtype=bool)
    if min_cells > 0:
        mask[evaluated] &= detected[evaluated] >= min_cells
    if max_cells > 0:
        mask[evaluated] &= detected[evaluated] <= max_cells
    if minimum_expression is not None:
        mask[evaluated] &= counts[evaluated] >= minimum_expression
    if maximum_expression is not None:
        mask[evaluated] &= counts[evaluated] <= maximum_expression
    retained_genes = int(mask.sum())
    # AnnData slicing returns a view; materialization is scientifically required before publishing an independent file.
    output = adata[:, mask].copy()
    report_parameters = {
        **expression.parameters(),
        "min_cells": min_cells,
        "max_cells": max_cells,
        "min_counts": minimum_expression,
        "max_counts": maximum_expression,
        "enable_min_counts": minimum_expression is not None,
        "enable_max_counts": maximum_expression is not None,
    }
    active_filters = [
        name
        for name, enabled in (
            ("min_cells", min_cells > 0),
            ("max_cells", max_cells > 0),
            ("min_counts", minimum_expression is not None),
            ("max_counts", maximum_expression is not None),
        )
        if enabled
    ]
    if retained_genes == 0:
        warnings.append(
            "The configured criteria removed every gene. The empty result is preserved as the explicit "
            "filtering outcome and will not be usable by most downstream analyses."
        )
    elif not active_filters:
        warnings.append("No gene-filter criterion was enabled; the output contains every input gene.")
    elif retained_genes == genes:
        warnings.append("All genes passed the configured filtering criteria.")
    finish_adata(output, "filter_genes", report_parameters, cells, genes, started_at, warnings=warnings)

    removed_genes = genes - retained_genes
    retention_percent = _percentage(retained_genes, genes)
    zero_total_expression_genes = int(((counts == 0) & evaluated).sum())
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellFilterGenes",
        title="Gene filtering summary",
        operation="filter_genes_report",
        methods=(
            f"Genes were filtered from {source_label} by retaining the intersection of all enabled total-expression "
            "and detected-cell criteria. This low-information filter is distinct from highly variable gene selection."
        ),
        results=(
            f"Of {genes:,} input genes, {retained_genes:,} ({retention_percent:.1f}%) were retained and "
            f"{removed_genes:,} were removed; all {cells:,} cells were preserved."
        ),
        key_results={
            "input_genes": genes,
            "retained_genes": retained_genes,
            "removed_genes": removed_genes,
            "retention_percent": retention_percent,
            "cells": cells,
            "source_features": int(matrix.shape[1]),
            "evaluated_current_genes": int(evaluated.sum()),
            "unevaluated_current_genes": int((~evaluated).sum()),
            "active_filters": active_filters,
            "zero_total_expression_genes": zero_total_expression_genes,
            "input_distributions": {
                "total_expression": summarize_numeric(counts),
                "detected_cells": summarize_numeric(detected),
            },
        },
        parameters=report_parameters,
        references=(SCANPY_REFERENCE, SCATER_REFERENCE),
        software_packages=QC_SOFTWARE_PACKAGES,
        warnings=warnings,
        limitations=(
            "Low-detection filtering can remove markers of rare populations; retain a full-gene Raw snapshot "
            "when downstream analyses require recoverable post-QC counts.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_filter_genes_code(expression, report_parameters),
    )
    return output, report, code


@register_operation("openbio.node.filtergenes")
def filter_genes(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    operation = "Filter Genes"
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(
        parameters,
        {
            "min_cells",
            "max_cells",
            "min_counts",
            "max_counts",
            "source",
            "enable_min_counts",
            "enable_max_counts",
        },
        operation=operation,
    )
    output, report, code = filter_genes_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_anndata_output(context, output))


def _qc_plot_metric_values(matrix, source_var, numpy, scipy_sparse):
    values = matrix.data if scipy_sparse.issparse(matrix) else numpy.asarray(matrix).ravel()
    values = numpy.asarray(values)
    if values.size and not numpy.isfinite(values).all():
        raise ValueError("Selected expression source contains non-finite expression values.")
    if scipy_sparse.issparse(matrix):
        total_expression = numpy.asarray(matrix.sum(axis=1)).ravel()
        detected_genes = numpy.asarray(matrix.count_nonzero(axis=1)).ravel()
    else:
        dense = numpy.asarray(matrix)
        total_expression = dense.sum(axis=1)
        detected_genes = numpy.count_nonzero(dense, axis=1)
    total_expression = numpy.asarray(total_expression, dtype=float)
    detected_genes = numpy.asarray(detected_genes, dtype=float)

    mitochondrial_percent = None
    if "mt" in source_var and bool(numpy.asarray(source_var["mt"], dtype=bool).any()):
        mitochondrial_mask = numpy.asarray(source_var["mt"], dtype=bool)
        mitochondrial_matrix = matrix[:, mitochondrial_mask]
        if scipy_sparse.issparse(mitochondrial_matrix):
            mitochondrial_expression = numpy.asarray(mitochondrial_matrix.sum(axis=1)).ravel()
        else:
            mitochondrial_expression = numpy.asarray(mitochondrial_matrix).sum(axis=1)
        mitochondrial_percent = numpy.divide(
            mitochondrial_expression * 100.0,
            total_expression,
            out=numpy.full(total_expression.shape, numpy.nan, dtype=float),
            where=total_expression != 0,
        )
    return total_expression, detected_genes, mitochondrial_percent


def _derive_qc_plot_metrics(
    matrix: Any,
    source_var: Any,
    source_label: str,
) -> tuple[Any, Any, Any | None, dict[str, str], list[str]]:
    science = dependencies.require_scientific_dependencies()
    warnings = validate_count_expression(matrix, source_label=source_label, require_nonnegative=False)
    total_expression, detected_genes, mitochondrial_percent = _qc_plot_metric_values(
        matrix,
        source_var,
        science.np,
        science.sparse,
    )
    metric_sources = {"total_expression": source_label, "detected_genes": source_label}
    if mitochondrial_percent is not None:
        metric_sources["mitochondrial_percent"] = source_label
        metric_sources["mitochondrial_panel_total_expression"] = source_label
    else:
        warnings.append("Mitochondrial annotations are unavailable; the mitochondrial panel was omitted.")
    return total_expression, detected_genes, mitochondrial_percent, metric_sources, warnings


def _resolve_qc_plot_view(value: Mapping[str, object] | None) -> dict[str, str]:
    if value is None:
        return {"view": "overview"}
    if not isinstance(value, Mapping):
        raise TypeError("QC plot view must be a DynamicCombo value.")
    mode = value.get("view")
    if not isinstance(mode, str):
        raise TypeError("QC plot view selection must be a string.")
    mode = mode.strip()
    if mode not in {"overview", "grouped"}:
        raise ValueError(f"Unsupported QC plot view: {mode!r}.")
    allowed = {"view", "groupby"} if mode == "grouped" else {"view"}
    unknown = sorted(str(key) for key in value if key not in allowed)
    if unknown:
        raise ValueError(f"QC plot {mode} view received inactive or unknown parameters: {unknown!r}.")
    if mode == "overview":
        return {"view": mode}
    groupby = value.get("groupby", "sample")
    if not isinstance(groupby, str):
        raise TypeError("QC plot groupby must be a string.")
    groupby = groupby.strip()
    if not groupby:
        raise ValueError("QC plot groupby cannot be empty.")
    return {"view": mode, "groupby": groupby}


def _qc_plot_group_details(obs, groupby, pandas, numpy):
    if groupby not in obs:
        raise ValueError(f"QC plot groupby column not found in obs: {groupby!r}.")
    groups = obs[groupby]
    if isinstance(groups.dtype, pandas.CategoricalDtype):
        codes = groups.cat.codes.to_numpy(copy=True)
        categories = list(groups.cat.categories)
    elif (
        pandas.api.types.is_object_dtype(groups.dtype)
        or isinstance(groups.dtype, pandas.StringDtype)
        or pandas.api.types.is_bool_dtype(groups.dtype)
    ):
        categories = []
        code_by_identity = {}
        codes = numpy.empty(len(groups), dtype=int)
        for row, value in enumerate(groups.tolist()):
            missing = pandas.isna(value)
            if not isinstance(missing, (bool, numpy.bool_)) or bool(missing):
                raise ValueError(f"QC plot obs[{groupby!r}] contains missing group labels.")
            identity = (type(value).__module__, type(value).__qualname__, repr(value))
            if identity not in code_by_identity:
                code_by_identity[identity] = len(categories)
                categories.append(value)
            codes[row] = code_by_identity[identity]
    else:
        raise TypeError(f"QC plot obs[{groupby!r}] must be categorical, boolean, or string-like.")
    if bool((codes < 0).any()):
        raise ValueError(f"QC plot obs[{groupby!r}] contains missing group labels.")
    if isinstance(groups.dtype, pandas.CategoricalDtype):
        declared_counts = numpy.bincount(codes, minlength=len(categories))
        represented = numpy.flatnonzero(declared_counts)
        remapped = numpy.full(len(categories), -1, dtype=int)
        remapped[represented] = numpy.arange(represented.size)
        codes = remapped[codes]
        categories = [categories[index] for index in represented]
        counts = declared_counts[represented]
    else:
        counts = numpy.bincount(codes, minlength=len(categories))
    labels = []
    for original in categories:
        value = original.item() if isinstance(original, numpy.generic) else original
        if isinstance(value, float) and not numpy.isfinite(value):
            raise ValueError(f"QC plot obs[{groupby!r}] contains non-finite group labels.")
        if not isinstance(value, (str, bool, int, float)):
            raise TypeError(f"QC plot obs[{groupby!r}] group labels must be scalar strings or finite numbers.")
        labels.append(str(value).strip())
    if any(not label for label in labels):
        raise ValueError(f"QC plot obs[{groupby!r}] contains blank group labels.")
    if len(set(labels)) != len(labels):
        raise ValueError(f"QC plot obs[{groupby!r}] contains labels that collide after string conversion.")
    if len(labels) > _QC_PLOT_MAX_GROUPS:
        raise ValueError(
            f"QC grouped plot supports at most {_QC_PLOT_MAX_GROUPS} visible groups; "
            f"obs[{groupby!r}] contains {len(labels)}."
        )
    ordered_rows = numpy.argsort(codes, kind="stable")
    positions = numpy.split(ordered_rows, numpy.cumsum(counts)[:-1])
    return labels, positions, {label: int(count) for label, count in zip(labels, counts, strict=True)}


def _qc_plot_figure(
    total_expression,
    detected_genes,
    mitochondrial_percent,
    view,
    groupby,
    group_labels,
    group_positions,
    Figure,
    numpy,
):
    if not numpy.isfinite(total_expression).any() or not numpy.isfinite(detected_genes).any():
        raise ValueError("QC plots require at least one finite total-expression and detected-gene value.")
    if view == "overview":
        figure = Figure(figsize=(10, 7), constrained_layout=True)
        axes = figure.subplots(2, 2)
        axes[0, 0].hist(total_expression[numpy.isfinite(total_expression)], bins=40, color="#246bfe")
        axes[0, 0].set_title("Total expression per cell")
        axes[0, 1].hist(detected_genes[numpy.isfinite(detected_genes)], bins=40, color="#17a673")
        axes[0, 1].set_title("Genes detected per cell")
        finite = numpy.isfinite(total_expression) & numpy.isfinite(detected_genes)
        axes[1, 0].scatter(total_expression[finite], detected_genes[finite], s=8, alpha=0.6, color="#6f4bf2")
        axes[1, 0].set_xlabel("Total expression")
        axes[1, 0].set_ylabel("Genes detected")
        if mitochondrial_percent is None:
            axes[1, 1].text(0.5, 0.5, "Mitochondrial metric unavailable", ha="center", va="center")
            axes[1, 1].set_axis_off()
        else:
            finite = numpy.isfinite(total_expression) & numpy.isfinite(mitochondrial_percent)
            axes[1, 1].scatter(
                total_expression[finite],
                mitochondrial_percent[finite],
                s=8,
                alpha=0.6,
                color="#e45d3a",
            )
            axes[1, 1].set_xlabel("Total expression")
            axes[1, 1].set_ylabel("Mitochondrial expression (%)")
        return figure

    metrics = [
        ("Total expression", total_expression, "#246bfe"),
        ("Genes detected", detected_genes, "#17a673"),
    ]
    if mitochondrial_percent is not None:
        metrics.append(("Mitochondrial expression (%)", mitochondrial_percent, "#e45d3a"))
    figure = Figure(figsize=(max(8, 4 * len(metrics)), 4.5), constrained_layout=True)
    axes = figure.subplots(1, len(metrics), squeeze=False)[0]
    tick_positions = numpy.arange(1, len(group_labels) + 1)
    for axis, (label, values, color) in zip(axes, metrics, strict=True):
        datasets = []
        positions = []
        for position, selected in zip(tick_positions, group_positions, strict=True):
            finite = values[selected]
            finite = finite[numpy.isfinite(finite)]
            if finite.size:
                datasets.append(finite)
                positions.append(position)
        if datasets:
            parts = axis.violinplot(datasets, positions=positions, showmedians=True)
            for body in parts["bodies"]:
                body.set_facecolor(color)
                body.set_edgecolor(color)
                body.set_alpha(0.75)
        axis.set_title(f"{label} by {groupby}")
        axis.set_xlabel(groupby)
        axis.set_ylabel(label)
        axis.set_xticks(tick_positions, group_labels, rotation=30, ha="right")
    return figure


def _qc_plots_code(expression: ExpressionSource, view: Mapping[str, str]) -> str:
    matrix_code = _expression_matrix_code(expression)
    layer_argument = repr(expression.layer_name) if expression.kind == "layer" else "None"
    use_raw = expression.kind == "raw"
    mode = view["view"]
    groupby = view.get("groupby")
    helper_source = "\n\n".join(
        dedent(inspect.getsource(helper)).strip()
        for helper in (_qc_plot_metric_values, _qc_plot_group_details, _qc_plot_figure)
    )
    implementation = dedent(
        f"""
        def qc_plots(adata):
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError("QC plots require at least one cell and one gene.")
            layer = {layer_argument}
            if layer is not None and layer not in adata.layers:
                raise ValueError(f"Expression layer not found: {{layer!r}}")
            if {use_raw!r} and adata.raw is None:
                raise ValueError("Raw expression is unavailable.")
            matrix = {matrix_code}
            source_var = adata.raw.var if {use_raw!r} else adata.var
            total_expression, detected_genes, mitochondrial_percent = _qc_plot_metric_values(
                matrix,
                source_var,
                np,
                sparse,
            )
            if mitochondrial_percent is not None and not np.isfinite(mitochondrial_percent).any():
                mitochondrial_percent = None
            group_labels = None
            group_positions = None
            if {mode!r} == "grouped":
                group_labels, group_positions, _ = _qc_plot_group_details(adata.obs, {groupby!r}, pd, np)
            return _qc_plot_figure(
                total_expression,
                detected_genes,
                mitochondrial_percent,
                {mode!r},
                {groupby!r},
                group_labels,
                group_positions,
                Figure,
                np,
            )
        """
    ).strip()
    imports = (
        "import numpy as np\nimport pandas as pd\nfrom matplotlib.figure import Figure\nfrom scipy import sparse\n"
        f"_QC_PLOT_MAX_GROUPS = {_QC_PLOT_MAX_GROUPS!r}"
    )
    return "\n\n".join((imports, helper_source, implementation))


def qc_plots_owned(
    adata: Any,
    *,
    view: Mapping[str, object] | None = None,
    source: DynamicExpressionSource | None = None,
) -> tuple[Any, Any, str]:
    parameters: dict[str, JSONValue] = {"source": source}
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    if cells == 0 or genes == 0:
        raise ValueError("QC plots require at least one cell and one gene.")

    resolved_view = _resolve_qc_plot_view(view)
    view_mode = resolved_view["view"]
    groupby = resolved_view.get("groupby")
    expression = QC_PLOTS_EXPRESSION_SOURCE.resolve(adata, _source_parameter(parameters))
    source_label = _expression_source_label(expression)
    matrix = expression.matrix(adata)
    source_var = _expression_var(adata, expression)
    total_counts, detected, pct_mt, metric_sources, warnings = _derive_qc_plot_metrics(
        matrix,
        source_var,
        source_label,
    )
    science = dependencies.require_scientific_dependencies()
    group_labels = None
    group_positions = None
    group_cell_counts = None
    if groupby is not None:
        group_labels, group_positions, group_cell_counts = _qc_plot_group_details(
            adata.obs,
            groupby,
            science.pd,
            science.np,
        )

    if pct_mt is not None:
        finite_pct_mt = science.np.isfinite(pct_mt)
        finite_pct_count = int(finite_pct_mt.sum())
        if finite_pct_count == 0:
            pct_mt = None
            metric_sources.pop("mitochondrial_percent", None)
            metric_sources.pop("mitochondrial_panel_total_expression", None)
            warnings.append(
                "Mitochondrial percentage annotations contained no finite values; the mitochondrial panel was omitted."
            )
        elif finite_pct_count != cells:
            warnings.append(
                f"{cells - finite_pct_count:,} cells with non-finite mitochondrial percentages were omitted "
                "from the mitochondrial panel and its distribution summary."
            )

    finite_total = total_counts[science.np.isfinite(total_counts)]
    finite_detected = detected[science.np.isfinite(detected)]
    if finite_total.size == 0 or finite_detected.size == 0:
        raise ValueError("QC plots require at least one finite total-expression and detected-gene value.")
    if finite_total.size != cells or finite_detected.size != cells:
        warnings.append("Non-finite QC observations were omitted from the affected plot panels.")
    figure = _qc_plot_figure(
        total_counts,
        detected,
        pct_mt,
        view_mode,
        groupby,
        group_labels,
        group_positions,
        science.Figure,
        science.np,
    )
    png = figure_to_png(figure)
    plot_parameters = {"view": resolved_view, **expression.parameters()}
    result = make_plot_result(
        title="Quality control plots",
        operation="qc_plots",
        parameters=plot_parameters,
        description=(
            "Cell count, feature, and mitochondrial quality metrics."
            if view_mode == "overview"
            else f"Cell-level quality metric distributions grouped by obs[{groupby!r}]."
        ),
        warnings=warnings,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        png=png,
    )
    distributions = {
        "total_expression": summarize_numeric(total_counts),
        "detected_genes": summarize_numeric(detected),
    }
    if pct_mt is not None:
        distributions["mitochondrial_percent"] = summarize_numeric(pct_mt)
    total_median = distributions["total_expression"]["median"]
    detected_median = distributions["detected_genes"]["median"]
    mito_median = distributions.get("mitochondrial_percent", {}).get("median")
    mito_text = (
        f" Median mitochondrial expression was {_display_number(mito_median)}%."
        if mito_median is not None
        else " Mitochondrial evidence was unavailable and its panel was omitted."
    )
    if view_mode == "overview":
        panels = [
            "total_expression_histogram",
            "detected_genes_histogram",
            "total_expression_vs_detected_genes",
            *(("total_expression_vs_mitochondrial_percent",) if pct_mt is not None else ()),
        ]
        grouped_results = {}
        methods = (
            "Cell-level total expression, detected-gene counts, and mitochondrial percentage were visualized as "
            f"distributions and pairwise QC scatter plots. Every metric was derived atomically from {source_label}; "
            "pre-existing observation-level QC columns were not used as an unverified cache."
        )
    else:
        panels = [
            "total_expression_by_group",
            "detected_genes_by_group",
            *(("mitochondrial_percent_by_group",) if pct_mt is not None else ()),
        ]
        group_distributions = {}
        for label, selected in zip(group_labels, group_positions, strict=True):
            group_distributions[label] = {
                "total_expression": summarize_numeric(total_counts[selected]),
                "detected_genes": summarize_numeric(detected[selected]),
            }
            if pct_mt is not None:
                group_distributions[label]["mitochondrial_percent"] = summarize_numeric(pct_mt[selected])
        grouped_results = {
            "view": view_mode,
            "groupby": groupby,
            "group_order": group_labels,
            "group_cell_counts": group_cell_counts,
            "group_distributions": group_distributions,
        }
        methods = (
            "Cell-level total expression, detected-gene counts, and available mitochondrial percentages were "
            f"visualized as grouped violin distributions over categorical obs[{groupby!r}]. Every metric was "
            f"derived atomically from {source_label}; pre-existing observation-level QC columns were not used as "
            "an unverified cache."
        )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellQCPlots",
        title="QC plot summary",
        operation="qc_plots_report",
        methods=methods,
        results=(
            f"The plot summarizes {cells:,} cells and {genes:,} genes. Median total expression was "
            f"{_display_number(total_median)} and the median number of detected genes was "
            f"{_display_number(detected_median)}.{mito_text}"
        ),
        key_results={
            "cells": cells,
            "genes": genes,
            "source_features": int(matrix.shape[1]),
            "metric_sources": metric_sources,
            "distributions": distributions,
            "panels": panels,
            **grouped_results,
        },
        parameters=plot_parameters,
        references=(SCANPY_REFERENCE, SCATER_REFERENCE, MITO_QC_REFERENCE),
        software_packages=QC_SOFTWARE_PACKAGES,
        warnings=warnings,
        limitations=QC_LIMITATIONS,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_qc_plots_code(expression, resolved_view),
    )
    return result, report, code


@register_operation("openbio.node.qcplots")
def qc_plots(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    operation = "QC Plots"
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(parameters, {"source", "view"}, operation=operation)
    result, report, code = qc_plots_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, result, kind=PLOT_KIND))


__all__ = [
    "calculate_qc",
    "calculate_qc_owned",
    "filter_cells",
    "filter_cells_owned",
    "filter_genes",
    "filter_genes_owned",
    "qc_plots",
    "qc_plots_owned",
]
