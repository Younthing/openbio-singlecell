from __future__ import annotations

import time
import warnings as python_warnings
from collections.abc import Mapping
from textwrap import dedent
from typing import Any

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, make_plot_result, matrix_totals_and_nonzero
from .expression_source import (
    _HVG_SPEC,
    _NORMALIZE_LAYER_SPEC,
    _NORMALIZE_TOTAL_SPEC,
    _PEARSON_RESIDUAL_SPEC,
    _SCALE_SPEC,
    ExpressionSource,
    ExpressionSourceSpec,
)
from .hvg_plotting import hvg_selection_plot_code, render_hvg_selection_plot
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
)
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

PREPROCESS_SOFTWARE_PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy")
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
NORMALIZE_TOTAL_SOURCE = _NORMALIZE_TOTAL_SPEC
NORMALIZE_LAYER_SOURCE = _NORMALIZE_LAYER_SPEC
PEARSON_SOURCE = _PEARSON_RESIDUAL_SPEC
HVG_SOURCE = _HVG_SPEC
SCALE_SOURCE = _SCALE_SPEC
SCANPY_REFERENCE = AnalysisReference(
    citation=(
        "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. "
        "Genome Biology. 2018;19:15."
    ),
    doi="10.1186/s13059-017-1382-0",
    url="https://doi.org/10.1186/s13059-017-1382-0",
    kind="software",
)
ANNDATA_REFERENCE = AnalysisReference(
    citation=(
        "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data "
        "matrices. Journal of Open Source Software. 2024;9(101):4371."
    ),
    doi="10.21105/joss.04371",
    url="https://doi.org/10.21105/joss.04371",
    kind="software",
)
ZHENG_REFERENCE = AnalysisReference(
    citation=(
        "Zheng GXY et al. Massively parallel digital transcriptional profiling of single cells. "
        "Nature Communications. 2017;8:14049."
    ),
    doi="10.1038/ncomms14049",
    url="https://doi.org/10.1038/ncomms14049",
    kind="method",
)
TRANSFORM_REFERENCE = AnalysisReference(
    citation=(
        "Ahlmann-Eltze C, Huber W. Comparison of transformations for single-cell RNA-seq data. "
        "Nature Methods. 2023;20:665-672."
    ),
    doi="10.1038/s41592-023-01814-1",
    url="https://doi.org/10.1038/s41592-023-01814-1",
    kind="method",
)
SATIJA_REFERENCE = AnalysisReference(
    citation=(
        "Satija R, Farrell JA, Gennert D, Schier AF, Regev A. Spatial reconstruction of single-cell "
        "gene expression data. Nature Biotechnology. 2015;33:495-502."
    ),
    doi="10.1038/nbt.3192",
    url="https://doi.org/10.1038/nbt.3192",
    kind="method",
)
STUART_REFERENCE = AnalysisReference(
    citation=("Stuart T et al. Comprehensive Integration of Single-Cell Data. Cell. 2019;177(7):1888-1902.e21."),
    doi="10.1016/j.cell.2019.05.031",
    url="https://doi.org/10.1016/j.cell.2019.05.031",
    kind="method",
)
LAUSE_REFERENCE = AnalysisReference(
    citation=(
        "Lause J, Berens P, Kobak D. Analytic Pearson residuals for normalization of single-cell "
        "RNA-seq UMI data. Genome Biology. 2021;22:258."
    ),
    doi="10.1186/s13059-021-02451-7",
    url="https://doi.org/10.1186/s13059-021-02451-7",
    kind="method",
)
HVG_REFERENCE_BY_FLAVOR = {
    "seurat": SATIJA_REFERENCE,
    "cell_ranger": ZHENG_REFERENCE,
    "seurat_v3": STUART_REFERENCE,
    "seurat_v3_paper": STUART_REFERENCE,
    "pearson_residuals": LAUSE_REFERENCE,
}

LOG_HVG_FLAVORS = frozenset({"seurat", "cell_ranger"})
COUNT_HVG_FLAVORS = frozenset({"seurat_v3", "seurat_v3_paper", "pearson_residuals"})
HVG_RESULT_COLUMNS = (
    "highly_variable",
    "highly_variable_algorithm",
    "highly_variable_forced",
    "highly_variable_rank",
    "highly_variable_nbatches",
    "highly_variable_intersection",
    "means",
    "dispersions",
    "dispersions_norm",
    "variances",
    "variances_norm",
    "residual_variances",
)


def _matrix_values(matrix: Any) -> Any:
    science = dependencies.require_scientific_dependencies()
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    return science.np.asarray(values)


def _matrix_description(matrix: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    total_entries = int(matrix.shape[0] * matrix.shape[1])
    if science.sparse.issparse(matrix):
        storage = f"sparse_{matrix.format}"
        nonzero_entries = int(matrix.count_nonzero())
        stored_values = science.np.asarray(matrix.data)
        value_min = float(stored_values.min()) if stored_values.size else 0.0
        value_max = float(stored_values.max()) if stored_values.size else 0.0
        if nonzero_entries < total_entries:
            value_min = min(value_min, 0.0)
            value_max = max(value_max, 0.0)
    else:
        storage = "dense"
        values = science.np.asarray(matrix)
        nonzero_entries = int(science.np.count_nonzero(values))
        value_min = float(values.min()) if values.size else None
        value_max = float(values.max()) if values.size else None
    return {
        "shape": [int(value) for value in matrix.shape],
        "dtype": str(matrix.dtype),
        "storage": storage,
        "total_entries": total_entries,
        "nonzero_entries": nonzero_entries,
        "value_min": value_min,
        "value_max": value_max,
    }


def _expression_source_label(expression: ExpressionSource) -> str:
    return f"layer:{expression.layer_name}" if expression.kind == "layer" else expression.kind


def _require_in_memory_nonempty(adata: Any, *, operation: str) -> None:
    if bool(adata.isbacked):
        raise ValueError(f"{operation} requires an in-memory AnnData.")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"{operation} requires at least one cell and one feature.")


def _require_unique_axes(adata: Any, *, operation: str) -> None:
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    if not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique feature identifiers.")


def _matrix_is_integer_like(matrix: Any) -> bool:
    science = dependencies.require_scientific_dependencies()
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix)
    for start in range(0, int(values.size), 1_000_000):
        block = science.np.asarray(values.flat[start : start + 1_000_000])
        if not bool(science.np.allclose(block, science.np.rint(block), rtol=0.0, atol=1e-8)):
            return False
    return True


def _dense_matrix_summary_memory_bounded(matrix: Any) -> dict[str, int | float | None]:
    """Summarize a dense result without allocating reporting arrays proportional to the result size."""

    science = dependencies.require_scientific_dependencies()
    values = science.np.asarray(matrix)
    finite_count = 0
    total = 0.0
    minimum: float | None = None
    maximum: float | None = None
    for start in range(0, int(values.size), 1_000_000):
        block = science.np.asarray(values.flat[start : start + 1_000_000], dtype=float)
        finite = block[science.np.isfinite(block)]
        if not finite.size:
            continue
        finite_count += int(finite.size)
        total += float(finite.sum(dtype=float))
        block_min = float(finite.min())
        block_max = float(finite.max())
        minimum = block_min if minimum is None else min(minimum, block_min)
        maximum = block_max if maximum is None else max(maximum, block_max)
    return {
        "n": int(values.size),
        "missing": int(values.size) - finite_count,
        "min": minimum,
        "q1": None,
        "median": None,
        "mean": None if finite_count == 0 else total / finite_count,
        "q3": None,
        "max": maximum,
    }


def _dense_matrix_all_finite(matrix: Any) -> bool:
    science = dependencies.require_scientific_dependencies()
    values = science.np.asarray(matrix)
    return all(
        bool(science.np.isfinite(science.np.asarray(values.flat[start : start + 1_000_000])).all())
        for start in range(0, int(values.size), 1_000_000)
    )


def _dense_matrix_boundary_count(matrix: Any, boundary: float) -> int:
    science = dependencies.require_scientific_dependencies()
    values = science.np.asarray(matrix)
    count = 0
    for start in range(0, int(values.size), 1_000_000):
        block = science.np.asarray(values.flat[start : start + 1_000_000])
        count += int(science.np.isclose(science.np.abs(block), boundary).sum())
    return count


def _validate_finite_numeric_matrix(matrix: Any, *, source_label: str) -> Any:
    science = dependencies.require_scientific_dependencies()
    if not hasattr(matrix, "shape") or len(matrix.shape) != 2:
        raise TypeError(f"{source_label} must be a two-dimensional numeric matrix.")
    if not (science.sparse.issparse(matrix) or isinstance(matrix, science.np.ndarray)):
        raise TypeError(f"{source_label} must be an in-memory NumPy or SciPy sparse matrix.")
    values = _matrix_values(matrix)
    if values.size and not bool(science.np.issubdtype(values.dtype, science.np.number)):
        raise TypeError(f"{source_label} must contain numeric values.")
    if values.size and not bool(science.np.isfinite(values).all()):
        raise ValueError(f"{source_label} contains non-finite expression values.")
    return values


def _validate_count_source(
    adata: Any,
    expression: ExpressionSource,
    *,
    operation: str,
) -> tuple[Any, Any, list[str]]:
    science = dependencies.require_scientific_dependencies()
    matrix = expression.matrix(adata)
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} source must align to the current AnnData shape.")
    values = _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
    warnings: list[str] = []
    if values.size and bool((values < 0).any()):
        warnings.append(
            f"{operation} source contains negative expression values. Scanpy can scale this signed matrix, "
            "but it is not a conventional non-negative count representation."
        )
    if not _matrix_is_integer_like(matrix):
        warnings.append(
            f"{operation} source contains non-integer values. Count totals are reported as expression sums, "
            "not verified UMI counts."
        )
    if values.size == 0 or not bool((values > 0).any()):
        warnings.append(f"{operation} source contains no positive expression values.")
    totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
    totals = science.np.asarray(totals, dtype=float)
    zero_total_cells = int((totals == 0).sum())
    if zero_total_cells:
        warnings.append(
            f"{operation} source contains {zero_total_cells} cells with zero total expression; Scanpy leaves "
            "those rows unchanged, so they cannot attain target_sum."
        )
    negative_total_cells = int((totals < 0).sum())
    if negative_total_cells:
        warnings.append(
            f"{operation} source contains {negative_total_cells} cells with negative total expression; the "
            "resulting scale-factor signs require expert interpretation."
        )
    return matrix, totals, warnings


def _validate_count_model_source(
    adata: Any,
    expression: ExpressionSource,
    *,
    operation: str,
    require_nonzero_totals: bool,
) -> tuple[Any, Any, Any, list[str]]:
    science = dependencies.require_scientific_dependencies()
    matrix = expression.matrix(adata)
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} source must align to the current AnnData shape.")
    values = _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
    warnings: list[str] = []
    if values.size and bool((values < 0).any()):
        warnings.append(
            f"{operation} source contains negative values; the selected expression is retained, but it is not "
            "a conventional non-negative count representation."
        )
    elif not _matrix_is_integer_like(matrix):
        warnings.append(
            f"{operation} source contains fractional non-negative values. The selected count-based method is "
            "executable, but the input cannot be described as verified integer UMI counts."
        )
    cell_totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
    gene_totals, _ = matrix_totals_and_nonzero(matrix, axis=0)
    cell_totals = science.np.asarray(cell_totals, dtype=float)
    gene_totals = science.np.asarray(gene_totals, dtype=float)
    zero_cells = int((cell_totals == 0).sum())
    if zero_cells:
        if require_nonzero_totals:
            raise ValueError(f"{operation} source contains {zero_cells} cells with zero total counts; filter them first.")
        warnings.append(f"{operation} source contains {zero_cells} cells with zero total expression; they were retained.")
    zero_genes = int((gene_totals == 0).sum())
    if require_nonzero_totals and zero_genes:
        raise ValueError(f"{operation} source contains {zero_genes} genes with zero total counts; filter them first.")
    return matrix, cell_totals, gene_totals, warnings


def _validate_logged_source(
    adata: Any,
    expression: ExpressionSource,
    *,
    operation: str,
) -> tuple[Any, Any, list[str]]:
    science = dependencies.require_scientific_dependencies()
    matrix = expression.matrix(adata)
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} source must align to the current AnnData shape.")
    values = _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
    warnings: list[str] = []
    if values.size and bool((values < 0).any()):
        warnings.append(
            f"{operation} source contains negative expression values; the selected expression is retained, "
            "but it is not a conventional normalized log count representation."
        )
    totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
    totals = science.np.asarray(totals, dtype=float)
    zero_cells = int((totals == 0).sum())
    if zero_cells:
        warnings.append(f"{operation} source contains {zero_cells} cells with zero total expression; they were retained.")
    return matrix, totals, warnings


def _number(value: JSONValue, *, name: str, positive: bool = False, nonnegative: bool = False) -> float:
    science = dependencies.require_scientific_dependencies()
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProtocolError(f"{name} must be a number.")
    number = float(value)
    if positive and (not bool(science.np.isfinite(number)) or number <= 0):
        raise ProtocolError(f"{name} must be finite and greater than zero.")
    if not positive and not bool(science.np.isfinite(number)):
        raise ProtocolError(f"{name} must be finite.")
    if nonnegative and number < 0:
        raise ProtocolError(f"{name} must be non-negative.")
    return number


def _positive_integer(value: JSONValue, *, name: str) -> int:
    if type(value) is not int or value < 1:
        raise ProtocolError(f"{name} must be an integer greater than zero.")
    return value


def _boolean(value: JSONValue, *, name: str) -> bool:
    if not isinstance(value, bool):
        raise ProtocolError(f"{name} must be a boolean.")
    return value


def _string(value: JSONValue, *, name: str, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise ProtocolError(f"{name} must be a string.")
    result = value.strip()
    if not allow_empty and not result:
        raise ProtocolError(f"{name} cannot be empty.")
    return result


def _source_parameter(parameters: dict[str, JSONValue], spec: ExpressionSourceSpec, adata: Any) -> ExpressionSource:
    value = parameters["source"]
    if not isinstance(value, dict):
        raise ProtocolError(f"{spec.description} must be a JSON object.")
    kind = value.get("source")
    expected = {"source", spec.layer_input_id} if kind == "layer" else {"source"}
    if set(value) != expected:
        raise ProtocolError(
            f"{spec.description} fields must be exactly {sorted(expected)!r}; got {sorted(value)!r}."
        )
    try:
        return spec.resolve(adata, value)
    except (TypeError, ValueError) as error:
        raise ProtocolError(str(error)) from error


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


def _validate_destination_layer(
    adata: Any,
    output_layer: JSONValue,
    overwrite_existing: JSONValue,
    *,
    operation: str,
) -> tuple[str, bool, bool]:
    output = _string(output_layer, name=f"{operation} output_layer")
    overwrite = _boolean(overwrite_existing, name=f"{operation} overwrite_existing")
    replaced = output in adata.layers
    if replaced and not overwrite:
        raise ValueError(f"{operation} output layer already exists: {output!r}")
    return output, overwrite, replaced


def _dense_array_gib(shape: tuple[int, int]) -> float:
    return float(int(shape[0]) * int(shape[1]) * 8 / 1024**3)


def _validate_dense_budget(shape: tuple[int, int], value: JSONValue, *, operation: str) -> tuple[float, float]:
    limit = _number(value, name=f"{operation} max_dense_gib", positive=True)
    estimate = _dense_array_gib(shape)
    if estimate > limit:
        raise ValueError(
            f"{operation} requires an estimated {estimate:.6g} GiB dense array, exceeding "
            f"max_dense_gib={limit:.6g}."
        )
    return estimate, limit


def _scale_feature_statistics(matrix: Any) -> tuple[Any, Any, Any]:
    science = dependencies.require_scientific_dependencies()
    if science.sparse.issparse(matrix):
        mean = science.np.asarray(matrix.mean(axis=0)).ravel()
        squared_mean = science.np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
        variance = (squared_mean - mean**2) * (matrix.shape[0] / (matrix.shape[0] - 1))
        variance = science.np.maximum(variance, 0.0)
    else:
        values = science.np.asarray(matrix, dtype=float)
        mean = values.mean(axis=0)
        variance = values.var(axis=0, ddof=1)
    std = science.np.sqrt(variance)
    return mean, std, science.np.isclose(std, 0.0)


def _log1p_code() -> str:
    return dedent(
        '''
        import warnings

        import numpy as np
        import scanpy as sc
        from scipy import sparse


        def log1p_expression(adata):
            operation = "Log1p"
            if adata.isbacked:
                raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError(f"{operation} requires at least one cell and one feature.")
            matrix = adata.X
            values = np.asarray(matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel())
            if values.size and not bool(np.isfinite(values).all()):
                raise ValueError("Log1p source contains non-finite expression values.")
            if values.size and bool((values <= -1).any()):
                raise ValueError("Log1p source contains values <= -1 outside the finite real log1p domain.")
            if values.size and bool((values < 0).any()):
                warnings.warn(
                    "Log1p source contains values in (-1, 0); the transform is finite but the matrix is "
                    "not conventional non-negative expression.",
                    UserWarning,
                    stacklevel=2,
                )
            sc.pp.log1p(adata, base=None, copy=False)
            output_values = np.asarray(
                adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).ravel()
            )
            if output_values.size and not np.isfinite(output_values).all():
                raise RuntimeError("Log1p produced non-finite expression values.")
            return adata
        '''
    ).strip()


@register_operation("openbio.node.log1p")
def log1p(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Log1p")
    require_parameters(parameters, set(), operation="Log1p")
    adata = read_anndata_input(inputs)
    if bool(adata.isbacked):
        raise ValueError("Log1p requires an in-memory AnnData.")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Log1p requires at least one cell and one feature.")
    matrix = adata.X
    values = _matrix_values(matrix)
    if values.size and not bool(dependencies.require_scientific_dependencies().np.isfinite(values).all()):
        raise ValueError("Log1p source contains non-finite expression values.")
    if values.size and bool((values <= -1).any()):
        raise ValueError("Log1p source contains values <= -1 outside the finite real log1p domain.")
    scientific_warnings = []
    if values.size and bool((values < 0).any()):
        scientific_warnings.append(
            "Log1p source contains values in (-1, 0); the transform is finite but the matrix is not "
            "conventional non-negative expression."
        )
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    input_description = _matrix_description(matrix)
    raw_present = adata.raw is not None
    layers_preserved = sorted(str(key) for key in adata.layers.keys() if key is not None)
    science = dependencies.require_scientific_dependencies()
    science.sc.pp.log1p(adata, base=None, copy=False)
    output_values = _matrix_values(adata.X)
    if output_values.size and not bool(science.np.isfinite(output_values).all()):
        raise RuntimeError("Log1p produced non-finite expression values.")
    report_parameters = {"source": "X", "base": "natural", "pseudocount": 1.0}
    finish_adata(
        adata,
        "log1p",
        report_parameters,
        cells,
        genes,
        started_at,
        warnings=scientific_warnings,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellLog1p",
        title="Log1p transformation summary",
        operation="log1p",
        methods=(
            f"Applied the natural-log transformation log(1+x) to active X for {cells:,} cells and "
            f"{genes:,} features with scanpy.pp.log1p on the worker's private AnnData."
        ),
        results=(
            f"Created a log1p expression representation in X for {cells:,} cells while preserving raw and "
            "all named layers unchanged."
        ),
        key_results={
            "cells": cells,
            "features": genes,
            "source": "X",
            "input_matrix": input_description,
            "output_matrix": _matrix_description(adata.X),
            "base": "natural",
            "pseudocount": 1.0,
            "raw_preserved": True,
            "raw_present": raw_present,
            "layers_preserved": layers_preserved,
        },
        parameters=report_parameters,
        references=(SCANPY_REFERENCE, TRANSFORM_REFERENCE, ANNDATA_REFERENCE),
        software_packages=PREPROCESS_SOFTWARE_PACKAGES,
        warnings=scientific_warnings,
        limitations=(
            "Log1p does not perform library-size normalization, Technical batch integration, or scaling.",
            "Count-based methods require an explicitly selected count expression source.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_log1p_code(),
    )
    return _records(context, adata, report, code)


def _generated_source_helpers() -> str:
    return dedent(
        """
        import warnings

        import numpy as np
        import scanpy as sc
        from anndata import AnnData
        from scipy import sparse


        def _openbio_validate_count_source(adata, matrix, operation):
            if tuple(matrix.shape) != tuple(adata.shape):
                raise ValueError(f"{operation} source must align to the current AnnData shape.")
            values = np.asarray(matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel())
            if values.size and not np.issubdtype(values.dtype, np.number):
                raise TypeError(f"{operation} source must contain numeric values.")
            if values.size and not bool(np.isfinite(values).all()):
                raise ValueError(f"{operation} source contains non-finite expression values.")
            scientific_warnings = []
            if values.size and bool((values < 0).any()):
                scientific_warnings.append(
                    f"{operation} source contains negative expression values. Scanpy can scale this signed "
                    "matrix, but it is not a conventional non-negative count representation."
                )
            if not bool(np.allclose(values, np.rint(values), rtol=0.0, atol=1e-8)):
                scientific_warnings.append(
                    f"{operation} source contains non-integer values. Count totals are reported as expression "
                    "sums, not verified UMI counts."
                )
            if values.size == 0 or not bool((values > 0).any()):
                scientific_warnings.append(f"{operation} source contains no positive expression values.")
            totals = np.asarray(matrix.sum(axis=1), dtype=float).ravel()
            zero_total_cells = int((totals == 0).sum())
            if zero_total_cells:
                scientific_warnings.append(
                    f"{operation} source contains {zero_total_cells} cells with zero total expression; "
                    "Scanpy leaves those rows unchanged, so they cannot attain target_sum."
                )
            negative_total_cells = int((totals < 0).sum())
            if negative_total_cells:
                scientific_warnings.append(
                    f"{operation} source contains {negative_total_cells} cells with negative total expression; "
                    "the resulting scale-factor signs require expert interpretation."
                )
            for message in scientific_warnings:
                warnings.warn(message, UserWarning, stacklevel=2)
            return totals
        """
    ).strip()


def _normalize_total_code(expression: ExpressionSource, target_sum: float) -> str:
    layer_name = expression.layer_name if expression.kind == "layer" else None
    execution = (
        "work = AnnData(X=matrix.copy()); "
        "sc.pp.normalize_total(work, target_sum=target_sum, exclude_highly_expressed=False, inplace=True); "
        "adata.X = work.X"
        if expression.kind == "layer"
        else "sc.pp.normalize_total(adata, target_sum=target_sum, exclude_highly_expressed=False, inplace=True)"
    )
    return (
        _generated_source_helpers()
        + "\n\n\n"
        + dedent(
            f"""
            def normalize_total_expression(adata):
                operation = "Normalize Total"
                if adata.isbacked:
                    raise ValueError(f"{{operation}} requires an in-memory AnnData; call adata.to_memory() first.")
                if adata.n_obs == 0 or adata.n_vars == 0:
                    raise ValueError(f"{{operation}} requires at least one cell and one feature.")
                target_sum = {target_sum!r}
                if not np.isfinite(target_sum) or target_sum <= 0:
                    raise ValueError(f"{{operation}} target_sum must be finite and greater than zero.")
                source_kind = {expression.kind!r}
                layer_name = {layer_name!r}
                if source_kind == "layer":
                    if layer_name not in adata.layers:
                        raise ValueError(f"{{operation}} source layer not found: {{layer_name!r}}")
                    matrix = adata.layers[layer_name]
                else:
                    matrix = adata.X
                _openbio_validate_count_source(adata, matrix, operation)
                {execution}
                normalized_values = np.asarray(
                    adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).ravel()
                )
                if normalized_values.size and not np.isfinite(normalized_values).all():
                    raise RuntimeError("Normalize Total produced non-finite expression values.")
                adata.uns.pop("log1p", None)
                return adata
            """
        ).strip()
    )


@register_operation("openbio.node.normalizetotal")
def normalize_total(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    adata = _input_anndata(
        inputs,
        parameters,
        operation="Normalize Total",
        expected_parameters={"target_sum", "source"},
    )
    _require_in_memory_nonempty(adata, operation="Normalize Total")
    target_sum = _number(parameters["target_sum"], name="Normalize Total target_sum", positive=True)
    expression = _source_parameter(parameters, NORMALIZE_TOTAL_SOURCE, adata)
    matrix, input_totals, report_warnings = _validate_count_source(
        adata,
        expression,
        operation="Normalize Total",
    )
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    input_description = _matrix_description(matrix)
    input_integer_like = _matrix_is_integer_like(matrix)
    raw_present = adata.raw is not None
    layers_preserved = sorted(str(key) for key in adata.layers.keys())
    stale_log1p_marker_removed = isinstance(adata.uns.get("log1p"), Mapping)

    if expression.kind == "X":
        # The worker owns this AnnData, so Scanpy can normalize X directly without a defensive AnnData copy.
        with python_warnings.catch_warnings(record=True) as caught:
            python_warnings.simplefilter("always")
            science.sc.pp.normalize_total(
                adata,
                target_sum=target_sum,
                exclude_highly_expressed=False,
                inplace=True,
            )
    else:
        # A matrix-only workspace preserves the selected layer while the normalized representation replaces X.
        work = science.ad.AnnData(X=matrix.copy())
        with python_warnings.catch_warnings(record=True) as caught:
            python_warnings.simplefilter("always")
            science.sc.pp.normalize_total(
                work,
                target_sum=target_sum,
                exclude_highly_expressed=False,
                inplace=True,
            )
        adata.X = work.X
    report_warnings.extend(f"Scanpy normalize_total warning: {item.message}" for item in caught)
    output_values = _matrix_values(adata.X)
    if output_values.size and not bool(science.np.isfinite(output_values).all()):
        raise RuntimeError("Normalize Total produced non-finite expression values.")
    adata.uns.pop("log1p", None)
    output_totals, _ = matrix_totals_and_nonzero(adata.X, axis=1)
    report_parameters = {
        "target_sum": target_sum,
        "exclude_highly_expressed": False,
        **expression.parameters(),
    }
    finish_adata(
        adata,
        "normalize_total",
        report_parameters,
        cells,
        genes,
        started_at,
        warnings=report_warnings,
    )
    source_label = _expression_source_label(expression)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellNormalizeTotal",
        title="Total-count normalization summary",
        operation="normalize_total",
        methods=(
            f"Each of {cells:,} cells was scaled from expert-selected source {source_label!r} with "
            f"scanpy.pp.normalize_total to a target total of {target_sum:,.6g}. The transformed matrix was "
            "written to X on the worker's private AnnData; raw and every named layer were preserved."
        ),
        results=(
            f"Total-count-normalized expression was produced for {cells:,} cells and {genes:,} features. "
            f"The median output cell total was {float(science.np.median(output_totals)):,.6g}."
        ),
        key_results={
            "cells": cells,
            "features": genes,
            "source": source_label,
            "input_matrix": input_description,
            "input_integer_like": input_integer_like,
            "input_cell_totals": summarize_numeric(input_totals),
            "zero_total_cells": int((science.np.asarray(input_totals) == 0).sum()),
            "negative_total_cells": int((science.np.asarray(input_totals) < 0).sum()),
            "target_sum": target_sum,
            "output_matrix": _matrix_description(adata.X),
            "output_cell_totals": summarize_numeric(output_totals),
            "stale_log1p_marker_removed": stale_log1p_marker_removed,
            "raw_preserved": True,
            "raw_present": raw_present,
            "layers_preserved": layers_preserved,
        },
        parameters=report_parameters,
        references=(SCANPY_REFERENCE, ZHENG_REFERENCE),
        software_packages=PREPROCESS_SOFTWARE_PACKAGES,
        warnings=report_warnings,
        limitations=(
            "Total-count scaling adjusts library depth only; it is not Technical batch integration or "
            "Sample-level inference.",
            "Count-based methods require an explicitly selected count expression source.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_normalize_total_code(expression, target_sum),
    )
    return _records(context, adata, report, code)


def _normalize_to_layer_code(
    expression: ExpressionSource,
    target_sum: float,
    transform: str,
    output_layer: str,
    overwrite_existing: bool,
) -> str:
    layer_name = expression.layer_name if expression.kind == "layer" else None
    return (
        _generated_source_helpers()
        + "\n\n\n"
        + dedent(
            f"""
            def normalize_expression_to_layer(adata):
                operation = "Normalize to Layer"
                if adata.isbacked:
                    raise ValueError(f"{{operation}} requires an in-memory AnnData; call adata.to_memory() first.")
                if adata.n_obs == 0 or adata.n_vars == 0:
                    raise ValueError(f"{{operation}} requires at least one cell and one feature.")
                source_kind = {expression.kind!r}
                source_layer = {layer_name!r}
                target_sum = {target_sum!r}
                transform = {transform!r}
                output_layer = {output_layer!r}
                overwrite_existing = {overwrite_existing!r}
                if not np.isfinite(target_sum) or target_sum <= 0:
                    raise ValueError(f"{{operation}} target_sum must be finite and greater than zero.")
                if transform not in {{"none", "log1p", "sqrt"}}:
                    raise ValueError(f"Unsupported Normalize to Layer transform: {{transform!r}}")
                if not output_layer:
                    raise ValueError("Normalize to Layer output_layer cannot be empty.")
                if output_layer in adata.layers and not overwrite_existing:
                    raise ValueError(f"Normalize to Layer output layer already exists: {{output_layer!r}}")
                if source_kind == "layer":
                    if source_layer not in adata.layers:
                        raise ValueError(f"{{operation}} source layer not found: {{source_layer!r}}")
                    matrix = adata.layers[source_layer]
                else:
                    matrix = adata.X
                _openbio_validate_count_source(adata, matrix, operation)
                work = AnnData(X=matrix.copy())
                sc.pp.normalize_total(
                    work,
                    target_sum=target_sum,
                    exclude_highly_expressed=False,
                    inplace=True,
                )
                normalized_values = np.asarray(
                    work.X.data if sparse.issparse(work.X) else np.asarray(work.X).ravel()
                )
                if normalized_values.size and not np.isfinite(normalized_values).all():
                    raise RuntimeError("Normalize to Layer produced non-finite normalized values.")
                if transform == "log1p":
                    if normalized_values.size and (normalized_values <= -1).any():
                        raise ValueError(
                            "Normalize to Layer log1p transform requires every normalized value to be greater than -1."
                        )
                    sc.pp.log1p(work, base=None)
                elif transform == "sqrt":
                    if normalized_values.size and (normalized_values < 0).any():
                        raise ValueError(
                            "Normalize to Layer sqrt transform requires non-negative normalized values."
                        )
                    work.X = work.X.sqrt() if sparse.issparse(work.X) else np.sqrt(work.X)
                output_values = np.asarray(
                    work.X.data if sparse.issparse(work.X) else np.asarray(work.X).ravel()
                )
                if output_values.size and not np.isfinite(output_values).all():
                    raise RuntimeError("Normalize to Layer produced non-finite output values.")
                adata.layers[output_layer] = work.X
                return adata
            """
        ).strip()
    )


@register_operation("openbio.node.normalizetolayer")
def normalize_to_layer(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    adata = _input_anndata(
        inputs,
        parameters,
        operation="Normalize to Layer",
        expected_parameters={"source", "target_sum", "transform", "output_layer", "overwrite_existing"},
    )
    operation = "Normalize to Layer"
    _require_in_memory_nonempty(adata, operation=operation)
    expression = _source_parameter(parameters, NORMALIZE_LAYER_SOURCE, adata)
    target_sum = _number(parameters["target_sum"], name=f"{operation} target_sum", positive=True)
    transform = _string(parameters["transform"], name=f"{operation} transform")
    if transform not in {"none", "log1p", "sqrt"}:
        raise ProtocolError(f"Unsupported {operation} transform: {transform!r}")
    output_layer, overwrite_existing, replaced_existing = _validate_destination_layer(
        adata,
        parameters["output_layer"],
        parameters["overwrite_existing"],
        operation=operation,
    )
    matrix, input_totals, report_warnings = _validate_count_source(
        adata,
        expression,
        operation=operation,
    )
    science = dependencies.require_scientific_dependencies()
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    input_description = _matrix_description(matrix)
    raw_present = adata.raw is not None
    # Transform a separate matrix before assigning the named output, preserving X and every other layer.
    work = science.ad.AnnData(X=matrix.copy())
    with python_warnings.catch_warnings(record=True) as caught:
        python_warnings.simplefilter("always")
        science.sc.pp.normalize_total(
            work,
            target_sum=target_sum,
            exclude_highly_expressed=False,
            inplace=True,
        )
    report_warnings.extend(f"Scanpy normalize_total warning: {item.message}" for item in caught)
    normalized_totals, _ = matrix_totals_and_nonzero(work.X, axis=1)
    normalized_values = _matrix_values(work.X)
    if normalized_values.size and not bool(science.np.isfinite(normalized_values).all()):
        raise RuntimeError("Normalize to Layer produced non-finite normalized values.")
    if transform == "log1p":
        if normalized_values.size and bool((normalized_values <= -1).any()):
            raise ValueError("Normalize to Layer log1p transform requires every normalized value to be greater than -1.")
        science.sc.pp.log1p(work, base=None, copy=False)
    elif transform == "sqrt":
        if normalized_values.size and bool((normalized_values < 0).any()):
            raise ValueError("Normalize to Layer sqrt transform requires non-negative normalized values.")
        work.X = work.X.sqrt() if science.sparse.issparse(work.X) else science.np.sqrt(work.X)
    output_values = _matrix_values(work.X)
    if output_values.size and not bool(science.np.isfinite(output_values).all()):
        raise RuntimeError("Normalize to Layer produced non-finite output values.")
    adata.layers[output_layer] = work.X
    report_parameters = {
        **expression.parameters(),
        "target_sum": target_sum,
        "exclude_highly_expressed": False,
        "transform": transform,
        "output_layer": output_layer,
        "overwrite_existing": overwrite_existing,
    }
    finish_adata(
        adata,
        "normalize_to_layer",
        report_parameters,
        cells,
        genes,
        started_at,
        warnings=report_warnings,
    )
    source_label = _expression_source_label(expression)
    transform_description = {
        "none": "no post-normalization transform",
        "log1p": "natural log(1+x)",
        "sqrt": "element-wise square root",
    }[transform]
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellNormalizeToLayer",
        title="Derived normalized layer summary",
        operation="normalize_to_layer",
        methods=(
            f"Scaled each cell from expert-selected source {source_label!r} to total {target_sum:,.6g} with "
            f"scanpy.pp.normalize_total, applied {transform_description}, and stored the result in AnnData "
            f"layer {output_layer!r} without replacing X or the Raw snapshot."
        ),
        results=(
            f"Created derived layer {output_layer!r} for {cells:,} cells and {genes:,} features using "
            f"transform {transform!r}."
        ),
        key_results={
            "cells": cells,
            "features": genes,
            "source": source_label,
            "input_matrix": input_description,
            "input_integer_like": _matrix_is_integer_like(matrix),
            "input_cell_totals": summarize_numeric(input_totals),
            "zero_total_cells": int((science.np.asarray(input_totals) == 0).sum()),
            "negative_total_cells": int((science.np.asarray(input_totals) < 0).sum()),
            "target_sum": target_sum,
            "intermediate_normalized_cell_totals": summarize_numeric(normalized_totals),
            "transform": transform,
            "output_layer": output_layer,
            "output_matrix": _matrix_description(adata.layers[output_layer]),
            "replaced_existing_layer": replaced_existing,
            "x_preserved": True,
            "source_preserved": expression.kind != "layer" or output_layer != expression.layer_name,
            "raw_preserved": True,
            "raw_present": raw_present,
            "canonical_counts_preserved": "counts" in adata.layers and output_layer != "counts",
        },
        parameters=report_parameters,
        references=(SCANPY_REFERENCE, TRANSFORM_REFERENCE, ANNDATA_REFERENCE),
        software_packages=PREPROCESS_SOFTWARE_PACKAGES,
        warnings=report_warnings,
        limitations=(
            "Only the intermediate linear normalized matrix has cell totals equal to target_sum; transformed "
            "layer sums do not.",
            "Count-based methods require an explicitly selected count expression source.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_normalize_to_layer_code(expression, target_sum, transform, output_layer, overwrite_existing),
    )
    return _records(context, adata, report, code)


def _generated_feature_helpers() -> str:
    return (
        _generated_source_helpers()
        + "\n\n\n"
        + dedent(
            """
            _OPENBIO_HVG_COLUMNS = (
                "highly_variable",
                "highly_variable_algorithm",
                "highly_variable_forced",
                "highly_variable_rank",
                "highly_variable_nbatches",
                "highly_variable_intersection",
                "means",
                "dispersions",
                "dispersions_norm",
                "variances",
                "variances_norm",
                "residual_variances",
            )


            def _openbio_require_feature_input(adata, operation, minimum_cells=2):
                if adata.isbacked:
                    raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
                if adata.n_obs < minimum_cells or adata.n_vars < 1:
                    raise ValueError(
                        f"{operation} requires at least {minimum_cells} cell(s) and one feature."
                    )
                if not adata.obs_names.is_unique:
                    raise ValueError(f"{operation} requires unique observation identifiers.")
                if not adata.var_names.is_unique:
                    raise ValueError(f"{operation} requires unique feature identifiers.")


            def _openbio_feature_source(adata, source_kind, source_layer, operation):
                if source_kind == "layer":
                    if source_layer not in adata.layers:
                        raise ValueError(f"{operation} source layer not found: {source_layer!r}")
                    matrix = adata.layers[source_layer]
                elif source_kind == "X":
                    matrix = adata.X
                else:
                    raise ValueError(f"Unsupported {operation} source: {source_kind!r}")
                if tuple(matrix.shape) != tuple(adata.shape):
                    raise ValueError(f"{operation} source must align to the current AnnData shape.")
                if not (sparse.issparse(matrix) or isinstance(matrix, np.ndarray)):
                    raise TypeError(f"{operation} source must be an in-memory NumPy or SciPy sparse matrix.")
                values = np.asarray(matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel())
                if values.size and not np.issubdtype(values.dtype, np.number):
                    raise TypeError(f"{operation} source must contain numeric values.")
                if values.size and not np.isfinite(values).all():
                    raise ValueError(f"{operation} source contains non-finite expression values.")
                return matrix, values


            def _openbio_validate_count_model(matrix, values, operation, require_nonzero_totals):
                if values.size and (values < 0).any():
                    warnings.warn(
                        f"{operation} source contains negative values; the selected expression is retained, but it is "
                        "not a conventional non-negative count representation.",
                        UserWarning,
                        stacklevel=2,
                    )
                elif not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-8):
                    warnings.warn(
                        f"{operation} source contains fractional non-negative values. The selected count-based "
                        "method is executable, but the input cannot be described as verified integer UMI counts.",
                        UserWarning,
                        stacklevel=2,
                    )
                cell_totals = np.asarray(matrix.sum(axis=1), dtype=float).ravel()
                gene_totals = np.asarray(matrix.sum(axis=0), dtype=float).ravel()
                zero_cells = int((cell_totals == 0).sum())
                if zero_cells:
                    if require_nonzero_totals:
                        raise ValueError(
                            f"{operation} source contains {zero_cells} cells with zero total counts; filter them first."
                        )
                    warnings.warn(
                        f"{operation} source contains {zero_cells} cells with zero total expression; they were retained.",
                        UserWarning,
                        stacklevel=2,
                    )
                zero_genes = int((gene_totals == 0).sum())
                if require_nonzero_totals and zero_genes:
                    raise ValueError(
                        f"{operation} source contains {zero_genes} genes with zero total counts; filter them first."
                    )


            def _openbio_validate_logged(matrix, values, operation):
                if values.size and (values < 0).any():
                    warnings.warn(
                        f"{operation} source contains negative expression values; the selected expression is retained, "
                        "but it is not a conventional normalized log count representation.",
                        UserWarning,
                        stacklevel=2,
                    )
                cell_totals = np.asarray(matrix.sum(axis=1), dtype=float).ravel()
                zero_cells = int((cell_totals == 0).sum())
                if zero_cells:
                    warnings.warn(
                        f"{operation} source contains {zero_cells} cells with zero total expression; they were retained.",
                        UserWarning,
                        stacklevel=2,
                    )


            def _openbio_validate_destination(adata, output_layer, overwrite_existing, operation):
                if not isinstance(output_layer, str) or not output_layer.strip():
                    raise ValueError(f"{operation} output_layer cannot be empty.")
                output_layer = output_layer.strip()
                if output_layer in adata.layers and not overwrite_existing:
                    raise ValueError(f"{operation} output layer already exists: {output_layer!r}")
                return output_layer


            def _openbio_validate_dense_budget(shape, max_dense_gib, operation):
                if not np.isfinite(max_dense_gib) or max_dense_gib <= 0:
                    raise ValueError(f"{operation} max_dense_gib must be finite and greater than zero.")
                estimate = int(shape[0]) * int(shape[1]) * 8 / 1024**3
                if estimate > max_dense_gib:
                    raise ValueError(
                        f"{operation} requires an estimated {estimate:.6g} GiB dense array, "
                        f"exceeding max_dense_gib={max_dense_gib:.6g}."
                    )


            def _openbio_validate_batch(adata, batch_key, operation, minimum_group_size):
                if not batch_key:
                    return
                if batch_key not in adata.obs:
                    raise ValueError(f"{operation} batch column not found in obs: {batch_key!r}")
                values = adata.obs[batch_key]
                if values.isna().any():
                    raise ValueError(f"{operation} batch column contains missing labels: {batch_key!r}")
                if any(isinstance(value, str) and not value.strip() for value in values.astype(object)):
                    raise ValueError(f"{operation} batch column contains blank labels: {batch_key!r}")
                counts = values.value_counts(sort=False, dropna=False)
                small = [(value, int(count)) for value, count in counts.items() if int(count) < minimum_group_size]
                if small:
                    raise ValueError(
                        f"{operation} backend requires at least {minimum_group_size} cells in every batch group; "
                        f"found {small!r}."
                    )
                singleton = [(value, int(count)) for value, count in counts.items() if int(count) == 1]
                if singleton:
                    warnings.warn(
                        f"Batch-aware selection includes singleton group(s) {singleton!r}; the selected flavor "
                        "can execute, but within-group variation is not independently supported.",
                        UserWarning,
                        stacklevel=2,
                    )
            """
        ).strip()
    )


def _python_float_literal(value: float | None) -> str:
    if value is None:
        return "None"
    if value == float("inf"):
        return 'float("inf")'
    return repr(value)


def _pearson_residuals_code(
    expression: ExpressionSource,
    *,
    theta: float,
    clip: float | None,
    output_layer: str,
    overwrite_existing: bool,
    max_dense_gib: float,
) -> str:
    source_layer = expression.layer_name if expression.kind == "layer" else None
    return (
        _generated_feature_helpers()
        + "\n\n\n"
        + dedent(
            f"""
            def pearson_residuals_to_layer(adata):
                operation = "Pearson Residuals"
                _openbio_require_feature_input(adata, operation, minimum_cells=1)
                source_kind = {expression.kind!r}
                source_layer = {source_layer!r}
                theta = {_python_float_literal(theta)}
                clip = {_python_float_literal(clip)}
                output_layer = {output_layer!r}
                overwrite_existing = {overwrite_existing!r}
                max_dense_gib = {max_dense_gib!r}
                if not (theta > 0 and (np.isfinite(theta) or theta == float("inf"))):
                    raise ValueError("Pearson Residuals theta must be greater than zero.")
                matrix, values = _openbio_feature_source(adata, source_kind, source_layer, operation)
                _openbio_validate_count_model(
                    matrix,
                    values,
                    operation,
                    require_nonzero_totals=True,
                )
                output_layer = _openbio_validate_destination(
                    adata,
                    output_layer,
                    overwrite_existing,
                    operation,
                )
                _openbio_validate_dense_budget(adata.shape, max_dense_gib, operation)
                work = AnnData(X=matrix.copy())
                normalized = sc.experimental.pp.normalize_pearson_residuals(
                    work,
                    theta=theta,
                    clip=clip,
                    check_values=False,
                    layer=None,
                    inplace=False,
                )
                residuals = np.asarray(normalized["X"])
                if not np.isfinite(residuals).all():
                    raise ValueError("Pearson Residuals produced non-finite values.")
                adata.layers[output_layer] = residuals
                return adata
            """
        ).strip()
    )


@register_operation("openbio.node.pearsonresidualstolayer")
def pearson_residuals_to_layer(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    expected = {
        "source",
        "theta",
        "clipping_mode",
        "custom_clip",
        "output_layer",
        "overwrite_existing",
        "max_dense_gib",
    }
    adata = _input_anndata(inputs, parameters, operation="Pearson Residuals", expected_parameters=expected)
    operation = "Pearson Residuals"
    _require_in_memory_nonempty(adata, operation=operation)
    _require_unique_axes(adata, operation=operation)
    expression = _source_parameter(parameters, PEARSON_SOURCE, adata)
    matrix, cell_totals, gene_totals, report_warnings = _validate_count_model_source(
        adata,
        expression,
        operation=operation,
        require_nonzero_totals=True,
    )
    theta = _number(parameters["theta"], name=f"{operation} theta", positive=True)
    clipping_mode = _string(parameters["clipping_mode"], name=f"{operation} clipping_mode")
    if clipping_mode not in {"sqrt_n_obs", "custom", "none"}:
        raise ProtocolError(f"Unsupported {operation} clipping_mode: {clipping_mode!r}")
    if clipping_mode == "sqrt_n_obs":
        scanpy_clip: float | None = float(dependencies.require_scientific_dependencies().np.sqrt(adata.n_obs))
        resolved_clip: float | str = scanpy_clip
    elif clipping_mode == "custom":
        scanpy_clip = _number(parameters["custom_clip"], name=f"{operation} custom_clip", nonnegative=True)
        resolved_clip = scanpy_clip
    else:
        scanpy_clip = None
        resolved_clip = "infinity"
    output_layer, overwrite_existing, replaced_existing = _validate_destination_layer(
        adata,
        parameters["output_layer"],
        parameters["overwrite_existing"],
        operation=operation,
    )
    estimated_dense_gib, max_dense_gib = _validate_dense_budget(
        tuple(adata.shape),
        parameters["max_dense_gib"],
        operation=operation,
    )
    science = dependencies.require_scientific_dependencies()
    input_sparse = bool(science.sparse.issparse(matrix))
    if input_sparse:
        report_warnings.append(
            "Analytic Pearson residuals convert sparse counts to a dense float64 matrix; the reported "
            "allocation is a lower bound because working arrays coexist."
        )
    if adata.n_obs == 1:
        report_warnings.append(
            "Pearson Residuals received one cell. Scanpy defines the residual transform, but per-gene sample "
            "variance with ddof=1 is unavailable and is reported as null."
        )
    if adata.n_vars == 1:
        report_warnings.append(
            "Pearson Residuals received one feature; the transform is finite but offers no multigene "
            "representation for downstream dimensionality reduction."
        )
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    input_description = _matrix_description(matrix)
    raw_present = adata.raw is not None
    # Scanpy's API consumes an AnnData; isolate the transform until the named output is assigned.
    work = science.ad.AnnData(X=matrix.copy())
    normalized = science.sc.experimental.pp.normalize_pearson_residuals(
        work,
        theta=theta,
        clip=scanpy_clip,
        check_values=False,
        layer=None,
        inplace=False,
    )
    residuals = science.np.asarray(normalized["X"])
    if not _dense_matrix_all_finite(residuals):
        raise ValueError(
            "Pearson Residuals produced non-finite values; verify that zero-total cells and genes were removed."
        )
    adata.layers[output_layer] = residuals
    residual_variance = (
        residuals.var(axis=0, ddof=1) if cells > 1 else science.np.full(genes, science.np.nan, dtype=float)
    )
    at_clip = 0 if scanpy_clip is None else _dense_matrix_boundary_count(residuals, scanpy_clip)
    residual_summary = _dense_matrix_summary_memory_bounded(residuals)
    report_parameters = {
        **expression.parameters(),
        "theta": theta,
        "clipping_mode": clipping_mode,
        "resolved_clip": resolved_clip,
        "check_values": False,
        "output_layer": output_layer,
        "overwrite_existing": overwrite_existing,
        "max_dense_gib": max_dense_gib,
    }
    finish_adata(
        adata,
        "pearson_residuals_to_layer",
        report_parameters,
        cells,
        genes,
        started_at,
        warnings=report_warnings,
    )
    source_label = _expression_source_label(expression)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellPearsonResidualsToLayer",
        title="Analytic Pearson residual transformation summary",
        operation="pearson_residuals_to_layer",
        methods=(
            f"Computed analytic Pearson residuals from expert-selected expression source {source_label!r} with "
            f"scanpy.experimental.pp.normalize_pearson_residuals using theta {theta!r} and clipping mode "
            f"{clipping_mode!r} (resolved bound {resolved_clip!r}). The dense residual matrix was stored in "
            f"layer {output_layer!r}; X, every other layer, and Raw were preserved."
        ),
        results=(
            f"Created residual representation {output_layer!r} for {cells:,} cells and {genes:,} genes. "
            f"Residual values ranged from {float(residuals.min()):.6g} to {float(residuals.max()):.6g}; "
            f"{at_clip:,} entries were at the resolved clipping boundary."
        ),
        key_results={
            "cells": cells,
            "features": genes,
            "source": source_label,
            "input_matrix": input_description,
            "input_integer_like": _matrix_is_integer_like(matrix),
            "input_cell_totals": summarize_numeric(cell_totals),
            "input_gene_totals": summarize_numeric(gene_totals),
            "zero_total_cells": 0,
            "zero_total_genes": 0,
            "theta": theta,
            "clipping_mode": clipping_mode,
            "resolved_clip": resolved_clip,
            "output_layer": output_layer,
            "output_matrix": _matrix_description(residuals),
            "residual_values": residual_summary,
            "residual_variances_ddof1": summarize_numeric(residual_variance),
            "entries_at_clipping_bound": at_clip,
            "estimated_dense_array_gib": estimated_dense_gib,
            "dense_budget_gib": max_dense_gib,
            "sparse_to_dense": input_sparse,
            "replaced_existing_layer": replaced_existing,
            "x_preserved": True,
            "source_preserved": expression.kind != "layer" or output_layer != expression.layer_name,
            "raw_preserved": True,
            "raw_present": raw_present,
        },
        parameters=report_parameters,
        references=(LAUSE_REFERENCE, SCANPY_REFERENCE, ANNDATA_REFERENCE),
        software_packages=PREPROCESS_SOFTWARE_PACKAGES,
        warnings=report_warnings,
        limitations=(
            "Analytic Pearson residuals are a model-based expression representation, not counts or "
            "sample-level inference.",
            "The Scanpy API is experimental; the dense allocation estimate covers one float64 result array "
            "and not all simultaneous intermediates.",
            "Residual value quartiles are omitted from the summary to avoid reporting-only full-matrix "
            "allocations; exact minimum, mean, and maximum are reported blockwise.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_pearson_residuals_code(
            expression,
            theta=theta,
            clip=scanpy_clip,
            output_layer=output_layer,
            overwrite_existing=overwrite_existing,
            max_dense_gib=max_dense_gib,
        ),
    )
    return _records(context, adata, report, code)


def _batch_group_sizes(
    adata: Any,
    batch_key: str,
    *,
    operation: str,
    minimum_group_size: int,
) -> tuple[list[dict[str, Any]], list[str]]:
    if not batch_key:
        return [], []
    if batch_key not in adata.obs:
        raise ValueError(f"{operation} batch column not found in obs: {batch_key!r}")
    values = adata.obs[batch_key]
    if bool(values.isna().any()):
        raise ValueError(f"{operation} batch column contains missing labels: {batch_key!r}")
    if any(isinstance(value, str) and not value.strip() for value in values.astype(object)):
        raise ValueError(f"{operation} batch column contains blank labels: {batch_key!r}")
    counts = values.value_counts(sort=False, dropna=False)
    groups = [{"label": str(value), "count": int(count)} for value, count in counts.items() if int(count) > 0]
    small = [group for group in groups if group["count"] < minimum_group_size]
    if small:
        raise ValueError(
            f"{operation} backend requires at least {minimum_group_size} cells in every batch group; "
            f"found {small!r}."
        )
    warnings: list[str] = []
    if len(groups) == 1:
        warnings.append(f"Batch-aware selection used {batch_key!r}, but the column has only one observed level.")
    singleton = [group for group in groups if group["count"] == 1]
    if singleton:
        warnings.append(
            f"Batch-aware selection includes singleton group(s) {singleton!r}; the selected flavor can "
            "execute, but within-group variation is not independently supported."
        )
    return groups, warnings


def _highly_variable_genes_code(
    expression: ExpressionSource,
    *,
    n_top_genes: int,
    flavor: str,
    batch_key: str,
    requested_keep: list[str],
    subset: bool,
    theta: float | None,
    clip: float | None,
    chunksize: int | None,
    span: float | None,
    n_bins: int | None,
    overwrite_existing: bool,
) -> str:
    source_layer = expression.layer_name if expression.kind == "layer" else None
    return (
        _generated_feature_helpers()
        + "\n\n\n"
        + dedent(
            f"""
            def select_highly_variable_genes(adata):
                operation = "Highly Variable Genes"
                source_kind = {expression.kind!r}
                source_layer = {source_layer!r}
                n_top_genes = {n_top_genes!r}
                flavor = {flavor!r}
                batch_key = {batch_key!r}
                requested_keep = {requested_keep!r}
                subset = {subset!r}
                theta = {_python_float_literal(theta)}
                clip = {_python_float_literal(clip)}
                chunksize = {chunksize!r}
                span = {span!r}
                n_bins = {n_bins!r}
                overwrite_existing = {overwrite_existing!r}
                minimum_cells = 1 if flavor == "pearson_residuals" else 2
                _openbio_require_feature_input(adata, operation, minimum_cells=minimum_cells)
                if n_top_genes < 1:
                    raise ValueError("Highly Variable Genes n_top_genes must be at least one.")
                if n_top_genes > adata.n_vars:
                    warnings.warn(
                        f"Highly Variable Genes requested {{n_top_genes}} genes from {{adata.n_vars}} features; "
                        "Scanpy will return at most the eligible feature count.",
                        UserWarning,
                        stacklevel=2,
                    )
                matrix, values = _openbio_feature_source(adata, source_kind, source_layer, operation)
                if flavor in {{"seurat", "cell_ranger"}}:
                    _openbio_validate_logged(matrix, values, operation)
                else:
                    _openbio_validate_count_model(
                        matrix,
                        values,
                        operation,
                        require_nonzero_totals=False,
                    )
                expressed_features = int((np.asarray(matrix.sum(axis=0), dtype=float).ravel() > 0).sum())
                if expressed_features < n_top_genes:
                    warnings.warn(
                        f"{{operation}} requested {{n_top_genes}} genes but only {{expressed_features}} features "
                        "have positive totals; the algorithm may return fewer.",
                        UserWarning,
                        stacklevel=2,
                    )
                _openbio_validate_batch(
                    adata,
                    batch_key,
                    operation,
                    1 if flavor == "pearson_residuals" else 2,
                )
                stale = [column for column in _OPENBIO_HVG_COLUMNS if column in adata.var]
                if (stale or "hvg" in adata.uns) and not overwrite_existing:
                    raise ValueError(
                        "Highly Variable Genes annotations already exist; enable overwrite_existing to replace them."
                    )
                adata.var.drop(columns=stale, inplace=True)
                adata.uns.pop("hvg", None)
                layer = source_layer if source_kind == "layer" else None
                if flavor == "pearson_residuals":
                    sc.experimental.pp.highly_variable_genes(
                        adata,
                        theta=theta,
                        clip=clip,
                        n_top_genes=n_top_genes,
                        batch_key=batch_key or None,
                        chunksize=chunksize,
                        flavor="pearson_residuals",
                        check_values=False,
                        layer=layer,
                        subset=False,
                        inplace=True,
                    )
                elif flavor in {{"seurat_v3", "seurat_v3_paper"}}:
                    sc.pp.highly_variable_genes(
                        adata,
                        layer=layer,
                        n_top_genes=n_top_genes,
                        span=span,
                        flavor=flavor,
                        subset=False,
                        inplace=True,
                        batch_key=batch_key or None,
                        check_values=False,
                    )
                else:
                    sc.pp.highly_variable_genes(
                        adata,
                        layer=layer,
                        n_top_genes=n_top_genes,
                        n_bins=n_bins,
                        flavor=flavor,
                        subset=False,
                        inplace=True,
                        batch_key=batch_key or None,
                        check_values=True,
                    )
                algorithm = adata.var["highly_variable"].astype(bool).to_numpy(copy=True)
                forced = adata.var_names.isin(requested_keep)
                adata.var["highly_variable_algorithm"] = algorithm
                adata.var["highly_variable_forced"] = forced
                adata.var["highly_variable"] = algorithm | forced
                if subset:
                    adata = adata[:, adata.var["highly_variable"].astype(bool)].copy()
                return adata
            """
        ).strip()
    )


@register_operation("openbio.node.highlyvariablegenes")
def highly_variable_genes(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    expected = {
        "n_top_genes",
        "flavor",
        "source",
        "batch_key",
        "always_keep_genes",
        "subset",
        "theta",
        "clipping_mode",
        "custom_clip",
        "chunksize",
        "span",
        "n_bins",
        "overwrite_existing",
    }
    operation = "Highly Variable Genes"
    adata = _input_anndata(inputs, parameters, operation=operation, expected_parameters=expected)
    _require_in_memory_nonempty(adata, operation=operation)
    _require_unique_axes(adata, operation=operation)
    n_top_genes = _positive_integer(parameters["n_top_genes"], name=f"{operation} n_top_genes")
    flavor = _string(parameters["flavor"], name=f"{operation} flavor")
    if flavor not in LOG_HVG_FLAVORS | COUNT_HVG_FLAVORS:
        raise ProtocolError(f"Unsupported {operation} flavor: {flavor!r}")
    if flavor != "pearson_residuals" and adata.n_obs < 2:
        raise ValueError(f"{operation} flavor {flavor!r} backend requires at least two cells.")
    batch_key = _string(parameters["batch_key"], name=f"{operation} batch_key", allow_empty=True)
    keep_text = _string(parameters["always_keep_genes"], name=f"{operation} always_keep_genes", allow_empty=True)
    subset = _boolean(parameters["subset"], name=f"{operation} subset")
    overwrite_existing = _boolean(parameters["overwrite_existing"], name=f"{operation} overwrite_existing")
    expression = _source_parameter(parameters, HVG_SOURCE, adata)
    source_label = _expression_source_label(expression)
    report_warnings: list[str] = []
    if n_top_genes > adata.n_vars:
        report_warnings.append(
            f"{operation} requested {n_top_genes} genes from {int(adata.n_vars)} features; Scanpy will return "
            "at most the eligible feature count."
        )
    if flavor in LOG_HVG_FLAVORS:
        matrix, cell_totals, source_warnings = _validate_logged_source(
            adata,
            expression,
            operation=operation,
        )
        gene_totals, gene_nonzero = matrix_totals_and_nonzero(matrix, axis=0)
        report_warnings.extend(source_warnings)
        count_validation = False
    else:
        matrix, cell_totals, gene_totals, source_warnings = _validate_count_model_source(
            adata,
            expression,
            operation=operation,
            require_nonzero_totals=False,
        )
        report_warnings.extend(source_warnings)
        _, gene_nonzero = matrix_totals_and_nonzero(matrix, axis=0)
        count_validation = True
    science = dependencies.require_scientific_dependencies()
    gene_totals = science.np.asarray(gene_totals, dtype=float)
    gene_nonzero = science.np.asarray(gene_nonzero, dtype=int)
    expressed_features = int((gene_totals > 0).sum())
    group_sizes, batch_warnings = _batch_group_sizes(
        adata,
        batch_key,
        operation=operation,
        minimum_group_size=1 if flavor == "pearson_residuals" else 2,
    )
    report_warnings.extend(batch_warnings)
    requested_keep = list(dict.fromkeys(gene.strip() for gene in keep_text.split(",") if gene.strip()))
    present_keep = [gene for gene in requested_keep if gene in adata.var_names]
    missing_keep = [gene for gene in requested_keep if gene not in adata.var_names]
    if missing_keep:
        report_warnings.append(f"Requested keep genes were not found: {missing_keep}")
    prior_columns = [column for column in HVG_RESULT_COLUMNS if column in adata.var]
    prior_uns = "hvg" in adata.uns
    if (prior_columns or prior_uns) and not overwrite_existing:
        raise ValueError(
            f"{operation} annotations already exist ({prior_columns!r}); enable overwrite_existing to replace them."
        )
    adata.var.drop(columns=prior_columns, inplace=True)
    adata.uns.pop("hvg", None)

    resolved_theta: float | None = None
    scanpy_clip: float | None = None
    resolved_clip: float | str | None = None
    resolved_chunksize: int | None = None
    resolved_span: float | None = None
    resolved_n_bins: int | None = None
    if flavor == "pearson_residuals":
        resolved_theta = _number(parameters["theta"], name=f"{operation} theta", positive=True)
        clipping_mode = _string(parameters["clipping_mode"], name=f"{operation} clipping_mode")
        if clipping_mode == "sqrt_n_obs":
            scanpy_clip = float(science.np.sqrt(adata.n_obs))
            resolved_clip = scanpy_clip
        elif clipping_mode == "custom":
            scanpy_clip = _number(parameters["custom_clip"], name=f"{operation} custom_clip", nonnegative=True)
            resolved_clip = scanpy_clip
        elif clipping_mode == "none":
            resolved_clip = "infinity"
        else:
            raise ProtocolError(f"Unsupported {operation} clipping_mode: {clipping_mode!r}")
        resolved_chunksize = _positive_integer(parameters["chunksize"], name=f"{operation} chunksize")
    elif flavor in {"seurat_v3", "seurat_v3_paper"}:
        resolved_span = _number(parameters["span"], name=f"{operation} span", positive=True)
        if resolved_span > 1:
            raise ValueError(f"{operation} span must be no greater than one.")
        clipping_mode = None
    else:
        resolved_n_bins = _positive_integer(parameters["n_bins"], name=f"{operation} n_bins")
        clipping_mode = None

    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    try:
        if flavor == "pearson_residuals":
            science.sc.experimental.pp.highly_variable_genes(
                adata,
                theta=resolved_theta,
                clip=scanpy_clip,
                n_top_genes=n_top_genes,
                batch_key=batch_key or None,
                chunksize=resolved_chunksize,
                flavor="pearson_residuals",
                check_values=False,
                layer=expression.scanpy_layer,
                subset=False,
                inplace=True,
            )
        elif flavor in {"seurat_v3", "seurat_v3_paper"}:
            science.sc.pp.highly_variable_genes(
                adata,
                layer=expression.scanpy_layer,
                n_top_genes=n_top_genes,
                span=resolved_span,
                flavor=flavor,
                subset=False,
                inplace=True,
                batch_key=batch_key or None,
                check_values=False,
            )
        else:
            science.sc.pp.highly_variable_genes(
                adata,
                layer=expression.scanpy_layer,
                n_top_genes=n_top_genes,
                n_bins=resolved_n_bins,
                flavor=flavor,
                subset=False,
                inplace=True,
                batch_key=batch_key or None,
                check_values=True,
            )
    except ImportError as error:
        if flavor in {"seurat_v3", "seurat_v3_paper"}:
            raise ImportError(
                "Highly Variable Genes seurat_v3 flavors require scikit-misc; install the recommended dependencies."
            ) from error
        raise
    algorithm_mask = adata.var["highly_variable"].astype(bool).to_numpy(copy=True)
    forced_mask = adata.var_names.isin(present_keep)
    final_mask = algorithm_mask | forced_mask
    adata.var["highly_variable_algorithm"] = algorithm_mask
    adata.var["highly_variable_forced"] = forced_mask
    adata.var["highly_variable"] = final_mask
    algorithm_count = int(algorithm_mask.sum())
    forced_added = int((~algorithm_mask & forced_mask).sum())
    final_count = int(final_mask.sum())
    metric_columns = [column for column in HVG_RESULT_COLUMNS if column in adata.var]
    metric_summaries = {
        column: summarize_numeric(adata.var[column].to_numpy())
        for column in metric_columns
        if column not in {"highly_variable", "highly_variable_algorithm", "highly_variable_forced"}
        and science.np.issubdtype(adata.var[column].dtype, science.np.number)
    }
    if "highly_variable_rank" in adata.var:
        rank = adata.var.loc[algorithm_mask, "highly_variable_rank"]
        top_features = [str(name) for name in rank.dropna().sort_values().index[: min(20, algorithm_count)]]
    else:
        top_features = [str(name) for name in adata.var_names[algorithm_mask][:20]]
    nbatches_distribution: dict[str, int] = {}
    if "highly_variable_nbatches" in adata.var:
        nbatches_distribution = {
            str(key): int(value)
            for key, value in adata.var["highly_variable_nbatches"].value_counts(dropna=False).sort_index().items()
        }
    intersection_count: int | None = None
    intersection_rate: float | None = None
    if "highly_variable_intersection" in adata.var:
        intersection_count = int(adata.var["highly_variable_intersection"].astype(bool).sum())
        intersection_rate = intersection_count / genes
    if subset:
        # AnnData slicing is a view; materialization is required to publish an independent subset artifact.
        adata = adata[:, final_mask].copy()
    report_parameters = {
        "n_top_genes": n_top_genes,
        "flavor": flavor,
        **expression.parameters(),
        "batch_key": batch_key,
        "always_keep_genes": requested_keep,
        "subset": subset,
        "theta": resolved_theta,
        "clipping_mode": clipping_mode,
        "resolved_clip": resolved_clip,
        "chunksize": resolved_chunksize,
        "span": resolved_span,
        "n_bins": resolved_n_bins,
        "check_values": flavor in LOG_HVG_FLAVORS,
        "overwrite_existing": overwrite_existing,
    }
    finish_adata(
        adata,
        "highly_variable_genes",
        report_parameters,
        cells,
        genes,
        started_at,
        warnings=report_warnings,
    )
    reference = HVG_REFERENCE_BY_FLAVOR[flavor]
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellHighlyVariableGenes",
        title="Highly variable gene selection summary",
        operation="highly_variable_genes",
        methods=(
            f"Selected up to {n_top_genes:,} algorithmic highly variable genes from source {source_label!r} "
            f"using flavor {flavor!r}; forced genes were applied after the algorithmic mask."
        ),
        results=(
            f"The algorithm selected {algorithm_count:,} genes; {forced_added:,} requested genes were added "
            f"to the final set of {final_count:,}. The returned AnnData has {int(adata.n_vars):,} features."
        ),
        key_results={
            "input_cells": cells,
            "input_features": genes,
            "output_cells": int(adata.n_obs),
            "output_features": int(adata.n_vars),
            "source": source_label,
            "source_matrix": _matrix_description(matrix),
            "count_model_validation": count_validation,
            "input_cell_totals": summarize_numeric(cell_totals),
            "input_gene_totals": summarize_numeric(gene_totals),
            "expressed_features": expressed_features,
            "cells_per_expressed_feature": summarize_numeric(gene_nonzero[gene_totals > 0]),
            "requested_algorithm_count": n_top_genes,
            "algorithm_selected_count": algorithm_count,
            "forced_requested_count": len(requested_keep),
            "forced_present_count": int(forced_mask.sum()),
            "forced_added_count": forced_added,
            "final_selected_count": final_count,
            "selection_rate": final_count / genes,
            "missing_forced_genes": missing_keep,
            "top_algorithm_features": top_features,
            "metric_summaries": metric_summaries,
            "batch_key": batch_key or None,
            "batch_group_sizes": group_sizes,
            "highly_variable_nbatches_distribution": nbatches_distribution,
            "highly_variable_intersection_count": intersection_count,
            "highly_variable_intersection_rate": intersection_rate,
            "prior_result_columns_replaced": prior_columns,
            "prior_hvg_uns_replaced": prior_uns,
            "raw_preserved": True,
            "subset_applied": subset,
        },
        parameters=report_parameters,
        references=(reference, SCANPY_REFERENCE, ANNDATA_REFERENCE),
        software_packages=PREPROCESS_SOFTWARE_PACKAGES,
        warnings=report_warnings,
        limitations=(
            "Highly variable genes are variation features, not marker genes or evidence of differential expression.",
            "Forced genes are user overrides and are excluded from the algorithm-selected count.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_highly_variable_genes_code(
            expression,
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key,
            requested_keep=requested_keep,
            subset=subset,
            theta=resolved_theta,
            clip=scanpy_clip,
            chunksize=resolved_chunksize,
            span=resolved_span,
            n_bins=resolved_n_bins,
            overwrite_existing=overwrite_existing,
        ),
    )
    return _records(context, adata, report, code)


def hvg_selection_plot_owned(adata: Any) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = render_hvg_selection_plot(adata)
    code = hvg_selection_plot_code()
    flavor = details["flavor"]
    counts = details["selection_counts"]
    warnings = list(details["warnings"])
    title = f"HVG Selection ({flavor})"
    plotted = make_plot_result(
        title=title,
        operation="hvg_selection_plot",
        parameters={},
        description=(
            f"Read-only rendering of stored {flavor!r} feature-selection evidence for "
            f"{details['input_features']:,} input features."
        ),
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        png=png,
    )
    reference = HVG_REFERENCE_BY_FLAVOR[flavor]
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellHVGSelectionPlot",
        title=title,
        operation="hvg_selection_plot",
        methods=(
            f"Read flavor {flavor!r} from adata.uns['hvg'] and plotted stored feature means against "
            f"{details['variability_columns']!r} without recomputing feature-selection statistics."
        ),
        results=(
            f"Classified {counts['selected']:,} selected and {counts['unselected']:,} unselected features; "
            "the algorithm-selected versus forced breakdown is unavailable."
            if counts['algorithm_selected'] is None else
            f"Classified {counts['algorithm_selected']:,} algorithm-selected, {counts['forced']:,} forced, "
            f"and {counts['unselected']:,} unselected features."
        ),
        key_results=details,
        parameters={},
        references=(reference, SCANPY_REFERENCE),
        software_packages=(*PREPROCESS_SOFTWARE_PACKAGES, "matplotlib"),
        warnings=warnings,
        limitations=(
            "This diagnostic displays stored feature-selection statistics and does not rerun or validate the "
            "biological suitability of the selected flavor.",
            "Forced features are user overrides, not algorithm-selected highly variable genes.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


@register_operation("openbio.node.hvgselectionplot")
def hvg_selection_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    operation = "HVG Selection Plot"
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(parameters, set(), operation=operation)
    plotted, report, code = hvg_selection_plot_owned(read_anndata_input(inputs))
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


def _scale_code(
    expression: ExpressionSource,
    *,
    zero_center: bool,
    max_value: float | None,
    output_layer: str,
    overwrite_existing: bool,
    max_dense_gib: float,
) -> str:
    source_layer = expression.layer_name if expression.kind == "layer" else None
    return (
        _generated_feature_helpers()
        + "\n\n\n"
        + dedent(
            f"""
            def scale_expression_to_layer(adata):
                operation = "Scale"
                _openbio_require_feature_input(adata, operation)
                source_kind = {expression.kind!r}
                source_layer = {source_layer!r}
                zero_center = {zero_center!r}
                max_value = {_python_float_literal(max_value)}
                output_layer = {output_layer!r}
                overwrite_existing = {overwrite_existing!r}
                max_dense_gib = {max_dense_gib!r}
                matrix, _ = _openbio_feature_source(adata, source_kind, source_layer, operation)
                output_layer = _openbio_validate_destination(
                    adata,
                    output_layer,
                    overwrite_existing,
                    operation,
                )
                if sparse.issparse(matrix) and zero_center:
                    _openbio_validate_dense_budget(adata.shape, max_dense_gib, operation)
                    warnings.warn(
                        "Zero-centering a sparse expression matrix densifies it; working arrays coexist.",
                        UserWarning,
                        stacklevel=2,
                    )
                work = AnnData(X=matrix.copy())
                sc.pp.scale(work, zero_center=zero_center, max_value=max_value, copy=False)
                if sparse.issparse(matrix) and not zero_center and not sparse.issparse(work.X):
                    raise RuntimeError("Scale unexpectedly densified sparse input with zero_center=False.")
                adata.layers[output_layer] = work.X
                return adata
            """
        ).strip()
    )


@register_operation("openbio.node.scale")
def scale(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    expected = {
        "source",
        "zero_center",
        "clipping_mode",
        "custom_max_value",
        "output_layer",
        "overwrite_existing",
        "max_dense_gib",
    }
    operation = "Scale"
    adata = _input_anndata(inputs, parameters, operation=operation, expected_parameters=expected)
    _require_in_memory_nonempty(adata, operation=operation)
    if adata.n_obs < 2:
        raise ValueError(f"{operation} requires at least two cells because Scanpy uses sample variance (ddof=1).")
    _require_unique_axes(adata, operation=operation)
    zero_center = _boolean(parameters["zero_center"], name=f"{operation} zero_center")
    expression = _source_parameter(parameters, SCALE_SOURCE, adata)
    matrix = expression.matrix(adata)
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} source must align to the current AnnData shape.")
    _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
    output_layer, overwrite_existing, replaced_existing = _validate_destination_layer(
        adata,
        parameters["output_layer"],
        parameters["overwrite_existing"],
        operation=operation,
    )
    clipping_mode = _string(parameters["clipping_mode"], name=f"{operation} clipping_mode")
    if clipping_mode not in {"custom", "none"}:
        raise ProtocolError(f"Unsupported {operation} clipping_mode: {clipping_mode!r}")
    resolved_max_value = (
        _number(parameters["custom_max_value"], name=f"{operation} custom_max_value")
        if clipping_mode == "custom"
        else None
    )
    science = dependencies.require_scientific_dependencies()
    max_dense_gib = _number(parameters["max_dense_gib"], name=f"{operation} max_dense_gib", positive=True)
    estimated_dense_gib = _dense_array_gib(tuple(adata.shape))
    input_sparse = bool(science.sparse.issparse(matrix))
    report_warnings: list[str] = []
    if resolved_max_value is not None and resolved_max_value < 0:
        report_warnings.append("A negative clipping bound was retained; output values follow Scanpy's clipping behavior.")
    if input_sparse and zero_center:
        _validate_dense_budget(tuple(adata.shape), parameters["max_dense_gib"], operation=operation)
        report_warnings.append(
            "Zero-centering a sparse expression matrix densifies it; the reported float64 allocation is a "
            "lower bound because working arrays coexist."
        )
    if input_sparse and not zero_center and resolved_max_value is not None:
        report_warnings.append(
            "With zero_center=False, Scanpy applies only an upper clipping bound; lower-tail values are not "
            "symmetrically truncated."
        )
    feature_mean, feature_std, constant_mask = _scale_feature_statistics(matrix)
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    input_description = _matrix_description(matrix)
    raw_present = adata.raw is not None
    # Scanpy scales X; isolate the transform until the named output is assigned.
    work = science.ad.AnnData(X=matrix.copy())
    science.sc.pp.scale(work, zero_center=zero_center, max_value=resolved_max_value, copy=False)
    if input_sparse and not zero_center and not science.sparse.issparse(work.X):
        raise RuntimeError("Scale unexpectedly densified sparse input with zero_center=False.")
    scaled = work.X
    adata.layers[output_layer] = scaled
    if resolved_max_value is None:
        boundary_count = 0
    else:
        scaled_values = _matrix_values(scaled)
        upper = science.np.isclose(scaled_values, resolved_max_value)
        if zero_center:
            boundary_count = int((upper | science.np.isclose(scaled_values, -resolved_max_value)).sum())
        else:
            boundary_count = int(upper.sum())
    report_parameters = {
        **expression.parameters(),
        "zero_center": zero_center,
        "clipping_mode": clipping_mode,
        "resolved_max_value": resolved_max_value,
        "ddof": 1,
        "output_layer": output_layer,
        "overwrite_existing": overwrite_existing,
        "max_dense_gib": max_dense_gib,
    }
    finish_adata(
        adata,
        "scale_to_layer",
        report_parameters,
        cells,
        genes,
        started_at,
        warnings=report_warnings,
    )
    source_label = _expression_source_label(expression)
    constant_names = [str(name) for name in adata.var_names[constant_mask][:20]]
    scaling_method = "centered z-scaling" if zero_center else "uncentered variance scaling"
    clipping_text = "disabled" if resolved_max_value is None else f"bounded at {resolved_max_value!r}"
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellScale",
        title="Gene scaling summary",
        operation="scale_to_layer",
        methods=(
            f"Applied {scaling_method} to source {source_label!r} with scanpy.pp.scale using per-gene sample "
            f"standard deviations (ddof=1); clipping was {clipping_text}. The result was written to layer "
            f"{output_layer!r}, preserving X, every other layer, and Raw."
        ),
        results=(
            f"Created scaled layer {output_layer!r} for {cells:,} cells and {genes:,} genes. "
            f"{int(constant_mask.sum()):,} constant genes followed Scanpy's unit-denominator behavior, and "
            f"{boundary_count:,} stored values were at a configured clipping boundary."
        ),
        key_results={
            "cells": cells,
            "features": genes,
            "source": source_label,
            "input_matrix": input_description,
            "feature_means_before_scaling": summarize_numeric(feature_mean),
            "feature_standard_deviations_ddof1": summarize_numeric(feature_std),
            "constant_feature_count": int(constant_mask.sum()),
            "constant_feature_examples": constant_names,
            "zero_center": zero_center,
            "clipping_mode": clipping_mode,
            "resolved_max_value": resolved_max_value,
            "values_at_clipping_boundary": boundary_count,
            "output_layer": output_layer,
            "output_matrix": _matrix_description(scaled),
            "estimated_dense_array_gib": estimated_dense_gib,
            "dense_budget_gib": max_dense_gib,
            "sparse_to_dense": input_sparse and zero_center,
            "sparse_preserved": input_sparse and not zero_center and science.sparse.issparse(scaled),
            "replaced_existing_layer": replaced_existing,
            "x_preserved": True,
            "source_preserved": expression.kind != "layer" or output_layer != expression.layer_name,
            "raw_preserved": True,
            "raw_present": raw_present,
        },
        parameters=report_parameters,
        references=(SCANPY_REFERENCE, ANNDATA_REFERENCE),
        software_packages=PREPROCESS_SOFTWARE_PACKAGES,
        warnings=report_warnings,
        limitations=(
            "Scaling changes feature units but does not normalize library size, integrate batches, or provide "
            "statistical inference.",
            "Constant-gene standard deviations are replaced with one internally by Scanpy; such genes are retained.",
            "The dense-memory estimate covers one float64 matrix and is lower than possible peak memory.",
        ),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_scale_code(
            expression,
            zero_center=zero_center,
            max_value=resolved_max_value,
            output_layer=output_layer,
            overwrite_existing=overwrite_existing,
            max_dense_gib=max_dense_gib,
        ),
    )
    return _records(context, adata, report, code)


__all__ = [
    "hvg_selection_plot",
    "hvg_selection_plot_owned",
    "highly_variable_genes",
    "log1p",
    "normalize_to_layer",
    "normalize_total",
    "pearson_residuals_to_layer",
    "scale",
]
