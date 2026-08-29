from __future__ import annotations

import io as bytes_io
import time
from typing import TYPE_CHECKING, Any

from . import dependencies
from .contracts import (
    PlotResult,
    SummaryResult,
    TableResult,
    make_analysis_source,
    record_history,
)

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


def matrix_totals_and_nonzero(matrix: Any, axis: int) -> tuple[Any, Any]:
    science = dependencies.require_scientific_dependencies()
    if science.sparse.issparse(matrix):
        totals = science.np.asarray(matrix.sum(axis=axis)).ravel()
        # ``getnnz`` counts stored entries, including explicit zeros.  Scanpy's
        # ``n_genes_by_counts`` / ``n_cells_by_counts`` semantics count values
        # that are actually non-zero, so use the value-aware sparse operation.
        # The one-shot worker owns this matrix, so let SciPy canonicalize
        # duplicate sparse entries in place instead of allocating a full copy.
        nonzero = science.np.asarray(matrix.count_nonzero(axis=axis)).ravel()
    else:
        array = science.np.asarray(matrix)
        totals = array.sum(axis=axis)
        nonzero = science.np.count_nonzero(array, axis=axis)
    return totals, nonzero


def validate_count_expression(
    matrix: Any,
    *,
    source_label: str,
    require_integers: bool = False,
    require_positive: bool = False,
    require_nonnegative: bool = True,
) -> list[str]:
    """Validate a count-dependent matrix with operation-specific strictness."""
    science = dependencies.require_scientific_dependencies()
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    values = science.np.asarray(values)
    if values.size and not bool(science.np.isfinite(values).all()):
        raise ValueError(f"{source_label} contains non-finite expression values.")
    warnings = []
    if values.size and bool((values < 0).any()):
        message = f"{source_label} contains negative expression values and is not an unnormalized count matrix."
        if require_nonnegative:
            raise ValueError(message)
        warnings.append(f"{message} Results are reported as expression sums rather than verified counts.")
    if values.size == 0 or not bool((values > 0).any()):
        message = f"{source_label} contains no positive expression values."
        if require_positive:
            raise ValueError(message)
        warnings.append(message)
    elif not bool(science.np.allclose(values, science.np.rint(values), rtol=0.0, atol=1e-8)):
        message = f"{source_label} contains non-integer values and is not an unnormalized count matrix."
        if require_integers:
            raise ValueError(message)
        warnings.append(f"{message} Count totals should be interpreted as expression sums.")
    return warnings


def finish_adata(
    adata: AnnData,
    operation: str,
    parameters: dict[str, Any],
    input_cells: int,
    input_genes: int,
    started_at: float,
    *,
    random_seed: int = 0,
    warnings: list[str] | None = None,
) -> AnnData:
    record_history(
        adata,
        operation,
        parameters,
        input_cells,
        input_genes,
        random_seed,
        time.perf_counter() - started_at,
        warnings,
    )
    return adata


def _result_fields(
    *,
    title: str,
    operation: str,
    parameters: dict[str, Any],
    description: str,
    warnings: list[str],
    input_cells: int,
    input_genes: int,
    started_at: float,
    random_seed: int = 0,
) -> dict[str, Any]:
    elapsed = time.perf_counter() - started_at
    source = make_analysis_source(
        operation,
        parameters,
        input_cells,
        input_genes,
        random_seed,
        elapsed,
    )
    return {
        "title": title,
        "parameters": parameters,
        "description": description,
        "warnings": warnings,
        "input_cells": input_cells,
        "input_genes": input_genes,
        "random_seed": random_seed,
        "elapsed_seconds": elapsed,
        "source": source,
    }


def make_summary_result(
    *,
    summary: dict[str, Any],
    title: str,
    operation: str,
    parameters: dict[str, Any],
    description: str,
    warnings: list[str],
    input_cells: int,
    input_genes: int,
    started_at: float,
    random_seed: int = 0,
) -> SummaryResult:
    return SummaryResult(
        summary=summary,
        **_result_fields(
            title=title,
            operation=operation,
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=input_genes,
            started_at=started_at,
            random_seed=random_seed,
        ),
    )


def make_table_result(
    *,
    table: DataFrame,
    title: str,
    operation: str,
    parameters: dict[str, Any],
    description: str,
    warnings: list[str],
    input_cells: int,
    input_genes: int,
    started_at: float,
    random_seed: int = 0,
) -> TableResult:
    return TableResult(
        table=table,
        **_result_fields(
            title=title,
            operation=operation,
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=input_genes,
            started_at=started_at,
            random_seed=random_seed,
        ),
    )


def make_plot_result(
    *,
    png: bytes,
    title: str,
    operation: str,
    parameters: dict[str, Any],
    description: str,
    warnings: list[str],
    input_cells: int,
    input_genes: int,
    started_at: float,
    random_seed: int = 0,
) -> PlotResult:
    return PlotResult(
        png=png,
        **_result_fields(
            title=title,
            operation=operation,
            parameters=parameters,
            description=description,
            warnings=warnings,
            input_cells=input_cells,
            input_genes=input_genes,
            started_at=started_at,
            random_seed=random_seed,
        ),
    )


def figure_to_png(figure: Any) -> bytes:
    science = dependencies.require_scientific_dependencies()
    science.FigureCanvasAgg(figure)
    buffer = bytes_io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    return buffer.getvalue()
