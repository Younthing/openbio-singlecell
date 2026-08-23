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
        nonzero = science.np.asarray(matrix.getnnz(axis=axis)).ravel()
    else:
        array = science.np.asarray(matrix)
        totals = array.sum(axis=axis)
        nonzero = science.np.count_nonzero(array, axis=axis)
    return totals, nonzero


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
    summary: Any,
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
