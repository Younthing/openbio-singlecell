from __future__ import annotations

import io as bytes_io
import time
from typing import TYPE_CHECKING, Any

from . import dependencies
from .contracts import ResultKind, SingleCellResult, make_analysis_source, record_history

if TYPE_CHECKING:
    from anndata import AnnData


def matrix_totals_and_nonzero(matrix: Any, axis: int) -> tuple[Any, Any]:
    if dependencies.sparse.issparse(matrix):
        totals = dependencies.np.asarray(matrix.sum(axis=axis)).ravel()
        nonzero = dependencies.np.asarray(matrix.getnnz(axis=axis)).ravel()
    else:
        array = dependencies.np.asarray(matrix)
        totals = array.sum(axis=axis)
        nonzero = dependencies.np.count_nonzero(array, axis=axis)
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


def make_result(
    *,
    kind: ResultKind,
    title: str,
    operation: str,
    parameters: dict[str, Any],
    description: str,
    warnings: list[str],
    input_cells: int,
    input_genes: int,
    started_at: float,
    random_seed: int = 0,
    table: Any = None,
    png: bytes | None = None,
    summary: Any = None,
) -> SingleCellResult:
    elapsed = time.perf_counter() - started_at
    source = make_analysis_source(
        operation,
        parameters,
        input_cells,
        input_genes,
        random_seed,
        elapsed,
    )
    return SingleCellResult(
        kind=kind,
        title=title,
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=input_cells,
        input_genes=input_genes,
        random_seed=random_seed,
        elapsed_seconds=elapsed,
        source=source,
        table=table,
        png=png,
        summary=summary,
    )


def figure_to_png(figure: Any) -> bytes:
    dependencies.FigureCanvasAgg(figure)
    buffer = bytes_io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    return buffer.getvalue()
