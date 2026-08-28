from __future__ import annotations

import time
import warnings as python_warnings
from collections.abc import Mapping
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, matrix_totals_and_nonzero, validate_count_expression
from .expression_source import DynamicExpressionSource, ExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, SummaryResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/preprocessing"
PREPROCESS_SOFTWARE_PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy")
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


def _require_in_memory_nonempty(adata: AnnData, *, operation: str) -> None:
    if bool(adata.isbacked):
        raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"{operation} requires at least one cell and one feature.")


def _require_unique_axes(adata: AnnData, *, operation: str) -> None:
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    if not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique feature identifiers.")


def _expression_source_label(expression: ExpressionSource) -> str:
    return f"layer:{expression.layer_name}" if expression.kind == "layer" else expression.kind


def _matrix_values(matrix: Any) -> Any:
    science = dependencies.require_scientific_dependencies()
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    return science.np.asarray(values)


def _matrix_description(matrix: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    total_entries = int(matrix.shape[0] * matrix.shape[1])
    if science.sparse.issparse(matrix):
        storage = f"sparse_{matrix.format}"
        nonzero_entries = int(matrix.copy().count_nonzero())
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


def _matrix_is_integer_like(matrix: Any) -> bool:
    science = dependencies.require_scientific_dependencies()
    values = _matrix_values(matrix)
    return bool(values.size == 0 or science.np.allclose(values, science.np.rint(values), rtol=0.0, atol=1e-8))


def _feature_expression_state(adata: AnnData, expression: ExpressionSource) -> tuple[str, str | None]:
    """Resolve the latest known state for each expression representation."""
    x_state = "logged_unverified" if isinstance(adata.uns.get("log1p"), Mapping) else "unknown"
    x_evidence = "AnnData uns['log1p'] marker without normalization provenance" if x_state != "unknown" else None
    layer_states: dict[str, tuple[str, str]] = {}
    for entry in _history_entries(adata):
        operation = entry.get("operation")
        parameters = entry.get("parameters")
        parameters = parameters if isinstance(parameters, Mapping) else {}

        if operation == "snapshot_expression":
            layer_states["counts"] = ("counts", "OpenBio Snapshot Expression history")
            if parameters.get("source") == "X":
                x_state, x_evidence = "counts", "OpenBio Snapshot Expression history"
        elif operation == "normalize_to_layer":
            output_layer = parameters.get("output_layer")
            transform = parameters.get("transform")
            if isinstance(output_layer, str):
                state = {"log1p": "logged", "none": "normalized", "sqrt": "transformed"}.get(
                    transform,
                    "derived",
                )
                layer_states[output_layer] = (state, "OpenBio Normalize to Layer history")
        elif operation == "pearson_residuals_to_layer":
            output_layer = parameters.get("output_layer")
            if isinstance(output_layer, str):
                layer_states[output_layer] = ("pearson_residuals", "OpenBio Pearson Residuals history")
        elif operation == "scale_to_layer":
            output_layer = parameters.get("output_layer")
            if isinstance(output_layer, str):
                layer_states[output_layer] = ("scaled", "OpenBio Scale history")

        if operation == "normalize_total":
            x_state, x_evidence = "normalized", "OpenBio Normalize Total history"
        elif operation == "log1p":
            if x_state == "normalized":
                x_state, x_evidence = "logged", "OpenBio Normalize Total followed by Log1p history"
            else:
                x_state, x_evidence = (
                    "logged_unverified",
                    "OpenBio Log1p history without proven prior normalization",
                )

    if expression.kind == "layer":
        return layer_states.get(expression.layer_name, ("unknown", None))
    return x_state, x_evidence


def _validate_finite_numeric_matrix(matrix: Any, *, source_label: str) -> Any:
    science = dependencies.require_scientific_dependencies()
    if tuple(getattr(matrix, "shape", ())) == () or len(matrix.shape) != 2:
        raise TypeError(f"{source_label} must be a two-dimensional numeric matrix.")
    if not (science.sparse.issparse(matrix) or isinstance(matrix, science.np.ndarray)):
        raise TypeError(f"{source_label} must be an in-memory NumPy or SciPy sparse matrix.")
    values = _matrix_values(matrix)
    if values.size and not bool(science.np.issubdtype(values.dtype, science.np.number)):
        raise TypeError(f"{source_label} must contain numeric values.")
    if values.size and not bool(science.np.isfinite(values).all()):
        raise ValueError(f"{source_label} contains non-finite expression values.")
    return values


def _validate_strict_count_source(
    adata: AnnData,
    expression: ExpressionSource,
    *,
    operation: str,
    reject_zero_genes: bool,
) -> tuple[Any, Any, Any, list[str], str, str | None]:
    science = dependencies.require_scientific_dependencies()
    matrix = expression.matrix(adata)
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} source must align to the current AnnData shape.")
    values = _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
    state, evidence = _feature_expression_state(adata, expression)
    warnings: list[str] = []
    if values.size and bool((values < 0).any()):
        raise ValueError(f"{operation} source contains negative values outside the non-negative count-model domain.")
    if values.size == 0 or not bool((values > 0).any()):
        raise ValueError(f"{operation} source contains no positive expression values.")
    if not _matrix_is_integer_like(matrix):
        warnings.append(
            f"{operation} source contains fractional non-negative values. The selected count-based method is "
            "executable, but the input cannot be described as verified integer UMI counts."
        )
    if state != "counts":
        detail = (
            f"provenance describes the source as {state} expression ({evidence})"
            if state != "unknown"
            else "OpenBio provenance cannot prove the selected expression state"
        )
        warnings.append(
            f"{detail}; honoring the expert-selected source without claiming verified raw UMI counts."
        )
    cell_totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
    gene_totals, _ = matrix_totals_and_nonzero(matrix, axis=0)
    cell_totals = science.np.asarray(cell_totals, dtype=float)
    gene_totals = science.np.asarray(gene_totals, dtype=float)
    zero_cells = int((cell_totals <= 0).sum())
    if zero_cells:
        raise ValueError(f"{operation} source contains {zero_cells} cells with zero total counts; filter them first.")
    zero_genes = int((gene_totals <= 0).sum())
    if reject_zero_genes and zero_genes:
        raise ValueError(f"{operation} source contains {zero_genes} genes with zero total counts; filter them first.")
    return matrix, cell_totals, gene_totals, warnings, state, evidence


def _validate_logged_source(
    adata: AnnData,
    expression: ExpressionSource,
    *,
    operation: str,
) -> tuple[Any, Any, list[str], str, str | None]:
    science = dependencies.require_scientific_dependencies()
    matrix = expression.matrix(adata)
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} source must align to the current AnnData shape.")
    values = _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
    if values.size and bool((values < 0).any()):
        raise ValueError(f"{operation} normalized log source contains negative expression values.")
    if values.size == 0 or not bool((values > 0).any()):
        raise ValueError(f"{operation} normalized log source contains no positive expression values.")
    state, evidence = _feature_expression_state(adata, expression)
    warnings: list[str] = []
    if state == "logged_unverified":
        warnings.append(
            "The selected source has a Scanpy log1p marker, but OpenBio provenance cannot prove that "
            "total-count normalization preceded logarithmization; confirm the external data contract."
        )
    elif state != "logged":
        detail = (
            f"provenance describes the source as {state} expression ({evidence})"
            if state != "unknown"
            else "OpenBio provenance cannot prove the source is normalized log1p expression"
        )
        warnings.append(
            f"{detail}; honoring the expert-selected source for the {operation} dispersion flavor."
        )
    totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
    totals = science.np.asarray(totals, dtype=float)
    zero_cells = int((totals <= 0).sum())
    if zero_cells:
        raise ValueError(f"{operation} source contains {zero_cells} cells with zero total expression.")
    return matrix, totals, warnings, state, evidence


def _positive_float(value: Any, *, name: str, allow_infinity: bool = False) -> float:
    science = dependencies.require_scientific_dependencies()
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a number.")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a number.") from error
    valid = number > 0 and (bool(science.np.isfinite(number)) or (allow_infinity and number == float("inf")))
    if not valid:
        suffix = " or positive infinity" if allow_infinity else ""
        raise ValueError(f"{name} must be finite and greater than zero{suffix}.")
    return number


def _nonnegative_float(value: Any, *, name: str) -> float:
    science = dependencies.require_scientific_dependencies()
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a number.")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a number.") from error
    if not bool(science.np.isfinite(number)) or number < 0:
        raise ValueError(f"{name} must be finite and non-negative.")
    return number


def _positive_integer(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer.")
    if value < 1:
        raise ValueError(f"{name} must be at least one.")
    return int(value)


def _resolve_clipping(
    mode: str,
    custom_value: Any,
    *,
    n_obs: int,
    value_name: str,
) -> tuple[str, float | None, float | str]:
    science = dependencies.require_scientific_dependencies()
    if not isinstance(mode, str):
        raise TypeError("clipping_mode must be a string.")
    mode = mode.strip()
    if mode not in {"sqrt_n_obs", "custom", "none"}:
        raise ValueError(f"Unsupported clipping_mode: {mode!r}")
    if mode == "sqrt_n_obs":
        return mode, None, float(science.np.sqrt(n_obs))
    if mode == "none":
        return mode, float("inf"), "infinity"
    custom = _nonnegative_float(custom_value, name=value_name)
    return mode, custom, custom


def _dense_array_gib(shape: tuple[int, int]) -> float:
    return float(int(shape[0]) * int(shape[1]) * 8 / 1024**3)


def _validate_dense_budget(shape: tuple[int, int], max_dense_gib: Any, *, operation: str) -> tuple[float, float]:
    limit = _positive_float(max_dense_gib, name=f"{operation} max_dense_gib")
    estimate = _dense_array_gib(shape)
    if estimate > limit:
        raise ValueError(
            f"{operation} requires an estimated {estimate:.6g} GiB dense array, exceeding max_dense_gib={limit:.6g}."
        )
    return estimate, limit


def _validate_destination_layer(
    adata: AnnData,
    expression: ExpressionSource,
    output_layer: Any,
    overwrite_existing: Any,
    *,
    operation: str,
) -> tuple[str, bool]:
    if not isinstance(output_layer, str):
        raise TypeError(f"{operation} output_layer must be a string.")
    output_layer = output_layer.strip()
    if not output_layer:
        raise ValueError(f"{operation} output_layer cannot be empty.")
    if output_layer == "counts":
        raise ValueError(f"{operation} cannot overwrite the reserved canonical counts layer.")
    if expression.kind == "layer" and output_layer == expression.layer_name:
        raise ValueError(f"{operation} output_layer must differ from the selected source layer.")
    if not isinstance(overwrite_existing, bool):
        raise TypeError(f"{operation} overwrite_existing must be a boolean.")
    replaced = output_layer in adata.layers
    if replaced and not overwrite_existing:
        raise ValueError(f"{operation} output layer already exists: {output_layer!r}")
    return output_layer, replaced


def _batch_group_sizes(
    adata: AnnData,
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
    groups = [{"label": value, "count": int(count)} for value, count in counts.items() if int(count) > 0]
    small = [group for group in groups if group["count"] < minimum_group_size]
    if small:
        raise ValueError(
            f"{operation} backend requires at least {minimum_group_size} cells in every batch group; "
            f"found {small!r}."
        )
    warnings = []
    singleton = [group for group in groups if group["count"] == 1]
    if singleton:
        warnings.append(
            f"Batch-aware selection includes singleton group(s) {singleton!r}; the selected flavor can execute, "
            "but within-group variation is not independently supported."
        )
    if len(groups) == 1:
        warnings.append(f"Batch-aware selection used {batch_key!r}, but the column has only one observed level.")
    return groups, warnings


def _validate_full_gene_raw_snapshot(adata: AnnData, *, operation: str) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    raw = adata.raw
    if raw is None:
        raise ValueError(f"{operation} subset=True requires a full-gene Raw snapshot of post-QC counts.")
    if int(raw.n_obs) != int(adata.n_obs) or not raw.obs_names.equals(adata.obs_names):
        raise ValueError(f"{operation} Raw snapshot observations are not aligned to the current AnnData.")
    if not bool(raw.var_names.is_unique):
        raise ValueError(f"{operation} Raw snapshot requires unique feature identifiers.")
    if int(raw.n_vars) < int(adata.n_vars) or not adata.var_names.isin(raw.var_names).all():
        raise ValueError(f"{operation} Raw snapshot does not preserve the full current feature set.")
    validate_count_expression(
        raw.X,
        source_label=f"{operation} Raw snapshot",
        require_integers=True,
        require_positive=True,
    )
    totals, _ = matrix_totals_and_nonzero(raw.X, axis=1)
    zero_cells = int((science.np.asarray(totals, dtype=float) <= 0).sum())
    if zero_cells:
        raise ValueError(f"{operation} Raw snapshot contains {zero_cells} cells with zero total counts.")
    return {
        "present": True,
        "cells": int(raw.n_obs),
        "features": int(raw.n_vars),
        "integer_like": True,
        "observations_aligned": True,
        "preserves_current_features": True,
    }


def _scale_feature_statistics(matrix: Any) -> tuple[Any, Any, Any]:
    science = dependencies.require_scientific_dependencies()
    n_obs = int(matrix.shape[0])
    if science.sparse.issparse(matrix):
        mean = science.np.asarray(matrix.mean(axis=0), dtype=float).ravel()
        mean_square = science.np.asarray(matrix.power(2).mean(axis=0), dtype=float).ravel()
        variance = science.np.maximum((mean_square - mean**2) * n_obs / (n_obs - 1), 0.0)
    else:
        values = science.np.asarray(matrix, dtype=float)
        mean = values.mean(axis=0)
        variance = values.var(axis=0, ddof=1)
    standard_deviation = science.np.sqrt(variance)
    constant = science.np.isclose(standard_deviation, 0.0, rtol=0.0, atol=1e-12)
    return mean, standard_deviation, constant


def _history_entries(adata: AnnData) -> list[Mapping[str, Any]]:
    metadata = adata.uns.get("openbio_singlecell")
    if not isinstance(metadata, Mapping):
        return []
    history = metadata.get("analysis_history")
    if not isinstance(history, Mapping):
        return []
    return [entry for entry in history.values() if isinstance(entry, Mapping)]


def _expression_state(adata: AnnData, expression: ExpressionSource) -> tuple[str, str | None]:
    return _feature_expression_state(adata, expression)


def _validate_count_source(
    adata: AnnData,
    expression: ExpressionSource,
    *,
    operation: str,
) -> tuple[Any, Any, list[str], str, str | None]:
    matrix = expression.matrix(adata)
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} source must align to the current AnnData shape.")
    values = _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
    state, evidence = _expression_state(adata, expression)
    warnings: list[str] = []
    if state != "counts":
        detail = (
            f"provenance describes the selected source as {state} expression ({evidence})"
            if state != "unknown"
            else "OpenBio provenance cannot prove the selected expression state"
        )
        warnings.append(
            f"{detail}; total-count scaling will honor the expert-selected matrix but must not be interpreted "
            "as verified raw UMI normalization."
        )
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
    science = dependencies.require_scientific_dependencies()
    totals = science.np.asarray(totals, dtype=float)
    zero_total_cells = int((totals == 0).sum())
    if zero_total_cells:
        warnings.append(
            f"{operation} source contains {zero_total_cells} cells with zero total expression; "
            "Scanpy leaves those rows unchanged, so they cannot attain target_sum."
        )
    negative_total_cells = int((totals < 0).sum())
    if negative_total_cells:
        warnings.append(
            f"{operation} source contains {negative_total_cells} cells with negative total expression; "
            "the resulting scale-factor signs require expert interpretation."
        )
    return matrix, totals, warnings, state, evidence


def _positive_target_sum(value: float, *, operation: str) -> float:
    science = dependencies.require_scientific_dependencies()
    if isinstance(value, bool):
        raise TypeError(f"{operation} target_sum must be a number.")
    try:
        target_sum = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{operation} target_sum must be a number.") from error
    if not bool(science.np.isfinite(target_sum)) or target_sum <= 0:
        raise ValueError(f"{operation} target_sum must be finite and greater than zero.")
    return target_sum


def _generated_state_helpers() -> str:
    return dedent(
        """
        import warnings
        from collections.abc import Mapping

        import numpy as np
        import scanpy as sc
        from anndata import AnnData
        from scipy import sparse


        def _openbio_history_entries(adata):
            metadata = adata.uns.get("openbio_singlecell")
            if not isinstance(metadata, Mapping):
                return []
            history = metadata.get("analysis_history")
            if not isinstance(history, Mapping):
                return []
            return [entry for entry in history.values() if isinstance(entry, Mapping)]


        def _openbio_record_operation(adata, operation, parameters):
            existing = adata.uns.get("openbio_singlecell")
            metadata = dict(existing) if isinstance(existing, Mapping) else {}
            existing_history = metadata.get("analysis_history")
            history = dict(existing_history) if isinstance(existing_history, Mapping) else {}
            indices = []
            for key in history:
                try:
                    indices.append(int(key))
                except (TypeError, ValueError):
                    continue
            index = max(indices, default=-1) + 1
            key = f"{index:06d}"
            while key in history:
                index += 1
                key = f"{index:06d}"
            history[key] = {"operation": operation, "parameters": dict(parameters)}
            metadata.setdefault("schema_version", 1)
            metadata.setdefault("version", "0.2.0")
            metadata.setdefault("display_name", "AnnData")
            metadata.setdefault("source", {})
            metadata.setdefault("random_seed", 0)
            metadata.setdefault("warnings", [])
            metadata["analysis_history"] = history
            adata.uns["openbio_singlecell"] = metadata


        def _openbio_expression_state(adata, source_kind, layer_name=None):
            if isinstance(adata.uns.get("log1p"), Mapping):
                x_state = "logged_unverified"
                x_evidence = "AnnData uns['log1p'] marker without normalization provenance"
            else:
                x_state = "unknown"
                x_evidence = None
            layer_states = {}
            for entry in _openbio_history_entries(adata):
                operation = entry.get("operation")
                parameters = entry.get("parameters")
                parameters = parameters if isinstance(parameters, Mapping) else {}
                if operation == "snapshot_expression":
                    layer_states["counts"] = ("counts", "OpenBio Snapshot Expression history")
                    if parameters.get("source") == "X":
                        x_state, x_evidence = "counts", "OpenBio Snapshot Expression history"
                elif operation == "normalize_to_layer":
                    output_layer = parameters.get("output_layer")
                    if isinstance(output_layer, str):
                        state = {
                            "log1p": "logged",
                            "none": "normalized",
                            "sqrt": "transformed",
                        }.get(parameters.get("transform"), "derived")
                        layer_states[output_layer] = (state, "OpenBio Normalize to Layer history")
                elif operation == "pearson_residuals_to_layer":
                    output_layer = parameters.get("output_layer")
                    if isinstance(output_layer, str):
                        layer_states[output_layer] = (
                            "pearson_residuals",
                            "OpenBio Pearson Residuals history",
                        )
                elif operation == "scale_to_layer":
                    output_layer = parameters.get("output_layer")
                    if isinstance(output_layer, str):
                        layer_states[output_layer] = ("scaled", "OpenBio Scale history")

                if operation == "normalize_total":
                    x_state, x_evidence = "normalized", "OpenBio Normalize Total history"
                elif operation == "log1p":
                    if x_state == "normalized":
                        x_state, x_evidence = (
                            "logged",
                            "OpenBio Normalize Total followed by Log1p history",
                        )
                    else:
                        x_state, x_evidence = (
                            "logged_unverified",
                            "OpenBio Log1p history without proven prior normalization",
                        )

            if source_kind == "layer":
                return layer_states.get(layer_name, ("unknown", None))
            return x_state, x_evidence


        def _openbio_validate_count_source(adata, matrix, source_kind, layer_name, operation):
            if tuple(matrix.shape) != tuple(adata.shape):
                raise ValueError(f"{operation} source must align to the current AnnData shape.")
            state, evidence = _openbio_expression_state(adata, source_kind, layer_name)
            values = np.asarray(matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel())
            if values.size and not np.issubdtype(values.dtype, np.number):
                raise TypeError(f"{operation} source must contain numeric values.")
            if values.size and not bool(np.isfinite(values).all()):
                raise ValueError(f"{operation} source contains non-finite expression values.")
            scientific_warnings = []
            if state != "counts":
                detail = (
                    f"provenance describes the selected source as {state} expression ({evidence})"
                    if state != "unknown"
                    else "OpenBio provenance cannot prove the selected expression state"
                )
                scientific_warnings.append(
                    f"{detail}; total-count scaling will honor the expert-selected matrix but must not be "
                    "interpreted as verified raw UMI normalization."
                )
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
    return (
        _generated_state_helpers()
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
                _openbio_validate_count_source(adata, matrix, source_kind, layer_name, operation)
                output = adata.copy()
                work = AnnData(X=matrix.copy())
                with warnings.catch_warnings(record=True):
                    warnings.simplefilter("always")
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
                    raise RuntimeError("Normalize Total produced non-finite expression values.")
                output.X = work.X.copy()
                output.uns.pop("log1p", None)
                _openbio_record_operation(
                    output,
                    "normalize_total",
                    {{
                        "target_sum": target_sum,
                        "exclude_highly_expressed": False,
                        "source": source_kind,
                        "source_layer": layer_name,
                    }},
                )
                return output
            """
        ).strip()
    )


def _log1p_code() -> str:
    return (
        _generated_state_helpers()
        + "\n\n\n"
        + dedent(
            """
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
                state, evidence = _openbio_expression_state(adata, "X")
                if values.size and bool((values < 0).any()):
                    warnings.warn(
                        "Log1p source contains values in (-1, 0); the transform is finite but the matrix is "
                        "not conventional non-negative expression.",
                        UserWarning,
                        stacklevel=2,
                    )
                if state in {"logged", "logged_unverified"}:
                    warnings.warn(
                        f"Log1p provenance indicates an already log-transformed source ({evidence}); applying "
                        "log1p again because the expert explicitly requested it.",
                        UserWarning,
                        stacklevel=2,
                    )
                elif state != "normalized":
                    warnings.warn(
                        f"Log1p source state is {state!r} ({evidence}); applying the requested mathematical "
                        "transform without claiming prior normalization.",
                        UserWarning,
                        stacklevel=2,
                    )
                output = adata.copy()
                sc.pp.log1p(output, base=None)
                output_values = np.asarray(
                    output.X.data if sparse.issparse(output.X) else np.asarray(output.X).ravel()
                )
                if output_values.size and not np.isfinite(output_values).all():
                    raise RuntimeError("Log1p produced non-finite expression values.")
                _openbio_record_operation(
                    output,
                    "log1p",
                    {"source": "X", "base": "natural", "pseudocount": 1.0},
                )
                return output
            """
        ).strip()
    )


def _normalize_to_layer_code(
    expression: ExpressionSource,
    target_sum: float,
    transform: str,
    output_layer: str,
    overwrite_existing: bool,
) -> str:
    layer_name = expression.layer_name if expression.kind == "layer" else None
    return (
        _generated_state_helpers()
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
                if output_layer == "counts":
                    raise ValueError("Normalize to Layer cannot write the reserved canonical counts layer.")
                if source_kind == "layer" and output_layer == source_layer:
                    raise ValueError("Normalize to Layer output_layer must differ from the selected source layer.")
                if output_layer in adata.layers and not overwrite_existing:
                    raise ValueError(f"Normalize to Layer output layer already exists: {{output_layer!r}}")
                if source_kind == "layer":
                    if source_layer not in adata.layers:
                        raise ValueError(f"{{operation}} source layer not found: {{source_layer!r}}")
                    matrix = adata.layers[source_layer]
                else:
                    matrix = adata.X
                _openbio_validate_count_source(adata, matrix, source_kind, source_layer, operation)
                output = adata.copy()
                work = AnnData(X=matrix.copy())
                with warnings.catch_warnings(record=True):
                    warnings.simplefilter("always")
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
                output.layers[output_layer] = work.X.copy()
                _openbio_record_operation(
                    output,
                    "normalize_to_layer",
                    {{
                        "source": source_kind,
                        "source_layer": source_layer,
                        "target_sum": target_sum,
                        "exclude_highly_expressed": False,
                        "transform": transform,
                        "output_layer": output_layer,
                        "overwrite_existing": overwrite_existing,
                    }},
                )
                return output
            """
        ).strip()
    )


class OpenBioSingleCellNormalizeTotal(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Normalize Total source",
        layer_input_id="source_layer",
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellNormalizeTotal",
            display_name="Normalize Total",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Float.Input("target_sum", default=10000.0, step=1000.0),
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
        target_sum: float = 10000.0,
        source: DynamicExpressionSource | None = None,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        _require_in_memory_nonempty(adata, operation="Normalize Total")
        target_sum = _positive_target_sum(target_sum, operation="Normalize Total")
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        matrix, input_totals, warnings, state, evidence = _validate_count_source(
            adata,
            expression,
            operation="Normalize Total",
        )
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        stale_log1p_marker_removed = isinstance(adata.uns.get("log1p"), Mapping)
        output = adata.copy()
        work = science.ad.AnnData(X=matrix.copy())
        with python_warnings.catch_warnings(record=True) as caught:
            python_warnings.simplefilter("always")
            science.sc.pp.normalize_total(
                work,
                target_sum=target_sum,
                exclude_highly_expressed=False,
                inplace=True,
            )
        for item in caught:
            warnings.append(f"Scanpy normalize_total warning: {item.message}")
        output_values = _matrix_values(work.X)
        if output_values.size and not bool(science.np.isfinite(output_values).all()):
            raise RuntimeError("Normalize Total produced non-finite expression values.")
        output.X = work.X.copy()
        output.uns.pop("log1p", None)
        output_totals, _ = matrix_totals_and_nonzero(output.X, axis=1)
        parameters = {
            "target_sum": target_sum,
            "exclude_highly_expressed": False,
            **expression.parameters(),
        }
        finish_adata(output, "normalize_total", parameters, cells, genes, started_at, warnings=warnings)
        source_label = _expression_source_label(expression)
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellNormalizeTotal",
            title="Total-count normalization summary",
            operation="normalize_total",
            methods=(
                f"Each of {cells:,} cells was scaled from expert-selected source {source_label!r} with "
                f"scanpy.pp.normalize_total to a target total of {target_sum:,.6g}. The transformed matrix was "
                "written to X on a copied AnnData; raw and every named layer were preserved."
            ),
            results=(
                f"Total-count-normalized expression was produced for {cells:,} cells and {genes:,} features. "
                f"The median output cell total was {float(science.np.median(output_totals)):,.6g}."
            ),
            key_results={
                "cells": cells,
                "features": genes,
                "source": source_label,
                "source_state": state,
                "source_state_evidence": evidence,
                "input_matrix": _matrix_description(matrix),
                "input_integer_like": _matrix_is_integer_like(matrix),
                "input_cell_totals": summarize_numeric(input_totals),
                "zero_total_cells": int((science.np.asarray(input_totals) == 0).sum()),
                "negative_total_cells": int((science.np.asarray(input_totals) < 0).sum()),
                "target_sum": target_sum,
                "output_matrix": _matrix_description(output.X),
                "output_cell_totals": summarize_numeric(output_totals),
                "stale_log1p_marker_removed": stale_log1p_marker_removed,
                "raw_preserved": True,
                "raw_present": adata.raw is not None,
                "layers_preserved": sorted(str(key) for key in adata.layers.keys() if key is not None),
            },
            parameters=parameters,
            references=(SCANPY_REFERENCE, ZHENG_REFERENCE),
            software_packages=PREPROCESS_SOFTWARE_PACKAGES,
            warnings=warnings,
            limitations=(
                "Total-count scaling adjusts library depth only; it is not Technical batch integration or "
                "Sample-level inference.",
                "The normalized X is derived expression and must not replace the Raw snapshot for count-dependent models.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            code=_normalize_total_code(expression, target_sum),
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellLog1p(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLog1p",
            display_name="Log1p",
            category=CATEGORY,
            inputs=[AnnDataType.Input("adata")],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(cls, adata: AnnData) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        _require_in_memory_nonempty(adata, operation="Log1p")
        matrix = adata.X
        values = _matrix_values(matrix)
        if values.size and not bool(science.np.isfinite(values).all()):
            raise ValueError("Log1p source contains non-finite expression values.")
        if values.size and bool((values <= -1).any()):
            raise ValueError("Log1p source contains values <= -1 outside the finite real log1p domain.")
        expression = ExpressionSource("X")
        state, evidence = _expression_state(adata, expression)
        warnings = []
        if values.size and bool((values < 0).any()):
            warnings.append(
                "Log1p source contains values in (-1, 0); the transform is finite but the matrix is not "
                "conventional non-negative expression."
            )
        if state in {"logged", "logged_unverified"}:
            warnings.append(
                f"Log1p provenance indicates an already log-transformed source ({evidence}); applying log1p "
                "again because the expert explicitly requested it."
            )
        elif state != "normalized":
            warnings.append(
                f"Log1p source state is {state!r} ({evidence}); applying the requested mathematical transform "
                "without claiming prior normalization."
            )
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        input_description = _matrix_description(matrix)
        output = adata.copy()
        science.sc.pp.log1p(output, base=None)
        output_values = _matrix_values(output.X)
        if output_values.size and not bool(science.np.isfinite(output_values).all()):
            raise RuntimeError("Log1p produced non-finite expression values.")
        parameters = {"source": "X", "base": "natural", "pseudocount": 1.0}
        finish_adata(output, "log1p", parameters, cells, genes, started_at, warnings=warnings)
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellLog1p",
            title="Log1p transformation summary",
            operation="log1p",
            methods=(
                f"Applied the natural-log transformation log(1+x) to active X for {cells:,} cells and "
                f"{genes:,} features with scanpy.pp.log1p on a copied AnnData."
            ),
            results=(
                f"Created a log1p expression representation in X for {cells:,} cells while preserving raw and "
                "all named layers unchanged."
            ),
            key_results={
                "cells": cells,
                "features": genes,
                "source": "X",
                "source_state": state,
                "source_state_evidence": evidence,
                "input_matrix": input_description,
                "output_matrix": _matrix_description(output.X),
                "base": "natural",
                "pseudocount": 1.0,
                "raw_preserved": True,
                "raw_present": adata.raw is not None,
                "layers_preserved": sorted(str(key) for key in adata.layers.keys() if key is not None),
            },
            parameters=parameters,
            references=(SCANPY_REFERENCE, TRANSFORM_REFERENCE, ANNDATA_REFERENCE),
            software_packages=PREPROCESS_SOFTWARE_PACKAGES,
            warnings=warnings,
            limitations=(
                "Log1p does not perform library-size normalization, Technical batch integration, or scaling.",
                "Log1p expression is not a count matrix and must not be used as the Raw snapshot for count models.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            code=_log1p_code(),
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellNormalizeToLayer(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Normalize source",
        layer_input_id="source_layer",
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellNormalizeToLayer",
            display_name="Normalize to Layer",
            category=CATEGORY,
            description="Normalize an expression matrix into a named layer without replacing adata.X.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Float.Input("target_sum", default=10000.0, step=1000.0),
                io.Combo.Input("transform", options=["none", "log1p", "sqrt"], default="log1p"),
                io.String.Input("output_layer", default="log1p_norm", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
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
        source: DynamicExpressionSource | None = None,
        target_sum: float = 10000.0,
        transform: str = "log1p",
        output_layer: str = "log1p_norm",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        _require_in_memory_nonempty(adata, operation="Normalize to Layer")
        target_sum = _positive_target_sum(target_sum, operation="Normalize to Layer")
        if not isinstance(transform, str):
            raise TypeError("Normalize to Layer transform must be a string.")
        transform = transform.strip()
        if transform not in {"none", "log1p", "sqrt"}:
            raise ValueError(f"Unsupported Normalize to Layer transform: {transform!r}")
        if not isinstance(output_layer, str):
            raise TypeError("Normalize to Layer output_layer must be a string.")
        output_layer = output_layer.strip()
        if not output_layer:
            raise ValueError("Normalize to Layer output_layer cannot be empty.")
        if output_layer == "counts":
            raise ValueError("Normalize to Layer cannot write the reserved canonical counts layer.")
        if not isinstance(overwrite_existing, bool):
            raise TypeError("Normalize to Layer overwrite_existing must be a boolean.")
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        if expression.kind == "layer" and output_layer == expression.layer_name:
            raise ValueError("Normalize to Layer output_layer must differ from the selected source layer.")
        if output_layer in adata.layers and not overwrite_existing:
            raise ValueError(f"Normalize to Layer output layer already exists: {output_layer!r}")
        matrix, input_totals, warnings, state, evidence = _validate_count_source(
            adata,
            expression,
            operation="Normalize to Layer",
        )

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        work = science.ad.AnnData(X=matrix.copy())
        with python_warnings.catch_warnings(record=True) as caught:
            python_warnings.simplefilter("always")
            science.sc.pp.normalize_total(
                work,
                target_sum=target_sum,
                exclude_highly_expressed=False,
                inplace=True,
            )
        for item in caught:
            warnings.append(f"Scanpy normalize_total warning: {item.message}")
        normalized_totals, _ = matrix_totals_and_nonzero(work.X, axis=1)
        normalized_values = _matrix_values(work.X)
        if normalized_values.size and not bool(science.np.isfinite(normalized_values).all()):
            raise RuntimeError("Normalize to Layer produced non-finite normalized values.")
        if transform == "log1p":
            if normalized_values.size and bool((normalized_values <= -1).any()):
                raise ValueError(
                    "Normalize to Layer log1p transform requires every normalized value to be greater than -1."
                )
            science.sc.pp.log1p(work, base=None)
        elif transform == "sqrt":
            if normalized_values.size and bool((normalized_values < 0).any()):
                raise ValueError("Normalize to Layer sqrt transform requires non-negative normalized values.")
            work.X = work.X.sqrt() if science.sparse.issparse(work.X) else science.np.sqrt(work.X)
        output_values = _matrix_values(work.X)
        if output_values.size and not bool(science.np.isfinite(output_values).all()):
            raise RuntimeError("Normalize to Layer produced non-finite output values.")
        replaced_existing = output_layer in adata.layers
        output.layers[output_layer] = work.X.copy()
        parameters = {
            **expression.parameters(),
            "target_sum": target_sum,
            "exclude_highly_expressed": False,
            "transform": transform,
            "output_layer": output_layer,
            "overwrite_existing": overwrite_existing,
        }
        finish_adata(output, "normalize_to_layer", parameters, cells, genes, started_at, warnings=warnings)
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
                "source_state": state,
                "source_state_evidence": evidence,
                "input_matrix": _matrix_description(matrix),
                "input_integer_like": _matrix_is_integer_like(matrix),
                "input_cell_totals": summarize_numeric(input_totals),
                "zero_total_cells": int((science.np.asarray(input_totals) == 0).sum()),
                "negative_total_cells": int((science.np.asarray(input_totals) < 0).sum()),
                "target_sum": target_sum,
                "intermediate_normalized_cell_totals": summarize_numeric(normalized_totals),
                "transform": transform,
                "output_layer": output_layer,
                "output_matrix": _matrix_description(output.layers[output_layer]),
                "replaced_existing_layer": replaced_existing,
                "x_preserved": True,
                "raw_preserved": True,
                "raw_present": adata.raw is not None,
                "canonical_counts_preserved": "counts" in adata.layers,
            },
            parameters=parameters,
            references=(SCANPY_REFERENCE, TRANSFORM_REFERENCE, ANNDATA_REFERENCE),
            software_packages=PREPROCESS_SOFTWARE_PACKAGES,
            warnings=warnings,
            limitations=(
                "Only the intermediate linear normalized matrix has cell totals equal to target_sum; transformed "
                "layer sums do not.",
                "The derived layer is not counts and must not replace the Raw snapshot in count-dependent models.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            code=_normalize_to_layer_code(
                expression,
                target_sum,
                transform,
                output_layer,
                overwrite_existing,
            ),
        )
        return io.NodeOutput(output, report, code)


def _generated_feature_modeling_helpers() -> str:
    return dedent(
        """
        import warnings
        from collections.abc import Mapping

        import numpy as np
        import scanpy as sc
        from anndata import AnnData
        from scipy import sparse


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


        def _openbio_feature_history(adata):
            metadata = adata.uns.get("openbio_singlecell")
            if not isinstance(metadata, Mapping):
                return []
            history = metadata.get("analysis_history")
            if not isinstance(history, Mapping):
                return []
            return [entry for entry in history.values() if isinstance(entry, Mapping)]


        def _openbio_feature_state(adata, source_kind, layer_name=None):
            x_state = "logged_unverified" if isinstance(adata.uns.get("log1p"), Mapping) else "unknown"
            x_evidence = (
                "AnnData uns['log1p'] marker without normalization provenance"
                if x_state != "unknown"
                else None
            )
            layer_states = {}
            for entry in _openbio_feature_history(adata):
                operation = entry.get("operation")
                parameters = entry.get("parameters")
                parameters = parameters if isinstance(parameters, Mapping) else {}
                if operation == "snapshot_expression":
                    layer_states["counts"] = ("counts", "OpenBio Snapshot Expression history")
                    if parameters.get("source") == "X":
                        x_state, x_evidence = "counts", "OpenBio Snapshot Expression history"
                elif operation == "normalize_to_layer":
                    destination = parameters.get("output_layer")
                    transform = parameters.get("transform")
                    if isinstance(destination, str):
                        state = {"log1p": "logged", "none": "normalized", "sqrt": "transformed"}.get(
                            transform, "derived"
                        )
                        layer_states[destination] = (state, "OpenBio Normalize to Layer history")
                elif operation == "pearson_residuals_to_layer":
                    destination = parameters.get("output_layer")
                    if isinstance(destination, str):
                        layer_states[destination] = (
                            "pearson_residuals",
                            "OpenBio Pearson Residuals history",
                        )
                elif operation == "scale_to_layer":
                    destination = parameters.get("output_layer")
                    if isinstance(destination, str):
                        layer_states[destination] = ("scaled", "OpenBio Scale history")

                if operation == "normalize_total":
                    x_state, x_evidence = "normalized", "OpenBio Normalize Total history"
                elif operation == "log1p":
                    if x_state == "normalized":
                        x_state, x_evidence = (
                            "logged",
                            "OpenBio Normalize Total followed by Log1p history",
                        )
                    else:
                        x_state, x_evidence = (
                            "logged_unverified",
                            "OpenBio Log1p history without proven prior normalization",
                        )
            if source_kind == "layer":
                return layer_states.get(layer_name, ("unknown", None))
            return x_state, x_evidence


        def _openbio_require_feature_input(adata, operation, minimum_features=1, minimum_cells=2):
            if adata.isbacked:
                raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
            if adata.n_obs < minimum_cells or adata.n_vars < minimum_features:
                raise ValueError(
                    f"{operation} requires at least {minimum_cells} cell(s) and {minimum_features} feature(s)."
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


        def _openbio_validate_strict_counts(
            adata,
            matrix,
            values,
            source_kind,
            source_layer,
            operation,
            reject_zero_genes,
        ):
            state, evidence = _openbio_feature_state(adata, source_kind, source_layer)
            if values.size and (values < 0).any():
                raise ValueError(
                    f"{operation} source contains negative values outside the non-negative count-model domain."
                )
            if values.size == 0 or not (values > 0).any():
                raise ValueError(f"{operation} source contains no positive expression values.")
            if not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-8):
                warnings.warn(
                    f"{operation} source contains fractional non-negative values. The selected count-based "
                    "method is executable, but the input cannot be described as verified integer UMI counts.",
                    UserWarning,
                    stacklevel=2,
                )
            if state not in {"counts", "unknown"}:
                detail = f"provenance describes the source as {state} expression ({evidence})"
                warnings.warn(
                    f"{detail}; honoring the expert-selected source without claiming verified raw UMI counts.",
                    UserWarning,
                    stacklevel=2,
                )
            cell_totals = np.asarray(matrix.sum(axis=1), dtype=float).ravel()
            gene_totals = np.asarray(matrix.sum(axis=0), dtype=float).ravel()
            zero_cells = int((cell_totals <= 0).sum())
            if zero_cells:
                raise ValueError(
                    f"{operation} source contains {zero_cells} cells with zero total counts; filter them first."
                )
            zero_genes = int((gene_totals <= 0).sum())
            if reject_zero_genes and zero_genes:
                raise ValueError(
                    f"{operation} source contains {zero_genes} genes with zero total counts; filter them first."
                )


        def _openbio_validate_logged(
            adata,
            matrix,
            values,
            source_kind,
            source_layer,
            operation,
        ):
            if values.size and (values < 0).any():
                raise ValueError(f"{operation} normalized log source contains negative expression values.")
            if values.size == 0 or not (values > 0).any():
                raise ValueError(f"{operation} normalized log source contains no positive expression values.")
            state, evidence = _openbio_feature_state(adata, source_kind, source_layer)
            if state == "logged_unverified":
                warnings.warn(
                    "The selected source has a Scanpy log1p marker, but OpenBio provenance cannot prove that "
                    "total-count normalization preceded logarithmization; confirm the external data contract.",
                    UserWarning,
                    stacklevel=2,
                )
            elif state != "logged":
                detail = (
                    f"provenance describes the source as {state} expression ({evidence})"
                    if state != "unknown"
                    else "OpenBio provenance cannot prove the source is normalized log1p expression"
                )
                warnings.warn(
                    f"{detail}; honoring the expert-selected source for the {operation} dispersion flavor.",
                    UserWarning,
                    stacklevel=2,
                )
            cell_totals = np.asarray(matrix.sum(axis=1), dtype=float).ravel()
            zero_cells = int((cell_totals <= 0).sum())
            if zero_cells:
                raise ValueError(f"{operation} source contains {zero_cells} cells with zero total expression.")


        def _openbio_validate_destination(
            adata,
            source_kind,
            source_layer,
            output_layer,
            overwrite_existing,
            operation,
        ):
            if not isinstance(output_layer, str) or not output_layer.strip():
                raise ValueError(f"{operation} output_layer cannot be empty.")
            output_layer = output_layer.strip()
            if output_layer == "counts":
                raise ValueError(f"{operation} cannot overwrite the reserved canonical counts layer.")
            if source_kind == "layer" and output_layer == source_layer:
                raise ValueError(f"{operation} output_layer must differ from the selected source layer.")
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
            small = [
                (value, int(count))
                for value, count in counts.items()
                if int(count) < minimum_group_size
            ]
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


        def _openbio_validate_raw_snapshot(adata, operation):
            raw = adata.raw
            if raw is None:
                raise ValueError(f"{operation} subset=True requires a full-gene Raw snapshot of post-QC counts.")
            if raw.n_obs != adata.n_obs or not raw.obs_names.equals(adata.obs_names):
                raise ValueError(f"{operation} Raw snapshot observations are not aligned to the current AnnData.")
            if not raw.var_names.is_unique:
                raise ValueError(f"{operation} Raw snapshot requires unique feature identifiers.")
            if raw.n_vars < adata.n_vars or not adata.var_names.isin(raw.var_names).all():
                raise ValueError(f"{operation} Raw snapshot does not preserve the full current feature set.")
            values = np.asarray(raw.X.data if sparse.issparse(raw.X) else np.asarray(raw.X).ravel())
            if values.size == 0 or (values < 0).any() or not (values > 0).any():
                raise ValueError(f"{operation} Raw snapshot is not a non-negative count matrix.")
            if not np.isfinite(values).all() or not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-8):
                raise ValueError(f"{operation} Raw snapshot is not an integer count matrix.")
            totals = np.asarray(raw.X.sum(axis=1), dtype=float).ravel()
            if (totals <= 0).any():
                raise ValueError(f"{operation} Raw snapshot contains cells with zero total counts.")
        """
    ).strip()


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
        _generated_feature_modeling_helpers()
        + "\n\n\n"
        + dedent(
            f"""
            def pearson_residuals_to_layer(adata):
                operation = "Pearson Residuals"
                _openbio_require_feature_input(adata, operation, minimum_features=1, minimum_cells=1)
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
                _openbio_validate_strict_counts(
                    adata,
                    matrix,
                    values,
                    source_kind,
                    source_layer,
                    operation,
                    reject_zero_genes=True,
                )
                output_layer = _openbio_validate_destination(
                    adata,
                    source_kind,
                    source_layer,
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
                output = adata.copy()
                output.layers[output_layer] = residuals.copy()
                return output
            """
        ).strip()
    )


def _highly_variable_genes_code(
    expression: ExpressionSource,
    *,
    n_top_genes: int,
    flavor: str,
    batch_key: str,
    requested_keep: list[str],
    subset: bool,
    theta: float,
    clip: float | None,
    chunksize: int,
    span: float,
    n_bins: int,
    overwrite_existing: bool,
) -> str:
    source_layer = expression.layer_name if expression.kind == "layer" else None
    return (
        _generated_feature_modeling_helpers()
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
                _openbio_require_feature_input(
                    adata,
                    operation,
                    minimum_features=1,
                    minimum_cells=minimum_cells,
                )
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
                    _openbio_validate_logged(
                        adata, matrix, values, source_kind, source_layer, operation
                    )
                else:
                    _openbio_validate_strict_counts(
                        adata,
                        matrix,
                        values,
                        source_kind,
                        source_layer,
                        operation,
                        reject_zero_genes=False,
                    )
                gene_totals = np.asarray(matrix.sum(axis=0), dtype=float).ravel()
                expressed_features = int((gene_totals > 0).sum())
                if expressed_features < n_top_genes:
                    warnings.warn(
                        f"{{operation}} requested {{n_top_genes}} genes but only "
                        f"{{expressed_features}} features have positive totals; the algorithm may return fewer.",
                        UserWarning,
                        stacklevel=2,
                    )
                minimum_group_size = 1 if flavor == "pearson_residuals" else 2
                _openbio_validate_batch(adata, batch_key, operation, minimum_group_size)
                if subset:
                    try:
                        _openbio_validate_raw_snapshot(adata, operation)
                    except (TypeError, ValueError) as error:
                        warnings.warn(
                            f"{{error}} Subsetting will proceed, but no verified full-gene Raw snapshot is retained.",
                            UserWarning,
                            stacklevel=2,
                        )
                stale = [column for column in _OPENBIO_HVG_COLUMNS if column in adata.var.columns]
                if (stale or "hvg" in adata.uns) and not overwrite_existing:
                    raise ValueError(
                        "Highly Variable Genes annotations already exist; enable overwrite_existing to replace them."
                    )
                output = adata.copy()
                output.var.drop(
                    columns=[column for column in _OPENBIO_HVG_COLUMNS if column in output.var],
                    inplace=True,
                )
                output.uns.pop("hvg", None)
                layer = source_layer if source_kind == "layer" else None
                if flavor == "pearson_residuals":
                    sc.experimental.pp.highly_variable_genes(
                        output,
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
                elif flavor in {"seurat_v3", "seurat_v3_paper"}:
                    sc.pp.highly_variable_genes(
                        output,
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
                        output,
                        layer=layer,
                        n_top_genes=n_top_genes,
                        n_bins=n_bins,
                        flavor=flavor,
                        subset=False,
                        inplace=True,
                        batch_key=batch_key or None,
                        check_values=True,
                    )
                algorithm = output.var["highly_variable"].astype(bool).to_numpy(copy=True)
                forced = output.var_names.isin(requested_keep)
                output.var["highly_variable_algorithm"] = algorithm
                output.var["highly_variable_forced"] = forced
                output.var["highly_variable"] = algorithm | forced
                if subset:
                    output = output[:, output.var["highly_variable"].astype(bool)].copy()
                return output
            """
        ).strip()
    )


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
        _generated_feature_modeling_helpers()
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
                    source_kind,
                    source_layer,
                    output_layer,
                    overwrite_existing,
                    operation,
                )
                if sparse.issparse(matrix) and zero_center:
                    _openbio_validate_dense_budget(adata.shape, max_dense_gib, operation)
                output = adata.copy()
                work = AnnData(X=matrix.copy())
                sc.pp.scale(work, zero_center=zero_center, max_value=max_value, copy=False)
                if sparse.issparse(matrix) and not zero_center and not sparse.issparse(work.X):
                    raise RuntimeError("Scale unexpectedly densified sparse input with zero_center=False.")
                output.layers[output_layer] = work.X.copy()
                return output
            """
        ).strip()
    )


class OpenBioSingleCellPearsonResidualsToLayer(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Pearson residual source",
        default="layer",
        layer_input_id="source_layer",
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPearsonResidualsToLayer",
            display_name="Pearson Residuals to Layer",
            category=CATEGORY,
            description="Compute analytic Pearson residuals and store them in an AnnData layer.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Float.Input("theta", default=100.0, step=10.0, advanced=True),
                io.Combo.Input(
                    "clipping_mode",
                    options=["sqrt_n_obs", "custom", "none"],
                    default="sqrt_n_obs",
                ),
                io.Float.Input("custom_clip", default=10.0, min=0.0, step=1.0, advanced=True),
                io.String.Input("output_layer", default="analytic_pearson_residuals", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_dense_gib", default=2.0, step=0.25, advanced=True),
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
        source: DynamicExpressionSource | None = None,
        theta: float = 100.0,
        clipping_mode: str = "sqrt_n_obs",
        custom_clip: float = 10.0,
        output_layer: str = "analytic_pearson_residuals",
        overwrite_existing: bool = False,
        max_dense_gib: float = 2.0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        operation = "Pearson Residuals"
        _require_in_memory_nonempty(adata, operation=operation)
        _require_unique_axes(adata, operation=operation)
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        matrix, cell_totals, gene_totals, report_warnings, state, evidence = _validate_strict_count_source(
            adata,
            expression,
            operation=operation,
            reject_zero_genes=True,
        )
        theta = _positive_float(theta, name=f"{operation} theta", allow_infinity=True)
        clipping_value = custom_clip
        clipping_mode, scanpy_clip, resolved_clip = _resolve_clipping(
            clipping_mode,
            clipping_value,
            n_obs=int(adata.n_obs),
            value_name=f"{operation} custom_clip",
        )
        output_layer, replaced_existing = _validate_destination_layer(
            adata,
            expression,
            output_layer,
            overwrite_existing,
            operation=operation,
        )
        estimated_dense_gib, max_dense_gib = _validate_dense_budget(
            tuple(adata.shape),
            max_dense_gib,
            operation=operation,
        )
        input_sparse = bool(science.sparse.issparse(matrix))
        if input_sparse:
            warning = (
                "Analytic Pearson residuals convert sparse counts to a dense float64 matrix; the reported "
                "allocation is a lower bound because working arrays coexist."
            )
            report_warnings.append(warning)
            python_warnings.warn(warning, UserWarning, stacklevel=2)
        if adata.n_obs == 1:
            report_warnings.append(
                "Pearson Residuals received one cell. Scanpy defines the residual transform, but per-gene "
                "sample variance with ddof=1 is unavailable and is reported as null."
            )
        if adata.n_vars == 1:
            report_warnings.append(
                "Pearson Residuals received one feature; the transform is finite but offers no multigene "
                "representation for downstream dimensionality reduction."
            )

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
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
        if not bool(science.np.isfinite(residuals).all()):
            raise ValueError(
                "Pearson Residuals produced non-finite values; verify that zero-total cells and genes were removed."
            )
        output = adata.copy()
        output.layers[output_layer] = residuals.copy()
        residual_variance = (
            residuals.var(axis=0, ddof=1)
            if cells > 1
            else science.np.full(genes, science.np.nan, dtype=float)
        )
        reported_theta: float | str = "infinity" if theta == float("inf") else theta
        if resolved_clip == "infinity":
            at_clip = 0
        else:
            at_clip = int(science.np.isclose(science.np.abs(residuals), float(resolved_clip)).sum())
        parameters = {
            **expression.parameters(),
            "theta": reported_theta,
            "clipping_mode": clipping_mode,
            "resolved_clip": resolved_clip,
            "check_values": False,
            "output_layer": output_layer,
            "overwrite_existing": overwrite_existing,
            "max_dense_gib": max_dense_gib,
        }
        finish_adata(
            output,
            "pearson_residuals_to_layer",
            parameters,
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
                f"Computed analytic Pearson residuals from expert-selected non-negative source {source_label!r} with "
                f"scanpy.experimental.pp.normalize_pearson_residuals using a shared negative-binomial theta "
                f"of {theta!r} and clipping mode {clipping_mode!r} (resolved bound {resolved_clip!r}). The dense "
                f"residual matrix was stored in layer {output_layer!r}; X, the source, and Raw were preserved."
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
                "source_state": state,
                "source_state_evidence": evidence,
                "input_matrix": _matrix_description(matrix),
                "input_integer_like": _matrix_is_integer_like(matrix),
                "input_cell_totals": summarize_numeric(cell_totals),
                "input_gene_totals": summarize_numeric(gene_totals),
                "zero_total_cells": 0,
                "zero_total_genes": 0,
                "theta": reported_theta,
                "clipping_mode": clipping_mode,
                "resolved_clip": resolved_clip,
                "output_layer": output_layer,
                "output_matrix": _matrix_description(residuals),
                "residual_values": summarize_numeric(residuals),
                "residual_variances_ddof1": summarize_numeric(residual_variance),
                "entries_at_clipping_bound": at_clip,
                "estimated_dense_array_gib": estimated_dense_gib,
                "dense_budget_gib": max_dense_gib,
                "sparse_to_dense": input_sparse,
                "replaced_existing_layer": replaced_existing,
                "x_preserved": True,
                "source_preserved": True,
                "raw_preserved": True,
                "raw_present": adata.raw is not None,
            },
            parameters=parameters,
            references=(LAUSE_REFERENCE, SCANPY_REFERENCE, ANNDATA_REFERENCE),
            software_packages=PREPROCESS_SOFTWARE_PACKAGES,
            warnings=report_warnings,
            limitations=(
                "Analytic Pearson residuals are a model-based expression representation, not counts or "
                "sample-level inference.",
                "The Scanpy API is experimental; the dense allocation estimate covers one float64 result array "
                "and not all simultaneous intermediates.",
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
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellHighlyVariableGenes(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Highly Variable Genes source",
        default="layer",
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellHighlyVariableGenes",
            display_name="Highly Variable Genes",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.Int.Input("n_top_genes", default=2000, min=1),
                io.Combo.Input(
                    "flavor",
                    options=["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper", "pearson_residuals"],
                    default="seurat",
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("batch_key", default=""),
                io.String.Input("always_keep_genes", default="", advanced=True),
                io.Boolean.Input("subset", default=False, advanced=True),
                io.Float.Input("theta", default=100.0, step=10.0, advanced=True),
                io.Combo.Input(
                    "clipping_mode",
                    options=["sqrt_n_obs", "custom", "none"],
                    default="sqrt_n_obs",
                    advanced=True,
                ),
                io.Float.Input("custom_clip", default=10.0, min=0.0, step=1.0, advanced=True),
                io.Int.Input("chunksize", default=1000, min=1, advanced=True),
                io.Float.Input("span", default=0.3, max=1.0, step=0.05, advanced=True),
                io.Int.Input("n_bins", default=20, min=1, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
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
        n_top_genes: int = 2000,
        flavor: str = "seurat",
        source: DynamicExpressionSource | None = None,
        batch_key: str = "",
        always_keep_genes: str = "",
        subset: bool = False,
        theta: float = 100.0,
        clipping_mode: str = "sqrt_n_obs",
        custom_clip: float = 10.0,
        chunksize: int = 1000,
        span: float = 0.3,
        n_bins: int = 20,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        operation = "Highly Variable Genes"
        _require_in_memory_nonempty(adata, operation=operation)
        _require_unique_axes(adata, operation=operation)
        n_top_genes = _positive_integer(n_top_genes, name=f"{operation} n_top_genes")
        if not isinstance(flavor, str):
            raise TypeError(f"{operation} flavor must be a string.")
        flavor = flavor.strip()
        supported_flavors = LOG_HVG_FLAVORS | COUNT_HVG_FLAVORS
        if flavor not in supported_flavors:
            raise ValueError(f"Unsupported {operation} flavor: {flavor!r}")
        if flavor != "pearson_residuals" and adata.n_obs < 2:
            raise ValueError(f"{operation} flavor {flavor!r} backend requires at least two cells.")
        if not isinstance(batch_key, str):
            raise TypeError(f"{operation} batch_key must be a string.")
        batch_key = batch_key.strip()
        if not isinstance(always_keep_genes, str):
            raise TypeError(f"{operation} always_keep_genes must be a comma-separated string.")
        if not isinstance(subset, bool):
            raise TypeError(f"{operation} subset must be a boolean.")
        if not isinstance(overwrite_existing, bool):
            raise TypeError(f"{operation} overwrite_existing must be a boolean.")

        report_warnings: list[str] = []
        if n_top_genes > adata.n_vars:
            report_warnings.append(
                f"{operation} requested {n_top_genes} genes from {int(adata.n_vars)} features; Scanpy will "
                "return at most the eligible feature count."
            )
        if adata.n_obs == 1:
            report_warnings.append(
                f"{operation} flavor {flavor!r} received one cell. The selected backend can execute, but "
                "population-level variability is not independently supported."
            )

        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        source_label = _expression_source_label(expression)
        if flavor in LOG_HVG_FLAVORS:
            matrix, cell_totals, state_warnings, state, evidence = _validate_logged_source(
                adata,
                expression,
                operation=operation,
            )
            report_warnings.extend(state_warnings)
            _, gene_nonzero = matrix_totals_and_nonzero(matrix, axis=0)
            gene_totals, _ = matrix_totals_and_nonzero(matrix, axis=0)
            count_validation = False
        else:
            matrix, cell_totals, gene_totals, count_warnings, state, evidence = _validate_strict_count_source(
                adata,
                expression,
                operation=operation,
                reject_zero_genes=False,
            )
            report_warnings.extend(count_warnings)
            _, gene_nonzero = matrix_totals_and_nonzero(matrix, axis=0)
            count_validation = True
        gene_totals = science.np.asarray(gene_totals, dtype=float)
        gene_nonzero = science.np.asarray(gene_nonzero, dtype=int)
        expressed_features = int((gene_totals > 0).sum())
        if expressed_features < n_top_genes:
            report_warnings.append(
                f"{operation} requested {n_top_genes} genes but only {expressed_features} features have positive "
                "totals; the algorithm may return fewer."
            )

        group_sizes, batch_warnings = _batch_group_sizes(
            adata,
            batch_key,
            operation=operation,
            minimum_group_size=1 if flavor == "pearson_residuals" else 2,
        )
        report_warnings.extend(batch_warnings)
        requested_keep = list(dict.fromkeys(gene.strip() for gene in always_keep_genes.split(",") if gene.strip()))
        input_var_names = adata.var_names
        present_keep = [gene for gene in requested_keep if gene in input_var_names]
        missing_keep = [gene for gene in requested_keep if gene not in input_var_names]
        if missing_keep:
            report_warnings.append(f"Requested keep genes were not found: {missing_keep}")

        raw_snapshot = {
            "present": adata.raw is not None,
            "cells": int(adata.raw.n_obs) if adata.raw is not None else None,
            "features": int(adata.raw.n_vars) if adata.raw is not None else None,
            "validated_full_gene_snapshot": False,
        }
        if subset:
            try:
                raw_snapshot = {
                    **_validate_full_gene_raw_snapshot(adata, operation=operation),
                    "validated_full_gene_snapshot": True,
                }
            except (TypeError, ValueError) as error:
                raw_snapshot["validation_issue"] = str(error)
                report_warnings.append(
                    f"{error} Subsetting will proceed, but no verified full-gene Raw snapshot is retained."
                )
        prior_columns = [column for column in HVG_RESULT_COLUMNS if column in adata.var.columns]
        prior_uns = "hvg" in adata.uns
        if (prior_columns or prior_uns) and not overwrite_existing:
            raise ValueError(
                f"{operation} annotations already exist ({prior_columns!r}); enable overwrite_existing to replace them."
            )

        resolved_theta: float | None = None
        scanpy_clip: float | None = None
        resolved_clip: float | str | None = None
        resolved_chunksize: int | None = None
        resolved_span: float | None = None
        resolved_n_bins: int | None = None
        if flavor == "pearson_residuals":
            resolved_theta = _positive_float(theta, name=f"{operation} theta", allow_infinity=True)
            clipping_mode, scanpy_clip, resolved_clip = _resolve_clipping(
                clipping_mode,
                custom_clip,
                n_obs=int(adata.n_obs),
                value_name=f"{operation} custom_clip",
            )
            if clipping_mode == "sqrt_n_obs":
                scanpy_clip = float(resolved_clip)
            resolved_chunksize = _positive_integer(chunksize, name=f"{operation} chunksize")
        elif flavor in {"seurat_v3", "seurat_v3_paper"}:
            resolved_span = _positive_float(span, name=f"{operation} span")
            if resolved_span > 1:
                raise ValueError(f"{operation} span must be no greater than one.")
        else:
            resolved_n_bins = _positive_integer(n_bins, name=f"{operation} n_bins")

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        output.var.drop(columns=[column for column in HVG_RESULT_COLUMNS if column in output.var], inplace=True)
        output.uns.pop("hvg", None)
        try:
            if flavor == "pearson_residuals":
                science.sc.experimental.pp.highly_variable_genes(
                    output,
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
                    output,
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
                    output,
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
                    "Highly Variable Genes seurat_v3 flavors require scikit-misc; install the plugin's "
                    "recommended scientific dependencies."
                ) from error
            raise

        algorithm_mask = output.var["highly_variable"].astype(bool).to_numpy(copy=True)
        forced_mask = output.var_names.isin(present_keep)
        final_mask = algorithm_mask | forced_mask
        output.var["highly_variable_algorithm"] = algorithm_mask
        output.var["highly_variable_forced"] = forced_mask
        output.var["highly_variable"] = final_mask
        algorithm_count = int(algorithm_mask.sum())
        forced_present_count = int(forced_mask.sum())
        forced_already_selected = int((algorithm_mask & forced_mask).sum())
        forced_added = int((~algorithm_mask & forced_mask).sum())
        final_count = int(final_mask.sum())
        metric_columns = [column for column in HVG_RESULT_COLUMNS if column in output.var]
        metric_summaries = {
            column: summarize_numeric(output.var[column].to_numpy())
            for column in metric_columns
            if column not in {"highly_variable", "highly_variable_algorithm", "highly_variable_forced"}
            and science.np.issubdtype(output.var[column].dtype, science.np.number)
        }
        if "highly_variable_rank" in output.var:
            rank = output.var.loc[algorithm_mask, "highly_variable_rank"]
            top_features = [str(name) for name in rank.dropna().sort_values().index[: min(20, algorithm_count)]]
        else:
            selected_names = output.var_names[algorithm_mask]
            top_features = [str(name) for name in selected_names[:20]]
        nbatches_distribution: dict[str, int] = {}
        if "highly_variable_nbatches" in output.var:
            nbatches_distribution = {
                str(key): int(value)
                for key, value in output.var["highly_variable_nbatches"].value_counts(dropna=False).sort_index().items()
            }
        intersection_count: int | None = None
        intersection_rate: float | None = None
        if "highly_variable_intersection" in output.var:
            intersection_count = int(output.var["highly_variable_intersection"].astype(bool).sum())
            intersection_rate = intersection_count / genes
        if subset:
            output = output[:, output.var["highly_variable"].astype(bool)].copy()
        parameters = {
            "n_top_genes": n_top_genes,
            "flavor": flavor,
            **expression.parameters(),
            "batch_key": batch_key,
            "always_keep_genes": requested_keep,
            "subset": subset,
            "theta": "infinity" if resolved_theta == float("inf") else resolved_theta,
            "clipping_mode": clipping_mode if flavor == "pearson_residuals" else None,
            "resolved_clip": resolved_clip,
            "chunksize": resolved_chunksize,
            "span": resolved_span,
            "n_bins": resolved_n_bins,
            "check_values": flavor in LOG_HVG_FLAVORS,
            "overwrite_existing": overwrite_existing,
        }
        finish_adata(
            output,
            "highly_variable_genes",
            parameters,
            cells,
            genes,
            started_at,
            warnings=report_warnings,
        )
        method_name = (
            "scanpy.experimental.pp.highly_variable_genes"
            if flavor == "pearson_residuals"
            else "scanpy.pp.highly_variable_genes"
        )
        expected_state = (
            "documented normalized log1p input"
            if flavor in LOG_HVG_FLAVORS
            else "documented non-negative count-like input"
        )
        grouping_text = (
            f"within each observed group in obs[{batch_key!r}] before Scanpy's cross-group rank merge"
            if batch_key
            else "without a grouping key"
        )
        reference = {
            "seurat": SATIJA_REFERENCE,
            "cell_ranger": ZHENG_REFERENCE,
            "seurat_v3": STUART_REFERENCE,
            "seurat_v3_paper": STUART_REFERENCE,
            "pearson_residuals": LAUSE_REFERENCE,
        }[flavor]
        software_packages = (
            (*PREPROCESS_SOFTWARE_PACKAGES, "scikit-misc")
            if flavor.startswith("seurat_v3")
            else PREPROCESS_SOFTWARE_PACKAGES
        )
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellHighlyVariableGenes",
            title="Highly variable gene selection summary",
            operation="highly_variable_genes",
            methods=(
                f"Selected up to {n_top_genes:,} algorithmic highly variable genes from {expected_state} in "
                f"source {source_label!r} using {method_name} flavor {flavor!r}, {grouping_text}. Forced genes "
                "were applied only after the algorithmic mask and were recorded separately."
            ),
            results=(
                f"The {flavor!r} algorithm selected {algorithm_count:,} genes; {forced_added:,} additional "
                f"requested genes were forced into the final set of {final_count:,}. The returned AnnData has "
                f"{int(output.n_vars):,} current features (subset={subset})."
            ),
            key_results={
                "input_cells": cells,
                "input_features": genes,
                "output_cells": int(output.n_obs),
                "output_features": int(output.n_vars),
                "source": source_label,
                "source_state": state,
                "source_state_evidence": evidence,
                "source_matrix": _matrix_description(matrix),
                "strict_count_validation": False,
                "count_model_validation": count_validation,
                "input_cell_totals": summarize_numeric(cell_totals),
                "input_gene_totals": summarize_numeric(gene_totals),
                "expressed_features": expressed_features,
                "cells_per_expressed_feature": summarize_numeric(gene_nonzero[gene_totals > 0]),
                "requested_algorithm_count": n_top_genes,
                "algorithm_selected_count": algorithm_count,
                "forced_requested_count": len(requested_keep),
                "forced_present_count": forced_present_count,
                "forced_already_algorithm_selected": forced_already_selected,
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
                "raw_snapshot": raw_snapshot,
                "raw_preserved": True,
                "subset_applied": subset,
            },
            parameters=parameters,
            references=(reference, SCANPY_REFERENCE, ANNDATA_REFERENCE),
            software_packages=software_packages,
            warnings=report_warnings,
            limitations=(
                "Highly variable genes are features with high model- or dispersion-based variation; they are not "
                "marker genes or evidence of differential expression.",
                "A grouping key changes recurrence/ranking across groups but does not perform batch integration "
                "or provide sample-level inference.",
                "Forced genes are user overrides and are explicitly excluded from the algorithm-selected count.",
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
                theta=resolved_theta if resolved_theta is not None else 100.0,
                clip=scanpy_clip,
                chunksize=resolved_chunksize if resolved_chunksize is not None else 1000,
                span=resolved_span if resolved_span is not None else 0.3,
                n_bins=resolved_n_bins if resolved_n_bins is not None else 20,
                overwrite_existing=overwrite_existing,
            ),
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellScale(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Scale source",
        default="layer",
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellScale",
            display_name="Scale",
            category=CATEGORY,
            description="Scale one expression source into a named layer without replacing X or Raw.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Boolean.Input("zero_center", default=True),
                io.Combo.Input("clipping_mode", options=["custom", "none"], default="custom"),
                io.Float.Input("custom_max_value", default=10.0, min=0.0, step=1.0, advanced=True),
                io.String.Input("output_layer", default="scaled", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
                io.Float.Input("max_dense_gib", default=2.0, step=0.25, advanced=True),
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
        source: DynamicExpressionSource | None = None,
        zero_center: bool = True,
        clipping_mode: str = "custom",
        custom_max_value: float = 10.0,
        output_layer: str = "scaled",
        overwrite_existing: bool = False,
        max_dense_gib: float = 2.0,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        operation = "Scale"
        _require_in_memory_nonempty(adata, operation=operation)
        if adata.n_obs < 2:
            raise ValueError(f"{operation} requires at least two cells because Scanpy uses sample variance (ddof=1).")
        _require_unique_axes(adata, operation=operation)
        if not isinstance(zero_center, bool):
            raise TypeError(f"{operation} zero_center must be a boolean.")
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        matrix = expression.matrix(adata)
        if tuple(matrix.shape) != tuple(adata.shape):
            raise ValueError(f"{operation} source must align to the current AnnData shape.")
        _validate_finite_numeric_matrix(matrix, source_label=f"{operation} source")
        output_layer, replaced_existing = _validate_destination_layer(
            adata,
            expression,
            output_layer,
            overwrite_existing,
            operation=operation,
        )
        if not isinstance(clipping_mode, str):
            raise TypeError(f"{operation} clipping_mode must be a string.")
        clipping_mode = clipping_mode.strip()
        if clipping_mode not in {"custom", "none"}:
            raise ValueError(f"Unsupported {operation} clipping_mode: {clipping_mode!r}")
        resolved_max_value = (
            _nonnegative_float(custom_max_value, name=f"{operation} custom_max_value")
            if clipping_mode == "custom"
            else None
        )
        max_dense_gib = _positive_float(max_dense_gib, name=f"{operation} max_dense_gib")
        estimated_dense_gib = _dense_array_gib(tuple(adata.shape))
        input_sparse = bool(science.sparse.issparse(matrix))
        report_warnings: list[str] = []
        if input_sparse and zero_center:
            _validate_dense_budget(tuple(adata.shape), max_dense_gib, operation=operation)
            warning = (
                "Zero-centering a sparse expression matrix densifies it; the reported float64 allocation is a "
                "lower bound because working arrays coexist."
            )
            report_warnings.append(warning)
            python_warnings.warn(warning, UserWarning, stacklevel=2)
        if input_sparse and not zero_center and resolved_max_value is not None:
            report_warnings.append(
                "With zero_center=False, Scanpy clips only the positive upper tail; negative values are not "
                "symmetrically truncated."
            )
        state, evidence = _feature_expression_state(adata, expression)
        if state in {"counts", "unknown"}:
            report_warnings.append(
                f"The selected source state is {state!r}; gene scaling is generally applied to a normalized "
                "continuous expression representation."
            )
        feature_mean, feature_std, constant_mask = _scale_feature_statistics(matrix)

        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata.copy()
        work = science.ad.AnnData(X=matrix.copy())
        science.sc.pp.scale(
            work,
            zero_center=zero_center,
            max_value=resolved_max_value,
            copy=False,
        )
        if input_sparse and not zero_center and not science.sparse.issparse(work.X):
            raise RuntimeError("Scale unexpectedly densified sparse input with zero_center=False.")
        scaled = work.X.copy()
        output.layers[output_layer] = scaled
        if resolved_max_value is None:
            boundary_count = 0
        else:
            scaled_values = _matrix_values(scaled)
            upper = science.np.isclose(scaled_values, resolved_max_value)
            if zero_center:
                lower = science.np.isclose(scaled_values, -resolved_max_value)
                boundary_count = int((upper | lower).sum())
            else:
                boundary_count = int(upper.sum())
        parameters = {
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
            output,
            "scale_to_layer",
            parameters,
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
                f"Applied {scaling_method} to source {source_label!r} with scanpy.pp.scale using per-gene "
                f"sample standard deviations (ddof=1); clipping was {clipping_text}. The result was written to "
                f"layer {output_layer!r}, preserving X, the selected source, other layers, and Raw."
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
                "source_state": state,
                "source_state_evidence": evidence,
                "input_matrix": _matrix_description(matrix),
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
                "source_preserved": True,
                "raw_preserved": True,
                "raw_present": adata.raw is not None,
            },
            parameters=parameters,
            references=(SCANPY_REFERENCE, ANNDATA_REFERENCE),
            software_packages=PREPROCESS_SOFTWARE_PACKAGES,
            warnings=report_warnings,
            limitations=(
                "Scaling changes feature units but does not normalize library size, integrate batches, or provide "
                "statistical inference.",
                "Constant-gene standard deviations are replaced with one internally by Scanpy; such genes are "
                "retained and disclosed rather than removed.",
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
        return io.NodeOutput(output, report, code)


PREPROCESS_NODE_CLASSES = [
    OpenBioSingleCellNormalizeTotal,
    OpenBioSingleCellLog1p,
    OpenBioSingleCellNormalizeToLayer,
    OpenBioSingleCellPearsonResidualsToLayer,
    OpenBioSingleCellHighlyVariableGenes,
    OpenBioSingleCellScale,
]
