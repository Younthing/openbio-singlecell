from __future__ import annotations

import hashlib
import inspect
import json
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from .contracts import TableResult

if TYPE_CHECKING:
    from pandas import DataFrame


MARKER_EVIDENCE_SCHEMA_VERSION = 2
MARKER_GENES_NODE_ID = "OpenBioSingleCellMarkerGenes"
FILTER_MARKER_GENES_NODE_ID = "OpenBioSingleCellFilterMarkerGenes"
MARKER_COLUMNS = [
    "group",
    "gene",
    "rank",
    "score",
    "log2_fold_change_approx",
    "p_value",
    "p_adjusted",
    "fraction_in_group",
    "fraction_reference",
]
MARKER_UNIVERSE_COLUMNS = ["gene", "universe_rank"]
MARKER_METHODS = ("wilcoxon", "t-test", "t-test_overestim_var")
_MARKER_NUMERIC_COLUMNS = MARKER_COLUMNS[2:]
_MARKER_UNIT_INTERVAL_COLUMNS = (
    "p_value",
    "p_adjusted",
    "fraction_in_group",
    "fraction_reference",
)
_SHA256_LENGTH = 64


def _canonical_sha256(payload: Any) -> str:
    encoder = json.JSONEncoder(
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    digest = hashlib.sha256()
    for chunk in encoder.iterencode(payload):
        digest.update(chunk.encode("utf-8"))
    return digest.hexdigest()


def _validate_identifier_values(values: Any, *, name: str) -> list[str]:
    rendered = []
    for value in values:
        if not isinstance(value, str):
            raise TypeError(f"{name} identifiers must be strings; received {type(value).__name__}.")
        if not value.strip():
            raise ValueError(f"{name} identifiers cannot be empty or whitespace-only.")
        if value != value.strip():
            raise ValueError(f"{name} identifiers cannot contain surrounding whitespace.")
        rendered.append(value)
    if len(set(rendered)) != len(rendered):
        raise ValueError(f"{name} identifiers must be unique.")
    return rendered


def _validate_integer(value: Any, *, name: str, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer.")
    minimum = 0 if allow_zero else 1
    if value < minimum:
        comparison = "non-negative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {comparison}.")
    return value


def _validate_positive_float(value: Any, *, name: str, np: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a finite number greater than zero.")
    normalized = float(value)
    if not bool(np.isfinite(normalized)) or normalized <= 0.0:
        raise ValueError(f"{name} must be finite and greater than zero.")
    return normalized


def _validate_finite_threshold(value: Any, *, name: str, np: Any, unit_interval: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a finite number.")
    normalized = float(value)
    if not bool(np.isfinite(normalized)):
        raise ValueError(f"{name} must be finite.")
    if unit_interval and not 0.0 <= normalized <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1 inclusive.")
    return normalized


def _validate_canonical_marker_table(frame: Any, *, np: Any, pd: Any) -> None:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Cluster marker evidence must be a pandas DataFrame.")
    if list(frame.columns) != MARKER_COLUMNS:
        raise ValueError(f"Cluster marker evidence must use the exact canonical columns in order: {MARKER_COLUMNS}.")

    for column in ("group", "gene"):
        values = frame[column].tolist()
        for value in values:
            if not isinstance(value, str):
                raise TypeError(f"Marker column {column!r} must contain strings only.")
            if not value.strip():
                raise ValueError(f"Marker column {column!r} cannot contain empty identifiers.")
            if value != value.strip():
                raise ValueError(f"Marker column {column!r} cannot contain surrounding whitespace.")

    if bool(frame.duplicated(subset=["group", "gene"]).any()):
        raise ValueError("Cluster marker evidence contains duplicate (group, gene) rows.")

    for column in _MARKER_NUMERIC_COLUMNS:
        series = frame[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"Marker column {column!r} must have a numeric dtype.")
        values = series.to_numpy(dtype=float)
        if values.size and not bool(np.isfinite(values).all()):
            raise ValueError(f"Marker column {column!r} must contain finite values.")

    ranks = frame["rank"].to_numpy(dtype=float)
    if ranks.size and (not bool(np.equal(ranks, np.floor(ranks)).all()) or bool((ranks < 1).any())):
        raise ValueError("Marker ranks must be positive integers.")
    for column in _MARKER_UNIT_INTERVAL_COLUMNS:
        values = frame[column].to_numpy(dtype=float)
        if values.size and bool(((values < 0.0) | (values > 1.0)).any()):
            raise ValueError(f"Marker column {column!r} must be between 0 and 1 inclusive.")
    p_values = frame["p_value"].to_numpy(dtype=float)
    adjusted = frame["p_adjusted"].to_numpy(dtype=float)
    tolerance = 1e-12 + 1e-7 * np.abs(p_values)
    if p_values.size and bool((adjusted + tolerance < p_values).any()):
        raise ValueError("Marker adjusted p-values cannot be smaller than raw p-values beyond numeric tolerance.")

    seen_groups: set[str] = set()
    active_group: str | None = None
    previous_rank = 0
    for group, rank in zip(frame["group"].tolist(), ranks, strict=True):
        if group != active_group:
            if group in seen_groups:
                raise ValueError("Marker rows for each group must form one contiguous block.")
            seen_groups.add(group)
            active_group = group
            previous_rank = 0
        integer_rank = int(rank)
        if integer_rank <= previous_rank:
            raise ValueError("Marker ranks must be strictly increasing within each group.")
        previous_rank = integer_rank


def _validate_universe_frame(frame: Any, *, np: Any, pd: Any) -> list[str]:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Tested-gene universe must be a pandas DataFrame.")
    if list(frame.columns) != MARKER_UNIVERSE_COLUMNS:
        raise ValueError(
            f"Tested-gene universe must use the exact canonical columns in order: {MARKER_UNIVERSE_COLUMNS}."
        )
    genes = _validate_identifier_values(frame["gene"].tolist(), name="Tested-gene universe")
    ranks = frame["universe_rank"]
    if pd.api.types.is_bool_dtype(ranks.dtype) or not pd.api.types.is_numeric_dtype(ranks.dtype):
        raise TypeError("Tested-gene universe ranks must have a numeric dtype.")
    numeric = ranks.to_numpy(dtype=float)
    if numeric.size and not bool(np.isfinite(numeric).all()):
        raise ValueError("Tested-gene universe ranks must be finite.")
    expected = np.arange(1, len(frame) + 1, dtype=float)
    if not bool(np.array_equal(numeric, expected)):
        raise ValueError("Tested-gene universe ranks must be consecutive, one-based, and order-stable.")
    return genes


def _all_values_integer_like(matrix: Any, *, np: Any, sparse: Any, chunk_size: int = 65_536) -> bool:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).reshape(-1)
    values = np.asarray(values)
    for start in range(0, int(values.size), chunk_size):
        chunk = values[start : start + chunk_size]
        if not bool(np.allclose(chunk, np.rint(chunk), rtol=0.0, atol=1e-8)):
            return False
    return True


def _frame_content_fingerprint(frame: Any, *, artifact: str, numeric_columns: set[str]) -> str:
    encoder = json.JSONEncoder(
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    digest = hashlib.sha256()

    def update(value: Any) -> None:
        for chunk in encoder.iterencode(value):
            digest.update(chunk.encode("utf-8"))

    digest.update(b'{"columns":')
    update(list(frame.columns))
    digest.update(b',"rows":[')
    first_row = True
    for row in frame.itertuples(index=False, name=None):
        encoded_row = []
        for column, value in zip(frame.columns, row, strict=True):
            if column in numeric_columns:
                encoded_row.append(float(value).hex())
            else:
                encoded_row.append(value)
        if not first_row:
            digest.update(b",")
        update(encoded_row)
        first_row = False
    digest.update(b'],"schema":')
    update(f"openbio-singlecell/{artifact}/v2")
    digest.update(b"}")
    return digest.hexdigest()


def _marker_table_content_fingerprint(frame: Any, *, artifact: str = "marker-table") -> str:
    return _frame_content_fingerprint(frame, artifact=artifact, numeric_columns=set(_MARKER_NUMERIC_COLUMNS))


def _universe_content_fingerprint(frame: Any) -> str:
    return _frame_content_fingerprint(
        frame,
        artifact="tested-gene-universe",
        numeric_columns={"universe_rank"},
    )


def marker_table_content_fingerprint(frame: Any) -> str:
    """Return the canonical SHA-256 identity of the current marker-table content."""
    return _marker_table_content_fingerprint(frame)


def marker_universe_content_fingerprint(frame: Any) -> str:
    """Return the canonical SHA-256 identity of the current tested-universe content."""
    return _universe_content_fingerprint(frame)


def _benjamini_hochberg(p_values: Any, *, np: Any, multipletests: Any) -> Any:
    values = np.asarray(p_values, dtype=float)
    if not values.size:
        return values.copy()
    adjusted = np.asarray(
        multipletests(values, alpha=0.05, method="fdr_bh", is_sorted=False, returnsorted=False)[1],
        dtype=float,
    )
    return np.maximum(adjusted, values)


def _working_memory_estimate(
    *,
    n_obs: int,
    n_vars: int,
    n_variable_genes: int,
    n_groups: int,
    actual_n_genes: int,
    sparse_expression: bool,
) -> dict[str, Any]:
    gib = float(1024**3)
    validation_bytes = n_obs * 256 + n_vars * 256 + 2 * 65_536 * 8
    if sparse_expression:
        selected_expression_bytes = n_obs * n_variable_genes * 16 + (n_obs + 1) * 8
        expression_policy = "worst-case sparse values/indices plus row pointers"
    else:
        selected_expression_bytes = n_obs * n_variable_genes * 8
        expression_policy = "float64 dense selected-expression copy"
    backend_arrays_bytes = n_groups * n_variable_genes * 8 * 8
    full_ranking_bytes = n_groups * n_vars * 768
    marker_output_bytes = n_groups * actual_n_genes * 256
    universe_output_bytes = n_vars * 160
    total_bytes = (
        validation_bytes
        + selected_expression_bytes
        + backend_arrays_bytes
        + full_ranking_bytes
        + marker_output_bytes
        + universe_output_bytes
    )
    return {
        "bounded_validation_buffers_gib": validation_bytes / gib,
        "selected_expression_temp_copy_gib": selected_expression_bytes / gib,
        "backend_group_by_variable_gene_arrays_gib": backend_arrays_bytes / gib,
        "complete_ranking_table_gib": full_ranking_bytes / gib,
        "marker_output_table_gib": marker_output_bytes / gib,
        "universe_output_table_gib": universe_output_bytes / gib,
        "estimated_peak_working_memory_gib": total_bytes / gib,
        "estimation_policy": (
            "Bounded extrema/integer-like validation buffers and streaming SHA-256 serialization; conservative "
            f"{expression_policy}; eight backend group-by-variable-gene numeric/prevalence arrays; 768 bytes per "
            "complete-ranking row, 256 bytes per exported-marker row, and 160 bytes per universe row."
        ),
    }


def _marker_analysis_fingerprint(
    *,
    groupby: str,
    method: str,
    source_kind: str,
    layer_name: str | None,
    n_genes: int,
    tie_correct: bool,
    genes: list[str],
    observation_ids: list[str],
    group_labels_by_observation: list[str],
) -> tuple[str, str]:
    universe_fingerprint = _canonical_sha256(
        {
            "schema": "openbio-singlecell/tested-gene-universe-identity/v2",
            "ordered_genes": genes,
        }
    )
    analysis_fingerprint = _canonical_sha256(
        {
            "schema": "openbio-singlecell/cluster-marker-analysis/v2",
            "groupby": groupby,
            "method": method,
            "source": source_kind,
            "layer_name": layer_name,
            "n_genes": n_genes,
            "tie_correct": tie_correct,
            "universe_fingerprint": universe_fingerprint,
            "ordered_observation_ids": observation_ids,
            "group_labels_by_observation": group_labels_by_observation,
            "fixed_policy": {
                "groups": "all",
                "reference": "rest",
                "rankby_abs": False,
                "pts": True,
                "corr_method": "benjamini-hochberg",
            },
        }
    )
    return analysis_fingerprint, universe_fingerprint


def _rank_marker_evidence_impl(
    adata: Any,
    *,
    groupby: str,
    method: str,
    source_kind: str,
    layer_name: str | None,
    n_genes: int,
    tie_correct: bool,
    max_output_rows: int,
    max_working_memory_gib: float,
    np: Any,
    pd: Any,
    sparse: Any,
    sc: Any,
    ad: Any,
    multipletests: Any,
) -> tuple[Any, Any, dict[str, Any]]:
    if getattr(adata, "isbacked", False):
        raise ValueError("Marker Genes requires an in-memory AnnData.")
    n_obs = int(adata.n_obs)
    n_vars = int(adata.n_vars)
    if n_obs < 1 or n_vars < 1:
        raise ValueError("Marker Genes requires a nonempty AnnData.")
    max_working_memory_gib = _validate_positive_float(
        max_working_memory_gib,
        name="Marker Genes max_working_memory_gib",
        np=np,
    )
    identifier_group_validation_gib = (n_obs * 256 + n_vars * 256 + 2 * 65_536 * 8) / float(1024**3)
    if identifier_group_validation_gib > max_working_memory_gib:
        raise ValueError(
            "Marker Genes requires an estimated "
            f"{identifier_group_validation_gib:.6g} GiB for identifier/group validation buffers, exceeding "
            f"max_working_memory_gib={max_working_memory_gib:.6g}."
        )
    if not bool(adata.obs_names.is_unique):
        raise ValueError("Marker Genes requires unique observation identifiers.")
    observation_ids = _validate_identifier_values(adata.obs_names, name="Observation")
    genes = _validate_identifier_values(adata.var_names, name="Selected-source gene")

    if not isinstance(groupby, str):
        raise TypeError("Marker Genes groupby must be a string.")
    groupby = groupby.strip()
    if not groupby:
        raise ValueError("Marker Genes groupby cannot be empty.")
    if groupby not in adata.obs:
        raise ValueError(f"Marker groupby column not found in obs: {groupby!r}")
    grouping = adata.obs[groupby]
    if not isinstance(grouping.dtype, pd.CategoricalDtype):
        raise TypeError("Marker Genes grouping must be categorical with explicit category order.")
    category_codes = np.asarray(grouping.cat.codes.to_numpy(copy=False), dtype=np.int64)
    if category_codes.size and int(category_codes.min()) < 0:
        raise ValueError("Marker Genes grouping cannot contain missing labels.")
    categories = list(grouping.cat.categories)
    category_counts = np.bincount(category_codes, minlength=len(categories))
    if bool((category_counts == 0).any()):
        unused = [str(category) for category, count in zip(categories, category_counts, strict=True) if count == 0]
        raise ValueError(f"Marker Genes grouping contains unused categories: {unused}.")
    group_labels = []
    for category in categories:
        rendered = str(category)
        if not rendered.strip():
            raise ValueError("Marker Genes group labels cannot render as empty strings.")
        if rendered != rendered.strip():
            raise ValueError("Marker Genes group labels cannot contain surrounding whitespace.")
        group_labels.append(rendered)
    if len(set(group_labels)) != len(group_labels):
        raise ValueError("Marker Genes group labels collide after string serialization.")
    if len(group_labels) < 2:
        raise ValueError("Marker Genes requires at least two observed groups.")
    labels_by_observation = [group_labels[int(code)] for code in category_codes]
    group_sizes = {
        label: int(count) for label, count in zip(group_labels, category_counts.tolist(), strict=True)
    }
    undersized = {label: size for label, size in group_sizes.items() if size < 2}
    if undersized:
        raise ValueError(f"Marker Genes requires at least two cells in every group: {undersized}.")
    total_cells = len(labels_by_observation)
    small_references = {label: total_cells - size for label, size in group_sizes.items() if total_cells - size < 2}
    if small_references:
        raise ValueError(f"Marker Genes requires at least two cells in every rest reference: {small_references}.")

    if method not in MARKER_METHODS:
        raise ValueError(f"Unsupported Marker Genes method: {method!r}; expected one of {list(MARKER_METHODS)}.")
    if not isinstance(tie_correct, bool):
        raise TypeError("Marker Genes tie_correct must be a boolean.")
    if method != "wilcoxon" and tie_correct:
        raise ValueError("Marker Genes tie_correct is only applicable to the Wilcoxon method; set it to false.")
    n_genes = _validate_integer(n_genes, name="Marker Genes n_genes", allow_zero=True)
    if n_genes > int(adata.n_vars):
        raise ValueError("Marker Genes n_genes cannot exceed the selected-source gene count.")
    max_output_rows = _validate_integer(max_output_rows, name="Marker Genes max_output_rows")
    actual_n_genes = int(adata.n_vars) if n_genes == 0 else n_genes
    expected_marker_rows = len(group_labels) * actual_n_genes
    expected_universe_rows = int(adata.n_vars)
    expected_output_rows = expected_marker_rows + expected_universe_rows
    if expected_output_rows > max_output_rows:
        raise ValueError(
            f"Marker Genes would produce {expected_output_rows:,} rows across marker and universe outputs "
            f"({expected_marker_rows:,} + {expected_universe_rows:,}), exceeding "
            f"max_output_rows={max_output_rows:,}."
        )

    if source_kind == "X":
        if layer_name is not None:
            raise ValueError("Marker Genes X source cannot carry a layer name.")
        matrix = adata.X
    elif source_kind == "layer":
        if not isinstance(layer_name, str) or not layer_name.strip():
            raise ValueError("Marker Genes layer source requires a nonempty layer name.")
        layer_name = layer_name.strip()
        if layer_name not in adata.layers:
            raise ValueError(f"Marker source layer not found: {layer_name!r}")
        matrix = adata.layers[layer_name]
    else:
        raise ValueError("Marker Genes expression source must be X or layer; Raw is unsupported.")
    if getattr(matrix, "shape", None) != (int(adata.n_obs), int(adata.n_vars)):
        raise ValueError("Marker Genes expression source is not aligned to AnnData observations and variables.")
    matrix_dtype = getattr(matrix, "dtype", None)
    if matrix_dtype is None or not np.issubdtype(matrix_dtype, np.number) or np.issubdtype(
        matrix_dtype, np.complexfloating
    ):
        raise TypeError("Marker Genes expression must contain real numeric values.")
    if sparse.issparse(matrix):
        minima_result = matrix.min(axis=0)
        maxima_result = matrix.max(axis=0)
        minima = np.asarray(
            minima_result.toarray() if sparse.issparse(minima_result) else minima_result,
            dtype=float,
        ).ravel()
        maxima = np.asarray(
            maxima_result.toarray() if sparse.issparse(maxima_result) else maxima_result,
            dtype=float,
        ).ravel()
    else:
        dense_matrix = np.asarray(matrix)
        minima = np.asarray(dense_matrix.min(axis=0), dtype=float).ravel()
        maxima = np.asarray(dense_matrix.max(axis=0), dtype=float).ravel()
    if not bool(np.isfinite(minima).all()) or not bool(np.isfinite(maxima).all()):
        raise ValueError("Marker Genes expression contains non-finite values.")
    contains_negative_values = bool((minima < 0).any())
    constant_mask = minima == maxima
    constant_indices = np.flatnonzero(constant_mask)
    variable_indices = np.flatnonzero(~constant_mask)

    memory = _working_memory_estimate(
        n_obs=int(adata.n_obs),
        n_vars=int(adata.n_vars),
        n_variable_genes=int(variable_indices.size),
        n_groups=len(group_labels),
        actual_n_genes=actual_n_genes,
        sparse_expression=bool(sparse.issparse(matrix)),
    )
    if memory["estimated_peak_working_memory_gib"] > max_working_memory_gib:
        raise ValueError(
            "Marker Genes requires an estimated "
            f"{memory['estimated_peak_working_memory_gib']:.6g} GiB working memory, exceeding "
            f"max_working_memory_gib={max_working_memory_gib:.6g}."
        )
    count_like_values = bool(_all_values_integer_like(matrix, np=np, sparse=sparse))

    analysis_fingerprint, universe_fingerprint = _marker_analysis_fingerprint(
        groupby=groupby,
        method=method,
        source_kind=source_kind,
        layer_name=layer_name,
        n_genes=n_genes,
        tie_correct=tie_correct,
        genes=genes,
        observation_ids=observation_ids,
        group_labels_by_observation=labels_by_observation,
    )
    variable_genes = [genes[int(index)] for index in variable_indices]
    expected_backend_rows = len(group_labels) * len(variable_genes)
    if variable_genes:
        selected_matrix = (
            matrix[:, variable_indices].copy()
            if sparse.issparse(matrix)
            else np.asarray(matrix)[:, variable_indices].copy()
        )
        work = ad.AnnData(
            X=selected_matrix,
            obs=pd.DataFrame(
                {
                    groupby: pd.Categorical(
                        labels_by_observation,
                        categories=group_labels,
                        ordered=bool(grouping.cat.ordered),
                    )
                },
                index=observation_ids,
            ),
            var=pd.DataFrame(index=variable_genes),
        )
        result_key = "_openbio_marker_genes"
        sc.tl.rank_genes_groups(
            work,
            groupby=groupby,
            groups="all",
            reference="rest",
            method=method,
            corr_method="benjamini-hochberg",
            n_genes=len(variable_genes),
            rankby_abs=False,
            pts=True,
            tie_correct=tie_correct if method == "wilcoxon" else False,
            use_raw=False,
            layer=None,
            key_added=result_key,
        )
        stored = work.uns.get(result_key)
        if not isinstance(stored, Mapping):
            raise RuntimeError("Marker Genes backend did not return a result mapping.")
        required_stored = {
            "params",
            "names",
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
            "pts",
            "pts_rest",
        }
        missing_stored = sorted(required_stored.difference(stored))
        if missing_stored:
            raise RuntimeError(f"Marker Genes backend result is incomplete: {missing_stored}.")
        backend_parameters = stored.get("params")
        if not isinstance(backend_parameters, Mapping):
            raise RuntimeError("Marker Genes backend did not return parameter metadata.")
        expected_backend_parameters = {
            "groupby": groupby,
            "reference": "rest",
            "method": method,
            "use_raw": False,
            "layer": None,
            "corr_method": "benjamini-hochberg",
        }
        for parameter, expected in expected_backend_parameters.items():
            if backend_parameters.get(parameter) != expected:
                raise RuntimeError(
                    f"Marker Genes backend reported inconsistent {parameter!r}: {backend_parameters.get(parameter)!r}."
                )

        ranked = sc.get.rank_genes_groups_df(work, group=None, key=result_key)
        if not isinstance(ranked, pd.DataFrame):
            raise RuntimeError("Marker Genes backend extractor did not return a pandas DataFrame.")
        expected_backend_columns = [
            "group",
            "names",
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
            "pct_nz_group",
            "pct_nz_reference",
        ]
        missing_columns = [column for column in expected_backend_columns if column not in ranked.columns]
        if missing_columns:
            raise RuntimeError(f"Marker Genes backend table is missing columns: {missing_columns}.")
        if len(ranked) != expected_backend_rows:
            raise RuntimeError(
                f"Marker Genes backend returned {len(ranked):,} rows; expected exactly "
                f"{expected_backend_rows:,} variable-gene rows."
            )
        for column in (
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
            "pct_nz_group",
            "pct_nz_reference",
        ):
            series = ranked[column]
            if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
                raise RuntimeError(f"Marker Genes backend column {column!r} must have a numeric dtype.")
            numeric = series.to_numpy(dtype=float)
            if numeric.size and not bool(np.isfinite(numeric).all()):
                raise RuntimeError(f"Marker Genes backend column {column!r} must contain finite values.")
    else:
        ranked = pd.DataFrame(
            columns=[
                "group",
                "names",
                "scores",
                "logfoldchanges",
                "pvals",
                "pvals_adj",
                "pct_nz_group",
                "pct_nz_reference",
            ]
        )

    gene_order = {gene: index for index, gene in enumerate(genes)}
    constant_positive = {genes[int(index)]: bool(minima[int(index)] > 0.0) for index in constant_indices}
    complete_blocks = []
    for group in group_labels:
        backend_group = ranked.loc[ranked["group"] == group]
        backend_genes = backend_group["names"].tolist()
        if any(not isinstance(gene, str) for gene in backend_genes):
            raise RuntimeError("Marker Genes backend returned non-string gene identifiers.")
        if len(backend_group) != len(variable_genes) or set(backend_genes) != set(variable_genes):
            raise RuntimeError(f"Marker Genes backend returned an incomplete variable-gene family for {group!r}.")
        records = [
            {
                "group": group,
                "gene": gene,
                "score": float(score),
                "log2_fold_change_approx": float(logfc),
                "p_value": float(p_value),
                "p_adjusted": 1.0,
                "fraction_in_group": float(fraction_in),
                "fraction_reference": float(fraction_reference),
            }
            for gene, score, logfc, p_value, fraction_in, fraction_reference in zip(
                backend_group["names"].tolist(),
                backend_group["scores"].tolist(),
                backend_group["logfoldchanges"].tolist(),
                backend_group["pvals"].tolist(),
                backend_group["pct_nz_group"].tolist(),
                backend_group["pct_nz_reference"].tolist(),
                strict=True,
            )
        ]
        records.extend(
            {
                "group": group,
                "gene": gene,
                "score": 0.0,
                "log2_fold_change_approx": 0.0,
                "p_value": 1.0,
                "p_adjusted": 1.0,
                "fraction_in_group": float(constant_positive[gene]),
                "fraction_reference": float(constant_positive[gene]),
            }
            for gene in (genes[int(index)] for index in constant_indices)
        )
        raw_p_values = [record["p_value"] for record in records]
        adjusted_values = _benjamini_hochberg(raw_p_values, np=np, multipletests=multipletests)
        for record, adjusted_value in zip(records, adjusted_values, strict=True):
            record["p_adjusted"] = float(adjusted_value)
        records.sort(key=lambda record: (-record["score"], gene_order[record["gene"]]))
        for rank, record in enumerate(records, start=1):
            record["rank"] = rank
        complete_blocks.append(pd.DataFrame.from_records(records, columns=MARKER_COLUMNS))

    complete_ranking = pd.concat(complete_blocks, ignore_index=True)
    try:
        _validate_canonical_marker_table(complete_ranking, np=np, pd=pd)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Marker Genes backend returned an invalid canonical table: {error}") from error
    for group in group_labels:
        returned_genes = complete_ranking.loc[complete_ranking["group"] == group, "gene"].tolist()
        if len(returned_genes) != len(genes) or set(returned_genes) != set(genes):
            raise RuntimeError(f"Marker Genes did not preserve the complete tested-gene family for {group!r}.")

    ranking_fingerprint = _marker_table_content_fingerprint(
        complete_ranking,
        artifact="complete-marker-ranking",
    )
    table = complete_ranking.groupby("group", sort=False, group_keys=False).head(actual_n_genes).reset_index(drop=True)
    _validate_canonical_marker_table(table, np=np, pd=pd)

    universe = pd.DataFrame(
        {
            "gene": genes,
            "universe_rank": np.arange(1, len(genes) + 1, dtype=int),
        },
        columns=MARKER_UNIVERSE_COLUMNS,
    )
    _validate_universe_frame(universe, np=np, pd=pd)
    table_content_fingerprint = _marker_table_content_fingerprint(table)
    universe_content_fingerprint = _universe_content_fingerprint(universe)
    warnings_list = []
    if count_like_values and not contains_negative_values:
        warnings_list.append(
            "The selected expression is integer-like/count-like. Marker statistics remain calculable, but library "
            "size, ties, prevalence, and approximate fold-change interpretation can differ from log-expression practice."
        )
    if contains_negative_values:
        warnings_list.append(
            "The selected expression contains negative values. Marker ranking remains calculable, but `>0` "
            "prevalence and approximate log2-fold-change fields do not have standard nonnegative log-abundance semantics."
        )
    if constant_indices.size:
        warnings_list.append(
            f"Retained {int(constant_indices.size):,} globally constant genes as neutral hypotheses "
            "(score/log2 fold change 0; p/p-adjusted 1) with observed prevalence."
        )
    small_groups = {label: size for label, size in group_sizes.items() if size < 20}
    if small_groups:
        warnings_list.append(
            f"Small groups can yield unstable Cluster marker evidence; groups with fewer than 20 cells: {small_groups}."
        )
    details = {
        "analysis_fingerprint": analysis_fingerprint,
        "ranking_fingerprint": ranking_fingerprint,
        "universe_fingerprint": universe_fingerprint,
        "table_content_fingerprint": table_content_fingerprint,
        "universe_content_fingerprint": universe_content_fingerprint,
        "group_labels": group_labels,
        "group_sizes": group_sizes,
        "source_gene_count": len(genes),
        "requested_n_genes": n_genes,
        "actual_n_genes_per_group": actual_n_genes,
        "ranking_truncated": actual_n_genes < len(genes),
        "constant_gene_count": int(constant_indices.size),
        "variable_gene_count": int(variable_indices.size),
        "backend_hypotheses_per_group": int(variable_indices.size),
        "bh_hypotheses_per_group": len(genes),
        "marker_output_rows": expected_marker_rows,
        "universe_output_rows": expected_universe_rows,
        "total_output_rows": expected_output_rows,
        "max_working_memory_gib": max_working_memory_gib,
        "working_memory": memory,
        "warnings": warnings_list,
    }
    return table, universe, details


def marker_provenance_parameters(
    *,
    artifact_role: str,
    groupby: str,
    method: str,
    source_kind: str,
    layer_name: str | None,
    n_genes: int,
    tie_correct: bool,
    max_output_rows: int,
    max_working_memory_gib: float,
    details: Mapping[str, Any],
) -> dict[str, Any]:
    if artifact_role == "marker_table":
        content_fingerprint = details["table_content_fingerprint"]
    elif artifact_role == "tested_gene_universe":
        content_fingerprint = details["universe_content_fingerprint"]
    else:
        raise ValueError(f"Unsupported marker artifact role: {artifact_role!r}.")
    parameters: dict[str, Any] = {
        "marker_evidence_schema_version": MARKER_EVIDENCE_SCHEMA_VERSION,
        "producer_node_id": MARKER_GENES_NODE_ID,
        "artifact_role": artifact_role,
        "analysis_fingerprint": details["analysis_fingerprint"],
        "ranking_fingerprint": details["ranking_fingerprint"],
        "universe_fingerprint": details["universe_fingerprint"],
        "content_fingerprint": content_fingerprint,
        "marker_groupby": groupby,
        "marker_method": method,
        "marker_source": source_kind,
        "group_labels": list(details["group_labels"]),
        "group_sizes": dict(details["group_sizes"]),
        "source_gene_count": int(details["source_gene_count"]),
        "requested_n_genes": n_genes,
        "actual_n_genes_per_group": int(details["actual_n_genes_per_group"]),
        "ranking_truncated": bool(details["ranking_truncated"]),
        "tie_correct": tie_correct,
        "max_output_rows": max_output_rows,
        "max_working_memory_gib": max_working_memory_gib,
        "constant_gene_count": int(details["constant_gene_count"]),
        "variable_gene_count": int(details["variable_gene_count"]),
        "bh_hypotheses_per_group": int(details["bh_hypotheses_per_group"]),
        "marker_output_rows": int(details["marker_output_rows"]),
        "universe_output_rows": int(details["universe_output_rows"]),
        "total_output_rows": int(details["total_output_rows"]),
        "groups": "all",
        "reference": "rest",
        "rankby_abs": False,
        "pts": True,
        "corr_method": "benjamini-hochberg",
    }
    if layer_name is not None:
        parameters["marker_layer_name"] = layer_name
    return parameters


def _source_parameters(result: TableResult, *, artifact_name: str) -> Mapping[str, Any]:
    if not isinstance(result.source, Mapping):
        raise ValueError(f"{artifact_name} provenance source must be a mapping.")
    source_parameters = result.source.get("parameters")
    if not isinstance(source_parameters, Mapping):
        raise ValueError(f"{artifact_name} provenance parameters must be a mapping.")
    if dict(source_parameters) != result.parameters:
        raise ValueError(f"{artifact_name} parameters disagree with its immutable analysis source.")
    return source_parameters


def validate_marker_artifact_pair(
    table: Any,
    universe: Any,
    *,
    np: Any,
    pd: Any,
    allowed_table_operations: tuple[str, ...] = ("marker_genes", "filter_marker_genes"),
) -> dict[str, Any]:
    if not isinstance(table, TableResult) or not isinstance(universe, TableResult):
        raise TypeError("Marker-evidence consumers require TableResult values for table and universe.")
    _validate_canonical_marker_table(table.table, np=np, pd=pd)
    genes = _validate_universe_frame(universe.table, np=np, pd=pd)
    if not genes:
        raise ValueError("Tested-gene universe cannot be empty.")

    table_parameters = _source_parameters(table, artifact_name="Marker table")
    universe_parameters = _source_parameters(universe, artifact_name="Tested-gene universe")
    table_operation = table.source.get("operation")
    known_table_producers = {
        "marker_genes": MARKER_GENES_NODE_ID,
        "filter_marker_genes": FILTER_MARKER_GENES_NODE_ID,
    }
    if (
        not isinstance(allowed_table_operations, tuple)
        or not allowed_table_operations
        or any(operation not in known_table_producers for operation in allowed_table_operations)
    ):
        raise ValueError("Marker consumer configured invalid allowed_table_operations.")
    if table_operation not in allowed_table_operations:
        raise ValueError(
            f"Marker consumer requires table operation in {list(allowed_table_operations)}; "
            f"received {table_operation!r}."
        )
    expected_table_producer = known_table_producers.get(table_operation)
    if expected_table_producer is None or table_parameters.get("producer_node_id") != expected_table_producer:
        raise ValueError("Marker table provenance does not identify an approved marker-evidence producer.")
    if universe.source.get("operation") != "marker_genes":
        raise ValueError("Tested-gene universe must come directly from Marker Genes.")
    if universe_parameters.get("producer_node_id") != MARKER_GENES_NODE_ID:
        raise ValueError("Tested-gene universe provenance has the wrong producer node.")
    if table_parameters.get("artifact_role") != "marker_table":
        raise ValueError("Marker table provenance has the wrong artifact role.")
    if universe_parameters.get("artifact_role") != "tested_gene_universe":
        raise ValueError("Tested-gene universe provenance has the wrong artifact role.")
    for parameters, artifact_name in (
        (table_parameters, "Marker table"),
        (universe_parameters, "Tested-gene universe"),
    ):
        if parameters.get("marker_evidence_schema_version") != MARKER_EVIDENCE_SCHEMA_VERSION:
            raise ValueError(f"{artifact_name} uses an unsupported marker-evidence schema.")
        for field in (
            "analysis_fingerprint",
            "ranking_fingerprint",
            "universe_fingerprint",
            "content_fingerprint",
        ):
            value = parameters.get(field)
            if (
                not isinstance(value, str)
                or len(value) != _SHA256_LENGTH
                or any(character not in "0123456789abcdef" for character in value)
            ):
                raise ValueError(f"{artifact_name} has an invalid SHA-256 {field}.")

    analysis_fingerprint = table_parameters["analysis_fingerprint"]
    if analysis_fingerprint != universe_parameters["analysis_fingerprint"]:
        raise ValueError("Marker table and tested-gene universe analysis fingerprints do not match.")
    ranking_fingerprint = table_parameters["ranking_fingerprint"]
    if ranking_fingerprint != universe_parameters["ranking_fingerprint"]:
        raise ValueError("Marker table and tested-gene universe ranking fingerprints do not match.")
    shared_fields = (
        "marker_groupby",
        "marker_method",
        "marker_source",
        "marker_layer_name",
        "group_labels",
        "group_sizes",
        "source_gene_count",
        "requested_n_genes",
        "actual_n_genes_per_group",
        "ranking_truncated",
        "tie_correct",
        "max_working_memory_gib",
        "constant_gene_count",
        "variable_gene_count",
        "bh_hypotheses_per_group",
        "marker_output_rows",
        "universe_output_rows",
        "total_output_rows",
        "groups",
        "reference",
        "rankby_abs",
        "pts",
        "corr_method",
    )
    mismatched_fields = [
        field for field in shared_fields if table_parameters.get(field) != universe_parameters.get(field)
    ]
    if mismatched_fields:
        raise ValueError(f"Marker table and tested-gene universe provenance disagree: {mismatched_fields}.")
    universe_fingerprint = _canonical_sha256(
        {
            "schema": "openbio-singlecell/tested-gene-universe-identity/v2",
            "ordered_genes": genes,
        }
    )
    if universe_parameters["universe_fingerprint"] != universe_fingerprint:
        raise ValueError("Tested-gene universe content does not match its SHA-256 identity.")
    if table_parameters["universe_fingerprint"] != universe_fingerprint:
        raise ValueError("Marker table provenance does not match the tested-gene universe identity.")
    table_content_fingerprint = _marker_table_content_fingerprint(table.table)
    if table_parameters["content_fingerprint"] != table_content_fingerprint:
        raise ValueError("Marker table content does not match its SHA-256 current-content fingerprint.")
    universe_content_fingerprint = _universe_content_fingerprint(universe.table)
    if universe_parameters["content_fingerprint"] != universe_content_fingerprint:
        raise ValueError("Tested-gene universe content does not match its SHA-256 current-content fingerprint.")
    if not set(table.table["gene"].tolist()).issubset(set(genes)):
        raise ValueError("Marker table contains genes outside the bound tested-gene universe.")
    if int(table.input_genes) != len(genes) or int(universe.input_genes) != len(genes):
        raise ValueError("Marker artifacts disagree with the tested-gene universe size.")

    group_labels = table_parameters.get("group_labels")
    if not isinstance(group_labels, list) or any(not isinstance(label, str) or not label for label in group_labels):
        raise ValueError("Marker table provenance must contain ordered string group labels.")
    if len(group_labels) < 2 or len(set(group_labels)) != len(group_labels):
        raise ValueError("Marker table provenance must contain at least two unique group labels.")
    group_sizes = table_parameters.get("group_sizes")
    if not isinstance(group_sizes, Mapping) or set(group_sizes) != set(group_labels):
        raise ValueError("Marker table provenance has invalid group-size keys.")
    if any(isinstance(size, bool) or not isinstance(size, int) or size < 2 for size in group_sizes.values()):
        raise ValueError("Marker table provenance group sizes must be integers of at least two.")
    source_gene_count = table_parameters.get("source_gene_count")
    if isinstance(source_gene_count, bool) or not isinstance(source_gene_count, int) or source_gene_count != len(genes):
        raise ValueError("Marker table provenance has an invalid source gene count.")
    if table_parameters.get("marker_source") not in {"X", "layer"}:
        raise ValueError("Marker table provenance has an invalid expression source.")
    if not isinstance(table_parameters.get("marker_groupby"), str) or not table_parameters["marker_groupby"].strip():
        raise ValueError("Marker table provenance has an invalid grouping column.")
    observed_groups = list(dict.fromkeys(table.table["group"].tolist()))
    if any(group not in group_labels for group in observed_groups):
        raise ValueError("Marker table contains a group absent from its provenance.")
    expected_observed_order = [group for group in group_labels if group in set(observed_groups)]
    if observed_groups != expected_observed_order:
        raise ValueError("Marker table group blocks do not follow the provenance category order.")
    method = table_parameters.get("marker_method")
    if method not in MARKER_METHODS:
        raise ValueError("Marker table provenance contains an unsupported marker method.")
    actual_n_genes = table_parameters.get("actual_n_genes_per_group")
    if isinstance(actual_n_genes, bool) or not isinstance(actual_n_genes, int) or actual_n_genes < 1:
        raise ValueError("Marker table provenance has an invalid actual_n_genes_per_group value.")
    if len(table.table) and int(table.table["rank"].max()) > actual_n_genes:
        raise ValueError("Marker table contains a rank outside its upstream ranking depth.")
    if table_operation == "marker_genes":
        expected_rows = len(group_labels) * actual_n_genes
        if len(table.table) != expected_rows:
            raise ValueError("Direct Marker Genes evidence has an incomplete group-by-rank grid.")
        for group in group_labels:
            ranks = table.table.loc[table.table["group"] == group, "rank"].tolist()
            if ranks != list(range(1, actual_n_genes + 1)):
                raise ValueError(f"Direct Marker Genes evidence has incomplete ranks for group {group!r}.")
    else:
        if table_parameters.get("upstream_producer_node_id") != MARKER_GENES_NODE_ID:
            raise ValueError("Filtered marker evidence must identify Marker Genes as its direct upstream producer.")
        upstream_content_fingerprint = table_parameters.get("upstream_content_fingerprint")
        if (
            not isinstance(upstream_content_fingerprint, str)
            or len(upstream_content_fingerprint) != _SHA256_LENGTH
            or any(character not in "0123456789abcdef" for character in upstream_content_fingerprint)
        ):
            raise ValueError("Filtered marker evidence has an invalid upstream content fingerprint.")
    return {
        "analysis_fingerprint": analysis_fingerprint,
        "ranking_fingerprint": ranking_fingerprint,
        "universe_fingerprint": universe_fingerprint,
        "table_content_fingerprint": table_content_fingerprint,
        "universe_content_fingerprint": universe_content_fingerprint,
        "table_operation": table_operation,
        "group_labels": list(group_labels),
        "marker_method": method,
        "ranking_truncated": bool(table_parameters.get("ranking_truncated")),
        "upstream_parameters": dict(table_parameters),
    }


def _filter_marker_table_impl(
    frame: Any,
    *,
    min_log2_fold_change: float,
    min_fraction_in_group: float,
    max_fraction_reference: float,
    max_p_adjusted: float,
    np: Any,
    pd: Any,
) -> tuple[Any, dict[str, Any]]:
    min_log2_fold_change = _validate_finite_threshold(
        min_log2_fold_change,
        name="min_log2_fold_change",
        np=np,
    )
    min_fraction_in_group = _validate_finite_threshold(
        min_fraction_in_group,
        name="min_fraction_in_group",
        np=np,
        unit_interval=True,
    )
    max_fraction_reference = _validate_finite_threshold(
        max_fraction_reference,
        name="max_fraction_reference",
        np=np,
        unit_interval=True,
    )
    max_p_adjusted = _validate_finite_threshold(
        max_p_adjusted,
        name="max_p_adjusted",
        np=np,
        unit_interval=True,
    )
    _validate_canonical_marker_table(frame, np=np, pd=pd)
    logfc_mask = frame["log2_fold_change_approx"].to_numpy(dtype=float) >= min_log2_fold_change
    fraction_in_mask = frame["fraction_in_group"].to_numpy(dtype=float) >= min_fraction_in_group
    fraction_reference_mask = frame["fraction_reference"].to_numpy(dtype=float) <= max_fraction_reference
    adjusted_p_mask = frame["p_adjusted"].to_numpy(dtype=float) <= max_p_adjusted
    retained_mask = logfc_mask & fraction_in_mask & fraction_reference_mask & adjusted_p_mask
    filtered = frame.loc[retained_mask, MARKER_COLUMNS].copy().reset_index(drop=True)
    _validate_canonical_marker_table(filtered, np=np, pd=pd)
    diagnostics = {
        "input_rows": len(frame),
        "retained_rows": len(filtered),
        "rejected_rows": int((~retained_mask).sum()),
        "failed_min_log2_fold_change": int((~logfc_mask).sum()),
        "failed_min_fraction_in_group": int((~fraction_in_mask).sum()),
        "failed_max_fraction_reference": int((~fraction_reference_mask).sum()),
        "failed_max_p_adjusted": int((~adjusted_p_mask).sum()),
        "content_fingerprint": _marker_table_content_fingerprint(filtered),
    }
    return filtered, diagnostics


def rank_marker_evidence(
    adata: Any,
    *,
    groupby: str,
    method: str,
    source_kind: str,
    layer_name: str | None,
    n_genes: int,
    tie_correct: bool,
    max_output_rows: int,
    max_working_memory_gib: float,
    np: Any,
    pd: Any,
    sparse: Any,
    sc: Any,
    ad: Any,
) -> tuple[DataFrame, DataFrame, dict[str, Any]]:
    from statsmodels.stats.multitest import multipletests

    return _rank_marker_evidence_impl(
        adata,
        groupby=groupby,
        method=method,
        source_kind=source_kind,
        layer_name=layer_name,
        n_genes=n_genes,
        tie_correct=tie_correct,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        np=np,
        pd=pd,
        sparse=sparse,
        sc=sc,
        ad=ad,
        multipletests=multipletests,
    )


def filter_marker_table(
    frame: Any,
    *,
    min_log2_fold_change: float,
    min_fraction_in_group: float,
    max_fraction_reference: float,
    max_p_adjusted: float,
    np: Any,
    pd: Any,
) -> tuple[DataFrame, dict[str, Any]]:
    return _filter_marker_table_impl(
        frame,
        min_log2_fold_change=min_log2_fold_change,
        min_fraction_in_group=min_fraction_in_group,
        max_fraction_reference=max_fraction_reference,
        max_p_adjusted=max_p_adjusted,
        np=np,
        pd=pd,
    )


def _validate_generated_marker_frames(frame: Any, universe: Any, *, np: Any, pd: Any) -> None:
    _validate_canonical_marker_table(frame, np=np, pd=pd)
    genes = _validate_universe_frame(universe, np=np, pd=pd)
    if not genes:
        raise ValueError("Tested-gene universe cannot be empty.")
    if not set(frame["gene"].tolist()).issubset(set(genes)):
        raise ValueError("Marker table contains genes outside the tested-gene universe.")


def marker_genes_code(
    *,
    groupby: str,
    method: str,
    source_kind: str,
    layer_name: str | None,
    n_genes: int,
    tie_correct: bool,
    max_output_rows: int,
    max_working_memory_gib: float,
) -> str:
    helpers = [
        _canonical_sha256,
        _validate_identifier_values,
        _validate_integer,
        _validate_positive_float,
        _validate_canonical_marker_table,
        _validate_universe_frame,
        _all_values_integer_like,
        _frame_content_fingerprint,
        _marker_table_content_fingerprint,
        _universe_content_fingerprint,
        _benjamini_hochberg,
        _working_memory_estimate,
        _marker_analysis_fingerprint,
        _rank_marker_evidence_impl,
    ]
    sources = "\n\n".join(inspect.getsource(helper).strip() for helper in helpers)
    return f"""from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests

MARKER_COLUMNS = {MARKER_COLUMNS!r}
MARKER_UNIVERSE_COLUMNS = {MARKER_UNIVERSE_COLUMNS!r}
MARKER_METHODS = {MARKER_METHODS!r}
_MARKER_NUMERIC_COLUMNS = {tuple(_MARKER_NUMERIC_COLUMNS)!r}
_MARKER_UNIT_INTERVAL_COLUMNS = {_MARKER_UNIT_INTERVAL_COLUMNS!r}

{sources}

def rank_marker_genes(adata):
    table, universe, _details = _rank_marker_evidence_impl(
        adata,
        groupby={groupby!r},
        method={method!r},
        source_kind={source_kind!r},
        layer_name={layer_name!r},
        n_genes={n_genes!r},
        tie_correct={tie_correct!r},
        max_output_rows={max_output_rows!r},
        max_working_memory_gib={max_working_memory_gib!r},
        np=np,
        pd=pd,
        sparse=sparse,
        sc=sc,
        ad=ad,
        multipletests=multipletests,
    )
    return table, universe
"""


def filter_marker_genes_code(
    *,
    min_log2_fold_change: float,
    min_fraction_in_group: float,
    max_fraction_reference: float,
    max_p_adjusted: float,
    expected_input_content_fingerprint: str,
    expected_universe_content_fingerprint: str,
) -> str:
    helpers = [
        _canonical_sha256,
        _validate_finite_threshold,
        _validate_canonical_marker_table,
        _validate_identifier_values,
        _validate_universe_frame,
        _frame_content_fingerprint,
        _marker_table_content_fingerprint,
        _universe_content_fingerprint,
        _validate_generated_marker_frames,
        _filter_marker_table_impl,
    ]
    sources = "\n\n".join(inspect.getsource(helper).strip() for helper in helpers)
    return f"""from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np
import pandas as pd

MARKER_COLUMNS = {MARKER_COLUMNS!r}
MARKER_UNIVERSE_COLUMNS = {MARKER_UNIVERSE_COLUMNS!r}
_MARKER_NUMERIC_COLUMNS = {tuple(_MARKER_NUMERIC_COLUMNS)!r}
_MARKER_UNIT_INTERVAL_COLUMNS = {_MARKER_UNIT_INTERVAL_COLUMNS!r}

{sources}

def filter_marker_genes(frame, universe):
    _validate_generated_marker_frames(frame, universe, np=np, pd=pd)
    if _marker_table_content_fingerprint(frame) != {expected_input_content_fingerprint!r}:
        raise ValueError("Marker table content does not match the source operation embedded in this code.")
    if _universe_content_fingerprint(universe) != {expected_universe_content_fingerprint!r}:
        raise ValueError("Tested-gene universe content does not match the source operation embedded in this code.")
    filtered, _ = _filter_marker_table_impl(
        frame,
        min_log2_fold_change={min_log2_fold_change!r},
        min_fraction_in_group={min_fraction_in_group!r},
        max_fraction_reference={max_fraction_reference!r},
        max_p_adjusted={max_p_adjusted!r},
        np=np,
        pd=pd,
    )
    return filtered, universe.copy()
"""


__all__ = [
    "FILTER_MARKER_GENES_NODE_ID",
    "MARKER_COLUMNS",
    "MARKER_EVIDENCE_SCHEMA_VERSION",
    "MARKER_GENES_NODE_ID",
    "MARKER_METHODS",
    "MARKER_UNIVERSE_COLUMNS",
    "filter_marker_genes_code",
    "filter_marker_table",
    "marker_genes_code",
    "marker_provenance_parameters",
    "marker_table_content_fingerprint",
    "marker_universe_content_fingerprint",
    "rank_marker_evidence",
    "validate_marker_artifact_pair",
]
