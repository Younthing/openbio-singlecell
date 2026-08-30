from __future__ import annotations

import hashlib
import inspect
import json
from collections.abc import Mapping, Sequence
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .analysis_reporting import collect_software_versions
from .analysis_utils import figure_to_png

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


POPULATION_CORRELATION_SCHEMA = "openbio-singlecell/population-centroid-correlation/v1"
POPULATION_CORRELATION_NODE_ID = "OpenBioSingleCellCellTypeCorrelation"

POPULATION_CORRELATION_REFERENCES = [
    {
        "citation": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19:15.",
        "doi": "10.1186/s13059-017-1382-0",
        "url": "https://doi.org/10.1186/s13059-017-1382-0",
        "kind": "software",
    },
    {
        "citation": "Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods. 2020;17:261-272.",
        "doi": "10.1038/s41592-019-0686-2",
        "url": "https://doi.org/10.1038/s41592-019-0686-2",
        "kind": "software",
    },
    {
        "citation": "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. Journal of Open Source Software. 2024;9:4371.",
        "doi": "10.21105/joss.04371",
        "url": "https://doi.org/10.21105/joss.04371",
        "kind": "software",
    },
]


def _strict_text(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    if not value or value != value.strip():
        raise ValueError(f"{label} must be nonempty and free of surrounding whitespace.")
    return value


def _missing_scalar(value: Any, pandas: Any, *, label: str) -> bool:
    missing = pandas.isna(value)
    if isinstance(missing, bool):
        return missing
    if getattr(missing, "ndim", 1) == 0:
        try:
            return bool(missing)
        except (TypeError, ValueError):
            pass
    raise ValueError(f"{label} values must be scalar.")


def _ordered_labels(series: Any, science: Any, *, key: str) -> tuple[list[str], list[str], list[str]]:
    values: list[str] = []
    observed: list[str] = []
    seen: set[str] = set()
    for value in series:
        if _missing_scalar(value, science.pd, label=f"Population column {key!r}"):
            raise ValueError(f"Population column {key!r} contains missing labels.")
        if not isinstance(value, str) or not value or value != value.strip():
            raise TypeError(
                f"Population column {key!r} must contain canonical nonempty string labels without coercion."
            )
        values.append(value)
        if value not in seen:
            seen.add(value)
            observed.append(value)
    unused: list[str] = []
    if isinstance(series.dtype, science.pd.CategoricalDtype):
        declared: list[str] = []
        for category in series.cat.categories:
            if not isinstance(category, str) or not category or category != category.strip():
                raise TypeError(f"Population categories for {key!r} must be canonical nonempty strings.")
            declared.append(category)
        observed = [category for category in declared if category in seen]
        unused = [category for category in declared if category not in seen]
    return values, observed, unused


def _axis_fingerprint(obs_names: Sequence[str], labels: Sequence[str]) -> str:
    payload = json.dumps(
        {"obs_names": list(obs_names), "population_labels": list(labels)},
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _representation_fingerprint(matrix: Any, numpy: Any) -> str:
    contiguous = numpy.ascontiguousarray(matrix)
    digest = hashlib.sha256()
    digest.update(json.dumps(list(contiguous.shape), separators=(",", ":")).encode("ascii"))
    digest.update(contiguous.dtype.str.encode("ascii"))
    digest.update(memoryview(contiguous).cast("B"))
    return digest.hexdigest()


def _representation_history(adata: AnnData, representation_key: str) -> list[dict[str, Any]]:
    root = adata.uns.get("openbio_singlecell")
    if not isinstance(root, Mapping):
        return []
    history = root.get("analysis_history")
    if not isinstance(history, Mapping):
        return []
    matches: list[dict[str, Any]] = []
    for history_id, entry in history.items():
        if not isinstance(entry, Mapping):
            continue
        parameters = entry.get("parameters")
        if not isinstance(parameters, Mapping):
            continue
        referenced = {
            value
            for key, value in parameters.items()
            if key in {"output_key", "representation_key", "use_rep", "embedding_key"} and isinstance(value, str)
        }
        if representation_key not in referenced:
            continue
        matches.append(
            {
                "history_id": str(history_id),
                "operation": str(entry.get("operation", "unknown")),
                "parameters": {str(key): value for key, value in parameters.items()},
            }
        )
    return matches[-5:]


def _validate_backend_dendrogram(
    info: Any,
    *,
    categories: Sequence[str],
    expected_matrix: Any,
    population_key: str,
    representation_key: str,
    correlation_method: str,
    linkage_method: str,
    science: Any,
) -> tuple[Any, list[str], Any]:
    if not isinstance(info, Mapping):
        raise RuntimeError("Scanpy dendrogram must return a result mapping when inplace=False.")
    required = {
        "linkage",
        "groupby",
        "use_rep",
        "cor_method",
        "linkage_method",
        "categories_ordered",
        "categories_idx_ordered",
        "correlation_matrix",
    }
    missing = sorted(required - set(info))
    if missing:
        raise RuntimeError(f"Scanpy dendrogram result is missing fields: {missing}.")
    matrix = science.np.asarray(info["correlation_matrix"], dtype=float)
    groups = len(categories)
    if matrix.shape != (groups, groups) or not bool(science.np.isfinite(matrix).all()):
        raise RuntimeError("Scanpy returned an invalid population-centroid correlation matrix.")
    if not bool(science.np.allclose(matrix, matrix.T, atol=1e-10, rtol=0.0)):
        raise RuntimeError("Scanpy population-centroid correlation matrix is not symmetric.")
    if not bool(science.np.allclose(science.np.diag(matrix), 1.0, atol=1e-10, rtol=0.0)):
        raise RuntimeError("Scanpy population-centroid correlation matrix must have a unit diagonal.")
    if bool(((matrix < -1.0 - 1e-12) | (matrix > 1.0 + 1e-12)).any()):
        raise RuntimeError("Scanpy population-centroid correlations must lie in [-1, 1].")
    if not bool(science.np.allclose(matrix, expected_matrix, atol=1e-10, rtol=1e-10)):
        raise RuntimeError("Scanpy correlation matrix disagrees with the independently recomputed centroid matrix.")
    ordered = list(info["categories_ordered"])
    indices = list(info["categories_idx_ordered"])
    if (
        not all(isinstance(value, str) for value in ordered)
        or len(ordered) != groups
        or set(ordered) != set(categories)
        or len(indices) != groups
        or sorted(indices) != list(range(groups))
        or ordered != [categories[index] for index in indices]
    ):
        raise RuntimeError("Scanpy returned an invalid dendrogram population order.")
    linkage = science.np.asarray(info["linkage"], dtype=float)
    if linkage.shape != (groups - 1, 4) or not bool(science.np.isfinite(linkage).all()):
        raise RuntimeError("Scanpy returned an invalid hierarchical-linkage matrix.")
    if (
        info["groupby"] != [population_key]
        or info["use_rep"] != representation_key
        or info["cor_method"] != correlation_method
        or info["linkage_method"] != linkage_method
    ):
        raise RuntimeError("Scanpy dendrogram metadata disagrees with requested method settings.")
    return matrix, ordered, linkage


def _render_correlation_png(
    matrix: Any,
    ordered: Sequence[str],
    *,
    color_map: str,
    show_numbers: bool,
    science: Any,
) -> bytes:
    import matplotlib

    if color_map not in matplotlib.colormaps:
        raise ValueError(f"Unknown Matplotlib colormap: {color_map!r}.")
    groups = len(ordered)
    side = min(18.0, max(6.0, 3.5 + groups * 0.42))
    figure = science.Figure(figsize=(side, side), constrained_layout=True)
    axis = figure.subplots()
    image = axis.imshow(matrix, cmap=color_map, vmin=-1.0, vmax=1.0, interpolation="nearest", aspect="equal")
    positions = list(range(groups))
    axis.set_xticks(positions, ordered, rotation=90)
    axis.set_yticks(positions, ordered)
    axis.set_xlabel("Population centroid")
    axis.set_ylabel("Population centroid")
    axis.set_title("Population centroid correlation")
    if show_numbers:
        for row in range(groups):
            for column in range(groups):
                value = float(matrix[row, column])
                color = "white" if abs(value) >= 0.55 else "black"
                axis.text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=7, color=color)
    colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    colorbar.set_label("Correlation")
    return figure_to_png(figure)


def _strongest_signed_pairs(table, *, direction: str, limit: int = 10):
    """Return deterministic signed extrema while retaining cutoff ties."""
    if direction == "positive":
        candidates = table.loc[table["correlation"] > 0].sort_values(
            ["correlation", "dendrogram_order_a", "dendrogram_order_b"],
            ascending=[False, True, True],
        )
        if len(candidates) > limit:
            cutoff = candidates.iloc[limit - 1]["correlation"]
            candidates = candidates.loc[candidates["correlation"] >= cutoff]
    elif direction == "negative":
        candidates = table.loc[table["correlation"] < 0].sort_values(
            ["correlation", "dendrogram_order_a", "dendrogram_order_b"],
            ascending=[True, True, True],
        )
        if len(candidates) > limit:
            cutoff = candidates.iloc[limit - 1]["correlation"]
            candidates = candidates.loc[candidates["correlation"] <= cutoff]
    else:
        raise ValueError(f"Unknown signed-extrema direction: {direction!r}.")
    return candidates[["population_a", "population_b", "correlation"]].to_dict(orient="records")


def analyze_population_centroid_correlation(
    adata: AnnData,
    *,
    population_key: str,
    representation_key: str,
    correlation_method: str,
    linkage_method: str,
    n_dimensions: int,
    annotation_status: str,
    color_map: str,
    show_numbers: bool,
    max_groups: int,
    max_output_rows: int,
) -> tuple[DataFrame, bytes, dict[str, Any]]:
    from . import dependencies

    science = dependencies.require_scientific_dependencies()
    population_key = _strict_text(population_key, label="Population key")
    representation_key = _strict_text(representation_key, label="Representation key")
    color_map = _strict_text(color_map, label="Color map")
    if correlation_method not in {"pearson", "spearman", "kendall"}:
        raise ValueError(f"Unsupported centroid correlation method: {correlation_method!r}.")
    if linkage_method not in {"complete", "average", "single", "weighted"}:
        raise ValueError(f"Unsupported centroid linkage method: {linkage_method!r}.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError(f"Unsupported annotation status: {annotation_status!r}.")
    if not isinstance(show_numbers, bool):
        raise TypeError("show_numbers must be a Boolean.")
    for value, name, minimum in (
        (n_dimensions, "n_dimensions", 0),
        (max_groups, "max_groups", 2),
        (max_output_rows, "max_output_rows", 1),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}.")
    try:
        obs_names = tuple(adata.obs_names)
        unique_obs = bool(adata.obs_names.is_unique)
        n_obs = int(adata.n_obs)
        n_vars = int(adata.n_vars)
    except AttributeError as error:
        raise TypeError("Population Centroid Correlation requires an AnnData input.") from error
    if not obs_names or not unique_obs:
        raise ValueError("Population Centroid Correlation requires nonempty unique obs_names.")
    if not all(isinstance(value, str) and value and value == value.strip() for value in obs_names):
        raise ValueError("Population Centroid Correlation requires canonical string obs_names.")
    if population_key not in adata.obs:
        raise ValueError(f"Population column not found in obs: {population_key!r}.")
    if representation_key not in adata.obsm:
        raise ValueError(f"Representation not found in obsm: {representation_key!r}.")
    labels, categories, unused_categories = _ordered_labels(
        adata.obs[population_key], science, key=population_key
    )
    if len(categories) < 2:
        raise ValueError("Population Centroid Correlation requires at least two observed populations.")
    if len(categories) > max_groups:
        raise ValueError(
            f"Population count {len(categories):,} exceeds max_groups={max_groups:,}; subset or raise the guard."
        )
    pairs = len(categories) * (len(categories) - 1) // 2
    if pairs > max_output_rows:
        raise ValueError(
            f"Correlation pair table would contain {pairs:,} rows, exceeding max_output_rows={max_output_rows:,}."
        )

    raw_representation = adata.obsm[representation_key]
    if science.sparse.issparse(raw_representation):
        raise TypeError("Population centroid representations must be a dense numeric matrix.")
    representation = science.np.asarray(raw_representation)
    if (
        representation.ndim != 2
        or representation.shape[0] != n_obs
        or representation.shape[1] < 2
        or science.np.issubdtype(representation.dtype, science.np.bool_)
        or science.np.issubdtype(representation.dtype, science.np.complexfloating)
        or not science.np.issubdtype(representation.dtype, science.np.number)
    ):
        raise ValueError("Population centroid representation must be a dense numeric cells-by-dimensions matrix.")
    if not bool(science.np.isfinite(representation).all()):
        raise ValueError("Population centroid representation contains NaN or infinity.")
    resolved_dimensions = representation.shape[1] if n_dimensions == 0 else n_dimensions
    if resolved_dimensions < 2 or resolved_dimensions > representation.shape[1]:
        raise ValueError(
            f"n_dimensions must resolve to 2..{representation.shape[1]}, observed {resolved_dimensions}."
        )
    selected = science.np.asarray(representation[:, :resolved_dimensions], dtype=float)
    informative_dimensions = int((science.np.ptp(selected, axis=0) > 0).sum())
    if informative_dimensions < 2:
        raise ValueError("Population centroid correlation requires at least two nonconstant dimensions.")

    counts = {category: labels.count(category) for category in categories}
    centroids = science.np.vstack([selected[science.np.asarray(labels) == category].mean(axis=0) for category in categories])
    expected_matrix = science.pd.DataFrame(centroids, index=categories).T.corr(method=correlation_method).to_numpy()
    if expected_matrix.shape != (len(categories), len(categories)) or not bool(
        science.np.isfinite(expected_matrix).all()
    ):
        raise ValueError(
            "Population centroid correlation is undefined; check for constant or collinear population centroids."
        )

    private_obs = science.pd.DataFrame(
        {population_key: science.pd.Categorical(labels, categories=categories, ordered=True)},
        index=list(obs_names),
    )
    work = science.ad.AnnData(science.np.zeros((n_obs, 1), dtype=science.np.float32), obs=private_obs)
    # Scanpy receives a private feature-only workspace so its dendrogram backend cannot alias the source representation.
    work.obsm[representation_key] = selected.copy()
    info = science.sc.tl.dendrogram(
        work,
        groupby=population_key,
        use_rep=representation_key,
        cor_method=correlation_method,
        linkage_method=linkage_method,
        optimal_ordering=True,
        inplace=False,
    )
    matrix, ordered, linkage = _validate_backend_dendrogram(
        info,
        categories=categories,
        expected_matrix=expected_matrix,
        population_key=population_key,
        representation_key=representation_key,
        correlation_method=correlation_method,
        linkage_method=linkage_method,
        science=science,
    )
    original_index = {category: index for index, category in enumerate(categories)}
    ordered_matrix = matrix[
        science.np.ix_([original_index[category] for category in ordered], [original_index[category] for category in ordered])
    ]
    records = []
    for left_index, left in enumerate(ordered):
        for right_index in range(left_index + 1, len(ordered)):
            right = ordered[right_index]
            records.append(
                {
                    "population_a": left,
                    "population_b": right,
                    "correlation": float(ordered_matrix[left_index, right_index]),
                    "dendrogram_order_a": left_index,
                    "dendrogram_order_b": right_index,
                    "cells_a": counts[left],
                    "cells_b": counts[right],
                    "correlation_method": correlation_method,
                    "representation_key": representation_key,
                    "n_dimensions": resolved_dimensions,
                }
            )
    columns = [
        "population_a",
        "population_b",
        "correlation",
        "dendrogram_order_a",
        "dendrogram_order_b",
        "cells_a",
        "cells_b",
        "correlation_method",
        "representation_key",
        "n_dimensions",
    ]
    table = science.pd.DataFrame.from_records(records, columns=columns)
    if len(table) != pairs or not bool(science.np.isfinite(table["correlation"].to_numpy(dtype=float)).all()):
        raise RuntimeError("Population centroid correlation pair table failed postconditions.")
    png = _render_correlation_png(
        ordered_matrix,
        ordered,
        color_map=color_map,
        show_numbers=show_numbers,
        science=science,
    )

    strongest_positive = _strongest_signed_pairs(table, direction="positive")
    strongest_negative = _strongest_signed_pairs(table, direction="negative")
    warnings = [
        "This is descriptive population-centroid similarity; it is not a p-value, formal population-difference test, lineage, or interaction analysis.",
        "Cells contribute directly to centroids, so Samples with more cells can have greater influence; Sample is not modeled as an inference unit.",
        "Technical batch is not adjusted here; interpretation depends on the declared upstream representation and its integration choices.",
        "Correlations across latent dimensions can change under representation rotation, scaling, subsetting, or retraining.",
    ]
    if annotation_status == "unknown":
        warnings.append("Population annotation status is unknown; labels must not be presented as curated identities.")
    elif annotation_status == "provisional":
        warnings.append("Population labels are provisional and require independent marker/evidence review.")
    if min(counts.values()) * 4 < max(counts.values()):
        warnings.append("Population cell counts differ by more than four-fold, so centroid support is strongly imbalanced.")
    if unused_categories:
        warnings.append("Unused declared population categories were excluded from the analysis.")
    history = _representation_history(adata, representation_key)
    results = (
        f"Centroids for {len(categories):,} observed populations were compared across the first "
        f"{resolved_dimensions:,} dimensions of {representation_key!r}, yielding {pairs:,} unique pairwise "
        f"{correlation_method} correlations and a {linkage_method}-linkage dendrogram."
    )
    parameters = {
        "population_key": population_key,
        "representation_key": representation_key,
        "correlation_method": correlation_method,
        "linkage_method": linkage_method,
        "n_dimensions": n_dimensions,
        "resolved_dimensions": resolved_dimensions,
        "annotation_status": annotation_status,
        "color_map": color_map,
        "show_numbers": show_numbers,
        "max_groups": max_groups,
        "max_output_rows": max_output_rows,
        "optimal_ordering": True,
    }
    summary = {
        "schema_version": POPULATION_CORRELATION_SCHEMA,
        "node_id": POPULATION_CORRELATION_NODE_ID,
        "status": "descriptive_population_similarity",
        "methods": (
            "For each population, arithmetic mean coordinates were computed in the selected representation. "
            f"Pairwise {correlation_method} correlations were calculated across representation dimensions; "
            f"Scanpy 1.12 dendrogram used 1-correlation distance, {linkage_method} linkage, and optimal ordering."
        ),
        "results": results,
        "key_results": {
            "input_cells": n_obs,
            "input_features": n_vars,
            "representation_shape": list(representation.shape),
            "resolved_dimensions": resolved_dimensions,
            "informative_dimensions": informative_dimensions,
            "population_count": len(categories),
            "pair_count": pairs,
            "population_cell_counts": counts,
            "declared_observed_order": categories,
            "dendrogram_order": ordered,
            "unused_declared_categories": unused_categories,
            "correlation_min": float(table["correlation"].min()),
            "correlation_max": float(table["correlation"].max()),
            "strongest_positive_pairs": strongest_positive,
            "strongest_negative_pairs": strongest_negative,
            "linkage_shape": list(linkage.shape),
            "observation_population_fingerprint_sha256": _axis_fingerprint(obs_names, labels),
            "representation_fingerprint_sha256": _representation_fingerprint(selected, science.np),
            "representation_history_matches": history,
            "annotation_status": annotation_status,
        },
        "parameters": parameters,
        "references": POPULATION_CORRELATION_REFERENCES,
        "software_versions": collect_software_versions(
            ("scanpy", "anndata", "numpy", "pandas", "scipy", "matplotlib")
        ),
        "warnings": warnings,
        "limitations": [
            "No statistical test or Sample-level uncertainty is estimated.",
            "Cell-count imbalance can weight population centroids toward Samples contributing more cells.",
            "Integrated representations can remove Technical-batch variation but can also alter biological geometry.",
            "Centroid correlation does not establish lineage, communication, causal similarity, or population identity.",
        ],
    }
    json.dumps(summary, allow_nan=False, ensure_ascii=False)
    return table, png, summary


def population_correlation_code(
    *,
    population_key: str,
    representation_key: str,
    correlation_method: str,
    linkage_method: str,
    n_dimensions: int,
    annotation_status: str,
    color_map: str,
    show_numbers: bool,
    max_groups: int,
    max_output_rows: int,
) -> str:
    helpers = (
        _strict_text,
        _missing_scalar,
        _ordered_labels,
        _axis_fingerprint,
        _representation_fingerprint,
        _representation_history,
        _validate_backend_dendrogram,
        _render_correlation_png,
        _strongest_signed_pairs,
    )
    helper_source = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    implementation = dedent(inspect.getsource(analyze_population_centroid_correlation)).strip()
    dependency_import = "    from . import dependencies\n\n    science = dependencies.require_scientific_dependencies()"
    if dependency_import not in implementation:
        raise RuntimeError("Population-correlation source extraction could not locate its dependency seam.")
    implementation = implementation.replace(
        dependency_import,
        "    science = _standalone_population_science_dependencies()",
        1,
    )
    return f'''from __future__ import annotations

import hashlib
import json
import platform
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from io import BytesIO
from types import SimpleNamespace
from typing import Any

POPULATION_CORRELATION_SCHEMA = {POPULATION_CORRELATION_SCHEMA!r}
POPULATION_CORRELATION_NODE_ID = {POPULATION_CORRELATION_NODE_ID!r}
POPULATION_CORRELATION_REFERENCES = {POPULATION_CORRELATION_REFERENCES!r}


def _standalone_population_science_dependencies():
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from matplotlib.figure import Figure
    from scipy import sparse

    return SimpleNamespace(ad=ad, np=np, pd=pd, sc=sc, sparse=sparse, Figure=Figure)


def figure_to_png(figure):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    FigureCanvasAgg(figure)
    buffer = BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    return buffer.getvalue()


def collect_software_versions(packages):
    """Collect the same strict version payload without importing OpenBio."""
    versions = {{
        "python": platform.python_version(),
        "openbio-singlecell": {PLUGIN_VERSION!r},
    }}
    for package in dict.fromkeys(packages):
        if package in versions:
            continue
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


{helper_source}


{implementation}


def population_centroid_correlation(adata):
    """Return (pair_table, png_bytes, strict_summary) without mutating adata."""
    return analyze_population_centroid_correlation(
        adata,
        population_key={population_key!r},
        representation_key={representation_key!r},
        correlation_method={correlation_method!r},
        linkage_method={linkage_method!r},
        n_dimensions={n_dimensions!r},
        annotation_status={annotation_status!r},
        color_map={color_map!r},
        show_numbers={show_numbers!r},
        max_groups={max_groups!r},
        max_output_rows={max_output_rows!r},
    )
'''


__all__ = [
    "POPULATION_CORRELATION_NODE_ID",
    "POPULATION_CORRELATION_REFERENCES",
    "POPULATION_CORRELATION_SCHEMA",
    "analyze_population_centroid_correlation",
    "population_correlation_code",
]
