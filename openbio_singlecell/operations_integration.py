from __future__ import annotations

import importlib
import inspect
import random
import re
import time
import warnings as python_warnings
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, make_table_result, matrix_totals_and_nonzero, validate_count_expression
from .artifact_codecs import read_anndata, write_table
from .artifact_envelope import result_metadata
from .expression_source import DynamicExpressionSource, ExpressionSource, ExpressionSourceSpec
from .expression_state import resolve_expression_state
from .graph_analysis import (
    assess_leiden_stability,
    graph_diagnostics,
    resolve_named_graph,
    run_leiden_partition,
    validate_graph_result_key,
    validate_leiden_settings,
    validate_random_seed,
)
from .operations_input import (
    ANNDATA_CODEC,
    ANNDATA_KIND,
    require_artifact_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
)
from .scvi_model import SCVI_GLOBAL_RNG_LOCK
from .worker_protocol import JSONValue, OperationContext, register_operation

if TYPE_CHECKING:
    from anndata import AnnData


INTEGRATION_SOFTWARE_PACKAGES = ("anndata", "numpy", "pandas", "scipy")
_LEIDEN_SWEEP_WORKLOAD_WARNING = 256
HARMONYPY_VERSION = "2.0.0"

HARMONY_REFERENCE = AnalysisReference(
    citation=(
        "Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell "
        "data with Harmony. Nature Methods. 2019;16:1289-1296."
    ),
    doi="10.1038/s41592-019-0619-0",
    url="https://doi.org/10.1038/s41592-019-0619-0",
    kind="method",
)
HARMONYPY_REFERENCE = AnalysisReference(
    citation="Slowikowski K et al. harmonypy method-author implementation and documentation.",
    url="https://github.com/slowkow/harmonypy",
    kind="software_documentation",
)
SCANPY_HARMONY_REFERENCE = AnalysisReference(
    citation="Scanpy external Harmony integration reference.",
    url="https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html",
    kind="software_documentation",
)
SCVI_REFERENCE = AnalysisReference(
    citation=(
        "Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell "
        "transcriptomics. Nature Methods. 2018;15:1053-1058."
    ),
    doi="10.1038/s41592-018-0229-2",
    url="https://doi.org/10.1038/s41592-018-0229-2",
    kind="method",
)
SCVI_TOOLS_REFERENCE = AnalysisReference(
    citation=(
        "Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell "
        "omics data. Nature Biotechnology. 2022;40:163-166."
    ),
    doi="10.1038/s41587-021-01206-w",
    url="https://doi.org/10.1038/s41587-021-01206-w",
    kind="software",
)
SCVI_API_REFERENCE = AnalysisReference(
    citation="scvi-tools SCVI public API and introductory workflow.",
    url="https://docs.scvi-tools.org/en/latest/api/reference/scvi.model.SCVI.html",
    kind="software_documentation",
)
INTEGRATION_BENCHMARK_REFERENCE = AnalysisReference(
    citation=(
        "Luecken MD, Buttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in "
        "single-cell genomics. Nature Methods. 2022;19:41-50."
    ),
    doi="10.1038/s41592-021-01336-8",
    url="https://doi.org/10.1038/s41592-021-01336-8",
    kind="practice",
)
LEIDEN_REFERENCE = AnalysisReference(
    citation="Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden. Scientific Reports. 2019;9:5233.",
    doi="10.1038/s41598-019-41695-z",
    url="https://doi.org/10.1038/s41598-019-41695-z",
    kind="method",
)
SCANPY_REFERENCE = AnalysisReference(
    citation="Wolf FA, Angerer P, Theis FJ. SCANPY. Genome Biology. 2018;19:15.",
    doi="10.1186/s13059-017-1382-0",
    url="https://doi.org/10.1186/s13059-017-1382-0",
    kind="software",
)
ARI_REFERENCE = AnalysisReference(
    citation="Hubert L, Arabie P. Comparing partitions. Journal of Classification. 1985;2:193-218.",
    doi="10.1007/BF01908075",
    url="https://doi.org/10.1007/BF01908075",
    kind="method",
)


def _require_in_memory_axes(adata: AnnData, *, operation: str, require_variables: bool = True) -> None:
    if bool(adata.isbacked):
        raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
    if int(adata.n_obs) < 2:
        raise ValueError(f"{operation} requires at least two cells.")
    if require_variables and int(adata.n_vars) < 1:
        raise ValueError(f"{operation} requires at least one feature.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    if require_variables and not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique feature identifiers.")


def _numeric_matrix(matrix: Any, *, description: str, shape: tuple[int, int]) -> Any:
    science = dependencies.require_scientific_dependencies()
    if tuple(matrix.shape) != shape:
        raise ValueError(f"{description} must have shape {shape}; received {tuple(matrix.shape)}.")
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    values = science.np.asarray(values)
    if values.dtype.kind not in "iuf":
        raise ValueError(f"{description} must be numeric.")
    if values.size and not bool(science.np.isfinite(values).all()):
        raise ValueError(f"{description} contains non-finite values.")
    return matrix


def _comma_separated_keys(value: str) -> list[str]:
    if not isinstance(value, str):
        raise TypeError("Observation-column lists must be strings.")
    keys = [item.strip() for item in value.split(",") if item.strip()]
    if len(keys) != len(set(keys)):
        raise ValueError("Observation-column lists cannot contain duplicate keys.")
    return keys


def _categorical_levels(adata: AnnData, keys: Sequence[str], *, role: str) -> dict[str, dict[str, int]]:
    levels: dict[str, dict[str, int]] = {}
    for key in keys:
        if key not in adata.obs:
            raise ValueError(f"{role} column not found in obs: {key!r}")
        series = adata.obs[key]
        if bool(series.isna().any()):
            raise ValueError(f"{role} column contains missing values: {key!r}")
        labels = series.astype(str)
        if bool(labels.str.strip().eq("").any()):
            raise ValueError(f"{role} column contains blank values: {key!r}")
        counts = labels.value_counts(sort=False)
        levels[key] = {str(label): int(count) for label, count in counts.items()}
    return levels


def _supports_keywords(callable_: Any, names: Sequence[str], *, description: str) -> None:
    try:
        signature = inspect.signature(callable_)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Cannot inspect the installed {description} interface.") from error
    accepts_kwargs = any(parameter.kind == inspect.Parameter.VAR_KEYWORD for parameter in signature.parameters.values())
    missing = [name for name in names if name not in signature.parameters and not accepts_kwargs]
    if missing:
        raise RuntimeError(f"Installed {description} is incompatible; missing parameters: {missing}.")


def _compatible_harmonypy_2_0(version: Any) -> bool:
    return isinstance(version, str) and re.fullmatch(r"2\.0\.\d+(?:\+[A-Za-z0-9.-]+)?", version) is not None


def _require_harmonypy_200() -> tuple[Any, str]:
    try:
        harmonypy = importlib.import_module("harmonypy")
    except (ImportError, OSError) as error:
        raise RuntimeError("Harmony Integration requires a compatible harmonypy 2.0.x installation.") from error
    module_version = getattr(harmonypy, "__version__", None)
    try:
        distribution_version = importlib_metadata.version("harmonypy")
    except importlib_metadata.PackageNotFoundError as error:
        raise RuntimeError("Harmony Integration cannot verify the installed harmonypy distribution.") from error
    if (
        module_version != distribution_version
        or not _compatible_harmonypy_2_0(module_version)
        or not _compatible_harmonypy_2_0(distribution_version)
    ):
        raise RuntimeError(
            "Harmony Integration requires agreeing module/distribution harmonypy 2.0.x identities; module reports "
            f"{module_version!r} and distribution reports "
            f"{distribution_version!r}."
        )
    return harmonypy, distribution_version


def _source_label(expression: ExpressionSource) -> str:
    return f"layer:{expression.layer_name}" if expression.kind == "layer" else expression.kind


def _validate_scvi_counts(
    adata: AnnData, expression: ExpressionSource
) -> tuple[Any, Any, Any, int, str, str | None, bool, int, list[str]]:
    science = dependencies.require_scientific_dependencies()
    matrix = _numeric_matrix(expression.matrix(adata), description="scVI count source", shape=tuple(adata.shape))
    state, evidence = resolve_expression_state(adata, expression)
    try:
        count_warnings = validate_count_expression(
            matrix,
            source_label="scVI count source",
            require_integers=False,
            require_positive=False,
        )
    except TypeError as error:
        raise ValueError("scVI count source must contain finite non-negative numeric values.") from error
    cell_totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
    gene_totals, _ = matrix_totals_and_nonzero(matrix, axis=0)
    cell_totals = science.np.asarray(cell_totals, dtype=float)
    gene_totals = science.np.asarray(gene_totals, dtype=float)
    zero_cells = int((cell_totals <= 0).sum())
    zero_genes = int((gene_totals <= 0).sum())
    if zero_cells:
        raise ValueError(f"scVI count source contains {zero_cells} cells with zero total counts.")
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    values = science.np.asarray(values)
    integer_like = bool(science.np.allclose(values, science.np.rint(values), rtol=0.0, atol=1e-8))
    nonzero = int(matrix.count_nonzero()) if science.sparse.issparse(matrix) else int(science.np.count_nonzero(matrix))
    return matrix, cell_totals, gene_totals, nonzero, state, evidence, integer_like, zero_genes, count_warnings


def _resolve_epochs(epochs: Mapping[str, object] | None) -> tuple[str, int | None]:
    if epochs is None:
        return "automatic", None
    if not isinstance(epochs, Mapping):
        raise TypeError("scVI epochs must be a DynamicCombo value.")
    mode = epochs.get("epochs")
    if mode == "automatic":
        return "automatic", None
    if mode != "fixed":
        raise ValueError(f"Unsupported scVI epochs mode: {mode!r}")
    value = epochs.get("max_epochs")
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("scVI fixed max_epochs must be a positive integer.")
    return "fixed", value


def _resolve_harmony_theta(theta: Mapping[str, object] | None) -> tuple[str, float | None]:
    if theta is None:
        return "automatic", None
    if not isinstance(theta, Mapping):
        raise TypeError("Harmony theta must be a DynamicCombo value.")
    mode = theta.get("theta")
    if mode == "automatic":
        return "automatic", None
    if mode != "custom":
        raise ValueError(f"Unsupported Harmony theta mode: {mode!r}")
    value = theta.get("theta_value")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("Harmony custom theta must be a finite non-negative number.")
    resolved = float(value)
    science = dependencies.require_scientific_dependencies()
    if not bool(science.np.isfinite(resolved)) or resolved < 0:
        raise ValueError("Harmony custom theta must be a finite non-negative number.")
    return "custom", resolved


def _history_values(value: Any, science: Any) -> Any:
    if hasattr(value, "to_numpy"):
        value = value.to_numpy()
    array = science.np.asarray(value, dtype=float).ravel()
    return array[science.np.isfinite(array)]


def _scvi_training_diagnostics(model: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    history = getattr(model, "history", None)
    metrics: dict[str, dict[str, float | int | None]] = {}
    lengths = []
    if isinstance(history, Mapping):
        for key in (
            "elbo_train",
            "elbo_validation",
            "reconstruction_loss_train",
            "reconstruction_loss_validation",
            "kl_local_train",
            "kl_local_validation",
        ):
            if key not in history:
                continue
            try:
                values = _history_values(history[key], science)
            except (TypeError, ValueError):
                continue
            lengths.append(int(values.size))
            if values.size:
                metrics[key] = {
                    "n": int(values.size),
                    "first": float(values[0]),
                    "last": float(values[-1]),
                    "best": float(values.min()),
                }
    return {
        "actual_epochs": max(lengths) if lengths else None,
        "metrics": metrics,
        "train_cells": len(getattr(model, "train_indices", ())) if hasattr(model, "train_indices") else None,
        "validation_cells": len(getattr(model, "validation_indices", ()))
        if hasattr(model, "validation_indices")
        else None,
        "test_cells": len(getattr(model, "test_indices", ())) if hasattr(model, "test_indices") else None,
        "device": str(getattr(model, "device", "unknown")),
    }


_LEIDEN_METRIC_COLUMNS = (
    "resolution",
    "key",
    "n_clusters",
    "min_cluster_size",
    "median_cluster_size",
    "max_cluster_size",
    "smallest_cluster_fraction",
    "largest_cluster_fraction",
    "singleton_clusters",
    "modularity",
    "stability_repeats",
    "stability_mean_ari",
    "stability_min_ari",
    "stability_max_ari",
    "adjacent_previous_resolution",
    "adjacent_resolution_ari",
)


def _parse_resolutions(value: str) -> list[float]:
    science = dependencies.require_scientific_dependencies()
    if not isinstance(value, str):
        raise TypeError("Leiden resolutions must be provided as a comma-separated string.")
    tokens = value.split(",")
    if not tokens or any(not token.strip() for token in tokens):
        raise ValueError("Leiden resolution list cannot be empty or contain empty tokens.")
    resolutions: list[float] = []
    for token in tokens:
        try:
            resolution = float(token.strip())
        except ValueError as error:
            raise ValueError(f"Invalid Leiden resolution: {token.strip()!r}.") from error
        if not bool(science.np.isfinite(resolution)):
            raise ValueError("Leiden resolutions must be finite.")
        if resolution < 0:
            raise ValueError("Leiden resolutions must be non-negative.")
        if resolution in resolutions:
            raise ValueError(f"duplicate Leiden resolution: {resolution!r}.")
        resolutions.append(resolution)
    return sorted(resolutions)


def _resolution_key(prefix: str, resolution: float) -> str:
    suffix = format(resolution, ".15g").replace("-", "m").replace(".", "_").replace("+", "p")
    return f"{prefix}_{suffix}"


def _leiden_resolution_sweep_code(parameters: Mapping[str, Any]) -> str:
    return dedent(
        f"""
        from collections.abc import Mapping

        import igraph
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import warnings
        from scipy import sparse
        from sklearn.metrics import adjusted_rand_score


        def leiden_resolution_sweep(adata):
            resolutions = {parameters["resolutions"]!r}
            keys = {parameters["keys"]!r}
            neighbors_key = {parameters["neighbors_key"]!r}
            n_iterations = {parameters["n_iterations"]!r}
            stability_repeats = {parameters["stability_repeats"]!r}
            random_seed = {parameters["random_seed"]!r}
            overwrite_existing = {parameters["overwrite_existing"]!r}
            if adata.isbacked or adata.n_obs < 2:
                raise ValueError("Leiden resolution sweep requires an in-memory AnnData with at least two observations.")
            if not adata.obs_names.is_unique:
                raise ValueError("Leiden resolution sweep requires unique observation identifiers.")
            if not isinstance(n_iterations, int) or isinstance(n_iterations, bool):
                raise TypeError("Leiden n_iterations must be an integer.")
            if not isinstance(stability_repeats, int) or isinstance(stability_repeats, bool):
                raise TypeError("Leiden sweep stability_repeats must be an integer.")
            if stability_repeats < 1:
                raise ValueError("Leiden sweep stability_repeats must be positive.")
            if not isinstance(random_seed, int) or isinstance(random_seed, bool):
                raise TypeError("Leiden random_seed must be an integer.")
            if random_seed < 0 or random_seed > 2**31 - stability_repeats:
                raise ValueError("Leiden random_seed range exceeds the supported seed boundary.")
            if stability_repeats == 1:
                warnings.warn("Leiden sweep stability was not assessed because stability_repeats=1.", UserWarning, stacklevel=2)
            if any(resolution == 0 for resolution in resolutions):
                warnings.warn("Leiden sweep includes resolution=0, which commonly produces very few communities.", UserWarning, stacklevel=2)
            if n_iterations == 0:
                warnings.warn("Leiden sweep n_iterations=0 requests no optimization iterations.", UserWarning, stacklevel=2)
            if len(resolutions) * stability_repeats > {_LEIDEN_SWEEP_WORKLOAD_WARNING!r}:
                warnings.warn("Leiden sweep requests more than 256 backend runs; execution may be resource intensive.", UserWarning, stacklevel=2)
            metadata = adata.uns.get(neighbors_key)
            if not isinstance(metadata, Mapping):
                raise ValueError("Named neighbor graph metadata is missing.")
            connectivities_key = metadata.get("connectivities_key")
            if not isinstance(connectivities_key, str) or not connectivities_key.strip():
                raise ValueError("Named neighbor graph has no valid connectivities pointer.")
            connectivities_key = connectivities_key.strip()
            if connectivities_key not in adata.obsp:
                raise ValueError("Named neighbor connectivity matrix is missing.")
            adjacency = adata.obsp[connectivities_key]
            if not sparse.issparse(adjacency):
                raise TypeError("Leiden connectivity matrix must be sparse.")
            if tuple(adjacency.shape) != (adata.n_obs, adata.n_obs):
                raise ValueError("Leiden connectivity matrix must have shape (n_obs, n_obs).")
            adjacency = adjacency.tocsr(copy=True)
            adjacency.sum_duplicates()
            adjacency.eliminate_zeros()
            values = np.asarray(adjacency.data)
            if values.size and values.dtype.kind not in "iuf":
                raise TypeError("Leiden connectivity weights must be numeric.")
            if values.size and not np.isfinite(values).all():
                raise ValueError("Leiden connectivity weights must be finite.")
            if values.size and (values < 0).any():
                raise ValueError("Leiden connectivity weights must be non-negative.")
            if not np.allclose(adjacency.diagonal(), 0.0, rtol=0.0, atol=1e-12):
                raise ValueError("Leiden connectivity matrix must have a zero diagonal.")
            delta = adjacency - adjacency.T
            delta.eliminate_zeros()
            if delta.nnz and np.max(np.abs(delta.data)) > 1e-8:
                raise ValueError("Leiden connectivity matrix must be symmetric and undirected.")
            adjacency = ((adjacency + adjacency.T) * 0.5).tocsr()
            adjacency.setdiag(0)
            adjacency.eliminate_zeros()
            if adjacency.nnz == 0:
                raise ValueError("Leiden connectivity graph contains no positive edges.")
            if any(key == neighbors_key for key in keys):
                raise ValueError(
                    "Leiden resolution output keys cannot equal the resolved neighbors_key because that uns key stores graph metadata."
                )
            collisions = [
                container
                for key in keys
                for container, present in ((f"obs[{{key!r}}]", key in adata.obs), (f"uns[{{key!r}}]", key in adata.uns))
                if present
            ]
            if collisions and not overwrite_existing:
                raise ValueError(f"Leiden resolution sweep outputs already exist: {{', '.join(collisions)}}")

            upper = sparse.triu(adjacency, k=1).tocoo()
            modularity_graph = igraph.Graph(
                n=adata.n_obs,
                edges=list(zip(upper.row.tolist(), upper.col.tolist(), strict=True)),
                directed=False,
            )
            modularity_graph.es["weight"] = [float(value) for value in upper.data]

            def validated_partition(target, key, resolution):
                labels = target.obs.get(key)
                if labels is None or not isinstance(labels.dtype, pd.CategoricalDtype):
                    raise RuntimeError("Leiden backend did not return categorical memberships.")
                if len(labels) != adata.n_obs or not labels.index.equals(adata.obs_names):
                    raise RuntimeError("Leiden memberships do not align to the graph observations.")
                if labels.isna().any():
                    raise RuntimeError("Leiden backend returned missing memberships.")
                observed_values = tuple(str(value) for value in labels.astype(object))
                observed = tuple(dict.fromkeys(observed_values))
                mapping = {{label: str(index) for index, label in enumerate(observed)}}
                membership = tuple(mapping[label] for label in observed_values)
                categories = tuple(str(index) for index in range(len(observed)))
                target.obs[key] = pd.Series(
                    pd.Categorical(membership, categories=categories),
                    index=target.obs_names,
                    name=key,
                )
                labels = target.obs[key]
                counts = labels.value_counts(sort=False)
                cluster_sizes = {{str(label): int(count) for label, count in counts.items()}}
                if sum(cluster_sizes.values()) != adata.n_obs or any(count <= 0 for count in cluster_sizes.values()):
                    raise RuntimeError("Leiden cluster sizes do not cover every observation exactly once.")
                result_metadata = target.uns.get(key)
                backend_value = result_metadata.get("modularity") if isinstance(result_metadata, Mapping) else None
                try:
                    backend_modularity = float(backend_value)
                except (TypeError, ValueError) as error:
                    raise RuntimeError("Leiden backend modularity is missing or invalid.") from error
                if not np.isfinite(backend_modularity):
                    raise RuntimeError("Leiden backend modularity is non-finite.")
                label_ids = {{label: index for index, label in enumerate(dict.fromkeys(membership))}}
                calculated_modularity = float(
                    modularity_graph.modularity(
                        [label_ids[label] for label in membership],
                        weights=modularity_graph.es["weight"],
                        resolution=resolution,
                        directed=False,
                    )
                )
                if not np.isfinite(calculated_modularity) or not np.isclose(
                    backend_modularity, calculated_modularity, rtol=1e-9, atol=1e-12
                ):
                    raise RuntimeError("Leiden backend modularity is inconsistent with the validated graph partition.")
                return membership, cluster_sizes, backend_modularity

            output = adata
            if overwrite_existing:
                for key in keys:
                    output.obs.drop(columns=[key], inplace=True, errors="ignore")
                    output.uns.pop(key, None)
            rows = []
            previous_resolution = None
            previous_membership = None
            for resolution, key in zip(resolutions, keys, strict=True):
                memberships = []
                cluster_sizes = None
                modularity = None
                for offset in range(stability_repeats):
                    run_key = key
                    if offset:
                        run_key = f"__openbio_leiden_stability_{{offset}}"
                        collision = 0
                        while run_key in output.obs or run_key in output.uns:
                            collision += 1
                            run_key = f"__openbio_leiden_stability_{{offset}}_{{collision}}"
                    try:
                        sc.tl.leiden(
                            output,
                            adjacency=adjacency,
                            resolution=resolution,
                            key_added=run_key,
                            random_state=random_seed + offset,
                            flavor="igraph",
                            n_iterations=n_iterations,
                            directed=False,
                            use_weights=True,
                            objective_function="modularity",
                        )
                        membership, run_sizes, run_modularity = validated_partition(output, run_key, resolution)
                    finally:
                        if offset:
                            output.obs.drop(columns=[run_key], inplace=True, errors="ignore")
                            output.uns.pop(run_key, None)
                    memberships.append(membership)
                    if offset == 0:
                        cluster_sizes = run_sizes
                        modularity = run_modularity
                pairwise = [
                    float(adjusted_rand_score(memberships[left], memberships[right]))
                    for left in range(len(memberships))
                    for right in range(left + 1, len(memberships))
                ]
                if pairwise and (not np.isfinite(pairwise).all() or any(not -1.0 <= value <= 1.0 for value in pairwise)):
                    raise RuntimeError("Leiden stability ARI diagnostics are invalid.")
                counts = list(cluster_sizes.values())
                adjacent_ari = (
                    None
                    if previous_membership is None
                    else float(adjusted_rand_score(previous_membership, memberships[0]))
                )
                row = {{
                    "resolution": resolution,
                    "key": key,
                    "n_clusters": len(cluster_sizes),
                    "min_cluster_size": min(counts),
                    "median_cluster_size": float(np.median(counts)),
                    "max_cluster_size": max(counts),
                    "smallest_cluster_fraction": float(min(counts) / adata.n_obs),
                    "largest_cluster_fraction": float(max(counts) / adata.n_obs),
                    "singleton_clusters": sum(count == 1 for count in counts),
                    "modularity": modularity,
                    "stability_repeats": stability_repeats,
                    "stability_mean_ari": float(sum(pairwise) / len(pairwise)) if pairwise else None,
                    "stability_min_ari": min(pairwise) if pairwise else None,
                    "stability_max_ari": max(pairwise) if pairwise else None,
                    "adjacent_previous_resolution": previous_resolution,
                    "adjacent_resolution_ari": adjacent_ari,
                }}
                rows.append(row)
                output.uns[key]["openbio_diagnostics"] = {{
                    "modularity": modularity,
                    "cluster_sizes": cluster_sizes,
                    "singleton_count": row["singleton_clusters"],
                    "stability": {{
                        "starts": stability_repeats,
                        "seeds": [random_seed + offset for offset in range(stability_repeats)],
                        "pairwise_ari": pairwise,
                        "mean_ari": row["stability_mean_ari"],
                        "min_ari": row["stability_min_ari"],
                        "max_ari": row["stability_max_ari"],
                    }},
                }}
                previous_resolution = resolution
                previous_membership = memberships[0]
            metrics = pd.DataFrame(rows, columns={list(_LEIDEN_METRIC_COLUMNS)!r})
            for column in {
            [
                "stability_mean_ari",
                "stability_min_ari",
                "stability_max_ari",
                "adjacent_previous_resolution",
                "adjacent_resolution_ari",
            ]!r
        }:
                metrics[column] = pd.array(metrics[column], dtype="Float64")
            return output, metrics
        """
    ).strip()


def _harmony_code(parameters: Mapping[str, Any]) -> str:
    return dedent(
        f"""
        def run_harmony_integration(adata):
            import inspect
            import re
            import warnings
            from collections.abc import Mapping
            from importlib import metadata as importlib_metadata

            import numpy as np
            import harmonypy

            technical_batch_keys = {parameters["technical_batch_keys"]!r}
            basis = {parameters["basis"]!r}
            adjusted_basis = {parameters["adjusted_basis"]!r}
            theta_mode = {parameters["theta_mode"]!r}
            theta_value = {parameters["theta"]!r}
            module_version = getattr(harmonypy, "__version__", None)
            try:
                distribution_version = importlib_metadata.version("harmonypy")
            except importlib_metadata.PackageNotFoundError as error:
                raise RuntimeError("Cannot verify the installed harmonypy distribution.") from error
            compatible = lambda value: isinstance(value, str) and re.fullmatch(
                r"2\\.0\\.\\d+(?:\\+[A-Za-z0-9.-]+)?", value
            ) is not None
            if module_version != distribution_version or not compatible(module_version):
                raise RuntimeError(
                    "Harmony requires agreeing module/distribution harmonypy 2.0.x identities; module reports "
                    f"{{module_version!r}} and distribution reports {{distribution_version!r}}."
                )
            if adata.isbacked or adata.n_obs < 2:
                raise ValueError("Harmony requires an in-memory nonempty AnnData with at least two cells.")
            if not adata.obs_names.is_unique:
                raise ValueError("Harmony requires unique observation identifiers.")
            if basis == adjusted_basis:
                raise ValueError("Harmony adjusted_basis must differ from basis so the original embedding is preserved.")
            if basis not in adata.obsm:
                raise ValueError(f"Harmony basis not found: {{basis!r}}")
            values = np.asarray(adata.obsm[basis])
            if values.ndim != 2 or values.shape[0] != adata.n_obs or values.shape[1] < 1:
                raise ValueError("Harmony basis must be n_obs by at least one component.")
            if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
                raise ValueError("Harmony basis must contain finite numeric values.")
            constant_components = np.flatnonzero(np.max(values, axis=0) == np.min(values, axis=0))
            if values.shape[1] == 1:
                warnings.warn(
                    "Harmony received a one-dimensional basis; correction is executable but geometric "
                    "interpretation is limited.",
                    UserWarning,
                    stacklevel=2,
                )
            if constant_components.size:
                warnings.warn(
                    f"Harmony basis contains {{int(constant_components.size)}} zero-variance component(s) at "
                    f"zero-based indices {{constant_components.tolist()}}; the backend result requires careful review.",
                    UserWarning,
                    stacklevel=2,
                )
            for key in technical_batch_keys:
                if key not in adata.obs or adata.obs[key].isna().any():
                    raise ValueError(f"Invalid Harmony Technical batch column: {{key!r}}")
                labels = adata.obs[key].astype(str)
                if labels.str.strip().eq("").any():
                    raise ValueError(f"Harmony Technical batch contains blank values: {{key!r}}")
                counts = labels.value_counts(sort=False)
                if labels.nunique() < 2:
                    warnings.warn(
                        f"Harmony correction covariate {{key!r}} has only one observed level and provides no "
                        "correction information; preserving the explicit expert choice.",
                        UserWarning,
                        stacklevel=2,
                    )
                if (counts < 2).any():
                    warnings.warn(
                        f"Harmony Technical batch {{key!r}} contains singleton levels and may be unstable.",
                        UserWarning,
                        stacklevel=2,
                    )
                small = {{str(level): int(count) for level, count in counts.items() if count < 10}}
                if small:
                    warnings.warn(
                        f"Small Harmony Technical-batch levels may be unstable for {{key!r}}: {{small}}.",
                        UserWarning,
                        stacklevel=2,
                    )
            output = adata
            if adjusted_basis in output.obsm and not {parameters["overwrite_existing"]!r}:
                raise ValueError(f"Harmony adjusted basis already exists: {{adjusted_basis!r}}")
            existing_metadata = output.uns.get("harmony")
            if existing_metadata is not None and not isinstance(existing_metadata, Mapping):
                raise ValueError("AnnData uns['harmony'] must be a mapping before Harmony metadata can be stored.")
            if (
                isinstance(existing_metadata, Mapping)
                and adjusted_basis in existing_metadata
                and not {parameters["overwrite_existing"]!r}
            ):
                raise ValueError(f"Harmony metadata already exists for adjusted basis: {{adjusted_basis!r}}")
            required = {parameters["required_backend_keywords"]!r}
            signature = inspect.signature(harmonypy.run_harmony)
            has_kwargs = any(item.kind == inspect.Parameter.VAR_KEYWORD for item in signature.parameters.values())
            missing = [name for name in required if name not in signature.parameters and not has_kwargs]
            if missing:
                raise RuntimeError(f"Installed harmonypy is incompatible; missing parameters: {{missing}}.")
            n_clusters = {parameters["n_clusters"]!r}
            if not 1 <= n_clusters <= adata.n_obs:
                raise ValueError("Harmony n_clusters must be between one and the number of cells, inclusive.")
            backend_basis = np.asarray(values, dtype=float).copy()
            run_kwargs = {{
                "lamb": {parameters["ridge_penalty"]!r},
                "sigma": np.full(n_clusters, {parameters["sigma"]!r}, dtype=float),
                "nclust": n_clusters,
                "tau": {parameters["tau"]!r},
                "block_size": 0.05,
                "max_iter_harmony": {parameters["max_iter_harmony"]!r},
                "max_iter_kmeans": {parameters["max_iter_kmeans"]!r},
                "epsilon_cluster": 0.001,
                "epsilon_harmony": 0.01,
                "verbose": False,
                "random_state": {parameters["random_seed"]!r},
                "ncores": 1,
            }}
            if theta_mode == "custom":
                run_kwargs["theta"] = theta_value
            elif theta_mode != "automatic" or theta_value is not None:
                raise RuntimeError("Generated Harmony theta configuration is invalid.")
            result = harmonypy.run_harmony(backend_basis, output.obs, technical_batch_keys, **run_kwargs)
            if not hasattr(result, "Z_corr"):
                raise RuntimeError("harmonypy result does not expose corrected coordinates in Z_corr.")
            raw_corrected = np.asarray(result.Z_corr)
            expected_shape = values.shape
            if raw_corrected.ndim != 2 or raw_corrected.shape != expected_shape:
                raise RuntimeError(
                    f"harmonypy 2.0.x Z_corr must preserve (n_obs, n_components)={{expected_shape}}; "
                    f"received {{raw_corrected.shape}}."
                )
            if raw_corrected.dtype.kind not in "iuf" or not np.isfinite(raw_corrected).all():
                raise RuntimeError("harmonypy Z_corr must contain finite real numeric coordinates.")
            corrected = np.asarray(raw_corrected, dtype=float).copy()
            output.obsm[adjusted_basis] = corrected
            harmony_metadata = dict(existing_metadata or {{}})
            harmony_metadata[adjusted_basis] = {{
                "technical_batch_keys": list(technical_batch_keys),
                "basis": basis,
                "theta_mode": theta_mode,
                "theta": theta_value,
                "random_seed": {parameters["random_seed"]!r},
            }}
            output.uns["harmony"] = harmony_metadata
            return output
        """
    ).strip()


def _scvi_code(parameters: Mapping[str, Any]) -> str:
    source_layer = parameters.get("counts_layer")
    check_val_every_n_epoch = parameters.get(
        "check_val_every_n_epoch",
        None if parameters.get("train_size") == 1 else 1,
    )
    return dedent(
        f"""
        def run_scvi_integration(adata):
            import inspect
            import random
            import warnings
            from collections.abc import Mapping

            import numpy as np
            import scvi
            import torch
            from scipy import sparse

            def resolve_source_state():
                x_state = "unknown"
                x_evidence = None
                layer_states = {{}}
                metadata = adata.uns.get("openbio_singlecell")
                history = metadata.get("analysis_history") if isinstance(metadata, Mapping) else None
                entries = history.values() if isinstance(history, Mapping) else ()
                for entry in entries:
                    if not isinstance(entry, Mapping):
                        continue
                    operation = entry.get("operation")
                    entry_parameters = entry.get("parameters")
                    entry_parameters = entry_parameters if isinstance(entry_parameters, Mapping) else {{}}
                    if operation == "snapshot_expression":
                        layer_states["counts"] = ("counts", "OpenBio Snapshot Expression history")
                        if entry_parameters.get("source") == "X":
                            x_state, x_evidence = "counts", "OpenBio Snapshot Expression history"
                    elif operation == "normalize_to_layer":
                        layer_name = entry_parameters.get("output_layer")
                        transform = entry_parameters.get("transform")
                        if isinstance(layer_name, str):
                            state = "logged" if transform == "log1p" else "transformed" if transform == "sqrt" else "normalized"
                            layer_states[layer_name] = (state, "OpenBio Normalize To Layer history")
                    elif operation == "pearson_residuals_to_layer":
                        layer_name = entry_parameters.get("output_layer")
                        if isinstance(layer_name, str):
                            layer_states[layer_name] = ("pearson_residuals", "OpenBio Pearson Residuals history")
                    elif operation == "scale_to_layer":
                        layer_name = entry_parameters.get("output_layer")
                        if isinstance(layer_name, str):
                            layer_states[layer_name] = ("scaled", "OpenBio Scale history")

                    if operation == "normalize_total":
                        x_state, x_evidence = "normalized", "OpenBio Normalize Total history"
                    elif operation == "log1p":
                        x_state, x_evidence = "logged", "OpenBio Log1p history"

                if source_kind == "layer":
                    return layer_states.get(source_layer, ("unknown", None))
                if x_state == "unknown" and isinstance(adata.uns.get("log1p"), Mapping):
                    return "logged", "AnnData uns['log1p'] marker"
                return x_state, x_evidence

            source_kind = {parameters["source"]!r}
            source_layer = {source_layer!r}
            technical_batch_key = {parameters["technical_batch_key"]!r}
            categorical_keys = {parameters["categorical_covariates"]!r}
            continuous_keys = {parameters["continuous_covariates"]!r}
            size_factor_key = {parameters["size_factor_key"]!r}
            output_key = {parameters["output_key"]!r}
            if adata.isbacked or adata.n_obs < 2 or adata.n_vars < 1:
                raise ValueError("scVI requires an in-memory nonempty AnnData with at least two cells.")
            if not adata.obs_names.is_unique or not adata.var_names.is_unique:
                raise ValueError("scVI requires unique observation and feature identifiers.")
            if source_kind == "layer" and source_layer not in adata.layers:
                raise ValueError(f"scVI count-source layer not found: {{source_layer!r}}")
            matrix = adata.layers[source_layer] if source_kind == "layer" else adata.X
            if tuple(matrix.shape) != tuple(adata.shape):
                raise ValueError("scVI count source must align to AnnData.")
            source_state, source_evidence = resolve_source_state()
            if source_state == "unknown":
                warnings.warn(
                    "The selected count-intended source has no recorded raw-count provenance; verify whether it "
                    "contains unnormalized UMI counts or another explicitly intended non-negative representation.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            elif source_state != "counts":
                warnings.warn(
                    f"The explicitly selected scVI source has recorded expression state {{source_state!r}} "
                    f"({{source_evidence or 'no evidence text'}}). The backend can train on finite non-negative "
                    "values, but count-likelihood interpretation may not apply.",
                    UserWarning,
                    stacklevel=2,
                )
            values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
            values = np.asarray(values)
            if values.dtype.kind not in "iuf" or not np.isfinite(values).all() or (values < 0).any():
                raise ValueError("scVI requires finite non-negative numeric values.")
            if not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-8):
                warnings.warn(
                    "scVI count source contains non-integer values and is not an unnormalized count matrix. "
                    "Count totals should be interpreted as expression sums.",
                    UserWarning,
                    stacklevel=2,
                )
            cell_totals = np.asarray(matrix.sum(axis=1)).ravel()
            gene_totals = np.asarray(matrix.sum(axis=0)).ravel()
            zero_cells = int((cell_totals <= 0).sum())
            zero_genes = int((gene_totals <= 0).sum())
            if zero_cells:
                raise ValueError(f"scVI count source contains {{zero_cells}} cells with zero total counts.")
            if zero_genes:
                warnings.warn(
                    f"The selected scVI source contains {{zero_genes}} all-zero gene(s); the audited backend can "
                    "retain them, but they provide no expression information.",
                    UserWarning,
                    stacklevel=2,
                )
            if {parameters["n_latent"]!r} >= min(adata.n_obs, adata.n_vars):
                warnings.warn(
                    f"scVI n_latent={parameters["n_latent"]!r} is at least one input-axis dimension; this "
                    "overcomplete architecture is executable but may be inefficient or weakly identified.",
                    UserWarning,
                    stacklevel=2,
                )
            if {parameters["n_layers"]!r} > 20:
                warnings.warn(
                    f"scVI n_layers={parameters["n_layers"]!r} requests an unusually deep architecture; resource "
                    "use and optimization stability require review.",
                    UserWarning,
                    stacklevel=2,
                )
            if {parameters["dropout_rate"]!r} == 1:
                warnings.warn(
                    "scVI dropout_rate=1 drops every affected activation; the backend accepts this endpoint, but "
                    "the resulting architecture is degenerate and requires explicit review.",
                    UserWarning,
                    stacklevel=2,
                )
            if {parameters["train_size"]!r} < 0.5:
                warnings.warn(
                    f"scVI train_size={parameters["train_size"]!r} uses fewer than half of cells for training; "
                    "this is executable but may reduce model stability.",
                    UserWarning,
                    stacklevel=2,
                )
            if {parameters["train_size"]!r} == 1:
                warnings.warn(
                    "scVI train_size=1 uses every cell for training and leaves no validation holdout; validation "
                    "metrics and early-stopping interpretation may be unavailable.",
                    UserWarning,
                    stacklevel=2,
                )
            categorical_all = [technical_batch_key, *categorical_keys]
            for key in categorical_all:
                if key not in adata.obs or adata.obs[key].isna().any():
                    raise ValueError(f"Invalid scVI categorical nuisance column: {{key!r}}")
                labels = adata.obs[key].astype(str)
                if labels.str.strip().eq("").any():
                    raise ValueError(f"scVI nuisance column contains blank values: {{key!r}}")
                if labels.nunique() < 2:
                    warnings.warn(
                        f"scVI nuisance column {{key!r}} has only one observed level and provides no correction "
                        "information; preserving the explicit expert choice.",
                        UserWarning,
                        stacklevel=2,
                    )
                small = {{str(level): int(count) for level, count in labels.value_counts(sort=False).items() if count < 10}}
                if small:
                    warnings.warn(
                        f"Small Technical-batch or categorical-covariate levels may be unstable for {{key!r}}: {{small}}.",
                        RuntimeWarning,
                        stacklevel=2,
                    )
            for key in continuous_keys:
                if key not in adata.obs:
                    raise ValueError(f"scVI continuous nuisance column not found in obs: {{key!r}}")
                numeric = np.asarray(adata.obs[key], dtype=float)
                if not np.isfinite(numeric).all():
                    raise ValueError(f"Invalid non-finite scVI continuous nuisance column: {{key!r}}")
                if np.ptp(numeric) <= 0:
                    warnings.warn(
                        f"scVI continuous nuisance column {{key!r}} is constant and provides no correction "
                        "information; preserving the explicit expert choice.",
                        UserWarning,
                        stacklevel=2,
                    )
            if size_factor_key:
                if size_factor_key not in adata.obs:
                    raise ValueError(f"scVI size-factor column not found in obs: {{size_factor_key!r}}")
                factors = np.asarray(adata.obs[size_factor_key], dtype=float)
                if not np.isfinite(factors).all() or (factors <= 0).any():
                    raise ValueError("scVI size factors must be finite, positive, and on the linear scale.")
            if output_key in adata.obsm and not {parameters["overwrite_existing"]!r}:
                raise ValueError(f"scVI output key already exists: {{output_key!r}}")
            setup_required = ("layer", "batch_key", "size_factor_key", "categorical_covariate_keys", "continuous_covariate_keys")
            setup_signature = inspect.signature(scvi.model.SCVI.setup_anndata)
            setup_kwargs = any(item.kind == inspect.Parameter.VAR_KEYWORD for item in setup_signature.parameters.values())
            missing = [name for name in setup_required if name not in setup_signature.parameters and not setup_kwargs]
            if missing:
                raise RuntimeError(f"Installed scvi-tools setup interface is incompatible: {{missing}}")
            constructor_required = ("n_layers", "n_latent", "gene_likelihood", "dispersion", "dropout_rate")
            constructor_signature = inspect.signature(scvi.model.SCVI)
            constructor_kwargs = any(
                item.kind == inspect.Parameter.VAR_KEYWORD for item in constructor_signature.parameters.values()
            )
            missing = [
                name for name in constructor_required if name not in constructor_signature.parameters and not constructor_kwargs
            ]
            if missing:
                raise RuntimeError(f"Installed scvi-tools model interface is incompatible: {{missing}}")
            # scvi-tools setup mutates registration state and the saved model must own that exact training view.
            training_adata = adata.copy()
            previous_seed = getattr(scvi.settings, "seed", None)
            python_state = random.getstate()
            numpy_state = np.random.get_state()
            torch_state = torch.random.get_rng_state()
            cuda_states = torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None
            mps_module = getattr(torch, "mps", None)
            mps_backend = getattr(getattr(torch, "backends", None), "mps", None)
            mps_available = (
                mps_module is not None
                and callable(getattr(mps_module, "get_rng_state", None))
                and callable(getattr(mps_module, "set_rng_state", None))
                and mps_backend is not None
                and callable(getattr(mps_backend, "is_available", None))
                and mps_backend.is_available()
            )
            mps_state = mps_module.get_rng_state() if mps_available else None
            try:
                with torch.inference_mode(False):
                    scvi.settings.seed = {parameters["random_seed"]!r}
                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore",
                            message=r".*does not contain unnormalized count data.*",
                            category=UserWarning,
                        )
                        scvi.model.SCVI.setup_anndata(
                            training_adata,
                            layer=source_layer if source_kind == "layer" else None,
                            batch_key=technical_batch_key,
                            size_factor_key=size_factor_key or None,
                            categorical_covariate_keys=categorical_keys or None,
                            continuous_covariate_keys=continuous_keys or None,
                        )
                    model = scvi.model.SCVI(
                        training_adata,
                        n_layers={parameters["n_layers"]!r},
                        n_latent={parameters["n_latent"]!r},
                        gene_likelihood={parameters["gene_likelihood"]!r},
                        dispersion={parameters["dispersion"]!r},
                        dropout_rate={parameters["dropout_rate"]!r},
                    )
                    train_required = (
                        "max_epochs",
                        "accelerator",
                        "devices",
                        "train_size",
                        "validation_size",
                        "batch_size",
                        "early_stopping",
                    )
                    train_signature = inspect.signature(model.train)
                    train_kwargs = any(
                        item.kind == inspect.Parameter.VAR_KEYWORD for item in train_signature.parameters.values()
                    )
                    missing = [
                        name for name in train_required if name not in train_signature.parameters and not train_kwargs
                    ]
                    if missing:
                        raise RuntimeError(f"Installed scvi-tools training interface is incompatible: {{missing}}")
                    model.train(
                        max_epochs={parameters["max_epochs"]!r},
                        accelerator={parameters["accelerator"]!r},
                        devices=1,
                        train_size={parameters["train_size"]!r},
                        validation_size=None,
                        batch_size={parameters["batch_size"]!r},
                        early_stopping={parameters["early_stopping"]!r},
                        early_stopping_monitor="elbo_validation",
                        early_stopping_patience=10,
                        early_stopping_min_delta=0.0,
                        check_val_every_n_epoch={check_val_every_n_epoch!r},
                    )
                    latent = np.asarray(model.get_latent_representation(), dtype=float)
            finally:
                scvi.settings.seed = previous_seed
                random.setstate(python_state)
                np.random.set_state(numpy_state)
                torch.random.set_rng_state(torch_state)
                if cuda_states is not None:
                    torch.cuda.set_rng_state_all(cuda_states)
                if mps_state is not None:
                    mps_module.set_rng_state(mps_state)
            if not bool(getattr(model, "is_trained", False)):
                raise RuntimeError("scVI backend did not report a trained model after training.")
            if latent.shape != (adata.n_obs, {parameters["n_latent"]!r}) or not np.isfinite(latent).all():
                raise RuntimeError("scVI returned an invalid latent representation.")
            constant_latent = np.flatnonzero(np.max(latent, axis=0) == np.min(latent, axis=0))
            if constant_latent.size:
                warnings.warn(
                    f"scVI returned {{int(constant_latent.size)}} constant latent dimension(s) at zero-based "
                    f"indices {{constant_latent.tolist()}}; downstream geometry may be degenerate.",
                    UserWarning,
                    stacklevel=2,
                )
            output = adata
            output.obsm[output_key] = latent.copy()
            return output, model
        """
    ).strip()


class OpenBioSingleCellHarmonyIntegration:
    @classmethod
    def execute(
        cls,
        adata: AnnData,
        technical_batch_keys: str = "batch",
        basis: str = "X_pca",
        theta: Mapping[str, object] | None = None,
        ridge_penalty: float = -1.0,
        sigma: float = 0.1,
        n_clusters: int = 0,
        tau: float = 0.0,
        adjusted_basis: str = "X_pca_harmony",
        overwrite_existing: bool = False,
        max_iter_harmony: int = 10,
        max_iter_kmeans: int = 4,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        _require_in_memory_axes(adata, operation="Harmony Integration", require_variables=False)
        keys = _comma_separated_keys(technical_batch_keys)
        if not keys:
            raise ValueError("Harmony requires at least one Technical batch key.")
        level_counts = _categorical_levels(adata, keys, role="Harmony Technical batch")
        singleton_levels = {
            key: {level: count for level, count in counts.items() if count < 2} for key, counts in level_counts.items()
        }
        singleton_levels = {key: counts for key, counts in singleton_levels.items() if counts}
        basis = basis.strip()
        adjusted_basis = adjusted_basis.strip()
        if not basis or not adjusted_basis:
            raise ValueError("Harmony basis and adjusted_basis cannot be empty.")
        if basis == adjusted_basis:
            raise ValueError("Harmony adjusted_basis must differ from basis so the original embedding is preserved.")
        if basis not in adata.obsm:
            raise ValueError(f"Harmony basis not found in obsm: {basis!r}")
        if science.sparse.issparse(adata.obsm[basis]):
            raise ValueError("Harmony PCA basis must be a dense numeric matrix.")
        basis_values = science.np.asarray(adata.obsm[basis])
        if basis_values.ndim != 2:
            raise ValueError("Harmony PCA basis must be a two-dimensional cell-by-component matrix.")
        _numeric_matrix(
            basis_values,
            description="Harmony PCA basis",
            shape=(int(adata.n_obs), int(basis_values.shape[1])),
        )
        if basis_values.shape[1] < 1:
            raise ValueError("Harmony PCA basis must contain at least one component.")
        constant_component_indices = science.np.flatnonzero(
            science.np.max(basis_values, axis=0) == science.np.min(basis_values, axis=0)
        )
        if adjusted_basis in adata.obsm and not overwrite_existing:
            raise ValueError(f"Harmony adjusted basis already exists: {adjusted_basis!r}")
        existing_metadata = adata.uns.get("harmony")
        if existing_metadata is not None and not isinstance(existing_metadata, Mapping):
            raise ValueError("AnnData uns['harmony'] must be a mapping before Harmony metadata can be stored.")
        if isinstance(existing_metadata, Mapping) and adjusted_basis in existing_metadata and not overwrite_existing:
            raise ValueError(f"Harmony metadata already exists for adjusted basis: {adjusted_basis!r}")
        theta_mode, resolved_theta = _resolve_harmony_theta(theta)
        if isinstance(tau, bool) or not isinstance(tau, (int, float)) or not science.np.isfinite(tau) or tau < 0:
            raise ValueError("Harmony tau must be finite and non-negative.")
        if (
            isinstance(ridge_penalty, bool)
            or not isinstance(ridge_penalty, (int, float))
            or not science.np.isfinite(ridge_penalty)
            or (ridge_penalty != -1 and ridge_penalty <= 0)
        ):
            raise ValueError("Harmony ridge_penalty must be -1 for automatic estimation or a positive value.")
        if (
            isinstance(sigma, bool)
            or not isinstance(sigma, (int, float))
            or not science.np.isfinite(sigma)
            or sigma <= 0
        ):
            raise ValueError("Harmony sigma must be finite and positive.")
        for name, value in (
            ("n_clusters", n_clusters),
            ("max_iter_harmony", max_iter_harmony),
            ("max_iter_kmeans", max_iter_kmeans),
            ("random_seed", random_seed),
        ):
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"Harmony {name} must be an integer.")
        if max_iter_harmony < 1 or max_iter_kmeans < 1 or not 0 <= random_seed <= 2**31 - 1:
            raise ValueError("Harmony iteration limits must be positive and random_seed in [0, 2**31 - 1].")
        if not isinstance(overwrite_existing, bool):
            raise TypeError("Harmony overwrite_existing must be a boolean.")
        if n_clusters < 0:
            raise ValueError("Harmony n_clusters must be zero for automatic selection or a positive integer.")
        resolved_clusters = max(1, min(round(int(adata.n_obs) / 30), 100)) if n_clusters == 0 else n_clusters
        if resolved_clusters > int(adata.n_obs):
            raise ValueError("Harmony n_clusters cannot exceed the number of cells.")

        harmonypy, harmonypy_version = _require_harmonypy_200()
        required_keywords = (
            "lamb",
            "sigma",
            "nclust",
            "tau",
            "block_size",
            "max_iter_harmony",
            "max_iter_kmeans",
            "epsilon_cluster",
            "epsilon_harmony",
            "verbose",
            "random_state",
            "ncores",
        )
        if theta_mode == "custom":
            required_keywords = (*required_keywords, "theta")
        _supports_keywords(harmonypy.run_harmony, required_keywords, description="harmonypy.run_harmony")
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata
        # harmonypy may mutate its array; keep only the basis-sized scratch/check buffers, never a whole AnnData copy.
        backend_basis = science.np.asarray(basis_values, dtype=float).copy()
        sigma_vector = science.np.full(resolved_clusters, float(sigma), dtype=float)
        run_kwargs = {
            "lamb": None if ridge_penalty == -1 else float(ridge_penalty),
            "sigma": sigma_vector,
            "nclust": resolved_clusters,
            "tau": float(tau),
            "block_size": 0.05,
            "max_iter_harmony": max_iter_harmony,
            "max_iter_kmeans": max_iter_kmeans,
            "epsilon_cluster": 1e-3,
            "epsilon_harmony": 1e-2,
            "verbose": False,
            "random_state": random_seed,
            "ncores": 1,
        }
        if theta_mode == "custom":
            run_kwargs["theta"] = resolved_theta
        result = harmonypy.run_harmony(backend_basis, output.obs, keys, **run_kwargs)
        if not hasattr(result, "Z_corr"):
            raise RuntimeError("harmonypy result does not expose corrected coordinates in Z_corr.")
        raw_corrected = science.np.asarray(result.Z_corr)
        expected_shape = tuple(int(value) for value in basis_values.shape)
        if raw_corrected.ndim != 2 or tuple(raw_corrected.shape) != expected_shape:
            raise RuntimeError(
                f"harmonypy {harmonypy_version} Z_corr must preserve "
                f"(n_obs, n_components)={expected_shape}; received {tuple(raw_corrected.shape)}."
            )
        if raw_corrected.dtype.kind not in "iuf" or not bool(science.np.isfinite(raw_corrected).all()):
            raise RuntimeError("harmonypy Z_corr must contain finite real numeric coordinates.")
        corrected = science.np.asarray(raw_corrected, dtype=float).copy()
        orientation = "n_obs_by_components_preserved"
        output.obsm[adjusted_basis] = corrected
        objective = science.np.asarray(getattr(result, "objective_harmony", ()), dtype=float).ravel()
        finite_objective = objective[science.np.isfinite(objective)]
        rounds = int(finite_objective.size)
        displacement = science.np.linalg.norm(corrected - basis_values, axis=1)
        parameters = {
            "technical_batch_keys": keys,
            "basis": basis,
            "adjusted_basis": adjusted_basis,
            "theta_mode": theta_mode,
            "theta": resolved_theta,
            "ridge_penalty": None if ridge_penalty == -1 else float(ridge_penalty),
            "sigma": float(sigma),
            "sigma_vector": sigma_vector.tolist(),
            "n_clusters": resolved_clusters,
            "tau": float(tau),
            "block_size": 0.05,
            "max_iter_harmony": max_iter_harmony,
            "max_iter_kmeans": max_iter_kmeans,
            "epsilon_cluster": 1e-3,
            "epsilon_harmony": 1e-2,
            "random_seed": random_seed,
            "ncores": 1,
            "harmonypy_version": harmonypy_version,
            "overwrite_existing": overwrite_existing,
            "required_backend_keywords": list(required_keywords),
        }
        harmony_metadata = dict(existing_metadata or {})
        harmony_metadata[adjusted_basis] = {
            "technical_batch_keys": keys,
            "basis": basis,
            "theta_mode": theta_mode,
            "theta": resolved_theta,
            "random_seed": random_seed,
        }
        output.uns["harmony"] = harmony_metadata
        warnings = [
            "Technical batch keys are user declarations; confounding can remove Sample or Condition biology.",
            "Harmony completion and objective change do not prove Technical-batch removal or biological conservation.",
        ]
        if basis_values.shape[1] == 1:
            warnings.append(
                "Harmony received a one-dimensional basis; correction is executable but geometric interpretation is limited."
            )
        if constant_component_indices.size:
            warnings.append(
                f"Harmony basis contains {int(constant_component_indices.size)} zero-variance component(s) at "
                f"zero-based indices {constant_component_indices.tolist()}; the backend result requires careful review."
            )
        noninformative_keys = [key for key, counts in level_counts.items() if len(counts) < 2]
        if noninformative_keys:
            warnings.append(
                "Harmony correction covariates with only one observed level provide no correction information: "
                f"{noninformative_keys}. They were preserved as explicit expert choices."
            )
        if singleton_levels:
            warnings.append(
                f"Harmony Technical-batch singleton levels may produce unstable correction: {singleton_levels}."
            )
        if theta_mode == "automatic":
            warnings.append(
                "Harmony theta was omitted so the installed harmonypy backend selected its official automatic default."
            )
        small_levels = {
            key: {level: count for level, count in counts.items() if count < 10} for key, counts in level_counts.items()
        }
        small_levels = {key: counts for key, counts in small_levels.items() if counts}
        if small_levels:
            warnings.append(f"Small Technical-batch levels may be unstable: {small_levels}.")
        finish_adata(
            output,
            "harmony_integration",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
            warnings=warnings,
        )
        code_parameters = {**parameters, "ridge_penalty": None if ridge_penalty == -1 else float(ridge_penalty)}
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellHarmonyIntegration",
            title="Harmony integration summary",
            operation="harmony_integration",
            methods=(
                f"Harmony corrected the declared PCA basis {basis!r} across user-declared Technical batch "
                f"covariates {keys!r} with harmonypy.run_harmony and stored a separate embedding at {adjusted_basis!r}. "
                f"Theta mode was {theta_mode!r}; automatic mode omits theta from the backend call."
            ),
            results=(
                f"A finite {cells:,} by {corrected.shape[1]:,} corrected embedding was produced. "
                "This structural result does not establish successful integration or biological-signal preservation."
            ),
            key_results={
                "cells": cells,
                "features": genes,
                "input_components": int(basis_values.shape[1]),
                "output_components": int(corrected.shape[1]),
                "technical_batch_levels": level_counts,
                "backend_orientation": orientation,
                "backend_z_corr_shape": list(expected_shape),
                "harmony_round_values": rounds,
                "objective_initial": float(finite_objective[0]) if rounds else None,
                "objective_final": float(finite_objective[-1]) if rounds else None,
                "displacement_norm": summarize_numeric(displacement),
                "constant_input_components": int(constant_component_indices.size),
                "constant_input_component_indices": constant_component_indices.tolist(),
                "original_basis_preserved": True,
                "corrected_finite": True,
            },
            parameters=parameters,
            references=(
                HARMONY_REFERENCE,
                HARMONYPY_REFERENCE,
                SCANPY_HARMONY_REFERENCE,
                INTEGRATION_BENCHMARK_REFERENCE,
            ),
            software_packages=(*INTEGRATION_SOFTWARE_PACKAGES, "harmonypy"),
            warnings=warnings,
            limitations=(
                "Harmony produces corrected low-dimensional coordinates, not corrected gene expression.",
                "Integration quality requires separate Technical-batch removal and biological-conservation diagnostics.",
                "The corrected embedding is exploratory and is not Sample-level Condition inference.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
            code=_harmony_code(code_parameters),
        )
        return output, report, code


class OpenBioSingleCellSCVIIntegration:
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="scVI counts source",
        default="layer",
        layer_input_id="counts_layer",
        layer_default="counts",
    )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        technical_batch_key: str = "batch",
        categorical_covariates: str = "",
        continuous_covariates: str = "",
        n_latent: int = 10,
        gene_likelihood: str = "zinb",
        n_layers: int = 1,
        dispersion: str = "gene",
        dropout_rate: float = 0.1,
        epochs: Mapping[str, object] | None = None,
        early_stopping: bool = False,
        train_size: float = 0.9,
        batch_size: int = 128,
        size_factor_key: str = "",
        accelerator: str = "auto",
        output_key: str = "X_scVI",
        overwrite_existing: bool = False,
        random_seed: int = 0,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        _require_in_memory_axes(adata, operation="scVI Integration")
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        (
            matrix,
            cell_totals,
            gene_totals,
            nonzero,
            state,
            evidence,
            integer_like,
            zero_genes,
            count_warnings,
        ) = _validate_scvi_counts(adata, expression)
        technical_batch_key = technical_batch_key.strip()
        if not technical_batch_key:
            raise ValueError("scVI Integration requires a Technical batch key.")
        categorical_keys = _comma_separated_keys(categorical_covariates)
        continuous_keys = _comma_separated_keys(continuous_covariates)
        size_factor_key = size_factor_key.strip()
        output_key = output_key.strip()
        if not output_key:
            raise ValueError("scVI output_key cannot be empty.")
        all_roles = [
            technical_batch_key,
            *categorical_keys,
            *continuous_keys,
            *([size_factor_key] if size_factor_key else []),
        ]
        if len(all_roles) != len(set(all_roles)):
            raise ValueError(
                "scVI observation columns cannot be reused across Technical-batch, covariate, and size-factor roles."
            )
        level_counts = _categorical_levels(
            adata,
            [technical_batch_key, *categorical_keys],
            role="scVI Technical nuisance",
        )
        continuous_summaries = {}
        for key in continuous_keys:
            if key not in adata.obs:
                raise ValueError(f"scVI continuous nuisance column not found in obs: {key!r}")
            try:
                values = science.np.asarray(adata.obs[key], dtype=float)
            except (TypeError, ValueError) as error:
                raise ValueError(f"scVI continuous nuisance column must be numeric: {key!r}") from error
            if not bool(science.np.isfinite(values).all()):
                raise ValueError(f"scVI continuous nuisance column must be finite: {key!r}")
            continuous_summaries[key] = summarize_numeric(values)
        if size_factor_key:
            if size_factor_key not in adata.obs:
                raise ValueError(f"scVI size-factor column not found in obs: {size_factor_key!r}")
            try:
                factors = science.np.asarray(adata.obs[size_factor_key], dtype=float)
            except (TypeError, ValueError) as error:
                raise ValueError("scVI size factors must be numeric and on the linear scale.") from error
            if not bool(science.np.isfinite(factors).all()) or bool((factors <= 0).any()):
                raise ValueError("scVI size factors must be finite, positive, and on the linear scale.")
        if isinstance(n_latent, bool) or not isinstance(n_latent, int) or n_latent < 1:
            raise ValueError("scVI n_latent must be a positive integer.")
        if gene_likelihood not in {"zinb", "nb", "poisson"}:
            raise ValueError(f"Unsupported scVI gene_likelihood: {gene_likelihood!r}")
        if dispersion not in {"gene", "gene-batch"}:
            raise ValueError(f"Unsupported scVI dispersion: {dispersion!r}")
        if isinstance(n_layers, bool) or not isinstance(n_layers, int) or n_layers < 1:
            raise ValueError("scVI n_layers must be a positive integer.")
        if (
            isinstance(dropout_rate, bool)
            or not isinstance(dropout_rate, (int, float))
            or not science.np.isfinite(dropout_rate)
            or not 0 <= dropout_rate <= 1
        ):
            raise ValueError("scVI dropout_rate must be finite in [0, 1].")
        epoch_mode, resolved_max_epochs = _resolve_epochs(epochs)
        if not science.np.isfinite(train_size) or not 0 < train_size <= 1:
            raise ValueError("scVI train_size must be finite in (0, 1].")
        if train_size == 1 and early_stopping:
            raise ValueError("scVI early_stopping requires a validation holdout, so train_size must be less than 1.")
        validation_interval = None if train_size == 1 else 1
        if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
            raise ValueError("scVI batch_size must be a positive integer.")
        if accelerator not in {"auto", "cpu", "gpu", "mps"}:
            raise ValueError(f"Unsupported scVI accelerator: {accelerator!r}")
        if isinstance(random_seed, bool) or not isinstance(random_seed, int) or not 0 <= random_seed <= 2**31 - 1:
            raise ValueError("scVI random_seed must be an integer in [0, 2**31 - 1].")
        if not isinstance(early_stopping, bool):
            raise TypeError("scVI early_stopping must be a boolean.")
        if not isinstance(overwrite_existing, bool):
            raise TypeError("scVI overwrite_existing must be a boolean.")
        if output_key in adata.obsm and not overwrite_existing:
            raise ValueError(f"scVI output key already exists: {output_key!r}")

        try:
            scvi = importlib.import_module("scvi")
            torch = importlib.import_module("torch")
        except (ImportError, OSError) as error:
            raise RuntimeError("scVI Integration requires the optional scvi-tools package.") from error
        model_class = scvi.model.SCVI
        _supports_keywords(
            model_class.setup_anndata,
            ("layer", "batch_key", "size_factor_key", "categorical_covariate_keys", "continuous_covariate_keys"),
            description="scvi.model.SCVI.setup_anndata",
        )
        _supports_keywords(
            model_class,
            ("n_layers", "n_latent", "gene_likelihood", "dispersion", "dropout_rate"),
            description="scvi.model.SCVI",
        )
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        # scvi-tools setup mutates registration state and the saved model must own that exact training view.
        training_adata = adata.copy()
        with SCVI_GLOBAL_RNG_LOCK:
            previous_seed = getattr(scvi.settings, "seed", None)
            python_state = random.getstate()
            numpy_state = science.np.random.get_state()
            torch_state = torch.random.get_rng_state()
            cuda_states = torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None
            mps_module = getattr(torch, "mps", None)
            mps_backend = getattr(getattr(torch, "backends", None), "mps", None)
            mps_available = (
                mps_module is not None
                and callable(getattr(mps_module, "get_rng_state", None))
                and callable(getattr(mps_module, "set_rng_state", None))
                and mps_backend is not None
                and callable(getattr(mps_backend, "is_available", None))
                and mps_backend.is_available()
            )
            mps_state = mps_module.get_rng_state() if mps_available else None
            try:
                with torch.inference_mode(False):
                    scvi.settings.seed = random_seed
                    with python_warnings.catch_warnings():
                        python_warnings.filterwarnings(
                            "ignore",
                            message=r".*does not contain unnormalized count data.*",
                            category=UserWarning,
                        )
                        model_class.setup_anndata(
                            training_adata,
                            layer=expression.scanpy_layer,
                            batch_key=technical_batch_key,
                            size_factor_key=size_factor_key or None,
                            categorical_covariate_keys=categorical_keys or None,
                            continuous_covariate_keys=continuous_keys or None,
                        )
                    model = model_class(
                        training_adata,
                        n_layers=n_layers,
                        n_latent=n_latent,
                        gene_likelihood=gene_likelihood,
                        dispersion=dispersion,
                        dropout_rate=float(dropout_rate),
                    )
                    _supports_keywords(
                        model.train,
                        (
                            "max_epochs",
                            "accelerator",
                            "devices",
                            "train_size",
                            "validation_size",
                            "batch_size",
                            "early_stopping",
                        ),
                        description="scvi.model.SCVI.train",
                    )
                    model.train(
                        max_epochs=resolved_max_epochs,
                        accelerator=accelerator,
                        devices=1,
                        train_size=float(train_size),
                        validation_size=None,
                        batch_size=batch_size,
                        early_stopping=early_stopping,
                        early_stopping_monitor="elbo_validation",
                        early_stopping_patience=10,
                        early_stopping_min_delta=0.0,
                        check_val_every_n_epoch=validation_interval,
                    )
                    latent = science.np.asarray(model.get_latent_representation(), dtype=float)
            finally:
                scvi.settings.seed = previous_seed
                random.setstate(python_state)
                science.np.random.set_state(numpy_state)
                torch.random.set_rng_state(torch_state)
                if cuda_states is not None:
                    torch.cuda.set_rng_state_all(cuda_states)
                if mps_state is not None:
                    mps_module.set_rng_state(mps_state)
        if not bool(getattr(model, "is_trained", False)):
            raise RuntimeError("scVI backend did not report a trained model after training.")
        if latent.shape != (cells, n_latent) or not bool(science.np.isfinite(latent).all()):
            raise RuntimeError("scVI returned an invalid latent representation.")
        with science.np.errstate(over="ignore", invalid="ignore"):
            latent_variance = science.np.var(latent, axis=0)
        constant_latent_indices = science.np.flatnonzero(
            science.np.max(latent, axis=0) == science.np.min(latent, axis=0)
        )
        diagnostics = _scvi_training_diagnostics(model)
        parameters = {
            **expression.parameters(),
            "count_source_state": state,
            "count_source_state_evidence": evidence,
            "technical_batch_key": technical_batch_key,
            "size_factor_key": size_factor_key,
            "categorical_covariates": categorical_keys,
            "continuous_covariates": continuous_keys,
            "n_latent": n_latent,
            "gene_likelihood": gene_likelihood,
            "n_layers": n_layers,
            "dispersion": dispersion,
            "dropout_rate": float(dropout_rate),
            "epochs_mode": epoch_mode,
            "max_epochs": resolved_max_epochs,
            "early_stopping": early_stopping,
            "early_stopping_monitor": "elbo_validation",
            "early_stopping_patience": 10,
            "early_stopping_min_delta": 0.0,
            "check_val_every_n_epoch": validation_interval,
            "train_size": float(train_size),
            "validation_size": None,
            "batch_size": batch_size,
            "accelerator": accelerator,
            "devices": 1,
            "output_key": output_key,
            "overwrite_existing": overwrite_existing,
            "random_seed": random_seed,
        }
        output = adata
        output.obsm[output_key] = latent.copy()
        warnings = [
            *count_warnings,
            "Technical batch and nuisance covariates are user declarations; confounding can remove Sample or Condition biology.",
            "Training completion and loss values do not prove Technical-batch removal or biological conservation.",
            "The trained scVI model is a session-only native artifact bound to the producing Python Worker.",
            "A fixed seed does not guarantee bitwise equality across devices, hardware, or software versions.",
        ]
        if state == "unknown":
            warnings.append(
                "The selected count-intended source has no recorded raw-count provenance; verify whether it contains "
                "unnormalized UMI counts or another explicitly intended non-negative representation."
            )
        elif state != "counts":
            warnings.append(
                f"The explicitly selected scVI source has recorded expression state {state!r} ({evidence or 'no evidence text'}). "
                "The backend can train on finite non-negative values, but count-likelihood interpretation may not apply."
            )
        if zero_genes:
            warnings.append(
                f"The selected scVI source contains {zero_genes} all-zero gene(s); the audited backend can retain them, "
                "but they provide no expression information."
            )
        if n_latent >= min(cells, genes):
            warnings.append(
                f"scVI n_latent={n_latent} is at least one input-axis dimension; this overcomplete architecture is "
                "executable but may be inefficient or weakly identified."
            )
        if n_layers > 20:
            warnings.append(
                f"scVI n_layers={n_layers} requests an unusually deep architecture; resource use and optimization "
                "stability require review."
            )
        if dropout_rate == 1:
            warnings.append(
                "scVI dropout_rate=1 drops every affected activation; the backend accepts this endpoint, but the "
                "resulting architecture is degenerate and requires explicit review."
            )
        if train_size < 0.5:
            warnings.append(
                f"scVI train_size={train_size:g} uses fewer than half of cells for training; this is executable but "
                "may reduce model stability."
            )
        if train_size == 1:
            warnings.append(
                "scVI train_size=1 uses every cell for training and leaves no validation holdout; validation metrics "
                "and early-stopping interpretation may be unavailable."
            )
        if constant_latent_indices.size:
            warnings.append(
                f"scVI returned {int(constant_latent_indices.size)} constant latent dimension(s) at zero-based indices "
                f"{constant_latent_indices.tolist()}; downstream geometry may be degenerate."
            )
        small_levels = {
            key: {level: count for level, count in counts.items() if count < 10} for key, counts in level_counts.items()
        }
        small_levels = {key: counts for key, counts in small_levels.items() if counts}
        if small_levels:
            warnings.append(f"Small Technical-batch or categorical-covariate levels may be unstable: {small_levels}.")
        noninformative_categorical = [key for key, counts in level_counts.items() if len(counts) < 2]
        if noninformative_categorical:
            warnings.append(
                "scVI nuisance columns with only one observed level provide no integration information: "
                f"{noninformative_categorical}. They were preserved as explicit expert choices."
            )
        noninformative_continuous = [
            key for key, values in continuous_summaries.items() if values["min"] == values["max"]
        ]
        if noninformative_continuous:
            warnings.append(
                "Constant scVI continuous nuisance columns provide no correction information: "
                f"{noninformative_continuous}. They were preserved as explicit expert choices."
            )
        if genes < 1000 or genes > 10000:
            warnings.append(
                f"The model used {genes:,} features; scvi-tools generally recommends approximately 1,000-10,000 HVGs."
            )
        if diagnostics["actual_epochs"] is None:
            warnings.append("The installed backend did not expose a usable training-history length.")
        finish_adata(
            output, "scvi_integration", parameters, cells, genes, started_at, random_seed=random_seed, warnings=warnings
        )
        trained_model = model
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellSCVIIntegration",
            title="scVI integration summary",
            operation="scvi_integration",
            methods=(
                f"scVI was trained on explicitly selected count-intended source {_source_label(expression)!r} "
                f"(resolved state {state!r}) with Technical batch {technical_batch_key!r}; the posterior-mean latent "
                f"representation was stored at {output_key!r}."
            ),
            results=(
                f"A trained session-only native model artifact and finite {cells:,} by {n_latent:,} latent representation were produced. "
                "This does not by itself establish successful integration or preservation of Condition biology."
            ),
            key_results={
                "cells": cells,
                "features": genes,
                "count_source": _source_label(expression),
                "count_source_state": state,
                "count_source_state_evidence": evidence,
                "count_source_integer_like": integer_like,
                "zero_total_genes": zero_genes,
                "count_nonzero": nonzero,
                "total_counts": float(cell_totals.sum()),
                "cell_totals": summarize_numeric(cell_totals),
                "gene_totals": summarize_numeric(gene_totals),
                "technical_nuisance_levels": level_counts,
                "continuous_nuisance_summaries": continuous_summaries,
                "training": diagnostics,
                "latent_shape": [cells, n_latent],
                "latent_variance": summarize_numeric(latent_variance),
                "constant_latent_dimensions": int(constant_latent_indices.size),
                "constant_latent_dimension_indices": constant_latent_indices.tolist(),
                "latent_finite": True,
                "model_session_only": True,
            },
            parameters=parameters,
            references=(SCVI_REFERENCE, SCVI_TOOLS_REFERENCE, SCVI_API_REFERENCE, INTEGRATION_BENCHMARK_REFERENCE),
            software_packages=(*INTEGRATION_SOFTWARE_PACKAGES, "scvi-tools", "torch", "lightning"),
            warnings=warnings,
            limitations=(
                "The latent representation is exploratory and is not corrected gene expression.",
                "Integration quality requires separate Technical-batch removal and biological-conservation diagnostics.",
                "This HVG-view representation is not full-gene or Sample-level Condition inference.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
            code=_scvi_code(parameters),
        )
        return output, trained_model, report, code


class OpenBioSingleCellLeidenResolutionSweep:
    @classmethod
    def execute(
        cls,
        adata: AnnData,
        resolutions: str = "0.25,0.5,1.0,2.0",
        key_prefix: str = "leiden",
        neighbors_key: str = "neighbors",
        n_iterations: int = 2,
        stability_repeats: int = 5,
        random_seed: int = 0,
        overwrite_existing: bool = False,
    ) -> tuple[Any, ...]:
        science = dependencies.require_scientific_dependencies()
        if not isinstance(key_prefix, str):
            raise TypeError("Leiden key_prefix must be a string.")
        key_prefix = key_prefix.strip()
        if not key_prefix:
            raise ValueError("Leiden key_prefix cannot be empty.")
        values = _parse_resolutions(resolutions)
        _, n_iterations, stability_repeats = validate_leiden_settings(values[0], n_iterations, stability_repeats)
        for resolution in values[1:]:
            validate_leiden_settings(resolution, n_iterations, stability_repeats)
        random_seed = validate_random_seed(random_seed, stability_repeats=stability_repeats)
        total_runs = len(values) * stability_repeats
        keys = [_resolution_key(key_prefix, resolution) for resolution in values]
        if len(keys) != len(set(keys)):
            raise ValueError("Leiden resolutions generate colliding output keys.")
        graph = resolve_named_graph(adata, neighbors_key, operation="Leiden resolution sweep")
        if graph.neighbors_key in keys:
            raise ValueError(
                "Leiden resolution output keys cannot equal the resolved neighbors_key because that uns key stores graph metadata."
            )
        for key in keys:
            validate_graph_result_key(
                adata,
                key,
                operation="Leiden resolution sweep",
                obs=True,
                uns=True,
                overwrite_existing=overwrite_existing,
            )

        fixed_policy = {
            "flavor": "igraph",
            "directed": False,
            "use_weights": True,
            "objective_function": "modularity",
        }
        parameters = {
            "resolutions": values,
            "keys": keys,
            "key_prefix": key_prefix,
            "neighbors_key": graph.neighbors_key,
            "n_iterations": n_iterations,
            "stability_repeats": stability_repeats,
            "random_seed": random_seed,
            "overwrite_existing": overwrite_existing,
            "total_runs": total_runs,
            "fixed_policy": fixed_policy,
        }
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output = adata
        if overwrite_existing:
            for key in keys:
                output.obs.drop(columns=[key], inplace=True, errors="ignore")
                output.uns.pop(key, None)

        from sklearn.metrics import adjusted_rand_score

        rows: list[dict[str, Any]] = []
        cluster_sizes_by_resolution: dict[str, dict[str, int]] = {}
        previous_resolution: float | None = None
        previous_membership: tuple[str, ...] | None = None
        for resolution, key in zip(values, keys, strict=True):
            partition = run_leiden_partition(
                output,
                resolution=resolution,
                key_added=key,
                graph=graph,
                random_seed=random_seed,
                n_iterations=n_iterations,
            )
            metadata = output.uns.get(key)
            backend_value = metadata.get("modularity") if isinstance(metadata, Mapping) else None
            if isinstance(backend_value, bool):
                raise RuntimeError("Leiden backend modularity is missing or invalid.")
            try:
                backend_modularity = float(backend_value)
            except (TypeError, ValueError) as error:
                raise RuntimeError("Leiden backend modularity is missing or invalid.") from error
            if not bool(science.np.isfinite(backend_modularity)):
                raise RuntimeError("Leiden backend modularity is non-finite.")
            if partition.modularity is None or not bool(
                science.np.isclose(backend_modularity, partition.modularity, rtol=1e-9, atol=1e-12)
            ):
                raise RuntimeError("Leiden backend modularity is inconsistent with the validated graph partition.")

            stability = assess_leiden_stability(
                output,
                graph,
                base_membership=partition.membership,
                resolution=resolution,
                random_seed=random_seed,
                n_iterations=n_iterations,
                repeats=stability_repeats,
            )
            if any(
                not bool(science.np.isfinite(value)) or not -1.0 <= value <= 1.0 for value in stability.pairwise_ari
            ):
                raise RuntimeError("Leiden stability ARI diagnostics are invalid.")
            cluster_sizes = dict(partition.cluster_sizes)
            counts = list(cluster_sizes.values())
            adjacent_ari = (
                None
                if previous_membership is None
                else float(adjusted_rand_score(previous_membership, partition.membership))
            )
            if adjacent_ari is not None and (
                not bool(science.np.isfinite(adjacent_ari)) or not -1.0 <= adjacent_ari <= 1.0
            ):
                raise RuntimeError("Leiden adjacent-resolution ARI is invalid.")
            row = {
                "resolution": resolution,
                "key": key,
                "n_clusters": len(cluster_sizes),
                "min_cluster_size": min(counts),
                "median_cluster_size": float(science.np.median(counts)),
                "max_cluster_size": max(counts),
                "smallest_cluster_fraction": float(min(counts) / cells),
                "largest_cluster_fraction": float(max(counts) / cells),
                "singleton_clusters": partition.singleton_count,
                "modularity": backend_modularity,
                "stability_repeats": stability.starts,
                "stability_mean_ari": stability.mean_ari,
                "stability_min_ari": stability.min_ari,
                "stability_max_ari": stability.max_ari,
                "adjacent_previous_resolution": previous_resolution,
                "adjacent_resolution_ari": adjacent_ari,
            }
            rows.append(row)
            cluster_sizes_by_resolution[key] = cluster_sizes
            output.uns[key]["openbio_diagnostics"] = {
                "modularity": backend_modularity,
                "cluster_sizes": cluster_sizes,
                "singleton_count": partition.singleton_count,
                "stability": {
                    "starts": stability.starts,
                    "seeds": list(stability.seeds),
                    "pairwise_ari": list(stability.pairwise_ari),
                    "mean_ari": stability.mean_ari,
                    "min_ari": stability.min_ari,
                    "max_ari": stability.max_ari,
                },
            }
            previous_resolution = resolution
            previous_membership = partition.membership

        finish_adata(
            output,
            "leiden_resolution_sweep",
            parameters,
            cells,
            genes,
            started_at,
            random_seed=random_seed,
        )
        warnings = []
        if stability_repeats == 1:
            warnings.append("Leiden sweep stability was not assessed because stability_repeats=1.")
        if any(resolution == 0 for resolution in values):
            warnings.append("Leiden sweep includes resolution=0, which commonly produces very few communities.")
        if n_iterations == 0:
            warnings.append("Leiden sweep n_iterations=0 requests no optimization iterations.")
        if total_runs > _LEIDEN_SWEEP_WORKLOAD_WARNING:
            warnings.append(f"Leiden sweep requests {total_runs} backend runs; execution may be resource intensive.")
        graph_report = graph_diagnostics(graph)
        graph_report["isolated_cells"] = graph_report["isolated_observation_count"]
        if graph_report["connected_components"] > 1:
            warnings.append(
                "The named graph has multiple connected components, which imposes a lower bound on community count."
            )
        if graph_report["isolated_observation_count"]:
            warnings.append("The named graph contains isolated observations; inspect their assigned singleton groups.")
        if any(row["singleton_clusters"] for row in rows):
            warnings.append("At least one surveyed partition contains singleton communities requiring review.")
        if any(row["largest_cluster_fraction"] > 0.8 for row in rows):
            warnings.append("At least one surveyed partition is dominated by a community containing over 80% of cells.")
        if any(row["stability_min_ari"] is not None and row["stability_min_ari"] < 0.9 for row in rows):
            warnings.append(
                "At least one resolution has minimum repeat-start ARI below 0.9; inspect optimization stability."
            )

        code = _leiden_resolution_sweep_code(parameters)
        table = science.pd.DataFrame(rows, columns=list(_LEIDEN_METRIC_COLUMNS))
        for column in (
            "stability_mean_ari",
            "stability_min_ari",
            "stability_max_ari",
            "adjacent_previous_resolution",
            "adjacent_resolution_ari",
        ):
            table[column] = science.pd.array(table[column], dtype="Float64")
        table_result = make_table_result(
            table=table,
            title="Leiden resolution metrics",
            operation="leiden_resolution_sweep",
            parameters=parameters,
            description=(
                "One row per declared resolution with partition sizes, resolution-aware modularity, "
                "repeat-start ARI, and adjacent-resolution ARI."
            ),
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        cluster_counts = [int(row["n_clusters"]) for row in rows]
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellLeidenResolutionSweep",
            title="Leiden resolution sweep",
            operation="leiden_resolution_sweep",
            methods=(
                f"Surveyed {len(values)} declared resolutions on one canonical weighted undirected graph with "
                f"Scanpy's igraph Leiden modularity objective, {n_iterations} iterations, and "
                f"{stability_repeats} total starts per resolution including the base start."
            ),
            results=(
                f"Across the declared graph granularities, community counts ranged from {min(cluster_counts)} "
                f"to {max(cluster_counts)}; inspect the structured size and ARI diagnostics alongside biological "
                "marker evidence and upstream graph choices."
            ),
            key_results={
                "resolution_count": len(values),
                "total_runs": total_runs,
                "stability_seeds": [random_seed + offset for offset in range(stability_repeats)],
                "resolution_metrics": rows,
                "cluster_sizes_by_resolution": cluster_sizes_by_resolution,
                "graph": graph_report,
            },
            parameters=parameters,
            references=(
                LEIDEN_REFERENCE,
                SCANPY_REFERENCE,
                *([ARI_REFERENCE] if stability_repeats > 1 else []),
            ),
            software_packages=(
                "scanpy",
                "anndata",
                "igraph",
                "numpy",
                "pandas",
                "scipy",
                "scikit-learn",
            ),
            warnings=warnings,
            limitations=(
                "Resolution and modularity are graph-dependent diagnostics, not automatic selection criteria.",
                "Optimization stability does not establish biological validity or curated cell identities.",
                "Graph communities are exploratory pooled-cell partitions, not Sample-level Condition inference.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            code=code,
            random_seed=random_seed,
        )
        return output, table_result, report, code


TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
TABLE_CODEC = "table-jsonl-v1"
SCVI_MODEL_KIND = "OPENBIO_SCVI_MODEL"
SCVI_MODEL_CODEC = "scvi-native-directory"


def _adata_input(inputs: dict[str, JSONValue], *, operation: str) -> Any:
    require_input_names(inputs, {"adata"}, operation=operation)
    root = require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    return read_anndata(root)


def _standard_records(context: OperationContext, output: Any, summary: Any, code: str) -> list[JSONValue]:
    return [
        write_anndata_output(context, output),
        {"type": "summary", "name": "summary", "value": result_metadata(summary)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.harmonyintegration")
def harmony_integration(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    expected = {
        "technical_batch_keys",
        "basis",
        "theta",
        "ridge_penalty",
        "sigma",
        "n_clusters",
        "tau",
        "adjusted_basis",
        "overwrite_existing",
        "max_iter_harmony",
        "max_iter_kmeans",
        "random_seed",
    }
    require_parameters(parameters, expected, operation="Harmony Integration")
    output, summary, code = OpenBioSingleCellHarmonyIntegration.execute(
        _adata_input(inputs, operation="Harmony Integration"), **parameters
    )
    return _standard_records(context, output, summary, code)


@register_operation("openbio.node.scviintegration")
def scvi_integration(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    expected = {
        "source",
        "technical_batch_key",
        "categorical_covariates",
        "continuous_covariates",
        "n_latent",
        "gene_likelihood",
        "n_layers",
        "dispersion",
        "dropout_rate",
        "epochs",
        "early_stopping",
        "train_size",
        "batch_size",
        "size_factor_key",
        "accelerator",
        "output_key",
        "overwrite_existing",
        "random_seed",
    }
    require_parameters(parameters, expected, operation="scVI Integration")
    output, model, summary, code = OpenBioSingleCellSCVIIntegration.execute(
        _adata_input(inputs, operation="scVI Integration"), **parameters
    )
    save = getattr(model, "save", None)
    if not callable(save):
        raise RuntimeError("scVI model does not expose the required official save API.")
    container = context.create_output_directory("model")
    model_root = container / "native"
    save(str(model_root), overwrite=False, save_anndata=True)
    if not model_root.is_dir():
        raise RuntimeError("scVI model.save did not create its native model directory.")
    return [
        write_anndata_output(context, output),
        {
            "type": "artifact",
            "name": "model",
            "kind": SCVI_MODEL_KIND,
            "codec": SCVI_MODEL_CODEC,
            "payload": model_root.relative_to(context.output_root).as_posix(),
        },
        {"type": "summary", "name": "summary", "value": result_metadata(summary)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.leidenresolutionsweep")
def leiden_resolution_sweep(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    expected = {
        "resolutions",
        "key_prefix",
        "neighbors_key",
        "n_iterations",
        "stability_repeats",
        "random_seed",
        "overwrite_existing",
    }
    require_parameters(parameters, expected, operation="Leiden Resolution Sweep")
    output, table, summary, code = OpenBioSingleCellLeidenResolutionSweep.execute(
        _adata_input(inputs, operation="Leiden Resolution Sweep"), **parameters
    )
    table_root = context.create_output_directory("resolution_metrics")
    write_table(table_root, table.table, result_metadata(table))
    return [
        write_anndata_output(context, output),
        {
            "type": "artifact",
            "name": "resolution_metrics",
            "kind": TABLE_KIND,
            "codec": TABLE_CODEC,
            "payload": table_root.relative_to(context.output_root).as_posix(),
        },
        {"type": "summary", "name": "summary", "value": result_metadata(summary)},
        {"type": "string", "name": "code", "value": code},
    ]


__all__ = ["harmony_integration", "leiden_resolution_sweep", "scvi_integration"]
