from __future__ import annotations

import hashlib
import inspect
import textwrap
from typing import Any

from .enrichment_resource import gene_set_resource_code, load_gene_set_resource
from .score_artifact import (
    SCORE_ARTIFACT_SCHEMA_VERSION,
    SCORE_ARTIFACTS_KEY,
    _canonical_json_sha256,
    _canonical_string_axis,
    build_score_artifact,
    score_frame_fingerprint,
    store_score_artifact,
)


def gene_set_resource_cache_fingerprint(
    path: str,
    identity: tuple[Any, ...],
) -> tuple[Any, ...]:
    """Bind Comfy cache invalidation to resource bytes as well as filesystem identity."""
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return ("openbio-gene-set-scoring-resource-v1", *identity, digest.hexdigest())


def _standalone_score_gene_sets(
    adata,
    *,
    method,
    gene_sets_path,
    requested_resource_path,
    resource_metadata_json,
    source_kind="X",
    layer_name=None,
    source_column="geneset",
    target_column="genesymbol",
    min_targets=5,
    output_key,
    overwrite_existing=False,
    expected_resource_sha256=None,
    n_top_features=0,
    batch_size=250000,
    max_output_rows=2000000,
    max_working_memory_gib=4.0,
    kernel="gaussian_normalized",
    maxdiff=True,
    absrnk=False,
    tau=1.0,
    panel="",
    ctrl_size=0,
    n_bins=25,
    random_seed=0,
    openbio_version="unknown",
):
    """Run one reviewed observation-level gene-set scoring method on a private expression view."""
    import hashlib
    import importlib
    import importlib.metadata
    import inspect as runtime_inspect
    import json
    import math
    import platform
    import warnings as runtime_warnings
    from collections.abc import Mapping

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import sparse

    supported_methods = {"aucell", "gsva", "panel"}
    if method not in supported_methods:
        raise ValueError(f"Gene-set scoring method must be one of {sorted(supported_methods)}; received {method!r}.")
    operation = {"aucell": "AUCell scores", "gsva": "GSVA scores", "panel": "Gene panel score"}[method]
    producer_node = {
        "aucell": "OpenBioSingleCellAUCellScores",
        "gsva": "OpenBioSingleCellGSVAScores",
        "panel": "OpenBioSingleCellGenePanelScores",
    }[method]
    method_display = {"aucell": "AUCell", "gsva": "GSVA", "panel": "Scanpy score_genes"}[method]

    if not hasattr(adata, "X") or not hasattr(adata, "obs_names") or not hasattr(adata, "var_names"):
        raise TypeError(f"{operation} requires an AnnData-like input with expression and named axes.")
    if int(getattr(adata, "n_obs", len(adata.obs_names))) <= 0:
        raise ValueError(f"{operation} requires at least one observation.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    observation_names, observation_axis_sha256 = _canonical_string_axis(
        adata.obs_names.tolist(), description=f"{operation} observation identifier"
    )
    if method == "gsva" and len(observation_names) < 2:
        raise ValueError("GSVA scores require at least two observations for cohort density estimation.")
    if source_kind not in {"X", "layer", "raw"}:
        raise ValueError(f"{operation} source_kind must be 'X', 'layer', or 'raw'.")
    if source_kind == "layer":
        if not isinstance(layer_name, str) or not layer_name or layer_name != layer_name.strip():
            raise ValueError(f"{operation} layer_name must be a nonblank, whitespace-canonical string.")
        if layer_name not in adata.layers:
            raise ValueError(f"{operation} expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
        feature_values = adata.var_names.tolist()
    elif source_kind == "raw":
        if getattr(adata, "raw", None) is None:
            raise ValueError(f"{operation} Raw snapshot was selected but adata.raw is unavailable.")
        matrix = adata.raw.X
        feature_values = adata.raw.var_names.tolist()
        layer_name = None
    else:
        matrix = adata.X
        feature_values = adata.var_names.tolist()
        layer_name = None
    feature_names, feature_axis_sha256 = _canonical_string_axis(
        feature_values, description=f"{operation} feature identifier"
    )
    if len(feature_names) < 2:
        raise ValueError(f"{operation} requires at least two uniquely named expression features.")
    if tuple(getattr(matrix, "shape", ())) != (len(observation_names), len(feature_names)):
        raise ValueError(f"{operation} selected expression matrix does not align to its named axes.")
    if not (sparse.issparse(matrix) or isinstance(matrix, np.ndarray)):
        raise TypeError(f"{operation} requires an in-memory NumPy or SciPy sparse expression matrix.")
    stored_values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    stored_values = np.asarray(stored_values)
    if stored_values.dtype.kind not in "iuf":
        raise TypeError(f"{operation} selected expression must be numeric.")
    if stored_values.size and not bool(np.isfinite(stored_values).all()):
        raise ValueError(f"{operation} selected expression contains non-finite values.")
    if sparse.issparse(matrix):
        canonical_sparse = sparse.csr_matrix(matrix, copy=True)
        canonical_sparse.sum_duplicates()
        canonical_sparse.eliminate_zeros()
        canonical_sparse.sort_indices()
        nonzero_per_observation = np.asarray(canonical_sparse.getnnz(axis=1)).ravel()
        semantic_values = canonical_sparse.data
    else:
        dense_view = np.asarray(matrix)
        nonzero_per_observation = np.count_nonzero(dense_view, axis=1)
        semantic_values = dense_view.ravel()
    zero_observations = int(np.count_nonzero(nonzero_per_observation == 0))
    warnings = []
    if zero_observations:
        warnings.append(
            f"{operation} retained {zero_observations} all-zero observations; their scores can reflect "
            "ties or background values rather than expression evidence."
        )
    value_min = float(np.min(semantic_values)) if semantic_values.size else 0.0
    value_max = float(np.max(semantic_values)) if semantic_values.size else 0.0
    integer_like = bool(np.allclose(semantic_values, np.rint(semantic_values), rtol=0.0, atol=1e-8))

    if method == "gsva" and kernel == "poisson_counts":
        if value_min < 0.0 or not integer_like:
            raise ValueError("Poisson GSVA requires finite nonnegative integer count values.")
    else:
        if method == "gsva" and kernel not in {"gaussian_normalized", "empirical"}:
            raise ValueError("GSVA kernel must be 'gaussian_normalized', 'poisson_counts', or 'empirical'.")
        if integer_like and value_min >= 0.0:
            warnings.append(
                f"{operation} is using nonnegative integer-like/count-like values. Normalized continuous expression "
                "is recommended for this configuration; sparsity, ties, and score scale can differ."
            )

    def expression_fingerprint():
        digest = hashlib.sha256()
        digest.update(b"openbio-singlecell/scoring-expression/v1\0")
        digest.update(observation_axis_sha256.encode("ascii"))
        digest.update(feature_axis_sha256.encode("ascii"))
        digest.update(str(tuple(matrix.shape)).encode("ascii"))
        if sparse.issparse(matrix):
            csr = sparse.csr_matrix(matrix, copy=True)
            csr.sum_duplicates()
            csr.eliminate_zeros()
            csr.sort_indices()
            digest.update(b"csr\0")
            digest.update(np.ascontiguousarray(csr.indptr, dtype="<i8").tobytes())
            digest.update(np.ascontiguousarray(csr.indices, dtype="<i8").tobytes())
            digest.update(np.ascontiguousarray(csr.data, dtype="<f8").tobytes())
        else:
            digest.update(b"dense\0")
            dense = np.asarray(matrix)
            for start in range(0, dense.shape[0], 1024):
                digest.update(np.ascontiguousarray(dense[start : start + 1024], dtype="<f8").tobytes())
        return digest.hexdigest()

    expression_sha256 = expression_fingerprint()
    if not isinstance(output_key, str) or not output_key or output_key != output_key.strip():
        raise ValueError(f"{operation} output_key must be a nonblank, whitespace-canonical string.")
    if not isinstance(overwrite_existing, bool):
        raise TypeError(f"{operation} overwrite_existing must be boolean.")
    storage = "obs" if method == "panel" else "obsm"
    storage_mapping = adata.obs if storage == "obs" else adata.obsm
    provenance_mapping = adata.uns.get(SCORE_ARTIFACTS_KEY, {})
    if not isinstance(provenance_mapping, Mapping):
        raise ValueError(f"{operation} existing score-artifact provenance must be a mapping.")
    collision = output_key in storage_mapping or output_key in provenance_mapping
    if collision and not overwrite_existing:
        raise ValueError(
            f"{operation} output {storage}[{output_key!r}] or its provenance already exists; "
            "enable overwrite_existing explicitly to replace it."
        )
    if isinstance(min_targets, bool) or not isinstance(min_targets, int) or min_targets < 1:
        raise TypeError(f"{operation} min_targets must be a positive integer.")

    resource = load_gene_set_resource(
        gene_sets_path,
        resource_metadata_json,
        source_column,
        target_column,
        feature_names,
        1 if method == "panel" else min_targets,
        expected_sha256=expected_resource_sha256,
        operation=operation,
    )
    warnings.extend(resource["warnings"])
    if not isinstance(requested_resource_path, str) or not requested_resource_path.strip():
        raise ValueError(f"{operation} requested_resource_path must be a nonblank string.")
    retained_sources = list(resource["retained_sources"])
    resource_contract_sha256 = _canonical_json_sha256(
        {
            "sha256": resource["sha256"],
            "metadata": resource["metadata"],
            "source_order": list(resource["source_order"]),
            "targets_before": {key: list(value) for key, value in resource["targets_before"].items()},
            "targets_in_universe": {key: list(value) for key, value in resource["targets_in_universe"].items()},
            "min_targets": 1 if method == "panel" else min_targets,
        }
    )
    set_accounting = resource["accounting"]["set_accounting"]
    accounting_preview_limit = 500
    resource_provenance = {
        "requested_path": requested_resource_path,
        "resolved_path": resource["resolved_path"],
        "size_bytes": resource["size_bytes"],
        "sha256": resource["sha256"],
        "resource_contract_sha256": resource_contract_sha256,
        "metadata": resource["metadata"],
        "extension": resource["accounting"]["extension"],
        "input_pairs": resource["accounting"]["input_pairs"],
        "unique_pairs": resource["accounting"]["unique_pairs"],
        "duplicate_pairs_removed": resource["accounting"]["duplicate_pairs_removed"],
        "source_count_before": resource["accounting"]["source_count_before"],
        "source_count_after_min_targets": resource["accounting"]["source_count_after_min_targets"],
        "sources_removed_by_min_targets_count": resource["accounting"]["sources_removed_by_min_targets_count"],
        "sources_removed_by_min_targets": resource["accounting"]["sources_removed_by_min_targets"][:100],
        "sources_removed_preview_truncated": len(resource["accounting"]["sources_removed_by_min_targets"]) > 100,
        "targets_before_universe_count": resource["accounting"]["targets_before_universe_count"],
        "targets_in_universe_count": resource["accounting"]["targets_in_universe_count"],
        "retained_target_rows_in_universe": resource["accounting"]["retained_target_rows_in_universe"],
        "set_accounting_preview": set_accounting[:accounting_preview_limit],
        "set_accounting_count": len(set_accounting),
        "set_accounting_preview_truncated": len(set_accounting) > accounting_preview_limit,
    }
    expression_provenance = {
        "source": source_kind,
        "layer_name": layer_name if source_kind == "layer" else None,
        "sparse": bool(sparse.issparse(matrix)),
        "observations": len(observation_names),
        "features": len(feature_names),
        "observation_axis_sha256": observation_axis_sha256,
        "feature_axis_sha256": feature_axis_sha256,
        "expression_content_sha256": expression_sha256,
        "value_min": value_min,
        "value_max": value_max,
        "integer_like_nonzero_values": integer_like,
        "zero_observations": zero_observations,
    }

    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
        raise TypeError(f"{operation} batch_size must be a positive integer.")
    if batch_size > 2**31 - 1:
        raise ValueError(f"{operation} batch_size exceeds the supported bound.")
    if isinstance(max_output_rows, bool) or not isinstance(max_output_rows, int) or max_output_rows < 1:
        raise TypeError(f"{operation} max_output_rows must be a positive integer.")
    if isinstance(max_working_memory_gib, bool) or not isinstance(max_working_memory_gib, (int, float)):
        raise TypeError(f"{operation} max_working_memory_gib must be numeric.")
    max_working_memory_gib = float(max_working_memory_gib)
    if not math.isfinite(max_working_memory_gib) or max_working_memory_gib <= 0.0:
        raise ValueError(f"{operation} max_working_memory_gib must be finite and positive.")

    backend_version = None
    backend_warnings = []
    generated_key = None
    score_frame = None
    method_parameters = {}
    if method in {"aucell", "gsva"}:
        if method == "aucell":
            if isinstance(n_top_features, bool) or not isinstance(n_top_features, int) or n_top_features < 0:
                raise TypeError("AUCell n_top_features must be a nonnegative integer.")
            if n_top_features == 1 or n_top_features > len(feature_names):
                raise ValueError(f"AUCell n_top_features must be 0 or between 2 and {len(feature_names)} inclusive.")
            resolved_n_up = (
                int(np.clip(np.ceil(0.05 * len(feature_names)), a_min=2, a_max=len(feature_names)))
                if n_top_features == 0
                else n_top_features
            )
            requested_output_rows = len(observation_names) * len(retained_sources)
            estimated_working_bytes = requested_output_rows * 8 + len(feature_names) * 8
            generated_key = "score_aucell"
            specific_kwargs = {"n_up": None if n_top_features == 0 else n_top_features}
            required_specific = {"n_up"}
            method_parameters = {
                "n_top_features": int(n_top_features),
                "resolved_n_up": int(resolved_n_up),
                "tie_policy": (
                    "decoupler 2.2 fixed seed-0 feature permutation followed by scipy ordinal ranking; "
                    "ties then follow the permuted feature order"
                ),
            }
        else:
            if kernel not in {"gaussian_normalized", "poisson_counts", "empirical"}:
                raise ValueError("GSVA kernel must be 'gaussian_normalized', 'poisson_counts', or 'empirical'.")
            if not isinstance(maxdiff, bool) or not isinstance(absrnk, bool):
                raise TypeError("GSVA maxdiff and absrnk must be boolean.")
            if isinstance(tau, bool) or not isinstance(tau, (int, float)):
                raise TypeError("GSVA tau must be numeric.")
            tau = float(tau)
            if not math.isfinite(tau) or tau <= 0.0:
                raise ValueError("GSVA tau must be finite and positive.")
            full_sets = [
                source
                for source in retained_sources
                if len(resource["targets_in_universe"][source]) >= len(feature_names)
            ]
            if full_sets:
                raise ValueError(
                    "GSVA resource sets cannot span the complete measured feature universe because the "
                    f"running-sum complement is empty; offending sets={full_sets[:20]}."
                )
            requested_output_rows = len(observation_names) * len(retained_sources)
            dense_conversion_bytes = len(observation_names) * len(feature_names) * 8
            estimated_working_bytes = dense_conversion_bytes
            generated_key = "score_gsva"
            kernel_argument = {
                "gaussian_normalized": "gaussian",
                "poisson_counts": "poisson",
                "empirical": None,
            }[kernel]
            specific_kwargs = {
                "layer": None,
                "kcdf": kernel_argument,
                "maxdiff": maxdiff,
                "absrnk": absrnk,
                "tau": tau,
            }
            required_specific = {"kcdf", "maxdiff", "absrnk", "tau"}
            method_parameters = {
                "kernel": kernel,
                "kcdf": kernel_argument,
                "maxdiff": maxdiff,
                "absrnk": absrnk,
                "tau": tau,
                "dense_conversion_bytes": dense_conversion_bytes,
                "cohort_dependent": True,
            }
        if requested_output_rows > max_output_rows:
            raise ValueError(
                f"{operation} would create {requested_output_rows:,} observation-by-set values, exceeding "
                f"max_output_rows={max_output_rows:,}."
            )
        guard_bytes = int(max_working_memory_gib * 1024**3)
        if estimated_working_bytes > guard_bytes:
            raise ValueError(
                f"{operation} working-memory preflight requires at least {estimated_working_bytes:,} bytes, "
                f"exceeding max_working_memory_gib={max_working_memory_gib:g} ({guard_bytes:,} bytes)."
            )
        try:
            backend_module = importlib.import_module("decoupler")
        except (ImportError, OSError) as exc:
            raise RuntimeError(f"{operation} requires decoupler 2.2, but it is unavailable ({exc}).") from exc
        backend_version = getattr(backend_module, "__version__", None)
        if not isinstance(backend_version, str) or not backend_version.strip():
            try:
                backend_version = importlib.metadata.version("decoupler")
            except importlib.metadata.PackageNotFoundError as exc:
                raise RuntimeError(f"{operation} cannot determine the installed decoupler version.") from exc
        try:
            version_parts = tuple(int(part) for part in backend_version.split(".")[:2])
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{operation} received invalid decoupler version {backend_version!r}.") from exc
        if version_parts != (2, 2):
            raise RuntimeError(
                f"{operation} requires the reviewed decoupler 2.2 public mt API; installed version is {backend_version!r}."
            )
        backend = getattr(getattr(backend_module, "mt", None), method, None)
        if not callable(backend):
            raise RuntimeError(f"{operation} requires callable decoupler.mt.{method}.")
        try:
            public_signature = runtime_inspect.signature(backend)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{operation} could not inspect decoupler.mt.{method}.") from exc
        required_public = {"data", "net", "tmin", "raw", "empty", "bsize", "verbose"}
        if not required_public.issubset(public_signature.parameters):
            raise RuntimeError(
                f"{operation} decoupler.mt.{method} is missing reviewed public parameters: "
                f"{sorted(required_public - set(public_signature.parameters))}."
            )
        public_has_kwargs = any(
            parameter.kind is runtime_inspect.Parameter.VAR_KEYWORD
            for parameter in public_signature.parameters.values()
        )
        missing_specific = required_specific - set(public_signature.parameters)
        if missing_specific:
            implementation = getattr(backend, "func", None)
            if not public_has_kwargs or not callable(implementation):
                raise RuntimeError(
                    f"{operation} decoupler.mt.{method} cannot verify method parameters: {sorted(missing_specific)}."
                )
            try:
                implementation_signature = runtime_inspect.signature(implementation)
            except (TypeError, ValueError) as exc:
                raise RuntimeError(f"{operation} could not inspect decoupler.mt.{method}.func.") from exc
            if not missing_specific.issubset(implementation_signature.parameters):
                raise RuntimeError(
                    f"{operation} decoupler.mt.{method}.func is missing parameters: "
                    f"{sorted(missing_specific - set(implementation_signature.parameters))}."
                )
        # Scientific workspace: decoupler writes result slots and may coerce its matrix input.
        work_matrix = matrix.copy()
        work = ad.AnnData(X=work_matrix)
        work.obs_names = pd.Index(observation_names, dtype="object")
        work.var_names = pd.Index(feature_names, dtype="object")
        backend(
            data=work,
            net=resource["network"],
            tmin=min_targets,
            raw=False,
            empty=False,
            bsize=batch_size,
            verbose=False,
            **specific_kwargs,
        )
        if generated_key not in work.obsm:
            raise RuntimeError(f"{operation} decoupler did not create obsm[{generated_key!r}].")
        raw_scores = work.obsm[generated_key]
        if not isinstance(raw_scores, pd.DataFrame):
            raise RuntimeError(f"{operation} decoupler output must be a named pandas DataFrame.")
        if not raw_scores.index.equals(pd.Index(observation_names)):
            raise RuntimeError(f"{operation} decoupler output changed the observation identity/order.")
        if len(raw_scores.columns) != len(set(raw_scores.columns)):
            raise RuntimeError(f"{operation} decoupler output contains duplicate score columns.")
        if set(raw_scores.columns) != set(retained_sources):
            raise RuntimeError(
                f"{operation} decoupler output set family differs from the independently retained resource: "
                f"missing={sorted(set(retained_sources) - set(raw_scores.columns))[:20]}, "
                f"unknown={sorted(set(raw_scores.columns) - set(retained_sources))[:20]}."
            )
        # Scientific result materialization: detach the canonical column order from the backend workspace.
        score_frame = raw_scores.loc[:, retained_sources].copy()
        try:
            score_values = score_frame.to_numpy(dtype=float, copy=True)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{operation} decoupler output scores must be numeric.") from exc
        if not bool(np.isfinite(score_values).all()):
            raise RuntimeError(f"{operation} decoupler output contains non-finite scores.")
        tolerance = 1e-12
        lower, upper = (0.0, 1.0) if method == "aucell" else (-1.0, 1.0)
        if bool((score_values < lower - tolerance).any()) or bool((score_values > upper + tolerance).any()):
            raise RuntimeError(f"{operation} decoupler scores fall outside the reviewed [{lower:g}, {upper:g}] range.")
        method_parameters.update(
            {
                "min_targets": min_targets,
                "batch_size": batch_size,
                "max_output_rows": max_output_rows,
                "max_working_memory_gib": max_working_memory_gib,
                "requested_output_values": requested_output_rows,
                "estimated_guarded_working_bytes": estimated_working_bytes,
                "decoupler_fixed_arguments": {
                    "raw": False,
                    "empty": False,
                    "verbose": False,
                },
            }
        )
    else:
        if not isinstance(panel, str) or not panel or panel != panel.strip():
            raise ValueError("Gene panel identifier must be a nonblank, whitespace-canonical string.")
        if panel not in resource["source_order"]:
            raise ValueError(f"Gene panel {panel!r} is not present in the supplied resource.")
        matched_genes = list(resource["targets_in_universe"][panel])
        requested_genes = list(resource["targets_before"][panel])
        if not matched_genes:
            raise ValueError(f"Gene panel {panel!r} has no exact match in the selected expression universe.")
        if isinstance(ctrl_size, bool) or not isinstance(ctrl_size, int) or ctrl_size < 0:
            raise TypeError("Gene panel ctrl_size must be a nonnegative integer.")
        resolved_ctrl_size = len(matched_genes) if ctrl_size == 0 else ctrl_size
        if resolved_ctrl_size < 1:
            raise ValueError("Gene panel resolved control size must be positive.")
        if isinstance(n_bins, bool) or not isinstance(n_bins, int) or n_bins < 2:
            raise TypeError("Gene panel n_bins must be an integer of at least 2.")
        if isinstance(random_seed, bool) or not isinstance(random_seed, int) or not 0 <= random_seed <= 2**31 - 1:
            raise TypeError("Gene panel random_seed must be an integer in [0, 2^31-1].")
        try:
            scanpy = importlib.import_module("scanpy")
        except (ImportError, OSError) as exc:
            raise RuntimeError(f"Gene panel scoring requires Scanpy 1.12, but it is unavailable ({exc}).") from exc
        try:
            backend_version = importlib.metadata.version("scanpy")
        except importlib.metadata.PackageNotFoundError as exc:
            raise RuntimeError("Gene panel scoring cannot determine the installed Scanpy version.") from exc
        try:
            scanpy_parts = tuple(int(part) for part in backend_version.split(".")[:3])
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"Gene panel scoring received invalid Scanpy version {backend_version!r}.") from exc
        if scanpy_parts[:2] != (1, 12) or scanpy_parts[2] < 3:
            raise RuntimeError(
                f"Gene panel scoring requires the reviewed Scanpy >=1.12.3,<1.13 interface; installed {backend_version!r}."
            )
        backend = getattr(getattr(scanpy, "tl", None), "score_genes", None)
        if not callable(backend):
            raise RuntimeError("Gene panel scoring requires callable scanpy.tl.score_genes.")
        required_parameters = {
            "adata",
            "gene_list",
            "ctrl_as_ref",
            "ctrl_size",
            "gene_pool",
            "n_bins",
            "score_name",
            "random_state",
            "copy",
            "use_raw",
            "layer",
        }
        try:
            signature = runtime_inspect.signature(backend)
        except (TypeError, ValueError) as exc:
            raise RuntimeError("Gene panel scoring could not inspect scanpy.tl.score_genes.") from exc
        if not required_parameters.issubset(signature.parameters):
            raise RuntimeError(
                "Gene panel scoring Scanpy interface is missing reviewed parameters: "
                f"{sorted(required_parameters - set(signature.parameters))}."
            )
        # Scientific workspace: Scanpy score_genes mutates obs and may inspect/coerce X.
        work = ad.AnnData(X=matrix.copy())
        work.obs_names = pd.Index(observation_names, dtype="object")
        work.var_names = pd.Index(feature_names, dtype="object")
        temporary_key = "__openbio_gene_panel_score__"
        global_rng_before = np.random.get_state()
        caught = []
        backend_changed_global_rng = False
        try:
            with runtime_warnings.catch_warnings(record=True) as caught:
                runtime_warnings.simplefilter("always")
                backend(
                    work,
                    gene_list=matched_genes,
                    ctrl_as_ref=False,
                    ctrl_size=resolved_ctrl_size,
                    gene_pool=work.var_names,
                    n_bins=n_bins,
                    score_name=temporary_key,
                    random_state=random_seed,
                    copy=False,
                    use_raw=False,
                    layer=None,
                )
        finally:
            global_rng_after = np.random.get_state()
            backend_changed_global_rng = not (
                global_rng_before[0] == global_rng_after[0]
                and np.array_equal(global_rng_before[1], global_rng_after[1])
                and global_rng_before[2:] == global_rng_after[2:]
            )
            np.random.set_state(global_rng_before)
        global_rng_restored = True
        backend_warnings = [str(item.message) for item in caught]
        if temporary_key not in work.obs:
            raise RuntimeError("Gene panel scoring Scanpy backend did not create the requested score column.")
        series = work.obs[temporary_key]
        if not series.index.equals(pd.Index(observation_names)):
            raise RuntimeError("Gene panel scoring backend changed the observation identity/order.")
        try:
            score_values = pd.to_numeric(series, errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise RuntimeError("Gene panel scoring backend returned a nonnumeric score.") from exc
        if not bool(np.isfinite(score_values).all()):
            raise RuntimeError("Gene panel scoring backend returned non-finite scores.")
        score_frame = pd.DataFrame({output_key: score_values}, index=pd.Index(observation_names))
        method_parameters = {
            "panel": panel,
            "requested_gene_count": len(requested_genes),
            "matched_gene_count": len(matched_genes),
            "matched_fraction": float(len(matched_genes) / len(requested_genes)),
            "requested_genes_preview": requested_genes[:200],
            "requested_genes_preview_truncated": len(requested_genes) > 200,
            "matched_genes_preview": matched_genes[:200],
            "matched_genes_preview_truncated": len(matched_genes) > 200,
            "excluded_genes_preview": [gene for gene in requested_genes if gene not in set(matched_genes)][:200],
            "excluded_gene_count": len(requested_genes) - len(matched_genes),
            "ctrl_size": ctrl_size,
            "resolved_ctrl_size": resolved_ctrl_size,
            "n_bins": n_bins,
            "random_seed": random_seed,
            "gene_pool": "complete selected expression feature universe",
            "gene_pool_count": len(feature_names),
            "scanpy_fixed_arguments": {
                "ctrl_as_ref": False,
                "copy": False,
                "use_raw": False,
                "layer": None,
            },
            "numpy_global_rng_state_restored": global_rng_restored,
            "scanpy_changed_numpy_global_rng_state_inside_isolation": backend_changed_global_rng,
        }
        if backend_warnings:
            warnings.extend(f"Scanpy score_genes: {message}" for message in backend_warnings[:20])
            if len(backend_warnings) > 20:
                warnings.append(
                    f"{len(backend_warnings) - 20} additional Scanpy warnings were omitted from this preview."
                )

    if score_frame is None:
        raise RuntimeError(f"{operation} failed to construct a canonical score frame.")
    score_identity = score_frame_fingerprint(score_frame, np=np, pd=pd)
    score_values = score_frame.to_numpy(dtype=float, copy=True)
    flat_values = score_values.ravel()
    quantiles = np.quantile(flat_values, [0.0, 0.25, 0.5, 0.75, 1.0])
    score_distribution_preview_limit = 200
    per_score_distributions = []
    for score_name in score_frame.columns[:score_distribution_preview_limit]:
        values = score_frame[score_name].to_numpy(dtype=float)
        score_quantiles = np.quantile(values, [0.0, 0.25, 0.5, 0.75, 1.0])
        per_score_distributions.append(
            {
                "score": str(score_name),
                "min": float(score_quantiles[0]),
                "q1": float(score_quantiles[1]),
                "median": float(score_quantiles[2]),
                "mean": float(values.mean()),
                "q3": float(score_quantiles[3]),
                "max": float(score_quantiles[4]),
            }
        )
    score_means = score_frame.mean(axis=0)
    highest = score_means.sort_values(ascending=False, kind="mergesort").head(5)
    lowest = score_means.sort_values(ascending=True, kind="mergesort").head(5)
    if resource["excluded_sources"]:
        warnings.append(
            f"Excluded {len(resource['excluded_sources'])} resource sets below the post-intersection target threshold."
        )
    if collision:
        warnings.append(f"Explicitly replaced existing {storage}[{output_key!r}] and/or its score provenance.")

    references = []
    if method == "aucell":
        references.extend(
            [
                {
                    "citation": (
                        "Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. "
                        "Nature Methods. 2017;14:1083-1086."
                    ),
                    "url": "https://doi.org/10.1038/nmeth.4463",
                    "kind": "method",
                    "doi": "10.1038/nmeth.4463",
                },
                {
                    "citation": (
                        "Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer "
                        "biological activities from omics data. Bioinformatics Advances. 2022;2:vbac016."
                    ),
                    "url": "https://doi.org/10.1093/bioadv/vbac016",
                    "kind": "software",
                    "doi": "10.1093/bioadv/vbac016",
                },
            ]
        )
    elif method == "gsva":
        references.extend(
            [
                {
                    "citation": (
                        "Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for "
                        "microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7."
                    ),
                    "url": "https://doi.org/10.1186/1471-2105-14-7",
                    "kind": "method",
                    "doi": "10.1186/1471-2105-14-7",
                },
                {
                    "citation": (
                        "Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer "
                        "biological activities from omics data. Bioinformatics Advances. 2022;2:vbac016."
                    ),
                    "url": "https://doi.org/10.1093/bioadv/vbac016",
                    "kind": "software",
                    "doi": "10.1093/bioadv/vbac016",
                },
            ]
        )
    else:
        references.extend(
            [
                {
                    "citation": (
                        "Tirosh I, et al. Dissecting the multicellular ecosystem of metastatic melanoma "
                        "by single-cell RNA-seq. Science. 2016;352:189-196."
                    ),
                    "url": "https://doi.org/10.1126/science.aad0501",
                    "kind": "method",
                    "doi": "10.1126/science.aad0501",
                },
                {
                    "citation": (
                        "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression "
                        "data analysis. Genome Biology. 2018;19:15."
                    ),
                    "url": "https://doi.org/10.1186/s13059-017-1382-0",
                    "kind": "software",
                    "doi": "10.1186/s13059-017-1382-0",
                },
            ]
        )
    references.extend(
        [
            {
                "citation": "Virshup I, et al. anndata: Annotated data. JOSS. 2024;9:4371.",
                "url": "https://doi.org/10.21105/joss.04371",
                "kind": "software",
                "doi": "10.21105/joss.04371",
            },
            {
                "citation": (
                    f"{resource['metadata'].get('name', 'Gene-set resource')} "
                    f"{resource['metadata'].get('version', '')} "
                    f"({resource['metadata'].get('date', 'date not supplied')}): "
                    f"{resource['metadata'].get('citation', 'citation not supplied')}"
                ),
                "url": f"urn:sha256:{resource['sha256']}",
                "kind": "resource",
                "doi": None,
            },
        ]
    )
    parameters = {
        "gene_sets_file": requested_resource_path,
        "resource_metadata": resource["metadata"],
        "source": source_kind,
        "layer_name": layer_name if source_kind == "layer" else None,
        "source_column": source_column,
        "target_column": target_column,
        "output_key": output_key,
        "overwrite_existing": overwrite_existing,
        **method_parameters,
    }
    artifact = build_score_artifact(
        frame=score_frame,
        feature_names=feature_names,
        producer_node=producer_node,
        method=method_display,
        storage=storage,
        score_key=output_key,
        resource=resource_provenance,
        expression=expression_provenance,
        parameters=parameters,
        references=references,
        np=np,
        pd=pd,
    )
    output = adata
    if method == "panel":
        output.obs[output_key] = score_frame[output_key]
    else:
        output.obsm[output_key] = score_frame
    store_score_artifact(output, artifact)
    stored_frame = (
        pd.DataFrame({output_key: output.obs[output_key]}, index=output.obs_names.copy())
        if method == "panel"
        else output.obsm[output_key]
    )
    stored_identity = score_frame_fingerprint(stored_frame, np=np, pd=pd)
    if stored_identity != score_identity:
        raise RuntimeError(f"{operation} stored score data differs from the validated backend result.")
    if not output.obs_names.equals(adata.obs_names) or not output.var_names.equals(adata.var_names):
        raise RuntimeError(f"{operation} changed the AnnData observation or active feature axis.")

    if method == "aucell":
        methods = (
            "decoupler 2.2 mt.aucell ranked each observation's explicitly selected expression and "
            "calculated an unweighted recovery-curve AUC for every supplied resource set meeting the exact "
            "post-intersection target threshold. No hypothesis test or p-value was computed."
        )
        limitations = [
            "AUCell scores are descriptive per observation and do not constitute Condition inference or p-values.",
            "Ties follow decoupler's fixed seed-0 feature permutation and ordinal feature order.",
            "Cells from one Sample are not biological replicates; downstream Condition inference must be Sample-level within a population.",
            "Gene identifiers were matched exactly without case folding, synonym translation, or online lookup.",
        ]
    elif method == "gsva":
        methods = (
            f"decoupler 2.2 mt.gsva computed {kernel!r} cohort-dependent, unweighted observation-level "
            f"gene-set scores with maxdiff={maxdiff}, absrnk={absrnk}, and tau={tau:g}. The complete selected "
            "observation cohort defined the density transformation; no hypothesis test or p-value was computed."
        )
        limitations = [
            "GSVA scores are relative to the observations included in this run and can change when the cohort is changed.",
            "GSVA scores are descriptive per observation and do not constitute Condition inference or p-values.",
            "Cells from one Sample are not biological replicates; downstream Condition inference must be Sample-level within a population.",
            "The in-memory decoupler 2.2 GSVA path can densify the complete expression matrix despite the public batch-size parameter.",
            "Gene identifiers were matched exactly without case folding, synonym translation, or online lookup.",
        ]
    else:
        methods = (
            f"Scanpy {backend_version} tl.score_genes calculated one predeclared panel score as mean matched-panel "
            "expression minus a seeded expression-bin-matched control mean. ctrl_as_ref=False, the complete selected "
            "feature universe as gene_pool, use_raw=False, copy=False, and layer=None were fixed explicitly."
        )
        limitations = [
            "The panel score is descriptive per cell and is not a formal enrichment p-value or causal activity state.",
            "Control-gene selection depends on the selected expression cohort and complete measured feature pool.",
            "Cells from one Sample are not biological replicates; downstream Condition inference must be Sample-level within a population.",
            "Gene identifiers were matched exactly without case folding, synonym translation, or online lookup.",
        ]
    mean_preview = ", ".join(f"{name} ({float(value):.3g})" for name, value in highest.items())
    results = (
        f"Scored {len(observation_names):,} observations across {score_frame.shape[1]:,} named score columns "
        f"using {len(feature_names):,} measured features and resource {resource['metadata'].get('name', requested_resource_path)!r} "
        f"version {resource['metadata'].get('version', 'not supplied')!r}. Scores ranged from {quantiles[0]:.4g} to "
        f"{quantiles[4]:.4g}; leading columns by mean score were {mean_preview or 'none'}. "
        "These are observation-level descriptive results, not a replicate-aware Condition contrast."
    )

    def package_version(distribution):
        try:
            return importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            return "not-installed"

    software_versions = {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
        "anndata": package_version("anndata"),
        "scanpy": package_version("scanpy"),
        "numpy": package_version("numpy"),
        "pandas": package_version("pandas"),
        "scipy": str(scipy.__version__),
    }
    if method in {"aucell", "gsva"}:
        software_versions["decoupler"] = str(backend_version)
    key_results = {
        "method": method_display,
        "statistical_test": None,
        "p_values_produced": False,
        "observation_count": len(observation_names),
        "feature_count": len(feature_names),
        "score_count": score_frame.shape[1],
        "score_names_preview": list(score_frame.columns[:200]),
        "score_names_preview_truncated": score_frame.shape[1] > 200,
        "score_range": {"min": float(quantiles[0]), "max": float(quantiles[4])},
        "score_distribution": {
            "min": float(quantiles[0]),
            "q1": float(quantiles[1]),
            "median": float(quantiles[2]),
            "mean": float(flat_values.mean()),
            "q3": float(quantiles[3]),
            "max": float(quantiles[4]),
        },
        "per_score_distributions_preview": per_score_distributions,
        "per_score_distributions_count": score_frame.shape[1],
        "per_score_distributions_preview_truncated": score_frame.shape[1] > score_distribution_preview_limit,
        "highest_mean_scores": [{"score": str(name), "mean": float(value)} for name, value in highest.items()],
        "lowest_mean_scores": [{"score": str(name), "mean": float(value)} for name, value in lowest.items()],
        "resource": resource_provenance,
        "expression": expression_provenance,
        "score_artifact_fingerprint_sha256": artifact["artifact_fingerprint_sha256"],
        "score_content_fingerprint_sha256": artifact["score_content_sha256"],
        "backend_warnings": backend_warnings[:20],
        "backend_warnings_count": len(backend_warnings),
    }
    summary = {
        "schema_version": 1,
        "node_id": producer_node,
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": limitations,
        "references": references,
        "software_versions": software_versions,
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return output, summary


def run_aucell_scores(
    adata: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str,
    source_kind: str,
    layer_name: str | None,
    source_column: str,
    target_column: str,
    min_targets: int,
    n_top_features: int,
    batch_size: int,
    max_output_rows: int,
    max_working_memory_gib: float,
    output_key: str,
    overwrite_existing: bool,
    expected_resource_sha256: str | None = None,
    openbio_version: str = "unknown",
) -> tuple[Any, dict[str, Any]]:
    return _standalone_score_gene_sets(
        adata,
        method="aucell",
        gene_sets_path=gene_sets_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        source_kind=source_kind,
        layer_name=layer_name,
        source_column=source_column,
        target_column=target_column,
        min_targets=min_targets,
        n_top_features=n_top_features,
        batch_size=batch_size,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        output_key=output_key,
        overwrite_existing=overwrite_existing,
        expected_resource_sha256=expected_resource_sha256,
        openbio_version=openbio_version,
    )


def run_gsva_scores(
    adata: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str,
    source_kind: str,
    layer_name: str | None,
    source_column: str,
    target_column: str,
    min_targets: int,
    kernel: str,
    maxdiff: bool,
    absrnk: bool,
    tau: float,
    batch_size: int,
    max_output_rows: int,
    max_working_memory_gib: float,
    output_key: str,
    overwrite_existing: bool,
    expected_resource_sha256: str | None = None,
    openbio_version: str = "unknown",
) -> tuple[Any, dict[str, Any]]:
    return _standalone_score_gene_sets(
        adata,
        method="gsva",
        gene_sets_path=gene_sets_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        source_kind=source_kind,
        layer_name=layer_name,
        source_column=source_column,
        target_column=target_column,
        min_targets=min_targets,
        kernel=kernel,
        maxdiff=maxdiff,
        absrnk=absrnk,
        tau=tau,
        batch_size=batch_size,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        output_key=output_key,
        overwrite_existing=overwrite_existing,
        expected_resource_sha256=expected_resource_sha256,
        openbio_version=openbio_version,
    )


def run_gene_panel_score(
    adata: Any,
    *,
    gene_sets_path: str,
    requested_resource_path: str,
    resource_metadata_json: str,
    panel: str,
    source_kind: str,
    layer_name: str | None,
    source_column: str,
    target_column: str,
    output_key: str,
    ctrl_size: int,
    n_bins: int,
    random_seed: int,
    overwrite_existing: bool,
    expected_resource_sha256: str | None = None,
    openbio_version: str = "unknown",
) -> tuple[Any, dict[str, Any]]:
    return _standalone_score_gene_sets(
        adata,
        method="panel",
        gene_sets_path=gene_sets_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        source_kind=source_kind,
        layer_name=layer_name,
        source_column=source_column,
        target_column=target_column,
        min_targets=1,
        panel=panel,
        output_key=output_key,
        ctrl_size=ctrl_size,
        n_bins=n_bins,
        random_seed=random_seed,
        overwrite_existing=overwrite_existing,
        expected_resource_sha256=expected_resource_sha256,
        openbio_version=openbio_version,
    )


def gene_set_scoring_code(*, function_name: str, parameters: dict[str, Any]) -> str:
    if function_name not in {"score_aucell", "score_gsva", "score_gene_panel"}:
        raise ValueError(f"Unsupported generated scoring function name: {function_name!r}.")
    helpers = "\n\n".join(
        textwrap.dedent(inspect.getsource(helper)).strip()
        for helper in (
            _canonical_json_sha256,
            _canonical_string_axis,
            score_frame_fingerprint,
            build_score_artifact,
            store_score_artifact,
        )
    )
    implementation = textwrap.dedent(inspect.getsource(_standalone_score_gene_sets)).strip()
    arguments = ",\n".join(f"        {name}={value!r}" for name, value in parameters.items())
    return f"""from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any

SCORE_ARTIFACTS_KEY = {SCORE_ARTIFACTS_KEY!r}
SCORE_ARTIFACT_SCHEMA_VERSION = {SCORE_ARTIFACT_SCHEMA_VERSION!r}

{gene_set_resource_code()}

load_gene_set_resource = _standalone_load_gene_set_resource

{helpers}

{implementation}


def {function_name}(adata):
    return _standalone_score_gene_sets(
        adata,
{arguments},
    )
"""


__all__ = [
    "gene_set_resource_cache_fingerprint",
    "gene_set_scoring_code",
    "run_aucell_scores",
    "run_gene_panel_score",
    "run_gsva_scores",
]
