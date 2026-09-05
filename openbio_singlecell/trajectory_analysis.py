from __future__ import annotations

import hashlib
import inspect
import json
import math
import textwrap
import threading
import warnings
from collections.abc import Mapping, Sequence
from typing import Any

from . import PLUGIN_VERSION
from .analysis_reporting import _package_version, collect_software_versions

DIFFMAP_NODE_ID = "OpenBioSingleCellDiffusionMap"
PAGA_NODE_ID = "OpenBioSingleCellPAGA"
DPT_NODE_ID = "OpenBioSingleCellDPT"

_TA_DIFFMAP_SCHEMA = "openbio-singlecell/diffusion-map/v1"
_TA_PAGA_SCHEMA = "openbio-singlecell/paga/v1"
_TA_DPT_SCHEMA = "openbio-singlecell/dpt/v1"
_TA_DIFFMAP_PROVENANCE_KEY = "openbio_diffusion_map"
_TA_PAGA_PROVENANCE_KEY = "openbio_paga"
_TA_DPT_PROVENANCE_KEY = "openbio_dpt"
_TA_PAGA_GROUP_WARNING_THRESHOLD = 256
_TA_MAX_PAGA_SUMMARY_EDGES = 2_000
_TA_MAX_GRAPH_COMPONENT_SIZES = 100
_TA_MAX_DIFFMAP_COMPONENT_SUMMARIES = 100
_TA_RUNTIME_LOCK = threading.RLock()


def _ta_import_science(operation):
    import importlib
    import re
    from importlib import metadata as package_metadata

    try:
        scanpy = importlib.import_module("scanpy")
        numpy = importlib.import_module("numpy")
        pandas = importlib.import_module("pandas")
        scipy_sparse = importlib.import_module("scipy.sparse")
        scipy_csgraph = importlib.import_module("scipy.sparse.csgraph")
    except ImportError as error:
        raise RuntimeError(f"{operation} requires Scanpy, AnnData, NumPy, pandas, and SciPy.") from error
    try:
        scanpy_version = package_metadata.version("scanpy")
    except package_metadata.PackageNotFoundError as error:
        raise RuntimeError(f"{operation} could not determine the installed Scanpy version.") from error
    match = re.match(r"^(\d+)\.(\d+)\.(\d+)", scanpy_version)
    if match is None:
        raise RuntimeError(f"{operation} could not parse Scanpy version {scanpy_version!r}.")
    release = tuple(int(value) for value in match.groups())
    if release < (1, 12, 3) or release >= (1, 13, 0):
        raise RuntimeError(f"{operation} supports the audited Scanpy >=1.12.3,<1.13 API; found {scanpy_version}.")
    return scanpy, numpy, pandas, scipy_sparse, scipy_csgraph, scanpy_version


def _ta_require_signature(callable_object, *, operation, expected_names):
    try:
        signature = inspect.signature(callable_object)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"{operation} could not inspect the Scanpy backend signature.") from error
    actual_names = tuple(signature.parameters)
    if actual_names != tuple(expected_names):
        raise RuntimeError(
            f"{operation} requires the audited Scanpy signature {tuple(expected_names)!r}; found {actual_names!r}."
        )


def _ta_clean_text(value, *, label):
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    if not value or value != value.strip():
        raise ValueError(f"{label} must be nonempty and free of surrounding whitespace.")
    return value


def _ta_integer(value, *, label, minimum, maximum):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer.")
    if value < minimum or value > maximum:
        raise ValueError(f"{label} must be between {minimum} and {maximum}.")
    return int(value)


def _ta_boolean(value, *, label):
    if not isinstance(value, bool):
        raise TypeError(f"{label} must be a boolean.")
    return value


def _ta_json_value(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not value == value or value in (float("inf"), float("-inf")):
            raise ValueError("A strict JSON value cannot contain NaN or infinity.")
        return value
    if hasattr(value, "item"):
        try:
            return _ta_json_value(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _ta_json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_ta_json_value(item) for item in value]
    raise TypeError(f"Unsupported strict JSON value type: {type(value).__name__}.")


def _ta_obs_names(adata, *, operation, minimum_cells):
    try:
        if bool(adata.isbacked):
            raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
        n_obs = int(adata.n_obs)
        obs_names = list(adata.obs_names)
    except AttributeError as error:
        raise TypeError(f"{operation} requires an AnnData input.") from error
    if n_obs < minimum_cells:
        raise ValueError(f"{operation} requires at least {minimum_cells} observations.")
    if len(obs_names) != n_obs or not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires a unique, observation-aligned cell axis.")
    for value in obs_names:
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{operation} requires canonical nonempty string observation identifiers.")
    return tuple(obs_names)


def _ta_hash_text_records(records):
    digest = hashlib.sha256()
    for record in records:
        encoded = str(record).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little", signed=False))
        digest.update(encoded)
    return digest.hexdigest()


def _ta_hash_sparse(matrix, *, numpy, preserve_zeros):
    canonical = matrix.tocsr(copy=True)
    canonical.sum_duplicates()
    canonical.sort_indices()
    if not preserve_zeros:
        canonical.eliminate_zeros()
        canonical.sort_indices()
    digest = hashlib.sha256()
    digest.update(b"csr-f8-preserve" if preserve_zeros else b"csr-f8-eliminate")
    digest.update(numpy.asarray(canonical.shape, dtype="<i8").tobytes())
    digest.update(numpy.asarray(canonical.indptr, dtype="<i8").tobytes())
    digest.update(numpy.asarray(canonical.indices, dtype="<i8").tobytes())
    digest.update(numpy.asarray(canonical.data, dtype="<f8").tobytes())
    return digest.hexdigest()


def _ta_hash_array(array, *, numpy):
    canonical = numpy.ascontiguousarray(array, dtype="<f8")
    digest = hashlib.sha256()
    digest.update(b"dense-f8")
    digest.update(numpy.asarray(canonical.shape, dtype="<i8").tobytes())
    digest.update(canonical.tobytes(order="C"))
    return digest.hexdigest()


def _ta_sparse_graph_matrix(matrix, *, n_obs, operation, name, numpy, scipy_sparse, preserve_zeros):
    if not scipy_sparse.issparse(matrix):
        raise TypeError(f"{operation} {name} matrix must be a SciPy sparse matrix.")
    if tuple(matrix.shape) != (n_obs, n_obs):
        raise ValueError(f"{operation} {name} matrix must have shape (n_obs, n_obs).")
    canonical = matrix.tocsr(copy=True)
    canonical.sum_duplicates()
    canonical.sort_indices()
    values = numpy.asarray(canonical.data)
    if values.size and (not numpy.issubdtype(values.dtype, numpy.number) or numpy.iscomplexobj(values)):
        raise TypeError(f"{operation} {name} matrix must contain real numeric values.")
    if values.size and not bool(numpy.isfinite(values).all()):
        raise ValueError(f"{operation} {name} matrix contains non-finite values.")
    if values.size and bool((values < 0).any()):
        raise ValueError(f"{operation} {name} matrix contains negative values.")
    coordinates = canonical.tocoo(copy=False)
    diagonal = coordinates.row == coordinates.col
    if bool(numpy.any(diagonal)):
        if not bool(numpy.allclose(coordinates.data[diagonal], 0.0, rtol=0.0, atol=1e-12)):
            raise ValueError(f"{operation} {name} matrix must have a zero diagonal.")
        off_diagonal = ~diagonal
        canonical = scipy_sparse.csr_matrix(
            (
                coordinates.data[off_diagonal],
                (coordinates.row[off_diagonal], coordinates.col[off_diagonal]),
            ),
            shape=canonical.shape,
        )
        canonical.sum_duplicates()
        canonical.sort_indices()
    if not preserve_zeros:
        canonical.eliminate_zeros()
        canonical.sort_indices()
    if canonical.nnz == 0:
        raise ValueError(f"{operation} {name} matrix contains no off-diagonal graph edges.")
    return canonical


def _ta_graph(adata, *, neighbors_key, operation, numpy, scipy_sparse, scipy_csgraph):
    obs_names = _ta_obs_names(adata, operation=operation, minimum_cells=2)
    neighbors_key = _ta_clean_text(neighbors_key, label=f"{operation} neighbors_key")
    metadata_record = adata.uns.get(neighbors_key)
    if not isinstance(metadata_record, Mapping):
        raise ValueError(f"{operation} neighbor metadata not found in uns[{neighbors_key!r}].")
    connectivity_key = metadata_record.get("connectivities_key")
    distance_key = metadata_record.get("distances_key")
    connectivity_key = _ta_clean_text(connectivity_key, label=f"{operation} connectivities_key")
    distance_key = _ta_clean_text(distance_key, label=f"{operation} distances_key")
    if connectivity_key == distance_key:
        raise ValueError(f"{operation} distance and connectivity pointers must name different matrices.")
    if connectivity_key not in adata.obsp:
        raise ValueError(f"{operation} connectivity matrix not found in obsp[{connectivity_key!r}].")
    if distance_key not in adata.obsp:
        raise ValueError(f"{operation} distance matrix not found in obsp[{distance_key!r}].")
    connectivities = _ta_sparse_graph_matrix(
        adata.obsp[connectivity_key],
        n_obs=len(obs_names),
        operation=operation,
        name="connectivity",
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        preserve_zeros=False,
    )
    delta = connectivities - connectivities.T
    delta.eliminate_zeros()
    if delta.nnz and float(numpy.max(numpy.abs(delta.data))) > 1e-8:
        raise ValueError(f"{operation} requires a symmetric undirected connectivity matrix.")
    connectivities = ((connectivities + connectivities.T) * 0.5).tocsr()
    connectivities.eliminate_zeros()
    connectivities.sort_indices()
    distances = _ta_sparse_graph_matrix(
        adata.obsp[distance_key],
        n_obs=len(obs_names),
        operation=operation,
        name="distance",
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        preserve_zeros=True,
    )
    n_components, component_labels = scipy_csgraph.connected_components(connectivities, directed=False)
    component_counts = numpy.bincount(component_labels, minlength=n_components)
    degree = numpy.asarray((connectivities > 0).sum(axis=1), dtype=int).ravel()
    weighted_degree = numpy.asarray(connectivities.sum(axis=1), dtype=float).ravel()
    parameters = metadata_record.get("params")
    if parameters is None:
        parameters = {}
    if not isinstance(parameters, Mapping):
        raise TypeError(f"{operation} uns[{neighbors_key!r}]['params'] must be a mapping when present.")
    parameters = _ta_json_value(dict(parameters))
    identity = {
        "neighbors_key": neighbors_key,
        "connectivities_key": connectivity_key,
        "distances_key": distance_key,
        "observation_fingerprint_sha256": _ta_hash_text_records(obs_names),
        "connectivities_fingerprint_sha256": _ta_hash_sparse(connectivities, numpy=numpy, preserve_zeros=False),
        "distances_fingerprint_sha256": _ta_hash_sparse(distances, numpy=numpy, preserve_zeros=True),
        "neighbor_parameters": parameters,
    }
    graph_fingerprint = hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")
    ).hexdigest()
    possible_edges = len(obs_names) * (len(obs_names) - 1) // 2
    return {
        "neighbors_key": neighbors_key,
        "connectivities_key": connectivity_key,
        "distances_key": distance_key,
        "connectivities": connectivities,
        "distances": distances,
        "obs_names": obs_names,
        "observation_fingerprint_sha256": identity["observation_fingerprint_sha256"],
        "graph_fingerprint_sha256": graph_fingerprint,
        "component_labels": component_labels,
        "component_sizes": tuple(int(value) for value in sorted(component_counts, reverse=True)),
        "degree": degree,
        "weighted_degree": weighted_degree,
        "positive_undirected_edges": int(connectivities.nnz // 2),
        "distance_directed_edges": int(distances.nnz),
        "density": float((connectivities.nnz // 2) / possible_edges) if possible_edges else 0.0,
        "neighbor_parameters": parameters,
    }


def _ta_assert_exact_graph_storage(adata, *, graph, operation, numpy, scipy_sparse):
    specifications = (
        ("connectivities", graph["connectivities_key"], False),
        ("distances", graph["distances_key"], True),
    )
    for graph_field, key, preserve_zeros in specifications:
        try:
            stored = _ta_sparse_graph_matrix(
                adata.obsp[key],
                n_obs=len(graph["obs_names"]),
                operation=operation,
                name=graph_field,
                numpy=numpy,
                scipy_sparse=scipy_sparse,
                preserve_zeros=preserve_zeros,
            )
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(f"{operation} backend corrupted the exact canonical named graph.") from error
        if _ta_hash_sparse(stored, numpy=numpy, preserve_zeros=preserve_zeros) != _ta_hash_sparse(
            graph[graph_field], numpy=numpy, preserve_zeros=preserve_zeros
        ):
            raise RuntimeError(f"{operation} backend modified the exact canonical named graph matrices.")


def _ta_numeric_summary(values, *, numpy):
    array = numpy.asarray(values, dtype=float).ravel()
    if not array.size or not bool(numpy.isfinite(array).all()):
        raise ValueError("Numeric summary requires at least one finite value.")
    quantiles = numpy.quantile(array, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "n": int(array.size),
        "min": float(quantiles[0]),
        "q1": float(quantiles[1]),
        "median": float(quantiles[2]),
        "mean": float(array.mean()),
        "q3": float(quantiles[3]),
        "max": float(quantiles[4]),
    }


def _ta_graph_diagnostics(graph, *, numpy):
    isolated = numpy.flatnonzero(graph["degree"] == 0)
    component_count = len(graph["component_sizes"])
    return {
        "neighbors_key": graph["neighbors_key"],
        "connectivities_key": graph["connectivities_key"],
        "distances_key": graph["distances_key"],
        "graph_fingerprint_sha256": graph["graph_fingerprint_sha256"],
        "observation_fingerprint_sha256": graph["observation_fingerprint_sha256"],
        "shape": [len(graph["obs_names"]), len(graph["obs_names"])],
        "positive_undirected_connectivity_edges": graph["positive_undirected_edges"],
        "stored_directed_distance_edges": graph["distance_directed_edges"],
        "connectivity_density": graph["density"],
        "connected_components": component_count,
        "component_sizes": list(graph["component_sizes"][:_TA_MAX_GRAPH_COMPONENT_SIZES]),
        "component_sizes_truncated": component_count > _TA_MAX_GRAPH_COMPONENT_SIZES,
        "component_sizes_limit": _TA_MAX_GRAPH_COMPONENT_SIZES,
        "isolated_observation_count": int(isolated.size),
        "isolated_observation_examples": [graph["obs_names"][int(index)] for index in isolated[:20]],
        "degree": _ta_numeric_summary(graph["degree"], numpy=numpy),
        "weighted_degree": _ta_numeric_summary(graph["weighted_degree"], numpy=numpy),
        "neighbor_parameters": graph["neighbor_parameters"],
    }


def _ta_warning_messages(records):
    messages = []
    for record in records:
        message = " ".join(str(record.message).split())
        rendered = f"{record.category.__name__}: {message}"
        if rendered not in messages:
            messages.append(rendered)
    return messages


def _ta_references(operation):
    scanpy = {
        "citation": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19:15.",
        "doi": "10.1186/s13059-017-1382-0",
        "url": "https://doi.org/10.1186/s13059-017-1382-0",
        "kind": "software",
    }
    if operation == "diffusion_map":
        return [
            {
                "citation": "Coifman RR, Lafon S, Lee AB, et al. Geometric diffusions as a tool for harmonic analysis and structure definition of data. PNAS. 2005;102:7426-7431.",
                "doi": "10.1073/pnas.0500334102",
                "url": "https://doi.org/10.1073/pnas.0500334102",
                "kind": "method",
            },
            {
                "citation": "Haghverdi L, Buettner F, Theis FJ. Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics. 2015;31:2989-2998.",
                "doi": "10.1093/bioinformatics/btv325",
                "url": "https://doi.org/10.1093/bioinformatics/btv325",
                "kind": "method",
            },
            scanpy,
        ]
    if operation == "paga":
        return [
            {
                "citation": "Wolf FA, Hamey FK, Plass M, et al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. Genome Biology. 2019;20:59.",
                "doi": "10.1186/s13059-019-1663-x",
                "url": "https://doi.org/10.1186/s13059-019-1663-x",
                "kind": "method",
            },
            scanpy,
        ]
    if operation == "dpt":
        return [
            {
                "citation": "Haghverdi L, Buttner M, Wolf FA, Buettner F, Theis FJ. Diffusion pseudotime robustly reconstructs lineage branching. Nature Methods. 2016;13:845-848.",
                "doi": "10.1038/nmeth.3971",
                "url": "https://doi.org/10.1038/nmeth.3971",
                "kind": "method",
            },
            scanpy,
        ]
    raise ValueError(f"Unknown trajectory reference family: {operation!r}.")


def _ta_summary(
    *, node_id, methods, results, key_results, parameters, warnings_list, limitations, references, versions
):
    summary = {
        "schema_version": 1,
        "node_id": node_id,
        "methods": str(methods),
        "results": str(results),
        "key_results": _ta_json_value(key_results),
        "parameters": _ta_json_value(parameters),
        "warnings": [str(value) for value in warnings_list],
        "limitations": [str(value) for value in limitations],
        "references": _ta_json_value(references),
        "software_versions": _ta_json_value(versions),
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def _ta_validate_dense_float(array, *, shape, label, numpy):
    values = numpy.asarray(array)
    if tuple(values.shape) != tuple(shape):
        raise RuntimeError(f"{label} has shape {tuple(values.shape)!r}; expected {tuple(shape)!r}.")
    if not numpy.issubdtype(values.dtype, numpy.floating) or numpy.iscomplexobj(values):
        raise RuntimeError(f"{label} must contain real floating-point values.")
    if not bool(numpy.isfinite(values).all()):
        raise RuntimeError(f"{label} contains non-finite values.")
    return values


def _ta_diffmap_output_fingerprint(coordinates, eigenvalues, *, numpy):
    payload = {
        "coordinates": _ta_hash_array(coordinates, numpy=numpy),
        "eigenvalues": _ta_hash_array(eigenvalues, numpy=numpy),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    ).hexdigest()


def _ta_transition_matrix(connectivities, *, numpy):
    density = numpy.asarray(connectivities.sum(axis=0), dtype=float).ravel()
    if not bool(numpy.isfinite(density).all()) or bool((density <= 0).any()):
        raise ValueError("Diffusion Map requires every cell to have positive connectivity degree.")
    inverse_density = 1.0 / density
    normalized = connectivities.multiply(inverse_density[:, None]).multiply(inverse_density[None, :]).tocsr()
    normalization = numpy.sqrt(numpy.asarray(normalized.sum(axis=0), dtype=float).ravel())
    if not bool(numpy.isfinite(normalization).all()) or bool((normalization <= 0).any()):
        raise ValueError("Diffusion Map transition normalization is undefined for the named graph.")
    inverse_normalization = 1.0 / normalization
    transition = normalized.multiply(inverse_normalization[:, None]).multiply(inverse_normalization[None, :])
    transition = transition.tocsr()
    transition.sum_duplicates()
    transition.eliminate_zeros()
    transition.sort_indices()
    return transition


def _ta_validate_diffmap_eigenpairs(connectivities, coordinates, eigenvalues, *, numpy):
    transition = _ta_transition_matrix(connectivities, numpy=numpy)
    vectors = numpy.asarray(coordinates, dtype=float)
    values = numpy.asarray(eigenvalues, dtype=float)
    if bool((values < -1.00001).any()) or bool((values > 1.00001).any()):
        raise ValueError("Diffusion Map eigenvalues fall outside the symmetric transition-operator spectrum.")
    gram = vectors.T @ vectors
    orthogonality_error = float(numpy.max(numpy.abs(gram - numpy.eye(vectors.shape[1]))))
    if orthogonality_error > 5e-5:
        raise ValueError("Diffusion Map coordinates are not an orthonormal transition-operator eigenbasis.")
    residuals = []
    for index, value in enumerate(values):
        vector = vectors[:, index]
        projected = numpy.asarray(transition @ vector, dtype=float).ravel()
        residual = float(numpy.linalg.norm(projected - value * vector))
        scale = max(float(numpy.linalg.norm(projected)), abs(float(value)), 1e-15)
        residuals.append(residual / scale)
    maximum_residual = max(residuals, default=0.0)
    if maximum_residual > 5e-5:
        raise ValueError("Diffusion Map coordinates/eigenvalues fail the named-graph eigenpair residual check.")
    return {
        "maximum_relative_eigenpair_residual": maximum_residual,
        "maximum_orthogonality_error": orthogonality_error,
        "validation_tolerance": 5e-5,
    }


def _ta_validate_diffmap_bundle(adata, *, graph, numpy):
    if "X_diffmap" not in adata.obsm or "diffmap_evals" not in adata.uns:
        raise ValueError("DPT requires the complete Diffusion Map bundle; run Diffusion Map first.")
    coordinates = numpy.asarray(adata.obsm["X_diffmap"])
    eigenvalues = numpy.asarray(adata.uns["diffmap_evals"])
    if coordinates.ndim != 2 or coordinates.shape[0] != len(graph["obs_names"]):
        raise ValueError("DPT Diffusion Map coordinates are not aligned to the current observations.")
    n_comps = int(coordinates.shape[1])
    if n_comps < 3:
        raise ValueError("DPT requires a Diffusion Map basis with at least three components.")
    coordinates = _ta_validate_dense_float(
        coordinates,
        shape=(len(graph["obs_names"]), n_comps),
        label="Diffusion Map coordinates",
        numpy=numpy,
    )
    eigenvalues = _ta_validate_dense_float(
        eigenvalues,
        shape=(n_comps,),
        label="Diffusion Map eigenvalues",
        numpy=numpy,
    )
    if bool(numpy.any(numpy.diff(eigenvalues.astype(float)) > 1e-6)):
        raise ValueError("DPT Diffusion Map eigenvalues are not in non-increasing order.")
    provenance = adata.uns.get(_TA_DIFFMAP_PROVENANCE_KEY)
    expected = {
        "neighbors_key": graph["neighbors_key"],
        "graph_fingerprint_sha256": graph["graph_fingerprint_sha256"],
        "observation_fingerprint_sha256": graph["observation_fingerprint_sha256"],
        "n_comps": n_comps,
        "output_fingerprint_sha256": _ta_diffmap_output_fingerprint(coordinates, eigenvalues, numpy=numpy),
    }
    if provenance is None:
        provenance = expected
    else:
        if not isinstance(provenance, Mapping) or provenance.get("schema") != _TA_DIFFMAP_SCHEMA:
            raise ValueError("DPT Diffusion Map provenance has an unsupported schema.")
        for key, value in expected.items():
            if provenance.get(key) != value:
                raise ValueError(f"DPT Diffusion Map provenance mismatch for {key!r}; recompute Diffusion Map.")
    validation = _ta_validate_diffmap_eigenpairs(graph["connectivities"], coordinates, eigenvalues, numpy=numpy)
    return coordinates, eigenvalues, dict(provenance), validation


def _ta_run_diffusion_map(
    adata,
    *,
    neighbors_key="neighbors",
    n_comps=15,
    overwrite_existing=False,
    random_seed=0,
    openbio_version="unknown",
):
    operation = "Diffusion Map"
    scanpy, numpy, _pandas, scipy_sparse, scipy_csgraph, scanpy_version = _ta_import_science(operation)
    _ta_require_signature(
        scanpy.tl.diffmap,
        operation=operation,
        expected_names=("adata", "n_comps", "neighbors_key", "random_state", "copy"),
    )
    graph = _ta_graph(
        adata,
        neighbors_key=neighbors_key,
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    n_comps = _ta_integer(n_comps, label="Diffusion Map n_comps", minimum=3, maximum=4096)
    if n_comps > len(graph["obs_names"]) - 1:
        raise ValueError("Diffusion Map n_comps must not exceed n_obs - 1 for Scanpy's sparse eigensolver.")
    random_seed = _ta_integer(random_seed, label="Diffusion Map random_seed", minimum=0, maximum=2**31 - 1)
    overwrite_existing = _ta_boolean(overwrite_existing, label="Diffusion Map overwrite_existing")
    downstream_dpt = [
        location
        for location, container, key in (
            ("obs['dpt_pseudotime']", adata.obs, "dpt_pseudotime"),
            ("uns['iroot']", adata.uns, "iroot"),
            (f"uns[{_TA_DPT_PROVENANCE_KEY!r}]", adata.uns, _TA_DPT_PROVENANCE_KEY),
        )
        if key in container
    ]
    if downstream_dpt:
        raise ValueError(
            "Diffusion Map refuses to overwrite or coexist with a downstream DPT bundle "
            f"({', '.join(downstream_dpt)}); remove/recompute DPT atomically after the new diffusion basis."
        )
    collisions = []
    for container_name, container, key in (
        ("obsm", adata.obsm, "X_diffmap"),
        ("uns", adata.uns, "diffmap_evals"),
        ("uns", adata.uns, _TA_DIFFMAP_PROVENANCE_KEY),
    ):
        if key in container:
            collisions.append(f"{container_name}[{key!r}]")
    if collisions and not overwrite_existing:
        raise ValueError(f"Diffusion Map output bundle already exists: {', '.join(collisions)}")
    output = adata
    output.obsp[graph["connectivities_key"]] = graph["connectivities"]
    output.obsp[graph["distances_key"]] = graph["distances"]
    if overwrite_existing:
        for container, key in (
            (output.obsm, "X_diffmap"),
            (output.uns, "diffmap_evals"),
            (output.uns, _TA_DIFFMAP_PROVENANCE_KEY),
        ):
            if key in container:
                del container[key]
    with _TA_RUNTIME_LOCK:
        rng_state = numpy.random.get_state()
        print_options = numpy.get_printoptions()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                backend_return = scanpy.tl.diffmap(
                    output,
                    n_comps=n_comps,
                    neighbors_key=graph["neighbors_key"],
                    random_state=random_seed,
                    copy=False,
                )
        finally:
            numpy.random.set_state(rng_state)
            numpy.set_printoptions(**print_options)
    if backend_return is not None:
        raise RuntimeError("Diffusion Map backend violated copy=False by returning a value.")
    _ta_assert_exact_graph_storage(
        output,
        graph=graph,
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    post_graph = _ta_graph(
        output,
        neighbors_key=graph["neighbors_key"],
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    if post_graph["graph_fingerprint_sha256"] != graph["graph_fingerprint_sha256"]:
        raise RuntimeError("Diffusion Map backend modified the validated named graph.")
    coordinates = _ta_validate_dense_float(
        output.obsm.get("X_diffmap"),
        shape=(len(graph["obs_names"]), n_comps),
        label="Diffusion Map coordinates",
        numpy=numpy,
    )
    eigenvalues = _ta_validate_dense_float(
        output.uns.get("diffmap_evals"),
        shape=(n_comps,),
        label="Diffusion Map eigenvalues",
        numpy=numpy,
    )
    if bool(numpy.any(numpy.diff(eigenvalues.astype(float)) > 1e-6)):
        raise RuntimeError("Diffusion Map backend returned eigenvalues outside non-increasing order.")
    try:
        scientific_validation = _ta_validate_diffmap_eigenpairs(
            graph["connectivities"], coordinates, eigenvalues, numpy=numpy
        )
    except ValueError as error:
        raise RuntimeError(f"Diffusion Map backend returned scientifically inconsistent eigenpairs: {error}") from error
    output_fingerprint = _ta_diffmap_output_fingerprint(coordinates, eigenvalues, numpy=numpy)
    versions = collect_software_versions(
        ("scanpy", "anndata", "numpy", "pandas", "scipy"),
        openbio_version=str(openbio_version),
    )
    provenance = {
        "schema": _TA_DIFFMAP_SCHEMA,
        "neighbors_key": graph["neighbors_key"],
        "connectivities_key": graph["connectivities_key"],
        "distances_key": graph["distances_key"],
        "graph_fingerprint_sha256": graph["graph_fingerprint_sha256"],
        "observation_fingerprint_sha256": graph["observation_fingerprint_sha256"],
        "n_comps": n_comps,
        "random_seed": random_seed,
        "output_fingerprint_sha256": output_fingerprint,
        "scanpy_version": scanpy_version,
    }
    output.uns[_TA_DIFFMAP_PROVENANCE_KEY] = provenance
    backend_warnings = _ta_warning_messages(caught)
    report_warnings = list(backend_warnings)
    if len(graph["component_sizes"]) > 1:
        report_warnings.append(
            "The named graph is disconnected; diffusion geometry spans disconnected components and is not a single rooted ordering."
        )
    component_summaries = [
        {"component_index": index, **_ta_numeric_summary(coordinates[:, index], numpy=numpy)}
        for index in range(1, min(n_comps, _TA_MAX_DIFFMAP_COMPONENT_SUMMARIES + 1))
    ]
    component_summaries_truncated = n_comps - 1 > _TA_MAX_DIFFMAP_COMPONENT_SUMMARIES
    if component_summaries_truncated:
        report_warnings.append(
            f"Informative component summaries were truncated to {_TA_MAX_DIFFMAP_COMPONENT_SUMMARIES} of {n_comps - 1}."
        )
    parameters = {
        "neighbors_key": graph["neighbors_key"],
        "n_comps": n_comps,
        "overwrite_existing": overwrite_existing,
        "random_seed": random_seed,
        "scanpy_fixed_arguments": {"copy": False},
    }
    key_results = {
        "graph": _ta_graph_diagnostics(graph, numpy=numpy),
        "outputs": {
            "coordinates_key": "X_diffmap",
            "coordinates_shape": [int(value) for value in coordinates.shape],
            "coordinates_dtype": str(coordinates.dtype),
            "eigenvalues_key": "diffmap_evals",
            "eigenvalues_shape": [int(value) for value in eigenvalues.shape],
            "eigenvalues_dtype": str(eigenvalues.dtype),
            "provenance_key": _TA_DIFFMAP_PROVENANCE_KEY,
            "output_fingerprint_sha256": output_fingerprint,
        },
        "eigenvalues": [float(value) for value in eigenvalues],
        "scientific_validation": scientific_validation,
        "informative_component_summaries": component_summaries,
        "informative_component_summary_count": n_comps - 1,
        "informative_component_summaries_truncated": component_summaries_truncated,
        "informative_component_summaries_limit": _TA_MAX_DIFFMAP_COMPONENT_SUMMARIES,
        "stationary_component_index": 0,
        "stationary_component_is_not_informative": True,
    }
    summary = _ta_summary(
        node_id=DIFFMAP_NODE_ID,
        methods=(
            f"Computed {n_comps} diffusion components with scanpy.tl.diffmap {scanpy_version} from the exact "
            f"named graph {graph['neighbors_key']!r}; the graph, axes, output bundle, seed, and backend signature "
            "were validated and fingerprinted."
        ),
        results=(
            f"The diffusion basis contains {n_comps - 1} informative components after the stationary column; "
            f"the named graph contained {len(graph['component_sizes'])} connected component(s)."
        ),
        key_results=key_results,
        parameters=parameters,
        warnings_list=report_warnings,
        limitations=(
            "Diffusion components are exploratory and their signs are arbitrary.",
            "The result does not select a biological root or estimate pseudotime or branches.",
            "No cluster, Sample-level, Condition, or causal inference was performed.",
        ),
        references=_ta_references("diffusion_map"),
        versions=versions,
    )
    return output, summary


def _ta_label_record(value, *, label, numpy):
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, bool):
        raise TypeError(f"{label} boolean labels are not supported.")
    if isinstance(value, str):
        if not value or value != value.strip():
            raise ValueError(f"{label} string labels must be nonempty and canonical.")
        return {"type": "string", "value": value, "display": value}
    if isinstance(value, int):
        return {"type": "integer", "value": int(value), "display": str(value)}
    if isinstance(value, float):
        if not bool(numpy.isfinite(value)):
            raise ValueError(f"{label} numeric labels must be finite.")
        return {"type": "number", "value": float(value), "display": repr(float(value))}
    raise TypeError(f"{label} labels must be canonical strings, integers, or finite numbers.")


def _ta_partition(adata, *, groupby, operation, pandas, numpy):
    groupby = _ta_clean_text(groupby, label=f"{operation} groupby")
    if groupby not in adata.obs:
        raise ValueError(f"{operation} groupby column not found in obs[{groupby!r}].")
    series = adata.obs[groupby]
    if not isinstance(series.dtype, pandas.CategoricalDtype):
        raise TypeError(f"{operation} obs[{groupby!r}] must be explicitly categorical.")
    if bool(series.isna().any()):
        raise ValueError(f"{operation} obs[{groupby!r}] contains missing labels.")
    category_records = [
        _ta_label_record(value, label=f"{operation} groupby", numpy=numpy) for value in series.cat.categories
    ]
    displays = [record["display"] for record in category_records]
    if len(displays) != len(set(displays)):
        raise ValueError(f"{operation} group categories collide after display rendering.")
    counts = series.value_counts(sort=False)
    used_mask = [int(counts.iloc[index]) > 0 for index in range(len(category_records))]
    unused = [record for record, used in zip(category_records, used_mask, strict=True) if not used]
    represented = [record for record, used in zip(category_records, used_mask, strict=True) if used]
    represented_counts = [int(value) for value, used in zip(counts.tolist(), used_mask, strict=True) if used]
    if len(represented) < 2:
        raise ValueError(f"{operation} requires at least two represented groups.")
    category_codes = numpy.asarray(series.cat.codes, dtype="<i8")
    digest = hashlib.sha256()
    category_payload = json.dumps(
        {"groupby": groupby, "categories": represented},
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    digest.update(len(category_payload).to_bytes(8, "little", signed=False))
    digest.update(category_payload)
    digest.update(category_codes.tobytes(order="C"))
    partition_fingerprint = digest.hexdigest()
    return {
        "groupby": groupby,
        "categories": represented,
        "counts": represented_counts,
        "unused_categories": unused,
        "codes": tuple(int(value) for value in category_codes),
        "partition_fingerprint_sha256": partition_fingerprint,
    }


def _ta_paga_matrix(matrix, *, n_groups, name, symmetric, numpy, scipy_sparse):
    if not scipy_sparse.issparse(matrix):
        raise RuntimeError(f"PAGA backend {name} must be a SciPy sparse matrix.")
    values = matrix.tocsr(copy=True)
    values.sum_duplicates()
    values.eliminate_zeros()
    values.sort_indices()
    if tuple(values.shape) != (n_groups, n_groups):
        raise RuntimeError(f"PAGA backend {name} has the wrong group-axis shape.")
    data = numpy.asarray(values.data)
    if data.size and (not numpy.issubdtype(data.dtype, numpy.number) or numpy.iscomplexobj(data)):
        raise RuntimeError(f"PAGA backend {name} must contain real numeric weights.")
    if data.size and (not bool(numpy.isfinite(data).all()) or bool((data < 0).any())):
        raise RuntimeError(f"PAGA backend {name} contains invalid weights.")
    if data.size and bool((data > 1.0 + 1e-8).any()):
        raise RuntimeError(f"PAGA backend {name} contains weights above the model-v1.2 bound of one.")
    if not bool(numpy.allclose(values.diagonal(), 0.0, rtol=0.0, atol=1e-12)):
        raise RuntimeError(f"PAGA backend {name} must have a zero diagonal.")
    if symmetric:
        delta = values - values.T
        delta.eliminate_zeros()
        if delta.nnz and float(numpy.max(numpy.abs(delta.data))) > 1e-8:
            raise RuntimeError(f"PAGA backend {name} must be symmetric.")
        values = ((values + values.T) * 0.5).tocsr()
        values.eliminate_zeros()
        values.sort_indices()
    return values


def _ta_expected_paga_v1_2(graph, partition, *, numpy, scipy_sparse, scipy_csgraph):
    codes = numpy.asarray(partition["codes"], dtype=int)
    n_groups = len(partition["categories"])
    directed = numpy.zeros((n_groups, n_groups), dtype=numpy.int64)
    inner = numpy.zeros(n_groups, dtype=numpy.int64)
    distances = graph["distances"].tocsr()
    for cell_index in range(distances.shape[0]):
        source = int(codes[cell_index])
        for position in range(distances.indptr[cell_index], distances.indptr[cell_index + 1]):
            target = int(codes[int(distances.indices[position])])
            if source == target:
                inner[source] += 1
            else:
                directed[source, target] += 1
    incident = inner + directed.sum(axis=1)
    pair_edges = directed + directed.T
    expected = numpy.zeros((n_groups, n_groups), dtype=float)
    denominator = int(codes.size) - 1
    for left in range(n_groups):
        for right in range(n_groups):
            observed = int(pair_edges[left, right])
            if observed == 0:
                continue
            null_edges = (
                int(incident[left]) * int(partition["counts"][right])
                + int(incident[right]) * int(partition["counts"][left])
            ) / denominator
            expected[left, right] = min(observed / null_edges, 1.0) if null_edges else 1.0
    connectivities = scipy_sparse.csr_matrix(expected)
    inverse = connectivities.copy()
    inverse.data = 1.0 / inverse.data
    minimum_tree = scipy_csgraph.minimum_spanning_tree(inverse).tocoo()
    tree = scipy_sparse.csr_matrix(
        (expected[minimum_tree.row, minimum_tree.col], (minimum_tree.row, minimum_tree.col)),
        shape=expected.shape,
    )
    return connectivities, tree


def _ta_validate_paga_bundle(output, *, graph, partition, numpy, scipy_sparse, scipy_csgraph):
    result = output.uns.get("paga")
    if not isinstance(result, Mapping) or result.get("groups") != partition["groupby"]:
        raise RuntimeError("PAGA backend did not return the requested canonical grouping key.")
    n_groups = len(partition["categories"])
    connectivities = _ta_paga_matrix(
        result.get("connectivities"),
        n_groups=n_groups,
        name="connectivities",
        symmetric=True,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    tree = _ta_paga_matrix(
        result.get("connectivities_tree"),
        n_groups=n_groups,
        name="connectivities_tree",
        symmetric=False,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    expected_connectivities, expected_tree = _ta_expected_paga_v1_2(
        graph,
        partition,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    if not bool(numpy.allclose(connectivities.toarray(), expected_connectivities.toarray(), rtol=1e-10, atol=1e-12)):
        raise RuntimeError("PAGA backend connectivities do not equal an independent model-v1.2 recomputation.")
    if not bool(numpy.allclose(tree.toarray(), expected_tree.toarray(), rtol=1e-10, atol=1e-12)):
        raise RuntimeError("PAGA backend spanning forest does not equal the model-v1.2 inverse-weight MST.")
    reciprocal = tree.multiply(tree.T)
    reciprocal.eliminate_zeros()
    if reciprocal.nnz:
        raise RuntimeError("PAGA backend tree must store one orientation per undirected edge.")
    tree_coo = tree.tocoo()
    if tree_coo.nnz:
        parent_weights = numpy.asarray(connectivities[tree_coo.row, tree_coo.col]).ravel()
        if bool((parent_weights <= 0).any()) or not bool(
            numpy.allclose(parent_weights, tree_coo.data, rtol=1e-8, atol=1e-12)
        ):
            raise RuntimeError("PAGA backend tree is not an exact weighted subgraph of connectivities.")
    full_components, _ = scipy_csgraph.connected_components(connectivities, directed=False)
    undirected_tree = (tree + tree.T).tocsr()
    tree_components, _ = scipy_csgraph.connected_components(undirected_tree, directed=False)
    if tree_components != full_components or tree.nnz != n_groups - full_components:
        raise RuntimeError("PAGA backend tree is not a spanning forest of the abstracted graph.")
    sizes_key = f"{partition['groupby']}_sizes"
    sizes = numpy.asarray(output.uns.get(sizes_key))
    if sizes.shape != (n_groups,) or not numpy.issubdtype(sizes.dtype, numpy.integer):
        raise RuntimeError("PAGA backend group-size sidecar is missing or invalid.")
    if [int(value) for value in sizes] != partition["counts"]:
        raise RuntimeError("PAGA backend group-size sidecar does not match the categorical partition.")
    post_partition = _ta_partition(
        output,
        groupby=partition["groupby"],
        operation="PAGA",
        pandas=__import__("pandas"),
        numpy=numpy,
    )
    if post_partition["partition_fingerprint_sha256"] != partition["partition_fingerprint_sha256"]:
        raise RuntimeError("PAGA backend changed the validated group partition.")
    return connectivities, tree, sizes, full_components


def _ta_paga_output_fingerprint(connectivities, tree, sizes, partition, *, numpy):
    payload = {
        "connectivities": _ta_hash_sparse(connectivities, numpy=numpy, preserve_zeros=False),
        "tree": _ta_hash_sparse(tree, numpy=numpy, preserve_zeros=False),
        "sizes": _ta_hash_array(sizes, numpy=numpy),
        "partition": partition["partition_fingerprint_sha256"],
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    ).hexdigest()


def _ta_run_paga(
    adata,
    *,
    groupby="leiden",
    neighbors_key="neighbors",
    overwrite_existing=False,
    openbio_version="unknown",
):
    operation = "PAGA"
    scanpy, numpy, pandas, scipy_sparse, scipy_csgraph, scanpy_version = _ta_import_science(operation)
    _ta_require_signature(
        scanpy.tl.paga,
        operation=operation,
        expected_names=("adata", "groups", "use_rna_velocity", "model", "neighbors_key", "copy"),
    )
    graph = _ta_graph(
        adata,
        neighbors_key=neighbors_key,
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    groupby = _ta_clean_text(groupby, label=f"{operation} groupby")
    if groupby in adata.obs and not isinstance(adata.obs[groupby].dtype, pandas.CategoricalDtype):
        grouping = adata.obs[[groupby]].copy()
        adata.strings_to_categoricals(grouping)
        adata.obs[groupby] = grouping[groupby]
    partition = _ta_partition(adata, groupby=groupby, operation=operation, pandas=pandas, numpy=numpy)
    overwrite_existing = _ta_boolean(overwrite_existing, label="PAGA overwrite_existing")
    sizes_key = f"{partition['groupby']}_sizes"
    collisions = [f"uns[{key!r}]" for key in ("paga", sizes_key, _TA_PAGA_PROVENANCE_KEY) if key in adata.uns]
    if collisions and not overwrite_existing:
        raise ValueError(f"PAGA output bundle already exists: {', '.join(collisions)}")
    output = adata
    previous_paga = adata.uns.get("paga")
    previous_sizes_key = None
    if "paga" in adata.uns:
        previous_groups = previous_paga.get("groups") if isinstance(previous_paga, Mapping) else None
        if not isinstance(previous_groups, str) or not previous_groups or previous_groups != previous_groups.strip():
            raise ValueError(
                "PAGA overwrite cannot identify the previous canonical group-size sidecar; remove the malformed "
                "uns['paga'] bundle explicitly before recomputing."
            )
        previous_sizes_key = f"{previous_groups}_sizes"
    output.obsp[graph["connectivities_key"]] = graph["connectivities"]
    output.obsp[graph["distances_key"]] = graph["distances"]
    if partition["unused_categories"]:
        output.obs[partition["groupby"]] = output.obs[partition["groupby"]].cat.remove_unused_categories()
        partition = _ta_partition(
            output,
            groupby=partition["groupby"],
            operation=operation,
            pandas=pandas,
            numpy=numpy,
        ) | {"unused_categories": partition["unused_categories"]}
    if overwrite_existing:
        for key in ("paga", sizes_key, previous_sizes_key, _TA_PAGA_PROVENANCE_KEY):
            if key is None:
                continue
            if key in output.uns:
                del output.uns[key]
    # PAGA needs a featureless AnnData workspace so backend-owned mutations do
    # not leak into the full scientific state beyond the validated PAGA bundle.
    scratch = output[:, []].copy()
    scratch.obs = output.obs[[partition["groupby"]]].copy()
    scratch.uns.clear()
    scratch.obsm.clear()
    scratch.varm.clear()
    scratch.varp.clear()
    scratch.layers.clear(keep_x=False)
    scratch.obsp.clear()
    scratch.uns[graph["neighbors_key"]] = {
        "connectivities_key": graph["connectivities_key"],
        "distances_key": graph["distances_key"],
        "params": graph["neighbor_parameters"],
    }
    scratch.obsp[graph["connectivities_key"]] = graph["connectivities"].copy()
    scratch.obsp[graph["distances_key"]] = graph["distances"].copy()
    with _TA_RUNTIME_LOCK:
        rng_state = numpy.random.get_state()
        print_options = numpy.get_printoptions()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                backend_return = scanpy.tl.paga(
                    scratch,
                    groups=partition["groupby"],
                    use_rna_velocity=False,
                    model="v1.2",
                    neighbors_key=graph["neighbors_key"],
                    copy=False,
                )
        finally:
            numpy.random.set_state(rng_state)
            numpy.set_printoptions(**print_options)
    if backend_return is not None:
        raise RuntimeError("PAGA backend violated copy=False by returning a value.")
    _ta_assert_exact_graph_storage(
        scratch,
        graph=graph,
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    scratch_paga = scratch.uns.get("paga")
    if isinstance(scratch_paga, Mapping):
        output.uns["paga"] = {
            key: value.copy() if hasattr(value, "copy") else value for key, value in scratch_paga.items()
        }
    if sizes_key in scratch.uns:
        output.uns[sizes_key] = numpy.asarray(scratch.uns[sizes_key]).copy()
    post_graph = _ta_graph(
        output,
        neighbors_key=graph["neighbors_key"],
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    if post_graph["graph_fingerprint_sha256"] != graph["graph_fingerprint_sha256"]:
        raise RuntimeError("PAGA backend modified the validated named graph.")
    connectivities, tree, sizes, abstract_components = _ta_validate_paga_bundle(
        output,
        graph=graph,
        partition=partition,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    output_fingerprint = _ta_paga_output_fingerprint(connectivities, tree, sizes, partition, numpy=numpy)
    versions = collect_software_versions(
        ("scanpy", "anndata", "numpy", "pandas", "scipy", "igraph"),
        openbio_version=str(openbio_version),
    )
    provenance = {
        "schema": _TA_PAGA_SCHEMA,
        "neighbors_key": graph["neighbors_key"],
        "graph_fingerprint_sha256": graph["graph_fingerprint_sha256"],
        "observation_fingerprint_sha256": graph["observation_fingerprint_sha256"],
        "groupby": partition["groupby"],
        "categories": partition["categories"],
        "group_counts": partition["counts"],
        "partition_fingerprint_sha256": partition["partition_fingerprint_sha256"],
        "model": "v1.2",
        "use_rna_velocity": False,
        "output_fingerprint_sha256": output_fingerprint,
        "scanpy_version": scanpy_version,
    }
    output.uns[_TA_PAGA_PROVENANCE_KEY] = provenance
    edges = []
    matrix = connectivities.tocoo()
    for left, right, weight in zip(matrix.row, matrix.col, matrix.data, strict=True):
        if int(left) >= int(right):
            continue
        edges.append(
            {
                "group_a": partition["categories"][int(left)],
                "group_b": partition["categories"][int(right)],
                "connectivity": float(weight),
                "in_spanning_forest": bool(tree[int(left), int(right)] or tree[int(right), int(left)]),
            }
        )
    edges.sort(
        key=lambda record: (
            -record["connectivity"],
            record["group_a"]["display"],
            record["group_b"]["display"],
        )
    )
    report_warnings = _ta_warning_messages(caught)
    if len(partition["categories"]) > _TA_PAGA_GROUP_WARNING_THRESHOLD:
        report_warnings.append(
            f"PAGA retains {len(partition['categories'])} represented groups; a large group-level graph "
            "can be costly to compute and difficult to interpret."
        )
    if partition["unused_categories"]:
        report_warnings.append(
            f"Removed {len(partition['unused_categories'])} globally unused categorical level(s) on the output copy."
        )
    if abstract_components > 1:
        report_warnings.append(f"The PAGA graph contains {abstract_components} disconnected group-level components.")
    edge_table_truncated = len(edges) > _TA_MAX_PAGA_SUMMARY_EDGES
    if edge_table_truncated:
        report_warnings.append(
            f"The PAGA edge table was truncated to the strongest {_TA_MAX_PAGA_SUMMARY_EDGES} of {len(edges)} edges."
        )
    parameters = {
        "groupby": partition["groupby"],
        "neighbors_key": graph["neighbors_key"],
        "overwrite_existing": overwrite_existing,
        "model": "v1.2",
        "use_rna_velocity": False,
        "copy": False,
        "group_count_warning_threshold": _TA_PAGA_GROUP_WARNING_THRESHOLD,
        "max_summary_edges": _TA_MAX_PAGA_SUMMARY_EDGES,
    }
    key_results = {
        "graph": _ta_graph_diagnostics(graph, numpy=numpy),
        "partition": {
            "groupby": partition["groupby"],
            "represented_groups": len(partition["categories"]),
            "categories": partition["categories"],
            "group_counts": partition["counts"],
            "unused_categories_removed": partition["unused_categories"],
            "partition_fingerprint_sha256": partition["partition_fingerprint_sha256"],
        },
        "abstract_graph": {
            "connectivities_shape": [int(value) for value in connectivities.shape],
            "connectivities_nnz": int(connectivities.nnz),
            "spanning_forest_shape": [int(value) for value in tree.shape],
            "spanning_forest_stored_edges": int(tree.nnz),
            "connected_components": int(abstract_components),
            "nonzero_undirected_edges": len(edges),
            "all_edges": edges[:_TA_MAX_PAGA_SUMMARY_EDGES],
            "strongest_edges": edges[:20],
            "edge_table_truncated": edge_table_truncated,
            "edge_table_limit": _TA_MAX_PAGA_SUMMARY_EDGES,
            "group_sizes_key": sizes_key,
            "provenance_key": _TA_PAGA_PROVENANCE_KEY,
            "output_fingerprint_sha256": output_fingerprint,
        },
    }
    summary = _ta_summary(
        node_id=PAGA_NODE_ID,
        methods=(
            f"Computed undirected partition-based graph abstraction with scanpy.tl.paga {scanpy_version}, "
            f"model='v1.2', over {len(partition['categories'])} represented categories in obs[{partition['groupby']!r}] "
            f"and named graph {graph['neighbors_key']!r}."
        ),
        results=(
            f"The PAGA abstraction contained {len(edges)} nonzero undirected group-pair edge(s) across "
            f"{abstract_components} connected component(s); the strongest reported edges and spanning forest are disclosed."
        ),
        key_results=key_results,
        parameters=parameters,
        warnings_list=report_warnings,
        limitations=(
            "PAGA connectivity weights are exploratory observed-to-expected connectivity evidence, not p-values or probabilities.",
            "The abstraction is conditional on the chosen cell partition and named neighbor graph.",
            "No direction, root, pseudotime, lineage proof, Sample-level, or Condition inference was performed.",
        ),
        references=_ta_references("paga"),
        versions=versions,
    )
    return output, summary


def _ta_select_root(adata, *, root_mode, root_cell_id, root_column, root_value, coordinates, n_dcs, numpy):
    root_mode = _ta_clean_text(root_mode, label="DPT root_mode")
    obs_names = tuple(adata.obs_names)
    if root_mode == "cell_id":
        root_cell_id = _ta_clean_text(root_cell_id, label="DPT root_cell_id")
        try:
            index = obs_names.index(root_cell_id)
        except ValueError as error:
            raise ValueError(f"DPT root cell {root_cell_id!r} was not found exactly in obs_names.") from error
        return index, {
            "mode": "cell_id",
            "root_cell_id": root_cell_id,
            "rationale": "The user nominated one exact biological root cell identifier.",
            "population": None,
        }
    if root_mode != "group_medoid":
        raise ValueError("DPT root_mode must be 'cell_id' or 'group_medoid'.")
    root_column = _ta_clean_text(root_column, label="DPT root_column")
    if root_column not in adata.obs:
        raise ValueError(f"DPT root column not found in obs[{root_column!r}].")
    series = adata.obs[root_column]
    if bool(series.isna().any()):
        raise ValueError(f"DPT obs[{root_column!r}] contains missing population labels.")
    target = _ta_label_record(root_value, label="DPT root_value", numpy=numpy)
    observed = [
        _ta_label_record(value, label="DPT root population", numpy=numpy) for value in series.astype(object).tolist()
    ]
    indices = [index for index, value in enumerate(observed) if value == target]
    if not indices:
        raise ValueError(
            f"DPT root population {target['display']!r} ({target['type']}) was not found exactly in obs[{root_column!r}]."
        )
    stable_indices = sorted(indices, key=lambda candidate: obs_names[candidate])
    informative = numpy.asarray(coordinates[stable_indices, 1:n_dcs], dtype=float)
    centroid = numpy.asarray(
        [
            math.fsum(float(value) for value in informative[:, column]) / len(stable_indices)
            for column in range(informative.shape[1])
        ],
        dtype=float,
    )
    distances = numpy.asarray(
        [
            math.fsum((float(value) - float(center)) ** 2 for value, center in zip(row, centroid, strict=True))
            for row in informative
        ],
        dtype=float,
    )
    minimum = float(distances.min())
    tie_tolerance = max(1e-15, abs(minimum) * 1e-12)
    tied = [
        stable_indices[position] for position, value in enumerate(distances) if float(value) <= minimum + tie_tolerance
    ]
    index = min(tied, key=lambda candidate: obs_names[candidate])
    return index, {
        "mode": "group_medoid",
        "root_cell_id": obs_names[index],
        "rationale": (
            f"The user nominated one root population; the selected cell is nearest its centroid over informative "
            f"Diffusion Map columns 1:{n_dcs}, with tolerance-aware lexical cell-ID tie breaking."
        ),
        "population": {
            "column": root_column,
            "value": target,
            "cell_count": len(indices),
            "squared_distance_to_centroid": minimum,
            "tolerance_tie_count": len(tied),
            "tie_tolerance": tie_tolerance,
        },
    }


def _ta_expected_dpt_pseudotime(coordinates, eigenvalues, *, root_index, n_dcs, reachable, numpy):
    squared = numpy.zeros(coordinates.shape[0], dtype=float)
    for component in range(n_dcs):
        value = float(eigenvalues[component])
        delta = float(coordinates[root_index, component]) - numpy.asarray(coordinates[:, component], dtype=float)
        if value < 0.9994:
            delta = (value / (1.0 - value)) * delta
        squared += delta**2
    distances = numpy.sqrt(squared)
    distances[~reachable] = numpy.inf
    maximum = float(distances[reachable].max())
    if not numpy.isfinite(maximum) or maximum <= 0:
        raise ValueError("DPT diffusion distances are degenerate and cannot define a relative ordering.")
    return distances / maximum


def _ta_dpt_ordering(pseudotime, obs_names):
    indices = sorted(range(len(obs_names)), key=lambda index: (float(pseudotime[index]), obs_names[index]))
    ordered_ids = [obs_names[index] for index in indices]
    return {
        "ordering_fingerprint_sha256": _ta_hash_text_records(ordered_ids),
        "first_cell_ids": ordered_ids[:20],
        "last_cell_ids": ordered_ids[-20:],
        "tie_breaker": "lexical cell ID after ascending pseudotime",
    }


def _ta_run_dpt(
    adata,
    *,
    neighbors_key="neighbors",
    root_mode="cell_id",
    root_cell_id="",
    root_column="",
    root_value="",
    n_dcs=10,
    overwrite_existing=False,
    openbio_version="unknown",
):
    operation = "Diffusion Pseudotime"
    scanpy, numpy, _pandas, scipy_sparse, scipy_csgraph, scanpy_version = _ta_import_science(operation)
    _ta_require_signature(
        scanpy.tl.dpt,
        operation=operation,
        expected_names=(
            "adata",
            "n_dcs",
            "n_branchings",
            "min_group_size",
            "allow_kendall_tau_shift",
            "neighbors_key",
            "copy",
        ),
    )
    graph = _ta_graph(
        adata,
        neighbors_key=neighbors_key,
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    coordinates, eigenvalues, diffusion_provenance, diffusion_validation = _ta_validate_diffmap_bundle(
        adata, graph=graph, numpy=numpy
    )
    n_dcs = _ta_integer(n_dcs, label="DPT n_dcs", minimum=1, maximum=4096)
    if n_dcs > coordinates.shape[1]:
        raise ValueError(
            f"DPT n_dcs={n_dcs} exceeds the graph-bound Diffusion Map basis ({coordinates.shape[1]} components)."
        )
    overwrite_existing = _ta_boolean(overwrite_existing, label="DPT overwrite_existing")
    root_index, root = _ta_select_root(
        adata,
        root_mode=root_mode,
        root_cell_id=root_cell_id,
        root_column=root_column,
        root_value=root_value,
        coordinates=coordinates,
        n_dcs=n_dcs,
        numpy=numpy,
    )
    reachable = graph["component_labels"] == graph["component_labels"][root_index]
    reachable_count = int(reachable.sum())
    unreachable_count = len(graph["obs_names"]) - reachable_count
    collisions = []
    for container_name, container, key in (
        ("obs", adata.obs, "dpt_pseudotime"),
        ("uns", adata.uns, "iroot"),
        ("uns", adata.uns, _TA_DPT_PROVENANCE_KEY),
    ):
        if key in container:
            collisions.append(f"{container_name}[{key!r}]")
    if collisions and not overwrite_existing:
        raise ValueError(f"DPT output bundle already exists: {', '.join(collisions)}")
    output = adata
    output.obsp[graph["connectivities_key"]] = graph["connectivities"]
    output.obsp[graph["distances_key"]] = graph["distances"]
    if overwrite_existing:
        for container, key in (
            (output.obs, "dpt_pseudotime"),
            (output.uns, "iroot"),
            (output.uns, _TA_DPT_PROVENANCE_KEY),
        ):
            if key in container:
                del container[key]
    graph_fingerprint = graph["graph_fingerprint_sha256"]
    diffusion_fingerprint = diffusion_provenance["output_fingerprint_sha256"]
    output.uns["iroot"] = int(root_index)
    with _TA_RUNTIME_LOCK:
        rng_state = numpy.random.get_state()
        print_options = numpy.get_printoptions()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                backend_return = scanpy.tl.dpt(
                    output,
                    n_dcs=n_dcs,
                    n_branchings=0,
                    min_group_size=0.01,
                    allow_kendall_tau_shift=True,
                    neighbors_key=graph["neighbors_key"],
                    copy=False,
                )
        finally:
            numpy.random.set_state(rng_state)
            numpy.set_printoptions(**print_options)
    if backend_return is not None:
        raise RuntimeError("DPT backend violated copy=False by returning a value.")
    _ta_assert_exact_graph_storage(
        output,
        graph=graph,
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    post_graph = _ta_graph(
        output,
        neighbors_key=graph["neighbors_key"],
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        scipy_csgraph=scipy_csgraph,
    )
    if post_graph["graph_fingerprint_sha256"] != graph_fingerprint:
        raise RuntimeError("DPT backend modified the validated named graph.")
    post_coordinates, post_eigenvalues, post_diffusion_provenance, _post_validation = _ta_validate_diffmap_bundle(
        output, graph=post_graph, numpy=numpy
    )
    if (
        post_diffusion_provenance["output_fingerprint_sha256"] != diffusion_fingerprint
        or not bool(numpy.array_equal(post_coordinates, coordinates))
        or not bool(numpy.array_equal(post_eigenvalues, eigenvalues))
    ):
        raise RuntimeError("DPT backend modified the graph-bound Diffusion Map bundle.")
    pseudotime = numpy.asarray(output.obs.get("dpt_pseudotime"))
    if pseudotime.shape != (len(graph["obs_names"]),):
        raise RuntimeError("DPT pseudotime must align to the current observations.")
    finite_pseudotime = _ta_validate_dense_float(
        pseudotime[reachable],
        shape=(reachable_count,),
        label="DPT pseudotime",
        numpy=numpy,
    )
    if not bool(numpy.isposinf(pseudotime[~reachable]).all()):
        raise RuntimeError("DPT must mark cells outside the root component with positive infinity.")
    if bool((finite_pseudotime < -1e-7).any()) or bool((finite_pseudotime > 1.0 + 1e-7).any()):
        raise RuntimeError("DPT backend returned pseudotime values outside [0, 1].")
    if not bool(numpy.isclose(pseudotime[root_index], 0.0, rtol=0.0, atol=1e-7)):
        raise RuntimeError("DPT backend did not assign zero pseudotime to the selected root cell.")
    expected_pseudotime = _ta_expected_dpt_pseudotime(
        coordinates,
        eigenvalues,
        root_index=root_index,
        n_dcs=n_dcs,
        reachable=reachable,
        numpy=numpy,
    )
    if not bool(numpy.allclose(pseudotime, expected_pseudotime, rtol=2e-6, atol=2e-7)):
        raise RuntimeError("DPT backend pseudotime does not equal the independently recomputed diffusion distance.")
    if output.uns.get("iroot") != root_index:
        raise RuntimeError("DPT backend changed the selected root index.")
    output_fingerprint = _ta_hash_array(pseudotime, numpy=numpy)
    versions = collect_software_versions(
        ("scanpy", "anndata", "numpy", "pandas", "scipy"),
        openbio_version=str(openbio_version),
    )
    provenance = {
        "schema": _TA_DPT_SCHEMA,
        "neighbors_key": graph["neighbors_key"],
        "graph_fingerprint_sha256": graph_fingerprint,
        "observation_fingerprint_sha256": graph["observation_fingerprint_sha256"],
        "diffusion_output_fingerprint_sha256": diffusion_fingerprint,
        "n_dcs": n_dcs,
        "n_branchings": 0,
        "root": root,
        "root_index": root_index,
        "output_fingerprint_sha256": output_fingerprint,
        "scanpy_version": scanpy_version,
    }
    output.uns[_TA_DPT_PROVENANCE_KEY] = provenance
    report_warnings = _ta_warning_messages(caught)
    if n_dcs == 1:
        report_warnings.append(
            "DPT n_dcs=1 uses only the stationary diffusion component; its executable ordering does not "
            "capture informative diffusion geometry."
        )
    if unreachable_count:
        report_warnings.append(
            f"DPT retained {unreachable_count} cells unreachable from the selected root with Scanpy's positive-infinity "
            "pseudotime sentinel; numeric summaries and ordering cover only the root component."
        )
    source_provenance_available = diffusion_provenance.get("schema") == _TA_DIFFMAP_SCHEMA
    if not source_provenance_available:
        report_warnings.append(
            "OpenBio Diffusion Map provenance was unavailable; the supplied coordinates and eigenvalues "
            "were checked against the selected graph, but their original computation settings are unverified."
        )
    parameters = {
        "neighbors_key": graph["neighbors_key"],
        "root_mode": root["mode"],
        "root_cell_id": root["root_cell_id"],
        "root_column": root["population"]["column"] if root["population"] else None,
        "root_value": root["population"]["value"] if root["population"] else None,
        "n_dcs": n_dcs,
        "overwrite_existing": overwrite_existing,
        "scanpy_fixed_arguments": {
            "n_branchings": 0,
            "min_group_size": 0.01,
            "allow_kendall_tau_shift": True,
            "copy": False,
        },
    }
    key_results = {
        "root": root,
        "root_index": root_index,
        "graph": _ta_graph_diagnostics(graph, numpy=numpy),
        "diffusion": {
            "n_available_components": int(coordinates.shape[1]),
            "source_provenance_available": source_provenance_available,
            "n_dcs_used": n_dcs,
            "graph_fingerprint_sha256": graph_fingerprint,
            "output_fingerprint_sha256": diffusion_fingerprint,
            "scientific_validation": diffusion_validation,
        },
        "pseudotime": {
            **_ta_numeric_summary(finite_pseudotime, numpy=numpy),
            "unreachable_cells": unreachable_count,
            "unreachable_value": "positive_infinity",
            "root_value": float(pseudotime[root_index]),
            "cells_at_or_before_0_05": int(numpy.sum(finite_pseudotime <= 0.05)),
            "cells_at_or_after_0_95": int(numpy.sum(finite_pseudotime >= 0.95)),
            "output_key": "dpt_pseudotime",
            "provenance_key": _TA_DPT_PROVENANCE_KEY,
            "output_fingerprint_sha256": output_fingerprint,
            "scientific_validation": {
                "formula": "Scanpy 1.12.3 DPT diffusion distance, normalized by the finite maximum",
                "comparison_rtol": 2e-6,
                "comparison_atol": 2e-7,
            },
            "stable_ordering": _ta_dpt_ordering(
                finite_pseudotime,
                tuple(name for name, included in zip(graph["obs_names"], reachable, strict=True) if included),
            ),
        },
    }
    summary = _ta_summary(
        node_id=DPT_NODE_ID,
        methods=(
            f"Computed single-root Diffusion Pseudotime with scanpy.tl.dpt {scanpy_version}, n_dcs={n_dcs}, "
            f"n_branchings=0, from exact root cell {root['root_cell_id']!r}. The named graph and upstream "
            "Diffusion Map eigenpairs were verified before and after execution."
        ),
        results=(
            f"Relative pseudotime ranged from {float(finite_pseudotime.min()):.6g} to {float(finite_pseudotime.max()):.6g} "
            f"across {reachable_count} reachable cells, with root {root['root_cell_id']!r} fixed at zero."
        ),
        key_results=key_results,
        parameters=parameters,
        warnings_list=report_warnings,
        limitations=(
            "DPT is a root-dependent relative ordering, not measured time or causal transition direction.",
            "The result is conditional on preprocessing, the named graph, diffusion basis, and user-nominated biological root.",
            "No branch detection, lineage-tracing proof, Sample-level, or Condition inference was performed.",
        ),
        references=_ta_references("dpt"),
        versions=versions,
    )
    return output, summary


_TA_COMMON_STANDALONE_HELPERS = (
    _ta_import_science,
    _ta_require_signature,
    _ta_clean_text,
    _ta_integer,
    _ta_boolean,
    _ta_json_value,
    _ta_obs_names,
    _ta_hash_text_records,
    _ta_hash_sparse,
    _ta_hash_array,
    _ta_sparse_graph_matrix,
    _ta_graph,
    _ta_assert_exact_graph_storage,
    _ta_numeric_summary,
    _ta_graph_diagnostics,
    _ta_warning_messages,
    _package_version,
    collect_software_versions,
    _ta_references,
    _ta_summary,
    _ta_validate_dense_float,
)

_TA_DIFFMAP_STANDALONE_HELPERS = (
    *_TA_COMMON_STANDALONE_HELPERS,
    _ta_diffmap_output_fingerprint,
    _ta_transition_matrix,
    _ta_validate_diffmap_eigenpairs,
    _ta_run_diffusion_map,
)

_TA_PAGA_STANDALONE_HELPERS = (
    *_TA_COMMON_STANDALONE_HELPERS,
    _ta_label_record,
    _ta_partition,
    _ta_paga_matrix,
    _ta_expected_paga_v1_2,
    _ta_validate_paga_bundle,
    _ta_paga_output_fingerprint,
    _ta_run_paga,
)

_TA_DPT_STANDALONE_HELPERS = (
    *_TA_COMMON_STANDALONE_HELPERS,
    _ta_diffmap_output_fingerprint,
    _ta_transition_matrix,
    _ta_validate_diffmap_eigenpairs,
    _ta_validate_diffmap_bundle,
    _ta_label_record,
    _ta_select_root,
    _ta_expected_dpt_pseudotime,
    _ta_dpt_ordering,
    _ta_run_dpt,
)


def analyze_diffusion_map(
    adata: Any,
    *,
    neighbors_key: str = "neighbors",
    n_comps: int = 15,
    overwrite_existing: bool = False,
    random_seed: int = 0,
) -> tuple[Any, dict[str, Any]]:
    return _ta_run_diffusion_map(
        adata,
        neighbors_key=neighbors_key,
        n_comps=n_comps,
        overwrite_existing=overwrite_existing,
        random_seed=random_seed,
        openbio_version=PLUGIN_VERSION,
    )


def analyze_paga(
    adata: Any,
    *,
    groupby: str = "leiden",
    neighbors_key: str = "neighbors",
    overwrite_existing: bool = False,
) -> tuple[Any, dict[str, Any]]:
    return _ta_run_paga(
        adata,
        groupby=groupby,
        neighbors_key=neighbors_key,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
    )


def analyze_dpt(
    adata: Any,
    *,
    neighbors_key: str = "neighbors",
    root_mode: str = "cell_id",
    root_cell_id: str = "",
    root_column: str = "",
    root_value: Any = "",
    n_dcs: int = 10,
    overwrite_existing: bool = False,
) -> tuple[Any, dict[str, Any]]:
    return _ta_run_dpt(
        adata,
        neighbors_key=neighbors_key,
        root_mode=root_mode,
        root_cell_id=root_cell_id,
        root_column=root_column,
        root_value=root_value,
        n_dcs=n_dcs,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
    )


def _trajectory_source(*, helpers: tuple[object, ...], declarations: str) -> str:
    helper_source = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in helpers)
    return f"""from __future__ import annotations

import hashlib
import inspect
import json
import math
import threading
import warnings
from collections.abc import Mapping, Sequence

{declarations}

{helper_source}
"""


def diffusion_map_code(*, neighbors_key: str, n_comps: int, overwrite_existing: bool, random_seed: int) -> str:
    return (
        _trajectory_source(
            helpers=_TA_DIFFMAP_STANDALONE_HELPERS,
            declarations=f"""_TA_DIFFMAP_SCHEMA = {_TA_DIFFMAP_SCHEMA!r}
_TA_DIFFMAP_PROVENANCE_KEY = {_TA_DIFFMAP_PROVENANCE_KEY!r}
_TA_DPT_PROVENANCE_KEY = {_TA_DPT_PROVENANCE_KEY!r}
_TA_MAX_GRAPH_COMPONENT_SIZES = {_TA_MAX_GRAPH_COMPONENT_SIZES!r}
_TA_MAX_DIFFMAP_COMPONENT_SUMMARIES = {_TA_MAX_DIFFMAP_COMPONENT_SUMMARIES!r}
_TA_RUNTIME_LOCK = threading.RLock()
DIFFMAP_NODE_ID = {DIFFMAP_NODE_ID!r}""",
        )
        + f"""

def run_diffusion_map(adata):
    \"\"\"Return (annotated_adata, strict_summary) without importing OpenBio helpers.\"\"\"
    return _ta_run_diffusion_map(
        adata,
        neighbors_key={neighbors_key!r},
        n_comps={n_comps!r},
        overwrite_existing={overwrite_existing!r},
        random_seed={random_seed!r},
        openbio_version={PLUGIN_VERSION!r},
    )
"""
    )


def paga_code(*, groupby: str, neighbors_key: str, overwrite_existing: bool) -> str:
    return (
        _trajectory_source(
            helpers=_TA_PAGA_STANDALONE_HELPERS,
            declarations=f"""_TA_PAGA_SCHEMA = {_TA_PAGA_SCHEMA!r}
_TA_PAGA_PROVENANCE_KEY = {_TA_PAGA_PROVENANCE_KEY!r}
_TA_PAGA_GROUP_WARNING_THRESHOLD = {_TA_PAGA_GROUP_WARNING_THRESHOLD!r}
_TA_MAX_PAGA_SUMMARY_EDGES = {_TA_MAX_PAGA_SUMMARY_EDGES!r}
_TA_MAX_GRAPH_COMPONENT_SIZES = {_TA_MAX_GRAPH_COMPONENT_SIZES!r}
_TA_RUNTIME_LOCK = threading.RLock()
PAGA_NODE_ID = {PAGA_NODE_ID!r}""",
        )
        + f"""

def run_paga(adata):
    \"\"\"Return (annotated_adata, strict_summary) without importing OpenBio helpers.\"\"\"
    return _ta_run_paga(
        adata,
        groupby={groupby!r},
        neighbors_key={neighbors_key!r},
        overwrite_existing={overwrite_existing!r},
        openbio_version={PLUGIN_VERSION!r},
    )
"""
    )


def dpt_code(
    *,
    neighbors_key: str,
    root_mode: str,
    root_cell_id: str,
    root_column: str,
    root_value: Any,
    n_dcs: int,
    overwrite_existing: bool,
) -> str:
    return (
        _trajectory_source(
            helpers=_TA_DPT_STANDALONE_HELPERS,
            declarations=f"""_TA_DIFFMAP_SCHEMA = {_TA_DIFFMAP_SCHEMA!r}
_TA_DPT_SCHEMA = {_TA_DPT_SCHEMA!r}
_TA_DIFFMAP_PROVENANCE_KEY = {_TA_DIFFMAP_PROVENANCE_KEY!r}
_TA_DPT_PROVENANCE_KEY = {_TA_DPT_PROVENANCE_KEY!r}
_TA_MAX_GRAPH_COMPONENT_SIZES = {_TA_MAX_GRAPH_COMPONENT_SIZES!r}
_TA_RUNTIME_LOCK = threading.RLock()
DPT_NODE_ID = {DPT_NODE_ID!r}""",
        )
        + f"""

def run_dpt(adata):
    \"\"\"Return (annotated_adata, strict_summary) without importing OpenBio helpers.\"\"\"
    return _ta_run_dpt(
        adata,
        neighbors_key={neighbors_key!r},
        root_mode={root_mode!r},
        root_cell_id={root_cell_id!r},
        root_column={root_column!r},
        root_value={root_value!r},
        n_dcs={n_dcs!r},
        overwrite_existing={overwrite_existing!r},
        openbio_version={PLUGIN_VERSION!r},
    )
"""
    )


__all__ = [
    "DIFFMAP_NODE_ID",
    "DPT_NODE_ID",
    "PAGA_NODE_ID",
    "analyze_diffusion_map",
    "analyze_dpt",
    "analyze_paga",
    "diffusion_map_code",
    "dpt_code",
    "paga_code",
]
