from __future__ import annotations

import inspect
import textwrap
from typing import Any


def _standalone_plain_json(value):
    from collections.abc import Mapping, Sequence

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if value == value and value not in (float("inf"), float("-inf")) else None
    if hasattr(value, "item"):
        try:
            return _standalone_plain_json(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _standalone_plain_json(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_standalone_plain_json(item) for item in value]
    if hasattr(value, "tolist"):
        return _standalone_plain_json(value.tolist())
    return str(value)


def _standalone_validate_settings(
    adata,
    *,
    random_seed,
    neighbors_key,
    key_added,
    posterior_samples,
    degree_correction,
    overwrite_existing,
    max_working_memory_gib,
):
    import math
    import re

    if bool(adata.isbacked):
        raise ValueError("Schist Nested-SBM Hierarchy requires an in-memory AnnData; call adata.to_memory() first.")
    if int(adata.n_obs) < 2:
        raise ValueError("Schist Nested-SBM Hierarchy requires at least two observations.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError("Schist Nested-SBM Hierarchy requires unique observation identifiers.")
    if isinstance(random_seed, bool) or not isinstance(random_seed, int):
        raise TypeError("Schist random_seed must be an integer.")
    if not 1 <= random_seed <= 2**31 - 1:
        raise ValueError("Schist random_seed must be between 1 and 2**31 - 1; zero leaves graph-tool unseeded.")
    if not isinstance(neighbors_key, str):
        raise TypeError("Schist neighbors_key must be a string.")
    neighbors_key = neighbors_key.strip()
    if not neighbors_key:
        raise ValueError("Schist neighbors_key cannot be empty.")
    if not isinstance(key_added, str):
        raise TypeError("Schist key_added must be a string.")
    key_added = key_added.strip()
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,63}", key_added):
        raise ValueError(
            "Schist key_added must start with a letter and contain at most 64 ASCII letters, digits, or underscores."
        )
    if isinstance(posterior_samples, bool) or not isinstance(posterior_samples, int):
        raise TypeError("Schist posterior_samples must be an integer.")
    if not 100 <= posterior_samples <= 2**31 - 1:
        raise ValueError("Schist posterior_samples must be at least 100 and no greater than 2**31 - 1.")
    if not isinstance(degree_correction, bool):
        raise TypeError("Schist degree_correction must be a boolean.")
    if not isinstance(overwrite_existing, bool):
        raise TypeError("Schist overwrite_existing must be a boolean.")
    if isinstance(max_working_memory_gib, bool):
        raise TypeError("Schist max_working_memory_gib must be a number.")
    try:
        max_working_memory_gib = float(max_working_memory_gib)
    except (TypeError, ValueError) as error:
        raise TypeError("Schist max_working_memory_gib must be a number.") from error
    if not math.isfinite(max_working_memory_gib) or max_working_memory_gib <= 0:
        raise ValueError("Schist max_working_memory_gib must be finite and greater than zero.")
    budget_bytes = int(max_working_memory_gib * (1024**3))
    lower_bound = 65_536 + int(adata.n_obs) * 256 + int(adata.n_vars) * 128
    if lower_bound > budget_bytes:
        raise MemoryError(
            "Schist conservative memory lower bound exceeds max_working_memory_gib before graph materialization "
            f"(estimated at least {lower_bound} bytes; budget {budget_bytes} bytes)."
        )
    return {
        "random_seed": random_seed,
        "neighbors_key": neighbors_key,
        "key_added": key_added,
        "posterior_samples": posterior_samples,
        "degree_correction": degree_correction,
        "overwrite_existing": overwrite_existing,
        "max_working_memory_gib": max_working_memory_gib,
        "budget_bytes": budget_bytes,
    }


def _standalone_output_family(adata, *, key_added, overwrite_existing):
    import re
    from collections.abc import Mapping

    obs_pattern = re.compile(rf"^{re.escape(key_added)}_level_([0-9]+)$")
    obsm_pattern = re.compile(rf"^CM_{re.escape(key_added)}_level_([0-9]+)$")
    obs_prefix = f"{key_added}_level_"
    obsm_prefix = f"CM_{key_added}_level_"
    confusing = [f"obs[{key!r}]" for key in adata.obs if key.startswith(obs_prefix) and not obs_pattern.fullmatch(key)]
    confusing.extend(
        f"obsm[{key!r}]" for key in adata.obsm if key.startswith(obsm_prefix) and not obsm_pattern.fullmatch(key)
    )
    if confusing:
        raise ValueError(
            "Schist output namespace has noncanonical prefix collisions that the backend could delete: "
            + ", ".join(confusing)
        )
    obs_keys = tuple(key for key in adata.obs if obs_pattern.fullmatch(key))
    obsm_keys = tuple(key for key in adata.obsm if obsm_pattern.fullmatch(key))
    schist_root = adata.uns.get("schist")
    if schist_root is not None and not isinstance(schist_root, Mapping):
        raise ValueError("Schist requires uns['schist'] to be mapping-like when it already exists.")
    uns_collision = isinstance(schist_root, Mapping) and key_added in schist_root
    collisions = [*(f"obs[{key!r}]" for key in obs_keys), *(f"obsm[{key!r}]" for key in obsm_keys)]
    if uns_collision:
        collisions.append(f"uns['schist'][{key_added!r}]")
    if collisions and not overwrite_existing:
        raise ValueError(
            "Schist output family already exists; enable overwrite_existing to replace the complete owned family: "
            + ", ".join(collisions)
        )
    return {
        "obs_keys": obs_keys,
        "obsm_keys": obsm_keys,
        "uns_collision": bool(uns_collision),
        "collisions": tuple(collisions),
    }


def _standalone_object_bytes(value, *, np, pd, sparse, seen):
    import sys
    from collections.abc import Mapping

    object_id = id(value)
    if object_id in seen:
        return 0
    seen.add(object_id)
    if sparse.issparse(value):
        total = int(getattr(value.data, "nbytes", 0))
        for name in ("indices", "indptr", "row", "col"):
            array = getattr(value, name, None)
            total += int(getattr(array, "nbytes", 0))
        return total
    if isinstance(value, np.ndarray):
        return int(value.nbytes)
    if isinstance(value, pd.DataFrame):
        return int(value.memory_usage(index=True, deep=True).sum())
    if isinstance(value, pd.Series):
        return int(value.memory_usage(index=True, deep=True))
    if isinstance(value, Mapping):
        return int(sys.getsizeof(value)) + sum(
            _standalone_object_bytes(key, np=np, pd=pd, sparse=sparse, seen=seen)
            + _standalone_object_bytes(item, np=np, pd=pd, sparse=sparse, seen=seen)
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple, set, frozenset)):
        return int(sys.getsizeof(value)) + sum(
            _standalone_object_bytes(item, np=np, pd=pd, sparse=sparse, seen=seen) for item in value
        )
    return int(sys.getsizeof(value))


def _standalone_adata_bytes(adata, *, np, pd, sparse):
    seen = set()
    total = _standalone_object_bytes(adata.X, np=np, pd=pd, sparse=sparse, seen=seen)
    for frame in (adata.obs, adata.var):
        total += _standalone_object_bytes(frame, np=np, pd=pd, sparse=sparse, seen=seen)
    for container in (adata.layers, adata.obsm, adata.varm, adata.obsp, adata.varp, adata.uns):
        total += _standalone_object_bytes(container, np=np, pd=pd, sparse=sparse, seen=seen)
    if adata.raw is not None:
        total += _standalone_object_bytes(adata.raw.X, np=np, pd=pd, sparse=sparse, seen=seen)
        total += _standalone_object_bytes(adata.raw.var, np=np, pd=pd, sparse=sparse, seen=seen)
        total += _standalone_object_bytes(adata.raw.varm, np=np, pd=pd, sparse=sparse, seen=seen)
    return int(total)


def _standalone_memory_estimate(adata, *, raw_graph_nnz, posterior_samples, budget_bytes, np, pd, sparse):
    n_obs = int(adata.n_obs)
    n_vars = int(adata.n_vars)
    depth_bound = max(2, n_obs.bit_length() + 1)
    ann_data_owned = _standalone_adata_bytes(adata, np=np, pd=pd, sparse=sparse)
    graph_conversion = int(max(1, raw_graph_nnz) * 2 * 40 + (n_obs + 1) * 16)
    graph_tool_safety = int(n_obs * 512 + max(1, raw_graph_nnz // 2) * 512)
    posterior_hierarchies = int(posterior_samples) * n_obs * depth_bound * 8
    worst_case_marginals = 16 * n_obs * n_obs * depth_bound
    hierarchy_outputs = n_obs * depth_bound * 96
    validation_and_report = max(1_048_576, n_obs * 256 + n_vars * 128)
    h5ad_roundtrip_and_output_copy = 2 * ann_data_owned
    components = {
        "ann_data_owned_bytes": ann_data_owned,
        "output_copy_and_h5ad_roundtrip_bytes": h5ad_roundtrip_and_output_copy,
        "canonical_graph_and_coordinate_buffers_bytes": graph_conversion,
        "graph_tool_opaque_safety_allowance_bytes": graph_tool_safety,
        "posterior_hierarchy_samples_bytes": posterior_hierarchies,
        "worst_case_dense_marginals_bytes": worst_case_marginals,
        "hierarchy_outputs_bytes": hierarchy_outputs,
        "validation_and_report_buffers_bytes": validation_and_report,
    }
    estimated = int(sum(components.values()))
    result = {
        "estimated_peak_bytes": estimated,
        "budget_bytes": int(budget_bytes),
        "estimated_peak_gib": float(estimated / (1024**3)),
        "budget_gib": float(budget_bytes / (1024**3)),
        "hierarchy_depth_bound": depth_bound,
        "graph_tool_allocations_are_opaque": True,
        "is_process_rss_ceiling": False,
        "components": components,
    }
    if estimated > budget_bytes:
        raise MemoryError(
            "Schist conservative working-memory estimate exceeds max_working_memory_gib before canonical graph, "
            f"AnnData copy, or backend work (estimated {estimated} bytes; budget {budget_bytes} bytes; "
            f"components={components})."
        )
    return result


def _standalone_numeric_summary(values, *, np):
    array = np.asarray(values, dtype=float).ravel()
    if not array.size or not bool(np.isfinite(array).all()):
        raise RuntimeError("Schist diagnostic values must be finite and nonempty.")
    quantiles = np.quantile(array, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "n": int(array.size),
        "min": float(quantiles[0]),
        "q1": float(quantiles[1]),
        "median": float(quantiles[2]),
        "mean": float(array.mean()),
        "q3": float(quantiles[3]),
        "max": float(quantiles[4]),
    }


def _standalone_graph_fingerprint(matrix, obs_names, *, np):
    import hashlib

    digest = hashlib.sha256()
    digest.update(b"openbio-schist-named-graph-v1\0")
    digest.update(int(matrix.shape[0]).to_bytes(8, "little", signed=False))
    for name in obs_names:
        encoded = str(name).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little", signed=False))
        digest.update(encoded)
    for array, dtype in (
        (matrix.indptr, "<i8"),
        (matrix.indices, "<i8"),
        (matrix.data, "<f8"),
    ):
        values = np.asarray(array, dtype=dtype)
        for start in range(0, int(values.size), 65_536):
            digest.update(values[start : start + 65_536].tobytes(order="C"))
    return digest.hexdigest()


def _standalone_validate_graph(adata, *, neighbors_key, posterior_samples, budget_bytes, np, pd, sparse):
    from collections.abc import Mapping

    from scipy.sparse.csgraph import connected_components

    metadata = adata.uns.get(neighbors_key)
    if not isinstance(metadata, Mapping):
        raise ValueError(f"Schist neighbor metadata not found in uns: {neighbors_key!r}.")
    connectivity_key = metadata.get("connectivities_key")
    if not isinstance(connectivity_key, str) or not connectivity_key.strip():
        raise ValueError(f"Schist uns[{neighbors_key!r}] has no valid connectivities_key.")
    connectivity_key = connectivity_key.strip()
    if connectivity_key not in adata.obsp:
        raise ValueError(f"Schist connectivity matrix not found in obsp: {connectivity_key!r}.")
    raw = adata.obsp[connectivity_key]
    if not sparse.issparse(raw):
        raise TypeError("Schist connectivity graph must be a SciPy sparse matrix.")
    if tuple(raw.shape) != (int(adata.n_obs), int(adata.n_obs)):
        raise ValueError("Schist connectivity graph must have shape (n_obs, n_obs).")
    raw_values = np.asarray(raw.data)
    if raw_values.size and (
        not bool(np.issubdtype(raw_values.dtype, np.number)) or bool(np.iscomplexobj(raw_values))
    ):
        raise TypeError("Schist connectivity graph must contain real numeric weights.")
    if raw_values.size and not bool(np.isfinite(raw_values).all()):
        raise ValueError("Schist connectivity graph contains non-finite weights.")
    if raw_values.size and bool((raw_values < 0).any()):
        raise ValueError("Schist connectivity graph contains negative weights.")
    memory = _standalone_memory_estimate(
        adata,
        raw_graph_nnz=int(raw.nnz),
        posterior_samples=posterior_samples,
        budget_bytes=budget_bytes,
        np=np,
        pd=pd,
        sparse=sparse,
    )
    matrix = sparse.csr_matrix(raw, dtype=np.float64, copy=True)
    matrix.sum_duplicates()
    matrix.sort_indices()
    matrix.eliminate_zeros()
    diagonal = np.asarray(matrix.diagonal(), dtype=float)
    if diagonal.size and not bool(np.allclose(diagonal, 0.0, rtol=0.0, atol=1e-12)):
        raise ValueError("Schist connectivity graph must have a zero diagonal.")
    delta = matrix - matrix.T
    delta.eliminate_zeros()
    if delta.nnz and float(np.max(np.abs(delta.data))) > 1e-8:
        raise ValueError("Schist connectivity graph must be symmetric within absolute tolerance 1e-8.")
    matrix = sparse.csr_matrix((matrix + matrix.T) * 0.5, dtype=np.float64)
    matrix.sum_duplicates()
    matrix.sort_indices()
    matrix.eliminate_zeros()
    if int(matrix.nnz) == 0:
        raise ValueError("Schist connectivity graph contains no positive off-diagonal edges.")
    degree = np.diff(matrix.indptr).astype(np.int64, copy=False)
    isolated = np.flatnonzero(degree == 0)
    if isolated.size:
        examples = [str(adata.obs_names[int(index)]) for index in isolated[:20]]
        raise ValueError(
            f"Schist connectivity graph contains {int(isolated.size)} isolated observations: {examples}."
        )
    n_components, labels = connected_components(matrix, directed=False, return_labels=True)
    component_sizes = np.bincount(labels, minlength=n_components)
    ordered_component_sizes = [int(value) for value in sorted(component_sizes, reverse=True)]
    positive_edges = int(matrix.nnz // 2)
    possible_edges = int(adata.n_obs) * (int(adata.n_obs) - 1) // 2
    graph = {
        "neighbors_key": neighbors_key,
        "connectivities_key": connectivity_key,
        "fingerprint_sha256": _standalone_graph_fingerprint(matrix, adata.obs_names, np=np),
        "shape": [int(adata.n_obs), int(adata.n_obs)],
        "storage": "canonical_symmetric_csr_float64",
        "symmetric_absolute_tolerance": 1e-8,
        "zero_diagonal_absolute_tolerance": 1e-12,
        "nonzero_entries": int(matrix.nnz),
        "positive_undirected_edges": positive_edges,
        "density": float(positive_edges / possible_edges) if possible_edges else 0.0,
        "connected_components": int(n_components),
        "component_sizes": ordered_component_sizes[:100],
        "component_sizes_total": len(ordered_component_sizes),
        "component_sizes_truncated": len(ordered_component_sizes) > 100,
        "component_sizes_record_limit": 100,
        "isolated_observation_count": 0,
        "degree": _standalone_numeric_summary(degree, np=np),
        "connectivity_weights": _standalone_numeric_summary(matrix.data, np=np),
        "upstream_neighbor_parameters": _standalone_plain_json(metadata.get("params", {})),
        "connectivity_weights_used_by_model": False,
    }
    return matrix, graph, memory


def _standalone_require_backend():
    import importlib
    import inspect as inspect_module

    guidance = (
        "Use the official dawe/schist v0.10.0 source or conda-forge Schist together with graph-tool; "
        "see https://github.com/dawe/schist and https://graph-tool.skewed.de/installation.html. "
        "The same-named PyPI knowledge-graph project is not compatible. On native Windows, use WSL or a supported "
        "Linux container/environment for graph-tool."
    )
    try:
        schist = importlib.import_module("schist")
    except (ImportError, OSError) as error:
        raise RuntimeError(f"Official single-cell Schist v0.10.0 is unavailable ({error}). {guidance}") from error
    module_path = str(getattr(schist, "__file__", "unknown"))
    detected_version = str(getattr(schist, "__version__", "unknown"))
    author = str(getattr(schist, "__author__", ""))
    if "Davide Cittaro" not in author or "Leonardo Morelli" not in author:
        raise RuntimeError(
            "The imported 'schist' module does not expose the official single-cell Schist authorship "
            f"(path={module_path!r}, version={detected_version!r}). {guidance}"
        )
    try:
        from packaging.version import InvalidVersion, Version
    except ImportError as error:
        raise RuntimeError("Version verification requires the public packaging library. " + guidance) from error
    try:
        parsed_version = Version(detected_version)
        expected_version = Version("0.10.0")
    except InvalidVersion as error:
        raise RuntimeError(
            f"Cannot verify the installed Schist version (path={module_path!r}, version={detected_version!r}). {guidance}"
        ) from error
    if parsed_version != expected_version:
        raise RuntimeError(
            "Schist compatibility is fail-closed at the audited exact release 0.10.0 "
            f"(path={module_path!r}, detected={detected_version!r}). {guidance}"
        )
    inference = getattr(schist, "inference", None)
    fit_model = getattr(inference, "fit_model", None)
    if not callable(fit_model):
        raise RuntimeError(
            f"Official Schist fit_model is unavailable (path={module_path!r}, version={detected_version!r}). {guidance}"
        )
    expected_parameters = (
        "adata",
        "nested",
        "assortative",
        "collect_marginals",
        "n_samples",
        "key_added",
        "adjacency",
        "neighbors_key",
        "constraint_key",
        "deg_corr",
        "directed",
        "use_weights",
        "bisection",
        "simple_init",
        "n_jobs",
        "n_iter",
        "beta",
        "save_model",
        "copy",
        "random_seed",
    )
    try:
        signature = inspect_module.signature(fit_model)
        detected_parameters = tuple(signature.parameters)
    except (TypeError, ValueError) as error:
        raise RuntimeError("Cannot inspect the installed Schist fit_model interface. " + guidance) from error
    if detected_parameters != expected_parameters:
        raise RuntimeError(
            "Installed Schist fit_model does not match the exact audited v0.10.0 interface "
            f"(detected parameters={detected_parameters!r}). {guidance}"
        )
    expected_defaults = {
        "adata": inspect_module.Parameter.empty,
        "nested": True,
        "assortative": False,
        "collect_marginals": True,
        "n_samples": 100,
        "key_added": None,
        "adjacency": None,
        "neighbors_key": "neighbors",
        "constraint_key": None,
        "deg_corr": True,
        "directed": False,
        "use_weights": False,
        "bisection": True,
        "simple_init": False,
        "n_jobs": -1,
        "n_iter": 10,
        "beta": 1.0,
        "save_model": None,
        "copy": False,
        "random_seed": None,
    }
    incompatible = [
        name
        for name, parameter in signature.parameters.items()
        if parameter.kind is not inspect_module.Parameter.POSITIONAL_OR_KEYWORD
        or parameter.default != expected_defaults[name]
    ]
    if incompatible:
        raise RuntimeError(
            "Installed Schist fit_model parameter kinds/defaults do not match the exact audited v0.10.0 interface "
            f"(incompatible={incompatible!r}). {guidance}"
        )
    try:
        graph_tool = importlib.import_module("graph_tool.all")
    except (ImportError, OSError) as error:
        raise RuntimeError(f"Schist requires graph-tool in a supported environment ({error}). {guidance}") from error
    graph_tool_version = str(getattr(graph_tool, "__version__", "unknown"))
    if not graph_tool_version or graph_tool_version == "unknown":
        raise RuntimeError("The imported graph-tool backend exposes no version information. " + guidance)
    return {
        "schist": schist,
        "fit_model": fit_model,
        "schist_version": detected_version,
        "schist_module_path": module_path,
        "schist_author": author,
        "graph_tool": graph_tool,
        "graph_tool_version": graph_tool_version,
        "graph_tool_module_path": str(getattr(graph_tool, "__file__", "unknown")),
    }


def _standalone_hash_update(hasher, value, *, np, pd, sparse):
    from collections.abc import Mapping

    def token(text):
        encoded = text.encode("utf-8", errors="backslashreplace")
        hasher.update(len(encoded).to_bytes(8, "little", signed=False))
        hasher.update(encoded)

    if value is None or isinstance(value, (str, bytes, bool, int, float, complex)):
        token(f"scalar:{type(value).__qualname__}:{value!r}")
        return
    if isinstance(value, np.generic):
        token(f"numpy-scalar:{value.dtype}:{value!r}")
        return
    if sparse.issparse(value):
        token(f"sparse:{type(value).__qualname__}:{value.format}:{value.shape}:{value.dtype}")
        matrix = value if value.format in {"csr", "csc", "coo"} else value.tocoo(copy=False)
        for name in ("data", "indices", "indptr", "row", "col"):
            array = getattr(matrix, name, None)
            if array is not None:
                token(name)
                _standalone_hash_update(hasher, np.asarray(array), np=np, pd=pd, sparse=sparse)
        return
    if isinstance(value, np.ndarray):
        token(f"ndarray:{value.shape}:{value.dtype}")
        if value.dtype.kind in "OSUV":
            for item in value.flat:
                _standalone_hash_update(hasher, item, np=np, pd=pd, sparse=sparse)
        else:
            iterator = np.nditer(
                value,
                flags=["external_loop", "buffered", "zerosize_ok", "refs_ok"],
                op_flags=["readonly"],
                order="C",
                buffersize=65_536,
            )
            for chunk in iterator:
                hasher.update(np.asarray(chunk).tobytes(order="C"))
        return
    if isinstance(value, pd.DataFrame):
        token("dataframe")
        _standalone_hash_update(hasher, value.index, np=np, pd=pd, sparse=sparse)
        for position, column in enumerate(value.columns):
            _standalone_hash_update(hasher, column, np=np, pd=pd, sparse=sparse)
            _standalone_hash_update(hasher, value.iloc[:, position], np=np, pd=pd, sparse=sparse)
        return
    if isinstance(value, pd.Series):
        token(f"series:{value.name!r}:{value.dtype}")
        _standalone_hash_update(hasher, value.index, np=np, pd=pd, sparse=sparse)
        if isinstance(value.dtype, pd.CategoricalDtype):
            _standalone_hash_update(hasher, value.cat.categories, np=np, pd=pd, sparse=sparse)
            _standalone_hash_update(hasher, value.cat.codes.to_numpy(), np=np, pd=pd, sparse=sparse)
        else:
            for item in value.array:
                _standalone_hash_update(hasher, item, np=np, pd=pd, sparse=sparse)
        return
    if isinstance(value, pd.Index):
        token(f"index:{type(value).__qualname__}:{value.dtype}:{value.name!r}")
        for item in value:
            _standalone_hash_update(hasher, item, np=np, pd=pd, sparse=sparse)
        return
    if isinstance(value, Mapping):
        token(f"mapping:{type(value).__qualname__}:{len(value)}")
        for key in sorted(value, key=lambda item: (type(item).__qualname__, repr(item))):
            _standalone_hash_update(hasher, key, np=np, pd=pd, sparse=sparse)
            _standalone_hash_update(hasher, value[key], np=np, pd=pd, sparse=sparse)
        return
    if isinstance(value, (list, tuple)):
        token(f"sequence:{type(value).__qualname__}:{len(value)}")
        for item in value:
            _standalone_hash_update(hasher, item, np=np, pd=pd, sparse=sparse)
        return
    token(f"fallback:{type(value).__module__}.{type(value).__qualname__}:{value!r}")


def _standalone_state_fingerprint(adata, *, key_added, np, pd, sparse):
    import hashlib
    import re
    from collections.abc import Mapping

    obs_pattern = re.compile(rf"^{re.escape(key_added)}_level_[0-9]+$")
    obsm_pattern = re.compile(rf"^CM_{re.escape(key_added)}_level_[0-9]+$")
    hasher = hashlib.sha256()
    hasher.update(b"openbio-schist-untouched-state-v1\0")
    _standalone_hash_update(hasher, adata.X, np=np, pd=pd, sparse=sparse)
    _standalone_hash_update(
        hasher,
        adata.obs[[column for column in adata.obs if not obs_pattern.fullmatch(column)]],
        np=np,
        pd=pd,
        sparse=sparse,
    )
    _standalone_hash_update(hasher, adata.var, np=np, pd=pd, sparse=sparse)
    for name, container in (
        ("layers", adata.layers),
        ("obsm", {key: value for key, value in adata.obsm.items() if not obsm_pattern.fullmatch(key)}),
        ("varm", adata.varm),
        ("obsp", adata.obsp),
        ("varp", adata.varp),
    ):
        _standalone_hash_update(hasher, name, np=np, pd=pd, sparse=sparse)
        _standalone_hash_update(hasher, container, np=np, pd=pd, sparse=sparse)
    uns = dict(adata.uns)
    schist_root = uns.get("schist")
    if isinstance(schist_root, Mapping):
        filtered = {key: value for key, value in schist_root.items() if key != key_added}
        if filtered:
            uns["schist"] = filtered
        else:
            uns.pop("schist", None)
    _standalone_hash_update(hasher, "uns", np=np, pd=pd, sparse=sparse)
    _standalone_hash_update(hasher, uns, np=np, pd=pd, sparse=sparse)
    if adata.raw is None:
        _standalone_hash_update(hasher, None, np=np, pd=pd, sparse=sparse)
    else:
        _standalone_hash_update(hasher, adata.raw.X, np=np, pd=pd, sparse=sparse)
        _standalone_hash_update(hasher, adata.raw.var, np=np, pd=pd, sparse=sparse)
        _standalone_hash_update(hasher, adata.raw.varm, np=np, pd=pd, sparse=sparse)
    return hasher.hexdigest()


def _standalone_clear_output_family(output, *, key_added, family):
    for key in family["obs_keys"]:
        del output.obs[key]
    for key in family["obsm_keys"]:
        del output.obsm[key]
    if family["uns_collision"]:
        del output.uns["schist"][key_added]


def _standalone_validate_backend_output(output, *, key_added, graph, settings, source_fingerprint, np, pd, sparse):
    import re
    from collections.abc import Mapping

    obs_pattern = re.compile(rf"^{re.escape(key_added)}_level_([0-9]+)$")
    obsm_pattern = re.compile(rf"^CM_{re.escape(key_added)}_level_([0-9]+)$")
    obs_levels = {int(match.group(1)): key for key in output.obs if (match := obs_pattern.fullmatch(key))}
    obsm_levels = {int(match.group(1)): key for key in output.obsm if (match := obsm_pattern.fullmatch(key))}
    if not obs_levels:
        raise RuntimeError("Schist backend returned no nested hierarchy levels.")
    expected_levels = list(range(max(obs_levels) + 1))
    if sorted(obs_levels) != expected_levels:
        raise RuntimeError("Schist backend hierarchy levels must be consecutive from level 0 through the root.")
    if sorted(obsm_levels) != expected_levels:
        raise RuntimeError("Schist backend must return one consecutive marginal matrix for every hierarchy level.")

    memberships = []
    cluster_counts = []
    cluster_sizes = []
    for level in expected_levels:
        series = output.obs[obs_levels[level]]
        if not series.index.equals(output.obs_names) or len(series) != int(output.n_obs):
            raise RuntimeError("Schist backend hierarchy membership does not align to observation order.")
        if not isinstance(series.dtype, pd.CategoricalDtype):
            raise RuntimeError("Schist backend hierarchy memberships must be categorical.")
        if bool(series.isna().any()):
            raise RuntimeError("Schist backend hierarchy memberships contain missing values.")
        categories = [str(value) for value in series.cat.categories]
        if categories != [str(index) for index in range(len(categories))]:
            raise RuntimeError("Schist backend hierarchy categories must be canonical consecutive strings.")
        codes = series.cat.codes.to_numpy(dtype=np.int64, copy=True)
        counts = np.bincount(codes, minlength=len(categories))
        if len(categories) < 1 or bool((counts <= 0).any()) or int(counts.sum()) != int(output.n_obs):
            raise RuntimeError("Schist backend hierarchy has unused categories or incomplete memberships.")
        memberships.append(codes)
        cluster_counts.append(len(categories))
        cluster_sizes.append(counts)
    if any(right > left for left, right in zip(cluster_counts, cluster_counts[1:], strict=False)):
        raise RuntimeError("Schist backend hierarchy must be ordered from finest to nonincreasing coarser levels.")
    if cluster_counts[-1] != 1:
        raise RuntimeError("Schist backend hierarchy must end in exactly one root block.")
    nested_to_next = []
    parent_maps = []
    for level in range(len(expected_levels) - 1):
        parent = np.full(cluster_counts[level], -1, dtype=np.int64)
        for fine, coarse in zip(memberships[level], memberships[level + 1], strict=True):
            if parent[fine] == -1:
                parent[fine] = coarse
            elif parent[fine] != coarse:
                raise RuntimeError("Schist backend hierarchy violates the nested parent mapping invariant.")
        if bool((parent < 0).any()):
            raise RuntimeError("Schist backend hierarchy has a fine block without a parent.")
        parent_maps.append(parent)
        nested_to_next.append(True)

    schist_root = output.uns.get("schist")
    bundle = schist_root.get(key_added) if isinstance(schist_root, Mapping) else None
    if not isinstance(bundle, Mapping):
        raise RuntimeError("Schist backend did not return uns['schist'][key_added].")
    stats = bundle.get("stats")
    params = bundle.get("params")
    blocks = bundle.get("blocks")
    if not isinstance(stats, Mapping) or not isinstance(params, Mapping) or not isinstance(blocks, Mapping):
        raise RuntimeError("Schist backend metadata requires mapping-like stats, params, and blocks entries.")
    try:
        entropy = float(stats.get("entropy"))
    except (TypeError, ValueError) as error:
        raise RuntimeError("Schist backend total entropy is missing or invalid.") from error
    modularity = np.asarray(stats.get("modularity"), dtype=float).ravel()
    level_entropy = np.asarray(stats.get("level_entropy"), dtype=float).ravel()
    if not np.isfinite(entropy) or entropy < 0:
        raise RuntimeError("Schist backend total entropy must be finite and non-negative.")
    if modularity.shape != (len(expected_levels),) or not bool(np.isfinite(modularity).all()):
        raise RuntimeError("Schist backend modularity must be finite and aligned to every hierarchy level.")
    if bool(((modularity < -1.0) | (modularity > 1.0)).any()):
        raise RuntimeError("Schist backend modularity values must be within [-1, 1].")
    if level_entropy.shape != (len(expected_levels),) or not bool(np.isfinite(level_entropy).all()):
        raise RuntimeError("Schist backend level entropy must be finite and aligned to every hierarchy level.")

    required_params = {
        "nested": True,
        "assortative": False,
        "neighbors_key": settings["neighbors_key"],
        "use_weights": False,
        "key_added": key_added,
        "n_samples": settings["posterior_samples"],
        "collect_marginals": True,
        "random_seed": settings["random_seed"],
        "deg_corr": settings["degree_correction"],
        "directed": False,
        "n_iter": 10,
        "beta": 1.0,
    }
    for name, expected in required_params.items():
        if name not in params or params[name] != expected:
            raise RuntimeError(f"Schist backend parameter metadata is inconsistent for {name!r}.")

    if set(blocks) != {str(level) for level in expected_levels}:
        raise RuntimeError("Schist backend block arrays must be present for every consecutive hierarchy level.")
    for level in expected_levels:
        array = np.asarray(blocks[str(level)])
        expected_length = int(output.n_obs) if level == 0 else cluster_counts[level - 1]
        if array.ndim != 1 or int(array.size) != expected_length or not np.issubdtype(array.dtype, np.integer):
            raise RuntimeError("Schist backend block arrays have invalid shape or dtype.")
        expected_codes = memberships[0] if level == 0 else parent_maps[level - 1]
        if not bool(np.array_equal(array.astype(np.int64, copy=False), expected_codes)):
            raise RuntimeError("Schist backend block arrays are inconsistent with hierarchy memberships.")

    level_records = []
    weak_mass = False
    for level in expected_levels:
        matrix = np.asarray(output.obsm[obsm_levels[level]])
        expected_shape = (int(output.n_obs), cluster_counts[level])
        if matrix.shape != expected_shape:
            raise RuntimeError("Schist backend marginal matrix shape does not match observations and used blocks.")
        if not np.issubdtype(matrix.dtype, np.number) or np.iscomplexobj(matrix):
            raise RuntimeError("Schist backend marginal matrices must be real numeric.")
        values = np.asarray(matrix, dtype=float)
        if not bool(np.isfinite(values).all()) or bool(((values < 0) | (values > 1)).any()):
            raise RuntimeError("Schist backend marginal matrices must contain finite probabilities in [0, 1].")
        row_sums = values.sum(axis=1)
        if bool((row_sums <= 0).any()) or bool((row_sums > 1.0 + 1e-8).any()):
            raise RuntimeError("Schist backend marginal rows must have positive represented mass no greater than one.")
        maximum = values.max(axis=1)
        weak_mass = weak_mass or bool((row_sums < 0.9).any())
        counts = cluster_sizes[level]
        ordered_counts = [(str(index), int(count)) for index, count in enumerate(counts)]
        size_quantiles = np.quantile(counts.astype(float), [0.0, 0.5, 1.0])
        level_records.append(
            {
                "level": level,
                "obs_key": obs_levels[level],
                "obsm_key": obsm_levels[level],
                "is_finest": level == 0,
                "is_root": level == expected_levels[-1],
                "cluster_count": cluster_counts[level],
                "cluster_sizes": dict(ordered_counts[:50]),
                "cluster_sizes_total": len(ordered_counts),
                "cluster_sizes_truncated": len(ordered_counts) > 50,
                "cluster_sizes_record_limit": 50,
                "minimum_cluster_size": int(size_quantiles[0]),
                "median_cluster_size": float(size_quantiles[1]),
                "maximum_cluster_size": int(size_quantiles[2]),
                "singleton_count": int((counts == 1).sum()),
                "largest_cluster_fraction": float(counts.max() / int(output.n_obs)),
                "modularity": float(modularity[level]),
                "level_entropy": float(level_entropy[level]),
                "marginal_shape": [int(value) for value in matrix.shape],
                "marginal_row_sum": _standalone_numeric_summary(row_sums, np=np),
                "maximum_membership_probability": _standalone_numeric_summary(maximum, np=np),
                "nested_to_next_level_passed": True if level < expected_levels[-1] else None,
            }
        )

    observed_fingerprint = _standalone_state_fingerprint(
        output,
        key_added=key_added,
        np=np,
        pd=pd,
        sparse=sparse,
    )
    if observed_fingerprint != source_fingerprint:
        raise RuntimeError("Schist backend changed input state outside its owned hierarchy output family.")
    graph_after = output.obsp.get(graph["connectivities_key"])
    if graph_after is None:
        raise RuntimeError("Schist backend removed the named input graph from obsp.")
    return {
        "entropy": entropy,
        "cluster_counts": cluster_counts,
        "level_records": level_records,
        "root_level": expected_levels[-1],
        "informative_level_count": sum(count > 1 for count in cluster_counts),
        "weak_marginal_mass": weak_mass,
        "obs_keys": [obs_levels[level] for level in expected_levels],
        "obsm_keys": [obsm_levels[level] for level in expected_levels],
    }


def _standalone_h5ad_roundtrip(output, *, key_added, expected_obs_keys, expected_obsm_keys):
    import tempfile
    from pathlib import Path

    import anndata as ad

    with tempfile.TemporaryDirectory(prefix="openbio-schist-") as directory:
        path = Path(directory) / "result.h5ad"
        output.write_h5ad(path, convert_strings_to_categoricals=False)
        restored = ad.read_h5ad(path)
        if tuple(restored.shape) != tuple(output.shape):
            raise RuntimeError("Schist output failed H5AD round-trip shape validation.")
        if any(key not in restored.obs for key in expected_obs_keys):
            raise RuntimeError("Schist output lost hierarchy memberships during H5AD round-trip.")
        if any(key not in restored.obsm for key in expected_obsm_keys):
            raise RuntimeError("Schist output lost marginal matrices during H5AD round-trip.")
        schist_root = restored.uns.get("schist")
        if not hasattr(schist_root, "get") or key_added not in schist_root:
            raise RuntimeError("Schist output lost backend metadata during H5AD round-trip.")


def _standalone_backend_build(backend):
    graph_tool = backend["graph_tool"]

    def resolved(name):
        value = getattr(graph_tool, name, None)
        if callable(value):
            try:
                value = value()
            except Exception as error:
                return f"unavailable:{type(error).__name__}"
        return _standalone_plain_json(value)

    return {
        "schist_module_path": backend["schist_module_path"],
        "schist_author": backend["schist_author"],
        "graph_tool_module_path": backend["graph_tool_module_path"],
        "graph_tool_openmp_enabled": resolved("openmp_enabled"),
        "graph_tool_detected_threads": resolved("openmp_get_num_threads"),
        "analysis_threads": 1,
    }


def _standalone_software_versions(*, openbio_version, backend):
    import platform
    from importlib import metadata

    def version(distribution):
        try:
            return metadata.version(distribution)
        except metadata.PackageNotFoundError:
            return "not-installed"

    return {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
        "schist": backend["schist_version"],
        "graph-tool": backend["graph_tool_version"],
        "scanpy": version("scanpy"),
        "anndata": version("anndata"),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "scipy": version("scipy"),
        "platform": platform.platform(),
        "python_implementation": platform.python_implementation(),
    }


def _standalone_summary(
    *,
    adata,
    settings,
    graph,
    memory,
    backend,
    hierarchy,
    overwritten,
    runtime_duration_seconds,
    openbio_version,
):
    import json

    warnings = []
    if graph["connected_components"] > 1:
        warnings.append(
            "The selected named graph is disconnected; components were preserved and modeled without artificial edges."
        )
    if overwritten:
        warnings.append("A pre-existing Schist-owned output family was removed on the private output copy before fitting.")
    if hierarchy["weak_marginal_mass"]:
        warnings.append(
            "At least one marginal row represents less than 0.9 total probability mass; Schist documents that sampled "
            "partitions with different block counts can leave represented mass below one."
        )
    if any(record["largest_cluster_fraction"] >= 0.9 and not record["is_root"] for record in hierarchy["level_records"]):
        warnings.append("A non-root hierarchy level is dominated by one block containing at least 90% of observations.")
    if any(record["singleton_count"] > 0 for record in hierarchy["level_records"]):
        warnings.append("At least one hierarchy level contains singleton blocks; interpret them cautiously.")

    parameters = {
        "neighbors_key": settings["neighbors_key"],
        "key_added": settings["key_added"],
        "posterior_samples": settings["posterior_samples"],
        "degree_correction": settings["degree_correction"],
        "overwrite_existing": settings["overwrite_existing"],
        "max_working_memory_gib": settings["max_working_memory_gib"],
        "random_seed": settings["random_seed"],
        "nested": True,
        "assortative": False,
        "collect_marginals": True,
        "constraint_key": None,
        "directed": False,
        "connectivity_weights_used_by_model": False,
        "bisection": True,
        "simple_init": False,
        "n_jobs": 1,
        "mcmc_sweeps_per_posterior_sample": 10,
        "beta": 1.0,
        "save_model": None,
        "backend_copy": False,
        "automatic_level_selection": False,
    }
    limitations = [
        "This is exploratory cell-neighborhood graph clustering, not a differential abundance or differential composition test.",
        "The inferred blocks are cluster evidence, not Curated annotation or biological ground truth.",
        "The hierarchy depends on the upstream graph and can reflect Sample or Technical batch structure present in that graph.",
        "Cells are graph vertices rather than independent biological replicates; no Sample-aware Condition inference is performed.",
        "Schist v0.10.0 exposes no ESS, R-hat, posterior trace, or formal convergence status for this workflow.",
        "Repeatability is limited to the same graph bytes, Schist and graph-tool builds, platform, parameters, and nonzero seed.",
        "The memory estimate is conservative but cannot guarantee a process RSS ceiling because graph-tool allocates opaque C++ state.",
    ]
    references = [
        {
            "citation": "Morelli L, Giansanti V, Cittaro D. Nested Stochastic Block Models applied to single cell data. BMC Bioinformatics. 2021;22:576.",
            "doi": "10.1186/s12859-021-04489-7",
            "url": "https://doi.org/10.1186/s12859-021-04489-7",
            "kind": "method",
        },
        {
            "citation": "Peixoto TP. Hierarchical Block Structures and High-Resolution Model Selection in Large Networks. Physical Review X. 2014;4:011047.",
            "doi": "10.1103/PhysRevX.4.011047",
            "url": "https://doi.org/10.1103/PhysRevX.4.011047",
            "kind": "method",
        },
        {
            "citation": "Peixoto TP. Nonparametric Bayesian inference of the microcanonical stochastic block model. Physical Review E. 2017;95:012317.",
            "doi": "10.1103/PhysRevE.95.012317",
            "url": "https://doi.org/10.1103/PhysRevE.95.012317",
            "kind": "method",
        },
        {
            "citation": "Peixoto TP. The graph-tool Python library. figshare. 2014.",
            "doi": "10.6084/m9.figshare.1164194",
            "url": "https://doi.org/10.6084/m9.figshare.1164194",
            "kind": "software",
        },
        {
            "citation": "Wolf FA, Angerer P, Theis FJ. SCANPY. Genome Biology. 2018;19:15.",
            "doi": "10.1186/s13059-017-1382-0",
            "url": "https://doi.org/10.1186/s13059-017-1382-0",
            "kind": "software",
        },
        {
            "citation": "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. JOSS. 2021;6:4371.",
            "doi": "10.21105/joss.04371",
            "url": "https://doi.org/10.21105/joss.04371",
            "kind": "software",
        },
        {
            "citation": "Schist v0.10.0 fit_model official source.",
            "url": "https://github.com/dawe/schist/blob/v0.10.0/schist/inference/_model.py",
            "kind": "software_documentation",
        },
    ]
    cluster_counts = hierarchy["cluster_counts"]
    results = (
        f"Schist v0.10.0 completed one unweighted undirected general nested-SBM fit on "
        f"{graph['positive_undirected_edges']} graph edges and returned {len(cluster_counts)} hierarchy level(s), "
        f"with block counts {cluster_counts} from finest through the one-block root. All registered postconditions passed; "
        "no hierarchy level was selected automatically."
    )
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellSchistNestedModel",
        "title": "Schist Nested-SBM Hierarchy",
        "operation": "schist_nested_sbm_hierarchy",
        "analysis_status": "exploratory_graph_clustering",
        "input_dimensions": {"observations": int(adata.n_obs), "variables": int(adata.n_vars)},
        "methods": (
            "The exact official Schist 0.10.0 fit_model interface fit a degree-corrected or uncorrected general "
            "nested stochastic block model to the explicitly named canonical symmetric sparse graph. The fixed policy "
            "used nested=True, assortative=False, collect_marginals=True, directed=False, use_weights=False, "
            "bisection=True, simple_init=False, n_jobs=1, n_iter=10, beta=1.0, save_model=None, copy=False, and a "
            "nonzero graph-tool seed. Every inferred level and its marginal matrix was retained."
        ),
        "results": results,
        "key_results": {
            "hierarchy_levels": hierarchy["level_records"],
            "cluster_counts_finest_to_root": cluster_counts,
            "total_entropy": hierarchy["entropy"],
            "informative_level_count": hierarchy["informative_level_count"],
            "root_level": hierarchy["root_level"],
            "output_storage": {
                "obs_membership_keys": hierarchy["obs_keys"],
                "obsm_marginal_keys": hierarchy["obsm_keys"],
                "uns_key": f"schist/{settings['key_added']}",
                "new_obsp_keys": [],
            },
        },
        "graph": graph,
        "parameters": parameters,
        "memory_preflight": memory,
        "execution": {
            "runtime_duration_seconds": float(runtime_duration_seconds),
            "duration_clock": "time.perf_counter",
            "backend_threads": 1,
        },
        "completion": {
            "backend_returned_and_postconditions_passed": True,
            "minimization_completed": True,
            "posterior_samples_collected": settings["posterior_samples"],
            "convergence_diagnostic_available": False,
            "convergence_explanation": (
                "Completion means the description-length minimization and fixed posterior sampling returned and passed "
                "storage invariants; Schist exposes no ESS, R-hat, trace, or formal convergence flag here."
            ),
            "h5ad_roundtrip_passed": True,
        },
        "warnings": warnings,
        "limitations": limitations,
        "references": references,
        "software_versions": _standalone_software_versions(openbio_version=openbio_version, backend=backend),
        "backend_build": _standalone_backend_build(backend),
    }
    summary = _standalone_plain_json(summary)
    json.dumps(summary, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)
    return summary


def _standalone_schist_nested_model(
    adata,
    *,
    random_seed=123,
    neighbors_key="neighbors",
    key_added="nsbm",
    posterior_samples=100,
    degree_correction=True,
    overwrite_existing=False,
    max_working_memory_gib=8.0,
    openbio_version="unknown",
):
    import time

    import numpy as np
    import pandas as pd
    from scipy import sparse

    started_at = time.perf_counter()

    settings = _standalone_validate_settings(
        adata,
        random_seed=random_seed,
        neighbors_key=neighbors_key,
        key_added=key_added,
        posterior_samples=posterior_samples,
        degree_correction=degree_correction,
        overwrite_existing=overwrite_existing,
        max_working_memory_gib=max_working_memory_gib,
    )
    family = _standalone_output_family(
        adata,
        key_added=settings["key_added"],
        overwrite_existing=settings["overwrite_existing"],
    )
    adjacency, graph, memory = _standalone_validate_graph(
        adata,
        neighbors_key=settings["neighbors_key"],
        posterior_samples=settings["posterior_samples"],
        budget_bytes=settings["budget_bytes"],
        np=np,
        pd=pd,
        sparse=sparse,
    )
    backend = _standalone_require_backend()
    source_fingerprint = _standalone_state_fingerprint(
        adata,
        key_added=settings["key_added"],
        np=np,
        pd=pd,
        sparse=sparse,
    )
    adjacency_fingerprint = _standalone_graph_fingerprint(adjacency, adata.obs_names, np=np)
    output = adata
    _standalone_clear_output_family(output, key_added=settings["key_added"], family=family)
    result = backend["fit_model"](
        output,
        nested=True,
        assortative=False,
        collect_marginals=True,
        n_samples=settings["posterior_samples"],
        key_added=settings["key_added"],
        adjacency=adjacency,
        neighbors_key=settings["neighbors_key"],
        constraint_key=None,
        deg_corr=settings["degree_correction"],
        directed=False,
        use_weights=False,
        bisection=True,
        simple_init=False,
        n_jobs=1,
        n_iter=10,
        beta=1.0,
        save_model=None,
        copy=False,
        random_seed=settings["random_seed"],
    )
    if result is not None:
        raise RuntimeError("Schist fit_model(copy=False) unexpectedly returned an AnnData object.")
    if _standalone_graph_fingerprint(adjacency, adata.obs_names, np=np) != adjacency_fingerprint:
        raise RuntimeError("Schist backend mutated the validated adjacency matrix.")
    if (
        _standalone_state_fingerprint(
            adata,
            key_added=settings["key_added"],
            np=np,
            pd=pd,
            sparse=sparse,
        )
        != source_fingerprint
    ):
        raise RuntimeError("Schist execution mutated state outside its owned output family.")
    hierarchy = _standalone_validate_backend_output(
        output,
        key_added=settings["key_added"],
        graph=graph,
        settings=settings,
        source_fingerprint=source_fingerprint,
        np=np,
        pd=pd,
        sparse=sparse,
    )
    _standalone_h5ad_roundtrip(
        output,
        key_added=settings["key_added"],
        expected_obs_keys=hierarchy["obs_keys"],
        expected_obsm_keys=hierarchy["obsm_keys"],
    )
    if (
        _standalone_state_fingerprint(
            output,
            key_added=settings["key_added"],
            np=np,
            pd=pd,
            sparse=sparse,
        )
        != source_fingerprint
    ):
        raise RuntimeError("H5AD round-trip validation mutated input state outside the Schist-owned output family.")
    summary = _standalone_summary(
        adata=adata,
        settings=settings,
        graph=graph,
        memory=memory,
        backend=backend,
        hierarchy=hierarchy,
        overwritten=bool(family["collisions"]),
        runtime_duration_seconds=time.perf_counter() - started_at,
        openbio_version=openbio_version,
    )
    return output, summary


_STANDALONE_HELPERS = (
    _standalone_plain_json,
    _standalone_validate_settings,
    _standalone_output_family,
    _standalone_object_bytes,
    _standalone_adata_bytes,
    _standalone_memory_estimate,
    _standalone_numeric_summary,
    _standalone_graph_fingerprint,
    _standalone_validate_graph,
    _standalone_require_backend,
    _standalone_hash_update,
    _standalone_state_fingerprint,
    _standalone_clear_output_family,
    _standalone_validate_backend_output,
    _standalone_h5ad_roundtrip,
    _standalone_backend_build,
    _standalone_software_versions,
    _standalone_summary,
    _standalone_schist_nested_model,
)


def run_schist_nested_model(
    adata: Any,
    *,
    random_seed: int = 123,
    neighbors_key: str = "neighbors",
    key_added: str = "nsbm",
    posterior_samples: int = 100,
    degree_correction: bool = True,
    overwrite_existing: bool = False,
    max_working_memory_gib: float = 8.0,
    openbio_version: str = "unknown",
) -> tuple[Any, dict[str, Any]]:
    return _standalone_schist_nested_model(
        adata,
        random_seed=random_seed,
        neighbors_key=neighbors_key,
        key_added=key_added,
        posterior_samples=posterior_samples,
        degree_correction=degree_correction,
        overwrite_existing=overwrite_existing,
        max_working_memory_gib=max_working_memory_gib,
        openbio_version=openbio_version,
    )


def schist_nested_model_code(
    *,
    random_seed: int,
    neighbors_key: str,
    key_added: str,
    posterior_samples: int,
    degree_correction: bool,
    overwrite_existing: bool,
    max_working_memory_gib: float,
    openbio_version: str,
) -> str:
    sources = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in _STANDALONE_HELPERS)
    return f"""from __future__ import annotations

{sources}


def run_schist_nsbm(adata):
    return _standalone_schist_nested_model(
        adata,
        random_seed={random_seed!r},
        neighbors_key={neighbors_key!r},
        key_added={key_added!r},
        posterior_samples={posterior_samples!r},
        degree_correction={degree_correction!r},
        overwrite_existing={overwrite_existing!r},
        max_working_memory_gib={max_working_memory_gib!r},
        openbio_version={openbio_version!r},
    )
"""


__all__ = [
    "run_schist_nested_model",
    "schist_nested_model_code",
]
