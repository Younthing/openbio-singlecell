from __future__ import annotations

import copy
import hashlib
import importlib
import inspect
import io
import json
import math
import platform
import random
import re
import sys
import threading
import warnings
from collections.abc import Mapping, Sequence
from contextlib import contextmanager
from importlib import metadata as importlib_metadata
from typing import Any

VELOCITY_PORTABLE_SCHEMA = "openbio-singlecell/velocity-portable-state/v1"
VELOCITY_SUMMARY_SCHEMA = "openbio-singlecell/velocity-summary/v1"
VELOCITY_ARTIFACT_TYPE = "OPENBIO_VELOCITY_STATE"
VELOCITY_STATE_KEY = "openbio_velocity_state"
VELOCITY_SCVELO_VERSION = "0.3.4"
VELOCITY_STAGES = (
    "prepared",
    "moments",
    "velocity_estimated",
    "dynamics_recovered",
    "velocity_graph",
)

VELOCITY_REFERENCES = [
    {
        "citation": (
            "Bergen V, Lange M, Peidli S, et al. Generalizing RNA velocity to transient cell states "
            "through dynamical modeling. Nature Biotechnology. 2020;38:1408-1414."
        ),
        "doi": "10.1038/s41587-020-0591-3",
        "url": "https://doi.org/10.1038/s41587-020-0591-3",
        "kind": "method",
    },
    {
        "citation": (
            "La Manno G, Soldatov R, Zeisel A, et al. RNA velocity of single cells. "
            "Nature. 2018;560:494-498."
        ),
        "doi": "10.1038/s41586-018-0414-6",
        "url": "https://doi.org/10.1038/s41586-018-0414-6",
        "kind": "method",
    },
    {
        "citation": "scVelo 0.3.4 public interface and signed release source.",
        "doi": None,
        "url": "https://github.com/theislab/scvelo/releases/tag/v0.3.4",
        "kind": "software_documentation",
    },
]

MOMENT_REFERENCES = [
    *VELOCITY_REFERENCES,
    {
        "citation": (
            "Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ. Diffusion pseudotime robustly "
            "reconstructs lineage branching. Nature Methods. 2016;13:845-848."
        ),
        "doi": "10.1038/nmeth.3971",
        "url": "https://doi.org/10.1038/nmeth.3971",
        "kind": "method",
    },
]

DYNAMICS_REFERENCES = [
    *VELOCITY_REFERENCES,
    {
        "citation": (
            "Bergen V, Soldatov RA, Kharchenko PV, Theis FJ. RNA velocity-current challenges and future "
            "perspectives. Molecular Systems Biology. 2021;17:e10282."
        ),
        "doi": "10.15252/msb.202110282",
        "url": "https://doi.org/10.15252/msb.202110282",
        "kind": "practice",
    },
]

_STATE_METADATA_KEYS = {
    "schema_version",
    "artifact_type",
    "stage",
    "producer_node_id",
    "scvelo_version",
    "parent_state_fingerprint_sha256",
    "parameters",
    "parameters_fingerprint_sha256",
    "neighbors_key",
    "neighbor_graph_fingerprint_sha256",
    "vkey",
    "velocity_mode",
    "dynamics_fit_fingerprint_sha256",
    "velocity_graph_positive_fingerprint_sha256",
    "velocity_graph_negative_fingerprint_sha256",
    "state_fingerprint_sha256",
}
_VELOCITY_BACKEND_LOCK = threading.RLock()

_EXPECTED_SIGNATURES = {
    "normalize_per_cell": (
        "data",
        "counts_per_cell_after",
        "counts_per_cell",
        "key_n_counts",
        "max_proportion_per_cell",
        "use_initial_size",
        "layers",
        "enforce",
        "copy",
    ),
    "moments": (
        "data",
        "n_neighbors",
        "n_pcs",
        "mode",
        "method",
        "use_rep",
        "use_highly_variable",
        "copy",
    ),
    "velocity": (
        "data",
        "vkey",
        "mode",
        "fit_offset",
        "fit_offset2",
        "filter_genes",
        "groups",
        "groupby",
        "groups_for_fit",
        "constrain_ratio",
        "use_raw",
        "use_latent_time",
        "perc",
        "min_r2",
        "min_likelihood",
        "r2_adjusted",
        "use_highly_variable",
        "diff_kinetics",
        "copy",
        "kwargs",
    ),
    "velocity_graph": (
        "data",
        "vkey",
        "xkey",
        "tkey",
        "basis",
        "n_neighbors",
        "n_recurse_neighbors",
        "random_neighbors_at_max",
        "sqrt_transform",
        "variance_stabilization",
        "gene_subset",
        "compute_uncertainties",
        "approx",
        "mode_neighbors",
        "copy",
        "n_jobs",
        "backend",
        "show_progress_bar",
    ),
    "recover_dynamics": (
        "data",
        "var_names",
        "n_top_genes",
        "max_iter",
        "assignment_mode",
        "t_max",
        "fit_time",
        "fit_scaling",
        "fit_steady_states",
        "fit_connected_states",
        "fit_basal_transcription",
        "use_raw",
        "load_pars",
        "return_model",
        "plot_results",
        "steady_state_prior",
        "add_key",
        "copy",
        "n_jobs",
        "backend",
        "show_progress_bar",
        "kwargs",
    ),
    "velocity_embedding_stream": (
        "adata",
        "basis",
        "vkey",
        "density",
        "smooth",
        "min_mass",
        "cutoff_perc",
        "arrow_color",
        "arrow_size",
        "arrow_style",
        "max_length",
        "integration_direction",
        "linewidth",
        "n_neighbors",
        "recompute",
        "color",
        "use_raw",
        "layer",
        "color_map",
        "colorbar",
        "palette",
        "size",
        "alpha",
        "perc",
        "X",
        "V",
        "X_grid",
        "V_grid",
        "sort_order",
        "groups",
        "components",
        "legend_loc",
        "legend_fontsize",
        "legend_fontweight",
        "xlabel",
        "ylabel",
        "title",
        "fontsize",
        "figsize",
        "dpi",
        "frameon",
        "show",
        "save",
        "ax",
        "ncols",
        "kwargs",
    ),
}


def _vs_science() -> tuple[Any, Any, Any, Any]:
    import anndata
    import numpy
    import pandas
    from scipy import sparse

    return anndata, numpy, pandas, sparse


def _vs_clean_text(value: Any, *, label: str, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise TypeError(f"Velocity {label} must be a string.")
    if value != value.strip():
        raise ValueError(f"Velocity {label} must not contain surrounding whitespace.")
    if not value and not allow_empty:
        raise ValueError(f"Velocity {label} cannot be empty.")
    return value


def _vs_output_key(value: Any, *, label: str) -> str:
    value = _vs_clean_text(value, label=label)
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", value) is None:
        raise ValueError(
            f"Velocity {label} must start with a letter and contain only letters, digits, and underscores."
        )
    if value in {"spliced", "unspliced", "Ms", "Mu"} or value.startswith("fit_"):
        raise ValueError(f"Velocity {label} collides with a reserved staged layer/fit namespace.")
    return value


def _vs_velocity_var_keys(vkey: str) -> tuple[str, ...]:
    """Return every var key written by scVelo 0.3.4's steady-state estimator."""
    return tuple(
        f"{vkey}{suffix}"
        for suffix in (
            "_offset",
            "_offset2",
            "_beta",
            "_gamma",
            "_qreg_ratio",
            "_r2",
            "_genes",
        )
    )


def _vs_integer(value: Any, *, label: str, minimum: int, maximum: int = 2**31 - 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"Velocity {label} must be an integer.")
    if not minimum <= value <= maximum:
        raise ValueError(f"Velocity {label} must be between {minimum} and {maximum}.")
    return value


def _vs_float(value: Any, *, label: str, minimum: float, strictly_greater: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"Velocity {label} must be numeric.")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"Velocity {label} must be finite.")
    if (strictly_greater and value <= minimum) or (not strictly_greater and value < minimum):
        relation = ">" if strictly_greater else ">="
        raise ValueError(f"Velocity {label} must be {relation} {minimum}.")
    return value


def _vs_boolean(value: Any, *, label: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"Velocity {label} must be a boolean.")
    return value


def _vs_json_value(value: Any) -> Any:
    if value is None or isinstance(value, str | bool | int):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "item"):
        try:
            return _vs_json_value(value.item())
        except (TypeError, ValueError):
            pass
    if hasattr(value, "tolist"):
        try:
            return _vs_json_value(value.tolist())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _vs_json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return [_vs_json_value(item) for item in value]
    return str(value)


def _vs_json_sha256(value: Any) -> str:
    encoded = json.dumps(
        _vs_json_value(value), ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _vs_package_version(package: str) -> str:
    try:
        return importlib_metadata.version(package)
    except importlib_metadata.PackageNotFoundError:
        return "not-installed"


def _vs_versions(scvelo_version: str) -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "openbio-singlecell": _vs_package_version("openbio-singlecell"),
        "scvelo": scvelo_version,
        "anndata": _vs_package_version("anndata"),
        "numpy": _vs_package_version("numpy"),
        "pandas": _vs_package_version("pandas"),
        "scipy": _vs_package_version("scipy"),
        "scanpy": _vs_package_version("scanpy"),
        "matplotlib": _vs_package_version("matplotlib"),
    }


def _vs_backend(scvelo_module: Any | None, operation: str) -> tuple[Any, str, Any]:
    if scvelo_module is None:
        try:
            import scvelo as scvelo_module
        except (ImportError, OSError) as exc:
            raise RuntimeError("Velocity analysis requires exact scVelo 0.3.4.") from exc
    version = getattr(scvelo_module, "__version__", None)
    if not isinstance(version, str) or not version:
        try:
            version = importlib_metadata.version("scvelo")
        except importlib_metadata.PackageNotFoundError as exc:
            raise RuntimeError("Velocity analysis could not establish the scVelo distribution version.") from exc
    if version != VELOCITY_SCVELO_VERSION:
        raise RuntimeError(
            f"Velocity analysis requires exact scVelo {VELOCITY_SCVELO_VERSION}; detected {version!r}."
        )
    if getattr(scvelo_module, "__name__", None) == "scvelo":
        try:
            distribution_version = importlib_metadata.version("scvelo")
        except importlib_metadata.PackageNotFoundError as exc:
            raise RuntimeError("Velocity analysis could not verify the scVelo distribution identity.") from exc
        if distribution_version != VELOCITY_SCVELO_VERSION or distribution_version != version:
            raise RuntimeError(
                "Velocity analysis requires matching module/distribution scVelo 0.3.4 identity; "
                f"module={version!r}, distribution={distribution_version!r}."
            )
    location = {
        "normalize_per_cell": ("pp", "normalize_per_cell"),
        "moments": ("pp", "moments"),
        "velocity": ("tl", "velocity"),
        "velocity_graph": ("tl", "velocity_graph"),
        "recover_dynamics": ("tl", "recover_dynamics"),
        "velocity_embedding_stream": ("pl", "velocity_embedding_stream"),
    }[operation]
    namespace = getattr(scvelo_module, location[0], None)
    function = getattr(namespace, location[1], None)
    if not callable(function):
        raise RuntimeError(f"scVelo 0.3.4 {location[0]}.{location[1]} is unavailable.")
    try:
        observed = tuple(inspect.signature(function).parameters)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"Could not inspect scVelo 0.3.4 {location[0]}.{location[1]}.") from exc
    expected = _EXPECTED_SIGNATURES[operation]
    if observed != expected:
        raise RuntimeError(
            f"scVelo 0.3.4 {location[0]}.{location[1]} interface is incompatible; "
            f"expected={expected}, observed={observed}."
        )
    return scvelo_module, version, function


@contextmanager
def _vs_preserve_global_state_unlocked(scvelo_module: Any):
    _, numpy, _, _ = _vs_science()
    numpy_state = numpy.random.get_state()
    python_state = random.getstate()
    print_options = numpy.get_printoptions()
    matplotlib_state = None
    if "matplotlib" in sys.modules:
        import matplotlib

        matplotlib_state = matplotlib.rcParams.copy()
    settings = getattr(scvelo_module, "settings", None)
    setting_names = (
        "verbosity",
        "presenter_view",
        "autoshow",
        "autosave",
        "figdir",
        "plot_prefix",
        "dpi",
        "dpi_save",
        "frameon",
        "vector_friendly",
        "file_format_figs",
        "transparent",
        "color_map",
    )
    setting_values = {
        name: copy.deepcopy(getattr(settings, name))
        for name in setting_names
        if settings is not None and hasattr(settings, name)
    }
    try:
        yield
    finally:
        numpy.random.set_state(numpy_state)
        random.setstate(python_state)
        numpy.set_printoptions(**print_options)
        if matplotlib_state is not None:
            import matplotlib

            matplotlib.rcParams.update(matplotlib_state)
        for name, value in setting_values.items():
            setattr(settings, name, value)


@contextmanager
def _vs_preserve_global_state(scvelo_module: Any):
    with _VELOCITY_BACKEND_LOCK, _vs_preserve_global_state_unlocked(scvelo_module):
        yield


@contextmanager
def _vs_recover_pandas3_compatibility(recover: Any):
    """Scope the two audited Pandas 3 compatibility fixes to one backend call."""

    _, numpy, pandas, _ = _vs_science()
    function_globals = getattr(recover, "__globals__", None)
    if not isinstance(function_globals, dict):
        yield ()
        return
    original_unique = function_globals.get("make_unique_list")
    original_read = function_globals.get("_read_pars")
    if original_unique is None and original_read is None:
        yield ()
        return
    if original_unique is None or original_read is None:
        raise RuntimeError("scVelo 0.3.4 Pandas compatibility seams are incomplete.")
    try:
        observed_unique = tuple(inspect.signature(original_unique).parameters)
        observed_read = tuple(inspect.signature(original_read).parameters)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("Could not inspect scVelo 0.3.4 Pandas compatibility seams.") from exc
    if observed_unique != ("key", "allow_array"):
        raise RuntimeError(
            "scVelo 0.3.4 make_unique_list compatibility seam changed; "
            f"expected=('key', 'allow_array'), observed={observed_unique}."
        )
    if observed_read != ("adata", "pars_names", "key"):
        raise RuntimeError(
            "scVelo 0.3.4 _read_pars compatibility seam changed; "
            f"expected=('adata', 'pars_names', 'key'), observed={observed_read}."
        )

    def compatible_make_unique_list(key: Any, allow_array: bool = False):
        if isinstance(key, list | tuple):
            values = numpy.asarray(key, dtype=object)
            if values.ndim == 1 and all(isinstance(value, str) for value in values):
                return pandas.unique(values)
        return original_unique(key, allow_array=allow_array)

    def compatible_read_pars(adata: Any, pars_names: Any = None, key: str = "fit"):
        return [
            numpy.array(values, copy=True)
            for values in original_read(adata, pars_names=pars_names, key=key)
        ]

    with _VELOCITY_BACKEND_LOCK:
        if (
            function_globals.get("make_unique_list") is not original_unique
            or function_globals.get("_read_pars") is not original_read
        ):
            raise RuntimeError("scVelo Pandas compatibility seams changed before the guarded dynamics call.")
        function_globals["make_unique_list"] = compatible_make_unique_list
        function_globals["_read_pars"] = compatible_read_pars
        try:
            yield ("stable_unique", "writable_parameter_arrays")
        finally:
            function_globals["make_unique_list"] = original_unique
            function_globals["_read_pars"] = original_read


@contextmanager
def _vs_velocity_pandas3_compatibility(scvelo_module: Any, *, enabled: bool):
    """Make only read-only divergence ndarrays writable for audited dynamical velocity."""

    if not enabled or getattr(scvelo_module, "__name__", None) != "scvelo":
        yield ()
        return
    try:
        utilities = importlib.import_module("scvelo.tools._em_model_utils")
    except (ImportError, OSError) as exc:
        raise RuntimeError("Could not import the audited scVelo dynamical utility module.") from exc
    original = getattr(utilities, "compute_divergence", None)
    if not callable(original):
        raise RuntimeError("scVelo 0.3.4 compute_divergence compatibility seam is unavailable.")
    expected = (
        "u",
        "s",
        "alpha",
        "beta",
        "gamma",
        "scaling",
        "t_",
        "u0_",
        "s0_",
        "tau",
        "tau_",
        "std_u",
        "std_s",
        "normalized",
        "mode",
        "assignment_mode",
        "var_scale",
        "kernel_width",
        "fit_steady_states",
        "connectivities",
        "constraint_time_increments",
        "reg_time",
        "reg_par",
        "min_confidence",
        "pval_steady",
        "steady_u",
        "steady_s",
        "noise_model",
        "time_connectivities",
        "clusters",
        "kwargs",
    )
    try:
        observed = tuple(inspect.signature(original).parameters)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("Could not inspect scVelo 0.3.4 compute_divergence.") from exc
    if observed != expected:
        raise RuntimeError(
            "scVelo 0.3.4 compute_divergence compatibility seam changed; "
            f"expected={expected}, observed={observed}."
        )

    def compatible_compute_divergence(*args: Any, **kwargs: Any):
        _, numpy, _, _ = _vs_science()

        def writable(value: Any) -> Any:
            if isinstance(value, numpy.ndarray) and not bool(value.flags.writeable):
                return numpy.array(value, copy=True)
            return value

        return original(
            *(writable(value) for value in args),
            **{key: writable(value) for key, value in kwargs.items()},
        )

    with _VELOCITY_BACKEND_LOCK:
        if getattr(utilities, "compute_divergence", None) is not original:
            raise RuntimeError("scVelo compute_divergence changed before the guarded velocity call.")
        utilities.compute_divergence = compatible_compute_divergence
        try:
            yield ("writable_divergence_inputs",)
        finally:
            utilities.compute_divergence = original


def _vs_axis(index: Any, *, label: str) -> tuple[str, ...]:
    if not bool(index.is_unique):
        raise ValueError(f"Velocity requires unique {label} identifiers.")
    values = []
    for position, value in enumerate(index.tolist()):
        if not isinstance(value, str) or not value or value != value.strip():
            raise TypeError(
                f"Velocity {label} identifier at position {position} must be a canonical nonempty string."
            )
        values.append(value)
    return tuple(values)


def _vs_matrix_values(matrix: Any, *, label: str, shape: tuple[int, int], finite: bool = True) -> Any:
    _, numpy, _, sparse = _vs_science()
    if tuple(getattr(matrix, "shape", ())) != shape:
        raise ValueError(f"Velocity {label} must have shape {shape}.")
    values = numpy.asarray(matrix.data) if sparse.issparse(matrix) else numpy.asarray(matrix).ravel()
    if numpy.issubdtype(values.dtype, numpy.bool_) or not numpy.issubdtype(values.dtype, numpy.number):
        raise TypeError(f"Velocity {label} must contain real numeric values.")
    if bool(numpy.iscomplexobj(values)):
        raise TypeError(f"Velocity {label} must contain real numeric values.")
    numeric = values.astype(float, copy=False)
    if finite and numeric.size and not bool(numpy.isfinite(numeric).all()):
        raise ValueError(f"Velocity {label} contains non-finite values.")
    return numeric


def _vs_matrix_fingerprint(matrix: Any, *, label: str) -> str:
    _, numpy, _, sparse = _vs_science()
    digest = hashlib.sha256()
    shape = tuple(int(value) for value in matrix.shape)
    digest.update(json.dumps({"label": label, "shape": shape}, sort_keys=True).encode("utf-8"))
    if sparse.issparse(matrix):
        canonical = matrix.tocsr(copy=True)
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
        digest.update(b"sparse-csr\0")
        digest.update(str(canonical.dtype).encode("ascii"))
        for array in (canonical.indptr, canonical.indices):
            digest.update(numpy.ascontiguousarray(array, dtype="<i8").tobytes(order="C"))
        data = numpy.asarray(canonical.data)
        digest.update(numpy.ascontiguousarray(data).tobytes(order="C"))
    else:
        array = numpy.asarray(matrix)
        digest.update(b"dense\0")
        digest.update(str(array.dtype).encode("ascii"))
        for start in range(0, array.shape[0], 256):
            digest.update(numpy.ascontiguousarray(array[start : start + 256]).tobytes(order="C"))
    return digest.hexdigest()


def _vs_series_fingerprint(frame: Any, columns: Sequence[str], *, label: str) -> str:
    _, _, pandas, _ = _vs_science()
    selected = frame.loc[:, list(columns)].copy()
    digest = hashlib.sha256()
    digest.update(label.encode("utf-8"))
    digest.update(json.dumps([str(value) for value in selected.index], ensure_ascii=False).encode("utf-8"))
    for column in selected.columns:
        digest.update(str(column).encode("utf-8"))
        digest.update(str(selected[column].dtype).encode("utf-8"))
        hashed = pandas.util.hash_pandas_object(selected[column], index=False, categorize=True)
        digest.update(hashed.to_numpy(dtype="<u8", copy=False).tobytes(order="C"))
    return digest.hexdigest()


def _vs_numeric_summary(values: Any) -> dict[str, int | float | None]:
    _, numpy, _, _ = _vs_science()
    array = numpy.asarray(values, dtype=float).ravel()
    finite = array[numpy.isfinite(array)]
    if not finite.size:
        return {
            "n": int(array.size),
            "missing": int(array.size),
            "minimum": None,
            "q25": None,
            "median": None,
            "mean": None,
            "q75": None,
            "maximum": None,
        }
    quantiles = numpy.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "n": int(array.size),
        "missing": int(array.size - finite.size),
        "minimum": float(quantiles[0]),
        "q25": float(quantiles[1]),
        "median": float(quantiles[2]),
        "mean": float(finite.mean()),
        "q75": float(quantiles[3]),
        "maximum": float(quantiles[4]),
    }


def _vs_row_sums(matrix: Any) -> Any:
    _, numpy, _, sparse = _vs_science()
    if sparse.issparse(matrix):
        return numpy.asarray(matrix.sum(axis=1), dtype=float).ravel()
    return numpy.asarray(matrix, dtype=float).sum(axis=1)


def _vs_graph_matrix(matrix: Any, *, n_obs: int, label: str, symmetric: bool) -> Any:
    _, numpy, _, sparse = _vs_science()
    if not sparse.issparse(matrix):
        raise TypeError(f"Velocity {label} graph matrix must be SciPy sparse.")
    if tuple(matrix.shape) != (n_obs, n_obs):
        raise ValueError(f"Velocity {label} graph matrix must have shape ({n_obs}, {n_obs}).")
    matrix = matrix.tocsr(copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    values = _vs_matrix_values(matrix, label=f"{label} graph", shape=(n_obs, n_obs))
    if values.size and bool((values < 0).any()):
        raise ValueError(f"Velocity {label} graph contains negative weights.")
    if not bool(numpy.allclose(matrix.diagonal(), 0.0, rtol=0.0, atol=1e-12)):
        raise ValueError(f"Velocity {label} graph must have a zero diagonal.")
    if symmetric:
        delta = matrix - matrix.T
        delta.eliminate_zeros()
        if delta.nnz and float(numpy.max(numpy.abs(delta.data))) > 1e-8:
            raise ValueError(f"Velocity {label} graph must be symmetric.")
    if matrix.nnz == 0:
        raise ValueError(f"Velocity {label} graph contains no edges.")
    return matrix


def _vs_resolve_graph(adata: Any, neighbors_key: str) -> dict[str, Any]:
    neighbors_key = _vs_clean_text(neighbors_key, label="neighbors_key")
    record = adata.uns.get(neighbors_key)
    if not isinstance(record, Mapping):
        raise ValueError(f"Velocity neighbor metadata not found in uns[{neighbors_key!r}].")
    connectivities_key = _vs_clean_text(record.get("connectivities_key"), label="connectivities_key")
    distances_key = _vs_clean_text(record.get("distances_key"), label="distances_key")
    if connectivities_key not in adata.obsp or distances_key not in adata.obsp:
        raise ValueError("Velocity named graph matrices are missing from obsp.")
    connectivities = _vs_graph_matrix(
        adata.obsp[connectivities_key], n_obs=int(adata.n_obs), label="connectivity", symmetric=True
    )
    distances = _vs_graph_matrix(
        adata.obsp[distances_key], n_obs=int(adata.n_obs), label="distance", symmetric=False
    )
    fingerprint = _vs_json_sha256(
        {
            "schema": "openbio-singlecell/velocity-neighbor-graph/v1",
            "neighbors_key": neighbors_key,
            "connectivities_key": connectivities_key,
            "distances_key": distances_key,
            "obs_names": list(_vs_axis(adata.obs_names, label="observation")),
            "connectivities_fingerprint_sha256": _vs_matrix_fingerprint(
                connectivities, label="connectivities"
            ),
            "distances_fingerprint_sha256": _vs_matrix_fingerprint(distances, label="distances"),
            "params": _vs_json_value(record.get("params", {})),
        }
    )
    return {
        "neighbors_key": neighbors_key,
        "connectivities_key": connectivities_key,
        "distances_key": distances_key,
        "connectivities": connectivities,
        "distances": distances,
        "fingerprint_sha256": fingerprint,
        "params": _vs_json_value(record.get("params", {})),
    }


@contextmanager
def _vs_alias_graph(adata: Any, graph: Mapping[str, Any]):
    had_uns = "neighbors" in adata.uns
    had_connectivities = "connectivities" in adata.obsp
    had_distances = "distances" in adata.obsp
    previous_uns = copy.deepcopy(adata.uns["neighbors"]) if had_uns else None
    previous_connectivities = adata.obsp["connectivities"] if had_connectivities else None
    previous_distances = adata.obsp["distances"] if had_distances else None
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": copy.deepcopy(graph["params"]),
    }
    adata.obsp["connectivities"] = graph["connectivities"].copy()
    adata.obsp["distances"] = graph["distances"].copy()
    try:
        yield
    finally:
        if not had_uns:
            del adata.uns["neighbors"]
        else:
            adata.uns["neighbors"] = previous_uns
        if not had_connectivities:
            del adata.obsp["connectivities"]
        else:
            adata.obsp["connectivities"] = previous_connectivities
        if not had_distances:
            del adata.obsp["distances"]
        else:
            adata.obsp["distances"] = previous_distances


@contextmanager
def _vs_optional_alias_graph(adata: Any, graph: Mapping[str, Any] | None):
    if graph is None:
        yield
    else:
        with _vs_alias_graph(adata, graph):
            yield


def _vs_dynamics_components(adata: Any) -> dict[str, Any]:
    var_columns = sorted(column for column in adata.var.columns if str(column).startswith("fit_"))
    layer_keys = sorted(key for key in adata.layers.keys() if str(key).startswith("fit_"))
    components = {
        "var_columns": var_columns,
        "var_fingerprint_sha256": (
            _vs_series_fingerprint(adata.var, var_columns, label="dynamics-var") if var_columns else None
        ),
        "layers": {
            key: _vs_matrix_fingerprint(adata.layers[key], label=f"dynamics-layer:{key}")
            for key in layer_keys
        },
        "loss_fingerprint_sha256": (
            _vs_matrix_fingerprint(adata.varm["loss"], label="dynamics-loss")
            if "loss" in adata.varm
            else None
        ),
        "recover_parameters": _vs_json_value(adata.uns.get("recover_dynamics")),
        "selected_fingerprint_sha256": (
            _vs_series_fingerprint(
                adata.var,
                [
                    "openbio_dynamics_selected",
                    "openbio_dynamics_fit_success",
                    "openbio_dynamics_velocity_gene",
                ],
                label="dynamics-selected",
            )
            if {
                "openbio_dynamics_selected",
                "openbio_dynamics_fit_success",
                "openbio_dynamics_velocity_gene",
            }.issubset(adata.var.columns)
            else None
        ),
    }
    components["fingerprint_sha256"] = _vs_json_sha256(components)
    return components


def _vs_state_components(adata: Any, state: Mapping[str, Any]) -> dict[str, Any]:
    stage = state["stage"]
    components: dict[str, Any] = {
        "obs_fingerprint_sha256": _vs_json_sha256(list(_vs_axis(adata.obs_names, label="observation"))),
        "var_fingerprint_sha256": _vs_json_sha256(list(_vs_axis(adata.var_names, label="feature"))),
        "X_fingerprint_sha256": _vs_matrix_fingerprint(adata.X, label="X"),
        "spliced_fingerprint_sha256": _vs_matrix_fingerprint(
            adata.layers["spliced"], label="spliced"
        ),
        "unspliced_fingerprint_sha256": _vs_matrix_fingerprint(
            adata.layers["unspliced"], label="unspliced"
        ),
    }
    if stage != "prepared":
        components["Ms_fingerprint_sha256"] = _vs_matrix_fingerprint(adata.layers["Ms"], label="Ms")
        components["Mu_fingerprint_sha256"] = _vs_matrix_fingerprint(adata.layers["Mu"], label="Mu")
    neighbors_key = state.get("neighbors_key")
    if neighbors_key is not None:
        components["neighbor_graph_fingerprint_sha256"] = _vs_resolve_graph(
            adata, str(neighbors_key)
        )["fingerprint_sha256"]
    vkey = state.get("vkey")
    if vkey is not None:
        components["velocity_layer_fingerprint_sha256"] = _vs_matrix_fingerprint(
            adata.layers[vkey], label=f"velocity:{vkey}"
        )
        velocity_var_columns = [
            key
            for key in _vs_velocity_var_keys(str(vkey))
            if key in adata.var
        ]
        components["velocity_var_fingerprint_sha256"] = _vs_series_fingerprint(
            adata.var, velocity_var_columns, label=f"velocity-var:{vkey}"
        )
        components["velocity_auxiliary_layer_fingerprints_sha256"] = {
            key: _vs_matrix_fingerprint(adata.layers[key], label=f"velocity-auxiliary:{key}")
            for key in (f"variance_{vkey}", f"{vkey}_u")
            if key in adata.layers
        }
        components["velocity_parameters_fingerprint_sha256"] = _vs_json_sha256(
            adata.uns.get(f"{vkey}_params")
        )
    if state.get("dynamics_fit_fingerprint_sha256") is not None:
        components["dynamics_fit_fingerprint_sha256"] = _vs_dynamics_components(adata)[
            "fingerprint_sha256"
        ]
    if stage == "velocity_graph":
        positive_graph = _vs_velocity_graph_matrix(
            adata.uns[f"{vkey}_graph"], n_obs=int(adata.n_obs), label="positive", sign="positive"
        )
        negative_graph = _vs_velocity_graph_matrix(
            adata.uns[f"{vkey}_graph_neg"], n_obs=int(adata.n_obs), label="negative", sign="negative"
        )
        components["velocity_graph_positive_fingerprint_sha256"] = _vs_matrix_fingerprint(
            positive_graph, label=f"{vkey}-graph-positive"
        )
        components["velocity_graph_negative_fingerprint_sha256"] = _vs_matrix_fingerprint(
            negative_graph, label=f"{vkey}-graph-negative"
        )
        components["self_transition_fingerprint_sha256"] = _vs_series_fingerprint(
            adata.obs, [f"{vkey}_self_transition"], label=f"{vkey}-self-transition"
        )
    return components


def _vs_write_state(
    adata: Any,
    *,
    stage: str,
    producer_node_id: str,
    scvelo_version: str,
    parameters: Mapping[str, Any],
    parent: Mapping[str, Any] | None,
    neighbors_key: str | None,
    vkey: str | None,
    velocity_mode: str | None,
    has_dynamics: bool,
) -> dict[str, Any]:
    if stage not in VELOCITY_STAGES:
        raise ValueError(f"Unsupported velocity stage: {stage!r}.")
    normalized_parameters = _vs_json_value(parameters)
    state = {
        "schema_version": VELOCITY_PORTABLE_SCHEMA,
        "artifact_type": VELOCITY_ARTIFACT_TYPE,
        "stage": stage,
        "producer_node_id": producer_node_id,
        "scvelo_version": scvelo_version,
        "parent_state_fingerprint_sha256": (
            None if parent is None else parent["state_fingerprint_sha256"]
        ),
        "parameters": normalized_parameters,
        "parameters_fingerprint_sha256": _vs_json_sha256(normalized_parameters),
        "neighbors_key": neighbors_key,
        "neighbor_graph_fingerprint_sha256": None,
        "vkey": vkey,
        "velocity_mode": velocity_mode,
        "dynamics_fit_fingerprint_sha256": None,
        "velocity_graph_positive_fingerprint_sha256": None,
        "velocity_graph_negative_fingerprint_sha256": None,
        "state_fingerprint_sha256": "",
    }
    components = _vs_state_components(adata, state)
    state["neighbor_graph_fingerprint_sha256"] = components.get(
        "neighbor_graph_fingerprint_sha256"
    )
    if has_dynamics:
        dynamics = _vs_dynamics_components(adata)["fingerprint_sha256"]
        state["dynamics_fit_fingerprint_sha256"] = dynamics
        components["dynamics_fit_fingerprint_sha256"] = dynamics
    if stage == "velocity_graph":
        state["velocity_graph_positive_fingerprint_sha256"] = components[
            "velocity_graph_positive_fingerprint_sha256"
        ]
        state["velocity_graph_negative_fingerprint_sha256"] = components[
            "velocity_graph_negative_fingerprint_sha256"
        ]
    state["state_fingerprint_sha256"] = _vs_json_sha256(
        {"metadata": {key: value for key, value in state.items() if key != "state_fingerprint_sha256"},
         "components": components}
    )
    adata.uns[VELOCITY_STATE_KEY] = copy.deepcopy(state)
    return state


def validate_portable_velocity_state(
    adata: Any,
    *,
    allowed_stages: Sequence[str] | None = None,
) -> dict[str, Any]:
    anndata, _, _, _ = _vs_science()
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("Velocity portable state must be an AnnData object.")
    if bool(adata.isbacked):
        raise ValueError("Velocity state must be in memory; call to_memory() first.")
    state = adata.uns.get(VELOCITY_STATE_KEY)
    if not isinstance(state, Mapping) or set(state) != _STATE_METADATA_KEYS:
        raise TypeError("Velocity input lacks the exact staged portable-state metadata.")
    state = copy.deepcopy(dict(state))
    if state["schema_version"] != VELOCITY_PORTABLE_SCHEMA or state["artifact_type"] != VELOCITY_ARTIFACT_TYPE:
        raise ValueError("Velocity portable-state identity is invalid.")
    stage = state["stage"]
    if stage not in VELOCITY_STAGES:
        raise ValueError("Velocity portable state has an unsupported stage.")
    if allowed_stages is not None and stage not in set(allowed_stages):
        raise ValueError(
            f"Velocity state stage {stage!r} is not valid here; expected one of {list(allowed_stages)}."
        )
    if state["scvelo_version"] != VELOCITY_SCVELO_VERSION:
        raise ValueError("Velocity state was not produced by exact scVelo 0.3.4.")
    if state["parameters_fingerprint_sha256"] != _vs_json_sha256(state["parameters"]):
        raise ValueError("Velocity state parameter provenance fingerprint is invalid.")
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise ValueError("Velocity state lacks canonical spliced/unspliced layers.")
    if stage != "prepared" and ("Ms" not in adata.layers or "Mu" not in adata.layers):
        raise ValueError("Velocity state lacks canonical moment layers.")
    components = _vs_state_components(adata, state)
    for metadata_key, component_key in (
        ("neighbor_graph_fingerprint_sha256", "neighbor_graph_fingerprint_sha256"),
        ("dynamics_fit_fingerprint_sha256", "dynamics_fit_fingerprint_sha256"),
        ("velocity_graph_positive_fingerprint_sha256", "velocity_graph_positive_fingerprint_sha256"),
        ("velocity_graph_negative_fingerprint_sha256", "velocity_graph_negative_fingerprint_sha256"),
    ):
        if state[metadata_key] != components.get(component_key):
            raise ValueError(f"Velocity state {metadata_key} failed its current-content check.")
    expected = _vs_json_sha256(
        {"metadata": {key: value for key, value in state.items() if key != "state_fingerprint_sha256"},
         "components": components}
    )
    if state["state_fingerprint_sha256"] != expected:
        raise ValueError("Velocity portable state failed its current-content fingerprint check.")
    json.dumps(state, ensure_ascii=False, allow_nan=False)
    return state


def _vs_summary(
    *,
    node_id: str,
    status: str,
    methods: str,
    results: str,
    key_results: Mapping[str, Any],
    parameters: Mapping[str, Any],
    references: Sequence[Mapping[str, Any]],
    software_versions: Mapping[str, str],
    report_warnings: Sequence[str],
    limitations: Sequence[str],
) -> dict[str, Any]:
    summary = {
        "schema_version": VELOCITY_SUMMARY_SCHEMA,
        "node_id": node_id,
        "status": status,
        "methods": methods.strip(),
        "results": results.strip(),
        "key_results": _vs_json_value(key_results),
        "parameters": _vs_json_value(parameters),
        "references": _vs_json_value(list(references)),
        "software_versions": dict(software_versions),
        "warnings": [str(value) for value in report_warnings],
        "limitations": [str(value) for value in limitations],
    }
    if not summary["methods"] or not summary["results"] or not summary["references"]:
        raise ValueError("Velocity summary requires methods, results, and references.")
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def prepare_velocity_abundances(
    adata: Any,
    *,
    spliced_layer: str = "spliced",
    unspliced_layer: str = "unspliced",
    min_shared_counts: int = 20,
    min_shared_cells: int = 0,
    normalization_target: str = "median_library",
    overwrite_existing: bool = False,
    scvelo_module: Any | None = None,
) -> tuple[Any, dict[str, Any]]:
    anndata, numpy, _, sparse = _vs_science()
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("Velocity preparation input must be an AnnData object.")
    if bool(adata.isbacked):
        raise ValueError("Velocity preparation requires in-memory AnnData.")
    if int(adata.n_obs) < 2 or int(adata.n_vars) < 1:
        raise ValueError("Velocity preparation requires at least two cells and one feature.")
    _vs_axis(adata.obs_names, label="observation")
    feature_names = _vs_axis(adata.var_names, label="feature")
    spliced_layer = _vs_clean_text(spliced_layer, label="spliced layer")
    unspliced_layer = _vs_clean_text(unspliced_layer, label="unspliced layer")
    if spliced_layer == unspliced_layer:
        raise ValueError("Velocity spliced and unspliced layers must be different.")
    if spliced_layer not in adata.layers or unspliced_layer not in adata.layers:
        raise ValueError("Velocity source spliced/unspliced layers are missing.")
    if adata.layers[spliced_layer] is adata.layers[unspliced_layer]:
        raise ValueError("Velocity spliced and unspliced layers must not alias the same matrix object.")
    min_shared_counts = _vs_integer(
        min_shared_counts, label="min_shared_counts", minimum=0
    )
    min_shared_cells = _vs_integer(min_shared_cells, label="min_shared_cells", minimum=0)
    if normalization_target != "median_library":
        raise ValueError("Velocity normalization_target must be 'median_library'.")
    overwrite_existing = _vs_boolean(overwrite_existing, label="overwrite_existing")
    if VELOCITY_STATE_KEY in adata.uns and not overwrite_existing:
        raise ValueError("Velocity preparation state already exists; enable overwrite_existing to replace it.")
    for canonical, source in (("spliced", spliced_layer), ("unspliced", unspliced_layer)):
        if canonical in adata.layers and canonical != source and not overwrite_existing:
            raise ValueError(f"Velocity canonical layer {canonical!r} already exists.")

    source_matrices = {
        "spliced": adata.layers[spliced_layer],
        "unspliced": adata.layers[unspliced_layer],
    }
    for name, matrix in source_matrices.items():
        values = _vs_matrix_values(
            matrix, label=f"raw {name} layer", shape=(int(adata.n_obs), int(adata.n_vars))
        )
        if values.size and bool((values < 0).any()):
            raise ValueError(f"Velocity raw {name} layer contains negative values.")
        if values.size and not bool(numpy.allclose(values, numpy.rint(values), rtol=0.0, atol=1e-8)):
            raise ValueError(f"Velocity raw {name} layer must contain integer-like molecule counts.")
        if not values.size or not bool((values > 0).any()):
            raise ValueError(f"Velocity raw {name} layer contains no positive counts.")

    spliced = source_matrices["spliced"]
    unspliced = source_matrices["unspliced"]
    if sparse.issparse(spliced) or sparse.issparse(unspliced):
        spliced_csr = sparse.csr_matrix(spliced)
        unspliced_csr = sparse.csr_matrix(unspliced)
        shared_presence = (spliced_csr > 0).multiply(unspliced_csr > 0)
        shared_counts = numpy.asarray(
            shared_presence.multiply(spliced_csr + unspliced_csr).sum(axis=0), dtype=float
        ).ravel()
        shared_cells = numpy.asarray(shared_presence.sum(axis=0), dtype=int).ravel()
    else:
        spliced_array = numpy.asarray(spliced)
        unspliced_array = numpy.asarray(unspliced)
        shared_presence = (spliced_array > 0) & (unspliced_array > 0)
        shared_counts = numpy.where(shared_presence, spliced_array + unspliced_array, 0).sum(axis=0)
        shared_cells = shared_presence.sum(axis=0)
    retained = numpy.ones(int(adata.n_vars), dtype=bool)
    if min_shared_counts > 0:
        retained &= shared_counts >= min_shared_counts
    if min_shared_cells > 0:
        retained &= shared_cells >= min_shared_cells
    if not bool(retained.any()):
        raise ValueError("Velocity filtering removed every feature.")

    output = adata.copy()
    output.layers["spliced"] = source_matrices["spliced"].copy()
    output.layers["unspliced"] = source_matrices["unspliced"].copy()
    if VELOCITY_STATE_KEY in output.uns:
        del output.uns[VELOCITY_STATE_KEY]
    output = output[:, retained].copy()
    before_totals = {
        key: _vs_row_sums(output.layers[key]) for key in ("spliced", "unspliced")
    }
    for key, totals in before_totals.items():
        zero = numpy.flatnonzero(totals <= 0)
        if zero.size:
            examples = [str(output.obs_names[index]) for index in zero[:10]]
            raise ValueError(
                f"Velocity filtered {key} layer has {zero.size} zero-library cell(s): {examples}."
            )
    baseline = output.copy()
    x_fingerprint_before = _vs_matrix_fingerprint(baseline.X, label="X-before-normalization")
    scvelo_module, scvelo_version, normalize = _vs_backend(scvelo_module, "normalize_per_cell")
    caught = []
    try:
        with _vs_preserve_global_state(scvelo_module), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            normalize(
                output,
                counts_per_cell_after=None,
                counts_per_cell=None,
                key_n_counts=None,
                max_proportion_per_cell=None,
                use_initial_size=True,
                layers=["spliced", "unspliced"],
                enforce=True,
                copy=False,
            )
    except Exception as exc:
        raise RuntimeError(
            f"scVelo 0.3.4 velocity-layer normalization failed ({type(exc).__name__}: {exc})."
        ) from exc
    normalized_layers = {}
    for key in ("spliced", "unspliced"):
        if key not in output.layers:
            raise RuntimeError(f"scVelo normalization removed canonical layer {key!r}.")
        normalized_layers[key] = output.layers[key].copy()
    output = baseline
    for key, matrix in normalized_layers.items():
        output.layers[key] = matrix
    if _vs_matrix_fingerprint(output.X, label="X-before-normalization") != x_fingerprint_before:
        raise RuntimeError("Velocity preparation failed to restore expression X exactly.")

    after_totals = {}
    for key in ("spliced", "unspliced"):
        values = _vs_matrix_values(
            output.layers[key], label=f"normalized {key}", shape=(int(output.n_obs), int(output.n_vars))
        )
        if values.size and bool((values < 0).any()):
            raise RuntimeError(f"scVelo returned negative normalized {key} values.")
        totals = _vs_row_sums(output.layers[key])
        target = float(numpy.median(before_totals[key]))
        if not bool(numpy.allclose(totals, target, rtol=1e-6, atol=1e-7)):
            raise RuntimeError(f"scVelo normalized {key} library totals do not match the declared median target.")
        after_totals[key] = totals

    removed_names = [feature_names[index] for index in numpy.flatnonzero(~retained)]
    parameters = {
        "spliced_layer": spliced_layer,
        "unspliced_layer": unspliced_layer,
        "min_shared_counts": min_shared_counts,
        "min_shared_cells": min_shared_cells,
        "normalization_target": normalization_target,
        "overwrite_existing": overwrite_existing,
    }
    state = _vs_write_state(
        output,
        stage="prepared",
        producer_node_id="OpenBioSingleCellVelocityFilterAndNormalize",
        scvelo_version=scvelo_version,
        parameters=parameters,
        parent=None,
        neighbors_key=None,
        vkey=None,
        velocity_mode=None,
        has_dynamics=False,
    )
    summary = _vs_summary(
        node_id="OpenBioSingleCellVelocityFilterAndNormalize",
        status="prepared_velocity_abundances",
        methods=(
            "Validated aligned raw integer-like spliced and unspliced molecule-count layers, retained genes "
            "meeting explicit shared-count/cell thresholds, and normalized only the persistent canonical velocity "
            "layers to their layer-specific median pre-normalization library size with scVelo 0.3.4."
        ),
        results=(
            f"Prepared {output.n_obs:,} cells and retained {output.n_vars:,} of {adata.n_vars:,} genes for "
            "downstream RNA-velocity modeling; no velocity, direction, time, or fate was estimated."
        ),
        key_results={
            "input_cells": int(adata.n_obs),
            "input_genes": int(adata.n_vars),
            "retained_genes": int(output.n_vars),
            "removed_genes": int((~retained).sum()),
            "removed_gene_examples": removed_names[:100],
            "removed_gene_identity_fingerprint_sha256": _vs_json_sha256(removed_names),
            "shared_count_distribution": _vs_numeric_summary(shared_counts),
            "shared_cell_distribution": _vs_numeric_summary(shared_cells),
            "library_totals_before": {
                key: _vs_numeric_summary(value) for key, value in before_totals.items()
            },
            "library_totals_after": {
                key: _vs_numeric_summary(value) for key, value in after_totals.items()
            },
            "state_fingerprint_sha256": state["state_fingerprint_sha256"],
            "canonical_layer_fingerprints_sha256": {
                key: _vs_matrix_fingerprint(output.layers[key], label=key)
                for key in ("spliced", "unspliced")
            },
            "expression_X_restored_exactly": True,
            "backend_side_effects_isolated_to_canonical_layers": True,
        },
        parameters=parameters,
        references=VELOCITY_REFERENCES,
        software_versions=_vs_versions(scvelo_version),
        report_warnings=[f"scVelo warning: {item.message}" for item in caught],
        limitations=[
            "Layer names and integer-like values support but cannot independently prove correct molecule quantification.",
            "This preparation stage does not estimate velocity, temporal order, lineage, transition probability, or fate.",
            "RNA-velocity conclusions remain sensitive to splicing quantification, sampling, kinetics, and batch effects.",
        ],
    )
    return output, summary


def compute_velocity_moments(
    prepared_adata: Any,
    *,
    neighbors_key: str = "neighbors",
    mode: str = "connectivities",
    max_dense_gib: float = 2.0,
    overwrite_existing: bool = False,
    scvelo_module: Any | None = None,
) -> tuple[Any, dict[str, Any]]:
    _, numpy, _, sparse = _vs_science()
    parent = validate_portable_velocity_state(prepared_adata, allowed_stages=("prepared",))
    neighbors_key = _vs_clean_text(neighbors_key, label="neighbors_key")
    if mode not in {"connectivities", "distances"}:
        raise ValueError("Velocity moments mode must be 'connectivities' or 'distances'.")
    max_dense_gib = _vs_float(max_dense_gib, label="max_dense_gib", minimum=0.0, strictly_greater=True)
    overwrite_existing = _vs_boolean(overwrite_existing, label="overwrite_existing")
    collisions = [key for key in ("Ms", "Mu") if key in prepared_adata.layers]
    if collisions and not overwrite_existing:
        raise ValueError(f"Velocity moment layers already exist: {collisions}.")
    dense_bytes = int(prepared_adata.n_obs) * int(prepared_adata.n_vars) * 2 * 4
    if dense_bytes > int(max_dense_gib * 1024**3):
        raise MemoryError(
            f"Velocity moments require at least {dense_bytes / 1024**3:.3f} GiB for dense Ms/Mu, "
            f"exceeding max_dense_gib={max_dense_gib}."
        )
    graph = _vs_resolve_graph(prepared_adata, neighbors_key)
    output = prepared_adata.copy()
    if VELOCITY_STATE_KEY in output.uns:
        del output.uns[VELOCITY_STATE_KEY]
    for key in collisions:
        del output.layers[key]
    abundance_before = {
        key: _vs_matrix_fingerprint(output.layers[key], label=f"before-moments:{key}")
        for key in ("spliced", "unspliced")
    }
    scvelo_module, scvelo_version, moments = _vs_backend(scvelo_module, "moments")
    try:
        with (
            _vs_preserve_global_state(scvelo_module),
            _vs_alias_graph(output, graph),
            warnings.catch_warnings(record=True) as caught,
        ):
            warnings.simplefilter("always")
            moments(
                output,
                n_neighbors=None,
                n_pcs=None,
                mode=mode,
                method="umap",
                use_rep=None,
                use_highly_variable=True,
                copy=False,
            )
    except Exception as exc:
        raise RuntimeError(f"scVelo 0.3.4 moments failed ({type(exc).__name__}: {exc}).") from exc
    warning_text = [str(item.message) for item in caught]
    if any("neighbor" in text.lower() and ("calculat" in text.lower() or "comput" in text.lower()) for text in warning_text):
        raise RuntimeError("scVelo attempted hidden neighbor construction during velocity moments.")
    if any("normaliz" in text.lower() for text in warning_text):
        raise RuntimeError("scVelo attempted hidden normalization during velocity moments.")
    for key, expected in abundance_before.items():
        if _vs_matrix_fingerprint(output.layers[key], label=f"before-moments:{key}") != expected:
            raise RuntimeError(f"scVelo moments unexpectedly modified canonical {key} abundances.")
    post_graph = _vs_resolve_graph(output, neighbors_key)
    if post_graph["fingerprint_sha256"] != graph["fingerprint_sha256"]:
        raise RuntimeError("scVelo moments modified the validated named neighbor graph.")
    diagnostics = {}
    for key in ("Ms", "Mu"):
        if key not in output.layers:
            raise RuntimeError(f"scVelo moments did not create layers[{key!r}].")
        matrix = output.layers[key]
        if sparse.issparse(matrix):
            raise RuntimeError(f"scVelo 0.3.4 moments layer {key!r} must be dense float32.")
        if numpy.asarray(matrix).dtype != numpy.dtype("float32"):
            raise RuntimeError(f"scVelo 0.3.4 moments layer {key!r} must use float32 dtype.")
        values = _vs_matrix_values(
            matrix, label=f"moment {key}", shape=(int(output.n_obs), int(output.n_vars))
        )
        if values.size and bool((values < 0).any()):
            raise RuntimeError(f"scVelo moments returned negative {key} values.")
        diagnostics[key] = {
            "shape": [int(output.n_obs), int(output.n_vars)],
            "dtype": str(numpy.asarray(matrix).dtype),
            "distribution": _vs_numeric_summary(values),
            "fingerprint_sha256": _vs_matrix_fingerprint(matrix, label=key),
        }
    parameters = {
        "neighbors_key": neighbors_key,
        "mode": mode,
        "max_dense_gib": max_dense_gib,
        "overwrite_existing": overwrite_existing,
        "n_neighbors": None,
        "n_pcs": None,
        "use_rep": None,
    }
    state = _vs_write_state(
        output,
        stage="moments",
        producer_node_id="OpenBioSingleCellVelocityMoments",
        scvelo_version=scvelo_version,
        parameters=parameters,
        parent=parent,
        neighbors_key=neighbors_key,
        vkey=None,
        velocity_mode=None,
        has_dynamics=False,
    )
    summary = _vs_summary(
        node_id="OpenBioSingleCellVelocityMoments",
        status="velocity_moments_computed",
        methods=(
            f"Computed first-order spliced and unspliced neighborhood moments with scVelo 0.3.4 using the exact "
            f"named graph {neighbors_key!r} and its {mode} weights; no graph construction, PCA, normalization, or "
            "feature selection was authorized."
        ),
        results=(
            f"Computed observation-by-feature Ms and Mu moment matrices for {output.n_obs:,} cells and "
            f"{output.n_vars:,} genes from the fingerprinted named graph."
        ),
        key_results={
            "input_cells": int(output.n_obs),
            "input_genes": int(output.n_vars),
            "neighbor_graph_fingerprint_sha256": graph["fingerprint_sha256"],
            "neighbor_graph": {
                "neighbors_key": neighbors_key,
                "connectivities_key": graph["connectivities_key"],
                "distances_key": graph["distances_key"],
                "connectivity_nnz": int(graph["connectivities"].nnz),
                "distance_nnz": int(graph["distances"].nnz),
            },
            "dense_bytes": dense_bytes,
            "dense_gib": float(dense_bytes / 1024**3),
            "moments": diagnostics,
            "state_fingerprint_sha256": state["state_fingerprint_sha256"],
        },
        parameters=parameters,
        references=MOMENT_REFERENCES,
        software_versions=_vs_versions(scvelo_version),
        report_warnings=[f"scVelo warning: {text}" for text in warning_text],
        limitations=[
            "Ms and Mu are local graph-smoothed abundances, not independent observations or velocities.",
            "Every later estimate is conditional on the exact representation and named neighbor graph used here.",
            "Dense moment matrices can be memory intensive; the declared preflight covers their minimum payload.",
        ],
    )
    return output, summary


def estimate_rna_velocity(
    velocity_adata: Any,
    *,
    mode: str = "stochastic",
    vkey: str = "velocity",
    min_r2: float = 0.01,
    min_likelihood: float = 0.001,
    overwrite_existing: bool = False,
    scvelo_module: Any | None = None,
) -> tuple[Any, dict[str, Any]]:
    _, numpy, _, sparse = _vs_science()
    if mode not in {"stochastic", "deterministic", "dynamical"}:
        raise ValueError("Velocity mode must be stochastic, deterministic, or dynamical.")
    allowed_stages = ("dynamics_recovered",) if mode == "dynamical" else ("moments",)
    parent = validate_portable_velocity_state(velocity_adata, allowed_stages=allowed_stages)
    vkey = _vs_output_key(vkey, label="vkey")
    min_r2 = _vs_float(min_r2, label="min_r2", minimum=0.0)
    min_likelihood = _vs_float(min_likelihood, label="min_likelihood", minimum=0.0)
    overwrite_existing = _vs_boolean(overwrite_existing, label="overwrite_existing")
    collision_locations = []
    for key in (vkey, f"variance_{vkey}", f"{vkey}_u"):
        if key in velocity_adata.layers:
            collision_locations.append(f"layers[{key!r}]")
    for key in _vs_velocity_var_keys(vkey):
        if key in velocity_adata.var:
            collision_locations.append(f"var[{key!r}]")
    if f"{vkey}_params" in velocity_adata.uns:
        collision_locations.append(f"uns[{f'{vkey}_params'!r}]")
    if collision_locations and not overwrite_existing:
        raise ValueError(f"Velocity estimate output already exists: {collision_locations}.")
    for key in ("Ms", "Mu"):
        values = _vs_matrix_values(
            velocity_adata.layers[key],
            label=f"{key} prerequisite",
            shape=(int(velocity_adata.n_obs), int(velocity_adata.n_vars)),
        )
        if values.size and bool((values < 0).any()):
            raise ValueError(f"Velocity {key} prerequisite contains negative values.")
    if mode == "dynamical":
        required_fit = {
            "fit_alpha",
            "fit_beta",
            "fit_gamma",
            "fit_t_",
            "fit_scaling",
            "fit_likelihood",
            "fit_r2",
            "openbio_dynamics_selected",
            "openbio_dynamics_fit_success",
        }
        missing_fit = sorted(required_fit - set(velocity_adata.var.columns))
        if missing_fit:
            raise ValueError(
                "Dynamical velocity requires an explicit complete recovered fit including fit_r2; "
                f"missing={missing_fit}. Run a steady-state velocity estimate before recovery when fit_r2 is absent."
            )
        recover_record = velocity_adata.uns.get("recover_dynamics")
        if not isinstance(recover_record, Mapping):
            raise ValueError("Dynamical velocity requires verified recover_dynamics parameters.")
        successful = velocity_adata.var["openbio_dynamics_fit_success"].to_numpy(dtype=bool)
        if int(successful.sum()) < 10:
            raise ValueError("Dynamical velocity requires at least 10 successfully recovered genes.")
        for key in required_fit - {"openbio_dynamics_selected", "openbio_dynamics_fit_success"}:
            values = velocity_adata.var[key].to_numpy(dtype=float)
            if bool(numpy.isinf(values).any()):
                raise ValueError(f"Dynamical velocity prerequisite {key!r} contains infinity.")
            if not bool(numpy.isfinite(values[successful]).all()):
                raise ValueError(f"Dynamical velocity prerequisite {key!r} is missing for successful genes.")
    model_graph = (
        _vs_resolve_graph(velocity_adata, str(parent["neighbors_key"]))
        if mode == "dynamical"
        else None
    )

    output = velocity_adata.copy()
    if VELOCITY_STATE_KEY in output.uns:
        del output.uns[VELOCITY_STATE_KEY]
    if overwrite_existing:
        for key in (vkey, f"variance_{vkey}", f"{vkey}_u"):
            if key in output.layers:
                del output.layers[key]
        for key in _vs_velocity_var_keys(vkey):
            if key in output.var:
                del output.var[key]
        if f"{vkey}_params" in output.uns:
            del output.uns[f"{vkey}_params"]
    obs_before = output.obs.copy(deep=True)
    var_before = output.var.copy(deep=True)
    scvelo_module, scvelo_version, velocity = _vs_backend(scvelo_module, "velocity")
    compatibility_actions: tuple[str, ...] = ()
    try:
        with (
            _vs_velocity_pandas3_compatibility(
                scvelo_module, enabled=mode == "dynamical"
            ) as compatibility_actions,
            _vs_preserve_global_state(scvelo_module),
            _vs_optional_alias_graph(output, model_graph),
            warnings.catch_warnings(record=True) as caught,
        ):
            warnings.simplefilter("always")
            velocity(
                output,
                vkey=vkey,
                mode=mode,
                fit_offset=False,
                fit_offset2=False,
                filter_genes=False,
                groups=None,
                groupby=None,
                groups_for_fit=None,
                constrain_ratio=None,
                use_raw=False,
                use_latent_time=None,
                perc=None,
                min_r2=min_r2,
                min_likelihood=min_likelihood,
                r2_adjusted=None,
                use_highly_variable=True,
                diff_kinetics=None,
                copy=False,
            )
    except Exception as exc:
        raise RuntimeError(f"scVelo 0.3.4 velocity estimation failed ({type(exc).__name__}: {exc}).") from exc
    generated_var = {
        key: output.var[key].copy(deep=True)
        for key in _vs_velocity_var_keys(vkey)
        if key in output.var
    }
    output.obs = obs_before
    output.var = var_before
    for key, values in generated_var.items():
        output.var[key] = values
    if model_graph is not None:
        post_graph = _vs_resolve_graph(output, str(parent["neighbors_key"]))
        if post_graph["fingerprint_sha256"] != model_graph["fingerprint_sha256"]:
            raise RuntimeError("Dynamical velocity modified the validated named neighbor graph.")
        if (
            _vs_dynamics_components(output)["fingerprint_sha256"]
            != parent["dynamics_fit_fingerprint_sha256"]
        ):
            raise RuntimeError("Dynamical velocity modified the validated recovered-dynamics fit.")
    if vkey not in output.layers:
        raise RuntimeError(f"scVelo velocity did not create layers[{vkey!r}].")
    params = output.uns.get(f"{vkey}_params")
    if not isinstance(params, Mapping) or params.get("mode") != mode:
        raise RuntimeError(
            f"scVelo velocity backend fell back or reported a different model; requested={mode!r}, "
            f"observed={params.get('mode') if isinstance(params, Mapping) else None!r}."
        )
    mask_key = f"{vkey}_genes"
    if mask_key not in output.var:
        raise RuntimeError(f"scVelo velocity did not create var[{mask_key!r}].")
    mask_series = output.var[mask_key]
    if not (mask_series.dtype == bool or str(mask_series.dtype) == "boolean") or bool(mask_series.isna().any()):
        raise RuntimeError("scVelo velocity-gene mask must be complete and boolean.")
    selected = mask_series.to_numpy(dtype=bool)
    selected_count = int(selected.sum())
    if selected_count < 10:
        raise ValueError(
            f"Velocity estimation selected only {selected_count} genes; at least 10 are required for reviewable evidence."
        )
    layer = output.layers[vkey]
    all_values = _vs_matrix_values(
        layer, label=f"velocity layer {vkey}", shape=(int(output.n_obs), int(output.n_vars)), finite=False
    )
    if bool(numpy.isinf(all_values).any()):
        raise RuntimeError("scVelo velocity layer contains infinity.")
    selected_values = (
        layer[:, selected].toarray() if sparse.issparse(layer) else numpy.asarray(layer)[:, selected]
    )
    selected_values = numpy.asarray(selected_values, dtype=float)
    if not bool(numpy.isfinite(selected_values).all()):
        raise RuntimeError("scVelo velocity layer is non-finite for selected velocity genes.")
    nonselected_missing = int(numpy.isnan(numpy.asarray(all_values, dtype=float)).sum())
    magnitudes = numpy.sqrt(numpy.mean(numpy.square(selected_values), axis=1))
    parameter_diagnostics = {}
    for key in (f"{vkey}_r2", f"{vkey}_beta", f"{vkey}_gamma"):
        if key in output.var:
            values = output.var[key].to_numpy(dtype=float)
            if bool(numpy.isinf(values).any()):
                raise RuntimeError(f"scVelo velocity parameter {key!r} contains infinity.")
            parameter_diagnostics[key] = _vs_numeric_summary(values)
    if mode == "dynamical":
        parameter_diagnostics["fit_likelihood"] = _vs_numeric_summary(
            output.var["fit_likelihood"].to_numpy(dtype=float)
        )
    parameters = {
        "mode": mode,
        "vkey": vkey,
        "min_r2": min_r2,
        "min_likelihood": min_likelihood,
        "overwrite_existing": overwrite_existing,
        "fit_offset": False,
        "fit_offset2": False,
        "filter_genes": False,
        "use_raw": False,
        "use_latent_time": None,
        "pandas3_compatibility_policy": (
            ["writable_divergence_inputs"] if mode == "dynamical" else []
        ),
    }
    state = _vs_write_state(
        output,
        stage="velocity_estimated",
        producer_node_id="OpenBioSingleCellEstimateVelocity",
        scvelo_version=scvelo_version,
        parameters=parameters,
        parent=parent,
        neighbors_key=parent["neighbors_key"],
        vkey=vkey,
        velocity_mode=mode,
        has_dynamics=parent["dynamics_fit_fingerprint_sha256"] is not None,
    )
    caught_text = [str(item.message) for item in caught]
    if any("falling back" in text.lower() for text in caught_text):
        raise RuntimeError("scVelo velocity emitted a model fallback warning.")
    summary = _vs_summary(
        node_id="OpenBioSingleCellEstimateVelocity",
        status="exploratory_model_dependent_velocity_estimate",
        methods=(
            f"Estimated per-gene RNA velocity with the explicit scVelo 0.3.4 {mode} kinetic model on validated "
            "Ms/Mu moments, with fixed offsets disabled, no grouping, no latent-time regularization, no raw-layer "
            "fallback, and no feature-axis filtering."
        ),
        results=(
            f"The {mode} model selected {selected_count:,} of {output.n_vars:,} genes and produced finite "
            "velocity values for every selected gene; this is exploratory vector-field evidence, not fate or time."
        ),
        key_results={
            "input_cells": int(output.n_obs),
            "input_genes": int(output.n_vars),
            "selected_velocity_genes": selected_count,
            "failed_or_unselected_genes": int(output.n_vars - selected_count),
            "nonselected_missing_velocity_values": nonselected_missing,
            "velocity_magnitude_distribution": _vs_numeric_summary(magnitudes),
            "parameter_diagnostics": parameter_diagnostics,
            "velocity_layer_fingerprint_sha256": _vs_matrix_fingerprint(layer, label=f"velocity:{vkey}"),
            "velocity_gene_fingerprint_sha256": _vs_series_fingerprint(
                output.var, [mask_key], label=f"velocity-mask:{mask_key}"
            ),
            "upstream_neighbor_graph_fingerprint_sha256": parent[
                "neighbor_graph_fingerprint_sha256"
            ],
            "dynamics_fit_fingerprint_sha256": state["dynamics_fit_fingerprint_sha256"],
            "state_fingerprint_sha256": state["state_fingerprint_sha256"],
            "pandas3_compatibility_actions_applied": list(compatibility_actions),
        },
        parameters=parameters,
        references=DYNAMICS_REFERENCES if mode == "dynamical" else VELOCITY_REFERENCES,
        software_versions=_vs_versions(scvelo_version),
        report_warnings=[
            *(
                [
                    "Applied and restored the audited scVelo 0.3.4 writable-divergence-input "
                    "compatibility seam for Pandas 3."
                ]
                if compatibility_actions
                else []
            ),
            *[f"scVelo warning: {text}" for text in caught_text],
        ],
        limitations=[
            "RNA velocity is conditional on a kinetic model and can be unstable under quantification error, bursts, batch effects, or violated steady/transient assumptions.",
            "Velocity vectors are not transition probabilities, fate probabilities, measured time, lineage proof, or replicate-aware Condition inference.",
            "Genes outside the selected mask are excluded from directional evidence; missing dynamical values there are explicitly counted.",
        ],
    )
    return output, summary


def recover_velocity_dynamics(
    velocity_adata: Any,
    *,
    gene_selection: str = "velocity_genes",
    n_top_genes: int = 0,
    max_iter: int = 10,
    n_jobs: int = 1,
    max_dense_gib: float = 2.0,
    overwrite_existing: bool = False,
    scvelo_module: Any | None = None,
) -> tuple[Any, dict[str, Any]]:
    _, numpy, _, sparse = _vs_science()
    parent = validate_portable_velocity_state(
        velocity_adata, allowed_stages=("moments", "velocity_estimated")
    )
    if parent["stage"] == "velocity_estimated" and parent["velocity_mode"] == "dynamical":
        raise ValueError("Recover Dynamics cannot consume an already dynamical velocity estimate.")
    if gene_selection not in {"velocity_genes", "all"}:
        raise ValueError("Velocity dynamics gene_selection must be 'velocity_genes' or 'all'.")
    n_top_genes = _vs_integer(n_top_genes, label="n_top_genes", minimum=0)
    max_iter = _vs_integer(max_iter, label="max_iter", minimum=1, maximum=10_000)
    n_jobs = _vs_integer(n_jobs, label="n_jobs", minimum=1, maximum=1024)
    max_dense_gib = _vs_float(max_dense_gib, label="max_dense_gib", minimum=0.0, strictly_greater=True)
    overwrite_existing = _vs_boolean(overwrite_existing, label="overwrite_existing")

    if gene_selection == "velocity_genes":
        upstream_vkey = parent.get("vkey")
        mask_key = f"{upstream_vkey}_genes" if upstream_vkey else None
        if mask_key is None or mask_key not in velocity_adata.var:
            raise ValueError(
                "Recover Dynamics gene_selection='velocity_genes' requires a verified upstream steady-state "
                "velocity-gene mask."
            )
        mask = velocity_adata.var[mask_key]
        if not (mask.dtype == bool or str(mask.dtype) == "boolean") or bool(mask.isna().any()):
            raise ValueError("Recover Dynamics requires a complete boolean upstream velocity-gene mask.")
        selected_positions = numpy.flatnonzero(mask.to_numpy(dtype=bool))
    else:
        selected_positions = numpy.arange(int(velocity_adata.n_vars), dtype=int)
    if n_top_genes > 0 and selected_positions.size > n_top_genes:
        matrix = velocity_adata.layers["Ms"][:, selected_positions]
        totals = (
            numpy.asarray(matrix.sum(axis=0), dtype=float).ravel()
            if sparse.issparse(matrix)
            else numpy.asarray(matrix, dtype=float).sum(axis=0)
        )
        order = numpy.lexsort((selected_positions, -totals))
        selected_positions = selected_positions[order[:n_top_genes]]
    selected_positions = numpy.sort(selected_positions)
    if selected_positions.size < 5:
        raise ValueError("Recover Dynamics requires at least five explicitly selected genes.")
    selected_genes = [str(velocity_adata.var_names[position]) for position in selected_positions]
    dense_bytes = int(velocity_adata.n_obs) * int(velocity_adata.n_vars) * 3 * 8
    if dense_bytes > int(max_dense_gib * 1024**3):
        raise MemoryError(
            f"Recover Dynamics requires at least {dense_bytes / 1024**3:.3f} GiB for fit_t/fit_tau/fit_tau_, "
            f"exceeding max_dense_gib={max_dense_gib}."
        )
    collision_locations = [
        f"var[{column!r}]" for column in velocity_adata.var.columns if str(column).startswith("fit_")
    ]
    collision_locations.extend(
        f"var[{column!r}]"
        for column in (
            "openbio_dynamics_selected",
            "openbio_dynamics_fit_success",
            "openbio_dynamics_velocity_gene",
        )
        if column in velocity_adata.var
    )
    collision_locations.extend(
        f"layers[{key!r}]" for key in velocity_adata.layers.keys() if str(key).startswith("fit_")
    )
    if "loss" in velocity_adata.varm:
        collision_locations.append("varm['loss']")
    if "recover_dynamics" in velocity_adata.uns:
        collision_locations.append("uns['recover_dynamics']")
    if collision_locations and not overwrite_existing:
        raise ValueError(f"Recover Dynamics output already exists: {collision_locations[:20]}.")
    graph = _vs_resolve_graph(velocity_adata, str(parent["neighbors_key"]))
    output = velocity_adata.copy()
    if VELOCITY_STATE_KEY in output.uns:
        del output.uns[VELOCITY_STATE_KEY]
    if overwrite_existing:
        for column in list(output.var.columns):
            if str(column).startswith("fit_") or column in {
                "openbio_dynamics_selected",
                "openbio_dynamics_fit_success",
                "openbio_dynamics_velocity_gene",
            }:
                del output.var[column]
        for key in list(output.layers.keys()):
            if str(key).startswith("fit_"):
                del output.layers[key]
        if "loss" in output.varm:
            del output.varm["loss"]
        if "recover_dynamics" in output.uns:
            del output.uns["recover_dynamics"]
    upstream_vkey = parent.get("vkey")
    if upstream_vkey and f"{upstream_vkey}_r2" in output.var and "fit_r2" not in output.var:
        output.var["fit_r2"] = output.var[f"{upstream_vkey}_r2"].to_numpy(dtype=float, copy=True)
    selected_mask = numpy.zeros(int(output.n_vars), dtype=bool)
    selected_mask[selected_positions] = True
    output.var["openbio_dynamics_selected"] = selected_mask
    if upstream_vkey and f"{upstream_vkey}_genes" in output.var:
        output.var["openbio_dynamics_velocity_gene"] = output.var[
            f"{upstream_vkey}_genes"
        ].to_numpy(dtype=bool, copy=True)
    else:
        output.var["openbio_dynamics_velocity_gene"] = numpy.zeros(int(output.n_vars), dtype=bool)
    scvelo_module, scvelo_version, recover = _vs_backend(scvelo_module, "recover_dynamics")
    compatibility_actions: tuple[str, ...] = ()
    try:
        with (
            _vs_recover_pandas3_compatibility(recover) as compatibility_actions,
            _vs_preserve_global_state(scvelo_module),
            _vs_alias_graph(output, graph),
            warnings.catch_warnings(record=True) as caught,
        ):
            warnings.simplefilter("always")
            recover(
                output,
                var_names=selected_genes,
                n_top_genes=None,
                max_iter=max_iter,
                assignment_mode="projection",
                t_max=None,
                fit_time=True,
                fit_scaling=True,
                fit_steady_states=True,
                fit_connected_states=True,
                fit_basal_transcription=False,
                use_raw=False,
                load_pars=False,
                return_model=False,
                plot_results=False,
                steady_state_prior=None,
                add_key="fit",
                copy=False,
                n_jobs=n_jobs,
                backend="loky",
                show_progress_bar=False,
            )
    except Exception as exc:
        raise RuntimeError(f"scVelo 0.3.4 recover_dynamics failed ({type(exc).__name__}: {exc}).") from exc
    post_graph = _vs_resolve_graph(output, str(parent["neighbors_key"]))
    if post_graph["fingerprint_sha256"] != graph["fingerprint_sha256"]:
        raise RuntimeError("Recover Dynamics modified the validated named neighbor graph.")
    required_parameters = (
        "fit_alpha",
        "fit_beta",
        "fit_gamma",
        "fit_t_",
        "fit_scaling",
        "fit_std_u",
        "fit_std_s",
        "fit_likelihood",
        "fit_u0",
        "fit_s0",
        "fit_pval_steady",
        "fit_steady_u",
        "fit_steady_s",
        "fit_variance",
        "fit_alignment_scaling",
    )
    missing = [key for key in required_parameters if key not in output.var]
    missing.extend(key for key in ("fit_t", "fit_tau", "fit_tau_") if key not in output.layers)
    if missing:
        raise RuntimeError(f"scVelo Recover Dynamics returned an incomplete fit bundle: {missing}.")
    recover_record = output.uns.get("recover_dynamics")
    if not isinstance(recover_record, Mapping) or not (
        recover_record.get("fit_connected_states") is True
        and recover_record.get("fit_basal_transcription") is False
        and recover_record.get("use_raw") is False
    ):
        raise RuntimeError("scVelo Recover Dynamics parameter metadata differs from the audited policy.")
    parameter_arrays = {}
    for key in required_parameters:
        values = output.var[key].to_numpy(dtype=float)
        if bool(numpy.isinf(values).any()):
            raise RuntimeError(f"Recovered parameter {key!r} contains infinity.")
        parameter_arrays[key] = values
    for key in ("fit_t", "fit_tau", "fit_tau_"):
        values = _vs_matrix_values(
            output.layers[key],
            label=f"recovered layer {key}",
            shape=(int(output.n_obs), int(output.n_vars)),
            finite=False,
        )
        if bool(numpy.isinf(values).any()):
            raise RuntimeError(f"Recovered layer {key!r} contains infinity.")
    success = selected_mask.copy()
    for key in ("fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_scaling", "fit_likelihood"):
        success &= numpy.isfinite(parameter_arrays[key])
    likelihood = parameter_arrays["fit_likelihood"]
    if bool((likelihood[success] < 0.0).any()):
        raise RuntimeError("Recovered fit_likelihood must be nonnegative for successful genes.")
    success_count = int(success.sum())
    if success_count == 0:
        raise ValueError("Recover Dynamics did not successfully fit any selected gene.")
    output.var["openbio_dynamics_fit_success"] = success
    # A recovered state owns the fit, not the stale steady-state velocity used to nominate genes.
    if upstream_vkey:
        for key in (upstream_vkey, f"variance_{upstream_vkey}", f"{upstream_vkey}_u"):
            if key in output.layers:
                del output.layers[key]
        for key in _vs_velocity_var_keys(str(upstream_vkey)):
            if key in output.var:
                del output.var[key]
        if f"{upstream_vkey}_params" in output.uns:
            del output.uns[f"{upstream_vkey}_params"]
    parameters = {
        "gene_selection": gene_selection,
        "n_top_genes": n_top_genes,
        "max_iter": max_iter,
        "n_jobs": n_jobs,
        "max_dense_gib": max_dense_gib,
        "overwrite_existing": overwrite_existing,
        "assignment_mode": "projection",
        "fit_connected_states": True,
        "fit_basal_transcription": False,
        "use_raw": False,
        "backend": "loky",
        "pandas3_compatibility_policy": ["stable_unique", "writable_parameter_arrays"],
    }
    state = _vs_write_state(
        output,
        stage="dynamics_recovered",
        producer_node_id="OpenBioSingleCellRecoverDynamics",
        scvelo_version=scvelo_version,
        parameters=parameters,
        parent=parent,
        neighbors_key=parent["neighbors_key"],
        vkey=None,
        velocity_mode=None,
        has_dynamics=True,
    )
    failed_genes = [
        str(output.var_names[index]) for index in numpy.flatnonzero(selected_mask & ~success)
    ]
    summary = _vs_summary(
        node_id="OpenBioSingleCellRecoverDynamics",
        status="exploratory_dynamical_parameters_recovered",
        methods=(
            "Fitted scVelo 0.3.4's explicit expectation-maximization dynamical splicing model to a "
            "preflighted deterministic gene set on validated Ms/Mu moments and the exact named graph, with "
            "projection assignment, connected-state fitting, fixed kinetic policy, and no raw or hidden selection fallback."
        ),
        results=(
            f"Recovered complete finite core kinetic parameters for {success_count:,} of "
            f"{selected_positions.size:,} selected genes; failed fits remain explicitly missing and are not rates measured directly."
        ),
        key_results={
            "input_cells": int(output.n_obs),
            "input_genes": int(output.n_vars),
            "selected_genes": int(selected_positions.size),
            "selected_gene_fingerprint_sha256": _vs_json_sha256(selected_genes),
            "selected_gene_examples": selected_genes[:100],
            "successful_fits": success_count,
            "failed_fits": int(selected_positions.size - success_count),
            "failed_gene_examples": failed_genes[:100],
            "fit_likelihood_distribution": _vs_numeric_summary(likelihood[selected_mask]),
            "fit_r2_distribution": (
                _vs_numeric_summary(output.var["fit_r2"].to_numpy(dtype=float)[selected_mask])
                if "fit_r2" in output.var
                else None
            ),
            "kinetic_parameter_distributions": {
                key: _vs_numeric_summary(values[selected_mask])
                for key, values in parameter_arrays.items()
                if key in {"fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_scaling"}
            },
            "dense_bytes_preflight": dense_bytes,
            "neighbor_graph_fingerprint_sha256": graph["fingerprint_sha256"],
            "dynamics_fit_fingerprint_sha256": state["dynamics_fit_fingerprint_sha256"],
            "state_fingerprint_sha256": state["state_fingerprint_sha256"],
            "fit_r2_source": (
                "upstream_steady_state_velocity_r2" if "fit_r2" in output.var else "not_available"
            ),
            "pandas3_compatibility_actions_applied": list(compatibility_actions),
        },
        parameters=parameters,
        references=DYNAMICS_REFERENCES,
        software_versions=_vs_versions(scvelo_version),
        report_warnings=[
            *(
                [
                    "Applied and restored the audited scVelo 0.3.4 stable-unique and writable-parameter "
                    "compatibility seams for Pandas 3."
                ]
                if compatibility_actions
                else []
            ),
            *[f"scVelo warning: {item.message}" for item in caught],
        ],
        limitations=[
            "Recovered kinetic parameters are model fits, not directly measured transcription, splicing, or degradation rates.",
            "Fit likelihood is model fit quality, not a p-value, differential-expression result, driver score, or fate probability.",
            "Recover Dynamics does not compute dynamical velocity, a velocity graph, latent time, lineage, or CellRank fate.",
        ],
    )
    return output, summary


def _vs_velocity_graph_matrix(
    matrix: Any,
    *,
    n_obs: int,
    label: str,
    sign: str,
) -> Any:
    _, numpy, _, sparse = _vs_science()
    if not sparse.issparse(matrix):
        raise TypeError(f"Velocity {label} directed graph must be SciPy sparse.")
    if tuple(matrix.shape) != (n_obs, n_obs):
        raise ValueError(f"Velocity {label} directed graph must have shape ({n_obs}, {n_obs}).")
    graph = matrix.tocsr(copy=True)
    graph.sum_duplicates()
    graph.eliminate_zeros()
    graph.sort_indices()
    values = _vs_matrix_values(graph, label=f"{label} directed graph", shape=(n_obs, n_obs))
    if not bool(numpy.allclose(graph.diagonal(), 0.0, rtol=0.0, atol=1e-12)):
        raise ValueError(f"Velocity {label} directed graph must have a zero diagonal.")
    if sign == "positive" and values.size and bool(((values < 0.0) | (values > 1.0)).any()):
        raise ValueError("Velocity positive directed-graph weights must lie in [0, 1].")
    if sign == "negative" and values.size and bool(((values < -1.0) | (values > 0.0)).any()):
        raise ValueError("Velocity negative directed-graph weights must lie in [-1, 0].")
    return graph


def build_velocity_graph(
    velocity_adata: Any,
    *,
    vkey: str = "velocity",
    xkey: str = "Ms",
    mode_neighbors: str = "distances",
    n_jobs: int = 1,
    overwrite_existing: bool = False,
    scvelo_module: Any | None = None,
) -> tuple[Any, dict[str, Any]]:
    _, numpy, _, sparse = _vs_science()
    parent = validate_portable_velocity_state(
        velocity_adata, allowed_stages=("velocity_estimated",)
    )
    vkey = _vs_output_key(vkey, label="vkey")
    if vkey != parent["vkey"]:
        raise ValueError(
            f"Velocity graph vkey {vkey!r} differs from the staged estimate {parent['vkey']!r}."
        )
    if xkey != "Ms":
        raise ValueError("Velocity graph xkey is fixed to the verified moment layer 'Ms'.")
    if mode_neighbors not in {"distances", "connectivities"}:
        raise ValueError("Velocity graph mode_neighbors must be 'distances' or 'connectivities'.")
    n_jobs = _vs_integer(n_jobs, label="n_jobs", minimum=1, maximum=1024)
    overwrite_existing = _vs_boolean(overwrite_existing, label="overwrite_existing")
    graph_key = f"{vkey}_graph"
    negative_key = f"{vkey}_graph_neg"
    self_key = f"{vkey}_self_transition"
    collisions = []
    if graph_key in velocity_adata.uns:
        collisions.append(f"uns[{graph_key!r}]")
    if negative_key in velocity_adata.uns:
        collisions.append(f"uns[{negative_key!r}]")
    if self_key in velocity_adata.obs:
        collisions.append(f"obs[{self_key!r}]")
    if collisions and not overwrite_existing:
        raise ValueError(f"Velocity graph output already exists: {collisions}.")
    mask_key = f"{vkey}_genes"
    mask = velocity_adata.var[mask_key]
    if not (mask.dtype == bool or str(mask.dtype) == "boolean") or bool(mask.isna().any()):
        raise ValueError("Velocity graph requires a complete boolean velocity-gene mask.")
    selected = mask.to_numpy(dtype=bool)
    if int(selected.sum()) < 10:
        raise ValueError("Velocity graph requires at least 10 selected velocity genes.")
    velocity_layer = velocity_adata.layers[vkey]
    selected_velocity = (
        velocity_layer[:, selected].toarray()
        if sparse.issparse(velocity_layer)
        else numpy.asarray(velocity_layer)[:, selected]
    )
    if not bool(numpy.isfinite(numpy.asarray(selected_velocity, dtype=float)).all()):
        raise ValueError("Velocity graph requires finite velocities for selected genes.")
    _vs_matrix_values(
        velocity_adata.layers[xkey],
        label=f"velocity graph expression layer {xkey}",
        shape=(int(velocity_adata.n_obs), int(velocity_adata.n_vars)),
    )
    source_graph = _vs_resolve_graph(velocity_adata, str(parent["neighbors_key"]))
    output = velocity_adata.copy()
    if VELOCITY_STATE_KEY in output.uns:
        del output.uns[VELOCITY_STATE_KEY]
    if overwrite_existing:
        for key in (graph_key, negative_key):
            if key in output.uns:
                del output.uns[key]
        if self_key in output.obs:
            del output.obs[self_key]
    scvelo_module, scvelo_version, velocity_graph = _vs_backend(scvelo_module, "velocity_graph")
    try:
        with (
            _vs_preserve_global_state(scvelo_module),
            _vs_alias_graph(output, source_graph),
            warnings.catch_warnings(record=True) as caught,
        ):
            warnings.simplefilter("always")
            velocity_graph(
                output,
                vkey=vkey,
                xkey=xkey,
                tkey=None,
                basis=None,
                n_neighbors=None,
                n_recurse_neighbors=1,
                random_neighbors_at_max=None,
                sqrt_transform=False,
                variance_stabilization=None,
                gene_subset=None,
                compute_uncertainties=False,
                approx=False,
                mode_neighbors=mode_neighbors,
                copy=False,
                n_jobs=n_jobs,
                backend="loky",
                show_progress_bar=False,
            )
    except Exception as exc:
        raise RuntimeError(f"scVelo 0.3.4 velocity_graph failed ({type(exc).__name__}: {exc}).") from exc
    post_source_graph = _vs_resolve_graph(output, str(parent["neighbors_key"]))
    if post_source_graph["fingerprint_sha256"] != source_graph["fingerprint_sha256"]:
        raise RuntimeError("scVelo velocity_graph modified the validated source neighbor graph.")
    if graph_key not in output.uns or negative_key not in output.uns:
        raise RuntimeError("scVelo velocity_graph did not return both positive and negative directed graphs.")
    positive = _vs_velocity_graph_matrix(
        output.uns[graph_key], n_obs=int(output.n_obs), label="positive", sign="positive"
    )
    negative = _vs_velocity_graph_matrix(
        output.uns[negative_key], n_obs=int(output.n_obs), label="negative", sign="negative"
    )
    overlap = positive.multiply(negative != 0)
    overlap.eliminate_zeros()
    if overlap.nnz:
        raise RuntimeError("Positive and negative velocity graphs must have disjoint directed-edge support.")
    source_support = source_graph[
        "distances" if mode_neighbors == "distances" else "connectivities"
    ]
    combined_support = (positive != 0).astype(int) + (negative != 0).astype(int)
    outside_support = combined_support - combined_support.multiply(source_support != 0)
    outside_support.eliminate_zeros()
    if outside_support.nnz:
        raise RuntimeError(
            "Velocity directed edges escaped the selected direct named-neighbor support."
        )
    output.uns[graph_key] = positive
    output.uns[negative_key] = negative
    params = output.uns.get(f"{vkey}_params")
    if not isinstance(params, Mapping) or params.get("mode") != parent["velocity_mode"]:
        raise RuntimeError("Velocity graph lost the exact upstream velocity-model parameters.")
    if params.get("mode_neighbors") != mode_neighbors:
        raise RuntimeError("Velocity graph backend reported a different neighbor mode.")
    if params.get("n_recurse_neighbors") != 1:
        raise RuntimeError("Velocity graph backend did not preserve the fixed direct-neighbor policy.")
    if self_key not in output.obs:
        raise RuntimeError("Velocity graph backend did not return self-transition evidence.")
    self_transition = output.obs[self_key].to_numpy(dtype=float)
    if not bool(numpy.isfinite(self_transition).all()) or bool(
        ((self_transition < 0.0) | (self_transition > 1.0)).any()
    ):
        raise RuntimeError("Velocity self-transition evidence must be finite and lie in [0, 1].")
    positive_out = numpy.asarray((positive != 0).sum(axis=1), dtype=int).ravel()
    positive_in = numpy.asarray((positive != 0).sum(axis=0), dtype=int).ravel()
    negative_out = numpy.asarray((negative != 0).sum(axis=1), dtype=int).ravel()
    parameters = {
        "vkey": vkey,
        "xkey": xkey,
        "mode_neighbors": mode_neighbors,
        "n_jobs": n_jobs,
        "overwrite_existing": overwrite_existing,
        "approx": False,
        "compute_uncertainties": False,
        "n_recurse_neighbors": 1,
        "sqrt_transform": False,
        "backend": "loky",
    }
    state = _vs_write_state(
        output,
        stage="velocity_graph",
        producer_node_id="OpenBioSingleCellVelocityGraph",
        scvelo_version=scvelo_version,
        parameters=parameters,
        parent=parent,
        neighbors_key=parent["neighbors_key"],
        vkey=vkey,
        velocity_mode=parent["velocity_mode"],
        has_dynamics=parent["dynamics_fit_fingerprint_sha256"] is not None,
    )
    summary = _vs_summary(
        node_id="OpenBioSingleCellVelocityGraph",
        status="directed_velocity_correlation_evidence",
        methods=(
            "Compared each selected-gene velocity vector with candidate cell-state changes using scVelo 0.3.4 "
            f"cosine correlations over the exact direct {mode_neighbors} neighbor support, full Ms space, fixed "
            "sqrt_transform=False, no approximation, and no uncertainty, recursive expansion, embedding, or hidden velocity fallback."
        ),
        results=(
            f"Constructed {positive.nnz:,} positive and {negative.nnz:,} negative directed cosine-alignment edges "
            f"for {output.n_obs:,} cells; these matrices are transition evidence, not probabilities or fates."
        ),
        key_results={
            "input_cells": int(output.n_obs),
            "selected_velocity_genes": int(selected.sum()),
            "velocity_mode": parent["velocity_mode"],
            "source_neighbor_graph_fingerprint_sha256": source_graph["fingerprint_sha256"],
            "positive_graph": {
                "shape": [int(output.n_obs), int(output.n_obs)],
                "nnz": int(positive.nnz),
                "weight_distribution": _vs_numeric_summary(positive.data),
                "out_degree_distribution": _vs_numeric_summary(positive_out),
                "in_degree_distribution": _vs_numeric_summary(positive_in),
                "zero_outdegree_cells": int((positive_out == 0).sum()),
                "fingerprint_sha256": state[
                    "velocity_graph_positive_fingerprint_sha256"
                ],
            },
            "negative_graph": {
                "shape": [int(output.n_obs), int(output.n_obs)],
                "nnz": int(negative.nnz),
                "weight_distribution": _vs_numeric_summary(negative.data),
                "out_degree_distribution": _vs_numeric_summary(negative_out),
                "zero_outdegree_cells": int((negative_out == 0).sum()),
                "fingerprint_sha256": state[
                    "velocity_graph_negative_fingerprint_sha256"
                ],
            },
            "self_transition_distribution": _vs_numeric_summary(self_transition),
            "dynamics_fit_fingerprint_sha256": state["dynamics_fit_fingerprint_sha256"],
            "state_fingerprint_sha256": state["state_fingerprint_sha256"],
        },
        parameters=parameters,
        references=VELOCITY_REFERENCES,
        software_versions=_vs_versions(scvelo_version),
        report_warnings=[f"scVelo warning: {item.message}" for item in caught],
        limitations=[
            "Cosine-alignment edges are directed kinetic-model evidence, not a row-stochastic transition matrix.",
            "Graph direction is conditional on quantification, moments, velocity-gene selection, kinetic model, and neighbor support.",
            "This node does not run a CellRank kernel or estimator and does not infer terminal states, fates, or lineage confidence.",
        ],
    )
    return output, summary


def rank_recovered_dynamics(
    fitted_adata: Any,
    *,
    top_n: int = 50,
    include_failed: bool = False,
    max_output_rows: int = 100_000,
) -> tuple[Any, dict[str, Any]]:
    _, numpy, pandas, _ = _vs_science()
    state = validate_portable_velocity_state(
        fitted_adata,
        allowed_stages=("dynamics_recovered", "velocity_estimated", "velocity_graph"),
    )
    if state["dynamics_fit_fingerprint_sha256"] is None:
        raise ValueError("Recovered Dynamics Fit Ranking requires a state with a verified dynamics fit.")
    top_n = _vs_integer(top_n, label="top_n", minimum=1)
    include_failed = _vs_boolean(include_failed, label="include_failed")
    max_output_rows = _vs_integer(max_output_rows, label="max_output_rows", minimum=1)
    required = {
        "openbio_dynamics_selected",
        "openbio_dynamics_fit_success",
        "openbio_dynamics_velocity_gene",
        "fit_likelihood",
        "fit_alpha",
        "fit_beta",
        "fit_gamma",
        "fit_scaling",
        "fit_t_",
    }
    missing = sorted(required - set(fitted_adata.var.columns))
    if missing:
        raise ValueError(f"Recovered Dynamics Fit Ranking is missing canonical fit fields: {missing}.")
    selected = fitted_adata.var["openbio_dynamics_selected"].to_numpy(dtype=bool)
    success = fitted_adata.var["openbio_dynamics_fit_success"].to_numpy(dtype=bool)
    velocity_gene = fitted_adata.var["openbio_dynamics_velocity_gene"].to_numpy(dtype=bool)
    if bool((success & ~selected).any()):
        raise ValueError("Recovered Dynamics fit-success mask includes unselected genes.")
    likelihood = fitted_adata.var["fit_likelihood"].to_numpy(dtype=float)
    if bool(numpy.isinf(likelihood).any()) or not bool(numpy.isfinite(likelihood[success]).all()):
        raise ValueError("Recovered Dynamics fit likelihood is invalid for successful genes.")
    original_position = numpy.arange(int(fitted_adata.n_vars), dtype=int)
    successful_positions = numpy.flatnonzero(success)
    successful_positions = successful_positions[
        numpy.lexsort((original_position[successful_positions], -likelihood[successful_positions]))
    ]
    returned_success = successful_positions[:top_n]
    failed_positions = numpy.flatnonzero(selected & ~success) if include_failed else numpy.asarray([], dtype=int)
    returned_positions = numpy.concatenate([returned_success, failed_positions])
    if returned_positions.size > max_output_rows:
        raise ValueError(
            f"Recovered Dynamics Fit Ranking would return {returned_positions.size:,} rows, "
            f"exceeding max_output_rows={max_output_rows:,}."
        )
    rows = []
    for rank_index, position in enumerate(returned_success, start=1):
        rows.append(
            {
                "rank": rank_index,
                "gene": str(fitted_adata.var_names[position]),
                "fit_status": "fitted",
                "fit_likelihood": float(likelihood[position]),
                "fit_r2": (
                    float(fitted_adata.var["fit_r2"].iloc[position])
                    if "fit_r2" in fitted_adata.var
                    and numpy.isfinite(float(fitted_adata.var["fit_r2"].iloc[position]))
                    else None
                ),
                "velocity_gene": bool(velocity_gene[position]),
                "fit_alpha": float(fitted_adata.var["fit_alpha"].iloc[position]),
                "fit_beta": float(fitted_adata.var["fit_beta"].iloc[position]),
                "fit_gamma": float(fitted_adata.var["fit_gamma"].iloc[position]),
                "fit_scaling": float(fitted_adata.var["fit_scaling"].iloc[position]),
                "fit_switch_time": float(fitted_adata.var["fit_t_"].iloc[position]),
            }
        )
    for position in failed_positions:
        rows.append(
            {
                "rank": None,
                "gene": str(fitted_adata.var_names[position]),
                "fit_status": "failed",
                "fit_likelihood": None,
                "fit_r2": None,
                "velocity_gene": bool(velocity_gene[position]),
                "fit_alpha": None,
                "fit_beta": None,
                "fit_gamma": None,
                "fit_scaling": None,
                "fit_switch_time": None,
            }
        )
    columns = [
        "rank",
        "gene",
        "fit_status",
        "fit_likelihood",
        "fit_r2",
        "velocity_gene",
        "fit_alpha",
        "fit_beta",
        "fit_gamma",
        "fit_scaling",
        "fit_switch_time",
    ]
    table = pandas.DataFrame(rows, columns=columns)
    table["rank"] = pandas.array(table["rank"], dtype="Int64")
    table["velocity_gene"] = table["velocity_gene"].astype(bool)
    parameters = {
        "top_n": top_n,
        "include_failed": include_failed,
        "max_output_rows": max_output_rows,
    }
    summary = _vs_summary(
        node_id="OpenBioSingleCellVelocityGeneRanking",
        status="descriptive_recovered_dynamics_fit_ranking",
        methods=(
            "Selected the canonical verified recover_dynamics fit bundle, ordered successful genes by descending "
            "fit_likelihood with original feature order as the stable tie-break, and optionally appended failed fits without ranks."
        ),
        results=(
            f"Returned {len(returned_success):,} ranked successful fit(s)"
            + (f" and {len(failed_positions):,} explicitly failed fit(s)." if include_failed else ".")
        ),
        key_results={
            "eligible_selected_genes": int(selected.sum()),
            "successful_fits": int(success.sum()),
            "failed_fits": int((selected & ~success).sum()),
            "returned_successful_fits": int(len(returned_success)),
            "returned_failed_fits": int(len(failed_positions)),
            "truncated": bool(len(successful_positions) > top_n),
            "fit_likelihood_distribution": _vs_numeric_summary(likelihood[selected]),
            "fit_r2_distribution": (
                _vs_numeric_summary(fitted_adata.var["fit_r2"].to_numpy(dtype=float)[selected])
                if "fit_r2" in fitted_adata.var
                else None
            ),
            "dynamics_fit_fingerprint_sha256": state["dynamics_fit_fingerprint_sha256"],
            "upstream_state_fingerprint_sha256": state["state_fingerprint_sha256"],
            "table_fingerprint_sha256": _vs_series_fingerprint(table, columns, label="dynamics-ranking"),
        },
        parameters=parameters,
        references=DYNAMICS_REFERENCES,
        software_versions=_vs_versions(str(state["scvelo_version"])),
        report_warnings=[],
        limitations=[
            "Fit likelihood ranks model fit quality; it is not a p-value, significance threshold, driver score, marker test, or Condition contrast.",
            "top_n is a reporting truncation and stable ties follow fitted feature order.",
            "Recovered kinetic parameters remain model-dependent rather than directly measured rates.",
        ],
    )
    return table, summary


def render_velocity_stream(
    graph_adata: Any,
    *,
    basis: str = "umap",
    color_key: str = "leiden",
    density: float = 2.0,
    smooth: float = 0.5,
    min_mass: float = 1.0,
    scvelo_module: Any | None = None,
) -> tuple[bytes, dict[str, Any]]:
    _, numpy, pandas, sparse = _vs_science()
    state = validate_portable_velocity_state(graph_adata, allowed_stages=("velocity_graph",))
    basis = _vs_clean_text(basis, label="basis")
    if basis.startswith("X_"):
        raise ValueError("Velocity stream basis must omit the 'X_' obsm prefix.")
    color_key = _vs_clean_text(color_key, label="color_key", allow_empty=True)
    density = _vs_float(density, label="density", minimum=0.0, strictly_greater=True)
    smooth = _vs_float(smooth, label="smooth", minimum=0.0, strictly_greater=True)
    min_mass = _vs_float(min_mass, label="min_mass", minimum=0.0)
    embedding_key = f"X_{basis}"
    if embedding_key not in graph_adata.obsm:
        raise ValueError(f"Velocity stream embedding not found in obsm[{embedding_key!r}].")
    embedding = graph_adata.obsm[embedding_key]
    if sparse.issparse(embedding):
        raise TypeError("Velocity stream embedding must be a dense numeric n_obs × 2 matrix.")
    embedding_values = numpy.asarray(embedding)
    if embedding_values.shape != (int(graph_adata.n_obs), 2):
        raise ValueError("Velocity stream embedding must have exact shape (n_obs, 2).")
    if not numpy.issubdtype(embedding_values.dtype, numpy.number) or bool(
        numpy.iscomplexobj(embedding_values)
    ):
        raise TypeError("Velocity stream embedding must be real numeric.")
    if not bool(numpy.isfinite(embedding_values).all()):
        raise ValueError("Velocity stream embedding contains non-finite coordinates.")
    color_semantics: dict[str, Any]
    if color_key:
        if color_key not in graph_adata.obs:
            raise ValueError(f"Velocity stream color column not found: {color_key!r}.")
        color = graph_adata.obs[color_key]
        if bool(color.isna().any()):
            raise ValueError("Velocity stream color column contains missing values.")
        if isinstance(color.dtype, pandas.CategoricalDtype):
            observed = list(dict.fromkeys(str(value) for value in color.astype(object)))
            color_semantics = {
                "kind": "categorical",
                "observed_categories": observed,
                "declared_categories": [str(value) for value in color.cat.categories],
            }
        elif pandas.api.types.is_numeric_dtype(color.dtype) and not pandas.api.types.is_bool_dtype(
            color.dtype
        ):
            numeric_color = color.to_numpy(dtype=float)
            if not bool(numpy.isfinite(numeric_color).all()):
                raise ValueError("Velocity stream numeric color column contains non-finite values.")
            color_semantics = {"kind": "numeric", "distribution": _vs_numeric_summary(numeric_color)}
        else:
            raise TypeError("Velocity stream color column must be categorical or real numeric.")
    else:
        color_semantics = {"kind": "none"}
    vkey = str(state["vkey"])
    mask = graph_adata.var[f"{vkey}_genes"].to_numpy(dtype=bool)
    velocity_layer = graph_adata.layers[vkey]
    selected_velocity = (
        velocity_layer[:, mask].toarray()
        if sparse.issparse(velocity_layer)
        else numpy.asarray(velocity_layer)[:, mask]
    )
    selected_velocity = numpy.asarray(selected_velocity, dtype=float)
    if not bool(numpy.isfinite(selected_velocity).all()):
        raise ValueError("Velocity stream requires finite selected-gene velocity vectors.")
    usable = numpy.linalg.norm(selected_velocity, axis=1) > 0
    work = graph_adata.copy()
    scvelo_module, scvelo_version, stream = _vs_backend(scvelo_module, "velocity_embedding_stream")
    import matplotlib.pyplot as pyplot
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    existing_figures = set(pyplot.get_fignums())
    figure = Figure(figsize=(8.0, 6.0), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    try:
        with _vs_preserve_global_state(scvelo_module), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            stream(
                work,
                basis=basis,
                vkey=vkey,
                density=density,
                smooth=smooth,
                min_mass=min_mass,
                color=color_key or None,
                show=False,
                save=None,
                ax=axis,
            )
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
        png = buffer.getvalue()
    except Exception as exc:
        raise RuntimeError(f"scVelo 0.3.4 velocity stream rendering failed ({type(exc).__name__}: {exc}).") from exc
    finally:
        for number in set(pyplot.get_fignums()) - existing_figures:
            pyplot.close(number)
        figure.clear()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError("Velocity stream renderer did not produce valid PNG bytes.")
    parameters = {
        "basis": basis,
        "color_key": color_key,
        "density": density,
        "smooth": smooth,
        "min_mass": min_mass,
        "vkey": vkey,
        "figure_size_inches": [8.0, 6.0],
        "dpi": 120,
    }
    summary = _vs_summary(
        node_id="OpenBioSingleCellVelocityStreamPlot",
        status="velocity_stream_visualization_rendered",
        methods=(
            "Rendered scVelo 0.3.4's smoothed embedding stream field from a private copy of an exact verified "
            "velocity-graph state, on one finite two-dimensional embedding with explicit density, smoothing, and mass policy."
        ),
        results=(
            f"Rendered a PNG stream-field visualization for {graph_adata.n_obs:,} cells on {embedding_key!r}; "
            "the streamlines are interpolation and not cell paths, lineage, time, transition probability, or fate."
        ),
        key_results={
            "input_cells": int(graph_adata.n_obs),
            "selected_velocity_genes": int(mask.sum()),
            "usable_nonzero_velocity_cells": int(usable.sum()),
            "usable_nonzero_velocity_fraction": float(usable.mean()),
            "embedding_key": embedding_key,
            "embedding_axis_1_distribution": _vs_numeric_summary(embedding_values[:, 0]),
            "embedding_axis_2_distribution": _vs_numeric_summary(embedding_values[:, 1]),
            "embedding_fingerprint_sha256": _vs_matrix_fingerprint(
                embedding_values, label=embedding_key
            ),
            "color_semantics": color_semantics,
            "velocity_mode": state["velocity_mode"],
            "velocity_graph_positive_fingerprint_sha256": state[
                "velocity_graph_positive_fingerprint_sha256"
            ],
            "velocity_graph_negative_fingerprint_sha256": state[
                "velocity_graph_negative_fingerprint_sha256"
            ],
            "upstream_state_fingerprint_sha256": state["state_fingerprint_sha256"],
            "private_embedding_cache_created": f"{vkey}_{basis}" in work.obsm,
            "png_sha256": hashlib.sha256(png).hexdigest(),
            "png_bytes": len(png),
        },
        parameters=parameters,
        references=VELOCITY_REFERENCES,
        software_versions=_vs_versions(scvelo_version),
        report_warnings=[f"scVelo warning: {item.message}" for item in caught],
        limitations=[
            "Streamlines are smoothed grid interpolation influenced by embedding geometry, density, smoothness, and mass thresholds.",
            "The image does not show individual-cell paths, branch confidence, transition probabilities, measured time, lineage, or fate.",
            "Biological direction must not be inferred from the picture without independent kinetic and experimental evidence.",
        ],
    )
    return png, summary


__all__ = [
    "VELOCITY_ARTIFACT_TYPE",
    "VELOCITY_PORTABLE_SCHEMA",
    "VELOCITY_SCVELO_VERSION",
    "VELOCITY_STAGES",
    "VELOCITY_STATE_KEY",
    "VELOCITY_SUMMARY_SCHEMA",
    "build_velocity_graph",
    "compute_velocity_moments",
    "estimate_rna_velocity",
    "prepare_velocity_abundances",
    "rank_recovered_dynamics",
    "recover_velocity_dynamics",
    "render_velocity_stream",
    "validate_portable_velocity_state",
]
