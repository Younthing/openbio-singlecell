from __future__ import annotations

import hashlib
import importlib
import inspect
import json
import math
import os
import re
import stat
import tempfile
import threading
import warnings
import weakref
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass, fields, is_dataclass
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any

ad: Any = None
np: Any = None
pd: Any = None
sparse: Any = None
KMeans: Any = None
euclidean_distances: Any = None


def _load_science() -> None:
    global KMeans, ad, euclidean_distances, np, pd, sparse
    if ad is not None:
        return
    try:
        import anndata as _ad
        import numpy as _np
        import pandas as _pd
        from scipy import sparse as _sparse
        from sklearn.cluster import KMeans as _KMeans
        from sklearn.metrics import euclidean_distances as _euclidean_distances
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            "cNMF scientific dependencies are unavailable; install the repository scientific requirements."
        ) from exc
    ad = _ad
    np = _np
    pd = _pd
    sparse = _sparse
    KMeans = _KMeans
    euclidean_distances = _euclidean_distances


CNMF_REQUIRED_VERSION = "1.7.1"
CNMF_ADAPTER_VERSION = 4
CNMF_LOCAL_NEIGHBORHOOD_SIZE = 0.30
CNMF_STATISTICS_DENSITY_THRESHOLD = 2.0
CNMF_RESOURCE_BUDGET_BYTES = 2 * 1024**3
CNMF_RESOURCE_SAFETY_FACTOR = 2.0
CNMF_AUDITED_PYPI_WHEEL_SHA256 = "a4b237169626f632fee04a31441365669b12a02c72a93b079686b0ce47a29ec9"
_RUN_NAME = "run"
_PREPARE_LOCK = threading.Lock()
_RLOCK_TYPE = type(threading.RLock())
_RUN_ARTIFACT_MARKER = ("OPENBIO_CNMF_RUN", 1, CNMF_ADAPTER_VERSION, CNMF_REQUIRED_VERSION)
_METRIC_RECORD_MARKER = ("OPENBIO_CNMF_METRIC", 1)
_POLICY_RECORD_MARKER = ("OPENBIO_CNMF_POLICY", 1)
_RESOURCE_RECORD_MARKER = ("OPENBIO_CNMF_RESOURCE", 1)
_METADATA_RECORD_MARKER = ("OPENBIO_CNMF_METADATA", 1)


@dataclass(frozen=True, slots=True)
class CNMFKMetric:
    _OPENBIO_RECORD_MARKER = _METRIC_RECORD_MARKER

    k: int
    stability: float
    prediction_error: float
    statistics_density_threshold: float
    expected_restarts: int
    completed_restarts: int
    combined_components: int

    def as_dict(self) -> dict[str, int | float]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class CNMFFixedPolicy:
    _OPENBIO_RECORD_MARKER = _POLICY_RECORD_MARKER

    densify: bool
    tpm_fn: str | None
    beta_loss: str
    genes_file: str | None
    alpha_usage: float
    alpha_spectra: float
    init: str
    max_nmf_iter: int
    worker_i: int
    total_workers: int
    factorize_skip_completed_runs: bool
    combine_skip_missing_files: bool
    statistics_density_threshold: float
    statistics_local_neighborhood_size: float
    statistics_show_clustering: bool
    statistics_skip_density_and_return_after_stats: bool
    statistics_close_clustergram_fig: bool
    statistics_refit_usage: bool
    statistics_normalize_tpm_spectra: bool
    statistics_build_ref: bool
    storage_mode: str


@dataclass(frozen=True, slots=True)
class CNMFResourceEstimate:
    _OPENBIO_RECORD_MARKER = _RESOURCE_RECORD_MARKER

    observation_count: int
    highvar_gene_count: int
    total_restarts: int
    iter_merged_spectra_bytes: int
    density_distance_work_bytes: int
    prediction_error_work_bytes: int
    statistics_stage_base_bytes: int
    safety_factor: float
    estimated_peak_bytes: int
    budget_bytes: int

    def as_dict(self) -> dict[str, int | float]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class CNMFRunMetadata:
    _OPENBIO_RECORD_MARKER = _METADATA_RECORD_MARKER

    adapter_version: int
    backend_name: str
    backend_version: str
    audited_pypi_wheel_sha256: str
    installed_distribution_attested: bool
    source_kind: str
    source_layer: str | None
    input_advisories: tuple[str, ...]
    input_cells: int
    input_genes: int
    current_features: int
    input_total_counts: float
    input_nonzero_entries: int
    input_fingerprint: str
    parameters_fingerprint: str
    candidate_ks: tuple[int, ...]
    n_iter: int
    num_highvar_genes: int
    realized_highvar_genes: tuple[str, ...]
    positive_variance_highvar_genes: int
    random_seed: int
    restart_seeds: tuple[tuple[int, int, int], ...]
    total_restarts: int
    requested_resource: CNMFResourceEstimate
    realized_resource: CNMFResourceEstimate
    fixed_policy: CNMFFixedPolicy
    backend_paths_fingerprint: str
    prepare_manifest_sha256: str
    factorization_manifest_sha256: str
    artifact_manifest_sha256: str
    artifact_hashes: tuple[tuple[str, str], ...]
    known_upstream_warnings: tuple[tuple[str, int], ...]


def _cleanup_temporary_directory(directory: tempfile.TemporaryDirectory[str]) -> None:
    directory.cleanup()


def _add_cleanup_note(primary: BaseException, cleanup: BaseException, *, context: str) -> None:
    primary.add_note(f"{context} cleanup also failed with {type(cleanup).__name__}: {cleanup}")


def _cleanup_preserving_primary(cleanup: Any, primary: BaseException, *, context: str) -> None:
    try:
        cleanup()
    except BaseException as cleanup_error:
        _add_cleanup_note(primary, cleanup_error, context=context)


def close_run_preserving_primary(run: Any, primary: BaseException, *, context: str) -> None:
    """Close an owned run without replacing an already-active primary exception."""

    _cleanup_preserving_primary(run.close, primary, context=context)


def _require_portable_record(value: Any, expected_type: type[Any], marker: tuple[str, int], *, label: str) -> None:
    value_type = type(value)
    if (
        vars(value_type).get("_OPENBIO_RECORD_MARKER") != marker
        or not is_dataclass(value)
        or not bool(getattr(getattr(value_type, "__dataclass_params__", None), "frozen", False))
        or hasattr(value, "__dict__")
        or tuple(field.name for field in fields(value)) != tuple(field.name for field in fields(expected_type))
    ):
        raise TypeError(f"OPENBIO_CNMF_RUN {label} is not the exact portable frozen record schema.")


def _is_plain_int(value: Any, *, minimum: int | None = None) -> bool:
    return type(value) is int and (minimum is None or value >= minimum)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"sha256:[0-9a-f]{64}", value) is not None


def _validate_run_records(metadata: Any, metrics: Any) -> None:
    _require_portable_record(metadata, CNMFRunMetadata, _METADATA_RECORD_MARKER, label="metadata")
    if not isinstance(metrics, tuple) or not metrics:
        raise TypeError("OPENBIO_CNMF_RUN metrics must be a nonempty tuple.")
    for metric in metrics:
        _require_portable_record(metric, CNMFKMetric, _METRIC_RECORD_MARKER, label="metric")
    _require_portable_record(
        metadata.requested_resource,
        CNMFResourceEstimate,
        _RESOURCE_RECORD_MARKER,
        label="requested resource",
    )
    _require_portable_record(
        metadata.realized_resource,
        CNMFResourceEstimate,
        _RESOURCE_RECORD_MARKER,
        label="realized resource",
    )
    _require_portable_record(metadata.fixed_policy, CNMFFixedPolicy, _POLICY_RECORD_MARKER, label="fixed policy")

    if metadata.adapter_version != CNMF_ADAPTER_VERSION:
        raise RuntimeError("OPENBIO_CNMF_RUN adapter version is incompatible.")
    if metadata.backend_name != "cnmf.cNMF" or metadata.backend_version != CNMF_REQUIRED_VERSION:
        raise RuntimeError("OPENBIO_CNMF_RUN metadata has an incompatible backend identity.")
    if metadata.audited_pypi_wheel_sha256 != CNMF_AUDITED_PYPI_WHEEL_SHA256:
        raise RuntimeError("OPENBIO_CNMF_RUN audited release evidence differs from this adapter.")
    if metadata.installed_distribution_attested is not False:
        raise RuntimeError("OPENBIO_CNMF_RUN must not claim installed-distribution attestation.")
    if metadata.source_kind not in {"X", "raw", "layer"}:
        raise RuntimeError("OPENBIO_CNMF_RUN source kind is invalid.")
    if metadata.source_kind in {"X", "raw"} and metadata.source_layer is not None:
        raise RuntimeError("OPENBIO_CNMF_RUN X/Raw source cannot name a layer.")
    if metadata.source_kind == "layer" and not isinstance(metadata.source_layer, str):
        raise RuntimeError("OPENBIO_CNMF_RUN layer source must name its layer.")
    if not isinstance(metadata.input_advisories, tuple) or any(
        not isinstance(message, str) or not message for message in metadata.input_advisories
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN input advisories must be a tuple of nonempty messages.")
    if not _is_plain_int(metadata.input_cells, minimum=1) or not _is_plain_int(metadata.input_genes, minimum=1):
        raise RuntimeError("OPENBIO_CNMF_RUN input dimensions are invalid.")
    if not _is_plain_int(metadata.current_features, minimum=0):
        raise RuntimeError("OPENBIO_CNMF_RUN Survey-time current feature count is invalid.")
    if metadata.source_kind != "raw" and metadata.current_features != metadata.input_genes:
        raise RuntimeError("OPENBIO_CNMF_RUN X/layer source features must match the current feature count.")
    if (
        type(metadata.input_total_counts) is not float
        or not math.isfinite(metadata.input_total_counts)
        or metadata.input_total_counts <= 0
        or not _is_plain_int(metadata.input_nonzero_entries, minimum=1)
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN input count accounting is invalid.")
    for label, value in (
        ("input", metadata.input_fingerprint),
        ("parameters", metadata.parameters_fingerprint),
        ("backend paths", metadata.backend_paths_fingerprint),
        ("prepare manifest", metadata.prepare_manifest_sha256),
        ("factorization manifest", metadata.factorization_manifest_sha256),
        ("artifact manifest", metadata.artifact_manifest_sha256),
    ):
        if not _is_sha256(value):
            raise RuntimeError(f"OPENBIO_CNMF_RUN {label} fingerprint is invalid.")
    if (
        not isinstance(metadata.candidate_ks, tuple)
        or not metadata.candidate_ks
        or any(not _is_plain_int(k, minimum=2) for k in metadata.candidate_ks)
        or metadata.candidate_ks != tuple(range(metadata.candidate_ks[0], metadata.candidate_ks[-1] + 1))
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN candidate K family is invalid.")
    if not _is_plain_int(metadata.n_iter, minimum=2) or not _is_plain_int(metadata.num_highvar_genes, minimum=1):
        raise RuntimeError("OPENBIO_CNMF_RUN restart/HVG request is invalid.")
    if (
        not isinstance(metadata.realized_highvar_genes, tuple)
        or not metadata.realized_highvar_genes
        or any(not isinstance(gene, str) or not gene for gene in metadata.realized_highvar_genes)
        or len(set(metadata.realized_highvar_genes)) != len(metadata.realized_highvar_genes)
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN realized HVG family is invalid.")
    if not _is_plain_int(metadata.positive_variance_highvar_genes, minimum=0) or (
        metadata.positive_variance_highvar_genes > len(metadata.realized_highvar_genes)
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN positive-variance HVG accounting is invalid.")
    if not _is_plain_int(metadata.random_seed, minimum=0) or metadata.random_seed > 2**31 - 1:
        raise RuntimeError("OPENBIO_CNMF_RUN random seed is invalid.")
    expected_total = len(metadata.candidate_ks) * metadata.n_iter
    if metadata.total_restarts != expected_total:
        raise RuntimeError("OPENBIO_CNMF_RUN total restart accounting is invalid.")
    if (
        not isinstance(metadata.restart_seeds, tuple)
        or len(metadata.restart_seeds) != expected_total
        or any(
            not isinstance(item, tuple) or len(item) != 3 or not all(_is_plain_int(value) for value in item)
            for item in metadata.restart_seeds
        )
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN restart-seed family is invalid.")
    expected_pairs = [(k, iteration) for k in metadata.candidate_ks for iteration in range(metadata.n_iter)]
    if [(k, iteration) for k, iteration, _ in metadata.restart_seeds] != expected_pairs:
        raise RuntimeError("OPENBIO_CNMF_RUN restart-seed ordering is invalid.")
    if asdict(metadata.fixed_policy) != asdict(_fixed_policy()):
        raise RuntimeError("OPENBIO_CNMF_RUN fixed policy differs from this adapter.")
    if asdict(metadata.requested_resource) != asdict(
        _resource_estimate(
            metadata.candidate_ks,
            metadata.n_iter,
            metadata.num_highvar_genes,
            metadata.input_cells,
        )
    ) or asdict(metadata.realized_resource) != asdict(
        _resource_estimate(
            metadata.candidate_ks,
            metadata.n_iter,
            len(metadata.realized_highvar_genes),
            metadata.input_cells,
        )
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN resource accounting differs from its dimensions.")
    if (
        not isinstance(metadata.artifact_hashes, tuple)
        or not metadata.artifact_hashes
        or any(
            not isinstance(item, tuple)
            or len(item) != 2
            or not isinstance(item[0], str)
            or not item[0]
            or Path(item[0]).is_absolute()
            or ".." in Path(item[0]).parts
            or not _is_sha256(item[1])
            for item in metadata.artifact_hashes
        )
        or len({path for path, _ in metadata.artifact_hashes}) != len(metadata.artifact_hashes)
        or metadata.artifact_manifest_sha256 != _manifest(metadata.artifact_hashes)
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN artifact-hash family is invalid.")
    if not isinstance(metadata.known_upstream_warnings, tuple) or any(
        not isinstance(item, tuple)
        or len(item) != 2
        or not isinstance(item[0], str)
        or not _is_plain_int(item[1], minimum=1)
        for item in metadata.known_upstream_warnings
    ):
        raise RuntimeError("OPENBIO_CNMF_RUN known-warning accounting is invalid.")

    if tuple(metric.k for metric in metrics) != metadata.candidate_ks:
        raise RuntimeError("OPENBIO_CNMF_RUN metrics do not match the candidate K family.")
    for metric in metrics:
        if (
            type(metric.stability) is not float
            or not math.isfinite(metric.stability)
            or not -1.0 <= metric.stability <= 1.0
            or type(metric.prediction_error) is not float
            or not math.isfinite(metric.prediction_error)
            or metric.prediction_error < 0
            or metric.statistics_density_threshold != CNMF_STATISTICS_DENSITY_THRESHOLD
            or metric.expected_restarts != metadata.n_iter
            or metric.completed_restarts != metadata.n_iter
            or metric.combined_components != metric.k * metadata.n_iter
        ):
            raise RuntimeError(f"OPENBIO_CNMF_RUN K={metric.k} metric contents are invalid.")


def _run_snapshot_fingerprint(metadata: Any, metrics: Any) -> str:
    return _json_fingerprint(
        {
            "artifact_marker": list(_RUN_ARTIFACT_MARKER),
            "metadata": asdict(metadata),
            "metrics": [asdict(metric) for metric in metrics],
        }
    )


class CNMFRun:
    """Owned, process-local, file-backed cNMF state for staged analysis."""

    _OPENBIO_PORTABLE_ARTIFACT = _RUN_ARTIFACT_MARKER

    __slots__ = (
        "_metadata",
        "_metrics",
        "_artifact_snapshot_sha256",
        "_backend",
        "_base_adata",
        "_temporary_directory",
        "_root",
        "_root_identity",
        "_lock",
        "_closed",
        "_finalizer",
        "__weakref__",
    )

    def __init__(
        self,
        *,
        metadata: CNMFRunMetadata,
        metrics: Sequence[CNMFKMetric],
        backend: Any,
        base_adata: Any,
        temporary_directory: tempfile.TemporaryDirectory[str],
    ) -> None:
        metrics_tuple = tuple(metrics)
        if tuple(metric.k for metric in metrics_tuple) != metadata.candidate_ks:
            raise ValueError("CNMFRun metrics must match the ordered candidate K values.")
        if backend is None or base_adata is None:
            raise ValueError("CNMFRun requires backend state and a private AnnData copy.")
        root = Path(temporary_directory.name).absolute()
        root_stat = _require_directory_identity(root, label="managed root")
        _validate_run_records(metadata, metrics_tuple)
        self._metadata = metadata
        self._metrics = metrics_tuple
        self._artifact_snapshot_sha256 = _run_snapshot_fingerprint(metadata, metrics_tuple)
        self._backend = backend
        self._base_adata = base_adata
        self._temporary_directory = temporary_directory
        self._root = root
        self._root_identity = _identity(root_stat)
        self._lock = threading.RLock()
        self._closed = False
        self._finalizer = weakref.finalize(self, _cleanup_temporary_directory, temporary_directory)

    def _assert_snapshot(self) -> None:
        _validate_run_records(self._metadata, self._metrics)
        current = _run_snapshot_fingerprint(self._metadata, self._metrics)
        if current != self._artifact_snapshot_sha256:
            raise RuntimeError("OPENBIO_CNMF_RUN metadata or metrics changed after Survey construction.")

    @property
    def metadata(self) -> CNMFRunMetadata:
        self._assert_snapshot()
        return self._metadata

    @property
    def metrics(self) -> tuple[CNMFKMetric, ...]:
        self._assert_snapshot()
        return self._metrics

    @property
    def candidate_ks(self) -> tuple[int, ...]:
        self._assert_snapshot()
        return self._metadata.candidate_ks

    @property
    def closed(self) -> bool:
        return self._closed

    @property
    def private_root(self) -> str:
        """Diagnostic-only path; callers must never reconnect a run from this string."""

        with self._lock:
            self._assert_live()
            return str(self._root)

    def _assert_live(self) -> None:
        self._assert_snapshot()
        if self._closed or self._backend is None or self._base_adata is None:
            raise RuntimeError("OPENBIO_CNMF_RUN is closed and can no longer be used.")
        current = _require_directory_identity(self._root, label="managed root")
        if _identity(current) != self._root_identity:
            raise RuntimeError("OPENBIO_CNMF_RUN managed root was replaced after creation.")

    def copy_base_adata(self) -> Any:
        with self._lock:
            self._assert_live()
            return self._base_adata.copy()

    @contextmanager
    def locked_backend(self) -> Iterator[Any]:
        with self._lock:
            self._assert_live()
            yield self._backend
            self._assert_live()

    def close(self) -> None:
        with self._lock:
            if self._closed and self._temporary_directory is None:
                return
            if not self._closed:
                self._closed = True
                self._backend = None
                self._base_adata = None
            directory = self._temporary_directory
            try:
                if directory is not None:
                    directory.cleanup()
            except BaseException:
                # Keep the owner/finalizer reachable so close() or finalization can retry.
                raise
            else:
                self._temporary_directory = None
                self._finalizer.detach()

    def __enter__(self) -> CNMFRun:
        self._assert_live()
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, traceback
        if exc is None:
            self.close()
            return
        close_run_preserving_primary(self, exc, context="OPENBIO_CNMF_RUN context-manager")

    def __reduce_ex__(self, protocol: int) -> Any:
        del protocol
        raise TypeError("OPENBIO_CNMF_RUN is process-local and cannot be serialized.")


def _validate_run_artifact(run: Any) -> None:
    run_type = type(run)
    if run_type is not CNMFRun:
        if vars(run_type).get("_OPENBIO_PORTABLE_ARTIFACT") != _RUN_ARTIFACT_MARKER:
            raise TypeError(
                "cNMF Consensus Programs requires a concrete or explicitly portable OPENBIO_CNMF_RUN artifact."
            )
        if vars(run_type).get("__slots__") != CNMFRun.__slots__ or hasattr(run, "__dict__"):
            raise TypeError("Portable OPENBIO_CNMF_RUN has an incompatible owned-state schema.")
        for name in ("metadata", "metrics", "candidate_ks", "closed", "private_root"):
            descriptor = vars(run_type).get(name)
            if not isinstance(descriptor, property) or descriptor.fset is not None:
                raise TypeError(f"Portable OPENBIO_CNMF_RUN {name} must be a getter-only property.")
        for name in ("_assert_snapshot", "_assert_live", "copy_base_adata", "locked_backend", "close"):
            if not callable(vars(run_type).get(name)):
                raise TypeError(f"Portable OPENBIO_CNMF_RUN is missing its class-owned {name} implementation.")
    try:
        metadata = object.__getattribute__(run, "_metadata")
        metrics = object.__getattribute__(run, "_metrics")
        snapshot = object.__getattribute__(run, "_artifact_snapshot_sha256")
        root = object.__getattribute__(run, "_root")
        root_identity = object.__getattribute__(run, "_root_identity")
        lock = object.__getattribute__(run, "_lock")
        closed = object.__getattribute__(run, "_closed")
        backend = object.__getattribute__(run, "_backend")
        base_adata = object.__getattribute__(run, "_base_adata")
        temporary = object.__getattribute__(run, "_temporary_directory")
        finalizer = object.__getattribute__(run, "_finalizer")
    except (AttributeError, TypeError) as exc:
        raise TypeError("OPENBIO_CNMF_RUN owned state is incomplete.") from exc
    _validate_run_records(metadata, metrics)
    if not _is_sha256(snapshot) or snapshot != _run_snapshot_fingerprint(metadata, metrics):
        raise RuntimeError("OPENBIO_CNMF_RUN metadata or metrics differ from their Survey construction snapshot.")
    if not isinstance(root, Path) or not (
        isinstance(root_identity, tuple)
        and len(root_identity) == 2
        and all(_is_plain_int(value) for value in root_identity)
    ):
        raise TypeError("OPENBIO_CNMF_RUN managed-root identity state is invalid.")
    if type(lock) is not _RLOCK_TYPE:
        raise TypeError("OPENBIO_CNMF_RUN lock owner is invalid.")
    if type(closed) is not bool:
        raise TypeError("OPENBIO_CNMF_RUN lifecycle state is invalid.")
    if closed or backend is None or base_adata is None or temporary is None:
        raise RuntimeError("OPENBIO_CNMF_RUN is closed and can no longer be used.")
    if not isinstance(temporary, tempfile.TemporaryDirectory) or not isinstance(finalizer, weakref.finalize):
        raise TypeError("OPENBIO_CNMF_RUN temporary-directory owner or finalizer is invalid.")
    current_root = _require_directory_identity(root, label="managed root")
    if _identity(current_root) != root_identity:
        raise RuntimeError("OPENBIO_CNMF_RUN managed root was replaced after creation.")


@dataclass(frozen=True, slots=True)
class _BackendAPI:
    module: Any
    cnmf_class: Any
    load_frame: Any
    version: str


def _identity(value: os.stat_result) -> tuple[int, int]:
    return int(value.st_dev), int(value.st_ino)


def _is_reparse(value: os.stat_result) -> bool:
    attributes = int(getattr(value, "st_file_attributes", 0))
    reparse_flag = int(getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400))
    return stat.S_ISLNK(value.st_mode) or bool(attributes & reparse_flag)


def _require_directory_identity(path: Path, *, label: str) -> os.stat_result:
    try:
        value = path.lstat()
    except FileNotFoundError as exc:
        raise RuntimeError(f"cNMF {label} is missing: {path}") from exc
    if _is_reparse(value) or not stat.S_ISDIR(value.st_mode):
        raise RuntimeError(f"cNMF {label} must be a real directory, not a link or reparse point: {path}")
    return value


def _integer(value: Any, *, name: str, minimum: int, maximum: int | None = None) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"cNMF {name} must be an integer.")
    try:
        converted = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"cNMF {name} must be an integer.") from exc
    try:
        equal = bool(converted == value)
    except Exception:
        equal = False
    if not equal:
        raise TypeError(f"cNMF {name} must be an integer.")
    if converted < minimum or (maximum is not None and converted > maximum):
        upper = f" and at most {maximum}" if maximum is not None else ""
        raise ValueError(f"cNMF {name} must be at least {minimum}{upper}.")
    return converted


def _finite_float(value: Any, *, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"cNMF {name} must be a number.")
    try:
        converted = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"cNMF {name} must be a number.") from exc
    if not np.isfinite(converted) or not minimum <= converted <= maximum:
        raise ValueError(f"cNMF {name} must be finite and between {minimum} and {maximum}.")
    return converted


def _require_parameters(callable_object: Any, names: Sequence[str], *, label: str, version: str) -> None:
    try:
        parameters = inspect.signature(callable_object).parameters
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"Cannot inspect cnmf {version} {label}.") from exc
    missing = [name for name in names if name not in parameters]
    if missing:
        raise RuntimeError(
            f"Installed cnmf {version} has an incompatible {label}; missing explicit parameters: {missing}."
        )


def _load_backend() -> _BackendAPI:
    try:
        distribution_version = importlib_metadata.version("cnmf")
    except importlib_metadata.PackageNotFoundError as exc:
        raise RuntimeError("The optional package 'cnmf==1.7.1' is required for cNMF nodes.") from exc
    if distribution_version != CNMF_REQUIRED_VERSION:
        raise RuntimeError(
            f"cNMF nodes require exact cnmf=={CNMF_REQUIRED_VERSION}; installed distribution is {distribution_version}."
        )
    try:
        module = importlib.import_module("cnmf")
    except (ImportError, OSError) as exc:
        raise RuntimeError(f"Cannot import cnmf=={CNMF_REQUIRED_VERSION}: {exc}") from exc
    module_version = str(getattr(module, "__version__", "unknown"))
    if module_version != CNMF_REQUIRED_VERSION:
        raise RuntimeError(
            f"cnmf module version {module_version!r} disagrees with required distribution {CNMF_REQUIRED_VERSION}."
        )
    cnmf_class = getattr(module, "cNMF", None)
    load_frame = getattr(module, "load_df_from_npz", None)
    if cnmf_class is None or load_frame is None:
        raise RuntimeError("cnmf==1.7.1 does not expose cNMF and load_df_from_npz.")
    _require_parameters(cnmf_class, ("output_dir", "name"), label="cNMF constructor", version=module_version)
    for method, names in {
        "prepare": (
            "counts_fn",
            "components",
            "n_iter",
            "densify",
            "tpm_fn",
            "seed",
            "beta_loss",
            "num_highvar_genes",
            "genes_file",
            "alpha_usage",
            "alpha_spectra",
            "init",
            "max_NMF_iter",
        ),
        "factorize": ("worker_i", "total_workers", "skip_completed_runs"),
        "combine": ("components", "skip_missing_files"),
        "consensus": (
            "k",
            "density_threshold",
            "local_neighborhood_size",
            "show_clustering",
            "build_ref",
            "skip_density_and_return_after_stats",
            "close_clustergram_fig",
            "refit_usage",
            "normalize_tpm_spectra",
            "norm_counts",
        ),
        "load_results": ("K", "density_threshold", "n_top_genes", "norm_usage"),
    }.items():
        target = getattr(cnmf_class, method, None)
        if target is None:
            raise RuntimeError(f"cnmf==1.7.1 is missing cNMF.{method}().")
        _require_parameters(target, names, label=f"cNMF.{method}", version=module_version)
    _require_parameters(load_frame, ("filename",), label="load_df_from_npz", version=module_version)
    return _BackendAPI(module=module, cnmf_class=cnmf_class, load_frame=load_frame, version=module_version)


def _fixed_policy() -> CNMFFixedPolicy:
    return CNMFFixedPolicy(
        densify=False,
        tpm_fn=None,
        beta_loss="frobenius",
        genes_file=None,
        alpha_usage=0.0,
        alpha_spectra=0.0,
        init="random",
        max_nmf_iter=1000,
        worker_i=0,
        total_workers=1,
        factorize_skip_completed_runs=False,
        combine_skip_missing_files=False,
        statistics_density_threshold=CNMF_STATISTICS_DENSITY_THRESHOLD,
        statistics_local_neighborhood_size=CNMF_LOCAL_NEIGHBORHOOD_SIZE,
        statistics_show_clustering=False,
        statistics_skip_density_and_return_after_stats=True,
        statistics_close_clustergram_fig=True,
        statistics_refit_usage=True,
        statistics_normalize_tpm_spectra=False,
        statistics_build_ref=False,
        storage_mode="owned-private-temporary-directory",
    )


def _resource_estimate(
    candidate_ks: Sequence[int], n_iter: int, highvar_genes: int, n_obs: int
) -> CNMFResourceEstimate:
    total_restarts = len(candidate_ks) * n_iter
    spectra_bytes = 2 * sum(candidate_ks) * n_iter * highvar_genes * 8
    density_bytes = 2 * max((k * n_iter) ** 2 for k in candidate_ks) * 8
    prediction_bytes = 3 * n_obs * highvar_genes * 8
    base_bytes = spectra_bytes + density_bytes + prediction_bytes
    return CNMFResourceEstimate(
        observation_count=n_obs,
        highvar_gene_count=highvar_genes,
        total_restarts=total_restarts,
        iter_merged_spectra_bytes=spectra_bytes,
        density_distance_work_bytes=density_bytes,
        prediction_error_work_bytes=prediction_bytes,
        statistics_stage_base_bytes=base_bytes,
        safety_factor=CNMF_RESOURCE_SAFETY_FACTOR,
        estimated_peak_bytes=math.ceil(base_bytes * CNMF_RESOURCE_SAFETY_FACTOR),
        budget_bytes=CNMF_RESOURCE_BUDGET_BYTES,
    )


def _validate_resource_budget(
    candidate_ks: Sequence[int], n_iter: int, highvar_genes: int, n_obs: int, *, stage: str
) -> CNMFResourceEstimate:
    del stage
    return _resource_estimate(candidate_ks, n_iter, highvar_genes, n_obs)


def _require_adata(adata: Any, *, operation: str, require_current_features: bool = True) -> None:
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"{operation} requires an AnnData input.")
    if bool(adata.isbacked):
        raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
    if adata.n_obs == 0 or (require_current_features and adata.n_vars == 0):
        raise ValueError(f"{operation} requires at least one cell and one selected-source feature.")
    _validate_axis_identifiers(adata.obs_names, axis_name="obs_names", operation=operation)
    if require_current_features:
        _validate_axis_identifiers(adata.var_names, axis_name="var_names", operation=operation)


def _validate_axis_identifiers(identifiers: Any, *, axis_name: str, operation: str) -> None:
    if len(identifiers) == 0:
        raise ValueError(f"{operation} requires a nonempty {axis_name} axis.")
    if not bool(identifiers.is_unique):
        raise ValueError(f"{operation} requires unique {axis_name} for strict alignment.")
    invalid = [
        str(value)
        for value in identifiers
        if not str(value) or any(character in str(value) for character in ("\x00", "\t", "\r", "\n"))
    ]
    if invalid:
        raise ValueError(
            f"{operation} requires nonempty {axis_name} without NUL, tab, carriage-return, or newline "
            "characters because cnmf==1.7.1 persists tab-delimited result axes."
        )


def _parse_source(source: str) -> tuple[str, str | None]:
    if source == "X":
        return "X", None
    if source == "raw":
        return "raw", None
    if isinstance(source, str) and source.startswith("layer:"):
        layer = source.removeprefix("layer:").strip()
        if layer:
            return "layer", layer
    raise ValueError("cNMF source must be 'X', 'raw', or 'layer:<name>'.")


def _selected_matrix(adata: Any, source_kind: str, source_layer: str | None) -> Any:
    if source_kind == "X":
        return adata.X
    if source_kind == "raw":
        if adata.raw is None:
            raise ValueError("cNMF source 'raw' was selected, but adata.raw is unavailable.")
        return adata.raw.X
    if source_layer not in adata.layers:
        raise ValueError(f"cNMF count layer not found: {source_layer!r}.")
    return adata.layers[source_layer]


def _source_axes(adata: Any, source_kind: str) -> tuple[Any, Any]:
    if source_kind == "raw":
        if adata.raw is None:
            raise ValueError("cNMF source 'raw' was selected, but adata.raw is unavailable.")
        return adata.raw.obs_names, adata.raw.var_names
    return adata.obs_names, adata.var_names


def _materialize_source_adata(adata: Any, source_kind: str) -> Any:
    if source_kind == "raw":
        if adata.raw is None:
            raise ValueError("cNMF source 'raw' was selected, but adata.raw is unavailable.")
        # Raw may have a wider feature axis; materializing that axis is part of the requested result semantics.
        return adata.raw.to_adata()
    # The caller is a one-shot worker and transfers ownership, so no upstream-protection copy is needed.
    return adata


def _matrix_values(matrix: Any) -> np.ndarray:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    return np.asarray(values)


def _matrix_totals(matrix: Any, axis: int) -> np.ndarray:
    return np.asarray(matrix.sum(axis=axis), dtype=float).reshape(-1)


def _validate_counts(
    adata: Any, source_kind: str, source_layer: str | None
) -> tuple[Any, Any, Any, np.ndarray, np.ndarray, int, float, tuple[str, ...]]:
    matrix = _selected_matrix(adata, source_kind, source_layer)
    source_obs, source_vars = _source_axes(adata, source_kind)
    if tuple(matrix.shape) != (len(source_obs), len(source_vars)):
        raise ValueError("cNMF expression source must align exactly to its declared observation/feature axes.")
    _validate_axis_identifiers(source_obs, axis_name="source obs_names", operation="cNMF Rank Survey")
    _validate_axis_identifiers(source_vars, axis_name="source var_names", operation="cNMF Rank Survey")
    advisories: list[str] = []
    values = _matrix_values(matrix)
    if values.dtype.kind not in "iuf" or not bool(np.isfinite(values).all()):
        raise ValueError("cNMF count source must contain finite numeric values.")
    if bool((values < 0).any()):
        raise ValueError("cNMF count source contains negative values.")
    if values.size and not bool(np.allclose(values, np.rint(values), rtol=0.0, atol=1e-8)):
        advisories.append(
            "The selected non-negative source is not integer-like. The backend can calculate it, but the default "
            "cNMF UMI-count interpretation may not apply."
        )
    cell_totals = _matrix_totals(matrix, 1)
    gene_totals = _matrix_totals(matrix, 0)
    zero_cells = int((cell_totals <= 0).sum())
    zero_genes = int((gene_totals <= 0).sum())
    if zero_cells:
        raise ValueError(f"cNMF count source contains {zero_cells} cells with zero total counts.")
    if zero_genes:
        advisories.append(
            f"The selected source contains {zero_genes} zero-total genes. The official backend can exclude them "
            "from HVG factorization and retain them as zero full-gene results; filtering them upstream may reduce work."
        )
    if sparse.issparse(matrix):
        canonical = matrix.tocsr(copy=True)
        canonical.sum_duplicates()
        nonzero = int(canonical.count_nonzero())
    else:
        nonzero = int(np.count_nonzero(np.asarray(matrix)))
    return (
        matrix,
        source_obs,
        source_vars,
        cell_totals,
        gene_totals,
        nonzero,
        float(cell_totals.sum()),
        tuple(advisories),
    )


def _update_hash_names(digest: Any, names: Any) -> None:
    for value in names:
        encoded = str(value).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)


def _axis_fingerprint(names: Any) -> str:
    digest = hashlib.sha256()
    _update_hash_names(digest, names)
    return f"sha256:{digest.hexdigest()}"


def _matrix_fingerprint(matrix: Any, obs_names: Any, var_names: Any) -> str:
    digest = hashlib.sha256()
    digest.update(str(tuple(int(value) for value in matrix.shape)).encode("ascii"))
    digest.update(str(matrix.dtype).encode("ascii"))
    _update_hash_names(digest, obs_names)
    _update_hash_names(digest, var_names)
    if sparse.issparse(matrix):
        canonical = matrix.tocsr(copy=True)
        canonical.sum_duplicates()
        canonical.sort_indices()
        for array in (canonical.indptr, canonical.indices, canonical.data):
            digest.update(np.ascontiguousarray(array).tobytes())
    else:
        digest.update(np.ascontiguousarray(np.asarray(matrix)).tobytes())
    return f"sha256:{digest.hexdigest()}"


def _json_fingerprint(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return f"sha256:{hashlib.sha256(encoded).hexdigest()}"


def _expected_paths(root: Path, name: str) -> dict[str, str]:
    run = root / name
    tmp = run / "cnmf_tmp"
    prefix = tmp / name
    public = run / name
    return {
        "normalized_counts": str(prefix) + ".norm_counts.h5ad",
        "nmf_replicate_parameters": str(prefix) + ".nmf_params.df.npz",
        "nmf_run_parameters": str(prefix) + ".nmf_idvrun_params.yaml",
        "nmf_genes_list": str(public) + ".overdispersed_genes.txt",
        "tpm": str(prefix) + ".tpm.h5ad",
        "tpm_stats": str(prefix) + ".tpm_stats.df.npz",
        "iter_spectra": str(prefix) + ".spectra.k_%d.iter_%d.df.npz",
        "iter_usages": str(prefix) + ".usages.k_%d.iter_%d.df.npz",
        "merged_spectra": str(prefix) + ".spectra.k_%d.merged.df.npz",
        "local_density_cache": str(prefix) + ".local_density_cache.k_%d.merged.df.npz",
        "consensus_spectra": str(prefix) + ".spectra.k_%d.dt_%s.consensus.df.npz",
        "consensus_spectra__txt": str(public) + ".spectra.k_%d.dt_%s.consensus.txt",
        "consensus_usages": str(prefix) + ".usages.k_%d.dt_%s.consensus.df.npz",
        "consensus_usages__txt": str(public) + ".usages.k_%d.dt_%s.consensus.txt",
        "consensus_stats": str(prefix) + ".stats.k_%d.dt_%s.df.npz",
        "clustering_plot": str(public) + ".clustering.k_%d.dt_%s.png",
        "gene_spectra_score": str(prefix) + ".gene_spectra_score.k_%d.dt_%s.df.npz",
        "gene_spectra_score__txt": str(public) + ".gene_spectra_score.k_%d.dt_%s.txt",
        "gene_spectra_tpm": str(prefix) + ".gene_spectra_tpm.k_%d.dt_%s.df.npz",
        "gene_spectra_tpm__txt": str(public) + ".gene_spectra_tpm.k_%d.dt_%s.txt",
        "starcat_spectra": str(prefix) + ".starcat_spectra.k_%d.dt_%s.df.npz",
        "starcat_spectra__txt": str(public) + ".starcat_spectra.k_%d.dt_%s.txt",
        "k_selection_plot": str(public) + ".k_selection.png",
        "k_selection_stats": str(public) + ".k_selection_stats.df.npz",
    }


def _normalized_lexical(path: str | os.PathLike[str]) -> str:
    return os.path.normcase(os.path.abspath(os.fspath(path)))


def _validate_backend_layout(model: Any, root: Path, name: str) -> str:
    output_dir = getattr(model, "output_dir", "")
    if (
        not isinstance(output_dir, (str, os.PathLike))
        or ".." in Path(os.fspath(output_dir)).parts
        or _normalized_lexical(output_dir) != _normalized_lexical(root)
    ):
        raise RuntimeError("cNMF backend output_dir no longer equals the owned private root.")
    if getattr(model, "name", None) != name:
        raise RuntimeError("cNMF backend run name no longer equals the owned private name.")
    paths = getattr(model, "paths", None)
    if not isinstance(paths, Mapping):
        raise RuntimeError("cNMF backend does not expose its public paths mapping.")
    expected = _expected_paths(root, name)
    for key, expected_value in expected.items():
        actual = paths.get(key)
        if (
            not isinstance(actual, (str, os.PathLike))
            or ".." in Path(os.fspath(actual)).parts
            or _normalized_lexical(actual) != _normalized_lexical(expected_value)
        ):
            raise RuntimeError(f"cNMF backend path {key!r} escaped or changed from its canonical private path.")
    payload = {key: _normalized_lexical(paths[key]) for key in sorted(expected)}
    return _json_fingerprint(payload)


def _safe_path(
    root: Path,
    candidate: str | os.PathLike[str],
    *,
    label: str,
    must_exist: bool,
    regular_file: bool = True,
) -> Path:
    raw = os.fspath(candidate)
    if not isinstance(raw, str) or not raw or "\x00" in raw or not os.path.isabs(raw) or ".." in Path(raw).parts:
        raise RuntimeError(f"cNMF {label} path must be an absolute path inside the managed root.")
    root_abs = Path(os.path.abspath(root))
    path = Path(os.path.abspath(raw))
    try:
        common = Path(os.path.commonpath((str(root_abs), str(path))))
    except ValueError as exc:
        raise RuntimeError(f"cNMF {label} path is outside the managed root.") from exc
    if _normalized_lexical(common) != _normalized_lexical(root_abs) or path == root_abs:
        raise RuntimeError(f"cNMF {label} path is outside the managed root.")
    _require_directory_identity(root_abs, label="managed root")
    current = root_abs
    relative = path.relative_to(root_abs)
    for part in relative.parts:
        if part in {"", ".", ".."}:
            raise RuntimeError(f"cNMF {label} path contains an unsafe component.")
        current = current / part
        try:
            value = current.lstat()
        except FileNotFoundError:
            if must_exist:
                raise RuntimeError(f"cNMF {label} file is missing: {current}") from None
            break
        if _is_reparse(value):
            raise RuntimeError(f"cNMF {label} path traverses a link or reparse point: {current}")
    if must_exist:
        try:
            value = path.lstat()
        except FileNotFoundError as exc:
            raise RuntimeError(f"cNMF {label} file is missing: {path}") from exc
        if _is_reparse(value) or (regular_file and not stat.S_ISREG(value.st_mode)):
            raise RuntimeError(f"cNMF {label} must be a regular non-linked file: {path}")
        resolved = path.resolve(strict=True)
        if _normalized_lexical(os.path.commonpath((str(root_abs), str(resolved)))) != _normalized_lexical(root_abs):
            raise RuntimeError(f"cNMF {label} resolved outside the managed root.")
    return path


def _model_path(
    model: Any,
    root: Path,
    key: str,
    *format_values: Any,
    must_exist: bool,
) -> Path:
    paths = getattr(model, "paths", None)
    if not isinstance(paths, Mapping) or key not in paths:
        raise RuntimeError(f"cNMF backend path mapping is missing {key!r}.")
    raw = paths[key]
    try:
        formatted = os.fspath(raw) % format_values if format_values else os.fspath(raw)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"cNMF backend path template {key!r} is malformed.") from exc
    return _safe_path(root, formatted, label=key, must_exist=must_exist)


def _stat_signature(value: os.stat_result) -> tuple[int, int, int, int]:
    return int(value.st_dev), int(value.st_ino), int(value.st_size), int(value.st_mtime_ns)


def _sha256_file(path: Path, root: Path, *, label: str, expected: str | None = None) -> str:
    safe = _safe_path(root, path, label=label, must_exist=True)
    digest = hashlib.sha256()
    with safe.open("rb") as handle:
        before = os.fstat(handle.fileno())
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
        after = os.fstat(handle.fileno())
    current = safe.lstat()
    if _stat_signature(before) != _stat_signature(after) or _stat_signature(after) != _stat_signature(current):
        raise RuntimeError(f"cNMF {label} changed while it was being hashed.")
    result = f"sha256:{digest.hexdigest()}"
    if expected is not None and result != expected:
        raise RuntimeError(f"cNMF {label} failed its immutable SHA-256 check.")
    return result


def _stable_load_frame(
    api: _BackendAPI,
    path: Path,
    root: Path,
    *,
    label: str,
    expected_sha256: str | None = None,
) -> pd.DataFrame:
    before = _sha256_file(path, root, label=label, expected=expected_sha256)
    try:
        value = api.load_frame(str(path))
    except Exception as exc:
        raise RuntimeError(f"Cannot load cNMF {label} table: {exc}") from exc
    after = _sha256_file(path, root, label=label, expected=before)
    if before != after:
        raise RuntimeError(f"cNMF {label} changed while it was being read.")
    if not isinstance(value, pd.DataFrame):
        raise RuntimeError(f"cNMF {label} is not a pandas DataFrame.")
    return value


def _stable_read_h5ad(
    path: Path,
    root: Path,
    *,
    label: str,
    expected_sha256: str | None = None,
) -> Any:
    before = _sha256_file(path, root, label=label, expected=expected_sha256)
    try:
        value = ad.read_h5ad(path)
    except Exception as exc:
        raise RuntimeError(f"Cannot load cNMF {label} H5AD: {exc}") from exc
    after = _sha256_file(path, root, label=label, expected=before)
    if after != before:
        raise RuntimeError(f"cNMF {label} changed while it was being read.")
    return value


def _stable_read_text(path: Path, root: Path, *, label: str) -> str:
    safe = _safe_path(root, path, label=label, must_exist=True)
    with safe.open("r", encoding="utf-8") as handle:
        before = os.fstat(handle.fileno())
        value = handle.read()
        after = os.fstat(handle.fileno())
    current = safe.lstat()
    if _stat_signature(before) != _stat_signature(after) or _stat_signature(after) != _stat_signature(current):
        raise RuntimeError(f"cNMF {label} changed while it was being read.")
    return value


def _stable_read_tsv(
    path: Path,
    root: Path,
    *,
    label: str,
    expected_sha256: str,
) -> pd.DataFrame:
    before = _sha256_file(path, root, label=label, expected=expected_sha256)
    try:
        frame = pd.read_csv(path, sep="\t", index_col=0)
    except Exception as exc:
        raise RuntimeError(f"Cannot load cNMF {label} text table: {exc}") from exc
    _sha256_file(path, root, label=label, expected=before)
    return frame


def _frame_values(frame: Any, *, description: str, nonnegative: bool) -> np.ndarray:
    try:
        values = np.asarray(frame.to_numpy(dtype=float), dtype=float)
    except (AttributeError, TypeError, ValueError) as exc:
        raise RuntimeError(f"cNMF {description} is not a numeric table.") from exc
    if not bool(np.isfinite(values).all()):
        raise RuntimeError(f"cNMF {description} contains non-finite values.")
    if nonnegative and bool((values < 0).any()):
        raise RuntimeError(f"cNMF {description} contains negative values.")
    return values


def _validate_axis(actual: Any, expected: Any, *, description: str) -> None:
    actual_values = [str(value) for value in actual]
    expected_values = [str(value) for value in expected]
    if len(set(actual_values)) != len(actual_values):
        raise RuntimeError(f"cNMF {description} identifiers are duplicated.")
    if actual_values != expected_values:
        raise RuntimeError(f"cNMF {description} axis does not exactly match the expected identifiers and order.")


def _variance_count(matrix: Any) -> int:
    if sparse.issparse(matrix):
        means = np.asarray(matrix.mean(axis=0), dtype=float).reshape(-1)
        squares = np.asarray(matrix.multiply(matrix).mean(axis=0), dtype=float).reshape(-1)
        variances = np.maximum(squares - means**2, 0.0)
    else:
        variances = np.asarray(matrix, dtype=float).var(axis=0)
    return int((variances > np.finfo(float).eps).sum())


def _validate_prepared_files(
    api: _BackendAPI,
    model: Any,
    root: Path,
    *,
    candidate_ks: tuple[int, ...],
    n_iter: int,
    input_obs: Any,
    input_vars: Any,
) -> tuple[tuple[str, ...], int, tuple[tuple[int, int, int], ...]]:
    genes_path = _model_path(model, root, "nmf_genes_list", must_exist=True)
    genes = tuple(
        line.strip() for line in _stable_read_text(genes_path, root, label="HVG list").splitlines() if line.strip()
    )
    if not genes or len(set(genes)) != len(genes):
        raise RuntimeError("cNMF realized high-variance gene list is empty or duplicated.")
    available = {str(value) for value in input_vars}
    unknown = sorted(set(genes) - available)
    if unknown:
        raise RuntimeError(f"cNMF realized high-variance genes are outside the input feature axis: {unknown[:5]}")
    normalized_path = _model_path(model, root, "normalized_counts", must_exist=True)
    normalized = _stable_read_h5ad(normalized_path, root, label="normalized counts")
    _validate_axis(normalized.obs_names, input_obs, description="normalized-count cell")
    _validate_axis(normalized.var_names, genes, description="normalized-count gene")
    values = _matrix_values(normalized.X)
    if not bool(np.isfinite(values).all()) or bool((values < 0).any()):
        raise RuntimeError("cNMF normalized counts contain non-finite or negative values.")
    positive_variance = _variance_count(normalized.X)
    tpm_path = _model_path(model, root, "tpm", must_exist=True)
    tpm = _stable_read_h5ad(tpm_path, root, label="TPM matrix")
    _validate_axis(tpm.obs_names, input_obs, description="TPM cell")
    _validate_axis(tpm.var_names, input_vars, description="TPM gene")
    tpm_values = _matrix_values(tpm.X)
    if not bool(np.isfinite(tpm_values).all()) or bool((tpm_values < 0).any()):
        raise RuntimeError("cNMF TPM matrix contains non-finite or negative values.")

    tpm_stats_path = _model_path(model, root, "tpm_stats", must_exist=True)
    tpm_stats = _stable_load_frame(api, tpm_stats_path, root, label="TPM statistics")
    _validate_axis(tpm_stats.index, input_vars, description="TPM-statistics gene")
    if list(tpm_stats.columns) != ["__mean", "__std"]:
        raise RuntimeError("cNMF TPM statistics must contain exact __mean and __std columns.")
    tpm_stats_values = _frame_values(tpm_stats, description="TPM statistics", nonnegative=False)
    if bool((tpm_stats_values[:, 1] < 0).any()):
        raise RuntimeError("cNMF TPM standard deviations must be non-negative.")

    replicate_path = _model_path(model, root, "nmf_replicate_parameters", must_exist=True)
    replicate = _stable_load_frame(api, replicate_path, root, label="replicate parameters")
    required_columns = ["n_components", "iter", "nmf_seed", "completed"]
    if list(replicate.columns) != required_columns:
        raise RuntimeError(f"cNMF replicate parameters must have exact columns {required_columns}.")
    expected_pairs = [(k, iteration) for k in candidate_ks for iteration in range(n_iter)]
    if len(replicate) != len(expected_pairs):
        raise RuntimeError("cNMF replicate-parameter row count does not match the complete restart family.")
    observed_pairs: list[tuple[int, int]] = []
    restart_seeds: list[tuple[int, int, int]] = []
    for row in replicate.itertuples(index=False):
        k = _integer(row.n_components, name="replicate K", minimum=2)
        iteration = _integer(row.iter, name="replicate iteration", minimum=0)
        seed = _integer(row.nmf_seed, name="replicate seed", minimum=1, maximum=2**31 - 2)
        if bool(row.completed):
            raise RuntimeError("A new private cNMF run unexpectedly reports a pre-completed restart.")
        observed_pairs.append((k, iteration))
        restart_seeds.append((k, iteration, seed))
    if observed_pairs != expected_pairs:
        raise RuntimeError("cNMF replicate parameters do not match the ordered complete K/restart family.")

    run_parameters_path = _model_path(model, root, "nmf_run_parameters", must_exist=True)
    text = _stable_read_text(run_parameters_path, root, label="NMF run parameters")
    try:
        import yaml

        run_parameters = yaml.safe_load(text)
    except Exception as exc:
        raise RuntimeError(f"Cannot parse cNMF NMF run parameters: {exc}") from exc
    expected_parameters = {
        "alpha_H": 0.0,
        "alpha_W": 0.0,
        "beta_loss": "frobenius",
        "init": "random",
        "l1_ratio": 0.0,
        "max_iter": 1000,
        "solver": "cd",
        "tol": 0.0001,
    }
    if run_parameters != expected_parameters:
        raise RuntimeError("cNMF persisted NMF run parameters differ from the fixed audited policy.")
    return genes, positive_variance, tuple(restart_seeds)


def _validate_iteration_files(
    api: _BackendAPI,
    model: Any,
    root: Path,
    candidate_ks: Sequence[int],
    n_iter: int,
    genes: Sequence[str],
) -> dict[int, int]:
    completed: dict[int, int] = {}
    for k in candidate_ks:
        count = 0
        for iteration in range(n_iter):
            path = _model_path(model, root, "iter_spectra", k, iteration, must_exist=True)
            frame = _stable_load_frame(api, path, root, label=f"K={k} restart={iteration} spectrum")
            if frame.shape != (k, len(genes)):
                raise RuntimeError(f"cNMF K={k} restart={iteration} spectrum has an invalid shape.")
            _validate_axis(frame.index, range(1, k + 1), description=f"K={k} restart={iteration} component")
            _validate_axis(frame.columns, genes, description=f"K={k} restart={iteration} gene")
            _frame_values(frame, description=f"K={k} restart={iteration} spectrum", nonnegative=True)
            count += 1
        completed[k] = count
    return completed


def _load_merged(
    api: _BackendAPI,
    model: Any,
    root: Path,
    k: int,
    n_iter: int,
    genes: Sequence[str],
    *,
    expected_sha256: str | None = None,
) -> pd.DataFrame:
    path = _model_path(model, root, "merged_spectra", k, must_exist=True)
    frame = _stable_load_frame(api, path, root, label=f"K={k} merged spectrum", expected_sha256=expected_sha256)
    expected_index = [f"iter{iteration}_topic{topic}" for iteration in range(n_iter) for topic in range(1, k + 1)]
    if frame.shape != (n_iter * k, len(genes)):
        raise RuntimeError(f"cNMF K={k} merged spectrum has an invalid shape.")
    _validate_axis(frame.index, expected_index, description=f"K={k} merged-component")
    _validate_axis(frame.columns, genes, description=f"K={k} merged-spectrum gene")
    _frame_values(frame, description=f"K={k} merged spectrum", nonnegative=True)
    return frame


def _known_warning(warning: warnings.WarningMessage) -> str | None:
    source = Path(str(warning.filename)).as_posix().lower()
    message = str(warning.message)
    if source.endswith("/cnmf/cnmf.py") and int(warning.lineno) in {727, 792}:
        pattern = re.compile(
            r"unclosed file <_io\.TextIOWrapper name='.*\.nmf_idvrun_params\.yaml' mode='r' encoding='[^']+'>"
        )
        if warning.category is ResourceWarning and pattern.fullmatch(message):
            return "cnmf-1.7.1-unclosed-yaml-reader"
    if source.endswith("/cnmf/cnmf.py") and int(warning.lineno) == 964:
        pattern = re.compile(
            r"unclosed file <_io\.TextIOWrapper name='.*\.overdispersed_genes\.txt' mode='r' encoding='[^']+'>"
        )
        if warning.category is ResourceWarning and pattern.fullmatch(message):
            return "cnmf-1.7.1-unclosed-hvg-reader"
    if source.endswith("/cnmf/cnmf.py") and int(warning.lineno) == 969:
        implicit_modification = getattr(ad, "ImplicitModificationWarning", None)
        if (
            implicit_modification is not None
            and issubclass(warning.category, implicit_modification)
            and message == "Setting element `.X` of view, initializing view as actual."
        ):
            return "cnmf-1.7.1-anndata-view-materialization"
    if source.endswith("/cnmf/cnmf.py") and int(warning.lineno) == 916:
        pandas4 = getattr(pd.errors, "Pandas4Warning", None)
        if (
            pandas4 is not None
            and issubclass(warning.category, pandas4)
            and message == "Starting with pandas version 4.0 all arguments of sum will be keyword-only."
        ):
            return "cnmf-1.7.1-pandas4-sum-positional"
    return None


def _invoke_upstream(stage: str, callable_object: Any, /, *args: Any, **kwargs: Any) -> tuple[Any, list[str]]:
    del stage
    caught: list[warnings.WarningMessage]
    result: Any = None
    error: BaseException | None = None
    with warnings.catch_warnings(record=True) as records:
        warnings.simplefilter("always")
        try:
            result = callable_object(*args, **kwargs)
        except BaseException as exc:
            error = exc
        caught = list(records)
    known: list[str] = []
    for warning in caught:
        label = _known_warning(warning)
        if label is not None:
            known.append(label)
            continue
        warnings.warn_explicit(
            warning.message,
            warning.category,
            warning.filename,
            warning.lineno,
            source=warning.source,
        )
    if error is not None:
        raise error
    return result, known


def _warning_counts(labels: Sequence[str]) -> tuple[tuple[str, int], ...]:
    return tuple((label, labels.count(label)) for label in dict.fromkeys(labels))


def _metric_from_statistics(
    statistics: Any,
    *,
    requested_k: int,
    n_iter: int,
    completed_restarts: int,
) -> CNMFKMetric:
    if (
        not isinstance(statistics, pd.DataFrame)
        or list(statistics.index)
        != [
            "k",
            "local_density_threshold",
            "silhouette",
            "prediction_error",
        ]
        or list(statistics.columns) != ["stats"]
    ):
        raise RuntimeError("cNMF statistics must be the exact documented four-row 'stats' DataFrame.")
    values = _frame_values(statistics, description=f"K={requested_k} statistics", nonnegative=False).reshape(-1)
    k_value = _integer(values[0], name="statistics K", minimum=2)
    threshold = float(values[1])
    stability = float(values[2])
    prediction_error = float(values[3])
    if k_value != requested_k:
        raise RuntimeError(f"cNMF statistics returned K={k_value}, expected K={requested_k}.")
    if threshold != CNMF_STATISTICS_DENSITY_THRESHOLD:
        raise RuntimeError("cNMF statistics did not record the fixed no-filter density threshold 2.0.")
    if not -1.0 <= stability <= 1.0:
        raise RuntimeError("cNMF stability must be finite and between -1 and 1.")
    if prediction_error < 0:
        raise RuntimeError("cNMF prediction error must be finite and non-negative.")
    return CNMFKMetric(
        k=requested_k,
        stability=stability,
        prediction_error=prediction_error,
        statistics_density_threshold=threshold,
        expected_restarts=n_iter,
        completed_restarts=completed_restarts,
        combined_components=n_iter * requested_k,
    )


def _relative_hash(root: Path, path: Path, *, label: str) -> tuple[str, str]:
    safe = _safe_path(root, path, label=label, must_exist=True)
    return safe.relative_to(root).as_posix(), _sha256_file(safe, root, label=label)


def _manifest(entries: Sequence[tuple[str, str]]) -> str:
    normalized = sorted((str(path), str(digest)) for path, digest in entries)
    return _json_fingerprint({path: digest for path, digest in normalized})


def _collect_prepare_hashes(model: Any, root: Path, counts_path: Path) -> list[tuple[str, str]]:
    result = [_relative_hash(root, counts_path, label="private count input")]
    for key in (
        "normalized_counts",
        "nmf_replicate_parameters",
        "nmf_run_parameters",
        "nmf_genes_list",
        "tpm",
        "tpm_stats",
    ):
        result.append(_relative_hash(root, _model_path(model, root, key, must_exist=True), label=key))
    return result


def _collect_iteration_hashes(
    model: Any, root: Path, candidate_ks: Sequence[int], n_iter: int
) -> list[tuple[str, str]]:
    result: list[tuple[str, str]] = []
    for k in candidate_ks:
        for iteration in range(n_iter):
            path = _model_path(model, root, "iter_spectra", k, iteration, must_exist=True)
            result.append(_relative_hash(root, path, label=f"K={k} restart={iteration} spectrum"))
    return result


def _collect_merged_hashes(model: Any, root: Path, candidate_ks: Sequence[int]) -> list[tuple[str, str]]:
    return [
        _relative_hash(
            root,
            _model_path(model, root, "merged_spectra", k, must_exist=True),
            label=f"K={k} merged spectrum",
        )
        for k in candidate_ks
    ]


def _verify_artifact_hashes(root: Path, entries: Sequence[tuple[str, str]]) -> str:
    for relative, expected in entries:
        relative_path = Path(relative)
        if relative_path.is_absolute() or ".." in relative_path.parts:
            raise RuntimeError("cNMF immutable artifact manifest contains an unsafe relative path.")
        _sha256_file(root / relative_path, root, label=f"immutable artifact {relative}", expected=expected)
    return _manifest(entries)


def _write_private_counts(root: Path, matrix: Any, obs_names: Any, var_names: Any, fingerprint: str) -> Path:
    counts_path = root / "input_counts.h5ad"
    private = ad.AnnData(
        # The matrix is worker-owned and only serialized here; sharing avoids a full defensive copy.
        X=matrix,
        obs=pd.DataFrame(index=pd.Index([str(value) for value in obs_names])),
        var=pd.DataFrame(index=pd.Index([str(value) for value in var_names])),
    )
    private.write_h5ad(counts_path)
    reopened = _stable_read_h5ad(counts_path, root, label="private count input")
    _validate_axis(reopened.obs_names, obs_names, description="private count-input cell")
    _validate_axis(reopened.var_names, var_names, description="private count-input gene")
    if _matrix_fingerprint(reopened.X, reopened.obs_names, reopened.var_names) != fingerprint:
        raise RuntimeError("cNMF private H5AD count input differs from the declared source after write/reopen.")
    return counts_path


def cnmf_rank_survey(
    adata: Any,
    *,
    source: str = "layer:counts",
    components_min: int = 3,
    components_max: int = 19,
    n_iter: int = 100,
    num_highvar_genes: int = 2000,
    random_seed: int = 123,
) -> tuple[CNMFRun, pd.DataFrame]:
    _load_science()
    source_kind, source_layer = _parse_source(source)
    _require_adata(
        adata,
        operation="cNMF Rank Survey",
        require_current_features=source_kind != "raw",
    )
    components_min = _integer(components_min, name="components_min", minimum=2)
    components_max = _integer(components_max, name="components_max", minimum=2)
    if components_min > components_max:
        raise ValueError("cNMF components_min cannot exceed components_max.")
    n_iter = _integer(n_iter, name="n_iter", minimum=2)
    num_highvar_genes = _integer(
        num_highvar_genes,
        name="num_highvar_genes",
        minimum=1,
    )
    random_seed = _integer(random_seed, name="random_seed", minimum=0, maximum=2**31 - 1)
    (
        matrix,
        source_obs,
        source_vars,
        _,
        _,
        nonzero,
        total_counts,
        input_advisories,
    ) = _validate_counts(adata, source_kind, source_layer)
    source_cells = len(source_obs)
    source_features = len(source_vars)
    current_features = int(adata.n_vars)
    input_fingerprint = _matrix_fingerprint(matrix, source_obs, source_vars)
    candidate_ks = tuple(range(components_min, components_max + 1))
    requested_resource = _validate_resource_budget(
        candidate_ks, n_iter, num_highvar_genes, source_cells, stage="requested-HVG preflight"
    )
    advisories = list(input_advisories)
    if components_max > source_cells:
        advisories.append(
            f"Maximum K={components_max} exceeds the {source_cells} input cells. The official random-init backend "
            "can calculate this overcomplete shape, but programs may be weakly identified and resource use grows quickly."
        )
    if num_highvar_genes > source_features:
        advisories.append(
            f"Requested num_highvar_genes={num_highvar_genes} exceeds the {source_features} source features; "
            "the official backend will retain at most the available eligible genes."
        )
    if requested_resource.estimated_peak_bytes > requested_resource.budget_bytes:
        advisories.append(
            "The conservative requested-HVG working-memory estimate exceeds 2 GiB. This is disclosed rather "
            "than blocked for expert execution; monitor memory and reduce K/restarts/HVGs if necessary."
        )
    fixed_policy = _fixed_policy()
    parameters = {
        "source": source_kind,
        "source_layer": source_layer,
        "source_features": source_features,
        "current_features": current_features,
        "components_min": components_min,
        "components_max": components_max,
        "candidate_ks": list(candidate_ks),
        "n_iter": n_iter,
        "num_highvar_genes": num_highvar_genes,
        "random_seed": random_seed,
        "fixed_policy": asdict(fixed_policy),
        "requested_resource": requested_resource.as_dict(),
    }
    parameter_fingerprint = _json_fingerprint(parameters)
    api = _load_backend()
    temporary = tempfile.TemporaryDirectory(prefix="openbio-cnmf-")
    root = Path(temporary.name).absolute()
    known_warnings: list[str] = []
    try:
        counts_path = _write_private_counts(root, matrix, source_obs, source_vars, input_fingerprint)
        model = api.cnmf_class(output_dir=str(root), name=_RUN_NAME)
        paths_fingerprint = _validate_backend_layout(model, root, _RUN_NAME)
        with _PREPARE_LOCK:
            legacy_rng_state = np.random.get_state()
            try:
                _, captured = _invoke_upstream(
                    "prepare",
                    model.prepare,
                    counts_fn=str(counts_path),
                    components=np.asarray(candidate_ks, dtype=int),
                    n_iter=n_iter,
                    densify=False,
                    tpm_fn=None,
                    seed=random_seed,
                    beta_loss="frobenius",
                    num_highvar_genes=num_highvar_genes,
                    genes_file=None,
                    alpha_usage=0.0,
                    alpha_spectra=0.0,
                    init="random",
                    max_NMF_iter=1000,
                )
            finally:
                np.random.set_state(legacy_rng_state)
        known_warnings.extend(captured)
        if _validate_backend_layout(model, root, _RUN_NAME) != paths_fingerprint:
            raise RuntimeError("cNMF backend paths changed during prepare.")
        genes, positive_variance, restart_seeds = _validate_prepared_files(
            api,
            model,
            root,
            candidate_ks=candidate_ks,
            n_iter=n_iter,
            input_obs=source_obs,
            input_vars=source_vars,
        )
        if components_max > len(genes):
            advisories.append(
                f"Maximum K={components_max} exceeds the {len(genes)} realized HVGs. The official random-init "
                "backend accepted this overcomplete shape; interpret program identifiability cautiously."
            )
        if components_max > positive_variance:
            advisories.append(
                f"Maximum K={components_max} exceeds the {positive_variance} positive-variance realized HVGs. "
                "Execution remains open to experts, but the factorization may be weakly identified."
            )
        persisted_counts = _stable_read_h5ad(counts_path, root, label="private count input")
        _validate_axis(persisted_counts.obs_names, source_obs, description="private count-input cell")
        _validate_axis(persisted_counts.var_names, source_vars, description="private count-input gene")
        if (
            _matrix_fingerprint(persisted_counts.X, persisted_counts.obs_names, persisted_counts.var_names)
            != input_fingerprint
        ):
            raise RuntimeError("cNMF private H5AD count input changed during backend preparation.")
        realized_resource = _validate_resource_budget(
            candidate_ks, n_iter, len(genes), source_cells, stage="realized-HVG post-prepare"
        )
        if realized_resource.estimated_peak_bytes > realized_resource.budget_bytes:
            advisories.append(
                "The conservative realized-HVG working-memory estimate exceeds 2 GiB. Execution remains open "
                "to experts and the complete estimate is recorded."
            )
        prepare_hashes = _collect_prepare_hashes(model, root, counts_path)
        prepare_manifest = _manifest(prepare_hashes)

        _, captured = _invoke_upstream(
            "factorize",
            model.factorize,
            worker_i=0,
            total_workers=1,
            skip_completed_runs=False,
        )
        known_warnings.extend(captured)
        if _validate_backend_layout(model, root, _RUN_NAME) != paths_fingerprint:
            raise RuntimeError("cNMF backend paths changed during factorization.")
        completed = _validate_iteration_files(api, model, root, candidate_ks, n_iter, genes)
        iteration_hashes = _collect_iteration_hashes(model, root, candidate_ks, n_iter)
        factorization_hashes = prepare_hashes + iteration_hashes
        factorization_manifest = _manifest(factorization_hashes)

        _, captured = _invoke_upstream(
            "combine",
            model.combine,
            components=np.asarray(candidate_ks, dtype=int),
            skip_missing_files=False,
        )
        known_warnings.extend(captured)
        if _validate_backend_layout(model, root, _RUN_NAME) != paths_fingerprint:
            raise RuntimeError("cNMF backend paths changed during combine.")
        for k in candidate_ks:
            _load_merged(api, model, root, k, n_iter, genes)
        merged_hashes = _collect_merged_hashes(model, root, candidate_ks)
        artifact_hashes = tuple(sorted(factorization_hashes + merged_hashes))
        artifact_manifest = _manifest(artifact_hashes)

        normalized_path = _model_path(model, root, "normalized_counts", must_exist=True)
        normalized_sha = dict(artifact_hashes)[normalized_path.relative_to(root).as_posix()]
        normalized = _stable_read_h5ad(normalized_path, root, label="normalized counts", expected_sha256=normalized_sha)
        metrics: list[CNMFKMetric] = []
        for k in candidate_ks:
            statistics, captured = _invoke_upstream(
                "statistics consensus",
                model.consensus,
                k=k,
                density_threshold=CNMF_STATISTICS_DENSITY_THRESHOLD,
                local_neighborhood_size=CNMF_LOCAL_NEIGHBORHOOD_SIZE,
                show_clustering=False,
                build_ref=False,
                skip_density_and_return_after_stats=True,
                close_clustergram_fig=True,
                refit_usage=True,
                normalize_tpm_spectra=False,
                norm_counts=normalized,
            )
            known_warnings.extend(captured)
            metrics.append(
                _metric_from_statistics(
                    statistics,
                    requested_k=k,
                    n_iter=n_iter,
                    completed_restarts=completed[k],
                )
            )
        if _verify_artifact_hashes(root, artifact_hashes) != artifact_manifest:
            raise RuntimeError("cNMF immutable artifact manifest changed during K-statistics calculation.")
        if _validate_backend_layout(model, root, _RUN_NAME) != paths_fingerprint:
            raise RuntimeError("cNMF backend paths changed during K-statistics calculation.")
        current_source_obs, current_source_vars = _source_axes(adata, source_kind)
        if (
            _matrix_fingerprint(
                _selected_matrix(adata, source_kind, source_layer),
                current_source_obs,
                current_source_vars,
            )
            != input_fingerprint
        ):
            raise RuntimeError("Input AnnData count source changed during cNMF Rank Survey.")
        metadata = CNMFRunMetadata(
            adapter_version=CNMF_ADAPTER_VERSION,
            backend_name="cnmf.cNMF",
            backend_version=api.version,
            audited_pypi_wheel_sha256=CNMF_AUDITED_PYPI_WHEEL_SHA256,
            installed_distribution_attested=False,
            source_kind=source_kind,
            source_layer=source_layer,
            input_advisories=tuple(advisories),
            input_cells=source_cells,
            input_genes=source_features,
            current_features=current_features,
            input_total_counts=total_counts,
            input_nonzero_entries=nonzero,
            input_fingerprint=input_fingerprint,
            parameters_fingerprint=parameter_fingerprint,
            candidate_ks=candidate_ks,
            n_iter=n_iter,
            num_highvar_genes=num_highvar_genes,
            realized_highvar_genes=genes,
            positive_variance_highvar_genes=positive_variance,
            random_seed=random_seed,
            restart_seeds=restart_seeds,
            total_restarts=requested_resource.total_restarts,
            requested_resource=requested_resource,
            realized_resource=realized_resource,
            fixed_policy=fixed_policy,
            backend_paths_fingerprint=paths_fingerprint,
            prepare_manifest_sha256=prepare_manifest,
            factorization_manifest_sha256=factorization_manifest,
            artifact_manifest_sha256=artifact_manifest,
            artifact_hashes=artifact_hashes,
            known_upstream_warnings=_warning_counts(known_warnings),
        )
        metrics_tuple = tuple(metrics)
        run = CNMFRun(
            metadata=metadata,
            metrics=metrics_tuple,
            backend=model,
            base_adata=_materialize_source_adata(adata, source_kind),
            temporary_directory=temporary,
        )
        table = pd.DataFrame.from_records(
            [metric.as_dict() for metric in metrics_tuple],
            columns=[
                "k",
                "stability",
                "prediction_error",
                "statistics_density_threshold",
                "expected_restarts",
                "completed_restarts",
                "combined_components",
            ],
        )
        return run, table
    except BaseException as primary:
        _cleanup_preserving_primary(
            temporary.cleanup,
            primary,
            context="cNMF Rank Survey managed-temporary rollback",
        )
        raise


def _immutable_hash_for(run: Any, relative_path: str) -> str:
    mapping = dict(run.metadata.artifact_hashes)
    try:
        return mapping[relative_path]
    except KeyError as exc:
        raise RuntimeError(f"cNMF immutable manifest is missing {relative_path!r}.") from exc


def _numeric_summary(values: Any) -> dict[str, int | float | None]:
    array = np.asarray(values, dtype=float).reshape(-1)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return {
            "n": int(array.size),
            "missing": int(array.size),
            "min": None,
            "q1": None,
            "median": None,
            "mean": None,
            "q3": None,
            "max": None,
        }
    quantiles = np.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "n": int(array.size),
        "missing": int(array.size - finite.size),
        "min": float(quantiles[0]),
        "q1": float(quantiles[1]),
        "median": float(quantiles[2]),
        "mean": float(finite.mean()),
        "q3": float(quantiles[3]),
        "max": float(quantiles[4]),
    }


def _independent_density(
    merged: pd.DataFrame,
    *,
    selected_k: int,
    n_iter: int,
    density_threshold: float,
    local_neighborhood_size: float,
) -> tuple[pd.Series, pd.Series]:
    values = _frame_values(merged, description="merged spectrum", nonnegative=True)
    norms = np.sqrt((values**2).sum(axis=1))
    if bool((norms <= 0).any()):
        raise RuntimeError("cNMF merged spectrum contains a zero-length component.")
    l2_values = values / norms[:, None]
    n_neighbors = int(local_neighborhood_size * n_iter)
    if n_neighbors < 1:
        raise ValueError(
            "cNMF local_neighborhood_size is too small for this run: "
            f"int({local_neighborhood_size} * {n_iter}) must be at least one."
        )
    distances = euclidean_distances(l2_values)
    order = np.argpartition(distances, n_neighbors + 1, axis=1)[:, : n_neighbors + 1]
    nearest = distances[np.arange(distances.shape[0])[:, None], order]
    density = pd.Series(nearest.sum(axis=1) / n_neighbors, index=merged.index, name="local_density")
    if not bool(np.isfinite(density.to_numpy()).all()) or bool((density.to_numpy() < 0).any()):
        raise RuntimeError("cNMF independent local density is not finite and non-negative.")
    retained = density < density_threshold
    if int(retained.sum()) < selected_k:
        raise RuntimeError(
            f"cNMF density filtering retained {int(retained.sum())} components, fewer than selected K={selected_k}."
        )
    retained_values = l2_values[retained.to_numpy(), :]
    model = KMeans(n_clusters=selected_k, n_init=10, random_state=1)
    labels = pd.Series(model.fit_predict(retained_values) + 1, index=merged.index[retained], name="cluster")
    label_values = {int(value) for value in labels.to_numpy()}
    expected_labels = set(range(1, selected_k + 1))
    if label_values != expected_labels:
        raise RuntimeError(f"cNMF independent consensus labels must contain every program 1..{selected_k}.")
    return density, labels


def _threshold_token(value: float) -> str:
    return str(value).replace(".", "_")


def _remove_mutable_outputs(model: Any, root: Path, selected_k: int, density_threshold: float) -> None:
    token = _threshold_token(density_threshold)
    specifications: list[tuple[str, tuple[Any, ...]]] = [("local_density_cache", (selected_k,))]
    for key in (
        "consensus_spectra",
        "consensus_spectra__txt",
        "consensus_usages",
        "consensus_usages__txt",
        "gene_spectra_score",
        "gene_spectra_score__txt",
        "gene_spectra_tpm",
        "gene_spectra_tpm__txt",
        "clustering_plot",
    ):
        specifications.append((key, (selected_k, token)))
    for key, values in specifications:
        path = _model_path(model, root, key, *values, must_exist=False)
        if path.exists():
            safe = _safe_path(root, path, label=f"stale {key}", must_exist=True)
            safe.unlink()


def _align_axis(frame: Any, target: Any, *, description: str) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        try:
            frame = pd.DataFrame(frame)
        except Exception as exc:
            raise RuntimeError(f"cNMF {description} is not a table.") from exc
    result = frame.copy()
    actual = [str(value) for value in result.index]
    expected = [str(value) for value in target]
    if len(set(actual)) != len(actual):
        raise RuntimeError(f"cNMF {description} identifiers are duplicated.")
    if set(actual) != set(expected) or len(actual) != len(expected):
        raise RuntimeError(f"cNMF {description} axis does not exactly match the input identifiers.")
    result.index = pd.Index(actual)
    result = result.loc[expected, :]
    result.index = target.copy()
    return result


def _program_number(value: Any) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise RuntimeError(f"cNMF result has invalid program label {value!r}.")
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, str) and value == value.strip() and value.isdigit() and str(int(value)) == value:
        return int(value)
    raise RuntimeError(f"cNMF result has invalid program label {value!r}.")


def _align_programs(frame: Any, k: int, *, description: str) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        try:
            frame = pd.DataFrame(frame)
        except Exception as exc:
            raise RuntimeError(f"cNMF {description} is not a table.") from exc
    numbers = [_program_number(value) for value in frame.columns]
    if len(set(numbers)) != len(numbers) or set(numbers) != set(range(1, k + 1)):
        raise RuntimeError(f"cNMF {description} program labels must be exactly 1..{k}.")
    result = frame.copy()
    result.columns = numbers
    result = result.loc[:, list(range(1, k + 1))]
    result.columns = [f"cNMF_{number}" for number in range(1, k + 1)]
    return result


def _result_paths(model: Any, root: Path, selected_k: int, density_threshold: float) -> dict[str, Path]:
    token = _threshold_token(density_threshold)
    return {
        key: Path(_model_path(model, root, key, selected_k, token, must_exist=True))
        for key in (
            "consensus_spectra",
            "consensus_spectra__txt",
            "consensus_usages",
            "consensus_usages__txt",
            "gene_spectra_score",
            "gene_spectra_score__txt",
            "gene_spectra_tpm",
            "gene_spectra_tpm__txt",
        )
    }


def _verify_base_input(run: Any, output: Any) -> None:
    if run.metadata.source_kind == "raw":
        matrix = output.X
    else:
        matrix = _selected_matrix(output, run.metadata.source_kind, run.metadata.source_layer)
    fingerprint = _matrix_fingerprint(matrix, output.obs_names, output.var_names)
    if fingerprint != run.metadata.input_fingerprint:
        raise RuntimeError("OPENBIO_CNMF_RUN private base AnnData no longer matches its input fingerprint.")


def cnmf_consensus_programs(
    run: Any,
    *,
    selected_k: int = 7,
    density_threshold: float = 2.0,
    local_neighborhood_size: float = CNMF_LOCAL_NEIGHBORHOOD_SIZE,
    n_top_genes: int = 100,
    overwrite_existing: bool = False,
) -> Any:
    _load_science()
    _validate_run_artifact(run)
    if run.metadata.adapter_version != CNMF_ADAPTER_VERSION:
        raise RuntimeError("OPENBIO_CNMF_RUN adapter version is incompatible with this consensus implementation.")
    selected_k = _integer(selected_k, name="selected_k", minimum=2)
    if selected_k not in tuple(run.metadata.candidate_ks):
        raise ValueError(
            f"cNMF selected_k={selected_k} was not surveyed; available values are {list(run.metadata.candidate_ks)}."
        )
    density_threshold = _finite_float(density_threshold, name="density_threshold", minimum=0.0, maximum=2.0)
    local_neighborhood_size = _finite_float(
        local_neighborhood_size, name="local_neighborhood_size", minimum=0.0, maximum=1.0
    )
    if local_neighborhood_size <= 0:
        raise ValueError("cNMF local_neighborhood_size must be greater than zero and at most 1.0.")
    density_neighbors = int(local_neighborhood_size * run.metadata.n_iter)
    if density_neighbors < 1:
        raise ValueError(
            "cNMF local_neighborhood_size is too small for this run: "
            f"int({local_neighborhood_size} * {run.metadata.n_iter}) must be at least one."
        )
    n_top_genes = _integer(n_top_genes, name="n_top_genes", minimum=1)
    if not isinstance(overwrite_existing, (bool, np.bool_)):
        raise TypeError("cNMF overwrite_existing must be a boolean.")
    overwrite_existing = bool(overwrite_existing)
    output = run.copy_base_adata()
    _verify_base_input(run, output)
    collisions = [
        key
        for container, key in (
            (output.obsm, "X_cnmf_usage"),
            (output.varm, "cnmf_gep_scores"),
            (output.varm, "cnmf_gep_tpm"),
            (output.uns, "cnmf"),
        )
        if key in container
    ]
    if collisions and not overwrite_existing:
        raise ValueError(f"cNMF result keys already exist: {collisions}; enable overwrite_existing to replace them.")
    api = _load_backend()
    if run.metadata.backend_name != "cnmf.cNMF" or run.metadata.backend_version != api.version:
        raise RuntimeError("OPENBIO_CNMF_RUN backend identity differs from exact installed cnmf==1.7.1.")
    known = [label for label, count in run.metadata.known_upstream_warnings for _ in range(count)]
    root = Path(run.private_root)
    with run.locked_backend() as model:
        if not isinstance(model, api.cnmf_class):
            raise RuntimeError("OPENBIO_CNMF_RUN backend object is not an instance of exact cnmf==1.7.1 cNMF.")
        paths_fingerprint = _validate_backend_layout(model, root, _RUN_NAME)
        if paths_fingerprint != run.metadata.backend_paths_fingerprint:
            raise RuntimeError("OPENBIO_CNMF_RUN backend paths failed their fingerprint check.")
        manifest_before = _verify_artifact_hashes(root, run.metadata.artifact_hashes)
        if manifest_before != run.metadata.artifact_manifest_sha256:
            raise RuntimeError("OPENBIO_CNMF_RUN immutable artifact manifest is inconsistent.")
        merged_path = _model_path(model, root, "merged_spectra", selected_k, must_exist=True)
        merged_relative = merged_path.relative_to(root).as_posix()
        merged = _load_merged(
            api,
            model,
            root,
            selected_k,
            run.metadata.n_iter,
            run.metadata.realized_highvar_genes,
            expected_sha256=_immutable_hash_for(run, merged_relative),
        )
        normalized_path = _model_path(model, root, "normalized_counts", must_exist=True)
        normalized_relative = normalized_path.relative_to(root).as_posix()
        normalized = _stable_read_h5ad(
            normalized_path,
            root,
            label="normalized counts",
            expected_sha256=_immutable_hash_for(run, normalized_relative),
        )
        _validate_axis(normalized.obs_names, output.obs_names, description="normalized-count cell")
        _validate_axis(
            normalized.var_names,
            run.metadata.realized_highvar_genes,
            description="normalized-count gene",
        )
        density, labels = _independent_density(
            merged,
            selected_k=selected_k,
            n_iter=run.metadata.n_iter,
            density_threshold=density_threshold,
            local_neighborhood_size=local_neighborhood_size,
        )
        _remove_mutable_outputs(model, root, selected_k, density_threshold)
        _, captured = _invoke_upstream(
            "final consensus",
            model.consensus,
            k=selected_k,
            density_threshold=density_threshold,
            local_neighborhood_size=local_neighborhood_size,
            show_clustering=False,
            build_ref=False,
            skip_density_and_return_after_stats=False,
            close_clustergram_fig=True,
            refit_usage=True,
            normalize_tpm_spectra=False,
            norm_counts=normalized,
        )
        known.extend(captured)
        if _validate_backend_layout(model, root, _RUN_NAME) != paths_fingerprint:
            raise RuntimeError("cNMF backend paths changed during final consensus.")
        density_cache_path = _model_path(model, root, "local_density_cache", selected_k, must_exist=True)
        density_cache = _stable_load_frame(api, density_cache_path, root, label="local-density cache")
        if list(density_cache.columns) != ["local_density"]:
            raise RuntimeError("cNMF local-density cache must contain one local_density column.")
        _validate_axis(density_cache.index, merged.index, description="local-density")
        cached_values = _frame_values(density_cache, description="local-density cache", nonnegative=True).reshape(-1)
        if not bool(np.allclose(cached_values, density.to_numpy(), rtol=1e-12, atol=1e-12)):
            raise RuntimeError("cNMF public local-density cache disagrees with the independent 1.7.1 verifier.")
        result_paths = {
            key: Path(path) for key, path in _result_paths(model, root, selected_k, density_threshold).items()
        }
        result_hashes_before = {
            path.relative_to(root).as_posix(): _sha256_file(path, root, label="consensus result")
            for path in result_paths.values()
        }

        def read_result_text(key: str) -> pd.DataFrame:
            path = result_paths[key]
            relative = path.relative_to(root).as_posix()
            return _stable_read_tsv(
                path,
                root,
                label=key,
                expected_sha256=result_hashes_before[relative],
            )

        file_usage = read_result_text("consensus_usages__txt")
        file_usage = file_usage.div(file_usage.sum(axis=1), axis=0)
        file_scores = read_result_text("gene_spectra_score__txt").T
        file_spectra = read_result_text("gene_spectra_tpm__txt").T
        file_top_genes = pd.DataFrame(
            {
                program: file_scores.sort_values(by=program, ascending=False).index[:n_top_genes].tolist()
                for program in file_scores.columns
            }
        )
        result_tuple, captured = _invoke_upstream(
            "load_results",
            model.load_results,
            K=selected_k,
            density_threshold=density_threshold,
            n_top_genes=n_top_genes,
            norm_usage=True,
        )
        known.extend(captured)
        for path in result_paths.values():
            relative = path.relative_to(root).as_posix()
            _sha256_file(path, root, label="consensus result", expected=result_hashes_before[relative])
        manifest_after = _verify_artifact_hashes(root, run.metadata.artifact_hashes)
        if manifest_after != manifest_before:
            raise RuntimeError("cNMF immutable artifacts changed during final consensus.")
    if not isinstance(result_tuple, tuple) or len(result_tuple) != 4:
        raise RuntimeError("cnmf==1.7.1 load_results must return (usage, spectra_scores, spectra_tpm, top_genes).")
    usage_raw, scores_raw, spectra_raw, top_genes_raw = result_tuple
    usage = _align_programs(
        _align_axis(usage_raw, output.obs_names, description="cell usage matrix"),
        selected_k,
        description="cell usage matrix",
    )
    scores = _align_programs(
        _align_axis(scores_raw, output.var_names, description="GEP score matrix"),
        selected_k,
        description="GEP score matrix",
    )
    spectra = _align_programs(
        _align_axis(spectra_raw, output.var_names, description="GEP TPM matrix"),
        selected_k,
        description="GEP TPM matrix",
    )
    usage_values = _frame_values(usage, description="normalized usage matrix", nonnegative=True)
    _frame_values(scores, description="GEP score matrix", nonnegative=False)
    _frame_values(spectra, description="GEP TPM matrix", nonnegative=True)
    row_sums = usage_values.sum(axis=1)
    if bool((row_sums <= 0).any()) or not bool(np.allclose(row_sums, 1.0, rtol=1e-6, atol=1e-8)):
        raise RuntimeError("cNMF normalized usage rows must be positive and sum to one.")
    top_genes = _align_programs(top_genes_raw, selected_k, description="top-gene table")
    expected_top = min(n_top_genes, int(output.n_vars))
    if len(top_genes) != expected_top:
        raise RuntimeError(f"cNMF returned {len(top_genes)} top genes per program; expected {expected_top}.")
    file_usage_aligned = _align_programs(
        _align_axis(file_usage, output.obs_names, description="persisted cell usage matrix"),
        selected_k,
        description="persisted cell usage matrix",
    )
    file_scores_aligned = _align_programs(
        _align_axis(file_scores, output.var_names, description="persisted GEP score matrix"),
        selected_k,
        description="persisted GEP score matrix",
    )
    file_spectra_aligned = _align_programs(
        _align_axis(file_spectra, output.var_names, description="persisted GEP TPM matrix"),
        selected_k,
        description="persisted GEP TPM matrix",
    )
    file_top_aligned = _align_programs(file_top_genes, selected_k, description="persisted top-gene table")
    for returned, persisted, description in (
        (usage, file_usage_aligned, "normalized usage"),
        (scores, file_scores_aligned, "GEP scores"),
        (spectra, file_spectra_aligned, "GEP TPM spectra"),
    ):
        if not bool(
            np.allclose(
                returned.to_numpy(dtype=float),
                persisted.to_numpy(dtype=float),
                rtol=1e-12,
                atol=1e-12,
            )
        ):
            raise RuntimeError(f"cNMF load_results {description} disagree with the persisted official text file.")
    if top_genes.astype(str).to_numpy().tolist() != file_top_aligned.astype(str).to_numpy().tolist():
        raise RuntimeError("cNMF load_results top genes disagree with the persisted official score ranking.")
    available_genes = {str(value) for value in output.var_names}
    top_records: dict[str, list[str]] = {}
    for program in top_genes.columns:
        genes = [str(value) for value in top_genes[program].tolist()]
        if any(gene not in available_genes for gene in genes):
            raise RuntimeError(f"cNMF top genes for {program} are outside the input feature axis.")
        if len(set(genes)) != len(genes):
            raise RuntimeError(f"cNMF top genes for {program} contain duplicate identifiers.")
        top_records[str(program)] = genes
    program_names = [f"cNMF_{number}" for number in range(1, selected_k + 1)]
    output.obsm["X_cnmf_usage"] = usage
    output.varm["cnmf_gep_scores"] = scores
    output.varm["cnmf_gep_tpm"] = spectra
    selected_metric = next(metric for metric in run.metrics if metric.k == selected_k)
    metric_records = [metric.as_dict() for metric in run.metrics]
    metric_columns = {key: [record[key] for record in metric_records] for key in metric_records[0]}
    components_before = int(len(merged))
    components_after = int(len(labels))
    filtered = components_before - components_after
    cluster_counts = {f"cNMF_{number}": int((labels.to_numpy() == number).sum()) for number in range(1, selected_k + 1)}
    output.uns["cnmf"] = {
        "schema_version": 2,
        "adapter_version": CNMF_ADAPTER_VERSION,
        "selected_k": selected_k,
        "candidate_ks": list(run.metadata.candidate_ks),
        "k_metrics": metric_columns,
        "density_threshold": density_threshold,
        "local_neighborhood_size": local_neighborhood_size,
        "density_neighbors": density_neighbors,
        "components_before_filtering": components_before,
        "components_after_filtering": components_after,
        "components_filtered": filtered,
        "filtering_fraction": filtered / components_before,
        "local_density_summary": _numeric_summary(density.to_numpy()),
        "cluster_component_counts": cluster_counts,
        "program_names": program_names,
        "top_genes": top_records,
        "observation_axis_fingerprint_sha256": _axis_fingerprint(output.obs_names),
        "feature_axis_fingerprint_sha256": _axis_fingerprint(output.var_names),
        "usage_fingerprint_sha256": _matrix_fingerprint(usage_values, output.obs_names, program_names),
        "gep_scores_fingerprint_sha256": _matrix_fingerprint(
            scores.to_numpy(dtype=float), output.var_names, program_names
        ),
        "top_genes_fingerprint_sha256": _json_fingerprint(top_records),
        "storage_keys": {
            "usage": "obsm:X_cnmf_usage",
            "gep_scores": "varm:cnmf_gep_scores",
            "gep_tpm": "varm:cnmf_gep_tpm",
        },
        "overwrote_existing": bool(collisions),
        "source": {
            "kind": run.metadata.source_kind,
            **({"layer": run.metadata.source_layer} if run.metadata.source_layer is not None else {}),
        },
        "survey": {
            "source_features": run.metadata.input_genes,
            "current_features": run.metadata.current_features,
            "input_advisories": list(run.metadata.input_advisories),
            "n_iter": run.metadata.n_iter,
            "num_highvar_genes": run.metadata.num_highvar_genes,
            "realized_highvar_genes": len(run.metadata.realized_highvar_genes),
            "positive_variance_highvar_genes": run.metadata.positive_variance_highvar_genes,
            "random_seed": run.metadata.random_seed,
            "input_fingerprint": run.metadata.input_fingerprint,
            "parameters_fingerprint": run.metadata.parameters_fingerprint,
            "backend": run.metadata.backend_name,
            "backend_version": run.metadata.backend_version,
            "audited_pypi_wheel_sha256": run.metadata.audited_pypi_wheel_sha256,
            "installed_distribution_attested": run.metadata.installed_distribution_attested,
            "total_restarts": run.metadata.total_restarts,
            "completed_restarts": sum(metric.completed_restarts for metric in run.metrics),
            "k_metrics": metric_columns,
            "requested_resource": run.metadata.requested_resource.as_dict(),
            "realized_resource": run.metadata.realized_resource.as_dict(),
            "fixed_policy": asdict(run.metadata.fixed_policy),
            "backend_paths_fingerprint": run.metadata.backend_paths_fingerprint,
            "prepare_manifest_sha256": run.metadata.prepare_manifest_sha256,
            "factorization_manifest_sha256": run.metadata.factorization_manifest_sha256,
            "artifact_manifest_sha256": run.metadata.artifact_manifest_sha256,
            "immutable_manifest_before_consensus": manifest_before,
            "immutable_manifest_after_consensus": manifest_after,
        },
        "selected_k_survey_stability": selected_metric.stability,
        "selected_k_survey_prediction_error": selected_metric.prediction_error,
        "known_upstream_warnings": [{"id": label, "count": count} for label, count in _warning_counts(known)],
    }
    return output


def equivalent_source() -> str:
    with Path(__file__).open("r", encoding="utf-8") as handle:
        return handle.read()


__all__ = [
    "CNMF_ADAPTER_VERSION",
    "CNMF_AUDITED_PYPI_WHEEL_SHA256",
    "CNMFFixedPolicy",
    "CNMFKMetric",
    "CNMF_LOCAL_NEIGHBORHOOD_SIZE",
    "CNMF_REQUIRED_VERSION",
    "CNMFResourceEstimate",
    "CNMFRun",
    "CNMFRunMetadata",
    "close_run_preserving_primary",
    "cnmf_consensus_programs",
    "cnmf_rank_survey",
    "equivalent_source",
]
