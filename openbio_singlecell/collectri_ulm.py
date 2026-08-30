from __future__ import annotations

import copy
import csv
import hashlib
import importlib
import inspect
import json
import math
import os
from collections.abc import Mapping, Sequence
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .analysis_reporting import _package_version, collect_software_versions
from .tf_activity_artifact import (
    TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION,
    TF_ACTIVITY_ARTIFACT_TYPE,
    TF_ACTIVITY_PRODUCER_NODE_ID,
    TF_ACTIVITY_PRODUCER_SCHEMA,
    TFActivityArtifact,
    _activity_frame_fingerprint,
    _canonical_axis,
    _canonical_json_sha256,
    _portable_tf_activity_payload,
    build_tf_activity_artifact,
    validate_tf_activity_artifact,
)

if TYPE_CHECKING:
    from anndata import AnnData


DECOUPLER_VERSION = "2.2.0"
COLLECTRI_SCORE_KEY = "collectri_ulm_score"
COLLECTRI_PADJ_KEY = "collectri_ulm_padj"
COLLECTRI_UNS_KEY = "openbio_collectri_ulm_activity"
COLLECTRI_SUMMARY_SCHEMA = "openbio-singlecell/collectri-ulm-summary/v1"
COLLECTRI_RESOURCE_METADATA_FIELDS = (
    "name",
    "version",
    "date",
    "organism",
    "identifier_namespace",
    "scope",
    "license",
    "citation",
)
COLLECTRI_REFERENCES = [
    {
        "citation": (
            "Müller-Dott S, et al. Expanding the coverage of regulons from high-confidence prior knowledge "
            "for accurate estimation of transcription factor activities. Nucleic Acids Research. "
            "2023;51:10934-10949."
        ),
        "doi": "10.1093/nar/gkad841",
        "url": "https://doi.org/10.1093/nar/gkad841",
        "kind": "method",
    },
    {
        "citation": (
            "Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological "
            "activities from omics data. Bioinformatics Advances. 2022;2:vbac016."
        ),
        "doi": "10.1093/bioadv/vbac016",
        "url": "https://doi.org/10.1093/bioadv/vbac016",
        "kind": "method",
    },
    {
        "citation": "CollecTRI version 2.0 publication resource, Zenodo record 8192729.",
        "doi": "10.5281/zenodo.8192729",
        "url": "https://doi.org/10.5281/zenodo.8192729",
        "kind": "software_documentation",
    },
    {
        "citation": "decoupler 2.2.0 CollecTRI and ULM public interfaces and tagged implementation.",
        "doi": None,
        "url": "https://github.com/scverse/decoupler/tree/v2.2.0/src/decoupler",
        "kind": "software_documentation",
    },
]


def _strict_axis(values: Sequence[Any], *, label: str) -> tuple[list[str], str]:
    canonical: list[str] = []
    for position, value in enumerate(values):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(
                f"CollecTRI ULM {label} at position {position} must be a nonblank, whitespace-canonical string."
            )
        canonical.append(value)
    if len(canonical) != len(set(canonical)):
        raise ValueError(f"CollecTRI ULM {label} values must be unique.")
    payload = json.dumps(canonical, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    return canonical, hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _strict_metadata(resource_metadata_json: str) -> dict[str, str]:
    if not isinstance(resource_metadata_json, str):
        raise TypeError("CollecTRI local resource metadata must be JSON text.")
    if len(resource_metadata_json.encode("utf-8")) > 65_536:
        raise ValueError("CollecTRI local resource metadata exceeds 65,536 UTF-8 bytes.")

    def reject_constant(value: str) -> None:
        raise ValueError(f"CollecTRI local resource metadata contains non-standard constant {value!r}.")

    def pairs_to_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if not isinstance(key, str) or not key or key != key.strip():
                raise ValueError("CollecTRI local resource metadata keys must be canonical nonblank strings.")
            if key in result:
                raise ValueError(f"CollecTRI local resource metadata contains duplicate key {key!r}.")
            result[key] = value
        return result

    try:
        value = json.loads(
            resource_metadata_json,
            object_pairs_hook=pairs_to_object,
            parse_constant=reject_constant,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"CollecTRI local resource metadata is invalid JSON ({exc.msg}).") from exc
    required = set(COLLECTRI_RESOURCE_METADATA_FIELDS)
    observed = set(value) if isinstance(value, dict) else set()
    if not isinstance(value, dict) or observed != required:
        raise ValueError(
            "CollecTRI local resource metadata must contain exactly the required fields; "
            f"missing={sorted(required - observed)}, unknown={sorted(observed - required)}."
        )
    result: dict[str, str] = {}
    for field in COLLECTRI_RESOURCE_METADATA_FIELDS:
        item = value[field]
        if not isinstance(item, str):
            raise TypeError(f"CollecTRI resource metadata field {field!r} must be a string.")
        if not item or item != item.strip():
            raise ValueError(f"CollecTRI resource metadata field {field!r} must be canonical and nonblank.")
        result[field] = item
    return result


def _file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collectri_resource_cache_fingerprint(path: str, identity: tuple[Any, ...]) -> tuple[Any, ...]:
    """Bind Comfy cache invalidation to local CollecTRI bytes as well as filesystem identity."""

    return ("openbio-collectri-resource-v1", *identity, _file_sha256(path))


def _canonical_network_sha256(network: Any) -> str:
    records = [
        [str(source), str(target), float(weight)]
        for source, target, weight in network[["source", "target", "weight"]].itertuples(index=False, name=None)
    ]
    payload = json.dumps(records, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _normalize_network(
    frame: Any,
    *,
    pandas: Any,
    complex_policy: str,
    official: bool,
) -> tuple[Any, dict[str, Any]]:
    if not isinstance(frame, pandas.DataFrame):
        raise TypeError("CollecTRI resource loader must return a pandas DataFrame.")
    required = {"source", "target", "weight"}
    official_columns = required | {"resources", "references", "sign_decision"}
    columns = set(frame.columns)
    if official and columns != official_columns:
        raise ValueError(
            "decoupler 2.2.0 CollecTRI output schema changed; "
            f"expected={sorted(official_columns)}, observed={sorted(columns)}."
        )
    if not official and columns not in {frozenset(required), frozenset(official_columns)}:
        raise ValueError(
            "Local CollecTRI CSV columns must be exactly source,target,weight or the six-column 2.2 table; "
            f"observed={list(frame.columns)!r}."
        )
    work = frame.copy(deep=True)
    input_rows = int(len(work))
    identifiers: dict[str, list[str]] = {}
    for column in ("source", "target"):
        canonical: list[str] = []
        for position, value in enumerate(work[column].tolist()):
            if not isinstance(value, str) or not value or value != value.strip():
                raise ValueError(
                    f"CollecTRI resource {column} row {position + 1} must be a canonical nonblank string."
                )
            canonical.append(value)
        identifiers[column] = canonical
    try:
        weights = work["weight"].to_numpy(dtype=float, copy=True)
    except (TypeError, ValueError) as exc:
        raise TypeError("CollecTRI resource weights must be numeric.") from exc
    if not bool(pandas.notna(weights).all()) or not bool(((weights == -1.0) | (weights == 1.0)).all()):
        raise ValueError("CollecTRI resource weights must be finite and exactly -1 or +1.")
    work["source"] = identifiers["source"]
    work["target"] = identifiers["target"]
    work["weight"] = weights
    if set(work.columns) == official_columns:
        for column in ("resources", "references", "sign_decision"):
            if bool(work[column].isna().any()):
                raise ValueError(f"CollecTRI resource column {column!r} contains missing values.")
    conflicts = work.groupby(["source", "target"], sort=False, observed=True)["weight"].nunique()
    if bool((conflicts > 1).any()):
        pair = conflicts[conflicts > 1].index[0]
        raise ValueError(f"CollecTRI resource contains conflicting weights for source-target pair {pair!r}.")
    duplicate_rows_removed = int(work.duplicated(["source", "target"]).sum())
    work = work.drop_duplicates(["source", "target"], keep="first")
    complexes_removed = 0
    if complex_policy == "remove":
        mask = work["source"].isin(["AP1", "NFKB"])
        complexes_removed = int(mask.sum())
        work = work.loc[~mask]
    elif complex_policy != "retain":
        raise ValueError("CollecTRI complex_policy must be 'retain' or 'remove'.")
    work = work[["source", "target", "weight"]].sort_values(
        ["source", "target"], kind="mergesort", ignore_index=True
    )
    if work.empty:
        raise ValueError("CollecTRI resource is empty after applying the complex policy.")
    if bool(work.duplicated(["source", "target"]).any()):
        raise RuntimeError("CollecTRI canonical network contains duplicate source-target pairs.")
    return work, {
        "input_rows": input_rows,
        "unique_source_target_rows": int(len(work) + complexes_removed),
        "duplicate_rows_removed": duplicate_rows_removed,
        "complex_edges_removed": complexes_removed,
    }


def load_local_collectri_network(
    path: str,
    resource_metadata_json: str,
    *,
    complex_policy: str,
    pandas: Any,
    expected_file_sha256: str | None = None,
) -> tuple[Any, dict[str, Any]]:
    metadata = _strict_metadata(resource_metadata_json)
    if not isinstance(path, str) or not path or path != path.strip():
        raise ValueError("CollecTRI local network path must be a canonical nonblank string.")
    resolved = os.path.realpath(path)
    if not os.path.isfile(resolved):
        raise FileNotFoundError(f"CollecTRI local network file not found: {resolved}")
    if os.path.splitext(resolved)[1].lower() != ".csv":
        raise ValueError("CollecTRI local network must use the .csv extension.")
    size_bytes = int(os.path.getsize(resolved))
    if size_bytes > 512 * 1024 * 1024:
        raise ValueError("CollecTRI local network exceeds the 512 MiB safety limit.")
    before_hash = _file_sha256(resolved)
    if expected_file_sha256 is not None:
        if (
            not isinstance(expected_file_sha256, str)
            or len(expected_file_sha256) != 64
            or any(character not in "0123456789abcdef" for character in expected_file_sha256)
        ):
            raise ValueError("Expected CollecTRI resource SHA-256 is invalid.")
        if before_hash != expected_file_sha256:
            raise ValueError(
                f"CollecTRI local resource changed: expected {expected_file_sha256}, observed {before_hash}."
            )
    try:
        with open(resolved, encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValueError("CollecTRI local network is empty.") from exc
            if not header or len(header) != len(set(header)):
                raise ValueError("CollecTRI local network header must be nonempty and unique.")
            for column in header:
                if not column or column != column.strip():
                    raise ValueError("CollecTRI local network column names must be canonical nonblank strings.")
            allowed = ["source", "target", "weight"]
            official = ["source", "target", "weight", "resources", "references", "sign_decision"]
            if header != allowed and header != official:
                raise ValueError(
                    "CollecTRI local CSV header must be exactly source,target,weight or "
                    "source,target,weight,resources,references,sign_decision in that order."
                )
            rows: list[list[str]] = []
            for row_number, row in enumerate(reader, start=1):
                if row_number > 2_000_000:
                    raise ValueError("CollecTRI local network exceeds the 2,000,000-row safety limit.")
                if len(row) != len(header):
                    raise ValueError(
                        f"CollecTRI local CSV row {row_number} has {len(row)} fields; expected {len(header)}."
                    )
                rows.append(row)
    except UnicodeError as exc:
        raise ValueError("CollecTRI local network must be valid UTF-8; an optional BOM is accepted.") from exc
    if not rows:
        raise ValueError("CollecTRI local network contains no data rows.")
    frame = pandas.DataFrame(rows, columns=header)
    network, accounting = _normalize_network(
        frame, pandas=pandas, complex_policy=complex_policy, official=False
    )
    after_hash = _file_sha256(resolved)
    if before_hash != after_hash or size_bytes != int(os.path.getsize(resolved)):
        raise ValueError("CollecTRI local network changed while it was being parsed.")
    provenance = {
        "mode": "local_network",
        "requested_path": path,
        "resolved_path": resolved,
        "size_bytes": size_bytes,
        "file_sha256": after_hash,
        "canonical_network_sha256": _canonical_network_sha256(network),
        "metadata": metadata,
        "network_access": False,
        "cache_state": "not_applicable_local_snapshot",
        "accounting": accounting,
    }
    return network, provenance


def _signature_shape(callable_value: Any) -> list[tuple[str, str, Any]]:
    result: list[tuple[str, str, Any]] = []
    for parameter in inspect.signature(callable_value).parameters.values():
        default = "<required>" if parameter.default is inspect.Parameter.empty else parameter.default
        result.append((parameter.name, parameter.kind.name, default))
    return result


def _validate_decoupler_ulm(decoupler: Any, *, require_collectri: bool) -> None:
    if getattr(decoupler, "__version__", None) != DECOUPLER_VERSION:
        raise RuntimeError(
            f"CollecTRI ULM requires decoupler exactly {DECOUPLER_VERSION}; "
            f"observed {getattr(decoupler, '__version__', None)!r}."
        )
    try:
        method = decoupler.mt.ulm
    except AttributeError as exc:
        raise RuntimeError("decoupler 2.2.0 public mt.ulm interface is unavailable.") from exc
    expected_public = [
        ("data", "POSITIONAL_OR_KEYWORD", "<required>"),
        ("net", "POSITIONAL_OR_KEYWORD", "<required>"),
        ("tmin", "POSITIONAL_OR_KEYWORD", 5),
        ("raw", "POSITIONAL_OR_KEYWORD", False),
        ("empty", "POSITIONAL_OR_KEYWORD", True),
        ("bsize", "POSITIONAL_OR_KEYWORD", 250000),
        ("verbose", "POSITIONAL_OR_KEYWORD", False),
        ("kwargs", "VAR_KEYWORD", "<required>"),
    ]
    expected_func = [
        ("mat", "POSITIONAL_OR_KEYWORD", "<required>"),
        ("adj", "POSITIONAL_OR_KEYWORD", "<required>"),
        ("tval", "POSITIONAL_OR_KEYWORD", True),
        ("verbose", "POSITIONAL_OR_KEYWORD", False),
    ]
    if _signature_shape(method) != expected_public or not hasattr(method, "func"):
        raise RuntimeError("decoupler 2.2.0 public mt.ulm signature changed.")
    if _signature_shape(method.func) != expected_func:
        raise RuntimeError("decoupler 2.2.0 mt.ulm.func signature changed.")
    if require_collectri:
        try:
            collectri = decoupler.op.collectri
        except AttributeError as exc:
            raise RuntimeError("decoupler 2.2.0 public op.collectri interface is unavailable.") from exc
        expected_collectri = [
            ("organism", "POSITIONAL_OR_KEYWORD", "human"),
            ("remove_complexes", "POSITIONAL_OR_KEYWORD", False),
            ("license", "POSITIONAL_OR_KEYWORD", "academic"),
            ("verbose", "POSITIONAL_OR_KEYWORD", False),
        ]
        if _signature_shape(collectri) != expected_collectri:
            raise RuntimeError("decoupler 2.2.0 public op.collectri signature changed.")


def _matrix_fingerprint(
    matrix: Any,
    *,
    observations: Sequence[str],
    features: Sequence[str],
    numpy: Any,
    sparse: Any,
) -> str:
    digest = hashlib.sha256()
    digest.update(b"openbio-singlecell/collectri-expression/v1\0")
    digest.update(
        json.dumps(
            {"observations": list(observations), "features": list(features)},
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    if sparse.issparse(matrix):
        canonical = sparse.csr_matrix(matrix, dtype=float, copy=True)
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
        for array, dtype in (
            (canonical.indptr, "<i8"),
            (canonical.indices, "<i8"),
            (canonical.data, "<f8"),
        ):
            digest.update(numpy.ascontiguousarray(array, dtype=dtype).tobytes(order="C"))
    else:
        digest.update(numpy.ascontiguousarray(matrix, dtype="<f8").tobytes(order="C"))
    return digest.hexdigest()


def _bh_adjust(pvalues: Any, *, numpy: Any) -> Any:
    values = numpy.asarray(pvalues, dtype=float)
    order = numpy.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = ranked * values.size / numpy.arange(1, values.size + 1, dtype=float)
    adjusted = numpy.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = numpy.clip(adjusted, 0.0, 1.0)
    result = numpy.empty_like(adjusted)
    result[order] = adjusted
    return result


def run_collectri_ulm(
    adata: AnnData,
    *,
    resource_mode: str = "local_network",
    network_path: str | None = None,
    resource_metadata_json: str = "{}",
    source_kind: str = "X",
    layer_name: str | None = None,
    affiliation_license: str = "academic",
    allow_network_access: bool = False,
    complex_policy: str = "retain",
    min_targets: int = 5,
    batch_size: int = 250_000,
    max_output_rows: int = 2_000_000,
    max_working_memory_gib: float = 4.0,
    overwrite_existing: bool = False,
    expected_resource_sha256: str | None = None,
    expected_canonical_network_sha256: str | None = None,
    decoupler_module: Any | None = None,
    openbio_version: str = PLUGIN_VERSION,
    _materialized_official_network_path: str | None = None,
    _resolved_official_resource_provenance: Mapping[str, Any] | None = None,
) -> tuple[AnnData, TFActivityArtifact, dict[str, Any]]:
    """Infer observation-level CollecTRI ULM activity in one worker-owned AnnData."""

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import sparse

    if resource_mode not in {"local_network", "official_collectri"}:
        raise ValueError("CollecTRI resource_mode must be 'local_network' or 'official_collectri'.")
    if affiliation_license not in {"academic", "commercial", "nonprofit"}:
        raise ValueError("CollecTRI affiliation_license must be academic, commercial, or nonprofit.")
    if not isinstance(allow_network_access, bool):
        raise TypeError("CollecTRI allow_network_access must be boolean.")
    if not isinstance(overwrite_existing, bool):
        raise TypeError("CollecTRI overwrite_existing must be boolean.")
    for value, label in (
        (min_targets, "min_targets"),
        (batch_size, "batch_size"),
        (max_output_rows, "max_output_rows"),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise TypeError(f"CollecTRI {label} must be a positive integer.")
    if not isinstance(max_working_memory_gib, (int, float)) or isinstance(max_working_memory_gib, bool):
        raise TypeError("CollecTRI max_working_memory_gib must be numeric.")
    max_working_memory_gib = float(max_working_memory_gib)
    if not math.isfinite(max_working_memory_gib) or max_working_memory_gib <= 0:
        raise ValueError("CollecTRI max_working_memory_gib must be finite and positive.")
    if int(getattr(adata, "n_obs", 0)) < 1 or int(getattr(adata, "n_vars", 0)) < 1:
        raise ValueError("CollecTRI ULM requires a nonempty AnnData input.")
    observations, observation_axis_sha256 = _strict_axis(
        adata.obs_names.tolist(), label="observation identifier"
    )
    if source_kind == "X":
        matrix = adata.X
        feature_values = adata.var_names.tolist()
    elif source_kind == "layer":
        if not isinstance(layer_name, str) or not layer_name or layer_name != layer_name.strip():
            raise ValueError("CollecTRI named layer must be a canonical nonblank string.")
        if layer_name not in adata.layers:
            raise ValueError(f"CollecTRI expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
        feature_values = adata.var_names.tolist()
    else:
        raise ValueError("CollecTRI source_kind must be 'X' or 'layer'; Raw snapshot is not accepted.")
    features, feature_axis_sha256 = _strict_axis(feature_values, label="feature identifier")
    if len(features) < 3:
        raise ValueError("CollecTRI ULM requires at least three expression features for positive residual df.")
    if tuple(getattr(matrix, "shape", ())) != (len(observations), len(features)):
        raise ValueError("CollecTRI selected expression matrix does not align to its named axes.")
    if not (isinstance(matrix, np.ndarray) or sparse.issparse(matrix)):
        raise TypeError("CollecTRI ULM requires an in-memory NumPy or SciPy sparse expression matrix.")
    stored = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    stored = np.asarray(stored)
    if stored.dtype.kind not in "iuf" or not bool(np.isfinite(stored).all()):
        raise ValueError("CollecTRI selected expression must be finite numeric data.")
    if sparse.issparse(matrix):
        canonical_sparse = sparse.csr_matrix(matrix, dtype=float, copy=True)
        canonical_sparse.sum_duplicates()
        canonical_sparse.eliminate_zeros()
        canonical_sparse.sort_indices()
        row_mean = np.asarray(canonical_sparse.sum(axis=1), dtype=float).ravel() / len(features)
        row_second_moment = (
            np.asarray(canonical_sparse.multiply(canonical_sparse).sum(axis=1), dtype=float).ravel()
            / len(features)
        )
        row_variance = row_second_moment - row_mean**2
        semantic_values = canonical_sparse.data
        value_min = float(min(0.0, semantic_values.min(initial=0.0)))
        value_max = float(max(0.0, semantic_values.max(initial=0.0)))
        integer_like = bool(
            np.allclose(semantic_values, np.rint(semantic_values), rtol=0.0, atol=1e-8)
        )
    else:
        dense_view = np.asarray(matrix, dtype=float)
        row_variance = np.var(dense_view, axis=1)
        value_min = float(dense_view.min())
        value_max = float(dense_view.max())
        integer_like = bool(np.allclose(dense_view, np.rint(dense_view), rtol=0.0, atol=1e-8))
    if not bool((row_variance > 0.0).all()):
        raise ValueError("CollecTRI ULM requires nonconstant expression across features for every observation.")
    warnings: list[str] = []
    if integer_like and value_min >= 0.0:
        warnings.append(
            "The explicitly selected expression is count-like. ULM remains numerically defined, but its normal-error "
            "interpretation and activity scale may not match the recommended normalized continuous input."
        )
    expression_sha256 = _matrix_fingerprint(
        matrix, observations=observations, features=features, numpy=np, sparse=sparse
    )

    decoupler = importlib.import_module("decoupler") if decoupler_module is None else decoupler_module
    _validate_decoupler_ulm(
        decoupler,
        require_collectri=(
            resource_mode == "official_collectri" and _materialized_official_network_path is None
        ),
    )
    if resource_mode == "local_network":
        if network_path is None:
            raise ValueError("CollecTRI local network mode requires network_path.")
        network, resource_provenance = load_local_collectri_network(
            network_path,
            resource_metadata_json,
            complex_policy=complex_policy,
            pandas=pd,
            expected_file_sha256=expected_resource_sha256,
        )
    else:
        if network_path is not None:
            raise ValueError("Official CollecTRI mode does not accept a local network_path.")
        if _materialized_official_network_path is not None:
            if not isinstance(_resolved_official_resource_provenance, Mapping):
                raise ValueError(
                    "Equivalent official CollecTRI execution requires the runtime-resolved resource provenance."
                )
            resolved_provenance = copy.deepcopy(dict(_resolved_official_resource_provenance))
            expected_provenance_fields = {
                "mode",
                "requested_path",
                "resolved_path",
                "size_bytes",
                "file_sha256",
                "canonical_network_sha256",
                "metadata",
                "affiliation_license_declaration",
                "affiliation_filter_enforced_by_decoupler_2_2",
                "network_access",
                "cache_state",
                "accounting",
            }
            if set(resolved_provenance) != expected_provenance_fields:
                raise ValueError("Equivalent official CollecTRI resource provenance schema is invalid.")
            materialized_metadata = json.dumps(
                resolved_provenance["metadata"],
                ensure_ascii=False,
                allow_nan=False,
                separators=(",", ":"),
            )
            network, _materialized_provenance = load_local_collectri_network(
                _materialized_official_network_path,
                materialized_metadata,
                complex_policy="retain",
                pandas=pd,
            )
            observed_materialized_hash = _canonical_network_sha256(network)
            if observed_materialized_hash != resolved_provenance["canonical_network_sha256"]:
                raise ValueError(
                    "Materialized official CollecTRI network does not match the runtime-resolved canonical SHA-256."
                )
            resource_provenance = resolved_provenance
        elif not allow_network_access:
            raise RuntimeError(
                "Official CollecTRI loading requires direct network access in decoupler 2.2.0; "
                "enable allow_network_access or use a fingerprinted local snapshot."
            )
        else:
            loaded = decoupler.op.collectri(
                organism="human",
                remove_complexes=complex_policy == "remove",
                license=affiliation_license,
                verbose=False,
            )
            network, accounting = _normalize_network(
                loaded, pandas=pd, complex_policy="retain", official=True
            )
            resource_provenance = {
                "mode": "official_collectri",
                "requested_path": None,
                "resolved_path": "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1",
                "size_bytes": None,
                "file_sha256": None,
                "canonical_network_sha256": _canonical_network_sha256(network),
                "metadata": {
                    "name": "CollecTRI",
                    "version": "2.0",
                    "date": "2023-03-29",
                    "organism": "human",
                    "identifier_namespace": "human gene symbols",
                    "scope": "post-decoupler-2.2.0 human table",
                    "license": "CC BY 4.0; incorporated-source licenses also apply",
                    "citation": "Müller-Dott et al. 2023; doi:10.1093/nar/gkad841",
                },
                "affiliation_license_declaration": affiliation_license,
                "affiliation_filter_enforced_by_decoupler_2_2": False,
                "network_access": True,
                "cache_state": "not_supported_by_decoupler_2_2",
                "accounting": accounting,
            }
    if expected_canonical_network_sha256 is not None:
        observed_hash = resource_provenance["canonical_network_sha256"]
        if observed_hash != expected_canonical_network_sha256:
            raise ValueError(
                "CollecTRI canonical network fingerprint changed: "
                f"expected {expected_canonical_network_sha256}, observed {observed_hash}."
            )
    feature_set = set(features)
    overlap = network.loc[network["target"].isin(feature_set)].copy()
    before_counts = network.groupby("source", sort=True, observed=True)["target"].nunique()
    overlap_counts = overlap.groupby("source", sort=True, observed=True)["target"].nunique()
    regulators = sorted(
        source for source, count in overlap_counts.items() if int(count) >= min_targets
    )
    if not regulators:
        raise ValueError(
            f"CollecTRI has no regulator with at least {min_targets} exact targets on the selected feature axis."
        )
    overlap = overlap.loc[overlap["source"].isin(regulators)].copy()
    overlap = overlap.sort_values(["source", "target"], kind="mergesort", ignore_index=True)
    output_rows = len(observations) * len(regulators)
    if output_rows > max_output_rows:
        raise ValueError(
            f"CollecTRI ULM would produce {output_rows:,} observation-regulator rows, exceeding "
            f"max_output_rows={max_output_rows:,}."
        )
    dense_expression_bytes = len(observations) * len(features) * 8
    memory_components = {
        "dense_expression_copy_bytes": dense_expression_bytes,
        "dense_adjacency_bytes": len(features) * len(regulators) * 8,
        "score_raw_pvalue_adjusted_peak_bytes": output_rows * 8 * 3,
        "canonical_output_frame_bytes": output_rows * 8 * 2,
        "bounded_result_ordering_index_bytes": output_rows * 8,
    }
    estimated_working_bytes = sum(memory_components.values())
    if estimated_working_bytes > int(max_working_memory_gib * 1024**3):
        raise MemoryError(
            f"CollecTRI ULM estimates {estimated_working_bytes / 1024**3:.3f} GiB working memory, "
            f"exceeding max_working_memory_gib={max_working_memory_gib:.3f}."
        )
    adjacency = np.zeros((len(features), len(regulators)), dtype=float)
    feature_positions = {feature: index for index, feature in enumerate(features)}
    regulator_positions = {regulator: index for index, regulator in enumerate(regulators)}
    for source, target, weight in overlap.itertuples(index=False, name=None):
        adjacency[feature_positions[target], regulator_positions[source]] = float(weight)
    adjacency_variance = np.var(adjacency, axis=0)
    if not bool((adjacency_variance > 0.0).all()):
        invalid = [regulators[index] for index in np.flatnonzero(adjacency_variance <= 0.0)]
        raise ValueError(
            "CollecTRI retained regulators have zero-variance complete feature-axis weights: "
            f"{invalid[:10]!r}."
        )
    dense_for_checks = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix, dtype=float)
    work = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=pd.Index(observations)),
        var=pd.DataFrame(index=pd.Index(features)),
    )
    work_expression_before = _matrix_fingerprint(
        work.X, observations=observations, features=features, numpy=np, sparse=sparse
    )
    work_obs_columns_before = list(work.obs.columns)
    work_var_columns_before = list(work.var.columns)
    work_layer_keys_before = list(work.layers.keys())
    work_uns_before = copy.deepcopy(dict(work.uns))
    work_obsp_keys_before = list(work.obsp.keys())
    work_varm_keys_before = list(work.varm.keys())
    work_varp_keys_before = list(work.varp.keys())
    work_raw_before = work.raw
    network_before = network.copy(deep=True)
    backend_return = decoupler.mt.ulm(
        data=work,
        net=network,
        tmin=min_targets,
        raw=False,
        empty=False,
        bsize=batch_size,
        verbose=False,
        tval=True,
    )
    if backend_return is not None:
        raise RuntimeError("decoupler 2.2.0 mt.ulm must mutate its private AnnData and return None.")
    if _matrix_fingerprint(
        work.X, observations=observations, features=features, numpy=np, sparse=sparse
    ) != work_expression_before:
        raise RuntimeError("decoupler mt.ulm modified the private expression matrix.")
    if list(work.obs_names) != observations or list(work.var_names) != features:
        raise RuntimeError("decoupler mt.ulm modified private AnnData axes.")
    if (
        list(work.obs.columns) != work_obs_columns_before
        or list(work.var.columns) != work_var_columns_before
        or list(work.layers.keys()) != work_layer_keys_before
        or dict(work.uns) != work_uns_before
        or list(work.obsp.keys()) != work_obsp_keys_before
        or list(work.varm.keys()) != work_varm_keys_before
        or list(work.varp.keys()) != work_varp_keys_before
        or work.raw is not work_raw_before
    ):
        raise RuntimeError(
            "decoupler mt.ulm modified unsupported private AnnData state: "
            f"obs={list(work.obs.columns)!r}, var={list(work.var.columns)!r}, "
            f"layers={list(work.layers.keys())!r}, uns={list(work.uns)!r}."
        )
    if set(work.obsm.keys()) != {"score_ulm", "padj_ulm"}:
        raise RuntimeError(
            "decoupler mt.ulm output keys changed; expected exactly score_ulm and padj_ulm."
        )
    try:
        pd.testing.assert_frame_equal(network, network_before, check_exact=True)
    except AssertionError as exc:
        raise RuntimeError("decoupler mt.ulm modified the caller-owned network table.") from exc
    scores = work.obsm["score_ulm"]
    adjusted = work.obsm["padj_ulm"]
    for frame, label in ((scores, "score_ulm"), (adjusted, "padj_ulm")):
        if not isinstance(frame, pd.DataFrame):
            raise RuntimeError(f"decoupler mt.ulm {label} output must be a pandas DataFrame.")
        if frame.index.tolist() != observations or frame.columns.tolist() != regulators:
            raise RuntimeError(f"decoupler mt.ulm {label} axes differ from the complete canonical family.")
        values = frame.to_numpy(dtype=float, copy=True)
        if not bool(np.isfinite(values).all()):
            raise RuntimeError(f"decoupler mt.ulm {label} contains non-finite values.")
    adjusted_values = adjusted.to_numpy(dtype=float, copy=True)
    if bool(((adjusted_values < 0.0) | (adjusted_values > 1.0)).any()):
        raise RuntimeError("decoupler mt.ulm adjusted p-values fall outside [0, 1].")
    expression_dense = np.asarray(dense_for_checks, dtype=float)
    centered_expression = expression_dense - expression_dense.mean(axis=1, keepdims=True)
    centered_adjacency = adjacency - adjacency.mean(axis=0, keepdims=True)
    covariance = centered_expression @ centered_adjacency / (len(features) - 1)
    expression_sd = expression_dense.std(axis=1, ddof=1)[:, None]
    adjacency_sd = adjacency.std(axis=0, ddof=1)[None, :]
    correlation = covariance / (expression_sd * adjacency_sd)
    degrees_freedom = len(features) - 2
    expected_scores = correlation * np.sqrt(
        degrees_freedom / ((1.0 - correlation + 2.2e-16) * (1.0 + correlation + 2.2e-16))
    )
    observed_scores = scores.to_numpy(dtype=float, copy=True)
    if not bool(np.allclose(observed_scores, expected_scores, rtol=2e-5, atol=2e-5)):
        raise RuntimeError("decoupler mt.ulm scores differ from independent complete-axis ULM calculation.")
    raw_pvalues = scipy.stats.t.sf(np.abs(expected_scores), degrees_freedom) * 2.0
    expected_adjusted = np.vstack([_bh_adjust(row, numpy=np) for row in raw_pvalues])
    if not bool(np.allclose(adjusted_values, expected_adjusted, rtol=5e-5, atol=5e-6)):
        raise RuntimeError("decoupler mt.ulm adjusted p-values differ from independent row-wise BH calculation.")
    network_accounting = []
    for regulator in sorted(before_counts.index.tolist()):
        network_accounting.append(
            {
                "regulator": regulator,
                "targets_before_feature_intersection": int(before_counts.loc[regulator]),
                "targets_on_feature_axis": int(overlap_counts.get(regulator, 0)),
                "retained_by_min_targets": regulator in set(regulators),
            }
        )
    expression_provenance = {
        "source": source_kind,
        "layer_name": layer_name if source_kind == "layer" else None,
        "observations": len(observations),
        "features": len(features),
        "observation_axis_sha256": observation_axis_sha256,
        "feature_axis_sha256": feature_axis_sha256,
        "expression_content_sha256": expression_sha256,
        "value_min": value_min,
        "value_max": value_max,
        "integer_like": integer_like,
        "sparse": bool(sparse.issparse(matrix)),
    }
    parameters = {
        "resource_mode": resource_mode,
        "source": source_kind,
        "layer_name": layer_name if source_kind == "layer" else None,
        "affiliation_license": affiliation_license,
        "allow_network_access": allow_network_access,
        "complex_policy": complex_policy,
        "min_targets": min_targets,
        "batch_size": batch_size,
        "max_output_rows": max_output_rows,
        "max_working_memory_gib": max_working_memory_gib,
        "overwrite_existing": overwrite_existing,
        "raw": False,
        "empty": False,
        "verbose": False,
        "tval": True,
    }
    provenance = {
        "method": "decoupler 2.2.0 ULM t statistic",
        "expression": expression_provenance,
        "resource": resource_provenance,
        "parameters": parameters,
        "retained_regulators": regulators,
        "retained_regulator_count": len(regulators),
        "retained_target_count": int(overlap["target"].nunique()),
        "retained_edge_count": int(len(overlap)),
        "network_accounting": network_accounting,
        "degrees_freedom": degrees_freedom,
        "multiple_testing_family": "within each observation across all retained regulators",
    }
    artifact = build_tf_activity_artifact(
        scores=scores,
        adjusted_pvalues=adjusted,
        provenance=provenance,
        numpy=np,
        pandas=pd,
        copy_frames=False,
    )
    output = adata
    collision = (
        COLLECTRI_SCORE_KEY in output.obsm
        or COLLECTRI_PADJ_KEY in output.obsm
        or COLLECTRI_UNS_KEY in output.uns
    )
    if collision and not overwrite_existing:
        raise ValueError(
            "CollecTRI owned output keys already exist; enable overwrite_existing explicitly to replace them."
        )
    output.obsm[COLLECTRI_SCORE_KEY] = scores
    output.obsm[COLLECTRI_PADJ_KEY] = adjusted
    output.uns[COLLECTRI_UNS_KEY] = {
        "artifact_fingerprint_sha256": artifact.fingerprint,
        "metadata": artifact.metadata,
        "provenance": artifact.provenance,
    }
    score_values = scores.to_numpy(dtype=float)
    strongest_positive = []
    strongest_negative = []
    for cell_position, regulator_position in np.dstack(
        np.unravel_index(np.argsort(score_values.ravel())[::-1], score_values.shape)
    )[0][:10]:
        strongest_positive.append(
            {
                "observation": observations[int(cell_position)],
                "regulator": regulators[int(regulator_position)],
                "score": float(score_values[cell_position, regulator_position]),
                "p_adjusted": float(adjusted_values[cell_position, regulator_position]),
            }
        )
    for cell_position, regulator_position in np.dstack(
        np.unravel_index(np.argsort(score_values.ravel()), score_values.shape)
    )[0][:10]:
        strongest_negative.append(
            {
                "observation": observations[int(cell_position)],
                "regulator": regulators[int(regulator_position)],
                "score": float(score_values[cell_position, regulator_position]),
                "p_adjusted": float(adjusted_values[cell_position, regulator_position]),
            }
        )
    warnings.extend(
        [
            "TF activities are exploratory cell-level scores; cells from one biological Sample are not independent replicates.",
            "ULM adjusted p-values test within-cell feature-weight regressions and are not Condition-level p-values.",
        ]
    )
    if resource_mode == "official_collectri":
        warnings.append(
            "decoupler 2.2.0 accepts the affiliation license argument but does not use it to filter the downloaded network."
        )
    software_packages = ["anndata", "decoupler", "numba", "numpy", "pandas", "scipy"]
    if resource_mode == "official_collectri":
        software_packages.append("requests")
    summary = {
        "schema_version": COLLECTRI_SUMMARY_SCHEMA,
        "node_id": "OpenBioSingleCellCollecTRIULM",
        "status": "exploratory_cell_level_tf_activity",
        "methods": (
            f"decoupler {DECOUPLER_VERSION} ULM regressed each observation's complete {len(features):,}-feature "
            f"axis on each signed CollecTRI weight vector ({degrees_freedom:,} df). Two-sided t p-values were "
            "Benjamini-Hochberg adjusted separately within each observation across all retained regulators."
        ),
        "results": (
            f"ULM produced finite activity scores for {len(observations):,} observations and "
            f"{len(regulators):,} regulators from {len(overlap):,} retained exact target edges."
        ),
        "key_results": {
            "observations": len(observations),
            "features": len(features),
            "regulators": len(regulators),
            "retained_targets": int(overlap["target"].nunique()),
            "retained_edges": int(len(overlap)),
            "degrees_freedom": degrees_freedom,
            "score_range": [float(score_values.min()), float(score_values.max())],
            "adjusted_pvalue_range": [float(adjusted_values.min()), float(adjusted_values.max())],
            "strongest_positive": strongest_positive,
            "strongest_negative": strongest_negative,
            "expression": expression_provenance,
            "resource": resource_provenance,
            "network_accounting": network_accounting,
            "memory_guard": {
                "estimated_working_bytes": estimated_working_bytes,
                "components": memory_components,
                "dense_bsize_limitation": "decoupler 2.2.0 bsize does not batch dense input",
            },
            "artifact_fingerprint_sha256": artifact.fingerprint,
        },
        "parameters": parameters,
        "references": copy.deepcopy(COLLECTRI_REFERENCES),
        "software_versions": collect_software_versions(
            software_packages, openbio_version=openbio_version
        ),
        "warnings": warnings,
        "limitations": [
            "Activity signs describe concordance with a signed regulon and do not establish direct binding or causality.",
            "Background-feature selection, normalization, scaling, identifier namespace, and network coverage affect ULM scores.",
            "The current AnnData contract cannot independently prove organism or feature-identifier namespace identity.",
            "Condition and Technical-batch effects require a separate replicate-aware Sample-level analysis.",
        ],
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return output, artifact, summary


def collectri_ulm_code(
    *,
    function_name: str = "infer_collectri_ulm",
    parameters: Mapping[str, Any],
    resolved_resource_provenance: Mapping[str, Any] | None = None,
) -> str:
    """Return importable equivalent source pinned to all runtime-resolved parameters and hashes."""

    if not isinstance(function_name, str) or not function_name.isidentifier():
        raise ValueError("Generated CollecTRI function_name must be a Python identifier.")
    portable_parameters = dict(parameters)
    for fixed_parameter in ("raw", "empty", "verbose", "tval"):
        portable_parameters.pop(fixed_parameter, None)
    official = portable_parameters.get("resource_mode") == "official_collectri"
    if official:
        if not isinstance(resolved_resource_provenance, Mapping):
            raise ValueError(
                "Generated official CollecTRI code requires runtime-resolved resource provenance."
            )
        portable_parameters["_resolved_official_resource_provenance"] = copy.deepcopy(
            dict(resolved_resource_provenance)
        )
    helper_sources = "\n\n".join(
        dedent(inspect.getsource(helper)).strip()
        for helper in (
            _canonical_json_sha256,
            _canonical_axis,
            _activity_frame_fingerprint,
            TFActivityArtifact,
            build_tf_activity_artifact,
            _portable_tf_activity_payload,
            validate_tf_activity_artifact,
            _strict_axis,
            _strict_metadata,
            _file_sha256,
            _canonical_network_sha256,
            _normalize_network,
            load_local_collectri_network,
            _signature_shape,
            _validate_decoupler_ulm,
            _matrix_fingerprint,
            _bh_adjust,
            _package_version,
            collect_software_versions,
            run_collectri_ulm,
        )
    )
    argument_lines = [
        f"        {name}={value!r}," for name, value in portable_parameters.items()
    ]
    materialized_argument = ", materialized_network_csv" if official else ""
    if official:
        argument_lines.append("        _materialized_official_network_path=materialized_network_csv,")
    argument_lines.append("        decoupler_module=decoupler_module,")
    arguments = "\n".join(argument_lines)
    generated_doc = (
        "The caller must supply a local post-decoupler network CSV whose canonical SHA-256 matches the runtime "
        "resource; this equivalent function performs no CollecTRI download."
        if official
        else "The local network bytes and canonical network SHA-256 are pinned; no download is performed."
    )
    return f'''from __future__ import annotations

import copy
import csv
import hashlib
import importlib
import inspect
import json
import math
import os
from collections.abc import Mapping, Sequence
from typing import Any

PLUGIN_VERSION = {PLUGIN_VERSION!r}
DECOUPLER_VERSION = {DECOUPLER_VERSION!r}
COLLECTRI_SCORE_KEY = {COLLECTRI_SCORE_KEY!r}
COLLECTRI_PADJ_KEY = {COLLECTRI_PADJ_KEY!r}
COLLECTRI_UNS_KEY = {COLLECTRI_UNS_KEY!r}
COLLECTRI_SUMMARY_SCHEMA = {COLLECTRI_SUMMARY_SCHEMA!r}
COLLECTRI_RESOURCE_METADATA_FIELDS = {COLLECTRI_RESOURCE_METADATA_FIELDS!r}
COLLECTRI_REFERENCES = {COLLECTRI_REFERENCES!r}
TF_ACTIVITY_ARTIFACT_TYPE = {TF_ACTIVITY_ARTIFACT_TYPE!r}
TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION = {TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION!r}
TF_ACTIVITY_PRODUCER_NODE_ID = {TF_ACTIVITY_PRODUCER_NODE_ID!r}
TF_ACTIVITY_PRODUCER_SCHEMA = {TF_ACTIVITY_PRODUCER_SCHEMA!r}

{helper_sources}


def {function_name}(adata{materialized_argument}, decoupler_module=None):
    """Return (AnnData copy, immutable TF activity artifact, strict summary).

    {generated_doc}
    """
    return run_collectri_ulm(
        adata,
{arguments}
    )
'''


__all__ = [
    "COLLECTRI_PADJ_KEY",
    "COLLECTRI_SCORE_KEY",
    "COLLECTRI_SUMMARY_SCHEMA",
    "COLLECTRI_UNS_KEY",
    "DECOUPLER_VERSION",
    "collectri_resource_cache_fingerprint",
    "collectri_ulm_code",
    "load_local_collectri_network",
    "run_collectri_ulm",
]
