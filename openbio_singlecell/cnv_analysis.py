from __future__ import annotations

import copy
import hashlib
import inspect
import json
import math
import platform
import re
import textwrap
import threading
import warnings
from collections.abc import Mapping, Sequence
from importlib import metadata
from typing import Any

from . import PLUGIN_VERSION

CNV_STATE_TYPE = "OPENBIO_CNV_STATE"
CNV_STATE_SCHEMA = "openbio-singlecell/cnv-state/v1"
INFER_CNV_NODE_ID = "OpenBioSingleCellInferCNV"
CNV_PCA_NODE_ID = "OpenBioSingleCellCNVPCA"
CNV_SCORE_NODE_ID = "OpenBioSingleCellCNVScore"
_CNV_PROVENANCE_KEY = "openbio_cnv_state"
_CNV_PCA_PROVENANCE_KEY = "openbio_cnv_pca"
_CNV_SCORE_PROVENANCE_KEY = "openbio_cnv_score"
_CNV_RUNTIME_LOCK = threading.RLock()


def _cnv_imports(operation):
    import importlib
    import re
    from importlib import metadata as package_metadata

    packages = {}
    for name in ("infercnvpy", "anndata", "numpy", "pandas", "scipy.sparse", "natsort"):
        try:
            packages[name] = importlib.import_module(name)
        except (ImportError, OSError) as error:
            raise RuntimeError(f"{operation} requires infercnvpy 0.6.1 and its scientific dependencies.") from error
    try:
        version = package_metadata.version("infercnvpy")
    except package_metadata.PackageNotFoundError as error:
        raise RuntimeError(f"{operation} could not determine the infercnvpy version.") from error
    match = re.match(r"^(\d+)\.(\d+)\.(\d+)", version)
    if match is None or tuple(int(value) for value in match.groups()) != (0, 6, 1):
        raise RuntimeError(f"{operation} supports the audited infercnvpy 0.6.1 API; found {version!r}.")
    return (
        packages["infercnvpy"],
        packages["anndata"],
        packages["numpy"],
        packages["pandas"],
        packages["scipy.sparse"],
        packages["natsort"],
        version,
    )


def _cnv_signature(callable_object, *, operation, expected_names):
    try:
        actual = tuple(inspect.signature(callable_object).parameters)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"{operation} could not inspect the infercnvpy backend signature.") from error
    if actual != tuple(expected_names):
        raise RuntimeError(
            f"{operation} requires the audited infercnvpy signature {tuple(expected_names)!r}; found {actual!r}."
        )


def _cnv_text(value, *, label):
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    if not value or value != value.strip():
        raise ValueError(f"{label} must be nonempty and free of surrounding whitespace.")
    return value


def _cnv_key(value, *, label, allow_x_prefix=False):
    value = _cnv_text(value, label=label)
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", value):
        raise ValueError(f"{label} must contain only letters, digits, and underscores and start with a letter.")
    if not allow_x_prefix and value.startswith("X_"):
        raise ValueError(f"{label} must be prefix-free; the matrix is stored under an automatic 'X_' prefix.")
    return value


def _cnv_int(value, *, label, minimum, maximum):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer.")
    if value < minimum or value > maximum:
        raise ValueError(f"{label} must be between {minimum} and {maximum}.")
    return int(value)


def _cnv_float(value, *, label, minimum, maximum=None):
    if isinstance(value, bool):
        raise TypeError(f"{label} must be numeric.")
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{label} must be numeric.") from error
    if not math.isfinite(value) or value < minimum or (maximum is not None and value > maximum):
        suffix = f" and at most {maximum}" if maximum is not None else ""
        raise ValueError(f"{label} must be finite, at least {minimum}{suffix}.")
    return value


def _cnv_bool(value, *, label):
    if not isinstance(value, bool):
        raise TypeError(f"{label} must be a boolean.")
    return value


def _cnv_json(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Strict JSON data cannot contain NaN or infinity.")
        return value
    if isinstance(value, Mapping):
        return {str(key): _cnv_json(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_cnv_json(item) for item in value]
    if hasattr(value, "tolist"):
        return _cnv_json(value.tolist())
    if hasattr(value, "item"):
        try:
            return _cnv_json(value.item())
        except (TypeError, ValueError):
            pass
    raise TypeError(f"Unsupported strict JSON value type: {type(value).__name__}.")


def _cnv_sha_json(value):
    payload = json.dumps(_cnv_json(value), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _cnv_axis_hash(values):
    digest = hashlib.sha256()
    for value in values:
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little", signed=False))
        digest.update(encoded)
    return digest.hexdigest()


def _cnv_matrix_hash(matrix, *, numpy, scipy_sparse):
    digest = hashlib.sha256()
    if scipy_sparse.issparse(matrix):
        values = scipy_sparse.csr_matrix(matrix, dtype=numpy.float64, copy=True)
        values.sum_duplicates()
        values.eliminate_zeros()
        values.sort_indices()
        digest.update(b"csr-f8")
        digest.update(numpy.asarray(values.shape, dtype="<i8").tobytes())
        for array, dtype in ((values.indptr, "<i8"), (values.indices, "<i8"), (values.data, "<f8")):
            canonical = numpy.ascontiguousarray(array, dtype=dtype)
            digest.update(canonical.tobytes(order="C"))
    else:
        values = numpy.ascontiguousarray(matrix, dtype="<f8")
        digest.update(b"dense-f8")
        digest.update(numpy.asarray(values.shape, dtype="<i8").tobytes())
        digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _cnv_axes(adata, *, operation):
    try:
        if bool(adata.isbacked):
            raise ValueError(f"{operation} requires an in-memory AnnData; call adata.to_memory() first.")
        obs_names = list(adata.obs_names)
        var_names = list(adata.var_names)
    except AttributeError as error:
        raise TypeError(f"{operation} requires an AnnData input.") from error
    if not obs_names or not var_names:
        raise ValueError(f"{operation} requires nonempty observation and feature axes.")
    if not bool(adata.obs_names.is_unique) or not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique observation and feature identifiers.")
    for values, label in ((obs_names, "observation"), (var_names, "feature")):
        if any(not isinstance(value, str) or not value or value != value.strip() for value in values):
            raise ValueError(f"{operation} requires canonical nonempty string {label} identifiers.")
    return tuple(obs_names), tuple(var_names)


def _cnv_source_matrix(adata, *, source_kind, layer_name, operation, numpy, scipy_sparse):
    if source_kind not in {"X", "layer"}:
        raise ValueError(f"{operation} source must be 'X' or 'layer'; Raw is not supported.")
    if source_kind == "layer":
        layer_name = _cnv_key(layer_name, label=f"{operation} layer_name", allow_x_prefix=True)
        if layer_name not in adata.layers:
            raise ValueError(f"{operation} expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
    else:
        layer_name = None
        matrix = adata.X
    if tuple(matrix.shape) != tuple(adata.shape):
        raise ValueError(f"{operation} expression source must align to the complete AnnData axes.")
    values = matrix.data if scipy_sparse.issparse(matrix) else numpy.asarray(matrix).ravel()
    values = numpy.asarray(values)
    if values.size and (not numpy.issubdtype(values.dtype, numpy.number) or numpy.iscomplexobj(values)):
        raise TypeError(f"{operation} expression source must contain real numeric values.")
    if values.size and not bool(numpy.isfinite(values).all()):
        raise ValueError(f"{operation} expression source must contain finite real numeric values.")
    return matrix, layer_name


def _cnv_csv(value, *, label):
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a comma-separated string.")
    values = []
    for position, item in enumerate(value.split(",")):
        item = item.strip()
        if not item:
            raise ValueError(f"{label} contains a blank item at position {position}.")
        if item in values:
            raise ValueError(f"{label} contains duplicate item {item!r}.")
        values.append(item)
    return tuple(values)


def _cnv_canonical_string_series(series, *, label, pandas):
    result = []
    for position, value in enumerate(series.tolist()):
        missing = pandas.isna(value)
        if not isinstance(missing, bool):
            try:
                missing = bool(missing)
            except (TypeError, ValueError) as error:
                raise TypeError(f"{label} row {position} must contain a scalar value.") from error
        if missing:
            raise ValueError(f"{label} contains missing values.")
        if not isinstance(value, str) or not value or value != value.strip():
            raise TypeError(f"{label} must contain canonical nonempty string labels without coercion.")
        result.append(value)
    return tuple(result)


def _cnv_reference_design(
    adata,
    *,
    reference_key,
    reference_categories,
    sample_key,
    minimum_reference_cells,
    pandas,
):
    reference_key = _cnv_key(reference_key, label="Infer CNV reference_key", allow_x_prefix=True)
    sample_key = _cnv_key(sample_key, label="Infer CNV sample_key", allow_x_prefix=True)
    if reference_key == sample_key:
        raise ValueError("Infer CNV reference_key and biological sample_key must be different columns.")
    for key, role in ((reference_key, "reference"), (sample_key, "biological Sample")):
        if key not in adata.obs:
            raise ValueError(f"Infer CNV {role} column not found in obs[{key!r}].")
    if not isinstance(adata.obs[reference_key].dtype, pandas.CategoricalDtype):
        raise TypeError("Infer CNV reference annotations must be explicitly categorical.")
    references = _cnv_canonical_string_series(
        adata.obs[reference_key], label="Infer CNV reference annotations", pandas=pandas
    )
    samples = _cnv_canonical_string_series(
        adata.obs[sample_key], label="Infer CNV biological Sample annotations", pandas=pandas
    )
    requested = _cnv_csv(reference_categories, label="Infer CNV reference_categories")
    available = tuple(str(value) for value in adata.obs[reference_key].cat.categories)
    if any(not value or value != value.strip() for value in available):
        raise ValueError("Infer CNV reference categories contain blank or noncanonical labels.")
    missing = [value for value in requested if value not in available]
    if missing:
        raise ValueError(f"Infer CNV reference categories were not found exactly: {missing!r}.")
    mask = tuple(value in requested for value in references)
    reference_count = sum(mask)
    if reference_count < minimum_reference_cells:
        raise ValueError(
            f"Infer CNV requires at least {minimum_reference_cells} reference cells; found {reference_count}."
        )
    if reference_count == len(references):
        raise ValueError("Infer CNV requires at least one non-reference cell in addition to the background cells.")
    category_counts = [
        {"category": value, "cells": sum(label == value for label in references)} for value in requested
    ]
    if any(record["cells"] == 0 for record in category_counts):
        raise RuntimeError("Infer CNV reference category accounting is inconsistent.")
    sample_counts = []
    for sample in dict.fromkeys(samples):
        sample_counts.append(
            {
                "sample": sample,
                "reference_cells": sum(flag and sample_value == sample for flag, sample_value in zip(mask, samples, strict=True)),
                "total_cells": sum(sample_value == sample for sample_value in samples),
            }
        )
    return {
        "reference_key": reference_key,
        "reference_categories": requested,
        "sample_key": sample_key,
        "reference_mask": mask,
        "reference_cells": reference_count,
        "non_reference_cells": len(references) - reference_count,
        "category_counts": category_counts,
        "sample_counts": sample_counts,
        "reference_samples": sum(record["reference_cells"] > 0 for record in sample_counts),
        "reference_labels": references,
        "sample_labels": samples,
    }


def _cnv_coordinates(
    adata,
    *,
    genome_assembly,
    exclude_chromosomes,
    window_size,
    step,
    numpy,
    pandas,
    natsort,
):
    genome_assembly = _cnv_text(genome_assembly, label="Infer CNV genome_assembly")
    required = ("chromosome", "start", "end")
    missing = [column for column in required if column not in adata.var]
    if missing:
        raise ValueError(f"Infer CNV genomic coordinate columns are missing: {missing!r}.")
    chromosomes = _cnv_canonical_string_series(
        adata.var["chromosome"], label="Infer CNV chromosome annotations", pandas=pandas
    )
    if any(not chromosome.startswith("chr") for chromosome in chromosomes):
        raise ValueError("Infer CNV 0.6.1 requires chromosome labels beginning with 'chr'.")
    coordinates = {}
    for column in ("start", "end"):
        series = adata.var[column]
        if bool(series.isna().any()) or not pandas.api.types.is_integer_dtype(series.dtype):
            raise TypeError(f"Infer CNV var[{column!r}] must contain complete integer coordinates.")
        values = numpy.asarray(series, dtype=numpy.int64)
        if bool((values < 0).any()):
            raise ValueError(f"Infer CNV var[{column!r}] cannot contain negative coordinates.")
        coordinates[column] = values
    if bool((coordinates["start"] >= coordinates["end"]).any()):
        raise ValueError("Infer CNV genomic coordinates require start < end for every feature.")
    coordinate_identities = list(zip(chromosomes, coordinates["start"], coordinates["end"], strict=True))
    if len(set(coordinate_identities)) != len(coordinate_identities):
        raise ValueError("Infer CNV genomic coordinates contain duplicate chromosome/start/end identities.")
    requested_exclusions = _cnv_csv(exclude_chromosomes, label="Infer CNV exclude_chromosomes")
    effective_exclusions = tuple(dict.fromkeys([*requested_exclusions, "chrM"]))
    counts = {}
    for chromosome in chromosomes:
        counts[chromosome] = counts.get(chromosome, 0) + 1
    retained_counts = {key: value for key, value in counts.items() if key not in effective_exclusions}
    excluded_counts = {key: value for key, value in counts.items() if key in effective_exclusions}
    if not retained_counts:
        raise ValueError("Infer CNV chromosome exclusions remove every annotated feature.")
    too_short = {key: value for key, value in retained_counts.items() if value < window_size}
    if too_short:
        raise ValueError(
            f"Infer CNV retained chromosomes must each contain at least window_size={window_size} genes; "
            f"too short: {too_short!r}. Exclude those chromosomes or reduce the window deliberately."
        )
    ordered = list(natsort.natsorted(retained_counts))
    windows = []
    offset = 0
    for chromosome in ordered:
        gene_count = retained_counts[chromosome]
        window_count = int(math.ceil((gene_count - window_size + 1) / step))
        windows.append(
            {
                "chromosome": chromosome,
                "genes": gene_count,
                "first_window_index": offset,
                "window_count": window_count,
            }
        )
        offset += window_count
    if offset < 2:
        raise ValueError("Infer CNV requires at least two genomic windows after exclusions.")
    coordinate_records = [
        [gene, chromosome, int(start), int(end)]
        for gene, chromosome, start, end in zip(
            adata.var_names, chromosomes, coordinates["start"], coordinates["end"], strict=True
        )
    ]
    return {
        "genome_assembly": genome_assembly,
        "requested_exclusions": requested_exclusions,
        "effective_exclusions": effective_exclusions,
        "implicit_chrM_added": "chrM" not in requested_exclusions,
        "chromosome_counts": counts,
        "retained_chromosome_counts": retained_counts,
        "excluded_chromosome_counts": excluded_counts,
        "retained_genes": sum(retained_counts.values()),
        "excluded_genes": sum(excluded_counts.values()),
        "windows": windows,
        "n_windows": offset,
        "coordinate_fingerprint_sha256": _cnv_sha_json(coordinate_records),
    }


def _cnv_numeric_summary(values, *, numpy):
    if hasattr(values, "data") and not isinstance(values, numpy.ndarray):
        array = numpy.asarray(values.data, dtype=float)
    else:
        array = numpy.asarray(values, dtype=float).ravel()
    if not array.size:
        return {
            "n": 0,
            "min": None,
            "q1": None,
            "median": None,
            "mean": None,
            "q3": None,
            "max": None,
        }
    if not bool(numpy.isfinite(array).all()):
        raise ValueError("CNV numeric summary requires finite values.")
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


def _cnv_absolute_matrix_summary(matrix, *, numpy, scipy_sparse):
    if not scipy_sparse.issparse(matrix):
        return _cnv_numeric_summary(numpy.abs(numpy.asarray(matrix, dtype=float)), numpy=numpy)
    total = int(matrix.shape[0] * matrix.shape[1])
    if total < 1:
        return _cnv_numeric_summary([], numpy=numpy)
    stored = numpy.abs(numpy.asarray(matrix.data, dtype=float))
    if stored.size and not bool(numpy.isfinite(stored).all()):
        raise ValueError("CNV absolute matrix summary requires finite values.")
    positive = numpy.sort(stored[stored > 0])
    zero_count = total - int(positive.size)

    def quantile(probability):
        position = (total - 1) * probability
        lower = int(math.floor(position))
        upper = int(math.ceil(position))

        def value_at(index):
            return 0.0 if index < zero_count else float(positive[index - zero_count])

        lower_value = value_at(lower)
        upper_value = value_at(upper)
        return lower_value + (upper_value - lower_value) * (position - lower)

    return {
        "n": total,
        "min": float(quantile(0.0)),
        "q1": float(quantile(0.25)),
        "median": float(quantile(0.5)),
        "mean": float(positive.sum() / total),
        "q3": float(quantile(0.75)),
        "max": float(quantile(1.0)),
    }


def _cnv_versions(*, openbio_version):
    versions = {"python": platform.python_version(), "openbio-singlecell": str(openbio_version)}
    for package in ("infercnvpy", "scanpy", "anndata", "numpy", "pandas", "scipy", "natsort"):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _cnv_references(operation):
    infercnv = {
        "citation": "Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016;352:189-196.",
        "doi": "10.1126/science.aad0501",
        "url": "https://doi.org/10.1126/science.aad0501",
        "kind": "method",
    }
    software = {
        "citation": "infercnvpy 0.6.1 public API and source.",
        "doi": None,
        "url": "https://github.com/icbi-lab/infercnvpy/tree/v0.6.1",
        "kind": "software_documentation",
    }
    if operation == "infer":
        return [infercnv, software]
    if operation == "pca":
        return [
            {
                "citation": "Lehoucq RB, Sorensen DC, Yang C. ARPACK Users' Guide: Solution of Large-Scale Eigenvalue Problems with Implicitly Restarted Arnoldi Methods. SIAM. 1998.",
                "doi": "10.1137/1.9780898719628",
                "url": "https://doi.org/10.1137/1.9780898719628",
                "kind": "method",
            },
            software,
        ]
    if operation == "score":
        return [infercnv, software]
    raise ValueError(f"Unknown CNV reference family: {operation!r}.")


def _cnv_summary(*, node_id, methods, results, key_results, parameters, warnings_list, limitations, references, versions):
    summary = {
        "schema_version": 1,
        "node_id": node_id,
        "methods": str(methods),
        "results": str(results),
        "key_results": _cnv_json(key_results),
        "parameters": _cnv_json(parameters),
        "warnings": [str(value) for value in warnings_list],
        "limitations": [str(value) for value in limitations],
        "references": _cnv_json(references),
        "software_versions": _cnv_json(versions),
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def _cnv_warning_messages(records):
    result = []
    for record in records:
        rendered = f"{record.category.__name__}: {' '.join(str(record.message).split())}"
        if rendered not in result:
            result.append(rendered)
    return result


def _cnv_artifact_fingerprints(adata, artifact_metadata, *, numpy, scipy_sparse):
    obs_names, var_names = _cnv_axes(adata, operation="CNV state validation")
    source_kind = artifact_metadata.get("source_kind")
    layer_name = artifact_metadata.get("layer_name")
    if source_kind == "layer":
        if not isinstance(layer_name, str) or layer_name not in adata.layers:
            raise ValueError("CNV state source layer is missing.")
        source = adata.layers[layer_name]
    elif source_kind == "X":
        source = adata.X
    else:
        raise ValueError("CNV state source metadata is invalid.")
    output_key = artifact_metadata.get("output_key")
    if not isinstance(output_key, str) or f"X_{output_key}" not in adata.obsm:
        raise ValueError("CNV state output matrix is missing.")
    cnv_matrix = adata.obsm[f"X_{output_key}"]
    reference_key = artifact_metadata.get("reference_key")
    sample_key = artifact_metadata.get("sample_key")
    if reference_key not in adata.obs or sample_key not in adata.obs:
        raise ValueError("CNV state reference/Sample annotations are missing.")
    required_coordinates = ("chromosome", "start", "end")
    if any(column not in adata.var for column in required_coordinates):
        raise ValueError("CNV state genomic coordinates are missing.")
    coordinate_rows = [
        [str(gene), str(chromosome), int(start), int(end)]
        for gene, chromosome, start, end in zip(
            var_names,
            adata.var["chromosome"],
            adata.var["start"],
            adata.var["end"],
            strict=True,
        )
    ]
    design_rows = [
        [str(cell), str(reference), str(sample)]
        for cell, reference, sample in zip(
            obs_names, adata.obs[reference_key], adata.obs[sample_key], strict=True
        )
    ]
    provenance = adata.uns.get(_CNV_PROVENANCE_KEY)
    if not isinstance(provenance, Mapping):
        raise ValueError("CNV state OpenBio provenance is missing or malformed.")
    return {
        "observation_axis_sha256": _cnv_axis_hash(obs_names),
        "feature_axis_sha256": _cnv_axis_hash(var_names),
        "source_matrix_sha256": _cnv_matrix_hash(source, numpy=numpy, scipy_sparse=scipy_sparse),
        "coordinate_sha256": _cnv_sha_json(coordinate_rows),
        "reference_design_sha256": _cnv_sha_json(design_rows),
        "cnv_matrix_sha256": _cnv_matrix_hash(cnv_matrix, numpy=numpy, scipy_sparse=scipy_sparse),
        "cnv_metadata_sha256": _cnv_sha_json(adata.uns.get(output_key)),
        "openbio_provenance_sha256": _cnv_sha_json(provenance),
    }


def _cnv_artifact_hash(metadata_record, fingerprints):
    payload = dict(metadata_record)
    payload.pop("artifact_fingerprint_sha256", None)
    payload["fingerprints"] = dict(fingerprints)
    return _cnv_sha_json(payload)


class _StandaloneCNVState:
    __slots__ = ("_adata", "_metadata")
    artifact_type = CNV_STATE_TYPE

    def __init__(self, adata, metadata_record, *, _owned=False):
        object.__setattr__(self, "_adata", adata if _owned else adata.copy())
        object.__setattr__(self, "_metadata", copy.deepcopy(dict(metadata_record)))

    def __setattr__(self, name, value):
        raise AttributeError("CNVState is immutable; create a new validated artifact instead.")

    @property
    def fingerprint(self):
        return str(self._metadata["artifact_fingerprint_sha256"])

    @property
    def metadata(self):
        return copy.deepcopy(self._metadata)

    def to_adata(self):
        return self._adata.copy()

    def _owned_adata(self):
        """Return the private payload for one-shot worker encoding/consumption only."""

        return self._adata


CNVState = _StandaloneCNVState


def _cnv_validate_state(artifact, *, numpy, scipy_sparse, _owned=False):
    if getattr(artifact, "artifact_type", None) != CNV_STATE_TYPE:
        raise TypeError("CNV analysis requires an OPENBIO_CNV_STATE artifact.")
    if not callable(getattr(artifact, "to_adata", None)):
        raise TypeError("CNV state artifact must expose a defensive to_adata() copy.")
    metadata_record = getattr(artifact, "metadata", None)
    if not isinstance(metadata_record, Mapping):
        raise TypeError("CNV state artifact metadata must be a mapping.")
    metadata_record = copy.deepcopy(dict(metadata_record))
    required = {
        "schema",
        "artifact_type",
        "producer_node_id",
        "stage",
        "source_kind",
        "layer_name",
        "output_key",
        "reference_key",
        "sample_key",
        "genome_assembly",
        "parameters",
        "fingerprints",
        "artifact_fingerprint_sha256",
    }
    if set(metadata_record) != required:
        raise ValueError(f"CNV state metadata fields differ from the exact schema: {sorted(set(metadata_record) ^ required)!r}.")
    if (
        metadata_record["schema"] != CNV_STATE_SCHEMA
        or metadata_record["artifact_type"] != CNV_STATE_TYPE
        or metadata_record["producer_node_id"] != INFER_CNV_NODE_ID
        or metadata_record["stage"] != "inferred"
    ):
        raise ValueError("CNV state identity or stage is invalid.")
    adata = artifact._owned_adata() if _owned else artifact.to_adata()
    if not _owned and adata is artifact.to_adata():
        raise RuntimeError("CNV state to_adata() must return a new defensive copy.")
    observed = _cnv_artifact_fingerprints(adata, metadata_record, numpy=numpy, scipy_sparse=scipy_sparse)
    if observed != metadata_record["fingerprints"]:
        raise ValueError("CNV state artifact content fingerprint validation failed.")
    expected_hash = _cnv_artifact_hash(metadata_record, observed)
    if metadata_record["artifact_fingerprint_sha256"] != expected_hash or getattr(artifact, "fingerprint", None) != expected_hash:
        raise ValueError("CNV state artifact envelope fingerprint validation failed.")
    return adata, metadata_record


def _cnv_validate_output_matrix(matrix, *, n_obs, n_windows, operation, numpy, scipy_sparse):
    if tuple(getattr(matrix, "shape", ())) != (n_obs, n_windows):
        raise RuntimeError(f"{operation} backend returned a matrix with the wrong shape.")
    values = matrix.data if scipy_sparse.issparse(matrix) else numpy.asarray(matrix).ravel()
    values = numpy.asarray(values)
    if values.size and (not numpy.issubdtype(values.dtype, numpy.number) or numpy.iscomplexobj(values)):
        raise RuntimeError(f"{operation} backend returned a non-real-numeric matrix.")
    if values.size and not bool(numpy.isfinite(values).all()):
        raise RuntimeError(f"{operation} backend returned non-finite values.")
    return matrix


def _cnv_validate_pca_scores(cnv_matrix, scores, *, random_seed, numpy, scipy_sparse):
    """Independently verify that scores are the leading uncentered singular coordinates."""
    scores64 = numpy.asarray(scores, dtype=float)
    gram = scores64.T @ scores64
    eigenvalues = numpy.diag(gram)
    if not bool(numpy.isfinite(eigenvalues).all()) or bool((eigenvalues <= 0).any()):
        raise RuntimeError("CNV PCA backend returned invalid component energies.")
    scale = numpy.sqrt(numpy.outer(eigenvalues, eigenvalues))
    off_diagonal = gram - numpy.diag(eigenvalues)
    if bool((numpy.abs(off_diagonal) > (1e-5 * numpy.maximum(scale, 1.0))).any()):
        raise RuntimeError("CNV PCA backend scores are not mutually orthogonal.")
    if bool((eigenvalues[:-1] + 1e-7 * numpy.maximum(eigenvalues[:-1], 1.0) < eigenvalues[1:]).any()):
        raise RuntimeError("CNV PCA backend components are not ordered by decreasing singular value.")
    for index, eigenvalue in enumerate(eigenvalues):
        vector = scores64[:, index]
        right = cnv_matrix.T @ vector
        observed = cnv_matrix @ right
        observed = numpy.asarray(observed, dtype=float).reshape(-1)
        expected = float(eigenvalue) * vector
        denominator = max(float(numpy.linalg.norm(observed)), float(numpy.linalg.norm(expected)), 1.0)
        residual = float(numpy.linalg.norm(observed - expected) / denominator)
        if residual > 2e-4:
            raise RuntimeError(
                f"CNV PCA backend component {index + 1} failed the uncentered eigenpair residual check."
            )
    try:
        expected_singular_values = scipy_sparse.linalg.svds(
            cnv_matrix,
            k=scores64.shape[1],
            which="LM",
            solver="arpack",
            return_singular_vectors=False,
            random_state=random_seed,
        )
    except (ArithmeticError, RuntimeError, ValueError) as error:
        raise RuntimeError("CNV PCA independent leading-spectrum validation failed to converge.") from error
    expected_singular_values = numpy.sort(numpy.asarray(expected_singular_values, dtype=float))[::-1]
    observed_singular_values = numpy.sqrt(eigenvalues)
    tolerance = max(float(expected_singular_values[0]) * 2e-4, 1e-6)
    if not bool(
        numpy.allclose(
            observed_singular_values,
            expected_singular_values,
            rtol=2e-4,
            atol=tolerance,
        )
    ):
        raise RuntimeError("CNV PCA backend scores disagree with an independent leading-spectrum calculation.")
    return {
        "eigenpair_residual_tolerance": 2e-4,
        "spectrum_relative_tolerance": 2e-4,
        "observed_singular_values": [float(value) for value in observed_singular_values],
    }


def _cnv_run_infer(
    adata,
    *,
    source_kind="layer",
    layer_name="log1p_norm",
    reference_key="cell_type",
    reference_categories="",
    sample_key="sample",
    genome_assembly="",
    window_size=100,
    step=10,
    lfc_clip=3.0,
    dynamic_threshold=1.5,
    exclude_chromosomes="chrX,chrY,chrM",
    output_key="cnv",
    minimum_reference_cells=20,
    chunksize=5000,
    n_jobs=1,
    max_output_gib=4.0,
    overwrite_existing=False,
    openbio_version="unknown",
    _worker_owned=False,
):
    operation = "Infer CNV"
    infercnvpy, anndata, numpy, pandas, scipy_sparse, natsort, backend_version = _cnv_imports(operation)
    _cnv_signature(
        infercnvpy.tl.infercnv,
        operation=operation,
        expected_names=(
            "adata",
            "reference_key",
            "reference_cat",
            "reference",
            "lfc_clip",
            "window_size",
            "step",
            "dynamic_threshold",
            "exclude_chromosomes",
            "chunksize",
            "n_jobs",
            "inplace",
            "layer",
            "key_added",
            "calculate_gene_values",
        ),
    )
    obs_names, var_names = _cnv_axes(adata, operation=operation)
    if len(obs_names) < 2:
        raise ValueError("Infer CNV requires at least two cells.")
    window_size = _cnv_int(window_size, label="Infer CNV window_size", minimum=2, maximum=10000)
    step = _cnv_int(step, label="Infer CNV step", minimum=1, maximum=window_size)
    lfc_clip = _cnv_float(lfc_clip, label="Infer CNV lfc_clip", minimum=1e-12)
    dynamic_threshold = _cnv_float(
        dynamic_threshold, label="Infer CNV dynamic_threshold", minimum=0.0
    )
    minimum_reference_cells = _cnv_int(
        minimum_reference_cells,
        label="Infer CNV minimum_reference_cells",
        minimum=2,
        maximum=len(obs_names) - 1,
    )
    chunksize = _cnv_int(chunksize, label="Infer CNV chunksize", minimum=1, maximum=2**31 - 1)
    n_jobs = _cnv_int(n_jobs, label="Infer CNV n_jobs", minimum=1, maximum=256)
    max_output_gib = _cnv_float(max_output_gib, label="Infer CNV max_output_gib", minimum=0.01)
    overwrite_existing = _cnv_bool(overwrite_existing, label="Infer CNV overwrite_existing")
    output_key = _cnv_key(output_key, label="Infer CNV output_key")
    matrix, layer_name = _cnv_source_matrix(
        adata,
        source_kind=source_kind,
        layer_name=layer_name,
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    design = _cnv_reference_design(
        adata,
        reference_key=reference_key,
        reference_categories=reference_categories,
        sample_key=sample_key,
        minimum_reference_cells=minimum_reference_cells,
        pandas=pandas,
    )
    coordinates = _cnv_coordinates(
        adata,
        genome_assembly=genome_assembly,
        exclude_chromosomes=exclude_chromosomes,
        window_size=window_size,
        step=step,
        numpy=numpy,
        pandas=pandas,
        natsort=natsort,
    )
    matrix_key = f"X_{output_key}"
    collisions = []
    for container_name, container, key in (
        ("obsm", adata.obsm, matrix_key),
        ("uns", adata.uns, output_key),
        ("uns", adata.uns, _CNV_PROVENANCE_KEY),
    ):
        if key in container:
            collisions.append(f"{container_name}[{key!r}]")
    if collisions and not overwrite_existing:
        raise ValueError(f"Infer CNV output bundle already exists: {', '.join(collisions)}")
    dense_chunk_bytes = min(chunksize, len(obs_names)) * len(var_names) * 8
    dense_output_bytes = len(obs_names) * coordinates["n_windows"] * 8
    estimated_bytes = dense_chunk_bytes + dense_output_bytes
    if estimated_bytes > max_output_gib * 1024**3:
        raise ValueError(
            f"Infer CNV estimated dense working/output memory ({estimated_bytes / 1024**3:.3f} GiB) "
            "exceeds max_output_gib."
        )
    source_fingerprint = _cnv_matrix_hash(matrix, numpy=numpy, scipy_sparse=scipy_sparse)
    working_matrix = matrix.copy() if hasattr(matrix, "copy") else numpy.array(matrix, copy=True)
    working = anndata.AnnData(
        X=working_matrix,
        obs=adata.obs[[design["reference_key"], design["sample_key"]]].copy(),
        var=adata.var[["chromosome", "start", "end"]].copy(),
    )
    working.obs_names = pandas.Index(numpy.asarray(obs_names, dtype=object), dtype=object)
    # infercnvpy 0.6.1 indexes a 2D convolution array into var_names. Pandas 3's ArrowStringArray rejects that
    # operation, so the private backend adapter uses an equivalent object-string Index and restores canonical axes.
    working.var_names = pandas.Index(numpy.asarray(var_names, dtype=object), dtype=object)
    working_design_fingerprint = _cnv_sha_json(
        {
            "obs_names": list(obs_names),
            "var_names": list(var_names),
            "reference": design["reference_labels"],
            "samples": design["sample_labels"],
            "coordinates": [
                [str(chromosome), int(start), int(end)]
                for chromosome, start, end in zip(
                    working.var["chromosome"], working.var["start"], working.var["end"], strict=True
                )
            ],
        }
    )
    with _CNV_RUNTIME_LOCK:
        rng_state = numpy.random.get_state()
        print_options = numpy.get_printoptions()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                backend_return = infercnvpy.tl.infercnv(
                    working,
                    reference_key=design["reference_key"],
                    reference_cat=list(design["reference_categories"]),
                    reference=None,
                    lfc_clip=lfc_clip,
                    window_size=window_size,
                    step=step,
                    dynamic_threshold=dynamic_threshold,
                    exclude_chromosomes=list(coordinates["effective_exclusions"]),
                    chunksize=chunksize,
                    n_jobs=n_jobs,
                    inplace=True,
                    layer=None,
                    key_added=output_key,
                    calculate_gene_values=False,
                )
        finally:
            numpy.random.set_state(rng_state)
            numpy.set_printoptions(**print_options)
    if backend_return is not None:
        raise RuntimeError("Infer CNV backend violated inplace=True by returning a value.")
    if _cnv_matrix_hash(matrix, numpy=numpy, scipy_sparse=scipy_sparse) != source_fingerprint:
        raise RuntimeError("Infer CNV backend modified the selected input expression source.")
    if _cnv_matrix_hash(working.X, numpy=numpy, scipy_sparse=scipy_sparse) != source_fingerprint:
        raise RuntimeError("Infer CNV backend modified its private expression input.")
    observed_working_design_fingerprint = _cnv_sha_json(
        {
            "obs_names": list(working.obs_names),
            "var_names": list(working.var_names),
            "reference": tuple(str(value) for value in working.obs[design["reference_key"]].tolist()),
            "samples": tuple(str(value) for value in working.obs[design["sample_key"]].tolist()),
            "coordinates": [
                [str(chromosome), int(start), int(end)]
                for chromosome, start, end in zip(
                    working.var["chromosome"], working.var["start"], working.var["end"], strict=True
                )
            ],
        }
    )
    if observed_working_design_fingerprint != working_design_fingerprint:
        raise RuntimeError("Infer CNV backend modified its private axes, reference/Sample design, or coordinates.")
    cnv_matrix = _cnv_validate_output_matrix(
        working.obsm.get(matrix_key),
        n_obs=len(obs_names),
        n_windows=coordinates["n_windows"],
        operation=operation,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    backend_metadata = working.uns.get(output_key)
    if not isinstance(backend_metadata, Mapping) or set(backend_metadata) != {"chr_pos"}:
        raise RuntimeError("Infer CNV backend returned malformed genomic-window metadata.")
    chr_pos = backend_metadata.get("chr_pos")
    if not isinstance(chr_pos, Mapping):
        raise RuntimeError("Infer CNV backend returned malformed chromosome start positions.")
    expected_chr_pos = {record["chromosome"]: record["first_window_index"] for record in coordinates["windows"]}
    observed_chr_pos = {str(key): int(value) for key, value in chr_pos.items()}
    if observed_chr_pos != expected_chr_pos:
        raise RuntimeError(
            f"Infer CNV backend chromosome/window metadata mismatch: expected {expected_chr_pos!r}, "
            f"found {observed_chr_pos!r}."
        )
    output = adata if _worker_owned else adata.copy()
    if overwrite_existing:
        for container, key in (
            (output.obsm, matrix_key),
            (output.uns, output_key),
            (output.uns, _CNV_PROVENANCE_KEY),
        ):
            if key in container:
                del container[key]
    # The infercnvpy working AnnData is private and dead after this transfer;
    # an owned worker can move its matrix reference without another full output.
    output.obsm[matrix_key] = cnv_matrix if _worker_owned else cnv_matrix.copy()
    output.uns[output_key] = copy.deepcopy(dict(backend_metadata))
    parameters = {
        "source_kind": source_kind,
        "layer_name": layer_name,
        "reference_key": design["reference_key"],
        "reference_categories": list(design["reference_categories"]),
        "sample_key": design["sample_key"],
        "genome_assembly": coordinates["genome_assembly"],
        "window_size": window_size,
        "step": step,
        "lfc_clip": lfc_clip,
        "dynamic_threshold": dynamic_threshold,
        "exclude_chromosomes_requested": list(coordinates["requested_exclusions"]),
        "exclude_chromosomes_effective": list(coordinates["effective_exclusions"]),
        "output_key": output_key,
        "minimum_reference_cells": minimum_reference_cells,
        "chunksize": chunksize,
        "n_jobs": n_jobs,
        "max_output_gib": max_output_gib,
        "overwrite_existing": overwrite_existing,
        "calculate_gene_values": False,
        "inplace": True,
        "backend_layer": None,
    }
    metadata_record = {
        "schema": CNV_STATE_SCHEMA,
        "artifact_type": CNV_STATE_TYPE,
        "producer_node_id": INFER_CNV_NODE_ID,
        "stage": "inferred",
        "source_kind": source_kind,
        "layer_name": layer_name,
        "output_key": output_key,
        "reference_key": design["reference_key"],
        "sample_key": design["sample_key"],
        "genome_assembly": coordinates["genome_assembly"],
        "parameters": parameters,
        "fingerprints": {},
        "artifact_fingerprint_sha256": "",
    }
    output.uns[_CNV_PROVENANCE_KEY] = {
        "schema": CNV_STATE_SCHEMA,
        "output_key": output_key,
        "genome_assembly": coordinates["genome_assembly"],
        "reference_key": design["reference_key"],
        "reference_categories": list(design["reference_categories"]),
        "sample_key": design["sample_key"],
        "windows": coordinates["windows"],
        "infercnvpy_version": backend_version,
    }
    fingerprints = _cnv_artifact_fingerprints(output, metadata_record, numpy=numpy, scipy_sparse=scipy_sparse)
    metadata_record["fingerprints"] = fingerprints
    metadata_record["artifact_fingerprint_sha256"] = _cnv_artifact_hash(metadata_record, fingerprints)
    artifact = _StandaloneCNVState(output, metadata_record, _owned=_worker_owned)
    _cnv_validate_state(
        artifact, numpy=numpy, scipy_sparse=scipy_sparse, _owned=_worker_owned
    )
    report_warnings = _cnv_warning_messages(caught)
    if coordinates["implicit_chrM_added"]:
        report_warnings.append(
            "infercnvpy 0.6.1 always omits chrM; it was added to the effective exclusion set and disclosed."
        )
    if design["reference_samples"] < 2:
        report_warnings.append(
            "Reference cells occur in fewer than two biological Samples; reference-specific technical effects cannot be assessed."
        )
    stored_values = cnv_matrix.data if scipy_sparse.issparse(cnv_matrix) else numpy.asarray(cnv_matrix).ravel()
    stored_values = numpy.asarray(stored_values, dtype=float)
    nonzero_count = int(numpy.count_nonzero(stored_values))
    total_values = int(len(obs_names) * coordinates["n_windows"])
    if nonzero_count == 0:
        report_warnings.append("All inferred CNV window values are zero after dynamic thresholding.")
    versions = _cnv_versions(openbio_version=openbio_version)
    key_results = {
        "artifact": {
            "artifact_type": CNV_STATE_TYPE,
            "schema": CNV_STATE_SCHEMA,
            "stage": "inferred",
            "fingerprint_sha256": metadata_record["artifact_fingerprint_sha256"],
        },
        "source": {
            "kind": source_kind,
            "layer_name": layer_name,
            "cells": len(obs_names),
            "genes": len(var_names),
            "fingerprint_sha256": source_fingerprint,
        },
        "reference": {
            "key": design["reference_key"],
            "categories": list(design["reference_categories"]),
            "reference_cells": design["reference_cells"],
            "non_reference_cells": design["non_reference_cells"],
            "category_counts": design["category_counts"],
            "sample_key": design["sample_key"],
            "sample_counts": design["sample_counts"],
            "reference_samples": design["reference_samples"],
        },
        "coordinates": {
            key: value
            for key, value in coordinates.items()
            if key not in {"requested_exclusions", "effective_exclusions", "implicit_chrM_added"}
        },
        "cnv_matrix": {
            "key": matrix_key,
            "shape": [len(obs_names), coordinates["n_windows"]],
            "storage": "sparse" if scipy_sparse.issparse(cnv_matrix) else "dense",
            "dtype": str(cnv_matrix.dtype),
            "nonzero_values": nonzero_count,
            "zero_fraction": float(1.0 - nonzero_count / total_values),
            "stored_value_summary": _cnv_numeric_summary(stored_values, numpy=numpy),
            "stored_absolute_value_summary": _cnv_numeric_summary(numpy.abs(stored_values), numpy=numpy),
            "absolute_value_summary_all_windows": _cnv_absolute_matrix_summary(
                cnv_matrix,
                numpy=numpy,
                scipy_sparse=scipy_sparse,
            ),
            "fingerprint_sha256": fingerprints["cnv_matrix_sha256"],
            "window_metadata": expected_chr_pos,
        },
        "memory": {
            "estimated_dense_chunk_bytes": dense_chunk_bytes,
            "estimated_dense_output_bytes": dense_output_bytes,
            "estimated_total_bytes": estimated_bytes,
        },
    }
    summary = _cnv_summary(
        node_id=INFER_CNV_NODE_ID,
        methods=(
            f"Inferred expression-derived CNV windows with infercnvpy {backend_version} from the explicitly "
            f"selected current-axis expression source, using {design['reference_cells']} nominated reference cells, "
            f"window_size={window_size}, step={step}, lfc_clip={lfc_clip}, and dynamic_threshold={dynamic_threshold}."
        ),
        results=(
            f"Generated {coordinates['n_windows']} genomic windows across {len(coordinates['windows'])} retained "
            f"chromosomes for {len(obs_names)} cells; {nonzero_count:,} of {total_values:,} values were nonzero "
            "after thresholding."
        ),
        key_results=key_results,
        parameters=parameters,
        warnings_list=report_warnings,
        limitations=(
            "Expression-derived CNV is indirect, experimental evidence and does not replace DNA copy-number measurement.",
            "Results depend on cell identity, technical effects, reference choice, genome assembly, gene coverage, smoothing, and thresholding.",
            "The analysis does not infer a tumor label or cutoff and performs no Sample-level Condition test.",
        ),
        references=_cnv_references("infer"),
        versions=versions,
    )
    return artifact, summary


def _cnv_run_pca(
    cnv_state,
    *,
    n_comps=30,
    output_key="X_cnv_pca",
    overwrite_existing=False,
    max_output_gib=2.0,
    random_seed=0,
    openbio_version="unknown",
    _worker_owned=False,
):
    operation = "CNV PCA"
    infercnvpy, _anndata, numpy, _pandas, scipy_sparse, _natsort, backend_version = _cnv_imports(operation)
    _cnv_signature(
        infercnvpy.tl.pca,
        operation=operation,
        expected_names=("adata", "svd_solver", "zero_center", "inplace", "use_rep", "key_added", "kwargs"),
    )
    adata, state_metadata = _cnv_validate_state(
        cnv_state,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        _owned=_worker_owned,
    )
    obs_names, _var_names = _cnv_axes(adata, operation=operation)
    use_rep = state_metadata["output_key"]
    cnv_matrix = adata.obsm[f"X_{use_rep}"]
    n_windows = int(cnv_matrix.shape[1])
    n_comps = _cnv_int(n_comps, label="CNV PCA n_comps", minimum=1, maximum=4096)
    if n_comps >= min(len(obs_names), n_windows):
        raise ValueError("CNV PCA with ARPACK requires n_comps < min(n_obs, n_windows).")
    output_key = _cnv_key(output_key, label="CNV PCA output_key", allow_x_prefix=True)
    if not output_key.startswith("X_") or len(output_key) <= 2:
        raise ValueError("CNV PCA output_key must use the explicit 'X_' AnnData representation namespace.")
    backend_key = output_key[2:]
    backend_key = _cnv_key(backend_key, label="CNV PCA backend key")
    overwrite_existing = _cnv_bool(overwrite_existing, label="CNV PCA overwrite_existing")
    max_output_gib = _cnv_float(max_output_gib, label="CNV PCA max_output_gib", minimum=0.01)
    random_seed = _cnv_int(random_seed, label="CNV PCA random_seed", minimum=0, maximum=2**31 - 1)
    collisions = [
        f"obsm[{output_key!r}]" if output_key in adata.obsm else None,
        f"uns[{_CNV_PCA_PROVENANCE_KEY!r}]" if _CNV_PCA_PROVENANCE_KEY in adata.uns else None,
    ]
    collisions = [value for value in collisions if value is not None]
    if collisions and not overwrite_existing:
        raise ValueError(f"CNV PCA output bundle already exists: {', '.join(collisions)}")
    estimated_bytes = len(obs_names) * n_comps * 4
    if estimated_bytes > max_output_gib * 1024**3:
        raise ValueError(
            f"CNV PCA dense output estimate ({estimated_bytes / 1024**3:.3f} GiB) exceeds max_output_gib."
        )
    input_state_fingerprint = state_metadata["artifact_fingerprint_sha256"]
    input_cnv_fingerprint = state_metadata["fingerprints"]["cnv_matrix_sha256"]
    output = adata if _worker_owned else adata.copy()
    if overwrite_existing:
        for container, key in (
            (output.obsm, output_key),
            (output.uns, _CNV_PCA_PROVENANCE_KEY),
        ):
            if key in container:
                del container[key]
    with _CNV_RUNTIME_LOCK:
        rng_state = numpy.random.get_state()
        print_options = numpy.get_printoptions()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                backend_return = infercnvpy.tl.pca(
                    output,
                    svd_solver="arpack",
                    zero_center=False,
                    inplace=True,
                    use_rep=use_rep,
                    key_added=backend_key,
                    n_comps=n_comps,
                    random_state=random_seed,
                    dtype="float32",
                )
        finally:
            numpy.random.set_state(rng_state)
            numpy.set_printoptions(**print_options)
    if backend_return is not None:
        raise RuntimeError("CNV PCA backend violated inplace=True by returning a value.")
    after_fingerprints = _cnv_artifact_fingerprints(
        output,
        state_metadata,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    if after_fingerprints != state_metadata["fingerprints"]:
        raise RuntimeError("CNV PCA backend modified the validated upstream CNV state.")
    scores = numpy.asarray(output.obsm.get(output_key))
    if scores.shape != (len(obs_names), n_comps):
        raise RuntimeError("CNV PCA backend returned the wrong score shape.")
    if not numpy.issubdtype(scores.dtype, numpy.floating) or numpy.iscomplexobj(scores):
        raise RuntimeError("CNV PCA backend scores must contain real floating-point values.")
    if not bool(numpy.isfinite(scores).all()):
        raise RuntimeError("CNV PCA backend scores contain non-finite values.")
    component_variances = numpy.var(scores.astype(float), axis=0)
    if bool((component_variances <= 0).any()):
        raise RuntimeError("CNV PCA backend returned one or more constant components.")
    scientific_validation = _cnv_validate_pca_scores(
        cnv_matrix,
        scores,
        random_seed=random_seed,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    output_fingerprint = _cnv_matrix_hash(scores, numpy=numpy, scipy_sparse=scipy_sparse)
    provenance = {
        "schema": "openbio-singlecell/cnv-pca/v1",
        "input_cnv_state_fingerprint_sha256": input_state_fingerprint,
        "input_cnv_matrix_fingerprint_sha256": input_cnv_fingerprint,
        "use_rep": use_rep,
        "output_key": output_key,
        "n_comps": n_comps,
        "svd_solver": "arpack",
        "zero_center": False,
        "dtype": "float32",
        "random_seed": random_seed,
        "output_fingerprint_sha256": output_fingerprint,
        "infercnvpy_version": backend_version,
    }
    output.uns[_CNV_PCA_PROVENANCE_KEY] = provenance
    report_warnings = _cnv_warning_messages(caught)
    parameters = {
        "n_comps": n_comps,
        "output_key": output_key,
        "overwrite_existing": overwrite_existing,
        "max_output_gib": max_output_gib,
        "random_seed": random_seed,
        "use_rep": use_rep,
        "svd_solver": "arpack",
        "zero_center": False,
        "dtype": "float32",
        "inplace": True,
    }
    key_results = {
        "input_cnv_state": {
            "artifact_fingerprint_sha256": input_state_fingerprint,
            "cnv_matrix_fingerprint_sha256": input_cnv_fingerprint,
            "genome_assembly": state_metadata["genome_assembly"],
            "reference_key": state_metadata["reference_key"],
            "sample_key": state_metadata["sample_key"],
            "parameters": state_metadata["parameters"],
        },
        "output": {
            "key": output_key,
            "shape": [len(obs_names), n_comps],
            "dtype": str(scores.dtype),
            "component_variances": [float(value) for value in component_variances],
            "component_summaries": [
                {"component": index + 1, **_cnv_numeric_summary(scores[:, index], numpy=numpy)}
                for index in range(n_comps)
            ],
            "independent_scientific_validation": scientific_validation,
            "output_fingerprint_sha256": output_fingerprint,
            "provenance_key": _CNV_PCA_PROVENANCE_KEY,
        },
        "memory": {"estimated_dense_output_bytes": estimated_bytes},
    }
    summary = _cnv_summary(
        node_id=CNV_PCA_NODE_ID,
        methods=(
            f"Reduced the validated {len(obs_names)} by {n_windows} inferred-CNV window matrix with "
            f"infercnvpy.tl.pca {backend_version}, n_comps={n_comps}, svd_solver='arpack', "
            "zero_center=False, and dtype='float32'."
        ),
        results=(
            f"Generated {n_comps} nonconstant CNV principal components for {len(obs_names)} cells under "
            f"{output_key!r}; all score and input-state fingerprints were verified."
        ),
        key_results=key_results,
        parameters=parameters,
        warnings_list=report_warnings,
        limitations=(
            "CNV principal components are exploratory and their signs are arbitrary.",
            "The representation is conditional on the upstream reference, genomic coordinates, smoothing, and thresholding.",
            "This node does not construct a graph, cluster cells, embed UMAP, score groups, classify tumors, or test Conditions.",
        ),
        references=_cnv_references("pca"),
        versions=_cnv_versions(openbio_version=openbio_version),
    )
    return output, summary


def _cnv_partition(adata, *, groupby, pandas, numpy):
    groupby = _cnv_key(groupby, label="CNV Score groupby", allow_x_prefix=True)
    if groupby not in adata.obs:
        raise ValueError(f"CNV Score groupby column not found in obs[{groupby!r}].")
    series = adata.obs[groupby]
    if not isinstance(series.dtype, pandas.CategoricalDtype):
        raise TypeError("CNV Score groupby must be explicitly categorical.")
    if bool(series.isna().any()):
        raise ValueError("CNV Score groupby contains missing labels.")
    categories = []
    displays = []
    counts = series.value_counts(sort=False)
    unused = []
    for index, value in enumerate(series.cat.categories):
        if hasattr(value, "item"):
            value = value.item()
        if isinstance(value, bool):
            raise TypeError("CNV Score boolean group labels are not supported.")
        if isinstance(value, str):
            if not value or value != value.strip():
                raise ValueError("CNV Score string group labels must be canonical and nonempty.")
            record = {"type": "string", "value": value, "display": value}
        elif isinstance(value, int):
            record = {"type": "integer", "value": int(value), "display": str(value)}
        elif isinstance(value, float) and bool(numpy.isfinite(value)):
            record = {"type": "number", "value": float(value), "display": repr(float(value))}
        else:
            raise TypeError("CNV Score group labels must be strings, integers, or finite numbers.")
        if record["display"] in displays:
            raise ValueError("CNV Score group categories collide after display rendering.")
        displays.append(record["display"])
        if int(counts.iloc[index]) > 0:
            categories.append(record)
        else:
            unused.append(record)
    if len(categories) < 2:
        raise ValueError("CNV Score requires at least two represented groups.")
    observed_records = []
    for value in series.astype(object).tolist():
        if hasattr(value, "item"):
            value = value.item()
        match = next(
            (
                record
                for record in categories
                if (
                    (record["type"] == "string" and isinstance(value, str) and record["value"] == value)
                    or (record["type"] == "integer" and isinstance(value, int) and not isinstance(value, bool) and record["value"] == value)
                    or (record["type"] == "number" and isinstance(value, float) and record["value"] == value)
                )
            ),
            None,
        )
        if match is None:
            raise RuntimeError("CNV Score categorical observation/category alignment is inconsistent.")
        observed_records.append(match)
    fingerprint = _cnv_sha_json({"groupby": groupby, "categories": categories, "assignments": observed_records})
    return {
        "groupby": groupby,
        "categories": categories,
        "unused_categories": unused,
        "assignments": observed_records,
        "fingerprint_sha256": fingerprint,
    }


def _cnv_run_score(
    cnv_state,
    adata,
    *,
    groupby,
    output_key="cnv_score",
    overwrite_existing=False,
    openbio_version="unknown",
    _worker_owned=False,
):
    operation = "CNV Score"
    infercnvpy, _anndata, numpy, pandas, scipy_sparse, _natsort, backend_version = _cnv_imports(operation)
    _cnv_signature(
        infercnvpy.tl.cnv_score,
        operation=operation,
        expected_names=("adata", "groupby", "use_rep", "key_added", "inplace", "obs_key"),
    )
    state_adata, state_metadata = _cnv_validate_state(
        cnv_state,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
        _owned=_worker_owned,
    )
    state_obs_names, state_var_names = _cnv_axes(state_adata, operation=operation)
    downstream_obs_names, downstream_var_names = _cnv_axes(adata, operation=operation)
    if downstream_obs_names != state_obs_names or downstream_var_names != state_var_names:
        raise ValueError("CNV Score downstream AnnData axes do not exactly match the typed CNV state.")
    downstream_fingerprints = _cnv_artifact_fingerprints(
        adata,
        state_metadata,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    if downstream_fingerprints != state_metadata["fingerprints"]:
        raise ValueError("CNV Score downstream AnnData content does not match the typed CNV state fingerprints.")
    partition = _cnv_partition(adata, groupby=groupby, pandas=pandas, numpy=numpy)
    output_key = _cnv_key(output_key, label="CNV Score output_key", allow_x_prefix=True)
    overwrite_existing = _cnv_bool(overwrite_existing, label="CNV Score overwrite_existing")
    collisions = [
        f"obs[{output_key!r}]" if output_key in adata.obs else None,
        f"uns[{_CNV_SCORE_PROVENANCE_KEY!r}]" if _CNV_SCORE_PROVENANCE_KEY in adata.uns else None,
    ]
    collisions = [value for value in collisions if value is not None]
    if collisions and not overwrite_existing:
        raise ValueError(f"CNV Score output bundle already exists: {', '.join(collisions)}")
    output = adata if _worker_owned else adata.copy()
    if partition["unused_categories"]:
        output.obs[partition["groupby"]] = output.obs[partition["groupby"]].cat.remove_unused_categories()
        retained_unused = partition["unused_categories"]
        partition = _cnv_partition(output, groupby=groupby, pandas=pandas, numpy=numpy)
        partition["unused_categories"] = retained_unused
    if overwrite_existing:
        for container, key in (
            (output.obs, output_key),
            (output.uns, _CNV_SCORE_PROVENANCE_KEY),
        ):
            if key in container:
                del container[key]
    use_rep = state_metadata["output_key"]
    cnv_matrix = output.obsm[f"X_{use_rep}"]
    input_cnv_fingerprint = state_metadata["fingerprints"]["cnv_matrix_sha256"]
    with _CNV_RUNTIME_LOCK:
        rng_state = numpy.random.get_state()
        print_options = numpy.get_printoptions()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                backend_return = infercnvpy.tl.cnv_score(
                    output,
                    groupby=partition["groupby"],
                    use_rep=use_rep,
                    key_added=output_key,
                    inplace=True,
                    obs_key=None,
                )
        finally:
            numpy.random.set_state(rng_state)
            numpy.set_printoptions(**print_options)
    if backend_return is not None:
        raise RuntimeError("CNV Score backend violated inplace=True by returning a value.")
    after_fingerprints = _cnv_artifact_fingerprints(
        output,
        state_metadata,
        numpy=numpy,
        scipy_sparse=scipy_sparse,
    )
    if after_fingerprints != state_metadata["fingerprints"]:
        raise RuntimeError("CNV Score backend modified the validated upstream CNV state.")
    scores = numpy.asarray(output.obs.get(output_key))
    if scores.shape != (output.n_obs,) or not numpy.issubdtype(scores.dtype, numpy.floating):
        raise RuntimeError("CNV Score backend returned an invalid per-cell score column.")
    if not bool(numpy.isfinite(scores).all()) or bool((scores < 0).any()):
        raise RuntimeError("CNV Score backend returned non-finite or negative scores.")
    table_rows = []
    for category in partition["categories"]:
        mask = numpy.asarray([record == category for record in partition["assignments"]], dtype=bool)
        count = int(mask.sum())
        selected = cnv_matrix[mask, :]
        if scipy_sparse.issparse(selected):
            absolute_sum = float(numpy.abs(selected).sum())
        else:
            dense = numpy.asarray(selected, dtype=float)
            absolute_sum = float(numpy.abs(dense).sum())
        absolute_summary = _cnv_absolute_matrix_summary(
            selected,
            numpy=numpy,
            scipy_sparse=scipy_sparse,
        )
        expected = absolute_sum / (count * int(cnv_matrix.shape[1]))
        observed = scores[mask]
        if not bool(numpy.allclose(observed, expected, rtol=1e-6, atol=1e-8)):
            raise RuntimeError("CNV Score backend disagrees with the independently recomputed group mean absolute CNV.")
        table_rows.append(
            {
                "group_type": category["type"],
                "group_value": category["value"],
                "group_display": category["display"],
                "cell_count": count,
                "cnv_score": float(expected),
                "absolute_value_summary": absolute_summary,
            }
        )
    table_rows.sort(key=lambda record: (-record["cnv_score"], record["group_display"]))
    table = pandas.DataFrame(
        [
            {
                "rank": index + 1,
                "group_type": record["group_type"],
                "group_value": record["group_value"],
                "group_display": record["group_display"],
                "cell_count": record["cell_count"],
                "cnv_score": record["cnv_score"],
                "absolute_min": record["absolute_value_summary"]["min"],
                "absolute_q1": record["absolute_value_summary"]["q1"],
                "absolute_median": record["absolute_value_summary"]["median"],
                "absolute_mean": record["absolute_value_summary"]["mean"],
                "absolute_q3": record["absolute_value_summary"]["q3"],
                "absolute_max": record["absolute_value_summary"]["max"],
            }
            for index, record in enumerate(table_rows)
        ],
        columns=(
            "rank",
            "group_type",
            "group_value",
            "group_display",
            "cell_count",
            "cnv_score",
            "absolute_min",
            "absolute_q1",
            "absolute_median",
            "absolute_mean",
            "absolute_q3",
            "absolute_max",
        ),
    )
    output_fingerprint = _cnv_matrix_hash(scores, numpy=numpy, scipy_sparse=scipy_sparse)
    provenance = {
        "schema": "openbio-singlecell/cnv-score/v1",
        "input_cnv_state_fingerprint_sha256": state_metadata["artifact_fingerprint_sha256"],
        "input_cnv_matrix_fingerprint_sha256": input_cnv_fingerprint,
        "partition_fingerprint_sha256": partition["fingerprint_sha256"],
        "groupby": partition["groupby"],
        "use_rep": use_rep,
        "output_key": output_key,
        "formula": "group mean(abs(inferred CNV)) across all member cells and windows",
        "independent_validation_rtol": 1e-6,
        "independent_validation_atol": 1e-8,
        "output_fingerprint_sha256": output_fingerprint,
        "infercnvpy_version": backend_version,
    }
    output.uns[_CNV_SCORE_PROVENANCE_KEY] = provenance
    report_warnings = _cnv_warning_messages(caught)
    if partition["unused_categories"]:
        report_warnings.append(
            f"Removed {len(partition['unused_categories'])} unused categorical level(s) on the output copy."
        )
    parameters = {
        "groupby": partition["groupby"],
        "use_rep": use_rep,
        "output_key": output_key,
        "overwrite_existing": overwrite_existing,
        "inplace": True,
        "obs_key": None,
        "independent_validation_rtol": 1e-6,
        "independent_validation_atol": 1e-8,
    }
    key_results = {
        "input_cnv_state": {
            "artifact_fingerprint_sha256": state_metadata["artifact_fingerprint_sha256"],
            "cnv_matrix_fingerprint_sha256": input_cnv_fingerprint,
            "genome_assembly": state_metadata["genome_assembly"],
            "reference_key": state_metadata["reference_key"],
            "sample_key": state_metadata["sample_key"],
            "parameters": state_metadata["parameters"],
            "downstream_adata_fingerprints_verified": True,
        },
        "partition": {
            "groupby": partition["groupby"],
            "represented_groups": len(partition["categories"]),
            "categories": partition["categories"],
            "unused_categories_removed": partition["unused_categories"],
            "fingerprint_sha256": partition["fingerprint_sha256"],
        },
        "group_scores": table_rows,
        "output": {
            "obs_key": output_key,
            "score_summary": _cnv_numeric_summary(scores, numpy=numpy),
            "output_fingerprint_sha256": output_fingerprint,
            "provenance_key": _CNV_SCORE_PROVENANCE_KEY,
        },
    }
    summary = _cnv_summary(
        node_id=CNV_SCORE_NODE_ID,
        methods=(
            f"Assigned infercnvpy.tl.cnv_score {backend_version} for each complete category in "
            f"obs[{partition['groupby']!r}] and independently verified mean(abs(X_{use_rep})) over all group cells "
            "and genomic windows."
        ),
        results=(
            f"Computed descriptive CNV scores for {len(table_rows)} groups across {output.n_obs} cells; "
            f"scores ranged from {float(table['cnv_score'].min()):.6g} to {float(table['cnv_score'].max()):.6g}."
        ),
        key_results=key_results,
        parameters=parameters,
        warnings_list=report_warnings,
        limitations=(
            "The score is a group mean and every member receives the same value; it is not an independent per-cell measurement.",
            "Scores are not absolute copy number, tumor probabilities, calibrated cutoffs, or hypothesis tests.",
            "Interpretation remains conditional on the grouping, reference, selected expression source, genome assembly, windows, and technical effects.",
            "No tumor label or Sample-level Condition inference was performed.",
        ),
        references=_cnv_references("score"),
        versions=_cnv_versions(openbio_version=openbio_version),
    )
    return output, table, summary


_CNV_INFER_SOURCE_HELPERS = (
    _cnv_imports,
    _cnv_signature,
    _cnv_text,
    _cnv_key,
    _cnv_int,
    _cnv_float,
    _cnv_bool,
    _cnv_json,
    _cnv_sha_json,
    _cnv_axis_hash,
    _cnv_matrix_hash,
    _cnv_axes,
    _cnv_source_matrix,
    _cnv_csv,
    _cnv_canonical_string_series,
    _cnv_reference_design,
    _cnv_coordinates,
    _cnv_numeric_summary,
    _cnv_absolute_matrix_summary,
    _cnv_versions,
    _cnv_references,
    _cnv_summary,
    _cnv_warning_messages,
    _cnv_artifact_fingerprints,
    _cnv_artifact_hash,
    _cnv_validate_state,
    _cnv_validate_output_matrix,
    _cnv_run_infer,
)

_CNV_PCA_SOURCE_HELPERS = (
    _cnv_imports,
    _cnv_signature,
    _cnv_text,
    _cnv_key,
    _cnv_int,
    _cnv_float,
    _cnv_bool,
    _cnv_json,
    _cnv_sha_json,
    _cnv_axis_hash,
    _cnv_matrix_hash,
    _cnv_axes,
    _cnv_numeric_summary,
    _cnv_versions,
    _cnv_references,
    _cnv_summary,
    _cnv_warning_messages,
    _cnv_artifact_fingerprints,
    _cnv_artifact_hash,
    _cnv_validate_state,
    _cnv_validate_pca_scores,
    _cnv_run_pca,
)

_CNV_SCORE_SOURCE_HELPERS = (
    _cnv_imports,
    _cnv_signature,
    _cnv_text,
    _cnv_key,
    _cnv_bool,
    _cnv_json,
    _cnv_sha_json,
    _cnv_axis_hash,
    _cnv_matrix_hash,
    _cnv_axes,
    _cnv_numeric_summary,
    _cnv_absolute_matrix_summary,
    _cnv_versions,
    _cnv_references,
    _cnv_summary,
    _cnv_warning_messages,
    _cnv_artifact_fingerprints,
    _cnv_artifact_hash,
    _cnv_validate_state,
    _cnv_partition,
    _cnv_run_score,
)


def analyze_infer_cnv(
    adata: Any,
    *,
    source_kind: str = "layer",
    layer_name: str | None = "log1p_norm",
    reference_key: str = "cell_type",
    reference_categories: str = "",
    sample_key: str = "sample",
    genome_assembly: str = "",
    window_size: int = 100,
    step: int = 10,
    lfc_clip: float = 3.0,
    dynamic_threshold: float = 1.5,
    exclude_chromosomes: str = "chrX,chrY,chrM",
    output_key: str = "cnv",
    minimum_reference_cells: int = 20,
    chunksize: int = 5000,
    n_jobs: int = 1,
    max_output_gib: float = 4.0,
    overwrite_existing: bool = False,
    _owned: bool = False,
) -> tuple[CNVState, dict[str, Any]]:
    return _cnv_run_infer(
        adata,
        source_kind=source_kind,
        layer_name=layer_name,
        reference_key=reference_key,
        reference_categories=reference_categories,
        sample_key=sample_key,
        genome_assembly=genome_assembly,
        window_size=window_size,
        step=step,
        lfc_clip=lfc_clip,
        dynamic_threshold=dynamic_threshold,
        exclude_chromosomes=exclude_chromosomes,
        output_key=output_key,
        minimum_reference_cells=minimum_reference_cells,
        chunksize=chunksize,
        n_jobs=n_jobs,
        max_output_gib=max_output_gib,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
        _worker_owned=_owned,
    )


def analyze_cnv_pca(
    cnv_state: Any,
    *,
    n_comps: int = 30,
    output_key: str = "X_cnv_pca",
    overwrite_existing: bool = False,
    max_output_gib: float = 2.0,
    random_seed: int = 0,
    _owned: bool = False,
) -> tuple[Any, dict[str, Any]]:
    return _cnv_run_pca(
        cnv_state,
        n_comps=n_comps,
        output_key=output_key,
        overwrite_existing=overwrite_existing,
        max_output_gib=max_output_gib,
        random_seed=random_seed,
        openbio_version=PLUGIN_VERSION,
        _worker_owned=_owned,
    )


def analyze_cnv_score(
    cnv_state: Any,
    adata: Any,
    *,
    groupby: str,
    output_key: str = "cnv_score",
    overwrite_existing: bool = False,
    _owned: bool = False,
) -> tuple[Any, Any, dict[str, Any]]:
    return _cnv_run_score(
        cnv_state,
        adata,
        groupby=groupby,
        output_key=output_key,
        overwrite_existing=overwrite_existing,
        openbio_version=PLUGIN_VERSION,
        _worker_owned=_owned,
    )


def _cnv_source(helpers, *, include_state_class):
    sources = [textwrap.dedent(inspect.getsource(function)).strip() for function in helpers]
    if include_state_class:
        class_source = textwrap.dedent(inspect.getsource(_StandaloneCNVState)).strip()
        insertion_index = sources.index(textwrap.dedent(inspect.getsource(_cnv_validate_state)).strip())
        sources.insert(insertion_index, class_source)
        sources.insert(insertion_index + 1, "CNVState = _StandaloneCNVState")
    implementation = "\n\n".join(sources)
    return f"""from __future__ import annotations

import copy
import hashlib
import inspect
import json
import math
import platform
import re
import threading
import warnings
from collections.abc import Mapping, Sequence
from importlib import metadata

CNV_STATE_TYPE = {CNV_STATE_TYPE!r}
CNV_STATE_SCHEMA = {CNV_STATE_SCHEMA!r}
INFER_CNV_NODE_ID = {INFER_CNV_NODE_ID!r}
CNV_PCA_NODE_ID = {CNV_PCA_NODE_ID!r}
CNV_SCORE_NODE_ID = {CNV_SCORE_NODE_ID!r}
_CNV_PROVENANCE_KEY = {_CNV_PROVENANCE_KEY!r}
_CNV_PCA_PROVENANCE_KEY = {_CNV_PCA_PROVENANCE_KEY!r}
_CNV_SCORE_PROVENANCE_KEY = {_CNV_SCORE_PROVENANCE_KEY!r}
_CNV_RUNTIME_LOCK = threading.RLock()

{implementation}
"""


def infer_cnv_code(**parameters: Any) -> str:
    return (
        _cnv_source(_CNV_INFER_SOURCE_HELPERS, include_state_class=True)
        + f"""

def run_infer_cnv(adata):
    \"\"\"Return (validated_cnv_state, strict_summary) without OpenBio imports.\"\"\"
    return _cnv_run_infer(
        adata,
        source_kind={parameters['source_kind']!r},
        layer_name={parameters['layer_name']!r},
        reference_key={parameters['reference_key']!r},
        reference_categories={parameters['reference_categories']!r},
        sample_key={parameters['sample_key']!r},
        genome_assembly={parameters['genome_assembly']!r},
        window_size={parameters['window_size']!r},
        step={parameters['step']!r},
        lfc_clip={parameters['lfc_clip']!r},
        dynamic_threshold={parameters['dynamic_threshold']!r},
        exclude_chromosomes={parameters['exclude_chromosomes']!r},
        output_key={parameters['output_key']!r},
        minimum_reference_cells={parameters['minimum_reference_cells']!r},
        chunksize={parameters['chunksize']!r},
        n_jobs={parameters['n_jobs']!r},
        max_output_gib={parameters['max_output_gib']!r},
        overwrite_existing={parameters['overwrite_existing']!r},
        openbio_version={PLUGIN_VERSION!r},
    )
"""
    )


def cnv_pca_code(**parameters: Any) -> str:
    return (
        _cnv_source(_CNV_PCA_SOURCE_HELPERS, include_state_class=False)
        + f"""

def run_cnv_pca(cnv_state):
    \"\"\"Return (annotated_adata, strict_summary) without OpenBio imports.\"\"\"
    return _cnv_run_pca(
        cnv_state,
        n_comps={parameters['n_comps']!r},
        output_key={parameters['output_key']!r},
        overwrite_existing={parameters['overwrite_existing']!r},
        max_output_gib={parameters['max_output_gib']!r},
        random_seed={parameters['random_seed']!r},
        openbio_version={PLUGIN_VERSION!r},
    )
"""
    )


def cnv_score_code(**parameters: Any) -> str:
    return (
        _cnv_source(_CNV_SCORE_SOURCE_HELPERS, include_state_class=False)
        + f"""

def run_cnv_score(cnv_state, adata):
    \"\"\"Return (annotated_adata, ranked_table, strict_summary) without OpenBio imports.\"\"\"
    return _cnv_run_score(
        cnv_state,
        adata,
        groupby={parameters['groupby']!r},
        output_key={parameters['output_key']!r},
        overwrite_existing={parameters['overwrite_existing']!r},
        openbio_version={PLUGIN_VERSION!r},
    )
"""
    )


__all__ = [
    "CNV_PCA_NODE_ID",
    "CNV_SCORE_NODE_ID",
    "CNV_STATE_SCHEMA",
    "CNV_STATE_TYPE",
    "CNVState",
    "INFER_CNV_NODE_ID",
    "analyze_cnv_pca",
    "analyze_cnv_score",
    "analyze_infer_cnv",
    "cnv_pca_code",
    "cnv_score_code",
    "infer_cnv_code",
]
