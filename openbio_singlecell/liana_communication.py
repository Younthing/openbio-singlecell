from __future__ import annotations

import copy
import csv
import hashlib
import importlib
import inspect
import json
import math
import os
import platform
import stat
import warnings
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from typing import Any

from .liana_result import LIANA_RESULT_COLUMNS, build_liana_result

LIANA_AUDITED_VERSION = "1.9.0"
LIANA_COMMUNICATION_NODE_ID = "OpenBioSingleCellLianaCommunication"
LIANA_RESOURCE_METADATA_FIELDS = (
    "name",
    "version",
    "release_date",
    "download_url",
    "organism",
    "gene_identifier_namespace",
    "scope",
    "license",
    "source_license_review",
    "citation",
)
LIANA_MAX_RESOURCE_BYTES = 512 * 1024**2
LIANA_MAX_RESOURCE_ROWS = 2_000_000
LIANA_REPORT_INTERACTION_LIMIT = 20
LIANA_METHOD_PROFILES = {
    "rank_aggregate": {
        "backend_columns": (
            "sample",
            "source",
            "target",
            "ligand_complex",
            "receptor_complex",
            "lr_means",
            "cellphone_pvals",
            "expr_prod",
            "scaled_weight",
            "lr_logfc",
            "spec_weight",
            "lrscore",
            "specificity_rank",
            "magnitude_rank",
        ),
        "magnitude_field": "magnitude_rank",
        "magnitude_direction": "lower ranks are prioritized",
        "specificity_field": "specificity_rank",
        "specificity_direction": "lower ranks are prioritized; this is not a p-value",
        "sort": (("magnitude_rank", True), ("specificity_rank", True)),
    },
    "cellphonedb": {
        "backend_columns": (
            "sample",
            "ligand",
            "ligand_complex",
            "ligand_means",
            "ligand_props",
            "receptor",
            "receptor_complex",
            "receptor_means",
            "receptor_props",
            "source",
            "target",
            "lr_means",
            "cellphone_pvals",
        ),
        "magnitude_field": "lr_means",
        "magnitude_direction": "higher expression means are prioritized",
        "specificity_field": "cellphone_pvals",
        "specificity_direction": "lower unadjusted cell-label permutation p-values are prioritized",
        "sort": (("lr_means", False), ("cellphone_pvals", True)),
    },
}
LIANA_BY_SAMPLE_SIGNATURE = ("adata", "sample_key", "key_added", "inplace", "verbose", "kwargs")
LIANA_METHOD_SIGNATURES = {
    "rank_aggregate": (
        "adata",
        "groupby",
        "resource_name",
        "expr_prop",
        "min_cells",
        "groupby_pairs",
        "base",
        "aggregate_method",
        "consensus_opts",
        "return_all_lrs",
        "key_added",
        "use_raw",
        "layer",
        "de_method",
        "n_perms",
        "seed",
        "n_jobs",
        "resource",
        "interactions",
        "mdata_kwargs",
        "spatial_key",
        "spatial_kwargs",
        "inplace",
        "verbose",
    ),
    "cellphonedb": (
        "adata",
        "groupby",
        "resource_name",
        "expr_prop",
        "min_cells",
        "groupby_pairs",
        "base",
        "supp_columns",
        "return_all_lrs",
        "key_added",
        "use_raw",
        "layer",
        "de_method",
        "n_perms",
        "seed",
        "n_jobs",
        "resource",
        "interactions",
        "spatial_key",
        "spatial_kwargs",
        "mdata_kwargs",
        "inplace",
        "verbose",
    ),
}


def _liana_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _liana_text(value: Any, *, label: str, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    result = value.strip()
    if not result and not allow_empty:
        raise ValueError(f"{label} cannot be empty.")
    if value != result:
        raise ValueError(f"{label} cannot contain surrounding whitespace.")
    if any(ord(character) < 32 or ord(character) == 127 for character in result):
        raise ValueError(f"{label} cannot contain ASCII control characters.")
    return result


def _liana_integer(value: Any, *, label: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer.")
    if value < minimum or value > maximum:
        raise ValueError(f"{label} must lie in [{minimum}, {maximum}].")
    return value


def _liana_number(
    value: Any,
    *,
    label: str,
    minimum: float,
    maximum: float | None = None,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be a real number.")
    result = float(value)
    if not math.isfinite(result) or result < minimum or (maximum is not None and result > maximum):
        interval = f"[{minimum}, {maximum}]" if maximum is not None else f">= {minimum}"
        raise ValueError(f"{label} must be finite and {interval}.")
    return result


def _liana_axis(values: Sequence[Any], *, label: str) -> tuple[list[str], str]:
    canonical = [_liana_text(value, label=label) for value in list(values)]
    if len(set(canonical)) != len(canonical):
        raise ValueError(f"LIANA requires unique {label} values.")
    digest = hashlib.sha256()
    digest.update(f"openbio-singlecell/liana-{label}-axis/v1\0".encode())
    for value in canonical:
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little", signed=False))
        digest.update(encoded)
    return canonical, digest.hexdigest()


def _liana_matrix_sha256(
    matrix: Any,
    *,
    observations: Sequence[str],
    features: Sequence[str],
    numpy: Any,
    sparse: Any,
) -> str:
    digest = hashlib.sha256()
    digest.update(b"openbio-singlecell/liana-expression/v1\0")
    digest.update(_liana_axis(observations, label="observation identifier")[1].encode("ascii"))
    digest.update(_liana_axis(features, label="feature identifier")[1].encode("ascii"))
    if sparse.issparse(matrix):
        canonical = sparse.csr_matrix(matrix, dtype=float, copy=True)
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
        digest.update(b"csr\0")
        digest.update(numpy.asarray(canonical.indptr, dtype="<i8").tobytes())
        digest.update(numpy.asarray(canonical.indices, dtype="<i8").tobytes())
        digest.update(numpy.asarray(canonical.data, dtype="<f8").tobytes())
    else:
        digest.update(b"dense\0")
        digest.update(numpy.ascontiguousarray(numpy.asarray(matrix, dtype="<f8")).tobytes(order="C"))
    return digest.hexdigest()


def _liana_canonical_labels(series: Any, *, label: str, pandas: Any) -> list[str]:
    if len(series) < 1:
        raise ValueError(f"LIANA {label} column cannot be empty.")
    result: list[str] = []
    for value in series.tolist():
        missing = pandas.isna(value)
        if getattr(missing, "shape", ()) != ():
            raise ValueError(f"LIANA {label} values must be scalar.")
        if bool(missing):
            raise ValueError(f"LIANA {label} column contains missing values.")
        result.append(_liana_text(value, label=f"{label} value"))
    return result


def _liana_history(adata: Any) -> list[Mapping[str, Any]]:
    metadata = adata.uns.get("openbio_singlecell")
    if not isinstance(metadata, Mapping):
        return []
    history = metadata.get("analysis_history")
    if not isinstance(history, Mapping):
        return []
    return [entry for entry in history.values() if isinstance(entry, Mapping)]


def _liana_expression_state(adata: Any, *, source_kind: str, layer_name: str | None) -> tuple[str, str | None]:
    x_state: tuple[str, str | None] = ("unknown", None)
    layer_states: dict[str, tuple[str, str]] = {}
    for entry in _liana_history(adata):
        operation = entry.get("operation")
        parameters = entry.get("parameters")
        parameters = parameters if isinstance(parameters, Mapping) else {}
        if operation == "snapshot_expression":
            layer_states["counts"] = ("counts", "OpenBio Snapshot Expression history")
            if parameters.get("source") == "X":
                x_state = ("counts", "OpenBio Snapshot Expression history")
        elif operation == "normalize_to_layer":
            output_layer = parameters.get("output_layer")
            transform = parameters.get("transform")
            if isinstance(output_layer, str):
                state = "logged" if transform == "log1p" else "normalized" if transform == "none" else "transformed"
                layer_states[output_layer] = (state, "OpenBio Normalize To Layer history")
        elif operation == "pearson_residuals_to_layer":
            output_layer = parameters.get("output_layer")
            if isinstance(output_layer, str):
                layer_states[output_layer] = ("pearson_residuals", "OpenBio Pearson Residuals history")
        elif operation == "scale_to_layer":
            output_layer = parameters.get("output_layer")
            if isinstance(output_layer, str):
                layer_states[output_layer] = ("scaled", "OpenBio Scale history")
        if operation == "normalize_total":
            x_state = ("normalized", "OpenBio Normalize Total history")
        elif operation == "log1p":
            x_state = ("logged", "OpenBio Log1p history")
    if source_kind == "layer":
        return layer_states.get(str(layer_name), ("unknown", None))
    if x_state[0] == "unknown" and isinstance(adata.uns.get("log1p"), Mapping):
        return "logged", "AnnData uns['log1p'] marker"
    return x_state


def _liana_validate_resource_metadata(
    value: str | Mapping[str, Any],
    *,
    organism: str,
) -> dict[str, str]:
    if isinstance(value, str):
        try:
            parsed = json.loads(
                value,
                object_pairs_hook=lambda pairs: _liana_unique_object(pairs),
                parse_constant=lambda token: (_ for _ in ()).throw(
                    ValueError(f"Non-finite JSON constant is not allowed: {token}")
                ),
            )
        except (json.JSONDecodeError, TypeError, ValueError) as exc:
            raise ValueError(f"LIANA resource_metadata_json is invalid strict JSON: {exc}") from exc
    else:
        parsed = copy.deepcopy(value)
    if not isinstance(parsed, Mapping) or set(parsed) != set(LIANA_RESOURCE_METADATA_FIELDS):
        raise ValueError(
            "LIANA resource metadata must contain exactly: " + ", ".join(LIANA_RESOURCE_METADATA_FIELDS) + "."
        )
    normalized = {
        field: _liana_text(parsed[field], label=f"LIANA resource metadata {field}")
        for field in LIANA_RESOURCE_METADATA_FIELDS
    }
    if normalized["organism"] != organism:
        raise ValueError("LIANA resource metadata organism must exactly match the declared organism.")
    if not normalized["download_url"].startswith(("https://", "http://")):
        raise ValueError("LIANA resource metadata download_url must be an explicit HTTP(S) source URL.")
    return normalized


def _liana_unique_object(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key is not allowed: {key!r}")
        result[key] = value
    return result


def _liana_resource_table(
    rows: Sequence[tuple[Any, Any]],
    *,
    features: Sequence[str],
    pandas: Any,
) -> tuple[Any, dict[str, Any]]:
    canonical_rows: list[tuple[str, str]] = []
    for row_index, (ligand, receptor) in enumerate(rows, start=1):
        canonical_rows.append(
            (
                _liana_text(ligand, label=f"LIANA resource ligand row {row_index}"),
                _liana_text(receptor, label=f"LIANA resource receptor row {row_index}"),
            )
        )
    if not canonical_rows:
        raise ValueError("LIANA resource contains no ligand-receptor rows.")
    unique_rows = sorted(set(canonical_rows))
    if len(unique_rows) > LIANA_MAX_RESOURCE_ROWS:
        raise ValueError(
            f"LIANA resource has {len(unique_rows):,} unique rows, exceeding the fixed "
            f"{LIANA_MAX_RESOURCE_ROWS:,}-row guard."
        )
    feature_set = set(features)
    retained_rows = [
        row for row in unique_rows if all(subunit in feature_set for partner in row for subunit in partner.split("_"))
    ]
    if not retained_rows:
        raise ValueError("LIANA resource has no complete ligand-receptor pair on the exact expression feature axis.")

    def rows_sha256(values: Sequence[tuple[str, str]], domain: bytes) -> str:
        digest = hashlib.sha256()
        digest.update(domain)
        for ligand, receptor in values:
            for value in (ligand, receptor):
                encoded = value.encode("utf-8")
                digest.update(len(encoded).to_bytes(8, "little", signed=False))
                digest.update(encoded)
        return digest.hexdigest()

    frame = pandas.DataFrame(retained_rows, columns=["ligand", "receptor"])
    accounting = {
        "parsed_rows": len(canonical_rows),
        "unique_rows": len(unique_rows),
        "duplicate_rows_collapsed": len(canonical_rows) - len(unique_rows),
        "retained_rows_on_exact_feature_axis": len(retained_rows),
        "excluded_rows_missing_subunits": len(unique_rows) - len(retained_rows),
        "canonical_resource_sha256": rows_sha256(unique_rows, b"openbio-singlecell/liana-resource-all/v1\0"),
        "effective_resource_sha256": rows_sha256(retained_rows, b"openbio-singlecell/liana-resource-effective/v1\0"),
    }
    return frame, accounting


def _liana_file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _liana_load_local_resource(
    path: str,
    *,
    features: Sequence[str],
    pandas: Any,
    expected_file_sha256: str | None,
) -> tuple[Any, dict[str, Any]]:
    path = os.path.realpath(_liana_text(path, label="LIANA local resource path"))
    if not path.lower().endswith(".csv"):
        raise ValueError("LIANA local resource must be a CSV file.")
    before = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(before.st_mode):
        raise ValueError("LIANA local resource path must resolve to a regular file.")
    if before.st_size < 1 or before.st_size > LIANA_MAX_RESOURCE_BYTES:
        raise ValueError(f"LIANA local resource size must lie in [1, {LIANA_MAX_RESOURCE_BYTES}] bytes.")
    raw_sha256 = _liana_file_sha256(path)
    if expected_file_sha256 is not None and raw_sha256 != expected_file_sha256:
        raise ValueError(
            f"LIANA local resource SHA-256 changed: expected {expected_file_sha256}, observed {raw_sha256}."
        )
    rows: list[tuple[str, str]] = []
    try:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, strict=True)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValueError("LIANA local resource is empty.") from exc
            if header != ["ligand", "receptor"]:
                raise ValueError("LIANA local resource header must be exactly ligand,receptor.")
            for logical_row, row in enumerate(reader, start=2):
                if len(row) != 2:
                    raise ValueError(
                        f"LIANA local resource logical row {logical_row} has {len(row)} fields; expected 2."
                    )
                rows.append((row[0], row[1]))
                if len(rows) > LIANA_MAX_RESOURCE_ROWS:
                    raise ValueError(f"LIANA local resource exceeds the fixed {LIANA_MAX_RESOURCE_ROWS:,}-row guard.")
    except (UnicodeDecodeError, csv.Error) as exc:
        raise ValueError(f"LIANA local resource is not strict UTF-8 CSV: {exc}") from exc
    after_sha256 = _liana_file_sha256(path)
    after = os.stat(path, follow_symlinks=False)
    if (
        raw_sha256 != after_sha256
        or before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
        or getattr(before, "st_ino", None) != getattr(after, "st_ino", None)
    ):
        raise RuntimeError("LIANA local resource changed while it was being read; execution stopped.")
    frame, accounting = _liana_resource_table(rows, features=features, pandas=pandas)
    accounting.update(
        {
            "raw_size_bytes": int(before.st_size),
            "raw_file_sha256": raw_sha256,
            "resolved_path": path,
        }
    )
    return frame, accounting


def _liana_validate_input(
    adata: Any,
    *,
    sample_key: str,
    condition_key: str,
    identity_key: str,
    annotation_status: str,
    source_kind: str,
    layer_name: str | None,
    min_cells: int,
    numpy: Any,
    pandas: Any,
    sparse: Any,
) -> tuple[
    Any,
    list[str],
    list[str],
    dict[str, str],
    dict[str, Any],
    dict[str, Any],
    list[str],
    dict[str, Any],
]:
    if not hasattr(adata, "obs") or not hasattr(adata, "var_names") or not hasattr(adata, "n_obs"):
        raise TypeError("LIANA Communication requires an in-memory AnnData object.")
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError("LIANA Communication does not support backed AnnData; load it into memory first.")
    sample_key = _liana_text(sample_key, label="LIANA Sample key")
    condition_key = _liana_text(condition_key, label="LIANA Condition key")
    identity_key = _liana_text(identity_key, label="LIANA identity key")
    if len({sample_key, condition_key, identity_key}) != 3:
        raise ValueError("LIANA Sample, Condition, and identity keys must be distinct.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError("LIANA annotation_status must be unknown, provisional, or curated.")
    for key, role in (
        (sample_key, "Sample"),
        (condition_key, "Condition"),
        (identity_key, "identity"),
    ):
        if key not in adata.obs:
            raise ValueError(f"LIANA {role} column not found in obs: {key!r}.")
    if int(adata.n_obs) < 1 or int(adata.n_vars) < 1:
        raise ValueError("LIANA Communication requires nonempty observation and feature axes.")
    observations, observation_axis_sha256 = _liana_axis(adata.obs_names.tolist(), label="observation identifier")
    features, feature_axis_sha256 = _liana_axis(adata.var_names.tolist(), label="feature identifier")

    if source_kind == "X":
        if layer_name is not None:
            raise ValueError("LIANA X expression source cannot carry a layer name.")
        matrix = adata.X
    elif source_kind == "layer":
        layer_name = _liana_text(layer_name, label="LIANA expression layer")
        if layer_name not in adata.layers:
            raise ValueError(f"LIANA expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
    else:
        raise ValueError("LIANA source_kind must be X or layer; Raw snapshots are not accepted as expression.")
    if tuple(getattr(matrix, "shape", ())) != (len(observations), len(features)):
        raise ValueError("LIANA selected expression matrix does not align to its named axes.")
    if not (isinstance(matrix, numpy.ndarray) or sparse.issparse(matrix)):
        raise TypeError("LIANA expression must be an in-memory NumPy array or SciPy sparse matrix.")
    stored = matrix.data if sparse.issparse(matrix) else numpy.asarray(matrix).ravel()
    stored = numpy.asarray(stored)
    if stored.dtype.kind not in "iuf" or not bool(numpy.isfinite(stored).all()):
        raise ValueError("LIANA expression must contain only finite real numeric values.")
    value_min = float(min(0.0, stored.min(initial=0.0)))
    integer_like = bool(numpy.allclose(stored, numpy.rint(stored), rtol=0.0, atol=1e-8))
    state, state_evidence = _liana_expression_state(adata, source_kind=source_kind, layer_name=layer_name)
    warnings_list: list[str] = []
    if state in {"counts", "pearson_residuals", "scaled", "transformed"}:
        warnings_list.append(
            f"The explicitly selected LIANA expression is described as {state!r} ({state_evidence}); the backend "
            "will use those finite values unchanged, but interaction scores do not have the recommended normalized "
            "log1p interpretation."
        )
    elif state == "unknown" and integer_like and value_min >= 0.0:
        warnings_list.append(
            "LIANA expression provenance is unknown and values are nonnegative integer-like/count-like. The "
            "explicit expert selection was honored, but library size and ties can dominate the result."
        )
    elif state == "unknown":
        warnings_list.append(
            "Expression provenance is unknown but values are non-count-like; LIANA assumes caller-supplied "
            "library-size-normalized log1p expression."
        )
    elif state == "normalized":
        warnings_list.append(
            "Expression is library-size normalized but no log1p transform is recorded; score magnitudes depend on "
            "this explicitly disclosed representation."
        )
    elif state != "logged":
        warnings_list.append(
            f"LIANA expression state {state!r} is nonstandard; the explicitly selected finite values were used unchanged."
        )
    if value_min < 0.0:
        warnings_list.append(
            "The selected LIANA expression contains negative values. The calculation remains available, but "
            "magnitude and expression-proportion semantics differ from nonnegative normalized/log1p practice."
        )
    warnings_list.append(
        "Full-gene completeness cannot be proven from the current AnnData feature axis; the selected feature set "
        "was used as declared by the caller."
    )

    sample_values = _liana_canonical_labels(adata.obs[sample_key], label="Sample", pandas=pandas)
    condition_values = _liana_canonical_labels(adata.obs[condition_key], label="Condition", pandas=pandas)
    identity_values = _liana_canonical_labels(adata.obs[identity_key], label="identity", pandas=pandas)
    sample_order = list(dict.fromkeys(sample_values))
    if len(sample_order) < 2:
        raise ValueError("LIANA Communication requires at least two biological Samples.")
    condition_map: dict[str, str] = {}
    for sample, condition in zip(sample_values, condition_values, strict=True):
        previous = condition_map.setdefault(sample, condition)
        if previous != condition:
            raise ValueError(f"LIANA Sample {sample!r} maps to more than one Condition.")
    counts: dict[tuple[str, str], int] = {}
    identity_order: dict[str, list[str]] = {sample: [] for sample in sample_order}
    for sample, identity in zip(sample_values, identity_values, strict=True):
        key = (sample, identity)
        counts[key] = counts.get(key, 0) + 1
        if identity not in identity_order[sample]:
            identity_order[sample].append(identity)
    eligible_by_sample: dict[str, list[str]] = {}
    strata = []
    for sample in sample_order:
        eligible_by_sample[sample] = []
        for identity in identity_order[sample]:
            count = counts[(sample, identity)]
            eligible = count >= min_cells
            if eligible:
                eligible_by_sample[sample].append(identity)
            strata.append(
                {
                    "sample": sample,
                    "condition": condition_map[sample],
                    "identity": identity,
                    "cells": count,
                    "eligible": eligible,
                }
            )
        if len(eligible_by_sample[sample]) < 2:
            raise ValueError(f"LIANA Sample {sample!r} has fewer than two identities with at least {min_cells} cells.")
    expression_sha256 = _liana_matrix_sha256(
        matrix,
        observations=observations,
        features=features,
        numpy=numpy,
        sparse=sparse,
    )
    expression = {
        "source_kind": source_kind,
        "layer_name": layer_name,
        "state": state,
        "state_evidence": state_evidence,
        "full_gene_completeness_verified": False,
        "full_gene_completeness_basis": "caller-selected current feature axis",
        "observation_axis_sha256": observation_axis_sha256,
        "feature_axis_sha256": feature_axis_sha256,
        "expression_content_sha256": expression_sha256,
        "cells": len(observations),
        "genes": len(features),
    }
    design = {
        "sample_order": sample_order,
        "condition_order": list(dict.fromkeys(condition_values)),
        "sample_to_condition": condition_map,
        "sample_count": len(sample_order),
        "condition_count": len(set(condition_values)),
        "identity_count": len(set(identity_values)),
        "sample_identity_strata": strata,
        "eligible_identity_count_by_sample": {sample: len(values) for sample, values in eligible_by_sample.items()},
        "excluded_strata": sum(not record["eligible"] for record in strata),
    }
    caller_identity = {
        "observations": observations,
        "features": features,
        "sample_values": sample_values,
        "condition_values": condition_values,
        "identity_values": identity_values,
        "expression_sha256": expression_sha256,
    }
    return (
        matrix,
        observations,
        features,
        condition_map,
        expression,
        {**design, "eligible_identities": eligible_by_sample},
        warnings_list,
        caller_identity,
    )


def _liana_current_caller_identity(
    adata: Any,
    *,
    sample_key: str,
    condition_key: str,
    identity_key: str,
    source_kind: str,
    layer_name: str | None,
    numpy: Any,
    pandas: Any,
    sparse: Any,
) -> dict[str, Any]:
    observations, _ = _liana_axis(adata.obs_names.tolist(), label="observation identifier")
    features, _ = _liana_axis(adata.var_names.tolist(), label="feature identifier")
    matrix = adata.X if source_kind == "X" else adata.layers[layer_name]
    return {
        "observations": observations,
        "features": features,
        "sample_values": _liana_canonical_labels(adata.obs[sample_key], label="Sample", pandas=pandas),
        "condition_values": _liana_canonical_labels(adata.obs[condition_key], label="Condition", pandas=pandas),
        "identity_values": _liana_canonical_labels(adata.obs[identity_key], label="identity", pandas=pandas),
        "expression_sha256": _liana_matrix_sha256(
            matrix,
            observations=observations,
            features=features,
            numpy=numpy,
            sparse=sparse,
        ),
    }


def _liana_signature(callable_object: Any, *, label: str, expected: Sequence[str]) -> None:
    if not callable(callable_object):
        raise RuntimeError(f"LIANA 1.9.0 public {label} is not callable.")
    try:
        signature = inspect.signature(callable_object)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"Cannot inspect LIANA 1.9.0 public {label} signature.") from exc
    observed = tuple(signature.parameters)
    if observed != tuple(expected):
        raise RuntimeError(
            f"LIANA 1.9.0 public {label} signature drifted: expected {tuple(expected)!r}, observed {observed!r}."
        )
    if expected and expected[-1] == "kwargs":
        parameter = signature.parameters["kwargs"]
        if parameter.kind is not inspect.Parameter.VAR_KEYWORD:
            raise RuntimeError(f"LIANA 1.9.0 public {label} kwargs parameter is not **kwargs.")


def _liana_runtime(
    *,
    method: str,
    resource_mode: str,
    liana_module: Any | None,
) -> tuple[Any, Any]:
    injected = liana_module is not None
    if liana_module is None:
        try:
            liana_module = importlib.import_module("liana")
        except (ImportError, OSError) as exc:
            raise RuntimeError(
                "LIANA Communication requires optional liana==1.9.0 with pandas<3; install that audited "
                f"environment before execution ({exc})."
            ) from exc
    version = getattr(liana_module, "__version__", None)
    if version != LIANA_AUDITED_VERSION:
        raise RuntimeError(f"LIANA Communication requires exact liana=={LIANA_AUDITED_VERSION}; observed {version!r}.")
    pandas = importlib.import_module("pandas")
    if not injected:
        pandas_version = getattr(pandas, "__version__", "unknown")
        try:
            pandas_major = int(str(pandas_version).split(".", 1)[0])
        except ValueError as exc:
            raise RuntimeError(f"Cannot interpret pandas version {pandas_version!r} for LIANA 1.9.0.") from exc
        if pandas_major >= 3:
            raise RuntimeError(
                f"LIANA 1.9.0 declares pandas<3 and this audited adapter rejects pandas {pandas_version}."
            )
    mt = getattr(liana_module, "mt", None)
    rs = getattr(liana_module, "rs", None)
    method_object = getattr(mt, method, None)
    if method_object is None:
        raise RuntimeError(f"LIANA 1.9.0 public mt.{method} is unavailable.")
    _liana_signature(
        getattr(method_object, "by_sample", None),
        label=f"mt.{method}.by_sample",
        expected=LIANA_BY_SAMPLE_SIGNATURE,
    )
    _liana_signature(
        method_object,
        label=f"mt.{method}",
        expected=LIANA_METHOD_SIGNATURES[method],
    )
    if resource_mode == "bundled_human":
        _liana_signature(
            getattr(rs, "select_resource", None),
            label="rs.select_resource",
            expected=("resource_name",),
        )
    return liana_module, method_object


def _liana_resolve_resource(
    liana_module: Any,
    *,
    resource_mode: str,
    resource_name: str,
    resource_path: str | None,
    resource_metadata_json: str | Mapping[str, Any],
    organism: str,
    features: Sequence[str],
    pandas: Any,
    expected_file_sha256: str | None,
    expected_canonical_resource_sha256: str | None,
    expected_effective_resource_sha256: str | None,
) -> tuple[Any, dict[str, Any]]:
    organism = _liana_text(organism, label="LIANA organism")
    resource_name = _liana_text(resource_name, label="LIANA resource name")
    metadata = _liana_validate_resource_metadata(resource_metadata_json, organism=organism)
    if metadata["name"] != resource_name:
        raise ValueError("LIANA resource metadata name must exactly match resource_name.")
    if resource_mode == "bundled_human":
        if organism != "Homo sapiens":
            raise ValueError("LIANA bundled_human resources require organism 'Homo sapiens'.")
        if resource_name.lower() == "mouseconsensus":
            raise ValueError("LIANA bundled_human mode rejects the mouseconsensus resource.")
        if resource_path is not None:
            raise ValueError("LIANA bundled_human resource does not accept a local resource path.")
        loaded = liana_module.rs.select_resource(resource_name=resource_name)
        if not isinstance(loaded, pandas.DataFrame) or list(loaded.columns) != ["ligand", "receptor"]:
            raise RuntimeError(
                "LIANA 1.9.0 rs.select_resource returned an unexpected schema; expected ligand,receptor."
            )
        resource, accounting = _liana_resource_table(
            list(loaded[["ligand", "receptor"]].itertuples(index=False, name=None)),
            features=features,
            pandas=pandas,
        )
        accounting.update(
            {
                "raw_size_bytes": None,
                "raw_file_sha256": None,
                "resolved_path": f"liana=={LIANA_AUDITED_VERSION}:rs.select_resource({resource_name!r})",
            }
        )
    elif resource_mode == "local_resource":
        if resource_path is None:
            raise ValueError("LIANA local_resource mode requires a resource CSV path.")
        resource, accounting = _liana_load_local_resource(
            resource_path,
            features=features,
            pandas=pandas,
            expected_file_sha256=expected_file_sha256,
        )
    else:
        raise ValueError("LIANA resource_mode must be bundled_human or local_resource.")
    if (
        expected_canonical_resource_sha256 is not None
        and accounting["canonical_resource_sha256"] != expected_canonical_resource_sha256
    ):
        raise ValueError(
            "LIANA canonical resource fingerprint changed: expected "
            f"{expected_canonical_resource_sha256}, observed {accounting['canonical_resource_sha256']}."
        )
    if (
        expected_effective_resource_sha256 is not None
        and accounting["effective_resource_sha256"] != expected_effective_resource_sha256
    ):
        raise ValueError(
            "LIANA effective resource fingerprint changed: expected "
            f"{expected_effective_resource_sha256}, observed {accounting['effective_resource_sha256']}."
        )
    provenance = {
        "mode": resource_mode,
        "name": resource_name,
        "metadata": metadata,
        "metadata_sha256": _liana_json_sha256(metadata),
        "accounting": accounting,
    }
    return resource, provenance


def _liana_resource_frame_sha256(frame: Any) -> str:
    digest = hashlib.sha256()
    digest.update(b"openbio-singlecell/liana-effective-resource-frame/v1\0")
    for ligand, receptor in frame[["ligand", "receptor"]].itertuples(index=False, name=None):
        for value in (ligand, receptor):
            encoded = value.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "little", signed=False))
            digest.update(encoded)
    return digest.hexdigest()


def _liana_validate_backend_output(
    table: Any,
    *,
    method: str,
    sample_order: Sequence[str],
    condition_map: Mapping[str, str],
    eligible_identities: Mapping[str, Sequence[str]],
    effective_resource: Any,
    max_output_rows: int,
    numpy: Any,
    pandas: Any,
) -> Any:
    if not isinstance(table, pandas.DataFrame):
        raise RuntimeError("LIANA 1.9.0 by_sample must return a pandas DataFrame with inplace=False.")
    expected_backend_columns = list(LIANA_METHOD_PROFILES[method]["backend_columns"])
    if list(table.columns) != expected_backend_columns:
        raise RuntimeError(
            f"LIANA 1.9.0 {method} returned an unexpected complete column family: expected "
            f"{expected_backend_columns!r}, observed {list(table.columns)!r}."
        )
    if len(table) < 1:
        raise ValueError("LIANA returned no candidate interactions after the declared filters.")
    if len(table) > max_output_rows:
        raise ValueError(f"LIANA returned {len(table):,} rows, exceeding max_output_rows={max_output_rows:,}.")
    work = table.copy(deep=True).reset_index(drop=True)
    text_columns = ["sample", "source", "target", "ligand_complex", "receptor_complex"]
    if method == "cellphonedb":
        text_columns.extend(["ligand", "receptor"])
    for column in text_columns:
        work[column] = [_liana_text(value, label=f"LIANA backend {column}") for value in work[column].tolist()]
    numeric_columns = [column for column in work.columns if column not in text_columns]
    for column in numeric_columns:
        if not pandas.api.types.is_numeric_dtype(work[column].dtype) or pandas.api.types.is_bool_dtype(
            work[column].dtype
        ):
            raise RuntimeError(f"LIANA backend column {column!r} is not numeric.")
        values = work[column].to_numpy(dtype=float, copy=True)
        if not bool(numpy.isfinite(values).all()):
            raise RuntimeError(f"LIANA backend column {column!r} contains non-finite values.")
        work[column] = values
    sample_set = set(work["sample"])
    if sample_set != set(sample_order):
        raise RuntimeError(
            "LIANA backend Sample family does not exactly match the validated input Samples: "
            f"expected {list(sample_order)!r}, observed {sorted(sample_set)!r}."
        )
    for sample, source, target in work[["sample", "source", "target"]].itertuples(index=False, name=None):
        eligible = set(eligible_identities[sample])
        if source not in eligible or target not in eligible:
            raise RuntimeError(
                "LIANA backend returned an interaction involving an ineligible Sample-by-identity stratum."
            )
    resource_pairs = set(effective_resource[["ligand", "receptor"]].itertuples(index=False, name=None))
    for pair in work[["ligand_complex", "receptor_complex"]].itertuples(index=False, name=None):
        if pair not in resource_pairs:
            raise RuntimeError("LIANA backend returned an interaction outside the pinned effective resource.")
    grain = ["sample", "source", "target", "ligand_complex", "receptor_complex"]
    if bool(work.duplicated(grain, keep=False).any()):
        raise RuntimeError("LIANA backend returned duplicate rows at the declared interaction grain.")
    work.insert(1, "condition", work["sample"].map(condition_map))
    if bool(work["condition"].isna().any()):
        raise RuntimeError("LIANA backend result could not be joined to the validated Condition mapping.")
    work["condition"] = [_liana_text(value, label="LIANA joined Condition") for value in work["condition"].tolist()]
    output_columns = list(LIANA_RESULT_COLUMNS[method])
    work = work[output_columns]
    sample_positions = {sample: position for position, sample in enumerate(sample_order)}
    work["__sample_order"] = work["sample"].map(sample_positions)
    sort_columns = ["__sample_order"]
    ascending = [True]
    for column, direction in LIANA_METHOD_PROFILES[method]["sort"]:
        sort_columns.append(column)
        ascending.append(direction)
    sort_columns.extend(["source", "target", "ligand_complex", "receptor_complex"])
    ascending.extend([True, True, True, True])
    work = work.sort_values(sort_columns, ascending=ascending, kind="mergesort", ignore_index=True)
    work = work.drop(columns="__sample_order")
    return work


def _liana_package_version(package: str) -> str:
    try:
        return importlib_metadata.version(package)
    except importlib_metadata.PackageNotFoundError:
        return "not-installed"


def _liana_references(method: str, resource: Mapping[str, Any]) -> list[dict[str, Any]]:
    references = [
        {
            "citation": (
                "Dimitrov D, et al. Comparison of methods and resources for cell-cell communication inference "
                "from single-cell RNA-Seq data. Nature Communications. 2022;13:3224."
            ),
            "doi": "10.1038/s41467-022-30755-0",
            "url": "https://doi.org/10.1038/s41467-022-30755-0",
            "kind": "method",
        },
        {
            "citation": (
                "Dimitrov D, et al. LIANA+ provides an all-in-one framework for cell-cell communication "
                "inference. Nature Cell Biology. 2024;26:1613-1622."
            ),
            "doi": "10.1038/s41556-024-01469-w",
            "url": "https://doi.org/10.1038/s41556-024-01469-w",
            "kind": "software",
        },
    ]
    if method == "rank_aggregate":
        references.append(
            {
                "citation": (
                    "Kolde R, et al. Robust rank aggregation for gene list integration and meta-analysis. "
                    "Bioinformatics. 2012;28:573-580."
                ),
                "doi": "10.1093/bioinformatics/btr709",
                "url": "https://doi.org/10.1093/bioinformatics/btr709",
                "kind": "method",
            }
        )
    else:
        references.append(
            {
                "citation": (
                    "Efremova M, et al. CellPhoneDB: inferring cell-cell communication from combined expression "
                    "of multi-subunit ligand-receptor complexes. Nature Protocols. 2020;15:1484-1506."
                ),
                "doi": "10.1038/s41596-020-0292-x",
                "url": "https://doi.org/10.1038/s41596-020-0292-x",
                "kind": "method",
            }
        )
    metadata = resource["metadata"]
    references.append(
        {
            "citation": metadata["citation"],
            "doi": None,
            "url": metadata["download_url"],
            "kind": "resource",
        }
    )
    return references


def _liana_run_impl(
    adata: Any,
    *,
    sample_key: str,
    condition_key: str,
    identity_key: str,
    annotation_status: str,
    organism: str,
    method: str,
    resource_mode: str,
    resource_name: str,
    resource_path: str | None,
    resource_metadata_json: str | Mapping[str, Any],
    source_kind: str,
    layer_name: str | None,
    expression_proportion: float,
    min_cells_per_identity_sample: int,
    permutations: int,
    random_seed: int,
    jobs: int,
    max_output_rows: int,
    max_working_memory_gib: float,
    openbio_version: str,
    expected_file_sha256: str | None = None,
    expected_canonical_resource_sha256: str | None = None,
    expected_effective_resource_sha256: str | None = None,
    liana_module: Any | None = None,
) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    sample_key = _liana_text(sample_key, label="LIANA Sample key")
    condition_key = _liana_text(condition_key, label="LIANA Condition key")
    identity_key = _liana_text(identity_key, label="LIANA identity key")
    organism = _liana_text(organism, label="LIANA organism")
    resource_name = _liana_text(resource_name, label="LIANA resource name")
    if method not in LIANA_METHOD_PROFILES:
        raise ValueError(f"Unsupported LIANA communication method: {method!r}.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError("LIANA annotation_status must be unknown, provisional, or curated.")
    expression_proportion = _liana_number(
        expression_proportion,
        label="LIANA expression_proportion",
        minimum=0.0,
        maximum=1.0,
    )
    min_cells_per_identity_sample = _liana_integer(
        min_cells_per_identity_sample,
        label="LIANA min_cells_per_identity_sample",
        minimum=2,
        maximum=2**31 - 1,
    )
    permutations = _liana_integer(
        permutations,
        label="LIANA permutations",
        minimum=1000,
        maximum=10_000_000,
    )
    random_seed = _liana_integer(
        random_seed,
        label="LIANA random_seed",
        minimum=1,
        maximum=2**32 - 1,
    )
    jobs = _liana_integer(jobs, label="LIANA jobs", minimum=1, maximum=1)
    max_output_rows = _liana_integer(
        max_output_rows,
        label="LIANA max_output_rows",
        minimum=1,
        maximum=2**31 - 1,
    )
    max_working_memory_gib = _liana_number(
        max_working_memory_gib,
        label="LIANA max_working_memory_gib",
        minimum=0.001,
        maximum=1024.0,
    )
    openbio_version = _liana_text(openbio_version, label="OpenBio version")
    numpy = importlib.import_module("numpy")
    pandas = importlib.import_module("pandas")
    sparse = importlib.import_module("scipy.sparse")
    liana_module, method_object = _liana_runtime(
        method=method,
        resource_mode=resource_mode,
        liana_module=liana_module,
    )
    (
        matrix,
        observations,
        features,
        condition_map,
        expression,
        design,
        warnings_list,
        caller_identity,
    ) = _liana_validate_input(
        adata,
        sample_key=sample_key,
        condition_key=condition_key,
        identity_key=identity_key,
        annotation_status=annotation_status,
        source_kind=source_kind,
        layer_name=layer_name,
        min_cells=min_cells_per_identity_sample,
        numpy=numpy,
        pandas=pandas,
        sparse=sparse,
    )
    resource, resource_provenance = _liana_resolve_resource(
        liana_module,
        resource_mode=resource_mode,
        resource_name=resource_name,
        resource_path=resource_path,
        resource_metadata_json=resource_metadata_json,
        organism=organism,
        features=features,
        pandas=pandas,
        expected_file_sha256=expected_file_sha256,
        expected_canonical_resource_sha256=expected_canonical_resource_sha256,
        expected_effective_resource_sha256=expected_effective_resource_sha256,
    )
    resource_frame_sha256 = _liana_resource_frame_sha256(resource)
    eligible_identities = design["eligible_identities"]
    tested_rows_upper_bound = sum(
        len(eligible_identities[sample]) ** 2 * len(resource) for sample in design["sample_order"]
    )
    if tested_rows_upper_bound > max_output_rows:
        raise ValueError(
            f"LIANA complete output family can contain up to {tested_rows_upper_bound:,} rows, exceeding "
            f"max_output_rows={max_output_rows:,}; narrow the explicit resource before analysis."
        )
    matrix_storage_bytes = int(getattr(matrix, "nbytes", 0))
    if sparse.issparse(matrix):
        matrix_storage_bytes = int(matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes)
    dense_expression_bytes = len(observations) * len(features) * 8
    memory_components = {
        "selected_expression_storage_bytes": matrix_storage_bytes,
        "backend_dense_expression_workspace_bytes": dense_expression_bytes * 2,
        "backend_and_canonical_result_bytes": tested_rows_upper_bound * len(LIANA_RESULT_COLUMNS[method]) * 16,
        "resource_copy_bytes": len(resource) * 2 * 64,
    }
    estimated_working_bytes = sum(memory_components.values())
    if estimated_working_bytes > int(max_working_memory_gib * 1024**3):
        raise MemoryError(
            f"LIANA estimates {estimated_working_bytes / 1024**3:.3f} GiB working memory, exceeding "
            f"max_working_memory_gib={max_working_memory_gib:.3f}."
        )
    work = adata
    work_matrix = work.X if source_kind == "X" else work.layers[layer_name]
    work_expression_before = _liana_matrix_sha256(
        work_matrix,
        observations=observations,
        features=features,
        numpy=numpy,
        sparse=sparse,
    )
    backend_resource = resource.copy(deep=True)
    backend_resource_sha256 = _liana_resource_frame_sha256(backend_resource)
    backend_kwargs: dict[str, Any] = {
        "groupby": identity_key,
        "resource_name": resource_name,
        "expr_prop": expression_proportion,
        "min_cells": min_cells_per_identity_sample,
        "groupby_pairs": None,
        "base": math.e,
        "return_all_lrs": False,
        "use_raw": False,
        "layer": layer_name if source_kind == "layer" else None,
        "de_method": "t-test",
        "n_perms": permutations,
        "seed": random_seed,
        "n_jobs": jobs,
        "resource": backend_resource,
        "interactions": None,
        "spatial_key": None,
        "spatial_kwargs": None,
        "mdata_kwargs": None,
    }
    if method == "rank_aggregate":
        backend_kwargs.update({"aggregate_method": "rra", "consensus_opts": None})
    else:
        backend_kwargs["supp_columns"] = None
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        backend_table = method_object.by_sample(
            work,
            sample_key=sample_key,
            key_added="__openbio_liana_private_scratch__",
            inplace=False,
            verbose=False,
            **backend_kwargs,
        )
    captured_backend_warnings = list(
        dict.fromkeys(f"LIANA backend {warning.category.__name__}: {warning.message}" for warning in caught_warnings)
    )
    warnings_list.extend(captured_backend_warnings)
    if _liana_resource_frame_sha256(backend_resource) != backend_resource_sha256 or not backend_resource.equals(
        resource
    ):
        raise RuntimeError("LIANA backend modified its private pinned resource input.")
    if list(work.obs_names) != observations or list(work.var_names) != features:
        raise RuntimeError("LIANA changed the worker-owned AnnData axes.")
    work_matrix_after = work.X if source_kind == "X" else work.layers[layer_name]
    if (
        _liana_matrix_sha256(
            work_matrix_after,
            observations=observations,
            features=features,
            numpy=numpy,
            sparse=sparse,
        )
        != work_expression_before
    ):
        raise RuntimeError("LIANA changed the worker-owned selected expression matrix.")
    if _liana_resource_frame_sha256(resource) != resource_frame_sha256:
        raise RuntimeError("LIANA backend modified the caller-owned canonical resource table.")
    if (
        _liana_current_caller_identity(
            adata,
            sample_key=sample_key,
            condition_key=condition_key,
            identity_key=identity_key,
            source_kind=source_kind,
            layer_name=layer_name,
            numpy=numpy,
            pandas=pandas,
            sparse=sparse,
        )
        != caller_identity
    ):
        raise RuntimeError("LIANA backend modified caller-owned AnnData scientific inputs.")
    table = _liana_validate_backend_output(
        backend_table,
        method=method,
        sample_order=design["sample_order"],
        condition_map=condition_map,
        eligible_identities=eligible_identities,
        effective_resource=resource,
        max_output_rows=max_output_rows,
        numpy=numpy,
        pandas=pandas,
    )
    if annotation_status != "curated":
        warnings_list.append(
            f"Identity labels are declared {annotation_status}; candidate interactions require annotation review."
        )
    accounting = resource_provenance["accounting"]
    if accounting["duplicate_rows_collapsed"]:
        warnings_list.append(
            f"Collapsed {accounting['duplicate_rows_collapsed']:,} exact duplicate resource rows before analysis."
        )
    if accounting["excluded_rows_missing_subunits"]:
        warnings_list.append(
            f"Excluded {accounting['excluded_rows_missing_subunits']:,} resource rows whose full complex subunits "
            "were absent from the exact expression feature axis."
        )
    warnings_list = list(dict.fromkeys(warnings_list))
    design_public = {key: copy.deepcopy(value) for key, value in design.items() if key != "eligible_identities"}
    profile = LIANA_METHOD_PROFILES[method]
    parameters = {
        "sample_key": sample_key,
        "condition_key": condition_key,
        "identity_key": identity_key,
        "annotation_status": annotation_status,
        "organism": organism,
        "method": method,
        "resource_mode": resource_mode,
        "resource_name": resource_name,
        "source_kind": source_kind,
        "layer_name": layer_name,
        "expression_proportion": expression_proportion,
        "min_cells_per_identity_sample": min_cells_per_identity_sample,
        "permutations": permutations,
        "random_seed": random_seed,
        "jobs": jobs,
        "return_all_lrs": False,
        "use_raw": False,
        "de_method": "t-test",
        "spatial_key": None,
        "max_output_rows": max_output_rows,
        "max_working_memory_gib": max_working_memory_gib,
    }
    result_record = {
        "rows": int(len(table)),
        "columns": list(table.columns),
        "complete_backend_family_returned": True,
        "tested_rows_upper_bound": tested_rows_upper_bound,
        "magnitude_field": profile["magnitude_field"],
        "magnitude_direction": profile["magnitude_direction"],
        "specificity_field": profile["specificity_field"],
        "specificity_direction": profile["specificity_direction"],
        "report_interaction_limit": LIANA_REPORT_INTERACTION_LIMIT,
        "report_interactions": table.head(LIANA_REPORT_INTERACTION_LIMIT).to_dict(orient="records"),
    }
    versions = {
        "python": platform.python_version(),
        "openbio-singlecell": openbio_version,
        "liana": str(liana_module.__version__),
        "numpy": str(numpy.__version__),
        "pandas": str(pandas.__version__),
        "scipy": _liana_package_version("scipy"),
        "anndata": _liana_package_version("anndata"),
        "scanpy": _liana_package_version("scanpy"),
    }
    limitations = [
        "Results are descriptive within-Sample expression-compatible ligand-receptor candidates; no "
        "replicate-aware Condition contrast is performed.",
        "Cell-label permutation p-values and aggregate ranks are not Condition-level inference and are not "
        "reported as multiple-testing-adjusted significance.",
        "Dissociated scRNA-seq transcript abundance does not establish protein abundance, secretion, binding, "
        "spatial contact, receptor activation, or signaling causality.",
        "Exact resource-to-feature string matching does not independently prove organism or identifier namespace; "
        "the declared resource metadata and upstream feature annotation must be reviewed together.",
        "Technical batch and paired-design effects are not modeled by this node.",
    ]
    summary = {
        "schema_version": 1,
        "node_id": LIANA_COMMUNICATION_NODE_ID,
        "status": "ok",
        "methods": (
            f"LIANA {LIANA_AUDITED_VERSION} public mt.{method}.by_sample ran separately within each biological "
            f"Sample using explicit {identity_key!r} identities, pinned normalized expression, and a "
            "license-reviewed SHA-256-bound ligand-receptor resource."
        ),
        "results": (
            f"Returned the complete validated family of {len(table):,} within-Sample candidate interactions "
            f"across {design_public['sample_count']:,} Samples; Condition labels are provenance only and were "
            "not tested."
        ),
        "key_results": {
            "design": design_public,
            "expression": expression,
            "resource": resource_provenance,
            "result": result_record,
            "memory_guard": {
                "components": memory_components,
                "estimated_working_bytes": estimated_working_bytes,
            },
        },
        "parameters": parameters,
        "references": _liana_references(method, resource_provenance),
        "software_versions": versions,
        "warnings": warnings_list,
        "limitations": limitations,
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    provenance = {
        "schema_version": 1,
        "method": method,
        "method_profile": {
            "magnitude_field": profile["magnitude_field"],
            "magnitude_direction": profile["magnitude_direction"],
            "specificity_field": profile["specificity_field"],
            "specificity_direction": profile["specificity_direction"],
        },
        "roles": {
            "sample_key": sample_key,
            "condition_key": condition_key,
            "identity_key": identity_key,
            "annotation_status": annotation_status,
        },
        "organism": organism,
        "expression": expression,
        "design": design_public,
        "resource": resource_provenance,
        "backend": {
            "package": "liana",
            "version": str(liana_module.__version__),
            "public_method": f"mt.{method}.by_sample",
            "captured_warnings": captured_backend_warnings,
        },
        "parameters": parameters,
        "result": result_record,
        "summary_sha256": _liana_json_sha256(summary),
    }
    json.dumps(provenance, ensure_ascii=False, allow_nan=False)
    return table, provenance, summary


def run_liana_communication(
    adata: Any,
    *,
    copy_table: bool = True,
    **parameters: Any,
) -> tuple[Any, dict[str, Any]]:
    table, provenance, summary = _liana_run_impl(adata, **parameters)
    numpy = importlib.import_module("numpy")
    pandas = importlib.import_module("pandas")
    artifact = build_liana_result(
        table=table,
        method=summary["parameters"]["method"],
        provenance=provenance,
        numpy=numpy,
        pandas=pandas,
        copy_table=copy_table,
    )
    return artifact, summary


def liana_resource_cache_fingerprint(path: str, resource_metadata_json: str) -> tuple[Any, ...]:
    canonical = os.path.realpath(path)
    stat_result = os.stat(canonical, follow_symlinks=False)
    if not stat.S_ISREG(stat_result.st_mode):
        raise ValueError("LIANA local resource cache input must resolve to a regular file.")
    if stat_result.st_size < 1 or stat_result.st_size > LIANA_MAX_RESOURCE_BYTES:
        raise ValueError(f"LIANA local resource cache input size must lie in [1, {LIANA_MAX_RESOURCE_BYTES}] bytes.")
    return (
        "openbio-liana-resource-v1",
        canonical,
        int(stat_result.st_size),
        int(stat_result.st_mtime_ns),
        _liana_file_sha256(canonical),
        _liana_json_sha256(
            _liana_validate_resource_metadata(
                resource_metadata_json,
                organism=_liana_validate_resource_metadata_unbound_organism(resource_metadata_json),
            )
        ),
    )


def _liana_validate_resource_metadata_unbound_organism(value: str) -> str:
    try:
        parsed = json.loads(value, object_pairs_hook=lambda pairs: _liana_unique_object(pairs))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        raise ValueError(f"LIANA resource_metadata_json is invalid strict JSON: {exc}") from exc
    if not isinstance(parsed, Mapping) or "organism" not in parsed:
        raise ValueError("LIANA resource metadata must declare organism.")
    return _liana_text(parsed["organism"], label="LIANA resource metadata organism")


def liana_communication_code(*, parameters: Mapping[str, Any]) -> str:
    if not isinstance(parameters, Mapping):
        raise TypeError("LIANA generated-code parameters must be a mapping.")
    payload = copy.deepcopy(dict(parameters))
    json.dumps(payload, ensure_ascii=False, allow_nan=False)
    helpers = (
        _liana_json_sha256,
        _liana_text,
        _liana_integer,
        _liana_number,
        _liana_axis,
        _liana_matrix_sha256,
        _liana_canonical_labels,
        _liana_history,
        _liana_expression_state,
        _liana_unique_object,
        _liana_validate_resource_metadata,
        _liana_resource_table,
        _liana_file_sha256,
        _liana_load_local_resource,
        _liana_validate_input,
        _liana_current_caller_identity,
        _liana_signature,
        _liana_runtime,
        _liana_resolve_resource,
        _liana_resource_frame_sha256,
        _liana_validate_backend_output,
        _liana_package_version,
        _liana_references,
        _liana_run_impl,
    )
    helper_source = "\n\n".join(inspect.getsource(helper) for helper in helpers)
    code = f'''from __future__ import annotations

import copy
import csv
import hashlib
import importlib
import inspect
import json
import math
import os
import platform
import stat
import warnings
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from typing import Any

LIANA_AUDITED_VERSION = {LIANA_AUDITED_VERSION!r}
LIANA_COMMUNICATION_NODE_ID = {LIANA_COMMUNICATION_NODE_ID!r}
LIANA_RESOURCE_METADATA_FIELDS = {LIANA_RESOURCE_METADATA_FIELDS!r}
LIANA_MAX_RESOURCE_BYTES = {LIANA_MAX_RESOURCE_BYTES!r}
LIANA_MAX_RESOURCE_ROWS = {LIANA_MAX_RESOURCE_ROWS!r}
LIANA_REPORT_INTERACTION_LIMIT = {LIANA_REPORT_INTERACTION_LIMIT!r}
LIANA_METHOD_PROFILES = {LIANA_METHOD_PROFILES!r}
LIANA_BY_SAMPLE_SIGNATURE = {LIANA_BY_SAMPLE_SIGNATURE!r}
LIANA_METHOD_SIGNATURES = {LIANA_METHOD_SIGNATURES!r}
LIANA_RESULT_COLUMNS = {LIANA_RESULT_COLUMNS!r}

{helper_source}

_OPENBIO_LIANA_PARAMETERS = {payload!r}

def run_liana_communication(adata, resource_path=None, liana_module=None):
    """Reproduce the OpenBio LIANA analysis and return its portable table and strict summary."""
    parameters = copy.deepcopy(_OPENBIO_LIANA_PARAMETERS)
    pinned_path = parameters.pop("resource_path")
    parameters["resource_path"] = pinned_path if resource_path is None else resource_path
    parameters["liana_module"] = liana_module
    table, _provenance, summary = _liana_run_impl(adata, **parameters)
    return table, summary
'''
    compile(code, "<openbio_liana_communication_code>", "exec")
    return code


__all__ = [
    "LIANA_AUDITED_VERSION",
    "LIANA_COMMUNICATION_NODE_ID",
    "LIANA_METHOD_PROFILES",
    "LIANA_RESOURCE_METADATA_FIELDS",
    "liana_communication_code",
    "liana_resource_cache_fingerprint",
    "run_liana_communication",
]
