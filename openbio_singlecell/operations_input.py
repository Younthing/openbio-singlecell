from __future__ import annotations

import gzip
import ntpath
import os
import posixpath
import time
import warnings
from collections.abc import Mapping, Sequence
from itertools import islice
from pathlib import Path
from textwrap import dedent
from typing import Any

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import validate_count_expression
from .artifact_codecs import PLOT_CODEC, TABLE_CODEC, read_anndata, read_table, write_anndata, write_plot, write_table
from .artifact_envelope import result_metadata
from .contracts import ensure_metadata
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

ANNDATA_KIND = "OPENBIO_ANNDATA"
ANNDATA_CODEC = "anndata-h5ad-v1"
TENX_FILE_CANDIDATES = {
    "matrix": ("matrix.mtx", "matrix.mtx.gz"),
    "barcodes": ("barcodes.tsv", "barcodes.tsv.gz"),
    "features": ("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"),
}
SLOT_ITEM_LIMIT = 64
SLOT_NAME_CHARACTER_LIMIT = 256
INDEX_EXAMPLE_LIMIT = 8
METADATA_WARNING_LIMIT = 32
ANNDATA_REFERENCE = AnalysisReference(
    citation=(
        "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data "
        "matrices. Journal of Open Source Software. 2024;9(101):4371."
    ),
    doi="10.21105/joss.04371",
    url="https://doi.org/10.21105/joss.04371",
    kind="software",
)
ANNDATA_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="AnnData documentation: annotated data matrix structure and storage slots.",
    url="https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html",
    kind="software_documentation",
)


def _require_exact_fields(value: dict[str, JSONValue], expected: set[str], description: str) -> None:
    if set(value) != expected:
        raise ProtocolError(f"{description} descriptor fields must be exactly {sorted(expected)!r}.")


def require_input_names(inputs: dict[str, JSONValue], expected: set[str], *, operation: str) -> None:
    if set(inputs) != expected:
        raise ProtocolError(f"{operation} inputs must be exactly {sorted(expected)!r}.")


def _require_input_descriptor(inputs: dict[str, JSONValue], name: str) -> dict[str, JSONValue]:
    if name not in inputs:
        raise ProtocolError(f"Worker inputs are missing {name!r}.")
    descriptor = inputs[name]
    if not isinstance(descriptor, dict):
        raise ProtocolError(f"Worker input {name!r} must be a descriptor object.")
    return descriptor


def _absolute_path(value: JSONValue, *, description: str, kind: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ProtocolError(f"{description} path must be a non-empty string.")
    if not os.path.isabs(value):
        raise ProtocolError(f"{description} path must be absolute.")
    path = Path(value).resolve(strict=True)
    if kind == "file" and not path.is_file():
        raise ProtocolError(f"{description} path must identify a file.")
    if kind == "directory" and not path.is_dir():
        raise ProtocolError(f"{description} path must identify a directory.")
    return path


def require_file_input(inputs: dict[str, JSONValue], name: str) -> tuple[Path, dict[str, JSONValue]]:
    descriptor = _require_input_descriptor(inputs, name)
    _require_exact_fields(descriptor, {"type", "path", "provenance"}, f"Worker input {name!r}")
    if descriptor["type"] != "file":
        raise ProtocolError(f"Worker input {name!r} must use the file descriptor type.")
    path = _absolute_path(descriptor["path"], description=f"Worker input {name!r}", kind="file")
    provenance = descriptor["provenance"]
    if not isinstance(provenance, dict):
        raise ProtocolError(f"Worker input {name!r} provenance must be a JSON object.")
    return path, provenance


def require_directory_input(inputs: dict[str, JSONValue], name: str) -> tuple[Path, dict[str, JSONValue]]:
    descriptor = _require_input_descriptor(inputs, name)
    _require_exact_fields(descriptor, {"type", "path", "provenance"}, f"Worker input {name!r}")
    if descriptor["type"] != "directory":
        raise ProtocolError(f"Worker input {name!r} must use the directory descriptor type.")
    path = _absolute_path(descriptor["path"], description=f"Worker input {name!r}", kind="directory")
    provenance = descriptor["provenance"]
    if not isinstance(provenance, dict):
        raise ProtocolError(f"Worker input {name!r} provenance must be a JSON object.")
    return path, provenance


def require_artifact_input(
    inputs: dict[str, JSONValue],
    name: str,
    *,
    kind: str,
    codec: str,
) -> Path:
    descriptor = _require_input_descriptor(inputs, name)
    _require_exact_fields(descriptor, {"type", "path", "kind", "codec"}, f"Worker input {name!r}")
    if descriptor["type"] != "artifact":
        raise ProtocolError(f"Worker input {name!r} must use the artifact descriptor type.")
    if descriptor["kind"] != kind or descriptor["codec"] != codec:
        raise ProtocolError(f"Worker input {name!r} must be {kind} with codec {codec}.")
    return _absolute_path(descriptor["path"], description=f"Worker input {name!r}", kind="directory")


def require_parameters(parameters: dict[str, JSONValue], expected: set[str], *, operation: str) -> None:
    if set(parameters) != expected:
        raise ProtocolError(f"{operation} parameters must be exactly {sorted(expected)!r}.")


def read_anndata_input(inputs: dict[str, JSONValue], name: str = "adata") -> Any:
    return read_anndata(require_artifact_input(inputs, name, kind=ANNDATA_KIND, codec=ANNDATA_CODEC))


def read_table_input(
    inputs: dict[str, JSONValue],
    name: str,
    *,
    kind: str,
    codec: str = TABLE_CODEC,
) -> tuple[Any, dict[str, Any]]:
    return read_table(require_artifact_input(inputs, name, kind=kind, codec=codec))


def _artifact_output(
    context: OperationContext,
    root: Path,
    *,
    name: str,
    kind: str,
    codec: str,
) -> dict[str, JSONValue]:
    return {
        "type": "artifact",
        "name": name,
        "kind": kind,
        "codec": codec,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def write_anndata_output(context: OperationContext, adata: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("adata")
    write_anndata(root, adata)
    return _artifact_output(context, root, name="adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC)


def write_table_output(
    context: OperationContext,
    result: Any,
    *,
    kind: str,
    name: str = "table",
    table: Any | None = None,
    json_list_columns: Sequence[str] = (),
) -> dict[str, JSONValue]:
    root = context.create_output_directory(name)
    write_table(
        root,
        result.table if table is None else table,
        result_metadata(result),
        json_list_columns=json_list_columns,
    )
    return _artifact_output(context, root, name=name, kind=kind, codec=TABLE_CODEC)


def write_plot_output(
    context: OperationContext,
    result: Any,
    *,
    kind: str,
    name: str = "plot",
) -> dict[str, JSONValue]:
    root = context.create_output_directory(name)
    write_plot(root, result.png, result_metadata(result))
    return _artifact_output(context, root, name=name, kind=kind, codec=PLOT_CODEC)


def analysis_outputs(report: Any, code: str, *outputs: JSONValue) -> list[JSONValue]:
    return [
        *outputs,
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


def _sanitize_embedded_provenance(value: Any) -> Any:
    if isinstance(value, (bytes, bytearray)):
        return _sanitize_embedded_provenance(bytes(value).decode("utf-8", errors="replace"))
    if not isinstance(value, (str, bytes, bytearray)) and hasattr(value, "tolist"):
        try:
            return _sanitize_embedded_provenance(value.tolist())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _sanitize_embedded_provenance(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_sanitize_embedded_provenance(item) for item in value]
    if isinstance(value, str):
        text = value.strip()
        if ntpath.isabs(text) or posixpath.isabs(text) or text.lower().startswith("file://"):
            return "<redacted-absolute-path>"
    return value


def _blank_index_count(index: Any) -> int:
    return sum(not str(value).strip() for value in index)


def _duplicate_index_count(index: Any) -> int:
    return int(len(index) - index.nunique())


def _validate_axis_names(adata: Any, *, make_var_names_unique: bool) -> dict[str, JSONValue]:
    blank_obs = _blank_index_count(adata.obs_names)
    duplicate_obs = _duplicate_index_count(adata.obs_names)
    blank_vars = _blank_index_count(adata.var_names)
    duplicate_vars = _duplicate_index_count(adata.var_names)
    if duplicate_vars and make_var_names_unique:
        adata.var_names_make_unique()
    return {
        "empty_observation_axis": int(adata.n_obs) == 0,
        "empty_variable_axis": int(adata.n_vars) == 0,
        "obs_names_unique": duplicate_obs == 0,
        "blank_obs_names": blank_obs,
        "duplicate_obs_name_occurrences": duplicate_obs,
        "blank_var_names": blank_vars,
        "var_names_unique_before_repair": duplicate_vars == 0,
        "duplicate_var_name_occurrences": duplicate_vars,
        "var_names_repaired": duplicate_vars if make_var_names_unique else 0,
        "var_names_unique_after_repair": bool(adata.var_names.is_unique),
    }


def _axis_advisories(audit: Mapping[str, Any], *, source_label: str = "H5AD input") -> list[str]:
    advisories = []
    if audit["empty_observation_axis"]:
        advisories.append(f"{source_label} has an empty observation axis.")
    if audit["empty_variable_axis"]:
        advisories.append(f"{source_label} has an empty variable axis.")
    if audit["blank_obs_names"]:
        advisories.append(f"{source_label} contains {audit['blank_obs_names']} blank observation name(s).")
    if audit["duplicate_obs_name_occurrences"]:
        advisories.append(
            f"{source_label} contains {audit['duplicate_obs_name_occurrences']} duplicate observation-name "
            "occurrence(s); later name-based alignment may be ambiguous."
        )
    if audit["blank_var_names"]:
        advisories.append(f"{source_label} contains {audit['blank_var_names']} blank variable name(s).")
    if not audit["var_names_unique_after_repair"]:
        advisories.append(
            f"{source_label} retains {audit['duplicate_var_name_occurrences']} duplicate variable-name "
            "occurrence(s) by explicit request."
        )
    return advisories


def _resolve_10x_mtx_files(directory: Path) -> dict[str, Path]:
    selected = {}
    for role, candidates in TENX_FILE_CANDIDATES.items():
        matches = [directory / filename for filename in candidates if (directory / filename).is_file()]
        if not matches:
            raise FileNotFoundError(f"10x directory is missing {role}; expected one of: {', '.join(candidates)}")
        if len(matches) > 1:
            raise ValueError(
                f"10x directory contains ambiguous {role} files: {', '.join(path.name for path in matches)}"
            )
        selected[role] = matches[0]
    return selected


def _validate_10x_options(var_names: object, make_unique: object, gex_only: object) -> tuple[str, bool, bool]:
    if var_names not in {"gene_symbols", "gene_ids"}:
        raise ProtocolError("10x var_names must be 'gene_symbols' or 'gene_ids'.")
    if not isinstance(make_unique, bool) or not isinstance(gex_only, bool):
        raise ProtocolError("10x make_unique and gex_only must be booleans.")
    return str(var_names), make_unique, gex_only


def _matrix_provenance(matrix: Any) -> dict[str, JSONValue]:
    science = dependencies.require_scientific_dependencies()
    storage = f"sparse_{matrix.format}" if science.sparse.issparse(matrix) else "dense"
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    values = science.np.asarray(values)
    finite = bool(values.size == 0 or science.np.isfinite(values).all())
    nonnegative = bool(values.size == 0 or (values >= 0).all())
    integer_like = bool(values.size == 0 or science.np.allclose(values, science.np.rint(values), rtol=0.0, atol=1e-8))
    return {
        "shape": [int(value) for value in matrix.shape],
        "dtype": str(matrix.dtype),
        "storage": storage,
        "finite": finite,
        "nonnegative": nonnegative,
        "integer_like": integer_like,
        "has_positive_values": bool(values.size and (values > 0).any()),
        "finite_nonnegative_integer_counts": finite and nonnegative and integer_like,
    }


def _validate_10x_payload(adata: Any, *, source_label: str) -> list[str]:
    science = dependencies.require_scientific_dependencies()
    if "gene_ids" not in adata.var:
        raise ValueError(f"{source_label} does not provide the required stable feature IDs in var['gene_ids'].")
    gene_ids = adata.var["gene_ids"]
    missing_gene_ids = int(gene_ids.isna().sum())
    blank_gene_ids = sum(not str(value).strip() for value in gene_ids if not bool(science.pd.isna(value)))
    duplicate_gene_ids = int(len(gene_ids) - gene_ids.astype(str).nunique())
    advisories = validate_count_expression(
        adata.X,
        source_label=f"{source_label} count matrix",
        require_integers=True,
        require_positive=False,
    )
    if missing_gene_ids or blank_gene_ids:
        advisories.append(
            f"{source_label} contains {missing_gene_ids + blank_gene_ids} missing or blank stable feature ID(s)."
        )
    if duplicate_gene_ids:
        advisories.append(f"{source_label} contains {duplicate_gene_ids} duplicate stable feature ID occurrence(s).")
    return advisories


def _read_10x_mtx(
    directory: Path,
    var_names: str,
    make_unique: bool,
    gex_only: bool,
) -> tuple[Any, dict[str, JSONValue]]:
    science = dependencies.require_scientific_dependencies()
    files = _resolve_10x_mtx_files(directory)
    opener = gzip.open if files["matrix"].name.lower().endswith(".gz") else open
    with opener(files["matrix"], "rb") as matrix_file:
        matrix = science.sparse.csr_matrix(science.mmread(matrix_file))
    barcodes = science.pd.read_csv(
        files["barcodes"],
        sep="\t",
        header=None,
        dtype=str,
        compression="infer",
        keep_default_na=False,
        skip_blank_lines=False,
    )
    features = science.pd.read_csv(
        files["features"],
        sep="\t",
        header=None,
        dtype=str,
        compression="infer",
        keep_default_na=False,
        skip_blank_lines=False,
    )
    if barcodes.shape[1] != 1:
        raise ValueError("The 10x barcodes file must contain exactly one column.")
    if features.shape[1] < 2:
        raise ValueError("The 10x features file must contain at least feature ID and feature name columns.")
    barcode_values = barcodes.iloc[:, 0].astype(str)
    gene_ids = features.iloc[:, 0].astype(str)
    gene_symbols = features.iloc[:, 1].astype(str)
    feature_types = features.iloc[:, 2].astype(str) if features.shape[1] > 2 else None
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            "10x matrix dimensions do not match the features and barcodes files: "
            f"matrix={matrix.shape}, features={len(features)}, barcodes={len(barcodes)}."
        )

    original_feature_count = int(len(features))
    feature_type_counts = (
        {str(key): int(value) for key, value in feature_types.value_counts(sort=False).items()}
        if feature_types is not None
        else {}
    )
    if gex_only and feature_types is not None:
        mask = feature_types.to_numpy() == "Gene Expression"
        matrix = matrix[mask, :]
        gene_ids = gene_ids[mask].reset_index(drop=True)
        gene_symbols = gene_symbols[mask].reset_index(drop=True)
        feature_types = feature_types[mask].reset_index(drop=True)
    selected_names = gene_ids if var_names == "gene_ids" else gene_symbols
    selected_index = science.pd.Index(selected_names.astype(str).to_numpy())
    duplicate_var_names = int(len(selected_index) - selected_index.nunique())
    if duplicate_var_names and make_unique:
        selected_index = science.ad.utils.make_index_unique(selected_index)
    var = science.pd.DataFrame(index=selected_index)
    var["gene_ids"] = gene_ids.astype(str).to_numpy()
    var["gene_symbols"] = gene_symbols.astype(str).to_numpy()
    if feature_types is not None:
        var["feature_types"] = feature_types.astype(str).to_numpy()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
        adata = science.ad.AnnData(
            X=matrix.transpose().tocsr(),
            obs=science.pd.DataFrame(index=barcode_values.to_numpy()),
            var=var,
        )
    axis_audit = _validate_axis_names(adata, make_var_names_unique=False)
    axis_audit.update(
        {
            "var_names_unique_before_repair": duplicate_var_names == 0,
            "duplicate_var_name_occurrences": duplicate_var_names,
            "var_names_repaired": duplicate_var_names if make_unique else 0,
            "var_names_unique_after_repair": bool(adata.var_names.is_unique),
        }
    )
    blank_barcodes = sum(not value.strip() for value in barcode_values)
    duplicate_barcodes = int(len(barcode_values) - barcode_values.nunique())
    blank_gene_ids = sum(not value.strip() for value in gene_ids)
    duplicate_gene_ids = int(len(gene_ids) - gene_ids.nunique())
    blank_symbols = sum(not value.strip() for value in gene_symbols)
    payload_advisories = _validate_10x_payload(adata, source_label="10x MTX input")
    details: dict[str, JSONValue] = {
        "reader_parameters": {"var_names": var_names, "make_unique": make_unique, "gex_only": gex_only},
        "input_features": original_feature_count,
        "retained_features": int(adata.n_vars),
        "cells": int(adata.n_obs),
        "feature_types_before_filter": feature_type_counts,
        "gex_filter_verified": feature_types is not None,
        "axis_names": axis_audit,
        "matrix": _matrix_provenance(adata.X),
        "warnings": [
            *(
                ["The legacy features file has no feature-type column; gex_only could not be verified."]
                if gex_only and feature_types is None
                else []
            ),
            *([f"The 10x barcodes file contains {blank_barcodes} blank barcode(s)."] if blank_barcodes else []),
            *(
                [
                    f"The 10x barcodes file contains {duplicate_barcodes} duplicate barcode occurrence(s); "
                    "later name-based cell alignment may be ambiguous."
                ]
                if duplicate_barcodes
                else []
            ),
            *(
                [f"The retained 10x features contain {blank_gene_ids} blank stable feature ID(s)."]
                if blank_gene_ids
                else []
            ),
            *(
                [f"The retained 10x features contain {duplicate_gene_ids} duplicate stable feature ID occurrence(s)."]
                if duplicate_gene_ids
                else []
            ),
            *([f"The retained 10x features contain {blank_symbols} blank feature name(s)."] if blank_symbols else []),
            *_axis_advisories(axis_audit, source_label="10x MTX input"),
            *payload_advisories,
        ],
    }
    return adata, details


def _discover_10x_study(root: Path) -> list[tuple[str, Path]]:
    children = sorted(root.iterdir(), key=lambda path: path.name.casefold())
    directories = [path for path in children if path.is_dir() or path.is_symlink()]
    if not directories:
        raise ValueError("The study directory does not contain any first-level Sample directories.")
    discovered = []
    errors = []
    for child in directories:
        if child.is_symlink():
            errors.append(f"{child.name}: symbolic-link Sample directories are not allowed")
            continue
        resolved = child.resolve(strict=True)
        if not resolved.is_relative_to(root):
            errors.append(f"{child.name}: Sample directory escapes the selected Study directory")
            continue
        try:
            _resolve_10x_mtx_files(resolved)
        except (OSError, ValueError) as error:
            errors.append(f"{child.name}: {error}")
            continue
        discovered.append((child.name, resolved))
    if errors:
        raise ValueError("Invalid 10x Study Sample inventory:\n" + "\n".join(f"- {item}" for item in errors))
    return discovered


def _validate_study_options(
    sample_key: object,
    var_names: object,
    join: object,
    make_unique: object,
    gex_only: object,
) -> tuple[str, str, str, bool, bool]:
    if not isinstance(sample_key, str) or not sample_key.strip():
        raise ProtocolError("Study sample_key cannot be empty.")
    if join not in {"inner", "outer"}:
        raise ProtocolError("Study join must be 'inner' or 'outer'.")
    validated_var_names, validated_make_unique, validated_gex_only = _validate_10x_options(
        var_names, make_unique, gex_only
    )
    return sample_key.strip(), validated_var_names, str(join), validated_make_unique, validated_gex_only


def _study_feature_metadata(samples: Mapping[str, Any]) -> dict[str, dict[str, str]]:
    science = dependencies.require_scientific_dependencies()
    observed: dict[str, tuple[str, str]] = {}
    conflicts = []
    metadata: dict[str, dict[str, str]] = {}
    for sample_name, sample in samples.items():
        if not sample.var_names.is_unique:
            raise ValueError(
                f"Sample {sample_name!r} has duplicate selected feature names; "
                "Study concatenation requires make_unique=True or gene_ids."
            )
        for var_name, gene_id in zip(sample.var_names.astype(str), sample.var["gene_ids"].astype(str), strict=True):
            previous = observed.get(var_name)
            if previous is not None and previous[1] != gene_id:
                conflicts.append(f"{var_name!r}: {previous[0]}={previous[1]!r}, {sample_name}={gene_id!r}")
            else:
                observed[var_name] = (sample_name, gene_id)
        for position, var_name in enumerate(sample.var_names.astype(str)):
            row = sample.var.iloc[position]
            target = metadata.setdefault(var_name, {})
            for column in ("gene_ids", "gene_symbols", "feature_types"):
                if column not in sample.var or bool(science.pd.isna(row[column])):
                    continue
                value = str(row[column])
                previous = target.get(column)
                if previous is not None and previous != value:
                    conflicts.append(
                        f"{var_name!r} has conflicting {column}: {previous!r} versus {value!r} in {sample_name}"
                    )
                else:
                    target[column] = value
    if conflicts:
        preview = "; ".join(conflicts[:8])
        suffix = f"; plus {len(conflicts) - 8} more" if len(conflicts) > 8 else ""
        raise ValueError(f"Selected features have conflicting identity metadata across Samples: {preview}{suffix}")
    return metadata


@register_operation("openbio.node.load10xstudy")
def load_10x_study(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"directory"}, operation="Load 10x Study")
    directory, provenance = require_directory_input(inputs, "directory")
    require_parameters(
        parameters,
        {"sample_key", "var_names", "join", "make_unique", "gex_only"},
        operation="Load 10x Study",
    )
    sample_key, var_names, join, make_unique, gex_only = _validate_study_options(
        parameters["sample_key"],
        parameters["var_names"],
        parameters["join"],
        parameters["make_unique"],
        parameters["gex_only"],
    )
    discovered = _discover_10x_study(directory)
    sample_provenance = provenance.get("samples")
    sample_names = [name for name, _ in discovered]
    if not isinstance(sample_provenance, dict) or list(sample_provenance) != sample_names:
        raise ProtocolError("Load 10x Study provenance samples must exactly match sorted discovered Samples.")
    if any(not isinstance(sample_provenance[name], dict) for name in sample_names):
        raise ProtocolError("Load 10x Study Sample provenance values must be JSON objects.")

    loaded = {}
    errors = []
    for sample_name, sample_directory in discovered:
        try:
            loaded[sample_name] = _read_10x_mtx(sample_directory, var_names, make_unique, gex_only)
        except (EOFError, OSError, TypeError, ValueError) as error:
            errors.append(f"{sample_name}: {error}")
    if errors:
        raise ValueError("Invalid 10x Study Sample payloads:\n" + "\n".join(f"- {item}" for item in errors))
    samples = {sample_name: value[0] for sample_name, value in loaded.items()}
    feature_metadata = _study_feature_metadata(samples)
    feature_sets = {sample_name: set(sample.var_names.astype(str)) for sample_name, sample in samples.items()}
    common_features = set.intersection(*feature_sets.values())
    union_features = set.union(*feature_sets.values())
    science = dependencies.require_scientific_dependencies()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
        adata = science.ad.concat(
            samples,
            join=join,
            merge="same",
            label=sample_key,
            index_unique="-",
            fill_value=0,
        )
    for column in ("gene_ids", "gene_symbols", "feature_types"):
        values = [feature_metadata[str(name)].get(column) for name in adata.var_names]
        if any(value is not None for value in values):
            adata.var[column] = values
    retained_features = common_features if join == "inner" else union_features
    alignment = {
        sample_name: {
            "input_features": len(features),
            "dropped_features": len(features - retained_features) if join == "inner" else 0,
            "zero_filled_features": len(retained_features - features) if join == "outer" else 0,
        }
        for sample_name, features in feature_sets.items()
    }
    study_warnings = []
    if join == "inner" and not common_features:
        study_warnings.append(
            "join='inner' produced an empty feature intersection; the returned Study has zero variables."
        )
    if join == "outer" and any(item["zero_filled_features"] for item in alignment.values()):
        study_warnings.append(
            "join='outer' represents features absent from a Sample as sparse zeros; these are structural "
            "missingness, not measured biological zeros."
        )
    if not adata.obs_names.is_unique:
        study_warnings.append(
            "Study concatenation retained duplicate observation names; later name-based cell alignment may be ambiguous."
        )
    if int(adata.n_obs) == 0:
        study_warnings.append("The concatenated Study has an empty observation axis.")
    root_path = provenance.get("path")
    if not isinstance(root_path, str) or not root_path:
        raise ProtocolError("Load 10x Study provenance path must be a non-empty string.")
    ensure_metadata(
        adata,
        display_name=directory.name or "10x study",
        source={
            "kind": "10x_study",
            "path": root_path,
            "reader_parameters": {
                "sample_key": sample_key,
                "var_names": var_names,
                "join": join,
                "make_unique": make_unique,
                "gex_only": gex_only,
            },
            "sample_semantics": (
                "Each first-level directory is asserted by the user to be one independent biological Sample."
            ),
            "sample_count": len(samples),
            "common_features": len(common_features),
            "union_features": len(union_features),
            "retained_features": len(retained_features),
            "feature_alignment": alignment,
            "index_unique_delimiter": "-",
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "warnings": study_warnings,
            "samples": {name: {**sample_provenance[name], **loaded[name][1]} for name in sample_names},
        },
    )
    return [write_anndata_output(context, adata)]


@register_operation("openbio.node.load10xh5")
def load_10x_h5(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"path"}, operation="Load 10x H5")
    path, provenance = require_file_input(inputs, "path")
    if not path.name.lower().endswith((".h5", ".hdf5")):
        raise ProtocolError("Load 10x H5 requires an .h5 or .hdf5 input file.")
    require_parameters(parameters, {"genome", "gex_only", "make_unique"}, operation="Load 10x H5")
    genome = parameters["genome"]
    gex_only = parameters["gex_only"]
    make_unique = parameters["make_unique"]
    if not isinstance(genome, str):
        raise ProtocolError("Load 10x H5 genome must be a string.")
    if not isinstance(gex_only, bool) or not isinstance(make_unique, bool):
        raise ProtocolError("Load 10x H5 gex_only and make_unique must be booleans.")
    normalized_genome = genome.strip()
    science = dependencies.require_scientific_dependencies()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
        adata = science.sc.read_10x_h5(
            path,
            genome=normalized_genome or None,
            gex_only=False,
            backup_url=None,
        )
        original_feature_count = int(adata.n_vars)
        feature_type_counts_before = (
            {str(key): int(value) for key, value in adata.var["feature_types"].value_counts(sort=False).items()}
            if "feature_types" in adata.var
            else {}
        )
        if gex_only and "feature_types" in adata.var:
            # The feature filter changes the scientific axes, so this slice must be materialized.
            adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    axis_audit = _validate_axis_names(adata, make_var_names_unique=make_unique)
    payload_advisories = _validate_10x_payload(adata, source_label="10x H5 input")
    feature_type_counts = (
        {str(key): int(value) for key, value in adata.var["feature_types"].value_counts(sort=False).items()}
        if "feature_types" in adata.var
        else {}
    )
    source: dict[str, Any] = dict(provenance)
    source.update(
        {
            "reader_parameters": {
                "genome": normalized_genome or None,
                "gex_only": gex_only,
                "make_unique": make_unique,
                "backup_url": None,
            },
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "input_features": original_feature_count,
            "retained_features": int(adata.n_vars),
            "filtered_features": original_feature_count - int(adata.n_vars),
            "feature_types_before_filter": feature_type_counts_before,
            "feature_types_retained": feature_type_counts,
            "axis_names": axis_audit,
            "matrix": _matrix_provenance(adata.X),
            "warnings": [
                *_axis_advisories(axis_audit, source_label="10x H5 input"),
                *payload_advisories,
            ],
        }
    )
    ensure_metadata(
        adata,
        display_name=path.stem,
        source={"kind": "10x_h5", **source},
    )
    return [write_anndata_output(context, adata)]


@register_operation("openbio.node.load10xmtx")
def load_10x_mtx(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"directory"}, operation="Load 10x MTX")
    directory, provenance = require_directory_input(inputs, "directory")
    require_parameters(parameters, {"var_names", "make_unique", "gex_only"}, operation="Load 10x MTX")
    var_names, make_unique, gex_only = _validate_10x_options(
        parameters["var_names"], parameters["make_unique"], parameters["gex_only"]
    )
    adata, details = _read_10x_mtx(directory, var_names, make_unique, gex_only)
    source: dict[str, Any] = dict(provenance)
    source.update(details)
    ensure_metadata(
        adata,
        display_name=directory.name or "10x",
        source={"kind": "10x_mtx", **source},
    )
    return [write_anndata_output(context, adata)]


def _bounded_name(value: Any) -> str:
    text = str(value)
    return text if len(text) <= SLOT_NAME_CHARACTER_LIMIT else f"{text[: SLOT_NAME_CHARACTER_LIMIT - 3]}..."


def _bounded_names(values: Any) -> dict[str, Any]:
    total = int(len(values))
    names = [_bounded_name(value) for value in islice(iter(values), SLOT_ITEM_LIMIT)]
    return {"names": names, "total": total, "truncated": total > SLOT_ITEM_LIMIT}


def _array_descriptor(value: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    if value is None:
        return {"present": False, "shape": None, "dtype": None, "python_type": None, "storage": "none"}
    shape = getattr(value, "shape", None)
    dtype = getattr(value, "dtype", None)
    if science.sparse.issparse(value):
        storage = f"sparse_{value.format}"
    elif getattr(value, "format", None) in {"csr", "csc"}:
        storage = f"sparse_backed_{value.format}"
    elif type(value).__module__.startswith("numpy"):
        storage = "dense"
    elif shape is not None:
        storage = "array_like"
    else:
        storage = "object"
    return {
        "present": True,
        "shape": [int(item) for item in shape] if shape is not None else None,
        "dtype": str(dtype) if dtype is not None else None,
        "python_type": f"{type(value).__module__}.{type(value).__qualname__}",
        "storage": storage,
    }


def _bounded_columns(frame: Any) -> dict[str, Any]:
    result = _bounded_names(frame.columns)
    result["dtypes"] = [str(dtype) for dtype in islice(iter(frame.dtypes), SLOT_ITEM_LIMIT)]
    return result


def _bounded_mapping(mapping: Any, *, exclude_none: bool = False) -> dict[str, Any]:
    total = sum(key is not None or not exclude_none for key in mapping.keys())
    keys = list(islice((key for key in mapping.keys() if key is not None or not exclude_none), SLOT_ITEM_LIMIT))
    descriptors = [_array_descriptor(mapping[key]) for key in keys]
    return {
        "names": [_bounded_name(key) for key in keys],
        "total": total,
        "truncated": total > SLOT_ITEM_LIMIT,
        "shapes": [item["shape"] for item in descriptors],
        "dtypes": [item["dtype"] for item in descriptors],
        "storage": [item["storage"] for item in descriptors],
    }


def _empty_column_inventory() -> dict[str, Any]:
    return {"names": [], "total": 0, "truncated": False, "dtypes": []}


def _empty_mapping_inventory() -> dict[str, Any]:
    return {"names": [], "total": 0, "truncated": False, "shapes": [], "dtypes": [], "storage": []}


def _index_descriptor(index: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    return {
        "name": str(index.name) if index.name is not None else None,
        "dtype": str(index.dtype),
        "is_unique": bool(index.is_unique),
        "missing": int(science.pd.isna(index).sum()),
        "blank": _blank_index_count(index),
        "examples": [_bounded_name(value) for value in islice(iter(index), INDEX_EXAMPLE_LIMIT)],
    }


def _metadata_warning_values(metadata: Mapping[str, Any] | None) -> list[Any] | None:
    if metadata is None:
        return None
    values = metadata.get("warnings")
    if isinstance(values, (list, tuple)):
        return list(values)
    if hasattr(values, "tolist"):
        converted = values.tolist()
        if isinstance(converted, list):
            return converted
    return None


def _anndata_summary_key_results(adata: Any) -> dict[str, Any]:
    raw = adata.raw
    metadata = adata.uns.get("openbio_singlecell")
    valid_metadata = isinstance(metadata, Mapping)
    source = metadata.get("source") if valid_metadata else None
    history = metadata.get("analysis_history") if valid_metadata else None
    metadata_warnings = _metadata_warning_values(metadata if valid_metadata else None)
    return {
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "object_state": {"is_view": bool(adata.is_view), "is_backed": bool(adata.isbacked)},
        "X": _array_descriptor(adata.X),
        "obs_index": _index_descriptor(adata.obs_names),
        "var_index": _index_descriptor(adata.var_names),
        "raw": {
            "present": raw is not None,
            "shape": [int(raw.n_obs), int(raw.n_vars)] if raw is not None else None,
            "X": _array_descriptor(raw.X if raw is not None else None),
            "var_index": _index_descriptor(raw.var_names) if raw is not None else None,
        },
        "obs_columns": _bounded_columns(adata.obs),
        "var_columns": _bounded_columns(adata.var),
        "layers": _bounded_mapping(adata.layers, exclude_none=True),
        "obsm": _bounded_mapping(adata.obsm),
        "varm": _bounded_mapping(adata.varm),
        "obsp": _bounded_mapping(adata.obsp),
        "varp": _bounded_mapping(adata.varp),
        "uns": _bounded_mapping(adata.uns),
        "raw_var_columns": _bounded_columns(raw.var) if raw is not None else _empty_column_inventory(),
        "raw_varm": _bounded_mapping(raw.varm) if raw is not None else _empty_mapping_inventory(),
        "openbio_metadata": {
            "present": metadata is not None,
            "valid_mapping": valid_metadata,
            "display_name": (
                _bounded_name(metadata.get("display_name")) if valid_metadata and metadata.get("display_name") else None
            ),
            "source_kind": (
                _bounded_name(source.get("kind")) if isinstance(source, Mapping) and source.get("kind") else None
            ),
            "history_entries": len(history) if isinstance(history, Mapping) else None,
            "warning_count": len(metadata_warnings) if metadata_warnings is not None else None,
        },
    }


def _anndata_summary_code() -> str:
    return dedent(
        f"""
        from collections.abc import Mapping
        from itertools import islice

        import pandas as pd
        from scipy import sparse

        SLOT_ITEM_LIMIT = {SLOT_ITEM_LIMIT}
        SLOT_NAME_CHARACTER_LIMIT = {SLOT_NAME_CHARACTER_LIMIT}
        INDEX_EXAMPLE_LIMIT = {INDEX_EXAMPLE_LIMIT}

        def _bounded_name(value):
            text = str(value)
            return text if len(text) <= SLOT_NAME_CHARACTER_LIMIT else f"{{text[: SLOT_NAME_CHARACTER_LIMIT - 3]}}..."

        def _bounded_names(values):
            total = int(len(values))
            return {{
                "names": [_bounded_name(value) for value in islice(iter(values), SLOT_ITEM_LIMIT)],
                "total": total,
                "truncated": total > SLOT_ITEM_LIMIT,
            }}

        def _array_descriptor(value):
            if value is None:
                return {{"present": False, "shape": None, "dtype": None, "python_type": None, "storage": "none"}}
            shape = getattr(value, "shape", None)
            dtype = getattr(value, "dtype", None)
            if sparse.issparse(value):
                storage = f"sparse_{{value.format}}"
            elif getattr(value, "format", None) in {{"csr", "csc"}}:
                storage = f"sparse_backed_{{value.format}}"
            elif type(value).__module__.startswith("numpy"):
                storage = "dense"
            elif shape is not None:
                storage = "array_like"
            else:
                storage = "object"
            return {{
                "present": True,
                "shape": [int(item) for item in shape] if shape is not None else None,
                "dtype": str(dtype) if dtype is not None else None,
                "python_type": f"{{type(value).__module__}}.{{type(value).__qualname__}}",
                "storage": storage,
            }}

        def _bounded_columns(frame):
            result = _bounded_names(frame.columns)
            result["dtypes"] = [str(dtype) for dtype in islice(iter(frame.dtypes), SLOT_ITEM_LIMIT)]
            return result

        def _bounded_mapping(mapping, *, exclude_none=False):
            total = sum(key is not None or not exclude_none for key in mapping.keys())
            keys = list(islice((key for key in mapping.keys() if key is not None or not exclude_none), SLOT_ITEM_LIMIT))
            descriptors = [_array_descriptor(mapping[key]) for key in keys]
            return {{
                "names": [_bounded_name(key) for key in keys],
                "total": total,
                "truncated": total > SLOT_ITEM_LIMIT,
                "shapes": [item["shape"] for item in descriptors],
                "dtypes": [item["dtype"] for item in descriptors],
                "storage": [item["storage"] for item in descriptors],
            }}

        def _index_descriptor(index):
            return {{
                "name": str(index.name) if index.name is not None else None,
                "dtype": str(index.dtype),
                "is_unique": bool(index.is_unique),
                "missing": int(pd.isna(index).sum()),
                "blank": sum(not str(value).strip() for value in index),
                "examples": [_bounded_name(value) for value in islice(iter(index), INDEX_EXAMPLE_LIMIT)],
            }}

        def _empty_columns():
            return {{"names": [], "total": 0, "truncated": False, "dtypes": []}}

        def _empty_mapping():
            return {{"names": [], "total": 0, "truncated": False, "shapes": [], "dtypes": [], "storage": []}}

        def summarize_anndata(adata):
            raw = adata.raw
            metadata = adata.uns.get("openbio_singlecell")
            valid_metadata = isinstance(metadata, Mapping)
            source = metadata.get("source") if valid_metadata else None
            history = metadata.get("analysis_history") if valid_metadata else None
            return {{
                "shape": [int(adata.n_obs), int(adata.n_vars)],
                "object_state": {{"is_view": bool(adata.is_view), "is_backed": bool(adata.isbacked)}},
                "X": _array_descriptor(adata.X),
                "obs_index": _index_descriptor(adata.obs_names),
                "var_index": _index_descriptor(adata.var_names),
                "raw": {{
                    "present": raw is not None,
                    "shape": [int(raw.n_obs), int(raw.n_vars)] if raw is not None else None,
                    "X": _array_descriptor(raw.X if raw is not None else None),
                    "var_index": _index_descriptor(raw.var_names) if raw is not None else None,
                }},
                "obs_columns": _bounded_columns(adata.obs),
                "var_columns": _bounded_columns(adata.var),
                "layers": _bounded_mapping(adata.layers, exclude_none=True),
                "obsm": _bounded_mapping(adata.obsm),
                "varm": _bounded_mapping(adata.varm),
                "obsp": _bounded_mapping(adata.obsp),
                "varp": _bounded_mapping(adata.varp),
                "uns": _bounded_mapping(adata.uns),
                "raw_var_columns": _bounded_columns(raw.var) if raw is not None else _empty_columns(),
                "raw_varm": _bounded_mapping(raw.varm) if raw is not None else _empty_mapping(),
                "openbio_metadata": {{
                    "present": metadata is not None,
                    "valid_mapping": valid_metadata,
                    "display_name": _bounded_name(metadata.get("display_name")) if valid_metadata and metadata.get("display_name") else None,
                    "source_kind": _bounded_name(source.get("kind")) if isinstance(source, Mapping) and source.get("kind") else None,
                    "history_entries": len(history) if isinstance(history, Mapping) else None,
                    "warning_count": len(metadata.get("warnings", [])) if valid_metadata and isinstance(metadata.get("warnings", []), (list, tuple)) else None,
                }},
            }}
        """
    )


@register_operation("openbio.node.anndatasummary")
def anndata_summary(
    _context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="AnnData Summary")
    require_parameters(parameters, set(), operation="AnnData Summary")
    adata = read_anndata_input(inputs)
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    key_results = _anndata_summary_key_results(adata)
    metadata = adata.uns.get("openbio_singlecell")
    display_name = (
        _bounded_name(metadata.get("display_name", "AnnData")) if isinstance(metadata, Mapping) else "AnnData"
    )
    warnings_list = []
    if cells == 0 or genes == 0:
        warnings_list.append("The AnnData object has an empty observation or variable axis.")
    if not adata.obs_names.is_unique:
        warnings_list.append("Observation names are not unique; cross-object cell alignment is ambiguous.")
    if not adata.var_names.is_unique:
        warnings_list.append("Variable names are not unique; name-based feature selection is ambiguous.")
    if adata.raw is not None and not adata.raw.var_names.is_unique:
        warnings_list.append("Raw variable names are not unique.")
    if metadata is not None and not isinstance(metadata, Mapping):
        warnings_list.append("openbio_singlecell metadata exists but is not a mapping.")
    metadata_warnings = _metadata_warning_values(metadata if isinstance(metadata, Mapping) else None)
    if metadata_warnings is not None:
        warnings_list.extend(_bounded_name(value) for value in metadata_warnings[:METADATA_WARNING_LIMIT])
        if len(metadata_warnings) > METADATA_WARNING_LIMIT:
            warnings_list.append(f"OpenBio metadata warnings were limited to the first {METADATA_WARNING_LIMIT} items.")
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellAnnDataSummary",
        title=f"{display_name} summary",
        operation="anndata_summary",
        methods=(
            "Inspected AnnData axis identity, storage state, raw snapshot, and bounded slot inventories on the "
            "worker's private file-backed artifact without modifying the input artifact."
        ),
        results=(
            f"The object contains {cells:,} observations and {genes:,} variables, with "
            f"{sum(key is not None for key in adata.layers.keys()):,} expression layer(s), "
            f"{len(adata.obsm):,} observation embedding(s), and "
            f"{'an' if adata.raw is not None else 'no'} attached raw snapshot."
        ),
        key_results=key_results,
        parameters={
            "slot_item_limit": SLOT_ITEM_LIMIT,
            "slot_name_character_limit": SLOT_NAME_CHARACTER_LIMIT,
            "index_example_limit": INDEX_EXAMPLE_LIMIT,
            "metadata_warning_limit": METADATA_WARNING_LIMIT,
        },
        references=[ANNDATA_REFERENCE, ANNDATA_DOCUMENTATION_REFERENCE],
        software_packages=["anndata", "numpy", "pandas", "scipy"],
        warnings=warnings_list,
        limitations=[
            "This structural inventory does not determine whether expression values are counts, normalized, "
            "logged, scaled, corrected, or otherwise transformed.",
            "It does not assess QC quality, annotation validity, biological interpretation, or whether the "
            "Study design supports statistical inference.",
            f"Each slot inventory is limited to {SLOT_ITEM_LIMIT} names and each displayed name to "
            f"{SLOT_NAME_CHARACTER_LIMIT} characters.",
            f"At most {METADATA_WARNING_LIMIT} existing OpenBio metadata warnings are copied into this report.",
        ],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=_anndata_summary_code(),
    )
    return analysis_outputs(report, code)


def load_h5ad(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    from anndata import read_h5ad

    require_input_names(inputs, {"path"}, operation="Load H5AD")
    path, provenance = require_file_input(inputs, "path")
    if path.suffix.lower() != ".h5ad":
        raise ProtocolError("Load H5AD requires an .h5ad input file.")
    require_parameters(parameters, {"make_var_names_unique"}, operation="Load H5AD")
    make_var_names_unique = parameters["make_var_names_unique"]
    if not isinstance(make_var_names_unique, bool):
        raise ProtocolError("Load H5AD make_var_names_unique must be a boolean.")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
        adata = read_h5ad(path, backed=None)
    axis_audit = _validate_axis_names(adata, make_var_names_unique=make_var_names_unique)
    existing_metadata = adata.uns.get("openbio_singlecell")
    embedded_source = existing_metadata.get("source", {}) if isinstance(existing_metadata, Mapping) else {}
    metadata_advisory = None
    try:
        ensure_metadata(adata)
    except ValueError as error:
        quarantined_metadata = _sanitize_embedded_provenance(existing_metadata)
        del adata.uns["openbio_singlecell"]
        ensure_metadata(adata)
        metadata_advisory = (
            "Embedded OpenBio metadata was incompatible with the current metadata schema and was quarantined "
            f"instead of used as workflow evidence: {error}"
        )
    else:
        quarantined_metadata = None

    source: dict[str, Any] = dict(provenance)
    source_warnings = _axis_advisories(axis_audit)
    if metadata_advisory:
        source_warnings.append(metadata_advisory)
    source.update(
        {
            "reader_parameters": {"backed": None, "make_var_names_unique": make_var_names_unique},
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "axis_names": axis_audit,
            "warnings": source_warnings,
        }
    )
    if embedded_source:
        source["embedded_source"] = _sanitize_embedded_provenance(embedded_source)
    if quarantined_metadata is not None:
        source["quarantined_embedded_openbio_metadata"] = quarantined_metadata
    ensure_metadata(
        adata,
        display_name=path.stem,
        source={"kind": "h5ad", **source},
    )
    return [write_anndata_output(context, adata)]


__all__ = [
    "ANNDATA_CODEC",
    "ANNDATA_KIND",
    "analysis_outputs",
    "anndata_summary",
    "load_10x_h5",
    "load_10x_mtx",
    "load_10x_study",
    "load_h5ad",
    "require_artifact_input",
    "require_directory_input",
    "require_file_input",
    "require_input_names",
    "require_parameters",
    "read_anndata_input",
    "read_table_input",
    "write_anndata_output",
    "write_plot_output",
    "write_table_output",
]
