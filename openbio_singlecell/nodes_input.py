from __future__ import annotations

import gzip
import ntpath
import os
import posixpath
import time
import warnings
from collections.abc import Mapping
from itertools import islice
from textwrap import dedent
from typing import TYPE_CHECKING, Any

import folder_paths
from comfy_api.latest import io

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import validate_count_expression
from .contracts import ensure_metadata
from .files import (
    input_file_fingerprint,
    input_file_provenance,
    resolve_10x_mtx_files,
    resolve_input_path,
    tenx_mtx_fingerprint,
    tenx_mtx_provenance,
)
from .node_types import AnnDataType, SummaryResultType

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/input"
DIAGNOSTIC_CATEGORY = "openbio/single-cell/diagnostics"
INPUT_FILE_UPLOAD_WIDGET = "OPENBIO_INPUT_FILE_UPLOAD_WIDGET"
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


def _file_upload_widget(*extensions: str, label: str, drop_label: str) -> dict[str, Any]:
    return {
        "widgetType": INPUT_FILE_UPLOAD_WIDGET,
        "allowed_extensions": list(extensions),
        "accept": ",".join(extensions),
        "upload_subfolder": "openbio-singlecell",
        "upload_label": label,
        "drop_label": drop_label,
    }


def _source_metadata(kind: str, provenance: dict[str, Any]) -> dict[str, Any]:
    return {"kind": kind, **provenance}


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


def _validate_axis_names(
    adata: AnnData,
    *,
    source_label: str,
    make_var_names_unique: bool,
) -> dict[str, Any]:
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


def _axis_advisories(audit: Mapping[str, Any], *, source_label: str) -> list[str]:
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


def _validate_10x_payload(adata: AnnData, *, source_label: str) -> list[str]:
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
        advisories.append(
            f"{source_label} contains {duplicate_gene_ids} duplicate stable feature ID occurrence(s)."
        )
    return advisories


def _matrix_provenance(matrix: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    if science.sparse.issparse(matrix):
        storage = f"sparse_{matrix.format}"
    else:
        storage = "dense"
    values = matrix.data if science.sparse.issparse(matrix) else science.np.asarray(matrix).ravel()
    values = science.np.asarray(values)
    finite = bool(values.size == 0 or science.np.isfinite(values).all())
    nonnegative = bool(values.size == 0 or (values >= 0).all())
    integer_like = bool(
        values.size == 0
        or science.np.allclose(values, science.np.rint(values), rtol=0.0, atol=1e-8)
    )
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


def _validate_10x_options(var_names: str, make_unique: bool, gex_only: bool) -> None:
    if var_names not in {"gene_symbols", "gene_ids"}:
        raise ValueError("10x var_names must be 'gene_symbols' or 'gene_ids'.")
    if not isinstance(make_unique, bool):
        raise TypeError("10x make_unique must be a boolean.")
    if not isinstance(gex_only, bool):
        raise TypeError("10x gex_only must be a boolean.")


def _read_10x_mtx(
    relative_directory: str,
    var_names: str,
    make_unique: bool,
    gex_only: bool,
) -> tuple[AnnData, dict[str, Any]]:
    _validate_10x_options(var_names, make_unique, gex_only)
    science = dependencies.require_scientific_dependencies()
    files = resolve_10x_mtx_files(relative_directory)

    opener = gzip.open if files["matrix"].lower().endswith(".gz") else open
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
    blank_barcodes = sum(not value.strip() for value in barcode_values)
    duplicate_barcodes = int(len(barcode_values) - barcode_values.nunique())

    gene_ids = features.iloc[:, 0].astype(str)
    gene_symbols = features.iloc[:, 1].astype(str)
    feature_types = features.iloc[:, 2].astype(str) if features.shape[1] > 2 else None
    if matrix.shape[0] != len(features) or matrix.shape[1] != len(barcodes):
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
    gex_filter_verified = feature_types is not None
    if gex_only and feature_types is not None:
        mask = feature_types.to_numpy() == "Gene Expression"
        matrix = matrix[mask, :]
        gene_ids = gene_ids[mask].reset_index(drop=True)
        gene_symbols = gene_symbols[mask].reset_index(drop=True)
        feature_types = feature_types[mask].reset_index(drop=True)

    blank_gene_ids = sum(not value.strip() for value in gene_ids)
    duplicate_gene_ids = int(len(gene_ids) - gene_ids.nunique())
    blank_symbols = sum(not value.strip() for value in gene_symbols)

    selected_names = gene_ids if var_names == "gene_ids" else gene_symbols
    selected_index = science.pd.Index(selected_names.astype(str).to_numpy())
    duplicate_var_names = int(len(selected_index) - selected_index.nunique())
    if duplicate_var_names and make_unique:
        selected_index = science.ad.utils.make_index_unique(selected_index)

    obs = science.pd.DataFrame(index=barcode_values.to_numpy())
    var = science.pd.DataFrame(index=selected_index)
    var["gene_ids"] = gene_ids.astype(str).to_numpy()
    var["gene_symbols"] = gene_symbols.astype(str).to_numpy()
    if feature_types is not None:
        var["feature_types"] = feature_types.astype(str).to_numpy()

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
        adata = science.ad.AnnData(X=matrix.transpose().tocsr(), obs=obs, var=var)
    axis_audit = _validate_axis_names(
        adata,
        source_label="10x MTX input",
        make_var_names_unique=False,
    )
    payload_advisories = _validate_10x_payload(adata, source_label="10x MTX input")
    axis_audit.update(
        {
            "var_names_unique_before_repair": duplicate_var_names == 0,
            "duplicate_var_name_occurrences": duplicate_var_names,
            "var_names_repaired": duplicate_var_names if make_unique else 0,
            "var_names_unique_after_repair": bool(adata.var_names.is_unique),
        }
    )
    details = {
        "reader_parameters": {
            "var_names": var_names,
            "make_unique": make_unique,
            "gex_only": gex_only,
        },
        "input_features": original_feature_count,
        "retained_features": int(adata.n_vars),
        "cells": int(adata.n_obs),
        "feature_types_before_filter": feature_type_counts,
        "gex_filter_verified": gex_filter_verified,
        "axis_names": axis_audit,
        "matrix": _matrix_provenance(adata.X),
        "warnings": [
            *(
                ["The legacy features file has no feature-type column; gex_only could not be verified."]
                if gex_only and not gex_filter_verified
                else []
            ),
            *(
                [f"The 10x barcodes file contains {blank_barcodes} blank barcode(s)."]
                if blank_barcodes
                else []
            ),
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
                [
                    f"The retained 10x features contain {duplicate_gene_ids} duplicate stable feature ID "
                    "occurrence(s)."
                ]
                if duplicate_gene_ids
                else []
            ),
            *(
                [f"The retained 10x features contain {blank_symbols} blank feature name(s)."]
                if blank_symbols
                else []
            ),
            *_axis_advisories(axis_audit, source_label="10x MTX input"),
            *payload_advisories,
        ],
    }
    return adata, details


def _discover_10x_study(relative_root: str) -> list[tuple[str, str]]:
    root = resolve_input_path(relative_root, kind="directory")
    input_root = os.path.realpath(folder_paths.get_input_directory())
    discovered = []
    with os.scandir(root) as entries:
        child_directories = sorted(
            (
                (entry.name, entry.path, entry.is_symlink())
                for entry in entries
                if entry.is_dir(follow_symlinks=False) or entry.is_symlink()
            ),
            key=lambda item: item[0].casefold(),
        )
    if not child_directories:
        raise ValueError("The study directory does not contain any first-level Sample directories.")

    errors = []
    for sample_name, child_path, is_symlink in child_directories:
        if is_symlink:
            errors.append(f"{sample_name}: symbolic-link Sample directories are not allowed")
            continue
        child = os.path.realpath(child_path)
        if not folder_paths.is_within_directory(root, child):
            errors.append(f"{sample_name}: Sample directory escapes the selected Study directory")
            continue
        relative_child = os.path.relpath(child, input_root).replace(os.sep, "/")
        try:
            resolve_10x_mtx_files(relative_child)
        except (ValueError, FileNotFoundError, OSError) as exc:
            errors.append(f"{sample_name}: {exc}")
            continue
        discovered.append((sample_name, relative_child))
    if errors:
        details = "\n".join(f"- {message}" for message in errors)
        raise ValueError(f"Invalid 10x Study Sample inventory:\n{details}")
    return discovered


def _portable_study_path(relative_root: str) -> str:
    root = resolve_input_path(relative_root, kind="directory")
    input_root = os.path.realpath(folder_paths.get_input_directory())
    return os.path.relpath(root, input_root).replace(os.sep, "/")


def _validate_study_options(sample_key: str, var_names: str, join: str, make_unique: bool, gex_only: bool) -> str:
    if not isinstance(sample_key, str) or not sample_key.strip():
        raise ValueError("Study sample_key cannot be empty.")
    if join not in {"inner", "outer"}:
        raise ValueError("Study join must be 'inner' or 'outer'.")
    _validate_10x_options(var_names, make_unique, gex_only)
    return sample_key.strip()


def _study_feature_metadata(samples: Mapping[str, AnnData]) -> dict[str, dict[str, str]]:
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
                conflicts.append(
                    f"{var_name!r}: {previous[0]}={previous[1]!r}, {sample_name}={gene_id!r}"
                )
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


class OpenBioSingleCellLoadH5AD(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoadH5AD",
            display_name="Load H5AD",
            category=CATEGORY,
            description="Load an AnnData H5AD file from the ComfyUI input directory.",
            inputs=[
                io.String.Input(
                    "path",
                    display_name="H5AD file",
                    default="openbio-singlecell/openbio_singlecell_demo.h5ad",
                    extra_dict=_file_upload_widget(
                        ".h5ad",
                        label="Choose H5AD file",
                        drop_label="Drop an .h5ad file here",
                    ),
                ),
                io.Boolean.Input("make_var_names_unique", default=True, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, path: str, make_var_names_unique: bool = True) -> bool | str:
        try:
            resolve_input_path(path, extensions=(".h5ad",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, path: str, make_var_names_unique: bool = True) -> Any:
        return input_file_fingerprint(path, (".h5ad",))

    @classmethod
    def execute(cls, path: str, make_var_names_unique: bool = True) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        resolved = resolve_input_path(path, extensions=(".h5ad",))
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
            adata = science.ad.read_h5ad(resolved, backed=None)
        axis_audit = _validate_axis_names(
            adata,
            source_label="H5AD input",
            make_var_names_unique=make_var_names_unique,
        )
        existing_metadata = adata.uns.get("openbio_singlecell")
        embedded_source = existing_metadata.get("source", {}) if isinstance(existing_metadata, Mapping) else {}
        metadata_advisory = None
        try:
            ensure_metadata(adata)
        except ValueError as exc:
            quarantined_metadata = _sanitize_embedded_provenance(existing_metadata)
            del adata.uns["openbio_singlecell"]
            ensure_metadata(adata)
            metadata_advisory = (
                "Embedded OpenBio metadata was incompatible with the current metadata schema and was quarantined "
                f"instead of used as workflow evidence: {exc}"
            )
        else:
            quarantined_metadata = None
        source = input_file_provenance(path, (".h5ad",))
        source_warnings = _axis_advisories(axis_audit, source_label="H5AD input")
        if metadata_advisory:
            source_warnings.append(metadata_advisory)
        source.update(
            {
                "reader_parameters": {
                    "backed": None,
                    "make_var_names_unique": make_var_names_unique,
                },
                "shape": [int(adata.n_obs), int(adata.n_vars)],
                "axis_names": axis_audit,
                "expression_state": "not_inferred",
                "warnings": source_warnings,
            }
        )
        if embedded_source:
            source["embedded_source"] = _sanitize_embedded_provenance(embedded_source)
        if quarantined_metadata is not None:
            source["quarantined_embedded_openbio_metadata"] = quarantined_metadata
        ensure_metadata(
            adata,
            display_name=os.path.splitext(os.path.basename(resolved))[0],
            source=_source_metadata("h5ad", source),
        )
        return io.NodeOutput(adata)


class OpenBioSingleCellLoad10xMTX(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xMTX",
            display_name="Load 10x MTX",
            category=CATEGORY,
            description="Load a 10x matrix, barcodes, and features directory under ComfyUI input.",
            inputs=[
                io.String.Input("directory", default="openbio-singlecell/10x"),
                io.Combo.Input("var_names", options=["gene_symbols", "gene_ids"], default="gene_symbols"),
                io.Boolean.Input("make_unique", default=True, advanced=True),
                io.Boolean.Input("gex_only", default=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, directory: str, **kwargs) -> bool | str:
        try:
            _validate_10x_options(
                kwargs.get("var_names", "gene_symbols"),
                kwargs.get("make_unique", True),
                kwargs.get("gex_only", True),
            )
            resolve_10x_mtx_files(directory)
        except (TypeError, ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, directory: str, **kwargs) -> Any:
        return tenx_mtx_fingerprint(directory)

    @classmethod
    def execute(
        cls,
        directory: str,
        var_names: str = "gene_symbols",
        make_unique: bool = True,
        gex_only: bool = True,
    ) -> io.NodeOutput:
        adata, details = _read_10x_mtx(directory, var_names, make_unique, gex_only)
        source = tenx_mtx_provenance(directory)
        source.update(details)
        ensure_metadata(
            adata,
            display_name=os.path.basename(os.path.normpath(directory)) or "10x",
            source=_source_metadata("10x_mtx", source),
        )
        return io.NodeOutput(adata)


class OpenBioSingleCellLoad10xStudy(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xStudy",
            display_name="Load 10x Study",
            category=CATEGORY,
            description=(
                "Load every first-level 10x MTX Sample directory and concatenate them without silently "
                "omitting invalid biological replicates."
            ),
            inputs=[
                io.String.Input("directory", default="openbio-singlecell/study"),
                io.String.Input("sample_key", default="sample"),
                io.Combo.Input("var_names", options=["gene_symbols", "gene_ids"], default="gene_symbols"),
                io.Combo.Input("join", options=["inner", "outer"], default="inner"),
                io.Boolean.Input("make_unique", default=True, advanced=True),
                io.Boolean.Input("gex_only", default=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, directory: str, **kwargs) -> bool | str:
        try:
            _validate_study_options(
                kwargs.get("sample_key", "sample"),
                kwargs.get("var_names", "gene_symbols"),
                kwargs.get("join", "inner"),
                kwargs.get("make_unique", True),
                kwargs.get("gex_only", True),
            )
            _discover_10x_study(directory)
        except (TypeError, ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, directory: str, **kwargs) -> Any:
        _validate_study_options(
            kwargs.get("sample_key", "sample"),
            kwargs.get("var_names", "gene_symbols"),
            kwargs.get("join", "inner"),
            kwargs.get("make_unique", True),
            kwargs.get("gex_only", True),
        )
        return tuple(
            (sample_name, tenx_mtx_fingerprint(sample_directory))
            for sample_name, sample_directory in _discover_10x_study(directory)
        )

    @classmethod
    def execute(
        cls,
        directory: str,
        sample_key: str = "sample",
        var_names: str = "gene_symbols",
        join: str = "inner",
        make_unique: bool = True,
        gex_only: bool = True,
    ) -> io.NodeOutput:
        sample_key = _validate_study_options(sample_key, var_names, join, make_unique, gex_only)

        science = dependencies.require_scientific_dependencies()
        discovered = _discover_10x_study(directory)
        loaded = {}
        load_errors = []
        for sample_name, sample_directory in discovered:
            try:
                loaded[sample_name] = _read_10x_mtx(sample_directory, var_names, make_unique, gex_only)
            except (OSError, ValueError, TypeError, EOFError) as exc:
                load_errors.append(f"{sample_name}: {exc}")
        if load_errors:
            details = "\n".join(f"- {message}" for message in load_errors)
            raise ValueError(f"Invalid 10x Study Sample payloads:\n{details}")
        samples = {sample_name: value[0] for sample_name, value in loaded.items()}
        feature_metadata = _study_feature_metadata(samples)
        feature_sets = {sample_name: set(sample.var_names.astype(str)) for sample_name, sample in samples.items()}
        common_features = set.intersection(*feature_sets.values())
        union_features = set.union(*feature_sets.values())
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
                "Study concatenation retained duplicate observation names; later name-based cell alignment may be "
                "ambiguous."
            )
        if int(adata.n_obs) == 0:
            study_warnings.append("The concatenated Study has an empty observation axis.")
        study_name = os.path.basename(os.path.normpath(directory)) or "10x study"
        ensure_metadata(
            adata,
            display_name=study_name,
            source=_source_metadata(
                "10x_study",
                {
                    "path": _portable_study_path(directory),
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
                    "samples": {
                        sample_name: {
                            **tenx_mtx_provenance(sample_directory),
                            **loaded[sample_name][1],
                        }
                        for sample_name, sample_directory in discovered
                    },
                },
            ),
        )
        return io.NodeOutput(adata)


class OpenBioSingleCellLoad10xH5(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xH5",
            display_name="Load 10x H5",
            category=CATEGORY,
            description="Load a 10x Genomics HDF5 matrix from the ComfyUI input directory.",
            inputs=[
                io.String.Input(
                    "path",
                    display_name="10x H5 file",
                    default="openbio-singlecell/filtered_feature_bc_matrix.h5",
                    extra_dict=_file_upload_widget(
                        ".h5",
                        ".hdf5",
                        label="Choose 10x H5 file",
                        drop_label="Drop a 10x .h5 or .hdf5 file here",
                    ),
                ),
                io.String.Input("genome", default="", advanced=True),
                io.Boolean.Input("gex_only", default=True),
                io.Boolean.Input("make_unique", default=True, advanced=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, path: str, **kwargs) -> bool | str:
        try:
            resolve_input_path(path, extensions=(".h5", ".hdf5"))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, path: str, **kwargs) -> Any:
        return input_file_fingerprint(path, (".h5", ".hdf5"))

    @classmethod
    def execute(
        cls,
        path: str,
        genome: str = "",
        gex_only: bool = True,
        make_unique: bool = True,
    ) -> io.NodeOutput:
        if not isinstance(genome, str):
            raise TypeError("10x H5 genome must be a string.")
        if not isinstance(gex_only, bool) or not isinstance(make_unique, bool):
            raise TypeError("10x H5 gex_only and make_unique must be booleans.")
        science = dependencies.require_scientific_dependencies()
        resolved = resolve_input_path(path, extensions=(".h5", ".hdf5"))
        normalized_genome = genome.strip()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
            adata = science.sc.read_10x_h5(
                resolved,
                genome=normalized_genome or None,
                gex_only=False,
                backup_url=None,
            )
            original_feature_count = int(adata.n_vars)
            feature_type_counts_before = (
                {
                    str(key): int(value)
                    for key, value in adata.var["feature_types"].value_counts(sort=False).items()
                }
                if "feature_types" in adata.var
                else {}
            )
            if gex_only and "feature_types" in adata.var:
                adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
        axis_audit = _validate_axis_names(
            adata,
            source_label="10x H5 input",
            make_var_names_unique=make_unique,
        )
        payload_advisories = _validate_10x_payload(adata, source_label="10x H5 input")
        feature_type_counts = (
            {str(key): int(value) for key, value in adata.var["feature_types"].value_counts(sort=False).items()}
            if "feature_types" in adata.var
            else {}
        )
        source = input_file_provenance(path, (".h5", ".hdf5"))
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
            display_name=os.path.splitext(os.path.basename(resolved))[0],
            source=_source_metadata("10x_h5", source),
        )
        return io.NodeOutput(adata)


def _bounded_name(value: Any) -> str:
    text = str(value)
    if len(text) <= SLOT_NAME_CHARACTER_LIMIT:
        return text
    return f"{text[: SLOT_NAME_CHARACTER_LIMIT - 3]}..."


def _bounded_names(values: Any) -> dict[str, Any]:
    total = int(len(values))
    names = [_bounded_name(value) for value in islice(iter(values), SLOT_ITEM_LIMIT)]
    return {"names": names, "total": total, "truncated": total > SLOT_ITEM_LIMIT}


def _array_descriptor(value: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    if value is None:
        return {
            "present": False,
            "shape": None,
            "dtype": None,
            "python_type": None,
            "storage": "none",
        }
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
    all_keys = mapping.keys()
    total = sum(key is not None or not exclude_none for key in all_keys)
    keys = list(islice((key for key in all_keys if key is not None or not exclude_none), SLOT_ITEM_LIMIT))
    result = {
        "names": [_bounded_name(key) for key in keys],
        "total": total,
        "truncated": total > SLOT_ITEM_LIMIT,
    }
    descriptors = [_array_descriptor(mapping[key]) for key in keys]
    result.update(
        {
            "shapes": [item["shape"] for item in descriptors],
            "dtypes": [item["dtype"] for item in descriptors],
            "storage": [item["storage"] for item in descriptors],
        }
    )
    return result


def _empty_column_inventory() -> dict[str, Any]:
    return {"names": [], "total": 0, "truncated": False, "dtypes": []}


def _empty_mapping_inventory() -> dict[str, Any]:
    return {
        "names": [],
        "total": 0,
        "truncated": False,
        "shapes": [],
        "dtypes": [],
        "storage": [],
    }


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


def _anndata_summary_key_results(adata: AnnData) -> dict[str, Any]:
    raw = adata.raw
    metadata = adata.uns.get("openbio_singlecell")
    valid_metadata = isinstance(metadata, Mapping)
    source = metadata.get("source") if valid_metadata else None
    history = metadata.get("analysis_history") if valid_metadata else None
    raw_descriptor: dict[str, Any] = {
        "present": raw is not None,
        "shape": [int(raw.n_obs), int(raw.n_vars)] if raw is not None else None,
        "X": _array_descriptor(raw.X if raw is not None else None),
        "var_index": _index_descriptor(raw.var_names) if raw is not None else None,
    }
    return {
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "object_state": {
            "is_view": bool(adata.is_view),
            "is_backed": bool(adata.isbacked),
        },
        "X": _array_descriptor(adata.X),
        "obs_index": _index_descriptor(adata.obs_names),
        "var_index": _index_descriptor(adata.var_names),
        "raw": raw_descriptor,
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
            "warning_count": (
                len(metadata.get("warnings", []))
                if valid_metadata and isinstance(metadata.get("warnings", []), (list, tuple))
                else None
            ),
        },
    }


def _anndata_summary_code() -> str:
    return dedent(
        f'''
        from collections.abc import Mapping
        from itertools import islice

        import pandas as pd
        from scipy import sparse


        SLOT_ITEM_LIMIT = {SLOT_ITEM_LIMIT}
        SLOT_NAME_CHARACTER_LIMIT = {SLOT_NAME_CHARACTER_LIMIT}
        INDEX_EXAMPLE_LIMIT = {INDEX_EXAMPLE_LIMIT}


        def _bounded_name(value):
            text = str(value)
            if len(text) <= SLOT_NAME_CHARACTER_LIMIT:
                return text
            return f"{{text[: SLOT_NAME_CHARACTER_LIMIT - 3]}}..."


        def _bounded_names(values):
            total = int(len(values))
            names = [_bounded_name(value) for value in islice(iter(values), SLOT_ITEM_LIMIT)]
            return {{"names": names, "total": total, "truncated": total > SLOT_ITEM_LIMIT}}


        def _array_descriptor(value):
            if value is None:
                return {{
                    "present": False, "shape": None, "dtype": None,
                    "python_type": None, "storage": "none",
                }}
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
            all_keys = mapping.keys()
            total = sum(key is not None or not exclude_none for key in all_keys)
            keys = list(islice(
                (key for key in all_keys if key is not None or not exclude_none),
                SLOT_ITEM_LIMIT,
            ))
            result = {{
                "names": [_bounded_name(key) for key in keys],
                "total": total,
                "truncated": total > SLOT_ITEM_LIMIT,
            }}
            descriptors = [_array_descriptor(mapping[key]) for key in keys]
            result.update({{
                "shapes": [item["shape"] for item in descriptors],
                "dtypes": [item["dtype"] for item in descriptors],
                "storage": [item["storage"] for item in descriptors],
            }})
            return result


        def _empty_column_inventory():
            return {{"names": [], "total": 0, "truncated": False, "dtypes": []}}


        def _empty_mapping_inventory():
            return {{
                "names": [], "total": 0, "truncated": False,
                "shapes": [], "dtypes": [], "storage": [],
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


        def summarize_anndata(adata):
            raw = adata.raw
            metadata = adata.uns.get("openbio_singlecell")
            valid_metadata = isinstance(metadata, Mapping)
            source = metadata.get("source") if valid_metadata else None
            history = metadata.get("analysis_history") if valid_metadata else None
            raw_descriptor = {{
                "present": raw is not None,
                "shape": [int(raw.n_obs), int(raw.n_vars)] if raw is not None else None,
                "X": _array_descriptor(raw.X if raw is not None else None),
                "var_index": _index_descriptor(raw.var_names) if raw is not None else None,
            }}
            return {{
                "shape": [int(adata.n_obs), int(adata.n_vars)],
                "object_state": {{"is_view": bool(adata.is_view), "is_backed": bool(adata.isbacked)}},
                "X": _array_descriptor(adata.X),
                "obs_index": _index_descriptor(adata.obs_names),
                "var_index": _index_descriptor(adata.var_names),
                "raw": raw_descriptor,
                "obs_columns": _bounded_columns(adata.obs),
                "var_columns": _bounded_columns(adata.var),
                "layers": _bounded_mapping(adata.layers, exclude_none=True),
                "obsm": _bounded_mapping(adata.obsm),
                "varm": _bounded_mapping(adata.varm),
                "obsp": _bounded_mapping(adata.obsp),
                "varp": _bounded_mapping(adata.varp),
                "uns": _bounded_mapping(adata.uns),
                "raw_var_columns": (
                    _bounded_columns(raw.var) if raw is not None else _empty_column_inventory()
                ),
                "raw_varm": _bounded_mapping(raw.varm) if raw is not None else _empty_mapping_inventory(),
                "openbio_metadata": {{
                    "present": metadata is not None,
                    "valid_mapping": valid_metadata,
                    "display_name": (
                        _bounded_name(metadata.get("display_name"))
                        if valid_metadata and metadata.get("display_name") else None
                    ),
                    "source_kind": (
                        _bounded_name(source.get("kind"))
                        if isinstance(source, Mapping) and source.get("kind") else None
                    ),
                    "history_entries": len(history) if isinstance(history, Mapping) else None,
                    "warning_count": (
                        len(metadata.get("warnings", []))
                        if valid_metadata and isinstance(metadata.get("warnings", []), (list, tuple))
                        else None
                    ),
                }},
            }}
        '''
    )


class OpenBioSingleCellAnnDataSummary(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAnnDataSummary",
            display_name="AnnData Summary",
            category=DIAGNOSTIC_CATEGORY,
            description="Report a bounded structural diagnostic of an existing AnnData object without loading data.",
            inputs=[AnnDataType.Input("adata")],
            outputs=[
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(cls, adata: AnnData) -> io.NodeOutput:
        start = time.perf_counter()
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
        if isinstance(metadata, Mapping) and isinstance(metadata.get("warnings"), (list, tuple)):
            metadata_warnings = metadata["warnings"]
            warnings_list.extend(
                _bounded_name(value) for value in metadata_warnings[:METADATA_WARNING_LIMIT]
            )
            if len(metadata_warnings) > METADATA_WARNING_LIMIT:
                warnings_list.append(
                    f"OpenBio metadata warnings were limited to the first {METADATA_WARNING_LIMIT} items."
                )
        code = _anndata_summary_code()
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellAnnDataSummary",
            title=f"{display_name} summary",
            operation="anndata_summary",
            methods=(
                "Inspected AnnData axis identity, storage state, raw snapshot, and bounded slot inventories "
                "without materializing expression matrices or modifying the input object."
            ),
            results=(
                f"The object contains {cells:,} observations and {genes:,} variables, "
                f"with {sum(key is not None for key in adata.layers.keys()):,} expression layer(s), "
                f"{len(adata.obsm):,} observation "
                f"embedding(s), and {'an' if adata.raw is not None else 'no'} attached raw snapshot."
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
            started_at=start,
            code=code,
        )
        return io.NodeOutput(report, code)


INPUT_NODE_CLASSES = [
    OpenBioSingleCellLoadH5AD,
    OpenBioSingleCellLoad10xMTX,
    OpenBioSingleCellLoad10xStudy,
    OpenBioSingleCellLoad10xH5,
    OpenBioSingleCellAnnDataSummary,
]
