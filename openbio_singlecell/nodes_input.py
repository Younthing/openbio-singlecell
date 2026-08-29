from __future__ import annotations

import os
from typing import Any

import folder_paths
from comfy_api.latest import io

from .files import (
    input_file_fingerprint,
    resolve_10x_mtx_files,
    resolve_input_path,
    tenx_mtx_fingerprint,
)
from .node_types import AnnDataType, SummaryResultType

CATEGORY = "openbio/single-cell/input"
DIAGNOSTIC_CATEGORY = "openbio/single-cell/diagnostics"
INPUT_FILE_UPLOAD_WIDGET = "OPENBIO_INPUT_FILE_UPLOAD_WIDGET"


def _file_upload_widget(*extensions: str, label: str, drop_label: str) -> dict[str, Any]:
    return {
        "widgetType": INPUT_FILE_UPLOAD_WIDGET,
        "allowed_extensions": list(extensions),
        "accept": ",".join(extensions),
        "upload_subfolder": "openbio-singlecell",
        "upload_label": label,
        "drop_label": drop_label,
    }


def _validate_10x_options(var_names: object, make_unique: object, gex_only: object) -> None:
    if var_names not in {"gene_symbols", "gene_ids"}:
        raise ValueError("10x var_names must be 'gene_symbols' or 'gene_ids'.")
    if not isinstance(make_unique, bool):
        raise TypeError("10x make_unique must be a boolean.")
    if not isinstance(gex_only, bool):
        raise TypeError("10x gex_only must be a boolean.")


def _discover_10x_study(relative_root: str) -> list[tuple[str, str]]:
    root = resolve_input_path(relative_root, kind="directory")
    input_root = os.path.realpath(folder_paths.get_input_directory())
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

    discovered = []
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


def _validate_study_options(
    sample_key: object,
    var_names: object,
    join: object,
    make_unique: object,
    gex_only: object,
) -> None:
    if not isinstance(sample_key, str) or not sample_key.strip():
        raise ValueError("Study sample_key cannot be empty.")
    if join not in {"inner", "outer"}:
        raise ValueError("Study join must be 'inner' or 'outer'.")
    _validate_10x_options(var_names, make_unique, gex_only)


class OpenBioSingleCellLoadH5AD(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoadH5AD",
            display_name="Load H5AD",
            category=CATEGORY,
            description="Load an AnnData H5AD file through a one-shot Python worker.",
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
    def validate_inputs(
        cls,
        path: str,
        make_var_names_unique: bool = True,
        **kwargs: Any,
    ) -> bool | str:
        try:
            resolve_input_path(path, extensions=(".h5ad",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(
        cls,
        path: str,
        make_var_names_unique: bool = True,
        **kwargs: Any,
    ) -> Any:
        return input_file_fingerprint(path, (".h5ad",))


class OpenBioSingleCellLoad10xMTX(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xMTX",
            display_name="Load 10x MTX",
            category=CATEGORY,
            description="Load a 10x matrix, barcodes, and features directory through a one-shot worker.",
            inputs=[
                io.String.Input("directory", default="openbio-singlecell/10x"),
                io.Combo.Input("var_names", options=["gene_symbols", "gene_ids"], default="gene_symbols"),
                io.Boolean.Input("make_unique", default=True, advanced=True),
                io.Boolean.Input("gex_only", default=True),
            ],
            outputs=[AnnDataType.Output(display_name="adata")],
        )

    @classmethod
    def validate_inputs(cls, directory: str, **kwargs: Any) -> bool | str:
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
    def fingerprint_inputs(cls, directory: str, **kwargs: Any) -> Any:
        return tenx_mtx_fingerprint(directory)


class OpenBioSingleCellLoad10xStudy(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xStudy",
            display_name="Load 10x Study",
            category=CATEGORY,
            description=(
                "Load every first-level 10x MTX Sample directory without silently omitting invalid replicates."
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
    def validate_inputs(cls, directory: str, **kwargs: Any) -> bool | str:
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
    def fingerprint_inputs(cls, directory: str, **kwargs: Any) -> Any:
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


class OpenBioSingleCellLoad10xH5(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLoad10xH5",
            display_name="Load 10x H5",
            category=CATEGORY,
            description="Load a 10x Genomics HDF5 matrix through a one-shot worker.",
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
    def validate_inputs(cls, path: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(path, extensions=(".h5", ".hdf5"))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, path: str, **kwargs: Any) -> Any:
        return input_file_fingerprint(path, (".h5", ".hdf5"))


class OpenBioSingleCellAnnDataSummary(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAnnDataSummary",
            display_name="AnnData Summary",
            category=DIAGNOSTIC_CATEGORY,
            description="Report a bounded structural diagnostic without returning large data to the main process.",
            inputs=[AnnDataType.Input("adata")],
            outputs=[
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )


INPUT_NODE_CLASSES = [
    OpenBioSingleCellLoadH5AD,
    OpenBioSingleCellLoad10xMTX,
    OpenBioSingleCellLoad10xStudy,
    OpenBioSingleCellLoad10xH5,
    OpenBioSingleCellAnnDataSummary,
]


__all__ = [
    "INPUT_NODE_CLASSES",
    "OpenBioSingleCellAnnDataSummary",
    "OpenBioSingleCellLoad10xH5",
    "OpenBioSingleCellLoad10xMTX",
    "OpenBioSingleCellLoad10xStudy",
    "OpenBioSingleCellLoadH5AD",
]
