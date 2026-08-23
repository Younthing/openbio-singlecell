from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from typing import Literal

import folder_paths

PLUGIN_DIRECTORY = "openbio-singlecell"
OUTPUT_EXTENSIONS = {"h5ad", "csv", "png"}
TENX_FILE_CANDIDATES = {
    "matrix": ("matrix.mtx", "matrix.mtx.gz"),
    "barcodes": ("barcodes.tsv", "barcodes.tsv.gz"),
    "features": ("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"),
}


@dataclass(frozen=True, slots=True)
class OutputTarget:
    path: str
    filename: str
    subfolder: str
    folder_type: Literal["output", "temp"]


def _same_path(left: str, right: str) -> bool:
    return os.path.normcase(os.path.realpath(left)) == os.path.normcase(os.path.realpath(right))


def resolve_input_path(
    relative_path: str,
    *,
    kind: Literal["file", "directory"] = "file",
    extensions: tuple[str, ...] | None = None,
) -> str:
    if not isinstance(relative_path, str) or not relative_path.strip():
        raise ValueError("Input path cannot be empty.")

    name, annotated_root = folder_paths.annotated_filepath(relative_path.strip())
    input_root = folder_paths.get_input_directory()
    if annotated_root is not None and not _same_path(annotated_root, input_root):
        raise ValueError("Only files under the ComfyUI input directory are allowed.")
    if os.path.isabs(name) or os.path.splitdrive(name)[0]:
        raise ValueError("Input path must be relative to the ComfyUI input directory.")

    target = os.path.abspath(os.path.join(input_root, name))
    if not folder_paths.is_within_directory(input_root, target):
        raise ValueError(f"Invalid input path: {relative_path!r}")
    target = os.path.realpath(target)

    if kind == "file" and not os.path.isfile(target):
        raise FileNotFoundError(f"Input file not found: {relative_path}")
    if kind == "directory" and not os.path.isdir(target):
        raise FileNotFoundError(f"Input directory not found: {relative_path}")
    if extensions is not None and not target.lower().endswith(tuple(ext.lower() for ext in extensions)):
        allowed = ", ".join(extensions)
        raise ValueError(f"Unsupported input file type; expected one of: {allowed}")
    return target


def file_fingerprint(path: str) -> tuple[str, int, int]:
    stat = os.stat(path)
    canonical_path = os.path.normcase(os.path.realpath(path))
    return canonical_path, int(stat.st_size), int(stat.st_mtime_ns)


def input_file_fingerprint(relative_path: str, extensions: tuple[str, ...]) -> tuple[str, int, int]:
    return file_fingerprint(resolve_input_path(relative_path, extensions=extensions))


def resolve_10x_mtx_files(relative_directory: str) -> dict[str, str]:
    directory = resolve_input_path(relative_directory, kind="directory")
    selected = {}
    for role, candidates in TENX_FILE_CANDIDATES.items():
        for filename in candidates:
            candidate = os.path.join(directory, filename)
            if os.path.isfile(candidate) and folder_paths.is_within_directory(directory, candidate):
                selected[role] = os.path.realpath(candidate)
                break
        else:
            expected = ", ".join(candidates)
            raise FileNotFoundError(f"10x directory is missing {role}; expected one of: {expected}")
    return selected


def tenx_mtx_fingerprint(relative_directory: str) -> tuple[tuple[str, str, int, int], ...]:
    selected = resolve_10x_mtx_files(relative_directory)
    return tuple((role, *file_fingerprint(selected[role])) for role in ("matrix", "barcodes", "features"))


def prepare_output_target(
    filename_prefix: str,
    extension: str,
    overwrite: bool = False,
    *,
    folder_type: Literal["output", "temp"] = "output",
) -> OutputTarget:
    extension = extension.lower().lstrip(".")
    if extension not in OUTPUT_EXTENSIONS:
        raise ValueError(f"Unsupported output extension: {extension!r}")
    if not isinstance(filename_prefix, str) or not filename_prefix.strip():
        raise ValueError("Filename prefix cannot be empty.")

    prefix = filename_prefix.strip()
    suffix = f".{extension}"
    if prefix.lower().endswith(suffix):
        prefix = prefix[: -len(suffix)]

    root = folder_paths.get_directory_by_type(folder_type)
    if root is None:
        raise ValueError(f"Unknown output folder type: {folder_type}")
    plugin_root = os.path.join(root, PLUGIN_DIRECTORY)
    if not folder_paths.is_within_directory(root, plugin_root):
        raise ValueError("openbio-singlecell output directory escapes the ComfyUI output root.")

    requested_folder = os.path.abspath(os.path.join(plugin_root, os.path.dirname(os.path.normpath(prefix))))
    if not folder_paths.is_within_directory(plugin_root, requested_folder):
        raise ValueError("Output filename prefix escapes the openbio-singlecell output directory.")

    full_folder, filename, counter, _, _ = folder_paths.get_save_image_path(prefix, plugin_root)
    if filename in {"", ".", ".."}:
        raise ValueError("Filename prefix must include a file name.")
    output_filename = f"{filename}.{extension}" if overwrite else f"{filename}_{counter:05}.{extension}"
    output_path = os.path.join(full_folder, output_filename)
    if not folder_paths.is_within_directory(plugin_root, output_path):
        raise ValueError("Output path escapes the openbio-singlecell output directory.")

    subfolder = os.path.relpath(full_folder, root)
    if subfolder == ".":
        subfolder = ""
    return OutputTarget(output_path, output_filename, subfolder.replace(os.sep, "/"), folder_type)


def preview_output_target(unique_id: str | int | None) -> OutputTarget:
    root = folder_paths.get_temp_directory()
    plugin_root = os.path.join(root, PLUGIN_DIRECTORY)
    if not folder_paths.is_within_directory(root, plugin_root):
        raise ValueError("openbio-singlecell preview directory escapes the ComfyUI temp root.")
    os.makedirs(plugin_root, exist_ok=True)
    digest = hashlib.sha256(str(unique_id or "preview").encode("utf-8")).hexdigest()[:20]
    filename = f"preview_{digest}.png"
    return OutputTarget(
        path=os.path.join(plugin_root, filename),
        filename=filename,
        subfolder=PLUGIN_DIRECTORY,
        folder_type="temp",
    )
