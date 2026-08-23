from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

PLUGIN_ROOT = Path(__file__).resolve().parents[1]


def resolve_comfy_root() -> Path:
    configured = os.environ.get("OPENBIO_COMFYUI_ROOT")
    if configured:
        candidates = [Path(configured)]
        source = "OPENBIO_COMFYUI_ROOT"
    else:
        candidates = []
        if PLUGIN_ROOT.parent.name.lower() == "custom_nodes":
            candidates.append(PLUGIN_ROOT.parent.parent)
        candidates.append(PLUGIN_ROOT.parent / "ComfyUI")
        source = None

    for candidate in candidates:
        if (candidate / "main.py").is_file():
            return candidate.resolve()

    if source is not None:
        raise RuntimeError(f"{source} does not point to a ComfyUI root containing main.py: {candidates[0]}")
    raise RuntimeError("ComfyUI root was not found. Set OPENBIO_COMFYUI_ROOT before running tests.")


COMFY_ROOT = resolve_comfy_root()
sys.path.insert(0, str(PLUGIN_ROOT))
sys.path.insert(0, str(COMFY_ROOT))

import folder_paths  # noqa: E402


@pytest.fixture
def comfy_directories(tmp_path):
    original_input = folder_paths.get_input_directory()
    original_output = folder_paths.get_output_directory()
    original_temp = folder_paths.get_temp_directory()
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    temp_dir = tmp_path / "temp"
    input_dir.mkdir()
    output_dir.mkdir()
    temp_dir.mkdir()
    folder_paths.set_input_directory(str(input_dir))
    folder_paths.set_output_directory(str(output_dir))
    folder_paths.set_temp_directory(str(temp_dir))
    try:
        yield input_dir, output_dir, temp_dir
    finally:
        folder_paths.set_input_directory(original_input)
        folder_paths.set_output_directory(original_output)
        folder_paths.set_temp_directory(original_temp)
