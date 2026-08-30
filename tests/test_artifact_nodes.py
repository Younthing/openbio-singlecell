from __future__ import annotations

import asyncio
import inspect
import sys
from dataclasses import replace

import pytest
from comfy_api.latest import io

import openbio_singlecell.extension as extension_module
from openbio_singlecell.extension import NODE_CLASSES
from openbio_singlecell.node_adapter import adapt_scientific_node
from openbio_singlecell.node_types import WorkerType
from openbio_singlecell.nodes_worker import OpenBioSingleCellPythonWorker, WorkerProfile

MAIN_PROCESS_NODE_IDS = {
    "OpenBioSingleCellCoreStudyParameters",
    "OpenBioSingleCellPythonWorker",
    "OpenBioSingleCellPreviewResult",
    "OpenBioSingleCellSaveH5AD",
    "OpenBioSingleCellExportCSV",
    "OpenBioSingleCellSavePNG",
    "OpenBioSingleCellPersistArtifact",
}


class _LazyProbeNode(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLazyProbe",
            display_name="Lazy Probe",
            category="openbio/tests",
            inputs=[io.String.Input("value", lazy=True)],
            outputs=[io.String.Output("value")],
        )

    @classmethod
    def check_lazy_status(cls, value):
        return ["value"] if value is None else []


def test_python_worker_is_registered_as_a_typed_source():
    schemas = {node.GET_SCHEMA().node_id: node.GET_SCHEMA() for node in NODE_CLASSES}

    schema = schemas["OpenBioSingleCellPythonWorker"]
    assert [item.id for item in schema.inputs] == ["python"]
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [("worker", WorkerType.io_type)]


def test_persist_artifact_is_registered_as_a_terminal_for_file_tickets():
    schemas = {node.GET_SCHEMA().node_id: node.GET_SCHEMA() for node in NODE_CLASSES}

    schema = schemas["OpenBioSingleCellPersistArtifact"]
    assert schema.display_name == "Persist Artifact"
    assert [item.id for item in schema.inputs] == ["artifact", "name"]
    assert schema.is_output_node is True
    assert schema.outputs == []


def test_python_worker_probes_the_exact_executable():
    output = asyncio.run(OpenBioSingleCellPythonWorker.execute(sys.executable))

    profile = output[0]
    assert isinstance(profile, WorkerProfile)
    assert profile.executable == OpenBioSingleCellPythonWorker.fingerprint_inputs(sys.executable)[0]
    assert profile.identity.python_version[:2] >= (3, 12)


def test_scientific_cache_fingerprint_changes_with_selected_python_path():
    profile = asyncio.run(OpenBioSingleCellPythonWorker.execute(sys.executable))[0]
    other = WorkerProfile(replace(profile.identity, executable=profile.executable + ".other"))
    log1p = next(node for node in NODE_CLASSES if node.GET_SCHEMA().node_id == "OpenBioSingleCellLog1p")

    assert log1p.fingerprint_inputs(worker=profile) != log1p.fingerprint_inputs(worker=other)


def test_every_scientific_node_is_async_and_has_optional_worker_socket():
    scientific_nodes = [
        node for node in NODE_CLASSES if node.GET_SCHEMA().node_id not in MAIN_PROCESS_NODE_IDS
    ]

    assert scientific_nodes
    for node in scientific_nodes:
        schema = node.GET_SCHEMA()
        worker_inputs = [item for item in schema.inputs if item.id == "worker"]
        assert len(worker_inputs) == 1, schema.node_id
        assert worker_inputs[0].get_io_type() == WorkerType.io_type
        assert worker_inputs[0].optional is True
        assert inspect.iscoroutinefunction(node.execute), schema.node_id


def test_scientific_adapter_forwards_lazy_status_without_the_worker_socket():
    adapted = adapt_scientific_node(_LazyProbeNode)

    assert asyncio.run(adapted.check_lazy_status(value=None, worker=None)) == ["value"]


def test_removed_scenic_binary_contract_is_not_registered():
    for node in NODE_CLASSES:
        schema = node.GET_SCHEMA()
        assert all("OPENBIO_SCENIC_BINARY" not in item.get_io_type() for item in schema.inputs)
        assert all(output.io_type != "OPENBIO_SCENIC_BINARY" for output in schema.outputs)


def test_extension_rejects_ram_pressure_cache_before_initializing_workers(monkeypatch):
    monkeypatch.setattr(extension_module.comfy_args, "cache_classic", False)
    monkeypatch.setattr(extension_module.comfy_args, "cache_none", False)
    monkeypatch.setattr(extension_module.comfy_args, "cache_lru", 0)
    monkeypatch.setattr(extension_module.comfy_args, "cache_ram", [])

    with pytest.raises(RuntimeError, match=r"RAM-pressure.*--cache-classic"):
        asyncio.run(extension_module.OpenBioSingleCellExtension().on_load())
