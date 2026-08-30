from __future__ import annotations

import asyncio
import gc
import sys
import textwrap
from pathlib import Path

import execution
import nodes
import numpy as np
import pytest
from comfy_api.latest import io

import openbio_singlecell.artifact_service as artifact_service_module
from openbio_singlecell.artifact_runtime import ArtifactTicket
from openbio_singlecell.artifact_service import (
    current_artifact_runtime,
    execute_artifact_node,
    initialize_artifact_service,
)
from openbio_singlecell.contracts import SummaryResult
from openbio_singlecell.node_adapter import adapt_scientific_node
from openbio_singlecell.node_types import AnnDataType, CNMFRunType, SCVIModelType, TableResultType
from openbio_singlecell.nodes_input import OpenBioSingleCellLoadH5AD
from openbio_singlecell.nodes_preprocess import OpenBioSingleCellLog1p
from openbio_singlecell.nodes_worker import WorkerProfile
from openbio_singlecell.worker_client import WorkerIdentity, _invoke_worker
from openbio_singlecell.worker_protocol import WorkerRequest


def _run(coroutine):
    return asyncio.run(coroutine)


class _ArgumentAdaptingNode(io.ComfyNode):
    @classmethod
    def define_schema(cls):
        return io.Schema(
            node_id="OpenBioSingleCellArgumentAdapterProbe",
            display_name="Argument Adapter Probe",
            category="openbio/tests",
            inputs=[],
            outputs=[io.String.Output("code")],
        )

    @classmethod
    def prepare_worker_arguments(cls, kwargs):
        nested = kwargs["nested"]
        return {"flattened": nested["value"]}


class _ClassicCacheProbeNode(io.ComfyNode):
    @classmethod
    def define_schema(cls):
        return io.Schema(
            node_id="OpenBioSingleCellClassicCacheProbe",
            display_name="Classic Cache Probe",
            category="openbio/tests",
            inputs=[io.String.Input("value", default="first")],
            outputs=[AnnDataType.Output("adata")],
        )


class _InputReferenceProbeNode(io.ComfyNode):
    @classmethod
    def define_schema(cls):
        return io.Schema(
            node_id="OpenBioSingleCellInputReferenceProbe",
            display_name="Input Reference Probe",
            category="openbio/tests",
            inputs=[AnnDataType.Input("adata")],
            outputs=[TableResultType.Output("table")],
        )


class _NativeProducerProbeNode(io.ComfyNode):
    node_id = ""
    output_type = AnnDataType

    @classmethod
    def define_schema(cls):
        return io.Schema(
            node_id=cls.node_id,
            display_name="Native Producer Probe",
            category="openbio/tests",
            inputs=[],
            outputs=[cls.output_type.Output("native")],
        )


class _SCVINativeProducerProbeNode(_NativeProducerProbeNode):
    node_id = "OpenBioSingleCellSCVINativeProducerProbe"
    output_type = SCVIModelType


class _CNMFNativeProducerProbeNode(_NativeProducerProbeNode):
    node_id = "OpenBioSingleCellCNMFNativeProducerProbe"
    output_type = CNMFRunType


class _NativeConsumerProbeNode(io.ComfyNode):
    @classmethod
    def define_schema(cls):
        return io.Schema(
            node_id="OpenBioSingleCellNativeConsumerProbe",
            display_name="Native Consumer Probe",
            category="openbio/tests",
            inputs=[],
            outputs=[io.String.Output("code")],
        )


class _Server:
    client_id = None

    @staticmethod
    def send_sync(*_args, **_kwargs):
        return None


def _worker_identity(marker: int) -> WorkerIdentity:
    executable = str(Path(sys.executable).resolve())
    stat = Path(executable).stat()
    return WorkerIdentity(
        executable=executable,
        size=stat.st_size,
        mtime_ns=stat.st_mtime_ns + marker,
        protocol_version=1,
        python_version=(3, 12, 0),
        dependencies=(
            ("anndata", "test"),
            ("numpy", "test"),
            ("pandas", "test"),
            ("scanpy", "test"),
            ("scipy", "test"),
        ),
    )


def test_artifact_service_applies_node_worker_argument_adapter(
    comfy_directories,
    monkeypatch,
):
    _, _, temp_dir = comfy_directories
    observed = {}

    async def fake_run_worker(identity, request_path, operation, *, inputs, parameters):
        observed.update(operation=operation, inputs=inputs, parameters=parameters)
        return [{"type": "string", "name": "code", "value": "ok"}]

    monkeypatch.setattr(artifact_service_module, "run_worker", fake_run_worker)
    _run(initialize_artifact_service(temp_dir, sys.executable))

    output = _run(
        execute_artifact_node(
            _ArgumentAdaptingNode,
            None,
            {"nested": {"value": "prepared"}},
        )
    )

    assert output[0] == "ok"
    assert observed == {
        "operation": "openbio.node.argumentadapterprobe",
        "inputs": {},
        "parameters": {"flattened": "prepared"},
    }


def test_artifact_service_rejects_non_string_worker_string_output(
    comfy_directories,
    monkeypatch,
):
    _, _, temp_dir = comfy_directories

    async def fake_run_worker(_identity, _request_path, _operation, *, inputs, parameters):
        assert inputs == {}
        assert parameters == {"flattened": "prepared"}
        return [{"type": "string", "name": "code", "value": {"not": "a string"}}]

    monkeypatch.setattr(artifact_service_module, "run_worker", fake_run_worker)
    _run(initialize_artifact_service(temp_dir, sys.executable))

    with pytest.raises(ValueError, match="string.*must contain a string"):
        _run(
            execute_artifact_node(
                _ArgumentAdaptingNode,
                None,
                {"nested": {"value": "prepared"}},
            )
        )

    assert not list(current_artifact_runtime().session_path.glob("*.partial"))


@pytest.mark.parametrize(("field", "invalid"), [("codec", []), ("payload", [])])
def test_artifact_service_rejects_non_string_worker_artifact_fields(
    comfy_directories,
    monkeypatch,
    field,
    invalid,
):
    _, _, temp_dir = comfy_directories

    async def fake_run_worker(_identity, request_path, _operation, *, inputs, parameters):
        assert inputs == {}
        assert parameters == {"value": "first"}
        payload = Path(request_path).parent / "payload"
        payload.mkdir()
        (payload / "data.bin").write_bytes(b"data")
        record = {
            "type": "artifact",
            "name": "adata",
            "kind": "OPENBIO_ANNDATA",
            "codec": "test-binary",
            "payload": "payload",
        }
        record[field] = invalid
        return [record]

    monkeypatch.setattr(artifact_service_module, "run_worker", fake_run_worker)
    gc.collect()
    _run(initialize_artifact_service(temp_dir, sys.executable))

    with pytest.raises(ValueError, match=rf"artifact.*{field}.*non-empty string"):
        _run(execute_artifact_node(_ClassicCacheProbeNode, None, {"value": "first"}))

    assert not list(current_artifact_runtime().session_path.glob("*.partial"))


def test_artifact_service_rejects_input_reference_with_the_wrong_wire_type(
    comfy_directories,
    monkeypatch,
):
    _, _, temp_dir = comfy_directories

    async def fake_run_worker(_identity, request_path, operation, *, inputs, parameters):
        if operation == "openbio.node.inputreferenceprobe":
            assert set(inputs) == {"adata"}
            assert parameters == {}
            return [{"type": "input_ref", "name": "table", "input": "adata"}]
        assert inputs == {}
        assert parameters == {"value": "first"}
        payload = Path(request_path).parent / "payload"
        payload.mkdir()
        (payload / "data.bin").write_bytes(b"data")
        return [
            {
                "type": "artifact",
                "name": "adata",
                "kind": "OPENBIO_ANNDATA",
                "codec": "test-binary",
                "payload": "payload",
            }
        ]

    monkeypatch.setattr(artifact_service_module, "run_worker", fake_run_worker)
    gc.collect()
    _run(initialize_artifact_service(temp_dir, sys.executable))
    adata = _run(execute_artifact_node(_ClassicCacheProbeNode, None, {"value": "first"}))[0]

    with pytest.raises(ValueError, match="input reference.*wire type"):
        _run(execute_artifact_node(_InputReferenceProbeNode, None, {"adata": adata}))

    del adata
    gc.collect()


@pytest.mark.parametrize(
    ("producer_node", "kind", "codec"),
    [
        (_SCVINativeProducerProbeNode, "OPENBIO_SCVI_MODEL", "scvi-native-directory"),
        (_CNMFNativeProducerProbeNode, "OPENBIO_CNMF_RUN", "cnmf-native-directory"),
    ],
)
def test_native_artifact_uses_its_producer_worker_and_rejects_a_different_worker(
    comfy_directories,
    monkeypatch,
    producer_node,
    kind,
    codec,
):
    gc.collect()
    _, _, temp_dir = comfy_directories
    default_identity = _worker_identity(0)
    producer_identity = _worker_identity(1)
    wrong_identity = _worker_identity(2)
    calls = []

    async def fake_probe(_python):
        return default_identity

    async def fake_run_worker(identity, request_path, operation, *, inputs, parameters):
        calls.append((operation, identity, inputs, parameters))
        if operation == "openbio.node.nativeconsumerprobe":
            return [{"type": "string", "name": "code", "value": "ok"}]
        payload = Path(request_path).parent / "native"
        payload.mkdir()
        (payload / "weights.bin").write_bytes(b"native")
        return [
            {
                "type": "artifact",
                "name": "native",
                "kind": kind,
                "codec": codec,
                "payload": "native",
            }
        ]

    monkeypatch.setattr(artifact_service_module, "probe_python", fake_probe)
    monkeypatch.setattr(artifact_service_module, "run_worker", fake_run_worker)
    _run(initialize_artifact_service(temp_dir, sys.executable))

    ticket = _run(execute_artifact_node(producer_node, WorkerProfile(producer_identity), {}))[0]
    consumed = _run(execute_artifact_node(_NativeConsumerProbeNode, None, {"native": ticket}))

    assert consumed[0] == "ok"
    assert [call[1] for call in calls] == [producer_identity, producer_identity]
    with pytest.raises(ValueError, match="does not match the native artifact producer"):
        _run(
            execute_artifact_node(
                _NativeConsumerProbeNode,
                WorkerProfile(wrong_identity),
                {"native": ticket},
            )
        )
    assert len(calls) == 2
    del consumed, ticket
    gc.collect()


def test_comfy_classic_cache_hit_skips_worker_and_eviction_releases_artifact(
    comfy_directories,
    monkeypatch,
):
    _, _, temp_dir = comfy_directories
    calls = 0

    async def fake_run_worker(_identity, request_path, _operation, *, inputs, parameters):
        nonlocal calls
        assert inputs == {}
        calls += 1
        payload = Path(request_path).parent / "payload"
        payload.mkdir()
        (payload / "data.bin").write_bytes(str(parameters["value"]).encode())
        return [
            {
                "type": "artifact",
                "name": "adata",
                "kind": "OPENBIO_ANNDATA",
                "codec": "test-binary",
                "payload": "payload",
            }
        ]

    monkeypatch.setattr(artifact_service_module, "run_worker", fake_run_worker)
    monkeypatch.setitem(
        nodes.NODE_CLASS_MAPPINGS,
        "OpenBioClassicCacheProbe",
        adapt_scientific_node(_ClassicCacheProbeNode),
    )

    async def exercise_classic_cache():
        gc.collect()
        await initialize_artifact_service(temp_dir, sys.executable)
        executor = execution.PromptExecutor(
            _Server(),
            cache_type=execution.CacheType.CLASSIC,
            cache_args={"ram": 0, "ram_inactive": 0},
        )
        first_prompt = {
            "1": {"class_type": "OpenBioClassicCacheProbe", "inputs": {"value": "first"}}
        }
        await executor.execute_async(first_prompt, "first", execute_outputs=["1"])
        runtime = current_artifact_runtime()
        first_ready = next(runtime.session_path.glob("*.ready"))

        await executor.execute_async(first_prompt, "hit", execute_outputs=["1"])
        cached = next(data for event, data in executor.status_messages if event == "execution_cached")
        assert cached["nodes"] == ["1"]
        assert calls == 1
        assert first_ready.is_dir()

        changed_prompt = {
            "1": {"class_type": "OpenBioClassicCacheProbe", "inputs": {"value": "changed"}}
        }
        await executor.execute_async(changed_prompt, "changed", execute_outputs=["1"])
        gc.collect()
        assert calls == 2
        assert not first_ready.exists()
        executor.reset()
        gc.collect()

    _run(exercise_classic_cache())


def test_cancelling_artifact_execution_removes_its_partial_run(
    comfy_directories,
    monkeypatch,
):
    _, _, temp_dir = comfy_directories
    started = asyncio.Event()
    partial_path = None

    async def blocking_worker(_identity, request_path, _operation, *, inputs, parameters):
        nonlocal partial_path
        assert inputs == {}
        assert parameters == {"flattened": "ignored"}
        partial_path = Path(request_path).parent
        (partial_path / "unfinished.bin").write_bytes(b"unfinished")
        started.set()
        await asyncio.Event().wait()
        raise AssertionError("unreachable")

    monkeypatch.setattr(artifact_service_module, "run_worker", blocking_worker)

    async def cancel_execution():
        gc.collect()
        await initialize_artifact_service(temp_dir, sys.executable)
        task = asyncio.create_task(
            execute_artifact_node(
                _ArgumentAdaptingNode,
                None,
                {"nested": {"value": "ignored"}},
            )
        )
        await started.wait()
        task.cancel()
        with pytest.raises(asyncio.CancelledError):
            await task

    _run(cancel_execution())
    assert partial_path is not None
    assert not partial_path.exists()
    assert not list(current_artifact_runtime().session_path.glob("*.partial"))
    assert current_artifact_runtime().close()


def test_comfy_ui_interrupt_cancels_worker_tree_and_removes_partial_run(
    comfy_directories,
    monkeypatch,
    tmp_path,
):
    _, _, temp_dir = comfy_directories
    started_marker = tmp_path / "worker-started.txt"
    descendant_marker = tmp_path / "descendant-finished.txt"
    worker_entry = tmp_path / "interrupt_worker.py"
    worker_entry.write_text(
        textwrap.dedent(
            """
            import argparse
            import json
            import subprocess
            import sys
            import time
            from pathlib import Path

            parser = argparse.ArgumentParser()
            parser.add_argument("--request", required=True)
            request_path = Path(parser.parse_args().request)
            request = json.loads(request_path.read_text(encoding="utf-8"))
            parameters = request["parameters"]

            time.sleep(0.25)
            child_code = (
                "import sys, time; from pathlib import Path; "
                "time.sleep(1.0); Path(sys.argv[1]).write_text('survived', encoding='utf-8')"
            )
            subprocess.Popen([sys.executable, "-I", "-c", child_code, parameters["descendant_marker"]])
            Path(parameters["started_marker"]).write_text("started", encoding="utf-8")
            time.sleep(2.0)
            response = {
                "version": request["version"],
                "request_id": request["request_id"],
                "ok": True,
                "outputs": [],
            }
            request_path.with_name("response.json").write_text(
                json.dumps(response), encoding="utf-8"
            )
            """
        ),
        encoding="utf-8",
    )

    async def process_worker(_identity, request_path, _operation, *, inputs, parameters):
        assert inputs == {}
        assert parameters == {"value": "first"}
        partial_path = Path(request_path).parent
        (partial_path / "unfinished.bin").write_bytes(b"unfinished")
        request = WorkerRequest.create(
            "openbio.test.interrupt",
            inputs={},
            parameters={
                "started_marker": str(started_marker),
                "descendant_marker": str(descendant_marker),
            },
        )
        await _invoke_worker(sys.executable, worker_entry, request_path, request)
        payload = partial_path / "payload"
        payload.mkdir()
        (payload / "data.bin").write_bytes(b"completed")
        return [
            {
                "type": "artifact",
                "name": "adata",
                "kind": "OPENBIO_ANNDATA",
                "codec": "test-binary",
                "payload": "payload",
            }
        ]

    monkeypatch.setattr(artifact_service_module, "run_worker", process_worker)
    monkeypatch.setitem(
        nodes.NODE_CLASS_MAPPINGS,
        "OpenBioInterruptProbe",
        adapt_scientific_node(_ClassicCacheProbeNode),
    )

    async def interrupt_prompt():
        gc.collect()
        await initialize_artifact_service(temp_dir, sys.executable)
        executor = execution.PromptExecutor(
            _Server(),
            cache_type=execution.CacheType.CLASSIC,
            cache_args={"ram": 0, "ram_inactive": 0},
        )
        prompt = {"1": {"class_type": "OpenBioInterruptProbe", "inputs": {"value": "first"}}}
        execution_task = asyncio.create_task(
            executor.execute_async(prompt, "interrupt", execute_outputs=["1"])
        )
        async with asyncio.timeout(5):
            while not started_marker.exists():
                await asyncio.sleep(0.02)
        interrupted_at = asyncio.get_running_loop().time()
        nodes.interrupt_processing()
        await asyncio.wait_for(execution_task, timeout=5)
        cancellation_latency = asyncio.get_running_loop().time() - interrupted_at
        await asyncio.sleep(1.1)
        return executor, cancellation_latency

    try:
        executor, cancellation_latency = _run(interrupt_prompt())
    finally:
        nodes.interrupt_processing(False)

    assert cancellation_latency < 1.0
    assert any(event == "execution_interrupted" for event, _data in executor.status_messages)
    assert not descendant_marker.exists()
    assert not list(current_artifact_runtime().session_path.glob("*.partial"))
    assert current_artifact_runtime().close()


def test_load_then_log1p_publishes_new_file_tickets_without_changing_input(
    comfy_directories,
    science,
):
    input_dir, _, temp_dir = comfy_directories
    source = science.ad.AnnData(
        X=science.sparse.csr_matrix(np.array([[0.0, 3.0], [1.0, 8.0]])),
    )
    source.obs_names = ["cell-a", "cell-b"]
    source.var_names = ["gene-a", "gene-b"]
    source_path = input_dir / "source.h5ad"
    source.write_h5ad(source_path, compression=None, convert_strings_to_categoricals=False)
    original_bytes = source_path.read_bytes()

    _run(initialize_artifact_service(temp_dir, sys.executable))
    loaded = _run(
        execute_artifact_node(
            OpenBioSingleCellLoadH5AD,
            None,
            {"path": "source.h5ad", "make_var_names_unique": True},
        )
    )[0]
    transformed_output = _run(
        execute_artifact_node(OpenBioSingleCellLog1p, None, {"adata": loaded})
    )
    transformed, summary, code = transformed_output.args

    assert isinstance(loaded, ArtifactTicket)
    assert isinstance(transformed, ArtifactTicket)
    assert transformed != loaded
    assert isinstance(summary, SummaryResult)
    assert summary.summary["node_id"] == "OpenBioSingleCellLog1p"
    assert "sc.pp.log1p" in code
    assert source_path.read_bytes() == original_bytes

    runtime = current_artifact_runtime()
    loaded_adata = science.ad.read_h5ad(runtime.resolve(loaded) / "data.h5ad")
    transformed_adata = science.ad.read_h5ad(runtime.resolve(transformed) / "data.h5ad")
    np.testing.assert_allclose(loaded_adata.X.toarray(), [[0.0, 3.0], [1.0, 8.0]])
    np.testing.assert_allclose(
        transformed_adata.X.toarray(),
        np.log1p([[0.0, 3.0], [1.0, 8.0]]),
    )
    assert loaded_adata.obs_names.tolist() == transformed_adata.obs_names.tolist()
    assert loaded_adata.var_names.tolist() == transformed_adata.var_names.tolist()


def test_worker_failure_leaves_no_ready_artifact(comfy_directories, science):
    input_dir, _, temp_dir = comfy_directories
    source = science.ad.AnnData(X=np.array([[0.0, -1.0]]))
    source_path = Path(input_dir, "invalid.h5ad")
    source.write_h5ad(source_path, compression=None, convert_strings_to_categoricals=False)

    _run(initialize_artifact_service(temp_dir, sys.executable))
    loaded = _run(
        execute_artifact_node(
            OpenBioSingleCellLoadH5AD,
            None,
            {"path": "invalid.h5ad", "make_var_names_unique": True},
        )
    )[0]
    runtime = current_artifact_runtime()
    ready_before = set(runtime.session_path.glob("*.ready"))

    try:
        _run(execute_artifact_node(OpenBioSingleCellLog1p, None, {"adata": loaded}))
    except RuntimeError as error:
        assert "values <= -1" in str(error)
    else:
        raise AssertionError("invalid Log1p input did not fail")

    assert set(runtime.session_path.glob("*.ready")) == ready_before
    assert not list(runtime.session_path.glob("*.partial"))


def test_two_async_branches_read_one_parent_and_publish_independent_runs(
    comfy_directories,
    science,
):
    input_dir, _, temp_dir = comfy_directories
    source = science.ad.AnnData(X=np.array([[0.0, 2.0], [3.0, 4.0]]))
    source.obs_names = ["a", "b"]
    source.var_names = ["x", "y"]
    source_path = input_dir / "branches.h5ad"
    source.write_h5ad(source_path, compression=None, convert_strings_to_categoricals=False)
    original_bytes = source_path.read_bytes()

    async def execute_branches():
        gc.collect()
        await initialize_artifact_service(temp_dir, sys.executable)
        parent = (
            await execute_artifact_node(
                OpenBioSingleCellLoadH5AD,
                None,
                {"path": "branches.h5ad", "make_var_names_unique": True},
            )
        )[0]
        left, right = await asyncio.gather(
            execute_artifact_node(OpenBioSingleCellLog1p, None, {"adata": parent}),
            execute_artifact_node(OpenBioSingleCellLog1p, None, {"adata": parent}),
        )
        return parent, left[0], right[0]

    parent, left, right = _run(execute_branches())
    runtime = current_artifact_runtime()

    assert left.run != right.run
    assert left != right
    assert runtime.resolve(parent).is_dir()
    np.testing.assert_allclose(
        science.ad.read_h5ad(runtime.resolve(parent) / "data.h5ad").X,
        source.X,
    )
    np.testing.assert_allclose(
        science.ad.read_h5ad(runtime.resolve(left) / "data.h5ad").X,
        science.ad.read_h5ad(runtime.resolve(right) / "data.h5ad").X,
    )
    assert source_path.read_bytes() == original_bytes
