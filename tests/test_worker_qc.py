from __future__ import annotations

import hashlib
import json
import uuid
from pathlib import Path

import pytest

from openbio_singlecell import operations_qc  # noqa: F401 - importing registers the QC operation allowlist.
from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_anndata, read_plot, write_anndata
from openbio_singlecell.nodes_qc import QC_NODE_CLASSES, OpenBioSingleCellQCPlots
from openbio_singlecell.worker_protocol import (
    OperationContext,
    ProtocolError,
    WorkerRequest,
    execute_request,
)

ANNDATA_KIND = "OPENBIO_ANNDATA"
ANNDATA_CODEC = "anndata-h5ad-v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifact_descriptor(root: Path) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": ANNDATA_KIND,
        "codec": ANNDATA_CODEC,
        "path": str(root.resolve()),
    }


def _context(root: Path) -> OperationContext:
    root.mkdir()
    return OperationContext.from_request_path(root / "request.json", str(uuid.uuid4()))


def _execute(
    context: OperationContext,
    operation: str,
    inputs: dict[str, object],
    parameters: dict[str, object],
):
    request = WorkerRequest.create(
        operation,
        inputs=inputs,
        parameters=parameters,
        request_id=context.request_id,
    )
    return execute_request(request, context).outputs


def test_qc_plots_schema_exposes_overview_and_grouped_views():
    schema = OpenBioSingleCellQCPlots.define_schema()

    assert [item.id for item in schema.inputs] == ["adata", "view", "source"]
    view = schema.inputs[1]
    assert [option.key for option in view.options] == ["overview", "grouped"]
    assert [item.id for item in view.options[0].inputs] == []
    assert [item.id for item in view.options[1].inputs] == ["groupby"]
    assert view.options[1].inputs[0].default == "sample"


def test_calculate_qc_publishes_new_anndata_without_mutating_input(tmp_path: Path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    matrix = science.sparse.csr_matrix([[10.0, 0.0, 2.0], [0.0, 4.0, 1.0]])
    adata = science.ad.AnnData(
        matrix,
        obs=science.pd.DataFrame(index=["cell_a", "cell_b"]),
        var=science.pd.DataFrame(index=["MT-A", "RPS1", "G1"]),
    )
    adata.layers["counts"] = matrix.copy()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = _sha256(input_path)

    records = _execute(
        _context(tmp_path / "calculate.partial"),
        "openbio.node.calculateqc",
        {"adata": _artifact_descriptor(input_root)},
        {
            "include_ribosomal": True,
            "include_hemoglobin": False,
            "mitochondrial_prefix": "MT-",
            "ribosomal_prefixes": "RPS,RPL",
            "hemoglobin_pattern": r"^HB[^(P)]",
            "percent_top": "1,3",
            "log1p": True,
            "source": {"source": "layer", "source_layer": "counts"},
        },
    )

    assert [record["name"] for record in records] == ["adata", "summary", "code"]
    assert records[0] == {
        "type": "artifact",
        "name": "adata",
        "kind": ANNDATA_KIND,
        "codec": ANNDATA_CODEC,
        "payload": "outputs/adata",
    }
    json.dumps(records, allow_nan=False)
    assert _sha256(input_path) == before

    output = read_anndata(tmp_path / "calculate.partial" / records[0]["payload"])
    science.np.testing.assert_allclose(output.obs["total_counts"], [12.0, 5.0])
    assert output.var["mt"].tolist() == [True, False, False]
    assert output.var["ribo"].tolist() == [False, True, False]
    assert records[1]["type"] == "summary"
    assert records[1]["value"]["summary"]["node_id"] == "OpenBioSingleCellCalculateQC"
    assert records[2]["type"] == "string"
    compile(records[2]["value"], "<calculate-qc-code>", "exec")


def test_filter_cells_materializes_only_the_selected_subset(tmp_path: Path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    adata = science.ad.AnnData(
        science.sparse.csr_matrix([[-1.5, 0.0], [0.25, 0.25], [1.0, 2.0]]),
        obs=science.pd.DataFrame(
            {"pct_counts_mt": [0.0, 50.0, 10.0]},
            index=["negative", "fractional", "positive"],
        ),
        var=science.pd.DataFrame(index=["gene_a", "gene_b"]),
    )
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = _sha256(input_path)

    records = _execute(
        _context(tmp_path / "cells.partial"),
        "openbio.node.filtercells",
        {"adata": _artifact_descriptor(input_root)},
        {
            "min_genes": 0,
            "max_genes": 0,
            "min_counts": 0.0,
            "max_counts": 0.0,
            "max_pct_mito": 20.0,
            "mito_column": "pct_counts_mt",
            "source": {"source": "X"},
            "enable_min_counts": True,
            "enable_max_counts": False,
            "enable_max_pct_mito": True,
        },
    )

    assert [record["name"] for record in records] == ["adata", "summary", "code"]
    assert _sha256(input_path) == before
    output = read_anndata(tmp_path / "cells.partial" / records[0]["payload"])
    assert output.obs_names.tolist() == ["positive"]
    summary = records[1]["value"]["summary"]
    assert summary["key_results"]["removed_cells"] == 2
    assert summary["parameters"]["enable_min_counts"] is True
    compile(records[2]["value"], "<filter-cells-code>", "exec")


def test_filter_genes_maps_broader_raw_axis_and_preserves_unmodified_input(tmp_path: Path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    base = science.ad.AnnData(
        science.sparse.csr_matrix([[2.0, 0.0, 5.0], [1.0, 3.0, 0.0]]),
        obs=science.pd.DataFrame(index=["cell_a", "cell_b"]),
        var=science.pd.DataFrame(index=["keep", "drop", "raw_only"]),
    )
    base.raw = base
    adata = base[:, ["keep", "drop"]].copy()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = _sha256(input_path)

    records = _execute(
        _context(tmp_path / "genes.partial"),
        "openbio.node.filtergenes",
        {"adata": _artifact_descriptor(input_root)},
        {
            "min_cells": 2,
            "max_cells": 0,
            "min_counts": 0.0,
            "max_counts": 0.0,
            "source": {"source": "raw"},
            "enable_min_counts": False,
            "enable_max_counts": False,
        },
    )

    assert [record["name"] for record in records] == ["adata", "summary", "code"]
    assert _sha256(input_path) == before
    output = read_anndata(tmp_path / "genes.partial" / records[0]["payload"])
    assert output.var_names.tolist() == ["keep"]
    assert output.raw.var_names.tolist() == ["keep", "drop", "raw_only"]
    summary = records[1]["value"]["summary"]
    assert summary["key_results"]["evaluated_current_genes"] == 2
    assert summary["key_results"]["removed_genes"] == 1
    compile(records[2]["value"], "<filter-genes-code>", "exec")


def test_qc_plots_publishes_png_artifact_from_selected_source(tmp_path: Path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    adata = science.ad.AnnData(
        science.sparse.csr_matrix([[1.0, 1.0], [1.0, 1.0]]),
        obs=science.pd.DataFrame(
            {"total_counts": [900.0, 900.0], "pct_counts_mt": [99.0, 99.0]},
            index=["mitochondrial", "nuclear"],
        ),
        var=science.pd.DataFrame({"mt": [True, False]}, index=["MT-G", "G"]),
    )
    adata.layers["counts"] = science.sparse.csr_matrix([[10.0, 0.0], [0.0, 20.0]])
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = _sha256(input_path)

    records = _execute(
        _context(tmp_path / "plots.partial"),
        "openbio.node.qcplots",
        {"adata": _artifact_descriptor(input_root)},
        {
            "view": {"view": "overview"},
            "source": {"source": "layer", "source_layer": "counts"},
        },
    )

    assert records[0] == {
        "type": "artifact",
        "name": "plot",
        "kind": "OPENBIO_SINGLE_CELL_PLOT",
        "codec": "plot-png-v1",
        "payload": "outputs/plot",
    }
    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert _sha256(input_path) == before
    png, plot_metadata = read_plot(tmp_path / "plots.partial" / records[0]["payload"])
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert plot_metadata["kind"] == "plot"
    summary = records[1]["value"]["summary"]
    assert summary["parameters"]["view"] == {"view": "overview"}
    assert summary["key_results"]["distributions"]["total_expression"]["median"] == 15.0
    assert set(summary["key_results"]["metric_sources"].values()) == {"AnnData layer 'counts'"}
    compile(records[2]["value"], "<qc-plots-code>", "exec")
    exec(records[2]["value"], {})


def test_qc_worker_boundary_is_strict_and_node_classes_are_schema_only(tmp_path: Path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, science.ad.AnnData(science.np.eye(2)))
    descriptor = {**_artifact_descriptor(input_root), "unexpected": True}
    context = _context(tmp_path / "invalid.partial")

    with pytest.raises(ProtocolError, match="descriptor fields"):
        _execute(
            context,
            "openbio.node.qcplots",
            {"adata": descriptor},
            {"view": {"view": "overview"}, "source": {"source": "X"}},
        )

    source = (Path(__file__).parents[1] / "openbio_singlecell" / "operations_qc.py").read_text(encoding="utf-8")
    assert "comfy_api" not in source
    assert "folder_paths" not in source
    assert "nodes_" not in source
    assert all("execute" not in node.__dict__ for node in QC_NODE_CLASSES)
