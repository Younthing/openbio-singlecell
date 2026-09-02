from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.nodes_embedding import OpenBioSingleCellPCAVariancePlot
from openbio_singlecell.operations_embedding import pca_variance_plot, pca_variance_plot_owned
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def test_pca_variance_plot_schema():
    schema = OpenBioSingleCellPCAVariancePlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellPCAVariancePlot"
    assert schema.display_name == "PCA Variance Plot"
    assert schema.category == "openbio/single-cell/dimension-reduction"
    assert [item.id for item in schema.inputs] == ["adata", "n_pcs"]
    assert schema.inputs[1].default == 0
    assert schema.inputs[1].min == 0
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_pca_variance_plot_reads_all_stored_ratios_without_recomputing_or_mutating(science, monkeypatch):
    adata = science.ad.AnnData(science.np.arange(20, dtype=float).reshape(5, 4))
    adata.uns["pca"] = {
        "params": {"zero_center": True},
        "variance_ratio": science.np.asarray([0.4, 0.25, 0.1], dtype=float),
    }
    before = adata.copy()

    def reject_recomputation(*args, **kwargs):
        raise AssertionError("PCA Variance Plot must not recompute PCA")

    monkeypatch.setattr(science.sc.pp, "pca", reject_recomputation)
    plotted, report, code = pca_variance_plot_owned(adata)

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellPCAVariancePlot"
    assert summary["parameters"] == {"n_pcs": 0}
    assert summary["key_results"]["variance_ratio_source"] == "adata.uns['pca']['variance_ratio']"
    assert summary["key_results"]["available_components"] == 3
    assert summary["key_results"]["plotted_components"] == 3
    assert summary["key_results"]["variance_ratio"] == [0.4, 0.25, 0.1]
    assert summary["key_results"]["cumulative_variance_ratio"] == [0.4, 0.65, 0.75]
    json.dumps(summary, allow_nan=False)

    compile(code, "<pca-variance-plot-code>", "exec")
    namespace: dict[str, object] = {}
    exec(code, namespace)
    assert namespace["plot_pca_variance"](adata) == plotted.png

    science.np.testing.assert_array_equal(adata.X, before.X)
    science.pd.testing.assert_frame_equal(adata.obs, before.obs)
    science.pd.testing.assert_frame_equal(adata.var, before.var)
    assert list(adata.uns) == list(before.uns)
    assert adata.uns["pca"]["params"] == before.uns["pca"]["params"]
    science.np.testing.assert_array_equal(
        adata.uns["pca"]["variance_ratio"], before.uns["pca"]["variance_ratio"]
    )


def test_pca_variance_plot_worker_round_trips_png_without_rewriting_input(tmp_path, science):
    adata = science.ad.AnnData(science.np.ones((4, 3), dtype=float))
    adata.uns["pca"] = {"variance_ratio": science.np.asarray([0.5, 0.2], dtype=float)}
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = pca_variance_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"n_pcs": 1},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert records[0]["kind"] == "OPENBIO_SINGLE_CELL_PLOT"
    assert records[0]["codec"] == "plot-png-v1"
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"n_pcs": 1}
    assert hashlib.sha256(input_path.read_bytes()).digest() == before


@pytest.mark.parametrize(
    ("case", "n_pcs", "error", "message"),
    [
        ("missing", 0, ValueError, "requires stored"),
        ("matrix", 0, ValueError, "one-dimensional"),
        ("boolean", 0, TypeError, "real numeric"),
        ("nonfinite", 0, ValueError, "non-finite"),
        ("negative", 0, ValueError, "negative"),
        ("total_over_one", 0, ValueError, "cannot exceed 1"),
        ("valid", 3, ValueError, "only 2 are stored"),
        ("valid", True, TypeError, "must be an integer"),
    ],
)
def test_pca_variance_plot_rejects_invalid_stored_evidence(science, case, n_pcs, error, message):
    adata = science.ad.AnnData(science.np.ones((3, 2), dtype=float))
    values = {
        "matrix": [[0.5, 0.2]],
        "boolean": [True, False],
        "nonfinite": [0.5, science.np.nan],
        "negative": [0.5, -0.1],
        "total_over_one": [0.8, 0.3],
        "valid": [0.5, 0.2],
    }
    if case != "missing":
        adata.uns["pca"] = {"variance_ratio": science.np.asarray(values[case])}

    with pytest.raises(error, match=message):
        pca_variance_plot_owned(adata, n_pcs=n_pcs)
