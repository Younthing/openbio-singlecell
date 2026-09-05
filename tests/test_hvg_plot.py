from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata
from openbio_singlecell.contracts import record_history
from openbio_singlecell.nodes_preprocess import OpenBioSingleCellHVGSelectionPlot
from openbio_singlecell.operations_preprocess import hvg_selection_plot_owned
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

FLAVOR_COLUMNS = [
    ("seurat", "dispersions", "dispersions_norm"),
    ("cell_ranger", "dispersions", "dispersions_norm"),
    ("seurat_v3", "variances", "variances_norm"),
    ("seurat_v3_paper", "variances", "variances_norm"),
    ("pearson_residuals", "variances", "residual_variances"),
]


def _hvg_adata(science, flavor: str, base_variability_column: str, selection_variability_column: str):
    adata = science.ad.AnnData(
        X=science.np.zeros((3, 6), dtype=float),
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(6)]),
    )
    adata.var["means"] = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
    adata.var[base_variability_column] = [0.2, 0.5, 0.4, 1.2, 0.8, 0.3]
    adata.var[selection_variability_column] = [0.1, 1.5, 0.2, 2.0, 0.5, 0.4]
    algorithm = science.np.asarray([True, False, False, True, False, False])
    forced = science.np.asarray([False, True, False, False, False, False])
    adata.var["highly_variable_algorithm"] = algorithm
    adata.var["highly_variable_forced"] = forced
    adata.var["highly_variable"] = algorithm | forced
    adata.uns["hvg"] = {"flavor": flavor}
    return adata


def test_hvg_selection_plot_schema_is_read_only_and_has_standard_outputs():
    schema = OpenBioSingleCellHVGSelectionPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellHVGSelectionPlot"
    assert schema.display_name == "HVG Selection Plot"
    assert schema.category == "openbio/single-cell/preprocessing"
    assert [item.id for item in schema.inputs] == ["adata"]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]
    assert [item.get_io_type() for item in schema.outputs] == [
        "OPENBIO_SINGLE_CELL_PLOT",
        "OPENBIO_SINGLE_CELL_SUMMARY",
        "STRING",
    ]


@pytest.mark.parametrize(("flavor", "base_variability_column", "selection_variability_column"), FLAVOR_COLUMNS)
def test_hvg_selection_plot_renders_stored_flavor_metrics_and_selection_classes(
    science,
    flavor,
    base_variability_column,
    selection_variability_column,
):
    adata = _hvg_adata(science, flavor, base_variability_column, selection_variability_column)
    original = adata.copy()

    plotted, report, code = hvg_selection_plot_owned(adata)

    assert plotted.kind == "plot"
    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["node_id"] == "OpenBioSingleCellHVGSelectionPlot"
    assert report.summary["key_results"] == {
        "flavor": flavor,
        "mean_column": "means",
        "variability_columns": [selection_variability_column, base_variability_column],
        "selection_counts": {
            "algorithm_selected": 2,
            "forced": 1,
            "unselected": 3,
        },
        "input_features": 6,
        "plotted_points_by_panel": {
            selection_variability_column: 6,
            base_variability_column: 6,
        },
        "omitted_features_by_panel": {
            selection_variability_column: 0,
            base_variability_column: 0,
        },
        "feature_axis_validation": "unverified",
        "warnings": [
            "HVG Selection Plot feature-axis completeness is unverified because no OpenBio "
            "highly_variable_genes history is available."
        ],
    }
    assert "Classified 2 algorithm-selected, 1 forced, and 3 unselected features" in report.summary["results"]
    assert "Displayed" not in report.summary["results"]
    json.dumps(report.summary, allow_nan=False)
    compile(code, "<hvg-selection-plot-code>", "exec")
    science.pd.testing.assert_frame_equal(adata.var, original.var)
    assert adata.uns == original.uns
    science.np.testing.assert_array_equal(adata.X, original.X)


@pytest.mark.parametrize("flavor", ["seurat", "cell_ranger"])
def test_hvg_selection_plot_omits_valid_nan_variability_by_panel(science, flavor):
    adata = _hvg_adata(science, flavor, "dispersions", "dispersions_norm")
    adata.var.loc["gene_2", "dispersions_norm"] = science.np.nan
    adata.var.loc["gene_4", "dispersions"] = science.np.nan
    record_history(adata, "highly_variable_genes", {}, adata.n_obs, adata.n_vars)

    plotted, report, _ = hvg_selection_plot_owned(adata)

    expected_omissions = {"dispersions_norm": 1, "dispersions": 1}
    expected_plotted = {"dispersions_norm": 5, "dispersions": 5}
    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["key_results"]["input_features"] == 6
    assert report.summary["key_results"]["plotted_points_by_panel"] == expected_plotted
    assert report.summary["key_results"]["omitted_features_by_panel"] == expected_omissions
    assert report.summary["key_results"]["warnings"] == plotted.warnings
    assert len(plotted.warnings) == 2
    assert plotted.warnings == report.summary["warnings"]
    assert all("omitted 1 feature" in warning for warning in plotted.warnings)


def test_hvg_selection_plot_rejects_panel_without_any_finite_points(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var["dispersions"] = science.np.nan

    with pytest.raises(ValueError, match="no finite points.*dispersions"):
        hvg_selection_plot_owned(adata)


def test_hvg_selection_plot_rejects_nonfinite_algorithm_selection_metric(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var.loc["gene_0", "dispersions_norm"] = science.np.nan

    with pytest.raises(ValueError, match="algorithm-selected.*finite.*dispersions_norm"):
        hvg_selection_plot_owned(adata)


@pytest.mark.parametrize("column", ["means", "dispersions", "dispersions_norm"])
def test_hvg_selection_plot_rejects_infinite_metric_evidence(science, column):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var.loc["gene_2", column] = science.np.inf

    with pytest.raises(ValueError, match=rf"finite numeric evidence.*{column}"):
        hvg_selection_plot_owned(adata)


@pytest.mark.parametrize(
    ("flavor", "base_variability_column", "selection_variability_column"),
    FLAVOR_COLUMNS[2:],
)
def test_hvg_selection_plot_rejects_nan_variability_for_non_log_flavors(
    science,
    flavor,
    base_variability_column,
    selection_variability_column,
):
    adata = _hvg_adata(science, flavor, base_variability_column, selection_variability_column)
    adata.var.loc["gene_2", base_variability_column] = science.np.nan

    with pytest.raises(ValueError, match=rf"finite numeric evidence.*{base_variability_column}"):
        hvg_selection_plot_owned(adata)


def test_hvg_selection_plot_rejects_missing_stored_flavor(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    del adata.uns["hvg"]

    with pytest.raises(ValueError, match="stored flavor metadata"):
        hvg_selection_plot_owned(adata)


def test_hvg_selection_plot_rejects_missing_flavor_metric(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    del adata.var["dispersions_norm"]

    with pytest.raises(ValueError, match="complete stored evidence.*dispersions_norm"):
        hvg_selection_plot_owned(adata)


@pytest.mark.parametrize(("with_history", "expected_validation"), [(False, "unverified"), (True, "verified")])
def test_hvg_selection_plot_allows_all_selected_evidence(
    science,
    with_history,
    expected_validation,
):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var["highly_variable_algorithm"] = True
    adata.var["highly_variable_forced"] = False
    adata.var["highly_variable"] = True
    if with_history:
        record_history(adata, "highly_variable_genes", {}, adata.n_obs, adata.n_vars)

    plotted, report, _ = hvg_selection_plot_owned(adata)

    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["key_results"]["selection_counts"]["unselected"] == 0
    assert report.summary["key_results"]["feature_axis_validation"] == expected_validation


def test_hvg_selection_plot_discloses_partial_feature_slice_and_preserves_generated_behavior(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    record_history(adata, "highly_variable_genes", {}, adata.n_obs, adata.n_vars)
    sliced = adata[:, :5].copy()
    assert bool((~sliced.var["highly_variable"].to_numpy()).any())

    plotted, report, code = hvg_selection_plot_owned(sliced)

    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["key_results"]["input_features"] == 5
    assert report.summary["key_results"]["feature_axis_validation"] == "different_from_recorded"
    assert any("6" in warning and "5" in warning for warning in plotted.warnings)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_hvg_selection"](sliced) == plotted.png


def test_hvg_selection_plot_accepts_native_scanpy_selection_without_inventing_origin(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var.drop(columns=["highly_variable_algorithm", "highly_variable_forced"], inplace=True)
    original = adata.copy()

    plotted, report, code = hvg_selection_plot_owned(adata)

    assert plotted.png.startswith(b"\x89PNG\r\n\x1a\n")
    assert report.summary["key_results"]["selection_counts"] == {
        "selected": 3,
        "algorithm_selected": None,
        "forced": None,
        "unselected": 3,
    }
    assert "3 selected and 3 unselected" in report.summary["results"]
    assert any("selection origin" in warning for warning in plotted.warnings)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_hvg_selection"](adata) == plotted.png
    science.pd.testing.assert_frame_equal(adata.var, original.var)
    assert adata.uns == original.uns


def test_hvg_selection_plot_marks_external_feature_axis_unverified(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")

    plotted, report, _ = hvg_selection_plot_owned(adata)

    assert report.summary["key_results"]["feature_axis_validation"] == "unverified"
    assert any("feature-axis completeness is unverified" in warning for warning in plotted.warnings)
    assert plotted.warnings == report.summary["warnings"]


def test_hvg_selection_plot_marks_complete_openbio_feature_axis_verified(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    record_history(adata, "highly_variable_genes", {}, adata.n_obs, adata.n_vars)

    plotted, report, _ = hvg_selection_plot_owned(adata)

    assert report.summary["key_results"]["feature_axis_validation"] == "verified"
    assert not plotted.warnings


def test_hvg_selection_plot_rejects_inconsistent_selection_masks(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var.loc["gene_2", "highly_variable"] = True

    with pytest.raises(ValueError, match="selection masks are inconsistent"):
        hvg_selection_plot_owned(adata)


def test_hvg_selection_plot_rejects_nonfinite_metric_evidence(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var.loc["gene_2", "means"] = science.np.nan

    with pytest.raises(ValueError, match="finite numeric evidence.*means"):
        hvg_selection_plot_owned(adata)


def test_hvg_selection_plot_rejects_nonboolean_selection_evidence(science):
    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    adata.var["highly_variable_algorithm"] = [1, 0, 0, 1, 0, 0]

    with pytest.raises(ValueError, match="boolean evidence.*highly_variable_algorithm"):
        hvg_selection_plot_owned(adata)


def test_hvg_selection_plot_does_not_change_matplotlib_backend(science):
    import matplotlib

    adata = _hvg_adata(science, "seurat", "dispersions", "dispersions_norm")
    original_backend = matplotlib.get_backend()
    matplotlib.use("svg", force=True)
    try:
        backend_before = matplotlib.get_backend()
        hvg_selection_plot_owned(adata)
        assert matplotlib.get_backend() == backend_before
    finally:
        matplotlib.use(original_backend, force=True)


def test_hvg_selection_plot_worker_round_trips_png_without_mutating_input(tmp_path, science):
    from openbio_singlecell.operations_preprocess import hvg_selection_plot

    adata = _hvg_adata(science, "pearson_residuals", "variances", "residual_variances")
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "hvg-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = hvg_selection_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        {},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert records[0]["kind"] == "OPENBIO_SINGLE_CELL_PLOT"
    assert records[0]["codec"] == "plot-png-v1"
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert metadata["kind"] == "plot"
    assert hashlib.sha256(input_path.read_bytes()).digest() == before
    json.dumps(records[1]["value"]["summary"], allow_nan=False)
    namespace = {}
    exec(records[2]["value"], namespace)
    assert namespace["plot_hvg_selection"](adata) == png
    assert adata.uns == {"hvg": {"flavor": "pearson_residuals"}}
