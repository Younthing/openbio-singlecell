from __future__ import annotations

import copy
import json

import pytest

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.nodes_data import OpenBioSingleCellSnapshotExpression
from openbio_singlecell.nodes_preprocess import (
    OpenBioSingleCellHighlyVariableGenes,
    OpenBioSingleCellLog1p,
    OpenBioSingleCellNormalizeTotal,
    OpenBioSingleCellScale,
)
from openbio_singlecell.nodes_qc import (
    OpenBioSingleCellCalculateQC,
    OpenBioSingleCellFilterCells,
    OpenBioSingleCellFilterGenes,
    OpenBioSingleCellQCPlots,
)


def output_value(node_output):
    return node_output.result[0]


def dense(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


@pytest.fixture
def adata(science):
    rng = science.np.random.default_rng(4)
    counts = rng.poisson(2.0, size=(30, 12)).astype(float)
    counts[:15, :3] += 5
    counts[15:, 3:6] += 5
    obs = science.pd.DataFrame(index=[f"cell_{index}" for index in range(30)])
    var = science.pd.DataFrame(index=["MT-G0", *[f"G{index}" for index in range(1, 12)]])
    value = science.ad.AnnData(science.sparse.csr_matrix(counts), obs=obs, var=var)
    ensure_metadata(value, display_name="synthetic", source={"kind": "test"})
    return value


def assert_adata_equal(actual, expected, science):
    science.np.testing.assert_array_equal(dense(actual.X, science), dense(expected.X, science))
    science.pd.testing.assert_frame_equal(actual.obs, expected.obs)
    science.pd.testing.assert_frame_equal(actual.var, expected.var)
    assert set(actual.layers) == set(expected.layers)
    for key in actual.layers:
        science.np.testing.assert_array_equal(
            dense(actual.layers[key], science),
            dense(expected.layers[key], science),
        )
    assert actual.uns == expected.uns


def test_qc_filter_and_preprocessing_chain_preserves_counts(adata, science):
    qc = output_value(OpenBioSingleCellCalculateQC.execute(adata, mitochondrial_prefix="MT-"))
    assert "total_counts" in qc.obs
    filtered_cells = output_value(
        OpenBioSingleCellFilterCells.execute(qc, min_genes=1, min_counts=1, enable_min_counts=True)
    )
    filtered = output_value(
        OpenBioSingleCellFilterGenes.execute(filtered_cells, min_cells=1, min_counts=1, enable_min_counts=True)
    )
    snapshotted = output_value(OpenBioSingleCellSnapshotExpression.execute(filtered))
    normalized = output_value(OpenBioSingleCellNormalizeTotal.execute(snapshotted, 10_000.0))
    before_log = dense(normalized.layers["counts"], science).copy()
    logged = output_value(OpenBioSingleCellLog1p.execute(normalized))
    variable = output_value(
        OpenBioSingleCellHighlyVariableGenes.execute(
            logged,
            n_top_genes=8,
            flavor="seurat",
            source={"source": "X"},
            subset=False,
        )
    )
    with pytest.warns(UserWarning, match="densifies"):
        scaled = output_value(
            OpenBioSingleCellScale.execute(
                variable,
                source={"source": "X"},
                clipping_mode="custom",
                custom_max_value=10.0,
            )
        )

    assert science.np.array_equal(dense(logged.layers["counts"], science), before_log)
    assert logged.raw is not None
    assert int(variable.var["highly_variable"].sum()) == 8
    assert science.sparse.issparse(scaled.X)
    assert not science.sparse.issparse(scaled.layers["scaled"])
    assert len(scaled.uns["openbio_singlecell"]["analysis_history"]) == 8


def test_modifying_nodes_do_not_mutate_upstream_anndata(adata, science):
    snapshot = adata.copy()
    snapshot_uns = copy.deepcopy(adata.uns)

    OpenBioSingleCellCalculateQC.execute(adata, mitochondrial_prefix="MT-")
    OpenBioSingleCellFilterCells.execute(adata, min_genes=1, min_counts=1, enable_min_counts=True)
    OpenBioSingleCellFilterGenes.execute(adata, min_cells=1, min_counts=1, enable_min_counts=True)
    normalized = output_value(OpenBioSingleCellNormalizeTotal.execute(adata, 10_000.0))
    logged = output_value(OpenBioSingleCellLog1p.execute(normalized))
    OpenBioSingleCellHighlyVariableGenes.execute(
        logged,
        n_top_genes=8,
        flavor="seurat",
        source={"source": "X"},
        subset=False,
    )
    with pytest.warns(UserWarning, match="densifies"):
        OpenBioSingleCellScale.execute(
            logged,
            source={"source": "X"},
            clipping_mode="custom",
            custom_max_value=10.0,
        )

    assert_adata_equal(adata, snapshot, science)
    assert adata.uns == snapshot_uns


def test_qc_plot_is_read_only_and_reports_missing_mitochondrial_annotations(adata, science):
    adata.var_names = [f"G{index}" for index in range(adata.n_vars)]
    snapshot = adata.copy()

    qc = output_value(
        OpenBioSingleCellCalculateQC.execute(
            adata,
            include_ribosomal=False,
            include_hemoglobin=False,
            mitochondrial_prefix="MT-",
            percent_top="",
        )
    )
    plot = output_value(OpenBioSingleCellQCPlots.execute(qc))

    assert qc.uns["openbio_singlecell"]["warnings"] == ["No genes matched mitochondrial prefix 'MT-'."]
    assert "pct_counts_mt" not in qc.obs
    assert plot.kind == "plot"
    assert plot.png
    assert any("panel was omitted" in warning for warning in plot.warnings)
    assert_adata_equal(adata, snapshot, science)


def test_qc_plot_rejects_empty_anndata(science):
    empty = science.ad.AnnData(science.np.empty((0, 2)))

    with pytest.raises(ValueError, match="at least one cell and one gene"):
        OpenBioSingleCellQCPlots.execute(empty)


def test_qc_nodes_expose_json_serializable_reports_and_equivalent_code(adata):
    qc, qc_report, qc_code = OpenBioSingleCellCalculateQC.execute(adata).result
    cells, cell_report, cell_code = OpenBioSingleCellFilterCells.execute(qc, min_genes=1).result
    genes, gene_report, gene_code = OpenBioSingleCellFilterGenes.execute(cells, min_cells=1).result
    plot, plot_report, plot_code = OpenBioSingleCellQCPlots.execute(genes).result

    assert qc.shape == cells.shape == genes.shape
    assert plot.kind == "plot"
    for node_id, report, code in (
        ("OpenBioSingleCellCalculateQC", qc_report, qc_code),
        ("OpenBioSingleCellFilterCells", cell_report, cell_code),
        ("OpenBioSingleCellFilterGenes", gene_report, gene_code),
        ("OpenBioSingleCellQCPlots", plot_report, plot_code),
    ):
        assert report.kind == "summary"
        assert report.summary["schema_version"] == 1
        assert report.summary["node_id"] == node_id
        assert report.summary["methods"]
        assert report.summary["results"]
        assert report.summary["key_results"]
        assert report.summary["references"]
        assert report.summary["software_versions"]["openbio-singlecell"] == "0.2.0"
        json.dumps(report.summary, allow_nan=False)
        assert code.endswith("\n")
        compile(code, f"<{node_id}-code>", "exec")

    _, mitochondrial_report, _ = OpenBioSingleCellFilterCells.execute(
        qc,
        max_pct_mito=100.0,
        enable_max_pct_mito=True,
    ).result
    mitochondrial_metric = mitochondrial_report.summary["key_results"]["mitochondrial_metric"]
    assert mitochondrial_metric["available"] is True
    assert mitochondrial_metric["used_for_filtering"] is True
    assert mitochondrial_metric["distribution"]["n"] == adata.n_obs
    assert mitochondrial_report.summary["key_results"]["zero_total_expression_cells"] == 0


def test_qc_count_source_is_explicit_and_validated(adata, science):
    expected_counts = adata.X.copy()
    adata.layers["counts"] = expected_counts
    adata.X = science.np.full(adata.shape, -1.0)

    with pytest.raises(ValueError, match="log1p annotations are undefined"):
        OpenBioSingleCellCalculateQC.execute(adata)

    signed = adata.copy()
    signed.X = expected_counts.toarray() if science.sparse.issparse(expected_counts) else expected_counts.copy()
    signed.X[0, 0] = -0.25
    _, signed_report, _ = OpenBioSingleCellCalculateQC.execute(signed, percent_top="").result
    assert any("negative expression" in warning for warning in signed_report.warnings)

    qc, report, _ = OpenBioSingleCellCalculateQC.execute(
        adata,
        source={"source": "layer", "source_layer": "counts"},
    ).result

    expected_totals = science.np.asarray(expected_counts.sum(axis=1)).ravel()
    science.np.testing.assert_allclose(qc.obs["total_counts"], expected_totals)
    assert report.summary["parameters"]["source"] == "layer"
    assert report.summary["parameters"]["source_layer"] == "counts"


def test_qc_signed_all_feature_percentages_preserve_math_and_disclose_zero_denominators(science):
    value = science.ad.AnnData(
        science.np.asarray(
            [
                [-2.0, 0.5],
                [-1.0, 1.0],
                [1.0, 1.0],
            ]
        ),
        obs=science.pd.DataFrame(index=["negative_total", "zero_total", "positive_total"]),
        var=science.pd.DataFrame(index=["gene_0", "gene_1"]),
    )

    output, report, code = OpenBioSingleCellCalculateQC.execute(
        value,
        include_ribosomal=False,
        include_hemoglobin=False,
        percent_top="3",
        log1p=False,
    ).result

    observed = science.np.asarray(output.obs["pct_counts_in_top_3_genes"], dtype=float)
    science.np.testing.assert_allclose(observed, [100.0, science.np.nan, 100.0], equal_nan=True)
    assert report.summary["key_results"]["undefined_all_feature_percent_cells"] == 1
    assert any("zero selected-source expression sum" in warning for warning in report.warnings)
    json.dumps(report.summary, allow_nan=False)

    namespace = {}
    exec(code, namespace)
    equivalent = namespace["calculate_qc_metrics"](value)
    science.np.testing.assert_allclose(
        equivalent.obs["pct_counts_in_top_3_genes"],
        observed,
        equal_nan=True,
    )


def test_qc_expression_threshold_enable_controls_accept_signed_fractional_and_zero(science):
    value = science.ad.AnnData(
        science.np.asarray(
            [
                [-1.5, 0.0, 0.0],
                [0.25, 0.25, 0.0],
                [1.0, 1.0, 1.0],
            ]
        ),
        obs=science.pd.DataFrame(
            {"pct_counts_mt": [0.0, 0.1, 0.0]},
            index=["negative", "fractional", "positive"],
        ),
        var=science.pd.DataFrame(index=["negative_sum", "fractional_sum", "positive_sum"]),
    )

    cells, cell_report, cell_code = OpenBioSingleCellFilterCells.execute(
        value,
        min_counts=0.0,
        enable_min_counts=True,
    ).result
    assert cells.obs_names.tolist() == ["fractional", "positive"]
    assert cell_report.summary["parameters"]["min_counts"] == 0.0
    assert cell_report.summary["parameters"]["enable_min_counts"] is True
    assert any("negative expression" in warning for warning in cell_report.warnings)
    namespace = {}
    exec(cell_code, namespace)
    assert namespace["filter_cells"](value).obs_names.equals(cells.obs_names)

    upper_cells, upper_report, _ = OpenBioSingleCellFilterCells.execute(
        value,
        max_counts=0.5,
        enable_max_counts=True,
    ).result
    assert upper_cells.obs_names.tolist() == ["negative", "fractional"]
    assert upper_report.summary["parameters"]["max_counts"] == 0.5

    mitochondrial_cells, mitochondrial_report, mitochondrial_code = OpenBioSingleCellFilterCells.execute(
        value,
        max_pct_mito=0.0,
        enable_max_pct_mito=True,
    ).result
    assert mitochondrial_cells.obs_names.tolist() == ["negative", "positive"]
    assert mitochondrial_report.summary["key_results"]["mitochondrial_metric"]["used_for_filtering"] is True
    namespace = {}
    exec(mitochondrial_code, namespace)
    assert namespace["filter_cells"](value).obs_names.equals(mitochondrial_cells.obs_names)

    genes, gene_report, gene_code = OpenBioSingleCellFilterGenes.execute(
        value,
        min_counts=0.0,
        enable_min_counts=True,
    ).result
    assert genes.var_names.tolist() == ["fractional_sum", "positive_sum"]
    assert gene_report.summary["parameters"]["enable_min_counts"] is True
    namespace = {}
    exec(gene_code, namespace)
    assert namespace["filter_genes"](value).var_names.equals(genes.var_names)

    disabled, disabled_report, disabled_code = OpenBioSingleCellFilterCells.execute(
        value,
        min_counts=science.np.nan,
        enable_min_counts=False,
    ).result
    assert disabled.obs_names.equals(value.obs_names)
    assert disabled_report.summary["parameters"]["min_counts"] is None
    json.dumps(disabled_report.summary, allow_nan=False)
    namespace = {}
    exec(disabled_code, namespace)
    assert namespace["filter_cells"](value).obs_names.equals(value.obs_names)

    with pytest.raises(ValueError, match="minimum cannot exceed"):
        OpenBioSingleCellFilterGenes.execute(
            value,
            min_counts=1.0,
            max_counts=0.0,
            enable_min_counts=True,
            enable_max_counts=True,
        )


def test_qc_sparse_explicit_zeros_are_not_counted_as_detected(science):
    matrix = science.sparse.csr_matrix(
        ([0.0, 5.0], ([0, 1], [0, 1])),
        shape=(2, 2),
    )
    value = science.ad.AnnData(
        matrix,
        obs=science.pd.DataFrame(index=["explicit_zero", "expressed"]),
        var=science.pd.DataFrame(index=["zero_gene", "expressed_gene"]),
    )
    assert value.X.nnz == 2

    qc, _, _ = OpenBioSingleCellCalculateQC.execute(value, percent_top="").result
    science.np.testing.assert_array_equal(qc.obs["n_genes_by_counts"], [0, 1])
    science.np.testing.assert_array_equal(qc.var["n_cells_by_counts"], [0, 1])

    cells, _, cells_code = OpenBioSingleCellFilterCells.execute(value, min_genes=1).result
    assert cells.obs_names.tolist() == ["expressed"]
    namespace = {}
    exec(cells_code, namespace)
    assert namespace["filter_cells"](value).obs_names.tolist() == ["expressed"]

    genes, _, genes_code = OpenBioSingleCellFilterGenes.execute(value, min_cells=1).result
    assert genes.var_names.tolist() == ["expressed_gene"]
    namespace = {}
    exec(genes_code, namespace)
    assert namespace["filter_genes"](value).var_names.tolist() == ["expressed_gene"]


def test_qc_sparse_duplicate_storage_is_not_mutated(science):
    matrix = science.sparse.csr_matrix(
        (
            science.np.asarray([1.0, 1.0, 2.0]),
            science.np.asarray([0, 0, 1]),
            science.np.asarray([0, 2, 3]),
        ),
        shape=(2, 2),
    )
    value = science.ad.AnnData(
        matrix,
        obs=science.pd.DataFrame(index=["duplicated", "expressed"]),
        var=science.pd.DataFrame(index=["duplicated_gene", "expressed_gene"]),
    )
    snapshot = (
        value.X.data.copy(),
        value.X.indices.copy(),
        value.X.indptr.copy(),
        value.X.has_canonical_format,
    )

    def assert_storage_unchanged():
        science.np.testing.assert_array_equal(value.X.data, snapshot[0])
        science.np.testing.assert_array_equal(value.X.indices, snapshot[1])
        science.np.testing.assert_array_equal(value.X.indptr, snapshot[2])
        assert value.X.has_canonical_format is snapshot[3]

    assert value.X.has_canonical_format is False
    OpenBioSingleCellCalculateQC.execute(value, percent_top="")
    assert_storage_unchanged()
    runtime_cells, _, cells_code = OpenBioSingleCellFilterCells.execute(value, min_genes=1).result
    assert_storage_unchanged()
    OpenBioSingleCellFilterGenes.execute(value, min_cells=1)
    assert_storage_unchanged()
    _, _, plot_code = OpenBioSingleCellQCPlots.execute(value).result
    assert_storage_unchanged()

    namespace = {}
    exec(plot_code, namespace)
    namespace["qc_plots"](value)
    assert_storage_unchanged()

    namespace = {}
    exec(cells_code, namespace)
    equivalent_cells = namespace["filter_cells"](value)
    assert_storage_unchanged()
    science.np.testing.assert_array_equal(equivalent_cells.X.data, runtime_cells.X.data)
    science.np.testing.assert_array_equal(equivalent_cells.X.indices, runtime_cells.X.indices)
    science.np.testing.assert_array_equal(equivalent_cells.X.indptr, runtime_cells.X.indptr)
    assert equivalent_cells.X.has_canonical_format is runtime_cells.X.has_canonical_format


def test_qc_filters_reject_contradictory_bounds_and_report_empty_results(adata):
    with pytest.raises(ValueError, match="minimum cannot exceed"):
        OpenBioSingleCellFilterCells.execute(adata, min_genes=10, max_genes=5)
    with pytest.raises(ValueError, match="minimum cannot exceed"):
        OpenBioSingleCellFilterGenes.execute(adata, min_cells=10, max_cells=5)
    empty_cells, cell_report, _ = OpenBioSingleCellFilterCells.execute(
        adata,
        min_genes=1_000_000,
    ).result
    empty_genes, gene_report, _ = OpenBioSingleCellFilterGenes.execute(
        adata,
        min_cells=1_000_000,
    ).result
    assert empty_cells.n_obs == 0
    assert empty_genes.n_vars == 0
    assert any("removed every cell" in warning for warning in cell_report.warnings)
    assert any("removed every gene" in warning for warning in gene_report.warnings)


def test_qc_filter_code_reproduces_empty_results(science):
    empty_counts = science.ad.AnnData(
        science.sparse.csr_matrix((2, 2), dtype=float),
        obs=science.pd.DataFrame(index=["cell_0", "cell_1"]),
        var=science.pd.DataFrame(index=["gene_0", "gene_1"]),
    )

    for node, function_name, keyword, empty_axis in (
        (OpenBioSingleCellFilterCells, "filter_cells", {"min_genes": 1}, "obs"),
        (OpenBioSingleCellFilterGenes, "filter_genes", {"min_cells": 1}, "var"),
    ):
        output, report, code = node.execute(empty_counts, **keyword).result
        assert len(getattr(output, empty_axis + "_names")) == 0
        assert any("removed every" in warning for warning in report.warnings)
        namespace = {}
        exec(code, namespace)
        equivalent = namespace[function_name](empty_counts)
        assert len(getattr(equivalent, empty_axis + "_names")) == 0


def test_qc_plot_reports_missing_mitochondrial_metric_as_unavailable(adata):
    adata.var_names = [f"G{index}" for index in range(adata.n_vars)]

    _, report, code = OpenBioSingleCellQCPlots.execute(adata).result

    assert "mitochondrial_percent" not in report.summary["key_results"]["distributions"]
    assert "total_expression_vs_mitochondrial_percent" not in report.summary["key_results"]["panels"]
    assert "Mitochondrial metric unavailable" in code


def test_qc_plot_derived_mitochondrial_percent_uses_one_expression_source(science):
    value = science.ad.AnnData(
        science.sparse.csr_matrix([[10.0, 0.0], [0.0, 10.0]]),
        obs=science.pd.DataFrame(
            {
                "total_counts": [100.0, 100.0],
                "n_genes_by_counts": [1, 1],
            },
            index=["mitochondrial", "nuclear"],
        ),
        var=science.pd.DataFrame({"mt": [True, False]}, index=["MT-G", "G"]),
    )

    _, report, code = OpenBioSingleCellQCPlots.execute(value).result

    assert report.summary["key_results"]["distributions"]["mitochondrial_percent"]["median"] == 50.0
    assert report.summary["key_results"]["metric_sources"]["mitochondrial_panel_total_expression"] == "AnnData.X"
    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](value)
    offsets = figure.axes[3].collections[0].get_offsets()
    science.np.testing.assert_allclose(offsets[:, 0], [10.0, 10.0])
    science.np.testing.assert_allclose(offsets[:, 1], [100.0, 0.0])


def test_qc_plot_omits_nonfinite_mitochondrial_evidence(adata, science):
    qc = OpenBioSingleCellCalculateQC.execute(adata).result[0]
    qc.obs["pct_counts_mt"] = science.np.nan

    _, report, code = OpenBioSingleCellQCPlots.execute(qc).result

    key_results = report.summary["key_results"]
    assert "mitochondrial_percent" not in key_results["distributions"]
    assert "total_expression_vs_mitochondrial_percent" not in key_results["panels"]
    assert any("no finite values" in warning for warning in report.warnings)
    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](qc)
    assert figure.axes[3].axison is False
    assert figure.axes[3].texts[0].get_text() == "Mitochondrial metric unavailable"


def test_qc_plot_reports_partially_nonfinite_mitochondrial_evidence(adata, science):
    qc = OpenBioSingleCellCalculateQC.execute(adata).result[0]
    qc.obs.loc[qc.obs_names[0], "pct_counts_mt"] = science.np.nan

    _, report, _ = OpenBioSingleCellQCPlots.execute(qc).result

    distribution = report.summary["key_results"]["distributions"]["mitochondrial_percent"]
    assert distribution["n"] == qc.n_obs
    assert distribution["missing"] == 1
    assert any("non-finite mitochondrial percentages" in warning for warning in report.warnings)


def test_qc_plot_code_reproduces_nonfinite_primary_metric_error(adata, science):
    _, _, code = OpenBioSingleCellQCPlots.execute(adata).result
    invalid = adata.copy()
    invalid.obs["total_counts"] = science.np.nan
    invalid.obs["n_genes_by_counts"] = science.np.inf

    with pytest.raises(ValueError, match="at least one finite"):
        OpenBioSingleCellQCPlots.execute(invalid)
    namespace = {}
    exec(code, namespace)
    with pytest.raises(ValueError, match="at least one finite"):
        namespace["qc_plots"](invalid)


def test_qc_equivalent_code_reproduces_primary_results(adata, science):
    qc, _, qc_code = OpenBioSingleCellCalculateQC.execute(adata, mitochondrial_prefix="mt-").result
    namespace = {}
    exec(qc_code, namespace)
    equivalent_qc = namespace["calculate_qc_metrics"](adata)
    science.pd.testing.assert_frame_equal(qc.obs, equivalent_qc.obs)
    science.pd.testing.assert_frame_equal(qc.var, equivalent_qc.var)

    filtered_cells, _, cells_code = OpenBioSingleCellFilterCells.execute(qc, min_genes=3).result
    namespace = {}
    exec(cells_code, namespace)
    equivalent_cells = namespace["filter_cells"](qc)
    assert equivalent_cells.obs_names.equals(filtered_cells.obs_names)

    filtered_genes, _, genes_code = OpenBioSingleCellFilterGenes.execute(filtered_cells, min_cells=2).result
    namespace = {}
    exec(genes_code, namespace)
    equivalent_genes = namespace["filter_genes"](filtered_cells)
    assert equivalent_genes.var_names.equals(filtered_genes.var_names)

    _, _, plot_code = OpenBioSingleCellQCPlots.execute(filtered_genes).result
    namespace = {}
    exec(plot_code, namespace)
    equivalent_figure = namespace["qc_plots"](filtered_genes)
    assert len(equivalent_figure.axes) == 4
    assert equivalent_figure.axes[0].get_title() == "Total expression per cell"
