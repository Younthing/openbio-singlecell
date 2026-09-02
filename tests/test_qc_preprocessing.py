from __future__ import annotations

import copy
import json

import pytest

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.operations_qc import (
    calculate_qc_owned,
    filter_cells_owned,
    filter_genes_owned,
    qc_plots_owned,
)


def output_value(node_output):
    return node_output[0]


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


def test_qc_filter_chain_preserves_counts(adata, science):
    expected = dense(adata.X, science).copy()
    qc = output_value(calculate_qc_owned(adata, mitochondrial_prefix="MT-"))
    assert "total_counts" in qc.obs
    filtered_cells = output_value(
        filter_cells_owned(qc, min_genes=1, min_counts=1, enable_min_counts=True)
    )
    filtered = output_value(
        filter_genes_owned(filtered_cells, min_cells=1, min_counts=1, enable_min_counts=True)
    )
    science.np.testing.assert_array_equal(dense(filtered.X, science), expected)
    assert len(filtered.uns["openbio_singlecell"]["analysis_history"]) == 3


def test_qc_owned_operations_mutate_only_worker_private_copy(adata, science):
    snapshot = adata.copy()
    snapshot_uns = copy.deepcopy(adata.uns)

    calculate_qc_owned(adata.copy(), mitochondrial_prefix="MT-")
    filter_cells_owned(adata.copy(), min_genes=1, min_counts=1, enable_min_counts=True)
    filter_genes_owned(adata.copy(), min_cells=1, min_counts=1, enable_min_counts=True)

    assert_adata_equal(adata, snapshot, science)
    assert adata.uns == snapshot_uns


def test_qc_plot_is_read_only_and_reports_missing_mitochondrial_annotations(adata, science):
    adata.var_names = [f"G{index}" for index in range(adata.n_vars)]
    snapshot = adata.copy()

    qc = output_value(
        calculate_qc_owned(
            adata.copy(),
            include_ribosomal=False,
            include_hemoglobin=False,
            mitochondrial_prefix="MT-",
            percent_top="",
        )
    )
    plot = output_value(qc_plots_owned(qc))

    assert qc.uns["openbio_singlecell"]["warnings"] == ["No genes matched mitochondrial prefix 'MT-'."]
    assert "pct_counts_mt" not in qc.obs
    assert plot.kind == "plot"
    assert plot.png
    assert any("panel was omitted" in warning for warning in plot.warnings)
    assert_adata_equal(adata, snapshot, science)


def test_qc_plot_rejects_empty_anndata(science):
    empty = science.ad.AnnData(science.np.empty((0, 2)))

    with pytest.raises(ValueError, match="at least one cell and one gene"):
        qc_plots_owned(empty)


def test_qc_nodes_expose_json_serializable_reports_and_equivalent_code(adata):
    qc, qc_report, qc_code = calculate_qc_owned(adata)
    cells, cell_report, cell_code = filter_cells_owned(qc, min_genes=1)
    genes, gene_report, gene_code = filter_genes_owned(cells, min_cells=1)
    plot, plot_report, plot_code = qc_plots_owned(genes)

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

    _, mitochondrial_report, _ = filter_cells_owned(
        qc,
        max_pct_mito=100.0,
        enable_max_pct_mito=True,
    )
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
        calculate_qc_owned(adata)

    signed = adata.copy()
    signed.X = expected_counts.toarray() if science.sparse.issparse(expected_counts) else expected_counts.copy()
    signed.X[0, 0] = -0.25
    _, signed_report, _ = calculate_qc_owned(signed, percent_top="")
    assert any("negative expression" in warning for warning in signed_report.warnings)

    qc, report, _ = calculate_qc_owned(
        adata,
        source={"source": "layer", "source_layer": "counts"},
    )

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

    output, report, code = calculate_qc_owned(
        value,
        include_ribosomal=False,
        include_hemoglobin=False,
        percent_top="3",
        log1p=False,
    )

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

    cells, cell_report, cell_code = filter_cells_owned(
        value,
        min_counts=0.0,
        enable_min_counts=True,
    )
    assert cells.obs_names.tolist() == ["fractional", "positive"]
    assert cell_report.summary["parameters"]["min_counts"] == 0.0
    assert cell_report.summary["parameters"]["enable_min_counts"] is True
    assert any("negative expression" in warning for warning in cell_report.warnings)
    namespace = {}
    exec(cell_code, namespace)
    assert namespace["filter_cells"](value).obs_names.equals(cells.obs_names)

    upper_cells, upper_report, _ = filter_cells_owned(
        value,
        max_counts=0.5,
        enable_max_counts=True,
    )
    assert upper_cells.obs_names.tolist() == ["negative", "fractional"]
    assert upper_report.summary["parameters"]["max_counts"] == 0.5

    mitochondrial_cells, mitochondrial_report, mitochondrial_code = filter_cells_owned(
        value,
        max_pct_mito=0.0,
        enable_max_pct_mito=True,
    )
    assert mitochondrial_cells.obs_names.tolist() == ["negative", "positive"]
    assert mitochondrial_report.summary["key_results"]["mitochondrial_metric"]["used_for_filtering"] is True
    namespace = {}
    exec(mitochondrial_code, namespace)
    assert namespace["filter_cells"](value).obs_names.equals(mitochondrial_cells.obs_names)

    genes, gene_report, gene_code = filter_genes_owned(
        value,
        min_counts=0.0,
        enable_min_counts=True,
    )
    assert genes.var_names.tolist() == ["fractional_sum", "positive_sum"]
    assert gene_report.summary["parameters"]["enable_min_counts"] is True
    namespace = {}
    exec(gene_code, namespace)
    assert namespace["filter_genes"](value).var_names.equals(genes.var_names)

    disabled, disabled_report, disabled_code = filter_cells_owned(
        value,
        min_counts=science.np.nan,
        enable_min_counts=False,
    )
    assert disabled.obs_names.equals(value.obs_names)
    assert disabled_report.summary["parameters"]["min_counts"] is None
    json.dumps(disabled_report.summary, allow_nan=False)
    namespace = {}
    exec(disabled_code, namespace)
    assert namespace["filter_cells"](value).obs_names.equals(value.obs_names)

    with pytest.raises(ValueError, match="minimum cannot exceed"):
        filter_genes_owned(
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

    qc, _, _ = calculate_qc_owned(value, percent_top="")
    science.np.testing.assert_array_equal(qc.obs["n_genes_by_counts"], [0, 1])
    science.np.testing.assert_array_equal(qc.var["n_cells_by_counts"], [0, 1])

    cells, _, cells_code = filter_cells_owned(value, min_genes=1)
    assert cells.obs_names.tolist() == ["expressed"]
    namespace = {}
    exec(cells_code, namespace)
    assert namespace["filter_cells"](value).obs_names.tolist() == ["expressed"]

    genes, _, genes_code = filter_genes_owned(value, min_cells=1)
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
    calculate_qc_owned(value.copy(), percent_top="")
    assert_storage_unchanged()
    runtime_cells, _, cells_code = filter_cells_owned(value.copy(), min_genes=1)
    assert_storage_unchanged()
    filter_genes_owned(value.copy(), min_cells=1)
    assert_storage_unchanged()
    _, _, plot_code = qc_plots_owned(value.copy())
    assert_storage_unchanged()

    namespace = {}
    exec(plot_code, namespace)
    namespace["qc_plots"](value.copy())
    assert_storage_unchanged()

    namespace = {}
    exec(cells_code, namespace)
    equivalent_cells = namespace["filter_cells"](value.copy())
    assert_storage_unchanged()
    science.np.testing.assert_array_equal(equivalent_cells.X.data, runtime_cells.X.data)
    science.np.testing.assert_array_equal(equivalent_cells.X.indices, runtime_cells.X.indices)
    science.np.testing.assert_array_equal(equivalent_cells.X.indptr, runtime_cells.X.indptr)
    assert equivalent_cells.X.has_canonical_format is runtime_cells.X.has_canonical_format


def test_qc_filters_reject_contradictory_bounds_and_report_empty_results(adata):
    with pytest.raises(ValueError, match="minimum cannot exceed"):
        filter_cells_owned(adata, min_genes=10, max_genes=5)
    with pytest.raises(ValueError, match="minimum cannot exceed"):
        filter_genes_owned(adata, min_cells=10, max_cells=5)
    empty_cells, cell_report, _ = filter_cells_owned(
        adata,
        min_genes=1_000_000,
    )
    empty_genes, gene_report, _ = filter_genes_owned(
        adata,
        min_cells=1_000_000,
    )
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

    for operation, function_name, keyword, empty_axis in (
        (filter_cells_owned, "filter_cells", {"min_genes": 1}, "obs"),
        (filter_genes_owned, "filter_genes", {"min_cells": 1}, "var"),
    ):
        output, report, code = operation(empty_counts, **keyword)
        assert len(getattr(output, empty_axis + "_names")) == 0
        assert any("removed every" in warning for warning in report.warnings)
        namespace = {}
        exec(code, namespace)
        equivalent = namespace[function_name](empty_counts)
        assert len(getattr(equivalent, empty_axis + "_names")) == 0


def test_qc_plot_reports_missing_mitochondrial_metric_as_unavailable(adata):
    adata.var_names = [f"G{index}" for index in range(adata.n_vars)]

    _, report, code = qc_plots_owned(adata)

    assert "mitochondrial_percent" not in report.summary["key_results"]["distributions"]
    assert "total_expression_vs_mitochondrial_percent" not in report.summary["key_results"]["panels"]
    assert "Mitochondrial metric unavailable" in code


def test_qc_plot_ignores_stale_obs_metrics_and_uses_one_explicit_expression_source(science):
    value = science.ad.AnnData(
        science.sparse.csr_matrix([[1.0, 1.0], [1.0, 1.0]]),
        obs=science.pd.DataFrame(
            {
                "total_counts": [900.0, 900.0],
                "n_genes_by_counts": [2, 2],
                "pct_counts_mt": [99.0, 99.0],
            },
            index=["mitochondrial", "nuclear"],
        ),
        var=science.pd.DataFrame({"mt": [True, False]}, index=["MT-G", "G"]),
    )
    value.layers["counts"] = science.sparse.csr_matrix([[10.0, 0.0], [0.0, 20.0]])
    snapshot = value.copy()

    _, report, code = qc_plots_owned(
        value,
        source={"source": "layer", "source_layer": "counts"},
    )

    key_results = report.summary["key_results"]
    assert key_results["distributions"]["total_expression"]["median"] == 15.0
    assert key_results["distributions"]["detected_genes"]["median"] == 1.0
    assert key_results["distributions"]["mitochondrial_percent"]["median"] == 50.0
    assert set(key_results["metric_sources"].values()) == {"AnnData layer 'counts'"}
    assert_adata_equal(value, snapshot, science)
    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](value)
    primary_offsets = figure.axes[2].collections[0].get_offsets()
    mito_offsets = figure.axes[3].collections[0].get_offsets()
    science.np.testing.assert_allclose(primary_offsets, [[10.0, 1.0], [20.0, 1.0]])
    science.np.testing.assert_allclose(mito_offsets, [[10.0, 100.0], [20.0, 0.0]])


def test_qc_grouped_plot_uses_categorical_groups_and_one_explicit_expression_source(science):
    value = science.ad.AnnData(
        science.np.full((4, 2), 999.0),
        obs=science.pd.DataFrame(
            {
                "sample": science.pd.Categorical(
                    ["S2", "S1", "S2", "S1"],
                    categories=["unused_before", "S1", "S2", "unused_after"],
                    ordered=True,
                ),
                "total_counts": [9_999.0] * 4,
                "pct_counts_mt": [99.0] * 4,
            },
            index=["cell_0", "cell_1", "cell_2", "cell_3"],
        ),
        var=science.pd.DataFrame({"mt": [True, False]}, index=["MT-G", "G"]),
    )
    value.layers["counts"] = science.sparse.csr_matrix(
        [[10.0, 0.0], [0.0, 20.0], [2.0, 3.0], [4.0, 1.0]]
    )
    snapshot = value.copy()

    plot, report, code = qc_plots_owned(
        value,
        view={"view": "grouped", "groupby": "sample"},
        source={"source": "layer", "source_layer": "counts"},
    )

    assert plot.png.startswith(b"\x89PNG\r\n\x1a\n")
    key_results = report.summary["key_results"]
    assert key_results["view"] == "grouped"
    assert key_results["groupby"] == "sample"
    assert key_results["group_order"] == ["S1", "S2"]
    assert key_results["group_cell_counts"] == {"S1": 2, "S2": 2}
    assert key_results["group_distributions"]["S1"]["total_expression"]["median"] == 12.5
    assert key_results["group_distributions"]["S2"]["total_expression"]["median"] == 7.5
    assert key_results["panels"] == [
        "total_expression_by_group",
        "detected_genes_by_group",
        "mitochondrial_percent_by_group",
    ]
    assert set(key_results["metric_sources"].values()) == {"AnnData layer 'counts'"}
    assert report.summary["parameters"]["view"] == {"view": "grouped", "groupby": "sample"}
    json.dumps(report.summary, allow_nan=False)
    assert_adata_equal(value, snapshot, science)

    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](value)
    assert [axis.get_title() for axis in figure.axes] == [
        "Total expression by sample",
        "Genes detected by sample",
        "Mitochondrial expression (%) by sample",
    ]
    assert [tick.get_text() for tick in figure.axes[0].get_xticklabels()] == ["S1", "S2"]
    assert max(path.vertices[:, 1].max() for body in figure.axes[0].collections for path in body.get_paths()) <= 20.0


@pytest.mark.parametrize("dtype", [object, "string"])
def test_qc_grouped_plot_orders_string_groups_by_first_appearance(dtype, science):
    value = science.ad.AnnData(
        science.np.eye(4),
        obs=science.pd.DataFrame(
            {"sample": science.pd.Series(["S2", "S1", "S2", "S3"], dtype=dtype).array},
            index=[f"cell_{index}" for index in range(4)],
        ),
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(4)]),
    )

    _, report, code = qc_plots_owned(value, view={"view": "grouped", "groupby": "sample"})

    assert report.summary["key_results"]["group_order"] == ["S2", "S1", "S3"]
    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](value)
    assert [tick.get_text() for tick in figure.axes[0].get_xticklabels()] == ["S2", "S1", "S3"]


def test_qc_grouped_plot_accepts_boolean_groups(science):
    value = science.ad.AnnData(
        science.np.eye(3),
        obs=science.pd.DataFrame({"selected": [True, False, True]}, index=[f"cell_{index}" for index in range(3)]),
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(3)]),
    )

    _, report, code = qc_plots_owned(value, view={"view": "grouped", "groupby": "selected"})

    assert report.summary["key_results"]["group_order"] == ["True", "False"]
    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](value)
    assert [tick.get_text() for tick in figure.axes[0].get_xticklabels()] == ["True", "False"]


@pytest.mark.parametrize(
    ("case", "error", "message"),
    [
        ("missing column", ValueError, "groupby column not found"),
        ("continuous numeric", TypeError, "must be categorical, boolean, or string-like"),
        ("missing categorical label", ValueError, "contains missing group labels"),
        ("missing object label", ValueError, "contains missing group labels"),
        ("blank label", ValueError, "contains blank group labels"),
        ("display collision", ValueError, "collide after string conversion"),
        ("non-finite object label", ValueError, "contains non-finite group labels"),
    ],
)
def test_qc_grouped_plot_validates_discrete_group_labels_in_runtime_and_code(case, error, message, science):
    if case == "missing column":
        obs = {"other": ["S1", "S2"]}
    elif case == "continuous numeric":
        obs = {"sample": [1.0, 2.0]}
    elif case == "missing categorical label":
        obs = {"sample": science.pd.Categorical(["S1", None], categories=["S1", "S2"])}
    elif case == "missing object label":
        obs = {"sample": science.pd.Series(["S1", None], dtype=object).array}
    elif case == "blank label":
        obs = {"sample": ["S1", " "]}
    elif case == "display collision":
        obs = {"sample": science.pd.Series([1, "1"], dtype=object).array}
    else:
        obs = {"sample": science.pd.Series([1, science.np.inf], dtype=object).array}
    value = science.ad.AnnData(
        science.np.eye(2),
        obs=science.pd.DataFrame(obs, index=["cell_0", "cell_1"]),
        var=science.pd.DataFrame(index=["gene_0", "gene_1"]),
    )
    valid = value.copy()
    valid.obs = science.pd.DataFrame({"sample": ["S1", "S2"]}, index=valid.obs_names)
    _, _, code = qc_plots_owned(valid, view={"view": "grouped", "groupby": "sample"})

    with pytest.raises(error, match=message):
        qc_plots_owned(value, view={"view": "grouped", "groupby": "sample"})
    namespace = {}
    exec(code, namespace)
    with pytest.raises(error, match=message):
        namespace["qc_plots"](value)


def test_qc_grouped_plot_caps_visible_groups_in_runtime_and_code(science):
    groups = [f"S{index:02}" for index in range(41)]
    value = science.ad.AnnData(
        science.np.eye(len(groups)),
        obs=science.pd.DataFrame({"sample": groups}, index=[f"cell_{index}" for index in range(len(groups))]),
        var=science.pd.DataFrame(index=[f"gene_{index}" for index in range(len(groups))]),
    )
    valid = value[:2].copy()
    _, _, code = qc_plots_owned(valid, view={"view": "grouped", "groupby": "sample"})

    with pytest.raises(ValueError, match="supports at most 40 visible groups"):
        qc_plots_owned(value, view={"view": "grouped", "groupby": "sample"})
    namespace = {}
    exec(code, namespace)
    with pytest.raises(ValueError, match="supports at most 40 visible groups"):
        namespace["qc_plots"](value)


def test_qc_plot_ignores_nonfinite_obs_mitochondrial_cache(adata, science):
    qc = calculate_qc_owned(adata)[0]
    qc.obs["pct_counts_mt"] = science.np.nan

    _, report, code = qc_plots_owned(qc)

    key_results = report.summary["key_results"]
    assert key_results["distributions"]["mitochondrial_percent"]["missing"] == 0
    assert "total_expression_vs_mitochondrial_percent" in key_results["panels"]
    assert set(key_results["metric_sources"].values()) == {"AnnData.X"}
    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](qc)
    assert figure.axes[3].axison is True


def test_qc_plot_reports_zero_library_mitochondrial_percentage_as_missing(science):
    value = science.ad.AnnData(
        science.sparse.csr_matrix([[0.0, 0.0], [0.0, 10.0]]),
        obs=science.pd.DataFrame(index=["empty", "expressed"]),
        var=science.pd.DataFrame({"mt": [True, False]}, index=["MT-G", "G"]),
    )

    _, report, code = qc_plots_owned(value)

    distribution = report.summary["key_results"]["distributions"]["mitochondrial_percent"]
    assert distribution["n"] == value.n_obs
    assert distribution["missing"] == 1
    assert any("non-finite mitochondrial percentages" in warning for warning in report.warnings)
    namespace = {}
    exec(code, namespace)
    figure = namespace["qc_plots"](value)
    science.np.testing.assert_allclose(figure.axes[3].collections[0].get_offsets(), [[10.0, 0.0]])


def test_qc_plot_code_reproduces_nonfinite_selected_source_error(science):
    value = science.ad.AnnData(science.np.asarray([[1.0, 0.0], [0.0, 1.0]]))
    _, _, code = qc_plots_owned(value)
    invalid = value.copy()
    invalid.X[0, 0] = science.np.nan

    with pytest.raises(ValueError, match="non-finite expression values"):
        qc_plots_owned(invalid)
    namespace = {}
    exec(code, namespace)
    with pytest.raises(ValueError, match="non-finite expression values"):
        namespace["qc_plots"](invalid)


def test_qc_equivalent_code_reproduces_primary_results(adata, science):
    qc, _, qc_code = calculate_qc_owned(adata, mitochondrial_prefix="mt-")
    namespace = {}
    exec(qc_code, namespace)
    equivalent_qc = namespace["calculate_qc_metrics"](adata)
    science.pd.testing.assert_frame_equal(qc.obs, equivalent_qc.obs)
    science.pd.testing.assert_frame_equal(qc.var, equivalent_qc.var)

    filtered_cells, _, cells_code = filter_cells_owned(qc, min_genes=3)
    namespace = {}
    exec(cells_code, namespace)
    equivalent_cells = namespace["filter_cells"](qc)
    assert equivalent_cells.obs_names.equals(filtered_cells.obs_names)

    filtered_genes, _, genes_code = filter_genes_owned(filtered_cells, min_cells=2)
    namespace = {}
    exec(genes_code, namespace)
    equivalent_genes = namespace["filter_genes"](filtered_cells)
    assert equivalent_genes.var_names.equals(filtered_genes.var_names)

    _, _, plot_code = qc_plots_owned(filtered_genes)
    namespace = {}
    exec(plot_code, namespace)
    equivalent_figure = namespace["qc_plots"](filtered_genes)
    assert len(equivalent_figure.axes) == 4
    assert equivalent_figure.axes[0].get_title() == "Total expression per cell"
