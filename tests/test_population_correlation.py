from __future__ import annotations

import copy
import json

import pytest

from openbio_singlecell.node_types import PlotResultType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_population import OpenBioSingleCellCellTypeCorrelation
from openbio_singlecell.population_correlation import _strongest_signed_pairs, analyze_population_centroid_correlation


@pytest.fixture
def correlation_adata(science):
    labels = science.pd.Categorical(
        ["A", "A", "A", "B", "B", "B", "C", "C", "C"],
        categories=["B", "A", "C", "unused"],
        ordered=True,
    )
    obs = science.pd.DataFrame(
        {"cell_type": labels, "sample": [f"s{index // 3}" for index in range(9)]},
        index=[f"c{index}" for index in range(9)],
    )
    var = science.pd.DataFrame(index=["g1", "g2"])
    adata = science.ad.AnnData(science.np.arange(18, dtype=float).reshape(9, 2), obs=obs, var=var)
    adata.obsm["X_pca"] = science.np.asarray(
        [
            [0.0, 1.0, 2.0, 3.0],
            [0.1, 1.2, 2.1, 3.2],
            [-0.1, 0.9, 1.8, 2.9],
            [0.0, 1.0, 2.5, 4.0],
            [0.2, 1.1, 2.7, 4.1],
            [-0.1, 0.8, 2.4, 3.8],
            [3.0, 2.2, 1.0, -0.2],
            [3.2, 2.0, 0.9, 0.0],
            [2.8, 2.1, 1.1, -0.1],
        ],
        dtype=float,
    )
    adata.uns["sentinel"] = {"nested": [1, 2, 3]}
    return adata


def _run(adata, **overrides):
    parameters = {
        "population_key": "cell_type",
        "representation_key": "X_pca",
        "correlation_method": "pearson",
        "linkage_method": "complete",
        "n_dimensions": 0,
        "annotation_status": "provisional",
        "color_map": "RdYlBu",
        "show_numbers": False,
        "max_groups": 200,
        "max_output_rows": 100000,
    }
    parameters.update(overrides)
    return analyze_population_centroid_correlation(adata, **parameters)


def test_population_correlation_schema_is_atomic():
    schema = OpenBioSingleCellCellTypeCorrelation.GET_SCHEMA()
    assert schema.display_name == "Population Centroid Correlation"
    assert [input_.id for input_ in schema.inputs] == [
        "adata",
        "population_key",
        "representation_key",
        "correlation_method",
        "linkage_method",
        "n_dimensions",
        "annotation_status",
        "color_map",
        "show_numbers",
        "max_groups",
        "max_output_rows",
    ]
    assert [(output.display_name, output.io_type) for output in schema.outputs] == [
        ("table", TableResultType.io_type),
        ("plot", PlotResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]


def test_population_correlation_matches_independent_centroids_and_is_immutable(correlation_adata, science):
    before = correlation_adata.copy()
    table, png, summary = _run(correlation_adata)

    assert table.shape == (3, 10)
    assert set(table["population_a"]) | set(table["population_b"]) == {"A", "B", "C"}
    assert table["cells_a"].tolist() == [3, 3, 3]
    assert table["cells_b"].tolist() == [3, 3, 3]
    labels = correlation_adata.obs["cell_type"]
    categories = ["B", "A", "C"]
    centroids = science.np.vstack(
        [correlation_adata.obsm["X_pca"][science.np.asarray(labels == category)].mean(axis=0) for category in categories]
    )
    expected = science.pd.DataFrame(centroids, index=categories).T.corr(method="pearson")
    for row in table.itertuples(index=False):
        assert row.correlation == pytest.approx(expected.loc[row.population_a, row.population_b], abs=1e-10)
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert summary["status"] == "descriptive_population_similarity"
    assert summary["key_results"]["declared_observed_order"] == ["B", "A", "C"]
    assert summary["key_results"]["unused_declared_categories"] == ["unused"]
    assert summary["key_results"]["population_cell_counts"] == {"B": 3, "A": 3, "C": 3}
    assert set(summary["key_results"]["dendrogram_order"]) == {"A", "B", "C"}
    assert any("provisional" in warning for warning in summary["warnings"])
    json.dumps(summary, allow_nan=False)

    science.np.testing.assert_array_equal(correlation_adata.X, before.X)
    science.np.testing.assert_array_equal(correlation_adata.obsm["X_pca"], before.obsm["X_pca"])
    science.pd.testing.assert_frame_equal(correlation_adata.obs, before.obs)
    assert correlation_adata.uns == before.uns


def test_population_correlation_node_and_generated_code_are_equivalent(correlation_adata, science):
    table_result, plot_result, report, code = OpenBioSingleCellCellTypeCorrelation.execute(
        correlation_adata,
        population_key="cell_type",
        representation_key="X_pca",
        correlation_method="spearman",
        linkage_method="average",
        n_dimensions=3,
        annotation_status="curated",
        color_map="viridis",
        show_numbers=True,
        max_groups=10,
        max_output_rows=100,
    ).result
    assert table_result.source["operation"] == "population_centroid_correlation"
    assert plot_result.png.startswith(b"\x89PNG")
    assert report.summary["parameters"]["resolved_dimensions"] == 3
    assert "from openbio_singlecell" not in code
    namespace = {}
    exec(code, namespace)
    reproduced_table, reproduced_png, reproduced_summary = namespace["population_centroid_correlation"](
        correlation_adata
    )
    science.pd.testing.assert_frame_equal(reproduced_table, table_result.table)
    assert reproduced_png == plot_result.png
    assert reproduced_summary == report.summary


@pytest.mark.parametrize("method", ["pearson", "spearman", "kendall"])
def test_population_correlation_supports_official_methods(correlation_adata, method):
    table, _, summary = _run(correlation_adata, correlation_method=method)
    assert table["correlation_method"].unique().tolist() == [method]
    assert summary["parameters"]["correlation_method"] == method


def _pair_table(science, correlations):
    return science.pd.DataFrame(
        {
            "population_a": [f"a{index}" for index in range(len(correlations))],
            "population_b": [f"b{index}" for index in range(len(correlations))],
            "correlation": correlations,
            "dendrogram_order_a": list(range(len(correlations))),
            "dendrogram_order_b": list(range(len(correlations))),
        }
    )


def test_population_correlation_signed_extrema_do_not_mislabel_absent_signs(science):
    positive = _pair_table(science, [0.9, 0.3, 0.0])
    negative = _pair_table(science, [-0.9, -0.3, 0.0])

    assert _strongest_signed_pairs(positive, direction="negative") == []
    assert _strongest_signed_pairs(negative, direction="positive") == []
    assert [row["correlation"] for row in _strongest_signed_pairs(positive, direction="positive")] == [0.9, 0.3]
    assert [row["correlation"] for row in _strongest_signed_pairs(negative, direction="negative")] == [-0.9, -0.3]


@pytest.mark.parametrize(
    ("direction", "correlations", "expected_count", "expected_cutoff"),
    [
        ("positive", [1.0 - index / 20 for index in range(9)] + [0.5, 0.5, 0.5, 0.1], 12, 0.5),
        ("negative", [-1.0 + index / 20 for index in range(9)] + [-0.5, -0.5, -0.5, -0.1], 12, -0.5),
    ],
)
def test_population_correlation_signed_extrema_retain_cutoff_ties(
    science, direction, correlations, expected_count, expected_cutoff
):
    result = _strongest_signed_pairs(_pair_table(science, correlations), direction=direction)

    assert len(result) == expected_count
    assert result[-1]["correlation"] == expected_cutoff


def _replace_population_labels(adata, labels):
    obs = adata.obs.copy()
    obs["cell_type"] = labels
    adata.obs = obs


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda adata: adata.obsm.__setitem__("X_pca", adata.obsm["X_pca"][:, :1]), "dense numeric"),
        (lambda adata: adata.obsm["X_pca"].__setitem__((0, 0), float("nan")), "NaN or infinity"),
        (lambda adata: _replace_population_labels(adata, ["A"] * adata.n_obs), "at least two"),
        (
            lambda adata: _replace_population_labels(
                adata,
                [" A", "A", "A", "B", "B", "B", "C", "C", "C"],
            ),
            "canonical",
        ),
    ],
)
def test_population_correlation_rejects_invalid_inputs(correlation_adata, mutate, message):
    mutate(correlation_adata)
    with pytest.raises((TypeError, ValueError), match=message):
        _run(correlation_adata)


def test_population_correlation_rejects_sparse_representation(correlation_adata, science):
    correlation_adata.obsm["X_pca"] = science.sparse.csr_matrix(correlation_adata.obsm["X_pca"])
    with pytest.raises(TypeError, match="dense"):
        _run(correlation_adata)


def test_population_correlation_rejects_complex_representation(correlation_adata):
    correlation_adata.obsm["X_pca"] = correlation_adata.obsm["X_pca"].astype(complex)
    with pytest.raises(ValueError, match="dense numeric"):
        _run(correlation_adata)


def test_population_correlation_rejects_undefined_constant_centroids(correlation_adata, science):
    correlation_adata.obsm["X_pca"] = science.np.tile(science.np.asarray([1.0, 1.0, 1.0, 1.0]), (9, 1))
    with pytest.raises(ValueError, match="nonconstant dimensions"):
        _run(correlation_adata)


def test_population_correlation_guards_and_visual_parameters(correlation_adata):
    with pytest.raises(ValueError, match="max_groups"):
        _run(correlation_adata, max_groups=2)
    with pytest.raises(ValueError, match="max_output_rows"):
        _run(correlation_adata, max_output_rows=2)
    with pytest.raises(ValueError, match="Unknown Matplotlib colormap"):
        _run(correlation_adata, color_map="not-a-colormap")
    with pytest.raises(ValueError, match="n_dimensions"):
        _run(correlation_adata, n_dimensions=1)


def test_population_correlation_rejects_malformed_scanpy_result(correlation_adata, science, monkeypatch):
    original = science.sc.tl.dendrogram

    def malformed(*args, **kwargs):
        result = copy.deepcopy(original(*args, **kwargs))
        result["correlation_matrix"][0, 1] = 0.123
        return result

    monkeypatch.setattr(science.sc.tl, "dendrogram", malformed)
    with pytest.raises(RuntimeError, match="not symmetric"):
        _run(correlation_adata)
