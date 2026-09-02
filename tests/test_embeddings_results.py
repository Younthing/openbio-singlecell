from __future__ import annotations

import copy

import pytest

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.marker_evidence import MARKER_COLUMNS
from openbio_singlecell.operations_embedding import leiden, neighbors, pca, umap
from openbio_singlecell.operations_preprocess import highly_variable_genes, normalize_to_layer
from openbio_singlecell.operations_results import embedding_plot_owned, marker_genes_owned
from tests.artifact_operation_harness import run_anndata_operation


def output_value(node_output):
    return node_output.result[0] if hasattr(node_output, "result") else node_output[0]


def dense(matrix, science):
    return matrix.toarray() if science.sparse.issparse(matrix) else science.np.asarray(matrix)


@pytest.fixture
def adata(science):
    rng = science.np.random.default_rng(4)
    counts = rng.poisson(2.0, size=(30, 12)).astype(float)
    counts[:15, :3] += 5
    counts[15:, 3:6] += 5
    obs = science.pd.DataFrame(
        {"group": science.pd.Categorical(["A"] * 15 + ["B"] * 15)},
        index=[f"cell_{index}" for index in range(30)],
    )
    var = science.pd.DataFrame(index=["MT-G0", *[f"G{index}" for index in range(1, 12)]])
    value = science.ad.AnnData(science.sparse.csr_matrix(counts), obs=obs, var=var)
    value.layers["counts"] = value.X.copy()
    value.layers["alternate"] = value.X.copy()
    value.raw = value.copy()
    ensure_metadata(value, display_name="synthetic", source={"kind": "test"})
    return value


def prepare_embedding_input(adata):
    logged = output_value(
        run_anndata_operation(
            normalize_to_layer,
            adata,
            {
                "source": {"source": "layer", "source_layer": "counts"},
                "target_sum": 10_000.0,
                "transform": "log1p",
                "output_layer": "log1p_norm",
                "overwrite_existing": False,
            },
        )
    )
    variable = output_value(
        run_anndata_operation(
            highly_variable_genes,
            logged,
            {
                "n_top_genes": 8,
                "flavor": "seurat",
                "source": {"source": "layer", "layer_name": "log1p_norm"},
                "batch_key": "",
                "always_keep_genes": "",
                "subset": False,
                "theta": 100.0,
                "clipping_mode": "sqrt_n_obs",
                "custom_clip": 10.0,
                "chunksize": 1_000,
                "span": 0.3,
                "n_bins": 20,
                "overwrite_existing": False,
            },
        )
    )
    return variable


def test_dimension_reduction_chain_publishes_distinct_artifacts(adata):
    scaled = prepare_embedding_input(adata)
    scaled_history = copy.deepcopy(scaled.uns["openbio_singlecell"]["analysis_history"])

    pca_output = run_anndata_operation(
        pca,
        scaled,
        {
            "n_comps": 3,
            "use_hvg": True,
            "source": {"source": "layer", "layer_name": "log1p_norm"},
            "overwrite_existing": False,
            "max_output_gib": 2.0,
            "random_seed": 0,
        },
    )
    pca_result = output_value(pca_output)
    assert "X_pca" not in scaled.obsm
    assert repr(scaled.uns["openbio_singlecell"]["analysis_history"]) == repr(scaled_history)

    neighbors_result = output_value(
        run_anndata_operation(
            neighbors,
            pca_result,
            {
                "use_rep": "X_pca",
                "n_dimensions": 3,
                "n_neighbors": 5,
                "metric": "cosine",
                "method": "umap",
                "key_added": "neighbors",
                "overwrite_existing": False,
                "random_seed": 0,
            },
        )
    )
    assert "connectivities" not in pca_result.obsp

    umap_result = output_value(
        run_anndata_operation(
            umap,
            neighbors_result,
            {
                "neighbors_key": "neighbors",
                "min_dist": 0.5,
                "spread": 1.0,
                "key_added": "X_umap",
                "overwrite_existing": False,
                "random_seed": 0,
            },
        )
    )
    assert "X_umap" not in neighbors_result.obsm

    leiden_result = output_value(
        run_anndata_operation(
            leiden,
            umap_result,
            {
                "resolution": 1.0,
                "key_added": "leiden",
                "neighbors_key": "neighbors",
                "n_iterations": 2,
                "stability_repeats": 5,
                "overwrite_existing": False,
                "random_seed": 0,
            },
        )
    )
    assert "leiden" not in umap_result.obs

    assert pca_result.obsm["X_pca"].shape == (adata.n_obs, 3)
    assert "connectivities" in neighbors_result.obsp
    assert umap_result.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert "leiden" in leiden_result.obs
    assert len(leiden_result.uns["openbio_singlecell"]["analysis_history"]) == 6


def test_embedding_preconditions_are_clear(adata):
    with pytest.raises(ValueError, match="Highly Variable Genes first"):
        run_anndata_operation(
            pca,
            adata,
            {
                "n_comps": 3,
                "use_hvg": True,
                "source": {"source": "layer", "layer_name": "counts"},
                "overwrite_existing": False,
                "max_output_gib": 2.0,
                "random_seed": 0,
            },
        )

    with pytest.raises(ValueError, match="cannot be empty"):
        run_anndata_operation(
            leiden,
            adata,
            {
                "resolution": 1.0,
                "key_added": "",
                "neighbors_key": "neighbors",
                "n_iterations": 2,
                "stability_repeats": 5,
                "overwrite_existing": False,
                "random_seed": 0,
            },
        )


def test_marker_layer_integration_has_stable_evidence_and_universe(adata):
    logged = prepare_embedding_input(adata)
    table, universe, summary, code = marker_genes_owned(
        logged,
        "group",
        "wilcoxon",
        {"source": "layer", "layer_name": "log1p_norm"},
        3,
        True,
        1_000,
    )

    assert table.kind == universe.kind == "table"
    assert summary.kind == "summary"
    assert list(table.table.columns) == MARKER_COLUMNS
    assert universe.table["gene"].tolist() == logged.var_names.tolist()
    assert set(table.table["group"]) == {"A", "B"}
    assert table.table["rank"].min() == 1
    assert table.parameters["analysis_fingerprint"] == universe.parameters["analysis_fingerprint"]
    compile(code, "<marker-code>", "exec")


def test_marker_preconditions_are_clear(adata):
    logged = prepare_embedding_input(adata)
    with pytest.raises(ValueError, match="groupby column not found"):
        marker_genes_owned(
            logged,
            "missing",
            "wilcoxon",
            {"source": "layer", "layer_name": "log1p_norm"},
            3,
            True,
            1_000,
        )
    with pytest.raises(ValueError, match="Marker source layer not found"):
        marker_genes_owned(
            logged,
            "group",
            "wilcoxon",
            {"source": "layer", "layer_name": "missing"},
            3,
            True,
            1_000,
        )
    with pytest.raises(ValueError, match="Unsupported marker source: 'raw'"):
        marker_genes_owned(logged, "group", "wilcoxon", {"source": "raw"}, 3, True, 1_000)
    with pytest.raises(ValueError, match="Unsupported Marker Genes method"):
        marker_genes_owned(
            logged,
            "group",
            "logreg",
            {"source": "layer", "layer_name": "log1p_norm"},
            3,
            False,
            1_000,
        )


def test_embedding_plot_is_read_only_for_categorical_and_numeric_colors(adata, science):
    adata.obsm["X_umap"] = science.np.arange(adata.n_obs * 2, dtype=float).reshape(adata.n_obs, 2)
    adata.obs["score"] = science.np.linspace(0.0, 1.0, adata.n_obs)
    snapshot_x = dense(adata.X, science).copy()
    snapshot_obs = adata.obs.copy(deep=True)
    snapshot_uns = copy.deepcopy(adata.uns)

    categorical = output_value(
        embedding_plot_owned(
            adata,
            embedding_key="X_umap",
            color={"color": "obs", "obs_key": "group", "color_mode": "auto"},
            point_size=8.0,
            continuous_color_map="viridis",
        )
    )
    numeric = output_value(
        embedding_plot_owned(
            adata,
            embedding_key="X_umap",
            color={"color": "obs", "obs_key": "score", "color_mode": "auto"},
            point_size=8.0,
            continuous_color_map="viridis",
        )
    )

    assert categorical.kind == "plot" and categorical.png
    assert numeric.kind == "plot" and numeric.png
    science.np.testing.assert_array_equal(dense(adata.X, science), snapshot_x)
    science.pd.testing.assert_frame_equal(adata.obs, snapshot_obs)
    assert adata.uns == snapshot_uns


def test_embedding_plot_preconditions_are_clear(adata, science):
    with pytest.raises(ValueError, match="coordinates not found"):
        embedding_plot_owned(
            adata,
            embedding_key="X_umap",
            color={"color": "obs", "obs_key": "group", "color_mode": "auto"},
            point_size=8.0,
            continuous_color_map="viridis",
        )

    adata.obsm["X_umap"] = science.np.ones((adata.n_obs, 1))
    with pytest.raises(ValueError, match="shape exactly"):
        embedding_plot_owned(
            adata,
            embedding_key="X_umap",
            color={"color": "obs", "obs_key": "group", "color_mode": "auto"},
            point_size=8.0,
            continuous_color_map="viridis",
        )

    adata.obsm["X_umap"] = science.np.ones((adata.n_obs, 2))
    with pytest.raises(ValueError, match="color column not found"):
        embedding_plot_owned(
            adata,
            embedding_key="X_umap",
            color={"color": "obs", "obs_key": "missing", "color_mode": "auto"},
            point_size=8.0,
            continuous_color_map="viridis",
        )
