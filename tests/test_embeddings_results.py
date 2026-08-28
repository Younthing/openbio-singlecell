from __future__ import annotations

import copy

import pytest

from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.nodes_embedding import (
    OpenBioSingleCellLeiden,
    OpenBioSingleCellNeighbors,
    OpenBioSingleCellPCA,
    OpenBioSingleCellUMAP,
)
from openbio_singlecell.nodes_preprocess import (
    OpenBioSingleCellHighlyVariableGenes,
    OpenBioSingleCellNormalizeToLayer,
)
from openbio_singlecell.nodes_results import (
    MARKER_COLUMNS,
    OpenBioSingleCellMarkerGenes,
    OpenBioSingleCellUMAPPlot,
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
        OpenBioSingleCellNormalizeToLayer.execute(
            adata,
            source={"source": "layer", "source_layer": "counts"},
            target_sum=10_000.0,
            transform="log1p",
            output_layer="log1p_norm",
        )
    )
    variable = output_value(
        OpenBioSingleCellHighlyVariableGenes.execute(
            logged,
            n_top_genes=8,
            flavor="seurat",
            source={"source": "layer", "layer_name": "log1p_norm"},
            subset=False,
        )
    )
    return variable


def test_dimension_reduction_chain_is_copy_on_write(adata):
    scaled = prepare_embedding_input(adata)
    scaled_history = copy.deepcopy(scaled.uns["openbio_singlecell"]["analysis_history"])

    pca = output_value(OpenBioSingleCellPCA.execute(scaled, n_comps=3, use_hvg=True, random_seed=0))
    assert "X_pca" not in scaled.obsm
    assert scaled.uns["openbio_singlecell"]["analysis_history"] == scaled_history

    neighbors = output_value(
        OpenBioSingleCellNeighbors.execute(pca, n_neighbors=5, n_dimensions=3, metric="cosine", random_seed=0)
    )
    assert "connectivities" not in pca.obsp

    umap = output_value(OpenBioSingleCellUMAP.execute(neighbors, min_dist=0.5, spread=1.0, random_seed=0))
    assert "X_umap" not in neighbors.obsm

    leiden = output_value(OpenBioSingleCellLeiden.execute(umap, resolution=1.0, key_added="leiden", random_seed=0))
    assert "leiden" not in umap.obs

    assert pca.obsm["X_pca"].shape == (adata.n_obs, 3)
    assert "connectivities" in neighbors.obsp
    assert umap.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert "leiden" in leiden.obs
    assert len(leiden.uns["openbio_singlecell"]["analysis_history"]) == 6


def test_embedding_preconditions_are_clear(adata):
    with pytest.raises(ValueError, match="Highly Variable Genes first"):
        OpenBioSingleCellPCA.execute(adata, n_comps=3, use_hvg=True, random_seed=0)

    with pytest.raises(ValueError, match="cannot be empty"):
        OpenBioSingleCellLeiden.execute(adata, resolution=1.0, key_added="", random_seed=0)


def test_marker_layer_integration_has_stable_evidence_and_universe(adata):
    logged = prepare_embedding_input(adata)
    table, universe, summary, code = OpenBioSingleCellMarkerGenes.execute(
        logged,
        "group",
        "wilcoxon",
        {"source": "layer", "layer_name": "log1p_norm"},
        3,
        True,
        1_000,
    ).result

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
        OpenBioSingleCellMarkerGenes.execute(
            logged,
            "missing",
            "wilcoxon",
            {"source": "layer", "layer_name": "log1p_norm"},
            3,
            True,
            1_000,
        )
    with pytest.raises(ValueError, match="Marker source layer not found"):
        OpenBioSingleCellMarkerGenes.execute(
            logged,
            "group",
            "wilcoxon",
            {"source": "layer", "layer_name": "missing"},
            3,
            True,
            1_000,
        )
    with pytest.raises(ValueError, match="Unsupported marker source: 'raw'"):
        OpenBioSingleCellMarkerGenes.execute(logged, "group", "wilcoxon", {"source": "raw"}, 3, True, 1_000)
    with pytest.raises(ValueError, match="Unsupported Marker Genes method"):
        OpenBioSingleCellMarkerGenes.execute(
            logged,
            "group",
            "logreg",
            {"source": "layer", "layer_name": "log1p_norm"},
            3,
            False,
            1_000,
        )


def test_umap_plot_is_read_only_for_categorical_and_numeric_colors(adata, science):
    adata.obsm["X_umap"] = science.np.arange(adata.n_obs * 2, dtype=float).reshape(adata.n_obs, 2)
    adata.obs["score"] = science.np.linspace(0.0, 1.0, adata.n_obs)
    snapshot_x = dense(adata.X, science).copy()
    snapshot_obs = adata.obs.copy(deep=True)
    snapshot_uns = copy.deepcopy(adata.uns)

    categorical = output_value(
        OpenBioSingleCellUMAPPlot.execute(adata, "X_umap", "group", "auto", 8.0, "viridis")
    )
    numeric = output_value(OpenBioSingleCellUMAPPlot.execute(adata, "X_umap", "score", "auto", 8.0, "viridis"))

    assert categorical.kind == "plot" and categorical.png
    assert numeric.kind == "plot" and numeric.png
    science.np.testing.assert_array_equal(dense(adata.X, science), snapshot_x)
    science.pd.testing.assert_frame_equal(adata.obs, snapshot_obs)
    assert adata.uns == snapshot_uns


def test_umap_plot_preconditions_are_clear(adata, science):
    with pytest.raises(ValueError, match="coordinates not found"):
        OpenBioSingleCellUMAPPlot.execute(adata, "X_umap", "group", "auto", 8.0, "viridis")

    adata.obsm["X_umap"] = science.np.ones((adata.n_obs, 1))
    with pytest.raises(ValueError, match="shape exactly"):
        OpenBioSingleCellUMAPPlot.execute(adata, "X_umap", "group", "auto", 8.0, "viridis")

    adata.obsm["X_umap"] = science.np.ones((adata.n_obs, 2))
    with pytest.raises(ValueError, match="color column not found"):
        OpenBioSingleCellUMAPPlot.execute(adata, "X_umap", "missing", "auto", 8.0, "viridis")
