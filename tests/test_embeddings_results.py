from __future__ import annotations

import copy
import json
import warnings
from pathlib import Path

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
    OpenBioSingleCellLog1p,
    OpenBioSingleCellNormalizeTotal,
    OpenBioSingleCellScale,
)
from openbio_singlecell.nodes_results import (
    MARKER_COLUMNS,
    OpenBioSingleCellMarkerGenes,
    OpenBioSingleCellUMAPPlot,
)

PLUGIN_ROOT = Path(__file__).resolve().parents[1]


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
    normalized = output_value(OpenBioSingleCellNormalizeTotal.execute(adata, 10_000.0))
    logged = output_value(OpenBioSingleCellLog1p.execute(normalized, False))
    variable = output_value(OpenBioSingleCellHighlyVariableGenes.execute(logged, 8, "seurat", False))
    with pytest.warns(UserWarning, match="densifies"):
        return output_value(OpenBioSingleCellScale.execute(variable, 10.0))


def full_example_node(node_type):
    workflow_path = PLUGIN_ROOT / "example_workflows" / "openbio_singlecell_full_analysis.json"
    workflow = json.loads(workflow_path.read_text(encoding="utf-8"))
    return next(node for node in workflow["nodes"] if node["type"] == node_type)


def test_dimension_reduction_chain_is_copy_on_write(adata):
    scaled = prepare_embedding_input(adata)
    scaled_history = copy.deepcopy(scaled.uns["openbio_singlecell"]["analysis_history"])

    pca = output_value(OpenBioSingleCellPCA.execute(scaled, n_comps=3, use_hvg=True, random_seed=0))
    assert "X_pca" not in scaled.obsm
    assert scaled.uns["openbio_singlecell"]["analysis_history"] == scaled_history

    neighbors = output_value(
        OpenBioSingleCellNeighbors.execute(pca, n_neighbors=5, n_pcs=3, metric="cosine", random_seed=0)
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
    assert len(leiden.uns["openbio_singlecell"]["analysis_history"]) == 8


def test_embedding_preconditions_are_clear(adata):
    with pytest.raises(ValueError, match="Highly Variable Genes first"):
        OpenBioSingleCellPCA.execute(adata, n_comps=3, use_hvg=True, random_seed=0)

    with pytest.raises(ValueError, match="cannot be empty"):
        OpenBioSingleCellLeiden.execute(adata, resolution=1.0, key_added="", random_seed=0)


@pytest.mark.parametrize(
    "source,layer_name",
    [("X", ""), ("raw", ""), ("layer", "alternate")],
)
def test_marker_sources_have_stable_columns(adata, source, layer_name):
    result = output_value(
        OpenBioSingleCellMarkerGenes.execute(adata, "group", "wilcoxon", source, layer_name, 3, True, 0)
    )

    assert result.kind == "table"
    assert list(result.table.columns) == MARKER_COLUMNS
    assert set(result.table["group"]) == {"A", "B"}
    assert result.table["rank"].min() == 1


def test_marker_preconditions_are_clear(adata):
    with pytest.raises(ValueError, match="groupby column not found"):
        OpenBioSingleCellMarkerGenes.execute(adata, "missing", "wilcoxon", "X", "", 3, True, 0)
    with pytest.raises(ValueError, match="Marker layer not found"):
        OpenBioSingleCellMarkerGenes.execute(adata, "group", "wilcoxon", "layer", "missing", 3, True, 0)

    without_raw = adata.copy()
    without_raw.raw = None
    with pytest.raises(ValueError, match="adata.raw is unavailable"):
        OpenBioSingleCellMarkerGenes.execute(without_raw, "group", "wilcoxon", "raw", "", 3, True, 0)


def test_marker_raw_provenance_uses_raw_gene_count(adata):
    subset = adata[:, :8].copy()
    result = output_value(OpenBioSingleCellMarkerGenes.execute(subset, "group", "wilcoxon", "raw", "", 3, True, 0))

    assert subset.n_vars == 8
    assert subset.raw.n_vars == adata.n_vars
    assert result.input_genes == adata.n_vars
    assert result.source["input_genes"] == adata.n_vars


def test_full_example_marker_logfc_is_finite_after_scaling_without_runtime_warnings(adata, science):
    log1p_values = full_example_node("OpenBioSingleCellLog1p")["widgets_values"]
    variable_values = full_example_node("OpenBioSingleCellHighlyVariableGenes")["widgets_values"]
    scale_values = full_example_node("OpenBioSingleCellScale")["widgets_values"]
    marker_values = full_example_node("OpenBioSingleCellMarkerGenes")["widgets_values"]

    normalized = output_value(OpenBioSingleCellNormalizeTotal.execute(adata, 10_000.0))
    logged = output_value(OpenBioSingleCellLog1p.execute(normalized, *log1p_values))
    variable = output_value(OpenBioSingleCellHighlyVariableGenes.execute(logged, *variable_values))
    with pytest.warns(UserWarning, match="densifies"):
        scaled = output_value(OpenBioSingleCellScale.execute(variable, *scale_values))
    scaled.obs["leiden"] = scaled.obs["group"].copy()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = output_value(OpenBioSingleCellMarkerGenes.execute(scaled, *marker_values))

    runtime_warnings = [warning for warning in caught if issubclass(warning.category, RuntimeWarning)]
    assert runtime_warnings == []
    assert result.parameters["source"] == "raw"
    assert science.np.isfinite(result.table["logFC"].to_numpy(dtype=float)).all()


@pytest.mark.parametrize(
    "method,expected_random_state",
    [("logreg", 17), ("wilcoxon", None)],
)
def test_marker_random_seed_is_only_forwarded_to_logreg(adata, science, monkeypatch, method, expected_random_state):
    received = {}

    def rank_genes_groups(*args, **kwargs):
        received.update(kwargs)

    def rank_genes_groups_df(*args, **kwargs):
        return science.pd.DataFrame(
            {
                "group": ["A", "B"],
                "names": ["G1", "G2"],
                "scores": [1.0, 0.5],
                "logfoldchanges": [1.0, 0.5],
                "pvals": [0.01, 0.02],
                "pvals_adj": [0.02, 0.04],
                "pct_nz_group": [0.8, 0.7],
                "pct_nz_reference": [0.2, 0.3],
            }
        )

    monkeypatch.setattr(science.sc.tl, "rank_genes_groups", rank_genes_groups)
    monkeypatch.setattr(science.sc.get, "rank_genes_groups_df", rank_genes_groups_df)
    OpenBioSingleCellMarkerGenes.execute(adata, "group", method, "X", "", 2, True, 17)

    assert received.get("random_state") == expected_random_state
    assert ("random_state" in received) is (method == "logreg")


def test_umap_plot_is_read_only_for_categorical_and_numeric_colors(adata, science):
    adata.obsm["X_umap"] = science.np.arange(adata.n_obs * 2, dtype=float).reshape(adata.n_obs, 2)
    adata.obs["score"] = science.np.linspace(0.0, 1.0, adata.n_obs)
    snapshot_x = dense(adata.X, science).copy()
    snapshot_obs = adata.obs.copy(deep=True)
    snapshot_uns = copy.deepcopy(adata.uns)

    categorical = output_value(OpenBioSingleCellUMAPPlot.execute(adata, "group", 8.0, "viridis"))
    numeric = output_value(OpenBioSingleCellUMAPPlot.execute(adata, "score", 8.0, "viridis"))

    assert categorical.kind == "plot" and categorical.png
    assert numeric.kind == "plot" and numeric.png
    science.np.testing.assert_array_equal(dense(adata.X, science), snapshot_x)
    science.pd.testing.assert_frame_equal(adata.obs, snapshot_obs)
    assert adata.uns == snapshot_uns


def test_umap_plot_preconditions_are_clear(adata, science):
    with pytest.raises(ValueError, match="run UMAP first"):
        OpenBioSingleCellUMAPPlot.execute(adata, "group", 8.0, "viridis")

    adata.obsm["X_umap"] = science.np.ones((adata.n_obs, 1))
    with pytest.raises(ValueError, match="at least two columns"):
        OpenBioSingleCellUMAPPlot.execute(adata, "group", 8.0, "viridis")

    adata.obsm["X_umap"] = science.np.ones((adata.n_obs, 2))
    with pytest.raises(ValueError, match="color column not found"):
        OpenBioSingleCellUMAPPlot.execute(adata, "missing", 8.0, "viridis")
