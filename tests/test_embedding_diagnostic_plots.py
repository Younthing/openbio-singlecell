from __future__ import annotations

import json

import pytest

from openbio_singlecell.nodes_embedding import (
    OpenBioSingleCellNeighborGraphDiagnosticsPlot,
    OpenBioSingleCellPCALoadingsPlot,
)
from openbio_singlecell.operations_embedding import (
    neighbor_graph_diagnostics_plot_owned,
    pca_loadings_plot_owned,
)

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def test_pca_loadings_plot_schema_and_known_component(science):
    schema = OpenBioSingleCellPCALoadingsPlot.define_schema()
    assert schema.category == "openbio/single-cell/dimension-reduction"
    assert [item.id for item in schema.inputs] == ["adata", "component", "n_genes"]
    assert schema.inputs[2].max == 2**31 - 1
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]

    adata = science.ad.AnnData(
        science.np.zeros((3, 4)),
        var=science.pd.DataFrame(index=["A", "B", "C", "D"]),
    )
    adata.varm["PCs"] = science.np.asarray(
        [[0.8, -0.1], [-0.6, 0.7], [0.2, -0.5], [0.1, 0.3]],
        dtype=float,
    )
    before = adata.copy()

    plotted, report, code = pca_loadings_plot_owned(adata, component=1, n_genes=3)

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["key_results"]["genes"] == ["A", "B", "C"]
    assert report.summary["key_results"]["loadings"] == [0.8, -0.6, 0.2]
    assert report.summary["key_results"]["component"] == 1
    json.dumps(report.summary, allow_nan=False)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_pca_loadings"](adata) == plotted.png
    science.np.testing.assert_array_equal(adata.varm["PCs"], before.varm["PCs"])


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda adata, np: adata.varm.__setitem__("PCs", np.ones((4, 1))), "only 1 components"),
        (lambda adata, np: adata.varm.__setitem__("PCs", np.full((4, 2), np.nan)), "non-finite"),
    ],
)
def test_pca_loadings_plot_rejects_invalid_component_evidence(science, mutation, message):
    adata = science.ad.AnnData(science.np.zeros((2, 4)))
    adata.var_names = ["A", "B", "C", "D"]
    adata.varm["PCs"] = science.np.ones((4, 2))
    mutation(adata, science.np)
    with pytest.raises(ValueError, match=message):
        pca_loadings_plot_owned(adata, component=2, n_genes=2)


def test_pca_loadings_plot_preserves_expert_gene_count_with_generated_parity(science):
    adata = science.ad.AnnData(science.np.zeros((2, 51)))
    adata.varm["PCs"] = science.np.ones((51, 1))

    plotted, report, code = pca_loadings_plot_owned(adata, component=1, n_genes=51)

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["key_results"]["requested_genes"] == 51
    assert report.summary["key_results"]["plotted_genes"] == 51
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_pca_loadings"](adata) == plotted.png


def test_neighbor_graph_diagnostics_plot_schema_and_stored_graph(science):
    schema = OpenBioSingleCellNeighborGraphDiagnosticsPlot.define_schema()
    assert schema.category == "openbio/single-cell/dimension-reduction"
    assert [item.id for item in schema.inputs] == ["adata", "neighbors_key", "max_working_memory_gib"]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]

    adata = science.ad.AnnData(science.np.zeros((4, 0)))
    adata.obs_names = ["c0", "c1", "c2", "c3"]
    connectivities = science.sparse.csr_matrix(
        [[0.0, 1.0, 0.5, 0.0], [1.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 2.0], [0.0, 0.0, 2.0, 0.0]]
    )
    distances = science.sparse.csr_matrix(
        [[0.0, 0.2, 0.8, 0.0], [0.2, 0.0, 0.0, 0.0], [0.8, 0.0, 0.0, 0.4], [0.0, 0.0, 0.4, 0.0]]
    )
    adata.obsp["connectivities"] = connectivities
    adata.obsp["distances"] = distances
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": 2, "metric": "euclidean"},
    }
    before_connectivities = connectivities.copy()

    plotted, report, code = neighbor_graph_diagnostics_plot_owned(adata)

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["connected_components"] == 1
    assert details["component_sizes"] == [4]
    assert details["degree_values"] == [2, 1, 2, 1]
    assert details["positive_distance_values"] == [0.2, 0.8, 0.2, 0.8, 0.4, 0.4]
    json.dumps(report.summary, allow_nan=False)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_neighbor_graph_diagnostics"](adata) == plotted.png
    assert (adata.obsp["connectivities"] != before_connectivities).nnz == 0


def test_neighbor_graph_diagnostics_rejects_missing_distance_bundle(science):
    adata = science.ad.AnnData(science.np.zeros((3, 0)))
    adata.obsp["connectivities"] = science.sparse.csr_matrix(
        [[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
    )
    adata.uns["neighbors"] = {"connectivities_key": "connectivities"}

    with pytest.raises(ValueError, match="distances_key"):
        neighbor_graph_diagnostics_plot_owned(adata)


def test_neighbor_graph_diagnostics_preflights_sparse_working_memory(science):
    adata = science.ad.AnnData(science.np.zeros((3, 0)))
    adata.obsp["connectivities"] = science.sparse.csr_matrix(
        [[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
    )
    adata.obsp["distances"] = science.sparse.csr_matrix(
        [[0.0, 0.2, 0.0], [0.2, 0.0, 0.3], [0.0, 0.3, 0.0]]
    )
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
    }

    with pytest.raises(MemoryError, match="working memory"):
        neighbor_graph_diagnostics_plot_owned(adata, max_working_memory_gib=1e-9)
