from __future__ import annotations

import copy

import pytest

from openbio_singlecell.contracts import ensure_metadata
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
    filtered_cells = output_value(OpenBioSingleCellFilterCells.execute(qc, 1, 0, 1, 0))
    filtered = output_value(OpenBioSingleCellFilterGenes.execute(filtered_cells, 1, 0, 1, 0))
    normalized = output_value(OpenBioSingleCellNormalizeTotal.execute(filtered, 10_000.0))
    before_log = dense(normalized.layers["counts"], science).copy()
    logged = output_value(OpenBioSingleCellLog1p.execute(normalized, True))
    variable = output_value(OpenBioSingleCellHighlyVariableGenes.execute(logged, 8, "seurat", False))
    with pytest.warns(UserWarning, match="densifies"):
        scaled = output_value(OpenBioSingleCellScale.execute(variable, 10.0))

    assert science.np.array_equal(dense(logged.layers["counts"], science), before_log)
    assert logged.raw is not None
    assert int(variable.var["highly_variable"].sum()) == 8
    assert not science.sparse.issparse(scaled.X)
    assert len(scaled.uns["openbio_singlecell"]["analysis_history"]) == 7


def test_modifying_nodes_do_not_mutate_upstream_anndata(adata, science):
    snapshot = adata.copy()
    snapshot_uns = copy.deepcopy(adata.uns)

    OpenBioSingleCellCalculateQC.execute(adata, mitochondrial_prefix="MT-")
    OpenBioSingleCellFilterCells.execute(adata, 1, 0, 1, 0)
    OpenBioSingleCellFilterGenes.execute(adata, 1, 0, 1, 0)
    OpenBioSingleCellNormalizeTotal.execute(adata, 10_000.0)
    OpenBioSingleCellLog1p.execute(adata, True)
    OpenBioSingleCellHighlyVariableGenes.execute(adata, 8, "seurat", False)
    with pytest.warns(UserWarning, match="densifies"):
        OpenBioSingleCellScale.execute(adata, 10.0)

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
    plot = output_value(OpenBioSingleCellQCPlots.execute(adata))

    assert qc.uns["openbio_singlecell"]["warnings"] == ["No genes matched mitochondrial prefix 'MT-'."]
    assert plot.kind == "plot"
    assert plot.png
    assert any("shown as zero" in warning for warning in plot.warnings)
    assert_adata_equal(adata, snapshot, science)


def test_qc_plot_rejects_empty_anndata(science):
    empty = science.ad.AnnData(science.np.empty((0, 2)))

    with pytest.raises(ValueError, match="at least one cell and one gene"):
        OpenBioSingleCellQCPlots.execute(empty)
