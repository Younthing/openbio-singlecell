from __future__ import annotations

import pytest

from openbio_singlecell import dependencies
from openbio_singlecell.contracts import ensure_metadata, record_history
from openbio_singlecell.nodes_input import (
    OpenBioSingleCellAnnDataSummary,
    OpenBioSingleCellLoad10xH5,
    OpenBioSingleCellLoad10xMTX,
    OpenBioSingleCellLoadH5AD,
)

pytestmark = pytest.mark.skipif(
    not dependencies.AVAILABLE,
    reason="Scientific dependencies are unavailable.",
)


def output_value(node_output):
    return node_output.result[0]


@pytest.fixture
def adata():
    counts = dependencies.sparse.csr_matrix(dependencies.np.arange(24).reshape(6, 4))
    obs = dependencies.pd.DataFrame(index=[f"cell_{index}" for index in range(6)])
    var = dependencies.pd.DataFrame(index=["MT-A", "B", "C", "D"])
    value = dependencies.ad.AnnData(counts, obs=obs, var=var)
    ensure_metadata(value, display_name="synthetic", source={"kind": "test"})
    return value


def test_h5ad_and_10x_mtx_loaders(comfy_directories, adata):
    input_dir, _, _ = comfy_directories
    h5ad_path = input_dir / "sample.h5ad"
    adata.write_h5ad(h5ad_path)

    loaded = output_value(OpenBioSingleCellLoadH5AD.execute("sample.h5ad", True))

    assert loaded.shape == adata.shape
    assert loaded.uns["openbio_singlecell"]["display_name"] == "synthetic"

    tenx = input_dir / "tenx"
    tenx.mkdir()
    scipy_io = pytest.importorskip("scipy.io")
    scipy_io.mmwrite(tenx / "matrix.mtx", adata.X.transpose())
    (tenx / "barcodes.tsv").write_text("\n".join(adata.obs_names) + "\n", encoding="utf-8")
    with open(tenx / "features.tsv", "w", encoding="utf-8") as handle:
        for name in adata.var_names:
            handle.write(f"id_{name}\t{name}\tGene Expression\n")

    loaded_10x = output_value(OpenBioSingleCellLoad10xMTX.execute("tenx", "gene_symbols", True, True))

    assert loaded_10x.shape == adata.shape
    assert dependencies.sparse.issparse(loaded_10x.X)
    assert loaded_10x.uns["openbio_singlecell"]["display_name"] == "tenx"


def test_h5ad_loader_boundaries(comfy_directories):
    input_dir, _, _ = comfy_directories
    assert "not found" in OpenBioSingleCellLoadH5AD.validate_inputs("missing.h5ad", True).lower()

    (input_dir / "corrupt.h5ad").write_bytes(b"not an hdf5 file")
    with pytest.raises(OSError):
        OpenBioSingleCellLoadH5AD.execute("corrupt.h5ad", True)

    empty = dependencies.ad.AnnData(dependencies.np.empty((0, 2)))
    empty.write_h5ad(input_dir / "empty.h5ad")
    with pytest.raises(ValueError, match="at least one cell and one gene"):
        OpenBioSingleCellLoadH5AD.execute("empty.h5ad", True)

    duplicate_names = dependencies.ad.AnnData(dependencies.np.eye(2))
    duplicate_names.var_names = ["gene", "gene"]
    duplicate_names.write_h5ad(input_dir / "duplicate.h5ad")
    with pytest.warns(UserWarning, match="Variable names are not unique"):
        loaded = output_value(OpenBioSingleCellLoadH5AD.execute("duplicate.h5ad", True))
    assert loaded.var_names.is_unique


def test_summary_normalizes_metadata(comfy_directories):
    input_dir, _, _ = comfy_directories
    value = dependencies.ad.AnnData(dependencies.np.eye(2))
    value.uns["openbio_singlecell"] = {
        "display_name": "sample",
        "warnings": 1,
        "analysis_history": 42,
        "random_seed": "bad",
    }
    value.write_h5ad(input_dir / "metadata.h5ad")

    loaded = output_value(OpenBioSingleCellLoadH5AD.execute("metadata.h5ad", True))
    summary = output_value(OpenBioSingleCellAnnDataSummary.execute(loaded))
    record_history(loaded, "test", {}, loaded.n_obs, loaded.n_vars)
    metadata = loaded.uns["openbio_singlecell"]

    assert summary.title == "sample summary"
    assert summary.summary["shape"] == [2, 2]
    assert summary.warnings == []
    assert summary.random_seed == 0
    assert metadata["warnings"] == []
    assert metadata["analysis_history"]["000000"]["operation"] == "test"


def test_10x_mtx_loader_boundaries(comfy_directories):
    input_dir, _, _ = comfy_directories
    missing = input_dir / "missing"
    missing.mkdir()
    (missing / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n0 0 0\n",
        encoding="utf-8",
    )
    validation = OpenBioSingleCellLoad10xMTX.validate_inputs("missing")
    assert "missing barcodes" in validation

    mismatched = input_dir / "mismatched"
    mismatched.mkdir()
    scipy_io = pytest.importorskip("scipy.io")
    scipy_io.mmwrite(mismatched / "matrix.mtx", dependencies.sparse.csr_matrix([[1], [2]]))
    (mismatched / "barcodes.tsv").write_text("cell\n", encoding="utf-8")
    (mismatched / "features.tsv").write_text("id\tgene\tGene Expression\n", encoding="utf-8")
    with pytest.raises(ValueError, match="dimensions do not match"):
        OpenBioSingleCellLoad10xMTX.execute("mismatched", "gene_symbols", True, True)


def test_10x_h5_loader(comfy_directories, adata):
    h5py = pytest.importorskip("h5py")
    input_dir, _, _ = comfy_directories
    path = input_dir / "matrix.h5"
    matrix = dependencies.sparse.csc_matrix(adata.X.transpose())
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=matrix.data)
        group.create_dataset("indices", data=matrix.indices)
        group.create_dataset("indptr", data=matrix.indptr)
        group.create_dataset("shape", data=matrix.shape)
        group.create_dataset("barcodes", data=dependencies.np.asarray(adata.obs_names, dtype="S"))
        features = group.create_group("features")
        features.create_dataset("id", data=dependencies.np.asarray([f"id_{v}" for v in adata.var_names], dtype="S"))
        features.create_dataset("name", data=dependencies.np.asarray(adata.var_names, dtype="S"))
        features.create_dataset(
            "feature_type",
            data=dependencies.np.asarray(["Gene Expression"] * adata.n_vars, dtype="S"),
        )
        features.create_dataset("genome", data=dependencies.np.asarray(["test"] * adata.n_vars, dtype="S"))
        features.create_dataset("_all_tag_keys", data=dependencies.np.asarray(["genome"], dtype="S"))

    loaded = output_value(OpenBioSingleCellLoad10xH5.execute("matrix.h5", "", True, True))

    assert loaded.shape == adata.shape
    assert loaded.uns["openbio_singlecell"]["display_name"] == "matrix"
