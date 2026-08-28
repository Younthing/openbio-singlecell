from __future__ import annotations

import json

import pytest

from openbio_singlecell import PLUGIN_VERSION, SCHEMA_VERSION
from openbio_singlecell.contracts import ensure_metadata, record_history
from openbio_singlecell.nodes_input import (
    OpenBioSingleCellAnnDataSummary,
    OpenBioSingleCellLoad10xH5,
    OpenBioSingleCellLoad10xMTX,
    OpenBioSingleCellLoad10xStudy,
    OpenBioSingleCellLoadH5AD,
)
from openbio_singlecell.payload import result_to_payload


def output_value(node_output):
    return node_output.result[0]


def write_10x_mtx(directory, science, matrix, barcodes, features):
    directory.mkdir(parents=True)
    scipy_io = pytest.importorskip("scipy.io")
    scipy_io.mmwrite(directory / "matrix.mtx", science.sparse.csr_matrix(matrix))
    (directory / "barcodes.tsv").write_text("\n".join(barcodes) + "\n", encoding="utf-8")
    (directory / "features.tsv").write_text(
        "".join("\t".join(feature) + "\n" for feature in features),
        encoding="utf-8",
    )


def write_10x_h5(path, science, matrix, barcodes, gene_ids, gene_names, feature_types=None):
    h5py = pytest.importorskip("h5py")
    feature_types = feature_types or ["Gene Expression"] * len(gene_ids)
    stored = science.sparse.csc_matrix(matrix)
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=stored.data)
        group.create_dataset("indices", data=stored.indices)
        group.create_dataset("indptr", data=stored.indptr)
        group.create_dataset("shape", data=stored.shape)
        group.create_dataset("barcodes", data=science.np.asarray(barcodes, dtype="S"))
        features = group.create_group("features")
        features.create_dataset("id", data=science.np.asarray(gene_ids, dtype="S"))
        features.create_dataset("name", data=science.np.asarray(gene_names, dtype="S"))
        features.create_dataset("feature_type", data=science.np.asarray(feature_types, dtype="S"))
        features.create_dataset("genome", data=science.np.asarray(["test"] * len(gene_ids), dtype="S"))
        features.create_dataset("_all_tag_keys", data=science.np.asarray(["genome"], dtype="S"))


def test_binary_file_loaders_use_the_web_upload_widget():
    expected = {
        OpenBioSingleCellLoadH5AD: [".h5ad"],
        OpenBioSingleCellLoad10xH5: [".h5", ".hdf5"],
    }

    for node_class, extensions in expected.items():
        path_input = next(input_ for input_ in node_class.GET_SCHEMA().inputs if input_.id == "path")
        assert path_input.extra_dict["widgetType"] == "OPENBIO_INPUT_FILE_UPLOAD_WIDGET"
        assert path_input.extra_dict["allowed_extensions"] == extensions
        assert path_input.extra_dict["upload_subfolder"] == "openbio-singlecell"


@pytest.fixture
def adata(science):
    counts = science.sparse.csr_matrix(science.np.arange(24).reshape(6, 4))
    obs = science.pd.DataFrame(index=[f"cell_{index}" for index in range(6)])
    var = science.pd.DataFrame(index=["MT-A", "B", "C", "D"])
    value = science.ad.AnnData(counts, obs=obs, var=var)
    ensure_metadata(value, display_name="synthetic", source={"kind": "test"})
    return value


def test_h5ad_and_10x_mtx_loaders(comfy_directories, adata, science):
    input_dir, _, _ = comfy_directories
    h5ad_path = input_dir / "sample.h5ad"
    adata.write_h5ad(h5ad_path)

    loaded = output_value(OpenBioSingleCellLoadH5AD.execute("sample.h5ad", True))

    assert loaded.shape == adata.shape
    assert loaded.uns["openbio_singlecell"]["display_name"] == "sample"
    h5ad_source = loaded.uns["openbio_singlecell"]["source"]
    assert h5ad_source["kind"] == "h5ad"
    assert h5ad_source["path"] == "sample.h5ad"
    assert "fingerprint" not in h5ad_source
    assert str(input_dir.resolve()).lower() not in json.dumps(h5ad_source).lower()

    tenx = input_dir / "tenx"
    write_10x_mtx(
        tenx,
        science,
        adata.X.transpose(),
        list(adata.obs_names),
        [(f"id_{name}", name, "Gene Expression") for name in adata.var_names],
    )

    loaded_10x = output_value(OpenBioSingleCellLoad10xMTX.execute("tenx", "gene_symbols", True, True))

    assert loaded_10x.shape == adata.shape
    assert science.sparse.issparse(loaded_10x.X)
    tenx_source = loaded_10x.uns["openbio_singlecell"]["source"]
    assert loaded_10x.uns["openbio_singlecell"]["display_name"] == "tenx"
    assert tenx_source["kind"] == "10x_mtx"
    assert tenx_source["path"] == "tenx"
    assert tenx_source["files"]["matrix"]["path"] == "tenx/matrix.mtx"
    assert "fingerprint" not in tenx_source
    assert str(input_dir.resolve()).lower() not in json.dumps(tenx_source).lower()


def test_h5ad_loader_boundaries(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    assert "not found" in OpenBioSingleCellLoadH5AD.validate_inputs("missing.h5ad", True).lower()

    (input_dir / "corrupt.h5ad").write_bytes(b"not an hdf5 file")
    with pytest.raises(OSError):
        OpenBioSingleCellLoadH5AD.execute("corrupt.h5ad", True)

    empty = science.ad.AnnData(science.np.empty((0, 2)))
    empty.write_h5ad(input_dir / "empty.h5ad")
    empty_loaded = output_value(OpenBioSingleCellLoadH5AD.execute("empty.h5ad", True))
    assert empty_loaded.shape == (0, 2)
    assert any(
        "empty observation axis" in warning
        for warning in empty_loaded.uns["openbio_singlecell"]["source"]["warnings"]
    )

    duplicate_names = science.ad.AnnData(science.np.eye(2))
    duplicate_names.var_names = ["gene", "gene"]
    duplicate_names.write_h5ad(input_dir / "duplicate.h5ad")
    loaded = output_value(OpenBioSingleCellLoadH5AD.execute("duplicate.h5ad", True))
    assert loaded.var_names.is_unique
    assert loaded.uns["openbio_singlecell"]["source"]["axis_names"]["var_names_repaired"] == 1

    preserved = output_value(OpenBioSingleCellLoadH5AD.execute("duplicate.h5ad", False))
    assert not preserved.var_names.is_unique
    assert preserved.uns["openbio_singlecell"]["source"]["warnings"]

    duplicate_cells = science.ad.AnnData(science.np.eye(2))
    duplicate_cells.obs_names = ["cell", "cell"]
    duplicate_cells.write_h5ad(input_dir / "duplicate_cells.h5ad")
    duplicate_cells_loaded = output_value(
        OpenBioSingleCellLoadH5AD.execute("duplicate_cells.h5ad", True)
    )
    duplicate_source = duplicate_cells_loaded.uns["openbio_singlecell"]["source"]
    assert not duplicate_cells_loaded.obs_names.is_unique
    assert any("duplicate observation-name" in warning for warning in duplicate_source["warnings"])


def test_summary_normalizes_metadata(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    value = science.ad.AnnData(science.np.eye(2))
    value.uns["openbio_singlecell"] = {
        "schema_version": SCHEMA_VERSION,
        "version": PLUGIN_VERSION,
        "display_name": "sample",
        "warnings": 1,
        "analysis_history": {},
        "random_seed": "bad",
    }
    value.write_h5ad(input_dir / "metadata.h5ad")

    loaded = output_value(OpenBioSingleCellLoadH5AD.execute("metadata.h5ad", True))
    summary = output_value(OpenBioSingleCellAnnDataSummary.execute(loaded))
    record_history(loaded, "test", {}, loaded.n_obs, loaded.n_vars)
    metadata = loaded.uns["openbio_singlecell"]

    assert summary.title == "metadata summary"
    assert summary.summary["key_results"]["shape"] == [2, 2]
    assert summary.summary["node_id"] == "OpenBioSingleCellAnnDataSummary"
    assert summary.warnings == []
    assert summary.random_seed == 0
    assert metadata["warnings"] == []
    assert metadata["analysis_history"]["000000"]["operation"] == "test"
    assert metadata["source"]["kind"] == "h5ad"
    assert metadata["source"]["path"] == "metadata.h5ad"


def test_h5ad_loader_redacts_absolute_paths_in_embedded_provenance(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    value = science.ad.AnnData(science.np.eye(2))
    ensure_metadata(
        value,
        source={
            "kind": "external",
            "windows_path": r"C:\private\donor.h5ad",
            "posix_path": "/private/donor.h5ad",
            "relative_path": "portable/donor.h5ad",
            "nested": {"file": "file:///private/donor.h5ad"},
        },
    )
    value.write_h5ad(input_dir / "embedded.h5ad")

    loaded = output_value(OpenBioSingleCellLoadH5AD.execute("embedded.h5ad"))
    embedded = loaded.uns["openbio_singlecell"]["source"]["embedded_source"]

    assert embedded["windows_path"] == "<redacted-absolute-path>"
    assert embedded["posix_path"] == "<redacted-absolute-path>"
    assert embedded["nested"]["file"] == "<redacted-absolute-path>"
    assert embedded["relative_path"] == "portable/donor.h5ad"

    array_value = science.ad.AnnData(science.np.eye(1))
    array_value.uns["openbio_singlecell"] = {
        "schema_version": SCHEMA_VERSION,
        "version": PLUGIN_VERSION,
        "display_name": "array source",
        "warnings": [],
        "analysis_history": {},
        "random_seed": 0,
        "source": {
            "paths": science.np.asarray(["/private/a.h5ad", r"C:\private\b.h5ad"]),
            "byte_paths": science.np.asarray([b"/private/c.h5ad", b"C:\\private\\d.h5ad"]),
        },
    }
    array_value.write_h5ad(input_dir / "embedded_array.h5ad")
    array_loaded = output_value(OpenBioSingleCellLoadH5AD.execute("embedded_array.h5ad"))
    assert array_loaded.uns["openbio_singlecell"]["source"]["embedded_source"]["paths"] == [
        "<redacted-absolute-path>",
        "<redacted-absolute-path>",
    ]
    assert array_loaded.uns["openbio_singlecell"]["source"]["embedded_source"]["byte_paths"] == [
        "<redacted-absolute-path>",
        "<redacted-absolute-path>",
    ]


def test_10x_mtx_loader_boundaries(comfy_directories, science):
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
    scipy_io.mmwrite(mismatched / "matrix.mtx", science.sparse.csr_matrix([[1], [2]]))
    (mismatched / "barcodes.tsv").write_text("cell\n", encoding="utf-8")
    (mismatched / "features.tsv").write_text("id\tgene\tGene Expression\n", encoding="utf-8")
    with pytest.raises(ValueError, match="dimensions do not match"):
        OpenBioSingleCellLoad10xMTX.execute("mismatched", "gene_symbols", True, True)


@pytest.mark.parametrize(
    ("name", "matrix", "barcodes", "features", "message"),
    [
        (
            "duplicate_barcodes",
            [[1, 2]],
            ["cell", "cell"],
            [("id", "gene", "Gene Expression")],
            "duplicate barcode",
        ),
        (
            "blank_barcode",
            [[1, 2]],
            ["cell", ""],
            [("id", "gene", "Gene Expression")],
            "blank barcode",
        ),
        (
            "blank_gene_id",
            [[1]],
            ["cell"],
            [("", "gene", "Gene Expression")],
            "blank stable feature ID",
        ),
        (
            "duplicate_gene_id",
            [[1], [2]],
            ["cell"],
            [("id", "A", "Gene Expression"), ("id", "B", "Gene Expression")],
            "duplicate stable feature ID",
        ),
        (
            "all_zero",
            [[0]],
            ["cell"],
            [("id", "gene", "Gene Expression")],
            "no positive expression",
        ),
    ],
)
def test_10x_mtx_discloses_identity_and_empty_count_advisories(
    comfy_directories, science, name, matrix, barcodes, features, message
):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(input_dir / name, science, matrix, barcodes, features)
    loaded = output_value(OpenBioSingleCellLoad10xMTX.execute(name, "gene_symbols", True, True))
    source = loaded.uns["openbio_singlecell"]["source"]
    assert any(message in warning for warning in source["warnings"])


@pytest.mark.parametrize(
    ("name", "matrix", "message"),
    [
        ("negative", [[-1]], "negative expression"),
        ("fractional", [[1.5]], "non-integer"),
    ],
)
def test_10x_mtx_rejects_values_incompatible_with_count_format(
    comfy_directories, science, name, matrix, message
):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(
        input_dir / name,
        science,
        matrix,
        ["cell"],
        [("id", "gene", "Gene Expression")],
    )
    with pytest.raises(ValueError, match=message):
        OpenBioSingleCellLoad10xMTX.execute(name, "gene_symbols", True, True)


def test_10x_mtx_repairs_only_duplicate_feature_names(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(
        input_dir / "duplicate_symbols",
        science,
        [[1], [2]],
        ["cell"],
        [("id_1", "gene", "Gene Expression"), ("id_2", "gene", "Gene Expression")],
    )

    repaired = output_value(
        OpenBioSingleCellLoad10xMTX.execute("duplicate_symbols", "gene_symbols", True, True)
    )
    preserved = output_value(
        OpenBioSingleCellLoad10xMTX.execute("duplicate_symbols", "gene_symbols", False, True)
    )

    assert repaired.obs_names.tolist() == ["cell"]
    assert repaired.var_names.tolist() == ["gene", "gene-1"]
    assert not preserved.var_names.is_unique
    assert repaired.uns["openbio_singlecell"]["source"]["axis_names"]["var_names_repaired"] == 1


def test_10x_mtx_rejects_ambiguous_file_roles(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_mtx(
        input_dir / "ambiguous",
        science,
        [[1]],
        ["cell"],
        [("id", "gene", "Gene Expression")],
    )
    (input_dir / "ambiguous" / "barcodes.tsv.gz").write_bytes(b"placeholder")
    message = OpenBioSingleCellLoad10xMTX.validate_inputs("ambiguous")
    assert "ambiguous barcodes" in message


def test_10x_study_never_silently_omits_invalid_samples(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[1]],
        ["cell_a"],
        [("id", "gene", "Gene Expression")],
    )
    broken = root / "Sample_B"
    broken.mkdir()
    (broken / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n0 0 0\n")

    validation = OpenBioSingleCellLoad10xStudy.validate_inputs("study")
    assert "Sample_B" in validation
    assert "missing barcodes" in validation
    with pytest.raises(ValueError, match="Sample_B"):
        OpenBioSingleCellLoad10xStudy.fingerprint_inputs("study")
    with pytest.raises(ValueError, match="Sample_B"):
        OpenBioSingleCellLoad10xStudy.execute("study")


def test_10x_study_reports_inner_and_outer_alignment(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[1], [2]],
        ["barcode"],
        [("id_common", "common", "Gene Expression"), ("id_a", "only_a", "Gene Expression")],
    )
    write_10x_mtx(
        root / "Sample_B",
        science,
        [[3], [4]],
        ["barcode"],
        [("id_common", "common", "Gene Expression"), ("id_b", "only_b", "Gene Expression")],
    )

    inner = output_value(OpenBioSingleCellLoad10xStudy.execute("study", join="inner"))
    outer = output_value(OpenBioSingleCellLoad10xStudy.execute("study", join="outer"))

    assert inner.shape == (2, 1)
    assert outer.shape == (2, 3)
    assert science.np.asarray(outer.X.toarray()).tolist() == [[1, 2, 0], [3, 0, 4]]
    assert inner.obs["sample"].astype(str).tolist() == ["Sample_A", "Sample_B"]
    assert inner.obs_names.tolist() == ["barcode-Sample_A", "barcode-Sample_B"]
    inner_source = inner.uns["openbio_singlecell"]["source"]
    outer_source = outer.uns["openbio_singlecell"]["source"]
    assert inner_source["sample_count"] == 2
    assert set(inner_source["samples"]) == {"Sample_A", "Sample_B"}
    assert inner_source["feature_alignment"]["Sample_A"]["dropped_features"] == 1
    assert outer_source["feature_alignment"]["Sample_A"]["zero_filled_features"] == 1
    assert outer_source["warnings"]
    assert outer.var["gene_ids"].tolist() == ["id_common", "id_a", "id_b"]
    assert outer.var["gene_symbols"].tolist() == ["common", "only_a", "only_b"]
    assert outer.var["feature_types"].tolist() == ["Gene Expression"] * 3
    assert str(input_dir.resolve()).lower() not in json.dumps(outer_source).lower()


def test_10x_study_returns_empty_inner_feature_intersection_with_advisory(
    comfy_directories,
    science,
):
    input_dir, _, _ = comfy_directories
    root = input_dir / "disjoint"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[1]],
        ["a"],
        [("id_a", "gene_a", "Gene Expression")],
    )
    write_10x_mtx(
        root / "Sample_B",
        science,
        [[2]],
        ["b"],
        [("id_b", "gene_b", "Gene Expression")],
    )

    output = output_value(OpenBioSingleCellLoad10xStudy.execute("disjoint", join="inner"))
    source = output.uns["openbio_singlecell"]["source"]
    assert output.shape == (2, 0)
    assert source["retained_features"] == 0
    assert any("empty feature intersection" in warning for warning in source["warnings"])


def test_10x_study_rejects_conflicting_symbol_identity(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "A",
        science,
        [[1]],
        ["a"],
        [("id_1", "shared_name", "Gene Expression")],
    )
    write_10x_mtx(
        root / "B",
        science,
        [[1]],
        ["b"],
        [("id_2", "shared_name", "Gene Expression")],
    )
    with pytest.raises(ValueError, match="conflicting identity metadata"):
        OpenBioSingleCellLoad10xStudy.execute("study")


def test_10x_study_aggregates_sample_payload_errors(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    root = input_dir / "study"
    write_10x_mtx(
        root / "Sample_A",
        science,
        [[-1]],
        ["a"],
        [("id_a", "gene_a", "Gene Expression")],
    )
    write_10x_mtx(
        root / "Sample_B",
        science,
        [[1, 2]],
        ["b", "b"],
        [("id_b", "gene_b", "Gene Expression")],
    )

    with pytest.raises(ValueError, match="Invalid 10x Study Sample payloads") as error:
        OpenBioSingleCellLoad10xStudy.execute("study")

    message = str(error.value)
    assert "Sample_A" in message and "negative expression" in message
    assert "Sample_B" not in message


def test_10x_h5_loader(comfy_directories, adata, science):
    input_dir, _, _ = comfy_directories
    path = input_dir / "matrix.h5"
    write_10x_h5(
        path,
        science,
        adata.X.transpose(),
        list(adata.obs_names),
        [f"id_{value}" for value in adata.var_names],
        list(adata.var_names),
    )

    loaded = output_value(OpenBioSingleCellLoad10xH5.execute("matrix.h5", "", True, True))

    assert loaded.shape == adata.shape
    assert loaded.uns["openbio_singlecell"]["display_name"] == "matrix"
    source = loaded.uns["openbio_singlecell"]["source"]
    assert source["kind"] == "10x_h5"
    assert source["path"] == "matrix.h5"
    assert "fingerprint" not in source
    assert str(input_dir.resolve()).lower() not in json.dumps(source).lower()


def test_10x_h5_rejects_duplicate_barcodes_and_invalid_counts(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_h5(
        input_dir / "duplicate.h5",
        science,
        [[1, 2]],
        ["cell", "cell"],
        ["id"],
        ["gene"],
    )
    duplicate = output_value(OpenBioSingleCellLoad10xH5.execute("duplicate.h5"))
    assert not duplicate.obs_names.is_unique
    assert any(
        "duplicate observation-name" in warning
        for warning in duplicate.uns["openbio_singlecell"]["source"]["warnings"]
    )

    write_10x_h5(
        input_dir / "fractional.h5",
        science,
        [[1.5]],
        ["cell"],
        ["id"],
        ["gene"],
    )
    with pytest.raises(ValueError, match="non-integer"):
        OpenBioSingleCellLoad10xH5.execute("fractional.h5")


def test_10x_h5_discloses_feature_filtering(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_10x_h5(
        input_dir / "multimodal.h5",
        science,
        [[1], [2]],
        ["cell"],
        ["gene_id", "antibody_id"],
        ["gene", "antibody"],
        ["Gene Expression", "Antibody Capture"],
    )

    loaded = output_value(OpenBioSingleCellLoad10xH5.execute("multimodal.h5", gex_only=True))
    source = loaded.uns["openbio_singlecell"]["source"]

    assert loaded.var_names.tolist() == ["gene"]
    assert source["input_features"] == 2
    assert source["retained_features"] == 1
    assert source["filtered_features"] == 1
    assert source["feature_types_before_filter"] == {
        "Gene Expression": 1,
        "Antibody Capture": 1,
    }


def test_anndata_summary_is_strict_bounded_read_only_and_code_equivalent(science):
    value = science.ad.AnnData(
        X=science.sparse.csr_matrix(science.np.eye(3)),
        obs=science.pd.DataFrame(index=["first", "second", "third"]),
        var=science.pd.DataFrame(index=["a", "b", "c"]),
    )
    value.layers["counts"] = value.X.copy()
    value.obsm["X_test"] = science.np.ones((3, 2))
    value.varm["loadings"] = science.np.ones((3, 2))
    value.obsp["connectivities"] = science.sparse.eye(3, format="csr")
    value.varp["similarity"] = science.sparse.eye(3, format="csr")
    value.uns["method"] = {"value": 1}
    value.raw = value.copy()
    value.obs_names = ["cell", "cell", "third"]
    before_x = value.X.copy()
    before_obs = value.obs.copy(deep=True)

    report, code = OpenBioSingleCellAnnDataSummary.execute(value).result
    key_results = report.summary["key_results"]

    json.dumps(report.summary, allow_nan=False)
    payload = result_to_payload(report)
    assert key_results["shape"] == [3, 3]
    assert key_results["X"]["storage"] == "sparse_csr"
    assert key_results["obs_index"]["is_unique"] is False
    assert key_results["raw"]["present"] is True
    assert key_results["layers"]["names"] == ["counts"]
    assert payload["summary"]["key_results"]["layers"]["names"] == ["counts"]
    assert key_results["obsp"]["names"] == ["connectivities"]
    assert key_results["varp"]["names"] == ["similarity"]
    assert any("not unique" in warning for warning in report.summary["warnings"])
    namespace = {}
    exec(compile(code, "<anndata-summary>", "exec"), namespace)
    assert namespace["summarize_anndata"](value) == key_results
    assert (value.X != before_x).nnz == 0
    assert value.obs.equals(before_obs)


def test_anndata_summary_bounds_names_without_mutating_invalid_metadata(science):
    value = science.ad.AnnData(science.np.eye(2))
    for index in range(70):
        value.uns[f"key_{index:02d}_" + "x" * 300] = index
    value.uns["openbio_singlecell"] = ["invalid"]
    before = list(value.uns["openbio_singlecell"])

    report, _ = OpenBioSingleCellAnnDataSummary.execute(value).result
    inventory = report.summary["key_results"]["uns"]

    assert inventory["total"] == 71
    assert len(inventory["names"]) == 64
    assert inventory["truncated"] is True
    assert max(map(len, inventory["names"])) <= 256
    assert report.summary["key_results"]["openbio_metadata"]["valid_mapping"] is False
    assert value.uns["openbio_singlecell"] == before


def test_anndata_summary_bounds_embedded_display_name_and_warnings(science):
    value = science.ad.AnnData(science.np.eye(2))
    value.uns["openbio_singlecell"] = {
        "display_name": "d" * 500,
        "warnings": [f"warning-{index}-" + "x" * 500 for index in range(40)],
        "source": {"kind": "source-" + "s" * 500},
        "analysis_history": {},
    }

    report, _ = OpenBioSingleCellAnnDataSummary.execute(value).result
    metadata = report.summary["key_results"]["openbio_metadata"]

    assert len(metadata["display_name"]) == 256
    assert len(metadata["source_kind"]) == 256
    assert len(report.title) == 264
    assert len(report.summary["warnings"]) == 33
    assert all(len(warning) <= 256 for warning in report.summary["warnings"][:-1])
    assert "limited to the first 32" in report.summary["warnings"][-1]
