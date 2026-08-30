from __future__ import annotations

import io

import h5py
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from PIL import Image
from scipy import sparse

from openbio_singlecell.artifact_codecs import (
    read_anndata,
    read_plot,
    read_table,
    write_anndata,
    write_plot,
    write_table,
)


def _assert_matrix_equal(actual, expected) -> None:
    actual_array = actual.toarray() if sparse.issparse(actual) else np.asarray(actual)
    expected_array = expected.toarray() if sparse.issparse(expected) else np.asarray(expected)
    np.testing.assert_array_equal(actual_array, expected_array)


@pytest.mark.parametrize("sparse_x", [False, True])
def test_anndata_round_trip_preserves_complete_state_without_categorizing_strings(tmp_path, sparse_x):
    counts = np.asarray([[1, 0, 3], [0, 2, 4]], dtype=np.float32)
    matrix = sparse.csr_matrix(counts) if sparse_x else counts.copy()
    obs = pd.DataFrame(
        {
            "batch": ["first", "second"],
            "quality": pd.array([1, None], dtype="Int64"),
        },
        index=pd.Index(["cell-a", "cell-b"], name="cell"),
    )
    var = pd.DataFrame(
        {"feature_type": ["gene", "gene", "control"]},
        index=pd.Index(["g1", "g2", "g3"], name="gene"),
    )
    adata = AnnData(matrix, obs=obs, var=var)
    adata.layers["counts"] = matrix.copy()
    adata.obsm["X_pca"] = np.asarray([[1.0, 2.0], [3.0, 4.0]])
    adata.obsp["connectivities"] = sparse.csr_matrix([[1.0, 0.25], [0.25, 1.0]])
    adata.uns["analysis"] = {
        "method": "test",
        "parameters": {"neighbors": 7},
        "analysis_history": {
            "snapshot": {"operation": "snapshot_expression"},
            "normalized": {"operation": "scale_to_layer"},
        },
        "set_accounting_preview": [
            {"set": "pathway-a", "features": ["g1", "g2"], "retained": 2},
            {"set": "pathway-b", "features": ["g3"], "retained": 1},
        ],
    }
    adata.raw = adata.copy()
    preview = adata.uns["analysis"]["set_accounting_preview"]

    payloads = write_anndata(tmp_path, adata)
    restored = read_anndata(tmp_path)

    assert adata.uns["analysis"]["set_accounting_preview"] is preview
    assert "__openbio_json_uns_v1__" not in adata.uns
    assert payloads == [{"path": "data.h5ad", "size": (tmp_path / "data.h5ad").stat().st_size}]
    with h5py.File(tmp_path / "data.h5ad", "r") as stored:
        matrix_dataset = stored["X/data"] if sparse_x else stored["X"]
        assert matrix_dataset.compression is None
    assert restored.shape == (2, 3)
    assert restored.obs_names.tolist() == ["cell-a", "cell-b"]
    assert restored.var_names.tolist() == ["g1", "g2", "g3"]
    assert not isinstance(restored.obs["batch"].dtype, pd.CategoricalDtype)
    assert not isinstance(restored.var["feature_type"].dtype, pd.CategoricalDtype)
    assert str(restored.obs["quality"].dtype) == "Int64"
    _assert_matrix_equal(restored.X, matrix)
    _assert_matrix_equal(restored.layers["counts"], matrix)
    _assert_matrix_equal(restored.raw.X, matrix)
    np.testing.assert_array_equal(restored.obsm["X_pca"], adata.obsm["X_pca"])
    _assert_matrix_equal(restored.obsp["connectivities"], adata.obsp["connectivities"])
    assert restored.uns["analysis"]["method"] == "test"
    assert restored.uns["analysis"]["parameters"]["neighbors"] == 7
    assert list(restored.uns["analysis"]["analysis_history"]) == ["snapshot", "normalized"]
    assert restored.uns["analysis"]["set_accounting_preview"] == [
        {"set": "pathway-a", "features": ["g1", "g2"], "retained": 2},
        {"set": "pathway-b", "features": ["g3"], "retained": 1},
    ]


def test_table_round_trip_preserves_order_nullable_categories_and_explicit_json_lists(tmp_path):
    table = pd.DataFrame(
        {
            "count": pd.array([4, None], dtype="Int64"),
            "score": pd.array([1.5, None], dtype="Float64"),
            "accepted": pd.array([True, None], dtype="boolean"),
            "label": pd.array(["first", None], dtype="string"),
            "cluster": pd.Categorical(
                ["beta", "alpha"],
                categories=["alpha", "beta", "unused"],
                ordered=True,
            ),
            "source": pd.Series(["manual", None], dtype=object, index=["row-b", "row-a"]),
            "motif_ids": pd.Series([["M1", "M2"], []], dtype=object, index=["row-b", "row-a"]),
        },
        index=pd.Index(["row-b", "row-a"], name="row_id"),
    )
    metadata = {
        "kind": "table",
        "title": "membership",
        "parameters": {"threshold": 0.25},
        "description": "Complete membership table.",
        "warnings": [],
        "input_cells": 2,
        "input_genes": 3,
        "random_seed": 0,
        "elapsed_seconds": 0.5,
        "source": {"operation": "test"},
    }

    payloads = write_table(tmp_path, table, metadata, json_list_columns={"motif_ids"})
    restored, restored_metadata = read_table(tmp_path)

    assert [payload["path"] for payload in payloads] == ["data.jsonl", "schema.json", "result.json"]
    pd.testing.assert_frame_equal(restored, table)
    assert restored.columns.tolist() == list(table.columns)
    assert restored.index.tolist() == ["row-b", "row-a"]
    assert restored.index.name == "row_id"
    assert restored_metadata == metadata


def test_plot_round_trip_validates_png_and_preserves_result_metadata(tmp_path):
    buffer = io.BytesIO()
    Image.new("RGB", (2, 1), color=(12, 34, 56)).save(buffer, format="PNG")
    png = buffer.getvalue()
    metadata = {"kind": "plot", "title": "embedding", "warnings": []}

    payloads = write_plot(tmp_path, png, metadata)
    restored_png, restored_metadata = read_plot(tmp_path)

    assert [payload["path"] for payload in payloads] == ["plot.png", "result.json"]
    assert restored_png == png
    assert restored_metadata == metadata


@pytest.mark.parametrize(
    ("table", "message"),
    [
        (pd.DataFrame([[1]], columns=[7]), "column names must be strings"),
        (pd.DataFrame([[1, 2]], columns=["value", "value"]), "column names must be unique"),
        (
            pd.DataFrame(
                [[1]],
                index=pd.MultiIndex.from_tuples([("sample", "cell")]),
                columns=["value"],
            ),
            "MultiIndex values",
        ),
        (
            pd.DataFrame([[1]], columns=pd.MultiIndex.from_tuples([("value", "mean")])),
            "MultiIndex columns",
        ),
    ],
)
def test_table_writer_rejects_nonportable_axes(tmp_path, table, message):
    with pytest.raises((TypeError, ValueError), match=message):
        write_table(tmp_path, table, {"kind": "table"})


@pytest.mark.parametrize(
    ("value", "json_list_columns", "message"),
    [
        ({"arbitrary": "object"}, (), "non-scalar"),
        (["not", "declared"], (), "non-scalar"),
        (float("nan"), (), "non-finite"),
        (float("inf"), (), "non-finite"),
        ([1.0, float("inf")], ("value",), "non-finite"),
    ],
)
def test_table_writer_rejects_arbitrary_objects_and_nonfinite_numbers(
    tmp_path,
    value,
    json_list_columns,
    message,
):
    table = pd.DataFrame({"value": pd.Series([value], dtype=object)})

    with pytest.raises((TypeError, ValueError), match=message):
        write_table(
            tmp_path,
            table,
            {"kind": "table"},
            json_list_columns=json_list_columns,
        )


@pytest.mark.parametrize("writer", ["table", "plot"])
def test_result_metadata_rejects_nonfinite_json_numbers(tmp_path, writer):
    metadata = {"kind": writer, "elapsed_seconds": float("nan")}

    with pytest.raises(ValueError, match="non-finite"):
        if writer == "table":
            write_table(tmp_path, pd.DataFrame({"value": [1]}), metadata)
        else:
            buffer = io.BytesIO()
            Image.new("RGB", (1, 1)).save(buffer, format="PNG")
            write_plot(tmp_path, buffer.getvalue(), metadata)


def test_plot_writer_rejects_bytes_with_only_a_png_signature(tmp_path):
    with pytest.raises(ValueError, match="valid PNG"):
        write_plot(tmp_path, b"\x89PNG\r\n\x1a\nnot-an-image", {"kind": "plot"})


def test_table_reader_rejects_numbers_that_overflow_to_infinity(tmp_path):
    write_table(tmp_path, pd.DataFrame({"value": [1.0]}), {"kind": "table"})
    (tmp_path / "result.json").write_text('{"elapsed_seconds":1e999}', encoding="utf-8")

    with pytest.raises(ValueError, match="non-finite"):
        read_table(tmp_path)
