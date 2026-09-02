from __future__ import annotations

import hashlib
import uuid
from pathlib import Path

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, read_table, write_anndata
from openbio_singlecell.contracts import ensure_metadata, record_history
from openbio_singlecell.operations_results import (
    embedding_plot,
    filter_marker_genes,
    marker_expression_plot,
    marker_genes,
    pca_metadata_associations,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse


def _artifact_descriptor(root: Path) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": "OPENBIO_ANNDATA",
        "codec": "anndata-h5ad-v1",
        "path": str(root.resolve()),
    }


def _table_descriptor(root: Path) -> dict[str, str]:
    return {
        "type": "artifact",
        "kind": "OPENBIO_SINGLE_CELL_TABLE",
        "codec": "table-jsonl-v1",
        "path": str(root.resolve()),
    }


def _marker_adata(science):
    values = science.np.log1p(
        science.np.asarray(
            [
                [9.0, 1.0, 2.0],
                [8.0, 2.0, 1.0],
                [7.0, 1.0, 2.0],
                [1.0, 9.0, 2.0],
                [2.0, 8.0, 1.0],
                [1.0, 7.0, 2.0],
            ]
        )
    )
    adata = science.ad.AnnData(
        X=science.sparse.csr_matrix(values),
        obs=science.pd.DataFrame(
            {"cluster": science.pd.Categorical(["A"] * 3 + ["B"] * 3)},
            index=[f"cell-{index}" for index in range(6)],
        ),
        var=science.pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    adata.layers["log1p_norm"] = adata.X.copy()
    ensure_metadata(adata, source={"kind": "test"})
    record_history(
        adata,
        "normalize_to_layer",
        {"source": "X", "output_layer": "log1p_norm", "transform": "log1p"},
        int(adata.n_obs),
        int(adata.n_vars),
    )
    return adata


def test_marker_genes_operation_writes_two_portable_tables_without_rewriting_input(tmp_path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, _marker_adata(science))
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "run.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = marker_genes(
        context,
        {"adata": _artifact_descriptor(input_root)},
        {
            "groupby": "cluster",
            "method": "wilcoxon",
            "source": {"source": "layer", "layer_name": "log1p_norm"},
            "n_genes": 0,
            "tie_correct": True,
            "max_output_rows": 1_000,
            "max_working_memory_gib": 4.0,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["table", "universe", "summary", "code"]
    assert [record["codec"] for record in records[:2]] == ["table-jsonl-v1", "table-jsonl-v1"]
    table, table_metadata = read_table(staging / records[0]["payload"])
    universe, universe_metadata = read_table(staging / records[1]["payload"])
    assert table_metadata["kind"] == universe_metadata["kind"] == "table"
    assert len(table) == 6
    assert universe["gene"].tolist() == ["g1", "g2", "g3"]
    assert hashlib.sha256(input_path.read_bytes()).digest() == before


def test_filter_marker_genes_operation_returns_the_original_universe_ticket_reference(tmp_path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, _marker_adata(science))
    ranking_root = tmp_path / "ranking.partial"
    ranking_root.mkdir()
    ranking_context = OperationContext(ranking_root, str(uuid.uuid4()))
    ranked = marker_genes(
        ranking_context,
        {"adata": _artifact_descriptor(input_root)},
        {
            "groupby": "cluster",
            "method": "wilcoxon",
            "source": {"source": "layer", "layer_name": "log1p_norm"},
            "n_genes": 0,
            "tie_correct": True,
            "max_output_rows": 1_000,
            "max_working_memory_gib": 4.0,
        },
    )
    table_root = ranking_root / ranked[0]["payload"]
    universe_root = ranking_root / ranked[1]["payload"]
    before = {
        path.relative_to(ranking_root): hashlib.sha256(path.read_bytes()).digest()
        for path in ranking_root.rglob("*")
        if path.is_file()
    }
    filter_root = tmp_path / "filter.partial"
    filter_root.mkdir()
    context = OperationContext(filter_root, str(uuid.uuid4()))

    records = filter_marker_genes(
        context,
        {"table": _table_descriptor(table_root), "universe": _table_descriptor(universe_root)},
        {
            "min_log2_fold_change": -100.0,
            "min_fraction_in_group": 0.0,
            "max_fraction_reference": 1.0,
            "max_p_adjusted": 1.0,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert records[0]["codec"] == "table-jsonl-v1"
    assert records[1] == {"type": "input_ref", "name": "universe", "input": "universe"}
    assert before == {
        path.relative_to(ranking_root): hashlib.sha256(path.read_bytes()).digest()
        for path in ranking_root.rglob("*")
        if path.is_file()
    }


def test_embedding_plot_operation_writes_png_without_rewriting_input(tmp_path, science):
    adata = _marker_adata(science)
    adata.obsm["X_umap"] = science.np.arange(12, dtype=float).reshape(6, 2)
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = embedding_plot(
        context,
        {"adata": _artifact_descriptor(input_root)},
        {
            "embedding_key": "X_umap",
            "x_dimension": 1,
            "y_dimension": 2,
            "color": {"color": "obs", "obs_key": "cluster", "color_mode": "auto"},
            "point_size": 8.0,
            "continuous_color_map": "viridis",
            "categorical_palette": "tab20",
            "sort_order": True,
            "missing_color": "lightgray",
            "legend_policy": "automatic",
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert records[0]["codec"] == "plot-png-v1"
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert metadata["kind"] == "plot"
    assert hashlib.sha256(input_path.read_bytes()).digest() == before


def test_marker_expression_operation_writes_png_artifact(tmp_path, science):
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, _marker_adata(science))
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = marker_expression_plot(
        context,
        {"adata": _artifact_descriptor(input_root)},
        {
            "genes": "g1,g2",
            "groupby": "cluster",
            "plot": {
                "plot": "dotplot",
                "standard_scale": "var",
                "expression_cutoff": 0.0,
                "mean_only_expressed": False,
            },
            "source": {"source": "layer", "layer_name": "log1p_norm"},
            "group_order": "observed",
            "random_seed": 0,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert records[0]["codec"] == "plot-png-v1"
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert metadata["parameters"]["genes"] == ["g1", "g2"]
    assert hashlib.sha256(input_path.read_bytes()).digest() == before


def test_pca_metadata_associations_operation_writes_portable_table(tmp_path, science):
    adata = _marker_adata(science)
    adata.obs["sample"] = ["s1", "s1", "s2", "s2", "s3", "s3"]
    adata.obs["age"] = [10.0, 10.0, 20.0, 20.0, 30.0, 30.0]
    adata.obsm["X_pca"] = science.np.asarray(
        [[0.0, 2.0], [0.2, 2.1], [1.0, 1.0], [1.2, 1.1], [2.0, 0.0], [2.2, 0.1]]
    )
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "pca.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = pca_metadata_associations(
        context,
        {"adata": _artifact_descriptor(input_root)},
        {
            "use_rep": "X_pca",
            "sample_key": "sample",
            "categorical_obs_keys": "",
            "continuous_obs_keys": "age",
            "alpha": 0.05,
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert records[0]["codec"] == "table-jsonl-v1"
    table, metadata = read_table(staging / records[0]["payload"])
    assert table["component"].tolist() == ["PC1", "PC2"]
    assert metadata["kind"] == "table"
    assert hashlib.sha256(input_path.read_bytes()).digest() == before


def test_results_nodes_are_schema_only_and_operations_do_not_import_comfy():
    package = Path(__file__).parents[1] / "openbio_singlecell"
    node_source = (package / "nodes_results.py").read_text(encoding="utf-8")
    operation_source = (package / "operations_results.py").read_text(encoding="utf-8")

    assert "def execute" not in node_source
    assert "dependencies" not in node_source
    assert "comfy_api" not in operation_source
    assert "folder_paths" not in operation_source
    assert "nodes_" not in operation_source
