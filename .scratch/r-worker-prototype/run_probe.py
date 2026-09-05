"""THROWAWAY: exercise a real R-backed node through the existing artifact service.

Run with D:/learn/ComfyUI/.venv/Scripts/python.exe; no production module is modified.
"""
import argparse
import asyncio
import gc
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT), str(ROOT.parent / "ComfyUI")]

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from comfy.cli_args import args as comfy_args
comfy_args.cpu = True
from comfy_api.latest import io

from openbio_singlecell import worker_client
from openbio_singlecell.artifact_codecs import read_anndata, read_table, write_anndata
from openbio_singlecell.artifact_service import (
    current_artifact_runtime, execute_artifact_node, initialize_artifact_service,
)
from openbio_singlecell.artifact_runtime import ArtifactTicket
from openbio_singlecell.contracts import ensure_metadata
from openbio_singlecell.node_types import AnnDataType, TableResultType, PlotResultType, analysis_outputs
from openbio_singlecell.nodes_data import OpenBioSingleCellSubsetObservations


class OpenBioSingleCellRPrototype(io.ComfyNode):
    @classmethod
    def define_schema(cls):
        return io.Schema(
            node_id="OpenBioSingleCellRPrototype", display_name="R PCA prototype", category="prototype",
            inputs=[AnnDataType.Input("adata"), io.String.Input("rscript"),
                    io.String.Input("r_bin"), io.String.Input("r_home")],
            outputs=analysis_outputs(AnnDataType.Output("adata"), TableResultType.Output("table"),
                                     PlotResultType.Output("plot")),
        )


def fixture(storage):
    identifiers = ["cell A", "cell/B", "细胞-三", "NA", "cell,E", "célula"]
    values = np.arange(24, dtype=np.float64).reshape(6, 4) / 4
    values[0, 0] = -0.999999999
    values[2, 1] = 0
    full = ad.AnnData(
        sparse.csr_matrix(values),
        obs=pd.DataFrame({
            "group": pd.Categorical(["naïve", "memory", None, "memory", "naïve", "memory"],
                                    categories=["memory", "naïve", "unused"], ordered=True),
            "nullable_integer": pd.array([1, None, 3, 4, 5, 6], dtype="Int64"),
            "nullable_boolean": pd.array([True, None, False, True, False, True], dtype="boolean"),
        }, index=identifiers),
        var=pd.DataFrame({"symbol": ["gene A", "基因B", "G3", "RAW_ONLY"]},
                         index=["G1", "G2", "G3", "RAW_ONLY"]),
    )
    full.raw = full
    value = full[:, ["G3", "G1", "G2"]].copy()
    if storage == "dense":
        value.X = value.X.toarray()
    value.layers["selected_snapshot"] = value.X.copy()
    value.obsm["old_embedding"] = np.arange(12, dtype=float).reshape(6, 2)
    value.obsp["old_graph"] = sparse.eye(6, format="csr")
    value.uns["nested"] = {"records": [{"label": "来源"}, None, {"values": [1, 2]}]}
    ensure_metadata(value, source={"kind": "prototype", "dataset": "r-file-bridge"})
    return value


def dense(value):
    return value.toarray() if sparse.issparse(value) else np.asarray(value)


async def main(arguments):
    output_root = Path(arguments.output).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    worker_client.WORKER_ENTRY = HERE / "entry.py"  # one prototype-only registration seam
    observations = []
    with tempfile.TemporaryDirectory(prefix="openbio-r-numeric-") as temporary:
        await initialize_artifact_service(temporary, sys.executable)
        runtime = current_artifact_runtime()
        for storage in ("dense", "csr"):
            original = fixture(storage)
            seed_lease = runtime.begin_run()
            input_root = seed_lease.staging_path / "outputs" / "adata"
            input_root.mkdir(parents=True)
            write_anndata(input_root, original)
            seed = seed_lease.publish({"adata": {"kind": "OPENBIO_ANNDATA", "codec": "anndata-h5ad-v1",
                                                "payload": "outputs/adata"}})
            input_ticket = seed["adata"]
            input_path = runtime.resolve(input_ticket) / "data.h5ad"
            before = hashlib.sha256(input_path.read_bytes()).hexdigest()
            returned = await execute_artifact_node(OpenBioSingleCellRPrototype, None, {
                "adata": input_ticket, "rscript": arguments.rscript,
                "r_bin": arguments.r_bin, "r_home": arguments.r_home,
            })
            adata_ticket, table_ticket, plot_ticket, report, code = returned.result
            assert all(isinstance(ticket, ArtifactTicket) for ticket in (adata_ticket, table_ticket, plot_ticket))
            assert "prcomp" in code and report.summary["key_results"]["code_language"] == "R"
            assert not {"totals", "scores", "obs_names"}.intersection(report.summary["key_results"])
            actual = read_anndata(runtime.resolve(adata_ticket))
            expected = dense(original.X)
            np.testing.assert_allclose(actual.obs["r_total"], expected.sum(axis=1), rtol=1e-14, atol=1e-14)
            centered = expected - expected.mean(axis=0)
            u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
            projected = u[:, :2] * singular_values[:2]
            r_scores = actual.obsm["X_r_pca"]
            np.testing.assert_allclose(r_scores @ r_scores.T, projected @ projected.T, atol=1e-10, rtol=1e-10)
            np.testing.assert_array_equal(dense(actual.X), expected)
            assert actual.X.dtype == original.X.dtype
            assert sparse.issparse(actual.X) == sparse.issparse(original.X)
            pd.testing.assert_frame_equal(actual.obs[original.obs.columns], original.obs)
            pd.testing.assert_frame_equal(actual.var, original.var)
            pd.testing.assert_frame_equal(actual.raw.var, original.raw.var)
            np.testing.assert_array_equal(dense(actual.raw.X), dense(original.raw.X))
            np.testing.assert_array_equal(dense(actual.layers["selected_snapshot"]), dense(original.layers["selected_snapshot"]))
            np.testing.assert_array_equal(actual.obsm["old_embedding"], original.obsm["old_embedding"])
            np.testing.assert_array_equal(dense(actual.obsp["old_graph"]), dense(original.obsp["old_graph"]))
            assert actual.uns["nested"] == original.uns["nested"]
            assert hashlib.sha256(input_path.read_bytes()).hexdigest() == before
            table, _ = read_table(runtime.resolve(table_ticket))
            assert table["obs_id"].tolist() == original.obs_names.tolist()
            downstream = await execute_artifact_node(OpenBioSingleCellSubsetObservations, None, {
                "adata": adata_ticket, "column": "group", "values": "memory",
                "invert": False, "missing_policy": "exclude",
            })
            selected = read_anndata(runtime.resolve(downstream.result[0]))
            assert selected.obs_names.tolist() == original.obs_names[[1, 3, 5]].tolist()
            np.testing.assert_array_equal(selected.obsm["X_r_pca"], r_scores[[1, 3, 5]])
            shutil.copyfile(runtime.resolve(adata_ticket) / "data.h5ad", output_root / f"{storage}-output.h5ad")
            shutil.copyfile(runtime.resolve(plot_ticket) / "plot.png", output_root / f"{storage}-r-pca.png")
            observations.append({
                "storage": storage, "shape": list(original.shape), "raw_shape": list(original.raw.shape),
                "r_versions": report.summary["key_results"]["versions"],
                "r_pid": report.summary["key_results"]["r_pid"],
                "wire_values_are_tickets": True, "input_file_unchanged": True,
                "numeric_matches_numpy": True, "axes_and_metadata_preserved": True,
                "raw_layers_graphs_nested_uns_preserved": True,
                "native_rds_roundtrip": report.summary["key_results"]["r_model_roundtrip"],
                "standard_python_downstream_cells": selected.n_obs,
                "r_warnings": report.summary["warnings"],
            })
            del downstream, returned, adata_ticket, table_ticket, plot_ticket, input_ticket, seed, seed_lease
            gc.collect()
        assert runtime.close()
    result = {"prototype": "Python Artifact node -> real Rscript -> Artifact node", "cases": observations}
    (output_root / "numeric-report.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--r-bin", required=True)
    parser.add_argument("--r-home", required=True)
    parser.add_argument("--output", default=str(HERE / "results"))
    asyncio.run(main(parser.parse_args()))
