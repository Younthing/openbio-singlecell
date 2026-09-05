"""THROWAWAY operation: reuse the Python worker, codecs and publisher unchanged."""
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy import sparse

from openbio_singlecell.analysis_reporting import AnalysisReference, make_analysis_report
from openbio_singlecell.analysis_utils import make_plot_result, make_table_result
from openbio_singlecell.operations_input import (
    analysis_outputs, read_anndata_input, write_anndata_output, write_plot_output, write_table_output,
)
from openbio_singlecell.worker_protocol import register_operation

SCRIPT = Path(__file__).with_name("r_analysis.R")


@register_operation("openbio.node.rprototype")
def r_prototype(context, inputs, parameters):
    started = time.perf_counter()
    adata = read_anndata_input(inputs)
    bridge = context.output_root / "r-bridge"
    bridge.mkdir()
    mmwrite(bridge / "expression.mtx", sparse.coo_matrix(adata.X), symmetry="general", precision=17)
    manifest = {
        "matrix": str(bridge / "expression.mtx"),
        "obs_names": adata.obs_names.tolist(), "var_names": adata.var_names.tolist(),
        "result": str(bridge / "result.json"), "plot": str(bridge / "plot.png"),
        "model": str(bridge / "model.rds"),
    }
    manifest_path = bridge / "request.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False), encoding="utf-8")
    environment = os.environ.copy()
    environment["R_HOME"] = parameters["r_home"]
    environment["PATH"] = os.pathsep.join([
        parameters["r_bin"], str(Path(parameters["rscript"]).parent), environment["PATH"],
    ])
    stderr_path = bridge / "r-stderr.log"
    with stderr_path.open("wb") as stderr:
        process = subprocess.run(
            [parameters["rscript"], "--vanilla", str(SCRIPT), str(manifest_path)],
            env=environment, stderr=stderr,
        )
    with stderr_path.open("rb") as stderr:
        stderr.seek(max(0, stderr_path.stat().st_size - 64 * 1024))
        error_tail = stderr.read().decode("utf-8", errors="replace")
    if process.returncode:
        raise RuntimeError(f"Rscript exited {process.returncode}: {error_tail}")
    result = json.loads((bridge / "result.json").read_text(encoding="utf-8"))
    assert result["obs_names"] == adata.obs_names.tolist()
    adata.obs["r_total"] = np.asarray(result["totals"], dtype=float)
    adata.obsm["X_r_pca"] = np.asarray(result["scores"], dtype=float)
    adata.uns["r_prototype"] = {"versions": result["versions"], "orientation": "cells_by_genes"}
    fields = dict(
        operation="r_prototype", parameters={}, warnings=[], input_cells=adata.n_obs,
        input_genes=adata.n_vars, started_at=started,
    )
    table = make_table_result(
        table=pd.DataFrame({"obs_id": adata.obs_names, "r_total": result["totals"]}),
        title="R cell totals", description="Real R Matrix::rowSums results, aligned by explicit cell IDs.", **fields,
    )
    plot = make_plot_result(
        png=(bridge / "plot.png").read_bytes(), title="R PCA", description="Real R graphics::plot PNG.", **fields,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellRPrototype", title="R file bridge prototype", operation="r_prototype",
        methods="Python decoded the input artifact; R read Matrix Market, ran stats::prcomp and Matrix::rowSums; Python attached only those results.",
        results=f"R analyzed {adata.n_obs} cells and {adata.n_vars} features.",
        key_results={"r_pid": result["r_pid"], "versions": result["versions"],
                     "shape": list(adata.shape), "n_components": 2,
                     "code_language": "R", "r_model_roundtrip": True},
        parameters={"center": True, "scale": False, "n_components": 2},
        references=[AnalysisReference(citation="R stats::prcomp documentation", url="https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html", kind="software")],
        software_packages=["anndata", "numpy", "scipy"], warnings=[error_tail] if error_tail else [],
        limitations=["Throwaway bridge: Matrix Market I/O and dense prcomp suit this small verification dataset only."],
        input_cells=adata.n_obs, input_genes=adata.n_vars, started_at=started,
        code=SCRIPT.read_text(encoding="utf-8"),
    )
    return analysis_outputs(
        report, code, write_anndata_output(context, adata),
        write_table_output(context, table, kind="OPENBIO_SINGLE_CELL_TABLE"),
        write_plot_output(context, plot, kind="OPENBIO_SINGLE_CELL_PLOT"),
    )
