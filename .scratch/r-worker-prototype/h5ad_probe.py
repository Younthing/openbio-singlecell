"""THROWAWAY: measure native R H5AD round-trip fidelity; no production changes.

Run from repository root with the project Python environment. R packages/runtime
are isolated in ~/.cache/openbio-r-prototype. Generated files stay beside this.
"""
from __future__ import annotations

import json
from importlib.metadata import version
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from openbio_singlecell.artifact_codecs import read_anndata, write_anndata

HERE = Path(__file__).resolve().parent
R_PREFIX = Path.home() / ".cache/openbio-r-prototype"


def fixture(object_strings=False):
    index = lambda values, name=None: pd.Index(values, name=name, dtype=object if object_strings else None)
    obs = pd.DataFrame({
        "group": pd.Categorical(["control", None, "treated", "control"],
                                categories=index(["treated", "control", "unused"]), ordered=True),
        "selected": np.array([True, False, True, False], dtype=bool),
        "n_counts": np.array([1, 2, 3, 4], dtype=np.int64),
        "nullable_count": pd.array([1, None, 3, 4], dtype="Int64"),
    }, index=index(["cell-d", "cell-b", "cell-a", "cell-c"], name="cell_id"))
    var = pd.DataFrame({"highly_variable": [True, False, True],
                        "n_cells": np.array([2, 3, 4], dtype=np.int64)},
                       index=index(["gene-z", "gene-a", "gene-m"], name="gene_id"))
    x = np.array([[1, 0, 2], [0, 3, 0], [4, 0, 5], [0, 6, 7]], dtype=np.float32)
    a = ad.AnnData(sparse.csr_matrix(x), obs=obs, var=var)
    a.layers["counts"] = sparse.csr_matrix(x.astype(np.int32))
    a.obsm["X_pca"] = np.arange(8, dtype=np.float32).reshape(4, 2) / 3
    a.obsp["connectivities"] = sparse.csr_matrix(np.array(
        [[0, .5, 0, 1], [.5, 0, 1, 0], [0, 1, 0, .25], [1, 0, .25, 0]], dtype=np.float32))
    a.varm["PCs"] = np.arange(6, dtype=np.float64).reshape(3, 2) / 7
    a.varp["gene_links"] = sparse.csc_matrix(np.eye(3, dtype=np.float64))
    raw_var = pd.DataFrame({"source": pd.Categorical(["a", "b", "a", "a", "b"], categories=index(["a", "b"]))},
                           index=index(["raw-e", "gene-a", "raw-d", "gene-z", "gene-m"], name="raw_gene_id"))
    a.raw = ad.AnnData(sparse.csr_matrix(np.arange(20, dtype=np.int32).reshape(4, 5)),
                      obs=obs.copy(), var=raw_var)
    a.uns["z_first"] = {"steps": [{"method": "fixture", "parameters": {"seed": 3}}, None],
                         "mapping": {"zeta": 2, "alpha": True}}
    a.uns["a_second"] = "opaque provenance preserved by the project codec"
    return a


def matrix_summary(value):
    return {"shape": list(value.shape), "dtype": str(value.dtype),
            "sparse": sparse.issparse(value), "format": getattr(value, "format", None),
            "values": (value.toarray() if sparse.issparse(value) else np.asarray(value)).tolist()}


def frame_summary(frame):
    return {"index": frame.index.tolist(), "index_name": frame.index.name,
            "columns": {key: {"dtype": str(col.dtype),
                              "values": [None if pd.isna(v) else v.item() if isinstance(v, np.generic) else v for v in col],
                              **({"categories": col.cat.categories.tolist(), "ordered": col.cat.ordered}
                                 if isinstance(col.dtype, pd.CategoricalDtype) else {})}
                        for key, col in frame.items()}}


def jsonable(value):
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def summary(a):
    return {"shape": list(a.shape), "X": matrix_summary(a.X),
            "obs": frame_summary(a.obs), "var": frame_summary(a.var),
            **{slot: {k: matrix_summary(getattr(a, slot)[k]) for k in getattr(a, slot).keys() if k is not None}
               for slot in ("layers", "obsm", "obsp", "varm", "varp")},
            "raw": None if a.raw is None else {"X": matrix_summary(a.raw.X), "var": frame_summary(a.raw.var)},
            "uns": jsonable(a.uns), "uns_key_order": list(a.uns)}


def differences(left, right, prefix=""):
    if isinstance(left, dict) and isinstance(right, dict):
        result = []
        for key in left.keys() | right.keys():
            at = f"{prefix}.{key}" if prefix else key
            if key not in left or key not in right:
                result.append({"path": at, "before": left.get(key), "after": right.get(key)})
            else:
                result.extend(differences(left[key], right[key], at))
        return result
    if left != right:
        return [{"path": prefix, "before": left, "after": right}]
    return []


def compare_files(source, names):
    baseline = summary(read_anndata(source))
    native_baseline = summary(ad.read_h5ad(source / "data.h5ad"))
    report = {"python": sys.version, "anndata": version("anndata"), "numpy": np.__version__,
              "pandas": pd.__version__, "baseline": baseline,
              "native_python_baseline": native_baseline, "outputs": {}}
    for name in names:
        root = HERE / name
        result = report["outputs"][name] = {}
        path = root / "data.h5ad"
        if not path.exists():
            result["missing_output"] = True
            continue
        with h5py.File(path) as handle:
            result["root_hdf5_keys"] = list(handle)
        for reader_name, reader in [("native_python", lambda: ad.read_h5ad(path)),
                                    ("project_codec", lambda: read_anndata(root))]:
            try:
                observed = summary(reader())
                expected = baseline if reader_name == "project_codec" else native_baseline
                result[reader_name] = {"summary": observed, "differences": differences(expected, observed)}
            except Exception as exc:
                result[reader_name] = {"error_type": type(exc).__name__, "error": str(exc)}
    return report


def run_variant(variant):
    prefix = "h5ad_" + variant
    source = HERE / (prefix + "_input")
    source.mkdir(exist_ok=True)
    (source / "data.h5ad").unlink(missing_ok=True)
    write_anndata(source, fixture(object_strings=variant == "object_strings"))
    env = os.environ.copy()
    env["R_HOME"] = str(R_PREFIX / "lib/R")
    env["PATH"] = os.pathsep.join([str(R_PREFIX / "Library/bin"), str(R_PREFIX / "lib/R/bin/x64"), env["PATH"]])
    names = [prefix + "_roundtrip", prefix + "_r_modified"]
    for name in names:
        (HERE / name / "data.h5ad").unlink(missing_ok=True)
    run = subprocess.run([str(R_PREFIX / "lib/R/bin/x64/Rscript.exe"), str(HERE / "h5ad_roundtrip.R"),
                          str(source / "data.h5ad"), str(HERE), prefix],
                         env=env, capture_output=True, text=True, timeout=120)
    (HERE / (prefix + "_run.log")).write_text(run.stdout + "\nSTDERR:\n" + run.stderr, encoding="utf-8")
    report = compare_files(source, names)
    report.update(r_exit=run.returncode, r_stdout=run.stdout, r_stderr=run.stderr)
    print(json.dumps({"variant": variant, "r_exit": run.returncode, "outputs": {
        name: {reader: [d["path"] for d in details["differences"]] if "differences" in details else details.get("error", details)
               for reader, details in value.items() if isinstance(details, dict)}
        for name, value in report["outputs"].items()}}, indent=2))
    return report


def main():
    if sys.argv[1:] == ["--compare-archived"]:
        report = compare_files(HERE / "h5ad_input", ["h5ad_roundtrip", "h5ad_r_modified"])
        report["provenance"] = "Original successful native-R run; see h5ad_run.log. This command only re-reads its preserved files."
        (HERE / "h5ad_initial_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        print("Compared preserved successful R outputs: h5ad_initial_report.json")
        return 0
    report = {variant: run_variant(variant) for variant in ("pandas_default", "object_strings")}
    (HERE / "h5ad_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0 if all(result["r_exit"] == 0 for result in report.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
