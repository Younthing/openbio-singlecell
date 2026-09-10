from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from openbio_singlecell.artifact_codecs import read_table, write_table
from openbio_singlecell.monocle2 import MONOCLE2_CODEC, MONOCLE2_KIND, read_cds, reproduction_code, run_analysis


def _runtime():
    executable = os.path.normcase(os.path.realpath(sys.executable))
    stat = Path(executable).stat()
    return {"executable": executable, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "r_home": "", "library_paths": "", "path_prefix": "", "r_version": "test", "packages": {}}


def _cds(tmp_path):
    root = tmp_path / "input-cds"
    root.mkdir()
    files = {"state.rds": b"R output placeholder at the subprocess boundary",
             "cell_metadata.csv": b"cell_id\nc0\nc1\nc2\n",
             "ordering_genes.csv": b"gene_id\ng0\ng1\n", "summary.json": b"{}",
             "cell_metadata.json": b'{"index":["c0","c1","c2"],"columns":{"cell_id":{"type":"string","values":["c0","c1","c2"]}}}',
             "ordering_genes.json": b'{"index":["g0","g1"],"columns":{"gene_id":{"type":"string","values":["g0","g1"]}}}'}
    for name, data in files.items():
        (root / name).write_bytes(data)
    (root / "metadata.json").write_text(json.dumps({
        "schema_version": 1, "kind": MONOCLE2_KIND, "codec": MONOCLE2_CODEC,
        "files": {name: hashlib.sha256(data).hexdigest() for name, data in files.items()},
        "n_obs": 3, "n_vars": 2, "source": {"expression": {"source": "raw"}},
    }), encoding="utf-8")
    return root


def _native_export(argv, **kwargs):
    request = json.loads(Path(argv[-1]).read_text(encoding="utf-8"))
    root = Path(request["output_dir"])
    (root / "summary.json").write_text(json.dumps({
        "operation": "export", "n_obs": 3, "n_vars": 2, "versions": {"monocle": "2.34.0"},
        "state_counts": {"1": 2, "2": 1}, "warnings": [],
    }), encoding="utf-8")
    (root / "cell_metadata.csv").write_text(
        "cell_id,Pseudotime,State\nc2,2.5,2\nc0,0,1\nc1,1.25,1\n", encoding="utf-8")
    (root / "embedding.csv").write_text(
        "cell_id,DDRTree1,DDRTree2\nc2,20,21\nc0,0,1\nc1,10,11\n", encoding="utf-8")
    (root / "ordering_genes.csv").write_text("gene_id\ng0\ng1\n", encoding="utf-8")
    (root / "cell_metadata.json").write_text(json.dumps({"index": ["c2", "c0", "c1"], "columns": {
        "cell_id": {"type": "string", "values": ["c2", "c0", "c1"]},
        "Pseudotime": {"type": "number", "values": [2.5, 0., 1.25]},
        "State": {"type": "category", "values": ["2", "1", "1"], "levels": ["1", "2"], "ordered": False},
    }}), encoding="utf-8")
    (root / "embedding.json").write_text(json.dumps({"index": ["c2", "c0", "c1"], "columns": {
        "cell_id": {"type": "string", "values": ["c2", "c0", "c1"]},
        "DDRTree1": {"type": "number", "values": [20., 0., 10.]},
        "DDRTree2": {"type": "number", "values": [21., 1., 11.]},
    }}), encoding="utf-8")
    return subprocess.CompletedProcess(argv, 0)


def test_export_aligns_results_and_preserves_original_expression_and_metadata(tmp_path, monkeypatch):
    original = ad.AnnData(sparse.csr_matrix([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]),
                         obs=pd.DataFrame({"label": pd.Categorical(["a", "b", None], ordered=True)},
                                          index=["c0", "c1", "c2"]))
    original.raw = original
    original.layers["counts"] = original.X.copy()
    original.uns["kept"] = {"value": "unchanged"}
    cds_path = _cds(tmp_path)
    monkeypatch.setattr(subprocess, "run", _native_export)
    parameters = {"embedding_key": "X_monocle2", "pseudotime_key": "monocle2_pseudotime",
                  "state_key": "monocle2_state", "overwrite_existing": False}
    destination = tmp_path / "output"
    destination.mkdir()
    result = run_analysis("export", _runtime(), destination, parameters, adata=original.copy(), cds_path=cds_path)
    output = result["adata"]
    assert output.obs["monocle2_pseudotime"].tolist() == [0., 1.25, 2.5]
    assert output.obs["monocle2_state"].tolist() == ["1", "1", "2"]
    np.testing.assert_array_equal(output.obsm["X_monocle2"], [[0, 1], [10, 11], [20, 21]])
    np.testing.assert_array_equal(output.X.toarray(), original.X.toarray())
    np.testing.assert_array_equal(output.raw.X.toarray(), original.raw.X.toarray())
    np.testing.assert_array_equal(output.layers["counts"].toarray(), original.layers["counts"].toarray())
    pd.testing.assert_series_equal(output.obs["label"], original.obs["label"])
    assert output.uns["kept"] == original.uns["kept"]
    code = reproduction_code("export", parameters, _runtime())
    assert "from openbio_singlecell" not in code
    namespace = {}
    exec(code, namespace)
    reproduced = namespace["run_monocle2"](original, cds_path, output_dir=tmp_path / "reproduced")
    pd.testing.assert_frame_equal(reproduced["adata"].obs, output.obs)
    assert "monocle2_pseudotime" not in original.obs


def test_cds_artifact_rejects_changed_native_payload(tmp_path):
    root = _cds(tmp_path)
    assert read_cds(root)["n_obs"] == 3
    (root / "state.rds").write_bytes(b"changed native state")
    with pytest.raises(ValueError, match="recorded content"):
        read_cds(root)


def test_export_requires_explicit_overwrite_and_matching_cell_axes(tmp_path, monkeypatch):
    adata = ad.AnnData(np.ones((3, 2)), obs=pd.DataFrame({"time": [99., 99., 99.]}, index=["c0", "c1", "c2"]))
    cds_path = _cds(tmp_path)
    monkeypatch.setattr(subprocess, "run", _native_export)
    parameters = {"embedding_key": "X_r", "pseudotime_key": "time", "state_key": "state",
                  "overwrite_existing": False}
    blocked = tmp_path / "blocked"
    blocked.mkdir()
    with pytest.raises(ValueError, match="overwrite_existing"):
        run_analysis("export", _runtime(), blocked, parameters, adata=adata, cds_path=cds_path)
    assert adata.obs["time"].tolist() == [99., 99., 99.]
    allowed = tmp_path / "allowed"
    allowed.mkdir()
    run_analysis("export", _runtime(), allowed, {**parameters, "overwrite_existing": True},
                 adata=adata, cds_path=cds_path)
    assert adata.obs["time"].tolist() == [0., 1.25, 2.5]
    mismatched = adata.copy()
    mismatched.obs_names = ["c0", "c1", "different"]
    wrong = tmp_path / "wrong"
    wrong.mkdir()
    with pytest.raises(ValueError, match="same cells"):
        run_analysis("export", _runtime(), wrong, {**parameters, "overwrite_existing": True},
                     adata=mismatched, cds_path=cds_path)


def test_native_evidence_preserves_labels_types_and_missing_metadata(tmp_path, monkeypatch):
    def run(argv, **kwargs):
        request = json.loads(Path(argv[-1]).read_text(encoding="utf-8"))
        root = Path(request["output_dir"])
        (root / "summary.json").write_text('{"operation":"differential_test"}', encoding="utf-8")
        (root / "table.csv").write_text(
            'gene_id,gene_short_name,score,feature_tag,status,pval,qval\n001,001,5,NA,OK,0.1,0.2\nNA,002,,,FAIL,,\n',
            encoding="utf-8",
        )
        (root / "table.json").write_text(json.dumps({"index": ["001", "NA"], "columns": {
            "gene_id": {"type": "string", "values": ["001", "NA"]},
            "gene_short_name": {"type": "string", "values": ["001", "002"]},
            "score": {"type": "number", "values": [5., None]},
            "feature_tag": {"type": "string", "values": ["NA", ""]},
            "status": {"type": "string", "values": ["OK", "FAIL"]},
            "pval": {"type": "number", "values": [0.1, None]},
            "qval": {"type": "number", "values": [0.2, None]},
        }}), encoding="utf-8")
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", run)
    cds_path = _cds(tmp_path)
    destination = tmp_path / "evidence"
    destination.mkdir()
    table = run_analysis("differential_test", _runtime(), destination, {}, cds_path=cds_path)["table"]
    assert table["gene_short_name"].tolist() == ["001", "002"]
    assert table["gene_id"].tolist() == ["001", "NA"]
    assert table["feature_tag"].tolist() == ["NA", ""]
    assert pd.api.types.is_float_dtype(table["score"])
    assert table["score"].iloc[0] == 5 and pd.isna(table["score"].iloc[1])
    assert table["status"].iloc[1] == "FAIL" and pd.isna(table["qval"].iloc[1])
    namespace = {}
    exec(reproduction_code("differential_test", {}, _runtime()), namespace)
    reproduced = namespace["run_monocle2"](cds_path=cds_path, output_dir=tmp_path / "code-result")
    pd.testing.assert_frame_equal(reproduced["table"], table)


def test_native_infinite_metadata_is_numeric_after_table_publication(tmp_path):
    table = pd.DataFrame({"native_numeric": pd.array([np.inf, -np.inf, None], dtype="Float64"),
                          "label": pd.array(["Inf", "001", "NA"], dtype="string")})
    write_table(tmp_path, table, {"kind": "table"})
    restored, _ = read_table(tmp_path)
    pd.testing.assert_frame_equal(restored, table)
    for line in (tmp_path / "data.jsonl").read_text(encoding="utf-8").splitlines():
        json.loads(line, parse_constant=lambda value: pytest.fail(f"Non-JSON numeric token: {value}"))
