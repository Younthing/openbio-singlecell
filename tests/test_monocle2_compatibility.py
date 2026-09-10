"""Optional scientific parity against an unmodified Monocle 2 installation."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite

from openbio_singlecell.monocle2 import MONOCLE2_PACKAGES, run_analysis
from openbio_singlecell.r_runtime import probe_r_runtime

pytestmark = pytest.mark.skipif(
    not (os.environ.get("OPENBIO_TEST_RSCRIPT") and os.environ.get("OPENBIO_REFERENCE_RSCRIPT")),
    reason="Set OPENBIO_TEST_RSCRIPT and OPENBIO_REFERENCE_RSCRIPT for modern/native Monocle 2 parity",
)
REFERENCE = Path(__file__).with_name("monocle2_compatibility_reference.R")


def test_modern_dependencies_preserve_native_monocle2_scientific_results(tmp_path):
    modern = os.environ["OPENBIO_TEST_RSCRIPT"]
    native = os.environ["OPENBIO_REFERENCE_RSCRIPT"]
    runtime = probe_r_runtime(rscript=modern, packages=MONOCLE2_PACKAGES)
    rng = np.random.default_rng(2026)
    progression = np.tile(np.linspace(0, 1, 60), 3)
    branch = np.repeat(["progenitor", "branch_a", "branch_b"], 60)
    latent = np.column_stack((progression + np.repeat([0, 1, 1], 60),
                             (branch == "branch_a") * progression, (branch == "branch_b") * progression))
    means = np.exp(1.5 + rng.normal(0, 0.8, (400, 3)) @ latent.T)
    data = ad.AnnData(sparse.csr_matrix(rng.negative_binomial(5, 5 / (5 + means)).T.astype(float)))
    data.obs_names = [f"cell_{i}" for i in range(180)]
    data.var_names = [f"gene_{i}" for i in range(400)]
    data.obs["branch"] = branch
    data.var["gene_short_name"] = data.var_names
    mmwrite(tmp_path / "expression.mtx", sparse.coo_matrix(data.X.T), precision=17, symmetry="general")
    data.obs.to_csv(tmp_path / "obs.csv")
    data.var.to_csv(tmp_path / "var.csv")
    reference_dir = tmp_path / "native"
    reference_dir.mkdir()

    def reference_process(executable, *args):
        result = subprocess.run([executable, "--vanilla", str(REFERENCE), *map(str, args)],
                                capture_output=True, text=True, timeout=180)
        assert result.returncode == 0, result.stdout + result.stderr

    def assert_semantically_equal(actual, expected, location="state"):
        if isinstance(expected, dict):
            assert actual.keys() == expected.keys(), location
            for key in expected:
                assert_semantically_equal(actual[key], expected[key], f"{location}.{key}")
        elif isinstance(expected, list):
            assert len(actual) == len(expected), location
            for index, (left, right) in enumerate(zip(actual, expected, strict=True)):
                assert_semantically_equal(left, right, f"{location}[{index}]")
        elif isinstance(expected, (int, float)):
            assert actual == pytest.approx(expected, rel=1e-8, abs=1e-10), location
        else:
            assert actual == expected, location

    def compare_state(result, stage):
        actual_path = tmp_path / f"{stage}-modern.json"
        expected_path = tmp_path / f"{stage}-native.json"
        reference_process(modern, "inspect", result["cds_path"] / "state.rds", actual_path)
        reference_process(native, "inspect", reference_dir / f"{stage}.rds", expected_path)
        assert_semantically_equal(json.loads(actual_path.read_text()), json.loads(expected_path.read_text()))

    def run(operation, parameters, *, previous=None, stage=None):
        destination = tmp_path / (stage or operation)
        destination.mkdir()
        return run_analysis(operation, runtime, destination, parameters,
                            adata=data.copy() if operation == "create" else None,
                            cds_path=None if previous is None else previous["cds_path"])

    reference_process(native, "reference", tmp_path, reference_dir)
    created = run("create", {
        "source": {"source": "X"}, "expression_family": "negbinomial.size",
        "lowerDetectionLimit": 0.1, "estimate_size_factors": True, "estimate_dispersions": True,
        "detect_genes": True, "size_factor_parameters": {}, "detect_parameters": {},
        "dispersion_parameters": {"modelFormulaStr": "~branch", "remove_outliers": False},
    })
    compare_state(created, "create")
    selected = run("ordering_genes", {"method": "explicit", "genes": data.var_names.tolist()}, previous=created)
    reduced = run("ddrtree", {"reduction_method": "DDRTree", "max_components": 2}, previous=selected)
    compare_state(reduced, "ddrtree")
    ordered = run("order_cells", {"reverse": False}, previous=reduced)
    compare_state(ordered, "order_cells")
    root_state = (reference_dir / "root_state.txt").read_text().strip()
    rooted = run("order_cells", {"root_state": root_state, "reverse": False}, previous=ordered, stage="reroot")
    compare_state(rooted, "reroot")
    assert rooted["summary"]["root_method"] == "state"
    assert np.isfinite(rooted["table"]["Pseudotime"]).all()
    for operation, parameters in (
        ("differential_test", {"fullModelFormulaStr": "~sm.ns(Pseudotime, df=3) + branch",
                               "reducedModelFormulaStr": "~branch", "relative_expr": False}),
        ("beam", {"branch_point": 1, "progenitor_method": "duplicate"}),
    ):
        result = run(operation, {"genes": data.var_names[:8].tolist(), "cores": 1, **parameters}, previous=rooted)
        expected = pd.read_csv(reference_dir / f"{operation}.csv", index_col=0)
        actual = result["table"].loc[expected.index, expected.columns]
        assert "OK" in actual["status"].values
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False, check_index_type=False,
                                      check_categorical=False,
                                      rtol=1e-8, atol=1e-10)
