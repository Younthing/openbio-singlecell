from __future__ import annotations

import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from openbio_singlecell.monocle2 import (
    MONOCLE2_PACKAGES,
    read_cds,
    reproduction_code,
    run_analysis,
    runtime_reproduction_code,
)
from openbio_singlecell.r_runtime import probe_r_runtime

pytestmark = pytest.mark.skipif(
    not os.environ.get("OPENBIO_TEST_RSCRIPT"), reason="Set OPENBIO_TEST_RSCRIPT for native Monocle 2",
)


@pytest.fixture(scope="module")
def native_runtime():
    parameters = {
        "rscript": os.environ["OPENBIO_TEST_RSCRIPT"],
        "r_home": os.environ.get("OPENBIO_TEST_R_HOME", ""),
        "library_paths": os.environ.get("OPENBIO_TEST_R_LIBRARIES", ""),
        "path_prefix": os.environ.get("OPENBIO_TEST_R_PATH", ""),
    }
    return parameters, probe_r_runtime(**parameters, packages=MONOCLE2_PACKAGES)


def test_runtime_reproduction_probes_the_actual_selected_packages(native_runtime):
    parameters, expected = native_runtime
    namespace = {}
    exec(runtime_reproduction_code(parameters), namespace)
    assert namespace["probe_monocle2_runtime"]() == expected


def test_reproduction_runs_portable_native_chain_and_preserves_export_input(tmp_path, native_runtime):
    from PIL import Image

    _, runtime = native_runtime
    rng = np.random.default_rng(2026)
    progression = np.tile(np.linspace(0, 1, 60), 3)
    branch = np.repeat(["progenitor", "branch_a", "branch_b"], 60)
    latent = np.column_stack((progression + np.repeat([0, 1, 1], 60),
                             (branch == "branch_a") * progression, (branch == "branch_b") * progression))
    means = np.exp(1.5 + rng.normal(0, 0.8, (400, 3)) @ latent.T)
    full = ad.AnnData(sparse.csr_matrix(rng.negative_binomial(5, 5 / (5 + means)).T.astype(float)))
    full.obs_names = [f"cell_{i}" for i in range(180)]
    full.var_names = [f"gene_{i}" for i in range(400)]
    full.obs["branch"] = pd.Categorical(branch, ordered=True)
    full.layers["counts"] = full.X.copy()
    full.uns["study"] = "亚群 Monocle 2"
    full.raw = full.copy()
    original = full[:, :200].copy()

    def pair(operation, parameters, *, previous=None, adata=None):
        stem = f"{operation}-{len(list(tmp_path.iterdir()))}"
        direct_dir = tmp_path / f"{stem}-direct"
        direct_dir.mkdir()
        direct = run_analysis(operation, runtime, direct_dir, parameters,
                              adata=None if adata is None else adata.copy(),
                              cds_path=None if previous is None else previous[0]["cds_path"])
        code = reproduction_code(operation, parameters, runtime)
        assert "from openbio_singlecell" not in code
        namespace = {}
        exec(code, namespace)
        runtime_override = {} if operation == "create" else {"r_runtime": json.dumps(runtime)}
        reproduced = namespace["run_monocle2"](
            adata=adata, cds_path=None if previous is None else previous[1]["cds_path"],
            output_dir=tmp_path / f"{stem}-reproduced", **runtime_override,
        )
        assert reproduced["summary"] == direct["summary"]
        for name in ("cell_metadata.csv", "ordering_genes.csv", "embedding.csv", "table.csv"):
            if (direct_dir / name).exists():
                pd.testing.assert_frame_equal(
                    pd.read_csv(reproduced["output_dir"] / name), pd.read_csv(direct_dir / name),
                )
        if "cds_path" in direct:
            # Native RDS bytes may include serialization details; compare its sealed semantic sidecars.
            direct_metadata, reproduced_metadata = read_cds(direct["cds_path"]), read_cds(reproduced["cds_path"])
            for field in ("kind", "codec", "n_obs", "n_vars", "operation", "versions", "source"):
                assert reproduced_metadata[field] == direct_metadata[field]
        if "table" in direct:
            pd.testing.assert_frame_equal(reproduced["table"], direct["table"])
        if "plot_path" in direct:
            with Image.open(direct["plot_path"]) as expected, Image.open(reproduced["plot_path"]) as actual:
                assert actual.format == "PNG"
                assert actual.size == expected.size == (480, 240)
                # Native pseudotime plots jitter points independently in each R process.
                assert np.asarray(actual).std() > 0
        return direct, reproduced

    created = pair("create", {
        "source": {"source": "raw"}, "gene_short_name_column": "", "expression_family": "negbinomial.size",
        "lowerDetectionLimit": 0.1, "estimate_size_factors": True, "estimate_dispersions": True,
        "detect_genes": True, "size_factor_parameters": {}, "dispersion_parameters": {}, "detect_parameters": {},
    }, adata=original)
    assert created[1]["summary"]["n_vars"] == 400
    selected = pair("ordering_genes", {"method": "explicit", "genes": full.var_names.tolist()}, previous=created)
    reduced = pair("ddrtree", {"reduction_method": "DDRTree", "max_components": 2}, previous=selected)
    ordered = pair("order_cells", {"root_state": None, "reverse": False, "num_paths": None}, previous=reduced)
    root_state = ordered[0]["table"].loc[branch == "progenitor", "State"].value_counts().index[0]
    rooted = pair("order_cells", {"root_state": str(root_state), "reverse": False, "num_paths": None}, previous=ordered)
    assert np.isfinite(rooted[1]["table"]["Pseudotime"]).all()
    assert rooted[1]["summary"]["root_method"] == "state"

    pair("differential_test", {
        "genes": full.var_names[:8].tolist(), "fullModelFormulaStr": "~sm.ns(Pseudotime, df=3)",
        "reducedModelFormulaStr": "~1", "relative_expr": True, "cores": 1,
    }, previous=rooted)
    pair("gene_trends", {
        "genes": full.var_names[:2].tolist(), "color_by": "branch", "trend_formula": "~sm.ns(Pseudotime, df=3)",
        "relative_expr": True, "ncol": 2, "width": 6, "height": 3, "dpi": 80,
    }, previous=rooted)
    exported = pair("export", {
        "embedding_key": "X_monocle2", "pseudotime_key": "monocle2_pseudotime", "state_key": "monocle2_state",
        "overwrite_existing": False,
    }, previous=rooted, adata=original)
    result = exported[1]["adata"]
    pd.testing.assert_frame_equal(result.obs, exported[0]["adata"].obs)
    np.testing.assert_array_equal(result.obsm["X_monocle2"], exported[0]["adata"].obsm["X_monocle2"])
    assert result is not original
    assert (result.X != original.X).nnz == 0
    assert (result.layers["counts"] != original.layers["counts"]).nnz == 0
    assert result.raw.shape == (180, 400)
    assert (result.raw.X != original.raw.X).nnz == 0
    assert result.uns == original.uns
    pd.testing.assert_frame_equal(original.obs, full.obs)
    assert len(original.obsm) == 0
