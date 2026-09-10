from __future__ import annotations

import json
import os
import uuid

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from openbio_singlecell.artifact_codecs import read_anndata, read_table, write_anndata, write_table
from openbio_singlecell.worker_protocol import OperationContext, ProtocolError


def _context(tmp_path):
    staging = tmp_path / f"run-{uuid.uuid4()}"
    staging.mkdir()
    return OperationContext(staging, str(uuid.uuid4()))


def test_prepare_publishes_native_cds_and_preserves_expert_expression_choice(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    adata = ad.AnnData(np.array([[0.5, -0.25], [1.75, 3.0]]))
    source = tmp_path / "adata"
    source.mkdir()
    write_anndata(source, adata)
    captured = {}

    def run_analysis(operation, runtime, output_dir, parameters, *, adata=None, cds_path=None):
        captured.update(operation=operation, runtime=runtime, parameters=parameters, adata=adata)
        state = output_dir / "state.rds"
        state.write_bytes(b"native-test-state")
        return {
            "output_dir": output_dir,
            "cds_path": state,
            "summary": {"n_obs": 2, "n_vars": 2, "warnings": ["native advisory"],
                        "versions": {"R": "4.4.3", "monocle": "2.34.0"}},
        }

    monkeypatch.setattr(operations, "run_analysis", run_analysis)
    monkeypatch.setattr(operations, "reproduction_code", lambda *args: "def run_monocle2():\n    pass\n")
    runtime = {"executable": "/expert/Rscript", "r_version": "4.4.3", "packages": {}}
    parameters = {
        "r_runtime": json.dumps(runtime), "source": {"source": "X"}, "gene_short_name_column": "",
        "expression_family": "tobit", "lower_detection_limit": -0.5,
        "estimate_size_factors": False, "estimate_dispersions": False, "detect_genes": False,
        "size_factor_parameters_json": "{}", "dispersion_parameters_json": "{}",
        "detect_parameters_json": '{"min_expr": -1}',
    }
    context = _context(tmp_path)
    records = operations.monocle2_prepare(context, {
        "adata": {"type": "artifact", "path": str(source.resolve()), "kind": "OPENBIO_ANNDATA",
                  "codec": "anndata-h5ad-v1"},
    }, parameters)

    assert [record["name"] for record in records] == ["cds", "summary", "code"]
    assert records[0]["codec"] == "monocle2-cds-rds-v1"
    assert (context.output_root / records[0]["payload"] / "state.rds").read_bytes() == b"native-test-state"
    assert captured["operation"] == "create"
    assert captured["parameters"]["lowerDetectionLimit"] == -0.5
    assert captured["parameters"]["detect_parameters"] == {"min_expr": -1}
    assert np.array_equal(captured["adata"].X, adata.X)
    assert records[1]["value"]["summary"]["software_versions"]["R:monocle"] == "2.34.0"
    assert "native advisory" in records[1]["value"]["summary"]["warnings"]


def test_runtime_probes_selected_environment_and_returns_reusable_descriptor(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    captured = {}
    runtime = {"executable": "/expert/Rscript", "r_version": "4.4.3", "packages": {
        "monocle": {"version": "2.34.0", "path": "/expert/library/monocle"},
    }}

    def probe(rscript, **kwargs):
        captured.update(rscript=rscript, **kwargs)
        return runtime

    monkeypatch.setattr(operations, "probe_r_runtime", probe, raising=False)
    records = operations.monocle2_runtime(_context(tmp_path), {}, {
        "rscript": "custom-Rscript", "r_home": "", "library_paths": "custom-library", "path_prefix": "",
    })
    assert [record["name"] for record in records] == ["r_runtime", "summary", "code"]
    assert json.loads(records[0]["value"]) == runtime
    assert captured["rscript"] == "custom-Rscript"
    assert captured["library_paths"] == "custom-library"
    assert "monocle" in captured["packages"]
    assert "def probe_monocle2_runtime" in records[-1]["value"]


def _cds_inputs(tmp_path):
    root = tmp_path / "cds-input"
    root.mkdir(exist_ok=True)
    return {"cds": {"type": "artifact", "path": str(root.resolve()),
                    "kind": "OPENBIO_MONOCLE2_CDS", "codec": "monocle2-cds-rds-v1"}}


def _capture_analysis(monkeypatch, *, table=None, plot_path=None, adata=None):
    from openbio_singlecell import operations_monocle2 as operations

    captured = {}

    def run(operation, runtime, output_dir, parameters, **kwargs):
        captured.update(operation=operation, runtime=runtime, parameters=parameters, **kwargs)
        result = {"output_dir": output_dir, "summary": {
            "n_obs": 3, "n_vars": 2, "warnings": [], "versions": {"monocle": "2.34.0"},
            "effective_parameters": parameters, "ordering_gene_count": 1,
        }}
        if operation in {"ordering_genes", "ddrtree", "order_cells"}:
            (output_dir / "state.rds").write_bytes(b"native-test-state")
            result["cds_path"] = output_dir
        if operation == "order_cells":
            result["summary"]["root_method"] = "state" if parameters["root_state"] else "automatic"
        if table is not None:
            result["table"] = table
        if plot_path is not None:
            result["plot_path"] = plot_path
        if adata is not None:
            result["adata"] = adata
        return result

    monkeypatch.setattr(operations, "run_analysis", run)
    monkeypatch.setattr(operations, "reproduction_code", lambda *args: "def run_monocle2():\n    pass\n")
    return captured


def test_ordering_genes_publishes_native_state_and_selection_table(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    captured = _capture_analysis(monkeypatch, table=pd.DataFrame({
        "gene_id": ["g1", "g2"], "use_for_ordering": [True, False],
    }))
    records = operations.monocle2_ordering_genes(_context(tmp_path), _cds_inputs(tmp_path), {
        "r_runtime": "{}", "method": "explicit", "genes_json": '["g1", "unknown"]',
        "var_column": "highly_variable", "min_mean_expression": 0.0, "dispersion_fold": 0.0,
    })
    assert [record["name"] for record in records] == ["cds", "table", "summary", "code"]
    assert captured["parameters"]["genes"] == ["g1", "unknown"]
    assert captured["parameters"]["mean_expression"] == 0.0
    assert captured["cds_path"] == (tmp_path / "cds-input").resolve()


def test_ddrtree_forwards_native_controls_and_advanced_overrides(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    captured = _capture_analysis(monkeypatch)
    records = operations.monocle2_ddrtree(_context(tmp_path), _cds_inputs(tmp_path), {
        "r_runtime": "{}", "num_components": 2, "norm_method": "none", "residual_formula": "~batch",
        "pseudo_expr": 0.0, "auto_param_selection": False, "scaling": False,
        "extra_parameters_json": '{"max_components": 3, "ncenter": 20}',
    })
    assert [record["name"] for record in records] == ["cds", "summary", "code"]
    assert captured["parameters"]["reduction_method"] == "DDRTree"
    assert captured["parameters"]["residualModelFormulaStr"] == "~batch"
    assert captured["parameters"]["extra_parameters"] == {"max_components": 3, "ncenter": 20}
    assert "random_seed" not in captured["parameters"]


@pytest.mark.parametrize("root,expected", [({"root": "automatic"}, None), ({"root": "state", "state": "3"}, "3")])
def test_order_cells_preserves_native_initial_and_explicit_root_choices(tmp_path, monkeypatch, root, expected):
    from openbio_singlecell import operations_monocle2 as operations

    captured = _capture_analysis(monkeypatch, table=pd.DataFrame({"cell_id": ["c1"], "Pseudotime": [0.0], "State": ["3"]}))
    records = operations.monocle2_order_cells(_context(tmp_path), _cds_inputs(tmp_path), {
        "r_runtime": "{}", "root": root, "reverse": True,
    })
    assert [record["name"] for record in records] == ["cds", "table", "summary", "code"]
    assert captured["parameters"] == {"root_state": expected, "reverse": True}
    if expected is None:
        assert "unoriented" in records[-2]["value"]["summary"]["results"].lower()


def test_explicit_root_mode_explains_the_missing_state_choice(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    _capture_analysis(monkeypatch, table=pd.DataFrame({"cell_id": ["c1"]}))
    with pytest.raises(ValueError, match="Choose a root State"):
        operations.monocle2_order_cells(_context(tmp_path), _cds_inputs(tmp_path), {
            "r_runtime": "{}", "root": {"root": "state", "state": ""}, "reverse": False,
        })


def _png(tmp_path):
    from PIL import Image

    path = tmp_path / "native-plot.png"
    Image.new("RGB", (2, 2), "white").save(path)
    return path


def test_trajectory_plot_publishes_native_png_with_expert_plot_options(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    captured = _capture_analysis(monkeypatch, plot_path=_png(tmp_path))
    records = operations.monocle2_trajectory_plot(_context(tmp_path), _cds_inputs(tmp_path), {
        "r_runtime": "{}", "color_by": "sample", "x_dimension": 1, "y_dimension": 3,
        "show_tree": False, "show_branch_points": False, "show_state_number": False, "cell_size": 0.2,
        "width": 8, "height": 6, "dpi": 150, "extra_parameters_json": '{"theta": 90}',
    })
    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert captured["parameters"]["y"] == 3
    assert captured["parameters"]["show_tree"] is False
    assert captured["parameters"]["extra_parameters"] == {"theta": 90}


def test_differential_test_preserves_failed_fits_and_user_formulas(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    captured = _capture_analysis(monkeypatch, table=pd.DataFrame({
        "gene_id": ["g1", "g2"], "status": ["OK", "FAIL"], "pval": [0.01, np.nan], "qval": [0.02, np.nan],
    }))
    context = _context(tmp_path)
    records = operations.monocle2_differential_test(context, _cds_inputs(tmp_path), {
        "r_runtime": "{}", "genes_json": "[]", "full_formula": "~sm.ns(Pseudotime, df=4)+batch",
        "reduced_formula": "~batch", "relative_expr": True, "cores": 1, "extra_parameters_json": "{}",
    })
    table, _metadata = read_table(context.output_root / records[0]["payload"])
    assert [record["name"] for record in records] == ["table", "summary", "code"]
    assert table["status"].tolist() == ["OK", "FAIL"]
    assert pd.isna(table.loc[1, "pval"])
    assert captured["parameters"]["fullModelFormulaStr"] == "~sm.ns(Pseudotime, df=4)+batch"
    assert captured["parameters"]["reducedModelFormulaStr"] == "~batch"
    assert captured["parameters"]["genes"] == []


def test_beam_forwards_explicit_branch_and_progenitor_method(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    captured = _capture_analysis(monkeypatch, table=pd.DataFrame({
        "gene_id": ["g1"], "status": ["OK"], "pval": [0.01], "qval": [0.01],
    }))
    records = operations.monocle2_beam(_context(tmp_path), _cds_inputs(tmp_path), {
        "r_runtime": "{}", "genes_json": '["g1"]', "branch_point": 2, "progenitor_method": "duplicate",
        "cores": 1, "extra_parameters_json": '{"relative_expr": false}',
    })
    assert [record["name"] for record in records] == ["table", "summary", "code"]
    assert captured["parameters"]["branch_point"] == 2
    assert captured["parameters"]["progenitor_method"] == "duplicate"
    assert captured["parameters"]["extra_parameters"] == {"relative_expr": False}


@pytest.mark.parametrize("explicit,expected", [("[]", ["g2", "g3"]), ('["chosen"]', ["chosen"])])
def test_gene_trends_uses_visible_table_top_n_or_explicit_genes(tmp_path, monkeypatch, explicit, expected):
    from openbio_singlecell import operations_monocle2 as operations

    table_root = tmp_path / "gene-table"
    table_root.mkdir()
    write_table(table_root, pd.DataFrame({
        "gene_id": ["g1", "g2", "g3"], "qval": [0.3, 0.01, 0.02], "pval": [0.1, 0.001, 0.002],
    }), {})
    inputs = _cds_inputs(tmp_path)
    inputs["table"] = {"type": "artifact", "path": str(table_root.resolve()),
                       "kind": "OPENBIO_SINGLE_CELL_TABLE", "codec": "table-jsonl-v1"}
    captured = _capture_analysis(monkeypatch, plot_path=_png(tmp_path))
    records = operations.monocle2_gene_trends(_context(tmp_path), inputs, {
        "r_runtime": "{}", "genes_json": explicit, "top_n": 2, "color_by": "State",
        "trend_formula": "~sm.ns(Pseudotime, df=4)", "relative_expr": True, "ncol": 2,
        "width": 10, "height": 7, "dpi": 150, "extra_parameters_json": "{}",
    })
    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    assert captured["parameters"]["genes"] == expected
    assert captured["parameters"]["trend_formula"] == "~sm.ns(Pseudotime, df=4)"


def test_export_publishes_owned_anndata_and_aligned_table(tmp_path, monkeypatch):
    from openbio_singlecell import operations_monocle2 as operations

    adata = ad.AnnData(np.array([[0.5, 2.0], [3.0, 4.5], [1.0, 3.0]]))
    adata.obs_names = ["c1", "c2", "c3"]
    source = tmp_path / "adata"
    source.mkdir()
    write_anndata(source, adata)
    original = (source / "data.h5ad").read_bytes()
    exported = adata.copy()
    exported.obs["monocle2_pseudotime"] = [0.0, 1.0, 2.0]
    captured = _capture_analysis(monkeypatch, adata=exported, table=pd.DataFrame({
        "cell_id": adata.obs_names, "Pseudotime": [0.0, 1.0, 2.0],
    }))
    inputs = _cds_inputs(tmp_path)
    inputs["adata"] = {"type": "artifact", "path": str(source.resolve()),
                       "kind": "OPENBIO_ANNDATA", "codec": "anndata-h5ad-v1"}
    context = _context(tmp_path)
    records = operations.monocle2_export(context, inputs, {
        "r_runtime": "{}", "embedding_key": "X_monocle2", "pseudotime_key": "monocle2_pseudotime",
        "state_key": "monocle2_state", "overwrite_existing": True,
    })
    actual = read_anndata(context.output_root / records[0]["payload"])
    assert [record["name"] for record in records] == ["adata", "table", "summary", "code"]
    assert np.array_equal(actual.X, adata.X)
    assert actual.obs["monocle2_pseudotime"].tolist() == [0.0, 1.0, 2.0]
    assert captured["parameters"]["overwrite_existing"] is True
    assert (source / "data.h5ad").read_bytes() == original


@pytest.mark.parametrize("extra", ['{"max_components": NaN}', '{"max_components": 2, "max_components": 3}'])
def test_advanced_parameters_use_the_workers_strict_json_boundary(tmp_path, monkeypatch, extra):
    from openbio_singlecell import operations_monocle2 as operations

    captured = _capture_analysis(monkeypatch)
    with pytest.raises(ProtocolError):
        operations.monocle2_ddrtree(_context(tmp_path), _cds_inputs(tmp_path), {
            "r_runtime": "{}", "num_components": 2, "norm_method": "log", "residual_formula": "",
            "pseudo_expr": 1.0, "auto_param_selection": True, "scaling": True,
            "extra_parameters_json": extra,
        })
    assert captured == {}


@pytest.mark.skipif(not os.environ.get("OPENBIO_TEST_RSCRIPT"), reason="Set OPENBIO_TEST_RSCRIPT for native Monocle 2")
def test_native_worker_operations_connect_r_trajectory_and_python_artifacts(tmp_path):
    from scipy import sparse

    from openbio_singlecell import operations_monocle2 as operations

    rng = np.random.default_rng(2026)
    n_cells, n_genes = 180, 400
    progression = np.tile(np.linspace(0, 1, 60), 3)
    branch = np.repeat(["progenitor", "branch_a", "branch_b"], 60)
    latent = np.column_stack((progression + np.repeat([0, 1, 1], 60),
                             (branch == "branch_a") * progression, (branch == "branch_b") * progression))
    means = np.exp(1.5 + rng.normal(0, 0.8, (n_genes, 3)) @ latent.T)
    counts = rng.negative_binomial(5, 5 / (5 + means)).T.astype(float)
    full = ad.AnnData(sparse.csr_matrix(counts))
    full.obs_names = [f"cell_{i}" for i in range(n_cells)]
    full.var_names = [f"gene_{i}" for i in range(n_genes)]
    full.obs["branch"] = pd.Categorical(branch)
    full.layers["counts"] = full.X.copy()
    full.raw = full.copy()
    adata = full[:, :200].copy()
    original = tmp_path / "original"
    original.mkdir()
    write_anndata(original, adata)
    source_hash = (original / "data.h5ad").read_bytes()
    adata_input = {"type": "artifact", "path": str(original.resolve()),
                   "kind": "OPENBIO_ANNDATA", "codec": "anndata-h5ad-v1"}

    runtime = operations.monocle2_runtime(_context(tmp_path), {}, {
        "rscript": os.environ["OPENBIO_TEST_RSCRIPT"], "r_home": os.environ.get("OPENBIO_TEST_R_HOME", ""),
        "library_paths": os.environ.get("OPENBIO_TEST_R_LIBRARIES", ""),
        "path_prefix": os.environ.get("OPENBIO_TEST_R_PATH", ""),
    })[0]["value"]

    def run(function, inputs, parameters):
        context = _context(tmp_path)
        records = function(context, inputs, {"r_runtime": runtime, **parameters})
        artifacts = {record["name"]: {"type": "artifact", "kind": record["kind"], "codec": record["codec"],
                     "path": str((context.output_root / record["payload"]).resolve())}
                     for record in records if record["type"] == "artifact"}
        return artifacts, records

    outputs, _ = run(operations.monocle2_prepare, {"adata": adata_input}, {
        "source": {"source": "raw"}, "gene_short_name_column": "", "expression_family": "negbinomial.size",
        "lower_detection_limit": 0.1, "estimate_size_factors": True, "estimate_dispersions": True,
        "detect_genes": True, "size_factor_parameters_json": "{}", "dispersion_parameters_json": "{}",
        "detect_parameters_json": "{}",
    })
    outputs, _ = run(operations.monocle2_ordering_genes, {"cds": outputs["cds"]}, {
        "method": "dispersion", "genes_json": "[]", "var_column": "highly_variable",
        "min_mean_expression": 0.0, "dispersion_fold": 0.0,
    })
    outputs, _ = run(operations.monocle2_ddrtree, {"cds": outputs["cds"]}, {
        "num_components": 2, "norm_method": "log", "residual_formula": "", "pseudo_expr": 1.0,
        "auto_param_selection": True, "scaling": True, "extra_parameters_json": "{}",
    })
    outputs, records = run(operations.monocle2_order_cells, {"cds": outputs["cds"]}, {
        "root": {"root": "automatic"}, "reverse": False,
    })
    assert "unoriented" in records[-2]["value"]["summary"]["results"]
    cells, _ = read_table(outputs["table"]["path"])
    root_state = cells.loc[branch == "progenitor", "State"].value_counts().index[0]
    outputs, _ = run(operations.monocle2_order_cells, {"cds": outputs["cds"]}, {
        "root": {"root": "state", "state": str(root_state)}, "reverse": False,
    })
    rooted = outputs["cds"]
    plot, _ = run(operations.monocle2_trajectory_plot, {"cds": rooted}, {
        "color_by": "Pseudotime", "x_dimension": 1, "y_dimension": 2, "show_tree": True,
        "show_branch_points": True, "show_state_number": True, "cell_size": 1.5,
        "width": 6, "height": 4, "dpi": 100, "extra_parameters_json": "{}",
    })
    assert "plot" in plot
    evidence, _ = run(operations.monocle2_differential_test, {"cds": rooted}, {
        "genes_json": json.dumps(full.var_names[:12].tolist()), "full_formula": "~sm.ns(Pseudotime, df=3)",
        "reduced_formula": "~1", "relative_expr": True, "cores": 1, "extra_parameters_json": "{}",
    })
    table, _ = read_table(evidence["table"]["path"])
    assert len(table) == 12
    beam, _ = run(operations.monocle2_beam, {"cds": rooted}, {
        "genes_json": json.dumps(full.var_names[:12].tolist()), "branch_point": 1,
        "progenitor_method": "duplicate", "cores": 1, "extra_parameters_json": "{}",
    })
    branch_table, _ = read_table(beam["table"]["path"])
    assert len(branch_table) == 12
    trends, _ = run(operations.monocle2_gene_trends, {"cds": rooted, "table": evidence["table"]}, {
        "genes_json": "[]", "top_n": 3, "color_by": "branch", "trend_formula": "~sm.ns(Pseudotime, df=3)",
        "relative_expr": True, "ncol": 3, "width": 8, "height": 3, "dpi": 100, "extra_parameters_json": "{}",
    })
    assert "plot" in trends
    exported, _ = run(operations.monocle2_export, {"cds": rooted, "adata": adata_input}, {
        "embedding_key": "X_monocle2", "pseudotime_key": "monocle2_pseudotime", "state_key": "monocle2_state",
        "overwrite_existing": False,
    })
    result = read_anndata(exported["adata"]["path"])
    assert result.shape == (n_cells, 200)
    assert result.raw.shape == (n_cells, n_genes)
    assert (result.X != adata.X).nnz == 0
    assert (result.layers["counts"] != adata.layers["counts"]).nnz == 0
    assert np.isfinite(result.obs["monocle2_pseudotime"]).all()
    assert result.obsm["X_monocle2"].shape == (n_cells, 2)
    assert (original / "data.h5ad").read_bytes() == source_hash
