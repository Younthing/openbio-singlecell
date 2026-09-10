from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from openbio_singlecell import nodes_monocle2 as nodes
from openbio_singlecell.nodes_monocle2 import OpenBioSingleCellMonocle2Prepare


def test_prepare_exposes_native_family_and_explicit_expression_processing():
    schema = OpenBioSingleCellMonocle2Prepare.define_schema()
    inputs = {item.id: item for item in schema.inputs}

    assert [(option.key, [item.id for item in option.inputs]) for option in inputs["source"].options] == [
        ("X", []), ("raw", []), ("layer", ["layer_name"]),
    ]
    assert inputs["expression_family"].options == ["negbinomial.size", "negbinomial", "tobit", "gaussianff"]
    assert inputs["gene_short_name_column"].default == ""
    assert inputs["lower_detection_limit"].default == 0.1
    assert inputs["lower_detection_limit"].min is None
    assert inputs["lower_detection_limit"].max is None
    for name in ("estimate_size_factors", "estimate_dispersions", "detect_genes"):
        assert inputs[name].default is True
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("cds", "OPENBIO_MONOCLE2_CDS"),
        ("summary", "OPENBIO_SINGLE_CELL_SUMMARY"),
        ("code", "STRING"),
    ]


def test_runtime_and_ordering_genes_keep_environment_and_gene_selection_explicit():
    runtime = nodes.OpenBioSingleCellMonocle2Runtime.define_schema()
    runtime_inputs = {item.id: item for item in runtime.inputs}
    assert runtime_inputs["rscript"].default == "Rscript"
    assert all(runtime_inputs[name].advanced for name in ("r_home", "library_paths", "path_prefix"))
    assert runtime.outputs[0].io_type == "STRING"
    assert runtime.outputs[0].display_name == "r_runtime"

    ordering = nodes.OpenBioSingleCellMonocle2OrderingGenes.define_schema()
    inputs = {item.id: item for item in ordering.inputs}
    assert inputs["method"].options == ["dispersion", "explicit", "var_column"]
    assert inputs["genes_json"].default == "[]"
    assert inputs["var_column"].default == "highly_variable"
    assert inputs["min_mean_expression"].default == 0.5
    assert inputs["dispersion_fold"].default == 1.0
    assert inputs["min_mean_expression"].min is None
    assert inputs["dispersion_fold"].max is None
    assert [item.display_name for item in ordering.outputs] == ["cds", "table", "summary", "code"]


def test_ddrtree_and_order_cells_expose_native_controls_without_selecting_a_biological_root():
    ddrtree = {item.id: item for item in nodes.OpenBioSingleCellMonocle2DDRTree.define_schema().inputs}
    assert ddrtree["norm_method"].options == ["log", "vstExprs", "none"]
    assert ddrtree["num_components"].default == 2
    assert ddrtree["num_components"].max is None
    assert ddrtree["extra_parameters_json"].default == "{}"
    assert "random_seed" not in ddrtree  # Native reduceDimension owns its fixed seed.
    assert "relative_expr" not in ddrtree  # Native reduceDimension ignores this formal argument.
    order = {item.id: item for item in nodes.OpenBioSingleCellMonocle2OrderCells.define_schema().inputs}
    assert [(option.key, [item.id for item in option.inputs]) for option in order["root"].options] == [
        ("automatic", []), ("state", ["state"]),
    ]
    assert order["root"].options[1].inputs[0].default == ""
    assert order["reverse"].default is False
    assert "num_paths" not in order  # Monocle 2 only applies this argument to ICA.


def test_gene_evidence_nodes_expose_native_formulas_and_visible_gene_selection():
    differential = {item.id: item for item in nodes.OpenBioSingleCellMonocle2DifferentialTest.define_schema().inputs}
    assert differential["full_formula"].default == "~sm.ns(Pseudotime, df=3)"
    assert differential["reduced_formula"].default == "~1"
    assert differential["genes_json"].default == "[]"
    assert differential["cores"].max is None
    beam = {item.id: item for item in nodes.OpenBioSingleCellMonocle2BEAM.define_schema().inputs}
    assert beam["branch_point"].default == 1
    assert beam["progenitor_method"].options == ["duplicate", "sequential_split"]
    trends = {item.id: item for item in nodes.OpenBioSingleCellMonocle2GeneTrends.define_schema().inputs}
    assert trends["table"].optional is True
    assert trends["genes_json"].default == "[]"
    assert trends["top_n"].default == 6
    assert trends["trend_formula"].default == "~ sm.ns(Pseudotime, df=3)"
    assert trends["top_n"].max is None


def test_trajectory_plot_and_export_expose_display_and_result_names():
    plot = {item.id: item for item in nodes.OpenBioSingleCellMonocle2TrajectoryPlot.define_schema().inputs}
    assert plot["color_by"].default == "State"
    assert plot["show_branch_points"].default is True
    assert plot["show_state_number"].default is False
    assert plot["cell_size"].max is None
    assert plot["width"].min is None
    export = {item.id: item for item in nodes.OpenBioSingleCellMonocle2Export.define_schema().inputs}
    assert export["embedding_key"].default == "X_monocle2"
    assert export["pseudotime_key"].default == "monocle2_pseudotime"
    assert export["state_key"].default == "monocle2_state"
    assert export["overwrite_existing"].default is False


def test_analysis_cache_tracks_r_identity_and_packaged_scripts(tmp_path, monkeypatch):
    executable = Path(sys.executable)
    runtime = {
        "executable": str(executable), "size": executable.stat().st_size,
        "mtime_ns": executable.stat().st_mtime_ns, "r_home": "", "library_paths": "", "path_prefix": "",
        "r_version": "4.4.3", "packages": {},
    }
    script = tmp_path / "monocle2.R"
    script.write_text("original", encoding="utf-8")
    monkeypatch.setattr(nodes, "MONOCLE2_SCRIPT", script)
    fingerprint = nodes.OpenBioSingleCellMonocle2Prepare.fingerprint_inputs(r_runtime=json.dumps(runtime))
    assert fingerprint != nodes.OpenBioSingleCellMonocle2Prepare.fingerprint_inputs(
        r_runtime=json.dumps({**runtime, "r_version": "4.5.0"}),
    )
    script.write_text("updated", encoding="utf-8")
    assert fingerprint != nodes.OpenBioSingleCellMonocle2Prepare.fingerprint_inputs(r_runtime=json.dumps(runtime))
    adapter = tmp_path / "monocle2_compat.R"
    adapter.write_text("original adapter", encoding="utf-8")
    monkeypatch.setattr(nodes, "MONOCLE2_COMPAT_SCRIPT", adapter, raising=False)
    fingerprint = nodes.OpenBioSingleCellMonocle2Prepare.fingerprint_inputs(r_runtime=json.dumps(runtime))
    adapter.write_text("updated adapter", encoding="utf-8")
    assert fingerprint != nodes.OpenBioSingleCellMonocle2Prepare.fingerprint_inputs(r_runtime=json.dumps(runtime))
    with pytest.raises(ValueError, match="changed"):
        nodes.OpenBioSingleCellMonocle2Prepare.fingerprint_inputs(r_runtime=json.dumps({**runtime, "size": 0}))
    for node in nodes.MONOCLE2_NODE_CLASSES:
        assert "fingerprint_inputs" in node.__dict__  # Adapter intentionally accepts direct overrides.
