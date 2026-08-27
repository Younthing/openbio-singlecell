from __future__ import annotations

import json
from pathlib import Path

from install import main as manager_install
from openbio_singlecell.extension import NODE_CLASSES
from scripts.generate_demo import KNOWN_MARKERS, build_demo, validate_existing
from tests.workflow_helpers import workflow_execute_kwargs

PLUGIN_ROOT = Path(__file__).resolve().parents[1]


def output_value(node_output):
    if len(node_output.result) == 1:
        return node_output.result[0]
    return node_output.result


def template_runner(template_name):
    workflow_path = PLUGIN_ROOT / "example_workflows" / f"{template_name}.json"
    workflow = json.loads(workflow_path.read_text(encoding="utf-8"))
    workflow_nodes = {node["type"]: node for node in workflow["nodes"]}
    node_classes = {node.GET_SCHEMA().node_id: node for node in NODE_CLASSES}

    def run(node_type, *input_values):
        node_class = node_classes[node_type]
        kwargs = workflow_execute_kwargs(node_class, workflow_nodes[node_type])
        return output_value(node_class.execute(*input_values, **kwargs))

    return run


def write_demo(input_dir, science):
    demo_path = input_dir / "openbio-singlecell" / "openbio_singlecell_demo.h5ad"
    demo_path.parent.mkdir()
    build_demo(science.np, science.pd, science.sparse, science.ad).write_h5ad(demo_path)


def test_demo_conditions_are_biological_cohorts_balanced_across_technical_batches(science):
    demo = build_demo(science.np, science.pd, science.sparse, science.ad)

    sample_metadata = demo.obs[["sample", "condition", "batch"]].astype(str).drop_duplicates()
    assert set(sample_metadata.itertuples(index=False, name=None)) == {
        ("sample_1", "control", "batch_1"),
        ("sample_2", "treated", "batch_1"),
        ("sample_3", "control", "batch_2"),
        ("sample_4", "treated", "batch_2"),
    }


def test_manager_install_generates_demo(tmp_path, monkeypatch, science):
    comfy_root = tmp_path / "ComfyUI"
    comfy_root.mkdir()
    (comfy_root / "main.py").touch()
    monkeypatch.delenv("OPENBIO_COMFYUI_ROOT", raising=False)
    monkeypatch.setenv("COMFYUI_FOLDERS_BASE_PATH", str(comfy_root))

    assert manager_install() == 0

    output = comfy_root / "input" / "openbio-singlecell" / "openbio_singlecell_demo.h5ad"
    valid, reason = validate_existing(output, science.ad, science.sparse)
    assert valid, reason
    persisted = science.ad.read_h5ad(output)
    history = persisted.uns["openbio_singlecell"]["analysis_history"]
    assert isinstance(history, dict)
    assert history == {}


def test_demo_validation_rejects_forged_same_shape_anndata(tmp_path, science):
    demo = build_demo(science.np, science.pd, science.sparse, science.ad)
    valid_path = tmp_path / "valid.h5ad"
    demo.write_h5ad(valid_path)
    assert demo.uns["openbio_singlecell"]["display_name"] == "OpenBio single-cell demo"
    assert demo.uns["openbio_singlecell"]["analysis_history"] == {}
    assert validate_existing(valid_path, science.ad, science.sparse) == (True, "valid")

    cases = [
        ("markers", "known marker genes"),
        ("groups", "exactly 200"),
        ("samples", "samples do not match"),
        ("conditions", "conditions do not match"),
        ("batches", "batches do not match"),
        ("fractional", "integer counts"),
        ("negative", "negative counts"),
        ("metadata", "metadata must be a mapping"),
        ("schema", "unsupported schema_version"),
        ("history", "analysis_history must be a mapping"),
    ]
    for issue, expected_reason in cases:
        forged = demo.copy()
        if issue == "markers":
            names = list(forged.var_names)
            names[names.index(next(iter(KNOWN_MARKERS["T cell"])))] = "FORGED_MARKER"
            forged.var_names = names
        elif issue == "groups":
            groups = forged.obs["cell_type"].astype(str).to_numpy().copy()
            groups[0] = next(group for group in KNOWN_MARKERS if group != groups[0])
            forged.obs["cell_type"] = science.pd.Categorical(groups, categories=list(KNOWN_MARKERS))
        elif issue == "samples":
            forged.obs["sample"] = science.pd.Categorical(["sample_1"] * forged.n_obs)
        elif issue == "conditions":
            forged.obs["condition"] = science.pd.Categorical(["control"] * forged.n_obs)
        elif issue == "batches":
            forged.obs["batch"] = science.pd.Categorical(["batch_1"] * forged.n_obs)
        elif issue == "fractional":
            forged.X = forged.X.astype(float)
            forged.X.data[0] = 0.5
        elif issue == "negative":
            forged.X = forged.X.copy()
            forged.X.data[0] = -1
        elif issue == "metadata":
            forged.uns["openbio_singlecell"] = science.np.asarray([], dtype=float)
        elif issue == "schema":
            forged.uns["openbio_singlecell"]["schema_version"] = 999
        elif issue == "history":
            forged.uns["openbio_singlecell"]["analysis_history"] = science.np.asarray([], dtype=float)

        path = tmp_path / f"forged_{issue}.h5ad"
        forged.write_h5ad(path)
        valid, reason = validate_existing(path, science.ad, science.sparse)
        assert not valid
        assert expected_reason in reason


def test_clustering_template_completes_demo_analysis(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_demo(input_dir, science)
    run = template_runner("Cell Clustering and Marker Discovery")

    loaded = run("OpenBioSingleCellLoadH5AD")
    normalized = run("OpenBioSingleCellNormalizeToLayer", loaded)
    variable = run("OpenBioSingleCellHighlyVariableGenes", normalized)
    highly_variable = int(variable.var["highly_variable"].sum())
    assert 0 < highly_variable < variable.n_vars

    pca = run("OpenBioSingleCellPCA", variable)
    neighbors = run("OpenBioSingleCellNeighbors", pca)
    umap = run("OpenBioSingleCellUMAP", neighbors)
    clustered = run("OpenBioSingleCellLeiden", umap)
    markers = run("OpenBioSingleCellMarkerGenes", clustered)
    filtered_markers = run("OpenBioSingleCellFilterMarkerGenes", markers)
    plot = run("OpenBioSingleCellUMAPPlot", clustered)
    summary = run("OpenBioSingleCellAnnDataSummary", clustered)

    assert markers.kind == "table" and not markers.table.empty
    assert science.np.isfinite(markers.table["logFC"].to_numpy(dtype=float)).all()
    assert filtered_markers.kind == "table" and not filtered_markers.table.empty
    assert plot.kind == "plot" and plot.png
    assert summary.kind == "summary" and summary.summary["shape"] == [clustered.n_obs, clustered.n_vars]


def test_quality_control_template_filters_demo_and_produces_consistent_plots(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_demo(input_dir, science)
    run = template_runner("Quality Control and Clean Counts")

    loaded = run("OpenBioSingleCellLoadH5AD")
    qc = run("OpenBioSingleCellCalculateQC", loaded)
    before_plot = run("OpenBioSingleCellQCPlots", qc)
    filtered_cells = run("OpenBioSingleCellFilterCells", qc)
    retained_plot = run("OpenBioSingleCellQCPlots", filtered_cells)
    filtered = run("OpenBioSingleCellFilterGenes", filtered_cells)
    summary = run("OpenBioSingleCellAnnDataSummary", filtered)

    assert (filtered.n_obs, filtered.n_vars) == (577, 500)
    assert before_plot.kind == "plot" and before_plot.png
    assert retained_plot.kind == "plot" and retained_plot.png
    assert summary.summary["shape"] == [577, 500]


def test_composition_template_produces_sample_level_tables(comfy_directories, science):
    input_dir, _, _ = comfy_directories
    write_demo(input_dir, science)
    run = template_runner("Sample Composition Comparison")

    sample_column, condition_column, batch_column, annotation_column, reference, comparison = run(
        "OpenBioSingleCellCoreStudyParameters"
    )
    loaded = run("OpenBioSingleCellLoadH5AD")
    composition = run(
        "OpenBioSingleCellSampleCompositionSummary",
        loaded,
        sample_column,
        condition_column,
        annotation_column,
    )
    contrast = run(
        "OpenBioSingleCellDifferentialCompositionTest",
        loaded,
        sample_column,
        condition_column,
        annotation_column,
        reference,
        comparison,
    )

    assert set(loaded.obs["condition"].astype(str)) == {"control", "treated"}
    assert batch_column == "batch"
    assert condition_column == "condition"
    assert reference == "control"
    assert comparison == "treated"
    assert composition.kind == "table" and len(composition.table) == 12
    assert contrast.kind == "table" and len(contrast.table) == 6
    assert set(contrast.table["scope"]) == {"global", "pairwise"}
