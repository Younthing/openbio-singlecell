from __future__ import annotations

import hashlib
import json
import uuid

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from openbio_singlecell.abundance_artifact_codecs import (
    COMPOSITION_MODEL_CODEC,
    COMPOSITION_MODEL_KIND,
    MILO_RESULT_CODEC,
    MILO_RESULT_KIND,
    write_composition_model_result,
    write_milo_result,
)
from openbio_singlecell.artifact_codecs import read_plot
from openbio_singlecell.composition_result import (
    SCCODA_RESULT_COLUMNS,
    TASCCODA_RESULT_COLUMNS,
    build_composition_model_result,
)
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.milo_analysis import MILO_TABLE_COLUMNS
from openbio_singlecell.milo_result import build_milo_result
from openbio_singlecell.nodes_abundance import (
    OpenBioSingleCellMiloDifferentialAbundancePlot,
    OpenBioSingleCellSccodaDifferentialCompositionPlot,
    OpenBioSingleCellTasccodaDifferentialCompositionPlot,
)
from openbio_singlecell.operations_abundance import (
    milo_differential_abundance_plot,
    milo_differential_abundance_plot_owned,
    sccoda_differential_composition_plot,
    sccoda_differential_composition_plot_owned,
    tasccoda_differential_composition_plot,
    tasccoda_differential_composition_plot_owned,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse, registered_operation_ids

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _milo_result():
    table = pd.DataFrame(
        [
            ["nhood_1", "c1", 2, 0.5, "control", "treated", -1.0, 4.0, 8.0, 0.01, 0.02, 0.03, "T", 1.0, "T", False],
            ["nhood_2", "c3", 2, 0.7, "control", "treated", 0.8, 5.0, 7.0, 0.02, 0.03, 0.04, "B", 1.0, "B", False],
        ],
        columns=MILO_TABLE_COLUMNS,
    )
    membership = sparse.csr_matrix([[1, 0], [1, 1], [0, 1]], dtype=np.int8)
    graph = membership.T @ membership
    graph.setdiag(0)
    graph.eliminate_zeros()
    return build_milo_result(
        table=table,
        membership=membership,
        graph=graph,
        representative_coordinates=np.array([[0.0, 1.0], [2.0, 3.0]]),
        observation_names=["c1", "c2", "c3"],
        neighborhood_names=["nhood_1", "nhood_2"],
        representation_key="X_pca",
        representation_sha256="a" * 64,
        condition_key="condition",
        annotation_key="cell_type",
        annotation_status="curated",
        n_neighbors=3,
        neighborhood_proportion=0.2,
        mixed_annotation_threshold=0.6,
        spatial_fdr_threshold=0.1,
        min_abs_log2_fold_change=0.0,
        random_seed=7,
    )


def test_abundance_plot_nodes_use_method_specific_typed_contracts():
    milo = OpenBioSingleCellMiloDifferentialAbundancePlot.define_schema()
    sccoda = OpenBioSingleCellSccodaDifferentialCompositionPlot.define_schema()
    tasccoda = OpenBioSingleCellTasccodaDifferentialCompositionPlot.define_schema()

    assert milo.node_id == "OpenBioSingleCellMiloDifferentialAbundancePlot"
    assert milo.display_name == "Milo Differential Abundance Plot"
    assert milo.category == "openbio/single-cell/differential-abundance"
    assert [item.id for item in milo.inputs] == ["result", "view"]
    assert [option.key for option in milo.inputs[1].options] == [
        "differential_evidence",
        "neighborhood_graph",
    ]

    assert sccoda.node_id == "OpenBioSingleCellSccodaDifferentialCompositionPlot"
    assert sccoda.display_name == "scCODA Differential Composition Plot"
    assert [item.id for item in sccoda.inputs] == ["result", "view"]
    assert [option.key for option in sccoda.inputs[1].options] == [
        "effect_forest",
        "posterior_distribution",
        "sampler_diagnostics",
    ]
    assert [item.id for item in sccoda.inputs[1].options[1].inputs] == ["cell_type"]

    assert tasccoda.node_id == "OpenBioSingleCellTasccodaDifferentialCompositionPlot"
    assert tasccoda.display_name == "tascCODA Differential Composition Plot"
    assert [item.id for item in tasccoda.inputs] == ["result", "view"]
    assert [option.key for option in tasccoda.inputs[1].options] == [
        "hierarchy_effect_forest",
        "derived_leaf_effects",
        "posterior_distribution",
        "sampler_diagnostics",
    ]
    assert [item.id for item in tasccoda.inputs[1].options[2].inputs] == ["effect_scope", "effect_name"]
    assert [item.display_name for item in milo.outputs] == ["plot", "summary", "code"]
    assert [item.display_name for item in sccoda.outputs] == ["plot", "summary", "code"]
    assert [item.display_name for item in tasccoda.outputs] == ["plot", "summary", "code"]


def test_milo_plot_renders_tamper_evident_neighborhood_evidence_without_reanalysis():
    result = _milo_result()
    before = result.portable()

    plotted, report, code = milo_differential_abundance_plot_owned(
        result,
        view={"view": "differential_evidence"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellMiloDifferentialAbundancePlot"
    assert summary["parameters"] == {"view": {"view": "differential_evidence"}}
    assert summary["key_results"]["neighborhoods"] == 2
    assert summary["key_results"]["spatial_fdr_selected"] == 2
    assert summary["key_results"]["annotation_key"] == "cell_type"
    assert summary["key_results"]["annotation_status"] == "curated"
    assert "neighborhood" in summary["methods"].lower()
    assert "cell-level" not in summary["results"].lower()
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<milo-plot-code>", "exec"), namespace)
    assert namespace["plot_milo_differential_abundance"](result) == plotted.png
    after = result.portable()
    pd.testing.assert_frame_equal(after["table"], before["table"])
    np.testing.assert_array_equal(after["membership"].toarray(), before["membership"].toarray())
    np.testing.assert_array_equal(after["graph"].toarray(), before["graph"].toarray())
    np.testing.assert_array_equal(after["coordinates"], before["coordinates"])
    assert after["metadata"] == before["metadata"]


def _sccoda_result(*, effect_center=0.4):
    table = pd.DataFrame(
        [
            ["treated_vs_control", "condition", "control", "treated", "B", "T", effect_center, 0.2, 0.6, 0.1, 0.95, True, 0.05, 0.9, 0.03, 10.0, 12.0, 0.2, False],
            ["treated_vs_control", "condition", "control", "treated", "T", "T", 0.0, 0.0, 0.0, 0.0, 0.0, False, 0.05, 0.9, 0.03, 8.0, 6.0, -0.4, True],
        ],
        columns=SCCODA_RESULT_COLUMNS,
    )
    return build_composition_model_result(
        table=table,
        method="sccoda",
        posterior={
            "intercept": np.array([[0.1, -0.1], [0.2, -0.2], [0.0, 0.0], [0.15, -0.15]]),
            "condition_effect": np.array([[0.2, 0.0], [0.4, 0.0], [0.5, 0.0], [0.3, 0.0]]),
        },
        sample_stats={
            "potential_energy": np.array([10.0, 10.5, 9.8, 10.1]),
            "num_steps": np.array([7.0, 8.0, 7.0, 9.0]),
            "step_size": np.array([0.1, 0.1, 0.1, 0.1]),
        },
        cell_types=["B", "T"],
        hierarchy=None,
        model_metadata={
            "condition_key": "condition",
            "reference_condition": "control",
            "comparison_condition": "treated",
            "reference_cell_type": "T",
            "annotation_key": "cell_type",
            "annotation_status": "curated",
            "model": {
                "family": "Dirichlet-multinomial",
                "method": "sccoda",
                "chain_count": 1,
                "posterior_draws": 4,
            },
            "diagnostics": {"rhat_available": False, "divergences_available": False},
            "selection": {"method": "posterior_expected_fdr"},
            "software_versions": {"pertpy": "1.3.0", "arviz": "0.22.0"},
        },
    )


def test_sccoda_plot_preserves_relative_compositional_effect_semantics_and_posterior_evidence():
    result = _sccoda_result()
    before = result.portable()

    plotted, report, code = sccoda_differential_composition_plot_owned(
        result,
        view={"view": "effect_forest"},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellSccodaDifferentialCompositionPlot"
    assert summary["parameters"] == {"view": {"view": "effect_forest"}}
    assert summary["key_results"]["method"] == "scCODA"
    assert summary["key_results"]["credible_effects"] == 1
    assert summary["key_results"]["reference_cell_type"] == "T"
    assert summary["key_results"]["posterior_draws"] == 4
    assert "relative" in summary["results"].lower()
    assert "absolute abundance" not in summary["results"].lower()
    json.dumps(summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<sccoda-plot-code>", "exec"), namespace)
    assert namespace["plot_sccoda_differential_composition"](result) == plotted.png
    after = result.portable()
    pd.testing.assert_frame_equal(after["table"], before["table"])
    for name in before["posterior"]:
        np.testing.assert_array_equal(after["posterior"][name], before["posterior"][name])
    assert after["metadata"] == before["metadata"]


def test_composition_effect_forest_allows_posterior_center_outside_hdi():
    result = _sccoda_result(effect_center=0.8)
    before = result.portable()["table"]

    plotted, _report, code = sccoda_differential_composition_plot_owned(result, view={"view": "effect_forest"})

    assert plotted.png.startswith(PNG_SIGNATURE)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_sccoda_differential_composition"](result) == plotted.png
    pd.testing.assert_frame_equal(result.portable()["table"], before)


def _tasccoda_result():
    rows = [
        ["hierarchy_node", "treated_vs_control", "condition", "control", "treated", "T", "B", "cell_type", 1, '["B"]', 0.4, 0.4, 0.2, 0.6, 0.1, 0.4, True, "direct_node_abs_median_strictly_greater_than_delta", pd.NA, pd.NA, pd.NA, False],
        ["hierarchy_node", "treated_vs_control", "condition", "control", "treated", "T", "T", "cell_type", 1, '["T"]', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, "direct_node_abs_median_strictly_greater_than_delta", pd.NA, pd.NA, pd.NA, True],
        ["derived_leaf", "treated_vs_control", "condition", "control", "treated", "T", "B", "cell_type", 1, '["B"]', 0.4, 0.4, 0.2, 0.6, 0.1, pd.NA, pd.NA, "derived_sum_of_selected_hierarchy_node_effects", 10.0, 12.0, 0.2, False],
        ["derived_leaf", "treated_vs_control", "condition", "control", "treated", "T", "T", "cell_type", 1, '["T"]', 0.0, 0.0, 0.0, 0.0, 0.0, pd.NA, pd.NA, "derived_sum_of_selected_hierarchy_node_effects", 8.0, 6.0, -0.4, True],
    ]
    return build_composition_model_result(
        table=pd.DataFrame(rows, columns=TASCCODA_RESULT_COLUMNS),
        method="tasccoda",
        posterior={
            "intercept": np.array([[0.1, -0.1], [0.2, -0.2], [0.0, 0.0], [0.15, -0.15]]),
            "hierarchy_node_effect": np.array([[0.2, 0.0], [0.4, 0.0], [0.5, 0.0], [0.3, 0.0]]),
            "derived_leaf_effect": np.array([[0.2, 0.0], [0.4, 0.0], [0.5, 0.0], [0.3, 0.0]]),
            "theta": np.array([0.4, 0.5, 0.45, 0.55]),
        },
        sample_stats={
            "potential_energy": np.array([10.0, 10.5, 9.8, 10.1]),
            "num_steps": np.array([7.0, 8.0, 7.0, 9.0]),
            "step_size": np.array([0.1, 0.1, 0.1, 0.1]),
        },
        cell_types=["B", "T"],
        hierarchy={
            "ancestor_keys": ["lineage"],
            "leaf_key": "cell_type",
            "levels_passed_to_pertpy": ["lineage", "cell_type"],
            "root": "root",
            "edges": [["root", "B"], ["root", "T"]],
            "node_names": ["B", "T"],
            "ancestor_matrix_shape": [2, 2],
            "fingerprint": "a" * 64,
            "reference_path": ["root", "T"],
            "declared_reference_path": ["root", "T"],
            "reference_nodes": ["T"],
            "derived_leaf_effects_are_propagated": True,
        },
        model_metadata={
            "condition_key": "condition",
            "reference_condition": "control",
            "comparison_condition": "treated",
            "reference_cell_type": "T",
            "annotation_key": "cell_type",
            "annotation_status": "curated",
            "model": {
                "family": "Dirichlet-multinomial",
                "method": "tasccoda",
                "chain_count": 1,
                "posterior_draws": 4,
            },
            "diagnostics": {"rhat_available": False, "divergences_available": False},
            "selection": {"method": "tree_adaptive_spike_and_slab_lasso"},
            "software_versions": {"pertpy": "1.3.0", "arviz": "0.22.0"},
        },
    )


def test_tasccoda_plot_counts_direct_hierarchy_evidence_without_promoting_derived_leaves():
    result = _tasccoda_result()

    plotted, report, code = tasccoda_differential_composition_plot_owned(
        result,
        view={"view": "hierarchy_effect_forest"},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    summary = report.summary
    assert summary["node_id"] == "OpenBioSingleCellTasccodaDifferentialCompositionPlot"
    assert summary["key_results"]["method"] == "tascCODA"
    assert summary["key_results"]["credible_effects"] == 1
    assert summary["key_results"]["plotted_effects"] == 2
    assert "propagated leaf effects" in summary["methods"]
    assert any("independent discoveries" in limitation for limitation in summary["limitations"])

    namespace: dict[str, object] = {}
    exec(compile(code, "<tasccoda-plot-code>", "exec"), namespace)
    assert namespace["plot_tasccoda_differential_composition"](result) == plotted.png


def test_milo_neighborhood_graph_uses_retained_membership_graph_and_coordinates():
    result = _milo_result()

    plotted, report, code = milo_differential_abundance_plot_owned(
        result,
        view={"view": "neighborhood_graph"},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["plotted_edges"] == 1
    assert details["membership_nonzero"] == 4
    assert details["representation_key"] == "X_pca"
    namespace: dict[str, object] = {}
    exec(compile(code, "<milo-neighborhood-graph-code>", "exec"), namespace)
    assert namespace["plot_milo_differential_abundance"](result) == plotted.png


@pytest.mark.parametrize(
    "view",
    [
        {"view": "posterior_distribution", "cell_type": ""},
        {"view": "sampler_diagnostics"},
    ],
)
def test_sccoda_posterior_views_render_retained_draws_without_refitting(view):
    result = _sccoda_result()

    plotted, report, code = sccoda_differential_composition_plot_owned(result, view=view)

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["key_results"]["view"] == view["view"]
    assert report.summary["key_results"]["posterior_draws"] == 4
    namespace: dict[str, object] = {}
    exec(compile(code, "<sccoda-posterior-view-code>", "exec"), namespace)
    assert namespace["plot_sccoda_differential_composition"](result) == plotted.png


@pytest.mark.parametrize(
    "view",
    [
        {"view": "derived_leaf_effects"},
        {"view": "posterior_distribution", "effect_scope": "hierarchy_node", "effect_name": ""},
        {"view": "posterior_distribution", "effect_scope": "derived_leaf", "effect_name": "B"},
        {"view": "sampler_diagnostics"},
    ],
)
def test_tasccoda_views_keep_direct_and_derived_effect_scopes_explicit(view):
    result = _tasccoda_result()

    plotted, report, code = tasccoda_differential_composition_plot_owned(result, view=view)

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["view"] == view["view"]
    if view["view"] == "derived_leaf_effects":
        assert details["derived_leaf_credibility"] == (
            "not_applicable_propagated_effects_are_not_direct_selections"
        )
    if view["view"] == "posterior_distribution":
        assert details["selected_effect_scope"] == view["effect_scope"]
    namespace: dict[str, object] = {}
    exec(compile(code, "<tasccoda-view-code>", "exec"), namespace)
    assert namespace["plot_tasccoda_differential_composition"](result) == plotted.png


def test_abundance_plots_reject_cross_method_and_inactive_view_fields():
    with pytest.raises(ValueError, match="exact scCODA"):
        sccoda_differential_composition_plot_owned(_tasccoda_result())
    with pytest.raises(ValueError, match="closed supported"):
        milo_differential_abundance_plot_owned(
            _milo_result(),
            view={"view": "differential_evidence", "effect_name": "B"},
        )
    with pytest.raises(ValueError, match="closed supported"):
        tasccoda_differential_composition_plot_owned(
            _tasccoda_result(),
            view={"view": "derived_leaf_effects", "effect_scope": "derived_leaf"},
        )


def _tree_hash(root):
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).digest()
        for path in root.rglob("*")
        if path.is_file()
    }


def test_milo_plot_worker_round_trips_png_without_rewriting_typed_evidence(tmp_path):
    input_root = tmp_path / "milo"
    input_root.mkdir()
    write_milo_result(input_root, _milo_result())
    before = _tree_hash(input_root)
    staging = tmp_path / "milo-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = milo_differential_abundance_plot(
        context,
        {
            "result": {
                "type": "artifact",
                "kind": MILO_RESULT_KIND,
                "codec": MILO_RESULT_CODEC,
                "path": str(input_root.resolve()),
            }
        },
        {"view": {"view": "neighborhood_graph"}},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"view": {"view": "neighborhood_graph"}}
    assert _tree_hash(input_root) == before


@pytest.mark.parametrize(
    ("operation", "build", "view"),
    [
        (sccoda_differential_composition_plot, _sccoda_result, {"view": "effect_forest"}),
        (
            tasccoda_differential_composition_plot,
            _tasccoda_result,
            {"view": "posterior_distribution", "effect_scope": "derived_leaf", "effect_name": "B"},
        ),
    ],
)
def test_composition_plot_workers_round_trip_png_without_rewriting_inference_data(
    tmp_path,
    operation,
    build,
    view,
):
    input_root = tmp_path / operation.__name__
    input_root.mkdir()
    write_composition_model_result(input_root, build())
    before = _tree_hash(input_root)
    staging = tmp_path / f"{operation.__name__}.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = operation(
        context,
        {
            "result": {
                "type": "artifact",
                "kind": COMPOSITION_MODEL_KIND,
                "codec": COMPOSITION_MODEL_CODEC,
                "path": str(input_root.resolve()),
            }
        },
        {"view": view},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"view": view}
    assert _tree_hash(input_root) == before


def test_abundance_plot_worker_operation_ids_are_registered():
    assert {
        "openbio.node.milodifferentialabundanceplot",
        "openbio.node.sccodadifferentialcompositionplot",
        "openbio.node.tasccodadifferentialcompositionplot",
    }.issubset(registered_operation_ids())


def test_abundance_plots_restore_matplotlib_and_numpy_global_state():
    import matplotlib
    import matplotlib.pyplot as pyplot

    original_linewidth = matplotlib.rcParams["lines.linewidth"]
    try:
        matplotlib.rcParams["lines.linewidth"] = 3.25
        expected_linewidth = matplotlib.rcParams["lines.linewidth"]
        expected_figures = pyplot.get_fignums()
        expected_numpy = np.random.get_state()

        milo_differential_abundance_plot_owned(_milo_result())
        sccoda_differential_composition_plot_owned(_sccoda_result())
        tasccoda_differential_composition_plot_owned(_tasccoda_result())

        assert matplotlib.rcParams["lines.linewidth"] == expected_linewidth
        assert pyplot.get_fignums() == expected_figures
        observed_numpy = np.random.get_state()
        assert observed_numpy[0] == expected_numpy[0]
        np.testing.assert_array_equal(observed_numpy[1], expected_numpy[1])
        assert observed_numpy[2:] == expected_numpy[2:]
    finally:
        matplotlib.rcParams["lines.linewidth"] = original_linewidth
