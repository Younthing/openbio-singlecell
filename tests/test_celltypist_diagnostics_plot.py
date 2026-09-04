from __future__ import annotations

import copy
import hashlib
import json
import uuid

import pytest

from openbio_singlecell.annotation_core import (
    _celltypist_axis_fingerprint,
    _celltypist_probability_fingerprint,
    _celltypist_selected_fingerprint,
)
from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata
from openbio_singlecell.nodes_annotation import OpenBioSingleCellCellTypistDiagnosticsPlot
from openbio_singlecell.operations_annotation import (
    celltypist_diagnostics_plot,
    celltypist_diagnostics_plot_owned,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _celltypist_output(science):
    obs_names = [f"cell_{index}" for index in range(8)]
    classes = ["A", "B", "C"]
    probabilities = science.np.asarray(
        [
            [0.90, 0.10, 0.05],
            [0.80, 0.20, 0.10],
            [0.70, 0.30, 0.10],
            [0.60, 0.40, 0.10],
            [0.20, 0.75, 0.10],
            [0.10, 0.85, 0.05],
            [0.15, 0.70, 0.20],
            [0.10, 0.20, 0.80],
        ],
        dtype=float,
    )
    labels = ["A", "A", "A", "A", "B", "B", "B", "C"]
    confidence = [float(probabilities[index, classes.index(label)]) for index, label in enumerate(labels)]
    obs = science.pd.DataFrame(
        {
            "leiden": science.pd.Categorical(
                ["0"] * 4 + ["1"] * 4,
                categories=["0", "1"],
                ordered=True,
            ),
            "celltypist_cell_type": science.pd.Categorical(labels, categories=classes),
            "celltypist_confidence": confidence,
            "celltypist_individual_label": science.pd.Categorical(labels, categories=classes),
            "celltypist_individual_probability": confidence,
        },
        index=obs_names,
    )
    adata = science.ad.AnnData(
        science.np.ones((8, 3)),
        obs=obs,
        var=science.pd.DataFrame(index=["G0", "G1", "G2"]),
    )
    adata.obsm["celltypist_probabilities"] = probabilities.copy()
    provenance = {
        "schema_version": 1,
        "operation": "celltypist_annotation",
        "annotation_status": "provisional",
        "expression": {
            "source": "layer",
            "layer_name": "log1p_norm",
            "declared_state": "verified_cp10k_log1p",
            "transform": "none",
            "target_sum": 10_000.0,
            "query_feature_count": 3,
        },
        "model": {
            "requested_identifier": "Immune_All_Low.pkl",
            "resolved_path": "C:/models/Immune_All_Low.pkl",
            "sha256": "a" * 64,
            "description": {},
            "class_count": 3,
            "feature_count": 3,
        },
        "feature_overlap": {"matched": 3},
        "classes": classes,
        "fixed_policy": {
            "mode": "best match",
            "p_thres": 0.5,
            "transpose_input": False,
            "gene_file": None,
            "cell_file": None,
            "use_GPU": False,
        },
        "majority_voting": False,
        "over_clustering_key": None,
        "min_prop": 0.0,
        "overwrite_existing": False,
        "columns": {
            "selected_label": "celltypist_cell_type",
            "selected_label_probability": "celltypist_confidence",
            "individual_label": "celltypist_individual_label",
            "individual_row_max_probability": "celltypist_individual_probability",
            "majority_label": None,
            "majority_support": None,
        },
        "matrices": {
            "probability_key": "celltypist_probabilities",
            "probability_shape": [8, 3],
            "decision_key": None,
            "decision_shape": None,
        },
        "artifact_replacement": {
            "prior_artifact_verified": False,
            "overwritten_obs_keys": [],
            "overwritten_obsm_keys": [],
            "removed_stale_obs_keys": [],
            "removed_stale_obsm_keys": [],
            "removed_stale_annotation_keys": [],
        },
        "software_versions": {},
        "warnings": [],
        "observation_axis_fingerprint_sha256": _celltypist_axis_fingerprint(obs_names),
        "probability_content_fingerprint_sha256": _celltypist_probability_fingerprint(
            probabilities,
            observation_ids=obs_names,
            classes=classes,
        ),
        "selected_annotation_fingerprint_sha256": _celltypist_selected_fingerprint(
            labels,
            confidence,
            observation_ids=obs_names,
        ),
    }
    adata.uns["celltypist"] = provenance
    adata.uns["openbio_singlecell"] = {
        "schema_version": 1,
        "annotations": {"celltypist_cell_type": copy.deepcopy(provenance)},
        "analysis_history": {
            "000000": {
                "operation": "celltypist_annotation",
                "parameters": {
                    "metadata_key": "celltypist",
                    "label_column": "celltypist_cell_type",
                    "probability_key": "celltypist_probabilities",
                    "model_sha256": "a" * 64,
                },
            }
        },
    }
    return adata


def test_celltypist_diagnostics_plot_schema_has_closed_same_contract_views():
    schema = OpenBioSingleCellCellTypistDiagnosticsPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellCellTypistDiagnosticsPlot"
    assert schema.display_name == "CellTypist Diagnostics Plot"
    assert schema.category == "openbio/single-cell/annotation"
    assert [item.id for item in schema.inputs] == ["adata", "metadata_key", "view"]
    assert [option.key for option in schema.inputs[2].options] == [
        "confidence_distributions",
        "probability_heatmap",
    ]
    assert [[item.id for item in option.inputs] for option in schema.inputs[2].options] == [
        ["max_labels"],
        ["groupby", "max_classes"],
    ]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_celltypist_confidence_plot_is_provisional_and_code_equivalent(science):
    adata = _celltypist_output(science)
    before = adata.copy()

    plotted, report, code = celltypist_diagnostics_plot_owned(
        adata,
        metadata_key="celltypist",
        view={"view": "confidence_distributions", "max_labels": 2},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["view"] == "confidence_distributions"
    assert details["annotation_status"] == "provisional"
    assert details["plotted_labels"] == ["A", "B"]
    assert details["label_counts"] == {"A": 4, "B": 3, "C": 1}
    assert details["plotted_cells"] == 8
    text = " ".join([report.summary["results"], *report.summary["limitations"]])
    assert "Provisional annotation" in text
    assert "Condition contrast" in text
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<celltypist-diagnostics-code>", "exec"), namespace)
    assert namespace["plot_celltypist_diagnostics"](adata) == plotted.png
    science.pd.testing.assert_frame_equal(adata.obs, before.obs)
    science.np.testing.assert_array_equal(
        adata.obsm["celltypist_probabilities"], before.obsm["celltypist_probabilities"]
    )
    assert adata.uns == before.uns


def test_celltypist_probability_heatmap_uses_stored_probabilities_and_explicit_grouping(science):
    adata = _celltypist_output(science)

    plotted, report, code = celltypist_diagnostics_plot_owned(
        adata,
        metadata_key="celltypist",
        view={"view": "probability_heatmap", "groupby": "leiden", "max_classes": 3},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["view"] == "probability_heatmap"
    assert details["groups"] == ["0", "1"]
    assert details["plotted_classes"] == ["A", "B", "C"]
    science.np.testing.assert_allclose(
        details["group_mean_probabilities"],
        [[0.75, 0.25, 0.0875], [0.1375, 0.625, 0.2875]],
    )
    namespace: dict[str, object] = {}
    exec(compile(code, "<celltypist-heatmap-code>", "exec"), namespace)
    assert namespace["plot_celltypist_diagnostics"](adata) == plotted.png


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda adata: adata.obsm["celltypist_probabilities"].__setitem__((0, 0), 0.4),
            "probability.*fingerprint",
        ),
        (lambda adata: setattr(adata, "obs_names", list(reversed(adata.obs_names))), "axis fingerprint"),
        (lambda adata: adata.uns["celltypist"].__setitem__("operation", "map_cluster_annotations"), "producer"),
    ],
)
def test_celltypist_diagnostics_rejects_tampered_producer_evidence(science, mutation, message):
    adata = _celltypist_output(science)
    mutation(adata)

    with pytest.raises((TypeError, ValueError), match=message):
        celltypist_diagnostics_plot_owned(adata)


def test_celltypist_diagnostics_rejects_inactive_view_parameters(science):
    adata = _celltypist_output(science)

    with pytest.raises(ValueError, match="inactive or unknown"):
        celltypist_diagnostics_plot_owned(
            adata,
            view={
                "view": "confidence_distributions",
                "max_labels": 2,
                "groupby": "leiden",
            },
        )


def test_celltypist_diagnostics_worker_preserves_the_input_h5ad(tmp_path, science):
    adata = _celltypist_output(science)
    input_root = tmp_path / "celltypist"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).digest()
    staging = tmp_path / "celltypist-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = celltypist_diagnostics_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        {
            "metadata_key": "celltypist",
            "view": {"view": "confidence_distributions", "max_labels": 2},
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"]["view"]["view"] == "confidence_distributions"
    assert hashlib.sha256(input_path.read_bytes()).digest() == before


def test_celltypist_diagnostics_restores_matplotlib_global_state(science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    matplotlib.rcParams["lines.linewidth"] = 3.25
    before_rc = matplotlib.rcParams["lines.linewidth"]
    before_figures = pyplot.get_fignums()

    celltypist_diagnostics_plot_owned(_celltypist_output(science))

    assert matplotlib.rcParams["lines.linewidth"] == before_rc
    assert pyplot.get_fignums() == before_figures
