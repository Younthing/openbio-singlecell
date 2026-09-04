from __future__ import annotations

import hashlib
import json
import time
import uuid

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.contracts import PlotResult
from openbio_singlecell.nodes_regulatory import (
    OpenBioSingleCellSCENICActivityPlot,
    OpenBioSingleCellSCENICBinarizationPlot,
    OpenBioSingleCellSCENICRegulonMembershipPlot,
    OpenBioSingleCellSCENICRegulonSpecificityPlot,
    OpenBioSingleCellTFActivityPlot,
    OpenBioSingleCellTFActivityRankingPlot,
)
from openbio_singlecell.operations_regulatory import (
    scenic_activity_plot_owned,
    scenic_binarization_plot_owned,
    scenic_membership_plot_owned,
    scenic_rss_plot_owned,
    tf_activity_plot_owned,
    tf_activity_ranking_plot_owned,
)
from openbio_singlecell.regulatory_artifact_codecs import (
    SCENIC_CODEC,
    TF_ACTIVITY_CODEC,
    write_scenic,
    write_tf_activity,
)
from openbio_singlecell.regulatory_plotting import table_evidence_fingerprint
from openbio_singlecell.scenic_artifact import (
    SCENIC_MEMBERSHIP_COLUMNS,
    SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
    SCENIC_RESULT_ARTIFACT_TYPE,
    SCENIC_RESULT_PRODUCER_NODE_ID,
    SCENIC_RESULT_PRODUCER_SCHEMA,
    SCENICResultArtifact,
)
from openbio_singlecell.scenic_binarization import SCENIC_THRESHOLD_COLUMNS
from openbio_singlecell.scenic_rss import SCENIC_RSS_COLUMNS
from openbio_singlecell.tf_activity_artifact import build_tf_activity_artifact
from openbio_singlecell.tf_activity_ranking import TF_RANKING_COLUMNS
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _adata(science):
    names = [f"cell_{index}" for index in range(6)]
    adata = science.ad.AnnData(
        science.np.ones((6, 3)),
        obs=science.pd.DataFrame(
            {"cell_type": science.pd.Categorical(["T", "T", "T", "B", "B", "B"], categories=["T", "B"])},
            index=names,
        ),
        var=science.pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    adata.obsm["X_umap"] = science.np.asarray(
        [[0.0, 0.0], [0.2, 0.1], [0.4, 0.2], [1.0, 1.0], [1.2, 1.1], [1.4, 1.2]]
    )
    return adata


def _tf_artifact(science, adata):
    scores = science.pd.DataFrame(
        {"TF1": [2.0, 1.5, 1.0, -1.0, -1.5, -2.0], "TF2": [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5]},
        index=adata.obs_names,
    )
    adjusted = science.pd.DataFrame(
        {"TF1": [0.01, 0.02, 0.03, 0.2, 0.3, 0.4], "TF2": [0.4, 0.3, 0.2, 0.03, 0.02, 0.01]},
        index=adata.obs_names,
    )
    return build_tf_activity_artifact(
        scores=scores,
        adjusted_pvalues=adjusted,
        provenance={"resource_sha256": "a" * 64, "observation_axis_sha256": "b" * 64},
        numpy=science.np,
        pandas=science.pd,
    )


def _scenic_artifact(science, adata):
    activity = science.pd.DataFrame(
        {"TF1(+)": [0.8, 0.7, 0.6, 0.2, 0.1, 0.0], "TF2(-)": [0.1, 0.2, 0.3, 0.7, 0.8, 0.9]},
        index=adata.obs_names,
    )
    membership = science.pd.DataFrame(
        [
            ["TF1(+)", "TF1", "activating", "ctx", "G1", 0.9, 1, ["M1"], 1],
            ["TF1(+)", "TF1", "activating", "ctx", "G2", 0.4, 1, ["M1"], 2],
            ["TF2(-)", "TF2", "repressing", "ctx", "G3", 0.7, 1, ["M2"], 1],
        ],
        columns=SCENIC_MEMBERSHIP_COLUMNS,
    )
    resources = [
        {
            "role": "tf_list", "file": "tfs.txt", "sha256": "1" * 64, "release": "fixture",
            "organism": "human", "genome_build": "GRCh38", "gene_namespace": "HGNC",
            "license": "fixture", "citation": "fixture",
        }
    ]
    accounting = {
        label: {"file": f"{label.replace(':', '_')}.dat", "sha256": "1" * 64, "bytes": 1}
        for label in ["manifest", "expression", "adjacency", "regulons", "aucell", "resource:0:tf_list"]
    }
    accounting["manifest"]["sha256"] = "2" * 64
    accounting["resource:0:tf_list"] = {"file": "tfs.txt", "sha256": "1" * 64, "bytes": 1}
    provenance = {
        "artifact_type": SCENIC_RESULT_ARTIFACT_TYPE,
        "artifact_schema_version": SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
        "producer_node_id": SCENIC_RESULT_PRODUCER_NODE_ID,
        "producer_schema": SCENIC_RESULT_PRODUCER_SCHEMA,
        "manifest_schema": "openbio-singlecell/pyscenic-external-run/v1",
        "manifest_sha256": "2" * 64,
        "run_id": "fixture-run",
        "organism": "human",
        "genome_build": "GRCh38",
        "gene_namespace": "HGNC",
        "expression_state": "normalized",
        "input_dimensions": {"cells": 6, "genes": 3},
        "expression_fingerprints": {
            "observation_ids_sha256": "3" * 64,
            "gene_ids_sha256": "4" * 64,
            "matrix_sha256": "5" * 64,
        },
        "container": None,
        "declared_software_versions": {
            "pyscenic": "0.12.1", "ctxcore": "0.2", "arboreto": "0.1", "python": "3.12"
        },
        "commands": {
            "grn": ["pyscenic", "grn", "expression.csv"],
            "ctx": ["pyscenic", "ctx", "adjacency.csv"],
            "aucell": ["pyscenic", "aucell", "expression.csv"],
        },
        "resources": resources,
        "bundle_accounting": accounting,
    }
    return SCENICResultArtifact(activity, membership, provenance)


def _table_result(science, contract: str):
    if contract == "tf_ranking":
        frame = science.pd.DataFrame(
            [
                ["T", "rest", "TF1", 4.0, 2.0, 0.001, 0.01, "higher_in_group", True, 1],
                ["B", "rest", "TF2", 3.0, 1.5, 0.002, 0.02, "higher_in_group", True, 1],
            ],
            columns=TF_RANKING_COLUMNS,
        )
        operation = "rank_tf_activities"
    elif contract == "scenic_rss":
        frame = science.pd.DataFrame(
            [["T", "TF1(+)", 0.9, 1, 3, 6], ["T", "TF2(-)", 0.2, 2, 3, 6], ["B", "TF2(-)", 0.85, 1, 3, 6], ["B", "TF1(+)", 0.25, 2, 3, 6]],
            columns=SCENIC_RSS_COLUMNS,
        )
        operation = "scenic_regulon_specificity"
    elif contract == "scenic_binarization":
        frame = science.pd.DataFrame(
            [["TF1(+)", 0.5, "pyscenic_0.12.1_hdt", 3, 6, 0.5], ["TF2(-)", 0.6, "manual_override", 2, 6, 1 / 3]],
            columns=SCENIC_THRESHOLD_COLUMNS,
        )
        operation = "scenic_activity_binarization"
    else:
        frame = science.pd.DataFrame(
            [
                ["TF1(+)", "TF1", "activating", "ctx", "G1", 0.9, 1, ["M1"], 1],
                ["TF1(+)", "TF1", "activating", "ctx", "G2", 0.4, 1, ["M1"], 2],
                ["TF2(-)", "TF2", "repressing", "ctx", "G3", 0.7, 1, ["M2"], 1],
            ],
            columns=SCENIC_MEMBERSHIP_COLUMNS,
        )
        operation = "scenic_final_regulon_membership"
    parameters = {
        "plot_evidence_schema": "openbio-singlecell/plot-evidence/v1",
        "plot_evidence_fingerprint_sha256": table_evidence_fingerprint(frame, contract=contract),
    }
    return make_table_result(
        table=frame,
        title=contract,
        operation=operation,
        parameters=parameters,
        description="fixture",
        warnings=[],
        input_cells=6,
        input_genes=3,
        started_at=time.perf_counter(),
    )


def _many_group_tf_result(science, groups: int, rows_per_group: int = 1):
    frame = science.pd.DataFrame(
        [
            [
                f"group_{index:03d}",
                "rest",
                f"TF{rank}",
                float(rows_per_group - rank + 1),
                0.5,
                0.01,
                0.02,
                "higher_in_group",
                True,
                rank,
            ]
            for index in range(groups)
            for rank in range(1, rows_per_group + 1)
        ],
        columns=TF_RANKING_COLUMNS,
    )
    parameters = {
        "report_p_adjusted": 0.05,
        "plot_evidence_schema": "openbio-singlecell/plot-evidence/v1",
        "plot_evidence_fingerprint_sha256": table_evidence_fingerprint(frame, contract="tf_ranking"),
    }
    return make_table_result(
        table=frame,
        title="tf_ranking",
        operation="rank_tf_activities",
        parameters=parameters,
        description="many-group fixture",
        warnings=[],
        input_cells=groups,
        input_genes=3,
        started_at=time.perf_counter(),
    )


def test_regulatory_plot_node_schemas_use_exact_artifact_types_and_closed_views():
    tf = OpenBioSingleCellTFActivityPlot.define_schema()
    scenic = OpenBioSingleCellSCENICActivityPlot.define_schema()
    assert [item.id for item in tf.inputs] == ["adata", "activities", "view"]
    assert [item.id for item in scenic.inputs] == ["adata", "scenic_result", "view"]
    assert [option.key for option in tf.inputs[2].options] == ["distributions", "grouped_heatmap", "embedding"]
    assert [option.key for option in scenic.inputs[2].options] == [
        "distributions",
        "grouped_heatmap",
        "embedding",
    ]
    assert tf.category == scenic.category == "openbio/single-cell/regulatory"
    assert [item.display_name for item in tf.outputs] == ["plot", "summary", "code"]

    cases = [
        (OpenBioSingleCellTFActivityRankingPlot, ["ranked_dot", "effect_significance"]),
        (OpenBioSingleCellSCENICRegulonSpecificityPlot, ["rss_heatmap", "top_regulons"]),
        (OpenBioSingleCellSCENICBinarizationPlot, ["thresholds", "active_fraction"]),
        (OpenBioSingleCellSCENICRegulonMembershipPlot, ["target_count", "top_targets"]),
    ]
    for node, views in cases:
        schema = node.define_schema()
        assert schema.category == "openbio/single-cell/regulatory"
        assert [item.id for item in schema.inputs] == ["table", "view"]
        assert [option.key for option in schema.inputs[1].options] == views
        assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


@pytest.mark.parametrize("kind", ["tf", "scenic"])
def test_regulatory_activity_plots_validate_typed_artifact_axes_and_report_group_means(science, kind):
    adata = _adata(science)
    artifact = _tf_artifact(science, adata) if kind == "tf" else _scenic_artifact(science, adata)
    snapshot = adata.copy()
    owned = tf_activity_plot_owned if kind == "tf" else scenic_activity_plot_owned
    kwargs = {"activities": artifact} if kind == "tf" else {"scenic_result": artifact}

    plotted, report, code = owned(
        adata,
        **kwargs,
        view={"view": "grouped_heatmap", "groupby": "cell_type", "max_activities": 10},
    )

    assert isinstance(plotted, PlotResult)
    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["groups"] == ["T", "B"]
    assert details["plotted_cells"] == 6
    assert details["plotted_activities"] == 2
    assert len(details["group_means"]) == 2
    json.dumps(report.summary, allow_nan=False)
    namespace: dict[str, object] = {}
    exec(compile(code, f"<{kind}-activity-plot>", "exec"), namespace)
    function = namespace[f"plot_{kind}_activity"]
    portable = artifact.portable() if kind == "tf" else artifact.to_portable()
    assert function(adata, artifact) == plotted.png
    assert function(adata, portable) == plotted.png
    science.np.testing.assert_array_equal(adata.X, snapshot.X)


@pytest.mark.parametrize(
    ("owned", "contract", "view", "function_name", "mark_key", "expected"),
    [
        (tf_activity_ranking_plot_owned, "tf_ranking", {"view": "ranked_dot", "max_per_group": 10}, "plot_tf_activity_ranking", "statistics", [4.0, 3.0]),
        (scenic_rss_plot_owned, "scenic_rss", {"view": "rss_heatmap", "max_regulons": 10}, "plot_scenic_regulon_specificity", "rss_values", [0.9, 0.2, 0.85, 0.25]),
        (scenic_binarization_plot_owned, "scenic_binarization", {"view": "thresholds", "max_regulons": 10}, "plot_scenic_binarization", "thresholds", [0.5, 0.6]),
        (scenic_membership_plot_owned, "scenic_membership", {"view": "target_count", "max_regulons": 10}, "plot_scenic_regulon_membership", "target_counts", [2, 1]),
    ],
)
def test_regulatory_table_plots_report_exact_marks_and_code_parity(
    science, owned, contract, view, function_name, mark_key, expected
):
    result = _table_result(science, contract)
    snapshot = result.table.copy(deep=True)

    plotted, report, code = owned(result, view=view)

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details[mark_key] == expected
    assert details["result_fingerprint_sha256"] == result.parameters["plot_evidence_fingerprint_sha256"]
    namespace: dict[str, object] = {}
    exec(compile(code, f"<{function_name}>", "exec"), namespace)
    assert namespace[function_name](result.table) == plotted.png
    science.pd.testing.assert_frame_equal(result.table, snapshot)


@pytest.mark.parametrize(
    ("groups", "rows_per_group", "message"),
    [(101, 1, "at most 100 groups"), (50, 3, "at most 100 plotted marks")],
)
def test_tf_activity_ranking_preflights_group_and_mark_bounds_before_figure_allocation(
    science,
    monkeypatch,
    groups,
    rows_per_group,
    message,
):
    result = _many_group_tf_result(science, groups, rows_per_group)

    def reject_figure(*_args, **_kwargs):
        raise AssertionError("Figure was allocated before TF-ranking mark preflight")

    monkeypatch.setattr("matplotlib.figure.Figure", reject_figure)
    with pytest.raises(ValueError, match=message):
        tf_activity_ranking_plot_owned(
            result,
            view={"view": "ranked_dot", "max_per_group": rows_per_group},
        )


def test_regulatory_plots_fail_closed_on_axis_and_table_tampering(science):
    adata = _adata(science)
    activities = _tf_artifact(science, adata)
    reordered = adata[[1, 0, 2, 3, 4, 5]].copy()
    with pytest.raises(ValueError, match="observation"):
        tf_activity_plot_owned(reordered, activities=activities)

    result = _table_result(science, "scenic_binarization")
    result.table.loc[0, "threshold"] = 0.9
    with pytest.raises(ValueError, match="fingerprint"):
        scenic_binarization_plot_owned(result)


@pytest.mark.parametrize("kind", ["tf", "scenic"])
def test_regulatory_activity_embedding_uses_first_stored_activity_by_default(science, kind):
    adata = _adata(science)
    artifact = _tf_artifact(science, adata) if kind == "tf" else _scenic_artifact(science, adata)
    owned = tf_activity_plot_owned if kind == "tf" else scenic_activity_plot_owned
    kwargs = {"activities": artifact} if kind == "tf" else {"scenic_result": artifact}

    plotted, report, _code = owned(
        adata,
        **kwargs,
        view={"view": "embedding", "embedding_key": "X_umap", "activity_name": ""},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    assert report.summary["parameters"]["view"]["activity_name"] in {"TF1", "TF1(+)"}
    assert report.summary["key_results"]["plotted_activities"] == 1


@pytest.mark.parametrize("kind", ["tf", "scenic"])
def test_regulatory_activity_plot_workers_use_typed_codecs_without_rewriting_inputs(
    tmp_path, science, kind
):
    from openbio_singlecell.operations_regulatory import (
        scenic_activity_plot_operation,
        tf_activity_plot_operation,
    )

    adata = _adata(science)
    artifact = _tf_artifact(science, adata) if kind == "tf" else _scenic_artifact(science, adata)
    adata_root = tmp_path / "adata"
    artifact_root = tmp_path / kind
    adata_root.mkdir()
    artifact_root.mkdir()
    write_anndata(adata_root, adata)
    if kind == "tf":
        write_tf_activity(artifact_root, artifact)
        operation = tf_activity_plot_operation
        input_name, artifact_kind, codec = "activities", "OPENBIO_TF_ACTIVITY", TF_ACTIVITY_CODEC
    else:
        write_scenic(artifact_root, artifact)
        operation = scenic_activity_plot_operation
        input_name, artifact_kind, codec = "scenic_result", "OPENBIO_SCENIC_RESULT", SCENIC_CODEC
    before_adata = hashlib.sha256((adata_root / ANNDATA_PAYLOAD).read_bytes()).digest()
    before_artifact = {
        str(path.relative_to(artifact_root)): hashlib.sha256(path.read_bytes()).digest()
        for path in artifact_root.rglob("*")
        if path.is_file()
    }
    staging = tmp_path / f"{kind}-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = operation(
        context,
        {
            "adata": {"type": "artifact", "kind": "OPENBIO_ANNDATA", "codec": "anndata-h5ad-v1", "path": str(adata_root.resolve())},
            input_name: {"type": "artifact", "kind": artifact_kind, "codec": codec, "path": str(artifact_root.resolve())},
        },
        {"view": {"view": "distributions", "max_activities": 20}},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert hashlib.sha256((adata_root / ANNDATA_PAYLOAD).read_bytes()).digest() == before_adata
    assert before_artifact == {
        str(path.relative_to(artifact_root)): hashlib.sha256(path.read_bytes()).digest()
        for path in artifact_root.rglob("*")
        if path.is_file()
    }


@pytest.mark.parametrize(
    ("contract", "operation_name", "view"),
    [
        ("tf_ranking", "tf_activity_ranking_plot_operation", {"view": "ranked_dot", "max_per_group": 10}),
        ("scenic_rss", "scenic_rss_plot_operation", {"view": "rss_heatmap", "max_regulons": 10}),
        ("scenic_binarization", "scenic_binarization_plot_operation", {"view": "active_fraction", "max_regulons": 10}),
        ("scenic_membership", "scenic_membership_plot_operation", {"view": "target_count", "max_regulons": 10}),
    ],
)
def test_regulatory_table_plot_workers_round_trip_each_public_adapter(
    tmp_path, science, contract, operation_name, view
):
    import openbio_singlecell.operations_regulatory as operations

    result = _table_result(science, contract)
    input_root = tmp_path / contract
    input_root.mkdir()
    list_columns = ("motif_ids",) if contract == "scenic_membership" else ()
    write_table(input_root, result.table, result_metadata(result), json_list_columns=list_columns)
    before = {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}
    staging = tmp_path / f"{contract}-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = getattr(operations, operation_name)(
        context,
        {"table": {"type": "artifact", "kind": "OPENBIO_SINGLE_CELL_TABLE", "codec": "table-jsonl-v1", "path": str(input_root.resolve())}},
        {"view": view},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert before == {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}
