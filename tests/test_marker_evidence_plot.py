from __future__ import annotations

import hashlib
import json
import time
import uuid
from dataclasses import replace

import pytest

from openbio_singlecell.analysis_utils import make_table_result
from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.marker_evidence import (
    MARKER_COLUMNS,
    MARKER_UNIVERSE_COLUMNS,
    _marker_analysis_fingerprint,
    _marker_table_content_fingerprint,
    marker_provenance_parameters,
    marker_table_content_fingerprint,
    marker_universe_content_fingerprint,
)
from openbio_singlecell.nodes_results import OpenBioSingleCellMarkerEvidencePlot
from openbio_singlecell.operations_results import marker_evidence_plot, marker_evidence_plot_owned
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _marker_inputs(science):
    obs_names = [f"cell_{index}" for index in range(8)]
    genes = ["G0", "G1", "G2", "G3"]
    labels = ["A"] * 4 + ["B"] * 4
    values = science.np.asarray(
        [
            [8.0, 1.0, 4.0, 0.0],
            [7.0, 0.0, 3.0, 0.0],
            [6.0, 1.0, 5.0, 1.0],
            [9.0, 0.0, 4.0, 0.0],
            [1.0, 8.0, 0.0, 4.0],
            [0.0, 7.0, 1.0, 3.0],
            [1.0, 9.0, 0.0, 5.0],
            [0.0, 6.0, 1.0, 4.0],
        ]
    )
    adata = science.ad.AnnData(
        values.copy(),
        obs=science.pd.DataFrame(
            {"cluster": science.pd.Categorical(labels, categories=["A", "B"], ordered=True)},
            index=obs_names,
        ),
        var=science.pd.DataFrame(index=genes),
    )
    adata.layers["log1p_norm"] = values.copy()
    table_frame = science.pd.DataFrame(
        [
            ("A", "G0", 1, 5.0, 2.0, 0.001, 0.004, 1.00, 0.50),
            ("A", "G2", 2, 4.0, 1.5, 0.002, 0.006, 1.00, 0.50),
            ("A", "G1", 3, 1.0, -1.0, 0.20, 0.30, 0.50, 1.00),
            ("A", "G3", 4, 0.5, -0.5, 0.40, 0.50, 0.25, 1.00),
            ("B", "G1", 1, 5.5, 2.1, 0.001, 0.004, 1.00, 0.50),
            ("B", "G3", 2, 4.5, 1.7, 0.002, 0.006, 1.00, 0.25),
            ("B", "G0", 3, 0.7, -1.2, 0.20, 0.30, 0.50, 1.00),
            ("B", "G2", 4, 0.4, -0.8, 0.40, 0.50, 0.50, 1.00),
        ],
        columns=MARKER_COLUMNS,
    )
    universe_frame = science.pd.DataFrame(
        {"gene": genes, "universe_rank": [1, 2, 3, 4]},
        columns=MARKER_UNIVERSE_COLUMNS,
    )
    analysis_fingerprint, universe_fingerprint = _marker_analysis_fingerprint(
        groupby="cluster",
        method="wilcoxon",
        source_kind="layer",
        layer_name="log1p_norm",
        n_genes=0,
        tie_correct=True,
        genes=genes,
        observation_ids=obs_names,
        group_labels_by_observation=labels,
    )
    details = {
        "analysis_fingerprint": analysis_fingerprint,
        "ranking_fingerprint": _marker_table_content_fingerprint(
            table_frame, artifact="complete-marker-ranking"
        ),
        "universe_fingerprint": universe_fingerprint,
        "table_content_fingerprint": marker_table_content_fingerprint(table_frame),
        "universe_content_fingerprint": marker_universe_content_fingerprint(universe_frame),
        "group_labels": ["A", "B"],
        "group_sizes": {"A": 4, "B": 4},
        "source_gene_count": 4,
        "actual_n_genes_per_group": 4,
        "ranking_truncated": False,
        "constant_gene_count": 0,
        "variable_gene_count": 4,
        "bh_hypotheses_per_group": 4,
        "marker_output_rows": 8,
        "universe_output_rows": 4,
        "total_output_rows": 12,
    }
    common = {
        "groupby": "cluster",
        "method": "wilcoxon",
        "source_kind": "layer",
        "layer_name": "log1p_norm",
        "n_genes": 0,
        "tie_correct": True,
        "max_output_rows": 10_000,
        "max_working_memory_gib": 4.0,
        "details": details,
    }
    table_parameters = marker_provenance_parameters(artifact_role="marker_table", **common)
    universe_parameters = marker_provenance_parameters(artifact_role="tested_gene_universe", **common)
    table = make_table_result(
        table=table_frame,
        title="Cluster marker evidence",
        operation="marker_genes",
        parameters=table_parameters,
        description="Exploratory Cluster marker evidence.",
        warnings=[],
        input_cells=8,
        input_genes=4,
        started_at=time.perf_counter(),
    )
    universe = make_table_result(
        table=universe_frame,
        title="Tested-gene universe",
        operation="marker_genes",
        parameters=universe_parameters,
        description="Exact tested-gene universe.",
        warnings=[],
        input_cells=8,
        input_genes=4,
        started_at=time.perf_counter(),
    )
    return adata, table, universe


def _replace_marker_ranks(_adata, table, _universe):
    table.table.loc[:, "rank"] = [2, 1, 3, 4, 1, 2, 3, 4]


def test_marker_evidence_plot_schema_uses_exact_evidence_and_expression_source():
    schema = OpenBioSingleCellMarkerEvidencePlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellMarkerEvidencePlot"
    assert schema.display_name == "Marker Evidence Plot"
    assert schema.category == "openbio/single-cell/marker-evidence"
    assert [item.id for item in schema.inputs] == [
        "adata",
        "table",
        "universe",
        "top_genes_per_group",
        "source",
    ]
    assert [option.key for option in schema.inputs[-1].options] == ["layer", "X"]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_marker_evidence_plot_uses_upstream_rank_without_reranking(science):
    adata, table, universe = _marker_inputs(science)
    before = adata.copy()
    table_before = table.table.copy(deep=True)
    universe_before = universe.table.copy(deep=True)

    plotted, report, code = marker_evidence_plot_owned(
        adata,
        table,
        universe,
        top_genes_per_group=2,
        source={"source": "layer", "layer_name": "log1p_norm"},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["genes"] == ["G0", "G2", "G1", "G3"]
    assert details["selected_upstream_ranks"] == {
        "A": [{"gene": "G0", "rank": 1}, {"gene": "G2", "rank": 2}],
        "B": [{"gene": "G1", "rank": 1}, {"gene": "G3", "rank": 2}],
    }
    assert details["display_semantics"] == {
        "color": "mean selected-source expression",
        "size": "fraction of cells with expression > 0",
    }
    assert details["size_legend_fractions"] == [0.25, 0.5, 1.0]
    assert details["size_legend_location"] == "outside lower center"
    assert details["status"] == "cluster_marker_evidence_plot"
    text = " ".join([report.summary["results"], *report.summary["limitations"]])
    assert "Condition contrast" in text
    assert "Curated annotation" in text
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<marker-evidence-plot-code>", "exec"), namespace)
    assert namespace["plot_marker_evidence"](adata, table.table, universe.table) == plotted.png
    science.np.testing.assert_array_equal(adata.layers["log1p_norm"], before.layers["log1p_norm"])
    science.pd.testing.assert_frame_equal(adata.obs, before.obs)
    assert adata.uns == before.uns
    science.pd.testing.assert_frame_equal(table.table, table_before)
    science.pd.testing.assert_frame_equal(universe.table, universe_before)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda adata, table, universe: table.source.__setitem__("operation", "condition_contrast"), "operation"),
        (_replace_marker_ranks, "rank"),
        (lambda adata, table, universe: adata.__setattr__("obs_names", list(reversed(adata.obs_names))), "axis"),
    ],
)
def test_marker_evidence_plot_rejects_wrong_or_tampered_evidence(science, mutation, message):
    adata, table, universe = _marker_inputs(science)
    mutation(adata, table, universe)

    with pytest.raises((TypeError, ValueError), match=message):
        marker_evidence_plot_owned(
            adata,
            table,
            universe,
            top_genes_per_group=2,
            source={"source": "layer", "layer_name": "log1p_norm"},
        )


def test_marker_evidence_plot_rejects_expression_source_different_from_marker_analysis(science):
    adata, table, universe = _marker_inputs(science)

    with pytest.raises(ValueError, match="expression source.*upstream"):
        marker_evidence_plot_owned(
            adata,
            table,
            universe,
            top_genes_per_group=2,
            source={"source": "X"},
        )


def test_marker_evidence_plot_rejects_mismatched_table_and_universe_axes(science):
    adata, table, universe = _marker_inputs(science)

    with pytest.raises(ValueError, match="table and universe input axes"):
        marker_evidence_plot_owned(
            adata,
            table,
            replace(universe, input_cells=7),
            top_genes_per_group=2,
            source={"source": "layer", "layer_name": "log1p_norm"},
        )


def test_marker_evidence_plot_worker_preserves_all_three_input_artifacts(tmp_path, science):
    adata, table, universe = _marker_inputs(science)
    adata_root = tmp_path / "adata"
    table_root = tmp_path / "table"
    universe_root = tmp_path / "universe"
    for root in (adata_root, table_root, universe_root):
        root.mkdir()
    write_anndata(adata_root, adata)
    write_table(table_root, table.table, result_metadata(table))
    write_table(universe_root, universe.table, result_metadata(universe))
    before_adata = hashlib.sha256((adata_root / ANNDATA_PAYLOAD).read_bytes()).digest()
    before_table = {path.name: hashlib.sha256(path.read_bytes()).digest() for path in table_root.iterdir()}
    before_universe = {
        path.name: hashlib.sha256(path.read_bytes()).digest() for path in universe_root.iterdir()
    }
    staging = tmp_path / "marker-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = marker_evidence_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(adata_root.resolve()),
            },
            "table": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
                "path": str(table_root.resolve()),
            },
            "universe": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
                "path": str(universe_root.resolve()),
            },
        },
        {
            "top_genes_per_group": 2,
            "source": {"source": "layer", "layer_name": "log1p_norm"},
        },
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"]["top_genes_per_group"] == 2
    assert hashlib.sha256((adata_root / ANNDATA_PAYLOAD).read_bytes()).digest() == before_adata
    assert before_table == {
        path.name: hashlib.sha256(path.read_bytes()).digest() for path in table_root.iterdir()
    }
    assert before_universe == {
        path.name: hashlib.sha256(path.read_bytes()).digest() for path in universe_root.iterdir()
    }


def test_marker_evidence_plot_restores_matplotlib_global_state(science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    adata, table, universe = _marker_inputs(science)
    matplotlib.rcParams["lines.linewidth"] = 3.25
    before_rc = matplotlib.rcParams["lines.linewidth"]
    before_figures = pyplot.get_fignums()

    marker_evidence_plot_owned(
        adata,
        table,
        universe,
        top_genes_per_group=2,
        source={"source": "layer", "layer_name": "log1p_norm"},
    )

    assert matplotlib.rcParams["lines.linewidth"] == before_rc
    assert pyplot.get_fignums() == before_figures
