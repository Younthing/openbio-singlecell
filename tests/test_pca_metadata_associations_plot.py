from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import read_plot, write_table
from openbio_singlecell.artifact_envelope import result_metadata
from openbio_singlecell.nodes_results import OpenBioSingleCellPCAMetadataAssociationsPlot
from openbio_singlecell.operations_results import (
    pca_metadata_associations_owned,
    pca_metadata_associations_plot,
    pca_metadata_associations_plot_owned,
)
from openbio_singlecell.worker_protocol import OperationContext, WorkerResponse

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _association_result(science):
    sample_names = [f"sample_{index}" for index in range(8)]
    conditions = ["control"] * 4 + ["treated"] * 4
    ages = [20.0, 28.0, 35.0, 43.0, 51.0, 60.0, 68.0, 77.0]
    pc1 = [0.0, 1.5, 0.5, 2.5, 5.0, 8.0, 6.5, 10.0]
    pc2 = [8.0, 6.0, 7.0, 5.0, 4.0, 2.0, 3.0, 1.0]
    adata = science.ad.AnnData(
        science.np.ones((8, 2)),
        obs=science.pd.DataFrame(
            {"sample": sample_names, "condition": conditions, "age": ages},
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=science.pd.DataFrame(index=["G0", "G1"]),
    )
    adata.obsm["X_pca"] = science.np.column_stack([pc1, pc2])
    return pca_metadata_associations_owned(
        adata,
        use_rep="X_pca",
        sample_key="sample",
        categorical_obs_keys="condition",
        continuous_obs_keys="age",
        alpha=0.05,
    )[0]


def _replace_effect_sizes(result):
    result.table.loc[:, "effect_size"] = -1.0


def _replace_adjusted_p_values(result):
    result.table.loc[:, "p_adjusted"] = 0.0


def test_pca_metadata_associations_plot_schema_has_closed_stored_evidence_views():
    schema = OpenBioSingleCellPCAMetadataAssociationsPlot.define_schema()

    assert schema.node_id == "OpenBioSingleCellPCAMetadataAssociationsPlot"
    assert schema.display_name == "PCA Metadata Associations Plot"
    assert schema.category == "openbio/single-cell/diagnostics"
    assert [item.id for item in schema.inputs] == ["table", "view"]
    assert [option.key for option in schema.inputs[1].options] == [
        "association_heatmap",
        "effect_sizes",
    ]
    assert [[item.id for item in option.inputs] for option in schema.inputs[1].options] == [
        [],
        ["max_associations"],
    ]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]


def test_pca_metadata_association_heatmap_uses_canonical_table_without_retesting(science):
    table = _association_result(science)
    before = table.table.copy(deep=True)

    plotted, report, code = pca_metadata_associations_plot_owned(
        table,
        view={"view": "association_heatmap"},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["view"] == "association_heatmap"
    assert details["components"] == ["PC1", "PC2"]
    assert details["metadata"] == ["condition", "age"]
    assert details["plotted_associations"] == 4
    assert details["effect_size_semantics"] == {
        "categorical": "eta_squared",
        "continuous": "rho_squared",
    }
    assert details["significance_annotation"] == "* BH-adjusted p ≤ α (α=0.05)"
    assert "not causal" in report.summary["results"]
    assert "Condition inference" in " ".join(report.summary["limitations"])
    assert "recomputed" not in report.summary["methods"].lower()
    assert "stored" in report.summary["methods"].lower()
    json.dumps(report.summary, allow_nan=False)

    namespace: dict[str, object] = {}
    exec(compile(code, "<pca-metadata-plot-code>", "exec"), namespace)
    assert namespace["plot_pca_metadata_associations"](table.table) == plotted.png
    science.pd.testing.assert_frame_equal(table.table, before)


def test_pca_metadata_effect_view_preserves_canonical_evidence_order(science):
    table = _association_result(science)

    plotted, report, code = pca_metadata_associations_plot_owned(
        table,
        view={"view": "effect_sizes", "max_associations": 3},
    )

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["view"] == "effect_sizes"
    assert details["plotted_associations"] == 3
    expected = [
        {"component": str(row.component), "metadata": str(row.metadata)}
        for row in table.table.iloc[:3].itertuples()
    ]
    assert [
        {"component": row["component"], "metadata": row["metadata"]}
        for row in details["displayed_rows"]
    ] == expected
    namespace: dict[str, object] = {}
    exec(compile(code, "<pca-metadata-effect-plot-code>", "exec"), namespace)
    assert namespace["plot_pca_metadata_associations"](table.table) == plotted.png


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda result: result.source.__setitem__("operation", "marker_genes"), "producer"),
        (_replace_effect_sizes, "effect size"),
        (_replace_adjusted_p_values, "adjusted p-value"),
    ],
)
def test_pca_metadata_plot_rejects_wrong_or_tampered_table(science, mutation, message):
    table = _association_result(science)
    mutation(table)

    with pytest.raises((TypeError, ValueError), match=message):
        pca_metadata_associations_plot_owned(table)


def test_pca_metadata_plot_rejects_inactive_view_fields(science):
    table = _association_result(science)

    with pytest.raises(ValueError, match="inactive or unknown"):
        pca_metadata_associations_plot_owned(
            table,
            view={"view": "association_heatmap", "max_associations": 3},
        )


def test_pca_metadata_plot_worker_preserves_the_input_table_artifact(tmp_path, science):
    table = _association_result(science)
    input_root = tmp_path / "pca-metadata"
    input_root.mkdir()
    write_table(input_root, table.table, result_metadata(table))
    before = {path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()}
    staging = tmp_path / "pca-metadata-plot.partial"
    staging.mkdir()
    context = OperationContext(staging, str(uuid.uuid4()))

    records = pca_metadata_associations_plot(
        context,
        {
            "table": {
                "type": "artifact",
                "kind": "OPENBIO_SINGLE_CELL_TABLE",
                "codec": "table-jsonl-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"view": {"view": "association_heatmap"}},
    )
    WorkerResponse.success(context.request_id, records)

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, metadata = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert metadata["parameters"] == {"view": {"view": "association_heatmap"}}
    assert before == {
        path.name: hashlib.sha256(path.read_bytes()).digest() for path in input_root.iterdir()
    }


def test_pca_metadata_plot_restores_matplotlib_global_state(science):
    import matplotlib
    import matplotlib.pyplot as pyplot

    matplotlib.rcParams["lines.linewidth"] = 3.25
    before_rc = matplotlib.rcParams["lines.linewidth"]
    before_figures = pyplot.get_fignums()

    pca_metadata_associations_plot_owned(_association_result(science))

    assert matplotlib.rcParams["lines.linewidth"] == before_rc
    assert pyplot.get_fignums() == before_figures
