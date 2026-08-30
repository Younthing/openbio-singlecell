from __future__ import annotations

import copy
import json

import pytest

from openbio_singlecell.node_types import SummaryResultType, TableResultType
from openbio_singlecell.nodes_abundance import (
    OpenBioSingleCellSampleCompositionSummary,
)
from openbio_singlecell.operations_abundance import sample_composition_owned


def _composition_adata(science):
    obs = science.pd.DataFrame(
        {
            "sample": ["s1", "s1", "s2", "s2", "s3"],
            "condition": ["control", "control", "control", "control", "treated"],
            "cell_type": science.pd.Categorical(
                ["A", "A", "B", "C", "A"],
                categories=["B", "A", "C", "unused"],
                ordered=True,
            ),
        },
        index=[f"cell_{index}" for index in range(5)],
    )
    adata = science.ad.AnnData(
        science.np.arange(10, dtype=float).reshape(5, 2),
        obs=obs,
        var=science.pd.DataFrame(index=["G1", "G2"]),
    )
    adata.layers["counts"] = science.np.arange(10, dtype=int).reshape(5, 2)
    adata.uns["fixture"] = {"preserved": True}
    return adata


def test_sample_composition_schema():
    schema = OpenBioSingleCellSampleCompositionSummary.define_schema()
    assert [item.id for item in schema.inputs] == [
        "adata",
        "sample_key",
        "condition_key",
        "annotation_key",
        "annotation_status",
        "max_output_rows",
    ]
    assert [(item.display_name, item.io_type) for item in schema.outputs] == [
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    inputs = {item.id: item for item in schema.inputs}
    assert inputs["condition_key"].default == "condition"
    assert inputs["annotation_status"].default == "unknown"
    assert inputs["max_output_rows"].default == 2_000_000



def test_sample_composition_complete_grid_summary_code_and_input_immutability(science):
    adata = _composition_adata(science)
    original_obs = adata.obs.copy(deep=True)
    original_x = adata.X.copy()
    original_counts = adata.layers["counts"].copy()
    original_uns = copy.deepcopy(adata.uns)

    result, report, code = sample_composition_owned(
        adata,
        annotation_status="provisional",
    )
    table = result.table
    assert table.columns.tolist() == [
        "sample",
        "condition",
        "annotation",
        "cell_count",
        "sample_total_cells",
        "proportion",
    ]
    assert table[["sample", "annotation"]].to_dict("records") == [
        {"sample": sample, "annotation": annotation}
        for sample in ["s1", "s2", "s3"]
        for annotation in ["B", "A", "C"]
    ]
    assert table["cell_count"].tolist() == [0, 2, 0, 1, 0, 1, 0, 1, 0]
    assert table["sample_total_cells"].tolist() == [2, 2, 2, 2, 2, 2, 1, 1, 1]
    assert table.groupby("sample", sort=False)["proportion"].sum().tolist() == [1.0, 1.0, 1.0]
    assert int(table["cell_count"].sum()) == adata.n_obs

    summary = report.summary
    json.dumps(summary, allow_nan=False)
    assert summary["node_id"] == "OpenBioSingleCellSampleCompositionSummary"
    assert summary["key_results"]["complete_grid_rows"] == 9
    assert summary["key_results"]["zero_count_rows"] == 5
    assert summary["key_results"]["unused_declared_annotations"] == ["unused"]
    assert summary["key_results"]["samples_per_condition"] == {"control": 2, "treated": 1}
    assert summary["key_results"]["annotation_status"] == "provisional"
    assert "no Condition contrast was performed" in summary["results"]
    assert set(summary["software_versions"]) == {
        "python",
        "openbio-singlecell",
        "anndata",
        "pandas",
        "numpy",
    }
    compile(code, "<sample-composition-code>", "exec")
    namespace = {}
    exec(code, namespace)
    reproduced, reproduced_summary = namespace["summarize_sample_composition"](adata)
    science.pd.testing.assert_frame_equal(reproduced, table)
    assert reproduced_summary == summary

    science.pd.testing.assert_frame_equal(adata.obs, original_obs)
    science.np.testing.assert_array_equal(adata.X, original_x)
    science.np.testing.assert_array_equal(adata.layers["counts"], original_counts)
    assert adata.uns == original_uns


@pytest.mark.parametrize("column", ["sample", "condition", "cell_type"])
def test_sample_composition_missing_metadata_fails_instead_of_dropping_cells(science, column):
    adata = _composition_adata(science)
    adata.obs[column] = adata.obs[column].astype(object)
    adata.obs.loc[adata.obs.index[0], column] = None

    with pytest.raises(ValueError, match="contains a missing value"):
        sample_composition_owned(adata)


def test_sample_composition_rejects_whitespace_collisions_and_multi_condition_samples(science):
    whitespace = _composition_adata(science)
    whitespace.obs["sample"] = whitespace.obs["sample"].astype(object)
    whitespace.obs.loc[whitespace.obs.index[0], "sample"] = " s1"
    with pytest.raises(ValueError, match="whitespace-padded"):
        sample_composition_owned(whitespace)

    collision = _composition_adata(science)
    collision.obs["sample"] = collision.obs["sample"].astype(object)
    collision.obs.iloc[0, collision.obs.columns.get_loc("sample")] = 1
    collision.obs.iloc[1, collision.obs.columns.get_loc("sample")] = "1"
    with pytest.raises(ValueError, match="collide after display normalization"):
        sample_composition_owned(collision)

    ambiguous = _composition_adata(science)
    ambiguous.obs.loc[ambiguous.obs.index[1], "condition"] = "treated"
    ambiguous.obs["cell_type"] = ambiguous.obs["cell_type"].astype(object)
    ambiguous.obs.loc[ambiguous.obs.index[1], "cell_type"] = None
    with pytest.raises(ValueError, match="one Condition per Sample"):
        sample_composition_owned(ambiguous)


def test_sample_composition_grid_guard_and_minimal_descriptive_input(science):
    adata = _composition_adata(science)
    with pytest.raises(ValueError, match=r"requires 9 rows.*max_output_rows=8"):
        sample_composition_owned(adata, max_output_rows=8)

    minimal = science.ad.AnnData(
        science.np.empty((2, 0)),
        obs=science.pd.DataFrame(
            {"sample": ["s1", "s1"], "condition": ["only", "only"], "cell_type": ["A", "A"]},
            index=["c1", "c2"],
        ),
        var=science.pd.DataFrame(index=[]),
    )
    result, report, _ = sample_composition_owned(minimal)
    assert result.table.to_dict("records") == [
        {
            "sample": "s1",
            "condition": "only",
            "annotation": "A",
            "cell_count": 2,
            "sample_total_cells": 2,
            "proportion": 1.0,
        }
    ]
    assert report.summary["key_results"]["samples"] == 1
    assert any("descriptive" in warning for warning in report.summary["warnings"])

