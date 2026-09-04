from __future__ import annotations

import hashlib
import json
import uuid

import pytest

from openbio_singlecell.artifact_codecs import ANNDATA_PAYLOAD, read_plot, write_anndata
from openbio_singlecell.contracts import record_history
from openbio_singlecell.nodes_schist import OpenBioSingleCellSchistHierarchyPlot
from openbio_singlecell.operations_schist import schist_hierarchy_plot, schist_hierarchy_plot_owned
from openbio_singlecell.schist_analysis import _standalone_hierarchy_evidence_fingerprint
from openbio_singlecell.worker_protocol import OperationContext

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _schist_result(science):
    adata = science.ad.AnnData(
        science.np.zeros((6, 2)),
        obs=science.pd.DataFrame(index=[f"cell-{index}" for index in range(6)]),
        var=science.pd.DataFrame(index=["gene-a", "gene-b"]),
    )
    memberships = (
        science.np.asarray([0, 0, 1, 1, 2, 2], dtype=int),
        science.np.asarray([0, 0, 0, 0, 1, 1], dtype=int),
        science.np.asarray([0, 0, 0, 0, 0, 0], dtype=int),
    )
    for level, values in enumerate(memberships):
        adata.obs[f"nsbm_level_{level}"] = science.pd.Categorical(values.astype(str))
    adata.obsm["CM_nsbm_level_0"] = science.np.asarray(
        [
            [0.9, 0.05, 0.05],
            [0.8, 0.1, 0.1],
            [0.05, 0.9, 0.05],
            [0.1, 0.8, 0.1],
            [0.05, 0.05, 0.9],
            [0.1, 0.1, 0.8],
        ]
    )
    adata.obsm["CM_nsbm_level_1"] = science.np.asarray(
        [[0.9, 0.1], [0.8, 0.2], [0.9, 0.1], [0.8, 0.2], [0.1, 0.9], [0.2, 0.8]]
    )
    adata.obsm["CM_nsbm_level_2"] = science.np.ones((6, 1))
    adata.uns["schist"] = {
        "nsbm": {
            "stats": {
                "entropy": 12.5,
                "modularity": science.np.asarray([0.42, 0.21, 0.0]),
                "level_entropy": science.np.asarray([5.0, 2.0, 0.0]),
            },
            "blocks": {
                "0": science.np.asarray([0, 0, 1, 1, 2, 2], dtype=int),
                "1": science.np.asarray([0, 0, 1], dtype=int),
                "2": science.np.asarray([0, 0], dtype=int),
            },
            "params": {"nested": True, "collect_marginals": True, "key_added": "nsbm"},
        }
    }
    adata.uns["schist"]["nsbm"]["openbio_evidence_schema"] = (
        "openbio-singlecell/schist-hierarchy-evidence/v1"
    )
    adata.uns["schist"]["nsbm"]["openbio_evidence_sha256"] = _standalone_hierarchy_evidence_fingerprint(
        adata,
        key_added="nsbm",
        np=science.np,
    )
    record_history(
        adata,
        "schist_nested_sbm_hierarchy",
        {"key_added": "nsbm"},
        input_cells=6,
        input_genes=2,
        random_seed=123,
    )
    return adata


def test_schist_hierarchy_plot_schema_and_stored_hierarchy(science):
    schema = OpenBioSingleCellSchistHierarchyPlot.define_schema()
    assert schema.category == "openbio/single-cell/clustering"
    assert [item.id for item in schema.inputs] == ["adata", "key_added"]
    assert [item.display_name for item in schema.outputs] == ["plot", "summary", "code"]

    adata = _schist_result(science)
    before = adata.copy()

    plotted, report, code = schist_hierarchy_plot_owned(adata, key_added="nsbm")

    assert plotted.png.startswith(PNG_SIGNATURE)
    details = report.summary["key_results"]
    assert details["cluster_counts_finest_to_root"] == [3, 2, 1]
    assert details["parent_edges"] == 5
    assert details["level_entropy"] == [5.0, 2.0, 0.0]
    assert details["modularity"] == [0.42, 0.21, 0.0]
    assert details["cluster_sizes_record_limit"] == 50
    assert details["cluster_sizes_truncated"] == [False, False, False]
    json.dumps(report.summary, allow_nan=False)
    namespace = {}
    exec(code, namespace)
    assert namespace["plot_schist_hierarchy"](adata) == plotted.png
    assert adata.uns["schist"]["nsbm"]["params"] == before.uns["schist"]["nsbm"]["params"]
    for key in ("entropy", "modularity", "level_entropy"):
        science.np.testing.assert_array_equal(
            adata.uns["schist"]["nsbm"]["stats"][key],
            before.uns["schist"]["nsbm"]["stats"][key],
        )
    for key in ("0", "1", "2"):
        science.np.testing.assert_array_equal(
            adata.uns["schist"]["nsbm"]["blocks"][key],
            before.uns["schist"]["nsbm"]["blocks"][key],
        )
    for key in adata.obsm:
        science.np.testing.assert_array_equal(adata.obsm[key], before.obsm[key])


def test_schist_hierarchy_plot_rejects_non_nested_membership(science):
    adata = _schist_result(science)
    adata.obs["nsbm_level_1"] = science.pd.Categorical(["0", "0", "0", "1", "1", "1"])

    with pytest.raises(ValueError, match="nested parent mapping"):
        schist_hierarchy_plot_owned(adata, key_added="nsbm")


def test_schist_hierarchy_plot_requires_matching_openbio_producer(science):
    adata = _schist_result(science)
    del adata.uns["openbio_singlecell"]

    with pytest.raises(ValueError, match="producer history"):
        schist_hierarchy_plot_owned(adata, key_added="nsbm")


def test_schist_hierarchy_plot_rejects_valid_range_evidence_tampering(science):
    adata = _schist_result(science)
    adata.uns["schist"]["nsbm"]["stats"]["modularity"][0] = 0.41

    with pytest.raises(ValueError, match="evidence fingerprint"):
        schist_hierarchy_plot_owned(adata, key_added="nsbm")


def test_schist_hierarchy_plot_worker_round_trips_without_rewriting_input(tmp_path, science):
    adata = _schist_result(science)
    input_root = tmp_path / "input"
    input_root.mkdir()
    write_anndata(input_root, adata)
    input_path = input_root / ANNDATA_PAYLOAD
    before = hashlib.sha256(input_path.read_bytes()).hexdigest()
    staging = tmp_path / "schist-plot.partial"
    staging.mkdir()
    context = OperationContext.from_request_path(staging / "request.json", str(uuid.uuid4()))

    records = schist_hierarchy_plot(
        context,
        {
            "adata": {
                "type": "artifact",
                "kind": "OPENBIO_ANNDATA",
                "codec": "anndata-h5ad-v1",
                "path": str(input_root.resolve()),
            }
        },
        {"key_added": "nsbm"},
    )

    assert [record["name"] for record in records] == ["plot", "summary", "code"]
    png, _ = read_plot(staging / records[0]["payload"])
    assert png.startswith(PNG_SIGNATURE)
    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == before
