from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import openbio_singlecell.milo_result as milo_result_module
from openbio_singlecell.abundance_artifact_codecs import (
    COMPOSITION_MODEL_CODEC,
    COMPOSITION_MODEL_KIND,
    MILO_RESULT_CODEC,
    MILO_RESULT_KIND,
    read_composition_model_result,
    read_milo_result,
    write_composition_model_result,
    write_milo_result,
)
from openbio_singlecell.artifact_codecs import TABLE_PAYLOAD
from openbio_singlecell.artifact_persist import persist_artifact
from openbio_singlecell.artifact_runtime import ArtifactRuntime
from openbio_singlecell.composition_result import (
    SCCODA_RESULT_COLUMNS,
    build_composition_model_result,
    validate_composition_model_result,
)
from openbio_singlecell.milo_analysis import MILO_TABLE_COLUMNS
from openbio_singlecell.milo_result import build_milo_result, validate_milo_result
from openbio_singlecell.node_types import CompositionModelResultType, MiloResultType, SummaryResultType, TableResultType
from openbio_singlecell.nodes_abundance import (
    OpenBioSingleCellMiloDifferentialAbundance,
    OpenBioSingleCellSccodaDifferentialComposition,
    OpenBioSingleCellTasccodaDifferentialComposition,
)


def _milo_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ["nhood_000001", "c1", 2, 0.5, "A", "B", -1.0, 4.0, 8.0, 0.01, 0.02, 0.03, "T", 1.0, "T", False],
            ["nhood_000002", "c3", 2, 0.7, "A", "B", 0.8, 5.0, 7.0, 0.02, 0.03, 0.04, "B", 1.0, "B", False],
        ],
        columns=MILO_TABLE_COLUMNS,
    )


def _milo_result():
    membership = sparse.csr_matrix([[1, 0], [1, 1], [0, 1]], dtype=np.int8)
    graph = membership.T @ membership
    graph.setdiag(0)
    graph.eliminate_zeros()
    return build_milo_result(
        table=_milo_table(),
        membership=membership,
        graph=graph,
        representative_coordinates=np.array([[0.0, 1.0], [2.0, 3.0]]),
        observation_names=["c1", "c2", "c3"],
        neighborhood_names=["nhood_000001", "nhood_000002"],
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


def test_milo_result_is_immutable_by_interface_and_validates_all_axes():
    result = _milo_result()
    table, membership, graph, coordinates, observations, neighborhoods, provenance, metadata = (
        validate_milo_result(result)
    )

    assert observations == ["c1", "c2", "c3"]
    assert neighborhoods == ["nhood_000001", "nhood_000002"]
    assert provenance["representation_key"] == "X_pca"
    assert provenance["annotation_status"] == "curated"
    assert provenance["n_neighbors"] == 3
    assert metadata["artifact_fingerprint_sha256"] == result.fingerprint
    assert table["neighborhood_size"].tolist() == [2, 2]
    assert membership.shape == (3, 2)
    assert graph.toarray().tolist() == [[0, 1], [1, 0]]
    assert coordinates.tolist() == [[0.0, 1.0], [2.0, 3.0]]

    table.iloc[0, table.columns.get_loc("spatial_fdr")] = 0.9
    membership[0, 0] = 0
    coordinates[0, 0] = 99.0
    assert result.table.loc[0, "spatial_fdr"] == 0.03
    assert result.membership[0, 0] == 1
    assert result.representative_coordinates[0, 0] == 0.0
    with pytest.raises(AttributeError, match="immutable"):
        result.fingerprint = "b" * 64


def test_milo_result_preflights_overlap_graph_before_sparse_multiplication(monkeypatch):
    monkeypatch.setattr(milo_result_module, "MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS", 0)

    with pytest.raises(MemoryError, match="overlap-graph pair contributions"):
        _milo_result()


def test_milo_result_codec_round_trips_and_rejects_tampered_table(tmp_path):
    result = _milo_result()

    descriptors = write_milo_result(tmp_path, result)
    restored = read_milo_result(tmp_path)

    assert MILO_RESULT_CODEC == "milo-result-h5ad-v1"
    assert {item["path"] for item in descriptors} == {
        "data.h5ad",
        "metadata.json",
        "table/data.jsonl",
        "table/result.json",
        "table/schema.json",
    }
    pd.testing.assert_frame_equal(restored.table, result.table)
    np.testing.assert_array_equal(restored.membership.toarray(), result.membership.toarray())
    np.testing.assert_array_equal(restored.graph.toarray(), result.graph.toarray())
    np.testing.assert_array_equal(restored.representative_coordinates, result.representative_coordinates)
    assert restored.fingerprint == result.fingerprint

    table_path = tmp_path / "table" / TABLE_PAYLOAD
    table_path.write_text(table_path.read_text(encoding="utf-8").replace("0.03", "0.93", 1), encoding="utf-8")
    with pytest.raises(ValueError, match="current-content fingerprint"):
        read_milo_result(tmp_path)


@pytest.mark.parametrize("part", ["membership", "graph", "coordinates", "observation_axis"])
def test_milo_result_rejects_tampered_payload_parts(part):
    portable = _milo_result().portable()
    if part == "membership":
        portable["membership"][0, 0] = 0
    elif part == "graph":
        portable["graph"][0, 1] = 0
    elif part == "coordinates":
        portable["coordinates"][0, 0] = 99.0
    else:
        portable["observation_names"][0] = "other"

    with pytest.raises(ValueError, match="current-content|axes"):
        validate_milo_result(portable, exact_type=False)


def test_abundance_producers_publish_typed_evidence_before_the_retained_table():
    expected = [
        ("result", MiloResultType.io_type),
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    assert [
        (output.display_name, output.io_type)
        for output in OpenBioSingleCellMiloDifferentialAbundance.define_schema().outputs
    ] == expected

    composition_expected = [
        ("result", CompositionModelResultType.io_type),
        ("table", TableResultType.io_type),
        ("summary", SummaryResultType.io_type),
        ("code", "STRING"),
    ]
    for node in (
        OpenBioSingleCellSccodaDifferentialComposition,
        OpenBioSingleCellTasccodaDifferentialComposition,
    ):
        assert [(output.display_name, output.io_type) for output in node.define_schema().outputs] == composition_expected


def _sccoda_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ["B_vs_A", "condition", "A", "B", "B", "T", 0.4, 0.4, 0.4, 0.0, 1.0, True, 0.05, 1.0, 0.0, 10.0, 12.0, 0.2, False],
            ["B_vs_A", "condition", "A", "B", "T", "T", 0.0, 0.0, 0.0, 0.0, 0.0, False, 0.05, 1.0, 0.0, 8.0, 6.0, -0.4, True],
        ],
        columns=SCCODA_RESULT_COLUMNS,
    )


def _sccoda_result():
    return build_composition_model_result(
        table=_sccoda_table(),
        method="sccoda",
        posterior={
            "intercept": np.array([[0.1, -0.1]] * 4),
            "condition_effect": np.array([[0.4, 0.0]] * 4),
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
            "reference_condition": "A",
            "comparison_condition": "B",
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
            "software_versions": {"pertpy": "1.3.0", "arviz": "1.3.0"},
        },
    )


def test_composition_model_result_retains_complete_posterior_draws_and_rejects_tampering():
    result = _sccoda_result()
    table, posterior, sample_stats, cell_types, hierarchy, model_metadata, metadata = (
        validate_composition_model_result(result)
    )

    pd.testing.assert_frame_equal(table, _sccoda_table())
    assert posterior["condition_effect"].shape == (1, 4, 2)
    assert sample_stats["potential_energy"].shape == (1, 4)
    assert cell_types == ["B", "T"]
    assert hierarchy is None
    assert model_metadata["model"]["posterior_draws"] == 4
    assert model_metadata["annotation_status"] == "curated"
    assert metadata["artifact_fingerprint_sha256"] == result.fingerprint

    portable = result.portable()
    portable["posterior"]["condition_effect"][0, 0, 0] = 99.0
    with pytest.raises(ValueError, match="posterior.*fingerprint"):
        validate_composition_model_result(portable, exact_type=False)

    axis_tampered = result.portable()
    axis_tampered["cell_types"].reverse()
    with pytest.raises(ValueError, match="axes|cell-type axis"):
        validate_composition_model_result(axis_tampered, exact_type=False)


def test_composition_model_codec_uses_portable_arviz_store_and_rejects_tampering(tmp_path):
    result = _sccoda_result()

    descriptors = write_composition_model_result(tmp_path, result)
    restored = read_composition_model_result(tmp_path)

    assert COMPOSITION_MODEL_CODEC == "composition-inferencedata-zarr-v1"
    assert any(item["path"].startswith("inference_data/") for item in descriptors)
    pd.testing.assert_frame_equal(restored.table, result.table)
    assert restored.fingerprint == result.fingerprint
    inference_data = restored.inference_data()
    assert set(inference_data.children) == {"posterior", "sample_stats"}
    assert inference_data["posterior"].to_dataset().sizes == {"chain": 1, "draw": 4, "cell_type": 2}

    table_path = tmp_path / "table" / TABLE_PAYLOAD
    table_path.write_text(table_path.read_text(encoding="utf-8").replace("0.4", "0.9", 1), encoding="utf-8")
    with pytest.raises(ValueError, match="current-content fingerprint"):
        read_composition_model_result(tmp_path)


def test_composition_model_codec_rejects_tampered_posterior_draw(tmp_path):
    import zarr

    write_composition_model_result(tmp_path, _sccoda_result())
    posterior = zarr.open_group(tmp_path / "inference_data", mode="r+")["posterior"]["condition_effect"]
    posterior[0, 0, 0] = 9.0

    with pytest.raises(ValueError, match="current-content fingerprint"):
        read_composition_model_result(tmp_path)


@pytest.mark.parametrize(
    ("build", "write", "kind", "codec"),
    [
        (_milo_result, write_milo_result, MILO_RESULT_KIND, MILO_RESULT_CODEC),
        (
            _sccoda_result,
            write_composition_model_result,
            COMPOSITION_MODEL_KIND,
            COMPOSITION_MODEL_CODEC,
        ),
    ],
)
def test_abundance_diagnostic_artifacts_are_portable_under_the_existing_persistence_policy(
    tmp_path, build, write, kind, codec
):
    identifiers = iter(("session", "run"))
    runtime = ArtifactRuntime(tmp_path / codec, uuid_factory=lambda: next(identifiers))
    lease = runtime.begin_run()
    payload = lease.staging_path / "result"
    payload.mkdir()
    write(payload, build())
    ticket = lease.publish({"result": {"kind": kind, "codec": codec, "payload": "result"}})["result"]
    output_root = tmp_path / f"output-{kind}"
    output_root.mkdir()

    persisted = persist_artifact(runtime, ticket, output_root, "diagnostic")

    assert (persisted.path / "manifest.json").is_file()
    assert runtime.resolve(ticket).is_dir()
