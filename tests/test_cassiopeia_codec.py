from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import openbio_singlecell.cassiopeia_codec as codec
from openbio_singlecell.cassiopeia_codec import (
    decode_characters,
    decode_tree,
    encode_characters,
    encode_tree,
    read_characters,
    read_tree,
    write_characters,
    write_tree,
)
from openbio_singlecell.cassiopeia_tree import _ARTIFACT_TOKEN, CassiopeiaCharacters, CassiopeiaTree


def _character_inputs():
    matrices = {
        "T/unsafe": pd.DataFrame([[2, 0], [1, -1]], index=["cell-b", "cell-a"], columns=["r2", "r1"]),
        "../other": pd.DataFrame([[1], [1]], index=["x2", "x1"], columns=["cut"]),
    }
    priors = {
        "T/unsafe": {0: {2: 0.25, 1: 0.75}, 1: {}},
        "../other": {0: {1: 1.0}},
    }
    state_maps = {
        "T/unsafe": {0: {2: "DEL-B", 1: "DEL-A"}, 1: {}},
        "../other": {0: {1: "INS"}},
    }
    metadata = {"fingerprint": "upstream-fingerprint", "nested": {"ok": [True, None, 3]}}
    return matrices, priors, state_maps, metadata


def test_characters_round_trip_uses_numbered_files_pair_arrays_and_preserves_axes(tmp_path):
    matrices, priors, state_maps, metadata = _character_inputs()
    destination = tmp_path / "characters"

    manifest = encode_characters(
        destination,
        matrices=matrices,
        priors=priors,
        state_maps=state_maps,
        metadata=metadata,
    )

    assert sorted(path.name for path in destination.iterdir()) == [
        "characters.json",
        "lineage_0000.json",
        "lineage_0000.npy",
        "lineage_0001.json",
        "lineage_0001.npy",
        "metadata.json",
    ]
    details = json.loads((destination / "lineage_0000.json").read_text(encoding="utf-8"))
    assert details["lineage_id"] == "T/unsafe"
    assert details["priors"] == [[0, [[1, 0.75], [2, 0.25]]], [1, []]]
    assert details["state_maps"] == [[0, [[1, "DEL-A"], [2, "DEL-B"]]], [1, []]]
    assert np.load(destination / "lineage_0000.npy", allow_pickle=False).dtype == np.dtype("int64")

    decoded = decode_characters(destination)

    assert decoded["lineage_ids"] == ("T/unsafe", "../other")
    for lineage_id, expected in matrices.items():
        assert_frame_equal(decoded["matrices"][lineage_id], expected.astype("int64"))
    assert decoded["priors"] == priors
    assert decoded["state_maps"] == state_maps
    assert decoded["metadata"] == metadata
    assert decoded["content_fingerprint"] == manifest["content_fingerprint"]
    reencoded = encode_characters(
        tmp_path / "characters-reencoded",
        matrices=decoded["matrices"],
        priors=decoded["priors"],
        state_maps=decoded["state_maps"],
        metadata=decoded["metadata"],
    )
    assert reencoded["content_fingerprint"] == manifest["content_fingerprint"]


def test_tree_round_trip_preserves_declared_node_edge_and_axis_order(tmp_path):
    matrix = pd.DataFrame([[1, 0], [2, -1]], index=["leaf-b", "leaf-a"], columns=["r2", "r1"])
    destination = tmp_path / "tree"

    manifest = encode_tree(
        destination,
        character_matrix=matrix,
        priors={0: {1: 0.4, 2: 0.6}, 1: {}},
        root="root",
        nodes=["root", "leaf-a", "leaf-b"],
        edges=[("root", "leaf-b", 2.0), ("root", "leaf-a", 1.0)],
        metadata={"fingerprint": "tree-fingerprint", "lineage_id": "T1"},
    )

    decoded = decode_tree(destination)

    assert decoded["nodes"] == ("root", "leaf-a", "leaf-b")
    assert decoded["edges"] == (("root", "leaf-b", 2.0), ("root", "leaf-a", 1.0))
    assert_frame_equal(decoded["character_matrix"], matrix.astype("int64"))
    assert decoded["priors"] == {0: {1: 0.4, 2: 0.6}, 1: {}}
    assert decoded["metadata"]["fingerprint"] == "tree-fingerprint"
    assert decoded["content_fingerprint"] == manifest["content_fingerprint"]
    reencoded = encode_tree(
        tmp_path / "tree-reencoded",
        character_matrix=decoded["character_matrix"],
        priors=decoded["priors"],
        root=decoded["root"],
        nodes=decoded["nodes"],
        edges=decoded["edges"],
        metadata=decoded["metadata"],
    )
    assert reencoded["content_fingerprint"] == manifest["content_fingerprint"]


@pytest.mark.parametrize("encoder", ["characters", "tree"])
def test_codecs_reject_object_character_matrices(tmp_path, encoder):
    matrix = pd.DataFrame(np.array([[object()]], dtype=object), index=["leaf"], columns=["r1"])

    with pytest.raises(TypeError, match="object dtype"):
        if encoder == "characters":
            encode_characters(
                tmp_path / "characters",
                matrices={"T1": matrix},
                priors={"T1": {0: {}}},
                state_maps={"T1": {0: {}}},
                metadata={},
            )
        else:
            encode_tree(
                tmp_path / "tree",
                character_matrix=matrix,
                priors={0: {}},
                root="leaf",
                nodes=["leaf"],
                edges=[],
                metadata={},
            )


def test_codecs_reject_nonfinite_priors_and_branch_lengths(tmp_path):
    matrix = pd.DataFrame([[1]], index=["leaf"], columns=["r1"])

    with pytest.raises(ValueError, match="prior.*finite"):
        encode_characters(
            tmp_path / "characters",
            matrices={"T1": matrix},
            priors={"T1": {0: {1: float("inf")}}},
            state_maps={"T1": {0: {1: "DEL"}}},
            metadata={},
        )
    with pytest.raises(ValueError, match="branch.*finite"):
        encode_tree(
            tmp_path / "tree",
            character_matrix=matrix,
            priors={0: {1: 1.0}},
            root="root",
            nodes=["root", "leaf"],
            edges=[("root", "leaf", float("nan"))],
            metadata={},
        )


def _resign_manifest(path: Path) -> None:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest.pop("content_fingerprint")
    encoded = json.dumps(manifest, allow_nan=False, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode(
        "utf-8"
    )
    manifest["content_fingerprint"] = hashlib.sha256(encoded).hexdigest()
    path.write_text(
        json.dumps(manifest, allow_nan=False, ensure_ascii=False, separators=(",", ":"), sort_keys=True),
        encoding="utf-8",
    )


def test_decode_rejects_manifest_path_escape_before_reading_outside_file(tmp_path):
    matrices, priors, state_maps, metadata = _character_inputs()
    destination = tmp_path / "characters"
    encode_characters(
        destination,
        matrices=matrices,
        priors=priors,
        state_maps=state_maps,
        metadata=metadata,
    )
    manifest_path = destination / "characters.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["lineages"][0]["matrix"] = "../outside.npy"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _resign_manifest(manifest_path)

    with pytest.raises(ValueError, match="escapes artifact directory"):
        decode_characters(destination)


def test_decode_rejects_object_npy_without_enabling_pickle(tmp_path):
    matrices, priors, state_maps, metadata = _character_inputs()
    destination = tmp_path / "characters"
    encode_characters(
        destination,
        matrices=matrices,
        priors=priors,
        state_maps=state_maps,
        metadata=metadata,
    )
    matrix_path = destination / "lineage_0000.npy"
    np.save(matrix_path, np.array([[object()]], dtype=object), allow_pickle=True)
    manifest_path = destination / "characters.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["lineages"][0]["matrix_sha256"] = hashlib.sha256(matrix_path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _resign_manifest(manifest_path)

    with pytest.raises(ValueError, match="Object arrays cannot be loaded"):
        decode_characters(destination)


class _BackendTree:
    def __init__(
        self,
        character_matrix=None,
        missing_state_indicator=-1,
        priors=None,
        tree=None,
        **_kwargs,
    ):
        self.character_matrix = character_matrix
        self.missing_state_indicator = missing_state_indicator
        self.priors = priors or {}
        self._graph = tree

    def copy(self):
        return copy.deepcopy(self)

    @property
    def root(self):
        return next(node for node in self._graph if self._graph.in_degree(node) == 0)

    def get_tree_topology(self):
        return self._graph.copy()


def test_worker_character_wrapper_round_trip(tmp_path):
    matrices, priors, state_maps, metadata = _character_inputs()
    artifact = CassiopeiaCharacters(matrices, priors, state_maps, metadata, _token=_ARTIFACT_TOKEN)

    write_characters(tmp_path / "characters", artifact)
    restored = read_characters(tmp_path / "characters")

    assert restored.fingerprint == artifact.fingerprint
    assert restored.lineage_ids == artifact.lineage_ids


def test_worker_character_reader_adopts_one_loaded_matrix_until_the_backend_needs_it(tmp_path, monkeypatch):
    matrix = pd.DataFrame([[1], [2]], index=["a", "b"], columns=["r1"])
    encode_characters(
        tmp_path / "characters",
        matrices={"T1": matrix},
        priors={"T1": {0: {1: 0.4, 2: 0.6}}},
        state_maps={"T1": {0: {1: "A", 2: "B"}}},
        metadata={},
    )
    loaded = []
    original_load = codec._load_matrix

    def capture_load(*args, **kwargs):
        result = original_load(*args, **kwargs)
        loaded.append(result)
        return result

    monkeypatch.setattr(codec, "_load_matrix", capture_load)

    restored = read_characters(tmp_path / "characters")
    owned_matrix, _, _ = restored.copy_lineage("T1", _owned=True)
    defensive_matrix, _, _ = restored.copy_lineage("T1")

    assert owned_matrix is loaded[0]
    assert defensive_matrix is not owned_matrix


def test_worker_tree_wrapper_round_trip(tmp_path, monkeypatch):
    matrix = pd.DataFrame([[1], [2]], index=["a", "b"], columns=["r1"])
    graph = nx.DiGraph()
    graph.add_edge("root", "a", length=1.0)
    graph.add_edge("root", "b", length=2.0)
    artifact = CassiopeiaTree(
        _BackendTree(character_matrix=matrix, priors={0: {1: 0.4, 2: 0.6}}, tree=graph),
        {"lineage_id": "T1"},
        _token=_ARTIFACT_TOKEN,
    )
    monkeypatch.setattr(
        "openbio_singlecell.cassiopeia_codec._require_cassiopeia",
        lambda: SimpleNamespace(data=SimpleNamespace(CassiopeiaTree=_BackendTree)),
    )

    write_tree(tmp_path / "tree", artifact)
    restored = read_tree(tmp_path / "tree")

    assert restored.fingerprint == artifact.fingerprint
    assert restored.topology_fingerprint == artifact.topology_fingerprint


def test_worker_tree_reader_adopts_the_decoded_matrix_and_backend_once(tmp_path, monkeypatch):
    matrix = pd.DataFrame([[1], [2]], index=["a", "b"], columns=["r1"])
    graph = nx.DiGraph()
    graph.add_edge("root", "a", length=1.0)
    graph.add_edge("root", "b", length=2.0)
    encode_tree(
        tmp_path / "tree",
        character_matrix=matrix,
        priors={0: {1: 0.4, 2: 0.6}},
        root="root",
        nodes=["root", "a", "b"],
        edges=[("root", "a", 1.0), ("root", "b", 2.0)],
        metadata={"lineage_id": "T1"},
    )
    loaded = []
    created = []
    copied = []
    original_load = codec._load_matrix

    def capture_load(*args, **kwargs):
        result = original_load(*args, **kwargs)
        loaded.append(result)
        return result

    class TrackingTree(_BackendTree):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            created.append(self)

        def copy(self):
            copied.append(self)
            return super().copy()

    monkeypatch.setattr(codec, "_load_matrix", capture_load)
    monkeypatch.setattr(
        "openbio_singlecell.cassiopeia_codec._require_cassiopeia",
        lambda: SimpleNamespace(data=SimpleNamespace(CassiopeiaTree=TrackingTree)),
    )

    restored = read_tree(tmp_path / "tree")
    owned_tree = restored.copy_tree(_owned=True)

    assert owned_tree is created[0]
    assert owned_tree.character_matrix is loaded[0]
    assert copied == []
    assert restored.copy_tree() is not owned_tree
    assert copied == [owned_tree]


def test_worker_codec_writers_do_not_request_defensive_owned_matrix_copies(tmp_path, monkeypatch):
    matrix = pd.DataFrame([[1], [2]], index=["a", "b"], columns=["r1"], dtype="int64")
    characters = CassiopeiaCharacters(
        {"T1": matrix},
        {"T1": {0: {1: 0.4, 2: 0.6}}},
        {"T1": {0: {1: "A", 2: "B"}}},
        {},
        _token=_ARTIFACT_TOKEN,
        _owned=True,
    )
    graph = nx.DiGraph()
    graph.add_edge("root", "a", length=1.0)
    graph.add_edge("root", "b", length=2.0)
    tree = CassiopeiaTree(
        _BackendTree(character_matrix=matrix, priors={0: {1: 0.4, 2: 0.6}}, tree=graph),
        {"lineage_id": "T1"},
        _token=_ARTIFACT_TOKEN,
        _owned=True,
    )
    copied = []
    original_to_numpy = pd.DataFrame.to_numpy

    def tracked_to_numpy(frame, *args, **kwargs):
        if kwargs.get("copy") is True:
            copied.append(frame)
        return original_to_numpy(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_numpy", tracked_to_numpy)

    write_characters(tmp_path / "characters-owned", characters)
    write_tree(tmp_path / "tree-owned", tree)

    assert copied == []
