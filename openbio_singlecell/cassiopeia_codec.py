from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .cassiopeia_tree import (
    _ARTIFACT_TOKEN,
    MISSING_STATE,
    CassiopeiaCharacters,
    CassiopeiaTree,
    _normalize_priors,
    _normalize_state_maps,
    _require_cassiopeia,
    _validate_character_frame,
    _validate_prior_closure,
)

CHARACTERS_CODEC = "openbio-cassiopeia-characters"
TREE_CODEC = "openbio-cassiopeia-tree"
CODEC_VERSION = 1


def _strict_value(value: Any) -> Any:
    if value is pd.NA:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Cassiopeia JSON values must be finite.")
        return value
    if isinstance(value, Mapping):
        result: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise TypeError("Cassiopeia JSON object keys must be strings.")
            result[key] = _strict_value(item)
        return result
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_strict_value(item) for item in value]
    raise TypeError(f"Cassiopeia metadata contains a non-JSON value: {type(value).__name__}.")


def _json_bytes(value: Any) -> bytes:
    return json.dumps(
        _strict_value(value),
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _json_hash(value: Any) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _write_json(path: Path, value: Any) -> str:
    payload = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(payload)
    return hashlib.sha256(payload).hexdigest()


def _reject_constant(value: str) -> None:
    raise ValueError(f"Cassiopeia JSON contains non-finite value {value}.")


def _finite_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError("Cassiopeia JSON contains a non-finite number.")
    return parsed


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Cassiopeia JSON contains duplicate key {key!r}.")
        result[key] = value
    return result


def _read_json(path: Path) -> Any:
    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=_reject_constant,
            parse_float=_finite_float,
            object_pairs_hook=_unique_object,
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"Cassiopeia JSON is unreadable: {path.name}.") from error


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_root(directory: str | Path) -> Path:
    root = Path(directory).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("Cassiopeia artifact path must be a directory.")
    return root


def _output_root(directory: str | Path) -> Path:
    root = Path(directory)
    try:
        root.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        if not root.is_dir() or any(root.iterdir()):
            raise
    return root


def _member(root: Path, relative: Any) -> Path:
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("Cassiopeia artifact member path is invalid or escapes artifact directory.")
    candidate = (root / relative).resolve(strict=False)
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise ValueError("Cassiopeia artifact member path escapes artifact directory.") from error
    if not candidate.is_file():
        raise ValueError(f"Cassiopeia artifact member is missing: {relative}.")
    return candidate


def _verified_member(root: Path, relative: Any, expected_sha256: Any) -> Path:
    if not isinstance(expected_sha256, str) or len(expected_sha256) != 64:
        raise ValueError("Cassiopeia artifact member hash is invalid.")
    path = _member(root, relative)
    if _sha256_file(path) != expected_sha256:
        raise ValueError(f"Cassiopeia artifact member failed its content fingerprint: {path.name}.")
    return path


def _exact_mapping(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise ValueError(f"Cassiopeia {label} schema is invalid.")
    return value


def _exact_strings(values: Any, label: str) -> list[str]:
    if not isinstance(values, list):
        raise ValueError(f"Cassiopeia {label} must be a JSON array.")
    result = list(values)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in result):
        raise ValueError(f"Cassiopeia {label} must contain exact non-empty strings.")
    if len(set(result)) != len(result):
        raise ValueError(f"Cassiopeia {label} must be unique.")
    return result


def _matrix(frame: Any, *, _owned: bool = False) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Cassiopeia character matrices must be pandas DataFrames.")
    if frame.to_numpy(copy=False).dtype.hasobject:
        raise TypeError("Cassiopeia character matrices cannot use object dtype.")
    return _validate_character_frame(frame, _owned=_owned)


def _save_matrix(path: Path, frame: pd.DataFrame) -> str:
    with path.open("xb") as stream:
        np.save(stream, frame.to_numpy(dtype=np.int64, copy=False), allow_pickle=False)
    return _sha256_file(path)


def _load_matrix(path: Path, rows: list[str], columns: list[str]) -> pd.DataFrame:
    values = np.load(path, allow_pickle=False)
    if values.dtype.hasobject:
        raise TypeError("Cassiopeia character matrices cannot use object dtype.")
    if values.ndim != 2 or values.shape != (len(rows), len(columns)):
        raise ValueError("Cassiopeia matrix shape does not match its declared axes.")
    return _matrix(pd.DataFrame(values, index=rows, columns=columns), _owned=True)


def _pairs(priors: Mapping[int, Mapping[int, float]]) -> list[list[Any]]:
    return [
        [character, [[state, probability] for state, probability in states.items()]]
        for character, states in priors.items()
    ]


def _checked_priors(value: Mapping[Any, Mapping[Any, Any]]) -> dict[int, dict[int, float]]:
    if not isinstance(value, Mapping):
        raise TypeError("Cassiopeia priors must be a mapping.")
    for character, states in value.items():
        if isinstance(character, bool) or not isinstance(character, (int, np.integer)):
            raise ValueError("Cassiopeia prior character keys must be integers.")
        if not isinstance(states, Mapping):
            raise ValueError("Cassiopeia prior states must be mappings.")
        for state, probability in states.items():
            if isinstance(state, bool) or not isinstance(state, (int, np.integer)):
                raise ValueError("Cassiopeia prior state keys must be integers.")
            if isinstance(probability, (bool, np.bool_)):
                raise ValueError("Cassiopeia prior probabilities must be finite values in (0, 1].")
    return _normalize_priors(value)


def _checked_state_maps(value: Mapping[Any, Mapping[Any, Any]]) -> dict[int, dict[int, str]]:
    if not isinstance(value, Mapping):
        raise TypeError("Cassiopeia state maps must be a mapping.")
    for character, states in value.items():
        if isinstance(character, bool) or not isinstance(character, (int, np.integer)):
            raise ValueError("Cassiopeia state-map character keys must be integers.")
        if not isinstance(states, Mapping):
            raise ValueError("Cassiopeia state-map states must be mappings.")
        if any(isinstance(state, bool) or not isinstance(state, (int, np.integer)) for state in states):
            raise ValueError("Cassiopeia state-map state keys must be integers.")
    return _normalize_state_maps(value)


def _map_pairs(state_maps: Mapping[int, Mapping[int, str]]) -> list[list[Any]]:
    return [
        [character, [[state, allele] for state, allele in states.items()]] for character, states in state_maps.items()
    ]


def _nested_pairs(value: Any, label: str) -> dict[int, dict[int, Any]]:
    if not isinstance(value, list):
        raise ValueError(f"Cassiopeia {label} must use pair arrays.")
    result: dict[int, dict[int, Any]] = {}
    for character_pair in value:
        if not isinstance(character_pair, list) or len(character_pair) != 2:
            raise ValueError(f"Cassiopeia {label} must use pair arrays.")
        character, state_pairs = character_pair
        if isinstance(character, bool) or not isinstance(character, int) or character in result:
            raise ValueError(f"Cassiopeia {label} character keys are invalid.")
        if not isinstance(state_pairs, list):
            raise ValueError(f"Cassiopeia {label} state values must use pair arrays.")
        states: dict[int, Any] = {}
        for state_pair in state_pairs:
            if not isinstance(state_pair, list) or len(state_pair) != 2:
                raise ValueError(f"Cassiopeia {label} state values must use pair arrays.")
            state, item = state_pair
            if isinstance(state, bool) or not isinstance(state, int) or state in states:
                raise ValueError(f"Cassiopeia {label} state keys are invalid.")
            states[state] = item
        result[character] = states
    return result


def _character_payload(
    matrices: Mapping[str, Any],
    priors: Mapping[str, Mapping[int, Mapping[int, float]]],
    state_maps: Mapping[str, Mapping[int, Mapping[int, str]]],
    *,
    _owned: bool = False,
) -> list[tuple[str, pd.DataFrame, dict[int, dict[int, float]], dict[int, dict[int, str]]]]:
    if not isinstance(matrices, Mapping) or not isinstance(priors, Mapping) or not isinstance(state_maps, Mapping):
        raise TypeError("Cassiopeia character payloads must be mappings.")
    lineage_ids = list(matrices)
    if set(lineage_ids) != set(priors) or set(lineage_ids) != set(state_maps):
        raise ValueError("Cassiopeia character lineage payloads are inconsistent.")
    result = []
    for lineage_id in lineage_ids:
        if not isinstance(lineage_id, str) or not lineage_id or lineage_id != lineage_id.strip():
            raise ValueError("Cassiopeia lineage IDs must be exact non-empty strings.")
        frame = _matrix(matrices[lineage_id], _owned=_owned)
        lineage_priors = _checked_priors(priors[lineage_id])
        lineage_maps = _checked_state_maps(state_maps[lineage_id])
        _validate_prior_closure(frame, lineage_priors, lineage_maps)
        result.append((lineage_id, frame, lineage_priors, lineage_maps))
    return result


def _manifest(codec: str, body: dict[str, Any]) -> dict[str, Any]:
    manifest = {"codec": codec, "version": CODEC_VERSION, **body}
    manifest["content_fingerprint"] = _json_hash(manifest)
    return manifest


def _validated_manifest(path: Path, codec: str, body_keys: set[str]) -> dict[str, Any]:
    manifest = _exact_mapping(
        _read_json(path),
        {"codec", "version", "content_fingerprint", *body_keys},
        "manifest",
    )
    if manifest["codec"] != codec or manifest["version"] != CODEC_VERSION:
        raise ValueError("Cassiopeia codec identity or version is unsupported.")
    fingerprint = manifest.pop("content_fingerprint")
    if not isinstance(fingerprint, str) or _json_hash(manifest) != fingerprint:
        raise ValueError("Cassiopeia artifact content fingerprint is invalid.")
    manifest["content_fingerprint"] = fingerprint
    return manifest


def encode_characters(
    directory: str | Path,
    *,
    matrices: Mapping[str, Any],
    priors: Mapping[str, Mapping[int, Mapping[int, float]]],
    state_maps: Mapping[str, Mapping[int, Mapping[int, str]]],
    metadata: Mapping[str, Any],
    _worker_owned: bool = False,
) -> dict[str, Any]:
    normalized = _character_payload(matrices, priors, state_maps, _owned=_worker_owned)
    metadata_value = _strict_value(metadata)
    if not isinstance(metadata_value, dict):
        raise TypeError("Cassiopeia character metadata must be a mapping.")
    root = _output_root(directory)
    metadata_hash = _write_json(root / "metadata.json", metadata_value)
    lineage_entries: list[dict[str, Any]] = []
    for index, (lineage_id, frame, lineage_priors, lineage_maps) in enumerate(normalized):
        stem = f"lineage_{index:04d}"
        matrix_name = f"{stem}.npy"
        details_name = f"{stem}.json"
        matrix_hash = _save_matrix(root / matrix_name, frame)
        details_hash = _write_json(
            root / details_name,
            {
                "lineage_id": lineage_id,
                "row_axis": list(frame.index),
                "column_axis": list(frame.columns),
                "priors": _pairs(lineage_priors),
                "state_maps": _map_pairs(lineage_maps),
            },
        )
        lineage_entries.append(
            {
                "lineage_id": lineage_id,
                "matrix": matrix_name,
                "matrix_sha256": matrix_hash,
                "details": details_name,
                "details_sha256": details_hash,
            }
        )
    manifest = _manifest(
        CHARACTERS_CODEC,
        {
            "lineages": lineage_entries,
            "metadata": "metadata.json",
            "metadata_sha256": metadata_hash,
        },
    )
    _write_json(root / "characters.json", manifest)
    return manifest


def decode_characters(directory: str | Path) -> dict[str, Any]:
    root = _artifact_root(directory)
    manifest = _validated_manifest(
        _member(root, "characters.json"), CHARACTERS_CODEC, {"lineages", "metadata", "metadata_sha256"}
    )
    if manifest["metadata"] != "metadata.json":
        raise ValueError("Cassiopeia character metadata must use its fixed filename.")
    metadata = _read_json(_verified_member(root, manifest["metadata"], manifest["metadata_sha256"]))
    if not isinstance(metadata, dict):
        raise ValueError("Cassiopeia character metadata must be a JSON object.")
    if not isinstance(manifest["lineages"], list):
        raise ValueError("Cassiopeia character manifest lineages must be an array.")
    matrices: dict[str, pd.DataFrame] = {}
    priors: dict[str, dict[int, dict[int, float]]] = {}
    state_maps: dict[str, dict[int, dict[int, str]]] = {}
    for index, raw_entry in enumerate(manifest["lineages"]):
        entry = _exact_mapping(
            raw_entry,
            {"lineage_id", "matrix", "matrix_sha256", "details", "details_sha256"},
            "lineage manifest",
        )
        matrix_path = _verified_member(root, entry["matrix"], entry["matrix_sha256"])
        details_path = _verified_member(root, entry["details"], entry["details_sha256"])
        if matrix_path.name != f"lineage_{index:04d}.npy" or details_path.name != f"lineage_{index:04d}.json":
            raise ValueError("Cassiopeia lineage payloads must use fixed numbered filenames.")
        details = _exact_mapping(
            _read_json(details_path),
            {"lineage_id", "row_axis", "column_axis", "priors", "state_maps"},
            "lineage details",
        )
        lineage_id = entry["lineage_id"]
        if (
            not isinstance(lineage_id, str)
            or not lineage_id
            or lineage_id != lineage_id.strip()
            or details["lineage_id"] != lineage_id
            or lineage_id in matrices
        ):
            raise ValueError("Cassiopeia lineage identity is invalid or duplicated.")
        rows = _exact_strings(details["row_axis"], "row axis")
        columns = _exact_strings(details["column_axis"], "column axis")
        frame = _load_matrix(matrix_path, rows, columns)
        lineage_priors = _checked_priors(_nested_pairs(details["priors"], "priors"))
        lineage_maps = _checked_state_maps(_nested_pairs(details["state_maps"], "state maps"))
        _validate_prior_closure(frame, lineage_priors, lineage_maps)
        matrices[lineage_id] = frame
        priors[lineage_id] = lineage_priors
        state_maps[lineage_id] = lineage_maps
    return {
        "lineage_ids": tuple(matrices),
        "matrices": matrices,
        "priors": priors,
        "state_maps": state_maps,
        "metadata": metadata,
        "content_fingerprint": manifest["content_fingerprint"],
    }


def _topology(
    root: str, nodes: Sequence[str], edges: Sequence[Sequence[Any]]
) -> tuple[list[str], list[tuple[str, str, float]]]:
    if isinstance(nodes, (str, bytes, bytearray)):
        raise ValueError("Cassiopeia topology nodes must be a sequence of identities.")
    node_axis = _exact_strings(list(nodes), "topology nodes")
    if not isinstance(root, str) or root not in node_axis:
        raise ValueError("Cassiopeia topology root is invalid.")
    normalized_edges: list[tuple[str, str, float]] = []
    parents = {node: [] for node in node_axis}
    children = {node: [] for node in node_axis}
    seen: set[tuple[str, str]] = set()
    for edge in edges:
        if not isinstance(edge, Sequence) or isinstance(edge, (str, bytes, bytearray)) or len(edge) != 3:
            raise ValueError("Cassiopeia topology edges must be parent, child, branch-length triples.")
        parent, child, raw_length = edge
        if parent not in parents or child not in parents or parent == child or (parent, child) in seen:
            raise ValueError("Cassiopeia topology contains an invalid or duplicate edge.")
        if isinstance(raw_length, bool):
            raise ValueError("Cassiopeia branch lengths must be finite and non-negative.")
        try:
            length = float(raw_length)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError("Cassiopeia branch lengths must be finite and non-negative.") from error
        if not math.isfinite(length) or length < 0.0:
            raise ValueError("Cassiopeia branch lengths must be finite and non-negative.")
        seen.add((parent, child))
        parents[child].append(parent)
        children[parent].append(child)
        normalized_edges.append((parent, child, length))
    if parents[root] or any(len(parents[node]) != 1 for node in node_axis if node != root):
        raise ValueError("Cassiopeia topology must give every non-root node exactly one parent.")
    if len(normalized_edges) != len(node_axis) - 1:
        raise ValueError("Cassiopeia topology is not a directed tree.")
    visited: set[str] = set()

    def visit(node: str) -> None:
        if node in visited:
            raise ValueError("Cassiopeia topology contains a cycle.")
        visited.add(node)
        for child in children[node]:
            visit(child)

    visit(root)
    if visited != set(node_axis):
        raise ValueError("Cassiopeia topology is disconnected from its root.")
    return node_axis, normalized_edges


def _tree_payload(
    character_matrix: Any,
    priors: Mapping[int, Mapping[int, float]],
    root: str,
    nodes: Sequence[str],
    edges: Sequence[Sequence[Any]],
    *,
    _owned: bool = False,
) -> tuple[pd.DataFrame, dict[int, dict[int, float]], list[str], list[tuple[str, str, float]]]:
    frame = _matrix(character_matrix, _owned=_owned)
    normalized_priors = _checked_priors(priors)
    observed_maps = {
        position: {
            int(state): f"observed_state_{int(state)}"
            for state in np.unique(frame.iloc[:, position].to_numpy())
            if int(state) > 0
        }
        for position in range(frame.shape[1])
    }
    _validate_prior_closure(frame, normalized_priors, observed_maps)
    node_axis, normalized_edges = _topology(root, nodes, edges)
    parents = {child for _, child, _ in normalized_edges}
    leaves = {node for node in node_axis if node not in {parent for parent, _, _ in normalized_edges}}
    if set(frame.index) != leaves or root in parents:
        raise ValueError("Cassiopeia character-matrix cells must exactly match topology leaves.")
    return frame, normalized_priors, node_axis, normalized_edges


def encode_tree(
    directory: str | Path,
    *,
    character_matrix: Any,
    priors: Mapping[int, Mapping[int, float]],
    root: str,
    nodes: Sequence[str],
    edges: Sequence[Sequence[Any]],
    metadata: Mapping[str, Any],
    _worker_owned: bool = False,
) -> dict[str, Any]:
    frame, normalized_priors, node_axis, normalized_edges = _tree_payload(
        character_matrix, priors, root, nodes, edges, _owned=_worker_owned
    )
    metadata_value = _strict_value(metadata)
    if not isinstance(metadata_value, dict):
        raise TypeError("Cassiopeia tree metadata must be a mapping.")
    target = _output_root(directory)
    hashes = {
        "matrix_sha256": _save_matrix(target / "character_matrix.npy", frame),
        "axes_sha256": _write_json(
            target / "axes.json", {"row_axis": list(frame.index), "column_axis": list(frame.columns)}
        ),
        "topology_sha256": _write_json(
            target / "topology.json",
            {"root": root, "nodes": node_axis, "edges": [list(edge) for edge in normalized_edges]},
        ),
        "priors_sha256": _write_json(target / "priors.json", _pairs(normalized_priors)),
        "metadata_sha256": _write_json(target / "metadata.json", metadata_value),
    }
    manifest = _manifest(
        TREE_CODEC,
        {
            "matrix": "character_matrix.npy",
            "axes": "axes.json",
            "topology": "topology.json",
            "priors": "priors.json",
            "metadata": "metadata.json",
            **hashes,
        },
    )
    _write_json(target / "tree.json", manifest)
    return manifest


def decode_tree(directory: str | Path) -> dict[str, Any]:
    root_path = _artifact_root(directory)
    manifest = _validated_manifest(
        _member(root_path, "tree.json"),
        TREE_CODEC,
        {
            "matrix",
            "matrix_sha256",
            "axes",
            "axes_sha256",
            "topology",
            "topology_sha256",
            "priors",
            "priors_sha256",
            "metadata",
            "metadata_sha256",
        },
    )
    fixed = {
        "matrix": "character_matrix.npy",
        "axes": "axes.json",
        "topology": "topology.json",
        "priors": "priors.json",
        "metadata": "metadata.json",
    }
    paths: dict[str, Path] = {}
    for key, filename in fixed.items():
        paths[key] = _verified_member(root_path, manifest[key], manifest[f"{key}_sha256"])
        if paths[key].name != filename:
            raise ValueError("Cassiopeia tree payload must use fixed filenames.")
    axes = _exact_mapping(_read_json(paths["axes"]), {"row_axis", "column_axis"}, "tree axes")
    rows = _exact_strings(axes["row_axis"], "row axis")
    columns = _exact_strings(axes["column_axis"], "column axis")
    frame = _load_matrix(paths["matrix"], rows, columns)
    topology = _exact_mapping(_read_json(paths["topology"]), {"root", "nodes", "edges"}, "topology")
    node_axis, normalized_edges = _topology(topology["root"], topology["nodes"], topology["edges"])
    normalized_priors = _checked_priors(_nested_pairs(_read_json(paths["priors"]), "priors"))
    frame, normalized_priors, node_axis, normalized_edges = _tree_payload(
        frame,
        normalized_priors,
        topology["root"],
        node_axis,
        normalized_edges,
        _owned=True,
    )
    metadata = _read_json(paths["metadata"])
    if not isinstance(metadata, dict):
        raise ValueError("Cassiopeia tree metadata must be a JSON object.")
    return {
        "character_matrix": frame,
        "priors": normalized_priors,
        "root": topology["root"],
        "nodes": tuple(node_axis),
        "edges": tuple(normalized_edges),
        "metadata": metadata,
        "content_fingerprint": manifest["content_fingerprint"],
    }


def write_characters(directory: str | Path, artifact: CassiopeiaCharacters) -> dict[str, Any]:
    if type(artifact) is not CassiopeiaCharacters:
        raise TypeError("Expected a CassiopeiaCharacters value.")
    payloads = {
        lineage_id: artifact.copy_lineage(lineage_id, _owned=True) for lineage_id in artifact.lineage_ids
    }
    return encode_characters(
        directory,
        matrices={lineage_id: payload[0] for lineage_id, payload in payloads.items()},
        priors={lineage_id: payload[1] for lineage_id, payload in payloads.items()},
        state_maps={lineage_id: payload[2] for lineage_id, payload in payloads.items()},
        metadata=artifact.metadata,
        _worker_owned=True,
    )


def read_characters(directory: str | Path) -> CassiopeiaCharacters:
    payload = decode_characters(directory)
    return CassiopeiaCharacters(
        payload["matrices"],
        payload["priors"],
        payload["state_maps"],
        payload["metadata"],
        _token=_ARTIFACT_TOKEN,
        _owned=True,
    )


def write_tree(directory: str | Path, artifact: CassiopeiaTree) -> dict[str, Any]:
    if type(artifact) is not CassiopeiaTree:
        raise TypeError("Expected a CassiopeiaTree value.")
    backend = artifact.copy_tree(_owned=True)
    graph = backend.get_tree_topology()
    return encode_tree(
        directory,
        character_matrix=backend.character_matrix,
        priors=backend.priors,
        root=backend.root,
        nodes=list(graph.nodes),
        edges=[
            (parent, child, float(attributes.get("length", 1.0)))
            for parent, child, attributes in graph.edges(data=True)
        ],
        metadata=artifact.metadata,
        _worker_owned=True,
    )


def read_tree(directory: str | Path) -> CassiopeiaTree:
    payload = decode_tree(directory)
    cassiopeia = _require_cassiopeia()
    try:
        import networkx as nx
    except ImportError as error:
        raise RuntimeError("Cassiopeia tree artifacts require networkx.") from error
    graph = nx.DiGraph()
    graph.add_nodes_from(payload["nodes"])
    graph.add_edges_from((parent, child, {"length": length}) for parent, child, length in payload["edges"])
    backend = cassiopeia.data.CassiopeiaTree(
        character_matrix=payload["character_matrix"],
        missing_state_indicator=MISSING_STATE,
        priors=payload["priors"],
        tree=graph,
    )
    return CassiopeiaTree(backend, payload["metadata"], _token=_ARTIFACT_TOKEN, _owned=True)


__all__ = [
    "CHARACTERS_CODEC",
    "TREE_CODEC",
    "decode_characters",
    "decode_tree",
    "encode_characters",
    "encode_tree",
    "read_characters",
    "read_tree",
    "write_characters",
    "write_tree",
]
