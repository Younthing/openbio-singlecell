"""Audited, portable Cassiopeia lineage-analysis core.

Only standard-library modules are imported eagerly. The complete source is the
independently executable implementation returned by every lineage ``code``
output; Cassiopeia itself is pinned, imported lazily, and never downloaded.
"""

from __future__ import annotations

import contextlib
import copy
import csv
import hashlib
import importlib
import inspect
import io
import json
import math
import os
import platform
import random
import sys
import threading
import warnings
from collections import Counter
from collections.abc import Iterator, Mapping, Sequence
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any

ALLELE_TABLE_EXTENSIONS = (".csv", ".tsv")
CASSIOPEIA_DISTRIBUTION = "cassiopeia-mt"
CASSIOPEIA_VERSION = "2.1.3"
CASSIOPEIA_SOURCE_TAG = "v2.1.3-mt"
CASSIOPEIA_SOURCE_COMMIT = "09868dd04c073d81dc60611f6f91cdbab17783f1"
CASSIOPEIA_SDIST_SHA256 = "5934043af6506a6439576816be49e170bc41dbe3766eaf49d14d05c3577a3e24"
CHARACTERS_ARTIFACT_TYPE = "OPENBIO_CASSIOPEIA_CHARACTERS"
CHARACTERS_SCHEMA_VERSION = 1
TREE_ARTIFACT_TYPE = "OPENBIO_CASSIOPEIA_TREE"
TREE_SCHEMA_VERSION = 2
PRODUCER = "openbio-singlecell"
MISSING_STATE = -1
UNCUT_STATE = 0
_INTERNAL_PREFIX = "__openbio_internal__"
_GROUP_PREFIX = "__openbio_state_group__"
_CAPTURE_LIMIT = 8192
_ARTIFACT_TOKEN = object()
_CASSIOPEIA_LOCK = threading.RLock()


def _numpy_pandas() -> tuple[Any, Any]:
    try:
        return importlib.import_module("numpy"), importlib.import_module("pandas")
    except ImportError as error:
        raise RuntimeError("Cassiopeia lineage analysis requires NumPy and pandas.") from error


def _json_scalar(value: Any) -> Any:
    np, pd = _numpy_pandas()
    if value is None or value is pd.NA:
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return None
    return str(value)


def _json_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_json_value(item) for item in value]
    if hasattr(value, "tolist"):
        try:
            return _json_value(value.tolist())
        except (TypeError, ValueError):
            pass
    return _json_scalar(value)


def _canonical_json(value: Any) -> str:
    return json.dumps(
        _json_value(value), allow_nan=False, ensure_ascii=False, separators=(",", ":"), sort_keys=True
    )


def _hash_payload(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _frame_payload(frame: Any) -> dict[str, Any]:
    return {
        "index": [_json_scalar(value) for value in frame.index.tolist()],
        "columns": [_json_scalar(value) for value in frame.columns.tolist()],
        "data": [[_json_scalar(value) for value in row] for row in frame.itertuples(index=False, name=None)],
    }


def _frame_hash(frame: Any) -> str:
    return _hash_payload(_frame_payload(frame))


def _normalize_priors(priors: Mapping[Any, Mapping[Any, Any]]) -> dict[int, dict[int, float]]:
    result: dict[int, dict[int, float]] = {}
    for character, states in priors.items():
        character_id = int(character)
        if character_id < 0 or character_id in result or not isinstance(states, Mapping):
            raise ValueError("Cassiopeia priors require unique non-negative integer character keys.")
        normalized: dict[int, float] = {}
        for state, probability in states.items():
            state_id = int(state)
            probability_value = float(probability)
            if state_id <= 0 or state_id in normalized:
                raise ValueError("Cassiopeia prior state keys must be unique positive integers.")
            if not math.isfinite(probability_value) or not 0.0 < probability_value <= 1.0:
                raise ValueError("Cassiopeia prior probabilities must be finite values in (0, 1].")
            normalized[state_id] = probability_value
        result[character_id] = dict(sorted(normalized.items()))
    return dict(sorted(result.items()))


def _normalize_state_maps(state_maps: Mapping[Any, Mapping[Any, Any]]) -> dict[int, dict[int, str]]:
    result: dict[int, dict[int, str]] = {}
    for character, states in state_maps.items():
        character_id = int(character)
        if character_id < 0 or character_id in result or not isinstance(states, Mapping):
            raise ValueError("Cassiopeia state maps require unique non-negative integer character keys.")
        normalized: dict[int, str] = {}
        for state, allele in states.items():
            state_id = int(state)
            if state_id <= 0 or state_id in normalized or not isinstance(allele, str) or not allele:
                raise ValueError("Cassiopeia state maps require positive integer states and non-empty alleles.")
            normalized[state_id] = allele
        result[character_id] = dict(sorted(normalized.items()))
    return dict(sorted(result.items()))


def _priors_payload(priors: Mapping[int, Mapping[int, float]]) -> list[list[Any]]:
    return [
        [int(character), [[int(state), float(probability)] for state, probability in sorted(states.items())]]
        for character, states in sorted(priors.items())
    ]


def _state_maps_payload(state_maps: Mapping[int, Mapping[int, str]]) -> list[list[Any]]:
    return [
        [int(character), [[int(state), str(allele)] for state, allele in sorted(states.items())]]
        for character, states in sorted(state_maps.items())
    ]


def _validate_character_frame(frame: Any, *, _owned: bool = False) -> Any:
    np, pd = _numpy_pandas()
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Cassiopeia character matrices must be pandas DataFrames.")
    if frame.empty or frame.shape[1] == 0:
        raise ValueError("Cassiopeia character matrices must contain cells and characters.")
    if not frame.index.is_unique or not frame.columns.is_unique:
        raise ValueError("Cassiopeia character axes must be unique.")
    cells = frame.index.tolist()
    columns = frame.columns.tolist()
    if any(not isinstance(value, str) or not value or value != value.strip() for value in cells + columns):
        raise ValueError("Cassiopeia character axes must be exact non-empty strings.")
    owned_values = frame.to_numpy(copy=False)
    if _owned and owned_values.dtype == np.dtype("int64"):
        if bool(np.any(owned_values < MISSING_STATE)):
            raise ValueError("Cassiopeia states must be -1, 0, or positive integers.")
        return frame
    values = frame.to_numpy(dtype=object, copy=True)
    integer_values = np.empty(values.shape, dtype=np.int64)
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = values[row_index, column_index]
            if isinstance(value, (bool, np.bool_)):
                raise ValueError("Cassiopeia states cannot be Boolean values.")
            try:
                integer = int(value)
                numeric = float(value)
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError("Cassiopeia states must be integers.") from error
            if not math.isfinite(numeric) or numeric != integer or integer < MISSING_STATE:
                raise ValueError("Cassiopeia states must be -1, 0, or positive integers.")
            integer_values[row_index, column_index] = integer
    return pd.DataFrame(integer_values, index=list(cells), columns=list(columns), dtype="int64")


def _validate_prior_closure(
    frame: Any,
    priors: Mapping[int, Mapping[int, float]],
    state_maps: Mapping[int, Mapping[int, str]],
) -> None:
    np, _ = _numpy_pandas()
    expected = set(range(frame.shape[1]))
    if set(priors) != expected or set(state_maps) != expected:
        raise ValueError("Cassiopeia prior/state-map axes do not close over the character matrix.")
    for position in range(frame.shape[1]):
        observed = {int(value) for value in np.unique(frame.iloc[:, position].to_numpy()) if int(value) > 0}
        if not observed.issubset(priors[position]) or not observed.issubset(state_maps[position]):
            raise ValueError("Cassiopeia priors and allele maps must cover every observed mutation state.")


class CassiopeiaCharacters:
    """Tamper-evident, defensive, process-local character artifact."""

    __slots__ = ("__matrices", "__priors", "__state_maps", "__metadata_json", "__fingerprint")

    def __init__(
        self,
        matrices: Mapping[str, Any],
        priors: Mapping[str, Mapping[Any, Mapping[Any, Any]]],
        state_maps: Mapping[str, Mapping[Any, Mapping[Any, Any]]],
        metadata: Mapping[str, Any],
        *,
        _token: object | None = None,
        _owned: bool = False,
    ) -> None:
        if _token is not _ARTIFACT_TOKEN:
            raise TypeError("CassiopeiaCharacters values can only be created by the audited preparation function.")
        lineage_ids = tuple(sorted(matrices))
        if set(priors) != set(lineage_ids) or set(state_maps) != set(lineage_ids):
            raise ValueError("Cassiopeia character artifact lineage payloads are inconsistent.")
        normalized_matrices: dict[str, Any] = {}
        normalized_priors: dict[str, dict[int, dict[int, float]]] = {}
        normalized_maps: dict[str, dict[int, dict[int, str]]] = {}
        lineage_hashes: dict[str, dict[str, str]] = {}
        for lineage_id in lineage_ids:
            if not isinstance(lineage_id, str) or not lineage_id or lineage_id != lineage_id.strip():
                raise ValueError("Cassiopeia lineage IDs must be exact non-empty strings.")
            frame = _validate_character_frame(matrices[lineage_id], _owned=_owned)
            lineage_priors = _normalize_priors(priors[lineage_id])
            lineage_maps = _normalize_state_maps(state_maps[lineage_id])
            _validate_prior_closure(frame, lineage_priors, lineage_maps)
            normalized_matrices[lineage_id] = frame
            normalized_priors[lineage_id] = lineage_priors
            normalized_maps[lineage_id] = lineage_maps
            lineage_hashes[lineage_id] = {
                "character_matrix_sha256": _frame_hash(frame),
                "priors_sha256": _hash_payload(_priors_payload(lineage_priors)),
                "state_maps_sha256": _hash_payload(_state_maps_payload(lineage_maps)),
            }
        normalized_metadata = _json_value(metadata)
        if not isinstance(normalized_metadata, dict):
            raise TypeError("Cassiopeia character metadata must be a mapping.")
        normalized_metadata.update(
            {
                "artifact_type": CHARACTERS_ARTIFACT_TYPE,
                "schema_version": CHARACTERS_SCHEMA_VERSION,
                "producer": PRODUCER,
                "lineage_ids": list(lineage_ids),
                "lineage_hashes": lineage_hashes,
            }
        )
        normalized_metadata.pop("fingerprint", None)
        fingerprint = _hash_payload(normalized_metadata)
        normalized_metadata["fingerprint"] = fingerprint
        self.__matrices = normalized_matrices
        self.__priors = normalized_priors
        self.__state_maps = normalized_maps
        self.__metadata_json = _canonical_json(normalized_metadata)
        self.__fingerprint = fingerprint

    @property
    def fingerprint(self) -> str:
        self._validate()
        return self.__fingerprint

    @property
    def lineage_ids(self) -> tuple[str, ...]:
        self._validate()
        return tuple(self.__matrices)

    @property
    def metadata(self) -> dict[str, Any]:
        self._validate()
        return json.loads(self.__metadata_json)

    @property
    def qc_records(self) -> tuple[dict[str, Any], ...]:
        return tuple(copy.deepcopy(self.metadata.get("qc_records", [])))

    def _validate(self) -> None:
        try:
            metadata = json.loads(self.__metadata_json)
        except (TypeError, json.JSONDecodeError) as error:
            raise ValueError("Cassiopeia character artifact metadata was corrupted.") from error
        if (
            metadata.get("artifact_type") != CHARACTERS_ARTIFACT_TYPE
            or metadata.get("schema_version") != CHARACTERS_SCHEMA_VERSION
            or metadata.get("producer") != PRODUCER
            or metadata.get("fingerprint") != self.__fingerprint
        ):
            raise ValueError("Cassiopeia character artifact provenance is invalid.")
        unsigned = dict(metadata)
        unsigned.pop("fingerprint", None)
        if _hash_payload(unsigned) != self.__fingerprint:
            raise ValueError("Cassiopeia character artifact metadata was tampered with.")
        if tuple(metadata.get("lineage_ids", [])) != tuple(self.__matrices):
            raise ValueError("Cassiopeia character artifact lineage identity was tampered with.")
        expected_hashes = metadata.get("lineage_hashes")
        if not isinstance(expected_hashes, dict):
            raise ValueError("Cassiopeia character artifact hashes are missing.")
        for lineage_id, frame in self.__matrices.items():
            lineage_hashes = expected_hashes.get(lineage_id, {})
            current = {
                "character_matrix_sha256": _frame_hash(frame),
                "priors_sha256": _hash_payload(_priors_payload(self.__priors[lineage_id])),
                "state_maps_sha256": _hash_payload(_state_maps_payload(self.__state_maps[lineage_id])),
            }
            if current != lineage_hashes:
                raise ValueError(f"Cassiopeia character artifact payload for lineage {lineage_id!r} was tampered with.")

    def copy_lineage(
        self, lineage_id: str, *, _owned: bool = False
    ) -> tuple[Any, dict[int, dict[int, float]], dict[int, dict[int, str]]]:
        self._validate()
        if lineage_id not in self.__matrices:
            available = ", ".join(self.__matrices) or "none"
            raise ValueError(f"Lineage {lineage_id!r} has no computable character payload; available: {available}.")
        if _owned:
            return self.__matrices[lineage_id], self.__priors[lineage_id], self.__state_maps[lineage_id]
        return (
            self.__matrices[lineage_id].copy(deep=True),
            copy.deepcopy(self.__priors[lineage_id]),
            copy.deepcopy(self.__state_maps[lineage_id]),
        )


class CassiopeiaTree:
    """Tamper-evident solved topology with defensive backend ownership."""

    __slots__ = ("__tree", "__metadata_json", "__fingerprint")

    def __init__(
        self,
        tree: Any,
        metadata: Mapping[str, Any],
        *,
        _token: object | None = None,
        _owned: bool = False,
    ) -> None:
        if _token is not _ARTIFACT_TOKEN:
            raise TypeError("CassiopeiaTree values can only be created by audited reconstruction.")
        if not hasattr(tree, "copy") or not callable(tree.copy):
            raise TypeError("Cassiopeia tree backend must provide the public copy method.")
        private_tree = tree if _owned else tree.copy()
        content = _tree_content(private_tree)
        normalized_metadata = _json_value(metadata)
        if not isinstance(normalized_metadata, dict):
            raise TypeError("Cassiopeia tree metadata must be a mapping.")
        normalized_metadata.update(
            {
                "artifact_type": TREE_ARTIFACT_TYPE,
                "schema_version": TREE_SCHEMA_VERSION,
                "producer": PRODUCER,
                "lineage_id": str(normalized_metadata.get("lineage_id", "")),
                "root": content["root"],
                "leaf_axis": content["leaves"],
                "character_axis": content["characters"],
                "topology_sha256": content["topology_sha256"],
                "character_matrix_sha256": content["character_matrix_sha256"],
                "priors_sha256": content["priors_sha256"],
                "node_count": content["node_count"],
                "edge_count": content["edge_count"],
                "max_depth": content["max_depth"],
            }
        )
        lineage_id = normalized_metadata["lineage_id"]
        if not lineage_id or lineage_id != lineage_id.strip():
            raise ValueError("Cassiopeia tree lineage identity is missing or invalid.")
        normalized_metadata.pop("fingerprint", None)
        fingerprint = _hash_payload(normalized_metadata)
        normalized_metadata["fingerprint"] = fingerprint
        self.__tree = private_tree
        self.__metadata_json = _canonical_json(normalized_metadata)
        self.__fingerprint = fingerprint

    def _validate(self) -> None:
        try:
            metadata = json.loads(self.__metadata_json)
        except (TypeError, json.JSONDecodeError) as error:
            raise ValueError("Cassiopeia tree artifact metadata was corrupted.") from error
        if (
            metadata.get("artifact_type") != TREE_ARTIFACT_TYPE
            or metadata.get("schema_version") != TREE_SCHEMA_VERSION
            or metadata.get("producer") != PRODUCER
            or metadata.get("fingerprint") != self.__fingerprint
        ):
            raise ValueError("Cassiopeia tree artifact provenance is invalid.")
        unsigned = dict(metadata)
        unsigned.pop("fingerprint", None)
        if _hash_payload(unsigned) != self.__fingerprint:
            raise ValueError("Cassiopeia tree artifact metadata was tampered with.")
        content = _tree_content(self.__tree)
        current = {
            "root": content["root"],
            "leaf_axis": content["leaves"],
            "character_axis": content["characters"],
            "topology_sha256": content["topology_sha256"],
            "character_matrix_sha256": content["character_matrix_sha256"],
            "priors_sha256": content["priors_sha256"],
            "node_count": content["node_count"],
            "edge_count": content["edge_count"],
            "max_depth": content["max_depth"],
        }
        if any(metadata.get(key) != value for key, value in current.items()):
            raise ValueError("Cassiopeia tree artifact topology, axes, or character payload was tampered with.")

    @property
    def fingerprint(self) -> str:
        self._validate()
        return self.__fingerprint

    @property
    def metadata(self) -> dict[str, Any]:
        self._validate()
        return json.loads(self.__metadata_json)

    @property
    def provenance(self) -> dict[str, Any]:
        return self.metadata

    @property
    def lineage_id(self) -> str:
        return str(self.metadata["lineage_id"])

    @property
    def root(self) -> str:
        return str(self.metadata["root"])

    @property
    def input_cells(self) -> int:
        return len(self.metadata["leaf_axis"])

    @property
    def character_count(self) -> int:
        return len(self.metadata["character_axis"])

    @property
    def topology_fingerprint(self) -> str:
        return str(self.metadata["topology_sha256"])

    def copy_tree(self, *, _owned: bool = False) -> Any:
        self._validate()
        return self.__tree if _owned else self.__tree.copy()


def _require_characters(value: Any) -> CassiopeiaCharacters:
    if type(value) is not CassiopeiaCharacters:
        raise TypeError("Expected an audited OPENBIO_CASSIOPEIA_CHARACTERS artifact.")
    value._validate()
    return value


def _require_exact_parameters(callable_object: Any, required: Sequence[str], label: str) -> None:
    try:
        parameters = inspect.signature(callable_object).parameters
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Could not inspect the pinned Cassiopeia {label} API.") from error
    missing = [name for name in required if name not in parameters]
    if missing:
        raise RuntimeError(f"Pinned Cassiopeia {label} API is incompatible; missing: {', '.join(missing)}.")


def _require_cassiopeia() -> Any:
    try:
        installed = importlib_metadata.version(CASSIOPEIA_DISTRIBUTION)
    except importlib_metadata.PackageNotFoundError as error:
        raise RuntimeError(
            f"Cassiopeia lineage nodes require {CASSIOPEIA_DISTRIBUTION}=={CASSIOPEIA_VERSION}; "
            "install the lineage optional dependency explicitly."
        ) from error
    if installed != CASSIOPEIA_VERSION:
        raise RuntimeError(
            f"Cassiopeia lineage nodes require exactly {CASSIOPEIA_DISTRIBUTION}=={CASSIOPEIA_VERSION}; "
            f"found {installed}."
        )
    cassiopeia = importlib.import_module("cassiopeia")
    if getattr(cassiopeia, "__version__", None) != CASSIOPEIA_VERSION:
        raise RuntimeError("The imported cassiopeia module does not match the pinned distribution version.")
    _require_exact_parameters(
        cassiopeia.pp.convert_alleletable_to_character_matrix,
        (
            "alleletable",
            "allele_rep_thresh",
            "missing_data_allele",
            "missing_data_state",
            "mutation_priors",
            "cut_sites",
            "collapse_duplicates",
        ),
        "allele-table conversion",
    )
    _require_exact_parameters(
        cassiopeia.pp.compute_empirical_indel_priors,
        ("allele_table", "grouping_variables", "cut_sites"),
        "empirical-prior",
    )
    _require_exact_parameters(
        cassiopeia.solver.VanillaGreedySolver,
        ("missing_data_classifier", "prior_transformation"),
        "VanillaGreedy constructor",
    )
    _require_exact_parameters(
        cassiopeia.solver.VanillaGreedySolver.solve,
        ("cassiopeia_tree", "collapse_mutationless_edges"),
        "VanillaGreedy solve",
    )
    _require_exact_parameters(
        cassiopeia.tl.compute_expansion_pvalues,
        ("tree", "min_clade_size", "min_depth", "copy"),
        "expansion test",
    )
    _require_exact_parameters(
        cassiopeia.tl.score_small_parsimony,
        ("cassiopeia_tree", "meta_item", "root", "infer_ancestral_states", "label_key"),
        "small parsimony",
    )
    return cassiopeia


@contextlib.contextmanager
def _captured_backend() -> Iterator[tuple[io.StringIO, list[warnings.WarningMessage]]]:
    stream = io.StringIO()
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
        warnings.simplefilter("always")
        yield stream, caught


def _captured_messages(stream: io.StringIO, caught: Sequence[warnings.WarningMessage]) -> list[str]:
    messages = [str(item.message) for item in caught]
    text = stream.getvalue().strip()
    if text:
        messages.append(text[:_CAPTURE_LIMIT] + ("…" if len(text) > _CAPTURE_LIMIT else ""))
    return list(dict.fromkeys(messages))


@contextlib.contextmanager
def _preserved_rng(seed: int = 0) -> Iterator[None]:
    np, _ = _numpy_pandas()
    with _CASSIOPEIA_LOCK:
        numpy_state = np.random.get_state()
        python_state = random.getstate()
        try:
            np.random.seed(int(seed))
            random.seed(int(seed))
            yield
        finally:
            np.random.set_state(numpy_state)
            random.setstate(python_state)


def _graph_structure(graph: Any) -> dict[str, Any]:
    if graph is None or not hasattr(graph, "nodes") or not hasattr(graph, "edges"):
        raise ValueError("Cassiopeia backend did not expose a public directed topology.")
    if hasattr(graph, "is_directed") and not bool(graph.is_directed()):
        raise ValueError("Cassiopeia topology must be directed.")
    nodes = list(graph.nodes)
    if not nodes:
        raise ValueError("Cassiopeia topology is empty.")
    if any(not isinstance(node, str) or not node or node != node.strip() for node in nodes):
        raise ValueError("Cassiopeia topology node identities must be exact non-empty strings.")
    if len(set(nodes)) != len(nodes):
        raise ValueError("Cassiopeia topology node identities must be unique.")
    edges: list[tuple[str, str, float]] = []
    parents: dict[str, list[str]] = {node: [] for node in nodes}
    children: dict[str, list[str]] = {node: [] for node in nodes}
    for parent, child, attributes in graph.edges(data=True):
        if parent not in children or child not in parents or parent == child:
            raise ValueError("Cassiopeia topology contains an invalid/self edge.")
        length_raw = attributes.get("length", 1.0)
        if isinstance(length_raw, bool):
            raise ValueError("Cassiopeia branch lengths cannot be Boolean.")
        length = float(length_raw)
        if not math.isfinite(length) or length < 0.0:
            raise ValueError("Cassiopeia branch lengths must be finite and non-negative.")
        parents[child].append(parent)
        children[parent].append(child)
        edges.append((str(parent), str(child), length))
    roots = [node for node in nodes if not parents[node]]
    if len(roots) != 1:
        raise ValueError("Cassiopeia topology must have exactly one root.")
    root = roots[0]
    if any(len(parents[node]) != 1 for node in nodes if node != root):
        raise ValueError("Cassiopeia topology must give every non-root node exactly one parent.")
    if len(edges) != len(nodes) - 1:
        raise ValueError("Cassiopeia topology is not a directed arborescence.")
    visiting: set[str] = set()
    visited: set[str] = set()
    depths: dict[str, int] = {}

    def visit(node: str, depth: int) -> None:
        if node in visiting:
            raise ValueError("Cassiopeia topology contains a cycle.")
        if node in visited:
            return
        visiting.add(node)
        depths[node] = depth
        for child in sorted(children[node]):
            visit(child, depth + 1)
        visiting.remove(node)
        visited.add(node)

    visit(root, 0)
    if visited != set(nodes):
        raise ValueError("Cassiopeia topology is disconnected from its root.")
    leaves = sorted(node for node in nodes if not children[node])
    canonical_edges = [[parent, child, length] for parent, child, length in sorted(edges)]
    topology_payload = {"root": root, "nodes": sorted(nodes), "edges": canonical_edges}
    return {
        "root": root,
        "nodes": sorted(nodes),
        "edges": canonical_edges,
        "parents": {node: tuple(parents[node]) for node in nodes},
        "children": {node: tuple(sorted(children[node])) for node in nodes},
        "leaves": leaves,
        "depths": depths,
        "max_depth": max(depths.values()),
        "topology_sha256": _hash_payload(topology_payload),
    }


def _tree_content(tree: Any) -> dict[str, Any]:
    if not hasattr(tree, "get_tree_topology") or not callable(tree.get_tree_topology):
        raise TypeError("Cassiopeia tree backend lacks the public get_tree_topology method.")
    graph = tree.get_tree_topology()
    structure = _graph_structure(graph)
    if str(getattr(tree, "root", "")) != structure["root"]:
        raise ValueError("Cassiopeia backend root disagrees with the public topology.")
    matrix = _validate_character_frame(getattr(tree, "character_matrix", None), _owned=True)
    cells = list(matrix.index)
    if set(cells) != set(structure["leaves"]):
        raise ValueError("Cassiopeia character-matrix cells do not exactly match topology leaves.")
    raw_priors = getattr(tree, "priors", None)
    if raw_priors is None:
        raw_priors = {}
    filled_priors = {position: dict(raw_priors.get(position, {})) for position in range(matrix.shape[1])}
    priors = _normalize_priors(filled_priors)
    observed_maps: dict[int, dict[int, str]] = {}
    for position in range(matrix.shape[1]):
        observed_maps[position] = {
            int(state): f"observed_state_{int(state)}"
            for state in sorted(set(matrix.iloc[:, position].tolist()))
            if int(state) > 0
        }
    _validate_prior_closure(matrix, priors, observed_maps)
    return {
        "root": structure["root"],
        "leaves": cells,
        "characters": list(matrix.columns),
        "topology_sha256": structure["topology_sha256"],
        "character_matrix_sha256": _frame_hash(matrix),
        "priors_sha256": _hash_payload(_priors_payload(priors)),
        "node_count": len(structure["nodes"]),
        "edge_count": len(structure["edges"]),
        "max_depth": structure["max_depth"],
    }


def _require_tree(value: Any) -> CassiopeiaTree:
    if type(value) is not CassiopeiaTree:
        raise TypeError("Expected an audited OPENBIO_CASSIOPEIA_TREE artifact.")
    value._validate()
    return value


def _canonicalize_solved_tree(
    cassiopeia: Any,
    tree: Any,
    matrix: Any,
    priors: Mapping[int, Mapping[int, float]],
    *,
    _worker_owned: bool = False,
) -> Any:
    nx = importlib.import_module("networkx")
    graph = tree.get_tree_topology()
    structure = _graph_structure(graph)
    expected_leaves = set(matrix.index)
    if set(structure["leaves"]) != expected_leaves:
        raise RuntimeError("VanillaGreedy output leaves do not match the selected character axis.")
    descendants: dict[str, tuple[str, ...]] = {}
    for node in sorted(structure["nodes"], key=lambda value: structure["depths"][value], reverse=True):
        children = structure["children"][node]
        if not children:
            descendants[node] = (node,)
        else:
            descendants[node] = tuple(sorted(leaf for child in children for leaf in descendants[child]))
    mapping: dict[str, str] = {}
    occupied = set(expected_leaves)
    for node in structure["nodes"]:
        if node in expected_leaves:
            continue
        digest = hashlib.sha256(_canonical_json(descendants[node]).encode("utf-8")).hexdigest()
        canonical = f"{_INTERNAL_PREFIX}{digest}"
        if canonical in occupied:
            raise RuntimeError("Canonical Cassiopeia internal-node identity collided with a leaf or subtree.")
        occupied.add(canonical)
        mapping[node] = canonical
    canonical_graph = nx.relabel_nodes(graph, mapping, copy=True)
    canonical_tree = cassiopeia.data.CassiopeiaTree(
        character_matrix=matrix if _worker_owned else matrix.copy(deep=True),
        missing_state_indicator=MISSING_STATE,
        priors=priors if _worker_owned else copy.deepcopy(dict(priors)),
        tree=canonical_graph,
    )
    _tree_content(canonical_tree)
    return canonical_tree


def _software_versions(openbio_version: str) -> dict[str, str]:
    versions = {"python": platform.python_version(), "openbio-singlecell": str(openbio_version)}
    for distribution in ("numpy", "pandas", CASSIOPEIA_DISTRIBUTION):
        try:
            versions[distribution] = importlib_metadata.version(distribution)
        except importlib_metadata.PackageNotFoundError:
            versions[distribution] = "not-installed"
    return versions


def _strict_summary(
    *,
    node_id: str,
    methods: str,
    results: str,
    key_results: Mapping[str, Any],
    parameters: Mapping[str, Any],
    warnings_list: Sequence[str],
    limitations: Sequence[str],
    references: Sequence[Mapping[str, Any]],
    openbio_version: str,
) -> dict[str, Any]:
    summary = {
        "schema_version": 1,
        "node_id": node_id,
        "methods": methods,
        "results": results,
        "key_results": _json_value(key_results),
        "parameters": _json_value(parameters),
        "warnings": [str(value) for value in warnings_list],
        "limitations": [str(value) for value in limitations],
        "references": [_json_value(value) for value in references],
        "software_versions": _software_versions(openbio_version),
    }
    json.dumps(summary, allow_nan=False)
    return summary


def _cassiopeia_reference() -> dict[str, Any]:
    return {
        "citation": "Jones MG et al. Inference of single-cell phylogenies from lineage tracing data using Cassiopeia. Genome Biology. 2020.",
        "doi": "10.1186/s13059-020-02000-8",
        "url": "https://doi.org/10.1186/s13059-020-02000-8",
        "kind": "method",
    }


def _software_reference() -> dict[str, Any]:
    return {
        "citation": f"Cassiopeia {CASSIOPEIA_VERSION} API; cassiopeia-mt redistribution tag {CASSIOPEIA_SOURCE_TAG}.",
        "url": f"https://github.com/andrecossa5/Cassiopeia/tree/{CASSIOPEIA_SOURCE_TAG}",
        "kind": "software",
    }


def _yang_reference() -> dict[str, Any]:
    return {
        "citation": "Yang D, Jones MG et al. Lineage tracing reveals the phylodynamics, plasticity, and paths of tumor evolution. Cell. 2022.",
        "doi": "10.1016/j.cell.2022.04.015",
        "url": "https://doi.org/10.1016/j.cell.2022.04.015",
        "kind": "method",
    }


def _standalone_source() -> str:
    try:
        source = Path(__file__).read_text(encoding="utf-8")
    except (NameError, OSError):
        try:
            source = inspect.getsource(sys.modules[__name__])
        except (KeyError, OSError, TypeError):
            source = "def run(*args, **kwargs):\n    raise RuntimeError('Save/import this generated source first.')\n"
    compile(source, "<openbio-cassiopeia-code>", "exec")
    return source.rstrip() + "\n"


def _parse_column_list(value: str, label: str) -> tuple[str, ...]:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a comma-separated string.")
    columns = tuple(item.strip() for item in value.split(","))
    if not columns or any(not item for item in columns):
        raise ValueError(f"{label} must contain exact non-empty column names.")
    if len(set(columns)) != len(columns):
        raise ValueError(f"{label} contains duplicate column names.")
    return columns


def _validate_fraction(value: float, label: str, *, lower_open: bool = False) -> float:
    numeric = float(value)
    valid_lower = numeric > 0.0 if lower_open else numeric >= 0.0
    if not math.isfinite(numeric) or not valid_lower or numeric > 1.0:
        bound = "(0, 1]" if lower_open else "[0, 1]"
        raise ValueError(f"{label} must be a finite value in {bound}.")
    return numeric


def _read_allele_table(
    path: str, first_column_as_index: bool, max_file_mib: int = 512
) -> tuple[Any, dict[str, Any]]:
    _, pd = _numpy_pandas()
    if not isinstance(path, str) or not path:
        raise ValueError("Allele-table path cannot be empty.")
    suffix = Path(path).suffix.lower()
    if suffix not in ALLELE_TABLE_EXTENSIONS:
        raise ValueError("Cassiopeia allele tables must be CSV or TSV files.")
    if int(max_file_mib) <= 0:
        raise ValueError("max_file_mib must be positive.")
    limit = int(max_file_mib) * 1024 * 1024
    with open(path, "rb") as handle:
        content = handle.read(limit + 1)
    if len(content) > limit:
        raise ValueError(f"Allele table exceeds max_file_mib={int(max_file_mib)}.")
    try:
        text = content.decode("utf-8-sig", errors="strict")
    except UnicodeDecodeError as error:
        raise ValueError("Allele table must be valid UTF-8 text.") from error
    delimiter = "," if suffix == ".csv" else "\t"
    rows = csv.reader(io.StringIO(text), delimiter=delimiter, strict=True)
    try:
        header = next(rows)
    except StopIteration as error:
        raise ValueError("Allele table is empty.") from error
    data_header = header[1:] if first_column_as_index else header
    if first_column_as_index and len(header) < 2:
        raise ValueError("Indexed allele table must contain at least one data column.")
    if any(not column for column in data_header):
        raise ValueError("Allele-table data headers cannot be blank.")
    if len(set(data_header)) != len(data_header):
        raise ValueError("Allele-table headers must be unique before pandas parsing.")
    row_count = 0
    for line_number, row in enumerate(rows, start=2):
        if not row:
            continue
        if len(row) != len(header):
            raise ValueError(f"Allele table row {line_number} has {len(row)} fields; expected {len(header)}.")
        row_count += 1
    if row_count == 0:
        raise ValueError("Allele table contains no records.")
    frame = pd.read_csv(
        io.StringIO(text),
        sep=delimiter,
        index_col=0 if first_column_as_index else None,
        dtype=object,
        keep_default_na=False,
        na_filter=False,
    ).reset_index(drop=True)
    if tuple(str(column) for column in frame.columns) != tuple(data_header):
        raise ValueError("Allele-table columns changed during parsing.")
    return frame, {
        "name": os.path.basename(path),
        "size_bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "format": suffix.lstrip("."),
        "record_count": row_count,
    }


def _canonical_row_order(frame: Any, columns: Sequence[str]) -> Any:
    keys = []
    for position, row in enumerate(frame.loc[:, list(columns)].itertuples(index=False, name=None)):
        key = tuple("<CSV-NULL>" if _json_scalar(value) is None else str(value) for value in row)
        keys.append((key, position))
    positions = [position for _, position in sorted(keys)]
    return frame.iloc[positions].reset_index(drop=True)


def _backend_character_labels(frame: Any, cut_sites: Sequence[str], threshold: float) -> list[str]:
    np, pd = _numpy_pandas()
    key_pairs: dict[str, tuple[str, str]] = {}
    allele_values: dict[str, list[str]] = {}
    for cell_id in frame["cellBC"].drop_duplicates().tolist():
        cell_frame = frame.loc[frame["cellBC"] == cell_id]
        for row in cell_frame.itertuples(index=False, name=None):
            record = dict(zip(frame.columns, row, strict=True))
            intbc = str(record["intBC"])
            for cut_site in cut_sites:
                key = f"{intbc}{cut_site}"
                pair = (intbc, cut_site)
                if key in key_pairs and key_pairs[key] != pair:
                    raise ValueError("Integration-barcode/cut-site concatenation is ambiguous for Cassiopeia conversion.")
                key_pairs[key] = pair
                value = record[cut_site]
                missing = pd.isna(value)
                normalized = "<CSV-NULL>" if isinstance(missing, (bool, np.bool_)) and missing else str(value)
                allele_values.setdefault(key, []).append(normalized)
    labels = []
    for key, values in allele_values.items():
        counts = Counter(values)
        if any(count / len(values) > threshold for count in counts.values()):
            continue
        intbc, cut_site = key_pairs[key]
        labels.append(f"{intbc}::{cut_site}")
    return labels


def _validate_identifier_column(frame: Any, column: str) -> None:
    for row_number, value in enumerate(frame[column].tolist(), start=2):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{column} row {row_number} must be an exact non-empty string.")
        if any(character in value for character in ("\x00", "\n", "\r", ",")):
            raise ValueError(f"{column} row {row_number} contains a forbidden control/comma character.")
        if value.startswith((_INTERNAL_PREFIX, _GROUP_PREFIX)):
            raise ValueError(f"{column} row {row_number} uses a reserved OpenBio node prefix.")


def _normalize_alleles(frame: Any, cut_sites: Sequence[str], missing_data_allele: str) -> Any:
    np, pd = _numpy_pandas()
    token = missing_data_allele
    if not isinstance(token, str):
        raise TypeError("missing_data_allele must be a string; blank selects CSV-null-only missingness.")
    if token and (token != token.strip() or "none" in token.lower()):
        raise ValueError("The missing allele must be exact and distinct from Cassiopeia's reserved uncut token.")
    normalized = frame.copy(deep=True)
    for cut_site in cut_sites:
        values = []
        for row_number, value in enumerate(normalized[cut_site].tolist(), start=2):
            if not isinstance(value, str):
                raise ValueError(f"Allele {cut_site} row {row_number} must be a string or CSV null.")
            if value == "" or (token and value == token):
                values.append(np.nan)
                continue
            if value != value.strip() or any(character in value for character in ("\x00", "\n", "\r")):
                raise ValueError(f"Allele {cut_site} row {row_number} contains whitespace/control data.")
            if value in {"NONE", "None"}:
                values.append("NONE")
                continue
            if "none" in value.lower():
                raise ValueError(
                    f"Allele {cut_site} row {row_number} contains reserved text 'None' but is not exactly NONE/None."
                )
            values.append(value)
        normalized[cut_site] = pd.Series(values, index=normalized.index, dtype=object)
    return normalized


def _allele_conflicts(frame: Any, cut_sites: Sequence[str]) -> list[str]:
    _, pd = _numpy_pandas()
    conflicts: list[str] = []
    for identity, group in frame.groupby(
        ["lineage_id", "cellBC", "intBC"], observed=True, sort=True, dropna=False
    ):
        for cut_site in cut_sites:
            keys = {"<CSV-NULL>" if pd.isna(value) else str(value) for value in group[cut_site].tolist()}
            if len(keys) > 1:
                conflicts.append(f"{identity!r}/{cut_site}")
    return conflicts


def _remap_character_payload(
    frame: Any,
    priors: Mapping[Any, Mapping[Any, Any]],
    state_maps: Mapping[Any, Mapping[Any, Any]],
    retained_columns: Sequence[str],
) -> tuple[Any, dict[int, dict[int, float]], dict[int, dict[int, str]]]:
    old_positions = {str(column): position for position, column in enumerate(frame.columns)}
    trimmed = frame.loc[:, list(retained_columns)].copy(deep=True)
    remapped_priors: dict[int, dict[int, float]] = {}
    remapped_maps: dict[int, dict[int, str]] = {}
    for new_position, column in enumerate(retained_columns):
        old_position = old_positions[str(column)]
        remapped_priors[new_position] = dict(priors.get(old_position, {}))
        remapped_maps[new_position] = dict(state_maps.get(old_position, {}))
    normalized_priors = _normalize_priors(remapped_priors)
    normalized_maps = _normalize_state_maps(remapped_maps)
    _validate_prior_closure(trimmed, normalized_priors, normalized_maps)
    return trimmed, normalized_priors, normalized_maps


def _convert_lineage(
    *,
    cassiopeia: Any,
    lineage_name: str,
    lineage_frame: Any,
    cut_sites: Sequence[str],
    normalized_groups: Sequence[str],
    empirical_priors: Any,
    threshold: float,
    minimum_cells: int,
    max_missing: float,
    max_uncut: float,
    min_unique: float,
    min_informative: float,
) -> tuple[
    dict[str, Any],
    tuple[Any, dict[int, dict[int, float]], dict[int, dict[int, str]]] | None,
    list[str],
]:
    np, pd = _numpy_pandas()
    canonical_columns = list(dict.fromkeys(["cellBC", "intBC", *cut_sites, *normalized_groups]))
    lineage_frame = _canonical_row_order(lineage_frame.reset_index(drop=True), canonical_columns)
    labels = _backend_character_labels(lineage_frame, cut_sites, threshold)
    original_cells = int(lineage_frame["cellBC"].nunique())
    qc_warnings: list[str] = []
    if original_cells < minimum_cells:
        qc_warnings.append("fewer_than_minimum_input_cells")
    conversion_before = _frame_hash(lineage_frame)
    with _captured_backend() as (convert_stream, convert_caught):
        character_matrix, lineage_priors, lineage_state_maps = (
            cassiopeia.pp.convert_alleletable_to_character_matrix(
                lineage_frame,
                allele_rep_thresh=threshold,
                missing_data_allele=None,
                missing_data_state=MISSING_STATE,
                mutation_priors=empirical_priors,
                cut_sites=list(cut_sites),
                collapse_duplicates=True,
            )
        )
    messages = _captured_messages(convert_stream, convert_caught)
    if _frame_hash(lineage_frame) != conversion_before:
        raise RuntimeError(f"Cassiopeia conversion mutated lineage {lineage_name!r} input.")
    if not isinstance(character_matrix, pd.DataFrame):
        raise RuntimeError("Cassiopeia conversion did not return a pandas character matrix.")
    if character_matrix.shape[1] != len(labels):
        raise RuntimeError("Cassiopeia conversion character identity does not match audited cassette ordering.")
    total_characters = int(character_matrix.shape[1])
    if total_characters == 0:
        record = {
            "lineage_id": lineage_name,
            "status": "unavailable",
            "qc_warnings": list(dict.fromkeys(qc_warnings)),
            "unavailability_reasons": [
                "no_characters_after_representation_filter",
                "no_mutation_bearing_characters",
            ],
            "input_cells": original_cells,
            "retained_cells": 0,
            "input_characters": 0,
            "retained_characters": 0,
            "informative_character_fraction": 0.0,
            "unique_profile_fraction": 0.0,
            "mean_missing_fraction": None,
            "mean_uncut_fraction_observed": None,
            "all_missing_cells": 0,
            "all_uncut_cells": 0,
            "missing_threshold_filtered_cells": 0,
            "uncut_threshold_filtered_cells": 0,
        }
        return record, None, messages
    if not character_matrix.index.is_unique or set(map(str, character_matrix.index)) != set(
        lineage_frame["cellBC"].unique()
    ):
        raise RuntimeError("Cassiopeia conversion changed cell identity.")
    character_matrix = character_matrix.copy(deep=True)
    character_matrix.index = [str(value) for value in character_matrix.index]
    character_matrix.columns = labels
    character_matrix = _validate_character_frame(character_matrix.sort_index(kind="stable"))
    keep_cells: list[str] = []
    all_missing_cells = 0
    all_uncut_cells = 0
    missing_filtered_cells = 0
    uncut_filtered_cells = 0
    cell_missing: list[float] = []
    cell_uncut: list[float] = []
    for cell_id, row in character_matrix.iterrows():
        values = row.to_numpy(dtype="int64", copy=True)
        missing_count = int(np.count_nonzero(values == MISSING_STATE))
        observed = values[values != MISSING_STATE]
        mutation_count = int(np.count_nonzero(observed > UNCUT_STATE))
        missing_fraction = missing_count / total_characters
        uncut_fraction = None
        if observed.size:
            uncut_fraction = float(np.count_nonzero(observed == UNCUT_STATE) / observed.size)
        cell_missing.append(float(missing_fraction))
        if uncut_fraction is not None:
            cell_uncut.append(uncut_fraction)
        if observed.size == 0:
            all_missing_cells += 1
        if mutation_count == 0:
            all_uncut_cells += 1
        if missing_fraction > max_missing:
            missing_filtered_cells += 1
            continue
        if uncut_fraction is not None and uncut_fraction > max_uncut:
            uncut_filtered_cells += 1
            continue
        keep_cells.append(str(cell_id))
    filtered = character_matrix.loc[keep_cells].copy(deep=True)
    if filtered.shape[0] < minimum_cells:
        qc_warnings.append("fewer_than_minimum_retained_cells")
    informative_columns = [
        str(column) for column in filtered.columns if bool(np.any(filtered[column].to_numpy(dtype="int64") > 0))
    ]
    informative_fraction = len(informative_columns) / total_characters
    if not informative_columns:
        qc_warnings.append("no_mutation_bearing_characters")
    if informative_fraction < min_informative:
        qc_warnings.append("below_minimum_informative_character_fraction")
    filtered = filtered.loc[:, informative_columns]
    unique_fraction = (
        float(filtered.drop_duplicates().shape[0] / filtered.shape[0]) if filtered.shape[0] and filtered.shape[1] else 0.0
    )
    if unique_fraction < min_unique:
        qc_warnings.append("below_minimum_unique_fraction")
    unavailable: list[str] = []
    if filtered.shape[0] == 0:
        unavailable.append("no_cells_after_declared_cell_filters")
    if not informative_columns:
        unavailable.append("no_mutation_bearing_characters")
    status = "unavailable" if unavailable else ("warning" if qc_warnings else "pass")
    record = {
        "lineage_id": lineage_name,
        "status": status,
        "qc_warnings": list(dict.fromkeys(qc_warnings)),
        "unavailability_reasons": list(dict.fromkeys(unavailable)),
        "input_cells": original_cells,
        "retained_cells": int(filtered.shape[0]),
        "input_characters": total_characters,
        "retained_characters": int(filtered.shape[1]),
        "informative_character_fraction": float(informative_fraction),
        "unique_profile_fraction": float(unique_fraction),
        "mean_missing_fraction": float(np.mean(cell_missing)) if cell_missing else None,
        "mean_uncut_fraction_observed": float(np.mean(cell_uncut)) if cell_uncut else None,
        "all_missing_cells": all_missing_cells,
        "all_uncut_cells": all_uncut_cells,
        "missing_threshold_filtered_cells": missing_filtered_cells,
        "uncut_threshold_filtered_cells": uncut_filtered_cells,
    }
    if status == "unavailable":
        return record, None, messages
    payload = _remap_character_payload(
        character_matrix.loc[keep_cells], lineage_priors, lineage_state_maps, informative_columns
    )
    return record, payload, messages


def prepare_cassiopeia_characters(
    allele_table_path: str,
    *,
    first_column_as_index: bool = True,
    lineage_column: str = "Tumor",
    cell_barcode_column: str = "cellBC",
    integration_barcode_column: str = "intBC",
    cut_site_columns: str = "r1,r2,r3",
    prior_grouping_columns: str = "Tumor,intBC",
    missing_data_allele: str = "",
    allele_representation_threshold: float = 0.98,
    minimum_cells: int = 2,
    maximum_missing_fraction: float = 0.8,
    maximum_uncut_fraction: float = 0.8,
    minimum_unique_fraction: float = 0.05,
    minimum_informative_character_fraction: float = 0.2,
    max_file_mib: int = 512,
    max_matrix_gib: float = 2.0,
    openbio_version: str = "standalone",
) -> tuple[CassiopeiaCharacters, Any, dict[str, Any], str]:
    _, pd = _numpy_pandas()
    cassiopeia = _require_cassiopeia()
    threshold = _validate_fraction(
        allele_representation_threshold, "allele_representation_threshold", lower_open=True
    )
    max_missing = _validate_fraction(maximum_missing_fraction, "maximum_missing_fraction")
    max_uncut = _validate_fraction(maximum_uncut_fraction, "maximum_uncut_fraction")
    min_unique = _validate_fraction(minimum_unique_fraction, "minimum_unique_fraction")
    min_informative = _validate_fraction(
        minimum_informative_character_fraction, "minimum_informative_character_fraction"
    )
    minimum_cells = int(minimum_cells)
    if minimum_cells < 1:
        raise ValueError("minimum_cells must be positive.")
    if not math.isfinite(float(max_matrix_gib)) or float(max_matrix_gib) <= 0:
        raise ValueError("max_matrix_gib must be positive and finite.")
    cut_sites = _parse_column_list(cut_site_columns, "cut_site_columns")
    prior_groups = _parse_column_list(prior_grouping_columns, "prior_grouping_columns")
    identifier_sources = (lineage_column, cell_barcode_column, integration_barcode_column)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in identifier_sources):
        raise ValueError("Identifier column names must be exact non-empty strings.")
    if len(set(identifier_sources)) != 3:
        raise ValueError("Lineage, cell-barcode, and integration-barcode columns must be distinct.")
    if set(cut_sites) & set(identifier_sources):
        raise ValueError("Cut-site columns cannot also be identifier columns.")
    policy_warnings: list[str] = []
    overlapping_prior_groups = sorted(set(cut_sites) & set(prior_groups))
    if overlapping_prior_groups:
        policy_warnings.append(
            "Empirical-prior grouping includes cut-site columns; this is an expert-selected estimand and may make "
            "mutation priors degenerate: " + ", ".join(overlapping_prior_groups) + "."
        )

    raw, source = _read_allele_table(allele_table_path, bool(first_column_as_index), int(max_file_mib))
    required = set(identifier_sources) | set(cut_sites) | set(prior_groups)
    missing_columns = sorted(required - set(raw.columns))
    if missing_columns:
        raise ValueError(f"Allele table is missing required columns: {', '.join(missing_columns)}.")
    rename = {
        lineage_column: "lineage_id",
        cell_barcode_column: "cellBC",
        integration_barcode_column: "intBC",
    }
    for source_name, target_name in rename.items():
        if target_name in raw.columns and source_name != target_name:
            raise ValueError(f"Allele-table normalization would collide with existing column {target_name!r}.")
    normalized_groups = tuple(rename.get(column, column) for column in prior_groups)
    if "lineage_id" not in normalized_groups or "intBC" not in normalized_groups:
        policy_warnings.append(
            "Empirical-prior grouping does not include both lineage and integration-barcode identity; the selected "
            "grouping was honored, but independence and allele-frequency assumptions require expert review."
        )
    frame = raw.rename(columns=rename).copy(deep=True)
    for column in ("lineage_id", "cellBC", "intBC"):
        _validate_identifier_column(frame, column)
    for column in dict.fromkeys(normalized_groups):
        if column not in {"lineage_id", "cellBC", "intBC"} and bool(frame[column].isna().any()):
            policy_warnings.append(
                f"Empirical-prior grouping column {column!r} contains missing values; backend grouping behavior was "
                "honored and the resulting estimand requires expert review."
            )
    frame = _normalize_alleles(frame, cut_sites, missing_data_allele)
    lineage_counts = frame.groupby("cellBC", observed=True, sort=True)["lineage_id"].nunique()
    cross_lineage_cells = sorted(str(value) for value in lineage_counts[lineage_counts > 1].index.tolist())
    if cross_lineage_cells:
        raise ValueError(
            "Cell barcode identity is not globally unique across lineages: " + ", ".join(cross_lineage_cells[:10])
        )
    conflicts = _allele_conflicts(frame, cut_sites)
    if conflicts:
        raise ValueError("Conflicting repeated cell/intBC/cut-site observations: " + ", ".join(conflicts[:10]))
    relevant_columns = list(dict.fromkeys(["lineage_id", "cellBC", "intBC", *normalized_groups, *cut_sites]))
    frame = _canonical_row_order(frame.loc[:, relevant_columns], relevant_columns)
    input_records = int(frame.shape[0])
    frame = frame.drop_duplicates(ignore_index=True)
    duplicate_records = input_records - int(frame.shape[0])
    potential_cells = int(frame["cellBC"].nunique())
    potential_intbcs = int(frame["intBC"].nunique())
    potential_bytes = potential_cells * potential_intbcs * len(cut_sites) * 16
    if potential_bytes > float(max_matrix_gib) * 1024**3:
        raise ValueError(
            f"Conservative character-matrix estimate {potential_bytes / 1024**3:.3f} GiB exceeds max_matrix_gib."
        )

    prior_frame = frame.copy(deep=True)
    for cut_site in cut_sites:
        prior_frame[cut_site] = prior_frame[cut_site].where(prior_frame[cut_site].notna(), "NONE")
    prior_before = _frame_hash(prior_frame)
    with _captured_backend() as (prior_stream, prior_caught):
        empirical_priors = cassiopeia.pp.compute_empirical_indel_priors(
            prior_frame, grouping_variables=list(normalized_groups), cut_sites=list(cut_sites)
        )
    if _frame_hash(prior_frame) != prior_before:
        raise RuntimeError("Cassiopeia empirical-prior computation mutated its input table.")
    backend_messages = _captured_messages(prior_stream, prior_caught)
    if not isinstance(empirical_priors, pd.DataFrame) or "freq" not in empirical_priors.columns:
        raise RuntimeError("Cassiopeia empirical-prior computation returned an incompatible table.")
    if not empirical_priors.index.is_unique:
        raise RuntimeError("Cassiopeia empirical-prior alleles are not unique.")
    for allele, probability in empirical_priors["freq"].items():
        if not isinstance(allele, str) or not allele or "none" in allele.lower():
            raise RuntimeError("Cassiopeia empirical priors contain an invalid mutation allele.")
        numeric = float(probability)
        if not math.isfinite(numeric) or not 0.0 < numeric <= 1.0:
            raise RuntimeError("Cassiopeia empirical priors contain invalid probabilities.")

    matrices: dict[str, Any] = {}
    priors_by_lineage: dict[str, dict[int, dict[int, float]]] = {}
    maps_by_lineage: dict[str, dict[int, dict[int, str]]] = {}
    qc_records: list[dict[str, Any]] = []
    for lineage_id, lineage_frame in frame.groupby("lineage_id", observed=True, sort=True):
        record, payload, messages = _convert_lineage(
            cassiopeia=cassiopeia,
            lineage_name=str(lineage_id),
            lineage_frame=lineage_frame,
            cut_sites=cut_sites,
            normalized_groups=normalized_groups,
            empirical_priors=empirical_priors,
            threshold=threshold,
            minimum_cells=minimum_cells,
            max_missing=max_missing,
            max_uncut=max_uncut,
            min_unique=min_unique,
            min_informative=min_informative,
        )
        qc_records.append(record)
        backend_messages.extend(messages)
        if payload is not None:
            matrix, lineage_priors, lineage_maps = payload
            matrices[str(lineage_id)] = matrix
            priors_by_lineage[str(lineage_id)] = lineage_priors
            maps_by_lineage[str(lineage_id)] = lineage_maps

    qc_table = pd.DataFrame.from_records(qc_records).sort_values("lineage_id", kind="stable").reset_index(drop=True)
    parameters = {
        "first_column_as_index": bool(first_column_as_index),
        "lineage_column": lineage_column,
        "cell_barcode_column": cell_barcode_column,
        "integration_barcode_column": integration_barcode_column,
        "cut_site_columns": list(cut_sites),
        "prior_grouping_columns": list(prior_groups),
        "missing_data_allele": missing_data_allele or None,
        "allele_representation_threshold": threshold,
        "minimum_cells": minimum_cells,
        "maximum_missing_fraction": max_missing,
        "maximum_uncut_fraction": max_uncut,
        "minimum_unique_fraction": min_unique,
        "minimum_informative_character_fraction": min_informative,
        "max_file_mib": int(max_file_mib),
        "max_matrix_gib": float(max_matrix_gib),
    }
    artifact = CassiopeiaCharacters(
        matrices,
        priors_by_lineage,
        maps_by_lineage,
        {
            "dependency": {
                "distribution": CASSIOPEIA_DISTRIBUTION,
                "version": CASSIOPEIA_VERSION,
                "source_tag": CASSIOPEIA_SOURCE_TAG,
                "source_commit": CASSIOPEIA_SOURCE_COMMIT,
                "sdist_sha256": CASSIOPEIA_SDIST_SHA256,
            },
            "source": source,
            "parameters": parameters,
            "qc_records": qc_records,
            "qc_table_sha256": _frame_hash(qc_table),
            "input_records": input_records,
            "identical_duplicate_records_removed": duplicate_records,
        },
        _token=_ARTIFACT_TOKEN,
    )
    passing = int(sum(record["status"] == "pass" for record in qc_records))
    warning_count = int(sum(record["status"] == "warning" for record in qc_records))
    unavailable = int(sum(record["status"] == "unavailable" for record in qc_records))
    available = passing + warning_count
    summary_warnings = list(dict.fromkeys((*policy_warnings, *backend_messages)))
    if duplicate_records:
        summary_warnings.append(
            f"Removed {duplicate_records} identical repeated allele-table records after conflict audit."
        )
    if warning_count:
        summary_warnings.append(
            f"{warning_count} computable lineage(s) did not meet at least one declared QC threshold; they remain "
            "available to expert downstream analysis with the flags disclosed."
        )
    if unavailable:
        summary_warnings.append(
            f"{unavailable} lineage(s) had no character payload satisfying backend/computability requirements."
        )
    results = (
        f"Prepared {available} computable lineage character payload(s) from {len(qc_records)} lineage(s) and "
        f"{potential_cells} globally identified cells; {warning_count} available lineage(s) carry QC warnings and "
        f"{unavailable} lineage(s) were computationally unavailable."
    )
    summary = _strict_summary(
        node_id="OpenBioSingleCellCassiopeiaLineageQC",
        methods=(
            "The bounded UTF-8 allele table was conflict-audited and canonically ordered. Cassiopeia public "
            "compute_empirical_indel_priors and convert_alleletable_to_character_matrix APIs encoded exact cut sites "
            "with -1 missing, 0 uncut, and positive mutation states. Missingness used all characters as denominator; "
            "uncut fractions used observed characters only. Declared missing/uncut cell filters were applied exactly; "
            "lineage-level minimum-cell, informative-character, and unique-profile thresholds were reported as QC "
            "warnings rather than downstream execution gates whenever a character payload remained computable."
        ),
        results=results,
        key_results={
            "lineages_total": len(qc_records),
            "lineages_passed": passing,
            "lineages_with_qc_warnings": warning_count,
            "lineages_available": available,
            "lineages_unavailable": unavailable,
            "lineages_failed": unavailable,
            "globally_unique_cells": potential_cells,
            "identical_duplicate_records_removed": duplicate_records,
            "artifact_fingerprint": artifact.fingerprint,
            "source_sha256": source["sha256"],
            "qc_records": qc_records,
        },
        parameters=parameters,
        warnings_list=summary_warnings,
        limitations=(
            "Lineage cells are observational descendants within clonal populations, not independent biological samples.",
            "Empirical priors depend on the declared lineage/integration grouping and allele table supplied.",
            "QC thresholds are protocol-dependent screening rules and do not establish biological validity.",
        ),
        references=(_cassiopeia_reference(), _software_reference()),
        openbio_version=openbio_version,
    )
    return artifact, qc_table, summary, _standalone_source()


def reconstruct_cassiopeia_tree(
    characters: CassiopeiaCharacters,
    lineage_id: str,
    *,
    prior_transformation: str = "negative_log",
    collapse_mutationless_edges: bool = False,
    openbio_version: str = "standalone",
    _worker_owned: bool = False,
) -> tuple[CassiopeiaTree, dict[str, Any], str]:
    artifact = _require_characters(characters)
    if not isinstance(lineage_id, str) or not lineage_id or lineage_id != lineage_id.strip():
        raise ValueError("lineage_id must be an exact non-empty canonical lineage ID.")
    transformations = {"negative_log", "inverse", "square_root_inverse"}
    if prior_transformation not in transformations:
        raise ValueError("prior_transformation must be negative_log, inverse, or square_root_inverse.")
    cassiopeia = _require_cassiopeia()
    artifact_metadata = artifact.metadata
    dependency = artifact_metadata.get("dependency", {})
    if (
        dependency.get("distribution") != CASSIOPEIA_DISTRIBUTION
        or dependency.get("version") != CASSIOPEIA_VERSION
        or dependency.get("source_commit") != CASSIOPEIA_SOURCE_COMMIT
    ):
        raise ValueError("Character artifact dependency provenance does not match the pinned Cassiopeia runtime.")
    upstream_fingerprint = artifact.fingerprint
    matrix, priors, _ = artifact.copy_lineage(lineage_id, _owned=_worker_owned)
    if matrix.shape[0] < 2:
        raise ValueError("VanillaGreedy reconstruction requires at least two character-matrix cells.")
    qc_record = next(
        (record for record in artifact.qc_records if record.get("lineage_id") == lineage_id),
        None,
    )
    backend_tree = cassiopeia.data.CassiopeiaTree(
        character_matrix=matrix if _worker_owned else matrix.copy(deep=True),
        missing_state_indicator=MISSING_STATE,
        priors=priors if _worker_owned else copy.deepcopy(priors),
    )
    constructor_parameters = inspect.signature(cassiopeia.solver.VanillaGreedySolver).parameters
    missing_classifier = constructor_parameters["missing_data_classifier"].default
    if not callable(missing_classifier) or getattr(missing_classifier, "__name__", "") != "assign_missing_average":
        raise RuntimeError("Pinned VanillaGreedy default is not the documented assign_missing_average callable.")
    solver = cassiopeia.solver.VanillaGreedySolver(
        missing_data_classifier=missing_classifier,
        prior_transformation=prior_transformation,
    )
    matrix_before = _frame_hash(backend_tree.character_matrix)
    priors_before = _hash_payload(_priors_payload(_normalize_priors(backend_tree.priors)))
    with _preserved_rng(0), _captured_backend() as (solver_stream, solver_caught):
        solver.solve(backend_tree, collapse_mutationless_edges=bool(collapse_mutationless_edges))
    messages = _captured_messages(solver_stream, solver_caught)
    if isinstance(qc_record, Mapping) and qc_record.get("status") == "warning":
        messages.append(
            "The selected lineage remained computable but did not meet all declared upstream QC thresholds: "
            + ", ".join(map(str, qc_record.get("qc_warnings", [])))
            + "."
        )
    if _frame_hash(backend_tree.character_matrix) != matrix_before:
        raise RuntimeError("VanillaGreedy mutated the selected character matrix.")
    if _hash_payload(_priors_payload(_normalize_priors(backend_tree.priors))) != priors_before:
        raise RuntimeError("VanillaGreedy mutated empirical priors.")
    artifact._validate()
    if artifact.fingerprint != upstream_fingerprint:
        raise RuntimeError("Source character artifact changed during reconstruction.")
    canonical_tree = _canonicalize_solved_tree(
        cassiopeia,
        backend_tree,
        matrix,
        priors,
        _worker_owned=_worker_owned,
    )
    content = _tree_content(canonical_tree)
    parameters = {
        "lineage_id": lineage_id,
        "prior_transformation": prior_transformation,
        "collapse_mutationless_edges": bool(collapse_mutationless_edges),
        "missing_data_classifier": "assign_missing_average",
        "missing_state_indicator": MISSING_STATE,
        "root_semantics": "solver_generated",
    }
    tree_artifact = CassiopeiaTree(
        canonical_tree,
        {
            "lineage_id": lineage_id,
            "dependency": copy.deepcopy(dependency),
            "upstream_characters_fingerprint": upstream_fingerprint,
            "upstream_character_matrix_sha256": matrix_before,
            "upstream_qc_status": qc_record.get("status") if isinstance(qc_record, Mapping) else None,
            "solver": {
                "class": "VanillaGreedySolver",
                **parameters,
            },
        },
        _token=_ARTIFACT_TOKEN,
        _owned=_worker_owned,
    )
    results = (
        f"VanillaGreedy reconstructed lineage {lineage_id!r} with {len(content['leaves'])} leaves, "
        f"{content['node_count']} total nodes, {content['edge_count']} edges, and maximum topological depth "
        f"{content['max_depth']}."
    )
    summary = _strict_summary(
        node_id="OpenBioSingleCellReconstructCassiopeiaTree",
        methods=(
            "One private computable character matrix was solved with Cassiopeia VanillaGreedy, "
            "the documented assign_missing_average callable, empirical mutation priors, and the declared prior "
            "transformation. The solver-generated root was retained. Internal nodes were relabeled from their exact "
            "descendant leaf sets, and the rooted arborescence, leaf/character axes, priors, and branch lengths were "
            "validated before creating a tamper-evident tree artifact."
        ),
        results=results,
        key_results={
            "lineage_id": lineage_id,
            "leaf_count": len(content["leaves"]),
            "character_count": len(content["characters"]),
            "node_count": content["node_count"],
            "edge_count": content["edge_count"],
            "maximum_depth": content["max_depth"],
            "root": content["root"],
            "root_semantics": "solver_generated",
            "topology_sha256": content["topology_sha256"],
            "tree_fingerprint": tree_artifact.fingerprint,
            "upstream_characters_fingerprint": upstream_fingerprint,
        },
        parameters=parameters,
        warnings_list=messages,
        limitations=(
            "VanillaGreedy is a heuristic reconstruction; a topology is not a confidence interval or proof of ancestry.",
            "Cells within this one lineage are not independent biological replicates.",
            "Solver tie resolution is encounter-order dependent; canonical input order and descendant-set node IDs make this run reproducible.",
        ),
        references=(_cassiopeia_reference(), _software_reference()),
        openbio_version=openbio_version,
    )
    return tree_artifact, summary, _standalone_source()


def _preorder_nodes(structure: Mapping[str, Any]) -> list[str]:
    order: list[str] = []

    def visit(node: str) -> None:
        order.append(node)
        for child in structure["children"][node]:
            visit(child)

    visit(str(structure["root"]))
    return order


def _descendant_leaf_counts(structure: Mapping[str, Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for node in sorted(structure["nodes"], key=lambda value: structure["depths"][value], reverse=True):
        children = structure["children"][node]
        counts[node] = 1 if not children else sum(counts[child] for child in children)
    return counts


def _benjamini_hochberg(values: Sequence[tuple[str, float]]) -> dict[str, float]:
    ordered = sorted(enumerate(values), key=lambda item: (item[1][1], item[0]))
    adjusted: dict[str, float] = {}
    running = 1.0
    family_size = len(ordered)
    for reverse_index in range(family_size - 1, -1, -1):
        original_index, (node, probability) = ordered[reverse_index]
        del original_index
        rank = reverse_index + 1
        running = min(running, probability * family_size / rank)
        adjusted[node] = min(1.0, float(running))
    return adjusted


def compute_cassiopeia_expansions(
    tree: CassiopeiaTree,
    *,
    minimum_clade_size: int = 10,
    minimum_depth: int = 1,
    fdr_threshold: float = 0.05,
    openbio_version: str = "standalone",
    _worker_owned: bool = False,
) -> tuple[Any, dict[str, Any], str]:
    _, pd = _numpy_pandas()
    artifact = _require_tree(tree)
    minimum_clade_size = int(minimum_clade_size)
    minimum_depth = int(minimum_depth)
    if minimum_clade_size < 1:
        raise ValueError("minimum_clade_size must be positive.")
    if minimum_depth < 0:
        raise ValueError("minimum_depth cannot be negative.")
    fdr_threshold = _validate_fraction(fdr_threshold, "fdr_threshold")
    cassiopeia = _require_cassiopeia()
    source_fingerprint = artifact.fingerprint
    backend_input = artifact.copy_tree(_owned=_worker_owned)
    input_content = _tree_content(backend_input)
    graph = backend_input.get_tree_topology()
    structure = _graph_structure(graph)
    leaf_counts = _descendant_leaf_counts(structure)
    order = _preorder_nodes(structure)
    eligibility: dict[str, bool] = {}
    reasons: dict[str, str | None] = {}
    for node in order:
        node_reasons: list[str] = []
        if node == structure["root"]:
            node_reasons.append("root_not_tested")
        if leaf_counts[node] < minimum_clade_size:
            node_reasons.append("below_minimum_clade_size")
        if structure["depths"][node] < minimum_depth:
            node_reasons.append("below_minimum_depth")
        eligibility[node] = not node_reasons
        reasons[node] = ";".join(node_reasons) if node_reasons else None
    with _preserved_rng(0), _captured_backend() as (stream, caught):
        backend_result = cassiopeia.tl.compute_expansion_pvalues(
            backend_input,
            min_clade_size=minimum_clade_size,
            min_depth=minimum_depth,
            copy=True,
        )
    messages = _captured_messages(stream, caught)
    if minimum_clade_size == 1:
        messages.append(
            "minimum_clade_size=1 includes singleton terminal clades; this expert-selected family is valid but "
            "usually uninformative because singleton expansion probabilities equal one."
        )
    if backend_result is None:
        raise RuntimeError("Cassiopeia expansion API returned None despite copy=True.")
    if _tree_content(backend_input) != input_content:
        raise RuntimeError("Cassiopeia expansion API mutated its source tree despite copy=True.")
    if _tree_content(backend_result) != input_content:
        raise RuntimeError("Cassiopeia expansion API changed tree topology, axes, characters, or priors.")
    raw_values: list[tuple[str, float]] = []
    pvalues: dict[str, float] = {}
    for node in order:
        try:
            backend_probability = float(backend_result.get_attribute(node, "expansion_pvalue"))
        except (AttributeError, KeyError, TypeError, ValueError) as error:
            raise RuntimeError(f"Cassiopeia expansion API omitted node {node!r} probability.") from error
        if not math.isfinite(backend_probability) or not 0.0 <= backend_probability <= 1.0:
            raise RuntimeError("Cassiopeia expansion API returned a non-finite/out-of-range probability.")
        if not eligibility[node]:
            if not math.isclose(backend_probability, 1.0, rel_tol=0.0, abs_tol=1e-15):
                raise RuntimeError("Cassiopeia expansion API assigned a probability to an ineligible node.")
            continue
        parent_values = structure["parents"][node]
        if len(parent_values) != 1:
            raise RuntimeError("Eligible expansion clade does not have exactly one parent.")
        parent = parent_values[0]
        population_size = leaf_counts[parent]
        clade_size = leaf_counts[node]
        child_count = len(structure["children"][parent])
        denominator = math.comb(population_size - 1, child_count - 1)
        if denominator <= 0:
            raise RuntimeError("Expansion null probability denominator is zero for an eligible clade.")
        expected = math.comb(population_size - clade_size, child_count - 1) / denominator
        if not math.isclose(backend_probability, expected, rel_tol=1e-12, abs_tol=1e-15):
            raise RuntimeError(
                f"Cassiopeia expansion probability for {node!r} disagrees with the audited coalescent formula."
            )
        pvalues[node] = expected
        raw_values.append((node, expected))
    adjusted = _benjamini_hochberg(raw_values)
    records: list[dict[str, Any]] = []
    for node in order:
        parent_values = structure["parents"][node]
        probability = pvalues.get(node)
        q_value = adjusted.get(node)
        records.append(
            {
                "node_id": node,
                "parent_id": parent_values[0] if parent_values else None,
                "depth": int(structure["depths"][node]),
                "leaf_count": int(leaf_counts[node]),
                "child_count": len(structure["children"][node]),
                "eligible": bool(eligibility[node]),
                "exclusion_reason": reasons[node],
                "p_value": probability,
                "q_value_bh": q_value,
                "significant_fdr": bool(q_value is not None and q_value <= fdr_threshold),
            }
        )
    table = pd.DataFrame.from_records(records)
    significant = [record for record in records if record["significant_fdr"]]
    eligible_records = [record for record in records if record["eligible"]]
    strongest = min(eligible_records, key=lambda record: (record["q_value_bh"], record["node_id"])) if eligible_records else None
    artifact._validate()
    if artifact.fingerprint != source_fingerprint:
        raise RuntimeError("Source Cassiopeia tree artifact changed during expansion testing.")
    parameters = {
        "minimum_clade_size": minimum_clade_size,
        "minimum_depth": minimum_depth,
        "fdr_threshold": fdr_threshold,
        "multiple_testing": "Benjamini-Hochberg across all eligible clades in this tree",
        "inference_unit": "nested clade within one inferred lineage tree",
    }
    results = (
        f"Tested {len(eligible_records)} eligible nested clade(s) in lineage {artifact.lineage_id!r}; "
        f"{len(significant)} clade(s) met Benjamini-Hochberg FDR ≤ {fdr_threshold:g}."
    )
    summary = _strict_summary(
        node_id="OpenBioSingleCellCassiopeiaExpansionTest",
        methods=(
            "Every topology node was enumerated before calling Cassiopeia compute_expansion_pvalues(copy=True). "
            "Eligibility used the declared integer descendant-leaf count and topological depth. Each eligible "
            "one-sided neutral-coalescent probability was independently recomputed as C(n-b,k-1)/C(n-1,k-1), "
            "then Benjamini-Hochberg correction was applied once to the complete eligible within-tree family."
        ),
        results=results,
        key_results={
            "lineage_id": artifact.lineage_id,
            "tree_fingerprint": source_fingerprint,
            "topology_sha256": artifact.topology_fingerprint,
            "topology_nodes_reported": len(records),
            "eligible_hypotheses": len(eligible_records),
            "fdr_significant_clades": len(significant),
            "strongest_clade": strongest,
        },
        parameters=parameters,
        warnings_list=(
            *messages,
            *( ["No clade met the declared FDR threshold."] if not significant else [] ),
        ),
        limitations=(
            "Clades are nested and arise from one inferred tree; they are not independent biological samples.",
            "The simple neutral-coalescent null does not model all sampling, editing, or fitness processes.",
            "Statistical significance is conditional on the reconstructed topology and declared hypothesis family.",
        ),
        references=(
            _yang_reference(),
            {
                "citation": "Benjamini Y, Hochberg Y. Controlling the false discovery rate. JRSS B. 1995.",
                "doi": "10.1111/j.2517-6161.1995.tb02031.x",
                "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
                "kind": "method",
            },
            _software_reference(),
        ),
        openbio_version=openbio_version,
    )
    return table, summary, _standalone_source()


def _value_nbytes(value: Any, seen: set[int]) -> int:
    if value is None or id(value) in seen:
        return 0
    seen.add(id(value))
    if isinstance(value, (str, bytes, bytearray)):
        return len(value)
    if hasattr(value, "data") and hasattr(value, "indices") and hasattr(value, "indptr"):
        return sum(_value_nbytes(getattr(value, name), seen) for name in ("data", "indices", "indptr"))
    if hasattr(value, "memory_usage"):
        try:
            usage = value.memory_usage(index=True, deep=True)
            return int(usage.sum() if hasattr(usage, "sum") else usage)
        except (TypeError, ValueError):
            pass
    if hasattr(value, "nbytes"):
        try:
            return int(value.nbytes)
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return sum(_value_nbytes(key, seen) + _value_nbytes(item, seen) for key, item in value.items())
    if isinstance(value, Sequence):
        return sum(_value_nbytes(item, seen) for item in value)
    return 64


def _adata_working_bytes(adata: Any, topology_nodes: int, state_count: int) -> int:
    seen: set[int] = set()
    total = _value_nbytes(getattr(adata, "X", None), seen)
    for attribute in ("obs", "var", "layers", "obsm", "varm", "obsp", "varp", "uns"):
        total += _value_nbytes(getattr(adata, attribute, None), seen)
    raw = getattr(adata, "raw", None)
    if raw is not None:
        total += _value_nbytes(getattr(raw, "X", None), seen)
    dynamic_program = int(topology_nodes) * max(1, int(state_count)) * 24
    topology_overhead = int(topology_nodes) * 1024
    return total * 2 + dynamic_program + topology_overhead


def _adata_snapshot(adata: Any, annotation_key: str) -> str:
    annotation = adata.obs[annotation_key]
    return _hash_payload(
        {
            "obs_names": [str(value) for value in adata.obs_names],
            "var_names": [str(value) for value in adata.var_names],
            "obs_columns": [str(value) for value in adata.obs.columns],
            "annotation_values": [_json_scalar(value) for value in annotation.tolist()],
            "annotation_categories": [str(value) for value in annotation.cat.categories.tolist()],
            "plasticity_provenance": _json_value(adata.uns.get("openbio_cassiopeia_plasticity")),
        }
    )


def _annotation_provenance(
    adata: Any, annotation_key: str, annotation_status: str, analysis_mode: str
) -> tuple[dict[str, Any] | None, list[str]]:
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError("annotation_status must be unknown, provisional, or curated.")
    if analysis_mode not in {"exploratory", "report_grade"}:
        raise ValueError("analysis_mode must be exploratory or report_grade.")
    messages: list[str] = []
    openbio = adata.uns.get("openbio_singlecell")
    provenance = None
    if openbio is not None:
        if not isinstance(openbio, Mapping):
            messages.append("Malformed OpenBio annotation registry was ignored; annotation provenance is unverified.")
        else:
            annotations = openbio.get("annotations", {})
            if annotations is not None and not isinstance(annotations, Mapping):
                messages.append(
                    "Malformed OpenBio annotation registry was ignored; annotation provenance is unverified."
                )
            else:
                candidate = annotations.get(annotation_key) if isinstance(annotations, Mapping) else None
                if candidate is not None and not isinstance(candidate, Mapping):
                    messages.append(
                        "Malformed annotation provenance was ignored; the caller-declared annotation status was used."
                    )
                elif isinstance(candidate, Mapping):
                    provenance = candidate
    if provenance is not None:
        registered = provenance.get("annotation_status", provenance.get("status"))
        if registered not in {"unknown", "provisional", "curated"}:
            messages.append("Registered annotation status is invalid and was treated as unverified provenance.")
        elif registered != annotation_status:
            messages.append(
                f"Caller-declared annotation status {annotation_status!r} disagrees with registered status "
                f"{registered!r}; both are disclosed and expert interpretation is required."
            )
    if analysis_mode == "report_grade":
        registered = provenance.get("annotation_status", provenance.get("status")) if provenance else None
        if annotation_status != "curated" or registered != "curated":
            messages.append(
                "Report-grade mode was requested without matching registered Curated annotation provenance; "
                "execution was allowed for expert use, but this provenance limitation must remain in the report."
            )
    if analysis_mode == "exploratory":
        messages.append(
            f"Exploratory analysis used annotation {annotation_key!r} with caller-declared status {annotation_status!r}."
        )
    if annotation_status == "curated":
        messages.append("Curated annotation status remains a caller/workflow provenance assertion.")
    return copy.deepcopy(dict(provenance)) if provenance is not None else None, messages


def _preprocess_plasticity_graph(
    graph: Any, leaf_states: Mapping[str, str], retained_states: set[str]
) -> tuple[Any, dict[str, Any], int]:
    structure = _graph_structure(graph)
    retained_leaves = sorted(
        leaf for leaf in structure["leaves"] if leaf in leaf_states and leaf_states[leaf] in retained_states
    )
    if len(retained_leaves) < 2:
        raise ValueError("Plasticity requires at least two retained tree leaves.")
    keep_nodes: set[str] = set()
    for leaf in retained_leaves:
        node = leaf
        while True:
            keep_nodes.add(node)
            parents = structure["parents"][node]
            if not parents:
                break
            node = parents[0]
    processed = graph.subgraph(keep_nodes).copy()
    original_root = structure["root"]
    while True:
        current = _graph_structure(processed)
        collapsible = [
            node
            for node in current["nodes"]
            if node != original_root and len(current["children"][node]) == 1
        ]
        if not collapsible:
            break
        node = sorted(collapsible, key=lambda value: current["depths"][value], reverse=True)[0]
        parent = current["parents"][node][0]
        child = current["children"][node][0]
        first_length = float(processed[parent][node].get("length", 1.0))
        second_length = float(processed[node][child].get("length", 1.0))
        processed.remove_node(node)
        processed.add_edge(parent, child, length=first_length + second_length)
    current = _graph_structure(processed)
    exhausted_nodes = 0
    for node in _preorder_nodes(current):
        children = list(current["children"][node])
        if len(children) < 3 or not all(not current["children"][child] for child in children):
            continue
        states = sorted({leaf_states[child] for child in children})
        if len(states) < 2:
            continue
        exhausted_nodes += 1
        for state in states:
            digest = hashlib.sha256(_canonical_json([node, state]).encode("utf-8")).hexdigest()
            group_node = f"{_GROUP_PREFIX}{digest}"
            if group_node in processed:
                raise RuntimeError("Deterministic exhausted-lineage node identity collision.")
            processed.add_edge(node, group_node, length=1.0)
            for child in sorted(child for child in children if leaf_states[child] == state):
                processed.remove_edge(node, child)
                processed.add_edge(group_node, child, length=0.0)
    final_structure = _graph_structure(processed)
    if set(final_structure["leaves"]) != set(retained_leaves):
        raise RuntimeError("Plasticity preprocessing changed retained leaf identity.")
    if final_structure["root"] != original_root:
        raise RuntimeError("Plasticity preprocessing changed the solver-generated root.")
    return processed, final_structure, exhausted_nodes


def _sankoff_subtree_scores(
    structure: Mapping[str, Any], leaf_states: Mapping[str, str], states: Sequence[str]
) -> tuple[dict[str, float], dict[str, int], int]:
    state_to_position = {state: position for position, state in enumerate(states)}
    costs: dict[str, list[float]] = {}
    edge_counts: dict[str, int] = {}
    nodes = sorted(structure["nodes"], key=lambda value: structure["depths"][value], reverse=True)
    for node in nodes:
        children = structure["children"][node]
        if not children:
            if node not in leaf_states or leaf_states[node] not in state_to_position:
                raise ValueError("A retained topology leaf lacks a retained categorical state.")
            own = state_to_position[leaf_states[node]]
            costs[node] = [0.0 if position == own else math.inf for position in range(len(states))]
            edge_counts[node] = 0
            continue
        node_costs = [0.0] * len(states)
        for child in children:
            child_costs = costs[child]
            child_minimum = min(child_costs)
            for position in range(len(states)):
                node_costs[position] += min(child_costs[position], child_minimum + 1.0)
        costs[node] = node_costs
        edge_counts[node] = sum(edge_counts[child] + 1 for child in children)
    subtree_scores: dict[str, float] = {}
    for node in structure["nodes"]:
        if edge_counts[node] == 0:
            continue
        score = min(costs[node]) / edge_counts[node]
        if not math.isfinite(score) or not 0.0 <= score <= 1.0:
            raise RuntimeError("EffectivePlasticity dynamic program produced an invalid score.")
        subtree_scores[node] = float(score)
    root_score = min(costs[structure["root"]])
    if not math.isfinite(root_score) or root_score != int(root_score):
        raise RuntimeError("Global small-parsimony score is not a finite integer.")
    return subtree_scores, edge_counts, int(root_score)


def _single_cell_effective_plasticity(
    structure: Mapping[str, Any], subtree_scores: Mapping[str, float]
) -> tuple[dict[str, float], dict[str, int]]:
    scores: dict[str, float] = {}
    path_counts: dict[str, int] = {}
    for leaf in structure["leaves"]:
        values: list[float] = []
        node = leaf
        while structure["parents"][node]:
            node = structure["parents"][node][0]
            if node in subtree_scores:
                values.append(float(subtree_scores[node]))
        if not values:
            raise RuntimeError("A retained leaf has no non-leaf ancestor plasticity score.")
        score = sum(values) / len(values)
        if not math.isfinite(score) or not 0.0 <= score <= 1.0:
            raise RuntimeError("Single-cell EffectivePlasticity is outside [0, 1].")
        scores[leaf] = float(score)
        path_counts[leaf] = len(values)
    return scores, path_counts


def add_cassiopeia_plasticity(
    adata: Any,
    tree: CassiopeiaTree,
    *,
    annotation_key: str = "cell_type",
    annotation_status: str = "unknown",
    analysis_mode: str = "exploratory",
    minimum_state_fraction: float = 0.025,
    output_key: str = "sc_effective_plasticity",
    overwrite_existing: bool = False,
    max_working_gib: float = 4.0,
    openbio_version: str = "standalone",
    _worker_owned: bool = False,
) -> tuple[Any, Any, dict[str, Any], str]:
    np, pd = _numpy_pandas()
    artifact = _require_tree(tree)
    if not hasattr(adata, "obs") or not hasattr(adata, "obs_names") or not hasattr(adata, "copy"):
        raise TypeError("Plasticity requires an in-memory AnnData-like object.")
    if not isinstance(annotation_key, str) or not annotation_key or annotation_key != annotation_key.strip():
        raise ValueError("annotation_key must be an exact non-empty string.")
    if annotation_key not in adata.obs:
        raise ValueError(f"AnnData obs is missing annotation column {annotation_key!r}.")
    if not isinstance(output_key, str) or not output_key or output_key != output_key.strip():
        raise ValueError("output_key must be an exact non-empty string.")
    if output_key == annotation_key:
        raise ValueError("output_key cannot overwrite the source annotation column.")
    status_key = f"{output_key}_status"
    minimum_state_fraction = _validate_fraction(minimum_state_fraction, "minimum_state_fraction")
    if not math.isfinite(float(max_working_gib)) or float(max_working_gib) <= 0:
        raise ValueError("max_working_gib must be positive and finite.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("AnnData observation and variable axes must be unique.")
    obs_names = list(adata.obs_names)
    var_names = list(adata.var_names)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in obs_names + var_names):
        raise ValueError("AnnData axes must be exact non-empty strings.")
    annotation = adata.obs[annotation_key]
    if not isinstance(annotation.dtype, pd.CategoricalDtype):
        raise ValueError("Plasticity annotation must use pandas categorical dtype.")
    categories = annotation.cat.categories.tolist()
    if any(not isinstance(value, str) or not value or value != value.strip() for value in categories):
        raise ValueError("Plasticity annotation categories must be exact non-empty strings.")
    _, provenance_messages = _annotation_provenance(
        adata, annotation_key, annotation_status, analysis_mode
    )
    provenance_registry = adata.uns.get("openbio_cassiopeia_plasticity", {})
    if provenance_registry is not None and not isinstance(provenance_registry, Mapping):
        raise ValueError("Existing Cassiopeia plasticity provenance must be a mapping.")
    ownership = provenance_registry.get(output_key) if isinstance(provenance_registry, Mapping) else None
    collisions = [column for column in (output_key, status_key) if column in adata.obs]
    if collisions and not overwrite_existing:
        raise ValueError("Plasticity output columns already exist and overwrite_existing is disabled.")
    if collisions and (
        not isinstance(ownership, Mapping)
        or ownership.get("producer") != PRODUCER
        or ownership.get("node_id") != "OpenBioSingleCellCassiopeiaPlasticity"
        or ownership.get("output_key") != output_key
    ):
        raise ValueError("Existing output columns are not owned by this audited plasticity node.")

    source_fingerprint = artifact.fingerprint
    backend_tree = artifact.copy_tree(_owned=_worker_owned)
    graph = backend_tree.get_tree_topology()
    structure = _graph_structure(graph)
    tree_leaves = set(structure["leaves"])
    missing_leaves = sorted(tree_leaves - set(obs_names))
    missing_annotation_leaves = sorted(
        leaf for leaf in tree_leaves & set(obs_names) if bool(pd.isna(annotation.loc[leaf]))
    )
    leaf_states = {
        leaf: str(annotation.loc[leaf])
        for leaf in structure["leaves"]
        if leaf in set(obs_names) and leaf not in set(missing_annotation_leaves)
    }
    if len(leaf_states) < 2:
        raise ValueError("Plasticity requires at least two tree leaves with available categorical annotations.")
    state_counts = Counter(leaf_states.values())
    tree_leaf_count = len(structure["leaves"])
    retained_states = {
        state for state, count in state_counts.items() if count / tree_leaf_count >= minimum_state_fraction
    }
    retained_leaf_count = sum(count for state, count in state_counts.items() if state in retained_states)
    if retained_leaf_count < 2:
        raise ValueError(
            "Plasticity requires at least two annotated tree leaves after the declared frequency filter."
        )
    if missing_leaves:
        provenance_messages.append(
            f"AnnData omitted {len(missing_leaves)} tree leaf/leaves; they were excluded from the descriptive "
            "plasticity estimand rather than treated as zero."
        )
    if missing_annotation_leaves:
        provenance_messages.append(
            f"Excluded {len(missing_annotation_leaves)} tree leaf/leaves with missing annotation state; their "
            "output scores remain null."
        )
    if len(retained_states) == 1:
        provenance_messages.append(
            "Only one annotation state remained after filtering; EffectivePlasticity is computable and equals zero, "
            "but there is no observed state diversity to interpret."
        )
    working_bytes = _adata_working_bytes(adata, len(structure["nodes"]), len(state_counts))
    if working_bytes > float(max_working_gib) * 1024**3:
        raise ValueError(
            f"Conservative plasticity working estimate {working_bytes / 1024**3:.3f} GiB exceeds max_working_gib."
        )
    input_snapshot = _adata_snapshot(adata, annotation_key)
    processed_graph, processed_structure, exhausted_nodes = _preprocess_plasticity_graph(
        graph, leaf_states, retained_states
    )
    retained_state_axis = tuple(sorted(retained_states))
    subtree_scores, edge_counts, global_parsimony = _sankoff_subtree_scores(
        processed_structure, leaf_states, retained_state_axis
    )
    cell_scores, path_counts = _single_cell_effective_plasticity(processed_structure, subtree_scores)

    cassiopeia = _require_cassiopeia()
    retained_leaf_axis = list(processed_structure["leaves"])
    cell_meta = pd.DataFrame(
        {
            "_openbio_state": pd.Categorical(
                [leaf_states[leaf] for leaf in retained_leaf_axis], categories=list(retained_state_axis)
            )
        },
        index=retained_leaf_axis,
    )
    backend_for_check = cassiopeia.data.CassiopeiaTree(
        cell_meta=cell_meta.copy(deep=True), tree=processed_graph.copy()
    )
    with _preserved_rng(0), _captured_backend() as (stream, caught):
        backend_parsimony = cassiopeia.tl.score_small_parsimony(
            backend_for_check,
            meta_item="_openbio_state",
            root=backend_for_check.root,
            infer_ancestral_states=True,
            label_key="_openbio_label",
        )
    backend_messages = _captured_messages(stream, caught)
    if isinstance(backend_parsimony, (bool, np.bool_)):
        raise RuntimeError("Cassiopeia small-parsimony backend returned a Boolean score.")
    try:
        backend_numeric = float(backend_parsimony)
    except (TypeError, ValueError, OverflowError) as error:
        raise RuntimeError("Cassiopeia small-parsimony backend returned a non-numeric score.") from error
    if not math.isfinite(backend_numeric) or backend_numeric != global_parsimony:
        raise RuntimeError("Audited dynamic program disagrees with Cassiopeia global small-parsimony score.")
    artifact._validate()
    if artifact.fingerprint != source_fingerprint:
        raise RuntimeError("Source Cassiopeia tree artifact changed during plasticity analysis.")
    if _adata_snapshot(adata, annotation_key) != input_snapshot:
        raise RuntimeError("Plasticity analysis mutated its AnnData input.")

    output = adata if _worker_owned else adata.copy()
    score_values: list[Any] = []
    status_values: list[str] = []
    table_records: list[dict[str, Any]] = []
    for cell_id in obs_names:
        raw_state = annotation.loc[cell_id]
        state = None if bool(pd.isna(raw_state)) else str(raw_state)
        if cell_id not in tree_leaves:
            status = "not_in_tree"
            score = None
            path_count = 0
            lineage_value = None
        elif state is None:
            status = "missing_annotation"
            score = None
            path_count = 0
            lineage_value = artifact.lineage_id
        elif state not in retained_states:
            status = "state_below_minimum_frequency"
            score = None
            path_count = 0
            lineage_value = artifact.lineage_id
        else:
            status = "included"
            score = float(cell_scores[cell_id])
            path_count = int(path_counts[cell_id])
            lineage_value = artifact.lineage_id
        score_values.append(pd.NA if score is None else score)
        status_values.append(status)
        table_records.append(
            {
                "cell_id": cell_id,
                "lineage_id": lineage_value,
                "annotation_state": state,
                "state_count_in_tree": int(state_counts.get(state, 0)),
                "state_fraction_in_tree": float(state_counts.get(state, 0) / tree_leaf_count),
                "status": status,
                "path_subtree_count": path_count,
                "sc_effective_plasticity": score,
            }
        )
    output.obs[output_key] = pd.Series(score_values, index=obs_names, dtype="Float64")
    status_categories = ["included", "missing_annotation", "state_below_minimum_frequency", "not_in_tree"]
    output.obs[status_key] = pd.Categorical(status_values, categories=status_categories)
    output_registry = output.uns.get("openbio_cassiopeia_plasticity", {})
    if not isinstance(output_registry, Mapping):
        raise ValueError("Output Cassiopeia plasticity provenance must be a mapping.")
    output_registry = copy.deepcopy(dict(output_registry))
    annotation_provenance = None
    openbio = output.uns.get("openbio_singlecell")
    if isinstance(openbio, Mapping) and isinstance(openbio.get("annotations"), Mapping):
        annotation_provenance = openbio["annotations"].get(annotation_key)
    provenance_record = {
        "producer": PRODUCER,
        "node_id": "OpenBioSingleCellCassiopeiaPlasticity",
        "output_key": output_key,
        "status_key": status_key,
        "tree_fingerprint": source_fingerprint,
        "topology_sha256": artifact.topology_fingerprint,
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "analysis_mode": analysis_mode,
        "annotation_provenance_sha256": _hash_payload(annotation_provenance),
        "minimum_state_fraction": minimum_state_fraction,
        "retained_states": list(retained_state_axis),
        "global_small_parsimony": global_parsimony,
    }
    output_registry[output_key] = provenance_record
    output.uns["openbio_cassiopeia_plasticity"] = output_registry
    if list(output.obs_names) != obs_names or list(output.var_names) != var_names:
        raise RuntimeError("Plasticity output changed AnnData observation or variable axes.")
    if not output.obs[annotation_key].equals(annotation):
        raise RuntimeError("Plasticity output changed the source annotation.")
    table = pd.DataFrame.from_records(table_records)
    included_values = [float(value) for value in cell_scores.values()]
    ordered_values = sorted(included_values)
    mean_score = float(sum(ordered_values) / len(ordered_values))
    median_score = float(np.median(np.asarray(ordered_values, dtype=float)))
    parameters = {
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "analysis_mode": analysis_mode,
        "minimum_state_fraction": minimum_state_fraction,
        "output_key": output_key,
        "overwrite_existing": bool(overwrite_existing),
        "max_working_gib": float(max_working_gib),
        "state_frequency_rule": "retain frequency >= threshold among original tree leaves",
        "subtree_denominator": "directed edge count",
        "single_cell_aggregation": "mean over non-leaf root-to-leaf ancestors",
    }
    results = (
        f"EffectivePlasticity was defined for {len(included_values)} retained cell(s) in lineage "
        f"{artifact.lineage_id!r} across {len(retained_state_axis)} retained states; the mean single-cell score was "
        f"{mean_score:.6g} and global Fitch-Hartigan parsimony was {global_parsimony}."
    )
    summary = _strict_summary(
        node_id="OpenBioSingleCellCassiopeiaPlasticity",
        methods=(
            "Annotation states with frequency at least the declared threshold among original tree leaves were retained. "
            "Rare-state leaves were pruned, non-root unifurcations collapsed, and terminal exhausted polytomies split "
            "into deterministic state-group nodes following the released study procedure. A bounded unit-cost Sankoff "
            "dynamic program computed Fitch-Hartigan minima for every rooted subtree; EffectivePlasticity divided each "
            "minimum by subtree edge count, and each retained cell received the mean over its non-leaf ancestors. The "
            "global minimum was independently checked with Cassiopeia score_small_parsimony under restored RNG state."
        ),
        results=results,
        key_results={
            "lineage_id": artifact.lineage_id,
            "tree_fingerprint": source_fingerprint,
            "tree_leaves": tree_leaf_count,
            "retained_cells": len(included_values),
            "excluded_rare_state_cells": sum(value == "state_below_minimum_frequency" for value in status_values),
            "tree_leaves_missing_from_anndata": len(missing_leaves),
            "tree_leaves_with_missing_annotation": len(missing_annotation_leaves),
            "ann_data_cells_not_in_tree": sum(value == "not_in_tree" for value in status_values),
            "original_state_count": len(state_counts),
            "retained_states": list(retained_state_axis),
            "processed_node_count": len(processed_structure["nodes"]),
            "processed_edge_count": len(processed_structure["edges"]),
            "exhausted_polytomies_split": exhausted_nodes,
            "global_small_parsimony": global_parsimony,
            "mean_sc_effective_plasticity": mean_score,
            "median_sc_effective_plasticity": median_score,
            "minimum_sc_effective_plasticity": min(ordered_values),
            "maximum_sc_effective_plasticity": max(ordered_values),
        },
        parameters=parameters,
        warnings_list=(*provenance_messages, *backend_messages),
        limitations=(
            "EffectivePlasticity is descriptive and conditional on one inferred tree and one supplied annotation.",
            "Cells and nested subtrees are not independent biological samples; no condition-level inference is performed.",
            "Rare-state filtering and terminal-polytomy preprocessing can materially change scores and are disclosed parameters.",
        ),
        references=(_yang_reference(), _cassiopeia_reference(), _software_reference()),
        openbio_version=openbio_version,
    )
    return output, table, summary, _standalone_source()
__all__ = [
    "ALLELE_TABLE_EXTENSIONS",
    "CassiopeiaCharacters",
    "CassiopeiaTree",
    "add_cassiopeia_plasticity",
    "compute_cassiopeia_expansions",
    "prepare_cassiopeia_characters",
    "reconstruct_cassiopeia_tree",
]
