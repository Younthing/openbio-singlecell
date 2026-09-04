from __future__ import annotations

import copy
import hashlib
import json
import math
import numbers
import re
from collections.abc import Mapping, Sequence
from typing import Any

MILO_RESULT_ARTIFACT_TYPE = "OPENBIO_MILO_RESULT"
MILO_RESULT_ARTIFACT_SCHEMA_VERSION = 1
MILO_RESULT_PRODUCER_NODE_ID = "OpenBioSingleCellMiloDifferentialAbundance"
MILO_RESULT_PRODUCER_SCHEMA = "openbio-singlecell/milo-result/v1"
MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS = 10_000_000
MILO_MAX_OVERLAP_WORKING_BYTES = 2 * 1024**3

MILO_RESULT_COLUMNS = (
    "neighborhood_id",
    "index_cell",
    "neighborhood_size",
    "kth_distance",
    "reference_condition",
    "comparison_condition",
    "log2_fold_change",
    "log_counts_per_million",
    "quasi_likelihood_f",
    "p_value",
    "p_adjusted_bh",
    "spatial_fdr",
    "majority_annotation",
    "majority_annotation_fraction",
    "reported_annotation",
    "is_mixed",
)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_TEXT_COLUMNS = {
    "neighborhood_id",
    "index_cell",
    "reference_condition",
    "comparison_condition",
    "majority_annotation",
    "reported_annotation",
}
_UNIT_COLUMNS = {"p_value", "p_adjusted_bh", "spatial_fdr", "majority_annotation_fraction"}


def _json_sha256(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _axis(values: Sequence[Any], *, label: str) -> tuple[list[str], str]:
    result: list[str] = []
    for position, value in enumerate(values):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"Milo {label} at position {position} must be a canonical nonblank string.")
        result.append(value)
    if not result or len(result) != len(set(result)):
        raise ValueError(f"Milo {label} must be nonempty and unique.")
    return result, _json_sha256(result)


def _csr(matrix: Any, *, label: str, rows: int, columns: int, binary: bool = False) -> Any:
    import numpy as np
    from scipy import sparse

    if not sparse.issparse(matrix):
        array = np.asarray(matrix)
        if array.ndim != 2:
            raise ValueError(f"Milo {label} must be two-dimensional.")
    result = sparse.csr_matrix(matrix, copy=True)
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    if result.shape != (rows, columns):
        raise ValueError(f"Milo {label} shape does not match its named axes.")
    values = np.asarray(result.data)
    if values.size and (
        not np.issubdtype(values.dtype, np.number)
        or np.issubdtype(values.dtype, np.bool_)
        or not bool(np.isfinite(values).all())
    ):
        raise TypeError(f"Milo {label} must contain finite non-boolean numbers.")
    if values.size and bool((values < 0).any()):
        raise ValueError(f"Milo {label} cannot contain negative values.")
    if binary and values.size and not bool(np.equal(values, 1).all()):
        raise ValueError("Milo neighborhood membership must be binary.")
    return result.astype(np.int64) if binary else result


def _matrix_sha256(matrix: Any, *, label: str) -> str:
    import numpy as np

    digest = hashlib.sha256(f"openbio-singlecell/{label}/v1\0".encode("ascii"))
    digest.update(json.dumps(list(matrix.shape), separators=(",", ":")).encode("ascii"))
    for array, dtype in ((matrix.indptr, "<i8"), (matrix.indices, "<i8"), (matrix.data, "<f8")):
        values = np.ascontiguousarray(array, dtype=dtype)
        digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _membership_overlap_graph(membership: Any) -> tuple[Any, dict[str, int]]:
    import numpy as np

    row_memberships = np.diff(membership.indptr)
    pair_contributions = sum(int(value) * (int(value) - 1) // 2 for value in row_memberships)
    membership_bytes = int(membership.data.nbytes + membership.indices.nbytes + membership.indptr.nbytes)
    estimated_working_bytes = membership_bytes * 3 + pair_contributions * 64
    if pair_contributions > MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS:
        raise MemoryError(
            "Milo overlap-graph pair contributions exceed the fixed artifact limit: "
            f"{pair_contributions:,} > {MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS:,}."
        )
    if estimated_working_bytes > MILO_MAX_OVERLAP_WORKING_BYTES:
        raise MemoryError(
            "Milo overlap-graph estimated working memory exceeds the fixed artifact limit: "
            f"{estimated_working_bytes / 1024**3:.3f} GiB > "
            f"{MILO_MAX_OVERLAP_WORKING_BYTES / 1024**3:.3f} GiB."
        )
    graph = (membership.T @ membership).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    graph.sort_indices()
    return graph, {
        "pair_contributions": pair_contributions,
        "estimated_working_bytes": estimated_working_bytes,
        "maximum_pair_contributions": MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS,
        "maximum_working_bytes": MILO_MAX_OVERLAP_WORKING_BYTES,
    }


def _coordinates(value: Any, *, neighborhoods: int) -> Any:
    import numpy as np

    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError("Milo representative coordinates must be numeric.") from error
    if result.shape != (neighborhoods, 2):
        raise ValueError("Milo representative coordinates must have one two-dimensional row per neighborhood.")
    if not bool(np.isfinite(result).all()):
        raise ValueError("Milo representative coordinates must be finite.")
    return np.array(result, dtype=float, copy=True, order="C")


def _coordinate_sha256(value: Any) -> str:
    import numpy as np

    digest = hashlib.sha256(b"openbio-singlecell/milo-representative-coordinates/v1\0")
    digest.update(json.dumps(list(value.shape), separators=(",", ":")).encode("ascii"))
    digest.update(np.ascontiguousarray(value, dtype="<f8").tobytes(order="C"))
    return digest.hexdigest()


def _table_identity(table: Any) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    if not isinstance(table, pd.DataFrame):
        raise TypeError("Milo result table must be a pandas DataFrame.")
    if tuple(table.columns) != MILO_RESULT_COLUMNS:
        raise ValueError("Milo result table columns do not match the current canonical schema.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError("Milo result table must use a zero-based RangeIndex.")
    if table.empty:
        raise ValueError("Milo result table cannot be empty.")

    digest = hashlib.sha256(b"openbio-singlecell/milo-result-table/v1\0")
    for column in MILO_RESULT_COLUMNS:
        digest.update(column.encode("utf-8") + b"\0")
        if column in _TEXT_COLUMNS:
            for value in table[column].tolist():
                if not isinstance(value, str) or not value or value != value.strip():
                    raise ValueError(f"Milo result column {column!r} contains a noncanonical string.")
                encoded = value.encode("utf-8")
                digest.update(len(encoded).to_bytes(8, "little") + encoded)
        elif column == "is_mixed":
            if not pd.api.types.is_bool_dtype(table[column].dtype):
                raise TypeError("Milo result is_mixed must use a boolean dtype.")
            digest.update(np.ascontiguousarray(table[column].to_numpy(dtype=np.uint8), dtype="u1").tobytes())
        elif column == "neighborhood_size":
            if not pd.api.types.is_integer_dtype(table[column].dtype):
                raise TypeError("Milo result neighborhood_size must use an integer dtype.")
            values = table[column].to_numpy(dtype=np.int64)
            if bool((values <= 0).any()):
                raise ValueError("Milo result neighborhood sizes must be positive.")
            digest.update(np.ascontiguousarray(values, dtype="<i8").tobytes())
        else:
            if not pd.api.types.is_numeric_dtype(table[column].dtype) or pd.api.types.is_bool_dtype(
                table[column].dtype
            ):
                raise TypeError(f"Milo result column {column!r} must be numeric and non-boolean.")
            values = table[column].to_numpy(dtype=float)
            if not bool(np.isfinite(values).all()):
                raise ValueError(f"Milo result column {column!r} must be finite.")
            if column in _UNIT_COLUMNS and bool(((values < 0.0) | (values > 1.0)).any()):
                raise ValueError(f"Milo result column {column!r} must lie in [0, 1].")
            if column == "kth_distance" and bool((values <= 0.0).any()):
                raise ValueError("Milo result kth distances must be positive.")
            if column in {"log_counts_per_million", "quasi_likelihood_f"} and bool((values < 0.0).any()):
                raise ValueError(f"Milo result column {column!r} cannot be negative.")
            digest.update(np.ascontiguousarray(values, dtype="<f8").tobytes())
    return {"rows": int(len(table)), "columns": list(MILO_RESULT_COLUMNS), "content_sha256": digest.hexdigest()}


class MiloResult:
    """Immutable-by-interface, tamper-evident Milo diagnostic evidence."""

    __slots__ = (
        "_coordinates",
        "_graph",
        "_membership",
        "_metadata",
        "_neighborhood_names",
        "_observation_names",
        "_provenance",
        "_table",
    )
    artifact_type = MILO_RESULT_ARTIFACT_TYPE

    @classmethod
    def _from_owned(cls, **payload: Any) -> MiloResult:
        result = object.__new__(cls)
        for key, value in payload.items():
            object.__setattr__(result, f"_{key}", value)
        return result

    def __setattr__(self, _name: str, _value: Any) -> None:
        raise AttributeError("MiloResult is immutable; create a new validated artifact instead.")

    @property
    def fingerprint(self) -> str:
        return str(self._metadata["artifact_fingerprint_sha256"])

    @property
    def table(self) -> Any:
        return self._table.copy(deep=True)

    @property
    def membership(self) -> Any:
        return self._membership.copy()

    @property
    def graph(self) -> Any:
        return self._graph.copy()

    @property
    def representative_coordinates(self) -> Any:
        return self._coordinates.copy()

    @property
    def observation_names(self) -> list[str]:
        return list(self._observation_names)

    @property
    def neighborhood_names(self) -> list[str]:
        return list(self._neighborhood_names)

    @property
    def provenance(self) -> dict[str, Any]:
        return copy.deepcopy(self._provenance)

    @property
    def metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._metadata)

    def portable(self) -> dict[str, Any]:
        return {
            "artifact_type": MILO_RESULT_ARTIFACT_TYPE,
            "table": self._table.copy(deep=True),
            "membership": self._membership.copy(),
            "graph": self._graph.copy(),
            "coordinates": self._coordinates.copy(),
            "observation_names": list(self._observation_names),
            "neighborhood_names": list(self._neighborhood_names),
            "provenance": copy.deepcopy(self._provenance),
            "metadata": copy.deepcopy(self._metadata),
        }


def build_milo_result(
    *,
    table: Any,
    membership: Any,
    graph: Any,
    representative_coordinates: Any,
    observation_names: Sequence[Any],
    neighborhood_names: Sequence[Any],
    representation_key: str,
    representation_sha256: str,
    condition_key: str,
    annotation_key: str,
    annotation_status: str,
    n_neighbors: int,
    neighborhood_proportion: float,
    mixed_annotation_threshold: float,
    spatial_fdr_threshold: float,
    min_abs_log2_fold_change: float,
    random_seed: int,
) -> MiloResult:
    observations, observation_hash = _axis(observation_names, label="observation axis")
    neighborhoods, neighborhood_hash = _axis(neighborhood_names, label="neighborhood axis")
    if not isinstance(representation_key, str) or not representation_key or representation_key != representation_key.strip():
        raise ValueError("Milo representation key must be a canonical nonblank string.")
    if not isinstance(representation_sha256, str) or _SHA256.fullmatch(representation_sha256) is None:
        raise ValueError("Milo representation fingerprint must be a lowercase SHA-256 digest.")
    for name, value in (("condition_key", condition_key), ("annotation_key", annotation_key)):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"Milo {name} must be a canonical nonblank string.")
    if annotation_status not in {"provisional", "curated"}:
        raise ValueError("Milo annotation status must be 'provisional' or 'curated'.")
    for name, value, minimum in (("n_neighbors", n_neighbors, 2), ("random_seed", random_seed, 0)):
        if isinstance(value, bool) or not isinstance(value, numbers.Integral) or int(value) < minimum:
            raise ValueError(f"Milo {name} must be an integer >= {minimum}.")
    numeric = {
        "neighborhood_proportion": neighborhood_proportion,
        "mixed_annotation_threshold": mixed_annotation_threshold,
        "spatial_fdr_threshold": spatial_fdr_threshold,
        "min_abs_log2_fold_change": min_abs_log2_fold_change,
    }
    if any(isinstance(value, bool) or not isinstance(value, numbers.Real) for value in numeric.values()):
        raise TypeError("Milo diagnostic thresholds must be real non-boolean numbers.")
    numeric = {name: float(value) for name, value in numeric.items()}
    if not all(math.isfinite(value) for value in numeric.values()):
        raise ValueError("Milo diagnostic thresholds must be finite.")
    if not 0.0 < numeric["neighborhood_proportion"] <= 1.0:
        raise ValueError("Milo neighborhood proportion must lie in (0, 1].")
    for name in ("mixed_annotation_threshold", "spatial_fdr_threshold"):
        if not 0.0 <= numeric[name] <= 1.0:
            raise ValueError(f"Milo {name} must lie in [0, 1].")
    if numeric["min_abs_log2_fold_change"] < 0.0:
        raise ValueError("Milo minimum absolute log2 fold change cannot be negative.")

    table_owned = table.copy(deep=True)
    table_identity = _table_identity(table_owned)
    if table_owned["neighborhood_id"].tolist() != neighborhoods:
        raise ValueError("Milo result table and neighborhood axes differ.")
    if len(table_owned) != len(neighborhoods):
        raise ValueError("Milo result table must contain one row per neighborhood.")
    if table_owned["index_cell"].duplicated().any() or not set(table_owned["index_cell"]).issubset(observations):
        raise ValueError("Milo neighborhood index cells must be unique members of the observation axis.")

    membership_owned = _csr(
        membership,
        label="neighborhood membership",
        rows=len(observations),
        columns=len(neighborhoods),
        binary=True,
    )
    sizes = membership_owned.sum(axis=0).A1.astype(int)
    if not bool((sizes == table_owned["neighborhood_size"].to_numpy(dtype=int)).all()):
        raise ValueError("Milo membership sizes do not match the canonical result table.")
    by_neighborhood = membership_owned.tocsc()
    signatures = [
        tuple(by_neighborhood.indices[by_neighborhood.indptr[index] : by_neighborhood.indptr[index + 1]].tolist())
        for index in range(len(neighborhoods))
    ]
    if len(signatures) != len(set(signatures)):
        raise ValueError("Milo neighborhood memberships must be unique.")

    expected_graph, _ = _membership_overlap_graph(membership_owned)
    graph_owned = _csr(
        graph,
        label="neighborhood graph",
        rows=len(neighborhoods),
        columns=len(neighborhoods),
    )
    if (graph_owned != expected_graph).nnz:
        raise ValueError("Milo neighborhood graph must equal membership overlap with a zero diagonal.")
    coordinates_owned = _coordinates(representative_coordinates, neighborhoods=len(neighborhoods))

    provenance = {
        "representation_key": representation_key,
        "representation_sha256": representation_sha256,
        "coordinate_dimensions": ["component_1", "component_2"],
        "condition_key": condition_key,
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "n_neighbors": int(n_neighbors),
        "neighborhood_proportion": numeric["neighborhood_proportion"],
        "mixed_annotation_threshold": numeric["mixed_annotation_threshold"],
        "spatial_fdr_threshold": numeric["spatial_fdr_threshold"],
        "min_abs_log2_fold_change": numeric["min_abs_log2_fold_change"],
        "random_seed": int(random_seed),
    }
    metadata: dict[str, Any] = {
        "schema_version": MILO_RESULT_ARTIFACT_SCHEMA_VERSION,
        "artifact_type": MILO_RESULT_ARTIFACT_TYPE,
        "producer_node_id": MILO_RESULT_PRODUCER_NODE_ID,
        "producer_schema": MILO_RESULT_PRODUCER_SCHEMA,
        "observation_count": len(observations),
        "observation_axis_sha256": observation_hash,
        "neighborhood_count": len(neighborhoods),
        "neighborhood_axis_sha256": neighborhood_hash,
        "table_identity": table_identity,
        "membership_sha256": _matrix_sha256(membership_owned, label="milo-membership"),
        "graph_sha256": _matrix_sha256(graph_owned, label="milo-neighborhood-graph"),
        "coordinates_sha256": _coordinate_sha256(coordinates_owned),
        "provenance_sha256": _json_sha256(provenance),
    }
    metadata["artifact_fingerprint_sha256"] = _json_sha256(metadata)
    result = MiloResult._from_owned(
        table=table_owned,
        membership=membership_owned,
        graph=graph_owned,
        coordinates=coordinates_owned,
        observation_names=observations,
        neighborhood_names=neighborhoods,
        provenance=provenance,
        metadata=metadata,
    )
    validate_milo_result(result, copy_result=False)
    return result


def validate_milo_result(
    result: Any,
    *,
    exact_type: bool = True,
    copy_result: bool = True,
) -> tuple[Any, Any, Any, Any, list[str], list[str], dict[str, Any], dict[str, Any]]:
    if exact_type:
        if type(result) is not MiloResult:
            raise TypeError("Expected an exact OPENBIO_MILO_RESULT artifact.")
        payload = result.portable() if copy_result else {
            "artifact_type": MILO_RESULT_ARTIFACT_TYPE,
            "table": object.__getattribute__(result, "_table"),
            "membership": object.__getattribute__(result, "_membership"),
            "graph": object.__getattribute__(result, "_graph"),
            "coordinates": object.__getattribute__(result, "_coordinates"),
            "observation_names": object.__getattribute__(result, "_observation_names"),
            "neighborhood_names": object.__getattribute__(result, "_neighborhood_names"),
            "provenance": object.__getattribute__(result, "_provenance"),
            "metadata": object.__getattribute__(result, "_metadata"),
        }
    elif isinstance(result, Mapping):
        payload = copy.deepcopy(dict(result)) if copy_result else dict(result)
    else:
        raise TypeError("Portable Milo result must be a mapping.")
    expected_payload = {
        "artifact_type",
        "table",
        "membership",
        "graph",
        "coordinates",
        "observation_names",
        "neighborhood_names",
        "provenance",
        "metadata",
    }
    if set(payload) != expected_payload or payload["artifact_type"] != MILO_RESULT_ARTIFACT_TYPE:
        raise ValueError("Milo portable result schema is invalid.")
    metadata = payload["metadata"]
    provenance = payload["provenance"]
    expected_metadata = {
        "schema_version",
        "artifact_type",
        "producer_node_id",
        "producer_schema",
        "observation_count",
        "observation_axis_sha256",
        "neighborhood_count",
        "neighborhood_axis_sha256",
        "table_identity",
        "membership_sha256",
        "graph_sha256",
        "coordinates_sha256",
        "provenance_sha256",
        "artifact_fingerprint_sha256",
    }
    if not isinstance(metadata, Mapping) or set(metadata) != expected_metadata:
        raise ValueError("Milo result metadata schema is invalid.")
    if (
        metadata["schema_version"] != MILO_RESULT_ARTIFACT_SCHEMA_VERSION
        or metadata["artifact_type"] != MILO_RESULT_ARTIFACT_TYPE
        or metadata["producer_node_id"] != MILO_RESULT_PRODUCER_NODE_ID
        or metadata["producer_schema"] != MILO_RESULT_PRODUCER_SCHEMA
    ):
        raise ValueError("Milo result producer or schema identity is invalid.")
    metadata_without_fingerprint = dict(metadata)
    fingerprint = metadata_without_fingerprint.pop("artifact_fingerprint_sha256")
    if fingerprint != _json_sha256(metadata_without_fingerprint):
        raise ValueError("Milo result metadata fingerprint is invalid.")
    if exact_type and result.fingerprint != fingerprint:
        raise ValueError("Milo result fingerprint property is inconsistent.")
    if not isinstance(provenance, Mapping) or set(provenance) != {
        "representation_key",
        "representation_sha256",
        "coordinate_dimensions",
        "condition_key",
        "annotation_key",
        "annotation_status",
        "n_neighbors",
        "neighborhood_proportion",
        "mixed_annotation_threshold",
        "spatial_fdr_threshold",
        "min_abs_log2_fold_change",
        "random_seed",
    }:
        raise ValueError("Milo result provenance schema is invalid.")
    if provenance["coordinate_dimensions"] != ["component_1", "component_2"]:
        raise ValueError("Milo coordinate dimensions are invalid.")
    if metadata["provenance_sha256"] != _json_sha256(dict(provenance)):
        raise ValueError("Milo result provenance fingerprint is invalid.")

    observations, observation_hash = _axis(payload["observation_names"], label="observation axis")
    neighborhoods, neighborhood_hash = _axis(payload["neighborhood_names"], label="neighborhood axis")
    if (
        metadata["observation_count"] != len(observations)
        or metadata["observation_axis_sha256"] != observation_hash
        or metadata["neighborhood_count"] != len(neighborhoods)
        or metadata["neighborhood_axis_sha256"] != neighborhood_hash
    ):
        raise ValueError("Milo result axes failed their current-content checks.")
    table = payload["table"]
    if _table_identity(table) != metadata["table_identity"]:
        raise ValueError("Milo result table failed its current-content fingerprint check.")
    if table["neighborhood_id"].tolist() != neighborhoods:
        raise ValueError("Milo result table and neighborhood axes differ.")
    if table["index_cell"].duplicated().any() or not set(table["index_cell"]).issubset(observations):
        raise ValueError("Milo neighborhood index cells are inconsistent with the observation axis.")
    membership = _csr(
        payload["membership"],
        label="neighborhood membership",
        rows=len(observations),
        columns=len(neighborhoods),
        binary=True,
    )
    if _matrix_sha256(membership, label="milo-membership") != metadata["membership_sha256"]:
        raise ValueError("Milo membership failed its current-content fingerprint check.")
    sizes = membership.sum(axis=0).A1.astype(int)
    if not bool((sizes == table["neighborhood_size"].to_numpy(dtype=int)).all()):
        raise ValueError("Milo membership sizes do not match the canonical result table.")
    expected_graph, _ = _membership_overlap_graph(membership)
    graph = _csr(
        payload["graph"],
        label="neighborhood graph",
        rows=len(neighborhoods),
        columns=len(neighborhoods),
    )
    if _matrix_sha256(graph, label="milo-neighborhood-graph") != metadata["graph_sha256"]:
        raise ValueError("Milo neighborhood graph failed its current-content fingerprint check.")
    if (graph != expected_graph).nnz:
        raise ValueError("Milo neighborhood graph does not match membership overlap.")
    coordinates = _coordinates(payload["coordinates"], neighborhoods=len(neighborhoods))
    if _coordinate_sha256(coordinates) != metadata["coordinates_sha256"]:
        raise ValueError("Milo representative coordinates failed their current-content fingerprint check.")
    return (
        table.copy(deep=True) if copy_result else table,
        membership.copy() if copy_result else membership,
        graph.copy() if copy_result else graph,
        coordinates.copy() if copy_result else coordinates,
        observations,
        neighborhoods,
        copy.deepcopy(dict(provenance)),
        copy.deepcopy(dict(metadata)),
    )


__all__ = [
    "MILO_RESULT_ARTIFACT_SCHEMA_VERSION",
    "MILO_RESULT_ARTIFACT_TYPE",
    "MILO_RESULT_COLUMNS",
    "MiloResult",
    "build_milo_result",
    "validate_milo_result",
]
