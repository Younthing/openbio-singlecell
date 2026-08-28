from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pandas import DataFrame


SCENIC_RESULT_ARTIFACT_TYPE = "openbio-singlecell/scenic-result"
SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION = 1
SCENIC_RESULT_PRODUCER_NODE_ID = "OpenBioSingleCellImportPySCENICResults"
SCENIC_RESULT_PRODUCER_SCHEMA = "openbio-singlecell/pyscenic-import/v1"

SCENIC_BINARY_ARTIFACT_TYPE = "openbio-singlecell/scenic-binary"
SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION = 1
SCENIC_BINARY_PRODUCER_NODE_ID = "OpenBioSingleCellSCENICActivityBinarization"
SCENIC_BINARY_PRODUCER_SCHEMA = "openbio-singlecell/scenic-binarization/v1"

SCENIC_MEMBERSHIP_COLUMNS = (
    "regulon",
    "transcription_factor",
    "regulation",
    "context",
    "target",
    "target_weight",
    "motif_evidence_count",
    "motif_ids",
    "target_rank",
)


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _canonical_axis(values: Any, *, name: str) -> tuple[str, ...]:
    try:
        result = tuple(values)
    except TypeError as exc:
        raise TypeError(f"{name} must be an ordered collection of strings.") from exc
    if not result:
        raise ValueError(f"{name} cannot be empty.")
    for position, value in enumerate(result):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{name} identifier at position {position} must be canonical and nonblank.")
    if len(result) != len(set(result)):
        raise ValueError(f"{name} identifiers must be unique.")
    return result


def _activity_bytes(activity: DataFrame, *, numpy: Any, pandas: Any) -> tuple[tuple[str, ...], tuple[str, ...], bytes]:
    if not isinstance(activity, pandas.DataFrame):
        raise TypeError("SCENIC activity must be a pandas DataFrame.")
    observations = _canonical_axis(activity.index.tolist(), name="SCENIC observation axis")
    regulons = _canonical_axis(activity.columns.tolist(), name="SCENIC regulon axis")
    values = numpy.asarray(activity.to_numpy(dtype=numpy.float64, copy=True), dtype="<f8", order="C")
    if values.shape != (len(observations), len(regulons)):
        raise ValueError("SCENIC activity shape does not match its named axes.")
    if not bool(numpy.isfinite(values).all()):
        raise ValueError("SCENIC activity values must be finite.")
    if bool(((values < 0.0) | (values > 1.0)).any()):
        raise ValueError("SCENIC AUCell activity values must lie in [0, 1].")
    return observations, regulons, values.tobytes(order="C")


def _canonical_membership(
    membership: DataFrame,
    *,
    regulons: tuple[str, ...],
    numpy: Any,
    pandas: Any,
) -> tuple[dict[str, Any], ...]:
    if not isinstance(membership, pandas.DataFrame):
        raise TypeError("SCENIC membership must be a pandas DataFrame.")
    if membership.columns.tolist() != list(SCENIC_MEMBERSHIP_COLUMNS):
        raise ValueError(
            "SCENIC membership columns must be exactly "
            f"{list(SCENIC_MEMBERSHIP_COLUMNS)!r}."
        )
    if membership.empty:
        raise ValueError("SCENIC membership cannot be empty.")
    regulon_order = {name: position for position, name in enumerate(regulons)}
    rows: list[dict[str, Any]] = []
    seen_edges: set[tuple[str, str]] = set()
    identity_by_regulon: dict[str, tuple[str, str, str]] = {}
    for position, row in membership.iterrows():
        del position
        text_fields = {}
        for column in ("regulon", "transcription_factor", "regulation", "context", "target"):
            value = row[column]
            if not isinstance(value, str) or not value or value != value.strip():
                raise ValueError(f"SCENIC membership {column} values must be canonical nonblank strings.")
            text_fields[column] = value
        if text_fields["regulon"] not in regulon_order:
            raise ValueError("SCENIC membership references a regulon outside the activity family.")
        if text_fields["regulation"] not in {"activating", "repressing"}:
            raise ValueError("SCENIC membership regulation must be activating or repressing.")
        suffix = "(+)" if text_fields["regulation"] == "activating" else "(-)"
        if text_fields["regulon"] != f"{text_fields['transcription_factor']}{suffix}":
            raise ValueError("SCENIC regulon name, TF, and regulation are inconsistent.")
        identity = (
            text_fields["transcription_factor"],
            text_fields["regulation"],
            text_fields["context"],
        )
        previous_identity = identity_by_regulon.setdefault(text_fields["regulon"], identity)
        if previous_identity != identity:
            raise ValueError("SCENIC membership has inconsistent identity within one regulon.")
        edge = (text_fields["regulon"], text_fields["target"])
        if edge in seen_edges:
            raise ValueError("SCENIC membership contains duplicate regulon-target edges.")
        seen_edges.add(edge)
        weight = row["target_weight"]
        if isinstance(weight, bool) or not isinstance(weight, (int, float, numpy.integer, numpy.floating)):
            raise TypeError("SCENIC target weights must be numeric.")
        weight = float(weight)
        if not bool(numpy.isfinite(weight)) or weight < 0.0:
            raise ValueError("SCENIC target weights must be finite and nonnegative.")
        motif_ids = row["motif_ids"]
        if isinstance(motif_ids, str):
            raise TypeError("SCENIC motif_ids must be an ordered collection, not one string.")
        motif_ids = _canonical_axis(motif_ids, name="SCENIC motif evidence")
        motif_count = row["motif_evidence_count"]
        if isinstance(motif_count, bool) or not isinstance(motif_count, (int, numpy.integer)):
            raise TypeError("SCENIC motif_evidence_count must be an integer.")
        if int(motif_count) != len(motif_ids):
            raise ValueError("SCENIC motif evidence count does not match motif_ids.")
        target_rank = row["target_rank"]
        if isinstance(target_rank, bool) or not isinstance(target_rank, (int, numpy.integer)):
            raise TypeError("SCENIC target_rank must be an integer.")
        if int(target_rank) < 1:
            raise ValueError("SCENIC target_rank must be positive.")
        rows.append(
            {
                **text_fields,
                "target_weight": weight,
                "motif_evidence_count": int(motif_count),
                "motif_ids": list(motif_ids),
                "target_rank": int(target_rank),
            }
        )
    if set(identity_by_regulon) != set(regulons):
        raise ValueError("SCENIC membership does not cover the complete activity regulon family.")
    expected = sorted(
        rows,
        key=lambda row: (
            regulon_order[row["regulon"]],
            -row["target_weight"],
            row["target"],
        ),
    )
    if rows != expected:
        raise ValueError("SCENIC membership rows are not in canonical regulon/weight/target order.")
    expected_rank: dict[str, int] = {}
    previous_weight: dict[str, float] = {}
    position_by_regulon: dict[str, int] = {}
    for row in rows:
        regulon = row["regulon"]
        position = position_by_regulon.get(regulon, 0) + 1
        position_by_regulon[regulon] = position
        if regulon not in previous_weight or row["target_weight"] != previous_weight[regulon]:
            expected_rank[regulon] = position
            previous_weight[regulon] = row["target_weight"]
        if row["target_rank"] != expected_rank[regulon]:
            raise ValueError("SCENIC target ranks do not match descending target weights.")
    return tuple(rows)


def _strict_provenance(
    provenance: Any,
    *,
    artifact_type: str,
    schema_version: int,
    producer_node_id: str,
    producer_schema: str,
) -> dict[str, Any]:
    if not isinstance(provenance, dict):
        raise TypeError("SCENIC artifact provenance must be a dictionary.")
    result = copy.deepcopy(provenance)
    required = {
        "artifact_type": artifact_type,
        "artifact_schema_version": schema_version,
        "producer_node_id": producer_node_id,
        "producer_schema": producer_schema,
    }
    for name, expected in required.items():
        if result.get(name) != expected:
            raise ValueError(f"SCENIC artifact provenance {name!r} is invalid.")
    if artifact_type == SCENIC_RESULT_ARTIFACT_TYPE:
        expected_fields = {
            *required,
            "manifest_schema",
            "manifest_sha256",
            "run_id",
            "organism",
            "genome_build",
            "gene_namespace",
            "expression_state",
            "input_dimensions",
            "expression_fingerprints",
            "container",
            "declared_software_versions",
            "commands",
            "resources",
            "bundle_accounting",
        }
        if set(result) != expected_fields:
            raise ValueError("SCENIC result provenance fields are invalid.")
        if result["manifest_schema"] != "openbio-singlecell/pyscenic-external-run/v1":
            raise ValueError("SCENIC result manifest schema is invalid.")
        for name in ("run_id", "organism", "genome_build", "gene_namespace"):
            value = result[name]
            if not isinstance(value, str) or not value or value != value.strip():
                raise ValueError(f"SCENIC result provenance {name!r} must be canonical and nonblank.")
        if (
            not isinstance(result["expression_state"], str)
            or not result["expression_state"]
            or result["expression_state"] != result["expression_state"].strip()
        ):
            raise ValueError("SCENIC result expression_state must be canonical and nonblank.")

        def checked_sha(value: Any, *, description: str) -> str:
            if (
                not isinstance(value, str)
                or len(value) != 64
                or value != value.lower()
                or any(character not in "0123456789abcdef" for character in value)
            ):
                raise ValueError(f"{description} must be a lowercase SHA-256 digest.")
            return value

        manifest_hash = checked_sha(result["manifest_sha256"], description="SCENIC manifest hash")
        dimensions = result["input_dimensions"]
        if not isinstance(dimensions, dict) or set(dimensions) != {"cells", "genes"}:
            raise ValueError("SCENIC result input_dimensions fields are invalid.")
        for name, value in dimensions.items():
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"SCENIC result dimension {name!r} must be a positive integer.")
        fingerprints = result["expression_fingerprints"]
        if not isinstance(fingerprints, dict) or set(fingerprints) != {
            "observation_ids_sha256",
            "gene_ids_sha256",
            "matrix_sha256",
        }:
            raise ValueError("SCENIC result expression fingerprint fields are invalid.")
        for name, value in fingerprints.items():
            checked_sha(value, description=f"SCENIC expression fingerprint {name}")
        container = result["container"]
        if container is not None:
            if not isinstance(container, dict) or set(container) != {"image", "digest"}:
                raise ValueError("SCENIC result container fields are invalid.")
            if not isinstance(container["image"], str) or not container["image"]:
                raise ValueError("SCENIC result container image must be nonblank.")
            digest = container["digest"]
            if not isinstance(digest, str) or not digest.startswith("sha256:"):
                raise ValueError("SCENIC result container digest is invalid.")
            checked_sha(digest[7:], description="SCENIC container digest")
        versions = result["declared_software_versions"]
        required_versions = {"pyscenic", "ctxcore", "arboreto", "python"}
        if not isinstance(versions, dict) or not required_versions.issubset(versions):
            raise ValueError("SCENIC result declared software-version fields are invalid.")
        if versions["pyscenic"] != "0.12.1":
            raise ValueError("SCENIC result requires declared pySCENIC exactly 0.12.1.")
        for name, value in versions.items():
            if (
                not isinstance(name, str)
                or not name
                or name != name.strip()
                or not isinstance(value, str)
                or not value
                or value != value.strip()
            ):
                raise ValueError("SCENIC result declared software versions must use canonical nonblank strings.")
        commands = result["commands"]
        if not isinstance(commands, dict) or set(commands) != {"grn", "ctx", "aucell"}:
            raise ValueError("SCENIC result command fields are invalid.")
        for stage, command in commands.items():
            if (
                not isinstance(command, list)
                or len(command) < 3
                or command[:2] != ["pyscenic", stage]
                or any(not isinstance(value, str) or not value for value in command)
            ):
                raise ValueError(f"SCENIC result {stage} command provenance is invalid.")
        resources = result["resources"]
        if not isinstance(resources, list) or not resources:
            raise ValueError("SCENIC result resource provenance must be a non-empty list.")
        resource_fields = {
            "role",
            "file",
            "sha256",
            "release",
            "organism",
            "genome_build",
            "gene_namespace",
            "license",
            "citation",
        }
        for position, resource in enumerate(resources):
            if not isinstance(resource, dict) or set(resource) != resource_fields:
                raise ValueError(f"SCENIC result resource {position} fields are invalid.")
            if resource["role"] not in {"tf_list", "ranking_database", "motif_annotations"}:
                raise ValueError(f"SCENIC result resource {position} role is invalid.")
            checked_sha(resource["sha256"], description=f"SCENIC resource {position} hash")
            for name in resource_fields - {"sha256"}:
                if not isinstance(resource[name], str) or not resource[name] or resource[name] != resource[name].strip():
                    raise ValueError(f"SCENIC result resource {position} {name} is invalid.")
            if any(resource[name] != result[name] for name in ("organism", "genome_build", "gene_namespace")):
                raise ValueError(f"SCENIC result resource {position} family provenance is inconsistent.")
        accounting = result["bundle_accounting"]
        expected_accounting = {
            "manifest",
            "expression",
            "adjacency",
            "regulons",
            "aucell",
            *(f"resource:{position}:{resource['role']}" for position, resource in enumerate(resources)),
        }
        if not isinstance(accounting, dict) or set(accounting) != expected_accounting:
            raise ValueError("SCENIC result bundle-accounting fields are invalid.")
        for label, item in accounting.items():
            if not isinstance(item, dict) or set(item) != {"file", "sha256", "bytes"}:
                raise ValueError(f"SCENIC result bundle accounting for {label!r} is invalid.")
            path = item["file"]
            if (
                not isinstance(path, str)
                or not path
                or "\\" in path
                or ":" in path
                or any(part in {"", ".", ".."} for part in path.split("/"))
            ):
                raise ValueError(f"SCENIC result bundle path for {label!r} is invalid.")
            checked_sha(item["sha256"], description=f"SCENIC bundle {label} hash")
            if isinstance(item["bytes"], bool) or not isinstance(item["bytes"], int) or item["bytes"] < 1:
                raise ValueError(f"SCENIC result bundle size for {label!r} is invalid.")
        if accounting["manifest"]["sha256"] != manifest_hash:
            raise ValueError("SCENIC result manifest hash differs from bundle accounting.")
        for position, resource in enumerate(resources):
            label = f"resource:{position}:{resource['role']}"
            if accounting[label]["sha256"] != resource["sha256"] or accounting[label]["file"] != resource["file"]:
                raise ValueError(f"SCENIC result resource {position} differs from bundle accounting.")
    elif artifact_type == SCENIC_BINARY_ARTIFACT_TYPE:
        expected_fields = {
            *required,
            "source_scenic_artifact_fingerprint_sha256",
            "source_manifest_sha256",
            "random_seed",
            "algorithm",
            "strict_comparison",
            "thresholds_sha256",
        }
        if set(result) != expected_fields:
            raise ValueError("SCENIC binary provenance fields are invalid.")
        for name in (
            "source_scenic_artifact_fingerprint_sha256",
            "source_manifest_sha256",
            "thresholds_sha256",
        ):
            value = result[name]
            if (
                not isinstance(value, str)
                or len(value) != 64
                or value != value.lower()
                or any(character not in "0123456789abcdef" for character in value)
            ):
                raise ValueError(f"SCENIC binary {name} is invalid.")
        seed = result["random_seed"]
        if isinstance(seed, bool) or not isinstance(seed, int) or not 1 <= seed <= 2**31 - 1:
            raise ValueError("SCENIC binary random_seed is invalid.")
        if result["algorithm"] != "pyscenic-0.12.1-hdt-compatible-single-worker":
            raise ValueError("SCENIC binary algorithm provenance is invalid.")
        if result["strict_comparison"] != "auc > threshold":
            raise ValueError("SCENIC binary comparison provenance is invalid.")
    _canonical_json(result)
    return result


def _result_fingerprint(
    observations: tuple[str, ...],
    regulons: tuple[str, ...],
    activity_data: bytes,
    membership_rows: tuple[dict[str, Any], ...],
    provenance: dict[str, Any],
) -> str:
    digest = hashlib.sha256()
    digest.update(_canonical_json({"observations": observations, "regulons": regulons}).encode("utf-8"))
    digest.update(b"\0")
    digest.update(activity_data)
    digest.update(b"\0")
    digest.update(_canonical_json(membership_rows).encode("utf-8"))
    digest.update(b"\0")
    digest.update(_canonical_json(provenance).encode("utf-8"))
    return digest.hexdigest()


@dataclass(frozen=True, slots=True, init=False)
class SCENICResultArtifact:
    _observations: tuple[str, ...] = field(repr=False)
    _regulons: tuple[str, ...] = field(repr=False)
    _activity_data: bytes = field(repr=False)
    _membership_json: str = field(repr=False)
    _provenance_json: str = field(repr=False)
    artifact_fingerprint_sha256: str

    def __init__(self, activity: DataFrame, membership: DataFrame, provenance: dict[str, Any]) -> None:
        import numpy as np
        import pandas as pd

        observations, regulons, activity_data = _activity_bytes(activity, numpy=np, pandas=pd)
        membership_rows = _canonical_membership(
            membership,
            regulons=regulons,
            numpy=np,
            pandas=pd,
        )
        checked_provenance = _strict_provenance(
            provenance,
            artifact_type=SCENIC_RESULT_ARTIFACT_TYPE,
            schema_version=SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
            producer_node_id=SCENIC_RESULT_PRODUCER_NODE_ID,
            producer_schema=SCENIC_RESULT_PRODUCER_SCHEMA,
        )
        if checked_provenance["input_dimensions"]["cells"] != len(observations):
            raise ValueError("SCENIC result provenance cell count differs from the activity axis.")
        fingerprint = _result_fingerprint(
            observations,
            regulons,
            activity_data,
            membership_rows,
            checked_provenance,
        )
        object.__setattr__(self, "_observations", observations)
        object.__setattr__(self, "_regulons", regulons)
        object.__setattr__(self, "_activity_data", activity_data)
        object.__setattr__(self, "_membership_json", _canonical_json(membership_rows))
        object.__setattr__(self, "_provenance_json", _canonical_json(checked_provenance))
        object.__setattr__(self, "artifact_fingerprint_sha256", fingerprint)

    @property
    def observation_names(self) -> tuple[str, ...]:
        return self._observations

    @property
    def regulon_names(self) -> tuple[str, ...]:
        return self._regulons

    @property
    def activity(self) -> DataFrame:
        import numpy as np
        import pandas as pd

        values = np.frombuffer(self._activity_data, dtype="<f8").reshape(
            len(self._observations), len(self._regulons)
        )
        return pd.DataFrame(values.copy(), index=self._observations, columns=self._regulons)

    @property
    def membership(self) -> DataFrame:
        import pandas as pd

        return pd.DataFrame(json.loads(self._membership_json), columns=SCENIC_MEMBERSHIP_COLUMNS)

    @property
    def provenance(self) -> dict[str, Any]:
        return json.loads(self._provenance_json)

    def to_portable(self) -> dict[str, Any]:
        return {
            "artifact_type": SCENIC_RESULT_ARTIFACT_TYPE,
            "artifact_schema_version": SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
            "observations": list(self._observations),
            "regulons": list(self._regulons),
            "activity": self.activity.to_numpy(dtype=float).tolist(),
            "membership": json.loads(self._membership_json),
            "provenance": json.loads(self._provenance_json),
            "artifact_fingerprint_sha256": self.artifact_fingerprint_sha256,
        }


def validate_scenic_result_artifact(
    value: Any,
    *,
    exact_type: bool,
    numpy: Any,
    pandas: Any,
) -> tuple[DataFrame, DataFrame, dict[str, Any], dict[str, Any]]:
    if isinstance(value, SCENICResultArtifact):
        activity = value.activity
        membership = value.membership
        provenance = value.provenance
        claimed_fingerprint = value.artifact_fingerprint_sha256
    elif not exact_type and isinstance(value, dict):
        if set(value) != {
            "artifact_type",
            "artifact_schema_version",
            "observations",
            "regulons",
            "activity",
            "membership",
            "provenance",
            "artifact_fingerprint_sha256",
        }:
            raise ValueError("Portable SCENIC result artifact fields are invalid.")
        if value["artifact_type"] != SCENIC_RESULT_ARTIFACT_TYPE:
            raise ValueError("Portable SCENIC result artifact type is invalid.")
        if value["artifact_schema_version"] != SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION:
            raise ValueError("Portable SCENIC result artifact schema is invalid.")
        observations = _canonical_axis(value["observations"], name="SCENIC observation axis")
        regulons = _canonical_axis(value["regulons"], name="SCENIC regulon axis")
        activity = pandas.DataFrame(value["activity"], index=observations, columns=regulons)
        membership = pandas.DataFrame(value["membership"], columns=SCENIC_MEMBERSHIP_COLUMNS)
        provenance = copy.deepcopy(value["provenance"])
        claimed_fingerprint = value["artifact_fingerprint_sha256"]
    else:
        raise TypeError("SCENIC result must be the typed immutable artifact produced by the audited importer.")
    observations, regulons, activity_data = _activity_bytes(activity, numpy=numpy, pandas=pandas)
    membership_rows = _canonical_membership(
        membership,
        regulons=regulons,
        numpy=numpy,
        pandas=pandas,
    )
    provenance = _strict_provenance(
        provenance,
        artifact_type=SCENIC_RESULT_ARTIFACT_TYPE,
        schema_version=SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
        producer_node_id=SCENIC_RESULT_PRODUCER_NODE_ID,
        producer_schema=SCENIC_RESULT_PRODUCER_SCHEMA,
    )
    if provenance["input_dimensions"]["cells"] != len(observations):
        raise ValueError("SCENIC result provenance cell count differs from the activity axis.")
    observed_fingerprint = _result_fingerprint(
        observations,
        regulons,
        activity_data,
        membership_rows,
        provenance,
    )
    if not isinstance(claimed_fingerprint, str) or claimed_fingerprint != observed_fingerprint:
        raise ValueError("SCENIC result artifact fingerprint does not match its canonical payload.")
    metadata = {
        "observations": observations,
        "regulons": regulons,
        "artifact_fingerprint_sha256": observed_fingerprint,
    }
    return activity.copy(deep=True), membership.copy(deep=True), copy.deepcopy(provenance), metadata


def _binary_fingerprint(
    observations: tuple[str, ...],
    regulons: tuple[str, ...],
    binary_data: bytes,
    provenance: dict[str, Any],
) -> str:
    digest = hashlib.sha256()
    digest.update(_canonical_json({"observations": observations, "regulons": regulons}).encode("utf-8"))
    digest.update(b"\0")
    digest.update(binary_data)
    digest.update(b"\0")
    digest.update(_canonical_json(provenance).encode("utf-8"))
    return digest.hexdigest()


@dataclass(frozen=True, slots=True, init=False)
class SCENICBinaryArtifact:
    _observations: tuple[str, ...] = field(repr=False)
    _regulons: tuple[str, ...] = field(repr=False)
    _binary_data: bytes = field(repr=False)
    _provenance_json: str = field(repr=False)
    artifact_fingerprint_sha256: str

    def __init__(self, binary: DataFrame, provenance: dict[str, Any]) -> None:
        import numpy as np
        import pandas as pd

        if not isinstance(binary, pd.DataFrame):
            raise TypeError("SCENIC binary activity must be a pandas DataFrame.")
        observations = _canonical_axis(binary.index.tolist(), name="SCENIC binary observation axis")
        regulons = _canonical_axis(binary.columns.tolist(), name="SCENIC binary regulon axis")
        values = binary.to_numpy(copy=True)
        if values.dtype != np.dtype(bool):
            if not bool(np.isin(values, [0, 1, False, True]).all()):
                raise ValueError("SCENIC binary activity values must be boolean.")
            values = values.astype(bool)
        values = np.asarray(values, dtype=np.bool_, order="C")
        checked_provenance = _strict_provenance(
            provenance,
            artifact_type=SCENIC_BINARY_ARTIFACT_TYPE,
            schema_version=SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION,
            producer_node_id=SCENIC_BINARY_PRODUCER_NODE_ID,
            producer_schema=SCENIC_BINARY_PRODUCER_SCHEMA,
        )
        binary_data = values.tobytes(order="C")
        fingerprint = _binary_fingerprint(observations, regulons, binary_data, checked_provenance)
        object.__setattr__(self, "_observations", observations)
        object.__setattr__(self, "_regulons", regulons)
        object.__setattr__(self, "_binary_data", binary_data)
        object.__setattr__(self, "_provenance_json", _canonical_json(checked_provenance))
        object.__setattr__(self, "artifact_fingerprint_sha256", fingerprint)

    @property
    def binary(self) -> DataFrame:
        import numpy as np
        import pandas as pd

        values = np.frombuffer(self._binary_data, dtype=np.bool_).reshape(
            len(self._observations), len(self._regulons)
        )
        return pd.DataFrame(values.copy(), index=self._observations, columns=self._regulons)

    @property
    def provenance(self) -> dict[str, Any]:
        return json.loads(self._provenance_json)

    def to_portable(self) -> dict[str, Any]:
        return {
            "artifact_type": SCENIC_BINARY_ARTIFACT_TYPE,
            "artifact_schema_version": SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION,
            "observations": list(self._observations),
            "regulons": list(self._regulons),
            "binary": self.binary.to_numpy(dtype=bool).tolist(),
            "provenance": self.provenance,
            "artifact_fingerprint_sha256": self.artifact_fingerprint_sha256,
        }


__all__ = [
    "SCENICBinaryArtifact",
    "SCENICResultArtifact",
    "SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION",
    "SCENIC_BINARY_ARTIFACT_TYPE",
    "SCENIC_BINARY_PRODUCER_NODE_ID",
    "SCENIC_BINARY_PRODUCER_SCHEMA",
    "SCENIC_MEMBERSHIP_COLUMNS",
    "SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION",
    "SCENIC_RESULT_ARTIFACT_TYPE",
    "SCENIC_RESULT_PRODUCER_NODE_ID",
    "SCENIC_RESULT_PRODUCER_SCHEMA",
    "validate_scenic_result_artifact",
]
