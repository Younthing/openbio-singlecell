from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pandas import DataFrame


TF_ACTIVITY_ARTIFACT_TYPE = "OPENBIO_TF_ACTIVITY"
TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION = 1
TF_ACTIVITY_PRODUCER_NODE_ID = "OpenBioSingleCellCollecTRIULM"
TF_ACTIVITY_PRODUCER_SCHEMA = "openbio-singlecell/collectri-ulm-activity/v1"


def _canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _canonical_axis(values: Sequence[Any], *, label: str) -> tuple[list[str], str]:
    canonical: list[str] = []
    for position, value in enumerate(values):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(
                f"TF activity {label} at position {position} must be a nonblank, whitespace-canonical string."
            )
        canonical.append(value)
    if len(canonical) != len(set(canonical)):
        raise ValueError(f"TF activity {label} values must be unique.")
    digest = _canonical_json_sha256(canonical)
    return canonical, digest


def _activity_frame_fingerprint(frame: Any, *, label: str, numpy: Any, pandas: Any) -> dict[str, Any]:
    if not isinstance(frame, pandas.DataFrame):
        raise TypeError(f"TF activity {label} must be a pandas DataFrame.")
    observations, observation_hash = _canonical_axis(frame.index.tolist(), label="observation identifier")
    regulators, regulator_hash = _canonical_axis(frame.columns.tolist(), label="regulator identifier")
    if not observations or not regulators:
        raise ValueError(f"TF activity {label} must contain at least one observation and regulator.")
    try:
        values = frame.to_numpy(dtype=float, copy=True)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"TF activity {label} must contain only numeric values.") from exc
    if values.shape != (len(observations), len(regulators)):
        raise RuntimeError(f"TF activity {label} shape does not match its named axes.")
    if not bool(numpy.isfinite(values).all()):
        raise ValueError(f"TF activity {label} contains non-finite values.")
    digest = hashlib.sha256()
    digest.update(f"openbio-singlecell/tf-activity-{label}/v1\0".encode("ascii"))
    digest.update(observation_hash.encode("ascii"))
    digest.update(regulator_hash.encode("ascii"))
    digest.update(numpy.ascontiguousarray(values, dtype="<f8").tobytes(order="C"))
    return {
        "observation_count": len(observations),
        "observation_axis_sha256": observation_hash,
        "regulator_count": len(regulators),
        "regulator_axis_sha256": regulator_hash,
        "content_sha256": digest.hexdigest(),
    }


class TFActivityArtifact:
    """Immutable-by-interface, tamper-evident CollecTRI ULM activity result."""

    __slots__ = ("_adjusted_pvalues", "_metadata", "_provenance", "_scores")
    artifact_type = TF_ACTIVITY_ARTIFACT_TYPE

    def __init__(
        self,
        *,
        scores: Any,
        adjusted_pvalues: Any,
        provenance: Mapping[str, Any],
        metadata: Mapping[str, Any],
    ) -> None:
        object.__setattr__(self, "_scores", scores.copy(deep=True))
        object.__setattr__(self, "_adjusted_pvalues", adjusted_pvalues.copy(deep=True))
        object.__setattr__(self, "_provenance", copy.deepcopy(dict(provenance)))
        object.__setattr__(self, "_metadata", copy.deepcopy(dict(metadata)))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("TFActivityArtifact is immutable; create a new validated artifact instead.")

    @property
    def fingerprint(self) -> str:
        return str(self._metadata["artifact_fingerprint_sha256"])

    @property
    def metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._metadata)

    @property
    def provenance(self) -> dict[str, Any]:
        return copy.deepcopy(self._provenance)

    @property
    def scores(self) -> DataFrame:
        return self._scores.copy(deep=True)

    @property
    def adjusted_pvalues(self) -> DataFrame:
        return self._adjusted_pvalues.copy(deep=True)

    def portable(self) -> dict[str, Any]:
        return {
            "artifact_type": TF_ACTIVITY_ARTIFACT_TYPE,
            "scores": self._scores.copy(deep=True),
            "adjusted_pvalues": self._adjusted_pvalues.copy(deep=True),
            "provenance": copy.deepcopy(self._provenance),
            "metadata": copy.deepcopy(self._metadata),
        }


def build_tf_activity_artifact(
    *,
    scores: Any,
    adjusted_pvalues: Any,
    provenance: Mapping[str, Any],
    numpy: Any,
    pandas: Any,
) -> TFActivityArtifact:
    if not isinstance(provenance, Mapping) or not provenance:
        raise ValueError("TF activity provenance must be a nonempty mapping.")
    score_identity = _activity_frame_fingerprint(scores, label="scores", numpy=numpy, pandas=pandas)
    pvalue_identity = _activity_frame_fingerprint(
        adjusted_pvalues, label="adjusted-pvalues", numpy=numpy, pandas=pandas
    )
    if score_identity["observation_axis_sha256"] != pvalue_identity["observation_axis_sha256"]:
        raise ValueError("TF activity score and adjusted-p-value observation axes differ.")
    if score_identity["regulator_axis_sha256"] != pvalue_identity["regulator_axis_sha256"]:
        raise ValueError("TF activity score and adjusted-p-value regulator axes differ.")
    adjusted_values = adjusted_pvalues.to_numpy(dtype=float, copy=True)
    if bool(((adjusted_values < 0.0) | (adjusted_values > 1.0)).any()):
        raise ValueError("TF activity adjusted p-values must lie in [0, 1].")
    provenance_copy = copy.deepcopy(dict(provenance))
    json.dumps(provenance_copy, ensure_ascii=False, allow_nan=False)
    metadata: dict[str, Any] = {
        "schema_version": TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION,
        "artifact_type": TF_ACTIVITY_ARTIFACT_TYPE,
        "producer_node_id": TF_ACTIVITY_PRODUCER_NODE_ID,
        "producer_schema": TF_ACTIVITY_PRODUCER_SCHEMA,
        "score_identity": score_identity,
        "adjusted_pvalue_identity": pvalue_identity,
        "provenance_sha256": _canonical_json_sha256(provenance_copy),
    }
    metadata["artifact_fingerprint_sha256"] = _canonical_json_sha256(metadata)
    result = TFActivityArtifact(
        scores=scores,
        adjusted_pvalues=adjusted_pvalues,
        provenance=provenance_copy,
        metadata=metadata,
    )
    validate_tf_activity_artifact(result)
    return result


def _portable_tf_activity_payload(result: Any, *, exact_type: bool) -> dict[str, Any]:
    if exact_type:
        if type(result) is not TFActivityArtifact:
            raise TypeError(
                "Rank TF Activities requires an exact OPENBIO_TF_ACTIVITY artifact from CollecTRI ULM."
            )
        return result.portable()
    if isinstance(result, Mapping):
        payload = copy.deepcopy(dict(result))
    elif getattr(result, "artifact_type", None) == TF_ACTIVITY_ARTIFACT_TYPE and callable(
        getattr(result, "portable", None)
    ):
        payload = result.portable()
    else:
        raise TypeError("Equivalent TF ranking requires a portable OPENBIO_TF_ACTIVITY artifact.")
    return payload


def validate_tf_activity_artifact(
    result: Any,
    *,
    exact_type: bool = True,
    numpy: Any | None = None,
    pandas: Any | None = None,
) -> tuple[Any, Any, dict[str, Any], dict[str, Any]]:
    """Validate producer identity, both matrices, provenance, and current-content fingerprints."""

    if numpy is None or pandas is None:
        import numpy as np
        import pandas as pd

        numpy = np
        pandas = pd
    payload = _portable_tf_activity_payload(result, exact_type=exact_type)
    expected_payload = {"artifact_type", "scores", "adjusted_pvalues", "provenance", "metadata"}
    if set(payload) != expected_payload:
        raise ValueError("TF activity portable artifact schema is invalid.")
    if payload["artifact_type"] != TF_ACTIVITY_ARTIFACT_TYPE:
        raise ValueError("TF activity artifact type identity is invalid.")
    metadata = payload["metadata"]
    provenance = payload["provenance"]
    if not isinstance(metadata, Mapping) or not isinstance(provenance, Mapping):
        raise ValueError("TF activity artifact metadata and provenance must be mappings.")
    expected_metadata = {
        "schema_version",
        "artifact_type",
        "producer_node_id",
        "producer_schema",
        "score_identity",
        "adjusted_pvalue_identity",
        "provenance_sha256",
        "artifact_fingerprint_sha256",
    }
    if set(metadata) != expected_metadata:
        raise ValueError("TF activity artifact metadata schema is invalid.")
    if metadata["schema_version"] != TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION:
        raise ValueError("TF activity artifact has an unsupported schema version.")
    if metadata["artifact_type"] != TF_ACTIVITY_ARTIFACT_TYPE:
        raise ValueError("TF activity metadata type identity is invalid.")
    if metadata["producer_node_id"] != TF_ACTIVITY_PRODUCER_NODE_ID:
        raise ValueError("TF activity artifact has the wrong producer node.")
    if metadata["producer_schema"] != TF_ACTIVITY_PRODUCER_SCHEMA:
        raise ValueError("TF activity artifact has an unsupported producer schema.")
    metadata_payload = {key: copy.deepcopy(value) for key, value in metadata.items() if key != "artifact_fingerprint_sha256"}
    if metadata["artifact_fingerprint_sha256"] != _canonical_json_sha256(metadata_payload):
        raise ValueError("TF activity artifact metadata fingerprint is invalid.")
    if exact_type and result.fingerprint != metadata["artifact_fingerprint_sha256"]:
        raise ValueError("TF activity artifact fingerprint property is inconsistent.")
    provenance_copy = copy.deepcopy(dict(provenance))
    json.dumps(provenance_copy, ensure_ascii=False, allow_nan=False)
    if metadata["provenance_sha256"] != _canonical_json_sha256(provenance_copy):
        raise ValueError("TF activity provenance fingerprint is invalid.")
    scores = payload["scores"]
    adjusted = payload["adjusted_pvalues"]
    score_identity = _activity_frame_fingerprint(scores, label="scores", numpy=numpy, pandas=pandas)
    adjusted_identity = _activity_frame_fingerprint(
        adjusted, label="adjusted-pvalues", numpy=numpy, pandas=pandas
    )
    if score_identity != metadata["score_identity"]:
        raise ValueError("TF activity score matrix failed its current-content fingerprint check.")
    if adjusted_identity != metadata["adjusted_pvalue_identity"]:
        raise ValueError("TF activity adjusted-p-value matrix failed its current-content fingerprint check.")
    if score_identity["observation_axis_sha256"] != adjusted_identity["observation_axis_sha256"]:
        raise ValueError("TF activity matrix observation axes differ.")
    if score_identity["regulator_axis_sha256"] != adjusted_identity["regulator_axis_sha256"]:
        raise ValueError("TF activity matrix regulator axes differ.")
    adjusted_values = adjusted.to_numpy(dtype=float, copy=True)
    if bool(((adjusted_values < 0.0) | (adjusted_values > 1.0)).any()):
        raise ValueError("TF activity adjusted p-values must lie in [0, 1].")
    return scores.copy(deep=True), adjusted.copy(deep=True), provenance_copy, copy.deepcopy(dict(metadata))


__all__ = [
    "TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION",
    "TF_ACTIVITY_ARTIFACT_TYPE",
    "TF_ACTIVITY_PRODUCER_NODE_ID",
    "TF_ACTIVITY_PRODUCER_SCHEMA",
    "TFActivityArtifact",
    "build_tf_activity_artifact",
    "validate_tf_activity_artifact",
]
