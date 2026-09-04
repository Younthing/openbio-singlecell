from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any

SCORE_ARTIFACTS_KEY = "openbio_singlecell_score_artifacts"
SCORE_ARTIFACT_SCHEMA_VERSION = 1


def _canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _canonical_string_axis(values: Sequence[Any], *, description: str) -> tuple[list[str], str]:
    normalized: list[str] = []
    for position, value in enumerate(values):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(
                f"{description} value at position {position} must be a nonblank, whitespace-canonical string."
            )
        normalized.append(value)
    if len(normalized) != len(set(normalized)):
        raise ValueError(f"{description} values must be unique.")
    return normalized, _canonical_json_sha256(normalized)


def score_frame_fingerprint(frame: Any, *, np: Any, pd: Any) -> dict[str, Any]:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Canonical score artifact data must be a pandas DataFrame.")
    observations, observation_fingerprint = _canonical_string_axis(
        frame.index.tolist(), description="Score observation identifier"
    )
    scores, score_axis_fingerprint = _canonical_string_axis(frame.columns.tolist(), description="Score name")
    if not observations:
        raise ValueError("Canonical score artifact data must contain at least one observation.")
    if not scores:
        raise ValueError("Canonical score artifact data must contain at least one named score.")
    if frame.shape != (len(observations), len(scores)):
        raise RuntimeError("Canonical score artifact matrix shape does not match its named axes.")
    digest = hashlib.sha256()
    digest.update(b"openbio-singlecell/score-frame/v1\0")
    digest.update(observation_fingerprint.encode("ascii"))
    digest.update(score_axis_fingerprint.encode("ascii"))
    digest.update(str(frame.shape).encode("ascii"))
    for start in range(0, len(observations), 1024):
        try:
            values = frame.iloc[start : start + 1024].to_numpy(dtype=float, copy=False)
        except (TypeError, ValueError) as exc:
            raise TypeError("Canonical score artifact data must contain only numeric values.") from exc
        if not bool(np.isfinite(values).all()):
            raise ValueError("Canonical score artifact data contains non-finite values.")
        digest.update(np.ascontiguousarray(values, dtype="<f8").tobytes(order="C"))
    return {
        "observation_count": len(observations),
        "observation_axis_sha256": observation_fingerprint,
        "score_count": len(scores),
        "score_names": scores,
        "score_axis_sha256": score_axis_fingerprint,
        "score_content_sha256": digest.hexdigest(),
    }


def build_score_artifact(
    *,
    frame: Any,
    feature_names: Sequence[Any],
    producer_node: str,
    method: str,
    storage: str,
    score_key: str,
    resource: Mapping[str, Any],
    expression: Mapping[str, Any],
    parameters: Mapping[str, Any],
    references: Sequence[Mapping[str, Any]],
    np: Any,
    pd: Any,
) -> dict[str, Any]:
    for value, description in (
        (producer_node, "Score producer node"),
        (method, "Score method"),
        (score_key, "Score key"),
    ):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{description} must be a nonblank, whitespace-canonical string.")
    if storage not in {"obsm", "obs"}:
        raise ValueError("Score artifact storage must be 'obsm' or 'obs'.")
    if not isinstance(resource, Mapping) or not resource:
        raise ValueError("Score artifact resource provenance must be a nonempty mapping.")
    if not isinstance(expression, Mapping) or not expression:
        raise ValueError("Score artifact expression provenance must be a nonempty mapping.")
    if not isinstance(parameters, Mapping) or not parameters:
        raise ValueError("Score artifact parameters must be a nonempty mapping.")
    if not isinstance(references, Sequence) or isinstance(references, (str, bytes)) or not references:
        raise ValueError("Score artifact references must be a nonempty sequence.")
    feature_axis, feature_axis_fingerprint = _canonical_string_axis(
        feature_names, description="Score feature identifier"
    )
    if not feature_axis:
        raise ValueError("Score artifact feature axis cannot be empty.")
    frame_identity = score_frame_fingerprint(frame, np=np, pd=pd)
    artifact: dict[str, Any] = {
        "schema_version": SCORE_ARTIFACT_SCHEMA_VERSION,
        "producer_node": producer_node,
        "method": method,
        "storage": storage,
        "score_key": score_key,
        **frame_identity,
        "feature_count": len(feature_axis),
        "feature_axis_sha256": feature_axis_fingerprint,
        "resource": copy.deepcopy(dict(resource)),
        "expression": copy.deepcopy(dict(expression)),
        "parameters": copy.deepcopy(dict(parameters)),
        "references": copy.deepcopy([dict(reference) for reference in references]),
    }
    # This is deliberately a self-authenticating integrity token, not a trust or signature claim.
    artifact["artifact_fingerprint_sha256"] = _canonical_json_sha256(artifact)
    json.dumps(artifact, allow_nan=False)
    return artifact


def store_score_artifact(output: Any, artifact: Mapping[str, Any]) -> None:
    if not hasattr(output, "uns"):
        raise TypeError("Score artifact output must expose AnnData.uns.")
    score_key = artifact.get("score_key") if isinstance(artifact, Mapping) else None
    if not isinstance(score_key, str) or not score_key:
        raise ValueError("Score artifact does not contain a valid score_key.")
    existing = output.uns.get(SCORE_ARTIFACTS_KEY, {})
    if not isinstance(existing, Mapping):
        raise ValueError(f"AnnData uns[{SCORE_ARTIFACTS_KEY!r}] must be a mapping when present.")
    updated = copy.deepcopy(dict(existing))
    updated[score_key] = copy.deepcopy(dict(artifact))
    output.uns[SCORE_ARTIFACTS_KEY] = updated


def validate_score_artifact(
    adata: Any,
    *,
    score_key: str,
    required_storage: str,
    allowed_methods: Sequence[str],
    np: Any,
    pd: Any,
) -> tuple[Any, dict[str, Any]]:
    if not isinstance(score_key, str) or not score_key or score_key != score_key.strip():
        raise ValueError("Score key must be a nonblank, whitespace-canonical string.")
    if required_storage not in {"obsm", "obs"}:
        raise ValueError("Required score storage must be 'obsm' or 'obs'.")
    artifacts = getattr(adata, "uns", {}).get(SCORE_ARTIFACTS_KEY)
    if not isinstance(artifacts, Mapping) or score_key not in artifacts:
        raise ValueError(
            f"Score {score_key!r} lacks canonical OpenBio score provenance at "
            f"uns[{SCORE_ARTIFACTS_KEY!r}]. Recompute it with a reviewed score node."
        )
    artifact = artifacts[score_key]
    if not isinstance(artifact, Mapping):
        raise ValueError(f"Score artifact for {score_key!r} must be a mapping.")
    required_fields = {
        "schema_version",
        "producer_node",
        "method",
        "storage",
        "score_key",
        "observation_count",
        "observation_axis_sha256",
        "score_count",
        "score_names",
        "score_axis_sha256",
        "score_content_sha256",
        "feature_count",
        "feature_axis_sha256",
        "resource",
        "expression",
        "parameters",
        "references",
        "artifact_fingerprint_sha256",
    }
    if set(artifact) != required_fields:
        raise ValueError(
            f"Score artifact for {score_key!r} has an invalid schema; "
            f"missing={sorted(required_fields - set(artifact))}, unknown={sorted(set(artifact) - required_fields)}."
        )
    artifact = copy.deepcopy(dict(artifact))
    if hasattr(artifact["score_names"], "tolist"):
        artifact["score_names"] = artifact["score_names"].tolist()
    if artifact["schema_version"] != SCORE_ARTIFACT_SCHEMA_VERSION:
        raise ValueError(f"Score artifact for {score_key!r} has an unsupported schema version.")
    if artifact["score_key"] != score_key or artifact["storage"] != required_storage:
        raise ValueError(f"Score artifact for {score_key!r} has inconsistent storage identity.")
    if artifact["method"] not in set(allowed_methods):
        raise ValueError(
            f"Score artifact method {artifact['method']!r} is not supported by this contrast; "
            f"expected one of {sorted(set(allowed_methods))}."
        )
    fingerprint_payload = {
        key: copy.deepcopy(value) for key, value in artifact.items() if key != "artifact_fingerprint_sha256"
    }
    if artifact["artifact_fingerprint_sha256"] != _canonical_json_sha256(fingerprint_payload):
        raise ValueError(f"Score artifact for {score_key!r} failed its provenance fingerprint check.")
    if required_storage == "obsm":
        if score_key not in adata.obsm:
            raise ValueError(f"Score artifact points to missing obsm[{score_key!r}].")
        frame = adata.obsm[score_key]
    else:
        if score_key not in adata.obs:
            raise ValueError(f"Score artifact points to missing obs[{score_key!r}].")
        frame = pd.DataFrame({score_key: adata.obs[score_key]}, index=adata.obs_names.copy())
    current = score_frame_fingerprint(frame, np=np, pd=pd)
    for field in (
        "observation_count",
        "observation_axis_sha256",
        "score_count",
        "score_names",
        "score_axis_sha256",
        "score_content_sha256",
    ):
        if current[field] != artifact[field]:
            raise ValueError(f"Score data for {score_key!r} differs from immutable provenance field {field!r}.")
    return frame, copy.deepcopy(dict(artifact))


__all__ = [
    "SCORE_ARTIFACTS_KEY",
    "SCORE_ARTIFACT_SCHEMA_VERSION",
    "build_score_artifact",
    "score_frame_fingerprint",
    "store_score_artifact",
    "validate_score_artifact",
]
