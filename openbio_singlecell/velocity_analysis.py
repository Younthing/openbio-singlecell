from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .velocity_portable import (
    VELOCITY_ARTIFACT_TYPE,
    VELOCITY_PORTABLE_SCHEMA,
    VELOCITY_SCVELO_VERSION,
    VELOCITY_SUMMARY_SCHEMA,
    build_velocity_graph,
    compute_velocity_moments,
    estimate_rna_velocity,
    prepare_velocity_abundances,
    rank_recovered_dynamics,
    recover_velocity_dynamics,
    render_velocity_stream,
    validate_portable_velocity_state,
)

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


VELOCITY_ARTIFACT_SCHEMA = "openbio-singlecell/velocity-state-artifact/v1"
_SUMMARY_KEYS = {
    "schema_version",
    "node_id",
    "status",
    "methods",
    "results",
    "key_results",
    "parameters",
    "references",
    "software_versions",
    "warnings",
    "limitations",
}
_METADATA_KEYS = {
    "schema_version",
    "artifact_type",
    "portable_schema",
    "stage",
    "producer_node_id",
    "state_fingerprint_sha256",
    "summary_fingerprint_sha256",
    "artifact_fingerprint_sha256",
}


def _canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_summary(summary: Any, *, state: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(summary, Mapping) or set(summary) != _SUMMARY_KEYS:
        raise ValueError("Velocity artifact summary has an invalid exact schema.")
    normalized = copy.deepcopy(dict(summary))
    if normalized["schema_version"] != VELOCITY_SUMMARY_SCHEMA:
        raise ValueError("Velocity artifact summary schema version is unsupported.")
    if normalized["node_id"] != state["producer_node_id"]:
        raise ValueError("Velocity artifact summary producer differs from its portable state.")
    for key in ("status", "methods", "results"):
        if not isinstance(normalized[key], str) or not normalized[key].strip():
            raise ValueError(f"Velocity artifact summary {key!r} must be nonempty text.")
    if not isinstance(normalized["references"], list) or not normalized["references"]:
        raise ValueError("Velocity artifact summary requires method/software references.")
    versions = normalized["software_versions"]
    if not isinstance(versions, Mapping) or versions.get("scvelo") != VELOCITY_SCVELO_VERSION:
        raise ValueError("Velocity artifact summary must prove exact scVelo 0.3.4 execution.")
    if normalized["parameters"] != state["parameters"]:
        raise ValueError("Velocity artifact summary parameters differ from staged state provenance.")
    key_results = normalized["key_results"]
    if not isinstance(key_results, Mapping) or key_results.get(
        "state_fingerprint_sha256"
    ) != state["state_fingerprint_sha256"]:
        raise ValueError("Velocity artifact summary lacks the exact staged-state fingerprint.")
    json.dumps(normalized, ensure_ascii=False, allow_nan=False)
    return normalized


class VelocityState:
    """Immutable-by-interface, defensive, tamper-evident staged RNA-velocity artifact."""

    __slots__ = ("_adata", "_metadata", "_summary")
    artifact_type = VELOCITY_ARTIFACT_TYPE

    def __init__(
        self,
        *,
        adata: AnnData,
        summary: Mapping[str, Any],
        metadata: Mapping[str, Any],
    ) -> None:
        object.__setattr__(self, "_adata", adata.copy())
        object.__setattr__(self, "_summary", copy.deepcopy(dict(summary)))
        object.__setattr__(self, "_metadata", copy.deepcopy(dict(metadata)))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("VelocityState is immutable; create a new validated stage instead.")

    @property
    def fingerprint(self) -> str:
        return str(self._metadata["artifact_fingerprint_sha256"])

    @property
    def stage(self) -> str:
        return str(self._metadata["stage"])

    @property
    def metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._metadata)

    @property
    def summary(self) -> dict[str, Any]:
        return copy.deepcopy(self._summary)

    def portable_adata(self) -> AnnData:
        return self._adata.copy()


def _build_velocity_state(adata: AnnData, summary: Mapping[str, Any]) -> VelocityState:
    state = validate_portable_velocity_state(adata)
    normalized_summary = _validate_summary(summary, state=state)
    metadata: dict[str, Any] = {
        "schema_version": VELOCITY_ARTIFACT_SCHEMA,
        "artifact_type": VELOCITY_ARTIFACT_TYPE,
        "portable_schema": VELOCITY_PORTABLE_SCHEMA,
        "stage": state["stage"],
        "producer_node_id": state["producer_node_id"],
        "state_fingerprint_sha256": state["state_fingerprint_sha256"],
        "summary_fingerprint_sha256": _canonical_json_sha256(normalized_summary),
    }
    metadata["artifact_fingerprint_sha256"] = _canonical_json_sha256(metadata)
    result = VelocityState(adata=adata, summary=normalized_summary, metadata=metadata)
    validate_velocity_state(result)
    return result


def validate_velocity_state(
    result: Any,
    *,
    allowed_stages: tuple[str, ...] | None = None,
) -> tuple[AnnData, dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Return defensive payload, strict summary, metadata, and portable-state metadata."""

    if type(result) is not VelocityState:
        raise TypeError(
            "Velocity node requires an exact OPENBIO_VELOCITY_STATE artifact from the staged OpenBio chain."
        )
    metadata = result.metadata
    if set(metadata) != _METADATA_KEYS:
        raise ValueError("Velocity artifact metadata has an invalid exact schema.")
    if metadata["schema_version"] != VELOCITY_ARTIFACT_SCHEMA:
        raise ValueError("Velocity artifact schema version is unsupported.")
    if (
        metadata["artifact_type"] != VELOCITY_ARTIFACT_TYPE
        or result.artifact_type != VELOCITY_ARTIFACT_TYPE
    ):
        raise ValueError("Velocity artifact type identity is invalid.")
    if metadata["portable_schema"] != VELOCITY_PORTABLE_SCHEMA:
        raise ValueError("Velocity artifact portable schema is unsupported.")
    payload = result.portable_adata()
    state = validate_portable_velocity_state(payload, allowed_stages=allowed_stages)
    if metadata["stage"] != state["stage"] or result.stage != state["stage"]:
        raise ValueError("Velocity artifact stage identity is invalid.")
    if metadata["producer_node_id"] != state["producer_node_id"]:
        raise ValueError("Velocity artifact producer identity is invalid.")
    if metadata["state_fingerprint_sha256"] != state["state_fingerprint_sha256"]:
        raise ValueError("Velocity artifact payload failed its current-content fingerprint check.")
    summary = _validate_summary(result.summary, state=state)
    if metadata["summary_fingerprint_sha256"] != _canonical_json_sha256(summary):
        raise ValueError("Velocity artifact summary failed its current-content fingerprint check.")
    fingerprint = metadata["artifact_fingerprint_sha256"]
    if not (
        isinstance(fingerprint, str)
        and len(fingerprint) == 64
        and all(character in "0123456789abcdef" for character in fingerprint)
        and fingerprint == result.fingerprint
    ):
        raise ValueError("Velocity artifact fingerprint identity is invalid.")
    fingerprint_payload = {
        key: value for key, value in metadata.items() if key != "artifact_fingerprint_sha256"
    }
    if fingerprint != _canonical_json_sha256(fingerprint_payload):
        raise ValueError("Velocity artifact metadata failed its provenance fingerprint check.")
    return payload, summary, metadata, state


def run_velocity_prepare(
    adata: AnnData,
    *,
    scvelo_module: Any | None = None,
    **parameters: Any,
) -> VelocityState:
    output, summary = prepare_velocity_abundances(
        adata, scvelo_module=scvelo_module, **parameters
    )
    return _build_velocity_state(output, summary)


def run_velocity_moments(
    velocity_state: VelocityState,
    *,
    scvelo_module: Any | None = None,
    **parameters: Any,
) -> VelocityState:
    adata, _, _, _ = validate_velocity_state(velocity_state, allowed_stages=("prepared",))
    output, summary = compute_velocity_moments(
        adata, scvelo_module=scvelo_module, **parameters
    )
    return _build_velocity_state(output, summary)


def run_velocity_estimate(
    velocity_state: VelocityState,
    *,
    scvelo_module: Any | None = None,
    **parameters: Any,
) -> VelocityState:
    adata, _, _, _ = validate_velocity_state(velocity_state)
    output, summary = estimate_rna_velocity(
        adata, scvelo_module=scvelo_module, **parameters
    )
    return _build_velocity_state(output, summary)


def run_velocity_recover(
    velocity_state: VelocityState,
    *,
    scvelo_module: Any | None = None,
    **parameters: Any,
) -> VelocityState:
    adata, _, _, _ = validate_velocity_state(velocity_state)
    output, summary = recover_velocity_dynamics(
        adata, scvelo_module=scvelo_module, **parameters
    )
    return _build_velocity_state(output, summary)


def run_velocity_graph(
    velocity_state: VelocityState,
    *,
    scvelo_module: Any | None = None,
    **parameters: Any,
) -> VelocityState:
    adata, _, _, _ = validate_velocity_state(
        velocity_state, allowed_stages=("velocity_estimated",)
    )
    output, summary = build_velocity_graph(
        adata, scvelo_module=scvelo_module, **parameters
    )
    return _build_velocity_state(output, summary)


def run_velocity_ranking(
    velocity_state: VelocityState,
    **parameters: Any,
) -> tuple[DataFrame, dict[str, Any]]:
    adata, _, _, _ = validate_velocity_state(velocity_state)
    return rank_recovered_dynamics(adata, **parameters)


def run_velocity_stream(
    velocity_state: VelocityState,
    *,
    scvelo_module: Any | None = None,
    **parameters: Any,
) -> tuple[bytes, dict[str, Any]]:
    adata, _, _, _ = validate_velocity_state(
        velocity_state, allowed_stages=("velocity_graph",)
    )
    return render_velocity_stream(adata, scvelo_module=scvelo_module, **parameters)


_CODE_OPERATIONS = {
    "prepare": ("run_velocity_filter_and_normalize", "prepare_velocity_abundances", True),
    "moments": ("run_velocity_moments", "compute_velocity_moments", True),
    "estimate": ("run_velocity_estimation", "estimate_rna_velocity", True),
    "graph": ("run_velocity_graph", "build_velocity_graph", True),
    "recover": ("run_velocity_dynamics_recovery", "recover_velocity_dynamics", True),
    "ranking": ("run_velocity_fit_ranking", "rank_recovered_dynamics", False),
    "stream": ("run_velocity_stream_plot", "render_velocity_stream", True),
}


def velocity_code(operation: str, **parameters: Any) -> str:
    """Return self-contained source with one baked, portable equivalent function."""

    if operation not in _CODE_OPERATIONS:
        raise ValueError(f"Unsupported velocity code operation: {operation!r}.")
    wrapper_name, implementation_name, accepts_backend = _CODE_OPERATIONS[operation]
    normalized = copy.deepcopy(dict(parameters))
    json.dumps(normalized, ensure_ascii=False, allow_nan=False)
    source = Path(__file__).with_name("velocity_portable.py").read_text(encoding="utf-8")
    backend_argument = ", scvelo_module=None" if accepts_backend else ""
    backend_forward = ", scvelo_module=scvelo_module" if accepts_backend else ""
    wrapper = (
        f"\n\ndef {wrapper_name}(adata{backend_argument}):\n"
        f"    \"\"\"Return the portable primary result and strict summary without mutating adata.\"\"\"\n"
        f"    return {implementation_name}(adata, **{normalized!r}{backend_forward})\n"
    )
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, f"<{wrapper_name}>", "exec")
    return code


__all__ = [
    "VELOCITY_ARTIFACT_SCHEMA",
    "VelocityState",
    "run_velocity_estimate",
    "run_velocity_graph",
    "run_velocity_moments",
    "run_velocity_prepare",
    "run_velocity_ranking",
    "run_velocity_recover",
    "run_velocity_stream",
    "validate_velocity_state",
    "velocity_code",
]
