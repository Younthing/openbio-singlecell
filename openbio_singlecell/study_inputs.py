from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

CORE_STUDY_PARAMETERS_SCHEMA_VERSION = 1

DEFAULT_CORE_STUDY_PARAMETERS_JSON = (
    '{"schema_version":1,"sample_column":"sample","condition_column":"group",'
    '"batch_column":"","annotation_column":"cell_type","reference":"","comparison":""}'
)

_CORE_STUDY_PARAMETERS_FIELDS = frozenset(
    {
        "schema_version",
        "sample_column",
        "condition_column",
        "batch_column",
        "annotation_column",
        "reference",
        "comparison",
    }
)


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite number {value!r} is not valid")


def _json_object(value: str, *, label: str) -> dict[str, Any]:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be provided as a JSON string.")
    try:
        payload = json.loads(value, parse_constant=_reject_json_constant)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{label} must be valid JSON: {exc}.") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} JSON must contain an object at the top level.")
    return payload


def _require_exact_fields(payload: Mapping[str, Any], expected: frozenset[str], *, label: str) -> None:
    actual = set(payload)
    missing = sorted(expected - actual)
    unexpected = sorted(actual - expected)
    if missing or unexpected:
        details = []
        if missing:
            details.append(f"missing fields {missing}")
        if unexpected:
            details.append(f"unexpected fields {unexpected}")
        raise ValueError(f"{label} has {' and '.join(details)}.")


def _schema_version(value: Any, *, expected: int, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value != expected:
        raise ValueError(f"{label} schema_version must be {expected}.")
    return value


def _text(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be a string.")
    return value


@dataclass(frozen=True, slots=True)
class CoreStudyParameters:
    """Validated values emitted individually by the Core Study Parameters node."""

    sample_column: str
    condition_column: str
    batch_column: str
    annotation_column: str
    reference: str
    comparison: str

    def __post_init__(self) -> None:
        for label, value in (
            ("sample_column", self.sample_column),
            ("condition_column", self.condition_column),
            ("batch_column", self.batch_column),
            ("annotation_column", self.annotation_column),
            ("reference", self.reference),
            ("comparison", self.comparison),
        ):
            _text(value, label=f"Core Study Parameters {label}")

    @classmethod
    def from_json(cls, value: str) -> CoreStudyParameters:
        payload = _json_object(value, label="Core Study Parameters")
        _require_exact_fields(
            payload,
            _CORE_STUDY_PARAMETERS_FIELDS,
            label="Core Study Parameters",
        )
        _schema_version(
            payload["schema_version"],
            expected=CORE_STUDY_PARAMETERS_SCHEMA_VERSION,
            label="Core Study Parameters",
        )

        return cls(
            sample_column=_text(
                payload["sample_column"],
                label="Core Study Parameters sample_column",
            ),
            condition_column=_text(
                payload["condition_column"],
                label="Core Study Parameters condition_column",
            ),
            batch_column=_text(
                payload["batch_column"],
                label="Core Study Parameters batch_column",
            ),
            annotation_column=_text(
                payload["annotation_column"],
                label="Core Study Parameters annotation_column",
            ),
            reference=_text(
                payload["reference"],
                label="Core Study Parameters reference",
            ),
            comparison=_text(
                payload["comparison"],
                label="Core Study Parameters comparison",
            ),
        )


__all__ = [
    "CoreStudyParameters",
    "DEFAULT_CORE_STUDY_PARAMETERS_JSON",
]
