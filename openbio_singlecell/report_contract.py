from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

CANONICAL_REPORT_FIELDS = (
    "schema_version",
    "node_id",
    "methods",
    "results",
    "key_results",
    "parameters",
    "warnings",
    "limitations",
    "references",
    "software_versions",
)


def validate_report_summary(summary: Any) -> None:
    """Validate the shared top-level shape of an analysis summary."""
    if not isinstance(summary, dict):
        raise TypeError("Analysis summary must be a mapping.")
    missing = [field for field in CANONICAL_REPORT_FIELDS if field not in summary]
    if missing:
        raise ValueError(f"Analysis summary is missing required fields: {', '.join(missing)}.")

    schema_version = summary["schema_version"]
    if not (
        isinstance(schema_version, int)
        and not isinstance(schema_version, bool)
        and schema_version > 0
        or isinstance(schema_version, str)
        and bool(schema_version.strip())
    ):
        raise ValueError("Analysis summary schema_version must be a positive integer or non-empty string.")
    for field in ("node_id", "methods", "results"):
        value = summary[field]
        if not isinstance(value, str) or not value.strip():
            raise TypeError(f"Analysis summary {field} must be a non-empty string.")
    for field in ("key_results", "parameters", "software_versions"):
        if not isinstance(summary[field], Mapping):
            raise TypeError(f"Analysis summary {field} must be a mapping.")
    for field in ("warnings", "limitations"):
        value = summary[field]
        if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
            raise TypeError(f"Analysis summary {field} must be a list of strings.")

    references = summary["references"]
    if not isinstance(references, list):
        raise TypeError("Analysis summary references must be a list.")
    if not references:
        raise ValueError("Analysis summary requires at least one reference.")
    for index, reference in enumerate(references):
        if not isinstance(reference, Mapping):
            raise TypeError(f"Analysis summary reference {index} must be a mapping.")
        missing_reference_fields = [field for field in ("citation", "url", "kind") if field not in reference]
        if missing_reference_fields:
            raise ValueError(
                f"Analysis summary reference {index} is missing required fields: "
                f"{', '.join(missing_reference_fields)}."
            )
        for field in ("citation", "url", "kind"):
            value = reference[field]
            if not isinstance(value, str) or not value.strip():
                raise TypeError(f"Analysis summary reference {index} {field} must be a non-empty string.")

    software_versions = summary["software_versions"]
    required_versions = ("python", "openbio-singlecell")
    missing_versions = [field for field in required_versions if field not in software_versions]
    if missing_versions:
        raise ValueError(f"Analysis summary is missing required software versions: {', '.join(missing_versions)}.")
    for field in required_versions:
        value = software_versions[field]
        if not isinstance(value, str) or not value.strip():
            raise TypeError(f"Analysis summary software version {field} must be a non-empty string.")

    try:
        json.dumps(summary, allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ValueError("Analysis summary must contain strict JSON-compatible values.") from error


__all__ = ["CANONICAL_REPORT_FIELDS", "validate_report_summary"]
