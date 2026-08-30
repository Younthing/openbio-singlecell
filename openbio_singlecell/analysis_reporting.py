from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING, Any, Literal

from . import PLUGIN_VERSION
from .analysis_utils import make_summary_result

if TYPE_CHECKING:
    from .contracts import SummaryResult


REPORT_SCHEMA_VERSION = 1
ReferenceKind = Literal["method", "software", "practice", "software_documentation"]


@dataclass(frozen=True, slots=True)
class AnalysisReference:
    citation: str
    url: str
    kind: ReferenceKind
    doi: str | None = None

    def __post_init__(self) -> None:
        if not self.citation.strip():
            raise ValueError("Analysis reference citation cannot be empty.")
        if not self.url.strip():
            raise ValueError("Analysis reference URL cannot be empty.")


def _json_value(value: Any) -> Any:
    from collections.abc import Mapping, Sequence

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return None
        return value
    if hasattr(value, "item"):
        try:
            return _json_value(value.item())
        except (TypeError, ValueError):
            pass
    if hasattr(value, "tolist"):
        try:
            return _json_value(value.tolist())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_json_value(item) for item in value]
    return str(value)


def _plain_json(value: Any) -> Any:
    import math
    from collections.abc import Mapping, Sequence

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "item"):
        try:
            return _plain_json(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _plain_json(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_plain_json(item) for item in value]
    if hasattr(value, "tolist"):
        return _plain_json(value.tolist())
    return str(value)


def _package_version(package: str) -> str:
    from importlib import metadata

    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return "not-installed"


def collect_software_versions(
    packages: Sequence[str], *, openbio_version: str | None = None
) -> dict[str, str]:
    import platform

    if openbio_version is None:
        openbio_version = globals().get("PLUGIN_VERSION") or _package_version("openbio-singlecell")
    versions = {
        "python": platform.python_version(),
        "openbio-singlecell": openbio_version,
    }
    for package in dict.fromkeys(packages):
        if package in versions:
            continue
        versions[package] = _package_version(package)
    return versions


def summarize_numeric(values: Any) -> dict[str, int | float | None]:
    import numpy

    array = numpy.asarray(values, dtype=float).ravel()
    finite = array[numpy.isfinite(array)]
    missing = int(array.size - finite.size)
    if finite.size == 0:
        return {
            "n": int(array.size),
            "missing": missing,
            "min": None,
            "q1": None,
            "median": None,
            "mean": None,
            "q3": None,
            "max": None,
        }
    quantiles = numpy.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "n": int(array.size),
        "missing": missing,
        "min": float(quantiles[0]),
        "q1": float(quantiles[1]),
        "median": float(quantiles[2]),
        "mean": float(finite.mean()),
        "q3": float(quantiles[3]),
        "max": float(quantiles[4]),
    }


def make_analysis_report(
    *,
    node_id: str,
    title: str,
    operation: str,
    methods: str,
    results: str,
    key_results: Mapping[str, Any],
    parameters: Mapping[str, Any],
    references: Sequence[AnalysisReference],
    software_packages: Sequence[str],
    warnings: Sequence[str],
    limitations: Sequence[str],
    input_cells: int,
    input_genes: int,
    started_at: float,
    code: str,
    random_seed: int = 0,
) -> tuple[SummaryResult, str]:
    if not node_id.strip():
        raise ValueError("Analysis report node_id cannot be empty.")
    if not methods.strip() or not results.strip():
        raise ValueError("Analysis report methods and results text cannot be empty.")
    if not references:
        raise ValueError("Analysis report requires at least one method, software, or practice reference.")
    normalized_code = code.strip()
    if not normalized_code:
        raise ValueError("Analysis report code cannot be empty.")

    summary = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "node_id": node_id,
        "methods": methods.strip(),
        "results": results.strip(),
        "key_results": _json_value(key_results),
        "parameters": _json_value(parameters),
        "warnings": [str(warning) for warning in warnings],
        "limitations": [str(limitation) for limitation in limitations],
        "references": [_json_value(asdict(reference)) for reference in references],
        "software_versions": collect_software_versions(
            software_packages, openbio_version=PLUGIN_VERSION
        ),
    }
    # Enforce the public promise at construction time instead of relying on the UI payload adapter.
    json.dumps(summary, allow_nan=False)
    report = make_summary_result(
        summary=summary,
        title=title,
        operation=operation,
        parameters=dict(parameters),
        description=results.strip(),
        warnings=[str(warning) for warning in warnings],
        input_cells=input_cells,
        input_genes=input_genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    return report, f"{normalized_code}\n"


__all__ = [
    "AnalysisReference",
    "REPORT_SCHEMA_VERSION",
    "collect_software_versions",
    "make_analysis_report",
    "summarize_numeric",
]
