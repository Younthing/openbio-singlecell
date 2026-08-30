from __future__ import annotations

import copy
import inspect
import json
import math
from collections.abc import Mapping
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .analysis_reporting import _package_version, collect_software_versions
from .scenic_artifact import (
    SCENIC_MEMBERSHIP_COLUMNS,
    SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
    SCENIC_RESULT_ARTIFACT_TYPE,
    SCENIC_RESULT_PRODUCER_NODE_ID,
    SCENIC_RESULT_PRODUCER_SCHEMA,
    _activity_bytes,
    _canonical_axis,
    _canonical_json,
    _canonical_membership,
    _result_fingerprint,
    _strict_provenance,
    validate_scenic_result_artifact,
)

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


SCENIC_RSS_SUMMARY_SCHEMA = "openbio-singlecell/scenic-rss-summary/v1"
SCENIC_RSS_COLUMNS = (
    "group",
    "regulon",
    "rss",
    "rank_within_group",
    "group_cells",
    "total_cells",
)
SCENIC_RSS_REFERENCES = [
    {
        "citation": (
            "Suo S, et al. Revealing the critical regulators of cell identity in the mouse cell atlas. "
            "Cell Reports. 2018;25:1436-1445.e3."
        ),
        "doi": "10.1016/j.celrep.2018.10.045",
        "url": "https://doi.org/10.1016/j.celrep.2018.10.045",
        "kind": "method",
    },
    {
        "citation": (
            "Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network "
            "analysis. Nature Protocols. 2020;15:2247-2276."
        ),
        "doi": "10.1038/s41596-020-0336-2",
        "url": "https://doi.org/10.1038/s41596-020-0336-2",
        "kind": "method",
    },
    {
        "citation": "pySCENIC 0.12.1 regulon_specificity_scores tagged implementation.",
        "doi": None,
        "url": "https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/rss.py",
        "kind": "software_documentation",
    },
]


def _canonical_annotation(series: Any, *, pandas: Any) -> tuple[list[str], list[str], list[str]]:
    labels: list[str] = []
    observed: list[str] = []
    identities_by_display: dict[str, tuple[str, str]] = {}
    for position, value in enumerate(series.tolist()):
        missing = pandas.isna(value)
        if not isinstance(missing, bool):
            try:
                missing = bool(missing)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"SCENIC annotation at row {position} must be scalar.") from exc
        if missing:
            raise ValueError("SCENIC annotation contains missing labels.")
        if isinstance(value, str):
            if not value or value != value.strip():
                raise ValueError("SCENIC annotation string labels must be canonical and nonblank.")
            display = value
            identity = ("str", value)
        elif isinstance(value, (bool, int, float)):
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError("SCENIC annotation numeric labels must be finite.")
            display = str(value)
            identity = (type(value).__name__, repr(value))
        else:
            raise TypeError("SCENIC annotation labels must be scalar strings, booleans, or finite numbers.")
        previous = identities_by_display.get(display)
        if previous is not None and previous != identity:
            raise ValueError(
                f"SCENIC annotation labels collide after display normalization at {display!r}."
            )
        identities_by_display[display] = identity
        labels.append(display)
        if display not in observed:
            observed.append(display)
    unused: list[str] = []
    if isinstance(series.dtype, pandas.CategoricalDtype):
        for category in series.cat.categories.tolist():
            if pandas.isna(category):
                continue
            display = str(category)
            if display not in observed:
                unused.append(display)
    if not observed:
        raise ValueError("SCENIC RSS requires at least one observed annotation group.")
    return labels, observed, unused


def compute_scenic_rss(
    adata: AnnData,
    scenic_result: Any,
    *,
    annotation_key: str = "cell_type",
    annotation_status: str = "unknown",
    max_output_rows: int = 100_000,
    openbio_version: str = PLUGIN_VERSION,
    _portable_artifact: bool = False,
    _worker_owned: bool = False,
) -> tuple[DataFrame, dict[str, Any]]:
    import numpy as np
    import pandas as pd
    import scipy

    if not isinstance(annotation_key, str) or not annotation_key or annotation_key != annotation_key.strip():
        raise ValueError("SCENIC RSS annotation_key must be a canonical nonblank string.")
    if annotation_key not in adata.obs:
        raise ValueError(f"SCENIC RSS annotation column not found in obs: {annotation_key!r}.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError("SCENIC RSS annotation_status must be unknown, provisional, or curated.")
    if isinstance(max_output_rows, bool) or not isinstance(max_output_rows, int) or max_output_rows < 1:
        raise TypeError("SCENIC RSS max_output_rows must be a positive integer.")
    activity, _membership, provenance, metadata = validate_scenic_result_artifact(
        scenic_result,
        exact_type=not _portable_artifact,
        numpy=np,
        pandas=pd,
        copy_result=not _portable_artifact,
    )
    observation_names = _canonical_axis(adata.obs_names.tolist(), name="AnnData observation axis")
    if activity.index.tolist() != list(observation_names):
        raise ValueError("SCENIC result observation IDs do not align exactly to AnnData observations.")
    labels, groups, unused_categories = _canonical_annotation(adata.obs[annotation_key], pandas=pd)
    regulons = activity.columns.tolist()
    required_rows = len(groups) * len(regulons)
    if required_rows > max_output_rows:
        raise ValueError(
            f"SCENIC RSS complete grid requires {required_rows:,} rows, exceeding "
            f"max_output_rows={max_output_rows:,}."
        )
    values = activity.to_numpy(dtype=float, copy=not _worker_owned)
    sums = values.sum(axis=0)
    degenerate = [regulon for regulon, total in zip(regulons, sums, strict=True) if total <= 0.0]
    if degenerate:
        raise ValueError(f"SCENIC RSS regulons must have positive total activity: {degenerate!r}.")
    label_array = np.asarray(labels, dtype=object)
    rows: list[dict[str, Any]] = []
    for group in groups:
        indicator = (label_array == group).astype(float)
        group_cells = int(indicator.sum())
        group_rows: list[dict[str, Any]] = []
        for position, regulon in enumerate(regulons):
            score = 1.0 - scipy.spatial.distance.jensenshannon(
                values[:, position] / sums[position],
                indicator / indicator.sum(),
            )
            score = float(score)
            if not math.isfinite(score):
                raise ValueError(f"SCENIC RSS is non-finite for group {group!r}, regulon {regulon!r}.")
            group_rows.append(
                {
                    "group": group,
                    "regulon": regulon,
                    "rss": score,
                    "group_cells": group_cells,
                    "total_cells": len(observation_names),
                }
            )
        group_rows.sort(key=lambda row: (-row["rss"], row["regulon"]))
        prior_score: float | None = None
        rank = 0
        for row_position, row in enumerate(group_rows, start=1):
            if row["rss"] != prior_score:
                rank = row_position
                prior_score = row["rss"]
            row["rank_within_group"] = rank
            rows.append(row)
    table = pd.DataFrame(rows, columns=SCENIC_RSS_COLUMNS)
    table["rank_within_group"] = table["rank_within_group"].astype(np.int64)
    if len(table) != required_rows or table.columns.tolist() != list(SCENIC_RSS_COLUMNS):
        raise RuntimeError("SCENIC RSS complete grid construction failed.")
    warnings = [
        "RSS is descriptive cell-identity concentration evidence, not a p-value, effect size, or Condition contrast.",
        "Cells are not biological replicates; Sample and Technical batch are not modeled.",
    ]
    if len(groups) == 1:
        warnings.append(
            "Only one annotation group was observed; the score is numerically defined but has no between-group "
            "specificity interpretation."
        )
    if annotation_status == "unknown":
        warnings.append("Annotation status is unknown; group labels must not be presented as curated identities.")
    elif annotation_status == "provisional":
        warnings.append("Annotation labels are provisional and require independent expert review.")
    top_by_group = [
        {
            "group": group,
            "group_cells": int((label_array == group).sum()),
            "leading_regulons": table.loc[table["group"] == group].head(10).to_dict(orient="records"),
        }
        for group in groups
    ]
    parameters = {
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "max_output_rows": max_output_rows,
        "formula": "1 - scipy.spatial.distance.jensenshannon(normalized_auc, normalized_group_indicator)",
    }
    summary = {
        "schema_version": SCENIC_RSS_SUMMARY_SCHEMA,
        "node_id": "OpenBioSingleCellSCENICRegulonSpecificity",
        "status": "descriptive_regulon_specificity_evidence",
        "methods": (
            "Computed the exact pySCENIC 0.12.1 RSS definition locally: one minus SciPy Jensen-Shannon distance "
            "between each normalized nonnegative regulon-AUC vector and each normalized group indicator."
        ),
        "results": (
            f"The complete table contains {len(table):,} group-regulon scores across {len(groups):,} observed "
            f"annotation groups and {len(regulons):,} regulons."
        ),
        "key_results": {
            "scientific_label": "Descriptive annotation-associated regulon specificity",
            "cells": len(observation_names),
            "groups": groups,
            "group_count": len(groups),
            "regulons": len(regulons),
            "complete_grid_rows": len(table),
            "rss_min": float(table["rss"].min()),
            "rss_max": float(table["rss"].max()),
            "unused_categorical_levels": unused_categories,
            "top_by_group": top_by_group,
            "scenic_artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
            "external_scenic_provenance": provenance,
        },
        "parameters": parameters,
        "references": copy.deepcopy(SCENIC_RSS_REFERENCES),
        "software_versions": {
            **collect_software_versions(
                ["anndata", "numpy", "pandas", "scipy"], openbio_version=openbio_version
            ),
            "pyscenic-method": "0.12.1",
        },
        "warnings": warnings,
        "limitations": [
            "RSS does not test differential regulon activity and produces no p-value or FDR.",
            "Group size and included-cell composition influence the cell-level activity distributions.",
            "Regulon inference and AUCell scores inherit external expression/resource provenance and do not prove binding or causality.",
        ],
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return table, summary


def scenic_rss_code(
    *,
    parameters: Mapping[str, Any],
    function_name: str = "scenic_regulon_specificity",
) -> str:
    if not isinstance(function_name, str) or not function_name.isidentifier():
        raise ValueError("Generated SCENIC RSS function_name must be a Python identifier.")
    portable_parameters = dict(parameters)
    portable_parameters["_portable_artifact"] = True
    helpers = (
        _canonical_json,
        _canonical_axis,
        _activity_bytes,
        _canonical_membership,
        _strict_provenance,
        _result_fingerprint,
        validate_scenic_result_artifact,
        _package_version,
        collect_software_versions,
        _canonical_annotation,
        compute_scenic_rss,
    )
    helper_source = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    arguments = "\n".join(f"        {name}={value!r}," for name, value in portable_parameters.items())
    return f'''from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from typing import Any

PLUGIN_VERSION = {PLUGIN_VERSION!r}
SCENIC_RESULT_ARTIFACT_TYPE = {SCENIC_RESULT_ARTIFACT_TYPE!r}
SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION = {SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION!r}
SCENIC_RESULT_PRODUCER_NODE_ID = {SCENIC_RESULT_PRODUCER_NODE_ID!r}
SCENIC_RESULT_PRODUCER_SCHEMA = {SCENIC_RESULT_PRODUCER_SCHEMA!r}
SCENIC_MEMBERSHIP_COLUMNS = {SCENIC_MEMBERSHIP_COLUMNS!r}
SCENIC_RSS_SUMMARY_SCHEMA = {SCENIC_RSS_SUMMARY_SCHEMA!r}
SCENIC_RSS_COLUMNS = {SCENIC_RSS_COLUMNS!r}
SCENIC_RSS_REFERENCES = {SCENIC_RSS_REFERENCES!r}
SCENICResultArtifact = ()

{helper_source}


def {function_name}(adata, scenic_result):
    """Return the complete exact pySCENIC-0.12.1 RSS grid and strict summary."""
    return compute_scenic_rss(
        adata,
        scenic_result,
{arguments}
    )
'''


__all__ = [
    "SCENIC_RSS_COLUMNS",
    "SCENIC_RSS_SUMMARY_SCHEMA",
    "compute_scenic_rss",
    "scenic_rss_code",
]
