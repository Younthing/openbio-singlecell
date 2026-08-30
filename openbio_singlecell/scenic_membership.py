from __future__ import annotations

import copy
import inspect
import json
from collections.abc import Mapping
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .pyscenic_import import _software_versions
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
    from pandas import DataFrame


SCENIC_MEMBERSHIP_SUMMARY_SCHEMA = "openbio-singlecell/scenic-membership-summary/v1"
SCENIC_MEMBERSHIP_REFERENCES = [
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
        "citation": (
            "Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. "
            "Nature Methods. 2017;14:1083-1086."
        ),
        "doi": "10.1038/nmeth.4463",
        "url": "https://doi.org/10.1038/nmeth.4463",
        "kind": "method",
    },
    {
        "citation": "pySCENIC 0.12.1 df2regulons tagged implementation and target-union semantics.",
        "doi": None,
        "url": "https://github.com/aertslab/pySCENIC/blob/0.12.1/src/pyscenic/transform.py",
        "kind": "software_documentation",
    },
]


def scenic_regulon_membership(
    scenic_result: Any,
    *,
    transcription_factor: str = "",
    max_output_rows: int = 100_000,
    openbio_version: str = PLUGIN_VERSION,
    _portable_artifact: bool = False,
) -> tuple[DataFrame, dict[str, Any]]:
    import numpy as np
    import pandas as pd

    if not isinstance(transcription_factor, str) or transcription_factor != transcription_factor.strip():
        raise ValueError("SCENIC transcription_factor must be a canonical string; blank means all TFs.")
    if isinstance(max_output_rows, bool) or not isinstance(max_output_rows, int) or max_output_rows < 1:
        raise TypeError("SCENIC membership max_output_rows must be a positive integer.")
    _activity, membership, provenance, metadata = validate_scenic_result_artifact(
        scenic_result,
        exact_type=not _portable_artifact,
        numpy=np,
        pandas=pd,
        copy_result=not _portable_artifact,
    )
    complete_rows = len(membership)
    complete_regulons = membership["regulon"].nunique()
    complete_tfs = membership["transcription_factor"].nunique()
    complete_targets = membership["target"].nunique()
    if transcription_factor:
        table = membership.loc[
            membership["transcription_factor"] == transcription_factor
        ].copy(deep=True)
    else:
        table = membership.copy(deep=True)
    if len(table) > max_output_rows:
        raise ValueError(
            f"SCENIC membership view requires {len(table):,} rows, exceeding "
            f"max_output_rows={max_output_rows:,}."
        )
    table = table.reset_index(drop=True)
    if table.columns.tolist() != list(SCENIC_MEMBERSHIP_COLUMNS):
        raise RuntimeError("SCENIC membership view lost its canonical schema.")
    warnings = [
        "This node projects imported final motif-pruned regulons; it does not rerun co-expression or cisTarget inference.",
        "Target membership and weights are inferred evidence, not proof of direct binding or causality.",
    ]
    if transcription_factor and table.empty:
        warnings.append(
            f"No imported final regulon matched transcription factor {transcription_factor!r}; returned a schema-stable empty table."
        )
    multi_regulon_tfs = (
        membership[["transcription_factor", "regulon"]]
        .drop_duplicates()
        .groupby("transcription_factor", sort=False)["regulon"]
        .nunique()
    )
    queried_regulons = int(table["regulon"].nunique()) if not table.empty else 0
    queried_targets = int(table["target"].nunique()) if not table.empty else 0
    parameters = {
        "transcription_factor": transcription_factor,
        "max_output_rows": max_output_rows,
    }
    summary = {
        "schema_version": SCENIC_MEMBERSHIP_SUMMARY_SCHEMA,
        "node_id": "OpenBioSingleCellSCENICTFModules",
        "status": "final_motif_pruned_regulon_membership_view",
        "methods": (
            "Projected the immutable, safely imported pySCENIC 0.12.1 final motif-pruned regulon membership. "
            "Repeated motif-supported targets had already been consolidated with pySCENIC's maximum-weight union rule."
        ),
        "results": (
            f"Returned {len(table):,} regulon-target edges across {queried_regulons:,} regulons and "
            f"{queried_targets:,} targets from a complete family of {complete_rows:,} edges."
        ),
        "key_results": {
            "scientific_label": "Final motif-pruned SCENIC regulon membership",
            "query": transcription_factor or "all_transcription_factors",
            "query_rows": len(table),
            "query_regulons": queried_regulons,
            "query_targets": queried_targets,
            "complete_rows": complete_rows,
            "complete_regulons": int(complete_regulons),
            "complete_transcription_factors": int(complete_tfs),
            "complete_targets": int(complete_targets),
            "transcription_factors_with_multiple_regulons": [
                str(name) for name, count in multi_regulon_tfs.items() if int(count) > 1
            ],
            "regulation_counts": (
                table["regulation"].value_counts(sort=False).astype(int).to_dict()
                if not table.empty
                else {}
            ),
            "leading_edges": table.head(20).to_dict(orient="records"),
            "scenic_artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
            "external_scenic_provenance": provenance,
        },
        "parameters": parameters,
        "references": copy.deepcopy(SCENIC_MEMBERSHIP_REFERENCES),
        "software_versions": {
            **_software_versions(["numpy", "pandas"], openbio_version=openbio_version),
            "pyscenic-method": "0.12.1",
        },
        "warnings": warnings,
        "limitations": [
            "The view performs no new statistical test and cannot support a Condition-effect claim.",
            "Motif support is prior-based evidence and does not establish physical TF binding in the assayed cells.",
            "The complete family depends on the external expression universe, TF list, ranking databases, and motif annotations.",
        ],
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return table, summary


def scenic_membership_code(
    *,
    parameters: Mapping[str, Any],
    function_name: str = "scenic_final_regulon_membership",
) -> str:
    if not isinstance(function_name, str) or not function_name.isidentifier():
        raise ValueError("Generated SCENIC membership function_name must be a Python identifier.")
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
        _software_versions,
        scenic_regulon_membership,
    )
    helper_source = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    arguments = "\n".join(f"        {name}={value!r}," for name, value in portable_parameters.items())
    return f'''from __future__ import annotations

import copy
import hashlib
import importlib.metadata
import json
import platform
from collections.abc import Mapping, Sequence
from typing import Any

PLUGIN_VERSION = {PLUGIN_VERSION!r}
SCENIC_RESULT_ARTIFACT_TYPE = {SCENIC_RESULT_ARTIFACT_TYPE!r}
SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION = {SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION!r}
SCENIC_RESULT_PRODUCER_NODE_ID = {SCENIC_RESULT_PRODUCER_NODE_ID!r}
SCENIC_RESULT_PRODUCER_SCHEMA = {SCENIC_RESULT_PRODUCER_SCHEMA!r}
SCENIC_MEMBERSHIP_COLUMNS = {SCENIC_MEMBERSHIP_COLUMNS!r}
SCENIC_MEMBERSHIP_SUMMARY_SCHEMA = {SCENIC_MEMBERSHIP_SUMMARY_SCHEMA!r}
SCENIC_MEMBERSHIP_REFERENCES = {SCENIC_MEMBERSHIP_REFERENCES!r}
SCENICResultArtifact = ()

{helper_source}


def {function_name}(scenic_result):
    """Return a bounded exact view of final motif-pruned regulon membership."""
    return scenic_regulon_membership(
        scenic_result,
{arguments}
    )
'''


__all__ = [
    "SCENIC_MEMBERSHIP_SUMMARY_SCHEMA",
    "scenic_membership_code",
    "scenic_regulon_membership",
]
