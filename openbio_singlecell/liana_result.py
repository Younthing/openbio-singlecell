from __future__ import annotations

import copy
import hashlib
import inspect
import json
import struct
from collections.abc import Mapping
from typing import Any

LIANA_RESULT_ARTIFACT_TYPE = "OPENBIO_LIANA_RESULT"
LIANA_RESULT_ARTIFACT_SCHEMA_VERSION = 1
LIANA_RESULT_PRODUCER_NODE_ID = "OpenBioSingleCellLianaCommunication"
LIANA_RESULT_PRODUCER_SCHEMA = "openbio-singlecell/liana-communication-result/v1"

LIANA_RESULT_COLUMNS = {
    "rank_aggregate": (
        "sample",
        "condition",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "lr_means",
        "cellphone_pvals",
        "expr_prod",
        "scaled_weight",
        "lr_logfc",
        "spec_weight",
        "lrscore",
        "specificity_rank",
        "magnitude_rank",
    ),
    "cellphonedb": (
        "sample",
        "condition",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "ligand",
        "receptor",
        "ligand_means",
        "ligand_props",
        "receptor_means",
        "receptor_props",
        "lr_means",
        "cellphone_pvals",
    ),
}
LIANA_RESULT_TEXT_COLUMNS = {
    "rank_aggregate": (
        "sample",
        "condition",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
    ),
    "cellphonedb": (
        "sample",
        "condition",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "ligand",
        "receptor",
    ),
}
LIANA_RESULT_GRAIN = (
    "sample",
    "source",
    "target",
    "ligand_complex",
    "receptor_complex",
)


def _liana_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _liana_result_text(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"LIANA result {label} values must be strings.")
    if value != value.strip() or not value:
        raise ValueError(f"LIANA result {label} values must be canonical nonblank strings.")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise ValueError(f"LIANA result {label} values cannot contain ASCII control characters.")
    return value


def _liana_result_table_identity(
    table: Any,
    *,
    method: str,
    numpy: Any,
    pandas: Any,
) -> dict[str, Any]:
    if method not in LIANA_RESULT_COLUMNS:
        raise ValueError(f"Unsupported LIANA result method: {method!r}.")
    if not isinstance(table, pandas.DataFrame):
        raise TypeError("LIANA result table must be a pandas DataFrame.")
    expected_columns = list(LIANA_RESULT_COLUMNS[method])
    if list(table.columns) != expected_columns:
        raise ValueError(
            f"LIANA result columns do not match the exact {method} schema: "
            f"expected {expected_columns!r}, observed {list(table.columns)!r}."
        )
    if len(table) < 1:
        raise ValueError("LIANA result table cannot be empty.")
    if not isinstance(table.index, pandas.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError("LIANA result table must use a canonical zero-based RangeIndex.")

    text_columns = set(LIANA_RESULT_TEXT_COLUMNS[method])
    numeric_columns = [column for column in expected_columns if column not in text_columns]
    for column in text_columns:
        for value in table[column].tolist():
            _liana_result_text(value, label=column)
    for column in numeric_columns:
        if not pandas.api.types.is_numeric_dtype(table[column].dtype) or pandas.api.types.is_bool_dtype(
            table[column].dtype
        ):
            raise TypeError(f"LIANA result {column} must be a non-boolean numeric column.")
        values = table[column].to_numpy(dtype=float, copy=True)
        if not bool(numpy.isfinite(values).all()):
            raise ValueError(f"LIANA result {column} contains non-finite values.")

    bounded_unit_columns = {"cellphone_pvals"}
    if method == "rank_aggregate":
        bounded_unit_columns.update({"specificity_rank", "magnitude_rank"})
    else:
        bounded_unit_columns.update({"ligand_props", "receptor_props"})
    for column in bounded_unit_columns:
        values = table[column].to_numpy(dtype=float, copy=True)
        if bool(((values < 0.0) | (values > 1.0)).any()):
            raise ValueError(f"LIANA result {column} values must lie in [0, 1].")
    nonnegative_columns = {"lr_means"}
    if method == "cellphonedb":
        nonnegative_columns.update({"ligand_means", "receptor_means"})
    for column in nonnegative_columns:
        if bool((table[column].to_numpy(dtype=float, copy=True) < 0.0).any()):
            raise ValueError(f"LIANA result {column} values cannot be negative.")
    if bool(table.duplicated(list(LIANA_RESULT_GRAIN), keep=False).any()):
        raise ValueError("LIANA result contains duplicate rows at its declared interaction grain.")
    sample_condition = table[["sample", "condition"]].drop_duplicates()
    if bool(sample_condition.duplicated("sample", keep=False).any()):
        raise ValueError("LIANA result Sample-to-Condition mapping is not one-to-one.")

    digest = hashlib.sha256()
    digest.update(b"openbio-singlecell/liana-result-table/v1\0")
    digest.update(method.encode("utf-8"))
    digest.update(b"\0")
    digest.update(json.dumps(expected_columns, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    digest.update(b"\0")
    for row in table.itertuples(index=False, name=None):
        for column, value in zip(expected_columns, row, strict=True):
            if column in text_columns:
                encoded = value.encode("utf-8")
                digest.update(b"s")
                digest.update(len(encoded).to_bytes(8, "little", signed=False))
                digest.update(encoded)
            else:
                digest.update(b"f")
                digest.update(struct.pack("<d", float(value)))
        digest.update(b"\n")
    return {
        "method": method,
        "row_count": int(len(table)),
        "column_count": len(expected_columns),
        "columns": expected_columns,
        "sample_count": int(table["sample"].nunique(dropna=False)),
        "condition_count": int(table["condition"].nunique(dropna=False)),
        "interaction_grain": list(LIANA_RESULT_GRAIN),
        "content_sha256": digest.hexdigest(),
    }


def _liana_result_provenance(
    provenance: Any,
    *,
    method: str,
    table_identity: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if not isinstance(provenance, Mapping):
        raise ValueError("LIANA result provenance must be a mapping.")
    expected = {
        "schema_version",
        "method",
        "method_profile",
        "roles",
        "organism",
        "expression",
        "design",
        "resource",
        "backend",
        "parameters",
        "result",
        "summary_sha256",
    }
    if set(provenance) != expected:
        raise ValueError("LIANA result provenance schema is invalid.")
    result = copy.deepcopy(dict(provenance))
    json.dumps(result, ensure_ascii=False, allow_nan=False)
    if result["schema_version"] != 1 or result["method"] != method:
        raise ValueError("LIANA result provenance schema version or method is invalid.")
    nested_fields = {
        "method_profile": {
            "magnitude_field",
            "magnitude_direction",
            "specificity_field",
            "specificity_direction",
        },
        "roles": {"sample_key", "condition_key", "identity_key", "annotation_status"},
        "expression": {
            "source_kind",
            "layer_name",
            "state",
            "state_evidence",
            "full_gene_completeness_verified",
            "full_gene_completeness_basis",
            "observation_axis_sha256",
            "feature_axis_sha256",
            "expression_content_sha256",
            "cells",
            "genes",
        },
        "design": {
            "sample_order",
            "condition_order",
            "sample_to_condition",
            "sample_count",
            "condition_count",
            "identity_count",
            "sample_identity_strata",
            "eligible_identity_count_by_sample",
            "excluded_strata",
        },
        "resource": {"mode", "name", "metadata", "metadata_sha256", "accounting"},
        "backend": {"package", "version", "public_method", "captured_warnings"},
        "parameters": {
            "sample_key",
            "condition_key",
            "identity_key",
            "annotation_status",
            "organism",
            "method",
            "resource_mode",
            "resource_name",
            "source_kind",
            "layer_name",
            "expression_proportion",
            "min_cells_per_identity_sample",
            "permutations",
            "random_seed",
            "jobs",
            "return_all_lrs",
            "use_raw",
            "de_method",
            "spatial_key",
            "max_output_rows",
            "max_working_memory_gib",
        },
        "result": {
            "rows",
            "columns",
            "complete_backend_family_returned",
            "tested_rows_upper_bound",
            "magnitude_field",
            "magnitude_direction",
            "specificity_field",
            "specificity_direction",
            "report_interaction_limit",
            "report_interactions",
        },
    }
    for key, fields in nested_fields.items():
        if not isinstance(result[key], Mapping) or set(result[key]) != fields:
            raise ValueError(f"LIANA result provenance {key} schema is invalid.")
    metadata_fields = {
        "name",
        "version",
        "release_date",
        "download_url",
        "organism",
        "gene_identifier_namespace",
        "scope",
        "license",
        "source_license_review",
        "citation",
    }
    resource_metadata = result["resource"]["metadata"]
    if not isinstance(resource_metadata, Mapping) or set(resource_metadata) != metadata_fields:
        raise ValueError("LIANA result resource metadata schema is invalid.")
    accounting_fields = {
        "parsed_rows",
        "unique_rows",
        "duplicate_rows_collapsed",
        "retained_rows_on_exact_feature_axis",
        "excluded_rows_missing_subunits",
        "canonical_resource_sha256",
        "effective_resource_sha256",
        "raw_size_bytes",
        "raw_file_sha256",
        "resolved_path",
    }
    accounting = result["resource"]["accounting"]
    if not isinstance(accounting, Mapping) or set(accounting) != accounting_fields:
        raise ValueError("LIANA result resource accounting schema is invalid.")
    if result["backend"]["package"] != "liana" or result["backend"]["version"] != "1.9.0":
        raise ValueError("LIANA result backend identity is invalid.")
    if result["parameters"]["method"] != method:
        raise ValueError("LIANA result parameter method is inconsistent.")
    if result["roles"]["annotation_status"] not in {"unknown", "provisional", "curated"}:
        raise ValueError("LIANA result annotation status is invalid.")
    if result["expression"]["full_gene_completeness_verified"] is not False:
        raise ValueError("LIANA result must not claim programmatic full-gene completeness proof.")
    if result["expression"]["full_gene_completeness_basis"] != "caller-selected current feature axis":
        raise ValueError("LIANA result full-gene completeness disclosure is invalid.")
    for key in (
        "observation_axis_sha256",
        "feature_axis_sha256",
        "expression_content_sha256",
    ):
        value = result["expression"][key]
        if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError(f"LIANA result expression {key} is not a canonical SHA-256.")
    for value in (
        result["resource"]["metadata_sha256"],
        accounting["canonical_resource_sha256"],
        accounting["effective_resource_sha256"],
        result["summary_sha256"],
    ):
        if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError("LIANA result provenance contains a noncanonical SHA-256.")
    if table_identity is not None:
        if result["result"]["rows"] != table_identity["row_count"]:
            raise ValueError("LIANA result provenance row count does not match its table.")
        if result["result"]["columns"] != table_identity["columns"]:
            raise ValueError("LIANA result provenance columns do not match its table.")
        if result["design"]["sample_count"] != table_identity["sample_count"]:
            raise ValueError("LIANA result provenance Sample count does not match its table.")
        if result["design"]["condition_count"] != table_identity["condition_count"]:
            raise ValueError("LIANA result provenance Condition count does not match its table.")
    return result


class LianaResult:
    """Immutable-by-interface, tamper-evident Sample-resolved LIANA result."""

    __slots__ = ("_metadata", "_provenance", "_table")
    artifact_type = LIANA_RESULT_ARTIFACT_TYPE

    def __init__(self, *, table: Any, provenance: Mapping[str, Any], metadata: Mapping[str, Any]) -> None:
        object.__setattr__(self, "_table", table.copy(deep=True))
        object.__setattr__(self, "_provenance", copy.deepcopy(dict(provenance)))
        object.__setattr__(self, "_metadata", copy.deepcopy(dict(metadata)))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("LianaResult is immutable; create a new validated artifact instead.")

    @property
    def fingerprint(self) -> str:
        return str(self._metadata["artifact_fingerprint_sha256"])

    @property
    def table(self) -> Any:
        return self._table.copy(deep=True)

    @property
    def provenance(self) -> dict[str, Any]:
        return copy.deepcopy(self._provenance)

    @property
    def metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._metadata)

    def portable(self) -> dict[str, Any]:
        return {
            "artifact_type": LIANA_RESULT_ARTIFACT_TYPE,
            "table": self._table.copy(deep=True),
            "provenance": copy.deepcopy(self._provenance),
            "metadata": copy.deepcopy(self._metadata),
        }


def build_liana_result(
    *,
    table: Any,
    method: str,
    provenance: Mapping[str, Any],
    numpy: Any,
    pandas: Any,
) -> LianaResult:
    if not isinstance(provenance, Mapping) or not provenance:
        raise ValueError("LIANA result provenance must be a nonempty mapping.")
    table_identity = _liana_result_table_identity(
        table,
        method=method,
        numpy=numpy,
        pandas=pandas,
    )
    provenance_copy = _liana_result_provenance(
        provenance,
        method=method,
        table_identity=table_identity,
    )
    metadata: dict[str, Any] = {
        "schema_version": LIANA_RESULT_ARTIFACT_SCHEMA_VERSION,
        "artifact_type": LIANA_RESULT_ARTIFACT_TYPE,
        "producer_node_id": LIANA_RESULT_PRODUCER_NODE_ID,
        "producer_schema": LIANA_RESULT_PRODUCER_SCHEMA,
        "method": method,
        "table_identity": table_identity,
        "provenance_sha256": _liana_json_sha256(provenance_copy),
    }
    metadata["artifact_fingerprint_sha256"] = _liana_json_sha256(metadata)
    result = LianaResult(table=table, provenance=provenance_copy, metadata=metadata)
    validate_liana_result(result, numpy=numpy, pandas=pandas)
    return result


def _portable_liana_result_payload(result: Any, *, exact_type: bool) -> dict[str, Any]:
    if exact_type:
        if type(result) is not LianaResult:
            raise TypeError("LIANA consumers require an exact OPENBIO_LIANA_RESULT artifact from LIANA Communication.")
        return result.portable()
    if isinstance(result, Mapping):
        return copy.deepcopy(dict(result))
    if getattr(result, "artifact_type", None) == LIANA_RESULT_ARTIFACT_TYPE and callable(
        getattr(result, "portable", None)
    ):
        return result.portable()
    raise TypeError("Equivalent LIANA consumers require a portable OPENBIO_LIANA_RESULT artifact.")


def validate_liana_result(
    result: Any,
    *,
    exact_type: bool = True,
    numpy: Any | None = None,
    pandas: Any | None = None,
) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    if numpy is None or pandas is None:
        import numpy as np
        import pandas as pd

        numpy = np
        pandas = pd
    payload = _portable_liana_result_payload(result, exact_type=exact_type)
    if set(payload) != {"artifact_type", "table", "provenance", "metadata"}:
        raise ValueError("LIANA portable result schema is invalid.")
    if payload["artifact_type"] != LIANA_RESULT_ARTIFACT_TYPE:
        raise ValueError("LIANA result artifact type identity is invalid.")
    provenance = payload["provenance"]
    metadata = payload["metadata"]
    if not isinstance(provenance, Mapping) or not isinstance(metadata, Mapping):
        raise ValueError("LIANA result provenance and metadata must be mappings.")
    expected_metadata = {
        "schema_version",
        "artifact_type",
        "producer_node_id",
        "producer_schema",
        "method",
        "table_identity",
        "provenance_sha256",
        "artifact_fingerprint_sha256",
    }
    if set(metadata) != expected_metadata:
        raise ValueError("LIANA result metadata schema is invalid.")
    if metadata["schema_version"] != LIANA_RESULT_ARTIFACT_SCHEMA_VERSION:
        raise ValueError("LIANA result has an unsupported artifact schema version.")
    if metadata["artifact_type"] != LIANA_RESULT_ARTIFACT_TYPE:
        raise ValueError("LIANA result metadata type identity is invalid.")
    if metadata["producer_node_id"] != LIANA_RESULT_PRODUCER_NODE_ID:
        raise ValueError("LIANA result has the wrong producer node.")
    if metadata["producer_schema"] != LIANA_RESULT_PRODUCER_SCHEMA:
        raise ValueError("LIANA result has an unsupported producer schema.")
    method = metadata["method"]
    if method not in LIANA_RESULT_COLUMNS or provenance.get("method") != method:
        raise ValueError("LIANA result method identity is invalid or inconsistent.")
    metadata_payload = {
        key: copy.deepcopy(value) for key, value in metadata.items() if key != "artifact_fingerprint_sha256"
    }
    if metadata["artifact_fingerprint_sha256"] != _liana_json_sha256(metadata_payload):
        raise ValueError("LIANA result metadata fingerprint is invalid.")
    if exact_type and result.fingerprint != metadata["artifact_fingerprint_sha256"]:
        raise ValueError("LIANA result fingerprint property is inconsistent.")
    provenance_copy = _liana_result_provenance(provenance, method=method)
    if metadata["provenance_sha256"] != _liana_json_sha256(provenance_copy):
        raise ValueError("LIANA result provenance fingerprint is invalid.")
    table = payload["table"]
    table_identity = _liana_result_table_identity(
        table,
        method=method,
        numpy=numpy,
        pandas=pandas,
    )
    if table_identity != metadata["table_identity"]:
        raise ValueError("LIANA result table failed its current-content fingerprint check.")
    provenance_copy = _liana_result_provenance(
        provenance_copy,
        method=method,
        table_identity=table_identity,
    )
    return table.copy(deep=True), provenance_copy, copy.deepcopy(dict(metadata))


def liana_result_portable_code() -> str:
    helpers = (
        _liana_json_sha256,
        _liana_result_text,
        _liana_result_table_identity,
        _liana_result_provenance,
        _portable_liana_result_payload,
        validate_liana_result,
    )
    source = "\n\n".join(inspect.getsource(helper) for helper in helpers)
    return f"""from __future__ import annotations

import copy
import hashlib
import json
import struct
from collections.abc import Mapping

LIANA_RESULT_ARTIFACT_TYPE = {LIANA_RESULT_ARTIFACT_TYPE!r}
LIANA_RESULT_ARTIFACT_SCHEMA_VERSION = {LIANA_RESULT_ARTIFACT_SCHEMA_VERSION!r}
LIANA_RESULT_PRODUCER_NODE_ID = {LIANA_RESULT_PRODUCER_NODE_ID!r}
LIANA_RESULT_PRODUCER_SCHEMA = {LIANA_RESULT_PRODUCER_SCHEMA!r}
LIANA_RESULT_COLUMNS = {LIANA_RESULT_COLUMNS!r}
LIANA_RESULT_TEXT_COLUMNS = {LIANA_RESULT_TEXT_COLUMNS!r}
LIANA_RESULT_GRAIN = {LIANA_RESULT_GRAIN!r}

{source}
"""


__all__ = [
    "LIANA_RESULT_ARTIFACT_SCHEMA_VERSION",
    "LIANA_RESULT_ARTIFACT_TYPE",
    "LIANA_RESULT_COLUMNS",
    "LIANA_RESULT_GRAIN",
    "LIANA_RESULT_PRODUCER_NODE_ID",
    "LIANA_RESULT_PRODUCER_SCHEMA",
    "LianaResult",
    "build_liana_result",
    "liana_result_portable_code",
    "validate_liana_result",
]
