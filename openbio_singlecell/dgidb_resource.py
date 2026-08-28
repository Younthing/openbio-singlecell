from __future__ import annotations

import copy
import hashlib
import inspect
import textwrap
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pandas import DataFrame


DGIDB_RESOURCE_EXTENSIONS = (".csv", ".tsv")
DGIDB_RESOURCE_COLUMNS = ("drug", "gene", "sources", "evidence")
DGIDB_METADATA_FIELDS = (
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
)
DGIDB_ACCOUNTING_FIELDS = (
    "extension",
    "requested_path",
    "resolved_path",
    "size_bytes",
    "logical_data_rows",
    "raw_pairs",
    "canonical_pairs",
    "duplicate_pairs_removed",
    "drug_count",
    "gene_count",
    "source_value_count",
    "evidence_value_count",
    "rows_without_source",
    "rows_without_evidence",
    "drug_order",
    "drug_column",
    "gene_column",
    "source_column",
    "evidence_column",
    "max_file_bytes",
    "max_rows",
)
DGIDB_ARTIFACT_METADATA_FIELDS = (
    "schema_version",
    "artifact_type",
    "producer_node_id",
    "producer_schema",
    "raw_file_sha256",
    "raw_size_bytes",
    "canonical_content_fingerprint_sha256",
    "metadata_fingerprint_sha256",
    "accounting_fingerprint_sha256",
    "canonical_rows",
    "artifact_fingerprint_sha256",
)
DGIDB_ARTIFACT_SCHEMA_VERSION = 1
DGIDB_ARTIFACT_TYPE = "OPENBIO_DGIDB_RESOURCE"
DGIDB_PRODUCER_NODE_ID = "OpenBioSingleCellDGIdbAnnotation"
DGIDB_PRODUCER_SCHEMA = "openbio-dgidb-resource/v1"
DGIDB_REQUIRED_ORGANISM = "Homo sapiens"
DGIDB_REQUIRED_NAMESPACE = "HGNC symbol"


def _standalone_dgidb_json_sha256(value):
    import hashlib
    import json

    payload = json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _standalone_validate_dgidb_metadata(value):
    import json
    from collections.abc import Mapping
    from datetime import date
    from urllib.parse import urlparse

    metadata_fields = (
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
    )

    if isinstance(value, str):
        if len(value.encode("utf-8")) > 65_536:
            raise ValueError("DGIdb resource metadata exceeds the 64 KiB guard.")

        def reject_constant(constant):
            raise ValueError(f"DGIdb resource metadata contains invalid JSON constant {constant!r}.")

        def unique_object(pairs):
            result = {}
            for key, item in pairs:
                if key in result:
                    raise ValueError(f"DGIdb resource metadata contains duplicate key {key!r}.")
                result[key] = item
            return result

        try:
            parsed = json.loads(
                value,
                object_pairs_hook=unique_object,
                parse_constant=reject_constant,
            )
        except json.JSONDecodeError as exc:
            raise ValueError(f"DGIdb resource metadata is invalid JSON: {exc.msg}.") from exc
    elif isinstance(value, Mapping):
        parsed = dict(value)
    else:
        raise TypeError("DGIdb resource metadata must be a JSON string or mapping.")

    if not isinstance(parsed, dict):
        raise TypeError("DGIdb resource metadata must decode to an object.")
    expected = set(metadata_fields)
    actual = set(parsed)
    if actual != expected:
        raise ValueError(
            "DGIdb resource metadata schema mismatch; "
            f"missing={sorted(expected - actual)}, unknown={sorted(actual - expected)}."
        )
    metadata = {}
    for field in metadata_fields:
        item = parsed[field]
        if not isinstance(item, str) or not item.strip():
            raise ValueError(f"DGIdb resource metadata field {field!r} must be a nonblank string.")
        if item != item.strip():
            raise ValueError(f"DGIdb resource metadata field {field!r} has surrounding whitespace.")
        if "\x00" in item:
            raise ValueError(f"DGIdb resource metadata field {field!r} contains a NUL character.")
        metadata[field] = item

    if "dgidb" not in metadata["name"].casefold():
        raise ValueError("DGIdb resource metadata name must identify DGIdb.")
    try:
        parsed_date = date.fromisoformat(metadata["release_date"])
    except ValueError as exc:
        raise ValueError("DGIdb release_date must be an ISO 8601 calendar date (YYYY-MM-DD).") from exc
    if parsed_date.isoformat() != metadata["release_date"]:
        raise ValueError("DGIdb release_date must use canonical YYYY-MM-DD form.")
    parsed_url = urlparse(metadata["download_url"])
    if parsed_url.scheme not in {"https", "http"} or not parsed_url.netloc:
        raise ValueError("DGIdb download_url must be an absolute HTTP(S) URL.")
    if metadata["organism"] != "Homo sapiens":
        raise ValueError("DGIdb downstream analysis currently requires organism 'Homo sapiens'.")
    if metadata["gene_identifier_namespace"] != "HGNC symbol":
        raise ValueError("DGIdb downstream analysis currently requires gene namespace 'HGNC symbol'.")
    return metadata


def _standalone_validate_dgidb_table(table):
    import json

    import pandas as pd

    columns = ["drug", "gene", "sources", "evidence"]
    if not isinstance(table, pd.DataFrame):
        raise TypeError("DGIdb resource table must be a pandas DataFrame.")
    if table.empty:
        raise ValueError("DGIdb resource table cannot be empty.")
    if list(table.columns) != columns:
        raise ValueError(f"DGIdb resource table columns must be exactly {columns} in order.")
    if list(table.index) != list(range(len(table))):
        raise ValueError("DGIdb resource table must use a canonical zero-based RangeIndex.")
    canonical = table.copy(deep=True)
    for column in ("drug", "gene"):
        values = canonical[column].tolist()
        for row_index, value in enumerate(values, start=1):
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"DGIdb resource {column} row {row_index} must be a nonblank string.")
            if value != value.strip():
                raise ValueError(f"DGIdb resource {column} row {row_index} has surrounding whitespace.")
            if any(character in value for character in ("\x00", "\r", "\n")):
                raise ValueError(f"DGIdb resource {column} row {row_index} contains a control character.")
    if canonical.duplicated(["drug", "gene"]).any():
        raise ValueError("DGIdb canonical resource contains duplicate drug-gene pairs.")

    for column in ("sources", "evidence"):
        for row_index, value in enumerate(canonical[column].tolist(), start=1):
            if not isinstance(value, str):
                raise TypeError(f"DGIdb resource {column} row {row_index} must be canonical JSON text.")
            try:
                values = json.loads(value)
            except json.JSONDecodeError as exc:
                raise ValueError(f"DGIdb resource {column} row {row_index} is invalid JSON.") from exc
            if not isinstance(values, list):
                raise ValueError(f"DGIdb resource {column} row {row_index} must encode a JSON list.")
            if any(
                not isinstance(item, str)
                or not item
                or item != item.strip()
                or any(character in item for character in ("\x00", "\r", "\n"))
                for item in values
            ):
                raise ValueError(f"DGIdb resource {column} row {row_index} contains an invalid value.")
            if len(values) != len(set(values)):
                raise ValueError(f"DGIdb resource {column} row {row_index} contains duplicate values.")
            rendered = json.dumps(values, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
            if value != rendered:
                raise ValueError(f"DGIdb resource {column} row {row_index} is not canonical JSON.")
    return canonical


def _standalone_dgidb_table_fingerprint(table):
    canonical = _standalone_validate_dgidb_table(table)
    rows = canonical.loc[:, ["drug", "gene", "sources", "evidence"]].values.tolist()
    return _standalone_dgidb_json_sha256(
        {
            "schema": "openbio-dgidb-canonical-table/v1",
            "columns": ["drug", "gene", "sources", "evidence"],
            "rows": rows,
        }
    )


def _standalone_validate_dgidb_accounting(accounting, table):
    import json
    from collections.abc import Mapping

    expected_fields = (
        "extension",
        "requested_path",
        "resolved_path",
        "size_bytes",
        "logical_data_rows",
        "raw_pairs",
        "canonical_pairs",
        "duplicate_pairs_removed",
        "drug_count",
        "gene_count",
        "source_value_count",
        "evidence_value_count",
        "rows_without_source",
        "rows_without_evidence",
        "drug_order",
        "drug_column",
        "gene_column",
        "source_column",
        "evidence_column",
        "max_file_bytes",
        "max_rows",
    )
    if not isinstance(accounting, Mapping):
        raise TypeError("DGIdb resource accounting must be a mapping.")
    actual = set(accounting)
    expected = set(expected_fields)
    if actual != expected:
        raise ValueError(
            "DGIdb resource accounting schema mismatch; "
            f"missing={sorted(expected - actual)}, unknown={sorted(actual - expected)}."
        )
    result = dict(accounting)
    for field in (
        "logical_data_rows",
        "raw_pairs",
        "canonical_pairs",
        "duplicate_pairs_removed",
        "drug_count",
        "gene_count",
        "source_value_count",
        "evidence_value_count",
        "rows_without_source",
        "rows_without_evidence",
    ):
        value = result[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"DGIdb resource accounting field {field!r} must be a nonnegative integer.")
    for field in ("size_bytes", "max_file_bytes", "max_rows"):
        value = result[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"DGIdb resource accounting field {field!r} must be a positive integer.")
    if result["extension"] not in {".csv", ".tsv"}:
        raise ValueError("DGIdb resource accounting extension is unsupported.")
    for field in ("requested_path", "resolved_path", "drug_column", "gene_column"):
        value = result[field]
        if not isinstance(value, str) or not value.strip() or value != value.strip():
            raise ValueError(
                f"DGIdb resource accounting field {field!r} must be nonblank whitespace-canonical text."
            )
    for field in ("source_column", "evidence_column"):
        value = result[field]
        if value is not None and (
            not isinstance(value, str) or not value.strip() or value != value.strip()
        ):
            raise ValueError(
                f"DGIdb resource accounting field {field!r} must be null or whitespace-canonical text."
            )
    active_columns = [
        item
        for item in (
            result["drug_column"],
            result["gene_column"],
            result["source_column"],
            result["evidence_column"],
        )
        if item is not None
    ]
    if len(active_columns) != len(set(active_columns)):
        raise ValueError("DGIdb resource accounting semantic columns must be distinct.")
    if not result["resolved_path"].casefold().endswith(result["extension"]):
        raise ValueError("DGIdb resource accounting extension differs from resolved_path.")
    drug_order = result["drug_order"]
    if not isinstance(drug_order, list) or any(not isinstance(item, str) for item in drug_order):
        raise ValueError("DGIdb resource accounting drug_order must be a list of strings.")
    if len(drug_order) != len(set(drug_order)):
        raise ValueError("DGIdb resource accounting drug_order must be unique.")
    observed_order = list(dict.fromkeys(table["drug"].tolist()))
    if drug_order != observed_order:
        raise ValueError("DGIdb resource accounting drug_order differs from canonical table order.")
    if result["canonical_pairs"] != len(table):
        raise ValueError("DGIdb canonical-pair accounting differs from the current resource table.")
    if result["raw_pairs"] - result["duplicate_pairs_removed"] != len(table):
        raise ValueError("DGIdb duplicate-pair accounting is internally inconsistent.")
    if result["logical_data_rows"] != result["raw_pairs"]:
        raise ValueError("DGIdb logical-row accounting differs from raw pair count.")
    if result["size_bytes"] > result["max_file_bytes"]:
        raise ValueError("DGIdb resource size exceeds its recorded max_file_bytes guard.")
    if result["logical_data_rows"] > result["max_rows"]:
        raise ValueError("DGIdb logical-row count exceeds its recorded max_rows guard.")
    if result["drug_count"] != len(drug_order):
        raise ValueError("DGIdb drug-count accounting differs from the current resource table.")
    if result["gene_count"] != int(table["gene"].nunique()):
        raise ValueError("DGIdb gene-count accounting differs from the current resource table.")
    source_value_count = sum(len(json.loads(value)) for value in table["sources"].tolist())
    evidence_value_count = sum(len(json.loads(value)) for value in table["evidence"].tolist())
    if result["source_value_count"] != source_value_count:
        raise ValueError("DGIdb source-value accounting differs from the current resource table.")
    if result["evidence_value_count"] != evidence_value_count:
        raise ValueError("DGIdb evidence-value accounting differs from the current resource table.")
    for column, count_field, missing_field in (
        ("sources", "source_value_count", "rows_without_source"),
        ("evidence", "evidence_value_count", "rows_without_evidence"),
    ):
        missing_rows = result[missing_field]
        if missing_rows > result["logical_data_rows"]:
            raise ValueError(f"DGIdb {missing_field} exceeds the logical-row count.")
        if result[count_field] > result["logical_data_rows"] - missing_rows:
            raise ValueError(f"DGIdb {count_field} exceeds the number of nonblank logical rows.")
        empty_canonical_pairs = sum(not json.loads(value) for value in table[column].tolist())
        if empty_canonical_pairs > missing_rows:
            raise ValueError(
                f"DGIdb {missing_field} is smaller than the number of canonical pairs lacking provenance."
            )
    for column_field, count_field, missing_field in (
        ("source_column", "source_value_count", "rows_without_source"),
        ("evidence_column", "evidence_value_count", "rows_without_evidence"),
    ):
        if result[column_field] is None and (
            result[count_field] != 0 or result[missing_field] != result["logical_data_rows"]
        ):
            raise ValueError(
                f"DGIdb {column_field} is null but its value/missing-row accounting is inconsistent."
            )
    return result


def _standalone_build_dgidb_artifact_metadata(table, metadata, accounting, raw_sha256, raw_size_bytes):
    if not (
        isinstance(raw_sha256, str)
        and len(raw_sha256) == 64
        and all(character in "0123456789abcdef" for character in raw_sha256)
    ):
        raise ValueError("DGIdb raw file SHA-256 must be 64 lowercase hexadecimal characters.")
    if isinstance(raw_size_bytes, bool) or not isinstance(raw_size_bytes, int) or raw_size_bytes < 1:
        raise ValueError("DGIdb raw file size must be a positive integer.")
    canonical = _standalone_validate_dgidb_table(table)
    validated_metadata = _standalone_validate_dgidb_metadata(metadata)
    validated_accounting = _standalone_validate_dgidb_accounting(accounting, canonical)
    if validated_accounting["size_bytes"] != raw_size_bytes:
        raise ValueError("DGIdb raw-size accounting differs from artifact metadata.")
    artifact = {
        "schema_version": 1,
        "artifact_type": "OPENBIO_DGIDB_RESOURCE",
        "producer_node_id": "OpenBioSingleCellDGIdbAnnotation",
        "producer_schema": "openbio-dgidb-resource/v1",
        "raw_file_sha256": raw_sha256,
        "raw_size_bytes": raw_size_bytes,
        "canonical_content_fingerprint_sha256": _standalone_dgidb_table_fingerprint(canonical),
        "metadata_fingerprint_sha256": _standalone_dgidb_json_sha256(validated_metadata),
        "accounting_fingerprint_sha256": _standalone_dgidb_json_sha256(validated_accounting),
        "canonical_rows": len(canonical),
    }
    artifact["artifact_fingerprint_sha256"] = _standalone_dgidb_json_sha256(artifact)
    return artifact


def _standalone_validate_portable_dgidb_resource(table, metadata, accounting, artifact_metadata):
    from collections.abc import Mapping

    expected_fields = {
        "schema_version",
        "artifact_type",
        "producer_node_id",
        "producer_schema",
        "raw_file_sha256",
        "raw_size_bytes",
        "canonical_content_fingerprint_sha256",
        "metadata_fingerprint_sha256",
        "accounting_fingerprint_sha256",
        "canonical_rows",
        "artifact_fingerprint_sha256",
    }
    if not isinstance(artifact_metadata, Mapping):
        raise TypeError("DGIdb artifact metadata must be a mapping.")
    actual = set(artifact_metadata)
    if actual != expected_fields:
        raise ValueError(
            "DGIdb artifact metadata schema mismatch; "
            f"missing={sorted(expected_fields - actual)}, unknown={sorted(actual - expected_fields)}."
        )
    supplied = dict(artifact_metadata)
    rebuilt = _standalone_build_dgidb_artifact_metadata(
        table,
        metadata,
        accounting,
        supplied["raw_file_sha256"],
        supplied["raw_size_bytes"],
    )
    if supplied != rebuilt:
        mismatched = sorted(key for key in expected_fields if supplied.get(key) != rebuilt.get(key))
        raise ValueError(f"DGIdb artifact is tampered or internally inconsistent: {mismatched}.")
    return (
        _standalone_validate_dgidb_table(table),
        _standalone_validate_dgidb_metadata(metadata),
        _standalone_validate_dgidb_accounting(accounting, table),
        rebuilt,
    )


def _standalone_load_dgidb_snapshot(
    path,
    requested_path,
    resource_metadata_json,
    drug_column="drug_claim_name",
    gene_column="gene_claim_name",
    source_column="interaction_claim_source",
    evidence_column="interaction_types",
    max_file_bytes=536870912,
    max_rows=2000000,
    expected_sha256=None,
):
    import csv
    import hashlib
    import json
    import os

    import pandas as pd

    operation = "Load DGIdb Resource"
    for name, value in (("max_file_bytes", max_file_bytes), ("max_rows", max_rows)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{operation} {name} must be a positive integer.")
    if not isinstance(requested_path, str) or not requested_path.strip():
        raise ValueError(f"{operation} requested_path must be nonblank text.")
    resolved_path = os.path.realpath(os.fspath(path))
    extension = os.path.splitext(resolved_path)[1].lower()
    if extension not in {".csv", ".tsv"}:
        raise ValueError(f"{operation} accepts only .csv or .tsv snapshots.")
    if not os.path.isfile(resolved_path):
        raise FileNotFoundError(f"{operation} snapshot not found: {requested_path!r}.")

    def canonical_column(value, label, *, optional=False):
        if optional and value in {None, ""}:
            return None
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{operation} {label} must be nonblank text.")
        if value != value.strip():
            raise ValueError(f"{operation} {label} has surrounding whitespace.")
        return value

    drug_column = canonical_column(drug_column, "drug_column")
    gene_column = canonical_column(gene_column, "gene_column")
    source_column = canonical_column(source_column, "source_column", optional=True)
    evidence_column = canonical_column(evidence_column, "evidence_column", optional=True)
    active_columns = [item for item in (drug_column, gene_column, source_column, evidence_column) if item is not None]
    if len(active_columns) != len(set(active_columns)):
        raise ValueError(f"{operation} semantic column names must be distinct.")

    metadata = _standalone_validate_dgidb_metadata(resource_metadata_json)
    stat_before = os.stat(resolved_path)
    size_bytes = int(stat_before.st_size)
    if size_bytes < 1:
        raise ValueError(f"{operation} snapshot is empty.")
    if size_bytes > max_file_bytes:
        raise ValueError(
            f"{operation} snapshot is {size_bytes:,} bytes, exceeding max_file_bytes={max_file_bytes:,}."
        )

    def stream_sha256():
        digest = hashlib.sha256()
        with open(resolved_path, "rb") as binary:
            for chunk in iter(lambda: binary.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    sha_before = stream_sha256()
    if expected_sha256 is not None:
        if not (
            isinstance(expected_sha256, str)
            and len(expected_sha256) == 64
            and all(character in "0123456789abcdef" for character in expected_sha256)
        ):
            raise ValueError(f"{operation} expected_sha256 must be 64 lowercase hexadecimal characters.")
        if sha_before != expected_sha256:
            raise ValueError(f"{operation} snapshot SHA-256 differs from the pinned expected bytes.")

    separator = "\t" if extension == ".tsv" else ","
    pairs = {}
    pair_order = []
    seen_drug_spellings = {}
    seen_gene_spellings = {}
    logical_rows = 0
    rows_without_source = 0
    rows_without_evidence = 0

    def identifier(raw, label, row_number, spellings):
        if not isinstance(raw, str):
            raise TypeError(f"{operation} {label} at row {row_number} must be text.")
        value = raw.strip()
        if not value:
            raise ValueError(f"{operation} {label} at row {row_number} is blank.")
        if any(character in value for character in ("\x00", "\r", "\n")):
            raise ValueError(f"{operation} {label} at row {row_number} contains a control character.")
        previous = spellings.get(value)
        if previous is not None and previous != raw:
            raise ValueError(
                f"{operation} {label} normalization collision at row {row_number}: {previous!r} and {raw!r}."
            )
        spellings[value] = raw
        return value

    def provenance_value(raw, label, row_number):
        if not isinstance(raw, str):
            raise TypeError(f"{operation} {label} at row {row_number} must be text.")
        value = raw.strip()
        if "\x00" in value:
            raise ValueError(f"{operation} {label} at row {row_number} contains a NUL character.")
        return value or None

    with open(resolved_path, encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter=separator, strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{operation} snapshot is missing a header row.") from exc
        if not header or any(not isinstance(item, str) or not item.strip() for item in header):
            raise ValueError(f"{operation} snapshot header contains a blank column name.")
        if len(header) != len(set(header)):
            raise ValueError(f"{operation} snapshot contains duplicate column names.")
        missing = [item for item in active_columns if item not in header]
        if missing:
            raise ValueError(f"{operation} snapshot is missing columns: {missing}.")
        indexes = {item: header.index(item) for item in active_columns}
        for row_number, row in enumerate(reader, start=2):
            logical_rows += 1
            if logical_rows > max_rows:
                raise ValueError(f"{operation} snapshot exceeds max_rows={max_rows:,}.")
            if len(row) != len(header):
                raise ValueError(
                    f"{operation} logical row {row_number} has {len(row)} fields; expected exactly {len(header)}."
                )
            drug = identifier(row[indexes[drug_column]], "drug identifier", row_number, seen_drug_spellings)
            gene = identifier(row[indexes[gene_column]], "gene identifier", row_number, seen_gene_spellings)
            source = (
                provenance_value(row[indexes[source_column]], "source value", row_number)
                if source_column is not None
                else None
            )
            evidence = (
                provenance_value(row[indexes[evidence_column]], "evidence value", row_number)
                if evidence_column is not None
                else None
            )
            rows_without_source += int(source is None)
            rows_without_evidence += int(evidence is None)
            pair = (drug, gene)
            if pair not in pairs:
                pairs[pair] = {"sources": [], "evidence": []}
                pair_order.append(pair)
            if source is not None and source not in pairs[pair]["sources"]:
                pairs[pair]["sources"].append(source)
            if evidence is not None and evidence not in pairs[pair]["evidence"]:
                pairs[pair]["evidence"].append(evidence)

    if logical_rows == 0:
        raise ValueError(f"{operation} snapshot contains no data rows.")
    stat_after = os.stat(resolved_path)
    sha_after = stream_sha256()
    if (
        int(stat_after.st_size) != size_bytes
        or int(stat_after.st_mtime_ns) != int(stat_before.st_mtime_ns)
        or sha_after != sha_before
    ):
        raise RuntimeError(f"{operation} snapshot changed while it was being parsed.")
    if expected_sha256 is not None and sha_after != expected_sha256:
        raise RuntimeError(f"{operation} snapshot no longer matches the pinned expected bytes.")

    records = []
    for drug, gene in pair_order:
        provenance = pairs[(drug, gene)]
        records.append(
            {
                "drug": drug,
                "gene": gene,
                "sources": json.dumps(
                    provenance["sources"], ensure_ascii=False, allow_nan=False, separators=(",", ":")
                ),
                "evidence": json.dumps(
                    provenance["evidence"], ensure_ascii=False, allow_nan=False, separators=(",", ":")
                ),
            }
        )
    table = pd.DataFrame.from_records(records, columns=["drug", "gene", "sources", "evidence"])
    table = _standalone_validate_dgidb_table(table)
    accounting = {
        "extension": extension,
        "requested_path": requested_path,
        "resolved_path": resolved_path,
        "size_bytes": size_bytes,
        "logical_data_rows": logical_rows,
        "raw_pairs": logical_rows,
        "canonical_pairs": len(table),
        "duplicate_pairs_removed": logical_rows - len(table),
        "drug_count": int(table["drug"].nunique()),
        "gene_count": int(table["gene"].nunique()),
        "source_value_count": sum(len(values["sources"]) for values in pairs.values()),
        "evidence_value_count": sum(len(values["evidence"]) for values in pairs.values()),
        "rows_without_source": rows_without_source,
        "rows_without_evidence": rows_without_evidence,
        "drug_order": list(dict.fromkeys(table["drug"].tolist())),
        "drug_column": drug_column,
        "gene_column": gene_column,
        "source_column": source_column,
        "evidence_column": evidence_column,
        "max_file_bytes": max_file_bytes,
        "max_rows": max_rows,
    }
    accounting = _standalone_validate_dgidb_accounting(accounting, table)
    artifact_metadata = _standalone_build_dgidb_artifact_metadata(
        table, metadata, accounting, sha_before, size_bytes
    )
    return {
        "table": table,
        "metadata": metadata,
        "accounting": accounting,
        "artifact_metadata": artifact_metadata,
    }


def _standalone_prepare_dgidb_for_universe(
    table,
    metadata,
    accounting,
    artifact_metadata,
    universe_genes,
    min_targets,
    operation,
):
    import pandas as pd

    table, metadata, accounting, artifact_metadata = _standalone_validate_portable_dgidb_resource(
        table, metadata, accounting, artifact_metadata
    )
    if isinstance(min_targets, bool) or not isinstance(min_targets, int) or min_targets < 1:
        raise ValueError(f"{operation} min_targets must be a positive integer.")
    if isinstance(universe_genes, str):
        raise TypeError(f"{operation} universe genes must be an ordered collection.")
    genes = list(universe_genes)
    if not genes:
        raise ValueError(f"{operation} exact tested-gene universe cannot be empty.")
    if any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in genes):
        raise ValueError(f"{operation} exact tested-gene universe contains an invalid identifier.")
    if len(genes) != len(set(genes)):
        raise ValueError(f"{operation} exact tested-gene universe must be unique.")
    universe_set = set(genes)
    source_order = list(accounting["drug_order"])
    targets_before = {
        drug: tuple(table.loc[table["drug"] == drug, "gene"].tolist()) for drug in source_order
    }
    targets_in_universe = {
        drug: tuple(gene for gene in targets_before[drug] if gene in universe_set) for drug in source_order
    }
    retained_sources = tuple(
        drug for drug in source_order if len(targets_in_universe[drug]) >= min_targets
    )
    retained_set = set(retained_sources)
    excluded_sources = tuple(drug for drug in source_order if drug not in retained_set)
    all_targets = {gene for drug in source_order for gene in targets_in_universe[drug]}
    if not all_targets:
        raise ValueError(f"{operation} DGIdb resource has zero HGNC overlap with the exact tested universe.")
    if not retained_sources:
        raise ValueError(
            f"{operation} has no DGIdb drug with at least {min_targets} targets after universe intersection."
        )
    network = pd.DataFrame(
        [(drug, gene) for drug in retained_sources for gene in targets_in_universe[drug]],
        columns=["source", "target"],
    )
    set_accounting = [
        {
            "source": drug,
            "source_order": index,
            "set_size_before_universe": len(targets_before[drug]),
            "set_size_in_universe": len(targets_in_universe[drug]),
            "retained_by_min_targets": drug in retained_set,
        }
        for index, drug in enumerate(source_order)
    ]
    generic_accounting = {
        "extension": accounting["extension"],
        "logical_data_rows": accounting["logical_data_rows"],
        "input_pairs": accounting["raw_pairs"],
        "unique_pairs": accounting["canonical_pairs"],
        "duplicate_pairs_removed": accounting["duplicate_pairs_removed"],
        "source_count_before": len(source_order),
        "source_count_after_min_targets": len(retained_sources),
        "sources_removed_by_min_targets_count": len(excluded_sources),
        "sources_removed_by_min_targets": list(excluded_sources),
        "targets_in_universe_count": len(all_targets),
        "set_accounting": set_accounting,
        "dgidb_artifact_fingerprint_sha256": artifact_metadata["artifact_fingerprint_sha256"],
        "dgidb_canonical_content_fingerprint_sha256": artifact_metadata[
            "canonical_content_fingerprint_sha256"
        ],
        "dgidb_source_value_count": accounting["source_value_count"],
        "dgidb_evidence_value_count": accounting["evidence_value_count"],
        "dgidb_download_url": metadata["download_url"],
        "dgidb_source_license_review": metadata["source_license_review"],
    }
    return {
        "network": network,
        "metadata": {
            "name": metadata["name"],
            "version": metadata["version"],
            "date": metadata["release_date"],
            "organism": metadata["organism"],
            "identifier_namespace": metadata["gene_identifier_namespace"],
            "scope": metadata["scope"],
            "license": metadata["license"],
            "citation": metadata["citation"],
        },
        "sha256": artifact_metadata["raw_file_sha256"],
        "size_bytes": artifact_metadata["raw_size_bytes"],
        "resolved_path": f"typed:{artifact_metadata['artifact_fingerprint_sha256']}",
        "accounting": generic_accounting,
        "source_order": tuple(source_order),
        "targets_before": targets_before,
        "targets_in_universe": targets_in_universe,
        "retained_sources": retained_sources,
        "excluded_sources": excluded_sources,
    }


def _standalone_build_dgidb_resource_summary(payload, openbio_version):
    import json
    import platform
    from importlib import metadata as importlib_metadata

    table, metadata, accounting, artifact = _standalone_validate_portable_dgidb_resource(
        payload["table"],
        payload["metadata"],
        payload["accounting"],
        payload["artifact_metadata"],
    )
    del table

    def package_version(name):
        try:
            return importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            return "not-installed"

    bounded_accounting = {key: value for key, value in accounting.items() if key != "drug_order"}
    bounded_accounting["drug_order_preview"] = accounting["drug_order"][:100]
    bounded_accounting["drug_order_preview_truncated"] = len(accounting["drug_order"]) > 100
    key_results = {
        "resource": {
            "metadata": metadata,
            "raw_file_sha256": artifact["raw_file_sha256"],
            "canonical_content_fingerprint_sha256": artifact[
                "canonical_content_fingerprint_sha256"
            ],
            "artifact_fingerprint_sha256": artifact["artifact_fingerprint_sha256"],
            "artifact_schema_version": artifact["schema_version"],
            "producer_schema": artifact["producer_schema"],
            "accounting": bounded_accounting,
        }
    }
    parameters = {
        "resource_path": accounting["requested_path"],
        "drug_column": accounting["drug_column"],
        "gene_column": accounting["gene_column"],
        "source_column": accounting["source_column"],
        "evidence_column": accounting["evidence_column"],
        "max_file_bytes": accounting["max_file_bytes"],
        "max_rows": accounting["max_rows"],
        "canonical_policy": {
            "trim_surrounding_identifier_whitespace": True,
            "normalization_collisions_rejected": True,
            "exact_duplicate_pairs_collapsed": True,
            "source_evidence_values_preserved": True,
            "online_fetch": False,
        },
    }
    warnings = [
        "Source-license review and biological suitability are caller attestations; SHA-256 proves exact bytes, not legal clearance or scientific validity.",
    ]
    if metadata["license"].casefold() in {"unknown", "none", "n/a", "na", "unreviewed"}:
        warnings.append(
            "The caller-declared license is a placeholder or explicitly unreviewed; downstream use and redistribution require an independent source-license review."
        )
    if not metadata["source_license_review"].casefold().startswith("reviewed"):
        warnings.append(
            "source_license_review does not declare a completed review. This caller attestation is retained for audit, not treated as programmatic legal clearance."
        )
    if accounting["duplicate_pairs_removed"]:
        warnings.append(
            f"Collapsed {accounting['duplicate_pairs_removed']:,} exact duplicate drug-gene rows while preserving unique source/evidence values."
        )
    if accounting["source_column"] is None:
        warnings.append("No source column was retained; contributing DGIdb sources cannot be stratified.")
    if accounting["evidence_column"] is None:
        warnings.append("No evidence column was retained; interaction evidence categories cannot be stratified.")
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellDGIdbAnnotation",
        "methods": (
            "Loaded one caller-described local DGIdb CSV/TSV snapshot without network access; validated the exact "
            "schema and ten-field release/license metadata, streamed SHA-256 before and after strict logical-row "
            "parsing, trimmed only surrounding identifier whitespace, rejected normalization collisions, and "
            "collapsed exact drug-gene duplicates while preserving source/evidence values."
        ),
        "results": (
            f"Validated DGIdb resource {metadata['version']!r} ({metadata['release_date']}) with "
            f"{accounting['canonical_pairs']:,} canonical drug-gene pairs, {accounting['drug_count']:,} drugs, "
            f"and {accounting['gene_count']:,} unique HGNC genes."
        ),
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": [
            "DGIdb edges are aggregated interaction claims and do not establish efficacy, therapeutic direction, dose, safety, population relevance, or a treatment recommendation.",
            "The DGIdb software license does not automatically license every contributing source; users must independently review and honor source-specific terms rather than relying on caller metadata as legal clearance.",
            "Only Homo sapiens HGNC symbols are accepted; this loader does not map aliases, normalize drug names, or infer clinical direction.",
        ],
        "references": [
            {
                "citation": "Cannon M, et al. DGIdb 5.0. Nucleic Acids Research. 2024;52:D1227-D1235.",
                "url": "https://doi.org/10.1093/nar/gkad1040",
                "doi": "10.1093/nar/gkad1040",
                "kind": "resource",
            },
            {
                "citation": "Heumos L, et al. pertpy: an end-to-end framework for perturbation analysis. Nature Methods. 2025.",
                "url": "https://doi.org/10.1038/s41592-025-02909-7",
                "doi": "10.1038/s41592-025-02909-7",
                "kind": "software",
            },
            {
                "citation": f"{metadata['name']} {metadata['version']}: {metadata['citation']}",
                "url": metadata["download_url"],
                "doi": None,
                "kind": "resource_snapshot",
            },
        ],
        "software_versions": {
            "python": platform.python_version(),
            "openbio-singlecell": str(openbio_version),
            "pandas": package_version("pandas"),
        },
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


class DGIdbResource:
    """Immutable-by-interface, tamper-evident local DGIdb resource artifact."""

    __slots__ = ("_accounting", "_artifact_metadata", "_metadata", "_table")
    artifact_type = DGIDB_ARTIFACT_TYPE

    def __init__(
        self,
        *,
        table: Any,
        metadata: Mapping[str, Any],
        accounting: Mapping[str, Any],
        artifact_metadata: Mapping[str, Any],
    ) -> None:
        object.__setattr__(self, "_table", table.copy(deep=True))
        object.__setattr__(self, "_metadata", copy.deepcopy(dict(metadata)))
        object.__setattr__(self, "_accounting", copy.deepcopy(dict(accounting)))
        object.__setattr__(self, "_artifact_metadata", copy.deepcopy(dict(artifact_metadata)))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("DGIdbResource is immutable; load a new validated artifact instead.")

    @property
    def fingerprint(self) -> str:
        return str(self._artifact_metadata["artifact_fingerprint_sha256"])

    @property
    def metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._metadata)

    @property
    def accounting(self) -> dict[str, Any]:
        return copy.deepcopy(self._accounting)

    @property
    def artifact_metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._artifact_metadata)

    def table(self) -> DataFrame:
        return self._table.copy(deep=True)


def validate_dgidb_resource(
    resource: Any,
) -> tuple[DataFrame, dict[str, Any], dict[str, Any], dict[str, Any]]:
    if type(resource) is not DGIdbResource:
        raise TypeError(
            "Drug analysis requires an exact OPENBIO_DGIDB_RESOURCE artifact from Load DGIdb Resource."
        )
    if resource.artifact_type != DGIDB_ARTIFACT_TYPE:
        raise ValueError("DGIdb resource type identity is invalid.")
    table, metadata, accounting, artifact = _standalone_validate_portable_dgidb_resource(
        resource.table(),
        resource.metadata,
        resource.accounting,
        resource.artifact_metadata,
    )
    if resource.fingerprint != artifact["artifact_fingerprint_sha256"]:
        raise ValueError("DGIdb resource fingerprint identity is invalid.")
    return table, metadata, accounting, artifact


def load_dgidb_resource(
    path: str,
    *,
    requested_path: str,
    resource_metadata_json: str,
    drug_column: str = "drug_claim_name",
    gene_column: str = "gene_claim_name",
    source_column: str = "interaction_claim_source",
    evidence_column: str = "interaction_types",
    max_file_bytes: int = 536_870_912,
    max_rows: int = 2_000_000,
    expected_sha256: str | None = None,
    openbio_version: str = "not-installed",
) -> tuple[DGIdbResource, dict[str, Any]]:
    payload = _standalone_load_dgidb_snapshot(
        path,
        requested_path,
        resource_metadata_json,
        drug_column=drug_column,
        gene_column=gene_column,
        source_column=source_column,
        evidence_column=evidence_column,
        max_file_bytes=max_file_bytes,
        max_rows=max_rows,
        expected_sha256=expected_sha256,
    )
    resource = DGIdbResource(**payload)
    validate_dgidb_resource(resource)
    summary = _standalone_build_dgidb_resource_summary(payload, openbio_version)
    return resource, summary


def dgidb_resource_cache_fingerprint(
    path: str,
    identity: tuple[Any, ...],
    *,
    max_file_bytes: int = 536_870_912,
) -> tuple[Any, ...]:
    """Bind caching to guarded snapshot bytes and reject mutation during fingerprinting."""
    import os

    if isinstance(max_file_bytes, bool) or not isinstance(max_file_bytes, int) or max_file_bytes < 1:
        raise ValueError("DGIdb cache max_file_bytes must be a positive integer.")
    stat_before = os.stat(path)
    size_bytes = int(stat_before.st_size)
    if size_bytes < 1:
        raise ValueError("DGIdb cache snapshot is empty.")
    if size_bytes > max_file_bytes:
        raise ValueError(
            f"DGIdb cache snapshot is {size_bytes:,} bytes, exceeding max_file_bytes={max_file_bytes:,}."
        )
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    stat_after = os.stat(path)
    if (
        int(stat_after.st_size) != size_bytes
        or int(stat_after.st_mtime_ns) != int(stat_before.st_mtime_ns)
    ):
        raise RuntimeError("DGIdb cache snapshot changed while it was being fingerprinted.")
    return ("openbio-dgidb-resource-v1", *identity, size_bytes, digest.hexdigest())


def dgidb_portable_validation_code(*, include_prepare: bool = False) -> str:
    functions = [
        _standalone_dgidb_json_sha256,
        _standalone_validate_dgidb_metadata,
        _standalone_validate_dgidb_table,
        _standalone_dgidb_table_fingerprint,
        _standalone_validate_dgidb_accounting,
        _standalone_build_dgidb_artifact_metadata,
        _standalone_validate_portable_dgidb_resource,
    ]
    if include_prepare:
        functions.append(_standalone_prepare_dgidb_for_universe)
    return "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)


def dgidb_resource_code(**parameters: Any) -> str:
    loader_parameters = dict(parameters)
    openbio_version = loader_parameters.pop("openbio_version", "not-installed")
    implementations = "\n\n".join(
        (
            dgidb_portable_validation_code(),
            textwrap.dedent(inspect.getsource(_standalone_load_dgidb_snapshot)).strip(),
            textwrap.dedent(inspect.getsource(_standalone_build_dgidb_resource_summary)).strip(),
        )
    )
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in loader_parameters.items())
    return f'''from __future__ import annotations

{implementations}


def load_dgidb_resource(resource_path):
    payload = _standalone_load_dgidb_snapshot(
        resource_path,
        {rendered},
    )
    summary = _standalone_build_dgidb_resource_summary(payload, {openbio_version!r})
    return payload["table"].copy(deep=True), summary
'''


__all__ = [
    "DGIDB_ARTIFACT_SCHEMA_VERSION",
    "DGIDB_ARTIFACT_TYPE",
    "DGIDB_METADATA_FIELDS",
    "DGIDB_PRODUCER_NODE_ID",
    "DGIDB_PRODUCER_SCHEMA",
    "DGIDB_RESOURCE_COLUMNS",
    "DGIDB_RESOURCE_EXTENSIONS",
    "DGIdbResource",
    "dgidb_portable_validation_code",
    "dgidb_resource_cache_fingerprint",
    "dgidb_resource_code",
    "load_dgidb_resource",
    "validate_dgidb_resource",
]
