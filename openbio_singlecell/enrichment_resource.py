from __future__ import annotations

import inspect
import textwrap
from typing import Any

GENE_SET_RESOURCE_EXTENSIONS = (".csv", ".tsv", ".gmt")
GENE_SET_RESOURCE_METADATA_FIELDS = (
    "name",
    "version",
    "date",
    "organism",
    "identifier_namespace",
    "scope",
    "license",
    "citation",
)


def _standalone_load_gene_set_resource(
    path,
    resource_metadata_json,
    source_column,
    target_column,
    universe_genes,
    min_targets,
    *,
    expected_sha256=None,
    operation="Gene-set resource",
):
    """Load one provenance-complete local gene-set resource.

    The returned mapping is owned by the caller.  Its network has canonical
    ``source``/``target`` columns and contains only sources that meet
    ``min_targets`` after exact intersection with ``universe_genes``.
    """
    import csv
    import hashlib
    import json
    import os
    from collections.abc import Sequence

    import pandas as pd

    metadata_fields = (
        "name",
        "version",
        "date",
        "organism",
        "identifier_namespace",
        "scope",
        "license",
        "citation",
    )
    supported_extensions = (".csv", ".tsv", ".gmt")

    if not isinstance(operation, str) or not operation.strip():
        raise ValueError("Gene-set resource operation label must be a nonblank string.")
    operation = operation.strip()

    def required_identifier(value, *, name):
        if not isinstance(value, str):
            raise TypeError(f"{operation} {name} must be a string.")
        if not value.strip():
            raise ValueError(f"{operation} {name} cannot be empty or whitespace-only.")
        if value != value.strip():
            raise ValueError(f"{operation} {name} cannot contain surrounding whitespace.")
        return value

    if isinstance(min_targets, bool) or not isinstance(min_targets, int) or min_targets < 1:
        raise ValueError(f"{operation} min_targets must be a positive integer.")
    source_column = required_identifier(source_column, name="source_column")
    target_column = required_identifier(target_column, name="target_column")
    if source_column == target_column:
        raise ValueError(f"{operation} source_column and target_column must be different.")

    if isinstance(universe_genes, (str, bytes)) or not isinstance(universe_genes, Sequence):
        raise TypeError(f"{operation} universe_genes must be an ordered sequence of identifiers.")
    ordered_universe = [required_identifier(value, name="tested-universe identifier") for value in universe_genes]
    if not ordered_universe:
        raise ValueError(f"{operation} tested-gene universe cannot be empty.")
    if len(ordered_universe) != len(set(ordered_universe)):
        raise ValueError(f"{operation} tested-gene universe identifiers must be unique.")
    universe_set = set(ordered_universe)

    if not isinstance(resource_metadata_json, str):
        raise TypeError(f"{operation} resource_metadata_json must be a string.")
    if len(resource_metadata_json.encode("utf-8")) > 65_536:
        raise ValueError(f"{operation} resource metadata exceeds 65,536 UTF-8 bytes.")

    def reject_constant(value):
        raise ValueError(f"{operation} resource metadata contains non-standard JSON constant {value!r}.")

    def object_from_pairs(pairs):
        result = {}
        for key, value in pairs:
            if not isinstance(key, str) or not key.strip():
                raise ValueError(f"{operation} resource metadata keys must be nonblank strings.")
            if key != key.strip():
                raise ValueError(f"{operation} resource metadata keys cannot contain surrounding whitespace.")
            if key in result:
                raise ValueError(f"{operation} resource metadata contains duplicate key {key!r}.")
            result[key] = value
        return result

    try:
        metadata = json.loads(
            resource_metadata_json,
            object_pairs_hook=object_from_pairs,
            parse_constant=reject_constant,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"{operation} resource metadata is invalid JSON ({exc.msg}).") from exc
    required_metadata = set(metadata_fields)
    if not isinstance(metadata, dict) or set(metadata) != required_metadata:
        actual = set(metadata) if isinstance(metadata, dict) else set()
        raise ValueError(
            f"{operation} resource metadata must contain exactly the required fields; "
            f"missing={sorted(required_metadata - actual)}, unknown={sorted(actual - required_metadata)}."
        )
    for field in metadata_fields:
        value = metadata[field]
        if not isinstance(value, str):
            raise TypeError(f"{operation} resource metadata field {field!r} must be a string.")
        if not value.strip():
            raise ValueError(f"{operation} resource metadata field {field!r} cannot be blank.")
        if value != value.strip():
            raise ValueError(f"{operation} resource metadata field {field!r} cannot contain surrounding whitespace.")

    if not isinstance(path, str) or not path.strip():
        raise ValueError(f"{operation} resource path cannot be empty.")
    if path != path.strip():
        raise ValueError(f"{operation} resource path cannot contain surrounding whitespace.")
    resolved_path = os.path.realpath(path)
    if not os.path.isfile(resolved_path):
        raise FileNotFoundError(f"{operation} resource file not found: {resolved_path}")
    extension = os.path.splitext(resolved_path)[1].lower()
    if extension not in supported_extensions:
        raise ValueError(
            f"{operation} resource extension must be one of {list(supported_extensions)}; received {extension!r}."
        )
    size_bytes = int(os.path.getsize(resolved_path))
    if size_bytes > 512 * 1024 * 1024:
        raise ValueError(f"{operation} resource exceeds the 512 MiB safety limit.")

    def file_sha256():
        digest = hashlib.sha256()
        with open(resolved_path, "rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    first_sha256 = file_sha256()
    if expected_sha256 is not None:
        if (
            not isinstance(expected_sha256, str)
            or len(expected_sha256) != 64
            or any(character not in "0123456789abcdef" for character in expected_sha256)
        ):
            raise ValueError(f"{operation} expected resource SHA-256 is invalid.")
        if first_sha256 != expected_sha256:
            raise ValueError(
                f"{operation} resource fingerprint changed: expected {expected_sha256}, observed {first_sha256}."
            )

    pairs = []
    logical_rows = 0
    if extension in {".csv", ".tsv"}:
        delimiter = "," if extension == ".csv" else "\t"
        try:
            with open(resolved_path, encoding="utf-8-sig", newline="") as stream:
                reader = csv.reader(stream, delimiter=delimiter)
                try:
                    header = next(reader)
                except StopIteration as exc:
                    raise ValueError(f"{operation} resource is empty.") from exc
                if not header:
                    raise ValueError(f"{operation} resource header is empty.")
                for column in header:
                    required_identifier(column, name="resource column name")
                if len(header) != len(set(header)):
                    raise ValueError(f"{operation} resource contains duplicate column names.")
                missing = [column for column in (source_column, target_column) if column not in header]
                if missing:
                    raise ValueError(f"{operation} resource is missing source/target columns: {missing}.")
                source_index = header.index(source_column)
                target_index = header.index(target_column)
                for logical_rows, row in enumerate(reader, start=1):
                    if logical_rows > 2_000_000:
                        raise ValueError(f"{operation} resource exceeds the 2,000,000-row safety limit.")
                    if len(row) != len(header):
                        raise ValueError(
                            f"{operation} resource logical data row {logical_rows} has {len(row)} fields; "
                            f"expected exactly {len(header)} from the header."
                        )
                    source = required_identifier(row[source_index], name="resource source identifier")
                    target = required_identifier(row[target_index], name="resource target identifier")
                    pairs.append((source, target))
        except UnicodeError as exc:
            raise ValueError(f"{operation} resource must be valid UTF-8 (optional BOM is accepted).") from exc
    else:
        try:
            with open(resolved_path, encoding="utf-8-sig", newline="") as stream:
                reader = csv.reader(stream, delimiter="\t")
                for logical_rows, row in enumerate(reader, start=1):
                    if logical_rows > 2_000_000:
                        raise ValueError(f"{operation} GMT exceeds the 2,000,000-set safety limit.")
                    if len(row) < 3:
                        raise ValueError(
                            f"{operation} GMT logical row {logical_rows} must contain a source, description, "
                            "and at least one target."
                        )
                    source = required_identifier(row[0], name="GMT source identifier")
                    for target_value in row[2:]:
                        target = required_identifier(target_value, name="GMT target identifier")
                        pairs.append((source, target))
                        if len(pairs) > 2_000_000:
                            raise ValueError(f"{operation} GMT exceeds the 2,000,000-edge safety limit.")
        except UnicodeError as exc:
            raise ValueError(f"{operation} GMT must be valid UTF-8 (optional BOM is accepted).") from exc

    second_sha256 = file_sha256()
    if first_sha256 != second_sha256:
        raise ValueError(f"{operation} resource changed while it was being parsed; rerun with stable bytes.")
    if int(os.path.getsize(resolved_path)) != size_bytes:
        raise ValueError(f"{operation} resource size changed while it was being parsed.")
    if not pairs:
        raise ValueError(f"{operation} resource contains no source-target pairs.")

    unique_pairs = list(dict.fromkeys(pairs))
    source_order = list(dict.fromkeys(source for source, _ in unique_pairs))
    targets_before = {source: [] for source in source_order}
    for source, target in unique_pairs:
        targets_before[source].append(target)
    targets_in_universe = {
        source: tuple(target for target in targets_before[source] if target in universe_set) for source in source_order
    }
    retained_sources = tuple(source for source in source_order if len(targets_in_universe[source]) >= min_targets)
    retained_source_set = set(retained_sources)
    excluded_sources = tuple(source for source in source_order if source not in retained_source_set)
    all_targets_in_universe = {target for source in source_order for target in targets_in_universe[source]}
    if not all_targets_in_universe:
        raise ValueError(f"{operation} resource has zero target overlap with the exact tested-gene universe.")
    if not retained_sources:
        raise ValueError(
            f"{operation} has no resource set with at least {min_targets} targets after universe intersection."
        )
    retained_pairs = [(source, target) for source in retained_sources for target in targets_in_universe[source]]
    network = pd.DataFrame(retained_pairs, columns=["source", "target"])
    if network.empty or bool(network.duplicated().any()):
        raise RuntimeError(f"{operation} failed to build a nonempty unique canonical network.")

    set_accounting = [
        {
            "source": source,
            "source_order": index,
            "set_size_before_universe": len(targets_before[source]),
            "set_size_in_universe": len(targets_in_universe[source]),
            "retained_by_min_targets": source in retained_source_set,
        }
        for index, source in enumerate(source_order)
    ]
    accounting = {
        "extension": extension,
        "logical_data_rows": logical_rows,
        "input_pairs": len(pairs),
        "unique_pairs": len(unique_pairs),
        "duplicate_pairs_removed": len(pairs) - len(unique_pairs),
        "source_count_before": len(source_order),
        "source_count_after_min_targets": len(retained_sources),
        "sources_removed_by_min_targets_count": len(excluded_sources),
        "sources_removed_by_min_targets": list(excluded_sources),
        "targets_before_universe_count": len({target for _, target in unique_pairs}),
        "targets_in_universe_count": len(all_targets_in_universe),
        "retained_target_rows_in_universe": len(retained_pairs),
        "universe_count": len(ordered_universe),
        "min_targets_after_universe_intersection": min_targets,
        "set_accounting": set_accounting,
    }
    json.dumps({"metadata": metadata, "accounting": accounting}, ensure_ascii=False, allow_nan=False)
    return {
        "network": network,
        "metadata": dict(metadata),
        "resolved_path": resolved_path,
        "size_bytes": size_bytes,
        "sha256": second_sha256,
        "source_order": tuple(source_order),
        "retained_sources": retained_sources,
        "excluded_sources": excluded_sources,
        "targets_before": {source: tuple(targets_before[source]) for source in source_order},
        "targets_in_universe": targets_in_universe,
        "accounting": accounting,
    }


def load_gene_set_resource(
    path: str,
    resource_metadata_json: str,
    source_column: str,
    target_column: str,
    universe_genes: list[str],
    min_targets: int,
    *,
    expected_sha256: str | None = None,
    operation: str = "Gene-set resource",
) -> dict[str, Any]:
    return _standalone_load_gene_set_resource(
        path,
        resource_metadata_json,
        source_column,
        target_column,
        universe_genes,
        min_targets,
        expected_sha256=expected_sha256,
        operation=operation,
    )


def gene_set_resource_code() -> str:
    """Return the standalone loader source for generated equivalent code."""
    return textwrap.dedent(inspect.getsource(_standalone_load_gene_set_resource)).strip()


__all__ = [
    "GENE_SET_RESOURCE_EXTENSIONS",
    "GENE_SET_RESOURCE_METADATA_FIELDS",
    "gene_set_resource_code",
    "load_gene_set_resource",
]
