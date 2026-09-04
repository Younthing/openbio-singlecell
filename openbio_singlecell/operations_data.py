from __future__ import annotations

import copy
import gzip
import re
import time
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from textwrap import dedent
from typing import Any

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata
from .expression_source import _SNAPSHOT_SPEC, ExpressionSource
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_file_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
)
from .worker_protocol import JSONValue, OperationContext, ProtocolError, register_operation

_GTF_ATTRIBUTE_PATTERN = re.compile(r'([^\s;]+)\s+"([^"]*)"')
_ENSEMBL_VERSION_PATTERN = re.compile(r"\.\d+(?=_PAR_Y$|$)")
CANONICAL_COUNTS_LAYER = "counts"
ORIGINAL_FEATURE_NAMES_COLUMN = "feature_names_before_gtf"
GTF_MATCH_COLUMN = "gtf_match_type"
GTF_MAPPED_COLUMN = "gtf_mapped"
SNAPSHOT_EXPRESSION_SOURCE = _SNAPSHOT_SPEC
ANNDATA_REFERENCE = AnalysisReference(
    citation=(
        "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data "
        "matrices. Journal of Open Source Software. 2024;9(101):4371."
    ),
    doi="10.21105/joss.04371",
    url="https://doi.org/10.21105/joss.04371",
    kind="software",
)
PANDAS_REFERENCE = AnalysisReference(
    citation="The pandas development team. pandas software documentation and citation guidance.",
    url="https://pandas.pydata.org/about/citing.html",
    kind="software",
)
GENCODE_REFERENCE = AnalysisReference(
    citation=(
        "Frankish A et al. GENCODE 2025: reference gene annotation for human and mouse. Nucleic Acids Research. 2025."
    ),
    doi="10.1093/nar/gkae1078",
    url="https://doi.org/10.1093/nar/gkae1078",
    kind="method",
)
ENSEMBL_REFERENCE = AnalysisReference(
    citation="Martin FJ et al. Ensembl 2025. Nucleic Acids Research. 2025.",
    doi="10.1093/nar/gkae1071",
    url="https://doi.org/10.1093/nar/gkae1071",
    kind="method",
)


def _required_string(value: JSONValue, *, name: str) -> str:
    if not isinstance(value, str):
        raise ProtocolError(f"{name} must be a string.")
    return value


def _required_bool(value: JSONValue, *, name: str) -> bool:
    if not isinstance(value, bool):
        raise ProtocolError(f"{name} must be a boolean.")
    return value


def _comma_separated_values(values: str) -> tuple[str, ...]:
    parsed = tuple(dict.fromkeys(value.strip() for value in values.split(",") if value.strip()))
    if not parsed:
        raise ValueError("Observation values must contain at least one non-empty value.")
    return parsed


def _stable_gene_id(value: str) -> str:
    return _ENSEMBL_VERSION_PATTERN.sub("", value)


@dataclass(frozen=True, slots=True)
class _GTFGeneMapping:
    exact: dict[str, str]
    normalized: dict[str, str]
    headers: dict[str, str]
    used_fallback_records: bool
    total_data_lines: int
    parsed_records: int
    unique_gene_ids: int
    malformed_lines: int
    missing_attribute_lines: int
    repeated_records: int


def _gtf_header_item(line: str) -> tuple[str, str] | None:
    text = line.lstrip("#").strip()
    if not text:
        return None
    if ":" in text:
        key, value = text.split(":", 1)
    elif line.startswith("#!") and " " in text:
        key, value = text.split(None, 1)
    else:
        return None
    key = re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")[:128]
    value = value.strip()[:4096]
    return (key, value) if key and value else None


def _read_gtf_gene_mapping(path: str) -> _GTFGeneMapping:
    headers: dict[str, str] = {}
    gene_exact_names: dict[str, set[str]] = {}
    gene_normalized_names: dict[str, set[str]] = {}
    fallback_exact_names: dict[str, set[str]] = {}
    fallback_normalized_names: dict[str, set[str]] = {}
    gene_seen_records: set[tuple[str, str]] = set()
    fallback_seen_records: set[tuple[str, str]] = set()
    gene_record_count = 0
    fallback_record_count = 0
    gene_repeated_records = 0
    fallback_repeated_records = 0
    total_data_lines = 0
    malformed_lines = 0
    missing_attribute_lines = 0
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                item = _gtf_header_item(line)
                if item is not None:
                    headers.setdefault(*item)
                continue
            total_data_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                malformed_lines += 1
                continue
            attributes = dict(_GTF_ATTRIBUTE_PATTERN.findall(fields[8]))
            gene_id = attributes.get("gene_id")
            gene_name = attributes.get("gene_name")
            if not gene_id or not gene_name:
                missing_attribute_lines += 1
                continue
            record = (gene_id.strip(), gene_name.strip())
            if not record[0] or not record[1]:
                missing_attribute_lines += 1
                continue
            if fields[2] == "gene":
                gene_record_count += 1
                if record in gene_seen_records:
                    gene_repeated_records += 1
                gene_seen_records.add(record)
                gene_exact_names.setdefault(record[0], set()).add(record[1])
                gene_normalized_names.setdefault(_stable_gene_id(record[0]), set()).add(record[1])
            else:
                fallback_record_count += 1
                if record in fallback_seen_records:
                    fallback_repeated_records += 1
                fallback_seen_records.add(record)
                fallback_exact_names.setdefault(record[0], set()).add(record[1])
                fallback_normalized_names.setdefault(_stable_gene_id(record[0]), set()).add(record[1])

    used_fallback = gene_record_count == 0 and fallback_record_count > 0
    if gene_record_count:
        exact_names = gene_exact_names
        normalized_names = gene_normalized_names
        parsed_records = gene_record_count
        repeated_records = gene_repeated_records
    else:
        exact_names = fallback_exact_names
        normalized_names = fallback_normalized_names
        parsed_records = fallback_record_count
        repeated_records = fallback_repeated_records
    if not parsed_records:
        raise ValueError("The GTF file does not contain gene entries with gene_id and gene_name attributes.")

    conflicts = [
        f"{gene_id!r}: {sorted(names)!r}"
        for gene_id, names in [*exact_names.items(), *normalized_names.items()]
        if len(names) > 1
    ]
    if conflicts:
        preview = "; ".join(conflicts[:8])
        suffix = f"; plus {len(conflicts) - 8} more" if len(conflicts) > 8 else ""
        raise ValueError(f"The GTF contains ambiguous gene_id-to-gene_name mappings: {preview}{suffix}")

    return _GTFGeneMapping(
        exact={gene_id: next(iter(names)) for gene_id, names in exact_names.items()},
        normalized={gene_id: next(iter(names)) for gene_id, names in normalized_names.items()},
        headers=headers,
        used_fallback_records=used_fallback,
        total_data_lines=total_data_lines,
        parsed_records=parsed_records,
        unique_gene_ids=len(exact_names),
        malformed_lines=malformed_lines,
        missing_attribute_lines=missing_attribute_lines,
        repeated_records=repeated_records,
    )


def _matrix_description(matrix: Any) -> dict[str, Any]:
    science = dependencies.require_scientific_dependencies()
    storage = f"sparse_{matrix.format}" if science.sparse.issparse(matrix) else "dense"
    return {
        "shape": [int(value) for value in matrix.shape],
        "dtype": str(matrix.dtype),
        "storage": storage,
    }


def _display_expression_source(expression: ExpressionSource) -> str:
    return f"layer:{expression.layer_name}" if expression.kind == "layer" else expression.kind


def _snapshot_expression_code(expression: ExpressionSource, overwrite_existing: bool) -> str:
    source_layer = expression.layer_name
    return dedent(
        f"""
        def create_raw_snapshot(adata):
            source_kind = {expression.kind!r}
            source_layer = {source_layer!r}
            if adata.isbacked:
                raise ValueError("Raw snapshot creation requires an in-memory AnnData; call adata.to_memory() first.")
            if source_kind == "layer" and source_layer not in adata.layers:
                raise ValueError(f"Raw snapshot source layer not found: {{source_layer!r}}")
            matrix = adata.X if source_kind == "X" else adata.layers[source_layer]
            existing = []
            if adata.raw is not None:
                existing.append("raw")
            if "counts" in adata.layers:
                existing.append("layers['counts']")
            if existing and not {overwrite_existing!r}:
                raise ValueError(f"Raw snapshot destinations already exist: {{', '.join(existing)}}")
            output = adata.copy()
            output.layers["counts"] = matrix.copy()
            current_x = output.X
            try:
                output.X = matrix.copy()
                output.raw = output
            finally:
                output.X = current_x
            return output
        """
    )


def _raw_snapshot_to_anndata_code() -> str:
    return dedent(
        """
        import copy


        def raw_snapshot_to_anndata(adata):
            if adata.raw is None:
                raise ValueError("Raw Snapshot to AnnData requires adata.raw; create a Raw snapshot first.")
            if adata.n_obs == 0:
                raise ValueError("Raw Snapshot to AnnData requires at least one observation; received zero cells.")
            if adata.raw.n_vars == 0:
                raise ValueError("Raw Snapshot to AnnData requires at least one Raw feature; adata.raw has zero variables.")
            metadata = copy.deepcopy(adata.uns.get("openbio_singlecell"))
            output = adata.raw.to_adata()
            if output.layers or output.raw is not None:
                raise RuntimeError("AnnData Raw materialization unexpectedly retained layers or nested Raw storage.")
            output.obsm.clear()
            output.obsp.clear()
            output.varm.clear()
            output.varp.clear()
            output.uns.clear()
            if metadata is not None:
                output.uns["openbio_singlecell"] = metadata
            return output
        """
    )


def _bounded_value_counts(values: Any, limit: int = 64) -> dict[str, Any]:
    missing_mask = values.isna()
    counts = values.loc[~missing_mask].astype("object").value_counts(dropna=True, sort=True)
    entries = []
    for value, count in counts.items():
        if hasattr(value, "item"):
            try:
                value = value.item()
            except (TypeError, ValueError):
                pass
        entries.append(
            {
                "value": value,
                "value_type": f"{type(value).__module__}.{type(value).__qualname__}",
                "is_missing": False,
                "count": int(count),
            }
        )
    missing_count = int(missing_mask.sum())
    if missing_count:
        entries.append({"value": None, "value_type": None, "is_missing": True, "count": missing_count})
    entries.sort(
        key=lambda item: (
            -item["count"],
            item["is_missing"],
            item["value_type"] or "",
            repr(item["value"]),
        )
    )
    return {
        "items": entries[:limit],
        "total_categories": len(entries),
        "truncated": len(entries) > limit,
    }


def _equal_annotation_values(left: Any, right: Any) -> bool:
    try:
        result = left == right
        return bool(result) if isinstance(result, bool) or type(result).__module__.startswith("numpy") else False
    except (TypeError, ValueError):
        return False


def _string_label_collision(values: Any, *, label: str) -> None:
    science = dependencies.require_scientific_dependencies()
    rendered: dict[str, list[Any]] = {}
    for value in values:
        if bool(science.pd.isna(value)):
            continue
        key = str(value)
        representatives = rendered.setdefault(key, [])
        if not any(_equal_annotation_values(value, existing) for existing in representatives):
            representatives.append(value)
    collisions = [key for key, representatives in rendered.items() if len(representatives) > 1]
    if collisions:
        raise ValueError(f"{label} contains distinct values that collide after string conversion: {collisions[:8]!r}")


def _subset_observations_code(
    column: str,
    selected_values: tuple[str, ...],
    invert: bool,
    missing_policy: str,
) -> str:
    return dedent(
        f"""
        import pandas as pd


        def _equal_values(left, right):
            try:
                result = left == right
                return bool(result) if isinstance(result, bool) or type(result).__module__.startswith("numpy") else False
            except (TypeError, ValueError):
                return False


        def subset_observations(adata):
            column = {column!r}.strip()
            selected_values = tuple({selected_values!r})
            invert = {invert!r}
            missing_policy = {missing_policy!r}
            if adata.isbacked:
                raise ValueError("Observation subsetting requires an in-memory AnnData; call adata.to_memory() first.")
            if not column:
                raise ValueError("Observation annotation name cannot be empty.")
            if column not in adata.obs:
                raise ValueError(f"Observation annotation not found in obs: {{column!r}}")
            if not selected_values:
                raise ValueError("Observation values must contain at least one non-empty value.")
            if missing_policy not in {{"exclude", "include", "error"}}:
                raise ValueError(f"Unsupported missing_policy: {{missing_policy!r}}")
            series = adata.obs[column]
            rendered = {{}}
            for value in series:
                if pd.isna(value):
                    continue
                key = str(value)
                representatives = rendered.setdefault(key, [])
                if not any(_equal_values(value, existing) for existing in representatives):
                    representatives.append(value)
            collisions = [key for key, representatives in rendered.items() if len(representatives) > 1]
            if collisions:
                raise ValueError(
                    "Observation annotation contains distinct values that collide after string conversion: "
                    f"{{collisions[:8]!r}}"
                )
            not_missing = series.notna()
            missing_count = int((~not_missing).sum())
            if missing_policy == "error" and missing_count:
                raise ValueError(f"Observation annotation contains {{missing_count}} missing value(s).")
            rendered_series = series.astype("string")
            matched = not_missing & rendered_series.isin(selected_values)
            mask = matched if not invert else not_missing & ~matched
            if missing_policy == "include":
                mask = mask | ~not_missing
            return adata[mask.to_numpy(), :].copy()
        """
    )


def _provenance_token(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _provenance_token(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_provenance_token(item) for item in value]
    if not isinstance(value, (str, bytes, bytearray)) and hasattr(value, "tolist"):
        try:
            return _provenance_token(value.tolist())
        except (TypeError, ValueError):
            pass
    if hasattr(value, "item"):
        try:
            return _provenance_token(value.item())
        except (TypeError, ValueError):
            pass
    return value


def _source_provenance(adata: Any) -> Mapping[str, Any] | None:
    metadata = adata.uns.get("openbio_singlecell")
    if not isinstance(metadata, Mapping):
        return None
    source = metadata.get("source")
    return source if isinstance(source, Mapping) and source else None


def _verify_shared_provenance(adata: Any, subset_adata: Any) -> str:
    target_source = _source_provenance(adata)
    subset_source = _source_provenance(subset_adata)
    if target_source is None or subset_source is None:
        return "not_available"
    if _provenance_token(target_source) != _provenance_token(subset_source):
        raise ValueError(
            "The inputs carry different OpenBio source provenance; annotation transfer requires a common source."
        )
    return "verified_equal"


def _annotation_provenance(adata: Any, column: str) -> Mapping[str, Any] | None:
    metadata = adata.uns.get("openbio_singlecell")
    if not isinstance(metadata, Mapping):
        return None
    annotations = metadata.get("annotations")
    if not isinstance(annotations, Mapping):
        return None
    provenance = annotations.get(column)
    return provenance if isinstance(provenance, Mapping) else None


def _prepare_annotation_provenance_transfer(
    adata: Any,
    subset_adata: Any,
    *,
    source_column: str,
    target_column: str,
) -> tuple[str, Any, dict[str, Any] | None]:
    source = _annotation_provenance(subset_adata, source_column)
    if source is None:
        return "not_available", None, None
    transferred = copy.deepcopy(dict(source))
    transferred["output_column"] = target_column
    metadata = adata.uns.get("openbio_singlecell")
    annotations = metadata.get("annotations") if isinstance(metadata, Mapping) else None
    if isinstance(annotations, Mapping) and target_column in annotations:
        if _provenance_token(annotations[target_column]) != _provenance_token(transferred):
            raise ValueError(
                f"Target annotation provenance for {target_column!r} differs from the transferred source provenance."
            )
        return "verified_equal", transferred.get("annotation_status"), None
    return "transferred", transferred.get("annotation_status"), transferred


def _store_annotation_provenance(output: Any, target_column: str, provenance: dict[str, Any] | None) -> None:
    if provenance is None:
        return
    metadata = output.uns["openbio_singlecell"]
    annotations = copy.deepcopy(metadata.get("annotations", {}))
    annotations[target_column] = provenance
    metadata["annotations"] = annotations
    output.uns["openbio_singlecell"] = metadata


def _merge_annotation_values(
    adata: Any,
    subset_adata: Any,
    *,
    source_column: str,
    target_column: str,
    conflict_policy: str,
) -> tuple[Any, dict[str, Any]]:
    science = dependencies.require_scientific_dependencies()
    source = subset_adata.obs[source_column].copy()
    source_dtype = str(source.dtype)
    source_ordered = bool(source.cat.ordered) if isinstance(source.dtype, science.pd.CategoricalDtype) else None
    aligned = source.reindex(adata.obs_names)
    source_present = aligned.notna()
    source_missing = int(source.isna().sum())
    target_exists = target_column in adata.obs
    target = (
        adata.obs[target_column].copy()
        if target_exists
        else science.pd.Series(science.pd.NA, index=adata.obs_names, dtype="object")
    )
    target_before_dtype = str(target.dtype) if target_exists else None
    target_is_categorical = target_exists and isinstance(target.dtype, science.pd.CategoricalDtype)
    target_ordered_before = bool(target.cat.ordered) if target_is_categorical else None

    target_present = target.notna()
    fill_mask = source_present & ~target_present
    identical_mask = science.pd.Series(False, index=adata.obs_names)
    conflict_mask = science.pd.Series(False, index=adata.obs_names)
    for name in adata.obs_names[source_present & target_present]:
        if _equal_annotation_values(aligned.loc[name], target.loc[name]):
            identical_mask.loc[name] = True
        else:
            conflict_mask.loc[name] = True
    conflict_count = int(conflict_mask.sum())
    if conflict_policy == "error" and conflict_count:
        examples = list(adata.obs_names[conflict_mask.to_numpy()][:8])
        raise ValueError(
            f"Observation annotation merge found {conflict_count} conflicting value(s); examples: {examples!r}"
        )

    resolved_object = target.astype("object").copy()
    resolved_object.loc[fill_mask] = aligned.loc[fill_mask].astype("object")
    overwritten = 0
    kept_target = 0
    if conflict_policy == "overwrite":
        resolved_object.loc[conflict_mask] = aligned.loc[conflict_mask].astype("object")
        overwritten = conflict_count
    elif conflict_policy == "keep_target":
        kept_target = conflict_count

    # The worker owns this AnnData instance, so annotation transfer mutates it directly.
    output = adata
    dtype_fallback_to_object = False
    source_is_categorical = isinstance(source.dtype, science.pd.CategoricalDtype)
    if not target_exists:
        output.obs[target_column] = aligned.copy()
    elif not bool(fill_mask.any()) and overwritten == 0:
        output.obs[target_column] = target.copy()
    elif source_is_categorical and target_is_categorical:
        categories = list(target.cat.categories)
        for category in source.cat.categories:
            if not any(_equal_annotation_values(category, existing) for existing in categories):
                categories.append(category)
        ordered = bool(
            target.cat.ordered and source.cat.ordered and target.cat.categories.equals(source.cat.categories)
        )
        output.obs[target_column] = science.pd.Categorical(
            resolved_object,
            categories=categories,
            ordered=ordered,
        )
    else:
        resolved_native = target.copy()
        try:
            resolved_native.loc[fill_mask] = aligned.loc[fill_mask]
            if conflict_policy == "overwrite":
                resolved_native.loc[conflict_mask] = aligned.loc[conflict_mask]
        except (TypeError, ValueError):
            dtype_fallback_to_object = True
            output.obs[target_column] = resolved_object
        else:
            if resolved_native.dtype != target.dtype:
                dtype_fallback_to_object = True
                output.obs[target_column] = resolved_object
            else:
                output.obs[target_column] = resolved_native

    statistics = {
        "source_missing": source_missing,
        "filled": int(fill_mask.sum()),
        "identical": int(identical_mask.sum()),
        "conflicts": conflict_count,
        "overwritten": overwritten,
        "kept_target": kept_target,
        "target_column_existed": target_exists,
        "source_dtype": source_dtype,
        "source_categorical_ordered": source_ordered,
        "target_dtype_before": target_before_dtype,
        "target_dtype_after": str(output.obs[target_column].dtype),
        "dtype_fallback_to_object": dtype_fallback_to_object,
        "target_categorical_ordered_before": target_ordered_before,
        "target_categorical_ordered_after": (
            bool(output.obs[target_column].cat.ordered)
            if isinstance(output.obs[target_column].dtype, science.pd.CategoricalDtype)
            else None
        ),
        "final_non_missing": int(output.obs[target_column].notna().sum()),
    }
    return output, statistics


def _merge_observation_annotations_code(
    source_column: str,
    target_column: str,
    conflict_policy: str,
) -> str:
    return dedent(
        f"""
        import copy
        from collections.abc import Mapping

        import pandas as pd


        def _equal_values(left, right):
            try:
                result = left == right
                return bool(result) if isinstance(result, bool) or type(result).__module__.startswith("numpy") else False
            except (TypeError, ValueError):
                return False


        def _provenance_token(value):
            if isinstance(value, Mapping):
                return {{str(key): _provenance_token(item) for key, item in value.items()}}
            if isinstance(value, (list, tuple)):
                return [_provenance_token(item) for item in value]
            if not isinstance(value, (str, bytes, bytearray)) and hasattr(value, "tolist"):
                try:
                    return _provenance_token(value.tolist())
                except (TypeError, ValueError):
                    pass
            if hasattr(value, "item"):
                try:
                    return _provenance_token(value.item())
                except (TypeError, ValueError):
                    pass
            return value


        def merge_observation_annotations(adata, subset_adata):
            source_column = {source_column!r}.strip()
            target_column = {target_column!r}.strip()
            conflict_policy = {conflict_policy!r}
            if adata.isbacked or subset_adata.isbacked:
                raise ValueError("Observation annotation transfer requires in-memory AnnData inputs; call to_memory() first.")
            if not source_column or not target_column:
                raise ValueError("Source and target observation annotation names cannot be empty.")
            if source_column not in subset_adata.obs:
                raise ValueError(f"Source annotation not found: {{source_column!r}}")
            if not adata.obs_names.is_unique or not subset_adata.obs_names.is_unique:
                raise ValueError("Annotation transfer requires unique obs_names in both inputs.")
            target_metadata = adata.uns.get("openbio_singlecell")
            subset_metadata = subset_adata.uns.get("openbio_singlecell")
            target_source = target_metadata.get("source") if isinstance(target_metadata, Mapping) else None
            subset_source = subset_metadata.get("source") if isinstance(subset_metadata, Mapping) else None
            if (
                isinstance(target_source, Mapping)
                and target_source
                and isinstance(subset_source, Mapping)
                and subset_source
                and _provenance_token(target_source) != _provenance_token(subset_source)
            ):
                raise ValueError(
                    "The inputs carry different OpenBio source provenance; annotation transfer requires a common source."
                )
            source_metadata = subset_adata.uns.get("openbio_singlecell")
            source_annotations = (
                source_metadata.get("annotations") if isinstance(source_metadata, Mapping) else None
            )
            source_provenance = (
                source_annotations.get(source_column) if isinstance(source_annotations, Mapping) else None
            )
            transferred_provenance = None
            if isinstance(source_provenance, Mapping):
                candidate = copy.deepcopy(dict(source_provenance))
                candidate["output_column"] = target_column
                target_metadata = adata.uns.get("openbio_singlecell")
                target_annotations = (
                    target_metadata.get("annotations") if isinstance(target_metadata, Mapping) else None
                )
                if isinstance(target_annotations, Mapping) and target_column in target_annotations:
                    if _provenance_token(target_annotations[target_column]) != _provenance_token(candidate):
                        raise ValueError(
                            f"Target annotation provenance for {{target_column!r}} differs from the "
                            "transferred source provenance."
                        )
                else:
                    transferred_provenance = candidate
            source = subset_adata.obs[source_column].copy()
            aligned = source.reindex(adata.obs_names)
            source_present = aligned.notna()
            target_exists = target_column in adata.obs
            target = (
                adata.obs[target_column].copy()
                if target_exists
                else pd.Series(pd.NA, index=adata.obs_names, dtype="object")
            )
            target_present = target.notna()
            fill_mask = source_present & ~target_present
            conflict_mask = pd.Series(False, index=adata.obs_names)
            for name in adata.obs_names[source_present & target_present]:
                if not _equal_values(aligned.loc[name], target.loc[name]):
                    conflict_mask.loc[name] = True
            conflict_count = int(conflict_mask.sum())
            if conflict_policy == "error" and conflict_count:
                raise ValueError(f"Observation annotation merge found {{conflict_count}} conflicting value(s).")
            resolved_object = target.astype("object").copy()
            resolved_object.loc[fill_mask] = aligned.loc[fill_mask].astype("object")
            if conflict_policy == "overwrite":
                resolved_object.loc[conflict_mask] = aligned.loc[conflict_mask].astype("object")
            output = adata.copy()
            source_is_categorical = isinstance(source.dtype, pd.CategoricalDtype)
            target_is_categorical = target_exists and isinstance(target.dtype, pd.CategoricalDtype)
            if not target_exists:
                output.obs[target_column] = aligned.copy()
            elif not bool(fill_mask.any()) and (
                conflict_policy != "overwrite" or not bool(conflict_mask.any())
            ):
                output.obs[target_column] = target.copy()
            elif source_is_categorical and target_is_categorical:
                categories = list(target.cat.categories)
                for category in source.cat.categories:
                    if not any(_equal_values(category, existing) for existing in categories):
                        categories.append(category)
                ordered = bool(
                    target.cat.ordered
                    and source.cat.ordered
                    and target.cat.categories.equals(source.cat.categories)
                )
                output.obs[target_column] = pd.Categorical(
                    resolved_object,
                    categories=categories,
                    ordered=ordered,
                )
            else:
                resolved_native = target.copy()
                try:
                    resolved_native.loc[fill_mask] = aligned.loc[fill_mask]
                    if conflict_policy == "overwrite":
                        resolved_native.loc[conflict_mask] = aligned.loc[conflict_mask]
                except (TypeError, ValueError):
                    output.obs[target_column] = resolved_object
                else:
                    output.obs[target_column] = (
                        resolved_native if resolved_native.dtype == target.dtype else resolved_object
                    )
            if transferred_provenance is not None:
                metadata = copy.deepcopy(output.uns.get("openbio_singlecell"))
                if metadata is None:
                    metadata = {{
                        "schema_version": 1,
                        "version": "not-installed",
                        "display_name": "AnnData",
                        "source": {{}},
                        "random_seed": 0,
                        "warnings": [],
                        "analysis_history": {{}},
                    }}
                annotations = copy.deepcopy(metadata.get("annotations", {{}}))
                annotations[target_column] = transferred_provenance
                metadata["annotations"] = annotations
                output.uns["openbio_singlecell"] = metadata
            return output
        """
    )


def _header_evidence(headers: Mapping[str, str], *keys: str) -> dict[str, Any]:
    for key in keys:
        value = headers.get(key)
        if value:
            rendered = str(value)
            return {"key": key, "value": rendered[:512], "truncated": len(rendered) > 512}
    return {"key": None, "value": None, "truncated": False}


def _bounded_gtf_headers(headers: Mapping[str, str], *, limit: int = 32) -> dict[str, str]:
    return {str(key)[:128]: str(value)[:512] for key, value in list(headers.items())[:limit]}


def _validate_gtf_identity_inputs(
    adata: Any,
    *,
    id_source: str,
    id_column: str,
    gene_name_column: str,
) -> tuple[str, str, list[str], list[str], dict[str, int]]:
    science = dependencies.require_scientific_dependencies()
    if id_source not in {"var_names", "var_column"}:
        raise ValueError(f"Unsupported gene ID source: {id_source!r}")
    id_column = id_column.strip()
    gene_name_column = gene_name_column.strip()
    if not gene_name_column:
        raise ValueError("Gene-name annotation column cannot be empty.")
    if id_source == "var_column" and not id_column:
        raise ValueError("Gene ID annotation column cannot be empty when id_source='var_column'.")
    if id_source == "var_column" and id_column not in adata.var:
        raise ValueError(f"Stable gene ID annotation not found in var: {id_column!r}")
    if id_source == "var_column" and id_column == gene_name_column:
        raise ValueError("Gene ID and gene-name annotation columns must be different.")
    if gene_name_column in {ORIGINAL_FEATURE_NAMES_COLUMN, GTF_MAPPED_COLUMN, GTF_MATCH_COLUMN}:
        raise ValueError(f"Gene-name annotation column is reserved: {gene_name_column!r}")
    if GTF_MAPPED_COLUMN in adata.var or GTF_MATCH_COLUMN in adata.var:
        raise ValueError("GTF mapping evidence columns already exist; refusing to overwrite them.")
    if id_source == "var_column" and ORIGINAL_FEATURE_NAMES_COLUMN in adata.var:
        raise ValueError(f"Reserved original feature-name annotation already exists: {ORIGINAL_FEATURE_NAMES_COLUMN!r}")

    source_series = (
        science.pd.Series(adata.var_names, index=adata.var_names) if id_source == "var_names" else adata.var[id_column]
    )
    if bool(source_series.isna().any()):
        raise ValueError("Stable gene IDs contain missing values.")
    gene_ids = [str(value).strip() for value in source_series]
    blank_count = sum(not value for value in gene_ids)
    if blank_count:
        raise ValueError(f"Stable gene IDs contain {blank_count} blank value(s).")
    normalized_ids = [_stable_gene_id(value) for value in gene_ids]
    existing_non_string_count = 0
    if gene_name_column in adata.var:
        existing_non_string_count = sum(
            not (isinstance(science.pd.isna(value), (bool, science.np.bool_)) and bool(science.pd.isna(value)))
            and not isinstance(value, str)
            for value in adata.var[gene_name_column]
        )
    return (
        id_column,
        gene_name_column,
        gene_ids,
        normalized_ids,
        {
            "duplicate_stable_id_occurrences": len(gene_ids) - len(set(gene_ids)),
            "version_normalized_id_collisions": len(normalized_ids) - len(set(normalized_ids)),
            "existing_non_string_gene_names": existing_non_string_count,
        },
    )


def _gtf_identity_operation(
    adata: Any,
    mapping: _GTFGeneMapping,
    *,
    id_source: str,
    id_column: str,
    gene_name_column: str,
) -> tuple[Any, dict[str, Any]]:
    science = dependencies.require_scientific_dependencies()
    id_column, gene_name_column, gene_ids, normalized_ids, identity_audit = _validate_gtf_identity_inputs(
        adata,
        id_source=id_source,
        id_column=id_column,
        gene_name_column=gene_name_column,
    )
    mapped_names: list[str | None] = []
    match_types = []
    exact_count = 0
    version_count = 0
    for gene_id, normalized_id in zip(gene_ids, normalized_ids, strict=True):
        if gene_id in mapping.exact:
            mapped_names.append(mapping.exact[gene_id])
            match_types.append("exact")
            exact_count += 1
        elif normalized_id in mapping.normalized:
            mapped_names.append(mapping.normalized[normalized_id])
            match_types.append("version_normalized")
            version_count += 1
        else:
            mapped_names.append(None)
            match_types.append("unmapped")

    existing = (
        adata.var[gene_name_column].copy()
        if gene_name_column in adata.var
        else science.pd.Series(science.pd.NA, index=adata.var_names, dtype="object")
    )
    existing_equal = 0
    existing_filled = 0
    conflicts = []
    final_names = existing.astype("object").tolist()
    for position, mapped_name in enumerate(mapped_names):
        if mapped_name is None:
            continue
        current = existing.iloc[position]
        missing_result = science.pd.isna(current)
        current_is_missing = bool(missing_result) if isinstance(missing_result, (bool, science.np.bool_)) else False
        if current_is_missing:
            final_names[position] = mapped_name
            existing_filled += 1
        elif _equal_annotation_values(current, mapped_name):
            final_names[position] = mapped_name
            existing_equal += 1
        else:
            conflicts.append(f"{gene_ids[position]!r}: existing={str(current)!r}, GTF={mapped_name!r}")

    # The worker owns this AnnData instance; only small annotation columns are replaced.
    output = adata
    if id_source == "var_column":
        output.var[ORIGINAL_FEATURE_NAMES_COLUMN] = output.var_names.astype(str).to_numpy()
        output.var_names = science.pd.Index(gene_ids)
    if all(
        isinstance(value, str)
        or (isinstance(science.pd.isna(value), (bool, science.np.bool_)) and bool(science.pd.isna(value)))
        for value in final_names
    ):
        output.var[gene_name_column] = science.pd.Categorical(final_names)
    else:
        output.var[gene_name_column] = final_names
    output.var[GTF_MAPPED_COLUMN] = science.np.asarray(
        [value is not None for value in mapped_names],
        dtype=bool,
    )
    output.var[GTF_MATCH_COLUMN] = science.pd.Categorical(
        match_types,
        categories=["exact", "version_normalized", "unmapped"],
        ordered=False,
    )
    final_symbols = []
    for value in final_names:
        missing_result = science.pd.isna(value)
        if isinstance(missing_result, (bool, science.np.bool_)) and bool(missing_result):
            continue
        final_symbols.append(value)
    unique_symbols = []
    for value in final_symbols:
        if not any(_equal_annotation_values(value, existing_value) for existing_value in unique_symbols):
            unique_symbols.append(value)
    mapped_count = exact_count + version_count
    statistics = {
        "features": int(adata.n_vars),
        "exact_matches": exact_count,
        "version_normalized_matches": version_count,
        "mapped_features": mapped_count,
        "unmapped_features": int(adata.n_vars) - mapped_count,
        "ambiguous_matches": 0,
        "existing_names_equal": existing_equal,
        "existing_names_filled": existing_filled,
        "existing_name_conflicts": len(conflicts),
        "existing_name_conflict_examples": conflicts[:8],
        "duplicate_gene_symbol_occurrences": len(final_symbols) - len(unique_symbols),
        **identity_audit,
    }
    return output, statistics


def _gtf_identity_code(id_source: str, id_column: str, gene_name_column: str) -> str:
    return dedent(
        f"""
        import gzip
        import re
        import warnings

        import numpy as np
        import pandas as pd


        ATTRIBUTE_PATTERN = re.compile(r'([^\\s;]+)\\s+"([^\"]*)"')
        VERSION_PATTERN = re.compile(r"\\.\\d+(?=_PAR_Y$|$)")


        def _stable_gene_id(value):
            return VERSION_PATTERN.sub("", value)


        def _is_missing(value):
            result = pd.isna(value)
            return bool(result) if isinstance(result, (bool, np.bool_)) else False


        def _equal_values(left, right):
            try:
                result = left == right
                return bool(result) if isinstance(result, (bool, np.bool_)) else False
            except (TypeError, ValueError):
                return False


        def _read_mapping(gtf_path):
            gene_exact_sets = {{}}
            gene_normalized_sets = {{}}
            fallback_exact_sets = {{}}
            fallback_normalized_sets = {{}}
            gene_record_count = 0
            fallback_record_count = 0
            opener = gzip.open if str(gtf_path).lower().endswith(".gz") else open
            with opener(gtf_path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip("\\n").split("\\t")
                    if len(fields) != 9:
                        continue
                    attributes = dict(ATTRIBUTE_PATTERN.findall(fields[8]))
                    gene_id = attributes.get("gene_id", "").strip()
                    gene_name = attributes.get("gene_name", "").strip()
                    if not gene_id or not gene_name:
                        continue
                    if fields[2] == "gene":
                        gene_record_count += 1
                        exact_sets = gene_exact_sets
                        normalized_sets = gene_normalized_sets
                    else:
                        fallback_record_count += 1
                        exact_sets = fallback_exact_sets
                        normalized_sets = fallback_normalized_sets
                    exact_sets.setdefault(gene_id, set()).add(gene_name)
                    normalized_sets.setdefault(_stable_gene_id(gene_id), set()).add(gene_name)
            if gene_record_count:
                exact_sets = gene_exact_sets
                normalized_sets = gene_normalized_sets
            else:
                exact_sets = fallback_exact_sets
                normalized_sets = fallback_normalized_sets
            if not gene_record_count and not fallback_record_count:
                raise ValueError("The GTF has no records with gene_id and gene_name.")
            conflicts = [key for key, names in [*exact_sets.items(), *normalized_sets.items()] if len(names) > 1]
            if conflicts:
                raise ValueError(f"The GTF contains ambiguous gene mappings: {{conflicts[:8]!r}}")
            exact = {{key: next(iter(names)) for key, names in exact_sets.items()}}
            normalized = {{key: next(iter(names)) for key, names in normalized_sets.items()}}
            return exact, normalized


        def annotate_gene_ids_from_gtf(adata, gtf_path):
            id_source = {id_source!r}
            id_column = {id_column!r}.strip()
            gene_name_column = {gene_name_column!r}.strip()
            if adata.isbacked:
                raise ValueError("GTF identity annotation requires an in-memory AnnData; call adata.to_memory() first.")
            if not gene_name_column:
                raise ValueError("Gene-name annotation column cannot be empty.")
            if id_source == "var_column" and not id_column:
                raise ValueError("Gene ID annotation column cannot be empty for var_column input.")
            if id_source == "var_column" and id_column not in adata.var:
                raise ValueError(f"Stable gene ID annotation not found: {{id_column!r}}")
            if id_source == "var_column" and id_column == gene_name_column:
                raise ValueError("Gene ID and gene-name columns must differ.")
            if gene_name_column in {{"feature_names_before_gtf", "gtf_mapped", "gtf_match_type"}}:
                raise ValueError(f"Gene-name annotation column is reserved: {{gene_name_column!r}}")
            if "gtf_mapped" in adata.var or "gtf_match_type" in adata.var:
                raise ValueError("GTF mapping evidence columns already exist.")
            if id_source == "var_column" and "feature_names_before_gtf" in adata.var:
                raise ValueError("Reserved original feature-name annotation already exists.")
            source = adata.var_names if id_source == "var_names" else adata.var[id_column]
            if pd.isna(source).any():
                raise ValueError("Stable gene IDs contain missing values.")
            gene_ids = [str(value).strip() for value in source]
            if any(not value for value in gene_ids):
                raise ValueError("Stable gene IDs contain blank values.")
            normalized_ids = [_stable_gene_id(value) for value in gene_ids]
            exact, normalized = _read_mapping(gtf_path)
            mapped_names = []
            match_types = []
            for gene_id, normalized_id in zip(gene_ids, normalized_ids, strict=True):
                if gene_id in exact:
                    mapped_names.append(exact[gene_id])
                    match_types.append("exact")
                elif normalized_id in normalized:
                    mapped_names.append(normalized[normalized_id])
                    match_types.append("version_normalized")
                else:
                    mapped_names.append(None)
                    match_types.append("unmapped")
            existing = (
                adata.var[gene_name_column].copy()
                if gene_name_column in adata.var
                else pd.Series(pd.NA, index=adata.var_names, dtype="object")
            )
            final_names = existing.astype("object").tolist()
            for position, mapped_name in enumerate(mapped_names):
                if mapped_name is None:
                    continue
                current = existing.iloc[position]
                if _is_missing(current):
                    final_names[position] = mapped_name
                elif _equal_values(current, mapped_name):
                    final_names[position] = mapped_name
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
                output = adata.copy()
            if id_source == "var_column":
                output.var["feature_names_before_gtf"] = output.var_names.astype(str).to_numpy()
                output.var_names = pd.Index(gene_ids)
            output.var[gene_name_column] = final_names
            output.var["gtf_mapped"] = np.asarray([value is not None for value in mapped_names], dtype=bool)
            output.var["gtf_match_type"] = pd.Categorical(
                match_types,
                categories=["exact", "version_normalized", "unmapped"],
                ordered=False,
            )
            return output
        """
    )


def _result_records(context: OperationContext, adata: Any, report: Any, code: str) -> list[dict[str, JSONValue]]:
    return analysis_outputs(report, code, write_anndata_output(context, adata))


@register_operation("openbio.node.snapshotexpression")
def snapshot_expression(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    require_input_names(inputs, {"adata"}, operation="Snapshot Expression")
    require_parameters(parameters, {"source", "overwrite_existing"}, operation="Snapshot Expression")
    overwrite_existing = _required_bool(parameters["overwrite_existing"], name="overwrite_existing")
    source = parameters["source"]
    if source is not None and not isinstance(source, dict):
        raise ProtocolError("source must be a DynamicCombo object or null.")
    adata = read_anndata_input(inputs)
    science = dependencies.require_scientific_dependencies()
    expression = SNAPSHOT_EXPRESSION_SOURCE.resolve(adata, source)
    matrix = expression.matrix(adata)
    existing_destinations = []
    if adata.raw is not None:
        existing_destinations.append("raw")
    if CANONICAL_COUNTS_LAYER in adata.layers:
        existing_destinations.append("layers['counts']")
    if existing_destinations and not overwrite_existing:
        raise ValueError(f"Raw snapshot destinations already exist: {', '.join(existing_destinations)}")

    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    # These two matrix copies are scientifically required: layer and Raw must be
    # independent snapshots while the worker-private current X stays unchanged.
    adata.layers[CANONICAL_COUNTS_LAYER] = matrix.copy()
    current_x = adata.X
    try:
        adata.X = matrix.copy()
        adata.raw = adata
    finally:
        adata.X = current_x
    parameters_report = {
        **expression.parameters(),
        "overwrite_existing": overwrite_existing,
    }
    warnings_list = []
    if cells == 0 or genes == 0:
        warnings_list.append(
            "The declared source has an empty observation or variable axis; an empty snapshot was created."
        )
    if not adata.obs_names.is_unique:
        warnings_list.append(
            "Observation names are not unique; the snapshot preserves them, and later name-based alignment may "
            "be ambiguous."
        )
    if not adata.var_names.is_unique:
        warnings_list.append(
            "Variable names are not unique; the snapshot preserves them, and later name-based feature selection "
            "may be ambiguous."
        )
    finish_adata(
        adata,
        "snapshot_expression",
        parameters_report,
        cells,
        genes,
        started_at,
        warnings=warnings_list,
    )
    nonzero = (
        int(matrix.count_nonzero())
        if science.sparse.issparse(matrix)
        else int(science.np.count_nonzero(science.np.asarray(matrix)))
    )
    code = _snapshot_expression_code(expression, overwrite_existing)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellSnapshotExpression",
        title="Expression and Raw snapshot",
        operation="snapshot_expression",
        methods=(
            f"Copied the user-selected expression source {_display_expression_source(expression)!r} into the "
            "repository's conventional AnnData counts layer and initialized raw from the same matrix while "
            "preserving current X."
        ),
        results=(
            f"Retained a user-selected Raw snapshot for {cells:,} cells and {genes:,} features "
            f"({nonzero:,} non-zero entries)."
        ),
        key_results={
            "cells": cells,
            "features": genes,
            "source": _display_expression_source(expression),
            "matrix": _matrix_description(matrix),
            "nonzero_entries": nonzero,
            "canonical_layer": CANONICAL_COUNTS_LAYER,
            "raw_created": True,
            "x_preserved": True,
            "overwritten_destinations": existing_destinations if overwrite_existing else [],
        },
        parameters=parameters_report,
        references=[ANNDATA_REFERENCE],
        software_packages=["anndata", "numpy", "pandas", "scipy"],
        warnings=warnings_list,
        limitations=["The snapshot preserves the selected values without interpreting their preprocessing state."],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=code,
    )
    return _result_records(context, adata, report, code)


@register_operation("openbio.node.rawsnapshottoanndata")
def raw_snapshot_to_anndata(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    require_input_names(inputs, {"adata"}, operation="Raw Snapshot to AnnData")
    require_parameters(parameters, set(), operation="Raw Snapshot to AnnData")
    adata = read_anndata_input(inputs)
    if adata.raw is None:
        raise ValueError("Raw Snapshot to AnnData requires adata.raw; create a Raw snapshot first.")
    if adata.n_obs == 0:
        raise ValueError("Raw Snapshot to AnnData requires at least one observation; received zero cells.")
    if adata.raw.n_vars == 0:
        raise ValueError("Raw Snapshot to AnnData requires at least one Raw feature; adata.raw has zero variables.")

    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    raw_genes = int(adata.raw.n_vars)
    metadata = copy.deepcopy(adata.uns.get("openbio_singlecell"))
    output = adata.raw.to_adata()
    if output.layers or output.raw is not None:
        raise RuntimeError("AnnData Raw materialization unexpectedly retained layers or nested Raw storage.")
    if not output.obs.equals(adata.obs) or not output.var.equals(adata.raw.var):
        raise RuntimeError("AnnData Raw materialization changed current observations or Raw feature annotations.")
    output.obsm.clear()
    output.obsp.clear()
    output.varm.clear()
    output.varp.clear()
    output.uns.clear()
    if metadata is not None:
        output.uns["openbio_singlecell"] = metadata
    finish_adata(
        output,
        "raw_snapshot_to_anndata",
        {},
        cells,
        genes,
        started_at,
    )
    code = _raw_snapshot_to_anndata_code()
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellRawSnapshotToAnnData",
        title="Raw snapshot materialization",
        operation="raw_snapshot_to_anndata",
        methods=(
            "Materialized the current observation subset with adata.raw.to_adata(), retained Raw X and var "
            "with current obs, and removed inherited analysis representations and non-OpenBio uns entries."
        ),
        results=(
            f"Created a clean AnnData with {cells:,} cells and {raw_genes:,} Raw features; the input current "
            f"feature axis contained {genes:,} features."
        ),
        key_results={
            "cells": cells,
            "input_features": genes,
            "raw_features": raw_genes,
            "output_features": int(output.n_vars),
            "matrix": _matrix_description(output.X),
            "current_obs_preserved": True,
            "raw_var_preserved": True,
            "derived_state_cleared": ["layers", "obsm", "obsp", "varm", "varp", "raw"],
            "uns_keys": list(output.uns.keys()),
        },
        parameters={},
        references=[ANNDATA_REFERENCE],
        software_packages=["anndata", "numpy", "pandas", "scipy"],
        warnings=[],
        limitations=[
            "Raw is an independent expression snapshot; materialization does not infer whether its values are "
            "counts, normalized, or transformed."
        ],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=code,
    )
    return _result_records(context, output, report, code)


@register_operation("openbio.node.subsetobservations")
def subset_observations(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    require_input_names(inputs, {"adata"}, operation="Subset Observations")
    require_parameters(
        parameters,
        {"column", "values", "invert", "missing_policy"},
        operation="Subset Observations",
    )
    column = _required_string(parameters["column"], name="column").strip()
    values = _required_string(parameters["values"], name="values")
    invert = _required_bool(parameters["invert"], name="invert")
    missing_policy = _required_string(parameters["missing_policy"], name="missing_policy")
    adata = read_anndata_input(inputs)
    if not column:
        raise ValueError("Observation annotation name cannot be empty.")
    if column not in adata.obs:
        raise ValueError(f"Observation annotation not found in obs: {column!r}")
    selected_values = _comma_separated_values(values)
    if missing_policy not in {"exclude", "include", "error"}:
        raise ValueError(f"Unsupported missing_policy: {missing_policy!r}")

    series = adata.obs[column]
    _string_label_collision(series, label="Observation annotation")
    not_missing = series.notna()
    missing_count = int((~not_missing).sum())
    if missing_policy == "error" and missing_count:
        raise ValueError(f"Observation annotation contains {missing_count} missing value(s).")
    rendered = series.astype("string")
    matched = not_missing & rendered.isin(selected_values)
    mask = matched if not invert else not_missing & ~matched
    if missing_policy == "include":
        mask = mask | ~not_missing

    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    # AnnData slicing is a view; materialization is required for an independent output artifact.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*names are not unique.*", category=UserWarning)
        output = adata[mask.to_numpy(), :].copy()
    available = {str(value) for value in rendered[not_missing].unique()}
    unmatched = [value for value in selected_values if value not in available]
    input_counts = {value: int((rendered == value).fillna(False).sum()) for value in selected_values}
    output_rendered = output.obs[column].astype("string")
    output_counts = {value: int((output_rendered == value).fillna(False).sum()) for value in selected_values}
    parameters_report = {
        "column": column,
        "values": list(selected_values),
        "invert": invert,
        "missing_policy": missing_policy,
    }
    warnings_list = []
    if unmatched:
        warnings_list.append(f"Requested observation values were not present: {unmatched!r}")
    if not bool(mask.any()):
        warnings_list.append(
            "The explicit observation selection retained zero cells; the empty AnnData result was returned."
        )
    if not adata.obs_names.is_unique:
        warnings_list.append(
            "Observation names are not unique; position-based subsetting is valid, but later cross-object "
            "annotation transfer will be ambiguous."
        )
    finish_adata(
        output,
        "subset_observations",
        parameters_report,
        cells,
        genes,
        started_at,
        warnings=warnings_list,
    )
    retained = int(output.n_obs)
    removed = cells - retained
    code = _subset_observations_code(column, selected_values, invert, missing_policy)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellSubsetObservations",
        title=f"Observation subset by {column}",
        operation="subset_observations",
        methods=(
            f"Selected observations by exact string membership in obs[{column!r}], with invert={invert!r} "
            f"and missing_policy={missing_policy!r}; AnnData slicing preserved all observation-aligned slots."
        ),
        results=(
            f"Retained {retained:,} of {cells:,} cells ({retained / cells:.1%}) and removed {removed:,} "
            f"using the user-declared annotation values {list(selected_values)!r}."
            if cells
            else "The empty input contained no observations."
        ),
        key_results={
            "input_cells": cells,
            "retained_cells": retained,
            "removed_cells": removed,
            "retained_fraction": float(retained / cells) if cells else None,
            "features_preserved": genes,
            "annotation_column": column,
            "requested_values": list(selected_values),
            "input_counts_by_requested_value": input_counts,
            "retained_counts_by_requested_value": output_counts,
            "missing_annotations": missing_count,
            "unmatched_requested_values": unmatched,
            "input_value_counts": _bounded_value_counts(series),
            "observation_order_preserved": True,
            "empty_output": retained == 0,
        },
        parameters=parameters_report,
        references=[ANNDATA_REFERENCE, PANDAS_REFERENCE],
        software_packages=["anndata", "numpy", "pandas", "scipy"],
        warnings=warnings_list,
        limitations=[
            "The comma-separated widget cannot represent a literal annotation value containing a comma.",
            "Selection by an annotation does not validate that provisional or curated biological labels are correct.",
        ],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=code,
    )
    return _result_records(context, output, report, code)


@register_operation("openbio.node.mergeobservationannotations")
def merge_observation_annotations(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    require_input_names(
        inputs,
        {"adata", "subset_adata"},
        operation="Merge Observation Annotations",
    )
    require_parameters(
        parameters,
        {"source_column", "target_column", "conflict_policy"},
        operation="Merge Observation Annotations",
    )
    source_column = _required_string(parameters["source_column"], name="source_column").strip()
    target_column = _required_string(parameters["target_column"], name="target_column").strip()
    conflict_policy = _required_string(parameters["conflict_policy"], name="conflict_policy")
    adata = read_anndata_input(inputs)
    subset_adata = read_anndata_input(inputs, "subset_adata")
    if not source_column or not target_column:
        raise ValueError("Source and target observation annotation names cannot be empty.")
    if source_column not in subset_adata.obs:
        raise ValueError(f"Source observation annotation not found in subset_adata.obs: {source_column!r}")
    if not adata.obs_names.is_unique or not subset_adata.obs_names.is_unique:
        raise ValueError("Merging observation annotations requires unique obs_names in both AnnData inputs.")
    if conflict_policy not in {"error", "keep_target", "overwrite"}:
        raise ValueError(f"Unsupported conflict_policy: {conflict_policy!r}")

    annotation_provenance_status, transferred_annotation_status, transferred_provenance = (
        _prepare_annotation_provenance_transfer(
            adata,
            subset_adata,
            source_column=source_column,
            target_column=target_column,
        )
    )
    source_only = subset_adata.obs_names.difference(adata.obs_names)
    shared_names = subset_adata.obs_names.intersection(adata.obs_names, sort=False)
    provenance_status = _verify_shared_provenance(adata, subset_adata)
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, statistics = _merge_annotation_values(
        adata,
        subset_adata,
        source_column=source_column,
        target_column=target_column,
        conflict_policy=conflict_policy,
    )
    warnings_list = []
    limitations = [
        "Exact shared observation identifiers are necessary but cannot by themselves prove that the two "
        "objects derive from the same biological source.",
        "Transferred annotations remain provisional or curated according to the declared source column; this "
        "operation does not validate their biological correctness.",
    ]
    if provenance_status == "not_available":
        warnings_list.append(
            "Compatible OpenBio source provenance was not available on both inputs; confirm that observation "
            "identifiers refer to the same cells."
        )
    if len(source_only):
        warnings_list.append(
            f"Ignored {len(source_only)} source observation(s) absent from the target; examples: "
            f"{list(source_only[:8])!r}."
        )
    if len(shared_names) == 0:
        warnings_list.append(
            "The inputs have no shared observation identifiers; the target annotation is unchanged or empty."
        )
    if not bool(subset_adata.obs[source_column].notna().any()):
        warnings_list.append("The source annotation contains no non-missing values; no target values were transferred.")
    if statistics["dtype_fallback_to_object"]:
        warnings_list.append(
            "Transferred values were incompatible with the target annotation dtype; the merged annotation "
            "was represented with pandas object dtype."
        )
    if (
        statistics["target_categorical_ordered_before"] is True
        and statistics["target_categorical_ordered_after"] is False
    ):
        warnings_list.append(
            "The target categorical order was removed because transferred values required a union of "
            "incompatible ordered categories."
        )
    parameters_report = {
        "source_column": source_column,
        "target_column": target_column,
        "conflict_policy": conflict_policy,
    }
    finish_adata(
        output,
        "merge_observation_annotations",
        parameters_report,
        cells,
        genes,
        started_at,
        warnings=warnings_list,
    )
    _store_annotation_provenance(output, target_column, transferred_provenance)
    code = _merge_observation_annotations_code(source_column, target_column, conflict_policy)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMergeObservationAnnotations",
        title=f"Transferred {source_column} to {target_column}",
        operation="merge_observation_annotations",
        methods=(
            f"Transferred obs[{source_column!r}] by exact unique obs_names into obs[{target_column!r}], "
            f"using conflict_policy={conflict_policy!r}; the target observation and variable axes were unchanged."
        ),
        results=(
            f"Among {len(shared_names):,} shared source cells, filled {statistics['filled']:,} missing target "
            f"annotations, observed {statistics['identical']:,} identical values and "
            f"{statistics['conflicts']:,} conflicts; {statistics['overwritten']:,} conflicts were overwritten."
        ),
        key_results={
            "target_cells": cells,
            "source_cells": int(subset_adata.n_obs),
            "shared_cells": int(len(shared_names)),
            "source_only_cells": int(len(source_only)),
            "target_only_cells": cells - int(len(shared_names)),
            "source_missing_values": statistics["source_missing"],
            **statistics,
            "source_column": source_column,
            "target_column": target_column,
            "conflict_policy": conflict_policy,
            "provenance_status": provenance_status,
            "annotation_provenance_status": annotation_provenance_status,
            "transferred_annotation_status": transferred_annotation_status,
            "final_value_counts": _bounded_value_counts(output.obs[target_column]),
            "expression_and_axes_preserved": True,
        },
        parameters=parameters_report,
        references=[ANNDATA_REFERENCE, PANDAS_REFERENCE],
        software_packages=["anndata", "numpy", "pandas", "scipy"],
        warnings=warnings_list,
        limitations=limitations,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=code,
    )
    return _result_records(context, output, report, code)


@register_operation("openbio.node.mapgeneidsfromgtf")
def map_gene_ids_from_gtf(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    require_input_names(inputs, {"adata", "gtf_path"}, operation="Map Gene IDs from GTF")
    require_parameters(
        parameters,
        {"id_source", "id_column", "gene_name_column"},
        operation="Map Gene IDs from GTF",
    )
    id_source = _required_string(parameters["id_source"], name="id_source")
    id_column = _required_string(parameters["id_column"], name="id_column").strip()
    gene_name_column = _required_string(
        parameters["gene_name_column"],
        name="gene_name_column",
    ).strip()
    gtf_path, provenance = require_file_input(inputs, "gtf_path")
    adata = read_anndata_input(inputs)
    _validate_gtf_identity_inputs(
        adata,
        id_source=id_source,
        id_column=id_column,
        gene_name_column=gene_name_column,
    )

    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    mapping = _read_gtf_gene_mapping(str(gtf_path))
    output, statistics = _gtf_identity_operation(
        adata,
        mapping,
        id_source=id_source,
        id_column=id_column,
        gene_name_column=gene_name_column,
    )
    provider_evidence = _header_evidence(
        mapping.headers,
        "provider",
        "annotation_source",
        "genebuild_provider",
    )
    release_evidence = _header_evidence(
        mapping.headers,
        "annotation_release",
        "release",
        "gencode_version",
        "ensembl_release",
    )
    genome_build_evidence = _header_evidence(
        mapping.headers,
        "genome_build",
        "genome_build_accession",
        "assembly",
        "genome_version",
    )
    provider = provider_evidence["value"]
    release = release_evidence["value"]
    genome_build = genome_build_evidence["value"]
    warnings_list = []
    if statistics["unmapped_features"]:
        warnings_list.append(
            f"{statistics['unmapped_features']} stable gene ID(s) were not mapped and remain in the object."
        )
    if statistics["version_normalized_matches"]:
        warnings_list.append(
            f"{statistics['version_normalized_matches']} match(es) required terminal version normalization; "
            "the expression matrix and GTF may use different annotation releases."
        )
    if mapping.used_fallback_records:
        warnings_list.append(
            "The GTF contained no usable gene feature rows; mapping fell back to non-gene records carrying "
            "gene_id and gene_name."
        )
    if mapping.malformed_lines or mapping.missing_attribute_lines:
        warnings_list.append(
            f"Skipped {mapping.malformed_lines} malformed line(s) and "
            f"{mapping.missing_attribute_lines} line(s) missing gene_id/gene_name."
        )
    if cells == 0 or genes == 0:
        warnings_list.append(
            "The input has an empty observation or variable axis; GTF annotations were applied to the "
            "available feature rows."
        )
    if adata.raw is not None:
        warnings_list.append(
            "The existing Raw snapshot was preserved with its independent feature namespace; current and Raw "
            "feature names were not required to match."
        )
    if statistics["mapped_features"] == 0:
        warnings_list.append(
            "None of the declared stable gene IDs matched the selected GTF; all features were retained and "
            "marked unmapped."
        )
    if statistics["duplicate_stable_id_occurrences"]:
        warnings_list.append(
            f"The declared stable IDs contain {statistics['duplicate_stable_id_occurrences']} duplicate "
            "occurrence(s); annotations were applied positionally."
        )
    if statistics["version_normalized_id_collisions"]:
        warnings_list.append(
            f"Terminal version removal yields {statistics['version_normalized_id_collisions']} duplicate ID "
            "occurrence(s); exact matches were preferred and all mappings are disclosed."
        )
    if statistics["existing_name_conflicts"]:
        warnings_list.append(
            f"Preserved {statistics['existing_name_conflicts']} existing gene-name annotation(s) that "
            "disagreed with the selected GTF."
        )
    if statistics["existing_non_string_gene_names"]:
        warnings_list.append(
            f"Preserved {statistics['existing_non_string_gene_names']} non-string existing gene-name value(s); "
            "GTF strings were filled only where values were missing."
        )
    limitations = []
    if not provider:
        limitations.append("The GTF header did not identify an annotation provider.")
    if not release:
        limitations.append("The GTF header did not identify an annotation release.")
    if not genome_build:
        limitations.append("The GTF header did not identify a genome build/assembly.")
    limitations.append(
        "Gene symbols are annotations that may repeat or change between releases; stable IDs remain the feature key."
    )
    if adata.raw is not None:
        limitations.append("Raw is an independent snapshot and retained its original feature annotations and index.")
    portable_path = provenance.get("path")
    if not isinstance(portable_path, str) or not portable_path:
        raise ProtocolError("GTF file provenance must contain a non-empty portable path.")
    parameters_report = {
        "gtf_path": portable_path,
        "id_source": id_source,
        "id_column": id_column,
        "gene_name_column": gene_name_column,
    }
    finish_adata(
        output,
        "map_gene_ids_from_gtf",
        parameters_report,
        cells,
        genes,
        started_at,
        warnings=warnings_list,
    )
    code = _gtf_identity_code(id_source, id_column, gene_name_column)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMapGeneIdsFromGTF",
        title="GTF gene identity annotation",
        operation="map_gene_ids_from_gtf",
        methods=(
            "Matched declared stable gene IDs to GTF gene records by exact identifier first, then by terminal "
            "numeric-version removal that preserves _PAR_Y; annotated gene symbols without filtering features."
        ),
        results=(
            f"Annotated {statistics['mapped_features']:,} of {genes:,} stable gene IDs "
            f"({statistics['exact_matches']:,} exact and "
            f"{statistics['version_normalized_matches']:,} version-normalized); "
            f"{statistics['unmapped_features']:,} features remained unmapped and were retained."
        ),
        key_results={
            **statistics,
            "mapped_fraction": float(statistics["mapped_features"] / genes) if genes else None,
            "id_source": id_source,
            "id_column": id_column if id_source == "var_column" else None,
            "gene_name_column": gene_name_column,
            "stable_ids_are_var_names": True,
            "features_filtered": 0,
            "raw_present": adata.raw is not None,
            "raw_feature_namespace_preserved": True,
            "gtf": {
                **provenance,
                "provider": provider,
                "provider_header": provider_evidence,
                "release": release,
                "release_header": release_evidence,
                "genome_build": genome_build,
                "genome_build_header": genome_build_evidence,
                "headers": _bounded_gtf_headers(mapping.headers),
                "header_count": len(mapping.headers),
                "headers_truncated": len(mapping.headers) > 32,
                "used_fallback_records": mapping.used_fallback_records,
                "total_data_lines": mapping.total_data_lines,
                "parsed_records": mapping.parsed_records,
                "unique_gene_ids": mapping.unique_gene_ids,
                "malformed_lines": mapping.malformed_lines,
                "missing_attribute_lines": mapping.missing_attribute_lines,
                "repeated_records": mapping.repeated_records,
            },
        },
        parameters=parameters_report,
        references=[GENCODE_REFERENCE, ENSEMBL_REFERENCE, ANNDATA_REFERENCE],
        software_packages=["anndata", "numpy", "pandas"],
        warnings=warnings_list,
        limitations=limitations,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        code=code,
    )
    return _result_records(context, output, report, code)


__all__ = [
    "map_gene_ids_from_gtf",
    "merge_observation_annotations",
    "raw_snapshot_to_anndata",
    "snapshot_expression",
    "subset_observations",
]
