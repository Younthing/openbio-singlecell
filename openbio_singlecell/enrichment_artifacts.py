from __future__ import annotations

import inspect
import math
import numbers
import textwrap
from collections.abc import Mapping
from typing import Any, Literal

from .contracts import TableResult
from .marker_evidence import (
    MARKER_COLUMNS,
    MARKER_UNIVERSE_COLUMNS,
    _canonical_sha256,
    _frame_content_fingerprint,
    validate_marker_artifact_pair,
)

ENRICHMENT_EVIDENCE_SCHEMA_VERSION = 1
GENERIC_RANKED_ARTIFACT_ROLE = "complete_ranked_evidence"
GENERIC_SELECTED_ARTIFACT_ROLE = "selected_gene_evidence"
GENERIC_UNIVERSE_ARTIFACT_ROLE = "tested_gene_universe"
EVIDENCE_SCOPES = ("cluster_marker_evidence", "condition_contrast")
_SHA256_LENGTH = 64


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == _SHA256_LENGTH
        and all(character in "0123456789abcdef" for character in value)
    )


def _declaration_value(value: Any) -> Any:
    """Preserve JSON-safe caller declarations without treating them as verified facts."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, numbers.Real):
        normalized = float(value)
        return normalized if math.isfinite(normalized) else str(normalized)
    if hasattr(value, "item"):
        try:
            return _declaration_value(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _declaration_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
        sequence = value.tolist() if hasattr(value, "tolist") else value
        return [_declaration_value(item) for item in sequence]
    return f"<{type(value).__module__}.{type(value).__qualname__}>"


def _json_cell(value: Any, *, operation: str | None = None) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, numbers.Real):
        normalized = float(value)
        if not math.isfinite(normalized):
            if operation is not None:
                raise ValueError(f"{operation} generic evidence contains non-finite numbers.")
            raise ValueError("Enrichment evidence cannot contain non-finite numeric values.")
        return {"float_hex": normalized.hex()}
    if hasattr(value, "item"):
        if operation is not None:
            return _json_cell(value.item(), operation=operation)
        try:
            return _json_cell(value.item())
        except (TypeError, ValueError):
            pass
    if operation is not None:
        raise TypeError(f"{operation} generic evidence contains unsupported values.")
    raise TypeError(
        "Enrichment evidence content fingerprints support only strings, booleans, integers, finite real "
        f"numbers, and null values; received {type(value).__name__}."
    )


def generic_enrichment_content_fingerprint(
    frame: Any,
    *,
    artifact_role: str,
    operation: str | None = None,
) -> str:
    """Return a canonical current-content fingerprint for a generic evidence table."""
    from pandas import DataFrame

    if not isinstance(frame, DataFrame):
        raise TypeError("Generic enrichment evidence must be a pandas DataFrame.")
    rows = [
        [_json_cell(value, operation=operation) for value in row]
        for row in frame.itertuples(index=False, name=None)
    ]
    return _canonical_sha256(
        {
            "schema": f"openbio-singlecell/{artifact_role}/v1",
            "columns": list(frame.columns),
            "rows": rows,
        }
    )


def enrichment_universe_identity_fingerprint(ordered_genes: list[str]) -> str:
    return _canonical_sha256(
        {
            "schema": "openbio-singlecell/tested-gene-universe-identity/v2",
            "ordered_genes": ordered_genes,
        }
    )


def enrichment_universe_content_fingerprint(frame: Any) -> str:
    """Return the marker-compatible current-content fingerprint of a universe table."""
    from pandas import DataFrame

    if not isinstance(frame, DataFrame):
        raise TypeError("Tested-gene universe must be a pandas DataFrame.")
    if list(frame.columns) != MARKER_UNIVERSE_COLUMNS:
        raise ValueError(
            f"Tested-gene universe must use the exact canonical columns in order: {MARKER_UNIVERSE_COLUMNS}."
        )
    return _frame_content_fingerprint(
        frame,
        artifact="tested-gene-universe",
        numeric_columns={"universe_rank"},
    )


def generic_ranking_fingerprint(
    frame: Any,
    *,
    comparison_column: str,
    gene_column: str,
    score_column: str,
) -> str:
    """Fingerprint every comparison's explicit row/tie order and finite score."""
    from pandas import DataFrame

    if not isinstance(frame, DataFrame):
        raise TypeError("Complete ranked evidence must be a pandas DataFrame.")
    required = [comparison_column, gene_column, score_column]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(f"Complete ranked evidence is missing columns: {missing}.")
    rows = []
    for ordinal, (comparison, gene, score) in enumerate(frame[required].itertuples(index=False, name=None), start=1):
        if not isinstance(comparison, str) or not comparison.strip() or comparison != comparison.strip():
            raise ValueError("Complete ranked evidence comparisons must be nonblank strings without whitespace.")
        if not isinstance(gene, str) or not gene.strip() or gene != gene.strip():
            raise ValueError("Complete ranked evidence genes must be nonblank strings without whitespace.")
        if isinstance(score, bool) or not isinstance(score, numbers.Real) or not math.isfinite(float(score)):
            raise ValueError("Complete ranked evidence scores must be finite real numbers.")
        rows.append([comparison, gene, float(score).hex(), ordinal])
    return _canonical_sha256(
        {
            "schema": "openbio-singlecell/complete-ranked-evidence-order/v1",
            "comparison_column": comparison_column,
            "gene_column": gene_column,
            "score_column": score_column,
            "rows": rows,
        }
    )


def _source_parameters(result: TableResult, *, artifact_name: str) -> Mapping[str, Any]:
    if not isinstance(result.source, Mapping):
        raise ValueError(f"{artifact_name} provenance source must be a mapping.")
    source_parameters = result.source.get("parameters")
    if not isinstance(source_parameters, Mapping):
        raise ValueError(f"{artifact_name} provenance parameters must be a mapping.")
    if dict(source_parameters) != result.parameters:
        raise ValueError(f"{artifact_name} parameters disagree with its immutable analysis source.")
    return source_parameters


def _validated_universe(
    frame: Any,
    *,
    np: Any,
    pd: Any,
    operation: str | None = None,
) -> list[str]:
    if not isinstance(frame, pd.DataFrame):
        if operation is not None:
            raise TypeError(f"{operation} universe must be a pandas DataFrame or table artifact.")
        raise TypeError("Tested-gene universe must be a pandas DataFrame.")
    if list(frame.columns) != MARKER_UNIVERSE_COLUMNS:
        if operation is not None:
            raise ValueError(f"{operation} universe must use canonical gene/universe_rank columns.")
        raise ValueError(
            f"Tested-gene universe must use the exact canonical columns in order: {MARKER_UNIVERSE_COLUMNS}."
        )
    genes = frame["gene"].tolist()
    if operation is not None:
        if not genes or len(genes) != len(set(genes)):
            raise ValueError(f"{operation} universe genes must be nonempty and unique.")
        for gene in genes:
            if not isinstance(gene, str) or not gene.strip() or gene != gene.strip():
                raise ValueError(f"{operation} universe genes must be nonblank strings without whitespace.")
    else:
        for gene in genes:
            if not isinstance(gene, str) or not gene.strip() or gene != gene.strip():
                raise ValueError("Tested-gene universe identifiers must be nonblank strings without whitespace.")
        if not genes or len(genes) != len(set(genes)):
            raise ValueError("Tested-gene universe identifiers must be nonempty and unique.")
    ranks = frame["universe_rank"]
    if pd.api.types.is_bool_dtype(ranks.dtype) or not pd.api.types.is_numeric_dtype(ranks.dtype):
        if operation is not None:
            raise TypeError(f"{operation} universe ranks must be numeric and non-boolean.")
        raise TypeError("Tested-gene universe ranks must be numeric and non-boolean.")
    if not bool(np.array_equal(ranks.to_numpy(dtype=float), np.arange(1, len(genes) + 1, dtype=float))):
        if operation is not None:
            raise ValueError(f"{operation} universe ranks must be consecutive and one-based.")
        raise ValueError("Tested-gene universe ranks must be consecutive, one-based, and order-stable.")
    return genes


def _selected_comparison_frame(
    frame: Any,
    *,
    comparison_column: str,
    selector: str,
    pd: Any,
    operation: str | None = None,
) -> tuple[Any, str, list[str]]:
    if comparison_column not in frame:
        if operation is not None:
            raise ValueError(f"{operation} table is missing comparison/gene columns.")
        raise ValueError(f"Enrichment evidence comparison column not found: {comparison_column!r}.")
    labels = frame[comparison_column].tolist()
    for label in labels:
        if not isinstance(label, str) or not label.strip() or label != label.strip():
            if operation is not None:
                raise ValueError(f"{operation} comparison identifiers are invalid.")
            raise ValueError("Enrichment evidence comparisons must be nonblank strings without whitespace.")
    ordered_labels = list(dict.fromkeys(labels))
    if operation is not None:
        selected = frame.loc[frame[comparison_column] == selector].copy().reset_index(drop=True)
        if selected.empty:
            raise ValueError(f"{operation} comparison {selector!r} is absent or empty.")
        return selected, selector, ordered_labels
    selector = selector.strip() if isinstance(selector, str) else selector
    if not isinstance(selector, str):
        raise TypeError("Enrichment comparison selector must be a string.")
    if not selector:
        if len(ordered_labels) != 1:
            raise ValueError(
                "Enrichment comparison selector is required when the evidence artifact contains multiple comparisons."
            )
        selector = ordered_labels[0]
    if selector not in ordered_labels:
        raise ValueError(f"Enrichment evidence contains no comparison {selector!r}.")
    selected = frame.loc[frame[comparison_column] == selector].copy().reset_index(drop=True)
    if not isinstance(selected, pd.DataFrame) or selected.empty:
        raise ValueError(f"Enrichment evidence comparison {selector!r} is empty.")
    return selected, selector, ordered_labels


def _declared_scope(
    parameters: Mapping[str, Any],
    *,
    artifact_name: str,
) -> tuple[dict[str, Any], list[str]]:
    """Normalize inference declarations for reporting; these fields are not an authority boundary."""
    warnings: list[str] = []
    scope = parameters.get("evidence_scope")
    if not isinstance(scope, str) or not scope.strip() or scope != scope.strip():
        warnings.append(
            f"Generic enrichment {artifact_name} does not declare a canonical evidence_scope; scope is reported "
            "as 'unspecified' and formal Condition interpretation is not validated."
        )
        scope = "unspecified"
    elif scope not in EVIDENCE_SCOPES:
        warnings.append(
            f"Generic enrichment {artifact_name} declares unfamiliar evidence_scope {scope!r}; the value is "
            "preserved but not endorsed by OpenBio."
        )
    inference_unit = parameters.get("inference_unit")
    direction = parameters.get("direction")
    if not isinstance(direction, str) or not direction.strip() or direction != direction.strip():
        warnings.append(
            f"Generic enrichment {artifact_name} does not declare a canonical ranking/selection direction; "
            "direction is reported as 'not_declared'."
        )
        direction = "not_declared"
    batch_handling = parameters.get("technical_batch_handling")
    replicate_aware = parameters.get("replicate_aware")
    declared_sample_level = (
        scope == "condition_contrast"
        and inference_unit == "Sample"
        and replicate_aware is True
        and isinstance(batch_handling, (str, Mapping))
        and bool(batch_handling)
        and parameters.get("condition_inference_level") != "cell_level"
    )
    if scope == "condition_contrast" and not declared_sample_level:
        warnings.append(
            f"Generic enrichment {artifact_name} does not declare a complete Sample-level, replicate-aware "
            "Condition design with Technical-batch handling; computation remains available, but formal Condition "
            "interpretation is not validated."
        )
    return {
        "evidence_scope": scope,
        "inference_unit": _declaration_value(inference_unit),
        "direction": direction,
        "replicate_aware": _declaration_value(replicate_aware),
        "technical_batch_handling": _declaration_value(batch_handling),
        "condition_inference_level": _declaration_value(
            parameters.get("condition_inference_level")
        ),
        "producer_declares_sample_level_condition_design": declared_sample_level,
    }, warnings


def validate_enrichment_artifact_pair(
    table: Any,
    universe: Any,
    *,
    purpose: Literal["ranked", "selected"],
    selector: str,
    gene_column: str,
    score_column: str | None,
    np: Any,
    pd: Any,
) -> dict[str, Any]:
    """Validate ranked/selected structure, current-content fingerprints, and paired axes."""
    if purpose not in {"ranked", "selected"}:
        raise ValueError(f"Unsupported enrichment artifact purpose: {purpose!r}.")
    if not isinstance(table, TableResult) or not isinstance(universe, TableResult):
        raise TypeError("Enrichment consumers require TableResult values for table and universe.")
    if not isinstance(gene_column, str) or not gene_column.strip() or gene_column != gene_column.strip():
        raise ValueError("Enrichment gene_column must be a nonblank string without surrounding whitespace.")
    if purpose == "ranked" and (
        not isinstance(score_column, str) or not score_column.strip() or score_column != score_column.strip()
    ):
        raise ValueError("Ranked enrichment score_column must be a nonblank string without whitespace.")

    table_parameters = _source_parameters(table, artifact_name="Enrichment evidence table")
    universe_parameters = _source_parameters(universe, artifact_name="Enrichment universe")
    marker_family = table_parameters.get("marker_evidence_schema_version") == 2
    if marker_family:
        allowed = ("marker_genes",) if purpose == "ranked" else ("filter_marker_genes",)
        provenance = validate_marker_artifact_pair(
            table,
            universe,
            np=np,
            pd=pd,
            allowed_table_operations=allowed,
        )
        if purpose == "ranked" and provenance["ranking_truncated"]:
            raise ValueError("Ranked enrichment requires a complete, non-truncated upstream ranking.")
        if gene_column != "gene":
            raise ValueError("Marker enrichment evidence requires gene_column='gene'.")
        if purpose == "ranked" and score_column not in MARKER_COLUMNS:
            raise ValueError(f"Marker ranked evidence has no score column {score_column!r}.")
        selected, comparison, comparison_labels = _selected_comparison_frame(
            table.table,
            comparison_column="group",
            selector=selector,
            pd=pd,
        )
        genes = _validated_universe(universe.table, np=np, pd=pd)
        if purpose == "ranked":
            selected_genes = selected["gene"].tolist()
            if len(selected_genes) != len(genes) or set(selected_genes) != set(genes):
                raise ValueError(
                    "Ranked enrichment requires exactly one row for every gene in the paired tested universe."
                )
            values = selected[score_column]
            if pd.api.types.is_bool_dtype(values.dtype) or not pd.api.types.is_numeric_dtype(values.dtype):
                raise TypeError("Ranked enrichment scores must be numeric and non-boolean.")
            scores = values.to_numpy(dtype=float)
            if not bool(np.isfinite(scores).all()) or bool(np.all(scores == scores[0])):
                raise ValueError("Ranked enrichment scores must be finite and nonconstant.")
            if scores.size > 1 and bool((scores[:-1] < scores[1:]).any()):
                raise ValueError("Ranked enrichment rows must preserve a non-increasing upstream score order.")
        return {
            **provenance,
            "artifact_family": "marker_v2",
            "purpose": purpose,
            "comparison": comparison,
            "comparison_labels": comparison_labels,
            "comparison_column": "group",
            "gene_column": "gene",
            "score_column": score_column,
            "selected_frame": selected,
            "universe_genes": genes,
            "evidence_scope": "cluster_marker_evidence",
            "inference_unit": "cells_exploratory",
            "direction": f"{comparison} versus rest; positive scores favor {comparison}",
            "replicate_aware": False,
            "technical_batch_handling": "not modeled by cluster-marker ranking",
            "provenance_warnings": [],
            "provenance_declarations": {
                "validation_scope": "marker_v2_structure_fingerprints_axes_and_direct_operation",
                "scientific_metadata_status": "cluster_marker_inference_is_exploratory",
                "formal_condition_interpretation_validated": False,
                "table_producer_node_id": provenance.get("table_producer_node_id"),
                "universe_producer_node_id": provenance.get("universe_producer_node_id"),
            },
        }

    provenance_warnings: list[str] = []
    declaration_records: dict[str, dict[str, Any]] = {}
    for parameters, name in (
        (table_parameters, "evidence table"),
        (universe_parameters, "tested-gene universe"),
    ):
        if parameters.get("enrichment_evidence_schema_version") != ENRICHMENT_EVIDENCE_SCHEMA_VERSION:
            raise ValueError(f"Generic enrichment {name} uses an unsupported evidence schema.")
        for field in (
            "analysis_fingerprint",
            "ranking_fingerprint",
            "universe_fingerprint",
            "content_fingerprint",
        ):
            if not _valid_sha256(parameters.get(field)):
                raise ValueError(f"Generic enrichment {name} has an invalid SHA-256 {field}.")
        producer = parameters.get("producer_node_id")
        producer_operation = parameters.get("producer_operation")
        source = table.source if name == "evidence table" else universe.source
        approved = parameters.get("enrichment_approved")
        source_operation = source.get("operation")
        declaration_records[name] = {
            "enrichment_approved": _declaration_value(approved),
            "producer_node_id": _declaration_value(producer),
            "producer_operation": _declaration_value(producer_operation),
            "source_operation": _declaration_value(source_operation),
        }
        if approved is not True:
            provenance_warnings.append(
                f"Generic enrichment {name} is not producer-declared enrichment_approved=True; structural and "
                "fingerprint checks passed, but producer suitability is not endorsed."
            )
        if not isinstance(producer, str) or not producer.strip() or producer != producer.strip():
            provenance_warnings.append(
                f"Generic enrichment {name} has no canonical producer node declaration."
            )
        if (
            not isinstance(producer_operation, str)
            or not producer_operation.strip()
            or producer_operation != producer_operation.strip()
        ):
            provenance_warnings.append(
                f"Generic enrichment {name} has no canonical producer operation declaration."
            )
        if source_operation != producer_operation:
            provenance_warnings.append(
                f"Generic enrichment {name} producer operation declaration disagrees with the table source "
                "operation; both are reported as unverified caller metadata."
            )

    expected_role = GENERIC_RANKED_ARTIFACT_ROLE if purpose == "ranked" else GENERIC_SELECTED_ARTIFACT_ROLE
    if table_parameters.get("artifact_role") != expected_role:
        raise ValueError(f"Generic enrichment table has wrong artifact_role for {purpose} enrichment.")
    if universe_parameters.get("artifact_role") != GENERIC_UNIVERSE_ARTIFACT_ROLE:
        raise ValueError("Generic enrichment universe has the wrong artifact_role.")
    for field in ("analysis_fingerprint", "ranking_fingerprint", "universe_fingerprint"):
        if table_parameters[field] != universe_parameters[field]:
            raise ValueError(f"Generic enrichment table and universe disagree on {field}.")
    if declaration_records["evidence table"]["producer_node_id"] != declaration_records[
        "tested-gene universe"
    ]["producer_node_id"]:
        provenance_warnings.append(
            "Generic enrichment table and universe declare different producer node identities; pairing is accepted "
            "only because their structural fingerprints and axes agree."
        )
    if table_parameters.get("evidence_scope") != universe_parameters.get("evidence_scope"):
        provenance_warnings.append(
            "Generic enrichment table and universe declare different evidence_scope values; both declarations are "
            "reported and formal inference is not validated."
        )
    if bool(table_parameters.get("ranking_truncated")) or bool(universe_parameters.get("ranking_truncated")):
        provenance_warnings.append(
            "Generic enrichment metadata declares a truncated upstream ranking. Current table/universe structure is "
            "validated, but upstream completeness is not programmatically proven."
        )
    if purpose == "selected":
        if table_parameters.get("selection_applied") is not True:
            provenance_warnings.append(
                "Selected enrichment evidence does not declare selection_applied=True; the supplied explicit gene "
                "set is analyzed, but its selection procedure is not verified."
            )
        if not _valid_sha256(table_parameters.get("upstream_content_fingerprint")):
            provenance_warnings.append(
                "Selected enrichment evidence has no valid direct-upstream content fingerprint."
            )
        upstream_producer = table_parameters.get("upstream_producer_node_id")
        if (
            not isinstance(upstream_producer, str)
            or not upstream_producer.strip()
            or upstream_producer != upstream_producer.strip()
        ):
            provenance_warnings.append(
                "Selected enrichment evidence has no canonical direct-upstream producer declaration."
            )

    genes = _validated_universe(universe.table, np=np, pd=pd)
    identity = enrichment_universe_identity_fingerprint(genes)
    if table_parameters["universe_fingerprint"] != identity:
        raise ValueError("Generic enrichment table does not match the tested-universe identity.")
    if universe_parameters["universe_fingerprint"] != identity:
        raise ValueError("Generic enrichment universe content does not match its identity fingerprint.")
    universe_content = enrichment_universe_content_fingerprint(universe.table)
    if universe_parameters["content_fingerprint"] != universe_content:
        raise ValueError("Generic enrichment universe content does not match its current-content fingerprint.")
    table_content = generic_enrichment_content_fingerprint(table.table, artifact_role=expected_role)
    if table_parameters["content_fingerprint"] != table_content:
        raise ValueError("Generic enrichment table content does not match its current-content fingerprint.")
    if int(table.input_genes) != len(genes) or int(universe.input_genes) != len(genes):
        raise ValueError("Generic enrichment artifacts report the wrong tested-universe size.")
    if int(table.input_cells) != int(universe.input_cells):
        raise ValueError("Generic enrichment artifacts disagree on upstream cell count.")

    if table_parameters.get("gene_column") != gene_column:
        provenance_warnings.append(
            "Generic enrichment gene_column differs from the producer declaration; the explicitly selected existing "
            "column is used and reported."
        )
    comparison_column = table_parameters.get("comparison_column")
    if not isinstance(comparison_column, str) or not comparison_column.strip():
        raise ValueError("Generic enrichment provenance has no comparison column.")
    if purpose == "ranked" and table_parameters.get("score_column") != score_column:
        provenance_warnings.append(
            "Generic ranked score_column differs from the producer declaration; the explicitly selected existing "
            "column is used and reported."
        )
    selected, comparison, comparison_labels = _selected_comparison_frame(
        table.table,
        comparison_column=comparison_column,
        selector=selector,
        pd=pd,
    )
    if gene_column not in selected:
        raise ValueError(f"Generic enrichment gene column not found: {gene_column!r}.")
    selected_genes = selected[gene_column].tolist()
    for gene in selected_genes:
        if not isinstance(gene, str) or not gene.strip() or gene != gene.strip():
            raise ValueError("Generic enrichment genes must be nonblank strings without whitespace.")
    if len(selected_genes) != len(set(selected_genes)):
        raise ValueError("Generic enrichment genes must be unique within the selected comparison.")
    if not set(selected_genes).issubset(set(genes)):
        raise ValueError("Generic enrichment table contains genes outside its paired universe.")
    if purpose == "ranked":
        if len(selected_genes) != len(genes) or set(selected_genes) != set(genes):
            raise ValueError("Generic ranked evidence must contain exactly every paired-universe gene.")
        values = selected[score_column]
        if pd.api.types.is_bool_dtype(values.dtype) or not pd.api.types.is_numeric_dtype(values.dtype):
            raise TypeError("Generic ranked evidence scores must be numeric and non-boolean.")
        scores = values.to_numpy(dtype=float)
        if not bool(np.isfinite(scores).all()) or bool(np.all(scores == scores[0])):
            raise ValueError("Generic ranked evidence scores must be finite and nonconstant.")
        if scores.size > 1 and bool((scores[:-1] < scores[1:]).any()):
            raise ValueError("Generic ranked evidence rows must preserve non-increasing score order.")
        observed_ranking = generic_ranking_fingerprint(
            table.table,
            comparison_column=comparison_column,
            gene_column=gene_column,
            score_column=score_column,
        )
        if table_parameters["ranking_fingerprint"] != observed_ranking:
            raise ValueError("Generic ranked evidence order/scores do not match the ranking fingerprint.")
    scope, table_scope_warnings = _declared_scope(
        table_parameters, artifact_name="evidence table"
    )
    universe_scope, universe_scope_warnings = _declared_scope(
        universe_parameters, artifact_name="tested-gene universe"
    )
    provenance_warnings.extend(table_scope_warnings)
    provenance_warnings.extend(universe_scope_warnings)
    if scope != universe_scope:
        provenance_warnings.append(
            "Generic enrichment table and universe inference declarations differ; the evidence-table declarations "
            "drive output labels, while neither declaration is programmatically endorsed."
        )
    provenance_warnings.insert(
        0,
        "Generic evidence producer identity, approval, scope, inference unit, replicate awareness, and "
        "Technical-batch handling are caller declarations; OpenBio validated structure, current-content "
        "fingerprints, and table/universe axes only.",
    )
    provenance_warnings = list(dict.fromkeys(provenance_warnings))
    provenance_declarations = {
        "validation_scope": "structure_current_content_fingerprints_and_axes_only",
        "scientific_metadata_status": "caller_declared_not_programmatically_verified",
        "formal_condition_interpretation_validated": False,
        "evidence_table": declaration_records["evidence table"],
        "tested_gene_universe": declaration_records["tested-gene universe"],
        "evidence_table_inference": scope,
        "tested_gene_universe_inference": universe_scope,
    }
    return {
        "artifact_family": "generic_v1",
        "purpose": purpose,
        "analysis_fingerprint": table_parameters["analysis_fingerprint"],
        "ranking_fingerprint": table_parameters["ranking_fingerprint"],
        "universe_fingerprint": identity,
        "table_content_fingerprint": table_content,
        "universe_content_fingerprint": universe_content,
        "table_operation": table.source["operation"],
        "comparison": comparison,
        "comparison_labels": comparison_labels,
        "comparison_column": comparison_column,
        "gene_column": gene_column,
        "score_column": score_column,
        "selected_frame": selected,
        "universe_genes": genes,
        "upstream_parameters": dict(table_parameters),
        "provenance_warnings": provenance_warnings,
        "provenance_declarations": provenance_declarations,
        **scope,
    }


def _standalone_validate_pinned_enrichment_pair(
    table_input,
    universe_input,
    *,
    artifact_family,
    purpose,
    comparison,
    comparison_column,
    gene_column,
    score_column,
    expected_analysis_fingerprint,
    expected_ranking_fingerprint,
    expected_universe_fingerprint,
    expected_table_content_fingerprint,
    expected_universe_content_fingerprint,
):
    """Revalidate current frame bytes against runtime-pinned enrichment identities."""
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd

    operation = "Pinned enrichment evidence"

    for name, value in (
        ("analysis", expected_analysis_fingerprint),
        ("ranking", expected_ranking_fingerprint),
        ("universe", expected_universe_fingerprint),
        ("table content", expected_table_content_fingerprint),
        ("universe content", expected_universe_content_fingerprint),
    ):
        if not _valid_sha256(value):
            raise ValueError(f"{operation} expected {name} SHA-256 is invalid.")

    def unwrap(value, *, name):
        if isinstance(value, pd.DataFrame):
            return value, None, None
        frame = getattr(value, "table", None)
        parameters = getattr(value, "parameters", None)
        source = getattr(value, "source", None)
        if not isinstance(frame, pd.DataFrame):
            raise TypeError(f"{operation} {name} must be a pandas DataFrame or table artifact.")
        if not isinstance(parameters, Mapping) or not isinstance(source, Mapping):
            raise ValueError(f"{operation} {name} has malformed provenance.")
        if not isinstance(source.get("parameters"), Mapping) or dict(source["parameters"]) != dict(parameters):
            raise ValueError(f"{operation} {name} parameters disagree with its immutable source.")
        return frame, dict(parameters), dict(source)

    table, table_parameters, _ = unwrap(table_input, name="table")
    universe, universe_parameters, _ = unwrap(universe_input, name="universe")
    universe_genes = _validated_universe(universe, np=np, pd=pd, operation=operation)
    observed_universe_identity = enrichment_universe_identity_fingerprint(universe_genes)
    if observed_universe_identity != expected_universe_fingerprint:
        raise ValueError(f"{operation} universe identity changed.")

    observed_universe_content = _frame_content_fingerprint(
        universe,
        artifact="tested-gene-universe",
        numeric_columns={"universe_rank"},
    )
    if observed_universe_content != expected_universe_content_fingerprint:
        raise ValueError(f"{operation} universe current content changed.")

    if comparison_column not in table or gene_column not in table:
        raise ValueError(f"{operation} table is missing comparison/gene columns.")
    selected, _, _ = _selected_comparison_frame(
        table,
        comparison_column=comparison_column,
        selector=comparison,
        pd=pd,
        operation=operation,
    )
    selected_genes = selected[gene_column].tolist()
    for gene in selected_genes:
        if not isinstance(gene, str) or not gene.strip() or gene != gene.strip():
            raise ValueError(f"{operation} selected genes are invalid.")
    if len(selected_genes) != len(set(selected_genes)):
        raise ValueError(f"{operation} selected comparison contains duplicate genes.")
    if not set(selected_genes).issubset(set(universe_genes)):
        raise ValueError(f"{operation} table contains genes outside the paired universe.")

    if artifact_family == "marker_v2":
        if list(table.columns) != MARKER_COLUMNS or comparison_column != "group" or gene_column != "gene":
            raise ValueError(f"{operation} marker table schema changed.")
        numeric_columns = set(MARKER_COLUMNS[2:])
        for column in numeric_columns:
            values = table[column]
            if pd.api.types.is_bool_dtype(values.dtype) or not pd.api.types.is_numeric_dtype(values.dtype):
                raise TypeError(f"{operation} marker numeric column {column!r} is invalid.")
            if not bool(np.isfinite(values.to_numpy(dtype=float)).all()):
                raise ValueError(f"{operation} marker numeric column {column!r} is non-finite.")
        observed_table_content = _frame_content_fingerprint(
            table,
            artifact="marker-table",
            numeric_columns=numeric_columns,
        )
    elif artifact_family == "generic_v1":
        expected_role = "complete_ranked_evidence" if purpose == "ranked" else "selected_gene_evidence"
        observed_table_content = generic_enrichment_content_fingerprint(
            table,
            artifact_role=expected_role,
            operation=operation,
        )
    else:
        raise ValueError(f"{operation} artifact family is unsupported: {artifact_family!r}.")
    if observed_table_content != expected_table_content_fingerprint:
        raise ValueError(f"{operation} table current content changed.")

    if purpose == "ranked":
        if score_column not in selected:
            raise ValueError(f"{operation} score column is missing.")
        scores = selected[score_column]
        if pd.api.types.is_bool_dtype(scores.dtype) or not pd.api.types.is_numeric_dtype(scores.dtype):
            raise TypeError(f"{operation} scores must be numeric and non-boolean.")
        numeric_scores = scores.to_numpy(dtype=float)
        if (
            not bool(np.isfinite(numeric_scores).all())
            or bool(np.all(numeric_scores == numeric_scores[0]))
            or (numeric_scores.size > 1 and bool((numeric_scores[:-1] < numeric_scores[1:]).any()))
        ):
            raise ValueError(f"{operation} complete ranking must be finite, nonconstant, and non-increasing.")
        if len(selected_genes) != len(universe_genes) or set(selected_genes) != set(universe_genes):
            raise ValueError(f"{operation} complete ranking does not contain the exact universe.")

    if table_parameters is not None:
        for field, expected in (
            ("analysis_fingerprint", expected_analysis_fingerprint),
            ("ranking_fingerprint", expected_ranking_fingerprint),
            ("universe_fingerprint", expected_universe_fingerprint),
            ("content_fingerprint", expected_table_content_fingerprint),
        ):
            if table_parameters.get(field) != expected:
                raise ValueError(f"{operation} table provenance {field} changed.")
    if universe_parameters is not None:
        for field, expected in (
            ("analysis_fingerprint", expected_analysis_fingerprint),
            ("ranking_fingerprint", expected_ranking_fingerprint),
            ("universe_fingerprint", expected_universe_fingerprint),
            ("content_fingerprint", expected_universe_content_fingerprint),
        ):
            if universe_parameters.get(field) != expected:
                raise ValueError(f"{operation} universe provenance {field} changed.")
    return selected, universe_genes


def pinned_enrichment_validation_code() -> str:
    functions = (
        _canonical_sha256,
        _valid_sha256,
        _json_cell,
        generic_enrichment_content_fingerprint,
        enrichment_universe_identity_fingerprint,
        _validated_universe,
        _selected_comparison_frame,
        _frame_content_fingerprint,
        _standalone_validate_pinned_enrichment_pair,
    )
    constants = (
        f"_SHA256_LENGTH = {_SHA256_LENGTH!r}",
        f"MARKER_COLUMNS = {MARKER_COLUMNS!r}",
        f"MARKER_UNIVERSE_COLUMNS = {MARKER_UNIVERSE_COLUMNS!r}",
    )
    return "\n\n".join(
        (
            "import hashlib\nimport json\nimport math\nimport numbers\nfrom typing import Any",
            "\n".join(constants),
            *(textwrap.dedent(inspect.getsource(function)).strip() for function in functions),
        )
    )


__all__ = [
    "ENRICHMENT_EVIDENCE_SCHEMA_VERSION",
    "EVIDENCE_SCOPES",
    "GENERIC_RANKED_ARTIFACT_ROLE",
    "GENERIC_SELECTED_ARTIFACT_ROLE",
    "GENERIC_UNIVERSE_ARTIFACT_ROLE",
    "enrichment_universe_content_fingerprint",
    "enrichment_universe_identity_fingerprint",
    "generic_enrichment_content_fingerprint",
    "generic_ranking_fingerprint",
    "pinned_enrichment_validation_code",
    "validate_enrichment_artifact_pair",
]
