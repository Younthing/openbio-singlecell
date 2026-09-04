from __future__ import annotations

import inspect
import textwrap
from typing import Any


def _standalone_map_cluster_annotations(
    adata,
    *,
    groupby,
    mapping_json,
    output_column,
    unmapped_policy="error",
    annotation_status="provisional",
    overwrite_existing=False,
    _return_diagnostics=False,
):
    """Apply one strict cluster-label map without depending on OpenBio runtime helpers."""
    import copy
    import json
    import math
    import numbers
    from importlib import metadata as importlib_metadata

    import numpy as np
    import pandas as pd

    operation = "Map Cluster Annotations"
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires an in-memory AnnData; call to_memory() first.")
    if int(getattr(adata, "n_obs", 0)) < 1 or int(getattr(adata, "n_vars", 0)) < 1:
        raise ValueError(f"{operation} requires non-empty observation and feature axes.")
    if not bool(adata.obs_names.is_unique) or not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique observation and feature identifiers.")
    if not isinstance(groupby, str) or not groupby.strip():
        raise ValueError("Cluster groupby column cannot be empty.")
    if not isinstance(output_column, str) or not output_column.strip():
        raise ValueError("Cluster annotation output column cannot be empty.")
    groupby = groupby.strip()
    output_column = output_column.strip()
    if groupby not in adata.obs:
        raise ValueError(f"Cluster groupby column not found in obs: {groupby!r}")
    if output_column == groupby:
        raise ValueError("Cluster annotation output_column cannot overwrite the groupby evidence column.")
    if not isinstance(overwrite_existing, bool):
        raise TypeError("Cluster annotation overwrite_existing must be a boolean.")
    output_existed = output_column in adata.obs
    if output_existed and not overwrite_existing:
        raise ValueError(
            f"Cluster annotation output column already exists: {output_column!r}; "
            "enable overwrite_existing to replace it."
        )
    if unmapped_policy not in {"error", "preserve_cluster_label", "set_missing"}:
        raise ValueError(f"Unsupported cluster annotation unmapped_policy: {unmapped_policy!r}")
    if annotation_status not in {"provisional", "curated"}:
        raise ValueError(f"Unsupported cluster annotation annotation_status: {annotation_status!r}")
    if not isinstance(mapping_json, str):
        raise TypeError("Cluster annotation mapping_json must be a string.")
    if len(mapping_json.encode("utf-8")) > 1_000_000:
        raise ValueError("Cluster annotation mapping_json exceeds the 1,000,000-byte safety limit.")

    def reject_constant(value):
        raise ValueError(f"Cluster annotation mapping contains non-standard JSON constant {value!r}.")

    def object_from_pairs(pairs):
        if len(pairs) > 10_000:
            raise ValueError("Cluster annotation mapping exceeds the 10,000-entry safety limit.")
        result = {}
        raw_seen = set()
        for raw_key, value in pairs:
            if raw_key in raw_seen:
                raise ValueError(f"Cluster annotation mapping contains duplicate key {raw_key!r}.")
            raw_seen.add(raw_key)
            key = raw_key.strip()
            if not key:
                raise ValueError("Cluster annotation mapping keys must be non-empty strings.")
            if key in result:
                raise ValueError(f"Cluster annotation mapping keys collide after trimming: {raw_key!r} -> {key!r}.")
            result[key] = value
        return result

    try:
        raw_mapping = json.loads(
            mapping_json,
            object_pairs_hook=object_from_pairs,
            parse_constant=reject_constant,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"Cluster annotation mapping is not valid JSON ({exc.msg}).") from exc
    if not isinstance(raw_mapping, dict) or not raw_mapping:
        raise ValueError("Cluster annotation mapping must be a non-empty JSON object.")
    mapping = {}
    for key, value in raw_mapping.items():
        if not isinstance(value, str):
            raise TypeError(f"Cluster annotation target for {key!r} must be a string, not {type(value).__name__}.")
        target = value.strip()
        if not target:
            raise ValueError(f"Cluster annotation target for {key!r} cannot be blank.")
        mapping[key] = target

    source = adata.obs[groupby]
    source_is_categorical = isinstance(source.dtype, pd.CategoricalDtype)
    source_ordered = bool(source.cat.ordered) if source_is_categorical else False
    canonical_signatures = {}

    def canonicalize(value, *, context):
        if not bool(pd.api.types.is_scalar(value)):
            raise TypeError(f"{operation} {context} contains a nonscalar source label.")
        missing = pd.isna(value)
        if not isinstance(missing, (bool, np.bool_)):
            raise TypeError(f"{operation} {context} contains a nonscalar source label.")
        if bool(missing):
            return None
        if isinstance(value, numbers.Real) and not math.isfinite(float(value)):
            raise ValueError(f"{operation} {context} contains a non-finite source label.")
        key = str(value).strip()
        if not key:
            raise ValueError(f"{operation} {context} contains a blank source label.")
        signature = (type(value).__module__, type(value).__qualname__, repr(value))
        previous = canonical_signatures.get(key)
        if previous is not None and previous != signature:
            raise ValueError(f"{operation} has distinct source labels that collapse to canonical key {key!r}.")
        canonical_signatures[key] = signature
        return key

    declared_levels = []
    if source_is_categorical:
        for category in source.cat.categories.tolist():
            key = canonicalize(category, context="declared categories")
            if key is None:
                raise ValueError("Cluster annotation categorical levels cannot be missing.")
            declared_levels.append(key)

    row_levels = []
    source_counts = {}
    missing_source_cells = 0
    for value in source.tolist():
        key = canonicalize(value, context="observations")
        row_levels.append(key)
        if key is None:
            missing_source_cells += 1
        else:
            source_counts[key] = source_counts.get(key, 0) + 1
    if not source_counts:
        raise ValueError(f"{operation} requires at least one non-missing observed source level.")

    observed_levels = (
        [level for level in declared_levels if level in source_counts]
        if source_is_categorical
        else sorted(source_counts)
    )
    unused_declared_levels = [level for level in declared_levels if level not in source_counts]
    known_levels = set(declared_levels) if source_is_categorical else set(observed_levels)
    unknown_mapping_keys = sorted(set(mapping) - known_levels)
    if unknown_mapping_keys:
        raise ValueError(
            "Cluster annotation mapping contains keys unknown to observed or declared source levels: "
            f"{unknown_mapping_keys}."
        )
    effective_mapping = {level: mapping[level] for level in observed_levels if level in mapping}
    if not effective_mapping:
        raise ValueError("Cluster annotation mapping must map at least one observed source level.")
    unused_mapping = {level: mapping[level] for level in unused_declared_levels if level in mapping}
    unmapped_levels = [level for level in observed_levels if level not in mapping]
    if unmapped_levels and unmapped_policy == "error":
        counts = {level: source_counts[level] for level in unmapped_levels}
        raise ValueError(f"Cluster annotation has observed unmapped levels: {counts}.")
    if unmapped_policy == "preserve_cluster_label":
        collisions = sorted(set(effective_mapping.values()) & set(unmapped_levels))
        if collisions:
            raise ValueError(f"Cluster annotation mapped targets collide with preserved source labels: {collisions}.")

    resolved_by_level = {}
    for level in observed_levels:
        if level in mapping:
            resolved_by_level[level] = mapping[level]
        elif unmapped_policy == "preserve_cluster_label":
            resolved_by_level[level] = level
        else:
            resolved_by_level[level] = None
    output_categories = []
    for level in observed_levels:
        target = resolved_by_level[level]
        if target is not None and target not in output_categories:
            output_categories.append(target)
    output_values = [None if level is None else resolved_by_level[level] for level in row_levels]
    output_categorical = pd.Categorical(output_values, categories=output_categories, ordered=False)

    reverse_targets = {}
    for level, target in effective_mapping.items():
        reverse_targets.setdefault(target, []).append(level)
    many_to_one = {target: levels for target, levels in reverse_targets.items() if len(levels) > 1}
    identity_mappings = [level for level, target in effective_mapping.items() if level == target]
    mapped_cells = int(sum(source_counts[level] for level in effective_mapping))
    unmapped_cells = int(sum(source_counts[level] for level in unmapped_levels))
    preserved_cells = unmapped_cells if unmapped_policy == "preserve_cluster_label" else 0
    set_missing_cells = unmapped_cells if unmapped_policy == "set_missing" else 0
    output_counts = {category: int(sum(value == category for value in output_values)) for category in output_categories}
    provenance = {
        "schema_version": 1,
        "operation": "map_cluster_annotations",
        "annotation_status": annotation_status,
        "curation_assertion": "caller_declared" if annotation_status == "curated" else None,
        "software_verified_review": False,
        "groupby": groupby,
        "output_column": output_column,
        "unmapped_policy": unmapped_policy,
        "overwrite_existing": overwrite_existing,
        "output_existed": output_existed,
        "effective_mapping": effective_mapping,
        "unused_declared_level_mapping": unused_mapping,
        "source_dtype": str(source.dtype),
        "source_ordered": source_ordered,
        "declared_source_categories": declared_levels,
        "observed_source_levels": observed_levels,
        "unused_declared_levels": unused_declared_levels,
        "source_counts": {level: int(source_counts[level]) for level in observed_levels},
        "missing_source_cells": missing_source_cells,
        "unmapped_source_levels": unmapped_levels,
        "mapped_cells": mapped_cells,
        "unmapped_cells": unmapped_cells,
        "preserved_cells": preserved_cells,
        "set_missing_cells": set_missing_cells,
        "output_categories": output_categories,
        "output_counts": output_counts,
        "output_ordered": False,
        "identity_mappings": identity_mappings,
        "many_to_one_merges": many_to_one,
    }
    diagnostics = {
        **provenance,
        "mapping_entries": len(mapping),
        "effective_mapping_entries": len(effective_mapping),
        "unused_mapping_entries": len(unused_mapping),
        "observed_level_count": len(observed_levels),
        "output_missing_cells": int(missing_source_cells + set_missing_cells),
    }

    existing_metadata = adata.uns.get("openbio_singlecell")
    if existing_metadata is not None:
        if not isinstance(existing_metadata, dict):
            raise ValueError("OpenBio metadata must be a mapping.")
        if existing_metadata.get("schema_version") != 1:
            raise ValueError(
                "OpenBio metadata has unsupported schema_version "
                f"{existing_metadata.get('schema_version')!r}; expected 1."
            )
    output = adata
    output.obs[output_column] = output_categorical
    metadata = (
        copy.deepcopy(existing_metadata)
        if existing_metadata is not None
        else {
            "schema_version": 1,
            "version": "not-installed",
            "display_name": "AnnData",
            "source": {},
            "random_seed": 0,
            "warnings": [],
            "analysis_history": {},
        }
    )
    if existing_metadata is None:
        try:
            metadata["version"] = importlib_metadata.version("openbio-singlecell")
        except importlib_metadata.PackageNotFoundError:
            pass
    annotations = metadata.get("annotations", {})
    if not isinstance(annotations, dict):
        raise ValueError("OpenBio annotation provenance must be a mapping.")
    annotations = copy.deepcopy(annotations)
    annotations[output_column] = provenance
    metadata["annotations"] = annotations
    output.uns["openbio_singlecell"] = metadata
    if not output.obs_names.equals(adata.obs_names) or not output.var_names.equals(adata.var_names):
        raise RuntimeError(f"{operation} changed observation or feature identity.")
    if not output.obs[groupby].equals(source):
        raise RuntimeError(f"{operation} changed the source grouping evidence.")
    stored_output = output.obs[output_column]
    if not isinstance(stored_output.dtype, pd.CategoricalDtype):
        raise RuntimeError(f"{operation} failed to store a categorical annotation.")
    if stored_output.cat.categories.tolist() != output_categories or bool(stored_output.cat.ordered):
        raise RuntimeError(f"{operation} stored the wrong annotation category order.")
    if int(stored_output.isna().sum()) != diagnostics["output_missing_cells"]:
        raise RuntimeError(f"{operation} stored the wrong number of missing annotation values.")
    if _return_diagnostics:
        return output, diagnostics
    return output


def _standalone_resolve_celltypist_model_path(model_identifier, models_module):
    """Resolve one explicit or locally installed CellTypist model without downloading."""
    from pathlib import Path

    if not isinstance(model_identifier, str) or not model_identifier.strip():
        raise ValueError("CellTypist model cannot be empty.")
    model_identifier = model_identifier.strip()
    candidate = Path(model_identifier).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    get_model_path = getattr(models_module, "get_model_path", None)
    if not callable(get_model_path):
        raise RuntimeError("CellTypist backend does not expose the official local model-path API.")
    try:
        resolved_model_path = Path(get_model_path(model_identifier)).expanduser().resolve()
    except Exception as exc:
        raise FileNotFoundError(
            f"CellTypist could not resolve local model {model_identifier!r}. The annotation node never "
            "downloads models. Explicitly run the network-enabled official command "
            f"celltypist.models.download_models(model={model_identifier!r}) before execution, or provide "
            "a trusted local .pkl path."
        ) from exc
    if not resolved_model_path.is_file():
        raise FileNotFoundError(
            f"CellTypist model {model_identifier!r} is not an existing path or installed local model; "
            f"expected local file {str(resolved_model_path)!r}. The annotation node never downloads models. "
            "Explicitly run the network-enabled official command "
            f"celltypist.models.download_models(model={model_identifier!r}) before execution, or provide "
            "a trusted local .pkl path."
        )
    return resolved_model_path


def _standalone_celltypist_software_versions(*, openbio_version):
    """Collect the version fields used by the standalone CellTypist report."""
    import platform
    from importlib import metadata as importlib_metadata

    versions = {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
    }
    for package in (
        "celltypist",
        "scikit-learn",
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
    ):
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _standalone_celltypist_summary(diagnostics):
    """Build the strict scientific CellTypist summary without mutating inputs."""
    import json
    import math
    from collections.abc import Mapping as RuntimeMapping

    import numpy as np

    if not isinstance(diagnostics, RuntimeMapping):
        raise TypeError("CellTypist summary diagnostics must be a mapping.")
    def plain_json(value):
        if value is None or isinstance(value, (str, bool, int)):
            return value
        if isinstance(value, float):
            return value if math.isfinite(value) else None
        if hasattr(value, "item"):
            try:
                return plain_json(value.item())
            except (TypeError, ValueError):
                pass
        if isinstance(value, RuntimeMapping):
            return {str(key): plain_json(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
            sequence = value.tolist() if hasattr(value, "tolist") else value
            return [plain_json(item) for item in sequence]
        return str(value)

    def summarize_numeric(values):
        array = np.asarray(values, dtype=float).ravel()
        finite = array[np.isfinite(array)]
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
        quantiles = np.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])
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

    required = {
        "cells",
        "expression_source",
        "expression_state",
        "expression_transform",
        "selected_label_counts",
        "selected_probabilities",
        "individual_probabilities",
        "majority_support",
        "warnings",
        "parameters",
        "software_versions",
    }
    missing = sorted(required - set(diagnostics))
    if missing:
        raise ValueError(f"CellTypist summary diagnostics are missing required fields: {missing}.")
    if not isinstance(diagnostics["software_versions"], RuntimeMapping):
        raise TypeError("CellTypist summary software_versions must be a mapping.")

    parameters = plain_json(diagnostics["parameters"])
    if not isinstance(parameters, dict):
        raise TypeError("CellTypist summary parameters must be a mapping.")
    majority_voting = parameters.get("majority_voting")
    if not isinstance(majority_voting, bool):
        raise TypeError("CellTypist summary majority_voting must be a boolean.")
    warnings = [str(value) for value in diagnostics["warnings"]]
    bounded_diagnostics = {
        key: plain_json(value)
        for key, value in diagnostics.items()
        if key
        not in {
            "selected_probabilities",
            "individual_probabilities",
            "majority_support",
            "parameters",
        }
    }
    key_results = {
        **bounded_diagnostics,
        "selected_probability_summary": summarize_numeric(diagnostics["selected_probabilities"]),
        "individual_probability_summary": summarize_numeric(diagnostics["individual_probabilities"]),
        "majority_support_summary": (
            summarize_numeric(diagnostics["majority_support"]) if majority_voting else None
        ),
        "annotation_status": "provisional",
        "selected_probability_semantics": (
            "probability assigned to the selected model class; missing for Heterogeneous"
            if majority_voting
            else "row-maximum class probability for the individual best match"
        ),
        "majority_support_semantics": (
            "dominant individual-label fraction within the explicit over-cluster" if majority_voting else None
        ),
    }
    cells = int(diagnostics["cells"])
    methods = (
        f"A pinned CellTypist best-match model was applied to {diagnostics['expression_source']} declared "
        f"as {diagnostics['expression_state']!r}; the matrix transform was "
        f"{diagnostics['expression_transform']!r}. Majority voting was "
        f"{'performed on explicit groups' if majority_voting else 'disabled'}."
    )
    model_class_count = len(diagnostics.get("class_order", []))
    heterogeneous_disclosure = (
        " plus the observed Heterogeneous sentinel"
        if "Heterogeneous" in diagnostics["selected_label_counts"]
        else ""
    )
    results = (
        f"CellTypist produced Provisional annotation for {cells:,} cells with a complete count map for "
        f"all {model_class_count:,} model classes{heterogeneous_disclosure}. These model-derived labels "
        "require marker and biological-context review before curation."
    )
    references = [
        {
            "citation": (
                "Domínguez Conde C, Xu C, Jarvis LB, et al. Cross-tissue immune cell analysis reveals "
                "tissue-specific features in humans. Science. 2022;376:eabl5197."
            ),
            "url": "https://doi.org/10.1126/science.abl5197",
            "kind": "method",
            "doi": "10.1126/science.abl5197",
        },
        {
            "citation": "CellTypist developers. celltypist.annotate official API documentation.",
            "url": "https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html",
            "kind": "software_documentation",
            "doi": None,
        },
        {
            "citation": "CellTypist developers. Official classifier and AnnotationResult implementation.",
            "url": "https://celltypist.readthedocs.io/en/stable/_modules/celltypist/classifier.html",
            "kind": "software_documentation",
            "doi": None,
        },
        {
            "citation": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.",
            "url": "https://doi.org/10.1186/s13059-017-1382-0",
            "kind": "software",
            "doi": "10.1186/s13059-017-1382-0",
        },
        {
            "citation": (
                "Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. "
                "Molecular Systems Biology. 2019;15:e8746."
            ),
            "url": "https://doi.org/10.15252/msb.20188746",
            "kind": "practice",
            "doi": "10.15252/msb.20188746",
        },
    ]
    limitations = [
        "The model vocabulary is closed; unseen or out-of-domain states can be forced to a known class.",
        "CellTypist sigmoid scores are model evidence, not calibrated biological certainty or ground truth.",
        "Majority voting, when enabled, can erase rare, transitional, or mixed states within a supplied group.",
        "Tissue, species, protocol, disease, Technical batch, training composition, and feature overlap can shift predictions.",
        "Cell-level Provisional annotation is not an independent-replicate Condition contrast; formal inference remains Sample-level.",
    ]
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellCellTypistAnnotation",
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": limitations,
        "references": references,
        "software_versions": plain_json(diagnostics["software_versions"]),
    }
    json.dumps(summary, allow_nan=False)
    return summary


def _celltypist_axis_fingerprint(observation_ids):
    import hashlib
    import json

    values = list(observation_ids)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in values):
        raise ValueError("CellTypist observation identifiers must be canonical strings.")
    payload = json.dumps(
        values,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _celltypist_probability_fingerprint(matrix, *, observation_ids, classes):
    import hashlib
    import json

    import numpy as np

    observations = list(observation_ids)
    class_axis = list(classes)
    values = np.asarray(matrix)
    if values.shape != (len(observations), len(class_axis)):
        raise ValueError("CellTypist probability matrix does not align to its observation and class axes.")
    if (
        not np.issubdtype(values.dtype, np.number)
        or np.issubdtype(values.dtype, np.bool_)
        or np.iscomplexobj(values)
    ):
        raise TypeError("CellTypist probabilities must contain real numeric values.")
    if not bool(np.isfinite(values).all()):
        raise ValueError("CellTypist probabilities must be finite.")
    header = json.dumps(
        {"observation_ids": observations, "classes": class_axis},
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    digest = hashlib.sha256(header)
    digest.update(np.ascontiguousarray(values, dtype="<f8").tobytes())
    return digest.hexdigest()


def _celltypist_selected_fingerprint(labels, probabilities, *, observation_ids):
    import hashlib
    import json
    import math

    observations = list(observation_ids)
    selected_labels = list(labels)
    selected_probabilities = list(probabilities)
    if len(selected_labels) != len(observations) or len(selected_probabilities) != len(observations):
        raise ValueError("CellTypist selected annotation does not align to the observation axis.")
    rows = []
    for label, probability in zip(selected_labels, selected_probabilities, strict=True):
        if not isinstance(label, str) or not label:
            raise ValueError("CellTypist selected labels must be nonempty strings.")
        value = float(probability)
        if not (math.isfinite(value) or math.isnan(value)):
            raise ValueError("CellTypist selected probabilities must be finite or missing.")
        rows.append([label, None if math.isnan(value) else value.hex()])
    payload = json.dumps(
        {"observation_ids": observations, "rows": rows},
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _standalone_celltypist_annotation(
    adata,
    *,
    source_kind,
    layer_name,
    expression_state,
    model_identifier,
    requested_model_identifier=None,
    expected_model_sha256=None,
    majority_voting=False,
    over_clustering_key="",
    min_prop=0.0,
    label_column="celltypist_cell_type",
    confidence_column="celltypist_confidence",
    probability_key="celltypist_probabilities",
    store_decision_matrix=False,
    decision_key="celltypist_decision_scores",
    metadata_key="celltypist",
    overwrite_existing=False,
    openbio_version="not-installed",
    celltypist_module=None,
    _return_diagnostics=False,
):
    """Run one pinned CellTypist model with explicit expression-state semantics."""
    import copy
    import hashlib
    import importlib
    import math
    from collections.abc import Mapping as RuntimeMapping
    from importlib import metadata as importlib_metadata

    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    operation = "CellTypist Annotation"
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires an in-memory AnnData; call to_memory() first.")
    if int(getattr(adata, "n_obs", 0)) < 1 or int(getattr(adata, "n_vars", 0)) < 1:
        raise ValueError(f"{operation} requires non-empty observation and feature axes.")
    if not bool(adata.obs_names.is_unique) or not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique observation and feature identifiers.")
    for observation in adata.obs_names.tolist():
        if not isinstance(observation, str) or not observation.strip() or observation != observation.strip():
            raise ValueError(f"{operation} requires non-empty observation identifiers without surrounding whitespace.")
    if source_kind not in {"X", "raw", "layer"}:
        raise ValueError(f"Unsupported CellTypist expression source: {source_kind!r}")
    if expression_state not in {"verified_counts", "verified_cp10k_log1p"}:
        raise ValueError(f"Unsupported CellTypist expression_state: {expression_state!r}")
    if not isinstance(model_identifier, str) or not model_identifier.strip():
        raise ValueError("CellTypist model cannot be empty.")
    model_identifier = model_identifier.strip()
    requested_model_identifier = (
        model_identifier if requested_model_identifier is None else str(requested_model_identifier).strip()
    )
    if not requested_model_identifier:
        raise ValueError("CellTypist requested model identifier cannot be empty.")
    if not isinstance(majority_voting, bool):
        raise TypeError("CellTypist majority_voting must be a boolean.")
    if isinstance(min_prop, bool) or not isinstance(min_prop, (int, float)):
        raise TypeError("CellTypist min_prop must be a finite number in [0, 1].")
    min_prop = float(min_prop)
    if not math.isfinite(min_prop) or not 0.0 <= min_prop <= 1.0:
        raise ValueError("CellTypist min_prop must be a finite number in [0, 1].")
    if not isinstance(over_clustering_key, str):
        raise TypeError("CellTypist over_clustering_key must be a string.")
    over_clustering_key = over_clustering_key.strip()
    if majority_voting:
        if int(adata.n_obs) <= 50:
            raise ValueError(
                "CellTypist explicit majority voting requires more than 50 cells; "
                "the official backend skips voting for smaller inputs."
            )
        if not over_clustering_key:
            raise ValueError("CellTypist majority voting requires an explicit categorical over_clustering_key.")
    elif over_clustering_key:
        raise ValueError("CellTypist over_clustering_key is only valid when majority_voting is enabled.")
    elif min_prop != 0.0:
        raise ValueError("CellTypist min_prop is only valid when majority_voting is enabled.")
    if not isinstance(store_decision_matrix, bool):
        raise TypeError("CellTypist store_decision_matrix must be a boolean.")
    if not isinstance(overwrite_existing, bool):
        raise TypeError("CellTypist overwrite_existing must be a boolean.")
    if not isinstance(openbio_version, str) or not openbio_version.strip():
        raise ValueError("CellTypist openbio_version cannot be empty.")
    openbio_version = openbio_version.strip()

    key_values = {
        "label_column": label_column,
        "confidence_column": confidence_column,
        "probability_key": probability_key,
        "decision_key": decision_key,
        "metadata_key": metadata_key,
    }
    for description, value in key_values.items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"CellTypist {description} cannot be empty.")
        key_values[description] = value.strip()
    label_column = key_values["label_column"]
    confidence_column = key_values["confidence_column"]
    probability_key = key_values["probability_key"]
    decision_key = key_values["decision_key"]
    metadata_key = key_values["metadata_key"]
    if metadata_key == "openbio_singlecell":
        raise ValueError("CellTypist metadata_key cannot use the reserved OpenBio metadata key.")
    individual_label_column = "celltypist_individual_label"
    individual_probability_column = "celltypist_individual_probability"
    majority_label_column = "celltypist_majority_label"
    majority_support_column = "celltypist_majority_support"
    obs_output_keys = [
        label_column,
        confidence_column,
        individual_label_column,
        individual_probability_column,
    ]
    if majority_voting:
        obs_output_keys.extend([majority_label_column, majority_support_column])
    if len(obs_output_keys) != len(set(obs_output_keys)):
        raise ValueError(f"CellTypist observation output keys must be distinct: {obs_output_keys}.")
    if majority_voting and over_clustering_key in obs_output_keys:
        raise ValueError("CellTypist output columns cannot overwrite the over-clustering evidence column.")
    if probability_key == decision_key:
        raise ValueError("CellTypist probability and decision matrix keys must be different.")
    obsm_output_keys = [probability_key] + ([decision_key] if store_decision_matrix else [])
    obs_collisions = [key for key in obs_output_keys if key in adata.obs]
    obsm_collisions = [key for key in obsm_output_keys if key in adata.obsm]
    uns_collisions = [metadata_key] if metadata_key in adata.uns else []
    prior_provenance = None
    prior_columns = None
    prior_matrices = None
    prior_owned_obs_keys = []
    prior_owned_obsm_keys = []
    stale_obs_keys = []
    stale_obsm_keys = []
    stale_annotation_keys = []
    prior_artifact_verified = False
    if obs_collisions or obsm_collisions or uns_collisions:
        if not overwrite_existing:
            raise ValueError(
                "CellTypist output keys already exist; enable overwrite_existing to replace them: "
                f"obs={obs_collisions}, obsm={obsm_collisions}, uns={uns_collisions}."
            )
        prior_provenance = adata.uns.get(metadata_key)
        if (
            not isinstance(prior_provenance, RuntimeMapping)
            or prior_provenance.get("schema_version") != 1
            or prior_provenance.get("operation") != "celltypist_annotation"
        ):
            raise ValueError(
                "CellTypist overwrite_existing requires prior schema-version-1 CellTypist provenance at "
                f"uns[{metadata_key!r}]; ownership of existing outputs could not be verified."
            )
        prior_columns = prior_provenance.get("columns")
        prior_matrices = prior_provenance.get("matrices")
        required_column_roles = {
            "selected_label",
            "selected_label_probability",
            "individual_label",
            "individual_row_max_probability",
            "majority_label",
            "majority_support",
        }
        required_matrix_roles = {
            "probability_key",
            "probability_shape",
            "decision_key",
            "decision_shape",
        }
        if not isinstance(prior_columns, RuntimeMapping) or set(prior_columns) != required_column_roles:
            raise ValueError(
                "CellTypist overwrite_existing cannot verify prior observation-column ownership from provenance."
            )
        if not isinstance(prior_matrices, RuntimeMapping) or set(prior_matrices) != required_matrix_roles:
            raise ValueError(
                "CellTypist overwrite_existing cannot verify prior matrix ownership from provenance."
            )
        required_prior_columns = (
            "selected_label",
            "selected_label_probability",
            "individual_label",
            "individual_row_max_probability",
        )
        if any(prior_columns[role] is None for role in required_prior_columns):
            raise ValueError(
                "CellTypist overwrite_existing found incomplete prior required observation-column ownership."
            )
        if (prior_columns["majority_label"] is None) != (prior_columns["majority_support"] is None):
            raise ValueError(
                "CellTypist overwrite_existing found inconsistent prior majority-column ownership."
            )
        if prior_matrices["probability_key"] is None:
            raise ValueError("CellTypist overwrite_existing found missing prior probability-matrix ownership.")
        if (prior_matrices["decision_key"] is None) != (prior_matrices["decision_shape"] is None):
            raise ValueError(
                "CellTypist overwrite_existing found inconsistent prior decision-matrix ownership."
            )

        def verified_matrix_shape(value, *, role):
            if hasattr(value, "tolist"):
                value = value.tolist()
            if (
                not isinstance(value, (list, tuple))
                or len(value) != 2
                or any(
                    isinstance(item, (bool, np.bool_))
                    or not isinstance(item, (int, np.integer))
                    or int(item) < 1
                    for item in value
                )
            ):
                raise ValueError(
                    f"CellTypist overwrite_existing found invalid prior matrix shape for {role!r}."
                )
            return tuple(int(item) for item in value)

        prior_probability_shape = verified_matrix_shape(
            prior_matrices["probability_shape"],
            role="probability_shape",
        )
        if prior_matrices["decision_key"] is not None:
            prior_decision_shape = verified_matrix_shape(
                prior_matrices["decision_shape"],
                role="decision_shape",
            )
        else:
            prior_decision_shape = None

        def verified_owned_keys(mapping, roles, *, description):
            owned = []
            for role in roles:
                value = mapping[role]
                if value is None:
                    continue
                if not isinstance(value, str) or not value.strip() or value != value.strip():
                    raise ValueError(
                        f"CellTypist overwrite_existing found invalid prior {description} ownership for {role!r}."
                    )
                owned.append(value)
            if len(owned) != len(set(owned)):
                raise ValueError(f"CellTypist overwrite_existing found duplicate prior {description} ownership.")
            return owned

        prior_owned_obs_keys = verified_owned_keys(
            prior_columns,
            (
                "selected_label",
                "selected_label_probability",
                "individual_label",
                "individual_row_max_probability",
                "majority_label",
                "majority_support",
            ),
            description="observation-column",
        )
        prior_owned_obsm_keys = verified_owned_keys(
            prior_matrices,
            ("probability_key", "decision_key"),
            description="matrix",
        )
        if prior_matrices["probability_key"] in adata.obsm and tuple(
            adata.obsm[prior_matrices["probability_key"]].shape
        ) != prior_probability_shape:
            raise ValueError(
                "CellTypist overwrite_existing found a prior probability matrix whose shape does not match provenance."
            )
        if (
            prior_matrices["decision_key"] is not None
            and prior_matrices["decision_key"] in adata.obsm
            and tuple(adata.obsm[prior_matrices["decision_key"]].shape) != prior_decision_shape
        ):
            raise ValueError(
                "CellTypist overwrite_existing found a prior decision matrix whose shape does not match provenance."
            )
        unowned_obs_collisions = sorted(set(obs_collisions) - set(prior_owned_obs_keys))
        unowned_obsm_collisions = sorted(set(obsm_collisions) - set(prior_owned_obsm_keys))
        if unowned_obs_collisions or unowned_obsm_collisions:
            raise ValueError(
                "CellTypist overwrite_existing refuses to replace outputs not owned by the verified prior "
                f"artifact: obs={unowned_obs_collisions}, obsm={unowned_obsm_collisions}."
            )
        stale_obs_keys = [
            key for key in prior_owned_obs_keys if key in adata.obs and key not in obs_output_keys
        ]
        stale_obsm_keys = [
            key for key in prior_owned_obsm_keys if key in adata.obsm and key not in obsm_output_keys
        ]
        prior_artifact_verified = True

    if source_kind == "raw":
        if adata.raw is None:
            raise ValueError("CellTypist expression source 'raw' was selected, but adata.raw is unavailable.")
        matrix = adata.raw.X
        feature_names = adata.raw.var_names
        source_label = "raw"
    elif source_kind == "layer":
        if not isinstance(layer_name, str) or not layer_name.strip():
            raise ValueError("CellTypist layer expression source requires a non-empty layer_name.")
        layer_name = layer_name.strip()
        if layer_name not in adata.layers:
            raise ValueError(f"CellTypist expression layer not found: {layer_name!r}")
        matrix = adata.layers[layer_name]
        feature_names = adata.var_names
        source_label = f"layers[{layer_name!r}]"
    else:
        if layer_name not in {None, ""}:
            raise ValueError("CellTypist X expression source cannot carry a layer_name.")
        matrix = adata.X
        feature_names = adata.var_names
        source_label = "X"
    if tuple(matrix.shape) != (int(adata.n_obs), int(len(feature_names))):
        raise ValueError(
            f"CellTypist source {source_label} has shape {tuple(matrix.shape)}, expected "
            f"({int(adata.n_obs)}, {int(len(feature_names))})."
        )
    if not bool(feature_names.is_unique):
        raise ValueError(f"CellTypist source {source_label} requires unique feature identifiers.")
    query_features = []
    for feature in feature_names.tolist():
        if not isinstance(feature, str) or not feature.strip() or feature != feature.strip():
            raise ValueError("CellTypist requires non-empty gene-symbol string feature identifiers.")
        query_features.append(feature)
    if np.iscomplexobj(matrix) or not np.issubdtype(matrix.dtype, np.number):
        raise TypeError(f"CellTypist source {source_label} must be a real numeric matrix.")
    stored_values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    stored_values = np.asarray(stored_values)
    if stored_values.size and not bool(np.isfinite(stored_values).all()):
        raise ValueError(f"CellTypist source {source_label} contains non-finite values.")
    if stored_values.size and bool((stored_values < 0).any()):
        raise ValueError(f"CellTypist source {source_label} contains negative values.")
    cell_totals = np.asarray(matrix.sum(axis=1)).ravel()
    if bool((cell_totals <= 0).any()):
        raise ValueError(f"CellTypist source {source_label} contains cells with zero total expression.")
    transform = "none"
    if expression_state == "verified_counts":
        if stored_values.size and not bool(np.allclose(stored_values, np.rint(stored_values), rtol=0.0, atol=1e-8)):
            raise ValueError(
                f"CellTypist source {source_label} was declared verified_counts but contains non-integer values."
            )
        transform = "normalize_total_10000_then_log1p"
    else:
        max_value = float(stored_values.max()) if stored_values.size else 0.0
        if max_value > math.log1p(10_000.0) + 1e-5:
            raise ValueError(f"CellTypist source {source_label} exceeds the possible CP10K/log1p maximum.")
        if sparse.issparse(matrix):
            restored = matrix.astype(float).tocsr(copy=True)
            restored.data = np.expm1(restored.data)
            restored_totals = np.asarray(restored.sum(axis=1)).ravel()
        else:
            restored_totals = np.expm1(np.asarray(matrix, dtype=float)).sum(axis=1)
        if not bool(np.allclose(restored_totals, 10_000.0, rtol=0.0, atol=1.0)):
            observed = [float(restored_totals.min()), float(restored_totals.max())]
            maximum_absolute_deviation = float(np.max(np.abs(restored_totals - 10_000.0)))
            raise ValueError(
                f"CellTypist source {source_label} was declared verified_cp10k_log1p, but inverse-transformed "
                f"cell totals span {observed} with maximum absolute deviation "
                f"{maximum_absolute_deviation:.6g}; every cell must be within 1 count of 10,000."
            )

    over_values = None
    group_order = []
    group_counts = {}
    if majority_voting:
        if over_clustering_key not in adata.obs:
            raise ValueError(f"CellTypist over-clustering column not found in obs: {over_clustering_key!r}")
        groups = adata.obs[over_clustering_key]
        if not isinstance(groups.dtype, pd.CategoricalDtype):
            raise TypeError("CellTypist over-clustering membership must be categorical.")
        if bool(groups.isna().any()):
            raise ValueError("CellTypist over-clustering membership cannot contain missing values.")
        canonical_seen = {}
        group_order = []
        for category in groups.cat.categories.tolist():
            if not bool(pd.api.types.is_scalar(category)) or bool(pd.isna(category)):
                raise TypeError("CellTypist over-clustering categories must be non-missing scalars.")
            key = str(category).strip()
            if not key:
                raise ValueError("CellTypist over-clustering categories cannot be blank.")
            signature = (type(category).__module__, type(category).__qualname__, repr(category))
            if key in canonical_seen and canonical_seen[key] != signature:
                raise ValueError(f"CellTypist over-clustering categories collapse to canonical key {key!r}.")
            canonical_seen[key] = signature
            if int((groups == category).sum()) > 0:
                group_order.append(key)
        over_values = np.asarray([str(value).strip() for value in groups.tolist()], dtype=object)
        group_counts = {key: int(np.count_nonzero(over_values == key)) for key in group_order}

    existing_metadata = adata.uns.get("openbio_singlecell")
    if existing_metadata is not None:
        if not isinstance(existing_metadata, dict):
            raise ValueError("OpenBio metadata must be a mapping.")
        if existing_metadata.get("schema_version") != 1:
            raise ValueError(
                "OpenBio metadata has unsupported schema_version "
                f"{existing_metadata.get('schema_version')!r}; expected 1."
            )
    if prior_artifact_verified and prior_columns["selected_label"] != label_column:
        existing_annotations = (
            existing_metadata.get("annotations", {}) if existing_metadata is not None else {}
        )
        if not isinstance(existing_annotations, RuntimeMapping):
            raise ValueError("OpenBio annotation provenance must be a mapping.")
        prior_selected_label = prior_columns["selected_label"]
        if prior_selected_label in existing_annotations:
            if existing_annotations[prior_selected_label] != prior_provenance:
                raise ValueError(
                    "CellTypist overwrite_existing refuses to remove the prior annotation entry because its "
                    "value does not match the verified prior CellTypist provenance."
                )
            stale_annotation_keys.append(prior_selected_label)

    if celltypist_module is None:
        try:
            celltypist_module = importlib.import_module("celltypist")
        except (ImportError, OSError) as exc:
            raise RuntimeError(f"{operation} requires the celltypist package, but it is unavailable ({exc}).") from exc
    models_module = getattr(celltypist_module, "models", None)
    if models_module is None or not hasattr(models_module, "Model"):
        raise RuntimeError("CellTypist backend does not expose the official models.Model API.")
    resolved_model_path = _standalone_resolve_celltypist_model_path(model_identifier, models_module)
    digest = hashlib.sha256()
    with resolved_model_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    model_sha256 = digest.hexdigest()
    if expected_model_sha256 is not None and model_sha256 != expected_model_sha256:
        raise ValueError(
            f"CellTypist model fingerprint changed: expected {expected_model_sha256}, observed {model_sha256}."
        )
    try:
        loaded_model = models_module.Model.load(str(resolved_model_path))
    except Exception as exc:
        raise RuntimeError(
            f"CellTypist failed to load model {str(resolved_model_path)!r}. CellTypist .pkl files are Python "
            "pickle artifacts that can execute code; load only official or explicitly user-trusted sources."
        ) from exc

    def normalized_strings(values, description):
        normalized = []
        for value in list(values):
            if not isinstance(value, str) or not value.strip() or value != value.strip():
                raise RuntimeError(f"CellTypist model has invalid {description} entries.")
            normalized.append(value)
        if not normalized or len(normalized) != len(set(normalized)):
            raise RuntimeError(f"CellTypist model {description} must be non-empty and unique.")
        return normalized

    model_classes = normalized_strings(getattr(loaded_model, "cell_types", []), "class")
    if "Heterogeneous" in model_classes:
        raise RuntimeError("CellTypist model class vocabulary uses reserved label 'Heterogeneous'.")
    model_features = normalized_strings(getattr(loaded_model, "features", []), "feature")
    matched_features = sorted(set(query_features) & set(model_features))
    if not matched_features:
        raise ValueError(
            f"CellTypist source {source_label} has zero feature overlap with model {str(resolved_model_path)!r}."
        )
    overlap = {
        "matched_features": len(matched_features),
        "query_features": len(query_features),
        "model_features": len(model_features),
        "query_fraction": float(len(matched_features) / len(query_features)),
        "model_fraction": float(len(matched_features) / len(model_features)),
    }

    # CellTypist preprocessing may normalize/log the prediction input; keep a
    # dedicated algorithm workspace so the selected source remains unchanged.
    prediction_matrix = matrix.copy()
    prediction_input = sc.AnnData(
        X=prediction_matrix,
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=pd.Index(query_features)),
    )
    if expression_state == "verified_counts":
        sc.pp.normalize_total(prediction_input, target_sum=10_000.0, inplace=True)
        sc.pp.log1p(prediction_input)
    transformed_values = (
        prediction_input.X.data if sparse.issparse(prediction_input.X) else np.asarray(prediction_input.X).ravel()
    )
    if transformed_values.size and not bool(np.isfinite(transformed_values).all()):
        raise RuntimeError("CellTypist expression preparation produced non-finite values.")

    predictions = celltypist_module.annotate(
        prediction_input,
        model=loaded_model,
        transpose_input=False,
        gene_file=None,
        cell_file=None,
        mode="best match",
        p_thres=0.5,
        majority_voting=majority_voting,
        over_clustering=over_values,
        use_GPU=False,
        min_prop=min_prop,
    )
    predicted = getattr(predictions, "predicted_labels", None)
    probability = getattr(predictions, "probability_matrix", None)
    decision = getattr(predictions, "decision_matrix", None)
    expected_index = adata.obs_names
    for table, description in (
        (predicted, "predicted_labels"),
        (probability, "probability_matrix"),
        (decision, "decision_matrix"),
    ):
        if not isinstance(table, pd.DataFrame):
            raise RuntimeError(f"CellTypist returned malformed {description}; expected a pandas DataFrame.")
        if not table.index.equals(expected_index):
            raise RuntimeError(f"CellTypist returned misaligned {description} observation identifiers.")
    if list(probability.columns) != model_classes or list(decision.columns) != model_classes:
        raise RuntimeError("CellTypist returned class columns that do not match the loaded model order.")
    try:
        probability_values = probability.to_numpy(dtype=float)
        decision_values = decision.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("CellTypist returned nonnumeric probability or decision values.") from exc
    expected_shape = (int(adata.n_obs), len(model_classes))
    if probability_values.shape != expected_shape or decision_values.shape != expected_shape:
        raise RuntimeError(
            f"CellTypist returned invalid matrix shapes; expected {expected_shape}, observed "
            f"probability={probability_values.shape}, decision={decision_values.shape}."
        )
    if not bool(np.isfinite(probability_values).all()) or bool(
        ((probability_values < 0.0) | (probability_values > 1.0)).any()
    ):
        raise RuntimeError("CellTypist probabilities must be finite and within [0, 1].")
    if not bool(np.isfinite(decision_values).all()):
        raise RuntimeError("CellTypist decision scores must be finite.")
    if "predicted_labels" not in predicted:
        raise RuntimeError("CellTypist result is missing the individual predicted_labels column.")
    individual_labels = []
    for value in predicted["predicted_labels"].tolist():
        if not isinstance(value, str) or value not in model_classes:
            raise RuntimeError(f"CellTypist returned invalid individual label {value!r}.")
        individual_labels.append(value)
    class_positions = {label: index for index, label in enumerate(model_classes)}
    row_positions = np.arange(int(adata.n_obs))
    individual_positions = np.asarray([class_positions[label] for label in individual_labels])
    individual_probabilities = probability_values[row_positions, individual_positions]
    row_maximum = probability_values.max(axis=1)
    if not bool(np.allclose(individual_probabilities, row_maximum, rtol=1e-12, atol=1e-12)):
        raise RuntimeError("CellTypist individual labels are not aligned with row-maximum class probabilities.")

    majority_labels = None
    majority_support = None
    agreement = None
    majority_tied_groups = []
    if majority_voting:
        required_columns = {"over_clustering", "majority_voting"}
        if not required_columns.issubset(predicted.columns):
            raise RuntimeError("CellTypist backend did not complete the requested majority voting.")
        returned_groups = np.asarray([str(value).strip() for value in predicted["over_clustering"].tolist()])
        if not bool(np.array_equal(returned_groups, over_values)):
            raise RuntimeError("CellTypist returned over-clustering memberships that do not match the input.")
        majority_labels = []
        for value in predicted["majority_voting"].tolist():
            if not isinstance(value, str) or value not in {*model_classes, "Heterogeneous"}:
                raise RuntimeError(f"CellTypist returned invalid majority label {value!r}.")
            majority_labels.append(value)
        majority_support = np.empty(int(adata.n_obs), dtype=float)
        individual_array = np.asarray(individual_labels, dtype=object)
        majority_array = np.asarray(majority_labels, dtype=object)
        for group in group_order:
            mask = over_values == group
            group_labels = individual_array[mask]
            counts = {label: int(np.count_nonzero(group_labels == label)) for label in model_classes}
            maximum = max(counts.values())
            winners = [label for label, count in counts.items() if count == maximum]
            if len(winners) > 1:
                majority_tied_groups.append(group)
            support = float(maximum / int(mask.sum()))
            returned = set(majority_array[mask].tolist())
            if len(returned) != 1:
                raise RuntimeError(f"CellTypist returned inconsistent majority labels within group {group!r}.")
            returned_label = next(iter(returned))
            if returned_label == "Heterogeneous":
                if not support < min_prop:
                    raise RuntimeError(
                        f"CellTypist returned Heterogeneous for group {group!r} despite support {support}."
                    )
            elif support < min_prop or counts[returned_label] != maximum:
                raise RuntimeError(
                    f"CellTypist returned invalid majority label {returned_label!r} for group {group!r}."
                )
            majority_support[mask] = support
        selected_labels = majority_labels
        selected_probabilities = np.asarray(
            [
                probability_values[index, class_positions[label]] if label in class_positions else np.nan
                for index, label in enumerate(majority_labels)
            ],
            dtype=float,
        )
        agreement = float(
            np.mean(np.asarray(individual_labels, dtype=object) == np.asarray(majority_labels, dtype=object))
        )
    else:
        selected_labels = individual_labels
        selected_probabilities = individual_probabilities.copy()
    selected_array = np.asarray(selected_probabilities, dtype=float)
    selected_missing = np.isnan(selected_array)
    expected_missing = np.asarray([label == "Heterogeneous" for label in selected_labels], dtype=bool)
    if not bool(np.array_equal(selected_missing, expected_missing)):
        raise RuntimeError("CellTypist selected confidence is missing outside the Heterogeneous label semantics.")
    selected_finite = selected_array[~selected_missing]
    if not bool(np.isfinite(selected_finite).all()) or bool(((selected_finite < 0.0) | (selected_finite > 1.0)).any()):
        raise RuntimeError("CellTypist selected confidence must be finite in [0, 1] when a model class is selected.")
    if majority_support is not None and (
        not bool(np.isfinite(majority_support).all())
        or bool(((majority_support <= 0.0) | (majority_support > 1.0)).any())
    ):
        raise RuntimeError("CellTypist majority support must be finite in (0, 1].")

    def plain_json(value):
        if value is None or isinstance(value, (str, bool, int)):
            return value
        if isinstance(value, float):
            return value if math.isfinite(value) else None
        if hasattr(value, "item"):
            try:
                return plain_json(value.item())
            except (TypeError, ValueError):
                pass
        if isinstance(value, RuntimeMapping):
            return {str(key): plain_json(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
            sequence = value.tolist() if hasattr(value, "tolist") else value
            return [plain_json(item) for item in sequence]
        return str(value)

    description = getattr(loaded_model, "description", {})
    model_metadata = plain_json(description if isinstance(description, RuntimeMapping) else {})
    model_provenance = {
        "requested_identifier": requested_model_identifier,
        "resolved_path": str(resolved_model_path),
        "sha256": model_sha256,
        "description": model_metadata,
        "class_count": len(model_classes),
        "feature_count": len(model_features),
    }
    selected_category_order = model_classes + (["Heterogeneous"] if "Heterogeneous" in selected_labels else [])
    selected_counts = {
        label: int(sum(value == label for value in selected_labels))
        for label in selected_category_order
    }
    individual_counts = {
        label: int(sum(value == label for value in individual_labels))
        for label in model_classes
    }
    warnings = [
        "CellTypist predictions are Provisional annotation and require marker, tissue, and study-context review before curation.",
        "CellTypist class scores are independent sigmoid outputs, are not multiclass probabilities summing to one, and are not calibrated biological certainty.",
        "The closed model vocabulary can force novel, low-quality, doublet, or out-of-domain cells to a known class.",
        "CellTypist .pkl models are Python pickle artifacts that can execute code when loaded; use only official or explicitly user-trusted sources. The SHA-256 identifies exact bytes but does not prove safety.",
    ]
    if len(matched_features) < len(model_features):
        warnings.append(
            "The query does not contain every model feature; incomplete feature overlap can weaken or shift predictions."
        )
    if majority_voting:
        warnings.append(
            "Majority voting smooths labels within caller-supplied groups and can hide rare, transitional, or mixed states."
        )
    if majority_tied_groups:
        warnings.append(
            f"{len(majority_tied_groups):,} over-clustering groups tied for the dominant individual label; "
            "the backend-selected tied class is retained and the ambiguity is disclosed."
        )
    parameters = {"source": source_kind}
    if source_kind == "layer":
        parameters["layer_name"] = layer_name
    parameters.update(
        {
            "expression_state": expression_state,
            "expression_transform": transform,
            "model": requested_model_identifier,
            "resolved_model_path": str(resolved_model_path),
            "model_sha256": model_sha256,
            "majority_voting": majority_voting,
            "over_clustering_key": over_clustering_key if majority_voting else "",
            "min_prop": min_prop,
            "label_column": label_column,
            "confidence_column": confidence_column,
            "individual_label_column": individual_label_column,
            "individual_probability_column": individual_probability_column,
            "majority_label_column": majority_label_column if majority_voting else None,
            "majority_support_column": majority_support_column if majority_voting else None,
            "probability_key": probability_key,
            "store_decision_matrix": store_decision_matrix,
            "decision_key": decision_key if store_decision_matrix else None,
            "metadata_key": metadata_key,
            "overwrite_existing": overwrite_existing,
            "target_sum": 10_000.0 if expression_state == "verified_counts" else None,
            "mode": "best match",
            "p_thres": 0.5,
            "transpose_input": False,
            "use_GPU": False,
        }
    )
    software_versions = _standalone_celltypist_software_versions(
        openbio_version=openbio_version,
    )
    provenance = {
        "schema_version": 1,
        "operation": "celltypist_annotation",
        "annotation_status": "provisional",
        "expression": {
            "source": source_kind,
            "layer_name": layer_name if source_kind == "layer" else None,
            "declared_state": expression_state,
            "transform": transform,
            "target_sum": 10_000.0,
            "query_feature_count": len(query_features),
        },
        "model": model_provenance,
        "feature_overlap": overlap,
        "classes": model_classes,
        "fixed_policy": {
            "mode": "best match",
            "p_thres": 0.5,
            "transpose_input": False,
            "gene_file": None,
            "cell_file": None,
            "use_GPU": False,
        },
        "majority_voting": majority_voting,
        "over_clustering_key": over_clustering_key if majority_voting else None,
        "min_prop": min_prop,
        "overwrite_existing": overwrite_existing,
        "columns": {
            "selected_label": label_column,
            "selected_label_probability": confidence_column,
            "individual_label": individual_label_column,
            "individual_row_max_probability": individual_probability_column,
            "majority_label": majority_label_column if majority_voting else None,
            "majority_support": majority_support_column if majority_voting else None,
        },
        "matrices": {
            "probability_key": probability_key,
            "probability_shape": list(probability_values.shape),
            "decision_key": decision_key if store_decision_matrix else None,
            "decision_shape": list(decision_values.shape) if store_decision_matrix else None,
        },
        "artifact_replacement": {
            "prior_artifact_verified": prior_artifact_verified,
            "overwritten_obs_keys": obs_collisions,
            "overwritten_obsm_keys": obsm_collisions,
            "removed_stale_obs_keys": stale_obs_keys,
            "removed_stale_obsm_keys": stale_obsm_keys,
            "removed_stale_annotation_keys": stale_annotation_keys,
        },
        "software_versions": software_versions,
        "warnings": warnings,
        "observation_axis_fingerprint_sha256": _celltypist_axis_fingerprint(adata.obs_names),
        "probability_content_fingerprint_sha256": _celltypist_probability_fingerprint(
            probability_values,
            observation_ids=adata.obs_names,
            classes=model_classes,
        ),
        "selected_annotation_fingerprint_sha256": _celltypist_selected_fingerprint(
            selected_labels,
            selected_probabilities,
            observation_ids=adata.obs_names,
        ),
    }
    diagnostics = {
        "cells": int(adata.n_obs),
        "source_features": len(query_features),
        "expression_source": source_label,
        "expression_state": expression_state,
        "expression_transform": transform,
        "model": model_provenance,
        "feature_overlap": overlap,
        "class_order": model_classes,
        "individual_label_counts": individual_counts,
        "selected_label_counts": selected_counts,
        "heterogeneous_cells": int(sum(value == "Heterogeneous" for value in selected_labels)),
        "selected_probabilities": selected_probabilities.tolist(),
        "individual_probabilities": individual_probabilities.tolist(),
        "probability_shape": list(probability_values.shape),
        "decision_shape": list(decision_values.shape),
        "majority_voting_completed": majority_voting,
        "over_clustering_key": over_clustering_key if majority_voting else None,
        "group_count": len(group_order),
        "group_order_preview": group_order[:100],
        "group_counts_preview": {group: group_counts[group] for group in group_order[:100]},
        "group_preview_truncated": len(group_order) > 100,
        "majority_tied_group_count": len(majority_tied_groups),
        "majority_tied_groups_preview": majority_tied_groups[:100],
        "majority_tied_groups_preview_truncated": len(majority_tied_groups) > 100,
        "majority_support": majority_support.tolist() if majority_support is not None else [],
        "individual_majority_agreement": agreement,
        "artifact_replacement": provenance["artifact_replacement"],
        "software_versions": software_versions,
        "warnings": warnings,
        "parameters": parameters,
        "provenance": provenance,
    }

    output = adata
    for key in stale_obs_keys:
        del output.obs[key]
    for key in stale_obsm_keys:
        del output.obsm[key]
    output.obs[individual_label_column] = pd.Categorical(
        individual_labels,
        categories=model_classes,
        ordered=False,
    )
    output.obs[individual_probability_column] = individual_probabilities
    output.obs[label_column] = pd.Categorical(
        selected_labels,
        categories=selected_category_order,
        ordered=False,
    )
    output.obs[confidence_column] = selected_probabilities
    if majority_voting:
        output.obs[majority_label_column] = pd.Categorical(
            majority_labels,
            categories=selected_category_order,
            ordered=False,
        )
        output.obs[majority_support_column] = majority_support
    output.obsm[probability_key] = probability_values.copy()
    if store_decision_matrix:
        output.obsm[decision_key] = decision_values.copy()
    output.uns[metadata_key] = provenance
    metadata = (
        copy.deepcopy(existing_metadata)
        if existing_metadata is not None
        else {
            "schema_version": 1,
            "version": "not-installed",
            "display_name": "AnnData",
            "source": {},
            "random_seed": 0,
            "warnings": [],
            "analysis_history": {},
        }
    )
    if existing_metadata is None:
        try:
            metadata["version"] = importlib_metadata.version("openbio-singlecell")
        except importlib_metadata.PackageNotFoundError:
            pass
    annotations = metadata.get("annotations", {})
    if not isinstance(annotations, dict):
        raise ValueError("OpenBio annotation provenance must be a mapping.")
    annotations = copy.deepcopy(annotations)
    for key in stale_annotation_keys:
        del annotations[key]
    annotations[label_column] = provenance
    metadata["annotations"] = annotations
    output.uns["openbio_singlecell"] = metadata
    if not output.obs_names.equals(adata.obs_names) or not output.var_names.equals(adata.var_names):
        raise RuntimeError(f"{operation} changed observation or feature identity.")
    if output.obsm[probability_key].shape != expected_shape:
        raise RuntimeError(f"{operation} stored a probability matrix with the wrong shape.")
    if output.obs[label_column].cat.categories.tolist() != selected_category_order:
        raise RuntimeError(f"{operation} stored the wrong selected-label category order.")
    if output.obs[individual_label_column].cat.categories.tolist() != model_classes:
        raise RuntimeError(f"{operation} stored the wrong individual-label category order.")
    if output.uns.get(metadata_key) != provenance:
        raise RuntimeError(f"{operation} failed to store the pinned annotation provenance.")
    if _return_diagnostics:
        return output, diagnostics
    return output


def run_map_cluster_annotations(
    adata: Any,
    *,
    groupby: str,
    mapping_json: str,
    output_column: str,
    unmapped_policy: str,
    annotation_status: str,
    overwrite_existing: bool,
) -> tuple[Any, dict[str, Any]]:
    return _standalone_map_cluster_annotations(
        adata,
        groupby=groupby,
        mapping_json=mapping_json,
        output_column=output_column,
        unmapped_policy=unmapped_policy,
        annotation_status=annotation_status,
        overwrite_existing=overwrite_existing,
        _return_diagnostics=True,
    )


def run_celltypist_annotation(
    adata: Any,
    *,
    source_kind: str,
    layer_name: str | None,
    expression_state: str,
    model_identifier: str,
    majority_voting: bool,
    over_clustering_key: str,
    min_prop: float,
    label_column: str,
    confidence_column: str,
    probability_key: str,
    store_decision_matrix: bool,
    decision_key: str,
    metadata_key: str,
    overwrite_existing: bool,
    openbio_version: str,
    celltypist_module: Any,
) -> tuple[Any, dict[str, Any]]:
    return _standalone_celltypist_annotation(
        adata,
        source_kind=source_kind,
        layer_name=layer_name,
        expression_state=expression_state,
        model_identifier=model_identifier,
        majority_voting=majority_voting,
        over_clustering_key=over_clustering_key,
        min_prop=min_prop,
        label_column=label_column,
        confidence_column=confidence_column,
        probability_key=probability_key,
        store_decision_matrix=store_decision_matrix,
        decision_key=decision_key,
        metadata_key=metadata_key,
        overwrite_existing=overwrite_existing,
        openbio_version=openbio_version,
        celltypist_module=celltypist_module,
        _return_diagnostics=True,
    )


def celltypist_model_fingerprint(model_identifier: str) -> tuple[Any, ...]:
    """Fingerprint an explicit or locally installed model without network access."""
    import hashlib
    import importlib
    import os
    from pathlib import Path

    if not isinstance(model_identifier, str) or not model_identifier.strip():
        return ("openbio-celltypist-model-v1", "invalid", "empty-model")
    requested = model_identifier.strip()
    candidate = Path(requested).expanduser()
    if candidate.is_file():
        resolved = candidate.resolve()
    else:
        try:
            celltypist_module = importlib.import_module("celltypist")
        except (ImportError, OSError):
            return ("openbio-celltypist-model-v1", "backend-unavailable", requested)
        models_module = getattr(celltypist_module, "models", None)
        if models_module is None or not hasattr(models_module, "Model"):
            return ("openbio-celltypist-model-v1", "backend-invalid", requested)
        try:
            resolved = _standalone_resolve_celltypist_model_path(requested, models_module)
        except FileNotFoundError:
            return ("openbio-celltypist-model-v1", "missing", requested)
    stat = resolved.stat()
    digest = hashlib.sha256()
    with resolved.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return (
        "openbio-celltypist-model-v1",
        os.path.normcase(os.path.realpath(resolved)),
        int(stat.st_size),
        int(stat.st_mtime_ns),
        digest.hexdigest(),
    )


def build_celltypist_summary(diagnostics: dict[str, Any]) -> dict[str, Any]:
    return _standalone_celltypist_summary(diagnostics)


def _standalone_source(function: Any) -> str:
    return textwrap.dedent(inspect.getsource(function)).strip()


def map_cluster_annotations_code(
    *,
    groupby: str,
    mapping_json: str,
    output_column: str,
    unmapped_policy: str,
    annotation_status: str,
    overwrite_existing: bool,
) -> str:
    return f"""from __future__ import annotations

{_standalone_source(_standalone_map_cluster_annotations)}


def map_cluster_annotations(adata):
    return _standalone_map_cluster_annotations(
        adata,
        groupby={groupby!r},
        mapping_json={mapping_json!r},
        output_column={output_column!r},
        unmapped_policy={unmapped_policy!r},
        annotation_status={annotation_status!r},
        overwrite_existing={overwrite_existing!r},
    )
"""


def celltypist_annotation_code(
    *,
    source_kind: str,
    layer_name: str | None,
    expression_state: str,
    resolved_model_path: str,
    requested_model_identifier: str,
    expected_model_sha256: str,
    majority_voting: bool,
    over_clustering_key: str,
    min_prop: float,
    label_column: str,
    confidence_column: str,
    probability_key: str,
    store_decision_matrix: bool,
    decision_key: str,
    metadata_key: str,
    overwrite_existing: bool,
    openbio_version: str,
) -> str:
    return f"""from __future__ import annotations

{_standalone_source(_standalone_resolve_celltypist_model_path)}


{_standalone_source(_celltypist_axis_fingerprint)}


{_standalone_source(_celltypist_probability_fingerprint)}


{_standalone_source(_celltypist_selected_fingerprint)}


{_standalone_source(_standalone_celltypist_annotation)}


{_standalone_source(_standalone_celltypist_software_versions)}


{_standalone_source(_standalone_celltypist_summary)}


def run_celltypist_annotation(adata):
    output, diagnostics = _standalone_celltypist_annotation(
        adata,
        source_kind={source_kind!r},
        layer_name={layer_name!r},
        expression_state={expression_state!r},
        model_identifier={resolved_model_path!r},
        requested_model_identifier={requested_model_identifier!r},
        expected_model_sha256={expected_model_sha256!r},
        majority_voting={majority_voting!r},
        over_clustering_key={over_clustering_key!r},
        min_prop={min_prop!r},
        label_column={label_column!r},
        confidence_column={confidence_column!r},
        probability_key={probability_key!r},
        store_decision_matrix={store_decision_matrix!r},
        decision_key={decision_key!r},
        metadata_key={metadata_key!r},
        overwrite_existing={overwrite_existing!r},
        openbio_version={openbio_version!r},
        _return_diagnostics=True,
    )
    summary = _standalone_celltypist_summary(diagnostics)
    return output, summary
"""


__all__ = [
    "_celltypist_axis_fingerprint",
    "_celltypist_probability_fingerprint",
    "_celltypist_selected_fingerprint",
    "build_celltypist_summary",
    "celltypist_annotation_code",
    "celltypist_model_fingerprint",
    "map_cluster_annotations_code",
    "run_celltypist_annotation",
    "run_map_cluster_annotations",
]
