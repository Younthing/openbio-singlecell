from __future__ import annotations

import inspect
import textwrap
from typing import TYPE_CHECKING, Any

from .dgidb_resource import (
    DGIdbResource,
    _standalone_dgidb_json_sha256,
    _standalone_validate_portable_dgidb_resource,
    dgidb_portable_validation_code,
    validate_dgidb_resource,
)

if TYPE_CHECKING:
    from anndata import AnnData


def _standalone_run_drug_score(
    adata,
    resource_table,
    resource_metadata,
    resource_accounting,
    resource_artifact_metadata,
    *,
    expected_resource_artifact_fingerprint,
    drug,
    source_kind="X",
    layer_name=None,
    output_key="drug_target_score",
    min_matched_targets=1,
    overwrite_existing=False,
    openbio_version="not-installed",
    pertpy_module=None,
):
    """Score one pinned DGIdb target set with Pertpy 1.3 mean scoring and independent verification."""
    import hashlib
    import importlib
    import inspect as runtime_inspect
    import json
    import platform
    from collections.abc import Mapping
    from importlib import metadata as importlib_metadata

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import sparse

    operation = "DGIdb Drug Score"
    resource_table, resource_metadata, resource_accounting, resource_artifact_metadata = (
        _standalone_validate_portable_dgidb_resource(
            resource_table,
            resource_metadata,
            resource_accounting,
            resource_artifact_metadata,
            copy_table=False,
        )
    )
    actual_resource_fingerprint = resource_artifact_metadata["artifact_fingerprint_sha256"]
    if actual_resource_fingerprint != expected_resource_artifact_fingerprint:
        raise ValueError(f"{operation} resource artifact differs from the pinned upstream resource.")
    if not isinstance(drug, str) or not drug or drug != drug.strip():
        raise ValueError(f"{operation} drug must be an exact nonblank whitespace-canonical identifier.")
    drug_rows = resource_table.loc[resource_table["drug"] == drug].copy()
    if drug_rows.empty:
        raise ValueError(f"{operation} drug {drug!r} is absent from the pinned DGIdb resource.")
    requested_targets = drug_rows["gene"].tolist()
    if any("|" in gene for gene in requested_targets):
        raise ValueError(f"{operation} HGNC target identifiers cannot contain Pertpy's '|' delimiter.")
    if isinstance(min_matched_targets, bool) or not isinstance(min_matched_targets, int) or min_matched_targets < 1:
        raise ValueError(f"{operation} min_matched_targets must be a positive integer.")
    if not isinstance(output_key, str) or not output_key or output_key != output_key.strip():
        raise ValueError(f"{operation} output_key must be nonblank whitespace-canonical text.")
    if not isinstance(overwrite_existing, bool):
        raise TypeError(f"{operation} overwrite_existing must be boolean.")
    if not hasattr(adata, "X") or not hasattr(adata, "obs_names") or not hasattr(adata, "var_names"):
        raise TypeError(f"{operation} requires an AnnData-like input.")
    if int(adata.n_obs) < 1:
        raise ValueError(f"{operation} requires nonempty observations.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique observation identifiers.")
    observation_names = adata.obs_names.tolist()
    active_feature_names = adata.var_names.tolist()
    if any(not isinstance(value, str) or not value or value != value.strip() for value in observation_names):
        raise ValueError(f"{operation} observation identifiers must be whitespace-canonical nonblank strings.")
    if source_kind not in {"X", "layer", "raw"}:
        raise ValueError(f"{operation} source_kind must be one of 'X', 'layer', or 'raw'.")
    if source_kind == "layer":
        if not bool(adata.var_names.is_unique):
            raise ValueError(f"{operation} requires unique selected-layer feature identifiers.")
        if any(not isinstance(value, str) or not value or value != value.strip() for value in active_feature_names):
            raise ValueError(
                f"{operation} selected-layer feature identifiers must be whitespace-canonical nonblank strings."
            )
        if not isinstance(layer_name, str) or not layer_name or layer_name != layer_name.strip():
            raise ValueError(f"{operation} layer_name must be nonblank whitespace-canonical text.")
        if layer_name not in adata.layers:
            raise ValueError(f"{operation} expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
        feature_names = active_feature_names
    elif source_kind == "raw":
        if layer_name is not None:
            raise ValueError(f"{operation} Raw source cannot carry layer_name.")
        if adata.raw is None:
            raise ValueError(f"{operation} Raw source was requested but adata.raw is absent.")
        if adata.raw.obs_names.tolist() != observation_names:
            raise ValueError(f"{operation} Raw observation axis does not align exactly to adata.obs_names.")
        feature_names = adata.raw.var_names.tolist()
        if not bool(adata.raw.var_names.is_unique):
            raise ValueError(f"{operation} requires unique Raw feature identifiers.")
        if any(not isinstance(value, str) or not value or value != value.strip() for value in feature_names):
            raise ValueError(f"{operation} Raw feature identifiers must be whitespace-canonical nonblank strings.")
        matrix = adata.raw.X
    else:
        if not bool(adata.var_names.is_unique):
            raise ValueError(f"{operation} requires unique X feature identifiers.")
        if any(not isinstance(value, str) or not value or value != value.strip() for value in active_feature_names):
            raise ValueError(f"{operation} X feature identifiers must be whitespace-canonical nonblank strings.")
        if layer_name is not None:
            raise ValueError(f"{operation} X source cannot carry layer_name.")
        matrix = adata.X
        feature_names = active_feature_names
    if not feature_names:
        raise ValueError(f"{operation} selected expression source requires at least one feature.")
    if tuple(getattr(matrix, "shape", ())) != (len(observation_names), len(feature_names)):
        raise ValueError(f"{operation} expression matrix does not align to its named axes.")
    if not (sparse.issparse(matrix) or isinstance(matrix, np.ndarray)):
        raise TypeError(f"{operation} requires an in-memory NumPy or SciPy sparse expression matrix.")
    stored_values = np.asarray(matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel())
    if stored_values.dtype.kind not in "iuf":
        raise TypeError(f"{operation} expression must be numeric.")
    if stored_values.size and not bool(np.isfinite(stored_values).all()):
        raise ValueError(f"{operation} expression contains non-finite values.")
    if sparse.issparse(matrix):
        canonical_matrix = sparse.csr_matrix(matrix, copy=True)
        canonical_matrix.sum_duplicates()
        canonical_matrix.eliminate_zeros()
        canonical_matrix.sort_indices()
        nonzero_per_observation = np.asarray(canonical_matrix.getnnz(axis=1)).ravel()
        semantic_values = canonical_matrix.data
    else:
        canonical_matrix = np.asarray(matrix)
        nonzero_per_observation = np.count_nonzero(canonical_matrix, axis=1)
        semantic_values = canonical_matrix.ravel()
    zero_observation_count = int((nonzero_per_observation == 0).sum())
    value_min = float(np.min(semantic_values)) if semantic_values.size else 0.0
    integer_like = bool(
        not semantic_values.size
        or np.allclose(semantic_values, np.rint(semantic_values), rtol=0.0, atol=1e-8)
    )

    def expression_state():
        x_state = "unknown"
        x_evidence = None
        layer_states = {}
        metadata_root = adata.uns.get("openbio_singlecell")
        history = metadata_root.get("analysis_history") if isinstance(metadata_root, Mapping) else None
        entries = history.values() if isinstance(history, Mapping) else ()
        for entry in entries:
            if not isinstance(entry, Mapping):
                continue
            history_operation = entry.get("operation")
            parameters = entry.get("parameters")
            parameters = parameters if isinstance(parameters, Mapping) else {}
            if history_operation == "snapshot_expression":
                layer_states["counts"] = ("counts", "OpenBio Snapshot Expression history")
                if parameters.get("source") == "X":
                    x_state, x_evidence = "counts", "OpenBio Snapshot Expression history"
            elif history_operation == "normalize_to_layer":
                output_layer = parameters.get("output_layer")
                transform = parameters.get("transform")
                if isinstance(output_layer, str):
                    state = {"log1p": "logged", "sqrt": "transformed", "none": "normalized"}.get(
                        transform, "normalized"
                    )
                    layer_states[output_layer] = (state, "OpenBio Normalize To Layer history")
            elif history_operation == "pearson_residuals_to_layer":
                output_layer = parameters.get("output_layer")
                if isinstance(output_layer, str):
                    layer_states[output_layer] = (
                        "pearson_residuals",
                        "OpenBio Pearson Residuals history",
                    )
            elif history_operation == "scale_to_layer":
                output_layer = parameters.get("output_layer")
                if isinstance(output_layer, str):
                    layer_states[output_layer] = ("scaled", "OpenBio Scale history")
            if history_operation == "normalize_total":
                x_state, x_evidence = "normalized", "OpenBio Normalize Total history"
            elif history_operation == "log1p":
                x_state, x_evidence = "logged", "OpenBio Log1p history"
        if source_kind == "raw":
            return "counts", "explicit AnnData Raw snapshot selection"
        if source_kind == "layer":
            return layer_states.get(layer_name, ("unknown", None))
        if x_state == "unknown" and isinstance(adata.uns.get("log1p"), Mapping):
            return "logged", "AnnData uns['log1p'] marker"
        return x_state, x_evidence

    resolved_state, state_evidence = expression_state()
    warnings = []
    if resolved_state in {"scaled", "pearson_residuals"}:
        warnings.append(
            f"The explicitly selected source is described as {resolved_state!r} ({state_evidence}); mean target "
            "values remain calculable but describe that transformed scale rather than normalized abundance."
        )
    elif resolved_state == "counts":
        warnings.append(
            "Explicitly selected count-like expression is being scored. Mean target expression is numerically "
            "defined, but library size and sequencing depth can dominate the descriptive score; normalized "
            "continuous expression is recommended for cross-cell interpretation."
        )
    elif resolved_state == "unknown" and integer_like and value_min >= 0.0:
        warnings.append(
            "Expression provenance is unknown and values are nonnegative count-like. Mean target expression is "
            "numerically defined, but library size and sequencing depth can dominate the descriptive score."
        )
    elif resolved_state == "unknown":
        warnings.append(
            "Expression provenance is unknown but values are non-count-like; scoring assumes caller-supplied normalized continuous expression."
        )
    elif resolved_state not in {"normalized", "logged", "transformed"}:
        warnings.append(
            f"Expression state {resolved_state!r} is nonstandard for this score; the explicitly selected finite "
            "values were averaged unchanged."
        )

    feature_set = set(feature_names)
    matched_targets = [gene for gene in requested_targets if gene in feature_set]
    unmatched_targets = [gene for gene in requested_targets if gene not in feature_set]
    if len(matched_targets) < min_matched_targets:
        raise ValueError(
            f"{operation} matched {len(matched_targets):,} targets, below min_matched_targets={min_matched_targets:,}."
        )
    if output_key in adata.obs and not overwrite_existing:
        raise ValueError(
            f"{operation} output obs[{output_key!r}] already exists; enable overwrite_existing explicitly."
        )

    def expression_fingerprint(target_matrix):
        digest = hashlib.sha256()
        digest.update(b"openbio-singlecell/drug-score-expression/v1\0")
        digest.update(json.dumps(observation_names, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
        digest.update(json.dumps(feature_names, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
        digest.update(source_kind.encode("utf-8"))
        digest.update(str(layer_name).encode("utf-8"))
        if sparse.issparse(target_matrix):
            csr = sparse.csr_matrix(target_matrix, copy=True)
            csr.sum_duplicates()
            csr.eliminate_zeros()
            csr.sort_indices()
            digest.update(b"csr\0")
            digest.update(np.ascontiguousarray(csr.indptr, dtype="<i8").tobytes())
            digest.update(np.ascontiguousarray(csr.indices, dtype="<i8").tobytes())
            digest.update(np.ascontiguousarray(csr.data, dtype="<f8").tobytes())
        else:
            digest.update(b"dense\0")
            dense = np.asarray(target_matrix)
            for start in range(0, dense.shape[0], 1024):
                digest.update(np.ascontiguousarray(dense[start : start + 1024], dtype="<f8").tobytes())
        return digest.hexdigest()

    input_expression_fingerprint = expression_fingerprint(matrix)
    if pertpy_module is None:
        try:
            pertpy_module = importlib.import_module("pertpy")
        except (ImportError, OSError) as exc:
            raise RuntimeError(f"{operation} requires Pertpy 1.3, but it is unavailable ({exc}).") from exc
    version = getattr(pertpy_module, "__version__", None)
    if not isinstance(version, str) or not version.strip():
        try:
            version = importlib_metadata.version("pertpy")
        except importlib_metadata.PackageNotFoundError as exc:
            raise RuntimeError(f"{operation} cannot determine the installed Pertpy version.") from exc
    try:
        version_parts = version.split(".")
        version_major = int(version_parts[0])
        version_minor = int(version_parts[1])
    except (IndexError, TypeError, ValueError) as exc:
        raise RuntimeError(f"{operation} received invalid Pertpy version {version!r}.") from exc
    if version_major != 1 or version_minor != 3:
        raise RuntimeError(f"{operation} requires Pertpy 1.3.x; installed version is {version!r}.")
    enrichment_class = getattr(getattr(pertpy_module, "tl", None), "Enrichment", None)
    if not callable(enrichment_class):
        raise RuntimeError(f"{operation} requires the public pertpy.tl.Enrichment interface.")
    enrichment = enrichment_class()
    score_method = getattr(enrichment, "score", None)
    if not callable(score_method):
        raise RuntimeError(f"{operation} requires pertpy.tl.Enrichment.score.")
    try:
        signature = runtime_inspect.signature(score_method)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{operation} could not inspect Pertpy Enrichment.score.") from exc
    required_parameters = {"adata", "layer", "targets", "nested", "method", "key_added"}
    if not required_parameters.issubset(signature.parameters):
        missing = sorted(required_parameters - set(signature.parameters))
        raise RuntimeError(f"{operation} Pertpy score is missing required parameters: {missing}.")
    if signature.parameters["adata"].kind not in {
        runtime_inspect.Parameter.POSITIONAL_ONLY,
        runtime_inspect.Parameter.POSITIONAL_OR_KEYWORD,
    }:
        raise RuntimeError(f"{operation} Pertpy score no longer accepts adata as one public positional input.")
    keyword_parameters = required_parameters - {"adata"}
    incompatible = sorted(
        name
        for name in keyword_parameters
        if signature.parameters[name].kind
        not in {
            runtime_inspect.Parameter.POSITIONAL_OR_KEYWORD,
            runtime_inspect.Parameter.KEYWORD_ONLY,
        }
    )
    if incompatible:
        raise RuntimeError(
            f"{operation} Pertpy score parameters are no longer keyword-compatible: {incompatible}."
        )

    # Scientific workspace: Pertpy writes several scratch slots and must not replace the selected source matrix.
    work = ad.AnnData(
        X=matrix.copy(),
        obs=pd.DataFrame(index=pd.Index(observation_names)),
        var=pd.DataFrame(index=pd.Index(feature_names)),
    )
    work_matrix_before = expression_fingerprint(work.X)
    scratch_key = "_openbio_dgidb_single_drug"
    score_method(
        work,
        layer=None,
        targets={drug: list(requested_targets)},
        nested=False,
        method="mean",
        key_added=scratch_key,
    )
    work_matrix_after = expression_fingerprint(work.X)
    if work_matrix_after != work_matrix_before:
        raise RuntimeError(f"{operation} Pertpy backend mutated the selected expression matrix.")
    if expression_fingerprint(matrix) != input_expression_fingerprint:
        raise RuntimeError(f"{operation} backend mutated the caller's input expression.")
    resource_table_after, _, _, resource_artifact_after = _standalone_validate_portable_dgidb_resource(
        resource_table,
        resource_metadata,
        resource_accounting,
        resource_artifact_metadata,
        copy_table=False,
    )
    if (
        resource_artifact_after["artifact_fingerprint_sha256"] != actual_resource_fingerprint
        or not resource_table_after.equals(resource_table)
    ):
        raise RuntimeError(f"{operation} resource changed during backend execution.")

    score_key = f"{scratch_key}_score"
    variables_key = f"{scratch_key}_variables"
    genes_key = f"{scratch_key}_genes"
    all_genes_key = f"{scratch_key}_all_genes"
    missing_keys = [key for key in (score_key, variables_key, genes_key, all_genes_key) if key not in work.uns]
    if missing_keys:
        raise RuntimeError(f"{operation} Pertpy output is missing scratch fields: {missing_keys}.")
    scores = np.asarray(work.uns[score_key], dtype=float)
    if scores.shape != (len(observation_names), 1) or not bool(np.isfinite(scores).all()):
        raise RuntimeError(f"{operation} Pertpy score must be one finite column aligned to observations.")
    variables = list(work.uns[variables_key])
    if variables != [drug]:
        raise RuntimeError(f"{operation} Pertpy variables differ from the requested one-drug target set.")
    try:
        backend_matched = work.uns[genes_key]["var"].loc[drug, "genes"].split("|")
        backend_all = work.uns[all_genes_key]["var"].loc[drug, "all_genes"].split("|")
    except (AttributeError, KeyError, TypeError) as exc:
        raise RuntimeError(f"{operation} Pertpy gene-accounting output is malformed.") from exc
    backend_matched = [] if backend_matched == [""] else backend_matched
    backend_all = [] if backend_all == [""] else backend_all
    matched_feature_order = [gene for gene in feature_names if gene in set(requested_targets)]
    if backend_matched != matched_feature_order or backend_all != requested_targets:
        raise RuntimeError(f"{operation} Pertpy target accounting differs from the exact pinned target set.")
    target_indexes = [feature_names.index(gene) for gene in matched_targets]
    expected_scores = (
        np.asarray(matrix[:, target_indexes].mean(axis=1)).ravel()
        if sparse.issparse(matrix)
        else np.asarray(matrix)[:, target_indexes].mean(axis=1)
    )
    observed_scores = scores[:, 0]
    if not bool(np.allclose(observed_scores, expected_scores, rtol=1e-12, atol=1e-12)):
        raise RuntimeError(f"{operation} Pertpy scores disagree with the independently computed target mean.")
    output = adata
    output.obs[output_key] = pd.Series(observed_scores, index=output.obs_names, dtype=float)

    source_values = []
    evidence_values = []
    for encoded in drug_rows["sources"].tolist():
        for item in json.loads(encoded):
            if item not in source_values:
                source_values.append(item)
    for encoded in drug_rows["evidence"].tolist():
        for item in json.loads(encoded):
            if item not in evidence_values:
                evidence_values.append(item)
    ordered_scores = list(zip(observation_names, observed_scores.tolist(), strict=True))
    lowest = sorted(ordered_scores, key=lambda item: (item[1], item[0]))[:10]
    highest = sorted(ordered_scores, key=lambda item: (-item[1], item[0]))[:10]
    score_fingerprint = _standalone_dgidb_json_sha256(
        {
            "schema": "openbio-drug-score/v1",
            "observations": observation_names,
            "scores": [float(value) for value in observed_scores],
        }
    )

    def package_version(name):
        try:
            return importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            return "not-installed"

    parameters = {
        "producer_node_id": "OpenBioSingleCellDrugScores",
        "drug": drug,
        "source": source_kind,
        "layer_name": layer_name,
        "output_key": output_key,
        "min_matched_targets": min_matched_targets,
        "overwrite_existing": overwrite_existing,
        "resource_artifact_fingerprint_sha256": actual_resource_fingerprint,
        "resource_raw_file_sha256": resource_artifact_metadata["raw_file_sha256"],
        "resource_version": resource_metadata["version"],
        "expression_fingerprint_sha256": input_expression_fingerprint,
        "score_fingerprint_sha256": score_fingerprint,
        "fixed_policy": {"method": "mean", "nested": False, "statistical_test": False},
    }
    key_results = {
        "drug": drug,
        "cells": len(observation_names),
        "features": len(feature_names),
        "requested_target_count": len(requested_targets),
        "matched_target_count": len(matched_targets),
        "unmatched_target_count": len(unmatched_targets),
        "requested_targets": requested_targets,
        "matched_targets": matched_targets,
        "unmatched_targets": unmatched_targets,
        "contributing_sources": source_values,
        "interaction_evidence": evidence_values,
        "score_distribution": {
            "minimum": float(np.min(observed_scores)),
            "first_quartile": float(np.quantile(observed_scores, 0.25)),
            "median": float(np.median(observed_scores)),
            "mean": float(np.mean(observed_scores)),
            "third_quartile": float(np.quantile(observed_scores, 0.75)),
            "maximum": float(np.max(observed_scores)),
            "standard_deviation": float(np.std(observed_scores, ddof=1)) if len(observed_scores) > 1 else 0.0,
        },
        "lowest_scores": [
            {"observation": observation, "score": float(value)} for observation, value in lowest
        ],
        "highest_scores": [
            {"observation": observation, "score": float(value)} for observation, value in highest
        ],
        "expression_state": resolved_state,
        "expression_state_evidence": state_evidence,
        "all_zero_observation_count": zero_observation_count,
        "resource": {
            "metadata": resource_metadata,
            "raw_file_sha256": resource_artifact_metadata["raw_file_sha256"],
            "canonical_content_fingerprint_sha256": resource_artifact_metadata[
                "canonical_content_fingerprint_sha256"
            ],
            "artifact_fingerprint_sha256": actual_resource_fingerprint,
        },
        "score_fingerprint_sha256": score_fingerprint,
    }
    warnings.extend(
        [
            "This per-cell mean target-expression score is descriptive and performs no Sample-level Condition test.",
            "A high target-expression score is not drug sensitivity, response, efficacy, mechanism direction, dose, safety, or a treatment recommendation.",
        ]
    )
    if unmatched_targets:
        warnings.append(
            f"{len(unmatched_targets):,} of {len(requested_targets):,} pinned DGIdb targets were absent from the exact expression feature axis."
        )
    if zero_observation_count:
        warnings.append(
            f"{zero_observation_count:,} observations were all zero on the selected expression source; their "
            "reported mean target-expression score is zero and may reflect technical absence rather than biology."
        )
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellDrugScores",
        "methods": (
            "For one exact drug identifier, Pertpy 1.3 Enrichment.score was run on a private AnnData copy with "
            "method='mean', nested=False, and the pinned DGIdb HGNC target set. The returned score and matched-target "
            "scratch fields were extracted, independently verified against the arithmetic mean, and only the finite "
            "one-column score was stored in obs."
        ),
        "results": (
            f"Scored {len(observation_names):,} cells for {drug!r} using {len(matched_targets):,} of "
            f"{len(requested_targets):,} pinned targets; mean score was {float(np.mean(observed_scores)):.6g}. "
            "No inferential test was performed."
        ),
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": [
            "Cells are observations, not independent biological replicates; downstream Condition comparisons require Sample-level inference within a defined population and Technical batch handling.",
            "DGIdb records interaction claims without therapeutic direction; target expression alone cannot establish benefit or signature reversal.",
            "Mean target expression gives equal weight to every matched target and does not model dose, target affinity, cell accessibility, direction, or off-target effects.",
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
                "citation": f"{resource_metadata['name']} {resource_metadata['version']}: {resource_metadata['citation']}",
                "url": resource_metadata["download_url"],
                "doi": None,
                "kind": "resource_snapshot",
            },
        ],
        "software_versions": {
            "python": platform.python_version(),
            "openbio-singlecell": str(openbio_version),
            "pertpy": str(version),
            "anndata": package_version("anndata"),
            "scanpy": package_version("scanpy"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return output, summary


def run_drug_score(
    adata: AnnData,
    resource: DGIdbResource,
    *,
    copy_resource: bool = True,
    **parameters: Any,
) -> tuple[AnnData, dict[str, Any]]:
    table, metadata, accounting, artifact = validate_dgidb_resource(
        resource, copy_payload=copy_resource
    )
    return _standalone_run_drug_score(
        adata,
        table,
        metadata,
        accounting,
        artifact,
        expected_resource_artifact_fingerprint=artifact["artifact_fingerprint_sha256"],
        **parameters,
    )


def drug_score_code(
    *,
    resource_metadata: dict[str, Any],
    resource_accounting: dict[str, Any],
    resource_artifact_metadata: dict[str, Any],
    **parameters: Any,
) -> str:
    implementations = "\n\n".join(
        (
            dgidb_portable_validation_code(),
            textwrap.dedent(inspect.getsource(_standalone_run_drug_score)).strip(),
        )
    )
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in parameters.items())
    return f'''from __future__ import annotations

{implementations}


def run_drug_score(adata, resource_table, pertpy_module=None):
    return _standalone_run_drug_score(
        adata,
        resource_table,
        {resource_metadata!r},
        {resource_accounting!r},
        {resource_artifact_metadata!r},
        expected_resource_artifact_fingerprint={resource_artifact_metadata["artifact_fingerprint_sha256"]!r},
        {rendered},
        pertpy_module=pertpy_module,
    )
'''


__all__ = ["drug_score_code", "run_drug_score"]
