from __future__ import annotations

import inspect
import textwrap
from typing import TYPE_CHECKING, Any

from .enrichment_artifacts import (
    _standalone_validate_pinned_enrichment_pair,
    pinned_enrichment_validation_code,
)
from .enrichment_resource import (
    _standalone_load_gene_set_resource,
    gene_set_resource_code,
)

if TYPE_CHECKING:
    from pandas import DataFrame


RANKED_GSEA_COLUMNS = [
    "comparison",
    "gene_set",
    "normalized_enrichment_score",
    "p_adjusted",
    "rank",
    "significant",
    "set_size_before_universe",
    "set_size_in_universe",
]


def _standalone_ranked_gsea_software_versions(
    *,
    openbio_version,
    decoupler_version,
    pandas_version,
    numpy_version,
):
    import platform

    return {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
        "decoupler": str(decoupler_version),
        "pandas": str(pandas_version),
        "numpy": str(numpy_version),
    }


def _standalone_ranked_gsea_summary(diagnostics):
    """Build one strict report-ready ranked-GSEA summary."""
    import json
    import math
    from collections.abc import Mapping

    if not isinstance(diagnostics, Mapping):
        raise TypeError("Ranked GSEA summary diagnostics must be a mapping.")
    required = {
        "parameters",
        "warnings",
        "software_versions",
        "resource",
        "ranking",
        "tested_gene_sets",
        "significant_gene_sets",
        "positive_gene_sets",
        "negative_gene_sets",
        "top_positive",
        "top_negative",
    }
    missing = sorted(required - set(diagnostics))
    if missing:
        raise ValueError(f"Ranked GSEA summary diagnostics are missing required fields: {missing}.")

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
        if isinstance(value, Mapping):
            return {str(key): plain_json(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
            sequence = value.tolist() if hasattr(value, "tolist") else value
            return [plain_json(item) for item in sequence]
        return str(value)

    parameters = plain_json(diagnostics["parameters"])
    key_results = {
        key: plain_json(value)
        for key, value in diagnostics.items()
        if key not in {"parameters", "warnings", "software_versions"}
    }
    ranking = key_results["ranking"]
    resource = key_results["resource"]
    metadata = resource["metadata"]
    scope = ranking["evidence_scope"]
    formal_permutation_inference_valid = bool(key_results["formal_permutation_inference_valid"])
    methods = (
        "A single complete, structurally validated and fingerprint-bound continuous feature ranking was tested "
        "with decoupler 2.2 "
        "mt.gsea using normalized enrichment scores and its returned Benjamini-Hochberg-adjusted empirical "
        "permutation p-values. Gene sets were intersected with the exact paired tested-gene universe before "
        "the size bounds; raw=False, empty=False, and verbose=False were fixed."
    )
    results = (
        f"For {ranking['comparison']!r}, tested {int(key_results['tested_gene_sets']):,} complete-family "
        f"gene sets using {int(key_results['permutations']):,} permutations and seed "
        f"{int(key_results['random_seed'])}. {int(key_results['significant_gene_sets']):,} sets had backend-returned "
        f"adjusted values <= 0.05; {int(key_results['positive_gene_sets']):,} NES values were positive and "
        f"{int(key_results['negative_gene_sets']):,} were negative. The supplied ranking scope was "
        f"{scope!r}; no raw p-values or leading-edge members are exposed by this backend result."
    )
    if not formal_permutation_inference_valid:
        results += (
            " Because seed 0 disables permutation shuffling in decoupler 2.2, adjusted-value threshold counts "
            "are computational flags and are not valid formal permutation evidence."
        )
    references = [
        {
            "citation": (
                "Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for "
                "interpreting genome-wide expression profiles. PNAS. 2005;102:15545-15550."
            ),
            "url": "https://doi.org/10.1073/pnas.0506580102",
            "kind": "method",
            "doi": "10.1073/pnas.0506580102",
        },
        {
            "citation": (
                "Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer "
                "biological activities from omics data. Bioinformatics Advances. 2022;2:vbac016."
            ),
            "url": "https://doi.org/10.1093/bioadv/vbac016",
            "kind": "software",
            "doi": "10.1093/bioadv/vbac016",
        },
        {
            "citation": "decoupler developers. decoupler.mt.gsea official API and implementation.",
            "url": "https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.gsea.html",
            "kind": "software_documentation",
            "doi": None,
        },
        {
            "citation": (
                "Benjamini Y, Hochberg Y. Controlling the false discovery rate. Journal of the Royal "
                "Statistical Society Series B. 1995;57:289-300."
            ),
            "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
            "kind": "method",
            "doi": "10.1111/j.2517-6161.1995.tb02031.x",
        },
        {
            "citation": (f"{metadata['name']} {metadata['version']} ({metadata['date']}): {metadata['citation']}"),
            "url": f"urn:sha256:{resource['sha256']}",
            "kind": "resource",
            "doi": None,
        },
    ]
    limitations = [
        "GSEA interpretation depends on the upstream ranking statistic, its direction, the exact tested universe, identifier compatibility, and resource coverage.",
        "decoupler mt.gsea returns adjusted empirical p-values and NES but does not expose raw p-values or leading-edge genes through this interface.",
        "Overlapping gene sets create dependent hypotheses; adjusted p-values do not make pathway labels independent biological truths.",
        "The producer ranking order is fingerprinted, but decoupler 2.2 deterministically reorders exact-score ties with its internal seed-0 tie breaker before enrichment; the caller random_seed controls permutations, not that tie break.",
        "Generic producer identity, approval, scope, inference unit, replicate awareness, and Technical-batch handling are caller declarations; this node validates structure and fingerprints, not upstream scientific authority.",
    ]
    if scope == "cluster_marker_evidence":
        limitations.append(
            "Cluster marker GSEA is exploratory annotation evidence; cells are not independent Sample replicates and the result is not Condition inference."
        )
    elif scope == "condition_contrast":
        limitations.append(
            "This enrichment node does not refit the upstream Condition model; validity depends on its Sample-level design and Technical batch audit."
        )
    else:
        limitations.append(
            "The upstream evidence scope is absent or unfamiliar, so the result is computationally valid but has no OpenBio-validated formal inference interpretation."
        )
    if int(key_results["permutations"]) < 1000:
        limitations.append(
            "Fewer than 1,000 requested permutations provide coarse empirical-p resolution and unstable tail estimates; this expert choice is reported rather than rejected."
        )
    if not formal_permutation_inference_valid:
        limitations.append(
            "With random_seed=0, decoupler 2.2 skips permutation shuffling; NES values remain computational outputs, but returned adjusted values do not support formal permutation inference."
        )
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellRankedGSEA",
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": [str(value) for value in diagnostics["warnings"]],
        "limitations": limitations,
        "references": references,
        "software_versions": plain_json(diagnostics["software_versions"]),
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def _standalone_run_ranked_gsea(
    table_input,
    universe_input,
    *,
    artifact_family,
    comparison,
    comparison_column,
    gene_column,
    score_column,
    evidence_scope,
    inference_unit,
    direction,
    replicate_aware,
    technical_batch_handling,
    upstream_parameters,
    expected_analysis_fingerprint,
    expected_ranking_fingerprint,
    expected_universe_fingerprint,
    expected_table_content_fingerprint,
    expected_universe_content_fingerprint,
    resource_path,
    requested_resource_path,
    resource_metadata_json,
    provenance_warnings=(),
    provenance_declarations=None,
    expected_resource_sha256=None,
    preloaded_resource=None,
    source_column="source",
    target_column="target",
    min_targets=15,
    max_targets=500,
    n_permutations=1000,
    random_seed=1,
    max_output_rows=100_000,
    openbio_version="not-installed",
    decoupler_module=None,
):
    """Run one complete ranked GSEA family with pinned inputs and resource bytes."""
    import hashlib
    import importlib
    import inspect as runtime_inspect
    import json
    from collections.abc import Mapping
    from importlib import metadata as importlib_metadata

    import numpy as np
    import pandas as pd

    operation = "Ranked GSEA"
    for name, value, minimum in (
        ("min_targets", min_targets, 1),
        ("max_targets", max_targets, 1),
        ("n_permutations", n_permutations, 2),
        ("random_seed", random_seed, 0),
        ("max_output_rows", max_output_rows, 1),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{operation} {name} must be an integer of at least {minimum}.")
    if min_targets > max_targets:
        raise ValueError(f"{operation} min_targets cannot exceed max_targets.")
    if not isinstance(upstream_parameters, Mapping):
        raise TypeError(f"{operation} upstream_parameters must be a mapping.")
    if not isinstance(provenance_warnings, (list, tuple)) or not all(
        isinstance(value, str) for value in provenance_warnings
    ):
        raise TypeError(f"{operation} provenance_warnings must be a list or tuple of strings.")
    runtime_provenance_warnings = list(provenance_warnings)
    if provenance_declarations is None:
        provenance_declarations = {
            "validation_scope": (
                "marker_v2_structure_fingerprints_axes_and_direct_operation"
                if artifact_family == "marker_v2"
                else "structure_current_content_fingerprints_and_axes_only"
            ),
            "scientific_metadata_status": (
                "cluster_marker_inference_is_exploratory"
                if artifact_family == "marker_v2"
                else "caller_declared_not_programmatically_verified"
            ),
            "formal_condition_interpretation_validated": False,
        }
    if not isinstance(provenance_declarations, Mapping):
        raise TypeError(f"{operation} provenance_declarations must be a mapping.")
    provenance_declarations = dict(provenance_declarations)
    if artifact_family == "generic_v1":
        generic_warning = (
            "Generic evidence producer identity, approval, scope, inference unit, replicate awareness, and "
            "Technical-batch handling are caller declarations; OpenBio validated structure, current-content "
            "fingerprints, and table/universe axes only."
        )
        if generic_warning not in runtime_provenance_warnings:
            runtime_provenance_warnings.insert(0, generic_warning)

    selected, universe_genes = _standalone_validate_pinned_enrichment_pair(
        table_input,
        universe_input,
        artifact_family=artifact_family,
        purpose="ranked",
        comparison=comparison,
        comparison_column=comparison_column,
        gene_column=gene_column,
        score_column=score_column,
        expected_analysis_fingerprint=expected_analysis_fingerprint,
        expected_ranking_fingerprint=expected_ranking_fingerprint,
        expected_universe_fingerprint=expected_universe_fingerprint,
        expected_table_content_fingerprint=expected_table_content_fingerprint,
        expected_universe_content_fingerprint=expected_universe_content_fingerprint,
    )
    scores = selected[score_column].to_numpy(dtype=float)
    genes = selected[gene_column].tolist()
    tie_blocks = []
    start = 0
    while start < len(scores):
        stop = start + 1
        while stop < len(scores) and scores[stop] == scores[start]:
            stop += 1
        if stop - start > 1:
            tie_blocks.append(
                {
                    "score": float(scores[start]),
                    "start_rank": start + 1,
                    "end_rank": stop,
                    "gene_count": stop - start,
                    "genes_preview": genes[start : min(stop, start + 10)],
                    "preview_truncated": stop - start > 10,
                }
            )
        start = stop

    if preloaded_resource is None:
        resource = _standalone_load_gene_set_resource(
            resource_path,
            resource_metadata_json,
            source_column,
            target_column,
            universe_genes,
            min_targets,
            expected_sha256=expected_resource_sha256,
            operation=operation,
        )
    elif not isinstance(preloaded_resource, Mapping):
        raise TypeError(f"{operation} preloaded_resource must be a validated mapping.")
    else:
        resource = preloaded_resource
    for identity_field in ("organism", "identifier_namespace"):
        upstream_identity = upstream_parameters.get(identity_field)
        if upstream_identity is not None and upstream_identity != resource["metadata"][identity_field]:
            runtime_provenance_warnings.append(
                f"{operation} resource {identity_field} conflicts with upstream evidence: "
                f"{resource['metadata'][identity_field]!r} != {upstream_identity!r}. Exact identifier matching "
                "was used, but biological compatibility is a caller responsibility."
            )

    whole_universe_sources = []
    max_size_sources = []
    retained_sources = []
    universe_count = len(universe_genes)
    for source in resource["retained_sources"]:
        size = len(resource["targets_in_universe"][source])
        if size == universe_count:
            whole_universe_sources.append(source)
        elif size > max_targets:
            max_size_sources.append(source)
        else:
            retained_sources.append(source)
    if not retained_sources:
        raise ValueError(
            f"{operation} has no nondegenerate resource set within [{min_targets}, {max_targets}] targets "
            "after exact universe intersection."
        )
    if len(retained_sources) > max_output_rows:
        raise ValueError(
            f"{operation} would return {len(retained_sources):,} complete-family rows, exceeding "
            f"max_output_rows={max_output_rows:,}."
        )
    retained_set = set(retained_sources)
    network = resource["network"].loc[resource["network"]["source"].isin(retained_set), ["source", "target"]].copy()
    if network.empty or set(network["source"].tolist()) != retained_set:
        raise RuntimeError(f"{operation} failed to retain the exact eligible resource-set family.")
    ranked_frame = pd.DataFrame([scores], index=[comparison], columns=genes)

    if decoupler_module is None:
        try:
            decoupler_module = importlib.import_module("decoupler")
        except (ImportError, OSError) as exc:
            raise RuntimeError(f"{operation} requires decoupler 2.2, but it is unavailable ({exc}).") from exc
    version = getattr(decoupler_module, "__version__", None)
    if not isinstance(version, str) or not version.strip():
        try:
            version = importlib_metadata.version("decoupler")
        except importlib_metadata.PackageNotFoundError as exc:
            raise RuntimeError(f"{operation} cannot determine the installed decoupler version.") from exc
    try:
        version_parts = version.split(".")
        version_major = int(version_parts[0])
        version_minor = int(version_parts[1])
    except (IndexError, TypeError, ValueError) as exc:
        raise RuntimeError(f"{operation} received invalid decoupler version {version!r}.") from exc
    if version_major != 2 or version_minor != 2:
        raise RuntimeError(f"{operation} requires decoupler 2.2.x; installed version is {version!r}.")
    gsea = getattr(getattr(decoupler_module, "mt", None), "gsea", None)
    if not callable(gsea):
        raise RuntimeError(f"{operation} requires the public decoupler.mt.gsea API.")
    try:
        signature = runtime_inspect.signature(gsea)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{operation} could not inspect decoupler.mt.gsea.") from exc
    public_keywords = {"data", "net", "tmin", "raw", "empty", "verbose"}
    backend_keywords = {"times", "seed"}
    supports_kwargs = any(
        parameter.kind is runtime_inspect.Parameter.VAR_KEYWORD for parameter in signature.parameters.values()
    )
    missing_public = sorted(public_keywords - set(signature.parameters))
    if missing_public:
        raise RuntimeError(
            f"{operation} decoupler.mt.gsea is missing required public parameters: {missing_public}."
        )
    missing_backend = backend_keywords - set(signature.parameters)
    if missing_backend and supports_kwargs:
        gsea_func = getattr(gsea, "func", None)
        if not callable(gsea_func):
            raise RuntimeError(
                f"{operation} cannot verify decoupler.mt.gsea passthrough parameters {sorted(missing_backend)}."
            )
        try:
            func_signature = runtime_inspect.signature(gsea_func)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{operation} could not inspect decoupler.mt.gsea.func.") from exc
        missing_backend -= set(func_signature.parameters)
    if missing_backend:
        missing = sorted(missing_backend)
        raise RuntimeError(f"{operation} decoupler.mt.gsea is missing required parameters: {missing}.")

    backend = gsea(
        data=ranked_frame,
        net=network,
        tmin=min_targets,
        raw=False,
        empty=False,
        verbose=False,
        times=n_permutations,
        seed=random_seed,
    )
    if not isinstance(backend, tuple) or len(backend) != 2:
        raise RuntimeError(f"{operation} decoupler.mt.gsea must return (scores, p_adjusted).")
    normalized, adjusted = backend
    for name, frame in (("scores", normalized), ("p_adjusted", adjusted)):
        if not isinstance(frame, pd.DataFrame):
            raise RuntimeError(f"{operation} backend {name} output must be a pandas DataFrame.")
        if frame.shape != (1, len(retained_sources)):
            raise RuntimeError(
                f"{operation} backend {name} shape {frame.shape} differs from the complete family "
                f"(1, {len(retained_sources)})."
            )
        if bool(frame.columns.duplicated().any()) or set(frame.columns.tolist()) != retained_set:
            raise RuntimeError(f"{operation} backend {name} sources differ from the retained resource family.")
    try:
        nes = normalized.loc[:, retained_sources].iloc[0].to_numpy(dtype=float)
        p_adjusted = adjusted.loc[:, retained_sources].iloc[0].to_numpy(dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError(f"{operation} backend returned nonnumeric or inaccessible values.") from exc
    if not bool(np.isfinite(nes).all()) or not bool(np.isfinite(p_adjusted).all()):
        raise RuntimeError(f"{operation} backend returned non-finite NES or adjusted p-values.")
    if bool(((p_adjusted < 0.0) | (p_adjusted > 1.0)).any()):
        raise RuntimeError(f"{operation} backend adjusted p-values must lie in [0, 1].")

    source_order = {source: index for index, source in enumerate(resource["source_order"])}
    rows = [
        {
            "comparison": comparison,
            "gene_set": source,
            "normalized_enrichment_score": float(score),
            "p_adjusted": float(padj),
            "significant": bool(padj <= 0.05),
            "set_size_before_universe": len(resource["targets_before"][source]),
            "set_size_in_universe": len(resource["targets_in_universe"][source]),
            "_source_order": source_order[source],
        }
        for source, score, padj in zip(retained_sources, nes.tolist(), p_adjusted.tolist(), strict=True)
    ]
    evidence = pd.DataFrame(rows).sort_values(
        ["p_adjusted", "normalized_enrichment_score", "_source_order"],
        ascending=[True, False, True],
        key=lambda column: column.abs() if column.name == "normalized_enrichment_score" else column,
        kind="mergesort",
        ignore_index=True,
    )
    ranks = []
    previous_key = None
    current_rank = 0
    for row in evidence.itertuples(index=False):
        scientific_key = (float(row.p_adjusted), float(row.normalized_enrichment_score))
        if scientific_key != previous_key:
            current_rank += 1
            previous_key = scientific_key
        ranks.append(current_rank)
    evidence["rank"] = ranks
    evidence = evidence[
        [
            "comparison",
            "gene_set",
            "normalized_enrichment_score",
            "p_adjusted",
            "rank",
            "significant",
            "set_size_before_universe",
            "set_size_in_universe",
        ]
    ]
    if len(evidence) != len(retained_sources) or bool(evidence["gene_set"].duplicated().any()):
        raise RuntimeError(f"{operation} failed to retain the complete unique tested family.")

    def result_sha256(frame):
        payload = {
            "schema": "openbio-singlecell/ranked-gsea-evidence/v1",
            "resource_sha256": resource["sha256"],
            "rows": [
                [
                    row.comparison,
                    row.gene_set,
                    float(row.normalized_enrichment_score).hex(),
                    float(row.p_adjusted).hex(),
                    int(row.rank),
                    bool(row.significant),
                    int(row.set_size_before_universe),
                    int(row.set_size_in_universe),
                ]
                for row in frame.itertuples(index=False)
            ],
        }
        encoder = json.JSONEncoder(
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        digest = hashlib.sha256()
        for chunk in encoder.iterencode(payload):
            digest.update(chunk.encode("utf-8"))
        return digest.hexdigest()

    top_positive = evidence.loc[evidence["normalized_enrichment_score"] > 0].head(10).to_dict("records")
    top_negative = (
        evidence.loc[evidence["normalized_enrichment_score"] < 0]
        .sort_values(
            ["p_adjusted", "normalized_enrichment_score", "gene_set"],
            ascending=[True, True, True],
            kind="mergesort",
        )
        .head(10)
        .to_dict("records")
    )
    resource_accounting = dict(resource["accounting"])
    set_accounting = resource_accounting.pop("set_accounting")
    removed_by_min = resource_accounting.pop("sources_removed_by_min_targets")
    resource_accounting.update(
        {
            "metadata": dict(resource["metadata"]),
            "resolved_path": resource["resolved_path"],
            "requested_path": requested_resource_path,
            "size_bytes": resource["size_bytes"],
            "sha256": resource["sha256"],
            "sources_removed_by_max_targets_count": len(max_size_sources),
            "sources_removed_by_min_targets_preview": removed_by_min[:100],
            "sources_removed_by_min_targets_preview_truncated": len(removed_by_min) > 100,
            "set_accounting_preview": set_accounting[:100],
            "set_accounting_preview_truncated": len(set_accounting) > 100,
            "sources_removed_by_max_targets_preview": max_size_sources[:100],
            "sources_removed_by_max_targets_preview_truncated": len(max_size_sources) > 100,
            "whole_universe_sources_removed_count": len(whole_universe_sources),
            "whole_universe_sources_removed_preview": whole_universe_sources[:100],
            "whole_universe_sources_removed_preview_truncated": len(whole_universe_sources) > 100,
            "tested_sources_preview": retained_sources[:100],
            "tested_sources_preview_truncated": len(retained_sources) > 100,
        }
    )
    parameters = {
        "producer_node_id": "OpenBioSingleCellRankedGSEA",
        "artifact_role": "ranked_gsea_evidence_table",
        "comparison": comparison,
        "comparison_column": comparison_column,
        "gene_column": gene_column,
        "score_column": score_column,
        "evidence_scope": evidence_scope,
        "inference_unit": inference_unit,
        "direction": direction,
        "replicate_aware": replicate_aware,
        "technical_batch_handling": technical_batch_handling,
        "input_provenance_declarations": provenance_declarations,
        "analysis_fingerprint": expected_analysis_fingerprint,
        "ranking_fingerprint": expected_ranking_fingerprint,
        "universe_fingerprint": expected_universe_fingerprint,
        "input_table_content_fingerprint": expected_table_content_fingerprint,
        "universe_content_fingerprint": expected_universe_content_fingerprint,
        "resource_path": requested_resource_path,
        "resource_sha256": resource["sha256"],
        "resource_metadata": dict(resource["metadata"]),
        "source_column": source_column,
        "target_column": target_column,
        "min_targets": min_targets,
        "max_targets": max_targets,
        "n_permutations": n_permutations,
        "random_seed": random_seed,
        "max_output_rows": max_output_rows,
        "adjusted_p_significance_threshold": 0.05,
        "fixed_policy": {
            "backend": "decoupler.mt.gsea",
            "raw": False,
            "empty": False,
            "verbose": False,
            "complete_family_returned": True,
            "raw_p_values_available": False,
            "leading_edge_available": False,
        },
        "result_fingerprint": result_sha256(evidence),
    }
    warnings = [
        *runtime_provenance_warnings,
        "Resource organism, namespace, scope, license, and citation are caller declarations; SHA-256 identifies exact bytes but not biological suitability.",
        "The backend reports adjusted empirical permutation p-values; no raw p-value or leading-edge column is available from this interface.",
    ]
    if n_permutations < 1000:
        warnings.append(
            f"Only {n_permutations:,} permutations were requested. Empirical-p resolution is coarse and tail "
            "estimates are unstable; at least 1,000 permutations is recommended for formal reporting."
        )
    if random_seed == 0:
        warnings.append(
            "random_seed=0 is accepted by decoupler 2.2 but disables permutation shuffling. NES values are "
            "computable, while backend-returned adjusted values are not valid formal permutation evidence."
        )
    if evidence_scope == "cluster_marker_evidence":
        warnings.append(
            "This is exploratory Cluster marker evidence for annotation review, not Sample-level Condition inference."
        )
    if tie_blocks:
        warnings.append(
            f"The ranking contains {len(tie_blocks):,} exact-score tie blocks. Producer order is fingerprinted, "
            "but decoupler 2.2 applies a deterministic internal seed-0 tie break before enrichment; random_seed "
            "controls permutations, not this tie break."
        )
    if resource["accounting"]["duplicate_pairs_removed"]:
        warnings.append(
            f"Collapsed {resource['accounting']['duplicate_pairs_removed']:,} exact duplicate resource edges."
        )
    if max_size_sources:
        warnings.append(f"Removed {len(max_size_sources):,} sets above max_targets after universe intersection.")
    if whole_universe_sources:
        warnings.append(
            f"Removed {len(whole_universe_sources):,} whole-universe sets because GSEA requires a non-set complement."
        )
    diagnostics = {
        "provenance": provenance_declarations,
        "ranking": {
            "comparison": comparison,
            "evidence_scope": evidence_scope,
            "inference_unit": inference_unit,
            "direction": direction,
            "replicate_aware": replicate_aware,
            "technical_batch_handling": technical_batch_handling,
            "gene_count": len(genes),
            "score_min": float(scores.min()),
            "score_max": float(scores.max()),
            "exact_tie_block_count": len(tie_blocks),
            "exact_tie_blocks_preview": tie_blocks[:20],
            "tie_preview_truncated": len(tie_blocks) > 20,
            "producer_row_order_fingerprinted": True,
            "backend_exact_tie_policy": "decoupler_2.2_deterministic_seed_0",
            "backend_tie_break_seed": 0,
            "analysis_fingerprint_sha256": expected_analysis_fingerprint,
            "ranking_fingerprint_sha256": expected_ranking_fingerprint,
            "universe_fingerprint_sha256": expected_universe_fingerprint,
            "table_content_fingerprint_sha256": expected_table_content_fingerprint,
            "universe_content_fingerprint_sha256": expected_universe_content_fingerprint,
        },
        "resource": resource_accounting,
        "tested_gene_sets": len(evidence),
        "significant_gene_sets": int(evidence["significant"].sum()),
        "positive_gene_sets": int((evidence["normalized_enrichment_score"] > 0).sum()),
        "negative_gene_sets": int((evidence["normalized_enrichment_score"] < 0).sum()),
        "zero_nes_gene_sets": int((evidence["normalized_enrichment_score"] == 0).sum()),
        "top_positive": top_positive,
        "top_negative": top_negative,
        "permutations": n_permutations,
        "random_seed": random_seed,
        "permutation_null_status": "not_shuffled_seed_zero" if random_seed == 0 else "shuffled",
        "formal_permutation_inference_valid": random_seed != 0,
        "backend": {
            "api": "decoupler.mt.gsea",
            "version": str(version),
            "return_semantics": "normalized_enrichment_score, Benjamini-Hochberg-adjusted empirical permutation p-value",
        },
        "parameters": parameters,
        "warnings": warnings,
        "software_versions": _standalone_ranked_gsea_software_versions(
            openbio_version=openbio_version,
            decoupler_version=str(version),
            pandas_version=pd.__version__,
            numpy_version=np.__version__,
        ),
    }
    json.dumps(diagnostics, ensure_ascii=False, allow_nan=False)
    return evidence, diagnostics


def run_ranked_gsea_evidence(table: Any, universe: Any, **kwargs: Any) -> tuple[DataFrame, dict[str, Any]]:
    return _standalone_run_ranked_gsea(table, universe, **kwargs)


def build_ranked_gsea_summary(diagnostics: dict[str, Any]) -> dict[str, Any]:
    return _standalone_ranked_gsea_summary(diagnostics)


def ranked_gsea_code(**parameters: Any) -> str:
    implementations = "\n\n".join(
        (
            pinned_enrichment_validation_code(),
            gene_set_resource_code(),
            textwrap.dedent(inspect.getsource(_standalone_ranked_gsea_software_versions)).strip(),
            textwrap.dedent(inspect.getsource(_standalone_ranked_gsea_summary)).strip(),
            textwrap.dedent(inspect.getsource(_standalone_run_ranked_gsea)).strip(),
        )
    )
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in parameters.items())
    return f"""from __future__ import annotations

{implementations}


def run_ranked_gsea(table, universe, decoupler_module=None):
    evidence, diagnostics = _standalone_run_ranked_gsea(
        table,
        universe,
        {rendered},
        decoupler_module=decoupler_module,
    )
    return evidence, _standalone_ranked_gsea_summary(diagnostics)
"""


__all__ = [
    "RANKED_GSEA_COLUMNS",
    "build_ranked_gsea_summary",
    "ranked_gsea_code",
    "run_ranked_gsea_evidence",
]
