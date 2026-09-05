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


GENERIC_ORA_COLUMNS = [
    "scope",
    "comparison",
    "gene_set",
    "selected_count",
    "universe_count",
    "set_size_before_universe",
    "set_size_in_universe",
    "overlap_count",
    "overlap_genes",
    "a",
    "b",
    "c",
    "d",
    "log_odds_ratio",
    "p_value",
    "p_adjusted",
    "rank",
    "positive_enrichment",
    "passes_min_overlap",
    "significant",
    "candidate",
    "candidate_status",
]


def _standalone_generic_ora_software_versions(
    *,
    openbio_version,
    decoupler_version,
    scipy_version,
    pandas_version,
    numpy_version,
):
    import platform

    return {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
        "decoupler": str(decoupler_version),
        "scipy": str(scipy_version),
        "pandas": str(pandas_version),
        "numpy": str(numpy_version),
    }


def _standalone_generic_ora_summary(diagnostics):
    """Build the strict report-ready selected-set ORA summary."""
    import json
    import math
    from collections.abc import Mapping

    if not isinstance(diagnostics, Mapping):
        raise TypeError("Gene-set ORA summary diagnostics must be a mapping.")
    required = {
        "parameters",
        "warnings",
        "software_versions",
        "resource",
        "selection",
        "tested_gene_sets",
        "candidate_gene_sets",
        "significant_gene_sets",
        "top_results",
    }
    missing = sorted(required - set(diagnostics))
    if missing:
        raise ValueError(f"Gene-set ORA summary diagnostics are missing required fields: {missing}.")

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
    selection = key_results["selection"]
    resource = key_results["resource"]
    metadata = resource["metadata"]
    methods = (
        "One explicit, structurally validated and fingerprint-bound gene set was tested against local resource sets "
        "intersected with its exact paired tested-gene universe. decoupler 2.2 mt.query_set used "
        "alternative='greater', n_bg equal to the universe size, Haldane-Anscombe correction 0.5, and the "
        "post-intersection min_targets. Every 2x2 contingency, one-sided Fisher p-value, corrected log odds "
        "ratio, and Benjamini-Hochberg adjusted p-value was independently reconstructed and verified."
    )
    results = (
        f"For {selection['comparison']!r}, {int(selection['selected_count']):,} selected genes were tested "
        f"within {int(selection['universe_count']):,} universe genes across "
        f"{int(key_results['tested_gene_sets']):,} retained resource sets. "
        f"{int(key_results['candidate_gene_sets']):,} sets met the positive-effect, overlap >= "
        f"{int(key_results['thresholds']['min_overlap'])}, and adjusted p <= "
        f"{float(key_results['thresholds']['max_p_adjusted']):g} descriptive gates. The complete tested "
        "resource family is retained in the table."
    )
    references = [
        {
            "citation": (
                "Fisher RA. On the interpretation of chi-square from contingency tables, and the "
                "calculation of P. Journal of the Royal Statistical Society. 1922;85:87-94."
            ),
            "url": "https://doi.org/10.2307/2340521",
            "kind": "method",
            "doi": "10.2307/2340521",
        },
        {
            "citation": (
                "Haldane JBS. The estimation and significance of the logarithm of a ratio of frequencies. "
                "Annals of Human Genetics. 1956;20:309-311."
            ),
            "url": "https://doi.org/10.1111/j.1469-1809.1955.tb01285.x",
            "kind": "method",
            "doi": "10.1111/j.1469-1809.1955.tb01285.x",
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
            "citation": (
                "Timmons JA, Szkop KJ, Gallagher IJ. Multiple sources of bias confound functional "
                "enrichment analysis. Genome Biology. 2015;16:186."
            ),
            "url": "https://doi.org/10.1186/s13059-015-0761-7",
            "kind": "practice",
            "doi": "10.1186/s13059-015-0761-7",
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
            "citation": "decoupler developers. decoupler.mt.query_set official API and implementation.",
            "url": "https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.query_set.html",
            "kind": "software_documentation",
            "doi": None,
        },
        {
            "citation": (
                f"{metadata.get('name', 'Gene-set resource')} {metadata.get('version', '')} "
                f"({metadata.get('date', 'date not supplied')}): {metadata.get('citation', 'citation not supplied')}"
            ),
            "url": f"urn:sha256:{resource['sha256']}",
            "kind": "resource",
            "doi": None,
        },
    ]
    limitations = [
        "ORA depends on the upstream selection rule, exact tested universe, identifier compatibility, and resource coverage; it cannot repair upstream confounding or pseudoreplication.",
        "Overlapping gene sets create dependent hypotheses, and data-dependent selection limits confirmatory interpretation.",
        "Haldane-corrected log odds ratios remain descriptive, especially for sparse or zero-overlap contingencies.",
        "Generic producer identity, approval, scope, inference unit, replicate awareness, and Technical-batch handling are caller declarations; this node validates structure and fingerprints, not upstream scientific authority.",
    ]
    if selection["evidence_scope"] == "cluster_marker_evidence":
        limitations.append(
            "Cluster marker ORA is exploratory annotation evidence, not a cell-type probability, Curated annotation, or Condition inference."
        )
    elif selection["evidence_scope"] == "condition_contrast":
        limitations.append(
            "This ORA node does not refit the upstream Condition model; interpretation depends on its Sample-level replicate design and Technical batch audit."
        )
    else:
        limitations.append(
            "The upstream evidence scope is absent or unfamiliar, so the result is computationally valid but has no OpenBio-validated formal inference interpretation."
        )
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellGeneSetOverrepresentation",
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


def _standalone_run_generic_ora(
    table_input,
    universe_input,
    *,
    artifact_family,
    comparison,
    comparison_column,
    gene_column,
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
    min_targets=3,
    min_overlap=2,
    max_p_adjusted=0.05,
    max_output_rows=100_000,
    openbio_version="not-installed",
    decoupler_module=None,
):
    """Run one selected-set ORA family and independently verify decoupler results."""
    import hashlib
    import importlib
    import inspect as runtime_inspect
    import json
    import math
    from collections.abc import Mapping
    from importlib import metadata as importlib_metadata

    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats

    operation = "Gene Set Overrepresentation"
    for name, value in (
        ("min_targets", min_targets),
        ("min_overlap", min_overlap),
        ("max_output_rows", max_output_rows),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{operation} {name} must be a positive integer.")
    if (
        isinstance(max_p_adjusted, bool)
        or not isinstance(max_p_adjusted, (int, float))
        or not math.isfinite(float(max_p_adjusted))
        or not 0.0 <= float(max_p_adjusted) <= 1.0
    ):
        raise ValueError(f"{operation} max_p_adjusted must be finite in [0, 1].")
    max_p_adjusted = float(max_p_adjusted)
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

    selected_frame, universe_genes = _standalone_validate_pinned_enrichment_pair(
        table_input,
        universe_input,
        artifact_family=artifact_family,
        purpose="selected",
        comparison=comparison,
        comparison_column=comparison_column,
        gene_column=gene_column,
        score_column=None,
        expected_analysis_fingerprint=expected_analysis_fingerprint,
        expected_ranking_fingerprint=expected_ranking_fingerprint,
        expected_universe_fingerprint=expected_universe_fingerprint,
        expected_table_content_fingerprint=expected_table_content_fingerprint,
        expected_universe_content_fingerprint=expected_universe_content_fingerprint,
    )
    if artifact_family == "marker_v2":
        thresholds = {
            "min_log2_fold_change": upstream_parameters.get("min_log2_fold_change"),
            "min_fraction_in_group": upstream_parameters.get("min_fraction_in_group"),
            "max_fraction_reference": upstream_parameters.get("max_fraction_reference"),
            "max_p_adjusted": upstream_parameters.get("max_p_adjusted"),
        }
        for name, value in thresholds.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ValueError(f"{operation} marker-filter provenance has invalid threshold {name!r}.")
            thresholds[name] = float(value)
        numeric = selected_frame[["log2_fold_change_approx", "fraction_in_group", "fraction_reference", "p_adjusted"]]
        values = numeric.to_numpy(dtype=float)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} selected marker evidence contains non-finite values.")
        valid = (
            (numeric["log2_fold_change_approx"].to_numpy(dtype=float) > 0.0)
            & (numeric["log2_fold_change_approx"].to_numpy(dtype=float) >= thresholds["min_log2_fold_change"])
            & (numeric["fraction_in_group"].to_numpy(dtype=float) >= thresholds["min_fraction_in_group"])
            & (numeric["fraction_reference"].to_numpy(dtype=float) <= thresholds["max_fraction_reference"])
            & (numeric["p_adjusted"].to_numpy(dtype=float) <= thresholds["max_p_adjusted"])
        )
        if not bool(valid.all()):
            raise ValueError(
                f"{operation} selected marker rows violate their direct filter provenance or positive direction."
            )

    selected_genes = selected_frame[gene_column].tolist()
    selected_set = set(selected_genes)
    universe_set = set(universe_genes)
    if not selected_set or not selected_set.issubset(universe_set):
        raise ValueError(f"{operation} selected genes must be a nonempty subset of the exact universe.")

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
        resource_identity = resource["metadata"].get(identity_field)
        if upstream_identity is not None and resource_identity and upstream_identity != resource_identity:
            runtime_provenance_warnings.append(
                f"{operation} resource {identity_field} conflicts with upstream evidence: "
                f"{resource['metadata'][identity_field]!r} != {upstream_identity!r}. Exact identifier matching "
                "was used, but biological compatibility is a caller responsibility."
            )
    retained_sources = list(resource["retained_sources"])
    if len(retained_sources) > max_output_rows:
        raise ValueError(
            f"{operation} would return {len(retained_sources):,} complete-family rows, exceeding "
            f"max_output_rows={max_output_rows:,}."
        )

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
    query_set = getattr(getattr(decoupler_module, "mt", None), "query_set", None)
    if not callable(query_set):
        raise RuntimeError(f"{operation} requires the public decoupler.mt.query_set API.")
    try:
        signature = runtime_inspect.signature(query_set)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{operation} could not inspect decoupler.mt.query_set.") from exc
    required_keywords = {"features", "net", "alternative", "n_bg", "ha_corr", "tmin", "verbose"}
    supports_kwargs = any(
        parameter.kind is runtime_inspect.Parameter.VAR_KEYWORD for parameter in signature.parameters.values()
    )
    if not supports_kwargs and not required_keywords.issubset(signature.parameters):
        missing = sorted(required_keywords - set(signature.parameters))
        raise RuntimeError(f"{operation} decoupler.mt.query_set is missing required parameters: {missing}.")

    ordered_selected = [gene for gene in universe_genes if gene in selected_set]
    backend = query_set(
        features=ordered_selected,
        net=resource["network"],
        alternative="greater",
        n_bg=len(universe_genes),
        ha_corr=0.5,
        tmin=min_targets,
        verbose=False,
    )
    if not isinstance(backend, pd.DataFrame):
        raise RuntimeError(f"{operation} backend result must be a pandas DataFrame.")
    required_columns = ["source", "stat", "pval", "padj"]
    missing_columns = [column for column in required_columns if column not in backend]
    if missing_columns:
        raise RuntimeError(f"{operation} backend result is missing columns: {missing_columns}.")
    if bool(backend["source"].duplicated().any()) or set(backend["source"].tolist()) != set(retained_sources):
        raise RuntimeError(f"{operation} backend source family differs from the complete retained resource.")
    indexed = backend.set_index("source")
    for column in ("stat", "pval", "padj"):
        if pd.api.types.is_bool_dtype(indexed[column].dtype) or not pd.api.types.is_numeric_dtype(
            indexed[column].dtype
        ):
            raise RuntimeError(f"{operation} backend column {column!r} must be numeric and non-boolean.")
    backend_values = indexed.loc[retained_sources, ["stat", "pval", "padj"]].to_numpy(dtype=float)
    if not bool(np.isfinite(backend_values).all()):
        raise RuntimeError(f"{operation} backend returned non-finite values.")
    if bool(((backend_values[:, 1:] < 0.0) | (backend_values[:, 1:] > 1.0)).any()):
        raise RuntimeError(f"{operation} backend returned p-values outside [0, 1].")

    raw_rows = []
    independent_p = []
    independent_log_odds = []
    universe_count = len(universe_genes)
    universe_position = {gene: index for index, gene in enumerate(universe_genes)}
    for source in retained_sources:
        resource_set = set(resource["targets_in_universe"][source])
        overlap = sorted(selected_set & resource_set, key=universe_position.__getitem__)
        a = len(overlap)
        b = len(resource_set - selected_set)
        c = len(selected_set - resource_set)
        d = universe_count - a - b - c
        if min(a, b, c, d) < 0 or a + b + c + d != universe_count:
            raise RuntimeError(f"{operation} constructed an invalid contingency for {source!r}.")
        log_odds = float(np.log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))))
        p_value = float(stats.fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
        independent_log_odds.append(log_odds)
        independent_p.append(p_value)
        raw_rows.append(
            {
                "scope": evidence_scope,
                "comparison": comparison,
                "gene_set": source,
                "selected_count": len(selected_set),
                "universe_count": universe_count,
                "set_size_before_universe": len(resource["targets_before"][source]),
                "set_size_in_universe": len(resource_set),
                "overlap_count": a,
                "overlap_genes": json.dumps(overlap, ensure_ascii=False, separators=(",", ":")),
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "log_odds_ratio": log_odds,
                "p_value": p_value,
            }
        )
    independent_p_array = np.asarray(independent_p, dtype=float)
    independent_log_odds_array = np.asarray(independent_log_odds, dtype=float)
    independent_adjusted = np.asarray(stats.false_discovery_control(independent_p_array, method="bh"), dtype=float)
    if not bool(np.allclose(backend_values[:, 0], independent_log_odds_array, rtol=1e-10, atol=1e-12)):
        raise RuntimeError(f"{operation} backend log odds disagree with independent Haldane contingencies.")
    if not bool(np.allclose(backend_values[:, 1], independent_p_array, rtol=1e-10, atol=1e-12)):
        raise RuntimeError(f"{operation} backend p-values disagree with independent Fisher tests.")
    if not bool(np.allclose(backend_values[:, 2], independent_adjusted, rtol=1e-10, atol=1e-12)):
        raise RuntimeError(f"{operation} backend adjusted p-values disagree with complete-family BH.")
    for row, adjusted in zip(raw_rows, independent_adjusted.tolist(), strict=True):
        row["p_adjusted"] = float(adjusted)

    source_order = {source: index for index, source in enumerate(resource["source_order"])}
    evidence = pd.DataFrame(raw_rows)
    evidence["_source_order"] = evidence["gene_set"].map(source_order)
    evidence = evidence.sort_values(
        ["p_adjusted", "p_value", "log_odds_ratio", "overlap_count", "_source_order"],
        ascending=[True, True, False, False, True],
        kind="mergesort",
        ignore_index=True,
    )
    ranks = []
    previous_key = None
    current_rank = 0
    for row in evidence.itertuples(index=False):
        scientific_key = (
            float(row.p_adjusted),
            float(row.p_value),
            float(row.log_odds_ratio),
            int(row.overlap_count),
        )
        if scientific_key != previous_key:
            current_rank += 1
            previous_key = scientific_key
        ranks.append(current_rank)
    evidence["rank"] = ranks
    evidence["positive_enrichment"] = evidence["log_odds_ratio"] > 0.0
    evidence["passes_min_overlap"] = evidence["overlap_count"] >= min_overlap
    evidence["significant"] = evidence["p_adjusted"] <= max_p_adjusted
    evidence["candidate"] = evidence["positive_enrichment"] & evidence["passes_min_overlap"] & evidence["significant"]
    statuses = []
    for row in evidence.itertuples(index=False):
        if not bool(row.positive_enrichment):
            statuses.append("nonpositive")
        elif not bool(row.passes_min_overlap):
            statuses.append("insufficient_overlap")
        elif not bool(row.significant):
            statuses.append("not_significant")
        else:
            statuses.append("candidate")
    evidence["candidate_status"] = statuses
    candidate_indices = evidence.index[evidence["candidate"]].tolist()
    if candidate_indices:
        best_rank = min(int(evidence.at[index, "rank"]) for index in candidate_indices)
        best = [index for index in candidate_indices if int(evidence.at[index, "rank"]) == best_rank]
        if len(best) > 1:
            evidence.loc[best, "candidate_status"] = "tied_best_candidate"
    evidence = evidence.drop(columns=["_source_order"])[
        [
            "scope",
            "comparison",
            "gene_set",
            "selected_count",
            "universe_count",
            "set_size_before_universe",
            "set_size_in_universe",
            "overlap_count",
            "overlap_genes",
            "a",
            "b",
            "c",
            "d",
            "log_odds_ratio",
            "p_value",
            "p_adjusted",
            "rank",
            "positive_enrichment",
            "passes_min_overlap",
            "significant",
            "candidate",
            "candidate_status",
        ]
    ]
    if len(evidence) != len(retained_sources) or bool(evidence["gene_set"].duplicated().any()):
        raise RuntimeError(f"{operation} failed to retain the complete unique tested family.")
    for row in evidence.itertuples(index=False):
        if row.a != row.overlap_count or row.a + row.b != row.set_size_in_universe:
            raise RuntimeError(f"{operation} contingency identities failed.")
        if row.a + row.c != row.selected_count or row.a + row.b + row.c + row.d != row.universe_count:
            raise RuntimeError(f"{operation} contingency identities failed.")

    def result_sha256(frame):
        payload = {
            "schema": "openbio-singlecell/generic-ora-evidence/v1",
            "resource_sha256": resource["sha256"],
            "rows": [
                [
                    row.scope,
                    row.comparison,
                    row.gene_set,
                    int(row.selected_count),
                    int(row.universe_count),
                    int(row.set_size_before_universe),
                    int(row.set_size_in_universe),
                    int(row.overlap_count),
                    row.overlap_genes,
                    int(row.a),
                    int(row.b),
                    int(row.c),
                    int(row.d),
                    float(row.log_odds_ratio).hex(),
                    float(row.p_value).hex(),
                    float(row.p_adjusted).hex(),
                    int(row.rank),
                    bool(row.positive_enrichment),
                    bool(row.passes_min_overlap),
                    bool(row.significant),
                    bool(row.candidate),
                    row.candidate_status,
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

    parameters = {
        "producer_node_id": "OpenBioSingleCellGeneSetOverrepresentation",
        "artifact_role": "generic_ora_evidence_table",
        "comparison": comparison,
        "comparison_column": comparison_column,
        "gene_column": gene_column,
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
        "min_overlap": min_overlap,
        "max_p_adjusted": max_p_adjusted,
        "max_output_rows": max_output_rows,
        "alternative": "greater",
        "n_bg": len(universe_genes),
        "ha_corr": 0.5,
        "multiple_testing": "Benjamini-Hochberg across every retained resource set",
        "complete_family_returned": True,
        "result_fingerprint": result_sha256(evidence),
    }
    warnings = [
        *resource.get("warnings", []),
        *runtime_provenance_warnings,
        "Resource organism, namespace, scope, license, and citation are caller declarations; SHA-256 identifies exact bytes but not biological suitability.",
        "Selection and ORA are not independent confirmatory tests; thresholds flag rows but never remove the tested family.",
    ]
    if evidence_scope == "cluster_marker_evidence":
        warnings.append(
            "This is exploratory Cluster marker evidence for annotation review, not Sample-level Condition inference or a Curated annotation."
        )
    if resource["accounting"]["duplicate_pairs_removed"]:
        warnings.append(
            f"Collapsed {resource['accounting']['duplicate_pairs_removed']:,} exact duplicate resource edges."
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
            "sources_removed_by_min_targets_preview": removed_by_min[:100],
            "sources_removed_by_min_targets_preview_truncated": len(removed_by_min) > 100,
            "set_accounting_preview": set_accounting[:100],
            "set_accounting_preview_truncated": len(set_accounting) > 100,
            "tested_sources_preview": retained_sources[:100],
            "tested_sources_preview_truncated": len(retained_sources) > 100,
        }
    )
    diagnostics = {
        "provenance": provenance_declarations,
        "selection": {
            "comparison": comparison,
            "evidence_scope": evidence_scope,
            "inference_unit": inference_unit,
            "direction": direction,
            "replicate_aware": replicate_aware,
            "technical_batch_handling": technical_batch_handling,
            "selected_count": len(selected_set),
            "universe_count": len(universe_genes),
            "selected_genes_preview": ordered_selected[:100],
            "selected_genes_preview_truncated": len(ordered_selected) > 100,
            "analysis_fingerprint_sha256": expected_analysis_fingerprint,
            "ranking_fingerprint_sha256": expected_ranking_fingerprint,
            "universe_fingerprint_sha256": expected_universe_fingerprint,
            "table_content_fingerprint_sha256": expected_table_content_fingerprint,
            "universe_content_fingerprint_sha256": expected_universe_content_fingerprint,
        },
        "resource": resource_accounting,
        "thresholds": {"min_overlap": min_overlap, "max_p_adjusted": max_p_adjusted},
        "tested_gene_sets": len(evidence),
        "positive_gene_sets": int(evidence["positive_enrichment"].sum()),
        "significant_gene_sets": int(evidence["significant"].sum()),
        "candidate_gene_sets": int(evidence["candidate"].sum()),
        "zero_overlap_gene_sets": int((evidence["overlap_count"] == 0).sum()),
        "top_results": evidence.head(10).to_dict("records"),
        "backend": {
            "api": "decoupler.mt.query_set",
            "version": str(version),
            "alternative": "greater",
            "n_bg": len(universe_genes),
            "ha_corr": 0.5,
            "tmin": min_targets,
            "scipy_version": scipy.__version__,
        },
        "parameters": parameters,
        "warnings": warnings,
        "software_versions": _standalone_generic_ora_software_versions(
            openbio_version=openbio_version,
            decoupler_version=str(version),
            scipy_version=scipy.__version__,
            pandas_version=pd.__version__,
            numpy_version=np.__version__,
        ),
    }
    json.dumps(diagnostics, ensure_ascii=False, allow_nan=False)
    return evidence, diagnostics


def run_generic_ora_evidence(table: Any, universe: Any, **kwargs: Any) -> tuple[DataFrame, dict[str, Any]]:
    return _standalone_run_generic_ora(table, universe, **kwargs)


def build_generic_ora_summary(diagnostics: dict[str, Any]) -> dict[str, Any]:
    return _standalone_generic_ora_summary(diagnostics)


def generic_ora_code(**parameters: Any) -> str:
    implementations = "\n\n".join(
        (
            pinned_enrichment_validation_code(),
            gene_set_resource_code(),
            textwrap.dedent(inspect.getsource(_standalone_generic_ora_software_versions)).strip(),
            textwrap.dedent(inspect.getsource(_standalone_generic_ora_summary)).strip(),
            textwrap.dedent(inspect.getsource(_standalone_run_generic_ora)).strip(),
        )
    )
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in parameters.items())
    return f"""from __future__ import annotations

{implementations}


def run_gene_set_overrepresentation(table, universe, decoupler_module=None):
    evidence, diagnostics = _standalone_run_generic_ora(
        table,
        universe,
        {rendered},
        decoupler_module=decoupler_module,
    )
    return evidence, _standalone_generic_ora_summary(diagnostics)
"""


__all__ = [
    "GENERIC_ORA_COLUMNS",
    "build_generic_ora_summary",
    "generic_ora_code",
    "run_generic_ora_evidence",
]
