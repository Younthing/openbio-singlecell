from __future__ import annotations

import inspect
import textwrap
from typing import TYPE_CHECKING, Any

from .dgidb_resource import (
    DGIdbResource,
    _standalone_prepare_dgidb_for_universe,
    _standalone_validate_portable_dgidb_resource,
    dgidb_portable_validation_code,
    validate_dgidb_resource,
)
from .enrichment_artifacts import (
    _standalone_validate_pinned_enrichment_pair,
    pinned_enrichment_validation_code,
)
from .generic_ora import (
    _standalone_generic_ora_software_versions,
    _standalone_run_generic_ora,
)
from .ranked_enrichment import (
    _standalone_ranked_gsea_software_versions,
    _standalone_run_ranked_gsea,
)

if TYPE_CHECKING:
    from pandas import DataFrame


DRUG_ORA_COLUMNS = [
    "scope",
    "comparison",
    "drug",
    "selected_count",
    "universe_count",
    "set_size_before_universe",
    "set_size_in_universe",
    "overlap_count",
    "overlap_genes",
    "resource_sources",
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

DRUG_GSEA_COLUMNS = [
    "comparison",
    "drug",
    "normalized_enrichment_score",
    "p_adjusted",
    "rank",
    "significant",
    "set_size_before_universe",
    "set_size_in_universe",
    "matched_targets",
    "resource_sources",
]


def _standalone_drug_evidence_fingerprint(frame, schema):
    import hashlib
    import json
    import math

    def plain(value):
        if value is None or isinstance(value, (str, bool, int)):
            return value
        if isinstance(value, float):
            if not math.isfinite(value):
                raise ValueError("Drug evidence fingerprint rejects non-finite values.")
            return value
        if hasattr(value, "item"):
            return plain(value.item())
        return str(value)

    payload = {
        "schema": schema,
        "columns": list(frame.columns),
        "rows": [[plain(value) for value in row] for row in frame.itertuples(index=False, name=None)],
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _standalone_resource_sources(resource_table):
    import json

    sources = {}
    for drug, group in resource_table.groupby("drug", sort=False):
        values = []
        for encoded in group["sources"].tolist():
            for item in json.loads(encoded):
                if item not in values:
                    values.append(item)
        sources[drug] = json.dumps(values, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    return sources


def _standalone_post_backend_pair_check(
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
    _standalone_validate_pinned_enrichment_pair(
        table_input,
        universe_input,
        artifact_family=artifact_family,
        purpose=purpose,
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


def _standalone_run_drug_ora(
    table_input,
    universe_input,
    resource_table,
    resource_metadata,
    resource_accounting,
    resource_artifact_metadata,
    *,
    expected_resource_artifact_fingerprint,
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
    provenance_warnings=(),
    provenance_declarations=None,
    min_targets=3,
    min_overlap=2,
    max_p_adjusted=0.05,
    max_output_rows=100000,
    openbio_version="not-installed",
    decoupler_module=None,
):
    import json

    _, universe_genes = _standalone_validate_pinned_enrichment_pair(
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
    resource_table, resource_metadata, resource_accounting, resource_artifact_metadata = (
        _standalone_validate_portable_dgidb_resource(
            resource_table,
            resource_metadata,
            resource_accounting,
            resource_artifact_metadata,
            copy_table=False,
        )
    )
    artifact_fingerprint = resource_artifact_metadata["artifact_fingerprint_sha256"]
    if artifact_fingerprint != expected_resource_artifact_fingerprint:
        raise ValueError("Drug Hypergeometric resource differs from the pinned upstream DGIdb artifact.")
    prepared = _standalone_prepare_dgidb_for_universe(
        resource_table,
        resource_metadata,
        resource_accounting,
        resource_artifact_metadata,
        universe_genes,
        min_targets,
        "Drug Hypergeometric",
    )
    typed_path = f"typed:{artifact_fingerprint}"
    evidence, diagnostics = _standalone_run_generic_ora(
        table_input,
        universe_input,
        artifact_family=artifact_family,
        comparison=comparison,
        comparison_column=comparison_column,
        gene_column=gene_column,
        evidence_scope=evidence_scope,
        inference_unit=inference_unit,
        direction=direction,
        replicate_aware=replicate_aware,
        technical_batch_handling=technical_batch_handling,
        upstream_parameters=upstream_parameters,
        provenance_warnings=provenance_warnings,
        provenance_declarations=provenance_declarations,
        expected_analysis_fingerprint=expected_analysis_fingerprint,
        expected_ranking_fingerprint=expected_ranking_fingerprint,
        expected_universe_fingerprint=expected_universe_fingerprint,
        expected_table_content_fingerprint=expected_table_content_fingerprint,
        expected_universe_content_fingerprint=expected_universe_content_fingerprint,
        resource_path=typed_path,
        requested_resource_path=typed_path,
        resource_metadata_json="{}",
        expected_resource_sha256=resource_artifact_metadata["raw_file_sha256"],
        preloaded_resource=prepared,
        source_column="source",
        target_column="target",
        min_targets=min_targets,
        min_overlap=min_overlap,
        max_p_adjusted=max_p_adjusted,
        max_output_rows=max_output_rows,
        openbio_version=openbio_version,
        decoupler_module=decoupler_module,
    )
    _standalone_post_backend_pair_check(
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
    _, _, _, resource_after = _standalone_validate_portable_dgidb_resource(
        resource_table,
        resource_metadata,
        resource_accounting,
        resource_artifact_metadata,
        copy_table=False,
    )
    if resource_after["artifact_fingerprint_sha256"] != artifact_fingerprint:
        raise RuntimeError("Drug Hypergeometric DGIdb artifact changed during backend execution.")

    sources = _standalone_resource_sources(resource_table)
    evidence = evidence.rename(columns={"gene_set": "drug"})
    evidence.insert(
        evidence.columns.get_loc("a"),
        "resource_sources",
        evidence["drug"].map(sources),
    )
    evidence = evidence.loc[:, [
        "scope",
        "comparison",
        "drug",
        "selected_count",
        "universe_count",
        "set_size_before_universe",
        "set_size_in_universe",
        "overlap_count",
        "overlap_genes",
        "resource_sources",
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
    ]].copy()
    result_fingerprint = _standalone_drug_evidence_fingerprint(evidence, "openbio-drug-ora/v1")
    parameters = dict(diagnostics["parameters"])
    parameters.update(
        {
            "producer_node_id": "OpenBioSingleCellDrugHypergeometric",
            "artifact_role": "drug_ora_evidence_table",
            "resource_artifact_fingerprint_sha256": artifact_fingerprint,
            "resource_canonical_content_fingerprint_sha256": resource_artifact_metadata[
                "canonical_content_fingerprint_sha256"
            ],
            "resource_metadata": resource_metadata,
            "resource_type": "OPENBIO_DGIDB_RESOURCE",
            "result_fingerprint": result_fingerprint,
        }
    )
    resource_diagnostics = dict(diagnostics["resource"])
    resource_diagnostics.update(
        {
            "metadata": resource_metadata,
            "artifact_fingerprint_sha256": artifact_fingerprint,
            "canonical_content_fingerprint_sha256": resource_artifact_metadata[
                "canonical_content_fingerprint_sha256"
            ],
            "source_license_review": resource_metadata["source_license_review"],
            "download_url": resource_metadata["download_url"],
        }
    )
    warnings = [
        warning
        for warning in diagnostics["warnings"]
        if not warning.startswith("Resource organism, namespace, scope, license")
    ]
    warnings.extend(
        [
            "The typed DGIdb artifact proves producer schema and content fingerprints; source-license review and biological suitability remain caller attestations.",
            "Drug-target overrepresentation is not therapeutic benefit, disease-signature reversal, efficacy, dose, safety, or a treatment recommendation.",
        ]
    )
    diagnostics.update(
        {
            "resource": resource_diagnostics,
            "tested_drugs": len(evidence),
            "positive_drugs": int(evidence["positive_enrichment"].sum()),
            "significant_drugs": int(evidence["significant"].sum()),
            "candidate_drugs": int(evidence["candidate"].sum()),
            "zero_overlap_drugs": int((evidence["overlap_count"] == 0).sum()),
            "top_drug_results": evidence.head(10).to_dict("records"),
            "parameters": parameters,
            "warnings": warnings,
            "result_fingerprint_sha256": result_fingerprint,
        }
    )
    json.dumps(diagnostics, ensure_ascii=False, allow_nan=False)
    return evidence, diagnostics


def _standalone_build_drug_ora_summary(diagnostics):
    import json
    import math
    from collections.abc import Mapping

    def plain(value):
        if value is None or isinstance(value, (str, bool, int)):
            return value
        if isinstance(value, float):
            return value if math.isfinite(value) else None
        if hasattr(value, "item"):
            return plain(value.item())
        if isinstance(value, Mapping):
            return {str(key): plain(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
            sequence = value.tolist() if hasattr(value, "tolist") else value
            return [plain(item) for item in sequence]
        return str(value)

    required = {
        "provenance",
        "selection",
        "resource",
        "thresholds",
        "tested_drugs",
        "positive_drugs",
        "significant_drugs",
        "candidate_drugs",
        "zero_overlap_drugs",
        "top_drug_results",
        "backend",
        "parameters",
        "warnings",
        "software_versions",
        "result_fingerprint_sha256",
    }
    if not isinstance(diagnostics, Mapping) or not required.issubset(diagnostics):
        missing = sorted(required - set(diagnostics)) if isinstance(diagnostics, Mapping) else sorted(required)
        raise ValueError(f"Drug Hypergeometric summary diagnostics are missing fields: {missing}.")
    key_results = {
        key: plain(diagnostics[key])
        for key in required
        if key not in {"parameters", "warnings", "software_versions"}
    }
    selection = key_results["selection"]
    resource = key_results["resource"]
    metadata = resource["metadata"]
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellDrugHypergeometric",
        "methods": (
            "One explicit selected-gene table whose structure and current-content fingerprints were validated, "
            "together with its fingerprint-paired exact tested universe, was "
            "tested against every eligible drug target set in the pinned DGIdb artifact. decoupler 2.2 mt.query_set "
            "used a one-sided greater alternative, exact background size, Haldane 0.5 correction, and min_targets "
            "after universe intersection; Fisher tables, log odds, raw p-values, and full-family Benjamini-Hochberg "
            "adjustment were independently reconstructed and verified."
        ),
        "results": (
            f"For {selection['comparison']!r}, tested {int(key_results['tested_drugs']):,} drugs from "
            f"DGIdb {metadata['version']!r}; {int(key_results['significant_drugs']):,} had adjusted p-values "
            f"within the declared threshold and {int(key_results['candidate_drugs']):,} also had positive log "
            "odds and met the declared overlap threshold."
        ),
        "key_results": key_results,
        "parameters": plain(diagnostics["parameters"]),
        "warnings": [str(value) for value in diagnostics["warnings"]],
        "limitations": [
            "Overrepresentation means selected genes overlap DGIdb target claims more than expected under the exact universe; it does not imply therapeutic benefit, reversal direction, dose, safety, efficacy, or recommendation.",
            "Selection and ORA are not independent confirmatory tests, overlapping drug target sets create dependent hypotheses, and BH adjustment does not make drug labels independent truths.",
            "Producer identity, approval, scope, inference unit, replicate awareness, and Technical-batch handling are caller declarations; OpenBio does not validate upstream scientific authority or formal Condition interpretation.",
            "DGIdb source coverage, evidence heterogeneity, HGNC matching, and source-specific licenses limit interpretation.",
        ],
        "references": [
            {"citation": "Fisher RA. On the interpretation of chi-square from contingency tables. JRSS. 1922;85:87-94.", "url": "https://doi.org/10.2307/2340521", "doi": "10.2307/2340521", "kind": "method"},
            {"citation": "Haldane JBS. The estimation and significance of the logarithm of a ratio of frequencies.", "url": "https://doi.org/10.1111/j.1469-1809.1955.tb01285.x", "doi": "10.1111/j.1469-1809.1955.tb01285.x", "kind": "method"},
            {"citation": "Benjamini Y, Hochberg Y. Controlling the false discovery rate.", "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x", "doi": "10.1111/j.2517-6161.1995.tb02031.x", "kind": "method"},
            {"citation": "Cannon M, et al. DGIdb 5.0. Nucleic Acids Research. 2024;52:D1227-D1235.", "url": "https://doi.org/10.1093/nar/gkad1040", "doi": "10.1093/nar/gkad1040", "kind": "resource"},
            {"citation": "Heumos L, et al. pertpy: an end-to-end framework for perturbation analysis. Nature Methods. 2025.", "url": "https://doi.org/10.1038/s41592-025-02909-7", "doi": "10.1038/s41592-025-02909-7", "kind": "software"},
            {"citation": "decoupler mt.query_set official documentation.", "url": "https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.query_set.html", "doi": None, "kind": "software_documentation"},
            {"citation": f"{metadata['name']} {metadata['version']}: {metadata['citation']}", "url": metadata["download_url"], "doi": None, "kind": "resource_snapshot"},
        ],
        "software_versions": plain(diagnostics["software_versions"]),
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def _standalone_run_drug_gsea(
    table_input,
    universe_input,
    resource_table,
    resource_metadata,
    resource_accounting,
    resource_artifact_metadata,
    *,
    expected_resource_artifact_fingerprint,
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
    provenance_warnings=(),
    provenance_declarations=None,
    min_targets=15,
    max_targets=500,
    n_permutations=1000,
    random_seed=1,
    max_output_rows=100000,
    openbio_version="not-installed",
    decoupler_module=None,
):
    import json

    _, universe_genes = _standalone_validate_pinned_enrichment_pair(
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
    resource_table, resource_metadata, resource_accounting, resource_artifact_metadata = (
        _standalone_validate_portable_dgidb_resource(
            resource_table,
            resource_metadata,
            resource_accounting,
            resource_artifact_metadata,
            copy_table=False,
        )
    )
    artifact_fingerprint = resource_artifact_metadata["artifact_fingerprint_sha256"]
    if artifact_fingerprint != expected_resource_artifact_fingerprint:
        raise ValueError("Drug GSEA resource differs from the pinned upstream DGIdb artifact.")
    prepared = _standalone_prepare_dgidb_for_universe(
        resource_table,
        resource_metadata,
        resource_accounting,
        resource_artifact_metadata,
        universe_genes,
        min_targets,
        "Drug GSEA",
    )
    typed_path = f"typed:{artifact_fingerprint}"
    evidence, diagnostics = _standalone_run_ranked_gsea(
        table_input,
        universe_input,
        artifact_family=artifact_family,
        comparison=comparison,
        comparison_column=comparison_column,
        gene_column=gene_column,
        score_column=score_column,
        evidence_scope=evidence_scope,
        inference_unit=inference_unit,
        direction=direction,
        replicate_aware=replicate_aware,
        technical_batch_handling=technical_batch_handling,
        upstream_parameters=upstream_parameters,
        provenance_warnings=provenance_warnings,
        provenance_declarations=provenance_declarations,
        expected_analysis_fingerprint=expected_analysis_fingerprint,
        expected_ranking_fingerprint=expected_ranking_fingerprint,
        expected_universe_fingerprint=expected_universe_fingerprint,
        expected_table_content_fingerprint=expected_table_content_fingerprint,
        expected_universe_content_fingerprint=expected_universe_content_fingerprint,
        resource_path=typed_path,
        requested_resource_path=typed_path,
        resource_metadata_json="{}",
        expected_resource_sha256=resource_artifact_metadata["raw_file_sha256"],
        preloaded_resource=prepared,
        source_column="source",
        target_column="target",
        min_targets=min_targets,
        max_targets=max_targets,
        n_permutations=n_permutations,
        random_seed=random_seed,
        max_output_rows=max_output_rows,
        openbio_version=openbio_version,
        decoupler_module=decoupler_module,
    )
    _standalone_post_backend_pair_check(
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
    _, _, _, resource_after = _standalone_validate_portable_dgidb_resource(
        resource_table,
        resource_metadata,
        resource_accounting,
        resource_artifact_metadata,
        copy_table=False,
    )
    if resource_after["artifact_fingerprint_sha256"] != artifact_fingerprint:
        raise RuntimeError("Drug GSEA DGIdb artifact changed during backend execution.")

    sources = _standalone_resource_sources(resource_table)
    evidence = evidence.rename(columns={"gene_set": "drug"})
    matched = {
        drug: json.dumps(
            list(prepared["targets_in_universe"][drug]),
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
        )
        for drug in evidence["drug"].tolist()
    }
    evidence["matched_targets"] = evidence["drug"].map(matched)
    evidence["resource_sources"] = evidence["drug"].map(sources)
    evidence = evidence.loc[:, [
        "comparison",
        "drug",
        "normalized_enrichment_score",
        "p_adjusted",
        "rank",
        "significant",
        "set_size_before_universe",
        "set_size_in_universe",
        "matched_targets",
        "resource_sources",
    ]].copy()
    result_fingerprint = _standalone_drug_evidence_fingerprint(evidence, "openbio-drug-gsea/v1")
    parameters = dict(diagnostics["parameters"])
    parameters.update(
        {
            "producer_node_id": "OpenBioSingleCellDrugGSEA",
            "artifact_role": "drug_gsea_evidence_table",
            "resource_artifact_fingerprint_sha256": artifact_fingerprint,
            "resource_canonical_content_fingerprint_sha256": resource_artifact_metadata[
                "canonical_content_fingerprint_sha256"
            ],
            "resource_metadata": resource_metadata,
            "resource_type": "OPENBIO_DGIDB_RESOURCE",
            "result_fingerprint": result_fingerprint,
        }
    )
    resource_diagnostics = dict(diagnostics["resource"])
    resource_diagnostics.update(
        {
            "metadata": resource_metadata,
            "artifact_fingerprint_sha256": artifact_fingerprint,
            "canonical_content_fingerprint_sha256": resource_artifact_metadata[
                "canonical_content_fingerprint_sha256"
            ],
            "source_license_review": resource_metadata["source_license_review"],
            "download_url": resource_metadata["download_url"],
        }
    )
    top_positive = (
        evidence.loc[evidence["normalized_enrichment_score"] > 0]
        .sort_values(
            ["p_adjusted", "normalized_enrichment_score", "drug"],
            ascending=[True, False, True],
            kind="mergesort",
        )
        .head(10)
        .to_dict("records")
    )
    top_negative = (
        evidence.loc[evidence["normalized_enrichment_score"] < 0]
        .sort_values(
            ["p_adjusted", "normalized_enrichment_score", "drug"],
            ascending=[True, True, True],
            kind="mergesort",
        )
        .head(10)
        .to_dict("records")
    )
    warnings = [
        warning
        for warning in diagnostics["warnings"]
        if not warning.startswith("Resource organism, namespace, scope, license")
    ]
    warnings.extend(
        [
            "The typed DGIdb artifact proves producer schema and content fingerprints; source-license review and biological suitability remain caller attestations.",
            "Positive or negative drug-target enrichment is not a connectivity-map signature-reversal result and conveys no efficacy, dose, safety, or recommendation.",
        ]
    )
    diagnostics.update(
        {
            "resource": resource_diagnostics,
            "tested_drugs": len(evidence),
            "significant_drugs": int(evidence["significant"].sum()),
            "positive_drugs": int((evidence["normalized_enrichment_score"] > 0).sum()),
            "negative_drugs": int((evidence["normalized_enrichment_score"] < 0).sum()),
            "zero_nes_drugs": int((evidence["normalized_enrichment_score"] == 0).sum()),
            "top_positive_drugs": top_positive,
            "top_negative_drugs": top_negative,
            "parameters": parameters,
            "warnings": warnings,
            "result_fingerprint_sha256": result_fingerprint,
        }
    )
    json.dumps(diagnostics, ensure_ascii=False, allow_nan=False)
    return evidence, diagnostics


def _standalone_build_drug_gsea_summary(diagnostics):
    import json
    import math
    from collections.abc import Mapping

    def plain(value):
        if value is None or isinstance(value, (str, bool, int)):
            return value
        if isinstance(value, float):
            return value if math.isfinite(value) else None
        if hasattr(value, "item"):
            return plain(value.item())
        if isinstance(value, Mapping):
            return {str(key): plain(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
            sequence = value.tolist() if hasattr(value, "tolist") else value
            return [plain(item) for item in sequence]
        return str(value)

    required = {
        "provenance",
        "ranking",
        "resource",
        "tested_drugs",
        "significant_drugs",
        "positive_drugs",
        "negative_drugs",
        "zero_nes_drugs",
        "top_positive_drugs",
        "top_negative_drugs",
        "permutations",
        "random_seed",
        "permutation_null_status",
        "formal_permutation_inference_valid",
        "backend",
        "parameters",
        "warnings",
        "software_versions",
        "result_fingerprint_sha256",
    }
    if not isinstance(diagnostics, Mapping) or not required.issubset(diagnostics):
        missing = sorted(required - set(diagnostics)) if isinstance(diagnostics, Mapping) else sorted(required)
        raise ValueError(f"Drug GSEA summary diagnostics are missing fields: {missing}.")
    key_results = {
        key: plain(diagnostics[key])
        for key in required
        if key not in {"parameters", "warnings", "software_versions"}
    }
    ranking = key_results["ranking"]
    resource = key_results["resource"]
    metadata = resource["metadata"]
    formal_permutation_inference_valid = bool(key_results["formal_permutation_inference_valid"])
    limitations = [
        "Positive or negative enrichment only locates DGIdb target claims along the supplied statistic; it is not a connectivity-map signature-reversal analysis and implies no therapeutic direction, efficacy, dose, safety, or recommendation.",
        "decoupler mt.gsea does not expose raw p-values or leading-edge genes through this interface; matched_targets reports resource coverage, not a leading edge.",
        "Producer identity, approval, scope, inference unit, replicate awareness, and Technical-batch handling are caller declarations; OpenBio does not validate upstream scientific authority or formal Condition interpretation.",
        "DGIdb source coverage, evidence heterogeneity, HGNC matching, overlapping target sets, and source-specific licenses limit interpretation.",
    ]
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
        "node_id": "OpenBioSingleCellDrugGSEA",
        "methods": (
            "One complete continuous ranking whose structure and current-content fingerprints were validated, "
            "together with its fingerprint-paired exact tested universe, was tested "
            "against every eligible target set in the pinned DGIdb artifact using decoupler 2.2 mt.gsea with "
            f"raw=False, empty=False, fixed exact rank orientation, caller seed, and {int(key_results['permutations']):,} requested permutations. "
            "The public backend returns normalized enrichment scores and full-family Benjamini-Hochberg-adjusted "
            "empirical p-values; raw p-values and leading-edge genes are unavailable."
        ),
        "results": (
            f"For {ranking['comparison']!r}, tested {int(key_results['tested_drugs']):,} drugs from DGIdb "
            f"{metadata['version']!r}; {int(key_results['significant_drugs']):,} had backend-returned adjusted values <= 0.05, "
            f"with {int(key_results['positive_drugs']):,} positive and {int(key_results['negative_drugs']):,} "
            "negative NES values."
            + (
                " Because seed 0 disables permutation shuffling in decoupler 2.2, adjusted-value threshold "
                "counts are computational flags and not valid formal permutation evidence."
                if not formal_permutation_inference_valid
                else ""
            )
        ),
        "key_results": key_results,
        "parameters": plain(diagnostics["parameters"]),
        "warnings": [str(value) for value in diagnostics["warnings"]],
        "limitations": limitations,
        "references": [
            {"citation": "Subramanian A, et al. Gene set enrichment analysis. PNAS. 2005;102:15545-15550.", "url": "https://doi.org/10.1073/pnas.0506580102", "doi": "10.1073/pnas.0506580102", "kind": "method"},
            {"citation": "Badia-i-Mompel P, et al. decoupleR. Bioinformatics Advances. 2022;2:vbac016.", "url": "https://doi.org/10.1093/bioadv/vbac016", "doi": "10.1093/bioadv/vbac016", "kind": "software"},
            {"citation": "Benjamini Y, Hochberg Y. Controlling the false discovery rate.", "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x", "doi": "10.1111/j.2517-6161.1995.tb02031.x", "kind": "method"},
            {"citation": "Cannon M, et al. DGIdb 5.0. Nucleic Acids Research. 2024;52:D1227-D1235.", "url": "https://doi.org/10.1093/nar/gkad1040", "doi": "10.1093/nar/gkad1040", "kind": "resource"},
            {"citation": "Heumos L, et al. pertpy: an end-to-end framework for perturbation analysis. Nature Methods. 2025.", "url": "https://doi.org/10.1038/s41592-025-02909-7", "doi": "10.1038/s41592-025-02909-7", "kind": "software"},
            {"citation": "decoupler mt.gsea official documentation.", "url": "https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.gsea.html", "doi": None, "kind": "software_documentation"},
            {"citation": f"{metadata['name']} {metadata['version']}: {metadata['citation']}", "url": metadata["download_url"], "doi": None, "kind": "resource_snapshot"},
        ],
        "software_versions": plain(diagnostics["software_versions"]),
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def run_drug_ora(
    table: Any,
    universe: Any,
    resource: DGIdbResource,
    *,
    copy_resource: bool = True,
    **parameters: Any,
) -> tuple[DataFrame, dict[str, Any]]:
    resource_table, metadata, accounting, artifact = validate_dgidb_resource(
        resource, copy_payload=copy_resource
    )
    evidence, diagnostics = _standalone_run_drug_ora(
        table,
        universe,
        resource_table,
        metadata,
        accounting,
        artifact,
        expected_resource_artifact_fingerprint=artifact["artifact_fingerprint_sha256"],
        **parameters,
    )
    return evidence, _standalone_build_drug_ora_summary(diagnostics)


def run_drug_gsea(
    table: Any,
    universe: Any,
    resource: DGIdbResource,
    *,
    copy_resource: bool = True,
    **parameters: Any,
) -> tuple[DataFrame, dict[str, Any]]:
    resource_table, metadata, accounting, artifact = validate_dgidb_resource(
        resource, copy_payload=copy_resource
    )
    evidence, diagnostics = _standalone_run_drug_gsea(
        table,
        universe,
        resource_table,
        metadata,
        accounting,
        artifact,
        expected_resource_artifact_fingerprint=artifact["artifact_fingerprint_sha256"],
        **parameters,
    )
    return evidence, _standalone_build_drug_gsea_summary(diagnostics)


def _drug_enrichment_implementations(*, kind: str) -> str:
    common = [
        pinned_enrichment_validation_code(),
        dgidb_portable_validation_code(include_prepare=True),
        textwrap.dedent(inspect.getsource(_standalone_drug_evidence_fingerprint)).strip(),
        textwrap.dedent(inspect.getsource(_standalone_resource_sources)).strip(),
        textwrap.dedent(inspect.getsource(_standalone_post_backend_pair_check)).strip(),
    ]
    if kind == "ora":
        common.extend(
            [
                textwrap.dedent(inspect.getsource(_standalone_generic_ora_software_versions)).strip(),
                textwrap.dedent(inspect.getsource(_standalone_run_generic_ora)).strip(),
                textwrap.dedent(inspect.getsource(_standalone_run_drug_ora)).strip(),
                textwrap.dedent(inspect.getsource(_standalone_build_drug_ora_summary)).strip(),
            ]
        )
    elif kind == "gsea":
        common.extend(
            [
                textwrap.dedent(inspect.getsource(_standalone_ranked_gsea_software_versions)).strip(),
                textwrap.dedent(inspect.getsource(_standalone_run_ranked_gsea)).strip(),
                textwrap.dedent(inspect.getsource(_standalone_run_drug_gsea)).strip(),
                textwrap.dedent(inspect.getsource(_standalone_build_drug_gsea_summary)).strip(),
            ]
        )
    else:
        raise ValueError(f"Unsupported drug enrichment code kind: {kind!r}.")
    return "\n\n".join(common)


def drug_ora_code(
    *,
    resource_metadata: dict[str, Any],
    resource_accounting: dict[str, Any],
    resource_artifact_metadata: dict[str, Any],
    **parameters: Any,
) -> str:
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in parameters.items())
    return f'''from __future__ import annotations

{_drug_enrichment_implementations(kind="ora")}


def run_drug_hypergeometric(table, universe, resource_table, decoupler_module=None):
    evidence, diagnostics = _standalone_run_drug_ora(
        table,
        universe,
        resource_table,
        {resource_metadata!r},
        {resource_accounting!r},
        {resource_artifact_metadata!r},
        expected_resource_artifact_fingerprint={resource_artifact_metadata["artifact_fingerprint_sha256"]!r},
        {rendered},
        decoupler_module=decoupler_module,
    )
    return evidence, _standalone_build_drug_ora_summary(diagnostics)
'''


def drug_gsea_code(
    *,
    resource_metadata: dict[str, Any],
    resource_accounting: dict[str, Any],
    resource_artifact_metadata: dict[str, Any],
    **parameters: Any,
) -> str:
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in parameters.items())
    return f'''from __future__ import annotations

{_drug_enrichment_implementations(kind="gsea")}


def run_drug_gsea(table, universe, resource_table, decoupler_module=None):
    evidence, diagnostics = _standalone_run_drug_gsea(
        table,
        universe,
        resource_table,
        {resource_metadata!r},
        {resource_accounting!r},
        {resource_artifact_metadata!r},
        expected_resource_artifact_fingerprint={resource_artifact_metadata["artifact_fingerprint_sha256"]!r},
        {rendered},
        decoupler_module=decoupler_module,
    )
    return evidence, _standalone_build_drug_gsea_summary(diagnostics)
'''


__all__ = [
    "DRUG_GSEA_COLUMNS",
    "DRUG_ORA_COLUMNS",
    "drug_gsea_code",
    "drug_ora_code",
    "run_drug_gsea",
    "run_drug_ora",
]
