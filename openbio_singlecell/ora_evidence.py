from __future__ import annotations

import inspect
import textwrap
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pandas import DataFrame


ORA_EVIDENCE_COLUMNS = [
    "group",
    "source",
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
    "p_adj_within_group",
    "p_adj_global",
    "rank_within_group",
    "positive_enrichment",
    "passes_min_overlap",
    "significant_within_group",
    "candidate_status",
]


def marker_ora_table_fingerprint(table: Any) -> str:
    import hashlib
    import json

    import pandas as pd

    if not isinstance(table, pd.DataFrame) or list(table.columns) != ORA_EVIDENCE_COLUMNS:
        raise ValueError("Marker ORA table fingerprint requires the exact canonical columns.")
    integer_columns = {
        "selected_count",
        "universe_count",
        "set_size_before_universe",
        "set_size_in_universe",
        "overlap_count",
        "a",
        "b",
        "c",
        "d",
        "rank_within_group",
    }
    float_columns = {"log_odds_ratio", "p_value", "p_adj_within_group", "p_adj_global"}
    boolean_columns = {
        "positive_enrichment",
        "passes_min_overlap",
        "significant_within_group",
    }
    rows = []
    for row in table.itertuples(index=False, name=None):
        normalized = []
        for column, value in zip(table.columns, row, strict=True):
            if column in boolean_columns:
                normalized.append(bool(value))
            elif column in integer_columns:
                normalized.append(int(value))
            elif column in float_columns:
                normalized.append(float(value).hex())
            else:
                normalized.append(str(value))
        rows.append(normalized)
    payload = json.dumps(
        {"schema": "openbio-singlecell/marker-ora-evidence/v1", "rows": rows},
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _standalone_marker_ora_software_versions(
    *,
    openbio_version,
    decoupler_version,
    scipy_version,
    pandas_version,
    numpy_version,
):
    """Collect the software-version map disclosed by standalone ORA reports."""
    import platform
    from importlib import metadata as importlib_metadata

    versions = {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
        "decoupler": str(decoupler_version),
        "scipy": str(scipy_version),
        "pandas": str(pandas_version),
        "numpy": str(numpy_version),
    }
    for package in ("scanpy", "anndata"):
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _standalone_marker_ora_summary(diagnostics):
    """Build the strict full Marker ORA scientific summary without side effects."""
    import json
    import math
    from collections.abc import Mapping

    if not isinstance(diagnostics, Mapping):
        raise TypeError("Marker ORA summary diagnostics must be a mapping.")
    required = {"parameters", "warnings", "software_versions", "resource", "thresholds"}
    missing = sorted(required - set(diagnostics))
    if missing:
        raise ValueError(f"Marker ORA summary diagnostics are missing required fields: {missing}.")
    if not isinstance(diagnostics["parameters"], Mapping):
        raise TypeError("Marker ORA summary parameters must be a mapping.")
    if not isinstance(diagnostics["software_versions"], Mapping):
        raise TypeError("Marker ORA summary software_versions must be a mapping.")

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
    warnings = [str(value) for value in diagnostics["warnings"]]
    key_results = {
        key: plain_json(value)
        for key, value in diagnostics.items()
        if key not in {"parameters", "warnings", "software_versions"}
    }
    thresholds = key_results["thresholds"]
    methods = (
        "For each explicit positive Cluster marker set, decoupler 2.x mt.query_set tested one-sided "
        "over-representation against resource sets intersected with the exact tested-gene universe. "
        "The analysis fixed Fisher's alternative='greater', n_bg to the universe size, and the "
        "Haldane-Anscombe correction to 0.5; contingencies and raw p-values were independently "
        "reconstructed, with Benjamini-Hochberg adjustment within each group and globally."
    )
    results = (
        f"Tested {int(key_results['tested_pairs']):,} group-resource pairs for "
        f"{int(key_results['group_count']):,} groups; {int(key_results['eligible_rows']):,} rows met the "
        f"descriptive positive, overlap >= {int(thresholds['min_overlap'])}, and within-group adjusted "
        f"p <= {float(thresholds['max_p_adjusted']):g} gates. "
        f"{int(key_results['unresolved_group_count']):,} groups were unresolved, including "
        f"{int(key_results['groups_without_selected_markers_count']):,} with no selected marker after "
        "upstream filtering. No annotation label was committed."
    )
    resource = key_results["resource"]
    resource_metadata = resource["metadata"]
    references = [
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
                f"{resource_metadata['name']} {resource_metadata['version']} "
                f"({resource_metadata['date']}): {resource_metadata['citation']}"
            ),
            "url": f"urn:sha256:{resource['sha256']}",
            "kind": "practice",
            "doi": None,
        },
    ]
    limitations = [
        "Enrichment depends on marker selection, clustering resolution, tested universe, identifier compatibility, and resource coverage.",
        "Overlapping resource sets and data-dependent marker selection limit confirmatory interpretation of adjusted p-values.",
        "A top ORA-supported candidate is not a probability of cell identity, ground truth, or Curated annotation.",
        "One group-level candidate hides within-cluster heterogeneity and novel or mixed cell states absent from the resource.",
        "Cluster marker evidence is not Sample-level Condition inference; cells are not independent biological replicates.",
    ]
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellMarkerORAEvidence",
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": limitations,
        "references": references,
        "software_versions": plain_json(diagnostics["software_versions"]),
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def _standalone_marker_ora_evidence(
    marker_input,
    universe_input,
    *,
    expected_analysis_fingerprint,
    expected_ranking_fingerprint,
    expected_universe_fingerprint,
    expected_marker_content_fingerprint,
    expected_universe_content_fingerprint,
    group_labels,
    upstream_parameters,
    resource_path,
    requested_resource_path=None,
    resource_metadata_json,
    expected_resource_sha256=None,
    source_column="source",
    target_column="target",
    min_targets=3,
    min_overlap=2,
    max_p_adjusted=0.05,
    openbio_version="not-installed",
    decoupler_module=None,
):
    """Run explicit-set ORA with independent contingency and p-value verification."""
    import csv
    import hashlib
    import importlib
    import inspect as runtime_inspect
    import json
    import math
    import os
    from collections.abc import Mapping
    from importlib import metadata as importlib_metadata

    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats

    marker_columns = [
        "group",
        "gene",
        "rank",
        "score",
        "log2_fold_change_approx",
        "p_value",
        "p_adjusted",
        "fraction_in_group",
        "fraction_reference",
    ]
    universe_columns = ["gene", "universe_rank"]
    evidence_columns = [
        "group",
        "source",
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
        "p_adj_within_group",
        "p_adj_global",
        "rank_within_group",
        "positive_enrichment",
        "passes_min_overlap",
        "significant_within_group",
        "candidate_status",
    ]
    operation = "Marker ORA Evidence"

    def sha256_json(payload):
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

    def unwrap_artifact(value, *, name):
        if isinstance(value, pd.DataFrame):
            return value, None, None, None
        frame = getattr(value, "table", None)
        parameters = getattr(value, "parameters", None)
        source = getattr(value, "source", None)
        input_genes = getattr(value, "input_genes", None)
        if not isinstance(frame, pd.DataFrame):
            raise TypeError(f"{operation} {name} must be a pandas DataFrame or table artifact.")
        if not isinstance(parameters, Mapping) or not isinstance(source, Mapping):
            raise ValueError(f"{operation} {name} has malformed provenance.")
        source_parameters = source.get("parameters")
        if not isinstance(source_parameters, Mapping) or dict(source_parameters) != dict(parameters):
            raise ValueError(f"{operation} {name} parameters disagree with its immutable source.")
        return frame, dict(parameters), dict(source), input_genes

    marker_frame, marker_parameters, marker_source, marker_input_genes = unwrap_artifact(
        marker_input, name="marker table"
    )
    universe_frame, universe_parameters, universe_source, universe_input_genes = unwrap_artifact(
        universe_input, name="tested-gene universe"
    )
    if list(marker_frame.columns) != marker_columns:
        raise ValueError(f"{operation} marker table must use the exact canonical columns in order.")
    if list(universe_frame.columns) != universe_columns:
        raise ValueError(f"{operation} universe must use the exact canonical columns in order.")

    for column in ("group", "gene"):
        for value in marker_frame[column].tolist():
            if not isinstance(value, str):
                raise TypeError(f"{operation} marker column {column!r} must contain strings only.")
            if not value.strip() or value != value.strip():
                raise ValueError(
                    f"{operation} marker column {column!r} requires nonblank identifiers without surrounding whitespace."
                )
    if bool(marker_frame.duplicated(subset=["group", "gene"]).any()):
        raise ValueError(f"{operation} marker table contains duplicate (group, gene) rows.")
    marker_numeric_columns = marker_columns[2:]
    for column in marker_numeric_columns:
        series = marker_frame[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} marker column {column!r} must have a numeric dtype.")
        values = series.to_numpy(dtype=float)
        if values.size and not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} marker column {column!r} must contain finite values.")
    marker_ranks = marker_frame["rank"].to_numpy(dtype=float)
    if marker_ranks.size and (
        not bool(np.equal(marker_ranks, np.floor(marker_ranks)).all()) or bool((marker_ranks < 1).any())
    ):
        raise ValueError(f"{operation} marker ranks must be positive integers.")
    for column in ("p_value", "p_adjusted", "fraction_in_group", "fraction_reference"):
        values = marker_frame[column].to_numpy(dtype=float)
        if values.size and bool(((values < 0.0) | (values > 1.0)).any()):
            raise ValueError(f"{operation} marker column {column!r} must lie in [0, 1].")
    seen_group_blocks = set()
    active_group = None
    previous_rank = 0
    for group, rank in zip(marker_frame["group"].tolist(), marker_ranks, strict=True):
        if group != active_group:
            if group in seen_group_blocks:
                raise ValueError(f"{operation} marker rows for each group must form one contiguous block.")
            seen_group_blocks.add(group)
            active_group = group
            previous_rank = 0
        integer_rank = int(rank)
        if integer_rank <= previous_rank:
            raise ValueError(f"{operation} marker ranks must be strictly increasing within each group.")
        previous_rank = integer_rank

    def valid_sha256(value):
        return (
            isinstance(value, str) and len(value) == 64 and all(character in "0123456789abcdef" for character in value)
        )

    expected_fingerprints = {
        "analysis": expected_analysis_fingerprint,
        "ranking": expected_ranking_fingerprint,
        "universe": expected_universe_fingerprint,
        "marker content": expected_marker_content_fingerprint,
        "universe content": expected_universe_content_fingerprint,
    }
    for description, fingerprint in expected_fingerprints.items():
        if not valid_sha256(fingerprint):
            raise ValueError(f"{operation} expected {description} fingerprint is invalid.")

    def frame_content_fingerprint(frame, *, artifact, numeric_columns):
        encoder = json.JSONEncoder(
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        digest = hashlib.sha256()

        def update(value):
            for chunk in encoder.iterencode(value):
                digest.update(chunk.encode("utf-8"))

        digest.update(b'{"columns":')
        update(list(frame.columns))
        digest.update(b',"rows":[')
        first_row = True
        for row in frame.itertuples(index=False, name=None):
            encoded_row = []
            for column, value in zip(frame.columns, row, strict=True):
                encoded_row.append(float(value).hex() if column in numeric_columns else value)
            if not first_row:
                digest.update(b",")
            update(encoded_row)
            first_row = False
        digest.update(b'],"schema":')
        update(f"openbio-singlecell/{artifact}/v2")
        digest.update(b"}")
        return digest.hexdigest()

    observed_marker_content_fingerprint = frame_content_fingerprint(
        marker_frame,
        artifact="marker-table",
        numeric_columns=set(marker_numeric_columns),
    )
    if observed_marker_content_fingerprint != expected_marker_content_fingerprint:
        raise ValueError(f"{operation} marker table content does not match its pinned SHA-256 fingerprint.")
    observed_universe_content_fingerprint = frame_content_fingerprint(
        universe_frame,
        artifact="tested-gene-universe",
        numeric_columns={"universe_rank"},
    )
    if observed_universe_content_fingerprint != expected_universe_content_fingerprint:
        raise ValueError(f"{operation} universe content does not match its pinned SHA-256 fingerprint.")
    if marker_parameters is not None:
        if marker_source.get("operation") != "filter_marker_genes":
            raise ValueError(f"{operation} requires a filtered marker-evidence artifact.")
        if marker_parameters.get("marker_evidence_schema_version") != 2:
            raise ValueError(f"{operation} marker table uses an unsupported marker-evidence schema.")
        if marker_parameters.get("producer_node_id") != "OpenBioSingleCellFilterMarkerGenes":
            raise ValueError(f"{operation} marker table has the wrong producer node.")
        if marker_parameters.get("artifact_role") != "marker_table":
            raise ValueError(f"{operation} marker table has the wrong artifact role.")
        if marker_parameters.get("analysis_fingerprint") != expected_analysis_fingerprint:
            raise ValueError(f"{operation} marker analysis fingerprint changed.")
        if marker_parameters.get("ranking_fingerprint") != expected_ranking_fingerprint:
            raise ValueError(f"{operation} marker ranking fingerprint changed.")
        if marker_parameters.get("universe_fingerprint") != expected_universe_fingerprint:
            raise ValueError(f"{operation} marker universe fingerprint changed.")
        if marker_parameters.get("content_fingerprint") != expected_marker_content_fingerprint:
            raise ValueError(f"{operation} marker current-content fingerprint changed.")
        upstream_content_fingerprint = marker_parameters.get("upstream_content_fingerprint")
        if not valid_sha256(upstream_content_fingerprint):
            raise ValueError(f"{operation} marker upstream-content fingerprint is invalid.")
    if universe_parameters is not None:
        if universe_source.get("operation") != "marker_genes":
            raise ValueError(f"{operation} universe must come directly from Marker Genes.")
        if universe_parameters.get("marker_evidence_schema_version") != 2:
            raise ValueError(f"{operation} universe uses an unsupported marker-evidence schema.")
        if universe_parameters.get("producer_node_id") != "OpenBioSingleCellMarkerGenes":
            raise ValueError(f"{operation} universe has the wrong producer node.")
        if universe_parameters.get("artifact_role") != "tested_gene_universe":
            raise ValueError(f"{operation} universe has the wrong artifact role.")
        if universe_parameters.get("analysis_fingerprint") != expected_analysis_fingerprint:
            raise ValueError(f"{operation} universe analysis fingerprint changed.")
        if universe_parameters.get("ranking_fingerprint") != expected_ranking_fingerprint:
            raise ValueError(f"{operation} universe ranking fingerprint changed.")
        if universe_parameters.get("universe_fingerprint") != expected_universe_fingerprint:
            raise ValueError(f"{operation} universe fingerprint changed.")
        if universe_parameters.get("content_fingerprint") != expected_universe_content_fingerprint:
            raise ValueError(f"{operation} universe current-content fingerprint changed.")

    universe_genes = []
    for value in universe_frame["gene"].tolist():
        if not isinstance(value, str) or not value.strip() or value != value.strip():
            raise ValueError(f"{operation} universe genes must be nonblank strings without surrounding whitespace.")
        universe_genes.append(value)
    if not universe_genes or len(universe_genes) != len(set(universe_genes)):
        raise ValueError(f"{operation} universe genes must be nonempty and unique.")
    universe_ranks = universe_frame["universe_rank"]
    if pd.api.types.is_bool_dtype(universe_ranks.dtype) or not pd.api.types.is_numeric_dtype(universe_ranks.dtype):
        raise TypeError(f"{operation} universe ranks must be numeric.")
    rank_values = universe_ranks.to_numpy(dtype=float)
    if not bool(np.array_equal(rank_values, np.arange(1, len(universe_genes) + 1, dtype=float))):
        raise ValueError(f"{operation} universe ranks must be consecutive and one-based.")
    observed_universe_fingerprint = sha256_json(
        {
            "schema": "openbio-singlecell/tested-gene-universe-identity/v2",
            "ordered_genes": universe_genes,
        }
    )
    if observed_universe_fingerprint != expected_universe_fingerprint:
        raise ValueError(f"{operation} universe content does not match its SHA-256 identity.")
    if marker_input_genes is not None and int(marker_input_genes) != len(universe_genes):
        raise ValueError(f"{operation} marker artifact reports the wrong universe size.")
    if universe_input_genes is not None and int(universe_input_genes) != len(universe_genes):
        raise ValueError(f"{operation} universe artifact reports the wrong universe size.")

    if not isinstance(group_labels, list) or any(not isinstance(group, str) or not group for group in group_labels):
        raise ValueError(f"{operation} requires ordered nonblank group labels.")
    if len(group_labels) < 2 or len(group_labels) != len(set(group_labels)):
        raise ValueError(f"{operation} requires at least two unique group labels.")
    if not isinstance(upstream_parameters, Mapping):
        raise TypeError(f"{operation} upstream_parameters must be a mapping.")
    if marker_parameters is not None and marker_parameters != dict(upstream_parameters):
        raise ValueError(f"{operation} marker provenance differs from the pinned upstream parameters.")
    if upstream_parameters.get("group_labels") != group_labels:
        raise ValueError(f"{operation} ordered group labels differ from the pinned upstream provenance.")
    actual_n_genes = upstream_parameters.get("actual_n_genes_per_group")
    if isinstance(actual_n_genes, bool) or not isinstance(actual_n_genes, int) or actual_n_genes < 1:
        raise ValueError(f"{operation} upstream ranking depth is invalid.")
    if marker_ranks.size and int(marker_ranks.max()) > actual_n_genes:
        raise ValueError(f"{operation} marker table contains a rank outside the upstream ranking depth.")
    observed_group_order = list(dict.fromkeys(marker_frame["group"].tolist()))
    expected_observed_order = [group for group in group_labels if group in set(observed_group_order)]
    if observed_group_order != expected_observed_order:
        raise ValueError(f"{operation} marker group blocks do not follow the pinned category order.")
    if bool(upstream_parameters.get("ranking_truncated")):
        raise ValueError(f"{operation} requires a non-truncated upstream ranking before marker filtering.")
    selected_by_group = {group: [] for group in group_labels}
    seen_pairs = set()
    universe_set = set(universe_genes)
    for row in marker_frame.itertuples(index=False):
        group = row.group
        gene = row.gene
        if not isinstance(group, str) or not group.strip() or group != group.strip() or group not in selected_by_group:
            raise ValueError(f"{operation} marker table contains an invalid group identifier {group!r}.")
        if not isinstance(gene, str) or not gene.strip() or gene != gene.strip():
            raise ValueError(f"{operation} marker genes must be nonblank strings without surrounding whitespace.")
        pair = (group, gene)
        if pair in seen_pairs:
            raise ValueError(f"{operation} marker table contains duplicate (group, gene) rows.")
        seen_pairs.add(pair)
        if gene not in universe_set:
            raise ValueError(f"{operation} selected marker {gene!r} is outside the tested-gene universe.")
        try:
            direction = float(row.log2_fold_change_approx)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{operation} marker log2 fold changes must be numeric.") from exc
        if not math.isfinite(direction) or direction <= 0.0:
            raise ValueError(f"{operation} requires strictly positive selected marker direction.")
        selected_by_group[group].append(gene)
    groups_without_selected_markers = [group for group, genes in selected_by_group.items() if not genes]
    for group, genes in selected_by_group.items():
        if len(genes) != len(set(genes)):
            raise ValueError(f"{operation} selected marker set is not unique for group {group!r}.")
    filter_thresholds = {
        "min_log2_fold_change": upstream_parameters.get("min_log2_fold_change"),
        "min_fraction_in_group": upstream_parameters.get("min_fraction_in_group"),
        "max_fraction_reference": upstream_parameters.get("max_fraction_reference"),
        "max_p_adjusted": upstream_parameters.get("max_p_adjusted"),
    }
    for field in (
        "min_log2_fold_change",
        "min_fraction_in_group",
        "max_fraction_reference",
        "max_p_adjusted",
    ):
        value = filter_thresholds[field]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError(f"{operation} filtered marker provenance has invalid {field!r}.")
        filter_thresholds[field] = float(value)
    selected_numeric = marker_frame[
        ["log2_fold_change_approx", "fraction_in_group", "fraction_reference", "p_adjusted"]
    ]
    if any(
        pd.api.types.is_bool_dtype(selected_numeric[column].dtype)
        or not pd.api.types.is_numeric_dtype(selected_numeric[column].dtype)
        for column in selected_numeric
    ):
        raise TypeError(f"{operation} filtered marker columns must remain numeric.")
    selected_values = selected_numeric.to_numpy(dtype=float)
    if not bool(np.isfinite(selected_values).all()):
        raise ValueError(f"{operation} filtered marker evidence contains non-finite values.")
    selected_mask = (
        (selected_numeric["log2_fold_change_approx"].to_numpy(dtype=float) >= filter_thresholds["min_log2_fold_change"])
        & (selected_numeric["fraction_in_group"].to_numpy(dtype=float) >= filter_thresholds["min_fraction_in_group"])
        & (selected_numeric["fraction_reference"].to_numpy(dtype=float) <= filter_thresholds["max_fraction_reference"])
        & (selected_numeric["p_adjusted"].to_numpy(dtype=float) <= filter_thresholds["max_p_adjusted"])
    )
    if not bool(selected_mask.all()):
        raise ValueError(f"{operation} marker rows do not satisfy their recorded filter thresholds.")

    for name, value in (("min_targets", min_targets), ("min_overlap", min_overlap)):
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
    if not isinstance(source_column, str) or not source_column.strip():
        raise ValueError(f"{operation} source_column cannot be empty.")
    if not isinstance(target_column, str) or not target_column.strip():
        raise ValueError(f"{operation} target_column cannot be empty.")
    source_column = source_column.strip()
    target_column = target_column.strip()
    if source_column == target_column:
        raise ValueError(f"{operation} source_column and target_column must be different.")

    if not isinstance(resource_metadata_json, str):
        raise TypeError(f"{operation} resource_metadata_json must be a string.")
    if len(resource_metadata_json.encode("utf-8")) > 65_536:
        raise ValueError(f"{operation} resource metadata exceeds 65,536 bytes.")

    def reject_constant(value):
        raise ValueError(f"{operation} resource metadata contains non-standard JSON constant {value!r}.")

    def metadata_from_pairs(pairs):
        result = {}
        for raw_key, value in pairs:
            if not isinstance(raw_key, str) or not raw_key.strip():
                raise ValueError(f"{operation} resource metadata keys must be nonblank strings.")
            key = raw_key.strip()
            if key in result:
                raise ValueError(f"{operation} resource metadata contains duplicate key {key!r}.")
            result[key] = value
        return result

    try:
        resource_metadata = json.loads(
            resource_metadata_json,
            object_pairs_hook=metadata_from_pairs,
            parse_constant=reject_constant,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"{operation} resource metadata is invalid JSON ({exc.msg}).") from exc
    required_metadata = {
        "name",
        "version",
        "date",
        "organism",
        "identifier_namespace",
        "scope",
        "license",
        "citation",
    }
    if not isinstance(resource_metadata, dict) or set(resource_metadata) != required_metadata:
        missing = sorted(required_metadata - set(resource_metadata)) if isinstance(resource_metadata, dict) else []
        unknown = sorted(set(resource_metadata) - required_metadata) if isinstance(resource_metadata, dict) else []
        raise ValueError(
            f"{operation} resource metadata must contain exactly the required fields; "
            f"missing={missing}, unknown={unknown}."
        )
    for key, value in list(resource_metadata.items()):
        if not isinstance(value, str) or not value.strip():
            raise TypeError(f"{operation} resource metadata field {key!r} must be a nonblank string.")
        resource_metadata[key] = value.strip()
    for optional_identity in ("organism", "identifier_namespace"):
        upstream_identity = upstream_parameters.get(optional_identity)
        if upstream_identity is not None and upstream_identity != resource_metadata[optional_identity]:
            raise ValueError(
                f"{operation} resource {optional_identity} conflicts with marker provenance: "
                f"{resource_metadata[optional_identity]!r} != {upstream_identity!r}."
            )

    if not isinstance(resource_path, str) or not resource_path.strip():
        raise ValueError(f"{operation} resource path cannot be empty.")
    if requested_resource_path is None:
        requested_resource_path = resource_path
    if not isinstance(requested_resource_path, str) or not requested_resource_path.strip():
        raise ValueError(f"{operation} requested resource path cannot be empty.")
    requested_resource_path = requested_resource_path.strip()
    if not isinstance(openbio_version, str) or not openbio_version.strip():
        raise ValueError(f"{operation} openbio_version cannot be empty.")
    openbio_version = openbio_version.strip()
    resource_path = os.path.realpath(resource_path)
    if not os.path.isfile(resource_path):
        raise FileNotFoundError(f"{operation} resource file not found: {resource_path}")
    if not resource_path.lower().endswith(".csv"):
        raise ValueError(f"{operation} resource must be a CSV file.")
    resource_size = int(os.path.getsize(resource_path))
    if resource_size > 512 * 1024 * 1024:
        raise ValueError(f"{operation} resource exceeds the 512 MiB safety limit.")
    digest = hashlib.sha256()
    with open(resource_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    resource_sha256 = digest.hexdigest()
    if expected_resource_sha256 is not None and resource_sha256 != expected_resource_sha256:
        raise ValueError(
            f"{operation} resource fingerprint changed: expected {expected_resource_sha256}, "
            f"observed {resource_sha256}."
        )
    resource_row_count = 0
    with open(resource_path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{operation} resource CSV is empty.") from exc
        for resource_row_count, row in enumerate(reader, start=1):
            if len(row) != len(header):
                raise ValueError(
                    f"{operation} resource CSV data row {resource_row_count} has {len(row)} fields; "
                    f"expected exactly {len(header)} from the header."
                )
            if resource_row_count > 2_000_000:
                raise ValueError(f"{operation} resource exceeds the 2,000,000-row safety limit.")
    if any(not value.strip() for value in header):
        raise ValueError(f"{operation} resource CSV contains a blank column name.")
    trimmed_header = [value.strip() for value in header]
    if len(trimmed_header) != len(set(trimmed_header)):
        raise ValueError(f"{operation} resource CSV contains duplicate column names after trimming.")
    if source_column not in header or target_column not in header:
        raise ValueError(
            f"{operation} resource CSV is missing source/target columns: "
            f"{[name for name in (source_column, target_column) if name not in header]}."
        )
    try:
        resource = pd.read_csv(
            resource_path,
            dtype=str,
            encoding="utf-8-sig",
            keep_default_na=False,
            na_filter=False,
        )
    except Exception as exc:
        raise ValueError(f"{operation} could not parse the resource CSV.") from exc
    if len(resource) != resource_row_count:
        raise ValueError(
            f"{operation} resource CSV parser row-count mismatch: validated {resource_row_count}, "
            f"parsed {len(resource)}."
        )
    if resource.empty:
        raise ValueError(f"{operation} resource CSV contains no data rows.")
    network_rows = []
    for source, target in zip(resource[source_column].tolist(), resource[target_column].tolist(), strict=True):
        if not isinstance(source, str) or not source.strip():
            raise ValueError(f"{operation} resource contains a blank source identifier.")
        if not isinstance(target, str) or not target.strip():
            raise ValueError(f"{operation} resource contains a blank target identifier.")
        network_rows.append((source.strip(), target.strip()))
    unique_pairs = list(dict.fromkeys(network_rows))
    duplicate_rows_removed = len(network_rows) - len(unique_pairs)
    source_sets_before = {}
    for source, target in unique_pairs:
        source_sets_before.setdefault(source, set()).add(target)
    source_sets_in_universe = {source: targets & universe_set for source, targets in source_sets_before.items()}
    all_resource_targets_in_universe = set().union(*source_sets_in_universe.values())
    retained_target_rows = int(sum(len(targets) for targets in source_sets_in_universe.values()))
    if retained_target_rows == 0:
        raise ValueError(f"{operation} resource has zero target overlap with the tested-gene universe.")
    surviving_sources = sorted(
        source for source, targets in source_sets_in_universe.items() if len(targets) >= min_targets
    )
    removed_sources = sorted(set(source_sets_before) - set(surviving_sources))
    if not surviving_sources:
        raise ValueError(
            f"{operation} has no resource source with at least {min_targets} targets after universe intersection."
        )
    pruned_pairs = [
        (source, target) for source in surviving_sources for target in sorted(source_sets_in_universe[source])
    ]
    network = pd.DataFrame(pruned_pairs, columns=["source", "target"])
    tested_resource_target_union = set(network["target"].tolist())

    if decoupler_module is None:
        try:
            decoupler_module = importlib.import_module("decoupler")
        except (ImportError, OSError) as exc:
            raise RuntimeError(
                f"{operation} requires the decoupler 2.x package, but it is unavailable ({exc})."
            ) from exc
    version = getattr(decoupler_module, "__version__", None)
    if not isinstance(version, str) or not version.strip():
        try:
            version = importlib_metadata.version("decoupler")
        except importlib_metadata.PackageNotFoundError as exc:
            raise RuntimeError(f"{operation} cannot determine the installed decoupler version.") from exc
    try:
        decoupler_major = int(version.split(".", 1)[0])
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{operation} received an invalid decoupler version {version!r}.") from exc
    if decoupler_major != 2:
        raise RuntimeError(f"{operation} requires decoupler 2.x mt.query_set; installed version is {version!r}.")
    query_set = getattr(getattr(decoupler_module, "mt", None), "query_set", None)
    if not callable(query_set):
        raise RuntimeError(f"{operation} requires the decoupler 2.x mt.query_set API.")
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

    raw_rows = []
    per_group_backend = {}
    universe_count = len(universe_genes)
    for group in group_labels:
        selected = set(selected_by_group[group])
        backend = query_set(
            features=sorted(selected),
            net=network.copy(),
            alternative="greater",
            n_bg=universe_count,
            ha_corr=0.5,
            tmin=min_targets,
            verbose=False,
        )
        if not isinstance(backend, pd.DataFrame):
            raise RuntimeError(f"{operation} backend returned a non-DataFrame result for group {group!r}.")
        required_columns = ["source", "stat", "pval", "padj"]
        missing_columns = [column for column in required_columns if column not in backend]
        if missing_columns:
            raise RuntimeError(f"{operation} backend result for group {group!r} is missing columns: {missing_columns}.")
        if bool(backend["source"].duplicated().any()):
            raise RuntimeError(f"{operation} backend returned duplicate sources for group {group!r}.")
        returned_sources = backend["source"].tolist()
        if any(not isinstance(source, str) for source in returned_sources) or set(returned_sources) != set(
            surviving_sources
        ):
            raise RuntimeError(
                f"{operation} backend source set for group {group!r} differs from the retained resource."
            )
        indexed = backend.set_index("source")
        for column in ("stat", "pval", "padj"):
            series = indexed[column]
            if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
                raise RuntimeError(
                    f"{operation} backend column {column!r} for group {group!r} must be numeric and non-boolean."
                )
        try:
            backend_values = indexed.loc[surviving_sources, ["stat", "pval", "padj"]].to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{operation} backend returned nonnumeric values for group {group!r}.") from exc
        if not bool(np.isfinite(backend_values).all()):
            raise RuntimeError(f"{operation} backend returned non-finite values for group {group!r}.")
        if bool(((backend_values[:, 1:] < 0.0) | (backend_values[:, 1:] > 1.0)).any()):
            raise RuntimeError(f"{operation} backend returned p-values outside [0, 1] for group {group!r}.")
        independent_p_values = []
        independent_statistics = []
        group_rows = []
        for source in surviving_sources:
            resource_set = source_sets_in_universe[source]
            overlap = sorted(selected & resource_set)
            a = len(overlap)
            b = len(resource_set - selected)
            c = len(selected - resource_set)
            d = universe_count - a - b - c
            if min(a, b, c, d) < 0 or a + b + c + d != universe_count:
                raise RuntimeError(
                    f"{operation} constructed an invalid contingency for group {group!r}, source {source!r}."
                )
            log_odds = float(np.log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))))
            p_value = float(stats.fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
            independent_statistics.append(log_odds)
            independent_p_values.append(p_value)
            group_rows.append(
                {
                    "group": group,
                    "source": source,
                    "selected_count": len(selected),
                    "universe_count": universe_count,
                    "set_size_before_universe": len(source_sets_before[source]),
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
        independent_p_values_array = np.asarray(independent_p_values, dtype=float)
        independent_statistics_array = np.asarray(independent_statistics, dtype=float)
        independent_adjusted = np.asarray(
            stats.false_discovery_control(independent_p_values_array, method="bh"), dtype=float
        )
        if not bool(
            np.allclose(
                backend_values[:, 0],
                independent_statistics_array,
                rtol=1e-10,
                atol=1e-12,
            )
        ):
            raise RuntimeError(f"{operation} backend log odds disagree with independent contingency verification.")
        if not bool(
            np.allclose(
                backend_values[:, 1],
                independent_p_values_array,
                rtol=1e-10,
                atol=1e-12,
            )
        ):
            raise RuntimeError(f"{operation} backend p-values disagree with independent Fisher exact tests.")
        if not bool(
            np.allclose(
                backend_values[:, 2],
                independent_adjusted,
                rtol=1e-10,
                atol=1e-12,
            )
        ):
            raise RuntimeError(f"{operation} backend adjusted p-values disagree with within-group BH correction.")
        for row, adjusted in zip(group_rows, independent_adjusted.tolist(), strict=True):
            row["p_adj_within_group"] = float(adjusted)
            raw_rows.append(row)
        per_group_backend[group] = {
            "sources": len(surviving_sources),
            "alternative": "greater",
            "n_bg": universe_count,
            "ha_corr": 0.5,
            "tmin": min_targets,
        }

    evidence = pd.DataFrame(raw_rows)
    evidence["p_adj_global"] = stats.false_discovery_control(evidence["p_value"].to_numpy(dtype=float), method="bh")
    group_position = {group: index for index, group in enumerate(group_labels)}
    evidence["_group_position"] = evidence["group"].map(group_position)
    evidence = evidence.sort_values(
        [
            "_group_position",
            "p_adj_within_group",
            "p_value",
            "log_odds_ratio",
            "overlap_count",
            "source",
        ],
        ascending=[True, True, True, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    ranks = []
    for group in group_labels:
        group_frame = evidence.loc[evidence["group"] == group]
        previous_key = None
        current_rank = 0
        for row in group_frame.itertuples(index=False):
            scientific_key = (
                float(row.p_adj_within_group),
                float(row.p_value),
                float(row.log_odds_ratio),
                int(row.overlap_count),
            )
            if scientific_key != previous_key:
                current_rank += 1
                previous_key = scientific_key
            ranks.append(current_rank)
    evidence["rank_within_group"] = ranks
    evidence["positive_enrichment"] = (evidence["selected_count"] > 0) & (evidence["log_odds_ratio"] > 0.0)
    evidence["passes_min_overlap"] = evidence["overlap_count"] >= min_overlap
    evidence["significant_within_group"] = evidence["p_adj_within_group"] <= max_p_adjusted
    statuses = []
    for row in evidence.itertuples(index=False):
        if int(row.selected_count) == 0:
            statuses.append("unresolved_no_selected_markers")
        elif not bool(row.positive_enrichment):
            statuses.append("nonpositive")
        elif not bool(row.passes_min_overlap):
            statuses.append("insufficient_overlap")
        elif not bool(row.significant_within_group):
            statuses.append("not_significant")
        else:
            statuses.append("eligible")
    evidence["candidate_status"] = statuses
    for group in group_labels:
        group_indices = evidence.index[evidence["group"] == group]
        eligible_indices = [index for index in group_indices if evidence.at[index, "candidate_status"] == "eligible"]
        if eligible_indices:
            best_rank = min(int(evidence.at[index, "rank_within_group"]) for index in eligible_indices)
            tied = [index for index in eligible_indices if int(evidence.at[index, "rank_within_group"]) == best_rank]
            if len(tied) > 1:
                evidence.loc[tied, "candidate_status"] = "tied_best"
    evidence = evidence.drop(columns=["_group_position"])
    evidence = evidence[evidence_columns]

    count_columns = [
        "selected_count",
        "universe_count",
        "set_size_before_universe",
        "set_size_in_universe",
        "overlap_count",
        "a",
        "b",
        "c",
        "d",
        "rank_within_group",
    ]
    for column in count_columns:
        evidence[column] = evidence[column].astype(int)
        if bool((evidence[column] < 0).any()):
            raise RuntimeError(f"{operation} produced a negative count in {column!r}.")
    float_columns = ["log_odds_ratio", "p_value", "p_adj_within_group", "p_adj_global"]
    numeric = evidence[float_columns].to_numpy(dtype=float)
    if not bool(np.isfinite(numeric).all()):
        raise RuntimeError(f"{operation} produced non-finite evidence values.")
    if bool(
        (
            (evidence[["p_value", "p_adj_within_group", "p_adj_global"]] < 0.0)
            | (evidence[["p_value", "p_adj_within_group", "p_adj_global"]] > 1.0)
        )
        .any()
        .any()
    ):
        raise RuntimeError(f"{operation} produced evidence p-values outside [0, 1].")
    if len(evidence) != len(group_labels) * len(surviving_sources):
        raise RuntimeError(f"{operation} produced an incomplete group-by-source evidence grid.")
    if bool(evidence.duplicated(subset=["group", "source"]).any()):
        raise RuntimeError(f"{operation} produced duplicate (group, source) rows.")
    for row in evidence.itertuples(index=False):
        if row.a != row.overlap_count or row.a + row.b != row.set_size_in_universe:
            raise RuntimeError(f"{operation} result contingency identities failed.")
        if row.a + row.c != row.selected_count or row.a + row.b + row.c + row.d != row.universe_count:
            raise RuntimeError(f"{operation} result contingency identities failed.")

    def numeric_summary(values):
        array = np.asarray(list(values), dtype=float)
        quantiles = np.quantile(array, [0.0, 0.25, 0.5, 0.75, 1.0])
        return {
            "n": int(array.size),
            "min": float(quantiles[0]),
            "q1": float(quantiles[1]),
            "median": float(quantiles[2]),
            "mean": float(array.mean()),
            "q3": float(quantiles[3]),
            "max": float(quantiles[4]),
        }

    selected_count_pairs = [(group, len(selected_by_group[group])) for group in group_labels]
    unmatched_pairs = [
        (group, len(set(selected_by_group[group]) - all_resource_targets_in_universe)) for group in group_labels
    ]
    unresolved = []
    ambiguous = []
    top_candidates = {}
    for group in group_labels:
        group_frame = evidence.loc[evidence["group"] == group]
        eligible = group_frame[group_frame["candidate_status"].isin(["eligible", "tied_best"])]
        if eligible.empty:
            unresolved.append(group)
        if int((group_frame["candidate_status"] == "tied_best").sum()) > 1:
            ambiguous.append(group)
        top_candidates[group] = [
            {
                "source": row.source,
                "rank": int(row.rank_within_group),
                "status": row.candidate_status,
                "log_odds_ratio": float(row.log_odds_ratio),
                "p_value": float(row.p_value),
                "p_adj_within_group": float(row.p_adj_within_group),
                "p_adj_global": float(row.p_adj_global),
                "overlap_count": int(row.overlap_count),
                "overlap_genes": json.loads(row.overlap_genes),
            }
            for row in group_frame.head(5).itertuples(index=False)
        ]
    group_preview_limit = 100
    preview_groups = group_labels[:group_preview_limit]
    source_before_sizes = [len(source_sets_before[source]) for source in sorted(source_sets_before)]
    source_after_sizes = [len(source_sets_in_universe[source]) for source in sorted(source_sets_before)]
    diagnostics = {
        "analysis_fingerprint_sha256": expected_analysis_fingerprint,
        "ranking_fingerprint_sha256": expected_ranking_fingerprint,
        "universe_fingerprint_sha256": expected_universe_fingerprint,
        "marker_content_fingerprint_sha256": expected_marker_content_fingerprint,
        "universe_content_fingerprint_sha256": expected_universe_content_fingerprint,
        "marker_method": upstream_parameters.get("marker_method"),
        "marker_groupby": upstream_parameters.get("marker_groupby"),
        "marker_filter_thresholds": filter_thresholds,
        "ranking_truncated": False,
        "group_count": len(group_labels),
        "groups_without_selected_markers_count": len(groups_without_selected_markers),
        "groups_without_selected_markers_preview": groups_without_selected_markers[:100],
        "groups_without_selected_markers_preview_truncated": len(groups_without_selected_markers) > 100,
        "group_preview_truncated": len(group_labels) > group_preview_limit,
        "selected_counts_by_group_preview": {group: dict(selected_count_pairs)[group] for group in preview_groups},
        "unmatched_selected_counts_by_group_preview": {group: dict(unmatched_pairs)[group] for group in preview_groups},
        "universe_count": universe_count,
        "resource": {
            "resolved_path": resource_path,
            "size_bytes": resource_size,
            "sha256": resource_sha256,
            "metadata": resource_metadata,
            "input_rows": len(network_rows),
            "unique_pairs": len(unique_pairs),
            "duplicate_rows_removed": duplicate_rows_removed,
            "source_count_before": len(source_sets_before),
            "source_count_after": len(surviving_sources),
            "sources_removed_by_min_targets_count": len(removed_sources),
            "sources_removed_by_min_targets_preview": removed_sources[:20],
            "sources_removed_preview_truncated": len(removed_sources) > 20,
            "targets_before_universe_count": len({target for _, target in unique_pairs}),
            "targets_in_universe_count": len(all_resource_targets_in_universe),
            "tested_targets_in_universe_count": len(tested_resource_target_union),
            "target_rows_in_universe": retained_target_rows,
            "set_size_before_universe": numeric_summary(source_before_sizes),
            "set_size_in_universe": numeric_summary(source_after_sizes),
        },
        "fixed_policy": {
            "alternative": "greater",
            "n_bg": universe_count,
            "ha_corr": 0.5,
            "min_targets_after_universe_intersection": min_targets,
            "within_group_correction": "Benjamini-Hochberg across retained sources",
            "global_correction": "Benjamini-Hochberg across all group-source pairs",
            "full_evidence_rows_retained": True,
        },
        "thresholds": {
            "min_overlap": min_overlap,
            "max_p_adjusted": max_p_adjusted,
        },
        "tested_pairs": len(evidence),
        "positive_rows": int(evidence["positive_enrichment"].sum()),
        "significant_within_group_rows": int(evidence["significant_within_group"].sum()),
        "eligible_rows": int(evidence["candidate_status"].isin(["eligible", "tied_best"]).sum()),
        "unresolved_group_count": len(unresolved),
        "unresolved_groups_preview": unresolved[:100],
        "unresolved_preview_truncated": len(unresolved) > 100,
        "ambiguous_group_count": len(ambiguous),
        "ambiguous_groups_preview": ambiguous[:100],
        "ambiguous_preview_truncated": len(ambiguous) > 100,
        "top_candidates_by_group_preview": {group: top_candidates[group] for group in preview_groups},
        "top_candidates_group_preview_truncated": len(group_labels) > group_preview_limit,
        "backend": {
            "decoupler_version": str(version),
            "api": "decoupler.mt.query_set",
            "calls": len(group_labels),
            "per_group_contract_preview": {group: per_group_backend[group] for group in preview_groups},
            "per_group_preview_truncated": len(group_labels) > group_preview_limit,
            "scipy_version": scipy.__version__,
        },
    }
    parameters = {
        "producer_node_id": "OpenBioSingleCellMarkerORAEvidence",
        "artifact_role": "marker_ora_evidence_table",
        "analysis_fingerprint": expected_analysis_fingerprint,
        "ranking_fingerprint": expected_ranking_fingerprint,
        "universe_fingerprint": expected_universe_fingerprint,
        "marker_content_fingerprint": expected_marker_content_fingerprint,
        "universe_content_fingerprint": expected_universe_content_fingerprint,
        "upstream_marker_content_fingerprint": upstream_parameters["upstream_content_fingerprint"],
        "resource_csv": requested_resource_path,
        "resource_sha256": resource_sha256,
        "resource_metadata": resource_metadata,
        "source_column": source_column,
        "target_column": target_column,
        "min_targets": min_targets,
        "min_overlap": min_overlap,
        "max_p_adjusted": max_p_adjusted,
        "alternative": "greater",
        "n_bg": universe_count,
        "ha_corr": 0.5,
        "within_group_correction": "benjamini-hochberg",
        "global_correction": "benjamini-hochberg",
        "table_content_fingerprint_sha256": marker_ora_table_fingerprint(evidence),
    }
    warnings = [
        "Marker ORA results are exploratory Cluster marker evidence for annotation review, not cell-type probabilities or Curated annotation.",
        "Resource organism, identifier namespace, scope, license, and citation are caller declarations; the file hash identifies bytes but not biological suitability.",
        "Marker selection and ORA reuse the same cells and are not independent confirmatory tests.",
    ]
    if "organism" not in upstream_parameters or "identifier_namespace" not in upstream_parameters:
        warnings.append(
            "Upstream marker artifacts do not declare organism and identifier namespace, so compatibility with the resource could not be independently verified."
        )
    if duplicate_rows_removed:
        warnings.append(
            f"Collapsed {duplicate_rows_removed:,} duplicate resource source-target rows before testing."
        )
    if removed_sources:
        warnings.append(
            f"Removed {len(removed_sources):,} resource sources below min_targets after universe intersection."
        )
    if unresolved:
        warnings.append(
            f"{len(unresolved):,} groups had no eligible ORA candidate under the descriptive gates."
        )
    if groups_without_selected_markers:
        warnings.append(
            f"{len(groups_without_selected_markers):,} groups had no marker surviving the upstream filter; "
            "their complete ORA rows are retained as unresolved evidence."
        )
    diagnostics["parameters"] = parameters
    diagnostics["warnings"] = warnings
    diagnostics["software_versions"] = _standalone_marker_ora_software_versions(
        openbio_version=openbio_version,
        decoupler_version=str(version),
        scipy_version=scipy.__version__,
        pandas_version=pd.__version__,
        numpy_version=np.__version__,
    )
    json.dumps(diagnostics, ensure_ascii=False, allow_nan=False)
    return evidence, diagnostics


def run_marker_ora_evidence(
    marker_input: Any,
    universe_input: Any,
    *,
    expected_analysis_fingerprint: str,
    expected_ranking_fingerprint: str,
    expected_universe_fingerprint: str,
    expected_marker_content_fingerprint: str,
    expected_universe_content_fingerprint: str,
    group_labels: list[str],
    upstream_parameters: dict[str, Any],
    resource_path: str,
    requested_resource_path: str,
    resource_metadata_json: str,
    source_column: str,
    target_column: str,
    min_targets: int,
    min_overlap: int,
    max_p_adjusted: float,
    openbio_version: str,
    decoupler_module: Any,
) -> tuple[DataFrame, dict[str, Any]]:
    return _standalone_marker_ora_evidence(
        marker_input,
        universe_input,
        expected_analysis_fingerprint=expected_analysis_fingerprint,
        expected_ranking_fingerprint=expected_ranking_fingerprint,
        expected_universe_fingerprint=expected_universe_fingerprint,
        expected_marker_content_fingerprint=expected_marker_content_fingerprint,
        expected_universe_content_fingerprint=expected_universe_content_fingerprint,
        group_labels=group_labels,
        upstream_parameters=upstream_parameters,
        resource_path=resource_path,
        requested_resource_path=requested_resource_path,
        resource_metadata_json=resource_metadata_json,
        source_column=source_column,
        target_column=target_column,
        min_targets=min_targets,
        min_overlap=min_overlap,
        max_p_adjusted=max_p_adjusted,
        openbio_version=openbio_version,
        decoupler_module=decoupler_module,
    )


def build_marker_ora_summary(diagnostics: dict[str, Any]) -> dict[str, Any]:
    """Return the strict, standalone Marker ORA scientific summary."""
    return _standalone_marker_ora_summary(diagnostics)


def marker_ora_evidence_code(
    *,
    expected_analysis_fingerprint: str,
    expected_ranking_fingerprint: str,
    expected_universe_fingerprint: str,
    expected_marker_content_fingerprint: str,
    expected_universe_content_fingerprint: str,
    group_labels: list[str],
    upstream_parameters: dict[str, Any],
    resource_path: str,
    requested_resource_path: str,
    resource_metadata_json: str,
    expected_resource_sha256: str,
    source_column: str,
    target_column: str,
    min_targets: int,
    min_overlap: int,
    max_p_adjusted: float,
    openbio_version: str,
) -> str:
    implementations = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (
            marker_ora_table_fingerprint,
            _standalone_marker_ora_software_versions,
            _standalone_marker_ora_summary,
            _standalone_marker_ora_evidence,
        )
    )
    return f"""from __future__ import annotations

ORA_EVIDENCE_COLUMNS = {ORA_EVIDENCE_COLUMNS!r}

{implementations}


def run_marker_ora_evidence(marker_table, universe):
    evidence, diagnostics = _standalone_marker_ora_evidence(
        marker_table,
        universe,
        expected_analysis_fingerprint={expected_analysis_fingerprint!r},
        expected_ranking_fingerprint={expected_ranking_fingerprint!r},
        expected_universe_fingerprint={expected_universe_fingerprint!r},
        expected_marker_content_fingerprint={expected_marker_content_fingerprint!r},
        expected_universe_content_fingerprint={expected_universe_content_fingerprint!r},
        group_labels={group_labels!r},
        upstream_parameters={upstream_parameters!r},
        resource_path={resource_path!r},
        requested_resource_path={requested_resource_path!r},
        resource_metadata_json={resource_metadata_json!r},
        expected_resource_sha256={expected_resource_sha256!r},
        source_column={source_column!r},
        target_column={target_column!r},
        min_targets={min_targets!r},
        min_overlap={min_overlap!r},
        max_p_adjusted={max_p_adjusted!r},
        openbio_version={openbio_version!r},
    )
    return evidence, _standalone_marker_ora_summary(diagnostics)
"""


__all__ = [
    "ORA_EVIDENCE_COLUMNS",
    "build_marker_ora_summary",
    "marker_ora_table_fingerprint",
    "marker_ora_evidence_code",
    "run_marker_ora_evidence",
]
