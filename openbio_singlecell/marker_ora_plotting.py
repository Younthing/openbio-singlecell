from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .ora_evidence import ORA_EVIDENCE_COLUMNS, marker_ora_table_fingerprint


def _standalone_marker_ora_evidence_plot(
    table,
    *,
    producer,
    view=None,
    _return_details=False,
):
    import io
    import json
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Marker ORA Evidence Plot"
    if view is None:
        mode = {"view": "enrichment_dot", "max_terms_per_group": 10}
    elif not isinstance(view, Mapping):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant not in {"enrichment_dot", "overlap_bar"}:
        raise ValueError(f"{operation} view must be enrichment_dot or overlap_bar.")
    if set(mode) != {"view", "max_terms_per_group"}:
        raise ValueError(f"{operation} {variant} received inactive or unknown parameters.")
    maximum = mode["max_terms_per_group"]
    if isinstance(maximum, (bool, np.bool_)) or not isinstance(maximum, (int, np.integer)):
        raise TypeError(f"{operation} max_terms_per_group must be an integer.")
    maximum = int(maximum)
    if not 1 <= maximum <= 30:
        raise ValueError(f"{operation} max_terms_per_group must be between 1 and 30.")
    mode["max_terms_per_group"] = maximum

    if not isinstance(table, pd.DataFrame) or list(table.columns) != ORA_EVIDENCE_COLUMNS:
        raise ValueError(f"{operation} requires the exact canonical Marker ORA evidence columns.")
    if table.empty:
        raise ValueError(f"{operation} requires at least one stored group-resource row.")
    if len(table) > 10_000:
        raise ValueError(f"{operation} supports at most 10,000 stored evidence rows.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError(f"{operation} requires the canonical zero-based row index.")
    if not isinstance(producer, Mapping) or set(producer) != {
        "operation",
        "parameters",
        "input_cells",
        "input_genes",
    }:
        raise ValueError(f"{operation} requires exact producer provenance.")
    if producer["operation"] != "marker_ora_evidence":
        raise ValueError(f"{operation} requires the Marker ORA Evidence producer.")
    parameters = producer["parameters"]
    required_parameters = {
        "producer_node_id",
        "artifact_role",
        "analysis_fingerprint",
        "ranking_fingerprint",
        "universe_fingerprint",
        "marker_content_fingerprint",
        "universe_content_fingerprint",
        "upstream_marker_content_fingerprint",
        "resource_csv",
        "resource_sha256",
        "resource_metadata",
        "source_column",
        "target_column",
        "min_targets",
        "min_overlap",
        "max_p_adjusted",
        "alternative",
        "n_bg",
        "ha_corr",
        "within_group_correction",
        "global_correction",
        "table_content_fingerprint_sha256",
    }
    if not isinstance(parameters, Mapping) or set(parameters) != required_parameters:
        raise ValueError(f"{operation} producer parameter schema is invalid.")
    if (
        parameters["producer_node_id"] != "OpenBioSingleCellMarkerORAEvidence"
        or parameters["artifact_role"] != "marker_ora_evidence_table"
        or parameters["alternative"] != "greater"
        or parameters["ha_corr"] != 0.5
        or parameters["within_group_correction"] != "benjamini-hochberg"
        or parameters["global_correction"] != "benjamini-hochberg"
    ):
        raise ValueError(f"{operation} producer scientific contract is invalid.")
    for field in (
        "analysis_fingerprint",
        "ranking_fingerprint",
        "universe_fingerprint",
        "marker_content_fingerprint",
        "universe_content_fingerprint",
        "upstream_marker_content_fingerprint",
        "resource_sha256",
        "table_content_fingerprint_sha256",
    ):
        value = parameters[field]
        if (
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
        ):
            raise ValueError(f"{operation} producer {field} must be a lowercase SHA-256 value.")
    for field in ("min_targets", "min_overlap", "n_bg"):
        value = parameters[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{operation} producer {field} must be a positive integer.")
    max_p_adjusted = parameters["max_p_adjusted"]
    if (
        isinstance(max_p_adjusted, bool)
        or not isinstance(max_p_adjusted, (int, float))
        or not np.isfinite(max_p_adjusted)
        or not 0 <= max_p_adjusted <= 1
    ):
        raise ValueError(f"{operation} producer max_p_adjusted must lie in [0, 1].")
    if producer["input_genes"] != parameters["n_bg"]:
        raise ValueError(f"{operation} tested-gene universe size disagrees with producer provenance.")
    if isinstance(producer["input_cells"], bool) or not isinstance(producer["input_cells"], int) or producer[
        "input_cells"
    ] < 1:
        raise ValueError(f"{operation} producer input_cells must be a positive integer.")

    for column in ("group", "source", "overlap_genes", "candidate_status"):
        if any(
            not isinstance(value, str) or (not value and column != "overlap_genes") or value != value.strip()
            for value in table[column].tolist()
        ):
            raise ValueError(f"{operation} {column!r} must contain canonical strings.")
    if bool(table.duplicated(["group", "source"]).any()):
        raise ValueError(f"{operation} contains duplicate group-resource rows.")
    groups = list(dict.fromkeys(table["group"].tolist()))
    active_group = None
    seen_groups = set()
    for group in table["group"]:
        if group != active_group:
            if group in seen_groups:
                raise ValueError(f"{operation} group blocks must be contiguous.")
            seen_groups.add(group)
            active_group = group
    sources_by_group = {
        group: set(table.loc[table["group"] == group, "source"].tolist()) for group in groups
    }
    first_sources = sources_by_group[groups[0]]
    if not first_sources or any(sources != first_sources for sources in sources_by_group.values()):
        raise ValueError(f"{operation} requires the complete group-by-resource evidence grid.")

    integer_columns = [
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
    for column in integer_columns:
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_integer_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} must contain canonical integers.")
        if bool((series < 0).any()):
            raise ValueError(f"{operation} {column!r} cannot contain negative values.")
    if bool((table["rank_within_group"] < 1).any()):
        raise ValueError(f"{operation} within-group ranks must be positive.")
    if not bool((table["universe_count"] == parameters["n_bg"]).all()):
        raise ValueError(f"{operation} universe counts disagree with producer provenance.")
    if not bool((table["set_size_before_universe"] >= table["set_size_in_universe"]).all()):
        raise ValueError(f"{operation} resource set sizes are inconsistent.")
    if not bool(
        (
            (table["a"] == table["overlap_count"])
            & (table["a"] + table["b"] == table["set_size_in_universe"])
            & (table["a"] + table["c"] == table["selected_count"])
            & (table["a"] + table["b"] + table["c"] + table["d"] == table["universe_count"])
        ).all()
    ):
        raise ValueError(f"{operation} contingency identities are invalid.")

    overlap_lists = []
    for raw, count in zip(table["overlap_genes"], table["overlap_count"], strict=True):
        try:
            genes = json.loads(raw)
        except (TypeError, json.JSONDecodeError) as error:
            raise ValueError(f"{operation} overlap_genes must be canonical JSON arrays.") from error
        if (
            not isinstance(genes, list)
            or len(genes) != int(count)
            or len(set(genes)) != len(genes)
            or any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in genes)
            or json.dumps(genes, ensure_ascii=False, separators=(",", ":")) != raw
        ):
            raise ValueError(f"{operation} overlap gene evidence disagrees with overlap_count.")
        overlap_lists.append(genes)

    for column in ("log_odds_ratio", "p_value", "p_adj_within_group", "p_adj_global"):
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} must contain real numeric values.")
        if not bool(np.isfinite(series.to_numpy(dtype=float)).all()):
            raise ValueError(f"{operation} {column!r} must contain finite values.")
    expected_log_odds = np.log(
        ((table["a"].to_numpy(dtype=float) + 0.5) * (table["d"].to_numpy(dtype=float) + 0.5))
        / ((table["b"].to_numpy(dtype=float) + 0.5) * (table["c"].to_numpy(dtype=float) + 0.5))
    )
    if not bool(
        np.allclose(table["log_odds_ratio"], expected_log_odds, rtol=1e-12, atol=1e-12)
    ):
        raise ValueError(f"{operation} log odds disagree with stored contingency evidence.")
    raw_p = table["p_value"].to_numpy(dtype=float)
    within_p = table["p_adj_within_group"].to_numpy(dtype=float)
    global_p = table["p_adj_global"].to_numpy(dtype=float)
    if bool(((raw_p < 0) | (raw_p > 1)).any()) or bool(((within_p < 0) | (within_p > 1)).any()) or bool(
        ((global_p < 0) | (global_p > 1)).any()
    ):
        raise ValueError(f"{operation} p-values must lie in [0, 1].")
    if bool((within_p + 1e-12 < raw_p).any()) or bool((global_p + 1e-12 < raw_p).any()):
        raise ValueError(f"{operation} adjusted p-values cannot be smaller than raw p-values.")
    for column in ("positive_enrichment", "passes_min_overlap", "significant_within_group"):
        if not pd.api.types.is_bool_dtype(table[column].dtype):
            raise TypeError(f"{operation} {column!r} must contain canonical booleans.")
    expected_positive = (table["selected_count"] > 0) & (table["log_odds_ratio"] > 0)
    expected_overlap = table["overlap_count"] >= parameters["min_overlap"]
    expected_significant = table["p_adj_within_group"] <= float(max_p_adjusted)
    if not table["positive_enrichment"].equals(expected_positive):
        raise ValueError(f"{operation} positive-enrichment flags are invalid.")
    if not table["passes_min_overlap"].equals(expected_overlap):
        raise ValueError(f"{operation} minimum-overlap flags are invalid.")
    if not table["significant_within_group"].equals(expected_significant):
        raise ValueError(f"{operation} within-group significance flags are invalid.")
    allowed_statuses = {
        "unresolved_no_selected_markers",
        "nonpositive",
        "insufficient_overlap",
        "not_significant",
        "eligible",
        "tied_best",
    }
    for row in table.itertuples(index=False):
        if row.candidate_status not in allowed_statuses:
            raise ValueError(f"{operation} candidate status is invalid: {row.candidate_status!r}.")
        expected_status = (
            "unresolved_no_selected_markers"
            if row.selected_count == 0
            else "nonpositive"
            if not row.positive_enrichment
            else "insufficient_overlap"
            if not row.passes_min_overlap
            else "not_significant"
            if not row.significant_within_group
            else None
        )
        if expected_status is None:
            if row.candidate_status not in {"eligible", "tied_best"}:
                raise ValueError(f"{operation} candidate status disagrees with its stored gates.")
        elif row.candidate_status != expected_status:
            raise ValueError(f"{operation} candidate status disagrees with its stored gates.")
    for group in groups:
        group_frame = table.loc[table["group"] == group]
        expected = group_frame.sort_values(
            ["p_adj_within_group", "p_value", "log_odds_ratio", "overlap_count", "source"],
            ascending=[True, True, False, False, True],
            kind="mergesort",
        )
        if not group_frame.equals(expected):
            raise ValueError(f"{operation} rows do not follow canonical within-group evidence order.")
        ranks = group_frame["rank_within_group"].tolist()
        if ranks != sorted(ranks) or ranks[0] != 1:
            raise ValueError(f"{operation} within-group ranks do not follow canonical order.")
    content_fingerprint = marker_ora_table_fingerprint(table)
    if parameters["table_content_fingerprint_sha256"] != content_fingerprint:
        raise ValueError(f"{operation} table content fingerprint does not match producer evidence.")

    selected_indices = [
        int(index)
        for group in groups
        for index in table.index[table["group"] == group][:maximum]
    ]
    selected = table.loc[selected_indices]
    selected_overlap = [overlap_lists[index] for index in selected_indices]
    labels = [f"{row.group} · {row.source}" for row in selected.itertuples(index=False)]
    height = max(4.5, min(18.0, 2.5 + 0.34 * len(selected)))
    figure = Figure(figsize=(9.5, height), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    positions = np.arange(len(selected))
    eligible = selected["candidate_status"].isin(["eligible", "tied_best"]).to_numpy(dtype=bool)
    if variant == "enrichment_dot":
        adjusted = selected["p_adj_within_group"].to_numpy(dtype=float)
        color = -np.log10(np.maximum(adjusted, np.finfo(float).tiny))
        points = axis.scatter(
            selected["log_odds_ratio"].to_numpy(dtype=float),
            positions,
            s=35.0 + 55.0 * selected["overlap_count"].to_numpy(dtype=float),
            c=color,
            cmap="viridis",
            edgecolors=np.where(eligible, "#d62728", "#333333"),
            linewidths=np.where(eligible, 1.4, 0.4),
        )
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.set_xlabel("Haldane–Anscombe log odds ratio")
        axis.set_title("Marker ORA stored enrichment evidence")
        colorbar = figure.colorbar(points, ax=axis)
        colorbar.set_label("−log10(within-group adjusted p)")
    else:
        set_sizes = selected["set_size_in_universe"].to_numpy(dtype=float)
        overlaps = selected["overlap_count"].to_numpy(dtype=float)
        axis.barh(positions, set_sizes, color="#d9d9d9", label="Set size in tested universe")
        axis.barh(positions, overlaps, color=np.where(eligible, "#E45756", "#4C78A8"), label="Overlap")
        axis.set_xlabel("Genes")
        axis.set_title("Marker ORA overlap within tested resource sets")
        axis.legend()
    axis.set_yticks(positions, labels)
    axis.invert_yaxis()
    axis.grid(axis="x", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    displayed_rows = []
    for row, genes in zip(selected.to_dict(orient="records"), selected_overlap, strict=True):
        row["overlap_genes"] = genes
        displayed_rows.append(row)
    details = {
        "view": variant,
        "view_parameters": mode,
        "groups": groups,
        "plotted_sources": list(dict.fromkeys(selected["source"].tolist())),
        "available_rows": len(table),
        "plotted_rows": len(selected),
        "eligible_rows": int(selected["candidate_status"].isin(["eligible", "tied_best"]).sum()),
        "displayed_rows": displayed_rows,
        "selection_policy": "first max_terms_per_group rows in canonical upstream rank order",
        "resource_identity": {
            "sha256": parameters["resource_sha256"],
            "metadata": parameters["resource_metadata"],
        },
        "thresholds": {
            "min_overlap": parameters["min_overlap"],
            "max_p_adjusted": float(max_p_adjusted),
        },
        "table_content_fingerprint_sha256": content_fingerprint,
        "rendering": {
            "figure_size_inches": [9.5, height],
            "dpi": 120,
            "maximum_available_rows": 10_000,
            "maximum_terms_per_group": 30,
        },
    }
    return (png, details) if _return_details else png


def marker_ora_evidence_plot_code(*, producer: dict[str, Any], view: dict[str, Any]) -> str:
    fingerprint = textwrap.dedent(inspect.getsource(marker_ora_table_fingerprint)).strip()
    implementation = textwrap.dedent(inspect.getsource(_standalone_marker_ora_evidence_plot)).strip()
    return f"""from __future__ import annotations

import hashlib
import json

ORA_EVIDENCE_COLUMNS = {ORA_EVIDENCE_COLUMNS!r}

{fingerprint}


{implementation}


def plot_marker_ora_evidence(table):
    return _standalone_marker_ora_evidence_plot(
        table,
        producer={producer!r},
        view={view!r},
    )
"""


__all__ = ["_standalone_marker_ora_evidence_plot", "marker_ora_evidence_plot_code"]
