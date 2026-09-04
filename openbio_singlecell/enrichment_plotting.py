from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .plot_evidence import PLOT_EVIDENCE_SCHEMA, table_evidence_fingerprint
from .score_artifact import (
    SCORE_ARTIFACT_SCHEMA_VERSION,
    SCORE_ARTIFACTS_KEY,
    _canonical_json_sha256,
    _canonical_string_axis,
    score_frame_fingerprint,
    validate_score_artifact,
)

_TABLE_CONTRACTS = {
    "pathway_score_contrast": {
        "operation": "pathway_score_sample_welch",
        "columns": [
            "population", "pathway", "condition_a", "condition_b", "n_samples_a", "n_samples_b",
            "n_cells_a", "n_cells_b", "mean_a", "mean_b", "sd_a", "sd_b", "mean_difference",
            "ci_low", "ci_high", "t_statistic", "degrees_of_freedom", "p_value", "p_adjusted", "rank",
        ],
        "item": "pathway",
    },
    "ranked_gsea": {
        "operation": "ranked_gsea",
        "columns": [
            "comparison", "gene_set", "normalized_enrichment_score", "p_adjusted", "rank", "significant",
            "set_size_before_universe", "set_size_in_universe",
        ],
        "item": "gene_set",
    },
    "gene_set_ora": {
        "operation": "gene_set_overrepresentation",
        "columns": [
            "scope", "comparison", "gene_set", "selected_count", "universe_count", "set_size_before_universe",
            "set_size_in_universe", "overlap_count", "overlap_genes", "a", "b", "c", "d",
            "log_odds_ratio", "p_value", "p_adjusted", "rank", "positive_enrichment", "passes_min_overlap",
            "significant", "candidate", "candidate_status",
        ],
        "item": "gene_set",
    },
    "drug_ora": {
        "operation": "drug_hypergeometric",
        "columns": [
            "scope", "comparison", "drug", "selected_count", "universe_count", "set_size_before_universe",
            "set_size_in_universe", "overlap_count", "overlap_genes", "resource_sources", "a", "b", "c", "d",
            "log_odds_ratio", "p_value", "p_adjusted", "rank", "positive_enrichment", "passes_min_overlap",
            "significant", "candidate", "candidate_status",
        ],
        "item": "drug",
    },
    "drug_gsea": {
        "operation": "drug_gsea",
        "columns": [
            "comparison", "drug", "normalized_enrichment_score", "p_adjusted", "rank", "significant",
            "set_size_before_universe", "set_size_in_universe", "matched_targets", "resource_sources",
        ],
        "item": "drug",
    },
}


def _standalone_score_activity_plot(
    adata,
    *,
    score_key="aucell_scores",
    view=None,
    _return_diagnostics=False,
):
    import hashlib
    import io
    import json
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Score/Activity Plot"
    if not isinstance(score_key, str) or not score_key or score_key != score_key.strip():
        raise ValueError(f"{operation} score_key must be canonical nonblank text.")
    if not hasattr(adata, "obs") or not hasattr(adata, "obsm") or not hasattr(adata, "uns"):
        raise TypeError(f"{operation} requires an AnnData-like input.")
    if int(adata.n_obs) < 1 or not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires a nonempty unique observation axis.")
    if view is None:
        mode = {"view": "distributions", "max_scores": 20}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant == "distributions":
        if set(mode) != {"view", "max_scores"}:
            raise ValueError(f"{operation} distributions accepts only max_scores.")
        max_scores = mode["max_scores"]
        groupby = None
    elif variant == "grouped_heatmap":
        if set(mode) != {"view", "groupby", "max_scores"}:
            raise ValueError(f"{operation} grouped_heatmap requires only groupby and max_scores.")
        groupby = mode["groupby"]
        if not isinstance(groupby, str) or not groupby or groupby != groupby.strip():
            raise ValueError(f"{operation} groupby must be canonical nonblank text.")
        max_scores = mode["max_scores"]
    elif variant == "embedding":
        if set(mode) != {"view", "embedding_key", "score_name"}:
            raise ValueError(f"{operation} embedding requires only embedding_key and score_name.")
        embedding_key = mode["embedding_key"]
        score_name = mode["score_name"]
        if not isinstance(embedding_key, str) or not embedding_key or embedding_key != embedding_key.strip():
            raise ValueError(f"{operation} embedding_key must be canonical nonblank text.")
        if not isinstance(score_name, str) or score_name != score_name.strip():
            raise ValueError(f"{operation} score_name must be canonical text when supplied.")
        max_scores = 1
        groupby = None
    else:
        raise ValueError(f"{operation} view must be a closed supported selection.")
    if variant != "embedding":
        if isinstance(max_scores, (bool, np.bool_)) or not isinstance(max_scores, (int, np.integer)):
            raise TypeError(f"{operation} max_scores must be an integer.")
        if not 1 <= int(max_scores) <= 50:
            raise ValueError(f"{operation} max_scores must be between 1 and 50.")
        mode["max_scores"] = int(max_scores)

    artifacts = adata.uns.get(SCORE_ARTIFACTS_KEY)
    if isinstance(artifacts, Mapping) and score_key in artifacts:
        artifact_value = artifacts[score_key]
        storage = artifact_value.get("storage") if isinstance(artifact_value, Mapping) else None
        frame, artifact = validate_score_artifact(
            adata,
            score_key=score_key,
            required_storage=storage,
            allowed_methods=("AUCell", "GSVA", "Scanpy score_genes"),
            np=np,
            pd=pd,
        )
        producer_by_method = {
            "AUCell": "OpenBioSingleCellAUCellScores",
            "GSVA": "OpenBioSingleCellGSVAScores",
            "Scanpy score_genes": "OpenBioSingleCellGenePanelScores",
        }
        if artifact["producer_node"] != producer_by_method[artifact["method"]]:
            raise ValueError(f"{operation} score producer does not match its exact method contract.")
        method = artifact["method"]
        producer_node = artifact["producer_node"]
        evidence_fingerprint = artifact["artifact_fingerprint_sha256"]
    else:
        if score_key not in adata.obs:
            raise ValueError(f"{operation} score {score_key!r} lacks supported stored producer evidence.")
        metadata = adata.uns.get("openbio_singlecell")
        history = metadata.get("analysis_history") if isinstance(metadata, Mapping) else None
        matches = []
        if isinstance(history, Mapping):
            for key, entry in history.items():
                if not isinstance(entry, Mapping) or entry.get("operation") != "dgidb_drug_score":
                    continue
                parameters = entry.get("parameters")
                if isinstance(parameters, Mapping) and parameters.get("output_key") == score_key:
                    matches.append((str(key), entry))
        if not matches:
            raise ValueError(f"{operation} score {score_key!r} lacks exact DGIdb Drug Score provenance.")
        entry = sorted(matches, key=lambda item: item[0])[-1][1]
        parameters = entry["parameters"]
        expected = parameters.get("score_fingerprint_sha256")
        values = pd.to_numeric(adata.obs[score_key], errors="raise").to_numpy(dtype=float)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} stored drug scores must be finite.")
        payload = {
            "schema": "openbio-drug-score/v1",
            "observations": [str(value) for value in adata.obs_names],
            "scores": [float(value) for value in values],
        }
        actual = hashlib.sha256(
            json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        if not isinstance(expected, str) or actual != expected:
            raise ValueError(f"{operation} stored drug score fingerprint does not match its producer evidence.")
        frame = pd.DataFrame({score_key: values}, index=adata.obs_names.copy())
        artifact = {"score_names": [score_key], "score_content_sha256": actual}
        method = "DGIdb mean target-expression score"
        producer_node = "OpenBioSingleCellDrugScores"
        evidence_fingerprint = actual

    if variant == "embedding":
        selected_name = score_name or str(frame.columns[0])
        if selected_name not in frame.columns:
            raise ValueError(f"{operation} score_name is unavailable in the stored score family: {selected_name!r}.")
        selected_names = [selected_name]
        mode["score_name"] = selected_name
    else:
        selected_names = list(frame.columns[: int(max_scores)])
    selected = frame.loc[:, selected_names]
    values = selected.to_numpy(dtype=float, copy=False)
    if selected.index.tolist() != adata.obs_names.tolist() or not bool(np.isfinite(values).all()):
        raise ValueError(f"{operation} stored score matrix does not align to the observation axis.")
    if int(adata.n_obs) > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 cells in one static plot.")
    groups = []
    group_means = []
    if variant == "distributions":
        figure = Figure(figsize=(max(7.0, 0.45 * len(selected_names) + 3.0), 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.boxplot([values[:, index] for index in range(len(selected_names))], tick_labels=selected_names, showfliers=False)
        axis.tick_params(axis="x", labelrotation=45)
        axis.set_ylabel("Stored descriptive score")
        axis.set_title(f"{method} distributions")
        axis.grid(axis="y", alpha=0.2)
        title = f"{method} distributions"
    elif variant == "grouped_heatmap":
        if groupby not in adata.obs:
            raise ValueError(f"{operation} groupby column is unavailable: {groupby!r}.")
        labels = adata.obs[groupby]
        if bool(labels.isna().any()):
            raise ValueError(f"{operation} groupby column contains missing labels.")
        if isinstance(labels.dtype, pd.CategoricalDtype):
            observed = set(str(value) for value in labels.astype(str))
            groups = [str(value) for value in labels.cat.categories if str(value) in observed]
        else:
            groups = list(dict.fromkeys(str(value) for value in labels))
        if any(not value or value != value.strip() for value in groups):
            raise ValueError(f"{operation} group labels must be canonical nonblank text.")
        if len(groups) > 100:
            raise ValueError(f"{operation} grouped_heatmap supports at most 100 observed groups.")
        label_values = np.asarray([str(value) for value in labels])
        group_means = [values[label_values == group].mean(axis=0).tolist() for group in groups]
        figure = Figure(figsize=(max(7.0, 0.45 * len(selected_names) + 3.0), max(4.5, 0.45 * len(groups) + 2.5)), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        image = axis.imshow(np.asarray(group_means), aspect="auto", cmap="coolwarm")
        axis.set_xticks(range(len(selected_names)), selected_names, rotation=45, ha="right")
        axis.set_yticks(range(len(groups)), groups)
        axis.set_xlabel("Stored score")
        axis.set_ylabel(groupby)
        axis.set_title(f"Mean {method} by {groupby}")
        figure.colorbar(image, ax=axis, label="Arithmetic cell mean")
        title = f"{method} grouped activity"
    else:
        if embedding_key not in adata.obsm:
            raise ValueError(f"{operation} embedding is unavailable: {embedding_key!r}.")
        coordinates = np.asarray(adata.obsm[embedding_key])
        if coordinates.ndim != 2 or coordinates.shape[0] != adata.n_obs or coordinates.shape[1] < 2:
            raise ValueError(f"{operation} embedding must align to observations and contain at least two dimensions.")
        if coordinates.dtype.kind not in "iuf" or not bool(np.isfinite(coordinates[:, :2]).all()):
            raise ValueError(f"{operation} embedding coordinates must be finite real values.")
        if adata.n_obs > 200_000:
            raise ValueError(f"{operation} embedding view supports at most 200,000 cells in one static plot.")
        figure = Figure(figsize=(7.0, 6.0), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        scatter = axis.scatter(
            coordinates[:, 0], coordinates[:, 1], c=values[:, 0], cmap="viridis", s=5, linewidths=0
        )
        axis.set_xlabel(f"{embedding_key} 1")
        axis.set_ylabel(f"{embedding_key} 2")
        axis.set_title(f"{selected_names[0]} on {embedding_key}")
        figure.colorbar(scatter, ax=axis, label=f"Stored {method} score")
        title = f"{method} score embedding"
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce PNG output.")
    details = {
        "view": variant,
        "view_parameters": mode,
        "score_key": score_key,
        "method": method,
        "producer_node": producer_node,
        "score_names": selected_names,
        "available_scores": int(frame.shape[1]),
        "plotted_scores": len(selected_names),
        "plotted_cells": int(frame.shape[0]),
        "score_means": values.mean(axis=0).tolist(),
        "groups": groups,
        "group_means": group_means,
        "evidence_fingerprint_sha256": evidence_fingerprint,
        "score_content_fingerprint_sha256": artifact["score_content_sha256"],
        "title": title,
    }
    return (png, details) if _return_diagnostics else png


def run_score_activity_plot(adata: Any, *, score_key: str, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _standalone_score_activity_plot(
        adata,
        score_key=score_key,
        view=view,
        _return_diagnostics=True,
    )


def score_activity_plot_code(*, score_key: str, view: dict[str, Any]) -> str:
    helpers = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (
            _canonical_json_sha256,
            _canonical_string_axis,
            score_frame_fingerprint,
            validate_score_artifact,
            _standalone_score_activity_plot,
        )
    )
    return f'''from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any

SCORE_ARTIFACTS_KEY = {SCORE_ARTIFACTS_KEY!r}
SCORE_ARTIFACT_SCHEMA_VERSION = {SCORE_ARTIFACT_SCHEMA_VERSION!r}

{helpers}


def plot_score_activity(adata):
    return _standalone_score_activity_plot(adata, score_key={score_key!r}, view={view!r})
'''


def _standalone_enrichment_table_plot(
    table,
    *,
    producer,
    contract,
    view=None,
    _return_diagnostics=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Enrichment Evidence Plot"
    if contract not in _TABLE_CONTRACTS:
        raise ValueError(f"{operation} contract is unsupported: {contract!r}.")
    spec = _TABLE_CONTRACTS[contract]
    if not isinstance(table, pd.DataFrame):
        raise TypeError(f"{operation} requires a pandas DataFrame.")
    if table.columns.tolist() != spec["columns"]:
        raise ValueError(f"{operation} requires canonical columns {spec['columns']!r} in that order.")
    if table.empty:
        raise ValueError(f"{operation} requires nonempty stored evidence.")
    if len(table) > 100_000:
        raise ValueError(f"{operation} supports at most 100,000 evidence rows.")
    if not isinstance(producer, Mapping) or producer.get("operation") != spec["operation"]:
        raise ValueError(f"{operation} requires the exact {spec['operation']!r} producer.")
    parameters = producer.get("parameters")
    if not isinstance(parameters, Mapping) or parameters.get("plot_evidence_schema") != PLOT_EVIDENCE_SCHEMA:
        raise ValueError(f"{operation} producer lacks the exact plot evidence schema.")
    expected_fingerprint = parameters.get("plot_evidence_fingerprint_sha256")
    actual_fingerprint = table_evidence_fingerprint(table, contract=contract)
    if not isinstance(expected_fingerprint, str) or expected_fingerprint != actual_fingerprint:
        raise ValueError(f"{operation} table fingerprint does not match producer evidence.")
    item_column = spec["item"]
    items = table[item_column].tolist()
    if any(not isinstance(item, str) or not item or item != item.strip() for item in items):
        raise ValueError(f"{operation} item labels must be canonical nonblank strings.")
    if len(items) != len(set(items)) and contract not in {"pathway_score_contrast"}:
        raise ValueError(f"{operation} item labels must be unique within this result contract.")

    if contract == "pathway_score_contrast":
        default_view = {"view": "forest", "max_pathways": 30}
        allowed = {"forest", "effect_significance"}
        max_key = "max_pathways"
    elif contract in {"ranked_gsea", "drug_gsea"}:
        default_view = {"view": "nes_lollipop", ("max_drugs" if contract == "drug_gsea" else "max_gene_sets"): 30}
        allowed = {"nes_lollipop", "significance_dot"}
        max_key = "max_drugs" if contract == "drug_gsea" else "max_gene_sets"
    elif contract == "drug_ora":
        default_view = {"view": "odds_ratio_dot", "max_drugs": 30}
        allowed = {"odds_ratio_dot", "overlap_bar"}
        max_key = "max_drugs"
    else:
        default_view = {"view": "odds_ratio_dot", "max_gene_sets": 30}
        allowed = {"odds_ratio_dot", "overlap_bar"}
        max_key = "max_gene_sets"
    if view is None:
        mode = default_view
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    if set(mode) != {"view", max_key} or mode.get("view") not in allowed:
        raise ValueError(f"{operation} view must be a closed supported selection for {contract!r}.")
    limit = mode[max_key]
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)):
        raise TypeError(f"{operation} {max_key} must be an integer.")
    if not 1 <= int(limit) <= 100:
        raise ValueError(f"{operation} {max_key} must be between 1 and 100.")
    mode[max_key] = int(limit)

    def finite(column):
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} column {column!r} must contain real numeric values.")
        values = series.to_numpy(dtype=float)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} column {column!r} must be finite.")
        return values

    def integers(column, *, minimum=0):
        values = finite(column)
        if not bool(np.equal(values, np.rint(values)).all()) or bool((values < minimum).any()):
            raise ValueError(f"{operation} column {column!r} must contain integers >= {minimum}.")
        return values.astype(np.int64)

    p_adjusted = finite("p_adjusted")
    if bool(((p_adjusted < 0.0) | (p_adjusted > 1.0)).any()):
        raise ValueError(f"{operation} adjusted p-values must lie in [0, 1].")
    ranks = integers("rank", minimum=1)
    if not bool(np.all(ranks[:-1] <= ranks[1:])):
        raise ValueError(f"{operation} rank axis must be in canonical nondecreasing order.")

    if contract == "pathway_score_contrast":
        for column in ("n_samples_a", "n_samples_b", "n_cells_a", "n_cells_b"):
            integers(column, minimum=1)
        for column in ("mean_a", "mean_b", "sd_a", "sd_b", "t_statistic", "degrees_of_freedom", "p_value"):
            finite(column)
        effects = finite("mean_difference")
        low, high = finite("ci_low"), finite("ci_high")
        if bool((low > effects).any()) or bool((effects > high).any()):
            raise ValueError(f"{operation} confidence intervals must contain their mean differences.")
        if not bool(np.allclose(effects, finite("mean_a") - finite("mean_b"), rtol=1e-9, atol=1e-12)):
            raise ValueError(f"{operation} mean differences disagree with stored group means.")
        marks = {}
    elif contract in {"ranked_gsea", "drug_gsea"}:
        finite("normalized_enrichment_score")
        for column in ("set_size_before_universe", "set_size_in_universe"):
            integers(column, minimum=1)
        if table["significant"].tolist() != [bool(value <= 0.05) for value in p_adjusted]:
            raise ValueError(f"{operation} significance flags disagree with adjusted p-values.")
        marks = {}
    else:
        finite("log_odds_ratio")
        for column in ("selected_count", "universe_count", "set_size_before_universe", "set_size_in_universe", "overlap_count", "a", "b", "c", "d"):
            integers(column, minimum=0)
        if not bool(np.equal(table["a"], table["overlap_count"]).all()):
            raise ValueError(f"{operation} overlap counts disagree with contingency table a.")
        if not bool(np.equal(table["a"] + table["b"], table["set_size_in_universe"]).all()):
            raise ValueError(f"{operation} set sizes disagree with contingency evidence.")
        if not bool(np.equal(table["a"] + table["c"], table["selected_count"]).all()):
            raise ValueError(f"{operation} selected counts disagree with contingency evidence.")
        if not bool(np.equal(table["a"] + table["b"] + table["c"] + table["d"], table["universe_count"]).all()):
            raise ValueError(f"{operation} universe counts disagree with contingency evidence.")
        marks = {}

    selected = table.iloc[: int(limit)].copy()
    labels = selected[item_column].tolist()
    positions = np.arange(len(selected), dtype=float)
    if contract == "pathway_score_contrast":
        marks = {
            "mean_differences": selected["mean_difference"].astype(float).tolist(),
            "ci_low": selected["ci_low"].astype(float).tolist(),
            "ci_high": selected["ci_high"].astype(float).tolist(),
        }
    elif contract in {"ranked_gsea", "drug_gsea"}:
        marks = {
            "normalized_enrichment_scores": selected["normalized_enrichment_score"].astype(float).tolist()
        }
    else:
        marks = {
            "log_odds_ratios": selected["log_odds_ratio"].astype(float).tolist(),
            "overlap_counts": selected["overlap_count"].astype(int).tolist(),
        }
    figure = Figure(figsize=(9.0, max(4.5, 0.34 * len(selected) + 2.5)), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    variant = mode["view"]
    if contract == "pathway_score_contrast" and variant == "forest":
        x = selected["mean_difference"].to_numpy(dtype=float)
        low = selected["ci_low"].to_numpy(dtype=float)
        high = selected["ci_high"].to_numpy(dtype=float)
        axis.errorbar(x, positions, xerr=np.vstack((x - low, high - x)), fmt="o", color="#4C78A8", capsize=3)
        axis.axvline(0.0, color="#777777", linewidth=1)
        xlabel, title = "Difference of Sample means", "Pathway score Condition contrast"
    elif contract == "pathway_score_contrast":
        x = selected["mean_difference"].to_numpy(dtype=float)
        color = -np.log10(np.maximum(selected["p_adjusted"].to_numpy(dtype=float), np.finfo(float).tiny))
        scatter = axis.scatter(x, positions, c=color, cmap="viridis", s=60)
        figure.colorbar(scatter, ax=axis, label="-log10 adjusted p")
        axis.axvline(0.0, color="#777777", linewidth=1)
        xlabel, title = "Difference of Sample means", "Pathway score effect and evidence"
    elif contract in {"ranked_gsea", "drug_gsea"}:
        x = selected["normalized_enrichment_score"].to_numpy(dtype=float)
        colors = np.where(x >= 0.0, "#E45756", "#4C78A8")
        if variant == "nes_lollipop":
            axis.hlines(positions, 0.0, x, color=colors, linewidth=2)
            axis.scatter(x, positions, color=colors, s=50)
        else:
            sizes = 25.0 + 25.0 * -np.log10(np.maximum(selected["p_adjusted"].to_numpy(dtype=float), np.finfo(float).tiny))
            axis.scatter(x, positions, color=colors, s=np.minimum(sizes, 250.0))
        axis.axvline(0.0, color="#777777", linewidth=1)
        xlabel, title = "Normalized enrichment score", "Drug GSEA evidence" if contract == "drug_gsea" else "Ranked GSEA evidence"
    elif variant == "odds_ratio_dot":
        x = selected["log_odds_ratio"].to_numpy(dtype=float)
        sizes = 30.0 + 20.0 * selected["overlap_count"].to_numpy(dtype=float)
        scatter = axis.scatter(x, positions, c=-np.log10(np.maximum(selected["p_adjusted"].to_numpy(dtype=float), np.finfo(float).tiny)), cmap="viridis", s=np.minimum(sizes, 250.0))
        figure.colorbar(scatter, ax=axis, label="-log10 adjusted p")
        axis.axvline(0.0, color="#777777", linewidth=1)
        xlabel, title = "Log odds ratio", "Drug overrepresentation evidence" if contract == "drug_ora" else "Gene-set overrepresentation evidence"
    else:
        x = selected["overlap_count"].to_numpy(dtype=float)
        axis.barh(positions, x, color="#4C78A8")
        xlabel, title = "Overlap count", "Drug target overlaps" if contract == "drug_ora" else "Gene-set overlaps"
    axis.set_yticks(positions, labels)
    axis.invert_yaxis()
    axis.set_xlabel(xlabel)
    axis.set_title(title)
    axis.grid(axis="x", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    details = {
        "contract": contract,
        "view": variant,
        "view_parameters": mode,
        "items": items,
        "plotted_item_names": labels,
        "total_items": len(table),
        "plotted_items": len(selected),
        "p_adjusted": p_adjusted.tolist(),
        "ranks": ranks.tolist(),
        "result_fingerprint_sha256": actual_fingerprint,
        "title": title,
        **marks,
    }
    return (png, details) if _return_diagnostics else png


def run_enrichment_table_plot(
    table: Any,
    *,
    producer: dict[str, Any],
    contract: str,
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_enrichment_table_plot(
        table,
        producer=producer,
        contract=contract,
        view=view,
        _return_diagnostics=True,
    )


def enrichment_table_plot_code(
    *,
    function_name: str,
    producer: dict[str, Any],
    contract: str,
    view: dict[str, Any],
) -> str:
    fingerprint = textwrap.dedent(inspect.getsource(table_evidence_fingerprint)).strip()
    renderer = textwrap.dedent(inspect.getsource(_standalone_enrichment_table_plot)).strip()
    return f'''from __future__ import annotations

import hashlib
import json
import math
from typing import Any

PLOT_EVIDENCE_SCHEMA = {PLOT_EVIDENCE_SCHEMA!r}
_TABLE_CONTRACTS = {_TABLE_CONTRACTS!r}

{fingerprint}


{renderer}


def {function_name}(table):
    return _standalone_enrichment_table_plot(
        table, producer={producer!r}, contract={contract!r}, view={view!r}
    )
'''


__all__ = [
    "PLOT_EVIDENCE_SCHEMA",
    "enrichment_table_plot_code",
    "run_enrichment_table_plot",
    "run_score_activity_plot",
    "score_activity_plot_code",
    "table_evidence_fingerprint",
]
