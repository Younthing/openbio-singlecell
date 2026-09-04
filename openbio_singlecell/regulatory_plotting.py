from __future__ import annotations

import inspect
import textwrap
from collections.abc import Mapping
from typing import Any

from .plot_evidence import PLOT_EVIDENCE_SCHEMA, table_evidence_fingerprint
from .scenic_artifact import (
    SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION,
    SCENIC_BINARY_ARTIFACT_TYPE,
    SCENIC_BINARY_PRODUCER_NODE_ID,
    SCENIC_BINARY_PRODUCER_SCHEMA,
    SCENIC_MEMBERSHIP_COLUMNS,
    SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
    SCENIC_RESULT_ARTIFACT_TYPE,
    SCENIC_RESULT_PRODUCER_NODE_ID,
    SCENIC_RESULT_PRODUCER_SCHEMA,
    _activity_bytes,
    _canonical_membership,
    _result_fingerprint,
    _strict_provenance,
    validate_scenic_result_artifact,
)
from .scenic_artifact import (
    _canonical_axis as _scenic_canonical_axis,
)
from .scenic_artifact import (
    _canonical_json as _scenic_canonical_json,
)
from .tf_activity_artifact import (
    TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION,
    TF_ACTIVITY_ARTIFACT_TYPE,
    TF_ACTIVITY_PRODUCER_NODE_ID,
    TF_ACTIVITY_PRODUCER_SCHEMA,
    _activity_frame_fingerprint,
    _canonical_json_sha256,
    _portable_tf_activity_payload,
    validate_tf_activity_artifact,
)
from .tf_activity_artifact import (
    _canonical_axis as _tf_canonical_axis,
)

_REGULATORY_TABLE_CONTRACTS = {
    "tf_ranking": {
        "operation": "rank_tf_activities",
        "columns": [
            "group", "reference", "regulator", "statistic", "mean_change", "p_value", "p_adjusted",
            "direction", "significant_adjusted", "rank_within_group",
        ],
    },
    "scenic_rss": {
        "operation": "scenic_regulon_specificity",
        "columns": ["group", "regulon", "rss", "rank_within_group", "group_cells", "total_cells"],
    },
    "scenic_binarization": {
        "operation": "scenic_activity_binarization",
        "columns": ["regulon", "threshold", "threshold_source", "active_cells", "total_cells", "active_proportion"],
    },
    "scenic_membership": {
        "operation": "scenic_final_regulon_membership",
        "columns": list(SCENIC_MEMBERSHIP_COLUMNS),
    },
}


def _standalone_regulatory_activity_plot(
    adata,
    activity,
    *,
    artifact_fingerprint,
    producer_node,
    method,
    view=None,
    _return_diagnostics=False,
):
    import io

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Regulatory Activity Plot"
    if not isinstance(activity, pd.DataFrame):
        raise TypeError(f"{operation} requires a canonical activity DataFrame.")
    if activity.empty or activity.index.tolist() != adata.obs_names.tolist():
        raise ValueError(f"{operation} activity observation axis does not align exactly to AnnData.")
    values = activity.to_numpy(dtype=float, copy=False)
    if not bool(np.isfinite(values).all()):
        raise ValueError(f"{operation} activity values must be finite.")
    if int(adata.n_obs) > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 cells in one static plot.")
    if not isinstance(artifact_fingerprint, str) or len(artifact_fingerprint) != 64:
        raise ValueError(f"{operation} requires a canonical artifact fingerprint.")
    if view is None:
        mode = {"view": "distributions", "max_activities": 20}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant == "distributions":
        if set(mode) != {"view", "max_activities"}:
            raise ValueError(f"{operation} distributions accepts only max_activities.")
        groupby = None
    elif variant == "grouped_heatmap":
        if set(mode) != {"view", "groupby", "max_activities"}:
            raise ValueError(f"{operation} grouped_heatmap requires only groupby and max_activities.")
        groupby = mode["groupby"]
        if not isinstance(groupby, str) or not groupby or groupby != groupby.strip():
            raise ValueError(f"{operation} groupby must be canonical nonblank text.")
    elif variant == "embedding":
        if set(mode) != {"view", "embedding_key", "activity_name"}:
            raise ValueError(f"{operation} embedding requires only embedding_key and activity_name.")
        embedding_key = mode["embedding_key"]
        activity_name = mode["activity_name"]
        if not isinstance(embedding_key, str) or not embedding_key or embedding_key != embedding_key.strip():
            raise ValueError(f"{operation} embedding_key must be canonical nonblank text.")
        if not isinstance(activity_name, str) or activity_name != activity_name.strip():
            raise ValueError(f"{operation} activity_name must be canonical text when supplied.")
        groupby = None
    else:
        raise ValueError(f"{operation} view must be a closed supported selection.")
    if variant == "embedding":
        selected_name = activity_name or str(activity.columns[0])
        if selected_name not in activity.columns:
            raise ValueError(f"{operation} activity_name is unavailable: {selected_name!r}.")
        names = [selected_name]
        mode["activity_name"] = selected_name
    else:
        limit = mode["max_activities"]
        if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)):
            raise TypeError(f"{operation} max_activities must be an integer.")
        if not 1 <= int(limit) <= 50:
            raise ValueError(f"{operation} max_activities must be between 1 and 50.")
        mode["max_activities"] = int(limit)
        names = activity.columns[: int(limit)].tolist()
    selected = activity.loc[:, names].to_numpy(dtype=float, copy=False)
    groups = []
    group_means = []
    if variant == "distributions":
        figure = Figure(figsize=(max(7.0, 0.45 * len(names) + 3.0), 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.boxplot([selected[:, index] for index in range(len(names))], tick_labels=names, showfliers=False)
        axis.tick_params(axis="x", labelrotation=45)
        axis.set_ylabel("Stored activity")
        axis.set_title(f"{method} activity distributions")
        axis.grid(axis="y", alpha=0.2)
        title = f"{method} activity distributions"
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
        group_means = [selected[label_values == group].mean(axis=0).tolist() for group in groups]
        figure = Figure(figsize=(max(7.0, 0.45 * len(names) + 3.0), max(4.5, 0.45 * len(groups) + 2.5)), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        image = axis.imshow(np.asarray(group_means), aspect="auto", cmap="coolwarm")
        axis.set_xticks(range(len(names)), names, rotation=45, ha="right")
        axis.set_yticks(range(len(groups)), groups)
        axis.set_xlabel("Stored activity")
        axis.set_ylabel(groupby)
        axis.set_title(f"Mean {method} activity by {groupby}")
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
            coordinates[:, 0], coordinates[:, 1], c=selected[:, 0], cmap="viridis", s=5, linewidths=0
        )
        axis.set_xlabel(f"{embedding_key} 1")
        axis.set_ylabel(f"{embedding_key} 2")
        axis.set_title(f"{names[0]} on {embedding_key}")
        figure.colorbar(scatter, ax=axis, label=f"Stored {method} activity")
        title = f"{method} activity embedding"
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    details = {
        "view": variant,
        "view_parameters": mode,
        "producer_node": producer_node,
        "method": method,
        "activity_names": names,
        "available_activities": int(activity.shape[1]),
        "plotted_activities": len(names),
        "plotted_cells": int(activity.shape[0]),
        "activity_means": selected.mean(axis=0).tolist(),
        "groups": groups,
        "group_means": group_means,
        "artifact_fingerprint_sha256": artifact_fingerprint,
        "title": title,
    }
    return (png, details) if _return_diagnostics else png


def run_tf_activity_plot(adata: Any, activities: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    scores, _adjusted, _provenance, metadata = validate_tf_activity_artifact(
        activities,
        exact_type=not isinstance(activities, Mapping),
        copy_result=False,
    )
    return _standalone_regulatory_activity_plot(
        adata,
        scores,
        artifact_fingerprint=metadata["artifact_fingerprint_sha256"],
        producer_node=metadata["producer_node_id"],
        method="CollecTRI ULM",
        view=view,
        _return_diagnostics=True,
    )


def run_scenic_activity_plot(adata: Any, scenic_result: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    import numpy as np
    import pandas as pd

    activity, _membership, _provenance, metadata = validate_scenic_result_artifact(
        scenic_result,
        exact_type=not isinstance(scenic_result, dict),
        numpy=np,
        pandas=pd,
        copy_result=False,
    )
    return _standalone_regulatory_activity_plot(
        adata,
        activity,
        artifact_fingerprint=metadata["artifact_fingerprint_sha256"],
        producer_node=SCENIC_RESULT_PRODUCER_NODE_ID,
        method="SCENIC AUCell",
        view=view,
        _return_diagnostics=True,
    )


def tf_activity_plot_code(*, view: dict[str, Any]) -> str:
    validator_parts = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (
            _canonical_json_sha256,
            _tf_canonical_axis,
            _activity_frame_fingerprint,
            _portable_tf_activity_payload,
            validate_tf_activity_artifact,
        )
    )
    renderer = textwrap.dedent(inspect.getsource(_standalone_regulatory_activity_plot)).strip()
    return f'''from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any

TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION = {TF_ACTIVITY_ARTIFACT_SCHEMA_VERSION!r}
TF_ACTIVITY_ARTIFACT_TYPE = {TF_ACTIVITY_ARTIFACT_TYPE!r}
TF_ACTIVITY_PRODUCER_NODE_ID = {TF_ACTIVITY_PRODUCER_NODE_ID!r}
TF_ACTIVITY_PRODUCER_SCHEMA = {TF_ACTIVITY_PRODUCER_SCHEMA!r}

{validator_parts}


{renderer}


def plot_tf_activity(adata, activities):
    scores, _adjusted, _provenance, metadata = validate_tf_activity_artifact(
        activities, exact_type=False, copy_result=False
    )
    return _standalone_regulatory_activity_plot(
        adata,
        scores,
        artifact_fingerprint=metadata["artifact_fingerprint_sha256"],
        producer_node=metadata["producer_node_id"],
        method="CollecTRI ULM",
        view={view!r},
    )
'''


def scenic_activity_plot_code(*, view: dict[str, Any]) -> str:
    validator_parts = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (
            _scenic_canonical_json,
            _scenic_canonical_axis,
            _activity_bytes,
            _canonical_membership,
            _strict_provenance,
            _result_fingerprint,
            validate_scenic_result_artifact,
        )
    )
    renderer = textwrap.dedent(inspect.getsource(_standalone_regulatory_activity_plot)).strip()
    return f'''from __future__ import annotations

import copy
import hashlib
import json

import numpy as np
import pandas as pd

SCENIC_RESULT_ARTIFACT_TYPE = {SCENIC_RESULT_ARTIFACT_TYPE!r}
SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION = {SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION!r}
SCENIC_RESULT_PRODUCER_NODE_ID = {SCENIC_RESULT_PRODUCER_NODE_ID!r}
SCENIC_RESULT_PRODUCER_SCHEMA = {SCENIC_RESULT_PRODUCER_SCHEMA!r}
SCENIC_BINARY_ARTIFACT_TYPE = {SCENIC_BINARY_ARTIFACT_TYPE!r}
SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION = {SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION!r}
SCENIC_BINARY_PRODUCER_NODE_ID = {SCENIC_BINARY_PRODUCER_NODE_ID!r}
SCENIC_BINARY_PRODUCER_SCHEMA = {SCENIC_BINARY_PRODUCER_SCHEMA!r}
SCENIC_MEMBERSHIP_COLUMNS = {SCENIC_MEMBERSHIP_COLUMNS!r}

class SCENICResultArtifact:
    pass

{validator_parts}


{renderer}


def plot_scenic_activity(adata, scenic_result):
    if not isinstance(scenic_result, dict) and callable(getattr(scenic_result, "to_portable", None)):
        scenic_result = scenic_result.to_portable()
    activity, _membership, _provenance, metadata = validate_scenic_result_artifact(
        scenic_result, exact_type=False, numpy=np, pandas=pd, copy_result=False
    )
    return _standalone_regulatory_activity_plot(
        adata,
        activity,
        artifact_fingerprint=metadata["artifact_fingerprint_sha256"],
        producer_node=SCENIC_RESULT_PRODUCER_NODE_ID,
        method="SCENIC AUCell",
        view={view!r},
    )
'''


def _standalone_regulatory_table_plot(
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

    operation = "Regulatory Evidence Plot"
    if contract not in _REGULATORY_TABLE_CONTRACTS:
        raise ValueError(f"{operation} contract is unsupported: {contract!r}.")
    spec = _REGULATORY_TABLE_CONTRACTS[contract]
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
    actual_fingerprint = table_evidence_fingerprint(table, contract=contract)
    if parameters.get("plot_evidence_fingerprint_sha256") != actual_fingerprint:
        raise ValueError(f"{operation} table fingerprint does not match producer evidence.")

    def finite(column):
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} column {column!r} must be real numeric.")
        values = series.to_numpy(dtype=float)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} column {column!r} must be finite.")
        return values

    def integers(column, minimum=0):
        values = finite(column)
        if not bool(np.equal(values, np.rint(values)).all()) or bool((values < minimum).any()):
            raise ValueError(f"{operation} column {column!r} must contain integers >= {minimum}.")
        return values.astype(np.int64)

    if contract == "tf_ranking":
        default = {"view": "ranked_dot", "max_per_group": 10}
        allowed = {"ranked_dot", "effect_significance"}
        max_key = "max_per_group"
    elif contract == "scenic_rss":
        default = {"view": "rss_heatmap", "max_regulons": 20}
        allowed = {"rss_heatmap", "top_regulons"}
        max_key = "max_regulons"
    elif contract == "scenic_binarization":
        default = {"view": "thresholds", "max_regulons": 30}
        allowed = {"thresholds", "active_fraction"}
        max_key = "max_regulons"
    else:
        default = {"view": "target_count", "max_regulons": 30}
        allowed = {"target_count", "top_targets"}
        max_key = "max_regulons" if (view or default).get("view") == "target_count" else "max_targets"
    if view is None:
        mode = default
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    expected_fields = {"view", max_key}
    if contract == "scenic_membership" and variant == "top_targets":
        expected_fields.add("regulon")
    if set(mode) != expected_fields or variant not in allowed:
        raise ValueError(f"{operation} view must be a closed supported selection for {contract!r}.")
    limit = mode[max_key]
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)):
        raise TypeError(f"{operation} {max_key} must be an integer.")
    if not 1 <= int(limit) <= 100:
        raise ValueError(f"{operation} {max_key} must be between 1 and 100.")
    mode[max_key] = int(limit)

    marks = {}
    plotted_marks = 0
    if contract == "tf_ranking":
        for column in ("group", "reference", "regulator", "direction"):
            if any(not isinstance(value, str) or not value or value != value.strip() for value in table[column]):
                raise ValueError(f"{operation} {column} labels must be canonical nonblank strings.")
        statistics, changes = finite("statistic"), finite("mean_change")
        p_adjusted = finite("p_adjusted")
        if bool(((p_adjusted < 0) | (p_adjusted > 1)).any()):
            raise ValueError(f"{operation} adjusted p-values must lie in [0, 1].")
        ranks = integers("rank_within_group", 1)
        for group in dict.fromkeys(table["group"]):
            group_ranks = ranks[table["group"].to_numpy() == group]
            if not bool(np.all(group_ranks[:-1] <= group_ranks[1:])):
                raise ValueError(f"{operation} rank order must be canonical within every group.")
        threshold = float(parameters.get("report_p_adjusted", 0.05))
        if table["significant_adjusted"].tolist() != [bool(value <= threshold) for value in p_adjusted]:
            raise ValueError(f"{operation} significance flags disagree with adjusted p-values.")
        groups = list(dict.fromkeys(table["group"].tolist()))
        if len(groups) > 100:
            raise ValueError(f"{operation} supports at most 100 groups in one static plot.")
        group_sizes = table.groupby("group", sort=False, observed=True).size().reindex(groups)
        planned_marks = sum(min(int(size), int(limit)) for size in group_sizes)
        if planned_marks > 100:
            raise ValueError(f"{operation} supports at most 100 plotted marks in one static plot.")
        selected = table.groupby("group", sort=False, observed=True).head(int(limit))
        plotted_marks = len(selected)
        if plotted_marks != planned_marks:
            raise RuntimeError(f"{operation} selected-row count disagrees with its plot preflight.")
        labels = [f"{row.group}: {row.regulator}" for row in selected.itertuples(index=False)]
        positions = np.arange(len(selected))
        x = selected["statistic" if variant == "ranked_dot" else "mean_change"].to_numpy(dtype=float)
        figure = Figure(figsize=(9.0, max(4.5, 0.34 * len(selected) + 2.5)), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        colors = np.where(x >= 0, "#E45756", "#4C78A8")
        sizes = 45.0 if variant == "ranked_dot" else np.minimum(250.0, 25.0 + 25.0 * -np.log10(np.maximum(selected["p_adjusted"], np.finfo(float).tiny)))
        axis.scatter(x, positions, color=colors, s=sizes)
        axis.axvline(0.0, color="#777777", linewidth=1)
        axis.set_xlabel("Test statistic" if variant == "ranked_dot" else "Arithmetic mean change")
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.set_title("Annotation-associated TF activity evidence")
        title = "TF activity ranking evidence"
        marks = {"statistics": statistics.tolist(), "mean_changes": changes.tolist(), "p_adjusted": p_adjusted.tolist()}
    elif contract == "scenic_rss":
        rss = finite("rss")
        if bool(((rss < 0) | (rss > 1)).any()):
            raise ValueError(f"{operation} RSS values must lie in [0, 1].")
        integers("rank_within_group", 1)
        integers("group_cells", 1)
        integers("total_cells", 1)
        for column in ("group", "regulon"):
            if any(not isinstance(value, str) or not value or value != value.strip() for value in table[column]):
                raise ValueError(f"{operation} {column} labels must be canonical nonblank strings.")
        if bool(table[["group", "regulon"]].duplicated().any()):
            raise ValueError(f"{operation} RSS group-regulon pairs must be unique.")
        groups = list(dict.fromkeys(table["group"].tolist()))
        if len(groups) > 100:
            raise ValueError(f"{operation} supports at most 100 groups in one static plot.")
        regulons = list(dict.fromkeys(table["regulon"].tolist()))[: int(limit)]
        expected_pairs = {(group, regulon) for group in groups for regulon in list(dict.fromkeys(table["regulon"]))}
        if set(zip(table["group"], table["regulon"], strict=True)) != expected_pairs:
            raise ValueError(f"{operation} RSS table must contain the complete group-by-regulon grid.")
        pivot = table.pivot(index="group", columns="regulon", values="rss").reindex(index=groups, columns=regulons)
        figure = Figure(figsize=(max(7.0, 0.5 * len(regulons) + 3.0), max(4.5, 0.45 * len(groups) + 2.5)), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        if variant == "rss_heatmap":
            image = axis.imshow(pivot.to_numpy(dtype=float), aspect="auto", vmin=0, vmax=1, cmap="viridis")
            axis.set_xticks(range(len(regulons)), regulons, rotation=45, ha="right")
            axis.set_yticks(range(len(groups)), groups)
            figure.colorbar(image, ax=axis, label="Regulon specificity score")
            plotted_marks = int(pivot.size)
        else:
            selected = table.loc[table["regulon"].isin(regulons)].groupby("group", sort=False, observed=True).head(5)
            labels = [f"{row.group}: {row.regulon}" for row in selected.itertuples(index=False)]
            axis.barh(range(len(selected)), selected["rss"], color="#4C78A8")
            axis.set_yticks(range(len(selected)), labels)
            axis.invert_yaxis()
            axis.set_xlabel("Regulon specificity score")
            plotted_marks = len(selected)
        axis.set_title("SCENIC regulon specificity")
        title = "SCENIC regulon specificity"
        marks = {"rss_values": rss.tolist(), "groups": groups, "regulons": regulons}
    elif contract == "scenic_binarization":
        thresholds, proportions = finite("threshold"), finite("active_proportion")
        active, totals = integers("active_cells", 0), integers("total_cells", 1)
        if bool(((proportions < 0) | (proportions > 1)).any()) or not bool(np.allclose(proportions, active / totals)):
            raise ValueError(f"{operation} active proportions disagree with stored cell counts.")
        if bool(table["regulon"].duplicated().any()):
            raise ValueError(f"{operation} regulon axis must be unique.")
        if any(
            value not in {"pyscenic_0.12.1_hdt", "manual_override"}
            for value in table["threshold_source"]
        ):
            raise ValueError(f"{operation} threshold sources are outside the canonical producer contract.")
        selected = table.iloc[: int(limit)]
        plotted_marks = len(selected)
        labels = selected["regulon"].tolist()
        x = selected["threshold" if variant == "thresholds" else "active_proportion"].to_numpy(dtype=float)
        figure = Figure(figsize=(9.0, max(4.5, 0.34 * len(selected) + 2.5)), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.barh(range(len(selected)), x, color="#4C78A8")
        axis.set_yticks(range(len(selected)), labels)
        axis.invert_yaxis()
        axis.set_xlabel("AUC threshold" if variant == "thresholds" else "Active-cell fraction")
        axis.set_title("SCENIC activity binarization evidence")
        title = "SCENIC activity binarization evidence"
        marks = {"thresholds": thresholds.tolist(), "active_proportions": proportions.tolist(), "active_cells": active.tolist()}
    else:
        if bool(table[["regulon", "target"]].duplicated().any()):
            raise ValueError(f"{operation} regulon-target edges must be unique.")
        weights = finite("target_weight")
        motif_counts = integers("motif_evidence_count", 1)
        target_ranks = integers("target_rank", 1)
        for position, row in enumerate(table.itertuples(index=False)):
            for value in (row.regulon, row.transcription_factor, row.regulation, row.context, row.target):
                if not isinstance(value, str) or not value or value != value.strip():
                    raise ValueError(f"{operation} membership labels must be canonical nonblank strings.")
            if row.regulation not in {"activating", "repressing"}:
                raise ValueError(f"{operation} regulation must be activating or repressing.")
            suffix = "(+)" if row.regulation == "activating" else "(-)"
            if row.regulon != f"{row.transcription_factor}{suffix}":
                raise ValueError(f"{operation} regulon identity disagrees with TF and regulation.")
            if not isinstance(row.motif_ids, (list, tuple)) or len(row.motif_ids) != motif_counts[position]:
                raise ValueError(f"{operation} motif evidence count disagrees with motif_ids.")
        for regulon in dict.fromkeys(table["regulon"]):
            ranks = target_ranks[table["regulon"].to_numpy() == regulon]
            if not bool(np.all(ranks[:-1] <= ranks[1:])):
                raise ValueError(f"{operation} target ranks must be canonical within each regulon.")
        regulons = list(dict.fromkeys(table["regulon"].tolist()))
        counts = table.groupby("regulon", sort=False, observed=True)["target"].nunique().reindex(regulons).astype(int)
        if variant == "target_count":
            labels = regulons[: int(limit)]
            x = counts.loc[labels].to_numpy(dtype=float)
            xlabel = "Unique targets"
            plotted_marks = len(labels)
        else:
            regulon = mode["regulon"]
            if not isinstance(regulon, str) or not regulon or regulon != regulon.strip() or regulon not in regulons:
                raise ValueError(f"{operation} top_targets requires one exact stored regulon.")
            selected = table.loc[table["regulon"] == regulon].head(int(limit))
            labels = selected["target"].tolist()
            x = selected["target_weight"].to_numpy(dtype=float)
            xlabel = "Stored target weight"
            plotted_marks = len(labels)
        figure = Figure(figsize=(9.0, max(4.5, 0.34 * len(labels) + 2.5)), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.barh(range(len(labels)), x, color="#4C78A8")
        axis.set_yticks(range(len(labels)), labels)
        axis.invert_yaxis()
        axis.set_xlabel(xlabel)
        axis.set_title("SCENIC final regulon membership")
        title = "SCENIC final regulon membership"
        marks = {"regulons": regulons, "target_counts": counts.tolist(), "target_weights": weights.tolist()}
    axis.grid(axis="x", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    details = {
        "contract": contract,
        "view": variant,
        "view_parameters": mode,
        "total_rows": len(table),
        "plotted_marks": int(plotted_marks),
        "result_fingerprint_sha256": actual_fingerprint,
        "title": title,
        **marks,
    }
    return (png, details) if _return_diagnostics else png


def run_regulatory_table_plot(table: Any, *, producer: dict[str, Any], contract: str, view: Any = None):
    return _standalone_regulatory_table_plot(
        table,
        producer=producer,
        contract=contract,
        view=view,
        _return_diagnostics=True,
    )


def regulatory_table_plot_code(*, function_name: str, producer: dict[str, Any], contract: str, view: dict[str, Any]):
    fingerprint = textwrap.dedent(inspect.getsource(table_evidence_fingerprint)).strip()
    renderer = textwrap.dedent(inspect.getsource(_standalone_regulatory_table_plot)).strip()
    return f'''from __future__ import annotations

import hashlib
import json
import math
from typing import Any

PLOT_EVIDENCE_SCHEMA = {PLOT_EVIDENCE_SCHEMA!r}
_REGULATORY_TABLE_CONTRACTS = {_REGULATORY_TABLE_CONTRACTS!r}

{fingerprint}


{renderer}


def {function_name}(table):
    return _standalone_regulatory_table_plot(
        table, producer={producer!r}, contract={contract!r}, view={view!r}
    )
'''


__all__ = [
    "regulatory_table_plot_code",
    "run_regulatory_table_plot",
    "run_scenic_activity_plot",
    "run_tf_activity_plot",
    "scenic_activity_plot_code",
    "table_evidence_fingerprint",
    "tf_activity_plot_code",
]
