from __future__ import annotations

import inspect
import textwrap
from typing import Any


def _standalone_mad_outlier_plot(adata, *, view=None, _return_diagnostics=False):
    import hashlib
    import io
    from collections.abc import Mapping

    import matplotlib
    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from scipy.stats import median_abs_deviation

    operation = "MAD Outlier Plot"
    if view is None:
        mode = {"view": "metric_distributions"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    if set(mode) != {"view"} or mode.get("view") not in {
        "metric_distributions",
        "sample_marked_fraction",
    }:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if not isinstance(adata, AnnData):
        raise TypeError(f"{operation} requires an in-memory AnnData object.")
    if adata.isbacked:
        raise ValueError(f"{operation} does not support backed AnnData; load it into memory first.")
    if adata.n_obs < 1:
        raise ValueError(f"{operation} requires at least one observation.")
    # ponytail: fixed static-PNG ceiling; add sampled/paginated rendering if million-cell plots become necessary.
    if adata.n_obs > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 observations in one static PNG.")

    evidence = adata.uns.get("openbio_mad_outliers")
    expected_evidence = {
        "schema_version",
        "producer_operation",
        "observation_axis_fingerprint_sha256",
        "metrics",
        "batch_key",
        "output_column",
        "nmads",
        "direction",
        "scale_mad",
        "minimum_group_size",
        "thresholds",
        "group_rates",
    }
    if not isinstance(evidence, Mapping) or set(evidence) != expected_evidence:
        raise ValueError(f"{operation} requires the exact stored Mark MAD Outliers diagnostic evidence schema.")
    if evidence["schema_version"] != 1 or evidence["producer_operation"] != "mark_mad_outliers":
        raise ValueError(f"{operation} evidence was not produced by Mark MAD Outliers.")

    def axis_fingerprint(names):
        digest = hashlib.sha256()
        for value in names:
            encoded = str(value).encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "little"))
            digest.update(encoded)
        return f"sha256:{digest.hexdigest()}"

    if evidence["observation_axis_fingerprint_sha256"] != axis_fingerprint(adata.obs_names):
        raise ValueError(f"{operation} observation axis fingerprint does not match the stored evidence.")
    metadata = adata.uns.get("openbio_singlecell")
    history = metadata.get("analysis_history") if isinstance(metadata, Mapping) else None
    if not isinstance(history, Mapping):
        raise ValueError(f"{operation} requires OpenBio analysis history from Mark MAD Outliers.")

    metrics_value = evidence["metrics"]
    if hasattr(metrics_value, "tolist"):
        metrics_value = metrics_value.tolist()
    if not isinstance(metrics_value, (list, tuple)) or not metrics_value:
        raise ValueError(f"{operation} stored metrics must be a non-empty ordered string axis.")
    metrics = list(metrics_value)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in metrics):
        raise ValueError(f"{operation} stored metrics must be canonical nonblank strings.")
    if len(set(metrics)) != len(metrics):
        raise ValueError(f"{operation} stored metric axis contains duplicates.")
    if len(metrics) > 8:
        raise ValueError(f"{operation} supports at most 8 metrics in one static PNG; found {len(metrics)}.")
    batch_key = evidence["batch_key"]
    output_column = evidence["output_column"]
    direction = evidence["direction"]
    scale_mad = evidence["scale_mad"]
    minimum_group_size = evidence["minimum_group_size"]
    nmads = evidence["nmads"]
    if not isinstance(batch_key, str) or not isinstance(output_column, str) or not output_column:
        raise ValueError(f"{operation} stored column identifiers are invalid.")
    if direction not in {"both", "upper", "lower"} or not isinstance(scale_mad, (bool, np.bool_)):
        raise ValueError(f"{operation} stored MAD rule is invalid.")
    if isinstance(nmads, (bool, np.bool_)) or not isinstance(nmads, (int, float)) or not np.isfinite(nmads) or nmads <= 0:
        raise ValueError(f"{operation} stored nmads must be finite and positive.")
    if isinstance(minimum_group_size, (bool, np.bool_)) or not isinstance(minimum_group_size, (int, np.integer)) or minimum_group_size < 1:
        raise ValueError(f"{operation} stored minimum_group_size must be a positive integer.")
    matching_history = []
    for entry in history.values():
        parameters = entry.get("parameters") if isinstance(entry, Mapping) else None
        if not isinstance(parameters, Mapping) or entry.get("operation") != "mark_mad_outliers":
            continue
        history_metrics = parameters.get("metrics")
        if hasattr(history_metrics, "tolist"):
            history_metrics = history_metrics.tolist()
        if (
            list(history_metrics) if isinstance(history_metrics, (list, tuple)) else None
        ) == metrics and all(
            parameters.get(key) == value
            for key, value in {
                "batch_key": batch_key,
                "output_column": output_column,
                "direction": direction,
                "scale_mad": bool(scale_mad),
                "minimum_group_size": int(minimum_group_size),
                "nmads": float(nmads),
            }.items()
        ):
            matching_history.append(entry)
    if not matching_history:
        raise ValueError(f"{operation} evidence does not match a Mark MAD Outliers provenance entry.")

    missing_metrics = [metric for metric in metrics if metric not in adata.obs]
    if missing_metrics:
        raise ValueError(f"{operation} metric columns are missing from obs: {missing_metrics!r}.")
    if output_column not in adata.obs:
        raise ValueError(f"{operation} marked-cell column is missing from obs: {output_column!r}.")
    marked = adata.obs[output_column]
    if not pd.api.types.is_bool_dtype(marked.dtype) or bool(marked.isna().any()):
        raise TypeError(f"{operation} marked-cell column must be complete boolean evidence.")
    if batch_key:
        if batch_key not in adata.obs:
            raise ValueError(f"{operation} Sample column is missing from obs: {batch_key!r}.")
        if bool(adata.obs[batch_key].isna().any()):
            raise ValueError(f"{operation} Sample column contains missing labels.")
        rendered = [str(value) for value in pd.unique(adata.obs[batch_key])]
        if len(set(rendered)) != len(rendered):
            raise ValueError(f"{operation} Sample labels collide after string rendering.")

    thresholds = evidence["thresholds"]
    rates = evidence["group_rates"]
    threshold_columns = [
        "metric",
        "group",
        "status",
        "n",
        "finite",
        "missing",
        "median",
        "mad",
        "lower_threshold",
        "upper_threshold",
        "flagged",
    ]
    rate_columns = ["group", "cells", "marked", "marked_fraction"]
    if not isinstance(thresholds, pd.DataFrame) or thresholds.columns.tolist() != threshold_columns:
        raise ValueError(f"{operation} requires the exact stored MAD threshold table schema.")
    if not isinstance(rates, pd.DataFrame) or rates.columns.tolist() != rate_columns:
        raise ValueError(f"{operation} requires the exact stored MAD group-rate table schema.")
    if rates.empty:
        raise ValueError(f"{operation} requires at least one stored Sample/group rate.")
    groups = rates["group"].tolist()
    if any(not isinstance(value, str) or not value or value != value.strip() for value in groups):
        raise ValueError(f"{operation} stored Sample/group axis must contain canonical strings.")
    if len(set(groups)) != len(groups):
        raise ValueError(f"{operation} stored Sample/group axis contains duplicates.")
    if len(groups) > 100:
        raise ValueError(f"{operation} supports at most 100 Samples/groups in one static PNG; found {len(groups)}.")
    expected_pairs = [(metric, group) for metric in metrics for group in groups]
    if list(thresholds[["metric", "group"]].itertuples(index=False, name=None)) != expected_pairs:
        raise ValueError(f"{operation} threshold rows must form the canonical metric-by-group grid.")

    if batch_key:
        group_masks = {
            group: np.asarray([str(value) == group for value in adata.obs[batch_key]], dtype=bool)
            for group in groups
        }
        if not bool(np.logical_or.reduce(list(group_masks.values())).all()) or any(
            int(mask.sum()) == 0 for mask in group_masks.values()
        ):
            raise ValueError(f"{operation} stored Sample/group axis does not partition the observations.")
    else:
        if groups != ["all_cells"]:
            raise ValueError(f"{operation} global evidence must use the single 'all_cells' group.")
        group_masks = {"all_cells": np.ones(adata.n_obs, dtype=bool)}

    union = np.zeros(adata.n_obs, dtype=bool)
    missing_values = 0
    threshold_records = []
    for row in thresholds.itertuples(index=False):
        values = adata.obs[row.metric]
        if pd.api.types.is_bool_dtype(values.dtype) or not pd.api.types.is_numeric_dtype(values.dtype):
            raise TypeError(f"{operation} metric {row.metric!r} must remain real numeric.")
        array = values.to_numpy(dtype=float, na_value=np.nan)
        selected = array[group_masks[row.group]]
        finite = np.isfinite(selected)
        missing_values += int((~finite).sum())
        integer_values = (row.n, row.finite, row.missing, row.flagged)
        if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0 for value in integer_values):
            raise ValueError(f"{operation} threshold count evidence must contain nonnegative integers.")
        if (int(row.n), int(row.finite), int(row.missing)) != (
            int(selected.size),
            int(finite.sum()),
            int((~finite).sum()),
        ):
            raise ValueError(f"{operation} threshold count evidence does not match the stored observations.")
        if row.status not in {"evaluated", "zero_mad", "skipped_insufficient_finite_values"}:
            raise ValueError(f"{operation} threshold status is invalid.")
        numbers = np.asarray([row.median, row.mad, row.lower_threshold, row.upper_threshold], dtype=float)
        if row.status == "skipped_insufficient_finite_values":
            if int(finite.sum()) >= int(minimum_group_size) or not bool(np.isnan(numbers).all()) or row.flagged != 0:
                raise ValueError(f"{operation} skipped-threshold evidence is inconsistent.")
            group_flags = np.zeros(selected.size, dtype=bool)
        else:
            if int(finite.sum()) < int(minimum_group_size) or not bool(np.isfinite(numbers[:2]).all()) or row.mad < 0:
                raise ValueError(f"{operation} evaluated-threshold evidence is inconsistent.")
            observed_median = float(np.median(selected[finite]))
            observed_mad = float(median_abs_deviation(selected[finite], scale="normal" if scale_mad else 1.0))
            if not np.isclose(row.median, observed_median, rtol=0.0, atol=1e-12) or not np.isclose(
                row.mad, observed_mad, rtol=0.0, atol=1e-12
            ):
                raise ValueError(f"{operation} stored threshold statistics do not match the metric values.")
            expected_lower = np.nan if direction == "upper" else observed_median - float(nmads) * observed_mad
            expected_upper = np.nan if direction == "lower" else observed_median + float(nmads) * observed_mad
            for stored, expected in ((row.lower_threshold, expected_lower), (row.upper_threshold, expected_upper)):
                if not (np.isnan(stored) and np.isnan(expected)) and not np.isclose(
                    stored, expected, rtol=0.0, atol=1e-12
                ):
                    raise ValueError(f"{operation} stored thresholds do not match the declared MAD rule.")
            lower = -np.inf if np.isnan(row.lower_threshold) else float(row.lower_threshold)
            upper = np.inf if np.isnan(row.upper_threshold) else float(row.upper_threshold)
            group_flags = finite & ((selected < lower) | (selected > upper))
            if int(group_flags.sum()) != int(row.flagged):
                raise ValueError(f"{operation} stored threshold flags do not match the metric values.")
        positions = np.flatnonzero(group_masks[row.group])
        union[positions[group_flags]] = True
        threshold_records.append(
            {
                key: (None if isinstance(value, float) and np.isnan(value) else value)
                for key, value in row._asdict().items()
            }
        )
    if not np.array_equal(union, marked.to_numpy(dtype=bool)):
        raise ValueError(f"{operation} marked-cell column does not equal the stored metric-threshold union.")
    rate_records = []
    for row in rates.itertuples(index=False):
        mask = group_masks[row.group]
        cells = int(mask.sum())
        count = int(union[mask].sum())
        fraction = count / cells
        if (
            isinstance(row.cells, (bool, np.bool_))
            or not isinstance(row.cells, (int, np.integer))
            or isinstance(row.marked, (bool, np.bool_))
            or not isinstance(row.marked, (int, np.integer))
            or not isinstance(row.marked_fraction, (int, float, np.integer, np.floating))
            or not np.isfinite(row.marked_fraction)
            or (int(row.cells), int(row.marked)) != (cells, count)
            or not np.isclose(row.marked_fraction, fraction, rtol=0.0, atol=1e-12)
        ):
            raise ValueError(f"{operation} stored Sample/group marked fractions do not match the observations.")
        rate_records.append(
            {"group": row.group, "cells": cells, "marked": count, "marked_fraction": float(fraction)}
        )

    if mode["view"] == "sample_marked_fraction" and not batch_key:
        raise ValueError(f"{operation} Sample marked fraction requires a non-empty producer batch_key.")
    if mode["view"] == "metric_distributions":
        if len(groups) > 12:
            raise ValueError(
                f"{operation} metric distributions support at most 12 Samples/groups in one static PNG; "
                f"found {len(groups)}."
            )
        figure = Figure(figsize=(9.0, min(20.0, 3.2 * len(metrics))), constrained_layout=True)
        FigureCanvasAgg(figure)
        axes = np.atleast_1d(figure.subplots(len(metrics), 1, squeeze=False)).reshape(-1)
        palette = matplotlib.colormaps["tab20"].resampled(max(1, len(groups)))
        for metric_index, (metric, axis) in enumerate(zip(metrics, axes, strict=True)):
            values = adata.obs[metric].to_numpy(dtype=float, na_value=np.nan)
            threshold_label_added = False
            for group_index, group in enumerate(groups):
                selected = values[group_masks[group]]
                finite = selected[np.isfinite(selected)]
                if finite.size:
                    axis.hist(
                        finite,
                        bins=min(30, max(5, int(np.ceil(np.sqrt(finite.size))))),
                        histtype="step",
                        linewidth=1.5,
                        color=palette(group_index),
                        label=group,
                    )
                row = thresholds.loc[(thresholds["metric"] == metric) & (thresholds["group"] == group)].iloc[0]
                for column in ("lower_threshold", "upper_threshold"):
                    threshold = float(row[column])
                    if np.isfinite(threshold):
                        label = "Stored MAD threshold" if metric_index == 0 and not threshold_label_added else None
                        axis.axvline(
                            threshold,
                            color=palette(group_index),
                            linestyle="--",
                            linewidth=1.0,
                            label=label,
                        )
                        threshold_label_added = True
            axis.set_title(metric)
            axis.set_xlabel("Metric value")
            axis.set_ylabel("Cells")
            axis.grid(axis="y", alpha=0.2)
            if metric_index == 0:
                axis.legend(title="Sample/group", bbox_to_anchor=(1.02, 1.0), loc="upper left")
        title = "MAD metric distributions and stored thresholds"
    else:
        figure = Figure(figsize=(max(6.4, min(20.0, 0.55 * len(groups) + 2.5)), 5.0))
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        positions = np.arange(len(groups), dtype=float)
        fractions = np.asarray([record["marked_fraction"] for record in rate_records], dtype=float)
        axis.bar(positions, fractions, color="#4C78A8")
        axis.set_xticks(positions, groups, rotation=45, ha="right")
        axis.set_ylim(0.0, 1.0)
        axis.set_ylabel("Marked-cell fraction")
        axis.set_title("MAD-marked fraction by Sample")
        axis.grid(axis="y", alpha=0.2)
        title = "MAD-marked fraction by Sample"
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")

    details = {
        "view": mode["view"],
        "view_parameters": mode,
        "metrics": metrics,
        "groups": groups,
        "batch_key": batch_key,
        "output_column": output_column,
        "plotted_observations": int(adata.n_obs),
        "plotted_metric_values": int(
            sum(np.isfinite(adata.obs[metric].to_numpy(dtype=float, na_value=np.nan)).sum() for metric in metrics)
        ),
        "missing_metric_values": int(missing_values),
        "plotted_thresholds": int(
            np.isfinite(thresholds[["lower_threshold", "upper_threshold"]].to_numpy(dtype=float)).sum()
        ),
        "threshold_line_legend": "Stored MAD threshold",
        "threshold_evidence": threshold_records,
        "group_rates": rate_records,
        "title": title,
    }
    if _return_diagnostics:
        return png, details
    return png


def run_mad_outlier_plot(adata: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _standalone_mad_outlier_plot(adata, view=view, _return_diagnostics=True)


def mad_outlier_plot_code(*, view: dict[str, Any]) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_mad_outlier_plot)).strip()
    return f"""from __future__ import annotations

{implementation}


def plot_mad_outliers(adata):
    return _standalone_mad_outlier_plot(adata, view={view!r})
"""


def _standalone_scrublet_diagnostics_plot(adata, *, view=None, _return_diagnostics=False):
    import hashlib
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Scrublet Diagnostics Plot"
    if view is None:
        mode = {"view": "score_distributions"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    if set(mode) != {"view"} or mode.get("view") not in {
        "score_distributions",
        "sample_predicted_fraction",
    }:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if not isinstance(adata, AnnData):
        raise TypeError(f"{operation} requires an in-memory AnnData object.")
    if adata.isbacked:
        raise ValueError(f"{operation} does not support backed AnnData; load it into memory first.")
    if adata.n_obs < 1:
        raise ValueError(f"{operation} requires at least one observation.")
    # ponytail: fixed static-PNG ceiling; add sampled/paginated rendering if million-cell plots become necessary.
    if adata.n_obs > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 observations in one static PNG.")

    evidence = adata.uns.get("openbio_scrublet_diagnostics")
    expected_evidence = {
        "schema_version",
        "producer_operation",
        "observation_axis_fingerprint_sha256",
        "batch_key",
        "score_column",
        "prediction_column",
        "group_statistics",
        "score_evidence_fingerprint_sha256",
    }
    if not isinstance(evidence, Mapping) or set(evidence) != expected_evidence:
        raise ValueError(f"{operation} requires the exact stored Scrublet diagnostic evidence schema.")
    if evidence["schema_version"] != 1 or evidence["producer_operation"] != "scrublet":
        raise ValueError(f"{operation} evidence was not produced by Run Scrublet.")

    def axis_fingerprint(names):
        digest = hashlib.sha256()
        for value in names:
            encoded = str(value).encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "little"))
            digest.update(encoded)
        return f"sha256:{digest.hexdigest()}"

    if evidence["observation_axis_fingerprint_sha256"] != axis_fingerprint(adata.obs_names):
        raise ValueError(f"{operation} observation axis fingerprint does not match the stored evidence.")
    batch_key = evidence["batch_key"]
    score_column = evidence["score_column"]
    prediction_column = evidence["prediction_column"]
    if not isinstance(batch_key, str) or score_column != "doublet_score" or prediction_column != "predicted_doublet":
        raise ValueError(f"{operation} stored column identifiers are invalid.")
    metadata = adata.uns.get("openbio_singlecell")
    history = metadata.get("analysis_history") if isinstance(metadata, Mapping) else None
    if not isinstance(history, Mapping) or not any(
        isinstance(entry, Mapping)
        and entry.get("operation") == "scrublet"
        and isinstance(entry.get("parameters"), Mapping)
        and entry["parameters"].get("batch_key") == batch_key
        for entry in history.values()
    ):
        raise ValueError(f"{operation} evidence does not match a Run Scrublet provenance entry.")
    if score_column not in adata.obs or prediction_column not in adata.obs:
        raise ValueError(f"{operation} requires the stored Scrublet observation columns.")
    scores_series = adata.obs[score_column]
    predictions_series = adata.obs[prediction_column]
    if pd.api.types.is_bool_dtype(scores_series.dtype) or not pd.api.types.is_numeric_dtype(scores_series.dtype):
        raise TypeError(f"{operation} observed scores must be real numeric.")
    scores = scores_series.to_numpy(dtype=float, na_value=np.nan)
    if not bool(np.isfinite(scores).all()):
        raise ValueError(f"{operation} observed scores must be finite.")
    if not pd.api.types.is_bool_dtype(predictions_series.dtype) or bool(predictions_series.isna().any()):
        raise TypeError(f"{operation} predictions must be complete boolean evidence.")
    predictions = predictions_series.to_numpy(dtype=bool)

    statistics = evidence["group_statistics"]
    expected_columns = [
        "group",
        "cells",
        "predicted_doublets",
        "predicted_fraction",
        "threshold",
        "simulated_count",
    ]
    if not isinstance(statistics, pd.DataFrame) or statistics.columns.tolist() != expected_columns:
        raise ValueError(f"{operation} requires the exact stored Scrublet group-statistics schema.")
    if statistics.empty:
        raise ValueError(f"{operation} requires at least one stored Sample/group statistic.")
    groups = statistics["group"].tolist()
    if any(not isinstance(value, str) or not value or value != value.strip() for value in groups):
        raise ValueError(f"{operation} stored Sample/group axis must contain canonical strings.")
    if len(set(groups)) != len(groups):
        raise ValueError(f"{operation} stored Sample/group axis contains duplicates.")
    if len(groups) > 100:
        raise ValueError(f"{operation} supports at most 100 Samples/groups in one static PNG; found {len(groups)}.")

    scrublet = adata.uns.get("scrublet")
    if not isinstance(scrublet, Mapping):
        raise ValueError(f"{operation} requires Scanpy's stored uns['scrublet'] result.")
    batches = scrublet.get("batches")
    if batches is None:
        payloads = {"all_cells": scrublet}
    elif isinstance(batches, Mapping):
        payloads = {str(label): payload for label, payload in batches.items() if isinstance(payload, Mapping)}
    else:
        raise ValueError(f"{operation} stored Scrublet batches mapping is invalid.")
    if set(payloads) != set(groups):
        raise ValueError(f"{operation} simulated-score groups do not match the diagnostic group axis.")

    if batch_key:
        if batch_key not in adata.obs or bool(adata.obs[batch_key].isna().any()):
            raise ValueError(f"{operation} requires the complete stored Sample column {batch_key!r}.")
        rendered = [str(value) for value in pd.unique(adata.obs[batch_key])]
        if len(set(rendered)) != len(rendered):
            raise ValueError(f"{operation} Sample labels collide after string rendering.")
        group_masks = {
            group: np.asarray([str(value) == group for value in adata.obs[batch_key]], dtype=bool)
            for group in groups
        }
        if not bool(np.logical_or.reduce(list(group_masks.values())).all()) or any(
            int(mask.sum()) == 0 for mask in group_masks.values()
        ):
            raise ValueError(f"{operation} stored Sample/group axis does not partition the observations.")
    else:
        if groups != ["all_cells"]:
            raise ValueError(f"{operation} global evidence must use the single 'all_cells' group.")
        group_masks = {"all_cells": np.ones(adata.n_obs, dtype=bool)}

    simulated_by_group = {}
    group_records = []
    total_simulated = 0
    for row in statistics.itertuples(index=False):
        payload = payloads[row.group]
        simulated = np.asarray(payload.get("doublet_scores_sim"))
        if (
            simulated.ndim != 1
            or simulated.size < 1
            or not bool(np.issubdtype(simulated.dtype, np.number))
            or bool(np.iscomplexobj(simulated))
            or not bool(np.isfinite(simulated).all())
        ):
            raise ValueError(f"{operation} requires finite one-dimensional simulated scores for {row.group!r}.")
        simulated = simulated.astype(float, copy=False)
        try:
            payload_threshold = float(payload.get("threshold"))
        except (TypeError, ValueError) as error:
            raise ValueError(f"{operation} requires a finite threshold for {row.group!r}.") from error
        if not np.isfinite(payload_threshold):
            raise ValueError(f"{operation} requires a finite threshold for {row.group!r}.")
        mask = group_masks[row.group]
        cells = int(mask.sum())
        predicted = int(predictions[mask].sum())
        fraction = predicted / cells
        integer_values = (row.cells, row.predicted_doublets, row.simulated_count)
        if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0 for value in integer_values):
            raise ValueError(f"{operation} stored group counts must be nonnegative integers.")
        if (
            (int(row.cells), int(row.predicted_doublets), int(row.simulated_count))
            != (cells, predicted, int(simulated.size))
            or not isinstance(row.predicted_fraction, (int, float, np.integer, np.floating))
            or not np.isfinite(row.predicted_fraction)
            or not np.isclose(row.predicted_fraction, fraction, rtol=0.0, atol=1e-12)
            or not isinstance(row.threshold, (int, float, np.integer, np.floating))
            or not np.isfinite(row.threshold)
            or not np.isclose(row.threshold, payload_threshold, rtol=0.0, atol=1e-12)
            or not np.array_equal(predictions[mask], scores[mask] > payload_threshold)
        ):
            raise ValueError(f"{operation} stored group statistics do not match the Scrublet evidence.")
        simulated_by_group[row.group] = simulated
        total_simulated += int(simulated.size)
        group_records.append(
            {
                "group": row.group,
                "cells": cells,
                "predicted_doublets": predicted,
                "predicted_fraction": float(fraction),
                "threshold": payload_threshold,
                "simulated_count": int(simulated.size),
            }
        )
    digest = hashlib.sha256(b"scrublet-score-evidence-v1")
    for values, dtype in ((scores, np.float64), (predictions, np.uint8)):
        array = np.asarray(values, dtype=dtype)
        digest.update(str(tuple(array.shape)).encode("ascii"))
        digest.update(np.ascontiguousarray(array).tobytes())
    for record in group_records:
        group = record["group"]
        encoded = group.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
        digest.update(np.asarray([record["threshold"]], dtype=np.float64).tobytes())
        simulated = np.asarray(simulated_by_group[group], dtype=np.float64)
        digest.update(str(tuple(simulated.shape)).encode("ascii"))
        digest.update(np.ascontiguousarray(simulated).tobytes())
    if evidence["score_evidence_fingerprint_sha256"] != f"sha256:{digest.hexdigest()}":
        raise ValueError(f"{operation} score evidence fingerprint does not match the stored values.")
    if total_simulated > 2_000_000:
        raise ValueError(
            f"{operation} supports at most 2,000,000 simulated scores in one static PNG; found {total_simulated:,}."
        )
    if mode["view"] == "sample_predicted_fraction" and not batch_key:
        raise ValueError(f"{operation} Sample predicted fraction requires a non-empty producer batch_key.")
    if mode["view"] == "score_distributions":
        if len(groups) > 12:
            raise ValueError(
                f"{operation} score distributions support at most 12 Samples/groups in one static PNG; "
                f"found {len(groups)}."
            )
        columns = min(3, len(groups))
        rows = int(np.ceil(len(groups) / columns))
        figure = Figure(figsize=(5.0 * columns, 3.6 * rows), constrained_layout=True)
        FigureCanvasAgg(figure)
        axes = np.atleast_1d(figure.subplots(rows, columns, squeeze=False)).reshape(-1)
        for axis, record in zip(axes, group_records, strict=False):
            group = record["group"]
            observed = scores[group_masks[group]]
            simulated = simulated_by_group[group]
            bins = min(40, max(8, int(np.ceil(np.sqrt(max(observed.size, simulated.size))))))
            axis.hist(observed, bins=bins, density=True, alpha=0.55, color="#4C78A8", label="Observed")
            axis.hist(simulated, bins=bins, density=True, histtype="step", linewidth=1.6, color="#F58518", label="Simulated")
            axis.axvline(record["threshold"], color="#E45756", linestyle="--", linewidth=1.4, label="Threshold")
            axis.set_title(group)
            axis.set_xlabel("Doublet score")
            axis.set_ylabel("Density")
            axis.grid(axis="y", alpha=0.2)
            axis.legend()
        for axis in axes[len(groups) :]:
            axis.set_visible(False)
        title = "Scrublet observed and simulated score distributions"
    else:
        figure = Figure(figsize=(max(6.4, min(20.0, 0.55 * len(groups) + 2.5)), 5.0))
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        positions = np.arange(len(groups), dtype=float)
        fractions = np.asarray([record["predicted_fraction"] for record in group_records], dtype=float)
        axis.bar(positions, fractions, color="#4C78A8")
        axis.set_xticks(positions, groups, rotation=45, ha="right")
        axis.set_ylim(0.0, 1.0)
        axis.set_ylabel("Predicted-doublet fraction")
        axis.set_title("Scrublet predicted-doublet fraction by Sample")
        axis.grid(axis="y", alpha=0.2)
        title = "Scrublet predicted-doublet fraction by Sample"
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")

    details = {
        "view": mode["view"],
        "view_parameters": mode,
        "groups": groups,
        "batch_key": batch_key,
        "score_column": score_column,
        "prediction_column": prediction_column,
        "plotted_observed_scores": int(adata.n_obs),
        "plotted_simulated_scores": total_simulated,
        "plotted_thresholds": len(groups),
        "group_statistics": group_records,
        "title": title,
    }
    if _return_diagnostics:
        return png, details
    return png


def run_scrublet_diagnostics_plot(adata: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _standalone_scrublet_diagnostics_plot(adata, view=view, _return_diagnostics=True)


def scrublet_diagnostics_plot_code(*, view: dict[str, Any]) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_scrublet_diagnostics_plot)).strip()
    return f"""from __future__ import annotations

{implementation}


def plot_scrublet_diagnostics(adata):
    return _standalone_scrublet_diagnostics_plot(adata, view={view!r})
"""


__all__ = [
    "mad_outlier_plot_code",
    "run_mad_outlier_plot",
    "run_scrublet_diagnostics_plot",
    "scrublet_diagnostics_plot_code",
]
