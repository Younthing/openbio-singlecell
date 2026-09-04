from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .cnv_analysis import (
    _CNV_PROVENANCE_KEY,
    CNV_STATE_SCHEMA,
    CNV_STATE_TYPE,
    INFER_CNV_NODE_ID,
    _cnv_artifact_fingerprints,
    _cnv_artifact_hash,
    _cnv_axes,
    _cnv_axis_hash,
    _cnv_json,
    _cnv_matrix_hash,
    _cnv_score_table_fingerprint,
    _cnv_sha_json,
    _cnv_validate_state,
)


def _standalone_cnv_heatmap_plot(
    cnv_state,
    *,
    groupby="",
    view=None,
    _return_details=False,
):
    import io

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from scipy import sparse

    operation = "CNV Heatmap Plot"
    adata, metadata = _cnv_validate_state(
        cnv_state,
        numpy=np,
        scipy_sparse=sparse,
    )
    provenance = adata.uns.get(_CNV_PROVENANCE_KEY)
    provenance_fields = {
        "schema",
        "output_key",
        "genome_assembly",
        "reference_key",
        "reference_categories",
        "sample_key",
        "windows",
        "infercnvpy_version",
    }
    if not isinstance(provenance, dict) or set(provenance) != provenance_fields:
        raise ValueError(f"{operation} requires exact OpenBio inferred-CNV provenance.")
    for key in ("schema", "output_key", "genome_assembly", "reference_key", "sample_key"):
        if provenance[key] != metadata[key if key != "schema" else "schema"]:
            raise ValueError(f"{operation} provenance differs from the typed state for {key!r}.")
    windows = provenance["windows"]
    if not isinstance(windows, list) or not windows:
        raise ValueError(f"{operation} requires nonempty genome-window metadata.")
    chromosomes = []
    chromosome_starts = []
    expected_start = 0
    for record in windows:
        if not isinstance(record, dict) or set(record) != {
            "chromosome",
            "genes",
            "first_window_index",
            "window_count",
        }:
            raise ValueError(f"{operation} genome-window records have an invalid exact schema.")
        chromosome = record["chromosome"]
        genes = record["genes"]
        first = record["first_window_index"]
        count = record["window_count"]
        if not isinstance(chromosome, str) or not chromosome.startswith("chr") or chromosome != chromosome.strip():
            raise ValueError(f"{operation} chromosome labels must be canonical chr-prefixed strings.")
        if chromosome in chromosomes:
            raise ValueError(f"{operation} chromosome window blocks must be unique and contiguous.")
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 1 for value in (genes, count)):
            raise ValueError(f"{operation} genome-window gene/count fields must be positive integers.")
        if isinstance(first, bool) or not isinstance(first, int) or first != expected_start:
            raise ValueError(f"{operation} genome-window blocks must form a zero-based contiguous axis.")
        chromosomes.append(chromosome)
        chromosome_starts.append(first)
        expected_start += count
    matrix = adata.obsm[f"X_{metadata['output_key']}"]
    if tuple(matrix.shape) != (int(adata.n_obs), expected_start):
        raise ValueError(f"{operation} CNV matrix does not align to the observation and genome-window axes.")
    stored_values = np.asarray(matrix.data if sparse.issparse(matrix) else matrix).ravel()
    if stored_values.size and (
        not np.issubdtype(stored_values.dtype, np.number)
        or np.iscomplexobj(stored_values)
        or not bool(np.isfinite(stored_values).all())
    ):
        raise ValueError(f"{operation} CNV matrix must contain finite real values.")
    chr_pos = adata.uns.get(metadata["output_key"])
    expected_chr_pos = dict(zip(chromosomes, chromosome_starts, strict=True))
    if not isinstance(chr_pos, dict) or set(chr_pos) != {"chr_pos"} or chr_pos["chr_pos"] != expected_chr_pos:
        raise ValueError(f"{operation} chromosome start coordinates do not match the genome-window axis.")
    if not isinstance(groupby, str):
        raise TypeError(f"{operation} groupby must be a string.")
    groupby = groupby.strip() or metadata["reference_key"]
    if groupby not in adata.obs:
        raise ValueError(f"{operation} groupby column not found in obs[{groupby!r}].")
    groups = adata.obs[groupby]
    if not isinstance(groups.dtype, pd.CategoricalDtype):
        raise TypeError(f"{operation} obs[{groupby!r}] must be explicitly categorical.")
    if bool(groups.isna().any()):
        raise ValueError(f"{operation} obs[{groupby!r}] contains missing labels.")
    observed_values = list(groups.astype(object))
    if groups.cat.ordered:
        group_values = [value for value in groups.cat.categories if value in set(observed_values)]
    else:
        group_values = list(dict.fromkeys(observed_values))
    group_labels = [str(value) for value in group_values]
    if any(not label or label != label.strip() for label in group_labels) or len(set(group_labels)) != len(group_labels):
        raise ValueError(f"{operation} group labels are blank or collide after display rendering.")
    group_records = [
        {"group": label, "cells": int(sum(value == group for value in observed_values))}
        for group, label in zip(group_values, group_labels, strict=True)
    ]
    mode = {"view": "group_mean", "max_groups": 50} if view is None else view
    variants = {"group_mean": {"view", "max_groups"}, "cells": {"view", "max_cells"}}
    if not isinstance(mode, dict) or mode.get("view") not in variants or set(mode) != variants.get(mode.get("view")):
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    limit_key = "max_groups" if mode["view"] == "group_mean" else "max_cells"
    limit = mode[limit_key]
    maximum = 256 if limit_key == "max_groups" else 10_000
    if isinstance(limit, bool) or not isinstance(limit, int) or not 2 <= limit <= maximum:
        raise ValueError(f"{operation} {limit_key} must be an integer between 2 and {maximum:,}.")
    if mode["view"] == "group_mean":
        if len(group_values) > limit:
            raise ValueError(f"{operation} has {len(group_values):,} groups, exceeding max_groups={limit:,}.")
        planned_rows = len(group_values)
    else:
        if int(adata.n_obs) > limit:
            raise ValueError(f"{operation} has {adata.n_obs:,} cells, exceeding max_cells={limit:,}; use group_mean.")
        planned_rows = int(adata.n_obs)
    planned_values = planned_rows * expected_start
    if planned_values > 5_000_000:
        raise ValueError(f"{operation} supports at most 5,000,000 displayed row-by-window values.")

    if mode["view"] == "group_mean":
        rows = []
        for value in group_values:
            mask = np.asarray([observed == value for observed in observed_values], dtype=bool)
            selected = matrix[mask, :]
            mean = np.asarray(selected.mean(axis=0), dtype=float).ravel()
            rows.append(mean)
        plotted = np.vstack(rows)
        row_labels = group_labels
        row_centers = np.arange(len(group_labels), dtype=float)
    else:
        positions = np.concatenate(
            [np.flatnonzero(np.asarray([observed == value for observed in observed_values], dtype=bool)) for value in group_values]
        )
        selected = matrix[positions, :]
        plotted = np.asarray(selected.toarray() if sparse.issparse(selected) else selected, dtype=float)
        row_labels = group_labels
        offsets = np.cumsum([0, *[record["cells"] for record in group_records]])
        row_centers = np.asarray([(offsets[index] + offsets[index + 1] - 1) / 2 for index in range(len(group_records))])
    if int(plotted.size) != planned_values:
        raise RuntimeError(f"{operation} materialized a matrix inconsistent with its display-size preflight.")
    if not bool(np.isfinite(plotted).all()):
        raise ValueError(f"{operation} displayed matrix contains non-finite values.")
    color_limit = float(np.max(np.abs(plotted))) if plotted.size else 0.0
    if color_limit == 0:
        color_limit = 1.0
    figure = Figure(
        figsize=(max(8.0, min(20.0, expected_start * 0.08 + 5.0)), max(4.5, min(20.0, plotted.shape[0] * 0.12 + 2.5))),
        constrained_layout=True,
    )
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    image = axis.imshow(
        plotted,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        vmin=-color_limit,
        vmax=color_limit,
    )
    chromosome_ends = [
        windows[index + 1]["first_window_index"] if index + 1 < len(windows) else expected_start
        for index in range(len(windows))
    ]
    chromosome_centers = [(start + end - 1) / 2 for start, end in zip(chromosome_starts, chromosome_ends, strict=True)]
    axis.set_xticks(chromosome_centers, chromosomes)
    axis.set_yticks(row_centers, row_labels)
    for boundary in chromosome_starts[1:]:
        axis.axvline(boundary - 0.5, color="black", linewidth=0.8)
    axis.set_xlabel(f"Genome-ordered inferred-CNV windows ({metadata['genome_assembly']})")
    axis.set_ylabel(groupby if mode["view"] == "group_mean" else f"Cells ordered by {groupby}")
    axis.set_title("Expression-derived inferred-CNV window heatmap")
    figure.colorbar(image, ax=axis, label="Centered inferred-CNV signal")
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "view": mode["view"],
        "view_parameters": dict(mode),
        "groupby": groupby,
        "groups": group_records,
        "row_count": int(plotted.shape[0]),
        "window_count": expected_start,
        "plotted_values": int(plotted.size),
        "chromosomes": chromosomes,
        "chromosome_start_indices": chromosome_starts,
        "genome_assembly": metadata["genome_assembly"],
        "reference_key": metadata["reference_key"],
        "sample_key": metadata["sample_key"],
        "color_limit": color_limit,
        "displayed_min": float(plotted.min()),
        "displayed_max": float(plotted.max()),
        "cnv_matrix_fingerprint_sha256": metadata["fingerprints"]["cnv_matrix_sha256"],
        "artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
        "title": "Expression-derived CNV window heatmap",
    }
    if _return_details:
        return png, details
    return png


def run_cnv_heatmap_plot(
    cnv_state: Any,
    *,
    groupby: str = "",
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_cnv_heatmap_plot(
        cnv_state,
        groupby=groupby,
        view=view,
        _return_details=True,
    )


def _cnv_plot_source(*functions: object) -> str:
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence

CNV_STATE_TYPE = {CNV_STATE_TYPE!r}
CNV_STATE_SCHEMA = {CNV_STATE_SCHEMA!r}
INFER_CNV_NODE_ID = {INFER_CNV_NODE_ID!r}
_CNV_PROVENANCE_KEY = {_CNV_PROVENANCE_KEY!r}

{implementation}
"""


def cnv_heatmap_plot_code(*, groupby: str, view: dict[str, Any]) -> str:
    source = _cnv_plot_source(
        _cnv_json,
        _cnv_sha_json,
        _cnv_axis_hash,
        _cnv_matrix_hash,
        _cnv_axes,
        _cnv_artifact_fingerprints,
        _cnv_artifact_hash,
        _cnv_validate_state,
        _standalone_cnv_heatmap_plot,
    )
    wrapper = f"""

def plot_cnv_heatmap(cnv_state):
    return _standalone_cnv_heatmap_plot(cnv_state, groupby={groupby!r}, view={view!r})
"""
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, "<plot_cnv_heatmap>", "exec")
    return code


def _standalone_cnv_score_plot(table, *, producer, view=None, _return_details=False):
    import io

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "CNV Score Plot"
    if not isinstance(producer, dict) or set(producer) != {
        "operation",
        "parameters",
        "input_cells",
        "input_genes",
        "random_seed",
    }:
        raise ValueError(f"{operation} requires exact producer provenance.")
    if producer["operation"] != "cnv_score" or not isinstance(producer["parameters"], dict):
        raise ValueError(f"{operation} requires the CNV Group Score producer.")
    parameters = producer["parameters"]
    required_parameters = {
        "groupby",
        "use_rep",
        "output_key",
        "overwrite_existing",
        "inplace",
        "independent_validation_rtol",
        "independent_validation_atol",
        "input_cnv_state_fingerprint_sha256",
        "input_cnv_matrix_fingerprint_sha256",
        "partition_fingerprint_sha256",
        "output_fingerprint_sha256",
        "table_fingerprint_sha256",
    }
    if set(parameters) != required_parameters:
        raise ValueError(f"{operation} producer parameters have an invalid exact schema.")
    for key in (
        "input_cnv_state_fingerprint_sha256",
        "input_cnv_matrix_fingerprint_sha256",
        "partition_fingerprint_sha256",
        "output_fingerprint_sha256",
        "table_fingerprint_sha256",
    ):
        value = parameters[key]
        if not isinstance(value, str) or len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
            raise ValueError(f"{operation} producer {key!r} is not a SHA-256 fingerprint.")
    columns = [
        "rank",
        "group_type",
        "group_value",
        "group_display",
        "cell_count",
        "cnv_score",
        "absolute_min",
        "absolute_q1",
        "absolute_median",
        "absolute_mean",
        "absolute_q3",
        "absolute_max",
    ]
    if not isinstance(table, pd.DataFrame) or list(table.columns) != columns or table.empty:
        raise ValueError(f"{operation} requires the nonempty canonical exact CNV score columns.")
    observed_fingerprint = _cnv_score_table_fingerprint(table)
    if observed_fingerprint != parameters["table_fingerprint_sha256"]:
        raise ValueError(f"{operation} table fingerprint does not match the producer evidence.")
    if len(table) > 256:
        raise ValueError(f"{operation} supports at most 256 scored groups.")
    ranks = table["rank"].to_numpy()
    if not np.issubdtype(ranks.dtype, np.integer) or not np.array_equal(ranks, np.arange(1, len(table) + 1)):
        raise ValueError(f"{operation} rank axis must be consecutive one-based integers.")
    counts = table["cell_count"].to_numpy()
    if not np.issubdtype(counts.dtype, np.integer) or bool((counts < 1).any()):
        raise ValueError(f"{operation} cell_count must contain positive integers.")
    if int(counts.sum()) != producer["input_cells"]:
        raise ValueError(f"{operation} cell counts do not cover the producer observation axis.")
    displays = table["group_display"].tolist()
    if any(not isinstance(value, str) or not value or value != value.strip() for value in displays) or len(set(displays)) != len(displays):
        raise ValueError(f"{operation} group display axis must contain unique canonical strings.")
    for row in table.itertuples(index=False):
        value = row.group_value
        if row.group_type == "string":
            valid = isinstance(value, str) and value == row.group_display
        elif row.group_type == "integer":
            valid = isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)) and str(int(value)) == row.group_display
        elif row.group_type == "number":
            valid = isinstance(value, (float, np.floating)) and np.isfinite(value) and repr(float(value)) == row.group_display
        else:
            valid = False
        if not valid:
            raise ValueError(f"{operation} group identity records are inconsistent.")
    numeric_columns = [
        "cnv_score",
        "absolute_min",
        "absolute_q1",
        "absolute_median",
        "absolute_mean",
        "absolute_q3",
        "absolute_max",
    ]
    values = table[numeric_columns].to_numpy(dtype=float)
    if not bool(np.isfinite(values).all()) or bool((values < 0).any()):
        raise ValueError(f"{operation} score and absolute-value evidence must be finite and nonnegative.")
    for row in table.itertuples(index=False):
        ordered = [row.absolute_min, row.absolute_q1, row.absolute_median, row.absolute_q3, row.absolute_max]
        if any(left > right for left, right in zip(ordered, ordered[1:], strict=False)):
            raise ValueError(f"{operation} absolute-value quantiles are not ordered.")
        if not row.absolute_min <= row.absolute_mean <= row.absolute_max:
            raise ValueError(f"{operation} absolute mean falls outside the stored range.")
        if not np.isclose(
            row.cnv_score,
            row.absolute_mean,
            rtol=float(parameters["independent_validation_rtol"]),
            atol=float(parameters["independent_validation_atol"]),
        ):
            raise ValueError(f"{operation} CNV score disagrees with the stored group absolute mean.")
    scores = table["cnv_score"].to_numpy(dtype=float)
    expected_order = sorted(range(len(table)), key=lambda index: (-scores[index], displays[index]))
    if expected_order != list(range(len(table))):
        raise ValueError(f"{operation} group rows do not preserve the canonical score ranking.")
    mode = {"view": "score_bar", "max_groups": 30} if view is None else view
    if not isinstance(mode, dict) or set(mode) != {"view", "max_groups"} or mode.get("view") not in {
        "score_bar",
        "absolute_interval",
    }:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    maximum = mode["max_groups"]
    if isinstance(maximum, bool) or not isinstance(maximum, int) or not 1 <= maximum <= 100:
        raise ValueError(f"{operation} max_groups must be an integer between 1 and 100.")
    selected = table.iloc[:maximum]
    groups = selected["group_display"].tolist()
    positions = np.arange(len(selected), dtype=float)
    figure = Figure(figsize=(8.0, max(4.0, min(18.0, 0.4 * len(selected) + 2.0))), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    if mode["view"] == "score_bar":
        axis.barh(positions, selected["cnv_score"].to_numpy(dtype=float), color="#4C78A8")
        axis.set_xlabel("Group mean absolute inferred-CNV signal")
        title = "Descriptive inferred-CNV group scores"
    else:
        medians = selected["absolute_median"].to_numpy(dtype=float)
        q1 = selected["absolute_q1"].to_numpy(dtype=float)
        q3 = selected["absolute_q3"].to_numpy(dtype=float)
        minimum = selected["absolute_min"].to_numpy(dtype=float)
        maximum_values = selected["absolute_max"].to_numpy(dtype=float)
        axis.errorbar(
            medians,
            positions,
            xerr=np.vstack([medians - q1, q3 - medians]),
            fmt="o",
            color="#4C78A8",
            capsize=3,
            label="median and IQR",
        )
        for position, left, right in zip(positions, minimum, maximum_values, strict=True):
            axis.plot([left, right], [position, position], color="#4C78A8", alpha=0.35, linewidth=1)
        axis.scatter(selected["absolute_mean"], positions, marker="D", color="#F58518", label="mean / score")
        axis.set_xlabel("Absolute inferred-CNV window signal")
        axis.legend()
        title = "Within-group absolute inferred-CNV summaries"
    axis.set_yticks(positions, groups)
    axis.invert_yaxis()
    axis.set_ylabel(parameters["groupby"])
    axis.set_title(title)
    axis.grid(axis="x", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": mode["view"],
        "view_parameters": dict(mode),
        "groupby": parameters["groupby"],
        "groups": groups,
        "ranks": [int(value) for value in selected["rank"]],
        "cell_counts": [int(value) for value in selected["cell_count"]],
        "cnv_scores": [float(value) for value in selected["cnv_score"]],
        "absolute_min": [float(value) for value in selected["absolute_min"]],
        "absolute_q1": [float(value) for value in selected["absolute_q1"]],
        "absolute_median": [float(value) for value in selected["absolute_median"]],
        "absolute_mean": [float(value) for value in selected["absolute_mean"]],
        "absolute_q3": [float(value) for value in selected["absolute_q3"]],
        "absolute_max": [float(value) for value in selected["absolute_max"]],
        "available_groups": len(table),
        "plotted_groups": len(selected),
        "table_fingerprint_sha256": observed_fingerprint,
        "input_cnv_state_fingerprint_sha256": parameters["input_cnv_state_fingerprint_sha256"],
        "input_cnv_matrix_fingerprint_sha256": parameters["input_cnv_matrix_fingerprint_sha256"],
        "partition_fingerprint_sha256": parameters["partition_fingerprint_sha256"],
        "score_output_fingerprint_sha256": parameters["output_fingerprint_sha256"],
        "title": title,
    }
    if _return_details:
        return png, details
    return png


def run_cnv_score_plot(
    table: Any,
    *,
    producer: dict[str, Any],
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_cnv_score_plot(
        table,
        producer=producer,
        view=view,
        _return_details=True,
    )


def cnv_score_plot_code(*, producer: dict[str, Any], view: dict[str, Any]) -> str:
    source = _cnv_plot_source(
        _cnv_json,
        _cnv_sha_json,
        _cnv_score_table_fingerprint,
        _standalone_cnv_score_plot,
    )
    wrapper = f"""

def plot_cnv_score(table):
    return _standalone_cnv_score_plot(table, producer={producer!r}, view={view!r})
"""
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, "<plot_cnv_score>", "exec")
    return code


__all__ = [
    "cnv_heatmap_plot_code",
    "cnv_score_plot_code",
    "run_cnv_heatmap_plot",
    "run_cnv_score_plot",
]
