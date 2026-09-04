from __future__ import annotations

import inspect
import textwrap

from .cell_cycle import _cell_cycle_axis_fingerprint, _cell_cycle_score_fingerprint


def _standalone_cell_cycle_score_plot(
    adata,
    *,
    output_prefix="cell_cycle",
    _return_details=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Cell Cycle Score Plot"
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires in-memory AnnData.")
    if int(getattr(adata, "n_obs", 0)) < 1 or not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires a non-empty unique observation axis.")
    if int(adata.n_obs) > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 cells in one static PNG.")
    if not isinstance(output_prefix, str) or not output_prefix or output_prefix != output_prefix.strip():
        raise ValueError(f"{operation} output_prefix must be a canonical nonempty string.")
    if any(character.isspace() for character in output_prefix):
        raise ValueError(f"{operation} output_prefix cannot contain whitespace.")
    score_columns = {
        "s_score": f"{output_prefix}_s_score",
        "g2m_score": f"{output_prefix}_g2m_score",
        "phase": f"{output_prefix}_phase",
    }
    missing = [column for column in score_columns.values() if column not in adata.obs]
    if missing:
        raise ValueError(f"{operation} requires the complete stored score bundle; missing {missing!r}.")

    openbio = adata.uns.get("openbio_singlecell")
    history = openbio.get("analysis_history") if isinstance(openbio, Mapping) else None
    matching = []
    if isinstance(history, Mapping):
        matching = [
            entry
            for entry in history.values()
            if isinstance(entry, Mapping)
            and entry.get("operation") == "cell_cycle_score"
            and isinstance(entry.get("parameters"), Mapping)
            and entry["parameters"].get("output_prefix") == output_prefix
        ]
    if not matching:
        raise ValueError(f"{operation} could not verify Cell Cycle Score producer history.")
    parameters = matching[-1]["parameters"]
    required_parameters = {
        "source_kind",
        "gene_set_source",
        "organism",
        "output_prefix",
        "observation_axis_fingerprint_sha256",
        "score_bundle_fingerprint_sha256",
    }
    if not required_parameters.issubset(parameters):
        raise ValueError(f"{operation} producer history lacks the complete Diagnostic evidence contract.")
    if parameters["source_kind"] not in {"X", "raw", "layer"}:
        raise ValueError(f"{operation} producer expression-source provenance is invalid.")

    observation_ids = list(adata.obs_names)
    axis_fingerprint = _cell_cycle_axis_fingerprint(observation_ids)
    if parameters["observation_axis_fingerprint_sha256"] != axis_fingerprint:
        raise ValueError(f"{operation} observation axis fingerprint does not match producer evidence.")
    s_score = pd.to_numeric(adata.obs[score_columns["s_score"]], errors="raise").to_numpy(dtype=float)
    g2m_score = pd.to_numeric(adata.obs[score_columns["g2m_score"]], errors="raise").to_numpy(dtype=float)
    if s_score.shape != (len(observation_ids),) or g2m_score.shape != (len(observation_ids),):
        raise ValueError(f"{operation} score columns are not observation aligned.")
    if not bool(np.isfinite(s_score).all()) or not bool(np.isfinite(g2m_score).all()):
        raise ValueError(f"{operation} scores must be finite.")
    phase_series = adata.obs[score_columns["phase"]]
    if not isinstance(phase_series.dtype, pd.CategoricalDtype) or phase_series.cat.categories.tolist() != [
        "G1",
        "S",
        "G2M",
    ]:
        raise ValueError(f"{operation} phase must retain the canonical categorical axis G1, S, G2M.")
    if bool(phase_series.isna().any()):
        raise ValueError(f"{operation} phase labels cannot be missing.")
    phases = phase_series.astype(str).tolist()
    score_fingerprint = _cell_cycle_score_fingerprint(
        observation_ids,
        s_score,
        g2m_score,
        phases,
    )
    if parameters["score_bundle_fingerprint_sha256"] != score_fingerprint:
        raise ValueError(f"{operation} score bundle fingerprint does not match producer evidence.")
    expected_phases = np.where(
        (s_score < 0.0) & (g2m_score < 0.0),
        "G1",
        np.where(g2m_score > s_score, "G2M", "S"),
    ).tolist()
    if phases != expected_phases:
        raise ValueError(f"{operation} phase labels disagree with the stored Scanpy score rule.")

    phase_order = ["G1", "S", "G2M"]
    phase_counts = {phase: phases.count(phase) for phase in phase_order}
    phase_proportions = {phase: phase_counts[phase] / len(phases) for phase in phase_order}
    colors = {"G1": "#9C9C9C", "S": "#4C78A8", "G2M": "#E45756"}
    figure = Figure(figsize=(11.0, 4.8), constrained_layout=True)
    FigureCanvasAgg(figure)
    score_axis, count_axis = figure.subplots(1, 2)
    for phase in phase_order:
        mask = np.asarray(phases, dtype=object) == phase
        score_axis.scatter(
            s_score[mask],
            g2m_score[mask],
            s=24.0,
            alpha=0.8,
            linewidths=0.0,
            color=colors[phase],
            label=phase,
        )
    lower = float(min(s_score.min(), g2m_score.min()))
    upper = float(max(s_score.max(), g2m_score.max()))
    score_axis.plot([lower, upper], [lower, upper], color="#333333", linewidth=0.8, linestyle="--")
    score_axis.axhline(0.0, color="#777777", linewidth=0.6)
    score_axis.axvline(0.0, color="#777777", linewidth=0.6)
    score_axis.set_xlabel("S score")
    score_axis.set_ylabel("G2/M score")
    score_axis.set_title("Stored cell-cycle score plane")
    score_axis.legend(title="Phase")
    score_axis.grid(alpha=0.15)
    count_axis.bar(phase_order, [phase_counts[phase] for phase in phase_order], color=[colors[p] for p in phase_order])
    count_axis.set_xlabel("Assigned phase")
    count_axis.set_ylabel("Cells")
    count_axis.set_title("Stored phase calls")
    count_axis.grid(axis="y", alpha=0.2)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")

    def numeric_summary(values):
        quantiles = np.quantile(values, [0.0, 0.25, 0.5, 0.75, 1.0])
        return {
            "n": int(len(values)),
            "missing": 0,
            "min": float(quantiles[0]),
            "q1": float(quantiles[1]),
            "median": float(quantiles[2]),
            "mean": float(np.mean(values)),
            "q3": float(quantiles[3]),
            "max": float(quantiles[4]),
        }

    details = {
        "output_prefix": output_prefix,
        "score_columns": score_columns,
        "phase_order": phase_order,
        "phase_counts": phase_counts,
        "phase_proportions": phase_proportions,
        "s_score_summary": numeric_summary(s_score),
        "g2m_score_summary": numeric_summary(g2m_score),
        "plotted_cells": len(observation_ids),
        "expression_source": {
            "source": parameters["source_kind"],
            "layer_name": parameters.get("layer_name") if parameters["source_kind"] == "layer" else None,
        },
        "gene_set_source": parameters["gene_set_source"],
        "organism": parameters["organism"],
        "observation_axis_fingerprint_sha256": axis_fingerprint,
        "score_bundle_fingerprint_sha256": score_fingerprint,
        "phase_rule": "G1 when both scores are negative; otherwise G2M when G2M > S, else S (ties resolve to S)",
        "rendering": {"figure_size_inches": [11.0, 4.8], "dpi": 120, "maximum_cells": 1_000_000},
    }
    return (png, details) if _return_details else png


def cell_cycle_score_plot_code(*, output_prefix: str) -> str:
    helpers = "\n\n".join(
        textwrap.dedent(inspect.getsource(helper)).strip()
        for helper in (_cell_cycle_axis_fingerprint, _cell_cycle_score_fingerprint)
    )
    implementation = textwrap.dedent(inspect.getsource(_standalone_cell_cycle_score_plot)).strip()
    return f"""from __future__ import annotations

import hashlib
import json

{helpers}


{implementation}


def plot_cell_cycle_scores(adata):
    return _standalone_cell_cycle_score_plot(adata, output_prefix={output_prefix!r})
"""


__all__ = ["_standalone_cell_cycle_score_plot", "cell_cycle_score_plot_code"]
