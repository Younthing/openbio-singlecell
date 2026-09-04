from __future__ import annotations

import inspect
import textwrap
from typing import Any


def _standalone_cnmf_rank_plot(table, *, producer, view=None, _return_diagnostics=False):
    import hashlib
    import io
    import json
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "cNMF Rank Plot"
    if view is None:
        mode = {"view": "stability_prediction_error"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    if set(mode) != {"view"} or mode.get("view") not in {
        "stability_prediction_error",
        "restart_completion",
    }:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if not isinstance(producer, Mapping) or set(producer) != {
        "operation",
        "parameters",
        "input_cells",
        "input_genes",
    }:
        raise ValueError(f"{operation} requires exact cNMF Rank Survey producer provenance.")
    if producer["operation"] != "cnmf_rank_survey":
        raise ValueError(f"{operation} requires a cNMF Rank Survey k_metrics artifact.")
    parameters = producer["parameters"]
    if not isinstance(parameters, Mapping):
        raise ValueError(f"{operation} producer parameters must be a mapping.")
    required_parameters = {
        "candidate_ks",
        "components_min",
        "components_max",
        "n_iter",
        "total_restarts",
        "source_features",
        "input_fingerprint",
        "k_metrics_fingerprint_sha256",
        "fixed_policy",
    }
    if not required_parameters.issubset(parameters):
        raise ValueError(f"{operation} producer provenance is missing cNMF rank-survey fields.")
    for name in ("input_cells", "input_genes"):
        value = producer[name]
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0:
            raise ValueError(f"{operation} {name} provenance must be a nonnegative integer.")
    fingerprint = parameters["input_fingerprint"]
    if not isinstance(fingerprint, str) or not fingerprint.startswith("sha256:") or len(fingerprint) != 71:
        raise ValueError(f"{operation} requires the cNMF input fingerprint provenance.")
    try:
        int(fingerprint[7:], 16)
    except ValueError as error:
        raise ValueError(f"{operation} requires the cNMF input fingerprint provenance.") from error
    metrics_fingerprint = parameters["k_metrics_fingerprint_sha256"]
    if (
        not isinstance(metrics_fingerprint, str)
        or not metrics_fingerprint.startswith("sha256:")
        or len(metrics_fingerprint) != 71
    ):
        raise ValueError(f"{operation} requires the cNMF rank metrics fingerprint provenance.")

    columns = [
        "k",
        "stability",
        "prediction_error",
        "statistics_density_threshold",
        "expected_restarts",
        "completed_restarts",
        "combined_components",
    ]
    if not isinstance(table, pd.DataFrame):
        raise TypeError(f"{operation} requires a pandas DataFrame.")
    if table.columns.tolist() != columns:
        raise ValueError(f"{operation} requires canonical columns {columns!r} in that order.")
    if table.empty:
        raise ValueError(f"{operation} requires at least one stored candidate rank.")
    # ponytail: fixed static-PNG ceiling; add pagination if surveys with more candidate ranks become useful.
    if len(table) > 100:
        raise ValueError(f"{operation} supports at most 100 candidate ranks in one static PNG; found {len(table)}.")
    integer_columns = ["k", "expected_restarts", "completed_restarts", "combined_components"]
    integers = {}
    for column in integer_columns:
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_integer_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} must contain canonical integers.")
        integers[column] = series.to_numpy(dtype=np.int64)
    ks = integers["k"]
    if bool((ks < 2).any()) or not bool(np.all(np.diff(ks) > 0)):
        raise ValueError(f"{operation} K values must be unique, strictly increasing integers no smaller than 2.")
    if not bool(np.array_equal(ks, np.arange(ks[0], ks[-1] + 1, dtype=np.int64))):
        raise ValueError(f"{operation} K values must form the complete candidate-rank interval.")
    candidate_value = parameters["candidate_ks"]
    if hasattr(candidate_value, "tolist"):
        candidate_value = candidate_value.tolist()
    if not isinstance(candidate_value, (list, tuple)) or list(candidate_value) != ks.tolist():
        raise ValueError(f"{operation} K axis does not match cNMF Rank Survey provenance.")
    if parameters["components_min"] != int(ks[0]) or parameters["components_max"] != int(ks[-1]):
        raise ValueError(f"{operation} K interval does not match cNMF Rank Survey provenance.")
    if parameters["source_features"] != producer["input_genes"]:
        raise ValueError(f"{operation} feature-axis provenance is inconsistent.")
    n_iter = parameters["n_iter"]
    total_restarts = parameters["total_restarts"]
    if (
        isinstance(n_iter, (bool, np.bool_))
        or not isinstance(n_iter, (int, np.integer))
        or n_iter < 2
        or isinstance(total_restarts, (bool, np.bool_))
        or not isinstance(total_restarts, (int, np.integer))
        or total_restarts != len(ks) * int(n_iter)
    ):
        raise ValueError(f"{operation} restart provenance is invalid.")
    if not bool(np.equal(integers["expected_restarts"], int(n_iter)).all()):
        raise ValueError(f"{operation} expected restart counts do not match producer provenance.")
    if not bool(np.array_equal(integers["completed_restarts"], integers["expected_restarts"])):
        raise ValueError(f"{operation} requires complete restart evidence for every candidate K.")
    if not bool(np.array_equal(integers["combined_components"], ks * integers["completed_restarts"])):
        raise ValueError(f"{operation} combined-component counts do not match K times completed restarts.")

    numeric = {}
    for column in ("stability", "prediction_error", "statistics_density_threshold"):
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} {column!r} must contain real numeric values.")
        values = series.to_numpy(dtype=float, na_value=np.nan)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} {column!r} values must be finite.")
        numeric[column] = values
    if bool(((numeric["stability"] < -1.0) | (numeric["stability"] > 1.0)).any()):
        raise ValueError(f"{operation} stability values must lie between -1 and 1.")
    if bool((numeric["prediction_error"] < 0.0).any()):
        raise ValueError(f"{operation} prediction error must be nonnegative.")
    if bool((numeric["statistics_density_threshold"] < 0.0).any()) or not bool(
        np.equal(numeric["statistics_density_threshold"], numeric["statistics_density_threshold"][0]).all()
    ):
        raise ValueError(f"{operation} statistics density threshold must be one finite nonnegative value.")
    fixed_policy = parameters["fixed_policy"]
    if not isinstance(fixed_policy, Mapping) or fixed_policy.get("statistics_density_threshold") != float(
        numeric["statistics_density_threshold"][0]
    ):
        raise ValueError(f"{operation} statistics density threshold does not match the fixed producer policy.")
    metric_rows = [
        [
            int(ks[index]),
            float(numeric["stability"][index]),
            float(numeric["prediction_error"][index]),
            float(numeric["statistics_density_threshold"][index]),
            int(integers["expected_restarts"][index]),
            int(integers["completed_restarts"][index]),
            int(integers["combined_components"][index]),
        ]
        for index in range(len(ks))
    ]
    encoded_metrics = json.dumps(metric_rows, separators=(",", ":"), allow_nan=False).encode("utf-8")
    if metrics_fingerprint != f"sha256:{hashlib.sha256(encoded_metrics).hexdigest()}":
        raise ValueError(f"{operation} rank metrics fingerprint does not match the canonical table.")

    if mode["view"] == "stability_prediction_error":
        figure = Figure(figsize=(9.0, 7.0), constrained_layout=True)
        FigureCanvasAgg(figure)
        stability_axis, error_axis = figure.subplots(2, 1, sharex=True)
        stability_axis.plot(ks, numeric["stability"], marker="o", color="#4C78A8")
        stability_axis.set_ylabel("Stability")
        stability_axis.set_title("cNMF rank stability")
        stability_axis.grid(alpha=0.2)
        error_axis.plot(ks, numeric["prediction_error"], marker="o", color="#F58518")
        error_axis.set_xlabel("K")
        error_axis.set_ylabel("Prediction error")
        error_axis.set_title("cNMF rank prediction error")
        error_axis.grid(alpha=0.2)
        error_axis.set_xticks(ks)
        title = "cNMF rank stability and prediction error"
    else:
        figure = Figure(figsize=(9.0, 5.0), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        width = 0.36
        positions = np.arange(len(ks), dtype=float)
        axis.bar(positions - width / 2, integers["expected_restarts"], width, label="Expected", color="#BAB0AC")
        axis.bar(positions + width / 2, integers["completed_restarts"], width, label="Completed", color="#4C78A8")
        axis.set_xticks(positions, ks)
        axis.set_xlabel("K")
        axis.set_ylabel("Restarts")
        axis.set_title("cNMF restart completion by candidate rank")
        axis.legend()
        axis.grid(axis="y", alpha=0.2)
        title = "cNMF restart completion"
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": mode["view"],
        "view_parameters": mode,
        "candidate_ks": ks.tolist(),
        "plotted_ranks": int(len(ks)),
        "stability": numeric["stability"].tolist(),
        "prediction_error": numeric["prediction_error"].tolist(),
        "statistics_density_threshold": float(numeric["statistics_density_threshold"][0]),
        "expected_restarts": integers["expected_restarts"].tolist(),
        "completed_restarts": integers["completed_restarts"].tolist(),
        "combined_components": integers["combined_components"].tolist(),
        "input_fingerprint": fingerprint,
        "k_metrics_fingerprint_sha256": metrics_fingerprint,
        "title": title,
    }
    if _return_diagnostics:
        return png, details
    return png


def run_cnmf_rank_plot(table: Any, *, producer: dict[str, Any], view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _standalone_cnmf_rank_plot(table, producer=producer, view=view, _return_diagnostics=True)


def cnmf_rank_plot_code(*, producer: dict[str, Any], view: dict[str, Any]) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_cnmf_rank_plot)).strip()
    return f"""from __future__ import annotations

{implementation}


def plot_cnmf_rank(k_metrics):
    return _standalone_cnmf_rank_plot(k_metrics, producer={producer!r}, view={view!r})
"""


def _standalone_cnmf_programs_plot(adata, *, view=None, _return_diagnostics=False):
    import hashlib
    import io
    import json
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "cNMF Programs Plot"
    if view is None:
        mode = {"view": "usage_distributions"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant == "usage_distributions":
        if set(mode) != {"view"}:
            raise ValueError(f"{operation} usage_distributions does not accept inactive view fields.")
    elif variant == "top_genes":
        if set(mode) != {"view", "genes_per_program"}:
            raise ValueError(f"{operation} top_genes requires only genes_per_program.")
        genes_per_program = mode["genes_per_program"]
        if isinstance(genes_per_program, (bool, np.bool_)) or not isinstance(
            genes_per_program, (int, np.integer)
        ):
            raise TypeError(f"{operation} genes_per_program must be an integer.")
        if not 1 <= genes_per_program <= 50:
            raise ValueError(f"{operation} genes_per_program must be between 1 and 50.")
        mode["genes_per_program"] = int(genes_per_program)
    else:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if not isinstance(adata, AnnData):
        raise TypeError(f"{operation} requires an in-memory AnnData object.")
    if adata.isbacked:
        raise ValueError(f"{operation} does not support backed AnnData; load it into memory first.")
    if adata.n_obs < 1 or adata.n_vars < 1:
        raise ValueError(f"{operation} requires at least one observation and one feature.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError(f"{operation} requires unique observation and feature identifiers.")

    state = adata.uns.get("cnmf")
    required_state = {
        "schema_version",
        "adapter_version",
        "selected_k",
        "program_names",
        "top_genes",
        "storage_keys",
        "observation_axis_fingerprint_sha256",
        "feature_axis_fingerprint_sha256",
        "usage_fingerprint_sha256",
        "gep_scores_fingerprint_sha256",
        "top_genes_fingerprint_sha256",
    }
    if not isinstance(state, Mapping) or not required_state.issubset(state):
        raise ValueError(f"{operation} requires the complete stored cNMF consensus-program evidence schema.")
    if state["schema_version"] != 2 or state["adapter_version"] != 4:
        raise ValueError(f"{operation} requires cNMF consensus-program schema 2 from adapter 4.")

    def update_names(digest, names):
        for value in names:
            encoded = str(value).encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "little"))
            digest.update(encoded)

    def axis_fingerprint(names):
        digest = hashlib.sha256()
        update_names(digest, names)
        return f"sha256:{digest.hexdigest()}"

    def matrix_fingerprint(matrix, rows, columns):
        values = np.asarray(matrix)
        digest = hashlib.sha256()
        digest.update(str(tuple(int(value) for value in values.shape)).encode("ascii"))
        digest.update(str(values.dtype).encode("ascii"))
        update_names(digest, rows)
        update_names(digest, columns)
        digest.update(np.ascontiguousarray(values).tobytes())
        return f"sha256:{digest.hexdigest()}"

    if state["observation_axis_fingerprint_sha256"] != axis_fingerprint(adata.obs_names):
        raise ValueError(f"{operation} observation axis fingerprint does not match the stored evidence.")
    if state["feature_axis_fingerprint_sha256"] != axis_fingerprint(adata.var_names):
        raise ValueError(f"{operation} feature axis fingerprint does not match the stored evidence.")
    selected_k = state["selected_k"]
    if isinstance(selected_k, (bool, np.bool_)) or not isinstance(selected_k, (int, np.integer)) or selected_k < 2:
        raise ValueError(f"{operation} selected_k evidence must be an integer no smaller than 2.")
    program_value = state["program_names"]
    if hasattr(program_value, "tolist"):
        program_value = program_value.tolist()
    if not isinstance(program_value, (list, tuple)):
        raise ValueError(f"{operation} program axis must be an ordered string sequence.")
    programs = list(program_value)
    expected_programs = [f"cNMF_{index}" for index in range(1, int(selected_k) + 1)]
    if programs != expected_programs:
        raise ValueError(f"{operation} program axis does not match selected_k.")
    if len(programs) > 50:
        raise ValueError(f"{operation} supports at most 50 programs in one static PNG; found {len(programs)}.")
    if adata.n_obs > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 cells in one static PNG.")
    # ponytail: validate the complete stored matrices; add chunked fingerprints if larger evidence becomes practical.
    if adata.n_obs * len(programs) > 50_000_000 or adata.n_vars * len(programs) > 50_000_000:
        raise ValueError(f"{operation} supports at most 50,000,000 values per stored program matrix.")
    if state["storage_keys"] != {
        "usage": "obsm:X_cnmf_usage",
        "gep_scores": "varm:cnmf_gep_scores",
        "gep_tpm": "varm:cnmf_gep_tpm",
    }:
        raise ValueError(f"{operation} cNMF storage-key schema is invalid.")

    metadata = adata.uns.get("openbio_singlecell")
    history = metadata.get("analysis_history") if isinstance(metadata, Mapping) else None
    if not isinstance(history, Mapping) or not any(
        isinstance(entry, Mapping)
        and entry.get("operation") == "cnmf_consensus_programs"
        and isinstance(entry.get("parameters"), Mapping)
        and entry["parameters"].get("selected_k") == int(selected_k)
        for entry in history.values()
    ):
        raise ValueError(f"{operation} evidence does not match cNMF Consensus Programs provenance.")

    usage = adata.obsm.get("X_cnmf_usage")
    scores = adata.varm.get("cnmf_gep_scores")
    if not isinstance(usage, pd.DataFrame) or usage.index.tolist() != adata.obs_names.tolist() or usage.columns.tolist() != programs:
        raise ValueError(f"{operation} usage matrix must be a canonical observation-by-program DataFrame.")
    if not isinstance(scores, pd.DataFrame) or scores.index.tolist() != adata.var_names.tolist() or scores.columns.tolist() != programs:
        raise ValueError(f"{operation} GEP score matrix must be a canonical feature-by-program DataFrame.")
    usage_values = usage.to_numpy()
    score_values = scores.to_numpy()
    for values, label in ((usage_values, "usage"), (score_values, "GEP score")):
        if not bool(np.issubdtype(values.dtype, np.number)) or bool(np.iscomplexobj(values)):
            raise TypeError(f"{operation} {label} matrix must contain real numeric values.")
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} {label} matrix must contain finite values.")
    if bool((usage_values < 0).any()) or not bool(
        np.allclose(usage_values.sum(axis=1), 1.0, rtol=1e-6, atol=1e-8)
    ):
        raise ValueError(f"{operation} normalized usage rows must be nonnegative and sum to one.")
    if state["usage_fingerprint_sha256"] != matrix_fingerprint(usage_values, adata.obs_names, programs):
        raise ValueError(f"{operation} usage fingerprint does not match the stored matrix.")
    if state["gep_scores_fingerprint_sha256"] != matrix_fingerprint(score_values, adata.var_names, programs):
        raise ValueError(f"{operation} GEP score fingerprint does not match the stored matrix.")

    top_value = state["top_genes"]
    if not isinstance(top_value, Mapping) or list(top_value) != programs:
        raise ValueError(f"{operation} top-gene evidence must align exactly to the program axis.")
    available = set(adata.var_names.astype(str))
    top_genes = {}
    for program in programs:
        genes = top_value[program]
        if hasattr(genes, "tolist"):
            genes = genes.tolist()
        if not isinstance(genes, (list, tuple)) or not genes:
            raise ValueError(f"{operation} top genes for {program} must be a non-empty ordered sequence.")
        genes = list(genes)
        if (
            any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in genes)
            or len(set(genes)) != len(genes)
            or any(gene not in available for gene in genes)
        ):
            raise ValueError(f"{operation} top genes for {program} are invalid or outside the feature axis.")
        top_genes[program] = genes
    encoded_top = json.dumps(top_genes, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    if state["top_genes_fingerprint_sha256"] != f"sha256:{hashlib.sha256(encoded_top).hexdigest()}":
        raise ValueError(f"{operation} top-gene fingerprint does not match the stored evidence.")

    usage_means = usage_values.mean(axis=0)
    usage_medians = np.median(usage_values, axis=0)
    top_records = []
    if variant == "usage_distributions":
        figure = Figure(figsize=(max(7.0, min(20.0, 0.6 * len(programs) + 3.0)), 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.boxplot([usage_values[:, index] for index in range(len(programs))], tick_labels=programs, showfliers=False)
        axis.tick_params(axis="x", labelrotation=45)
        axis.set_ylabel("Normalized usage")
        axis.set_title("cNMF program usage distributions")
        axis.grid(axis="y", alpha=0.2)
        title = "cNMF program usage distributions"
    else:
        if len(programs) > 20:
            raise ValueError(
                f"{operation} top-gene view supports at most 20 programs in one static PNG; found {len(programs)}."
            )
        count = mode["genes_per_program"]
        columns = min(3, len(programs))
        rows = int(np.ceil(len(programs) / columns))
        figure = Figure(figsize=(5.2 * columns, max(3.8, 0.34 * count + 2.2) * rows), constrained_layout=True)
        FigureCanvasAgg(figure)
        axes = np.atleast_1d(figure.subplots(rows, columns, squeeze=False)).reshape(-1)
        score_by_gene = {str(gene): index for index, gene in enumerate(adata.var_names)}
        for axis, program in zip(axes, programs, strict=False):
            genes = top_genes[program][:count]
            program_index = programs.index(program)
            values = np.asarray([score_values[score_by_gene[gene], program_index] for gene in genes], dtype=float)
            positions = np.arange(len(genes), dtype=float)
            axis.barh(positions, values, color="#4C78A8")
            axis.set_yticks(positions, genes)
            axis.invert_yaxis()
            axis.set_xlabel("GEP score")
            axis.set_title(program)
            axis.grid(axis="x", alpha=0.2)
            top_records.append(
                {"program": program, "genes": genes, "gep_scores": values.tolist()}
            )
        for axis in axes[len(programs) :]:
            axis.set_visible(False)
        title = "cNMF top genes by program"
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": variant,
        "view_parameters": mode,
        "programs": programs,
        "plotted_programs": len(programs),
        "plotted_cells": int(adata.n_obs),
        "usage_means": usage_means.tolist(),
        "usage_medians": usage_medians.tolist(),
        "plotted_top_genes": sum(len(record["genes"]) for record in top_records),
        "top_gene_evidence": top_records,
        "title": title,
    }
    if _return_diagnostics:
        return png, details
    return png


def run_cnmf_programs_plot(adata: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _standalone_cnmf_programs_plot(adata, view=view, _return_diagnostics=True)


def cnmf_programs_plot_code(*, view: dict[str, Any]) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_cnmf_programs_plot)).strip()
    return f"""from __future__ import annotations

{implementation}


def plot_cnmf_programs(adata):
    return _standalone_cnmf_programs_plot(adata, view={view!r})
"""


__all__ = [
    "cnmf_programs_plot_code",
    "cnmf_rank_plot_code",
    "run_cnmf_programs_plot",
    "run_cnmf_rank_plot",
]
