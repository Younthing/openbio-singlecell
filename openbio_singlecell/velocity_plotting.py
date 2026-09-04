from __future__ import annotations

from typing import Any

from .velocity_analysis import (
    _VELOCITY_COMMON_STANDALONE_CONSTANTS,
    _VELOCITY_COMMON_STANDALONE_HELPERS,
    _velocity_source,
)
from .velocity_portable import _vs_series_fingerprint, validate_portable_velocity_state


def _standalone_velocity_dynamics_plot(
    velocity_adata,
    *,
    gene="",
    view=None,
    _return_details=False,
):
    import io

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from scipy import sparse

    operation = "Velocity Dynamics Plot"
    state = validate_portable_velocity_state(
        velocity_adata,
        allowed_stages=("dynamics_recovered", "velocity_estimated", "velocity_graph"),
    )
    if state["dynamics_fit_fingerprint_sha256"] is None:
        raise ValueError(f"{operation} requires a state with verified recovered-dynamics evidence.")
    mode = {"view": "phase_portrait"} if view is None else view
    if not isinstance(mode, dict) or set(mode) != {"view"} or mode.get("view") not in {
        "phase_portrait",
        "fit_loss",
    }:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    if not isinstance(gene, str):
        raise TypeError(f"{operation} gene must be a string.")
    gene = gene.strip()
    required = {
        "openbio_dynamics_fit_success",
        "fit_likelihood",
        "fit_r2",
        "fit_alpha",
        "fit_beta",
        "fit_gamma",
        "fit_scaling",
        "fit_t_",
    }
    missing = sorted(required - set(velocity_adata.var.columns))
    if missing:
        raise ValueError(f"{operation} is missing canonical kinetic fields: {missing!r}.")
    success = velocity_adata.var["openbio_dynamics_fit_success"]
    if success.dtype != bool or bool(success.isna().any()):
        raise TypeError(f"{operation} fit-success axis must be complete boolean evidence.")
    successful = np.flatnonzero(success.to_numpy(dtype=bool))
    if successful.size == 0:
        raise ValueError(f"{operation} has no successful recovered fit to display.")
    likelihood = velocity_adata.var["fit_likelihood"].to_numpy(dtype=float)
    if not bool(np.isfinite(likelihood[successful]).all()) or bool((likelihood[successful] < 0).any()):
        raise ValueError(f"{operation} successful fit likelihoods must be finite and nonnegative.")
    if gene:
        if gene not in velocity_adata.var_names:
            raise ValueError(f"{operation} gene was not found exactly: {gene!r}.")
        gene_index = int(velocity_adata.var_names.get_loc(gene))
        selection = "explicit"
    else:
        gene_index = int(successful[np.lexsort((successful, -likelihood[successful]))[0]])
        gene = str(velocity_adata.var_names[gene_index])
        selection = "highest_stored_fit_likelihood"
    if not bool(success.iloc[gene_index]):
        raise ValueError(f"{operation} gene {gene!r} does not have a successful recovered fit.")
    kinetic = {}
    for key in ("fit_likelihood", "fit_r2", "fit_alpha", "fit_beta", "fit_gamma", "fit_scaling", "fit_t_"):
        value = float(velocity_adata.var[key].iloc[gene_index])
        if not np.isfinite(value):
            raise ValueError(f"{operation} gene {gene!r} has non-finite {key!r} evidence.")
        kinetic[key] = value
    if kinetic["fit_likelihood"] < 0 or any(kinetic[key] <= 0 for key in ("fit_alpha", "fit_beta", "fit_gamma", "fit_scaling")):
        raise ValueError(f"{operation} gene {gene!r} has invalid stored kinetic parameters.")

    details = {
        "view": mode["view"],
        "view_parameters": dict(mode),
        "gene": gene,
        "gene_index": gene_index,
        "gene_selection": selection,
        **kinetic,
        "dynamics_fit_fingerprint_sha256": state["dynamics_fit_fingerprint_sha256"],
        "upstream_state_fingerprint_sha256": state["state_fingerprint_sha256"],
    }
    if mode["view"] == "phase_portrait":
        if int(velocity_adata.n_obs) > 200_000:
            raise ValueError(f"{operation} phase portrait supports at most 200,000 cells in one static PNG.")
        for key in ("Ms", "Mu", "fit_t"):
            if key not in velocity_adata.layers:
                raise ValueError(f"{operation} is missing canonical layer {key!r}.")
            if tuple(velocity_adata.layers[key].shape) != tuple(velocity_adata.shape):
                raise ValueError(f"{operation} layer {key!r} is not aligned to the cell-by-gene axes.")
        def column(key):
            values = velocity_adata.layers[key][:, gene_index]
            return np.asarray(values.toarray() if sparse.issparse(values) else values, dtype=float).ravel()

        spliced = column("Ms")
        unspliced = column("Mu")
        fitted_time = column("fit_t")
        if not bool(np.isfinite(spliced).all() and np.isfinite(unspliced).all() and np.isfinite(fitted_time).all()):
            raise ValueError(f"{operation} phase-portrait axes contain non-finite stored values.")
        figure = Figure(figsize=(7.0, 5.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        marks = axis.scatter(spliced, unspliced, c=fitted_time, cmap="viridis", s=18, alpha=0.8, linewidths=0)
        axis.set_xlabel("Spliced moment (Ms)")
        axis.set_ylabel("Unspliced moment (Mu)")
        axis.set_title(f"Stored phase portrait: {gene}")
        axis.grid(alpha=0.2)
        figure.colorbar(marks, ax=axis, label="Stored fitted time (fit_t)")
        title = f"Velocity dynamics phase portrait: {gene}"
        details.update(
            {
                "plotted_cells": int(velocity_adata.n_obs),
                "spliced_range": [float(spliced.min()), float(spliced.max())],
                "unspliced_range": [float(unspliced.min()), float(unspliced.max())],
                "fitted_time_range": [float(fitted_time.min()), float(fitted_time.max())],
            }
        )
    else:
        if "loss" not in velocity_adata.varm:
            raise ValueError(f"{operation} is missing canonical varm['loss'] fit evidence.")
        loss_matrix = np.asarray(velocity_adata.varm["loss"])
        if loss_matrix.ndim != 2 or loss_matrix.shape[0] != int(velocity_adata.n_vars):
            raise ValueError(f"{operation} loss matrix is not aligned to the feature axis.")
        loss = np.asarray(loss_matrix[gene_index], dtype=float).ravel()
        if loss.size == 0 or loss.size > 10_000 or not bool(np.isfinite(loss).all()):
            raise ValueError(f"{operation} selected-gene loss evidence must contain 1 to 10,000 finite points.")
        figure = Figure(figsize=(7.0, 4.8), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        iterations = np.arange(1, loss.size + 1)
        axis.plot(iterations, loss, marker="o", color="#4C78A8")
        axis.set_xlabel("Stored fit-loss step")
        axis.set_ylabel("Fit loss")
        axis.set_title(f"Stored recovered-dynamics loss: {gene}")
        axis.grid(alpha=0.2)
        title = f"Velocity dynamics fit loss: {gene}"
        details.update(
            {
                "plotted_loss_points": int(loss.size),
                "initial_loss": float(loss[0]),
                "final_loss": float(loss[-1]),
                "minimum_loss": float(loss.min()),
            }
        )
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details["title"] = title
    if _return_details:
        return png, details
    return png


def run_velocity_dynamics_plot(
    velocity_adata: Any,
    *,
    gene: str = "",
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_velocity_dynamics_plot(
        velocity_adata,
        gene=gene,
        view=view,
        _return_details=True,
    )


def velocity_dynamics_plot_code(*, gene: str, view: dict[str, Any]) -> str:
    helpers = (
        *_VELOCITY_COMMON_STANDALONE_HELPERS,
        validate_portable_velocity_state,
        _standalone_velocity_dynamics_plot,
    )
    source = _velocity_source(
        helpers=helpers,
        constants=(*_VELOCITY_COMMON_STANDALONE_CONSTANTS, "_STATE_METADATA_KEYS"),
        backend_operation=None,
    )
    wrapper = f"""

def plot_velocity_dynamics(velocity_adata):
    return _standalone_velocity_dynamics_plot(velocity_adata, gene={gene!r}, view={view!r})
"""
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, "<plot_velocity_dynamics>", "exec")
    return code


def _standalone_velocity_gene_ranking_plot(table, *, producer, view=None, _return_details=False):
    import io

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Velocity Gene Ranking Plot"
    if not isinstance(producer, dict) or set(producer) != {
        "operation",
        "parameters",
        "input_cells",
        "input_genes",
        "random_seed",
    }:
        raise ValueError(f"{operation} requires exact producer provenance.")
    if producer["operation"] != "rank_recovered_dynamics" or not isinstance(producer["parameters"], dict):
        raise ValueError(f"{operation} requires the Recovered Dynamics Fit Ranking producer.")
    parameters = producer["parameters"]
    required_parameters = {
        "top_n",
        "include_failed",
        "max_output_rows",
        "dynamics_fit_fingerprint_sha256",
        "upstream_state_fingerprint_sha256",
        "table_fingerprint_sha256",
    }
    if set(parameters) != required_parameters:
        raise ValueError(f"{operation} producer parameters have an invalid exact schema.")
    for key in (
        "dynamics_fit_fingerprint_sha256",
        "upstream_state_fingerprint_sha256",
        "table_fingerprint_sha256",
    ):
        value = parameters[key]
        if not isinstance(value, str) or len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
            raise ValueError(f"{operation} producer {key!r} is not a SHA-256 fingerprint.")
    mode = {"view": "fit_quality", "max_genes": 30} if view is None else view
    if not isinstance(mode, dict) or set(mode) != {"view", "max_genes"} or mode.get("view") not in {
        "fit_quality",
        "kinetic_parameters",
    }:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    max_genes = mode["max_genes"]
    if isinstance(max_genes, bool) or not isinstance(max_genes, int) or not 1 <= max_genes <= 100:
        raise ValueError(f"{operation} max_genes must be an integer between 1 and 100.")
    columns = [
        "rank",
        "gene",
        "fit_status",
        "fit_likelihood",
        "fit_r2",
        "velocity_gene",
        "fit_alpha",
        "fit_beta",
        "fit_gamma",
        "fit_scaling",
        "fit_switch_time",
    ]
    if not isinstance(table, pd.DataFrame) or list(table.columns) != columns:
        raise ValueError(f"{operation} requires the canonical exact ranking columns.")
    if len(table) == 0 or len(table) > int(parameters["max_output_rows"]):
        raise ValueError(f"{operation} ranking row count conflicts with the producer bound.")
    observed_fingerprint = _vs_series_fingerprint(table, columns, label="dynamics-ranking-plot")
    if observed_fingerprint != parameters["table_fingerprint_sha256"]:
        raise ValueError(f"{operation} table fingerprint does not match the producer evidence.")
    genes = table["gene"].tolist()
    if any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in genes) or len(set(genes)) != len(genes):
        raise ValueError(f"{operation} gene axis must contain unique canonical strings.")
    if bool(table["fit_status"].isna().any()) or not set(table["fit_status"]).issubset({"fitted", "failed"}):
        raise ValueError(f"{operation} fit_status must use the canonical fitted/failed values.")
    fitted = table.loc[table["fit_status"] == "fitted"].copy()
    failed = table.loc[table["fit_status"] == "failed"]
    if fitted.empty or list(table["fit_status"]) != ["fitted"] * len(fitted) + ["failed"] * len(failed):
        raise ValueError(f"{operation} requires ranked fitted rows before any failed-fit rows.")
    ranks = fitted["rank"].to_numpy(dtype=int)
    if not np.array_equal(ranks, np.arange(1, len(fitted) + 1)):
        raise ValueError(f"{operation} fitted rows must preserve the canonical consecutive rank axis.")
    if bool(failed["rank"].notna().any()):
        raise ValueError(f"{operation} failed-fit rows cannot carry ranks.")
    likelihood = fitted["fit_likelihood"].to_numpy(dtype=float)
    if not bool(np.isfinite(likelihood).all()) or bool((likelihood < 0).any()) or bool(np.any(np.diff(likelihood) > 0)):
        raise ValueError(f"{operation} fit likelihood must be finite, nonnegative, and rank-ordered.")
    if len(fitted) > int(parameters["top_n"]):
        raise ValueError(f"{operation} fitted rows exceed the producer top_n bound.")
    if not bool(parameters["include_failed"]) and not failed.empty:
        raise ValueError(f"{operation} contains failed fits excluded by producer policy.")
    if any(not isinstance(value, (bool, np.bool_)) for value in table["velocity_gene"]):
        raise TypeError(f"{operation} velocity_gene must contain complete booleans.")
    numeric_columns = [
        "fit_likelihood",
        "fit_r2",
        "fit_alpha",
        "fit_beta",
        "fit_gamma",
        "fit_scaling",
        "fit_switch_time",
    ]
    if any(bool(failed[column].notna().any()) for column in numeric_columns):
        raise ValueError(f"{operation} failed-fit rows must retain scientific nulls for numeric evidence.")
    for column in numeric_columns:
        values = fitted[column].to_numpy(dtype=float, na_value=np.nan)
        if column == "fit_r2":
            if bool(np.isinf(values).any()):
                raise ValueError(f"{operation} fit_r2 cannot contain infinity.")
        elif not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} fitted {column!r} values must be finite.")
    for column in ("fit_alpha", "fit_beta", "fit_gamma", "fit_scaling"):
        if bool((fitted[column].to_numpy(dtype=float) <= 0).any()):
            raise ValueError(f"{operation} fitted {column!r} values must be positive.")
    selected = fitted.iloc[:max_genes]
    selected_genes = selected["gene"].tolist()
    selected_ranks = [int(value) for value in selected["rank"]]
    selected_likelihood = [float(value) for value in selected["fit_likelihood"]]
    details = {
        "view": mode["view"],
        "view_parameters": dict(mode),
        "genes": selected_genes,
        "ranks": selected_ranks,
        "fit_likelihood": selected_likelihood,
        "plotted_genes": len(selected),
        "available_ranked_genes": len(fitted),
        "stored_failed_fits": len(failed),
        "table_fingerprint_sha256": observed_fingerprint,
        "dynamics_fit_fingerprint_sha256": parameters["dynamics_fit_fingerprint_sha256"],
        "upstream_state_fingerprint_sha256": parameters["upstream_state_fingerprint_sha256"],
    }
    positions = np.arange(len(selected), dtype=float)
    if mode["view"] == "fit_quality":
        r2 = selected["fit_r2"].to_numpy(dtype=float, na_value=np.nan)
        figure = Figure(figsize=(9.0, max(4.5, min(18.0, 0.35 * len(selected) + 2.0))), constrained_layout=True)
        FigureCanvasAgg(figure)
        likelihood_axis, r2_axis = figure.subplots(1, 2, sharey=True)
        likelihood_axis.barh(positions, selected_likelihood, color="#4C78A8")
        likelihood_axis.set_yticks(positions, selected_genes)
        likelihood_axis.invert_yaxis()
        likelihood_axis.set_xlabel("Stored fit likelihood")
        likelihood_axis.set_title("Recovered-fit rank")
        likelihood_axis.grid(axis="x", alpha=0.2)
        finite_r2 = np.isfinite(r2)
        r2_axis.scatter(r2[finite_r2], positions[finite_r2], color="#F58518")
        r2_axis.axvline(0.0, color="black", linewidth=0.8, alpha=0.5)
        r2_axis.set_xlabel("Stored fit R²")
        r2_axis.set_title("Fit evidence")
        r2_axis.grid(axis="x", alpha=0.2)
        title = "Recovered dynamics fit-quality ranking"
        details["fit_r2"] = [None if not np.isfinite(value) else float(value) for value in r2]
        details["plotted_fit_r2"] = int(finite_r2.sum())
    else:
        parameter_columns = ["fit_alpha", "fit_beta", "fit_gamma", "fit_scaling"]
        values = selected[parameter_columns].to_numpy(dtype=float)
        logged = np.log10(values)
        figure = Figure(figsize=(7.0, max(4.5, min(18.0, 0.35 * len(selected) + 2.0))), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        image = axis.imshow(logged, aspect="auto", interpolation="nearest", cmap="viridis")
        axis.set_xticks(np.arange(len(parameter_columns)), [name.removeprefix("fit_") for name in parameter_columns])
        axis.set_yticks(positions, selected_genes)
        axis.set_xlabel("Stored kinetic parameter")
        axis.set_ylabel("Upstream rank order")
        axis.set_title("Recovered kinetic parameters (log10)")
        figure.colorbar(image, ax=axis, label="log10 fitted value")
        title = "Recovered dynamics kinetic-parameter evidence"
        details["kinetic_parameter_columns"] = parameter_columns
        details["kinetic_parameter_values"] = values.tolist()
        details["kinetic_parameter_log10_values"] = logged.tolist()
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details["title"] = title
    if _return_details:
        return png, details
    return png


def run_velocity_gene_ranking_plot(
    table: Any,
    *,
    producer: dict[str, Any],
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_velocity_gene_ranking_plot(
        table,
        producer=producer,
        view=view,
        _return_details=True,
    )


def velocity_gene_ranking_plot_code(*, producer: dict[str, Any], view: dict[str, Any]) -> str:
    helpers = (
        *_VELOCITY_COMMON_STANDALONE_HELPERS,
        _standalone_velocity_gene_ranking_plot,
    )
    source = _velocity_source(
        helpers=helpers,
        constants=_VELOCITY_COMMON_STANDALONE_CONSTANTS,
        backend_operation=None,
    )
    wrapper = f"""

def plot_velocity_gene_ranking(table):
    return _standalone_velocity_gene_ranking_plot(table, producer={producer!r}, view={view!r})
"""
    code = f"{source.rstrip()}\n{wrapper}"
    compile(code, "<plot_velocity_gene_ranking>", "exec")
    return code


__all__ = [
    "run_velocity_dynamics_plot",
    "run_velocity_gene_ranking_plot",
    "velocity_dynamics_plot_code",
    "velocity_gene_ranking_plot_code",
]
