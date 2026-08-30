from __future__ import annotations

import time
from typing import Any

from .analysis_utils import make_plot_result, make_summary_result, make_table_result
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_artifact_input,
    require_input_names,
    require_parameters,
    write_plot_output,
    write_table_output,
)
from .staged_state_codec import (
    VELOCITY_STATE_CODEC,
    read_velocity_state,
    write_velocity_state,
)
from .velocity_analysis import (
    run_velocity_estimate,
    run_velocity_graph,
    run_velocity_moments,
    run_velocity_prepare,
    run_velocity_ranking,
    run_velocity_recover,
    run_velocity_stream,
    validate_velocity_state,
    velocity_code,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

VELOCITY_STATE_KIND = "OPENBIO_VELOCITY_STATE"


def _anndata_input(inputs: dict[str, JSONValue], *, operation: str) -> Any:
    require_input_names(inputs, {"adata"}, operation=operation)
    return read_anndata_input(inputs)


def _state_input(inputs: dict[str, JSONValue], *, operation: str) -> Any:
    require_input_names(inputs, {"velocity_state"}, operation=operation)
    root = require_artifact_input(
        inputs,
        "velocity_state",
        kind=VELOCITY_STATE_KIND,
        codec=VELOCITY_STATE_CODEC,
    )
    return read_velocity_state(root)


def _summary_result(
    summary: dict[str, Any],
    *,
    title: str,
    operation: str,
    started_at: float,
    cells: int,
    genes: int,
) -> Any:
    key_results = summary["key_results"]
    return make_summary_result(
        summary=summary,
        title=title,
        operation=operation,
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=int(key_results.get("input_cells", cells)),
        input_genes=int(key_results.get("input_genes", genes)),
        started_at=started_at,
    )


def _state_records(
    context: OperationContext,
    state: Any,
    report: Any,
    code: str,
) -> list[dict[str, JSONValue]]:
    root = context.create_output_directory("velocity_state")
    write_velocity_state(root, state)
    return analysis_outputs(
        report,
        code,
        {
            "type": "artifact",
            "name": "velocity_state",
            "kind": VELOCITY_STATE_KIND,
            "codec": VELOCITY_STATE_CODEC,
            "payload": root.relative_to(context.output_root).as_posix(),
        },
    )


@register_operation("openbio.node.velocityfilterandnormalize")
def velocity_filter_and_normalize(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Prepare Velocity Abundances"
    require_parameters(
        parameters,
        {
            "spliced_layer",
            "unspliced_layer",
            "min_shared_counts",
            "min_shared_cells",
            "normalization_target",
            "overwrite_existing",
        },
        operation=operation,
    )
    adata = _anndata_input(inputs, operation=operation)
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    state = run_velocity_prepare(adata, _owned=True, **parameters)
    summary = state.summary
    report = _summary_result(
        summary,
        title="Velocity abundance preparation summary",
        operation="prepare_velocity_abundances",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    return _state_records(context, state, report, velocity_code("prepare", **parameters))


@register_operation("openbio.node.velocitymoments")
def velocity_moments(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Velocity Moments"
    require_parameters(
        parameters,
        {"neighbors_key", "mode", "max_dense_gib", "overwrite_existing"},
        operation=operation,
    )
    state = _state_input(inputs, operation=operation)
    adata, _summary, _metadata, _portable = validate_velocity_state(state, allowed_stages=("prepared",), _owned=True)
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    result = run_velocity_moments(state, _owned=True, **parameters)
    summary = result.summary
    report = _summary_result(
        summary,
        title="Velocity moments summary",
        operation="velocity_moments",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    return _state_records(context, result, report, velocity_code("moments", **parameters))


@register_operation("openbio.node.estimatevelocity")
def estimate_velocity(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Estimate RNA Velocity"
    require_parameters(
        parameters,
        {"mode", "vkey", "min_r2", "min_likelihood", "overwrite_existing"},
        operation=operation,
    )
    state = _state_input(inputs, operation=operation)
    adata, _summary, _metadata, _portable = validate_velocity_state(state, _owned=True)
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    result = run_velocity_estimate(state, _owned=True, **parameters)
    summary = result.summary
    report = _summary_result(
        summary,
        title="RNA velocity estimate summary",
        operation="estimate_rna_velocity",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    return _state_records(context, result, report, velocity_code("estimate", **parameters))


@register_operation("openbio.node.recoverdynamics")
def recover_dynamics(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Recover RNA Velocity Dynamics"
    require_parameters(
        parameters,
        {
            "gene_selection",
            "n_top_genes",
            "max_iter",
            "n_jobs",
            "max_dense_gib",
            "overwrite_existing",
        },
        operation=operation,
    )
    state = _state_input(inputs, operation=operation)
    adata, _summary, _metadata, _portable = validate_velocity_state(state, _owned=True)
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    result = run_velocity_recover(state, _owned=True, **parameters)
    summary = result.summary
    report = _summary_result(
        summary,
        title="Recovered RNA velocity dynamics summary",
        operation="recover_velocity_dynamics",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    return _state_records(context, result, report, velocity_code("recover", **parameters))


@register_operation("openbio.node.velocitygraph")
def velocity_graph(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Velocity Directed-Correlation Graph"
    require_parameters(
        parameters,
        {"vkey", "xkey", "mode_neighbors", "n_jobs", "overwrite_existing"},
        operation=operation,
    )
    state = _state_input(inputs, operation=operation)
    adata, _summary, _metadata, _portable = validate_velocity_state(
        state, allowed_stages=("velocity_estimated",), _owned=True
    )
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    result = run_velocity_graph(state, _owned=True, **parameters)
    summary = result.summary
    report = _summary_result(
        summary,
        title="Velocity directed-correlation graph summary",
        operation="build_velocity_graph",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    return _state_records(context, result, report, velocity_code("graph", **parameters))


@register_operation("openbio.node.velocitygeneranking")
def velocity_gene_ranking(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Recovered Dynamics Fit Ranking"
    require_parameters(
        parameters,
        {"top_n", "include_failed", "max_output_rows"},
        operation=operation,
    )
    state = _state_input(inputs, operation=operation)
    adata, _summary, _metadata, _portable = validate_velocity_state(state, _owned=True)
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    table, summary = run_velocity_ranking(state, _owned=True, **parameters)
    # Failed fits are scientific nulls, not non-finite JSON numbers. The table
    # is bounded and may be materialized with pandas nullable dtypes for JSONL.
    table = table.convert_dtypes()
    table_result = make_table_result(
        title="Recovered dynamics fit ranking",
        operation="rank_recovered_dynamics",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        table=table,
    )
    report = _summary_result(
        summary,
        title="Recovered dynamics fit ranking summary",
        operation="rank_recovered_dynamics",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    return analysis_outputs(
        report,
        velocity_code("ranking", **parameters),
        write_table_output(context, table_result, kind="OPENBIO_SINGLE_CELL_TABLE"),
    )


@register_operation("openbio.node.velocitystreamplot")
def velocity_stream_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[dict[str, JSONValue]]:
    operation = "Velocity Stream Plot"
    require_parameters(
        parameters,
        {"basis", "color_key", "density", "smooth", "min_mass"},
        operation=operation,
    )
    state = _state_input(inputs, operation=operation)
    adata, _summary, _metadata, _portable = validate_velocity_state(
        state, allowed_stages=("velocity_graph",), _owned=True
    )
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    started_at = time.perf_counter()
    png, summary = run_velocity_stream(state, _owned=True, **parameters)
    plot_result = make_plot_result(
        title=f"RNA velocity stream on {parameters['basis']}",
        operation="render_velocity_stream",
        parameters=summary["parameters"],
        description=summary["results"],
        warnings=summary["warnings"],
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        png=png,
    )
    report = _summary_result(
        summary,
        title="Velocity stream visualization summary",
        operation="render_velocity_stream",
        started_at=started_at,
        cells=cells,
        genes=genes,
    )
    return analysis_outputs(
        report,
        velocity_code("stream", **parameters),
        write_plot_output(context, plot_result, kind="OPENBIO_SINGLE_CELL_PLOT"),
    )


__all__ = [
    "VELOCITY_STATE_KIND",
    "estimate_velocity",
    "recover_dynamics",
    "velocity_filter_and_normalize",
    "velocity_gene_ranking",
    "velocity_graph",
    "velocity_moments",
    "velocity_stream_plot",
]
