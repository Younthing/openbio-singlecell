from __future__ import annotations

import math
import time
from collections.abc import Mapping
from typing import Any

from .analysis_utils import finish_adata, make_summary_result
from .cell_cycle import analyze_cell_cycle_score, cell_cycle_code
from .expression_source import _CELL_CYCLE_SPEC, DynamicExpressionSource
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
)
from .trajectory_analysis import (
    analyze_diffusion_map,
    analyze_dpt,
    analyze_paga,
    diffusion_map_code,
    dpt_code,
    paga_code,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

CELL_CYCLE_EXPRESSION_SOURCE = _CELL_CYCLE_SPEC


def parse_dpt_root_value(value_type: object, value: object) -> str | int | float:
    if value_type not in {"string", "integer", "number"}:
        raise ValueError("DPT root_value_type must be 'string', 'integer', or 'number'.")
    if not isinstance(value, str):
        raise TypeError("DPT root_value must be entered as text and interpreted by root_value_type.")
    if not value or value != value.strip():
        raise ValueError("DPT root_value must be nonempty and free of surrounding whitespace.")
    if value_type == "string":
        return value
    if value_type == "integer":
        try:
            parsed = int(value, 10)
        except ValueError as error:
            raise ValueError("DPT integer root_value must use canonical base-10 syntax.") from error
        if str(parsed) != value:
            raise ValueError("DPT integer root_value must use canonical base-10 syntax.")
        return parsed
    try:
        parsed_number = float(value)
    except ValueError as error:
        raise ValueError("DPT number root_value must be a finite numeric literal.") from error
    if not math.isfinite(parsed_number):
        raise ValueError("DPT number root_value must be a finite numeric literal.")
    return parsed_number


def _input_adata(inputs: dict[str, JSONValue]) -> Any:
    return read_anndata_input(inputs)


def _records(context: OperationContext, result: tuple[Any, Any, str]) -> list[JSONValue]:
    output, report, code = result
    return analysis_outputs(report, code, write_anndata_output(context, output))


def cell_cycle_score_owned(
    adata: Any,
    *,
    source: DynamicExpressionSource | None = None,
    gene_set_source: Mapping[str, object] | None = None,
    organism: str = "human",
    output_prefix: str = "cell_cycle",
    overwrite_existing: bool = False,
    random_seed: int = 0,
) -> tuple[Any, Any, str]:
    expression = CELL_CYCLE_EXPRESSION_SOURCE.resolve(adata, source)
    if gene_set_source is None:
        gene_set_source = {"gene_set_source": "regev_human_97"}
    if not isinstance(gene_set_source, Mapping):
        raise TypeError("Cell-cycle gene_set_source must be a DynamicCombo value.")
    mode = gene_set_source.get("gene_set_source")
    if mode not in {"regev_human_97", "custom"}:
        raise ValueError(f"Unsupported cell-cycle gene-set source: {mode!r}.")
    custom_s = gene_set_source.get("s_genes") if mode == "custom" else None
    custom_g2m = gene_set_source.get("g2m_genes") if mode == "custom" else None
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, summary = analyze_cell_cycle_score(
        adata,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        gene_set_source=str(mode),
        organism=organism,
        s_genes=custom_s,
        g2m_genes=custom_g2m,
        output_prefix=output_prefix,
        overwrite_existing=overwrite_existing,
        random_seed=random_seed,
    )
    parameters = dict(summary["parameters"])
    warnings = [str(warning) for warning in summary["warnings"]]
    finish_adata(
        output,
        "cell_cycle_score",
        parameters,
        cells,
        genes,
        started_at,
        random_seed=random_seed,
        warnings=warnings,
    )
    report = make_summary_result(
        summary=summary,
        title="Cell-cycle program score summary",
        operation="cell_cycle_score",
        parameters=parameters,
        description=str(summary["results"]),
        warnings=warnings,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    code = cell_cycle_code(
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        gene_set_source=str(mode),
        organism=organism,
        s_genes=custom_s,
        g2m_genes=custom_g2m,
        output_prefix=output_prefix,
        overwrite_existing=overwrite_existing,
        random_seed=random_seed,
    )
    return output, report, code


def diffusion_map_owned(
    adata: Any,
    *,
    neighbors_key: str = "neighbors",
    n_comps: int = 15,
    overwrite_existing: bool = False,
    random_seed: int = 0,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, summary = analyze_diffusion_map(
        adata,
        neighbors_key=neighbors_key,
        n_comps=n_comps,
        overwrite_existing=overwrite_existing,
        random_seed=random_seed,
    )
    parameters = dict(summary["parameters"])
    warnings = [str(value) for value in summary["warnings"]]
    finish_adata(
        output,
        "diffusion_map",
        parameters,
        cells,
        genes,
        started_at,
        random_seed=random_seed,
        warnings=warnings,
    )
    report = make_summary_result(
        summary=summary,
        title="Diffusion Map summary",
        operation="diffusion_map",
        parameters=parameters,
        description=str(summary["results"]),
        warnings=warnings,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    code = diffusion_map_code(
        neighbors_key=str(parameters["neighbors_key"]),
        n_comps=int(parameters["n_comps"]),
        overwrite_existing=bool(parameters["overwrite_existing"]),
        random_seed=int(parameters["random_seed"]),
    )
    return output, report, code


def paga_owned(
    adata: Any,
    *,
    groupby: str = "leiden",
    neighbors_key: str = "neighbors",
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, summary = analyze_paga(
        adata,
        groupby=groupby,
        neighbors_key=neighbors_key,
        overwrite_existing=overwrite_existing,
    )
    parameters = dict(summary["parameters"])
    warnings = [str(value) for value in summary["warnings"]]
    finish_adata(output, "paga", parameters, cells, genes, started_at, warnings=warnings)
    report = make_summary_result(
        summary=summary,
        title="PAGA connectivity summary",
        operation="paga",
        parameters=parameters,
        description=str(summary["results"]),
        warnings=warnings,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    code = paga_code(
        groupby=str(parameters["groupby"]),
        neighbors_key=str(parameters["neighbors_key"]),
        overwrite_existing=bool(parameters["overwrite_existing"]),
    )
    return output, report, code


def dpt_owned(
    adata: Any,
    *,
    neighbors_key: str = "neighbors",
    root_mode: Mapping[str, object] | None = None,
    n_dcs: int = 10,
    overwrite_existing: bool = False,
) -> tuple[Any, Any, str]:
    if root_mode is None:
        root_mode = {"root_mode": "cell_id", "root_cell_id": ""}
    if not isinstance(root_mode, Mapping):
        raise TypeError("DPT root_mode must be a DynamicCombo value.")
    mode = root_mode.get("root_mode")
    if mode not in {"cell_id", "group_medoid"}:
        raise ValueError(f"Unsupported DPT root mode: {mode!r}.")
    root_cell_id = root_mode.get("root_cell_id", "") if mode == "cell_id" else ""
    root_column = root_mode.get("root_column", "") if mode == "group_medoid" else ""
    root_value: Any = ""
    if mode == "group_medoid":
        root_value = parse_dpt_root_value(root_mode.get("root_value_type"), root_mode.get("root_value"))
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, summary = analyze_dpt(
        adata,
        neighbors_key=neighbors_key,
        root_mode=str(mode),
        root_cell_id=root_cell_id,
        root_column=root_column,
        root_value=root_value,
        n_dcs=n_dcs,
        overwrite_existing=overwrite_existing,
    )
    parameters = dict(summary["parameters"])
    warnings = [str(value) for value in summary["warnings"]]
    finish_adata(output, "dpt", parameters, cells, genes, started_at, warnings=warnings)
    report = make_summary_result(
        summary=summary,
        title="Diffusion Pseudotime summary",
        operation="dpt",
        parameters=parameters,
        description=str(summary["results"]),
        warnings=warnings,
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
    )
    code = dpt_code(
        neighbors_key=str(parameters["neighbors_key"]),
        root_mode=str(parameters["root_mode"]),
        root_cell_id=str(parameters["root_cell_id"]),
        root_column=str(parameters["root_column"] or ""),
        root_value=root_value,
        n_dcs=int(parameters["n_dcs"]),
        overwrite_existing=bool(parameters["overwrite_existing"]),
    )
    return output, report, code


@register_operation("openbio.node.cellcyclescore")
def cell_cycle_score(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Cell Cycle Score")
    require_parameters(
        parameters,
        {"source", "gene_set_source", "organism", "output_prefix", "overwrite_existing", "random_seed"},
        operation="Cell Cycle Score",
    )
    return _records(context, cell_cycle_score_owned(_input_adata(inputs), **parameters))


@register_operation("openbio.node.diffusionmap")
def diffusion_map(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Diffusion Map")
    require_parameters(
        parameters,
        {"neighbors_key", "n_comps", "overwrite_existing", "random_seed"},
        operation="Diffusion Map",
    )
    return _records(context, diffusion_map_owned(_input_adata(inputs), **parameters))


@register_operation("openbio.node.paga")
def paga(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="PAGA")
    require_parameters(parameters, {"groupby", "neighbors_key", "overwrite_existing"}, operation="PAGA")
    return _records(context, paga_owned(_input_adata(inputs), **parameters))


@register_operation("openbio.node.dpt")
def dpt(
    context: OperationContext, inputs: dict[str, JSONValue], parameters: dict[str, JSONValue]
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Diffusion Pseudotime")
    require_parameters(
        parameters,
        {"neighbors_key", "root_mode", "n_dcs", "overwrite_existing"},
        operation="Diffusion Pseudotime",
    )
    return _records(context, dpt_owned(_input_adata(inputs), **parameters))


__all__ = [
    "CELL_CYCLE_EXPRESSION_SOURCE",
    "cell_cycle_score",
    "cell_cycle_score_owned",
    "diffusion_map",
    "diffusion_map_owned",
    "dpt",
    "dpt_owned",
    "paga",
    "paga_owned",
    "parse_dpt_root_value",
]
