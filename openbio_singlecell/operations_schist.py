from __future__ import annotations

import time
from typing import Any

from . import PLUGIN_VERSION
from .analysis_utils import finish_adata, make_summary_result
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
)
from .schist_analysis import run_schist_nested_model, schist_nested_model_code
from .worker_protocol import JSONValue, OperationContext, register_operation


def schist_nested_model_owned(
    adata: Any,
    *,
    random_seed: int = 123,
    neighbors_key: str = "neighbors",
    key_added: str = "nsbm",
    posterior_samples: int = 100,
    degree_correction: bool = True,
    overwrite_existing: bool = False,
    max_working_memory_gib: float = 8.0,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    cells, genes = int(adata.n_obs), int(adata.n_vars)
    output, summary = run_schist_nested_model(
        adata,
        random_seed=random_seed,
        neighbors_key=neighbors_key,
        key_added=key_added,
        posterior_samples=posterior_samples,
        degree_correction=degree_correction,
        overwrite_existing=overwrite_existing,
        max_working_memory_gib=max_working_memory_gib,
        openbio_version=PLUGIN_VERSION,
    )
    code = schist_nested_model_code(
        random_seed=random_seed,
        neighbors_key=neighbors_key,
        key_added=key_added,
        posterior_samples=posterior_samples,
        degree_correction=degree_correction,
        overwrite_existing=overwrite_existing,
        max_working_memory_gib=max_working_memory_gib,
        openbio_version=PLUGIN_VERSION,
    )
    finish_adata(
        output,
        "schist_nested_sbm_hierarchy",
        dict(summary["parameters"]),
        cells,
        genes,
        started_at,
        random_seed=random_seed,
        warnings=list(summary["warnings"]),
    )
    report = make_summary_result(
        summary=summary,
        title="Schist Nested-SBM hierarchy summary",
        operation="schist_nested_sbm_hierarchy",
        parameters=dict(summary["parameters"]),
        description=str(summary["results"]),
        warnings=list(summary["warnings"]),
        input_cells=cells,
        input_genes=genes,
        started_at=started_at,
        random_seed=random_seed,
    )
    return output, report, code


@register_operation("openbio.node.schistnestedmodel")
def schist_nested_model(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Schist Nested-SBM Hierarchy")
    require_parameters(
        parameters,
        {
            "random_seed",
            "neighbors_key",
            "key_added",
            "posterior_samples",
            "degree_correction",
            "overwrite_existing",
            "max_working_memory_gib",
        },
        operation="Schist Nested-SBM Hierarchy",
    )
    output, report, code = schist_nested_model_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_anndata_output(context, output))


__all__ = ["schist_nested_model", "schist_nested_model_owned"]
