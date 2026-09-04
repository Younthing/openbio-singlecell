from __future__ import annotations

import time
from typing import Any

from . import PLUGIN_VERSION
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata, make_plot_result, make_summary_result
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
)
from .schist_analysis import run_schist_nested_model, schist_nested_model_code
from .schist_plotting import _standalone_schist_hierarchy_plot, schist_hierarchy_plot_code
from .worker_protocol import JSONValue, OperationContext, register_operation

PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
SCHIST_REFERENCE = AnalysisReference(
    citation=(
        "Morelli L, Giansanti V, Cittaro D. Nested Stochastic Block Models applied to single cell data. "
        "BMC Bioinformatics. 2021;22:576."
    ),
    doi="10.1186/s12859-021-04489-7",
    url="https://doi.org/10.1186/s12859-021-04489-7",
    kind="method",
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)


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


def schist_hierarchy_plot_owned(adata: Any, *, key_added: str = "nsbm") -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = _standalone_schist_hierarchy_plot(
        adata,
        key_added=key_added,
        _return_details=True,
    )
    parameters = {"key_added": details["key_added"]}
    warnings = []
    if details["hierarchy_truncated"]:
        warnings.append(
            "The hierarchy panel displays the largest stored blocks within a 500-block rendering bound; "
            "all levels and aggregate metrics remain represented."
        )
    code = schist_hierarchy_plot_code(**parameters)
    plotted = make_plot_result(
        png=png,
        title=f"Schist hierarchy: {details['key_added']}",
        operation="schist_hierarchy_plot",
        parameters=parameters,
        description="Read-only hierarchy, modularity, entropy, and stored block-relationship diagnostics.",
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellSchistHierarchyPlot",
        title=f"Schist hierarchy: {details['key_added']}",
        operation="schist_hierarchy_plot",
        methods=(
            "Validated OpenBio Schist producer provenance, cell memberships, posterior marginal shapes, stored "
            "block-parent mappings, and per-level statistics, then rendered the retained nested hierarchy without "
            "refitting or choosing a preferred level."
        ),
        results=(
            f"Displayed {details['hierarchy_levels']:,} stored hierarchy level(s) with block counts "
            f"{details['cluster_counts_finest_to_root']} from finest to root."
        ),
        key_results=details,
        parameters=parameters,
        references=[SCHIST_REFERENCE, MATPLOTLIB_REFERENCE],
        software_packages=["anndata", "matplotlib", "numpy"],
        warnings=warnings,
        limitations=[
            "The hierarchy is exploratory graph-clustering evidence, not Curated annotation or biological ground truth.",
            "The plot describes one stored fit and does not establish convergence or select an optimal hierarchy level.",
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


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


@register_operation("openbio.node.schisthierarchyplot")
def schist_hierarchy_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Schist Hierarchy Plot")
    require_parameters(parameters, {"key_added"}, operation="Schist Hierarchy Plot")
    plotted, report, code = schist_hierarchy_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


__all__ = [
    "schist_hierarchy_plot",
    "schist_hierarchy_plot_owned",
    "schist_nested_model",
    "schist_nested_model_owned",
]
