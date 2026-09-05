from __future__ import annotations

import math
import time
from collections.abc import Mapping
from typing import Any

from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata, make_plot_result, make_summary_result
from .cell_cycle import analyze_cell_cycle_score, cell_cycle_code
from .cell_cycle_plotting import _standalone_cell_cycle_score_plot, cell_cycle_score_plot_code
from .expression_source import _CELL_CYCLE_SPEC, _MARKER_PLOT_SPEC, DynamicExpressionSource
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
)
from .trajectory_analysis import (
    analyze_diffusion_map,
    analyze_dpt,
    analyze_paga,
    diffusion_map_code,
    dpt_code,
    paga_code,
)
from .trajectory_plotting import (
    diffusion_spectrum_plot_code,
    dpt_gene_trend_plot_code,
    paga_plot_code,
    run_diffusion_spectrum_plot,
    run_dpt_gene_trend_plot,
    run_paga_plot,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

CELL_CYCLE_EXPRESSION_SOURCE = _CELL_CYCLE_SPEC
DPT_GENE_TREND_EXPRESSION_SOURCE = _MARKER_PLOT_SPEC
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"

CELL_CYCLE_REFERENCE = AnalysisReference(
    citation=(
        "Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic "
        "melanoma by single-cell RNA-seq. Science. 2016;352:189-196."
    ),
    doi="10.1126/science.aad0501",
    url="https://doi.org/10.1126/science.aad0501",
    kind="method",
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)
PAGA_REFERENCE = AnalysisReference(
    citation=(
        "Wolf FA, Hamey FK, Plass M, et al. PAGA: graph abstraction reconciles clustering with trajectory "
        "inference through a topology preserving map of single cells. Genome Biology. 2019;20:59."
    ),
    doi="10.1186/s13059-019-1663-x",
    url="https://doi.org/10.1186/s13059-019-1663-x",
    kind="method",
)
DIFFUSION_MAP_REFERENCE = AnalysisReference(
    citation=(
        "Haghverdi L, Buettner F, Theis FJ. Diffusion maps for high-dimensional single-cell analysis of "
        "differentiation data. Bioinformatics. 2015;31:2989-2998."
    ),
    doi="10.1093/bioinformatics/btv325",
    url="https://doi.org/10.1093/bioinformatics/btv325",
    kind="method",
)
DPT_REFERENCE = AnalysisReference(
    citation=(
        "Haghverdi L, Buttner M, Wolf FA, Buettner F, Theis FJ. Diffusion pseudotime robustly reconstructs "
        "lineage branching. Nature Methods. 2016;13:845-848."
    ),
    doi="10.1038/nmeth.3971",
    url="https://doi.org/10.1038/nmeth.3971",
    kind="method",
)


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


def cell_cycle_score_plot_owned(
    adata: Any,
    *,
    output_prefix: str = "cell_cycle",
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = _standalone_cell_cycle_score_plot(
        adata,
        output_prefix=output_prefix,
        _return_details=True,
    )
    parameters = {"output_prefix": details["output_prefix"]}
    warnings = [
        "Cell-cycle scores and phases are supervised per-cell annotations, not a continuous biological clock."
    ]
    code = cell_cycle_score_plot_code(output_prefix=details["output_prefix"])
    title = "Cell-cycle score and phase diagnostics"
    plotted = make_plot_result(
        png=png,
        title=title,
        operation="cell_cycle_score_plot",
        parameters=parameters,
        description=(
            f"Read-only S/G2M score plane and phase counts for {details['plotted_cells']:,} scored cells."
        ),
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCellCycleScorePlot",
        title=title,
        operation="cell_cycle_score_plot",
        methods=(
            "Validated the Cell Cycle Score producer history, observation-axis and score-bundle fingerprints, "
            "canonical S/G2M score columns, categorical phase axis, and Scanpy phase rule. Rendered the stored "
            "score plane and phase counts without reading expression or rescoring genes."
        ),
        results=(
            f"Displayed {details['plotted_cells']:,} cells with stored phase calls "
            + ", ".join(f"{phase}={details['phase_counts'][phase]:,}" for phase in details["phase_order"])
            + "."
        ),
        key_results=details,
        parameters=parameters,
        references=[CELL_CYCLE_REFERENCE, MATPLOTLIB_REFERENCE],
        software_packages=["anndata", "matplotlib", "numpy", "pandas"],
        warnings=warnings,
        limitations=[
            "Phase labels are supervised descriptive annotations, not synchronized time, lineage, or a "
            "proliferation-rate measurement.",
            "Cells are not independent biological replicates; this plot performs no Sample-level Condition inference.",
            "Scores depend on the upstream explicit expression source, observed program genes, control-gene "
            "selection, binning, and random seed.",
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


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


def paga_plot_owned(
    adata: Any,
    *,
    min_connectivity: float = 0.0,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_paga_plot(adata, min_connectivity=min_connectivity)
    parameters = {"min_connectivity": details["min_connectivity"]}
    code = paga_plot_code(min_connectivity=details["min_connectivity"])
    description = (
        f"Rendered {details['plotted_edges']:,} of {details['available_edges']:,} stored undirected PAGA "
        f"edges across {len(details['groups']):,} partition groups."
    )
    warnings = [
        "The circular display layout is deterministic presentation geometry and is not a learned trajectory embedding."
    ]
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="paga_plot",
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellPAGAPlot",
        title="PAGA plot summary",
        operation="paga_plot",
        methods=(
            "Validated the exact OpenBio PAGA producer, named neighbor graph and observation fingerprints, "
            "categorical partition axis, stored connectivity and spanning-forest matrices, group sizes, and output "
            "fingerprint before applying a fixed circular display layout without recomputing PAGA."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key != "title"},
        parameters=parameters,
        references=(PAGA_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("anndata", "numpy", "pandas", "scipy", "scanpy", "matplotlib"),
        warnings=warnings,
        limitations=(
            "PAGA connectivity is exploratory observed-to-expected graph evidence, not a probability or p-value.",
            "Circular node positions carry no distance, direction, root, pseudotime, lineage, fate, or causal meaning.",
            "Edge filtering changes only displayed marks and performs no Sample-level Condition inference.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


def diffusion_spectrum_plot_owned(adata: Any) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_diffusion_spectrum_plot(adata)
    parameters: dict[str, object] = {}
    code = diffusion_spectrum_plot_code()
    description = (
        f"Rendered all {details['plotted_components']:,} stored graph-bound diffusion eigenvalues, including "
        "the separately identified stationary component."
    )
    warnings = [
        "The stationary component is displayed for completeness but is not an informative diffusion coordinate."
    ]
    if not details["source_provenance_available"]:
        warnings.append(
            "OpenBio Diffusion Map provenance was unavailable; graph eigenpairs were validated, "
            "but the original computation settings are unverified."
        )
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="diffusion_spectrum_plot",
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellDiffusionSpectrumPlot",
        title="Diffusion spectrum plot summary",
        operation="diffusion_spectrum_plot",
        methods=(
            "Validated the selected named graph, component axis, non-increasing spectrum, orthonormality, "
            "and graph eigenpair residuals before read-only Matplotlib rendering; checked recorded OpenBio "
            "provenance when available."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key != "title"},
        parameters=parameters,
        references=(DIFFUSION_MAP_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("anndata", "numpy", "scipy", "scanpy", "matplotlib"),
        warnings=warnings,
        limitations=(
            "The spectrum is conditional on the selected named neighbor graph and preprocessing choices.",
            "Diffusion component signs are arbitrary, and an eigenvalue does not define a fate, branch, or measured time.",
            "The plot performs no root selection, pseudotime estimation, Sample-level Condition test, or causal inference.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


def dpt_gene_trend_plot_owned(
    adata: Any,
    *,
    genes: str = "",
    source: DynamicExpressionSource | None = None,
    n_bins: int = 20,
) -> tuple[Any, Any, str]:
    expression = DPT_GENE_TREND_EXPRESSION_SOURCE.resolve(adata, source)
    started_at = time.perf_counter()
    png, details = run_dpt_gene_trend_plot(
        adata,
        genes=genes,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        n_bins=n_bins,
    )
    parameters = {
        "genes": details["genes"],
        **details["source"],
        "n_bins": details["requested_bins"],
    }
    code = dpt_gene_trend_plot_code(
        genes=",".join(details["genes"]),
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        n_bins=details["requested_bins"],
    )
    description = (
        f"Rendered mean and interquartile expression summaries for {len(details['genes']):,} selected gene(s) "
        f"across {details['nonempty_bins']:,} nonempty bins of stored root-dependent diffusion pseudotime."
    )
    warnings = [
        "Displayed bin means and interquartile ranges are descriptive cell-level summaries, not a fitted temporal model.",
        *details["warnings"],
    ]
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="dpt_gene_trend_plot",
        parameters=parameters,
        description=description,
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=details["source_gene_count"],
        started_at=started_at,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellDPTGeneTrendPlot",
        title="DPT gene trend plot summary",
        operation="dpt_gene_trend_plot",
        methods=(
            "Validated the exact OpenBio DPT producer, selected root, observation/pseudotime and named-graph "
            "fingerprints, graph-bound Diffusion Map bundle, explicit expression source and feature axis, selected "
            "gene values, and fixed pseudotime bins before read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key != "title"},
        parameters=parameters,
        references=(DPT_REFERENCE, MATPLOTLIB_REFERENCE),
        software_packages=("anndata", "numpy", "pandas", "scipy", "scanpy", "matplotlib"),
        warnings=warnings,
        limitations=(
            "DPT is a root-dependent relative ordering, not measured time or causal transition direction.",
            "Binned expression summaries are descriptive and do not identify drivers, fates, branches, or regulatory effects.",
            "Cells are not independent biological replicates; no Sample-level Condition inference was performed.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=details["source_gene_count"],
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


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


@register_operation("openbio.node.cellcyclescoreplot")
def cell_cycle_score_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Cell Cycle Score Plot")
    require_parameters(parameters, {"output_prefix"}, operation="Cell Cycle Score Plot")
    plotted, report, code = cell_cycle_score_plot_owned(_input_adata(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


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


@register_operation("openbio.node.pagaplot")
def paga_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="PAGA Plot")
    require_parameters(parameters, {"min_connectivity"}, operation="PAGA Plot")
    plotted, report, code = paga_plot_owned(_input_adata(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.diffusionspectrumplot")
def diffusion_spectrum_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Diffusion Spectrum Plot")
    require_parameters(parameters, set(), operation="Diffusion Spectrum Plot")
    plotted, report, code = diffusion_spectrum_plot_owned(_input_adata(inputs))
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.dptgenetrendplot")
def dpt_gene_trend_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    operation = "DPT Gene Trend Plot"
    require_input_names(inputs, {"adata"}, operation=operation)
    require_parameters(parameters, {"genes", "source", "n_bins"}, operation=operation)
    source = parameters["source"]
    expected = {
        "X": {"source"},
        "raw": {"source"},
        "layer": {"source", "layer_name"},
    }
    if not isinstance(source, dict) or source.get("source") not in expected or set(source) != expected.get(source.get("source")):
        raise ValueError(f"{operation} source must be a closed supported DynamicCombo selection.")
    plotted, report, code = dpt_gene_trend_plot_owned(_input_adata(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


__all__ = [
    "CELL_CYCLE_EXPRESSION_SOURCE",
    "cell_cycle_score",
    "cell_cycle_score_owned",
    "cell_cycle_score_plot",
    "cell_cycle_score_plot_owned",
    "diffusion_map",
    "diffusion_map_owned",
    "diffusion_spectrum_plot",
    "diffusion_spectrum_plot_owned",
    "dpt",
    "dpt_gene_trend_plot",
    "dpt_gene_trend_plot_owned",
    "dpt_owned",
    "paga",
    "paga_plot",
    "paga_plot_owned",
    "paga_owned",
    "parse_dpt_root_value",
]
