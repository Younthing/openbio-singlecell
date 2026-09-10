from __future__ import annotations

import json
import time
from typing import Any

from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import make_plot_result, make_table_result
from .monocle2 import (
    MONOCLE2_CODEC,
    MONOCLE2_KIND,
    MONOCLE2_PACKAGES,
    reproduction_code,
    run_analysis,
    runtime_reproduction_code,
)
from .operations_input import (
    analysis_outputs,
    read_anndata_input,
    read_table_input,
    require_artifact_input,
    require_input_names,
    require_parameters,
    write_anndata_output,
    write_plot_output,
    write_table_output,
)
from .r_runtime import probe_r_runtime
from .worker_protocol import (
    JSONValue,
    OperationContext,
    ProtocolError,
    _object_without_duplicates,
    _reject_constant,
    register_operation,
)

MONOCLE2_REFERENCE = AnalysisReference(
    citation="Monocle 2 documentation: constructing and interpreting single-cell trajectories with DDRTree.",
    url="https://cole-trapnell-lab.github.io/monocle-release/docs/",
    kind="software_documentation",
)


@register_operation("openbio.node.monocle2runtime")
def monocle2_runtime(context, inputs, parameters):
    require_input_names(inputs, set(), operation="Monocle 2 Runtime")
    require_parameters(parameters, {"rscript", "r_home", "library_paths", "path_prefix"},
                       operation="Monocle 2 Runtime")
    started_at = time.perf_counter()
    runtime = probe_r_runtime(**parameters, packages=MONOCLE2_PACKAGES)
    versions = {name: package["version"] for name, package in runtime["packages"].items()}
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMonocle2Runtime", title="Monocle 2 Runtime", operation="monocle2_runtime",
        methods="Probed the selected Rscript and actual R package identities in its child-process environment.",
        results=f"R {runtime['r_version']} loaded {len(versions)} package identities.",
        key_results={"r_version": runtime["r_version"], "packages": versions}, parameters=parameters,
        references=[MONOCLE2_REFERENCE], software_packages=[], warnings=[], limitations=[],
        input_cells=0, input_genes=0, started_at=started_at, code=runtime_reproduction_code(parameters),
    )
    report.summary["software_versions"].update({f"R:{name}": version for name, version in versions.items()})
    report.summary["software_versions"]["R"] = runtime["r_version"]
    return analysis_outputs(report, code, {
        "type": "string", "name": "r_runtime", "value": json.dumps(runtime, ensure_ascii=False, allow_nan=False),
    })


def _json_object(value: JSONValue, name: str) -> dict[str, Any]:
    try:
        parsed = json.loads(value, parse_constant=_reject_constant, object_pairs_hook=_object_without_duplicates)
    except (TypeError, json.JSONDecodeError) as error:
        raise ProtocolError(f"{name} must contain a JSON object.") from error
    if not isinstance(parsed, dict):
        raise ProtocolError(f"{name} must contain a JSON object.")
    return parsed


def _genes(value: JSONValue) -> list[str]:
    try:
        genes = json.loads(value)
    except (TypeError, json.JSONDecodeError) as error:
        raise ProtocolError("genes_json must be a JSON array of feature IDs.") from error
    if not isinstance(genes, list) or any(not isinstance(gene, str) for gene in genes):
        raise ProtocolError("genes_json must be a JSON array of feature IDs.")
    return genes


def _report(operation, title, node_id, native_parameters, runtime, result, started_at):
    details = dict(result["summary"])
    warnings = list(details["warnings"])
    if operation == "ddrtree":
        details["backend_fixed_seed"] = 2016
        warnings.append("Monocle 2 reduceDimension owns random seeding and sets its seed to 2016 internally.")
    results = f"Processed {details['n_obs']:,} cells and {details['n_vars']:,} features with native Monocle 2."
    if operation == "order_cells" and native_parameters["root_state"] is None:
        results += " Pseudotime remains biologically unoriented; the native graph endpoint supplies its direction."
    limitations = []
    if operation in {"differential_test", "beam"}:
        results += f" Retained all {len(result['table']):,} native gene-test rows, including failed fits."
        limitations.append("Cell-level trajectory associations do not establish replicate-aware Condition effects.")
    report, code = make_analysis_report(
        node_id=node_id,
        title=title,
        operation=f"monocle2_{operation}",
        methods=f"Executed native Monocle 2 {operation} through Rscript using the user's expression and parameters.",
        results=results,
        key_results=details,
        parameters=native_parameters,
        references=[MONOCLE2_REFERENCE],
        software_packages=["anndata", "numpy", "pandas", "scipy"],
        warnings=warnings,
        limitations=limitations,
        input_cells=details["n_obs"],
        input_genes=details["n_vars"],
        started_at=started_at,
        code=reproduction_code(operation, native_parameters, runtime),
    )
    report.summary["software_versions"].update(
        {f"R:{package}": version for package, version in details["versions"].items()}
    )
    return report, code


def _execute_cds(context, inputs, parameters, native_parameters, *, operation, title, suffix, outputs):
    started_at = time.perf_counter()
    runtime = _json_object(parameters["r_runtime"], "r_runtime")
    cds_path = require_artifact_input(inputs, "cds", kind=MONOCLE2_KIND, codec=MONOCLE2_CODEC)
    output_dir = context.create_output_directory("cds" if "cds" in outputs else "native")
    result = run_analysis(operation, runtime, output_dir, native_parameters, cds_path=cds_path,
                          adata=read_anndata_input(inputs) if "adata" in inputs else None)
    report, code = _report(operation, title, f"OpenBioSingleCellMonocle2{suffix}",
                           native_parameters, runtime, result, started_at)
    records = []
    if "cds" in outputs:
        records.append({"type": "artifact", "name": "cds", "kind": MONOCLE2_KIND, "codec": MONOCLE2_CODEC,
                        "payload": output_dir.relative_to(context.output_root).as_posix()})
    if "adata" in outputs:
        records.append(write_anndata_output(context, result["adata"]))
    if "table" in outputs:
        table = make_table_result(table=result["table"].convert_dtypes(), title=title, operation=f"monocle2_{operation}",
            parameters=native_parameters, description=report.description, warnings=report.warnings,
            input_cells=report.input_cells, input_genes=report.input_genes, started_at=started_at)
        records.append(write_table_output(context, table, kind="OPENBIO_SINGLE_CELL_TABLE"))
    if "plot" in outputs:
        plot = make_plot_result(png=result["plot_path"].read_bytes(), title=title, operation=f"monocle2_{operation}",
            parameters=native_parameters, description=report.description, warnings=report.warnings,
            input_cells=report.input_cells, input_genes=report.input_genes, started_at=started_at)
        records.append(write_plot_output(context, plot, kind="OPENBIO_SINGLE_CELL_PLOT"))
    return analysis_outputs(report, code, *records)


@register_operation("openbio.node.monocle2prepare")
def monocle2_prepare(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Monocle 2 Prepare")
    require_parameters(parameters, {
        "r_runtime", "source", "gene_short_name_column", "expression_family", "lower_detection_limit",
        "estimate_size_factors", "estimate_dispersions", "detect_genes", "size_factor_parameters_json",
        "dispersion_parameters_json", "detect_parameters_json",
    }, operation="Monocle 2 Prepare")
    started_at = time.perf_counter()
    runtime = _json_object(parameters["r_runtime"], "r_runtime")
    native_parameters = {
        "source": parameters["source"],
        "gene_short_name_column": parameters["gene_short_name_column"],
        "expression_family": parameters["expression_family"],
        "lowerDetectionLimit": parameters["lower_detection_limit"],
        "estimate_size_factors": parameters["estimate_size_factors"],
        "estimate_dispersions": parameters["estimate_dispersions"],
        "detect_genes": parameters["detect_genes"],
        "size_factor_parameters": _json_object(parameters["size_factor_parameters_json"], "size_factor_parameters_json"),
        "dispersion_parameters": _json_object(parameters["dispersion_parameters_json"], "dispersion_parameters_json"),
        "detect_parameters": _json_object(parameters["detect_parameters_json"], "detect_parameters_json"),
    }
    output_dir = context.create_output_directory("cds")
    result = run_analysis("create", runtime, output_dir, native_parameters, adata=read_anndata_input(inputs))
    report, code = _report("create", "Monocle 2 Prepare", "OpenBioSingleCellMonocle2Prepare",
                           native_parameters, runtime, result, started_at)
    return analysis_outputs(report, code, {
        "type": "artifact", "name": "cds", "kind": MONOCLE2_KIND, "codec": MONOCLE2_CODEC,
        "payload": output_dir.relative_to(context.output_root).as_posix(),
    })


@register_operation("openbio.node.monocle2orderinggenes")
def monocle2_ordering_genes(context, inputs, parameters):
    require_input_names(inputs, {"cds"}, operation="Monocle 2 Ordering Genes")
    require_parameters(parameters, {
        "r_runtime", "method", "genes_json", "var_column", "min_mean_expression", "dispersion_fold",
    }, operation="Monocle 2 Ordering Genes")
    return _execute_cds(context, inputs, parameters, {
        "method": parameters["method"], "genes": _genes(parameters["genes_json"]),
        "var_column": parameters["var_column"], "mean_expression": parameters["min_mean_expression"],
        "dispersion_fold": parameters["dispersion_fold"],
    }, operation="ordering_genes", title="Monocle 2 Ordering Genes", suffix="OrderingGenes", outputs=("cds", "table"))


@register_operation("openbio.node.monocle2ddrtree")
def monocle2_ddrtree(context, inputs, parameters):
    require_input_names(inputs, {"cds"}, operation="Monocle 2 DDRTree")
    require_parameters(parameters, {
        "r_runtime", "num_components", "norm_method", "residual_formula", "pseudo_expr",
        "auto_param_selection", "scaling", "extra_parameters_json",
    }, operation="Monocle 2 DDRTree")
    return _execute_cds(context, inputs, parameters, {
        "reduction_method": "DDRTree", "max_components": parameters["num_components"],
        "norm_method": parameters["norm_method"], "residualModelFormulaStr": parameters["residual_formula"] or None,
        "pseudo_expr": parameters["pseudo_expr"],
        "auto_param_selection": parameters["auto_param_selection"], "scaling": parameters["scaling"],
        "extra_parameters": _json_object(parameters["extra_parameters_json"], "extra_parameters_json"),
    }, operation="ddrtree", title="Monocle 2 DDRTree", suffix="DDRTree", outputs=("cds",))


@register_operation("openbio.node.monocle2ordercells")
def monocle2_order_cells(context, inputs, parameters):
    require_input_names(inputs, {"cds"}, operation="Monocle 2 Order Cells")
    require_parameters(parameters, {"r_runtime", "root", "reverse"}, operation="Monocle 2 Order Cells")
    root = parameters["root"]
    if not isinstance(root, dict) or root.get("root") not in {"automatic", "state"}:
        raise ProtocolError("Monocle 2 root must select automatic or state.")
    if root["root"] == "state":
        require_parameters(root, {"root", "state"}, operation="Monocle 2 root selection")
        root_state = root["state"]
        if root_state is None or isinstance(root_state, str) and not root_state.strip():
            raise ValueError("Choose a root State after inspecting the initial trajectory, or select automatic orientation.")
    else:
        require_parameters(root, {"root"}, operation="Monocle 2 root selection")
        root_state = None
    return _execute_cds(context, inputs, parameters, {
        "root_state": root_state, "reverse": parameters["reverse"],
    }, operation="order_cells", title="Monocle 2 Order Cells", suffix="OrderCells", outputs=("cds", "table"))


@register_operation("openbio.node.monocle2trajectoryplot")
def monocle2_trajectory_plot(context, inputs, parameters):
    require_input_names(inputs, {"cds"}, operation="Monocle 2 Trajectory Plot")
    require_parameters(parameters, {
        "r_runtime", "color_by", "x_dimension", "y_dimension", "show_tree", "show_branch_points", "show_state_number",
        "cell_size", "width", "height", "dpi", "extra_parameters_json",
    }, operation="Monocle 2 Trajectory Plot")
    native = {name: parameters[name] for name in (
        "color_by", "show_tree", "show_branch_points", "show_state_number", "cell_size", "width", "height", "dpi",
    )}
    native.update(x=parameters["x_dimension"], y=parameters["y_dimension"],
                  extra_parameters=_json_object(parameters["extra_parameters_json"], "extra_parameters_json"))
    return _execute_cds(context, inputs, parameters, native, operation="trajectory_plot", title="Monocle 2 Trajectory Plot",
                        suffix="TrajectoryPlot", outputs=("plot",))


@register_operation("openbio.node.monocle2differentialtest")
def monocle2_differential_test(context, inputs, parameters):
    require_input_names(inputs, {"cds"}, operation="Monocle 2 Differential Test")
    require_parameters(parameters, {
        "r_runtime", "genes_json", "full_formula", "reduced_formula", "relative_expr", "cores", "extra_parameters_json",
    }, operation="Monocle 2 Differential Test")
    return _execute_cds(context, inputs, parameters, {
        "genes": _genes(parameters["genes_json"]), "fullModelFormulaStr": parameters["full_formula"],
        "reducedModelFormulaStr": parameters["reduced_formula"], "relative_expr": parameters["relative_expr"],
        "cores": parameters["cores"],
        "extra_parameters": _json_object(parameters["extra_parameters_json"], "extra_parameters_json"),
    }, operation="differential_test", title="Monocle 2 Differential Test", suffix="DifferentialTest", outputs=("table",))


@register_operation("openbio.node.monocle2beam")
def monocle2_beam(context, inputs, parameters):
    require_input_names(inputs, {"cds"}, operation="Monocle 2 BEAM")
    require_parameters(parameters, {
        "r_runtime", "genes_json", "branch_point", "progenitor_method", "cores", "extra_parameters_json",
    }, operation="Monocle 2 BEAM")
    return _execute_cds(context, inputs, parameters, {
        "genes": _genes(parameters["genes_json"]), "branch_point": parameters["branch_point"],
        "progenitor_method": parameters["progenitor_method"], "cores": parameters["cores"],
        "extra_parameters": _json_object(parameters["extra_parameters_json"], "extra_parameters_json"),
    }, operation="beam", title="Monocle 2 BEAM", suffix="BEAM", outputs=("table",))


@register_operation("openbio.node.monocle2genetrends")
def monocle2_gene_trends(context, inputs, parameters):
    require_input_names(inputs, {"cds", "table"} if "table" in inputs else {"cds"}, operation="Monocle 2 Gene Trends")
    require_parameters(parameters, {
        "r_runtime", "genes_json", "top_n", "color_by", "trend_formula", "relative_expr", "ncol",
        "width", "height", "dpi", "extra_parameters_json",
    }, operation="Monocle 2 Gene Trends")
    genes = _genes(parameters["genes_json"])
    if not genes and "table" in inputs:
        table, _metadata = read_table_input(inputs, "table", kind="OPENBIO_SINGLE_CELL_TABLE")
        columns = [column for column in ("qval", "pval") if column in table.columns]
        if columns:
            table = table.sort_values(columns, kind="stable")
        genes = table.head(parameters["top_n"])["gene_id"].tolist()
    native = {name: parameters[name] for name in (
        "color_by", "trend_formula", "relative_expr", "ncol", "width", "height", "dpi",
    )}
    native.update(genes=genes, extra_parameters=_json_object(parameters["extra_parameters_json"], "extra_parameters_json"))
    return _execute_cds(context, inputs, parameters, native, operation="gene_trends", title="Monocle 2 Gene Trends",
                        suffix="GeneTrends", outputs=("plot",))


@register_operation("openbio.node.monocle2export")
def monocle2_export(context, inputs, parameters):
    require_input_names(inputs, {"adata", "cds"}, operation="Monocle 2 Export")
    require_parameters(parameters, {
        "r_runtime", "embedding_key", "pseudotime_key", "state_key", "overwrite_existing",
    }, operation="Monocle 2 Export")
    native = {name: parameters[name] for name in ("embedding_key", "pseudotime_key", "state_key", "overwrite_existing")}
    return _execute_cds(context, inputs, parameters, native, operation="export", title="Monocle 2 Export",
                        suffix="Export", outputs=("adata", "table"))
