from __future__ import annotations

import inspect
import time
from typing import Any

from . import dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import make_plot_result, make_table_result
from .artifact_codecs import read_anndata, read_table, write_plot, write_table
from .artifact_envelope import result_metadata, table_from_metadata
from .contracts import TableResult
from .expression_source import DynamicExpressionSource, ExpressionSource, ExpressionSourceSpec
from .expression_state import resolve_expression_state
from .marker_evidence import (
    MARKER_COLUMNS,
    filter_marker_genes_code,
    filter_marker_table,
    marker_genes_code,
    marker_provenance_parameters,
    rank_marker_evidence,
    validate_marker_artifact_pair,
)
from .operations_input import (
    ANNDATA_CODEC,
    ANNDATA_KIND,
    require_artifact_input,
    require_input_names,
    require_parameters,
)
from .result_plotting import (
    MARKER_PLOT_RNG_LOCK,
    _marker_expression_plot_impl,
    _plot_umap_impl,
    marker_expression_plot_code,
    marker_plot_expression_state,
    umap_embedding_provenance,
    umap_plot_code,
)
from .worker_protocol import JSONValue, OperationContext, register_operation

TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
TABLE_CODEC = "table-jsonl-v1"
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"
PLOT_CODEC = "plot-png-v1"
PCA_METADATA_COLUMNS = [
    "component",
    "metadata",
    "metadata_type",
    "n_cells_total",
    "n_cells_missing",
    "n_samples_total",
    "n_samples_analyzed",
    "n_samples_missing",
    "n_levels",
    "statistic_name",
    "statistic",
    "estimate_name",
    "estimate",
    "effect_size_name",
    "effect_size",
    "p_value",
    "p_adjusted",
    "significant",
]

MARKER_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="Marker source",
    default="layer",
    layer_default="log1p_norm",
)
MARKER_PLOT_EXPRESSION_SOURCE = ExpressionSourceSpec(
    description="Marker plot source",
    default="layer",
    include_raw=True,
    layer_default="log1p_norm",
)

BH_REFERENCE = AnalysisReference(
    citation=(
        "Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful "
        "Approach to Multiple Testing. Journal of the Royal Statistical Society Series B. 1995;57:289-300."
    ),
    url="https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
    doi="10.1111/j.2517-6161.1995.tb02031.x",
    kind="method",
)
SCIPY_REFERENCE = AnalysisReference(
    citation=(
        "Virtanen P et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. "
        "Nature Methods. 2020;17:261-272."
    ),
    url="https://doi.org/10.1038/s41592-019-0686-2",
    doi="10.1038/s41592-019-0686-2",
    kind="software",
)
NUMPY_REFERENCE = AnalysisReference(
    citation="Harris CR et al. Array programming with NumPy. Nature. 2020;585:357-362.",
    url="https://doi.org/10.1038/s41586-020-2649-2",
    doi="10.1038/s41586-020-2649-2",
    kind="software",
)
PANDAS_REFERENCE = AnalysisReference(
    citation=(
        "McKinney W. Data Structures for Statistical Computing in Python. "
        "Proceedings of the 9th Python in Science Conference. 2010:56-61."
    ),
    url="https://doi.org/10.25080/Majora-92bf1922-00a",
    doi="10.25080/Majora-92bf1922-00a",
    kind="software",
)
ANNDATA_REFERENCE = AnalysisReference(
    citation=(
        "Virshup I et al. anndata: Access and store annotated data matrices. "
        "Journal of Open Source Software. 2024;9:4371."
    ),
    url="https://doi.org/10.21105/joss.04371",
    doi="10.21105/joss.04371",
    kind="software",
)
SCANPY_REFERENCE = AnalysisReference(
    citation="Wolf FA, Angerer P, Theis FJ. SCANPY. Genome Biology. 2018;19:15.",
    url="https://doi.org/10.1186/s13059-017-1382-0",
    doi="10.1186/s13059-017-1382-0",
    kind="software",
)
SCANPY_MARKER_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy rank_genes_groups and rank_genes_groups_df API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html",
    kind="software_documentation",
)
SCANPY_FILTER_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy filter_rank_genes_groups API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.filter_rank_genes_groups.html",
    kind="software_documentation",
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    url="https://doi.org/10.1109/MCSE.2007.55",
    doi="10.1109/MCSE.2007.55",
    kind="software",
)
UMAP_REFERENCE = AnalysisReference(
    citation=(
        "McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension "
        "Reduction. arXiv. 2018:1802.03426."
    ),
    url="https://doi.org/10.48550/arXiv.1802.03426",
    doi="10.48550/arXiv.1802.03426",
    kind="method",
)
SCANPY_EMBEDDING_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy embedding plotting API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding.html",
    kind="software_documentation",
)
SEABORN_REFERENCE = AnalysisReference(
    citation="Waskom ML. seaborn: statistical data visualization. Journal of Open Source Software. 2021;6:3021.",
    url="https://doi.org/10.21105/joss.03021",
    doi="10.21105/joss.03021",
    kind="software",
)
SCANPY_DOTPLOT_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy dotplot API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.dotplot.html",
    kind="software_documentation",
)
SCANPY_MATRIXPLOT_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy matrixplot API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.matrixplot.html",
    kind="software_documentation",
)
SCANPY_TRACKSPLOT_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy tracksplot API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.tracksplot.html",
    kind="software_documentation",
)
SCANPY_VIOLIN_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy violin API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.violin.html",
    kind="software_documentation",
)
SCANPY_DENDROGRAM_DOCUMENTATION_REFERENCE = AnalysisReference(
    citation="Scanpy dendrogram API documentation (accessed 2026-08-28).",
    url="https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.tl.dendrogram.html",
    kind="software_documentation",
)
WILCOXON_REFERENCE = AnalysisReference(
    citation="Wilcoxon F. Individual comparisons by ranking methods. Biometrics Bulletin. 1945;1(6):80-83.",
    url="https://doi.org/10.2307/3001968",
    doi="10.2307/3001968",
    kind="method",
)
WELCH_REFERENCE = AnalysisReference(
    citation=(
        "Welch BL. The generalization of Student's problem when several different population variances are "
        "involved. Biometrika. 1947;34:28-35."
    ),
    url="https://doi.org/10.1093/biomet/34.1-2.28",
    doi="10.1093/biomet/34.1-2.28",
    kind="method",
)
SQUAIR_REFERENCE = AnalysisReference(
    citation=(
        "Squair JW et al. Confronting false discoveries in single-cell differential expression. "
        "Nature Communications. 2021;12:5692."
    ),
    url="https://doi.org/10.1038/s41467-021-25960-2",
    doi="10.1038/s41467-021-25960-2",
    kind="practice",
)
STATSMODELS_REFERENCE = AnalysisReference(
    citation=(
        "Seabold S, Perktold J. Statsmodels: Econometric and Statistical Modeling with Python. "
        "Proceedings of the 9th Python in Science Conference. 2010:92-96."
    ),
    url="https://doi.org/10.25080/Majora-92bf1922-011",
    doi="10.25080/Majora-92bf1922-011",
    kind="software",
)
SPEARMAN_REFERENCE = AnalysisReference(
    citation=(
        "Spearman C. The proof and measurement of association between two things. "
        "The American Journal of Psychology. 1904;15:72-101."
    ),
    url="https://doi.org/10.2307/1412159",
    doi="10.2307/1412159",
    kind="method",
)
EFFECT_SIZE_REFERENCE = AnalysisReference(
    citation=(
        "Lakens D. Calculating and reporting effect sizes to facilitate cumulative science: "
        "a practical primer for t-tests and ANOVAs. Frontiers in Psychology. 2013;4:863."
    ),
    url="https://doi.org/10.3389/fpsyg.2013.00863",
    doi="10.3389/fpsyg.2013.00863",
    kind="method",
)
PSEUDOREPLICATION_REFERENCE = AnalysisReference(
    citation=(
        "Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias "
        "in single-cell studies. Nature Communications. 2021;12:738."
    ),
    url="https://doi.org/10.1038/s41467-021-21038-1",
    doi="10.1038/s41467-021-21038-1",
    kind="practice",
)


def _table_record(context: OperationContext, name: str, result: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory(name)
    write_table(root, result.table, result_metadata(result))
    return {
        "type": "artifact",
        "name": name,
        "kind": TABLE_KIND,
        "codec": TABLE_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _plot_record(context: OperationContext, result: Any) -> dict[str, JSONValue]:
    root = context.create_output_directory("plot")
    write_plot(root, result.png, result_metadata(result))
    return {
        "type": "artifact",
        "name": "plot",
        "kind": PLOT_KIND,
        "codec": PLOT_CODEC,
        "payload": root.relative_to(context.output_root).as_posix(),
    }


def _table_input(inputs: dict[str, JSONValue], name: str) -> TableResult:
    root = require_artifact_input(inputs, name, kind=TABLE_KIND, codec=TABLE_CODEC)
    table, metadata = read_table(root)
    return table_from_metadata(metadata, table)


def marker_genes_owned(
    adata: Any,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    source: DynamicExpressionSource | None = None,
    n_genes: int = 100,
    tie_correct: bool = True,
    max_output_rows: int = 1_000_000,
    max_working_memory_gib: float = 4.0,
) -> tuple[Any, Any, Any, str]:
    science = dependencies.require_scientific_dependencies()
    expression = MARKER_EXPRESSION_SOURCE.resolve(adata, source)
    expression_state, expression_evidence = resolve_expression_state(adata, expression)
    started_at = time.perf_counter()
    marker_table, universe_table, details = rank_marker_evidence(
        adata,
        groupby=groupby,
        method=method,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        n_genes=n_genes,
        tie_correct=tie_correct,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        expression_state=expression_state,
        expression_evidence=expression_evidence,
        np=science.np,
        pd=science.pd,
        sparse=science.sparse,
        sc=science.sc,
        ad=science.ad,
    )
    table_parameters = marker_provenance_parameters(
        artifact_role="marker_table",
        groupby=groupby.strip(),
        method=method,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        n_genes=n_genes,
        tie_correct=tie_correct,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        details=details,
    )
    universe_parameters = marker_provenance_parameters(
        artifact_role="tested_gene_universe",
        groupby=groupby.strip(),
        method=method,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        n_genes=n_genes,
        tie_correct=tie_correct,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
        details=details,
    )
    warnings_list = list(details["warnings"])
    table_result = make_table_result(
        title=f"Cluster marker evidence by {groupby.strip()}",
        operation="marker_genes",
        parameters=table_parameters,
        description="Exploratory all-groups-versus-rest Cluster marker evidence.",
        warnings=warnings_list,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        table=marker_table,
    )
    universe_result = make_table_result(
        title="Tested-gene universe for Cluster marker evidence",
        operation="marker_genes",
        parameters=universe_parameters,
        description="Complete ordered gene universe tested by Marker Genes.",
        warnings=warnings_list,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        table=universe_table,
    )
    score_meaning = (
        f"normal-approximation standardized Wilcoxon rank-sum statistic (tie correction {tie_correct!r})"
        if method == "wilcoxon"
        else "Welch-style t statistic"
        if method == "t-test"
        else "Scanpy variance-overestimating t statistic"
    )
    top_markers = {
        group: marker_table.loc[marker_table["group"] == group]
        .head(5)[["gene", "rank", "score", "log2_fold_change_approx", "p_adjusted"]]
        .to_dict(orient="records")
        for group in details["group_labels"]
    }
    code = marker_genes_code(
        groupby=groupby.strip(),
        method=method,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        n_genes=n_genes,
        tie_correct=tie_correct,
        max_output_rows=max_output_rows,
        max_working_memory_gib=max_working_memory_gib,
    )
    method_reference = WILCOXON_REFERENCE if method == "wilcoxon" else WELCH_REFERENCE
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMarkerGenes",
        title=f"Cluster marker evidence by {groupby.strip()}",
        operation="marker_genes",
        methods=(
            f"Scanpy {method!r} ranked every observed {groupby.strip()!r} group against all remaining cells "
            f"using {expression.kind} expression interpreted as {details['expression_interpretation']}. "
            "Globally variable genes were evaluated by Scanpy; globally constant genes were retained as "
            "neutral hypotheses. Positive-score ranking and prevalence were computed, then statsmodels "
            "Benjamini-Hochberg correction was applied separately to the complete tested-gene family for "
            "each group."
        ),
        results=(
            f"Returned {len(marker_table):,} rows of exploratory Cluster marker evidence across "
            f"{len(details['group_labels']):,} groups and preserved the complete ordered universe of "
            f"{details['source_gene_count']:,} tested genes. The score is a {score_meaning}; "
            "log2_fold_change_approx is Scanpy's approximation from mean log expression."
        ),
        key_results={
            "analysis_fingerprint_sha256": details["analysis_fingerprint"],
            "ranking_fingerprint_sha256": details["ranking_fingerprint"],
            "universe_fingerprint_sha256": details["universe_fingerprint"],
            "marker_table_content_fingerprint_sha256": details["table_content_fingerprint"],
            "universe_content_fingerprint_sha256": details["universe_content_fingerprint"],
            "group_sizes": details["group_sizes"],
            "source_gene_count": details["source_gene_count"],
            "requested_n_genes": details["requested_n_genes"],
            "actual_n_genes_per_group": details["actual_n_genes_per_group"],
            "ranking_truncated": details["ranking_truncated"],
            "constant_gene_count": details["constant_gene_count"],
            "variable_gene_count": details["variable_gene_count"],
            "bh_hypotheses_per_group": details["bh_hypotheses_per_group"],
            "marker_output_rows": details["marker_output_rows"],
            "universe_output_rows": details["universe_output_rows"],
            "total_output_rows": details["total_output_rows"],
            "working_memory_preflight": details["working_memory"],
            "score_semantics": score_meaning,
            "numeric_diagnostics": {
                column: summarize_numeric(marker_table[column].to_numpy(dtype=float))
                for column in MARKER_COLUMNS[2:]
            },
            "top_markers_by_group": top_markers,
        },
        parameters=table_parameters,
        references=[
            method_reference,
            BH_REFERENCE,
            SQUAIR_REFERENCE,
            SCANPY_REFERENCE,
            SCANPY_MARKER_DOCUMENTATION_REFERENCE,
            ANNDATA_REFERENCE,
            NUMPY_REFERENCE,
            PANDAS_REFERENCE,
            SCIPY_REFERENCE,
            STATSMODELS_REFERENCE,
        ],
        software_packages=["scanpy", "anndata", "numpy", "pandas", "scipy", "statsmodels"],
        warnings=warnings_list,
        limitations=[
            "Cluster marker evidence is exploratory and depends on clustering resolution, preprocessing, "
            "Technical batch structure, group size, reference composition, and the tested-gene universe.",
            "Cells from the same Sample are not independent replicates; these cell-level p-values must not be "
            "used as a formal Condition contrast or interpreted as cell-type truth.",
            "Benjamini-Hochberg adjustment is per group and does not correct for trying multiple clusterings, "
            "resolutions, rankings, or downstream thresholds.",
            "Prevalence fractions describe values greater than zero in the selected logged representation, "
            "not average expression or proof of exclusivity.",
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return table_result, universe_result, report, code


def filter_marker_genes_owned(
    table: TableResult,
    universe: TableResult,
    min_log2_fold_change: float = 1.0,
    min_fraction_in_group: float = 0.25,
    max_fraction_reference: float = 0.5,
    max_p_adjusted: float = 0.05,
) -> tuple[TableResult, TableResult, Any, str]:
    science = dependencies.require_scientific_dependencies()
    if not isinstance(table, TableResult) or not isinstance(universe, TableResult):
        raise TypeError("Filter Marker Genes requires TableResult values for table and universe.")
    started_at = time.perf_counter()
    frame, diagnostics = filter_marker_table(
        table.table,
        min_log2_fold_change=min_log2_fold_change,
        min_fraction_in_group=min_fraction_in_group,
        max_fraction_reference=max_fraction_reference,
        max_p_adjusted=max_p_adjusted,
        np=science.np,
        pd=science.pd,
    )
    provenance = validate_marker_artifact_pair(
        table,
        universe,
        np=science.np,
        pd=science.pd,
        allowed_table_operations=("marker_genes",),
    )
    upstream = provenance["upstream_parameters"]
    parameters = {
        "marker_evidence_schema_version": upstream["marker_evidence_schema_version"],
        "producer_node_id": "OpenBioSingleCellFilterMarkerGenes",
        "upstream_producer_node_id": upstream["producer_node_id"],
        "artifact_role": "marker_table",
        "analysis_fingerprint": provenance["analysis_fingerprint"],
        "ranking_fingerprint": provenance["ranking_fingerprint"],
        "universe_fingerprint": provenance["universe_fingerprint"],
        "content_fingerprint": diagnostics["content_fingerprint"],
        "upstream_content_fingerprint": provenance["table_content_fingerprint"],
        "marker_groupby": upstream["marker_groupby"],
        "marker_method": provenance["marker_method"],
        "marker_source": upstream["marker_source"],
        "group_labels": provenance["group_labels"],
        "group_sizes": upstream["group_sizes"],
        "source_gene_count": upstream["source_gene_count"],
        "requested_n_genes": upstream["requested_n_genes"],
        "actual_n_genes_per_group": upstream["actual_n_genes_per_group"],
        "ranking_truncated": provenance["ranking_truncated"],
        "tie_correct": upstream["tie_correct"],
        "max_output_rows": upstream["max_output_rows"],
        "max_working_memory_gib": upstream["max_working_memory_gib"],
        "constant_gene_count": upstream["constant_gene_count"],
        "variable_gene_count": upstream["variable_gene_count"],
        "bh_hypotheses_per_group": upstream["bh_hypotheses_per_group"],
        "marker_output_rows": upstream["marker_output_rows"],
        "universe_output_rows": upstream["universe_output_rows"],
        "total_output_rows": upstream["total_output_rows"],
        "groups": "all",
        "reference": "rest",
        "rankby_abs": False,
        "pts": True,
        "corr_method": "benjamini-hochberg",
        "expression_state": upstream["expression_state"],
        "expression_interpretation": upstream["expression_interpretation"],
        "min_log2_fold_change": min_log2_fold_change,
        "min_fraction_in_group": min_fraction_in_group,
        "max_fraction_reference": max_fraction_reference,
        "max_p_adjusted": max_p_adjusted,
        "comparison_semantics": "inclusive",
    }
    for optional_name in ("marker_layer_name", "expression_evidence"):
        if optional_name in upstream:
            parameters[optional_name] = upstream[optional_name]
    warning_list = [
        "Marker thresholds are exploratory analyst choices, not universal biological cutoffs or validation."
    ]
    if not len(frame):
        warning_list.append("No marker rows passed all four validated thresholds.")
    filtered = make_table_result(
        title=f"Filtered {table.title}",
        operation="filter_marker_genes",
        parameters=parameters,
        description="Filtered Cluster marker evidence using four inclusive analyst-selected thresholds.",
        warnings=warning_list,
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        table=frame,
    )
    retained_by_group = {group: int((frame["group"] == group).sum()) for group in provenance["group_labels"]}
    input_by_group = {
        group: int((table.table["group"] == group).sum()) for group in provenance["group_labels"]
    }
    retained_counts = list(retained_by_group.values())
    groups_without_markers = [group for group, count in retained_by_group.items() if count == 0]
    code = filter_marker_genes_code(
        min_log2_fold_change=min_log2_fold_change,
        min_fraction_in_group=min_fraction_in_group,
        max_fraction_reference=max_fraction_reference,
        max_p_adjusted=max_p_adjusted,
        expected_input_content_fingerprint=provenance["table_content_fingerprint"],
        expected_universe_content_fingerprint=provenance["universe_content_fingerprint"],
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellFilterMarkerGenes",
        title=f"Filtered {table.title}",
        operation="filter_marker_genes",
        methods=(
            f"Validated upstream {provenance['marker_method']!r} Cluster marker evidence grouped by "
            f"{upstream['marker_groupby']!r} from {upstream['marker_source']!r} expression, then applied a "
            "deterministic pandas mask using inclusive comparisons: log2_fold_change_approx >= "
            "min_log2_fold_change, fraction_in_group >= min_fraction_in_group, fraction_reference <= "
            "max_fraction_reference, and p_adjusted <= max_p_adjusted. Upstream ranks and row order were "
            "preserved. This differs from Scanpy's strict AnnData-mutating filter helper."
        ),
        results=(
            f"Retained {len(frame):,} of {len(table.table):,} validated Cluster marker-evidence rows "
            f"({len(frame) / len(table.table):.1%} retained)."
            if len(table.table)
            else "Received a valid empty Cluster marker-evidence table and returned a valid empty table."
        ),
        key_results={
            **diagnostics,
            "retained_fraction": len(frame) / len(table.table) if len(table.table) else 0.0,
            "input_rows_by_group": input_by_group,
            "retained_rows_by_group": retained_by_group,
            "groups_without_retained_markers": groups_without_markers,
            "retained_count_min": min(retained_counts, default=0),
            "retained_count_median": float(science.np.median(retained_counts)) if retained_counts else 0.0,
            "retained_count_max": max(retained_counts, default=0),
            "upstream_marker_method": provenance["marker_method"],
            "upstream_groupby": upstream["marker_groupby"],
            "upstream_expression_source": upstream["marker_source"],
            "upstream_group_sizes": upstream["group_sizes"],
            "upstream_requested_n_genes": upstream["requested_n_genes"],
            "upstream_ranking_truncated": provenance["ranking_truncated"],
            "analysis_fingerprint_sha256": provenance["analysis_fingerprint"],
            "ranking_fingerprint_sha256": provenance["ranking_fingerprint"],
            "universe_fingerprint_sha256": provenance["universe_fingerprint"],
            "upstream_marker_content_fingerprint_sha256": provenance["table_content_fingerprint"],
            "filtered_table_content_fingerprint_sha256": diagnostics["content_fingerprint"],
            "validated_numeric_ranges": {
                column: summarize_numeric(table.table[column].to_numpy(dtype=float))
                for column in (
                    "log2_fold_change_approx",
                    "fraction_in_group",
                    "fraction_reference",
                    "p_adjusted",
                )
            },
        },
        parameters=parameters,
        references=[
            BH_REFERENCE,
            SQUAIR_REFERENCE,
            SCANPY_REFERENCE,
            SCANPY_FILTER_DOCUMENTATION_REFERENCE,
            NUMPY_REFERENCE,
            PANDAS_REFERENCE,
            STATSMODELS_REFERENCE,
        ],
        software_packages=["scanpy", "numpy", "pandas", "statsmodels"],
        warnings=warning_list,
        limitations=[
            "Thresholds are analyst choices and repeated tuning can introduce selection bias not represented "
            "by the upstream adjusted p-values.",
            "Filtering cannot recover genes omitted by a truncated upstream ranking.",
            "Passing rows remain exploratory Cluster marker evidence; they do not validate a Curated "
            "annotation, prove marker exclusivity, or support a formal Condition contrast.",
            "Cell-level adjusted p-values do not make cells from the same Sample independent replicates.",
        ],
        input_cells=table.input_cells,
        input_genes=table.input_genes,
        started_at=started_at,
        code=code,
    )
    return filtered, universe, report, code


def umap_plot_owned(
    adata: Any,
    embedding_key: str = "X_umap",
    color: str = "leiden",
    color_mode: str = "auto",
    point_size: float = 10.0,
    continuous_color_map: str = "viridis",
    categorical_palette: str = "tab20",
    sort_order: bool = True,
    missing_color: str = "lightgray",
    legend_policy: str = "automatic",
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = _plot_umap_impl(
        adata,
        embedding_key=embedding_key,
        color=color,
        color_mode=color_mode,
        point_size=point_size,
        continuous_color_map=continuous_color_map,
        categorical_palette=categorical_palette,
        sort_order=sort_order,
        missing_color=missing_color,
        legend_policy=legend_policy,
    )
    embedding_key = details["embedding_key"]
    provenance, provenance_warnings = umap_embedding_provenance(adata, embedding_key)
    warnings = [*details["warnings"], *provenance_warnings]
    parameters = {
        "embedding_key": embedding_key,
        "color": details["color"]["column"] or "",
        "color_mode": color_mode.strip() if isinstance(color_mode, str) else color_mode,
        "point_size": details["rendering"]["point_size_points_squared"],
        "continuous_color_map": continuous_color_map.strip()
        if isinstance(continuous_color_map, str)
        else continuous_color_map,
        "categorical_palette": categorical_palette.strip()
        if isinstance(categorical_palette, str)
        else categorical_palette,
        "sort_order": bool(sort_order),
        "missing_color": missing_color,
        "legend_policy": legend_policy.strip() if isinstance(legend_policy, str) else legend_policy,
    }
    code = umap_plot_code(parameters=parameters, details=details)
    title = f"UMAP colored by {parameters['color']}" if parameters["color"] else "UMAP"
    result = make_plot_result(
        title=title,
        operation="umap_plot",
        parameters=parameters,
        description=(
            f"Read-only rendering of {details['plotted_observations']:,} cells from {embedding_key!r} "
            f"with {details['color']['effective_mode']} observation coloring."
        ),
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        png=png,
    )
    uses_scanpy_metadata = bool(provenance["scanpy_metadata_verified"] or provenance["openbio_history_verified"])
    software_packages = ["matplotlib", "numpy", "pandas", "scipy", "anndata"]
    if uses_scanpy_metadata:
        software_packages.extend(["scanpy", "umap-learn"])
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellUMAPPlot",
        title=title,
        operation="umap_plot",
        methods=(
            f"Validated the stored {embedding_key!r} array as an observation-aligned, finite, real-valued "
            f"embedding with {details['available_coordinate_dimensions']:,} available dimensions, then "
            "rendered its first two dimensions for every observation with Matplotlib Agg. "
            f"Color mode resolved to {details['color']['effective_mode']!r}; point area was fixed at "
            f"{parameters['point_size']:.6g} points squared, alpha at 0.8, and marker linewidth at zero."
        ),
        results=(
            f"Displayed all {details['plotted_observations']:,} input cells. Coordinate ranges were "
            f"{details['coordinate_ranges']!r}; color missingness, ordering, ranges, and legend state are "
            "reported in key_results."
        ),
        key_results={**details, "embedding_provenance": provenance},
        parameters=parameters,
        references=[
            UMAP_REFERENCE,
            SCANPY_EMBEDDING_DOCUMENTATION_REFERENCE,
            MATPLOTLIB_REFERENCE,
            ANNDATA_REFERENCE,
            NUMPY_REFERENCE,
            PANDAS_REFERENCE,
            SCIPY_REFERENCE,
        ],
        software_packages=software_packages,
        warnings=warnings,
        limitations=[
            "UMAP is a nonlinear exploratory embedding; axis orientation, sign, origin, and absolute scale "
            "have no independent biological meaning.",
            "Apparent separation or mixing does not by itself establish a Curated annotation, Technical-batch "
            "correction, preservation of biology, or a Condition effect.",
            "Cell-level coloring is descriptive and does not create a replicate-aware Condition contrast; "
            "the Sample remains the inference unit.",
            "Occlusion, point size, category colors, continuous normalization, and software/font versions can "
            "change visual emphasis without changing the stored coordinates.",
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return result, report, code


def marker_expression_plot_owned(
    adata: Any,
    genes: str = "",
    groupby: str = "leiden",
    plot: dict[str, object] | None = None,
    source: DynamicExpressionSource | None = None,
    group_order: str = "observed",
    random_seed: int = 0,
) -> tuple[Any, Any, str]:
    import matplotlib

    matplotlib.use("Agg", force=True)
    science = dependencies.require_scientific_dependencies()
    expression = MARKER_PLOT_EXPRESSION_SOURCE.resolve(adata, source)
    if expression.kind == "raw":
        current_x_state, current_x_evidence = resolve_expression_state(adata, ExpressionSource("X"))
        expression_state, expression_evidence = "unknown", None
    else:
        expression_state, expression_evidence = resolve_expression_state(adata, expression)
    started_at = time.perf_counter()
    png, details = _marker_expression_plot_impl(
        adata,
        genes=genes,
        groupby=groupby,
        plot=plot,
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        group_order=group_order,
        random_seed=random_seed,
        _science=science,
        _rng_lock=MARKER_PLOT_RNG_LOCK,
    )
    if expression.kind == "raw":
        expression_state, expression_evidence = marker_plot_expression_state(
            adata,
            source_kind="raw",
            genes=details["genes"],
            current_x_state=current_x_state,
            current_x_evidence=current_x_evidence,
            np=science.np,
            sparse=science.sparse,
        )
    plot_type = details["plot"]["plot"]
    parameters = {
        "genes": details["genes"],
        "groupby": details["groupby"],
        "plot": details["plot"],
        **expression.parameters(),
        "group_order": details["group_order_policy"],
        "resolved_group_order": details["resolved_group_order"],
        "random_seed": details["rendering"]["random_seed"],
    }
    warnings = list(details["warnings"])
    if expression_state == "unknown":
        warnings.append(
            "The selected expression representation has no verifiable transformation history. Values were "
            "plotted unchanged without assuming counts, normalization, or logarithmization; confirm the scale."
        )
    elif expression_state == "counts":
        warnings.append(
            "The selected expression representation is proven count-scale data; library size can dominate "
            "cell-level visual differences unless normalization is performed upstream."
        )
    if expression_state == "scaled":
        warnings.append(
            "The selected representation is scaled and may contain negative values; displayed magnitudes are "
            "not expression abundance on the original normalized scale."
        )
    if plot_type == "dotplot":
        semantics = (
            "dot color is the selected mean expression and dot size is the fraction of cells strictly above "
            f"{details['plot']['expression_cutoff']:.6g}"
        )
    elif plot_type == "matrixplot":
        semantics = "color is group mean expression"
    elif plot_type == "tracksplot":
        semantics = "tracks show the cell-level expression distribution within each group"
    else:
        semantics = (
            f"violin width uses density_norm={details['plot']['density_norm']!r} and the y-axis is "
            f"{details['plot']['y_scale']!r}"
        )
    code = marker_expression_plot_code(
        genes=",".join(details["genes"]),
        groupby=details["groupby"],
        plot=details["plot"],
        source_kind=expression.kind,
        layer_name=expression.layer_name,
        group_order=details["group_order_policy"],
        random_seed=details["rendering"]["random_seed"],
        resolved_group_order=details["resolved_group_order"],
    )
    title = f"{plot_type}: {', '.join(details['genes'])} by {details['groupby']}"
    plotted = make_plot_result(
        title=title,
        operation="marker_expression_plot",
        parameters=parameters,
        description=(
            f"Read-only Scanpy {plot_type} rendering of {len(details['genes']):,} selected genes across "
            f"{len(details['resolved_group_order']):,} groups; {semantics}."
        ),
        warnings=warnings,
        input_cells=int(adata.n_obs),
        input_genes=details["source_gene_count"],
        started_at=started_at,
        png=png,
        random_seed=details["rendering"]["random_seed"],
    )
    documentation_reference = {
        "dotplot": SCANPY_DOTPLOT_DOCUMENTATION_REFERENCE,
        "matrixplot": SCANPY_MATRIXPLOT_DOCUMENTATION_REFERENCE,
        "tracksplot": SCANPY_TRACKSPLOT_DOCUMENTATION_REFERENCE,
        "violin": SCANPY_VIOLIN_DOCUMENTATION_REFERENCE,
    }[plot_type]
    references = [
        SCANPY_REFERENCE,
        documentation_reference,
        MATPLOTLIB_REFERENCE,
        ANNDATA_REFERENCE,
        NUMPY_REFERENCE,
        PANDAS_REFERENCE,
        SCIPY_REFERENCE,
    ]
    software_packages = ["scanpy", "matplotlib", "anndata", "numpy", "pandas", "scipy"]
    if plot_type == "violin":
        references.append(SEABORN_REFERENCE)
        software_packages.append("seaborn")
    if details["dendrogram"]["enabled"]:
        references.append(SCANPY_DENDROGRAM_DOCUMENTATION_REFERENCE)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellMarkerExpressionPlot",
        title=title,
        operation="marker_expression_plot",
        methods=(
            f"Selected the exact ordered panel {details['genes']!r} from {expression.kind!r} expression "
            f"({expression_state!r} state), validated finite values and complete categorical groups, and "
            f"rendered a ScanpyV1 {plot_type}. {semantics.capitalize()}. Group order followed "
            f"{details['group_order_policy']!r}; dendrogram order, when requested, used complete linkage on "
            "Pearson correlations of group means computed from this same gene panel and expression source. "
            f"Rendering ran under an isolated NumPy RNG state with random_seed={random_seed}."
        ),
        results=(
            f"Displayed {len(details['genes']):,} genes across {len(details['resolved_group_order']):,} groups "
            f"and {int(adata.n_obs):,} cells as descriptive Cluster marker evidence. The selected panel range "
            f"was {details['selected_panel_range']!r}; exact group counts, order, and group-by-gene descriptive "
            "statistics are available in key_results."
        ),
        key_results={
            **details,
            "expression_state": expression_state,
            "expression_evidence": expression_evidence,
            "display_semantics": semantics,
        },
        parameters=parameters,
        references=references,
        software_packages=software_packages,
        warnings=warnings,
        limitations=[
            "This visualization is descriptive Cluster marker evidence; it is not a replicate-aware Condition "
            "contrast and does not establish a Curated annotation or cell-type truth.",
            "Cells from the same Sample are not independent biological replicates, and cell count is not a "
            "substitute for Sample-level replication.",
            "The selected panel, expression representation, group definition, scaling, cutoff, and density "
            "settings determine the visible pattern and can change qualitative interpretation.",
            "Dendrogram order, when enabled, summarizes only the selected genes and should not be interpreted "
            "as a genome-wide or lineage relationship.",
        ],
        input_cells=int(adata.n_obs),
        input_genes=details["source_gene_count"],
        started_at=started_at,
        code=code,
        random_seed=details["rendering"]["random_seed"],
    )
    return plotted, report, code


def _pca_metadata_associations_impl(
    adata,
    use_rep="X_pca",
    sample_key="sample",
    categorical_obs_keys="",
    continuous_obs_keys="",
    alpha=0.05,
):
    """Compute Sample-level PCA diagnostics without plugin-specific dependencies."""
    import numpy as np
    import pandas as pd
    from scipy import sparse, stats

    columns = [
        "component",
        "metadata",
        "metadata_type",
        "n_cells_total",
        "n_cells_missing",
        "n_samples_total",
        "n_samples_analyzed",
        "n_samples_missing",
        "n_levels",
        "statistic_name",
        "statistic",
        "estimate_name",
        "estimate",
        "effect_size_name",
        "effect_size",
        "p_value",
        "p_adjusted",
        "significant",
    ]

    def parse_keys(value, label):
        if not isinstance(value, str):
            raise TypeError(f"{label} must be a comma-separated string.")
        return list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))

    def render_labels(series, label, *, allow_missing):
        rendered = []
        identities = {}
        missing_count = 0
        for value in series.tolist():
            missing_value = pd.isna(value)
            if not isinstance(missing_value, (bool, np.bool_)):
                raise ValueError(f"{label} values must be scalar.")
            if bool(missing_value):
                rendered.append(pd.NA)
                missing_count += 1
                continue
            text = str(value).strip()
            if not text:
                rendered.append(pd.NA)
                missing_count += 1
                continue
            identity = (type(value).__module__, type(value).__qualname__, repr(value))
            identities.setdefault(text, set()).add(identity)
            rendered.append(text)
        collisions = sorted(text for text, values in identities.items() if len(values) > 1)
        if collisions:
            raise ValueError(f"{label} has distinct values that collapse after string conversion: {collisions!r}.")
        if missing_count and not allow_missing:
            raise ValueError(f"{label} contains {missing_count} missing or blank value(s).")
        return pd.Series(pd.array(rendered, dtype="string"), index=series.index), missing_count

    def is_constant(values):
        return values.size == 0 or bool(np.all(values == values[0]))

    def spearman_associations(x, y):
        """Return SciPy Spearman rho and two-sided p-values for all columns in y."""
        results = [stats.spearmanr(x, y[:, index], alternative="two-sided") for index in range(y.shape[1])]
        statistics = np.asarray([result.statistic for result in results], dtype=float)
        p_values = np.asarray([result.pvalue for result in results], dtype=float)
        return statistics, p_values

    if not isinstance(use_rep, str) or not use_rep.strip():
        raise ValueError("PCA metadata association representation cannot be empty.")
    use_rep = use_rep.strip()
    if use_rep not in adata.obsm:
        raise ValueError(f"PCA metadata association representation not found in obsm: {use_rep!r}")
    if not isinstance(sample_key, str) or not sample_key.strip():
        raise ValueError("PCA metadata association Sample column cannot be empty.")
    sample_key = sample_key.strip()
    if sample_key not in adata.obs:
        raise ValueError(f"PCA metadata association Sample column not found in obs: {sample_key!r}")
    if not adata.obs_names.is_unique:
        raise ValueError("PCA metadata associations require unique observation names.")
    if isinstance(alpha, (bool, np.bool_)):
        raise TypeError("PCA metadata association alpha must be numeric.")
    try:
        alpha = float(alpha)
    except (TypeError, ValueError, OverflowError) as error:
        raise TypeError("PCA metadata association alpha must be numeric.") from error
    if not np.isfinite(alpha) or not 0.0 <= alpha <= 1.0:
        raise ValueError("PCA metadata association alpha must be finite and between zero and one, inclusive.")
    preflight_warnings = []
    if alpha in {0.0, 1.0}:
        preflight_warnings.append(
            f"PCA metadata association alpha={alpha:g} is an executable boundary value; "
            "interpret the significant flag according to the disclosed inclusive threshold."
        )

    categorical_keys = parse_keys(categorical_obs_keys, "categorical_obs_keys")
    continuous_keys = parse_keys(continuous_obs_keys, "continuous_obs_keys")
    if not categorical_keys and not continuous_keys:
        raise ValueError("PCA metadata associations require at least one categorical or continuous obs column.")
    overlap = sorted(set(categorical_keys).intersection(continuous_keys))
    if overlap:
        raise ValueError(f"PCA metadata columns cannot be both categorical and continuous: {overlap!r}")
    selected_keys = [*categorical_keys, *continuous_keys]
    missing_columns = [key for key in selected_keys if key not in adata.obs]
    if missing_columns:
        raise ValueError(f"PCA metadata association columns not found in obs: {missing_columns!r}")
    if sample_key in selected_keys:
        raise ValueError("The Sample column cannot also be tested as PCA metadata.")

    representation = adata.obsm[use_rep]
    if sparse.issparse(representation):
        representation = representation.toarray()
    elif hasattr(representation, "to_numpy"):
        representation = representation.to_numpy()
    if np.iscomplexobj(representation):
        raise ValueError("PCA metadata association representation must be real-valued.")
    try:
        representation = np.asarray(representation, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError("PCA metadata association representation must be numeric.") from error
    if representation.ndim != 2 or representation.shape[0] != int(adata.n_obs):
        raise ValueError("PCA metadata association representation must be two-dimensional and aligned to observations.")
    if representation.shape[0] == 0 or representation.shape[1] == 0:
        raise ValueError("PCA metadata association representation cannot be empty.")
    if not bool(np.isfinite(representation).all()):
        raise ValueError("PCA metadata association representation contains non-finite values.")

    sample_labels, _ = render_labels(adata.obs[sample_key], f"Sample column {sample_key!r}", allow_missing=False)
    components = [f"PC{index}" for index in range(1, representation.shape[1] + 1)]
    score_frame = pd.DataFrame(representation, index=adata.obs_names, columns=components)
    score_frame["_sample"] = sample_labels.to_numpy(dtype=str)
    grouped_scores = score_frame.groupby("_sample", sort=False, observed=True)
    sample_scores = grouped_scores[components].mean()
    cells_per_sample_series = grouped_scores.size()
    n_samples_total = int(sample_scores.shape[0])
    if n_samples_total < 3:
        raise ValueError("PCA metadata associations require at least three independent Samples.")
    cells_per_sample = {str(key): int(value) for key, value in cells_per_sample_series.items()}

    rows = []
    warnings_list = list(preflight_warnings)
    skipped = []
    metadata_diagnostics = {}

    for key in continuous_keys:
        series = adata.obs[key]
        if not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"Continuous PCA metadata {key!r} must have a numeric dtype.")
        if pd.api.types.is_bool_dtype(series.dtype):
            warnings_list.append(
                f"Continuous PCA metadata {key!r} is boolean and is being analyzed numerically as 0/1."
            )
        try:
            values = pd.to_numeric(series, errors="raise").to_numpy(dtype=float, na_value=np.nan)
        except (TypeError, ValueError) as error:
            raise TypeError(f"Continuous PCA metadata {key!r} must contain numeric values.") from error
        if bool(np.isinf(values).any()):
            raise ValueError(f"Continuous PCA metadata {key!r} contains infinite values.")
        n_cells_missing = int(np.isnan(values).sum())
        value_frame = pd.DataFrame({"_sample": sample_labels.to_numpy(dtype=str), "_value": values})
        sample_values = (
            value_frame.groupby("_sample", sort=False, observed=True)["_value"].mean().reindex(sample_scores.index)
        )
        valid = sample_values.notna()
        n_samples_analyzed = int(valid.sum())
        n_samples_missing = n_samples_total - n_samples_analyzed
        metadata_diagnostics[key] = {
            "metadata_type": "continuous",
            "n_cells_total": int(adata.n_obs),
            "n_cells_missing": n_cells_missing,
            "n_samples_total": n_samples_total,
            "n_samples_analyzed": n_samples_analyzed,
            "n_samples_missing": n_samples_missing,
            "test": "scipy.stats.spearmanr",
            "alternative": "two-sided",
        }
        if n_samples_analyzed < 3:
            reason = "Fewer than three Samples remained after missing values were excluded."
            skipped.extend({"component": component, "metadata": key, "reason": reason} for component in components)
            warnings_list.append(f"Continuous metadata {key!r} was skipped: {reason}")
            continue
        x = sample_values.loc[valid].to_numpy(dtype=float)
        if is_constant(x):
            reason = "Sample-level metadata values were constant."
            skipped.extend({"component": component, "metadata": key, "reason": reason} for component in components)
            warnings_list.append(f"Continuous metadata {key!r} was skipped: {reason}")
            continue
        valid_component_indices = []
        for component_index, component in enumerate(components):
            y = sample_scores.loc[valid, component].to_numpy(dtype=float)
            if is_constant(y):
                skipped.append(
                    {"component": component, "metadata": key, "reason": "Sample-mean PC scores were constant."}
                )
            else:
                valid_component_indices.append(component_index)
        statistics = np.asarray([], dtype=float)
        p_values = np.asarray([], dtype=float)
        if valid_component_indices:
            y = sample_scores.loc[valid, [components[index] for index in valid_component_indices]].to_numpy(dtype=float)
            statistics, p_values = spearman_associations(x, y)
        for component_index, statistic, p_value in zip(
            valid_component_indices,
            statistics,
            p_values,
            strict=True,
        ):
            component = components[component_index]
            statistic = float(statistic)
            p_value = float(p_value)
            if not np.isfinite(statistic) or not np.isfinite(p_value):
                skipped.append({"component": component, "metadata": key, "reason": "Spearman result was non-finite."})
                continue
            rows.append(
                {
                    "component": component,
                    "metadata": key,
                    "metadata_type": "continuous",
                    "n_cells_total": int(adata.n_obs),
                    "n_cells_missing": n_cells_missing,
                    "n_samples_total": n_samples_total,
                    "n_samples_analyzed": n_samples_analyzed,
                    "n_samples_missing": n_samples_missing,
                    "n_levels": None,
                    "statistic_name": "spearman_rho",
                    "statistic": statistic,
                    "estimate_name": "spearman_rho",
                    "estimate": statistic,
                    "effect_size_name": "rho_squared",
                    "effect_size": float(statistic**2),
                    "p_value": p_value,
                }
            )

    for key in categorical_keys:
        labels, n_cells_missing = render_labels(adata.obs[key], f"Categorical PCA metadata {key!r}", allow_missing=True)
        sample_values_dict = {}
        for sample in sample_scores.index:
            sample_mask = sample_labels == sample
            observed = labels.loc[sample_mask.to_numpy()].dropna().unique().tolist()
            if len(observed) > 1:
                raise ValueError(
                    f"Categorical PCA metadata {key!r} must be constant within each Sample; "
                    f"Sample {sample!r} contains {observed!r}."
                )
            sample_values_dict[sample] = observed[0] if observed else pd.NA
        sample_values = pd.Series(
            pd.array(list(sample_values_dict.values()), dtype="string"), index=sample_scores.index
        )
        valid = sample_values.notna()
        n_samples_analyzed = int(valid.sum())
        n_samples_missing = n_samples_total - n_samples_analyzed
        levels = sample_values.loc[valid].drop_duplicates().tolist()
        group_sizes = {str(level): int((sample_values.loc[valid] == level).sum()) for level in levels}
        metadata_diagnostics[key] = {
            "metadata_type": "categorical",
            "n_cells_total": int(adata.n_obs),
            "n_cells_missing": int(n_cells_missing),
            "n_samples_total": n_samples_total,
            "n_samples_analyzed": n_samples_analyzed,
            "n_samples_missing": n_samples_missing,
            "n_levels": len(levels),
            "group_sizes": group_sizes,
        }
        if len(levels) < 2:
            reason = "Fewer than two observed Sample-level groups remained."
            skipped.extend({"component": component, "metadata": key, "reason": reason} for component in components)
            warnings_list.append(f"Categorical metadata {key!r} was skipped: {reason}")
            continue
        singleton_groups = {level: size for level, size in group_sizes.items() if size == 1}
        if singleton_groups:
            warnings_list.append(
                f"Categorical metadata {key!r} has singleton Sample group(s) {singleton_groups!r}; "
                "ANOVA estimates are weakly supported."
            )
        below_three = {level: size for level, size in group_sizes.items() if size < 3}
        if below_three and not singleton_groups:
            warnings_list.append(
                f"Categorical metadata {key!r} has fewer than three Samples in group(s) {below_three!r}; "
                "ANOVA estimates are weakly supported."
            )
        if all(size == 1 for size in group_sizes.values()):
            reason = "Every observed category was represented by a single Sample."
            skipped.extend({"component": component, "metadata": key, "reason": reason} for component in components)
            warnings_list.append(f"Categorical metadata {key!r} was skipped: {reason}")
            continue
        for component in components:
            y = sample_scores.loc[valid, component].to_numpy(dtype=float)
            if is_constant(y):
                skipped.append(
                    {"component": component, "metadata": key, "reason": "Sample-mean PC scores were constant."}
                )
                continue
            category = sample_values.loc[valid]
            groups = [y[(category == level).to_numpy()] for level in levels]
            if all(is_constant(group) for group in groups):
                skipped.append(
                    {
                        "component": component,
                        "metadata": key,
                        "reason": "Every categorical group had constant Sample-mean PC scores.",
                    }
                )
                continue
            grand_mean = float(y.mean())
            ss_total = float(np.square(y - grand_mean).sum())
            if not np.isfinite(ss_total) or ss_total <= 0.0:
                skipped.append(
                    {"component": component, "metadata": key, "reason": "Total Sample-level variance was zero."}
                )
                continue
            ss_between = float(sum(group.size * (float(group.mean()) - grand_mean) ** 2 for group in groups))
            result = stats.f_oneway(*groups)
            statistic = float(result.statistic)
            p_value = float(result.pvalue)
            effect_size = float(ss_between / ss_total)
            if not all(np.isfinite(value) for value in (statistic, p_value, effect_size)):
                skipped.append({"component": component, "metadata": key, "reason": "ANOVA result was non-finite."})
                continue
            rows.append(
                {
                    "component": component,
                    "metadata": key,
                    "metadata_type": "categorical",
                    "n_cells_total": int(adata.n_obs),
                    "n_cells_missing": int(n_cells_missing),
                    "n_samples_total": n_samples_total,
                    "n_samples_analyzed": n_samples_analyzed,
                    "n_samples_missing": n_samples_missing,
                    "n_levels": len(levels),
                    "statistic_name": "anova_f",
                    "statistic": statistic,
                    "estimate_name": None,
                    "estimate": None,
                    "effect_size_name": "eta_squared",
                    "effect_size": effect_size,
                    "p_value": p_value,
                }
            )

    if not rows:
        raise ValueError("PCA metadata associations produced no valid hypotheses after Sample-level validation.")
    p_values = np.asarray([row["p_value"] for row in rows], dtype=float)
    adjusted = np.asarray(stats.false_discovery_control(p_values, method="bh"), dtype=float)
    if not bool(np.isfinite(adjusted).all()):
        raise RuntimeError("PCA metadata association adjusted p-values were non-finite.")
    for row, p_adjusted in zip(rows, adjusted, strict=True):
        row["p_adjusted"] = float(p_adjusted)
        row["significant"] = bool(p_adjusted <= alpha)
    table = pd.DataFrame(rows, columns=columns)
    table = table.sort_values(
        ["p_adjusted", "effect_size", "component", "metadata"],
        ascending=[True, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    for optional_column in ("n_levels", "estimate_name", "estimate"):
        table[optional_column] = table[optional_column].astype(object)
        table.loc[table[optional_column].isna(), optional_column] = None
    numeric_columns = [
        "n_cells_total",
        "n_cells_missing",
        "n_samples_total",
        "n_samples_analyzed",
        "n_samples_missing",
        "statistic",
        "effect_size",
        "p_value",
        "p_adjusted",
    ]
    if not all(bool(np.isfinite(table[column].to_numpy(dtype=float)).all()) for column in numeric_columns):
        raise RuntimeError("PCA metadata association table contains non-finite numeric results.")
    estimates = table["estimate"].dropna().to_numpy(dtype=float)
    if estimates.size and not bool(np.isfinite(estimates).all()):
        raise RuntimeError("PCA metadata association table contains non-finite estimates.")
    if skipped:
        warnings_list.append(
            f"Skipped {len(skipped)} degenerate PC-by-metadata hypothesis/hypotheses; reasons are recorded in summary."
        )
    details = {
        "use_rep": use_rep,
        "sample_key": sample_key,
        "categorical_obs_keys": categorical_keys,
        "continuous_obs_keys": continuous_keys,
        "alpha": alpha,
        "components": components,
        "n_samples_total": n_samples_total,
        "cells_per_sample": cells_per_sample,
        "metadata_diagnostics": metadata_diagnostics,
        "skipped_hypotheses": skipped,
        "warnings": warnings_list,
    }
    return table, details


def _pca_metadata_associations_code(
    use_rep: str,
    sample_key: str,
    categorical_obs_keys: list[str],
    continuous_obs_keys: list[str],
    alpha: float,
) -> str:
    implementation = inspect.getsource(_pca_metadata_associations_impl).strip()
    return (
        f"{implementation}\n\n"
        "def pca_metadata_associations(adata):\n"
        "    import warnings\n"
        "    table, details = _pca_metadata_associations_impl(\n"
        "        adata,\n"
        f"        use_rep={use_rep!r},\n"
        f"        sample_key={sample_key!r},\n"
        f"        categorical_obs_keys={','.join(categorical_obs_keys)!r},\n"
        f"        continuous_obs_keys={','.join(continuous_obs_keys)!r},\n"
        f"        alpha={alpha!r},\n"
        "    )\n"
        "    for message in details['warnings']:\n"
        "        warnings.warn(message, UserWarning, stacklevel=2)\n"
        "    return table\n"
    )


def pca_metadata_associations_owned(
    adata: Any,
    use_rep: str = "X_pca",
    sample_key: str = "sample",
    categorical_obs_keys: str = "",
    continuous_obs_keys: str = "",
    alpha: float = 0.05,
) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    table, details = _pca_metadata_associations_impl(
        adata,
        use_rep=use_rep,
        sample_key=sample_key,
        categorical_obs_keys=categorical_obs_keys,
        continuous_obs_keys=continuous_obs_keys,
        alpha=alpha,
    )
    use_rep = details["use_rep"]
    sample_key = details["sample_key"]
    categorical_keys = details["categorical_obs_keys"]
    continuous_keys = details["continuous_obs_keys"]
    alpha = details["alpha"]
    warnings_list = details["warnings"]
    parameters = {
        "use_rep": use_rep,
        "sample_key": sample_key,
        "categorical_obs_keys": categorical_keys,
        "continuous_obs_keys": continuous_keys,
        "alpha": alpha,
        "sample_aggregation": "unweighted_mean",
        "continuous_test": "scipy_spearmanr_two_sided",
        "categorical_test": "one_way_anova",
        "categorical_effect_size": "eta_squared",
        "p_adjust_method": "benjamini_hochberg",
        "p_adjust_scope": "all_valid_pc_metadata_pairs",
    }
    result = make_table_result(
        title=f"Sample-level metadata diagnostics for {use_rep}",
        operation="pca_metadata_associations",
        parameters=parameters,
        description=(
            "Exploratory Sample-level PCA diagnostics using two-sided Spearman tests for continuous metadata "
            "and one-way ANOVA for categorical metadata."
        ),
        warnings=warnings_list,
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        table=table,
    )
    strongest = table.iloc[0]
    significant_count = int(table["significant"].sum())
    requested_hypotheses = len(details["components"]) * (len(categorical_keys) + len(continuous_keys))
    skipped_hypothesis_count = len(details["skipped_hypotheses"])
    top_associations = table.head(10).to_dict(orient="records")
    code = _pca_metadata_associations_code(use_rep, sample_key, categorical_keys, continuous_keys, alpha)
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellPCAMetadataAssociations",
        title=f"Sample-level PCA metadata diagnostics for {use_rep}",
        operation="pca_metadata_associations",
        methods=(
            f"Aggregated {use_rep!r} PC scores to one unweighted mean per Sample from obs[{sample_key!r}]. "
            "Continuous metadata were averaged per Sample and screened with SciPy's two-sided Spearman "
            "rank-correlation test; categorical metadata were required to be constant within each Sample and "
            "screened with standard one-way ANOVA, with eta squared reported as effect size. "
            "Benjamini-Hochberg correction was applied once across every valid PC-by-metadata hypothesis."
        ),
        results=(
            f"Screened {len(table):,} valid Sample-level hypotheses across "
            f"{details['n_samples_total']:,} independent Samples; {significant_count:,} had BH-adjusted "
            f"p <= {alpha:.6g}. The strongest screened pair was {strongest['component']} with "
            f"{strongest['metadata']!r} ({strongest['effect_size_name']}="
            f"{float(strongest['effect_size']):.4g}, adjusted p="
            f"{float(strongest['p_adjusted']):.4g}, analyzed Samples="
            f"{int(strongest['n_samples_analyzed']):,}). These are exploratory diagnostics, not causal "
            "findings or formal Condition inference."
        ),
        key_results={
            "input_cells": int(adata.n_obs),
            "input_genes": int(adata.n_vars),
            "samples": details["n_samples_total"],
            "components": len(details["components"]),
            "cells_per_sample": details["cells_per_sample"],
            "metadata_diagnostics": details["metadata_diagnostics"],
            "requested_hypotheses": requested_hypotheses,
            "valid_hypotheses": len(table),
            "skipped_hypothesis_count": skipped_hypothesis_count,
            "skipped_hypotheses": details["skipped_hypotheses"],
            "significant_hypotheses": significant_count,
            "top_associations": top_associations,
        },
        parameters=parameters,
        references=[
            SPEARMAN_REFERENCE,
            EFFECT_SIZE_REFERENCE,
            BH_REFERENCE,
            PSEUDOREPLICATION_REFERENCE,
            SCIPY_REFERENCE,
            NUMPY_REFERENCE,
            PANDAS_REFERENCE,
            ANNDATA_REFERENCE,
        ],
        software_packages=["anndata", "numpy", "pandas", "scipy"],
        warnings=warnings_list,
        limitations=[
            "This is an exploratory PCA diagnostic. It does not establish causality and is not a formal "
            "Condition contrast or gene-level inferential analysis.",
            "Sample-mean PC positions can reflect cell-type composition, unequal cell recovery, preprocessing, "
            "or confounding; equal Sample weighting prevents cell count from becoming replicate weight but does "
            "not remove those influences.",
            "Standard one-way ANOVA assumes independent Samples, normal within-group populations, and equal "
            "population variances; small groups provide weak support for those assumptions.",
            "Spearman and ANOVA p-values are limited by the number of independent Samples, and correlated PCs or "
            "metadata make the screened hypotheses dependent even after Benjamini-Hochberg adjustment.",
            "SciPy's standard Spearman p-value is an asymptotic calculation whose accuracy is limited for very "
            "small Sample counts; this node is an exploratory screen and does not perform a permutation test.",
        ],
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return result, report, code


@register_operation("openbio.node.markergenes")
def marker_genes(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Marker Genes")
    require_parameters(
        parameters,
        {
            "groupby",
            "method",
            "source",
            "n_genes",
            "tie_correct",
            "max_output_rows",
            "max_working_memory_gib",
        },
        operation="Marker Genes",
    )
    root = require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    table, universe, report, code = marker_genes_owned(read_anndata(root), **parameters)
    return [
        _table_record(context, "table", table),
        _table_record(context, "universe", universe),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.filtermarkergenes")
def filter_marker_genes(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"table", "universe"}, operation="Filter Marker Genes")
    require_parameters(
        parameters,
        {
            "min_log2_fold_change",
            "min_fraction_in_group",
            "max_fraction_reference",
            "max_p_adjusted",
        },
        operation="Filter Marker Genes",
    )
    table, _universe, report, code = filter_marker_genes_owned(
        _table_input(inputs, "table"),
        _table_input(inputs, "universe"),
        **parameters,
    )
    return [
        _table_record(context, "table", table),
        {"type": "input_ref", "name": "universe", "input": "universe"},
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.umapplot")
def umap_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="UMAP Plot")
    require_parameters(
        parameters,
        {
            "embedding_key",
            "color",
            "color_mode",
            "point_size",
            "continuous_color_map",
            "categorical_palette",
            "sort_order",
            "missing_color",
            "legend_policy",
        },
        operation="UMAP Plot",
    )
    root = require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    plotted, report, code = umap_plot_owned(read_anndata(root), **parameters)
    return [
        _plot_record(context, plotted),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.markerexpressionplot")
def marker_expression_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="Marker Expression Plot")
    require_parameters(
        parameters,
        {"genes", "groupby", "plot", "source", "group_order", "random_seed"},
        operation="Marker Expression Plot",
    )
    root = require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    plotted, report, code = marker_expression_plot_owned(read_anndata(root), **parameters)
    return [
        _plot_record(context, plotted),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


@register_operation("openbio.node.pcametadataassociations")
def pca_metadata_associations(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="PCA Metadata Associations")
    require_parameters(
        parameters,
        {"use_rep", "sample_key", "categorical_obs_keys", "continuous_obs_keys", "alpha"},
        operation="PCA Metadata Associations",
    )
    root = require_artifact_input(inputs, "adata", kind=ANNDATA_KIND, codec=ANNDATA_CODEC)
    table, report, code = pca_metadata_associations_owned(read_anndata(root), **parameters)
    return [
        _table_record(context, "table", table),
        {"type": "summary", "name": "summary", "value": result_metadata(report)},
        {"type": "string", "name": "code", "value": code},
    ]


__all__ = [
    "MARKER_EXPRESSION_SOURCE",
    "MARKER_PLOT_EXPRESSION_SOURCE",
    "PCA_METADATA_COLUMNS",
    "filter_marker_genes",
    "filter_marker_genes_owned",
    "marker_genes",
    "marker_genes_owned",
    "marker_expression_plot",
    "marker_expression_plot_owned",
    "pca_metadata_associations",
    "pca_metadata_associations_owned",
    "umap_plot",
    "umap_plot_owned",
]
