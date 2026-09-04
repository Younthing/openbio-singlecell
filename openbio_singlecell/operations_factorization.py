from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict
from typing import TYPE_CHECKING, Any

from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, make_plot_result, make_table_result, matrix_totals_and_nonzero
from .artifact_envelope import table_from_metadata
from .cnmf_native_codec import CNMF_NATIVE_CODEC, checkout_cnmf_run, write_cnmf_run
from .cnmf_run import CNMFRun
from .cnmf_standalone import (
    CNMF_AUDITED_PYPI_WHEEL_SHA256,
    CNMF_LOCAL_NEIGHBORHOOD_SIZE,
    CNMF_REQUIRED_VERSION,
    close_run_preserving_primary,
    equivalent_source,
)
from .cnmf_standalone import (
    cnmf_consensus_programs as _cnmf_consensus_programs_science,
)
from .cnmf_standalone import (
    cnmf_rank_survey as _cnmf_rank_survey_science,
)
from .expression_source import _CNMF_SPEC, DynamicExpressionSource, ExpressionSource
from .factorization_plotting import (
    cnmf_programs_plot_code,
    cnmf_rank_plot_code,
    run_cnmf_programs_plot,
    run_cnmf_rank_plot,
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
from .worker_protocol import JSONValue, OperationContext, register_operation

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/factorization"
CNMF_SOFTWARE_PACKAGES = (
    "cnmf",
    "anndata",
    "scanpy",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "pyyaml",
)

CNMF_REFERENCE = AnalysisReference(
    citation=(
        "Kotliar D, Veres A, Nagy MA, et al. Identifying gene expression programs of cell-type "
        "identity and cellular activity with single-cell RNA-Seq. eLife. 2019;8:e43803."
    ),
    doi="10.7554/eLife.43803",
    url="https://doi.org/10.7554/eLife.43803",
    kind="method",
)
CNMF_REPOSITORY_REFERENCE = AnalysisReference(
    citation="Kotliar et al. cNMF method-author repository and implementation.",
    url="https://github.com/dylkot/cNMF",
    kind="software_documentation",
)
CNMF_GUIDE_REFERENCE = AnalysisReference(
    citation="Kotliar et al. cNMF stepwise guide.",
    url="https://github.com/dylkot/cNMF/blob/master/Stepwise_Guide.md",
    kind="software_documentation",
)
CNMF_PYPI_REFERENCE = AnalysisReference(
    citation=(
        "cNMF 1.7.1 Python distribution (MIT); audited PyPI wheel SHA-256 "
        f"{CNMF_AUDITED_PYPI_WHEEL_SHA256}. The installed artifact is not attested by this runtime check."
    ),
    url="https://pypi.org/project/cnmf/1.7.1/",
    kind="software_documentation",
)
ANNDATA_REFERENCE = AnalysisReference(
    citation=(
        "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated "
        "data matrices. Journal of Open Source Software. 2024;9(101):4371."
    ),
    doi="10.21105/joss.04371",
    url="https://doi.org/10.21105/joss.04371",
    kind="software",
)
NUMPY_REFERENCE = AnalysisReference(
    citation="Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature. 2020;585:357-362.",
    doi="10.1038/s41586-020-2649-2",
    url="https://doi.org/10.1038/s41586-020-2649-2",
    kind="software",
)
PANDAS_REFERENCE = AnalysisReference(
    citation="The pandas development team. pandas-dev/pandas: Pandas. Zenodo.",
    doi="10.5281/zenodo.3509134",
    url="https://doi.org/10.5281/zenodo.3509134",
    kind="software",
)
SCIPY_REFERENCE = AnalysisReference(
    citation="Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0. Nature Methods. 2020;17:261-272.",
    doi="10.1038/s41592-019-0686-2",
    url="https://doi.org/10.1038/s41592-019-0686-2",
    kind="software",
)
SCIKIT_LEARN_REFERENCE = AnalysisReference(
    citation="Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn. JMLR. 2011;12:2825-2830.",
    url="https://jmlr.org/papers/v12/pedregosa11a.html",
    kind="software",
)
CNMF_REFERENCES = (
    CNMF_REFERENCE,
    CNMF_REPOSITORY_REFERENCE,
    CNMF_GUIDE_REFERENCE,
    CNMF_PYPI_REFERENCE,
    ANNDATA_REFERENCE,
    NUMPY_REFERENCE,
    PANDAS_REFERENCE,
    SCIPY_REFERENCE,
    SCIKIT_LEARN_REFERENCE,
)
MATPLOTLIB_REFERENCE = AnalysisReference(
    citation="Hunter JD. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007;9:90-95.",
    doi="10.1109/MCSE.2007.55",
    url="https://doi.org/10.1109/MCSE.2007.55",
    kind="software",
)


def _source_string(expression: ExpressionSource) -> str:
    if expression.kind == "X":
        return "X"
    if expression.kind == "raw":
        return "raw"
    if expression.kind == "layer":
        return f"layer:{expression.layer_name}"
    raise ValueError(f"Unsupported cNMF expression source: {expression.kind!r}.")


def _source_label(expression: ExpressionSource) -> str:
    return f"layer:{expression.layer_name}" if expression.kind == "layer" else expression.kind


def _rank_metrics_fingerprint(table: Any) -> str:
    rows = [
        [
            int(row.k),
            float(row.stability),
            float(row.prediction_error),
            float(row.statistics_density_threshold),
            int(row.expected_restarts),
            int(row.completed_restarts),
            int(row.combined_components),
        ]
        for row in table.itertuples(index=False)
    ]
    payload = json.dumps(rows, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return f"sha256:{hashlib.sha256(payload).hexdigest()}"


def _survey_code(
    expression: ExpressionSource,
    *,
    components_min: int,
    components_max: int,
    n_iter: int,
    num_highvar_genes: int,
    random_seed: int,
) -> str:
    source = equivalent_source().rstrip()
    return (
        f"{source}\n\n"
        "_cnmf_rank_survey_equivalent = cnmf_rank_survey\n\n"
        "def cnmf_rank_survey(adata):\n"
        "    return _cnmf_rank_survey_equivalent(\n"
        "        adata,\n"
        f"        source={_source_string(expression)!r},\n"
        f"        components_min={components_min},\n"
        f"        components_max={components_max},\n"
        f"        n_iter={n_iter},\n"
        f"        num_highvar_genes={num_highvar_genes},\n"
        f"        random_seed={random_seed},\n"
        "    )\n"
    )


def _consensus_code(
    *,
    selected_k: int,
    density_threshold: float,
    local_neighborhood_size: float,
    n_top_genes: int,
    overwrite_existing: bool,
) -> str:
    source = equivalent_source().rstrip()
    return (
        f"{source}\n\n"
        "_cnmf_consensus_programs_equivalent = cnmf_consensus_programs\n\n"
        "def cnmf_consensus_programs(run):\n"
        "    return _cnmf_consensus_programs_equivalent(\n"
        "        run,\n"
        f"        selected_k={selected_k},\n"
        f"        density_threshold={density_threshold!r},\n"
        f"        local_neighborhood_size={local_neighborhood_size!r},\n"
        f"        n_top_genes={n_top_genes},\n"
        f"        overwrite_existing={overwrite_existing!r},\n"
        "    )\n"
    )


def _warning_text(records: tuple[tuple[str, int], ...]) -> list[str]:
    messages: list[str] = []
    for label, count in records:
        if label == "cnmf-1.7.1-unclosed-yaml-reader":
            messages.append(
                f"cnmf==1.7.1 emitted its known unclosed YAML-reader ResourceWarning {count} time(s); "
                "only the exact audited warning was isolated."
            )
        elif label == "cnmf-1.7.1-unclosed-hvg-reader":
            messages.append(
                f"cnmf==1.7.1 emitted its known unclosed high-variance-gene reader ResourceWarning {count} "
                "time(s); only the exact audited warning was isolated."
            )
        elif label == "cnmf-1.7.1-anndata-view-materialization":
            messages.append(
                f"cnmf==1.7.1 emitted its known AnnData view-materialization warning {count} time(s) during "
                "final usage refitting; only the exact audited warning was isolated."
            )
        elif label == "cnmf-1.7.1-pandas4-sum-positional":
            messages.append(
                f"cnmf==1.7.1 emitted its known pandas-4 positional-sum compatibility warning {count} time(s); "
                "only the exact audited warning was isolated."
            )
    return messages


class OpenBioSingleCellCNMFRankSurvey:
    EXPRESSION_SOURCE = _CNMF_SPEC

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        components_min: int = 3,
        components_max: int = 19,
        n_iter: int = 100,
        num_highvar_genes: int = 2000,
        random_seed: int = 123,
    ) -> tuple[CNMFRun, Any, Any, str]:
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        started_at = time.perf_counter()
        run, metric_table = _cnmf_rank_survey_science(
            adata,
            source=_source_string(expression),
            components_min=components_min,
            components_max=components_max,
            n_iter=n_iter,
            num_highvar_genes=num_highvar_genes,
            random_seed=random_seed,
        )
        try:
            metadata = run.metadata
            matrix = expression.matrix(adata)
            cell_totals, _ = matrix_totals_and_nonzero(matrix, axis=1)
            gene_totals, _ = matrix_totals_and_nonzero(matrix, axis=0)
            metric_records = [metric.as_dict() for metric in run.metrics]
            parameters = {
                **expression.parameters(),
                "components_min": metadata.candidate_ks[0],
                "components_max": metadata.candidate_ks[-1],
                "candidate_ks": list(metadata.candidate_ks),
                "n_iter": metadata.n_iter,
                "num_highvar_genes": metadata.num_highvar_genes,
                "random_seed": metadata.random_seed,
                "total_restarts": metadata.total_restarts,
                "source_features": metadata.input_genes,
                "current_features": metadata.current_features,
                "input_fingerprint": metadata.input_fingerprint,
                "k_metrics_fingerprint_sha256": _rank_metrics_fingerprint(metric_table),
                "fixed_policy": asdict(metadata.fixed_policy),
                "requested_resource": metadata.requested_resource.as_dict(),
            }
            warning_messages = [
                "OPENBIO_CNMF_RUN is a session-only native-directory artifact bound to its producing Python "
                "Worker. The Rank Survey worker closes its live run after encoding; each consumer receives a "
                "private writable checkout that is removed when that one-shot worker exits.",
                "The audited PyPI wheel SHA-256 is reference evidence only; this runtime did not attest the installed artifact.",
            ]
            warning_messages.extend(_warning_text(metadata.known_upstream_warnings))
            warning_messages.extend(metadata.input_advisories)
            if metadata.n_iter < 100:
                warning_messages.append(
                    "This survey uses fewer than the commonly documented 100 NMF restarts per K; stability evidence may be weak."
                )
            code = _survey_code(
                expression,
                components_min=metadata.candidate_ks[0],
                components_max=metadata.candidate_ks[-1],
                n_iter=metadata.n_iter,
                num_highvar_genes=metadata.num_highvar_genes,
                random_seed=metadata.random_seed,
            )
            results_text = (
                f"Completed {len(metadata.candidate_ks)} candidate ranks "
                f"({metadata.candidate_ks[0]}–{metadata.candidate_ks[-1]}) with {metadata.n_iter} complete "
                "NMF restarts per rank. Stability and prediction error were disclosed for analyst review; "
                "no rank was selected automatically and no final programs were calculated."
            )
            table_result = make_table_result(
                table=metric_table,
                title="cNMF rank-survey metrics",
                operation="cnmf_rank_survey",
                parameters=parameters,
                description=results_text,
                warnings=warning_messages,
                input_cells=metadata.input_cells,
                input_genes=metadata.input_genes,
                started_at=started_at,
                random_seed=metadata.random_seed,
            )
            report, code = make_analysis_report(
                node_id="OpenBioSingleCellCNMFRankSurvey",
                title="cNMF rank-survey summary",
                operation="cnmf_rank_survey",
                methods=(
                    f"The method authors' exact cnmf=={CNMF_REQUIRED_VERSION} file-backed implementation was run "
                    f"on declared non-negative expression source {_source_label(expression)!r}. A unique owned "
                    "temporary directory "
                    f"held {metadata.total_restarts} strict CPU/single-worker restarts, combined spectra, and no-filter "
                    "consensus statistics. Every official file path, axis, restart family, SHA-256 manifest, and "
                    "backend signature was validated."
                ),
                results=results_text,
                key_results={
                    "cells": metadata.input_cells,
                    "genes": metadata.input_genes,
                    "source": _source_label(expression),
                    "source_features": metadata.input_genes,
                    "current_features": metadata.current_features,
                    "input_advisories": list(metadata.input_advisories),
                    "total_counts": metadata.input_total_counts,
                    "nonzero_entries": metadata.input_nonzero_entries,
                    "cell_totals": summarize_numeric(cell_totals),
                    "gene_totals": summarize_numeric(gene_totals),
                    "realized_highvar_genes": len(metadata.realized_highvar_genes),
                    "positive_variance_highvar_genes": metadata.positive_variance_highvar_genes,
                    "candidate_ks": list(metadata.candidate_ks),
                    "k_metrics": metric_records,
                    "total_restarts": metadata.total_restarts,
                    "completed_restarts": sum(metric.completed_restarts for metric in run.metrics),
                    "requested_resource": metadata.requested_resource.as_dict(),
                    "realized_resource": metadata.realized_resource.as_dict(),
                    "input_fingerprint": metadata.input_fingerprint,
                    "parameters_fingerprint": metadata.parameters_fingerprint,
                    "backend": metadata.backend_name,
                    "backend_version": metadata.backend_version,
                    "backend_paths_fingerprint": metadata.backend_paths_fingerprint,
                    "prepare_manifest_sha256": metadata.prepare_manifest_sha256,
                    "factorization_manifest_sha256": metadata.factorization_manifest_sha256,
                    "artifact_manifest_sha256": metadata.artifact_manifest_sha256,
                    "integrity_manifest_artifact_count": len(metadata.artifact_hashes),
                    "cleanup_ownership": {
                        "owner": "one-shot Python Worker",
                        "live_run": "closed after native artifact encoding",
                        "consumer": "private writable checkout removed after the operation",
                        "node_graph": "session artifact retained while Classic cache references its ticket",
                        "cached_directory_may_persist": True,
                    },
                    "execution_mode": "CPU, single-worker, owned private file-backed run",
                    "audited_pypi_wheel_sha256": metadata.audited_pypi_wheel_sha256,
                    "installed_distribution_attested": metadata.installed_distribution_attested,
                    "known_upstream_warnings": [
                        {"id": label, "count": count} for label, count in metadata.known_upstream_warnings
                    ],
                },
                parameters=parameters,
                references=CNMF_REFERENCES,
                software_packages=CNMF_SOFTWARE_PACKAGES,
                warnings=warning_messages,
                limitations=(
                    "Rank selection remains an analyst decision supported by stability, reconstruction error, and biological interpretability.",
                    "The survey pools cells and may reflect Sample or Technical batch effects; it is not replicate-aware Condition inference.",
                    "The Worker-bound native run is session-only and cannot be persisted; export final annotated AnnData for reuse.",
                ),
                input_cells=metadata.input_cells,
                input_genes=metadata.input_genes,
                started_at=started_at,
                code=code,
                random_seed=metadata.random_seed,
            )
            return run, table_result, report, code
        except BaseException as primary:
            close_run_preserving_primary(
                run,
                primary,
                context="cNMF Rank Survey node-report rollback",
            )
            raise


class OpenBioSingleCellCNMF:
    @classmethod
    def execute(
        cls,
        run: CNMFRun,
        selected_k: int = 7,
        density_threshold: float = 2.0,
        local_neighborhood_size: float = CNMF_LOCAL_NEIGHBORHOOD_SIZE,
        n_top_genes: int = 100,
        overwrite_existing: bool = False,
    ) -> tuple[Any, Any, str]:
        if not isinstance(run, CNMFRun):
            raise TypeError("cNMF Consensus Programs requires an OPENBIO_CNMF_RUN from cNMF Rank Survey.")
        started_at = time.perf_counter()
        output = _cnmf_consensus_programs_science(
            run,
            selected_k=selected_k,
            density_threshold=density_threshold,
            local_neighborhood_size=local_neighborhood_size,
            n_top_genes=n_top_genes,
            overwrite_existing=overwrite_existing,
        )
        state = output.uns["cnmf"]
        survey = state["survey"]
        source = state["source"]
        source_label = f"layer:{source['layer']}" if source["kind"] == "layer" else source["kind"]
        parameters = {
            "source": source_label,
            "source_features": survey["source_features"],
            "current_features": survey["current_features"],
            "selected_k": state["selected_k"],
            "density_threshold": state["density_threshold"],
            "local_neighborhood_size": state["local_neighborhood_size"],
            "density_neighbors": state["density_neighbors"],
            "n_top_genes": n_top_genes,
            "overwrite_existing": bool(overwrite_existing),
            "survey_candidate_ks": state["candidate_ks"],
            "survey_n_iter": survey["n_iter"],
            "survey_num_highvar_genes": survey["num_highvar_genes"],
            "survey_realized_highvar_genes": survey["realized_highvar_genes"],
            "survey_random_seed": survey["random_seed"],
            "survey_backend": survey["backend"],
            "survey_backend_version": survey["backend_version"],
            "survey_execution_mode": "CPU, single-worker, owned private file-backed run",
            "survey_fixed_policy": survey["fixed_policy"],
            "survey_k_metrics": survey["k_metrics"],
            "survey_total_restarts": survey["total_restarts"],
            "survey_completed_restarts": survey["completed_restarts"],
            "survey_requested_resource": survey["requested_resource"],
            "survey_realized_resource": survey["realized_resource"],
            "input_fingerprint": survey["input_fingerprint"],
        }
        warning_messages: list[str] = []
        warning_messages.extend(str(message) for message in survey["input_advisories"])
        warning_messages.append(
            "The immutable OPENBIO_CNMF_RUN input remains owned by the Artifact Runtime. Consensus writes only to "
            "a private checkout, which the one-shot worker removes before returning."
        )
        if state["density_threshold"] >= 2.0:
            warning_messages.append("Density threshold 2.0 retained the documented no-filter initial consensus.")
        warning_records = tuple((item["id"], item["count"]) for item in state["known_upstream_warnings"])
        warning_messages.extend(_warning_text(warning_records))
        finish_adata(
            output,
            "cnmf_consensus_programs",
            parameters,
            run.metadata.input_cells,
            run.metadata.input_genes,
            started_at,
            random_seed=run.metadata.random_seed,
            warnings=warning_messages,
        )
        code = _consensus_code(
            selected_k=state["selected_k"],
            density_threshold=state["density_threshold"],
            local_neighborhood_size=state["local_neighborhood_size"],
            n_top_genes=n_top_genes,
            overwrite_existing=bool(overwrite_existing),
        )
        usage_summaries = {
            program: summarize_numeric(output.obsm["X_cnmf_usage"][program].to_numpy())
            for program in state["program_names"]
        }
        results_text = (
            f"Resolved {state['selected_k']} consensus gene-expression programs for {output.n_obs:,} cells. "
            f"Local-density filtering retained {state['components_after_filtering']:,} of "
            f"{state['components_before_filtering']:,} restart components "
            f"({state['components_filtered']:,} removed); program usages remain continuous mixtures."
        )
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellCNMF",
            title="cNMF consensus-program summary",
            operation="cnmf_consensus_programs",
            methods=(
                f"The method authors' exact cnmf=={CNMF_REQUIRED_VERSION} final consensus was applied in a private "
                f"writable checkout of the integrity-checked native run at analyst-selected K={state['selected_k']}. "
                "Local density and KMeans "
                "labels were independently reconstructed from the Survey-snapshotted merged spectrum, the exact upstream "
                "cache/result family was isolated, and the official four-tuple load_results output was aligned by identifiers."
            ),
            results=results_text,
            key_results={
                "cells": int(output.n_obs),
                "genes": int(output.n_vars),
                "source": source_label,
                "source_features": survey["source_features"],
                "current_features": survey["current_features"],
                "selected_k": state["selected_k"],
                "selected_k_survey_stability": state["selected_k_survey_stability"],
                "selected_k_survey_prediction_error": state["selected_k_survey_prediction_error"],
                "components_before_filtering": state["components_before_filtering"],
                "components_after_filtering": state["components_after_filtering"],
                "components_filtered": state["components_filtered"],
                "filtering_fraction": state["filtering_fraction"],
                "local_density": state["local_density_summary"],
                "cluster_component_counts": state["cluster_component_counts"],
                "local_neighborhood_size": state["local_neighborhood_size"],
                "density_neighbors": state["density_neighbors"],
                "usage_by_program": usage_summaries,
                "top_genes_by_program": {
                    program: genes[: min(10, len(genes))] for program, genes in state["top_genes"].items()
                },
                "storage_keys": state["storage_keys"],
                "overwrote_existing": state["overwrote_existing"],
                "input_fingerprint": survey["input_fingerprint"],
                "input_advisories": survey["input_advisories"],
                "backend": survey["backend"],
                "backend_version": survey["backend_version"],
                "audited_pypi_wheel_sha256": survey["audited_pypi_wheel_sha256"],
                "installed_distribution_attested": survey["installed_distribution_attested"],
                "survey_fixed_policy": survey["fixed_policy"],
                "survey_k_metrics": [metric.as_dict() for metric in run.metrics],
                "survey_total_restarts": survey["total_restarts"],
                "survey_completed_restarts": survey["completed_restarts"],
                "survey_requested_resource": survey["requested_resource"],
                "survey_realized_resource": survey["realized_resource"],
                "immutable_manifest_before_consensus": survey["immutable_manifest_before_consensus"],
                "immutable_manifest_after_consensus": survey["immutable_manifest_after_consensus"],
                "known_upstream_warnings": state["known_upstream_warnings"],
                "cleanup_ownership": {
                    "owner": "one-shot Python Worker",
                    "live_run": "closed after consensus output encoding",
                    "consumer": "private writable checkout removed after the operation",
                    "node_graph": "session artifact retained while Classic cache references its ticket",
                    "cached_directory_may_persist": True,
                },
            },
            parameters=parameters,
            references=CNMF_REFERENCES,
            software_packages=CNMF_SOFTWARE_PACKAGES,
            warnings=warning_messages,
            limitations=(
                "K and the local-density threshold are analyst decisions supported, not proven, by diagnostics.",
                "Gene-expression programs require biological annotation and external validation; usages are continuous mixtures, not cell clusters.",
                "Pooled-cell programs may reflect Sample or Technical batch effects and are not replicate-aware Condition inference.",
                "Non-negative factorization cannot directly represent repression, and low nonzero usages may reflect overfitting.",
            ),
            input_cells=run.metadata.input_cells,
            input_genes=run.metadata.input_genes,
            started_at=started_at,
            code=code,
            random_seed=run.metadata.random_seed,
        )
        return output, report, code


def cnmf_rank_plot_owned(k_metrics: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    from .contracts import TableResult

    if not isinstance(k_metrics, TableResult):
        raise TypeError("cNMF Rank Plot requires the k_metrics TableResult from cNMF Rank Survey.")
    source_parameters = k_metrics.source.get("parameters") if isinstance(k_metrics.source, dict) else None
    provenance_keys = (
        "candidate_ks",
        "components_min",
        "components_max",
        "n_iter",
        "total_restarts",
        "source_features",
        "input_fingerprint",
        "k_metrics_fingerprint_sha256",
    )
    if not isinstance(source_parameters, dict) or any(
        k_metrics.parameters.get(key) != source_parameters.get(key) for key in provenance_keys
    ):
        raise ValueError("cNMF Rank Plot requires internally consistent k_metrics producer provenance.")
    top_policy = k_metrics.parameters.get("fixed_policy")
    source_policy = source_parameters.get("fixed_policy")
    if (
        not isinstance(top_policy, dict)
        or not isinstance(source_policy, dict)
        or top_policy.get("statistics_density_threshold")
        != source_policy.get("statistics_density_threshold")
    ):
        raise ValueError("cNMF Rank Plot requires internally consistent k_metrics producer provenance.")
    producer = {
        "operation": k_metrics.source.get("operation"),
        "parameters": k_metrics.parameters,
        "input_cells": k_metrics.input_cells,
        "input_genes": k_metrics.input_genes,
    }
    started_at = time.perf_counter()
    png, details = run_cnmf_rank_plot(k_metrics.table, producer=producer, view=view)
    parameters = {"view": details["view_parameters"]}
    code = cnmf_rank_plot_code(producer=producer, view=details["view_parameters"])
    description = (
        f"Rendered stored cNMF rank evidence for {details['plotted_ranks']:,} complete candidate K values. "
        "No rank was selected or recommended."
    )
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": [],
        "input_cells": k_metrics.input_cells,
        "input_genes": k_metrics.input_genes,
        "started_at": started_at,
        "random_seed": k_metrics.random_seed,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="cnmf_rank_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCNMFRankPlot",
        title="cNMF rank plot summary",
        operation="cnmf_rank_plot",
        methods=(
            "Validated the exact cNMF Rank Survey producer, canonical K-metrics schema, source fingerprint, "
            "candidate-rank axis, and complete restart family before read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=(*CNMF_REFERENCES, MATPLOTLIB_REFERENCE),
        software_packages=(*CNMF_SOFTWARE_PACKAGES, "matplotlib"),
        warnings=[],
        limitations=(
            "Rank selection remains an analyst decision supported by stability, reconstruction error, and biological interpretability.",
            "The plot is diagnostic and pooled-cell; it is not replicate-aware Condition inference.",
        ),
        input_cells=k_metrics.input_cells,
        input_genes=k_metrics.input_genes,
        started_at=started_at,
        code=code,
        random_seed=k_metrics.random_seed,
    )
    return plotted, report, code


def cnmf_programs_plot_owned(adata: Any, *, view: Any = None) -> tuple[Any, Any, str]:
    started_at = time.perf_counter()
    png, details = run_cnmf_programs_plot(adata, view=view)
    parameters = {"view": details["view_parameters"]}
    code = cnmf_programs_plot_code(view=details["view_parameters"])
    description = (
        f"Rendered stored cNMF {details['view'].replace('_', ' ')} for {details['plotted_programs']:,} "
        f"continuous gene-expression programs across {details['plotted_cells']:,} cells."
    )
    common = {
        "parameters": parameters,
        "description": description,
        "warnings": [],
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "started_at": started_at,
    }
    plotted = make_plot_result(
        png=png,
        title=details["title"],
        operation="cnmf_programs_plot",
        **common,
    )
    report, code = make_analysis_report(
        node_id="OpenBioSingleCellCNMFProgramsPlot",
        title="cNMF programs plot summary",
        operation="cnmf_programs_plot",
        methods=(
            "Validated cNMF Consensus Programs provenance, observation and feature axes, normalized usage, GEP "
            "scores, stored top-gene ordering, and all content fingerprints before read-only Matplotlib rendering."
        ),
        results=description,
        key_results={key: value for key, value in details.items() if key not in {"view_parameters", "title"}},
        parameters=parameters,
        references=(*CNMF_REFERENCES, MATPLOTLIB_REFERENCE),
        software_packages=(*CNMF_SOFTWARE_PACKAGES, "matplotlib"),
        warnings=[],
        limitations=(
            "Gene-expression programs require biological annotation and external validation; usages are continuous mixtures, not cell clusters.",
            "Pooled-cell programs may reflect Sample or Technical batch effects and are not replicate-aware Condition inference.",
        ),
        input_cells=int(adata.n_obs),
        input_genes=int(adata.n_vars),
        started_at=started_at,
        code=code,
    )
    return plotted, report, code


CNMF_RUN_KIND = "OPENBIO_CNMF_RUN"
TABLE_KIND = "OPENBIO_SINGLE_CELL_TABLE"
PLOT_KIND = "OPENBIO_SINGLE_CELL_PLOT"


@register_operation("openbio.node.cnmfranksurvey")
def cnmf_rank_survey_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="cNMF Rank Survey")
    require_parameters(
        parameters,
        {
            "source",
            "components_min",
            "components_max",
            "n_iter",
            "num_highvar_genes",
            "random_seed",
        },
        operation="cNMF Rank Survey",
    )
    adata = read_anndata_input(inputs)
    run, table, summary, code = OpenBioSingleCellCNMFRankSurvey.execute(adata, **parameters)
    try:
        run_root = context.create_output_directory("run")
        # Encode the worker-owned run directly; the codec deliberately avoids copy_base_adata().
        write_cnmf_run(run_root, run)
        records = analysis_outputs(
            summary,
            code,
            {
                "type": "artifact",
                "name": "run",
                "kind": CNMF_RUN_KIND,
                "codec": CNMF_NATIVE_CODEC,
                "payload": run_root.relative_to(context.output_root).as_posix(),
            },
            write_table_output(context, table, name="k_metrics", kind=TABLE_KIND),
        )
    except BaseException as primary:
        close_run_preserving_primary(run, primary, context="cNMF Rank Survey artifact encoding")
        raise
    run.close()
    return records


@register_operation("openbio.node.cnmfrankplot")
def cnmf_rank_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"k_metrics"}, operation="cNMF Rank Plot")
    require_parameters(parameters, {"view"}, operation="cNMF Rank Plot")
    table, metadata = read_table_input(inputs, "k_metrics", kind=TABLE_KIND)
    result = table_from_metadata(metadata, table)
    plotted, report, code = cnmf_rank_plot_owned(result, **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


@register_operation("openbio.node.cnmf")
def cnmf_operation(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"run"}, operation="cNMF Consensus Programs")
    require_parameters(
        parameters,
        {
            "selected_k",
            "density_threshold",
            "local_neighborhood_size",
            "n_top_genes",
            "overwrite_existing",
        },
        operation="cNMF Consensus Programs",
    )
    run_root = require_artifact_input(
        inputs,
        "run",
        kind=CNMF_RUN_KIND,
        codec=CNMF_NATIVE_CODEC,
    )
    with checkout_cnmf_run(run_root, checkout_parent=context.output_root) as run:
        # The consensus output is a scientifically new AnnData materialized from the saved Survey base.
        output, summary, code = OpenBioSingleCellCNMF.execute(run, **parameters)
        return analysis_outputs(summary, code, write_anndata_output(context, output))


@register_operation("openbio.node.cnmfprogramsplot")
def cnmf_programs_plot(
    context: OperationContext,
    inputs: dict[str, JSONValue],
    parameters: dict[str, JSONValue],
) -> list[JSONValue]:
    require_input_names(inputs, {"adata"}, operation="cNMF Programs Plot")
    require_parameters(parameters, {"view"}, operation="cNMF Programs Plot")
    plotted, report, code = cnmf_programs_plot_owned(read_anndata_input(inputs), **parameters)
    return analysis_outputs(report, code, write_plot_output(context, plotted, kind=PLOT_KIND))


# Stable public names match operation IDs while avoiding collisions with the standalone science functions.
cnmf_rank_survey = cnmf_rank_survey_operation
cnmf = cnmf_operation


__all__ = [
    "CNMF_RUN_KIND",
    "OpenBioSingleCellCNMF",
    "OpenBioSingleCellCNMFRankSurvey",
    "cnmf",
    "cnmf_programs_plot",
    "cnmf_programs_plot_owned",
    "cnmf_rank_plot",
    "cnmf_rank_plot_owned",
    "cnmf_rank_survey",
]
