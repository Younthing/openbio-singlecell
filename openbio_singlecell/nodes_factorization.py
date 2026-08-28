from __future__ import annotations

import time
from dataclasses import asdict
from typing import TYPE_CHECKING

from comfy_api.latest import io

from .analysis_reporting import AnalysisReference, make_analysis_report, summarize_numeric
from .analysis_utils import finish_adata, make_table_result, matrix_totals_and_nonzero
from .cnmf_run import CNMFRun
from .cnmf_standalone import (
    CNMF_AUDITED_PYPI_WHEEL_SHA256,
    CNMF_LOCAL_NEIGHBORHOOD_SIZE,
    CNMF_REQUIRED_VERSION,
    close_run_preserving_primary,
    cnmf_consensus_programs,
    cnmf_rank_survey,
    equivalent_source,
)
from .expression_source import DynamicExpressionSource, ExpressionSource, ExpressionSourceSpec
from .node_types import AnnDataType, CNMFRunType, SummaryResultType, TableResultType

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


class OpenBioSingleCellCNMFRankSurvey(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="cNMF count source",
        default="layer",
        include_raw=True,
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNMFRankSurvey",
            display_name="cNMF Rank Survey",
            category=CATEGORY,
            description="Run the complete official cNMF rank family in an owned private file-backed artifact.",
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Int.Input("components_min", default=3, min=2),
                io.Int.Input("components_max", default=19, min=2),
                io.Int.Input("n_iter", default=100, min=2),
                io.Int.Input("num_highvar_genes", default=2000, min=1),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                CNMFRunType.Output(display_name="run"),
                TableResultType.Output(display_name="k_metrics"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

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
    ) -> io.NodeOutput:
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        started_at = time.perf_counter()
        run, metric_table = cnmf_rank_survey(
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
                "fixed_policy": asdict(metadata.fixed_policy),
                "requested_resource": metadata.requested_resource.as_dict(),
            }
            warning_messages = [
                "OPENBIO_CNMF_RUN is process-local and non-serializable. Equivalent Python callers should call "
                "run.close() after the final consensus; node-graph storage is released by the finalizer only after "
                "the engine evicts all cached references, so the private directory may persist until cache eviction "
                "or process exit.",
                "The audited PyPI wheel SHA-256 is reference evidence only; this runtime did not attest the installed artifact.",
            ]
            warning_messages.extend(_warning_text(metadata.known_upstream_warnings))
            warning_messages.extend(metadata.input_advisories)
            if metadata.source_state == "unknown":
                warning_messages.append(
                    "OpenBio provenance cannot establish the selected source's measurement scale or transformation "
                    "history; confirm the external data contract before interpretation."
                )
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
                    "source_state": metadata.source_state,
                    "source_state_evidence": metadata.source_state_evidence,
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
                        "owner": "OPENBIO_CNMF_RUN",
                        "explicit_python": "run.close() or context-manager exit",
                        "node_graph": "finalizer after engine cache/reference release",
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
                    "Numeric validation alone cannot establish an externally supplied matrix's measurement scale or transformation provenance.",
                    "The survey pools cells and may reflect Sample or Technical batch effects; it is not replicate-aware Condition inference.",
                    "The typed run is process-local; only final annotated AnnData should be persisted across application restarts.",
                ),
                input_cells=metadata.input_cells,
                input_genes=metadata.input_genes,
                started_at=started_at,
                code=code,
                random_seed=metadata.random_seed,
            )
            return io.NodeOutput(run, table_result, report, code)
        except BaseException as primary:
            close_run_preserving_primary(
                run,
                primary,
                context="cNMF Rank Survey node-report rollback",
            )
            raise


class OpenBioSingleCellCNMF(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCNMF",
            display_name="cNMF Consensus Programs",
            category=CATEGORY,
            description="Apply one surveyed rank/density choice to a live owned official cNMF run.",
            inputs=[
                CNMFRunType.Input("run"),
                io.Int.Input("selected_k", default=7, min=2),
                io.Float.Input("density_threshold", default=2.0, min=0.0, max=2.0, step=0.05),
                io.Float.Input(
                    "local_neighborhood_size",
                    default=CNMF_LOCAL_NEIGHBORHOOD_SIZE,
                    min=0.001,
                    max=1.0,
                    step=0.05,
                    advanced=True,
                ),
                io.Int.Input("n_top_genes", default=100, min=1, max=2**31 - 1, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        run: CNMFRun,
        selected_k: int = 7,
        density_threshold: float = 2.0,
        local_neighborhood_size: float = CNMF_LOCAL_NEIGHBORHOOD_SIZE,
        n_top_genes: int = 100,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        if not isinstance(run, CNMFRun):
            raise TypeError("cNMF Consensus Programs requires an OPENBIO_CNMF_RUN from cNMF Rank Survey.")
        started_at = time.perf_counter()
        output = cnmf_consensus_programs(
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
        if survey["source_state"] == "unknown":
            warning_messages.append(
                "OpenBio provenance cannot establish the selected source's measurement scale or transformation "
                "history; confirm the external data contract before interpretation."
            )
        warning_messages.append(
            "The reusable OPENBIO_CNMF_RUN remains owned by the upstream node. Equivalent Python callers should "
            "close it after their final decision; a node graph releases its private directory after cached references "
            "are evicted, which may be as late as process exit."
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
                f"The method authors' exact cnmf=={CNMF_REQUIRED_VERSION} final consensus was applied to the live "
                f"integrity-checked staged run at analyst-selected K={state['selected_k']}. Local density and KMeans "
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
                "source_state": survey["source_state"],
                "source_state_evidence": survey["source_state_evidence"],
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
                    "owner": "OPENBIO_CNMF_RUN",
                    "explicit_python": "run.close() or context-manager exit",
                    "node_graph": "finalizer after engine cache/reference release",
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
        return io.NodeOutput(output, report, code)


FACTORIZATION_NODE_CLASSES = [OpenBioSingleCellCNMFRankSurvey, OpenBioSingleCellCNMF]


__all__ = [
    "FACTORIZATION_NODE_CLASSES",
    "OpenBioSingleCellCNMF",
    "OpenBioSingleCellCNMFRankSurvey",
]
