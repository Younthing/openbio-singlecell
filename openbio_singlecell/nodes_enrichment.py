from __future__ import annotations

import hashlib
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import PLUGIN_VERSION, dependencies
from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .dgidb_resource import (
    DGIDB_RESOURCE_EXTENSIONS,
    DGIdbResource,
    dgidb_resource_cache_fingerprint,
    dgidb_resource_code,
    load_dgidb_resource,
    validate_dgidb_resource,
)
from .drug_enrichment import drug_gsea_code, drug_ora_code, run_drug_gsea, run_drug_ora
from .drug_score import drug_score_code, run_drug_score
from .enrichment_artifacts import validate_enrichment_artifact_pair
from .expression_source import (
    DynamicExpressionSource,
    ExpressionSourceSpec,
)
from .files import input_file_fingerprint, resolve_input_path
from .gene_set_scoring import (
    gene_set_resource_cache_fingerprint,
    gene_set_scoring_code,
    run_aucell_scores,
    run_gene_panel_score,
    run_gsva_scores,
)
from .generic_ora import build_generic_ora_summary, generic_ora_code, run_generic_ora_evidence
from .node_types import AnnDataType, DGIdbResourceType, SummaryResultType, TableResultType
from .pathway_score_contrast import pathway_score_contrast_code, run_pathway_score_contrast
from .ranked_enrichment import build_ranked_gsea_summary, ranked_gsea_code, run_ranked_gsea_evidence

if TYPE_CHECKING:
    from anndata import AnnData

    from .contracts import TableResult


CATEGORY = "openbio/single-cell/enrichment"
GENE_SET_EXTENSIONS = (".csv", ".tsv", ".gmt")


class OpenBioSingleCellAUCellScores(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="AUCell expression source",
        include_raw=True,
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellAUCellScores",
            display_name="AUCell Scores",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gene_sets_file", default=""),
                io.String.Input("resource_metadata_json", default="{}"),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.Int.Input("min_targets", default=5, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("n_top_features", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("batch_size", default=250_000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_output_rows", default=2_000_000, min=1, max=2**31 - 1, advanced=True),
                io.Float.Input("max_working_memory_gib", default=4.0, min=0.001, max=1024.0, advanced=True),
                io.String.Input("output_key", default="aucell_scores", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        return gene_set_resource_cache_fingerprint(path, identity)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        gene_sets_file: str = "",
        resource_metadata_json: str = "{}",
        source: DynamicExpressionSource | None = None,
        source_column: str = "geneset",
        target_column: str = "genesymbol",
        min_targets: int = 5,
        n_top_features: int = 0,
        batch_size: int = 250_000,
        max_output_rows: int = 2_000_000,
        max_working_memory_gib: float = 4.0,
        output_key: str = "aucell_scores",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        resolved_path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        output, summary = run_aucell_scores(
            adata,
            gene_sets_path=resolved_path,
            requested_resource_path=gene_sets_file,
            resource_metadata_json=resource_metadata_json,
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            source_column=source_column,
            target_column=target_column,
            min_targets=min_targets,
            n_top_features=n_top_features,
            batch_size=batch_size,
            max_output_rows=max_output_rows,
            max_working_memory_gib=max_working_memory_gib,
            output_key=output_key,
            overwrite_existing=overwrite_existing,
            openbio_version=PLUGIN_VERSION,
        )
        report = make_summary_result(
            summary=summary,
            title="AUCell observation scores",
            operation="aucell_scores",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        code = gene_set_scoring_code(
            function_name="score_aucell",
            parameters={
                "method": "aucell",
                "gene_sets_path": resolved_path,
                "requested_resource_path": gene_sets_file,
                "resource_metadata_json": resource_metadata_json,
                "source_kind": expression.kind,
                "layer_name": expression.layer_name,
                "source_column": source_column,
                "target_column": target_column,
                "min_targets": min_targets,
                "output_key": output_key,
                "overwrite_existing": overwrite_existing,
                "expected_resource_sha256": summary["key_results"]["resource"]["sha256"],
                "n_top_features": n_top_features,
                "batch_size": batch_size,
                "max_output_rows": max_output_rows,
                "max_working_memory_gib": max_working_memory_gib,
                "openbio_version": PLUGIN_VERSION,
            },
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellGSVAScores(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="GSVA expression source",
        include_raw=True,
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellGSVAScores",
            display_name="GSVA Scores",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gene_sets_file", default=""),
                io.String.Input("resource_metadata_json", default="{}"),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.Int.Input("min_targets", default=10, min=1, max=2**31 - 1, advanced=True),
                io.Combo.Input(
                    "kernel",
                    options=["gaussian_normalized", "poisson_counts", "empirical"],
                    default="gaussian_normalized",
                ),
                io.Boolean.Input("maxdiff", default=True, advanced=True),
                io.Boolean.Input("absrnk", default=False, advanced=True),
                io.Float.Input("tau", default=1.0, min=0.001, max=1000.0, advanced=True),
                io.Float.Input("max_working_memory_gib", default=4.0, min=0.001, max=1024.0, advanced=True),
                io.Int.Input("batch_size", default=250_000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_output_rows", default=2_000_000, min=1, max=2**31 - 1, advanced=True),
                io.String.Input("output_key", default="gsva_scores", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        return gene_set_resource_cache_fingerprint(path, identity)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        gene_sets_file: str = "",
        resource_metadata_json: str = "{}",
        source: DynamicExpressionSource | None = None,
        source_column: str = "geneset",
        target_column: str = "genesymbol",
        min_targets: int = 10,
        kernel: str = "gaussian_normalized",
        maxdiff: bool = True,
        absrnk: bool = False,
        tau: float = 1.0,
        max_working_memory_gib: float = 4.0,
        batch_size: int = 250_000,
        max_output_rows: int = 2_000_000,
        output_key: str = "gsva_scores",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        resolved_path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        output, summary = run_gsva_scores(
            adata,
            gene_sets_path=resolved_path,
            requested_resource_path=gene_sets_file,
            resource_metadata_json=resource_metadata_json,
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            source_column=source_column,
            target_column=target_column,
            min_targets=min_targets,
            kernel=kernel,
            maxdiff=maxdiff,
            absrnk=absrnk,
            tau=tau,
            batch_size=batch_size,
            max_output_rows=max_output_rows,
            max_working_memory_gib=max_working_memory_gib,
            output_key=output_key,
            overwrite_existing=overwrite_existing,
            openbio_version=PLUGIN_VERSION,
        )
        report = make_summary_result(
            summary=summary,
            title="GSVA observation scores",
            operation="gsva_scores",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        code = gene_set_scoring_code(
            function_name="score_gsva",
            parameters={
                "method": "gsva",
                "gene_sets_path": resolved_path,
                "requested_resource_path": gene_sets_file,
                "resource_metadata_json": resource_metadata_json,
                "source_kind": expression.kind,
                "layer_name": expression.layer_name,
                "source_column": source_column,
                "target_column": target_column,
                "min_targets": min_targets,
                "output_key": output_key,
                "overwrite_existing": overwrite_existing,
                "expected_resource_sha256": summary["key_results"]["resource"]["sha256"],
                "kernel": kernel,
                "maxdiff": maxdiff,
                "absrnk": absrnk,
                "tau": tau,
                "batch_size": batch_size,
                "max_output_rows": max_output_rows,
                "max_working_memory_gib": max_working_memory_gib,
                "openbio_version": PLUGIN_VERSION,
            },
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellGenePanelScores(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Gene panel expression source",
        include_raw=True,
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellGenePanelScores",
            display_name="Gene Panel Score",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("gene_sets_file", default=""),
                io.String.Input("resource_metadata_json", default="{}"),
                io.String.Input("panel", default=""),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("source_column", default="geneset", advanced=True),
                io.String.Input("target_column", default="genesymbol", advanced=True),
                io.String.Input("output_key", default="panel_score", advanced=True),
                io.Int.Input("ctrl_size", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input("n_bins", default=25, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=0, min=0, max=2**31 - 1, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        return gene_set_resource_cache_fingerprint(path, identity)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        gene_sets_file: str = "",
        resource_metadata_json: str = "{}",
        panel: str = "",
        source: DynamicExpressionSource | None = None,
        source_column: str = "geneset",
        target_column: str = "genesymbol",
        output_key: str = "panel_score",
        ctrl_size: int = 0,
        n_bins: int = 25,
        random_seed: int = 0,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        resolved_path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        output, summary = run_gene_panel_score(
            adata,
            gene_sets_path=resolved_path,
            requested_resource_path=gene_sets_file,
            resource_metadata_json=resource_metadata_json,
            panel=panel,
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            source_column=source_column,
            target_column=target_column,
            output_key=output_key,
            ctrl_size=ctrl_size,
            n_bins=n_bins,
            random_seed=random_seed,
            overwrite_existing=overwrite_existing,
            openbio_version=PLUGIN_VERSION,
        )
        report = make_summary_result(
            summary=summary,
            title=f"Gene panel score: {panel}",
            operation="gene_panel_score",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
        )
        code = gene_set_scoring_code(
            function_name="score_gene_panel",
            parameters={
                "method": "panel",
                "gene_sets_path": resolved_path,
                "requested_resource_path": gene_sets_file,
                "resource_metadata_json": resource_metadata_json,
                "source_kind": expression.kind,
                "layer_name": expression.layer_name,
                "source_column": source_column,
                "target_column": target_column,
                "min_targets": 1,
                "output_key": output_key,
                "overwrite_existing": overwrite_existing,
                "expected_resource_sha256": summary["key_results"]["resource"]["sha256"],
                "panel": panel,
                "ctrl_size": ctrl_size,
                "n_bins": n_bins,
                "random_seed": random_seed,
                "openbio_version": PLUGIN_VERSION,
            },
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellPathwayScoreTTest(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPathwayScoreTTest",
            display_name="Pathway Score T-Test",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("annotation_key", default="cell_type"),
                io.String.Input("population", default=""),
                io.String.Input("condition_a", default=""),
                io.String.Input("condition_b", default=""),
                io.String.Input("score_key", default="aucell_scores"),
                io.Int.Input("min_cells_per_sample_population", default=10, min=1, max=2**31 - 1),
                io.Int.Input("min_samples_per_condition", default=3, min=2, max=2**31 - 1),
                io.String.Input("technical_batch_key", default="", advanced=True),
                io.Float.Input("confidence_level", default=0.95, min=0.5, max=0.999, advanced=True),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                    advanced=True,
                ),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        condition_key: str = "condition",
        annotation_key: str = "cell_type",
        population: str = "",
        condition_a: str = "",
        condition_b: str = "",
        score_key: str = "aucell_scores",
        min_cells_per_sample_population: int = 10,
        min_samples_per_condition: int = 3,
        technical_batch_key: str = "",
        confidence_level: float = 0.95,
        annotation_status: str = "unknown",
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, summary = run_pathway_score_contrast(
            adata,
            sample_key=sample_key,
            condition_key=condition_key,
            annotation_key=annotation_key,
            population=population,
            condition_a=condition_a,
            condition_b=condition_b,
            score_key=score_key,
            min_cells_per_sample_population=min_cells_per_sample_population,
            min_samples_per_condition=min_samples_per_condition,
            technical_batch_key=technical_batch_key,
            confidence_level=confidence_level,
            annotation_status=annotation_status,
            openbio_version=PLUGIN_VERSION,
        )
        result = make_table_result(
            title=f"{condition_a} vs {condition_b} pathway scores in {population}",
            operation="pathway_score_sample_welch",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title=f"Pathway score contrast: {condition_a} vs {condition_b}",
            operation="pathway_score_sample_welch",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        code = pathway_score_contrast_code(
            sample_key=sample_key,
            condition_key=condition_key,
            annotation_key=annotation_key,
            population=population,
            condition_a=condition_a,
            condition_b=condition_b,
            score_key=score_key,
            min_cells_per_sample_population=min_cells_per_sample_population,
            min_samples_per_condition=min_samples_per_condition,
            technical_batch_key=technical_batch_key,
            confidence_level=confidence_level,
            annotation_status=annotation_status,
            openbio_version=PLUGIN_VERSION,
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellRankedGSEA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRankedGSEA",
            display_name="Ranked GSEA",
            category=CATEGORY,
            inputs=[
                TableResultType.Input("table"),
                TableResultType.Input("universe"),
                io.String.Input("gene_sets_file", default="openbio-singlecell/gene_sets.csv"),
                io.String.Input("resource_metadata_json", default="{}"),
                io.String.Input("comparison", default=""),
                io.String.Input("gene_column", default="gene", advanced=True),
                io.String.Input("score_column", default="score", advanced=True),
                io.String.Input("source_column", default="source", advanced=True),
                io.String.Input("target_column", default="target", advanced=True),
                io.Int.Input("min_targets", default=15, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_targets", default=500, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("n_permutations", default=1000, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input(
                    "max_output_rows",
                    default=100_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        digest = hashlib.sha256()
        with open(path, "rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        return ("openbio-ranked-gsea-resource-v1", *identity, digest.hexdigest())

    @classmethod
    def execute(
        cls,
        table: TableResult,
        universe: TableResult,
        gene_sets_file: str = "openbio-singlecell/gene_sets.csv",
        resource_metadata_json: str = "{}",
        comparison: str = "",
        gene_column: str = "gene",
        score_column: str = "score",
        source_column: str = "source",
        target_column: str = "target",
        min_targets: int = 15,
        max_targets: int = 500,
        n_permutations: int = 1000,
        random_seed: int = 123,
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        provenance = validate_enrichment_artifact_pair(
            table,
            universe,
            purpose="ranked",
            selector=comparison,
            gene_column=gene_column,
            score_column=score_column,
            np=science.np,
            pd=science.pd,
        )
        resolved_resource = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        runtime_parameters = {
            "artifact_family": provenance["artifact_family"],
            "comparison": provenance["comparison"],
            "comparison_column": provenance["comparison_column"],
            "gene_column": provenance["gene_column"],
            "score_column": provenance["score_column"],
            "evidence_scope": provenance["evidence_scope"],
            "inference_unit": provenance["inference_unit"],
            "direction": provenance["direction"],
            "replicate_aware": provenance["replicate_aware"],
            "technical_batch_handling": provenance["technical_batch_handling"],
            "upstream_parameters": provenance["upstream_parameters"],
            "provenance_warnings": provenance["provenance_warnings"],
            "provenance_declarations": provenance["provenance_declarations"],
            "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
            "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
            "expected_universe_fingerprint": provenance["universe_fingerprint"],
            "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
            "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
            "resource_path": resolved_resource,
            "requested_resource_path": gene_sets_file,
            "resource_metadata_json": resource_metadata_json,
            "source_column": source_column,
            "target_column": target_column,
            "min_targets": min_targets,
            "max_targets": max_targets,
            "n_permutations": n_permutations,
            "random_seed": random_seed,
            "max_output_rows": max_output_rows,
            "openbio_version": PLUGIN_VERSION,
        }
        evidence, diagnostics = run_ranked_gsea_evidence(
            table,
            universe,
            **runtime_parameters,
            decoupler_module=None,
        )
        summary = build_ranked_gsea_summary(diagnostics)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        result = make_table_result(
            title=f"Ranked GSEA: {provenance['comparison']}",
            operation="ranked_gsea",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
            random_seed=random_seed,
            table=evidence,
        )
        report = make_summary_result(
            summary=summary,
            title=f"Ranked GSEA summary: {provenance['comparison']}",
            operation="ranked_gsea",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        code = ranked_gsea_code(
            **runtime_parameters,
            expected_resource_sha256=diagnostics["resource"]["sha256"],
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellGeneSetOverrepresentation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellGeneSetOverrepresentation",
            display_name="Gene Set Overrepresentation",
            category=CATEGORY,
            inputs=[
                TableResultType.Input("table"),
                TableResultType.Input("universe"),
                io.String.Input("gene_sets_file", default="openbio-singlecell/gene_sets.csv"),
                io.String.Input("resource_metadata_json", default="{}"),
                io.String.Input("comparison", default=""),
                io.String.Input("gene_column", default="gene", advanced=True),
                io.String.Input("source_column", default="source", advanced=True),
                io.String.Input("target_column", default="target", advanced=True),
                io.Int.Input("min_targets", default=3, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("min_overlap", default=2, min=1, max=2**31 - 1),
                io.Float.Input(
                    "max_p_adjusted",
                    default=0.05,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
                io.Int.Input(
                    "max_output_rows",
                    default=100_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        digest = hashlib.sha256()
        with open(path, "rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        return ("openbio-generic-ora-resource-v1", *identity, digest.hexdigest())

    @classmethod
    def execute(
        cls,
        table: TableResult,
        universe: TableResult,
        gene_sets_file: str = "openbio-singlecell/gene_sets.csv",
        resource_metadata_json: str = "{}",
        comparison: str = "",
        gene_column: str = "gene",
        source_column: str = "source",
        target_column: str = "target",
        min_targets: int = 3,
        min_overlap: int = 2,
        max_p_adjusted: float = 0.05,
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        provenance = validate_enrichment_artifact_pair(
            table,
            universe,
            purpose="selected",
            selector=comparison,
            gene_column=gene_column,
            score_column=None,
            np=science.np,
            pd=science.pd,
        )
        resolved_resource = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        runtime_parameters = {
            "artifact_family": provenance["artifact_family"],
            "comparison": provenance["comparison"],
            "comparison_column": provenance["comparison_column"],
            "gene_column": provenance["gene_column"],
            "evidence_scope": provenance["evidence_scope"],
            "inference_unit": provenance["inference_unit"],
            "direction": provenance["direction"],
            "replicate_aware": provenance["replicate_aware"],
            "technical_batch_handling": provenance["technical_batch_handling"],
            "upstream_parameters": provenance["upstream_parameters"],
            "provenance_warnings": provenance["provenance_warnings"],
            "provenance_declarations": provenance["provenance_declarations"],
            "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
            "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
            "expected_universe_fingerprint": provenance["universe_fingerprint"],
            "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
            "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
            "resource_path": resolved_resource,
            "requested_resource_path": gene_sets_file,
            "resource_metadata_json": resource_metadata_json,
            "source_column": source_column,
            "target_column": target_column,
            "min_targets": min_targets,
            "min_overlap": min_overlap,
            "max_p_adjusted": max_p_adjusted,
            "max_output_rows": max_output_rows,
            "openbio_version": PLUGIN_VERSION,
        }
        evidence, diagnostics = run_generic_ora_evidence(
            table,
            universe,
            **runtime_parameters,
            decoupler_module=None,
        )
        summary = build_generic_ora_summary(diagnostics)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        result = make_table_result(
            title=f"Gene-set overrepresentation: {provenance['comparison']}",
            operation="gene_set_overrepresentation",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
            table=evidence,
        )
        report = make_summary_result(
            summary=summary,
            title=f"Gene-set overrepresentation summary: {provenance['comparison']}",
            operation="gene_set_overrepresentation",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
        )
        code = generic_ora_code(
            **runtime_parameters,
            expected_resource_sha256=diagnostics["resource"]["sha256"],
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellDGIdbAnnotation(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDGIdbAnnotation",
            display_name="Load DGIdb Resource",
            category=CATEGORY,
            description="Load and validate one reviewed, pinned local DGIdb snapshot without network access.",
            inputs=[
                io.String.Input("resource_file", default="openbio-singlecell/dgidb.tsv"),
                io.String.Input("resource_metadata_json", default="{}"),
                io.String.Input("drug_column", default="drug_claim_name", advanced=True),
                io.String.Input("gene_column", default="gene_claim_name", advanced=True),
                io.String.Input("source_column", default="interaction_claim_source", advanced=True),
                io.String.Input("evidence_column", default="interaction_types", advanced=True),
                io.Int.Input(
                    "max_file_bytes",
                    default=536_870_912,
                    min=1,
                    max=2**63 - 1,
                    advanced=True,
                ),
                io.Int.Input("max_rows", default=2_000_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=[
                DGIdbResourceType.Output(display_name="resource"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, resource_file: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(resource_file, extensions=DGIDB_RESOURCE_EXTENSIONS)
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, resource_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(resource_file, extensions=DGIDB_RESOURCE_EXTENSIONS)
        identity = input_file_fingerprint(resource_file, DGIDB_RESOURCE_EXTENSIONS)
        return dgidb_resource_cache_fingerprint(
            path,
            identity,
            max_file_bytes=kwargs.get("max_file_bytes", 536_870_912),
        )

    @classmethod
    def execute(
        cls,
        resource_file: str = "openbio-singlecell/dgidb.tsv",
        resource_metadata_json: str = "{}",
        drug_column: str = "drug_claim_name",
        gene_column: str = "gene_claim_name",
        source_column: str = "interaction_claim_source",
        evidence_column: str = "interaction_types",
        max_file_bytes: int = 536_870_912,
        max_rows: int = 2_000_000,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        resolved = resolve_input_path(resource_file, extensions=DGIDB_RESOURCE_EXTENSIONS)
        resource, summary = load_dgidb_resource(
            resolved,
            requested_path=resource_file,
            resource_metadata_json=resource_metadata_json,
            drug_column=drug_column,
            gene_column=gene_column,
            source_column=source_column,
            evidence_column=evidence_column,
            max_file_bytes=max_file_bytes,
            max_rows=max_rows,
            openbio_version=PLUGIN_VERSION,
        )
        _, _, accounting, artifact = validate_dgidb_resource(resource)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        report = make_summary_result(
            summary=summary,
            title=f"DGIdb resource {summary['key_results']['resource']['metadata']['version']}",
            operation="load_dgidb_resource",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=0,
            input_genes=int(accounting["gene_count"]),
            started_at=started_at,
        )
        code = dgidb_resource_code(
            requested_path=resource_file,
            resource_metadata_json=resource_metadata_json,
            drug_column=drug_column,
            gene_column=gene_column,
            source_column=source_column,
            evidence_column=evidence_column,
            max_file_bytes=max_file_bytes,
            max_rows=max_rows,
            expected_sha256=artifact["raw_file_sha256"],
            openbio_version=PLUGIN_VERSION,
        )
        return io.NodeOutput(resource, report, code)


class OpenBioSingleCellDrugScores(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="Explicit drug-score expression source",
        include_raw=True,
        layer_default="log1p_norm",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDrugScores",
            display_name="DGIdb Single-Drug Target Score",
            category=CATEGORY,
            description="Compute one descriptive per-cell mean-expression score for one pinned DGIdb target set.",
            inputs=[
                AnnDataType.Input("adata"),
                DGIdbResourceType.Input("resource"),
                io.String.Input("drug", default=""),
                cls.EXPRESSION_SOURCE.input(),
                io.String.Input("output_key", default="drug_target_score"),
                io.Int.Input("min_matched_targets", default=1, min=1, max=2**31 - 1, advanced=True),
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
        adata: AnnData,
        resource: DGIdbResource,
        drug: str = "",
        source: DynamicExpressionSource | None = None,
        output_key: str = "drug_target_score",
        min_matched_targets: int = 1,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        runtime_parameters = {
            "drug": drug,
            "source_kind": expression.kind,
            "layer_name": expression.layer_name,
            "output_key": output_key,
            "min_matched_targets": min_matched_targets,
            "overwrite_existing": overwrite_existing,
            "openbio_version": PLUGIN_VERSION,
        }
        output, summary = run_drug_score(
            adata,
            resource,
            **runtime_parameters,
            pertpy_module=None,
        )
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        finish_adata(
            output,
            "dgidb_drug_score",
            parameters,
            int(adata.n_obs),
            int(adata.n_vars),
            started_at,
            warnings=warnings,
        )
        report = make_summary_result(
            summary=summary,
            title=f"DGIdb target-expression score: {drug}",
            operation="dgidb_drug_score",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        resource_table, resource_metadata, resource_accounting, resource_artifact = (
            validate_dgidb_resource(resource)
        )
        del resource_table
        code = drug_score_code(
            resource_metadata=resource_metadata,
            resource_accounting=resource_accounting,
            resource_artifact_metadata=resource_artifact,
            **runtime_parameters,
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellDrugHypergeometric(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDrugHypergeometric",
            display_name="Drug Hypergeometric Enrichment",
            category=CATEGORY,
            description="Test one explicit selected-gene set against every eligible pinned DGIdb drug target set.",
            inputs=[
                TableResultType.Input("table"),
                TableResultType.Input("universe"),
                DGIdbResourceType.Input("resource"),
                io.String.Input("comparison", default=""),
                io.String.Input("gene_column", default="gene", advanced=True),
                io.Int.Input("min_targets", default=3, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("min_overlap", default=2, min=1, max=2**31 - 1),
                io.Float.Input("max_p_adjusted", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Int.Input(
                    "max_output_rows",
                    default=100_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        table: TableResult,
        universe: TableResult,
        resource: DGIdbResource,
        comparison: str = "",
        gene_column: str = "gene",
        min_targets: int = 3,
        min_overlap: int = 2,
        max_p_adjusted: float = 0.05,
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        provenance = validate_enrichment_artifact_pair(
            table,
            universe,
            purpose="selected",
            selector=comparison,
            gene_column=gene_column,
            score_column=None,
            np=science.np,
            pd=science.pd,
        )
        runtime_parameters = {
            "artifact_family": provenance["artifact_family"],
            "comparison": provenance["comparison"],
            "comparison_column": provenance["comparison_column"],
            "gene_column": provenance["gene_column"],
            "evidence_scope": provenance["evidence_scope"],
            "inference_unit": provenance["inference_unit"],
            "direction": provenance["direction"],
            "replicate_aware": provenance["replicate_aware"],
            "technical_batch_handling": provenance["technical_batch_handling"],
            "upstream_parameters": provenance["upstream_parameters"],
            "provenance_warnings": provenance["provenance_warnings"],
            "provenance_declarations": provenance["provenance_declarations"],
            "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
            "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
            "expected_universe_fingerprint": provenance["universe_fingerprint"],
            "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
            "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
            "min_targets": min_targets,
            "min_overlap": min_overlap,
            "max_p_adjusted": max_p_adjusted,
            "max_output_rows": max_output_rows,
            "openbio_version": PLUGIN_VERSION,
        }
        evidence, summary = run_drug_ora(table, universe, resource, **runtime_parameters)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        result = make_table_result(
            title=f"Drug target overrepresentation: {provenance['comparison']}",
            operation="drug_hypergeometric",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
            table=evidence,
        )
        report = make_summary_result(
            summary=summary,
            title=f"Drug target overrepresentation summary: {provenance['comparison']}",
            operation="drug_hypergeometric",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
        )
        resource_table, resource_metadata, resource_accounting, resource_artifact = (
            validate_dgidb_resource(resource)
        )
        del resource_table
        code = drug_ora_code(
            resource_metadata=resource_metadata,
            resource_accounting=resource_accounting,
            resource_artifact_metadata=resource_artifact,
            **runtime_parameters,
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellDrugGSEA(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellDrugGSEA",
            display_name="Drug GSEA",
            category=CATEGORY,
            description="Test one complete pinned ranking against every eligible pinned DGIdb drug target set.",
            inputs=[
                TableResultType.Input("table"),
                TableResultType.Input("universe"),
                DGIdbResourceType.Input("resource"),
                io.String.Input("comparison", default=""),
                io.String.Input("gene_column", default="gene", advanced=True),
                io.String.Input("score_column", default="score", advanced=True),
                io.Int.Input("min_targets", default=15, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_targets", default=500, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("n_permutations", default=1000, min=2, max=2**31 - 1, advanced=True),
                io.Int.Input("random_seed", default=123, min=0, max=2**31 - 1, advanced=True),
                io.Int.Input(
                    "max_output_rows",
                    default=100_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        table: TableResult,
        universe: TableResult,
        resource: DGIdbResource,
        comparison: str = "",
        gene_column: str = "gene",
        score_column: str = "score",
        min_targets: int = 15,
        max_targets: int = 500,
        n_permutations: int = 1000,
        random_seed: int = 123,
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        provenance = validate_enrichment_artifact_pair(
            table,
            universe,
            purpose="ranked",
            selector=comparison,
            gene_column=gene_column,
            score_column=score_column,
            np=science.np,
            pd=science.pd,
        )
        runtime_parameters = {
            "artifact_family": provenance["artifact_family"],
            "comparison": provenance["comparison"],
            "comparison_column": provenance["comparison_column"],
            "gene_column": provenance["gene_column"],
            "score_column": provenance["score_column"],
            "evidence_scope": provenance["evidence_scope"],
            "inference_unit": provenance["inference_unit"],
            "direction": provenance["direction"],
            "replicate_aware": provenance["replicate_aware"],
            "technical_batch_handling": provenance["technical_batch_handling"],
            "upstream_parameters": provenance["upstream_parameters"],
            "provenance_warnings": provenance["provenance_warnings"],
            "provenance_declarations": provenance["provenance_declarations"],
            "expected_analysis_fingerprint": provenance["analysis_fingerprint"],
            "expected_ranking_fingerprint": provenance["ranking_fingerprint"],
            "expected_universe_fingerprint": provenance["universe_fingerprint"],
            "expected_table_content_fingerprint": provenance["table_content_fingerprint"],
            "expected_universe_content_fingerprint": provenance["universe_content_fingerprint"],
            "min_targets": min_targets,
            "max_targets": max_targets,
            "n_permutations": n_permutations,
            "random_seed": random_seed,
            "max_output_rows": max_output_rows,
            "openbio_version": PLUGIN_VERSION,
        }
        evidence, summary = run_drug_gsea(table, universe, resource, **runtime_parameters)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        result = make_table_result(
            title=f"Drug GSEA: {provenance['comparison']}",
            operation="drug_gsea",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
            random_seed=random_seed,
            table=evidence,
        )
        report = make_summary_result(
            summary=summary,
            title=f"Drug GSEA summary: {provenance['comparison']}",
            operation="drug_gsea",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
            random_seed=random_seed,
        )
        resource_table, resource_metadata, resource_accounting, resource_artifact = (
            validate_dgidb_resource(resource)
        )
        del resource_table
        code = drug_gsea_code(
            resource_metadata=resource_metadata,
            resource_accounting=resource_accounting,
            resource_artifact_metadata=resource_artifact,
            **runtime_parameters,
        )
        return io.NodeOutput(result, report, code)


ENRICHMENT_NODE_CLASSES = [
    OpenBioSingleCellAUCellScores,
    OpenBioSingleCellGSVAScores,
    OpenBioSingleCellGenePanelScores,
    OpenBioSingleCellPathwayScoreTTest,
    OpenBioSingleCellRankedGSEA,
    OpenBioSingleCellGeneSetOverrepresentation,
    OpenBioSingleCellDGIdbAnnotation,
    OpenBioSingleCellDrugScores,
    OpenBioSingleCellDrugHypergeometric,
    OpenBioSingleCellDrugGSEA,
]


__all__ = [
    "ENRICHMENT_NODE_CLASSES",
    "OpenBioSingleCellAUCellScores",
    "OpenBioSingleCellDGIdbAnnotation",
    "OpenBioSingleCellDrugGSEA",
    "OpenBioSingleCellDrugHypergeometric",
    "OpenBioSingleCellDrugScores",
    "OpenBioSingleCellGenePanelScores",
    "OpenBioSingleCellGeneSetOverrepresentation",
    "OpenBioSingleCellGSVAScores",
    "OpenBioSingleCellPathwayScoreTTest",
    "OpenBioSingleCellRankedGSEA",
]
