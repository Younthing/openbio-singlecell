from __future__ import annotations

import hashlib
from typing import Any

from comfy_api.latest import io

from .dgidb_resource import DGIDB_RESOURCE_EXTENSIONS, dgidb_resource_cache_fingerprint
from .expression_source import _AUCELL_SPEC, _DRUG_SCORE_SPEC, _GENE_PANEL_SPEC, _GSVA_SPEC
from .files import input_file_fingerprint, resolve_input_path
from .gene_set_scoring import gene_set_resource_cache_fingerprint
from .node_types import AnnDataType, DGIdbResourceType, TableResultType, analysis_outputs

CATEGORY = "openbio/single-cell/enrichment"
GENE_SET_EXTENSIONS = (".csv", ".tsv", ".gmt")


def _validate_gene_set_file(path: str) -> bool | str:
    try:
        resolve_input_path(path, extensions=GENE_SET_EXTENSIONS)
    except (ValueError, FileNotFoundError, OSError) as exc:
        return str(exc)
    return True


class OpenBioSingleCellAUCellScores(io.ComfyNode):
    EXPRESSION_SOURCE = _AUCELL_SPEC

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
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        return _validate_gene_set_file(gene_sets_file)

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        return gene_set_resource_cache_fingerprint(path, identity)


class OpenBioSingleCellGSVAScores(io.ComfyNode):
    EXPRESSION_SOURCE = _GSVA_SPEC

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
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        return _validate_gene_set_file(gene_sets_file)

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        return gene_set_resource_cache_fingerprint(path, identity)


class OpenBioSingleCellGenePanelScores(io.ComfyNode):
    EXPRESSION_SOURCE = _GENE_PANEL_SPEC

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
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        return _validate_gene_set_file(gene_sets_file)

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        return gene_set_resource_cache_fingerprint(path, identity)


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
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


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
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        return _validate_gene_set_file(gene_sets_file)

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        digest = hashlib.sha256()
        with open(path, "rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        return ("openbio-ranked-gsea-resource-v1", *identity, digest.hexdigest())


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
                io.Float.Input("max_p_adjusted", default=0.05, min=0.0, max=1.0, step=0.01),
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )

    @classmethod
    def validate_inputs(cls, gene_sets_file: str, **kwargs: Any) -> bool | str:
        return _validate_gene_set_file(gene_sets_file)

    @classmethod
    def fingerprint_inputs(cls, gene_sets_file: str, **kwargs: Any) -> Any:
        path = resolve_input_path(gene_sets_file, extensions=GENE_SET_EXTENSIONS)
        identity = input_file_fingerprint(gene_sets_file, GENE_SET_EXTENSIONS)
        digest = hashlib.sha256()
        with open(path, "rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        return ("openbio-generic-ora-resource-v1", *identity, digest.hexdigest())


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
                io.Int.Input("max_file_bytes", default=536_870_912, min=1, max=2**63 - 1, advanced=True),
                io.Int.Input("max_rows", default=2_000_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(DGIdbResourceType.Output(display_name="resource")),
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


class OpenBioSingleCellDrugScores(io.ComfyNode):
    EXPRESSION_SOURCE = _DRUG_SCORE_SPEC

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
            outputs=analysis_outputs(AnnDataType.Output(display_name="adata")),
        )


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
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


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
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


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
