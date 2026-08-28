from __future__ import annotations

import hashlib
import time
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import PLUGIN_VERSION, dependencies
from .analysis_reporting import AnalysisReference, make_analysis_report
from .analysis_utils import finish_adata, make_summary_result, make_table_result
from .annotation_core import (
    build_celltypist_summary,
    celltypist_annotation_code,
    celltypist_model_fingerprint,
    map_cluster_annotations_code,
    run_celltypist_annotation,
    run_map_cluster_annotations,
)
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .files import input_file_fingerprint, resolve_input_path
from .marker_evidence import validate_marker_artifact_pair
from .node_types import AnnDataType, SummaryResultType, TableResultType
from .ora_evidence import build_marker_ora_summary, marker_ora_evidence_code, run_marker_ora_evidence

if TYPE_CHECKING:
    from .contracts import TableResult

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/annotation"

CELLTYPIST_METHOD_REFERENCE = AnalysisReference(
    citation=(
        "Domínguez Conde C, Xu C, Jarvis LB, et al. Cross-tissue immune cell analysis reveals "
        "tissue-specific features in humans. Science. 2022;376:eabl5197."
    ),
    doi="10.1126/science.abl5197",
    url="https://doi.org/10.1126/science.abl5197",
    kind="method",
)
CELLTYPIST_API_REFERENCE = AnalysisReference(
    citation="CellTypist developers. celltypist.annotate official API documentation.",
    url="https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html",
    kind="software_documentation",
)
CELLTYPIST_SOURCE_REFERENCE = AnalysisReference(
    citation="CellTypist developers. Official classifier and AnnotationResult implementation.",
    url="https://celltypist.readthedocs.io/en/stable/_modules/celltypist/classifier.html",
    kind="software_documentation",
)
SCANPY_REFERENCE = AnalysisReference(
    citation="Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.",
    doi="10.1186/s13059-017-1382-0",
    url="https://doi.org/10.1186/s13059-017-1382-0",
    kind="software",
)
ANNOTATION_PRACTICE_REFERENCE = AnalysisReference(
    citation=(
        "Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. "
        "Molecular Systems Biology. 2019;15:e8746."
    ),
    doi="10.15252/msb.20188746",
    url="https://doi.org/10.15252/msb.20188746",
    kind="practice",
)
PANDAS_MAPPING_REFERENCE = AnalysisReference(
    citation="pandas developers. pandas.Series.map and categorical data documentation.",
    url="https://pandas.pydata.org/docs/reference/api/pandas.Series.map.html",
    kind="software_documentation",
)
PANDAS_REFERENCE = AnalysisReference(
    citation="McKinney W. Data Structures for Statistical Computing in Python. SciPy 2010.",
    doi="10.25080/Majora-92bf1922-00a",
    url="https://doi.org/10.25080/Majora-92bf1922-00a",
    kind="software",
)
ANNDATA_REFERENCE = AnalysisReference(
    citation="Virshup I, et al. anndata: Annotated data. Journal of Open Source Software. 2024;9:4371.",
    doi="10.21105/joss.04371",
    url="https://doi.org/10.21105/joss.04371",
    kind="software",
)
JSON_REFERENCE = AnalysisReference(
    citation="Bray T. The JavaScript Object Notation (JSON) Data Interchange Format. RFC 8259.",
    doi="10.17487/RFC8259",
    url="https://www.rfc-editor.org/info/rfc8259/",
    kind="software_documentation",
)
DECOUPLER_REFERENCE = AnalysisReference(
    citation=(
        "Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer "
        "biological activities from omics data. Bioinformatics Advances. 2022;2:vbac016."
    ),
    doi="10.1093/bioadv/vbac016",
    url="https://doi.org/10.1093/bioadv/vbac016",
    kind="software",
)
DECOUPLER_QUERY_SET_REFERENCE = AnalysisReference(
    citation="decoupler developers. decoupler.mt.query_set official API and implementation.",
    url="https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.query_set.html",
    kind="software_documentation",
)
FISHER_REFERENCE = AnalysisReference(
    citation=(
        "Fisher RA. On the interpretation of chi-square from contingency tables, and the "
        "calculation of P. Journal of the Royal Statistical Society. 1922;85:87-94."
    ),
    doi="10.2307/2340521",
    url="https://doi.org/10.2307/2340521",
    kind="method",
)
HALDANE_REFERENCE = AnalysisReference(
    citation=(
        "Haldane JBS. The estimation and significance of the logarithm of a ratio of frequencies. "
        "Annals of Human Genetics. 1956;20:309-311."
    ),
    doi="10.1111/j.1469-1809.1955.tb01285.x",
    url="https://doi.org/10.1111/j.1469-1809.1955.tb01285.x",
    kind="method",
)
BH_REFERENCE = AnalysisReference(
    citation=(
        "Benjamini Y, Hochberg Y. Controlling the false discovery rate. Journal of the Royal "
        "Statistical Society Series B. 1995;57:289-300."
    ),
    doi="10.1111/j.2517-6161.1995.tb02031.x",
    url="https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
    kind="method",
)
ENRICHMENT_BIAS_REFERENCE = AnalysisReference(
    citation=(
        "Timmons JA, Szkop KJ, Gallagher IJ. Multiple sources of bias confound functional "
        "enrichment analysis. Genome Biology. 2015;16:186."
    ),
    doi="10.1186/s13059-015-0761-7",
    url="https://doi.org/10.1186/s13059-015-0761-7",
    kind="practice",
)


class OpenBioSingleCellCellTypistAnnotation(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="CellTypist expression source",
        default="X",
        include_raw=True,
        layer_input_id="layer_name",
        layer_default="counts",
    )

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCellTypistAnnotation",
            display_name="CellTypist Annotation",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                cls.EXPRESSION_SOURCE.input(),
                io.Combo.Input(
                    "expression_state",
                    options=["verified_cp10k_log1p", "verified_counts"],
                    default="verified_cp10k_log1p",
                ),
                io.String.Input("model", default=""),
                io.Boolean.Input("majority_voting", default=False),
                io.String.Input("over_clustering_key", default="", advanced=True),
                io.Float.Input("min_prop", default=0.0, min=0.0, max=1.0, step=0.05, advanced=True),
                io.String.Input("label_column", default="celltypist_cell_type", advanced=True),
                io.String.Input("confidence_column", default="celltypist_confidence", advanced=True),
                io.String.Input("probability_key", default="celltypist_probabilities", advanced=True),
                io.Boolean.Input("store_decision_matrix", default=False, advanced=True),
                io.String.Input("decision_key", default="celltypist_decision_scores", advanced=True),
                io.String.Input("metadata_key", default="celltypist", advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def fingerprint_inputs(cls, model: str = "", **kwargs: Any) -> Any:
        return celltypist_model_fingerprint(model)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        source: DynamicExpressionSource | None = None,
        expression_state: str = "verified_cp10k_log1p",
        model: str = "",
        majority_voting: bool = False,
        over_clustering_key: str = "",
        min_prop: float = 0.0,
        label_column: str = "celltypist_cell_type",
        confidence_column: str = "celltypist_confidence",
        probability_key: str = "celltypist_probabilities",
        store_decision_matrix: bool = False,
        decision_key: str = "celltypist_decision_scores",
        metadata_key: str = "celltypist",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output, diagnostics = run_celltypist_annotation(
            adata,
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            expression_state=expression_state,
            model_identifier=model,
            majority_voting=majority_voting,
            over_clustering_key=over_clustering_key,
            min_prop=min_prop,
            label_column=label_column,
            confidence_column=confidence_column,
            probability_key=probability_key,
            store_decision_matrix=store_decision_matrix,
            decision_key=decision_key,
            metadata_key=metadata_key,
            overwrite_existing=overwrite_existing,
            openbio_version=PLUGIN_VERSION,
            celltypist_module=None,
        )
        parameters = dict(diagnostics["parameters"])
        warnings = list(diagnostics["warnings"])
        finish_adata(
            output,
            "celltypist_annotation",
            parameters,
            cells,
            genes,
            started_at,
            warnings=warnings,
        )
        code = celltypist_annotation_code(
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            expression_state=expression_state,
            resolved_model_path=diagnostics["model"]["resolved_path"],
            requested_model_identifier=diagnostics["model"]["requested_identifier"],
            expected_model_sha256=diagnostics["model"]["sha256"],
            majority_voting=majority_voting,
            over_clustering_key=over_clustering_key.strip(),
            min_prop=float(min_prop),
            label_column=label_column.strip(),
            confidence_column=confidence_column.strip(),
            probability_key=probability_key.strip(),
            store_decision_matrix=store_decision_matrix,
            decision_key=decision_key.strip(),
            metadata_key=metadata_key.strip(),
            overwrite_existing=overwrite_existing,
            openbio_version=PLUGIN_VERSION,
        )
        summary = build_celltypist_summary(diagnostics)
        report = make_summary_result(
            summary=summary,
            title="CellTypist provisional annotation summary",
            operation="celltypist_annotation",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
        )
        return io.NodeOutput(output, report, code)


class OpenBioSingleCellMarkerORAEvidence(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMarkerORAEvidence",
            display_name="Marker ORA Evidence",
            category=CATEGORY,
            inputs=[
                TableResultType.Input("table"),
                TableResultType.Input("universe"),
                io.String.Input("resource_csv", default=""),
                io.String.Input("resource_metadata_json", default="{}"),
                io.String.Input("source_column", default="source", advanced=True),
                io.String.Input("target_column", default="target", advanced=True),
                io.Int.Input("min_targets", default=3, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("min_overlap", default=2, min=1, max=2**31 - 1, advanced=True),
                io.Float.Input(
                    "max_p_adjusted",
                    default=0.05,
                    min=0.0,
                    max=1.0,
                    step=0.01,
                ),
            ],
            outputs=[
                TableResultType.Output(display_name="table"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, resource_csv: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(resource_csv, extensions=(".csv",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, resource_csv: str, **kwargs: Any) -> Any:
        path = resolve_input_path(resource_csv, extensions=(".csv",))
        identity = input_file_fingerprint(resource_csv, (".csv",))
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return ("openbio-marker-ora-resource-v1", *identity, digest.hexdigest())

    @classmethod
    def execute(
        cls,
        table: TableResult,
        universe: TableResult,
        resource_csv: str = "",
        resource_metadata_json: str = "{}",
        source_column: str = "source",
        target_column: str = "target",
        min_targets: int = 3,
        min_overlap: int = 2,
        max_p_adjusted: float = 0.05,
    ) -> io.NodeOutput:
        science = dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        provenance = validate_marker_artifact_pair(
            table,
            universe,
            np=science.np,
            pd=science.pd,
            allowed_table_operations=("filter_marker_genes",),
        )
        if provenance["table_operation"] != "filter_marker_genes":
            raise RuntimeError("Marker ORA Evidence accepted a non-filtered marker artifact unexpectedly.")
        resolved_resource = resolve_input_path(resource_csv, extensions=(".csv",))
        evidence, diagnostics = run_marker_ora_evidence(
            table,
            universe,
            expected_analysis_fingerprint=provenance["analysis_fingerprint"],
            expected_ranking_fingerprint=provenance["ranking_fingerprint"],
            expected_universe_fingerprint=provenance["universe_fingerprint"],
            expected_marker_content_fingerprint=provenance["table_content_fingerprint"],
            expected_universe_content_fingerprint=provenance["universe_content_fingerprint"],
            group_labels=provenance["group_labels"],
            upstream_parameters=provenance["upstream_parameters"],
            resource_path=resolved_resource,
            requested_resource_path=resource_csv.strip(),
            resource_metadata_json=resource_metadata_json,
            source_column=source_column,
            target_column=target_column,
            min_targets=min_targets,
            min_overlap=min_overlap,
            max_p_adjusted=max_p_adjusted,
            openbio_version=PLUGIN_VERSION,
            decoupler_module=None,
        )
        summary = build_marker_ora_summary(diagnostics)
        parameters = dict(summary["parameters"])
        warnings = list(summary["warnings"])
        result = make_table_result(
            table=evidence,
            title="Marker ORA evidence",
            operation="marker_ora_evidence",
            parameters=parameters,
            description=(
                "Explicit selected-set over-representation evidence for Cluster marker review; no AnnData labels were assigned."
            ),
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
        )
        code = marker_ora_evidence_code(
            expected_analysis_fingerprint=provenance["analysis_fingerprint"],
            expected_ranking_fingerprint=provenance["ranking_fingerprint"],
            expected_universe_fingerprint=provenance["universe_fingerprint"],
            expected_marker_content_fingerprint=provenance["table_content_fingerprint"],
            expected_universe_content_fingerprint=provenance["universe_content_fingerprint"],
            group_labels=provenance["group_labels"],
            upstream_parameters=provenance["upstream_parameters"],
            resource_path=resolved_resource,
            requested_resource_path=resource_csv.strip(),
            resource_metadata_json=resource_metadata_json,
            expected_resource_sha256=diagnostics["resource"]["sha256"],
            source_column=source_column.strip(),
            target_column=target_column.strip(),
            min_targets=min_targets,
            min_overlap=min_overlap,
            max_p_adjusted=float(max_p_adjusted),
            openbio_version=PLUGIN_VERSION,
        )
        report = make_summary_result(
            summary=summary,
            title="Marker ORA evidence summary",
            operation="marker_ora_evidence",
            parameters=parameters,
            description=summary["results"],
            warnings=warnings,
            input_cells=table.input_cells,
            input_genes=table.input_genes,
            started_at=started_at,
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellMapClusterAnnotations(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellMapClusterAnnotations",
            display_name="Map Cluster Annotations",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("groupby", default="leiden"),
                io.String.Input("mapping_json", default="{}"),
                io.String.Input("output_column", default="cell_type", advanced=True),
                io.Combo.Input(
                    "unmapped_policy",
                    options=["error", "preserve_cluster_label", "set_missing"],
                    default="error",
                ),
                io.Combo.Input(
                    "annotation_status",
                    options=["provisional", "curated"],
                    default="provisional",
                ),
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
        groupby: str = "leiden",
        mapping_json: str = "{}",
        output_column: str = "cell_type",
        unmapped_policy: str = "error",
        annotation_status: str = "provisional",
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        dependencies.require_scientific_dependencies()
        started_at = time.perf_counter()
        cells, genes = int(adata.n_obs), int(adata.n_vars)
        output, diagnostics = run_map_cluster_annotations(
            adata,
            groupby=groupby,
            mapping_json=mapping_json,
            output_column=output_column,
            unmapped_policy=unmapped_policy,
            annotation_status=annotation_status,
            overwrite_existing=overwrite_existing,
        )
        normalized_mapping = {
            **diagnostics["effective_mapping"],
            **diagnostics["unused_declared_level_mapping"],
        }
        parameters = {
            "groupby": groupby.strip(),
            "mapping": normalized_mapping,
            "output_column": output_column.strip(),
            "unmapped_policy": unmapped_policy,
            "annotation_status": annotation_status,
            "overwrite_existing": overwrite_existing,
        }
        warnings = []
        if diagnostics["missing_source_cells"]:
            warnings.append(
                f"{diagnostics['missing_source_cells']:,} cells had missing source labels; true missing values were preserved."
            )
        if unmapped_policy == "preserve_cluster_label" and diagnostics["unmapped_cells"]:
            warnings.append("The output mixes mapped biological labels with preserved raw cluster identifiers.")
        if annotation_status == "curated":
            warnings.append(
                "Curated status is a caller declaration; the software did not verify expert review, marker evidence, or ontology validity."
            )
        finish_adata(
            output,
            "map_cluster_annotations",
            parameters,
            cells,
            genes,
            started_at,
            warnings=warnings,
        )
        code = map_cluster_annotations_code(
            groupby=groupby.strip(),
            mapping_json=mapping_json,
            output_column=output_column.strip(),
            unmapped_policy=unmapped_policy,
            annotation_status=annotation_status,
            overwrite_existing=overwrite_existing,
        )
        mapped_clusters = len(diagnostics["effective_mapping"])
        output_counts_text = ", ".join(
            f"{category!r}={count:,}" for category, count in diagnostics["output_counts"].items()
        ) or "none"
        many_to_one_count = len(diagnostics["many_to_one_merges"])
        report, code = make_analysis_report(
            node_id="OpenBioSingleCellMapClusterAnnotations",
            title="Cluster annotation mapping summary",
            operation="map_cluster_annotations",
            methods=(
                f"A strict JSON map was applied to observed {groupby.strip()!r} levels. "
                f"Levels absent from the map were handled by policy {unmapped_policy!r}; output categories followed the declared "
                "categorical source order or lexical observed-level order."
            ),
            results=(
                f"Mapped {mapped_clusters:,} of {diagnostics['observed_level_count']:,} observed clusters, covering "
                f"{diagnostics['mapped_cells']:,} of {cells:,} cells. Policy {unmapped_policy!r} preserved "
                f"{diagnostics['preserved_cells']:,} cells and set {diagnostics['set_missing_cells']:,} cells missing; "
                f"{diagnostics['missing_source_cells']:,} cells were already missing in the source and "
                f"{diagnostics['output_missing_cells']:,} cells are missing in the output. The output has "
                f"{len(diagnostics['output_categories']):,} categories with counts {output_counts_text}, and the map "
                f"contains {many_to_one_count:,} many-to-one target merges. The output column existed before execution: "
                f"{'yes' if diagnostics['output_existed'] else 'no'}; overwrite_existing was "
                f"{'true' if overwrite_existing else 'false'}. "
                f"The result is recorded as {'caller-declared Curated annotation' if annotation_status == 'curated' else 'Provisional annotation'}."
            ),
            key_results={
                **diagnostics,
                "annotation_provenance": output.uns["openbio_singlecell"]["annotations"][output_column.strip()],
                "mapping_unknown_keys": [],
                "output_dtype": "category",
                "output_ordered": False,
            },
            parameters=parameters,
            references=(
                ANNOTATION_PRACTICE_REFERENCE,
                PANDAS_MAPPING_REFERENCE,
                PANDAS_REFERENCE,
                ANNDATA_REFERENCE,
                JSON_REFERENCE,
            ),
            software_packages=("pandas", "anndata"),
            warnings=warnings,
            limitations=(
                "Every cell in one source cluster receives the same label, so within-cluster heterogeneity is hidden.",
                "Cluster resolution and preprocessing affect the source partition and therefore the applied annotation.",
                "The mapping operation does not inspect Cluster marker evidence, reference expression, or ontology identifiers.",
                "Provisional annotation requires review before use as Curated annotation in Sample-level Condition contrast.",
            ),
            input_cells=cells,
            input_genes=genes,
            started_at=started_at,
            code=code,
        )
        return io.NodeOutput(output, report, code)


ANNOTATION_NODE_CLASSES = [
    OpenBioSingleCellCellTypistAnnotation,
    OpenBioSingleCellMarkerORAEvidence,
    OpenBioSingleCellMapClusterAnnotations,
]


__all__ = [
    "ANNOTATION_NODE_CLASSES",
    "OpenBioSingleCellCellTypistAnnotation",
    "OpenBioSingleCellMapClusterAnnotations",
    "OpenBioSingleCellMarkerORAEvidence",
]
