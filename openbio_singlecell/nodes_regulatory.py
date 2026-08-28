from __future__ import annotations

import time
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import PLUGIN_VERSION
from .analysis_utils import make_summary_result, make_table_result
from .collectri_ulm import (
    collectri_resource_cache_fingerprint,
    collectri_ulm_code,
    run_collectri_ulm,
)
from .expression_source import (
    DynamicExpressionSource,
    ExpressionSourceSpec,
)
from .files import input_file_fingerprint, resolve_input_path
from .node_types import (
    AnnDataType,
    SCENICBinaryArtifactType,
    SCENICResultArtifactType,
    SummaryResultType,
    TableResultType,
    TFActivityArtifactType,
)
from .pyscenic_import import (
    import_pyscenic_bundle,
    pyscenic_bundle_cache_fingerprint,
    pyscenic_import_code,
)
from .scenic_artifact import SCENICResultArtifact
from .scenic_binarization import binarize_scenic_activity, scenic_binarization_code
from .scenic_membership import scenic_membership_code, scenic_regulon_membership
from .scenic_rss import compute_scenic_rss, scenic_rss_code
from .tf_activity_artifact import TFActivityArtifact
from .tf_activity_ranking import TF_RANKING_METHODS, rank_tf_activities, rank_tf_activities_code

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/regulatory"


class OpenBioSingleCellCollecTRIULM(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="CollecTRI normalized expression source",
        default="layer",
        include_raw=False,
        layer_default="log1p_norm",
    )

    @staticmethod
    def _resolve_resource(resource: Mapping[str, object] | None) -> dict[str, object]:
        if resource is None:
            resource = {
                "resource": "local_network",
                "network_csv": "",
                "resource_metadata_json": "{}",
            }
        if not isinstance(resource, Mapping):
            raise TypeError("CollecTRI resource must be a DynamicCombo value.")
        mode = resource.get("resource")
        if mode == "local_network":
            expected = {"resource", "network_csv", "resource_metadata_json"}
            if set(resource) != expected:
                raise ValueError(
                    "CollecTRI local resource fields are invalid; "
                    f"missing={sorted(expected - set(resource))}, unknown={sorted(set(resource) - expected)}."
                )
            network_csv = resource.get("network_csv")
            metadata_json = resource.get("resource_metadata_json")
            if not isinstance(network_csv, str):
                raise TypeError("CollecTRI network_csv must be a string.")
            if not isinstance(metadata_json, str):
                raise TypeError("CollecTRI resource_metadata_json must be a string.")
            return {
                "resource_mode": "local_network",
                "network_csv": network_csv,
                "resource_metadata_json": metadata_json,
                "affiliation_license": "academic",
                "allow_network_access": False,
            }
        if mode == "official_collectri":
            expected = {"resource", "affiliation_license", "allow_network_access"}
            if set(resource) != expected:
                raise ValueError(
                    "CollecTRI official resource fields are invalid; "
                    f"missing={sorted(expected - set(resource))}, unknown={sorted(set(resource) - expected)}."
                )
            affiliation_license = resource.get("affiliation_license")
            allow_network_access = resource.get("allow_network_access")
            if not isinstance(affiliation_license, str):
                raise TypeError("CollecTRI affiliation_license must be a string.")
            if not isinstance(allow_network_access, bool):
                raise TypeError("CollecTRI allow_network_access must be boolean.")
            return {
                "resource_mode": "official_collectri",
                "network_csv": None,
                "resource_metadata_json": "{}",
                "affiliation_license": affiliation_license,
                "allow_network_access": allow_network_access,
            }
        raise ValueError(f"Unsupported CollecTRI resource mode: {mode!r}.")

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellCollecTRIULM",
            display_name="CollecTRI ULM Activities",
            category=CATEGORY,
            description=(
                "Infer exploratory observation-level TF activities from normalized expression and one explicitly "
                "licensed, fingerprinted signed CollecTRI network."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.DynamicCombo.Input(
                    "resource",
                    options=[
                        io.DynamicCombo.Option(
                            "local_network",
                            [
                                io.String.Input("network_csv", default=""),
                                io.String.Input("resource_metadata_json", default="{}", multiline=True),
                            ],
                        ),
                        io.DynamicCombo.Option(
                            "official_collectri",
                            [
                                io.Combo.Input(
                                    "affiliation_license",
                                    options=["academic", "commercial", "nonprofit"],
                                    default="academic",
                                ),
                                io.Boolean.Input("allow_network_access", default=False),
                            ],
                        ),
                    ],
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.Combo.Input("complex_policy", options=["retain", "remove"], default="retain"),
                io.Int.Input("min_targets", default=5, min=1, max=2**31 - 1),
                io.Int.Input("batch_size", default=250_000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_output_rows", default=2_000_000, min=1, max=2**31 - 1, advanced=True),
                io.Float.Input(
                    "max_working_memory_gib", default=4.0, min=0.001, max=1024.0, advanced=True
                ),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                TFActivityArtifactType.Output(display_name="activities"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, resource: Mapping[str, object] | None = None, **kwargs: Any) -> bool | str:
        try:
            resolved = cls._resolve_resource(resource)
            if resolved["resource_mode"] == "local_network":
                resolve_input_path(str(resolved["network_csv"]), extensions=(".csv",))
        except (TypeError, ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, resource: Mapping[str, object] | None = None, **kwargs: Any) -> Any:
        resolved = cls._resolve_resource(resource)
        if resolved["resource_mode"] == "official_collectri":
            return (
                "openbio-official-collectri-v2",
                resolved["affiliation_license"],
                resolved["allow_network_access"],
            )
        requested_path = str(resolved["network_csv"])
        path = resolve_input_path(requested_path, extensions=(".csv",))
        identity = input_file_fingerprint(requested_path, (".csv",))
        return collectri_resource_cache_fingerprint(path, identity)

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        resource: Mapping[str, object] | None = None,
        source: DynamicExpressionSource | None = None,
        complex_policy: str = "retain",
        min_targets: int = 5,
        batch_size: int = 250_000,
        max_output_rows: int = 2_000_000,
        max_working_memory_gib: float = 4.0,
        overwrite_existing: bool = False,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        resolved_resource = cls._resolve_resource(resource)
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        requested_path = resolved_resource["network_csv"]
        network_path = (
            resolve_input_path(str(requested_path), extensions=(".csv",))
            if resolved_resource["resource_mode"] == "local_network"
            else None
        )
        output, activities, summary = run_collectri_ulm(
            adata,
            resource_mode=str(resolved_resource["resource_mode"]),
            network_path=network_path,
            resource_metadata_json=str(resolved_resource["resource_metadata_json"]),
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            affiliation_license=str(resolved_resource["affiliation_license"]),
            allow_network_access=bool(resolved_resource["allow_network_access"]),
            complex_policy=complex_policy,
            min_targets=min_targets,
            batch_size=batch_size,
            max_output_rows=max_output_rows,
            max_working_memory_gib=max_working_memory_gib,
            overwrite_existing=overwrite_existing,
            openbio_version=PLUGIN_VERSION,
        )
        report = make_summary_result(
            summary=summary,
            title="CollecTRI ULM activity summary",
            operation="collectri_ulm",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        code = collectri_ulm_code(
            parameters={
                "resource_mode": str(resolved_resource["resource_mode"]),
                "network_path": network_path,
                "resource_metadata_json": str(resolved_resource["resource_metadata_json"]),
                "source_kind": expression.kind,
                "layer_name": expression.layer_name,
                "affiliation_license": str(resolved_resource["affiliation_license"]),
                "allow_network_access": bool(resolved_resource["allow_network_access"]),
                "complex_policy": complex_policy,
                "min_targets": min_targets,
                "batch_size": batch_size,
                "max_output_rows": max_output_rows,
                "max_working_memory_gib": max_working_memory_gib,
                "overwrite_existing": overwrite_existing,
                "expected_resource_sha256": summary["key_results"]["resource"]["file_sha256"],
                "expected_canonical_network_sha256": summary["key_results"]["resource"][
                    "canonical_network_sha256"
                ],
                "openbio_version": PLUGIN_VERSION,
            },
            resolved_resource_provenance=summary["key_results"]["resource"],
        )
        return io.NodeOutput(output, activities, report, code)


class OpenBioSingleCellRankTFActivities(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellRankTFActivities",
            display_name="Rank TF Activities",
            category=CATEGORY,
            description=(
                "Characterize a validated cell-level TF activity artifact across annotations without rerunning "
                "activity inference. This is exploratory annotation evidence, not Condition inference."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                TFActivityArtifactType.Input("activities"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                ),
                io.String.Input("reference", default="rest"),
                io.Combo.Input(
                    "method", options=list(TF_RANKING_METHODS), default="t-test_overestim_var"
                ),
                io.Float.Input(
                    "report_p_adjusted", default=0.05, min=0.0, max=1.0, step=0.01, advanced=True
                ),
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
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
        activities: TFActivityArtifact,
        annotation_key: str = "cell_type",
        annotation_status: str = "unknown",
        reference: str = "rest",
        method: str = "t-test_overestim_var",
        report_p_adjusted: float = 0.05,
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, summary = rank_tf_activities(
            adata,
            activities,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            reference=reference,
            method=method,
            report_p_adjusted=report_p_adjusted,
            max_output_rows=max_output_rows,
            openbio_version=PLUGIN_VERSION,
        )
        result = make_table_result(
            title=f"TF activities by {annotation_key}",
            operation="rank_tf_activities",
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
            title=f"TF activity ranking by {annotation_key}",
            operation="rank_tf_activities",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        code = rank_tf_activities_code(
            parameters={
                "annotation_key": annotation_key,
                "annotation_status": annotation_status,
                "reference": reference,
                "method": method,
                "report_p_adjusted": report_p_adjusted,
                "max_output_rows": max_output_rows,
                "openbio_version": PLUGIN_VERSION,
            }
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellImportPySCENICResults(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellImportPySCENICResults",
            display_name="Import pySCENIC Results",
            category=CATEGORY,
            description=(
                "Safely import one complete hash-bound external pySCENIC 0.12.1 CSV bundle without loading "
                "pySCENIC, evaluating code, downloading resources, or executing a shell command."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("run_manifest_json", default=""),
                io.Boolean.Input("overwrite", default=False, advanced=True),
                io.Int.Input(
                    "max_file_bytes", default=2_147_483_647, min=1, max=2**63 - 1, advanced=True
                ),
                io.Int.Input(
                    "max_adjacency_edges", default=10_000_000, min=1, max=2**31 - 1, advanced=True
                ),
                io.Int.Input(
                    "max_regulon_edges", default=2_000_000, min=1, max=2**31 - 1, advanced=True
                ),
                io.Int.Input(
                    "max_dense_bytes", default=1_073_741_824, min=1, max=2**63 - 1, advanced=True
                ),
            ],
            outputs=[
                AnnDataType.Output(display_name="adata"),
                SCENICResultArtifactType.Output(display_name="scenic_result"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, run_manifest_json: str, **kwargs: Any) -> bool | str:
        try:
            resolve_input_path(run_manifest_json, extensions=(".json",))
        except (ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(
        cls,
        run_manifest_json: str,
        max_file_bytes: int = 2_147_483_647,
        **kwargs: Any,
    ) -> Any:
        path = resolve_input_path(run_manifest_json, extensions=(".json",))
        return pyscenic_bundle_cache_fingerprint(
            path,
            max_file_bytes=max_file_bytes,
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        run_manifest_json: str = "",
        overwrite: bool = False,
        max_file_bytes: int = 2_147_483_647,
        max_adjacency_edges: int = 10_000_000,
        max_regulon_edges: int = 2_000_000,
        max_dense_bytes: int = 1_073_741_824,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        manifest_path = resolve_input_path(run_manifest_json, extensions=(".json",))
        output, scenic_result, summary = import_pyscenic_bundle(
            adata,
            manifest_path,
            overwrite=overwrite,
            max_file_bytes=max_file_bytes,
            max_adjacency_edges=max_adjacency_edges,
            max_regulon_edges=max_regulon_edges,
            max_dense_bytes=max_dense_bytes,
            openbio_version=PLUGIN_VERSION,
        )
        report = make_summary_result(
            summary=summary,
            title="Imported pySCENIC regulatory evidence",
            operation="import_pyscenic_results",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        code = pyscenic_import_code(
            manifest_path=manifest_path,
            overwrite=overwrite,
            max_file_bytes=max_file_bytes,
            max_adjacency_edges=max_adjacency_edges,
            max_regulon_edges=max_regulon_edges,
            max_dense_bytes=max_dense_bytes,
        )
        return io.NodeOutput(output, scenic_result, report, code)


class OpenBioSingleCellSCENICRegulonSpecificity(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCENICRegulonSpecificity",
            display_name="SCENIC Regulon Specificity",
            category=CATEGORY,
            description=(
                "Compute descriptive pySCENIC 0.12.1 regulon-specificity scores from a validated immutable "
                "SCENIC artifact; this is annotation evidence, not a Condition test."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                SCENICResultArtifactType.Input("scenic_result"),
                io.String.Input("annotation_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                ),
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
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
        scenic_result: SCENICResultArtifact,
        annotation_key: str = "cell_type",
        annotation_status: str = "unknown",
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, summary = compute_scenic_rss(
            adata,
            scenic_result,
            annotation_key=annotation_key,
            annotation_status=annotation_status,
            max_output_rows=max_output_rows,
            openbio_version=PLUGIN_VERSION,
        )
        result = make_table_result(
            title=f"SCENIC regulon specificity by {annotation_key}",
            operation="scenic_regulon_specificity",
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
            title=f"SCENIC regulon specificity by {annotation_key}",
            operation="scenic_regulon_specificity",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
        )
        code = scenic_rss_code(
            parameters={
                "annotation_key": annotation_key,
                "annotation_status": annotation_status,
                "max_output_rows": max_output_rows,
                "openbio_version": PLUGIN_VERSION,
            }
        )
        return io.NodeOutput(result, report, code)


class OpenBioSingleCellSCENICActivityBinarization(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCENICActivityBinarization",
            display_name="SCENIC Activity Binarization",
            category=CATEGORY,
            description=(
                "Apply the deterministic one-process pySCENIC 0.12.1 HDT-compatible heuristic to one "
                "validated activity artifact and return immutable binary evidence plus complete thresholds."
            ),
            inputs=[
                SCENICResultArtifactType.Input("scenic_result"),
                io.Int.Input("random_seed", default=1, min=1, max=2**31 - 1),
                TableResultType.Input("threshold_overrides", optional=True),
                io.Int.Input(
                    "max_dense_bytes", default=1_073_741_824, min=1, max=2**63 - 1, advanced=True
                ),
            ],
            outputs=[
                SCENICBinaryArtifactType.Output(display_name="binary"),
                TableResultType.Output(display_name="thresholds"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def execute(
        cls,
        scenic_result: SCENICResultArtifact,
        random_seed: int = 1,
        threshold_overrides: Any | None = None,
        max_dense_bytes: int = 1_073_741_824,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        overrides_frame = (
            threshold_overrides.table if hasattr(threshold_overrides, "table") else threshold_overrides
        )
        binary, threshold_table, summary = binarize_scenic_activity(
            scenic_result,
            random_seed=random_seed,
            threshold_overrides=overrides_frame,
            max_dense_bytes=max_dense_bytes,
            openbio_version=PLUGIN_VERSION,
        )
        provenance = scenic_result.provenance
        dimensions = provenance["input_dimensions"]
        threshold_result = make_table_result(
            title="SCENIC activity thresholds",
            operation="scenic_activity_binarization",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(dimensions["cells"]),
            input_genes=int(dimensions["genes"]),
            started_at=started_at,
            random_seed=random_seed,
            table=threshold_table,
        )
        report = make_summary_result(
            summary=summary,
            title="SCENIC activity binarization",
            operation="scenic_activity_binarization",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(dimensions["cells"]),
            input_genes=int(dimensions["genes"]),
            started_at=started_at,
            random_seed=random_seed,
        )
        code = scenic_binarization_code(
            parameters={
                "random_seed": random_seed,
                "max_dense_bytes": max_dense_bytes,
                "openbio_version": PLUGIN_VERSION,
            }
        )
        return io.NodeOutput(binary, threshold_result, report, code)


class OpenBioSingleCellSCENICTFModules(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCENICTFModules",
            display_name="SCENIC Final Regulon Membership",
            category=CATEGORY,
            description=(
                "Project the final motif-pruned regulon membership from a validated immutable SCENIC result; "
                "this node does not rerun pre-cisTarget adjacency modules."
            ),
            inputs=[
                SCENICResultArtifactType.Input("scenic_result"),
                io.String.Input("transcription_factor", default=""),
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
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
        scenic_result: SCENICResultArtifact,
        transcription_factor: str = "",
        max_output_rows: int = 100_000,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        table, summary = scenic_regulon_membership(
            scenic_result,
            transcription_factor=transcription_factor,
            max_output_rows=max_output_rows,
            openbio_version=PLUGIN_VERSION,
        )
        dimensions = scenic_result.provenance["input_dimensions"]
        result = make_table_result(
            title=f"SCENIC final regulons: {transcription_factor or 'all TFs'}",
            operation="scenic_final_regulon_membership",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(dimensions["cells"]),
            input_genes=int(dimensions["genes"]),
            started_at=started_at,
            table=table,
        )
        report = make_summary_result(
            summary=summary,
            title=f"SCENIC final regulons: {transcription_factor or 'all TFs'}",
            operation="scenic_final_regulon_membership",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(dimensions["cells"]),
            input_genes=int(dimensions["genes"]),
            started_at=started_at,
        )
        code = scenic_membership_code(
            parameters={
                "transcription_factor": transcription_factor,
                "max_output_rows": max_output_rows,
                "openbio_version": PLUGIN_VERSION,
            }
        )
        return io.NodeOutput(result, report, code)


REGULATORY_NODE_CLASSES = [
    OpenBioSingleCellCollecTRIULM,
    OpenBioSingleCellRankTFActivities,
    OpenBioSingleCellImportPySCENICResults,
    OpenBioSingleCellSCENICRegulonSpecificity,
    OpenBioSingleCellSCENICActivityBinarization,
    OpenBioSingleCellSCENICTFModules,
]


__all__ = [
    "REGULATORY_NODE_CLASSES",
    "OpenBioSingleCellCollecTRIULM",
    "OpenBioSingleCellImportPySCENICResults",
    "OpenBioSingleCellRankTFActivities",
    "OpenBioSingleCellSCENICActivityBinarization",
    "OpenBioSingleCellSCENICRegulonSpecificity",
    "OpenBioSingleCellSCENICTFModules",
]
