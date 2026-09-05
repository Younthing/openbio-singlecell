from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from comfy_api.latest import io

from .collectri_ulm import collectri_resource_cache_fingerprint
from .expression_source import _COLLECTRI_SPEC
from .files import input_file_fingerprint, resolve_input_path
from .node_types import (
    AnnDataType,
    PlotResultType,
    SCENICResultArtifactType,
    TableResultType,
    TFActivityArtifactType,
    analysis_outputs,
)
from .pyscenic_import import pyscenic_bundle_cache_fingerprint

CATEGORY = "openbio/single-cell/regulatory"
TF_RANKING_METHODS = ("t-test_overestim_var", "t-test", "wilcoxon")


class OpenBioSingleCellCollecTRIULM(io.ComfyNode):
    EXPRESSION_SOURCE = _COLLECTRI_SPEC

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
                io.Float.Input("max_working_memory_gib", default=4.0, min=0.001, max=1024.0, advanced=True),
                io.Boolean.Input("overwrite_existing", default=False, advanced=True),
            ],
            outputs=analysis_outputs(
                AnnDataType.Output(display_name="adata"),
                TFActivityArtifactType.Output(display_name="activities"),
            ),
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
    def prepare_worker_arguments(cls, kwargs: dict[str, Any]) -> dict[str, Any]:
        prepared = dict(kwargs)
        resolved = cls._resolve_resource(prepared.pop("resource", None))
        prepared.update({key: value for key, value in resolved.items() if value is not None})
        return prepared


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
                io.Combo.Input("method", options=list(TF_RANKING_METHODS), default="t-test_overestim_var"),
                io.Float.Input("report_p_adjusted", default=0.05, min=0.0, max=1.0, step=0.01, advanced=True),
                io.Int.Input("max_output_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


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
                io.Int.Input("max_file_bytes", default=2_147_483_647, min=1, max=2**63 - 1, advanced=True),
                io.Int.Input("max_adjacency_edges", default=10_000_000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_regulon_edges", default=2_000_000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input("max_dense_bytes", default=1_073_741_824, min=1, max=2**63 - 1, advanced=True),
            ],
            outputs=analysis_outputs(
                AnnDataType.Output(display_name="adata"),
                SCENICResultArtifactType.Output(display_name="scenic_result"),
            ),
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
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


class OpenBioSingleCellSCENICActivityBinarization(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSCENICActivityBinarization",
            display_name="SCENIC Activity Binarization",
            category=CATEGORY,
            description=(
                "Apply the deterministic one-process pySCENIC 0.12.1 HDT-compatible heuristic to one "
                "validated activity artifact and return complete thresholds plus a strict summary."
            ),
            inputs=[
                SCENICResultArtifactType.Input("scenic_result"),
                io.Int.Input("random_seed", default=1, min=0, max=2**31 - 1),
                TableResultType.Input("threshold_overrides", optional=True),
                io.Int.Input("max_dense_bytes", default=1_073_741_824, min=1, max=2**63 - 1, advanced=True),
            ],
            outputs=analysis_outputs(TableResultType.Output(display_name="thresholds")),
        )

    @classmethod
    def prepare_worker_arguments(cls, kwargs: dict[str, Any]) -> dict[str, Any]:
        prepared = dict(kwargs)
        if prepared.get("threshold_overrides") is None:
            prepared.pop("threshold_overrides", None)
        return prepared


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
            outputs=analysis_outputs(TableResultType.Output(display_name="table")),
        )


def _activity_plot_schema(*, node_id: str, display_name: str, artifact_input: Any) -> io.Schema:
    return io.Schema(
        node_id=node_id,
        display_name=display_name,
        category=CATEGORY,
        inputs=[
            AnnDataType.Input("adata"),
            artifact_input,
            io.DynamicCombo.Input(
                "view",
                options=[
                    io.DynamicCombo.Option(
                        "distributions",
                        [io.Int.Input("max_activities", default=20, min=1, max=50)],
                    ),
                    io.DynamicCombo.Option(
                        "grouped_heatmap",
                        [
                            io.String.Input("groupby", default="cell_type"),
                            io.Int.Input("max_activities", default=30, min=1, max=50),
                        ],
                    ),
                    io.DynamicCombo.Option(
                        "embedding",
                        [
                            io.String.Input("embedding_key", default="X_umap"),
                            io.String.Input("activity_name", default=""),
                        ],
                    ),
                ],
            ),
        ],
        outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
    )


class OpenBioSingleCellTFActivityPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return _activity_plot_schema(
            node_id="OpenBioSingleCellTFActivityPlot",
            display_name="TF Activity Plot",
            artifact_input=TFActivityArtifactType.Input("activities"),
        )


class OpenBioSingleCellSCENICActivityPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return _activity_plot_schema(
            node_id="OpenBioSingleCellSCENICActivityPlot",
            display_name="SCENIC Activity Plot",
            artifact_input=SCENICResultArtifactType.Input("scenic_result"),
        )


def _table_plot_schema(
    *,
    node_id: str,
    display_name: str,
    options: list[io.DynamicCombo.Option],
) -> io.Schema:
    return io.Schema(
        node_id=node_id,
        display_name=display_name,
        category=CATEGORY,
        inputs=[
            TableResultType.Input("table"),
            io.DynamicCombo.Input("view", options=options),
        ],
        outputs=analysis_outputs(PlotResultType.Output(display_name="plot")),
    )


class OpenBioSingleCellTFActivityRankingPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return _table_plot_schema(
            node_id="OpenBioSingleCellTFActivityRankingPlot",
            display_name="TF Activity Ranking Plot",
            options=[
                io.DynamicCombo.Option(name, [io.Int.Input("max_per_group", default=10, min=1, max=100)])
                for name in ("ranked_dot", "effect_significance")
            ],
        )


class OpenBioSingleCellSCENICRegulonSpecificityPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return _table_plot_schema(
            node_id="OpenBioSingleCellSCENICRegulonSpecificityPlot",
            display_name="SCENIC Regulon Specificity Plot",
            options=[
                io.DynamicCombo.Option(name, [io.Int.Input("max_regulons", default=20, min=1, max=100)])
                for name in ("rss_heatmap", "top_regulons")
            ],
        )


class OpenBioSingleCellSCENICBinarizationPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return _table_plot_schema(
            node_id="OpenBioSingleCellSCENICBinarizationPlot",
            display_name="SCENIC Binarization Plot",
            options=[
                io.DynamicCombo.Option(name, [io.Int.Input("max_regulons", default=30, min=1, max=100)])
                for name in ("thresholds", "active_fraction")
            ],
        )


class OpenBioSingleCellSCENICRegulonMembershipPlot(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return _table_plot_schema(
            node_id="OpenBioSingleCellSCENICRegulonMembershipPlot",
            display_name="SCENIC Regulon Membership Plot",
            options=[
                io.DynamicCombo.Option(
                    "target_count",
                    [io.Int.Input("max_regulons", default=30, min=1, max=100)],
                ),
                io.DynamicCombo.Option(
                    "top_targets",
                    [
                        io.String.Input("regulon", default=""),
                        io.Int.Input("max_targets", default=30, min=1, max=100),
                    ],
                ),
            ],
        )


REGULATORY_NODE_CLASSES = [
    OpenBioSingleCellCollecTRIULM,
    OpenBioSingleCellRankTFActivities,
    OpenBioSingleCellImportPySCENICResults,
    OpenBioSingleCellSCENICRegulonSpecificity,
    OpenBioSingleCellSCENICActivityBinarization,
    OpenBioSingleCellSCENICTFModules,
    OpenBioSingleCellTFActivityPlot,
    OpenBioSingleCellTFActivityRankingPlot,
    OpenBioSingleCellSCENICActivityPlot,
    OpenBioSingleCellSCENICRegulonSpecificityPlot,
    OpenBioSingleCellSCENICBinarizationPlot,
    OpenBioSingleCellSCENICRegulonMembershipPlot,
]


__all__ = [
    "REGULATORY_NODE_CLASSES",
    "OpenBioSingleCellCollecTRIULM",
    "OpenBioSingleCellImportPySCENICResults",
    "OpenBioSingleCellRankTFActivities",
    "OpenBioSingleCellTFActivityPlot",
    "OpenBioSingleCellTFActivityRankingPlot",
    "OpenBioSingleCellSCENICActivityPlot",
    "OpenBioSingleCellSCENICActivityBinarization",
    "OpenBioSingleCellSCENICBinarizationPlot",
    "OpenBioSingleCellSCENICRegulonMembershipPlot",
    "OpenBioSingleCellSCENICRegulonSpecificity",
    "OpenBioSingleCellSCENICRegulonSpecificityPlot",
    "OpenBioSingleCellSCENICTFModules",
]
