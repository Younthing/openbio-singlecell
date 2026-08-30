from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from comfy_api.latest import io

from .expression_source import ExpressionSourceSpec
from .files import resolve_input_path
from .liana_communication import liana_resource_cache_fingerprint
from .node_types import AnnDataType, LianaResultType, PlotResultType, SummaryResultType

CATEGORY = "openbio/single-cell/cell-communication"
METHODS = ["cellphonedb", "rank_aggregate"]


class OpenBioSingleCellLianaCommunication(io.ComfyNode):
    EXPRESSION_SOURCE = ExpressionSourceSpec(
        description="LIANA expression source",
        include_raw=False,
        layer_default="log1p_norm",
    )

    @staticmethod
    def _resolve_resource(resource: Mapping[str, object] | None) -> dict[str, str | None]:
        if resource is None:
            resource = {
                "resource": "bundled_human",
                "resource_name": "consensus",
                "resource_metadata_json": "{}",
            }
        if not isinstance(resource, Mapping):
            raise TypeError("LIANA resource must be a DynamicCombo value.")
        mode = resource.get("resource")
        if mode == "bundled_human":
            if set(resource) != {"resource", "resource_name", "resource_metadata_json"}:
                raise ValueError("LIANA bundled_human resource payload has inactive or missing fields.")
            return {
                "resource_mode": "bundled_human",
                "resource_name": str(resource["resource_name"]),
                "resource_csv": None,
                "resource_metadata_json": str(resource["resource_metadata_json"]),
            }
        if mode == "local_resource":
            if set(resource) != {"resource", "resource_name", "resource_csv", "resource_metadata_json"}:
                raise ValueError("LIANA local_resource payload has inactive or missing fields.")
            return {
                "resource_mode": "local_resource",
                "resource_name": str(resource["resource_name"]),
                "resource_csv": str(resource["resource_csv"]),
                "resource_metadata_json": str(resource["resource_metadata_json"]),
            }
        raise ValueError(f"Unsupported LIANA resource mode: {mode!r}.")

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLianaCommunication",
            display_name="LIANA Communication",
            category=CATEGORY,
            description=(
                "Run exact LIANA 1.9.0 candidate communication separately within each biological Sample using "
                "one explicit, license-reviewed and SHA-256-bound ligand-receptor resource."
            ),
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("sample_key", default="sample"),
                io.String.Input("condition_key", default="condition"),
                io.String.Input("identity_key", default="cell_type"),
                io.Combo.Input(
                    "annotation_status",
                    options=["unknown", "provisional", "curated"],
                    default="unknown",
                ),
                io.String.Input("organism", default="Homo sapiens"),
                io.Combo.Input("method", options=METHODS, default="rank_aggregate"),
                io.DynamicCombo.Input(
                    "resource",
                    options=[
                        io.DynamicCombo.Option(
                            "bundled_human",
                            [
                                io.String.Input("resource_name", default="consensus"),
                                io.String.Input("resource_metadata_json", default="{}", multiline=True),
                            ],
                        ),
                        io.DynamicCombo.Option(
                            "local_resource",
                            [
                                io.String.Input("resource_name", default="local_ligand_receptor_snapshot"),
                                io.String.Input("resource_csv", default=""),
                                io.String.Input("resource_metadata_json", default="{}", multiline=True),
                            ],
                        ),
                    ],
                ),
                cls.EXPRESSION_SOURCE.input(),
                io.Float.Input("expression_proportion", default=0.1, min=0.0, max=1.0, step=0.05),
                io.Int.Input("min_cells_per_identity_sample", default=5, min=2, max=2**31 - 1),
                io.Int.Input("permutations", default=1000, min=1000, max=10_000_000),
                io.Int.Input("random_seed", default=1337, min=1, max=2**32 - 1),
                io.Int.Input("jobs", default=1, min=1, max=1, advanced=True),
                io.Int.Input("max_output_rows", default=2_000_000, min=1, max=2**31 - 1, advanced=True),
                io.Float.Input(
                    "max_working_memory_gib",
                    default=4.0,
                    min=0.001,
                    max=1024.0,
                    advanced=True,
                ),
            ],
            outputs=[
                LianaResultType.Output(display_name="result"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def validate_inputs(cls, resource: Mapping[str, object] | None = None, **kwargs: Any) -> bool | str:
        try:
            resolved = cls._resolve_resource(resource)
            if resolved["resource_mode"] == "local_resource":
                resolve_input_path(str(resolved["resource_csv"]), extensions=(".csv",))
        except (TypeError, ValueError, FileNotFoundError, OSError) as exc:
            return str(exc)
        return True

    @classmethod
    def fingerprint_inputs(cls, resource: Mapping[str, object] | None = None, **kwargs: Any) -> Any:
        resolved = cls._resolve_resource(resource)
        if resolved["resource_mode"] == "bundled_human":
            return (
                "openbio-liana-bundled-v1",
                resolved["resource_name"],
                resolved["resource_metadata_json"],
            )
        path = resolve_input_path(str(resolved["resource_csv"]), extensions=(".csv",))
        return liana_resource_cache_fingerprint(path, str(resolved["resource_metadata_json"]))

    @classmethod
    def prepare_worker_arguments(cls, kwargs: dict[str, Any]) -> dict[str, Any]:
        prepared = dict(kwargs)
        resolved = cls._resolve_resource(prepared.pop("resource", None))
        prepared.update({key: value for key, value in resolved.items() if value is not None})
        return prepared


class OpenBioSingleCellLianaDotPlot(io.ComfyNode):
    @staticmethod
    def _resolve_selection(selection: Mapping[str, object] | None) -> tuple[str, float]:
        if selection is None:
            selection = {"selection": "rank_aggregate", "max_specificity_rank": 0.05}
        if not isinstance(selection, Mapping):
            raise TypeError("LIANA plot selection must be a DynamicCombo value.")
        mode = selection.get("selection")
        if mode == "rank_aggregate":
            if set(selection) != {"selection", "max_specificity_rank"}:
                raise ValueError("LIANA rank_aggregate selection payload has inactive or missing fields.")
            value = selection["max_specificity_rank"]
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError("LIANA max_specificity_rank must be a real number.")
            return "rank_aggregate", float(value)
        if mode == "cellphonedb":
            if set(selection) != {"selection", "max_cellphone_pvalue"}:
                raise ValueError("LIANA cellphonedb selection payload has inactive or missing fields.")
            value = selection["max_cellphone_pvalue"]
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError("LIANA max_cellphone_pvalue must be a real number.")
            return "cellphonedb", float(value)
        raise ValueError(f"Unsupported LIANA plot selection branch: {mode!r}.")

    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellLianaDotPlot",
            display_name="LIANA Dot Plot",
            category=CATEGORY,
            description=(
                "Render a deterministic Sample-faceted view of a validated typed LIANA result without rerunning "
                "communication inference."
            ),
            inputs=[
                LianaResultType.Input("result"),
                io.String.Input("source_labels", default=""),
                io.String.Input("target_labels", default=""),
                io.DynamicCombo.Input(
                    "selection",
                    options=[
                        io.DynamicCombo.Option(
                            "rank_aggregate",
                            [
                                io.Float.Input(
                                    "max_specificity_rank",
                                    default=0.05,
                                    min=0.0,
                                    max=1.0,
                                    step=0.01,
                                )
                            ],
                        ),
                        io.DynamicCombo.Option(
                            "cellphonedb",
                            [
                                io.Float.Input(
                                    "max_cellphone_pvalue",
                                    default=0.05,
                                    min=0.0,
                                    max=1.0,
                                    step=0.01,
                                )
                            ],
                        ),
                    ],
                ),
                io.Int.Input("top_n", default=20, min=1, max=2**31 - 1),
                io.Float.Input("figure_width", default=12.0, min=1.0, max=30.0, step=1.0, advanced=True),
                io.Float.Input("figure_height", default=8.0, min=1.0, max=30.0, step=1.0, advanced=True),
                io.Int.Input("max_plot_rows", default=100_000, min=1, max=2**31 - 1, advanced=True),
                io.Int.Input(
                    "max_image_pixels",
                    default=40_000_000,
                    min=1,
                    max=2**31 - 1,
                    advanced=True,
                ),
            ],
            outputs=[
                PlotResultType.Output(display_name="plot"),
                SummaryResultType.Output(display_name="summary"),
                io.String.Output("code"),
            ],
        )

    @classmethod
    def prepare_worker_arguments(cls, kwargs: dict[str, Any]) -> dict[str, Any]:
        prepared = dict(kwargs)
        method, threshold = cls._resolve_selection(prepared.pop("selection", None))
        prepared.update(selection_method=method, selection_threshold=threshold)
        return prepared


COMMUNICATION_NODE_CLASSES = [
    OpenBioSingleCellLianaCommunication,
    OpenBioSingleCellLianaDotPlot,
]

__all__ = [
    "COMMUNICATION_NODE_CLASSES",
    "OpenBioSingleCellLianaCommunication",
    "OpenBioSingleCellLianaDotPlot",
]
