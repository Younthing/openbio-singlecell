from __future__ import annotations

import time
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io

from . import PLUGIN_VERSION
from .analysis_utils import make_plot_result, make_summary_result
from .expression_source import DynamicExpressionSource, ExpressionSourceSpec
from .files import resolve_input_path
from .liana_communication import (
    liana_communication_code,
    liana_resource_cache_fingerprint,
    run_liana_communication,
)
from .liana_plot import liana_dot_plot_code, render_liana_dot_plot
from .liana_result import LianaResult, validate_liana_result
from .node_types import (
    AnnDataType,
    LianaResultType,
    PlotResultType,
    SummaryResultType,
)

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/cell-communication"
METHODS = ["cellphonedb", "rank_aggregate"]
def _label_list(value: str, *, label: str) -> list[str]:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a comma-delimited string.")
    if not value.strip():
        return []
    parts = value.split(",")
    labels = [part.strip() for part in parts]
    if any(not item for item in labels):
        raise ValueError(f"{label} contains an empty comma-delimited item.")
    if len(labels) != len(set(labels)):
        raise ValueError(f"{label} contains duplicate labels.")
    return labels


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
            expected = {"resource", "resource_name", "resource_metadata_json"}
            if set(resource) != expected:
                raise ValueError("LIANA bundled_human resource payload has inactive or missing fields.")
            return {
                "resource_mode": "bundled_human",
                "resource_name": str(resource["resource_name"]),
                "resource_csv": None,
                "resource_metadata_json": str(resource["resource_metadata_json"]),
            }
        if mode == "local_resource":
            expected = {"resource", "resource_name", "resource_csv", "resource_metadata_json"}
            if set(resource) != expected:
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
        return liana_resource_cache_fingerprint(
            path,
            str(resolved["resource_metadata_json"]),
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        sample_key: str = "sample",
        condition_key: str = "condition",
        identity_key: str = "cell_type",
        annotation_status: str = "unknown",
        organism: str = "Homo sapiens",
        method: str = "rank_aggregate",
        resource: Mapping[str, object] | None = None,
        source: DynamicExpressionSource | None = None,
        expression_proportion: float = 0.1,
        min_cells_per_identity_sample: int = 5,
        permutations: int = 1000,
        random_seed: int = 1337,
        jobs: int = 1,
        max_output_rows: int = 2_000_000,
        max_working_memory_gib: float = 4.0,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        resolved_resource = cls._resolve_resource(resource)
        expression = cls.EXPRESSION_SOURCE.resolve(adata, source)
        resource_path = (
            resolve_input_path(str(resolved_resource["resource_csv"]), extensions=(".csv",))
            if resolved_resource["resource_mode"] == "local_resource"
            else None
        )
        result, summary = run_liana_communication(
            adata,
            sample_key=sample_key,
            condition_key=condition_key,
            identity_key=identity_key,
            annotation_status=annotation_status,
            organism=organism,
            method=method,
            resource_mode=str(resolved_resource["resource_mode"]),
            resource_name=str(resolved_resource["resource_name"]),
            resource_path=resource_path,
            resource_metadata_json=str(resolved_resource["resource_metadata_json"]),
            source_kind=expression.kind,
            layer_name=expression.layer_name,
            expression_proportion=expression_proportion,
            min_cells_per_identity_sample=min_cells_per_identity_sample,
            permutations=permutations,
            random_seed=random_seed,
            jobs=jobs,
            max_output_rows=max_output_rows,
            max_working_memory_gib=max_working_memory_gib,
            openbio_version=PLUGIN_VERSION,
        )
        report = make_summary_result(
            summary=summary,
            title="Sample-resolved LIANA communication summary",
            operation="liana_communication",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(adata.n_obs),
            input_genes=int(adata.n_vars),
            started_at=started_at,
            random_seed=random_seed,
        )
        accounting = summary["key_results"]["resource"]["accounting"]
        code = liana_communication_code(
            parameters={
                "sample_key": sample_key,
                "condition_key": condition_key,
                "identity_key": identity_key,
                "annotation_status": annotation_status,
                "organism": organism,
                "method": method,
                "resource_mode": str(resolved_resource["resource_mode"]),
                "resource_name": str(resolved_resource["resource_name"]),
                "resource_path": resource_path,
                "resource_metadata_json": str(resolved_resource["resource_metadata_json"]),
                "source_kind": expression.kind,
                "layer_name": expression.layer_name,
                "expression_proportion": expression_proportion,
                "min_cells_per_identity_sample": min_cells_per_identity_sample,
                "permutations": permutations,
                "random_seed": random_seed,
                "jobs": jobs,
                "max_output_rows": max_output_rows,
                "max_working_memory_gib": max_working_memory_gib,
                "openbio_version": PLUGIN_VERSION,
                "expected_file_sha256": accounting["raw_file_sha256"],
                "expected_canonical_resource_sha256": accounting["canonical_resource_sha256"],
                "expected_effective_resource_sha256": accounting["effective_resource_sha256"],
            }
        )
        return io.NodeOutput(result, report, code)


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
    def execute(
        cls,
        result: LianaResult,
        source_labels: str = "",
        target_labels: str = "",
        selection: Mapping[str, object] | None = None,
        top_n: int = 20,
        figure_width: float = 12.0,
        figure_height: float = 8.0,
        max_plot_rows: int = 100_000,
        max_image_pixels: int = 40_000_000,
    ) -> io.NodeOutput:
        started_at = time.perf_counter()
        _table, provenance, _metadata = validate_liana_result(result)
        selection_method, selection_threshold = cls._resolve_selection(selection)
        if selection_method != provenance["method"]:
            raise ValueError(
                f"LIANA plot selection branch {selection_method!r} does not match the typed result method "
                f"{provenance['method']!r}."
            )
        source_label_list = _label_list(source_labels, label="LIANA source_labels")
        target_label_list = _label_list(target_labels, label="LIANA target_labels")
        png, summary = render_liana_dot_plot(
            result,
            source_labels=source_label_list,
            target_labels=target_label_list,
            selection_threshold=selection_threshold,
            top_n=top_n,
            figure_width=figure_width,
            figure_height=figure_height,
            max_plot_rows=max_plot_rows,
            max_image_pixels=max_image_pixels,
            openbio_version=PLUGIN_VERSION,
        )
        plotted = make_plot_result(
            title=f"LIANA {provenance['method'].replace('_', ' ')} by-Sample dot plot",
            operation="liana_dot_plot",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(provenance["expression"]["cells"]),
            input_genes=int(provenance["expression"]["genes"]),
            started_at=started_at,
            png=png,
        )
        report = make_summary_result(
            summary=summary,
            title="LIANA dot-plot selection summary",
            operation="liana_dot_plot",
            parameters=summary["parameters"],
            description=summary["results"],
            warnings=summary["warnings"],
            input_cells=int(provenance["expression"]["cells"]),
            input_genes=int(provenance["expression"]["genes"]),
            started_at=started_at,
        )
        code = liana_dot_plot_code(
            parameters={
                "source_labels": source_label_list,
                "target_labels": target_label_list,
                "selection_threshold": selection_threshold,
                "top_n": top_n,
                "figure_width": figure_width,
                "figure_height": figure_height,
                "max_plot_rows": max_plot_rows,
                "max_image_pixels": max_image_pixels,
                "openbio_version": PLUGIN_VERSION,
            }
        )
        return io.NodeOutput(plotted, report, code)


COMMUNICATION_NODE_CLASSES = [
    OpenBioSingleCellLianaCommunication,
    OpenBioSingleCellLianaDotPlot,
]


__all__ = [
    "COMMUNICATION_NODE_CLASSES",
    "OpenBioSingleCellLianaCommunication",
    "OpenBioSingleCellLianaDotPlot",
]
