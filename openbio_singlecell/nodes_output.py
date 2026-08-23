from __future__ import annotations

from typing import TYPE_CHECKING, Any

from comfy_api.latest import io, ui

from .contracts import PlotResult, SingleCellResult, SummaryResult, TableResult
from .files import OutputTarget, prepare_output_target, preview_output_target
from .node_types import AnnDataType, PlotResultType, SummaryResultType, TableResultType
from .payload import result_to_payload

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/output"


def _saved_result(target: OutputTarget) -> ui.SavedResult:
    folder_type = io.FolderType.output if target.folder_type == "output" else io.FolderType.temp
    return ui.SavedResult(target.filename, target.subfolder, folder_type)


def _preview_identity(unique_id: str | int | None, extra_pnginfo: Any) -> str:
    workflow_id = None
    if isinstance(extra_pnginfo, dict):
        workflow = extra_pnginfo.get("workflow")
        if isinstance(workflow, dict):
            candidate = workflow.get("id")
            if isinstance(candidate, (str, int)) and not isinstance(candidate, bool):
                candidate = str(candidate).strip()
                if candidate:
                    workflow_id = candidate

    workflow_key = workflow_id or "workflow-unknown"
    node_key = str(unique_id).strip() if unique_id is not None else ""
    node_key = node_key or "preview"
    return f"{len(workflow_key)}:{workflow_key}{len(node_key)}:{node_key}"


def preview_result(
    result: SingleCellResult,
    unique_id: str | int | None = None,
    extra_pnginfo: Any = None,
) -> dict[str, Any]:
    if not isinstance(result, (SummaryResult, TableResult, PlotResult)):
        raise TypeError("Expected an OpenBio summary, table, or plot result value.")
    output: dict[str, Any] = {"openbio_singlecell": [result_to_payload(result)]}
    if isinstance(result, PlotResult):
        target = preview_output_target(_preview_identity(unique_id, extra_pnginfo))
        with open(target.path, "wb") as handle:
            handle.write(result.png)
        output["images"] = [_saved_result(target)]
    return output


def save_h5ad(adata: AnnData, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    target = prepare_output_target(filename_prefix, "h5ad", overwrite)
    adata.write_h5ad(target.path)
    return target


def export_csv(table: TableResult, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    if not isinstance(table, TableResult):
        raise TypeError("Expected an OPENBIO_SINGLE_CELL_TABLE value.")
    target = prepare_output_target(filename_prefix, "csv", overwrite)
    table.table.to_csv(target.path, index=False)
    return target


def save_png(plot: PlotResult, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    if not isinstance(plot, PlotResult):
        raise TypeError("Expected an OPENBIO_SINGLE_CELL_PLOT value.")
    target = prepare_output_target(filename_prefix, "png", overwrite)
    with open(target.path, "wb") as handle:
        handle.write(plot.png)
    return target


class OpenBioSingleCellPreviewResult(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPreviewResult",
            display_name="Preview Result",
            category=CATEGORY,
            inputs=[io.MultiType.Input("result", [SummaryResultType, TableResultType, PlotResultType])],
            outputs=[],
            hidden=[io.Hidden.unique_id, io.Hidden.extra_pnginfo],
            is_output_node=True,
        )

    @classmethod
    def execute(cls, result: SingleCellResult) -> io.NodeOutput:
        return io.NodeOutput(ui=preview_result(result, cls.hidden.unique_id, cls.hidden.extra_pnginfo))


class OpenBioSingleCellSaveH5AD(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSaveH5AD",
            display_name="Save H5AD",
            category=CATEGORY,
            inputs=[
                AnnDataType.Input("adata"),
                io.String.Input("filename_prefix", default="adata"),
                io.Boolean.Input("overwrite", default=False, advanced=True),
            ],
            outputs=[],
            is_output_node=True,
            not_idempotent=True,
        )

    @classmethod
    def execute(cls, adata: AnnData, filename_prefix: str = "adata", overwrite: bool = False) -> io.NodeOutput:
        target = save_h5ad(adata, filename_prefix, overwrite)
        return io.NodeOutput(ui={"files": [_saved_result(target)]})


class OpenBioSingleCellExportCSV(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellExportCSV",
            display_name="Export CSV",
            category=CATEGORY,
            inputs=[
                TableResultType.Input("table"),
                io.String.Input("filename_prefix", default="table"),
                io.Boolean.Input("overwrite", default=False, advanced=True),
            ],
            outputs=[],
            is_output_node=True,
            not_idempotent=True,
        )

    @classmethod
    def execute(
        cls,
        table: TableResult,
        filename_prefix: str = "table",
        overwrite: bool = False,
    ) -> io.NodeOutput:
        target = export_csv(table, filename_prefix, overwrite)
        return io.NodeOutput(ui={"files": [_saved_result(target)]})


class OpenBioSingleCellSavePNG(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSavePNG",
            display_name="Save PNG",
            category=CATEGORY,
            inputs=[
                PlotResultType.Input("plot"),
                io.String.Input("filename_prefix", default="plot"),
                io.Boolean.Input("overwrite", default=False, advanced=True),
            ],
            outputs=[],
            is_output_node=True,
            not_idempotent=True,
        )

    @classmethod
    def execute(
        cls,
        plot: PlotResult,
        filename_prefix: str = "plot",
        overwrite: bool = False,
    ) -> io.NodeOutput:
        target = save_png(plot, filename_prefix, overwrite)
        return io.NodeOutput(ui={"images": [_saved_result(target)]})


OUTPUT_NODE_CLASSES = [
    OpenBioSingleCellPreviewResult,
    OpenBioSingleCellSaveH5AD,
    OpenBioSingleCellExportCSV,
    OpenBioSingleCellSavePNG,
]
