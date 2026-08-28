from __future__ import annotations

import csv
import warnings
from io import BytesIO
from typing import TYPE_CHECKING, Any

from comfy_api.latest import io, ui

from .contracts import PlotResult, SingleCellResult, SummaryResult, TableResult
from .files import OutputTarget, atomic_write_output, prepare_output_target, preview_output_target
from .node_types import AnnDataType, PlotResultType, SummaryResultType, TableResultType
from .payload import result_to_payload

if TYPE_CHECKING:
    from anndata import AnnData


CATEGORY = "openbio/single-cell/output"
H5AD_COMPRESSIONS = ("gzip", "lzf", "none")
CSV_MISSING_VALUE = "NA"


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


def _validate_png_bytes(png: bytes) -> None:
    if not isinstance(png, bytes):
        raise TypeError("PNG payload must be bytes.")
    try:
        from PIL import Image

        with warnings.catch_warnings():
            warnings.simplefilter("error", Image.DecompressionBombWarning)
            with Image.open(BytesIO(png)) as image:
                if image.format != "PNG":
                    raise ValueError(f"Detected {image.format or 'unknown'} image data instead of PNG.")
                image.verify()
            with Image.open(BytesIO(png)) as image:
                if image.format != "PNG":
                    raise ValueError(f"Detected {image.format or 'unknown'} image data instead of PNG.")
                image.load()
                if image.width < 1 or image.height < 1:
                    raise ValueError("PNG dimensions must both be positive.")
    except Exception as error:
        raise ValueError("PlotResult contains an invalid, truncated, or unsafe PNG payload.") from error


def _write_png(target: OutputTarget, png: bytes, *, overwrite: bool) -> None:
    _validate_png_bytes(png)

    def writer(staged_path: str) -> None:
        with open(staged_path, "wb") as handle:
            handle.write(png)

    def validator(staged_path: str) -> None:
        with open(staged_path, "rb") as handle:
            staged_png = handle.read()
        if staged_png != png:
            raise RuntimeError("Staged PNG bytes differ from the supplied PlotResult.")
        _validate_png_bytes(staged_png)

    atomic_write_output(target, writer, overwrite=overwrite, validator=validator)


def _h5ad_compression(value: str) -> str | None:
    if not isinstance(value, str) or value not in H5AD_COMPRESSIONS:
        allowed = ", ".join(H5AD_COMPRESSIONS)
        raise ValueError(f"H5AD compression must be one of: {allowed}.")
    return None if value == "none" else value


def _csv_label(value: Any, *, role: str) -> str:
    if value is None:
        raise ValueError(f"{role} cannot be None; name it explicitly before CSV export.")
    label = str(value)
    if not label.strip():
        raise ValueError(f"{role} cannot be empty; name it explicitly before CSV export.")
    if "\x00" in label:
        raise ValueError(f"{role} cannot contain a NUL character.")
    return label


def _prepare_csv_frame(frame: Any) -> tuple[Any, list[str]]:
    if int(frame.columns.nlevels) != 1:
        raise ValueError("CSV export requires a single-level column axis; flatten MultiIndex columns first.")
    column_labels = [_csv_label(value, role="CSV column label") for value in frame.columns]
    if len(column_labels) != len(set(column_labels)):
        raise ValueError("CSV column labels must be unique after string conversion.")

    index_labels = []
    for level, name in enumerate(frame.index.names):
        if name is None or (isinstance(name, str) and not name.strip()):
            label = "__openbio_row_id__" if frame.index.nlevels == 1 else f"__openbio_row_id_{level}__"
        else:
            label = _csv_label(name, role=f"CSV index level {level} name")
        index_labels.append(label)
    if len(index_labels) != len(set(index_labels)):
        raise ValueError("CSV index labels must be unique.")
    collisions = sorted(set(index_labels).intersection(column_labels))
    if collisions:
        raise ValueError(
            "CSV index labels collide with data columns; rename them before export: " + ", ".join(collisions)
        )

    export_frame = frame.copy(deep=False)
    export_frame.columns = column_labels
    return export_frame, index_labels


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
        _write_png(target, result.png, overwrite=True)
        output["images"] = [_saved_result(target)]
    return output


def save_h5ad(
    adata: AnnData,
    filename_prefix: str,
    overwrite: bool = False,
    compression: str = "gzip",
) -> OutputTarget:
    from anndata import AnnData as AnnDataClass
    from anndata import read_h5ad

    if not isinstance(adata, AnnDataClass):
        raise TypeError("Expected an OPENBIO_ANNDATA value.")
    selected_compression = _h5ad_compression(compression)
    target = prepare_output_target(filename_prefix, "h5ad", overwrite)
    expected_shape = tuple(int(value) for value in adata.shape)
    expected_obs_names = adata.obs_names.copy()
    expected_var_names = adata.var_names.copy()

    def writer(staged_path: str) -> None:
        adata.write_h5ad(staged_path, compression=selected_compression)

    def validator(staged_path: str) -> None:
        if tuple(int(value) for value in adata.shape) != expected_shape:
            raise RuntimeError("Input AnnData shape changed while it was being written.")
        if not adata.obs_names.equals(expected_obs_names) or not adata.var_names.equals(expected_var_names):
            raise RuntimeError("Input AnnData axis identity changed while it was being written.")
        restored = read_h5ad(staged_path, backed="r")
        try:
            if tuple(int(value) for value in restored.shape) != expected_shape:
                raise RuntimeError("Written H5AD shape does not match the source AnnData.")
            if not restored.obs_names.equals(expected_obs_names):
                raise RuntimeError("Written H5AD observation identifiers do not match the source AnnData.")
            if not restored.var_names.equals(expected_var_names):
                raise RuntimeError("Written H5AD variable identifiers do not match the source AnnData.")
        finally:
            restored.file.close()

    atomic_write_output(target, writer, overwrite=overwrite, validator=validator)
    return target


def export_csv(table: TableResult, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    if not isinstance(table, TableResult):
        raise TypeError("Expected an OPENBIO_SINGLE_CELL_TABLE value.")
    frame, index_labels = _prepare_csv_frame(table.table)
    target = prepare_output_target(filename_prefix, "csv", overwrite)
    expected_header = index_labels + list(frame.columns)
    expected_index = table.table.index.copy()
    expected_columns = table.table.columns.copy()

    def writer(staged_path: str) -> None:
        index_label: str | list[str] = index_labels[0] if len(index_labels) == 1 else index_labels
        frame.to_csv(
            staged_path,
            index=True,
            index_label=index_label,
            encoding="utf-8",
            lineterminator="\n",
            na_rep=CSV_MISSING_VALUE,
        )

    def validator(staged_path: str) -> None:
        if not table.table.index.equals(expected_index) or not table.table.columns.equals(expected_columns):
            raise RuntimeError("Input TableResult axes changed while it was being exported.")
        with open(staged_path, encoding="utf-8", newline="") as handle:
            rows = csv.reader(handle)
            header = next(rows, None)
            if header != expected_header:
                raise RuntimeError("Staged CSV header does not match the explicit row/column identity contract.")
            row_count = sum(1 for _ in rows)
        if row_count != len(frame):
            raise RuntimeError("Staged CSV row count does not match the source table.")

    atomic_write_output(target, writer, overwrite=overwrite, validator=validator)
    return target


def save_png(plot: PlotResult, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    if not isinstance(plot, PlotResult):
        raise TypeError("Expected an OPENBIO_SINGLE_CELL_PLOT value.")
    target = prepare_output_target(filename_prefix, "png", overwrite)
    _write_png(target, plot.png, overwrite=overwrite)
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
                io.Combo.Input("compression", options=list(H5AD_COMPRESSIONS), default="gzip", advanced=True),
            ],
            outputs=[],
            is_output_node=True,
            not_idempotent=True,
        )

    @classmethod
    def execute(
        cls,
        adata: AnnData,
        filename_prefix: str = "adata",
        overwrite: bool = False,
        compression: str = "gzip",
    ) -> io.NodeOutput:
        target = save_h5ad(adata, filename_prefix, overwrite, compression)
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
