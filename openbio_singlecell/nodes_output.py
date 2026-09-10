from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import Any

import folder_paths
from comfy_api.latest import io, ui

from .artifact_codecs import (
    ANNDATA_CODEC,
    ANNDATA_PAYLOAD,
    PLOT_CODEC,
    PLOT_PAYLOAD,
    RESULT_METADATA,
    TABLE_CODEC,
    TABLE_PAYLOAD,
    TABLE_SCHEMA,
)
from .artifact_persist import persist_artifact
from .artifact_runtime import ArtifactTicket
from .artifact_service import current_artifact_runtime
from .contracts import SummaryResult
from .files import OutputTarget, atomic_write_output, prepare_output_target, preview_output_target
from .node_types import (
    AnnDataType,
    AugurResultType,
    CassiopeiaCharactersType,
    CassiopeiaTreeType,
    CNVStateType,
    CompositionModelResultType,
    DGIdbResourceType,
    LianaResultType,
    MiloResultType,
    Monocle2CDSType,
    PlotResultType,
    PseudobulkType,
    SCENICResultArtifactType,
    SummaryResultType,
    TableResultType,
    TFActivityArtifactType,
    VelocityStateType,
)
from .payload import MAX_COLUMNS, MAX_ROWS, artifact_metadata_to_payload, result_to_payload

CATEGORY = "openbio/single-cell/output"
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


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"Artifact JSON contains a non-finite number: {value}")


def _strict_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_json_constant)
    if not isinstance(value, dict):
        raise ValueError(f"Artifact JSON must be an object: {path.name}")
    return value


def _ticket_root(ticket: ArtifactTicket, *, kind: str, codec: str) -> Path:
    if not isinstance(ticket, ArtifactTicket):
        raise TypeError(f"Expected a {kind} ArtifactTicket.")
    if ticket.kind != kind or ticket.codec != codec:
        raise ValueError(f"Expected {kind} with codec {codec}; received {ticket.kind} with {ticket.codec}.")
    return current_artifact_runtime().resolve(ticket)


def _copy_artifact_file(source: Path, target: OutputTarget, *, overwrite: bool) -> None:
    source_size = source.stat().st_size

    def writer(staged_path: str) -> None:
        with source.open("rb") as source_handle, open(staged_path, "wb") as staged_handle:
            shutil.copyfileobj(source_handle, staged_handle, length=1024 * 1024)

    def validator(staged_path: str) -> None:
        if source.stat().st_size != source_size or Path(staged_path).stat().st_size != source_size:
            raise RuntimeError("Artifact source changed while it was being copied.")

    atomic_write_output(target, writer, overwrite=overwrite, validator=validator)


def _table_rows(root: Path):
    schema = _strict_json(root / TABLE_SCHEMA)
    columns = schema.get("columns")
    index = schema.get("index")
    if not isinstance(columns, list) or not isinstance(index, dict):
        raise ValueError("Table artifact schema is invalid.")
    column_names = [item.get("name") if isinstance(item, dict) else None for item in columns]
    if not all(isinstance(name, str) for name in column_names):
        raise ValueError("Table artifact columns are invalid.")
    index_name = index.get("name")
    index_label = index_name if isinstance(index_name, str) and index_name.strip() else "__openbio_row_id__"
    if index_label in column_names:
        raise ValueError(
            f"CSV index label collides with a data column; rename it before export: {index_label}"
        )

    def rows():
        with (root / TABLE_PAYLOAD).open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                row = json.loads(line, parse_constant=_reject_json_constant)
                if not isinstance(row, list) or len(row) != len(column_names) + 1:
                    raise ValueError(f"Table artifact row {line_number} is invalid.")
                yield row

    return index_label, column_names, rows()


def _csv_value(value: Any) -> Any:
    if value is None:
        return CSV_MISSING_VALUE
    if isinstance(value, (list, dict)):
        return json.dumps(value, allow_nan=False, ensure_ascii=False, separators=(",", ":"))
    return value


def preview_result(
    result: SummaryResult | ArtifactTicket,
    unique_id: str | int | None = None,
    extra_pnginfo: Any = None,
) -> dict[str, Any]:
    if isinstance(result, SummaryResult):
        return {"openbio_singlecell": [result_to_payload(result)]}
    if not isinstance(result, ArtifactTicket):
        raise TypeError("Expected an OpenBio summary or table/plot ArtifactTicket.")
    if result.kind == "OPENBIO_SINGLE_CELL_TABLE":
        root = _ticket_root(result, kind=result.kind, codec=TABLE_CODEC)
        index_name, columns, rows = _table_rows(root)
        preview_rows = []
        index_values = []
        total_rows = 0
        for row in rows:
            if total_rows < MAX_ROWS:
                index_values.append(row[0])
                preview_rows.append(row[1 : MAX_COLUMNS + 1])
            total_rows += 1
        payload = artifact_metadata_to_payload(
            _strict_json(root / RESULT_METADATA),
            columns=columns,
            rows=preview_rows,
            total_rows=total_rows,
            index_name=index_name,
            index_values=index_values,
        )
        return {"openbio_singlecell": [payload]}
    if result.kind == "OPENBIO_SINGLE_CELL_PLOT":
        root = _ticket_root(result, kind=result.kind, codec=PLOT_CODEC)
        target = preview_output_target(_preview_identity(unique_id, extra_pnginfo))
        _copy_artifact_file(root / PLOT_PAYLOAD, target, overwrite=True)
        payload = artifact_metadata_to_payload(_strict_json(root / RESULT_METADATA))
        return {"openbio_singlecell": [payload], "images": [_saved_result(target)]}
    raise ValueError("Preview accepts only OPENBIO_SINGLE_CELL_TABLE or OPENBIO_SINGLE_CELL_PLOT tickets.")


def save_h5ad(adata: ArtifactTicket, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    root = _ticket_root(adata, kind="OPENBIO_ANNDATA", codec=ANNDATA_CODEC)
    target = prepare_output_target(filename_prefix, "h5ad", overwrite)
    _copy_artifact_file(root / ANNDATA_PAYLOAD, target, overwrite=overwrite)
    return target


def export_csv(table: ArtifactTicket, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    root = _ticket_root(table, kind="OPENBIO_SINGLE_CELL_TABLE", codec=TABLE_CODEC)
    target = prepare_output_target(filename_prefix, "csv", overwrite)
    index_label, columns, rows = _table_rows(root)

    def writer(staged_path: str) -> None:
        with open(staged_path, "w", encoding="utf-8", newline="") as handle:
            output = csv.writer(handle, lineterminator="\n")
            output.writerow([index_label, *columns])
            for row in rows:
                output.writerow([_csv_value(value) for value in row])

    atomic_write_output(target, writer, overwrite=overwrite)
    return target


def save_png(plot: ArtifactTicket, filename_prefix: str, overwrite: bool = False) -> OutputTarget:
    root = _ticket_root(plot, kind="OPENBIO_SINGLE_CELL_PLOT", codec=PLOT_CODEC)
    target = prepare_output_target(filename_prefix, "png", overwrite)
    _copy_artifact_file(root / PLOT_PAYLOAD, target, overwrite=overwrite)
    return target


class OpenBioSingleCellPreviewResult(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPreviewResult",
            display_name="Preview Result",
            category=CATEGORY,
            description="Preview a summary, table, or plot result without creating a durable output.",
            inputs=[io.MultiType.Input("result", [SummaryResultType, TableResultType, PlotResultType])],
            outputs=[],
            hidden=[io.Hidden.unique_id, io.Hidden.extra_pnginfo],
            is_output_node=True,
        )

    @classmethod
    def execute(cls, result: SummaryResult | ArtifactTicket) -> io.NodeOutput:
        return io.NodeOutput(ui=preview_result(result, cls.hidden.unique_id, cls.hidden.extra_pnginfo))


class OpenBioSingleCellSaveH5AD(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellSaveH5AD",
            display_name="Save H5AD",
            category=CATEGORY,
            description="Save an AnnData artifact as an H5AD file in ComfyUI output storage.",
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
    def execute(
        cls,
        adata: ArtifactTicket,
        filename_prefix: str = "adata",
        overwrite: bool = False,
    ) -> io.NodeOutput:
        target = save_h5ad(adata, filename_prefix, overwrite)
        return io.NodeOutput(ui={"files": [_saved_result(target)]})


class OpenBioSingleCellExportCSV(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellExportCSV",
            display_name="Export CSV",
            category=CATEGORY,
            description="Export a complete table artifact, including row identity, as CSV in ComfyUI output storage.",
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
        table: ArtifactTicket,
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
            description="Save a plot artifact as a PNG file in ComfyUI output storage.",
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
        plot: ArtifactTicket,
        filename_prefix: str = "plot",
        overwrite: bool = False,
    ) -> io.NodeOutput:
        target = save_png(plot, filename_prefix, overwrite)
        return io.NodeOutput(ui={"images": [_saved_result(target)]})


class OpenBioSingleCellPersistArtifact(io.ComfyNode):
    @classmethod
    def define_schema(cls) -> io.Schema:
        return io.Schema(
            node_id="OpenBioSingleCellPersistArtifact",
            display_name="Persist Artifact",
            category=CATEGORY,
            description="Copy a temporary portable file artifact into durable ComfyUI output storage.",
            inputs=[
                io.MultiType.Input(
                    "artifact",
                    [
                        AnnDataType,
                        TableResultType,
                        PlotResultType,
                        PseudobulkType,
                        CNVStateType,
                        Monocle2CDSType,
                        VelocityStateType,
                        TFActivityArtifactType,
                        SCENICResultArtifactType,
                        AugurResultType,
                        DGIdbResourceType,
                        LianaResultType,
                        CassiopeiaCharactersType,
                        CassiopeiaTreeType,
                        MiloResultType,
                        CompositionModelResultType,
                    ],
                ),
                io.String.Input("name", default="artifact"),
            ],
            outputs=[],
            is_output_node=True,
            not_idempotent=True,
        )

    @classmethod
    def execute(cls, artifact: ArtifactTicket, name: str = "artifact") -> io.NodeOutput:
        persisted = persist_artifact(
            current_artifact_runtime(),
            artifact,
            folder_paths.get_output_directory(),
            name,
        )
        return io.NodeOutput(ui={"text": [str(persisted.path)]})


OUTPUT_NODE_CLASSES = [
    OpenBioSingleCellPreviewResult,
    OpenBioSingleCellSaveH5AD,
    OpenBioSingleCellExportCSV,
    OpenBioSingleCellSavePNG,
    OpenBioSingleCellPersistArtifact,
]
