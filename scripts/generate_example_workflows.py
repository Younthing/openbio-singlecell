from __future__ import annotations

import argparse
import json
import os
import sys
import uuid
from collections.abc import Sequence
from pathlib import Path
from typing import Any

PLUGIN_ROOT = Path(__file__).resolve().parents[1]


def _resolve_comfy_root() -> Path:
    configured = os.environ.get("OPENBIO_COMFYUI_ROOT")
    if configured:
        candidates = [Path(configured)]
        source = "OPENBIO_COMFYUI_ROOT"
    else:
        candidates = []
        if PLUGIN_ROOT.parent.name.lower() == "custom_nodes":
            candidates.append(PLUGIN_ROOT.parent.parent)
        candidates.append(PLUGIN_ROOT.parent / "ComfyUI")
        source = None

    for candidate in candidates:
        if (candidate / "main.py").is_file():
            return candidate.resolve()

    if source is not None:
        raise RuntimeError(f"{source} does not point to a ComfyUI root containing main.py: {candidates[0]}")
    raise RuntimeError("ComfyUI root was not found. Set OPENBIO_COMFYUI_ROOT before generating workflows.")


sys.path.insert(0, str(_resolve_comfy_root()))
sys.path.insert(0, str(PLUGIN_ROOT))

from comfy_api.latest import io  # noqa: E402

from openbio_singlecell.extension import NODE_CLASSES  # noqa: E402

EXAMPLE_DIRECTORY = PLUGIN_ROOT / "example_workflows"
DEMO_PATH = "openbio-singlecell/openbio_singlecell_demo.h5ad"
BEST_PRACTICE_PATH = "openbio-singlecell/anndata_qc.h5ad"
DEMO_STUDY_PARAMETERS_JSON = json.dumps(
    {
        "schema_version": 1,
        "sample_column": "sample",
        "condition_column": "condition",
        "batch_column": "batch",
        "annotation_column": "cell_type",
        "reference": "control",
        "comparison": "treated",
    },
    separators=(",", ":"),
)
BEST_PRACTICE_STUDY_PARAMETERS_JSON = json.dumps(
    {
        "schema_version": 1,
        "sample_column": "sample",
        "condition_column": "group",
        "batch_column": "batch",
        "annotation_column": "celltypist_cell_type",
        "reference": "Normal",
        "comparison": "nonDM_ED",
    },
    separators=(",", ":"),
)
BEST_PRACTICE_MARKER_PANEL = (
    "PECAM1,VWF,KDR,COL1A1,DCN,RGS5,CSPG4,ACTA2,TAGLN,SOX10,S100B,"
    "CD68,LYZ,CD3D,CD3E,NKG7,MS4A1,CD79A"
)
GROUP_COLORS = ["#3f789e", "#2f806d", "#76558f", "#9a6138", "#526d82", "#775b52"]

SCHEMAS = {node.GET_SCHEMA().node_id: node.GET_SCHEMA() for node in NODE_CLASSES}


def _is_widget_input(item: Any) -> bool:
    return isinstance(item, io.WidgetInput) or item.get_io_type() == "COMFY_DYNAMICCOMBO_V3"


def _is_wire_input(item: Any, wired_inputs: frozenset[str]) -> bool:
    return not _is_widget_input(item) or item.id in wired_inputs


def _widget_entries(item: Any, overrides: dict[str, Any]) -> list[tuple[str, Any, bool]]:
    if item.get_io_type() != "COMFY_DYNAMICCOMBO_V3":
        return [(item.id, overrides.get(item.id, item.default), item.advanced is True)]

    selected = overrides.get(item.id, item.options[0].key)
    option = next((option for option in item.options if option.key == selected), None)
    if option is None:
        raise ValueError(f"{item.id} has no dynamic option named {selected!r}")

    entries = [(item.id, selected, item.advanced is True)]
    for nested in option.inputs:
        name = f"{item.id}.{nested.id}"
        value = overrides.get(name, nested.default)
        entries.append((name, value, nested.advanced is True))
    return entries


def _node(
    alias: str,
    node_type: str,
    pos: tuple[int, int],
    widgets: dict[str, Any] | None = None,
    size: tuple[int, int] | None = None,
    *,
    wired_inputs: Sequence[str] = (),
) -> dict[str, Any]:
    schema = SCHEMAS[node_type]
    wired_input_ids = frozenset(wired_inputs)
    schema_input_ids = {item.id for item in schema.inputs}
    unknown_wired_inputs = wired_input_ids - schema_input_ids
    if unknown_wired_inputs:
        raise ValueError(f"{node_type} has no inputs named {sorted(unknown_wired_inputs)}")

    wire_inputs = [item for item in schema.inputs if _is_wire_input(item, wired_input_ids)]
    widget_inputs = [item for item in schema.inputs if _is_widget_input(item)]
    overrides = widgets or {}
    widget_ids = {item.id for item in widget_inputs}
    widget_ids.update(
        f"{item.id}.{nested.id}"
        for item in widget_inputs
        if item.get_io_type() == "COMFY_DYNAMICCOMBO_V3"
        for option in item.options
        for nested in option.inputs
    )
    unknown = set(overrides) - widget_ids
    if unknown:
        raise ValueError(f"{node_type} has no widget inputs named {sorted(unknown)}")

    widget_entries = [entry for item in widget_inputs for entry in _widget_entries(item, overrides)]
    active_widget_ids = {name for name, _, _ in widget_entries}
    inactive = set(overrides) - active_widget_ids
    if inactive:
        raise ValueError(f"{node_type} widget inputs are inactive for the selected dynamic option: {sorted(inactive)}")
    if size is None:
        visible_widget_count = sum(
            not advanced and name.partition(".")[0] not in wired_input_ids for name, _, advanced in widget_entries
        )
        size = (320, max(80, 52 + 24 * visible_widget_count))
    if node_type == "OpenBioSingleCellPreviewResult":
        size = (300, 220)

    return {
        "alias": alias,
        "type": node_type,
        "pos": list(pos),
        "size": list(size),
        "wire_inputs": [(item.id, item.get_io_type(), _is_widget_input(item)) for item in wire_inputs],
        "outputs": [(output.display_name, output.io_type) for output in schema.outputs],
        "widgets": [(name, value) for name, value, _ in widget_entries],
    }


def _group(group_id: int, title: str, bounding: tuple[int, int, int, int]) -> dict[str, Any]:
    return {
        "id": group_id,
        "title": title,
        "bounding": list(bounding),
        "color": GROUP_COLORS[(group_id - 1) % len(GROUP_COLORS)],
        "font_size": 24,
        "flags": {},
    }


def _workflow(
    stem: str,
    node_specs: list[dict[str, Any]],
    connections: list[tuple[str, str, str, str]],
    groups: list[dict[str, Any]],
    *,
    scale: float,
    offset: tuple[int, int],
) -> dict[str, Any]:
    aliases = {spec["alias"]: index for index, spec in enumerate(node_specs, start=1)}
    if len(aliases) != len(node_specs):
        raise ValueError(f"{stem} contains duplicate node aliases")

    nodes = []
    for order, spec in enumerate(node_specs):
        nodes.append(
            {
                "id": order + 1,
                "type": spec["type"],
                "pos": spec["pos"],
                "size": spec["size"],
                "flags": {},
                "order": order,
                "mode": 0,
                "inputs": [
                    {
                        "name": name,
                        "type": wire_type,
                        **({"widget": {"name": name}} if is_widget else {}),
                        "link": None,
                    }
                    for name, wire_type, is_widget in spec["wire_inputs"]
                ],
                "outputs": [{"name": name, "type": wire_type, "links": []} for name, wire_type in spec["outputs"]],
                "properties": {"Node name for S&R": spec["type"]},
                "widgets_values": [value for _, value in spec["widgets"]],
                "widgets_values_named": dict(spec["widgets"]),
            }
        )

    links = []
    for link_id, (source_alias, source_output, target_alias, target_input) in enumerate(connections, start=1):
        source = nodes[aliases[source_alias] - 1]
        target = nodes[aliases[target_alias] - 1]
        source_slot = next(index for index, item in enumerate(source["outputs"]) if item["name"] == source_output)
        target_slot = next(index for index, item in enumerate(target["inputs"]) if item["name"] == target_input)
        source_type = source["outputs"][source_slot]["type"]
        accepted_types = target["inputs"][target_slot]["type"].split(",")
        if source_type not in accepted_types:
            raise ValueError(
                f"Cannot connect {source_alias}.{source_output} ({source_type}) "
                f"to {target_alias}.{target_input} ({accepted_types})"
            )
        if target["inputs"][target_slot]["link"] is not None:
            raise ValueError(f"{target_alias}.{target_input} is already connected")

        source["outputs"][source_slot]["links"].append(link_id)
        target["inputs"][target_slot]["link"] = link_id
        links.append(
            {
                "id": link_id,
                "origin_id": source["id"],
                "origin_slot": source_slot,
                "target_id": target["id"],
                "target_slot": target_slot,
                "type": source_type,
            }
        )

    for node in nodes:
        disconnected = [item["name"] for item in node["inputs"] if item["link"] is None]
        if disconnected:
            raise ValueError(f"{stem}: {node['type']} has disconnected inputs {disconnected}")

    return {
        "id": str(uuid.uuid5(uuid.NAMESPACE_URL, f"openbio-singlecell:{stem}")),
        "revision": 0,
        "version": 1,
        "state": {
            "lastGroupId": max((group["id"] for group in groups), default=0),
            "lastNodeId": len(nodes),
            "lastLinkId": len(links),
            "lastRerouteId": 0,
        },
        "config": {},
        "groups": groups,
        "nodes": nodes,
        "links": links,
        "extra": {
            "ds": {"scale": scale, "offset": list(offset)},
            "frontendVersion": "1.52.3",
        },
    }


def _quality_control() -> tuple[str, dict[str, Any]]:
    stem = "Quality Control and Clean Counts"
    nodes = [
        _node("load", "OpenBioSingleCellLoadH5AD", (40, 40), {"path": DEMO_PATH}),
        _node("qc", "OpenBioSingleCellCalculateQC", (400, 40)),
        _node("before_plot", "OpenBioSingleCellQCPlots", (800, -360)),
        _node("before_preview", "OpenBioSingleCellPreviewResult", (1160, -430)),
        _node(
            "before_save",
            "OpenBioSingleCellSavePNG",
            (1160, -150),
            {"filename_prefix": "qc_before_filtering"},
        ),
        _node(
            "filter_cells",
            "OpenBioSingleCellFilterCells",
            (800, 180),
            {
                "min_genes": 90,
                "max_genes": 0,
                "min_counts": 0,
                "max_counts": 0,
                "max_pct_mito": 20.0,
                "mito_column": "pct_counts_mt",
            },
        ),
        _node("filter_genes", "OpenBioSingleCellFilterGenes", (1160, 180), {"min_cells": 3}),
        _node("after_plot", "OpenBioSingleCellQCPlots", (1540, -360)),
        _node("after_preview", "OpenBioSingleCellPreviewResult", (1900, -430)),
        _node(
            "after_save",
            "OpenBioSingleCellSavePNG",
            (1900, -150),
            {"filename_prefix": "qc_retained_cells"},
        ),
        _node("summary", "OpenBioSingleCellAnnDataSummary", (1540, 80)),
        _node("summary_preview", "OpenBioSingleCellPreviewResult", (1900, 30)),
        _node(
            "save_adata",
            "OpenBioSingleCellSaveH5AD",
            (1540, 340),
            {"filename_prefix": "qc_clean_counts"},
        ),
    ]
    connections = [
        ("load", "adata", "qc", "adata"),
        ("qc", "adata", "before_plot", "adata"),
        ("before_plot", "plot", "before_preview", "result"),
        ("before_plot", "plot", "before_save", "plot"),
        ("qc", "adata", "filter_cells", "adata"),
        ("filter_cells", "adata", "filter_genes", "adata"),
        ("filter_cells", "adata", "after_plot", "adata"),
        ("after_plot", "plot", "after_preview", "result"),
        ("after_plot", "plot", "after_save", "plot"),
        ("filter_genes", "adata", "summary", "adata"),
        ("summary", "summary", "summary_preview", "result"),
        ("filter_genes", "adata", "save_adata", "adata"),
    ]
    groups = [
        _group(1, "1 · Load and calculate QC metrics", (0, -20, 760, 220)),
        _group(2, "2 · Inspect counts before filtering", (760, -500, 740, 460)),
        _group(3, "3 · Apply QC filters", (760, 110, 740, 300)),
        _group(4, "4 · Review and export clean counts", (1500, -500, 740, 950)),
    ]
    return stem, _workflow(stem, nodes, connections, groups, scale=0.62, offset=(60, 430))


def _cell_clustering() -> tuple[str, dict[str, Any]]:
    stem = "Cell Clustering and Marker Discovery"
    nodes = [
        _node("load", "OpenBioSingleCellLoadH5AD", (40, 70), {"path": DEMO_PATH}),
        _node(
            "normalize",
            "OpenBioSingleCellNormalizeToLayer",
            (400, 70),
            {
                "source": "X",
                "target_sum": 10000.0,
                "transform": "log1p",
                "output_layer": "log1p_norm",
            },
        ),
        _node(
            "hvg",
            "OpenBioSingleCellHighlyVariableGenes",
            (760, 40),
            {
                "n_top_genes": 200,
                "flavor": "seurat",
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "batch_key": "",
                "always_keep_genes": "",
                "subset": False,
            },
        ),
        _node(
            "pca",
            "OpenBioSingleCellPCA",
            (1140, 70),
            {
                "n_comps": 30,
                "use_hvg": True,
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "random_seed": 0,
            },
        ),
        _node(
            "neighbors",
            "OpenBioSingleCellNeighbors",
            (1500, 40),
            {
                "n_neighbors": 15,
                "n_pcs": 30,
                "metric": "cosine",
                "use_rep": "X_pca",
                "key_added": "",
                "random_seed": 0,
            },
        ),
        _node("umap", "OpenBioSingleCellUMAP", (1860, 70)),
        _node("leiden", "OpenBioSingleCellLeiden", (2220, 70)),
        _node(
            "markers",
            "OpenBioSingleCellMarkerGenes",
            (2600, -360),
            {
                "groupby": "leiden",
                "method": "wilcoxon",
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "n_genes": 100,
                "pts": True,
                "random_seed": 0,
            },
        ),
        _node(
            "markers_all_csv",
            "OpenBioSingleCellExportCSV",
            (2960, -450),
            {"filename_prefix": "cluster_markers_all"},
        ),
        _node("marker_filter", "OpenBioSingleCellFilterMarkerGenes", (2960, -180)),
        _node("markers_preview", "OpenBioSingleCellPreviewResult", (3320, -260)),
        _node(
            "markers_filtered_csv",
            "OpenBioSingleCellExportCSV",
            (3320, 20),
            {"filename_prefix": "cluster_markers_filtered"},
        ),
        _node(
            "umap_plot",
            "OpenBioSingleCellUMAPPlot",
            (2600, 220),
            {"color": "leiden", "point_size": 10.0, "color_map": "viridis"},
        ),
        _node("umap_preview", "OpenBioSingleCellPreviewResult", (2960, 150)),
        _node(
            "umap_png",
            "OpenBioSingleCellSavePNG",
            (2960, 430),
            {"filename_prefix": "clusters_umap"},
        ),
        _node("summary", "OpenBioSingleCellAnnDataSummary", (2600, 590)),
        _node("summary_preview", "OpenBioSingleCellPreviewResult", (2960, 560)),
        _node(
            "save_adata",
            "OpenBioSingleCellSaveH5AD",
            (2600, 800),
            {"filename_prefix": "clustered_adata"},
        ),
    ]
    connections = [
        ("load", "adata", "normalize", "adata"),
        ("normalize", "adata", "hvg", "adata"),
        ("hvg", "adata", "pca", "adata"),
        ("pca", "adata", "neighbors", "adata"),
        ("neighbors", "adata", "umap", "adata"),
        ("umap", "adata", "leiden", "adata"),
        ("leiden", "adata", "markers", "adata"),
        ("markers", "table", "markers_all_csv", "table"),
        ("markers", "table", "marker_filter", "table"),
        ("marker_filter", "table", "markers_preview", "result"),
        ("marker_filter", "table", "markers_filtered_csv", "table"),
        ("leiden", "adata", "umap_plot", "adata"),
        ("umap_plot", "plot", "umap_preview", "result"),
        ("umap_plot", "plot", "umap_png", "plot"),
        ("leiden", "adata", "summary", "adata"),
        ("summary", "summary", "summary_preview", "result"),
        ("leiden", "adata", "save_adata", "adata"),
    ]
    groups = [
        _group(1, "1 · Load raw integer counts in adata.X", (0, 0, 360, 250)),
        _group(2, "2 · Normalize to a layer and select variable genes", (360, -20, 760, 330)),
        _group(3, "3 · Embed and cluster cells", (1100, -20, 1480, 330)),
        _group(4, "4 · Discover and filter cluster markers", (2560, -520, 1100, 660)),
        _group(5, "5 · Review and export", (2560, 150, 1100, 780)),
    ]
    return stem, _workflow(stem, nodes, connections, groups, scale=0.48, offset=(40, 520))


def _sample_composition() -> tuple[str, dict[str, Any]]:
    stem = "Sample Composition Comparison"
    nodes = [
        _node(
            "study_parameters",
            "OpenBioSingleCellCoreStudyParameters",
            (40, -400),
            {"study_parameters_json": DEMO_STUDY_PARAMETERS_JSON},
            size=(440, 340),
        ),
        _node("load", "OpenBioSingleCellLoadH5AD", (560, -260), {"path": DEMO_PATH}),
        _node(
            "composition",
            "OpenBioSingleCellSampleCompositionSummary",
            (1100, -300),
            wired_inputs=("sample_key", "group_key", "annotation_key"),
        ),
        _node("composition_preview", "OpenBioSingleCellPreviewResult", (1480, -390)),
        _node(
            "composition_csv",
            "OpenBioSingleCellExportCSV",
            (1480, -110),
            {"filename_prefix": "sample_composition"},
        ),
        _node(
            "contrast",
            "OpenBioSingleCellDifferentialCompositionTest",
            (1100, 240),
            wired_inputs=(
                "sample_key",
                "group_key",
                "annotation_key",
                "control_group",
                "comparison_groups",
            ),
        ),
        _node("contrast_preview", "OpenBioSingleCellPreviewResult", (1480, 150)),
        _node(
            "contrast_csv",
            "OpenBioSingleCellExportCSV",
            (1480, 430),
            {"filename_prefix": "differential_composition"},
        ),
        _node("summary", "OpenBioSingleCellAnnDataSummary", (1100, 750)),
        _node("summary_preview", "OpenBioSingleCellPreviewResult", (1480, 690)),
    ]
    connections = [
        ("load", "adata", "composition", "adata"),
        ("study_parameters", "sample_column", "composition", "sample_key"),
        ("study_parameters", "condition_column", "composition", "group_key"),
        ("study_parameters", "annotation_column", "composition", "annotation_key"),
        ("composition", "table", "composition_preview", "result"),
        ("composition", "table", "composition_csv", "table"),
        ("load", "adata", "contrast", "adata"),
        ("study_parameters", "sample_column", "contrast", "sample_key"),
        ("study_parameters", "condition_column", "contrast", "group_key"),
        ("study_parameters", "annotation_column", "contrast", "annotation_key"),
        ("study_parameters", "reference", "contrast", "control_group"),
        ("study_parameters", "comparison", "contrast", "comparison_groups"),
        ("contrast", "table", "contrast_preview", "result"),
        ("contrast", "table", "contrast_csv", "table"),
        ("load", "adata", "summary", "adata"),
        ("summary", "summary", "summary_preview", "result"),
    ]
    groups = [
        _group(1, "1 · Select reusable study parameters", (0, -470, 500, 410)),
        _group(2, "2 · Load AnnData with sample metadata", (520, -330, 420, 300)),
        _group(3, "3 · Summarize and compare sample composition", (1060, -450, 780, 1080)),
        _group(4, "4 · Inspect the input AnnData", (1060, 670, 780, 370)),
    ]
    return stem, _workflow(stem, nodes, connections, groups, scale=0.58, offset=(80, 500))


def _scvi_integration() -> tuple[str, dict[str, Any]]:
    stem = "scVI Batch Integration and Contrast"
    nodes = [
        _node("load", "OpenBioSingleCellLoadH5AD", (40, 100), {"path": DEMO_PATH}),
        _node("qc", "OpenBioSingleCellCalculateQC", (400, 80)),
        _node(
            "filter_cells",
            "OpenBioSingleCellFilterCells",
            (760, 50),
            {
                "min_genes": 90,
                "max_genes": 0,
                "min_counts": 0,
                "max_counts": 0,
                "max_pct_mito": 20.0,
                "mito_column": "pct_counts_mt",
            },
        ),
        _node("filter_genes", "OpenBioSingleCellFilterGenes", (1120, 60), {"min_cells": 3}),
        _node(
            "normalize",
            "OpenBioSingleCellNormalizeToLayer",
            (1500, 60),
            {
                "source": "X",
                "target_sum": 10000.0,
                "transform": "log1p",
                "output_layer": "log1p_norm",
            },
        ),
        _node(
            "scvi",
            "OpenBioSingleCellSCVIIntegration",
            (1840, 50),
            {"source": "X", "batch_key": "batch"},
        ),
        _node(
            "neighbors",
            "OpenBioSingleCellNeighbors",
            (2220, 50),
            {
                "n_neighbors": 15,
                "n_pcs": 10,
                "metric": "cosine",
                "use_rep": "X_scVI",
                "key_added": "",
                "random_seed": 0,
            },
        ),
        _node("umap", "OpenBioSingleCellUMAP", (2580, 70)),
        _node("leiden", "OpenBioSingleCellLeiden", (2940, 70)),
        _node(
            "de",
            "OpenBioSingleCellSCVIDifferentialExpression",
            (2220, -450),
            {
                "groupby": "cell_type",
                "group1": "T cell",
                "group2": "B cell",
                "subset_column": "",
                "subset_value": "",
                "mode": "vanilla",
                "delta": 0.25,
            },
        ),
        _node("de_preview", "OpenBioSingleCellPreviewResult", (2580, -540)),
        _node(
            "de_csv",
            "OpenBioSingleCellExportCSV",
            (2580, -260),
            {"filename_prefix": "scvi_differential_expression"},
        ),
        _node(
            "batch_plot",
            "OpenBioSingleCellUMAPPlot",
            (3320, -330),
            {"color": "batch", "point_size": 10.0, "color_map": "viridis"},
        ),
        _node("batch_preview", "OpenBioSingleCellPreviewResult", (3680, -420)),
        _node(
            "batch_png",
            "OpenBioSingleCellSavePNG",
            (3680, -140),
            {"filename_prefix": "scvi_umap_batch"},
        ),
        _node(
            "cell_type_plot",
            "OpenBioSingleCellUMAPPlot",
            (3320, 100),
            {"color": "cell_type", "point_size": 10.0, "color_map": "viridis"},
        ),
        _node("cell_type_preview", "OpenBioSingleCellPreviewResult", (3680, 20)),
        _node(
            "cell_type_png",
            "OpenBioSingleCellSavePNG",
            (3680, 300),
            {"filename_prefix": "scvi_umap_cell_type"},
        ),
        _node(
            "save_adata",
            "OpenBioSingleCellSaveH5AD",
            (3320, 520),
            {"filename_prefix": "scvi_integrated_adata"},
        ),
    ]
    connections = [
        ("load", "adata", "qc", "adata"),
        ("qc", "adata", "filter_cells", "adata"),
        ("filter_cells", "adata", "filter_genes", "adata"),
        ("filter_genes", "adata", "normalize", "adata"),
        ("normalize", "adata", "scvi", "adata"),
        ("scvi", "adata", "neighbors", "adata"),
        ("neighbors", "adata", "umap", "adata"),
        ("umap", "adata", "leiden", "adata"),
        ("scvi", "adata", "de", "adata"),
        ("scvi", "model", "de", "model"),
        ("de", "table", "de_preview", "result"),
        ("de", "table", "de_csv", "table"),
        ("leiden", "adata", "batch_plot", "adata"),
        ("batch_plot", "plot", "batch_preview", "result"),
        ("batch_plot", "plot", "batch_png", "plot"),
        ("leiden", "adata", "cell_type_plot", "adata"),
        ("cell_type_plot", "plot", "cell_type_preview", "result"),
        ("cell_type_plot", "plot", "cell_type_png", "plot"),
        ("leiden", "adata", "save_adata", "adata"),
    ]
    groups = [
        _group(1, "1 · Load and clean counts", (0, 0, 1460, 350)),
        _group(2, "2 · Normalize without replacing raw adata.X", (1460, 0, 360, 350)),
        _group(3, "3 · Train scVI from raw integer counts", (1800, -20, 400, 370)),
        _group(4, "4 · Build integrated neighborhoods", (2180, 0, 1100, 350)),
        _group(5, "5 · Contrast cell types with the trained model", (2180, -590, 760, 540)),
        _group(6, "6 · Review and export integrated results", (3280, -450, 760, 1100)),
    ]
    return stem, _workflow(stem, nodes, connections, groups, scale=0.44, offset=(40, 520))


def _single_cell_best_practice() -> tuple[str, dict[str, Any]]:
    stem = "Single-Cell Best Practice"
    nodes = [
        _node(
            "study_parameters",
            "OpenBioSingleCellCoreStudyParameters",
            (40, -470),
            {"study_parameters_json": BEST_PRACTICE_STUDY_PARAMETERS_JSON},
            size=(440, 340),
        ),
        _node("load", "OpenBioSingleCellLoadH5AD", (40, 0), {"path": BEST_PRACTICE_PATH}),
        _node("qc_initial", "OpenBioSingleCellCalculateQC", (520, 0)),
        _node("qc_before_plot", "OpenBioSingleCellQCPlots", (900, -450)),
        _node("qc_before_preview", "OpenBioSingleCellPreviewResult", (1260, -520)),
        _node(
            "qc_before_png",
            "OpenBioSingleCellSavePNG",
            (1260, -240),
            {"filename_prefix": "best_practice_qc_before_filtering"},
        ),
        _node(
            "hard_filter_cells",
            "OpenBioSingleCellFilterCells",
            (900, 0),
            {
                "min_genes": 200,
                "max_genes": 0,
                "min_counts": 0,
                "max_counts": 0,
                "max_pct_mito": 0.0,
                "mito_column": "pct_counts_mt",
            },
        ),
        _node(
            "hard_filter_genes",
            "OpenBioSingleCellFilterGenes",
            (1260, 0),
            {"min_cells": 3},
        ),
        _node("qc_refreshed", "OpenBioSingleCellCalculateQC", (1620, 0)),
        _node(
            "scrublet",
            "OpenBioSingleCellScrublet",
            (1980, 0),
            {"batch_key": "sample", "random_seed": 0},
            wired_inputs=("batch_key",),
        ),
        _node("filter_doublets", "OpenBioSingleCellFilterDoublets", (2340, 0)),
        _node(
            "mad_general",
            "OpenBioSingleCellMarkMADOutliers",
            (2700, -110),
            {
                "metrics": "total_counts,n_genes_by_counts,pct_counts_in_top_20_genes",
                "batch_key": "sample",
                "nmads": 5.0,
                "direction": "both",
                "output_column": "qc_mad_general",
                "scale_mad": False,
            },
            wired_inputs=("batch_key",),
        ),
        _node(
            "mad_mito",
            "OpenBioSingleCellMarkMADOutliers",
            (3060, -110),
            {
                "metrics": "pct_counts_mt",
                "batch_key": "sample",
                "nmads": 3.0,
                "direction": "upper",
                "output_column": "qc_mad_mito",
                "scale_mad": False,
            },
            wired_inputs=("batch_key",),
        ),
        _node(
            "keep_general",
            "OpenBioSingleCellSubsetObservations",
            (3420, -180),
            {"column": "qc_mad_general", "values": "True", "invert": True},
        ),
        _node(
            "keep_mito",
            "OpenBioSingleCellSubsetObservations",
            (3420, 130),
            {"column": "qc_mad_mito", "values": "True", "invert": True},
        ),
        _node(
            "mito_cap",
            "OpenBioSingleCellFilterCells",
            (3780, 0),
            {
                "min_genes": 0,
                "max_genes": 0,
                "min_counts": 0,
                "max_counts": 0,
                "max_pct_mito": 20.0,
                "mito_column": "pct_counts_mt",
            },
        ),
        _node("qc_after_plot", "OpenBioSingleCellQCPlots", (4140, -450)),
        _node("qc_after_preview", "OpenBioSingleCellPreviewResult", (4500, -520)),
        _node(
            "qc_after_png",
            "OpenBioSingleCellSavePNG",
            (4500, -240),
            {"filename_prefix": "best_practice_qc_retained_singlets"},
        ),
        _node("qc_summary", "OpenBioSingleCellAnnDataSummary", (4140, 100)),
        _node("qc_summary_preview", "OpenBioSingleCellPreviewResult", (4500, 70)),
        _node(
            "snapshot",
            "OpenBioSingleCellSnapshotExpression",
            (4860, 0),
            {"destination": "layer_and_raw", "layer_name": "counts"},
        ),
        _node(
            "celltypist",
            "OpenBioSingleCellCellTypistAnnotation",
            (5220, 0),
            {
                "model": "Adult_Human_Vascular.pkl",
                "use_raw": True,
                "majority_voting": True,
                "label_column": "celltypist_cell_type",
                "confidence_column": "celltypist_confidence",
            },
        ),
        _node(
            "normalize",
            "OpenBioSingleCellNormalizeToLayer",
            (5580, 0),
            {
                "source": "layer",
                "source.source_layer": "counts",
                "target_sum": 10000.0,
                "transform": "log1p",
                "output_layer": "log1p_norm",
            },
        ),
        _node(
            "hvg_full",
            "OpenBioSingleCellHighlyVariableGenes",
            (5940, -120),
            {
                "n_top_genes": 5000,
                "flavor": "seurat",
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "batch_key": "sample",
                "always_keep_genes": "",
                "subset": False,
            },
            wired_inputs=("batch_key",),
        ),
        _node(
            "hvg_subset",
            "OpenBioSingleCellHighlyVariableGenes",
            (6300, 0),
            {
                "n_top_genes": 5000,
                "flavor": "seurat",
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "batch_key": "sample",
                "always_keep_genes": "",
                "subset": True,
            },
            wired_inputs=("batch_key",),
        ),
        _node(
            "scvi",
            "OpenBioSingleCellSCVIIntegration",
            (6660, 0),
            {
                "source": "layer",
                "source.counts_layer": "counts",
                "batch_key": "batch",
                "n_latent": 30,
                "gene_likelihood": "nb",
                "max_epochs": 0,
                "early_stopping": True,
                "output_key": "X_scVI",
                "random_seed": 0,
            },
            wired_inputs=("batch_key",),
        ),
        _node(
            "neighbors",
            "OpenBioSingleCellNeighbors",
            (7020, 0),
            {
                "n_neighbors": 15,
                "n_pcs": 30,
                "metric": "cosine",
                "use_rep": "X_scVI",
                "key_added": "",
                "random_seed": 0,
            },
        ),
        _node(
            "umap",
            "OpenBioSingleCellUMAP",
            (7380, 0),
            {"min_dist": 0.3, "spread": 1.0, "neighbors_key": "neighbors", "random_seed": 0},
        ),
        _node(
            "leiden",
            "OpenBioSingleCellLeiden",
            (7740, 0),
            {"resolution": 1.0, "key_added": "leiden_scvi", "neighbors_key": "neighbors", "random_seed": 0},
        ),
        _node(
            "hvg_markers",
            "OpenBioSingleCellMarkerGenes",
            (8100, -760),
            {
                "groupby": "leiden_scvi",
                "method": "wilcoxon",
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "n_genes": 50,
                "pts": True,
                "random_seed": 0,
            },
        ),
        _node("hvg_marker_filter", "OpenBioSingleCellFilterMarkerGenes", (8460, -760)),
        _node("hvg_markers_preview", "OpenBioSingleCellPreviewResult", (8820, -850)),
        _node(
            "hvg_markers_csv",
            "OpenBioSingleCellExportCSV",
            (8820, -570),
            {"filename_prefix": "best_practice_hvg_cluster_markers"},
        ),
        _node(
            "umap_batch",
            "OpenBioSingleCellUMAPPlot",
            (8100, -330),
            {"color": "batch", "point_size": 3.0, "color_map": "viridis"},
        ),
        _node("umap_batch_preview", "OpenBioSingleCellPreviewResult", (8460, -410)),
        _node(
            "umap_batch_png",
            "OpenBioSingleCellSavePNG",
            (8460, -130),
            {"filename_prefix": "best_practice_umap_batch"},
        ),
        _node(
            "umap_condition",
            "OpenBioSingleCellUMAPPlot",
            (8100, 90),
            {"color": "group", "point_size": 3.0, "color_map": "viridis"},
        ),
        _node("umap_condition_preview", "OpenBioSingleCellPreviewResult", (8460, 10)),
        _node(
            "umap_condition_png",
            "OpenBioSingleCellSavePNG",
            (8460, 290),
            {"filename_prefix": "best_practice_umap_condition"},
        ),
        _node(
            "umap_clusters",
            "OpenBioSingleCellUMAPPlot",
            (8100, 510),
            {"color": "leiden_scvi", "point_size": 3.0, "color_map": "viridis"},
        ),
        _node("umap_clusters_preview", "OpenBioSingleCellPreviewResult", (8460, 430)),
        _node(
            "umap_clusters_png",
            "OpenBioSingleCellSavePNG",
            (8460, 710),
            {"filename_prefix": "best_practice_umap_clusters"},
        ),
        _node(
            "umap_celltypist",
            "OpenBioSingleCellUMAPPlot",
            (8100, 930),
            {"color": "celltypist_cell_type", "point_size": 3.0, "color_map": "viridis"},
        ),
        _node("umap_celltypist_preview", "OpenBioSingleCellPreviewResult", (8460, 850)),
        _node(
            "umap_celltypist_png",
            "OpenBioSingleCellSavePNG",
            (8460, 1130),
            {"filename_prefix": "best_practice_umap_celltypist"},
        ),
        _node(
            "merge_clusters",
            "OpenBioSingleCellMergeObservationAnnotations",
            (9180, 0),
            {
                "source_column": "leiden_scvi",
                "target_column": "leiden_scvi",
                "prefix": "",
                "overwrite": False,
            },
        ),
        _node(
            "full_markers",
            "OpenBioSingleCellMarkerGenes",
            (9540, -640),
            {
                "groupby": "leiden_scvi",
                "method": "wilcoxon",
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "n_genes": 100,
                "pts": True,
                "random_seed": 0,
            },
        ),
        _node("full_marker_filter", "OpenBioSingleCellFilterMarkerGenes", (9900, -640)),
        _node("full_markers_preview", "OpenBioSingleCellPreviewResult", (10260, -730)),
        _node(
            "full_markers_csv",
            "OpenBioSingleCellExportCSV",
            (10260, -450),
            {"filename_prefix": "best_practice_full_gene_cluster_markers"},
        ),
        _node(
            "marker_panel",
            "OpenBioSingleCellMarkerExpressionPlot",
            (9540, -80),
            {
                "genes": BEST_PRACTICE_MARKER_PANEL,
                "groupby": "leiden_scvi",
                "plot_type": "dotplot",
                "source": "layer",
                "source.layer_name": "log1p_norm",
                "standard_scale": "var",
                "dendrogram": False,
                "log": False,
            },
        ),
        _node("marker_panel_preview", "OpenBioSingleCellPreviewResult", (9900, -150)),
        _node(
            "marker_panel_png",
            "OpenBioSingleCellSavePNG",
            (9900, 130),
            {"filename_prefix": "best_practice_reference_marker_panel"},
        ),
        _node(
            "composition",
            "OpenBioSingleCellSampleCompositionSummary",
            (9540, 390),
            wired_inputs=("sample_key", "group_key", "annotation_key"),
        ),
        _node("composition_preview", "OpenBioSingleCellPreviewResult", (9900, 330)),
        _node(
            "composition_csv",
            "OpenBioSingleCellExportCSV",
            (9900, 610),
            {"filename_prefix": "best_practice_celltypist_composition"},
        ),
        _node(
            "contrast_conditions",
            "OpenBioSingleCellSubsetObservations",
            (10620, 0),
            {"column": "group", "values": "Normal,nonDM_ED", "invert": False},
            wired_inputs=("column",),
        ),
        _node(
            "pseudobulk",
            "OpenBioSingleCellPseudobulk",
            (10980, 0),
            {
                "sample_key": "sample",
                "groupby": "celltypist_cell_type",
                "source": "layer",
                "source.layer_name": "counts",
                "mode": "sum",
                "min_cells": 10,
                "min_counts": 1000,
            },
            wired_inputs=("sample_key", "groupby"),
        ),
        _node(
            "selected_population",
            "OpenBioSingleCellSubsetObservations",
            (11340, 0),
            {"column": "celltypist_cell_type", "values": "smc_pc_intermediate", "invert": False},
            wired_inputs=("column",),
        ),
        _node(
            "deseq2",
            "OpenBioSingleCellPseudobulkDESeq2",
            (11700, 0),
            {
                "design": "~group",
                "contrast_column": "group",
                "baseline": "Normal",
                "comparison": "nonDM_ED",
            },
            wired_inputs=("contrast_column", "baseline", "comparison"),
        ),
        _node("deseq2_preview", "OpenBioSingleCellPreviewResult", (12060, -100)),
        _node(
            "deseq2_csv",
            "OpenBioSingleCellExportCSV",
            (12060, 180),
            {"filename_prefix": "best_practice_smc_pc_intermediate_nonDM_ED_vs_Normal"},
        ),
        _node("full_summary", "OpenBioSingleCellAnnDataSummary", (12780, -340)),
        _node("full_summary_preview", "OpenBioSingleCellPreviewResult", (13140, -430)),
        _node(
            "save_full",
            "OpenBioSingleCellSaveH5AD",
            (13140, -150),
            {"filename_prefix": "best_practice_full_gene_results"},
        ),
        _node("integrated_summary", "OpenBioSingleCellAnnDataSummary", (12780, 300)),
        _node("integrated_summary_preview", "OpenBioSingleCellPreviewResult", (13140, 210)),
        _node(
            "save_integrated",
            "OpenBioSingleCellSaveH5AD",
            (13140, 490),
            {"filename_prefix": "best_practice_hvg_scvi_results"},
        ),
    ]
    connections = [
        ("load", "adata", "qc_initial", "adata"),
        ("qc_initial", "adata", "qc_before_plot", "adata"),
        ("qc_before_plot", "plot", "qc_before_preview", "result"),
        ("qc_before_plot", "plot", "qc_before_png", "plot"),
        ("qc_initial", "adata", "hard_filter_cells", "adata"),
        ("hard_filter_cells", "adata", "hard_filter_genes", "adata"),
        ("hard_filter_genes", "adata", "qc_refreshed", "adata"),
        ("qc_refreshed", "adata", "scrublet", "adata"),
        ("study_parameters", "sample_column", "scrublet", "batch_key"),
        ("scrublet", "adata", "filter_doublets", "adata"),
        ("filter_doublets", "adata", "mad_general", "adata"),
        ("study_parameters", "sample_column", "mad_general", "batch_key"),
        ("mad_general", "adata", "mad_mito", "adata"),
        ("study_parameters", "sample_column", "mad_mito", "batch_key"),
        ("mad_mito", "adata", "keep_general", "adata"),
        ("keep_general", "adata", "keep_mito", "adata"),
        ("keep_mito", "adata", "mito_cap", "adata"),
        ("mito_cap", "adata", "qc_after_plot", "adata"),
        ("qc_after_plot", "plot", "qc_after_preview", "result"),
        ("qc_after_plot", "plot", "qc_after_png", "plot"),
        ("mito_cap", "adata", "qc_summary", "adata"),
        ("qc_summary", "summary", "qc_summary_preview", "result"),
        ("mito_cap", "adata", "snapshot", "adata"),
        ("snapshot", "adata", "celltypist", "adata"),
        ("celltypist", "adata", "normalize", "adata"),
        ("normalize", "adata", "hvg_full", "adata"),
        ("study_parameters", "sample_column", "hvg_full", "batch_key"),
        ("hvg_full", "adata", "hvg_subset", "adata"),
        ("study_parameters", "sample_column", "hvg_subset", "batch_key"),
        ("hvg_subset", "adata", "scvi", "adata"),
        ("study_parameters", "batch_column", "scvi", "batch_key"),
        ("scvi", "adata", "neighbors", "adata"),
        ("neighbors", "adata", "umap", "adata"),
        ("umap", "adata", "leiden", "adata"),
        ("leiden", "adata", "hvg_markers", "adata"),
        ("hvg_markers", "table", "hvg_marker_filter", "table"),
        ("hvg_marker_filter", "table", "hvg_markers_preview", "result"),
        ("hvg_marker_filter", "table", "hvg_markers_csv", "table"),
        ("leiden", "adata", "umap_batch", "adata"),
        ("umap_batch", "plot", "umap_batch_preview", "result"),
        ("umap_batch", "plot", "umap_batch_png", "plot"),
        ("leiden", "adata", "umap_condition", "adata"),
        ("umap_condition", "plot", "umap_condition_preview", "result"),
        ("umap_condition", "plot", "umap_condition_png", "plot"),
        ("leiden", "adata", "umap_clusters", "adata"),
        ("umap_clusters", "plot", "umap_clusters_preview", "result"),
        ("umap_clusters", "plot", "umap_clusters_png", "plot"),
        ("leiden", "adata", "umap_celltypist", "adata"),
        ("umap_celltypist", "plot", "umap_celltypist_preview", "result"),
        ("umap_celltypist", "plot", "umap_celltypist_png", "plot"),
        ("hvg_full", "adata", "merge_clusters", "adata"),
        ("leiden", "adata", "merge_clusters", "subset_adata"),
        ("merge_clusters", "adata", "full_markers", "adata"),
        ("full_markers", "table", "full_marker_filter", "table"),
        ("full_marker_filter", "table", "full_markers_preview", "result"),
        ("full_marker_filter", "table", "full_markers_csv", "table"),
        ("merge_clusters", "adata", "marker_panel", "adata"),
        ("marker_panel", "plot", "marker_panel_preview", "result"),
        ("marker_panel", "plot", "marker_panel_png", "plot"),
        ("merge_clusters", "adata", "composition", "adata"),
        ("study_parameters", "sample_column", "composition", "sample_key"),
        ("study_parameters", "condition_column", "composition", "group_key"),
        ("study_parameters", "annotation_column", "composition", "annotation_key"),
        ("composition", "table", "composition_preview", "result"),
        ("composition", "table", "composition_csv", "table"),
        ("merge_clusters", "adata", "contrast_conditions", "adata"),
        ("study_parameters", "condition_column", "contrast_conditions", "column"),
        ("contrast_conditions", "adata", "pseudobulk", "adata"),
        ("study_parameters", "sample_column", "pseudobulk", "sample_key"),
        ("study_parameters", "annotation_column", "pseudobulk", "groupby"),
        ("pseudobulk", "adata", "selected_population", "adata"),
        ("study_parameters", "annotation_column", "selected_population", "column"),
        ("selected_population", "adata", "deseq2", "adata"),
        ("study_parameters", "condition_column", "deseq2", "contrast_column"),
        ("study_parameters", "reference", "deseq2", "baseline"),
        ("study_parameters", "comparison", "deseq2", "comparison"),
        ("deseq2", "table", "deseq2_preview", "result"),
        ("deseq2", "table", "deseq2_csv", "table"),
        ("merge_clusters", "adata", "full_summary", "adata"),
        ("full_summary", "summary", "full_summary_preview", "result"),
        ("merge_clusters", "adata", "save_full", "adata"),
        ("leiden", "adata", "integrated_summary", "adata"),
        ("integrated_summary", "summary", "integrated_summary_preview", "result"),
        ("leiden", "adata", "save_integrated", "adata"),
    ]
    groups = [
        _group(1, "1 · Study contract and audited input", (0, -560, 860, 900)),
        _group(2, "2 · Hard QC, Scrublet, then sample-wise MAD flags", (860, -600, 3900, 1000)),
        _group(3, "3 · Preserve full counts/raw and add provisional CellTypist labels", (4760, -200, 780, 520)),
        _group(4, "4 · Full-gene log1p plus 5,000 batch-aware Seurat HVGs", (5480, -260, 1140, 600)),
        _group(5, "5 · Train scVI on the HVG count view and cluster X_scVI", (6580, -220, 1480, 560)),
        _group(6, "6 · Integrated diagnostics + HVG Scanpy marker evidence (not scVI DE)", (8040, -920, 1160, 2350)),
        _group(7, "7 · Merge clusters into full genes; review two marker evidence branches", (9160, -820, 1420, 1700)),
        _group(
            8,
            "8 · One provisional population: smc_pc_intermediate · nonDM_ED vs Normal pseudobulk",
            (10580, -220, 1820, 650),
        ),
        _group(9, "9 · Save full-gene and HVG-scVI objects separately", (12680, -520, 820, 1250)),
    ]
    return stem, _workflow(stem, nodes, connections, groups, scale=0.18, offset=(40, 520))


def build_workflows() -> dict[str, dict[str, Any]]:
    workflows = [
        _quality_control(),
        _cell_clustering(),
        _sample_composition(),
        _scvi_integration(),
        _single_cell_best_practice(),
    ]
    return {f"{stem}.json": workflow for stem, workflow in workflows}


def _serialized(workflow: dict[str, Any]) -> str:
    return json.dumps(workflow, ensure_ascii=False, indent=2) + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate schema-checked OpenBio single-cell example workflows.")
    parser.add_argument(
        "--check", action="store_true", help="Fail if committed workflows differ from generated output."
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    workflows = build_workflows()
    expected_paths = {EXAMPLE_DIRECTORY / filename for filename in workflows}

    if args.check:
        errors = []
        for filename, workflow in workflows.items():
            path = EXAMPLE_DIRECTORY / filename
            expected = _serialized(workflow)
            if not path.is_file():
                errors.append(f"missing: {path}")
            elif path.read_text(encoding="utf-8") != expected:
                errors.append(f"stale: {path}")
        unexpected = sorted(set(EXAMPLE_DIRECTORY.glob("*.json")) - expected_paths)
        errors.extend(f"unexpected: {path}" for path in unexpected)
        if errors:
            raise SystemExit("\n".join(errors))
        print(f"Verified {len(workflows)} generated workflows.")
        return 0

    EXAMPLE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    for filename, workflow in workflows.items():
        path = EXAMPLE_DIRECTORY / filename
        path.write_text(_serialized(workflow), encoding="utf-8")
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
