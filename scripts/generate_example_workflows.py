from __future__ import annotations

import argparse
import json
import uuid
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from openbio_singlecell.extension import NODE_CLASSES

PLUGIN_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIRECTORY = PLUGIN_ROOT / "example_workflows"
DEMO_PATH = "openbio-singlecell/openbio_singlecell_demo.h5ad"
GROUP_COLORS = ["#3f789e", "#2f806d", "#76558f", "#9a6138", "#526d82", "#775b52"]

SCHEMAS = {node.GET_SCHEMA().node_id: node.GET_SCHEMA() for node in NODE_CLASSES}
WIRE_TYPES = {output.io_type for schema in SCHEMAS.values() for output in schema.outputs}


def _is_wire_input(item: Any) -> bool:
    return bool(set(item.get_io_type().split(",")) & WIRE_TYPES)


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
) -> dict[str, Any]:
    schema = SCHEMAS[node_type]
    wire_inputs = [item for item in schema.inputs if _is_wire_input(item)]
    widget_inputs = [item for item in schema.inputs if not _is_wire_input(item)]
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
        raise ValueError(
            f"{node_type} widget inputs are inactive for the selected dynamic option: {sorted(inactive)}"
        )
    if size is None:
        visible_widget_count = sum(not advanced for _, _, advanced in widget_entries)
        size = (320, max(80, 52 + 24 * visible_widget_count))
    if node_type == "OpenBioSingleCellPreviewResult":
        size = (300, 220)

    return {
        "alias": alias,
        "type": node_type,
        "pos": list(pos),
        "size": list(size),
        "wire_inputs": [(item.id, item.get_io_type()) for item in wire_inputs],
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
                "inputs": [{"name": name, "type": wire_type, "link": None} for name, wire_type in spec["wire_inputs"]],
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
        _node("load", "OpenBioSingleCellLoadH5AD", (40, 60), {"path": DEMO_PATH}),
        _node(
            "composition",
            "OpenBioSingleCellSampleCompositionSummary",
            (400, -310),
            {"sample_key": "sample", "group_key": "batch", "annotation_key": "cell_type"},
        ),
        _node("composition_preview", "OpenBioSingleCellPreviewResult", (760, -400)),
        _node(
            "composition_csv",
            "OpenBioSingleCellExportCSV",
            (760, -120),
            {"filename_prefix": "sample_composition"},
        ),
        _node(
            "contrast",
            "OpenBioSingleCellDifferentialCompositionTest",
            (400, 160),
            {
                "sample_key": "sample",
                "group_key": "batch",
                "annotation_key": "cell_type",
                "control_group": "batch_1",
                "comparison_groups": "batch_2",
                "pseudocount": 0.001,
            },
        ),
        _node("contrast_preview", "OpenBioSingleCellPreviewResult", (760, 60)),
        _node(
            "contrast_csv",
            "OpenBioSingleCellExportCSV",
            (760, 340),
            {"filename_prefix": "differential_composition"},
        ),
        _node("summary", "OpenBioSingleCellAnnDataSummary", (400, 650)),
        _node("summary_preview", "OpenBioSingleCellPreviewResult", (760, 600)),
    ]
    connections = [
        ("load", "adata", "composition", "adata"),
        ("composition", "table", "composition_preview", "result"),
        ("composition", "table", "composition_csv", "table"),
        ("load", "adata", "contrast", "adata"),
        ("contrast", "table", "contrast_preview", "result"),
        ("contrast", "table", "contrast_csv", "table"),
        ("load", "adata", "summary", "adata"),
        ("summary", "summary", "summary_preview", "result"),
    ]
    groups = [
        _group(1, "1 · Load annotated AnnData", (0, 0, 360, 240)),
        _group(2, "2 · Summarize sample composition", (360, -450, 760, 430)),
        _group(
            3,
            "3 · Replace batch / batch_1 / batch_2 with your study condition",
            (360, 20, 760, 500),
        ),
        _group(4, "4 · Inspect AnnData metadata", (360, 560, 760, 340)),
    ]
    return stem, _workflow(stem, nodes, connections, groups, scale=0.7, offset=(90, 430))


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


def build_workflows() -> dict[str, dict[str, Any]]:
    workflows = [
        _quality_control(),
        _cell_clustering(),
        _sample_composition(),
        _scvi_integration(),
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
