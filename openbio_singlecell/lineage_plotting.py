from __future__ import annotations

import inspect
import textwrap


def _standalone_lineage_scalar(value):
    import math

    import numpy as np
    import pandas as pd

    if value is None or value is pd.NA:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return None
    return str(value)


def _standalone_lineage_table_fingerprint(table):
    import hashlib
    import json

    payload = {
        "index": [_standalone_lineage_scalar(value) for value in table.index.tolist()],
        "columns": [_standalone_lineage_scalar(value) for value in table.columns.tolist()],
        "data": [
            [_standalone_lineage_scalar(value) for value in row] for row in table.itertuples(index=False, name=None)
        ],
    }
    encoded = json.dumps(payload, allow_nan=False, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def lineage_table_fingerprint(table) -> str:
    return _standalone_lineage_table_fingerprint(table)


def _standalone_cassiopeia_lineage_qc_plot(
    characters,
    table,
    *,
    producer,
    _return_details=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Cassiopeia Lineage QC Plot"
    metadata = characters.metadata
    if not isinstance(metadata, Mapping) or metadata.get("artifact_type") != "OPENBIO_CASSIOPEIA_CHARACTERS":
        raise TypeError(f"{operation} requires an audited Cassiopeia character artifact.")
    expected_producer = {"operation", "parameters", "input_cells", "input_genes", "random_seed"}
    if not isinstance(producer, Mapping) or set(producer) != expected_producer:
        raise ValueError(f"{operation} requires exact QC table producer provenance.")
    if producer["operation"] != "prepare_cassiopeia_characters":
        raise ValueError(f"{operation} requires the canonical Cassiopeia character-preparation table.")
    parameters = producer["parameters"]
    if not isinstance(parameters, Mapping):
        raise ValueError(f"{operation} producer parameters are invalid.")
    if parameters.get("characters_fingerprint") != characters.fingerprint:
        raise ValueError(f"{operation} table is not bound to this character artifact.")
    expected_hash = metadata.get("qc_table_sha256")
    if not isinstance(expected_hash, str) or parameters.get("qc_table_sha256") != expected_hash:
        raise ValueError(f"{operation} QC table fingerprint provenance is inconsistent.")
    columns = [
        "lineage_id",
        "status",
        "qc_warnings",
        "unavailability_reasons",
        "input_cells",
        "retained_cells",
        "input_characters",
        "retained_characters",
        "informative_character_fraction",
        "unique_profile_fraction",
        "mean_missing_fraction",
        "mean_uncut_fraction_observed",
        "all_missing_cells",
        "all_uncut_cells",
        "missing_threshold_filtered_cells",
        "uncut_threshold_filtered_cells",
    ]
    if not isinstance(table, pd.DataFrame) or list(table.columns) != columns:
        raise ValueError(f"{operation} table does not match the canonical lineage-QC schema.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError(f"{operation} table must use the canonical zero-based row index.")
    if table.empty or len(table) > 500:
        raise ValueError(f"{operation} supports 1 through 500 lineage rows in one static plot.")
    if _standalone_lineage_table_fingerprint(table) != expected_hash:
        raise ValueError(f"{operation} QC table content was changed after production.")

    lineage_ids = table["lineage_id"].tolist()
    if (
        any(not isinstance(value, str) or not value or value != value.strip() for value in lineage_ids)
        or lineage_ids != sorted(lineage_ids)
        or len(lineage_ids) != len(set(lineage_ids))
    ):
        raise ValueError(f"{operation} lineage identity must be unique canonical sorted text.")
    statuses = table["status"].tolist()
    allowed_statuses = ("pass", "warning", "unavailable")
    if any(value not in allowed_statuses for value in statuses):
        raise ValueError(f"{operation} contains an unsupported QC status.")
    available = [
        lineage_id for lineage_id, status in zip(lineage_ids, statuses, strict=True) if status != "unavailable"
    ]
    if tuple(available) != tuple(characters.lineage_ids):
        raise ValueError(f"{operation} table lineage identity disagrees with the character artifact.")
    integer_columns = [
        "input_cells",
        "retained_cells",
        "input_characters",
        "retained_characters",
        "all_missing_cells",
        "all_uncut_cells",
        "missing_threshold_filtered_cells",
        "uncut_threshold_filtered_cells",
    ]
    numeric = {}
    for column in integer_columns:
        values = table[column].to_numpy()
        if pd.api.types.is_bool_dtype(table[column].dtype) or not pd.api.types.is_integer_dtype(table[column].dtype):
            raise TypeError(f"{operation} column {column!r} must contain nonnegative integers.")
        numeric[column] = values.astype(np.int64, copy=False)
        if bool((numeric[column] < 0).any()):
            raise ValueError(f"{operation} column {column!r} cannot be negative.")
    for column in (
        "informative_character_fraction",
        "unique_profile_fraction",
        "mean_missing_fraction",
        "mean_uncut_fraction_observed",
    ):
        if pd.api.types.is_bool_dtype(table[column].dtype) or not pd.api.types.is_numeric_dtype(table[column].dtype):
            raise TypeError(f"{operation} column {column!r} must be numeric.")
        values = table[column].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if bool(((finite < 0.0) | (finite > 1.0)).any()) or bool(np.isinf(values).any()):
            raise ValueError(f"{operation} finite {column!r} values must lie in [0, 1].")
        numeric[column] = values
    for row_index, lineage_id in enumerate(lineage_ids):
        for column in ("qc_warnings", "unavailability_reasons"):
            if not isinstance(table.at[row_index, column], list):
                raise TypeError(f"{operation} column {column!r} must retain JSON lists.")
        if statuses[row_index] != "unavailable":
            matrix, _, _ = characters.copy_lineage(lineage_id)
            if matrix.shape != (
                int(numeric["retained_cells"][row_index]),
                int(numeric["retained_characters"][row_index]),
            ):
                raise ValueError(f"{operation} retained axes disagree for lineage {lineage_id!r}.")

    row_count = len(table)
    figure_height = max(5.0, min(18.0, 2.8 + 0.24 * row_count))
    figure = Figure(figsize=(15.0, figure_height), constrained_layout=True)
    FigureCanvasAgg(figure)
    cells_axis, fractions_axis, status_axis = figure.subplots(1, 3, gridspec_kw={"width_ratios": [1.4, 2.2, 0.8]})
    positions = np.arange(row_count)
    cells_axis.barh(positions, numeric["input_cells"], color="#D9E2EC", label="Input")
    cells_axis.barh(positions, numeric["retained_cells"], color="#3E7CB1", label="Retained")
    cells_axis.set_xlabel("Cells")
    cells_axis.set_title("Cell retention")
    cells_axis.grid(axis="x", alpha=0.2)
    cells_axis.legend()
    labels = lineage_ids if row_count <= 60 else ["" for _ in lineage_ids]
    cells_axis.set_yticks(positions, labels)
    cells_axis.invert_yaxis()

    fraction_columns = [
        "informative_character_fraction",
        "unique_profile_fraction",
        "mean_missing_fraction",
        "mean_uncut_fraction_observed",
    ]
    fraction_values = np.column_stack([numeric[column] for column in fraction_columns])
    image = fractions_axis.imshow(
        np.ma.masked_invalid(fraction_values),
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
    )
    fractions_axis.set_xticks(
        np.arange(4),
        ["Informative\ncharacters", "Unique\nprofiles", "Mean\nmissing", "Mean observed\nuncut"],
    )
    fractions_axis.set_yticks(positions, labels)
    fractions_axis.set_title("Character and profile fractions")
    figure.colorbar(image, ax=fractions_axis, label="Fraction")

    status_counts = {status: statuses.count(status) for status in allowed_statuses}
    status_axis.bar(
        np.arange(3),
        [status_counts[status] for status in allowed_statuses],
        color=["#3A923A", "#D59B21", "#C44E52"],
    )
    status_axis.set_xticks(np.arange(3), ["Pass", "Warning", "Unavailable"], rotation=30, ha="right")
    status_axis.set_ylabel("Lineages")
    status_axis.set_title("QC status")
    status_axis.grid(axis="y", alpha=0.2)

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")
    details = {
        "characters_fingerprint": characters.fingerprint,
        "qc_table_sha256": expected_hash,
        "lineages_plotted": row_count,
        "lineages_available": len(available),
        "status_counts": status_counts,
        "input_cells_total": int(numeric["input_cells"].sum()),
        "retained_cells_total": int(numeric["retained_cells"].sum()),
        "lineage_labels_rendered": row_count <= 60,
    }
    return (png, details) if _return_details else png


def run_cassiopeia_lineage_qc_plot(characters, table, *, producer):
    return _standalone_cassiopeia_lineage_qc_plot(
        characters,
        table,
        producer=producer,
        _return_details=True,
    )


def cassiopeia_lineage_qc_plot_code(*, producer) -> str:
    functions = (
        _standalone_lineage_scalar,
        _standalone_lineage_table_fingerprint,
        _standalone_cassiopeia_lineage_qc_plot,
    )
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

{implementation}


def plot_cassiopeia_lineage_qc(characters, table):
    return _standalone_cassiopeia_lineage_qc_plot(characters, table, producer={producer!r})
"""


def _standalone_cassiopeia_tree_evidence(tree):
    import math
    from collections.abc import Mapping

    operation = "Cassiopeia Tree Plot"
    metadata = tree.metadata
    if not isinstance(metadata, Mapping) or metadata.get("artifact_type") != "OPENBIO_CASSIOPEIA_TREE":
        raise TypeError(f"{operation} requires an audited Cassiopeia tree artifact.")
    backend = tree.copy_tree()
    if not hasattr(backend, "get_tree_topology") or not callable(backend.get_tree_topology):
        raise TypeError(f"{operation} tree backend does not expose its stored topology.")
    graph = backend.get_tree_topology()
    if graph is None or not hasattr(graph, "nodes") or not hasattr(graph, "edges"):
        raise ValueError(f"{operation} stored topology is unavailable.")
    if hasattr(graph, "is_directed") and not bool(graph.is_directed()):
        raise ValueError(f"{operation} requires a directed topology.")
    nodes = list(graph.nodes)
    if (
        not nodes
        or len(nodes) > 2_000
        or len(nodes) != len(set(nodes))
        or any(not isinstance(node, str) or not node or node != node.strip() for node in nodes)
    ):
        raise ValueError(f"{operation} requires 1 through 2,000 uniquely named topology nodes.")
    parents = {node: [] for node in nodes}
    children = {node: [] for node in nodes}
    lengths = {}
    for parent, child, attributes in graph.edges(data=True):
        if parent not in children or child not in parents or parent == child:
            raise ValueError(f"{operation} topology contains an invalid edge.")
        length = float(attributes.get("length", 1.0))
        if not math.isfinite(length) or length < 0.0:
            raise ValueError(f"{operation} branch lengths must be finite and nonnegative.")
        parents[child].append(parent)
        children[parent].append(child)
        lengths[(parent, child)] = length
    roots = [node for node in nodes if not parents[node]]
    if len(roots) != 1 or any(len(parents[node]) != 1 for node in nodes if node != roots[0]):
        raise ValueError(f"{operation} topology must be one rooted arborescence.")
    if len(lengths) != len(nodes) - 1:
        raise ValueError(f"{operation} topology edge count is inconsistent.")
    root = roots[0]
    if str(metadata.get("root")) != root or str(getattr(backend, "root", "")) != root:
        raise ValueError(f"{operation} root identity disagrees with stored provenance.")
    visiting = set()
    visited = set()
    depth = {}

    def visit(node, value):
        if node in visiting:
            raise ValueError(f"{operation} topology contains a cycle.")
        if node in visited:
            return
        visiting.add(node)
        depth[node] = value
        children[node].sort()
        for child in children[node]:
            visit(child, value + 1)
        visiting.remove(node)
        visited.add(node)

    visit(root, 0)
    if visited != set(nodes):
        raise ValueError(f"{operation} topology is disconnected.")
    leaves = [node for node in nodes if not children[node]]
    if set(leaves) != set(metadata.get("leaf_axis", [])):
        raise ValueError(f"{operation} leaf identity disagrees with stored provenance.")
    leaf_order = []
    y = {}

    def position(node):
        if not children[node]:
            y[node] = float(len(leaf_order))
            leaf_order.append(node)
        else:
            for child in children[node]:
                position(child)
            y[node] = sum(y[child] for child in children[node]) / len(children[node])

    position(root)
    return {
        "metadata": metadata,
        "root": root,
        "nodes": nodes,
        "children": children,
        "parents": parents,
        "lengths": lengths,
        "depth": depth,
        "y": y,
        "leaves": leaf_order,
    }


def _standalone_cassiopeia_tree_plot(
    tree,
    adata,
    *,
    annotation_key,
    _return_details=False,
):
    import io

    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    operation = "Cassiopeia Tree Plot"
    evidence = _standalone_cassiopeia_tree_evidence(tree)
    if not hasattr(adata, "obs") or not hasattr(adata, "obs_names"):
        raise TypeError(f"{operation} requires an AnnData-like annotation input.")
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires in-memory AnnData.")
    if not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires unique AnnData observation identities.")
    obs_names = list(adata.obs_names)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in obs_names):
        raise ValueError(f"{operation} AnnData observation identities must be exact nonblank strings.")
    if not isinstance(annotation_key, str) or not annotation_key or annotation_key != annotation_key.strip():
        raise ValueError(f"{operation} annotation_key must be exact nonblank text.")
    if annotation_key not in adata.obs:
        raise ValueError(f"{operation} AnnData is missing annotation column {annotation_key!r}.")
    annotation = adata.obs[annotation_key]
    if not isinstance(annotation.dtype, pd.CategoricalDtype):
        raise ValueError(f"{operation} annotation must use pandas categorical dtype.")
    missing_leaves = sorted(set(evidence["leaves"]) - set(obs_names))
    if missing_leaves:
        raise ValueError(f"{operation} AnnData is missing {len(missing_leaves)} exact tree leaf identities.")
    values = []
    for leaf in evidence["leaves"]:
        value = annotation.loc[leaf]
        if bool(pd.isna(value)):
            raise ValueError(f"{operation} tree leaf {leaf!r} has a missing annotation.")
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(f"{operation} leaf annotations must be exact nonblank strings.")
        values.append(value)
    observed = set(values)
    categories = [str(value) for value in annotation.cat.categories if str(value) in observed]
    if set(categories) != observed:
        raise ValueError(f"{operation} annotation categories do not close over tree leaves.")
    palette = [
        "#4C78A8",
        "#F58518",
        "#54A24B",
        "#E45756",
        "#72B7B2",
        "#B279A2",
        "#FF9DA6",
        "#9D755D",
        "#BAB0AC",
        "#1F77B4",
        "#FF7F0E",
        "#2CA02C",
        "#D62728",
        "#9467BD",
        "#8C564B",
        "#E377C2",
        "#7F7F7F",
        "#BCBD22",
        "#17BECF",
        "#AEC7E8",
    ]
    if len(categories) > len(palette):
        raise ValueError(f"{operation} supports at most {len(palette)} annotation categories in one legend.")
    colors = dict(zip(categories, palette, strict=False))
    leaf_count = len(evidence["leaves"])
    figure = Figure(figsize=(12.0, max(5.0, min(20.0, 2.8 + 0.24 * leaf_count))), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    for parent in evidence["nodes"]:
        children = evidence["children"][parent]
        if not children:
            continue
        child_y = [evidence["y"][child] for child in children]
        axis.plot(
            [evidence["depth"][parent], evidence["depth"][parent]],
            [min(child_y), max(child_y)],
            color="#7A7A7A",
            linewidth=0.8,
        )
        for child in children:
            axis.plot(
                [evidence["depth"][parent], evidence["depth"][child]],
                [evidence["y"][child], evidence["y"][child]],
                color="#7A7A7A",
                linewidth=0.8,
            )
    leaf_x = [evidence["depth"][leaf] for leaf in evidence["leaves"]]
    leaf_y = [evidence["y"][leaf] for leaf in evidence["leaves"]]
    axis.scatter(
        leaf_x, leaf_y, c=[colors[value] for value in values], s=30, edgecolors="white", linewidths=0.4, zorder=3
    )
    show_labels = leaf_count <= 80
    if show_labels:
        for leaf in evidence["leaves"]:
            axis.annotate(
                leaf,
                (evidence["depth"][leaf], evidence["y"][leaf]),
                xytext=(5, 0),
                textcoords="offset points",
                va="center",
                fontsize=7,
            )
    axis.set_xlabel("Topological depth")
    axis.set_ylabel("Leaves")
    axis.set_yticks([])
    axis.set_title(f"Cassiopeia lineage {tree.lineage_id}: {annotation_key}")
    axis.grid(axis="x", alpha=0.15)
    axis.legend(
        handles=[Line2D([0], [0], marker="o", linestyle="", color=colors[value], label=value) for value in categories],
        title=annotation_key,
        loc="best",
    )
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")
    details = {
        "lineage_id": tree.lineage_id,
        "tree_fingerprint": tree.fingerprint,
        "topology_sha256": tree.topology_fingerprint,
        "root": evidence["root"],
        "node_count": len(evidence["nodes"]),
        "edge_count": len(evidence["nodes"]) - 1,
        "leaf_count": leaf_count,
        "maximum_depth": max(evidence["depth"].values()),
        "annotation_key": annotation_key,
        "annotation_categories": categories,
        "ann_data_cells_not_in_tree": len(set(obs_names) - set(evidence["leaves"])),
        "leaf_labels_rendered": show_labels,
    }
    return (png, details) if _return_details else png


def run_cassiopeia_tree_plot(tree, adata, *, annotation_key):
    return _standalone_cassiopeia_tree_plot(
        tree,
        adata,
        annotation_key=annotation_key,
        _return_details=True,
    )


def cassiopeia_tree_plot_code(*, annotation_key: str) -> str:
    functions = (_standalone_cassiopeia_tree_evidence, _standalone_cassiopeia_tree_plot)
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

{implementation}


def plot_cassiopeia_tree(tree, adata):
    return _standalone_cassiopeia_tree_plot(tree, adata, annotation_key={annotation_key!r})
"""


def _standalone_cassiopeia_expansion_plot(
    tree,
    table,
    *,
    producer,
    _return_details=False,
):
    import io
    import math
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Cassiopeia Expansion Plot"
    evidence = _standalone_cassiopeia_tree_evidence(tree)
    expected_producer = {"operation", "parameters", "input_cells", "input_genes", "random_seed"}
    if not isinstance(producer, Mapping) or set(producer) != expected_producer:
        raise ValueError(f"{operation} requires exact expansion-table producer provenance.")
    if producer["operation"] != "cassiopeia_clade_expansion_test":
        raise ValueError(f"{operation} requires the canonical Cassiopeia clade-expansion table.")
    parameters = producer["parameters"]
    if not isinstance(parameters, Mapping):
        raise ValueError(f"{operation} producer parameters are invalid.")
    identity = {
        "lineage_id": tree.lineage_id,
        "tree_fingerprint": tree.fingerprint,
        "topology_sha256": tree.topology_fingerprint,
    }
    if any(parameters.get(key) != value for key, value in identity.items()):
        raise ValueError(f"{operation} table is not bound to this exact tree topology.")
    expected_hash = parameters.get("table_sha256")
    if not isinstance(expected_hash, str):
        raise ValueError(f"{operation} table fingerprint provenance is missing.")
    columns = [
        "node_id",
        "parent_id",
        "depth",
        "leaf_count",
        "child_count",
        "eligible",
        "exclusion_reason",
        "p_value",
        "q_value_bh",
        "significant_fdr",
    ]
    if not isinstance(table, pd.DataFrame) or list(table.columns) != columns:
        raise ValueError(f"{operation} table does not match the canonical clade-expansion schema.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError(f"{operation} table must use the canonical zero-based row index.")
    if _standalone_lineage_table_fingerprint(table) != expected_hash:
        raise ValueError(f"{operation} table content was changed after production.")
    if len(table) != len(evidence["nodes"]):
        raise ValueError(f"{operation} table must retain every tested topology node.")
    node_ids = table["node_id"].tolist()
    if (
        any(not isinstance(value, str) or not value or value != value.strip() for value in node_ids)
        or len(node_ids) != len(set(node_ids))
        or set(node_ids) != set(evidence["nodes"])
    ):
        raise ValueError(f"{operation} node identity does not match the tree topology.")
    if any(
        pd.api.types.is_bool_dtype(table[column].dtype) or not pd.api.types.is_integer_dtype(table[column].dtype)
        for column in ("depth", "leaf_count", "child_count")
    ):
        raise TypeError(f"{operation} topology counts must be integer columns.")
    if not pd.api.types.is_bool_dtype(table["eligible"].dtype) or not pd.api.types.is_bool_dtype(
        table["significant_fdr"].dtype
    ):
        raise TypeError(f"{operation} eligibility and FDR calls must be Boolean columns.")
    for column in ("p_value", "q_value_bh"):
        if pd.api.types.is_bool_dtype(table[column].dtype) or not pd.api.types.is_numeric_dtype(table[column].dtype):
            raise TypeError(f"{operation} probability columns must be numeric.")
        values = table[column].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if bool(np.isinf(values).any()) or bool(((finite < 0.0) | (finite > 1.0)).any()):
            raise ValueError(f"{operation} finite probabilities must lie in [0, 1].")

    leaf_counts = {}
    for node in sorted(evidence["nodes"], key=lambda value: evidence["depth"][value], reverse=True):
        children = evidence["children"][node]
        leaf_counts[node] = 1 if not children else sum(leaf_counts[child] for child in children)
    records = table.set_index("node_id", drop=False)
    fdr_threshold = float(parameters.get("fdr_threshold", -1.0))
    if not math.isfinite(fdr_threshold) or not 0.0 <= fdr_threshold <= 1.0:
        raise ValueError(f"{operation} producer FDR threshold is invalid.")
    for node in evidence["nodes"]:
        row = records.loc[node]
        parent = evidence["parents"][node][0] if evidence["parents"][node] else None
        observed_parent = row["parent_id"]
        observed_parent = None if pd.isna(observed_parent) else observed_parent
        if (
            observed_parent != parent
            or int(row["depth"]) != evidence["depth"][node]
            or int(row["leaf_count"]) != leaf_counts[node]
            or int(row["child_count"]) != len(evidence["children"][node])
        ):
            raise ValueError(f"{operation} topology columns disagree at node {node!r}.")
        p_value = row["p_value"]
        q_value = row["q_value_bh"]
        eligible = bool(row["eligible"])
        if eligible != (not pd.isna(p_value) and not pd.isna(q_value)):
            raise ValueError(f"{operation} eligible rows must retain both raw and adjusted probabilities.")
        expected_call = eligible and float(q_value) <= fdr_threshold
        if bool(row["significant_fdr"]) != expected_call:
            raise ValueError(f"{operation} FDR calls disagree with the stored threshold.")

    figure = Figure(
        figsize=(14.0, max(5.0, min(18.0, 2.8 + 0.18 * len(evidence["leaves"])))),
        constrained_layout=True,
    )
    FigureCanvasAgg(figure)
    tree_axis, evidence_axis = figure.subplots(1, 2, gridspec_kw={"width_ratios": [1.3, 1.0]})
    for parent in evidence["nodes"]:
        children = evidence["children"][parent]
        if not children:
            continue
        child_y = [evidence["y"][child] for child in children]
        tree_axis.plot(
            [evidence["depth"][parent], evidence["depth"][parent]],
            [min(child_y), max(child_y)],
            color="#8A8A8A",
            linewidth=0.8,
        )
        for child in children:
            tree_axis.plot(
                [evidence["depth"][parent], evidence["depth"][child]],
                [evidence["y"][child], evidence["y"][child]],
                color="#8A8A8A",
                linewidth=0.8,
            )
    row_by_node = {row.node_id: row for row in table.itertuples(index=False)}
    node_colors = []
    node_sizes = []
    for node in evidence["nodes"]:
        row = row_by_node[node]
        node_colors.append("#D64A4A" if row.significant_fdr else ("#3E7CB1" if row.eligible else "#B8B8B8"))
        node_sizes.append(18.0 + 80.0 * leaf_counts[node] / len(evidence["leaves"]))
    tree_axis.scatter(
        [evidence["depth"][node] for node in evidence["nodes"]],
        [evidence["y"][node] for node in evidence["nodes"]],
        c=node_colors,
        s=node_sizes,
        edgecolors="white",
        linewidths=0.4,
        zorder=3,
    )
    tree_axis.set_xlabel("Topological depth")
    tree_axis.set_ylabel("Leaves")
    tree_axis.set_yticks([])
    tree_axis.set_title("Expansion calls on stored topology")
    tree_axis.grid(axis="x", alpha=0.15)

    eligible_table = table.loc[table["eligible"]].copy()
    if not eligible_table.empty:
        q_values = eligible_table["q_value_bh"].to_numpy(dtype=float)
        y_values = -np.log10(np.maximum(q_values, 1e-300))
        evidence_axis.scatter(
            eligible_table["depth"].to_numpy(dtype=int),
            y_values,
            s=20.0 + 120.0 * eligible_table["leaf_count"].to_numpy(dtype=float) / len(evidence["leaves"]),
            c=["#D64A4A" if value else "#3E7CB1" for value in eligible_table["significant_fdr"]],
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
        )
    evidence_axis.axhline(
        -math.log10(max(fdr_threshold, 1e-300)),
        color="#555555",
        linestyle="--",
        linewidth=0.9,
    )
    evidence_axis.set_xlabel("Topological depth")
    evidence_axis.set_ylabel("−log10(BH adjusted p-value)")
    evidence_axis.set_title("Eligible clade evidence")
    evidence_axis.grid(alpha=0.2)

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")
    details = {
        "lineage_id": tree.lineage_id,
        "tree_fingerprint": tree.fingerprint,
        "topology_sha256": tree.topology_fingerprint,
        "table_sha256": expected_hash,
        "topology_nodes_plotted": len(evidence["nodes"]),
        "leaves_plotted": len(evidence["leaves"]),
        "eligible_clades": int(table["eligible"].sum()),
        "significant_clades": int(table["significant_fdr"].sum()),
        "fdr_threshold": fdr_threshold,
    }
    return (png, details) if _return_details else png


def run_cassiopeia_expansion_plot(tree, table, *, producer):
    return _standalone_cassiopeia_expansion_plot(
        tree,
        table,
        producer=producer,
        _return_details=True,
    )


def cassiopeia_expansion_plot_code(*, producer) -> str:
    functions = (
        _standalone_lineage_scalar,
        _standalone_lineage_table_fingerprint,
        _standalone_cassiopeia_tree_evidence,
        _standalone_cassiopeia_expansion_plot,
    )
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

{implementation}


def plot_cassiopeia_expansion(tree, table):
    return _standalone_cassiopeia_expansion_plot(tree, table, producer={producer!r})
"""


def _standalone_cassiopeia_plasticity_plot(
    adata,
    tree,
    table,
    *,
    producer,
    output_key,
    _return_details=False,
):
    import io
    import math
    from collections.abc import Mapping, Sequence

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.figure import Figure

    operation = "Cassiopeia Plasticity Plot"
    evidence = _standalone_cassiopeia_tree_evidence(tree)
    expected_producer = {"operation", "parameters", "input_cells", "input_genes", "random_seed"}
    if not isinstance(producer, Mapping) or set(producer) != expected_producer:
        raise ValueError(f"{operation} requires exact plasticity-table producer provenance.")
    if producer["operation"] != "cassiopeia_effective_plasticity":
        raise ValueError(f"{operation} requires the canonical Cassiopeia EffectivePlasticity table.")
    parameters = producer["parameters"]
    if not isinstance(parameters, Mapping):
        raise ValueError(f"{operation} producer parameters are invalid.")
    if not isinstance(output_key, str) or not output_key or output_key != output_key.strip():
        raise ValueError(f"{operation} output_key must be exact nonblank text.")
    identity = {
        "lineage_id": tree.lineage_id,
        "tree_fingerprint": tree.fingerprint,
        "topology_sha256": tree.topology_fingerprint,
        "output_key": output_key,
    }
    if any(parameters.get(key) != value for key, value in identity.items()):
        raise ValueError(f"{operation} table is not bound to this exact tree and output key.")
    expected_hash = parameters.get("table_sha256")
    if not isinstance(expected_hash, str):
        raise ValueError(f"{operation} table fingerprint provenance is missing.")
    if not hasattr(adata, "obs") or not hasattr(adata, "obs_names") or not hasattr(adata, "uns"):
        raise TypeError(f"{operation} requires an AnnData-like plasticity result.")
    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires in-memory AnnData.")
    if not bool(adata.obs_names.is_unique) or not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique AnnData axes.")
    obs_names = list(adata.obs_names)
    if any(not isinstance(value, str) or not value or value != value.strip() for value in obs_names):
        raise ValueError(f"{operation} observation identities must be exact nonblank strings.")
    if producer["input_cells"] != int(adata.n_obs) or producer["input_genes"] != int(adata.n_vars):
        raise ValueError(f"{operation} AnnData axes do not match table producer provenance.")
    registry = adata.uns.get("openbio_cassiopeia_plasticity")
    provenance = registry.get(output_key) if isinstance(registry, Mapping) else None
    expected_provenance = {
        "producer",
        "node_id",
        "output_key",
        "status_key",
        "tree_fingerprint",
        "topology_sha256",
        "annotation_key",
        "annotation_status",
        "analysis_mode",
        "annotation_provenance_sha256",
        "minimum_state_fraction",
        "retained_states",
        "global_small_parsimony",
    }
    if not isinstance(provenance, Mapping) or set(provenance) != expected_provenance:
        raise ValueError(f"{operation} requires exact stored EffectivePlasticity provenance.")
    if (
        provenance["producer"] != "openbio-singlecell"
        or provenance["node_id"] != "OpenBioSingleCellCassiopeiaPlasticity"
        or provenance["output_key"] != output_key
        or provenance["tree_fingerprint"] != tree.fingerprint
        or provenance["topology_sha256"] != tree.topology_fingerprint
    ):
        raise ValueError(f"{operation} AnnData provenance is not bound to this tree and score column.")
    annotation_key = provenance["annotation_key"]
    status_key = provenance["status_key"]
    if annotation_key not in adata.obs or output_key not in adata.obs or status_key not in adata.obs:
        raise ValueError(f"{operation} AnnData is missing stored annotation, score, or status columns.")
    annotation = adata.obs[annotation_key]
    if not isinstance(annotation.dtype, pd.CategoricalDtype):
        raise ValueError(f"{operation} annotation must retain pandas categorical dtype.")
    if pd.api.types.is_bool_dtype(adata.obs[output_key].dtype) or not pd.api.types.is_numeric_dtype(
        adata.obs[output_key].dtype
    ):
        raise TypeError(f"{operation} score column must be numeric.")
    columns = [
        "cell_id",
        "lineage_id",
        "annotation_state",
        "state_count_in_tree",
        "state_fraction_in_tree",
        "status",
        "path_subtree_count",
        "sc_effective_plasticity",
    ]
    if not isinstance(table, pd.DataFrame) or list(table.columns) != columns:
        raise ValueError(f"{operation} table does not match the canonical EffectivePlasticity schema.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError(f"{operation} table must use the canonical zero-based row index.")
    if len(table) != int(adata.n_obs) or table["cell_id"].tolist() != obs_names:
        raise ValueError(f"{operation} table cell axis must exactly match AnnData order.")
    if _standalone_lineage_table_fingerprint(table) != expected_hash:
        raise ValueError(f"{operation} table content was changed after production.")
    allowed_statuses = ("included", "missing_annotation", "state_below_minimum_frequency", "not_in_tree")
    table_status = table["status"].tolist()
    if any(value not in allowed_statuses for value in table_status):
        raise ValueError(f"{operation} table contains an unsupported cell status.")
    if [str(value) for value in adata.obs[status_key].tolist()] != table_status:
        raise ValueError(f"{operation} AnnData status column disagrees with the canonical table.")
    tree_leaves = set(evidence["leaves"])
    retained_states = provenance["retained_states"]
    if hasattr(retained_states, "tolist"):
        retained_states = retained_states.tolist()
    if (
        not isinstance(retained_states, Sequence)
        or isinstance(retained_states, (str, bytes, bytearray))
        or any(not isinstance(value, str) or not value or value != value.strip() for value in retained_states)
        or list(retained_states) != sorted(set(retained_states))
    ):
        raise ValueError(f"{operation} retained annotation-state identity is invalid.")
    retained_states = list(retained_states)
    included_scores = []
    included_states = []
    score_by_leaf = {}
    for row_index, cell_id in enumerate(obs_names):
        row = table.iloc[row_index]
        raw_annotation = annotation.iloc[row_index]
        annotation_value = None if pd.isna(raw_annotation) else str(raw_annotation)
        table_annotation = None if pd.isna(row["annotation_state"]) else row["annotation_state"]
        if annotation_value != table_annotation:
            raise ValueError(f"{operation} annotation values disagree at cell {cell_id!r}.")
        expected_lineage = None if cell_id not in tree_leaves else tree.lineage_id
        table_lineage = None if pd.isna(row["lineage_id"]) else row["lineage_id"]
        if table_lineage != expected_lineage:
            raise ValueError(f"{operation} lineage membership disagrees at cell {cell_id!r}.")
        adata_score = adata.obs[output_key].iloc[row_index]
        table_score = row["sc_effective_plasticity"]
        if pd.isna(adata_score) != pd.isna(table_score):
            raise ValueError(f"{operation} score missingness disagrees at cell {cell_id!r}.")
        if not pd.isna(adata_score) and not math.isclose(
            float(adata_score), float(table_score), rel_tol=0.0, abs_tol=1e-12
        ):
            raise ValueError(f"{operation} stored scores disagree at cell {cell_id!r}.")
        if row["status"] == "included":
            score = float(table_score)
            if (
                cell_id not in tree_leaves
                or annotation_value not in retained_states
                or not math.isfinite(score)
                or not 0 <= score <= 1
            ):
                raise ValueError(f"{operation} included score evidence is invalid at cell {cell_id!r}.")
            included_scores.append(score)
            included_states.append(annotation_value)
            score_by_leaf[cell_id] = score
        elif not pd.isna(table_score):
            raise ValueError(f"{operation} excluded cells must retain null scores.")
    if not included_scores:
        raise ValueError(f"{operation} requires at least one included score.")
    observed_states = set(included_states)
    categories = [str(value) for value in annotation.cat.categories if str(value) in observed_states]
    if set(categories) != observed_states:
        raise ValueError(f"{operation} annotation categories do not close over included scores.")

    figure = Figure(
        figsize=(14.0, max(5.0, min(18.0, 2.8 + 0.18 * len(evidence["leaves"])))),
        constrained_layout=True,
    )
    FigureCanvasAgg(figure)
    tree_axis, distribution_axis = figure.subplots(1, 2, gridspec_kw={"width_ratios": [1.3, 1.0]})
    for parent in evidence["nodes"]:
        children = evidence["children"][parent]
        if not children:
            continue
        child_y = [evidence["y"][child] for child in children]
        tree_axis.plot(
            [evidence["depth"][parent], evidence["depth"][parent]],
            [min(child_y), max(child_y)],
            color="#8A8A8A",
            linewidth=0.8,
        )
        for child in children:
            tree_axis.plot(
                [evidence["depth"][parent], evidence["depth"][child]],
                [evidence["y"][child], evidence["y"][child]],
                color="#8A8A8A",
                linewidth=0.8,
            )
    normalization = Normalize(vmin=0.0, vmax=max(included_scores) if max(included_scores) > 0 else 1.0)
    color_map = "viridis"
    available_leaves = [leaf for leaf in evidence["leaves"] if leaf in score_by_leaf]
    missing_score_leaves = [leaf for leaf in evidence["leaves"] if leaf not in score_by_leaf]
    tree_axis.scatter(
        [evidence["depth"][leaf] for leaf in available_leaves],
        [evidence["y"][leaf] for leaf in available_leaves],
        c=[score_by_leaf[leaf] for leaf in available_leaves],
        cmap=color_map,
        norm=normalization,
        s=32,
        edgecolors="white",
        linewidths=0.4,
        zorder=3,
    )
    if missing_score_leaves:
        tree_axis.scatter(
            [evidence["depth"][leaf] for leaf in missing_score_leaves],
            [evidence["y"][leaf] for leaf in missing_score_leaves],
            color="#C7C7C7",
            s=25,
            edgecolors="white",
            linewidths=0.4,
            zorder=3,
        )
    figure.colorbar(ScalarMappable(norm=normalization, cmap=color_map), ax=tree_axis, label="EffectivePlasticity")
    tree_axis.set_xlabel("Topological depth")
    tree_axis.set_ylabel("Leaves")
    tree_axis.set_yticks([])
    tree_axis.set_title("Stored single-cell scores on topology")
    tree_axis.grid(axis="x", alpha=0.15)

    grouped = [
        [score for score, state in zip(included_scores, included_states, strict=True) if state == category]
        for category in categories
    ]
    distribution_axis.boxplot(grouped, tick_labels=categories, showfliers=False)
    for position, values in enumerate(grouped, start=1):
        offsets = np.linspace(-0.14, 0.14, len(values)) if len(values) > 1 else np.asarray([0.0])
        distribution_axis.scatter(position + offsets, sorted(values), color="#3E7CB1", s=24, alpha=0.8, zorder=3)
    distribution_axis.set_xlabel(annotation_key)
    distribution_axis.set_ylabel("EffectivePlasticity")
    distribution_axis.set_ylim(bottom=0.0)
    distribution_axis.set_title("Scores by retained annotation state")
    distribution_axis.tick_params(axis="x", rotation=30)
    distribution_axis.grid(axis="y", alpha=0.2)

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} did not produce a valid PNG payload.")
    status_counts = {status: table_status.count(status) for status in allowed_statuses}
    details = {
        "lineage_id": tree.lineage_id,
        "tree_fingerprint": tree.fingerprint,
        "topology_sha256": tree.topology_fingerprint,
        "table_sha256": expected_hash,
        "output_key": output_key,
        "status_key": status_key,
        "annotation_key": annotation_key,
        "annotation_categories": categories,
        "tree_leaves": len(evidence["leaves"]),
        "included_cells": len(included_scores),
        "tree_leaves_without_scores": len(missing_score_leaves),
        "status_counts": status_counts,
        "minimum_score": float(min(included_scores)),
        "maximum_score": float(max(included_scores)),
        "mean_score": float(sum(included_scores) / len(included_scores)),
    }
    return (png, details) if _return_details else png


def run_cassiopeia_plasticity_plot(adata, tree, table, *, producer, output_key):
    return _standalone_cassiopeia_plasticity_plot(
        adata,
        tree,
        table,
        producer=producer,
        output_key=output_key,
        _return_details=True,
    )


def cassiopeia_plasticity_plot_code(*, producer, output_key: str) -> str:
    functions = (
        _standalone_lineage_scalar,
        _standalone_lineage_table_fingerprint,
        _standalone_cassiopeia_tree_evidence,
        _standalone_cassiopeia_plasticity_plot,
    )
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

{implementation}


def plot_cassiopeia_plasticity(adata, tree, table):
    return _standalone_cassiopeia_plasticity_plot(
        adata,
        tree,
        table,
        producer={producer!r},
        output_key={output_key!r},
    )
"""


__all__ = [
    "cassiopeia_expansion_plot_code",
    "cassiopeia_lineage_qc_plot_code",
    "cassiopeia_plasticity_plot_code",
    "cassiopeia_tree_plot_code",
    "lineage_table_fingerprint",
    "run_cassiopeia_expansion_plot",
    "run_cassiopeia_lineage_qc_plot",
    "run_cassiopeia_plasticity_plot",
    "run_cassiopeia_tree_plot",
]
