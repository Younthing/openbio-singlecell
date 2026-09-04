from __future__ import annotations

import inspect
import textwrap
from typing import Any


def _standalone_marker_evidence_plot(
    adata,
    table,
    universe,
    *,
    producer,
    source_kind="layer",
    layer_name="log1p_norm",
    top_genes_per_group=5,
    _return_details=False,
):
    """Render an upstream-ranked marker panel without recomputing marker evidence."""
    import hashlib
    import io
    import json
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from scipy import sparse

    operation = "Marker Evidence Plot"
    marker_columns = [
        "group",
        "gene",
        "rank",
        "score",
        "log2_fold_change_approx",
        "p_value",
        "p_adjusted",
        "fraction_in_group",
        "fraction_reference",
    ]
    universe_columns = ["gene", "universe_rank"]
    numeric_columns = set(marker_columns[2:])

    def canonical_sha256(payload):
        encoded = json.dumps(
            payload,
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def frame_fingerprint(frame, *, artifact, numeric):
        encoder = json.JSONEncoder(
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        digest = hashlib.sha256()

        def update(value):
            for chunk in encoder.iterencode(value):
                digest.update(chunk.encode("utf-8"))

        digest.update(b'{"columns":')
        update(list(frame.columns))
        digest.update(b',"rows":[')
        first = True
        for row in frame.itertuples(index=False, name=None):
            encoded_row = [
                float(value).hex() if column in numeric else value
                for column, value in zip(frame.columns, row, strict=True)
            ]
            if not first:
                digest.update(b",")
            update(encoded_row)
            first = False
        digest.update(b'],"schema":')
        update(f"openbio-singlecell/{artifact}/v2")
        digest.update(b"}")
        return digest.hexdigest()

    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires in-memory AnnData.")
    if int(getattr(adata, "n_obs", 0)) < 1 or int(getattr(adata, "n_vars", 0)) < 1:
        raise ValueError(f"{operation} requires non-empty observation and feature axes.")
    if int(adata.n_obs) > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 cells in one static PNG.")
    if not bool(adata.obs_names.is_unique) or not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique observation and feature identifiers.")
    if not isinstance(table, pd.DataFrame) or list(table.columns) != marker_columns:
        raise ValueError(f"{operation} requires the exact canonical Cluster marker evidence columns.")
    if not isinstance(universe, pd.DataFrame) or list(universe.columns) != universe_columns:
        raise ValueError(f"{operation} requires the exact canonical tested-gene universe columns.")
    if table.empty:
        raise ValueError(f"{operation} cannot render an empty marker table.")
    if bool(table.duplicated(["group", "gene"]).any()):
        raise ValueError(f"{operation} marker evidence contains duplicate group/gene rows.")
    for column in ("group", "gene"):
        values = table[column].tolist()
        if any(not isinstance(value, str) or not value or value != value.strip() for value in values):
            raise ValueError(f"{operation} marker {column} values must be canonical strings.")
    for column in marker_columns[2:]:
        series = table[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(series.dtype):
            raise TypeError(f"{operation} marker {column!r} must be numeric.")
        if not bool(np.isfinite(series.to_numpy(dtype=float)).all()):
            raise ValueError(f"{operation} marker {column!r} must be finite.")
    ranks = table["rank"].to_numpy(dtype=float)
    if not bool(np.equal(ranks, np.floor(ranks)).all()) or bool((ranks < 1).any()):
        raise ValueError(f"{operation} marker ranks must be positive integers.")
    previous_group = None
    seen_groups = set()
    previous_rank = 0
    for group, rank in zip(table["group"].tolist(), ranks, strict=True):
        if group != previous_group:
            if group in seen_groups:
                raise ValueError(f"{operation} marker group blocks must be contiguous.")
            seen_groups.add(group)
            previous_group = group
            previous_rank = 0
        if int(rank) <= previous_rank:
            raise ValueError(f"{operation} marker ranks must increase within each group.")
        previous_rank = int(rank)
    for column in ("p_value", "p_adjusted", "fraction_in_group", "fraction_reference"):
        values = table[column].to_numpy(dtype=float)
        if bool(((values < 0.0) | (values > 1.0)).any()):
            raise ValueError(f"{operation} marker {column!r} must lie in [0, 1].")

    universe_genes = universe["gene"].tolist()
    if (
        not universe_genes
        or any(not isinstance(value, str) or not value or value != value.strip() for value in universe_genes)
        or len(set(universe_genes)) != len(universe_genes)
    ):
        raise ValueError(f"{operation} tested-gene universe requires unique canonical strings.")
    universe_ranks = universe["universe_rank"]
    if pd.api.types.is_bool_dtype(universe_ranks.dtype) or not pd.api.types.is_numeric_dtype(
        universe_ranks.dtype
    ):
        raise TypeError(f"{operation} universe ranks must be numeric.")
    if not bool(
        np.array_equal(
            universe_ranks.to_numpy(dtype=float),
            np.arange(1, len(universe) + 1, dtype=float),
        )
    ):
        raise ValueError(f"{operation} universe ranks must be consecutive and one-based.")
    if not set(table["gene"]).issubset(set(universe_genes)):
        raise ValueError(f"{operation} marker genes must belong to the tested-gene universe.")

    expected_producer_keys = {
        "table_operation",
        "universe_operation",
        "table_parameters",
        "universe_parameters",
        "input_cells",
        "input_genes",
    }
    if not isinstance(producer, Mapping) or set(producer) != expected_producer_keys:
        raise ValueError(f"{operation} requires exact marker producer provenance.")
    if producer["table_operation"] not in {"marker_genes", "filter_marker_genes"}:
        raise ValueError(f"{operation} requires a marker_genes or filter_marker_genes table operation.")
    if producer["universe_operation"] != "marker_genes":
        raise ValueError(f"{operation} requires a tested universe from marker_genes.")
    table_parameters = producer["table_parameters"]
    universe_parameters = producer["universe_parameters"]
    if not isinstance(table_parameters, Mapping) or not isinstance(universe_parameters, Mapping):
        raise ValueError(f"{operation} producer parameters must be mappings.")
    required_parameters = {
        "marker_evidence_schema_version",
        "producer_node_id",
        "artifact_role",
        "analysis_fingerprint",
        "ranking_fingerprint",
        "universe_fingerprint",
        "content_fingerprint",
        "marker_groupby",
        "marker_method",
        "marker_source",
        "group_labels",
        "group_sizes",
        "source_gene_count",
        "requested_n_genes",
        "actual_n_genes_per_group",
        "tie_correct",
    }
    if not required_parameters.issubset(table_parameters) or not required_parameters.issubset(
        universe_parameters
    ):
        raise ValueError(f"{operation} marker provenance is missing required schema fields.")
    expected_table_node = {
        "marker_genes": "OpenBioSingleCellMarkerGenes",
        "filter_marker_genes": "OpenBioSingleCellFilterMarkerGenes",
    }[producer["table_operation"]]
    if (
        table_parameters["marker_evidence_schema_version"] != 2
        or universe_parameters["marker_evidence_schema_version"] != 2
        or table_parameters["producer_node_id"] != expected_table_node
        or universe_parameters["producer_node_id"] != "OpenBioSingleCellMarkerGenes"
        or table_parameters["artifact_role"] != "marker_table"
        or universe_parameters["artifact_role"] != "tested_gene_universe"
    ):
        raise ValueError(f"{operation} marker producer/schema contract is invalid.")
    shared = (
        "analysis_fingerprint",
        "ranking_fingerprint",
        "universe_fingerprint",
        "marker_groupby",
        "marker_method",
        "marker_source",
        "marker_layer_name",
        "group_labels",
        "group_sizes",
        "source_gene_count",
        "requested_n_genes",
        "actual_n_genes_per_group",
        "tie_correct",
    )
    if any(table_parameters.get(key) != universe_parameters.get(key) for key in shared):
        raise ValueError(f"{operation} marker table and universe provenance disagree.")
    if producer["input_cells"] != int(adata.n_obs) or producer["input_genes"] != len(universe_genes):
        raise ValueError(f"{operation} input axis sizes do not match marker provenance.")
    table_fingerprint = frame_fingerprint(table, artifact="marker-table", numeric=numeric_columns)
    universe_content = frame_fingerprint(
        universe,
        artifact="tested-gene-universe",
        numeric={"universe_rank"},
    )
    if table_parameters["content_fingerprint"] != table_fingerprint:
        raise ValueError(f"{operation} marker table current-content fingerprint does not match provenance.")
    if universe_parameters["content_fingerprint"] != universe_content:
        raise ValueError(f"{operation} universe current-content fingerprint does not match provenance.")
    universe_fingerprint = canonical_sha256(
        {
            "schema": "openbio-singlecell/tested-gene-universe-identity/v2",
            "ordered_genes": universe_genes,
        }
    )
    if (
        table_parameters["universe_fingerprint"] != universe_fingerprint
        or universe_parameters["universe_fingerprint"] != universe_fingerprint
    ):
        raise ValueError(f"{operation} universe identity does not match provenance.")

    if source_kind not in {"X", "layer"}:
        raise ValueError(f"{operation} requires an explicit X or layer expression source.")
    expected_source = table_parameters["marker_source"]
    expected_layer = table_parameters.get("marker_layer_name")
    if source_kind != expected_source or (source_kind == "layer" and layer_name != expected_layer):
        raise ValueError(f"{operation} expression source must match the upstream Marker Genes source.")
    if source_kind == "X":
        if layer_name is not None:
            raise ValueError(f"{operation} X source cannot carry a layer name.")
        matrix = adata.X
    else:
        if not isinstance(layer_name, str) or not layer_name or layer_name != layer_name.strip():
            raise ValueError(f"{operation} layer source requires a canonical layer name.")
        if layer_name not in adata.layers:
            raise ValueError(f"{operation} expression layer not found: {layer_name!r}.")
        matrix = adata.layers[layer_name]
    if list(adata.var_names) != universe_genes:
        raise ValueError(f"{operation} feature axis does not match the tested-gene universe.")
    if tuple(matrix.shape) != (int(adata.n_obs), len(universe_genes)):
        raise ValueError(f"{operation} expression matrix is not aligned to its axes.")
    matrix_values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    matrix_values = np.asarray(matrix_values)
    if (
        not np.issubdtype(matrix_values.dtype, np.number)
        or np.issubdtype(matrix_values.dtype, np.bool_)
        or np.iscomplexobj(matrix_values)
    ):
        raise TypeError(f"{operation} expression must contain real numeric values.")
    if matrix_values.size and not bool(np.isfinite(matrix_values).all()):
        raise ValueError(f"{operation} expression contains non-finite values.")

    groupby = table_parameters["marker_groupby"]
    if not isinstance(groupby, str) or not groupby or groupby not in adata.obs:
        raise ValueError(f"{operation} grouping column does not match marker provenance.")
    grouping = adata.obs[groupby]
    if not isinstance(grouping.dtype, pd.CategoricalDtype) or bool(grouping.isna().any()):
        raise ValueError(f"{operation} grouping must remain a complete categorical column.")
    raw_categories = grouping.cat.categories.tolist()
    group_labels = [str(value) for value in raw_categories]
    if group_labels != list(table_parameters["group_labels"]):
        raise ValueError(f"{operation} group axis does not match marker provenance.")
    codes = grouping.cat.codes.to_numpy(dtype=int)
    labels_by_observation = [group_labels[int(code)] for code in codes]
    group_sizes = {
        label: int(np.count_nonzero(codes == index)) for index, label in enumerate(group_labels)
    }
    if group_sizes != dict(table_parameters["group_sizes"]):
        raise ValueError(f"{operation} group sizes do not match marker provenance.")
    analysis_fingerprint = canonical_sha256(
        {
            "schema": "openbio-singlecell/cluster-marker-analysis/v2",
            "groupby": groupby,
            "method": table_parameters["marker_method"],
            "source": source_kind,
            "layer_name": layer_name if source_kind == "layer" else None,
            "n_genes": table_parameters["requested_n_genes"],
            "tie_correct": table_parameters["tie_correct"],
            "universe_fingerprint": universe_fingerprint,
            "ordered_observation_ids": list(adata.obs_names),
            "group_labels_by_observation": labels_by_observation,
            "fixed_policy": {
                "groups": "all",
                "reference": "rest",
                "rankby_abs": False,
                "pts": True,
                "corr_method": "benjamini-hochberg",
            },
        }
    )
    if table_parameters["analysis_fingerprint"] != analysis_fingerprint:
        raise ValueError(f"{operation} observation axis or group-label alignment does not match provenance.")

    if isinstance(top_genes_per_group, (bool, np.bool_)) or not isinstance(
        top_genes_per_group, (int, np.integer)
    ):
        raise TypeError(f"{operation} top_genes_per_group must be an integer.")
    top_genes_per_group = int(top_genes_per_group)
    if not 1 <= top_genes_per_group <= 20:
        raise ValueError(f"{operation} top_genes_per_group must be between 1 and 20.")
    rows_by_group = {
        group: table.loc[table["group"] == group].iloc[:top_genes_per_group]
        for group in group_labels
    }
    selected_upstream_ranks = {
        group: [
            {"gene": str(row.gene), "rank": int(row.rank)}
            for row in rows.itertuples(index=False)
        ]
        for group, rows in rows_by_group.items()
    }
    genes = list(
        dict.fromkeys(
            row["gene"]
            for group in group_labels
            for row in selected_upstream_ranks[group]
        )
    )
    if not genes:
        raise ValueError(f"{operation} found no upstream-ranked genes to plot.")
    if len(group_labels) > 50 or len(genes) > 50:
        raise ValueError(
            f"{operation} supports at most 50 groups and 50 unique genes in one static PNG; "
            f"received {len(group_labels)} groups and {len(genes)} genes."
        )
    indices = [universe_genes.index(gene) for gene in genes]
    selected = matrix[:, indices]
    panel = selected.toarray() if sparse.issparse(selected) else np.asarray(selected).copy()
    panel = np.asarray(panel, dtype=float)
    means = np.vstack([panel[codes == index].mean(axis=0) for index in range(len(group_labels))])
    fractions = np.vstack(
        [(panel[codes == index] > 0.0).mean(axis=0) for index in range(len(group_labels))]
    )

    width = min(20.0, max(7.0, 3.0 + 0.55 * len(genes)))
    height = min(16.0, max(4.5, 2.5 + 0.42 * len(group_labels)))
    figure = Figure(figsize=(width, height), constrained_layout=True)
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    x_values = np.tile(np.arange(len(genes)), len(group_labels))
    y_values = np.repeat(np.arange(len(group_labels)), len(genes))
    points = axis.scatter(
        x_values,
        y_values,
        c=means.ravel(),
        s=20.0 + 280.0 * fractions.ravel(),
        cmap="viridis",
        edgecolors="#333333",
        linewidths=0.25,
    )
    axis.set_xticks(np.arange(len(genes)), genes, rotation=45, ha="right")
    axis.set_yticks(np.arange(len(group_labels)), group_labels)
    axis.invert_yaxis()
    axis.set_xlabel("Upstream-ranked genes")
    axis.set_ylabel(groupby)
    axis.set_title("Cluster marker evidence: expression and detected fraction")
    axis.grid(alpha=0.15)
    colorbar = figure.colorbar(points, ax=axis)
    colorbar.set_label("Mean selected-source expression")
    size_legend_fractions = [0.25, 0.5, 1.0]
    size_handles = [
        axis.scatter(
            [],
            [],
            s=20.0 + 280.0 * fraction,
            facecolors="none",
            edgecolors="#333333",
            linewidths=0.6,
        )
        for fraction in size_legend_fractions
    ]
    figure.legend(
        size_handles,
        [f"{fraction:.0%}" for fraction in size_legend_fractions],
        title="Detected fraction (> 0)",
        loc="outside lower center",
        ncols=3,
        frameon=True,
    )
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "genes": genes,
        "groups": group_labels,
        "groupby": groupby,
        "group_cell_counts": group_sizes,
        "top_genes_per_group": top_genes_per_group,
        "selected_upstream_ranks": selected_upstream_ranks,
        "plotted_groups": len(group_labels),
        "plotted_genes": len(genes),
        "plotted_cells": int(adata.n_obs),
        "expression_source": {
            "source": source_kind,
            "layer_name": layer_name if source_kind == "layer" else None,
        },
        "table_operation": producer["table_operation"],
        "analysis_fingerprint_sha256": analysis_fingerprint,
        "marker_table_content_fingerprint_sha256": table_fingerprint,
        "universe_content_fingerprint_sha256": universe_content,
        "display_semantics": {
            "color": "mean selected-source expression",
            "size": "fraction of cells with expression > 0",
        },
        "size_legend_fractions": size_legend_fractions,
        "size_legend_location": "outside lower center",
        "mean_expression": means.tolist(),
        "detected_fraction": fractions.tolist(),
        "rendering": {
            "figure_size_inches": [width, height],
            "dpi": 120,
            "maximum_groups": 50,
            "maximum_unique_genes": 50,
            "maximum_cells": 1_000_000,
        },
    }
    return (png, details) if _return_details else png


def marker_evidence_plot_code(
    *,
    producer: dict[str, Any],
    source_kind: str,
    layer_name: str | None,
    top_genes_per_group: int,
) -> str:
    implementation = textwrap.dedent(inspect.getsource(_standalone_marker_evidence_plot)).strip()
    return f"""from __future__ import annotations

{implementation}


def plot_marker_evidence(adata, table, universe):
    return _standalone_marker_evidence_plot(
        adata,
        table,
        universe,
        producer={producer!r},
        source_kind={source_kind!r},
        layer_name={layer_name!r},
        top_genes_per_group={top_genes_per_group!r},
    )
"""


__all__ = ["_standalone_marker_evidence_plot", "marker_evidence_plot_code"]
