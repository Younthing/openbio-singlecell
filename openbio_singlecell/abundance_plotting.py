from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .composition_result import (
    _HIERARCHY_FIELDS as _COMPOSITION_HIERARCHY_FIELDS,
)
from .composition_result import (
    _MODEL_METADATA_FIELDS as _COMPOSITION_MODEL_METADATA_FIELDS,
)
from .composition_result import (
    _POSTERIOR_VARIABLES as _COMPOSITION_POSTERIOR_VARIABLES,
)
from .composition_result import (
    _SAMPLE_STATS as _COMPOSITION_SAMPLE_STATS,
)
from .composition_result import (
    COMPOSITION_MODEL_ARTIFACT_SCHEMA_VERSION,
    COMPOSITION_MODEL_ARTIFACT_TYPE,
    SCCODA_RESULT_COLUMNS,
    TASCCODA_RESULT_COLUMNS,
    CompositionModelResult,
    validate_composition_model_result,
)
from .composition_result import (
    _array_identities as _composition_array_identities,
)
from .composition_result import (
    _axis as _composition_axis,
)
from .composition_result import (
    _canonical_arrays as _composition_canonical_arrays,
)
from .composition_result import (
    _canonical_json as _composition_canonical_json,
)
from .composition_result import (
    _json_sha256 as _composition_json_sha256,
)
from .composition_result import (
    _sample_stat_arrays as _composition_sample_stat_arrays,
)
from .composition_result import (
    _table_identity as _composition_table_identity,
)
from .composition_result import (
    _validate_hierarchy as _composition_validate_hierarchy,
)
from .composition_result import (
    _validate_model_metadata as _composition_validate_model_metadata,
)
from .milo_result import (
    MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS,
    MILO_MAX_OVERLAP_WORKING_BYTES,
    MILO_RESULT_ARTIFACT_SCHEMA_VERSION,
    MILO_RESULT_ARTIFACT_TYPE,
    MILO_RESULT_COLUMNS,
    MILO_RESULT_PRODUCER_NODE_ID,
    MILO_RESULT_PRODUCER_SCHEMA,
    MiloResult,
    _axis,
    _coordinate_sha256,
    _coordinates,
    _csr,
    _json_sha256,
    _matrix_sha256,
    _membership_overlap_graph,
    _table_identity,
    validate_milo_result,
)

_MILO_TEXT_COLUMNS = {
    "neighborhood_id",
    "index_cell",
    "reference_condition",
    "comparison_condition",
    "majority_annotation",
    "reported_annotation",
}
_MILO_UNIT_COLUMNS = {"p_value", "p_adjusted_bh", "spatial_fdr", "majority_annotation_fraction"}


def _standalone_milo_differential_abundance_plot(result, *, view=None, _return_diagnostics=False):
    import io

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Milo Differential Abundance Plot"
    if view is None:
        mode = {"view": "differential_evidence"}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    if set(mode) != {"view"} or mode.get("view") not in {"differential_evidence", "neighborhood_graph"}:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    portable = result.portable() if callable(getattr(result, "portable", None)) else result
    table, membership, graph, coordinates, observations, neighborhoods, provenance, metadata = validate_milo_result(
        portable,
        exact_type=False,
        copy_result=False,
    )
    if len(neighborhoods) > 100_000:
        raise ValueError(f"{operation} supports at most 100,000 neighborhoods in one static PNG.")
    reference_conditions = table["reference_condition"].unique().tolist()
    comparison_conditions = table["comparison_condition"].unique().tolist()
    if (
        len(reference_conditions) != 1
        or len(comparison_conditions) != 1
        or reference_conditions[0] == comparison_conditions[0]
    ):
        raise ValueError(f"{operation} table must retain one explicit pairwise Condition direction.")

    effects = table["log2_fold_change"].to_numpy(dtype=float)
    spatial_fdr = table["spatial_fdr"].to_numpy(dtype=float)
    sizes = membership.sum(axis=0).A1.astype(int)
    selected = (spatial_fdr <= float(provenance["spatial_fdr_threshold"])) & (
        np.abs(effects) >= float(provenance["min_abs_log2_fold_change"])
    )
    annotations = [str(value) for value in pd.unique(table["majority_annotation"])]
    if mode["view"] == "differential_evidence" and len(annotations) > 100:
        raise ValueError(f"{operation} differential_evidence supports at most 100 majority annotations.")
    annotation_positions = {value: index for index, value in enumerate(annotations)}
    colors = [f"C{annotation_positions[value] % 10}" for value in table["majority_annotation"]]
    point_sizes = 20.0 + 90.0 * np.sqrt(sizes / max(1, int(sizes.max())))

    if mode["view"] == "differential_evidence":
        figure = Figure(figsize=(8.0, 5.8), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        y = -np.log10(np.maximum(spatial_fdr, 1e-300))
        axis.scatter(effects, y, s=point_sizes, c=colors, alpha=0.75, linewidths=0)
        if bool(selected.any()):
            axis.scatter(
                effects[selected],
                y[selected],
                s=point_sizes[selected] + 20,
                facecolors="none",
                edgecolors="#111111",
                linewidths=0.9,
            )
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        threshold = float(provenance["spatial_fdr_threshold"])
        axis.axhline(-np.log10(max(threshold, 1e-300)), color="#555555", linestyle="--", linewidth=0.9)
        minimum_effect = float(provenance["min_abs_log2_fold_change"])
        if minimum_effect > 0:
            axis.axvline(minimum_effect, color="#777777", linestyle="--", linewidth=0.8)
            axis.axvline(-minimum_effect, color="#777777", linestyle="--", linewidth=0.8)
        for annotation in annotations:
            axis.scatter([], [], color=f"C{annotation_positions[annotation] % 10}", label=annotation)
        axis.set_xlabel(
            f"Neighborhood log2 fold change ({table['comparison_condition'].iloc[0]} − "
            f"{table['reference_condition'].iloc[0]})"
        )
        axis.set_ylabel("−log10(Milo Spatial FDR)")
        axis.set_title("Milo neighborhood differential-abundance evidence")
        axis.grid(alpha=0.2)
        axis.legend(title=f"Majority {provenance['annotation_key']}", bbox_to_anchor=(1.02, 1.0), loc="upper left")
        title = "Milo neighborhood differential-abundance evidence"
        plotted_edges = 0
    else:
        undirected_edges = int(graph.nnz // 2)
        if len(neighborhoods) > 10_000 or undirected_edges > 50_000:
            raise ValueError(
                f"{operation} neighborhood_graph supports at most 10,000 neighborhoods and 50,000 overlap edges."
            )
        edge_rows, edge_columns = graph.nonzero()
        edge_pairs = [
            (int(row), int(column))
            for row, column in zip(edge_rows, edge_columns, strict=True)
            if row < column
        ]
        figure = Figure(figsize=(8.0, 6.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        maximum_overlap = max((float(graph[row, column]) for row, column in edge_pairs), default=1.0)
        for row, column in edge_pairs:
            overlap = float(graph[row, column])
            axis.plot(
                coordinates[[row, column], 0],
                coordinates[[row, column], 1],
                color="#A8A8A8",
                alpha=0.25 + 0.45 * overlap / maximum_overlap,
                linewidth=0.4 + 1.2 * overlap / maximum_overlap,
                zorder=1,
            )
        image = axis.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            c=effects,
            cmap="coolwarm",
            s=point_sizes,
            alpha=0.85,
            linewidths=np.where(selected, 1.0, 0.0),
            edgecolors="#111111",
            zorder=2,
        )
        axis.set_xlabel(f"{provenance['representation_key']} component 1")
        axis.set_ylabel(f"{provenance['representation_key']} component 2")
        axis.set_title("Milo retained neighborhood-overlap graph")
        axis.grid(alpha=0.15)
        figure.colorbar(image, ax=axis, label="Neighborhood log2 fold change")
        title = "Milo neighborhood graph"
        plotted_edges = len(edge_pairs)

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": mode["view"],
        "view_parameters": mode,
        "neighborhoods": len(neighborhoods),
        "observations": len(observations),
        "plotted_edges": plotted_edges,
        "spatial_fdr_selected": int(selected.sum()),
        "positive_selected": int((selected & (effects > 0)).sum()),
        "negative_selected": int((selected & (effects < 0)).sum()),
        "reference_condition": str(table["reference_condition"].iloc[0]),
        "comparison_condition": str(table["comparison_condition"].iloc[0]),
        "condition_key": provenance["condition_key"],
        "annotation_key": provenance["annotation_key"],
        "annotation_status": provenance["annotation_status"],
        "representation_key": provenance["representation_key"],
        "spatial_fdr_threshold": float(provenance["spatial_fdr_threshold"]),
        "minimum_absolute_log2_fold_change": float(provenance["min_abs_log2_fold_change"]),
        "membership_nonzero": int(membership.nnz),
        "artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
        "title": title,
    }
    if _return_diagnostics:
        return png, details
    return png


def run_milo_differential_abundance_plot(result: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    if type(result) is not MiloResult:
        raise TypeError("Milo Differential Abundance Plot requires an exact OPENBIO_MILO_RESULT artifact.")
    return _standalone_milo_differential_abundance_plot(result, view=view, _return_diagnostics=True)


def milo_differential_abundance_plot_code(*, view: dict[str, Any]) -> str:
    functions = (
        _json_sha256,
        _axis,
        _csr,
        _matrix_sha256,
        _membership_overlap_graph,
        _coordinates,
        _coordinate_sha256,
        _table_identity,
        validate_milo_result,
        _standalone_milo_differential_abundance_plot,
    )
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence

MILO_RESULT_ARTIFACT_TYPE = {MILO_RESULT_ARTIFACT_TYPE!r}
MILO_RESULT_ARTIFACT_SCHEMA_VERSION = {MILO_RESULT_ARTIFACT_SCHEMA_VERSION!r}
MILO_RESULT_PRODUCER_NODE_ID = {MILO_RESULT_PRODUCER_NODE_ID!r}
MILO_RESULT_PRODUCER_SCHEMA = {MILO_RESULT_PRODUCER_SCHEMA!r}
MILO_RESULT_COLUMNS = {MILO_RESULT_COLUMNS!r}
MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS = {MILO_MAX_OVERLAP_PAIR_CONTRIBUTIONS!r}
MILO_MAX_OVERLAP_WORKING_BYTES = {MILO_MAX_OVERLAP_WORKING_BYTES!r}
_TEXT_COLUMNS = {_MILO_TEXT_COLUMNS!r}
_UNIT_COLUMNS = {_MILO_UNIT_COLUMNS!r}

{implementation}


def plot_milo_differential_abundance(result):
    return _standalone_milo_differential_abundance_plot(result, view={view!r})
"""


def _standalone_composition_differential_plot(
    result,
    *,
    expected_method,
    view=None,
    _return_diagnostics=False,
):
    import io

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    display_method = "scCODA" if expected_method == "sccoda" else "tascCODA"
    operation = f"{display_method} Differential Composition Plot"
    default_view = "effect_forest" if expected_method == "sccoda" else "hierarchy_effect_forest"
    if view is None:
        mode = {"view": default_view}
    elif not isinstance(view, dict):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if expected_method == "sccoda":
        if variant in {"effect_forest", "sampler_diagnostics"}:
            expected_view = {"view"}
        elif variant == "posterior_distribution":
            expected_view = {"view", "cell_type"}
        else:
            expected_view = set()
    else:
        if variant in {"hierarchy_effect_forest", "derived_leaf_effects", "sampler_diagnostics"}:
            expected_view = {"view"}
        elif variant == "posterior_distribution":
            expected_view = {"view", "effect_scope", "effect_name"}
        else:
            expected_view = set()
    if set(mode) != expected_view:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")

    portable = result.portable() if callable(getattr(result, "portable", None)) else result
    table, posterior, sample_stats, cell_types, hierarchy, model_metadata, metadata = (
        validate_composition_model_result(portable, exact_type=False, copy_result=False)
    )
    if metadata["method"] != expected_method:
        raise ValueError(f"{operation} requires the exact {display_method} typed result.")
    if len(table) > 100_000 or len(cell_types) > 10_000:
        raise ValueError(f"{operation} result exceeds the fixed static-plot evidence limit.")
    draw_count = int(metadata["draw_count"])
    if draw_count > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 retained posterior draws.")
    reference = model_metadata["reference_cell_type"]
    reference_constraint = table["reference_constraint"]
    if bool(reference_constraint.isna().any()) or any(
        not isinstance(value, (bool, np.bool_)) for value in reference_constraint
    ):
        raise TypeError(f"{operation} table column 'reference_constraint' must contain complete Boolean values.")

    def finite_numeric(frame, columns):
        values = {}
        for column in columns:
            series = frame[column]
            if pd.api.types.is_bool_dtype(series.dtype):
                raise TypeError(f"{operation} table column {column!r} must be numeric and non-Boolean.")
            try:
                numeric = series.to_numpy(dtype=float, na_value=np.nan)
            except (TypeError, ValueError) as error:
                raise TypeError(f"{operation} table column {column!r} must be numeric and non-Boolean.") from error
            if not bool(np.isfinite(numeric).all()):
                raise ValueError(f"{operation} table column {column!r} must contain finite values.")
            values[column] = numeric
        return values

    if expected_method == "sccoda":
        credible = table["credible_effect"]
        if bool(credible.isna().any()) or any(not isinstance(value, (bool, np.bool_)) for value in credible):
            raise TypeError(f"{operation} table column 'credible_effect' must contain complete Boolean values.")
        numeric = finite_numeric(
            table,
            (
            "model_coefficient",
            "hdi_lower",
            "hdi_upper",
            "posterior_sd",
            "inclusion_probability",
            "estimated_fdr",
            "inclusion_probability_threshold",
            "realized_expected_fdr",
            "expected_count_reference",
            "expected_count_comparison",
            "compositional_log2_fold_change",
            ),
        )
        expected_reference = numeric["expected_count_reference"]
        expected_comparison = numeric["expected_count_comparison"]
    else:
        direct = table.loc[table["effect_scope"] == "hierarchy_node"]
        leaves = table.loc[table["effect_scope"] == "derived_leaf"]
        if direct.empty or leaves.empty or len(direct) + len(leaves) != len(table):
            raise ValueError(f"{operation} requires complete hierarchy-node and derived-leaf scopes.")
        direct_credible = direct["credible_effect"]
        if bool(direct_credible.isna().any()) or any(
            not isinstance(value, (bool, np.bool_)) for value in direct_credible
        ):
            raise TypeError(f"{operation} hierarchy-node credible_effect values must be complete Booleans.")
        if not bool(leaves["credible_effect"].isna().all()):
            raise ValueError(f"{operation} derived-leaf credible_effect values must be null by contract.")
        if not bool(leaves["selection_delta"].isna().all()):
            raise ValueError(f"{operation} derived-leaf selection_delta values must be null by contract.")
        for column in ("expected_count_reference", "expected_count_comparison", "compositional_log2_fold_change"):
            if not bool(direct[column].isna().all()):
                raise ValueError(f"{operation} hierarchy-node {column} values must be null by contract.")
        if direct["selection_basis"].tolist() != [
            "direct_node_abs_median_strictly_greater_than_delta"
        ] * len(direct) or leaves["selection_basis"].tolist() != [
            "derived_sum_of_selected_hierarchy_node_effects"
        ] * len(leaves):
            raise ValueError(f"{operation} direct/derived selection-basis semantics are inconsistent.")
        if any(
            not isinstance(value, str) or not value or value != value.strip() for value in table["hierarchy_level"]
        ):
            raise ValueError(f"{operation} hierarchy levels must be canonical nonblank keys.")
        numeric = finite_numeric(
            table,
            (
            "descendant_leaf_count",
            "model_effect",
            "posterior_median",
            "hdi_lower",
            "hdi_upper",
            "posterior_sd",
            ),
        )
        selection_delta = finite_numeric(direct, ("selection_delta",))["selection_delta"]
        if bool((selection_delta < 0).any()):
            raise ValueError(f"{operation} hierarchy-node selection deltas cannot be negative.")
        leaf_numeric = finite_numeric(
            leaves,
            (
            "expected_count_reference",
            "expected_count_comparison",
            "compositional_log2_fold_change",
            ),
        )
        expected_reference = leaf_numeric["expected_count_reference"]
        expected_comparison = leaf_numeric["expected_count_comparison"]
    if bool((numeric["posterior_sd"] < 0).any()):
        raise ValueError(f"{operation} posterior standard deviations cannot be negative.")
    if bool((expected_reference < 0).any()) or bool((expected_comparison < 0).any()):
        raise ValueError(f"{operation} expected composition counts cannot be negative.")

    selected_effect_name = None
    selected_effect_scope = None
    credible_effects = int(
        table["credible_effect"].sum()
        if expected_method == "sccoda"
        else table.loc[table["effect_scope"] == "hierarchy_node", "credible_effect"].sum()
    )
    if variant in {"effect_forest", "hierarchy_effect_forest"}:
        selected = table if expected_method == "sccoda" else table.loc[table["effect_scope"] == "hierarchy_node"]
        if len(selected) > 100:
            raise ValueError(f"{operation} forest view supports at most 100 modeled effects.")
        labels = (
            selected["cell_type"].astype(str).tolist()
            if expected_method == "sccoda"
            else selected["effect_name"].astype(str).tolist()
        )
        centers = (
            selected["model_coefficient"].to_numpy(dtype=float)
            if expected_method == "sccoda"
            else selected["posterior_median"].to_numpy(dtype=float)
        )
        lower = selected["hdi_lower"].to_numpy(dtype=float)
        upper = selected["hdi_upper"].to_numpy(dtype=float)
        if bool((lower > upper).any()):
            raise ValueError(f"{operation} HDI lower bounds cannot exceed upper bounds.")
        credible = selected["credible_effect"].to_numpy(dtype=bool)
        figure = Figure(figsize=(8.5, max(4.8, min(22.0, 0.35 * len(selected) + 2.5))), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        positions = np.arange(len(selected), dtype=float)
        for position, center, low, high, is_credible in zip(
            positions, centers, lower, upper, credible, strict=True
        ):
            color = "#D64A4A" if is_credible and center > 0 else "#3D6FB6" if is_credible else "#8A8A8A"
            axis.plot([low, high], [position, position], color=color, marker="|", markersize=6)
            axis.plot(center, position, "o", color=color)
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.set_xlabel("Model coefficient with posterior HDI" if expected_method == "sccoda" else "Hierarchy-node posterior effect with HDI")
        axis.set_title(f"{display_method} relative compositional effects")
        axis.grid(axis="x", alpha=0.2)
        title = f"{display_method} relative effect forest"
        plotted_effects = len(selected)
    elif variant == "derived_leaf_effects":
        selected = table.loc[table["effect_scope"] == "derived_leaf"]
        if len(selected) > 100:
            raise ValueError(f"{operation} derived-leaf view supports at most 100 cell types.")
        labels = selected["effect_name"].astype(str).tolist()
        effects = selected["compositional_log2_fold_change"].to_numpy(dtype=float)
        colors = np.where(effects > 0, "#D64A4A", np.where(effects < 0, "#3D6FB6", "#8A8A8A"))
        figure = Figure(figsize=(8.5, max(4.8, min(22.0, 0.35 * len(selected) + 2.5))), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        positions = np.arange(len(selected), dtype=float)
        axis.barh(positions, effects, color=colors)
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.set_xlabel(
            f"Derived compositional log2 fold change ({model_metadata['comparison_condition']} − "
            f"{model_metadata['reference_condition']})"
        )
        axis.set_title("tascCODA propagated leaf effects")
        axis.grid(axis="x", alpha=0.2)
        title = "tascCODA derived leaf effects"
        plotted_effects = len(selected)
    elif variant == "posterior_distribution":
        if expected_method == "sccoda":
            requested = mode["cell_type"]
            if not isinstance(requested, str) or requested != requested.strip():
                raise ValueError(f"{operation} cell_type must be canonical text or empty for automatic selection.")
            candidates = [
                value
                for value in table.loc[table["credible_effect"] & ~table["reference_constraint"], "cell_type"]
                if value != reference
            ]
            candidates.extend(value for value in cell_types if value != reference and value not in candidates)
            if requested:
                if requested not in cell_types:
                    raise ValueError(f"{operation} cell_type {requested!r} is absent from the retained posterior axis.")
                selected_effect_name = requested
            elif candidates:
                selected_effect_name = candidates[0]
            else:
                selected_effect_name = reference
            selected_effect_scope = "cell_type"
            draws = posterior["condition_effect"][0, :, cell_types.index(selected_effect_name)]
        else:
            scope = mode["effect_scope"]
            requested = mode["effect_name"]
            if scope not in {"hierarchy_node", "derived_leaf"}:
                raise ValueError(f"{operation} effect_scope must be hierarchy_node or derived_leaf.")
            if not isinstance(requested, str) or requested != requested.strip():
                raise ValueError(f"{operation} effect_name must be canonical text or empty for automatic selection.")
            if scope == "hierarchy_node":
                axis_names = hierarchy["node_names"]
                reference_names = set(hierarchy["reference_nodes"])
                variable = "hierarchy_node_effect"
            else:
                axis_names = cell_types
                reference_names = {reference}
                variable = "derived_leaf_effect"
            candidates = [name for name in axis_names if name not in reference_names]
            if requested:
                if requested not in axis_names:
                    raise ValueError(f"{operation} effect_name {requested!r} is absent from the retained {scope} axis.")
                selected_effect_name = requested
            elif candidates:
                selected_effect_name = candidates[0]
            else:
                selected_effect_name = axis_names[0]
            selected_effect_scope = scope
            draws = posterior[variable][0, :, axis_names.index(selected_effect_name)]
        figure = Figure(figsize=(9.5, 4.8), constrained_layout=True)
        FigureCanvasAgg(figure)
        trace_axis, density_axis = figure.subplots(1, 2)
        draw_indices = np.arange(draw_count, dtype=int)
        trace_axis.plot(draw_indices, draws, color="#4C78A8", linewidth=0.8)
        trace_axis.axhline(0.0, color="#555555", linewidth=0.8)
        trace_axis.set_xlabel("Retained posterior draw")
        trace_axis.set_ylabel("Relative compositional effect")
        trace_axis.set_title(f"{selected_effect_name} posterior trace")
        trace_axis.grid(alpha=0.2)
        density_axis.hist(draws, bins=min(50, max(5, int(np.ceil(np.sqrt(draw_count))))), color="#4C78A8", alpha=0.8)
        density_axis.axvline(0.0, color="#555555", linewidth=0.8)
        density_axis.set_xlabel("Relative compositional effect")
        density_axis.set_ylabel("Posterior draws")
        density_axis.set_title(f"{selected_effect_name} posterior distribution")
        density_axis.grid(axis="y", alpha=0.2)
        title = f"{display_method} posterior distribution: {selected_effect_name}"
        plotted_effects = 1
    else:
        figure = Figure(figsize=(9.0, 8.0), constrained_layout=True)
        FigureCanvasAgg(figure)
        axes = np.atleast_1d(figure.subplots(3, 1, squeeze=False)).reshape(-1)
        draw_indices = np.arange(draw_count, dtype=int)
        labels = {
            "potential_energy": "Potential energy",
            "num_steps": "NUTS steps",
            "step_size": "NUTS step size",
        }
        for axis, name in zip(axes, ("potential_energy", "num_steps", "step_size"), strict=True):
            values = sample_stats[name][0]
            axis.plot(draw_indices, values, color="#4C78A8", linewidth=0.8)
            axis.set_ylabel(labels[name])
            axis.grid(alpha=0.2)
        axes[-1].set_xlabel("Retained posterior draw")
        axes[0].set_title(f"{display_method} retained one-chain sampler diagnostics")
        title = f"{display_method} sampler diagnostics"
        plotted_effects = 0

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": variant,
        "view_parameters": mode,
        "method": display_method,
        "reference_condition": model_metadata["reference_condition"],
        "comparison_condition": model_metadata["comparison_condition"],
        "condition_key": model_metadata["condition_key"],
        "reference_cell_type": reference,
        "annotation_key": model_metadata["annotation_key"],
        "annotation_status": model_metadata["annotation_status"],
        "cell_types": cell_types,
        "posterior_draws": draw_count,
        "chain_count": int(metadata["chain_count"]),
        "rhat_available": bool(model_metadata["diagnostics"]["rhat_available"]),
        "divergences_available": bool(model_metadata["diagnostics"]["divergences_available"]),
        "credible_effects": credible_effects,
        "derived_leaf_credibility": (
            "not_applicable_propagated_effects_are_not_direct_selections"
            if expected_method == "tasccoda"
            else None
        ),
        "plotted_effects": plotted_effects,
        "selected_effect_scope": selected_effect_scope,
        "selected_effect_name": selected_effect_name,
        "artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
        "title": title,
    }
    if variant == "posterior_distribution":
        details["posterior_distribution"] = {
            "minimum": float(np.min(draws)),
            "median": float(np.median(draws)),
            "maximum": float(np.max(draws)),
        }
    if variant == "sampler_diagnostics":
        details["sample_stat_ranges"] = {
            name: {"minimum": float(values.min()), "maximum": float(values.max())}
            for name, values in sample_stats.items()
        }
    if _return_diagnostics:
        return png, details
    return png


def _run_composition_differential_plot(
    result: Any,
    *,
    expected_method: str,
    view: Any = None,
) -> tuple[bytes, dict[str, Any]]:
    if type(result) is not CompositionModelResult:
        raise TypeError(
            "Differential Composition Plot requires an exact OPENBIO_COMPOSITION_MODEL_RESULT artifact."
        )
    return _standalone_composition_differential_plot(
        result,
        expected_method=expected_method,
        view=view,
        _return_diagnostics=True,
    )


def run_sccoda_differential_composition_plot(result: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _run_composition_differential_plot(result, expected_method="sccoda", view=view)


def run_tasccoda_differential_composition_plot(result: Any, *, view: Any = None) -> tuple[bytes, dict[str, Any]]:
    return _run_composition_differential_plot(result, expected_method="tasccoda", view=view)


def _composition_differential_plot_code(
    *,
    expected_method: str,
    function_name: str,
    view: dict[str, Any],
) -> str:
    functions = (
        _composition_json_sha256,
        _composition_axis,
        _composition_canonical_json,
        _composition_table_identity,
        _composition_canonical_arrays,
        _composition_sample_stat_arrays,
        _composition_array_identities,
        _composition_validate_model_metadata,
        _composition_validate_hierarchy,
        validate_composition_model_result,
        _standalone_composition_differential_plot,
    )
    implementation = "\n\n".join(textwrap.dedent(inspect.getsource(function)).strip() for function in functions)
    return f"""from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence

COMPOSITION_MODEL_ARTIFACT_TYPE = {COMPOSITION_MODEL_ARTIFACT_TYPE!r}
COMPOSITION_MODEL_ARTIFACT_SCHEMA_VERSION = {COMPOSITION_MODEL_ARTIFACT_SCHEMA_VERSION!r}
SCCODA_RESULT_COLUMNS = {SCCODA_RESULT_COLUMNS!r}
TASCCODA_RESULT_COLUMNS = {TASCCODA_RESULT_COLUMNS!r}
_POSTERIOR_VARIABLES = {_COMPOSITION_POSTERIOR_VARIABLES!r}
_SAMPLE_STATS = {_COMPOSITION_SAMPLE_STATS!r}
_MODEL_METADATA_FIELDS = {_COMPOSITION_MODEL_METADATA_FIELDS!r}
_HIERARCHY_FIELDS = {_COMPOSITION_HIERARCHY_FIELDS!r}

{implementation}


def {function_name}(result):
    return _standalone_composition_differential_plot(
        result,
        expected_method={expected_method!r},
        view={view!r},
    )
"""


def sccoda_differential_composition_plot_code(*, view: dict[str, Any]) -> str:
    return _composition_differential_plot_code(
        expected_method="sccoda",
        function_name="plot_sccoda_differential_composition",
        view=view,
    )


def tasccoda_differential_composition_plot_code(*, view: dict[str, Any]) -> str:
    return _composition_differential_plot_code(
        expected_method="tasccoda",
        function_name="plot_tasccoda_differential_composition",
        view=view,
    )


__all__ = [
    "milo_differential_abundance_plot_code",
    "run_milo_differential_abundance_plot",
    "run_sccoda_differential_composition_plot",
    "run_tasccoda_differential_composition_plot",
    "sccoda_differential_composition_plot_code",
    "tasccoda_differential_composition_plot_code",
]
