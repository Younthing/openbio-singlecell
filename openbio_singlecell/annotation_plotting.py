from __future__ import annotations

import inspect
import textwrap
from typing import Any

from .annotation_core import (
    _celltypist_axis_fingerprint,
    _celltypist_probability_fingerprint,
    _celltypist_selected_fingerprint,
)


def _standalone_celltypist_diagnostics_plot(
    adata,
    *,
    metadata_key="celltypist",
    view=None,
    _return_details=False,
):
    import io
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "CellTypist Diagnostics Plot"
    if view is None:
        mode = {"view": "confidence_distributions", "max_labels": 30}
    elif not isinstance(view, Mapping):
        raise TypeError(f"{operation} view must be a DynamicCombo value.")
    else:
        mode = dict(view)
    variant = mode.get("view")
    if variant == "confidence_distributions":
        if set(mode) != {"view", "max_labels"}:
            raise ValueError(f"{operation} confidence view received inactive or unknown parameters.")
        limit_name = "max_labels"
    elif variant == "probability_heatmap":
        if set(mode) != {"view", "groupby", "max_classes"}:
            raise ValueError(f"{operation} probability heatmap received inactive or unknown parameters.")
        limit_name = "max_classes"
    else:
        raise ValueError(f"{operation} view must be a closed supported DynamicCombo selection.")
    limit = mode[limit_name]
    if isinstance(limit, (bool, np.bool_)) or not isinstance(limit, (int, np.integer)):
        raise TypeError(f"{operation} {limit_name} must be an integer.")
    limit = int(limit)
    if not 1 <= limit <= 100:
        raise ValueError(f"{operation} {limit_name} must be between 1 and 100.")
    mode[limit_name] = limit

    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires in-memory AnnData.")
    if int(getattr(adata, "n_obs", 0)) < 1 or not bool(adata.obs_names.is_unique):
        raise ValueError(f"{operation} requires a non-empty unique observation axis.")
    if int(adata.n_obs) > 1_000_000:
        raise ValueError(f"{operation} supports at most 1,000,000 cells in one static PNG.")
    if not isinstance(metadata_key, str) or not metadata_key or metadata_key != metadata_key.strip():
        raise ValueError(f"{operation} metadata_key must be a canonical nonempty string.")
    provenance = adata.uns.get(metadata_key)
    required_provenance = {
        "schema_version",
        "operation",
        "annotation_status",
        "expression",
        "model",
        "feature_overlap",
        "classes",
        "fixed_policy",
        "majority_voting",
        "over_clustering_key",
        "min_prop",
        "overwrite_existing",
        "columns",
        "matrices",
        "artifact_replacement",
        "software_versions",
        "warnings",
        "observation_axis_fingerprint_sha256",
        "probability_content_fingerprint_sha256",
        "selected_annotation_fingerprint_sha256",
    }
    if not isinstance(provenance, Mapping) or not required_provenance.issubset(provenance):
        raise ValueError(f"{operation} requires the complete CellTypist producer evidence schema.")
    if (
        provenance["schema_version"] != 1
        or provenance["operation"] != "celltypist_annotation"
        or provenance["annotation_status"] != "provisional"
    ):
        raise ValueError(f"{operation} producer must be schema-version-1 provisional CellTypist annotation.")
    columns = provenance["columns"]
    matrices = provenance["matrices"]
    if not isinstance(columns, Mapping) or set(columns) != {
        "selected_label",
        "selected_label_probability",
        "individual_label",
        "individual_row_max_probability",
        "majority_label",
        "majority_support",
    }:
        raise ValueError(f"{operation} observation-column schema is invalid.")
    if not isinstance(matrices, Mapping) or set(matrices) != {
        "probability_key",
        "probability_shape",
        "decision_key",
        "decision_shape",
    }:
        raise ValueError(f"{operation} matrix schema is invalid.")
    classes_value = provenance["classes"]
    if hasattr(classes_value, "tolist"):
        classes_value = classes_value.tolist()
    if not isinstance(classes_value, (list, tuple)):
        raise ValueError(f"{operation} model class axis must be an ordered sequence.")
    classes = list(classes_value)
    if (
        not classes
        or len(classes) > 500
        or any(not isinstance(label, str) or not label or label != label.strip() for label in classes)
        or len(set(classes)) != len(classes)
    ):
        raise ValueError(f"{operation} model class axis is invalid.")
    observation_ids = list(adata.obs_names)
    if provenance["observation_axis_fingerprint_sha256"] != _celltypist_axis_fingerprint(observation_ids):
        raise ValueError(f"{operation} observation axis fingerprint does not match producer evidence.")

    selected_label_key = columns["selected_label"]
    confidence_key = columns["selected_label_probability"]
    individual_label_key = columns["individual_label"]
    individual_probability_key = columns["individual_row_max_probability"]
    required_columns = [
        selected_label_key,
        confidence_key,
        individual_label_key,
        individual_probability_key,
    ]
    if any(not isinstance(key, str) or not key or key not in adata.obs for key in required_columns):
        raise ValueError(f"{operation} required CellTypist observation columns are missing.")
    selected_series = adata.obs[selected_label_key]
    individual_series = adata.obs[individual_label_key]
    if not isinstance(selected_series.dtype, pd.CategoricalDtype) or not isinstance(
        individual_series.dtype, pd.CategoricalDtype
    ):
        raise ValueError(f"{operation} CellTypist labels must retain categorical axes.")
    selected_categories = selected_series.cat.categories.tolist()
    expected_selected_categories = classes + (
        ["Heterogeneous"] if "Heterogeneous" in selected_series.astype(str).tolist() else []
    )
    if selected_categories != expected_selected_categories or individual_series.cat.categories.tolist() != classes:
        raise ValueError(f"{operation} label category axes do not match model provenance.")
    selected_labels = selected_series.astype(str).tolist()
    individual_labels = individual_series.astype(str).tolist()
    confidence = pd.to_numeric(adata.obs[confidence_key], errors="raise").to_numpy(dtype=float)
    individual_probability = pd.to_numeric(
        adata.obs[individual_probability_key], errors="raise"
    ).to_numpy(dtype=float)
    if confidence.shape != (len(observation_ids),) or individual_probability.shape != (
        len(observation_ids),
    ):
        raise ValueError(f"{operation} confidence vectors are not observation aligned.")
    expected_missing = np.asarray([label == "Heterogeneous" for label in selected_labels], dtype=bool)
    if not bool(np.array_equal(np.isnan(confidence), expected_missing)):
        raise ValueError(f"{operation} selected confidence missingness does not match Heterogeneous semantics.")
    finite_confidence = confidence[~expected_missing]
    if (
        not bool(np.isfinite(finite_confidence).all())
        or bool(((finite_confidence < 0.0) | (finite_confidence > 1.0)).any())
        or not bool(np.isfinite(individual_probability).all())
        or bool(((individual_probability < 0.0) | (individual_probability > 1.0)).any())
    ):
        raise ValueError(f"{operation} confidence values must lie in [0, 1].")

    probability_key = matrices["probability_key"]
    if not isinstance(probability_key, str) or probability_key not in adata.obsm:
        raise ValueError(f"{operation} stored probability matrix is missing.")
    probabilities = np.asarray(adata.obsm[probability_key])
    expected_shape = (len(observation_ids), len(classes))
    if tuple(matrices["probability_shape"]) != expected_shape or probabilities.shape != expected_shape:
        raise ValueError(f"{operation} probability matrix does not match its stored axes.")
    if (
        not np.issubdtype(probabilities.dtype, np.number)
        or np.issubdtype(probabilities.dtype, np.bool_)
        or np.iscomplexobj(probabilities)
    ):
        raise TypeError(f"{operation} probabilities must be real numeric values.")
    probabilities = np.asarray(probabilities, dtype=float)
    if not bool(np.isfinite(probabilities).all()) or bool(
        ((probabilities < 0.0) | (probabilities > 1.0)).any()
    ):
        raise ValueError(f"{operation} probabilities must be finite in [0, 1].")
    probability_fingerprint = _celltypist_probability_fingerprint(
        probabilities,
        observation_ids=observation_ids,
        classes=classes,
    )
    if provenance["probability_content_fingerprint_sha256"] != probability_fingerprint:
        raise ValueError(f"{operation} probability content fingerprint does not match producer evidence.")
    selected_fingerprint = _celltypist_selected_fingerprint(
        selected_labels,
        confidence,
        observation_ids=observation_ids,
    )
    if provenance["selected_annotation_fingerprint_sha256"] != selected_fingerprint:
        raise ValueError(f"{operation} selected annotation fingerprint does not match producer evidence.")
    class_position = {label: index for index, label in enumerate(classes)}
    individual_positions = np.asarray([class_position[label] for label in individual_labels], dtype=int)
    if not bool(
        np.allclose(
            individual_probability,
            probabilities[np.arange(len(observation_ids)), individual_positions],
            rtol=0.0,
            atol=1e-12,
        )
    ) or not bool(
        np.allclose(individual_probability, probabilities.max(axis=1), rtol=0.0, atol=1e-12)
    ):
        raise ValueError(f"{operation} individual predictions disagree with the stored probability matrix.")
    expected_selected_probability = np.asarray(
        [
            probabilities[index, class_position[label]] if label in class_position else np.nan
            for index, label in enumerate(selected_labels)
        ],
        dtype=float,
    )
    if not bool(
        np.allclose(confidence, expected_selected_probability, rtol=0.0, atol=1e-12, equal_nan=True)
    ):
        raise ValueError(f"{operation} selected confidence disagrees with the stored probability matrix.")

    openbio = adata.uns.get("openbio_singlecell")
    annotations = openbio.get("annotations") if isinstance(openbio, Mapping) else None
    history = openbio.get("analysis_history") if isinstance(openbio, Mapping) else None

    def plain_value(value):
        if isinstance(value, Mapping):
            return {str(key): plain_value(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
            sequence = value.tolist() if hasattr(value, "tolist") else value
            return [plain_value(item) for item in sequence]
        if hasattr(value, "item"):
            return plain_value(value.item())
        return value

    annotation_provenance = annotations.get(selected_label_key) if isinstance(annotations, Mapping) else None
    if not isinstance(annotation_provenance, Mapping) or plain_value(annotation_provenance) != plain_value(
        provenance
    ):
        raise ValueError(f"{operation} selected annotation does not match OpenBio producer provenance.")
    model = provenance["model"]
    model_sha256 = model.get("sha256") if isinstance(model, Mapping) else None
    if not isinstance(history, Mapping) or not any(
        isinstance(entry, Mapping)
        and entry.get("operation") == "celltypist_annotation"
        and isinstance(entry.get("parameters"), Mapping)
        and entry["parameters"].get("metadata_key") == metadata_key
        and entry["parameters"].get("label_column") == selected_label_key
        and entry["parameters"].get("probability_key") == probability_key
        and entry["parameters"].get("model_sha256") == model_sha256
        for entry in history.values()
    ):
        raise ValueError(f"{operation} could not verify the CellTypist producer history.")

    label_counts = {label: selected_labels.count(label) for label in selected_categories}
    warnings = []
    if variant == "confidence_distributions":
        finite_labels = [label for label in selected_categories if label != "Heterogeneous" and label_counts[label]]
        ranked_labels = sorted(
            finite_labels,
            key=lambda label: (-label_counts[label], selected_categories.index(label)),
        )
        plotted_labels = ranked_labels[:limit]
        if not plotted_labels:
            raise ValueError(f"{operation} has no finite selected-label confidence distributions to plot.")
        values_by_label = [
            confidence[np.asarray(selected_labels, dtype=object) == label] for label in plotted_labels
        ]
        figure = Figure(
            figsize=(8.5, max(4.5, min(16.0, 2.5 + 0.38 * len(plotted_labels)))),
            constrained_layout=True,
        )
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        axis.boxplot(
            values_by_label,
            tick_labels=plotted_labels,
            orientation="horizontal",
            showfliers=False,
        )
        axis.set_xlim(0.0, 1.0)
        axis.set_xlabel("Selected-class sigmoid score")
        axis.set_ylabel("Provisional label")
        axis.set_title("CellTypist selected-label confidence distributions")
        axis.grid(axis="x", alpha=0.2)
        groups = []
        plotted_classes = []
        group_means = []
        title = "CellTypist confidence distributions"
        if len(ranked_labels) > len(plotted_labels):
            warnings.append(
                f"Displayed {len(plotted_labels)} of {len(ranked_labels)} observed finite-confidence labels "
                "by descending cell count with model-class order for ties."
            )
    else:
        groupby = mode["groupby"]
        if not isinstance(groupby, str) or not groupby or groupby != groupby.strip():
            raise ValueError(f"{operation} groupby must be a canonical nonempty string.")
        if groupby not in adata.obs:
            raise ValueError(f"{operation} groupby column not found: {groupby!r}.")
        grouping = adata.obs[groupby]
        if isinstance(grouping.dtype, pd.CategoricalDtype):
            if bool(grouping.isna().any()):
                raise ValueError(f"{operation} groupby cannot contain missing values.")
            raw_groups = grouping.cat.categories.tolist()
            codes = grouping.cat.codes.to_numpy(dtype=int)
            if any(not bool((codes == index).any()) for index in range(len(raw_groups))):
                raise ValueError(f"{operation} groupby cannot contain unused categories.")
            groups = [str(value) for value in raw_groups]
        else:
            if bool(grouping.isna().any()):
                raise ValueError(f"{operation} groupby cannot contain missing values.")
            raw_values = grouping.tolist()
            raw_groups = list(dict.fromkeys(raw_values))
            groups = [str(value) for value in raw_groups]
            position = {value: index for index, value in enumerate(raw_groups)}
            codes = np.asarray([position[value] for value in raw_values], dtype=int)
        if any(not label or label != label.strip() for label in groups) or len(set(groups)) != len(groups):
            raise ValueError(f"{operation} group labels must be unique canonical strings after serialization.")
        if not groups or len(groups) > 100:
            raise ValueError(f"{operation} supports between 1 and 100 groups in a probability heatmap.")
        global_means = probabilities.mean(axis=0)
        class_order = sorted(range(len(classes)), key=lambda index: (-global_means[index], index))
        selected_positions = class_order[:limit]
        plotted_classes = [classes[index] for index in selected_positions]
        group_means_array = np.vstack(
            [probabilities[codes == index][:, selected_positions].mean(axis=0) for index in range(len(groups))]
        )
        group_means = group_means_array.tolist()
        figure = Figure(
            figsize=(
                max(7.0, min(20.0, 3.0 + 0.48 * len(plotted_classes))),
                max(4.5, min(16.0, 2.5 + 0.38 * len(groups))),
            ),
            constrained_layout=True,
        )
        FigureCanvasAgg(figure)
        axis = figure.subplots()
        image = axis.imshow(group_means_array, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
        axis.set_xticks(np.arange(len(plotted_classes)), plotted_classes, rotation=45, ha="right")
        axis.set_yticks(np.arange(len(groups)), groups)
        axis.set_xlabel("CellTypist model class")
        axis.set_ylabel(groupby)
        axis.set_title("Mean CellTypist class scores by explicit group")
        colorbar = figure.colorbar(image, ax=axis)
        colorbar.set_label("Mean sigmoid score")
        title = "CellTypist probability heatmap"
        plotted_labels = []
        if len(classes) > len(plotted_classes):
            warnings.append(
                f"Displayed {len(plotted_classes)} of {len(classes)} model classes by descending global mean "
                "sigmoid score with model-class order for ties."
            )

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120, bbox_inches="tight")
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    details = {
        "view": variant,
        "view_parameters": mode,
        "metadata_key": metadata_key,
        "annotation_status": "provisional",
        "model_classes": classes,
        "model_class_count": len(classes),
        "label_counts": label_counts,
        "plotted_labels": plotted_labels,
        "groups": groups,
        "plotted_classes": plotted_classes,
        "group_mean_probabilities": group_means,
        "plotted_cells": len(observation_ids),
        "probability_shape": list(probabilities.shape),
        "observation_axis_fingerprint_sha256": provenance["observation_axis_fingerprint_sha256"],
        "probability_content_fingerprint_sha256": probability_fingerprint,
        "selected_annotation_fingerprint_sha256": selected_fingerprint,
        "model_sha256": model_sha256,
        "selection_policy": (
            "descending selected-label cell count; model class order for ties"
            if variant == "confidence_distributions"
            else "descending global mean sigmoid score; model class order for ties"
        ),
        "render_limits": {"maximum_cells": 1_000_000, "maximum_classes": 500},
        "title": title,
        "warnings": warnings,
    }
    return (png, details) if _return_details else png


def celltypist_diagnostics_plot_code(*, metadata_key: str, view: dict[str, Any]) -> str:
    helpers = "\n\n".join(
        textwrap.dedent(inspect.getsource(helper)).strip()
        for helper in (
            _celltypist_axis_fingerprint,
            _celltypist_probability_fingerprint,
            _celltypist_selected_fingerprint,
        )
    )
    implementation = textwrap.dedent(inspect.getsource(_standalone_celltypist_diagnostics_plot)).strip()
    return f"""from __future__ import annotations

{helpers}


{implementation}


def plot_celltypist_diagnostics(adata):
    return _standalone_celltypist_diagnostics_plot(
        adata,
        metadata_key={metadata_key!r},
        view={view!r},
    )
"""


__all__ = ["_standalone_celltypist_diagnostics_plot", "celltypist_diagnostics_plot_code"]
