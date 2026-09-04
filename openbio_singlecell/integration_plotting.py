from __future__ import annotations

import inspect
from typing import Any

from .scvi_model import (
    SCVI_EPOCH_METRICS,
    SCVI_NATIVE_DIAGNOSTIC_SCHEMA,
    SCVI_TRAINING_DIAGNOSTIC_FILENAME,
    SCVI_TRAINING_DIAGNOSTIC_SCHEMA,
    _file_sha256,
    _fingerprint_json,
    _json_safe,
    read_scvi_training_diagnostics,
    validate_scvi_training_diagnostics,
)

HARMONY_DIAGNOSTIC_SCHEMA = "openbio-singlecell/harmony-diagnostics/v1"


def make_harmony_diagnostic_evidence(
    input_basis: Any,
    adjusted_values: Any,
    observation_names: Any,
    objective_harmony: Any,
    objective_kmeans: Any,
    kmeans_rounds: Any,
    *,
    basis: str,
    adjusted_basis: str,
    technical_batch_keys: list[str],
    technical_batch_values: Any,
    theta_mode: str,
    theta: float | None,
    harmonypy_version: str,
    random_seed: int,
    max_iter_harmony: int,
    max_iter_kmeans: int,
) -> dict[str, Any]:
    import hashlib
    import json

    import numpy as np

    def array_fingerprint(value):
        array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest = hashlib.sha256()
        digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
        digest.update(memoryview(array).cast("B"))
        return digest.hexdigest()

    def json_fingerprint(value):
        encoded = json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def objective_sequence(value, name):
        objective = np.asarray(value)
        if objective.ndim != 1 or objective.size == 0:
            raise RuntimeError(f"harmonypy {name} must be a non-empty one-dimensional sequence.")
        if objective.dtype.kind not in "iuf" or np.issubdtype(objective.dtype, np.bool_) or np.iscomplexobj(objective):
            raise RuntimeError(f"harmonypy {name} must contain real numeric values.")
        objective = np.asarray(objective, dtype=float)
        if not bool(np.isfinite(objective).all()):
            raise RuntimeError(f"harmonypy {name} must contain only finite values.")
        return objective

    harmony_objective = objective_sequence(objective_harmony, "objective_harmony")
    kmeans_objective = objective_sequence(objective_kmeans, "objective_kmeans")
    rounds_array = np.asarray(kmeans_rounds)
    if rounds_array.ndim != 1 or rounds_array.size == 0:
        raise RuntimeError("harmonypy kmeans_rounds must be a non-empty one-dimensional sequence.")
    if rounds_array.dtype.kind not in "iu" or np.issubdtype(rounds_array.dtype, np.bool_):
        raise RuntimeError("harmonypy kmeans_rounds must contain integers.")
    rounds = [int(value) for value in rounds_array]
    if not 1 <= len(rounds) <= max_iter_harmony:
        raise RuntimeError("harmonypy kmeans_rounds violates the configured outer-iteration bound.")
    if any(value < 1 or value > max_iter_kmeans for value in rounds):
        raise RuntimeError("harmonypy kmeans_rounds violates the configured clustering-iteration bound.")
    if harmony_objective.size != len(rounds) + 1:
        raise RuntimeError("harmonypy objective_harmony must contain initialization plus one value per outer round.")
    if kmeans_objective.size != 1 + sum(rounds):
        raise RuntimeError("harmonypy objective_kmeans must contain initialization plus every clustering update.")
    boundaries = [0, *np.cumsum(rounds).astype(int).tolist()]
    if not bool(np.array_equal(harmony_objective, kmeans_objective[boundaries])):
        raise RuntimeError("harmonypy objective sequences disagree at their shared outer-round boundaries.")
    if not isinstance(basis, str) or not basis or not isinstance(adjusted_basis, str) or not adjusted_basis:
        raise RuntimeError("Harmony diagnostic basis identifiers must be non-empty strings.")
    input_array = np.asarray(input_basis)
    adjusted_array = np.asarray(adjusted_values)
    if input_array.shape != adjusted_array.shape or adjusted_array.ndim != 2:
        raise RuntimeError("Harmony diagnostic coordinate matrices must have the same two-dimensional shape.")
    names = list(observation_names)
    if len(names) != adjusted_array.shape[0] or any(not isinstance(name, str) or not name for name in names):
        raise RuntimeError("Harmony diagnostic observation identity is invalid.")
    if not isinstance(technical_batch_values, dict) or list(technical_batch_values) != technical_batch_keys:
        raise RuntimeError("Harmony diagnostic Technical-batch provenance is invalid.")
    technical_batches = []
    for index, key in enumerate(technical_batch_keys):
        labels = list(technical_batch_values[key])
        if len(labels) != len(names) or any(not isinstance(label, str) or not label for label in labels):
            raise RuntimeError("Harmony diagnostic Technical-batch labels are invalid.")
        technical_batches.append(
            {"index": index, "key": key, "labels_fingerprint_sha256": json_fingerprint(labels)}
        )

    payload = {
        "adjusted_basis": adjusted_basis,
        "basis": basis,
        "technical_batches": technical_batches,
        "theta_mode": theta_mode,
        "theta": theta,
        "harmonypy_version": harmonypy_version,
        "random_seed": random_seed,
        "max_iter_harmony": max_iter_harmony,
        "max_iter_kmeans": max_iter_kmeans,
        "harmony_objective": [
            {"iteration": index, "value": float(value)}
            for index, value in enumerate(harmony_objective)
        ],
        "kmeans_objective": [
            {"update": index, "value": float(value)}
            for index, value in enumerate(kmeans_objective)
        ],
        "outer_rounds": [
            {"round": index + 1, "kmeans_updates": value} for index, value in enumerate(rounds)
        ],
        "coordinate_shape": {
            "observations": int(adjusted_array.shape[0]),
            "components": int(adjusted_array.shape[1]),
        },
        "input_basis_fingerprint_sha256": array_fingerprint(input_array),
        "adjusted_basis_fingerprint_sha256": array_fingerprint(adjusted_array),
        "observation_axis_fingerprint_sha256": json_fingerprint(names),
    }
    return {
        "schema_version": HARMONY_DIAGNOSTIC_SCHEMA,
        "payload": payload,
        "payload_sha256": json_fingerprint(payload),
    }


def _standalone_harmony_convergence_plot(
    adata,
    *,
    adjusted_basis="X_pca_harmony",
    max_image_pixels=40_000_000,
    _return_diagnostics=False,
):
    import hashlib
    import io
    import json
    from collections.abc import Mapping

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "Harmony Convergence Plot"
    if not isinstance(adjusted_basis, str) or not adjusted_basis or adjusted_basis != adjusted_basis.strip():
        raise ValueError(f"{operation} adjusted_basis must be a canonical non-empty string.")
    if isinstance(max_image_pixels, (bool, np.bool_)) or not isinstance(max_image_pixels, (int, np.integer)):
        raise TypeError(f"{operation} max_image_pixels must be an integer.")
    max_image_pixels = int(max_image_pixels)
    if max_image_pixels < 1:
        raise ValueError(f"{operation} max_image_pixels must be positive.")

    metadata = adata.uns.get("harmony")
    if not isinstance(metadata, Mapping) or adjusted_basis not in metadata:
        raise ValueError(f"{operation} requires retained diagnostics for adjusted basis {adjusted_basis!r}.")
    evidence = metadata[adjusted_basis]
    if not isinstance(evidence, Mapping) or set(evidence) != {"schema_version", "payload", "payload_sha256"}:
        raise ValueError(f"{operation} Harmony diagnostic evidence schema is invalid.")
    if evidence["schema_version"] != HARMONY_DIAGNOSTIC_SCHEMA:
        raise ValueError(f"{operation} Harmony diagnostic evidence schema version is unsupported.")
    payload = evidence["payload"]
    expected_payload_fields = {
        "adjusted_basis",
        "basis",
        "technical_batches",
        "theta_mode",
        "theta",
        "harmonypy_version",
        "random_seed",
        "max_iter_harmony",
        "max_iter_kmeans",
        "harmony_objective",
        "kmeans_objective",
        "outer_rounds",
        "coordinate_shape",
        "input_basis_fingerprint_sha256",
        "adjusted_basis_fingerprint_sha256",
        "observation_axis_fingerprint_sha256",
    }
    if not isinstance(payload, Mapping) or set(payload) != expected_payload_fields:
        raise ValueError(f"{operation} Harmony diagnostic payload schema is invalid.")

    def plain_json(value):
        if value is None or isinstance(value, (str, bool, int, float)):
            return value.item() if hasattr(value, "item") else value
        if hasattr(value, "tolist"):
            return plain_json(value.tolist())
        if isinstance(value, Mapping):
            return {str(key): plain_json(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [plain_json(item) for item in value]
        raise TypeError(f"{operation} Harmony diagnostic payload contains a non-JSON value.")

    normalized = plain_json(payload)
    encoded = json.dumps(
        normalized,
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    payload_sha256 = hashlib.sha256(encoded).hexdigest()
    if evidence["payload_sha256"] != payload_sha256:
        raise ValueError(f"{operation} Harmony diagnostic payload fingerprint does not match its contents.")
    if normalized["adjusted_basis"] != adjusted_basis:
        raise ValueError(f"{operation} adjusted-basis identity does not match the requested representation.")
    basis = normalized["basis"]
    if not isinstance(basis, str) or not basis or basis != basis.strip() or basis not in adata.obsm:
        raise ValueError(f"{operation} source-basis provenance is unavailable.")
    key_records = normalized["technical_batches"]
    if not isinstance(key_records, list) or not key_records:
        raise ValueError(f"{operation} Technical-batch provenance is invalid.")
    technical_batch_keys = []
    for expected_index, record in enumerate(key_records):
        if (
            not isinstance(record, dict)
            or set(record) != {"index", "key", "labels_fingerprint_sha256"}
            or record["index"] != expected_index
        ):
            raise ValueError(f"{operation} Technical-batch provenance order is invalid.")
        key = record["key"]
        if not isinstance(key, str) or not key or key != key.strip() or key in technical_batch_keys:
            raise ValueError(f"{operation} Technical-batch provenance keys are invalid.")
        if key not in adata.obs:
            raise ValueError(f"{operation} Technical-batch provenance column is unavailable: {key!r}.")
        labels = adata.obs[key].astype(str).tolist()
        labels_encoded = json.dumps(
            labels,
            allow_nan=False,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        if hashlib.sha256(labels_encoded).hexdigest() != record["labels_fingerprint_sha256"]:
            raise ValueError(f"{operation} Technical-batch labels no longer match their diagnostic fingerprint.")
        technical_batch_keys.append(key)
    theta_mode = normalized["theta_mode"]
    theta = normalized["theta"]
    if theta_mode not in {"automatic", "custom"} or (
        theta_mode == "automatic" and theta is not None
    ) or (
        theta_mode == "custom"
        and (isinstance(theta, bool) or not isinstance(theta, (int, float)) or not np.isfinite(theta) or theta < 0)
    ):
        raise ValueError(f"{operation} retained theta provenance is invalid.")
    harmonypy_version = normalized["harmonypy_version"]
    if not isinstance(harmonypy_version, str) or not harmonypy_version.startswith("2.0."):
        raise ValueError(f"{operation} retained harmonypy version provenance is invalid.")
    random_seed = normalized["random_seed"]
    if isinstance(random_seed, bool) or not isinstance(random_seed, int) or not 0 <= random_seed <= 2**31 - 1:
        raise ValueError(f"{operation} retained random-seed provenance is invalid.")
    if adjusted_basis not in adata.obsm:
        raise ValueError(f"{operation} adjusted representation is unavailable: {adjusted_basis!r}.")

    def array_fingerprint(value):
        array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest = hashlib.sha256()
        digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
        digest.update(memoryview(array).cast("B"))
        return digest.hexdigest()

    input_values = np.asarray(adata.obsm[basis])
    adjusted_values = np.asarray(adata.obsm[adjusted_basis])
    if (
        input_values.ndim != 2
        or input_values.dtype.kind not in "iuf"
        or adjusted_values.dtype.kind not in "iuf"
        or not bool(np.isfinite(input_values).all())
        or not bool(np.isfinite(adjusted_values).all())
    ):
        raise ValueError(f"{operation} coordinate provenance requires finite real matrices.")
    shape = {"observations": int(adjusted_values.shape[0]), "components": int(adjusted_values.shape[1])}
    if input_values.shape != adjusted_values.shape or shape != normalized["coordinate_shape"]:
        raise ValueError(f"{operation} stored coordinate shape does not match its diagnostic provenance.")
    fingerprints = (
        normalized["input_basis_fingerprint_sha256"],
        normalized["adjusted_basis_fingerprint_sha256"],
        normalized["observation_axis_fingerprint_sha256"],
    )
    if any(
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
        for value in fingerprints
    ):
        raise ValueError(f"{operation} retained SHA-256 fingerprints are invalid.")
    if array_fingerprint(input_values) != normalized["input_basis_fingerprint_sha256"]:
        raise ValueError(f"{operation} source basis no longer matches its diagnostic fingerprint.")
    if array_fingerprint(adjusted_values) != normalized["adjusted_basis_fingerprint_sha256"]:
        raise ValueError(f"{operation} adjusted basis no longer matches its diagnostic fingerprint.")
    names_encoded = json.dumps(
        [str(value) for value in adata.obs_names],
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    if hashlib.sha256(names_encoded).hexdigest() != normalized["observation_axis_fingerprint_sha256"]:
        raise ValueError(f"{operation} observation axis no longer matches its diagnostic fingerprint.")

    def objective_records(name, index_name):
        records = normalized[name]
        if not isinstance(records, list) or not records:
            raise ValueError(f"{operation} {name} must be a non-empty record sequence.")
        values = []
        for expected_index, record in enumerate(records):
            if not isinstance(record, dict) or set(record) != {index_name, "value"}:
                raise ValueError(f"{operation} {name} record schema is invalid.")
            if record[index_name] != expected_index:
                raise ValueError(f"{operation} {name} indices are not canonical.")
            value = record["value"]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError(f"{operation} {name} must contain finite real values.")
            values.append(float(value))
        return np.asarray(values, dtype=float)

    harmony_objective = objective_records("harmony_objective", "iteration")
    kmeans_objective = objective_records("kmeans_objective", "update")
    round_records = normalized["outer_rounds"]
    if not isinstance(round_records, list) or not round_records:
        raise ValueError(f"{operation} outer_rounds must be a non-empty record sequence.")
    rounds = []
    for expected_round, record in enumerate(round_records, start=1):
        if not isinstance(record, dict) or set(record) != {"round", "kmeans_updates"}:
            raise ValueError(f"{operation} outer_rounds record schema is invalid.")
        if record["round"] != expected_round:
            raise ValueError(f"{operation} outer_rounds identities are not canonical.")
        value = record["kmeans_updates"]
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"{operation} outer_rounds update counts must be integers.")
        rounds.append(value)
    max_harmony = normalized["max_iter_harmony"]
    max_kmeans = normalized["max_iter_kmeans"]
    if (
        isinstance(max_harmony, bool)
        or not isinstance(max_harmony, int)
        or max_harmony < 1
        or isinstance(max_kmeans, bool)
        or not isinstance(max_kmeans, int)
        or max_kmeans < 1
    ):
        raise ValueError(f"{operation} retained iteration bounds are invalid.")
    if not 1 <= len(rounds) <= max_harmony or any(value < 1 or value > max_kmeans for value in rounds):
        raise ValueError(f"{operation} retained iteration counts violate their configured bounds.")
    if harmony_objective.size != len(rounds) + 1 or kmeans_objective.size != 1 + sum(rounds):
        raise ValueError(f"{operation} retained objective sequences do not match the recorded round identities.")
    outer_boundaries = [0, *np.cumsum(rounds).astype(int).tolist()]
    if not bool(np.array_equal(harmony_objective, kmeans_objective[outer_boundaries])):
        raise ValueError(f"{operation} retained objective sequences disagree at outer-round boundaries.")
    harmony_iterations = list(range(int(harmony_objective.size)))
    kmeans_updates = list(range(int(kmeans_objective.size)))
    # ponytail: one static curve supports 100k stored points; use pagination if real training ever exceeds this.
    if harmony_objective.size + kmeans_objective.size > 100_000:
        raise ValueError(f"{operation} supports at most 100,000 objective points in one static PNG.")

    width, height, dpi = 8.0, 6.0, 120
    rendered_pixels = int(width * dpi) * int(height * dpi)
    if rendered_pixels > max_image_pixels:
        raise ValueError(
            f"{operation} requires {rendered_pixels:,} pixels, exceeding max_image_pixels={max_image_pixels:,}."
        )
    figure = Figure(figsize=(width, height))
    FigureCanvasAgg(figure)
    harmony_axis, kmeans_axis = figure.subplots(2, 1)
    harmony_axis.plot(harmony_iterations, harmony_objective, color="#4c78a8", marker="o", linewidth=2.0)
    harmony_axis.set_title("Harmony outer-loop objective")
    harmony_axis.set_xlabel("Outer round (0 = initialization)")
    harmony_axis.set_ylabel("Objective")
    harmony_axis.grid(alpha=0.25)
    kmeans_axis.plot(kmeans_updates, kmeans_objective, color="#f58518", linewidth=1.8)
    interior_boundaries = outer_boundaries[1:-1]
    for boundary_index, boundary in enumerate(interior_boundaries):
        kmeans_axis.axvline(
            boundary,
            color="#777777",
            linestyle="--",
            linewidth=0.8,
            alpha=0.7,
            label="Outer-round boundary" if boundary_index == 0 else None,
        )
    kmeans_axis.set_title("Harmony clustering objective")
    kmeans_axis.set_xlabel("Clustering update (0 = initialization)")
    kmeans_axis.set_ylabel("Objective")
    kmeans_axis.grid(alpha=0.25)
    if interior_boundaries:
        kmeans_axis.legend(frameon=False, loc="best")
    figure.tight_layout()
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi)
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")

    details = {
        "adjusted_basis": adjusted_basis,
        "basis": basis,
        "technical_batch_keys": technical_batch_keys,
        "harmonypy_version": harmonypy_version,
        "harmony_objective_points": int(harmony_objective.size),
        "kmeans_objective_points": int(kmeans_objective.size),
        "outer_rounds": len(rounds),
        "kmeans_rounds": rounds,
        "outer_round_boundaries": outer_boundaries,
        "boundary_line_legend": "Outer-round boundary" if interior_boundaries else None,
        "objective_harmony_initial": float(harmony_objective[0]),
        "objective_harmony_final": float(harmony_objective[-1]),
        "objective_kmeans_initial": float(kmeans_objective[0]),
        "objective_kmeans_final": float(kmeans_objective[-1]),
        "payload_sha256": payload_sha256,
        "input_cells": int(adata.n_obs),
        "input_genes": int(adata.n_vars),
        "rendered_pixels": rendered_pixels,
        "warnings": [],
    }
    return (png, details) if _return_diagnostics else png


def render_harmony_convergence_plot(
    adata: Any,
    *,
    adjusted_basis: str = "X_pca_harmony",
    max_image_pixels: int = 40_000_000,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_harmony_convergence_plot(
        adata,
        adjusted_basis=adjusted_basis,
        max_image_pixels=max_image_pixels,
        _return_diagnostics=True,
    )


def harmony_convergence_plot_code(*, adjusted_basis: str, max_image_pixels: int) -> str:
    implementation = inspect.getsource(_standalone_harmony_convergence_plot)
    return (
        "from __future__ import annotations\n\n"
        f"HARMONY_DIAGNOSTIC_SCHEMA = {HARMONY_DIAGNOSTIC_SCHEMA!r}\n\n\n"
        f"{implementation}\n\n\n"
        "def plot_harmony_convergence(adata):\n"
        "    return _standalone_harmony_convergence_plot(\n"
        f"        adata, adjusted_basis={adjusted_basis!r}, max_image_pixels={max_image_pixels!r}\n"
        "    )\n"
    )


def _standalone_scvi_training_plot(
    diagnostics,
    *,
    max_image_pixels=40_000_000,
    _return_diagnostics=False,
):
    import io

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    operation = "scVI Training Plot"
    validated = validate_scvi_training_diagnostics(diagnostics)
    if isinstance(max_image_pixels, (bool, np.bool_)) or not isinstance(max_image_pixels, (int, np.integer)):
        raise TypeError(f"{operation} max_image_pixels must be an integer.")
    max_image_pixels = int(max_image_pixels)
    if max_image_pixels < 1:
        raise ValueError(f"{operation} max_image_pixels must be positive.")
    metric_families = (
        ("elbo", "ELBO"),
        ("reconstruction_loss", "Reconstruction loss"),
        ("kl_local", "Local KL divergence"),
        ("kl_global", "Global KL divergence"),
    )
    metrics = validated["metrics"]
    selected_families = [
        (family, label)
        for family, label in metric_families
        if any(f"{family}_{split}" in metrics for split in ("train", "validation"))
    ]
    panels = [family for family, _label in selected_families]
    total_points = sum(len(payload["records"]) for payload in metrics.values())
    # ponytail: one static figure supports 100k epoch-metric points; use pagination if training ever exceeds this.
    if total_points > 100_000:
        raise ValueError(f"{operation} supports at most 100,000 epoch-metric points in one static PNG.")
    width, height, dpi = 8.0, 3.0 * len(panels), 120
    rendered_pixels = int(width * dpi) * int(height * dpi)
    if rendered_pixels > max_image_pixels:
        raise ValueError(
            f"{operation} requires {rendered_pixels:,} pixels, exceeding max_image_pixels={max_image_pixels:,}."
        )
    figure = Figure(figsize=(width, height))
    FigureCanvasAgg(figure)
    axes = np.atleast_1d(figure.subplots(len(panels), 1)).tolist()
    metric_summaries = {}
    plotted_metrics = []
    for axis, (family, label) in zip(axes, selected_families, strict=True):
        for split, color in (("train", "#4c78a8"), ("validation", "#f58518")):
            metric = f"{family}_{split}"
            if metric not in metrics:
                continue
            records = metrics[metric]["records"]
            epochs = [record["epoch"] for record in records]
            values = [record["value"] for record in records]
            axis.plot(epochs, values, color=color, linewidth=2.0, marker="o", markersize=3.0, label=split.title())
            plotted_metrics.append(metric)
            metric_summaries[metric] = {
                "points": len(records),
                "first_epoch": epochs[0],
                "last_epoch": epochs[-1],
                "first": values[0],
                "last": values[-1],
                "best": min(values),
            }
        axis.set_title(label)
        axis.set_xlabel("Epoch (native zero-based identity)")
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    figure.tight_layout()
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi)
    png = buffer.getvalue()
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError(f"{operation} renderer did not produce a valid PNG payload.")
    warnings = []
    if validated["validation_cells"] == 0:
        warnings.append("The fitted model retained no validation holdout or validation metric sequence.")
    details = {
        "actual_epochs": validated["actual_epochs"],
        "plotted_panels": panels,
        "plotted_metrics": plotted_metrics,
        "metric_summaries": metric_summaries,
        "train_cells": validated["train_cells"],
        "validation_cells": validated["validation_cells"],
        "test_cells": validated["test_cells"],
        "training_observations": validated["training_observations"],
        "fitted_features": validated["fitted_features"],
        "device": validated["device"],
        "total_plotted_points": total_points,
        "rendered_pixels": rendered_pixels,
        "warnings": warnings,
    }
    return (png, details) if _return_diagnostics else png


def render_scvi_training_plot(
    diagnostics: Any,
    *,
    max_image_pixels: int = 40_000_000,
) -> tuple[bytes, dict[str, Any]]:
    return _standalone_scvi_training_plot(
        diagnostics,
        max_image_pixels=max_image_pixels,
        _return_diagnostics=True,
    )


def scvi_training_plot_code(*, max_image_pixels: int) -> str:
    json_safe = inspect.getsource(_json_safe)
    fingerprint = inspect.getsource(_fingerprint_json)
    file_fingerprint = inspect.getsource(_file_sha256)
    validator = inspect.getsource(validate_scvi_training_diagnostics)
    reader = inspect.getsource(read_scvi_training_diagnostics)
    implementation = inspect.getsource(_standalone_scvi_training_plot)
    return (
        "from __future__ import annotations\n\n"
        "import copy\n"
        "import hashlib\n"
        "import json\n"
        "import math\n"
        "from collections.abc import Mapping, Sequence\n\n"
        f"SCVI_NATIVE_DIAGNOSTIC_SCHEMA = {SCVI_NATIVE_DIAGNOSTIC_SCHEMA!r}\n"
        f"SCVI_TRAINING_DIAGNOSTIC_FILENAME = {SCVI_TRAINING_DIAGNOSTIC_FILENAME!r}\n"
        f"SCVI_TRAINING_DIAGNOSTIC_SCHEMA = {SCVI_TRAINING_DIAGNOSTIC_SCHEMA!r}\n"
        f"SCVI_EPOCH_METRICS = {SCVI_EPOCH_METRICS!r}\n\n\n"
        f"{json_safe}\n\n\n"
        f"{fingerprint}\n\n\n"
        f"{file_fingerprint}\n\n\n"
        f"{validator}\n\n\n"
        f"{reader}\n\n\n"
        f"{implementation}\n\n\n"
        "def plot_scvi_training(native_model_root):\n"
        "    native_diagnostics = read_scvi_training_diagnostics(native_model_root)\n"
        "    return _standalone_scvi_training_plot(\n"
        "        native_diagnostics[\"training_diagnostics\"],\n"
        f"        max_image_pixels={max_image_pixels!r},\n"
        "    )\n"
    )


__all__ = [
    "HARMONY_DIAGNOSTIC_SCHEMA",
    "harmony_convergence_plot_code",
    "make_harmony_diagnostic_evidence",
    "render_harmony_convergence_plot",
    "render_scvi_training_plot",
    "scvi_training_plot_code",
]
