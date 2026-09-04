from __future__ import annotations

import copy
import hashlib
import inspect
import json
import math
import random
import sys
import threading
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


SCVI_MODEL_ARTIFACT_SCHEMA = "openbio-singlecell/scvi-model/v4"
SCVI_TRAINING_DIAGNOSTIC_SCHEMA = "openbio-singlecell/scvi-training-diagnostics/v1"
SCVI_NATIVE_DIAGNOSTIC_SCHEMA = "openbio-singlecell/scvi-native-diagnostics/v1"
SCVI_TRAINING_DIAGNOSTIC_FILENAME = "openbio-training-diagnostics.json"
SCVI_EPOCH_METRICS = (
    "elbo_train",
    "elbo_validation",
    "reconstruction_loss_train",
    "reconstruction_loss_validation",
    "kl_local_train",
    "kl_local_validation",
    "kl_global_train",
    "kl_global_validation",
)
SCVI_GLOBAL_RNG_LOCK = threading.RLock()


@dataclass(frozen=True, slots=True)
class SCVIDifferentialRun:
    """One raw scvi-tools contrast plus evidence needed for validation/reporting."""

    table: DataFrame
    comparison_evidence: Mapping[str, Any]
    backend_parameters: Mapping[str, Any]


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("scVI model evidence cannot contain NaN or infinity.")
        return value
    if hasattr(value, "item"):
        try:
            return _json_safe(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_json_safe(item) for item in value]
    return str(value)


def _fingerprint_json(value: Any) -> str:
    payload = json.dumps(_json_safe(value), allow_nan=False, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_scvi_training_diagnostics(diagnostics: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(diagnostics, Mapping):
        raise TypeError("scVI training diagnostics must be a mapping.")
    expected = {
        "schema_version",
        "actual_epochs",
        "metrics",
        "train_cells",
        "validation_cells",
        "test_cells",
        "training_observations",
        "fitted_features",
        "device",
    }
    if set(diagnostics) != expected:
        raise ValueError(f"scVI training diagnostics fields must be exactly {sorted(expected)!r}.")
    if diagnostics["schema_version"] != SCVI_TRAINING_DIAGNOSTIC_SCHEMA:
        raise ValueError("scVI training diagnostic schema version is unsupported.")
    integer_fields = (
        "actual_epochs",
        "train_cells",
        "validation_cells",
        "test_cells",
        "training_observations",
        "fitted_features",
    )
    integers: dict[str, int] = {}
    for field in integer_fields:
        value = diagnostics[field]
        if hasattr(value, "item"):
            value = value.item()
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"scVI training diagnostic {field} must be an integer.")
        if value < (
            1
            if field in {"actual_epochs", "train_cells", "training_observations", "fitted_features"}
            else 0
        ):
            raise ValueError(f"scVI training diagnostic {field} is outside its valid range.")
        integers[field] = value
    if integers["train_cells"] + integers["validation_cells"] + integers["test_cells"] != integers[
        "training_observations"
    ]:
        raise ValueError("scVI training split counts must cover every training observation exactly once.")
    device = diagnostics["device"]
    if not isinstance(device, str) or not device or device != device.strip():
        raise ValueError("scVI training diagnostic device must be a canonical non-empty string.")
    metrics = diagnostics["metrics"]
    if not isinstance(metrics, Mapping):
        raise TypeError("scVI training diagnostic metrics must be a mapping.")
    required = {"elbo_train", "reconstruction_loss_train", "kl_local_train"}
    allowed = {*required, "kl_global_train"}
    if integers["validation_cells"]:
        required.update({"elbo_validation", "reconstruction_loss_validation", "kl_local_validation"})
        allowed.update({*required, "kl_global_validation"})
    missing = sorted(required.difference(metrics))
    unexpected = sorted(set(metrics).difference(allowed))
    if missing or unexpected:
        raise ValueError(
            "scVI training diagnostic metrics must contain every required sequence and only supported optional "
            f"sequences; missing={missing!r}, unexpected={unexpected!r}."
        )

    normalized_metrics: dict[str, dict[str, list[dict[str, int | float]]]] = {}
    epoch_axes: dict[str, list[int]] = {}
    for metric in SCVI_EPOCH_METRICS:
        if metric not in metrics:
            continue
        metric_payload = metrics[metric]
        if not isinstance(metric_payload, Mapping) or set(metric_payload) != {"records"}:
            raise ValueError(f"scVI training metric {metric!r} must contain exactly one records field.")
        records = metric_payload["records"]
        if not isinstance(records, Sequence) or isinstance(records, (str, bytes, bytearray)) or not records:
            raise ValueError(f"scVI training metric {metric!r} records must be a non-empty sequence.")
        normalized_records: list[dict[str, int | float]] = []
        epochs: list[int] = []
        for record in records:
            if not isinstance(record, Mapping) or set(record) != {"epoch", "value"}:
                raise ValueError(f"scVI training metric {metric!r} record schema is invalid.")
            epoch = record["epoch"]
            value = record["value"]
            if hasattr(epoch, "item"):
                epoch = epoch.item()
            if hasattr(value, "item"):
                value = value.item()
            if isinstance(epoch, bool) or not isinstance(epoch, int) or epoch < 0:
                raise ValueError(f"scVI training metric {metric!r} epochs must be non-negative integers.")
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ValueError(f"scVI training metric {metric!r} values must be finite real numbers.")
            epochs.append(epoch)
            normalized_records.append({"epoch": epoch, "value": float(value)})
        if epochs != list(range(len(epochs))):
            raise ValueError(f"scVI training metric {metric!r} must retain every zero-based epoch exactly once.")
        epoch_axes[metric] = epochs
        normalized_metrics[metric] = {"records": normalized_records}
    train_axes = [
        epoch_axes[metric]
        for metric in SCVI_EPOCH_METRICS
        if metric.endswith("_train") and metric in epoch_axes
    ]
    if any(axis != train_axes[0] for axis in train_axes[1:]):
        raise ValueError("scVI training metric epoch axes must match.")
    validation_axes = [
        epoch_axes[metric]
        for metric in SCVI_EPOCH_METRICS
        if metric.endswith("_validation") and metric in epoch_axes
    ]
    if validation_axes and any(axis != train_axes[0] for axis in validation_axes):
        raise ValueError("scVI validation metric epoch axes must match the training epoch axis.")
    if integers["actual_epochs"] != len(train_axes[0]):
        raise ValueError("scVI actual_epochs must equal the complete retained training epoch count.")
    return {
        "schema_version": SCVI_TRAINING_DIAGNOSTIC_SCHEMA,
        "actual_epochs": integers["actual_epochs"],
        "metrics": normalized_metrics,
        "train_cells": integers["train_cells"],
        "validation_cells": integers["validation_cells"],
        "test_cells": integers["test_cells"],
        "training_observations": integers["training_observations"],
        "fitted_features": integers["fitted_features"],
        "device": device,
    }


def build_scvi_training_diagnostics(model: Any) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    history = getattr(model, "history", None)
    if not isinstance(history, Mapping):
        raise RuntimeError("scVI did not expose its required training-history mapping.")
    validation_cells = len(model.validation_indices)
    required = {"elbo_train", "reconstruction_loss_train", "kl_local_train"}
    if validation_cells:
        required.update({"elbo_validation", "reconstruction_loss_validation", "kl_local_validation"})
    elif any(metric.endswith("_validation") and metric in history for metric in SCVI_EPOCH_METRICS):
        raise RuntimeError("scVI training history contains validation metrics without validation observations.")
    missing = sorted(required.difference(history))
    if missing:
        raise RuntimeError(f"scVI training history is missing supported epoch metrics: {missing}.")
    retained = required | {metric for metric in ("kl_global_train", "kl_global_validation") if metric in history}
    metrics = {}
    for metric in SCVI_EPOCH_METRICS:
        if metric not in retained:
            continue
        frame = history[metric]
        if not isinstance(frame, pd.DataFrame) or frame.columns.tolist() != [metric]:
            raise RuntimeError(f"scVI training history {metric!r} must be a one-column pandas DataFrame.")
        if frame.index.name != "epoch":
            raise RuntimeError(f"scVI training history {metric!r} must use its native epoch index.")
        records = []
        for raw_epoch, raw_value in zip(frame.index.tolist(), frame[metric].tolist(), strict=True):
            epoch = raw_epoch.item() if hasattr(raw_epoch, "item") else raw_epoch
            value = raw_value.item() if hasattr(raw_value, "item") else raw_value
            if isinstance(epoch, bool) or not isinstance(epoch, int):
                raise RuntimeError(f"scVI training history {metric!r} contains a non-integer epoch identity.")
            if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
                raise RuntimeError(f"scVI training history {metric!r} contains a non-numeric value.")
            records.append({"epoch": epoch, "value": float(value)})
        metrics[metric] = {"records": records}
    train_epochs = [record["epoch"] for record in metrics["elbo_train"]["records"]]
    diagnostics = {
        "schema_version": SCVI_TRAINING_DIAGNOSTIC_SCHEMA,
        "actual_epochs": len(train_epochs),
        "metrics": metrics,
        "train_cells": len(model.train_indices),
        "validation_cells": validation_cells,
        "test_cells": len(model.test_indices),
        "training_observations": int(model.adata.n_obs),
        "fitted_features": int(model.adata.n_vars),
        "device": str(model.device),
    }
    try:
        return validate_scvi_training_diagnostics(diagnostics)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"scVI exposed invalid training diagnostics: {error}") from error


def _file_sha256(path: Any) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def write_scvi_training_diagnostics(
    root: Any,
    diagnostics: Mapping[str, Any],
    *,
    training_parameters: Mapping[str, Any],
    software_versions: Mapping[str, Any],
) -> dict[str, Any]:
    from pathlib import Path

    root = Path(root).resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError("scVI native artifact root must be a directory.")
    model_path = root / "model.pt"
    if not model_path.is_file() or model_path.is_symlink():
        raise ValueError("scVI native artifact requires a regular model.pt file before diagnostics are written.")
    sidecar = root / SCVI_TRAINING_DIAGNOSTIC_FILENAME
    if sidecar.exists() or sidecar.is_symlink():
        raise FileExistsError(f"scVI native diagnostics already exist: {sidecar.name}.")
    if not isinstance(training_parameters, Mapping) or not isinstance(software_versions, Mapping):
        raise TypeError("scVI native diagnostic provenance must use mappings.")
    normalized_diagnostics = validate_scvi_training_diagnostics(diagnostics)
    normalized_parameters = _json_safe(dict(training_parameters))
    normalized_versions = _json_safe(dict(software_versions))
    if not isinstance(normalized_parameters, dict) or not isinstance(normalized_versions, dict):
        raise TypeError("scVI native diagnostic provenance must remain JSON mappings.")
    if not all(isinstance(key, str) and isinstance(value, str) for key, value in normalized_versions.items()):
        raise ValueError("scVI native diagnostic software versions must map package names to version strings.")
    payload = {
        "training_diagnostics": normalized_diagnostics,
        "training_parameters": normalized_parameters,
        "software_versions": normalized_versions,
        "model_identity": {
            "path": "model.pt",
            "size": model_path.stat().st_size,
            "sha256": _file_sha256(model_path),
        },
    }
    envelope = {
        "schema_version": SCVI_NATIVE_DIAGNOSTIC_SCHEMA,
        "payload": payload,
        "payload_sha256": _fingerprint_json(payload),
    }
    sidecar.write_text(
        json.dumps(envelope, allow_nan=False, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return copy.deepcopy(payload)


def read_scvi_training_diagnostics(root: Any) -> dict[str, Any]:
    from pathlib import Path

    root = Path(root).resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError("scVI native artifact root must be a directory.")
    sidecar = root / SCVI_TRAINING_DIAGNOSTIC_FILENAME
    if not sidecar.is_file() or sidecar.is_symlink():
        raise ValueError("scVI native artifact is missing its regular OpenBio training-diagnostic sidecar.")

    def reject_constant(value: str) -> None:
        raise ValueError(f"scVI native diagnostics contain invalid JSON constant {value!r}.")

    try:
        envelope = json.loads(sidecar.read_text(encoding="utf-8"), parse_constant=reject_constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("scVI native training diagnostics are unreadable.") from error
    if not isinstance(envelope, Mapping) or set(envelope) != {"schema_version", "payload", "payload_sha256"}:
        raise ValueError("scVI native training diagnostic envelope schema is invalid.")
    if envelope["schema_version"] != SCVI_NATIVE_DIAGNOSTIC_SCHEMA:
        raise ValueError("scVI native training diagnostic envelope version is unsupported.")
    payload = envelope["payload"]
    if not isinstance(payload, Mapping) or set(payload) != {
        "training_diagnostics",
        "training_parameters",
        "software_versions",
        "model_identity",
    }:
        raise ValueError("scVI native training diagnostic payload schema is invalid.")
    if envelope["payload_sha256"] != _fingerprint_json(payload):
        raise ValueError("scVI native training diagnostic payload fingerprint does not match its contents.")
    diagnostics = validate_scvi_training_diagnostics(payload["training_diagnostics"])
    parameters = payload["training_parameters"]
    versions = payload["software_versions"]
    identity = payload["model_identity"]
    if not isinstance(parameters, Mapping):
        raise TypeError("scVI native training parameters must be a mapping.")
    if not isinstance(versions, Mapping) or not all(
        isinstance(key, str) and isinstance(value, str) for key, value in versions.items()
    ):
        raise TypeError("scVI native software versions must map package names to version strings.")
    if not isinstance(identity, Mapping) or set(identity) != {"path", "size", "sha256"}:
        raise ValueError("scVI native model identity schema is invalid.")
    if identity["path"] != "model.pt" or isinstance(identity["size"], bool) or not isinstance(identity["size"], int):
        raise ValueError("scVI native model identity fields are invalid.")
    if (
        not isinstance(identity["sha256"], str)
        or len(identity["sha256"]) != 64
        or any(character not in "0123456789abcdef" for character in identity["sha256"])
    ):
        raise ValueError("scVI native model identity SHA-256 is invalid.")
    model_path = root / "model.pt"
    if not model_path.is_file() or model_path.is_symlink():
        raise ValueError("scVI native artifact model.pt is missing or invalid.")
    if model_path.stat().st_size != identity["size"] or _file_sha256(model_path) != identity["sha256"]:
        raise ValueError("scVI native model identity does not match its retained training diagnostics.")
    return {
        "training_diagnostics": diagnostics,
        "training_parameters": copy.deepcopy(dict(parameters)),
        "software_versions": copy.deepcopy(dict(versions)),
        "model_identity": copy.deepcopy(dict(identity)),
    }


def _matrix_fingerprint(matrix: Any) -> str:
    """Hash one registered matrix without densifying sparse input."""
    import numpy as np

    digest = hashlib.sha256()
    shape = tuple(int(value) for value in matrix.shape)
    digest.update(json.dumps(shape, separators=(",", ":")).encode("ascii"))
    if callable(getattr(matrix, "tocsr", None)):
        canonical = matrix.tocsr(copy=True)
        canonical.sum_duplicates()
        canonical.sort_indices()
        canonical.eliminate_zeros()
        digest.update(b"csr")
        arrays = (canonical.indptr, canonical.indices, canonical.data)
    else:
        array = np.asarray(matrix)
        if array.ndim != 2:
            raise ValueError("The registered scVI expression matrix must be two-dimensional.")
        digest.update(b"dense")
        arrays = (array,)
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(contiguous.dtype.str.encode("ascii"))
        digest.update(memoryview(contiguous).cast("B"))
    return digest.hexdigest()


def _registered_count_source(adata: Any, parameters: Mapping[str, Any]) -> tuple[Any, str]:
    source = parameters["source"]
    if source == "X":
        return adata.X, "X"
    if source == "layer":
        layer = parameters.get("counts_layer")
        if not isinstance(layer, str) or not layer:
            raise ValueError("Layer-backed scVI provenance requires a nonempty counts_layer.")
        if layer not in adata.layers:
            raise ValueError(f"Registered scVI AnnData lacks declared count layer {layer!r}.")
        return adata.layers[layer], f"layers[{layer!r}]"
    raise ValueError(f"Unsupported registered scVI count source in artifact provenance: {source!r}.")


def _model_state_fingerprint(model: Any) -> str | None:
    owner = model if callable(getattr(model, "state_dict", None)) else getattr(model, "module", None)
    state_dict = getattr(owner, "state_dict", None)
    if not callable(state_dict):
        return None
    try:
        state = state_dict()
    except Exception:
        return None
    if not isinstance(state, Mapping):
        return None
    digest = hashlib.sha256()
    for key in sorted(state, key=str):
        value = state[key]
        digest.update(str(key).encode("utf-8"))
        try:
            array = value.detach().cpu().contiguous().numpy()
            digest.update(str(array.dtype).encode("ascii"))
            digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
            digest.update(memoryview(array).cast("B"))
        except Exception:
            digest.update(repr(value).encode("utf-8"))
    return digest.hexdigest()


def _registered_obs_fingerprint(adata: Any, parameters: Mapping[str, Any]) -> str:
    import pandas as pd

    keys = [
        parameters.get("technical_batch_key"),
        parameters.get("size_factor_key"),
        *(parameters.get("categorical_covariates") or ()),
        *(parameters.get("continuous_covariates") or ()),
    ]
    keys = list(dict.fromkeys(key for key in keys if isinstance(key, str) and key))
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise ValueError(f"Registered scVI AnnData lacks declared setup metadata columns: {missing}.")
    payload: dict[str, Any] = {}
    for key in keys:
        series = adata.obs[key]
        values = []
        for value in series:
            if _missing_scalar(value, pd, label=f"registered scVI metadata {key!r}"):
                values.append({"missing": True})
            else:
                scalar = _python_scalar(value)
                values.append({"type": type(scalar).__name__, "value": _json_safe(scalar)})
        payload[key] = {"dtype": str(series.dtype), "values": values}
    return _fingerprint_json(payload)


def _validate_axis(values: Sequence[Any], *, label: str) -> tuple[str, ...]:
    normalized: list[str] = []
    for index, value in enumerate(values):
        if not isinstance(value, str):
            raise TypeError(f"{label}[{index}] must be a string identifier.")
        if not value or value != value.strip():
            raise ValueError(f"{label}[{index}] must be nonempty and free of surrounding whitespace.")
        normalized.append(value)
    return tuple(normalized)


def _strict_text(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    if not value or value != value.strip():
        raise ValueError(f"{label} must be nonempty and free of surrounding whitespace.")
    return value


def _missing_scalar(value: Any, pandas: Any, *, label: str) -> bool:
    missing = pandas.isna(value)
    if isinstance(missing, bool):
        return missing
    if getattr(missing, "ndim", 1) == 0:
        try:
            return bool(missing)
        except (TypeError, ValueError):
            pass
    raise ValueError(f"{label} values must be scalar.")


def _python_scalar(value: Any) -> Any:
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    return value


def _batch_token(value: Any, pandas: Any, *, label: str) -> tuple[str, Any]:
    if _missing_scalar(value, pandas, label=label):
        raise ValueError(f"{label} cannot be missing for compared cells.")
    value = _python_scalar(value)
    if isinstance(value, bool):
        raise TypeError(f"{label} cannot use Boolean categories.")
    if isinstance(value, int):
        return "integer", value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{label} must contain finite categories.")
        return "float", value.hex()
    if isinstance(value, str):
        if not value or value != value.strip():
            raise ValueError(f"{label} string categories must be nonempty and free of surrounding whitespace.")
        return "string", value
    raise TypeError(f"{label} categories must be strings or finite numbers, not {type(value).__name__}.")


def _ordered_batch_levels(series: Any, pandas: Any, *, label: str) -> tuple[list[Any], dict[tuple[str, Any], Any]]:
    declared = list(series.cat.categories) if isinstance(series.dtype, pandas.CategoricalDtype) else list(series)
    ordered: list[Any] = []
    by_token: dict[tuple[str, Any], Any] = {}
    for raw in declared:
        if _missing_scalar(raw, pandas, label=label):
            continue
        token = _batch_token(raw, pandas, label=label)
        if token not in by_token:
            value = _python_scalar(raw)
            by_token[token] = value
            ordered.append(value)
    return ordered, by_token


def _public_keyword_support(function: Any, keywords: Sequence[str], *, description: str) -> None:
    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Cannot inspect the installed {description} API.") from error
    missing = [keyword for keyword in keywords if keyword not in signature.parameters]
    if missing:
        raise RuntimeError(f"Installed {description} does not support required arguments: {missing}.")


def _require_audited_scvi_version(scvi_module: Any) -> str:
    """Restrict model-DE semantics to the audited scvi-tools 1.5 release family."""
    from importlib import metadata

    version_value = getattr(scvi_module, "__version__", None)
    if version_value is None:
        try:
            version_value = metadata.version("scvi-tools")
        except metadata.PackageNotFoundError as exc:
            raise RuntimeError("Cannot determine the installed scvi-tools version.") from exc
    version = str(version_value)
    numeric = version.split("+", 1)[0].split(".post", 1)[0].split(".dev", 1)[0]
    try:
        major_minor = tuple(int(part) for part in numeric.split(".")[:2])
    except ValueError as exc:
        raise RuntimeError(f"Cannot parse installed scvi-tools version {version!r}.") from exc
    if major_minor != (1, 5):
        raise RuntimeError(f"scVI Model DE Evidence requires audited scvi-tools >=1.5,<1.6; detected {version!r}.")
    return version


class SCVIModel:
    """A trained, process-local scVI model bound to an immutable registered expression view."""

    __slots__ = (
        "_artifact_fingerprint",
        "_count_fingerprint",
        "_diagnostics",
        "_feature_fingerprint",
        "_lock",
        "_model",
        "_model_class",
        "_observation_fingerprint",
        "_obs_name_set",
        "_obs_names",
        "_registered_adata",
        "_registered_count_source",
        "_registered_obs_fingerprint",
        "_state_fingerprint",
        "_training_parameters",
        "_var_names",
    )

    def __init__(
        self,
        model: Any,
        registered_adata: AnnData,
        training_parameters: Mapping[str, Any],
        diagnostics: Mapping[str, Any] | None = None,
    ) -> None:
        if not callable(getattr(model, "differential_expression", None)):
            raise TypeError("SCVIModel requires a scVI model with differential-expression support.")
        if not callable(getattr(model, "deregister_manager", None)):
            raise TypeError("SCVIModel requires a scVI model with AnnData manager lifecycle support.")
        if getattr(model, "adata", None) is not registered_adata:
            raise ValueError("The scVI model must be attached to the registered AnnData object.")
        if not bool(getattr(model, "is_trained", False)):
            raise ValueError("SCVIModel requires a trained scVI model.")
        if not isinstance(training_parameters, Mapping):
            raise TypeError("SCVIModel training_parameters must be a mapping.")
        required_provenance = ("source",)
        missing_provenance = [key for key in required_provenance if key not in training_parameters]
        if missing_provenance:
            raise ValueError(
                "SCVIModel training_parameters are missing required current provenance fields: "
                f"{missing_provenance!r}."
            )
        if diagnostics is not None and not isinstance(diagnostics, Mapping):
            raise TypeError("SCVIModel diagnostics must be a mapping or None.")

        try:
            obs_names = _validate_axis(tuple(registered_adata.obs_names), label="registered obs_names")
            var_names = _validate_axis(tuple(registered_adata.var_names), label="registered var_names")
            obs_names_are_unique = bool(registered_adata.obs_names.is_unique)
            var_names_are_unique = bool(registered_adata.var_names.is_unique)
        except AttributeError as error:
            raise TypeError("SCVIModel requires a registered AnnData object.") from error
        if not obs_names or not var_names:
            raise ValueError("The registered scVI AnnData must contain observations and variables.")
        if not obs_names_are_unique or not var_names_are_unique:
            raise ValueError("The registered scVI AnnData requires unique obs_names and var_names.")

        self._model = model
        self._registered_adata = registered_adata
        self._training_parameters = copy.deepcopy(dict(training_parameters))
        self._diagnostics = (
            validate_scvi_training_diagnostics(diagnostics) if diagnostics else {}
        )
        _json_safe(self._training_parameters)
        _json_safe(self._diagnostics)
        self._lock = threading.RLock()
        self._obs_names = obs_names
        self._obs_name_set = frozenset(obs_names)
        self._var_names = var_names
        model_type = type(model)
        self._model_class = f"{model_type.__module__}.{model_type.__qualname__}"
        self._observation_fingerprint = _fingerprint_json(obs_names)
        self._feature_fingerprint = _fingerprint_json(var_names)
        registered_counts, self._registered_count_source = _registered_count_source(
            registered_adata, self._training_parameters
        )
        self._count_fingerprint = _matrix_fingerprint(registered_counts)
        self._registered_obs_fingerprint = _registered_obs_fingerprint(registered_adata, self._training_parameters)
        self._state_fingerprint = _model_state_fingerprint(model)
        self._artifact_fingerprint = _fingerprint_json(
            {
                "schema": SCVI_MODEL_ARTIFACT_SCHEMA,
                "model_class": self._model_class,
                "observations": self._observation_fingerprint,
                "features": self._feature_fingerprint,
                "registered_counts": self._count_fingerprint,
                "registered_setup_metadata": self._registered_obs_fingerprint,
                "fitted_state": self._state_fingerprint,
                "training_parameters": self._training_parameters,
                "diagnostics": self._diagnostics,
            }
        )

    @property
    def registered_adata(self) -> AnnData:
        with self._lock:
            self._validate_registered_state()
            return self._registered_adata.copy()

    @property
    def training_parameters(self) -> Mapping[str, Any]:
        return MappingProxyType(copy.deepcopy(self._training_parameters))

    @property
    def model_class(self) -> str:
        return self._model_class

    @property
    def diagnostics(self) -> Mapping[str, Any]:
        with self._lock:
            self._validate_registered_state()
            return MappingProxyType(copy.deepcopy(self._diagnostics))

    @property
    def obs_names(self) -> tuple[str, ...]:
        return self._obs_names

    @property
    def var_names(self) -> tuple[str, ...]:
        return self._var_names

    @property
    def evidence(self) -> Mapping[str, Any]:
        with self._lock:
            self._validate_registered_state()
            return MappingProxyType(
                {
                    "artifact_schema_version": SCVI_MODEL_ARTIFACT_SCHEMA,
                    "artifact_fingerprint_sha256": self._artifact_fingerprint,
                    "model_class": self._model_class,
                    "fitted_state_fingerprint_sha256": self._state_fingerprint,
                    "registered_count_fingerprint_sha256": self._count_fingerprint,
                    "registered_count_source": self._registered_count_source,
                    "registered_setup_metadata_fingerprint_sha256": self._registered_obs_fingerprint,
                    "observation_axis_fingerprint_sha256": self._observation_fingerprint,
                    "feature_axis_fingerprint_sha256": self._feature_fingerprint,
                    "training_observations": len(self._obs_names),
                    "fitted_features": len(self._var_names),
                    "training_parameters": copy.deepcopy(self._training_parameters),
                    "training_diagnostics": copy.deepcopy(self._diagnostics),
                    "process_local": True,
                    "serializable_fitted_weights": False,
                }
            )

    def differential_expression_evidence(
        self,
        adata: AnnData,
        *,
        groupby: str,
        group1: str,
        group2: str,
        subset_column: str | None,
        subset_value: str | None,
        mode: str,
        delta: float | None,
        fdr_target: float | None,
        batch_handling: str,
        n_samples_overall: int,
        random_seed: int,
    ) -> SCVIDifferentialRun:
        """Run one validated two-population contrast using the public scvi-tools API."""
        import numpy as np
        import pandas as pd

        groupby = _strict_text(groupby, label="scVI groupby")
        group1 = _strict_text(group1, label="scVI group1")
        group2 = _strict_text(group2, label="scVI group2")
        if group1 == group2:
            raise ValueError("scVI group1 and group2 must be different.")
        if (subset_column is None) != (subset_value is None):
            raise ValueError("scVI subset_column and subset_value must either both be absent or both be present.")
        if subset_column is not None:
            subset_column = _strict_text(subset_column, label="scVI subset_column")
            subset_value = _strict_text(subset_value, label="scVI subset_value")
        if mode not in {"change", "vanilla"}:
            raise ValueError(f"Unsupported scVI model-DE mode: {mode!r}.")
        if mode == "change":
            if isinstance(delta, bool) or not isinstance(delta, (int, float)) or not math.isfinite(float(delta)):
                raise TypeError("scVI change-mode delta must be a finite positive number.")
            if float(delta) <= 0:
                raise ValueError("scVI change-mode delta must be positive.")
            if (
                isinstance(fdr_target, bool)
                or not isinstance(fdr_target, (int, float))
                or not math.isfinite(float(fdr_target))
            ):
                raise TypeError("scVI change-mode fdr_target must be a finite number between zero and one.")
            if not 0 < float(fdr_target) < 1:
                raise ValueError("scVI change-mode fdr_target must be strictly between zero and one.")
            delta = float(delta)
            fdr_target = float(fdr_target)
        elif delta is not None or fdr_target is not None:
            raise ValueError("scVI vanilla mode does not accept delta or fdr_target parameters.")
        if batch_handling not in {"shared_technical_batches", "observed_technical_batches"}:
            raise ValueError(f"Unsupported scVI batch handling: {batch_handling!r}.")
        if isinstance(n_samples_overall, bool) or not isinstance(n_samples_overall, int) or n_samples_overall < 1:
            raise ValueError("scVI n_samples_overall must be a positive integer.")
        if isinstance(random_seed, bool) or not isinstance(random_seed, int) or not 0 <= random_seed <= 2**31 - 1:
            raise ValueError("scVI random_seed must be an integer in [0, 2**31 - 1].")

        parameters = self._training_parameters
        technical_batch_key = parameters.get("technical_batch_key")
        nuisance_keys = {
            key
            for key in (
                technical_batch_key,
                parameters.get("size_factor_key"),
                *(parameters.get("categorical_covariates") or ()),
                *(parameters.get("continuous_covariates") or ()),
            )
            if isinstance(key, str) and key
        }
        if groupby in nuisance_keys:
            raise ValueError(
                f"scVI groupby {groupby!r} was declared as a Technical/nuisance covariate during model training; "
                "it cannot simultaneously define the biological comparison."
            )
        if batch_handling == "shared_technical_batches" and (
            not isinstance(technical_batch_key, str) or not technical_batch_key
        ):
            raise ValueError(
                "Shared-Technical-batch scVI evidence requires a primary technical_batch_key in model provenance."
            )

        copy_columns = [groupby]
        if subset_column is not None and subset_column != technical_batch_key:
            copy_columns.append(subset_column)
        copy_columns = list(dict.fromkeys(copy_columns))

        with SCVI_GLOBAL_RNG_LOCK, self._lock:
            self._validate_registered_state()
            analysis_adata = self._copy_compatible_adata(adata, copy_columns)
            if isinstance(technical_batch_key, str) and technical_batch_key and technical_batch_key not in analysis_adata.obs:
                raise RuntimeError(
                    f"The registered scVI AnnData no longer contains Technical batch {technical_batch_key!r}."
                )
            total_cells = int(analysis_adata.n_obs)

            if subset_column is not None:
                if subset_column not in analysis_adata.obs:
                    raise ValueError(f"scVI subset observation column not found: {subset_column!r}.")
                subset_series = analysis_adata.obs[subset_column]
                scope_mask: list[bool] = []
                for value in subset_series:
                    if _missing_scalar(value, pd, label=f"scVI subset column {subset_column!r}"):
                        scope_mask.append(False)
                        continue
                    value = _python_scalar(value)
                    if not isinstance(value, str) or not value or value != value.strip():
                        raise TypeError(
                            f"scVI subset column {subset_column!r} must contain canonical string labels or missing values."
                        )
                    scope_mask.append(value == subset_value)
                if not any(scope_mask):
                    raise ValueError(f"scVI subset contains no observations: {subset_column}={subset_value!r}.")
                analysis_adata = analysis_adata[np.asarray(scope_mask, dtype=bool)].copy()

            group_series = analysis_adata.obs[groupby]
            if isinstance(group_series.dtype, pd.CategoricalDtype):
                declared_group_categories = []
                for category in group_series.cat.categories:
                    category = _python_scalar(category)
                    if not isinstance(category, str) or not category or category != category.strip():
                        raise TypeError(
                            f"scVI groupby {groupby!r} declared categories must be canonical strings."
                        )
                    declared_group_categories.append(category)
            else:
                declared_group_categories = []
            group1_mask: list[bool] = []
            group2_mask: list[bool] = []
            missing_group_mask: list[bool] = []
            observed_groups: set[str] = set()
            for value in group_series:
                if _missing_scalar(value, pd, label=f"scVI groupby {groupby!r}"):
                    group1_mask.append(False)
                    group2_mask.append(False)
                    missing_group_mask.append(True)
                    continue
                value = _python_scalar(value)
                if not isinstance(value, str) or not value or value != value.strip():
                    raise TypeError(f"scVI groupby {groupby!r} must contain canonical string labels or missing values.")
                observed_groups.add(value)
                if not isinstance(group_series.dtype, pd.CategoricalDtype) and value not in declared_group_categories:
                    declared_group_categories.append(value)
                group1_mask.append(value == group1)
                group2_mask.append(value == group2)
                missing_group_mask.append(False)
            missing_groups = [value for value in (group1, group2) if value not in observed_groups]
            if missing_groups:
                raise ValueError(f"scVI groups not found in the selected population scope: {missing_groups}.")
            group1_array = np.asarray(group1_mask, dtype=bool)
            group2_array = np.asarray(group2_mask, dtype=bool)
            missing_group_array = np.asarray(missing_group_mask, dtype=bool)
            group1_cells = int(group1_array.sum())
            group2_cells = int(group2_array.sum())
            if group1_cells < 2 or group2_cells < 2:
                raise ValueError(
                    "scVI model-DE requires at least two selected cells in each population; "
                    f"observed {group1!r}={group1_cells}, {group2!r}={group2_cells}."
                )

            batch_evidence, shared_batches = self._resolve_batch_evidence(
                analysis_adata,
                technical_batch_key=technical_batch_key,
                group1=group1,
                group2=group2,
                group1_mask=group1_array,
                group2_mask=group2_array,
                batch_handling=batch_handling,
                n_samples_overall=n_samples_overall,
                pandas=pd,
            )

            backend_kwargs: dict[str, Any] = {
                "adata": analysis_adata,
                "groupby": groupby,
                "group1": [group1],
                "group2": group2,
                "mode": mode,
                "batch_size": None,
                "all_stats": True,
                "batch_correction": batch_handling == "shared_technical_batches",
                "silent": True,
                "weights": "uniform",
                "filter_outlier_cells": False,
                "importance_weighting_kwargs": None,
                "dataloader": None,
                "n_samples_overall": n_samples_overall,
                "use_permutation": False,
                "pseudocounts": None,
            }
            if batch_handling == "shared_technical_batches":
                backend_kwargs["batchid1"] = shared_batches
                backend_kwargs["batchid2"] = shared_batches
            if mode == "change":
                backend_kwargs.update(
                    {
                        "delta": delta,
                        "fdr_target": fdr_target,
                        "change_fn": None,
                        "m1_domain_fn": None,
                        "test_mode": "two",
                    }
                )

            # Lazy scvi/torch imports may initialize global RNG-backed subsystems. Preserve
            # Python/NumPy state across that import as well as across posterior sampling.
            pre_import_python_state = random.getstate()
            pre_import_numpy_state = np.random.get_state()
            try:
                import scvi
                import torch
                from scvi.model.base import DifferentialComputation

                scvi_tools_version = _require_audited_scvi_version(scvi)
                _public_keyword_support(
                    self._model.differential_expression,
                    (
                        "adata",
                        "groupby",
                        "group1",
                        "group2",
                        "mode",
                        "batch_size",
                        "all_stats",
                        "batch_correction",
                        "batchid1",
                        "batchid2",
                        "delta",
                        "fdr_target",
                        "silent",
                        "weights",
                        "filter_outlier_cells",
                        "importance_weighting_kwargs",
                        "dataloader",
                    ),
                    description="scvi.model.SCVI.differential_expression",
                )
                _public_keyword_support(
                    DifferentialComputation.get_bayes_factors,
                    ("n_samples_overall", "use_permutation", "pseudocounts"),
                    description="scvi.model.base.DifferentialComputation.get_bayes_factors",
                )
                if mode == "change":
                    _public_keyword_support(
                        DifferentialComputation.get_bayes_factors,
                        ("change_fn", "m1_domain_fn", "test_mode"),
                        description="scvi.model.base.DifferentialComputation.get_bayes_factors",
                    )
                table = self._run_backend_with_scoped_rng(
                    backend_kwargs, random_seed=random_seed, numpy=np, scvi=scvi, torch=torch
                )
            finally:
                random.setstate(pre_import_python_state)
                np.random.set_state(pre_import_numpy_state)
            selected_cells = int(analysis_adata.n_obs)
            comparison_evidence = {
                "groupby": groupby,
                "group1": group1,
                "group2": group2,
                "population_scope": (
                    {"kind": "all"}
                    if subset_column is None
                    else {"kind": "obs_value", "column": subset_column, "value": subset_value}
                ),
                "input_cells": total_cells,
                "selected_scope_cells": selected_cells,
                "excluded_by_scope_cells": total_cells - selected_cells,
                "group1_cells": group1_cells,
                "group2_cells": group2_cells,
                "missing_group_label_cells": int(missing_group_array.sum()),
                "other_group_label_cells": selected_cells
                - group1_cells
                - group2_cells
                - int(missing_group_array.sum()),
                "declared_group_categories": declared_group_categories,
                "observed_group_categories": [
                    value for value in declared_group_categories if value in observed_groups
                ],
                "unused_declared_group_categories": [
                    value for value in declared_group_categories if value not in observed_groups
                ],
                **batch_evidence,
            }
            disclosed_backend = {
                key: _json_safe(value) for key, value in backend_kwargs.items() if key != "adata"
            }
            disclosed_backend["scvi_tools_version"] = scvi_tools_version
            return SCVIDifferentialRun(
                table=table,
                comparison_evidence=MappingProxyType(comparison_evidence),
                backend_parameters=MappingProxyType(disclosed_backend),
            )

    def _resolve_batch_evidence(
        self,
        analysis_adata: AnnData,
        *,
        technical_batch_key: Any,
        group1: str,
        group2: str,
        group1_mask: Any,
        group2_mask: Any,
        batch_handling: str,
        n_samples_overall: int,
        pandas: Any,
    ) -> tuple[dict[str, Any], list[Any]]:
        evidence = {
            "technical_batch_key": technical_batch_key or None,
            "batch_handling": batch_handling,
            "registered_technical_batches": [],
            "group1_observed_technical_batches": [],
            "group2_observed_technical_batches": [],
            "shared_technical_batches_used": [],
            "excluded_nonshared_technical_batches": [],
            "group_by_technical_batch_counts": [],
        }
        if not isinstance(technical_batch_key, str) or not technical_batch_key:
            return evidence, []

        label = f"Technical batch {technical_batch_key!r}"
        batch_series = analysis_adata.obs[technical_batch_key]
        registered_series = self._registered_adata.obs[technical_batch_key]
        registered_levels, registered_by_token = _ordered_batch_levels(registered_series, pandas, label=label)
        registered_tokens = [_batch_token(value, pandas, label=label) for value in registered_levels]
        group_tokens: dict[str, set[tuple[str, Any]]] = {group1: set(), group2: set()}
        counts: dict[tuple[str, tuple[str, Any]], int] = {}
        for row, raw_batch in enumerate(batch_series):
            if not (group1_mask[row] or group2_mask[row]):
                continue
            token = _batch_token(raw_batch, pandas, label=label)
            if token not in registered_by_token:
                raise RuntimeError("Selected Technical-batch metadata is incompatible with model registration.")
            group = group1 if group1_mask[row] else group2
            group_tokens[group].add(token)
            counts[(group, token)] = counts.get((group, token), 0) + 1
        group1_levels = [registered_by_token[token] for token in registered_tokens if token in group_tokens[group1]]
        group2_levels = [registered_by_token[token] for token in registered_tokens if token in group_tokens[group2]]
        shared_tokens = [
            token for token in registered_tokens if token in group_tokens[group1] and token in group_tokens[group2]
        ]
        shared_levels = [registered_by_token[token] for token in shared_tokens]
        excluded_levels = [
            registered_by_token[token]
            for token in registered_tokens
            if (token in group_tokens[group1]) ^ (token in group_tokens[group2])
        ]
        count_rows = [
            {
                "technical_batch": _json_safe(registered_by_token[token]),
                "group1_cells": counts.get((group1, token), 0),
                "group2_cells": counts.get((group2, token), 0),
            }
            for token in registered_tokens
            if token in group_tokens[group1] or token in group_tokens[group2]
        ]
        if batch_handling == "shared_technical_batches" and not shared_levels:
            raise ValueError(
                "The two scVI populations share no observed Technical-batch level. Counterfactual decoding would be "
                "unsupported; use Sample-level replicate-aware inference for Condition questions."
            )
        if batch_handling == "shared_technical_batches" and n_samples_overall < 2 * len(shared_levels):
            raise ValueError(
                "scVI n_samples_overall must provide at least two posterior samples per shared Technical batch."
            )
        evidence.update(
            {
                "registered_technical_batches": [_json_safe(value) for value in registered_levels],
                "group1_observed_technical_batches": [_json_safe(value) for value in group1_levels],
                "group2_observed_technical_batches": [_json_safe(value) for value in group2_levels],
                "shared_technical_batches_used": [_json_safe(value) for value in shared_levels],
                "excluded_nonshared_technical_batches": [_json_safe(value) for value in excluded_levels],
                "group_by_technical_batch_counts": count_rows,
            }
        )
        return evidence, shared_levels

    def _run_backend_with_scoped_rng(
        self, backend_kwargs: Mapping[str, Any], *, random_seed: int, numpy: Any, scvi: Any, torch: Any
    ) -> DataFrame:
        analysis_adata = backend_kwargs.get("adata")
        if analysis_adata is None:
            raise RuntimeError("scVI model-DE backend call requires one exact analysis AnnData for scoped cleanup.")
        previous_scvi_seed = getattr(scvi.settings, "seed", None)
        python_state = random.getstate()
        numpy_state = numpy.random.get_state()
        torch_state = torch.random.get_rng_state()
        cuda_states = torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None
        mps_module = getattr(torch, "mps", None)
        mps_backend = getattr(getattr(torch, "backends", None), "mps", None)
        mps_available = (
            mps_module is not None
            and callable(getattr(mps_module, "get_rng_state", None))
            and callable(getattr(mps_module, "set_rng_state", None))
            and mps_backend is not None
            and callable(getattr(mps_backend, "is_available", None))
            and bool(mps_backend.is_available())
        )
        mps_state = mps_module.get_rng_state() if mps_available else None
        try:
            with torch.inference_mode(False):
                scvi.settings.seed = random_seed
                return self._model.differential_expression(**backend_kwargs)
        finally:
            differential_error = sys.exception()
            try:
                self._deregister_analysis_manager(analysis_adata)
            except Exception:
                if differential_error is None:
                    raise
            finally:
                scvi.settings.seed = previous_scvi_seed
                random.setstate(python_state)
                numpy.random.set_state(numpy_state)
                torch.random.set_rng_state(torch_state)
                if cuda_states is not None:
                    torch.cuda.set_rng_state_all(cuda_states)
                if mps_state is not None:
                    mps_module.set_rng_state(mps_state)

    def _deregister_analysis_manager(self, analysis_adata: AnnData) -> None:
        """Remove exactly one transferred scvi-tools 1.5 manager through public APIs."""
        get_manager = getattr(self._model, "get_anndata_manager", None)
        register_manager = getattr(self._model, "register_manager", None)
        if callable(get_manager) and callable(register_manager):
            manager = get_manager(analysis_adata, required=False)
            if manager is None:
                return
            # scvi-tools 1.5 registers transferred fields only in the per-instance
            # store, while deregister_manager(adata) resolves through the class store.
            # Publishing this exact manager through the public class API makes the
            # documented targeted deregistration usable without touching either
            # private manager dictionary or any unrelated manager.
            register_manager(manager)
        self._model.deregister_manager(analysis_adata)

    def _copy_compatible_adata(self, adata: AnnData, required_obs_columns: Sequence[str]) -> AnnData:
        try:
            downstream_obs_names = _validate_axis(tuple(adata.obs_names), label="downstream obs_names")
            downstream_var_names = _validate_axis(tuple(adata.var_names), label="downstream var_names")
            obs_names_are_unique = bool(adata.obs_names.is_unique)
            var_names_are_unique = bool(adata.var_names.is_unique)
        except AttributeError as error:
            raise TypeError("SCVIModel requires a downstream AnnData object.") from error
        if not downstream_obs_names:
            raise ValueError("Downstream AnnData must contain at least one observation.")
        if not obs_names_are_unique or not var_names_are_unique:
            raise ValueError("Downstream AnnData requires unique obs_names and var_names.")
        if downstream_var_names != self._var_names:
            raise ValueError("Downstream AnnData var_names are incompatible with the trained scVI model.")
        unknown = [name for name in downstream_obs_names if name not in self._obs_name_set]
        if unknown:
            raise ValueError(
                "Downstream AnnData obs_names are incompatible with the trained scVI model; "
                f"first unknown observation: {unknown[0]!r}."
            )
        missing_columns = [column for column in required_obs_columns if column not in adata.obs]
        if missing_columns:
            raise ValueError(f"scVI differential-expression observation columns not found: {missing_columns}")

        analysis_adata = self._registered_adata[list(downstream_obs_names), :].copy()
        for column in required_obs_columns:
            analysis_adata.obs[column] = adata.obs[column].reindex(analysis_adata.obs_names).copy()
        return analysis_adata

    def _validate_registered_state(self) -> None:
        if getattr(self._model, "adata", None) is not self._registered_adata:
            raise RuntimeError("The scVI model is no longer attached to its registered AnnData object.")
        if not bool(getattr(self._model, "is_trained", False)):
            raise RuntimeError("The scVI model is no longer in a trained state.")
        if tuple(self._registered_adata.obs_names) != self._obs_names:
            raise RuntimeError("The registered scVI AnnData obs_names changed after training.")
        if tuple(self._registered_adata.var_names) != self._var_names:
            raise RuntimeError("The registered scVI AnnData var_names changed after training.")
        registered_counts, source = _registered_count_source(self._registered_adata, self._training_parameters)
        if source != self._registered_count_source or _matrix_fingerprint(registered_counts) != self._count_fingerprint:
            raise RuntimeError("The registered scVI count matrix changed after training.")
        if _registered_obs_fingerprint(self._registered_adata, self._training_parameters) != self._registered_obs_fingerprint:
            raise RuntimeError("The registered scVI setup metadata changed after training.")
        current_state = _model_state_fingerprint(self._model)
        if self._state_fingerprint is not None and current_state != self._state_fingerprint:
            raise RuntimeError("The fitted scVI model weights changed after artifact construction.")


__all__ = [
    "SCVI_EPOCH_METRICS",
    "SCVI_GLOBAL_RNG_LOCK",
    "SCVI_MODEL_ARTIFACT_SCHEMA",
    "SCVI_NATIVE_DIAGNOSTIC_SCHEMA",
    "SCVI_TRAINING_DIAGNOSTIC_FILENAME",
    "SCVI_TRAINING_DIAGNOSTIC_SCHEMA",
    "SCVIDifferentialRun",
    "SCVIModel",
    "build_scvi_training_diagnostics",
    "read_scvi_training_diagnostics",
    "validate_scvi_training_diagnostics",
    "write_scvi_training_diagnostics",
]
