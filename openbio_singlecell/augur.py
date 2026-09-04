from __future__ import annotations

import copy
import hashlib
import inspect
import json
import math
import re
import warnings as python_warnings
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .analysis_reporting import collect_software_versions

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame


AUGUR_ARTIFACT_TYPE = "OPENBIO_AUGUR_RESULT"
AUGUR_ARTIFACT_SCHEMA_VERSION = 2
AUGUR_PRODUCER_NODE_ID = "OpenBioSingleCellAugur"
AUGUR_RESULTS_NODE_ID = "OpenBioSingleCellAugurResults"
AUGUR_PLOT_NODE_ID = "OpenBioSingleCellAugurPlot"
AUGUR_PRODUCER_SCHEMA = "openbio-singlecell/augur-result/v2"
AUGUR_SUMMARY_SCHEMA = "openbio-singlecell/augur-summary/v2"
AUGUR_VIEWS = ("priorities", "cross_validation", "feature_importance", "predictions")
AUGUR_PLOT_VIEWS = ("priorities", "cross_validation", "feature_importance", "roc")
AUGUR_CLASSIFIERS = ("random_forest_classifier", "logistic_regression_classifier")

PRIORITY_COLUMNS = (
    "rank",
    "population",
    "mean_augur_score",
    "mean_auc",
    "mean_accuracy",
    "mean_precision",
    "mean_f1",
    "mean_recall",
    "cells_control",
    "cells_treatment",
    "samples_control",
    "samples_treatment",
)
CROSS_VALIDATION_COLUMNS = ("population", "subsample", "fold", "auc")
FEATURE_IMPORTANCE_COLUMNS = ("population", "subsample", "fold", "gene", "importance")
PREDICTION_COLUMNS = (
    "population",
    "subsample",
    "fold",
    "observation",
    "true_label",
    "prediction_score",
)

_PRIORITY_FLOAT_COLUMNS = (
    "mean_augur_score",
    "mean_auc",
    "mean_accuracy",
    "mean_precision",
    "mean_f1",
    "mean_recall",
)
_PRIORITY_INTEGER_COLUMNS = (
    "rank",
    "cells_control",
    "cells_treatment",
    "samples_control",
    "samples_treatment",
)

AUGUR_REFERENCES = [
    {
        "citation": (
            "Skinnider MA, Squair JW, Kathe C, et al. Cell type prioritization in single-cell data. "
            "Nature Biotechnology. 2021;39:30-34."
        ),
        "doi": "10.1038/s41587-020-0605-1",
        "url": "https://doi.org/10.1038/s41587-020-0605-1",
        "kind": "method",
    },
    {
        "citation": (
            "Squair JW, Gautier M, Kathe C, et al. Prioritization of cell types responsive to biological "
            "perturbations with Augur. Nature Protocols. 2021;16:3836-3873."
        ),
        "doi": "10.1038/s41596-021-00561-x",
        "url": "https://doi.org/10.1038/s41596-021-00561-x",
        "kind": "method",
    },
    {
        "citation": (
            "Heumos L, Schaar AC, Lance C, et al. pertpy: an end-to-end framework for perturbation analysis. "
            "Nature Methods. 2025."
        ),
        "doi": "10.1038/s41592-025-02909-7",
        "url": "https://doi.org/10.1038/s41592-025-02909-7",
        "kind": "software",
    },
    {
        "citation": "Pertpy 1.3.0 Augur public interface and implementation.",
        "doi": None,
        "url": "https://github.com/scverse/pertpy/blob/v1.3.0/src/pertpy/tools/_augur.py",
        "kind": "software_documentation",
    },
]


def _canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _strict_text(value: Any, *, label: str, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise TypeError(f"Augur {label} must be a string.")
    if value != value.strip():
        raise ValueError(f"Augur {label} must not contain surrounding whitespace.")
    if not value and not allow_empty:
        raise ValueError(f"Augur {label} cannot be empty.")
    return value


def _canonical_labels(series: Any, *, label: str, pandas: Any) -> tuple[list[str], list[str]]:
    values: list[str] = []
    observed: list[str] = []
    seen: set[str] = set()
    for position, value in enumerate(series.tolist()):
        missing = pandas.isna(value)
        if not isinstance(missing, bool):
            try:
                missing = bool(missing)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Augur {label} at row {position} must be scalar.") from exc
        if missing:
            raise ValueError(f"Augur {label} contains missing values.")
        if not isinstance(value, str) or not value or value != value.strip():
            raise TypeError(
                f"Augur {label} must contain canonical nonempty string labels without coercion."
            )
        values.append(value)
        if value not in seen:
            seen.add(value)
            observed.append(value)
    return values, observed


def _count_matrix_fingerprint(
    matrix: Any,
    *,
    observation_names: Sequence[str],
    feature_names: Sequence[str],
    sparse: Any,
    numpy: Any,
) -> str:
    digest = hashlib.sha256()
    header = {
        "schema": "openbio-singlecell/augur-count-matrix/v1",
        "shape": [len(observation_names), len(feature_names)],
        "observations": list(observation_names),
        "features": list(feature_names),
    }
    digest.update(
        json.dumps(header, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True).encode("utf-8")
    )
    if sparse.issparse(matrix):
        canonical = sparse.csr_matrix(matrix, dtype=numpy.int64, copy=True)
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
        for array in (canonical.indptr, canonical.indices, canonical.data):
            values = numpy.ascontiguousarray(array, dtype="<i8")
            digest.update(int(values.size).to_bytes(8, "little", signed=False))
            digest.update(values.tobytes(order="C"))
    else:
        values = numpy.ascontiguousarray(matrix, dtype="<i8")
        digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _observation_fingerprint(
    *,
    observation_names: Sequence[str],
    samples: Sequence[str],
    populations: Sequence[str],
    conditions: Sequence[str],
    technical_batches: Sequence[str] | None,
) -> str:
    return _canonical_json_sha256(
        {
            "schema": "openbio-singlecell/augur-observations/v1",
            "observation_names": list(observation_names),
            "samples": list(samples),
            "populations": list(populations),
            "conditions": list(conditions),
            "technical_batches": None if technical_batches is None else list(technical_batches),
        }
    )


def _append_string(digest: Any, value: str) -> None:
    encoded = value.encode("utf-8")
    digest.update(len(encoded).to_bytes(8, "little", signed=False))
    digest.update(encoded)


def _frame_fingerprint(frame: Any, *, view: str, numpy: Any, pandas: Any) -> str:
    if not isinstance(frame, pandas.DataFrame):
        raise TypeError(f"Augur {view} table must be a pandas DataFrame.")
    expected = {
        "priorities": PRIORITY_COLUMNS,
        "cross_validation": CROSS_VALIDATION_COLUMNS,
        "feature_importance": FEATURE_IMPORTANCE_COLUMNS,
        "predictions": PREDICTION_COLUMNS,
    }[view]
    if tuple(frame.columns) != expected:
        raise ValueError(
            f"Augur {view} table columns differ from the canonical schema: {list(frame.columns)!r}."
        )
    digest = hashlib.sha256()
    digest.update(f"openbio-singlecell/augur-{view}/v1\0".encode("ascii"))
    digest.update(int(len(frame)).to_bytes(8, "little", signed=False))
    for column in expected:
        _append_string(digest, column)
        if column in {"population", "gene", "observation"}:
            for value in frame[column].tolist():
                if not isinstance(value, str) or not value or value != value.strip():
                    raise ValueError(f"Augur {view} column {column!r} contains a noncanonical identifier.")
                _append_string(digest, value)
        elif column in {
            "rank",
            "cells_control",
            "cells_treatment",
            "samples_control",
            "samples_treatment",
            "subsample",
            "fold",
            "true_label",
        }:
            if not pandas.api.types.is_integer_dtype(frame[column].dtype):
                raise TypeError(f"Augur {view} column {column!r} must use an integer dtype.")
            try:
                values = frame[column].to_numpy(dtype=numpy.int64, copy=True)
            except (TypeError, ValueError, OverflowError) as exc:
                raise TypeError(f"Augur {view} column {column!r} must be integer-valued.") from exc
            original = frame[column].to_numpy(copy=False)
            try:
                original_numeric = numpy.asarray(original, dtype=float)
            except (TypeError, ValueError) as exc:
                raise TypeError(f"Augur {view} column {column!r} must be integer-valued.") from exc
            if not bool(numpy.isfinite(original_numeric).all()) or not bool(
                numpy.equal(original_numeric, values.astype(float)).all()
            ):
                raise ValueError(f"Augur {view} column {column!r} must contain finite exact integers.")
            digest.update(numpy.ascontiguousarray(values, dtype="<i8").tobytes(order="C"))
        else:
            if (
                not pandas.api.types.is_numeric_dtype(frame[column].dtype)
                or pandas.api.types.is_bool_dtype(frame[column].dtype)
            ):
                raise TypeError(f"Augur {view} column {column!r} must use a numeric non-boolean dtype.")
            try:
                values = frame[column].to_numpy(dtype=float, copy=True)
            except (TypeError, ValueError) as exc:
                raise TypeError(f"Augur {view} column {column!r} must be numeric.") from exc
            if not bool(numpy.isfinite(values).all()):
                raise ValueError(f"Augur {view} column {column!r} contains non-finite values.")
            digest.update(numpy.ascontiguousarray(values, dtype="<f8").tobytes(order="C"))
    return digest.hexdigest()


def _validate_canonical_tables(
    tables: Mapping[str, Any],
    *,
    parameters: Mapping[str, Any],
    numpy: Any,
    pandas: Any,
) -> dict[str, str]:
    if set(tables) != set(AUGUR_VIEWS):
        raise ValueError(
            f"Augur artifact table family is invalid; expected={list(AUGUR_VIEWS)}, observed={sorted(tables)}."
        )
    priorities = tables["priorities"]
    cross_validation = tables["cross_validation"]
    feature_importance = tables["feature_importance"]
    predictions = tables["predictions"]
    fingerprints = {
        view: _frame_fingerprint(tables[view], view=view, numpy=numpy, pandas=pandas) for view in AUGUR_VIEWS
    }

    if priorities.empty:
        raise ValueError("Augur priorities table cannot be empty.")
    populations = priorities["population"].tolist()
    if len(populations) != len(set(populations)):
        raise ValueError("Augur priorities must contain exactly one row per population.")
    for column in _PRIORITY_FLOAT_COLUMNS:
        values = priorities[column].to_numpy(dtype=float)
        if bool(((values < 0.0) | (values > 1.0)).any()):
            raise ValueError(f"Augur priority metric {column!r} must lie in [0, 1].")
    for column in _PRIORITY_INTEGER_COLUMNS:
        minimum = 1
        if bool((priorities[column].to_numpy(dtype=int) < minimum).any()):
            raise ValueError(f"Augur priority support column {column!r} must be >= {minimum}.")
    expected_ranks = (
        priorities["mean_augur_score"].rank(method="min", ascending=False).astype(int).tolist()
    )
    if priorities["rank"].astype(int).tolist() != expected_ranks:
        raise ValueError("Augur priority ranks do not match deterministic minimum ranks of mean Augur scores.")
    expected_priority_order = priorities.sort_values(
        ["rank", "population"], ascending=[True, True], kind="mergesort"
    )["population"].tolist()
    if populations != expected_priority_order:
        raise ValueError("Augur priorities are not in canonical rank/population order.")

    n_subsamples = parameters.get("n_subsamples")
    folds = parameters.get("folds")
    if isinstance(n_subsamples, bool) or not isinstance(n_subsamples, int) or n_subsamples < 1:
        raise ValueError("Augur artifact parameters contain an invalid n_subsamples.")
    if isinstance(folds, bool) or not isinstance(folds, int) or folds < 2:
        raise ValueError("Augur artifact parameters contain invalid folds.")
    if cross_validation.empty:
        raise ValueError("Augur cross-validation table cannot be empty.")
    if cross_validation.duplicated(["population", "subsample", "fold"]).any():
        raise ValueError("Augur cross-validation keys must be unique.")
    if set(cross_validation["population"]) != set(populations):
        raise ValueError("Augur cross-validation populations differ from priorities.")
    if bool(((cross_validation["auc"] < 0.0) | (cross_validation["auc"] > 1.0)).any()):
        raise ValueError("Augur cross-validation AUC must lie in [0, 1].")
    expected_cv_keys = {
        (population, subsample, fold)
        for population in populations
        for subsample in range(n_subsamples)
        for fold in range(folds)
    }
    observed_cv_keys = set(
        cross_validation[["population", "subsample", "fold"]].itertuples(index=False, name=None)
    )
    if observed_cv_keys != expected_cv_keys:
        raise ValueError("Augur cross-validation rows do not form the complete population/subsample/fold family.")
    expected_cv_order = sorted(observed_cv_keys)
    if list(cross_validation[["population", "subsample", "fold"]].itertuples(index=False, name=None)) != expected_cv_order:
        raise ValueError("Augur cross-validation rows are not in canonical order.")

    if feature_importance.empty:
        raise ValueError("Augur feature-importance table cannot be empty.")
    if feature_importance.duplicated(["population", "subsample", "fold", "gene"]).any():
        raise ValueError("Augur feature-importance keys must be unique.")
    if set(feature_importance["population"]) != set(populations):
        raise ValueError("Augur feature-importance populations differ from priorities.")
    feature_groups = set(
        feature_importance[["population", "subsample", "fold"]].itertuples(index=False, name=None)
    )
    if feature_groups != expected_cv_keys:
        raise ValueError("Augur feature importance does not cover every cross-validation fold.")
    expected_feature_order = sorted(
        feature_importance[["population", "subsample", "fold", "gene"]].itertuples(index=False, name=None)
    )
    if list(
        feature_importance[["population", "subsample", "fold", "gene"]].itertuples(index=False, name=None)
    ) != expected_feature_order:
        raise ValueError("Augur feature-importance rows are not in canonical order.")
    classifier = parameters.get("classifier")
    if classifier == "random_forest_classifier":
        importance = feature_importance["importance"].to_numpy(dtype=float)
        if bool(((importance < 0.0) | (importance > 1.0)).any()):
            raise ValueError("Random-forest Augur feature importances must lie in [0, 1].")
        sums = feature_importance.groupby(["population", "subsample", "fold"], observed=True)[
            "importance"
        ].sum()
        if not bool(numpy.allclose(sums.to_numpy(dtype=float), 1.0, rtol=1e-7, atol=1e-9)):
            raise ValueError("Random-forest Augur feature importances must sum to one within every fitted fold.")
    elif classifier != "logistic_regression_classifier":
        raise ValueError(f"Augur artifact contains unsupported classifier {classifier!r}.")

    if predictions.empty:
        raise ValueError("Augur prediction-evidence table cannot be empty.")
    prediction_keys = ["population", "subsample", "fold", "observation"]
    if predictions.duplicated(prediction_keys).any():
        raise ValueError("Augur prediction-evidence keys must be unique.")
    if set(predictions["population"]) != set(populations):
        raise ValueError("Augur prediction-evidence populations differ from priorities.")
    if bool((predictions.groupby("observation", observed=True)["population"].nunique() > 1).any()):
        raise ValueError("An Augur prediction observation cannot belong to more than one population.")
    prediction_groups = set(
        predictions[["population", "subsample", "fold"]].itertuples(index=False, name=None)
    )
    if prediction_groups != expected_cv_keys:
        raise ValueError("Augur prediction evidence does not cover every cross-validation fold.")
    expected_prediction_order = sorted(
        predictions[prediction_keys].itertuples(index=False, name=None)
    )
    if list(predictions[prediction_keys].itertuples(index=False, name=None)) != expected_prediction_order:
        raise ValueError("Augur prediction-evidence rows are not in canonical order.")
    if bool((~predictions["true_label"].isin([0, 1])).any()):
        raise ValueError("Augur prediction true labels must be binary integers 0 or 1.")
    prediction_scores = predictions["prediction_score"].to_numpy(dtype=float)
    if bool(((prediction_scores < 0.0) | (prediction_scores > 1.0)).any()):
        raise ValueError("Augur prediction scores must lie in [0, 1].")
    subsample_size = parameters.get("subsample_size")
    if isinstance(subsample_size, bool) or not isinstance(subsample_size, int) or subsample_size < 2:
        raise ValueError("Augur artifact parameters contain an invalid subsample_size.")
    for (population, subsample), group in predictions.groupby(
        ["population", "subsample"], observed=True, sort=False
    ):
        if len(group) != 2 * subsample_size or group["observation"].duplicated().any():
            raise ValueError(
                "Augur prediction evidence must cover each balanced subsample observation exactly once."
            )
        label_counts = group["true_label"].value_counts().to_dict()
        if label_counts != {0: subsample_size, 1: subsample_size}:
            raise ValueError(
                f"Augur prediction evidence is not balanced for population {population!r}, "
                f"subsample {int(subsample)}."
            )
    auc_by_key = {
        key: float(value)
        for key, value in cross_validation.set_index(["population", "subsample", "fold"])["auc"].items()
    }
    for key, group in predictions.groupby(
        ["population", "subsample", "fold"], observed=True, sort=False
    ):
        labels = group["true_label"].to_numpy(dtype=int)
        scores = group["prediction_score"].to_numpy(dtype=float)
        positives = labels == 1
        positive_count = int(positives.sum())
        negative_count = int((~positives).sum())
        if not positive_count or not negative_count:
            raise ValueError("Every Augur prediction fold must contain both true labels.")
        ranks = pandas.Series(scores).rank(method="average").to_numpy(dtype=float)
        observed_auc = float(
            (ranks[positives].sum() - positive_count * (positive_count + 1) / 2)
            / (positive_count * negative_count)
        )
        if not math.isclose(observed_auc, auc_by_key[key], rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError(
                "Augur prediction evidence recomputed AUC does not match cross-validation AUC "
                f"for population/subsample/fold {key!r}."
            )
    return fingerprints


class AugurResult:
    """Immutable-by-interface and tamper-evident canonical Augur result artifact."""

    __slots__ = ("_metadata", "_summary", "_tables")
    artifact_type = AUGUR_ARTIFACT_TYPE

    def __init__(self, *, tables: Mapping[str, Any], summary: Mapping[str, Any], metadata: Mapping[str, Any]) -> None:
        object.__setattr__(self, "_tables", {view: tables[view].copy(deep=True) for view in AUGUR_VIEWS})
        object.__setattr__(self, "_summary", copy.deepcopy(dict(summary)))
        object.__setattr__(self, "_metadata", copy.deepcopy(dict(metadata)))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("AugurResult is immutable; create a new validated artifact instead.")

    @property
    def fingerprint(self) -> str:
        return str(self._metadata["artifact_fingerprint_sha256"])

    @property
    def metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._metadata)

    @property
    def summary(self) -> dict[str, Any]:
        return copy.deepcopy(self._summary)

    def table(self, view: str) -> DataFrame:
        if view not in AUGUR_VIEWS:
            raise ValueError(f"Unsupported Augur result view: {view!r}.")
        return self._tables[view].copy(deep=True)

    def portable_tables(self) -> dict[str, DataFrame]:
        return {view: self._tables[view].copy(deep=True) for view in AUGUR_VIEWS}


def _build_augur_artifact(
    *,
    tables: Mapping[str, Any],
    summary: Mapping[str, Any],
    numpy: Any,
    pandas: Any,
) -> tuple[dict[str, DataFrame], dict[str, Any], dict[str, Any]]:
    summary_copy = copy.deepcopy(dict(summary))
    parameters = summary_copy.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ValueError("Augur summary parameters must be a mapping.")
    table_fingerprints = _validate_canonical_tables(
        tables, parameters=parameters, numpy=numpy, pandas=pandas
    )
    key_results = summary_copy.get("key_results")
    if not isinstance(key_results, dict):
        raise ValueError("Augur summary key_results must be a mapping.")
    key_results["table_fingerprints_sha256"] = dict(table_fingerprints)
    analysis_fingerprint = key_results.get("analysis_fingerprint_sha256")
    if not isinstance(analysis_fingerprint, str) or len(analysis_fingerprint) != 64:
        raise ValueError("Augur summary is missing its analysis fingerprint.")
    result_content_fingerprint = _canonical_json_sha256(
        {
            "schema": AUGUR_PRODUCER_SCHEMA,
            "analysis_fingerprint_sha256": analysis_fingerprint,
            "table_fingerprints_sha256": table_fingerprints,
        }
    )
    key_results["result_content_fingerprint_sha256"] = result_content_fingerprint
    json.dumps(summary_copy, ensure_ascii=False, allow_nan=False)
    summary_fingerprint = _canonical_json_sha256(summary_copy)
    metadata: dict[str, Any] = {
        "schema_version": AUGUR_ARTIFACT_SCHEMA_VERSION,
        "artifact_type": AUGUR_ARTIFACT_TYPE,
        "producer_node_id": AUGUR_PRODUCER_NODE_ID,
        "producer_schema": AUGUR_PRODUCER_SCHEMA,
        "table_fingerprints_sha256": dict(table_fingerprints),
        "result_content_fingerprint_sha256": result_content_fingerprint,
        "summary_fingerprint_sha256": summary_fingerprint,
    }
    metadata["artifact_fingerprint_sha256"] = _canonical_json_sha256(metadata)
    return validate_augur_portable(dict(tables), summary_copy, metadata)


def _build_augur_result(
    *,
    tables: Mapping[str, Any],
    summary: Mapping[str, Any],
    numpy: Any,
    pandas: Any,
) -> AugurResult:
    owned_tables, owned_summary, metadata = _build_augur_artifact(
        tables=tables,
        summary=summary,
        numpy=numpy,
        pandas=pandas,
    )
    result = AugurResult(tables=owned_tables, summary=owned_summary, metadata=metadata)
    validate_augur_result(result)
    return result


def validate_augur_portable(
    tables: Mapping[str, DataFrame],
    summary: Mapping[str, Any],
    metadata: Mapping[str, Any],
) -> tuple[dict[str, DataFrame], dict[str, Any], dict[str, Any]]:
    """Validate a decoded portable Augur table family without constructing a live result object."""

    import numpy as np
    import pandas as pd

    if not isinstance(tables, Mapping) or not isinstance(summary, Mapping) or not isinstance(metadata, Mapping):
        raise TypeError("Augur portable artifact requires table, summary, and metadata mappings.")
    tables = dict(tables)
    summary = dict(summary)
    metadata = dict(metadata)
    expected_metadata = {
        "schema_version",
        "artifact_type",
        "producer_node_id",
        "producer_schema",
        "table_fingerprints_sha256",
        "result_content_fingerprint_sha256",
        "summary_fingerprint_sha256",
        "artifact_fingerprint_sha256",
    }
    if set(metadata) != expected_metadata:
        raise ValueError(
            "Augur artifact metadata schema mismatch; "
            f"missing={sorted(expected_metadata - set(metadata))}, unknown={sorted(set(metadata) - expected_metadata)}."
        )
    if metadata["schema_version"] != AUGUR_ARTIFACT_SCHEMA_VERSION:
        raise ValueError("Augur artifact has an unsupported schema version.")
    if metadata["artifact_type"] != AUGUR_ARTIFACT_TYPE:
        raise ValueError("Augur artifact type identity is invalid.")
    if metadata["producer_node_id"] != AUGUR_PRODUCER_NODE_ID:
        raise ValueError("Augur artifact has the wrong producer node.")
    if metadata["producer_schema"] != AUGUR_PRODUCER_SCHEMA:
        raise ValueError("Augur artifact has an unsupported producer schema.")
    fingerprint = metadata["artifact_fingerprint_sha256"]
    if not (
        isinstance(fingerprint, str)
        and len(fingerprint) == 64
        and all(character in "0123456789abcdef" for character in fingerprint)
    ):
        raise ValueError("Augur artifact fingerprint identity is invalid.")

    expected_summary = {
        "schema_version",
        "node_id",
        "status",
        "methods",
        "results",
        "key_results",
        "parameters",
        "references",
        "software_versions",
        "warnings",
        "limitations",
    }
    if set(summary) != expected_summary:
        raise ValueError("Augur artifact summary schema is invalid.")
    if summary["schema_version"] != AUGUR_SUMMARY_SCHEMA or summary["node_id"] != AUGUR_PRODUCER_NODE_ID:
        raise ValueError("Augur artifact summary producer identity is invalid.")
    if summary["status"] != "exploratory_cell_level_cross_validation":
        raise ValueError("Augur artifact summary has an invalid scientific status.")
    if not isinstance(summary["references"], list) or not summary["references"]:
        raise ValueError("Augur artifact summary must contain method/software references.")
    versions = summary["software_versions"]
    if not isinstance(versions, Mapping) or not _is_compatible_pertpy_1_3(versions.get("pertpy")):
        raise ValueError("Augur artifact summary must report a compatible Pertpy 1.3.x execution.")
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    if metadata["summary_fingerprint_sha256"] != _canonical_json_sha256(summary):
        raise ValueError("Augur artifact summary failed its current-content fingerprint check.")

    parameters = summary["parameters"]
    table_fingerprints = _validate_canonical_tables(
        tables, parameters=parameters, numpy=np, pandas=pd
    )
    if metadata["table_fingerprints_sha256"] != table_fingerprints:
        raise ValueError("Augur artifact tables failed their current-content fingerprint checks.")
    key_results = summary["key_results"]
    if not isinstance(key_results, Mapping):
        raise ValueError("Augur artifact summary key_results must be a mapping.")
    encoding = key_results.get("prediction_label_encoding")
    expected_encoding_keys = {
        "negative_label",
        "negative_condition",
        "positive_label",
        "positive_condition",
    }
    if not isinstance(encoding, Mapping) or set(encoding) != expected_encoding_keys:
        raise ValueError("Augur artifact prediction label encoding is invalid.")
    if encoding["negative_label"] != 0 or encoding["positive_label"] != 1:
        raise ValueError("Augur artifact prediction labels must be encoded as 0 and 1.")
    if not all(
        isinstance(encoding[key], str) and encoding[key]
        for key in ("negative_condition", "positive_condition")
    ):
        raise ValueError("Augur artifact prediction conditions must be nonempty strings.")
    if {encoding["negative_condition"], encoding["positive_condition"]} != {
        parameters.get("control"),
        parameters.get("treatment"),
    }:
        raise ValueError("Augur artifact prediction label encoding differs from its selected Conditions.")
    if key_results.get("table_fingerprints_sha256") != table_fingerprints:
        raise ValueError("Augur summary and artifact table fingerprints disagree.")
    analysis_fingerprint = key_results.get("analysis_fingerprint_sha256")
    result_content_fingerprint = _canonical_json_sha256(
        {
            "schema": AUGUR_PRODUCER_SCHEMA,
            "analysis_fingerprint_sha256": analysis_fingerprint,
            "table_fingerprints_sha256": table_fingerprints,
        }
    )
    if (
        metadata["result_content_fingerprint_sha256"] != result_content_fingerprint
        or key_results.get("result_content_fingerprint_sha256") != result_content_fingerprint
    ):
        raise ValueError("Augur artifact result-content fingerprint is invalid.")
    fingerprint_payload = {key: value for key, value in metadata.items() if key != "artifact_fingerprint_sha256"}
    if fingerprint != _canonical_json_sha256(fingerprint_payload):
        raise ValueError("Augur artifact metadata failed its provenance fingerprint check.")
    return tables, summary, metadata


def validate_augur_result(result: Any) -> tuple[dict[str, DataFrame], dict[str, Any], dict[str, Any]]:
    """Validate class, producer schema, strict summary, tables, and all current-content fingerprints."""

    if type(result) is not AugurResult:
        raise TypeError("Augur Results requires an exact OPENBIO_AUGUR_RESULT artifact from OpenBio Augur.")
    metadata = result.metadata
    if result.fingerprint != metadata.get("artifact_fingerprint_sha256"):
        raise ValueError("Augur artifact fingerprint identity is invalid.")
    return validate_augur_portable(result.portable_tables(), result.summary, metadata)


def _require_pertpy() -> Any:
    try:
        import pertpy
    except (ImportError, OSError) as exc:
        raise RuntimeError("Augur requires a compatible Pertpy 1.3.x installation.") from exc
    return pertpy


def _pertpy_version(pertpy_module: Any) -> str:
    module_version = getattr(pertpy_module, "__version__", None)
    if isinstance(module_version, str) and module_version:
        return module_version
    try:
        return importlib_metadata.version("pertpy")
    except importlib_metadata.PackageNotFoundError as exc:
        raise RuntimeError("Augur could not establish the installed Pertpy distribution version.") from exc


def _is_compatible_pertpy_1_3(version: Any) -> bool:
    return isinstance(version, str) and re.fullmatch(r"1\.3\.\d+(?:\+[A-Za-z0-9.-]+)?", version) is not None


def _validate_pertpy_1_3_interface(pertpy_module: Any) -> tuple[Any, str]:
    version = _pertpy_version(pertpy_module)
    if not _is_compatible_pertpy_1_3(version):
        raise RuntimeError(f"Augur requires a compatible Pertpy 1.3.x release; detected {version!r}.")
    augur_class = getattr(getattr(pertpy_module, "tl", None), "Augur", None)
    if not callable(augur_class):
        raise RuntimeError(f"Pertpy {version} tl.Augur public interface is unavailable.")
    expected = {
        "constructor": (
            "estimator",
            "n_estimators",
            "max_depth",
            "max_features",
            "penalty",
            "random_state",
        ),
        "load": (
            "self",
            "input",
            "layer",
            "meta",
            "label_col",
            "cell_type_col",
            "condition_label",
            "treatment_label",
        ),
        "predict": (
            "self",
            "adata",
            "n_subsamples",
            "subsample_size",
            "folds",
            "min_cells",
            "feature_perc",
            "var_quantile",
            "span",
            "filter_negative_residuals",
            "n_threads",
            "augur_mode",
            "select_variance_features",
            "key_added",
            "random_state",
            "zero_division",
        ),
        "run_cross_validation": (
            "self",
            "subsample",
            "subsample_idx",
            "folds",
            "random_state",
            "zero_division",
        ),
    }
    callables = {
        "constructor": augur_class,
        "load": getattr(augur_class, "load", None),
        "predict": getattr(augur_class, "predict", None),
        "run_cross_validation": getattr(augur_class, "run_cross_validation", None),
    }
    for name, callable_value in callables.items():
        if not callable(callable_value):
            raise RuntimeError(f"Pertpy {version} Augur.{name} is unavailable.")
        try:
            observed = tuple(inspect.signature(callable_value).parameters)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"Could not inspect Pertpy {version} Augur.{name} interface.") from exc
        if observed != expected[name]:
            raise RuntimeError(
                f"Pertpy {version} Augur.{name} interface is incompatible; expected={expected[name]}, observed={observed}."
            )
    return augur_class, version


def _with_prediction_evidence(
    results: Any,
    subsample: Any,
    *,
    folds: int,
    random_state: int | None,
) -> dict[str, Any]:
    import numpy as np
    from sklearn.model_selection import StratifiedKFold

    if not isinstance(results, Mapping):
        raise RuntimeError("Pertpy Augur run_cross_validation returned a non-mapping result.")
    estimators = results.get("estimator")
    if not isinstance(estimators, Sequence) or len(estimators) != folds:
        raise RuntimeError("Pertpy Augur did not retain one fitted estimator per fold.")
    matrix = subsample.X.toarray() if hasattr(subsample.X, "toarray") else np.asarray(subsample.X)
    labels = np.asarray(subsample.obs["y_"])
    splitter = StratifiedKFold(n_splits=folds, random_state=random_state, shuffle=True)
    evidence = []
    for fold, ((_, test_positions), estimator) in enumerate(
        zip(splitter.split(matrix, labels), estimators, strict=True)
    ):
        classes = np.asarray(getattr(estimator, "classes_", []))
        positive_columns = np.flatnonzero(classes == 1)
        if classes.size != 2 or positive_columns.size != 1 or not callable(getattr(estimator, "predict_proba", None)):
            raise RuntimeError("Pertpy Augur fitted classifier cannot provide binary fold probabilities.")
        probabilities = np.asarray(estimator.predict_proba(matrix[test_positions]), dtype=float)
        if probabilities.shape != (len(test_positions), 2):
            raise RuntimeError("Pertpy Augur fitted classifier returned malformed fold probabilities.")
        scores = probabilities[:, int(positive_columns[0])]
        evidence.extend(
            {
                "fold": fold,
                "observation": str(subsample.obs_names[position]),
                "true_label": int(labels[position]),
                "prediction_score": float(score),
            }
            for position, score in zip(test_positions, scores, strict=True)
        )
    owned = dict(results)
    owned["openbio_prediction_evidence"] = evidence
    return owned


def _prediction_retaining_augur_class(base_class: Any) -> Any:
    class PredictionRetainingAugur(base_class):
        def run_cross_validation(
            self: Any,
            subsample: Any,
            *,
            subsample_idx: int,
            folds: int,
            random_state: int | None,
            zero_division: int | str,
        ) -> dict[str, Any]:
            results = super().run_cross_validation(
                subsample,
                subsample_idx=subsample_idx,
                folds=folds,
                random_state=random_state,
                zero_division=zero_division,
            )
            return _with_prediction_evidence(
                results,
                subsample,
                folds=folds,
                random_state=random_state,
            )

    return PredictionRetainingAugur


def _validate_integer(value: Any, *, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"Augur {name} must be an integer >= {minimum}.")
    return value


def _backend_float_column(series: Any, *, label: str, numpy: Any, pandas: Any) -> Any:
    if (
        not pandas.api.types.is_numeric_dtype(series.dtype)
        or pandas.api.types.is_bool_dtype(series.dtype)
    ):
        raise RuntimeError(f"Pertpy Augur {label} must use a numeric non-boolean dtype.")
    try:
        values = series.to_numpy(dtype=float, copy=True)
    except (TypeError, ValueError, OverflowError) as exc:
        raise RuntimeError(f"Pertpy Augur {label} must be numeric.") from exc
    if not bool(numpy.isfinite(values).all()):
        raise RuntimeError(f"Pertpy Augur {label} contains non-finite values.")
    return values


def _backend_integer_column(series: Any, *, label: str, numpy: Any, pandas: Any) -> Any:
    values = _backend_float_column(series, label=label, numpy=numpy, pandas=pandas)
    rounded = numpy.rint(values)
    if not bool(numpy.equal(values, rounded).all()):
        raise RuntimeError(f"Pertpy Augur {label} must contain exact integers.")
    if values.size and (
        float(values.min()) < numpy.iinfo(numpy.int64).min
        or float(values.max()) > numpy.iinfo(numpy.int64).max
    ):
        raise RuntimeError(f"Pertpy Augur {label} exceeds the supported int64 range.")
    return rounded.astype(numpy.int64)


def _canonicalize_backend_results(
    backend_results: Any,
    *,
    eligible_populations: Sequence[str],
    support: Mapping[str, Mapping[str, int]],
    n_subsamples: int,
    subsample_size: int,
    folds: int,
    classifier: str,
    max_result_rows: int,
    max_result_mib: float,
    numpy: Any,
    pandas: Any,
) -> dict[str, DataFrame]:
    if not isinstance(backend_results, Mapping):
        raise RuntimeError("Pertpy Augur returned a non-mapping results object.")
    required_keys = {"summary_metrics", "full_results", "feature_importances", *eligible_populations}
    if set(backend_results) != required_keys:
        raise RuntimeError(
            "Pertpy Augur result family changed; "
            f"missing={sorted(required_keys - set(backend_results))}, unknown={sorted(set(backend_results) - required_keys)}."
        )
    summary_metrics = backend_results["summary_metrics"]
    full_results = backend_results["full_results"]
    feature_importances = backend_results["feature_importances"]
    for value, name in (
        (summary_metrics, "summary_metrics"),
        (full_results, "full_results"),
        (feature_importances, "feature_importances"),
    ):
        if not isinstance(value, pandas.DataFrame):
            raise RuntimeError(f"Pertpy Augur {name} must be a pandas DataFrame in 1.3.0.")
        if not value.columns.is_unique:
            raise RuntimeError(f"Pertpy Augur {name} contains duplicate columns.")

    required_metrics = (
        "mean_augur_score",
        "mean_auc",
        "mean_accuracy",
        "mean_precision",
        "mean_f1",
        "mean_recall",
    )
    if not summary_metrics.index.is_unique or set(summary_metrics.index) != set(required_metrics):
        raise RuntimeError(
            "Pertpy Augur summary metric schema changed; "
            f"expected={sorted(required_metrics)}, observed={sorted(map(str, summary_metrics.index))}."
        )
    if set(summary_metrics.columns) != set(eligible_populations):
        raise RuntimeError("Pertpy Augur summary populations differ from the preflight eligible family.")
    priority = summary_metrics.loc[list(required_metrics), list(eligible_populations)].T.copy()
    priority.index.name = "population"
    priority = priority.reset_index()
    for column in required_metrics:
        priority[column] = _backend_float_column(
            priority[column], label=f"summary metric {column!r}", numpy=numpy, pandas=pandas
        )
    priority["rank"] = priority["mean_augur_score"].rank(method="min", ascending=False).astype(int)
    for column in ("cells_control", "cells_treatment", "samples_control", "samples_treatment"):
        priority[column] = [int(support[population][column]) for population in priority["population"]]
    priority = priority[list(PRIORITY_COLUMNS)].sort_values(
        ["rank", "population"], ascending=[True, True], kind="mergesort", ignore_index=True
    )

    if set(full_results.columns) != {"idx", "augur_score", "folds", "cell_type"}:
        raise RuntimeError("Pertpy Augur full_results schema changed from the audited 1.3.0 interface.")
    cross_validation = full_results.rename(
        columns={"cell_type": "population", "idx": "subsample", "folds": "fold", "augur_score": "auc"}
    )[list(CROSS_VALIDATION_COLUMNS)].copy()
    cross_populations, _ = _canonical_labels(
        cross_validation["population"], label="backend cross-validation population", pandas=pandas
    )
    cross_validation["population"] = cross_populations
    for column in ("subsample", "fold"):
        cross_validation[column] = _backend_integer_column(
            cross_validation[column],
            label=f"full_results column {column!r}",
            numpy=numpy,
            pandas=pandas,
        )
    cross_validation["auc"] = _backend_float_column(
        cross_validation["auc"], label="full_results AUC", numpy=numpy, pandas=pandas
    )
    cross_validation = cross_validation.sort_values(
        ["population", "subsample", "fold"], kind="mergesort", ignore_index=True
    )

    if set(feature_importances.columns) != {
        "genes",
        "feature_importances",
        "subsample_idx",
        "fold",
        "cell_type",
    }:
        raise RuntimeError("Pertpy Augur feature_importances schema changed from the audited 1.3.0 interface.")
    feature_importance = feature_importances.rename(
        columns={
            "cell_type": "population",
            "subsample_idx": "subsample",
            "genes": "gene",
            "feature_importances": "importance",
        }
    )[list(FEATURE_IMPORTANCE_COLUMNS)].copy()
    feature_populations, _ = _canonical_labels(
        feature_importance["population"], label="backend feature-importance population", pandas=pandas
    )
    feature_genes, _ = _canonical_labels(
        feature_importance["gene"], label="backend feature-importance gene", pandas=pandas
    )
    feature_importance["population"] = feature_populations
    feature_importance["gene"] = feature_genes
    for column in ("subsample", "fold"):
        feature_importance[column] = _backend_integer_column(
            feature_importance[column],
            label=f"feature_importances column {column!r}",
            numpy=numpy,
            pandas=pandas,
        )
    feature_importance["importance"] = _backend_float_column(
        feature_importance["importance"],
        label="feature_importances importance",
        numpy=numpy,
        pandas=pandas,
    )
    feature_importance = feature_importance.sort_values(
        ["population", "subsample", "fold", "gene"], kind="mergesort", ignore_index=True
    )

    prediction_rows = []
    for population in eligible_populations:
        population_results = backend_results[population]
        if not isinstance(population_results, Sequence) or isinstance(population_results, str | bytes):
            raise RuntimeError(f"Pertpy Augur fold results for population {population!r} must be a sequence.")
        if len(population_results) != n_subsamples:
            raise RuntimeError(
                f"Pertpy Augur fold results for population {population!r} do not cover every subsample."
            )
        for subsample, fold_results in enumerate(population_results):
            if not isinstance(fold_results, Mapping):
                raise RuntimeError("Pertpy Augur fold result must be a mapping.")
            evidence = fold_results.get("openbio_prediction_evidence")
            if not isinstance(evidence, Sequence) or isinstance(evidence, str | bytes):
                raise RuntimeError("Pertpy Augur did not retain fold-level prediction evidence.")
            for row in evidence:
                if not isinstance(row, Mapping) or set(row) != {
                    "fold",
                    "observation",
                    "true_label",
                    "prediction_score",
                }:
                    raise RuntimeError("Pertpy Augur prediction-evidence schema changed.")
                prediction_rows.append({"population": population, "subsample": subsample, **row})
    predictions = pandas.DataFrame(prediction_rows, columns=PREDICTION_COLUMNS)
    prediction_populations, _ = _canonical_labels(
        predictions["population"], label="backend prediction population", pandas=pandas
    )
    prediction_observations, _ = _canonical_labels(
        predictions["observation"], label="backend prediction observation", pandas=pandas
    )
    predictions["population"] = prediction_populations
    predictions["observation"] = prediction_observations
    for column in ("subsample", "fold", "true_label"):
        predictions[column] = _backend_integer_column(
            predictions[column], label=f"prediction-evidence column {column!r}", numpy=numpy, pandas=pandas
        )
    predictions["prediction_score"] = _backend_float_column(
        predictions["prediction_score"],
        label="prediction-evidence prediction_score",
        numpy=numpy,
        pandas=pandas,
    )
    predictions = predictions.sort_values(
        ["population", "subsample", "fold", "observation"], kind="mergesort", ignore_index=True
    )

    tables = {
        "priorities": priority,
        "cross_validation": cross_validation,
        "feature_importance": feature_importance,
        "predictions": predictions,
    }
    total_rows = sum(len(table) for table in tables.values())
    if total_rows > max_result_rows:
        raise ValueError(
            f"Augur canonical results contain {total_rows:,} rows, exceeding max_result_rows={max_result_rows:,}."
        )
    memory_bytes = sum(int(table.memory_usage(index=True, deep=True).sum()) for table in tables.values())
    budget_bytes = int(max_result_mib * 1024 * 1024)
    if memory_bytes > budget_bytes:
        raise MemoryError(
            f"Augur canonical results use {memory_bytes / (1024 * 1024):,.1f} MiB, "
            f"exceeding max_result_mib={max_result_mib:,.1f}."
        )
    _validate_canonical_tables(
        tables,
        parameters={
            "n_subsamples": n_subsamples,
            "subsample_size": subsample_size,
            "folds": folds,
            "classifier": classifier,
        },
        numpy=numpy,
        pandas=pandas,
    )
    observed_auc_means = cross_validation.groupby("population", observed=True)["auc"].mean()
    for record in priority.itertuples(index=False):
        expected_auc = float(observed_auc_means.loc[record.population])
        if not (
            math.isclose(float(record.mean_auc), expected_auc, rel_tol=1e-10, abs_tol=1e-12)
            and math.isclose(float(record.mean_augur_score), expected_auc, rel_tol=1e-10, abs_tol=1e-12)
        ):
            raise RuntimeError(
                "Pertpy Augur summary mean AUC/Augur score does not equal the complete full_results fold mean "
                f"for population {record.population!r}."
            )
    return tables


def _run_augur_artifact_owned(
    adata: AnnData,
    *,
    sample_key: str = "sample",
    population_key: str = "cell_type",
    condition_key: str = "condition",
    control: str,
    treatment: str,
    classifier: str = "random_forest_classifier",
    source_kind: str = "raw",
    layer_name: str | None = None,
    annotation_status: str = "unknown",
    technical_batch_key: str = "",
    n_subsamples: int = 50,
    subsample_size: int = 20,
    folds: int = 3,
    n_threads: int = 1,
    random_seed: int = 123,
    max_result_rows: int = 10_000_000,
    max_result_mib: float = 1024.0,
    pertpy_module: Any | None = None,
) -> tuple[dict[str, DataFrame], dict[str, Any], dict[str, Any]]:
    """Run one audited two-condition classifier prioritization into an owned portable payload."""

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from scipy import sparse

    if not isinstance(adata, AnnData):
        raise TypeError("Augur input must be an AnnData object.")
    if getattr(adata, "isbacked", False):
        raise ValueError("Augur requires an in-memory AnnData; call to_memory() first.")
    if adata.n_obs < 1 or adata.n_vars < 1:
        raise ValueError("Augur requires at least one cell and one current feature.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Augur requires unique current cell and feature identifiers.")
    if not all(isinstance(value, str) and value and value == value.strip() for value in adata.obs_names):
        raise ValueError("Augur requires canonical nonempty string cell identifiers.")
    sample_key = _strict_text(sample_key, label="Sample key")
    population_key = _strict_text(population_key, label="population key")
    condition_key = _strict_text(condition_key, label="Condition key")
    technical_batch_key = _strict_text(
        technical_batch_key, label="Technical batch key", allow_empty=True
    )
    control = _strict_text(control, label="control Condition")
    treatment = _strict_text(treatment, label="treatment Condition")
    if control == treatment:
        raise ValueError("Augur control and treatment Conditions must differ.")
    role_keys = [sample_key, population_key, condition_key]
    if technical_batch_key:
        role_keys.append(technical_batch_key)
    sample_role_reuse = [
        role
        for role, key in (
            ("population", population_key),
            ("Condition", condition_key),
            ("Technical batch", technical_batch_key),
        )
        if key and sample_key == key
    ]
    if population_key == condition_key:
        raise ValueError("Augur population and Condition keys must differ so each population can contain both arms.")
    missing_columns = [key for key in dict.fromkeys(role_keys) if key not in adata.obs]
    if missing_columns:
        raise ValueError(f"Augur observation columns were not found: {missing_columns}.")
    if classifier not in AUGUR_CLASSIFIERS:
        raise ValueError(f"Augur supports classifier-only choices {list(AUGUR_CLASSIFIERS)}; observed {classifier!r}.")
    if annotation_status not in {"unknown", "provisional", "curated"}:
        raise ValueError(f"Unsupported Augur annotation status: {annotation_status!r}.")
    n_subsamples = _validate_integer(n_subsamples, name="n_subsamples", minimum=1)
    subsample_size = _validate_integer(subsample_size, name="subsample_size", minimum=2)
    folds = _validate_integer(folds, name="folds", minimum=2)
    n_threads = _validate_integer(n_threads, name="n_threads", minimum=1)
    random_seed = _validate_integer(random_seed, name="random_seed", minimum=1)
    max_result_rows = _validate_integer(max_result_rows, name="max_result_rows", minimum=1)
    if folds > subsample_size:
        raise ValueError("Augur folds cannot exceed per-Condition subsample_size.")
    if isinstance(max_result_mib, bool) or not isinstance(max_result_mib, int | float):
        raise TypeError("Augur max_result_mib must be numeric.")
    max_result_mib = float(max_result_mib)
    if not math.isfinite(max_result_mib) or max_result_mib <= 0:
        raise ValueError("Augur max_result_mib must be finite and > 0.")

    if source_kind not in {"raw", "X", "layer"}:
        raise ValueError(f"Unsupported Augur count source: {source_kind!r}.")
    if source_kind == "raw":
        if layer_name is not None:
            raise ValueError("Augur Raw source cannot carry a layer name.")
        if adata.raw is None:
            raise ValueError("Augur selected Raw, but this AnnData has no Raw snapshot.")
        if not adata.raw.obs_names.equals(adata.obs_names):
            raise ValueError("Augur Raw snapshot observations are not aligned to current cells.")
        source_matrix = adata.raw.X
        source_var = adata.raw.var.copy()
        source_label = "raw.X"
    elif source_kind == "X":
        if layer_name is not None:
            raise ValueError("Augur X source cannot carry a layer name.")
        source_matrix = adata.X
        source_var = adata.var.copy()
        source_label = "X"
    else:
        layer_name = _strict_text(layer_name, label="count layer")
        if layer_name not in adata.layers:
            raise ValueError(f"Augur count layer not found: {layer_name!r}.")
        source_matrix = adata.layers[layer_name]
        source_var = adata.var.copy()
        source_label = f"layers[{layer_name!r}]"
    if not source_var.index.is_unique or not all(
        isinstance(value, str) and value and value == value.strip() for value in source_var.index
    ):
        raise ValueError("Augur selected source requires unique canonical feature identifiers.")
    values = np.asarray(source_matrix.data) if sparse.issparse(source_matrix) else np.asarray(source_matrix).ravel()
    if np.issubdtype(values.dtype, np.bool_) or not np.issubdtype(values.dtype, np.number):
        raise TypeError("Augur selected expression source must be numeric and non-boolean.")
    numeric = values.astype(float, copy=False)
    if numeric.size and not bool(np.isfinite(numeric).all()):
        raise ValueError("Augur selected expression source contains non-finite values.")
    if numeric.size and bool((numeric < 0).any()):
        raise ValueError("Augur selected expression source contains negative values.")
    if not numeric.size or not bool((numeric > 0).any()):
        raise ValueError("Augur selected expression source contains no positive values.")
    integer_like = bool(np.allclose(numeric, np.rint(numeric), rtol=0.0, atol=1e-8))
    if integer_like and float(numeric.max(initial=0.0)) > np.iinfo(np.int64).max:
        raise OverflowError("Augur Raw counts exceed the supported int64 range.")
    target_dtype = np.int64 if integer_like else np.float64
    counts = (
        source_matrix.astype(target_dtype).tocsr()
        if sparse.issparse(source_matrix)
        else np.asarray(source_matrix, dtype=target_dtype)
    )

    sample_labels, sample_order = _canonical_labels(adata.obs[sample_key], label="Sample labels", pandas=pd)
    population_labels, population_order = _canonical_labels(
        adata.obs[population_key], label="population labels", pandas=pd
    )
    condition_labels, condition_order = _canonical_labels(
        adata.obs[condition_key], label="Condition labels", pandas=pd
    )
    technical_labels = None
    technical_order: list[str] = []
    if technical_batch_key:
        technical_labels, technical_order = _canonical_labels(
            adata.obs[technical_batch_key], label="Technical batch labels", pandas=pd
        )
    if control not in condition_order or treatment not in condition_order:
        raise ValueError(
            f"Augur requested Conditions were not both observed; available={condition_order}, "
            f"control={control!r}, treatment={treatment!r}."
        )

    sample_conditions: dict[str, set[str]] = {}
    sample_batches: dict[str, set[str]] = {}
    sample_cells: dict[str, int] = {}
    for index, sample in enumerate(sample_labels):
        condition = condition_labels[index]
        sample_conditions.setdefault(sample, set()).add(condition)
        sample_cells[sample] = sample_cells.get(sample, 0) + 1
        if technical_labels is not None:
            sample_batches.setdefault(sample, set()).add(technical_labels[index])

    selected_positions = [
        index for index, condition in enumerate(condition_labels) if condition in {control, treatment}
    ]
    selected_samples_by_condition = {
        condition: [
            sample
            for sample in sample_order
            if any(
                sample_labels[index] == sample and condition_labels[index] == condition
                for index in selected_positions
            )
        ]
        for condition in (control, treatment)
    }
    sample_condition_mapping_valid = all(len(values) == 1 for values in sample_conditions.values())
    sample_batch_mapping_valid = all(len(values) == 1 for values in sample_batches.values())

    batch_audit: dict[str, Any]
    selected_batch_order: list[str] = []
    if technical_labels is None:
        batch_audit = {
            "technical_batch_key": None,
            "status": "not_provided",
            "sample_batch_mapping_valid": None,
            "perfect_condition_confounding": None,
            "condition_batch_imbalanced": None,
            "condition_batch_sample_counts": {},
            "batch_condition_records": [],
            "population_batch_records": [],
        }
    else:
        batch_condition_records = []
        batch_conditions: dict[str, set[str]] = {}
        for batch in technical_order:
            batch_positions = [
                index for index in selected_positions if technical_labels[index] == batch
            ]
            conditions = [
                condition
                for condition in (control, treatment)
                if any(condition_labels[index] == condition for index in batch_positions)
            ]
            if batch_positions:
                selected_batch_order.append(batch)
                batch_conditions[batch] = set(conditions)
                batch_condition_records.append(
                    {
                        "technical_batch": batch,
                        "conditions": conditions,
                        "sample_counts": {
                            condition: len(
                                {
                                    sample_labels[index]
                                    for index in batch_positions
                                    if condition_labels[index] == condition
                                }
                            )
                            for condition in (control, treatment)
                        },
                    }
                )
        perfect_confounding = bool(batch_conditions) and all(
            len(conditions) == 1 for conditions in batch_conditions.values()
        )
        condition_batch_sample_counts = {
            condition: {
                batch: len(
                    {
                        sample_labels[index]
                        for index in selected_positions
                        if condition_labels[index] == condition and technical_labels[index] == batch
                    }
                )
                for batch in batch_conditions
            }
            for condition in (control, treatment)
        }
        condition_batch_imbalanced = any(
            condition_batch_sample_counts[control][batch]
            != condition_batch_sample_counts[treatment][batch]
            for batch in batch_conditions
        )
        batch_audit = {
            "technical_batch_key": technical_batch_key,
            "status": "audited",
            "sample_batch_mapping_valid": sample_batch_mapping_valid,
            "perfect_condition_confounding": perfect_confounding,
            "condition_batch_imbalanced": condition_batch_imbalanced,
            "condition_batch_sample_counts": condition_batch_sample_counts,
            "batch_condition_records": batch_condition_records,
            "population_batch_records": [],
        }

    support: dict[str, dict[str, int]] = {}
    skipped: list[dict[str, Any]] = []
    eligible: list[str] = []
    for population in population_order:
        population_positions = [
            index for index in selected_positions if population_labels[index] == population
        ]
        arm_positions = {
            condition: [index for index in population_positions if condition_labels[index] == condition]
            for condition in (control, treatment)
        }
        arm_samples = {
            condition: list(dict.fromkeys(sample_labels[index] for index in arm_positions[condition]))
            for condition in (control, treatment)
        }
        record = {
            "cells_control": len(arm_positions[control]),
            "cells_treatment": len(arm_positions[treatment]),
            "samples_control": len(arm_samples[control]),
            "samples_treatment": len(arm_samples[treatment]),
        }
        support[population] = record
        reasons = []
        if len(population_positions) < n_subsamples:
            reasons.append("below_pertpy_min_cells_n_subsamples")
        if len(arm_positions[control]) < subsample_size:
            reasons.append("control_below_subsample_size")
        if len(arm_positions[treatment]) < subsample_size:
            reasons.append("treatment_below_subsample_size")
        if technical_labels is not None:
            population_batch_counts = {
                condition: {
                    batch: len(
                        {
                            sample_labels[index]
                            for index in arm_positions[condition]
                            if technical_labels[index] == batch
                        }
                    )
                    for batch in selected_batch_order
                }
                for condition in (control, treatment)
            }
            population_shared_batches = [
                batch
                for batch in selected_batch_order
                if population_batch_counts[control][batch] > 0
                and population_batch_counts[treatment][batch] > 0
            ]
            population_perfect_confounding = not population_shared_batches
            population_batch_imbalanced = any(
                population_batch_counts[control][batch]
                != population_batch_counts[treatment][batch]
                for batch in selected_batch_order
            )
            population_batch_record = {
                "population": population,
                "condition_batch_sample_counts": population_batch_counts,
                "shared_batches": population_shared_batches,
                "perfect_condition_confounding": population_perfect_confounding,
                "condition_batch_imbalanced": population_batch_imbalanced,
                "eligible_by_support": not reasons,
            }
            batch_audit["population_batch_records"].append(population_batch_record)
        if reasons:
            skipped.append({"population": population, "reasons": reasons, **record})
        else:
            eligible.append(population)
    if not eligible:
        raise ValueError(
            "Augur found no eligible population after the required cell/subsample support checks."
        )

    expected_cross_rows = len(eligible) * n_subsamples * folds
    estimated_features_per_fold = max(1, math.ceil(int(source_var.shape[0]) * 0.25))
    estimated_feature_rows = expected_cross_rows * estimated_features_per_fold
    estimated_prediction_rows = len(eligible) * n_subsamples * 2 * subsample_size
    estimated_total_rows = len(eligible) + expected_cross_rows + estimated_feature_rows + estimated_prediction_rows
    if estimated_total_rows > max_result_rows:
        raise ValueError(
            f"Augur estimated result size is {estimated_total_rows:,} rows, exceeding "
            f"max_result_rows={max_result_rows:,}; reduce populations/subsamples/folds or raise the guard."
        )
    estimated_bytes = (
        len(eligible) * 512
        + expected_cross_rows * 96
        + estimated_feature_rows * 160
        + estimated_prediction_rows * 192
    )
    if estimated_bytes > int(max_result_mib * 1024 * 1024):
        raise MemoryError(
            f"Augur estimated canonical result memory is {estimated_bytes / (1024 * 1024):,.1f} MiB, "
            f"exceeding max_result_mib={max_result_mib:,.1f}."
        )

    analysis_positions = [
        index for index in selected_positions if population_labels[index] in set(eligible)
    ]
    # Pertpy owns and mutates this scientifically selected cell subset during classifier preparation.
    analysis_counts = counts[analysis_positions].copy()
    analysis_obs_names = [str(adata.obs_names[index]) for index in analysis_positions]
    analysis_samples = [sample_labels[index] for index in analysis_positions]
    analysis_populations = [population_labels[index] for index in analysis_positions]
    analysis_conditions = [condition_labels[index] for index in analysis_positions]
    analysis_batches = (
        [technical_labels[index] for index in analysis_positions] if technical_labels is not None else None
    )
    selected_source_fingerprint = _count_matrix_fingerprint(
        counts,
        observation_names=[str(value) for value in adata.obs_names],
        feature_names=[str(value) for value in source_var.index],
        sparse=sparse,
        numpy=np,
    )
    analysis_count_fingerprint = _count_matrix_fingerprint(
        analysis_counts,
        observation_names=analysis_obs_names,
        feature_names=[str(value) for value in source_var.index],
        sparse=sparse,
        numpy=np,
    )
    analysis_observation_fingerprint = _observation_fingerprint(
        observation_names=analysis_obs_names,
        samples=analysis_samples,
        populations=analysis_populations,
        conditions=analysis_conditions,
        technical_batches=analysis_batches,
    )
    work_obs_data = {
        population_key: pd.Categorical(analysis_populations, categories=eligible, ordered=True),
        condition_key: pd.Categorical(
            analysis_conditions, categories=[control, treatment], ordered=True
        ),
    }
    if sample_key not in work_obs_data:
        work_obs_data[sample_key] = analysis_samples
    if (
        technical_batch_key
        and analysis_batches is not None
        and technical_batch_key not in work_obs_data
    ):
        work_obs_data[technical_batch_key] = analysis_batches
    work_obs = pd.DataFrame(work_obs_data, index=analysis_obs_names)
    work = AnnData(X=analysis_counts, obs=work_obs, var=source_var.copy())

    parameters = {
        "sample_key": sample_key,
        "population_key": population_key,
        "condition_key": condition_key,
        "control": control,
        "treatment": treatment,
        "classifier": classifier,
        "source": source_kind,
        "layer_name": layer_name,
        "source_label": source_label,
        "annotation_status": annotation_status,
        "technical_batch_key": technical_batch_key or None,
        "n_subsamples": n_subsamples,
        "subsample_size": subsample_size,
        "folds": folds,
        "n_threads": n_threads,
        "random_seed": random_seed,
        "max_result_rows": max_result_rows,
        "max_result_mib": max_result_mib,
        "min_cells": None,
        "resolved_min_cells": n_subsamples,
        "feature_perc": 0.5,
        "var_quantile": 0.5,
        "span": 0.75,
        "filter_negative_residuals": False,
        "augur_mode": "default",
        "select_variance_features": True,
        "key_added": "openbio_augur",
        "zero_division": 0,
        "estimator_hyperparameters": "Pertpy 1.3.0 defaults",
    }
    analysis_fingerprint = _canonical_json_sha256(
        {
            "schema": "openbio-singlecell/augur-analysis/v1",
            "parameters": parameters,
            "selected_source_fingerprint_sha256": selected_source_fingerprint,
            "analysis_count_fingerprint_sha256": analysis_count_fingerprint,
            "analysis_observation_fingerprint_sha256": analysis_observation_fingerprint,
            "eligible_populations": eligible,
            "skipped_populations": skipped,
        }
    )

    pertpy_module = _require_pertpy() if pertpy_module is None else pertpy_module
    augur_class, pertpy_version = _validate_pertpy_1_3_interface(pertpy_module)
    try:
        augur = _prediction_retaining_augur_class(augur_class)(classifier, random_state=random_seed)
    except Exception as exc:
        raise RuntimeError(f"Pertpy 1.3.0 Augur initialization failed ({type(exc).__name__}: {exc}).") from exc
    with python_warnings.catch_warnings(record=True) as caught:
        python_warnings.simplefilter("always")
        try:
            loaded = augur.load(
                work,
                label_col=condition_key,
                cell_type_col=population_key,
                condition_label=control,
                treatment_label=treatment,
            )
            if not isinstance(loaded, AnnData) or not {"label", "y_"}.issubset(loaded.obs):
                raise RuntimeError("Pertpy Augur load did not expose canonical binary label encoding.")
            encoded_conditions: dict[int, str] = {}
            for condition, encoded in zip(loaded.obs["label"], loaded.obs["y_"], strict=True):
                if isinstance(encoded, bool) or not isinstance(encoded, int | np.integer) or int(encoded) not in {0, 1}:
                    raise RuntimeError("Pertpy Augur load returned a non-binary label encoding.")
                code = int(encoded)
                condition = str(condition)
                if code in encoded_conditions and encoded_conditions[code] != condition:
                    raise RuntimeError("Pertpy Augur load returned an inconsistent label encoding.")
                encoded_conditions[code] = condition
            if set(encoded_conditions) != {0, 1} or set(encoded_conditions.values()) != {control, treatment}:
                raise RuntimeError("Pertpy Augur load label encoding differs from the selected Conditions.")
            prediction_label_encoding = {
                "negative_label": 0,
                "negative_condition": encoded_conditions[0],
                "positive_label": 1,
                "positive_condition": encoded_conditions[1],
            }
            updated, backend_results = augur.predict(
                loaded,
                n_subsamples=n_subsamples,
                subsample_size=subsample_size,
                folds=folds,
                min_cells=None,
                feature_perc=0.5,
                var_quantile=0.5,
                span=0.75,
                filter_negative_residuals=False,
                n_threads=n_threads,
                augur_mode="default",
                select_variance_features=True,
                key_added="openbio_augur",
                random_state=random_seed,
                zero_division=0,
            )
        except Exception as exc:
            raise RuntimeError(f"Pertpy 1.3.0 Augur execution failed ({type(exc).__name__}: {exc}).") from exc
    if not isinstance(updated, AnnData):
        raise RuntimeError("Pertpy 1.3.0 Augur returned a non-AnnData work result.")
    tables = _canonicalize_backend_results(
        backend_results,
        eligible_populations=eligible,
        support=support,
        n_subsamples=n_subsamples,
        subsample_size=subsample_size,
        folds=folds,
        classifier=classifier,
        max_result_rows=max_result_rows,
        max_result_mib=max_result_mib,
        numpy=np,
        pandas=pd,
    )

    priorities = tables["priorities"]
    cross_validation_auc_distribution = []
    for population in eligible:
        auc_values = tables["cross_validation"].loc[
            tables["cross_validation"]["population"] == population, "auc"
        ].to_numpy(dtype=float)
        cross_validation_auc_distribution.append(
            {
                "population": population,
                "observations": int(auc_values.size),
                "mean": float(auc_values.mean()),
                "standard_deviation": float(auc_values.std(ddof=1)),
                "minimum": float(auc_values.min()),
                "q025": float(np.quantile(auc_values, 0.025)),
                "q25": float(np.quantile(auc_values, 0.25)),
                "median": float(np.median(auc_values)),
                "q75": float(np.quantile(auc_values, 0.75)),
                "q975": float(np.quantile(auc_values, 0.975)),
                "maximum": float(auc_values.max()),
            }
        )
    top_rank = int(priorities["rank"].min())
    top_records = priorities.loc[
        priorities["rank"] == top_rank, ["population", "rank", "mean_augur_score", "mean_auc"]
    ].to_dict(orient="records")
    tie_records = []
    for score, group in priorities.groupby("mean_augur_score", sort=False):
        if len(group) > 1:
            tie_records.append(
                {
                    "mean_augur_score": float(score),
                    "rank": int(group["rank"].iloc[0]),
                    "populations": sorted(group["population"].tolist()),
                }
            )
    sample_records = [
        {
            "sample": sample,
            "conditions": [
                condition for condition in condition_order if condition in sample_conditions[sample]
            ],
            "technical_batches": [
                batch for batch in technical_order if batch in sample_batches.get(sample, set())
            ],
            "cells": sample_cells[sample],
            "selected": bool(sample_conditions[sample] & {control, treatment}),
        }
        for sample in sample_order
    ]
    minimum_two_samples_per_arm = all(
        len(selected_samples_by_condition[condition]) >= 2
        for condition in (control, treatment)
    )
    sample_audit = {
        "sample_key": sample_key,
        "identity_status": (
            "role_reused_unverified" if sample_role_reuse else "distinct_column_user_declared"
        ),
        "reused_roles": sample_role_reuse,
        "sample_condition_mapping_valid": sample_condition_mapping_valid,
        "minimum_two_samples_per_arm": minimum_two_samples_per_arm,
        "selected_sample_counts": {
            control: len(selected_samples_by_condition[control]),
            treatment: len(selected_samples_by_condition[treatment]),
        },
        "sample_records": sample_records,
    }
    report_warnings = [
        "Augur 1.3.0 cross-validation splits cells, not biological Samples; results are exploratory prioritization and not Sample-level Condition inference.",
        "Pertpy 1.3.0 seeds cell subsampling by subsample index; random_seed controls estimator/fold randomness but does not change that deterministic subsample sequence.",
        "Feature importance is classifier- and sampled-feature-dependent; it is not causal evidence or a differential-expression effect.",
    ]
    report_warnings.extend(f"Pertpy warning: {warning.message}" for warning in caught)
    if not integer_like:
        report_warnings.append(
            "The explicitly selected expression source is nonnegative but not integer-like; it was preserved "
            "without coercion. Confirm that this input is appropriate for Pertpy Augur normalization."
        )
    if not sample_condition_mapping_valid:
        report_warnings.append(
            "At least one Sample appears in more than one Condition. Augur remains a cell-level prioritization; "
            "the repeated-measure mapping is disclosed but not modeled."
        )
    if sample_role_reuse:
        report_warnings.append(
            f"sample_key reuses the {sample_role_reuse} role(s); Sample identity is unverified and the reported "
            "unique-label support must not be interpreted as biological replication."
        )
    if not minimum_two_samples_per_arm:
        report_warnings.append(
            "At least one selected Condition has fewer than two biological Samples; the cell-level Augur result "
            "is not independently replicated."
        )
    if technical_labels is None:
        report_warnings.append(
            "No Technical batch column was supplied; perfect Condition/batch confounding could not be audited."
        )
    else:
        if not sample_batch_mapping_valid:
            report_warnings.append(
                "At least one Sample spans multiple Technical batches; the mapping is disclosed but not modeled."
            )
        if batch_audit["perfect_condition_confounding"]:
            report_warnings.append(
                "Condition is perfectly confounded with Technical batch in the selected cells; Augur cannot "
                "separate those sources of classifier performance."
            )
        elif batch_audit["condition_batch_imbalanced"]:
            report_warnings.append(
                "Biological Sample counts across Technical batches differ between Conditions; Augur does not model "
                "this partial batch imbalance."
            )
    if technical_labels is not None and any(
        record["eligible_by_support"] and record["condition_batch_imbalanced"]
        for record in batch_audit["population_batch_records"]
    ):
        report_warnings.append(
            "At least one eligible population has unequal biological Sample support across Technical batches; "
            "population priorities can reflect this partial batch imbalance."
        )
    if technical_labels is not None and any(
        record["eligible_by_support"] and record["perfect_condition_confounding"]
        for record in batch_audit["population_batch_records"]
    ):
        report_warnings.append(
            "At least one eligible population has perfect Condition/Technical-batch confounding; its priority "
            "cannot distinguish biological Condition from Technical batch."
        )
    if annotation_status == "unknown":
        report_warnings.append(
            "Population annotation status is unknown; population labels must not be presented as curated identities."
        )
    elif annotation_status == "provisional":
        report_warnings.append("Population labels are provisional and require independent evidence review.")
    if skipped:
        report_warnings.append(
            f"{len(skipped)} population(s) were excluded by declared cell/Sample support checks before Pertpy."
        )
    if len(selected_samples_by_condition[control]) != len(selected_samples_by_condition[treatment]):
        report_warnings.append("Selected Conditions contain unequal numbers of biological Samples.")
    results_text = (
        f"Augur produced exploratory cell-level cross-validation priorities for {len(eligible):,} population(s); "
        f"the highest rank contains {len(top_records):,} population(s) with mean classifier AUC "
        f"{float(priorities['mean_augur_score'].max()):.3f}."
    )
    versions = collect_software_versions(
        ("pertpy", "scikit-learn", "anndata", "numpy", "pandas", "scipy", "scikit-misc")
    )
    versions["pertpy"] = pertpy_version
    summary = {
        "schema_version": AUGUR_SUMMARY_SCHEMA,
        "node_id": AUGUR_PRODUCER_NODE_ID,
        "status": "exploratory_cell_level_cross_validation",
        "methods": (
            f"Using explicitly selected {source_label} expression, cells from two explicit Conditions were audited "
            "by caller-declared Sample labels, population, and optional Technical batch. These audits were "
            "descriptive and did "
            "not gate the cell-level calculation. Populations with sufficient Pertpy cell support were "
            f"prioritized with Pertpy {pertpy_version} Augur {classifier}, {n_subsamples} balanced cell "
            f"subsamples of {subsample_size} cells per Condition, and {folds}-fold cell-level cross-validation. "
            "True labels and positive-class probabilities from each fitted test fold were retained as Diagnostic "
            "evidence and checked against the reported fold AUC. "
            "Original Augur variance-feature selection used feature_perc=0.5, var_quantile=0.5, span=0.75, "
            "and filter_negative_residuals=False."
        ),
        "results": results_text,
        "key_results": {
            "input_cells": int(adata.n_obs),
            "input_current_features": int(adata.n_vars),
            "selected_source_features": int(source_var.shape[0]),
            "raw_snapshot_features": int(adata.raw.n_vars) if source_kind == "raw" else None,
            "selected_condition_cells": len(selected_positions),
            "analyzed_cells": len(analysis_positions),
            "condition_cell_counts": {
                control: sum(condition_labels[index] == control for index in selected_positions),
                treatment: sum(condition_labels[index] == treatment for index in selected_positions),
            },
            "eligible_populations": eligible,
            "eligible_population_count": len(eligible),
            "skipped_populations": skipped,
            "priority_results": priorities[
                ["rank", "population", "mean_augur_score", "mean_auc", "mean_accuracy", "mean_f1"]
            ].to_dict(orient="records"),
            "top_priority": top_records,
            "priority_ties": tie_records,
            "cross_validation_rows": len(tables["cross_validation"]),
            "cross_validation_auc_distribution": cross_validation_auc_distribution,
            "feature_importance_rows": len(tables["feature_importance"]),
            "prediction_evidence_rows": len(tables["predictions"]),
            "prediction_label_encoding": prediction_label_encoding,
            "sample_audit": sample_audit,
            "technical_batch_audit": batch_audit,
            "annotation_status": annotation_status,
            "count_source": {
                "source": source_label,
                "selection": "explicit",
                "integer_like": integer_like,
            },
            "selected_source_fingerprint_sha256": selected_source_fingerprint,
            "raw_snapshot_fingerprint_sha256": (
                selected_source_fingerprint if source_kind == "raw" else None
            ),
            "analysis_count_fingerprint_sha256": analysis_count_fingerprint,
            "analysis_observation_fingerprint_sha256": analysis_observation_fingerprint,
            "analysis_fingerprint_sha256": analysis_fingerprint,
        },
        "parameters": parameters,
        "references": copy.deepcopy(AUGUR_REFERENCES),
        "software_versions": versions,
        "warnings": report_warnings,
        "limitations": [
            "Cross-validation is cell-level and allows cells from one biological Sample to occur in training and test folds.",
            "Augur AUC measures label separability/responsiveness priority, not effect direction, a p-value, or statistical significance.",
            "Sample support and batch audits improve reviewability but do not make this a replicate-aware Condition analysis.",
            "Cell abundance, library properties, donor composition, population labeling, and Technical batch can influence classifier performance.",
            "Priorities and feature importance depend on the classifier, random feature subsets, fold assignment, and Pertpy 1.3.0 implementation.",
            "Prediction rows are out-of-fold within each subsample; one cell can appear again in another subsample.",
            "The selected expression matrix is structurally validated; source choice remains part of the workflow.",
        ],
    }
    return _build_augur_artifact(tables=tables, summary=summary, numpy=np, pandas=pd)


def run_augur_artifact(
    adata: AnnData,
    **kwargs: Any,
) -> tuple[dict[str, DataFrame], dict[str, Any], dict[str, Any]]:
    """Return worker-owned canonical tables, strict summary, and artifact metadata."""

    return _run_augur_artifact_owned(adata, **kwargs)


def run_augur_analysis(adata: AnnData, **kwargs: Any) -> AugurResult:
    """Build the legacy in-process value for direct Python consumers."""

    tables, summary, metadata = _run_augur_artifact_owned(adata, **kwargs)
    result = AugurResult(tables=tables, summary=summary, metadata=metadata)
    validate_augur_result(result)
    return result


def run_augur_portable(adata: AnnData, **kwargs: Any) -> tuple[dict[str, DataFrame], dict[str, Any]]:
    """Portable equivalent returning the owned canonical DataFrames plus the strict summary."""

    tables, summary, _ = _run_augur_artifact_owned(adata, **kwargs)
    return tables, summary


def select_augur_view(result: AugurResult, *, view: str) -> tuple[DataFrame, dict[str, Any]]:
    if view not in AUGUR_VIEWS:
        raise ValueError(f"Unsupported Augur result view: {view!r}; expected one of {list(AUGUR_VIEWS)}.")
    tables, parent_summary, metadata = validate_augur_result(result)
    selected = tables[view].copy(deep=True)
    table_fingerprint = metadata["table_fingerprints_sha256"][view]
    summary = copy.deepcopy(parent_summary)
    summary["selected_view"] = {
        "node_id": AUGUR_RESULTS_NODE_ID,
        "view": view,
        "rows": int(selected.shape[0]),
        "columns": int(selected.shape[1]),
        "column_names": list(selected.columns),
        "table_fingerprint_sha256": table_fingerprint,
        "artifact_fingerprint_sha256": metadata["artifact_fingerprint_sha256"],
        "operation": "pure_artifact_view",
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return selected, summary


def render_augur_plot(
    tables: Mapping[str, DataFrame],
    parent_summary: Mapping[str, Any],
    metadata: Mapping[str, Any],
    *,
    view: str,
    top_n: int,
    max_plot_rows: int,
    max_image_pixels: int,
) -> tuple[bytes, dict[str, Any]]:
    import io

    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    if view not in AUGUR_PLOT_VIEWS:
        raise ValueError(f"Unsupported Augur plot view: {view!r}.")
    for value, label in (
        (top_n, "top_n"),
        (max_plot_rows, "max_plot_rows"),
        (max_image_pixels, "max_image_pixels"),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"Augur plot {label} must be an integer >= 1.")
    if top_n > 200:
        raise ValueError("Augur plot top_n cannot exceed 200.")
    owned_tables, owned_parent, owned_metadata = validate_augur_portable(tables, parent_summary, metadata)
    source_view = "predictions" if view == "roc" else view
    source = owned_tables[source_view]
    if len(source) > max_plot_rows:
        raise ValueError(
            f"Augur {source_view} plot evidence has {len(source):,} rows, exceeding max_plot_rows={max_plot_rows:,}."
        )
    width, height, dpi = 10.0, 6.0, 120
    nominal_pixels = int(width * dpi) * int(height * dpi)
    if nominal_pixels > max_image_pixels:
        raise MemoryError(
            f"Augur plot requires {nominal_pixels:,} nominal pixels, exceeding max_image_pixels={max_image_pixels:,}."
        )
    figure = Figure(figsize=(width, height))
    axis = figure.subplots()

    def plot_label(value: str) -> str:
        return value if len(value) <= 60 else f"{value[:57]}..."
    displayed = min(top_n, len(owned_tables["priorities"]))
    available = len(owned_tables["priorities"])
    displayed_kind = "ranked populations"
    roc_curve_auc: dict[str, float] = {}
    if view == "priorities":
        selected = owned_tables["priorities"].head(top_n).iloc[::-1]
        axis.barh([plot_label(value) for value in selected["population"]], selected["mean_augur_score"], color="#4C78A8")
        axis.axvline(0.5, color="#666666", linestyle="--", linewidth=1)
        axis.set(xlabel="Mean cross-validation AUC", ylabel="Population", xlim=(0.0, 1.02))
        axis.set_xticks(np.linspace(0.0, 1.0, 6))
        axis.set_title("Augur population priorities")
    elif view == "cross_validation":
        populations = owned_tables["priorities"].head(top_n)["population"].tolist()[::-1]
        distributions = [
            owned_tables["cross_validation"].loc[
                owned_tables["cross_validation"]["population"] == population, "auc"
            ].to_numpy(dtype=float)
            for population in populations
        ]
        positions = list(range(1, len(populations) + 1))
        boxes = axis.boxplot(distributions, positions=positions, orientation="horizontal", patch_artist=True)
        for box in boxes["boxes"]:
            box.set_facecolor("#72B7B2")
        for position, values in zip(positions, distributions, strict=True):
            offsets = np.linspace(-0.12, 0.12, len(values)) if len(values) > 1 else np.zeros(1)
            axis.scatter(values, position + offsets, color="#1F5A85", s=16, zorder=3)
        axis.axvline(0.5, color="#666666", linestyle="--", linewidth=1)
        axis.set_yticks(positions, labels=[plot_label(value) for value in populations])
        axis.set(xlabel="Fold AUC", ylabel="Population", xlim=(-0.02, 1.02))
        axis.set_xticks(np.linspace(0.0, 1.0, 6))
        axis.set_title("Augur cross-validation distribution")
    elif view == "feature_importance":
        aggregated = (
            owned_tables["feature_importance"]
            .groupby(["population", "gene"], observed=True, as_index=False)["importance"]
            .mean()
        )
        aggregated["magnitude"] = aggregated["importance"].abs()
        selected = aggregated.sort_values(
            ["magnitude", "population", "gene"], ascending=[False, True, True], kind="mergesort"
        ).head(top_n).iloc[::-1]
        labels = [
            plot_label(f"{population} · {gene}")
            for population, gene in selected[["population", "gene"]].itertuples(index=False, name=None)
        ]
        colors = ["#4C78A8" if value >= 0 else "#E45756" for value in selected["importance"]]
        axis.barh(labels, selected["importance"], color=colors)
        axis.axvline(0.0, color="#666666", linewidth=1)
        axis.set(xlabel="Mean fold importance", ylabel="Population · gene")
        axis.set_title("Augur feature importance")
        displayed = len(selected)
        available = len(aggregated)
        displayed_kind = "population-gene pairs"
    elif view == "roc":
        populations = owned_tables["priorities"].head(top_n)["population"].tolist()
        mean_auc = owned_tables["priorities"].set_index("population")["mean_auc"]
        for population in populations:
            population_predictions = owned_tables["predictions"].loc[
                owned_tables["predictions"]["population"] == population
            ]
            fold_curves = []
            for _, fold_predictions in population_predictions.groupby(
                ["subsample", "fold"], observed=True, sort=False
            ):
                labels = fold_predictions["true_label"].to_numpy(dtype=int)
                scores = fold_predictions["prediction_score"].to_numpy(dtype=float)
                order = np.argsort(-scores, kind="stable")
                ordered_labels = labels[order]
                ordered_scores = scores[order]
                threshold_ends = np.flatnonzero(np.r_[np.diff(ordered_scores) != 0, True])
                true_positives = np.cumsum(ordered_labels)[threshold_ends]
                false_positives = 1 + threshold_ends - true_positives
                fold_tpr = np.r_[0.0, true_positives / true_positives[-1]]
                fold_fpr = np.r_[0.0, false_positives / false_positives[-1]]
                unique_fpr = np.unique(fold_fpr)
                fold_curves.append(
                    (
                        unique_fpr,
                        np.asarray([fold_tpr[fold_fpr == value].max() for value in unique_fpr]),
                    )
                )
            fpr = np.unique(np.concatenate([curve[0] for curve in fold_curves]))
            tpr = np.mean(
                np.vstack(
                    [
                        fold_tpr[np.searchsorted(fold_fpr, fpr, side="right") - 1]
                        for fold_fpr, fold_tpr in fold_curves
                    ]
                ),
                axis=0,
            )
            curve_auc = float(np.sum(tpr[:-1] * np.diff(fpr)))
            if not math.isclose(curve_auc, float(mean_auc[population]), rel_tol=1e-10, abs_tol=1e-12):
                raise RuntimeError(f"Augur mean ROC curve AUC disagrees with fold evidence for {population!r}.")
            roc_curve_auc[population] = curve_auc
            axis.plot(
                np.r_[0.0, fpr],
                np.r_[0.0, tpr],
                linewidth=2,
                drawstyle="steps-post",
                label=plot_label(f"{population} (mean AUC {mean_auc[population]:.3f})"),
            )
        axis.plot([0.0, 1.0], [0.0, 1.0], color="#666666", linestyle="--", linewidth=1)
        positive_condition = str(owned_parent["key_results"]["prediction_label_encoding"]["positive_condition"])
        axis.set(
            xlabel=f"False-positive rate (positive: {plot_label(positive_condition)})",
            ylabel="True-positive rate",
            xlim=(-0.02, 1.02),
            ylim=(-0.02, 1.02),
        )
        axis.set_xticks(np.linspace(0.0, 1.0, 6))
        axis.set_yticks(np.linspace(0.0, 1.0, 6))
        axis.set_title("Augur mean fold ROC")
        axis.legend(loc="lower right", fontsize="small")
    else:
        raise NotImplementedError(f"Augur plot view {view!r} is not implemented.")
    figure.tight_layout()
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    tight_bbox = figure.get_tightbbox(canvas.get_renderer())
    tight_pixels = math.ceil((tight_bbox.width + 0.2) * dpi) * math.ceil((tight_bbox.height + 0.2) * dpi)
    if tight_pixels > max_image_pixels:
        raise MemoryError(
            f"Augur plot requires {tight_pixels:,} tight-layout pixels, exceeding "
            f"max_image_pixels={max_image_pixels:,}."
        )
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    png = buffer.getvalue()
    truncated = available > displayed
    warnings = [str(value) for value in owned_parent["warnings"]]
    if truncated:
        warnings.append(f"The plot displays {displayed} of {available} {displayed_kind}.")
    summary = {
        "schema_version": "openbio-singlecell/augur-plot-summary/v1",
        "node_id": AUGUR_PLOT_NODE_ID,
        "status": "diagnostic_visualization",
        "methods": f"Rendered the {view!r} view from validated immutable Augur diagnostic evidence.",
        "results": f"Displayed {displayed:,} item(s) in the Augur {view.replace('_', ' ')} view.",
        "key_results": {
            "view": view,
            "displayed_items": displayed,
            "available_items": available,
            "displayed_kind": displayed_kind,
            "source_rows": int(len(source)),
            "image_pixels": tight_pixels,
            **({"roc_curve_auc": roc_curve_auc} if view == "roc" else {}),
            "artifact_fingerprint_sha256": owned_metadata["artifact_fingerprint_sha256"],
        },
        "parameters": {
            "view": view,
            "top_n": top_n,
            "max_plot_rows": max_plot_rows,
            "max_image_pixels": max_image_pixels,
        },
        "references": copy.deepcopy(owned_parent["references"]),
        "software_versions": copy.deepcopy(owned_parent["software_versions"]),
        "warnings": warnings,
        "limitations": [
            "This visualization is diagnostic evidence and does not add Condition-level inference.",
            "Only the retained Augur result is plotted; the classifier analysis is not rerun.",
            *(
                [
                    "The displayed ROC is the deterministic mean of fold curves; it is diagnostic, not a new "
                    "population-level estimate."
                ]
                if view == "roc"
                else []
            ),
        ],
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return png, summary


def augur_plot_code(*, view: str, top_n: int, max_plot_rows: int, max_image_pixels: int) -> str:
    helpers = (
        _canonical_json_sha256,
        _append_string,
        _frame_fingerprint,
        _validate_canonical_tables,
        _is_compatible_pertpy_1_3,
        validate_augur_portable,
        render_augur_plot,
    )
    implementation = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    return f'''from __future__ import annotations

import copy
import hashlib
import json
import math
import re
from collections.abc import Mapping

AUGUR_ARTIFACT_TYPE = {AUGUR_ARTIFACT_TYPE!r}
AUGUR_ARTIFACT_SCHEMA_VERSION = {AUGUR_ARTIFACT_SCHEMA_VERSION!r}
AUGUR_PRODUCER_NODE_ID = {AUGUR_PRODUCER_NODE_ID!r}
AUGUR_PLOT_NODE_ID = {AUGUR_PLOT_NODE_ID!r}
AUGUR_PRODUCER_SCHEMA = {AUGUR_PRODUCER_SCHEMA!r}
AUGUR_SUMMARY_SCHEMA = {AUGUR_SUMMARY_SCHEMA!r}
AUGUR_VIEWS = {AUGUR_VIEWS!r}
AUGUR_PLOT_VIEWS = {AUGUR_PLOT_VIEWS!r}
PRIORITY_COLUMNS = {PRIORITY_COLUMNS!r}
CROSS_VALIDATION_COLUMNS = {CROSS_VALIDATION_COLUMNS!r}
FEATURE_IMPORTANCE_COLUMNS = {FEATURE_IMPORTANCE_COLUMNS!r}
PREDICTION_COLUMNS = {PREDICTION_COLUMNS!r}
_PRIORITY_FLOAT_COLUMNS = {_PRIORITY_FLOAT_COLUMNS!r}
_PRIORITY_INTEGER_COLUMNS = {_PRIORITY_INTEGER_COLUMNS!r}

{implementation}


def plot_augur_result(result):
    """Return (PNG bytes, strict plot summary) without rerunning Augur."""
    if getattr(result, "artifact_type", None) != AUGUR_ARTIFACT_TYPE:
        raise TypeError("Augur Plot requires an OPENBIO_AUGUR_RESULT artifact.")
    portable_tables = getattr(result, "portable_tables", None)
    if not callable(portable_tables):
        raise TypeError("Augur Plot requires portable typed artifact evidence.")
    metadata = result.metadata
    if result.fingerprint != metadata.get("artifact_fingerprint_sha256"):
        raise ValueError("Augur artifact fingerprint identity is invalid.")
    return render_augur_plot(
        portable_tables(),
        result.summary,
        metadata,
        view={view!r},
        top_n={top_n!r},
        max_plot_rows={max_plot_rows!r},
        max_image_pixels={max_image_pixels!r},
    )
'''


def augur_code(**parameters: Any) -> str:
    portable_parameters = dict(parameters)
    portable_parameters["source_kind"] = portable_parameters.pop("source")
    if portable_parameters.get("technical_batch_key") is None:
        portable_parameters["technical_batch_key"] = ""
    portable_parameters.pop("source_label", None)
    portable_parameters.pop("min_cells", None)
    portable_parameters.pop("resolved_min_cells", None)
    portable_parameters.pop("feature_perc", None)
    portable_parameters.pop("var_quantile", None)
    portable_parameters.pop("span", None)
    portable_parameters.pop("filter_negative_residuals", None)
    portable_parameters.pop("augur_mode", None)
    portable_parameters.pop("select_variance_features", None)
    portable_parameters.pop("key_added", None)
    portable_parameters.pop("zero_division", None)
    portable_parameters.pop("estimator_hyperparameters", None)
    helpers = (
        _canonical_json_sha256,
        _strict_text,
        _canonical_labels,
        _count_matrix_fingerprint,
        _observation_fingerprint,
        _append_string,
        _frame_fingerprint,
        _validate_canonical_tables,
        AugurResult,
        _build_augur_artifact,
        _build_augur_result,
        validate_augur_portable,
        validate_augur_result,
        _require_pertpy,
        _pertpy_version,
        _is_compatible_pertpy_1_3,
        _validate_pertpy_1_3_interface,
        _with_prediction_evidence,
        _prediction_retaining_augur_class,
        _validate_integer,
        _backend_float_column,
        _backend_integer_column,
        _canonicalize_backend_results,
        _run_augur_artifact_owned,
        run_augur_artifact,
        run_augur_analysis,
        run_augur_portable,
    )
    implementation = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    return f'''from __future__ import annotations

import copy
import hashlib
import inspect
import json
import math
import platform
import re
import warnings as python_warnings
from collections.abc import Mapping, Sequence
from importlib import metadata as importlib_metadata
from typing import Any

AUGUR_ARTIFACT_TYPE = {AUGUR_ARTIFACT_TYPE!r}
AUGUR_ARTIFACT_SCHEMA_VERSION = {AUGUR_ARTIFACT_SCHEMA_VERSION!r}
AUGUR_PRODUCER_NODE_ID = {AUGUR_PRODUCER_NODE_ID!r}
AUGUR_PRODUCER_SCHEMA = {AUGUR_PRODUCER_SCHEMA!r}
AUGUR_SUMMARY_SCHEMA = {AUGUR_SUMMARY_SCHEMA!r}
AUGUR_VIEWS = {AUGUR_VIEWS!r}
AUGUR_CLASSIFIERS = {AUGUR_CLASSIFIERS!r}
PRIORITY_COLUMNS = {PRIORITY_COLUMNS!r}
CROSS_VALIDATION_COLUMNS = {CROSS_VALIDATION_COLUMNS!r}
FEATURE_IMPORTANCE_COLUMNS = {FEATURE_IMPORTANCE_COLUMNS!r}
PREDICTION_COLUMNS = {PREDICTION_COLUMNS!r}
_PRIORITY_FLOAT_COLUMNS = {_PRIORITY_FLOAT_COLUMNS!r}
_PRIORITY_INTEGER_COLUMNS = {_PRIORITY_INTEGER_COLUMNS!r}
AUGUR_REFERENCES = {AUGUR_REFERENCES!r}


def collect_software_versions(packages):
    """Collect the same strict version payload without importing OpenBio."""
    versions = {{
        "python": platform.python_version(),
        "openbio-singlecell": {PLUGIN_VERSION!r},
    }}
    for package in dict.fromkeys(packages):
        if package in versions:
            continue
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


{implementation}


def run_augur_prioritization(adata, pertpy_module=None):
    """Return (canonical_table_mapping, strict_summary) without mutating adata."""
    return run_augur_portable(
        adata,
        **{portable_parameters!r},
        pertpy_module=pertpy_module,
    )
'''


def augur_results_code(*, view: str) -> str:
    if view not in AUGUR_VIEWS:
        raise ValueError(f"Unsupported Augur result view: {view!r}.")
    return dedent(
        f'''\
        def select_augur_results(result):
            """Return (selected_dataframe, inherited_strict_summary) without rerunning Augur."""
            from openbio_singlecell.augur import select_augur_view

            return select_augur_view(result, view={view!r})
        '''
    )


__all__ = [
    "AUGUR_ARTIFACT_SCHEMA_VERSION",
    "AUGUR_ARTIFACT_TYPE",
    "AUGUR_CLASSIFIERS",
    "AUGUR_PRODUCER_NODE_ID",
    "AUGUR_PRODUCER_SCHEMA",
    "AUGUR_PLOT_VIEWS",
    "AUGUR_PLOT_NODE_ID",
    "AUGUR_REFERENCES",
    "AUGUR_RESULTS_NODE_ID",
    "AUGUR_SUMMARY_SCHEMA",
    "AUGUR_VIEWS",
    "AugurResult",
    "CROSS_VALIDATION_COLUMNS",
    "FEATURE_IMPORTANCE_COLUMNS",
    "PREDICTION_COLUMNS",
    "PRIORITY_COLUMNS",
    "augur_code",
    "augur_plot_code",
    "augur_results_code",
    "run_augur_analysis",
    "run_augur_artifact",
    "run_augur_portable",
    "render_augur_plot",
    "select_augur_view",
    "validate_augur_portable",
    "validate_augur_result",
]
