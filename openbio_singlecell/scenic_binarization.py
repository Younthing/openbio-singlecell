from __future__ import annotations

import collections
import copy
import hashlib
import inspect
import json
import math
from collections.abc import Mapping
from textwrap import dedent
from typing import TYPE_CHECKING, Any

from . import PLUGIN_VERSION
from .pyscenic_import import _software_versions
from .scenic_artifact import (
    SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION,
    SCENIC_BINARY_ARTIFACT_TYPE,
    SCENIC_BINARY_PRODUCER_NODE_ID,
    SCENIC_BINARY_PRODUCER_SCHEMA,
    SCENIC_MEMBERSHIP_COLUMNS,
    SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION,
    SCENIC_RESULT_ARTIFACT_TYPE,
    SCENIC_RESULT_PRODUCER_NODE_ID,
    SCENIC_RESULT_PRODUCER_SCHEMA,
    SCENICBinaryArtifact,
    SCENICResultArtifact,
    _activity_bytes,
    _binary_fingerprint,
    _canonical_axis,
    _canonical_json,
    _canonical_membership,
    _result_fingerprint,
    _strict_provenance,
    validate_scenic_result_artifact,
)

if TYPE_CHECKING:
    from pandas import DataFrame


SCENIC_BINARIZATION_SUMMARY_SCHEMA = "openbio-singlecell/scenic-binarization-summary/v1"
SCENIC_THRESHOLD_COLUMNS = (
    "regulon",
    "threshold",
    "threshold_source",
    "active_cells",
    "total_cells",
    "active_proportion",
)
HDT_DIP_SIMULATIONS = 1000
SCENIC_BINARIZATION_REFERENCES = [
    {
        "citation": (
            "Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network "
            "analysis. Nature Protocols. 2020;15:2247-2276."
        ),
        "doi": "10.1038/s41596-020-0336-2",
        "url": "https://doi.org/10.1038/s41596-020-0336-2",
        "kind": "method",
    },
    {
        "citation": (
            "Hartigan JA, Hartigan PM. The dip test of unimodality. The Annals of Statistics. 1985;13:70-84."
        ),
        "doi": "10.1214/aos/1176346577",
        "url": "https://doi.org/10.1214/aos/1176346577",
        "kind": "method",
    },
    {
        "citation": (
            "pySCENIC 0.12.1 (commit ce41b61b6570490949bd12b514e9f6de46d19c1f) "
            "binarization.py derive_threshold HDT branch and bundled diptest.py implementations."
        ),
        "doi": None,
        "url": (
            "https://github.com/aertslab/pySCENIC/tree/"
            "ce41b61b6570490949bd12b514e9f6de46d19c1f/src/pyscenic"
        ),
        "kind": "software_documentation",
    },
]


# GPL-3.0-or-later compatibility rewrite derived from pySCENIC 0.12.1 commit
# ce41b61b6570490949bd12b514e9f6de46d19c1f, src/pyscenic/diptest.py
# (_gcm_, _lcm_, _touch_diffs_, dip, diptst). That upstream file credits Johannes Bauer's
# tatome/dip_test commit a0e3d448a4b266f54ec63a5b3d5be351fbd1db1c and BenjaminDoran/unidip.
# See THIRD_PARTY_NOTICES.md. Generated standalone code carries the same attribution.
def _gcm(cdf: Any, indices: Any, *, numpy: Any) -> tuple[Any, Any]:
    work_cdf = cdf
    work_indices = indices
    gcm = [work_cdf[0]]
    touchpoints = [0]
    while len(work_cdf) > 1:
        distances = work_indices[1:] - work_indices[0]
        slopes = (work_cdf[1:] - work_cdf[0]) / distances
        minimum_slope = slopes.min()
        minimum_position = numpy.where(slopes == minimum_slope)[0][0] + 1
        gcm.extend(work_cdf[0] + distances[:minimum_position] * minimum_slope)
        touchpoints.append(touchpoints[-1] + minimum_position)
        work_cdf = work_cdf[minimum_position:]
        work_indices = work_indices[minimum_position:]
    return numpy.asarray(gcm), numpy.asarray(touchpoints)


def _lcm(cdf: Any, indices: Any, *, numpy: Any) -> tuple[Any, Any]:
    values, touchpoints = _gcm(1 - cdf[::-1], indices.max() - indices[::-1], numpy=numpy)
    return 1 - values[::-1], len(cdf) - 1 - touchpoints[::-1]


def _touch_differences(part1: Any, part2: Any, touchpoints: Any, *, numpy: Any) -> tuple[float, Any]:
    differences = numpy.abs(part2[touchpoints] - part1[touchpoints])
    return float(differences.max()), differences


def _dip_function(data: Any, *, numpy: Any, is_histogram: bool = False, just_dip: bool = False) -> Any:
    if is_histogram:
        histogram = data
        indices = numpy.arange(len(histogram))
    else:
        counts = collections.Counter(data)
        indices = numpy.sort(list(counts.keys()))
        histogram = numpy.asarray([counts[index] for index in indices])
    if len(indices) <= 4 or indices[0] == indices[-1]:
        left: list[Any] = []
        right = [1]
        dip = 0.0
        return dip if just_dip else (dip, (None, indices, left, None, right, None))
    cdf = numpy.cumsum(histogram, dtype=float)
    cdf /= cdf[-1]
    work_indices = indices
    work_histogram = numpy.asarray(histogram, dtype=float) / numpy.sum(histogram)
    work_cdf = cdf
    maximum_difference = 0.0
    left = [0]
    right = [1]
    while True:
        left_part, left_touchpoints = _gcm(
            work_cdf - work_histogram, work_indices, numpy=numpy
        )
        right_part, right_touchpoints = _lcm(work_cdf, work_indices, numpy=numpy)
        left_distance, left_differences = _touch_differences(
            left_part, right_part, left_touchpoints, numpy=numpy
        )
        right_distance, right_differences = _touch_differences(
            left_part, right_part, right_touchpoints, numpy=numpy
        )
        if right_distance > left_distance:
            right_index = right_touchpoints[right_distance == right_differences][-1]
            left_index = left_touchpoints[left_touchpoints <= right_index][-1]
            distance = right_distance
        else:
            left_index = left_touchpoints[left_distance == left_differences][0]
            right_index = right_touchpoints[right_touchpoints >= left_index][0]
            distance = left_distance
        left_difference = numpy.abs(left_part[: left_index + 1] - work_cdf[: left_index + 1]).max()
        right_difference = numpy.abs(
            right_part[right_index:] - work_cdf[right_index:] + work_histogram[right_index:]
        ).max()
        if distance <= maximum_difference or right_index == 0 or left_index == len(work_cdf):
            dip = max(
                numpy.abs(cdf[: len(left)] - left).max(),
                numpy.abs(cdf[-len(right) - 1 : -1] - right).max(),
            )
            return dip / 2 if just_dip else (
                dip / 2,
                (cdf, indices, left, left_part, right, right_part),
            )
        maximum_difference = max(maximum_difference, left_difference, right_difference)
        work_cdf = work_cdf[left_index : right_index + 1]
        work_indices = work_indices[left_index : right_index + 1]
        work_histogram = work_histogram[left_index : right_index + 1]
        left[len(left) :] = left_part[1 : left_index + 1]
        right[:0] = right_part[right_index:-1]


def _dip_test(data: Any, *, numpy: Any, simulations: int = 1000) -> tuple[float, float | None]:
    dip, (_, indices, _left, _left_part, _right, _right_part) = _dip_function(
        data, numpy=numpy
    )
    uniforms = numpy.random.uniform(size=simulations * indices.shape[0]).reshape(
        simulations, indices.shape[0]
    )
    uniform_dips = numpy.apply_along_axis(
        lambda row: _dip_function(row, numpy=numpy, just_dip=True), 1, uniforms
    )
    p_value = (
        None
        if uniform_dips.sum() == 0
        else float((numpy.less(dip, uniform_dips).sum() + 1) / (float(simulations) + 1))
    )
    return float(dip), p_value


def _derive_threshold_exact(
    auc_matrix: DataFrame,
    regulon_name: str,
    seed: int,
    *,
    numpy: Any,
    scipy_stats: Any,
    minimize_scalar: Any,
    mixture: Any,
) -> float:
    # GPL-3.0-or-later compatibility transcription of the HDT branch in
    # pySCENIC 0.12.1 src/pyscenic/binarization.py::derive_threshold at commit
    # ce41b61b6570490949bd12b514e9f6de46d19c1f. Removed NumPy aliases are replaced by
    # np.sort/builtin float; seed, dip-test, GMM, KDE, and minimization semantics remain.
    data = auc_matrix[regulon_name].to_numpy(dtype=float, copy=True)
    if seed:
        numpy.random.seed(seed=seed)
    _dip, p_value = _dip_test(numpy.sort(data), numpy=numpy)
    is_bimodal = p_value is not None and p_value <= 0.05
    if not is_bimodal:
        return float(data.mean() + 2.0 * data.std())
    gmm = mixture.GaussianMixture(
        n_components=2, covariance_type="full", random_state=seed
    ).fit(data.reshape(-1, 1))
    result = minimize_scalar(
        fun=scipy_stats.gaussian_kde(data), bounds=sorted(gmm.means_), method="bounded"
    ).x
    return float(result[0] if hasattr(result, "__len__") else result)


def _threshold_overrides_frame(value: Any, *, regulons: list[str], numpy: Any, pandas: Any) -> dict[str, float]:
    if value is None:
        return {}
    if not isinstance(value, pandas.DataFrame):
        raise TypeError("SCENIC threshold_overrides must be a pandas DataFrame when provided.")
    if value.columns.tolist() != ["regulon", "threshold"]:
        raise ValueError("SCENIC threshold_overrides columns must be exactly ['regulon', 'threshold'].")
    result: dict[str, float] = {}
    known = set(regulons)
    for row in value.itertuples(index=False, name=None):
        regulon, threshold = row
        if not isinstance(regulon, str) or not regulon or regulon != regulon.strip():
            raise ValueError("SCENIC override regulon IDs must be canonical nonblank strings.")
        if regulon not in known:
            raise ValueError(f"SCENIC threshold override references unknown regulon {regulon!r}.")
        if regulon in result:
            raise ValueError(f"SCENIC threshold overrides contain duplicate regulon {regulon!r}.")
        if isinstance(threshold, bool) or not isinstance(
            threshold, (int, float, numpy.integer, numpy.floating)
        ):
            raise TypeError("SCENIC override thresholds must be numeric.")
        threshold = float(threshold)
        if not math.isfinite(threshold):
            raise ValueError("SCENIC override thresholds must be finite.")
        result[regulon] = threshold
    return result


def binarize_scenic_activity(
    scenic_result: Any,
    *,
    random_seed: int = 1,
    threshold_overrides: DataFrame | None = None,
    max_dense_bytes: int = 1_073_741_824,
    openbio_version: str = PLUGIN_VERSION,
    _portable_artifact: bool = False,
    _threshold_backend: Any | None = None,
) -> tuple[SCENICBinaryArtifact, DataFrame, dict[str, Any]]:
    import numpy as np
    import pandas as pd
    import scipy
    import sklearn

    if isinstance(random_seed, bool) or not isinstance(random_seed, int) or not 1 <= random_seed <= 2**31 - 1:
        raise ValueError("SCENIC binarization random_seed must be an integer in [1, 2^31-1].")
    if isinstance(max_dense_bytes, bool) or not isinstance(max_dense_bytes, int) or max_dense_bytes < 1:
        raise TypeError("SCENIC binarization max_dense_bytes must be a positive integer.")
    if isinstance(scenic_result, SCENICResultArtifact):
        declared_cells = len(scenic_result.observation_names)
        declared_regulon_names = list(scenic_result.regulon_names)
    elif _portable_artifact and isinstance(scenic_result, dict):
        observations = scenic_result.get("observations")
        regulons = scenic_result.get("regulons")
        declared_cells = len(observations) if isinstance(observations, (list, tuple)) else 0
        declared_regulon_names = list(regulons) if isinstance(regulons, (list, tuple)) else []
    else:
        declared_cells = 0
        declared_regulon_names = []
    declared_regulons = len(declared_regulon_names)
    overrides = _threshold_overrides_frame(
        threshold_overrides,
        regulons=declared_regulon_names,
        numpy=np,
        pandas=pd,
    )
    matrix_working_bytes = declared_cells * declared_regulons * 96
    dip_test_working_bytes = (
        declared_cells * HDT_DIP_SIMULATIONS * 16
        if declared_regulons > len(overrides)
        else 0
    )
    estimated_peak_dense_bytes = matrix_working_bytes + dip_test_working_bytes
    if estimated_peak_dense_bytes > max_dense_bytes:
        raise ValueError(
            f"SCENIC binarization requires an estimated {estimated_peak_dense_bytes:,} dense working bytes, exceeding "
            f"max_dense_bytes={max_dense_bytes:,}."
        )
    activity, _membership, source_provenance, source_metadata = validate_scenic_result_artifact(
        scenic_result,
        exact_type=not _portable_artifact,
        numpy=np,
        pandas=pd,
    )
    cells, regulon_count = activity.shape
    if (cells, regulon_count) != (declared_cells, declared_regulons):
        raise ValueError("SCENIC binarization declared artifact shape differs from its validated activity matrix.")
    regulons = activity.columns.tolist()
    if regulons != declared_regulon_names:
        raise ValueError("SCENIC binarization declared regulon axis differs from the validated activity matrix.")
    backend = _derive_threshold_exact if _threshold_backend is None else _threshold_backend
    if not callable(backend):
        raise TypeError("SCENIC threshold backend must be callable.")
    private_activity = activity.copy(deep=True)
    activity_before = private_activity.copy(deep=True)
    global_state = np.random.get_state()
    thresholds: dict[str, float] = {}
    low_information: list[str] = []
    try:
        for regulon in regulons:
            values = private_activity[regulon].to_numpy(dtype=float, copy=False)
            distinct = int(np.unique(values).size)
            if distinct <= 4:
                low_information.append(regulon)
            if regulon in overrides:
                threshold = overrides[regulon]
            else:
                if _threshold_backend is None:
                    threshold = backend(
                        private_activity,
                        regulon,
                        random_seed,
                        numpy=np,
                        scipy_stats=scipy.stats,
                        minimize_scalar=scipy.optimize.minimize_scalar,
                        mixture=sklearn.mixture,
                    )
                else:
                    threshold = backend(private_activity, regulon, random_seed)
            if isinstance(threshold, bool) or not isinstance(
                threshold, (int, float, np.integer, np.floating)
            ):
                raise RuntimeError(f"SCENIC threshold backend returned a nonnumeric value for {regulon!r}.")
            threshold = float(threshold)
            if not math.isfinite(threshold):
                raise RuntimeError(f"SCENIC threshold backend returned a non-finite value for {regulon!r}.")
            thresholds[regulon] = threshold
        try:
            pd.testing.assert_frame_equal(private_activity, activity_before, check_exact=True)
        except AssertionError as exc:
            raise RuntimeError("SCENIC threshold backend modified its private activity input.") from exc
    finally:
        np.random.set_state(global_state)
    if list(thresholds) != regulons:
        raise RuntimeError("SCENIC threshold family/order differs from the complete regulon family.")
    threshold_series = pd.Series(thresholds, index=regulons, dtype=float)
    binary = activity.gt(threshold_series, axis="columns")
    if binary.shape != activity.shape or binary.columns.tolist() != regulons:
        raise RuntimeError("SCENIC binary output axes differ from the source artifact.")
    threshold_rows = []
    for regulon in regulons:
        active_cells = int(binary[regulon].sum())
        threshold_rows.append(
            {
                "regulon": regulon,
                "threshold": thresholds[regulon],
                "threshold_source": "manual_override" if regulon in overrides else "pyscenic_0.12.1_hdt",
                "active_cells": active_cells,
                "total_cells": cells,
                "active_proportion": float(active_cells / cells),
            }
        )
    threshold_table = pd.DataFrame(threshold_rows, columns=SCENIC_THRESHOLD_COLUMNS)
    provenance = {
        "artifact_type": SCENIC_BINARY_ARTIFACT_TYPE,
        "artifact_schema_version": SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION,
        "producer_node_id": SCENIC_BINARY_PRODUCER_NODE_ID,
        "producer_schema": SCENIC_BINARY_PRODUCER_SCHEMA,
        "source_scenic_artifact_fingerprint_sha256": source_metadata[
            "artifact_fingerprint_sha256"
        ],
        "source_manifest_sha256": source_provenance["manifest_sha256"],
        "random_seed": random_seed,
        "algorithm": "pyscenic-0.12.1-hdt-compatible-single-worker",
        "strict_comparison": "auc > threshold",
        "thresholds_sha256": hashlib.sha256(
            _canonical_json(threshold_rows).encode("utf-8")
        ).hexdigest(),
    }
    binary_artifact = SCENICBinaryArtifact(binary, provenance)
    warnings = [
        "Binarization is a data-dependent descriptive heuristic; thresholds are not p-values.",
        "A cell is active only when AUC is strictly greater than the threshold; equality is off.",
        "Thresholds depend on the complete included-cell composition and may change when cells are added or removed.",
    ]
    if low_information:
        warnings.append(
            f"{len(low_information)} regulon(s) had at most four distinct AUC values and therefore followed the "
            "official unimodal mean-plus-two-standard-deviations branch: {low_information!r}."
        )
    out_of_range_overrides = [
        regulon for regulon, threshold in overrides.items() if threshold < 0.0 or threshold > 1.0
    ]
    if out_of_range_overrides:
        warnings.append(
            "Manual thresholds outside the AUCell [0, 1] range were retained as explicit expert choices; they "
            "necessarily produce all-on or all-off calls for the named regulons: "
            f"{out_of_range_overrides!r}."
        )
    parameters = {
        "random_seed": random_seed,
        "num_workers": 1,
        "method": "hdt",
        "dip_test_simulations": HDT_DIP_SIMULATIONS,
        "strict_comparison": ">",
        "manual_override_regulons": list(overrides),
        "out_of_range_override_regulons": out_of_range_overrides,
        "max_dense_bytes": max_dense_bytes,
    }
    summary = {
        "schema_version": SCENIC_BINARIZATION_SUMMARY_SCHEMA,
        "node_id": SCENIC_BINARY_PRODUCER_NODE_ID,
        "status": "descriptive_binary_regulon_activity",
        "methods": (
            "Applied a NumPy-2-compatible exact transcription of pySCENIC 0.12.1 HDT threshold derivation in one "
            "process with a nonzero seed and isolated legacy NumPy RNG state, followed by the official strict AUC > "
            "threshold rule."
        ),
        "results": (
            f"Binarized {cells:,} cells across {regulon_count:,} regulons; "
            f"{int(binary.to_numpy(dtype=bool).sum()):,} cell-regulon pairs were active."
        ),
        "key_results": {
            "scientific_label": "Exploratory binary regulon activity",
            "cells": cells,
            "regulons": regulon_count,
            "active_calls": int(binary.to_numpy(dtype=bool).sum()),
            "inactive_calls": int(binary.size - binary.to_numpy(dtype=bool).sum()),
            "threshold_min": float(threshold_table["threshold"].min()),
            "threshold_max": float(threshold_table["threshold"].max()),
            "manual_override_count": len(overrides),
            "out_of_range_override_regulons": out_of_range_overrides,
            "estimated_peak_dense_bytes": estimated_peak_dense_bytes,
            "low_information_regulons": low_information,
            "leading_activity_calls": threshold_table.sort_values(
                ["active_proportion", "regulon"], ascending=[False, True]
            ).head(20).to_dict(orient="records"),
            "source_scenic_artifact_fingerprint_sha256": source_metadata[
                "artifact_fingerprint_sha256"
            ],
            "binary_artifact_fingerprint_sha256": binary_artifact.artifact_fingerprint_sha256,
            "external_scenic_provenance": source_provenance,
        },
        "parameters": parameters,
        "references": copy.deepcopy(SCENIC_BINARIZATION_REFERENCES),
        "software_versions": {
            **_software_versions(
                ["numpy", "pandas", "scipy", "scikit-learn"], openbio_version=openbio_version
            ),
            "pyscenic-method": "0.12.1",
        },
        "warnings": warnings,
        "limitations": [
            "Binary calls discard continuous AUCell information and are not statistical significance decisions.",
            "Binarization does not establish TF binding, causality, cell identity, or a Condition effect.",
            "The output inherits expression/resource and regulon-inference limitations from the source SCENIC artifact.",
        ],
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return binary_artifact, threshold_table, summary


def scenic_binarization_code(
    *,
    parameters: Mapping[str, Any],
    function_name: str = "binarize_scenic_regulon_activity",
) -> str:
    if not isinstance(function_name, str) or not function_name.isidentifier():
        raise ValueError("Generated SCENIC binarization function_name must be a Python identifier.")
    portable_parameters = dict(parameters)
    portable_parameters["_portable_artifact"] = True
    helpers = (
        _canonical_json,
        _canonical_axis,
        _activity_bytes,
        _canonical_membership,
        _strict_provenance,
        _result_fingerprint,
        validate_scenic_result_artifact,
        _binary_fingerprint,
        SCENICBinaryArtifact,
        _software_versions,
        _gcm,
        _lcm,
        _touch_differences,
        _dip_function,
        _dip_test,
        _derive_threshold_exact,
        _threshold_overrides_frame,
        binarize_scenic_activity,
    )
    helper_source = "\n\n".join(dedent(inspect.getsource(helper)).strip() for helper in helpers)
    arguments = "\n".join(f"        {name}={value!r}," for name, value in portable_parameters.items())
    return f'''from __future__ import annotations

# Portions of this standalone function are GPL-3.0-or-later compatibility rewrites of
# pySCENIC 0.12.1 commit ce41b61b6570490949bd12b514e9f6de46d19c1f:
# src/pyscenic/diptest.py and the HDT branch of
# src/pyscenic/binarization.py::derive_threshold. The upstream dip-test file credits
# Johannes Bauer's tatome/dip_test commit a0e3d448a4b266f54ec63a5b3d5be351fbd1db1c
# and BenjaminDoran/unidip. See pySCENIC's GPL-3.0-or-later license and the source
# distribution's THIRD_PARTY_NOTICES.md.

import collections
import copy
import hashlib
import importlib.metadata
import json
import math
import platform
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

PLUGIN_VERSION = {PLUGIN_VERSION!r}
SCENIC_RESULT_ARTIFACT_TYPE = {SCENIC_RESULT_ARTIFACT_TYPE!r}
SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION = {SCENIC_RESULT_ARTIFACT_SCHEMA_VERSION!r}
SCENIC_RESULT_PRODUCER_NODE_ID = {SCENIC_RESULT_PRODUCER_NODE_ID!r}
SCENIC_RESULT_PRODUCER_SCHEMA = {SCENIC_RESULT_PRODUCER_SCHEMA!r}
SCENIC_BINARY_ARTIFACT_TYPE = {SCENIC_BINARY_ARTIFACT_TYPE!r}
SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION = {SCENIC_BINARY_ARTIFACT_SCHEMA_VERSION!r}
SCENIC_BINARY_PRODUCER_NODE_ID = {SCENIC_BINARY_PRODUCER_NODE_ID!r}
SCENIC_BINARY_PRODUCER_SCHEMA = {SCENIC_BINARY_PRODUCER_SCHEMA!r}
SCENIC_MEMBERSHIP_COLUMNS = {SCENIC_MEMBERSHIP_COLUMNS!r}
SCENIC_BINARIZATION_SUMMARY_SCHEMA = {SCENIC_BINARIZATION_SUMMARY_SCHEMA!r}
SCENIC_THRESHOLD_COLUMNS = {SCENIC_THRESHOLD_COLUMNS!r}
HDT_DIP_SIMULATIONS = {HDT_DIP_SIMULATIONS!r}
SCENIC_BINARIZATION_REFERENCES = {SCENIC_BINARIZATION_REFERENCES!r}
SCENICResultArtifact = ()

{helper_source}


def {function_name}(scenic_result, threshold_overrides=None):
    """Return immutable binary activity, complete thresholds, and a strict summary."""
    return binarize_scenic_activity(
        scenic_result,
        threshold_overrides=threshold_overrides,
{arguments}
    )
'''


__all__ = [
    "HDT_DIP_SIMULATIONS",
    "SCENIC_BINARIZATION_SUMMARY_SCHEMA",
    "SCENIC_THRESHOLD_COLUMNS",
    "binarize_scenic_activity",
    "scenic_binarization_code",
]
