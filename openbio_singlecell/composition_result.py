from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from typing import Any

COMPOSITION_MODEL_ARTIFACT_TYPE = "OPENBIO_COMPOSITION_MODEL_RESULT"
COMPOSITION_MODEL_ARTIFACT_SCHEMA_VERSION = 1

SCCODA_RESULT_COLUMNS = (
    "contrast",
    "condition_key",
    "reference_condition",
    "comparison_condition",
    "cell_type",
    "reference_cell_type",
    "model_coefficient",
    "hdi_lower",
    "hdi_upper",
    "posterior_sd",
    "inclusion_probability",
    "credible_effect",
    "estimated_fdr",
    "inclusion_probability_threshold",
    "realized_expected_fdr",
    "expected_count_reference",
    "expected_count_comparison",
    "compositional_log2_fold_change",
    "reference_constraint",
)

TASCCODA_RESULT_COLUMNS = (
    "effect_scope",
    "contrast",
    "condition_key",
    "reference_condition",
    "comparison_condition",
    "reference_cell_type",
    "effect_name",
    "hierarchy_level",
    "descendant_leaf_count",
    "descendant_leaves_json",
    "model_effect",
    "posterior_median",
    "hdi_lower",
    "hdi_upper",
    "posterior_sd",
    "selection_delta",
    "credible_effect",
    "selection_basis",
    "expected_count_reference",
    "expected_count_comparison",
    "compositional_log2_fold_change",
    "reference_constraint",
)

_POSTERIOR_VARIABLES = {
    "sccoda": {"intercept": "cell_type", "condition_effect": "cell_type"},
    "tasccoda": {
        "intercept": "cell_type",
        "hierarchy_node_effect": "hierarchy_node",
        "derived_leaf_effect": "cell_type",
        "theta": None,
    },
}
_SAMPLE_STATS = ("potential_energy", "num_steps", "step_size")
_MODEL_METADATA_FIELDS = {
    "condition_key",
    "reference_condition",
    "comparison_condition",
    "reference_cell_type",
    "annotation_key",
    "annotation_status",
    "model",
    "diagnostics",
    "selection",
    "software_versions",
}
_HIERARCHY_FIELDS = {
    "ancestor_keys",
    "leaf_key",
    "levels_passed_to_pertpy",
    "root",
    "edges",
    "node_names",
    "ancestor_matrix_shape",
    "fingerprint",
    "reference_path",
    "declared_reference_path",
    "reference_nodes",
    "derived_leaf_effects_are_propagated",
}


def _json_sha256(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _axis(values: Sequence[Any], *, label: str) -> tuple[list[str], str]:
    result = []
    for position, value in enumerate(values):
        if not isinstance(value, str) or not value or value != value.strip():
            raise ValueError(
                f"Composition-model {label} at position {position} must be a canonical nonblank string."
            )
        result.append(value)
    if not result or len(result) != len(set(result)):
        raise ValueError(f"Composition-model {label} must be nonempty and unique.")
    return result, _json_sha256(result)


def _canonical_json(value: Any, *, label: str) -> Any:
    try:
        payload = copy.deepcopy(value)
        json.dumps(payload, ensure_ascii=False, allow_nan=False)
    except (TypeError, ValueError) as error:
        raise TypeError(f"Composition-model {label} must be strict JSON data.") from error
    return payload


def _table_identity(table: Any, *, method: str, cell_types: list[str], hierarchy: Any) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    if not isinstance(table, pd.DataFrame):
        raise TypeError("Composition-model result table must be a pandas DataFrame.")
    expected = SCCODA_RESULT_COLUMNS if method == "sccoda" else TASCCODA_RESULT_COLUMNS
    if tuple(table.columns) != expected:
        raise ValueError(f"{method} result table columns do not match the current canonical schema.")
    if not isinstance(table.index, pd.RangeIndex) or table.index.start != 0 or table.index.step != 1:
        raise ValueError("Composition-model result table must use a zero-based RangeIndex.")
    if table.empty:
        raise ValueError("Composition-model result table cannot be empty.")

    if method == "sccoda":
        if table["cell_type"].tolist() != cell_types:
            raise ValueError("scCODA result rows do not match the cell-type axis.")
        reference = table["reference_cell_type"].tolist()
        if len(set(reference)) != 1 or reference[0] not in cell_types:
            raise ValueError("scCODA compositional reference is invalid.")
        expected_reference = [cell_type == reference[0] for cell_type in cell_types]
        if table["reference_constraint"].tolist() != expected_reference:
            raise ValueError("scCODA reference-constraint rows do not match its declared reference.")
        for column in ("credible_effect", "reference_constraint"):
            if not pd.api.types.is_bool_dtype(table[column].dtype):
                raise TypeError(f"scCODA result {column!r} must use a boolean dtype.")
        for column in (
            "inclusion_probability",
            "estimated_fdr",
            "inclusion_probability_threshold",
            "realized_expected_fdr",
        ):
            values = table[column].to_numpy(dtype=float)
            if bool(((values < 0.0) | (values > 1.0)).any()):
                raise ValueError(f"scCODA result {column!r} must lie in [0, 1].")
    else:
        direct = table.loc[table["effect_scope"] == "hierarchy_node"]
        leaves = table.loc[table["effect_scope"] == "derived_leaf"]
        if len(direct) + len(leaves) != len(table):
            raise ValueError("tascCODA effect_scope contains unsupported values.")
        if direct["effect_name"].tolist() != hierarchy["node_names"]:
            raise ValueError("tascCODA direct rows do not match the hierarchy-node axis.")
        if leaves["effect_name"].tolist() != cell_types:
            raise ValueError("tascCODA derived rows do not match the cell-type axis.")

    digest = hashlib.sha256(f"openbio-singlecell/{method}-result-table/v1\0".encode("ascii"))
    for row in table.itertuples(index=False, name=None):
        for value in row:
            missing = pd.isna(value)
            if isinstance(missing, (bool, np.bool_)) and bool(missing):
                digest.update(b"n")
            elif isinstance(value, (bool, np.bool_)):
                digest.update(b"b1" if bool(value) else b"b0")
            elif isinstance(value, str):
                if not value or value != value.strip():
                    raise ValueError("Composition-model result strings must be canonical and nonblank.")
                encoded = value.encode("utf-8")
                digest.update(b"s" + len(encoded).to_bytes(8, "little") + encoded)
            elif isinstance(value, (int, float, np.integer, np.floating)):
                numeric = float(value)
                if not math.isfinite(numeric):
                    raise ValueError("Composition-model result numeric values must be finite.")
                digest.update(b"f" + np.asarray([numeric], dtype="<f8").tobytes())
            else:
                raise TypeError(
                    f"Composition-model result contains unsupported value type {type(value).__name__}."
                )
    return {"rows": int(len(table)), "columns": list(expected), "content_sha256": digest.hexdigest()}


def _canonical_arrays(
    values: Mapping[str, Any],
    *,
    expected: Mapping[str, str | None],
    draw_count: int | None,
    axis_sizes: Mapping[str, int],
    add_chain_axis: bool,
    copy_arrays: bool,
    label: str,
) -> tuple[dict[str, Any], int]:
    import numpy as np

    if not isinstance(values, Mapping) or set(values) != set(expected):
        raise ValueError(f"Composition-model {label} variables do not match the current exact schema.")
    result = {}
    resolved_draws = draw_count
    for name, axis in expected.items():
        try:
            array = np.asarray(values[name], dtype=float)
        except (TypeError, ValueError) as error:
            raise TypeError(f"Composition-model {label} variable {name!r} must be numeric.") from error
        expected_tail = () if axis is None else (axis_sizes[axis],)
        if add_chain_axis:
            if array.ndim != len(expected_tail) + 1 or array.shape[1:] != expected_tail:
                raise ValueError(f"Composition-model {label} variable {name!r} has invalid dimensions.")
            array = array[None, ...]
        elif array.ndim != len(expected_tail) + 2 or array.shape[0] != 1 or array.shape[2:] != expected_tail:
            raise ValueError(f"Composition-model {label} variable {name!r} has invalid chain/draw dimensions.")
        current_draws = int(array.shape[1])
        if current_draws < 1 or (resolved_draws is not None and current_draws != resolved_draws):
            raise ValueError(f"Composition-model {label} variables must share a nonempty draw axis.")
        resolved_draws = current_draws
        if not bool(np.isfinite(array).all()):
            raise ValueError(f"Composition-model {label} variable {name!r} contains non-finite values.")
        result[name] = (
            np.array(array, dtype=float, copy=True, order="C")
            if copy_arrays
            else np.asarray(array, dtype=float, order="C")
        )
    assert resolved_draws is not None
    return result, resolved_draws


def _sample_stat_arrays(
    values: Mapping[str, Any],
    *,
    draw_count: int,
    add_chain_axis: bool,
    copy_arrays: bool,
) -> dict[str, Any]:
    import numpy as np

    if not isinstance(values, Mapping) or set(values) != set(_SAMPLE_STATS):
        raise ValueError("Composition-model sample-stat variables do not match the current exact schema.")
    result = {}
    for name in _SAMPLE_STATS:
        try:
            array = np.asarray(values[name], dtype=float)
        except (TypeError, ValueError) as error:
            raise TypeError(f"Composition-model sample statistic {name!r} must be numeric.") from error
        expected_shape = (draw_count,) if add_chain_axis else (1, draw_count)
        if array.shape != expected_shape:
            raise ValueError(
                f"Composition-model sample statistic {name!r} must contain one value per posterior draw."
            )
        if not bool(np.isfinite(array).all()):
            raise ValueError(f"Composition-model sample statistic {name!r} contains non-finite values.")
        if name in {"num_steps", "step_size"} and bool((array <= 0.0).any()):
            raise ValueError(f"Composition-model sample statistic {name!r} must be positive.")
        canonical = array[None, ...] if add_chain_axis else array
        result[name] = (
            np.array(canonical, dtype=float, copy=True, order="C")
            if copy_arrays
            else np.asarray(canonical, dtype=float, order="C")
        )
    return result


def _array_identities(values: Mapping[str, Any], *, namespace: str) -> dict[str, Any]:
    import numpy as np

    result = {}
    for name, value in values.items():
        digest = hashlib.sha256(f"openbio-singlecell/{namespace}/{name}/v1\0".encode("ascii"))
        digest.update(json.dumps(list(value.shape), separators=(",", ":")).encode("ascii"))
        digest.update(np.ascontiguousarray(value, dtype="<f8").tobytes(order="C"))
        result[name] = {"shape": list(value.shape), "content_sha256": digest.hexdigest()}
    return result


def _validate_model_metadata(value: Any, *, method: str, draw_count: int) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != _MODEL_METADATA_FIELDS:
        raise ValueError("Composition-model metadata has an invalid exact schema.")
    result = _canonical_json(dict(value), label="model metadata")
    for name in (
        "condition_key",
        "reference_condition",
        "comparison_condition",
        "reference_cell_type",
        "annotation_key",
    ):
        if not isinstance(result[name], str) or not result[name] or result[name] != result[name].strip():
            raise ValueError(f"Composition-model metadata {name!r} must be canonical nonblank text.")
    if result["reference_condition"] == result["comparison_condition"]:
        raise ValueError("Composition-model metadata Conditions must differ.")
    if result["annotation_status"] not in {"provisional", "curated"}:
        raise ValueError("Composition-model annotation status is invalid.")
    model = result["model"]
    diagnostics = result["diagnostics"]
    selection = result["selection"]
    software_versions = result["software_versions"]
    if not all(isinstance(value, dict) and value for value in (model, diagnostics, selection, software_versions)):
        raise ValueError("Composition-model model, diagnostics, selection, and software metadata must be objects.")
    if software_versions.get("pertpy") != "1.3.0" or not isinstance(software_versions.get("arviz"), str):
        raise ValueError("Composition-model software metadata must identify Pertpy 1.3.0 and ArviZ.")
    if model.get("chain_count") != 1 or model.get("posterior_draws") != draw_count:
        raise ValueError("Composition-model metadata does not match the retained chain/draw axes.")
    if model.get("method") != method or selection.get("method") != {
        "sccoda": "posterior_expected_fdr",
        "tasccoda": "tree_adaptive_spike_and_slab_lasso",
    }[method]:
        raise ValueError("Composition-model method metadata is inconsistent with its typed result.")
    if diagnostics.get("rhat_available") is not False or diagnostics.get("divergences_available") is not False:
        raise ValueError(f"{method} retained diagnostics must disclose unavailable R-hat and divergences.")
    return result


def _validate_hierarchy(value: Any, *, method: str) -> tuple[dict[str, Any] | None, list[str]]:
    if method == "sccoda":
        if value is not None:
            raise ValueError("scCODA result cannot carry tascCODA hierarchy metadata.")
        return None, []
    if not isinstance(value, Mapping) or set(value) != _HIERARCHY_FIELDS:
        raise ValueError("tascCODA hierarchy metadata has an invalid exact schema.")
    result = _canonical_json(dict(value), label="hierarchy metadata")
    nodes, _ = _axis(result["node_names"], label="hierarchy-node axis")
    if result["derived_leaf_effects_are_propagated"] is not True:
        raise ValueError("tascCODA hierarchy must disclose propagated derived leaf effects.")
    if not set(result["reference_nodes"]).issubset(nodes):
        raise ValueError("tascCODA reference path contains unknown hierarchy nodes.")
    return result, nodes


class CompositionModelResult:
    """Immutable-by-interface posterior diagnostic evidence for scCODA or tascCODA."""

    __slots__ = (
        "_cell_types",
        "_hierarchy",
        "_metadata",
        "_method",
        "_model_metadata",
        "_posterior",
        "_sample_stats",
        "_table",
    )
    artifact_type = COMPOSITION_MODEL_ARTIFACT_TYPE

    @classmethod
    def _from_owned(cls, **payload: Any) -> CompositionModelResult:
        result = object.__new__(cls)
        for key, value in payload.items():
            object.__setattr__(result, f"_{key}", value)
        return result

    def __setattr__(self, _name: str, _value: Any) -> None:
        raise AttributeError("CompositionModelResult is immutable; create a new validated artifact instead.")

    @property
    def fingerprint(self) -> str:
        return str(self._metadata["artifact_fingerprint_sha256"])

    @property
    def method(self) -> str:
        return self._method

    @property
    def table(self) -> Any:
        return self._table.copy(deep=True)

    @property
    def posterior(self) -> dict[str, Any]:
        return {name: value.copy() for name, value in self._posterior.items()}

    @property
    def sample_stats(self) -> dict[str, Any]:
        return {name: value.copy() for name, value in self._sample_stats.items()}

    @property
    def cell_types(self) -> list[str]:
        return list(self._cell_types)

    @property
    def hierarchy(self) -> dict[str, Any] | None:
        return copy.deepcopy(self._hierarchy)

    @property
    def model_metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._model_metadata)

    @property
    def metadata(self) -> dict[str, Any]:
        return copy.deepcopy(self._metadata)

    def portable(self) -> dict[str, Any]:
        return {
            "artifact_type": COMPOSITION_MODEL_ARTIFACT_TYPE,
            "method": self._method,
            "table": self._table.copy(deep=True),
            "posterior": self.posterior,
            "sample_stats": self.sample_stats,
            "cell_types": list(self._cell_types),
            "hierarchy": copy.deepcopy(self._hierarchy),
            "model_metadata": copy.deepcopy(self._model_metadata),
            "metadata": copy.deepcopy(self._metadata),
        }

    def inference_data(self) -> Any:
        try:
            import arviz as az
        except ImportError as error:
            raise RuntimeError("Composition-model InferenceData requires the installed composition dependencies.") from error
        dims = {
            name: ([] if axis is None else [axis])
            for name, axis in _POSTERIOR_VARIABLES[self._method].items()
        }
        coords = {"cell_type": list(self._cell_types)}
        if self._method == "tasccoda":
            coords["hierarchy_node"] = list(self._hierarchy["node_names"])
        return az.from_dict(
            {"posterior": self.posterior, "sample_stats": self.sample_stats},
            sample_dims=["chain", "draw"],
            coords=coords,
            dims=dims,
        )


def build_composition_model_result(
    *,
    table: Any,
    method: str,
    posterior: Mapping[str, Any],
    sample_stats: Mapping[str, Any],
    cell_types: Sequence[Any],
    hierarchy: Mapping[str, Any] | None,
    model_metadata: Mapping[str, Any],
) -> CompositionModelResult:
    if method not in _POSTERIOR_VARIABLES:
        raise ValueError("Composition-model method must be 'sccoda' or 'tasccoda'.")
    cell_type_axis, cell_type_hash = _axis(cell_types, label="cell-type axis")
    hierarchy_owned, hierarchy_nodes = _validate_hierarchy(hierarchy, method=method)
    axis_sizes = {"cell_type": len(cell_type_axis), "hierarchy_node": len(hierarchy_nodes)}
    posterior_owned, draw_count = _canonical_arrays(
        posterior,
        expected=_POSTERIOR_VARIABLES[method],
        draw_count=None,
        axis_sizes=axis_sizes,
        add_chain_axis=True,
        copy_arrays=True,
        label="posterior",
    )
    sample_stats_owned = _sample_stat_arrays(
        sample_stats,
        draw_count=draw_count,
        add_chain_axis=True,
        copy_arrays=True,
    )
    model_metadata_owned = _validate_model_metadata(model_metadata, method=method, draw_count=draw_count)
    if model_metadata_owned["reference_cell_type"] not in cell_type_axis:
        raise ValueError("Composition-model reference cell type is absent from its cell-type axis.")
    reference_index = cell_type_axis.index(model_metadata_owned["reference_cell_type"])
    if method == "sccoda" and not bool((posterior_owned["condition_effect"][:, :, reference_index] == 0.0).all()):
        raise ValueError("scCODA posterior violates the reference-cell zero-effect constraint.")
    if method == "tasccoda":
        reference_indices = [hierarchy_nodes.index(name) for name in hierarchy_owned["reference_nodes"]]
        if reference_indices and not bool(
            (posterior_owned["hierarchy_node_effect"][:, :, reference_indices] == 0.0).all()
        ):
            raise ValueError("tascCODA posterior violates the reference-path zero-effect constraint.")
    table_owned = table.copy(deep=True)
    table_identity = _table_identity(
        table_owned,
        method=method,
        cell_types=cell_type_axis,
        hierarchy=hierarchy_owned,
    )
    posterior_identity = _array_identities(posterior_owned, namespace=f"{method}-posterior")
    sample_stats_identity = _array_identities(sample_stats_owned, namespace=f"{method}-sample-stats")
    metadata: dict[str, Any] = {
        "schema_version": COMPOSITION_MODEL_ARTIFACT_SCHEMA_VERSION,
        "artifact_type": COMPOSITION_MODEL_ARTIFACT_TYPE,
        "producer_node_id": (
            "OpenBioSingleCellSccodaDifferentialComposition"
            if method == "sccoda"
            else "OpenBioSingleCellTasccodaDifferentialComposition"
        ),
        "producer_schema": f"openbio-singlecell/{method}-composition-model-result/v1",
        "method": method,
        "chain_count": 1,
        "draw_count": draw_count,
        "cell_type_count": len(cell_type_axis),
        "cell_type_axis_sha256": cell_type_hash,
        "hierarchy_node_count": len(hierarchy_nodes),
        "hierarchy_node_axis_sha256": _json_sha256(hierarchy_nodes),
        "table_identity": table_identity,
        "posterior_identity": posterior_identity,
        "sample_stats_identity": sample_stats_identity,
        "hierarchy_sha256": _json_sha256(hierarchy_owned),
        "model_metadata_sha256": _json_sha256(model_metadata_owned),
    }
    metadata["artifact_fingerprint_sha256"] = _json_sha256(metadata)
    result = CompositionModelResult._from_owned(
        method=method,
        table=table_owned,
        posterior=posterior_owned,
        sample_stats=sample_stats_owned,
        cell_types=cell_type_axis,
        hierarchy=hierarchy_owned,
        model_metadata=model_metadata_owned,
        metadata=metadata,
    )
    validate_composition_model_result(result, copy_result=False)
    return result


def validate_composition_model_result(
    result: Any,
    *,
    exact_type: bool = True,
    copy_result: bool = True,
) -> tuple[Any, dict[str, Any], dict[str, Any], list[str], Any, dict[str, Any], dict[str, Any]]:
    if exact_type:
        if type(result) is not CompositionModelResult:
            raise TypeError("Expected an exact OPENBIO_COMPOSITION_MODEL_RESULT artifact.")
        payload = {
            "artifact_type": COMPOSITION_MODEL_ARTIFACT_TYPE,
            "method": object.__getattribute__(result, "_method"),
            "table": object.__getattribute__(result, "_table"),
            "posterior": object.__getattribute__(result, "_posterior"),
            "sample_stats": object.__getattribute__(result, "_sample_stats"),
            "cell_types": object.__getattribute__(result, "_cell_types"),
            "hierarchy": object.__getattribute__(result, "_hierarchy"),
            "model_metadata": object.__getattribute__(result, "_model_metadata"),
            "metadata": object.__getattribute__(result, "_metadata"),
        }
    elif isinstance(result, Mapping):
        payload = dict(result)
    else:
        raise TypeError("Portable composition-model result must be a mapping.")
    expected_payload = {
        "artifact_type",
        "method",
        "table",
        "posterior",
        "sample_stats",
        "cell_types",
        "hierarchy",
        "model_metadata",
        "metadata",
    }
    if set(payload) != expected_payload or payload["artifact_type"] != COMPOSITION_MODEL_ARTIFACT_TYPE:
        raise ValueError("Composition-model portable result schema is invalid.")
    method = payload["method"]
    if method not in _POSTERIOR_VARIABLES:
        raise ValueError("Composition-model method identity is invalid.")
    metadata = payload["metadata"]
    expected_metadata = {
        "schema_version",
        "artifact_type",
        "producer_node_id",
        "producer_schema",
        "method",
        "chain_count",
        "draw_count",
        "cell_type_count",
        "cell_type_axis_sha256",
        "hierarchy_node_count",
        "hierarchy_node_axis_sha256",
        "table_identity",
        "posterior_identity",
        "sample_stats_identity",
        "hierarchy_sha256",
        "model_metadata_sha256",
        "artifact_fingerprint_sha256",
    }
    if not isinstance(metadata, Mapping) or set(metadata) != expected_metadata:
        raise ValueError("Composition-model artifact metadata schema is invalid.")
    expected_producer = (
        "OpenBioSingleCellSccodaDifferentialComposition"
        if method == "sccoda"
        else "OpenBioSingleCellTasccodaDifferentialComposition"
    )
    if (
        metadata["schema_version"] != COMPOSITION_MODEL_ARTIFACT_SCHEMA_VERSION
        or metadata["artifact_type"] != COMPOSITION_MODEL_ARTIFACT_TYPE
        or metadata["method"] != method
        or metadata["producer_node_id"] != expected_producer
        or metadata["producer_schema"] != f"openbio-singlecell/{method}-composition-model-result/v1"
        or metadata["chain_count"] != 1
    ):
        raise ValueError("Composition-model artifact producer or schema identity is invalid.")
    metadata_without_fingerprint = dict(metadata)
    fingerprint = metadata_without_fingerprint.pop("artifact_fingerprint_sha256")
    if fingerprint != _json_sha256(metadata_without_fingerprint):
        raise ValueError("Composition-model artifact metadata fingerprint is invalid.")
    if exact_type and result.fingerprint != fingerprint:
        raise ValueError("Composition-model artifact fingerprint property is inconsistent.")

    cell_types, cell_type_hash = _axis(payload["cell_types"], label="cell-type axis")
    hierarchy, hierarchy_nodes = _validate_hierarchy(payload["hierarchy"], method=method)
    if (
        metadata["cell_type_count"] != len(cell_types)
        or metadata["cell_type_axis_sha256"] != cell_type_hash
        or metadata["hierarchy_node_count"] != len(hierarchy_nodes)
        or metadata["hierarchy_node_axis_sha256"] != _json_sha256(hierarchy_nodes)
        or metadata["hierarchy_sha256"] != _json_sha256(hierarchy)
    ):
        raise ValueError("Composition-model artifact axes or hierarchy failed current-content checks.")
    draw_count = metadata["draw_count"]
    if not isinstance(draw_count, int) or isinstance(draw_count, bool) or draw_count < 1:
        raise ValueError("Composition-model draw count is invalid.")
    axis_sizes = {"cell_type": len(cell_types), "hierarchy_node": len(hierarchy_nodes)}
    posterior, observed_draws = _canonical_arrays(
        payload["posterior"],
        expected=_POSTERIOR_VARIABLES[method],
        draw_count=draw_count,
        axis_sizes=axis_sizes,
        add_chain_axis=False,
        copy_arrays=False,
        label="posterior",
    )
    if observed_draws != draw_count or _array_identities(
        posterior, namespace=f"{method}-posterior"
    ) != metadata["posterior_identity"]:
        raise ValueError("Composition-model posterior failed its current-content fingerprint check.")
    sample_stats = _sample_stat_arrays(
        payload["sample_stats"],
        draw_count=draw_count,
        add_chain_axis=False,
        copy_arrays=False,
    )
    if _array_identities(sample_stats, namespace=f"{method}-sample-stats") != metadata["sample_stats_identity"]:
        raise ValueError("Composition-model sample statistics failed their current-content fingerprint check.")
    model_metadata = _validate_model_metadata(payload["model_metadata"], method=method, draw_count=draw_count)
    if metadata["model_metadata_sha256"] != _json_sha256(model_metadata):
        raise ValueError("Composition-model metadata failed its current-content fingerprint check.")
    reference = model_metadata["reference_cell_type"]
    if reference not in cell_types:
        raise ValueError("Composition-model reference cell type is absent from its axis.")
    if method == "sccoda" and not bool((posterior["condition_effect"][:, :, cell_types.index(reference)] == 0).all()):
        raise ValueError("scCODA posterior violates the reference-cell zero-effect constraint.")
    if method == "tasccoda":
        reference_indices = [hierarchy_nodes.index(name) for name in hierarchy["reference_nodes"]]
        if reference_indices and not bool((posterior["hierarchy_node_effect"][:, :, reference_indices] == 0).all()):
            raise ValueError("tascCODA posterior violates the reference-path zero-effect constraint.")
    table = payload["table"]
    if _table_identity(table, method=method, cell_types=cell_types, hierarchy=hierarchy) != metadata["table_identity"]:
        raise ValueError("Composition-model result table failed its current-content fingerprint check.")
    return (
        table.copy(deep=True) if copy_result else table,
        {name: value.copy() for name, value in posterior.items()} if copy_result else posterior,
        {name: value.copy() for name, value in sample_stats.items()} if copy_result else sample_stats,
        cell_types,
        copy.deepcopy(hierarchy),
        copy.deepcopy(model_metadata),
        copy.deepcopy(dict(metadata)),
    )


__all__ = [
    "COMPOSITION_MODEL_ARTIFACT_SCHEMA_VERSION",
    "COMPOSITION_MODEL_ARTIFACT_TYPE",
    "CompositionModelResult",
    "SCCODA_RESULT_COLUMNS",
    "TASCCODA_RESULT_COLUMNS",
    "build_composition_model_result",
    "validate_composition_model_result",
]
