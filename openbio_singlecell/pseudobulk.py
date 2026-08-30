from __future__ import annotations

import inspect
import textwrap
from typing import TYPE_CHECKING, Any

from .analysis_reporting import _package_version, _plain_json, collect_software_versions

if TYPE_CHECKING:
    from anndata import AnnData


PSEUDOBULK_ARTIFACT_TYPE = "OPENBIO_SINGLE_CELL_PSEUDOBULK"
PSEUDOBULK_SCHEMA_VERSION = 1
PSEUDOBULK_QC_COLUMNS = ["openbio_n_cells", "openbio_total_counts"]


def _standalone_artifact_fingerprint(adata, metadata):
    """Hash canonical artifact metadata, identities, profile metadata, and integer counts."""
    import hashlib
    import json

    import numpy as np
    from scipy import sparse

    metadata_without_fingerprint = dict(metadata)
    metadata_without_fingerprint.pop("artifact_fingerprint", None)
    header = {
        "schema": "openbio-singlecell/pseudobulk-artifact/v1",
        "metadata": _plain_json(metadata_without_fingerprint),
        "obs_columns": [str(value) for value in adata.obs.columns.tolist()],
        "obs_names": [str(value) for value in adata.obs_names.tolist()],
        "obs_rows": _plain_json(adata.obs.to_numpy(dtype=object).tolist()),
        "var_columns": [str(value) for value in adata.var.columns.tolist()],
        "var_names": [str(value) for value in adata.var_names.tolist()],
        "var_rows": _plain_json(adata.var.to_numpy(dtype=object).tolist()),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
    }
    encoder = json.JSONEncoder(ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)
    digest = hashlib.sha256()
    for chunk in encoder.iterencode(header):
        digest.update(chunk.encode("utf-8"))
    matrix = adata.X.tocsr() if sparse.issparse(adata.X) else np.asarray(adata.X)
    for row_index in range(int(adata.n_obs)):
        row = matrix.getrow(row_index).toarray().ravel() if sparse.issparse(matrix) else matrix[row_index]
        canonical = np.asarray(row, dtype="<i8")
        digest.update(int(canonical.size).to_bytes(8, "little", signed=False))
        digest.update(canonical.tobytes(order="C"))
    return digest.hexdigest()


class _StandalonePseudobulkArtifact:
    """Immutable-by-interface Sample-by-population raw-count artifact."""

    __slots__ = ("_adata", "_metadata")
    artifact_type = "OPENBIO_SINGLE_CELL_PSEUDOBULK"

    def __init__(self, adata, metadata):
        import copy

        object.__setattr__(self, "_adata", adata.copy())
        object.__setattr__(self, "_metadata", copy.deepcopy(dict(metadata)))

    @classmethod
    def _from_owned(cls, adata, metadata):
        """Build from worker-owned state without a defensive full-matrix copy."""

        import copy

        artifact = object.__new__(cls)
        object.__setattr__(artifact, "_adata", adata)
        object.__setattr__(artifact, "_metadata", copy.deepcopy(dict(metadata)))
        return artifact

    def __setattr__(self, name, value):
        raise AttributeError("PseudobulkArtifact is immutable; create a new validated artifact instead.")

    @property
    def fingerprint(self):
        return str(self._metadata["artifact_fingerprint"])

    @property
    def metadata(self):
        import copy

        return copy.deepcopy(self._metadata)

    def to_adata(self):
        """Return a defensive copy of the private pseudobulk AnnData."""
        return self._adata.copy()

    def _owned_adata(self):
        """Return worker-private state for trusted internal codec/analysis paths."""

        return self._adata


PseudobulkArtifact = _StandalonePseudobulkArtifact


def _standalone_validate_pseudobulk_artifact(artifact, *, copy_result=True):
    """Validate the complete artifact; public callers receive a defensive data copy by default."""
    import json
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from scipy import sparse

    operation = "Pseudobulk artifact validation"
    if getattr(artifact, "artifact_type", None) != "OPENBIO_SINGLE_CELL_PSEUDOBULK":
        raise TypeError(f"{operation} requires an OPENBIO_SINGLE_CELL_PSEUDOBULK artifact.")
    if not callable(getattr(artifact, "to_adata", None)):
        raise TypeError(f"{operation} artifact does not expose a defensive to_adata() copy.")
    metadata = getattr(artifact, "metadata", None)
    if not isinstance(metadata, Mapping):
        raise TypeError(f"{operation} metadata must be a mapping.")
    metadata = dict(metadata)
    required_metadata = {
        "schema_version",
        "artifact_type",
        "producer_node_id",
        "artifact_fingerprint",
        "aggregation_mode",
        "role_keys",
        "declared_covariates",
        "annotation_status",
        "annotation_provenance_operation",
        "inference_status",
        "declarations",
        "role_aliases",
        "formal_interpretation_invalid",
        "formal_interpretation_invalid_reasons",
        "count_source",
        "profile_filters",
        "input_dimensions",
        "profile_records",
        "missing_combinations",
        "profile_order",
        "obs_columns",
        "backend",
    }
    missing = sorted(required_metadata - set(metadata))
    unknown = sorted(set(metadata) - required_metadata)
    if missing or unknown:
        raise ValueError(f"{operation} metadata schema mismatch; missing={missing}, unknown={unknown}.")
    if metadata["schema_version"] != 2:
        raise ValueError(f"{operation} requires schema_version=2.")
    if metadata["artifact_type"] != "OPENBIO_SINGLE_CELL_PSEUDOBULK":
        raise ValueError(f"{operation} metadata has the wrong artifact type.")
    if metadata["producer_node_id"] != "OpenBioSingleCellPseudobulk":
        raise ValueError(f"{operation} metadata has the wrong producer node.")
    if metadata["aggregation_mode"] != "sum":
        raise ValueError(f"{operation} requires summed values.")
    fingerprint = metadata["artifact_fingerprint"]
    if not (
        isinstance(fingerprint, str)
        and len(fingerprint) == 64
        and all(character in "0123456789abcdef" for character in fingerprint)
    ):
        raise ValueError(f"{operation} fingerprint is invalid.")
    if getattr(artifact, "fingerprint", None) != fingerprint:
        raise ValueError(f"{operation} fingerprint property disagrees with metadata.")
    owned_accessor = getattr(artifact, "_owned_adata", None)
    adata = artifact.to_adata() if copy_result or not callable(owned_accessor) else owned_accessor()
    if not isinstance(adata, AnnData):
        raise TypeError(f"{operation} to_adata() must return an AnnData object.")
    if getattr(adata, "isbacked", False):
        raise ValueError(f"{operation} requires an in-memory artifact.")
    if adata.n_obs < 1 or adata.n_vars < 1:
        raise ValueError(f"{operation} requires nonempty profiles and genes.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError(f"{operation} requires unique profile and gene identifiers.")
    if any(not isinstance(value, str) or not value for value in adata.var_names.tolist()):
        raise ValueError(f"{operation} requires nonblank string gene identifiers.")
    layer_keys = list(adata.layers.keys())
    has_only_anndata_none_alias = not layer_keys or (
        layer_keys == [None]
        and getattr(adata.layers[None], "shape", None) == getattr(adata.X, "shape", None)
        and (
            (sparse.csr_matrix(adata.layers[None]) - sparse.csr_matrix(adata.X)).nnz == 0
            if sparse.issparse(adata.layers[None]) or sparse.issparse(adata.X)
            else bool(np.array_equal(np.asarray(adata.layers[None]), np.asarray(adata.X)))
        )
    )
    if (
        adata.raw is not None
        or not has_only_anndata_none_alias
        or len(adata.obsm)
        or len(adata.varm)
        or len(adata.obsp)
        or len(adata.varp)
        or len(adata.uns)
    ):
        raise ValueError(f"{operation} must not carry unvalidated Raw/layer/multidimensional/uns state.")
    expected_columns = metadata["obs_columns"]
    if not isinstance(expected_columns, list) or list(adata.obs.columns) != expected_columns:
        raise ValueError(f"{operation} observation columns differ from immutable provenance.")
    roles = metadata["role_keys"]
    if not isinstance(roles, Mapping) or set(roles) != {
        "sample",
        "population",
        "condition",
        "technical_batch",
    }:
        raise ValueError(f"{operation} role_keys are malformed.")
    for role in ("sample", "population", "condition"):
        if not isinstance(roles[role], str) or not roles[role].strip() or roles[role] != roles[role].strip():
            raise ValueError(f"{operation} {role} role must be a nonblank canonical column name.")
    if roles["technical_batch"] is not None and (
        not isinstance(roles["technical_batch"], str)
        or not roles["technical_batch"].strip()
        or roles["technical_batch"] != roles["technical_batch"].strip()
    ):
        raise ValueError(f"{operation} technical_batch role must be null or a nonblank canonical column name.")
    required_obs = [
        roles["sample"],
        roles["population"],
        roles["condition"],
        "openbio_n_cells",
        "openbio_total_counts",
    ]
    if roles["technical_batch"]:
        required_obs.append(roles["technical_batch"])
    declared_covariates = metadata["declared_covariates"]
    if not isinstance(declared_covariates, Mapping) or set(declared_covariates) != {
        "categorical",
        "continuous",
    }:
        raise ValueError(f"{operation} declared_covariates are malformed.")
    for kind in ("categorical", "continuous"):
        keys = declared_covariates[kind]
        if (
            not isinstance(keys, list)
            or any(not isinstance(key, str) or not key.strip() or key != key.strip() for key in keys)
            or len(keys) != len(set(keys))
        ):
            raise ValueError(f"{operation} declared {kind} covariates are malformed.")
    required_obs.extend(declared_covariates["categorical"])
    required_obs.extend(declared_covariates["continuous"])
    if any(column not in adata.obs for column in required_obs):
        raise ValueError(f"{operation} required observation roles are missing.")
    role_entries = [
        (roles["sample"], "sample"),
        (roles["population"], "population"),
        (roles["condition"], "condition"),
    ]
    if roles["technical_batch"]:
        role_entries.append((roles["technical_batch"], "technical_batch"))
    role_entries.extend((key, "categorical_covariate") for key in declared_covariates["categorical"])
    role_entries.extend((key, "continuous_covariate") for key in declared_covariates["continuous"])
    if any(key in {"openbio_n_cells", "openbio_total_counts"} for key, _ in role_entries):
        raise ValueError(f"{operation} declared roles collide with reserved profile-QC columns.")
    roles_by_key = {}
    for key, role in role_entries:
        roles_by_key.setdefault(key, []).append(role)
    expected_role_aliases = {key: value for key, value in roles_by_key.items() if len(value) > 1}
    if metadata["role_aliases"] != expected_role_aliases:
        raise ValueError(f"{operation} role_aliases disagree with declared metadata roles.")
    annotation_status = metadata["annotation_status"]
    annotation_provenance_operation = metadata["annotation_provenance_operation"]
    inference_status = metadata["inference_status"]
    if annotation_status not in {"curated", "provisional", "unknown"}:
        raise ValueError(f"{operation} annotation_status is invalid.")
    if annotation_provenance_operation not in {
        None,
        "map_cluster_annotations",
        "celltypist_annotation",
    }:
        raise ValueError(f"{operation} annotation provenance operation is invalid.")
    if (annotation_status == "unknown") != (annotation_provenance_operation is None):
        raise ValueError(f"{operation} annotation status/provenance operation are inconsistent.")
    allowed_inference = {
        "formal_sample_level_condition_inference_ready",
        "formal_interpretation_invalid",
        "exploratory_sample_level_condition_evidence",
    }
    if inference_status not in allowed_inference:
        raise ValueError(f"{operation} inference_status is invalid.")
    declarations = metadata["declarations"]
    expected_declaration_keys = {
        "count_source",
        "feature_completeness",
        "sample_identity",
        "population_annotation",
    }
    if (
        not isinstance(declarations, Mapping)
        or set(declarations) != expected_declaration_keys
        or any(not isinstance(value, str) or not value for value in declarations.values())
    ):
        raise ValueError(f"{operation} caller declarations are malformed.")
    formal_invalid = metadata["formal_interpretation_invalid"]
    formal_invalid_reasons = metadata["formal_interpretation_invalid_reasons"]
    if not isinstance(formal_invalid, bool) or (
        not isinstance(formal_invalid_reasons, list)
        or any(not isinstance(reason, str) or not reason for reason in formal_invalid_reasons)
        or len(formal_invalid_reasons) != len(set(formal_invalid_reasons))
    ):
        raise ValueError(f"{operation} formal-interpretation disclosure is malformed.")
    if formal_invalid != bool(formal_invalid_reasons):
        raise ValueError(f"{operation} formal-interpretation flag disagrees with its reasons.")
    expected_status = (
        "exploratory_sample_level_condition_evidence"
        if "exploratory_mode_requested" in formal_invalid_reasons
        else "formal_interpretation_invalid"
        if formal_invalid
        else "formal_sample_level_condition_inference_ready"
    )
    if inference_status != expected_status:
        raise ValueError(f"{operation} inference status disagrees with interpretation disclosures.")
    count_source = metadata["count_source"]
    if not isinstance(count_source, Mapping) or set(count_source) != {
        "kind",
        "layer_name",
        "label",
        "feature_axis",
    }:
        raise ValueError(f"{operation} count_source provenance is malformed.")
    source_kind = count_source["kind"]
    layer_name = count_source["layer_name"]
    if source_kind not in {"raw", "X", "layer"}:
        raise ValueError(f"{operation} count_source kind is invalid.")
    if source_kind == "layer":
        if not isinstance(layer_name, str) or not layer_name.strip() or layer_name != layer_name.strip():
            raise ValueError(f"{operation} layer count_source requires a canonical layer name.")
        expected_source_label = f"layers[{layer_name!r}]"
    else:
        if layer_name is not None:
            raise ValueError(f"{operation} non-layer count_source cannot name a layer.")
        expected_source_label = "raw.X" if source_kind == "raw" else "X"
    if (
        count_source["label"] != expected_source_label
        or count_source["feature_axis"] != ("raw.var_names" if source_kind == "raw" else "adata.var_names")
    ):
        raise ValueError(f"{operation} count_source is inconsistent.")
    profile_filters = metadata["profile_filters"]
    if not isinstance(profile_filters, Mapping) or set(profile_filters) != {"min_cells", "min_counts"}:
        raise ValueError(f"{operation} profile_filters are malformed.")
    for name, minimum in (("min_cells", 1), ("min_counts", 0)):
        value = profile_filters[name]
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or int(value) < minimum:
            raise ValueError(f"{operation} profile filter {name!r} is invalid.")
    input_dimensions = metadata["input_dimensions"]
    if not isinstance(input_dimensions, Mapping) or set(input_dimensions) != {"cells", "genes", "samples"}:
        raise ValueError(f"{operation} input_dimensions are malformed.")
    for name in ("cells", "genes", "samples"):
        value = input_dimensions[name]
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or int(value) < 1:
            raise ValueError(f"{operation} input dimension {name!r} is invalid.")
    if int(input_dimensions["genes"]) != int(adata.n_vars):
        raise ValueError(f"{operation} gene dimension differs from immutable provenance.")
    backend = metadata["backend"]
    if not isinstance(backend, Mapping) or set(backend) != {
        "package",
        "version",
        "interface",
        "cell_qc_column",
        "synthetic_empty_profiles_removed",
    }:
        raise ValueError(f"{operation} backend provenance is malformed.")
    if (
        backend["package"] != "decoupler"
        or not isinstance(backend["version"], str)
        or not backend["version"].strip()
        or backend["interface"] != "decoupler.pp.pseudobulk"
        or backend["cell_qc_column"] not in {"psbulk_cells", "psbulk_n_cells"}
        or isinstance(backend["synthetic_empty_profiles_removed"], bool)
        or not isinstance(backend["synthetic_empty_profiles_removed"], (int, np.integer))
        or int(backend["synthetic_empty_profiles_removed"]) < 0
    ):
        raise ValueError(f"{operation} backend provenance is inconsistent.")
    for column in (roles["sample"], roles["population"], roles["condition"]):
        values = adata.obs[column]
        if values.isna().any() or any(not isinstance(value, str) or not value for value in values.tolist()):
            raise ValueError(f"{operation} column {column!r} contains invalid identifiers.")
    if bool(adata.obs.duplicated(subset=[roles["sample"], roles["population"]]).any()):
        raise ValueError(f"{operation} contains duplicate Sample-by-population profiles.")
    for key in ([roles["technical_batch"]] if roles["technical_batch"] else []) + declared_covariates["categorical"]:
        values = adata.obs[key]
        if values.isna().any() or any(not isinstance(value, str) or not value for value in values.tolist()):
            raise ValueError(f"{operation} categorical column {key!r} contains invalid labels.")
    for key in declared_covariates["continuous"]:
        if pd.api.types.is_bool_dtype(adata.obs[key].dtype):
            raise ValueError(f"{operation} continuous column {key!r} must be numeric and non-boolean.")
        try:
            continuous = pd.to_numeric(adata.obs[key], errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{operation} continuous column {key!r} must be numeric.") from exc
        if not bool(np.isfinite(continuous).all()):
            raise ValueError(f"{operation} continuous column {key!r} must be finite.")
    sample_level_columns = [roles["condition"]]
    if roles["technical_batch"]:
        sample_level_columns.append(roles["technical_batch"])
    sample_level_columns.extend(declared_covariates["categorical"])
    sample_level_columns.extend(declared_covariates["continuous"])
    for key in dict.fromkeys(sample_level_columns):
        if bool(adata.obs.groupby(roles["sample"], sort=False, observed=True)[key].nunique(dropna=False).gt(1).any()):
            raise ValueError(f"{operation} Sample-level column {key!r} is inconsistent within Sample.")

    counts = adata.X
    if getattr(counts, "shape", None) != (adata.n_obs, adata.n_vars):
        raise ValueError(f"{operation} count matrix shape is invalid.")
    if sparse.issparse(counts):
        values = np.asarray(counts.data)
    else:
        values = np.asarray(counts).ravel()
    if np.issubdtype(values.dtype, np.bool_) or not np.issubdtype(values.dtype, np.integer):
        raise TypeError(f"{operation} counts must have a non-boolean integer dtype.")
    if values.size and (bool((values < 0).any()) or int(values.max()) > int(np.iinfo(np.int64).max)):
        raise ValueError(f"{operation} counts must be non-negative and fit int64.")
    count_rows = counts.tocsr() if sparse.issparse(counts) else np.asarray(counts)
    totals = []
    for row_index in range(int(adata.n_obs)):
        row = count_rows.getrow(row_index).data if sparse.issparse(count_rows) else count_rows[row_index]
        total = int(sum(int(value) for value in np.asarray(row).ravel()))
        if not 0 < total <= int(np.iinfo(np.int64).max):
            raise ValueError(f"{operation} requires positive profile libraries within int64 range.")
        totals.append(total)
    n_cells_values = adata.obs["openbio_n_cells"].tolist()
    total_count_values = adata.obs["openbio_total_counts"].tolist()
    if any(
        isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) < 1
        for value in n_cells_values
    ) or any(
        isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) <= 0
        for value in total_count_values
    ):
        raise ValueError(f"{operation} profile QC columns must contain positive integers.")
    n_cells = [int(value) for value in n_cells_values]
    total_counts = [int(value) for value in total_count_values]
    if total_counts != totals:
        raise ValueError(f"{operation} profile QC columns disagree with the count matrix.")

    profile_records = metadata["profile_records"]
    if not isinstance(profile_records, list) or not profile_records:
        raise ValueError(f"{operation} profile_records must be a nonempty list.")
    record_fields = {
        "sample",
        "population",
        "condition",
        "n_cells",
        "total_counts",
        "retained",
        "removal_reasons",
    }
    represented_pairs = set()
    sample_conditions = {}
    retained_records = []
    represented_cell_count = 0
    for record in profile_records:
        if not isinstance(record, Mapping) or set(record) != record_fields:
            raise ValueError(f"{operation} profile record schema is malformed.")
        if any(not isinstance(record[key], str) or not record[key] for key in ("sample", "population", "condition")):
            raise ValueError(f"{operation} profile record identities are invalid.")
        pair = (record["sample"], record["population"])
        if pair in represented_pairs:
            raise ValueError(f"{operation} profile_records contain duplicate Sample-by-population identities.")
        represented_pairs.add(pair)
        prior_condition = sample_conditions.setdefault(record["sample"], record["condition"])
        if prior_condition != record["condition"]:
            raise ValueError(f"{operation} profile_records contain conflicting Sample-to-Condition mappings.")
        for key, minimum in (("n_cells", 1), ("total_counts", 1)):
            value = record[key]
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or int(value) < minimum:
                raise ValueError(f"{operation} profile record {key!r} is invalid.")
        reasons = []
        if int(record["n_cells"]) < int(profile_filters["min_cells"]):
            reasons.append("below_min_cells")
        if int(record["total_counts"]) < int(profile_filters["min_counts"]):
            reasons.append("below_min_counts")
        if not isinstance(record["retained"], bool) or record["removal_reasons"] != reasons:
            raise ValueError(f"{operation} profile record QC decision is inconsistent.")
        if record["retained"] != (not reasons):
            raise ValueError(f"{operation} profile record retained flag is inconsistent.")
        represented_cell_count += int(record["n_cells"])
        if record["retained"]:
            retained_records.append(record)
    if represented_cell_count != int(input_dimensions["cells"]):
        raise ValueError(f"{operation} contributing-cell accounting differs from immutable provenance.")
    if len(sample_conditions) != int(input_dimensions["samples"]):
        raise ValueError(f"{operation} Sample count differs from immutable provenance.")
    if len(retained_records) != int(adata.n_obs):
        raise ValueError(f"{operation} retained profile records do not match the count rows.")
    for row_index, record in enumerate(retained_records):
        if (
            adata.obs.iloc[row_index][roles["sample"]] != record["sample"]
            or adata.obs.iloc[row_index][roles["population"]] != record["population"]
            or adata.obs.iloc[row_index][roles["condition"]] != record["condition"]
            or n_cells[row_index] != int(record["n_cells"])
            or total_counts[row_index] != int(record["total_counts"])
        ):
            raise ValueError(f"{operation} retained profile rows differ from immutable provenance.")

    missing_combinations = metadata["missing_combinations"]
    if not isinstance(missing_combinations, list):
        raise ValueError(f"{operation} missing_combinations must be a list.")
    missing_pairs = set()
    for record in missing_combinations:
        if not isinstance(record, Mapping) or set(record) != {"sample", "population", "condition"}:
            raise ValueError(f"{operation} missing-combination schema is malformed.")
        sample = record["sample"]
        population = record["population"]
        condition = record["condition"]
        if (
            not isinstance(sample, str)
            or not sample
            or not isinstance(population, str)
            or not population
            or condition != sample_conditions.get(sample)
        ):
            raise ValueError(f"{operation} missing-combination identity is invalid.")
        pair = (sample, population)
        if pair in missing_pairs or pair in represented_pairs:
            raise ValueError(f"{operation} missing-combination identities overlap or repeat.")
        missing_pairs.add(pair)
    populations = {population for _, population in represented_pairs | missing_pairs}
    complete_grid = {(sample, population) for sample in sample_conditions for population in populations}
    if represented_pairs | missing_pairs != complete_grid:
        raise ValueError(f"{operation} represented/missing profiles do not form the disclosed complete grid.")
    if int(backend["synthetic_empty_profiles_removed"]) != len(missing_pairs):
        raise ValueError(f"{operation} synthetic-empty profile accounting is inconsistent.")
    strata_before = {}
    strata_after = {}
    for record in profile_records:
        stratum = (record["population"], record["condition"])
        strata_before[stratum] = strata_before.get(stratum, 0) + 1
        if record["retained"]:
            strata_after[stratum] = strata_after.get(stratum, 0) + 1
    profile_order = metadata["profile_order"]
    observed_order = [
        [sample, population]
        for sample, population in zip(
            adata.obs[roles["sample"]].tolist(),
            adata.obs[roles["population"]].tolist(),
            strict=True,
        )
    ]
    if profile_order != observed_order:
        raise ValueError(f"{operation} profile identity/order differs from immutable provenance.")
    expected_profile_names = [json.dumps(pair, ensure_ascii=False, separators=(",", ":")) for pair in profile_order]
    if adata.obs_names.tolist() != expected_profile_names:
        raise ValueError(f"{operation} profile row identifiers are not in the canonical Sample/population form.")
    observed_fingerprint = _standalone_artifact_fingerprint(adata, metadata)
    if observed_fingerprint != fingerprint:
        raise ValueError(
            f"{operation} content fingerprint changed: expected {fingerprint}, observed {observed_fingerprint}."
        )
    json.dumps(_plain_json(metadata), ensure_ascii=False, allow_nan=False)
    return adata, metadata


def _standalone_pseudobulk_summary(diagnostics):
    """Build the complete strict pseudobulk aggregation report."""
    import json
    from collections.abc import Mapping

    if not isinstance(diagnostics, Mapping):
        raise TypeError("Pseudobulk summary diagnostics must be a mapping.")
    for field in ("parameters", "warnings", "software_versions"):
        if field not in diagnostics:
            raise ValueError(f"Pseudobulk summary diagnostics are missing {field!r}.")
    key_results = {
        str(key): _plain_json(value)
        for key, value in diagnostics.items()
        if key not in {"parameters", "warnings", "software_versions"}
    }
    parameters = _plain_json(diagnostics["parameters"])
    warnings = [str(value) for value in diagnostics["warnings"]]
    methods = (
        f"Used the analyst-selected {parameters['count_source']['label']!r} axis and validated finite, "
        "non-negative integer-like values without requiring workflow history or equality to Raw; "
        f"Decoupler 2.x pp.pseudobulk summed the declared count source by Sample "
        f"{parameters['sample_key']!r} and population {parameters['population_key']!r}. OpenBio exact integer sums "
        "were authoritative; backend values were checked exactly through 2^53 and for disclosed float64 rounding "
        "equivalence above that limit. pp.filter_samples behavior was verified against its returned QC, synthetic "
        "empty combinations were excluded, and exact profile QC "
        f"required at least {parameters['min_cells']} cells and {parameters['min_counts']} total counts. Sample, "
        "not cell, is the inferential unit for downstream Condition contrasts."
    )
    results = (
        f"Declared counts from {int(key_results['input_cells']):,} cells and "
        f"{int(key_results['input_genes']):,} selected-axis genes "
        f"were summed into {int(key_results['profiles_before_qc']):,} represented Sample-by-population profiles "
        f"across {int(key_results['input_samples']):,} independent Samples; profile QC retained "
        f"{int(key_results['profiles_retained']):,} and removed {int(key_results['profiles_removed']):,}. "
        "This operation created count-model input and did not test a Condition effect."
    )
    references = [
        {
            "citation": "decoupler developers. decoupler.pp.pseudobulk and filter_samples official documentation.",
            "url": "https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.pp.pseudobulk.html",
            "kind": "software_documentation",
            "doi": None,
        },
        {
            "citation": (
                "Badia-I-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological "
                "activities from omics data. Bioinformatics Advances. 2022;2:vbac016."
            ),
            "url": "https://doi.org/10.1093/bioadv/vbac016",
            "kind": "software",
            "doi": "10.1093/bioadv/vbac016",
        },
        {
            "citation": (
                "Squair JW, et al. Confronting false discoveries in single-cell differential expression. "
                "Nature Communications. 2021;12:5692."
            ),
            "url": "https://doi.org/10.1038/s41467-021-25960-2",
            "kind": "practice",
            "doi": "10.1038/s41467-021-25960-2",
        },
        {
            "citation": (
                "Crowell HL, et al. muscat detects subpopulation-specific state transitions from multi-sample "
                "multi-condition single-cell transcriptomics data. Nature Communications. 2020;11:6077."
            ),
            "url": "https://doi.org/10.1038/s41467-020-19894-4",
            "kind": "practice",
            "doi": "10.1038/s41467-020-19894-4",
        },
        {
            "citation": "Virshup I, et al. anndata: Annotated data. Journal of Open Source Software. 2021;6:4371.",
            "url": "https://doi.org/10.21105/joss.04371",
            "kind": "software",
            "doi": "10.21105/joss.04371",
        },
    ]
    limitations = [
        "Source selection, biological raw-count identity, feature completeness, population annotation, and biological Sample identity are analyst declarations; numeric validation and workflow metadata cannot prove them.",
        "Profile min_cells and min_counts thresholds are dataset-dependent QC choices rather than universal validity cutoffs.",
        "Missing Sample-by-population combinations remain absent and may reflect biology, sampling depth, annotation, or QC; they are not zero-count replicates.",
        "Aggregation does not address paired/repeated Samples, technical replicates, interactions, or other designs requiring a dedicated advanced model.",
        "Cell totals describe aggregation depth and must not be used as the replicate count for Condition inference.",
    ]
    summary = {
        "schema_version": 1,
        "node_id": "OpenBioSingleCellPseudobulk",
        "methods": methods,
        "results": results,
        "key_results": key_results,
        "parameters": parameters,
        "warnings": warnings,
        "limitations": limitations,
        "references": references,
        "software_versions": _plain_json(diagnostics["software_versions"]),
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False)
    return summary


def _standalone_run_pseudobulk(
    adata,
    *,
    sample_key,
    population_key,
    condition_key,
    technical_batch_key="",
    categorical_covariate_keys="",
    continuous_covariate_keys="",
    source_kind="raw",
    layer_name=None,
    inference_mode="formal",
    min_cells=10,
    min_counts=1000,
    openbio_version="not-installed",
    decoupler_module=None,
):
    """Validate and aggregate one analyst-selected count source into a typed artifact."""
    import importlib
    import inspect as runtime_inspect
    import json
    import math
    from collections.abc import Mapping
    from importlib import metadata as importlib_metadata

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from scipy import sparse

    operation = "Pseudobulk"

    def required_key(value, description):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{operation} {description} cannot be empty.")
        return value.strip()

    def parse_keys(value, description):
        if not isinstance(value, str):
            raise TypeError(f"{operation} {description} must be a comma-separated string.")
        keys = [item.strip() for item in value.split(",") if item.strip()]
        if len(keys) != len(set(keys)):
            raise ValueError(f"{operation} {description} contains duplicate keys.")
        return keys

    def canonical_labels(series, description):
        labels = []
        identities = {}
        for position, original in enumerate(series.tolist()):
            if not pd.api.types.is_scalar(original):
                raise TypeError(f"{operation} {description} value at position {position} must be scalar.")
            if bool(pd.isna(original)):
                raise ValueError(f"{operation} {description} contains missing values.")
            value = original.item() if isinstance(original, np.generic) else original
            if isinstance(value, str):
                if not value or value != value.strip():
                    raise ValueError(f"{operation} {description} contains blank or padded labels.")
                label = value
            elif isinstance(value, bool):
                raise TypeError(f"{operation} {description} boolean labels are not supported.")
            elif isinstance(value, int):
                label = str(value)
            elif isinstance(value, float):
                if not math.isfinite(value):
                    raise ValueError(f"{operation} {description} contains non-finite numeric labels.")
                label = str(value)
            else:
                raise TypeError(
                    f"{operation} {description} labels must be strings or finite integer/float scalars; "
                    f"observed {type(value).__name__}."
                )
            identity = (type(value).__module__, type(value).__qualname__, repr(value))
            identities.setdefault(label, set()).add(identity)
            labels.append(label)
        collisions = sorted(label for label, values in identities.items() if len(values) > 1)
        if collisions:
            raise ValueError(
                f"{operation} {description} has distinct values that collapse after string conversion: "
                f"{collisions[:8]}."
            )
        return labels

    def matrix_values(matrix):
        return np.asarray(matrix.data) if sparse.issparse(matrix) else np.asarray(matrix).ravel()

    def backend_count_agreement(observed, exact):
        value = observed.item() if isinstance(observed, np.generic) else observed
        if isinstance(value, bool):
            return False, False
        if isinstance(value, int):
            return value == exact, False
        if not isinstance(value, float) or not math.isfinite(value) or value < 0.0 or not value.is_integer():
            return False, False
        if exact <= 2**53:
            return value == exact, False
        return value == float(exact), True

    if not isinstance(adata, AnnData):
        raise TypeError(f"{operation} input must be an AnnData object.")
    if getattr(adata, "isbacked", False):
        raise ValueError(f"{operation} requires an in-memory AnnData; call to_memory() first.")
    if adata.n_obs < 1:
        raise ValueError(f"{operation} requires at least one cell.")
    if not adata.obs_names.is_unique:
        raise ValueError(f"{operation} requires unique cell and gene identifiers.")
    openbio_metadata = adata.uns.get("openbio_singlecell")
    if not isinstance(openbio_metadata, Mapping):
        openbio_metadata = {}

    sample_key = required_key(sample_key, "Sample key")
    population_key = required_key(population_key, "population key")
    condition_key = required_key(condition_key, "Condition key")
    if not isinstance(technical_batch_key, str):
        raise TypeError(f"{operation} Technical batch key must be a string.")
    technical_batch_key = technical_batch_key.strip()
    categorical_keys = parse_keys(categorical_covariate_keys, "categorical covariate keys")
    continuous_keys = parse_keys(continuous_covariate_keys, "continuous covariate keys")
    role_entries = [
        (sample_key, "sample"),
        (population_key, "population"),
        (condition_key, "condition"),
    ]
    if technical_batch_key:
        role_entries.append((technical_batch_key, "technical_batch"))
    role_entries.extend((key, "categorical_covariate") for key in categorical_keys)
    role_entries.extend((key, "continuous_covariate") for key in continuous_keys)
    roles_by_key = {}
    for key, role in role_entries:
        roles_by_key.setdefault(key, []).append(role)
    role_aliases = {key: roles for key, roles in roles_by_key.items() if len(roles) > 1}
    role_keys = [key for key, _ in role_entries]
    reserved = {"openbio_n_cells", "openbio_total_counts"}
    if reserved & set(role_keys):
        raise ValueError(f"{operation} metadata roles collide with reserved profile-QC columns.")
    missing_columns = [key for key in role_keys if key not in adata.obs]
    if missing_columns:
        raise ValueError(f"{operation} observation columns were not found: {missing_columns}.")

    if source_kind not in {"raw", "X", "layer"}:
        raise ValueError(f"{operation} unsupported count source {source_kind!r}.")
    if source_kind == "layer":
        layer_name = required_key(layer_name, "count layer")
        if layer_name not in adata.layers:
            raise ValueError(f"{operation} count layer not found: {layer_name!r}.")
        source_matrix = adata.layers[layer_name]
        source_var = adata.var
        source_label = f"layers[{layer_name!r}]"
        feature_axis = "adata.var_names"
    elif source_kind == "X":
        if layer_name is not None:
            raise ValueError(f"{operation} X source cannot carry a layer name.")
        source_matrix = adata.X
        source_var = adata.var
        source_label = "X"
        feature_axis = "adata.var_names"
    else:
        if layer_name is not None:
            raise ValueError(f"{operation} Raw source cannot carry a layer name.")
        if adata.raw is None:
            raise ValueError(f"{operation} Raw source was selected but adata.raw is absent.")
        if not adata.raw.obs_names.equals(adata.obs_names):
            raise ValueError(f"{operation} Raw observations are not aligned to current cells.")
        source_matrix = adata.raw.X
        source_var = adata.raw.var
        source_label = "raw.X"
        feature_axis = "raw.var_names"
    if len(source_var.index) < 1 or not source_var.index.is_unique:
        raise ValueError(f"{operation} requires unique cell and gene identifiers on the selected source axis.")
    if getattr(source_matrix, "shape", None) != (adata.n_obs, len(source_var.index)):
        raise ValueError(f"{operation} selected count source is not aligned to its observation/feature axes.")
    values = matrix_values(source_matrix)
    if np.issubdtype(values.dtype, np.bool_) or not np.issubdtype(values.dtype, np.number):
        raise TypeError(f"{operation} raw counts must be numeric and non-boolean.")
    numeric = values.astype(float, copy=False)
    if numeric.size and not bool(np.isfinite(numeric).all()):
        raise ValueError(f"{operation} raw counts contain non-finite values.")
    if numeric.size and bool((numeric < 0).any()):
        raise ValueError(f"{operation} raw counts contain negative values.")
    if numeric.size and not bool(np.allclose(numeric, np.rint(numeric), rtol=0.0, atol=1e-8)):
        raise ValueError(f"{operation} requires integer-like raw counts; transformed values are unsupported.")
    input_total = float(numeric.sum(dtype=np.float64))
    if not math.isfinite(input_total) or input_total <= 0:
        raise ValueError(f"{operation} raw count matrix contains no positive counts.")
    if values.size:
        max_value = int(values.max()) if np.issubdtype(values.dtype, np.integer) else float(numeric.max())
        if max_value > int(np.iinfo(np.int64).max):
            raise OverflowError(f"{operation} raw count exceeds int64 range.")
    sample_labels = canonical_labels(adata.obs[sample_key], "Sample labels")
    population_labels = canonical_labels(adata.obs[population_key], "population labels")
    condition_labels = canonical_labels(adata.obs[condition_key], "Condition labels")
    categorical_values = {
        key: canonical_labels(adata.obs[key], f"categorical covariate {key!r}") for key in categorical_keys
    }
    technical_values = (
        canonical_labels(adata.obs[technical_batch_key], "Technical batch labels") if technical_batch_key else None
    )
    continuous_values = {}
    for key in continuous_keys:
        series = adata.obs[key]
        if series.isna().any() or pd.api.types.is_bool_dtype(series.dtype):
            raise ValueError(f"{operation} continuous covariate {key!r} must be complete numeric values.")
        try:
            array = pd.to_numeric(series, errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{operation} continuous covariate {key!r} must be numeric.") from exc
        if not bool(np.isfinite(array).all()):
            raise ValueError(f"{operation} continuous covariate {key!r} must be finite.")
        continuous_values[key] = array.tolist()

    sample_order = list(dict.fromkeys(sample_labels))
    population_order = list(dict.fromkeys(population_labels))
    condition_order = list(dict.fromkeys(condition_labels))
    sample_metadata = {}
    conflicts = {}
    for index, sample in enumerate(sample_labels):
        row = {
            condition_key: condition_labels[index],
            **({technical_batch_key: technical_values[index]} if technical_values is not None else {}),
            **{key: categorical_values[key][index] for key in categorical_keys},
            **{key: continuous_values[key][index] for key in continuous_keys},
        }
        if sample not in sample_metadata:
            sample_metadata[sample] = row
            continue
        differing = [key for key, value in row.items() if sample_metadata[sample][key] != value]
        if differing:
            conflicts.setdefault(sample, set()).update(differing)
    if conflicts:
        rendered = {sample: sorted(keys) for sample, keys in sorted(conflicts.items())}
        raise ValueError(
            f"{operation} requires one Condition/batch/covariate mapping per Sample; conflicts={rendered}."
        )

    annotations = openbio_metadata.get("annotations") if isinstance(openbio_metadata, Mapping) else None
    annotation_provenance = annotations.get(population_key) if isinstance(annotations, Mapping) else None
    annotation_status = "unknown"
    annotation_operation = None
    if isinstance(annotation_provenance, Mapping) and annotation_provenance.get("schema_version") == 1:
        candidate_operation = annotation_provenance.get("operation")
        provenance_matches_column = False
        if candidate_operation == "map_cluster_annotations":
            provenance_matches_column = annotation_provenance.get("output_column") == population_key
        elif candidate_operation == "celltypist_annotation":
            columns = annotation_provenance.get("columns")
            provenance_matches_column = isinstance(columns, Mapping) and columns.get("selected_label") == population_key
        candidate_status = annotation_provenance.get("annotation_status")
        if (
            provenance_matches_column
            and isinstance(candidate_status, str)
            and candidate_status.lower() in {"curated", "provisional"}
        ):
            annotation_status = candidate_status.lower()
            annotation_operation = candidate_operation
    if inference_mode not in {"formal", "exploratory"}:
        raise ValueError(f"{operation} inference_mode must be 'formal' or 'exploratory'.")
    formal_interpretation_invalid_reasons = []
    if inference_mode == "exploratory":
        formal_interpretation_invalid_reasons.append("exploratory_mode_requested")
    elif annotation_status != "curated":
        formal_interpretation_invalid_reasons.append("population_annotation_not_curated")
    if role_aliases:
        formal_interpretation_invalid_reasons.append("metadata_role_aliases")
    formal_interpretation_invalid = bool(formal_interpretation_invalid_reasons)
    inference_status = (
        "exploratory_sample_level_condition_evidence"
        if inference_mode == "exploratory"
        else "formal_interpretation_invalid"
        if formal_interpretation_invalid
        else "formal_sample_level_condition_inference_ready"
    )
    declarations = {
        "count_source": "caller_selected_and_declared_raw_counts",
        "feature_completeness": "caller_declared_selected_axis",
        "sample_identity": "caller_declared_independent_biological_sample",
        "population_annotation": "upstream_metadata_or_unknown",
    }
    for name, value, minimum in (("min_cells", min_cells, 1), ("min_counts", min_counts, 0)):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{operation} {name} must be an integer >= {minimum}.")
    if not isinstance(openbio_version, str) or not openbio_version.strip():
        raise ValueError(f"{operation} openbio_version cannot be empty.")

    pair_order = list(dict.fromkeys(zip(sample_labels, population_labels, strict=True)))
    positions_by_pair = {}
    for index, pair in enumerate(zip(sample_labels, population_labels, strict=True)):
        positions_by_pair.setdefault(pair, []).append(index)
    exact_profile_totals = {pair: 0 for pair in pair_order}
    if sparse.issparse(source_matrix):
        pairs_by_row = list(zip(sample_labels, population_labels, strict=True))
        coordinate_matrix = source_matrix.tocoo(copy=False)
        for row_index, value in zip(coordinate_matrix.row, coordinate_matrix.data, strict=True):
            pair = pairs_by_row[int(row_index)]
            exact_profile_totals[pair] += int(value)
    else:
        source_values = np.asarray(source_matrix)
        for pair, positions in positions_by_pair.items():
            exact_profile_totals[pair] = sum(int(value) for value in source_values[positions].ravel())
    int64_max = int(np.iinfo(np.int64).max)
    for (sample, population), exact_total in exact_profile_totals.items():
        if exact_total > int64_max:
            raise OverflowError(
                f"Pseudobulk aggregate library for Sample={sample!r}, population={population!r} "
                "exceeds the supported int64 range."
            )
    counts = (
        source_matrix.astype(np.int64).tocsr()
        if sparse.issparse(source_matrix)
        else np.asarray(source_matrix, dtype=np.int64)
    )
    obs_data = {
        sample_key: sample_labels,
        population_key: population_labels,
        condition_key: condition_labels,
    }
    if technical_batch_key:
        obs_data[technical_batch_key] = technical_values
    obs_data.update(categorical_values)
    for key, values_for_key in continuous_values.items():
        obs_data.setdefault(key, values_for_key)
    work = AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(obs_data, index=adata.obs_names.copy()),
        var=source_var.copy(),
    )
    if decoupler_module is None:
        try:
            decoupler_module = importlib.import_module("decoupler")
        except (ImportError, OSError) as exc:
            raise RuntimeError(
                "Pseudobulk requires Decoupler 2.x; install a compatible decoupler release before execution."
            ) from exc
    decoupler_version_value = getattr(decoupler_module, "__version__", None)
    if decoupler_version_value is None:
        try:
            decoupler_version_value = importlib_metadata.version("decoupler")
        except importlib_metadata.PackageNotFoundError:
            decoupler_version_value = "unknown"
    decoupler_version = str(decoupler_version_value)
    try:
        major_version = int(decoupler_version.split(".", 1)[0])
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"Pseudobulk could not parse Decoupler version {decoupler_version!r}.") from exc
    if major_version != 2:
        raise RuntimeError(f"Pseudobulk requires Decoupler 2.x; detected {decoupler_version!r}.")
    pseudobulk_function = getattr(getattr(decoupler_module, "pp", None), "pseudobulk", None)
    if not callable(pseudobulk_function):
        raise RuntimeError("Pseudobulk requires the Decoupler 2.x pp.pseudobulk public interface.")
    filter_samples_function = getattr(getattr(decoupler_module, "pp", None), "filter_samples", None)
    if not callable(filter_samples_function):
        raise RuntimeError("Pseudobulk requires the Decoupler 2.x pp.filter_samples public interface.")
    required_parameters = {
        "adata",
        "sample_col",
        "groups_col",
        "layer",
        "raw",
        "empty",
        "mode",
        "skip_checks",
        "bsize",
        "verbose",
    }
    try:
        signature_parameters = set(runtime_inspect.signature(pseudobulk_function).parameters)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("Pseudobulk could not inspect Decoupler pp.pseudobulk.") from exc
    if not required_parameters.issubset(signature_parameters):
        missing = sorted(required_parameters - signature_parameters)
        raise RuntimeError(
            f"Pseudobulk Decoupler {decoupler_version} pp.pseudobulk interface is incompatible; missing={missing}."
        )
    try:
        filter_parameters = set(runtime_inspect.signature(filter_samples_function).parameters)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("Pseudobulk could not inspect Decoupler pp.filter_samples.") from exc
    required_filter_parameters = {"adata", "min_cells", "min_counts", "inplace"}
    if not required_filter_parameters.issubset(filter_parameters):
        missing = sorted(required_filter_parameters - filter_parameters)
        raise RuntimeError(
            f"Pseudobulk Decoupler {decoupler_version} pp.filter_samples interface is incompatible; missing={missing}."
        )
    try:
        backend = pseudobulk_function(
            work,
            sample_col=sample_key,
            groups_col=population_key,
            layer=None,
            raw=False,
            empty=False,
            mode="sum",
            skip_checks=False,
            bsize=250_000,
            verbose=False,
        )
    except Exception as exc:
        raise RuntimeError(f"Pseudobulk Decoupler aggregation failed ({type(exc).__name__}: {exc}).") from exc
    if not isinstance(backend, AnnData):
        raise RuntimeError("Pseudobulk Decoupler backend returned a non-AnnData result.")
    if not backend.var_names.equals(source_var.index):
        raise RuntimeError("Pseudobulk Decoupler backend changed the selected feature identity or order.")
    if sample_key not in backend.obs or population_key not in backend.obs:
        raise RuntimeError("Pseudobulk Decoupler backend omitted Sample/population identity columns.")

    expected_rows = []
    profile_records = []
    for sample, population in pair_order:
        positions = positions_by_pair[(sample, population)]
        subset = counts[positions]
        exact_total = exact_profile_totals[(sample, population)]
        summed = np.asarray(subset.sum(axis=0, dtype=np.int64)).ravel().astype(np.int64, copy=False)
        total = int(sum(int(value) for value in summed.tolist()))
        if total != exact_total:
            raise RuntimeError(
                f"Pseudobulk aggregate for Sample={sample!r}, population={population!r} "
                "did not preserve the exact non-negative profile total."
            )
        if total <= 0:
            raise ValueError(
                f"Pseudobulk represented profile Sample={sample!r}, population={population!r} has a zero library."
            )
        expected_rows.append(summed)
        reasons = []
        if len(positions) < min_cells:
            reasons.append("below_min_cells")
        if total < min_counts:
            reasons.append("below_min_counts")
        profile_records.append(
            {
                "sample": sample,
                "population": population,
                "condition": sample_metadata[sample][condition_key],
                "n_cells": len(positions),
                "total_counts": total,
                "retained": not reasons,
                "removal_reasons": reasons,
            }
        )
    expected_matrix = np.vstack(expected_rows).astype(np.int64, copy=False)

    backend_samples = canonical_labels(backend.obs[sample_key], "backend Sample labels")
    backend_populations = canonical_labels(backend.obs[population_key], "backend population labels")
    backend_pairs = list(zip(backend_samples, backend_populations, strict=True))
    if len(backend_pairs) != len(set(backend_pairs)):
        raise RuntimeError("Pseudobulk Decoupler backend returned duplicate Sample-by-population profiles.")
    backend_by_pair = {pair: index for index, pair in enumerate(backend_pairs)}
    cell_column = None
    for candidate in ("psbulk_cells", "psbulk_n_cells"):
        if candidate in backend.obs:
            cell_column = candidate
            break
    if cell_column is None or "psbulk_counts" not in backend.obs:
        raise RuntimeError("Pseudobulk Decoupler backend omitted its cell/count QC fields.")
    backend_cell_series = pd.to_numeric(backend.obs[cell_column], errors="raise")
    backend_count_series = pd.to_numeric(backend.obs["psbulk_counts"], errors="raise")
    backend_cells = backend_cell_series.to_numpy(dtype=float)
    backend_counts = backend_count_series.to_numpy(dtype=float)
    backend_cell_scalars = backend_cell_series.tolist()
    backend_count_scalars = backend_count_series.tolist()
    if "psbulk_cells" in backend.obs and "psbulk_n_cells" in backend.obs:
        documented_cells = pd.to_numeric(backend.obs["psbulk_n_cells"], errors="raise").to_numpy(dtype=float)
        if not bool(np.array_equal(backend_cells, documented_cells)):
            raise RuntimeError("Pseudobulk Decoupler backend cell-QC aliases disagree.")
    # The Decoupler result is worker-owned. Reuse it for its non-mutating filter API instead of
    # retaining a second full AnnData merely to protect this private intermediate.
    filter_input = backend
    if "psbulk_cells" not in filter_input.obs and "psbulk_n_cells" in filter_input.obs:
        filter_input.obs["psbulk_cells"] = filter_input.obs["psbulk_n_cells"].to_numpy(copy=True)
    try:
        backend_retained_names = filter_samples_function(
            filter_input,
            min_cells=min_cells,
            min_counts=min_counts,
            inplace=False,
        )
    except Exception as exc:
        raise RuntimeError(f"Pseudobulk Decoupler profile filtering failed ({type(exc).__name__}: {exc}).") from exc
    if backend_retained_names is None:
        raise RuntimeError("Pseudobulk Decoupler pp.filter_samples(inplace=False) returned no profile identities.")
    backend_retained_names = [str(value) for value in np.asarray(backend_retained_names).ravel().tolist()]
    independently_retained_names = (
        backend.obs_names[(backend_cells >= min_cells) & (backend_counts >= min_counts)].astype(str).tolist()
    )
    if backend_retained_names != independently_retained_names:
        raise RuntimeError("Pseudobulk Decoupler profile-filter identities disagree with independent QC checks.")
    backend_matrix = backend.X.tocsr() if sparse.issparse(backend.X) else np.asarray(backend.X)
    backend_precision_limited_profiles = []
    for expected_index, pair in enumerate(pair_order):
        if pair not in backend_by_pair:
            raise RuntimeError(f"Pseudobulk Decoupler backend omitted represented profile {pair!r}.")
        backend_index = backend_by_pair[pair]
        observed = (
            backend_matrix.getrow(backend_index).toarray().ravel()
            if sparse.issparse(backend_matrix)
            else np.asarray(backend_matrix[backend_index]).ravel()
        )
        count_agreements = [
            backend_count_agreement(observed_value, int(exact_value))
            for observed_value, exact_value in zip(observed, expected_matrix[expected_index], strict=True)
        ]
        if not all(matches for matches, _ in count_agreements):
            raise RuntimeError(f"Pseudobulk Decoupler backend counts disagree for profile {pair!r}.")
        record = profile_records[expected_index]
        cells_match, cells_precision_limited = backend_count_agreement(
            backend_cell_scalars[backend_index],
            record["n_cells"],
        )
        total_match, total_precision_limited = backend_count_agreement(
            backend_count_scalars[backend_index],
            record["total_counts"],
        )
        if not cells_match or not total_match:
            raise RuntimeError(f"Pseudobulk Decoupler backend QC counts disagree for profile {pair!r}.")
        if cells_precision_limited or total_precision_limited or any(
            precision_limited for _, precision_limited in count_agreements
        ):
            backend_precision_limited_profiles.append(
                {
                    "sample": record["sample"],
                    "population": record["population"],
                    "total_counts": record["total_counts"],
                }
            )
    extra_pairs = [pair for pair in backend_pairs if pair not in positions_by_pair]
    for pair in extra_pairs:
        backend_index = backend_by_pair[pair]
        observed = (
            backend_matrix.getrow(backend_index).toarray().ravel()
            if sparse.issparse(backend_matrix)
            else np.asarray(backend_matrix[backend_index]).ravel()
        )
        if backend_cells[backend_index] != 0 or backend_counts[backend_index] != 0 or bool(np.any(observed != 0)):
            raise RuntimeError(
                f"Pseudobulk Decoupler backend returned an unexpected nonempty synthetic profile {pair!r}."
            )

    retained_indices = [index for index, record in enumerate(profile_records) if record["retained"]]
    if not retained_indices:
        raise ValueError("Pseudobulk profile QC removed every represented Sample-by-population profile.")
    before_strata = {}
    after_strata = {}
    for record in profile_records:
        key = (record["population"], record["condition"])
        before_strata[key] = before_strata.get(key, 0) + 1
        if record["retained"]:
            after_strata[key] = after_strata.get(key, 0) + 1
    lost_strata = [list(key) for key in before_strata if after_strata.get(key, 0) == 0]
    output_rows = expected_matrix[retained_indices]
    output_records = [profile_records[index] for index in retained_indices]
    output_obs_data = {
        sample_key: [record["sample"] for record in output_records],
        population_key: [record["population"] for record in output_records],
        condition_key: [record["condition"] for record in output_records],
    }
    if technical_batch_key:
        output_obs_data[technical_batch_key] = [
            sample_metadata[record["sample"]][technical_batch_key] for record in output_records
        ]
    for key in categorical_keys:
        output_obs_data[key] = [sample_metadata[record["sample"]][key] for record in output_records]
    for key in continuous_keys:
        output_obs_data.setdefault(
            key,
            [sample_metadata[record["sample"]][key] for record in output_records],
        )
    output_obs_data["openbio_n_cells"] = [record["n_cells"] for record in output_records]
    output_obs_data["openbio_total_counts"] = [record["total_counts"] for record in output_records]
    profile_names = [
        json.dumps([record["sample"], record["population"]], ensure_ascii=False, separators=(",", ":"))
        for record in output_records
    ]
    output = AnnData(
        X=sparse.csr_matrix(output_rows) if sparse.issparse(source_matrix) else output_rows,
        obs=pd.DataFrame(output_obs_data, index=profile_names),
        var=source_var.copy(),
    )

    all_pairs = [(sample, population) for sample in sample_order for population in population_order]
    missing_combinations = [
        {
            "sample": sample,
            "population": population,
            "condition": sample_metadata[sample][condition_key],
        }
        for sample, population in all_pairs
        if (sample, population) not in positions_by_pair
    ]
    stratum_records = []
    for population, condition in sorted(before_strata):
        matching = [
            record
            for record in profile_records
            if record["population"] == population and record["condition"] == condition
        ]
        retained_matching = [record for record in matching if record["retained"]]
        stratum_records.append(
            {
                "population": population,
                "condition": condition,
                "profiles_before_qc": len(matching),
                "profiles_after_qc": len(retained_matching),
                "contributing_cells_before_qc": sum(record["n_cells"] for record in matching),
                "contributing_cells_after_qc": sum(record["n_cells"] for record in retained_matching),
                "total_counts_before_qc": sum(record["total_counts"] for record in matching),
                "total_counts_after_qc": sum(record["total_counts"] for record in retained_matching),
            }
        )

    def numeric_summary(values):
        array = np.asarray(list(values), dtype=float)
        quantiles = np.quantile(array, [0.0, 0.25, 0.5, 0.75, 1.0])
        return {
            "n": int(array.size),
            "min": float(quantiles[0]),
            "q1": float(quantiles[1]),
            "median": float(quantiles[2]),
            "mean": float(array.mean()),
            "q3": float(quantiles[3]),
            "max": float(quantiles[4]),
        }

    parameters = {
        "sample_key": sample_key,
        "population_key": population_key,
        "condition_key": condition_key,
        "technical_batch_key": technical_batch_key or None,
        "categorical_covariate_keys": categorical_keys,
        "continuous_covariate_keys": continuous_keys,
        "count_source": {
            "kind": source_kind,
            "layer_name": layer_name if source_kind == "layer" else None,
            "label": source_label,
            "feature_axis": feature_axis,
        },
        "declarations": declarations,
        "role_aliases": role_aliases,
        "annotation_status": annotation_status,
        "annotation_provenance_operation": annotation_operation,
        "inference_mode": inference_mode,
        "inference_status": inference_status,
        "formal_interpretation_invalid": formal_interpretation_invalid,
        "formal_interpretation_invalid_reasons": formal_interpretation_invalid_reasons,
        "min_cells": min_cells,
        "min_counts": min_counts,
        "mode": "sum",
        "synthetic_empty_profiles": False,
        "summation_policy": {
            "output_dtype": "int64",
            "overflow_check": "exact_nonnegative_profile_total_before_backend",
            "backend_exact_integer_limit": 2**53,
            "backend_above_limit_check": "float64_rounding_equivalence_with_disclosure",
        },
        "backend_policy": {
            "raw": False,
            "empty": False,
            "skip_checks": False,
            "bsize": 250_000,
            "verbose": False,
        },
    }
    artifact_metadata = {
        "schema_version": 2,
        "artifact_type": "OPENBIO_SINGLE_CELL_PSEUDOBULK",
        "producer_node_id": "OpenBioSingleCellPseudobulk",
        "aggregation_mode": "sum",
        "role_keys": {
            "sample": sample_key,
            "population": population_key,
            "condition": condition_key,
            "technical_batch": technical_batch_key or None,
        },
        "declared_covariates": {
            "categorical": categorical_keys,
            "continuous": continuous_keys,
        },
        "annotation_status": annotation_status,
        "annotation_provenance_operation": annotation_operation,
        "inference_status": parameters["inference_status"],
        "declarations": declarations,
        "role_aliases": role_aliases,
        "formal_interpretation_invalid": formal_interpretation_invalid,
        "formal_interpretation_invalid_reasons": formal_interpretation_invalid_reasons,
        "count_source": parameters["count_source"],
        "profile_filters": {"min_cells": min_cells, "min_counts": min_counts},
        "input_dimensions": {
            "cells": int(adata.n_obs),
            "genes": int(len(source_var.index)),
            "samples": len(sample_order),
        },
        "profile_records": profile_records,
        "missing_combinations": missing_combinations,
        "profile_order": [[record["sample"], record["population"]] for record in output_records],
        "obs_columns": list(output.obs.columns),
        "backend": {
            "package": "decoupler",
            "version": decoupler_version,
            "interface": "decoupler.pp.pseudobulk",
            "cell_qc_column": cell_column,
            "synthetic_empty_profiles_removed": len(extra_pairs),
        },
    }
    artifact_fingerprint = _standalone_artifact_fingerprint(output, artifact_metadata)
    artifact_metadata["artifact_fingerprint"] = artifact_fingerprint
    artifact = _StandalonePseudobulkArtifact._from_owned(output, artifact_metadata)
    _standalone_validate_pseudobulk_artifact(artifact, copy_result=False)

    warnings = [
        "Profile min_cells and min_counts are dataset-dependent QC choices; sensitivity analysis may be warranted.",
        f"Count source {source_label!r}, raw-count identity, and feature completeness are caller declarations; software validated numeric count invariants but did not require workflow history or Raw equivalence.",
        f"Column {sample_key!r} is caller-declared as independent biological Sample identity; software checked within-Sample metadata consistency but cannot prove biological independence.",
    ]
    if backend_precision_limited_profiles:
        warnings.append(
            f"Decoupler returned float-backed aggregate/QC values for {len(backend_precision_limited_profiles):,} "
            "profile(s) containing integers above 2^53. Backend agreement is rounding-equivalent rather than exact; "
            "the emitted int64 counts and profile-QC decisions use OpenBio's authoritative exact sums."
        )
    if inference_mode == "exploratory":
        warnings.append(
            f"Exploratory mode was requested for population annotation status {annotation_status!r}; formal interpretation is invalid."
        )
    elif annotation_status != "curated":
        warnings.append(
            f"Population annotation status is {annotation_status!r}; aggregation remains executable, but formal interpretation is invalid until the analyst justifies the population labels."
        )
    if role_aliases:
        warnings.append(
            f"Metadata keys were declared for multiple roles; aggregation remains executable, but formal interpretation is invalid unless the aliases are intentional: {role_aliases}."
        )
    if missing_combinations:
        warnings.append(
            f"{len(missing_combinations):,} Sample-by-population combinations had no cells and remain absent rather than zero-filled."
        )
    removed_count = len(profile_records) - len(output_records)
    if removed_count:
        warnings.append(f"Profile QC removed {removed_count:,} represented Sample-by-population profiles.")
    if lost_strata:
        warnings.append(
            f"Profile QC removed an entire represented population/Condition stratum; the artifact remains usable for other estimable selections: {lost_strata[:20]}."
        )
    fragile = [
        {"population": population, "condition": condition, "samples": count}
        for (population, condition), count in sorted(after_strata.items())
        if count < 3
    ]
    if fragile:
        warnings.append(
            "Some retained population/Condition strata contain fewer than three declared independent Samples; "
            f"downstream inference may be fragile or non-estimable depending on the complete design: {fragile[:20]}."
        )
    if not technical_batch_key:
        warnings.append(
            "No Technical batch key was declared; downstream models cannot adjust an undeclared batch effect."
        )
    software_versions = collect_software_versions(
        ("decoupler", "anndata", "numpy", "pandas", "scipy"),
        openbio_version=str(openbio_version),
    )
    software_versions["decoupler"] = str(decoupler_version)
    diagnostics = {
        "artifact_fingerprint_sha256": artifact_fingerprint,
        "inference_status": parameters["inference_status"],
        "declarations": declarations,
        "role_aliases": role_aliases,
        "formal_interpretation_invalid": formal_interpretation_invalid,
        "formal_interpretation_invalid_reasons": formal_interpretation_invalid_reasons,
        "annotation_status": annotation_status,
        "input_cells": int(adata.n_obs),
        "input_genes": int(len(source_var.index)),
        "input_samples": len(sample_order),
        "population_count": len(population_order),
        "condition_count": len(condition_order),
        "profiles_before_qc": len(profile_records),
        "profiles_retained": len(output_records),
        "profiles_removed": removed_count,
        "lost_strata": lost_strata,
        "contributing_cells_retained": sum(record["n_cells"] for record in output_records),
        "counts_retained": sum(record["total_counts"] for record in output_records),
        "missing_combination_count": len(missing_combinations),
        "missing_combinations_preview": missing_combinations[:100],
        "missing_combinations_preview_truncated": len(missing_combinations) > 100,
        "backend_synthetic_empty_profiles_removed": len(extra_pairs),
        "backend_count_verification": {
            "authoritative_sum": "openbio_exact_python_integer_then_int64",
            "backend_matrix_dtype": str(backend_matrix.dtype),
            "backend_qc_dtype": str(backend_count_series.dtype),
            "exact_integer_limit": 2**53,
            "precision_limited_profile_count": len(backend_precision_limited_profiles),
            "precision_limited_profiles_preview": backend_precision_limited_profiles[:100],
            "precision_limited_profiles_preview_truncated": len(backend_precision_limited_profiles) > 100,
        },
        "retained_profile_records_preview": output_records[:100],
        "retained_profile_records_preview_truncated": len(output_records) > 100,
        "removed_profile_records_preview": [record for record in profile_records if not record["retained"]][:100],
        "removed_profile_records_preview_truncated": removed_count > 100,
        "strata_preview": stratum_records[:100],
        "strata_preview_truncated": len(stratum_records) > 100,
        "retained_cell_count_distribution": numeric_summary(record["n_cells"] for record in output_records),
        "retained_library_size_distribution": numeric_summary(record["total_counts"] for record in output_records),
        "backend": artifact_metadata["backend"],
        "parameters": parameters,
        "warnings": warnings,
        "software_versions": software_versions,
    }
    json.dumps(_plain_json(diagnostics), ensure_ascii=False, allow_nan=False)
    return artifact, diagnostics


def run_pseudobulk(
    adata: AnnData,
    **kwargs: Any,
) -> tuple[PseudobulkArtifact, dict[str, Any]]:
    return _standalone_run_pseudobulk(adata, **kwargs)


def build_pseudobulk_summary(diagnostics: dict[str, Any]) -> dict[str, Any]:
    return _standalone_pseudobulk_summary(diagnostics)


def validate_pseudobulk_artifact(artifact: Any) -> tuple[AnnData, dict[str, Any]]:
    return _standalone_validate_pseudobulk_artifact(artifact)


def pseudobulk_code(**parameters: Any) -> str:
    implementations = "\n\n".join(
        textwrap.dedent(inspect.getsource(function)).strip()
        for function in (
            _plain_json,
            _package_version,
            collect_software_versions,
            _standalone_artifact_fingerprint,
            _StandalonePseudobulkArtifact,
            _standalone_validate_pseudobulk_artifact,
            _standalone_pseudobulk_summary,
            _standalone_run_pseudobulk,
        )
    )
    rendered = ",\n        ".join(f"{key}={value!r}" for key, value in parameters.items())
    return f"""from __future__ import annotations

{implementations}


def run_pseudobulk(adata):
    artifact, diagnostics = _standalone_run_pseudobulk(
        adata,
        {rendered}
    )
    return artifact, _standalone_pseudobulk_summary(diagnostics)
"""


__all__ = [
    "PSEUDOBULK_ARTIFACT_TYPE",
    "PSEUDOBULK_QC_COLUMNS",
    "PSEUDOBULK_SCHEMA_VERSION",
    "PseudobulkArtifact",
    "build_pseudobulk_summary",
    "pseudobulk_code",
    "run_pseudobulk",
    "validate_pseudobulk_artifact",
]
