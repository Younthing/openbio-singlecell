from __future__ import annotations

import inspect
import textwrap
import threading
from typing import Any

from .milo_result import MILO_RESULT_COLUMNS

MILO_TABLE_COLUMNS = list(MILO_RESULT_COLUMNS)

_MILO_RNG_LOCK = threading.RLock()


def _standalone_milo_differential_abundance_impl(
    adata,
    *,
    sample_key,
    condition_key,
    reference_condition,
    comparison_condition,
    technical_batch_key="",
    categorical_covariate_keys_json="[]",
    continuous_covariate_keys_json="[]",
    annotation_key="cell_type",
    annotation_status="provisional",
    representation_key="X_pca",
    n_neighbors=30,
    neighborhood_proportion=0.1,
    mixed_annotation_threshold=0.6,
    spatial_fdr_threshold=0.1,
    min_abs_log2_fold_change=0.0,
    random_seed=123,
    openbio_version="unknown",
):
    """Run one strict pairwise Sample-level Milo contrast without OpenBio runtime helpers."""
    import hashlib
    import importlib
    import inspect as runtime_inspect
    import json
    import math
    import numbers
    import platform
    import random
    from importlib import metadata as importlib_metadata

    import numpy as np
    import pandas as pd
    from scipy import sparse

    operation = "Milo Differential Abundance"
    feature_key = "rna"

    def canonical_text(value, label, *, allow_empty=False):
        if not isinstance(value, str):
            raise TypeError(f"{operation} {label} must be a string.")
        if value != value.strip():
            raise ValueError(f"{operation} {label} cannot contain surrounding whitespace.")
        if not value and not allow_empty:
            raise ValueError(f"{operation} {label} cannot be empty.")
        return value

    def numeric_value(value, label, *, integer=False, minimum=None, maximum=None, maximum_inclusive=True):
        expected = numbers.Integral if integer else numbers.Real
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, expected):
            kind = "integer" if integer else "number"
            article = "an" if integer else "a"
            raise TypeError(f"{operation} {label} must be {article} {kind}.")
        converted = int(value) if integer else float(value)
        if not integer and not math.isfinite(converted):
            raise ValueError(f"{operation} {label} must be finite.")
        if minimum is not None and converted < minimum:
            raise ValueError(f"{operation} {label} must be at least {minimum}.")
        if maximum is not None:
            invalid = converted > maximum if maximum_inclusive else converted >= maximum
            if invalid:
                relation = "at most" if maximum_inclusive else "less than"
                raise ValueError(f"{operation} {label} must be {relation} {maximum}.")
        return converted

    def parse_key_list(value, label):
        if not isinstance(value, str):
            raise TypeError(f"{operation} {label} must be a JSON string array.")
        if len(value.encode("utf-8")) > 65_536:
            raise ValueError(f"{operation} {label} exceeds the 65,536-byte safety limit.")

        def reject_constant(constant):
            raise ValueError(f"{operation} {label} contains non-standard JSON constant {constant!r}.")

        try:
            parsed = json.loads(value, parse_constant=reject_constant)
        except json.JSONDecodeError as error:
            raise ValueError(f"{operation} {label} is not valid JSON ({error.msg}).") from error
        if not isinstance(parsed, list):
            raise TypeError(f"{operation} {label} must decode to a JSON array.")
        if len(parsed) > 32:
            raise ValueError(f"{operation} {label} supports at most 32 covariates.")
        result = []
        for position, item in enumerate(parsed):
            item = canonical_text(item, f"{label}[{position}]")
            if item in result:
                raise ValueError(f"{operation} {label} contains duplicate column {item!r}.")
            result.append(item)
        return result

    def string_values(series, label):
        values = []
        for position, value in enumerate(series.tolist()):
            missing = pd.isna(value)
            if not isinstance(missing, (bool, np.bool_)):
                raise TypeError(f"{operation} {label} row {position} is nonscalar.")
            if bool(missing):
                raise ValueError(f"{operation} {label} contains missing values.")
            if not isinstance(value, str):
                raise TypeError(
                    f"{operation} {label} must contain strings; row {position} is {type(value).__name__}."
                )
            if not value or value != value.strip():
                raise ValueError(f"{operation} {label} contains an empty or non-canonical value at row {position}.")
            values.append(value)
        return values

    def observed_levels(series, values):
        observed = set(values)
        if isinstance(series.dtype, pd.CategoricalDtype):
            declared = []
            for category in series.cat.categories.tolist():
                category = canonical_text(category, "categorical level")
                if category in observed:
                    declared.append(category)
            if len(declared) != len(observed):
                raise RuntimeError(f"{operation} failed to preserve categorical level identities.")
            return declared
        return list(dict.fromkeys(values))

    def safe_column(base, occupied):
        candidate = base
        suffix = 1
        while candidate in occupied:
            candidate = f"{base}_{suffix}"
            suffix += 1
        occupied.add(candidate)
        return candidate

    def matrix_payload(matrix, label, *, rows=None, columns=None, binary=False, nonnegative=False):
        if sparse.issparse(matrix):
            canonical = sparse.csr_matrix(matrix, copy=True)
            canonical.sum_duplicates()
            canonical.eliminate_zeros()
            canonical.sort_indices()
            values = np.asarray(canonical.data)
        else:
            array = np.asarray(matrix)
            if array.ndim != 2:
                raise ValueError(f"{operation} {label} must be two-dimensional.")
            canonical = array
            values = array.ravel()
        if len(canonical.shape) != 2:
            raise ValueError(f"{operation} {label} must be two-dimensional.")
        if rows is not None and int(canonical.shape[0]) != int(rows):
            raise RuntimeError(f"{operation} {label} has the wrong row axis.")
        if columns is not None and int(canonical.shape[1]) != int(columns):
            raise RuntimeError(f"{operation} {label} has the wrong column axis.")
        if values.size and not bool(np.issubdtype(values.dtype, np.number)):
            raise TypeError(f"{operation} {label} must be numeric.")
        if values.size and not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} {label} contains non-finite values.")
        if nonnegative and values.size and bool((values < 0).any()):
            raise ValueError(f"{operation} {label} contains negative values.")
        if binary and values.size and not bool(np.isin(values, [0, 1]).all()):
            raise ValueError(f"{operation} {label} must be binary.")
        return canonical

    def dense_array(matrix):
        return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)

    def numeric_summary(values):
        array = np.asarray(values, dtype=float)
        if array.size == 0 or not bool(np.isfinite(array).all()):
            raise RuntimeError(f"{operation} could not summarize an empty or non-finite diagnostic array.")
        return {
            "min": float(array.min()),
            "median": float(np.median(array)),
            "max": float(array.max()),
        }

    def representation_fingerprint(matrix, obs_names):
        digest = hashlib.sha256()
        digest.update(b"openbio-singlecell/milo-representation/v1\0")
        digest.update(json.dumps(list(matrix.shape), separators=(",", ":")).encode("ascii"))
        digest.update(json.dumps(list(obs_names), ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
        if sparse.issparse(matrix):
            canonical = sparse.csr_matrix(matrix, dtype=np.float64, copy=True)
            canonical.sum_duplicates()
            canonical.eliminate_zeros()
            canonical.sort_indices()
            digest.update(np.asarray(canonical.data, dtype="<f8").tobytes())
            digest.update(np.asarray(canonical.indices, dtype="<i8").tobytes())
            digest.update(np.asarray(canonical.indptr, dtype="<i8").tobytes())
        else:
            digest.update(np.ascontiguousarray(np.asarray(matrix, dtype="<f8")).tobytes())
        return digest.hexdigest()

    def require_parameters(callable_object, required, label):
        try:
            parameters = runtime_inspect.signature(callable_object).parameters
        except (TypeError, ValueError) as error:
            raise RuntimeError(f"{operation} cannot inspect the {label} interface.") from error
        missing = [name for name in required if name not in parameters]
        if missing:
            raise RuntimeError(
                f"{operation} requires the audited Pertpy/Scanpy interface for {label}; missing parameters: {missing}."
            )

    def distribution_version(distribution, module=None):
        try:
            value = importlib_metadata.version(distribution)
        except importlib_metadata.PackageNotFoundError:
            value = None
        if not value and module is not None:
            value = getattr(module, "__version__", None)
        if not isinstance(value, str) or not value.strip():
            raise RuntimeError(f"{operation} could not determine the installed {distribution} version.")
        return value.strip()

    def parsed_release(value, label):
        pieces = []
        for token in value.split("."):
            digits = "".join(character for character in token if character.isdigit())
            if not digits:
                break
            pieces.append(int(digits))
            if len(pieces) == 3:
                break
        if len(pieces) < 2:
            raise RuntimeError(f"{operation} cannot interpret {label} version {value!r}.")
        return tuple([*pieces, 0, 0][:3])

    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires an in-memory AnnData; call to_memory() first.")
    if int(getattr(adata, "n_obs", 0)) < 1 or int(getattr(adata, "n_vars", 0)) < 1:
        raise ValueError(f"{operation} requires non-empty observation and feature axes.")
    if not bool(adata.obs_names.is_unique) or not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique observation and feature identifiers.")
    obs_names = [canonical_text(value, "observation identifier") for value in adata.obs_names.tolist()]
    for value in adata.var_names.tolist():
        canonical_text(value, "feature identifier")

    sample_key = canonical_text(sample_key, "sample_key")
    condition_key = canonical_text(condition_key, "condition_key")
    reference_condition = canonical_text(reference_condition, "reference_condition")
    comparison_condition = canonical_text(comparison_condition, "comparison_condition")
    if reference_condition == comparison_condition:
        raise ValueError(f"{operation} reference_condition and comparison_condition must differ.")
    technical_batch_key = canonical_text(technical_batch_key, "technical_batch_key", allow_empty=True)
    annotation_key = canonical_text(annotation_key, "annotation_key")
    representation_key = canonical_text(representation_key, "representation_key")
    if annotation_status not in {"provisional", "curated"}:
        raise ValueError(f"{operation} annotation_status must be 'provisional' or 'curated'.")
    categorical_keys = parse_key_list(categorical_covariate_keys_json, "categorical_covariate_keys_json")
    continuous_keys = parse_key_list(continuous_covariate_keys_json, "continuous_covariate_keys_json")

    n_neighbors = numeric_value(n_neighbors, "n_neighbors", integer=True, minimum=2)
    random_seed = numeric_value(random_seed, "random_seed", integer=True, minimum=0, maximum=2**31 - 1)
    neighborhood_proportion = numeric_value(
        neighborhood_proportion,
        "neighborhood_proportion",
        minimum=0.0,
        maximum=1.0,
    )
    if neighborhood_proportion <= 0:
        raise ValueError(f"{operation} neighborhood_proportion must be greater than zero.")
    mixed_annotation_threshold = numeric_value(
        mixed_annotation_threshold,
        "mixed_annotation_threshold",
        minimum=0.0,
        maximum=1.0,
    )
    spatial_fdr_threshold = numeric_value(
        spatial_fdr_threshold,
        "spatial_fdr_threshold",
        minimum=0.0,
        maximum=1.0,
    )
    min_abs_log2_fold_change = numeric_value(
        min_abs_log2_fold_change,
        "min_abs_log2_fold_change",
        minimum=0.0,
    )

    categorical_roles = ([technical_batch_key] if technical_batch_key else []) + categorical_keys
    observation_roles = [sample_key, condition_key, annotation_key, *categorical_roles, *continuous_keys]
    duplicates = sorted({key for key in observation_roles if observation_roles.count(key) > 1})
    if duplicates:
        raise ValueError(f"{operation} observation roles must use distinct columns; repeated: {duplicates}.")
    missing_columns = [key for key in observation_roles if key not in adata.obs]
    if missing_columns:
        raise ValueError(f"{operation} observation columns not found: {missing_columns}.")
    if representation_key not in adata.obsm:
        raise ValueError(f"{operation} representation not found in obsm: {representation_key!r}.")

    all_sample_values = string_values(adata.obs[sample_key], f"obs[{sample_key!r}]")
    all_condition_values = string_values(adata.obs[condition_key], f"obs[{condition_key!r}]")
    sample_to_condition = {}
    for sample, condition in zip(all_sample_values, all_condition_values, strict=True):
        previous = sample_to_condition.setdefault(sample, condition)
        if previous != condition:
            raise ValueError(
                f"{operation} Sample {sample!r} maps to multiple Conditions: {previous!r} and {condition!r}."
            )
    present_conditions = list(dict.fromkeys(all_condition_values))
    missing_conditions = [
        condition for condition in (reference_condition, comparison_condition) if condition not in present_conditions
    ]
    if missing_conditions:
        raise ValueError(f"{operation} requested Conditions are absent: {missing_conditions}.")
    selected_mask = np.asarray(
        [condition in {reference_condition, comparison_condition} for condition in all_condition_values],
        dtype=bool,
    )
    selected_cells = int(selected_mask.sum())
    if selected_cells <= n_neighbors:
        raise ValueError(
            f"{operation} requires n_neighbors < selected cells ({n_neighbors} >= {selected_cells})."
        )
    if round(selected_cells * neighborhood_proportion) < 1:
        raise ValueError(f"{operation} neighborhood_proportion samples fewer than one candidate vertex.")

    # Milo needs an independently materialized two-Condition subset for graph and neighborhood mutation.
    selected = adata[selected_mask].copy()
    selected_obs_names = [obs_names[index] for index, keep in enumerate(selected_mask) if keep]
    if selected.obs_names.tolist() != selected_obs_names:
        raise RuntimeError(f"{operation} failed to preserve selected observation identity and order.")
    selected_sample_values = [all_sample_values[index] for index, keep in enumerate(selected_mask) if keep]
    selected_condition_values = [all_condition_values[index] for index, keep in enumerate(selected_mask) if keep]
    annotation_values = string_values(selected.obs[annotation_key], f"obs[{annotation_key!r}]")
    annotation_levels = [str(value) for value in pd.get_dummies(selected.obs[annotation_key]).columns.tolist()]
    if set(annotation_levels) != set(annotation_values):
        raise RuntimeError(f"{operation} failed to preserve annotation category identities.")

    sample_order = list(dict.fromkeys(selected_sample_values))
    samples_by_condition = {
        reference_condition: [sample for sample in sample_order if sample_to_condition[sample] == reference_condition],
        comparison_condition: [sample for sample in sample_order if sample_to_condition[sample] == comparison_condition],
    }
    insufficient = {condition: len(samples) for condition, samples in samples_by_condition.items() if len(samples) < 3}
    if insufficient:
        raise ValueError(
            f"{operation} requires at least 3 independent Samples per Condition; observed {insufficient}."
        )

    sample_cell_counts = {sample: selected_sample_values.count(sample) for sample in sample_order}
    condition_cell_counts = {
        reference_condition: selected_condition_values.count(reference_condition),
        comparison_condition: selected_condition_values.count(comparison_condition),
    }

    role_values = {}
    role_levels = {}
    for key in categorical_roles:
        values = string_values(selected.obs[key], f"obs[{key!r}]")
        levels = observed_levels(selected.obs[key], values)
        if len(levels) < 2:
            raise ValueError(f"{operation} categorical nuisance {key!r} must have at least two observed levels.")
        role_values[key] = values
        role_levels[key] = levels
    for key in continuous_keys:
        series = selected.obs[key]
        if bool(pd.api.types.is_bool_dtype(series.dtype)) or not bool(pd.api.types.is_numeric_dtype(series.dtype)):
            raise TypeError(f"{operation} continuous nuisance {key!r} must be numeric and non-boolean.")
        values = np.asarray(series, dtype=float)
        if not bool(np.isfinite(values).all()):
            raise ValueError(f"{operation} continuous nuisance {key!r} contains missing or non-finite values.")
        role_values[key] = values.tolist()

    sample_metadata = {sample: {condition_key: sample_to_condition[sample]} for sample in sample_order}
    selected_positions_by_sample = {}
    for position, sample in enumerate(selected_sample_values):
        selected_positions_by_sample.setdefault(sample, []).append(position)
    for key in [*categorical_roles, *continuous_keys]:
        for sample, positions in selected_positions_by_sample.items():
            unique = []
            for position in positions:
                value = role_values[key][position]
                if value not in unique:
                    unique.append(value)
            if len(unique) != 1:
                raise ValueError(
                    f"{operation} Sample {sample!r} maps to multiple values for nuisance {key!r}: {unique}."
                )
            sample_metadata[sample][key] = unique[0]

    representation = matrix_payload(
        selected.obsm[representation_key],
        f"obsm[{representation_key!r}]",
        rows=selected_cells,
    )
    if int(representation.shape[1]) < 2:
        raise ValueError(f"{operation} representation must contain at least two dimensions.")
    representation_dense_variance = np.asarray(
        representation.toarray().var(axis=0) if sparse.issparse(representation) else np.var(representation, axis=0)
    )
    if not bool((representation_dense_variance > 0).any()):
        raise ValueError(f"{operation} representation is degenerate across selected cells.")
    representation_sha256 = representation_fingerprint(representation, selected_obs_names)

    occupied = set(selected.obs.columns)
    sample_alias = safe_column("openbio_sample", occupied)
    condition_alias = safe_column("openbio_condition", occupied)
    selected.obs[sample_alias] = pd.Categorical(selected_sample_values, categories=sample_order, ordered=True)
    selected.obs[condition_alias] = pd.Categorical(
        ["reference" if value == reference_condition else "comparison" for value in selected_condition_values],
        categories=["reference", "comparison"],
        ordered=True,
    )

    encoded_categorical = []
    for position, key in enumerate(categorical_roles):
        base = "openbio_batch" if technical_batch_key and key == technical_batch_key else f"openbio_cat_{position}"
        alias = safe_column(base, occupied)
        levels = role_levels[key]
        encoded = [f"level_{levels.index(value)}" for value in role_values[key]]
        safe_levels = [f"level_{index}" for index in range(len(levels))]
        selected.obs[alias] = pd.Categorical(encoded, categories=safe_levels, ordered=True)
        encoded_categorical.append(
            {
                "source_key": key,
                "alias": alias,
                "levels": levels,
                "safe_levels": safe_levels,
                "reference_level": levels[0],
            }
        )

    encoded_continuous = []
    for position, key in enumerate(continuous_keys):
        alias = safe_column(f"openbio_cont_{position}", occupied)
        selected.obs[alias] = np.asarray(role_values[key], dtype=float)
        encoded_continuous.append({"source_key": key, "alias": alias})

    formula_terms = [condition_alias, *[item["alias"] for item in encoded_categorical]]
    formula_terms.extend(item["alias"] for item in encoded_continuous)
    design_formula = "~ " + " + ".join(formula_terms)
    model_contrast = f"{condition_alias}comparison-{condition_alias}reference"

    def build_design(rows):
        matrix_columns = [f"{condition_alias}reference", f"{condition_alias}comparison"]
        columns = [
            np.asarray([1.0 if sample_metadata[sample][condition_key] == reference_condition else 0.0 for sample in rows]),
            np.asarray(
                [1.0 if sample_metadata[sample][condition_key] == comparison_condition else 0.0 for sample in rows]
            ),
        ]
        categorical_evidence = []
        for item in encoded_categorical:
            key = item["source_key"]
            levels = item["levels"]
            if len(levels) >= len(rows):
                raise ValueError(
                    f"{operation} categorical nuisance {key!r} cannot use one distinct level per Sample."
                )
            for level_index, level in enumerate(levels[1:], start=1):
                matrix_columns.append(f"{item['alias']}level_{level_index}")
                columns.append(np.asarray([1.0 if sample_metadata[sample][key] == level else 0.0 for sample in rows]))
            categorical_evidence.append(
                {
                    "source_key": key,
                    "reference_level": levels[0],
                    "levels": list(levels),
                }
            )
        continuous_evidence = []
        for item in encoded_continuous:
            key = item["source_key"]
            values = np.asarray([sample_metadata[sample][key] for sample in rows], dtype=float)
            if bool(np.allclose(values, values[0], rtol=0.0, atol=0.0)):
                raise ValueError(f"{operation} continuous nuisance {key!r} is constant across Samples.")
            matrix_columns.append(item["alias"])
            columns.append(values)
            continuous_evidence.append(
                {
                    "source_key": key,
                    "min": float(values.min()),
                    "max": float(values.max()),
                }
            )
        matrix = np.column_stack(columns).astype(float, copy=False)
        rank = int(np.linalg.matrix_rank(matrix))
        residual_df = int(len(rows) - rank)
        if rank != int(matrix.shape[1]):
            raise ValueError(
                f"{operation} design is rank deficient ({rank} < {matrix.shape[1]}); "
                "Condition and nuisance covariates are confounded or redundant."
            )
        if residual_df <= 0:
            raise ValueError(
                f"{operation} design has no residual degrees of freedom ({len(rows)} Samples, rank {rank})."
            )
        contrast_vector = [-1.0, 1.0, *([0.0] * (matrix.shape[1] - 2))]
        return {
            "matrix": matrix,
            "columns": matrix_columns,
            "rank": rank,
            "residual_df": residual_df,
            "contrast_vector": contrast_vector,
            "categorical_covariates": categorical_evidence,
            "continuous_covariates": continuous_evidence,
        }

    build_design(sample_order)

    try:
        pertpy = importlib.import_module("pertpy")
    except (ImportError, OSError) as error:
        raise RuntimeError(
            f"{operation} requires pertpy 1.3.x with the Milo edgeR stack, rpy2, R, edgeR, limma and statmod."
        ) from error
    try:
        scanpy = importlib.import_module("scanpy")
    except (ImportError, OSError) as error:
        raise RuntimeError(f"{operation} requires Scanpy 1.12.x.") from error
    pertpy_version = distribution_version("pertpy", pertpy)
    if not ((1, 3, 0) <= parsed_release(pertpy_version, "Pertpy") < (1, 4, 0)):
        raise RuntimeError(
            f"{operation} supports the audited Pertpy >=1.3,<1.4 interface; found {pertpy_version}."
        )
    scanpy_version = distribution_version("scanpy", scanpy)
    if not ((1, 12, 0) <= parsed_release(scanpy_version, "Scanpy") < (1, 13, 0)):
        raise RuntimeError(
            f"{operation} supports the audited Scanpy >=1.12,<1.13 interface; found {scanpy_version}."
        )

    try:
        milo = pertpy.tl.Milo()
    except (AttributeError, TypeError) as error:
        raise RuntimeError(f"{operation} could not construct pertpy.tl.Milo from the audited interface.") from error
    require_parameters(milo.load, ["input", "feature_key"], "Milo.load")
    require_parameters(
        milo.make_nhoods,
        ["data", "neighbors_key", "feature_key", "prop", "seed", "copy"],
        "Milo.make_nhoods",
    )
    require_parameters(milo.count_nhoods, ["data", "sample_col", "feature_key"], "Milo.count_nhoods")
    require_parameters(
        milo.da_nhoods,
        [
            "mdata",
            "design",
            "model_contrasts",
            "subset_samples",
            "add_intercept",
            "feature_key",
            "reml",
            "max_iter",
            "tol",
            "solver",
        ],
        "Milo.da_nhoods",
    )
    require_parameters(milo.annotate_nhoods, ["mdata", "anno_col", "feature_key"], "Milo.annotate_nhoods")
    require_parameters(
        scanpy.pp.neighbors,
        [
            "adata",
            "n_neighbors",
            "n_pcs",
            "distances",
            "use_rep",
            "knn",
            "method",
            "transformer",
            "metric",
            "metric_kwds",
            "random_state",
            "key_added",
            "copy",
        ],
        "scanpy.pp.neighbors",
    )

    try:
        robjects = importlib.import_module("rpy2.robjects")

        def r_version(expression):
            result = robjects.r(expression)
            value = result[0]
            if not isinstance(value, str) or not value.strip():
                raise RuntimeError
            return value.strip()

        r_versions = {
            "R": r_version("as.character(getRversion())"),
            "edgeR": r_version('as.character(packageVersion("edgeR"))'),
            "limma": r_version('as.character(packageVersion("limma"))'),
            "statmod": r_version('as.character(packageVersion("statmod"))'),
        }
    except (ImportError, OSError, AttributeError, IndexError, TypeError, RuntimeError) as error:
        raise RuntimeError(
            f"{operation} could not verify R, edgeR, limma and statmod versions; "
            "install a working rpy2/R/Bioconductor edgeR stack before running Milo."
        ) from error

    mdata = milo.load(selected, feature_key=feature_key)
    try:
        rna = mdata[feature_key]
    except (KeyError, TypeError) as error:
        raise RuntimeError(f"{operation} Milo.load did not return a {feature_key!r} cell modality.") from error
    if rna.obs_names.tolist() != selected_obs_names:
        raise RuntimeError(f"{operation} Milo.load changed the selected observation identity or order.")
    loaded_representation = matrix_payload(
        rna.obsm.get(representation_key),
        f"loaded obsm[{representation_key!r}]",
        rows=selected_cells,
        columns=int(representation.shape[1]),
    )
    if representation_fingerprint(loaded_representation, selected_obs_names) != representation_sha256:
        raise RuntimeError(f"{operation} Milo.load changed the selected representation.")

    graph_key = "_openbio_milo_graph"
    graph_suffix = 1
    while (
        graph_key in rna.uns
        or f"{graph_key}_distances" in rna.obsp
        or f"{graph_key}_connectivities" in rna.obsp
    ):
        graph_key = f"_openbio_milo_graph_{graph_suffix}"
        graph_suffix += 1

    with _MILO_RNG_LOCK:
        python_rng_state = random.getstate()
        numpy_rng_state = np.random.get_state()
        try:
            neighbor_return = scanpy.pp.neighbors(
                rna,
                n_neighbors=n_neighbors,
                n_pcs=None,
                distances=None,
                use_rep=representation_key,
                knn=True,
                method="umap",
                transformer="pynndescent",
                metric="euclidean",
                metric_kwds={},
                random_state=random_seed,
                key_added=graph_key,
                copy=False,
            )
            if neighbor_return is not None:
                raise RuntimeError(f"{operation} Scanpy neighbors violated copy=False in-place semantics.")
            nhood_return = milo.make_nhoods(
                mdata,
                neighbors_key=graph_key,
                feature_key=feature_key,
                prop=neighborhood_proportion,
                seed=random_seed,
                copy=False,
            )
            if nhood_return is not None:
                raise RuntimeError(f"{operation} Milo.make_nhoods violated copy=False in-place semantics.")
        finally:
            random.setstate(python_rng_state)
            np.random.set_state(numpy_rng_state)

    try:
        graph_metadata = rna.uns[graph_key]
        distances_key = f"{graph_key}_distances"
        connectivities_key = f"{graph_key}_connectivities"
        graph_distances = rna.obsp[distances_key]
        graph_connectivities = rna.obsp[connectivities_key]
    except KeyError as error:
        raise RuntimeError(f"{operation} Scanpy did not create the complete private neighbor graph bundle.") from error
    if not isinstance(graph_metadata, dict) or not isinstance(graph_metadata.get("params"), dict):
        raise RuntimeError(f"{operation} private neighbor graph metadata is malformed.")
    if (
        graph_metadata.get("distances_key") != distances_key
        or graph_metadata.get("connectivities_key") != connectivities_key
    ):
        raise RuntimeError(f"{operation} private neighbor graph bundle key metadata is inconsistent.")
    graph_params = graph_metadata["params"]
    if (
        graph_params.get("use_rep") != representation_key
        or graph_params.get("n_neighbors") != n_neighbors
        or graph_params.get("method") != "umap"
        or graph_params.get("metric") != "euclidean"
        or graph_params.get("random_state") != random_seed
    ):
        raise RuntimeError(f"{operation} private neighbor graph metadata does not match the requested representation/k.")
    graph_distances = matrix_payload(
        graph_distances,
        "neighbor distances",
        rows=selected_cells,
        columns=selected_cells,
        nonnegative=True,
    )
    graph_connectivities = matrix_payload(
        graph_connectivities,
        "neighbor connectivities",
        rows=selected_cells,
        columns=selected_cells,
        nonnegative=True,
    )
    if rna.uns.get("nhood_neighbors_key") != graph_key:
        raise RuntimeError(f"{operation} Milo neighborhoods are not linked to the private neighbor graph.")

    if "nhoods" not in rna.obsm:
        raise RuntimeError(f"{operation} Milo.make_nhoods did not create obsm['nhoods'].")
    membership = matrix_payload(rna.obsm["nhoods"], "neighborhood membership", rows=selected_cells, binary=True)
    membership = sparse.csr_matrix(membership)
    neighborhood_count = int(membership.shape[1])
    if neighborhood_count < 2:
        raise RuntimeError(f"{operation} requires at least two refined neighborhoods; observed {neighborhood_count}.")
    neighborhood_sizes = np.asarray(membership.sum(axis=0)).ravel().astype(int)
    if bool((neighborhood_sizes <= 0).any()):
        raise RuntimeError(f"{operation} produced an empty neighborhood.")
    membership_by_neighborhood = membership.tocsc()
    membership_signatures = []
    for column in range(neighborhood_count):
        start, stop = membership_by_neighborhood.indptr[column : column + 2]
        membership_signatures.append(tuple(membership_by_neighborhood.indices[start:stop].tolist()))
    if len(set(membership_signatures)) != neighborhood_count:
        raise RuntimeError(f"{operation} requires unique neighborhood memberships after refinement.")
    for column in ("nhood_ixs_random", "nhood_ixs_refined", "nhood_kth_distance"):
        if column not in rna.obs:
            raise RuntimeError(f"{operation} Milo.make_nhoods did not create obs[{column!r}].")
    random_flags = np.asarray(rna.obs["nhood_ixs_random"])
    refined_flags = np.asarray(rna.obs["nhood_ixs_refined"])
    if not bool(np.isin(random_flags, [0, 1]).all()) or not bool(np.isin(refined_flags, [0, 1]).all()):
        raise RuntimeError(f"{operation} neighborhood index flags must be binary.")
    expected_random_indices = int(np.round(selected_cells * neighborhood_proportion))
    if int((random_flags == 1).sum()) != expected_random_indices:
        raise RuntimeError(
            f"{operation} random candidate index count does not match round(n_obs * neighborhood_proportion)."
        )
    refined_positions = np.flatnonzero(refined_flags == 1)
    if int(refined_positions.size) != neighborhood_count:
        raise RuntimeError(f"{operation} refined index count does not match the neighborhood axis.")
    expected_membership = sparse.csr_matrix(graph_connectivities)[:, refined_positions].copy()
    expected_membership.data = np.ones_like(expected_membership.data)
    expected_membership.eliminate_zeros()
    if (membership != expected_membership).nnz:
        raise RuntimeError(f"{operation} neighborhood membership does not match the private graph and refined indices.")
    expected_index_cells = [selected_obs_names[position] for position in refined_positions]
    kth_all = np.asarray(rna.obs["nhood_kth_distance"], dtype=float)
    kth_distances = kth_all[refined_positions]
    if not bool(np.isfinite(kth_distances).all()) or bool((kth_distances <= 0).any()):
        raise RuntimeError(f"{operation} requires finite positive kth-neighbor distances for SpatialFDR.")
    if sparse.issparse(graph_distances):
        graph_kth_distances = np.asarray(graph_distances[refined_positions].max(axis=1).toarray()).ravel()
    else:
        graph_kth_distances = np.asarray(graph_distances)[refined_positions].max(axis=1)
    if not bool(np.allclose(kth_distances, graph_kth_distances, rtol=0.0, atol=1e-12)):
        raise RuntimeError(f"{operation} kth-neighbor distances do not match the private graph and refined indices.")

    counted = milo.count_nhoods(mdata, sample_col=sample_alias, feature_key=feature_key)
    if counted is None:
        raise RuntimeError(f"{operation} Milo.count_nhoods did not return the documented MuData result.")
    mdata = counted
    try:
        rna = mdata[feature_key]
        sample_adata = mdata["milo"]
    except (KeyError, TypeError) as error:
        raise RuntimeError(f"{operation} count_nhoods did not return rna and milo modalities.") from error
    if rna.obs_names.tolist() != selected_obs_names:
        raise RuntimeError(f"{operation} count_nhoods changed cell identity or order.")
    counted_membership = matrix_payload(
        rna.obsm.get("nhoods"),
        "count_nhoods neighborhood membership",
        rows=selected_cells,
        columns=neighborhood_count,
        binary=True,
    )
    if (sparse.csr_matrix(counted_membership) != membership).nnz:
        raise RuntimeError(f"{operation} count_nhoods changed neighborhood membership.")
    if int(sample_adata.n_vars) != neighborhood_count:
        raise RuntimeError(f"{operation} Sample-count neighborhood axis does not match membership.")
    if not bool(sample_adata.obs_names.is_unique) or not bool(sample_adata.var_names.is_unique):
        raise RuntimeError(f"{operation} Sample-count axes must be unique.")
    counted_samples = [canonical_text(value, "counted Sample identifier") for value in sample_adata.obs_names.tolist()]
    if set(counted_samples) != set(sample_order) or len(counted_samples) != len(sample_order):
        raise RuntimeError(f"{operation} Sample-count rows do not match the selected independent Samples.")
    count_matrix = matrix_payload(
        sample_adata.X,
        "Sample-by-neighborhood counts",
        rows=len(sample_order),
        columns=neighborhood_count,
        nonnegative=True,
    )
    count_values = count_matrix.data if sparse.issparse(count_matrix) else np.asarray(count_matrix).ravel()
    if count_values.size and not bool(np.allclose(count_values, np.rint(count_values), rtol=0.0, atol=1e-8)):
        raise RuntimeError(f"{operation} Sample-by-neighborhood counts must be integers.")
    count_dense = dense_array(count_matrix).astype(float, copy=True)
    if bool((count_dense.sum(axis=1) <= 0).any()):
        raise RuntimeError(f"{operation} every selected Sample must have a nonzero neighborhood library.")
    sample_dummy = pd.get_dummies(pd.Categorical(selected_sample_values, categories=sample_order)).to_numpy(dtype=int)
    expected_counts = np.asarray(membership.T.dot(sample_dummy)).T
    expected_by_sample = {sample: expected_counts[position] for position, sample in enumerate(sample_order)}
    aligned_expected = np.vstack([expected_by_sample[sample] for sample in counted_samples])
    if not bool(np.array_equal(count_dense, aligned_expected)):
        raise RuntimeError(f"{operation} Sample counts do not equal the validated cell-neighborhood memberships.")
    for column in ("index_cell", "kth_distance"):
        if column not in sample_adata.var:
            raise RuntimeError(f"{operation} count_nhoods did not preserve neighborhood {column!r}.")
    counted_index_cells = string_values(sample_adata.var["index_cell"], "neighborhood index-cell identity")
    if counted_index_cells != expected_index_cells:
        raise RuntimeError(f"{operation} count_nhoods changed neighborhood index-cell identity.")
    if not bool(pd.api.types.is_numeric_dtype(sample_adata.var["kth_distance"].dtype)) or bool(
        pd.api.types.is_bool_dtype(sample_adata.var["kth_distance"].dtype)
    ):
        raise TypeError(f"{operation} count_nhoods kth-neighbor distances must be numeric.")
    counted_kth = np.asarray(sample_adata.var["kth_distance"], dtype=float)
    if not bool(np.array_equal(counted_kth, kth_distances)):
        raise RuntimeError(f"{operation} count_nhoods changed kth-neighbor distances.")

    design = build_design(counted_samples)
    original_neighborhood_axis = sample_adata.var_names.tolist()
    da_return = milo.da_nhoods(
        mdata,
        design=design_formula,
        model_contrasts=model_contrast,
        add_intercept=False,
        feature_key=feature_key,
        reml=True,
        max_iter=50,
        tol=1e-5,
        solver="edger",
    )
    if da_return is not None:
        raise RuntimeError(f"{operation} Milo.da_nhoods violated in-place semantics.")
    if sample_adata.obs_names.tolist() != counted_samples or sample_adata.var_names.tolist() != original_neighborhood_axis:
        raise RuntimeError(f"{operation} da_nhoods changed Sample or neighborhood identity.")
    fitted_counts = matrix_payload(
        sample_adata.X,
        "post-fit Sample-by-neighborhood counts",
        rows=len(sample_order),
        columns=neighborhood_count,
        nonnegative=True,
    )
    if not bool(np.array_equal(dense_array(fitted_counts), count_dense)):
        raise RuntimeError(f"{operation} da_nhoods changed Sample-by-neighborhood counts.")
    annotation_return = milo.annotate_nhoods(mdata, anno_col=annotation_key, feature_key=feature_key)
    if annotation_return is not None:
        raise RuntimeError(f"{operation} Milo.annotate_nhoods violated in-place semantics.")
    if sample_adata.var_names.tolist() != original_neighborhood_axis:
        raise RuntimeError(f"{operation} annotate_nhoods changed neighborhood identity.")

    required_result_columns = [
        "logFC",
        "logCPM",
        "F",
        "PValue",
        "FDR",
        "SpatialFDR",
        "nhood_annotation",
        "nhood_annotation_frac",
    ]
    missing_results = [column for column in required_result_columns if column not in sample_adata.var]
    if missing_results:
        raise RuntimeError(f"{operation} backend result is missing columns: {missing_results}.")
    result_arrays = {}
    for column in required_result_columns[:6]:
        series = sample_adata.var[column]
        if not bool(pd.api.types.is_numeric_dtype(series.dtype)) or bool(pd.api.types.is_bool_dtype(series.dtype)):
            raise TypeError(f"{operation} result column {column!r} must be numeric.")
        values = np.asarray(series, dtype=float)
        if values.shape != (neighborhood_count,) or not bool(np.isfinite(values).all()):
            raise RuntimeError(f"{operation} result column {column!r} must be complete and finite.")
        result_arrays[column] = values
    if bool((result_arrays["F"] < 0).any()):
        raise RuntimeError(f"{operation} quasi-likelihood F statistics cannot be negative.")
    for column in ("PValue", "FDR", "SpatialFDR"):
        if bool(((result_arrays[column] < 0) | (result_arrays[column] > 1)).any()):
            raise RuntimeError(f"{operation} result column {column!r} must lie in [0, 1].")

    p_values = result_arrays["PValue"]
    ordinary_order = np.argsort(p_values, kind="mergesort")
    ordinary_sorted = p_values[ordinary_order]
    ordinary_adjusted_sorted = np.minimum.accumulate(
        (ordinary_sorted * neighborhood_count / np.arange(1, neighborhood_count + 1))[::-1]
    )[::-1]
    ordinary_expected = np.empty_like(p_values)
    ordinary_expected[ordinary_order] = np.minimum(ordinary_adjusted_sorted, 1.0)
    if not bool(np.allclose(result_arrays["FDR"], ordinary_expected, rtol=1e-10, atol=1e-12)):
        raise RuntimeError(f"{operation} failed independent ordinary-BH FDR verification.")

    spatial_weights = 1.0 / kth_distances
    spatial_order = np.argsort(p_values)
    spatial_p_sorted = p_values[spatial_order]
    spatial_weight_sorted = spatial_weights[spatial_order]
    spatial_adjusted_sorted = np.minimum.accumulate(
        (spatial_weights.sum() * spatial_p_sorted / np.cumsum(spatial_weight_sorted))[::-1]
    )[::-1]
    spatial_expected = np.empty_like(p_values)
    spatial_expected[spatial_order] = np.minimum(spatial_adjusted_sorted, 1.0)
    if not bool(np.allclose(result_arrays["SpatialFDR"], spatial_expected, rtol=1e-10, atol=1e-12)):
        raise RuntimeError(f"{operation} failed independent SpatialFDR verification.")

    annotation_dummy = pd.get_dummies(rna.obs[annotation_key])
    if [str(value) for value in annotation_dummy.columns.tolist()] != annotation_levels:
        raise RuntimeError(f"{operation} annotation category order changed before neighborhood annotation.")
    annotation_counts = np.asarray(membership.T.dot(annotation_dummy.to_numpy(dtype=float)))
    annotation_totals = annotation_counts.sum(axis=1)
    if bool((annotation_totals <= 0).any()):
        raise RuntimeError(f"{operation} produced a neighborhood without annotated cells.")
    expected_fractions = annotation_counts / annotation_totals[:, None]
    expected_majority_positions = expected_fractions.argmax(axis=1)
    expected_majority = [annotation_levels[position] for position in expected_majority_positions]
    expected_majority_fraction = expected_fractions.max(axis=1)
    majority = string_values(sample_adata.var["nhood_annotation"], "neighborhood majority annotation")
    if not bool(pd.api.types.is_numeric_dtype(sample_adata.var["nhood_annotation_frac"].dtype)) or bool(
        pd.api.types.is_bool_dtype(sample_adata.var["nhood_annotation_frac"].dtype)
    ):
        raise TypeError(f"{operation} neighborhood majority annotation fractions must be numeric.")
    majority_fraction = np.asarray(sample_adata.var["nhood_annotation_frac"], dtype=float)
    if (
        majority_fraction.shape != (neighborhood_count,)
        or not bool(np.isfinite(majority_fraction).all())
        or bool(((majority_fraction < 0) | (majority_fraction > 1)).any())
    ):
        raise RuntimeError(f"{operation} neighborhood majority annotation fractions must be complete and lie in [0, 1].")
    if majority != expected_majority or not bool(
        np.allclose(majority_fraction, expected_majority_fraction, rtol=0.0, atol=1e-12)
    ):
        raise RuntimeError(f"{operation} neighborhood majority annotation does not match membership evidence.")
    if "frac_annotation" not in sample_adata.varm:
        raise RuntimeError(f"{operation} annotate_nhoods did not preserve the full annotation fraction matrix.")
    stored_fractions_raw = np.asarray(sample_adata.varm["frac_annotation"])
    if not bool(np.issubdtype(stored_fractions_raw.dtype, np.number)) or bool(
        np.issubdtype(stored_fractions_raw.dtype, np.bool_)
    ):
        raise TypeError(f"{operation} annotation fraction matrix must be numeric and non-boolean.")
    stored_fractions = np.asarray(stored_fractions_raw, dtype=float)
    stored_labels = [str(value) for value in sample_adata.uns.get("annotation_labels", [])]
    if (
        sample_adata.uns.get("annotation_obs") != annotation_key
        or stored_labels != annotation_levels
        or stored_fractions.shape != expected_fractions.shape
    ):
        raise RuntimeError(f"{operation} annotation fraction axes are malformed.")
    if not bool(np.isfinite(stored_fractions).all()) or bool(((stored_fractions < 0) | (stored_fractions > 1)).any()):
        raise RuntimeError(f"{operation} annotation fraction matrix must be complete and lie in [0, 1].")
    if not bool(np.allclose(stored_fractions, expected_fractions, rtol=0.0, atol=1e-12)):
        raise RuntimeError(f"{operation} annotation fraction matrix does not match membership evidence.")

    is_mixed = majority_fraction < mixed_annotation_threshold
    reported_annotation = ["Mixed" if mixed else label for label, mixed in zip(majority, is_mixed, strict=True)]
    neighborhood_ids = [f"nhood_{index + 1:06d}" for index in range(neighborhood_count)]
    table = pd.DataFrame(
        {
            "neighborhood_id": neighborhood_ids,
            "index_cell": expected_index_cells,
            "neighborhood_size": neighborhood_sizes.astype(int),
            "kth_distance": kth_distances.astype(float),
            "reference_condition": [reference_condition] * neighborhood_count,
            "comparison_condition": [comparison_condition] * neighborhood_count,
            "log2_fold_change": result_arrays["logFC"],
            "log_counts_per_million": result_arrays["logCPM"],
            "quasi_likelihood_f": result_arrays["F"],
            "p_value": result_arrays["PValue"],
            "p_adjusted_bh": result_arrays["FDR"],
            "spatial_fdr": result_arrays["SpatialFDR"],
            "majority_annotation": majority,
            "majority_annotation_fraction": majority_fraction,
            "reported_annotation": reported_annotation,
            "is_mixed": is_mixed.astype(bool),
        },
        columns=[
            "neighborhood_id",
            "index_cell",
            "neighborhood_size",
            "kth_distance",
            "reference_condition",
            "comparison_condition",
            "log2_fold_change",
            "log_counts_per_million",
            "quasi_likelihood_f",
            "p_value",
            "p_adjusted_bh",
            "spatial_fdr",
            "majority_annotation",
            "majority_annotation_fraction",
            "reported_annotation",
            "is_mixed",
        ],
    )
    if table["neighborhood_id"].duplicated().any() or table["index_cell"].duplicated().any():
        raise RuntimeError(f"{operation} canonical neighborhood identities must be unique.")

    calls = (table["spatial_fdr"] <= spatial_fdr_threshold) & (
        table["log2_fold_change"].abs() >= min_abs_log2_fold_change
    )
    comparison_calls = calls & (table["log2_fold_change"] > 0)
    reference_calls = calls & (table["log2_fold_change"] < 0)

    def result_records(mask, *, ascending):
        subset = table.loc[mask].sort_values(
            ["spatial_fdr", "log2_fold_change", "neighborhood_id"],
            ascending=[True, ascending, True],
            kind="mergesort",
        )
        records = []
        for row in subset.head(10).itertuples(index=False):
            records.append(
                {
                    "neighborhood_id": row.neighborhood_id,
                    "index_cell": row.index_cell,
                    "log2_fold_change": float(row.log2_fold_change),
                    "spatial_fdr": float(row.spatial_fdr),
                    "majority_annotation": row.majority_annotation,
                    "majority_annotation_fraction": float(row.majority_annotation_fraction),
                    "neighborhood_size": int(row.neighborhood_size),
                }
            )
        return records

    def python_module_version(distribution, module_name):
        module = None
        try:
            module = importlib.import_module(module_name)
        except (ImportError, OSError):
            pass
        return distribution_version(distribution, module)

    software_versions = {
        "python": platform.python_version(),
        "openbio-singlecell": str(openbio_version),
        "pertpy": pertpy_version,
        "scanpy": scanpy_version,
        "anndata": python_module_version("anndata", "anndata"),
        "mudata": python_module_version("mudata", "mudata"),
        "numpy": python_module_version("numpy", "numpy"),
        "pandas": python_module_version("pandas", "pandas"),
        "scipy": python_module_version("scipy", "scipy"),
        "scikit-learn": python_module_version("scikit-learn", "sklearn"),
        "pynndescent": python_module_version("pynndescent", "pynndescent"),
        "rpy2": python_module_version("rpy2", "rpy2"),
        **r_versions,
    }

    excluded_conditions = [condition for condition in present_conditions if condition not in {reference_condition, comparison_condition}]
    warnings = []
    if excluded_conditions:
        warnings.append(
            f"Excluded {int((~selected_mask).sum())} cells from non-compared Conditions {excluded_conditions} before graph construction."
        )
    if annotation_status == "provisional":
        warnings.append(
            "Neighborhood labels are provisional annotations and must not be presented as curated cell-type truth."
        )
    condition_cell_imbalance_ratio = float(max(condition_cell_counts.values()) / min(condition_cell_counts.values()))
    sample_cell_imbalance_ratio = float(max(sample_cell_counts.values()) / min(sample_cell_counts.values()))
    if condition_cell_imbalance_ratio >= 4.0 or sample_cell_imbalance_ratio >= 4.0:
        warnings.append(
            "Selected cell counts are at least 4-fold imbalanced across Conditions or Samples; "
            "review representation coverage and neighborhood-count sensitivity even though biological replication is adequate."
        )
    parameters = {
        "sample_key": sample_key,
        "condition_key": condition_key,
        "reference_condition": reference_condition,
        "comparison_condition": comparison_condition,
        "technical_batch_key": technical_batch_key,
        "categorical_covariate_keys": categorical_keys,
        "continuous_covariate_keys": continuous_keys,
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "representation_key": representation_key,
        "n_neighbors": n_neighbors,
        "neighborhood_proportion": neighborhood_proportion,
        "mixed_annotation_threshold": mixed_annotation_threshold,
        "spatial_fdr_threshold": spatial_fdr_threshold,
        "min_abs_log2_fold_change": min_abs_log2_fold_change,
        "random_seed": random_seed,
        "neighbor_metric": "euclidean",
        "neighbor_connectivity_method": "umap",
        "neighbor_transformer": "pynndescent",
        "solver": "edger",
        "normalization": "TMM",
        "test": "robust edgeR quasi-likelihood F",
        "multiple_testing_primary": "k-distance weighted SpatialFDR",
        "cell_count_imbalance_warning_ratio": 4.0,
    }
    sample_metadata_rows = []
    for sample in counted_samples:
        sample_metadata_rows.append(
            {
                "sample_id": sample,
                "condition": sample_metadata[sample][condition_key],
                "categorical_covariates": {
                    key: sample_metadata[sample][key] for key in categorical_roles
                },
                "continuous_covariates": {
                    key: float(sample_metadata[sample][key]) for key in continuous_keys
                },
            }
        )
    condition_sample_coverage = {}
    for condition in (reference_condition, comparison_condition):
        positions = [index for index, sample in enumerate(counted_samples) if sample_to_condition[sample] == condition]
        condition_sample_coverage[condition] = numeric_summary((count_dense[positions] > 0).sum(axis=0))
    analysis_summary = {
        "tested_neighborhoods": neighborhood_count,
        "comparison_enriched_calls": int(comparison_calls.sum()),
        "reference_enriched_calls": int(reference_calls.sum()),
        "mixed_annotation_neighborhoods": int(is_mixed.sum()),
        "ordinary_bh_is_primary": False,
        "spatial_fdr_threshold": spatial_fdr_threshold,
        "min_abs_log2_fold_change": min_abs_log2_fold_change,
    }
    results_text = (
        f"Milo tested {neighborhood_count} overlapping neighborhoods using {len(sample_order)} independent Samples "
        f"for {comparison_condition} versus {reference_condition}. At SpatialFDR <= {spatial_fdr_threshold:g} "
        f"and |log2 fold change| >= {min_abs_log2_fold_change:g}, "
        f"{analysis_summary['comparison_enriched_calls']} neighborhoods were enriched in {comparison_condition} and "
        f"{analysis_summary['reference_enriched_calls']} were enriched in {reference_condition}."
    )
    methods_text = (
        "Pairwise Milo overlapping K-nearest-neighbor differential abundance used independent Sample neighborhood "
        "cell counts, TMM-normalized negative-binomial edgeR robust quasi-likelihood F tests, and k-distance weighted "
        "SpatialFDR over the complete neighborhood universe. Neighborhood labels are descriptive majority annotations "
        "with an explicit mixed-label threshold."
    )
    summary = {
        "schema_version": 1,
        "schema_id": "openbio-singlecell/milo-differential-abundance/v1",
        "node_id": "OpenBioSingleCellMiloDifferentialAbundance",
        "analysis_status": "sample_level_condition_inference",
        "methods": methods_text,
        "results": results_text,
        "comparison": {
            "condition_key": condition_key,
            "reference": reference_condition,
            "comparison": comparison_condition,
            "positive_log2_fold_change": f"{comparison_condition} enriched versus {reference_condition}",
        },
        "input_evidence": {
            "input_cells": int(adata.n_obs),
            "selected_cells": selected_cells,
            "excluded_cells": int(adata.n_obs - selected_cells),
            "excluded_conditions": excluded_conditions,
            "independent_samples": len(sample_order),
            "sample_order": sample_order,
            "sample_cell_counts": sample_cell_counts,
            "condition_cell_counts": condition_cell_counts,
            "condition_cell_imbalance_ratio": condition_cell_imbalance_ratio,
            "sample_cell_imbalance_ratio": sample_cell_imbalance_ratio,
            "annotation_key": annotation_key,
            "annotation_status": annotation_status,
            "representation": {
                "key": representation_key,
                "shape": [int(representation.shape[0]), int(representation.shape[1])],
                "sha256": representation_sha256,
                "state": "precomputed representation used to derive the private comparison-specific graph",
            },
        },
        "design_evidence": {
            "sample_unit": "independent biological Sample",
            "sample_key": sample_key,
            "samples_by_condition": {
                reference_condition: len(samples_by_condition[reference_condition]),
                comparison_condition: len(samples_by_condition[comparison_condition]),
            },
            "sample_order": counted_samples,
            "sample_metadata": sample_metadata_rows,
            "formula": design_formula,
            "columns": design["columns"],
            "rank": design["rank"],
            "residual_degrees_of_freedom": design["residual_df"],
            "model_contrast": model_contrast,
            "contrast_vector": design["contrast_vector"],
            "categorical_covariates": design["categorical_covariates"],
            "continuous_covariates": design["continuous_covariates"],
            "minimum_samples_per_condition_policy": 3,
            "full_rank_verified": True,
            "contrast_estimable": True,
        },
        "graph_evidence": {
            "neighbors_key_scope": "private per-execution key",
            "n_neighbors": n_neighbors,
            "metric": "euclidean",
            "connectivity_method": "umap",
            "transformer": "pynndescent",
            "random_seed": random_seed,
            "requested_candidate_vertices": int(round(selected_cells * neighborhood_proportion)),
            "refined_neighborhoods": neighborhood_count,
            "neighborhood_size": numeric_summary(neighborhood_sizes),
            "kth_distance": numeric_summary(kth_distances),
            "samples_represented_per_neighborhood": numeric_summary((count_dense > 0).sum(axis=0)),
            "samples_represented_per_neighborhood_by_condition": condition_sample_coverage,
            "selected_cell_coverage": float((np.asarray(membership.sum(axis=1)).ravel() > 0).mean()),
        },
        "key_results": {
            **analysis_summary,
            "comparison_enriched": result_records(comparison_calls, ascending=False),
            "reference_enriched": result_records(reference_calls, ascending=True),
        },
        "parameters": parameters,
        "warnings": warnings,
        "limitations": [
            "Milo findings are overlapping local graph regions, not independent cell populations or cell-type-wide tests.",
            "The result depends on the supplied representation, k, neighborhood proportion and random seed.",
            "Three independent Samples per Condition is a conservative minimum policy, not a power guarantee.",
            "Cell counts do not increase the biological replicate n; the model rows are independent Samples.",
            "Paired, longitudinal or repeated-donor Samples require a dedicated repeated-measures model and are outside this independent-Sample interface.",
            "Neighborhood annotation is descriptive and does not enter the abundance model.",
            "Batch correction in the representation and Sample-level nuisance adjustment cannot guarantee removal of all technical effects.",
            "Observed Condition associations are not causal effects without an appropriate experimental or causal design.",
            "SpatialFDR is the primary Milo multiplicity control; ordinary BH is retained only as a diagnostic.",
        ],
        "references": [
            {
                "citation": (
                    "Dann E, Henderson NC, Teichmann SA, Morgan MD, Marioni JC. Differential abundance testing "
                    "on single-cell data using K-nearest neighbor graphs. Nature Biotechnology. 2022;40:245-253."
                ),
                "doi": "10.1038/s41587-021-01033-z",
                "url": "https://doi.org/10.1038/s41587-021-01033-z",
                "kind": "method",
            },
            {
                "citation": "Heumos L, et al. Pertpy: an end-to-end framework for perturbation analysis. Nature Methods. 2026.",
                "doi": "10.1038/s41592-025-02909-7",
                "url": "https://doi.org/10.1038/s41592-025-02909-7",
                "kind": "software",
            },
            {
                "citation": "Robinson MD, McCarthy DJ, Smyth GK. edgeR. Bioinformatics. 2010;26:139-140.",
                "doi": "10.1093/bioinformatics/btp616",
                "url": "https://doi.org/10.1093/bioinformatics/btp616",
                "kind": "software",
            },
            {
                "citation": "Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis. Genome Biology. 2010;11:R25.",
                "doi": "10.1186/gb-2010-11-3-r25",
                "url": "https://doi.org/10.1186/gb-2010-11-3-r25",
                "kind": "method",
            },
            {
                "citation": "Lun ATL, Chen Y, Smyth GK. It's DE-licious: a recipe for differential expression analyses of RNA-seq experiments using quasi-likelihood methods in edgeR. 2016.",
                "doi": "10.1007/978-1-4939-3578-9_19",
                "url": "https://doi.org/10.1007/978-1-4939-3578-9_19",
                "kind": "method",
            },
            {
                "citation": "Pertpy developers. pertpy.tools.Milo official API documentation.",
                "url": "https://pertpy.readthedocs.io/en/latest/api/tools/pertpy.tools.Milo.html",
                "kind": "software_documentation",
            },
            {
                "citation": (
                    "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. "
                    "Genome Biology. 2018;19:15."
                ),
                "doi": "10.1186/s13059-017-1382-0",
                "url": "https://doi.org/10.1186/s13059-017-1382-0",
                "kind": "software",
            },
            {
                "citation": (
                    "McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for "
                    "Dimension Reduction. Journal of Open Source Software. 2018;3:861."
                ),
                "doi": "10.21105/joss.00861",
                "url": "https://doi.org/10.21105/joss.00861",
                "kind": "method",
            },
            {
                "citation": (
                    "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful "
                    "approach to multiple testing. Journal of the Royal Statistical Society B. 1995;57:289-300."
                ),
                "doi": "10.1111/j.2517-6161.1995.tb02031.x",
                "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
                "kind": "method",
            },
        ],
        "software_versions": software_versions,
    }
    row_memberships = np.diff(membership.indptr)
    overlap_pair_contributions = sum(int(value) * (int(value) - 1) // 2 for value in row_memberships)
    membership_bytes = int(membership.data.nbytes + membership.indices.nbytes + membership.indptr.nbytes)
    overlap_working_bytes = membership_bytes * 3 + overlap_pair_contributions * 64
    # ponytail: fixed artifact bound; make this a node parameter only if larger Milo graphs become a demonstrated need.
    maximum_overlap_pair_contributions = 10_000_000
    maximum_overlap_working_bytes = 2 * 1024**3
    if overlap_pair_contributions > maximum_overlap_pair_contributions:
        raise MemoryError(
            "Milo overlap-graph pair contributions exceed the fixed artifact limit: "
            f"{overlap_pair_contributions:,} > {maximum_overlap_pair_contributions:,}."
        )
    if overlap_working_bytes > maximum_overlap_working_bytes:
        raise MemoryError(
            "Milo overlap-graph estimated working memory exceeds the fixed artifact limit: "
            f"{overlap_working_bytes / 1024**3:.3f} GiB > "
            f"{maximum_overlap_working_bytes / 1024**3:.3f} GiB."
        )
    neighborhood_graph = (membership.T @ membership).tocsr()
    neighborhood_graph.setdiag(0)
    neighborhood_graph.eliminate_zeros()
    neighborhood_graph.sort_indices()
    summary["graph_evidence"]["overlap_graph_preflight"] = {
        "pair_contributions": overlap_pair_contributions,
        "estimated_working_bytes": overlap_working_bytes,
        "maximum_pair_contributions": maximum_overlap_pair_contributions,
        "maximum_working_bytes": maximum_overlap_working_bytes,
    }
    evidence = {
        "membership": membership,
        "graph": neighborhood_graph,
        "representative_coordinates": np.asarray(
            dense_array(representation[refined_positions, :2]),
            dtype=float,
        ),
        "observation_names": list(selected_obs_names),
        "neighborhood_names": list(neighborhood_ids),
        "representation_key": representation_key,
        "representation_sha256": representation_sha256,
        "condition_key": condition_key,
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "n_neighbors": n_neighbors,
        "neighborhood_proportion": neighborhood_proportion,
        "mixed_annotation_threshold": mixed_annotation_threshold,
        "spatial_fdr_threshold": spatial_fdr_threshold,
        "min_abs_log2_fold_change": min_abs_log2_fold_change,
        "random_seed": random_seed,
    }
    json.dumps(summary, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)
    return table, summary, evidence


def _standalone_milo_differential_abundance(adata, **parameters):
    """Isolate all process-global RNG state used by the audited third-party stack."""
    import random

    import numpy as np

    with _MILO_RNG_LOCK:
        python_rng_state = random.getstate()
        numpy_rng_state = np.random.get_state()
        try:
            return _standalone_milo_differential_abundance_impl(adata, **parameters)
        finally:
            random.setstate(python_rng_state)
            np.random.set_state(numpy_rng_state)


def run_milo_differential_abundance(
    adata: Any,
    *,
    sample_key: str,
    condition_key: str,
    reference_condition: str,
    comparison_condition: str,
    technical_batch_key: str = "",
    categorical_covariate_keys_json: str = "[]",
    continuous_covariate_keys_json: str = "[]",
    annotation_key: str = "cell_type",
    annotation_status: str = "provisional",
    representation_key: str = "X_pca",
    n_neighbors: int = 30,
    neighborhood_proportion: float = 0.1,
    mixed_annotation_threshold: float = 0.6,
    spatial_fdr_threshold: float = 0.1,
    min_abs_log2_fold_change: float = 0.0,
    random_seed: int = 123,
    openbio_version: str = "unknown",
) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    return _standalone_milo_differential_abundance(
        adata,
        sample_key=sample_key,
        condition_key=condition_key,
        reference_condition=reference_condition,
        comparison_condition=comparison_condition,
        technical_batch_key=technical_batch_key,
        categorical_covariate_keys_json=categorical_covariate_keys_json,
        continuous_covariate_keys_json=continuous_covariate_keys_json,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        representation_key=representation_key,
        n_neighbors=n_neighbors,
        neighborhood_proportion=neighborhood_proportion,
        mixed_annotation_threshold=mixed_annotation_threshold,
        spatial_fdr_threshold=spatial_fdr_threshold,
        min_abs_log2_fold_change=min_abs_log2_fold_change,
        random_seed=random_seed,
        openbio_version=openbio_version,
    )


def milo_differential_abundance_code(
    *,
    sample_key: str,
    condition_key: str,
    reference_condition: str,
    comparison_condition: str,
    technical_batch_key: str,
    categorical_covariate_keys_json: str,
    continuous_covariate_keys_json: str,
    annotation_key: str,
    annotation_status: str,
    representation_key: str,
    n_neighbors: int,
    neighborhood_proportion: float,
    mixed_annotation_threshold: float,
    spatial_fdr_threshold: float,
    min_abs_log2_fold_change: float,
    random_seed: int,
    openbio_version: str,
) -> str:
    implementation_source = textwrap.dedent(inspect.getsource(_standalone_milo_differential_abundance_impl)).strip()
    wrapper_source = textwrap.dedent(inspect.getsource(_standalone_milo_differential_abundance)).strip()
    return f"""from __future__ import annotations

import threading

_MILO_RNG_LOCK = threading.RLock()


{implementation_source}


{wrapper_source}


def run_milo_differential_abundance(adata):
    return _standalone_milo_differential_abundance(
        adata,
        sample_key={sample_key!r},
        condition_key={condition_key!r},
        reference_condition={reference_condition!r},
        comparison_condition={comparison_condition!r},
        technical_batch_key={technical_batch_key!r},
        categorical_covariate_keys_json={categorical_covariate_keys_json!r},
        continuous_covariate_keys_json={continuous_covariate_keys_json!r},
        annotation_key={annotation_key!r},
        annotation_status={annotation_status!r},
        representation_key={representation_key!r},
        n_neighbors={n_neighbors!r},
        neighborhood_proportion={neighborhood_proportion!r},
        mixed_annotation_threshold={mixed_annotation_threshold!r},
        spatial_fdr_threshold={spatial_fdr_threshold!r},
        min_abs_log2_fold_change={min_abs_log2_fold_change!r},
        random_seed={random_seed!r},
        openbio_version={openbio_version!r},
    )
"""


__all__ = [
    "MILO_TABLE_COLUMNS",
    "milo_differential_abundance_code",
    "run_milo_differential_abundance",
]
