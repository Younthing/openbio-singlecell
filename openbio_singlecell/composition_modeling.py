from __future__ import annotations

import inspect
import textwrap
import threading
from typing import Any

SCCODA_TABLE_COLUMNS = [
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
]

TASCCODA_TABLE_COLUMNS = [
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
]

_COMPOSITION_RNG_LOCK = threading.RLock()


class _PinnedPertpyBackend:
    """Narrow, fail-closed adapter for the audited Pertpy 1.3.0 composition API."""

    def __init__(self, method: str, openbio_version: str):
        import importlib
        import inspect as runtime_inspect
        from importlib import metadata as importlib_metadata

        required_modules = [
            "pertpy",
            "numpy",
            "pandas",
            "anndata",
            "jax",
            "jaxlib",
            "numpyro",
            "arviz",
            "mudata",
            "patsy",
        ]
        if method == "tasccoda":
            required_modules.extend(["ete4", "toytree"])
        try:
            pertpy_version = importlib_metadata.version("pertpy")
        except importlib_metadata.PackageNotFoundError as error:
            raise RuntimeError(f"{method} could not determine the installed Pertpy version.") from error
        if pertpy_version != "1.3.0":
            raise RuntimeError(
                f"{method} requires exact Pertpy 1.3.0; installed version is {pertpy_version!r}. "
                "Use an isolated environment with the audited version."
            )
        modules = {}
        for module_name in required_modules:
            try:
                modules[module_name] = importlib.import_module(module_name)
            except (ImportError, OSError) as error:
                extra = "pertpy[tcoda]" if method == "tasccoda" else "pertpy with its scCODA dependencies"
                raise RuntimeError(
                    f"{method} requires Pertpy 1.3.0 and the audited optional stack ({extra}); "
                    f"failed to import {module_name!r}."
                ) from error
        distribution_names = {
            "pertpy": "pertpy",
            "numpy": "numpy",
            "pandas": "pandas",
            "anndata": "anndata",
            "jax": "jax",
            "jaxlib": "jaxlib",
            "numpyro": "numpyro",
            "arviz": "arviz",
            "mudata": "mudata",
            "patsy": "patsy",
            "ete4": "ete4",
            "toytree": "toytree",
        }
        versions = {"openbio-singlecell": openbio_version}
        for key, distribution in distribution_names.items():
            if key not in required_modules:
                continue
            try:
                value = importlib_metadata.version(distribution)
            except importlib_metadata.PackageNotFoundError as error:
                raise RuntimeError(f"{method} could not determine the installed {distribution} version.") from error
            if not isinstance(value, str) or not value.strip():
                raise RuntimeError(f"{method} received an empty version for {distribution}.")
            versions[key] = value.strip()
        self.method = method
        self.modules = modules
        self.pertpy = modules["pertpy"]
        self.az = modules["arviz"]
        self.np = importlib.import_module("numpy")
        self.pd = importlib.import_module("pandas")
        self.inspect = runtime_inspect
        self.software_versions = versions

    def _require_parameters(self, callable_object: Any, required: list[str], label: str) -> None:
        try:
            parameters = self.inspect.signature(callable_object).parameters
        except (TypeError, ValueError) as error:
            raise RuntimeError(f"{self.method} cannot inspect the audited Pertpy {label} interface.") from error
        missing = [name for name in required if name not in parameters]
        if missing:
            raise RuntimeError(
                f"{self.method} requires the Pertpy 1.3.0 {label} interface; missing parameters: {missing}."
            )

    def _prepare_model(self, context: dict[str, Any], *, hierarchical: bool):
        model = self.pertpy.tl.Tasccoda() if hierarchical else self.pertpy.tl.Sccoda()
        load_required = [
            "adata",
            "type",
            "cell_type_identifier",
            "sample_identifier",
            "covariate_obs",
        ]
        if hierarchical:
            load_required.extend(["levels_orig", "add_level_name", "key_added"])
        else:
            load_required.append("generate_sample_level")
        self._require_parameters(model.load, load_required, f"{type(model).__name__}.load")
        prepare_required = [
            "data",
            "formula",
            "reference_cell_type",
            "automatic_reference_absence_threshold",
            "modality_key",
        ]
        if hierarchical:
            prepare_required.extend(["tree_key", "pen_args"])
        self._require_parameters(model.prepare, prepare_required, f"{type(model).__name__}.prepare")
        self._require_parameters(
            model.run_nuts,
            ["data", "modality_key", "rng_key", "num_samples", "num_warmup", "copy"],
            f"{type(model).__name__}.run_nuts",
        )
        summary_parameters = ["sample_adata"]
        if not hierarchical:
            summary_parameters.append("est" + "_fdr")
        self._require_parameters(
            model.summary_prepare,
            summary_parameters,
            f"{type(model).__name__}.summary_prepare",
        )
        self._require_parameters(
            model.make_arviz,
            ["data", "modality_key", "rng_key", "num_prior_samples", "use_posterior_predictive"],
            f"{type(model).__name__}.make_arviz",
        )
        load_kwargs = {
            "adata": context["working_adata"],
            "type": "cell_level",
            "cell_type_identifier": context["internal_annotation_key"],
            "sample_identifier": context["internal_sample_key"],
            "covariate_obs": context["design"]["model_columns"],
        }
        tree_key = None
        if hierarchical:
            tree_key = "__openbio_tasccoda_tree"
            suffix = 1
            while tree_key in context["working_adata"].uns:
                tree_key = f"__openbio_tasccoda_tree_{suffix}"
                suffix += 1
            load_kwargs.update(
                {
                    "levels_orig": [*context["hierarchy_keys"], context["internal_annotation_key"]],
                    "add_level_name": True,
                    "key_added": tree_key,
                }
            )
        else:
            load_kwargs["generate_sample_level"] = True
        mdata = model.load(**load_kwargs)
        if mdata is None or "coda" not in mdata.mod:
            raise RuntimeError(f"{context['operation']} Pertpy.load did not return a 'coda' modality.")
        sample_adata = mdata["coda"]
        observed_samples = [str(value) for value in sample_adata.obs_names.tolist()]
        observed_cell_types = [str(value) for value in sample_adata.var_names.tolist()]
        if set(observed_samples) != set(context["sample_internal_order"]):
            raise RuntimeError(f"{context['operation']} Pertpy.load changed the biological Sample identities.")
        if set(observed_cell_types) != set(context["cell_types"]):
            raise RuntimeError(f"{context['operation']} Pertpy.load changed the modeled cell-type identities.")
        sample_adata = sample_adata[context["sample_internal_order"], context["cell_types"]].copy()
        mdata.mod["coda"] = sample_adata
        backend_counts = self.np.asarray(sample_adata.X)
        if backend_counts.shape != context["count_matrix"].shape or not self.np.array_equal(
            backend_counts,
            context["count_matrix"],
        ):
            raise RuntimeError(f"{context['operation']} Pertpy.load produced a different Sample-by-cell-type table.")
        prepare_kwargs = {
            "data": mdata,
            "modality_key": "coda",
            "formula": context["design"]["formula"],
            "reference_cell_type": context["resolved_reference"],
            "automatic_reference_absence_threshold": 0.05,
        }
        if hierarchical:
            prepare_kwargs.update(
                {
                    "tree_key": tree_key,
                    "pen_args": {
                        "lambda_0": 50.0,
                        "lambda_1": 5.0,
                        "theta": 0.5,
                        "phi": context["aggregation_bias"],
                    },
                }
            )
        prepared = model.prepare(**prepare_kwargs)
        if prepared is not None:
            mdata = prepared
        if "coda" not in mdata.mod:
            raise RuntimeError(f"{context['operation']} Pertpy.prepare removed the 'coda' modality.")
        sample_adata = mdata["coda"]
        prepared_counts = self.np.asarray(sample_adata.X, dtype=float)
        expected_prepared = context["count_matrix"].astype(float)
        expected_prepared[expected_prepared == 0] = 0.5
        if not self.np.array_equal(prepared_counts, expected_prepared):
            raise RuntimeError(
                f"{context['operation']} Pertpy.prepare did not apply only the audited zero pseudocount 0.5."
            )
        model.run_nuts(
            mdata,
            modality_key="coda",
            rng_key=context["random_seed"],
            num_samples=context["num_samples"],
            num_warmup=context["num_warmup"],
            copy=False,
        )
        params = sample_adata.uns.get("scCODA_params")
        if not isinstance(params, dict):
            raise RuntimeError(f"{context['operation']} Pertpy omitted scCODA_params after fitting.")
        mcmc = params.get("mcmc")
        if not isinstance(mcmc, dict) or not isinstance(mcmc.get("samples"), dict):
            raise RuntimeError(f"{context['operation']} Pertpy omitted stored NUTS samples.")
        arviz_data = model.make_arviz(
            mdata,
            modality_key="coda",
            rng_key=context["random_seed"],
            num_prior_samples=0,
            use_posterior_predictive=False,
        )
        return model, mdata, sample_adata, params, mcmc, arviz_data, backend_counts

    def _diagnostics(self, mcmc: dict[str, Any], arviz_data: Any, variable: str, focal_index: int):
        try:
            posterior = arviz_data["posterior"].to_dataset()
            focal = posterior[variable].isel(covariate=focal_index).to_dataset(name=variable)
            ess_bulk = self.az.ess(focal, method="bulk").to_array().values
            ess_tail = self.az.ess(focal, method="tail").to_array().values
            mcse_mean = self.az.mcse(focal, method="mean").to_array().values
        except Exception as error:
            raise RuntimeError(f"{self.method} could not calculate audited within-chain ArviZ diagnostics.") from error
        extra = mcmc.get("extra_fields")
        if not isinstance(extra, dict):
            raise RuntimeError(f"{self.method} Pertpy omitted stored NUTS extra fields.")
        required_extra = {
            "potential_energy": "potential_energy",
            "num_steps": "num_steps",
            "step_size": "adapt_state_step_size",
        }
        diagnostics = {
            "acceptance_rate": float(self.np.asarray(mcmc.get("acceptance_rate"))),
            "ess_bulk": self.np.asarray(ess_bulk, dtype=float).reshape(-1),
            "ess_tail": self.np.asarray(ess_tail, dtype=float).reshape(-1),
            "mcse_mean": self.np.asarray(mcse_mean, dtype=float).reshape(-1),
            "rhat_available": False,
            "divergences_available": False,
        }
        for output_key, stored_key in required_extra.items():
            if stored_key not in extra:
                raise RuntimeError(f"{self.method} Pertpy omitted NUTS extra field {stored_key!r}.")
            diagnostics[output_key] = self.np.asarray(extra[stored_key], dtype=float).reshape(-1)
        return diagnostics

    def _runtime_details(self) -> dict[str, Any]:
        jax = self.modules["jax"]
        try:
            x64_enabled = bool(jax.config.read("jax_enable_x64"))
            devices = [
                {
                    "platform": str(device.platform),
                    "device_kind": str(device.device_kind),
                }
                for device in jax.devices()
            ]
        except Exception as error:
            raise RuntimeError(f"{self.method} could not inspect the audited JAX runtime.") from error
        if not x64_enabled:
            raise RuntimeError(f"{self.method} requires Pertpy's audited JAX 64-bit mode.")
        if not devices:
            raise RuntimeError(f"{self.method} found no JAX execution device.")
        return {"jax_enable_x64": True, "jax_devices": devices}

    def fit_sccoda(self, context: dict[str, Any]) -> dict[str, Any]:
        model, _mdata, sample_adata, params, mcmc, arviz_data, backend_counts = self._prepare_model(
            context,
            hierarchical=False,
        )
        summaries = model.summary_prepare(sample_adata, est_fdr=context["flat_expected_fdr"])
        if not isinstance(summaries, tuple) or len(summaries) != 2:
            raise RuntimeError("scCODA summary_prepare did not return the audited intercept/effect pair.")
        covariates = list(params.get("covariate_names", []))
        if covariates != context["design"]["model_columns"]:
            raise RuntimeError("scCODA Pertpy design columns differ from the audited exact contrast.")
        focal_index = covariates.index(context["design"]["focal_column"])
        samples = mcmc["samples"]
        beta = self.np.asarray(samples.get("beta"), dtype=float)
        alpha = self.np.asarray(samples.get("alpha"), dtype=float)
        expected_beta_shape = (context["num_samples"], len(covariates), len(context["cell_types"]))
        expected_alpha_shape = (context["num_samples"], len(context["cell_types"]))
        if beta.shape != expected_beta_shape or alpha.shape != expected_alpha_shape:
            raise RuntimeError(
                f"scCODA stored posterior shapes are beta={beta.shape}, alpha={alpha.shape}; "
                f"expected {expected_beta_shape} and {expected_alpha_shape}."
            )
        return {
            "cell_types": list(context["cell_types"]),
            "resolved_reference": params.get("reference_cell_type"),
            "posterior_intercept_samples": alpha,
            "posterior_focal_effect_samples": beta[:, focal_index, :],
            "backend_count_matrix": backend_counts,
            "backend_sample_ids": list(context["sample_internal_order"]),
            "upstream_contract": {
                "model_type": params.get("model_type"),
                "select_type": params.get("select_type"),
                "covariate_names": covariates,
                "reference_index": int(params.get("reference_index")),
                "automatic_reference_absence_threshold": float(params.get("automatic_reference_absence_threshold")),
                "pseudocount": 0.5,
                "algorithm": mcmc.get("algorithm"),
                "num_samples": int(mcmc.get("num_samples")),
                "num_warmup": int(mcmc.get("num_warmup")),
                "chain_count": int(mcmc.get("num_chains")),
            },
            "diagnostics": self._diagnostics(mcmc, arviz_data, "beta", focal_index),
            "backend_runtime": self._runtime_details(),
            "software_versions": dict(self.software_versions),
        }

    def fit_tasccoda(self, context: dict[str, Any]) -> dict[str, Any]:
        model, _mdata, sample_adata, params, mcmc, arviz_data, backend_counts = self._prepare_model(
            context,
            hierarchical=True,
        )
        summaries = model.summary_prepare(sample_adata)
        if not isinstance(summaries, tuple) or len(summaries) != 3:
            raise RuntimeError(
                "tascCODA summary_prepare did not return the audited intercept/effect/node-effect triple."
            )
        covariates = list(params.get("covariate_names", []))
        if covariates != context["design"]["model_columns"]:
            raise RuntimeError("tascCODA Pertpy design columns differ from the audited exact contrast.")
        focal_index = covariates.index(context["design"]["focal_column"])
        backend_nodes = list(params.get("node_names", []))
        canonical_nodes = context["hierarchy"]["node_names"]
        if set(backend_nodes) != set(canonical_nodes) or len(backend_nodes) != len(canonical_nodes):
            raise RuntimeError("tascCODA Pertpy hierarchy nodes differ from the canonical hierarchy manifest.")
        node_permutation = [backend_nodes.index(node) for node in canonical_nodes]
        backend_leaves = [str(value) for value in sample_adata.var_names.tolist()]
        if set(backend_leaves) != set(context["cell_types"]):
            raise RuntimeError("tascCODA Pertpy leaf axis differs from the canonical annotation leaf set.")
        leaf_permutation = [backend_leaves.index(leaf) for leaf in context["cell_types"]]
        backend_ancestor = self.np.asarray(params.get("ancestor_matrix"), dtype=float)
        canonical_ancestor = backend_ancestor[self.np.ix_(leaf_permutation, node_permutation)]
        expected_ancestor = context["hierarchy"]["ancestor_matrix"]
        if canonical_ancestor.shape != expected_ancestor.shape or not self.np.array_equal(
            canonical_ancestor,
            expected_ancestor,
        ):
            raise RuntimeError("tascCODA Pertpy ancestor matrix differs from the canonical hierarchy contract.")
        backend_reference_nodes = list(params.get("reference_nodes", []))
        if set(backend_reference_nodes) != set(context["hierarchy"]["reference_nodes"]):
            raise RuntimeError("tascCODA Pertpy reference ancestry differs from the canonical reference path.")
        reference_indices = [canonical_nodes.index(node) for node in context["hierarchy"]["reference_nodes"]]
        backend_reference_indices = [backend_nodes.index(node) for node in backend_reference_nodes]
        pen = params.get("sslasso_pen_args")
        if not isinstance(pen, dict):
            raise RuntimeError("tascCODA Pertpy omitted the tree-adaptive penalty parameters.")
        backend_nonreference = [index for index in range(len(backend_nodes)) if index not in backend_reference_indices]
        canonical_nonreference = [index for index in range(len(canonical_nodes)) if index not in reference_indices]
        full_node_leaves = self.np.full(len(backend_nodes), self.np.nan)
        full_scaled = self.np.full(len(backend_nodes), self.np.nan)
        full_node_leaves[backend_nonreference] = self.np.asarray(pen.get("node_leaves"), dtype=float)
        full_scaled[backend_nonreference] = self.np.asarray(pen.get("lambda_1_scaled"), dtype=float)
        canonical_node_leaves = full_node_leaves[node_permutation][canonical_nonreference]
        canonical_scaled = full_scaled[node_permutation][canonical_nonreference]
        samples = mcmc["samples"]
        node_effect = self.np.asarray(samples.get("b_tilde"), dtype=float)
        alpha = self.np.asarray(samples.get("alpha"), dtype=float)
        theta = self.np.asarray(samples.get("theta"), dtype=float).reshape(-1)
        expected_node_shape = (context["num_samples"], len(covariates), len(backend_nodes))
        expected_alpha_shape = (context["num_samples"], len(context["cell_types"]))
        if node_effect.shape != expected_node_shape or alpha.shape != expected_alpha_shape:
            raise RuntimeError(
                f"tascCODA stored posterior shapes are b_tilde={node_effect.shape}, alpha={alpha.shape}; "
                f"expected {expected_node_shape} and {expected_alpha_shape}."
            )
        if theta.shape != (context["num_samples"],):
            raise RuntimeError(f"tascCODA stored theta shape is {theta.shape}; expected {(context['num_samples'],)}.")
        return {
            "cell_types": list(context["cell_types"]),
            "resolved_reference": params.get("reference_cell_type"),
            "posterior_intercept_samples": alpha[:, leaf_permutation],
            "posterior_focal_node_effect_samples": node_effect[:, focal_index, :][:, node_permutation],
            "posterior_theta_samples": theta,
            "backend_count_matrix": backend_counts,
            "backend_sample_ids": list(context["sample_internal_order"]),
            "upstream_contract": {
                "model_type": params.get("model_type"),
                "select_type": params.get("select_type"),
                "covariate_names": covariates,
                "reference_index": reference_indices,
                "reference_nodes": list(context["hierarchy"]["reference_nodes"]),
                "node_names": canonical_nodes,
                "ancestor_matrix": canonical_ancestor,
                "automatic_reference_absence_threshold": float(params.get("automatic_reference_absence_threshold")),
                "pseudocount": 0.5,
                "algorithm": mcmc.get("algorithm"),
                "num_samples": int(mcmc.get("num_samples")),
                "num_warmup": int(mcmc.get("num_warmup")),
                "chain_count": int(mcmc.get("num_chains")),
                "lambda_0": float(pen.get("lambda_0")),
                "lambda_1": float(pen.get("lambda_1")),
                "theta": float(pen.get("theta")),
                "phi": float(pen.get("phi")),
                "node_leaves_nonreference": canonical_node_leaves,
                "lambda_1_scaled_nonreference": canonical_scaled,
            },
            "diagnostics": self._diagnostics(mcmc, arviz_data, "b_tilde", focal_index),
            "backend_runtime": self._runtime_details(),
            "software_versions": dict(self.software_versions),
        }


def _standalone_composition_model_impl(
    adata,
    *,
    method,
    sample_key,
    annotation_key,
    annotation_status,
    condition_key,
    reference_condition,
    comparison_condition,
    adjustment_covariate_keys_json="[]",
    reference_cell_type="automatic",
    flat_expected_fdr=None,
    hierarchy_keys_json="[]",
    aggregation_bias=0.0,
    num_samples=10000,
    num_warmup=1000,
    random_seed=0,
    openbio_version="unknown",
    _backend=None,
):
    """Validate and fit one flat or hierarchical Sample-level composition contrast."""
    import hashlib
    import json
    import math
    import numbers

    import numpy as np
    import pandas as pd

    operation = "scCODA Differential Composition" if method == "sccoda" else "tascCODA Differential Composition"
    if method not in {"sccoda", "tasccoda"}:
        raise ValueError(f"Unknown composition method {method!r}.")

    def canonical_text(value, label, *, allow_empty=False):
        if not isinstance(value, str):
            raise TypeError(f"{operation} {label} must be a string.")
        if value != value.strip():
            raise ValueError(f"{operation} {label} cannot contain surrounding whitespace.")
        if not value and not allow_empty:
            raise ValueError(f"{operation} {label} cannot be empty.")
        return value

    def numeric_value(value, label, *, integer=False, minimum=None, maximum=None, strict_maximum=False):
        expected = numbers.Integral if integer else numbers.Real
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, expected):
            kind = "integer" if integer else "number"
            raise TypeError(f"{operation} {label} must be a {kind}.")
        converted = int(value) if integer else float(value)
        if not integer and not math.isfinite(converted):
            raise ValueError(f"{operation} {label} must be finite.")
        if minimum is not None and converted < minimum:
            raise ValueError(f"{operation} {label} must be at least {minimum}.")
        if maximum is not None:
            invalid = converted >= maximum if strict_maximum else converted > maximum
            if invalid:
                relation = "less than" if strict_maximum else "at most"
                raise ValueError(f"{operation} {label} must be {relation} {maximum}.")
        return converted

    def parse_key_list(value, label, *, required=False):
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
        if required and not parsed:
            raise ValueError(f"{operation} {label} must contain at least one key.")
        if len(parsed) > 32:
            raise ValueError(f"{operation} {label} supports at most 32 keys.")
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
            values.append(canonical_text(value, f"{label} row {position}"))
        return values

    def sample_identity(value, position):
        missing = pd.isna(value)
        if not isinstance(missing, (bool, np.bool_)):
            raise TypeError(f"{operation} sample identifier row {position} is nonscalar.")
        if bool(missing):
            raise ValueError(f"{operation} sample identifiers contain missing values.")
        if isinstance(value, str):
            display = canonical_text(value, f"sample identifier row {position}")
            return ("string", display), display
        if isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{operation} sample identifiers cannot be booleans.")
        if isinstance(value, numbers.Integral):
            return ("integer", int(value)), int(value)
        if isinstance(value, numbers.Real) and math.isfinite(float(value)):
            return ("number", float(value)), float(value)
        raise TypeError(f"{operation} sample identifier row {position} must be a finite scalar or string.")

    def json_scalar(value):
        return value.item() if isinstance(value, np.generic) else value

    sample_key = canonical_text(sample_key, "sample_key")
    annotation_key = canonical_text(annotation_key, "annotation_key")
    annotation_status = canonical_text(annotation_status, "annotation_status")
    condition_key = canonical_text(condition_key, "condition_key")
    reference_condition = canonical_text(reference_condition, "reference_condition")
    comparison_condition = canonical_text(comparison_condition, "comparison_condition")
    reference_cell_type = canonical_text(reference_cell_type, "reference_cell_type")
    openbio_version = canonical_text(openbio_version, "openbio_version")
    if annotation_status not in {"provisional", "curated"}:
        raise ValueError(f"{operation} annotation_status must be 'provisional' or 'curated'.")
    if reference_condition == comparison_condition:
        raise ValueError(f"{operation} reference and comparison Conditions must be different.")
    adjustment_keys = parse_key_list(adjustment_covariate_keys_json, "adjustment_covariate_keys_json")
    hierarchy_keys = parse_key_list(
        hierarchy_keys_json,
        "hierarchy_keys_json",
        required=method == "tasccoda",
    )
    if method == "sccoda" and hierarchy_keys:
        raise ValueError("scCODA does not accept hierarchy keys; use the tascCODA node.")
    role_keys = [sample_key, annotation_key, condition_key, *adjustment_keys, *hierarchy_keys]
    if len(role_keys) != len(set(role_keys)):
        raise ValueError(f"{operation} metadata roles must use distinct observation columns.")
    num_samples = numeric_value(num_samples, "num_samples", integer=True, minimum=1)
    num_warmup = numeric_value(num_warmup, "num_warmup", integer=True, minimum=1)
    random_seed = numeric_value(random_seed, "random_seed", integer=True, minimum=0, maximum=2**31 - 1)
    aggregation_bias = numeric_value(aggregation_bias, "aggregation_bias")
    if method == "sccoda":
        flat_expected_fdr = numeric_value(
            flat_expected_fdr,
            "estimated" + "_fdr",
            minimum=0.0,
            maximum=1.0,
        )
        if flat_expected_fdr in {0.0, 1.0}:
            raise ValueError(f"{operation} posterior expected-FDR target must be strictly between 0 and 1.")
    elif flat_expected_fdr is not None:
        raise ValueError("tascCODA does not use an estimated-FDR parameter.")

    if bool(getattr(adata, "isbacked", False)):
        raise ValueError(f"{operation} requires an in-memory AnnData; call to_memory() first.")
    if int(getattr(adata, "n_obs", 0)) < 1 or int(getattr(adata, "n_vars", 0)) < 1:
        raise ValueError(f"{operation} requires non-empty observation and feature axes.")
    if not bool(adata.obs_names.is_unique) or not bool(adata.var_names.is_unique):
        raise ValueError(f"{operation} requires unique observation and feature identifiers.")
    required_obs = list(dict.fromkeys(role_keys))
    missing_obs = [key for key in required_obs if key not in adata.obs]
    if missing_obs:
        raise ValueError(f"{operation} observation columns not found: {missing_obs}.")

    condition_values = string_values(adata.obs[condition_key], condition_key)
    annotation_values = string_values(adata.obs[annotation_key], annotation_key)
    hierarchy_values = {key: string_values(adata.obs[key], key) for key in hierarchy_keys}
    for key, values in hierarchy_values.items():
        series = adata.obs[key]
        if isinstance(series.dtype, pd.CategoricalDtype):
            declared = [canonical_text(value, f"{key} category") for value in series.cat.categories.tolist()]
            unused = [value for value in declared if value not in set(values)]
            if unused:
                raise ValueError(f"{operation} hierarchy key {key!r} has unused categories: {unused}.")
    sample_identities = []
    sample_displays = []
    for position, value in enumerate(adata.obs[sample_key].tolist()):
        identity, display = sample_identity(value, position)
        sample_identities.append(identity)
        sample_displays.append(display)
    sample_order = list(dict.fromkeys(sample_identities))
    sample_internal = {identity: f"sample_{position:06d}" for position, identity in enumerate(sample_order)}

    cell_metadata = pd.DataFrame(
        {
            "sample_internal": [sample_internal[value] for value in sample_identities],
            "sample": sample_displays,
            "condition": condition_values,
            "cell_type": annotation_values,
            **{key: values for key, values in hierarchy_values.items()},
        },
        index=adata.obs_names.copy(),
    )
    for key in adjustment_keys:
        cell_metadata[key] = adata.obs[key].to_numpy(copy=True)

    constant_keys = ["condition", *adjustment_keys]
    for key in constant_keys:
        per_sample = cell_metadata.groupby("sample_internal", sort=False, observed=True)[key].nunique(dropna=False)
        conflicting = per_sample[per_sample != 1]
        if not conflicting.empty:
            raise ValueError(
                f"{operation} requires {key!r} to be unique within each biological Sample; "
                f"conflicting internal Samples: {conflicting.index.tolist()}."
            )
    sample_metadata = cell_metadata.drop_duplicates("sample_internal", keep="first").copy()
    observed_conditions = set(sample_metadata["condition"].tolist())
    declared_conditions = {reference_condition, comparison_condition}
    if observed_conditions != declared_conditions:
        raise ValueError(
            f"{operation} requires exactly the declared two Condition levels {sorted(declared_conditions)!r}; "
            f"observed {sorted(observed_conditions)!r}. Subset upstream or correct the explicit contrast."
        )
    samples_per_condition = sample_metadata.groupby("condition", observed=True)["sample_internal"].nunique()
    insufficient = {condition: int(samples_per_condition.get(condition, 0)) for condition in declared_conditions}
    if any(count < 2 for count in insufficient.values()):
        raise ValueError(
            f"{operation} requires at least two biological Samples per Condition; observed {insufficient}."
        )

    def ordered_levels(series, values, label):
        observed = set(values)
        if isinstance(series.dtype, pd.CategoricalDtype):
            declared = [canonical_text(value, f"{label} category") for value in series.cat.categories.tolist()]
            unused = [value for value in declared if value not in observed]
            if unused:
                raise ValueError(f"{operation} {label} has unused declared categories: {unused}.")
            if len(declared) != len(observed):
                raise RuntimeError(f"{operation} could not preserve {label} category identities.")
            return declared
        return sorted(observed)

    def safe_column(base, occupied):
        candidate = base
        suffix = 1
        while candidate in occupied:
            candidate = f"{base}_{suffix}"
            suffix += 1
        occupied.add(candidate)
        return candidate

    def finite_array(value, label, *, shape=None):
        array = np.asarray(value, dtype=float)
        if shape is not None and tuple(array.shape) != tuple(shape):
            raise RuntimeError(f"{operation} backend {label} has shape {array.shape}; expected {shape}.")
        if array.size == 0 or not bool(np.isfinite(array).all()):
            raise RuntimeError(f"{operation} backend {label} must be non-empty and finite.")
        return array

    def numeric_summary(value, label, *, positive=False, nonnegative=False):
        array = finite_array(value, label).reshape(-1)
        if positive and bool((array <= 0).any()):
            raise RuntimeError(f"{operation} backend {label} must be positive.")
        if nonnegative and bool((array < 0).any()):
            raise RuntimeError(f"{operation} backend {label} must be nonnegative.")
        return {
            "min": float(array.min()),
            "median": float(np.median(array)),
            "max": float(array.max()),
        }

    cell_types = ordered_levels(adata.obs[annotation_key], annotation_values, annotation_key)
    if len(cell_types) < 2:
        raise ValueError(f"{operation} requires at least two modeled cell types.")
    sample_internal_order = [sample_internal[identity] for identity in sample_order]
    count_frame = (
        cell_metadata.groupby(["sample_internal", "cell_type"], sort=False, observed=True)
        .size()
        .rename("count")
        .reindex(
            pd.MultiIndex.from_product(
                [sample_internal_order, cell_types],
                names=["sample_internal", "cell_type"],
            ),
            fill_value=0,
        )
        .reset_index()
    )
    count_matrix = count_frame.pivot(index="sample_internal", columns="cell_type", values="count").reindex(
        index=sample_internal_order,
        columns=cell_types,
    )
    count_matrix = count_matrix.to_numpy(dtype=np.int64)
    if bool((count_matrix < 0).any()) or not bool(np.equal(count_matrix, np.floor(count_matrix)).all()):
        raise RuntimeError(f"{operation} constructed an invalid count matrix.")
    sample_totals = count_matrix.sum(axis=1)
    if bool((sample_totals <= 0).any()):
        raise ValueError(f"{operation} requires a positive cell total for every Sample.")
    count_digest = hashlib.sha256()
    count_digest.update(b"openbio-singlecell/composition-counts/v1\0")
    count_digest.update(json.dumps(sample_internal_order, separators=(",", ":")).encode("ascii"))
    count_digest.update(json.dumps(cell_types, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    count_digest.update(np.ascontiguousarray(count_matrix, dtype="<i8").tobytes())
    count_fingerprint = count_digest.hexdigest()

    occupied = set(adata.obs.columns)
    internal_sample_key = safe_column("__openbio_sample", occupied)
    internal_annotation_key = safe_column("__openbio_cell_type", occupied)
    focal_column = safe_column("__openbio_condition_comparison", occupied)
    sample_condition = dict(zip(sample_metadata["sample_internal"], sample_metadata["condition"], strict=True))
    focal_by_sample = {
        sample: 1.0 if condition == comparison_condition else 0.0 for sample, condition in sample_condition.items()
    }
    design_columns = [focal_column]
    design_values = [np.array([focal_by_sample[sample] for sample in sample_internal_order], dtype=float)]
    adjustment_report = []
    adjustment_cell_columns = {}
    for position, key in enumerate(adjustment_keys):
        original = adata.obs[key]
        sample_values = sample_metadata.set_index("sample_internal")[key].reindex(sample_internal_order)
        if bool(sample_values.isna().any()):
            raise ValueError(f"{operation} adjustment covariate {key!r} contains missing Sample values.")
        if pd.api.types.is_bool_dtype(original.dtype):
            raise TypeError(f"{operation} adjustment covariate {key!r} cannot be boolean.")
        if pd.api.types.is_numeric_dtype(original.dtype):
            values = sample_values.to_numpy(dtype=float)
            if not bool(np.isfinite(values).all()):
                raise ValueError(f"{operation} numeric adjustment covariate {key!r} must be finite.")
            if bool(np.allclose(values, values[0])):
                raise ValueError(f"{operation} numeric adjustment covariate {key!r} is constant.")
            unique_values = np.unique(values)
            paired_numeric = len(unique_values) * 2 == len(values) and all(
                {
                    sample_condition[sample_internal_order[index]]
                    for index, observed in enumerate(values)
                    if observed == value
                }
                == declared_conditions
                for value in unique_values
            )
            identifier_name = any(
                token in key.casefold() for token in ["subject", "patient", "donor", "pair", "individual"]
            )
            if paired_numeric and (len(unique_values) >= 3 or identifier_name):
                raise ValueError(
                    f"{operation} adjustment covariate {key!r} encodes a paired/repeated-Sample pattern. "
                    "This one-row-per-biological-Sample model has no random-effect or repeated-measures contract."
                )
            center = float(values.mean())
            column = safe_column(f"__openbio_adjustment_{position}", occupied)
            encoded = values - center
            design_columns.append(column)
            design_values.append(encoded)
            value_by_sample = dict(zip(sample_internal_order, encoded, strict=True))
            adjustment_cell_columns[column] = (
                cell_metadata["sample_internal"].map(value_by_sample).to_numpy(dtype=float)
            )
            adjustment_report.append(
                {
                    "key": key,
                    "kind": "continuous_mean_centered",
                    "center": center,
                    "model_columns": [column],
                }
            )
            continue
        if not isinstance(original.dtype, pd.CategoricalDtype) or not bool(original.cat.ordered):
            raise TypeError(
                f"{operation} categorical adjustment covariate {key!r} must be an ordered pandas Categorical; "
                "the first declared category is the reference."
            )
        categories = [canonical_text(value, f"{key} category") for value in original.cat.categories.tolist()]
        observed_adjustment = string_values(original, key)
        unused = [value for value in categories if value not in set(observed_adjustment)]
        if unused:
            raise ValueError(f"{operation} adjustment covariate {key!r} has unused categories: {unused}.")
        sample_strings = [canonical_text(value, f"{key} Sample value") for value in sample_values.tolist()]
        if len(set(sample_strings)) < 2:
            raise ValueError(f"{operation} categorical adjustment covariate {key!r} is constant.")
        level_conditions = {}
        for level in categories:
            indices = [index for index, value in enumerate(sample_strings) if value == level]
            level_conditions[level] = {sample_condition[sample_internal_order[index]] for index in indices}
        paired_categorical = len(categories) * 2 == len(sample_internal_order) and all(
            conditions == declared_conditions for conditions in level_conditions.values()
        )
        identifier_name = any(
            token in key.casefold() for token in ["subject", "patient", "donor", "pair", "individual"]
        )
        if paired_categorical and (len(categories) >= 3 or identifier_name):
            raise ValueError(
                f"{operation} adjustment covariate {key!r} encodes a paired/repeated-Sample pattern. "
                "This one-row-per-biological-Sample model has no random-effect or repeated-measures contract."
            )
        created = []
        for level_position, level in enumerate(categories[1:], start=1):
            column = safe_column(f"__openbio_adjustment_{position}_{level_position}", occupied)
            encoded = np.array([float(value == level) for value in sample_strings], dtype=float)
            design_columns.append(column)
            design_values.append(encoded)
            value_by_sample = dict(zip(sample_internal_order, encoded, strict=True))
            adjustment_cell_columns[column] = (
                cell_metadata["sample_internal"].map(value_by_sample).to_numpy(dtype=float)
            )
            created.append(column)
        adjustment_report.append(
            {
                "key": key,
                "kind": "ordered_categorical_treatment",
                "reference_category": categories[0],
                "categories": categories,
                "model_columns": created,
            }
        )
    design_matrix = np.column_stack([np.ones(len(sample_internal_order), dtype=float), *design_values])
    design_rank = int(np.linalg.matrix_rank(design_matrix))
    if design_rank != design_matrix.shape[1]:
        raise ValueError(
            f"{operation} design is rank deficient or perfectly confounded "
            f"(rank {design_rank}, {design_matrix.shape[1]} columns)."
        )
    residual_df = int(design_matrix.shape[0] - design_rank)
    if residual_df < 1:
        raise ValueError(f"{operation} design has no residual degrees of freedom after Sample-level adjustment.")
    condition_number = float(np.linalg.cond(design_matrix))
    if not math.isfinite(condition_number):
        raise ValueError(f"{operation} design has a non-finite condition number.")
    design_formula = "1 + " + " + ".join(design_columns)

    relative_counts = count_matrix / sample_totals[:, None]
    zero_fractions = np.mean(count_matrix == 0, axis=0)
    means = np.mean(relative_counts, axis=0)
    dispersions = np.divide(
        np.var(relative_counts, axis=0),
        means,
        out=np.full(len(cell_types), np.inf, dtype=float),
        where=means > 0,
    )
    reference_candidates = [
        {
            "cell_type": cell_type,
            "zero_fraction": float(zero_fractions[index]),
            "dispersion": float(dispersions[index]),
            "eligible": bool(zero_fractions[index] < 0.05),
        }
        for index, cell_type in enumerate(cell_types)
    ]
    if reference_cell_type == "automatic":
        eligible = [index for index, value in enumerate(zero_fractions) if value < 0.05]
        if not eligible:
            raise ValueError(
                f"{operation} automatic reference found no cell type with zero fraction strictly below 0.05."
            )
        resolved_reference_index = min(eligible, key=lambda index: (dispersions[index], index))
        resolved_reference = cell_types[resolved_reference_index]
    else:
        if reference_cell_type not in cell_types:
            raise ValueError(
                f"{operation} reference_cell_type {reference_cell_type!r} is not a modeled cell type: {cell_types}."
            )
        resolved_reference = reference_cell_type
        resolved_reference_index = cell_types.index(resolved_reference)
    if bool((count_matrix[:, resolved_reference_index] <= 0).all()):
        raise ValueError(f"{operation} reference cell type {resolved_reference!r} has no observed cells.")

    hierarchy = None
    if method == "tasccoda":
        paths_by_leaf = {}
        for row in cell_metadata[[*hierarchy_keys, "cell_type"]].itertuples(index=False, name=None):
            path = tuple(row[:-1])
            leaf = row[-1]
            previous = paths_by_leaf.setdefault(leaf, path)
            if previous != path:
                raise ValueError(
                    f"{operation} leaf {leaf!r} maps to more than one hierarchy path: {previous!r} and {path!r}."
                )
        if set(paths_by_leaf) != set(cell_types):
            raise ValueError(f"{operation} hierarchy leaf set does not equal the modeled annotation leaf set.")
        roots = {path[0] for path in paths_by_leaf.values()}
        if len(roots) != 1:
            raise ValueError(f"{operation} hierarchy must have one root; observed {sorted(roots)!r}.")
        for level_index in range(1, len(hierarchy_keys)):
            parents = {}
            for path in paths_by_leaf.values():
                child = path[level_index]
                parent = path[level_index - 1]
                old_parent = parents.setdefault(child, parent)
                if old_parent != parent:
                    raise ValueError(
                        f"{operation} hierarchy child {child!r} has multiple parents: {old_parent!r}, {parent!r}."
                    )
        root_value = next(iter(roots))
        root_name = f"{hierarchy_keys[0]}_{root_value}"
        node_records = []
        seen_nodes = {root_name}
        for level_index, key in enumerate(hierarchy_keys[1:], start=1):
            for value in sorted({path[level_index] for path in paths_by_leaf.values()}):
                if level_index + 1 < len(hierarchy_keys):
                    immediate_children = {
                        path[level_index + 1] for path in paths_by_leaf.values() if path[level_index] == value
                    }
                else:
                    immediate_children = {leaf for leaf, path in paths_by_leaf.items() if path[level_index] == value}
                if len(immediate_children) < 2:
                    # Pertpy 1.3.0 collapses unary hierarchy nodes before constructing the ancestor matrix.
                    continue
                name = f"{key}_{value}"
                if name in seen_nodes:
                    raise ValueError(f"{operation} hierarchy produces duplicate prefixed node name {name!r}.")
                seen_nodes.add(name)
                descendants = sorted(leaf for leaf, path in paths_by_leaf.items() if path[level_index] == value)
                node_records.append(
                    {
                        "name": name,
                        "level": key,
                        "descendant_leaves": descendants,
                    }
                )
        for leaf in cell_types:
            if leaf in seen_nodes:
                raise ValueError(f"{operation} hierarchy leaf {leaf!r} collides with an internal node name.")
            seen_nodes.add(leaf)
            node_records.append({"name": leaf, "level": annotation_key, "descendant_leaves": [leaf]})
        node_names = [record["name"] for record in node_records]
        ancestor_matrix = np.array(
            [[float(leaf in record["descendant_leaves"]) for record in node_records] for leaf in cell_types],
            dtype=float,
        )
        reference_path_values = paths_by_leaf[resolved_reference]
        declared_reference_path = [root_name]
        declared_reference_path.extend(
            f"{hierarchy_keys[index]}_{reference_path_values[index]}" for index in range(1, len(hierarchy_keys))
        )
        declared_reference_path.append(resolved_reference)
        reference_path = [node for node in declared_reference_path if node in {root_name, *node_names}]
        reference_nodes = [node for node in reference_path if node in node_names]
        edges = []
        for leaf, path in sorted(paths_by_leaf.items()):
            names = [root_name]
            names.extend(f"{hierarchy_keys[index]}_{path[index]}" for index in range(1, len(hierarchy_keys)))
            names.append(leaf)
            names = [name for name in names if name in {root_name, *node_names}]
            edges.extend(zip(names[:-1], names[1:], strict=True))
        edges = sorted(set(edges))
        hierarchy_digest = hashlib.sha256()
        hierarchy_digest.update(b"openbio-singlecell/composition-hierarchy/v1\0")
        hierarchy_digest.update(json.dumps(edges, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
        hierarchy = {
            "ancestor_keys": hierarchy_keys,
            "levels_passed_to_pertpy": [*hierarchy_keys, annotation_key],
            "root": root_name,
            "edges": [list(edge) for edge in edges],
            "node_records": node_records,
            "node_names": node_names,
            "ancestor_matrix": ancestor_matrix,
            "declared_reference_path": declared_reference_path,
            "reference_path": reference_path,
            "reference_nodes": reference_nodes,
            "fingerprint": hierarchy_digest.hexdigest(),
        }

    working_adata = adata.copy()
    working_adata.obs[internal_sample_key] = [sample_internal[value] for value in sample_identities]
    working_adata.obs[internal_annotation_key] = annotation_values
    working_adata.obs[focal_column] = cell_metadata["sample_internal"].map(focal_by_sample).to_numpy(dtype=float)
    for column, values in adjustment_cell_columns.items():
        working_adata.obs[column] = values

    context = {
        "method": method,
        "operation": operation,
        "working_adata": working_adata,
        "cell_metadata": cell_metadata,
        "sample_metadata": sample_metadata,
        "sample_key": sample_key,
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "condition_key": condition_key,
        "reference_condition": reference_condition,
        "comparison_condition": comparison_condition,
        "adjustment_keys": adjustment_keys,
        "design": {
            "formula": design_formula,
            "focal_column": focal_column,
            "model_columns": design_columns,
            "matrix": design_matrix,
            "rank": design_rank,
            "residual_degrees_of_freedom": residual_df,
            "condition_number": condition_number,
            "adjustments": adjustment_report,
        },
        "hierarchy_keys": hierarchy_keys,
        "hierarchy": hierarchy,
        "reference_cell_type": reference_cell_type,
        "resolved_reference": resolved_reference,
        "resolved_reference_index": resolved_reference_index,
        "reference_candidates": reference_candidates,
        "cell_types": cell_types,
        "sample_internal_order": sample_internal_order,
        "sample_display_order": [
            json_scalar(sample_metadata.set_index("sample_internal").loc[sample, "sample"])
            for sample in sample_internal_order
        ],
        "count_matrix": count_matrix,
        "count_fingerprint": count_fingerprint,
        "sample_totals": sample_totals,
        "internal_sample_key": internal_sample_key,
        "internal_annotation_key": internal_annotation_key,
        "flat_expected_fdr": flat_expected_fdr,
        "aggregation_bias": aggregation_bias,
        "num_samples": num_samples,
        "num_warmup": num_warmup,
        "random_seed": random_seed,
        "openbio_version": openbio_version,
    }
    if _backend is None:
        _backend = _PinnedPertpyBackend(method, openbio_version)
    backend_method = getattr(_backend, f"fit_{method}", None)
    if not callable(backend_method):
        raise TypeError(f"{operation} backend must provide fit_{method}(context).")
    backend_result = backend_method(context)
    if not isinstance(backend_result, dict):
        raise TypeError(f"{operation} backend must return a dictionary.")
    if backend_result.get("cell_types") != cell_types:
        raise RuntimeError(f"{operation} backend changed or reordered the modeled cell-type axis.")
    if backend_result.get("resolved_reference") != resolved_reference:
        raise RuntimeError(f"{operation} backend resolved a different compositional reference.")
    backend_counts = np.asarray(backend_result.get("backend_count_matrix"))
    if backend_counts.shape != count_matrix.shape or not bool(np.array_equal(backend_counts, count_matrix)):
        raise RuntimeError(f"{operation} backend Sample-by-cell-type counts differ from the audited count table.")
    if backend_result.get("backend_sample_ids") != sample_internal_order:
        raise RuntimeError(f"{operation} backend changed or reordered the biological Sample axis.")
    contract = backend_result.get("upstream_contract")
    if not isinstance(contract, dict):
        raise RuntimeError(f"{operation} backend omitted its audited upstream contract.")
    expected_contract = {
        "model_type": "classic" if method == "sccoda" else "tree_agg",
        "select_type": "spikeslab" if method == "sccoda" else "sslasso",
        "covariate_names": design_columns,
        "algorithm": "NUTS",
        "num_samples": num_samples,
        "num_warmup": num_warmup,
        "chain_count": 1,
    }
    for key, expected in expected_contract.items():
        if contract.get(key) != expected:
            raise RuntimeError(f"{operation} backend contract {key!r} is {contract.get(key)!r}; expected {expected!r}.")
    if float(contract.get("automatic_reference_absence_threshold", math.nan)) != 0.05:
        raise RuntimeError(f"{operation} backend did not use the fixed 0.05 reference-candidate threshold.")
    if float(contract.get("pseudocount", math.nan)) != 0.5:
        raise RuntimeError(f"{operation} backend did not use Pertpy's fixed zero pseudocount 0.5.")
    if method == "sccoda" and int(contract.get("reference_index", -1)) != resolved_reference_index:
        raise RuntimeError(f"{operation} backend reference index differs from the audited cell-type axis.")

    software_versions = backend_result.get("software_versions")
    if not isinstance(software_versions, dict) or software_versions.get("pertpy") != "1.3.0":
        raise RuntimeError(f"{operation} requires exact Pertpy version 1.3.0.")
    required_software = {
        "openbio-singlecell",
        "pertpy",
        "numpy",
        "pandas",
        "anndata",
        "jax",
        "jaxlib",
        "numpyro",
        "arviz",
        "mudata",
        "patsy",
    }
    if method == "tasccoda":
        required_software.update({"ete4", "toytree"})
    missing_versions = sorted(required_software.difference(software_versions))
    if missing_versions or any(
        not isinstance(software_versions[key], str) or not software_versions[key].strip() for key in required_software
    ):
        raise RuntimeError(f"{operation} backend omitted software versions: {missing_versions}.")
    if software_versions["openbio-singlecell"] != openbio_version:
        raise RuntimeError(f"{operation} backend reported a different OpenBio version than the executing node.")
    backend_runtime = backend_result.get("backend_runtime")
    if not isinstance(backend_runtime, dict) or backend_runtime.get("jax_enable_x64") is not True:
        raise RuntimeError(f"{operation} backend did not confirm JAX 64-bit mode.")
    jax_devices = backend_runtime.get("jax_devices")
    if (
        not isinstance(jax_devices, list)
        or not jax_devices
        or any(
            not isinstance(item, dict)
            or not isinstance(item.get("platform"), str)
            or not item["platform"].strip()
            or not isinstance(item.get("device_kind"), str)
            or not item["device_kind"].strip()
            for item in jax_devices
        )
    ):
        raise RuntimeError(f"{operation} backend did not report its JAX execution devices.")

    diagnostics_raw = backend_result.get("diagnostics")
    if not isinstance(diagnostics_raw, dict):
        raise RuntimeError(f"{operation} backend omitted diagnostics.")
    acceptance_rate = float(diagnostics_raw.get("acceptance_rate", math.nan))
    if not math.isfinite(acceptance_rate) or not 0 <= acceptance_rate <= 1:
        raise RuntimeError(f"{operation} backend acceptance rate must be finite and in [0, 1].")
    if diagnostics_raw.get("rhat_available") is not False:
        raise RuntimeError(f"{operation} cannot claim R-hat availability for the audited one-chain fit.")
    if diagnostics_raw.get("divergences_available") is not False:
        raise RuntimeError(f"{operation} cannot claim divergence availability for Pertpy 1.3.0.")
    diagnostic_status = "limited_single_chain" if 0.6 <= acceptance_rate <= 0.95 else "warning"
    diagnostics = {
        "status": diagnostic_status,
        "chain_count": 1,
        "acceptance_rate": acceptance_rate,
        "acceptance_policy": {"minimum": 0.6, "maximum": 0.95, "inside_range": diagnostic_status != "warning"},
        "ess_bulk": numeric_summary(diagnostics_raw.get("ess_bulk"), "bulk ESS", positive=True),
        "ess_tail": numeric_summary(diagnostics_raw.get("ess_tail"), "tail ESS", positive=True),
        "mcse_mean": numeric_summary(diagnostics_raw.get("mcse_mean"), "mean MCSE", nonnegative=True),
        "potential_energy": numeric_summary(diagnostics_raw.get("potential_energy"), "potential energy"),
        "num_steps": numeric_summary(diagnostics_raw.get("num_steps"), "NUTS step count", positive=True),
        "step_size": numeric_summary(diagnostics_raw.get("step_size"), "NUTS step size", positive=True),
        "rhat_available": False,
        "rhat": None,
        "divergences_available": False,
        "divergences": None,
    }

    contrast = f"{comparison_condition} vs {reference_condition}"
    warnings = []
    if annotation_status == "provisional":
        warnings.append("Annotation status is provisional; treat the composition result as exploratory.")
    if reference_cell_type == "automatic":
        warnings.append(
            "The compositional reference was selected from these data and is not proven biologically invariant."
        )
    if min(insufficient.values()) < 3:
        warnings.append(
            "Only two biological Samples are available in at least one Condition; power and stability are limited."
        )
    if diagnostic_status == "warning":
        warnings.append("Pertpy's mean NUTS acceptance rate is outside the audited [0.6, 0.95] range.")
    if method == "tasccoda":
        warning_leaf_counts = hierarchy["ancestor_matrix"].sum(axis=0).astype(float)
        warning_scale = 1.0 / (
            1.0
            + np.exp(
                np.clip(
                    -aggregation_bias * (warning_leaf_counts / len(cell_types) - 0.5),
                    -700.0,
                    700.0,
                )
            )
        )
        if float(np.min(warning_scale)) < 0.01 or float(np.max(warning_scale)) > 0.99:
            warnings.append(
                "The declared aggregation bias nearly saturates at least one hierarchy penalty scale; "
                "assess sensitivity to alternative finite values."
            )

    references = [
        {
            "citation": (
                "Büttner M, Ostner J, Müller CL, Theis FJ, Schubert B. scCODA is a Bayesian model for "
                "compositional single-cell data analysis. Nature Communications. 2021;12:6876."
            ),
            "doi": "10.1038/s41467-021-27150-6",
            "url": "https://doi.org/10.1038/s41467-021-27150-6",
            "kind": "method",
        },
        {
            "citation": f"Pertpy developers. Official {'scCODA' if method == 'sccoda' else 'tascCODA'} API, version 1.3.0.",
            "url": (
                "https://pertpy.readthedocs.io/en/latest/api/tools/pertpy.tools.Sccoda.html"
                if method == "sccoda"
                else "https://pertpy.readthedocs.io/en/latest/api/tools/pertpy.tools.Tasccoda.html"
            ),
            "kind": "software_documentation",
        },
        {
            "citation": "Heumos L, et al. Pertpy: an end-to-end framework for perturbation analysis. Nature Methods.",
            "doi": "10.1038/s41592-025-02909-7",
            "url": "https://doi.org/10.1038/s41592-025-02909-7",
            "kind": "software",
        },
        {
            "citation": "Hoffman MD, Gelman A. The No-U-Turn Sampler. Journal of Machine Learning Research. 2014.",
            "url": "https://jmlr.org/papers/v15/hoffman14a.html",
            "kind": "sampler",
        },
        {
            "citation": "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. bioRxiv.",
            "doi": "10.1101/2021.12.16.473007",
            "url": "https://doi.org/10.1101/2021.12.16.473007",
            "kind": "software",
        },
        {
            "citation": "Bradbury J, et al. JAX: composable transformations of Python+NumPy programs. 2018.",
            "url": "https://github.com/jax-ml/jax",
            "kind": "software",
        },
        {
            "citation": "Phan D, Pradhan N, Jankowiak M. Composable effects for flexible and accelerated probabilistic programming in NumPyro. 2019.",
            "doi": "10.48550/arXiv.1912.11554",
            "url": "https://doi.org/10.48550/arXiv.1912.11554",
            "kind": "software",
        },
        {
            "citation": "Kumar R, Carroll C, Hartikainen A, Martin O. ArviZ: exploratory analysis of Bayesian models. JOSS. 2019.",
            "doi": "10.21105/joss.01143",
            "url": "https://doi.org/10.21105/joss.01143",
            "kind": "software",
        },
    ]

    posterior_intercepts = finite_array(
        backend_result.get("posterior_intercept_samples"),
        "posterior intercept samples",
        shape=(num_samples, len(cell_types)),
    )
    mean_total = float(sample_totals.mean())
    intercept = np.mean(posterior_intercepts, axis=0)
    expected_reference = np.exp(intercept - np.max(intercept))
    expected_reference = expected_reference / expected_reference.sum() * mean_total

    if method == "sccoda":
        focal_samples = finite_array(
            backend_result.get("posterior_focal_effect_samples"),
            "posterior focal-effect samples",
            shape=(num_samples, len(cell_types)),
        )
        if not bool(np.allclose(focal_samples[:, resolved_reference_index], 0.0, atol=1e-12, rtol=0.0)):
            raise RuntimeError(f"{operation} backend violated the reference-cell zero-effect constraint.")
        included = np.abs(focal_samples) > 1e-3
        inclusion_probability = included.mean(axis=0)
        conditional_mean = np.array(
            [
                float(focal_samples[included[:, index], index].mean()) if included[:, index].any() else 0.0
                for index in range(len(cell_types))
            ]
        )
        candidate_probabilities = np.sort(np.unique(inclusion_probability[inclusion_probability > 0]))
        inclusion_threshold = 1.0
        threshold_expected_fdr = 0.0
        threshold_fallback = True
        for candidate in candidate_probabilities:
            selected = inclusion_probability >= candidate
            candidate_fdr = float(np.mean(1.0 - inclusion_probability[selected]))
            if candidate_fdr < flat_expected_fdr:
                inclusion_threshold = math.floor(float(candidate) * 1000.0) / 1000.0
                threshold_expected_fdr = candidate_fdr
                threshold_fallback = False
                break
        credible = inclusion_probability >= inclusion_threshold
        realized_expected_fdr = float(np.mean(1.0 - inclusion_probability[credible])) if credible.any() else 0.0
        model_coefficient = np.where(credible, conditional_mean, 0.0)
        expected_comparison = np.exp(intercept + model_coefficient - np.max(intercept + model_coefficient))
        expected_comparison = expected_comparison / expected_comparison.sum() * mean_total
        log2_fold_change = np.log2(expected_comparison / expected_reference)
        lower = np.zeros(len(cell_types), dtype=float)
        upper = np.zeros(len(cell_types), dtype=float)
        posterior_sd = np.zeros(len(cell_types), dtype=float)
        for index in range(len(cell_types)):
            selected_samples = focal_samples[included[:, index], index]
            if selected_samples.size:
                lower[index], upper[index] = np.quantile(selected_samples, [0.03, 0.97])
                posterior_sd[index] = float(np.std(selected_samples, ddof=1)) if selected_samples.size > 1 else 0.0
        table = pd.DataFrame(
            {
                "contrast": contrast,
                "condition_key": condition_key,
                "reference_condition": reference_condition,
                "comparison_condition": comparison_condition,
                "cell_type": cell_types,
                "reference_cell_type": resolved_reference,
                "model_coefficient": model_coefficient,
                "hdi_lower": lower,
                "hdi_upper": upper,
                "posterior_sd": posterior_sd,
                "inclusion_probability": inclusion_probability,
                "credible_effect": credible.astype(bool),
                "estimated" + "_fdr": float(flat_expected_fdr),
                "inclusion_probability_threshold": inclusion_threshold,
                "realized_expected_fdr": realized_expected_fdr,
                "expected_count_reference": expected_reference,
                "expected_count_comparison": expected_comparison,
                "compositional_log2_fold_change": log2_fold_change,
                "reference_constraint": [index == resolved_reference_index for index in range(len(cell_types))],
            },
            columns=SCCODA_TABLE_COLUMNS,
        )
        credible_records = table.loc[table["credible_effect"], ["cell_type", "model_coefficient"]].to_dict("records")
        selection = {
            "method": "posterior_expected_fdr",
            "requested_expected_fdr": float(flat_expected_fdr),
            "strict_expected_fdr_inequality": True,
            "nonzero_tolerance": 1e-3,
            "threshold_floor_decimal_places": 3,
            "inclusion_probability_threshold": inclusion_threshold,
            "candidate_expected_fdr_before_floor": threshold_expected_fdr,
            "realized_expected_fdr_after_floor": realized_expected_fdr,
            "fallback_threshold_used": threshold_fallback,
            "credible_effect_count": int(credible.sum()),
        }
        method_name = "scCODA flat Dirichlet-multinomial composition model"
        results = {
            "row_count": int(len(table)),
            "credible_effect_count": int(credible.sum()),
            "key_effects": credible_records,
            "writing_summary": (
                f"Under Pertpy scCODA 1.3.0 for {contrast} with {resolved_reference!r} as the compositional "
                f"reference, {int(credible.sum())} of {len(cell_types)} cell types met the posterior expected-FDR "
                f"{flat_expected_fdr:g} credibility rule."
            ),
        }
    else:
        node_names = hierarchy["node_names"]
        ancestor_matrix = hierarchy["ancestor_matrix"]
        reference_nodes = hierarchy["reference_nodes"]
        reference_indices = [node_names.index(name) for name in reference_nodes]
        backend_node_names = contract.get("node_names")
        if backend_node_names != node_names:
            raise RuntimeError(f"{operation} backend hierarchy node names differ from the canonical manifest.")
        if contract.get("reference_nodes") != reference_nodes:
            raise RuntimeError(f"{operation} backend reference path differs from the canonical hierarchy path.")
        if contract.get("reference_index") != reference_indices:
            raise RuntimeError(f"{operation} backend reference-node indices differ from the canonical hierarchy path.")
        backend_ancestor = np.asarray(contract.get("ancestor_matrix"), dtype=float)
        if backend_ancestor.shape != ancestor_matrix.shape or not bool(
            np.array_equal(backend_ancestor, ancestor_matrix)
        ):
            raise RuntimeError(f"{operation} backend ancestor matrix differs from the canonical hierarchy contract.")
        if not bool(np.isin(backend_ancestor, [0.0, 1.0]).all()):
            raise RuntimeError(f"{operation} backend ancestor matrix is not binary.")
        fixed_values = {
            "lambda_0": 50.0,
            "lambda_1": 5.0,
            "theta": 0.5,
            "phi": aggregation_bias,
        }
        for key, expected in fixed_values.items():
            observed = float(contract.get(key, math.nan))
            if not math.isfinite(observed) or observed != expected:
                raise RuntimeError(f"{operation} backend prior {key!r} is {observed!r}; expected {expected!r}.")
        node_leaf_counts = ancestor_matrix.sum(axis=0).astype(float)
        nonreference_indices = [index for index in range(len(node_names)) if index not in reference_indices]
        backend_node_leaves = finite_array(
            contract.get("node_leaves_nonreference"),
            "non-reference node leaf counts",
            shape=(len(nonreference_indices),),
        )
        if not bool(np.array_equal(backend_node_leaves, node_leaf_counts[nonreference_indices])):
            raise RuntimeError(f"{operation} backend node-leaf counts differ from the canonical ancestor matrix.")
        scale = 1.0 / (
            1.0
            + np.exp(
                np.clip(
                    -aggregation_bias * (node_leaf_counts / len(cell_types) - 0.5),
                    -700.0,
                    700.0,
                )
            )
        )
        lambda_1_scaled = 2.0 * 5.0 * scale
        if not bool(np.isfinite(lambda_1_scaled).all()):
            raise RuntimeError(f"{operation} produced non-finite hierarchy-scaled slab penalties.")
        backend_scaled = finite_array(
            contract.get("lambda_1_scaled_nonreference"),
            "non-reference scaled lambda_1",
            shape=(len(nonreference_indices),),
        )
        if not bool(np.allclose(backend_scaled, lambda_1_scaled[nonreference_indices], atol=1e-12, rtol=1e-12)):
            raise RuntimeError(f"{operation} backend hierarchy-scaled slab penalties violate the pinned Pertpy rule.")
        node_samples = finite_array(
            backend_result.get("posterior_focal_node_effect_samples"),
            "posterior focal node-effect samples",
            shape=(num_samples, len(node_names)),
        )
        for index in reference_indices:
            if not bool(np.allclose(node_samples[:, index], 0.0, atol=1e-12, rtol=0.0)):
                raise RuntimeError(f"{operation} backend violated the reference-path zero-effect constraint.")
        theta_samples = finite_array(
            backend_result.get("posterior_theta_samples"),
            "posterior theta samples",
            shape=(num_samples,),
        )
        if not bool(np.allclose(theta_samples, 0.5, atol=1e-12, rtol=0.0)):
            raise RuntimeError(f"{operation} backend did not hold theta fixed at 0.5.")
        theta_realized = float(np.median(theta_samples))
        node_median = np.median(node_samples, axis=0)
        node_lower = np.quantile(node_samples, 0.03, axis=0)
        node_upper = np.quantile(node_samples, 0.97, axis=0)
        node_sd = np.std(node_samples, axis=0, ddof=1)
        node_delta = np.zeros(len(node_names), dtype=float)
        for index in nonreference_indices:
            scaled = float(lambda_1_scaled[index])
            posterior_slab_probability = (theta_realized * scaled / 2.0) / (
                (theta_realized * scaled / 2.0) + ((1.0 - theta_realized) * 50.0 / 2.0)
            )
            node_delta[index] = 1.0 / (50.0 - scaled) * math.log(1.0 / posterior_slab_probability - 1.0)
        if not bool(np.isfinite(node_delta).all()) or bool((node_delta < 0).any()):
            raise RuntimeError(f"{operation} produced invalid tascCODA Delta selection boundaries.")
        node_credible = np.abs(node_median) > node_delta
        node_credible[reference_indices] = False
        selected_node_effect = np.where(node_credible, node_median, 0.0)
        derived_leaf_effect = ancestor_matrix @ selected_node_effect
        raw_leaf_samples = node_samples @ ancestor_matrix.T
        leaf_median = np.median(raw_leaf_samples, axis=0)
        leaf_lower = np.quantile(raw_leaf_samples, 0.03, axis=0)
        leaf_upper = np.quantile(raw_leaf_samples, 0.97, axis=0)
        leaf_sd = np.std(raw_leaf_samples, axis=0, ddof=1)
        if not bool(np.allclose(derived_leaf_effect, ancestor_matrix @ selected_node_effect, atol=1e-12, rtol=0.0)):
            raise RuntimeError(f"{operation} failed its independent node-to-leaf propagation check.")
        expected_comparison = np.exp(intercept + derived_leaf_effect - np.max(intercept + derived_leaf_effect))
        expected_comparison = expected_comparison / expected_comparison.sum() * mean_total
        log2_fold_change = np.log2(expected_comparison / expected_reference)
        record_by_name = {record["name"]: record for record in hierarchy["node_records"]}
        rows = []
        for index, node in enumerate(node_names):
            record = record_by_name[node]
            rows.append(
                {
                    "effect_scope": "hierarchy_node",
                    "contrast": contrast,
                    "condition_key": condition_key,
                    "reference_condition": reference_condition,
                    "comparison_condition": comparison_condition,
                    "reference_cell_type": resolved_reference,
                    "effect_name": node,
                    "hierarchy_level": record["level"],
                    "descendant_leaf_count": len(record["descendant_leaves"]),
                    "descendant_leaves_json": json.dumps(
                        record["descendant_leaves"], ensure_ascii=False, separators=(",", ":")
                    ),
                    "model_effect": float(selected_node_effect[index]),
                    "posterior_median": float(node_median[index]),
                    "hdi_lower": float(node_lower[index]),
                    "hdi_upper": float(node_upper[index]),
                    "posterior_sd": float(node_sd[index]),
                    "selection_delta": float(node_delta[index]),
                    "credible_effect": bool(node_credible[index]),
                    "selection_basis": "direct_node_abs_median_strictly_greater_than_delta",
                    "expected_count_reference": pd.NA,
                    "expected_count_comparison": pd.NA,
                    "compositional_log2_fold_change": pd.NA,
                    "reference_constraint": node in reference_nodes,
                }
            )
        for index, leaf in enumerate(cell_types):
            rows.append(
                {
                    "effect_scope": "derived_leaf",
                    "contrast": contrast,
                    "condition_key": condition_key,
                    "reference_condition": reference_condition,
                    "comparison_condition": comparison_condition,
                    "reference_cell_type": resolved_reference,
                    "effect_name": leaf,
                    "hierarchy_level": annotation_key,
                    "descendant_leaf_count": 1,
                    "descendant_leaves_json": json.dumps([leaf], ensure_ascii=False, separators=(",", ":")),
                    "model_effect": float(derived_leaf_effect[index]),
                    "posterior_median": float(leaf_median[index]),
                    "hdi_lower": float(leaf_lower[index]),
                    "hdi_upper": float(leaf_upper[index]),
                    "posterior_sd": float(leaf_sd[index]),
                    "selection_delta": pd.NA,
                    "credible_effect": pd.NA,
                    "selection_basis": "derived_sum_of_selected_hierarchy_node_effects",
                    "expected_count_reference": float(expected_reference[index]),
                    "expected_count_comparison": float(expected_comparison[index]),
                    "compositional_log2_fold_change": float(log2_fold_change[index]),
                    "reference_constraint": leaf == resolved_reference,
                }
            )
        table = pd.DataFrame(rows, columns=TASCCODA_TABLE_COLUMNS)
        selection = {
            "method": "tree_adaptive_spike_and_slab_lasso",
            "lambda_0": 50.0,
            "lambda_1": 5.0,
            "theta": 0.5,
            "aggregation_bias": aggregation_bias,
            "scaled_lambda_1": [float(value) for value in lambda_1_scaled],
            "delta": [float(value) for value in node_delta],
            "strict_delta_inequality": True,
            "direct_credible_node_count": int(node_credible.sum()),
        }
        method_name = "tascCODA tree-aggregated Dirichlet-multinomial composition model"
        direct_records = [
            {"hierarchy_node": node_names[index], "model_effect": float(selected_node_effect[index])}
            for index in range(len(node_names))
            if node_credible[index]
        ]
        results = {
            "row_count": int(len(table)),
            "direct_hierarchy_node_rows": len(node_names),
            "derived_leaf_rows": len(cell_types),
            "direct_credible_node_count": int(node_credible.sum()),
            "key_effects": direct_records,
            "writing_summary": (
                f"Under Pertpy tascCODA 1.3.0 for {contrast} with {resolved_reference!r} and its ancestry "
                f"constrained as the reference path, {int(node_credible.sum())} of {len(node_names)} modeled "
                "hierarchy nodes met the direct Delta credibility rule; leaf effects are derived consequences."
            ),
        }
        references.append(
            {
                "citation": (
                    "Ostner J, Carcy B, Müller CL. tascCODA: Bayesian tree-aggregated analysis of compositional "
                    "amplicon and single-cell data. Frontiers in Genetics. 2021;12:766313."
                ),
                "doi": "10.3389/fgene.2021.766313",
                "url": "https://doi.org/10.3389/fgene.2021.766313",
                "kind": "method",
            }
        )

    count_records = []
    display_by_internal = dict(zip(sample_internal_order, context["sample_display_order"], strict=True))
    condition_by_internal = dict(zip(sample_metadata["sample_internal"], sample_metadata["condition"], strict=True))
    for sample_index, sample in enumerate(sample_internal_order):
        for cell_index, cell_type in enumerate(cell_types):
            count_records.append(
                {
                    "sample": display_by_internal[sample],
                    "condition": condition_by_internal[sample],
                    "cell_type": cell_type,
                    "count": int(count_matrix[sample_index, cell_index]),
                }
            )
    parameters = {
        "sample_key": sample_key,
        "annotation_key": annotation_key,
        "annotation_status": annotation_status,
        "condition_key": condition_key,
        "reference_condition": reference_condition,
        "comparison_condition": comparison_condition,
        "adjustment_covariate_keys": adjustment_keys,
        "reference_cell_type": reference_cell_type,
        "num_samples": num_samples,
        "num_warmup": num_warmup,
        "random_seed": random_seed,
    }
    if method == "sccoda":
        parameters["estimated" + "_fdr"] = flat_expected_fdr
    else:
        parameters["hierarchy_keys"] = hierarchy_keys
        parameters["aggregation_bias"] = aggregation_bias
    summary = {
        "schema_version": 1,
        "node_id": (
            "OpenBioSingleCellSccodaDifferentialComposition"
            if method == "sccoda"
            else "OpenBioSingleCellTasccodaDifferentialComposition"
        ),
        "operation": "sccoda_differential_composition" if method == "sccoda" else "tasccoda_differential_composition",
        "method": method_name,
        "status": diagnostic_status,
        "inference_scope": "formal" if annotation_status == "curated" else "exploratory",
        "parameters": parameters,
        "question": {
            "condition_key": condition_key,
            "contrast": contrast,
            "reference_condition": reference_condition,
            "comparison_condition": comparison_condition,
            "sample_is_experimental_unit": True,
            "backend_runtime": backend_runtime,
            "adjustment_covariates": adjustment_keys,
        },
        "input": {
            "cells": int(adata.n_obs),
            "features_not_modeled": int(adata.n_vars),
            "samples": len(sample_internal_order),
            "samples_per_condition": {key: insufficient[key] for key in sorted(insufficient)},
            "cell_types": cell_types,
            "count_table_shape": [len(sample_internal_order), len(cell_types)],
            "count_table_fingerprint": count_fingerprint,
            "sample_total_cells": numeric_summary(sample_totals, "Sample totals", positive=True),
            "zero_fraction": float(np.mean(count_matrix == 0)),
            "annotation_key": annotation_key,
            "annotation_status": annotation_status,
            "complete_sample_cell_type_counts": count_records,
        },
        "design": {
            "formula": design_formula,
            "columns": ["Intercept", *design_columns],
            "focal_column": focal_column,
            "adjustments": adjustment_report,
            "dimensions": [int(value) for value in design_matrix.shape],
            "rank": design_rank,
            "residual_degrees_of_freedom": residual_df,
            "condition_number": condition_number,
            "full_rank": True,
            "expected_count_baseline": {
                "condition": reference_condition,
                "continuous_adjustments": "mean-centered value 0",
                "categorical_adjustments": "declared first category",
            },
        },
        "reference": {
            "requested": reference_cell_type,
            "resolved": resolved_reference,
            "automatic_absence_threshold": 0.05,
            "candidate_rule_uses_strict_inequality": True,
            "candidates": reference_candidates,
        },
        "model": {
            "family": "Dirichlet-multinomial",
            "method": method,
            "sampler": "NUTS",
            "posterior_draws": num_samples,
            "warmup_draws": num_warmup,
            "chain_count": 1,
            "random_seed": random_seed,
            "zero_pseudocount": 0.5,
            "automatic_reference_absence_threshold": 0.05,
            "sample_is_experimental_unit": True,
            **(
                {}
                if method == "sccoda"
                else {
                    "lambda_0": 50.0,
                    "lambda_1": 5.0,
                    "theta": 0.5,
                    "aggregation_bias": aggregation_bias,
                }
            ),
        },
        "selection": selection,
        "diagnostics": diagnostics,
        "results": results,
        "warnings": warnings,
        "limitations": [
            "The one-chain Pertpy fit cannot establish between-chain mixing; R-hat is unavailable.",
            "Pertpy 1.3.0 does not retain divergence counts for this audited runner.",
            "Composition effects are relative to the declared cell-type reference and are not absolute abundance effects.",
            "Cells are subsamples within biological Samples and do not increase the replicate count.",
            "The model has no random-effect or repeated-measures contract.",
            "Observed Condition associations are not causal without an appropriate experimental or causal design.",
        ],
        "references": references,
        "software_versions": {key: software_versions[key] for key in sorted(software_versions)},
    }
    if method == "tasccoda":
        summary["hierarchy"] = {
            "ancestor_keys": hierarchy_keys,
            "leaf_key": annotation_key,
            "levels_passed_to_pertpy": hierarchy["levels_passed_to_pertpy"],
            "root": hierarchy["root"],
            "edges": hierarchy["edges"],
            "node_names": hierarchy["node_names"],
            "ancestor_matrix_shape": [int(value) for value in hierarchy["ancestor_matrix"].shape],
            "fingerprint": hierarchy["fingerprint"],
            "reference_path": hierarchy["reference_path"],
            "declared_reference_path": hierarchy["declared_reference_path"],
            "reference_nodes": hierarchy["reference_nodes"],
            "derived_leaf_effects_are_propagated": True,
        }
        summary["limitations"].extend(
            [
                "Hierarchy-node credibility belongs only to direct modeled nodes; derived leaf rows are not independent discoveries.",
                "The result depends on the declared hierarchy and the fixed tree-adaptive prior, including aggregation bias.",
            ]
        )
    json.dumps(summary, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)
    return table, summary


def _run_composition_model_with_preserved_rng(*args: Any, **kwargs: Any) -> tuple[Any, dict[str, Any]]:
    import random

    import numpy as np

    with _COMPOSITION_RNG_LOCK:
        python_state = random.getstate()
        numpy_state = np.random.get_state()
        try:
            return _standalone_composition_model_impl(*args, **kwargs)
        finally:
            random.setstate(python_state)
            np.random.set_state(numpy_state)


def run_sccoda_differential_composition(
    adata,
    *,
    sample_key="sample",
    annotation_key="cell_type",
    annotation_status="provisional",
    condition_key="condition",
    reference_condition="",
    comparison_condition="",
    adjustment_covariate_keys_json="[]",
    reference_cell_type="automatic",
    estimated_fdr=0.05,
    num_samples=10000,
    num_warmup=1000,
    random_seed=0,
    openbio_version="unknown",
    _backend=None,
):
    return _run_composition_model_with_preserved_rng(
        adata,
        method="sccoda",
        sample_key=sample_key,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        condition_key=condition_key,
        reference_condition=reference_condition,
        comparison_condition=comparison_condition,
        adjustment_covariate_keys_json=adjustment_covariate_keys_json,
        reference_cell_type=reference_cell_type,
        flat_expected_fdr=estimated_fdr,
        num_samples=num_samples,
        num_warmup=num_warmup,
        random_seed=random_seed,
        openbio_version=openbio_version,
        _backend=_backend,
    )


def run_tasccoda_differential_composition(
    adata,
    *,
    sample_key="sample",
    annotation_key="cell_type",
    annotation_status="provisional",
    hierarchy_keys_json="[]",
    condition_key="condition",
    reference_condition="",
    comparison_condition="",
    adjustment_covariate_keys_json="[]",
    reference_cell_type="automatic",
    aggregation_bias=0.0,
    num_samples=10000,
    num_warmup=1000,
    random_seed=0,
    openbio_version="unknown",
    _backend=None,
):
    return _run_composition_model_with_preserved_rng(
        adata,
        method="tasccoda",
        sample_key=sample_key,
        annotation_key=annotation_key,
        annotation_status=annotation_status,
        condition_key=condition_key,
        reference_condition=reference_condition,
        comparison_condition=comparison_condition,
        adjustment_covariate_keys_json=adjustment_covariate_keys_json,
        reference_cell_type=reference_cell_type,
        hierarchy_keys_json=hierarchy_keys_json,
        aggregation_bias=aggregation_bias,
        num_samples=num_samples,
        num_warmup=num_warmup,
        random_seed=random_seed,
        openbio_version=openbio_version,
        _backend=_backend,
    )


_PINNED_BACKEND_SOURCE_CLASS = _PinnedPertpyBackend


def _pinned_backend_source(method: str) -> str:
    method_names = [
        "__init__",
        "_require_parameters",
        "_prepare_model",
        "_diagnostics",
        "_runtime_details",
        "fit_sccoda" if method == "sccoda" else "fit_tasccoda",
    ]
    methods = []
    for name in method_names:
        source = textwrap.dedent(inspect.getsource(getattr(_PINNED_BACKEND_SOURCE_CLASS, name))).strip()
        methods.append(textwrap.indent(source, "    "))
    return (
        "class _PinnedPertpyBackend:\n"
        '    """Narrow, fail-closed adapter for the audited Pertpy 1.3.0 composition API."""\n\n' + "\n\n".join(methods)
    )


def _composition_code(method: str, parameters: dict[str, Any]) -> str:
    backend_source = _pinned_backend_source(method)
    implementation = textwrap.dedent(inspect.getsource(_standalone_composition_model_impl))
    rng_wrapper = textwrap.dedent(inspect.getsource(_run_composition_model_with_preserved_rng))
    function_name = (
        "run_sccoda_differential_composition" if method == "sccoda" else "run_tasccoda_differential_composition"
    )
    ordered_names = [
        "sample_key",
        "annotation_key",
        "annotation_status",
        "condition_key",
        "reference_condition",
        "comparison_condition",
        "adjustment_covariate_keys_json",
        "reference_cell_type",
    ]
    if method == "sccoda":
        ordered_names.append("estimated_fdr")
    else:
        ordered_names.extend(["hierarchy_keys_json", "aggregation_bias"])
    ordered_names.extend(["num_samples", "num_warmup", "random_seed", "openbio_version"])
    missing = [name for name in ordered_names if name not in parameters]
    if missing:
        raise TypeError(f"Missing code-generation parameters: {missing}.")
    parameter_lines = []
    for name in ordered_names:
        internal_name = "flat_expected_fdr" if name == "estimated_fdr" else name
        parameter_lines.append(f"        {internal_name}={parameters[name]!r},")
    parameter_lines_text = "\n".join(parameter_lines)
    public_source = f"""def {function_name}(adata, *, _backend=None):
    return _run_composition_model_with_preserved_rng(
        adata,
        method={method!r},
{parameter_lines_text}
        _backend=_backend,
    )"""
    return "\n\n".join(
        [
            "from __future__ import annotations\n\nimport threading\nfrom typing import Any",
            (
                "_COMPOSITION_RNG_LOCK = threading.RLock()\n\n"
                + (
                    f"SCCODA_TABLE_COLUMNS = {SCCODA_TABLE_COLUMNS!r}"
                    if method == "sccoda"
                    else f"TASCCODA_TABLE_COLUMNS = {TASCCODA_TABLE_COLUMNS!r}"
                )
            ),
            backend_source,
            implementation,
            rng_wrapper,
            public_source,
        ]
    )


def sccoda_differential_composition_code(**parameters: Any) -> str:
    return _composition_code("sccoda", parameters)


def tasccoda_differential_composition_code(**parameters: Any) -> str:
    return _composition_code("tasccoda", parameters)


__all__ = [
    "SCCODA_TABLE_COLUMNS",
    "TASCCODA_TABLE_COLUMNS",
    "run_sccoda_differential_composition",
    "run_tasccoda_differential_composition",
    "sccoda_differential_composition_code",
    "tasccoda_differential_composition_code",
]
