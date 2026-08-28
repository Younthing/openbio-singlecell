# Harmony Integration: module design review

## Current module

At the start of the original audit, the node copied AnnData, checked only that one `batch_key` and one `obsm` basis existed, called `scanpy.external.pp.harmony_integrate`, stored one adjusted basis, appended generic history, and returned only AnnData. The refactored node now calls the reviewed harmonypy 2.0.x interface directly and returns AnnData, `summary`, and `code`; this follow-up removes advisory shape/size gates and fixes source-basis preservation reporting.

There is no dedicated Harmony test. Registration checks only prove that the node exists. No example workflow currently exercises Harmony. There is no consumer coupled to the Harmony node as an object; generic Neighbors consumes the adjusted `obsm` key by string.

## Correctness and interface problems

1. `max_iter_kmeans=100` overrides the current harmonypy default of 4 without a documented scientific reason.
2. Random initialization is left to an implicit upstream default; there is no seed in history or output.
3. Official Scanpy accepts multiple batch covariates, while the current interface admits only one string.
4. The input basis is not checked for numeric dtype, finite values, row alignment, useful dimension, or duplicate observation identifiers.
5. Missing/null/blank batch labels, a one-level column, tiny levels, and collision with the output key are not handled.
6. `batch_key` language invites confusion between Technical batch, Sample, and Condition. The node neither protects nor warns about biological variables.
7. The original and corrected embeddings are not compared even at a descriptive level, and convergence/objective evidence is discarded by the Scanpy wrapper.
8. The implementation does not isolate BLAS resource use or declare CPU/thread policy.
9. There is no strict `summary` or equivalent `code` output.

## Decision: keep as one atomic embedding correction and deepen it

Keep Harmony separate from PCA and Neighbors. Its atomic result is one corrected PCA embedding. Merging it with PCA would prevent reuse of a reviewed PCA result; merging it with Neighbors would bind integration to one graph parameterization and obscure the corrected coordinates.

Do not add a persistent Harmony object wire. Downstream callers need the corrected coordinates, not mutable soft clusters or backend caches. Objective history is small reporting metadata and can be materialized into `uns`/summary. The deletion test supports this seam: deleting the node would force each caller to repeat covariate validation, upstream-version adaptation, orientation checks, seed/resource policy, output collision handling, and report wording.

Use the reviewed harmonypy 2.0.x public contract directly behind an internal adapter so objective histories remain
available. Keep 2.0.0 as the release pin; at runtime require agreeing module/distribution 2.0.x identities,
runtime-check the required public keywords, require
raw `result.Z_corr.shape == input_basis.shape == (n_obs, n_components)`, and store it without transposition. Do not
infer orientation from whichever shape appears plausible, including the ambiguous square case. The historical Scanpy
wrapper's `.T` is a 0.x adapter and must not be used as evidence for 2.0.0 behavior. There is no supported fallback to
0.x or an unversioned future interface.

## Target interface

- Visible: `adata`, comma-separated `technical_batch_keys` (default `batch`), and PCA `basis` (default `X_pca`).
- Advanced scientific: `theta` (automatic official default or explicit scalar/per-key values), `lambda` mode/value, positive `sigma`, `nclust` auto/custom, and `tau`; no arbitrary UI maxima.
- Advanced reproducibility/state: `adjusted_basis`, `overwrite_existing`, `max_iter_harmony=10`, `max_iter_kmeans=4`, `random_seed=0`.
- Hidden fixed policy: non-mutating copy, explicit finite/axis validation, `verbose=False`, CPU harmonypy backend, single-process call, and preferably `ncores=1` for a reproducible/resource-bounded baseline.
- Outputs: copied `adata`, strict JSON `summary`, equivalent Python source `code`.

If per-covariate `theta`/lambda would make the normal interface shallow and error-prone, use a conditional automatic/custom control and accept either one scalar for all keys or exactly one comma-separated value per key. Do not expose every harmonypy implementation tuning knob merely because it exists. The selected fixed defaults still appear in summary and generated code.

## Domain contract

Rename the conceptual input to `technical_batch_keys`. Saved workflow migration may map the old `batch_key` value, but the migration and report must call it a user declaration, not prove that the column is technical.

An optional protected-column check may accept `sample_key` and `condition_key` only if those ordinary string inputs are already available from Core Study Parameters; any overlap with correction keys is an error. Their absence does not authorize the report to claim Sample or Condition preservation. The summary always warns that confounding can cause biological signal removal and that the corrected embedding is exploratory, not a Condition contrast.

Harmony's deep module interface is the declared `.obsm` basis plus `.obs` covariates; current `X` and the variable axis are outside this seam. Require at least one finite aligned component, warn on one-dimensional or zero-variance bases, and verify the backend's finite exact-shape output. Resolve automatic cluster count with a floor of one, accept custom `1 <= nclust <= n_obs`, and pass a length-`nclust` sigma vector to avoid harmonypy 2.0.0's scalar-`sigma` one-cluster defect. `nclust > n_obs` remains a safety failure because the upstream without-replacement initializer cannot terminate after exhausting cells.

Reject `basis == adjusted_basis` even when overwrite is requested. Overwrite applies only to a distinct existing destination; otherwise the node cannot satisfy its stated preservation contract. Verify source coordinates remain unchanged before reporting `original_basis_preserved=true`.

## Summary contract

`summary` includes:

- methods: PCA-basis Harmony correction, technical covariates, public backend call, and resolved correction parameters;
- results: corrected basis location/shape and cautious descriptive changes only;
- key results: input/output dimensions, Technical-batch level counts, rounds/objectives when supported, finite checks, and displacement quantiles;
- parameters: every visible/advanced and fixed hidden value, including seed and thread policy;
- warnings: user-declared Technical batch, possible confounding/over-correction, no guarantee of biological conservation, no corrected expression matrix, and downstream graph dependence;
- references and runtime software versions.

Strict JSON forbids NumPy scalars, arrays, NaN, and Infinity. If an objective is unavailable or non-finite, use `null` plus a warning rather than a false value.

## Equivalent-code feasibility

Equivalent source is straightforward. The generated function should copy AnnData, parse/validate Technical-batch keys and PCA coordinates, validate overwrite policy, import a capability-checked harmonypy 2.0.x release, call `run_harmony` with every resolved parameter explicit, require raw `Z_corr` to be finite real numeric with exact shape `(n_obs, n_components)`, store it unchanged in `obsm`, and return AnnData. It should include the same distribution/module-version and runtime-interface checks and must not depend on ComfyUI or plugin history.

The code may omit plugin-only elapsed-time/history bookkeeping. It may optionally return a plain diagnostics dictionary, but the primary returned AnnData must match runtime scientific state. Code and runtime must use the same backend path; generating Scanpy code while runtime calls harmonypy directly would not be equivalent. Single-level and singleton correction strata are warnings in both paths, not hard gates.

Generated code must also reproduce one-dimensional/zero-variance and small-level warnings, the same `nclust`/sigma normalization, the distinct source/destination rule, and the 2.0.x (not exact-2.0.0 wording) compatibility contract.

## Failure and verification plan

- Compare dense float32/float64 PCA inputs with a direct harmonypy call under the same explicit seed/parameters.
- Assert original basis preservation, exact output shape/order, no input mutation, and overwrite refusal.
- Reject absent/malformed/non-finite bases, empty/duplicate observation axes, absent/null/blank Technical-batch metadata, duplicate keys, a destination equal to its source, `nclust` outside `[1,n_obs]`, and protected Sample/Condition overlap; do not require a nonempty or unique variable axis. Disclose one-dimensional/zero-variance bases and one-level/singleton correction strata without rejecting them.
- Cover multiple Technical-batch covariates, unbalanced/small levels, scalar/per-key penalties, and invalid penalty lengths/ranges.
- Fake-backend tests pin every explicit argument, current defaults, the capability-checked 2.0.x family, orientation,
  objective extraction, actual version reporting, and actionable optional-dependency errors. Include non-square and
  non-symmetric square `Z_corr` cases so runtime and generated code prove that they preserve the 2.0 orientation
  rather than transpose or guess it, plus warning-only one-level/singleton strata.
- Run twice with the same seed/thread policy and compare outputs within documented numerical tolerance; report that cross-platform bitwise equality is not guaranteed.
- Validate strict JSON, cautious claims, citations/software versions, generated-code compilation, and runtime/code state equivalence.
- Exercise `nclust=1`, `nclust=n_obs`, the automatic small-`n_obs` floor, and fail-before-backend `nclust>n_obs`; assert runtime/generated pass identical sigma vectors and preserve the original basis.
