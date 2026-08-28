# Diffusion Map: module design review

## Decision

**Retain and deepen `OpenBioSingleCellDiffusionMap` as an atomic named-graph-to-diffusion-basis module.** Do not merge it with DPT: the basis is independently inspectable/reusable, whereas DPT adds a biological root and oriented ordering.

Inputs are `adata`, visible `neighbors_key="neighbors"`, visible `n_comps=15`, advanced `overwrite_existing=False`, and advanced `random_seed=0`. Outputs are `adata`, strict JSON `summary`, and equivalent `code`. Keep kernel/eigensolver internals hidden at the audited Scanpy policy.

## Deep graph/state contract

Resolve the named graph through the shared Batch05 graph module and validate its metadata, sparse matrices, alignment, finite nonnegative weights, symmetry policy, connected components, and observation fingerprint. Reject duplicate observations, empty/tiny inputs, invalid component bounds, and output-bundle collisions before Scanpy. Preserve all input matrices and unrelated state.

After execution, require `X_diffmap.shape == (n_obs, n_comps)` and `diffmap_evals.shape == (n_comps,)`, finite real values, and a non-increasing eigenvalue sequence within numerical tolerance. Store OpenBio provenance beside the Scanpy bundle: graph key/fingerprint, obs fingerprint, n_comps, seed, Scanpy version and output fingerprint. DPT must verify this record instead of trusting a coincidentally named matrix.

`summary` contains graph diagnostics, component sizes/isolates, exact output shapes/dtypes, all eigenvalues, informative-component ranges/quantiles, stationary-component disclosure, seed, references, versions and limitations. It must say that no root, pseudotime, branch, cluster, `Sample`, or `Condition` analysis occurred.

`code` defines a standalone function that performs the same graph/collision checks, calls `sc.tl.diffmap(..., neighbors_key=..., random_state=..., copy=False)`, validates the bundle, writes portable provenance, and returns `(output_adata, summary_dict)`.

## Migration and tests

The old schema maps directly, with `overwrite_existing=False` added. Workflows containing DPT without Diffusion Map must be migrated by inserting this node on the same graph; migrations must not silently run a different graph.

Tests cover default/custom graph keys, malformed pointers/matrices, disconnected/isolated graphs, too few cells/components, collisions/overwrite, repeatability, input immutability, postcondition tampering, graph/provenance fingerprint mismatch, strict JSON, compiled source, and runtime/code equivalence.

## 2026-08-28 adversarial implementation audit

The first implementation was not release-ready. It fingerprinted an averaged approximately symmetric connectivity matrix while Scanpy consumed the unmodified matrix; two distinct computation inputs could therefore share one graph fingerprint and yield different diffusion bases. Shape, finiteness, and eigenvalue ordering also allowed a signature-compatible backend to return unrelated zero coordinates and receive trusted provenance. Replacing a Diffusion Map could leave an existing DPT bundle silently stale, and generated source embedded all three trajectory analyses.

Release closure requires one canonical matrix to be both fingerprinted and consumed, independent transition-operator eigenpair residual checks, fail-closed downstream DPT invalidation, a Diffusion-Map-only standalone dependency closure, and malicious-backend plus migration regression tests. The legacy positional schema must be migrated atomically before registration can be treated as safe.

## 2026-08-28 implemented closure

The validated canonical connectivity and distance matrices are now written onto the output copy before the public Scanpy call, so the computation input, postcondition input, and graph fingerprint are identical. The implementation independently reconstructs Scanpy 1.12.3's density-normalized symmetric transition operator and verifies eigenvalue bounds, eigenvector orthonormality, and every relative eigenpair residual before writing provenance. A downstream DPT bundle, including an orphaned `iroot`, blocks Diffusion Map replacement rather than being left stale. NumPy RNG and print state are restored on success and failure.

Postconditions compare the exact stored sparse connectivity/distance fingerprints against the consumed canonical matrices before applying any symmetry tolerance again. Thus a backend cannot hide an equal-and-opposite asymmetric edit behind the input canonicalization average.

Standalone source now contains only the Diffusion Map dependency closure. Whole-workflow migration recognizes exact old positional `n_comps,neighbors_key,random_seed`, rewrites the current order with `overwrite_existing=false`, preserves slot-zero fan-out and connected/exposed controls, appends `summary/code`, and is idempotent. Partial/mixed schemas and malformed links fail during all-graph preflight before any workflow mutation. Real Scanpy, asymmetric-within-policy canonicalization, malicious orthonormal-but-wrong eigenpairs, stale-DPT rejection, strict JSON, generated-code parity, object/array-link, subgraph, exposure, and idempotence regressions cover the closure.

Report size is explicitly bounded: all at-most-4,096 eigenvalues remain disclosed, while per-component distribution summaries are capped at 100 with total/truncation metadata; named-graph component-size lists are likewise capped at 100 and isolated-cell examples at 20.
