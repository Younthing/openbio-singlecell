# t-SNE: module design review

## Decision

**Keep and deepen `OpenBioSingleCellTSNE` as an atomic representation-to-2D visualization module.** It should remain separate from PCA/scVI/Harmony because those representations are reusable and have distinct model contracts. It should remain separate from Neighbors and UMAP because t-SNE consumes a representation directly rather than the stored Scanpy graph.

The seam is:

```text
one explicit cell representation -> one named 2D t-SNE coordinate representation
```

It must not silently compute PCA, build a graph, cluster cells, or present map islands as biological populations.

## Proposed interface and P0 correction

- `adata`, required;
- `use_rep`, visible and explicit, default `X_pca` for the conventional PCA flow;
- `n_dimensions`, visible, default `0` meaning all available dimensions;
- `perplexity`, visible, default `30`;
- `metric`, visible bounded combo, default `euclidean`;
- `early_exaggeration`, advanced, default `12`;
- `learning_rate`, advanced, default `1000` to preserve the audited Scanpy adapter;
- `key_added`, advanced, default `X_tsne`;
- `overwrite_existing`, advanced, default `False`;
- `random_seed`, advanced, default `0`.

`n_dimensions=0` maps to `n_pcs=None`; positive values slice after checking the actual representation width. This fixes 10D scVI and removes inaccurate universal “PC” terminology. Migration should retain old positive `n_pcs` only when it fits; oversized old values should migrate to all dimensions with an explicit warning.

Keep `n_components=2`, `use_fast_tsne=False`, `n_jobs=None`, and `copy=False` hidden. Exposing MulticoreTSNE or arbitrary job counts creates backend/performance surface without a second supported scientific adapter.

## Deep implementation

The module should resolve and validate the representation itself instead of allowing Scanpy's implicit PCA branch. It should validate finite numeric values, axis alignment, unique observations, dimension selection, `0 < perplexity < n_obs`, positive advanced parameters, and output collisions before calculation.

After execution, validate a finite `n_obs × 2` result and the matching uns parameter record. Preserve all biological matrices, annotations, graph states, raw, and unrelated embeddings. Disclose but do not silently “fix” duplicate rows or degenerate coordinate ranges.

This is a deep module because a compact interface hides Scanpy/scikit version-default differences, representation slicing, collision semantics, stochasticity, and result diagnostics with locality.

## Summary contract

Strict JSON `summary` should include:

- t-SNE, Scanpy, and scikit-learn references;
- representation key, source shape, available/used dimensions, and metric;
- observation count, perplexity and perplexity-to-observation ratio;
- early exaggeration, learning rate, 2D output, backend `sklearn`, hidden `use_fast_tsne=False`, thread policy, and seed;
- output/parameter keys, overwrite status, finite coordinate checks, per-axis range/median/IQR, and duplicate-coordinate count;
- warnings for small datasets, extreme perplexity, degenerate inputs/output, and stochastic/non-convex sensitivity;
- Scanpy, scikit-learn, AnnData, NumPy, SciPy, and OpenBio versions;
- limitations about arbitrary axes, unreliable long-range distance/area, visualization-only use, and absence of `Technical batch` integration or replicate-aware `Condition` inference.

Do not report KL divergence unless it is actually obtained from the fitted backend. Scanpy discards the estimator and does not store that diagnostic in its AnnData result.

## Equivalent code and tests

`code` should be a standalone function that resolves the explicit representation, performs all preflight checks and collision handling, calls `sc.tl.tsne` with every supported/hidden parameter explicit, validates output, and returns a copied AnnData. It must reproduce warnings and must not depend on ComfyUI or OpenBio internals.

Tests should cover PCA and 10D scVI inputs, all/subset dimensions, missing/misaligned/nonfinite representations, perplexity bounds, seed repeatability, metric and advanced parameters, custom/default key asymmetry, collisions/overwrite, input preservation, strict JSON, and generated-code equivalence.

Runtime and generated source use the same representation-only precondition, allow positive perplexity below one, and warn about extreme ratios without changing the analyst's value. Fixed UI ceilings do not substitute for actual dimensional bounds.
