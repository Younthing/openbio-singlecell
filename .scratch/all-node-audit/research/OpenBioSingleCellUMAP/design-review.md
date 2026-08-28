# UMAP: module design review

## Decision

**Keep `OpenBioSingleCellUMAP` as an atomic graph-to-two-dimensional-embedding module.** Do not merge it with Neighbors: graph construction is reused by UMAP, graph drawing, and Leiden, while UMAP parameter/seed changes should not rebuild the scientific graph. Do not merge it with plotting, because coordinates are reusable data and plots have separate presentation choices.

The seam is:

```text
one named validated neighbor graph -> one named 2D UMAP coordinate representation
```

The deletion test confirms depth: removing this module would spread graph-key resolution, UMAP prerequisites, collision checks, reproducibility disclosure, and coordinate validation into every workflow/plot consumer.

## Proposed interface

- `adata`, required;
- `neighbors_key`, visible because it selects the scientific graph, default `neighbors`;
- `min_dist`, visible, default `0.5`;
- `spread`, visible, default `1.0`;
- `key_added`, advanced, default `X_umap` or an interface-normalized suffix with unambiguous storage semantics;
- `overwrite_existing`, advanced, default `False`;
- `random_seed`, advanced, default `0`.

Keep `n_components=2`, `init_pos="spectral"`, `maxiter=None`, `alpha=1`, `gamma=1`, `negative_sample_rate=5`, `a=None`, `b=None`, CPU `method="umap"`, and `copy=False` hidden. They define the supported visualization adapter, not ordinary decisions at this seam. If initialization is later exposed, use a bounded mode plus validated coordinate key; do not restore an unconstrained string that may silently imply PAGA prerequisites.

`key_added` needs a clear wrapper convention because Scanpy's default uses `X_umap` in obsm but `umap` in uns, while a custom key is used identically in both. The summary must report both resolved keys, not merely the input string.

## Deep implementation

Before calling Scanpy, the implementation should validate AnnData, unique observations, finite parameter values, graph metadata indirection, graph shape/nonnegativity/finiteness, edge presence, and output collisions. It should calculate connected components and isolated cells before embedding so failures and limitations are reported locally.

After Scanpy returns, verify a finite `n_obs × 2` coordinate array and expected uns parameter record. Preserve X, raw, layers, obs, var, the neighbor graph, and unrelated embeddings. The node must never rebuild neighbors or infer a representation from X.

This creates leverage through a small interface while graph validation, Scanpy storage asymmetry, hidden optimization policy, and diagnostics stay within the module.

## Summary contract

`summary` must be strict JSON and include:

- UMAP and Scanpy method references;
- resolved neighbors uns/connectivity keys and graph shape/nnz;
- connected-component count/sizes and isolated observations;
- `min_dist`, `spread`, derived/default `a/b` policy, 2D output, spectral initialization, optimizer defaults, CPU method, and seed;
- coordinate key and parameter key, shape, finite status, per-axis range/median/IQR, and overwrite status;
- warnings for disconnected graphs, small components, parameter sensitivity, and stochastic/version sensitivity;
- Scanpy, umap-learn, scikit-learn, AnnData, NumPy, SciPy, and OpenBio versions;
- limitations that the embedding is exploratory, axes are arbitrary, visual density/separation is parameter-dependent, and it is neither `Technical batch` integration nor `Sample`-level inference.

No claim about preserved “global biology,” cluster truth, or condition separation is justified by the picture alone.

## Equivalent code and tests

`code` should define a standalone function that validates the named graph and collision policy, copies AnnData, calls `sc.tl.umap` with all hidden/default values explicit, validates the `n_obs × 2` output, and returns the copy. It must not call Neighbors, plotting functions, ComfyUI, or OpenBio history helpers.

Interface tests should cover default/custom graph keys, missing/malformed graphs, disconnected/isolated graphs, invalid `min_dist`/`spread`, output collisions and overwrite, seed repeatability, input preservation, custom key storage semantics, strict JSON, and generated-code equality for X/layers/raw/uns/obsm/obsp.

The shared graph seam must expose an explicit `require_undirected`/diagonal policy. UMAP selects the permissive directed policy and preserves the supplied graph, while Leiden and force-directed drawing retain their reviewed undirected contract. Runtime and generated code accept zero-variable graph-only AnnData and emit the same unusual-parameter warnings.
