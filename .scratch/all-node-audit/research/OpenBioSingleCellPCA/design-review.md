# PCA: module design review

## Current module

At the start of the original audit, the node copied `AnnData`, resolved `X` or a named layer, optionally requested `var["highly_variable"]`, and called `scanpy.pp.pca` with fixed `svd_solver="arpack"`, canonical result keys, and a seed while returning only the modified copy. The refactored node already returns AnnData, strict `summary`, and equivalent `code`; this follow-up tightens its open-expert disclosure and removes non-method UI bounds.

The packaged clustering path correctly sends the reusable normalized/logarithmized `log1p_norm` layer into this node and does not require the Scale node solely to reach PCA.

## Correctness and interface problems

1. `n_comps` is not checked against the selected matrix. ARPACK requires `n_comps < min(n_obs, n_selected_vars)`; the UI's generic upper bound permits a late dependency error.
2. `use_hvg=True` checks only that the column exists. It does not require boolean dtype, any selected features, enough selected features, or disclose constant selected features.
3. The source default is `X`, even though this repository normally preserves counts in `X` and constructs `log1p_norm` for PCA. A default call can therefore analyze counts while only the example canvas is safe.
4. Expression state is not disclosed. Counts or untransformed size-normalized values should trigger a strong methodological warning, while unknown provenance needs its own limitation; none is a mathematical PCA execution error.
5. Non-finite values, empty axes, duplicate identifiers, backed/lazy arrays, and unsupported matrix types are not translated into clear interface errors.
6. Repeated execution silently replaces `obsm["X_pca"]`, `varm["PCs"]`, and `uns["pca"]`; a partial pre-existing bundle is also silently mixed/replaced.
7. Although sparse centered input is handled efficiently by Scanpy/ARPACK, persistent dense score/loadings allocation is not estimated or limited.
8. The node does not summarize explained variance, cumulative variance, or loadings, and exposes neither strict `summary` nor equivalent `code` outputs.

## Decision: keep and enhance; do not merge

Keep PCA as one atomic analysis: transform one deliberately selected expression representation into one canonical linear component representation. Do not merge it with normalization, scaling, HVG selection, Neighbors, integration, clustering, or metadata association. Each has a separately meaningful result and a different expression or independence contract.

The deletion test supports this seam. Without a PCA module, source/state validation, sparse solver policy, component bounds, overwrite semantics, output memory accounting, variance reporting, and code generation would reappear in every consumer. Conversely, merging PCA with NormalizeToLayer or HVG would force one preprocessing recipe on integration and residual workflows, while merging it with Neighbors would hide the choice of how many computed PCs a graph consumes.

PCA Metadata Associations remains a downstream read-only diagnostic module. It can be rerun against different metadata without recomputing PCA and must use Sample as its independent unit; folding it into PCA would couple an unsupervised representation to one study design and one family of hypotheses.

## Target interface

- Visible scientific inputs: `adata`, expression source (`X` or named current-axis layer; default suggestion `layer/log1p_norm`), `use_hvg`, and `n_comps`.
- Advanced state/resource inputs: `random_seed`, `overwrite_existing`, and a positive analyst-declared dense-output ceiling in GiB, without an arbitrary wrapper floor or ceiling.
- Outputs: copied `adata`, strict JSON `summary`, equivalent Python source text `code`.

Keep `zero_center=True`, `svd_solver="arpack"`, `dtype="float32"`, canonical output keys, `chunked=False`, `copy=False` for the internal Scanpy call, and `mask_var` routing as hidden policy. Exposing solver selection would make callers learn compatibility and numerical details without a demonstrated second production policy. `zero_center=False` is truncated SVD and should become an explicitly named method/node if needed rather than a PCA toggle. Chunked IncrementalPCA similarly deserves a deliberate scalability design because it ignores the seed/solver and densifies sparse chunks.

Do not expose a custom result key in this pass. The rest of the workflow uses the canonical `X_pca` contract. Multiple PCA variants remain possible through separate branches or explicit overwrite; custom Scanpy key conventions would expand every downstream interface.

## Validation and state contract

- Require an in-memory `AnnData` with at least two observations, at least two current-axis features, unique observation and variable names, and a dense/CSR/CSC finite numeric selected matrix. Reject backed or unsupported lazy matrices with an actionable `to_memory()` instruction.
- Resolve the source before validation. Never replace `X`, the source layer, another layer, or Raw snapshot.
- Accept any explicitly selected finite numeric source on which centered PCA is defined. Prefer provenance proving logarithmized, scaled, analytic-Pearson-residual, or another explicitly transformed expression state. Proven counts and size-normalized-but-untransformed values continue only with a prominent methodological warning; unknown state has a separate warning. The report must never upgrade an advisory state audit into programmatic proof of appropriate preprocessing.
- When `use_hvg=True`, require a boolean, non-missing `var["highly_variable"]`, at least two selected features, and apply it only to the current variable axis. When false, explicitly pass `mask_var=None` so an incidental external HVG column is not silently used by Scanpy defaults.
- Require `1 <= n_comps < min(n_obs, n_selected_vars)` for ARPACK. Count selected features with zero finite variance and report them; fail if every selected feature is constant.
- Preflight the lower-bound persistent dense output size using float32 cell scores plus float64 full-axis loadings and variance arrays. Fail before copying/computation when it exceeds the chosen ceiling and warn that solver workspace and the AnnData copy raise peak memory above this bound.
- Do not impose fixed UI maxima on `n_comps` or the memory budget. The positive scalar domains plus the exact ARPACK shape relation and the user's own budget are the complete execution contract.
- Treat `obsm["X_pca"]`, `varm["PCs"]`, and `uns["pca"]` as one bundle. If any exists, fail unless `overwrite_existing=True`; on permitted overwrite remove all three before calling Scanpy.
- After the call, validate score, loading, variance, and ratio shapes/dtypes/finiteness and the zero-loading rows for excluded features. Record canonical keys and expression provenance in namespaced history without altering Scanpy's own `uns["pca"]` contract.

## Report contract

`summary` includes:

- methods text naming centered PCA, explicit ARPACK solver, selected source/state, feature-mask policy, output dtype, and seed;
- results text stating the number of PCs and cumulative variance captured, without calling PCs clusters, corrected data, or Condition effects;
- key results for input cells/current genes/selected genes/constant genes, sparse input, requested/computed PC count, per-PC variance and variance ratio, cumulative ratios, leading absolute positive/negative loading identifiers, output shapes/keys, and dense persistent-output lower bound;
- resolved visible and hidden parameters, including separate JSON-native `expression_state` and `expression_state_evidence` fields, warnings for unknown external state or near-zero components, limitations about linearity, sign indeterminacy, PC-count choice, and cross-platform floating-point reproducibility;
- method/practice/software references and dynamic Python/openbio-singlecell/Scanpy/AnnData/scikit-learn/NumPy/SciPy versions.

`code` is a self-contained function with resolved arguments. It copies the input, routes the same source, performs the same provenance/mask/finite/ARPACK/overwrite/memory validations, calls `scanpy.pp.pca` with explicit `layer`, `mask_var`, `zero_center=True`, `svd_solver="arpack"`, `random_state`, `dtype="float32"`, `chunked=False`, and canonical keys, validates outputs, and returns an equivalent `AnnData`. It does not recreate plugin history.

## Cohesion and coupling assessment

The enhanced PCA node is a deep in-process module: a small scientific interface hides expression-state gating, current-axis feature selection, solver compatibility, sparse centering, shape/rank constraints, bundled output ownership, resource safeguards, variance/loadings reporting, citations, and version capture. A new adapter seam is not justified because only the public Scanpy implementation is supported.

## Verification plan

- Compare dense, CSR, and CSC `X`/layer inputs to direct Scanpy 1.12.3 ARPACK results for scores, loadings, variance, and variance ratios, allowing PCA sign equivalence where appropriate.
- Verify normalized-log, scaled, and Pearson-residual proven states pass; proven counts and untransformed normalized values run with the documented expert-choice warning; unknown external state emits its limitation.
- Cover boolean HVG selection, `use_hvg=False` in the presence of an incidental mask, non-boolean/missing/empty masks, excluded-feature zero loadings, and full-axis loading alignment.
- Exercise the exact ARPACK boundary (`n_comps=min_dim-1` succeeds; `min_dim` fails), one observation, too few features, all-constant/partly constant features, non-finite data, duplicate axes, backed data, and unsupported matrices.
- Verify the persistent dense-output estimate and a deliberately low ceiling fail before computation while ordinary sparse PCA does not materialize a cells-by-genes centered matrix.
- Verify default calls choose `layer/log1p_norm` in the UI, input/source/Raw remain unchanged, and canonical results are copy-on-write.
- Verify absent, partial, and complete pre-existing PCA bundles obey overwrite policy and never leave mixed old/new fields.
- Parse `summary` with strict JSON, compile `code`, and compare the generated function's primary object and failure modes with the node.
- Assert schema controls do not impose arbitrary component or memory maxima/floors, while direct calls still enforce positive values and the exact ARPACK/output-budget relations.
