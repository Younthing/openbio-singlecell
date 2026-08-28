# t-SNE: official usage and scientific practice

## Audited runtime and primary sources

The environment contains Scanpy 1.12.3 and scikit-learn 1.9.0. The installed Scanpy signature is:

```python
scanpy.tl.tsne(
    adata, n_pcs=None, *, n_components=2, use_rep=None,
    perplexity=30, metric="euclidean", early_exaggeration=12,
    learning_rate=1000, random_state=0, use_fast_tsne=False,
    n_jobs=None, key_added=None, copy=False,
)
```

Primary references:

- [Scanpy `tl.tsne` documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.tsne.html)
- [scikit-learn 1.9 `TSNE` documentation](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)
- [van der Maaten and Hinton, “Visualizing Data using t-SNE,” JMLR 9:2579–2605](https://www.jmlr.org/beta/papers/v9/vandermaaten08a.html)
- [Scanpy paper, DOI 10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
- [scikit-learn paper, JMLR 12:2825–2830](https://jmlr.org/papers/v12/pedregosa11a.html)

Scanpy 1.12.3 uses scikit-learn by default. `use_fast_tsne=True` is deprecated and MulticoreTSNE support is slated for removal, so the supported wrapper should keep it false.

## Representation resolution and current P0

Scanpy t-SNE uses the same `_choose_representation` logic as Neighbors. With an explicit `.obsm` `use_rep` and integer `n_pcs`, it slices the first `n_pcs` dimensions and raises if that representation is narrower. With `n_pcs=None`, it uses all explicit `.obsm` dimensions.

The pre-refactor OpenBio node defaulted to `use_rep="X_pca"` and always passed `n_pcs=50`, which made narrower learned latent representations fail. The refactored interface uses representation-neutral “dimensions”; `n_dimensions=0` maps to `n_pcs=None` and therefore uses all dimensions of an explicit `.obsm` representation, while positive values are shape-checked.

If `use_rep=None`, Scanpy may choose X or silently calculate PCA depending on feature count. That hidden data-dependent branch is undesirable in a reproducible workflow wrapper; the representation should be explicit.

## Parameter behavior

- `perplexity` controls the effective neighborhood scale in the high-dimensional affinity distribution. scikit-learn requires it to be strictly less than the number of observations; 5–50 is a common exploration range, not a universal optimum.
- `metric` defines geometry in the selected representation. Scanpy defaults to Euclidean.
- `early_exaggeration=12` affects initial cluster spacing.
- Scanpy fixes `learning_rate=1000`, whereas current scikit-learn defaults to `"auto"`; the executed value must be attributed to Scanpy and stated explicitly.
- `random_state=0` makes the stochastic non-convex optimization repeatable within a fixed environment. Different initializations or versions can produce different maps.
- `n_components=2` is appropriate for this visualization node.
- Scanpy passes `n_jobs` to scikit-learn; `None` uses Scanpy settings. Thread count is a performance control, not a scientific workflow parameter here.

With `key_added=None`, Scanpy stores `obsm["X_tsne"]` and `uns["tsne"]`. With a custom key it uses the exact custom key for both. Existing values are silently overwritten upstream.

## Preconditions and failure modes

- Explicit representation key must exist, align to `n_obs`, have at least one dimension, and contain finite numeric values.
- Requested dimensions must not exceed representation width.
- `0 < perplexity < n_obs`; very small datasets or a default perplexity of 30 need a clear preflight error.
- Unique observation identifiers are required for safe downstream joins.
- Constant/duplicate rows and inappropriate metrics can lead to degenerate affinities or unstable layouts and should be disclosed.
- `early_exaggeration` and `learning_rate` must be finite and positive. scikit-learn warns that too-high learning rates can create a ball-like map and too-low rates can compress points into a dense cloud.
- Output keys need collision policy.
- Backed data and non-finite output coordinates should fail explicitly.

## General scientific practice

scikit-learn recommends reducing very high-dimensional data to a moderate dimension (for example about 50) before t-SNE for speed and noise suppression. That does not justify demanding 50 dimensions from a learned 10D scVI space: the learned latent space is already a model-derived low-dimensional representation.

t-SNE is an exploratory visualization. Its objective is non-convex, axes are arbitrary, global distances and apparent cluster sizes are not directly interpretable, and visually separated islands do not establish cell types or `Condition` effects. Results should be checked across plausible perplexities and seeds. Clustering should use a high-dimensional graph/representation, not the two-dimensional t-SNE coordinates.

## Reporting expectations

Report the representation key, available/used dimensions, metric, perplexity, early exaggeration, learning rate, seed, backend, resolved output/parameter keys, coordinate shape and ranges, duplicate/degeneracy diagnostics, and software versions. State that no neighbor graph, `Technical batch` correction, clustering, `Sample` aggregation, or replicate-aware `Condition` inference was performed.

## Open expert-tool boundary

The backend requirement is `0 < perplexity < n_obs`; sub-unit positive perplexities are executable and must not be excluded by the workflow schema. Extreme but valid perplexity-to-cell ratios, repeated rows, and degenerate layouts are warnings/diagnostics. An explicit `.obsm` representation does not require any gene variables. Representation alignment, real finite values, metric-specific undefined rows, and finite two-dimensional output remain hard.
