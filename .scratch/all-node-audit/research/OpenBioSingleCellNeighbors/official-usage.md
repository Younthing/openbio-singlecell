# Neighbors: official usage and scientific practice

## Audited runtime

The repository currently calls `scanpy.pp.neighbors` from Scanpy 1.12.3. The inspected environment contains AnnData 0.13.2, NumPy 2.4.4, SciPy 1.17.1, scikit-learn 1.9.0, umap-learn 0.5.12, and igraph 1.0.0.

The installed signature is:

```python
scanpy.pp.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=None,
    *,
    distances=None,
    use_rep=None,
    knn=True,
    method="umap",
    transformer=None,
    metric=None,
    metric_kwds={},
    random_state=0,
    key_added=None,
    copy=False,
)
```

Primary references:

- [Scanpy `pp.neighbors` documentation](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html)
- [Scanpy 1.12.3 implementation](https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/neighbors/__init__.py)
- [Scanpy paper, DOI 10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
- [UMAP reference implementation paper, DOI 10.21105/joss.00861](https://doi.org/10.21105/joss.00861)

## Official operation

Scanpy performs two conceptually separate stages: a nearest-neighbor search creates distances, then `method` converts those results into connectivities. `transformer` chooses the search implementation; `method` chooses the connectivity kernel. They must not be reported as one interchangeable “neighbor method.”

With default storage, Scanpy writes:

- `adata.uns["neighbors"]`, containing parameter and key metadata;
- `adata.obsp["distances"]`;
- `adata.obsp["connectivities"]`.

With `key_added="latent"`, it instead writes `uns["latent"]`, `obsp["latent_distances"]`, and `obsp["latent_connectivities"]`. Consumers such as UMAP and Leiden receive the **uns key** (`neighbors_key="latent"`), not the obsp key.

`n_neighbors=15`, `knn=True`, the UMAP connectivity kernel, and Euclidean metric are Scanpy defaults. The refactored OpenBio node also defaults to Euclidean while exposing metric as an explicit scientific geometry choice. Scanpy describes roughly 2–100 neighbors as a general practical range: smaller values emphasize local neighborhoods and larger values emphasize broader structure.

## Representation resolution and the dimensionality P0

Scanpy's installed `_choose_representation` has exact, non-interchangeable branches:

1. explicit `use_rep` matching `.obsm` plus integer `n_pcs`: slice `obsm[use_rep][:, :n_pcs]`;
2. explicit `.obsm` plus `n_pcs=None`: use **all available dimensions**;
3. `use_rep="X"`: use X; `n_pcs` does not slice X;
4. `use_rep=None`: use X for a small feature axis, otherwise use or silently calculate PCA.

If `n_pcs` exceeds an explicit `.obsm` representation's width, Scanpy raises `ValueError`. This is not a harmless upper bound.

The pre-refactor OpenBio interface always passed `n_pcs=50`, which made common 10-dimensional scVI representations fail. The refactored representation-neutral `n_dimensions=0` mode passes `n_pcs=None` and uses every available dimension; a positive explicit value remains shape-checked against the selected representation width.

For learned latent representations such as scVI or Harmony, general practice is to use all model-produced dimensions unless the model contract or analyst explicitly selects a subset. Calling every latent dimension a “PC” is also inaccurate.

## Search, kernel, metric, and randomness

- `knn=True` yields a hard k-nearest-neighbor graph. `knn=False` uses a Gaussian-width interpretation and can require much more memory; Scanpy warns for at least 10,000 observations.
- `method="umap"` uses UMAP fuzzy-simplicial-set connectivities. `gauss` uses an adaptive Gaussian kernel; `jaccard` uses a shared-neighbor/Jaccard kernel.
- `transformer=None` lets Scanpy choose exact search for small data and PyNNDescent for larger data. A seed matters for approximate search; exact search may not consume it.
- `metric` defines geometry in the selected representation. Euclidean is the official default. Cosine can be useful for direction-based latent geometry but is not universally preferable. Correlation is undefined or unstable for constant rows and all metrics require finite numeric input.
- A fixed `random_state` is necessary for reproducible approximate search, but software/backend versions and the selected representation remain part of reproducibility.

The node should not expose arbitrary metric callables or `metric_kwds` through a workflow JSON interface. A bounded documented metric set is serializable and testable.

## Preconditions and failure modes

- The selected representation must exist, be two-dimensional, align exactly to `n_obs`, contain at least one dimension, and contain only finite numeric values.
- Observation identifiers should be unique because downstream graph-derived annotations are aligned by observation identity.
- `n_neighbors` must be at least 2. Values below `n_obs` are usual practice, but requests at or above `n_obs` are executable; Scanpy may resolve an oversized request to `1 + floor(n_obs / 2)`. The wrapper therefore runs the expert choice and explicitly reports requested and effective values rather than failing or claiming the request was executed unchanged.
- Requested dimensions must be positive and no larger than the selected representation width; “all dimensions” should map to `n_pcs=None`.
- Existing `uns`/`obsp` output keys are silently overwritten by Scanpy. Collision policy therefore belongs in the wrapper.
- The returned distance graph is directed/sparse while connectivities are normally symmetric; disconnected components and zero-degree observations require disclosure because they affect embeddings and clustering.
- Backed/view inputs and malformed pre-existing graph metadata should be rejected before calculation.

## Reporting expectations

A reproducible report must state the declared and resolved representation, available and actually used dimensions, metric, requested/effective neighbor count, search transformer policy, connectivity kernel, `knn`, seed, and exact output keys. Graph diagnostics should include distance/connectivity nonzero counts, degree/weighted-degree summaries, connected-component count and sizes, isolated observations, symmetry checks, and whether existing outputs were replaced.

The graph is an exploratory representation of cells, not a `Sample`-level result, `Technical batch` correction, cluster label, or replicate-aware `Condition` inference.

## Open expert-tool boundary

`n_neighbors >= 2` remains the documented lower bound. A request equal to `n_obs` is executable, and a larger request is resolved by Scanpy to an effective neighborhood size; these expert choices must run with requested/effective values and a warning rather than fail. The selected representation may live entirely in `.obsm`, so a nonempty gene axis is not a prerequisite unless `use_rep="X"`. Fixed UI ceilings on dimensions or neighbors are not scientific validity checks; actual representation width and explicit resource failures are the hard bounds.
