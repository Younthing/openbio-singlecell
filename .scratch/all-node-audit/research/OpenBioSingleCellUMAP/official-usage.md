# UMAP: official usage and scientific practice

## Audited runtime and primary sources

The installed runtime is Scanpy 1.12.3 with umap-learn 0.5.12, scikit-learn 1.9.0, AnnData 0.13.2, NumPy 2.4.4, and SciPy 1.17.1. The installed Scanpy signature is:

```python
scanpy.tl.umap(
    adata,
    *, min_dist=0.5, spread=1.0, n_components=2,
    maxiter=None, alpha=1.0, gamma=1.0, negative_sample_rate=5,
    init_pos="spectral", random_state=0, a=None, b=None,
    method="umap", key_added=None, neighbors_key="neighbors", copy=False,
)
```

Primary references:

- [Scanpy `tl.umap` documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html)
- [umap-learn API guide](https://umap-learn.readthedocs.io/en/latest/api.html)
- [umap-learn parameter guide](https://umap-learn.readthedocs.io/en/latest/parameters.html)
- [umap-learn reproducibility guide](https://umap-learn.readthedocs.io/en/latest/reproducibility.html)
- [UMAP software paper, DOI 10.21105/joss.00861](https://doi.org/10.21105/joss.00861)
- [UMAP method preprint, DOI 10.48550/arXiv.1802.03426](https://doi.org/10.48550/arxiv.1802.03426)

## Official operation and storage

Scanpy UMAP embeds an **existing neighborhood graph**; it does not choose the high-dimensional representation or recompute neighbors. `neighbors_key` points to `adata.uns[neighbors_key]`, whose `connectivities_key` points to the actual sparse graph in `.obsp`.

With `key_added=None`, results are stored in `obsm["X_umap"]` and `uns["umap"]`. With a custom `key_added`, both the coordinate and parameter records use that exact custom key. Scanpy silently replaces existing values, so a wrapper needs explicit collision policy.

The current OpenBio node correctly keeps `n_components=2`, `maxiter=None`, `alpha=1`, `gamma=1`, `negative_sample_rate=5`, `init_pos="spectral"`, `a/b=None`, and CPU `method="umap"` hidden. It exposes Scanpy's `min_dist=0.5` and `spread=1.0`; note that bare umap-learn defaults `min_dist` to 0.1, so reports must attribute the executed default to Scanpy.

## Parameter interpretation

- `min_dist` is the effective minimum separation of nearby embedded points. Smaller values make local groups visually tighter; it is not a biological cluster-separation parameter.
- `spread` sets the embedding's overall effective scale together with `min_dist`. Values are meaningful only as a pair; `spread` must be positive and unusual but executable `min_dist > spread` configurations should be clearly warned.
- `init_pos="spectral"` initializes from the graph spectrum. `random`, an existing `.obsm` key, and PAGA initialization are also supported upstream, but each adds a prerequisite contract.
- `random_state` controls stochastic optimization. umap-learn documents that fixed seeds improve exact reproducibility but disable non-deterministic parallel optimization and can reduce performance.
- `a` and `b` are normally derived from `min_dist` and `spread`; exposing them alongside those higher-level controls would create a contradictory shallow interface.
- GPU `method="rapids"` in Scanpy is deprecated in favor of rapids-singlecell and is a different dependency/backend contract.

## Preconditions and failure modes

- `uns[neighbors_key]` must exist and be a mapping with a valid `connectivities_key`.
- The referenced `.obsp` matrix must exist, be square `n_obs × n_obs`, numeric, finite, nonnegative, and contain meaningful edges.
- Observation identifiers should be unique.
- Disconnected graphs are allowed by UMAP but can create separately placed components and warnings; component count and isolated observations materially affect interpretation.
- `min_dist` must be finite and nonnegative; `spread` must be finite and positive.
- Existing output keys must not be overwritten silently.
- PAGA initialization requires prior PAGA calculation and positions; an arbitrary `.obsm` initializer must be finite, aligned, and two-dimensional.
- A fixed seed cannot make results portable across different input graphs, package versions, or backends.

## General scientific practice

UMAP is chiefly a nonlinear exploratory visualization. Local neighbor relationships reflect the supplied graph, while apparent gaps, density, cluster area, and long-range separation can change with graph construction, seed, and UMAP parameters. A two-dimensional UMAP must not be used as the input for clustering when the high-dimensional neighbor graph is available, and visual separation is not evidence for a `Condition` effect.

For reproducibility, publish the graph provenance (`use_rep`, dimensions, metric, and neighbor count from the Neighbors summary), UMAP parameters, seed, output key, and software versions. Comparing several reasonable seeds/parameters is a robustness check; selecting a visually pleasing map is not statistical validation.

## Reporting expectations

The report should disclose the resolved graph and connectivity key, graph component diagnostics, `min_dist`, `spread`, hidden optimization defaults, initialization, seed, output key, coordinate shape/ranges and finite status, and overwrite status. It should explicitly state that UMAP did not change the graph, integrate a `Technical batch`, discover clusters, or perform `Sample`-level `Condition` inference.

## Open expert-tool boundary

Finite `min_dist >= 0` and `spread > 0` are hard numeric requirements. `min_dist > spread` is unusual but executable and is therefore a warning, not an error; values above one are also valid expert inputs. UMAP can consume finite nonnegative directed/asymmetric connectivity matrices and matrices with self-loops. Their symmetry and diagonal are reported as graph diagnostics instead of being silently changed or rejected. A gene axis is irrelevant because UMAP consumes only the named observation graph.
