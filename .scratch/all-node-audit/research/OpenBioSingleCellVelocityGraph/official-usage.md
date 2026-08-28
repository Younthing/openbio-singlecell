# Velocity Graph: official usage and scientific practice

## Official API

For scVelo 0.3.4:

```python
scvelo.tl.velocity_graph(
    data, vkey="velocity", xkey="Ms", tkey=None, basis=None,
    n_neighbors=None, n_recurse_neighbors=None,
    random_neighbors_at_max=None, sqrt_transform=None,
    variance_stabilization=None, gene_subset=None,
    compute_uncertainties=None, approx=None,
    mode_neighbors="distances", copy=False, n_jobs=None,
    backend="loky", show_progress_bar=True,
)
```

It compares each velocity vector with candidate state changes using cosine correlation and writes sparse positive and negative directed graphs under `uns[f"{vkey}_graph"]` and `uns[f"{vkey}_graph_neg"]`, plus parameters. It depends on the velocity layer, the expression/moment layer selected by `xkey`, the velocity-gene mask and the neighbor graph.

The 0.3.4 constructor resolves `n_recurse_neighbors=None` to 2 for distance-mode graphs when `n_neighbors` is also `None`, so the apparent public default expands candidates to neighbors-of-neighbors. A wrapper claiming exact direct named-graph support must pass `n_recurse_neighbors=1` explicitly and verify the backend records 1. It should also disclose a fixed `sqrt_transform=False`; leaving that argument `None` would automatically transform stochastic velocities and make graph policy depend on the upstream model.

Sources:

- [scVelo `velocity_graph`](https://scvelo.readthedocs.io/en/stable/scvelo.tl.velocity_graph.html)
- [scVelo getting started](https://scvelo.readthedocs.io/en/stable/getting_started.html)
- Bergen et al. 2020, DOI [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)
- La Manno et al. 2018, DOI [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6)

## Required graph artifact boundary

The current node exposes only `n_jobs`; it does not establish which velocity/moment/neighbor state is being used or validate the returned directed matrices. A named sparse matrix alone is insufficient: observation order, vkey, xkey, source neighbor graph, gene subset and model provenance determine its meaning.

The result is directed transition **evidence** from cosine alignment, not a row-stochastic transition-probability matrix. CellRank's `VelocityKernel.compute_transition_matrix` performs a separate model/normalization step, and GPCCA is a separate estimator. These must not be smuggled into this node or represented by the same artifact.

Report matrix shape/nnz, positive/negative edge distributions, zero-outdegree cells, model and source graph identities, compute policy, and versions. Direction is kinetic-model dependent and does not prove fate or lineage.
