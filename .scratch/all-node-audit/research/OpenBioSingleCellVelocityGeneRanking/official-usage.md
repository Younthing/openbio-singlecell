# RNA Velocity Gene Ranking: official usage and scientific practice

## What the current node actually computes

The current implementation does not call scVelo's `rank_velocity_genes`. It sorts `adata.var["fit_likelihood"]` from `recover_dynamics` and returns the first rows. This is a global dynamical-model fit-quality ranking. The display name “RNA Velocity Gene Ranking” is therefore ambiguous and can be mistaken for genes characterizing a cell group or for statistical differential evidence.

Related scVelo 0.3.4 public interfaces are:

```python
scvelo.tl.rank_velocity_genes(data, vkey="velocity", n_genes=100, groupby=None, ...)
scvelo.tl.rank_dynamical_genes(data, n_genes=100, groupby=None, copy=False)
```

Those routines answer group-specific questions and require a grouping key. They are not equivalent to sorting global `fit_likelihood`.

Sources:

- [scVelo API: dynamical and velocity genes](https://scvelo.readthedocs.io/en/stable/api.html)
- [scVelo dynamical modeling tutorial](https://scvelo.readthedocs.io/en/latest/DynamicalModeling.html)
- Bergen et al. 2020, DOI [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)

## Scientific contract

`fit_likelihood` is meaningful only for a verified `recover_dynamics` result. Require exact feature alignment, finite likelihood for eligible genes and explicit handling of genes that failed fitting. Stable ties must follow the fitted feature order. A `top_n` truncation is a reporting choice, not a significance threshold.

The table should name the metric and include fit status and relevant diagnostics. It must not call genes “differential,” “markers,” “drivers,” or “significant,” and it must not be used as `Condition` evidence.
