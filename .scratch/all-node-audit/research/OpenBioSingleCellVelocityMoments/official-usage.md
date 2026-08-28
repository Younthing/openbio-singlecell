# Velocity Moments: official usage and scientific practice

## Official API

For scVelo 0.3.4:

```python
scvelo.pp.moments(
    data, n_neighbors=30, n_pcs=None, mode="connectivities",
    method="umap", use_rep=None, use_highly_variable=True, copy=False,
)
```

It computes first-order neighborhood moments of canonical `spliced` and `unspliced` layers and writes dense float32 `Ms` and `Mu`. It expects a valid neighbor graph. In 0.3.4 it can still automatically compute/recompute neighbors when the requested neighbor count exceeds the stored count, but that behavior emits a deprecation warning and is scheduled for removal; official direction is to compute neighbors first with Scanpy.

Sources:

- [scVelo `moments`](https://scvelo.readthedocs.io/en/stable/scvelo.pp.moments.html)
- [scVelo `neighbors`](https://scvelo.readthedocs.io/en/stable/scvelo.pp.neighbors.html)
- [scVelo 0.3.4 release](https://github.com/theislab/scvelo/releases/tag/v0.3.4)
- Bergen et al. 2020, DOI [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)
- Haghverdi et al. 2016, DOI [10.1038/nmeth.3971](https://doi.org/10.1038/nmeth.3971)

## Graph contract and current risks

Moments are defined by the exact neighbor weights, not merely `n_pcs`/`n_neighbors`. The current node supplies `None` for zeros and lets scVelo choose or reuse global defaults. It exposes no graph key and can therefore consume a stale/unrelated default graph or trigger hidden graph construction. That couples this atomic smoothing step to representation selection and graph building.

The supported wrapper must consume an explicitly named, already validated observation graph and a prepared velocity state with canonical normalized layers. It must not normalize raw layers, compute PCA, build neighbors or change the gene axis implicitly. Because scVelo 0.3.4 reads canonical neighbor storage, a private working copy may alias the requested named graph to canonical keys for the call; the public artifact retains the original name/fingerprint.

`Ms`/`Mu` can be dense and large. Estimate memory before execution and expose a guard. Verify both matrices are finite, nonnegative, observation/feature aligned and derived from the exact requested graph.

Moments are local smoothed abundances, not velocities or independent observations. Report graph identity, mode, dimensions, memory, layer distributions and versions.
