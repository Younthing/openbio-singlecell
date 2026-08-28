# CNV Structure: official usage and current composite audit

## Official infercnvpy operations

The repository does not lock infercnvpy. Audit against official infercnvpy **0.6.1** and require an exact compatible optional dependency smoke before implementing any replacement. The current node calls five distinct public operations with hidden defaults:

```python
infercnvpy.tl.pca(adata, use_rep="cnv", key_added="cnv_pca", ...)
infercnvpy.pp.neighbors(adata, use_rep="cnv_pca", key_added="cnv_neighbors", ...)
infercnvpy.tl.leiden(adata, neighbors_key="cnv_neighbors", key_added="cnv_leiden", ...)
infercnvpy.tl.umap(adata, neighbors_key="cnv_neighbors", key_added="cnv_umap", ...)
infercnvpy.tl.cnv_score(adata, groupby="cnv_leiden", use_rep="cnv", key_added="cnv_score", ...)
```

Sources:

- [infercnvpy API](https://infercnvpy.readthedocs.io/en/stable/api.html)
- [infercnvpy PCA](https://infercnvpy.readthedocs.io/en/stable/generated/infercnvpy.tl.pca.html)
- [infercnvpy neighbors](https://infercnvpy.readthedocs.io/en/latest/generated/infercnvpy.pp.neighbors.html)
- [infercnvpy Leiden](https://infercnvpy.readthedocs.io/en/latest/generated/infercnvpy.tl.leiden.html)
- [infercnvpy UMAP](https://infercnvpy.readthedocs.io/en/latest/generated/infercnvpy.tl.umap.html)
- [infercnvpy CNV score](https://infercnvpy.readthedocs.io/en/stable/generated/infercnvpy.tl.cnv_score.html)
- [infercnvpy 0.6.1 release](https://github.com/icbi-lab/infercnvpy/releases/tag/v0.6.1)
- Tirosh et al., *Science* 2016, DOI [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)
- McInnes, Healy and Melville, *arXiv* 2018 (UMAP), DOI [10.48550/arXiv.1802.03426](https://doi.org/10.48550/arXiv.1802.03426)
- Traag, Waltman and van Eck, *Scientific Reports* 2019 (Leiden), DOI [10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)

## Correctness/cohesion findings

These are not one atomic analysis. PCA chooses representation dimensionality; neighbors chooses graph geometry; Leiden chooses a resolution/stochastic partition; UMAP is a visualization; `cnv_score` assigns the **group mean absolute inferred-CNV value** to every member. Every important parameter is currently hidden, there is no seed/collision policy, and the input only assumes default `X_cnv` keys without verifying an InferCNV state.

The resulting score is cluster-dependent and descriptive. It is not an absolute copy number, tumor probability, calibrated cutoff or independent per-cell measurement. Clustering and UMAP can change while the upstream CNV matrix remains identical. Official tutorials inspect CNV clusters/scores and then make a domain-informed tumor/normal assignment; the package does not provide an automatic universal classification threshold.

This composite cannot produce an honest equivalent atomic function or report which choice drove which result. It should not be retained as an executing recipe.

## Public compatibility status

The stable composite ID remains solely as a non-executing migration surface. Its public schema must explicitly mark
deprecation and identify the atomic CNV PCA, graph, embedding, grouping, and score replacement path. Because the
shim performs none of those official operations, it has no legitimate scientific `summary` or equivalent `code`;
each active replacement reports its own method, references, software, and result.
