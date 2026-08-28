# Velocity Filter and Normalize: official usage and scientific practice

## Version audit and official interface

The repository does not currently declare or lock a scVelo optional dependency. Research was therefore performed against the latest signed stable release, scVelo 0.3.4 (2026-02-24); implementation must add and test one exact compatible optional stack rather than accept arbitrary installed versions.

In 0.3.4 the public interface is:

```python
scvelo.pp.filter_and_normalize(
    data, min_counts=None, min_counts_u=None, min_cells=None,
    min_cells_u=None, min_shared_counts=None, min_shared_cells=None,
    retain_genes=None, layers_normalize=None, copy=False, **kwargs,
)
```

It calls `filter_genes` and `normalize_per_cell`. Despite a stale docstring phrase, the 0.3.4 implementation no longer log-transforms and no longer performs highly-variable-gene selection; the 0.3.4 release removed `log1p` and `filter_genes_dispersion` from this recipe.

Primary sources:

- [scVelo 0.3.4 release](https://github.com/theislab/scvelo/releases/tag/v0.3.4)
- [scVelo `filter_and_normalize`](https://scvelo.readthedocs.io/en/stable/scvelo.pp.filter_and_normalize.html)
- [scVelo `filter_genes`](https://scvelo.readthedocs.io/en/stable/scvelo.pp.filter_genes.html)
- [scVelo getting started](https://scvelo.readthedocs.io/en/stable/getting_started.html)
- Bergen et al., *Nature Biotechnology* 2020, DOI [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)
- La Manno et al., *Nature* 2018, DOI [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6)

## P0 incompatibility and input contract

The current node passes `n_top_genes` and `subset_highly_variable` to 0.3.4. They are not parameters of `filter_and_normalize`; they flow through `**kwargs` to `normalize_per_cell`, which does not accept them. The current default execution therefore fails on the audited release.

RNA velocity requires observation- and feature-aligned, nonnegative raw spliced and unspliced molecule/count layers produced by an appropriate quantifier (for example velocyto). A layer name is not proof of raw count state. Reject missing, negative, nonfinite, non-integer-like, mismatched or already normalized layers. Spliced and unspliced must remain separate and may not alias the same object/key.

scVelo's normalization detection samples values and can infer state from dtype/value shape. A production wrapper must not rely on this heuristic. It must validate declared count provenance, call normalization with an explicit layer list and enforced policy, then verify library-size/post-filter results. Versioned 0.3.4 source inspection also shows that `normalize_per_cell` unconditionally prepends `X` to the requested layer list. Therefore `layers=["spliced", "unspliced"]` does **not** by itself mean “layers only”; a wrapper that promises not to redefine expression `X` must run on a private copy, snapshot and restore `X` and any `n_counts` collision, and prove that only the canonical velocity layers persistently changed.

Filtering changes the feature axis. It must occur on a dedicated velocity working state rather than destroy the repository's full-gene `Raw snapshot` or the expression object used for count-model inference. `X`, raw, embeddings, and expression provenance cannot be silently redefined as splicing kinetics.

## Scientific disclosure

This stage prepares velocity abundances; it does not estimate velocity, direction, pseudotime or fate. Report input/output cells/genes, layer identities, per-layer totals, removed genes and threshold reasons, normalization target/policy, zero-library cells, versions and warnings. The scVelo models assume splicing kinetics and are sensitive to quantification, sampling and biological steady/transient-state assumptions.
