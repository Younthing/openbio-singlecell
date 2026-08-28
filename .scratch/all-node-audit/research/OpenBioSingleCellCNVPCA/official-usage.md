# CNV PCA: official usage and scientific practice

## Audited API

The replacement is audited against infercnvpy 0.6.1. Its public API is:

```python
infercnvpy.tl.pca(
    adata,
    svd_solver="arpack",
    zero_center=False,
    inplace=True,
    use_rep="cnv",
    key_added="cnv_pca",
    **kwargs,
)
```

The wrapper requires `obsm[f"X_{use_rep}"]`, calls Scanpy PCA on that matrix, and writes
`obsm[f"X_{key_added}"]`. Its defaults deliberately use `zero_center=False`, which avoids densifying a sparse CNV
window matrix. The result is a lower-dimensional exploratory representation; it is not an expression PCA, a copy
number estimate, a tumor classifier, or an inferential test.

Primary sources:

- [infercnvpy `tl.pca`](https://infercnvpy.readthedocs.io/en/stable/generated/infercnvpy.tl.pca.html)
- [infercnvpy 0.6.1 source](https://github.com/icbi-lab/infercnvpy/blob/v0.6.1/src/infercnvpy/tl/_pca.py)
- [Scanpy PCA](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.pca.html)
- Tirosh et al., *Science* 2016, DOI [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)
- Lehoucq, Sorensen and Yang, *ARPACK Users' Guide*, SIAM 1998, DOI
  [10.1137/1.9780898719628](https://doi.org/10.1137/1.9780898719628)

## Scientific and output contract

The input must be a validated `OPENBIO_CNV_STATE` at stage `inferred`; a coincidentally named `X_cnv` matrix is
insufficient. Validate its observation axis, finite real sparse/dense CNV matrix, window metadata, source/reference
provenance, and complete artifact fingerprint before PCA. Require `n_comps < min(n_obs, n_windows)` for ARPACK,
guard output memory and collisions, and use an explicit seed.

Call infercnvpy with `svd_solver="arpack"`, `zero_center=False`, `inplace=True`, explicit `use_rep`, `key_added`,
`n_comps`, `random_state`, and floating output dtype. Verify an observation-aligned finite nonconstant score matrix
with exact component count, unchanged input artifact, and repeatability. Report variance of every component
independently because the thin wrapper does not store a Scanpy PCA metadata bundle for an array input. Component
signs are arbitrary and the representation remains conditional on the upstream reference, genomic windows and
smoothing choices.
