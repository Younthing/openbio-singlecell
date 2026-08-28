# Scale: official usage research

Researched: 2026-08-28

## Official Scanpy interface and examples

The repository pins Scanpy `>=1.12.3,<1.13`; the installed version is 1.12.3. The public AnnData interface is:

```python
scanpy.pp.scale(
    adata,
    *,
    zero_center=True,
    max_value=None,
    copy=False,
    layer=None,
    obsm=None,
    mask_obs=None,
)
```

Scanpy defines the operation as scaling variables to unit variance and, by default, zero mean. Its current PBMC tutorial keeps transformed expression in layers and calls:

```python
adata.layers["scaled"] = adata.X.toarray()
scanpy.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], layer="scaled")
scanpy.pp.scale(adata, max_value=10, layer="scaled")
```

The legacy workflow uses `scanpy.pp.scale(adata, max_value=10)` after normalization, logarithmization, HVG selection, and optional regression.

- Official API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.scale.html
- Official PBMC workflow: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html
- Official PCA API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.pca.html

## Statistical semantics

Scanpy 1.12.3 computes each feature's variance and standard deviation with Bessel's correction (`ddof=1`, denominator `n_obs - 1`). This matters for very small datasets and slightly changes where clipping occurs. With `zero_center=True`, it subtracts the feature mean and divides by the standard deviation. With `zero_center=False`, it divides by the standard deviation without subtracting the mean; this is sparse-preserving variance scaling, not a zero-mean z-score.

Constant features are retained. Scanpy replaces a zero standard deviation with one; centered constant features become zero, while uncentered constant features retain their original values. This behavior is explicitly documented as subject to possible future change and must therefore be versioned and reported.

`max_value` is applied after scaling. The Scanpy 1.12.3 implementation clips centered output symmetrically to `[-max_value, max_value]`. When `zero_center=False`, it clips only the upper tail and logs a caution. `None` means no clipping. A UI must not use numeric zero as a sentinel for `None`, because `max_value=0` is a valid, materially different call.

For AnnData input, Scanpy writes the scaled matrix to `X` or the selected layer and adds feature-level means and standard deviations to `.var`. The wrapper must avoid confusing generic `.var["mean"]`/`.var["std"]` from multiple source layers; output-layer-specific provenance or namespaced statistics are safer.

## Sparse/dense and memory behavior

The official API states that `zero_center=False` permits efficient sparse handling. Zero-centering a CSR/CSC matrix emits a warning and densifies it. This can allocate at least `n_obs * n_vars * sizeof(float)` bytes, with a larger peak during copying and calculation. Integer input is cast to floating point. A node that copies the entire AnnData and then adds a dense layer can retain the sparse source and dense output simultaneously, so preflight must account for more than the final layer alone.

Explicit scaling is not required solely to obtain centered PCA from sparse data: Scanpy's default PCA solvers support implicit zero-centering without materializing a centered dense matrix. Therefore Scale should remain an explicit scientific choice rather than an unavoidable stage in every clustering workflow.

## Expression state and scientific practice

Common Scanpy workflows scale normalized/logarithmized expression after feature selection, optionally clipping extreme standardized values. Scaling raw counts directly does not replace library-size normalization or variance stabilization. Pearson residuals are already a centered, variance-stabilized representation under their own model, so applying Scale to them is a separate choice that can alter their relative variance and must be disclosed.

Scale is a feature transformation for representation learning or visualization. It does not remove a Technical batch, estimate a Condition effect, or produce inferential evidence. Whether to scale, whether to center, and whether to clip are method choices; thresholds such as ten standard deviations are common tutorial values, not universal biological constants.

## References

- Scanpy `pp.scale` official documentation: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.scale.html
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

## Required disclosure

`summary` should report the source expression state, dimensions and sparsity, centering choice, `ddof=1`, clipping mode and effective bound(s), input feature mean/standard-deviation distributions, constant-feature count, clipped-value count/rate, output range/dtype/density, sparse-to-dense transition and allocation estimate, destination, warnings/limitations, Scanpy/AnnData references, and Python/openbio-singlecell/Scanpy/AnnData/NumPy/SciPy versions. A report-ready result must distinguish centered z-scaling from uncentered variance scaling and must describe `max_value=10` as a user-selected truncation threshold rather than a biological standard.

`code` should provide an equivalent self-contained function that uses the same source, centering, clipping, storage, and memory rules; preserves the input expression state; and returns the same scaled destination layer without plugin-only history.

## Open-expert boundary audit

No scientific-preference gate is required beyond warnings already described above. The one-observation rejection is retained because Scanpy 1.12.3 uses `ddof=1` and cannot produce finite variance scaling; finite numeric/aligned input, destination collision checks, and the user-declared dense-memory ceiling are computation and safety boundaries. Expression-state suitability, clipping conventions, and whether scaling is useful remain warnings/limitations rather than hard validation.
