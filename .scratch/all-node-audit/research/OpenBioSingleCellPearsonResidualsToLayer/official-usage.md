# Pearson Residuals to Layer: official usage research

Researched: 2026-08-28

## Official Scanpy interface and example

The repository pins Scanpy `>=1.12.3,<1.13`; the installed version is 1.12.3. Its public experimental interface is:

```python
scanpy.experimental.pp.normalize_pearson_residuals(
    adata,
    *,
    theta=100,
    clip=None,
    check_values=True,
    layer=None,
    obsm=None,
    inplace=True,
    copy=False,
)
```

The function expects raw UMI counts. With `inplace=False`, it returns a dictionary whose `"X"` entry is the residual matrix and whose other fields disclose `theta`, `clip`, and the source on which the residuals were computed. This maps directly to a non-destructive layer workflow:

```python
result = scanpy.experimental.pp.normalize_pearson_residuals(
    adata,
    layer="counts",
    theta=100,
    clip=None,
    check_values=True,
    inplace=False,
)
adata.layers["analytic_pearson_residuals"] = result["X"]
```

- Official API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.normalize_pearson_residuals.html
- Official tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html
- Experimental-module status: https://scanpy.readthedocs.io/en/stable/api/experimental.html

## Method and parameter semantics

For count matrix \(x_{ij}\), the method fits a negative-binomial offset null model with expected value determined by the cell and gene totals. It reports the Pearson residual

\[
r_{ij} = \frac{x_{ij} - \mu_{ij}}{\sqrt{\mu_{ij} + \mu_{ij}^{2}/\theta}}.
\]

`theta` is a gene-shared overdispersion parameter: smaller positive values imply greater overdispersion, `theta=100` is the Scanpy and Lause et al. recommendation, and `theta=np.inf` gives the Poisson limit. `theta` must be positive.

Clipping is part of the scientific method, not merely a performance option:

- `clip=None` means the default symmetric interval `[-sqrt(n_obs), sqrt(n_obs)]`;
- a finite scalar `c >= 0` means `[-c, c]`;
- `clip=np.inf` means no clipping.

Consequently, a UI value of `0` must not be overloaded to mean “automatic”: zero is a valid Scanpy threshold that makes every residual zero. The interface needs an explicit automatic/custom/no-clipping choice and must report both the requested mode and the resolved numeric bound.

Analytic Pearson residuals replace library-size normalization and log transformation for this branch of a workflow. They should not be computed on a matrix that has already been normalized, logged, scaled, or batch-corrected.

## Sparse, dense, layer, and Raw behavior

The official tutorial warns that sparse UMI counts become a dense residual matrix. It recommends selecting highly variable genes from counts first, using the chunked Pearson-residual HVG method, and computing the full residual matrix only on the reduced gene set. The output therefore has an unavoidable allocation of at least `n_obs * n_vars * sizeof(float)` bytes, with a higher peak while the expected-value and intermediate matrices coexist.

`AnnData.layers` entries have the same observation-by-variable shape as `X`, so storing the residuals in a named layer preserves the selected count source and keeps the expression state explicit. A Raw snapshot may contain a different, full-gene variable axis after feature selection; it is therefore not a suitable direct source for a current-axis residual layer. The node should support current `X` or a named layer, not `raw`.

- AnnData layers: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.layers.html
- AnnData Raw slicing behavior: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.raw.html

Scanpy 1.12.3 checks integer likeness only by warning. Its implementation divides by the cell total and model standard deviation; an all-zero dataset, a zero-total cell, or a zero-total gene yields undefined values and remains a hard precondition. Non-negativity is part of the declared count-model definition and is also hard. Fractional non-negative values, however, are numerically executable: they proceed with a count-assumption warning, and provenance suggesting a transformed state is likewise disclosed rather than trusted as an infallible veto. Singleton axes are accepted when the backend returns a finite residual matrix; sample-variance summaries are then unavailable and must be null/disclosed rather than blocking the transformation.

## Scientific practice and references

The method was developed and benchmarked for single-cell RNA-seq UMI counts. The report must call the result a residual representation, not corrected counts and not condition-level inference. It is an experimental Scanpy interface, so its runtime version is important provenance.

- Lause J, Berens P, Kobak D. Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data. *Genome Biology*. 2021;22:258. https://doi.org/10.1186/s13059-021-02451-7
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

## Required disclosure

`summary` should report the count source, input dimensions and sparsity, total counts, zero-total validation, `theta`, clipping mode and resolved bound, output layer, output dtype/density and memory estimate, residual distribution, clipped-value count/rate, per-gene residual-variance distribution, warnings/limitations, the Lause/Scanpy/AnnData references, and Python/openbio-singlecell/Scanpy/AnnData/NumPy/SciPy versions. A report-ready result should state that analytic Pearson residuals were computed from raw UMI counts under a shared-theta negative-binomial offset model and should not imply biological correction or inferential significance.

`code` should provide an equivalent self-contained function that resolves the same count source and clipping rule, reproduces the same hard numeric/model checks and expert warnings, preserves the source state, finite-postvalidates the backend result, and writes the residual matrix to the disclosed destination layer.
