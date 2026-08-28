# Highly Variable Genes: official usage research

Researched: 2026-08-28

## Official Scanpy interfaces

The repository pins Scanpy `>=1.12.3,<1.13`; the installed version is 1.12.3. Stable and experimental HVG methods are exposed through two public functions:

```python
scanpy.pp.highly_variable_genes(
    adata,
    *,
    layer=None,
    n_top_genes=None,
    min_disp=0.5,
    max_disp=np.inf,
    min_mean=0.0125,
    max_mean=3,
    span=0.3,
    n_bins=20,
    flavor="seurat",
    subset=False,
    inplace=True,
    batch_key=None,
    filter_unexpressed_genes=None,
    check_values=True,
)

scanpy.experimental.pp.highly_variable_genes(
    adata,
    *,
    theta=100,
    clip=None,
    n_top_genes=None,
    batch_key=None,
    chunksize=1000,
    flavor="pearson_residuals",
    check_values=True,
    layer=None,
    subset=False,
    inplace=True,
)
```

- Stable API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html
- Pearson-residual API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.highly_variable_genes.html
- Official PBMC workflow: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html
- Official Pearson-residual tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html

## Flavor-specific input and ranking contract

The selected flavor determines both the method and the required expression state. A single generic “expression matrix” contract is not scientifically sufficient.

| Flavor | Public function | Required source | Method and ranking | Additional dependency |
| --- | --- | --- | --- | --- |
| `seurat` | `scanpy.pp.highly_variable_genes` | normalized, logarithmized expression | mean-binned normalized dispersion; with `n_top_genes`, cutoff parameters are ignored | none |
| `cell_ranger` | `scanpy.pp.highly_variable_genes` | normalized, logarithmized expression | Cell Ranger-style dispersion normalization/ranking | none |
| `seurat_v3` | `scanpy.pp.highly_variable_genes` | raw counts | regularized standard deviation and normalized variance; with batches, median rank first, number of batches as tie-breaker | `scikit-misc` |
| `seurat_v3_paper` | `scanpy.pp.highly_variable_genes` | raw counts | same within-batch statistic; with batches, number of batches first, median rank as tie-breaker, matching `SelectIntegrationFeatures` | `scikit-misc` |
| `pearson_residuals` | `scanpy.experimental.pp.highly_variable_genes` | raw UMI counts | genes ranked by variance of analytic Pearson residuals | none; experimental Scanpy interface |

The standard API explicitly says that `seurat` and `cell_ranger` expect logarithmized data, whereas the v3 flavors expect counts. The experimental API explicitly expects raw counts. These expectations must be disclosed prominently, but provenance metadata is not infallible and integer likeness is not a backend precondition: finite non-negative fractional count-like inputs may execute with warnings. Actual non-finite values, negative values for count-model flavors, zero totals that make the selected statistic undefined, and backend-incompatible axis/group sizes remain hard numeric or backend gates.

Official examples illustrate both branches:

```python
# Dispersion-based branch
scanpy.pp.normalize_total(adata, target_sum=1e4)
scanpy.pp.log1p(adata)
scanpy.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# Count-based v3 branch
scanpy.pp.highly_variable_genes(
    adata,
    layer="counts",
    flavor="seurat_v3",
    n_top_genes=2000,
)

# Count-based Pearson-residual branch
scanpy.experimental.pp.highly_variable_genes(
    adata,
    layer="counts",
    flavor="pearson_residuals",
    n_top_genes=2000,
)
```

For Pearson residuals, `theta=100` and `clip=None` (resolved to `sqrt(n_obs)`) are model parameters. `chunksize` is a performance parameter controlling how many genes are processed at once; smaller chunks reduce memory at the cost of speed. Unlike full residual normalization, this method need not materialize the entire dense residual matrix.

## Batch-aware selection

When `batch_key` is provided, Scanpy selects features within each group and merges the rankings. For dispersion flavors and Pearson residuals, genes are prioritized by the number of groups in which they are HVGs, with the documented flavor-specific statistic breaking ties. `seurat_v3` and `seurat_v3_paper` differ only in the ordering of median rank versus group recurrence when a batch key is present.

The output may include `highly_variable_nbatches` and `highly_variable_intersection`. These are selection diagnostics, not evidence that unwanted variation was removed. In this repository's domain language, `Sample` is a biological replicate and `Technical batch` is a separate covariate. The report must name the actual column used and must not relabel a Sample or Condition as a Technical batch. HVG selection is exploratory feature construction, not condition-level inference.

## Output fields and axis effects

With `inplace=True` and `subset=False`, Scanpy annotates `.var` with `highly_variable` and method-specific metrics:

- dispersion flavors: `means`, `dispersions`, `dispersions_norm`;
- v3 flavors: `means`, `variances`, `variances_norm`, `highly_variable_rank`;
- Pearson residuals: `means`, `variances`, `residual_variances`, `highly_variable_rank`;
- batch-aware runs: `highly_variable_nbatches`, `highly_variable_intersection` where supported.

`subset=True` then slices the current variable axis and every aligned current-axis matrix/layer. AnnData's official Raw behavior is different: variable slicing leaves `.raw` unchanged. Therefore a full-gene Raw snapshot can preserve the accepted post-QC counts while the working object is reduced, but subsetting without that snapshot loses full-gene recovery from current matrices.

That recovery concern is a best-practice warning, not a prerequisite for AnnData subsetting. An expert may deliberately subset without `.raw`; the node must disclose the loss of an in-object full-gene snapshot rather than veto the operation. Likewise, `n_top_genes` larger than available/expressed features is a request for “up to” that many genes: Scanpy may return all eligible genes or fewer and warn, so the wrapper reports requested versus actual counts instead of rejecting the request in advance.

- AnnData Raw: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.raw.html
- AnnData layers: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.layers.html

Stable HVG methods support dense, CSR, and CSC arrays. Pearson-residual HVG supports dense and SciPy sparse matrices and provides `chunksize` specifically to limit working memory. Backed objects and unsupported lazy arrays should not be assumed to work merely because AnnData can store them.

## Scientific practice and references

HVG selection is a dimension-reduction heuristic for downstream representation learning. The selected set depends on the expression state, flavor, number of requested genes, grouping key, and any manually forced genes. It is not a list of differentially expressed genes and does not constitute condition evidence.

- Satija R, Farrell JA, Gennert D, Schier AF, Regev A. Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology*. 2015;33:495-502. https://doi.org/10.1038/nbt.3192
- Zheng GXY et al. Massively parallel digital transcriptional profiling of single cells. *Nature Communications*. 2017;8:14049. https://doi.org/10.1038/ncomms14049
- Stuart T et al. Comprehensive Integration of Single-Cell Data. *Cell*. 2019;177(7):1888-1902.e21. https://doi.org/10.1016/j.cell.2019.05.031
- Lause J, Berens P, Kobak D. Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data. *Genome Biology*. 2021;22:258. https://doi.org/10.1186/s13059-021-02451-7
- Wolf FA, Angerer P, Theis FJ. SCANPY. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Virshup I et al. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

## Required disclosure

`summary` should report the flavor and applicable method reference, expression source/state and count validation, requested top-gene count, algorithm-selected count, manually forced count, final selected count, missing forced genes, input/output dimensions, subset behavior and Raw-snapshot availability, grouping column and group sizes, recurrence/intersection diagnostics, the applicable flavor-specific parameters and metric distributions, warnings/limitations, and runtime Python/openbio-singlecell/Scanpy/AnnData/NumPy/Pandas/SciPy plus `scikit-misc` when used. A report-ready result must distinguish selected features from biological markers and disclose any manual overrides that increased the final set beyond `n_top_genes`.

`code` should provide an equivalent self-contained function that dispatches to the same public Scanpy flavor, reproduces hard numeric/backend checks and expression-state warnings, applies the same grouping and forced-gene rules, and reproduces annotation or explicit subsetting without plugin-only history.
