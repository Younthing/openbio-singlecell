# PCA: official usage research

Researched: 2026-08-28

## Repository and installed software context

The repository pins Scanpy `>=1.12.3,<1.13`; the audited environment contains Scanpy 1.12.3, scikit-learn 1.9.0, SciPy 1.17.1, AnnData 0.13.2, NumPy 2.4.4, and pandas 3.0.5. The public function used by this node is:

```python
scanpy.pp.pca(
    data,
    n_comps=None,
    *,
    layer=None,
    obsm=None,
    zero_center=True,
    svd_solver=None,
    chunked=False,
    chunk_size=None,
    random_state=0,
    return_info=False,
    mask_var=_empty,
    use_highly_variable=None,
    dtype="float32",
    key_added=None,
    copy=False,
)
```

- Scanpy 1.12 PCA API: https://scanpy.readthedocs.io/en/1.12.x/generated/scanpy.pp.pca.html
- Current Scanpy PCA API and storage contract: https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pp.pca.html
- scikit-learn PCA API: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
- SciPy ARPACK/SVD interface: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html

Scanpy's official example is the direct in-place call:

```python
scanpy.pp.pca(
    adata,
    n_comps=50,
    layer="log1p_norm",
    mask_var="highly_variable",
    zero_center=True,
    svd_solver="arpack",
    random_state=0,
    dtype="float32",
)
```

The repository's packaged clustering workflow follows the same scientific order: construct a normalized/logarithmized `log1p_norm` layer, identify highly variable genes, and run PCA from that layer without first materializing a centered dense scaled matrix.

## Expression-state and feature-mask contract

PCA is a linear variance decomposition. It does not normalize library size, log-transform counts, select features, or remove a Technical batch. Scanpy centers features for ordinary PCA but does not scale each feature to unit variance. Therefore the caller must deliberately choose a meaningful expression state.

For the standard single-cell clustering branch, normalized/logarithmized expression over highly variable genes is the safe default. Centered PCA is also meaningful on an explicitly scaled layer or on analytic Pearson residuals, although those are distinct preprocessing choices and must be disclosed. Proven raw counts must not be silently presented as the recommended workflow: library size and a few highly expressed genes can dominate their variance. PCA is nevertheless mathematically defined on such a matrix, so an expert's explicit source choice is accepted with a prominent state warning rather than blocked. Unknown provenance likewise continues with an explicit limitation.

`mask_var="highly_variable"` requires a boolean variable-axis annotation. Merely having a column of integer 0/1 values is not sufficient in Scanpy 1.12.3. When a mask is used, PCA scores are computed from the selected variables, while `varm["PCs"]` remains aligned to the complete current variable axis; excluded features receive zero loadings. Feature selection is a PCA input decision, not a physical subsetting requirement.

The current node accepts `X` or a current-axis named layer and deliberately excludes Raw. The repository defines Raw snapshot as the retained post-QC full-gene count state. It is not the default transformed PCA input and may have a different variable axis after feature selection.

## Centering, sparse input, solver, and component bounds

With `zero_center=True`, Scanpy performs ordinary centered PCA. Since Scanpy 1.5, supported sparse CSR/CSC inputs can be zero-centered implicitly instead of first creating a full dense cells-by-genes matrix. In Scanpy 1.12, efficient centered sparse PCA is supported with `arpack` and `covariance_eigh`.

The solver choices are not interchangeable implementation trivia:

- `arpack` uses truncated SVD, supports dense and SciPy sparse matrices, and uses `random_state`. scikit-learn requires strictly `0 < n_components < min(n_samples, n_features)` for this solver.
- `covariance_eigh` is deterministic and suited to tall-and-skinny inputs, but materializes a feature-by-feature covariance matrix and is less numerically stable because forming the covariance doubles the condition number.
- `randomized` is approximate and stochastic. Scanpy notes that it may not be exactly reproducible across platforms; it is not the efficient centered sparse path documented for the pinned implementation.
- `zero_center=False` dispatches to truncated SVD rather than ordinary PCA and changes the method's interpretation. It should not be hidden behind a setting still reported simply as PCA.
- `chunked=True` uses IncrementalPCA, ignores the solver and random-state choices, and densifies sparse chunks. It is a separate scalability/method tradeoff rather than an invisible optimization.

For a selected matrix with `n_obs` rows and `n_selected_vars` columns, the fixed ARPACK design must validate

```text
1 <= n_comps < min(n_obs, n_selected_vars)
```

before calling Scanpy. Centering also limits the mathematical rank to at most `n_obs - 1`; constant selected features should be counted and disclosed because they can produce numerically negligible components even when the shape bound passes.

Sparse input does not imply sparse outputs. `obsm["X_pca"]` is a dense cells-by-components array (requested as float32), `varm["PCs"]` is a dense current-genes-by-components loadings array (observed as float64 with the pinned stack), and variance arrays are dense. A workflow wrapper should preflight and disclose the lower-bound size of these persistent outputs.

## Randomness and reproducibility

scikit-learn documents that `random_state` is used by the ARPACK and randomized solvers. Passing an integer makes repeated calls reproducible in the same compatible software stack. The seed, solver, output dtype, and dynamic software versions must be reported. Reproducibility is not a promise of bitwise identity across arbitrary BLAS, SciPy, architecture, or package-version changes.

## Official output contract

With the canonical Scanpy-v1 keys (`key_added=None`), the call writes:

- `.obsm["X_pca"]`: cell scores, shape `n_obs × n_comps`;
- `.varm["PCs"]`: full current-axis feature loadings, shape `n_vars × n_comps`, with zero rows for excluded features;
- `.uns["pca"]["variance"]`: explained variance/eigenvalue for each component;
- `.uns["pca"]["variance_ratio"]`: fraction of total selected-matrix variance explained by each component;
- `.uns["pca"]["params"]`: the Scanpy parameter record.

scikit-learn states that components are ordered by decreasing explained variance and estimates explained variance with `n_samples - 1` degrees of freedom. A truncated set's variance ratios do not generally sum to one. Loadings have an arbitrary sign: flipping one loading vector and its corresponding scores leaves the same PCA solution, so tests and reports must not attach biological meaning to sign alone.

Repeated calls overwrite the three canonical result locations. A wrapper must treat `X_pca`, `PCs`, and `pca` as one result bundle and require explicit overwrite rather than combine fields from different runs.

## Scientific practice and references

PCA is an exploratory representation and dimensionality-reduction method. It is not clustering, Technical-batch correction, feature selection, or Condition inference. The number of retained PCs is a downstream modeling choice; explained variance and sensitivity analyses are more defensible than describing the default 50 as biologically optimal.

- Pearson K. On lines and planes of closest fit to systems of points in space. *Philosophical Magazine*. 1901;2:559-572. https://doi.org/10.1080/14786440109462720
- Pedregosa F et al. Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research*. 2011;12:2825-2830. https://www.jmlr.org/papers/v12/pedregosa11a.html
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*. 2019;15:e8746. https://doi.org/10.15252/msb.20188746
- Virshup I et al. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9:4371. https://doi.org/10.21105/joss.04371

Luecken and Theis describe PCA as general-purpose linear summarization, commonly after selecting informative genes, and recommend inspecting variance explained rather than confusing PCA with a visualization-specific nonlinear embedding.

## Required disclosure

`summary` should report the expression source and the resolved state **and state evidence as separate machine-readable parameters**, full and selected dimensions, mask name and number of selected/constant features, requested and computed components, solver/centering/dtype/seed, sparse input and dense-output memory estimate, exact output keys, per-PC variance and variance ratio, cumulative variance ratio, top absolute loading features per leading PC, warnings and limitations, and dynamic Python/openbio-singlecell/Scanpy/AnnData/scikit-learn/NumPy/SciPy versions. References must include the PCA method plus the actual Scanpy/scikit-learn software stack. Provenance remains advisory: `unknown` or a count-like state changes warnings and interpretation, not whether finite centered PCA can execute.

`code` should be a self-contained equivalent function that resolves the same `X`/layer source, audits expression provenance where available, validates boolean HVG selection, finite values, ARPACK bounds, overwrite policy, and output budget, calls public `scanpy.pp.pca` with all hidden choices explicit, and returns an equivalent copied `AnnData` without plugin-only history.

The only numeric wrapper bounds are the ARPACK component relation, positive component count, and a positive user-declared output-memory budget. Fixed UI ceilings such as 4,096 components or 128 GiB, and an arbitrary 0.01-GiB floor, are not Scanpy or ARPACK validity conditions and must not hide otherwise valid expert inputs.
