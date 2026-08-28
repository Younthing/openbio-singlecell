# Log1p: official usage research

Researched: 2026-08-28

## Official Scanpy interface and example

Scanpy documents `scanpy.pp.log1p(X, *, base=None, copy=False, chunked=None, chunk_size=None, layer=None, obsm=None)`. It computes

```text
X <- log(X + 1)
```

using the natural logarithm when `base=None`. The public function accepts an AnnData, dense array, or sparse matrix; for AnnData it can transform `X`, a named layer, or an `obsm` entry. `copy`, `chunked`, and `chunk_size` control execution and storage rather than the scientific question. The AnnData implementation records the logarithm base in `adata.uns["log1p"]` and warns when that marker already exists.

The official Scanpy preprocessing tutorial uses `log1p` after total-count normalization:

```python
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

For OpenBio, the first line is not part of this node: it is the explicit responsibility of `OpenBioSingleCellSnapshotExpression`.

- Scanpy `log1p` API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html
- Scanpy official preprocessing tutorial: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering.html#normalization
- Official Scanpy implementation: https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/preprocessing/_simple.py

## AnnData `raw` is a separate operation

`set_raw` is not a Scanpy `log1p` parameter. AnnData documents raw creation separately as `adata.raw = adata.copy()`, storing the then-current `X` and `var`. Observation slicing also slices `raw`, while variable slicing leaves it intact.

- AnnData raw API: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.raw.html
- AnnData paper: https://doi.org/10.21105/joss.04371

The repository gives **Raw snapshot** a narrower scientific meaning: post-QC, full-gene counts retained before highly variable gene selection. Assigning the log-transformed object to `raw` violates that meaning and can silently send transformed values to later count-dependent consumers. `Log1p` must never create, replace, or clear `raw`.

## Scientific input state and practice

Log1p is normally applied to non-negative, depth-normalized single-cell expression. The pseudocount keeps zeros at zero and makes fold-like differences less dominated by very large values, but the transformation does not perform library-size normalization, Technical batch integration, variance scaling, or Sample-level inference. A transformed matrix is a derived representation for exploratory feature selection, visualization, and related methods whose assumptions match it; it is not a count matrix for generative count models or pseudobulk count inference.

The mathematical domain for a finite real `log(1+x)` result is finite `x > -1`. Non-negative depth-normalized expression is the usual scientific input, but values in `(-1, 0)` are executable expert choices and therefore trigger a warning rather than rejection. Applying the function twice changes the scientific scale, but provenance is not infallible: proven or suspected prior transformation also triggers an explicit warning and disclosure rather than overriding the user's declared operation. Values `<= -1` remain a hard numeric error because they produce infinity or NaN.

Natural log is the repository convention. A different base changes every nonzero value but is rarely a routine scientific tuning parameter; fixing `base=None` and reporting it yields a smaller, clearer interface. Chunking/backed execution is an implementation/performance concern and should not be exposed unless the node adopts a separately verified backed-data contract.

## Edge conditions and failure behavior

- Reject zero cells, zero genes, non-finite values, and values `<= -1`; warn for finite negative values in `(-1, 0)`.
- Warn when provenance or `uns["log1p"]` indicates an already transformed matrix and record that the user explicitly requested another log transform.
- Reject backed AnnData with an actionable `to_memory()` instruction under the copy-on-write node contract.
- Preserve observation/variable order, annotations, every layer, and existing `raw`; if `raw` is absent, it must remain absent.
- Preserve sparse storage without densification where Scanpy supports it.
- Warn or disclose a limitation for absent, transformed, or otherwise inconsistent expression-state provenance; values alone do not prove normalization history.

## Methods and software references

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Ahlmann-Eltze C, Huber W. Comparison of transformations for single-cell RNA-seq data. *Nature Methods*. 2023;20:665-672. https://doi.org/10.1038/s41592-023-01814-1
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

The Scanpy publication is the software citation, the transformation comparison provides current method context, and the AnnData publication/API supports the storage distinction.

## Required `summary` and `code` disclosure

`summary` should state that active `X` was transformed by natural `log(1+x)`, report dimensions, storage type/dtype, nonzero count, input and output finite-value ranges/distributions, the resolved base, evidence about prior normalization/log transformation, and confirmation that layers and existing `raw` were preserved. Warnings and limitations must say that log1p output is not counts and that normalization state cannot always be inferred. Capture Python, openbio-singlecell, Scanpy, AnnData, NumPy, and SciPy versions (plus Pandas when used for reporting) and include the software/method references.

`code` should define a standalone copy-on-write function that validates the finite mathematical domain, warns on possible double transformation, calls `scanpy.pp.log1p` with the resolved natural base, postvalidates finite output, and leaves `raw` and layers unchanged. It must not recreate plugin history.
