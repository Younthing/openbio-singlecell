# Normalize to Layer: official usage research

Researched: 2026-08-28

## Composition of documented Scanpy and AnnData interfaces

There is no single upstream Scanpy function named “normalize to a new layer.” The node composes three public interfaces:

1. `scanpy.pp.normalize_total` scales each cell by total expression and accepts `target_sum` plus a source `layer`.
2. `scanpy.pp.log1p` computes `log(1+x)` and can transform `X` or a named layer.
3. `AnnData.layers` is a dictionary-like collection of matrices aligned to current `X` dimensions; assigning a key creates or replaces that layer.

- Scanpy `normalize_total`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html
- Scanpy `log1p`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html
- AnnData layers: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.layers.html
- AnnData getting-started layer example: https://anndata.readthedocs.io/en/stable/tutorials/notebooks/getting-started.html#layers

The upstream `layer=` parameters modify the selected layer; they do not provide separate source and destination names. To preserve the selected counts and write a new derived layer, the documented primitives can be composed through a private working AnnData:

```python
output = adata.copy()
work = anndata.AnnData(X=selected_counts.copy())
sc.pp.normalize_total(work, target_sum=10_000)
sc.pp.log1p(work)
output.layers["log1p_norm"] = work.X.copy()
```

This is a composition of public interfaces rather than an upstream one-call recipe. The official Scanpy clustering tutorial likewise presents total-count normalization followed by log1p as a common preprocessing path, while explicitly preserving counts before those transforms. OpenBio delegates that count checkpoint to `SnapshotExpression`.

For `transform="sqrt"`, NumPy documents `numpy.sqrt` as an element-wise square root. It requires non-negative expression for real-valued output and is a distinct declared transformation, not an alias for log1p.

- NumPy `sqrt`: https://numpy.org/doc/stable/reference/generated/numpy.sqrt.html
- Scanpy official preprocessing tutorial: https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering.html#normalization

## Scientific input and output states

The selected source must be aligned count-like expression in current `X` or a layer, normally the post-QC canonical `layers["counts"]`. `raw` is intentionally excluded: AnnData raw may retain more variables than current `X`, whereas every layer must match the current observation-by-variable axes.

The result is a derived layer. With `transform="none"`, it contains linear total-count-normalized expression. With `log1p`, it contains natural `log(1+x)` of that normalized matrix. With `sqrt`, it contains the element-wise square root of the normalized matrix. Only the intermediate linear normalized values have per-cell sums equal to `target_sum`; after log1p or square root, final layer sums must not be reported as normalization targets.

The source must be an aligned finite numeric matrix. Non-negative integer-like UMI input with positive cell totals is typical, but that practice is not a universal execution precondition: Scanpy leaves zero-total rows unchanged and can scale finite signed/non-integer matrices. Those conditions and provenance inconsistent with counts therefore become warnings and exact disclosures. The chosen post-transform still has a mathematical domain: `sqrt` requires non-negative intermediate values and `log1p` requires every intermediate value to be greater than `-1`; violations are hard errors because they would create non-finite output.

This operation must leave current `X`, canonical counts, and `raw` unchanged. It is safe to fan the derived layer out to multiple exploratory/representation consumers, ensuring they use exactly the same preprocessing rather than recomputing slightly different matrices.

## Edge conditions and failure behavior

- Reject empty axes, non-numeric/non-finite or misaligned expression, and post-normalization values outside the selected transform's real finite domain. Warn for negative source values/totals, all-zero input, zero-total cells, non-integer values, and inconsistent/unknown count provenance.
- Require a finite `target_sum > 0`; do not use numeric zero as a sentinel for the dataset-derived Scanpy default.
- Reject unsupported transform strings during direct execution instead of silently treating them as `none`.
- Require a nonempty output key; reject the reserved `counts` name and an output key equal to the selected source layer.
- By default reject an existing output layer. If explicit advanced overwrite is supported, report the replaced state.
- Reject missing/blank source layers and backed AnnData with an actionable `to_memory()` instruction.
- Preserve sparse storage where the selected Scanpy/NumPy operation supports it and do not mutate the input.
- Preserve observation/variable order and all annotations, layers, and `raw` except for the one explicitly created/replaced derived layer.

## Methods and software references

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Ahlmann-Eltze C, Huber W. Comparison of transformations for single-cell RNA-seq data. *Nature Methods*. 2023;20:665-672. https://doi.org/10.1038/s41592-023-01814-1
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371
- Zheng GXY et al. Massively parallel digital transcriptional profiling of single cells. *Nature Communications*. 2017;8:14049. https://doi.org/10.1038/ncomms14049

The first and third are the required software/storage citations. The transformation comparison and the official Scanpy-linked Zheng example provide method context for total-count normalization and transformed derived representations.

## Required `summary` and `code` disclosure

`summary` should report the selected source, output-layer key, overwrite decision, dimensions, input storage/dtype and integer-like status, input per-cell total distribution, target total, intermediate normalized-total distribution, exact transform/formula, final layer value distribution/nonzero count, and confirmation that active `X`, `raw`, canonical counts, and unrelated layers were preserved. It must distinguish the linear normalization target from final transformed sums and warn that the layer is not counts. Include transform-appropriate references and dynamic Python/openbio-singlecell/Scanpy/AnnData/NumPy/SciPy/Pandas versions.

`code` should define a standalone function using public packages, resolve `X` or one named layer, reproduce the same warnings/domain/collision policy and finite postvalidation, compute through a private working matrix/object, assign only the requested derived layer on a copied AnnData, and reproduce the primary scientific result without plugin history.
