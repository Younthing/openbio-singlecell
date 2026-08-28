# Use Expression Layer: official usage research

Researched: 2026-08-28

## Official AnnData and Scanpy semantics

`AnnData.layers` is a dictionary-like collection of matrices with exactly the same observation and variable dimensions as `X`. The official example reads a named layer as `adata.layers["spliced"]` and assigns a matrix to a layer with `adata.layers["spliced"] = ...`.

- AnnData layers API: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.layers.html
- AnnData object and aligned-axis slicing: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html

Assigning `adata.X = adata.layers[name].copy()` is therefore a valid AnnData storage operation, but it is not a biological method and AnnData does not define it as the preferred way to select an expression representation for an analysis. Many Scanpy operations expose representation selection at the consuming call, for example `layer` in normalization and highly-variable-gene interfaces or `use_raw` in reporting interfaces. The repository already represents the same choice as an `ExpressionSourceSpec` (`X`, a layer, or `raw`) at consuming nodes.

- Scanpy normalization API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html
- Scanpy highly variable genes API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html
- Scanpy expression-frame accessor: https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.obs_df.html

For a third-party routine that only reads `X`, the scientifically explicit pattern is to resolve the selected source in the wrapper and construct an isolated working `AnnData` whose `X` is that matrix. Only the intended result annotations are then copied back. This keeps the caller's expression state unchanged and makes the selected representation part of the analysis interface.

## Method and software references

This utility has no method paper because it only materializes one AnnData storage slot into another. The relevant software reference is:

- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. *Journal of Open Source Software*. 2021;6(63):4371. https://doi.org/10.21105/joss.04371

## Scientific implications

- The meaning of `X` is not self-describing. Replacing it makes every downstream node depend on workflow order unless each result also records the source layer.
- Copying a whole `AnnData` and then copying one expression matrix can be expensive for large single-cell datasets.
- The operation leaves `raw` untouched, so `X` and `raw.X` can silently refer to different scientific states.
- Expression-source choice is result-defining and belongs in the interface of the analysis that consumes it, not in an upstream global-state switch.
- If this compatibility node remains temporarily registered, its report must identify the layer, matrix shape/dtype/storage, the fact that `X` was replaced, and the fact that `raw` was not changed. It must be described as storage adaptation, not analysis.

## Open expert-boundary re-review (2026-08-28)

AnnData layers are shape-aligned even when either axis is empty. Materializing a named existing layer into `X` is
therefore defined for empty objects. Empty shape is disclosed as an advisory, not rejected. Missing layer names and
backed-object mutation remain true operation preconditions.
