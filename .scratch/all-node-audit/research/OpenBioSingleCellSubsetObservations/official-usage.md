# Subset Observations: official usage research

Researched: 2026-08-28

## Official AnnData usage

AnnData supports observation-axis indexing such as `adata[mask, :]`. Its documentation states that aligned elements are subset with the object: observation annotations, multidimensional observation annotations, pairwise observation matrices, layers, and `X` remain aligned. Indexing returns a view; calling `.copy()` materializes an independent result.

- AnnData object, indexing, views, and aligned slicing: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html
- AnnData raw behavior: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.raw.html

The appropriate implementation pattern for an immutable workflow node is therefore:

```python
series = adata.obs[column]
mask = ...  # one boolean value per observation
output = adata[mask.to_numpy(), :].copy()
```

AnnData `raw` is sliced on the observation axis too, so a correct observation subset keeps the same cells aligned across `X`, layers, `obsm`, `obsp`, annotations, and `raw`.

## Official pandas matching semantics

`Series.isin(values)` performs exact membership matching against a list-like collection. It does not implement substring or regular-expression matching. A direct inversion of the boolean result will include missing annotation values because a missing value is not a member of the requested set; the workflow must therefore define missing-value behavior separately.

- pandas `Series.isin`: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.isin.html
- pandas missing-data guide: https://pandas.pydata.org/docs/user_guide/missing_data.html

For the existing text widget, values can be compared to `series.astype("string")` after recording missingness. This supports string, categorical, numeric, and boolean annotations with one disclosed string-matching rule. Comma-separated input cannot represent a literal category containing a comma; that limitation must be reported until a structured-list widget exists.

## Scientific practice and reference

Selecting a cell population by an annotation is a dataset transformation, not a statistical test. It must not imply that an automated label is ground truth. The repository's **Provisional annotation** and **Curated annotation** terms remain distinct, and the report must name the actual column used.

Software reference:

- Virshup I et al. anndata: Annotated data. *Journal of Open Source Software*. 2021;6(63):4371. https://doi.org/10.21105/joss.04371

## Required disclosure

Report the annotation column, exact requested values, invert flag, missing policy, input/retained/removed cells and proportions, per-value counts, missing count, values not found, observation-order preservation, and the fact that all observation-aligned AnnData slots were subset together. An empty result is an error rather than a valid silent output.

## Open expert-boundary re-review (2026-08-28)

Official AnnData boolean slicing is defined when the mask selects zero observations. A wholly unmatched positive
selection or any final empty subset is therefore a valid explicit result and will be returned with warnings and zero
counts. The selected values, unmatched values, and empty-output status are disclosed in `summary`. Ambiguous string
coercion and an explicitly selected `missing_policy="error"` remain hard because the requested selection cannot then
be applied unambiguously.
