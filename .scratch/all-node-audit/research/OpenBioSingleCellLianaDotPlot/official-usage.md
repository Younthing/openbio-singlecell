# OpenBioSingleCellLianaDotPlot — official usage research

## Audited baseline

The current node reads an untyped table from `adata.uns`, guesses a method from a separate input, and calls
`liana.pl.dotplot` with `filterby` and `filter_lambda`. The deliberately audited LIANA 1.9.0 public signature uses
one `filter_fun` callable; those two keywords are unsupported. LIANA 1.10.0 is newer as of 2026-08-28 but is not
accepted without a separate audit. The node is therefore API-incompatible.

Its single `significance_threshold` is also scientifically wrong across methods. For CellPhoneDB,
`cellphone_pvals` is a cell-label-permutation specificity p-value; for rank aggregate, `specificity_rank` is an
aggregate rank quantity, not a generic p-value. Applying the same `<=0.01` label to both misstates the result. The
communication analysis is now Sample-resolved, so the official `dotplot_by_sample` is the appropriate renderer.

## Official LIANA 1.9 plotting usage

Official interfaces:

- https://liana-py.readthedocs.io/en/stable/generated/liana.plotting.dotplot.html
- https://liana-py.readthedocs.io/en/stable/generated/liana.plotting.dotplot_by_sample.html
- https://liana-py.readthedocs.io/en/stable/notebooks/basic_usage.html
- https://liana-py.readthedocs.io/en/stable/notebooks/mofatalk.html

For a direct pooled table, LIANA supports:

```python
p = li.pl.dotplot(
    liana_res=table.copy(deep=True),
    colour="magnitude_rank",
    size="specificity_rank",
    inverse_colour=True,
    inverse_size=True,
    source_labels=["CD4T"],
    target_labels=["B"],
    top_n=20,
    orderby="magnitude_rank",
    orderby_ascending=True,
    filter_fun=lambda frame: frame["specificity_rank"] <= 0.05,
    figure_size=(12, 8),
)
```

For `.by_sample` output, the documented call is:

```python
p = li.pl.dotplot_by_sample(
    liana_res=selected_table.copy(deep=True),
    sample_key="sample",
    colour="magnitude_rank",
    size="specificity_rank",
    inverse_colour=True,
    inverse_size=True,
    source_labels=["CD4T"],
    target_labels=["B"],
    figure_size=(12, 8),
)
```

`dotplot_by_sample` does not expose `top_n` or `filter_fun`; deterministic selection must happen before rendering.
Despite a `Figure` return annotation in 1.9, `return_fig=True` actually returns a plotnine `ggplot`; `.draw()` yields
the Matplotlib figure used for PNG export. LIANA 1.9's `_prep_liana_res` internally copies a supplied DataFrame before
inversion, so the previous claim that the public renderer necessarily mutates its caller's frame was incorrect.
OpenBio still passes a defensive copy and verifies pre/post fingerprints as part of its own stable contract.

## Method-specific visual semantics

For rank aggregate, lower `magnitude_rank` and `specificity_rank` are prioritized. Inversion is a visual transform,
not a changed scientific value; axes/legends and reports must state it. For CellPhoneDB, larger `lr_means` is stronger
expression magnitude while smaller `cellphone_pvals` is greater permutation specificity. The plot must not call an
unadjusted CellPhoneDB p-value FDR or Condition significance.

The plot is evidence selection. Every filter, sender/receiver subset, metric, direction, cutoff, top-N scope, input
and plotted row count, and omitted Sample must be reportable. A Sample-facet figure describes within-Sample
candidate interactions and does not compare Conditions statistically.

## References to emit

- Dimitrov D, et al. Comparison of methods and resources for cell-cell communication inference from single-cell
  RNA-Seq data. *Nature Communications*. 2022;13:3224. https://doi.org/10.1038/s41467-022-30755-0
- Dimitrov D, et al. LIANA+ provides an all-in-one framework for cell-cell communication inference. *Nature Cell
  Biology*. 2024;26:1613-1622. https://doi.org/10.1038/s41556-024-01469-w
- Efremova M, et al. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit
  ligand-receptor complexes. *Nature Protocols*. 2020;15:1484-1506.
  https://doi.org/10.1038/s41596-020-0292-x

## Report and generated-code implications

The report must disclose method-native colour/size/order fields and directions, any visual inversion, exact selection
and top-N scope, Sample facets, sender/receiver filters, plotted/input rows and interactions, resource/analysis
provenance, non-inferential limitations, references, and dynamic versions. Equivalent code must accept the portable
table/metadata from the explicit typed result, select deterministically on a copy, call the exact LIANA 1.9.0
plotting interface, draw its plotnine result, and return PNG bytes plus identical strict summary.
