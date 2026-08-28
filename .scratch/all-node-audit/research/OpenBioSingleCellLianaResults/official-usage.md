# OpenBioSingleCellLianaResults — official usage research

## Audited baseline

This node does not call a LIANA method. It copies a DataFrame from an arbitrary `adata.uns[result_key]` and wraps it
as a table without validating LIANA version, method, Sample design, resource, expression state, column semantics, or
score direction. Its existence is required only because the current communication node hides its primary result in
AnnData. This is a shallow pass-through rather than an atomic scientific analysis.

## Official LIANA result retrieval

The deliberately audited LIANA 1.9.0 methods accept `inplace=False` and return their result DataFrame directly. The
Sample-aware interface also returns a DataFrame directly, although its implementation still writes scratch state
and may categorical-cast the Sample column on the object it receives; callers must therefore use a private copy:

```python
import liana as li

work = adata.copy()
table = li.mt.rank_aggregate.by_sample(
    work,
    sample_key="sample",
    groupby="cell_type",
    resource=resource,
    use_raw=False,
    layer="log1p_norm",
    inplace=False,
)
```

Official documentation:

- https://liana-py.readthedocs.io/en/stable/generated/liana.method.AggregateClass.html
- https://liana-py.readthedocs.io/en/stable/generated/liana.method.rank_aggregate.__call__.html
- https://liana-py.readthedocs.io/en/stable/notebooks/mofatalk.html

LIANA 1.9's wrapper creates `adata.uns[key_added]` while iterating regardless of `inplace`; `inplace=False` controls
the final return but does not make the call read-only. When `inplace=True`, LIANA conventionally leaves the final
long table at that key, default `liana_res`. That storage option
is convenient in notebooks but is not a reason for a workflow system to discard the direct return value and expose
an untyped string key. A result table's interpretation depends on the producing method: rank aggregate uses lower
`magnitude_rank`/`specificity_rank`; CellPhoneDB uses higher `lr_means` and lower permutation
`cellphone_pvals`. A generic extractor cannot reconstruct those semantics safely.

## References

- Dimitrov D, et al. Comparison of methods and resources for cell-cell communication inference from single-cell
  RNA-Seq data. *Nature Communications*. 2022;13:3224. https://doi.org/10.1038/s41467-022-30755-0
- Dimitrov D, et al. LIANA+ provides an all-in-one framework for cell-cell communication inference. *Nature Cell
  Biology*. 2024;26:1613-1622. https://doi.org/10.1038/s41556-024-01469-w

## Consequence for this node

The correct official usage removes the need for this node. `OpenBioSingleCellLianaCommunication` should return a
validated, method-tagged, Sample-resolved typed result directly with its own `summary` and `code`. Plotting consumes
that explicit result. An arbitrary `uns` extractor would remain unsafe even if `summary`/`code` were bolted onto it.

## Public compatibility status

The old ID is retained only to explain why a hidden `adata.uns` table cannot be trusted or reconstructed. Its public
schema must be marked deprecated and non-executing so catalogue consumers cannot mistake it for an active LIANA
result operation. The direct typed Communication result is the only active method-bearing surface and owns the
LIANA references, versions, `summary`, and `code`.
