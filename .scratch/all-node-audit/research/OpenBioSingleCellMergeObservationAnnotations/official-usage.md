# Merge Observation Annotations: official usage research

Researched: 2026-08-28

## Official AnnData and pandas semantics

AnnData stores per-observation annotations in `adata.obs`, a pandas `DataFrame` indexed by `obs_names`. Observation identity is therefore the correct seam for transferring one annotation between AnnData objects whose variable axes may differ.

- AnnData `obs` API: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.obs.html
- AnnData object and axis semantics: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html

Pandas `Series.reindex(target_index)` conforms values to a declared index and inserts missing values where labels are absent. Official duplicate-label guidance explains that non-unique labels make reindexing/alignment ambiguous. Both AnnData inputs must consequently have unique `obs_names`; silently calling `obs_names_make_unique()` would change cell identity and is not acceptable.

- pandas `Series.reindex`: https://pandas.pydata.org/docs/reference/api/pandas.Series.reindex.html
- pandas duplicate labels: https://pandas.pydata.org/docs/user_guide/duplicates.html
- pandas missing-value semantics: https://pandas.pydata.org/docs/user_guide/missing_data.html

An order-independent transfer follows the index rather than positional order:

```python
source = subset_adata.obs[source_column]
aligned = source.reindex(adata.obs_names)
output = adata.copy()
# classify missing, same, and conflicting values before assignment
output.obs[target_column] = resolved_values
```

The variable axes should not be required to match. A valid repository workflow transfers annotations from an HVG/model branch back to a full-gene branch while retaining the same cells.

## Scientific practice

Matching barcode text alone cannot prove that two AnnData objects derive from the same biological source. Unique exact identifiers are necessary but not sufficient. When both inputs contain OpenBio file/study provenance, the module should compare it; otherwise the report must state the limitation and require the user to confirm shared origin.

Transferring an annotation does not validate it. The target column name and report must preserve the repository distinction between **Provisional annotation** and **Curated annotation**. Conflict behavior must be explicit because overwriting a curated annotation with a provisional label is materially different from filling missing values.

## Method and software references

This is an index-alignment/data-provenance operation, not a statistical method. Relevant references are:

- Virshup I et al. anndata: Annotated data. *Journal of Open Source Software*. 2021;6(63):4371. https://doi.org/10.21105/joss.04371
- pandas official citation guidance: https://pandas.pydata.org/about/citing.html

## Required disclosure

Report source/target/shared/source-only/target-only observation counts; source missing values; newly filled, identical, conflicting, overwritten, and target-preserved values; source/target columns; conflict policy; final non-missing count and category counts; dtype changes; provenance evidence/limitation; and the fact that expression matrices and axes were not changed.

## Open expert-boundary re-review (2026-08-28)

Exact unique shared observation identifiers are sufficient to perform deterministic alignment. File/source provenance
can increase confidence but cannot prove biological identity and must not gate an expert-requested merge. Conflicting
or unavailable provenance, source-only observations, no shared observations, and all-missing source annotations are
therefore disclosed as warnings; source-only rows are ignored because they have no target slot. Duplicate labels and
the caller's explicit `conflict_policy="error"` remain hard because assignment would be ambiguous or contrary to the
declared policy.
