# Legacy Decoupler pseudobulk contrast: official usage research

Researched: 2026-08-28

## Finding

There is no current official Decoupler API corresponding to this node's analysis. The node is implemented against Decoupler 1.x functions that were deprecated during the 2.0 rewrite and are absent from current Decoupler 2.2.0. The current official changelog explicitly says:

- `get_pseudobulk` was renamed to `pp.pseudobulk`, and profile QC moved to `pp.filter_samples`;
- `get_contrast`, `get_top_targets`, and `format_contrast_results` were deprecated;
- PyDESeq2 should be used instead.

Therefore the old call chain is not a stable current method to preserve or silently reimplement. Current supported formal inference should be `OpenBioSingleCellPseudobulk` followed by a typed pseudobulk EdgeR or PyDESeq2 contrast.

Primary current sources:

- Decoupler changelog: https://decoupler.readthedocs.io/en/stable/changelog.html
- Decoupler 2.2.0 changelog/source captured at tag `v2.2.0`, commit `7e8e957cbbb2230079cdd3c13a0ac114a67b9655`: https://github.com/scverse/decoupler/blob/v2.2.0/CHANGELOG.md
- Current Decoupler pseudobulk API: https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.pp.pseudobulk.html
- Current official pseudobulk/DE vignette: https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_psbk.html
- Current Pertpy differential-expression tutorial using EdgeR/PyDESeq2: https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html

## Historical Decoupler 1.9.2 API used by the node

For migration/audit only, Decoupler 1.9.2 exposed:

```python
decoupler.get_pseudobulk(
    adata,
    sample_col,
    groups_col,
    obs=None,
    layer=None,
    use_raw=False,
    mode="sum",
    min_cells=10,
    min_counts=1000,
    dtype=np.float32,
    skip_checks=False,
    min_prop=None,
    min_smpls=None,
    remove_empty=True,
) -> anndata.AnnData

decoupler.get_contrast(
    adata,
    group_col,
    condition_col,
    condition,
    reference=None,
    method="t-test",
) -> tuple[pandas.DataFrame, pandas.DataFrame]

decoupler.format_contrast_results(logFCs, pvals) -> pandas.DataFrame
```

The historical implementation of `get_contrast` does not fit a pseudobulk count GLM. For every `group_col` population, it subsets pseudobulk profiles, retains only the Condition field, checks that the requested Condition and reference each have at least two profiles, then calls:

```python
scanpy.tl.rank_genes_groups(
    sub_adata,
    groupby=condition_col,
    groups=[condition, reference],
    reference=reference,
    method=method,
)
```

It accepts Scanpy rank-genes methods such as `t-test` and `wilcoxon`, returns only log-fold-change and unadjusted-p-value matrices, and its formatter applies per-contrast FDR adjustment. The current OpenBio node further normalizes each pseudobulk profile with `scanpy.pp.normalize_total`, log-transforms it with `scanpy.pp.log1p`, then runs this historical test.

Historical primary sources:

- Decoupler 1.9.2 API index: https://decoupler.readthedocs.io/en/v1.9.2/api.html
- Decoupler 1.9.2 `get_pseudobulk` documentation: https://decoupler.readthedocs.io/en/v1.9.2/generated/decoupler.get_pseudobulk.html
- Decoupler 1.9.2 official pseudobulk notebook: https://decoupler.readthedocs.io/en/v1.9.2/notebooks/pseudobulk.html
- Captured 1.9.2 implementation, commit `14396390202af82516f5857db03f16ea5f8a48ab`: https://github.com/scverse/decoupler/blob/v1.9.2/decoupler/utils_anndata.py
- Scanpy `rank_genes_groups` API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html

The captured source is the authoritative description of this legacy computation. It also reveals important behavior that a workflow cannot safely infer from the node name: insufficient groups are printed to stderr and skipped rather than producing a structured failure; `reference=None` becomes comparison to all other levels; and no design matrix or nuisance covariate is accepted.

## Why the historical method is not the supported formal replacement

Aggregating by Sample means the historical t-test/Wilcoxon is not cell-level pseudoreplication **if** every profile truly is one independent Sample. It is nevertheless a poor general formal Condition node for this repository because:

- it discards raw count mean-variance/library-size structure through normalize-total/log1p and does not use count-model offsets;
- it has no design matrix and cannot adjust Technical batch or any Sample-level covariate;
- it iterates every population behind one call instead of making the tested population and estimand explicit;
- skipped/missing strata can disappear from output without a structured invalid-analysis result;
- `reference="rest"` mixes heterogeneous Conditions and is not the requested two-level contrast contract;
- historical log-fold-change calculation and test semantics depend on Scanpy's rank-genes implementation and transformed matrix;
- the API and result formatter have been retired by their maintainer, whose migration guidance is PyDESeq2.

For a deliberately simple two-group transformed-expression exploratory test, a newly named generic Sample-profile t-test/Wilcoxon node could someday be designed and cited on its own merits. It must not masquerade as the current Decoupler-supported pseudobulk contrast and should not be created without a concrete user need.

## Current supported replacements

The formal pipeline is:

```text
cell-level raw counts
  -> OpenBioSingleCellPseudobulk
       (sum by independent Sample and population)
  -> select exactly one Curated population
  -> OpenBioSingleCellPseudobulkEdgeR
       or OpenBioSingleCellPseudobulkDESeq2 (PyDESeq2)
```

Both replacements must preserve Sample as the inferential unit, model raw summed counts, retain Technical batch separately as an additive covariate, validate a full-rank estimable design, filter weak genes on the modeled population/Sample set, disclose replicate counts, and emit a complete result table plus `summary` and `code`.

Relevant method references:

- Robinson MD, McCarthy DJ, Smyth GK. edgeR. *Bioinformatics*. 2010;26:139-140. https://doi.org/10.1093/bioinformatics/btp616
- Lun ATL, Chen Y, Smyth GK. EdgeR quasi-likelihood methods. *Methods in Molecular Biology*. 2016;1418:391-416. https://doi.org/10.1007/978-1-4939-3578-9_19
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics*. 2023;39:btad547. https://doi.org/10.1093/bioinformatics/btad547
- Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. https://doi.org/10.1038/s41467-021-25960-2

## Migration-time validation

A saved legacy workflow contains only generic AnnData wiring and string keys. It does not prove that its selected matrix is raw counts, that Sample is an independent biological replicate, that Condition/Technical batch are constant within Sample, or that a particular population was intended for one formal contrast. Migration must therefore be conservative:

- do not translate a legacy connection directly into EdgeR/PyDESeq2 unless the typed pseudobulk provenance and roles are explicit;
- do not preserve `reference="rest"` as a formal pairwise contrast;
- do not silently choose EdgeR versus PyDESeq2;
- do not carry the old normalized/logged matrix into a count engine;
- do not convert multiple historical group results into a claim that one design was jointly modeled;
- emit an actionable manual replacement path naming required Sample, population, Condition, optional Technical batch, count source, reference and comparison choices.

## Summary/code requirement under deprecation

The retired node should not remain an executable analysis node merely to attach `summary` and `code`; that would institutionalize an unsupported overlapping method. A non-executing workflow-load shim is not an analysis/processing node and must fail with migration guidance before producing scientific outputs. Consequently it should not emit a fabricated analysis summary or equivalent source.

The two supported active replacement engines each emit their own truthful `summary` and `code`. If migration tooling produces a machine-readable diagnostic, label it `migration_diagnostic`, not the scientific `summary`, and never represent it as an executed result.

## Public compatibility status

The stable node ID now denotes only a non-executing workflow-load shim, not an official Decoupler analysis. Its
public schema must therefore be marked deprecated and visibly describe the manual typed-pseudobulk replacement.
The shim has no method result to cite or report; citations, software versions, `summary`, and `code` belong to the
chosen active count engine after the user has supplied an explicit Sample-level design.
