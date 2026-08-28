# OpenBioSingleCellMarkerGenes — official usage research

## Scope and current environment

This record covers **Cluster marker evidence**: ranking expression features that characterize discovered cell groups and help expert annotation. It does not cover a Condition contrast, replicate-aware differential expression, automated annotation, or proof of cell-type identity.

Environment inspected on 2026-08-28:

- `scanpy 1.12.3`;
- `anndata 0.13.2`;
- `numpy 2.4.4`;
- `pandas 3.0.5`;
- `scipy 1.17.1`;
- `statsmodels 0.14.6`;
- `scikit-learn 1.9.0`.

The current node accepts `groupby`, four Scanpy methods, an X/Raw/layer expression source, `n_genes`, `pts`, and `random_seed`. It copies AnnData, converts a non-categorical grouping column to string categories, runs `scanpy.tl.rank_genes_groups` with every group compared with `rest`, extracts `scanpy.get.rank_genes_groups_df(group=None)`, and returns one `TableResult`. It has no `summary` or `code` output.

## Primary sources

- Scanpy `rank_genes_groups`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
- Scanpy `rank_genes_groups_df`: https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.rank_genes_groups_df.html
- Scanpy 1.12.3 ranking implementation: https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/tools/_rank_genes_groups.py
- Scanpy 1.12.3 DataFrame extractor: https://github.com/scverse/scanpy/blob/1.12.3/src/scanpy/get/get.py
- scikit-learn `LogisticRegression`: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html
- SciPy `ttest_ind_from_stats`: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind_from_stats.html
- statsmodels `multipletests`: https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
- Wilcoxon F. Individual comparisons by ranking methods. *Biometrics Bulletin*. 1945;1(6):80–83. https://doi.org/10.2307/3001968
- Student. The probable error of a mean. *Biometrika*. 1908;6(1):1–25. https://doi.org/10.2307/2331554
- Welch BL. The generalization of Student's problem when several different population variances are involved. *Biometrika*. 1947;34:28–35. https://doi.org/10.1093/biomet/34.1-2.28
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B*. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Ntranos V et al. A discriminative learning approach to differential expression analysis for single-cell RNA-seq. *Nature Methods*. 2019;16:163–166. https://doi.org/10.1038/s41592-018-0303-9
- Squair JW et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. https://doi.org/10.1038/s41467-021-25960-2
- Wolf FA, Angerer P, Theis FJ. SCANPY. *Genome Biology*. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
- Repository ADR: `docs/adr/0001-separate-marker-evidence-from-condition-inference.md`.

## Scanpy 1.12.3 interface

The installed signature is:

```python
scanpy.tl.rank_genes_groups(
    adata,
    groupby,
    *,
    mask_var=None,
    use_raw=None,
    groups="all",
    reference="rest",
    n_genes=None,
    rankby_abs=False,
    pts=False,
    key_added=None,
    copy=False,
    method=None,
    corr_method="benjamini-hochberg",
    tie_correct=False,
    layer=None,
    **kwds,
)
```

Official documentation says the function expects logarithmized data. OpenBio therefore defaults to a named
`log1p_norm` layer and reports the selected scale, but mutable transformation history is not a runtime credential.
Explicit current-axis `X` or layer selection is sufficient; finite count-like, normalized-but-unlogged, scaled,
residual, negative, or unknown representations remain calculable and are disclosed as nonrecommended expert choices.

`groups="all"` and `reference="rest"` rank each group against the union of every other observed group. If `groups` is restricted while `reference="rest"`, the reference still uses all groups, not only the returned subset. `n_genes=None` returns all tested genes; a positive integer limits the stored/exported ranking, not the family of hypotheses used for the per-group p-value adjustment. `rankby_abs=False` ranks positive markers before negative features. `pts=True` stores the nonzero fractions inside the group and in its rest reference.

`corr_method="benjamini-hochberg"` adjusts p-values separately for the tested gene family of each group. It controls an expected false-discovery proportion under its assumptions; it is not a posterior probability that an individual gene is a marker. Multiplicity across exploratory clustering choices, several resolutions, or repeated downstream filters is not repaired by the per-ranking BH calculation.

The hardened adapter removes globally constant genes only from the Scanpy backend call, because their test statistic is undefined and Scanpy otherwise emits warning-prone non-finite values. It then restores each constant gene as an explicit neutral hypothesis for every group: score and approximate log2 fold change are zero, raw and adjusted p-values are one, and both prevalence values reflect the observed `> 0` state. Benjamini-Hochberg adjustment is recomputed with statsmodels over the **complete** original gene family per group, including these neutral hypotheses. Deterministic score-descending, original-gene-order tie breaking occurs before any requested ranking truncation.

## Method-specific semantics and result incompatibility

Scanpy supports `t-test`, `t-test_overestim_var`, `wilcoxon`, and `logreg`, but these methods do not share one complete result contract.

- `wilcoxon` is Scanpy's rank-sum path. Its score is a standardized ranking statistic and its p-value uses a large-sample normal approximation, not an exact small-sample test. `tie_correct` is used only here and defaults to false. Droplet expression has many tied zero values, so the chosen tie policy must be explicit; enabling it improves statistical fidelity at additional runtime cost.
- `t-test` calls SciPy's unequal-variance `ttest_ind_from_stats(..., equal_var=False)`, i.e. a Welch-style comparison. The score is a t statistic.
- `t-test_overestim_var` uses the same unequal-variance calculation but substitutes the focal-group sample count for the rest sample count. Scanpy source calls this a variance-overestimation “hack” for small groups; it is not a separately validated biological model and must be named exactly in reports.
- `logreg` fits scikit-learn multinomial/binary logistic regression and ranks classifier coefficients. It is a discriminative feature-ranking method, not the same inferential test as the other three. Solver, penalty, regularization, class weighting, convergence, feature scaling, and the multiclass parameterization materially affect coefficients.

The official 1.12.3 extractor is decisive: `scanpy.get.rank_genes_groups_df` returns `names` and `scores` only when `method == "logreg"`; the `logfoldchanges`, `pvals`, and `pvals_adj` fields are returned only for non-logistic methods. In the binary case, Scanpy's logistic implementation yields one coefficient vector rather than symmetric per-class tables. Therefore a uniform table promising group-wise `logFC`, `p`, and `p_adj` cannot truthfully include `logreg`. Filling those columns with NaN is not an equivalent result and makes the current downstream filter reject every logistic row.

The cohesive marker-evidence node should support `wilcoxon`, `t-test`, and `t-test_overestim_var`. Remove `logreg` from this interface. A future discriminative-feature node would need its own coefficient-oriented table, explicit model settings, convergence diagnostics, and multiclass interpretation; it must not masquerade as a p-value/fold-change marker test.

## Stored and tabular result contract

For non-logistic methods, `adata.uns[key]` contains structured per-group arrays for names, scores, approximate log2 fold changes, raw p-values, adjusted p-values, parameters, and—when requested—`pts`/`pts_rest` DataFrames. `scanpy.get.rank_genes_groups_df(adata, group=None, key=...)` returns every stored group in long form.

The official documentation warns that `logfoldchanges` is an approximation derived from mean log expression. It is not a model coefficient and should be labeled `log2_fold_change_approx`. Sparse and dense implementations can differ slightly. Fractions mean the fraction of cells with a value greater than zero in the selected logarithmized expression representation; they are not average expression and become uninterpretable after centering/scaling.

The proposed stable primary table should contain:

- `group`, `gene`, and one-based within-group `rank`;
- `score` with a method-specific interpretation recorded in summary;
- `log2_fold_change_approx`;
- `p_value` and `p_adjusted`;
- `fraction_in_group` and `fraction_reference`.

Every non-key column must be numeric, finite where the method promises it, and aligned to the same `(group, gene)` row. Gene identifiers and group labels must be unambiguous strings without collisions after serialization. Duplicate `(group, gene)` rows are an error.

## Scientific preconditions and caveats

- Require an in-memory, nonempty AnnData with unique observation and selected-source feature identifiers.
- Require only a finite real numeric representation aligned to the selected current axes for safe execution. The report must warn for states other than logged expression, for count-like values, and for negative/scaled representations because fold-change and `>0` prevalence semantics change. These are interpretation disclosures, not history-backed code gates.
- Require a complete categorical grouping with at least two observed groups and no unused categories, missing labels, or string-rendering collisions. Each focal group and its rest reference need at least two cells for the supported inferential paths. Report all group sizes and warn on small groups rather than claiming stable evidence.
- Accept `n_genes=0` to mean all genes or require `1 <= n_genes <= selected_source_n_vars`; disclose whether the ranking was truncated. Retaining the current practical default of 100 bounds the primary table for interactive annotation, while exhaustive downstream filtering requires an explicit zero/all choice and a preflight output-row budget.
- Set `pts=True`, `groups="all"`, `reference="rest"`, `rankby_abs=False`, and `corr_method="benjamini-hochberg"` as fixed, disclosed policy for positive cluster-characterizing evidence. Expose `tie_correct` as an advanced Wilcoxon control and default it to true for tied single-cell values.
- Retain every tested gene in a separate ordered universe artifact even when the displayed ranking is truncated. Constant genes remain in that universe and in the complete multiplicity family instead of causing the batch to fail.
- Results depend on clustering resolution, preprocessing, gene universe, group size, Technical batch structure, and the selected reference. Marker ranking does not validate the clusters or determine a Curated annotation.
- Cell-level tests treat cells as observations. They are acceptable here only as exploratory Cluster marker evidence. Squair et al. and Scanpy explicitly warn that cell-level comparisons inflate p-values when used for Condition inference because cells from the same Sample are not independent. Formal Condition contrasts must use Sample-level replicate-aware methods.

## Reporting and equivalent code

Primary outputs are the stable marker-evidence `TableResult` and the complete ordered tested-gene-universe `TableResult`, read-only with respect to input AnnData. They share an analysis fingerprint and a cryptographic fingerprint of the complete untruncated ranking. Each artifact additionally carries a separately recomputable current-content fingerprint, so a consumer detects numeric table tampering as well as accidental pairing of identical axes/parameters produced from different expression results.

`summary` should report the selected expression source and state/evidence, source feature count, group/reference policy, group sizes, method and score meaning, tie policy when relevant, correction family, requested/actual rows per group, top bounded marker examples, distributions/ranges of scores, approximate log2 fold changes, adjusted p-values, and prevalence, plus missing/non-finite counts. It must describe the result as exploratory Cluster marker evidence and explicitly deny Condition inference or cell-type truth.

References should include the active method paper, BH, Scanpy, statsmodels, and the replicate-awareness caveat. Runtime versions should include Python, openbio-singlecell, Scanpy, AnnData, NumPy, pandas, SciPy, and statsmodels; scikit-learn is included only if a future method actually uses it.

`code` should be a self-contained Python function that validates the same expression/group contracts, builds a minimal temporary AnnData containing only the selected expression, grouping, and globally variable genes, calls public `scanpy.tl.rank_genes_groups` with all fixed and visible values explicit, restores neutral constant-gene hypotheses, performs complete-family BH, and returns `(table, universe)` as equivalent pandas DataFrames. It should not depend on ComfyUI, plugin history bookkeeping, or Scanpy private functions.

## Performance and identity contract added by independent audit

`max_output_rows` guards the **sum** of exported marker-table rows and tested-universe rows. Advanced `max_working_memory_gib` first guards identifier/group validation before materializing observation-scale Python lists, then guards a disclosed conservative peak estimate before backend execution. The estimate includes bounded extrema and integer-like scan buffers, observation/gene identifiers and serialized group membership, a float64 dense selected-expression copy or worst-case sparse values/indices/row pointers, eight group-by-variable-gene backend statistic/name/prevalence arrays, the complete internal ranking, the truncated marker output, and the universe output. Content fingerprints use streaming canonical JSON hashing rather than an additional complete Python/JSON row copy. Sparse execution may use less memory than that worst-case bound, but it is subject to the same preflight.

Canonical SHA-256 identities use sorted compact JSON without NaN. DataFrame content payloads preserve column and row order; string values remain exact and numeric values use `float.hex()` tokens. The complete ranking, current marker table, and current universe use distinct schema tags. Observation, gene, and group identifiers must be nonempty strings without surrounding whitespace; gene and observation identifiers must also be unique.
