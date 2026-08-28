# PCA Metadata Associations: official usage research

Researched: 2026-08-28

## Audited decoupler interfaces and version break

The pre-refactor node called this decoupler 1.x function:

```python
decoupler.get_metadata_associations(
    data,
    obs_keys=None,
    obsm_key=None,
    use_X=False,
    layer=None,
    uns_key=None,
    inplace=False,
    alpha=0.05,
    method="fdr_bh",
    verbose=False,
)
```

- Archived decoupler 1.9.2 API: https://decoupler.readthedocs.io/en/v1.9.2/generated/decoupler.get_metadata_associations.html
- Fixed 1.9.2 source: https://github.com/scverse/decoupler/blob/v1.9.2/decoupler/utils_anndata.py#L942-L1035
- decoupler changelog: https://decoupler.readthedocs.io/en/stable/changelog.html

This is not the current public API. decoupler 2.x rewrote the package and moved the capability to:

```python
decoupler.tl.rankby_obsm(
    adata,
    key,
    uns_key="rank_obsm",
    obs_keys=None,
)
```

- Current official API: https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.tl.rankby_obsm.html
- Fixed 2.1.2 source: https://github.com/scverse/decoupler/blob/v2.1.2/src/decoupler/tl/_rankby_obsm.py
- Current official pseudobulk example: https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_psbk.html

The changelog explicitly warns that errors are expected when 1.x calls are used with decoupler 2.x. That made the pre-refactor node unusable in the supported base environment and left its major-version dependency unspecified. The refactored node no longer calls decoupler: it implements the reviewed Sample-level Spearman/one-way-ANOVA/global-BH family directly with SciPy and exposes the exact policy in its report.

## What the two official implementations actually test

decoupler 1.9.2 fits a separate formula `PC ~ metadata` for every PC/metadata pair using statsmodels OLS. Object/category metadata are wrapped as categorical terms; exact float/int dtypes are treated as continuous. It obtains a type-II ANOVA F-test and reports `eta_sq`, `pval`, and `p_adj`. Adjustment is performed separately within each PC over the selected metadata columns. `alpha` is passed to `multipletests`, but for ordinary `fdr_bh` the returned adjusted p-values are independent of alpha and the function discards the rejection flag, so alpha does not change the returned table.

decoupler 2.x changed the methods: numeric metadata use Spearman correlation, categorical metadata use one-way ANOVA, and all PC/metadata p-values are adjusted together using Benjamini-Hochberg. It returns/stores `obsm`, `obs`, `stat`, `pval`, and `padj`. It does not return categorical effect sizes, missingness counts, sample sizes, or an alpha-specific significance flag.

The current official example runs `rankby_obsm` on a pseudobulk object whose observations are biological samples, after PCA is computed on sample-level profiles. It does not present a cell-level PCA matrix with thousands of cells from the same Sample as thousands of independent biological replicates.

## Required distinction between continuous and categorical metadata

Metadata type changes both the hypothesis and the reported effect:

- Continuous metadata: a rank association such as Spearman's rho tests monotonic association. The signed rho is an effect-size/estimate and the p-value tests a null of no monotonic association under the function's assumptions.
- Categorical metadata: a one-way ANOVA F-test asks whether group means differ. Eta squared, `SS_between / SS_total`, is a unitless variance-explained effect size. The F statistic itself is not an effect size and has no direction.

Inferring type solely from storage dtype is unsafe. Integer-coded groups can be misclassified as continuous, while pandas nullable numeric/category/string dtypes may not match legacy exact-dtype checks. A workflow interface should make categorical and continuous key lists separate and explicit, reject overlap, and report the chosen type for every variable.

Official statistical interfaces relevant to an implementation are:

```python
scipy.stats.spearmanr(a, b, *, axis=0, nan_policy="propagate", alternative="two-sided")
scipy.stats.f_oneway(*samples, axis=0, equal_var=True, nan_policy="propagate")
scipy.stats.false_discovery_control(ps, *, axis=0, method="bh")
```

- Spearman correlation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html
- One-way ANOVA: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.f_oneway.html
- False-discovery control: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html
- statsmodels multiple-testing reference used by decoupler 1.x: https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html

Method and software references:

- Spearman C. The proof and measurement of association between two things. *The American Journal of Psychology*. 1904;15:72-101. https://doi.org/10.2307/1412159
- Lakens D. Calculating and reporting effect sizes to facilitate cumulative science: a practical primer for t-tests and ANOVAs. *Frontiers in Psychology*. 2013;4:863. https://doi.org/10.3389/fpsyg.2013.00863
- Virtanen P et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*. 2020;17:261-272. https://doi.org/10.1038/s41592-019-0686-2
- Badia-i-Mompel P et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances*. 2022;2:vbac016. https://doi.org/10.1093/bioadv/vbac016

`f_oneway` assumes independent observations, normal populations, and equal population variances for the standard test. Spearman correlation also assumes paired independent observational units for the usual p-value interpretation. Missing values must be handled per hypothesis and their removal disclosed. Constant numeric variables, one-level categorical variables, groups too small to estimate within-group variation, and constant PC scores make tests undefined or uninformative and must not silently produce a reportable association.

## Sample independence and the repository domain

The repository defines **Sample** as an independent biological specimen and the replicate unit for condition-level inference. Cells nested within the same Sample share biological and technical influences and are not independent biological replicates. Testing cell-level PC scores against Sample, Condition, or Technical batch labels with ordinary ANOVA/correlation inflates degrees of freedom and can generate pseudoreplicated p-values.

Zimmerman, Espeland, and Langefeld document that cells from the same individual are subsamples/pseudoreplicates and that treating them as independent inflates type-I error. The safe diagnostic design is therefore:

1. require an explicit Sample column with complete labels;
2. aggregate each PC to one unweighted mean score per Sample;
3. aggregate explicitly continuous metadata to one mean per Sample;
4. require categorical metadata to have exactly one non-missing value within each Sample;
5. run each association using Samples, not cells, as observations;
6. disclose cell counts and missing cells per Sample but never use cell count as replicate weight.

This aggregation makes the hypothesis a sample-level association of mean PC position. It is still exploratory. It can reflect cell-type composition, unequal cell recovery, preprocessing choices, or confounding and is not a causal effect or a formal gene-level Condition contrast. If a question concerns one Curated annotation, the object should be subset to that population before the PCA/diagnostic path.

- Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*. 2021;12:738. https://doi.org/10.1038/s41467-021-21038-1
- Squair JW et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. https://doi.org/10.1038/s41467-021-25960-2
- Repository ADR: `docs/adr/0001-separate-marker-evidence-from-condition-inference.md`

## Multiple testing and missingness

The family of hypotheses must be defined before adjustment. For this node, the transparent default family is all valid tested PC-by-metadata pairs from one invocation. Benjamini-Hochberg adjusted p-values and `significant = p_adjusted <= alpha` should be computed once across that family. The method, alpha, family size, skipped/failed pairs, and adjustment scope must be reported.

Benjamini-Hochberg controls the false discovery rate under its stated independence/positive-dependence conditions; PCs and metadata tests can be correlated, so adjusted p-values remain a screening device rather than confirmatory proof.

- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Missingness should be removed separately for each metadata variable after sample aggregation. The output must include total Samples, analyzed Samples, excluded Samples, missing cells, and group sizes. Silent formula-level row deletion is not adequate reporting. Finite `alpha` is accepted on the closed probability interval `0 <= alpha <= 1`; boundary values are expert choices that trigger a warning. Alpha is a decision threshold, not a parameter that changes ordinary BH adjusted p-values.

## Stable table and disclosure contract

The primary table should have one row per valid PC/metadata pair with stable columns:

```text
component, metadata, metadata_type,
n_cells_total, n_cells_missing,
n_samples_total, n_samples_analyzed, n_samples_missing, n_levels,
statistic_name, statistic,
estimate_name, estimate,
effect_size_name, effect_size,
p_value, p_adjusted, significant
```

For continuous variables, `statistic_name="spearman_rho"`, `statistic` and signed `estimate` are rho, and `effect_size_name="rho_squared"` with `effect_size=rho**2`. For categorical variables, `statistic_name="anova_f"`, estimate fields are null, and `effect_size_name="eta_squared"`. Strict JSON must convert undefined numeric values to null rather than NaN/Infinity.

`summary` should identify Sample as the independent unit, aggregation rules, representation and component count, explicit categorical/continuous variables, per-variable sample/group/missingness counts, test methods/effect sizes, global hypothesis count, adjustment method/alpha, the strongest associations by adjusted p-value and effect size, warnings/skipped tests, exploratory limitations, references, and dynamic Python/openbio-singlecell/AnnData/NumPy/Pandas/SciPy versions.

`code` should be a self-contained equivalent function that validates the same representation, Sample labels, metadata roles, within-Sample constancy, missingness, minimum replication, and finite values; aggregates without cell-count weighting; computes the same tests/effect sizes/global adjustment; and returns the same stable pandas table without plugin-only result wrappers.

## Open expert-tool boundary

Finite alpha is accepted on the closed interval `[0, 1]`; boundary values mean no/virtually all finite p-values pass and are warned, while values outside the probability range remain invalid. Explicit boolean continuous metadata are analyzed as 0/1 ranks with a warning that categorical declaration may be more interpretable. A categorical singleton group is weakly supported but can yield a finite ANOVA when another group has replication, so it is warned rather than rejected. Constant/one-level/insufficient metadata variables are skipped with machine-readable reasons while other valid variables continue; the call fails only when no finite hypothesis remains.
