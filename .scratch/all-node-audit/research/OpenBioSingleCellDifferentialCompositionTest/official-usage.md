# OpenBioSingleCellDifferentialCompositionTest — official usage research

Research date: 2026-08-28

## Scope

The current node is a collection of marginal rank tests on Sample-level proportions:

- one Kruskal–Wallis test per annotation across all selected Conditions;
- one two-sided Mann–Whitney U test per annotation and comparison Condition against the control;
- one Benjamini–Hochberg adjustment over the annotation-level omnibus p-values;
- a separate Benjamini–Hochberg adjustment over all annotation-by-comparison pairwise p-values.

It is not a joint compositional model, regression model, paired analysis, mixed model, or covariate-adjusted **Condition contrast**.

## Official SciPy semantics

### Kruskal–Wallis

Official call:

~~~python
from scipy import stats

result = stats.kruskal(
    control_proportions,
    treatment_a_proportions,
    treatment_b_proportions,
    nan_policy="raise",
)
~~~

SciPy documents Kruskal–Wallis as an H-test for two or more independent samples. Its p-value uses a chi-square approximation; the documentation gives a typical rule that each group have at least five measurements. Rejection is omnibus and does not identify which groups differ. For repeated measurements SciPy points to Friedman rather than Kruskal–Wallis.

Source: [SciPy kruskal reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kruskal.html).

### Mann–Whitney U

Official call:

~~~python
result = stats.mannwhitneyu(
    treatment_proportions,
    control_proportions,
    alternative="two-sided",
    method="auto",
    nan_policy="raise",
)
~~~

SciPy defines Mann–Whitney as a test on two independent samples. The general null concerns equality of the underlying distributions; interpreting it only as a median test requires stronger distributional assumptions. For paired/related samples SciPy directs users to the Wilcoxon signed-rank test.

With method="auto", current SciPy selects an exact calculation only when a sample is small and there are no ties; otherwise it uses an asymptotic calculation with tie correction. Proportion data commonly contain tied zeros and repeated fractions, so the method is often asymptotic even at small Sample counts. Exact mode does not correct for ties, and SciPy recommends a permutation method for small tied samples.

Sources:

- [SciPy mannwhitneyu reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)
- [SciPy wilcoxon reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)

## Experimental-unit and design limits

The vector passed to either rank test must contain one value per independent **Sample**. Cells are subsamples within a Sample, not independent replicates. The current aggregation does use one proportion per Sample, which is better than cell-level pseudoreplication, but the node cannot establish that the chosen sample_key truly identifies independent specimens.

The tests are appropriate only for a narrow design:

- one categorical Condition factor;
- independent, unpaired Samples;
- no repeated subject, donor blocking, pairing, family structure, or longitudinal correlation;
- no Technical-batch, sex, age, center, exposure, or other covariate requiring adjustment;
- adequate Sample counts for the selected p-value approximation;
- an exploratory marginal question about each annotation's proportion distribution.

Neither test can adjust a confounded Condition/Technical-batch design. Pairing cannot be repaired by keeping the same Sample identifier because that merely collapses distinct observations; it requires a paired or repeated-measures method. A one-line warning is not enough to turn an unsupported design into a valid formal contrast.

The official pertpy scCODA tutorial explicitly treats each statistical sample as one aggregated row and checks that covariates are unique within sample. The propeller manual calls its sample variable the biological replicate.

Sources:

- [pertpy scCODA tutorial](https://pertpy.readthedocs.io/en/1.0.1/tutorials/notebooks/sccoda.html)
- [speckle propeller manual](https://bioconductor.org/packages/release/bioc/manuals/speckle/man/speckle.pdf)

## Multiplicity

The current statsmodels usage is syntactically correct:

~~~python
from statsmodels.stats.multitest import multipletests

p_adjusted = multipletests(p_values, method="fdr_bh")[1]
~~~

Source: [statsmodels multipletests reference](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html).

However, a correction is meaningful only relative to an explicitly defined hypothesis family. The current combined output calls both values p_adj without saying that:

- omnibus BH covers one Kruskal–Wallis p-value per annotation;
- pairwise BH separately covers every selected annotation-by-control contrast;
- the two families are not adjusted together;
- with two Conditions, the omnibus and pairwise rows largely duplicate the same marginal question.

Benjamini–Hochberg addresses multiplicity among the supplied p-values. It does not adjust for pairing, confounding, incorrect Sample units, selection of comparisons after seeing results, or compositional closure. The original guarantee assumes independent tests, with extensions for certain positive dependence; closed proportions create structured dependence and a negative-correlation mechanism. At minimum the family and dependence limitation must be reported.

Primary reference: Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. JRSS B. 1995;57:289–300. [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).

## Compositional closure

Within every Sample, annotation proportions sum to one. An increase in one captured population forces the sum of all other proportions downward. The scCODA paper gives the direct single-cell consequence: independent univariate population tests can call shifts that are induced by this closure rather than distinct biological changes.

The current Kruskal–Wallis and Mann–Whitney calls operate on each annotation separately. BH correction does not transform those marginal hypotheses into a joint compositional model.

The current effect columns have additional interpretation problems:

- difference is the difference of arithmetic mean proportions, while the rank tests target distributional/rank differences;
- log2_fold_change is a log2 ratio of mean proportions after adding an arbitrary 0.001 proportion-scale pseudocount;
- the pseudocount affects only the reported effect, not either p-value;
- the result is relative captured composition, not absolute abundance;
- the generic name log2_fold_change can be mistaken for gene-expression or absolute cell-count fold change.

Sources:

- Aitchison J. The Statistical Analysis of Compositional Data. JRSS B. 1982;44:139–177. [doi:10.1111/j.2517-6161.1982.tb01195.x](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x)
- Büttner M et al. scCODA is a Bayesian model for compositional single-cell data analysis. Nature Communications. 2021;12:6876. [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)

## Composition-aware and replicate-aware alternatives

### scCODA

scCODA aggregates a Sample-by-cell-type count matrix and jointly models its columns with a hierarchical Dirichlet–multinomial model. It supports formula covariates, selects credible effects by posterior inclusion probability at a requested FDR, and requires a reference cell type or an automatic reference choice. Effects are relative to that population and reference choice can materially change interpretation.

Official pertpy usage:

~~~python
import pertpy as pt

model = pt.tl.Sccoda()
mdata = model.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_type",
    sample_identifier="sample",
    covariate_obs=["condition", "batch"],
)
mdata = model.prepare(
    mdata,
    formula="condition + batch",
    reference_cell_type="automatic",
)
model.run_nuts(mdata, rng_key=123, num_warmup=1000, num_samples=10000)
model.set_fdr(mdata, est_fdr=0.05)
~~~

The official effect table contains covariate, cell type, final parameter, posterior interval, inclusion probability, expected sample, and log2 fold change. The reference population, formula, sampling diagnostics, and FDR threshold are essential report fields.

Sources:

- [pertpy scCODA tutorial](https://pertpy.readthedocs.io/en/1.0.1/tutorials/notebooks/sccoda.html)
- [pertpy Sccoda source/reference](https://pertpy.readthedocs.io/en/1.0.1/_modules/pertpy/tools/_coda/_sccoda.html)
- Büttner M et al. [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)

The repository already registers OpenBioSingleCellSccodaDifferentialComposition, but that Adapter has not been certified by this Batch07 review. In particular, the official loader may skip a nonunique covariate; a production wrapper should fail instead. Migration must therefore wait for that node's own audit and tests.

### propeller

propeller computes Sample-level proportions, applies an arcsine-square-root or logit transformation, and uses limma empirical-Bayes linear modeling. Its lower-level t-test/ANOVA interfaces accept a design matrix and can model additional covariates. It leverages biological replication and applies BH across cell types. It is a useful replicate-aware alternative for suitable designs, but it remains a per-cell-type transformed-proportion method rather than scCODA's joint Dirichlet–multinomial model.

Source: Phipson B et al. [doi:10.1093/bioinformatics/btac582](https://doi.org/10.1093/bioinformatics/btac582).

### Milo

Milo tests abundance over overlapping graph neighborhoods and is useful when discrete annotations obscure continuous cell states. It answers a different question and is not a drop-in replacement for a discrete Curated-annotation composition table.

Source: Dann E et al. Differential abundance testing on single-cell data using k-nearest neighbor graphs. Nature Biotechnology. 2022;40:245–253. [doi:10.1038/s41587-021-01033-z](https://doi.org/10.1038/s41587-021-01033-z).

## Current implementation deviations

- The display name promises a Differential Composition Test even though the method is an exploratory marginal, unadjusted-design screen.
- group_key and group output use a discouraged synonym for Condition.
- Missing metadata cells are silently removed before analysis.
- No minimum Sample count is enforced or reported.
- The bundled demonstration has only two Samples per Condition, below SciPy's typical five-measurement guidance for the Kruskal–Wallis chi-square approximation.
- Independent/unpaired design, pairing, repeated subjects, covariates, and confounding are neither represented nor rejected.
- SciPy's Mann–Whitney method selection, continuity correction, ties, and exact/asymptotic status are not recorded.
- A Kruskal–Wallis ValueError is converted to NaN without distinguishing an all-identical neutral annotation from another failure.
- Global and pairwise rows with different estimands and correction families are concatenated into one shallow schema.
- Global testing is redundant for two Conditions; pairwise testing proceeds regardless of omnibus result for more than two.
- The arbitrary pseudocount and relative mean-proportion effect are under-disclosed.
- Annotation curation status, closure, absolute-abundance limitation, and Sample balance are absent.
- There is no summary or code output and no focused correctness test.

## Method and software references

- Kruskal WH, Wallis WA. Use of Ranks in One-Criterion Variance Analysis. JASA. 1952;47:583–621. [doi:10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441)
- Mann HB, Whitney DR. On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other. Annals of Mathematical Statistics. 1947;18:50–60. [doi:10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491)
- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate. JRSS B. 1995;57:289–300. [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- Virtanen P et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods. 2020;17:261–272. [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)
- Seabold S, Perktold J. statsmodels: Econometric and statistical modeling with Python. Proceedings of the 9th Python in Science Conference. 2010. [doi:10.25080/Majora-92bf1922-011](https://doi.org/10.25080/Majora-92bf1922-011)
- Büttner M et al. scCODA is a Bayesian model for compositional single-cell data analysis. [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)
- Phipson B et al. propeller. [doi:10.1093/bioinformatics/btac582](https://doi.org/10.1093/bioinformatics/btac582)

## Public compatibility status

The audited marginal-test implementation is retired rather than modernized in place. Its stable ID remains solely
so a saved workflow can surface migration guidance. The public schema must mark that compatibility surface as
deprecated and non-executing; it must not imply that the old rank tests are an endorsed current method. Active
descriptive and model-based replacements own their own references, versions, `summary`, and `code` outputs.
