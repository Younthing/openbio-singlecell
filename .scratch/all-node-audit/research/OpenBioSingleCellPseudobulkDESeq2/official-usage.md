# Pseudobulk PyDESeq2 contrast: official usage research

Researched: 2026-08-28

## Open expert boundary correction (normative)

PyDESeq2 0.5.4 does not define a universal per-Condition replicate minimum for constructing or fitting `DeseqDataSet`. Its `min_replicates=7` parameter controls whether Cook-outlier samples are eligible for replacement/refitting; it is not a minimum for the negative-binomial model itself. The hard boundary for this node is therefore invalid counts/axes, duplicate modeled Sample rows, absent contrast levels, a rank-deficient or saturated design with no residual degrees of freedom, no retained genes, unsafe required resources, or a backend fit/test failure.

The supported-version adapter resolves result-defining public defaults explicitly instead of inheriting patch-version drift. For PyDESeq2 0.5.4 it passes `control_genes=None`, `min_mu=0.5`, `min_disp=1e-8`, `max_disp=10.0`, `min_replicates=7`, and `beta_tol=1e-8` to `DeseqDataSet` through Pertpy, and passes `lfc_null=0.0`, `alt_hypothesis=None`, `prior_LFC_var=None`, `cooks_filter=True`, `independent_filter=True`, `alpha=<requested alpha>`, `quiet=False`, and the resolved `n_cpus` to `DeseqStats` through `test_contrasts`. Pertpy 1.3.0 itself hardcodes `refit_cooks=True`; the adapter verifies that supported interface and reports the resolved value rather than passing a duplicate keyword. `lfc_shrink=None` remains the OpenBio no-shrink policy and no shrink call is made.

PyDESeq2 stores the effective dispersion ceiling as `max(max_disp, n_samples)`, so `max_disp=10.0` is the explicit constructor input rather than always the realized ceiling. The report distinguishes those values and the post-fit check verifies the fitted `dds` attributes for fit/size-factor type, control-gene scope, `min_mu`, `min_disp`, effective `max_disp`, `refit_cooks`, `min_replicates`, `beta_tol`, quiet and low-memory policy.

Likewise, `size_factors_fit_type="ratio"` is the explicit request, not an unconditional claim about the realized normalization. PyDESeq2 0.5.4 automatically switches to iterative size-factor fitting, with a `UserWarning`, when every retained gene contains at least one zero across modeled Samples. OpenBio determines that official predicate on the fitted count matrix, reports requested and realized methods separately, retains the backend warning, and describes the realized method in report prose. This fallback is scientifically disclosed rather than rejected.

Post-fit validation must cover the backend-owned `model.dds`, not only the caller-owned/model adapter AnnData. The fitted `DeseqDataSet` must preserve the exact modeled count values, Sample and gene axes, every pre-existing `obs`/`var` field, and the exact numeric design matrix with its row/column identity. Backend-added fit metadata are allowed. Any mutation of those pre-fit inputs is a backend contract failure before results are reported; runtime and equivalent generated source enforce the same snapshot.

The object passed to Pertpy is deliberately minimal: copied counts, exact Sample/gene axes, OpenBio-namespaced axis-identity canaries, and the separately supplied numeric design. Artifact metadata remain in the fingerprint-validated preparation/report layer instead of being copied wholesale into `DeseqDataSet`, because PyDESeq2 legitimately owns and overwrites names such as `obs["size_factors"]`, `obs["replaceable"]`, and `var["dispersions"]`. The `dds` snapshot applies to this true minimal backend input, so backend-owned additions are allowed without sacrificing count/axis/design/canary immutability.

Low or asymmetric replication, incomplete nuisance overlap, non-Curated annotation, role aliases, and caller-declared Sample identity remain interpretation limitations. They are disclosed as warnings and machine-readable `formal_interpretation_invalid` reasons while a computable expert-selected model is allowed to run. Fewer than seven matching design replicates is disclosed separately as a Cook replacement/refit limitation and does not by itself invalidate interpretation. Exact Condition/nuisance aliasing remains a hard rank failure. OpenBio does not treat workflow history, Raw equality, or metadata naming as proof of biological independence or count provenance.

## Engine identity and scientific scope

Despite the node ID, the current implementation runs **PyDESeq2**, not the R/Bioconductor DESeq2 package. PyDESeq2 implements the DESeq2 negative-binomial workflow in Python and is separately versioned/cited. The stable node label, summary, result provenance, generated code and software versions must all say “PyDESeq2”; cite R DESeq2 for the statistical method and PyDESeq2 for the executed implementation.

For single-cell pseudobulk, each modeled row must be one independent biological **Sample** within exactly one selected population. Cells are summed within Sample/population; they are not independent model rows. Technical batch remains a separate Sample-level additive nuisance covariate. The tested pairwise estimand is:

```text
log2 fold change = comparison Condition / reference Condition
```

Positive log2 fold change means higher expression in the comparison Condition after adjustment for the declared additive covariates. Curated annotation is preferred for formal interpretation; provisional or unknown status remains executable but invalidates formal interpretation in the report.

General pseudobulk-practice sources:

- Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. https://doi.org/10.1038/s41467-021-25960-2
- Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*. 2021;12:738. https://doi.org/10.1038/s41467-021-21038-1
- Crowell HL, Soneson C, Germain P-L, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nature Communications*. 2020;11:6077. https://doi.org/10.1038/s41467-020-19894-4

## Current PyDESeq2 0.5.4 API

The current stable and project-pinned PyDESeq2 release researched was 0.5.4. Its main dataset signature includes:

```python
pydeseq2.dds.DeseqDataSet(
    *,
    adata=None,
    counts=None,
    metadata=None,
    design="~condition",          # formula or numeric design matrix
    fit_type="parametric",
    size_factors_fit_type="ratio",
    control_genes=None,
    min_mu=0.5,
    min_disp=1e-8,
    max_disp=10.0,
    refit_cooks=True,
    min_replicates=7,
    beta_tol=1e-8,
    n_cpus=None,                  # deprecated in favor of inference in current API
    inference=None,
    quiet=False,
    low_memory=False,
)
```

Local API verification used the project environment's Pertpy 1.3.0/PyDESeq2 0.5.4 with a full-rank intercept-plus-Condition design containing one reference and two comparison Samples (`residual_df=1`). The real adapter completed 250 Wald-test rows. A second real OpenBio end-to-end run (`X` -> Decoupler pseudobulk artifact -> Pertpy/PyDESeq2) completed all 120 tested genes with the same 1-vs-2 design and machine-reported `formal_interpretation_invalid=true`. PyDESeq2 warned that the parametric trend can fall back to `mean` and that dispersion-prior estimation with residual df below three is fragile; it did not reject the replicate pattern. This supports execution plus captured backend/low-replication warnings rather than a fabricated per-Condition hard minimum.

The current result/test interface includes:

```python
pydeseq2.ds.DeseqStats(
    dds,
    contrast,
    alpha=0.05,
    cooks_filter=True,
    independent_filter=True,
    prior_LFC_var=None,
    lfc_null=0.0,
    alt_hypothesis=None,
    inference=None,
    quiet=False,
    n_cpus=None,
)
```

Official minimal usage is structurally:

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(
    counts=counts,              # rows are biological samples; columns are genes
    metadata=sample_metadata,
    design="~ technical_batch + condition",
    refit_cooks=True,
)
dds.deseq2()

stat_res = DeseqStats(
    dds,
    contrast=["condition", "treated", "control"],
    alpha=0.05,
    cooks_filter=True,
    independent_filter=True,
)
stat_res.summary()
result = stat_res.results_df
```

PyDESeq2 requires raw non-negative integer counts with rows aligned to metadata samples. It estimates size factors, gene-wise dispersions, a dispersion trend and shrunk dispersion estimates, fits negative-binomial GLMs, and performs Wald tests. `DeseqStats` returns `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`. `alpha` participates in independent-filter optimization, not merely display. Cook's-distance outlier filtering and independent filtering can intentionally produce missing p-values or adjusted p-values; those must be counted/reported and serialized as JSON null rather than silently dropped.

- PyDESeq2 `DeseqDataSet` API: https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.dds.DeseqDataSet.html
- PyDESeq2 `DeseqStats` API: https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.ds.DeseqStats.html
- Official minimal pipeline: https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_minimal_pydeseq2_pipeline.html
- PyDESeq2 releases: https://github.com/owkin/PyDESeq2/releases
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics*. 2023;39:btad547. https://doi.org/10.1093/bioinformatics/btad547

## DESeq2 method reference

The current R/Bioconductor release researched was DESeq2 1.52.0 under Bioconductor 3.23. The foundational method models counts with negative-binomial GLMs, estimates size factors and dispersions with information sharing, and tests coefficients/contrasts. The R package is not executed by this node and its version must not be reported as runtime software unless actually used.

- DESeq2 Bioconductor package page/version and vignette: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8

## Current Pertpy 1.3.0 adapter

Pertpy 1.3.0 exposes:

```python
model = pertpy.tools.PyDESeq2(
    adata,
    design,       # formula string or numeric design matrix
    mask=None,
    layer=None,
)
model.fit(**deseq_dataset_kwargs)
contrast = model.contrast(column, baseline, group_to_compare)
result = model.test_contrasts(
    contrast,
    alpha=0.05,
    lfc_shrink=None,
    cooks_filter=True,
    independent_filter=True,
)
```

Captured 1.3.0 behavior relevant to reproducibility:

- it validates integer counts;
- it densifies sparse count input because of PyDESeq2 compatibility, which creates a material memory requirement that should be preflighted/reported;
- its base class converts a formula to a numeric Formulaic design matrix;
- `fit` constructs `DefaultInference`, passes the numeric design matrix to `DeseqDataSet`, fixes `refit_cooks=True`, and calls `dds.deseq2()`;
- `_test_single_contrast` constructs `DeseqStats`, calls `summary()`, and maps columns to `variable`, `baseMean`, `log_fc`, `lfcSE`, `stat`, `p_value`, and `adj_p_value`;
- LFC shrinkage is opt-in and accepts a named fitted coefficient, not an arbitrary numeric contrast. Pertpy explicitly does not infer a shrink coefficient from its numerical contrast vector.

The current official Pertpy DE tutorial first constructs sum pseudobulks with Decoupler and filters profiles/features, then constructs PyDESeq2 with an additive design. Relevant example:

```python
pdata = dc.pp.pseudobulk(
    adata,
    sample_col="Patient",
    groups_col="Cluster",
    layer="counts",
    mode="sum",
)
dc.pp.filter_samples(pdata, inplace=True)
# Select one population and weak-expression filter on its modeled samples.
pds2 = pt.tl.PyDESeq2(pdata, design="~Efficacy+Treatment")
pds2.fit()
contrast = pds2.contrast("Treatment", "pre", "post")
table = pds2.test_contrasts(contrast, alpha=0.05)
```

- Pertpy DE tutorial: https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html
- Pertpy PyDESeq2 API: https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.PyDESeq2.html
- Pertpy 1.3.0 adapter source, commit `e2976e13944bbdb3e30db25bb4b231f730e83ea5`: https://github.com/scverse/pertpy/blob/1.3.0/src/pertpy/tools/_differential_gene_expression/_pydeseq2.py
- Pertpy linear model/contrast source: https://github.com/scverse/pertpy/blob/1.3.0/src/pertpy/tools/_differential_gene_expression/_base.py
- Heumos L, Schaar AC, Lance C, et al. Pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025. https://doi.org/10.1038/s41592-025-02909-7

## Count filtering, design and contrast assumptions

Input must be the audited typed pseudobulk-count artifact. Select exactly one population, then require:

- one unique raw integer count library per independent Sample;
- nonzero library totals, unique genes and finite non-negative integer counts;
- exactly represented reference and comparison Conditions; low/asymmetric replication is warned and invalidates formal interpretation, while the complete design's rank and residual degrees of freedom decide computability;
- complete Sample-level Technical batch/covariates, with roles separate from Sample/Condition/population;
- a deterministic additive design such as `~ technical_batch + covariates + condition` aligned exactly to libraries;
- full column rank, positive residual degrees of freedom, no exact Condition/nuisance aliasing, and a nonzero estimable comparison-minus-reference contrast;
- weak-expression filtering computed on this population and modeled sample set before fit.

The same Decoupler `pp.filter_by_expr`/edgeR-derived filter policy may be shared with EdgeR to make the tested universe deliberate and auditable. PyDESeq2's later `independent_filter=True` is separate: it chooses a mean-expression threshold for multiple-testing power and can leave `padj` null. The report must disclose both prefit gene filtering and result-stage independent filtering.

Reject missing/blank levels, reference equal to comparison, levels lost after profile QC, duplicate Sample rows, invalid counts, rank deficiency, zero residual df, exact Condition-batch aliasing, all-zero/no retained genes, or non-estimable contrast. Constant nuisance terms may be omitted as redundant with the intercept if that omission is explicit in diagnostics. Partial imbalance remains executable with an invalidity warning; never switch the contrast or treat cells as residual replication.

A free-form formula plus a simple pairwise helper is unsafe for general interactions/transforms. The general node should expose additive metadata roles; arbitrary interactions, repeated measures/pairing, splines, custom numeric contrasts and nested technical replicates require a future typed advanced-design module.

## LFC shrinkage policy

Report unshrunk maximum-likelihood log2 fold changes in this core pairwise node. Pertpy's numeric contrast may combine coefficients after covariate encoding, while PyDESeq2 shrinkage targets a named coefficient. Guessing a coefficient can shrink a different estimand. Do not silently call `lfc_shrink`.

A future coefficient-specific shrinkage node may consume the fitted design artifact and require an exact named coefficient/shrink method. Until then, summary warnings should state that reported LFCs are unshrunk and may be unstable for low-count genes/small samples.

## Output and reporting contract

Return every fitted/tested gene in a canonical table:

```text
gene
base_mean
log2_fold_change
log2_fold_change_standard_error
wald_statistic
p_value
p_adjusted
```

`p_value` may be null for Cook's-filtered outliers and `p_adjusted` may be null after independent filtering; retain those rows and count each reason. BH/FDR calls use non-null adjusted p-values and declared `alpha`; display/top limits never change the tested/corrected universe.

Required `summary` fields:

- `methods`: Sample-level sum pseudobulk, selected population, raw-count source, prefit filter, size-factor/dispersion/NB GLM/Wald pipeline, Cook's policy, independent filtering, BH adjustment, and no LFC shrinkage;
- `results`: Sample counts by Condition, design columns/rank/residual df, genes before/after prefilter, fitted/tested rows, Cook's-null and independent-filter-null counts, size-factor/dispersion summaries, significant up/down counts and bounded top genes;
- `key_results`: report-ready comparison-minus-reference prose including population and independent replicate counts, without causal language;
- `parameters`: roles, filters, alpha/call thresholds, exact design/contrast vector and every fixed model policy;
- `warnings`: caller-declared Sample/source semantics, low/asymmetric replication, Cook-refit eligibility below seven, role aliases, exact/partial confounding diagnostics, provisional/unknown status, memory densification, null-result reasons, unshrunk effect sizes and model warnings;
- `references`: DESeq2 method, PyDESeq2 implementation, Pertpy adapter, BH and pseudobulk practice;
- `software_versions`: openbio-singlecell, Pertpy, PyDESeq2, Decoupler if filtering used, AnnData, Formulaic/Formulaic-Contrasts, NumPy, pandas, SciPy and Python. Do not report an R/DESeq2 runtime version unless R was actually called.

The `code` endpoint must provide executable equivalent function source using public APIs and the same typed-input/population/count/design/filter validation, explicit contrast direction, model controls, canonical renaming and null preservation. It must return the same scientific table and not depend on ComfyUI or plugin-history helpers.
