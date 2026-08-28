# Pseudobulk EdgeR contrast: official usage research

Researched: 2026-08-28

## Open expert boundary correction (normative)

The current edgeR User's Guide says multiple replicates in every Condition are desirable but not universally required; edgeR can handle a Condition with one library when the overall design still supplies biological replication and positive residual degrees of freedom. For this Pertpy QL path, the hard boundary is therefore an actually non-estimable model: invalid counts/axes, duplicate modeled Sample rows, absent contrast levels, a rank-deficient design, zero residual degrees of freedom for dispersion estimation, no retained genes, or backend failure. A fixed `>=3 Samples per Condition` gate is not an edgeR API prerequisite.

Low or asymmetric replication, incomplete batch/covariate overlap, non-Curated annotation, role aliases, and caller-declared Sample identity remain serious interpretation limitations. They must not block an otherwise computable expert analysis; they produce warnings, machine-readable invalidity reasons, and `formal_interpretation_invalid=true`. Exact Condition/nuisance aliasing that makes the design rank deficient still fails because the requested coefficient cannot be estimated. OpenBio does not claim that metadata or workflow history proves biological independence.

## Method scope and statistical unit

EdgeR models raw gene-count libraries with negative-binomial generalized linear models. In single-cell pseudobulk use, each library must represent one independent biological **Sample** within one selected population. Cells contribute counts to a Sample library but are not model rows or independent replicates. Technical batch is a separate Sample-level nuisance covariate.

The node's tested estimand should be one unambiguous pairwise Condition effect within exactly one population:

```text
log2 fold change = comparison Condition / reference Condition
```

A positive value means higher expression in the comparison Condition after adjustment for the declared additive covariates. Curated annotation is preferred for formal population-specific interpretation. Provisional or unknown annotations remain executable but force an interpretation warning and invalidity flag.

General-practice evidence that Sample-level aggregation avoids cell pseudoreplication:

- Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021;12:5692. https://doi.org/10.1038/s41467-021-25960-2
- Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*. 2021;12:738. https://doi.org/10.1038/s41467-021-21038-1
- Crowell HL, Soneson C, Germain P-L, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nature Communications*. 2020;11:6077. https://doi.org/10.1038/s41467-020-19894-4

## Current EdgeR release and canonical quasi-likelihood workflow

The current Bioconductor 3.23 release researched was edgeR 4.10.3 for R 4.6. The official guide's count-GLM workflow is:

```r
y <- edgeR::DGEList(counts = counts, samples = sample_metadata)
keep <- edgeR::filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::normLibSizes(y)  # calcNormFactors is retained as an alias
y <- edgeR::estimateDisp(y, design, robust = TRUE)
fit <- edgeR::glmQLFit(y, design, robust = TRUE)
qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
result <- edgeR::topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
```

Key official assumptions and practices:

- `DGEList` receives unnormalized non-negative integer counts with genes as rows and Sample libraries as columns;
- weakly expressed genes are filtered using `filterByExpr`, which considers library sizes and the design/group structure;
- TMM normalization (`normLibSizes`/historical `calcNormFactors`) adjusts effective library sizes; it does not transform the count matrix into analysis input supplied by the user;
- the design matrix must be numeric, full column rank, aligned to libraries, and leave positive residual degrees of freedom;
- dispersions are estimated from biological replication;
- `glmQLFit` plus `glmQLFTest` is the quasi-likelihood pipeline; edgeR guidance favors QL tests for stricter error-rate control, particularly with small replicate counts;
- the contrast must be the same length/order as the design columns and estimable under the design;
- `topTags(..., adjust.method="BH")` controls false discovery rate over the tested gene universe using Benjamini-Hochberg.

`robust=TRUE` makes dispersion/QL fitting less sensitive to hypervariable genes and is an appropriate fixed formal-analysis policy when supported by the installed edgeR release. It must be explicit and disclosed, not left to a changing wrapper default.

Primary official sources:

- EdgeR Bioconductor package page/version: https://bioconductor.org/packages/release/bioc/html/edgeR.html
- EdgeR User's Guide: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
- EdgeR reference manual (`filterByExpr`, `normLibSizes`, `estimateDisp`, `glmQLFit`, `glmQLFTest`, `topTags`): https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf
- Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*. 2010;26:139-140. https://doi.org/10.1093/bioinformatics/btp616
- Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*. 2010;11:R25. https://doi.org/10.1186/gb-2010-11-3-r25
- Lun ATL, Chen Y, Smyth GK. It's DE-licious: a recipe for differential expression analyses of RNA-seq experiments using quasi-likelihood methods in edgeR. *Methods in Molecular Biology*. 2016;1418:391-416. https://doi.org/10.1007/978-1-4939-3578-9_19
- Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research*. 2016;5:1438. https://doi.org/10.12688/f1000research.8987.2
- Chen Y, Lun ATL, Smyth GK. edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. *Nucleic Acids Research*. 2025. https://doi.org/10.1093/nar/gkaf018
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS B*. 1995;57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

## Current Pertpy 1.3.0 adapter

Pertpy 1.3.0 exposes a Python/R bridge with the constructor inherited from `LinearModelBase`:

```python
model = pertpy.tools.EdgeR(
    adata,
    design,       # formula string or numeric design matrix
    mask=None,
    layer=None,
)
model.fit(**glm_ql_fit_kwargs)
contrast = model.contrast(column, baseline, group_to_compare)
result = model.test_contrasts(contrast)
```

The captured 1.3.0 implementation validates integer counts and calls, in order:

```text
edgeR::DGEList
edgeR::calcNormFactors
edgeR::estimateDisp
edgeR::glmQLFit
edgeR::glmQLFTest
edgeR::topTags(n=Inf, adjust.method="BH", sort.by="PValue")
```

It uses `rpy2` and requires a valid R installation plus edgeR (its error also names BiocParallel and RhpcBLASctl). Its result columns are `variable`, `log_fc`, `logCPM`, `F`, `p_value`, and `adj_p_value`. The adapter does **not** call `filterByExpr`; feature filtering must occur before construction. `fit(**kwargs)` passes controls to `glmQLFit`, but `estimateDisp` is called without those kwargs. The OpenBio report/code must describe the actual wrapper pipeline rather than claim controls that it did not execute.

Pertpy's current official DE tutorial uses Decoupler pseudobulk and sample/expression filtering before constructing EdgeR. It is the relevant official Python example:

```python
pdata = dc.pp.pseudobulk(
    adata,
    sample_col="Patient",
    groups_col="Cluster",
    layer="counts",
    mode="sum",
)
dc.pp.filter_samples(pdata, inplace=True)
# Select one population and filter weak genes for the modeled libraries.
edgr = pt.tl.EdgeR(pdata, design="~Efficacy+Treatment")
edgr.fit()
contrast = edgr.contrast("Treatment", "pre", "post")
table = edgr.test_contrasts(contrast)
```

- Pertpy DE tutorial: https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html
- Pertpy EdgeR source captured at tag `1.3.0`, commit `e2976e13944bbdb3e30db25bb4b231f730e83ea5`: https://github.com/scverse/pertpy/blob/1.3.0/src/pertpy/tools/_differential_gene_expression/_edger.py
- Pertpy linear-model/contrast source: https://github.com/scverse/pertpy/blob/1.3.0/src/pertpy/tools/_differential_gene_expression/_base.py
- Heumos L, Schaar AC, Lance C, et al. Pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. 2025. https://doi.org/10.1038/s41592-025-02909-7

The repository's currently installed environment has Pertpy 1.3.0 but no `rpy2`; therefore the present EdgeR node cannot execute there. Dependency preflight must report Python/Pertpy/rpy2/R/edgeR versions or the exact missing component before a long model attempt.

## Count and design assumptions

Input must be the typed output of the audited pseudobulk aggregator. For one selected population require:

- one unique count library per independent Sample; no duplicate Sample rows and no cell-level rows;
- raw, finite, non-negative integer counts with nonzero library totals and unique genes;
- exactly the requested reference and comparison Condition levels represented; low/asymmetric replication is warned and invalidates formal interpretation, while execution depends on the complete design retaining positive residual degrees of freedom;
- complete Sample-level Technical batch/covariates with no string collisions;
- a deterministic additive design, for example `~ technical_batch + covariates + condition`, whose rows exactly match count libraries;
- full column rank, positive residual degrees of freedom, no exact Condition/batch/covariate aliasing, and a nonzero estimable comparison-minus-reference contrast;
- expression filtering derived on this selected population and modeled Sample set, not from pooled cells or all populations.

A free-form formula plus a separate `contrast_column/baseline/comparison` is unsafe: interactions, transformations, removed levels, and coding conventions can make the helper contrast differ from the analyst's intended estimand. The general node should expose metadata roles and build a simple additive design. Arbitrary interactions, pairing/random effects, donor blocking across repeated samples, nested technical replicates, or custom numeric contrasts need a separate advanced-design module with a typed design/contrast artifact.

Reject missing/blank levels, reference equal to comparison, unrepresented levels after QC, rank deficiency, zero residual df, all-zero genes after selection, no genes after filtering, exact Condition/nuisance aliasing, or non-estimable contrasts. Constant nuisance terms may be omitted because they add no column beyond the intercept, but the omission and declaration are reported. Partial imbalance or incomplete overlap that leaves a full-rank design remains executable with `formal_interpretation_invalid`; do not silently switch the contrast or use cells as residual replication.

## Output and reporting contract

The canonical table should retain every tested gene and use stable engine-neutral names:

```text
gene
log2_fold_change
log_counts_per_million
quasi_likelihood_f
p_value
p_adjusted
```

Sort order is presentational; gene identity and the tested universe must remain explicit. Values must be finite where the backend defines them, with JSON null rather than NaN/Infinity in summaries. `p_adjusted` is BH FDR over the reported tested-gene universe. Significance calls additionally require a declared `fdr_threshold` and optional absolute log2-fold-change threshold; changing display row limits must never change correction.

Required `summary` fields include:

- `methods`: Sample-level sum pseudobulk, population, raw-count source, expression filter, TMM effective-library normalization, negative-binomial dispersion estimation, robust QL fit/test policy actually executed, and BH correction;
- `results`: Sample counts by Condition, total/retained genes, design columns/rank/residual df, library-size/normalization summaries, null counts, FDR-significant up/down counts, and a bounded top-gene record;
- `key_results`: report-ready comparison-minus-reference prose with population and Sample replication, without causal language;
- `parameters`: all visible roles/thresholds plus fixed filter/model/correction controls and exact contrast vector;
- `warnings`: caller-declared Sample/source semantics, low/asymmetric replication, exact and partial confounding diagnostics, role aliases, provisional/unknown annotation, filtered universe, outliers/model warnings, and absent effect-size shrinkage;
- `references`: edgeR/QL/TMM/BH, Pertpy adapter, pseudobulk-practice and AnnData references;
- `software_versions`: openbio-singlecell, Pertpy, Decoupler if filtering used, rpy2, Python, R, edgeR, AnnData, NumPy, pandas, SciPy, and relevant R dependencies.

The `code` endpoint must contain executable equivalent function source using public APIs, the same typed-input/metadata/design/filter checks, the same explicit population subset and contrast direction, and the same backend calls/renaming. It must return the canonical result table and reproduce primary scientific state. It cannot be a schematic snippet or rely on ComfyUI/plugin-history helpers.
