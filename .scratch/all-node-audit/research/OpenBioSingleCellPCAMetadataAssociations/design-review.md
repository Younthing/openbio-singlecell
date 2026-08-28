# PCA Metadata Associations: module design review

## Current module

At the start of this audit, the node accepted an arbitrary `obsm` key and comma-separated observation columns, imported decoupler at runtime, and called the decoupler 1.x `get_metadata_associations` function. It wrapped the returned frame as a `TableResult` with the generic description "ANOVA associations between an embedding and selected observation metadata." The refactored node instead owns an explicit Sample-level SciPy implementation and returns the reviewed table, `summary`, and `code` contracts below.

There are no direct tests of this node's statistics, data types, missingness, output columns, dependency version, or independence assumptions.

## Correctness and interface problems

1. The called function was removed from the decoupler 2.x top-level API and replaced by `decoupler.tl.rankby_obsm`; the repository declares no decoupler dependency, and the supported environment has no decoupler installation. The node cannot run as shipped.
2. The old decoupler implementation treats every row as independent. On cell-level `AnnData`, repeated Sample/Condition/Technical-batch labels create pseudoreplication and invalid inferential p-values.
3. A single `obs_keys` list delegates categorical/continuous classification to fragile storage-dtype heuristics. Integer-coded categories and nullable dtypes can be analyzed by the wrong method or rejected.
4. The description says ANOVA for every variable, although decoupler 1.x uses OLS/ANOVA for continuous covariates and decoupler 2.x uses Spearman correlation for them.
5. Missing values are silently dropped by formula fitting. The analyzed Sample count, excluded values, per-group sizes, and failed hypotheses are absent.
6. `alpha` appears scientific but does not change the ordinary BH-adjusted p-values returned by the old function, and the old function discards its rejection flag. The node does not add `significant` itself.
7. The free-form adjustment string exposes unsupported/misspelled methods and does not disclose whether adjustment is within PC or global. decoupler 1.x and 2.x use different scopes.
8. The old table reports eta squared but the node does not guarantee stable columns, effect interpretation, sorting, finite JSON values, or software/method references.
9. Constant PCs/covariates, one-level factors, too-small groups, non-finite representations, duplicate Sample identifiers, metadata varying within Sample, and too few independent Samples are not validated clearly.
10. The node returns no strict `summary` and no equivalent `code`.

## Decision: keep the diagnostic concept, replace the implementation, and enhance

Keep the node as a read-only exploratory diagnostic, but replace the undeclared/version-fragile decoupler call with a small, explicit SciPy implementation matching the current official conceptual split: Spearman association for continuous variables and one-way ANOVA for categorical variables. This avoids pinning a large enrichment package for three public SciPy statistics, gives the repository ownership of a stable table, and makes sample aggregation/effect-size reporting testable at the module interface.

Do not delete the diagnostic: examining whether high-variance PCs align with Technical batch, Condition, Sample-level QC, or other covariates is useful for deciding what to inspect or include in a later model. Do not merge it into PCA, integration, or condition inference. PCA is unsupervised and reusable; metadata hypotheses depend on study design; integration must not be selected automatically from a p-value; formal Condition contrast belongs to replicate-aware downstream modules under ADR-0001.

Do not generalize the name to arbitrary embedding associations in this pass. The module's interface and report should say PCA, validate a PC score matrix, and use `PC1...PCn` labels. A future generic representation-diagnostic module would need its own naming, variance/effect interpretation, and evidence.

## Target interface

- Visible scientific inputs: `adata`, PCA representation key (default `X_pca`), `sample_key` (default `sample`), comma-separated `categorical_obs_keys`, and comma-separated `continuous_obs_keys`.
- Advanced decision input: finite `alpha` on the closed interval `0 <= alpha <= 1`; boundary values are warned.
- Hidden fixed policy: unweighted per-Sample mean PC aggregation; per-Sample mean for continuous cell metadata; exact within-Sample constancy for categorical metadata; two-sided Spearman correlation; standard one-way ANOVA; eta-squared calculation; and global Benjamini-Hochberg adjustment. Singleton categorical levels are warned and retained when the test remains finite.
- Outputs: stable `TableResult` primary table, strict JSON `summary`, equivalent Python source text `code`.

Remove the free-form `p_adjust_method` for this cohesive node and fix global BH as the reviewed policy. If users later need dependence-robust BY or family-wise correction, add a constrained Combo only after a concrete workflow requires it. The adjustment method and scope remain fully disclosed despite being hidden.

Separate categorical and continuous inputs are intentional interface information, not avoidable implementation detail: they define different hypotheses. `sample_key` is required because Sample is the repository's biological replicate. Cell-level inference mode must not be offered.

## Validation and statistical contract

- Require a two-dimensional dense/CSR/CSC finite representation in `obsm`, aligned to `n_obs`, with at least one component and at least one cell. Convert only the relatively small score matrix to a dense numeric array.
- Require a present Sample column with complete, nonblank labels and at least three unique Samples overall. Aggregate every PC by Sample with equal Sample weight and report cell counts per Sample.
- Parse/deduplicate each metadata list; require at least one key total, present columns, no Sample key among tested metadata, and no key in both type lists.
- For a categorical key, require at most one distinct non-missing value within each Sample because mixed values would make the declared Sample-level factor ambiguous. Drop Samples whose value is entirely missing. Fewer than two observed levels or all-singleton groups skip that variable with a recorded reason; singleton or fewer-than-three-Sample levels otherwise warn and run when ANOVA remains finite. Do not turn cell types or other within-Sample annotations into sample-level factors.
- For a continuous key, accept only truly numeric finite-or-missing values, aggregate non-missing cell values to one mean per Sample, and drop Samples with no value. Fewer than three analyzed Samples or a constant Sample-level value skips that variable with a recorded reason; the whole call fails only when no finite hypothesis remains across all declared variables.
- Analyze missingness separately per metadata key. Include total/analyzed/missing Sample counts, total/missing cell counts, level/group sizes, and any skipped reason in the report. Fail the whole call if no hypothesis remains valid rather than return an apparently successful empty table.
- For every component/continuous key pair, compute two-sided `scipy.stats.spearmanr`; report signed rho as statistic/estimate and rho squared as variance-explained effect size.
- For every component/categorical key pair, compute standard `scipy.stats.f_oneway`; compute eta squared from between/total sums of squares; report F and eta squared, with estimate fields null. Fail or skip constant PC/degenerate results according to one documented rule and never serialize NaN/Infinity as JSON.
- Apply `scipy.stats.false_discovery_control(..., method="bh")` once to all finite p-values from the invocation. Set `significant = p_adjusted <= alpha`; sort stably by adjusted p-value, then descending effect size, component, and metadata.
- The table and summary must call this an exploratory Sample-level PCA diagnostic. They must not describe a significant row as a causal Technical-batch effect, a biological discovery, or formal Condition inference.

## Report contract

The primary table uses the exact stable column order documented in `official-usage.md`. `TableResult.parameters` records representation, Sample key, explicit metadata roles, alpha, aggregation policy, tests, and global adjustment.

`summary` includes:

- methods text stating that PC and continuous metadata were averaged per Sample, categorical values had to be Sample-constant, Samples were equally weighted, continuous pairs used two-sided Spearman, categorical pairs used one-way ANOVA/eta squared, and BH was global;
- results text naming the strongest associations with component, metadata, effect size, adjusted p-value, and number of Samples, while calling them exploratory diagnostics;
- key results for input cells/Samples/components, cells per Sample, requested and analyzed variables, per-variable missing/group sizes, total/valid/skipped hypothesis counts, significant count, and top rows by adjusted p-value/effect;
- warnings for small groups and limitations covering observational confounding, cell-composition effects, PCA/preprocessing dependence, test assumptions, correlated hypotheses, limited power with few Samples, and non-causal interpretation;
- Zimmerman pseudoreplication, Spearman/ANOVA/BH method, SciPy, and AnnData references plus dynamic Python/openbio-singlecell/AnnData/NumPy/Pandas/SciPy versions. decoupler is not listed as used software after replacement.

`code` is a self-contained function with resolved arguments. It performs the same parsing, representation/Sample/metadata/missingness/replication validation, equal-Sample aggregation, Spearman/ANOVA/eta-squared tests, global BH correction, stable columns/sorting, and returns an equivalent pandas DataFrame. It does not recreate plugin result metadata.

## Cohesion and coupling assessment

This redesign makes the node a deep in-process module: callers declare only the PCA representation, Sample identity, and metadata roles, while the implementation owns pseudoreplication prevention, aggregation, statistical dispatch, effect sizes, missingness, multiplicity, stable reporting, and citations. The decoupler adapter seam is removed because there is no second implementation to vary and the old adapter leaked more version/statistical complexity than it hid.

## Verification plan

- Construct synthetic Sample-level shifts where categorical metadata changes mean PC scores and continuous metadata has monotonic positive/negative association; compare F, rho, eta squared, p-values, and global BH output to direct SciPy calculations.
- Duplicate every cell within each Sample and verify p-values/effects do not change, proving cells are not treated as new independent replicates; vary cells per Sample and verify equal Sample weighting.
- Compare against a naive cell-level calculation in a regression test demonstrating that the node uses the smaller Sample degrees of freedom.
- Test categorical/object/string/bool values and explicitly declared integer-coded categories; continuous integer/float/nullable numeric values; overlapping/duplicate/missing keys; and metadata type declarations that do not match actual values.
- Cover mixed categorical values within one Sample, blank/missing Sample labels, entirely missing Sample metadata, constant continuous variables, one-level categories, one-Sample levels, small-level warnings, too few Samples, constant PC columns, non-finite scores, wrong representation shape, sparse representations, and no valid hypotheses.
- Verify per-variable missingness and group-size disclosure, stable columns/order/types, null rather than non-finite JSON values, global rather than per-PC adjustment, and alpha affecting `significant` but not ordinary BH adjusted p-values.
- Verify input AnnData is read-only, the primary `TableResult` is stable, `summary` parses under strict JSON, and generated `code` compiles and reproduces the table and failure modes.

The implementation keeps Sample identity, within-Sample categorical constancy, axis alignment, and finite result requirements hard. Advisory replication and metadata-quality failures are localized per declared variable and returned as skips/warnings, preserving valid hypotheses. Runtime and generated source are the same implementation and therefore must reproduce warnings and skip reasons exactly.
