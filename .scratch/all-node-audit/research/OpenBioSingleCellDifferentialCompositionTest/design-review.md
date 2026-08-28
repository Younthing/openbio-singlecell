# OpenBioSingleCellDifferentialCompositionTest — design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. Current Sample-level
> compositional methods remain; no graph migration or compatibility shim is retained.

Research date: 2026-08-28

## Decision

**Retire this production node after an explicit workflow migration. Do not merge it into Sample Composition Summary, and do not silently relabel its current results as a formal Condition contrast.**

The current Implementation does use one proportion per Sample, but its permanent Interface cannot express the assumptions that determine validity: independent/unpaired design, pairing or repeated subjects, covariates, Technical-batch confounding, composition-aware reference choice, and adequate replicate counts. It also tests each closed proportion separately. The broad name and convenient p-values create more apparent depth than the Module provides.

The deletion test favors retirement. Removing this node removes an under-specified inferential claim. The valid scientific work does not disappear into callers: descriptive aggregation remains in Sample Composition Summary, while formal inference belongs in a separately reviewed compositional-model Adapter such as scCODA or a replicate-aware design model such as propeller.

Do not automatically certify the repository's current scCODA node through this decision. It requires its own official audit, packaging check, validation, convergence diagnostics, summary, code, and tests before it becomes the recommended production path.

## Why a rename alone is insufficient

Renaming the display to "Exploratory Unadjusted Composition Screen" would improve honesty, but the existing Interface would still:

- accept designs it cannot analyze;
- run an asymptotic omnibus test at two Samples per Condition in the packaged example;
- mix omnibus and pairwise hypotheses and correction families;
- report a mean-proportion pseudocount fold change beside rank tests;
- encourage annotation-by-annotation interpretation despite compositional closure.

An exploratory rank screen can be scientifically defensible for independent, unpaired, one-factor Samples when its exact hypotheses and limits are explicit. That would be a new, narrower Module with a new node ID and its own research review. It should not inherit the misleading old ID or be created by default in this refactor.

## Correct Seam placement

The modules should be separated as follows:

1. Sample Composition Summary owns validation, the complete Sample-by-annotation zero grid, cell counts, proportions, and descriptive reporting.
2. A formal compositional Condition-contrast Adapter owns model design, covariates, reference-population semantics, inference, multiplicity, diagnostics, and method-specific reporting.
3. Plot modules consume a canonical composition table but perform no inference.
4. Workflow migration owns the analyst-visible choice of formal model; it must not guess reference populations or formulas.

The private zero-grid constructor is a useful Seam shared by descriptive and formal Adapters. A generic "test proportions" function is not: it would move design validity and compositional interpretation into every caller and reduce Locality.

## Proposed formal replacement Interface

The existing OpenBioSingleCellSccodaDifferentialComposition is the nearest candidate, subject to its own audit. Its eventual atomic Interface should include:

Visible inputs:

- adata;
- sample_key;
- annotation_key;
- covariate_keys;
- formula;
- reference_cell_type;
- estimated_fdr.

Advanced inputs:

- num_samples;
- num_warmup;
- random_seed;
- explicit convergence/diagnostic thresholds and resource guards if the backend supports them.

Outputs:

- table;
- summary;
- code.

Fixed policy:

- aggregate a complete Sample-by-annotation integer count matrix, including zeros;
- require unique observation identifiers and complete Sample/annotation/covariate metadata;
- require every sample-level covariate to have exactly one value per Sample; never accept pertpy's skip-and-continue behavior;
- preserve the Sample as inference unit;
- require a Curated annotation for formal claims, or prominently mark the run exploratory;
- record formula coding and the annotation reference separately from the biological Condition baseline;
- preserve caller AnnData;
- report sampler/convergence status before presenting credible effects.

The formal table should stabilize the official scCODA effect information into at least:

| Column | Meaning |
| --- | --- |
| covariate_term | Exact encoded design term |
| annotation | Curated population label |
| reference_annotation | scCODA reference population |
| final_parameter | FDR-thresholded posterior effect parameter |
| posterior_hdi_low | Lower reported posterior interval |
| posterior_hdi_high | Upper reported posterior interval |
| posterior_sd | Posterior standard deviation |
| inclusion_probability | Posterior inclusion probability |
| credible_effect | Whether the requested estimated FDR selects the effect |
| expected_sample_cells | Model's expected count at the documented mean sampling depth |
| log2_fold_change_relative | Model-reported relative compositional fold change |

Column names and interval level must be pinned against the audited pertpy version. Intercepts should be a separately identified table section rather than mixed with effects without a scope field.

## Formal summary and code contract

The replacement summary must contain:

- methods and report-ready key results;
- Samples per Condition and complete grid/zero diagnostics;
- exact formula, design matrix columns, Condition coding, and Technical-batch/covariate terms;
- annotation status and reference population, including whether selected automatically;
- mean and range of cells per Sample;
- estimated FDR, posterior inclusion cutoff, credible-effect count, and bounded strongest effects;
- NUTS warmup/draws, seed, divergences, effective sample size, R-hat or every diagnostic exposed by the backend;
- relative-composition and absolute-abundance limitations;
- method/software references and dynamically collected versions for openbio-singlecell, pertpy, anndata, pandas, NumPy, and the actual inference stack;
- strict finite JSON or explicit nulls, never NaN/Infinity.

The equivalent code must define one callable that reproduces aggregation, strict sample-covariate validation, preparation, reference choice, RNG controls, sampling, FDR selection, stable table extraction, diagnostics, and caller immutability. It can omit Comfy result wrapping and timing. It cannot omit the complete zero grid, design coding, or convergence checks.

## Retirement migration

1. Remove OpenBioSingleCellDifferentialCompositionTest from new-workflow registration after a graph migration path exists.
2. Remove it from Sample Composition Comparison and change that example to a descriptive Sample composition workflow. Do not insert scCODA automatically until the scCODA Adapter passes its own review.
3. Replace README claims of a generic composition comparison with explicit descriptive wording and a link to the reviewed formal model.
4. For saved workflows, show an actionable migration notice explaining:
   - old control_group is a Condition baseline;
   - scCODA reference_cell_type is an annotation reference and is not the same setting;
   - comparison_groups cannot determine a model formula or covariate structure;
   - pseudocount has no direct formal-model equivalent.
5. Require the analyst to select a formula, covariates, Condition coding, and reference population. Automatic conversion would manufacture scientific decisions.
6. Remove the legacy CSV/export links or reconnect them to the selected formal table only after schema review.

If release policy requires one executable deprecation cycle, use a fail-closed compatibility Interface:

- retain the legacy positional inputs only for graph loading;
- rename the display to "Deprecated: Exploratory Unadjusted Composition Rank Screen";
- add a visible acknowledgement defaulting to false that states the design is independent, unpaired, one-factor, and exploratory;
- refuse execution unless the acknowledgement is true and every selected Condition has at least five Samples;
- append summary and code after table;
- emit a mandatory retirement and closure warning on every successful run;
- never call the result a formal differential-composition analysis.

This temporary path is a migration aid, not the recommended refactor.

## Temporary compatibility table/summary/code

If the fail-closed cycle is implemented, simplify the legacy analysis rather than reproduce the mixed table:

- exactly two Conditions: one two-sided Mann–Whitney test per annotation, with an explicit tie-aware p-value method;
- three or more Conditions: one Kruskal–Wallis omnibus test per annotation, with no automatic pairwise tests in the same Module;
- one BH family across the annotation tests actually emitted;
- no generic pseudocount log2 fold change.

Temporary table columns:

| Column | Meaning |
| --- | --- |
| annotation | Annotation tested |
| test | Mann–Whitney U or Kruskal–Wallis |
| condition_levels | Ordered selected Conditions |
| n_samples_by_condition | Strict-JSON mapping |
| median_proportion_by_condition | Strict-JSON mapping |
| mean_proportion_by_condition | Strict-JSON mapping |
| zero_samples_by_condition | Strict-JSON mapping |
| statistic | Rank-test statistic |
| p_value | Unadjusted p-value |
| p_adjusted | BH-adjusted value over the emitted annotation family |
| adjustment_family_size | Number of annotations tested |
| p_value_method | Exact, asymptotic, or permutation policy |
| exploratory_only | Always true |

For two Conditions, a rank-biserial/common-language effect can be included if its direction and tie convention are pinned. For more than two, do not invent a directional omnibus effect. Condition-specific medians and distributions are descriptive.

The temporary summary must state the exact independent/unpaired assertion, Sample counts, tie/zero frequencies, test method, correction family, absence of covariates, closure warning, annotation status, and retirement path. Its results may say "exploratory rank screen"; they may not say a population is differentially abundant in absolute terms.

The temporary code must reproduce the complete zero grid, all validation, fixed branch (two-condition Mann–Whitney versus multi-condition Kruskal–Wallis), explicit SciPy p-value method, one BH family, and strict table schema. It must compile and match values/failures without mutating AnnData.

## Required retirement and replacement tests

### Migration tests

- every packaged workflow contains no legacy node after migration;
- old saved graphs receive the actionable migration notice;
- control Condition is never mapped to scCODA reference annotation;
- graph migration preserves unrelated widgets and links;
- registry/examples/README remain consistent.

### Formal replacement tests

- complete Sample-by-annotation zeros and cell-count conservation;
- Sample-to-covariate one-to-one validation, including the pertpy skip case;
- missing/blank/colliding identifiers and duplicated cells fail;
- Curated versus Provisional annotation reporting;
- formula/design coding and Condition baseline pinning;
- reference-population validation and automatic-reference disclosure;
- fixed RNG reproducibility and caller immutability;
- malformed/missing effect tables fail rather than becoming NaN;
- sampler diagnostic failures prevent report-ready claims;
- strict JSON and actual software/reference lists;
- generated-code table/diagnostic equivalence;
- output and working-memory preflight.

### Temporary compatibility tests, only if retained

- default acknowledgement fails closed;
- paired/repeated/covariate claims are explicitly unsupported;
- fewer than five Samples in any Condition fails;
- two Conditions emit only Mann–Whitney rows;
- three or more emit only omnibus Kruskal–Wallis rows;
- tied zeros use and report the pinned p-value method;
- all-identical proportions yield a defined neutral result or a clear failure, never unexplained NaN;
- BH is computed over exactly one disclosed family;
- no pseudocount fold-change column survives;
- summary/code/input-immutability contracts hold.

## Final cohesion assessment

Sample Composition Summary is a deep descriptive Module. A reviewed scCODA or propeller Adapter can be a deep inferential Module. The current Differential Composition Test sits between them: it duplicates aggregation, exposes several scientific switches, and still leaves the decisive design and closure obligations to the caller. Retirement yields the cleaner, safer module structure.

## Registry disposition

Retain the legacy ID only with `is_deprecated=True`, an explicit “Migration Required” or “Retired” display name,
and a description naming Sample Composition Summary versus scCODA/tascCODA as distinct choices. Registration tests
must verify those catalogue signals and fail-before-read behavior. The shim remains exempt from `summary`/`code`
because it performs no analysis.
