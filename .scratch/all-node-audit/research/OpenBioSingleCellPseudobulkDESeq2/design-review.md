# Pseudobulk PyDESeq2 contrast: module design review

## Current module

`OpenBioSingleCellPseudobulkDESeq2` currently accepts generic AnnData, a free-form design formula (default `~group`), optional layer, and `contrast_column`, `baseline`, and `comparison`. It constructs Pertpy `PyDESeq2`, fits, builds a helper contrast, tests it, and returns only the upstream table.

The pre-change refactored engine consumed a typed pseudobulk artifact and established a reproducible design/filter/report contract, but treated Curated/formal provenance, role exclusivity, and a fixed replicate recommendation as execution authorization. The corrected implementation moves those items to disclosure.

The corrected boundary keeps rank, residual-degree, input-shape/count, tested-universe, memory/resource and backend failures hard. Weak/asymmetric replication, caller-declared Sample identity, annotation status, role aliases, constant nuisance terms and incomplete (but not exact) nuisance overlap remain executable with warnings and machine-readable `formal_interpretation_invalid` reasons.

The fixed supported-version policy is now explicit at both PyDESeq2 construction and test seams: `control_genes=None`, `min_mu=0.5`, `min_disp=1e-8`, `max_disp=10.0`, `min_replicates=7`, `beta_tol=1e-8`, `lfc_null=0.0`, `alt_hypothesis=None`, `prior_LFC_var=None`, `cooks_filter=True`, `independent_filter=True`, and `quiet=False`. The resolved CPU count is passed to both inference stages. Pertpy's adapter-owned `refit_cooks=True` is interface-checked and reported. These values are hidden because this atomic node does not expose alternative estimands/outlier policies, but none are left as undocumented PyDESeq2 patch defaults.

The `max_disp=10.0` item is the fixed constructor input; PyDESeq2 realizes it as `max(10.0, modeled Samples)`. Both requested and effective values are reported, and fitted `dds` attributes are checked rather than inferred from a successful call.

The same requested/realized distinction applies to size factors. Although the fixed request is `size_factors_fit_type="ratio"`, PyDESeq2 automatically realizes `iterative` when every retained gene has at least one modeled zero. That supported fallback remains executable, is detected from the fitted count state, and is disclosed in methods, parameters, warnings and runtime/generated parity tests.

Input immutability is checked on the object PyDESeq2 actually fits. After testing, `model.dds` must retain the exact pre-fit count matrix, Sample/gene axes, pre-existing observation/feature metadata, and design matrix including row/column identity; backend-added fit columns and arrays remain permitted. Checking only `model.adata` is insufficient because that is an adapter copy distinct from the fitted `DeseqDataSet`.

To avoid false failures and name collisions, “pre-existing metadata” at this seam means only OpenBio's namespaced Sample/gene identity canaries in a minimal backend AnnData. Full artifact `obs`/`var` metadata are already protected by the typed artifact fingerprint and remain available to preparation/reporting, but are not passed into a backend that owns common column names such as `size_factors`, `replaceable`, and `dispersions`.

## Correctness findings

### P0: the input type does not establish Sample-level raw counts

Generic AnnData allows a cell-level matrix, normalized pseudobulk values, or an unrelated layer to reach the negative-binomial model. This can produce pseudoreplication or invalid count input. Consume only `OPENBIO_SINGLE_CELL_PSEUDOBULK`, then select exactly one population and verify one row per independent Sample.

### P0: free-form design and simple helper contrast are semantically under-specified

An arbitrary formula can encode interactions/transformations while the UI separately requests one metadata-column pairwise contrast. The helper's numeric vector is not guaranteed to represent an intended marginal/conditional effect. The current result does not disclose encoded design columns, factor references, rank, residual degrees of freedom or contrast vector.

Use a bounded additive role interface for the general node. Reserve custom formulas/interactions/repeated-measures designs for a future advanced typed design/contrast artifact.

### Population and replicate strength must be reported without inventing an API minimum

Selecting one population and enforcing one modeled row per declared Sample remain structural requirements. PyDESeq2 has no universal per-Condition fit minimum: `min_replicates=7` controls Cook replacement/refitting. Low/asymmetric replication and non-Curated annotation therefore warn and invalidate formal interpretation but do not block a full-rank positive-residual design.

### Estimability/confounding and feature filtering are absent

The node does not reject missing levels, rank deficiency, zero residual df, Condition/Technical-batch confounding, invalid/constant covariates or non-estimable contrasts. It does not prefilter weak genes on the selected population/model rows. PyDESeq2 independent filtering occurs later and is not a substitute for a documented prefit universe policy.

### Fixed PyDESeq2 behavior is hidden

Pertpy densifies sparse input, uses a numeric design matrix, fixes `refit_cooks=True`, constructs a default parallel inference object, and by default uses Cook's filtering and independent filtering. The node reports none of these facts, null reasons, memory implications, versions or alpha semantics.

### The name can falsely imply R DESeq2

The class/node ID says DESeq2 while execution is PyDESeq2. Saved workflow identity may keep the ID, but display name, description, method report, generated import/code, citations and versions must identify PyDESeq2 exactly.

### Output/report/source contracts are missing

The upstream result schema is not normalized; rows with null p/padj lack explanation. There is no complete tested-universe record, design/sample context, summary, references, software versions, report-ready prose or equivalent code.

## Decision: keep and deepen as a focused PyDESeq2 engine

Retain the node ID for saved-workflow continuity but rename the user-facing operation to “Pseudobulk PyDESeq2 Contrast.” Its one atomic job is: select one population from typed pseudobulk counts, construct/preflight one bounded additive Sample-level design and pairwise contrast, apply a declared weak-expression filter, run PyDESeq2 once, and return a canonical complete table plus summary/code.

Do not merge with aggregation; the typed pseudobulk artifact is reusable and inspectable. Do not merge with EdgeR; normalization, dispersion, outlier, filtering and test outputs differ and need separate citations/versioning. Share only pure preparation/report primitives.

`OpenBioSingleCellDecouplerPseudobulkContrast` should migrate to this pipeline or EdgeR by an explicit analyst choice, never as an automatic internal alias.

## Target interface

Required inputs:

- `pseudobulk`: typed `OPENBIO_SINGLE_CELL_PSEUDOBULK`;
- `population`: exactly one represented population;
- `condition_key`: normally inherited/locked from artifact metadata;
- `reference_condition`;
- `comparison_condition`.

Visible optional inputs:

- `technical_batch_key`;
- `categorical_covariate_keys`;
- `continuous_covariate_keys`;
- `fdr_threshold`/`alpha`, default 0.05, explicitly documented as both reporting threshold and independent-filter target;
- `min_abs_log2_fold_change`, default 0, for reporting calls only.

Advanced prefit-expression filter inputs:

- `min_count`, `min_total_count`, `large_n`, and `min_prop`, matching the executed current Decoupler/edgeR-derived policy.

Performance-only advanced input:

- `n_cpus`, bounded to available positive CPUs and disclosed; it must not alter scientific defaults.

Fixed hidden scientific policy:

- one population/one row per Sample;
- additive intercept design ordered Technical batch, categorical covariates, continuous covariates, Condition;
- comparison-minus-reference contrast;
- `fit_type="parametric"`, `size_factors_fit_type="ratio"`, no control-gene subset;
- `control_genes=None`, `min_mu=0.5`, `min_disp=1e-8`, `max_disp=10.0`, `min_replicates=7`, `beta_tol=1e-8`;
- `refit_cooks=True`, `cooks_filter=True`, `independent_filter=True`, `lfc_null=0.0`, `alt_hypothesis=None`, `prior_LFC_var=None`;
- unshrunk MLE log2 fold changes (`lfc_shrink=None`);
- no top-N truncation of the table and no input mutation.

Outputs:

- `table`: canonical complete PyDESeq2 result;
- `summary`: strict JSON text;
- `code`: executable equivalent function source.

Any fixed PyDESeq2 controls must be checked against the installed supported version and passed explicitly where the public API permits. If Pertpy hardcodes a control, the adapter/version check and report record that exact behavior.

## Shared preparation seam with EdgeR

A private pure helper should consume the typed artifact and roles, then return copied one-population counts/metadata, deterministic Sample/gene order, encoded additive design, exact contrast vector, rank/residual-df/confounding diagnostics, Condition replicate counts, and the weak-expression mask/report.

This seam owns common invariants only. It must not choose a model, normalization, Cook's policy, result columns or engine citations. EdgeR and PyDESeq2 adapters remain independently replaceable/testable.

The filter is computed from exactly the selected population and modeled two-Condition Sample rows. Record genes before/after, thresholds and retained identities. PyDESeq2 independent-filter decisions are recorded separately after testing.

## Design and input invariants

- Artifact fingerprint/schema/count state is valid and annotation status permits formal use (otherwise explicit exploratory result only).
- Selected population exists; after choosing reference/comparison there is one unique row per Sample.
- Counts are finite, non-negative integers with positive libraries; genes and Samples are unique.
- Both Conditions are present/distinct; warn and mark formal interpretation invalid below three Samples in either level, while rank/residual checks determine computability.
- Role aliases are recorded and warned; categorical/continuous fields are complete and correctly typed. Constant nuisance terms are omitted as intercept-redundant and disclosed, while aliases that create a rank-deficient design still fail.
- Reference coding is deterministic and Condition reference equals `reference_condition`.
- Design row identity exactly matches count rows, rank equals number of columns, and residual df is positive.
- Exact confounding/aliased columns and non-estimable/zero contrasts are rejected, never auto-dropped.
- At least one gene survives filtering and the dense-memory estimate fits a fixed safety budget before Pertpy densifies.

Unsupported designs—paired/repeated observations, donor blocking across multiple rows, random effects, interactions, splines, nested technical replicates and arbitrary numeric user contrasts—fail with an actionable advanced-design message. Sample identity must never be included as a fixed effect when there is only one row per Sample.

## Canonical result contract and postconditions

Return one row per fitted gene:

```text
gene: string
base_mean: float >= 0
log2_fold_change: nullable float
log2_fold_change_standard_error: nullable float >= 0
wald_statistic: nullable float
p_value: nullable float in [0, 1]
p_adjusted: nullable float in [0, 1]
```

Nulls are preserved, not dropped or coerced to significance-safe numbers. Classify/report their backend cause where available: Cook's outlier filtering, independent filtering, or model failure. Validate unique exact gene set, value domains and full row accounting. Stable presentation sort is p-value then gene with nulls last; calls use non-null `p_adjusted <= alpha` plus the optional magnitude threshold.

## Summary contract

Use `methods`, `results`, `key_results`, `parameters`, `warnings`, `references`, and `software_versions`. Include:

- artifact fingerprint, caller-declared count source/Sample identity, annotation status, selected population, `formal_interpretation_invalid`, and invalidity reasons;
- exact reference/comparison direction and independent Sample counts;
- design rendering/columns/references/rank/residual df/contrast vector;
- genes before/after prefilter and threshold policy;
- size-factor and dispersion summaries, fit type and size-factor method;
- fitted rows, Cook's-p-null and independent-filter-padj-null counts, and complete universe accounting;
- FDR/magnitude-significant up/down counts and bounded top-gene records;
- fixed refit/Cook/independent-filter/no-shrink policy, alpha and CPU/memory disclosure;
- low/asymmetric-replicate, role-alias, confounding/overlap, provisional/unknown, low-count/unshrunk-LFC and backend warnings; Cook-refit eligibility below seven is a separate warning and not an invalidity reason by itself;
- DESeq2 method, PyDESeq2 implementation, Pertpy adapter, BH and pseudobulk-practice citations;
- runtime Python/openbio-singlecell/Pertpy/PyDESeq2/Decoupler/AnnData/Formulaic-Contrasts/numerical-stack versions, with no fictitious R DESeq2 version.

Report prose is associative, for example: “Within Curated population P, PyDESeq2 compared C versus R across nC and nR independent Samples while adjusting for B; K genes met FDR <= alpha.” If none meet threshold, state none. Cell totals may describe aggregation depth but never n.

## Equivalent-code design

Generate one executable public-package function that takes portable pseudobulk AnnData plus explicit artifact metadata/parameters, performs the same count/population/Sample/role/design/filter/memory validation, constructs the same numeric design and comparison-minus-reference vector, calls Pertpy PyDESeq2 with explicit supported controls, and maps/preserves the full result identically.

Runtime/code parity covers:

- Sample and gene identities/order;
- design columns/rank/references/contrast;
- prefilter universe;
- PyDESeq2 fit, Cook, independent-filter, alpha and no-shrink controls;
- canonical rows/columns/nulls and call counts.

Generated source may omit plugin history and timing only. It cannot switch to R DESeq2, use a simplified two-column metadata formula, omit covariates/checks, or drop null rows. It reports package-version incompatibility clearly.

## Dependency and memory policy

Check supported Pertpy and PyDESeq2 interfaces before the fit and include installed versions in any incompatibility error. No silent fallback to EdgeR, Scanpy, t-test or Wilcoxon is allowed. Because Pertpy densifies sparse counts, compute `n_samples * n_genes * dtype_size` plus a conservative fit overhead budget before construction and fail with actionable population/filter guidance when unsafe; do not mutate/densify the source artifact.

## Verification plan for the future implementation batch

- Typed socket rejects generic/cell AnnData and transformed count state.
- Multi-population fixture requires one population and proves one model row per Sample.
- Balanced and low/asymmetric additive Technical-batch/covariate designs produce deterministic columns/vector/direction; absent levels, exact confounding, rank deficiency, zero residual df and invalid covariates fail before backend.
- A one-versus-two Condition fixture executes with positive residual df and exact runtime/generated parity, reports low replication and Cook-refit limits, and sets formal interpretation invalid; one-versus-one remains a saturated-design hard failure.
- Filter mask is population/model-row specific; independent-filter nulls remain a distinct reported stage.
- Fake Pertpy/PyDESeq2 adapter asserts numeric design, `refit_cooks`, fit/size-factor policies, alpha, Cook/independent filtering, no shrink, CPU and canonical mapping.
- Runtime/generated parity tests mutate only backend-owned `model.dds` counts, Sample/gene axes, pre-existing `obs`/`var` metadata, or design and require the same actionable contract failure; the count case includes an integer above float64's exact-integer range so `+1` cannot be hidden by lossy fingerprint coercion. A real PyDESeq2 0.5.4 smoke verifies the unmodified fitted snapshot.
- Null p/padj rows survive with reason counts; significant calls use only valid adjusted values and never truncate correction universe.
- Dense-memory preflight is deterministic and source artifact remains unchanged.
- User-facing label/report says PyDESeq2; no R runtime version appears.
- Optional real 0.5.4 smoke test, strict JSON without NaN/Infinity, citation/version completeness, generated-code compilation and primary table/design/filter/null parity.
