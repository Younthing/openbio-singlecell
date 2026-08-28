# Pseudobulk EdgeR contrast: module design review

## Current module

`OpenBioSingleCellPseudobulkEdgeR` currently accepts generic AnnData, a free-form design formula (default `~group`), optional count layer, and `contrast_column`, `baseline`, and `comparison`. It constructs Pertpy EdgeR, calls `fit`, builds a helper contrast, calls `test_contrasts`, and returns only the upstream table.

The pre-change refactored EdgeR QL engine consumed the typed artifact and performed deterministic design/filter/backend validation, but treated three Samples per Condition, Curated/formal provenance, and all role distinctions as execution authorization. The corrected implementation treats these as interpretation disclosures rather than edgeR computation prerequisites.

The corrected design preserves hard failures only when the requested calculation is undefined or the backend cannot estimate it. Weak/asymmetric replication, caller-declared Sample identity, annotation status, role aliases, and incomplete (but not exact) nuisance overlap are disclosed as warnings and `formal_interpretation_invalid` reasons while execution proceeds.

## Correctness findings

### P0: generic AnnData permits cells or transformed values as model rows

Nothing proves that each row is one independent Sample or that the matrix contains summed integer counts. A cell-level AnnData can be connected directly, producing pseudoreplicated inference, while a log-normalized pseudobulk matrix fails late or is misinterpreted. Consume only the concrete `OPENBIO_SINGLE_CELL_PSEUDOBULK` artifact.

### P0: arbitrary formula and simple helper contrast are not one coherent estimand

The design may contain interactions or transformations while the UI asks for one column/baseline/comparison triple. With treatment coding, absent levels, interactions, or aliased columns, the simple helper contrast may not represent the requested marginal or conditional effect. The node does not report the encoded design columns or contrast vector.

Replace this with a deliberately bounded additive-role interface. Build and disclose the formula/design and an explicit comparison-minus-reference vector. Complex formulas/custom contrasts belong in a future advanced-design artifact, not an unchecked string field on the general node.

### No population selection or one-Sample-one-row enforcement

An aggregate containing multiple populations can enter one model, making the same Sample occur multiple times and treating correlated population profiles as independent libraries. The node must select exactly one visible `population` before modeling, then verify one retained row per Sample.

### Estimability remains hard; scientific-strength diagnostics become disclosure

Condition representation, rank, residual degrees of freedom and contrast estimability are computation-critical and remain preflight errors. Exact Condition/Technical-batch aliasing is rank deficiency and remains hard. Low/asymmetric replication, constant nuisance terms, role aliases, and incomplete overlap can still yield a valid numeric design; report them instead of rejecting. Constant terms are deterministically omitted as redundant with the intercept and listed in diagnostics.

### Weak-expression filtering is absent

The Pertpy adapter starts at `DGEList` and never runs edgeR `filterByExpr`. Testing a large low-count universe reduces power and can compromise dispersion behavior. Filter once on the selected population and modeled Sample/design set, disclose every threshold and the before/after universe, and compute BH only over the retained tested genes.

### Backend defaults and dependency state are hidden

Pertpy calls TMM normalization, dispersion estimation, QL fitting, QL testing and BH topTags, but the node discloses none of them. It does not explicitly request robust QL fitting and does not reveal Python/R bridge or R package versions. In the audited local environment `rpy2` is absent, so execution currently fails only at runtime.

### Output is upstream-shaped and under-reported

`variable`, `log_fc`, `logCPM`, `F`, `p_value`, `adj_p_value` are Pertpy implementation names rather than a stable OpenBio contract. There is no complete tested universe disclosure, Sample/design context, key result prose, references, software versions, strict JSON `summary`, or equivalent `code`.

## Decision: keep as a focused typed-pseudobulk EdgeR QL engine

Retain and deepen this node. Its atomic operation is: select one population from a validated pseudobulk-count artifact, construct/preflight one bounded additive Sample-level design and pairwise contrast, filter weak genes, run one EdgeR QL contrast, and return the complete canonical table plus report/source.

Do not merge it into `OpenBioSingleCellPseudobulk`: aggregation is a reusable upstream scientific state transition and should not be rerun for each model. Do not merge it with PyDESeq2: the engines use different normalization, dispersion/outlier and test semantics and must remain separately citable/reproducible. They may share a private typed-artifact selector/design-validator/filter-report helper.

## Target interface

Required inputs:

- `pseudobulk`: `OPENBIO_SINGLE_CELL_PSEUDOBULK` only;
- `population`: exactly one artifact population;
- `condition_key`: normally inherited/locked from the artifact, visible only if the artifact explicitly permits one of several validated Sample-level Condition fields;
- `reference_condition`;
- `comparison_condition`.

Visible optional inputs:

- `technical_batch_key`: an artifact-validated Sample-level field, default inherited/none;
- `categorical_covariate_keys`: a small explicit list of additive nuisance fields;
- `continuous_covariate_keys`: a small explicit list of additive numeric fields;
- `fdr_threshold`, default 0.05, for calls/reporting only;
- `min_abs_log2_fold_change`, default 0, for calls/reporting only.

Advanced expression-filter controls, matching the executed public filter policy:

- `min_count`, `min_total_count`, `large_n`, `min_prop` with documented edgeR/Decoupler semantics.

Fixed hidden policy:

- one population and one row per Sample;
- additive intercept design ordered as Technical batch, declared covariates, then Condition;
- comparison-minus-reference contrast;
- TMM effective-library normalization;
- QL fitting/testing with `robust=True` where the audited adapter can execute it, otherwise direct edgeR calls or an explicit incompatibility failure;
- BH adjustment over all retained tested genes;
- no LFC shrinkage and no top-N truncation of the scientific table;
- no mutation of the input artifact.

Outputs:

- `table`: canonical complete gene result table;
- `summary`: strict JSON text;
- `code`: executable equivalent function source.

## Shared low-coupling preflight seam

EdgeR and PyDESeq2 should call one private, pure preparation module that accepts the typed artifact and explicit roles, then returns:

- copied counts and Sample metadata for one population;
- deterministic Sample/gene identities;
- encoded design matrix with term/column metadata;
- comparison-minus-reference contrast vector;
- design rank, residual df and confounding/estimability diagnostics;
- modeled Condition replicate counts;
- an expression-filter mask and full threshold/universe record.

This helper owns common validity, not model execution. It must not know EdgeR result columns or PyDESeq2 fitting objects. The engines own package adapters, model-specific warnings, canonical table mapping, citations, versions and equivalent source.

The expression filter can use current `decoupler.pp.filter_by_expr` on the exact selected libraries/design grouping, or a validated direct equivalent. Pin/check the supported Decoupler interface because the filter is part of the scientific method. Never reuse a mask computed across all populations or cells.

## Design construction and validation

Role aliases are recorded and warned rather than rejected before design construction; if an alias creates duplicate/aliased design columns, the rank check still fails. After subsetting to the two requested Conditions:

- require distinct represented levels; warn and mark formal interpretation invalid below three Samples in either level, but let full-rank/positive-residual checks determine computability;
- reject duplicate Sample rows, missing/blank metadata, unknown levels, and string-coercion collisions;
- categorical/continuous terms must be complete and correctly typed; constant terms are omitted and warned because they add no estimable column; sample-unique or otherwise saturated terms ordinarily fail the later rank/residual check;
- encode an intercept plus declared additive terms, record categorical reference levels, and force Condition reference to `reference_condition`;
- require rank equal to design-column count and `n_samples - rank > 0`;
- identify exact confounding/aliased columns and reject; identify incomplete overlap/imbalance separately and warn without blocking;
- construct the pairwise Condition vector against the encoded matrix and verify it is nonzero, finite and estimable;
- report design row order, column names, rank, residual df and vector.

This node does not support paired/repeated Sample measurements, donor random effects, interactions, splines, nested technical replicates, or arbitrary continuous-condition contrasts. Fail with an actionable migration message to a future advanced node rather than approximate them.

## Canonical result and report contract

Return one row for every retained/tested gene:

```text
gene: string
log2_fold_change: float
log_counts_per_million: float
quasi_likelihood_f: float
p_value: float
p_adjusted: float
```

Validate unique genes, exact retained-gene set, expected numeric domains (`0 <= p <= 1`, `0 <= p_adjusted <= 1`, non-negative F), and no unaccounted row loss. Stable sorting can use p-value then gene, but the table is never truncated to “significant” genes. Report call thresholds do not affect model fit or BH universe.

The strict summary schema contains `methods`, `results`, `key_results`, `parameters`, `warnings`, `references`, and `software_versions`. In addition to method/version disclosure it must include:

- artifact fingerprint, caller-declared count source/Sample identity, annotation status, `formal_interpretation_invalid`, and invalidity reasons;
- population and exact comparison direction;
- independent Sample counts by Condition and total modeled rows;
- design formula-like rendering, columns, rank, residual df, categorical references and contrast vector;
- genes before/after expression filtering and all filter thresholds;
- TMM normalization-factor/library-size summaries;
- complete null/nonfinite accounting;
- FDR-significant up/down counts under declared call thresholds and bounded top-gene records;
- warnings for low/asymmetric replication, provisional/unknown labels, role aliases, incomplete nuisance overlap, filter sensitivity, observed confounding diagnostics and lack of effect-size shrinkage;
- EdgeR/QL/TMM/BH, Pertpy adapter and pseudobulk-practice citations;
- exact Python, Pertpy, Decoupler, rpy2, R, edgeR and numerical-stack versions.

Report prose says “associated with” rather than causal “caused by,” and never turns cell count into n. If no gene meets the call threshold, say so; do not describe the smallest p-value as significant.

## Equivalent-code design

Generated source should define a portable function accepting the pseudobulk AnnData plus explicit artifact metadata/parameters. It performs the same population/Sample/count checks, deterministic role encoding, rank/estimability checks and expression filter, then executes the same public Pertpy/edgeR bridge or direct public rpy2 edgeR sequence. Every package call uses explicit resolved controls.

Runtime and generated code must agree on:

- included Samples/genes and their order;
- design columns, references and contrast direction/vector;
- filter mask and TMM/QL/robust/BH policy;
- canonical column names and full result rows.

If Pertpy cannot pass a required control to the correct edgeR stage, either use an audited direct edgeR adapter consistently in both paths or disclose the exact wrapper behavior; never claim direct-canonical behavior while executing something else. Generated code may omit plugin history/timing, not scientific validation or result-state checks.

## Dependency and failure policy

Preflight imports and reports components in order: Pertpy, rpy2, R runtime, edgeR and required R dependencies. The error includes detected versions and an install/configuration action. There is no silent fallback to a t-test, Wilcoxon test, Scanpy rank-genes routine, PyDESeq2, or cell-level method because that changes the requested analysis.

## Verification plan for the future implementation batch

- Socket/registration tests accept only typed pseudobulk and reject generic/cell AnnData.
- Multi-population artifacts require one explicit valid population; selected data has one row per Sample.
- Count tests reject floats/logs/negative/non-finite/zero libraries, duplicate Samples/genes and altered artifact fingerprints.
- Design fixtures cover balanced and low/asymmetric replication, Technical-batch adjusted, role aliases, constant and multiple categorical/continuous covariates, absent levels, rank deficiency, zero residual df, exact Condition-batch confounding and non-estimable contrasts.
- A one-versus-two Condition fixture runs when the design has positive residual df, emits low-replication warnings/invalidity, and preserves runtime/generated parity; one-versus-one remains a zero-residual hard failure.
- Direction fixture proves positive log2FC is comparison/reference.
- Filter fixtures prove selection occurs on one population/modeled rows and BH covers exactly retained genes.
- Fake Pertpy/rpy2 adapter asserts DGEList/TMM/dispersion/robust QL/QLF/BH controls and canonical renaming; optional real EdgeR numerical smoke test runs only with an explicitly compatible environment.
- Complete result validation, empty-significance reporting, strict JSON without NaN/Infinity, bounded summaries, citations/version completeness, generated-code compilation and primary table/design/filter parity.
