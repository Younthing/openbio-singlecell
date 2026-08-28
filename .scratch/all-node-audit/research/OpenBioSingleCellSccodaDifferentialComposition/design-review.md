# OpenBioSingleCellSccodaDifferentialComposition: pre-change design review

## Decision

**Retain and substantially enhance this node as the flat scCODA Sample-level differential-composition module. Do not merge it with tascCODA, a descriptive composition summary, a generic proportion test, marker analysis, or cell-level differential expression. Do not delete it.**

The atomic scientific question is:

> For one explicitly directed two-condition comparison, do any annotated cell types have a credible relative compositional change across biological Samples under the flat scCODA model, optionally adjusting for declared Sample-level nuisance covariates?

scCODA and tascCODA may share one private adapter for validated cell-level AnnData to Sample-by-cell-type counts, deterministic design coding, reference resolution, NUTS execution, diagnostics, versions, and strict JSON. They must remain separate public nodes. tascCODA's hierarchy changes the parameterization, selection prior, result rows, tuning requirements, interpretation, dependencies, and migration safety. A `method = flat | hierarchical` switch would hide those differences and make the module shallow.

The deletion test favors retention: removing the scCODA module would force every caller to reimplement Sample aggregation, compositional reference semantics, posterior expected-FDR selection, NUTS disclosure, and canonical effect reporting. Those are substantial, reusable hidden decisions.

Repository constraints applied from `CONTEXT.md` and `docs/agents/domain.md` are: Sample means biological replicate, Condition contrast is Sample-level and replicate-aware, and Technical batch is not silently promoted to replication. ADR `docs/adr/0001-separate-marker-evidence-from-condition-inference.md` also requires marker evidence to remain separate from formal Condition inference. The local issue-tracker instructions were read for process context; this research-only task intentionally creates no issue or ledger entry.

## Problems in the current seam

The current wrapper exposes upstream mechanics while omitting the scientific contract:

- `covariate_keys` plus arbitrary `formula` can describe multiple effects or disagree, so a run is not one atomic comparison;
- no input identifies reference and comparison Condition levels, making direction data-order dependent;
- Pertpy can silently skip a covariate that varies within Sample;
- no check establishes that Sample is the biological replicate or that each Sample maps to one Condition;
- no minimum replication, extra-level, design-rank, confounding, count-table, or reference-candidate validation exists;
- intercept, focal effects, nuisance effects, and possibly node effects are concatenated into a heterogeneous table;
- the only output omits the posterior selection threshold, diagnostics, method references, software versions, limitations, and executable equivalent code;
- a successful one-chain Pertpy fit could be mistaken for proven convergence;
- dependency behavior is unpinned even though OpenBio does not currently declare the Pertpy/JAX stack.

This is low information hiding: a user must understand Pertpy loader logs, Patsy coding, scCODA selection internals, and result-frame structure to use a small-looking node correctly.

## Proposed public interface

Recommended outputs:

```text
table, summary, code
```

Recommended inputs, conceptually:

| Input | Visibility | Contract |
|---|---|---|
| `adata` | primary | Cell-level AnnData. It is read-only; implementation works on a copy. |
| `sample_key` | primary | One canonical biological experimental-unit identifier. No composite string joining inside Pertpy. |
| `annotation_key` | primary | Categorical cell-type annotation used to count cells. Annotation provenance/status is disclosed. |
| `annotation_status` | primary | Required `curated` or `provisional`; provisional labels make the result explicitly exploratory. |
| `condition_key` | primary | Sample-constant Condition column. |
| `reference_condition` | primary | Exact baseline level. |
| `comparison_condition` | primary | Exact compared level; table direction is comparison relative to reference. |
| `reference_cell_type` | primary | Exact annotation label or `automatic`; the resolved reference is always reported. |
| `adjustment_covariates` | advanced | Optional Sample-level nuisance keys, default empty, with deterministic coding and disclosed bases. No arbitrary transformations/interactions. |
| `estimated_fdr` | primary | Posterior expected-FDR target; default `0.05`, strictly between 0 and 1. |
| `num_samples` | advanced | NUTS posterior draws; default `10000`, positive. |
| `num_warmup` | advanced | NUTS adaptation draws; default `1000`, positive under project policy. |
| `random_seed` | advanced | Nonnegative integer, default `0` to match pinned Pertpy; one seed controls initialization/JAX sampling and is disclosed. |

`reference_condition` and `comparison_condition` must not default silently to lexical or observed order. If UI defaults are necessary, they should remain unresolved sentinels that force a user choice after inspecting available levels.

For adjustment covariates, a comma-separated widget may remain a transport format, but the module contract is an ordered list of unique keys. Numeric covariates must be finite. Categorical nuisance variables should be ordered categoricals with their declared first category as base, or use an explicit structured coding policy; fallback to first observed Sample is forbidden. The summary records all generated design columns and bases.

## Parameters to hide

The following are implementation policy, not ordinary public choices:

- Pertpy cell-level mode, `generate_sample_level=True`, modality names `"rna"`/`"coda"`, and internal count-table layout;
- Patsy formula text and safe internal column names;
- fixed intercept handling and treatment coding implementation;
- the zero pseudocount of `0.5` used by pinned Pertpy;
- automatic-reference absence threshold `0.05`;
- NUTS rather than a sampler switch;
- NumPyro's pinned kernel defaults such as target acceptance and maximum tree depth, unless a later diagnostic-specific design separately reviews them;
- backend plumbing, tree keys, upstream frame names, and table-normalization mechanics.

Hidden does not mean undisclosed. Every result-defining resolved value appears in `summary` and generated `code`. If users later need interactions, continuous focal effects, multi-level omnibus inference, multi-chain custom NumPyro execution, or a different sampler, those are separate reviewed capabilities rather than strings smuggled through `formula`.

## Deep execution seam

```text
validate AnnData axes and categorical annotation
  -> materialize one Sample metadata row per biological unit
  -> validate Sample-constant Condition/adjustments and exactly two selected levels
  -> aggregate integer cell counts by Sample x cell type
  -> validate replication, totals, cell-type support, reference, and design rank
  -> create collision-safe internal design and exact focal coefficient
  -> run version-gated Pertpy 1.3.0 scCODA NUTS on a copy
  -> extract posterior selection, focal effect rows, and available diagnostics
  -> enforce finite/stable output postconditions
  -> construct canonical table, strict summary JSON, and equivalent source atomically
```

The private composition adapter is worthwhile because it hides a changing third-party MuData API and makes validation testable with a deterministic fake. It must not become a public plug-in socket or accept arbitrary model objects; there is only one supported versioned scCODA path.

## Input and design preconditions

All inexpensive checks occur before importing/initializing JAX or fitting the model.

### AnnData and count construction

- nonempty observations and variables; unique observation identifiers;
- required `.obs` keys exist and have nonmissing, finite/string-safe values;
- `annotation_key` has at least two modeled cell types and no missing labels;
- `annotation_status` is explicit; `provisional` is allowed for exploration but forces `inference_scope="exploratory"` and provisional-label language in every result disclosure;
- each cell belongs to exactly one canonical Sample;
- grouping produces a nonnegative integer Sample-by-cell-type count table with a positive total for every Sample;
- the table preserves deterministic Sample and cell-type ordering and records a fingerprint;
- no silent removal of rare cell types or empty Samples; any explicit upstream subset is outside this node and is reported by its workflow provenance.

The input expression matrix is not used. This node counts annotation labels; its `input_genes` metadata is incidental and should not appear as if genes were modeled.

### Experimental-unit and contrast checks

- every Sample maps to exactly one Condition and one value for every adjustment covariate;
- selected reference and comparison levels both exist and are distinct;
- no third Condition level enters the fit; fail with an upstream-subset instruction;
- at least two biological Samples per compared level are required by OpenBio policy for a formal contrast; fewer is not rescued by more cells;
- below a documented stronger replication threshold, emit a prominent low-replication limitation rather than claiming adequate power;
- technical replicates are not counted as independent Samples unless upstream aggregation defines the biological unit correctly;
- repeated/paired observations are rejected unless the final implementation documents a valid supported fixed-effect design; scCODA here has no random-effect contract.

The two-Sample minimum is an OpenBio reporting safeguard derived from the repository's Sample-level domain model, not a claim that the original Bayesian method guarantees validity at that boundary.

### Design checks

- build one explicit focal treatment coefficient for comparison versus reference;
- include only declared adjustments; no interactions or data-dependent transformations;
- verify finite design values, nonzero variance where applicable, matrix dimensions, rank, and condition number;
- reject perfect confounding/rank deficiency and any inability to map exactly one upstream effect row back to the focal contrast;
- disclose deterministic categorical bases and continuous covariate scaling; preferably require upstream scaling rather than hide it here.

### Reference checks

- explicit reference label exists in the modeled cell-type axis and is not empty across all Samples;
- automatic reference uses the pinned `0.05` zero-fraction threshold, records all eligible candidates and dispersion values, and fails when none is eligible;
- warn that automatic selection is data-dependent and does not establish biological invariance.

## Canonical `table` contract

The primary table contains one row per modeled cell type for the exact focal Condition coefficient. It does not mix intercept or nuisance rows. Suggested stable snake-case fields are:

| Field | Meaning |
|---|---|
| `contrast` | Canonical `comparison_condition vs reference_condition` label. |
| `condition_key` | Source Condition field. |
| `reference_condition` / `comparison_condition` | Exact selected levels. |
| `cell_type` | Exact modeled annotation label. |
| `reference_cell_type` | Resolved compositional reference. |
| `model_coefficient` | scCODA final relative coefficient after posterior selection. |
| `hdi_lower` / `hdi_upper` | Pertpy-reported conditional HDI bounds. |
| `posterior_sd` | Pertpy-reported posterior SD. |
| `inclusion_probability` | Posterior inclusion probability from pinned Pertpy. |
| `credible_effect` | Boolean under the realized posterior expected-FDR threshold. |
| `estimated_fdr` | Requested posterior expected-FDR target. |
| `inclusion_probability_threshold` | Realized threshold used for selection. |
| `expected_count_reference` / `expected_count_comparison` | Model-derived compositional expected counts at the documented baseline covariate profile. |
| `compositional_log2_fold_change` | Renormalized relative change; not an absolute abundance fold change. |

If Pertpy exposes only one `Expected Sample` per effect row, the adapter must derive and validate the two named expected-count quantities from the pinned model parameters or narrow the schema honestly. It must not relabel one ambiguous field as both groups.

The reference row remains present with a fixed coefficient and explicit `reference_constraint=true` if that helps downstream joins. Missing/undefined numeric values use nullable columns, not strings or NaN hidden inside object dtype. Ordering is deterministic: credible effects first only if explicitly documented, otherwise annotation category order.

Intercept and adjustment-effect tables may be embedded under `summary.model_details` as JSON-safe records or omitted from the public artifact when not needed for the focal report. They must not be presented as additional discoveries.

## `summary` JSON contract

`summary` is valid strict JSON: no NaN/Infinity, NumPy/pandas scalars, tuples used as object keys, or backend objects. Suggested top-level content:

- `schema_version`, node/method identifiers, status, and UTC completion metadata;
- `question`: Condition key, exact direction, adjustment set, and statement that Sample is the unit;
- `input`: cell count, Sample count, per-group Samples, cell types, count-table shape/fingerprint, Sample total distribution, zero fraction, annotation key/status/provenance;
- `design`: internal generated formula representation, exact columns, focal column, categorical bases, dimensions, rank, condition number, and validation results;
- `reference`: requested/resolved cell type, automatic threshold, candidate diagnostics, zero fraction, and dispersion;
- `model`: Dirichlet-multinomial scCODA, zero pseudocount, NUTS settings, seed, backend/device/x64, posterior expected-FDR target, and realized threshold;
- `selection`: exact pinned cutoff rule and realized posterior expected-FDR quantity recomputed after Pertpy's three-decimal threshold flooring;
- `diagnostics`: chain count, official mean-acceptance check, within-chain ESS/MCSE where calculable, potential-energy/step-size/step-count summaries, finite checks, `divergences_available=false`, and a status that never asserts multi-chain convergence;
- `results`: row count, credible-effect count, ordered key-result records, strongest positive/negative relative effects, and explicit coefficient versus compositional-LFC distinction;
- `warnings` and `limitations`: low replication, provisional annotation, automatic reference, single-chain evidence, compositional/relative interpretation, no absolute abundance, no cell-level replication, no random effects;
- `software`: exact package/backend versions;
- `references`: machine-readable list of method/package citations and official URLs.

Writing-ready prose must be conditional on diagnostics. For example: “Under scCODA 1.3.0 with B relative to A and cell type R as reference, K cell types met the posterior expected-FDR 0.05 credibility rule.” It must not use “significant” as a substitute, claim convergence, or convert a relative composition result into an absolute cell-count conclusion.

## Equivalent `code` contract

The text output defines one importable function that accepts the same scientific inputs and returns the same canonical table plus summary dictionary. It must:

1. perform the same Sample/covariate/annotation/count/rank/reference preflight;
2. build the same collision-safe exact contrast rather than accept an arbitrary formula;
3. pin or check the same Pertpy contract and work on `adata.copy()`;
4. call `Sccoda.load`, `prepare`, `run_nuts`, and expected-FDR selection with every hidden parameter explicit;
5. extract and validate the same focal row and diagnostic fields;
6. construct the same stable table schema and strict summary content;
7. include method/package references and actual software versions.

Equivalent means scientifically reproducible, not merely a four-line upstream example. Wall-clock timing and plugin-only transport metadata may differ. The function should fail under the same malformed inputs/backend outputs as runtime code, and it must not mutate the caller's AnnData on failure.

## Diagnostics and reporting policy

Pertpy 1.3.0 uses one NumPyro chain by default. The result can pass numeric postconditions but cannot establish between-chain convergence. Proposed status meanings:

- `limited_single_chain`: finite fit, Pertpy mean acceptance in `[0.6, 0.95]`, and no declared within-chain diagnostic warning; still explicitly limited;
- `warning`: fit returned but official acceptance warning, poor/undefined ESS/MCSE, extreme tree-step/energy behavior, low replication, or another disclosed concern exists;
- `failed`: nonfinite posterior/output, missing focal coefficient, malformed shapes, or an upstream error; emit no partial result bundle.

Do not invent a zero divergence count: the pinned runner does not retain that field. If a future version captures divergences or supports tested independent chains, it requires a contract/version update and fixtures.

`summary_prepare(kind="stats")` does not by itself return ESS or R-hat. If the implementation uses `make_arviz`/ArviZ for diagnostics, test exact chain/draw dimensions and preserve the fact that R-hat is undefined for one chain. Numeric diagnostic thresholds beyond Pertpy's official acceptance policy must be documented as OpenBio policy, versioned, and reported.

## Boundary and failure semantics

Preflight failures include missing/duplicate axes; missing/multi-valued Sample metadata; invalid Sample definition; fewer than two conditions/cell types/Samples per group; extra Condition levels; missing or nonfinite adjustment values; rank deficiency/confounding; malformed counts; invalid reference; no automatic-reference candidate; invalid FDR/draw/warm-up/seed; missing/incompatible Pertpy extras; and unsupported annotation provenance under the chosen formal-inference policy.

Post-fit failures include absent modalities, altered Sample/cell-type axes, missing focal coefficient, duplicated rows, nonfinite posterior statistics, impossible probabilities/thresholds, malformed diagnostic arrays, or a table/summary/code bundle that cannot be serialized atomically.

All validation and model work occurs on temporary objects. Outputs are committed only after table, summary, and code all pass postconditions. Errors name the failed invariant, observed Sample/group/reference values where safe, installed Pertpy version, and an actionable correction.

## Workflow migration boundary

The legacy schema cannot generally be migrated automatically:

- arbitrary `formula` and `covariate_keys` do not identify exact reference/comparison levels;
- the legacy schema does not record whether the counted annotation is curated or provisional;
- implicit Patsy bases may depend on ordering;
- the meaning of the old heterogeneous table changes to a focal-only canonical table;
- two outputs are added;
- connected/exposed legacy widgets can carry dynamic formulas or covariate sets that cannot be proved equivalent.

Recommended migration behavior is full-graph atomic preflight followed by an actionable rejection for every genuine legacy node: recreate it with explicit `condition_key`, `reference_condition`, `comparison_condition`, optional adjustments, and review downstream table consumers. A very narrow automatic migration would be acceptable only if serialized data proves one categorical covariate, both exact levels and their direction, a stable reference, no formula exposure/connection, and compatible table consumers; typical legacy workflows do not contain enough information, so rejection is safer.

Current-schema nodes are no-ops. Partial output bundles, mixed old/new widgets, divergent named/positional values, malformed input/output endpoints, and proxy/direct exposure conflicts reject before any graph mutation. Root and subgraph traversal, object/array links, global link-ID allocation, and idempotence follow repository migration policy, but implementation belongs in the later migration phase, not this research change.

## Verification plan

### Unit and schema tests

- exact input/output names, order, defaults, visibility, and registration;
- deterministic contrast coding for string, categorical, and safe escaped names/levels;
- Sample-by-cell-type aggregation and ordering on dense/sparse AnnData;
- Sample uniqueness, missing metadata, extra levels, replication boundaries, and technical-unit misuse messages;
- categorical adjustment bases, numeric adjustments, rank/condition-number reporting, and confounding rejection;
- explicit and automatic cell-type references, candidate threshold, and no-candidate failure;
- exact pinned posterior expected-FDR cutoff semantics, including strict inequality, three-decimal flooring, no-candidate fallback, and stable canonical columns;
- coefficient versus compositional-LFC disclosure, including a zero coefficient with nonzero renormalized LFC;
- strict JSON and complete version/reference lists;
- generated code compiles and matches runtime validation/table/summary on deterministic fake results;
- caller AnnData remains unchanged after success and every failure.

### Adapter and diagnostic tests

- a fake Pertpy adapter verifies exact `load`/`prepare`/`run_nuts` arguments and returned MuData axes;
- requested covariate cannot be silently skipped;
- focal coefficient mapping is unique and nuisance/intercept rows never leak into the primary table;
- one-chain status, undefined R-hat, available ESS/MCSE, acceptance warnings, missing divergence field, nonfinite posterior, and malformed extra-field shapes;
- dependency/version gate and clear installation error for missing `pertpy[jax]`.

### Integration and migration tests

- version-pinned CPU smoke test with very short NUTS settings, marked/skipped when the optional stack is absent; scientific defaults remain 10,000/1,000;
- strict legacy rejection leaves root/subgraphs byte-equivalent;
- connected/exposed formula/covariate widgets, partial outputs, malformed links, proxy exposures, and mixed schemas reject atomically;
- current schema is zero-migration and repeated migration is idempotent;
- example workflow demonstrates Sample-level Condition wiring and consumes `table`, `summary`, and `code` appropriately.

## Acceptance criteria

The refactor is scientifically and architecturally complete only when:

1. one node run represents exactly one directed Sample-level Condition contrast;
2. Sample/covariate uniqueness and design rank are proven before inference;
3. the primary table contains only the focal flat-scCODA effect rows with stable semantics;
4. posterior expected-FDR selection is reported accurately as credibility, not frequentist significance;
5. single-chain diagnostic limitations and unavailable divergences cannot be hidden;
6. all result-defining inputs, hidden fixed policy, references, and actual versions appear in strict `summary` and equivalent `code`;
7. tascCODA remains a separate hierarchy-aware node, with only a private validated-count adapter shared;
8. legacy ambiguity is rejected atomically rather than guessed.

The official and primary-source evidence for these decisions is recorded in [`official-usage.md`](./official-usage.md), especially the pinned Pertpy loader, design, NUTS, selection, and output-field behavior.
