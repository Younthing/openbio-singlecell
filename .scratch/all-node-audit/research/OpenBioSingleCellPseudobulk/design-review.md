# Pseudobulk aggregation: module design review

## Current module and open-boundary correction

The pre-change typed implementation emitted an immutable pseudobulk artifact plus `summary` and standalone `code`, used Decoupler sum aggregation, and exposed an explicit Dynamic Expression Source, but it over-validated Snapshot history, aligned full-gene Raw, X/layer equality to Raw, Curated formal provenance, global role disjointness, and complete-stratum policy. The corrected implementation removes those authorization gates because they do not prove the corresponding biological claims.

The corrected boundary treats source, feature completeness, annotation status, Sample identity, and metadata roles as analyst declarations. It retains numeric/axis/overflow checks and immutable content fingerprints, reports declaration limitations in strict JSON, and sets `formal_interpretation_invalid` when reporting practice is weak without blocking a computable aggregation.

The overflow postcondition is exact rather than conservative. For each non-negative Sample-by-population profile, the adapter first sums source entries as unbounded Python integers. It rejects only an actual total above signed-int64 maximum, then performs the per-gene int64 sum and verifies that its exact total agrees. This keeps overflow a necessary safety boundary while accepting sparse or zero-heavy profiles whose `maximum * cell_count` bound is a false positive.

The internal exact sum is authoritative because Decoupler 2.2 emits float64 aggregates/QC. Backend cross-checks are exact only while expected integers are at most `2^53`; larger values can be checked only for float64 rounding equivalence and are explicitly counted/warned in `summary`. This precision limitation does not change the exact int64 artifact or exact profile-filter decision.

## Correctness findings from the pre-change audit

### P0: the type system cannot prevent cell-level pseudoreplication

Both this output and the three current contrast inputs use the same generic AnnData socket. A raw single-cell object, normalized object, or unrelated aggregate can therefore be wired directly into a count-model node. Metadata names and matrix appearance are not adequate proof of the statistical unit.

The output must be a concrete `OPENBIO_SINGLE_CELL_PSEUDOBULK` artifact carrying validated count state, Sample/population identity, Condition and Technical-batch roles, aggregation/filter policy, annotation status, and provenance. EdgeR/PyDESeq2 must consume that type only.

### P0: mean aggregation is incompatible with the advertised count-model family

Mean profiles are not raw integer library counts. They invalidate the negative-binomial count assumptions checked by EdgeR and PyDESeq2. A generic summary-statistic node could expose mean in the future, but this formal pseudobulk node must fix aggregation to sum and remove `mode` from its public scientific interface.

### Sample and Technical batch integrity is not enforced

The current node does not require a Condition key, does not verify Condition is constant within Sample, and has no separate Technical batch role. Pertpy preserves an obs field only when it happens to be constant within an aggregate, so an inconsistent field can disappear instead of producing a scientific error. The adapter must preflight these invariants and name every conflicting Sample.

Technical batch must not be concatenated with Sample to create extra observations. It is preserved once per Sample for downstream adjustment. If Technical batch varies inside Sample, fail this simple design rather than manufacturing technical profiles as biological replicates.

### The expression-state contract must be explicit without pretending to prove biology

The node resolves one explicit Raw/X/layer source and validates finite, non-negative, integer-like counts on that source's own aligned feature axis. It must not require history, Raw equality, or a current/full-gene match: these are neither algorithmic prerequisites nor reliable evidence of biological completeness. The artifact records `caller_declared_raw_counts`, the selected feature axis, and that provenance/completeness were not software-verified.

### Filtering and missing profiles are under-disclosed

`min_cells` and `min_counts` are defensible visible QC choices, but the present output does not say which Sample-population profiles were removed or how Condition replication changed. A missing Sample-by-population profile is not equivalent to an all-zero library and must not be synthesized. Every retained/removed profile and before/after replicate count belongs in `summary`.

### The adapter diverges from the current official family example

Pertpy's DE tutorial now uses `decoupler.pp.pseudobulk` plus `decoupler.pp.filter_samples`; `PseudobulkSpace` is a generic perturbation-space utility. Decoupler 2.2.0 itself has a QC field-name mismatch between docstring and implementation. A private adapter can handle that incompatibility and emit stable canonical fields; the inconsistency must not leak into node wires.

### No report or equivalent source

The node returns only data. It lacks method references, replicate/profile disclosure, software versions, a report-ready aggregation summary, and code that a researcher can rerun outside ComfyUI.

## Decision: keep and deepen as the only Sample-by-population count aggregator

Retain and enhance `OpenBioSingleCellPseudobulk`. It should perform exactly one atomic operation: validate a declared raw-count source, establish Sample/metadata integrity, sum cells into Sample-by-population profiles, apply profile-level QC, and package the result as a typed immutable scientific artifact.

Do not merge it with EdgeR or PyDESeq2. Aggregation is shared by multiple count engines, is useful for QC before choosing an engine, and has distinct failure/reporting semantics. Conversely, do not keep the redundant aggregation inside `OpenBioSingleCellDecouplerPseudobulkContrast`; migrate that legacy node to this explicit seam.

## Target interface

Required inputs:

- `adata`: cell-level AnnData;
- `sample_key`: independent biological Sample identity;
- `population_key`: cell population used to define within-population profiles;
- `condition_key`: Sample-level Condition metadata.

Visible optional/advanced inputs:

- `count_source`: explicit Raw, `X`, or layer selector; numeric validation does not infer or upgrade biological provenance;
- `technical_batch_key`: optional, separate Sample-level nuisance variable;
- `annotation_status`: detected from upstream metadata where available and otherwise `unknown`; it controls disclosure, not execution authorization;
- `min_cells`: visible dataset-dependent profile threshold;
- `min_counts`: visible dataset-dependent profile-library threshold.

Fixed hidden policy:

- `mode="sum"`;
- finite, non-negative integer-like selected values and no `skip_checks` escape hatch; suitability as biological raw counts remains caller-declared;
- no empty-grid expansion or synthetic zero profiles;
- exact selected-source feature alignment and no feature selection performed by this node; Raw and current axes need not match;
- deterministic profile order and canonical QC field names;
- no mutation of the input AnnData.

Outputs:

- `pseudobulk`: concrete `OPENBIO_SINGLE_CELL_PSEUDOBULK` artifact;
- `summary`: strict JSON text;
- `code`: executable equivalent Python function source.

## Typed artifact boundary

The artifact should wrap a private AnnData rather than expose a bare mutable object as proof of validity. At minimum its immutable metadata records:

- schema version and fingerprint;
- Sample, population, Condition and optional Technical-batch keys and role definitions;
- annotation status (`Curated`, `Provisional`, or `Unknown`) and a machine-readable interpretation-validity flag;
- caller-declared count source/expression state, selected feature axis, and explicit `history_required=false`/`raw_equivalence_required=false` disclosure;
- fixed sum mode and profile filters;
- original cell/gene/Sample totals;
- retained and removed profile identities/reasons;
- deterministic mapping from each profile to exactly one Sample and population.

The contained AnnData uses one row per retained Sample-by-population profile, `X` as integer summed counts, stable canonical `obs` columns such as `openbio_n_cells` and `openbio_total_counts`, and the original gene axis. Upstream-specific `psbulk_*` names may be preserved privately but are not the public contract.

The artifact is justified by the deletion test: without it, every engine must rediscover whether a generic matrix is really summed Sample-level counts, and the dangerous generic-AnnData wire returns. It is not a generic Python-object/model socket.

## Invariants and failure behavior

Preflight before external package construction:

- materialize a view as a copy; reject backed/empty objects and duplicate obs/var identifiers;
- resolve exactly one matrix and validate numeric, finite, non-negative, integer-like values, nonzero total, shape, and its own exact feature alignment;
- validate Sample/population/Condition labels as non-missing and non-blank before string conversion, detecting conversion collisions;
- require one Condition and at most one Technical batch per Sample; identify all offending Samples;
- reject Samples represented by no cells implicitly rather than fabricate them;
- allow formal or exploratory aggregation for any declared annotation status; non-Curated formal use is marked `formal_interpretation_invalid` with a prominent immutable warning;
- allow role aliases that do not make aggregation structurally ambiguous, disclose them, and defer actual design non-estimability to the downstream engine. Reserved QC-column collisions remain errors.

Postconditions:

- exactly one unique row per represented Sample-population pair;
- all aggregate values equal exact columnwise sums of contributing cells and remain finite integer counts;
- stable gene identity/order and nonzero cell/library totals;
- retained/removed profile sets exactly match the configured thresholds;
- no population/Condition stratum silently loses all profiles;
- report sparse Condition/population strata and fewer than three Samples as fragile without asserting a universal engine minimum; downstream rank/residual checks decide estimability.

Failures are actionable and occur before a downstream engine sees data. Do not coerce metadata, round non-integer values, repair Sample conflicts, substitute missing categories, or treat cell count as replicate count.

## Summary contract

Use the repository report schema: `methods`, `results`, `key_results`, `parameters`, `warnings`, `references`, and `software_versions`. Required payload includes:

- source matrix, selected feature axis, numeric count validation, and caller-declared biological count state;
- input/output dimensions and unique Sample count;
- per-population/per-Condition Sample counts before and after QC;
- contributing cell and total-count distributions;
- bounded retained/removed profile records with complete aggregate counts;
- annotation status, requested mode, `formal_interpretation_invalid`, and machine-readable invalidity reasons;
- explicit statement that Sample, not cell, is the statistical unit;
- min-cell/min-count values and their dataset-specific nature;
- Condition/Technical-batch integrity checks;
- Decoupler, general pseudobulk-practice, and AnnData references;
- runtime openbio-singlecell, Decoupler, AnnData, NumPy, pandas and SciPy versions.

Core prose may say: “Raw counts from N cells were summed into M Sample-by-population profiles across S independent biological Samples; profile QC removed R profiles. This step created count-model inputs and did not itself test Condition effects.” It must not claim a differential result.

## Equivalent-code design

The `code` endpoint emits one self-contained function using public AnnData/Decoupler interfaces. It must:

1. resolve and copy the exact Raw/X/layer source without requiring Snapshot history or equality to Raw;
2. run the same identifiers/count/metadata checks;
3. call `dc.pp.pseudobulk` with explicit `sample_col`, `groups_col`, `layer`, `raw=False`, `empty=False`, `mode="sum"`, and `skip_checks=False`;
4. accept either Decoupler 2.2.0 cell-QC spelling internally and rename to the stable OpenBio field;
5. apply identical profile QC and postconditions;
6. construct and return the same typed scientific artifact or, if artifact classes are intentionally unavailable outside the plugin, an explicitly documented portable AnnData-plus-metadata representation with identical primary state.

Do not generate a simplified `groupby().sum()` snippet that omits sparse/count validation, metadata integrity, removed-profile accounting, or feature alignment. Plugin execution history and timing may be omitted; all scientific defaults and versions remain visible.

## Cohesion and coupling assessment

The deep seam is `cell AnnData -> validated typed Sample-by-population counts`. Internal helpers may be shared with EdgeR/PyDESeq2 for JSON-safe reporting and label validation, but engine design/preflight does not belong here. Decoupler version quirks are local-substitutable implementation details hidden behind one adapter.

This gives high cohesion: every rule establishes or verifies aggregate profile identity. It lowers coupling because engines consume one stable contract rather than knowledge of Pertpy/Decoupler aggregation fields, count layers, and Sample semantics.

## Verification plan for the future implementation batch

- Dense, CSR and CSC exact integer sums from X and a named layer; source selection must be demonstrable.
- Accept independently valid Raw, X, and named-layer count sources, including Raw/current feature-axis differences; reject negative/non-finite/non-integer/empty data, duplicate selected-axis identifiers, blank labels, and string-collision labels.
- Condition and Technical-batch conflict fixtures name all conflicting Samples and never create extra replicates.
- Missing Sample-population combinations remain absent; no empty zero profiles are created.
- Profile thresholds remove exactly the expected rows and disclose before/after Condition replication.
- Curated, Provisional, and Unknown status remain executable; formal interpretation validity/reasons are explicit and immutable.
- Snapshot history and X/layer equality to Raw are ignored and are not boundary tests.
- Typed-socket registration prevents generic cell AnnData wiring to EdgeR/PyDESeq2.
- Fake Decoupler variants emitting `psbulk_cells` or documented `psbulk_n_cells` normalize to one contract; unsupported API fails with installed version.
- Input immutability, deterministic row/gene order, strict JSON without NaN/Infinity, complete citations/versions, generated-code compilation, and runtime/code primary-state parity.
- Runtime and generated source both accept a non-negative zero-heavy profile whose exact total fits signed int64 even when `maximum * cell_count` does not, and both reject a profile whose exact total truly exceeds signed int64. These are source-selection tests only; they do not assert Raw/current equality or analysis-history provenance.
