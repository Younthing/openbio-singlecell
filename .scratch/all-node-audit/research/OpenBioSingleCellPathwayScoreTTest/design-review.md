# OpenBioSingleCellPathwayScoreTTest — design review

## Decision

**Replace the implementation in place and retain the ID as a Sample-level pathway-score Welch contrast.** The old
cell-level behavior is P0 and must be deleted, not compatibility-preserved. The name still describes the corrected
method, but the interface is intentionally breaking.

## Atomic interface

One execution tests one population and one two-Condition contrast across the complete pathway-score family.

Inputs, in order:

1. `adata` carrying a canonical, provenance-bound score DataFrame;
2. `sample_key`, default `sample`;
3. `condition_key`, default `condition`;
4. `annotation_key`, default `cell_type`;
5. required `population`, `condition_a`, and `condition_b`;
6. `score_key` for a named DataFrame in `obsm`;
7. `min_cells_per_sample_population`, default 10;
8. `min_samples_per_condition`, default 3 with hard minimum 2; two-Sample arms remain executable but warning-grade;
9. optional advanced `technical_batch_key` used only for audit;
10. advanced `confidence_level`, default 0.95;
11. `annotation_status` (`unknown`, `provisional`, `curated`) for disclosure.

Outputs:

```text
table, summary, code
```

Hidden fixed choices are arithmetic Sample means, two-sided Welch test, `nan_policy="raise"`, mean difference
`condition_a - condition_b`, and BH across all tested pathway columns. No equal-variance switch, per-cell mode,
one-sided alternative, arbitrary p-filter, or Technical batch adjustment is exposed.

## Canonical result table

One deterministic row per pathway:

```text
population, pathway, condition_a, condition_b,
n_samples_a, n_samples_b, n_cells_a, n_cells_b,
mean_a, mean_b, sd_a, sd_b, mean_difference,
ci_low, ci_high, t_statistic, degrees_of_freedom,
p_value, p_adjusted, rank
```

All numeric outputs are finite and family-complete. Pathway columns with non-finite input or zero variance that
cannot yield a finite test result cause a bounded hard failure rather than silent family shrinkage. Ties share rank.

## Provenance and summary/code contract

Validate unique observation IDs; nonblank scalar Sample/Condition/population labels; Sample-to-Condition mapping;
audit rather than reject Sample-to-Technical-batch multiplicity and Condition/batch confounding; exact selected
Conditions; score index equality; finite score matrix; score
producer/resource/method fingerprints; and no output truncation. Pin the ordered observations, selected cells,
Sample aggregation frame, and complete pathway family in fingerprints.

`summary` includes methods/results prose suitable for a report, per-arm support, exclusion reasons, top effects and
alternatives, score provenance, Technical batch audit, warnings/limitations, Welch/BH/upstream-method references,
and dynamic package versions. Strict JSON forbids NaN/Infinity. `code` defines an equivalent function returning
`(DataFrame, summary_dict)` and shares the canonical aggregation/statistics/summary logic.

## Migration and tests

Migration must require new Sample, population, and Condition choices; old workflows cannot be auto-converted because
their `group_column` may be a cell-level label and they lack a Sample key. Preserve links only for inspection and add
an actionable migration note.

Tests: hand-calculated Sample means/Welch df/CI/BH; strong within-Sample correlation fixture demonstrating cell-level
inflation; Sample-to-Condition conflicts; two-Sample warning and truly insufficient Samples/cells; missing
populations; categorical unused levels; score-index/provenance tampering; non-finite/constant pathways; Technical
batch multiplicity/confounding warning plus formal-invalid status; complete-family BH;
annotation-status disclosure; deterministic ordering/ties; immutability; strict JSON; compiled code; and exact
runtime/generated table and summary equality.

## Cohesion assessment

The revised module owns one Condition contrast on an already computed score family. Upstream score computation and
resource choice remain separate, while aggregation, replicate validation, Welch inference, correction, and reporting
are localized behind one scientific interface.
