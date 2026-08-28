# Deepen the canonical report contract

Type: task
Status: resolved

## Question

How can one small report Interface validate every analysis summary while preserving method-specific detail and avoiding
circular imports or a rewrite of all callers?

## Scope

- `SummaryResult`, `make_summary_result`, and `make_analysis_report`
- Milo differential abundance
- scCODA and tascCODA differential composition
- scVI differential expression

## Evidence

- `.scratch/all-node-audit/research/OpenBioSingleCellMiloDifferentialAbundance/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellMiloDifferentialAbundance/design-review.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellSccodaDifferentialComposition/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellSccodaDifferentialComposition/design-review.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellTasccodaDifferentialComposition/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellTasccodaDifferentialComposition/design-review.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellSCVIDifferentialExpression/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellSCVIDifferentialExpression/design-review.md`

## Acceptance criteria

- The shared seam rejects missing or wrongly typed canonical fields.
- Method-specific fields may remain as additive JSON-compatible detail.
- All existing analysis/report nodes cross the same validation seam.
- The four known drifting node implementations return the canonical types and required version disclosure.

## Comments

- 2026-08-29: Claimed first because the current `Any`-typed seam permits semantic drift even when every port and JSON
  serialization test passes.
- 2026-08-29: Added the dependency-free `report_contract` Module behind `SummaryResult`; it validates the ten canonical
  fields, field types, references, mandatory Python/OpenBio versions, and strict JSON while allowing scientific detail.
- 2026-08-29: Repaired Milo, scCODA/tascCODA, and scVI-DE without compatibility aliases. Targeted affected suites pass
  (`162 passed, 1 skipped`); the first full run found only two stale internal callers, both updated at the canonical seam.

## Answer

Keep the existing `SummaryResult` Interface and deepen it with one package-internal, standard-library-only validator.
All analysis/report nodes now receive the invariant automatically, while builders retain responsibility only for
normalizing scientific values and composing reports. Four drifting summaries now use string `methods`/`results`, a
mapping `key_results`, complete references and runtime versions, with method-specific evidence preserved as additive
fields rather than alternate schema branches.
