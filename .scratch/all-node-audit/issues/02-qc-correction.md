# QC decisions and artifact correction

Type: task
Status: resolved

Blocked by: 01

## Scope

- Mark MAD Outliers
- Run Scrublet
- Filter Predicted Doublets

## Initial audit decisions

- Correct SciPy MAD normal scaling and make per-Sample thresholds/reporting explicit.
- Require an explicit count expression source for Scrublet and use Sample/capture as the default grouping concept.
- Keep doublet prediction and destructive filtering as separate atomic nodes.
- Reject ambiguous string/object prediction columns rather than coercing them with `astype(bool)`.

## Comments

- 2026-08-28: Initial audit complete; official-usage and design documents pending before code changes.
- 2026-08-28: Landed separate official/method usage research and design-review records for all three nodes before implementation. Decisions: enhance all three; retain annotation, stochastic detection, and destructive filtering as separate atomic boundaries.
- 2026-08-28: Implemented strict count-source and boolean contracts, per-Sample MAD/Scrublet reporting, strict JSON summaries, equivalent source outputs, dependency disclosure, and legacy widget-order regression coverage. Ruff passed and `tests/test_correction.py -W error` completed with 14 passing tests.

## Answer

All three nodes remain separate atomic operations. MAD thresholds now use the documented SciPy scaling direction and disclose per-group thresholds; Scrublet consumes an explicit integer count source and fails clearly when a group cannot support the requested PCA or automatic threshold; filtering accepts only complete boolean predictions. Each node now emits its primary result plus strict `summary` and equivalent `code` outputs, with references, versions, warnings, limitations, and group-level key results.
