# Integrate the complete plotting surface

Type: task
Status: resolved
Blocked by: 02, 03, 04, 05, 06, 07, 08, 09, 10

## Scope

- Registry counts/categories/output contracts.
- Representative example workflows, README, changelog, release metadata, and visual QA.
- Complete regression and final two-axis review.

## Acceptance criteria

- Every required node is registered exactly once and every worker operation matches.
- Examples and docs use current schemas; all validation commands pass; review has no unresolved findings.

## Comments

- 2026-09-04: Claimed as the final integration ticket.
- 2026-09-04: Central registry/type subtask resolved at 144 unique nodes. Exact 27-category counts, all standard
  Plot `plot/summary/code` outputs, Milo/composition typed outputs and consumers, Persist Artifact multitypes,
  Worker operation parity, and 29 Expression-source defaults are covered; registration + expression tests are
  22 passed and Ruff is GREEN. Docs/workflows/release/full regression/review remain ticket-level work.
- 2026-09-04: Representative workflow generation, README/changelog/release contracts, 23-view visual QA, and final
  integration checks are complete. Visual QA found seven legend/annotation defects; all were fixed and re-rendered.
- 2026-09-04: Two review cycles fixed genuine tascCODA null semantics, table/hierarchy fingerprints, standalone
  code parity, Cell Cycle tie behavior, renderer-side Pseudobulk analysis, and allocation-order/resource bounds.
  Fresh v3 Standards and Spec reviews both report zero findings. Final suite: 1,734 passed, 4 optional skips; Ruff,
  workflow generation, release/registration checks, and `git diff --check` pass.

## Answer

- The registry now contains 144 unique current nodes, including 45 new domain-specific companion Plot nodes. All
  pure Plot nodes expose `plot`, `summary`, and `code`; every adapted node has one matching Worker operation.
- Milo and composition results are portable typed artifacts accepted by Persist Artifact and declared in release
  metadata. Existing example workflows were regenerated against current schemas and now demonstrate core
  clustering/marker and scVI diagnostic plots.
- README and changelog document the categorized plot surface and the rerun requirement for artifacts created before
  the expanded Harmony/scVI/Milo/composition/Augur evidence contracts.
- Full regression, resource/tamper tests, actual PNG inspection, and independent Standards/Spec reviews leave no
  unresolved finding.
