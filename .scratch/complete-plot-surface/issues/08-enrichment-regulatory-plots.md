# Add enrichment and regulatory plots

Type: task
Status: resolved
Blocked by: 01

## Scope

- Stored Score/Activity, pathway contrast, GSEA, ORA, drug enrichment, TF ranking, SCENIC RSS, binarization, and membership Plot adapters.

## Acceptance criteria

- Each adapter validates its producer/schema and reports the exact score, effect, enrichment, overlap, or threshold encoded by its visual marks.

## Answer

- Added enrichment companions for stored Score/Activity, Pathway Score contrast, Ranked GSEA, gene-set ORA, Drug Hypergeometric, and Drug GSEA evidence.
- Added regulatory companions for typed CollecTRI activity, TF activity ranking, typed SCENIC activity, regulon specificity, binarization, and final regulon membership.
- Every public adapter stays in its scientific category, uses a closed view union, returns `plot`, strict `summary`, and equivalent standalone `code`, and has a one-shot Worker operation.
- Canonical result tables now carry a portable current-content plot-evidence fingerprint; activity adapters validate their existing typed artifact fingerprints and exact observation axes. The existing score-artifact validator now normalizes H5AD-restored `score_names` on a private validation copy.

## Comments

- 2026-09-04: Claimed; implementation starts after the shared convention is frozen.
- 2026-09-04: Completed 12 domain-specific Plot adapters and all requested views. Focused domain regression: 195 passed. New adapter suite: 37 passed. Ruff and `git diff --check` passed for the owned surface.
