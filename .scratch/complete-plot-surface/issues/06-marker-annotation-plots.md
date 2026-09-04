# Add marker and annotation plots

Type: task
Status: resolved
Blocked by: 01

## Scope

- Marker Evidence auto-panel, Marker ORA Evidence, CellTypist Diagnostics, and Cell Cycle Score plots.

## Acceptance criteria

- Cluster marker evidence never becomes a Condition contrast or Curated annotation.
- Expression-derived views use one explicit Expression source.

## Answer

- `Marker Evidence Plot` consumes the exact AnnData, marker table, and bound tested-gene universe; validates
  producer/content/axis provenance; selects the first N stored ranks per group without reranking; and renders
  mean explicit-source expression plus detected-cell fraction.
- `CellTypist Diagnostics Plot` renders closed confidence-distribution and grouped probability-heatmap views
  from stored Provisional annotation. CellTypist now retains observation-axis, probability-content, and selected-
  annotation fingerprints so the plot rejects stale or changed evidence.
- `Cell Cycle Score Plot` renders the stored S-versus-G2M score plane and phase counts. Cell Cycle Score now
  records observation-axis and score-bundle fingerprints in producer history; the plot validates both plus the
  canonical Scanpy phase rule and never reads expression or reruns scoring.
- `Marker ORA Evidence Plot` keeps this separate annotation-domain contract and renders closed enrichment-dot
  and overlap-bar views from stored contingencies, overlaps, log odds, ranks, and adjusted p-values. The producer
  now fingerprints the complete canonical evidence table; the plot never opens the source set file or retests.
- All four nodes return `plot`, strict `summary`, and standalone equivalent `code` through one-shot Worker
  operations, preserve input artifacts, bound static rendering, and keep Cluster marker evidence, Provisional
  annotation, Curated annotation, and Condition inference distinct.

## Comments

- 2026-09-04: Claimed; implementation starts after the shared convention is frozen.
- 2026-09-04: RED was recorded as missing plot interfaces/fingerprint helpers. Focused schema, owned-operation,
  tamper, immutability, code-parity, Worker artifact, and Matplotlib-state tests are GREEN (36 plot tests).
  Relevant producer regressions and Ruff are GREEN.
- 2026-09-04: The combined results/annotation/cell-cycle regression run is 184 passed with two expected central
  integration failures: the expression-source registry still expects 27 sources and has not yet added Marker
  Evidence Plot to the default-layer set. Central expected-count maintenance is intentionally ticket 11 scope.
- 2026-09-04: Completeness audit added the separate Marker ORA plot after confirming the generic enrichment ORA
  renderer intentionally rejects `marker_ora_evidence`; annotation regression is 74 passed and Ruff is GREEN.
- 2026-09-04: Representative Marker Evidence visual QA found and then verified the fix for its missing detected-
  fraction scale: a 25%/50%/100% size legend now sits outside the data axes, without occlusion or clipping.
