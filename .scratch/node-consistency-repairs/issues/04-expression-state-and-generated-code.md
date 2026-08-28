# Centralize expression state and clean generated preprocessing code

Type: task
Status: resolved
Blocked by: 03

## Question

Which expression-state rules can live behind one runtime Interface, and which minimal standalone helpers must remain in
generated code so that equivalent functions are portable without mutating OpenBio history?

## Evidence

- `.scratch/all-node-audit/research/OpenBioSingleCellCNMF/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellCNMF/design-review.md`
- Preprocessing node records under `.scratch/all-node-audit/research/`

## Acceptance criteria

- cNMF recognizes `scale_to_layer` as scaled state consistently with the canonical resolver.
- cNMF runtime and standalone code agree with the authoritative resolver for `scale_to_layer`.
- Generated preprocessing functions reproduce primary scientific outputs but do not add or invent `analysis_history`.
- Generated-code and runtime validation remain behaviorally equivalent for scientific inputs and results.

## Comments

- 2026-08-29: Claimed after the named-graph regression suite passed. This ticket is split into two tracer bullets:
  cNMF `scale_to_layer` classification, then removal of generated preprocessing history writes.
- 2026-08-29: cNMF runtime/generated parity passes in the full factorization suite (`38 passed`). Basic preprocessing
  generated code no longer contains a history writer and its suite passes (`34 passed`).
- 2026-08-29: Deleted the zero-caller `read_10x_mtx` pass-through; the internal audited loader remains the sole Module
  used by both public input nodes.

## Answer

Fix the concrete classifier drift in the standalone cNMF implementation, and keep generated preprocessing functions
scientifically equivalent but provenance-read-only. Do not make emitted source impersonate a plugin runtime by
inventing metadata defaults or signing history. Broader resolver consolidation is split into ticket 06 because the
current consumers use two scientifically distinct logged-state vocabularies that require one explicit design decision.
