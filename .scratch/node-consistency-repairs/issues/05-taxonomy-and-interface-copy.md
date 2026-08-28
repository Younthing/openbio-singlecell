# Align node taxonomy and interface copy

Type: task
Status: resolved
Blocked by: 04

## Question

Which node categories and concise descriptions best expose the already-established domain distinctions without merging
atomic scientific analyses?

## Acceptance criteria

- Leiden and its resolution sweep share a clustering taxonomy.
- Cluster marker evidence is not presented as a formal Condition contrast, consistent with ADR-0001.
- Touched schemas describe their scientific role and important source/assumption parameters without exposing internal
  implementation knobs.
- Registration, workflow, and release manifests remain exact.

## Comments

- 2026-08-29: Claimed after the expression/generated-code repairs passed. Category changes are presentation-only and
  intentionally do not merge nodes or alter schemas, execution parameters, outputs, or scientific implementations.
- 2026-08-29: Final verification passed: the full Python suite (`1164 passed, 4 skipped`), 14 frontend tests, five
  generated-workflow checks, Ruff, Python byte-compilation, 10 JavaScript syntax checks, and `git diff --check`.
  Registration, packaged-workflow, and release-manifest assertions are included in the full suite.

## Answer

Classify each node by its atomic scientific claim rather than by source file or whether it renders a plot. Leiden is
`clustering`; MarkerGenes and FilterMarkerGenes are `marker-evidence`; AnnDataSummary and Population Centroid
Correlation are `diagnostics`; CellCycleScore is `annotation`. Retain all six as cohesive nodes: none becomes simpler
or more correct by merging it with data loading, representation construction, Condition inference, or trajectory
analysis.
