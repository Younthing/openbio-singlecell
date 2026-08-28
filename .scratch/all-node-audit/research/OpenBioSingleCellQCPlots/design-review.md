# QC Plots: module design review

## Current module

This is a read-only reporting module with one coherent purpose: visualize cell-level QC evidence. It already returns a concrete plot artifact and does not mutate `AnnData`.

## Problems found during the audit

- When mitochondrial annotations are unavailable, it plots zeros. A warning exists, but the image can still be misread as observed zero mitochondrial content.
- Counts and detected genes were recalculated from implicit `X` when official QC columns were absent.
- It returns no structured distribution summary or code.

## Decision: enhance, do not merge

Keep separate from metric calculation. A plot is a report adapter over available QC evidence and may be rerun without modifying data. Merging would make a read-only visualization unexpectedly transform `AnnData`.

## Target interface

- Inputs: `adata`; one explicit expression source used for the complete plotted metric family.
- Outputs: `plot`, structured `summary`, equivalent Python `code`.
- Invariants: at least one cell and one gene; finite selected-source values; selected source exists.
- Warnings: missing mitochondrial evidence, non-count-like selected sources, non-finite derived observations dropped from display.
- Provenance rule: total expression, detected genes, mitochondrial numerator, denominator, and scatter-plot totals are
  derived atomically from the selected expression source. Pre-existing `obs` QC columns are not treated as a cache
  because they carry no immutable matrix binding.

## Depth and dependencies

The module hides metric selection/fallback, robust distribution summarization, figure construction, PNG conversion, report construction, and code rendering. Dependencies are in-process; no adapter seam is needed.

## Open expert boundary decision

Raw/current equality and analysis-history provenance are not gates. The selected representation owns totals, detected features, and any derived mitochondrial numerator/denominator. If the selected representation lacks compatible mitochondrial annotations, that panel is omitted and disclosed. Signed or fractional finite evidence is drawn with warnings; the only hard data-domain boundary is that each required primary panel needs at least one finite value.

## 2026-08-29 consistency repair

The post-refactor audit found that the implementation still preferred each existing `obs` metric independently, so a
single plot could combine different expression states. The node now ignores those unbound cache columns and derives the
complete metric family once from the explicit source. Runtime and generated code call the same source-embeddable pure
metric function, use zero-library `NaN` mitochondrial percentages matching Scanpy semantics, and preserve the input
AnnData.
