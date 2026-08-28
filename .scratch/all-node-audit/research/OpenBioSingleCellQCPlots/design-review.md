# QC Plots: module design review

## Current module

This is a read-only reporting module with one coherent purpose: visualize cell-level QC evidence. It already returns a concrete plot artifact and does not mutate `AnnData`.

## Problems

- When mitochondrial annotations are unavailable, it plots zeros. A warning exists, but the image can still be misread as observed zero mitochondrial content.
- Counts and detected genes are recalculated from implicit `X` even when official QC columns may already represent the intended count source.
- It returns no structured distribution summary or code.

## Decision: enhance, do not merge

Keep separate from metric calculation. A plot is a report adapter over available QC evidence and may be rerun without modifying data. Merging would make a read-only visualization unexpectedly transform `AnnData`.

## Target interface

- Inputs: `adata`; expression source used only when standard QC columns are absent.
- Outputs: `plot`, structured `summary`, equivalent Python `code`.
- Invariants: at least one cell and one gene; finite plottable values; selected fallback source exists.
- Warnings: missing mitochondrial evidence, fallback-derived metrics, non-finite observations dropped from display.
- Provenance rule: an existing metric is used as a whole column. When mitochondrial percentage must be derived, its numerator, denominator, and scatter-plot total all come from the same selected expression source; historical `obs` totals are never mixed into that calculation.

## Depth and dependencies

The module hides metric selection/fallback, robust distribution summarization, figure construction, PNG conversion, report construction, and code rendering. Dependencies are in-process; no adapter seam is needed.

## Open expert boundary decision

Raw/current equality and analysis-history provenance are not gates. The selected representation owns fallback totals, detected features, and any derived mitochondrial numerator/denominator. If the selected representation lacks compatible mitochondrial annotations, that panel is omitted and disclosed. Signed or fractional finite evidence is drawn with warnings; the only hard data-domain boundary is that each required primary panel needs at least one finite value.
