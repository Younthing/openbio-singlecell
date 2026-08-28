# Filter Genes: module design review

## Current module

The module performs one atomic low-information feature filter using an intersection of count and detection bounds. Its copy-on-write behavior is correct.

## Problems

- Count semantics are implicit in `X`.
- Disabled criteria and contradictory ranges are not explicitly handled in the report/interface.
- The key effect size—how many genes were removed—is hidden in the resulting shape and history.
- No reproducible source is exposed.

## Decision: enhance, do not merge

Keep separate from cell filtering because cells and genes are different scientific objects with different reporting implications. Keep separate from highly variable gene selection because the Raw snapshot/full-gene state must remain available for marker evidence and condition contrasts.

## Target interface

- Inputs: `adata`, count expression source, visible minimum thresholds, visible maximum thresholds.
- Outputs: filtered `adata`, structured `summary`, equivalent Python `code`.
- Invariants: selected source exists; active min/max ranges are consistent; source statistics can be aligned without ambiguous duplicate feature identity.
- Warnings: no active filter; all genes retained; all genes removed; finite signed/fractional source is not verified counts; Raw/current feature coverage; rare-feature loss should be considered for the study.

## Depth and dependencies

Sparse/dense reduction, mask composition, result summarization, citations, versions, and code generation are hidden in one in-process implementation. No adapter seam is warranted.

The node does not inspect mutable history or require a Raw/current binding. For an explicitly selected broader Raw axis, it computes in Raw and maps statistics onto current features by exact unique identifier. Current features absent from that selected source are retained and counted as unevaluated. A zero-gene output remains a valid explicit filtering outcome and is reported as likely unusable downstream rather than prohibited.

## Freeze-review interface correction

Keep the existing numeric ids and Python positional prefix, then add explicit `enable_min_counts` and `enable_max_counts` controls. Enabled total-expression bounds are finite floats and may be negative, fractional, or exactly zero. Inactive values are ignored. Direct legacy calls with omitted enable arguments infer the old non-zero rule, while serialized workflows require central migration to materialize the two booleans. Generated code consumes the resolved controls, not the historical zero sentinel.
