# Filter Cells: module design review

## Current module

The module atomically keeps cells satisfying the intersection of up to five criteria. This is more useful than Scanpy's one-threshold-per-call interface and keeps the complete filtering decision in one provenance entry.

## Problems

- It always derives counts and detected genes from `X`, which may no longer contain counts.
- Zero disables a bound, but a fully disabled filter produces no warning.
- It exposes only `adata`; callers cannot directly report retained/removed cells or reconstruct the exact operation.
- It does not reject contradictory min/max pairs early.

## Decision: enhance, do not merge

Keep one atomic cell-filtering node. Do not merge with metric calculation: measurement should precede and inform a decision. Do not split one node per threshold: that would scatter one scientific filtering decision across several shallow modules and provenance entries.

## Target interface

- Inputs: `adata`, count expression source, visible min/max count and detected-gene bounds, visible mitochondrial limit, advanced mitochondrial column.
- Outputs: filtered `adata`, structured `summary`, equivalent Python `code`.
- Invariants: non-empty input matrix; selected source exists; active minimum does not exceed active maximum; mitochondrial column exists and is decidable when used.
- Warnings: no active filter; all cells retained; all cells removed; finite signed/fractional source is not verified counts; thresholds should be reviewed per Sample/tissue.

## Depth and dependencies

The in-process implementation hides sparse/dense reductions, combined masking, result metrics, reporting, and code rendering. The interface remains one call with one coherent filtering decision; no external adapter seam is needed.

An empty output is a valid explicit selection result and is not a runtime error. The report carries `retained_cells=0`, the removal rate, and a downstream-usability warning. The node never consults analysis history or binds a selected Raw source to current `X`.

## Freeze-review interface correction

- Keep the legacy numeric input ids and Python positional prefix, but expose `enable_min_counts`, `enable_max_counts`, and `enable_max_pct_mito` as result-defining controls.
- Enabled total-expression bounds accept finite floats without a sign restriction. Disabled values are ignored rather than validated. Enabled min/max bounds must still be ordered because otherwise the declared interval is contradictory.
- A missing enable argument in direct legacy Python calls infers the former non-zero activation rule; the public schema supplies explicit booleans. Saved workflows must be migrated to explicit booleans before execution.
- Runtime and generated code must both count sparse non-zero values without canonicalizing the returned AnnData as an incidental side effect.
