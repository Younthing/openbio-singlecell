# Export CSV: module design review

## Current module

The node accepts a typed `TableResult`, selects a contained output path, and calls
`DataFrame.to_csv(..., index=False)` directly on that final path.

## Correctness and scientific-practice findings

- `index=False` silently discards the row axis. Many scientific result tables encode genes,
  cells, samples, contrasts, or compound IDs in that axis, so the current export can lose the
  only row identifier.
- Direct writing can truncate an existing file or expose a partial CSV after failure.
- Unnamed or colliding index headers can make later import ambiguous.
- CSV cannot preserve pandas dtype/category metadata; the adapter must not imply a lossless
  scientific object round trip.

## Decision: enhance and keep separate

Keep CSV export as a narrow terminal adapter. It should not merge with table-producing
analysis nodes: the same typed result may be previewed, filtered, or exported independently,
and computation must remain free of filesystem side effects.

This node is not an analysis/processing node. It intentionally has no `summary` or `code`
output; it persists an upstream disclosure artifact without creating a new scientific result.

## Target interface

- Visible inputs: typed `table`, `filename_prefix`.
- Advanced input: `overwrite`.
- Graph outputs: none; UI file reference only.
- Hidden fixed behavior: index preservation, deterministic index-label policy, UTF-8, LF,
  stable missing-value spelling, output containment, secure staging, fsync, atomic commit,
  non-overwrite collision protection, and cleanup.

## Invariants and tests

Require `TableResult`, a single-level uniquely named column axis, and a representable one- or
multi-level row index. Reject duplicate emitted header labels, MultiIndex columns, and
NUL-containing column/index names. Tests should cover named, unnamed, and multi-level row
indices; Unicode; missing values; header collisions; row and column order; failure atomicity;
and concurrent non-overwrite collision behavior.
