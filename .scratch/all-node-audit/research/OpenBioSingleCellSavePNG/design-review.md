# Save PNG: module design review

## Current module

The node accepts `PlotResult`, selects a contained `.png` path, and writes the byte payload
directly to the destination. `PlotResult` construction checks only the PNG signature.

## Correctness and cohesion findings

- Signature bytes do not prove that the PNG is complete or decodable.
- Direct destination writing can replace a valid image with a corrupt or partial file.
- Naming/path policy, byte validation, and durable commit belong inside this one adapter; plot
  creation and scientific interpretation do not.

## Decision: enhance and keep separate

Keep this terminal adapter distinct from Preview Result. Durable output uses a user-selected
name and collision policy, while previews use deterministic temporary replacement. Keep it
distinct from plotting analyses so those nodes remain pure and reusable.

This is not an analysis/processing node. It has no graph outputs and does not add `summary`
or `code`; it writes an already disclosed plot artifact.

## Target interface

- Visible inputs: typed `plot`, `filename_prefix`.
- Advanced input: `overwrite`.
- Graph outputs: none; UI image reference only.
- Hidden fixed behavior: full PNG verification, byte preservation, output containment,
  same-directory secure staging, flush/fsync, atomic commit, race-safe non-overwrite, cleanup.

## Invariants and tests

Reject non-`PlotResult` values and signature-only, truncated, non-PNG, or oversized malformed
payloads before publication. Verify bytes before touching the destination. Tests should use
real PNGs and cover exact-byte persistence, failure cleanup, prior-file preservation,
incrementing names, overwrite, and a late non-overwrite collision.
