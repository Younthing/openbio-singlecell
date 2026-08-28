# Load 10x H5: module design review

## Current module

The node is a 10x HDF5 adapter. It resolves a local file, delegates parsing to the official Scanpy reader, optionally makes both axis names unique, checks non-empty dimensions, and attaches source metadata.

## Correctness and interface problems

- The current uniqueness option calls `obs_names_make_unique()`, silently suffixing duplicate barcodes and hiding ambiguous cell identity.
- The returned matrix is assumed to be counts but is not checked for finite, non-negative, integer-like values.
- Naming repair, selected genome, and feature filtering are not fully disclosed in provenance.
- Input validation does not distinguish a valid generic HDF5 file from the required 10x feature-barcode schema until execution; execution errors should be translated with file/format context.

## Decision: enhance, do not merge

Keep the dedicated H5 adapter and share validation helpers with the MTX adapter. A universal 10x loader would need format conditionals and mutually irrelevant options (`directory` versus file, compression/layout versus legacy genome), reducing interface depth.

This is ingestion rather than analysis, so the only output remains `adata`. No forced `summary` or `code` ports are added; full ingestion details belong in source metadata.

## Target interface

- Visible: uploaded `path` and `gex_only`, because modality retention changes the output feature set.
- Advanced: `genome` for legacy multi-genome files and `make_unique` for displayed variable-name repair.
- Hidden fixed behavior: official Scanpy reader, `backup_url=None`, in-memory output, barcode rejection rather than repair, count validation, metadata preservation, sandboxed path resolution, and fingerprinting.
- Output: count-state `adata` with JSON-safe source provenance.

## Invariants and errors

Require a supported extension, valid 10x H5 schema, at least one cell and retained feature, unique non-empty barcodes, valid counts, and consistent feature annotations. If a legacy multi-genome file needs `genome`, preserve Scanpy's requirement but report it as an actionable domain error. If variable names collide, repair only when requested and disclose the number changed.

## Depth and test surface

The module hides secure file resolution, 10x schema parsing, legacy genome handling, modality selection, identity/count checks, variable-index policy, and provenance behind four inputs. Tests should cover modern and legacy fixtures, modality filtering, duplicate barcode rejection, variable-name repair, invalid counts/schema, and source metadata.

## Open expert-boundary decision (2026-08-28)

Relax empty-axis, all-zero, blank-label, and duplicate-label checks to warnings recorded in source provenance. Retain
only reader/file errors and the declared 10x count-format checks as hard failures. This makes H5 and MEX adapters
consistent without merging their distinct decoding implementations.
