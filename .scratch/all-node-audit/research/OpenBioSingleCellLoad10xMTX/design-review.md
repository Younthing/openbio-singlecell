# Load 10x MTX: module design review

## Current module

The node is a format adapter that resolves a 10x MEX directory, reads Matrix Market and TSV files, constructs a sparse `AnnData`, and attaches provenance. A shared private reader is also used by the Study loader.

## Correctness and interface problems

- The current `make_unique` path calls both `obs_names_make_unique()` and `var_names_make_unique()`. Official Scanpy behavior only repairs the variable index; changing duplicate barcodes hides an invalid cell-identity input.
- Matrix values are not currently checked for finite, non-negative, integer-like counts.
- Blank/duplicate barcodes and invalid feature IDs are not rejected explicitly.
- The behavior for missing feature-type metadata is implicit, so `gex_only=True` can look verified when it could not be applied.

## Decision: enhance, do not merge

Keep this adapter. Do not merge it with `Load 10x H5`: the two formats need different file resolution, integrity checks, and format-specific parameters. Share private count/index validation to create locality without broadening either public interface.

This node performs ingestion, not scientific analysis. It therefore keeps a single `adata` output and does not receive forced `summary`/`code` ports. Its source metadata is the portable ingestion provenance.

## Target interface

- Visible: `directory`, `var_names` (`gene_symbols` or `gene_ids`), and `gex_only`; these choices determine retained features and the variable index.
- Advanced: `make_unique`, because it controls only a technical variable-index repair.
- Hidden fixed behavior: auto-detection of supported compressed/uncompressed 10x filenames, sparse CSR output, cell-by-feature orientation, barcode rejection rather than repair, original ID/symbol preservation, count validation, path sandboxing, and file fingerprinting.
- Output: count-state `adata` with source metadata.

## Invariants and errors

Require one unambiguous file triplet, matching dimensions, at least one retained cell and feature, valid integer-like counts, unique non-empty barcodes, and valid feature identifiers. If `gex_only` removes every feature, fail explicitly. If gene symbols collide and uniqueness repair is disabled, preserve them only with a clear provenance warning; gene IDs remain available for unambiguous mapping.

## Depth and test surface

The module hides format discovery, decompression, Matrix Market orientation, legacy/v3 feature parsing, modality selection, index construction, validation, and provenance behind four inputs. Dependencies are local/in-process. Tests should use small real-format fixtures and assert orientation, values, identifiers, filtering, compressed/uncompressed parity, source metadata, and rejection of malformed inputs.

## Open expert-boundary decision (2026-08-28)

Revise the earlier identity gate: uniqueness and non-emptiness are downstream usability properties, not prerequisites
for parsing a dimensionally valid MEX payload. Preserve questionable identifiers and disclose counts rather than
blocking expert ingestion. Keep hard errors only for path/file-role ambiguity, malformed tables, dimension mismatch,
and values incompatible with the declared 10x count format. All-zero or zero-dimensional payloads are advisories.
