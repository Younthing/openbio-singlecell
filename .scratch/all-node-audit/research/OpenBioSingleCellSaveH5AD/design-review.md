# Save H5AD: module design review

## Current module

The node resolves a safe path and calls `adata.write_h5ad()` directly on the final filename.
It is a terminal output node with no graph outputs.

## Correctness and interface findings

- Direct final-path writing can destroy a valid prior file or leave a partial HDF5 container
  when serialization fails.
- The current code reports success without reopening the file, so a malformed or axis-mismatched
  artifact could be published unnoticed.
- The counter lookup and later write form a non-overwrite race: two executions can select the
  same name.
- Compression is a meaningful storage/performance choice, but exotic HDF5 filters should stay
  hidden and unsupported because they weaken portability.

## Decision: enhance and keep as a terminal adapter

Keep the node. H5AD persistence is one atomic I/O responsibility and should not be merged
with preprocessing or loaders. The node must preserve the supplied AnnData rather than infer
or alter its expression state.

This is not a scientific analysis/processing node. It has no primary data output and does not
add `summary` or `code`; analysis disclosures are already embedded in OpenBio metadata and
upstream report artifacts.

## Target interface

- Visible inputs: `adata`, `filename_prefix`.
- Advanced inputs: `overwrite`; `compression` limited to `gzip`, `lzf`, or `none`.
- Graph outputs: none; UI file reference only.
- Hidden fixed behavior: output-root containment, extension normalization, same-directory
  secure staging, HDF5 close, file fsync, backed read verification, atomic commit, race-safe
  non-overwrite, and cleanup.

## Invariants and tests

Require a real AnnData object and a contained destination. The round-trip verifier compares
shape, ordered `obs_names`, and ordered `var_names`; it does not reinterpret biology or load
the full matrix. Tests should inject writer/verifier/commit failures and prove both absence of
partial new files and byte-for-byte preservation of an old file. A real round trip must retain
OpenBio history and supported AnnData slots.
