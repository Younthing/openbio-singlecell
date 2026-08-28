# Load H5AD: module design review

## Current module

The node is a local-file adapter at the H5AD-to-`AnnData` seam. It resolves an uploaded path, calls the official AnnData reader, optionally makes variable names unique, validates non-empty dimensions, and attaches source metadata.

## Correctness and interface problems

- Duplicate or blank observation identifiers are currently accepted, even though downstream joins and cell-level annotations depend on stable cell identity.
- Variable-name repair is performed without disclosing how many names changed or preserving the pre-repair uniqueness state in provenance.
- The interface does not state that the whole object is loaded into memory; large H5AD files can therefore fail for resource reasons that look like data errors.
- It would be scientifically wrong to infer that `X` is counts merely because the source is H5AD. The adapter must remain expression-state agnostic.

## Decision: enhance, keep as an adapter

Keep the node and strengthen validation/provenance. Do not merge it with the 10x loaders: H5AD preserves a complete, potentially processed AnnData object, whereas 10x loaders construct a count-state object from a feature-barcode matrix. Combining them would make one shallow conditional interface with incompatible invariants.

This is an input adapter, not an analysis or transformation node. Its primary and only port remains `adata`; the audit-wide `summary` and `code` ports are not forced here. Reproducibility belongs in the loaded object's source metadata.

## Target interface

- Visible input: uploaded H5AD `path`.
- Advanced input: `make_var_names_unique`, because it is a compatibility/naming policy rather than a biological question.
- Hidden fixed behavior: in-memory `backed=None`; preservation of all H5AD slots; sandboxed path resolution; stable file fingerprinting; no expression/count interpretation.
- Output: loaded `adata` with bounded, JSON-safe OpenBio source metadata.

## Invariants and errors

Require an existing `.h5ad` under the input directory, successful native-format decoding, non-empty axes, and non-empty unique observation identifiers. If variable identifiers are duplicated, either apply the explicit suffix policy or preserve and warn according to the user's setting; never silently change observation identifiers. Surface memory/read errors with the file context and do not return a partially loaded object.

## Depth and test surface

The module hides path safety, fingerprinting, AnnData I/O, identity checks, optional variable-index repair, and provenance behind two inputs. AnnData/HDF5 are local dependencies; no additional adapter seam is justified. Tests should cross the node interface and assert slot preservation, identifier policy, provenance, and failure modes.

## Open expert-boundary decision (2026-08-28)

Revise the earlier invariants: file/path safety and native H5AD decoding remain hard errors, while empty axes, blank
axis labels, and duplicate observation/variable labels are load-time advisories. Their consequences belong in source
metadata so an expert can inspect the object or repair it explicitly. Do not use embedded analysis history as a gate
and do not add history-tamper tests. This keeps the adapter deep at the file-decoding seam without turning it into a
workflow-policy validator.
