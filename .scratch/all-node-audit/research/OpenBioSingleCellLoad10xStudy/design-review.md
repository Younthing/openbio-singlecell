# Load 10x Study: module design review

## Current module

The node discovers first-level directories, loads each with the shared 10x MTX adapter, and concatenates a mapping with `anndata.concat(..., label=sample_key, index_unique="-")`.

## Critical correctness problem

`_discover_10x_study` currently catches errors from a child directory and silently continues. A misspelled, incomplete, corrupt, or unsupported Sample directory therefore disappears from the returned Study without an error. This is a P0 scientific-integrity bug because Sample is the biological replicate and the omission can alter downstream Condition contrasts.

Additional gaps are the absence of an explicit empty-intersection error, incomplete disclosure of inner/outer feature alignment, and inherited MTX barcode/count validation weaknesses.

## Decision: enhance; keep distinct from the single-Sample loader

Keep the node as a Study-assembly adapter and use the single-Sample MTX reader as a private implementation seam. Do not merge the public nodes: loading one Sample and assembling a multi-Sample Study have different invariants, error modes, and scientific meaning. Deleting the Study adapter would spread directory inventory and concatenation policy across workflows.

As an input adapter, it exposes only `adata`; it does not receive forced `summary`/`code` ports. The exact Sample inventory and join decision must be captured in source metadata.

## Target interface

- Visible: Study `directory`, `sample_key`, `var_names`, `join`, and `gex_only`. Sample identity and feature-set alignment are result-defining.
- Advanced: `make_unique`, limited to the variable index in each Sample.
- Hidden fixed behavior: strict first-level inventory; deterministic case-insensitive directory order; shared MTX validation; `axis="obs"`; `merge="same"`; Sample-key mapping; global index construction with `"-"`; path safety and per-file fingerprints.
- Output: one count-state `adata` with an explicit Sample column and complete source provenance.

## Invariants and errors

Require at least one child directory and require every child to be a valid Sample. Aggregate all invalid-child diagnostics in one error so users can repair the whole inventory. Require distinct non-empty Sample names, a non-empty `sample_key`, valid counts and identifiers in each Sample, and at least one retained feature after the selected join. Never downgrade an invalid Sample to a warning.

For `outer`, record the per-Sample missing-feature counts and warn about sparse zero fill. For `inner`, record how many features were dropped from each Sample. The adapter does not claim that directories are biological replicates beyond the user's explicit directory organization; it makes that assertion visible and auditable.

## Depth and test surface

The module earns depth by hiding strict discovery, N Sample loads, feature alignment, index construction, and aggregate provenance behind a compact interface. Dependencies are local/in-process. Tests must prove deterministic inventory, no silent omission, multiple simultaneous diagnostics, Sample labels/indexes, inner/outer semantics, and complete provenance.

## Open expert-boundary decision (2026-08-28)

Remove the special hard gate for an empty inner feature intersection and return the valid zero-variable AnnData with
a prominent provenance warning. Keep strict Sample-inventory decoding and feature-identity conflict checks: those are
required to know what is being concatenated and to avoid conflating different stable features. Do not use analysis
history or Raw equality as evidence at this adapter seam.
