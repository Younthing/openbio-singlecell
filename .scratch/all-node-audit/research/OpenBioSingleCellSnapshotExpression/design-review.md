# Snapshot Expression: module design review

## Current module

The node copies current `X` into an arbitrary layer, `raw`, or both. It is currently a generic storage copier even though the repository defines **Raw snapshot** narrowly as post-QC, full-gene counts saved before highly variable gene selection.

## Correctness and interface problems

- The default `destination="layer"` does not create a Raw snapshot despite the node's role in the workflow.
- Any current `X`, including normalized/log/scaled values, can be stored without count validation.
- Existing `raw` or a same-named layer is overwritten silently, erasing provenance.
- User-selectable destination and layer names let equivalent workflows invent incompatible conventions.
- It does not validate nonempty/finite/non-negative/nonzero input or unique cell/feature identifiers.
- It cannot detect a too-late snapshot after feature subsetting, and it does not disclose when full-gene status cannot be proven.
- Other nodes also write `.raw` or implicitly create count layers, so ownership of the expression-state checkpoint is diffuse.

## Decision: keep, narrow, and enhance

Keep this as the single module that creates the repository-defined Raw snapshot. Do not merge it into normalization or HVG selection: the checkpoint is an independently auditable data-state decision, and downstream branches may consume either canonical counts or the full-feature `raw` state.

This module should become the only production node allowed to create or overwrite `.raw`. Nodes such as log transformation must not expose `set_raw`; package adapters that need a modified raw-like object should create it only inside a private working copy. Count-layer creation in normalization should be migrated to this checkpoint so the canonical convention has one owner.

## Target interface

- Inputs: `adata`; explicit expression `source`, default `X`, restricted to `X` or a layer; advanced `overwrite_existing=False`.
- Fixed hidden policy: create both `.raw` and `layers["counts"]`; `counts` is the canonical, non-configurable layer name.
- Outputs: copied `adata`, structured `summary`, equivalent `code`, in that order.
- Legacy migration: saved `destination="layer_and_raw", layer_name="counts"` maps exactly to the new default. Other legacy destinations should receive a clear migration error rather than silently changing meaning.

The source is visible because it determines which matrix becomes the canonical count state. The canonical key and copy semantics are hidden because varying them adds coupling without changing the scientific question. Overwrite is advanced because it is a provenance/safety operation, not a routine tuning choice.

## Invariants and failure modes

- Require a nonempty matrix with finite, non-negative values and at least one positive count.
- Require unique `obs_names` and `var_names`; do not synthesize identities.
- Detect integer-like values. Accept non-integer corrected counts only with a prominent warning and report limitation.
- Reject existing `.raw` or `layers["counts"]` by default. With explicit overwrite, report exactly which states were replaced.
- Inspect OpenBio analysis history: if an HVG operation with `subset=True` already occurred, fail because full-gene status is known to be lost.
- If provenance cannot establish whether an external object is full-gene/post-QC, proceed only with an explicit limitation; the module cannot infer that fact from shape or values.
- Preserve current `X` exactly, including when the source is another layer. Preserve observation/variable order and all annotations.

## Report contract

`summary` reports input and saved dimensions, total counts, nonzero entries, source, dtype/storage, integer-like status, canonical layer/raw creation, overwrite decisions, provenance evidence, warnings, AnnData reference, and dynamic Python/openbio-singlecell/AnnData/NumPy/Pandas/SciPy versions as used. The core results sentence should be suitable for a methods/results report, for example: “A Raw snapshot containing N post-QC cells and G full-gene features was retained from the declared count source.” When full-gene status is not provable, that sentence must be qualified.

`code` returns a new AnnData, performs the same validation, creates the same canonical layer and `raw`, preserves `X`, and is self-contained apart from public scientific packages. It does not need plugin-only history bookkeeping.

## Cohesion and coupling assessment

The narrowed module has high cohesion: one call establishes one canonical expression checkpoint. Fixing the storage convention behind a small interface gives downstream nodes leverage and makes count provenance local. The old arbitrary destination/layer interface is removed because it spreads convention knowledge across callers.

## Verification plan

- Dense, CSR, and CSC count matrices from `X` and a selected layer; returned `X` remains byte/value equivalent to input.
- Negative, non-finite, all-zero, and duplicate observation/variable identifiers fail; non-integer corrected counts warn and are disclosed.
- Existing `raw`/canonical layer fail unless overwrite is explicit; overwrite reports both replaced states.
- A known post-HVG-subset history fails, while unprovable external full-gene provenance yields a limitation.
- Subsequent observation slicing narrows `raw`, whereas variable slicing retains the Raw snapshot's full feature axis.
- Strict-JSON summary and compiling generated code reproduce `X`, `layers["counts"]`, `raw.X`, and feature identity.
- Repository-level scan proves no other production analysis node creates or overwrites `.raw`.

## Open expert-boundary decision (2026-08-28)

Remove the HVG-history gate and all matrix/axis scientific gates that AnnData storage itself does not require. Convert
non-finite, negative, non-integer, all-zero, empty-axis, and duplicate-axis observations into summary warnings. Keep
source existence/alignment (provided by AnnData layers), in-memory mutation, fixed destination convention, and
explicit overwrite as hard preconditions. Report `history_used=false`, `full_gene_status=user_declared`, and the
selected-source audit. Do not retain history-tamper or Raw-binding tests.
