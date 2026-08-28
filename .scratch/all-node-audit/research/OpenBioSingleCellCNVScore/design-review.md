# CNV Score: module design review

## Decision

**Add `OpenBioSingleCellCNVScore` as the atomic descriptive-score replacement for the final step of the retired CNV
Structure composite.** It consumes one validated inferred CNV artifact, one downstream AnnData carrying the
partition, and one explicit categorical partition key. The separate AnnData input is required because PCA,
Neighbors and Leiden extend an AnnData copy rather than mutating the immutable inferred-state artifact. Exact axes
and all CNV/source/reference/coordinate fingerprints must match the artifact before the annotation is used. It does
not require or imply CNV PCA, neighbors, Leiden or UMAP; any biologically meaningful categorical grouping is allowed
and disclosed.

Inputs are `cnv_state`, `adata`; visible `groupby`; advanced `output_key="cnv_score"` and
`overwrite_existing=False`. Outputs are the annotated `adata`, a group-level `table`, strict JSON `summary`, and
equivalent `code`. Fix infercnvpy `use_rep` from the artifact, `inplace=True`, and `obs_key=None`.

Validate artifact identity/tamper evidence, exact downstream AnnData state/axis agreement, and a complete categorical partition with at least two
represented groups, exact typed labels and output collision before execution. Remove unused categories only on the
private output copy and disclose them. Call the public backend, then independently recompute
`mean(abs(X_cnv[group_cells, :]))` for every complete group and require exact/tight agreement with every emitted
cell value. The group table contains typed group identity/display, cell count, score, and finite absolute-CNV
summaries, sorted by score without inventing a cutoff. Bind output/provenance fingerprints to the input artifact and
partition.

`summary` contains report-ready Methods/Results, every group score, cell/support accounting, source/reference/window
provenance, warnings, citations, versions, and strong descriptive limitations. Standalone `code` accepts
`(cnv_state, adata)` and performs the same artifact/bundle validation and official call without importing OpenBio
helpers. Tests cover mismatched downstream state, typed/missing/colliding labels, unused levels, dense/sparse
matrices, tampering/collisions, independent formula, malformed backend output,
input immutability, real 0.6.1 smoke, strict JSON, compilation and parity.

During the atomic legacy expansion the old Structure node object is converted in place to CNV Score, so every
existing slot-zero AnnData consumer and link ID now originates from the Score `adata` output without a transient or
duplicate legacy node. One new typed branch connects the original migrated InferCNV state to `cnv_state`; the final
UMAP AnnData connects separately to `adata`. Migration fixes `groupby="cnv_leiden"`,
`output_key="cnv_score"`, and `overwrite_existing=False`, and appends `table + summary + code`. Multiple legacy
Structure nodes may reuse one unambiguous typed state through distinct collision-free links. Boundary, reroute,
shared-definition, malformed-link, or non-Infer producers fail before any graph mutation. Exact current CNV Score
schemas are no-ops.
