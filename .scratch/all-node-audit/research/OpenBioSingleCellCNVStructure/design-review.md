# CNV Structure: module design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. Its atomic CNV analysis/view
> nodes remain, with no compatibility ID or automatic workflow migration.

## Decision

**Retire the executing `OpenBioSingleCellCNVStructure` composite and keep its ID only as a nonexecuting migration shim.** It fails the atomicity/cohesion requirement by hiding five independently meaningful analyses. A `summary`/`code` pair on the composite would document a recipe, not repair the interface.

## Replacement graph

Migration expands the old node into:

1. a new, separately researched `OpenBioSingleCellCNVPCA` consuming `OPENBIO_CNV_STATE` and producing a portable AnnData with `X_cnv_pca`, `summary`, and `code`;
2. existing audited `Neighbors` with `use_rep="X_cnv_pca"`, `key_added="cnv_neighbors"`;
3. existing audited `Leiden` with `neighbors_key="cnv_neighbors"`, `key_added="cnv_leiden"`;
4. existing audited `UMAP` with the same graph and `key_added="X_cnv_umap"`;
5. a new, separately researched `OpenBioSingleCellCNVScore` taking both the original typed CNV state and the
   downstream AnnData carrying an explicit grouping, verifying their exact shared CNV fingerprints, and returning
   `adata`, a group score table, `summary`, and `code`.

The two new IDs must receive their own exactly two pre-change documents before implementation. No replacement should infer a tumor label or cutoff.

## Shim, migration and tests

The nonexecuting shim accepts the exact serialized legacy input (`adata`) only so workflows can load; the old node exposed no scientific parameters. It has no executable analysis outputs and always raises an actionable message identifying the five replacements. Packaged workflows are structurally migrated; no run may silently execute old hidden defaults. Because the shim is no longer an analysis node, it has no fake `summary` or `code` outputs. Each replacement has precise typed state/graph inputs and its own result plus required `summary`/`code` endpoints.

Migration must preserve the infercnvpy historical defaults explicitly (PCA solver/centering, neighbor count/metric, Leiden resolution/seed, UMAP parameters/seed, score grouping) or, when an old backend default was never serialized, record that exact reproduction is impossible and require review. Tests assert the shim never executes, migration produces ordered wired nodes, all parameters are explicit, no legacy ID remains in packaged workflows, and replacement outputs independently satisfy report/code equivalence.

The implemented migration performs one whole-workflow preflight and then one commit. A Structure is expandable only
when its exact legacy AnnData input directly traces, in the same graph scope, to an eligible legacy InferCNV node;
all of that InferCNV output's consumers must be paired Structures. It inserts CNV PCA, reviewed Neighbors
(`X_cnv_pca`, all available dimensions, 15, Euclidean, UMAP connectivities, `cnv_neighbors`, seed 0), reviewed
Leiden (resolution 1, `cnv_leiden`, `cnv_neighbors`, convergence iterations, five audited starts, seed 0), and
reviewed UMAP (`cnv_neighbors`, min_dist 0.5, spread 1, `X_cnv_umap`, seed 0), then converts Structure to CNV
Score. The migration note explicitly records ARPACK/uncentered PCA and that the unrecorded component count plus
version-dependent historical Leiden flavor cannot be reproduced exactly; `n_comps=30` and the audited igraph
five-start implementation require review.

Object and array links, positional/named widgets, direct/proxied retained controls, nested single-host subgraphs,
global node/link allocation, typed-state fan-out, and downstream consumer links are preserved. Subgraph boundary
links, reroutes/link extensions, orphan/multi-host/duplicate definitions, disabled/bypassed nodes, connected unknown
controls, partial/mixed schemas, broken links, missing hints, empty references, and non-paired producers reject the
whole graph byte-for-byte. A second run returns zero; exact current InferCNV/CNV PCA/CNV Score/shim schemas are
no-ops.

## Registry disposition

Keep the old ID deserializable with `is_deprecated=True`, a visibly retired display name, and an actionable
non-executing description. Registration tests must verify these public signals and the existing fail-before-read
contract. This does not collapse the five atomic replacements back into a composite and does not add fabricated
reporting ports to the shim.
