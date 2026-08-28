# CNV PCA: module design review

## Decision

**Add `OpenBioSingleCellCNVPCA` as the atomic dimensionality-reduction replacement for the first step of the retired
CNV Structure composite.** It consumes a typed inferred CNV state and produces one portable AnnData copy containing
the named CNV PCA representation, plus `summary` and equivalent `code`. Graph construction, clustering, UMAP and
scoring remain separate nodes.

Inputs are `cnv_state`; visible `n_comps=30`; advanced `output_key="X_cnv_pca"`,
`overwrite_existing=False`, `max_output_gib=2.0`, and `random_seed=0`. Keep `svd_solver="arpack"`,
`zero_center=False`, `dtype="float32"`, and infercnvpy `inplace=True` fixed and disclosed. The public output key maps
to infercnvpy's prefix-free `key_added` only after validating an exact `X_` namespace.

The implementation validates/takes one defensive artifact copy, fingerprints the input CNV matrix and axes, calls
the exact public backend on a private AnnData, restores any global RNG state, and independently checks shape,
floating dtype, finite values, nonconstant dimensions and unchanged CNV state. It stores portable PCA provenance
binding the output to the CNV artifact and output fingerprint. `summary` reports upstream source/reference/window
identity, dimensions, per-component variance, memory, seed, warnings, references and versions. `code` is standalone
from OpenBio helpers and returns `(adata_with_cnv_pca, summary_dict)`.

Tests cover typed/tampered artifacts, dense/sparse CNV matrices, bounds, collisions/overwrite, memory, backend
arguments and tampering, input ownership, RNG isolation, real 0.6.1 smoke, strict JSON, compiled code and parity.

Workflow migration now inserts this node only while atomically expanding a legacy CNV Structure whose input is a
direct same-graph edge from an InferCNV node upgraded in the same reviewed plan. It fixes
`output_key="X_cnv_pca"`, `overwrite_existing=False`, `max_output_gib=2.0`, and `random_seed=0`; ARPACK and
uncentered PCA remain fixed in the runtime and are written into the migration disclosure. The old composite did
not serialize Scanpy's data-dependent component count, so migration chooses the reviewed node default
`n_comps=30` and marks that choice for review rather than claiming exact historical equivalence. Current CNV PCA
schemas are migration no-ops; partial output or widget bundles reject the whole workflow before insertion.
