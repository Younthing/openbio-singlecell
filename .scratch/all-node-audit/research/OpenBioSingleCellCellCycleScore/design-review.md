# Cell Cycle Score: module design review

## Decision

**Retain and deepen `OpenBioSingleCellCellCycleScore`.** It should perform one atomic operation: score a declared cell-cycle program on one validated transformed expression source. Do not merge it with normalization, regression, PCA, or trajectory inference.

The seam is:

```text
full-gene transformed expression + versioned organism-specific phase programs
    -> per-cell S score, G2/M score, and phase annotation
```

## Target interface

Inputs:

- `adata`, required;
- `source`, the repository typed expression-source selector, visible and defaulting to `log1p_norm`, while also
  allowing explicit Raw, X, or other layer choices without a workflow-history attestation;
- `gene_set_source`, bounded combo `regev_human_97|custom`, default `regev_human_97`;
- `organism`, bounded combo `human|mouse|other`, default `human`; the bundled resource is identified as a human
  symbol program, while a deliberate nonhuman use is allowed as a warned expert override with no orthology claim;
- `s_genes` and `g2m_genes`, visible only for `custom`, exact comma-separated identifiers;
- `output_prefix`, advanced, default `cell_cycle`;
- `overwrite_existing`, advanced, default `False`;
- `random_seed`, advanced, default `0`.

Remove `use_raw`. A slot boolean is a shallow proxy for expression state. Keep Scanpy's binning and control policy hidden and fixed for the audited version. Bundle the exact 43/54 resource with a resource version, checksum, organism, identifier namespace, and citations; never download it at execution time.

The official-practice targets are at least 20 genes and half of each bundled phase, or at least 10 observed genes per
custom phase. These targets describe evidence quality, not whether Scanpy can execute. Require only one exact
observed gene per phase; report every missing gene and warn below the targets. Reject duplicate identifiers within a
single declared set because duplicate weighting is structurally ambiguous, but allow cross-phase overlap with an
explicit warning because the score remains defined.

Outputs:

- `adata`, a copy with `<prefix>_s_score`, `<prefix>_g2m_score`, and `<prefix>_phase`;
- `summary`, strict JSON;
- `code`, equivalent Python source.

Call Scanpy on an isolated working object, then copy the three verified columns to the requested namespace. This hides Scanpy's fixed output names and prevents accidental replacement of unrelated columns.

## State, summary, and equivalent code

The deep implementation validates only the source properties needed for correct execution: canonical aligned axes
and a real finite matrix. It profiles, but does not certify, count-like/nonnegative/signed/constant state. Full-gene
log-normalized expression remains the recommended report template; Raw, an HVG view, and unverified transformed
state are permitted expert selections with cautious warnings. No mutable analysis-history or Raw-binding record is an
execution dependency. Input `.X`, `.raw`, layers, and existing annotations remain unchanged.

`summary` must include resource identity/checksum, organism and identifier namespace, source descriptor/state, requested/found/missing phase genes and coverage, control size, fixed bin policy, seed, score quantiles/ranges, phase counts/proportions, output keys, methods/references, and Scanpy/AnnData/NumPy/pandas/OpenBio versions. Limitations must say that phase is a supervised expression label and not a continuous clock or replicate-aware comparison.

`code` defines a self-contained function that selects the same Raw/X/layer matrix, validates axes/value shape and
gene availability, calls `sc.tl.score_genes_cell_cycle` on an isolated AnnData with explicit phase lists,
`use_raw=False` and the seed, verifies finite aligned outputs, namespaces them, and returns
`(output_adata, summary_dict)`. It contains the exact bundled list/checksum, carries the same warnings, does not
import OpenBio, and does not access ComfyUI history.

## Migration and tests

Legacy workflows map the old seven-gene defaults to `regev_human_97` and `use_raw=True` to the canonical transformed layer only when provenance proves it exists; otherwise migration adds a blocking warning. Explicit old custom lists migrate to `custom` and retain their exact strings. Existing output collisions require an explicit migration choice rather than silent overwrite.

Tests must cover exact bundled-list identity/checksum, human/custom paths, nonhuman bundled-resource warning, duplicate
identifiers, overlap and low-coverage advisories, Raw/X/layer plus count-like/signed/feature-restricted choices,
dense/sparse matrices, retained structural/non-finite failures, deterministic seed and global RNG isolation,
collision behavior, input preservation, phase/score postconditions, strict JSON without nonfinite values,
self-contained compiled code, and runtime/code equivalence.

## 2026-08-29 repair record

- New category: `openbio/single-cell/annotation`.
- Classification rationale: the atomic analysis assigns supervised per-cell phase scores and a phase label from a declared gene program; it does not infer a continuous trajectory or lineage.
- Merge/delete decision: retain the standalone scoring node. Do not merge it with normalization, regression, PCA, or trajectory reconstruction and do not delete it, because its gene-program validation and annotation outputs form one cohesive operation.
