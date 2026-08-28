# Unify named-graph base validation

Type: task
Status: resolved
Blocked by: 02

## Question

How should graph consumers share structural validation while retaining method-specific requirements and standalone
generated-code equivalence?

## Evidence

- `.scratch/all-node-audit/research/OpenBioSingleCellDiffusionMap/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellDiffusionMap/design-review.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellPAGA/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellPAGA/design-review.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellDPT/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellDPT/design-review.md`

## Acceptance criteria

- Explicitly stored sparse zeros are canonicalized before diagonal/self-loop checks in every path.
- Base graph failures are consistent across Leiden/UMAP, trajectory, and Schist consumers.
- Method-specific validation remains local to the method that needs it.

## Comments

- 2026-08-29: Claimed after coherent QC sourcing passed. The first tracer bullet targets the observable trajectory
  discrepancy for explicitly stored zero diagonal coordinates while preserving legitimate off-diagonal zero distances.
- 2026-08-29: Canonicalization now rejects only nonzero diagonal values, removes zero diagonal coordinates, and retains
  off-diagonal zeros when distance structure requires them. Runtime and generated code use the same embedded helper.
- 2026-08-29: The trajectory suite passes (`29 passed`), including fingerprint equivalence, input storage immutability,
  generated-code parity, preserved off-diagonal zero distances, and rejected nonzero self-loops.

## Answer

Keep one narrow sparse-matrix canonicalization policy inside the existing trajectory graph Module: validate numeric
values, remove semantically empty diagonal coordinates, optionally preserve off-diagonal zeros, then apply
method-specific graph checks. Do not merge the complete Leiden/trajectory/Schist resolvers: their distance, directed,
isolate, memory, and fingerprint policies genuinely differ and would create a wide policy Interface.
