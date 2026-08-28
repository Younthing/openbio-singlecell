# Keep QC plot metrics on one expression source

Type: task
Status: resolved
Blocked by: 01

## Question

When pre-existing `obs` QC columns are partial or stale, what public behavior guarantees that total counts, detected
genes, and mitochondrial percentages shown together all describe the node's explicit source selection?

## Evidence

- `.scratch/all-node-audit/research/OpenBioSingleCellQCPlots/official-usage.md`
- `.scratch/all-node-audit/research/OpenBioSingleCellQCPlots/design-review.md`

## Acceptance criteria

- A partial metric family cannot be mixed with freshly computed values.
- Mitochondrial unavailability is disclosed rather than represented as biological zero.
- Sparse and dense inputs behave equivalently, and the input AnnData is unchanged.

## Comments

- 2026-08-29: Claimed after the report-contract seam passed its affected suites. The explicit selected expression
  source will be the sole source of the plotted metric family; stale or partial `obs` QC columns are not provenance.
- 2026-08-29: Runtime and generated code now derive totals, detected genes, and mitochondrial percentage atomically
  from the selected matrix and its matching `var['mt']`; zero-library percentages are missing rather than zero.
- 2026-08-29: The QC suite passes (`18 passed`), including stale complete `obs` columns, sparse layer sourcing,
  zero-library disclosure, generated-code parity, and input immutability.

## Answer

Do not infer ownership from existing QC column names or mutable history. The selected expression source is the sole
fact source for the plot; derive the complete metric family in one private Module call, disclose that source for every
metric, and omit mitochondrial evidence when the selected feature metadata cannot define it.
