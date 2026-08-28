# OpenBioSingleCellSCENICActivityBinarization — design review

## Decision

**Enhance and retain as an atomic, explicitly exploratory activity-binarization module.** Consume the verified
SCENIC artifact, remove arbitrary AnnData-key coupling, fix a nonzero deterministic seed and one-worker execution,
isolate global RNG, and return typed binary evidence plus thresholds. Do not merge with RSS or TF-module reporting.

## Atomic interface

Inputs, in order:

1. immutable `scenic_result` (`SCENICResultArtifact`);
2. `random_seed`, default 1 and constrained to `1..2^31-1`;
3. optional typed `threshold_overrides` table with exact `regulon,threshold` columns (no free-form dictionary);
4. `max_dense_bytes` guard.

Outputs, in order:

```text
binary, thresholds, summary, code
```

`binary` is an immutable `SCENICBinaryArtifact` containing exact observation/regulon IDs and a defensive boolean
matrix. `thresholds` is a table with
`regulon,threshold,threshold_source,active_cells,total_cells,active_proportion`. No AnnData input, `binary_key`, or
`threshold_key` is needed. A later explicit adapter may attach the artifact to AnnData if a workflow requires that
storage representation.

An in-process single-worker equivalent, default HDT selection, and the upstream strict `>` rule are fixed
implementation details. Overrides are matched exactly to known regulons, finite, defensively copied, and fully
reported. Expert thresholds outside the AUCell score range `[0,1]` are allowed because they have a well-defined
all-on/all-off consequence, but every such regulon is prominently warned and disclosed. Thresholds need not lie
inside a regulon's observed range: the official unimodal rule can legitimately exceed the observed maximum and call
every cell off.

## Invariants and reporting

- The artifact has unique nonblank observation/regulon IDs and finite `[0,1]` AUC values. Zero-sum, constant, or
  at-most-four-distinct-value regulons take the official unimodal branch and are named in warnings; a zero-sum
  regulon therefore receives threshold zero and all cells off rather than being rejected;
  KDE/mixture-degenerate regulons fail with named diagnostics unless explicitly overridden.
- Threshold and binary families equal the complete artifact regulon family; shapes, order, booleans, finite
  thresholds, and active counts are independently verified.
- The wrapper snapshots/restores NumPy global RNG even on error; seed 0 and `None` are rejected; input and artifacts
  are immutable.
- Dense matrix bytes are guarded before artifact materialization or long-form reporting. The conservative estimate
  covers private/canonical cell-by-regulon copies and, whenever any threshold is derived rather than overridden, the
  fixed 1,000-simulation cell-by-simulation dip-test allocation. The binary artifact remains matrix-shaped to avoid
  an avoidable cells-times-regulons public table.
- `summary` follows the strict report schema, labels binarization a heuristic descriptive transformation, reports
  thresholds and bounded active/inactive patterns, provenance, warnings/limitations, references, and dynamic
  versions for OpenBio, pySCENIC, NumPy, pandas, SciPy, and scikit-learn. The fixed 1,000 dip-test simulations and
  conservative peak-byte estimate are disclosed. JSON rejects NaN/Infinity.
- `code` returns `(binary_artifact, thresholds_table, summary_dict)` and reproduces seed/RNG isolation, overrides,
  complete verification, and reporting.

## Scientific boundary

On/off calls are not p-values, Cell identity labels, direct TF-binding evidence, or Condition contrasts. The
threshold is learned across all included cells, so changing the cell composition can change every call. That
dependence and the external SCENIC resource/run provenance must be disclosed.

## Migration and tests

Saved workflows remove `adata`, `activity_key`, `binary_key`, and `threshold_key`; connect `scenic_result`; change a
legacy seed of zero to an explicit nonzero choice; and connect `binary`, `thresholds`, `summary`, and `code`. Existing
unseeded cached calls are not reproducible and must not be relabeled as audited.

Tests must cover fixed synthetic unimodal/bimodal distributions, strict equality-at-threshold behavior, seed 0
rejection, repeated determinism, global RNG restoration on success/error, overrides, constant/small/NaN/Inf inputs,
backend malformed outputs, complete families, artifact ownership, row/memory guards, strict JSON, report wording,
and exact runtime/generated equivalence.

## Cohesion assessment

The module owns one transformation and hides pySCENIC's global RNG, multiprocessing, threshold derivation, override,
validation, and reporting details. Typed matrix output keeps the interface small while preserving all calls.
