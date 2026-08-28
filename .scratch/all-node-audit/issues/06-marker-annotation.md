# Cluster marker evidence and annotation

Type: task
Status: resolved

Blocked by: 05

## Initial audit focus

Separate Cluster marker evidence from Condition contrast, remove the broken logreg table contract, keep provisional and curated annotation distinct, and prevent read-only plot mutation.

## Comments

- 2026-08-28: Initial inventory and code audit complete; per-node pre-change documents pending.
- 2026-08-28: Landed and read the required official-usage and design-review records for all eight registered or
  compatibility node identities before
  implementation. The reviewed boundary keeps marker ranking, threshold selection, evidence plotting, provisional
  model annotation, explicit-set ORA evidence, and reviewed label commitment as separate atomic operations.
- 2026-08-28: The marker-to-ORA seam will use two ordinary typed table artifacts: a canonical marker-evidence table
  and the exact tested-gene-universe table. Both carry the same SHA-256 analysis fingerprint; Filter Marker Genes
  must preserve the universe artifact unchanged, and ORA rejects mismatched pairs. This avoids an opaque generic
  object while preventing a same-size/different-gene background from silently changing the null model.
- 2026-08-28: Marker Genes/Filter Marker Genes passed independent final audit after exact dense/sparse constant-gene
  handling, complete-family BH, direct-only filtering, five provenance/content fingerprints, streaming canonical
  hashing, and two-stage working-memory guards were verified. The final marker audit reported 109 passing tests under
  `-W error` with no remaining P0/P1/P2 finding.
- 2026-08-28: UMAP Plot and Marker Expression Plot now emit `plot, summary, code`, construct only the selected private
  data, disclose actual expression/embedding evidence, and isolate seeded violin jitter. Their focused test suite
  passed 26 tests under `-W error`.
- 2026-08-28: The unsafe all-in-one Marker ORA Annotation executor is now a non-executing compatibility shim. The new
  Marker ORA Evidence node accepts a direct Filter Marker Genes table and its paired Marker Genes universe, uses the
  decoupler 2.x `mt.query_set` contract, retains a full evidence grid, and leaves label commitment to Map Cluster
  Annotations.
- 2026-08-28: Batch06 workflow migration is frozen and independently re-run: 53/53 tests pass, mixed/partial schemas
  fail before mutation, and all five packaged workflows are current-format no-ops. CellTypist/ORA generated-summary
  equivalence remains the final implementation gate before this issue can close.
- 2026-08-28: The final annotation audit closed every remaining report/resource boundary: CellTypist 1.7.1 and
  decoupler 2.2.0 runtime/generated outputs and full summaries match; ORA rejects ragged CSV records before pandas
  coercion and fingerprints exact resource bytes for caching; Map's writing-ready results disclose all unmapped,
  missing, category-count, many-to-one, and overwrite outcomes. The stable annotation suite passes 56 tests under
  `-W error`; the integrated Batch06 Python suite passes 148 tests, web migrations pass 53 tests, Ruff/compile pass,
  and all five packaged workflows are migration no-ops. Independent audit reports no remaining P0/P1/P2 finding.

## Answer

Batch06 is complete. Marker ranking and filtering now form a tamper-evident, memory-guarded evidence/universe pair;
plots are read-only report nodes; CellTypist produces provisional, provenance-bound annotation; ORA is explicit-set
evidence rather than automatic label commitment; Map Cluster Annotations is the separate reviewed commitment step.
Every active analysis node exposes its primary result plus strict `summary` and equivalent `code`, while the unsafe
legacy ORA executor remains a non-executing migration shim. Official/API practice, statistical interpretation,
software versions, references, limitations, migration atomicity, and workflow examples are all covered by focused
and cross-node regression evidence.
