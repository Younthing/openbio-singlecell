# OpenBioSingleCellLianaCommunication — design review

## Decision

**Enhance and retain as one atomic Sample-resolved candidate-communication module.** Exact-lock the audited LIANA
1.9.0 interface (1.10.x is newer but separately unaudited), call the
official `.by_sample` interface, make organism/resource/expression/identity roles explicit, return a typed result directly,
and add strict reporting/code. Do not perform a pooled analysis or claim Condition inference.

Merge the useful output role of `OpenBioSingleCellLianaResults` into this node: writing a hidden `adata.uns` key and
requiring a second extractor is a shallow, tightly coupled interface. Keep plotting separate because figure
selection and rendering have a different reason to change.

## Atomic interface

Inputs, in order:

1. `adata`;
2. `sample_key`, default `sample`;
3. `condition_key`, default `condition`;
4. `identity_key`, default `cell_type`;
5. `annotation_status` (`unknown`, `provisional`, or `curated`), default `unknown`;
6. `organism`, explicit;
7. `method` (`rank_aggregate` or `cellphonedb`);
8. `resource` dynamic choice:
   - `bundled_human`: `resource_name`, default `consensus`; or
   - `local_resource`: explicit `resource_name` plus `resource_csv`;
   both branches require strict `resource_metadata_json` because the public bundled selector does not disclose
   release/license/source-license facts;
9. explicit finite expression `source` (`X` or named layer; no Raw snapshot); normalized full-gene/log1p input is
   the documented recommendation, while state and feature completeness remain unverified report fields;
10. `expression_proportion`, default 0.1;
11. `min_cells_per_identity_sample`, default 5;
12. `permutations`, default 1000, and `random_seed`, default 1337;
13. advanced fixed-by-default `jobs=1` and `max_output_rows`/memory guard.

Outputs, in order:

```text
result, summary, code
```

`result` is a typed, immutable-by-interface and tamper-evident `LianaResult` with stable leading columns
`sample,condition,source,target,ligand_complex,receptor_complex` followed by the exact audited method-native columns.
CellPhoneDB preserves its additional `ligand`/`receptor` subunit and group-statistic fields; rank aggregation does not
fabricate them. Its table is built from the direct `.by_sample(..., inplace=False)` return, defensively owned, and never
stored in or read from an arbitrary caller-owned `adata.uns[result_key]`. Because LIANA 1.9 still mutates the object
passed to `.by_sample` even with `inplace=False`, the adapter always uses and discards a private AnnData copy.

`use_raw=False`, `return_all_lrs=False`, `verbose=False`, and one job are fixed implementation details. The result
metadata records method, native score names/directions, sample/condition/identity roles, resource fingerprint, and
expression state so plotting cannot guess them from columns.

## Invariants and reporting

- Role keys are nonblank and distinct. Cell IDs and current feature identifiers are unique; selected expression is
  finite, aligned, and matched exactly to resource subunit strings. Full-gene normalized/log1p input remains the
  recommendation, while nonstandard state and unverified completeness are warnings. String overlap does not
  independently prove organism or identifier namespace, so the declared resource metadata and upstream feature
  annotation remain an explicit report limitation. Neither history nor Raw/current equality is a hard gate.
- Sample, Condition, and identity labels are complete nonblank scalars; each Sample maps to one Condition; at least
  two Samples are required. Every Sample-by-identity cell count and eligibility decision is independently computed.
- Bundled mode is human-only and rejects `mouseconsensus`. Both resource modes require exact name/version/date/
  organism/namespace/scope/license/source-license-review/citation metadata; local resources require exact
  `ligand,receptor` columns. Raw/local bytes where applicable, canonical rows, and metadata receive SHA-256; no
  online translation.
- Exact LIANA 1.9.0, its declared pandas `<3` runtime range, the `by_sample` wrapper, and the underlying method-call
  parameter names are checked at runtime.
  Output Sample identities, common keys, method-native
  columns/ranges, duplicates, row count, and finite values are verified; Condition is joined by validated mapping.
- No generic `significant` flag is created. Rank-aggregate and CellPhoneDB score meanings/directions are preserved.
- `summary` follows the strict report schema, provides bounded report-ready per-Sample candidate interactions and
  design/resource accounting, explicitly denies Condition inference, and includes dynamic OpenBio, LIANA, Scanpy,
  AnnData, NumPy, pandas, SciPy, and plot-independent dependency versions. JSON rejects NaN/Infinity.
- `code` returns the portable `(table, summary_dict)` equivalent and reproduces role/resource/state validation, the
  exact `.by_sample` call on a private copy, output verification, deterministic sorting, and report construction
  without downloads or caller mutation.

## Scientific boundary

The node estimates expression-compatible sender-receiver candidates separately within each Sample. Permutation
p-values and aggregate ranks are within-Sample cell-label evidence, not replicate-aware Condition tests. Technical
batch, paired designs, and cross-Sample Condition effects require another explicitly Sample-level method. Reports
must state that transcript co-expression is not proof of secretion, binding, spatial contact, or signaling causality.

## Migration and tests

The exact legacy-to-current matrix is:

| Legacy contract | Current contract | Safe translation |
| --- | --- | --- |
| input 0 `adata` | input 0 `adata` | Preserve the link only as the expression object for a fresh run. Never reuse its hidden LIANA table. |
| input 1 `groupby` | input 3 `identity_key` | The string may be offered as the identity role, but execution remains blocked until Sample, Condition, and annotation status are reviewed. |
| input 2 `method` | input 6 `method` | Preserve only exact `rank_aggregate` or `cellphonedb`; the result must be recomputed with LIANA 1.9.0. |
| input 3 `resource_name` | input 7 `resource.resource_name` | May enter the `bundled_human` branch only after explicit human organism plus the ten reviewed metadata fields are supplied. No metadata is invented. |
| input 4 dynamic `source` (`X`, layer, or Raw) | input 8 dynamic `source` (`X` or layer) | Legacy scientific semantics remain too ambiguous for automatic migration; runtime execution of an explicit current `X`/layer choice does not require history or a Raw binding. The legacy Raw branch has no automatic mapping. |
| input 5 `result_key` | removed | Delete; hidden `uns` storage is not part of the new contract. |
| absent | inputs 1–2, 4–5 and 9–15 | `sample_key`, `condition_key`, `annotation_status`, `organism`, resource metadata, thresholds, seed, and guards require explicit reviewed values; defaults are not evidence. |
| output 0 mutated `adata` | output 0 `OPENBIO_LIANA_RESULT` | This is a semantic type change. Rewire only consumers explicitly migrated to the new typed result after a fresh run. |
| absent | outputs 1 `summary`, 2 `code` | Trailing additions for the rebuilt node; never synthesize them for a legacy cached result. |

Therefore a genuine legacy Communication node is an atomic **migration rejection**, not an automatic widget reorder.
Current-schema nodes are no-ops. A later migration may proceed only with an explicit reviewed payload supplying every
new role/resource field, then must replace the node, remove a directly connected legacy Results shim, rewire only
type-compatible consumers, and preserve unrelated graph state atomically. Old pooled cached results cannot be
upgraded to Sample-resolved evidence.

Tests must cover both methods against exact LIANA 1.9.0 signatures with an isolated pandas 2.3 smoke, explicit
rejection of 1.10.x and pandas 3.x, by-Sample calls and its known private-copy mutation, Sample-to-Condition conflicts,
missing/colliding labels, sparse/dense recommended and nonstandard finite input disclosure, too-small strata, human
bundled and local non-human resources, hashes/metadata/complexes, seed/permutations, native result schemas and score
directions, malformed backend rows, row guards, input immutability, strict JSON, scientific wording, and exact
runtime/generated equivalence.

## Cohesion assessment

The module hides version drift, per-Sample partitioning, role and resource validation, method-specific schema checks,
and report construction behind a compact scientific interface. Returning the table directly deletes an artificial
storage/extraction coupling while preserving plotting as a real downstream seam.
