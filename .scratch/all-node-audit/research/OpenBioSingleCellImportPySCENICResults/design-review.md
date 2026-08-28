# OpenBioSingleCellImportPySCENICResults — design review

## Decision

**Enhance and make this the sole supported pySCENIC boundary.** It is an atomic import-and-validation module, not a
second execution engine. Replace live pySCENIC helper imports with a safe parser and return a typed immutable result
artifact consumed by RSS, binarization, and TF-module reporting.

## Atomic interface

Inputs, in order:

1. `adata`, the exact full-gene finite numeric expression object exported to the external workflow;
2. `run_manifest_json`, the bundle manifest with schema `openbio-singlecell/pyscenic-external-run/v1`;
3. advanced `overwrite=false` for the fixed convenience activity key;
4. `max_file_bytes`, `max_adjacency_edges`, `max_regulon_edges`, and `max_dense_bytes` guards.

Outputs, in order:

```text
adata, scenic_result, summary, code
```

`scenic_result` is an immutable typed `SCENICResultArtifact` owning:

- exact ordered observation IDs and cell-by-regulon AUCell DataFrame;
- canonical motif-pruned regulon membership with regulon, TF, sign/context, target, target weight, and motif evidence;
- verified external manifest, file/resource hashes, and software/container provenance;
- expression-state, organism, genome build, and gene-namespace fingerprints.

The convenience AnnData output stores only the verified activity frame at a fixed owned key (for example
`scenic_auc`); it does not invoke `add_scenic_metadata`, create arbitrary keys, or expose an overwrite-prone
`activity_key` string. Downstream nodes consume `scenic_result`, not hidden AnnData state.

## Manifest and parser invariants

- The manifest has a strict allowlisted schema and finite JSON values. It records pySCENIC exactly 0.12.1, exact
  declared arboreto/ctxcore/Python versions plus any additional canonical package-version entries, three CLI
  argument arrays, an explicit expression-state label/fingerprint, and
  every consumed bundle-member SHA-256 plus organism/build/namespace/license/citation. An immutable container image
  digest, explicit `grn`/`aucell` seeds, workers, method, mode, and AUC threshold are preferred; official 0.12.1
  defaults or a non-container execution are accepted and prominently warned rather than treated as malformed output.
- Every declared member path is a normalized relative path confined beneath the manifest directory. Actual bytes are
  streamed into SHA-256 and must match before parsing. Symlinks on consumed paths, path traversal,
  path/size/mtime identity changes, network URLs, and implicit downloads are prohibited. Unrelated undeclared files
  are ignored and never treated as content-bound provenance.
- The GRN adjacency is a strict official `TF,target,importance` CSV with finite nonnegative importance, expression-axis
  membership, unique TF-target edges, and its own row guard. The TF list has unique canonical identifiers; motif
  annotations must contain pySCENIC's five required tabular fields with strict row widths and finite numeric values;
  ranking databases must carry Feather-v2 (`ARROW1`) framing. These checks validate safe import structure without
  importing the obsolete external pySCENIC/ctxcore runtime.
- AUCell CSV cell IDs must be unique and in the exact AnnData order; AUC is finite, bounded to `[0,1]`, named, and
  defensively copied. The exact ordered gene axis and a canonical sparse/dense value fingerprint bind the current
  AnnData object to the externally exported finite numeric expression matrix. The bundle expression export is itself
  a strict official cell-by-gene CSV; its axes and values are parsed and independently required to match that
  fingerprint. The first/index column header may be blank or any canonical label because pandas ignores that label
  when loading the index. Non-raw or negative-valued states are retained as explicit expert choices and warned.
- Motif CSV header depth and required columns are exact. Target-pair lists use `ast.literal_eval`; Context accepts
  only the official textual grammar `frozenset(<set-literal>)` and applies `ast.literal_eval` to the inner set.
  Explicit byte/item/depth limits apply; general `eval`, pickle, object arrays, and arbitrary imports are prohibited.
- Consolidated final regulon names equal AUC columns exactly. pySCENIC's final name is `TF(+)` or `TF(-)` according to
  context. Repeated motif-supported TF-target pairs are unioned at maximum weight and their motif IDs retained;
  finite zero weights are retained and warned, while negative weights, exact duplicate motif rows, ambiguous
  regulation context, or targets outside the bound gene axis fail.
- The input AnnData is immutable; output observation and feature identity/order are unchanged.

## Reporting and code

`summary` follows the standard report schema and includes method stages, report-ready bounded regulon/activity
results, complete resource/manifest accounting, warnings, limitations, references, and dynamic versions for
OpenBio, Python, AnnData, NumPy, pandas, and the externally declared pySCENIC stack. External versions are
clearly labeled declared-and-hash-bound rather than locally imported. JSON rejects NaN/Infinity.

`code` defines a standalone safe bundle importer returning `(output_adata, scenic_result, summary_dict)`. It reproduces
hashes, schema checks, literal parsing, exact identity joins, and report construction and never downloads resources.

## Scientific boundary

Import validates provenance and structure; it does not rerun, repair, normalize, or reinterpret the external model.
It cannot prove that an external computation was honest, only that declared artifacts are internally consistent and
content-bound. SCENIC evidence remains exploratory and cell-level, not a Condition contrast.

## Migration and tests

Saved import workflows replace separate loom/regulon paths with the manifest, remove `activity_key`, and connect
`scenic_result`, `summary`, and `code`. Old loom/regulon files without a complete manifest fail closed; metadata must
not be invented from filenames. The
old Run node's temporary results are unrecoverable and cannot be migrated.

Tests must cover a valid official-shaped AUCell CSV/regulon/manifest fixture, every hash mismatch,
short/extra/reordered/duplicate cells,
regulon/AUC mismatch, malformed multi-index CSV, malicious literal strings, oversized nesting/edge counts, NaN/Inf,
unknown targets, resource/build mismatches, manifest version drift, omitted optional provenance/defaults,
transformed numeric expression, overwrite policy, defensive ownership, no
pySCENIC import/eval/network access, input immutability, strict JSON, dynamic versions, and exact generated/runtime
equivalence.

## Cohesion assessment

The importer hides complex untrusted-file parsing, identity binding, provenance verification, and artifact
normalization behind two scientific inputs plus bounded safety controls. The immutable artifact is a real seam with
three independent consumers,
providing leverage and keeping pySCENIC's obsolete runtime out of downstream modules.
