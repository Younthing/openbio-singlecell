# OpenBioSingleCellCollecTRIULM — design review

## Decision

**Enhance and retain as one atomic CollecTRI ULM observation-scoring module.** Migrate to decoupler 2.2, make the
resource/license/expression state explicit, and introduce an immutable TF-activity artifact. Do not merge it with TF
ranking: scoring a declared expression matrix against a signed network and comparing scores across annotations have
different statistical assumptions and reasons to change.

## Atomic interface

Inputs, in order:

1. `adata`;
2. `resource` dynamic choice:
   - `local_network`: `network_csv` plus strict `resource_metadata_json`; or
   - `official_collectri`: the directly published `human` resource, affiliation `license`, and
     `allow_network_access=false` by default; non-human runs require a local translated snapshot because 2.2 adds an
     independently changing, uncached HCOP request with hidden translation defaults;
3. explicit normalized expression `source` (`X` or named layer; no Raw snapshot);
4. `complex_policy` (`retain` or `remove`), visible because it changes regulator identity;
5. `min_targets`, default 5, applied after exact measured-feature intersection;
6. advanced `batch_size`, default 250,000, and a bounded memory guard;
7. `overwrite=false` for the fixed owned AnnData keys.

Outputs, in order:

```text
adata, activities, summary, code
```

`activities` is a typed immutable `TFActivityArtifact` containing exact observation IDs, regulator IDs, score and
adjusted-p-value DataFrames, expression-state provenance, canonical network fingerprint, resource metadata, and
backend versions. The convenience AnnData output owns fixed keys such as `collectri_ulm_score` and
`collectri_ulm_padj`; arbitrary key strings are removed from the public interface. The artifact is the supported seam
for `OpenBioSingleCellRankTFActivities`.

`raw=False`, `empty=False`, `verbose=False`, and `tval=True` are fixed implementation details. `organism` is read
from local metadata in local-resource mode. The current AnnData contract carries no mandatory canonical organism or
identifier-namespace identity, so compatibility cannot be independently proven from AnnData alone: exact overlap is
enforced and the caller declaration/verification limitation is reported. No implicit gene-symbol translation or
case folding is allowed in local mode.

In official mode the affiliation `license` value is validated and forwarded, but the report explicitly records that
decoupler 2.2.0 does not use it to filter the fixed Zenodo download. `allow_network_access=false` fails before the
loader; `true` records a direct network retrieval and `cache_state="not_supported_by_decoupler_2_2"`. Generated
code never repeats that download: an official-mode equivalent function requires a caller-materialized local network
whose canonical hash equals the runtime-resolved hash.

## Invariants and reporting

- Input observation/feature IDs are unique, nonblank, and unchanged; the input is immutable.
- Selected expression values must be finite and axis-aligned. Normalized continuous values remain the recommended
  default; count-like, residual, scaled, transformed, or unknown states remain executable expert choices and receive
  specific interpretation warnings. Mutable analysis history never acts as an execution credential.
- Resource columns and metadata are strict; exact duplicate source-target rows are deterministically collapsed and
  counted, conflicting duplicate weights fail closed, and CollecTRI weights are exactly `-1` or `+1`; runtime overlap
  and `tmin` pruning are independently verified.
- Local resource bytes and the canonical normalized edge table receive SHA-256 hashes. Official-loader mode records
  the resolved network hash, direct network access, and the absence of a decoupler 2.2 cache contract.
- The complete ordered expression-feature axis receives a SHA-256 because ULM zero-fills every non-target feature
  and uses `n_features - 2` degrees of freedom. Every observation vector and every retained regulator's complete
  zero-filled weight vector must have nonzero variance before the backend runs.
- The memory guard separately accounts for the full dense feature-by-regulator adjacency, expression conversion
  peak, and two full observation-by-regulator outputs. In decoupler 2.2 `bsize` does not batch dense AnnData or
  DataFrame input and must not be presented as a dense-memory control.
- Backend output has exactly the retained regulators and input observations in canonical order, finite scores,
  adjusted p-values in `[0,1]`, and no extra or missing rows/columns.
- `summary` follows the repository report schema, calls the result exploratory TF-activity evidence, provides
  report-ready strongest positive/negative results and network accounting, and includes dynamic OpenBio, decoupler,
  AnnData, NumPy, pandas, SciPy, and Numba versions, plus Requests for official-loader mode. JSON rejects
  NaN/Infinity.
- `code` defines an importable equivalent function that accepts the same local resource/artifact inputs and returns
  `(output_adata, activity_artifact, summary_dict)` without downloads or hidden global state.

## Scientific boundary

The module estimates cell-level activities. It does not rank annotations, assign cell identities, infer a Condition
effect, or use cells as biological replicates. `padj_ulm` is not a Sample-level Condition p-value. Reports must state
that Technical batch and Condition belong to a later replicate-aware model.

## Migration and tests

Saved workflows replace `split_complexes` with `complex_policy`, remove `activity_key`/`pvalue_key`, select a
licensed and fingerprinted resource mode, and connect the new `activities`, `summary`, and `code` ports. A legacy
workflow must never reinterpret `split_complexes=True` as `remove_complexes=True` without an explicit migration
message because the semantics differ.

Tests must cover the exact decoupler 2.2 signatures/output keys, rejection of 1.x/missing interfaces, license and
network-access policy, local resource hashes, complex policy, post-overlap `tmin`, sparse/dense recommended and
explicitly selected nonstandard input with cautious disclosure, signed toy-network scores, output tampering,
overwrite behavior, memory guards, input
immutability, artifact defensive ownership, strict JSON, dynamic versions, and exact runtime/generated equivalence.

## Cohesion assessment

Resource loading, version drift, network normalization, expression validation, ULM batching, output verification,
and reporting are substantial hidden implementation behind a compact scoring interface. The artifact creates a real
seam because both ranking and future plotting can consume it without learning AnnData storage keys.
