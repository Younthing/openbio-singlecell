# cNMF Rank Survey: module design review

## Decision

Retain `OpenBioSingleCellCNMFRankSurvey` as the first atomic analysis in the staged cNMF workflow, preserve its public schema, and replace the former OmicVerse in-memory adapter with a narrow file-backed adapter for the method authors' exact `cnmf==1.7.1` release.

Inputs remain `adata`, `source`, `components_min`, `components_max`, `n_iter`, `num_highvar_genes`, and `random_seed`. `source` now openly supports current `X`, a named layer, or explicit Raw through the existing dynamic expression-source interface. Outputs remain `run`, `k_metrics`, `summary`, and `code`. The `run` wire remains the concrete typed `OPENBIO_CNMF_RUN`; no string path or generic Python object is accepted.

## Why the staged seam remains justified

Rank surveying is expensive and produces reusable state. Final K and density are analyst decisions normally made after inspecting the complete family of metrics. Keeping the survey and consensus separate avoids recomputing all restarts for each decision while keeping each node cohesive:

- Survey owns count validation, official preparation, complete restart execution, combination, and K metrics.
- Consensus owns one chosen K/density decision, final consensus, result loading, and AnnData annotation.

Deleting the typed seam would force expensive factorization and completeness validation back into every downstream attempt. The seam is therefore scientifically meaningful, not accidental infrastructure.

## Deep module boundary

The public interface is small; the internal module hides:

1. exact package-version and signature guards;
2. private H5AD materialization and on-disk axis/fingerprint verification;
3. unique managed-temporary creation and ownership transfer;
4. official prepare/factorize/combine/statistics calls;
5. path containment, symlink/reparse, regular-file, and replacement checks;
6. restart/merged/normalized artifact validation;
7. narrow known-warning compatibility handling;
8. strict report and equivalent-source construction.

The filesystem is a local-substitutable internal seam. Users do not choose paths, names, caches, workers, or resume behavior. The `cnmf` package is the true external dependency and is guarded at the adapter edge.

## Managed run contract

`CNMFRun` owns the `TemporaryDirectory` object rather than merely storing a path. Its invariants are:

- metadata and metrics are private-slot values exposed through getter-only properties, and a construction-time canonical snapshot fingerprint binds their complete typed contents;
- every public access, lock acquisition, and downstream consumption revalidates current metadata/metrics against that snapshot, including source state/evidence/advisories and every K metric;
- portable generated-code interoperability requires an exact class-owned artifact marker/schema and complete snapshot/lifecycle/content validation; generic duck-typed proxies and instance-forged markers are rejected;
- the canonical root and official run directory are unique, private, live, and non-symlinked;
- every backend path resolves strictly under that root;
- `locked_backend()` verifies liveness and serializes mutable upstream cache/result operations;
- `close()` is explicit and idempotent in equivalent Python, removes the whole owned tree, clears backend reachability, and causes later access to fail;
- `weakref.finalize` is cleanup fallback only;
- node-graph cleanup occurs when the engine releases the cached typed run; Consensus does not auto-close because one run may feed multiple sibling decisions, and the retained-until-cache-release limitation is reported;
- a cleanup exception leaves the run closed but keeps the temporary owner and finalizer attached, so fallback cleanup remains possible;
- when cleanup is secondary to an existing construction, context-body, or report exception, it is attached with `add_note()` and cannot mask the primary failure; standalone explicit `close()` cleanup failure remains outward and retryable;
- construction exceptions clean immediately; ownership transfers to `CNMFRun` only after all survey postconditions pass;
- the artifact cannot be pickled or transported across processes.

The node never returns or accepts a bare resource path. This prevents stale, shared, substituted, or user-controlled results from crossing the scientific boundary.

## Execution and validation design

The survey performs these steps transactionally:

1. Resolve the declared X/layer/Raw source and materialize a private source-aligned AnnData without mutating input AnnData. Raw uses its own complete observation/feature axes and becomes the owned output feature axis; it is not intersected with current variables.
2. Validate only backend-required numeric/non-negative/axis/dimension invariants; disclose non-integer/transformed-source, small-restart, requested-HVG, and resource-estimate concerns as expert advisories.
3. Create the unique temporary root and private counts H5AD; reopen and verify exact axes/content fingerprint.
4. Import exact `cnmf==1.7.1` and guard `cNMF`, `prepare`, `factorize`, `combine`, `consensus`, and `load_results` signatures.
5. Construct `cNMF(output_dir=root, name=name)` and call `prepare` with every resolved parameter explicit.
6. Validate replicate-parameter/seeds and normalized-count artifacts.
7. Call one complete `factorize(worker_i=0,total_workers=1,skip_completed_runs=False)`.
8. Verify every expected per-restart spectrum.
9. Call strict `combine(...,skip_missing_files=False)` and verify every merged spectrum.
10. Load normalized counts from the official file and call statistics consensus once per K with `build_ref=False`.
11. Validate the full metric family, bind read-only metadata/metrics to a complete canonical artifact snapshot, transfer temporary ownership, and return.

Partial output is never a valid run. Resume/skip controls are deliberately hidden because this node has no authenticated checkpoint contract.

## Expert-open validation boundary

Runtime hard gates are restricted to conditions required for calculation or safe staged reuse:

- exact backend version/signatures and complete official files;
- in-memory nonempty AnnData with unique aligned identifiers;
- for explicit Raw, existing nonempty unique Raw observation/feature axes; Raw/current axis equality, overlap, values, history, and historical integrity are deliberately not gates;
- numeric, finite, non-negative input and nonzero cell totals required by TPM/final normalized usages;
- K at least two and restart count at least two, the minimum domain for the reported silhouette family; nonempty realized backend features;
- strict typed-run liveness, private-path containment, immutable hashes, and result axes;
- the collision/overwrite contract.

The following are report warnings, not runtime gates: non-integer-but-non-negative input, provenance-confirmed transformed input, zero-total genes, fewer than 100 restarts, K above cells or realized/positive-variance features, requested HVGs above available genes, a conservative memory estimate above the 2 GiB reference envelope, and the analyst's eventual K/density interpretation. There are no arbitrary K, restart, or HVG upper limits. This keeps the tool open to expert use without weakening backend, identity, or filesystem safety.

## Warning compatibility adapter

Warnings are captured only around the upstream calls known to emit them. A warning is consumed only when category, exact stable message, and source file/function boundary match the audited 1.7.1 defect. Unknown warnings are re-emitted with their original category/location and therefore fail under `-W error`. Known warning counts and descriptions are included in summary and run metadata.

## Reporting and source output

The strict summary contains the entire K family and restart accounting, not only a preferred row. References identify Kotliar et al., the method-author repository/guide, PyPI 1.7.1, and scientific packages; software identifies `cnmf`, AnnData, NumPy, pandas, SciPy, and scikit-learn. The report does not cite OmicVerse as the backend.

Generated code repeats the same file-backed lifecycle and validation path. It returns an independently usable owned run rather than a soon-deleted directory. Runtime and generated parity compares metrics, official artifact axes/fingerprints, metadata, liveness, and cleanup—not private plugin history.

The summary records `source="raw"`, Raw `source_features`, current features as context, and the absence of a Raw/current equivalence claim. Consensus from that run returns an AnnData on Raw's complete feature axis so its gene-program matrices remain shape- and identifier-aligned.

## Migration matrix

The already-reviewed legacy-to-staged workflow migration remains unchanged because public schemas and the `OPENBIO_CNMF_RUN` wire are unchanged:

- legacy `adata`, count source, K bounds, restart count, HVG count, and seed route to Rank Survey;
- Rank Survey `run` routes to Consensus `run`;
- legacy selected K, density threshold, neighborhood size, top-gene count, and overwrite policy route to Consensus;
- removed worker/directory/name/GPU/hard-classification controls are not reintroduced and connected unsafe legacy controls continue to fail closed;
- explicit Raw is valid in the current node interface; the already-reviewed legacy `raw` migration remains an explicit warning rather than silently asserting that an older workflow intended the new full-Raw-axis contract;
- bypassed legacy nodes continue to preserve only safe AnnData fan-out.

No workflow migration file change is required for the backend replacement.

## Verification plan

- exact fake-backend call sequence and every explicit argument;
- real tiny `cnmf==1.7.1` CPU smoke through prepare, four restarts, combine, and K statistics;
- dense/CSR/CSC and X/layer/explicit-Raw sources, with input immutability; the positive Raw case uses a Raw feature axis distinct from current variables and verifies runtime/generated full-axis parity;
- hard failure for non-finite/negative data, zero-total cells, invalid axes, incomplete jobs, and impossible backend dimensions; advisory-only coverage for non-integer/transformed sources, zero-total genes, small restart families, excess requested HVGs, and resource estimates;
- missing/replaced/misaligned per-restart, merged-spectrum, normalized-count, parameter, and seed files;
- metadata/metric rebinding and low-level mutation rejection, plus portable generated-artifact identity and fake-proxy rejection;
- live-directory, explicit-close, idempotent-close, destructor fallback, retryable cleanup failure, and primary-exception-preserving construction/context/report cleanup tests;
- absolute escape, `..`, symlink/reparse, and out-of-root backend-path rejection;
- known-warning capture plus proof an unrelated warning is not swallowed under `-W error`;
- strict JSON and standalone generated-code parity;
- registration and existing cNMF migration tests.

## Cohesion assessment

The result is high-cohesion because the node performs one analysis: a complete rank survey. It is low-coupling because downstream code sees a typed run with scientific metadata, while file layout, lifecycle, external quirks, and cache details stay behind the adapter. Retaining the node is preferable to merging it with consensus or general preprocessing.

## 2026-08-29 expression-state consistency repair

The standalone cNMF resolver now recognizes `scale_to_layer` provenance as `scaled`, matching the package resolver and
triggering the existing cautious expert-source disclosure. Runtime and generated rank-survey functions are
contract-tested to return the same state and evidence.
