# OpenBioSingleCellFilterMarkerGenes — design review

## Decision

**Keep and enhance this node as an atomic, table-native marker-evidence selection module. Do not merge it into MarkerGenes and do not replace it with Scanpy's AnnData-mutating filter.**

Ranking and selection have different costs and decisions. MarkerGenes performs the expensive statistical ranking under one expression/group contract; FilterMarkerGenes lets an analyst revise transparent thresholds without recomputing statistics. Merging them would couple exploratory display/selection choices to the statistical analysis and reduce reuse. The Scanpy helper stores NaN placeholders in AnnData and lacks the adjusted-p criterion, so it is not an equivalent adapter for the desired primary result.

## Current interface problems

- The node accepts any `TableResult` with four coincidentally named columns; it does not require MarkerGenes provenance or the full canonical schema.
- `pandas.to_numeric(errors="coerce")` converts malformed values to NaN and silently removes them. This hides the current broken `logreg` contract by returning an empty table.
- It does not validate finite thresholds or finite/ranged input statistics.
- It does not check duplicate `(group, gene)` rows, positive ranks, stable ordering, or group-label consistency.
- Its terse legacy column names obscure approximate fold change and adjusted-p semantics.
- It reports neither per-criterion exclusion nor groups left with zero markers.
- The difference between inclusive OpenBio comparisons and strict Scanpy filtering is undisclosed.
- It has no required `summary` or `code` output.

There are currently no focused threshold, malformed-table, or copy/ordering tests for this node. Registration and demo coverage only prove that the node exists and that one default workflow path executes; they do not verify its scientific interface.

## Proposed interface

Visible inputs:

- canonical marker-evidence `table` produced directly by MarkerGenes;
- its bound complete tested-gene `universe`;
- `min_log2_fold_change=1.0`;
- `min_fraction_in_group=0.25`;
- `max_fraction_reference=0.5`;
- `max_p_adjusted=0.05`.

No advanced parameters are needed. Comparison operators are fixed and inclusive. Do not expose missing-value policy, arbitrary query strings, absolute fold change, sorting, top-N, or Scanpy storage keys. Missing statistics are contract failures, not a configurable filtering choice.

Outputs:

- filtered marker-evidence `table`;
- the original tested-gene `universe` unchanged;
- strict JSON `summary`;
- equivalent Python `code`.

The four thresholds remain visible because each changes the scientific selection. Their defaults should be described as conventional starting values inherited from current workflows, not generally valid biological constants.

## Deep module and seam placement

Use the shared internal canonical marker-artifact-pair validator proposed for MarkerGenes. This is a real seam with multiple adapters/callers: MarkerGenes creates the pair, while FilterMarkerGenes and ORA evidence consume it. The validator owns column names, dtypes/ranges, tolerance-aware adjusted-p invariants, key uniqueness, rank invariants, direct/filtered producer policy, shared complete-ranking identity, ordered-universe identity, and recomputation of each artifact's current-content fingerprint. This creates locality for future table-contract changes.

The filter module becomes deeper by hiding provenance/schema validation, numeric/range checks, four criterion masks, overlap-aware diagnostics, per-group retention accounting, ordering preservation, report construction, and code rendering behind five simple inputs.

The deletion test supports a separate node: without it, every workflow or UI consumer would duplicate threshold masks, handling of invalid tables, and disclosure of why markers were excluded. A generic “filter table” node would move marker semantics to every caller and lose leverage.

## Cohesion and coupling

All four visible thresholds contribute to one selection predicate over a marker-evidence table. The node is intentionally coupled to the canonical MarkerGenes table-plus-universe interface and its provenance, but not to AnnData, expression matrices, clustering implementations, plots, annotation databases, or Condition models. It accepts only a direct MarkerGenes table: Filter→Filter chaining is rejected so each threshold revision starts from the same ranked evidence and remains independently auditable.

Keep filtering read-only and deterministic. Preserve original row order, upstream ranks, complete-ranking fingerprint, and the universe artifact; do not rerank retained genes because rank remains the upstream all-gene comparison. The filtered table gets a new current-content fingerprint and records the direct upstream table's content fingerprint. A filtered rank could be added as a differently named derived column only if a concrete downstream consumer needs it; no such seam exists now.

## Summary and code contract

`summary` must state the exact inclusive predicate and report input/retained counts globally and per group, individual criterion failures, groups with zero retained rows, upstream truncation/provenance, and threshold limitations. It should distinguish “no passing rows” from invalid/missing upstream statistics.

The report must call the output filtered Cluster marker evidence. It must not say that thresholds validate cell types, identify a Curated annotation automatically, control all exploratory multiplicity, or support a Condition contrast.

`code` takes marker and universe pandas DataFrames, validates both canonical schemas, gene membership, adjusted-p invariants, and source-operation content identities embedded in the code, creates four named masks, applies their conjunction without reordering, resets the index, and returns `(filtered, universe_copy)`. Its two table values and malformed-input failures must match the node; plugin timing/result wrappers may be omitted.

## Failure and migration plan

- Reject non-`TableResult` inputs, any non-direct-MarkerGenes source operation (including FilterMarkerGenes), legacy/incomplete columns, mismatched ranking/universe identities, recomputed current-content mismatches, duplicate keys, invalid ranks, nonnumeric/non-finite values, tolerance-breaking `p_adjusted < p_value`, out-of-range p/fraction fields, and empty/whitespace-ambiguous group/gene identifiers.
- Reject bool/non-numeric/non-finite thresholds and percentage/p-value thresholds outside `[0,1]` before examining rows.
- A direct MarkerGenes table is never structurally empty because its group-by-rank grid is complete. For a valid nonempty input with zero passing rows, return a canonical empty output with a prominent result statement. That empty filtered output is terminal and cannot be fed back into this node.
- Migrate the legacy column names in the same change as MarkerGenes, downstream plots, workflow generators, serialized examples, and tests.
- Add a regression proving logistic-shaped NaN statistics fail loudly rather than silently producing an empty result.
- Test input/universe immutability, inclusive equality at every threshold, per-group diagnostics, overlapping failures, preserved order/ranks, direct-only enforcement, content-number tampering, same-axis/different-ranking cross-pair rejection, strict JSON, and generated-code two-table equivalence.
- Do not add a compatibility branch accepting both legacy and canonical names indefinitely; that would create two public interfaces and spread migration complexity across every consumer.
