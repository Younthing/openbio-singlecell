# OpenBioSingleCellMarkerGenes — design review

## Decision

**Keep and enhance this node as the atomic Cluster marker evidence ranking module. Remove `logreg` from its interface; do not merge marker ranking with filtering, plotting, annotation, or Condition inference.**

The atomic analysis is one complete all-groups-versus-rest ranking on one explicitly selected expression representation. Logarithmized abundance remains the report recommendation, not a history-backed runtime gate. Filtering is a cheap, repeatable selection over an already computed table. Plotting is read-only presentation. Curated annotation requires expert judgment. A Condition contrast uses Samples as replicates and answers a different question.

## Current interface problems

The current module is a shallow adapter around two Scanpy calls and hides several scientific invariants:

- it accepts the repository Raw snapshot even though Raw is counts and Scanpy expects logarithmized data;
- it does not validate expression state, finite values, unique genes, group completeness, observed category count, label collisions, or minimum group size;
- `groups="all"`, `reference="rest"`, `corr_method="benjamini-hochberg"`, `rankby_abs=False`, and `tie_correct=False` are hidden without disclosure;
- optional `pts=False` contradicts the stable table and downstream filter, which still expose prevalence columns but fill them with NaN;
- `logreg` has a different official result interface, yet missing fold-change and p-value fields are silently converted to NaN; binary logistic regression also does not emit symmetric per-group evidence;
- `random_seed` is irrelevant for the three statistical-test methods and exists only because of the incompatible logistic mode;
- `n_genes=100` silently truncates evidence before downstream threshold review;
- the output names `logFC`, `p`, and `p_adj` are terse, and the approximate nature of Scanpy fold change is undisclosed;
- it returns no required `summary` or `code`.

## Proposed interface

Visible inputs:

- `adata`;
- `groupby="leiden"`;
- `method="wilcoxon"`, limited to `wilcoxon`, `t-test`, and `t-test_overestim_var`;
- explicit expression source limited to X or layer, default layer `log1p_norm`;
- `n_genes=100`, where zero means all selected-source genes.

Advanced input:

- `tie_correct=True`, active only for Wilcoxon.
- `max_output_rows`, a bounded performance guard checked against marker-output rows plus universe-output rows before constructing either artifact.
- `max_working_memory_gib=4.0`, an advanced two-stage guard covering identifier/group validation before observation-scale list materialization, bounded matrix scans, the selected-expression temporary copy, group-by-gene backend/prevalence arrays, complete internal ranking, streaming fingerprints, and both outputs.

Fixed and disclosed policy:

- all observed groups;
- each group compared with `rest`;
- positive-score ranking (`rankby_abs=False`);
- `pts=True`;
- Benjamini–Hochberg adjustment;
- internal collision-safe Scanpy result key;
- copy-on-analysis and no AnnData output.

Remove `random_seed`, because no retained method is stochastic. Remove Raw from the source selector because the repository's Raw snapshot is counts. Remove `pts` as a public toggle because prevalence is part of the node's promised result. Do not expose arbitrary `groups`, reference groups, logistic kwargs, correction algorithms, or Scanpy storage keys in the first refactor; those would enlarge the interface into several different analyses.

The current tests cover X/Raw/layer routing, stable legacy columns, Raw feature-count provenance, finite fold changes in the packaged logarithmized-layer path, and logreg-only seed forwarding. They do not cover expression-state rejection, malformed/missing groups, group-size limits, logistic table incompleteness, numeric/range invariants, strict JSON, or generated-code equivalence. Those missing cases define the replacement test surface.

Outputs:

- primary marker-evidence `table`;
- complete ordered tested-gene `universe`;
- strict JSON `summary`;
- equivalent Python `code`.

The two table artifacts share the analysis identity, complete-untruncated-ranking fingerprint, and ordered-universe identity. Each also carries a recomputable current-content fingerprint. This makes the pair cryptographically depend on the scientific ranking result without hashing the entire input expression matrix, and detects numeric table mutation or same-axis cross-run pairing.

The table migration should rename `logFC` to `log2_fold_change_approx`, `p` to `p_value`, `p_adj` to `p_adjusted`, `pct_in_group` to `fraction_in_group`, and `pct_rest` to `fraction_reference`. Update downstream marker filters, plots, workflow templates, examples, and tests atomically.

## Deep module and seam placement

The external seam is “validated grouping plus logarithmized expression in, stable marker-evidence and tested-universe tables plus audit record out.” The module becomes deep by hiding expression-state resolution, count-like counterevidence detection, categorical-label normalization, group/reference validation, minimal temporary AnnData construction, constant-gene neutralization, complete-family BH, stable long-table extraction, fingerprints, performance preflight, finite/alignment checks, report construction, and code rendering behind that small interface.

Introduce one internal marker-table seam shared with `FilterMarkerGenes`: a schema/provenance validator for the canonical marker-evidence DataFrame. This is a real seam because ranking produces the table and filtering consumes it. Keep Scanpy structured arrays and temporary AnnData keys inside the ranking implementation; exposing them would leak adapter details and reduce locality.

The deletion test favors retention. Deleting this module would duplicate expression validation, group/reference setup, statistical policy, structured-array extraction, schema stabilization, and reporting wherever marker evidence is needed. Deleting `logreg` from this interface removes complexity without forcing required marker-test behavior into callers.

## Cohesion and coupling

All retained visible parameters define one ranking question. Upstream coupling to clustering and expression preprocessing is scientifically necessary and should be disclosed through `groupby`, expression state, group sizes, and source provenance—not recreated as inputs. The node must remain uncoupled from Sample/Condition labels, downstream thresholds, marker plots, external marker databases, provisional annotations, and curated annotations.

`tie_correct` is cohesive only for Wilcoxon; the interface should either make it conditionally active or reject a non-default value for t-test modes. Hidden fixed policy earns leverage because most workflows need the same positive all-groups-versus-rest evidence and reports can state it exactly.

## Summary and code contract

The primary table must be deterministic for a fixed input/software stack and contain complete method-appropriate numeric columns. Globally constant genes are neutral hypotheses, not batch-fatal errors: they receive score/log2 fold change zero, p/p-adjusted one, and true `>0` prevalence, while statsmodels BH is recomputed over the complete gene family per group. `summary` contains report-ready methods/results, group sizes, expression source/state or explicit user-verification assumption, ranking truncation, constant/variable counts, score and approximate-fold-change semantics, complete-family correction size, output/memory budgets, fingerprints, top bounded examples, missing/non-finite diagnostics, parameters, warnings, limitations, citations, and software versions.

The summary must use the domain term Cluster marker evidence and must not call the result a Condition contrast, independent-replicate differential expression, a Curated annotation, or biological truth. If p-values are reported, the warning about cell non-independence and their exploratory role is mandatory.

`code` returns `(table, universe)` as numerically equivalent canonical DataFrames and reproduces structural/numeric expression and group failures, nonstandard-state warnings, minimal variable-gene backend execution, constant-gene neutralization, complete-family BH, both budgets, and all fixed Scanpy settings. It may omit OpenBio result wrapping and timing, but it cannot omit state disclosure, tie policy, result-column validation, or log-fold-change labeling.

## Failure and migration plan

- Reject backed/empty inputs, wrong-shape/non-numeric/non-finite expression, duplicate or whitespace-ambiguous identifiers, and unavailable layers. Proven or apparent non-log/count/scaled/residual/negative states continue with explicit warnings and never a claim that logging was verified.
- Reject missing/empty group keys, missing labels, fewer than two observed groups, unused categories, string label collisions, and groups/rest references with fewer than two cells.
- Reject invalid methods, invalid `n_genes`, and irrelevant tie settings rather than silently forwarding them.
- Assert input AnnData immutability, stable group ordering, exact rows per group, unique `(group, gene)` rows, finite promised statistics, prevalence in `[0,1]`, p-values in `[0,1]`, tolerance-aware `p_adjusted >= p_value`, and one-based ranks.
- Add explicit tests proving `logreg` is unavailable and explaining migration. Do not preserve its NaN-shaped legacy output.
- Migrate legacy workflows using Raw to a full-gene `log1p_norm` layer. Existing packaged workflows already use that layer and Wilcoxon, so their scientific path is retained.
- Migrate table-column consumers and serialized workflows together. A temporary alias layer inside `FilterMarkerGenes` would prolong two table interfaces and weaken the seam; prefer one coordinated migration.
- Verify strict JSON, generated-code two-table equivalence, dense/sparse parity for all three methods, exact-extrema detection of all-zero/constant-positive neutral genes (without a tolerance that can erase real small nonzero values), complete-family BH, memory/output preflight before backend, content tampering, same-axis/different-ranking pair rejection, and small clusters.

## 2026-08-29 repair record

- New category: `openbio/single-cell/marker-evidence`.
- Classification rationale: this atomic analysis ranks **Cluster marker evidence** for exploratory annotation review; it is not a replicate-aware **Condition contrast**. This preserves the boundary defined by [ADR-0001](../../../../docs/adr/0001-separate-marker-evidence-from-condition-inference.md) and the terms in [CONTEXT.md](../../../../CONTEXT.md).
- Merge/delete decision: retain the standalone ranking node. Do not merge it with Condition inference, threshold filtering, plotting, or annotation, and do not delete it; each neighboring operation has a different decision, cost, or claim.
