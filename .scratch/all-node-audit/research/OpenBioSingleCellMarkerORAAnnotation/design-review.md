# MarkerORAAnnotation: module design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. Current workflows use
> `OpenBioSingleCellMarkerORAEvidence` followed by explicit reviewed annotation commitment.

## Decision

**Replace the current compound annotation behavior with an enhanced atomic Marker ORA evidence module. Preserve the
capability, but remove automatic AnnData label mutation. Do not merge marker selection, ORA statistics, and label
commit into one node.**

In keep/merge/delete/enhance terms:

- **keep** the scientifically useful operation: testing explicit cluster marker sets against curated cell-type marker
  sets;
- **enhance** its contracts, statistical implementation, evidence output, provenance, failure behavior, `summary`,
  and `code`;
- **delete** the existing internal sequence that performs per-cell ORA, a second cell-level t-test, and unconditional
  top-label assignment;
- **do not merge** ORA with `MarkerGenes`/`FilterMarkerGenes`; they define which genes are Cluster marker evidence;
- **do not merge** ORA with `MapClusterAnnotations`; the former supplies evidence and the latter records a reviewed
  Provisional or Curated annotation decision.

If a breaking node identity is acceptable, introduce `OpenBioSingleCellMarkerORAEvidence` and deprecate
`OpenBioSingleCellMarkerORAAnnotation`. If the existing ID must be retained, its display name and schema still need a
versioned breaking migration because the present output type encodes the wrong module responsibility.

The deletion test supports retaining the ORA module: removing it entirely would redistribute resource validation,
universe construction, contingency accounting, multiple-testing scope, overlap evidence, and provenance into callers.
Those are substantial details worth hiding behind a small interface. The current annotation behavior does not pass
the deletion test: `MapClusterAnnotations` already owns the label-commit seam more coherently.

## Current depth, cohesion, and locality problems

The node's apparent single output hides four operations:

```text
expression-source coercion
  -> per-cell top-expression ORA
  -> group-versus-rest source ranking
  -> one-label-per-cluster mutation
```

That is a shallow module with poor locality. To understand one `ora_annotation` value, a caller must know the Raw
snapshot conversion, decoupler 1.x top-5% and 20,000-background defaults, seeded tie breaking, discarded per-cell
p-values, the second conservative t-test, ranking sort order, and that significance/direction were ignored. None of
this is visible in the output.

It is also tightly coupled to volatile dependency names and storage layout. The current implementation writes
decoupler's intermediate `.obsm` matrices into a copied AnnData, extracts them back out, and then retains those
matrices even though the public contract never explains them. It changes Raw snapshot dtype/storage for backend
compatibility and reports the current `adata.n_vars` even when `.raw` supplies a different feature axis.

There are no direct tests, no compatible decoupler version constraint, no resource shipped at the default path, and
no evidence artifact from which a reviewer could reproduce or challenge the label.

## Design it twice

### Design A: compatibility-shaped observation scoring

One possible refactor would keep AnnData input, declare an expression source, run `dc.mt.ora`, extract
`score_ora`, call `dc.tl.rankby_group`, filter positive/significant rows, write provisional labels, and add
`summary`/`code`.

This would resemble the official tutorial, but it has a large public surface and weak scientific depth:

- expression source/state, `n_up`, `n_bm`, background, tie policy, ORA tail, group-ranking method, reference,
  significance threshold, and output storage all become coupled choices;
- sparse cell rows make top-expression set membership tie-sensitive;
- two statistical layers obscure what the final group label means;
- group-level p-values treat cells as observations and can look more inferentially precise than warranted;
- overlap genes and the cluster's explicit marker evidence remain indirect;
- the current decoupler 2.2.0 `mt.ora` release has a top-feature boundary defect documented in the research audit.

This design is acceptable only as a separately named **cell-level marker-set scoring** capability. It should not be
the implementation of a Cluster marker evidence annotation node.

### Design B: explicit cluster-marker set ORA

The preferred design accepts already filtered positive marker evidence plus the exact tested-gene universe, validates
one marker resource, and runs `dc.mt.query_set` once per cluster. It returns the full cluster-by-resource evidence
table and no biological label mutation.

This is the deeper module:

- callers choose marker ranking/filter criteria in the nodes that own those choices;
- ORA owns exactly one null model and one contingency-table analysis;
- the universe and every overlap are inspectable;
- there is no expression matrix, Raw snapshot, normalization, per-cell tie breaking, or second cell-level test in the
  interface;
- the backend seam is narrow and can be tested against known contingency tables;
- candidate review and label commit remain downstream.

Choose Design B. It has fewer public knobs, stronger information hiding, more leverage per input, and a smaller test
surface while preserving the scientifically useful capability.

## Required upstream contract before this module is implementable

The current MarkerGenes `TableResult` is not self-sufficient: it contains only the returned top `n_genes` and an
`input_genes` count. An exact ORA universe cannot be reconstructed from that count.

Before implementing this refactor, establish one of these explicit interfaces, in priority order:

1. a typed Cluster marker evidence result containing the selected marker table, exact tested-gene universe, universe
   SHA-256, upstream method/reference/direction/filter parameters, and identifier namespace; or
2. a second `universe` TableResult output from MarkerGenes with one unique gene per row and a shared analysis
   fingerprint, consumed alongside the filtered marker table.

Do not infer the universe from an unrelated AnnData connection. `adata.var_names` may represent HVGs while the
marker test used Raw snapshot genes; anndata layers share `.var` but `.raw` can have a different variable axis. A
count match is not identity. Do not use the resource target union or a numeric genome-size default.

The ORA node should reject generic/unfiltered tables unless they satisfy the versioned marker evidence contract. It
must not quietly treat the first 100 ranked genes—including weak, negative, or nonsignificant rows—as a selected set.

## Target public interface

Recommended conceptual schema:

```text
marker_table       # filtered Cluster marker evidence
universe           # exact tested-gene universe from the same analysis
resource_csv
resource_metadata  # typed artifact or strict metadata object
min_targets=3
min_overlap=2
max_p_adj=0.05
    -> table, summary, code
```

The precise typed artifact may evolve, but the scientific information cannot be optional.

### Exposed parameters

| Input | Visibility | Contract |
|---|---|---|
| `marker_table` | primary | Required selected positive marker evidence with `group`, `gene`, upstream method and filter provenance. |
| `universe` | primary | Required exact unique gene universe from the same marker analysis; identity/fingerprint must match. |
| `resource_csv` | primary | Required marker-set resource; no nonexistent packaged default. |
| `resource_metadata` | primary | Required name/version or retrieval date, organism/taxon, identifier namespace, scope, license, and citation. Prefer a typed resource artifact over many strings. |
| `min_targets` | advanced | Default 3 for continuity with the official cell-marker examples; applies after universe intersection. |
| `min_overlap` | advanced | Default 2 for candidate-eligibility only; does not remove tested evidence rows. It prevents a one-gene match from becoming an automatic candidate flag. |
| `max_p_adj` | primary | Default 0.05; controls only candidate/significant flags and report counts, never deletion of full results. |

`source_column` and `target_column` may remain advanced migration inputs for arbitrary CSVs, but a validated resource
artifact should normalize them to `source` and `target` and hide those names. If raw CSV remains the immediate
interface, validate nonblank distinct column names before reading the data.

### Hidden fixed policy

Keep the following implementation choices hidden and disclose them in `summary`/`code`:

- exact one-sided over-representation alternative: `greater`;
- Haldane-Anscombe correction: 0.5;
- contingency orientation and use of a single explicit universe;
- exact duplicate collapse before `min_targets`;
- resource/universe intersections and strict missing-identifier handling;
- BH within each cluster across retained sources, plus a globally adjusted evidence column;
- stable output sorting and display-only lexical tiebreaker;
- full evidence retention regardless of significance;
- no random seed because explicit set ORA is deterministic;
- no expression normalization, layer choice, Raw snapshot access, or AnnData mutation.

Do not expose `alternative`, `ha_corr`, `n_bg`, `seed`, `use_raw`, group-ranking method, `reference="rest"`, or a
numeric background. Exposing them would either permit a different method under the same node name or preserve
obsolete implementation coupling.

## Atomic execution seam

```text
validate marker-evidence and universe identity
  -> resolve, fingerprint, and validate resource plus metadata
  -> canonicalize exact set members and account for missing/duplicates
  -> intersect resource with universe and apply min_targets
  -> call narrow decoupler query_set adapter once per group
  -> reconstruct/validate contingency and overlap evidence
  -> add declared adjustment families and candidate-status fields
  -> build immutable table result, strict summary, and equivalent code
```

An internal decoupler adapter is justified because decoupler is a volatile optional dependency and its 1.x/2.x
interface and method semantics differ. The adapter should expose one private operation such as
`query_feature_set(selected, network, universe_size, min_targets) -> frame`. Production uses decoupler 2.x; a fake
and independent SciPy gold calculations make the module testable. Do not make the adapter a public plugin socket.

All validation and resource parsing should occur before invoking the adapter. Results should be accumulated in local
objects, validated completely, and committed to output wrappers only at the end. Failure therefore cannot leave a
partially annotated AnnData or partial resource state.

## Canonical evidence-table contract

One row per tested `(group, source)` pair, with stable columns:

```text
group
source
selected_count
universe_count
set_size_before_universe
set_size_in_universe
overlap_count
overlap_genes
a
b
c
d
log_odds_ratio
p_value
p_adj_within_group
p_adj_global
rank_within_group
positive_enrichment
passes_min_overlap
significant_within_group
candidate_status
```

`overlap_genes` should be a deterministic sorted JSON array in a machine-readable table cell if the result contract
supports nested values; otherwise use a delimiter-safe serialized JSON string, not an ambiguous semicolon join.
Every numeric column must be finite or represented as true missing with an explained warning. Counts must be integer,
nonnegative, and satisfy the contingency identities.

The row order is deterministic: source-independent group order from the input contract, then
`p_adj_within_group`, `p_value`, descending `log_odds_ratio`, descending overlap count, and lexical source solely as a
display tiebreaker. `rank_within_group` may share a rank for equal scientific keys; lexical order must not turn a tie
into different biological ranks.

Suggested `candidate_status` is descriptive and nonexclusive at the cluster level:

- `eligible`: positive, within-group adjusted-significant, and `overlap_count >= min_overlap`;
- `not_significant`, `nonpositive`, or `insufficient_overlap` for failed gates;
- if multiple rows share the complete best scientific ranking key, all are `tied_best` and the cluster is ambiguous;
- if no row is eligible, the cluster is unresolved.

Do not create a `confidence` column. Odds ratio, p-value, adjusted p-value, and overlap are distinct evidence fields,
not probabilities of cell identity.

## Cluster-label output semantics

The node may summarize the best evidence, but it must not write `.obs["ora_annotation"]` or output a mapping that
looks curated. The summary should use phrases such as “top ORA-supported candidate” and preserve alternatives.

A label is eligible for review only when enrichment is positive, overlap meets the explicit threshold, and the
within-cluster adjusted p-value passes `max_p_adj`. Even then:

- it remains a Provisional annotation candidate;
- exact ties remain an array of candidates;
- near alternatives and their evidence remain visible;
- a resource can lack the correct cell type, so a top rank is not proof;
- unresolved clusters remain unresolved rather than inheriting the first source.

`MapClusterAnnotations` is the downstream decision module. It accepts a caller-reviewed mapping and explicitly
records whether the user declares it provisional or curated. This seam permits combining ORA evidence with marker
plots, domain knowledge, other models, and nomenclature review without rerunning statistics.

## Resource, identifier, missing, duplicate, and tie policy

### Marker evidence

- Require `group` and `gene`; reject missing/blank values rather than silently drop scientific evidence.
- Require unique `(group, gene)` members. Exact duplicates are an upstream contract error, not extra evidence.
- Allow a gene in multiple groups.
- Require every selected gene in the exact universe; list offending values on failure within a bounded preview.
- Require positive marker direction or a producer declaration that selection already enforced direction.
- Reject a truncated top-N input without explicit selection/filter provenance.

### Universe

- Require a nonempty unique string identifier set and stable namespace/species metadata.
- Reject missing/blank/duplicated identifiers and fingerprint mismatch with marker evidence.
- Preserve the exact identity set in provenance through a bounded preview plus count and SHA-256; do not place an
  unbounded list in the summary JSON when a typed artifact owns it.

### Resource

- Reject missing/blank identifiers and declared metadata conflicts.
- Trim surrounding whitespace only; do not case-fold or translate IDs.
- Collapse exact `(source, target)` duplicates before set-size calculations and report counts.
- Intersect targets with the universe, retaining before/after counts per source.
- Apply `min_targets` after that intersection and report removed sources.
- Globally zero overlap or zero surviving sources is a hard failure. Per-group zero overlap is unresolved evidence.

### Ties

- The explicit set test has no expression-rank ties.
- Upstream marker ties belong to the marker selection contract; a top-N producer must include all cutoff ties or
  report a declared tie policy.
- ORA result ties remain multiple candidates. Stable lexical sorting is presentation only.

## Statistical contract

For every group use one common universe and explicit sets:

```text
a = selected and in resource set
b = not selected and in resource set
c = selected and not in resource set
d = neither selected nor in resource set, within universe
```

Call decoupler 2.x `mt.query_set` with `alternative="greater"`, `n_bg=len(universe)`, `ha_corr=0.5`, and
`tmin=min_targets`. Independently reconstruct `a:b:c:d` and overlap genes and assert that the returned corrected log
odds ratio matches the contingency table within a documented numeric tolerance. This guards against wrong universe,
orientation, or dependency drift.

Keep raw p-values. Use decoupler's per-call BH result as `p_adj_within_group`. Calculate `p_adj_global` across every
tested group/source pair with SciPy's BH implementation and report both families. Candidate flags use
`p_adj_within_group` because the immediate question is which resource labels have evidence for one cluster; reports
must say that global significance may differ.

Do not interpret the marker-ranking p-values and ORA p-values as independent layers of confirmatory inference. Marker
selection is data dependent, cells are reused across group-vs-rest marker comparisons, and marker sets overlap. The
result is exploratory Cluster marker evidence for annotation review.

## `summary` contract

Return strict JSON through the audit-wide analysis-report schema. Required contents:

- `schema_version`, node ID/display name, operation, run status, and exact resolved parameters;
- a methods paragraph naming explicit selected-set ORA, contingency orientation, one-sided Fisher exact,
  Haldane-Anscombe 0.5, within-group/global BH families, and candidate gates;
- references to official decoupler APIs/source, decoupler paper, Fisher, Haldane, Benjamini-Hochberg, the actual
  supplied marker resource, and the background-bias paper;
- upstream marker provenance: operation/method/reference/direction, all filtering criteria, selected count by group,
  tested-universe count/fingerprint, and producer versions;
- resource provenance: requested/resolved path, SHA-256, metadata/citation/license, input/valid/duplicate/missing rows,
  source/target counts, overlap losses, before/after set-size distributions, and sources removed by `min_targets`;
- key results: tested pair count, positive/significant/eligible row counts, per-group selected/unmatched counts,
  unresolved/ambiguous groups, and bounded top candidates including effect size, raw/adjusted p-values, overlap
  count/genes, and closest alternatives;
- warnings/limitations: identifier/resource suitability is caller-declared, enrichment cannot validate upstream
  markers, marker sets overlap, clustering/selection/resource context changes results, within-cluster heterogeneity is
  hidden, and cell evidence is not Sample-level Condition inference;
- dynamic Python, openbio-singlecell, decoupler, SciPy, pandas, NumPy, Scanpy, and anndata versions.

No NaN/Infinity, NumPy scalar, tuple key, or unbounded full result table may leak into JSON. Full evidence belongs in
the table output; summary arrays must have explicit totals and truncation flags.

## Equivalent `code` contract

The text output should define one self-contained importable function that returns the same evidence DataFrame plus a
plain summary dictionary. It must:

- contain resolved scientific defaults and explicit resource metadata;
- validate the marker/universe provenance and all identifier invariants;
- hash and parse the resource with the same missing/duplicate/intersection accounting;
- require a compatible decoupler 2.x `mt.query_set` API and record `importlib.metadata.version("decoupler")`;
- call `query_set` per group with the exact alternative/background/correction/minimum size;
- reconstruct contingency/overlap fields, both BH families, stable shared ranks, and candidate statuses;
- perform the same result postconditions and strict JSON normalization;
- avoid ComfyUI, plugin-only result classes, Raw snapshot mutation, expression preprocessing, and hidden network access.

It must not be a narrative pseudocode snippet. Tests should compile and run it on fixtures and compare every
scientific table column and normalized summary field, excluding only runtime/timestamp wrapper metadata.

## Failure behavior

Fail before the adapter is called for:

- missing/wrong table contracts, mismatched producer/universe fingerprints, missing/duplicate/blank identifiers,
  selected genes outside the universe, no groups, or no selected genes in a group;
- invalid/nonregular resource path, unsupported extension, missing/identical source-target columns, invalid metadata,
  hash/read failure, duplicate-collapse ambiguity, zero universe overlap, or no source surviving `min_targets`;
- organism or identifier-namespace conflict between marker evidence, universe, and resource;
- invalid `min_targets`, `min_overlap`, or adjusted-p threshold;
- unavailable decoupler, unsupported major/API, or a backend result missing required rows/columns or containing invalid
  counts/nonfinite/out-of-range p-values;
- contingency mismatch, BH mismatch, duplicate group/source rows, or nondeterministic scientific ties disguised as
  one result.

Errors should name the node operation, relevant group/resource/source where safe, expected invariant, observed count,
and installed decoupler version. Never silently fall back from the explicit-set design to old `run_ora`, current
observation-level `mt.ora`, a fixed 20,000 background, or unadjusted/source-first labeling.

No output wrapper is constructed until the complete evidence and report pass postconditions. Input tables/resources
remain unchanged on every failure.

## Migration plan

The old and target schemas are scientifically incompatible, so migration must be explicit:

1. Add a compatible decoupler optional extra, preferably `decoupler>=2.1.4,<3`, because `mt.query_set` with explicit
   `alternative` is required. Record the exact imported version; do not support both 1.x and 2.x semantics under one
   result schema.
2. Extend MarkerGenes to provide the exact tested universe and shared producer fingerprint; make FilterMarkerGenes
   preserve that artifact/provenance.
3. Introduce the evidence table schema and `summary`/`code`; verify it before changing workflow nodes.
4. Replace `min_n` with `min_targets`. Remove `adata`, `groupby`, `use_raw`, `seed`, and `output_column` from the ORA
   analysis interface. They cannot be safely auto-translated into an explicit marker set/universe.
5. Do not migrate an existing `ora_annotation` column into a reviewed annotation. Preserve old data files as legacy
   provenance only; require rerun and review.
6. Rewrite example workflows to connect MarkerGenes/FilterMarkerGenes plus the universe artifact to ORA evidence,
   then route a human-reviewed mapping to MapClusterAnnotations.
7. If workflow migration cannot synthesize the new required inputs, load the legacy node in a clearly blocked
   `migration required` state with actionable wiring instructions. Do not guess selected markers from current
   expression or silently disconnect downstream nodes.
8. Preserve output slot identity only within a node ID whose output types remain compatible. Because `adata -> table`
   is a type change, a new `OpenBioSingleCellMarkerORAEvidence` ID plus deprecation of the old ID is safer than
   pretending the change is append-only.

A one-release compatibility shim may keep the old node loadable only to explain migration. It should not execute the
scientifically defective algorithm or pin users indefinitely to decoupler 1.x.

## Verification obligations

- Schema: exact typed inputs/outputs, no generic object ports, `table, summary, code` order, and no AnnData mutation
  port in the evidence node.
- Versioning: decoupler absent, 1.9.x, compatible 2.x, missing `mt.query_set`, malformed adapter output, and dynamic
  software version reporting.
- Gold statistics: hand-calculated 2x2 tables covering positive, null, zero cells, Haldane correction, one-sided
  Fisher p-values, within-group BH, and global BH; compare with SciPy and decoupler.
- Universe: exact identity/fingerprint match, same count/different genes rejection, selected gene outside universe,
  resource targets outside universe, fixed-20,000 regression guard, empty/one-gene universes, and Raw/HVG mismatch
  fixtures.
- Resource: missing path/columns/metadata, null/blank IDs, whitespace, exact duplicates, same gene across sources,
  organism/namespace mismatch, zero overlap, per-source pruning, and min-target boundary after intersection.
- Marker input: filtered positive table, unfiltered/top-N rejection, missing group/gene, duplicate within group,
  same gene across groups, empty group, upstream cutoff ties, unmatched selected markers, and stable group order.
- Results: every contingency identity, exact overlap arrays, positive/nonpositive direction, minimum-overlap gate,
  significant/nonsignificant gates, exact ties retained, unresolved groups, shared ranks, stable row order under input
  permutations, and no label confidence fiction.
- Atomicity: inputs unchanged; adapter is not called after preflight failure; no partial result or annotation state;
  bounded errors and summaries for large resources.
- Reporting: strict JSON, complete parameters/references/hashes/versions, adjustment-family names, bounded top
  candidates/alternatives, and required Provisional/Cluster-marker-evidence limitations.
- Generated code: compiles, runs independently, produces equivalent evidence/summary, and rejects the same malformed
  fixtures.
- Workflow migration: root graphs and subgraphs, named and positional widget state, connected legacy output,
  idempotence, and explicit blocked migration when required inputs cannot be synthesized.

## Cohesion and coupling assessment

After the redesign the module has one reason to change: the definition and reporting of explicit-set marker ORA. Its
interface exposes only inputs that change that scientific test. decoupler API drift, set canonicalization,
contingency accounting, p-value correction, resource hashing, and report serialization are hidden behind the module,
improving depth and locality.

The only necessary upstream coupling is a versioned Cluster marker evidence/universe contract; that is semantic
coupling, not incidental storage coupling. The only downstream coupling is the stable evidence table consumed by
review/reporting. Label application remains a separate decision seam, and expression/Raw snapshot state disappears
from this node entirely.

## References

- [Official `decoupler.mt.query_set` API](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.mt.query_set.html)
- [Official versioned `query_set` source](https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/mt/_query_set.py)
- [Official network validation/pruning source](https://github.com/scverse/decoupler/blob/v2.2.0/src/decoupler/pp/net.py)
- [Official decoupler changelog](https://decoupler.readthedocs.io/en/stable/changelog.html)
- [Official current single-cell tutorial](https://decoupler.readthedocs.io/en/stable/notebooks/scell/rna_sc.html)
- Badia-i-Mompel P, et al. *Bioinformatics Advances*. 2022;2:vbac016.
  [doi:10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016).
- Fisher RA. *Journal of the Royal Statistical Society*. 1922;85:87-94.
  [doi:10.2307/2340521](https://doi.org/10.2307/2340521).
- Haldane JBS. *Annals of Human Genetics*. 1956;20:309-311.
  [doi:10.1111/j.1469-1809.1955.tb01285.x](https://doi.org/10.1111/j.1469-1809.1955.tb01285.x).
- Benjamini Y, Hochberg Y. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300.
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
- Timmons JA, Szkop KJ, Gallagher IJ. *Genome Biology*. 2015;16:186.
  [doi:10.1186/s13059-015-0761-7](https://doi.org/10.1186/s13059-015-0761-7).
- Franzén O, Gan L-M, Björkegren JLM. *Database*. 2019;2019:baz046.
  [doi:10.1093/database/baz046](https://doi.org/10.1093/database/baz046).

## Final replacement reporting seam

The replacement implementation no longer assembles a short generated diagnostic report beside a separate runtime
analysis report. Its validated core canonicalizes parameters, warnings, resource provenance, and software versions;
one pure builder constructs the complete strict-JSON scientific summary for both callers. Generated execution returns
the evidence table and that exact summary. This keeps the legacy node's migration boundary intact while making the
active evidence node's reproducibility promise testable without ComfyUI timing/history state.

## Registry disposition

Keep the stable ID for deserialization, but set `is_deprecated=True`, put “Migration Required” or “Retired” in the
display name, and make the schema description identify the non-executing replacement path. Tests must assert all
three public signals as well as fail-before-read execution. This registry-only clarification does not add outputs or
restore the unsafe analysis.
