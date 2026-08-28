# OpenBioSingleCellMarkerORAEvidence — design review

## Decision

**Create this node as the active, atomic explicit-set ORA evidence module. Keep the old
`OpenBioSingleCellMarkerORAAnnotation` ID only as a registered, non-executable migration shim for one compatibility
period. Do not mutate AnnData or commit labels here.**

The legacy node combined expression coercion, per-cell ORA, a second cell-level group test, and unconditional first
label assignment. Those operations have different inputs and uncertainties and cannot be repaired under its
`adata -> adata` contract. The new ID makes the breaking `marker evidence + universe -> evidence table` type change
honest. `MapClusterAnnotations` remains the reviewed decision seam.

The deeper legacy assessment is in `../MarkerORAAnnotation/design-review.md`; this document fixes the public seam of
the new node before code is added.

## Atomic interface

Inputs, in order:

1. canonical filtered marker-evidence `table` (`OPENBIO_SINGLE_CELL_TABLE`);
2. exact tested-gene `universe` (`OPENBIO_SINGLE_CELL_TABLE`);
3. explicit `resource_csv` under the ComfyUI input boundary;
4. strict `resource_metadata_json` containing nonblank `name`, `version`, `date`, `organism`,
   `identifier_namespace`, `scope`, `license`, and `citation` strings;
5. advanced `source_column` and `target_column` names;
6. `min_targets` applied after resource/universe intersection, while report accounting keeps all in-universe
   resource targets distinct from the subset belonging to tested/surviving sources;
7. `min_overlap` used only for candidate eligibility;
8. `max_p_adjusted` used only for descriptive significance/candidate flags.

Outputs, in order:

```text
table, summary, code
```

No random seed, expression source, AnnData, groupby, output label column, arbitrary background size, alternative,
or correction constant is exposed. Those inputs either belong upstream or would change the named method.

## Upstream evidence seam

Marker Genes emits a canonical marker table and complete universe table under provenance schema v2. Filter Marker
Genes validates that direct pair and returns the same universe object. ORA accepts only that direct Filter output
paired with the direct Marker universe; unfiltered Marker tables and Filter-to-Filter chains are rejected. ORA must
revalidate:

- producer/schema/artifact roles;
- the shared `analysis_fingerprint` for axes/parameters/group membership;
- the shared `ranking_fingerprint` for the complete, non-truncated all-group by all-gene numeric ranking;
- the shared `universe_fingerprint` for ordered tested-gene identity;
- each artifact's recomputed `content_fingerprint` for its actual current DataFrame;
- the Filter artifact's `upstream_content_fingerprint`, which binds it to its direct Marker table;
- unique nonblank identifiers and consecutive universe ranks;
- unique `(group, gene)` marker rows;
- selected genes are inside the universe;
- positive marker direction and non-truncated/filter provenance suitable for explicit-set testing.

Current-content fingerprints use canonical schema-v2 row serialization with exact strings and hexadecimal floats.
Two independent table wires keep generic preview/CSV tooling available. Their cryptographic binding prevents
numeric marker tampering, same-count/different-gene universes, and cross-pairing artifacts from different complete
rankings. A filtered subset cannot reconstruct the complete ranking, so ORA cross-checks its pinned ranking identity
across the two direct artifacts rather than claiming to recompute it.

## Resource seam

Resolve the CSV through the existing input-directory boundary, fingerprint bytes with SHA-256, parse without silent
type coercion, and validate metadata before importing decoupler. Exact duplicate pairs may be collapsed and counted;
blank identifiers, conflicting metadata, zero universe overlap, or no source surviving `min_targets` fail.

Identifiers are trimmed only. Case folding, symbol translation, synonym expansion, organism inference, network
downloads, and bundled marker defaults are prohibited. Resource suitability remains a caller-reviewed declaration
and is disclosed as a limitation.

## Canonical evidence table

One deterministic row per tested `(group, source)` pair:

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

`overlap_genes` is a deterministic JSON array string. Counts obey all contingency identities. Numeric values are
finite and p-values lie in `[0,1]`. Rows sort by stable group order, within-group adjusted/raw p-value, descending
effect and overlap; lexical source is display-only. Complete scientific ties share rank and remain multiple
candidates. Full rows are never dropped by `max_p_adjusted` or `min_overlap`.

## Summary and code

`summary` follows the repository analysis-report schema and must include report-ready methods/results, exact
parameters, upstream fingerprints and selection policy, resource path/hash/metadata/accounting, group/resource
counts, overlaps, unresolved/ambiguous groups, bounded top alternatives, Fisher/Haldane/BH/decoupler references,
strict JSON values, installed software versions, and explicit Provisional/Cluster-marker-evidence limitations.

`code` defines one importable function that performs the same input/resource validation, recomputes marker/universe
current-content and ordered-universe identities, cross-checks pinned analysis/ranking provenance when artifacts are
provided, checks the decoupler 2.x API/version, performs per-group `query_set` calls, independently verifies
contingencies, both BH families, ordering/tie/candidate logic, and the scientific summary. It may omit ComfyUI
wrappers and wall-clock history but cannot omit fingerprints, resource metadata, or postconditions.

## Failure and compatibility design

All preflight validation completes before the optional backend is called. No partial table is returned after a
malformed backend result. Missing/incompatible decoupler, wrong API, wrong row set, non-finite results, contingency
or BH mismatch, identifier/resource errors, and invalid thresholds fail with bounded actionable messages.

The old `OpenBioSingleCellMarkerORAAnnotation` remains registered with its legacy inputs and `adata` output so saved
workflows remain inspectable. Its execution always raises a migration-required error explaining the new wiring. The
workflow migrator may add an idempotent note but must not change its ID, ports, or links. It must never execute the
defective 1.x algorithm.

No old ORA workflow can be automatically translated because it contains neither an explicit marker set nor the
exact tested universe, and its output type is incompatible. No example workflow should add the new node until a
licensed, versioned, fixed-hash resource can be supplied.

## Verification obligations

- hand-calculated contingency, Fisher, Haldane, within-group BH and global BH fixtures;
- decoupler absent/1.x/missing API/malformed output cases through a narrow fake adapter;
- analysis/ranking/universe/current-content fingerprints, numeric tampering, cross-pair, same-size/different-gene,
  selected-outside-universe and direct-producer mismatch failures;
- resource metadata, columns, missing/duplicates, intersection, minimum-size and hash accounting;
- ties, unresolved groups, minimum-overlap and significance status without row deletion;
- input immutability, deterministic ordering, strict JSON, bounded summaries and equivalent generated code;
- registration of active Evidence plus deprecated shim, shim execution failure, and migration port/link preservation.

## Cohesion assessment

The module has one reason to change: the definition and reporting of explicit selected-set over-representation.
Resource parsing, optional-dependency drift, contingency reconstruction, correction families, provenance and result
validation are hidden behind a small scientific interface. Marker ranking/filtering, plotting, model annotation,
curated mapping and Condition inference remain separate modules, preserving locality and interpretability.

## Final summary-builder decision

The validated ORA core now owns one canonical parameter map, warning list, and software-version map. A pure,
standalone-embeddable builder turns those diagnostics into the complete analysis-report schema. Both the node and the
generated function use that builder, so methods/results, dynamic resource reference and SHA-256, references,
limitations, warnings, target accounting, and software versions cannot drift between two assembly paths. The
generated function returns `(evidence_dataframe, summary_dict)` and is tested for exact summary equality, strict JSON,
input immutability, empty selected groups, tampering, and malformed-backend failure equivalence.
