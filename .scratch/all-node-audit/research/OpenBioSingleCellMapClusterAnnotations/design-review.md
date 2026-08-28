# Map Cluster Annotations: module design review

## Current module

`OpenBioSingleCellMapClusterAnnotations` is a small in-process transformation: parse a JSON dictionary, stringify
both sides, map one `obs` column, optionally preserve unmapped values, infer a categorical, overwrite the requested
column on an `AnnData.copy()`, and append generic analysis history. The schema is:

```text
adata, groupby="leiden", mapping_json="{}", output_column="cell_type", keep_unmapped=True
    -> adata
```

There are no direct behavior, schema, migration, provenance, category-order, nonmutation, or failure tests.

## Correctness and interface problems

1. JSON duplicate names are silently last-write-wins, and non-standard `NaN`/`Infinity` constants are accepted by
   Python's default decoder.
2. Every JSON value is stringified. `null`, booleans, numbers, arrays, and nested objects silently become labels.
   Empty/blank keys and values are accepted.
3. Distinct native source values can collapse after `astype("string")`; no ambiguity check protects the mapping.
4. Mapping keys absent from the data are silently ignored. Declared-but-unused categorical levels are not
   distinguished from truly unknown keys.
5. `keep_unmapped` hides two materially different outcomes behind a boolean: retain raw cluster IDs in an annotation
   column or turn unmatched cells into missing values. The default mixes curated-looking labels and raw numeric
   cluster IDs without a warning.
6. Under preservation, a mapped target can collide with a preserved source label and silently merge unrelated cells.
7. Source missing values, observed unmapped levels, unused declared levels, and mapping-generated missing values are
   not separately counted.
8. `pd.Categorical` infers the output categories. The public order is not tied to the source categorical order and
   is not documented or tested.
9. Duplicate mapping targets legitimately merge clusters but the merge is not disclosed. Cluster counts before and
   after mapping are absent.
10. Any existing output column is overwritten without opt-in. `output_column == groupby` destroys the source cluster
    assignment in the returned object.
11. The full copy prevents direct input mutation, but backed input behavior/cost is unspecified and there are no
    invariants proving X, layers, raw, embeddings, graphs, axes, and existing annotations are preserved.
12. Generic history stores parameters but no machine-readable Provisional-versus-Curated status or annotation
    provenance suitable for downstream review.
13. The node name/result can be read as establishing cell-type truth even though it examines neither Cluster marker
    evidence nor external references.
14. There is no strict `summary` or equivalent `code` output.

## Decision: keep, deepen, and enhance

Keep the node and deepen it at its existing seam. It is the correct atomic operation for applying a caller-reviewed,
cluster-level label map to one observation annotation. The deletion test supports keeping it: without this module,
strict JSON parsing, label canonicalization, mapping coverage, categorical order, collision safety, provenance, and
reporting would recur across workflows or ad-hoc scripts.

Do not merge it with Leiden. Clustering is an unsupervised reusable partition; annotation changes with evidence,
reviewer judgment, nomenclature, and desired granularity. Do not merge it with Marker ORA or marker ranking: those
produce Cluster marker evidence, while this module records a decision made from evidence. Do not merge it with
CellTypist: automated labels are Provisional annotation and have model/confidence provenance. Do not merge it with
cross-AnnData annotation transfer: this node maps values within one aligned object and needs no identity join.

No adapter is warranted. JSON, pandas, and AnnData are in-process dependencies, and there is only one implementation.
Private helpers for parsing/canonicalization may form internal seams, but the node's interface is the test surface.

## Target interface

Keep the first three scientific widgets recognizable and replace the ambiguous boolean with constrained policies:

- `adata`
- `groupby="leiden"`: source cluster column
- `mapping_json`: non-empty strict JSON object from source label to annotation label
- `output_column="cell_type"`: destination annotation column
- `unmapped_policy = error | preserve_cluster_label | set_missing`, default `error`
- `annotation_status = provisional | curated`, default `provisional`
- advanced `overwrite_existing=False`

Outputs:

```text
adata, summary, code
```

`annotation_status` is scientific provenance and must remain exposed. `curated` means only that the caller asserts
expert review; the implementation cannot verify that assertion. `unmapped_policy` is also scientific because it
changes which cells receive biological labels. `overwrite_existing` is advanced operational policy. Do not expose
string-canonicalization flags, categorical ordering, missing-source handling, JSON decoder choices, mapping-size
limits, or copy mode; these are fixed implementation invariants that make the module deep.

Do not add marker thresholds, ontology lookup, fuzzy matching, case folding, regex mapping, per-cell overrides, or
automatic "best label" selection. Each would introduce a different analysis or an external seam.

## Fixed validation and mapping contract

### Input and collision validation

- Require an in-memory, non-empty AnnData and a present nonblank `groupby` column with at least one non-missing
  observed source level.
- Require nonblank `output_column`; always reject `output_column == groupby` so cluster evidence remains auditable.
- If the output column exists, fail unless `overwrite_existing=True`. Overwrite replaces the whole column; this node
  does not merge old and new annotations.
- Parse with `object_pairs_hook` to reject duplicate JSON members and `parse_constant` to reject `NaN`/`Infinity`.
  Require a non-empty object, bound text/entry counts, and require every target value to be a nonblank string.
- Trim source keys and target labels, reject blank values, and reject raw or normalized duplicate keys. Reject
  nonscalar source values and distinct source values that collapse to one canonical string.
- Allow duplicate targets because they intentionally merge clusters. Record the reverse target-to-source mapping.

### Coverage, unknown, missing, and unused semantics

- Canonicalize observed values separately from declared categorical levels.
- Reject mapping keys unknown to both observed and declared source levels.
- Accept mapping keys for declared-but-unused levels only as reported unused entries; they do not create output
  categories or count as an effective mapping.
- Require at least one observed source level to be explicitly mapped.
- Default `unmapped_policy="error"` rejects any observed unmapped level and lists its cell count.
- `preserve_cluster_label` retains canonical source labels for unmapped observed levels, warns that the output mixes
  biological annotations with cluster IDs, and rejects collisions between preserved labels and mapped targets.
- `set_missing` turns cells in unmapped observed levels into true missing values and reports them.
- Missing source values always propagate as true missing values under a fixed hidden policy. They are never converted
  to text or treated as an unmapped category. Report their count and warn when nonzero.

### Deterministic categorical order

Do not use JSON object order; RFC 8259 defines objects as unordered. Determine source-level order as follows:

1. categorical source: declared `cat.categories` order, canonicalized and filtered to observed levels;
2. noncategorical source: lexically sorted canonical observed labels;
3. resolve each source level through the map or unmapped policy;
4. deduplicate targets in first occurrence along that source-level order;
5. build an explicit unordered `pd.Categorical` with exactly those realized categories.

This is invariant to cell row order, keeps Scanpy Leiden's declared category order in the common path, supports
many-to-one merges, drops unused source/target categories from the result, and never invents a biological order.

## Atomicity, copy semantics, and provenance

The node changes exactly one observation annotation plus OpenBio provenance. It must preserve observation/variable
order and identity, shape, X, layers, Raw snapshot, obsm/varm, graphs, all unrelated `obs`/`var`, and input `uns`.
Create a full output copy only after all validation succeeds so failures leave no partially modified object. The input
must remain byte/value-equivalent at the public scientific attributes.

In addition to normal analysis history, store a bounded machine-readable annotation entry under OpenBio metadata,
keyed by `output_column`, containing:

- `annotation_status` and `curation_assertion="caller_declared"` when status is curated;
- source/output columns, unmapped and overwrite policies;
- normalized effective mapping and unused declared-level mappings;
- source and output category order plus ordered flags;
- observed level/cell counts, missing/unmapped/mapped counts, and many-to-one merges.

On explicit overwrite, replace that column's current annotation-provenance entry while generic analysis history keeps
the chronological record. Do not claim the mapping was expert-reviewed unless the caller selected `curated`, and even
then state that the software did not verify the review.

## Summary and code contract

`summary` must be strict JSON and include method/package references and dynamic software versions. Its core report
should state, in writing-ready language:

- how many cells and observed clusters were processed;
- how many clusters/cells were explicitly mapped, preserved, set missing, or already missing;
- the output annotation categories and cell counts;
- which source clusters were intentionally merged into one label;
- whether the result is Provisional annotation or caller-declared Curated annotation;
- whether an existing column was overwritten.

Key results should include source dtype/ordered/categories, declared unused levels, mapping entry counts, unused
mapping keys, per-source and per-output counts, missingness, output category order/dtype, collision checks, and the
annotation-provenance entry. Warnings/limitations must say that one label is assigned to every cell in a source
cluster, within-cluster heterogeneity is hidden, cluster resolution/preprocessing affects the result, no marker or
ontology evidence was checked, and provisional labels are not suitable as Curated annotation for formal Condition
contrast without review.

References should include Luecken/Theis annotation practice, pandas mapping/categorical documentation and citation,
AnnData documentation and citation, and RFC 8259/Python JSON documentation. Record Python, openbio-singlecell,
pandas, and AnnData versions dynamically.

`code` should be a self-contained function with resolved arguments. It must implement the same strict JSON parser,
label/coverage/collision validation, deterministic category construction, missing/unmapped policy, overwrite rule,
input nonmutation, and annotation provenance, and return a new AnnData. It need not recreate ComfyUI result wrappers
or timestamped plugin history.

## Failure and migration plan

The old widget layout is:

```text
groupby, mapping_json, output_column, keep_unmapped
```

Migration must rewrite the fourth value rather than let a boolean enter a Combo:

- `keep_unmapped=true` -> `unmapped_policy="preserve_cluster_label"`
- `keep_unmapped=false` -> `unmapped_policy="set_missing"`
- add `annotation_status="provisional"`
- add `overwrite_existing=false`

Append `summary` and `code` outputs after the existing AnnData slot. Migrate named widget values as well as positional
arrays and make the transform idempotent in root graphs and subgraphs. A connected or exposed legacy boolean cannot
be converted dynamically into a string policy; fail migration with a clear message rather than silently disconnect or
reinterpret it. Existing workflows that relied on silent output overwrite will now fail until the user explicitly
opts in. Strict parsing may also reject previously accepted nested/non-string/duplicate/unknown mappings; these are
intentional safety failures and should name the offending entry.

## Verification plan

- Exact mapping for string, integer, boolean, categorical, and ordered-categorical source labels; reject scalar
  string-conversion collisions and nonscalar values.
- Invalid/empty/non-object JSON; duplicate raw and normalized keys; `NaN`/`Infinity`; null/numeric/boolean/nested or
  blank target labels; bounded input and entry counts.
- Unknown keys, declared-but-unused keys, observed unmapped keys under all three policies, all-missing source, partial
  missing source, and zero effective mappings.
- Many-to-one targets, identity mappings, mapped-target/preserved-source collisions, deterministic order under cell
  row permutation, declared category order, and noncategorical lexical order.
- Existing output with default error and explicit overwrite; unconditional rejection of `output_column == groupby`;
  provenance replacement/history behavior.
- Input unchanged; output shape/axes/X/layers/raw/obsm/graphs/unrelated annotations preserved; backed input rejected.
- Strict JSON summary, bounded key results, dynamic versions/references, Provisional-versus-Curated wording, and no
  truth/causal/formal-inference claims.
- Generated code compiles and matches success output, categories/order/missingness, provenance, nonmutation, and major
  failure modes.
- Workflow migration for positional/named widgets, connected/exposed legacy boolean rejection, subgraphs, output-slot
  preservation, and idempotence.

## Cohesion and coupling assessment

The redesigned module has a small interface but owns substantial implementation complexity: interoperable parsing,
canonical identity, categorical semantics, mapping coverage, collision prevention, deterministic ordering, copy
safety, domain status, provenance, and reporting. That depth gives callers leverage without coupling annotation to
clustering, evidence generation, automated classifiers, ontology services, or Condition inference. The interface is
also the complete test surface; private parser/order helpers remain internal seams rather than new public adapters.
