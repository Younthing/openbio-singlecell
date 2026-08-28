# OpenBioSingleCellRankTFActivities — design review

## Decision

**Enhance and retain as an atomic exploratory annotation-ranking module.** Migrate to decoupler 2.2, consume the
typed TF-activity artifact, preserve the complete regulator family, and remove generic Condition-like semantics. Do
not merge with CollecTRI ULM because ranking can operate on an already validated score artifact and has separate
statistical assumptions.

## Atomic interface

Inputs, in order:

1. `adata`, used only for aligned annotation labels;
2. immutable `activities` (`TFActivityArtifact`);
3. `annotation_key`, default `cell_type`;
4. `annotation_status` (`unknown`, `provisional`, or `curated`), default `unknown`;
5. `reference`, default `rest`, with an explicit single label or comma-separated exact label list allowed;
6. `method` (`t-test_overestim_var`, `t-test`, or `wilcoxon`);
7. advanced `report_p_adjusted`, default 0.05, used for flags/summary only;
8. `max_output_rows` guard.

Outputs, in order:

```text
table, summary, code
```

The table contains the complete canonical columns
`group,reference,regulator,statistic,mean_change,p_value,p_adjusted,direction,significant_adjusted,rank_within_group`
for every tested group-regulator pair. It is deterministically ordered by observed annotation order and then
`p_adjusted`, `p_value`, descending statistic, and regulator ID. There is no `top_n` table truncation and no raw
p-value filter; bounded leading hits belong only in `summary`. `rank_within_group` uses competition/minimum rank on
the scientific tuple `(p_adjusted, p_value, -statistic)`; regulator ID stabilizes display order but does not invent a
different scientific rank. `significant_adjusted` is inclusive (`p_adjusted <= report_p_adjusted`).

`dc.pp.get_obsm` and `dc.tl.rankby_group` are hidden implementation details. The node independently verifies that
the activity artifact's observation IDs align exactly to `adata`; callers never supply or learn an `obsm` key.

## Invariants and reporting

- Annotation labels are nonmissing, nonblank scalar values after collision-safe display normalization; globally
  unused categorical levels are excluded and disclosed.
- Artifact score and adjusted-p-value matrices are defensively owned, finite, cell-aligned, and contain unique
  nonblank regulator names.
- Tested and reference groups are nonempty and meet method-specific minimum sizes; a group equal to its explicit
  reference is not tested. For a reference list, every member label is excluded from canonical tested groups even
  though decoupler 2.2 emits overlapping backend rows for those labels. Comma-list entries are stripped, nonblank,
  distinct, all known, and collision-safe before their union is built.
- Backend output must equal the independently expected 2.2 backend grid with no duplicates, missing rows, or unknown
  rows; only the scientifically valid non-reference group grid is released. Canonical statistics and mean changes
  are finite, p-values and adjusted p-values lie in `[0,1]`, and BH is recomputed for verification. Both t-test
  modes are two-sided Welch tests; `wilcoxon` is specifically SciPy's asymptotic `ranksums` without tie correction.
  Direction follows the signed statistic, while mean-change/statistic sign disagreement is explicitly counted.
- `summary` uses the standard strict report schema, explicitly labels the result `Cluster/annotation-associated TF
  activity evidence`, reports complete-family accounting and bounded key results, carries resource provenance from
  the artifact, and emits dynamic versions without NaN/Infinity.
- `code` defines an equivalent function returning `(table, summary_dict)` and reproduces validation, complete output,
  deterministic sorting, and report construction.

## Scientific boundary

This node must not accept or advertise a Condition contrast. Cells are the descriptive units and are not biological
replicates. `annotation_status=provisional` produces a visible warning; even curated labels do not turn the cell-level
test into Sample-level inference. Technical batch adjustment, repeated measures, and Condition effects are out of
scope.

## Migration and tests

Saved workflows connect `activities`, rename `groupby` to `annotation_key`, add annotation status, remove `top_n` and
`max_pvalue`, and add `summary`/`code`. A migration notice explains that the new table is complete and that
`report_p_adjusted` does not filter it.

Tests must cover all three official methods, exact 2.2 result fields, rest/single/list references, positive and
negative ranks, complete BH families, ties and ordering, missing/colliding labels, unused categories, too-small
groups, artifact/AnnData ID mismatch, malformed backend rows, no input mutation, row guards, strict JSON, dynamic
versions, exploratory wording, and exact runtime/generated equivalence.

## Cohesion assessment

The module owns one decision: how a validated cell-by-regulator activity matrix is characterized across annotation
groups. Backend adaptation, role-safe labels, full-family validation, sorting, and reporting remain local; scoring
and Condition inference remain on their own seams.
