# OpenBioSingleCellDrugHypergeometric — design review

## Decision

**Enhance and retain as an atomic DGIdb explicit-set ORA evidence module.** Delete internal Scanpy marker ranking and
resource construction. Share the hardened ORA implementation with generic/Marker ORA; retain this public identity
for drug-specific resource validation and interpretation warnings.

## Atomic interface

Inputs, in order:

1. canonical selected-evidence `table` (generic producer/approval metadata is declarative);
2. paired exact tested-gene `universe`;
3. typed DGIdb `resource`;
4. required `comparison`/`group` selector if needed;
5. post-universe `min_targets`;
6. `min_overlap` and `max_p_adjusted` for flags only;
7. output-row guard.

Outputs:

```text
table, summary, code
```

Remove AnnData, expression source, `groupby`, marker threshold/rank method, resource download, identifier translation,
and result filtering. Fix one-sided greater Fisher/hypergeometric semantics, Haldane 0.5, exact background, and BH
over the complete retained drug family.

## Evidence contract

Hard-validate selected/universe roles and analysis/ranking/universe/current-content fingerprints, selected genes
inside universe, DGIdb producer/content/version/license/hash, human HGNC namespace, and current-content tampering.
Generic producer/approval/direct-upstream/truncation and inference-scope metadata are user declarations; disclose
missing, unfamiliar, contradictory, or cell-level values with cautious warnings rather than rejecting the arithmetic.
Return one row per retained
drug with exact set/selection/universe sizes, overlap JSON, `a,b,c,d`, log odds ratio, raw/BH p, rank, flags, and scope.

`summary` includes report-ready method/results, upstream design, complete drug family/accounting, top alternatives
and ties, resource/evidence sources, warnings/limitations and clinical disclaimer, references, and dynamic versions.
`code` consumes portable tables/metadata and returns `(DataFrame, summary_dict)` with identical core and report.

## Sample and Technical batch semantics

The node does no expression model. Sample-level Condition inference and Technical-batch handling/audit remain required
for a formal biological claim, but generic metadata cannot prove them. Unsuitable or absent declarations make the
result exploratory and trigger a warning rather than a computation error. Cluster marker evidence is explicitly
non-confirmatory. No result may be described as therapeutic response or recommendation.

## Migration and tests

Legacy links cannot be auto-migrated because the node synthesized its own rank and never emitted the exact universe.
Migration requires direct selected+universe evidence and a reviewed resource artifact. Do not infer these from the
old AnnData.

Tests reuse ORA hand contingencies plus fresh-DGIdb regression, structural-pair/tamper/cross-universe failures,
custom/missing/cell-level provenance disclosure, DGIdb namespace/version/license/hash, duplicate targets, zero overlap, complete
family/BH, ties/flags, no network/no Scanpy ranking, guards, immutability, strict JSON, compiled generated code, and
runtime/generated equivalence.

## Cohesion assessment

One module now answers one question: which pinned DGIdb target sets are overrepresented in one explicit evidence set?
Ranking, filtering, resource lifecycle, and clinical interpretation stay outside, while shared statistical complexity
is localized.
