# OpenBioSingleCellDrugGSEA — design review

## Decision

**Enhance and retain as an atomic DGIdb ranked-enrichment evidence module.** Delete internal marker ranking and
resource construction. Share the complete-rank GSEA core with Ranked GSEA; keep this node's public interface for the
typed DGIdb resource and required clinical disclaimer.

## Atomic interface

Inputs, in order:

1. canonical complete ranked-evidence `table`;
2. paired exact `universe`;
3. typed DGIdb `resource`;
4. required `comparison`/`group` selector if necessary;
5. advanced `gene_column`/`score_column`;
6. post-universe `min_targets`/`max_targets`;
7. `n_permutations` (hard minimum 2; values below 1000 warning-grade);
8. nonnegative `random_seed`; zero remains executable but disables decoupler 2.2 shuffling and invalidates formal
   empirical-null interpretation;
9. output-row guard.

Outputs:

```text
table, summary, code
```

Remove AnnData, expression source, `groupby`, internal Wilcoxon, marker threshold, resource URL/cache, absolute-rank
switch, and output filtering. One invocation processes one complete rank. Fixed backend choices are decoupler 2.2
public `mt.gsea`, exact rank orientation, full drug family, and accurately named adjusted p-values.

## Evidence contract

Hard-validate table/universe roles, complete unique finite nonconstant ranking, all content/ranking/universe
fingerprints, DGIdb artifact fingerprint/release/license/HGNC namespace, and set intersection and size family.
Generic producer/approval/direct-upstream/truncation/scope/direction/Sample/batch fields are unverified
declarations; preserve and caution on missing, unfamiliar, contradictory, or cell-level values without blocking an
otherwise defined GSEA. Return all tested drugs with comparison, NES, adjusted p, rank/significance,
set-size accounting, and deterministic resource order. No clinical direction field may be inferred from NES.

`summary` includes report-ready method/results, upstream Sample/Technical-batch provenance, complete rank/resource
accounting, top positive/negative ties, warnings/limitations and explicit non-efficacy/non-signature-reversal text,
references, and dynamic versions. `code` consumes portable tables/metadata and returns `(DataFrame, summary_dict)`
using the same core.

## Sample and Technical batch semantics

Formal Condition ranks should use Sample-level inference and Technical-batch handling. Generic declarations do not
prove that design: absent or unsuitable values make interpretation explicitly exploratory and emit a warning rather
than a hard error. Cluster-marker ranks remain exploratory. The node itself has neither cells nor Samples and cannot
repair upstream pseudoreplication.

## Migration and tests

No automatic migration is scientifically possible from the old `adata + groupby` input because it did not preserve
an external complete rank/universe or resource release. Stop with an actionable wiring note.

Tests: official decoupler 2.2 output; cautious seed-zero/low-permutation execution; complete versus truncated rank;
duplicate/non-finite/tied
genes; cross-pair/fingerprint tampering; custom/missing/cell-level producer declarations with cautious summaries;
DGIdb release/license/hash/namespace; after-universe set bounds; positive and
negative hand fixtures; complete family/adjusted-p naming; no network/no Scanpy ranking; output guards; immutability;
strict JSON; compiled code; and runtime/generated equality.

## Cohesion assessment

The module now answers one ranked drug-target enrichment question. Resource lifecycle and ranking are clean upstream
seams; shared GSEA implementation creates locality without conflating drug-specific interpretation with generic sets.
