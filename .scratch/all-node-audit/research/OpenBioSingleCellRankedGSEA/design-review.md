# OpenBioSingleCellRankedGSEA — design review

## Decision

**Enhance and retain as an atomic one-ranking GSEA evidence module.** Remove internal marker ranking and generic-table
guessing. Share a deep internal ranked-enrichment/resource validator with Drug GSEA, but retain distinct public nodes
because the drug resource has specialized provenance and clinical limitations.

## Atomic interface

Inputs, in order:

1. canonical complete ranked-evidence `table`;
2. exact tested-gene `universe` artifact paired to that ranking;
3. local `gene_sets_file` and strict `resource_metadata_json`;
4. required `comparison`/`group` selector when the artifact contains several complete rankings;
5. advanced canonical `gene_column` and `score_column`;
6. resource `source_column`, `target_column`;
7. `min_targets` and `max_targets`, both applied after universe intersection;
8. `n_permutations` (hard minimum 2; values below 1000 warning-grade);
9. nonnegative `random_seed`; zero is executable but prominently disclosed as disabling decoupler 2.2 shuffling;
10. advanced output-row guard.

Outputs:

```text
table, summary, code
```

No AnnData, expression source, marker-generation method, significance cutoff, absolute-rank switch, gene-translation
option, or p-value display filter belongs here. One invocation analyzes one complete ranking. `raw=False`,
`empty=False`, and `verbose=False` are hidden.

## Evidence contract

Require the structural artifact role/schema, exactly one row per universe gene for the selected ranking, finite
nonconstant scores, and matching recomputed content/universe/ranking fingerprints. Exact duplicate resource pairs
collapse; identifiers are not normalized. Resource sets surviving intersection and size bounds define the complete
tested family. A generic table's producer identity, `enrichment_approved`, scope/direction, truncation declaration,
inference unit, replicate awareness, and Technical-batch handling are caller/producer assertions, not an authority
boundary. Preserve them verbatim when JSON-safe and emit cautious warnings when absent, unfamiliar, contradictory,
or unsuitable for formal Condition interpretation; do not hard-fail the calculation solely on those declarations.

Canonical table:

```text
comparison, gene_set, normalized_enrichment_score,
p_adjusted, rank, significant, set_size_before_universe,
set_size_in_universe
```

All retained sets are returned; significance is descriptive and never filters rows. The adjusted p-value column is
named accurately. Stable sorting is adjusted p, descending absolute NES, then original resource-set order; complete
ties share rank.

`summary` includes a report-ready method/result narrative, evidence scope and upstream design, set/universe/resource
accounting, bounded top positive/negative results and ties, permutations/seed, warnings/limitations, references, and
dynamic software versions. `code` returns `(DataFrame, summary_dict)` and performs identical validation/statistics.

## Sample and Technical batch semantics

This node performs no replicate-level model. Formal Condition rankings should already use Sample as inference unit
and model Technical batch as nuisance; cluster-marker rankings remain exploratory. Generic producer and inference
metadata are disclosed as unverified declarations. Cell-level, missing, or unfamiliar declarations therefore trigger
an explicit “not validated for formal Condition inference” warning rather than rejection. Rank completeness,
table/universe axis integrity, finite scores, and current-content fingerprints remain hard requirements because they
define the computation itself.

## Migration and tests

The old marker-table-only interface cannot be auto-migrated because no exact paired universe/completeness proof was
required. Migrations add a note requiring a compatible direct ranking+universe pair and resource metadata. Preserve
an old seed 0 only with the explicit warning that decoupler 2.2 does not shuffle that empirical null.

Tests: official 2.2 tuple/output semantics; absence/1.x APIs; seed-0 and low-permutation cautious execution;
hand-ranked positive/negative sets;
permutation/BH output naming; incomplete/duplicate/non-finite/constant ranks; cross-paired/tampered fingerprints;
after-universe size filtering; malformed/changed resource; no-overlap/whole-universe sets; ties/order; custom or
missing producer/inference declarations with cautious summaries; output guard; strict JSON; immutable inputs; code
compilation; and runtime/generated equivalence.

## Cohesion assessment

The module owns exactly one ranked-set enrichment method. Upstream ranking and downstream visualization stay outside;
backend quirks, resource validation, completeness proof, set accounting, and report generation remain local.
