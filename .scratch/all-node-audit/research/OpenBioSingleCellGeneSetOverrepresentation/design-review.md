# OpenBioSingleCellGeneSetOverrepresentation — design review

## Decision

**Enhance and retain as a generic explicit-selected-set ORA evidence module.** Merge its statistical implementation
with the already hardened Marker ORA core through an internal seam, but do not merge the public identities: this node
accepts generic compatible evidence while Marker ORA enforces a specialized Cluster-marker direct-pair contract.

## Atomic interface

Inputs, in order:

1. canonical selected-evidence `table` with structural/fingerprint metadata (producer approval is a declaration);
2. its exact paired tested-gene `universe` artifact;
3. local `gene_sets_file` and strict `resource_metadata_json`;
4. required `comparison`/`group` selector if multiple selected sets are present;
5. advanced resource `source_column`, `target_column`, and post-intersection `min_targets`;
6. `min_overlap` and `max_p_adjusted` used only for flags;
7. advanced output-row guard.

Outputs:

```text
table, summary, code
```

Remove AnnData, `universe_source`, selection column/operator/threshold, gene translation, alternative, arbitrary
background size, and multiple-testing method. Filtering belongs upstream. Fix one-sided greater Fisher semantics,
Haldane correction 0.5, exact background, and BH across the full source family.

## Evidence and table contract

Hard-validate schema/artifact roles, analysis/ranking/universe/current-content fingerprints, unique identifiers,
selected genes inside the universe, and one selected evidence set. Producer identity/approval, direct-upstream,
truncation, scope, direction, inference unit, replicate awareness, and batch handling are generic caller declarations:
preserve them in the summary and warn when absent, unfamiliar, contradictory, or cell-level. They do not establish an
authority boundary and must not independently block a structurally valid ORA.

Return one row per retained resource set with selected/universe/set sizes, overlap JSON, `a,b,c,d`, log odds ratio,
raw p, BH p, deterministic rank, overlap/significance/candidate flags, and scope/comparison. Full rows are retained.
Resource metadata and streamed SHA-256 are part of the result fingerprint.

`summary` is the strict scientific report with method/results, upstream inference and selection, contingency and
resource accounting, bounded alternatives/ties, warnings/limitations, references, and dynamic versions. `code`
uses the same pure core/summary builder and returns `(DataFrame, summary_dict)`.

## Sample and Technical batch semantics

For a Condition contrast, Sample-level inference and a disclosed Technical-batch model/audit remain the recommended
scientific design. The ORA node cannot verify those declarations, does no expression modeling, and cannot adjust
Technical batch itself. Missing or unsuitable declarations produce a prominent “formal Condition interpretation not
validated” warning; Cluster marker declarations remain explicitly exploratory.

## Migration and tests

Legacy workflows cannot be auto-converted: current AnnData genes are not proof of the upstream tested universe and
the node currently performs its own selection. Keep links inspectable and emit a migration-required note pointing
to an explicit selected-set+universe pair. Do not guess selection or resource metadata.

Reuse hand-calculated tests from Marker ORA plus generic producer/scope fixtures: same-size/different-gene universe,
cross-pair and numeric tampering, custom/missing/cell-level producer declarations with cautious summaries, selection outside universe, complete family,
zero overlap, duplicates, resource metadata/hash changes, set pruning, Fisher/Haldane/BH verification, ties/flags,
output guards, immutability, strict JSON, generated-code compilation, and exact runtime/generated equivalence.

## Cohesion assessment

The public module has one reason to change: the definition/reporting of generic explicit-set ORA. Shared statistical
and resource machinery stays in one deep internal module; upstream filtering and biological interpretation remain
separate.
