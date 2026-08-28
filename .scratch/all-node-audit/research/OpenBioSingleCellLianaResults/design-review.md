# OpenBioSingleCellLianaResults — design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. The LIANA analysis returns its
> typed primary result directly; no hidden-state extractor remains.

## Decision

**Merge into `OpenBioSingleCellLianaCommunication` and retire as a non-executing compatibility shim.** The producing
analysis returns its primary typed result directly; retaining a second hidden-key extractor creates coupling without
scientific leverage. Do not repurpose this stable ID as a generic AnnData `uns` reader.

## Compatibility-shim interface

Keep the old node ID so saved workflows resolve. Preserve its old `adata,result_key -> table` schema only long enough
to emit a fail-closed migration error. Execution must not read `adata.uns`, guess a method from columns, mutate the
input, or fabricate provenance. The error instructs users to:

1. rerun the enhanced Sample-resolved `OpenBioSingleCellLianaCommunication`;
2. connect its direct `result` output to downstream consumers;
3. remove this node and the `result_key` setting;
4. use `OpenBioSingleCellLianaDotPlot` with the typed result.

As an explicitly retired shim, it has no active scientific result and is exempt from fake `summary`/`code` outputs.

## Why not retain it as an adapter

A valid LIANA table is not identified by being a pandas DataFrame in `uns`. Correct interpretation requires the
method, native score directions, Sample/Condition/identity roles, expression state, resource fingerprint, LIANA
version, and complete-family accounting. Reconstructing those facts from a key or column names is impossible. A
typed direct result carries them at the actual producing seam.

## Migration and tests

The exact legacy-to-current matrix is:

| Legacy contract | Current disposition | Safe translation |
| --- | --- | --- |
| input 0 `adata` | no active input after node removal | Trace it only to identify the legacy Communication producer; never read `uns`. |
| input 1 `result_key` | removed | Delete; no key identifies an audited typed result. |
| output 0 generic `TABLE_RESULT` | no Results output after node removal | A directly connected new DotPlot is rewired to the freshly rerun Communication `OPENBIO_LIANA_RESULT`. Any arbitrary table consumer is type-incompatible and blocks automatic migration. |
| executing extractor | registered non-executing shim until migration | Before workflow migration, loading remains inspectable but every execution fails with the exact actionable message. |

Packaged workflows remove the node only as part of the same atomic rebuild of its Communication producer. Old cached
hidden tables are not accepted as audited evidence because they were pooled and lack resource/Sample provenance. A
standalone Results node, a missing producer, fan-in from an unknown AnnData producer, mixed schemas, or a generic
table consumer is an explicit migration rejection. No automatic table conversion is performed.

Tests must prove that the shim fails before reading any AnnData slot, the message names the replacement and direct
port, registration marks deprecation, packaged workflows contain no active instance, and no output/report/code is
fabricated.

## Cohesion assessment

The deletion test shows that removing this node does not spread complexity: the producer already owns and returns
the typed result. Retirement improves locality by keeping method semantics and provenance where the result is created.

## Registry disposition

Retain the stable ID with `is_deprecated=True`, a “Retired” or “Migration Required” display name, and a description
that names the direct typed result replacement. Registration tests must assert those signals and fail-before-read
execution. This change is intentionally catalogue-only and does not reintroduce a generic `uns` adapter or
fabricate `summary`/`code` outputs.
