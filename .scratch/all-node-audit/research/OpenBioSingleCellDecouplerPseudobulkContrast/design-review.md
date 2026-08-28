# Legacy Decoupler pseudobulk contrast: module design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. The supported surface contains
> only current Sample-level contrast nodes; no decoupler 1.x compatibility entry point is retained.

## Current module

`OpenBioSingleCellDecouplerPseudobulkContrast` currently accepts cell-level generic AnnData plus Sample, population and Condition keys; a reference/comparison; `t-test` or `wilcoxon`; count source; pseudobulk QC thresholds; and normalization target. In one node it:

1. aggregates with historical `decoupler.get_pseudobulk`;
2. stores counts;
3. normalize-total/log1p transforms the profiles;
4. calls historical `decoupler.get_contrast` across populations;
5. calls historical `format_contrast_results`;
6. returns one table.

Current Decoupler 2.2.0 contains none of the three called top-level analysis functions. Decoupler is not declared/installed in the audited local environment, so the node also has no executable supported dependency contract.

## Correctness findings

### P0: the implementation targets a retired API family

Decoupler 2.0 rewrote its functions, removed/deprecated this contrast path, and explicitly directs users to PyDESeq2. Keeping a private emulation under the Decoupler name would misrepresent current official support and burden the codebase with obsolete semantics.

### P0: Technical batch cannot be modeled

The old test consumes only population and Condition. It cannot include Technical batch or other Sample-level covariates, detect full-rank design/contrast estimability, or express an adjusted Condition effect. This violates the repository's Sample/Technical-batch separation for general formal Condition inference.

### Aggregation, transformation, population iteration and testing are over-coupled

One invocation makes several decisions with different reuse/change rates. It repeats aggregation for every method choice, hides every population contrast behind one table, and prevents the same typed profiles from feeding EdgeR and PyDESeq2. It also exposes normalization target as if arbitrary library scaling were part of one stable Decoupler inferential API.

### The statistical/report contract is weak

The historical method uses Sample profiles rather than cells, but it performs a Scanpy rank-genes test on normalized/log profiles rather than a count GLM with library-size offsets and dispersion modeling. Groups with fewer than two profiles are skipped to stderr. The node does not enforce Curated population status, one Condition/Technical batch per Sample, exact replicate disclosure, count-state integrity, covariate adjustment, test-universe semantics, report references/versions, or equivalent source.

### Saved-workflow compatibility cannot validate the missing scientific choices

An old node stores keys and parameters, not proof that Sample is biological replication or that its chosen layer is raw counts. Automatically replacing it with a formal engine would require inventing population selection, design roles and engine choice and could change the analysis silently.

## Decision: retire the analysis; use a non-executing shim only for bounded migration

Do not enhance or reimplement the old analysis. Remove it from the normal node catalogue and replace its registered legacy ID, for at most one documented transition release, with a non-executing schema shim whose only purpose is to let saved workflows load and show an actionable error. Remove the shim in the next major release after migration coverage is complete.

The shim must not import Decoupler, run t-tests/Wilcoxon tests, return stale cached-looking tables, or silently select a replacement. Its message directs the analyst to:

1. run `OpenBioSingleCellPseudobulk` on a declared raw-count source with independent Sample, population, Condition and optional Technical-batch roles;
2. choose exactly one population;
3. choose either EdgeR QL or PyDESeq2 and explicit reference/comparison Conditions;
4. review profile QC, replicate counts, design rank/confounding and expression filtering.

If project policy prefers immediate deletion over a one-release load shim, deletion is scientifically safe; the shim is solely a recoverability aid and does not preserve method compatibility.

## Why this should not be merged into another active node

No valid scientific capability is lost:

- sum aggregation moves to the deeper typed `OpenBioSingleCellPseudobulk` seam;
- count-model Condition inference belongs to the existing EdgeR and PyDESeq2 nodes;
- Technical-batch/design validation is shared privately between those engines;
- a transformed-profile t-test/Wilcoxon has no current concrete workflow consumer or official Decoupler API and therefore fails the deletion test.

Do not add a “method = t-test/wilcoxon/EdgeR/DESeq2” mega-node. That would conflate incompatible inputs, normalization, models, result fields and citations and would make summary/code parity harder to verify.

## Legacy shim interface and behavior

Preserve only enough historical input/output schema metadata for workflow deserialization. On execution, raise one structured migration exception containing:

- legacy node ID and detected saved parameters;
- the reason: Decoupler `get_contrast`/formatter retired and not a supported formal design;
- exact replacement node IDs;
- unresolved required choices: raw count source, independent Sample key, population, Condition, optional Technical batch, reference, comparison, and engine;
- a statement that the old normalized/logged profiles cannot feed a count engine;
- relevant Decoupler changelog URL.

Do not expose active `summary` or `code` endpoints on the shim. The repository requirement applies to executed analysis/processing nodes; pretending a failed migration is an analysis result would be worse than an explicit schema boundary. If the frontend needs structured details, use a migration-specific diagnostic channel owned by migration infrastructure.

## Migration policy

Static workflow migration may safely transform only UI/category metadata and insert replacement placeholders. It must not construct an executable formal workflow unless all semantic roles come from a validated typed upstream artifact and the analyst has already selected the engine/population/contrast.

Explicitly reject automatic mappings in these cases:

- old `reference` is empty/`rest` or comparison/reference is absent/equal;
- multiple populations were implicitly analyzed;
- count layer provenance is missing or matrix state is normalized/logged;
- one Sample maps to multiple Condition or Technical-batch values;
- Technical batch is absent but known confounding cannot be ruled out;
- a legacy test method has no like-for-like replacement;
- fewer than two independent Samples per requested Condition/population remain after profile QC;
- Curated population status is unavailable for a formal claim.

No migration should copy old `target_sum`, old rank-test p-values, or old formatted columns into the new engine's parameters/result contract.

## Documentation/report disposition

The deprecation notice should cite:

- current Decoupler 2.x changelog/source and the retired 1.9.2 implementation;
- the current pseudobulk/DE official tutorial;
- EdgeR or PyDESeq2 sources only in the respective replacement's own report;
- Sample-level pseudobulk practice evidence.

An existing previously executed legacy table remains historical output and should be labeled with the recorded old software/method if it is displayed. Never retroactively assign current Decoupler/PyDESeq2 citations or versions to it.

## Cohesion and coupling assessment

The current node is shallow: callers must coordinate count source, aggregation, profile QC, transformation, population grouping, test selection and formatting in one call, while its only external dependency has removed the method. Deleting the executor increases cohesion of all surviving modules. Aggregation becomes one reusable state transition; each count engine becomes one model/contrast; migration handling becomes one non-scientific compatibility boundary.

The one-release shim is intentionally shallow because it contains no domain logic. Its value is recoverability of saved graphs, not reuse or scientific leverage. It should have a scheduled removal criterion rather than accumulate new features.

## Verification plan for the future implementation batch

- Current Decoupler 2.x environment never attempts missing legacy imports.
- Legacy workflow fixture deserializes through the shim and receives deterministic actionable migration guidance.
- Shim never executes a test and exposes no scientific `summary`/`code` outputs.
- Node is absent from new-workflow discovery/catalogue while the legacy ID remains resolvable only for the transition window.
- Migration tests refuse `rest`, ambiguous populations, transformed/unknown count sources, missing semantic roles and engine choice.
- No automatic mapping treats Technical batch as Sample or chooses EdgeR/PyDESeq2.
- Replacement pipeline tests prove typed aggregation plus one-population engine wiring, replicate/design disclosure and required `summary`/`code` parity.
- Final-major-removal test and release note remove the shim after the announced boundary.

## Registry disposition

For the compatibility window, set `is_deprecated=True`, label the display name “Migration Required” or “Retired”,
and expose an actionable non-executing description. A registration test must cover those public signals for the
exact shim ID. This is a catalogue correction only: it neither selects an inferential engine nor creates a fake
scientific output contract.
