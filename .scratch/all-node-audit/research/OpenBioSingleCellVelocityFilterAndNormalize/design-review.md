# Velocity Filter and Normalize: module design review

## Decision

**Retain the stable node ID but reshape it into an atomic “prepare velocity abundances” module.** The current public name can remain for migration, while the display name should state the actual 0.3.4 behavior. Remove the invalid HVG arguments. Do not merge it with moments or model fitting.

The deletion test justifies a typed `OPENBIO_VELOCITY_STATE`: without it, every downstream node must rediscover whether generic AnnData layers are raw/prepared, aligned, filtered, normalized, and bound to the same cells/features.

## Target interface and artifact

Inputs: `adata`; visible `spliced_layer="spliced"`; visible `unspliced_layer="unspliced"`; visible `min_shared_counts=20`; advanced `min_shared_cells=0`; advanced `normalization_target="median_library"`; advanced `overwrite_existing=False`. Remove `n_top_genes` and free-form `subset_highly_variable`.

Outputs: `velocity_state` (`OPENBIO_VELOCITY_STATE`), `summary`, `code`.

The process-local artifact owns a defensive AnnData copy and immutable schema/provenance containing stage `prepared`, exact obs/var fingerprints, source layer identities and count-state evidence, thresholds, normalization policy, retained/removed gene identities, canonical `spliced`/`unspliced` layer fingerprints, and scVelo version. It exposes defensive copies only; later stages return a new artifact.

The implementation validates raw integer-like layers, copies them to canonical names, subsets genes by explicit shared thresholds, and normalizes the canonical velocity layers with `scv.pp.normalize_per_cell(..., layers=["spliced", "unspliced"], enforce=True)`. Because the audited 0.3.4 implementation always prepends `X`, the wrapper must snapshot and restore private-copy `X` plus any pre-existing `obs["n_counts"]`, then verify exact `X` restoration and that only canonical velocity layers persistently changed. It must verify feature alignment, finite values and nonzero libraries. Do not call the misleading all-in-one recipe and do not mutate/normalize expression `.X` as a public side effect.

`summary` includes layer/count diagnostics, retained/removed gene counts and bounded identities, thresholds, normalization totals before/after, artifact fingerprints, references, limitations and scVelo/AnnData/NumPy/SciPy/OpenBio versions. `code` returns a portable `(prepared_adata, summary_dict)` with identical canonical layers and metadata; that is the external equivalent of the typed artifact.

## Migration and tests

Legacy `n_top_genes`/`subset_highly_variable` values cannot be honored by 0.3.4 and are removed with a migration warning. If users require feature ranking, it must become a separately researched node. Existing custom layer names map directly.

Tests cover dense/CSR/CSC counts, aliases, missing/misaligned/negative/nonfinite/noninteger/already-normalized layers, gene filtering exactness, zero libraries, feature/obs duplicates, artifact defensive ownership and stage/fingerprints, input immutability, strict JSON, source compilation, fake/real 0.3.4 smoke, and runtime/code parity.
