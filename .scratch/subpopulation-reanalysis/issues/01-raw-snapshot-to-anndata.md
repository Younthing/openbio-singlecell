# Raw Snapshot to AnnData

Type: task
Status: resolved

Implement the confirmed Raw materialization node and its public behavior tests.

## Comments

- The user explicitly confirmed that `adata.raw.to_adata()` is the intended restoration primitive.

## Answer

Added `OpenBioSingleCellRawSnapshotToAnnData`. It materializes the selected Raw snapshot through
`adata.raw.to_adata()`, preserves the current observations and OpenBio provenance, removes inherited analysis
slots, reports that Raw state is not inferred, and exposes equivalent code.
