# Annotation round trip

Type: task
Status: resolved

Preserve reviewed annotation provenance when a subtype column is transferred back to its parent AnnData.

## Comments

- Keep the existing Merge Observation Annotations interface.

## Answer

`Merge Observation Annotations` now transfers compatible annotation provenance to the target column, redirects its
stored output-column identity, rejects conflicting target annotation provenance, and rejects explicitly different
OpenBio source provenance while retaining the existing missing-provenance advisory.
