# OpenBioSingleCellDGIdbAnnotation — design review

## Decision

**Retire the legacy AnnData-annotation behavior and repurpose this registered ID as the atomic DGIdb resource
loader/validator.** Display it as “Load DGIdb Resource”. This is a deliberately breaking correction: a hidden remote
annotation side effect cannot serve as a reproducible resource seam.

The stable ID is kept so workflows remain discoverable, but `adata -> adata` links cannot be auto-migrated. If ID
semantics are required to remain literal, retain a non-executing shim and introduce a new resource ID; under no
circumstance keep the old download behavior as an automatic fallback.

## Atomic interface

Inputs, in order:

1. local DGIdb `.tsv`/`.csv` snapshot under the input boundary;
2. strict `resource_metadata_json` containing nonblank `name`, `version`, `release_date`, `download_url`, `organism`,
   `gene_identifier_namespace`, `scope`, `license`, `source_license_review`, and `citation`;
3. advanced exact `drug_column`, default `drug_claim_name`;
4. advanced exact `gene_column`, default `gene_claim_name`;
5. optional source/evidence column names retained as provenance;
6. advanced maximum file bytes and row count.

Outputs:

```text
resource, summary, code
```

`resource` is an immutable typed artifact containing a canonical long DataFrame, strict metadata, streamed file
SHA-256, canonical-content fingerprint, deterministic drug order, duplicate/source accounting, and producer schema.
Generated code returns the portable `(DataFrame, summary_dict)` equivalent instead of the process-local wrapper.

No AnnData, online URL fetch, cache path, “latest” version, annotation column, expression source, drug filter, or
statistical parameter is exposed. Parsing and fingerprinting form one deep local-substitutable module shared by all
drug evidence nodes.

## Validation and disclosure

Reject blank/non-scalar identifiers, delimiter/field-count corruption, malformed metadata, ambiguous namespace,
empty data, duplicate column names, and rows beyond guards. License and source-review text is required for disclosure
but is a caller declaration: placeholder or unreviewed status emits a warning instead of pretending the program can
verify a legal review. Trim surrounding identifier
whitespace only; reject normalization collisions; collapse exact drug-gene duplicates and preserve contributing
source/evidence values in deterministic arrays. Human/HGNC is required for the current downstream implementation.

`summary` is strict JSON and reports all resource/accounting/license/version/hash facts plus citations, warnings, and
the DGIdb research-only/non-efficacy limitation. `code` repeats byte/content pinning and returns the same summary.

## Sample, Technical batch, and scientific boundary

This adapter performs no expression analysis or inference, so Sample and Technical batch do not enter its interface.
Downstream nodes remain responsible for Sample-level Condition semantics. A DGIdb edge is an aggregated interaction
claim; it does not establish drug efficacy, direction, dose, safety, or relevance to a population.

## Migration and tests

Legacy annotation links must stop with an actionable migration note: load a reviewed snapshot, wire the resource to
one evidence node, and use a separate generic merge/annotation operation only if a var annotation is truly needed.
Do not download or invent metadata during migration.

Tests: Pertpy fresh-object regression showing `.dictionary` is absent; no-network enforcement; quoted fields and
strict row widths; missing/blank/collision identifiers; metadata fields; reviewed/unreviewed/placeholder license
declaration warnings; exact duplicate aggregation;
source/evidence preservation; raw and canonical SHA changes including same-size/mtime changes; file/row guards;
deterministic order; artifact tampering; strict JSON; code compilation; and runtime/generated table+summary equality.

## Cohesion assessment

The new module has one reason to change: the accepted DGIdb snapshot contract. This makes the resource a real seam
with two consumers (runtime evidence and generated-code reproduction), while drug scoring and enrichment remain
independent modules.
