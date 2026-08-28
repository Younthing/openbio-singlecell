# OpenBioSingleCellRunPySCENIC — design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. The supported seam is an
> external current pySCENIC run followed by strict `OpenBioSingleCellImportPySCENICResults` import.

## Decision

**Retire as a non-executing compatibility shim.** The supported workflow is an explicitly versioned external
pySCENIC 0.12.1 container (or equivalent isolated environment) followed by
`OpenBioSingleCellImportPySCENICResults`. Do not silently patch NumPy aliases, pin the entire plugin to an obsolete
NumPy line, or launch a privileged container engine from this node.

This is a compatibility-driven retirement, not a rejection of SCENIC. The analysis remains supported at the
external artifact/import seam where its resource and environment provenance can be validated.

## Compatibility-shim interface

Keep the existing node ID so saved workflows resolve. Preserve enough of the old input schema to produce a precise
migration error, but `execute` always fails before reading expression, downloading resources, creating temporary
files, or starting a subprocess. The message must state:

1. pySCENIC 0.12.1 is incompatible with the plugin's Python 3.12/NumPy 2 runtime;
2. the former node deleted critical external artifacts and cannot meet the audit contract;
3. run the official three-stage container workflow from the `code` template/documentation;
4. retain the required manifest and outputs;
5. import them through `OpenBioSingleCellImportPySCENICResults`.

As an explicitly retired shim it has no scientific primary result and is exempt from fabricated `summary`/`code`
ports. Registration and workflow migration must mark it deprecated/non-executing. Do not reuse the stable ID for a
mere expression exporter: the display name promises a completed pySCENIC run.

## Supported external seam

The importer accepts one versioned `openbio-singlecell/pyscenic-external-run/v1` bundle manifest. Relative paths in
that manifest identify the exported expression input, adjacency table, motif-enrichment CSV, AUCell CSV, TF list,
every ranking database, and motif annotation file; every byte stream is SHA-256 verified before a result is parsed.
The manifest schema fixes the scientific state and resources without making Docker, a scheduler, or a remote cluster
part of the plugin interface. The external workflow may vary internally as long as its outputs and manifest satisfy
that interface; this creates a deep, testable import module instead of embedding environment orchestration in every
caller.

The recommended documented export source is the Raw snapshot: post-QC full-gene counts, unique observation/gene IDs,
explicit organism and gene namespace. No default HVG restriction or log-normalized layer is silently substituted.
An expert finite numeric state is importable only when its exact state label and value fingerprint are declared and
warned. Resource bundle compatibility and immutable hashes are checked before the external run and again at import.

## Migration and tests

Packaged workflows must remove this node and replace it with the external-run/import boundary. Old cached outputs
are invalid because the temporary adjacency/regulon/AUCell files and their hashes were not retained. Migration must
not infer a manifest from filenames or trust the old `ScenicNetwork.provenance` as a completed-run record.

Tests must prove that every old input combination fails before optional imports, network access, subprocess calls,
temporary-file creation, or input mutation; the error names the importer and required manifest; schema registration
marks deprecation; packaged workflows contain no executing instance; and no fake report is emitted.

## Cohesion assessment

Deleting in-process execution removes an interface that exposed host orchestration, three scientific stages,
resources, storage, and import as one operation. The replacement seam concentrates validation and interpretation in
the importer while leaving execution to the environment designed to run pySCENIC.

## Registry disposition

Keep the stable deserialization ID with `is_deprecated=True`, a visibly retired display name, and an actionable
non-executing description. Tests must verify these public catalogue signals in addition to proving that execution
fails before reading data or touching the environment. Deprecation does not change the external artifact/import
architecture or add fake result ports.
