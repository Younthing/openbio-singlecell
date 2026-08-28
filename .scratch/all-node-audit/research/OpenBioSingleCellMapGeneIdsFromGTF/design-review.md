# Map Gene IDs from GTF: module design review

## Current module

The node streams gene records from a plain/gzipped GTF, maps current `var_names` to `gene_name`, optionally drops unmatched features, optionally makes resulting names unique, and stores the old index in a user-named column.

## Correctness and interface problems

1. Default Scanpy 10x input uses gene symbols as `var_names` and stores stable IDs in `.var["gene_ids"]`, so the current default assumption commonly produces zero matches.
2. `_stable_gene_id(value.split(".", 1)[0])` removes the GENCODE `_PAR_Y` suffix along with the version.
3. `setdefault` silently accepts the first name when normalized IDs conflict.
4. Replacing stable IDs with possibly duplicated/changing symbols, then `make_unique()`, manufactures technical names rather than preserving gene identity.
5. `drop_unmapped` mixes annotation with feature filtering and changes the expression axis.
6. If `.raw` exists, only current `.var` is changed, leaving incompatible feature namespaces.
7. File provenance and annotation provider/release/build are not reported.

## Decision: keep, redefine, and enhance

Retain the GTF seam, but redefine it as **canonicalize and annotate gene identity from GTF**. The atomic result is stable feature identity plus GTF-derived annotation; it must not filter features. Run it after loading and before QC/Raw snapshot creation.

The module should keep stable gene IDs as the principal feature key and write GTF `gene_name` to a variable annotation such as `.var["gene_symbols"]`. Repeated symbols are valid annotations. Add `.var["gtf_mapped"]` rather than dropping unmatched features. If source IDs live in a var column while the existing index holds display names, preserve the prior index in a fixed provenance annotation before promoting stable IDs; do not expose arbitrary rename mechanics.

## Target interface

- Inputs: `adata`, `gtf_path`.
- Visible identity choice: `id_source = var_names | var_column`, default `var_column` for 10x interoperability.
- Advanced: `id_column="gene_ids"` when `id_source=var_column`; `gene_name_column="gene_symbols"`.
- Removed: `original_id_column`, `drop_unmapped`, and `make_unique`.
- Fixed hidden outputs: stable IDs remain/become the principal index; `gtf_mapped` records mapping evidence; no feature is removed.
- Outputs: annotated `adata`, `summary`, `code`.

The source and destination columns are exposed because they identify the external schema. Exact-then-version-normalized matching, `_PAR_Y` preservation, ambiguity detection, streaming, and provenance parsing are hidden implementation policy: callers should not tune correctness rules.

## Matching and validation contract

- Reject an empty dataset, blank/missing source column, blank or duplicate stable IDs, and any existing `.raw`; this operation must precede Raw snapshot creation.
- Parse plain or gzip GTF incrementally. Prefer `feature == "gene"`. If a community GTF lacks gene rows, optionally fall back to other records carrying both `gene_id` and `gene_name`, and disclose that fallback.
- Parse/report header metadata when available and fingerprint the resolved file by portable path, size, and modification time. Missing release/build information is a report limitation.
- Try exact IDs first. Only exact failures may use terminal numeric-version removal that preserves `_PAR_Y`.
- Allow repeated identical GTF records. Detect one normalized stable ID associated with different names as ambiguous; never select the first silently.
- Treat version-normalized matches as valid but warn that matrix and annotation releases may differ.
- Permit multiple stable IDs to share a symbol; report the duplication without manufacturing suffixes.
- If `gene_name_column` already exists, allow equal values and fill missing values. Reject non-missing conflicts by default rather than overwrite them.
- Require at least one successful mapping. Leave unmatched features and expression matrices unchanged with `gtf_mapped=False`.

## Report contract

`summary` includes input/output feature counts; identifier source; exact, version-normalized, unmapped, and ambiguous matches; gene-symbol duplicate counts; existing-name equal/fill/conflict counts; GTF parsed/skipped/malformed/conflicting record counts; provider/release/build/fingerprint; warnings and limitations; GENCODE/Ensembl/AnnData references; and dynamic Python/openbio-singlecell/AnnData/pandas/NumPy versions. A report-ready result states how many of G stable IDs received GTF gene-symbol annotation without implying that unmapped features were discarded.

`code` must be self-contained, include the GTF parser and `_PAR_Y`-safe version helper, validate the same invariants, avoid ComfyUI-private path resolution, annotate without filtering, and return a new AnnData equivalent to the primary result.

## Cohesion and coupling assessment

This redesign creates a deep annotation module: a small identity interface hides file parsing, exact/version-aware reconciliation, ambiguity checks, axis preservation, provenance, and reporting. Feature filtering, arbitrary regex edits, and package-specific symbol-index adaptation remain outside the seam. That separation keeps stable identity local and prevents symbol formatting from coupling every downstream analysis.

## Verification plan

- Plain and gzipped GTF, comments/header metadata, primary gene rows, and the disclosed no-gene-row fallback.
- `var_names` and `.var["gene_ids"]` sources; exact, version-normalized, unmapped, and zero-match datasets.
- `ENSG….7_PAR_Y` normalizes only its numeric version and retains `_PAR_Y`.
- Repeated same-ID/same-name records are tolerated; same normalized ID/different names are rejected; different stable IDs sharing one symbol are retained.
- Existing gene-symbol annotations cover equal, missing-fill, and conflicting values; duplicate/blank stable IDs fail.
- Existing Raw snapshot, invalid paths/extensions, empty/malformed GTF, and absent release/build metadata follow the stated failure/limitation behavior.
- Expression matrices and feature count remain unchanged, summary is strict JSON, and self-contained generated code compiles and reproduces the annotations and stable identity.

## Open expert-boundary decision (2026-08-28)

Revise the earlier pre-Raw/coverage/uniqueness contract. Allow Raw, zero-match, empty-axis, duplicate stable IDs,
version-normalization collisions, and pre-existing symbol conflicts; preserve Raw and existing conflicting values and
record each condition in `summary`. Keep path/parse ambiguity, missing declared ID source, blank/missing IDs, reserved
column overwrite, and invalid destination-column relationships as hard errors. This keeps annotation atomic and gives
experts the complete evidence without imposing a preferred annotation release. Do not add Raw-binding tests.
