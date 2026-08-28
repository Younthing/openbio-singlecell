# Normalization, feature modeling, and cNMF

Type: task
Status: resolved

Blocked by: 03

## Initial audit focus

Remove implicit counts/raw side effects, distinguish requested/effective sentinel values, enforce flavor-specific expression semantics, prevent silent layer overwrite, and isolate cNMF run artifacts.

## Comments

- 2026-08-28: Initial inventory and code audit complete; per-node pre-change documents pending.
- 2026-08-28: Landed official-usage and design-review records for all eight preprocessing/cNMF nodes before implementation. Normalization and feature-modeling nodes now enforce representation-specific state, protected layer destinations, copy-on-write behavior, flavor-specific scientific validation, strict reports, and equivalent source.
- 2026-08-28: Split cNMF into Rank Survey and Consensus Programs through the concrete process-local `OPENBIO_CNMF_RUN` contract. Removed hard RFC/argmax labels and shared filesystem/worker controls; added complete-restart checks, rank metrics, explicit density/neighborhood decisions, conservative 2 GiB resource preflight, canonical aligned continuous outputs, and safe v1/legacy canvas migration including subgraphs, bypass, exposed widgets, and reroutes.

## Answer

The six normalization and feature-modeling nodes remain audited and complete. Expression state is resolved per representation, so normalized/logged/scaled/residual matrices cannot silently re-enter count-dependent methods; generated code reproduces the same state and warning contracts. HVG selection reports algorithmic, forced, final, grouped-intersection, and Raw-subset outcomes, while Scale and Pearson residuals preserve source states and preflight densification.

The two staged cNMF nodes are temporarily reopened for migration to the method authors' exact `cnmf==1.7.1` file-backed API. Their public schemas, typed `OPENBIO_CNMF_RUN` seam, and reviewed workflow migration remain unchanged, but the earlier OmicVerse implementation and its `300 passed, 1 skipped` closure claim are superseded. Batch04 will return to resolved only after managed-temporary lifecycle, official-file completeness/tamper checks, generated-code parity, exact-warning isolation, and a real 1.7.1 Survey-to-Consensus smoke are independently verified.

- 2026-08-28: Final release dependency resolution invalidated the earlier optional-backend assumption. The published
  `omicverse==2.3.1` distribution requires `anndata<0.12` and `pandas<3`, while this audited package requires
  AnnData 0.13 and supports Pandas 3, so no standard `.[cnmf]` extra can install the tested OmicVerse adapter.
  Both cNMF nodes were returned to provisional status rather than advertising an un-installable backend. The
  method-author `cnmf==1.7.1` wheel resolves with the repository's exact core stack and explicitly includes the
  NumPy-2 clustering compatibility update; an exact file-backed adapter, real smoke, regenerated standalone code,
  and dependency/notice changes are now required before Batch04 can be closed again.

- 2026-08-28: Reclosed after replacing OmicVerse with the method-author `cnmf==1.7.1` file-backed API. Both nodes
  own a managed temporary workspace, validate the complete official file set and exact axes before reuse, preserve
  an explicit Raw/X/layer expert choice without history or Raw-binding gates, and emit strict summaries plus
  standalone equivalent source. A real 1.7.1 Survey-to-Consensus smoke, focused 37/37 tests under default and
  `-W error`, registration, Ruff, and an independent P0/P1/P2=0 review passed. The release extra, manifest, README,
  notices, and typed `OPENBIO_CNMF_RUN` contract now name the installable method-author distribution.
