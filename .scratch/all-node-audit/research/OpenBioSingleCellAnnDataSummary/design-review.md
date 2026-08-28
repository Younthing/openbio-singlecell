# AnnData Summary: module design review

## Current module

The node is a read-only report node whose primary result is already a `SummaryResult`. It currently returns shape, `repr(adata)`, and bounded names for a subset of AnnData mappings.

## Correctness and reporting gaps

- The current payload predates the standard analysis-report schema and lacks method/reference/software-version fields and a report-ready results narrative.
- It omits `X` dtype/storage, index uniqueness, `raw`, `obsp`, `varp`, view/backed state, and detailed source provenance.
- A free-form `repr` is version-dependent and less useful than deterministic structured fields.
- It does not expose equivalent source through a `code` port.
- It can tempt readers to infer biological meaning from structure unless limitations explicitly state that no expression-state or biological validity assessment was performed.

## Decision: enhance, do not merge

Keep a standalone structural report. It is useful at any workflow seam and is orthogonal to QC or biological analysis. Merging it into every loader/transform would duplicate inventory logic and couple those nodes to presentation.

Because the primary output is already the required `summary` artifact, expose it exactly once and add `code`; do not add a duplicate summary port. The report must use the audit-wide strict schema and include AnnData documentation/paper references and runtime versions.

## Target interface

- Input: `adata` only.
- Visible parameters: none; the node makes no scientific choice.
- Advanced parameters: none.
- Hidden fixed behavior: read-only traversal; deterministic ordering; a fixed bounded-name limit (64) with explicit totals/truncation; strict JSON normalization; no full matrix materialization; no biological interpretation.
- Outputs: enriched `summary`, then equivalent Python `code`.

## Invariants and report semantics

Accept any valid AnnData shape, including an empty object, because emptiness is itself structural evidence to report. Never mutate the input or force a backed matrix into memory. Access only metadata and cheap shape/dtype/storage properties. If a custom object cannot expose a property safely, record it as unavailable with a warning rather than scanning or coercing the matrix.

The key result is the data model inventory: dimensions, representation/storage, axis identity, available annotations/representations, raw state, and provenance. State explicitly that this does not validate counts, normalization, QC quality, annotation truth, or inferential design.

## Depth and test surface

The module hides bounded traversal of all AnnData slots, dtype/shape normalization, provenance extraction, narrative construction, references, version capture, and code rendering behind a one-input interface. Dependencies are in-process. Tests should assert strict JSON, deterministic bounds, sparse/dense/backed/view behavior without materialization, complete slot coverage, non-unique indexes, provenance, and generated-code equivalence.

## 2026-08-29 repair record

- New category: `openbio/single-cell/diagnostics`.
- Classification rationale: the atomic analysis is a read-only structural inventory of an existing AnnData object; it neither ingests data nor creates an input artifact, so `input` is misleading.
- Merge/delete decision: retain the standalone diagnostic. Do not merge it into loaders, QC, or every transformation and do not delete it, because the same bounded inventory is useful at any workflow seam.
