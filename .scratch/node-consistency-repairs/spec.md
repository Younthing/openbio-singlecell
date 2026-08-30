# Node consistency repairs

## Objective

Remove the concrete semantic inconsistencies found after the current-only 94-node refactor while keeping every
scientific node atomic. The public node Interface remains the node schema plus `execute`; analysis/report nodes keep
their primary output(s), `summary`, and `code` outputs. Historical schemas and implicit Raw fallbacks remain out of
scope.

## Confirmed public seams

- `SummaryResult.summary` is the report artifact seam. Every analysis/report node must return the same canonical
  report fields and JSON-compatible value types, regardless of which internal builder it uses.
- A node's explicit expression-source input is the expression-source seam. All quantities combined in one result or
  plot must use that selected source.
- A named graph selected through a node input is the graph seam. Structurally equivalent sparse matrices must receive
  the same base validation before method-specific rules are applied.
- Generated `code` is the reproducibility seam for the scientific result. It must reproduce the equivalent function
  without inventing or mutating OpenBio-only history metadata.

These seams were specified by the user's required `summary`/`code` outputs and explicit Raw selection policy; tests
observe only those public interfaces.

## Required pre-change evidence

The official-usage and design-review documents created during `.scratch/all-node-audit/research/` remain the required
pre-change record for every affected node. Each repair ticket names the relevant records and must update the design
review when its final behavior changes.

## Repair order

1. Deepen the report-contract Module and repair Milo, scCODA/tascCODA, and scVI differential-expression summaries.
2. Make QC plot metrics a coherent family from one explicit source.
3. Make named-graph base validation consistent, beginning with stored sparse zero-diagonal entries.
4. Remove inferred expression state and generated-history behavior, then delete proven dead wrappers.
5. Align node taxonomy and user-facing descriptions after behavioral seams are stable.

## Acceptance criteria

- Invalid analysis summaries fail at the shared report seam; all 85 analysis/report nodes satisfy the canonical
  report contract.
- Milo, scCODA, tascCODA, and scVI differential-expression summaries expose string `methods`/`results`, mapping
  `key_results`, references, warnings, limitations, parameters, and Python/OpenBio software versions.
- QC plots never combine metrics resolved from different expression sources.
- All named-graph consumers agree on base sparse-graph validity; explicit stored zeros do not create a path-specific
  failure.
- Runtime, summaries, artifacts, models, and generated standalone code contain no inferred expression-state model.
- No unused compatibility or pass-through helper remains in the touched paths.
- Targeted tests are written red-first at the confirmed seams, then the full Python/frontend/workflow/lint checks pass.

## Deferred architectural opportunity

A typed `GeneSetResource` Adapter could replace six repeated file/metadata input groups, but it changes several public
node Interfaces and needs its own research and migration-free design decision. It is not bundled into these correctness
repairs.

Issue 06 resolved the former expression-state consolidation idea by deleting the redundant abstraction instead.
