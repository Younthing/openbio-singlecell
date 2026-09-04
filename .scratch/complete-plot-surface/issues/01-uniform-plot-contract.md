# Establish the companion Plot contract

Type: task
Status: resolved

## Scope

- Freeze the schema/category/output conventions from the spec.
- Reuse existing PlotResult, report, Worker, Preview, Save PNG, and Persist seams.
- Add only proven shared validation/render helpers needed by at least two adapters.
- Move Sample Composition Plot to `differential-abundance` and lock the convention in registration tests.

## Acceptance criteria

- No universal Plot node, DSL, new dependency, default producer plot, or compatibility alias.
- Public seams and helper behavior have red-green tests.

## Comments

- 2026-09-04: Claimed after the user approved domain-specific companion Plot nodes and complete evidence-backed coverage.
- 2026-09-04: Resolved before domain implementation; 29 focused contract/category tests passed.

## Answer

- Companion plots use exact scientific inputs, remain in the producer's scientific category, and expose only
  closed same-contract view choices followed by `plot`, `summary`, and `code` outputs.
- Existing PlotResult, reporting, File-artifact codecs, Preview/Save PNG/Persist, and one-shot Worker seams are
  reused directly. No universal plotting DSL, generic renderer registry, dependency, producer-side default plot,
  or compatibility alias was introduced.
- Sample Composition Plot now belongs to `differential-abundance`; registration tests reserve `visualization` for
  genuinely cross-domain AnnData views.
