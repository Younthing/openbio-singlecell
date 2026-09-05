# Expert validation boundaries

## Request

Fix annotation merges blocked by provenance differences and audit the other nodes for similarly overconservative
validation. The user explicitly wants expert-facing wrappers to be as permissive as Scanpy and the selected backend.
This supersedes the source-provenance rejection and empty-Raw-axis rejection rules in the subpopulation-reanalysis spec.

## Behavior

- Honor the explicit expression source and algorithm parameters. Recommendations, count-likeness, provenance
  confidence, sampling recommendations, and interpretation preferences must not become additional runtime gates.
- Preserve checks needed for defined computation, structural/axis alignment, file safety, actual backend
  compatibility, and prevention of silent data loss. Do not silently change inputs or parameters to make them pass.
- Annotation values follow the existing error/keep_target/overwrite policy. Different input provenance is advisory.
  Preserve annotation provenance without falsely assigning one source's declaration to a mixed column.
- Fix additional restrictions only when native backend behavior or its implementation establishes that execution
  is supported. Audit all node families and record the limits of verification.
- Keep equivalent generated code consistent with runtime behavior.

## Verification

Use the existing public scientific-operation/node seams and equivalent generated-code tests. Work in vertical
red-green slices, then run relevant regressions, repository checks, and independent standards/spec review.

## Fixed point

`1b8d52b208c9c52aa8c84f5f21c607b8af68ab6d`; the worktree was clean before implementation.
