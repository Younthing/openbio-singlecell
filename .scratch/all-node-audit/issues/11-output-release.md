# Output adapters and final release audit

Type: task
Status: resolved

Blocked by: 10

## Initial audit focus

Atomic writes, CSV index identity, summary JSON/code persistence, preview of the complete reporting schema, workflow regeneration, release-manifest contract, and final registry audit.

## Comments

- 2026-08-28: Historical initial checkpoint: inventory and code audit were complete; per-node documents had not yet landed.
- 2026-08-28: Completed the two required pre-change documents for each of Preview Result, Save H5AD, Export CSV,
  and Save PNG. All four are terminal/UI adapters rather than analyses, so their design reviews explicitly retain
  zero graph outputs and do not invent `summary`/`code` values; they preserve disclosures produced upstream.
- 2026-08-28: Refactored all four adapters through one contained same-directory staging/validation/atomic-commit
  seam. Durable non-overwrite commits fail closed on a late collision; validation/writer failures clean staging
  files and preserve an old destination. PNG publication now requires Pillow verification plus full decoding and
  preserves exact bytes. H5AD exposes only portable gzip/lzf/none compression, defaults to gzip, and reopens a
  staged file backed/read-only to verify exact shape and axis identity. CSV now preserves the complete row index,
  assigns deterministic names to unnamed levels, rejects ambiguous/colliding headers and MultiIndex columns, and
  fixes UTF-8/LF/NA serialization.
- 2026-08-28: Output-focused evidence is 30 passing tests under `-W error`, including all three declared H5AD
  compression modes and real compressed H5AD
  round-trip, named/unnamed/multi-level row identity, corrupt PNG rejection, writer and validator failures, exact
  old-file preservation, staging cleanup, and late non-overwrite races. Focused Ruff, `py_compile`, release checks,
  and `git diff --check` pass. Pillow is now a declared core dependency with release manifest and license notice.
  At that checkpoint full workflow regeneration and repository-wide release checks were deferred until Batch09/10
  froze; the Answer below records the completed final gate.

## Answer

The output and release surface is closed. Preview, H5AD, CSV, and PNG remain cohesive terminal/UI adapters with
failure-atomic publication and no fabricated scientific outputs. Five generated workflows and covers, README,
release manifest, dependency/license notices, frontend payload, registry metadata, and the schema-driven generator
match the final interfaces. All six non-executing compatibility shims are publicly deprecated with actionable
replacement descriptions; the 87 active reporting nodes retain primary result + strict `summary` + equivalent
`code`, while nine adapters and six shims have explicit documented exemptions.

Final evidence: 102/102 node research directories contain exactly `official-usage.md` and `design-review.md`;
repository-wide Python tests pass 1151 with 4 expected skips under `-W error`; the migration suite passes 140/140;
release/registration/example/demo tests pass 37/37; Ruff, workflow regeneration check, JavaScript syntax check, and
`git diff --check` pass.
