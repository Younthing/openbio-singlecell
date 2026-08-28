# Current-only compatibility cleanup

## Objective

Expose only the current OpenBio Single-Cell node interfaces. Saved workflows using removed node IDs or historical
schemas are intentionally unsupported and receive no automatic rewrite, compatibility alias, or load-time shim.

## Remove

- Every registered node whose only purpose is legacy workflow loading, migration guidance, or deprecated mechanical
  adaptation.
- The browser-side workflow migration Module, its graph-load hook, and its migration-only tests.
- Optional-dependency aliases, release metadata, documentation, and test expectations that exist only for former
  names or interfaces.
- Legacy interface branches and helpers with no current caller after those seams are deleted.

## Retain

- Current analysis nodes and current typed-result/configuration/output Adapters that have an independent present-day
  use.
- Versioned scientific artifact/file formats, explicit schema-version validation, backend capability checks, and
  dependency version constraints required to interpret current results safely.
- Historical research records under `.scratch/all-node-audit/`, updated where necessary to state the final removal
  disposition instead of presenting compatibility retention as current policy.

## Acceptance criteria

- No deprecated or migration-only node is registered.
- No frontend code mutates a workflow during graph loading.
- Removed IDs, legacy migration entry points, and compatibility-only dependency aliases are absent from shipped code,
  generated workflows, release metadata, and current-facing README text.
- Registry and release tests lock the current-only node set and reject reintroduction of deprecated schemas.
- Full Python tests, current frontend tests, Ruff, generated-workflow verification, JavaScript syntax checks, and
  `git diff --check` pass.

