# Remove backward-compatibility surfaces

Type: task
Status: resolved

## Question

Which public nodes, aliases, browser hooks, migration helpers, tests, and release references exist solely to support
historical workflows, and can they be deleted without weakening a current scientific Module or Adapter?

## Answer

Eight node IDs, the browser graph-load migration path, hidden direct-call aliases, old analysis-history readers, and
the `ora` extra existed only for historical interfaces and were deleted. The resulting registry contains 94 current
nodes and no deprecated entries. Current scientific artifact/file-schema validation and backend capability checks
remain because they validate present inputs and execution dependencies rather than historical OpenBio workflows.

## Comments

- 2026-08-28: Claimed after the user explicitly removed the backward-compatibility requirement. The cleanup uses a
  current-only interface: old node IDs and schemas are unsupported rather than migrated or represented by shims.
- 2026-08-29: Removed eight migration/compatibility-only node IDs, browser graph-load migration, direct-call aliases,
  old analysis-history readers, and the compatibility-only `ora` extra. The registry is locked to 94 current nodes
  with zero deprecated entries. Retained only current artifact/file-schema validation and backend capability checks.
- 2026-08-29: Verification passed: 1137 Python tests (4 skipped), 14 frontend tests, Ruff, JavaScript syntax,
  five generated workflows, and `git diff --check`.
- 2026-08-29: Final independent audit found and closed two residual defaults: UMAP provenance now requires an explicit
  matching `key_added`, and `SCVIModel` requires explicit current `source`, `count_source_state`, and
  `count_source_state_evidence`. The full Python suite remained at 1137 passed and 4 skipped after the change.
