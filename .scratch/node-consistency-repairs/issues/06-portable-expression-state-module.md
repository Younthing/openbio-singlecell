# Remove inferred expression state

Type: task
Status: resolved
Blocked by: 05

## Question

Should OpenBio infer an expression state from workflow history, AnnData storage slots, or mutable metadata?

## Acceptance criteria

- Runtime, summaries, generated code, typed artifacts, and model metadata contain no inferred expression state or
  state-evidence fields.
- Explicit X, Raw, and layer source selection remains unchanged.
- Raw remains a user-owned AnnData snapshot without OpenBio count/normalization semantics or new gates.
- Scientific operations retain only the numeric and structural validation they actually need.
- CellTypist input mode and the pySCENIC external-manifest declaration remain because they are explicit behavior or
  import contracts, not inferred state.

## Answer

Delete the abstraction. The visible workflow and explicit source already express user intent; reconstructing a
second state model from mutable history adds duplicated knowledge and contradictory behavior. Scientific callers now
validate their actual input requirements directly, while Raw follows ordinary AnnData semantics.

## Comments

- 2026-08-31: The user rejected both inferred `state` and `evidence` as redundant and explicitly required OpenBio to
  preserve user intuition instead of policing `adata.raw`.
- 2026-08-31: Removed the shared resolver, caller-local/generated copies, inferred metadata and state-only warnings;
  retained explicit sources, real numeric gates, CellTypist behavior selection, and pySCENIC manifest declarations.
- 2026-08-31: Verification passed: 1,371 Python tests (4 skipped), 14 frontend tests, Ruff, five generated workflows,
  import checks, JavaScript syntax, and `git diff --check`.
