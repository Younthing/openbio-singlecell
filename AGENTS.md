# Engineering workflow

## Core policy

- Use Ponytail as the default implementation filter, after understanding the real flow: reuse existing code, then stdlib/native capabilities, then installed dependencies, and only then write the minimum necessary code.

- Optimize for the simplest design, not the smallest textual diff. Readability and separation of concerns take priority over backward compatibility or migration support unless the task explicitly requires them.

- A larger focused change is acceptable when it is necessary to make the requested behavior clearer. It is not permission for unrelated cleanup.

- Separate responsibilities at evidence-backed seams. Do not create abstractions merely to make the code look layered.

## Feature routing

Choose the shortest flow that resolves the actual uncertainty:

```text

Behavior and boundaries are clear

  -> tdd -> code-review

Requirements, domain language, or boundaries are materially unclear

  -> grilling + domain-modeling (the grill-with-docs flow)

  -> shared understanding -> tdd -> code-review

A decision requires running or seeing something to settle it

  -> prototype one explicit question -> confirm the conclusion with the user

  -> implement production behavior with tdd -> code-review

```

- Do not grill or prototype a well-scoped change.

- Prototype code is evidence, not production code; carry forward the conclusion, then implement it under production constraints.

- At the start of implementation, capture the review fixed point. Use a user-specified base when provided; otherwise use the current `HEAD` and account for any pre-existing worktree changes.

- `code-review` must evaluate both repository standards and the originating request/spec. Review the complete change, including staged, unstaged, and relevant untracked files; do not commit merely to make the review possible.

- Report review findings outside the requested scope, but do not fix them unless the user expands the task.

## Codebase discovery

- Start with `rg`, direct source reading, existing tests, `CONTEXT.md`, and relevant ADRs.

- Use GitNexus only when the change crosses modules, callers or dependants are unclear, the architecture is unfamiliar, or the blast radius is uncertain.

- Use `query` to find execution flows, `context` for a symbol's callers/callees, `impact` before changing a shared symbol, and `detect-changes` to inspect a broad diff.

- Treat the graph as navigation evidence, not ground truth; verify relevant results in source with `rg`.

- Keep GitNexus opt-in. Refresh with `gitnexus analyze --index-only .` only when its index is stale; do not install its generated rules, skills, or hooks.

## Before changing code

Before editing, be able to:

- State the requested behavior change in one or two sentences.

- Point to the closest existing implementation or repository convention.

- Name the files that actually need to change.

- Separate required work from merely desirable improvements.

If a material product or domain decision remains and cannot be discovered from project evidence, ask the user. Do not ask for confirmation when the requested behavior and existing convention already make the answer clear.

## Scope discipline

- Every changed file and hunk must be necessary for the requested task.

- Do not fix unrelated problems or perform unrelated cleanup, formatting, renaming, or dependency updates.

- Reuse an existing pattern instead of introducing a new abstraction.

- Do not generalize for hypothetical future use cases.

- Add a dependency only when the task cannot reasonably be completed with the codebase, standard library, native platform, or existing dependencies.

- Do not add backward-compatibility shims, migration paths, or dual implementations unless explicitly required. This does not authorize silent data loss or destructive schema/data changes.

## Defensive code

- Scientific nodes are expert-facing wrappers. Follow the selected backend's executable input and parameter
  behavior; do not promote its recommendations or warnings into additional hard errors. Expression-state guesses,
  provenance confidence, recommended sample sizes, plotting preferences, and documentation completeness are
  advisories. Verify a proposed scientific hard gate against backend behavior before adding it.

- Add validation or guards only for realistic states supported by project evidence or an explicit trust boundary.

- Do not add checks merely because a value could theoretically be null, missing, malformed, or unexpected.

- Preserve and expose invariant violations instead of silently masking them.

- Do not add fallback behavior without a concrete requirement.

- Never simplify away security, data-loss prevention, accessibility, or required trust-boundary validation.

## Tests

- In TDD, work in vertical red-green slices: one behavior test followed by only enough implementation to pass it.

- Test behavior through the narrowest existing public seam, not implementation details.

- Ask the user to choose a seam only when that choice would materially change the public interface or design; use an evident existing seam without ceremony.

- Use the narrowest existing test command that validates the change.

- Add a test only when it directly protects changed behavior or reproduces the fixed bug.

- Do not expand the matrix for hypothetical edge cases or rewrite tests to match a preferred implementation.

## Diff review

Before finishing, verify:

- Every changed hunk can be justified from the original request; if reverted, a named requirement would stop being satisfied.

- The solution follows repository conventions rather than generic best practices.

- A simpler solution would not satisfy the same behavior equally well.

- No unrelated improvement survived into the diff.

Ask: **“Did I solve the requested problem, or did I start improving the project?”** If it is the latter, remove the unrelated improvements.

## Agent skills

### Issue tracker

Issues and specs are tracked as local Markdown under `.scratch/<feature-slug>/`. See `docs/agents/issue-tracker.md`.

### Domain docs

This is a single-context repository. See `docs/agents/domain.md`.
