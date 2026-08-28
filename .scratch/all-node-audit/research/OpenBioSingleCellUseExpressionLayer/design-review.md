# Use Expression Layer: module design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. Expression source is selected
> explicitly by each consuming analysis; no deserialization shim or workflow migration is retained.

## Current module

The node copies its input, copies one named layer into `X`, and returns the result. Its interface is small, but the implementation is an almost direct assignment. More importantly, it moves the expression-source decision away from the analysis that uses the matrix and turns that decision into implicit workflow state.

## Correctness and interface problems

- Downstream nodes see only `X` and cannot know which layer produced it without reconstructing workflow history.
- A later node may overwrite `X` again, so correctness depends on ordering rather than a local invariant.
- `raw` and all other layers remain unchanged, allowing inconsistent representations to coexist without an explicit contract.
- The full-object copy plus matrix copy spends substantial memory for no scientific computation.
- Adding more switches such as in-place mode or copy policy would enlarge a shallow interface instead of creating depth.

## Decision: migrate consumers, then delete

The deletion test shows that the module does not earn its seam: deleting it does not remove the need to select an expression representation; that choice simply reappears at every real analysis. The correct locality is the consumer module. Each consumer should expose one explicit expression-source interface and hide any conversion required by its implementation.

This node must not be removed until all migration prerequisites are satisfied:

1. Inventory every registered node and generated workflow that reads `X`; classify whether it requires counts, normalized values, scaled values, or an arbitrary expression representation.
2. Add an explicit `ExpressionSourceSpec` to every consumer for which the representation is selectable. Count-dependent consumers must validate and report the selected source.
3. When a wrapped package only accepts `X`, create a private working `AnnData` inside that consumer and transfer only the intended results back to a copy of the original object.
4. Migrate saved/example workflows so no connection relies on `UseExpressionLayer` as a state switch; add a repository test that no supported workflow contains it.
5. Test every migrated consumer directly through its interface for `X`, layer, and `raw` behavior as applicable, including sparse input and preservation of the caller's other expression states.
6. Only after those checks pass, remove the node from registration and document the saved-workflow migration in release notes.

## Compatibility-period contract

If saved-workflow compatibility requires one deprecation cycle, retain the current positional inputs only and rename the display label to `Deprecated: Materialize Layer into X`. Do not add more behavior. Append `summary` and `code` outputs after the primary `adata` output.

The temporary implementation must reject an empty/missing layer, copy rather than mutate the input, preserve cell/feature axes and `raw`, and disclose a deprecation warning. `summary` must report source layer, shape, dtype, sparse/dense storage, whether a previous `X` was replaced, that no biological transformation was performed, the AnnData reference, and dynamically collected Python/openbio-singlecell/AnnData/NumPy/Pandas versions as applicable. `code` must reproduce the copy and replacement exactly.

## Parameter policy

- Visible while deprecated: `adata`, `layer_name`.
- Hidden: copy semantics are fixed; mutation and performance controls are not exposed.
- Removed rather than generalized: destination slot, in-place mode, fallback source, and automatic layer discovery.

## Cohesion and coupling assessment

The current node is internally cohesive only as a storage assignment, but it increases coupling across the whole workflow by changing the implicit meaning of `X`. Moving source resolution into each consuming module creates a deeper interface: callers name the scientific source once, while validation, temporary-object adaptation, and provenance remain local to the consumer. Final decision is **delete after migration**, not merge into another standalone state-switch node.

## Verification plan

- Repository inventory proves that every surviving `X` consumer either has fixed, documented `X` semantics or an explicit expression source.
- Consumer tests cover selected `X`, layer, and `raw` paths as applicable and confirm that private work objects do not alter caller expression state.
- Compatibility-node tests, if temporarily retained, cover missing/blank layers, sparse and dense matrices, input immutability, unchanged `raw`, strict-JSON summary, and compiling/equivalent code.
- Workflow and registry tests prove that no supported workflow references the node before its registration is removed.

## Open expert-boundary decision (2026-08-28)

During the compatibility period, remove the non-empty-axis gate and report an empty-materialization warning. No new
options are added. This preserves the small adapter interface while allowing every AnnData state for which the public
layer assignment is defined.
