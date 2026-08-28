# Log1p: module design review

## Current module

The module copies AnnData, applies `scanpy.pp.log1p` to `X`, and by default assigns the transformed copy to `output.raw`. Its visible `set_raw=True` switch couples a numerical transformation to a separate storage checkpoint.

## Correctness and interface problems

- The default produces a log-normalized `raw`, contradicting the repository's **Raw snapshot** definition and risking invalid count-dependent downstream analysis.
- It has no finite/non-negative check, no explicit double-log guard, and no disclosure when normalization provenance is absent.
- `set_raw` is not part of the upstream Scanpy method; it is an OpenBio side effect that broadens the interface without strengthening the transformation.
- Existing `raw` can be overwritten silently; absent `raw` is manufactured even when no full-gene count checkpoint was intended.
- No report exposes the formula, resolved base, value changes, state limitations, references, or runtime versions.
- No equivalent source is returned.

## Decision: keep and narrow to one atomic transformation

Keep the module as the atomic active-matrix log transformation and remove raw ownership. Do not merge it permanently into `NormalizeTotal`: log1p is independently useful after other normalization methods, while some callers legitimately require linear normalized values. `NormalizeToLayer` may compose both operations for a reusable derived layer, but its output-state contract is different. Shared numerical validation/code rendering can live behind an internal seam so the implementations do not drift.

The module remains useful under the deletion test because removing it would reproduce log-state validation, double-transform protection, sparse handling, reporting, and equivalent-code behavior in callers that intentionally transform `X`. Its interface becomes deeper by deleting the unrelated storage switch.

## Target interface and parameter policy

- Input: `adata` only; active `X` has fixed, explicit meaning.
- Outputs: copied `adata`, structured `summary`, equivalent Python `code`, in that order.
- Fixed hidden policy: natural logarithm (`base=None`), pseudocount one as defined by `log1p`, copy-on-write, no chunk/backed mutation, no layer or `obsm` routing, and no creation/replacement/deletion of `raw` or layers.
- Removed: user-facing `set_raw`. Snapshot ownership belongs exclusively to `OpenBioSingleCellSnapshotExpression`.

The base is fixed because natural log is the established repository convention and an extra numeric control offers little routine leverage. Source selection is not added here: a generic source-plus-destination interface would duplicate `NormalizeToLayer` and turn this atomic module into another expression-state router.

Saved-workflow migration must remove the `set_raw` widget. During a compatibility period, the Python method may accept a deprecated optional argument outside the schema: `False` maps to the new behavior, while `True` must raise an actionable instruction to insert `SnapshotExpression` before preprocessing. It must never silently create a transformed `raw`.

## Invariants and failure modes

- Require nonempty, in-memory AnnData and finite active `X` strictly greater than `-1`, the real-valued mathematical domain of `log1p`.
- Warn for values in `(-1, 0)` and whenever `uns["log1p"]` or OpenBio history indicates prior transformation; provenance informs interpretation but does not veto an explicit expert operation.
- If normalization/log state cannot be established, proceed with a limitation; matrix values are insufficient proof.
- Preserve cell/feature identity and order, annotations, every layer, and the exact presence/content of existing `raw`.
- Preserve sparsity and do not mutate the input object.

## `summary` and `code` contract

The report states the exact natural-log formula, cells/features transformed, sparse/dense storage and dtype, nonzero entries, bounded input/output value distributions, state evidence, and raw/layer preservation. Its result sentence must call the output a log-transformed expression representation rather than normalized counts. Include Scanpy/AnnData method/software references and dynamic Python/openbio-singlecell/Scanpy/AnnData/NumPy/SciPy/Pandas versions.

Generated code performs the same domain validation and warnings, copies the object, calls the public Scanpy function with `base=None`, postvalidates finite output, and returns the transformed AnnData without plugin-only history.

## Cohesion, coupling, and seam placement

The external seam exposes one operation with no storage decision. Scanpy and AnnData are in-process dependencies and need no adapter. Raw snapshot construction remains local to one separate module, eliminating an ordering-dependent convention from this interface. The resulting module has high cohesion: every invariant and result concerns exactly one `log(1+x)` transformation of active `X`.

## Verification plan

- Dense and sparse non-negative matrices match `scanpy.pp.log1p` and do not mutate the caller.
- Existing `raw` and all layers remain value/identity equivalent; absent `raw` remains absent.
- Empty, non-finite, `x <= -1`, and backed inputs fail clearly; `-1 < x < 0` and proven already-logged inputs execute with prominent warnings.
- Unproven normalization state is disclosed without inventing certainty.
- Deprecated `set_raw=False` compatibility, if retained, matches the target operation; `True` raises and never writes `raw`.
- Strict-JSON `summary` and compiling generated `code` reproduce the primary scientific state and invariant failures.
