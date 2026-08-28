# Pearson Residuals to Layer: module design review

## Current module

The node resolves `X` or a named layer, calls `scanpy.experimental.pp.normalize_pearson_residuals(..., inplace=False)`, copies the input `AnnData`, and stores the returned dense matrix in an output layer. The transformation is scientifically atomic and already avoids replacing the source matrix.

## Correctness and interface problems

1. The expression-source control defaults to `X` even though the operation requires raw UMI counts; the nested layer suggestion `counts` does not protect callers who accept the default.
2. Scanpy's `check_values=True` only warns on non-integer values. The wrapper currently accepts normalized or logged matrices and records no explicit count-suitability evidence.
3. `clip=0.0` is translated to `None`. This conflates a valid zero threshold with Scanpy's automatic `sqrt(n_obs)` rule and makes the resolved method invisible.
4. Sparse input becomes dense, but the interface has no preflight memory estimate or resource ceiling. A full-gene residual layer can exhaust memory even though the Scanpy tutorial recommends feature selection first.
5. A pre-existing output layer is silently overwritten, which can destroy a prior analysis state.
6. Zero-total cells/genes, empty axes, non-finite values, backed mode, and non-unique axes are not rejected before the model calculation.
7. The node discards Scanpy's returned method settings and exposes neither a report-ready `summary` nor equivalent `code`.

## Decision: keep and enhance; do not merge

Keep this node as the atomic conversion from one explicit count state to one residual representation. Do not merge it with highly variable gene selection or PCA. Pearson-residual HVG selection is deliberately chunked and should precede the dense transform on large data; PCA is a separate representation analysis with its own component and solver choices. Keeping these seams lets callers inspect the selected features, reuse the residual layer, and avoid allocating residuals when only HVG scores are needed.

## Target interface

- Inputs: `adata`; explicit `X`/named-layer count source; `theta`; clipping mode (`sqrt_n_obs`, `custom`, `none`); conditional custom clipping value; output layer.
- Advanced technical inputs: overwrite-existing opt-in and a dense-allocation ceiling in GiB.
- Outputs: copied `adata` with a new residual layer, structured `summary`, equivalent Python `code`.
- Visible scientific choice: expression source. `theta` and clipping remain exposed as advanced scientific controls because they change the model, even though the defaults are strongly recommended.
- Hidden policy: `check_values=True`, `inplace=False`, copy-on-write, exact finite/integer validation, resolved automatic bound, matrix conversion, and report/version construction.

The output-layer name, overwrite permission, and memory ceiling are advanced because they control storage and execution rather than the biological question. `raw` is deliberately absent from the source choices because its variable axis may differ from the current object and cannot safely populate a current-axis layer.

## Validation and state contract

- Require an in-memory nonempty `AnnData` with unique observation/variable identifiers. One-cell or one-gene inputs are allowed when Scanpy returns finite residuals; warn that `ddof=1` residual variance is unavailable for a singleton cell axis.
- Require a finite, non-negative source matrix with a positive grand total and no zero-total cell or gene. Fractional non-negative values and provenance inconsistent with counts proceed with explicit warnings; use `check_values=False` after this local audit so backend warnings cannot escape `-W error` or diverge from generated code.
- Require finite `theta > 0`, or explicitly supported positive infinity for the Poisson limit.
- Resolve automatic clipping to `sqrt(n_obs)` for reporting while still passing the semantically equivalent `None` to Scanpy; require custom clipping to be finite and non-negative; map “none” to positive infinity.
- Estimate one dense array conservatively at eight bytes per entry, disclose that peak use is higher because intermediates coexist, warn for every sparse-to-dense transition, and fail above the configured ceiling before copying or calling Scanpy.
- Reject a blank output-layer name, an output equal to the selected source layer, and an existing destination unless overwrite was explicitly requested.
- Preserve `X`, every source layer, `raw`, observations, variables, and unrelated annotations. The only scientific matrix change is the destination layer.
- Reject backed mode with an actionable instruction to materialize via `adata.to_memory()`; the output is necessarily dense and in-memory.

## Report contract

`summary` contains:

- methods text naming the negative-binomial offset model, shared `theta`, and clipping rule;
- results text suitable for a Methods/Results draft without calling the residuals counts or inference;
- key results for dimensions, input sparsity/count totals, resolved clipping bound, residual quantiles/range, clipped entries, residual-variance quantiles, output density/dtype, and estimated allocation;
- all resolved parameters, warnings and experimental-interface limitations;
- Lause, Scanpy, and AnnData references plus dynamic Python/openbio-singlecell/Scanpy/AnnData/NumPy/SciPy versions.

`code` is a self-contained function that copies the input, validates the same model and memory invariants, emits the same expert warnings, calls the public experimental Scanpy function with explicit `theta`, resolved `clip`, `check_values=False`, and `inplace=False`, refuses unintended overwrite, finite-postvalidates and writes the same destination layer, and returns the equivalent `AnnData`. It does not reproduce plugin-only history bookkeeping.

## Cohesion and coupling assessment

The node can be a deep in-process module: a small scientific interface hides source routing, strict count validation, clipping resolution, dense-allocation safeguards, Scanpy's return shape, layer preservation, reporting, citations, and version capture. No adapter seam is justified. Combining it with HVG selection or PCA would reduce locality because each operation has different memory behavior, outputs, and independently meaningful parameters.

## Verification plan

- Dense and CSR count sources from `X` and a named layer produce equivalent residual layers without mutating input/source matrices.
- Automatic, custom-zero, custom-positive, and no-clipping modes match direct Scanpy results; invalid `theta`/clip fail.
- Negative, non-finite, all-zero, zero-cell-total, and zero-gene-total inputs fail before Scanpy; fractional and provenance-inconsistent inputs run with warnings and are never called proven raw UMI counts.
- Existing/blank/same-as-source destinations, backed data, empty axes, and duplicate identifiers fail clearly.
- Sparse input reports densification and allocation; a deliberately low ceiling fails before allocation.
- `summary` is strict JSON, contains the resolved `sqrt(n_obs)` value and versions, and `code` compiles and reproduces the primary result and error invariants.

## Final audit amendment

The count-layer suggestion must also be the selected interface default; otherwise the visible default remains `X` despite the layer-preserving workflow. Legacy direct Python calls that omit the dynamic source may fall back to `X` only when the suggested layer is absent; normal UI calls carry their explicit selection. The Scanpy call must run on an isolated working `AnnData`, because the public experimental function materializes an `AnnData` view even with `inplace=False`. This preserves the caller object as well as its scientific matrices, and the generated function must use the same isolation seam.
