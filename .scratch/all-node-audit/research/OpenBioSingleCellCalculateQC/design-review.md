# Calculate QC Metrics: module design review

## Current module

The module performs one coherent operation: annotate configurable QC gene sets and calculate observation/variable QC metrics. It copies `AnnData`, so the external seam is deterministic and side-effect-free for callers.

## Problems

- The interface implicitly reads `adata.X` and labels its values as counts. Once `X` has been normalized or logged, this is scientifically misleading even though Scanpy accepts the matrix.
- The implementation records parameters and warnings in hidden history but exposes neither a report-ready result nor reproducible code.
- The key result—what genes matched and what the QC distributions look like—is not surfaced.

## Decision: enhance, do not merge

Keep this as an atomic analysis node. Merging it with filtering would couple measurement to a dataset-specific decision and make it impossible to inspect QC before choosing thresholds. `QC Plots` remains a separate read-only report because users may recalculate or inspect metrics without filtering.

## Target interface

- Inputs: `adata`; expression source; include flags; organism/annotation-dependent patterns; `percent_top`; `log1p`.
- Outputs: transformed `adata`, structured `summary`, equivalent Python `code`.
- Invariants: at least one cell and gene; non-empty mitochondrial prefix; enabled patterns must be non-empty; raw/count suitability is disclosed; input is not mutated.
- Errors: missing selected layer/raw, empty enabled patterns, empty matrix, invalid `percent_top`.
- Availability: when a configured QC feature set matches no genes, its zero-valued Scanpy observation columns are removed and the report marks the evidence unavailable; an annotation mismatch must not be presented as a biological zero.

Frequently chosen include flags remain visible. Prefix/patterns, `percent_top`, and log-derived annotations are advanced. The expression source stays visible because it changes the scientific meaning.

## Depth and dependencies

Dependencies are in-process (`AnnData`, Scanpy, NumPy/Pandas). The module hides gene-set construction, Scanpy routing, unsupported top-rank handling, result summarization, citations, and version capture behind one node interface. No adapter seam is justified.

## Open expert boundary decision

- Do not inspect `analysis_history`, prove a full-gene snapshot, or require Raw to equal current `X`/`var`; the selected source is the user's declaration.
- Raw with a broader variable axis is handled in an isolated Raw-feature working object. Observation QC metrics return to the caller's current object; variable metrics are mapped by unambiguous feature identity and mapping coverage is reported.
- Signed and fractional finite inputs produce suitability warnings and explicit “expression sums, not verified counts” wording. The node blocks only non-finite data or a transform/output domain that cannot be computed safely.
- These are report disclosures, not tests that simulate history tampering or assert a Raw provenance binding.

For wrapper-derived `percent_top` columns above the selected feature count, use 100% only when the selected expression sum is non-zero. Preserve a missing value for a zero denominator and report how many cells were undefined. Runtime and generated code must use the same rule.
