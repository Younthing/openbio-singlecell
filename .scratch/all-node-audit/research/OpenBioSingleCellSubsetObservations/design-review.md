# Subset Observations: module design review

## Current module

The node parses comma-separated values, stringifies an `obs` column, performs exact membership matching, optionally inverts the mask, and returns an observation-axis copy. This is a coherent single transformation: it changes the observation set and nothing else.

## Correctness and interface problems

- With `invert=True`, missing annotation values are currently retained implicitly.
- Missing requested values, partially unmatched requests, and an empty output pass silently.
- String coercion is not disclosed, and mixed labels that collapse to the same string can be ambiguous.
- The comma-separated widget cannot express a category containing a comma.
- The result lacks a quantitative selection report and equivalent code.

## Decision: keep and enhance; do not merge with QC filters

Keep this node separate from cell-QC filters. QC filters apply numerical acceptance thresholds and report technical quality decisions; this node selects a scientific population from a declared annotation. Combining them would weaken both report semantics and make the interface shallower.

## Target interface

- Preserve legacy positional inputs: `adata`, `column`, `values`, `invert`.
- Append advanced `missing_policy` with `exclude | include | error`, default `exclude`.
- Keep exact matching fixed and hidden. Do not add regex, substring, case-folding, or fuzzy matching modes.
- Outputs: subset `adata`, `summary`, `code`.

`column`, `values`, and `invert` remain visible because they define the selected population. Missing policy is advanced but explicit because it can alter the retained cells. Copy semantics, observation-order preservation, and alignment behavior are fixed implementation invariants.

## Resolved matching semantics

Record missingness before conversion, then compare exact string representations:

```python
not_missing = series.notna()
matched = not_missing & series.astype("string").isin(selected_values)
mask = matched if not invert else not_missing & ~matched
```

For `missing_policy="include"`, add missing observations to either positive or inverted selection. For `error`, fail if any annotation is missing. The default `exclude` never retains missing observations, including under inversion. Detect distinct non-missing labels that collide after string conversion and fail rather than select ambiguously.

## Invariants and failure modes

- Require a nonempty column name, existing annotation, and at least one nonempty requested value.
- Preserve input observation order; never sort by category or index.
- Fail when positive selection matches none or whenever the final output is empty.
- If all requested values are absent under inversion, keep the eligible observations but warn; if only some values are absent, report each absent value.
- Duplicate `obs_names` do not make this position-based slice incorrect, so do not reject them. Warn that the result is unsafe for later cross-object annotation merging.
- Preserve the variable axis exactly and rely on AnnData slicing to align `X`, layers, `obs`, `obsm`, `obsp`, and `raw`.

## Report contract

`summary` includes input/retained/removed counts and rates, requested values, per-value input and selected counts, missing observations, unmatched values, invert and missing-policy settings, output shape, duplicate-identifier warning, AnnData/pandas references, and dynamic software versions. It must describe the selected annotation as a user-declared population and make no claim that labels are biologically validated.

`code` parses the same comma-separated form, enforces the same missing and collision policies, preserves order, copies the slice, and returns the equivalent primary result.

## Cohesion and coupling assessment

The node remains atomic and high-cohesion because only the observation set changes. A small exact-match interface hides missingness handling, aligned-slot slicing, validation, and reporting. Annotation transfer remains a separate module because it does not change the observation set.

## Verification plan

- String, categorical, numeric, and boolean annotation columns, including string-representation collisions.
- Missing annotations under all three policies and both positive/inverted selection.
- Partially and wholly unmatched values, empty output, and labels containing commas as a disclosed input-format limitation.
- Duplicate `obs_names` produce the documented warning without changing position-based slicing.
- Observation order and aligned `X`, layers, `obsm`, `obsp`, and `raw` contents remain correct.
- Inputs remain unchanged; summary is strict JSON; generated code compiles and returns an equivalent AnnData.

## Open expert-boundary decision (2026-08-28)

Replace the no-match and empty-result gates with warnings. This node's atomic job is to apply the expert's exact mask,
not to decide whether the resulting population is large enough or scientifically useful. Keep column existence,
matching-mode ambiguity, and explicit missing-policy errors as true interface constraints. Runtime and generated code
must both return the same zero-observation AnnData.
