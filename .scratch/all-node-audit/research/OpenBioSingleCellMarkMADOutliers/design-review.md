# Mark MAD Outliers: module design review

## Current module

The node performs one coherent annotation operation: apply the same MAD rule to one or more existing cell-level QC metrics and store the union as a boolean observation column. It copies `AnnData` and does not remove cells.

## Correctness and interface problems

- `median_abs_deviation(..., scale=1.4826)` divides by 1.4826, so the current optional normal scaling is reversed.
- The default grouping key is the technical `batch`, although co-varying QC distributions should normally be inspected within biological `Sample`/capture.
- NaN values, missing Sample labels, non-numeric columns, small groups, and zero MAD are not handled or disclosed.
- Only the union flag is exposed; users cannot audit which metric/Sample threshold caused a cell to be marked.
- There is no structured report or equivalent source.

## Decision: enhance, do not merge

Keep the node separate from destructive subsetting. Marking is an exploratory QC decision that should be inspected before removal. Do not merge it with Calculate QC: measurement and dataset-specific cutoff decisions are different operations.

## Target interface

- Inputs: `adata`, comma-separated metrics, grouping key (default `sample`, blank for global), `nmads`, one direction, output column, minimum finite group size, and optional normal-consistent scaling.
- Outputs: annotated `adata`, structured `summary`, equivalent `code`.
- Visible scientific choices: metrics, Sample/group key, `nmads`, and direction.
- Advanced controls: output column, minimum finite group size, and normal-consistent scaling.
- Errors: empty data, missing/non-numeric metrics, missing grouping column/labels, invalid output name, invalid group size, and unsupported direction.
- Missing metric values remain unflagged and are reported. Groups below the minimum finite size are skipped with a warning rather than assigned unstable thresholds.

## Report contract

Report the union count/rate; flag counts by metric and group; every median, MAD and effective threshold; skipped/zero-MAD groups; missing observations; resolved scaling convention; references; and SciPy/NumPy/Pandas/AnnData/plugin versions.

## Depth and coupling

The node hides robust-statistic calculation, grouped threshold bookkeeping, warnings, union construction, report generation, and source rendering behind a small interface. It consumes generic QC columns and produces a generic boolean annotation, so it stays decoupled from particular downstream filters.

## Open expert boundary decision

The lower bound for `minimum_group_size` is one. Singleton and other zero-MAD groups are computationally valid and are evaluated when the expert requests them; their degeneracy is disclosed in `summary`. Missing grouping labels remain a hard ambiguity because the node cannot assign those observations to a declared threshold population.

The shared group-label helper may be reused internally, but error text is part of this module's interface. MAD failures must refer only to grouped threshold reporting; caller-specific Scrublet wording must stay in the Scrublet module.
