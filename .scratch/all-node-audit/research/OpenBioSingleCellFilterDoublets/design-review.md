# Filter Predicted Doublets: module design review

## Current module

This is a small but valid atomic transformation: consume an upstream boolean decision and remove marked observations while preserving all variables.

## Problems

- `astype(bool)` interprets every non-empty string, including `"False"`, as true.
- Missing column, missing predictions, empty names, empty input, and an all-removed result have no domain-specific validation.
- No removal summary, upstream-evidence limitation, references, versions, or equivalent source is exposed.

## Decision: enhance, do not merge

Keep separate from Scrublet. This preserves inspectability, permits other doublet callers, and prevents a stochastic analysis node from unexpectedly deleting observations. The module is intentionally shallow scientifically but owns an important destructive boundary and strict data contract.

## Target interface

- Input: `adata` and an advanced prediction-column name.
- Outputs: filtered `adata`, structured `summary`, equivalent `code`.
- Require a non-empty existing column with strict boolean dtype and no missing values.
- Reject empty input; allow and warn on an all-removed output, and warn when no cells are marked.
- Copy on subset and preserve observation identifiers/order and every gene.

## Report contract

Report input/retained/removed counts and percentages, the prediction column and its dtype, unchanged gene count, warning/limitation that calls depend on an upstream method/threshold, the AnnData subsetting reference, and runtime versions. The upstream caller's own summary carries its scientific method citation; this generic consumer must not falsely relabel predictions from another tool as Scrublet. Generated code must reproduce the same boolean and missing-value checks and the same zero-observation primary result. The structured report, rather than a hard gate, carries the warning that downstream analysis will generally need retained cells.
