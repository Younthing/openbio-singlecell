# Core Study Parameters: module design review

## Current module

The node converts one UI-owned, schema-versioned JSON value into six reusable strings. `CoreStudyParameters.from_json` centralizes exact-field, type, version, and whitespace validation.

## Decision: keep unchanged as a configuration adapter

The module is cohesive and deep enough for its role: callers learn one versioned payload while parsing and validation remain local. Deleting it would duplicate column/contrast widgets across Study-level workflows; merging it with any analysis would couple configuration reuse to one method.

It performs no dataset access, transformation, or scientific analysis. Therefore it must not gain artificial `summary`/`code` ports. Its six string outputs are the primary configuration values.

## Interface and parameter policy

- Visible input: `study_parameters_json`, rendered by the dedicated Core Study Parameters widget. All six scientific names/levels remain visible inside the JSON because silently chosen Sample, Condition, or contrast values would change interpretation.
- Advanced inputs: none. Hiding any field would make downstream scientific assumptions implicit.
- Hidden behavior: exact schema-field set, strict JSON decoding, `schema_version=1` dispatch, type checks, and trimming checks.
- Outputs: `sample_column`, `condition_column`, `batch_column`, `annotation_column`, `reference`, and `comparison` in stable order.

`batch_column`, `reference`, and `comparison` may be empty when a particular downstream node does not need them. Dataset- and method-specific requirements remain the responsibility of the consuming node, which has the necessary context to issue an accurate error.

## Cohesion, coupling, and invariants

The node owns configuration syntax and vocabulary, not validation against an `AnnData`. This seam keeps configuration reusable and prevents coupling to input order or dataset state. Require exact versioned JSON and trimmed strings; do not infer missing fields, rename columns, swap reference/comparison, or conflate Sample with Technical batch.

Dependencies are pure in-process Python/JSON. The interface is the test surface: tests should cover the canonical payload, stable output ordering, all schema/type/whitespace errors, and intentionally empty optional strings.

## Open expert-boundary decision (2026-08-28)

Keep exact schema/version/type validation but remove the whitespace-policy gate. Trimming would silently change an
expert's exact column or Condition-level identity, while rejection imposes a style rule unrelated to JSON decoding.
Preserving strings verbatim gives consumers one unambiguous value and keeps scientific semantics at their own seam.
