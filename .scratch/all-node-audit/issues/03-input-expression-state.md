# Input, study, and expression-state nodes

Type: task
Status: resolved

Blocked by: 02

## Initial audit focus

Strict multi-Sample loading, Raw snapshot semantics, stable feature identifiers, collision policies, expression-source provenance, and data preparation report/code outputs.

## Comments

- 2026-08-28: Landed the required official-usage and design-review documents for all six input/configuration/structure nodes before implementation. Enhanced four loaders, kept Core Study Parameters unchanged as a configuration adapter, and upgraded AnnData Summary to the standard `summary + code` contract. Strict 10x Study inventory, stable feature metadata restoration, count/identity checks, portable provenance, and bounded structural reporting are covered by 25 warning-strict tests and an independent review with no remaining P0/P1/P2 findings.

- 2026-08-28: Initial inventory and code audit complete; per-node pre-change documents pending.

- 2026-08-28: Landed both required pre-change documents for the six expression-state/data nodes, implemented their strict scientific contracts and report/code outputs, and added frontend migration for the two schemas whose appended widgets changed old UI-workflow semantics. Independent final review drove fixes for semantic scalar collisions, nullable dtype preservation, typed value counts, and unified missing-value reporting.

## Answer

All twelve input/configuration/expression-state nodes have been audited and completed. Loaders now enforce stable observation/feature identity, count semantics, complete Sample inventories where applicable, portable provenance, and bounded strict structural reporting. Core Study Parameters remains a configuration adapter; AnnData Summary emits the standard enriched summary and equivalent code.

Expression-state ownership is now explicit. Snapshot Expression is the sole canonical creator of the user-declared Raw snapshot and `layers["counts"]`; it preserves `X`, requires explicit overwrite for occupied destinations, and reports that snapshot timing and full-gene completeness are not programmatically proven. Analysis history is not a late-snapshot gate, and explicit Raw selection is sufficient without a Raw/current binding test. Generic layer materialization and regex feature-name replacement remain only as deprecated compatibility nodes pending consumer migration. Observation subsetting has explicit missing-value behavior, annotation transfer has explicit conflict/provenance/dtype policy, and GTF annotation preserves stable IDs while adding gene symbols without filtering features.

Every retained processing/reporting node in this batch exposes strict JSON `summary` and executable equivalent-function `code`. Validation passed with Ruff, 27 warning-strict data tests, 17 frontend tests, five generated workflow checks, and the full warning-strict Python suite (`214 passed`). The independent final audit reported no remaining P0/P1/P2 findings.
