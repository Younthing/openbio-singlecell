# OpenBioSingleCellGenePanelScores — design review

## Final open-expert boundary

Expression history is descriptive only. The selected source defines the calculation; nonrecommended state becomes a
warning, not a gate. No history-tamper or Raw/current binding contract is introduced.

## Decision

**Enhance and narrow to one named panel per execution.** The current all-panels loop is a shallow batch operation
with dynamic outputs and partial-success semantics. One-panel scoring makes the operation atomic, gives one stable
output key, and lets failure/reporting describe the exact scientific object.

## Atomic interface

Inputs, in order:

1. `adata`;
2. local `gene_sets_file` and strict `resource_metadata_json`;
3. required `panel` identifier;
4. explicit expression `source` (`X`, named layer, or Raw snapshot), with normalized expression recommended and
   Raw/count-like scoring disclosed as an expert choice;
5. advanced resource `source_column`/`target_column`;
6. nonblank `output_key`, default `panel_score`;
7. advanced `ctrl_size`, where `0` resolves to the matched panel size and a positive value is explicit;
8. advanced `n_bins=25`;
9. `random_seed`;
10. explicit overwrite policy.

Outputs:

```text
adata, summary, code
```

Fix `ctrl_as_ref=False`, `gene_pool` to the complete selected expression universe, `use_raw=False`, and `copy=False`
on a private work copy. Resource parser and fingerprinting are shared internal implementation. Do not expose a list
of panels, output prefix, Condition, group, statistical threshold, or control-gene file until a second justified
adapter exists.

## Validation and result contract

Require unique observation/feature IDs, finite selected expression, one exact panel match, at least one matched gene,
valid bin/control settings, and enough usable gene-pool diversity for Scanpy. Duplicate resource edges collapse and
are counted; conflicting identifiers or blank values fail. Store one finite numeric column at
`output.obs[output_key]`; preserve all other input state and never partially add columns.

`summary` includes the strict resource identity/hash/metadata, requested and matched target accounting, exact source
and gene pool, resolved algorithm settings, score distribution and bounded extremes, methods/results prose,
warnings/limitations, references, and dynamic software versions. `code` defines an equivalent one-panel function
returning `(AnnData, summary_dict)` and must use the same resource and summary builders.

## Scientific boundary

The node produces per-cell descriptive scores only. Sample is the inference unit for downstream Condition analysis;
Technical batch is nuisance and is neither adjusted nor erased here. The report must explicitly prohibit treating
cell count as replicate count or interpreting a positive score as a probability/causal pathway state.

## Migration and tests

Legacy all-panel workflows cannot be translated without a panel choice. Migrations should preserve the node for
inspection and add a migration note requiring one panel/resource-metadata selection. Existing output-prefix links
cannot be guessed. Add `summary` and `code` only after this breaking choice is explicit.

Tests cover official 1.12.3 arguments, future-default-resistant `ctrl_as_ref=False`, one-panel selection, resolved
control size, deterministic seeded scores, exact source use, Raw/count-like warnings, missing/partial targets,
malformed and
changed resources, identifier collisions, output collision policy, sparse/dense inputs, immutability/atomic failure,
strict JSON, code compilation, and exact runtime/generated summary and score equivalence.

## Cohesion assessment

The revised module owns exactly one panel-score definition. Parsing, matching, control selection, Scanpy drift,
postconditions, and reporting remain hidden behind a compact interface. Batch orchestration belongs to the workflow,
not inside the scientific module.
