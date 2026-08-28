# OpenBioSingleCellAUCellScores — design review

## Final open-expert boundary

The expression-state resolver is a reporting seam, not an authority seam. Any finite explicitly selected source is
executed, while nonrecommended scale is disclosed as a warning; mutable history and Raw/current equality are not
runtime gates. This final boundary supersedes any stricter validation language below.

## Decision

**Enhance and retain as one atomic AUCell observation-scoring module.** Migrate to decoupler 2.2 `mt.aucell`, remove
the inert seed, and add strict resource/state/output contracts. Do not merge it with GSVA or
panel scoring: those methods have different score definitions, parameters, ranges, and reporting limitations.

## Atomic interface

Inputs, in order:

1. `adata`;
2. local `gene_sets_file` (`.csv`, `.tsv`, or `.gmt`);
3. strict `resource_metadata_json` with nonblank `name`, `version`, `date`, `organism`, `identifier_namespace`,
   `scope`, `license`, and `citation`;
4. explicit expression `source` (`X`, named layer, or Raw snapshot), with normalized continuous expression
   recommended and Raw/count-like use disclosed as an expert choice;
5. `source_column` and `target_column` for tabular resources;
6. `min_targets` after exact expression-universe intersection;
7. optional `n_top_features`, where `0` means the documented 5% rule;
8. advanced bounded `batch_size` and `max_output_rows`/memory guard;
9. nonblank `output_key`, default `aucell_scores`.

Outputs, in order:

```text
adata, summary, code
```

The score frame is stored as a named pandas DataFrame in `output.obsm[output_key]`, indexed exactly by observation
IDs with one column per retained set. Existing keys require an explicit overwrite policy; silently replacing a score
matrix is prohibited.

No `random_seed`, group, Sample, Condition, Technical batch, p-value threshold, or network weight is exposed. AUCell
2.2 is unweighted and non-inferential. `raw=False`, `empty=False`, and `verbose=False` are fixed implementation
details. A shared internal resource parser is an internal seam, not another public input type.

## Invariants and reporting

- Input and output observation IDs/order are identical; the input is immutable.
- Feature IDs, set IDs, and targets are nonblank scalar strings; duplicate feature IDs fail; exact duplicate resource
  edges are collapsed and counted.
- Resource targets are intersected without case folding or online translation. Empty sets and sets below
  `min_targets` are disclosed; no surviving set is a hard error.
- Values are finite and scores lie in `[0,1]`; output shape and column identity equal the independently computed
  retained resource family.
- Resource bytes are streamed into SHA-256 and metadata are strict JSON. File content, not path/mtime alone, controls
  fingerprints and cache invalidation.
- `summary` follows the repository analysis-report schema and includes methods, report-ready key results, resource
  accounting, limitations, references, and dynamic versions for OpenBio, decoupler, Scanpy, AnnData, NumPy, pandas,
  and SciPy. JSON must reject NaN/Infinity.
- `code` defines an importable equivalent function returning `(AnnData, summary_dict)` and contains the same fixed
  method arguments, resource hash/metadata checks, postconditions, and summary builder.

## Scientific boundary

This node characterizes each observation. It does not compare Conditions, treat cells as replicates, adjust for
Technical batch, assign a population, or establish pathway activity as causal truth. Reports must say that any
Condition inference needs a separate Sample-level method within a defined population. If the input uses a Technical
batch-integrated representation elsewhere, that is irrelevant here because AUCell reads declared expression only.

## Migration and tests

Saved workflows require port migration: add metadata, remove seed, preserve an explicit source choice, and add
`summary`/`code`. The old default
`openbio-singlecell/gene_sets.csv` must not be silently retained unless it becomes an explicitly versioned,
licensed, fixed-hash asset. Do not auto-invent metadata.

Tests must cover the decoupler 2.2 signature and output key, rejection of 1.x/missing APIs, exact hand-ranked scores,
ties/feature-order disclosure, `n_up`, set pruning, malformed GMT/CSV, identifier collisions, resource SHA changes,
normalized and explicit Raw/count-like sources with cautious disclosure, sparse input, observation preservation,
overwrite behavior, memory/output
guards, input immutability, strict JSON, dynamic versions, and exact runtime/generated summary equality.

## Cohesion assessment

The module has one reason to change: AUCell's definition and observation-level result contract. Resource parsing,
backend drift, batching, validation, and reporting are hidden implementation complexity, creating leverage behind a
small interface. Statistical comparison stays downstream.
