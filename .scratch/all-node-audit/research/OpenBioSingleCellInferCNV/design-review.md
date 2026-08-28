# Infer CNV: module design review

## Decision

**Retain and deeply refactor `OpenBioSingleCellInferCNV` as one expression-to-smoothed-CNV operation.** Keep representation, neighbors, clustering, UMAP and group scoring downstream. Introduce a typed `OPENBIO_CNV_STATE` because generic AnnData keys cannot prove expression source, reference, coordinates or output identity.

## Target interface

Inputs:

- `adata`;
- typed `source`, explicit current-axis expression selector defaulting to `log1p_norm`; normalized/log-transformed full-gene expression is the documented recommendation, not a history-backed runtime gate;
- `reference_key`, visible;
- `reference_categories`, visible and required nonblank exact labels;
- `sample_key="sample"`, visible, naming the biological-replicate column used to disclose reference support;
- `genome_assembly`, visible and required, recording the assembly matched by the coordinate columns;
- `window_size=100`, visible;
- `step=10`, advanced;
- `lfc_clip=3.0`, `dynamic_threshold=1.5`, advanced;
- `exclude_chromosomes="chrX,chrY,chrM"`, advanced, making infercnvpy 0.6.1's otherwise implicit mitochondrial exclusion explicit;
- `output_key="cnv"`, advanced;
- `minimum_reference_cells=20`, advanced;
- `chunksize=5000`, `n_jobs=1`, `max_output_gib=4.0`, advanced;
- `overwrite_existing=False`, advanced.

Fix `calculate_gene_values=False`; that high-memory output requires a separate deliberate capability. Do not offer all-cell reference or direct untyped arrays in this node.

Outputs: `cnv_state` (`OPENBIO_CNV_STATE`), `summary`, `code`. The artifact owns a defensive AnnData copy and records stage `inferred`, obs/var/source/coordinate/reference/output fingerprints, exact coordinate assembly metadata, parameters, retained chromosome/window mapping and infercnvpy version. Portable code returns the equivalent standalone immutable `(cnv_state, summary_dict)` artifact protocol without importing OpenBio helpers.

Preflight label integrity before conversion, one-to-one selected-source/current-feature alignment, coordinate values/order, reference presence/count/Sample coverage, finite real numeric values, output collisions and memory. Expression-state evidence and apparent count/scaled inputs produce warnings only; mutable analysis history and a Raw/current binding are never execution credentials. Feature-axis completeness is user-declared and reported as unverified. Postflight requires a finite numeric sparse/dense `n_obs × n_windows` `X_<key>`, valid window metadata, unchanged input expression, and exact artifact fingerprints.

`summary` includes reference counts by category and Sample, source/coordinate provenance, gene/chromosome/window accounting, CNV value/absolute-value quantiles, zero fraction, parameters, experimental warnings, citations and infercnvpy/Scanpy/AnnData/NumPy/SciPy/pandas/OpenBio versions.

## Migration and tests

Old empty reference defaults become a blocking migration prompt, not all-cell reference. Old `window_size=250` is preserved for existing workflows with a change note; new runs use official 100. A missing explicit transformed source maps only when state provenance proves the canonical layer.

Tests cover dense/sparse selected sources and cautious disclosure for nonstandard expression; typed label collisions/missing categories/too few references; reference distribution; coordinate nulls/types/order/assembly mismatch/duplicates; sex chromosome exclusion; memory/collision handling; exact fake backend call; malformed outputs; artifact ownership/fingerprints; real 0.6.1 smoke; strict JSON, compilation and runtime/code parity. They do not delete analysis history or assert a Raw/current provenance binding.

The implemented workflow migration is deliberately narrower than the runtime interface. A legacy node is eligible
only when it carries the exact, user-authored
`properties.openbio_cnv_migration={schema_version:1, source, sample_key, genome_assembly}` record. `source` is
exactly `{kind:"X"}` or `{kind:"layer",layer_name:"..."}`; Raw, missing/extra fields, blank or noncanonical
Sample/assembly/layer values, and a pre-existing review record block the entire workflow before mutation. This hint
records a migration choice but does **not** prove normalized log-transformed full-gene state; the migrated runtime
still reports the source choice and any available state evidence without treating history or Raw as proof.

`reference_key`, `reference_categories`, and `window_size` come only from the legacy serialization, never from the
hint. Empty reference labels block migration; legacy comma trimming/deduplication is canonicalized using the old
parser's semantics. The historical `window_size=250` is retained (and any eligible legacy window must be at least
10 so the fixed historical `step=10` remains valid). On success the executable hint is deleted and replaced by the
informational `openbio_singlecell_migration_review` record, including source/Sample/assembly, preserved window, and
the mandatory runtime-provenance warning. The output becomes `OPENBIO_CNV_STATE + summary + code`. Exact current
schemas are no-ops; mixed, partial, unconnected, rerouted, boundary-crossing, ambiguous, or non-paired legacy graphs
fail atomically.
