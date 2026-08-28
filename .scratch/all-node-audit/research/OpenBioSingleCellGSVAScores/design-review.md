# OpenBioSingleCellGSVAScores — design review

## Final open-expert boundary

State provenance is report context, never execution authority. For Poisson, value-domain compatibility is enforced
independently of the chosen X/Raw/layer container; other kernels warn for nonrecommended scale. This supersedes any
stricter source-container or history language below.

## Decision

**Enhance and retain as one atomic GSVA observation-scoring module.** Migrate to decoupler 2.2 `mt.gsva`, restore
official defaults, expose kernel semantics, and add a mandatory dense-memory guard. Do not merge with AUCell: GSVA
is cohort-dependent, density based, signed, and materially more memory intensive.

## Atomic interface

Inputs, in order:

1. `adata`;
2. local `gene_sets_file`;
3. strict `resource_metadata_json`;
4. explicit expression `source`;
5. advanced `source_column`, `target_column`, and `min_targets`;
6. `kernel` with `gaussian_normalized`, `poisson_counts`, or `empirical` choices;
7. advanced `maxdiff=True`, `absrnk=False`, and positive finite `tau=1.0`;
8. `max_working_memory_gib` and bounded `batch_size`;
9. nonblank `output_key`, default `gsva_scores`, plus explicit overwrite policy.

Outputs:

```text
adata, summary, code
```

Kernel and data state are coupled in reporting rather than by an unnecessary authority gate: Gaussian and empirical
recommend finite normalized values, while an explicit expert-selected count-like/Raw source remains executable with
a prominent warning. Poisson requires finite nonnegative integer values, independent of whether they are explicitly
selected from X, Raw, or a named layer; only violating that arithmetic kernel domain remains a hard error.
`random_seed` is removed because decoupler 2.2 GSVA is
deterministic. `raw=False`, `empty=False`, and `verbose=False` are hidden.

## Output and provenance contract

Store a named DataFrame at `obsm[output_key]` with exact observation order and independently validated retained set
columns. Validate finite values and the documented `[-1,1]` method range within numeric tolerance. The input is
immutable. Resource parsing/fingerprinting follows one shared internal module and records exact duplicate collapse,
universe intersection, per-set sizes, exclusions, path and streamed SHA-256.

`summary` must include report-ready methods/results, score distribution and extreme-set previews, cohort dependence,
kernel/state rationale, dense-byte estimate versus guard, set/target accounting, warnings/limitations, references,
and dynamic versions. It must be strict JSON. `code` must define the equivalent function returning
`(AnnData, summary_dict)`, with identical resource validation, memory calculation, decoupler call, postconditions,
and scientific summary.

## Scientific boundary

No Condition, Sample, population, or Technical batch parameter belongs to GSVA scoring. The score module does not
perform inference. A downstream contrast must treat Sample as the replicate and Technical batch as nuisance; a
cell-level t-test is explicitly invalid. GSVA scores also change with the observation cohort, so the report must pin
ordered observation identity/fingerprint and warn against comparing scores computed in separate runs.

## Migration and verification

Migration removes `random_seed`, renames old `mx_diff` to the implementation-accurate `maxdiff`, changes the unsafe
`abs_rnk=True` default to `False`, adds kernel/resource metadata/memory/overwrite inputs, and adds `summary`/`code`.
Unversioned default resources cannot be auto-migrated.

Tests: 2.2 signature/output key; 1.x/missing API failure; Gaussian/empirical expert count-state warnings and strict
Poisson state validation; official
default regression fixture; signed scores and column identity; cohort/order pinning; set pruning; malformed resources
and changing bytes; sparse/dense memory estimates and guard edges; zero/empty/non-finite observations; overwrite and
immutability; strict JSON; compiled generated code; and exact runtime/generated summary equivalence.

## Cohesion assessment

The interface exposes only choices that change the named GSVA method or resource/state interpretation. Backend
batching, densification, resource parsing, validation, and reporting stay inside the module, concentrating the
complexity that would otherwise leak to every caller.
