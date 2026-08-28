# OpenBioSingleCellTasccodaDifferentialComposition: pre-change design review

## Decision

**Retain and substantially enhance this node as the hierarchy-aware tascCODA Sample-level differential-composition module. Do not merge it with flat scCODA behind a method switch, do not merge it with descriptive composition summaries or generic tests, and do not delete it.**

Its atomic scientific question is:

> For one explicitly directed two-condition comparison and one validated, pre-specified cell-type hierarchy, which hierarchy nodes have credible relative compositional effects under tascCODA, and what leaf-level compositional changes are implied by those selected node effects?

The hierarchy and aggregation bias change the estimand and selection structure. The primary discoveries are hierarchy nodes, while leaf effects are propagated consequences. A flat/hierarchical switch would make input requirements, output rows, tuning, dependency installation, migration, and interpretation conditional throughout one shallow interface.

Share only a private, well-tested composition adapter with the scCODA node: canonical Sample aggregation, Sample metadata validation, exact contrast construction, reference-candidate checks, version-gated NUTS execution, diagnostic normalization, software provenance, strict JSON, and atomic result assembly. tascCODA owns hierarchy validation, reference-path constraints, tree-adaptive penalty policy, node selection, and node-to-leaf reporting.

Repository constraints applied from `CONTEXT.md` and `docs/agents/domain.md` are: Sample means biological replicate, Condition contrast is Sample-level and replicate-aware, and Technical batch is not silently promoted to replication. ADR `docs/adr/0001-separate-marker-evidence-from-condition-inference.md` keeps marker evidence separate from this formal Condition inference. The local issue-tracker instructions were read for process context; this research-only task intentionally creates no issue or ledger entry.

## Why the current module is unsafe

- `covariate_keys` plus arbitrary `formula` does not define one comparison or its direction.
- `hierarchy_levels` has no public statement about whether it includes the leaf; upstream requires a root-to-leaf list.
- no invariant proves that every leaf maps to one path, every child to one parent, and the leaf set equals `annotation_key`.
- `phi` wrongly excludes negative values, hiding the part of the official parameter space that favors higher-level aggregation.
- `lambda_1=3.5` silently departs from the paper/Pertpy 1.3.0 default of `5`.
- `estimated_fdr` is scientifically inert for tascCODA's spike-and-slab LASSO, yet appears to promise controllable FDR.
- intercept, derived leaf effects, and selected node effects are concatenated into a heterogeneous table.
- the result omits hierarchy fingerprint, reference ancestors, Delta selection, `phi` sensitivity, sampler diagnostics, dependency versions, citations, limitations, and equivalent code.
- Pertpy's one-chain result can be misread as proven convergence.

These are not isolated UI defects. They distribute essential method knowledge across callers and make legacy serialized values impossible to interpret safely.

## Proposed public interface

Recommended outputs:

```text
table, summary, code
```

Recommended inputs:

| Input | Visibility | Contract |
|---|---|---|
| `adata` | primary | Cell-level AnnData; never mutated. |
| `sample_key` | primary | Canonical biological experimental-unit identifier. |
| `annotation_key` | primary | Exact cell-type leaf label. |
| `annotation_status` | primary | Required `curated` or `provisional`; provisional leaves force an explicitly exploratory result. |
| `hierarchy_keys` | primary | Ordered ancestor columns from root through immediate parent; implementation appends `annotation_key` as leaf. Nonempty. |
| `condition_key` | primary | Sample-constant Condition field. |
| `reference_condition` | primary | Exact baseline level. |
| `comparison_condition` | primary | Exact compared level; all directions are comparison relative to reference. |
| `reference_cell_type` | primary | Exact leaf or `automatic`; report resolved leaf and reference path. |
| `aggregation_bias` | primary or clearly result-defining advanced | Signed finite `phi`, default `0.0`; negative favors larger subtrees, positive favors leaves. |
| `adjustment_covariates` | advanced | Optional Sample-level nuisance fields, default empty, with deterministic, disclosed coding; no arbitrary transformations/interactions. |
| `num_samples` | advanced | NUTS posterior draws, default `10000`. |
| `num_warmup` | advanced | NUTS adaptation draws, default `1000`. |
| `random_seed` | advanced | Nonnegative seed, default `0` to match pinned Pertpy, for initialization/JAX sampling. |

Remove `estimated_fdr` entirely. It has no effect on the pinned tascCODA selection method. Remove public `lambda_1` unless a later, separate prior-sensitivity capability is designed; the current unexplained value is more dangerous than useful.

Condition levels must be explicit, not inferred from observation/category order. `hierarchy_keys` must not accept the leaf again; reject a list containing `annotation_key` with a message that the leaf is appended automatically. A structured hierarchy artifact could be considered later, but ordered `.obs` columns are sufficient for the present atomic node when strictly validated.

## Hidden, fixed, and disclosed policy

Hide as versioned adapter policy:

- Pertpy `type="cell_level"`, MuData modality names, `add_level_name=True`, and internal tree key;
- Patsy formula string and collision-safe internal design-column names;
- fixed zero pseudocount `0.5` and automatic-reference absence threshold `0.05`;
- fixed NUTS inference path and pinned kernel defaults;
- tascCODA `lambda_0=50`, `lambda_1=5`, and `theta=0.5` under Pertpy 1.3.0;
- ancestor-matrix construction, reference-node indices, and upstream frame layout;
- dependency-specific tree object types (`ete4`/`toytree`).

All fixed values remain visible in `summary` and `code`. They are hidden from routine tuning because changing them creates a prior-sensitivity analysis, not because they are unimportant.

Do not expose Pertpy's alternative `dendrogram_key` inside this inference node. A data-derived dendrogram and a curated lineage hierarchy have different provenance and validation questions. If dendrogram-derived hierarchy is needed, an upstream hierarchy-construction/review node should materialize the same explicit ancestor-column contract before tascCODA.

`aggregation_bias` cannot be hidden: it changes whether the prior favors high-level or leaf effects, and the original paper recommends sensitivity/cross-validation across values to avoid misspecification. The node fits one pre-specified value atomically. A future `tascCODA phi sensitivity` node could compare several values with held-out predictive scoring, but must be separate from this single-fit result.

## Deep execution seam

```text
validate AnnData/Sample/Condition/annotation metadata
  -> aggregate deterministic Sample x leaf counts
  -> validate ordered root-to-leaf paths and construct canonical tree
  -> fingerprint levels, nodes, leaves, edges, and ancestor matrix
  -> validate exact contrast, replication, adjustments, rank, and reference path
  -> run pinned Pertpy 1.3.0 tascCODA NUTS on temporary copies
  -> extract node-selection evidence, propagated leaf effects, and diagnostics
  -> prove node/leaf/ancestor axes and numeric postconditions
  -> emit one scoped canonical table, strict summary JSON, and equivalent source atomically
```

Tree validation belongs inside this module because the external library's permissive object construction is not the scientific boundary OpenBio promises. The internal adapter should return typed canonical data—count frame, Sample frame, tree manifest, exact design manifest—not a raw MuData that leaks Pertpy layout to callers.

## Sample, count, and design preconditions

The shared checks match the flat scCODA node:

- nonempty AnnData; unique observation identifiers; required `.obs` keys present and complete;
- each cell maps to one canonical Sample and one leaf label;
- `annotation_status` is explicit; provisional leaves are permitted only with `inference_scope="exploratory"` and provisional-label wording throughout the table and summary;
- every Sample maps to exactly one Condition, hierarchy-independent adjustment value, and biological experimental unit;
- nonnegative integer Sample-by-leaf counts and positive Sample totals;
- exactly the two selected Condition levels enter the model;
- reference/comparison are distinct and each has at least two biological Samples under OpenBio formal-contrast policy;
- technical replicates do not masquerade as biological replication;
- repeated/paired designs reject unless an explicitly supported fixed-effect policy is implemented; this node has no random effects;
- deterministic nuisance coding, finite design, exact focal coefficient, full rank, and reported condition number;
- at least two modeled leaves and a valid explicit/automatic reference.

Low Sample replication, uneven cells per Sample, sparse leaves, provisional annotation, and strong imbalance remain prominent limitations even when preflight passes. More cells within a Sample do not increase the number of independent experimental units.

## Hierarchy preconditions and manifest

`hierarchy_keys` describes ancestors only. The canonical upstream sequence is `[*hierarchy_keys, annotation_key]`. Preflight must establish:

1. at least one ancestor level and one leaf level;
2. unique, nonreserved column names and no repeated `annotation_key`;
3. no missing/empty hierarchy value on any modeled cell;
4. each leaf has exactly one complete ancestor path;
5. for every adjacent pair of levels, a child has exactly one parent;
6. one root after normalization, all nodes reachable, no cycles or orphan branches;
7. leaf labels exactly equal the modeled count-table columns, including ordering after canonicalization;
8. internal-node names remain unique after level prefixing and cannot collide with leaves/reference suffixes;
9. every internal node has descendants and every declared category is represented;
10. the reference leaf and each ancestor are in the computed reference-node set;
11. deterministic node/edge/leaf order and a stable hash of the hierarchy contract.

The manifest in `summary` should include ordered level keys, root, node/edge/leaf counts, branching-degree summary, ordered leaf list/hash, canonical edge-list hash, ancestor-matrix shape/hash, and reference path. It should also record whether hierarchy labels came from curated knowledge or data-driven clustering; uncertainty in that hierarchy is a scientific limitation.

A leaf count at or above the original paper's typical `p < 100` regime should trigger a scope/resource warning and disclose expected computational risk. A hard ceiling should be based on measured OpenBio resource tests rather than invented from the paper.

## `aggregation_bias` and sensitivity policy

Validate `aggregation_bias` as finite and signed; do not impose the current `min=0`. No official hard magnitude bound was found, so the implementation should not invent one silently. It must compute the resulting node-wise penalty scales before inference, reject nonfinite values, and warn when extreme `phi` saturates the scaling so strongly that the effective prior is nearly binary. The summary explains:

- negative: less penalty on nodes representing larger subtrees, favoring aggregation nearer the root;
- zero: equal node penalty under the pinned scaling;
- positive: increasing preference for detailed/leaf effects.

This direction must be covered by exact tests against Pertpy's scaling formula. The report records the scaled `lambda_1` range actually produced across nonreference nodes, because the same scalar `phi` interacts with tree size.

For confirmatory use, `phi` should be pre-specified from biological expectations or selected in a separate held-out/CV workflow. If the user reports one exploratory fit, mark it `phi_sensitivity_not_assessed`. The node itself must not fit multiple values, select the one with the most discoveries, or claim false-discovery control.

Keep `lambda_0=50`, `lambda_1=5`, and `theta=0.5` fixed and versioned. The output explicitly notes that Pertpy 1.3.0 fixes `theta`, rather than sampling the Beta prior in the original paper, following the documented HMC/NUTS collapse defect. A dependency below 1.1.1 is a hard error; exact 1.3.0 pin/testing is recommended.

## Canonical `table` contract

One table can carry both necessary scopes only when each row is unambiguous. Recommended stable union:

| Field | Node row | Leaf row |
|---|---|---|
| `effect_scope` | `hierarchy_node` | `derived_leaf` |
| `contrast` and Condition fields | exact focal comparison | same |
| `reference_cell_type` | resolved leaf | same |
| `effect_name` | canonical prefixed node | exact leaf label |
| `hierarchy_level` | node level | leaf level / annotation key |
| `descendant_leaf_count` | positive integer | `1` |
| `descendant_leaves` | deterministic JSON array/string artifact | the leaf itself |
| `model_effect` | selected node final parameter | propagated sum of selected ancestor/node effects |
| `posterior_median` | node median | leaf beta median only if upstream meaning is validated; otherwise null |
| `hdi_lower` / `hdi_upper` | node HDI | leaf HDI |
| `posterior_sd` | node SD | leaf SD |
| `selection_delta` | node-specific Delta | null |
| `credible_effect` | direct tascCODA node decision | null, not inherited as an independent discovery |
| `selection_basis` | `sslasso_delta` | `ancestor_node_sum` |
| `expected_count_reference` / `expected_count_comparison` | null | derived compositional expectations if exactly recoverable |
| `compositional_log2_fold_change` | null | relative renormalized leaf change |
| `reference_constraint` | true for reference nodes | true for reference leaf |

Primary ordering should place credible `hierarchy_node` rows first in deterministic hierarchy order, then noncredible node rows if the complete model table is retained, then `derived_leaf` rows. Alternatively, return only credible node rows plus all derived leaves, but the omission policy must be explicit. Never repeat one selected ancestral effect as multiple “credible cell-type discoveries.”

Intercept and nuisance-covariate frames do not belong in the primary table. Store only JSON-safe model details when needed. If two separate table outputs are later preferred, that changes the required output contract; this audit honors the requested single `table` endpoint and makes scope explicit within it.

## `summary` JSON contract

`summary` is strict JSON with no NaN/Infinity, NumPy/pandas scalars, tuple keys, or tree/backend objects. It should contain:

- `schema_version`, node/method identifiers, run status, and UTC completion metadata;
- `question`: Sample unit, exact Condition direction, adjustments, annotation key, and hierarchy-aware estimand;
- `input`: cells, Samples, group sizes, leaves, Sample-count distribution, zeros, pseudocount, annotation status/provenance;
- `design`: exact generated columns/focal term, categorical bases, dimensions, rank, condition number, and validation results;
- `hierarchy`: ordered keys, root, node/edge/leaf counts, branching summary, fingerprints, reference path, lineage source/status, and validation manifest;
- `model`: tascCODA/Dirichlet-multinomial, signed `aggregation_bias`, scaled penalty range, fixed `lambda_0=50`, `lambda_1=5`, `theta=0.5`, NUTS draws/warm-up/seed/backend;
- `selection`: `sslasso_delta`, explicit statement that no adjustable FDR was applied, credible-node count, and node-specific Delta disclosure;
- `diagnostics`: chain count, single-chain limitation, mean acceptance, within-chain ESS/MCSE where valid, potential-energy/steps/step size, finite checks, and `divergences_available=false`;
- `results`: ordered writing-ready credible node records with level/descendants/effect/HDI/Delta, plus separately labeled derived leaf consequences;
- `sensitivity`: whether `phi` was pre-specified, whether held-out/CV sensitivity was assessed elsewhere, and the original-paper misspecification caveat;
- `warnings`/`limitations`: relative-only inference, low replication, annotation/hierarchy uncertainty, moderate-leaf regime, no random effects, single chain, and no FDR-control claim;
- `software`: exact OpenBio/Pertpy/JAX/NumPyro/ArviZ/MuData/AnnData/Patsy/pandas/NumPy/ete4/toytree/Python versions and device;
- `references`: tascCODA, scCODA, NUTS, Pertpy, compositional-data, and official API/source citations.

Writing prose should take the form: “Under tascCODA with B relative to A, reference leaf R, and aggregation bias φ, K hierarchy nodes were credible under the pinned spike-and-slab LASSO/Delta rule; node N represents descendants [...].” It must not say “FDR-adjusted significant”, assert absolute population changes, call every descendant independently credible, or assert convergence from one chain.

If diagnostics are warning/failed, key-result prose is qualified or suppressed. The JSON may contain model evidence on warning but must never return a success status with nonfinite values.

## Equivalent `code` contract

The `code` endpoint defines one importable function that reproduces the same table and summary. It must:

1. validate Sample, contrast, adjustment, count, annotation, and complete hierarchy invariants;
2. append the leaf key itself and construct/fingerprint the same canonical root-to-leaf tree;
3. build the same collision-safe exact comparison and verify rank;
4. gate the same Pertpy 1.3.0 dependency contract and copy input data;
5. call `Tasccoda.load`, `prepare` with explicit `phi/lambda_0/lambda_1/theta/tree_key`, and NUTS with explicit draws/warm-up/seed;
6. avoid `estimated_fdr` entirely;
7. validate node/leaf/reference axes and produce the same scoped table;
8. extract the same available diagnostics and strict JSON with versions/references/limitations.

An illustrative upstream snippet alone is not equivalent code. Runtime and code output must reject the same malformed hierarchy or mixed Sample metadata and must not partially mutate the input or return only the leaf frame.

## Diagnostic and failure policy

Use the same stable statuses as the flat node:

- `limited_single_chain`: finite fit with Pertpy mean acceptance inside `[0.6, 0.95]` and no declared within-chain warning, while explicitly not claiming between-chain convergence;
- `warning`: fit returned but acceptance/ESS/MCSE/energy/tree-step/resource/replication/annotation/hierarchy/phi-sensitivity concern exists;
- `failed`: invalid inputs, nonfinite posterior, malformed node/leaf/reference axes, missing focal effect, or result-bundle serialization failure.

R-hat is undefined for Pertpy's one-chain default. Divergence count is unavailable because Pertpy 1.3.0 does not retain the field. Do not convert either into a reassuring value. If ArviZ is used, test chain/draw dimensions and only report diagnostics it can validly compute.

Hierarchy-specific post-fit checks include:

- upstream node names exactly match the canonical manifest after prefixing;
- ancestor-matrix dimensions are leaves by nodes and its reference constraints match the reference path;
- every derived leaf effect equals the selected node-effect sum indicated by the ancestor matrix within numeric tolerance;
- node `Delta`, credibility, final parameter, and median relations follow the pinned rule;
- no FDR/inclusion-probability field is fabricated;
- all HDI/SD/effect values are finite or explicitly nullable under a documented reference constraint.

All outputs are staged and committed together after checks. Errors name the hierarchy level/path/node or Sample/design invariant that failed and include the installed Pertpy version.

## Workflow migration boundary

Legacy tascCODA nodes must be rejected atomically rather than automatically rewritten:

- old `formula`/`covariate_keys` do not encode exact Condition levels/direction;
- the old schema does not record curated versus provisional annotation status;
- old `hierarchy_levels` may include or omit the leaf, and no stored manifest proves its contract;
- old `phi` could never be negative due to the widget constraint;
- old `lambda_1=3.5` is not equivalent to the proposed fixed `5`;
- old `estimated_fdr` is inert and cannot map to any new selection control;
- old table consumers may assume concatenated intercept/leaf/node rows, while the new stable table has explicit scopes;
- `summary` and `code` outputs are added.

The migration error should instruct users to recreate the node with explicit reference/comparison Conditions, ancestor-only hierarchy keys, reviewed `aggregation_bias`, and downstream table-scope review. Connected or exposed legacy formula, hierarchy, `lambda_1`, `estimated_fdr`, or reordered widgets make rejection mandatory.

Migration preflight covers the whole graph before mutation. Current schema is no-op; mixed old/new inputs, extra legacy keys, divergent named/positional values, partial output bundles, malformed endpoints/links, and proxy/direct exposure conflicts reject. Root/subgraphs, object/array links, and idempotence follow repository migration policy during the later implementation phase.

## Verification plan

### Schema and contract tests

- exact input/output names/order/defaults/visibility and registration;
- absence of `estimated_fdr` and `lambda_1`; signed finite `aggregation_bias` accepts negative/zero/positive values;
- hierarchy list excludes the leaf and appends `annotation_key` exactly once;
- strict JSON, complete method/software references, and code compilation/parity;
- caller AnnData is unchanged after success/failure.

### Sample/design tests

- deterministic Sample-by-leaf counts for dense/sparse data;
- Sample-constant Condition/adjustments, missing values, extra levels, minimum replication, technical-unit and repeated-design errors;
- exact safe Condition contrast, categorical nuisance bases, rank/condition-number/confounding checks;
- explicit/automatic reference leaf, 0.05 candidate threshold, no candidate, and complete reference path.

### Hierarchy tests

- valid bifurcating and multifurcating examples;
- missing paths, leaf-to-multiple-path, child-to-multiple-parent, multiple roots, orphan/unused category, empty internal node, duplicate/colliding names, wrong order, repeated leaf key, leaf-set mismatch, and reference-path corruption;
- deterministic edge/ancestor-matrix fingerprint under stable category order;
- exact negative/zero/positive `phi` scaling direction and scaled penalty disclosure;
- extreme finite `phi` saturation warning and nonfinite scaled-penalty rejection;
- warning at the documented moderate-feature regime without inventing an unsupported hard cutoff.

### Adapter/output/diagnostic tests

- deterministic fake Pertpy verifies exact `load` levels, `add_level_name`, tree key, fixed penalties, NUTS arguments, and absence of FDR calls;
- node Delta/median/final/credibility rule and node-to-leaf propagation equality;
- stable union table with nullable scope-specific fields; leaf rows never receive an independent credibility flag;
- coefficient versus renormalized compositional-LFC semantics;
- one-chain limitation, undefined R-hat, available ESS/MCSE, acceptance warning, unavailable divergences, nonfinite posterior, malformed axes, and theta/version gate;
- dependency error identifies `pertpy[tcoda]` and rejects versions earlier than the pinned contract;
- short version-pinned CPU smoke fit, optional/skipped when extras are absent.

### Migration and workflow tests

- every legacy schema rejects without changing any root/subgraph byte;
- formula/hierarchy/phi/lambda/FDR exposures or connections, proxy widgets, partial outputs, malformed object/array links, and mixed schemas reject in full-graph preflight;
- current schema produces zero migrations and repeated runs are idempotent;
- example workflow supplies curated hierarchy ancestry, explicit Sample-level Condition direction, and consumes `table`, `summary`, and `code` with node-versus-leaf scope awareness.

## Acceptance criteria

The refactor is complete only when:

1. one run is one explicit Sample-level two-condition hierarchy-aware comparison;
2. the hierarchy contract and fingerprint are proven before inference;
3. negative, zero, and positive aggregation bias have correct, disclosed semantics;
4. fixed prior values match pinned Pertpy 1.3.0 and the broken historical theta behavior is version-gated;
5. no inert FDR control remains anywhere in the schema, code, table, or prose;
6. direct credible node effects and derived leaf consequences cannot be confused;
7. one-chain and `phi`-sensitivity limitations are impossible to hide;
8. strict `summary` and equivalent `code` disclose every result-defining input, hidden policy, version, reference, diagnostic, and limitation;
9. flat scCODA remains a separate public node while only the validated internal composition adapter is shared;
10. ambiguous legacy state is rejected atomically rather than guessed.

The official and primary-source evidence behind these choices is recorded in [`official-usage.md`](./official-usage.md), including the root-to-leaf loader contract, tree-adaptive prior, signed `phi`, fixed penalty defaults, inactive FDR input, selection fields, dependency extras, and one-chain diagnostics.
