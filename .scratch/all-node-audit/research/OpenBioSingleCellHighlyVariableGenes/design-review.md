# Highly Variable Genes: module design review

## Current module

The node performs one feature-selection analysis, dispatching standard flavors to `scanpy.pp.highly_variable_genes` and the Pearson flavor to `scanpy.experimental.pp.highly_variable_genes`. It copies the input, optionally marks forced genes, and optionally subsets the variable axis.

## Correctness and interface problems

1. All flavors share one source default even though `seurat`/`cell_ranger` require normalized logarithmized expression and `seurat_v3`/`seurat_v3_paper`/`pearson_residuals` require counts. Switching flavor can silently apply the method to the wrong state.
2. Count flavors rely on Scanpy's warning-only integer check; normalized or log-transformed input can pass and produce a result.
3. Pearson `theta`, clipping, and `chunksize`, Seurat-v3 `span`, and dispersion `n_bins` are not representable or disclosed. In particular, the Pearson clipping rule is a result-defining model choice.
4. `always_keep_genes` changes the scientific selection but the result does not distinguish algorithm-selected, forced-added, and final counts. A final set larger than `n_top_genes` is therefore unexplained.
5. `subset=True` can remove full-gene matrices without requiring the repository's full-gene Raw snapshot invariant.
6. Missing/null grouping labels, one-cell groups, a grouping column with only one level, `n_top_genes` larger than eligible features, prior HVG annotations, backed mode, and non-unique axes are not handled explicitly.
7. The fixed Scanpy `.var` fields are silently replaced on repeated runs, and there is no overwrite policy.
8. The node exposes neither a structured `summary` nor equivalent `code`.

## Decision: keep, unify dispatch, and enhance

Keep one node for all flavors. Although their statistics differ, they answer the same atomic question and produce the same feature-selection contract. A single deep module centralizes flavor-to-expression-state validation, batch-aware ranking disclosure, forced-gene policy, subsetting safeguards, and reporting. Splitting by flavor would duplicate those rules and invite drift.

Do not merge HVG selection with normalization, Pearson-residual transformation, PCA, or integration. Those operations consume different expression states and have independently meaningful results. In particular, chunked Pearson HVG selection should be usable before the dense Pearson transform, and Scanpy PCA can consume the `highly_variable` mask without physical subsetting.

## Target interface

- Visible inputs: `adata`, flavor, flavor-appropriate expression source, `n_top_genes`, and optional grouping key.
- Advanced scientific inputs: forced-gene list; flavor-conditional `theta` and clipping for Pearson residuals; `span` for v3; `n_bins` for dispersion flavors.
- Advanced performance/state inputs: Pearson `chunksize`, `subset`, and overwrite-existing opt-in.
- Outputs: annotated or explicitly subset copied `adata`, structured `summary`, equivalent Python `code`.

The flavor control should own conditional subinputs so the suggested source changes with the method: normalized/logarithmized layer for `seurat` and `cell_ranger`, count layer for v3 and Pearson residuals. The source remains visible and validated; the node must not silently reinterpret one matrix based on its name. Stable defaults (`n_bins=20`, `span=0.3`, `theta=100`, automatic `sqrt(n_obs)` clipping, `chunksize=1000`) may be advanced, but they must be explicit in generated code and summary.

`check_values`, `inplace`, Scanpy function dispatch, cutoff parameters ignored by rank-based `n_top_genes`, and method-specific output-column routing are hidden implementation policy. `filter_unexpressed_genes` should also be fixed according to documented Scanpy behavior rather than exposed as a second feature-filtering decision.

## Validation and scientific contract

- Require an in-memory object with non-empty axes, unique observation/variable identifiers, and a finite numeric selected matrix. Enforce the minimum observation/group size only for flavors whose backend actually fails on the supplied shape; Pearson-residual HVG is defined for singleton inputs/groups in Scanpy 1.12.3.
- For `seurat` and `cell_ranger`, non-negative normalized/logarithmized expression is the documented expectation. Proven/unknown state mismatches are warnings because provenance is evidence, not proof; numeric conditions that make the backend non-finite may still fail.
- For v3 and Pearson flavors, require finite non-negative values and positive totals where the count statistic divides by them. Fractional values and provenance inconsistent with counts proceed with warnings and cannot be reported as verified raw counts.
- Require positive `n_top_genes`, but allow it to exceed `n_vars` or the expressed-feature count. Report requested versus algorithm-selected counts and captured backend warnings rather than pre-empting Scanpy's documented “return all eligible genes” behavior.
- If grouped, require a present column with no missing/blank labels. Apply a per-group minimum only when required by the selected flavor backend; otherwise warn for singleton groups. Report whether the column represents Sample or Technical batch; never infer that role from its name. Warn when only one level exists because grouping then has no effect.
- Resolve forced genes against stable `var_names`; deduplicate the request, report missing identifiers, add only present genes after algorithmic selection, and preserve a separate boolean/summary count for forced additions. Never claim they were selected by the method.
- If `subset=True`, audit the full-gene Raw snapshot when present and warn when absent/misaligned/incomplete; subsetting itself remains executable. Preserve valid Raw while slicing current `X`, layers, and variable annotations.
- Fail when prior HVG result columns exist unless overwrite was explicitly requested. Remove stale flavor-specific HVG columns before a permitted rerun so metrics from different methods cannot coexist misleadingly.
- Translate a missing `scikit-misc` dependency for v3 flavors into an actionable installation error.

## Report contract

`summary` includes:

- methods text naming the selected flavor, exact public Scanpy function, expected source state, grouping/ranking rule, and any forced-gene override;
- results text stating the algorithm-selected and final feature counts without presenting HVGs as markers or condition evidence;
- key results for input/output dimensions, requested/actual/forced counts, selection rate, missing forced genes, flavor-specific score distributions and top-ranked feature identifiers, group sizes, `highly_variable_nbatches` distribution/intersection, and Raw preservation when subset;
- all resolved flavor-specific and hidden parameters, warnings and limitations;
- only the applicable method citation(s), plus Scanpy/AnnData references and dynamic software versions, including `scikit-misc` for a v3 run.

`code` is a self-contained function with resolved parameters. It copies the input, validates the same source-state/group/count/overwrite/Raw invariants, calls the correct public Scanpy function with explicit arguments and `subset=False`, marks forced genes separately, optionally slices only after preserving the Raw invariant, and returns an equivalent `AnnData`. It does not recreate plugin history.

## Cohesion and coupling assessment

This is a deep in-process module when flavor dispatch and validation stay behind one interface. Callers learn one feature-selection result, while the implementation owns five method variants, two source-state families, batch merging, optional dependencies, forced-gene attribution, subsetting safety, citations, and version capture. Normalization, scaling, integration, and downstream PCA remain on separate seams because deleting this node would otherwise spread method-selection and provenance logic across all of them.

## Verification plan

- Each flavor is compared with a direct Scanpy call on the declared source; provenance mismatches and fractional count-like values are limited/disclosed while hard numeric/backend failures remain errors.
- `seurat_v3` and `seurat_v3_paper` agree without grouping and exercise their documented different group-ranking priorities with grouping.
- Pearson `theta`, automatic/custom/no clipping, and different chunk sizes preserve selections as expected and report resolved settings.
- Dense/CSR/CSC sources, missing layers, duplicate axes, empty axes, zero totals, non-finite data, positive oversized `n_top_genes`, and missing optional dependency follow the contract.
- Grouped selection covers missing labels, too-small groups, one-level warning, recurrence distribution, and Sample/Technical-batch wording.
- Forced genes cover duplicates, missing identifiers, already-selected genes, added genes, and final counts greater than `n_top_genes`.
- `subset=True` warns on absent/misaligned Raw, preserves a valid full-gene Raw snapshot, and slices all current aligned matrices consistently.
- Repeated runs honor overwrite policy, remove stale metrics, return strict-JSON summary, and generated code reproduces annotations, axis shape, and primary values.

## Final audit amendment

For the default `seurat` flavor, the normalized-log layer suggestion must be the selected interface default. Legacy direct Python calls that omit the dynamic source may fall back to `X` only when the suggested layer is absent; normal UI calls carry their explicit selection. A standard external Scanpy object with `uns["log1p"]` is acceptable log-transformation evidence even when OpenBio cannot prove the preceding normalization; this uncertainty is a report limitation, not a reason to reject an otherwise valid official workflow.

For grouped Pearson-residual HVG selection, `sqrt_n_obs` is resolved once from the full input and passed to Scanpy as an explicit numeric clip. This keeps the reported parameter, generated code, and computation identical instead of inheriting Scanpy 1.12.3's internal first-group resolution of `clip=None`. Grouped reports also disclose the all-group intersection count/rate, and `top_algorithm_features` is restricted to the actual algorithm mask even when forced genes enlarge the final set.
