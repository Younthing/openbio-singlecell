# Expert validation audit

Reviewed against `1b8d52b208c9c52aa8c84f5f21c607b8af68ab6d`, following the request in `spec.md`.

## Method

Inventory runtime errors and UI bounds across the active node families; trace shared validators and generated code.
Compare questionable restrictions with installed backend implementations and small native executions. Add a failing
behavior test before each correction at the existing operation/node seam and verify equivalent generated code.
This is an audit of additional wrapper restrictions, not a claim to expose every option of every backend API.

## Corrections

| Family | Extra restriction removed | Evidence and scope |
| --- | --- | --- |
| Data | Input provenance equality; annotation-record equality; empty Raw axes | Native AnnData materializes empty axes. Merge follows explicit value-conflict policy, records both annotation histories for mixed columns, and preserves metadata when no values change. Missing Raw and ambiguous cell identifiers still fail. |
| QC/filtering | Mitochondrial percentage restricted to 0–100; contradictory bounds rejected | Signed expression produces percentages outside that range; explicit comparisons remain defined. Contradictory bounds produce an empty selection consistently in runtime and generated code. |
| HVG/Pearson/Scrublet | Signed/fractional/zero-total recommendations promoted to universal count gates | Native Scanpy 1.12.3 probes produce valid outputs for the retained expert examples. Individual undefined arithmetic or invalid backend results still fail. |
| Expression destinations | Reserved `counts` output name; selected-layer overwrite rejected despite explicit overwrite | Named destinations and explicit overwrite are honored, with truthful source-preservation reporting. Existing destinations remain protected unless overwrite is enabled. |
| Marker ranking/filtering | Extra tie-correction gate, excessive requested gene-count rejection, pre-categorical-only grouping, arbitrary filter-cutoff ranges | Native Scanpy ignores irrelevant tie correction, caps gene count, and converts repeated string labels. Filtering applies the supplied finite comparisons; actual probability values and paired evidence remain validated. |
| Embedding/integration | Cosine zero-row ban; scVI zero-count and role-reuse bans | Native Scanpy Neighbors/t-SNE and real scVI CPU training support the tested cases. Input axes, finite values, and actual output contracts remain checked. |
| PAGA/DPT | 256-group ceiling, pre-categorical-only PAGA groups, mandatory OpenBio Diffusion Map provenance, connected-only DPT, at least two diffusion components | Native Scanpy supports 257 groups, repeated string groups, ordinary stored diffusion maps, one component, and disconnected DPT. DPT preserves the native `+inf` unreachable-cell sentinel; summaries use reachable cells. |
| CellTypist | Exact CP10K total, integer-only declared counts, zero-total cells, overly tight log bound, prior ownership required for explicitly named overwrite | Installed backend warnings and accepted range are respected. The chosen transformation stays explicit. Ownership checks remain for automatic deletion of stale outputs. |
| Scoring/correlation | Zero-expression rows; incorrect requirement for two dimensions varying across cells | Native Scanpy panel scoring and Decoupler support the tested zero-row inputs. All-zero empirical GSVA is allowed; unsupported Gaussian behavior is not replaced with another kernel. Correlation validates actual centroid results. |
| Resource metadata | Mandatory complete release/citation metadata; fixed DGIdb name/date/URL/species labels | Metadata does not determine numerical execution. Optional metadata and custom string fields are retained, absent documentation disclosed, and actual content fingerprints used for references. Declared identity mismatches remain distinct from missing declarations. |
| Pathway contrast | Reused metadata columns rejected before computation | Role aliases are recorded, while the actual Welch design and degrees of freedom remain validated. |
| Augur | Signed-expression and seed-zero bans; count-like data coerced to integers | Selected source values and dtype are retained. Count-likeness describes the input and never authorizes changing it. Source fingerprints also distinguish fractional values. |
| Composition/Milo | Universal minimum per-arm sample counts; Bayesian paired-name/rank/residual-degree rules | Installed Pertpy separates Bayesian composition from edgeR design requirements. Low support and confounding are interpretation advisories; Milo keeps actual edgeR estimability requirements. |
| CNV/Schist/LIANA | Conservative chromosome/window/reference/group, sampling, and per-Sample identity thresholds | Installed backend implementations support the added examples. Schist seed zero remains rejected because the pinned backend actually fails; LIANA rank aggregation retains its real group-versus-rest requirement. |
| Velocity/SCENIC | Fractional/zero-library abundance, small selected gene sets, seed-zero preferences | Native scVelo conventions are retained and effective per-layer normalization targets are reported, including the zero-median target of one. SCENIC seed-zero behavior is disclosed. |
| Generic plots | Zero-containing log violin, two-group dendrograms, 50-gene aesthetic caps, infinite observation colors, duplicate/overflow palette colors | Native Scanpy renders the supported selections. Infinite color values are masked without changing stored values; explicit palette mappings remain reported. User-selected expression and panel sizes are preserved within actual resource limits. |
| Diagnostic plots | HVG history gene-count gate and required OpenBio HVG origins; mandatory OpenBio diffusion provenance; posterior point required inside HDI; disconnected/one-component/one-bin DPT trends | Native feature slices and standard stored Scanpy evidence work. Unknown HVG selection origins are reported as unknown; posterior points and intervals are rendered independently; DPT trends bin reachable cells only. |

## Retained boundaries and limits

- File/path safety, immutable artifact ownership, explicit overwrite protection, ambiguous identifiers, shape/axis
  alignment, and paired evidence/content fingerprints remain enforced.
- Mathematical/backend constraints remain local to the selected operation. Examples include undefined finite
  correlation/variance results, missing DPT root, single-group native PAGA failure, real LOESS span constraints,
  unsupported GSVA kernels, and integer pseudobulk/DE artifact contracts.
- MAD accepts zero or larger finite nonnegative multipliers without an arbitrary UI maximum; Scale accepts native
  negative clipping bounds. These choices are applied as supplied, not substituted with recommended defaults.
- Current graph pipelines still require graphs compatible with their undirected, zero-diagonal modularity/diagnostic
  contracts. This change does not add directed/self-loop algorithm variants.
- Full companion evidence plots (for example CellTypist, Schist, cNMF and correction diagnostics) still consume their
  defined evidence bundles. Supporting arbitrary external equivalents would require an import contract, not simply
  deleting an identity check.
- Optional backend families were checked by installed source and existing public-seam tests. Native smoke tests were
  used where available; a stubbed backend test alone is not treated as proof of native scientific support.
- Input adapters, save/preview nodes, study inputs, artifact transport and codecs were also inspected for semantic
  policy gates. Their remaining validations protect parsing, identity, file ownership, or typed artifact transport.

## Verification

Family regressions include real backend examples and generated-code parity (including PNG equality for changed
plots). The final frozen-tree run completed on 2026-09-05:

- `..\ComfyUI\.venv\Scripts\python.exe -m pytest -q --tb=short --disable-warnings -ra`:
  **1833 passed, 4 skipped**, 215 warnings, 450.54 seconds.
- The four skipped optional native smoke tests require the explicit Pertpy composition environment flag, exact
  LIANA 1.9.0, the supported Schist/graph-tool environment, or the audited scVelo wheel. Their existing contract and
  generated-code tests ran; the full-suite result does not claim these four optional integrations were executed.
- Frontend Node tests: **25 passed**. Example workflow generator: **6 verified**.
- Full-repository Ruff and `git diff --check`: passed.
- Independent Standards and Spec reviews: all findings closed. Review corrections covered actual Velocity
  normalization targets, the remaining MAD UI ceiling, DPT-to-Embedding infinity masking, and Augur numeric/dtype
  preservation in both backend input and source fingerprints.
- The original multi-step annotation reproducer now succeeds: T subtype merge, disjoint B subtype merge, then
  correction of the T subtype with explicit overwrite; runtime and generated outputs match.

No dependencies were added and no deployment was performed.
