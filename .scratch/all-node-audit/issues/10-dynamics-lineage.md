# Trajectory, velocity, CNV, and lineage

Type: task
Status: resolved

Blocked by: 09

## Initial audit focus

Required graph/diffusion prerequisites, reproducible roots, scVelo workflow consolidation, infercnvpy expression/reference validation, protocol-specific lineage QC, and lineage multiple-testing correction.

## Comments

- 2026-08-28: Initial inventory and code audit complete; per-node pre-change documents pending.
- 2026-08-28: Pre-change research complete for all 17 Batch10 nodes (34 documents, exactly `official-usage.md` plus `design-review.md` per node). Decisions: retain/deepen 16 nodes; retire the executing `OpenBioSingleCellCNVStructure` composite into a nonexecuting migration shim and split its five analyses; keep `OpenBioSingleCellVelocityGeneRanking` but rename its display/claim to recovered-dynamics fit ranking. P0 findings include the truncated/default-wrong CellCycle program and counts-backed `use_raw`; DPT missing the DiffMap prerequisite and using row-order-dependent lossy group roots; scVelo 0.3.4 rejecting the current deprecated HVG arguments, possible dynamical-to-stochastic fallback, and absent package lock; no latent-time or CellRank kernel/estimator artifact boundary; infercnvpy missing a lock and explicit transformed full-gene/reference/coordinate contract; all-missing Cassiopeia cells passing QC, private Cassiopeia API use and unresolved package identity; tree reconstruction not consuming the QC character artifact or declaring root; expansion tests lacking family-wide BH; and EffectivePlasticity using node/leaf rather than edge denominators and omitting the published 2.5% state filter. No runtime, tests, dependencies, workflows, ledger or issue Status were changed.
- 2026-08-28: Historical implementation checkpoint, superseded by the later open-expert-boundary entry:
  `OpenBioSingleCellCellCycleScore` was deepened after its two documents landed. The bundled
  Regev/Seurat resource is now the immutable 43 S + 54 G2/M human-symbol program with a checked SHA-256; custom
  organism-specific exact programs remain available. The node requires a provenance-backed full-gene
  log-normalized source, rejects the count-backed Raw snapshot and HVG views, isolates Scanpy/global RNG state,
  namespaces three verified outputs, and emits a strict report plus equivalent code. Thirteen focused tests pass
  with warnings as errors; Ruff and bytecode compilation pass.
- 2026-08-28: Release-status audit restored Cell Cycle Score to provisional despite its frozen core: its design
  review requires a legacy seven-gene/custom-program and `use_raw` migration with provenance-aware fail-closed
  behavior, but the shared workflow migration layer does not yet contain that rule. The central migration owner
  will close it serially with the remaining Batch08-10 breaking schemas.
- 2026-08-28: Independent adversarial review of the first Diffusion Map/PAGA/DPT implementation found release
  blockers: no atomic legacy workflow migration; graph fingerprints that did not bind Scanpy's actual matrix;
  structural but not scientific backend validation; PAGA mutation of unrelated string annotations; an unreachable
  numeric group-root interface and row-order-sensitive medoid; stale downstream/sidecar overwrite state; and
  oversized/shared generated-code reports. The three design reviews now record the exact closure criteria. The
  ledger remains provisional until real-Scanpy, malicious-backend, migration, strict-report, and generated-code
  parity regressions pass.
- 2026-08-28: The three trajectory release blockers from that adversarial review are closed. Scanpy now consumes
  the exact canonical matrices used for graph identity; Diffusion Map eigenpairs, PAGA model-v1.2 weights/MST,
  and DPT diffusion distances are independently recomputed before provenance. PAGA runs on a minimal scratch
  AnnData and cleans obsolete size state; DPT exposes strict typed roots and uses stable cell-ID medoids; every
  standalone source is a node-specific dependency closure and reports are bounded. Batch10 migration performs
  all-graph preflight, preserves root/subgraph object/array links and widget exposures, upgrades old DiffMap/PAGA,
  and rewrites legacy DPT to `group_medoid` with a reused or inserted same-key Diffusion Map. Unsafe boundary,
  missing, mismatched, partial, or broken cases fail without mutation. Real Scanpy, malicious fake backend,
  strict-report/code-parity, process-state, overwrite, and migration/idempotence regressions pass.
- 2026-08-28: `OpenBioSingleCellDiffusionMap`, `OpenBioSingleCellPAGA`, and `OpenBioSingleCellDPT` were deepened
  as three separate graph-state transitions after their six documents landed. All three resolve and fingerprint
  the exact named distance/connectivity bundle and canonical observation axis. Diffusion Map now guards the full
  output bundle, verifies the spectrum, restores Scanpy's NumPy print/RNG globals, and writes graph-bound
  provenance. PAGA fixes model v1.2 and undirected mode, requires a complete categorical partition, handles
  unused levels only on the output copy, validates the official group-size sidecar and Scanpy's one-orientation
  spanning-forest encoding, and discloses every bounded group edge. DPT no longer computes or trusts an implicit
  basis: it requires matching OpenBio Diffusion Map provenance and one connected component, and uses either one
  exact cell ID or a deterministic nominated-population medoid over informative diffusion components. Each node
  now emits `adata + summary + code`; generated code is standalone from OpenBio helpers. Sixteen focused real
  Scanpy tests and the combined CellCycle/registration suite pass with warnings as errors; Ruff and bytecode
  compilation pass. The PAGA documents were corrected before implementation to record the canonical asymmetric
  storage of `connectivities_tree` and the official `{groupby}_sizes` sidecar.
- 2026-08-28: Historical implementation checkpoint, superseded by the later open-expert-boundary closure: CNV core
  redesign completed after all eight original/new-node research documents were present. Infer
  CNV is now one infercnvpy 0.6.1 expression-to-window operation over a provenance-backed full-gene logged source,
  exact integer Raw-count axis proof, explicit normal reference and biological Sample support, declared genome
  assembly, strict coordinates/window accounting, bounded memory, immutable tamper-evident `OPENBIO_CNV_STATE`,
  strict report and standalone code. `OpenBioSingleCellCNVStructure` is a nonexecuting zero-output migration shim;
  its atomic replacements are CNV PCA and CNV Group Score. PCA fixes audited ARPACK/uncentered settings and verifies
  every returned component by independent eigenpair and leading-spectrum calculations. Score accepts both the
  original typed state and downstream annotated AnnData, requires exact shared fingerprints/axes, then independently
  recomputes every infercnvpy group score without inventing tumor labels or cutoffs. The optional `.[cnv]` extra pins
  `infercnvpy==0.6.1`; the Pandas-3 object-index backend adapter and chrM behavior are documented. Eight focused
  tests plus registration/release checks pass with warnings as errors, including real 0.6.1 smoke, malicious
  backends, tampering, dense/sparse inputs, strict JSON and three generated-code parity paths. Legacy Infer CNV and
  CNV Structure workflow migration remains a release blocker, so the four ledger entries stay provisional.
- 2026-08-28: The seven RNA-velocity nodes were rebuilt as one explicit staged chain after all fourteen research
  documents were reread and corrected against the signed scVelo 0.3.4 wheel. `OPENBIO_VELOCITY_STATE` is now an
  immutable-by-interface, defensive and current-content-fingerprinted artifact spanning prepared abundances,
  exact named-graph moments, one explicit kinetic estimate or recovered dynamics, and directed velocity-graph
  evidence. Preparation validates raw integer-like layers, applies exact shared-gene filters and restores `X`;
  moments cannot reconstruct PCA/neighbors; model fallback fails closed; recovery fixes one deterministic gene
  family and bounded dense outputs; graph matrices are sparse/sign/range/support checked; ranking is a pure fit
  view; stream plotting uses a private copy/figure and makes no fate or lineage claim. Every node emits its primary
  result, strict JSON summary and self-contained equivalent source. Exact version/signature guards, fake/malicious
  stages, sparse/dense inputs, ownership/tamper checks, memory limits, global restoration and generated-code parity
  pass. A real CPU smoke from the audited wheel completed preparation, moments, deterministic velocity, recovered
  dynamics, dynamical velocity, directed graph, ranking and PNG rendering; the two recover and one dynamical
  velocity Pandas-3 compatibility seams were documented, signature-guarded, lock-scoped and restored. Focused
  velocity tests pass with warnings as errors (9 normal plus one gated real smoke), and registration passes 7/7.
  No dependency, manifest, workflow, example or web file was changed in this velocity batch.
- 2026-08-28: The CNV release blocker is closed with a fail-closed whole-graph migration. Legacy InferCNV now
  upgrades only from an exact user-authored source/Sample/assembly hint, retains nonempty reference labels and the
  historical window, emits typed CNV state/report/code, deletes the executable hint, and leaves runtime provenance
  validation mandatory. A directly paired legacy CNV Structure expands atomically to CNV PCA, reviewed Neighbors,
  Leiden, UMAP, and CNV Group Score; the original typed state and final annotated AnnData reach Score separately,
  while all prior Structure consumers retain their AnnData links. Historical defaults are explicit and the
  unrecorded PCA component count/version-dependent Leiden flavor are marked for review. Root, nested single-host
  subgraph, object/array link, proxy, fan-out, global-ID and idempotence cases pass; boundary, reroute, shared,
  unconnected, ambiguous, partial/mixed, disabled and malformed cases are mutation-free rejections. The full web
  migration suite passes 74/74, focused CNV runtime tests 8/8 with warnings as errors, registration 7/7, plus Node
  syntax, Ruff and bytecode compilation.
- 2026-08-28: Velocity core evidence is complete, but the seven ledger entries remain provisional until one atomic
  legacy workflow migration lands. Required cases are: remove only default/legacy-invalid HVG controls during
  preparation; map zero Moments dimensions to explicit graph reuse and block nonzero hidden-recompute requests;
  preserve explicitly serialized Estimate mode while inserting/rewiring the typed staged chain; map Recover to
  `velocity_genes` only when a verified upstream mask is established (otherwise block or require explicit `all`);
  rewire Graph/Ranking/Stream consumers to the correct typed stages; map only canonical `fit_likelihood` and block
  nondefault ranking columns; and rename Stream `groupby` to `color_key`. Migration must preflight every affected
  node/link before mutation, preserve widgets and IDs, be idempotent, and fail closed on unsafe partial chains.
- 2026-08-28: The four Cassiopeia cores were rebuilt as four atomic modules after their eight node-local documents
  were corrected against the stable public API, the pinned fork source, and the Jones/Yang method implementations.
  Character preparation/QC now owns bounded exact-byte parsing, declared cut sites and prior population, conflict
  and global-cell-identity checks, correct missing/observed denominators, all-missing/all-uncut handling, canonical
  state order, and an immutable-by-interface defensive `OPENBIO_CASSIOPEIA_CHARACTERS` artifact. VanillaGreedy
  consumes only that artifact, fixes the official missing-data callable, exposes only the two topology-affecting
  policies, isolates RNG, canonicalizes solver-generated internal IDs, and emits a tamper-evident defensive
  `OPENBIO_CASSIOPEIA_TREE`. Expansion independently verifies the exact coalescent probability for every eligible
  clade and applies stable BH to the complete within-tree family. EffectivePlasticity reproduces the published
  rare-state, pruning, unifurcation and exhausted-polytomy preprocessing, uses subtree edges as the denominator,
  computes every cell path score with a bounded dynamic program, and cross-checks official small parsimony without
  leaking global RNG. Each node emits its primary result, strict JSON summary and standalone equivalent source.
  Sixteen fake/adversarial/generated-code tests pass with warnings as errors; the last independent registration plus
  release freeze check passed 12/12. A real QC -> tree -> expansion -> plasticity smoke passes with `cassiopeia-mt==2.1.3`
  on WSL Ubuntu 24.04 / Python 3.12.3, including repeat-topology fingerprint and generated-code compilation.
- 2026-08-28: Cassiopeia dependency evidence defines a narrower release boundary than package metadata. The exact
  distribution is `cassiopeia-mt==2.1.3`, tag `v2.1.3-mt`, commit
  `09868dd04c073d81dc60611f6f91cdbab17783f1`, sdist SHA-256
  `5934043af6506a6439576816be49e170bc41dbe3766eaf49d14d05c3577a3e24`. Metadata declares Python `>=3.10,<4`,
  NumPy `>=1.22,<3`, and pandas `>=2.2.1`, but native Windows/Python 3.13 installation fails both a transitive
  `pysam` source build and the Cassiopeia sdist's own Windows build path. Verified release support is therefore
  Linux/Python 3.12; Windows must use WSL/Linux, while macOS and Python 3.13 remain unverified. The proposed
  fail-closed extra is `lineage = ["cassiopeia-mt==2.1.3; platform_system == 'Linux'"]`. Notice/manifest metadata
  must name the distribution separately from import package `cassiopeia`, include the exact source identity and MIT
  license, and declare `OPENBIO_CASSIOPEIA_CHARACTERS` plus `OPENBIO_CASSIOPEIA_TREE` as process-local wire types.
- 2026-08-28: Exact Cassiopeia legacy/current migration matrix is frozen for the central migration owner; no web
  migration was edited in this batch. (1) LineageQC changes output `[table]` to
  `[characters, table, summary, code]`, so every old output-0 table link must move to output 1. Directly reviewed
  widget renames are `tumor_column -> lineage_column`, `minimum_cells_for_summary -> minimum_cells`,
  `maximum_uncut_fraction -> maximum_uncut_fraction`, `allele_representation_threshold ->
  allele_representation_threshold`, and `percent_unique_threshold -> minimum_unique_fraction`. Migration is blocked
  unless a workflow-owned review payload supplies exact `cut_site_columns`, `prior_grouping_columns`,
  `missing_data_allele`, `maximum_missing_fraction`, `minimum_informative_character_fraction`, and both resource
  limits; `cut_sites_per_intbc`, private-filter `minimum_intbc_fraction`, `lineage_size_threshold`, and the
  wrong-denominator `percent_unsaturated_threshold` have no equivalent. (2) Reconstruct changes eight file widgets
  to `[characters, lineage_id, prior_transformation, collapse_mutationless_edges]` and appends summary/code after
  tree output 0. It requires one atomically inserted/reused current QC node with that complete reviewed payload;
  then `tumor -> lineage_id`, file/column/allele settings move upstream, and legacy solver behavior becomes
  `negative_log`, `False`, root `solver_generated`. Missing cut-site/missing/conflict semantics or the ambiguous
  `mutation_family_column` block migration. (3) Expansion keeps tree and `minimum_depth`, changes fractional clade
  cutoff to integer `minimum_clade_size`, changes raw-p cutoff to `fdr_threshold`, and appends summary/code after
  table output 0. The two changed statistical controls have no semantics-preserving mapping and require explicit
  review. (4) Plasticity keeps adata/tree and maps `annotation_key` directly; `output_key` is accepted only after
  collision review, `summary_key` is deleted, and new status/mode/fraction/overwrite/resource controls are required.
  Legacy values are mathematically invalid and must be marked stale and recomputed; absent provenance defaults only
  to `unknown/exploratory`, never Curated/report-grade, and no existing column may be overwritten automatically.
  The eventual whole-graph migration must preflight before mutation, be atomic and idempotent, preserve array/object
  links, subgraphs, exposures and global IDs, and leave every unprovable graph unchanged. Until that migration and
  its adversarial tests land, all four ledger entries remain provisional.
- 2026-08-28: The Cassiopeia execution gates were re-audited for the expert-tool policy. Hard failures are now limited
  to computation/type/identity/resource requirements, exact dependency/signature contracts, backend disagreement or
  mutation, and invalid tree/artifact structure. A computable lineage that misses minimum-cell, uniqueness or
  informative-character QC thresholds remains in `OPENBIO_CASSIOPEIA_CHARACTERS` with `status="warning"`; tree
  reconstruction carries those flags into its report. All-missing/all-uncut cells follow the explicitly selected
  fraction filters without an extra hidden rejection. Expert empirical-prior grouping is honored even when it omits
  lineage/intBC identity or includes a cut-site column, with assumptions prominently warned. Expansion permits a
  positive singleton cutoff and warns about its uninformative family. Plasticity excludes absent/null annotation
  leaves with null status, permits a one-state zero-plasticity result, and treats missing/malformed/mismatched Curated
  provenance in report-grade mode as a report warning rather than an execution gate. Conflicting allele identities,
  duplicate global cell identity, non-categorical phenotype type, fewer than two retained leaves, invalid probability
  domains, output collisions, and forged/tampered artifacts remain blocking because computation or data ownership is
  otherwise undefined. Updated focused tests pass 16/16 with warnings as errors, and the pinned Linux real-backend
  chain still passes.
- 2026-08-28: Cell Cycle Score was reopened under the expert-tool policy before migration freeze. It now accepts an
  explicit Raw, X, or layer source without Snapshot/normalization history, current-to-Raw feature equality, or a
  certified log/full-gene state. Count-like, signed, constant, low-coverage and cross-phase-overlap cases are retained
  with report warnings; finite real aligned values, at least one observed gene per phase, unambiguous within-set
  identifiers, collision safety and verified backend output remain hard. Generated code is self-contained. Legacy
  `use_raw=true|false` maps directly to Raw|X when serialized; dynamic connected booleans remain untranslatable.
- 2026-08-28: The Cell Cycle plus seven-stage velocity migration is frozen. Current Raw/X/layer and all typed velocity
  stages are exact no-ops; safe legacy chains preserve serialized model/source/view choices and unsafe hidden
  recomputation or incompatible wire semantics fail before mutation. Main-thread Cell Cycle/trajectory/registration
  verification passed 46 tests under `-W error`; the complete central migration suite passed 98/98. These eight
  ledger entries are resolved; Cassiopeia remains the only provisional family in Batch10.

## Answer

All 17 Batch10 node identities are complete. The four Cassiopeia nodes now have a graph-wide preflighted, atomic,
idempotent migration that preserves array/object links, subgraphs, direct/proxy exposures and global IDs. Reviewed
legacy QC preparation is reused or inserted only with exact static evidence; Expansion requires explicit integer/FDR
review, and invalid legacy plasticity values are marked stale for recomputation. Connected widget types are checked
across target input, link and origin output/subgraph boundary. Current interfaces are no-ops. The final independent
review reported P0/P1/P2=0; Cassiopeia targeted tests passed 7/7 and the then-complete central suite passed 110/110.
No analysis-history or Raw-binding gate/test was introduced.
