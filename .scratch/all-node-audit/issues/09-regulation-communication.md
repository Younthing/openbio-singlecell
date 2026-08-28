# Regulation and candidate cell communication

Type: task
Status: resolved

Blocked by: 08

## Initial audit focus

decoupler 2.x migration, pySCENIC resource bundles and result contracts, stochastic binarization provenance, LIANA method-specific result semantics, and Sample-aware communication analysis.

## Comments

- 2026-08-28: Historical initial checkpoint: inventory and code audit were complete; per-node documents had not yet landed.
- 2026-08-28: Pre-change research complete for all 10 Batch09 node IDs; every node now has exactly
  `official-usage.md` and `design-review.md`. Audited targets are decoupler 2.2.0, external pySCENIC 0.12.1, and
  LIANA 1.9.0. Decisions: enhance CollecTRI ULM and TF ranking around a typed activity artifact; retire the
  incompatible in-process pySCENIC runner and make a safe manifest/hash-bound importer the sole supported seam; fix
  RSS orientation with an exact NumPy/SciPy implementation; make binarization nonzero-seeded/RNG-isolated; report
  imported motif-pruned regulons rather than recomputing pre-cisTarget modules; make LIANA Sample-resolved with an
  explicit organism/resource and direct table output; retire the hidden-`uns` results extractor; and make the dot
  plot consume typed results with method-specific score semantics.
- 2026-08-28: P0 research findings: both regulatory decoupler nodes call removed 1.x interfaces; pySCENIC 0.12.1
  fails in the configured NumPy 2.4 runtime (`np.object`/`np.float`), while its legacy motif CSV loader uses `eval`
  on external content; the current RSS table swaps group and regulon meanings; pooled LIANA ignores biological
  Sample replication; and the dot plot passes unsupported LIANA 1.9 keywords while treating aggregate ranks as a
  generic significance p-value. No runtime, test, dependency, workflow, or ledger change was made in this research
  pass.
- 2026-08-28: Completed the CollecTRI ULM and Rank TF Activities pair. CollecTRI now uses the exact decoupler 2.2.0
  `op.collectri`/`mt.ulm` public contracts, accepts a strict local signed-network snapshot or an explicit human
  official-loader policy, fingerprints bytes/canonical edges/the complete expression feature axis, guards dense
  adjacency/output memory, independently verifies ULM t statistics and observation-wise BH, and emits an immutable
  `OPENBIO_TF_ACTIVITY` artifact plus AnnData/summary/standalone code. Official generated code requires a
  caller-materialized hash-matching network and never repeats the download. Rank consumes only that typed artifact,
  uses exact 2.2.0 `pp.get_obsm`/`tl.rankby_group`, validates rest/single/list comparator semantics, independently
  recomputes all three supported tests and complete per-group BH families, and emits the full bounded canonical table,
  strict summary, and standalone code without rerunning inference. Focused verification: 40 tests passed under
  `-W error`, including a real decoupler 2.2.0 sparse toy smoke and adversarial backend mutation checks; Ruff,
  `py_compile`, and diff whitespace checks passed for the owned surface.
- 2026-08-28: Historical implementation checkpoint, superseded by the final open-expert-boundary closure: LIANA core
  froze provisionally. Communication became the sole executor: it exact-locks public
  `liana==1.9.0`/pandas `<3` signatures, runs both audited methods by biological Sample on a private AnnData copy,
  at this checkpoint bound full-gene normalized expression and caller roles, and consumed either an audited bundled-human selector or a
  strict local CSV, and records license/citation plus raw/canonical/effective SHA-256 accounting. It returns a
  defensive `OPENBIO_LIANA_RESULT` artifact, strict summary, and standalone equivalent code. Results is a
  fail-before-read migration shim; DotPlot consumes only the typed artifact, applies method-specific deterministic
  selection, and renders without rerunning inference. Focused fake/adversarial/generated tests pass 22/22 under
  `-W error`; the optional installed-package test skips in the normal environment, while an isolated real
  LIANA 1.9.0 + pandas 2.3.3 smoke passed both methods and DotPlot. The three design reviews contain exact
  legacy-to-current mapping/rejection matrices. At that checkpoint the ledger was not yet closed because the shared
  migration owner still had to implement atomic Communication→Results→DotPlot rewrites/rejections, current-schema
  no-ops, and idempotence; the Answer below records their completion. Release integration also still needed an audited LIANA optional extra, notices, and artifact manifest entry;
  no shared dependency or migration file was changed in this pass.
- 2026-08-28: pySCENIC five-node core frozen provisionally. The in-process runner is a fail-before-read retired shim;
  the importer is the sole supported boundary and accepts only a dedicated `openbio-singlecell/pyscenic-external-run/v1`
  manifest with strict finite/duplicate-key JSON, exact pySCENIC 0.12.1 declaration, bound official CLI arrays,
  streamed SHA-256/TOCTOU checks, finite numeric expression and adjacency
  CSVs, TF/motif/Feather-v2 resource validation, safe `frozenset`/target literal parsing, and exact cell/gene/regulon
  identities. In the final open-expert gate review, non-Raw finite expression, an omitted container digest, official
  CLI defaults/unseeded execution, zero adjacency importance, unrelated directory files, and one-group RSS are
  accepted with explicit warnings; expert thresholds outside `[0,1]` and zero-sum binarization are likewise retained
  with disclosed all-on/all-off consequences. Safety, shape/axis, nonfinite, hash/resource-family, typed provenance,
  and memory/row gates remain fail-closed. It emits immutable `OPENBIO_SCENIC_RESULT`, a defensive AnnData activity copy, strict summary, and
  standalone equivalent importer. RSS, deterministic one-process HDT-compatible binarization, and final
  motif-pruned regulon membership consume only this artifact and emit bounded primary results plus strict
  summary/code; binarization additionally emits `OPENBIO_SCENIC_BINARY`. The dead `OPENBIO_SCENIC_NETWORK` public
  type/module was removed because it has no active producer or consumer. Focused verification passed 21/21 under
  `-W error`, including real pySCENIC 0.12.1 CLI-parser, RSS, and HDT threshold smoke tests, fake/adversarial mutation and
  malformed-bundle cases, input ownership, strict JSON, and exact generated/runtime equivalence. Ruff, `py_compile`,
  and owned-path diff checks pass. `research/batch09-regulatory-migration-matrix.md` records exact fail-closed
  legacy/current matrices for these five nodes and CollecTRI/Rank TF, plus dependency and release-artifact guidance.
  At that checkpoint graph-wide atomic migration/block rules had not yet shipped; the Answer below supersedes this
  historical note.

## Answer

All ten Batch09 stable node IDs are researched, implemented, and closed. Active regulatory, external-import, and
LIANA nodes expose typed primary artifacts plus strict `summary` and equivalent `code`; Run pySCENIC and LIANA
Results are publicly deprecated fail-before-read compatibility shims. Exact current schemas are no-ops. Every legacy,
hybrid, partial, hidden, or semantically untranslatable signature is discovered through the whole serialized graph
and rejected before unrelated mutation; safe migrations share collision-free, finite-safe-integer allocation with
complete bidirectional link/reference closure and exact capacity reservation.

The central adversarial migration suite passes 140/140, including nested graph containers, sentinels, numeric
subgraph hosts, unsafe IDs/counters/state shapes, mixed counter metadata, dangling references, mutation-free failure,
and idempotence. This closure adds no analysis-history or Raw/current binding gate: nodes expose only their audited
source branches, and wherever Raw is exposed elsewhere an explicit Raw selection is sufficient; scientific
confidence/completeness concerns are reported.
