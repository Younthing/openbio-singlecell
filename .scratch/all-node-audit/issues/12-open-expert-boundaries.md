# Open expert-tool validation boundaries

Type: task
Status: resolved

## Initial audit focus

Re-audit every active node after the user clarified that OpenBio is an open expert tool. A hard error is justified
only when execution cannot be defined safely or faithfully (for example: incompatible type/shape, axis mismatch,
non-finite values unsupported by the method, invalid algorithm choice, unsafe path, or a backend's documented hard
precondition). Methodological preferences, provenance confidence, recommended thresholds, data-quality concerns, and
reporting best practices must remain user-selectable and be represented as warnings and explicit `summary`
limitations instead of runtime gates.

## Acceptance criteria

- Inventory hard-error conditions for all active analysis/processing nodes and classify each as computationally
  necessary or advisory.
- Downgrade advisory conditions to warnings plus cautious, machine-readable summary fields without silently changing
  the user's inputs.
- Do not turn report-template claims into programmatic claims that the implementation cannot actually prove.
- Preserve structural, alignment, tamper, path-safety, and backend-compatibility checks that prevent silent wrong
  results or unsafe execution.
- Cover representative relaxed expert choices and retained essential failures in runtime and generated-code tests.

## Comments

- 2026-08-28: Claimed immediately after the user added the open-tool boundary to the active goal. Augur was corrected
  first: explicit `Raw` selection is sufficient; mutable analysis history is not a runtime gate, and full-gene
  completeness is reported as user-declared rather than programmatically verified.
- 2026-08-28: Augur's non-integer selected source, low Sample support, repeated-measure Sample mapping, multi-batch
  Sample mapping, and global/population Condition-batch confounding were downgraded to warnings plus machine-readable
  audits. X/layer are now genuine explicit expert choices rather than aliases permitted only when byte-equal to Raw.
  Required two-arm population structure, Pertpy cell/subsample/fold support, finite nonnegative numeric state,
  canonical aligned axes, resource guards, and exact audited backend interface remain hard. Focused + registration
  regression: 62 passed with warnings treated as errors; Ruff clean.
- 2026-08-28: Augur now also permits an expert to reuse the selected `sample_key` for another observation role.
  The calculation remains defined, but `sample_audit.identity_status=role_reused_unverified`, `reused_roles`, and a
  report warning prevent those labels from being described as proven biological replicates. Pertpy is release-pinned
  at 1.3.0 but runtime-compatible across capability-checked 1.3.x patch releases; incompatible public signatures
  still fail. The user explicitly declined analysis-history tamper and Raw-binding boundary tests: choosing Raw is
  the complete source-selection action. Augur focused regression: 37 passed under `-W error`; Ruff clean.
- 2026-08-28: PCA no longer blocks caller-selected expression merely because an OpenBio state record declares
  counts, normalized values, or cannot prove the state. It emits a reproducibility/method warning and performs PCA
  on the explicit source; numeric/finite matrix shape, dimensions, collision/resource safety, and backend
  postconditions remain hard. Embedding regression: 33 passed under `-W error`; Ruff clean.
- 2026-08-28: Harmony accepts one observed or singleton Technical-batch strata with explicit warnings, and accepts
  capability-checked `harmonypy` 2.0.x patch releases while retaining exact call-signature, version-consistency,
  orientation, and output-shape checks. scVI likewise reports one-level/singleton categorical covariates and finite
  constant continuous covariates as advisories rather than blocking expert execution; blank/non-finite covariates,
  axis mismatch, and backend incompatibility remain hard. Integration regression: 36 passed under `-W error`;
  generated source compiled and Ruff was clean.
- 2026-08-28: Augur and Population Centroid Correlation `code` outputs were deepened from package-import wrappers to
  self-contained equivalent source. The generated programs carry their validation, backend checks, canonical result
  construction, references, software-version reporting, and strict summary logic without importing OpenBio. This
  closes an output-contract defect without adding a scientific gate. Focused regressions: 37 + 18 passed under
  `-W error`; Ruff and source compilation clean.
- 2026-08-28: Cell Cycle Score no longer requires mutable Snapshot/normalization history, current-to-Raw feature-axis
  equality, a certified log state, or recommended gene-coverage thresholds. Explicit Raw/X/layer sources run when
  axes and finite real values are valid; count-like, signed, constant, feature-completeness-unverified, low-coverage,
  and cross-phase-overlap choices are preserved and reported. At least one observed gene per phase, unique identifiers
  within each declared set, output safety, and complete finite backend results remain hard. Its `code` is now
  self-contained and Raw was covered through the public node plus generated source. Focused + adjacent registration
  regression: 46 passed under `-W error`; Ruff and diff checks clean.
- 2026-08-28: Batch08's complete 13-node pass downgraded count-like/Raw recommendations for applicable scoring
  kernels, low GSEA permutations/seed zero, low-but-computable Welch replication, batch confounding, resource review,
  and scientific producer/approval claims to warnings and explicit machine-readable interpretation fields. Poisson
  GSVA integer support, Welch's two-observation mathematical minimum, exact evidence/universe pairing, resource/file
  safety, current-content fingerprints and verified backend arithmetic remain hard. Main-thread verification: 143
  Batch08 + registration tests under `-W error`, with Ruff, pycompile and diff checks clean.

## Answer

The whole active-node inventory now follows the open expert-tool boundary. Hard failures are limited to conditions
needed for defined computation, semantic/axis alignment, backend compatibility, memory/path safety, output ownership,
or an exact typed-artifact contract. Recommendations about expression state, count-likeness, feature completeness,
replication, balance/confounding, thresholds, seeds, annotation confidence, resource approval, and report-grade use
remain executable whenever the method is defined; they are preserved in warnings, audits, key results, and
limitations without silently rewriting the user's choices.

Explicit Raw selection is the complete Raw-source action. No analysis-history tamper test, late-history gate, or
Raw/current binding/equality test was added. Raw full-gene/post-QC status is user-declared and disclosed as unverified
where it cannot be proven. The final 1151-test Python suite, 140-test graph-migration suite, strict report probes,
generated-code compilation/parity coverage, and clean Ruff/release checks close this boundary.
