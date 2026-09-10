# Monocle 2 subpopulation trajectories

## Confirmed request

Add production R analysis nodes and a subpopulation best-practice workflow using the accepted Python one-shot
Worker -> Rscript -> result files -> existing Artifact publisher route. The user explicitly selected **Monocle 2 /
DDRTree**, not Monocle 3. Preserve the repository's expert-facing policy: native executable choices remain available,
recommendations are advisories, and no inferred counts-state or sample-size gates are added.

## Workflow and node seams

- Reuse Load H5AD, Subset Observations and Raw Snapshot to AnnData for selection and an explicit full-feature source.
- Monocle 2 Runtime selects Rscript and optional child-process library/home/PATH settings, probes actual required
  packages, and returns a compact JSON STRING descriptor. Runtime/package/script identity participates in caching.
- Monocle 2 Prepare constructs a native CellDataSet from an explicit X/layer/Raw source. Family, size-factor,
  dispersion and gene-detection processing are explicit choices. Preserve feature IDs and observation metadata.
- Monocle 2 Ordering Genes supports explicit gene IDs, a feature-metadata flag or the native dispersion table.
  It marks ordering genes without removing other genes from the CellDataSet; thresholds remain user parameters.
- Monocle 2 DDRTree performs native reduceDimension with explicit native parameters and optional advanced arguments.
- Monocle 2 Order Cells supports native initial orientation, explicit root_state after initial ordering, and reverse.
  State identifies trajectory segments, not cell types. Initial/root-free pseudotime must be described as unoriented.
- Monocle 2 Trajectory Plot exposes native State, Pseudotime, annotation and gene coloring as appropriate.
- Monocle 2 Differential Test exposes native full/reduced formulas for pseudotime-associated genes; preserve native
  status/NA results. Monocle 2 BEAM provides optional user-selected branch-point analysis.
- Monocle 2 Gene Trends plots user-selected genes or a visible top-N selection from a supplied result table.
- Monocle 2 Export to AnnData attaches only declared cell-level results and DDRTree coordinates to the original
  selected AnnData, preserving its expression/feature axis/Raw/layers and unrelated analysis metadata.
- Every scientific node has the conventional primary outputs, summary and fully reproducible code. R source must be
  shipped as a package resource and included in the emitted reproduction path; do not copy the throwaway prototype
  into production or require Python-in-R conversion libraries.

## File and runtime contract

Use a dedicated OPENBIO_MONOCLE2_CDS artifact, codec monocle2-cds-rds-v1, containing native state.rds, compact metadata
and ordinary result sidecars. Native RDS save/read is verified between fresh R processes. Existing RunLease owns
atomic publication, cancellation and cleanup. Large expression/embedding data stay in files, not JSON summaries.
No universal multi-language Worker framework, vendor patches, live R objects in ComfyUI, or automatic scientific
package installation is part of this feature. R errors and warnings remain visible.

Runtime controls are explicit and record actual versions/capabilities, rather than enforcing the recommended tested
versions as an arbitrary gate. R package or bundled-script changes must invalidate relevant node caches. Expert
advanced arguments can override visible native defaults when this is documented and effective arguments are reported;
they cannot replace the owned CDS object or output paths.

## Example

Package one schema-generated Monocle 2 subpopulation workflow. Start from an annotated parent with user-verified
count Raw, subset a selected population, create the CDS, select ordering genes within that population, recompute
DDRTree, inspect State and known annotations, then explicitly reroot before presenting final pseudotime. Include
trajectory diagnostics, pseudotime differential evidence, gene trends, and explicit H5AD/CSV/PNG/CDS persistence.
Do not infer the biological root or treat pseudotime association/BEAM as replicate-aware Condition inference.
Selection, root State and environment paths are user-editable; the template should clearly identify required edits.

## Verification

Vertical TDD at the existing public operation/node and emitted-code seams. Real Monocle 2 smoke data must exercise
DDRTree, root selection, RDS roundtrip, downstream Python attachment, plots and gene evidence. Retain structural,
axis, finite-serialization, file-ownership and genuine backend checks. Run focused tests, full repository checks,
workflow generation, release packaging checks and independent standards/spec review.

The validated starting environment is local WSL R 4.4.3, monocle 2.34.0, DDRTree 0.1.6, igraph 2.0.3 and dplyr 1.1.4.
New igraph/dplyr removals break current upstream Monocle 2. This native environment remains the reference;
the subsequently authorized [local dependency adapter](../monocle2-compat/spec.md) supports modern dependencies
without modifying the installed packages.

## Fixed point

57eb64e; production worktree was clean. Earlier R prototype evidence is on codex/r-worker-prototype at faf1fef.
