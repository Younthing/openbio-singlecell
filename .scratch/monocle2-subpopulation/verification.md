# Monocle 2 verification

Implemented the user-confirmed Monocle 2 / DDRTree scope against 57eb64e. The earlier R prototype remains separate
on codex/r-worker-prototype. No vendor source or production Python dependency was changed for this feature.

## Coverage

- Ten nodes: Runtime, Prepare, Ordering Genes, DDRTree, Order Cells, Trajectory Plot, Differential Test, BEAM,
  Gene Trends and Export to AnnData.
- Existing Python one-shot Worker, cancellation and Artifact publication; RDS native state plus typed JSON and CSV
  sidecars. Labels such as `001`, literal `NA`, empty strings, missing values, categories and native numeric fields
  retain their meaning. Floating infinities use valid quoted JSON tokens in float-typed table columns and restore
  through the existing table reader; metadata JSON remains strict.
- No inferred-counts or sample-size gate. Native advanced arguments are recorded, fixed DDRTree identity is retained,
  and upstream-ignored controls are not exposed as effective widgets.
- Initial native orientation is explicitly described as unoriented. Explicit empty root selection gets an actionable
  message. Native ordering is invalidated when its graph is recomputed; unrelated input annotations remain preserved.
- Original AnnData expression, feature axis, Raw, layers and metadata remain intact when attaching declared results.
- A 26-node subpopulation template includes initial State/annotation inspection, explicit rerooting, pseudotime,
  differential evidence, trends, H5AD/CSV/PNG saves and CDS persistence. The parent/population/environment and final
  root State are required user choices. The cover is generated artwork; analysis plots are actual Monocle output.

## Execution evidence

1. Native R driver: all 13 validation sections passed, including separate-process RDS stages, rerooting, BEAM,
   recomputation invalidation and typed feature metadata.
2. Final native Python/R checks: **29 passed in 61.00 seconds**, using isolated Python 3.12.2 with AnnData **0.13.2**,
   pandas **3.0.5**, and actual R **4.4.3** / monocle **2.34.0** / DDRTree **0.1.6** / igraph **2.0.3** / dplyr **1.1.4**.
   This covers runtime, codec, public operations and independent generated-code execution, not just compilation.
3. Full Windows regression: **1865 passed, 8 skipped**, with two failures caused only by the pre-existing hardcoded
   count of 29 expression-source nodes. The new Prepare makes this 30; those assertions and its counts-layer default
   classification were updated. Follow-up expression-source, release, registration and Monocle schema checks:
   **34 passed**. No production code changed after the full run.
4. The four new optional R checks skipped on Windows were executed successfully in the native Linux/WSL run above.
   The other four skips are the existing optional composition, LIANA, Schist and scVelo smoke tests.
5. Frontend checks: **25 passed**. Schema-generated workflows: **7 verified**. A real built wheel includes both
   `r/probe_runtime.R` and `r/monocle2.R` with bytes matching the source. Ruff and diff whitespace checks passed.
6. Independent Standards and Spec reviews have **zero open Monocle findings** after the typed-table correction.

The native integration fixture exercises 180 cells and 400 raw genes while its Python AnnData retains only 200 current
features. It produces a branched DDRTree, explicitly reroots, executes DE and BEAM, renders native plots, and verifies
the original H5AD bytes and matrix/layer/Raw contents after attachment. This is functional validation, not a large-data
performance benchmark. Native display jitter can legitimately change individual PNG pixels; scientific tables and
coordinates were compared independently of that rendering variation.

## Reproduce native checks

```powershell
wsl -d Ubuntu-24.04 -- env OPENBIO_TEST_RSCRIPT=/root/.cache/openbio-monocle2/bin/Rscript /root/.cache/openbio-monocle2-python/bin/python -m pytest --noconftest --confcutdir=tests tests/test_r_runtime.py tests/test_monocle2_codec.py tests/test_monocle2_operations.py tests/test_monocle2_reproduction.py -q
```

See runtime-setup/README.md for the isolated R environment and original native smoke data. Python and R run within
the same OS; native Windows Monocle 2 was not validated. The existing base WSL Python/R environments were retained.

Concurrent workspace changes to cache_policy.py, its tests, cache wording in README/ADR, and architecture diagrams
were not authored by this task and were excluded from the Monocle review. Their hunks were preserved. The full
regression ran against the shared worktree, including those concurrent edits.
