# Monocle 2 dependency adapter verification

Verified on 2026-09-05 against HEAD `57eb64e` and the uncommitted Monocle 2 feature present at turn start.
The adapter is implemented in production and has not been committed. Unrelated cache-pressure and architecture
changes were preserved; see [spec and fixed-point snapshot](spec.md).

## Environments

Both R environments run inside the existing Ubuntu-24.04 WSL distribution, using R 4.4.3, Monocle 2.34.0 and
DDRTree 0.1.6. Only the installed igraph and dplyr R package versions differ:

| Runtime | Rscript | igraph | dplyr |
| --- | --- | --- | --- |
| Unmodified native reference | `/root/.cache/openbio-monocle2/bin/Rscript` | 2.0.3 | 1.1.4 |
| Modern dependencies | `/root/.cache/openbio-monocle2-modern/bin/Rscript` | 2.3.3 | 1.2.1 |

The modern environment was cloned and only those two packages were rebuilt from official CRAN source. The original
environment remains unchanged. The shorter README conda command, with no igraph/dplyr pins, was separately verified
with a successful dry-run resolving the same R/Monocle/DDRTree/igraph/dplyr versions; that fresh conda environment was
not installed. See [exact environment setup and native failures](runtime/README.md).

## Red-green evidence

- Before adaptation, actual modern native calls failed at `estimateDispersions` (`group_by_`), `orderCells`
  (`dfs(neimode=...)`), and BEAM (`nei`). The old native reference succeeded for all three.
- The existing cache behavior test failed when only the adapter contents changed; it passes after adding the
  adapter's hash to the existing node fingerprint.
- Independent emitted Python code initially failed because its R bundle lacked `monocle2_compat.R`; it now ships
  the adapter and completes the native analysis chain.
- The full modern operation test exposed `plot_cell_trajectory`'s additional defunct `select_` calls. The adapter now
  preserves both named numeric column positions and string selections; the original failing test passes.
- Native formula syntax such as a backticked `branch label` remains executable. A compact-RDS regression protects
  against serializing the private adapter through native state-container parent environments or the fitted closure.

## Final checks

1. **Scientific parity: 1 passed in 33.16 seconds** on the final adapter. The independent reference script calls
   installed native Monocle directly without loading this adapter. The deterministic 180-cell / 400-gene branching
   fixture compares size factors, grouped dispersion estimates and coefficients, DDRTree coordinates, graph edges
   and weights, cell projections, closest vertices, DFS parents, roots, State and Pseudotime before and after rerooting,
   and eight-gene DE / BEAM results. Numeric tolerance is `rtol=1e-8`, `atol=1e-10`; categorical values and structure
   match exactly. Every RDS is read by a fresh R process without loading the adapter. User-selected grouped formulas,
   covariate formulas, `remove_outliers=False` and `relative_expr=False` are covered.
2. **Modern public operations and standalone reproduction: 18 passed in 59.30 seconds.** This includes trajectory
   plots, gene trends, native RDS stages, typed artifacts and downstream AnnData attachment. The earlier combined
   runtime/codec/operation/reproduction run had 28 passed and one plotting failure; that failure was fixed and the
   complete operation/reproduction modules were rerun successfully. Runtime and codec checks passed in that run.
3. **Native R checks:** the complete modern driver passed all 13 sections. The dedicated adapter check passed all
   four sections in both reference and modern environments, including namespace isolation, quoted annotation names
   and compact saved state.
4. **Windows focused regressions: 60 passed, 2 optional native checks skipped in 8.43 seconds.** Covers node schemas,
   cache invalidation, codecs, public operations, runtime behavior, release contract, examples and a real built wheel.
   The wheel contains byte-identical probe, driver and adapter R resources. Native checks run under WSL as above.
5. Ruff and Git whitespace checks passed. Independent **Standards: 0 findings; Spec: 0 findings**, including a second
   review of trajectory plotting and RDS environment handling. Reviewers checked native container consumers and found
   only direct `$` / `[[` access, with no inherited-environment lookup dependency.

## Reproduce

Run from the repository root:

```powershell
wsl -d Ubuntu-24.04 -- env OPENBIO_TEST_RSCRIPT=/root/.cache/openbio-monocle2-modern/bin/Rscript OPENBIO_REFERENCE_RSCRIPT=/root/.cache/openbio-monocle2/bin/Rscript /root/.cache/openbio-monocle2-python/bin/python -m pytest --noconftest --confcutdir=tests tests/test_monocle2_compatibility.py -q
wsl -d Ubuntu-24.04 -- env OPENBIO_TEST_RSCRIPT=/root/.cache/openbio-monocle2-modern/bin/Rscript /root/.cache/openbio-monocle2-python/bin/python -m pytest --noconftest --confcutdir=tests tests/test_r_runtime.py tests/test_monocle2_codec.py tests/test_monocle2_operations.py tests/test_monocle2_reproduction.py -q
wsl -d Ubuntu-24.04 -- /root/.cache/openbio-monocle2-modern/bin/Rscript --vanilla tests/test_monocle2_r_driver.R
wsl -d Ubuntu-24.04 -- /root/.cache/openbio-monocle2-modern/bin/Rscript --vanilla tests/test_monocle2_compat.R
D:\learn\ComfyUI\.venv\Scripts\python.exe -m pytest tests/test_monocle2_nodes.py tests/test_monocle2_codec.py tests/test_monocle2_operations.py tests/test_r_runtime.py tests/test_r_package_resources.py tests/test_release.py tests/test_examples.py -q
```

The full repository suite was not repeated for this localized adapter change. Earlier feature-wide verification is
recorded [separately](../monocle2-subpopulation/verification.md). These results demonstrate equivalence on the tested
branching fixture and package versions, not every possible dataset or future upstream release. Python and R must
run in the same OS; native Windows Monocle 2 remains unverified. Native plot jitter can change individual pixels.
