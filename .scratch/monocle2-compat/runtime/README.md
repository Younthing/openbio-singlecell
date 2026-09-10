# Modern Monocle 2 dependency runtime

This experiment keeps the existing reference environment at
`/root/.cache/openbio-monocle2` unchanged. It clones that environment to
`/root/.cache/openbio-monocle2-modern` in the existing `Ubuntu-24.04` WSL
distribution, then replaces only the clone's igraph and dplyr R packages.
No system R, Python, shell configuration, or global PATH is changed.

## Reproduce the isolated setup

Run from the Windows repository root:

```powershell
wsl -d Ubuntu-24.04 --exec /root/miniconda3/bin/conda create -y --prefix /root/.cache/openbio-monocle2-modern --clone /root/.cache/openbio-monocle2
wsl -d Ubuntu-24.04 --exec env PATH=/root/.cache/openbio-monocle2-modern/bin:/usr/bin:/bin PKG_CONFIG_PATH=/root/.cache/openbio-monocle2-modern/lib/pkgconfig MAKEFLAGS=-j4 /root/.cache/openbio-monocle2-modern/bin/Rscript --vanilla /mnt/d/learn/openbio-singlecell/.scratch/monocle2-compat/runtime/install_modern.R
```

`install_modern.R` installs official CRAN source tarballs into the cloned
environment's R library. The existing machine has GCC/G++ 13 and the cloned
environment supplies its Fortran compiler. Compiler and pkg-config paths are
set only for that installation process.

The current CRAN package pages were checked on 2026-09-05:

- [igraph 2.3.3](https://cran.r-project.org/web/packages/igraph/index.html),
  published 2026-06-26. The indexed igraph changelog still showed 2.2.3, so
  the live package page and actual source tarball determine this test version.
- [dplyr 1.2.1](https://dplyr.tidyverse.org/news/index.html), published
  2026-04-03. dplyr 1.2.0 made the underscored verbs used by Monocle 2 defunct.

An initial conda installation of `r-igraph=2.2.3 r-dplyr=1.2.1` with
`--freeze-installed` failed because the igraph binary requires a newer
libxml2 ABI. No conda transaction was applied. Source installation keeps the
reference environment's remaining package versions constant. Accordingly,
`conda list` still describes the cloned R packages' original conda records;
use R's `packageVersion()` and the runtime probe for the installed versions.

## Native failure probes

The probe uses the existing deterministic, branched CDS fixture from the
Monocle 2 runtime smoke test. It calls the installed native Monocle functions
directly and does not source the production driver or compatibility adapter.
Each operation catches its own error so all failures are visible in one run.

```powershell
wsl -d Ubuntu-24.04 --exec /root/.cache/openbio-monocle2/bin/Rscript --vanilla /mnt/d/learn/openbio-singlecell/.scratch/monocle2-compat/runtime/reproduce_native_failures.R /mnt/d/learn/openbio-singlecell/.scratch/monocle2-subpopulation/runtime-setup/smoke-output/cds.rds
wsl -d Ubuntu-24.04 --exec /root/.cache/openbio-monocle2-modern/bin/Rscript --vanilla /mnt/d/learn/openbio-singlecell/.scratch/monocle2-compat/runtime/reproduce_native_failures.R /mnt/d/learn/openbio-singlecell/.scratch/monocle2-subpopulation/runtime-setup/smoke-output/cds.rds
```

The reference environment successfully runs `estimateDispersions`,
`orderCells`, and `BEAM`. Its loaded versions remain monocle 2.34.0,
DDRTree 0.1.6, igraph 2.0.3, dplyr 1.1.4, and Matrix 1.7-6.

The modern installation completed successfully. Both environments use
R 4.4.3. Comparing every installed R package confirms that only igraph and
dplyr versions differ:

| Package | Reference | Modern |
| --- | --- | --- |
| igraph | 2.0.3 | 2.3.3 |
| dplyr | 1.1.4 | 1.2.1 |

All three native failures are reproduced with the modern environment:

| Operation | Actual native error |
| --- | --- |
| `estimateDispersions(cds)` | `group_by_()` is now defunct; use `group_by()` |
| `orderCells(cds)` | `dfs()` argument `neimode` is now defunct; use `mode` |
| `BEAM(cds[1:12, ], branch_point=1, progenitor_method="duplicate", cores=1)` | `nei()` is now defunct; use `.nei()` |

After the modern source installations, the reference probe was rerun and
all three operations still succeeded. Its loaded package versions remained
unchanged. Local full output is in ignored `install.log`, `native-old.log`,
and `native-modern.log` files.
