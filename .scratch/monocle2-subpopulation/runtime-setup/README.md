# Isolated Monocle 2 runtime validation

The local Windows computer already had the `Ubuntu-24.04` WSL distribution.
Only `/root/.cache/openbio-monocle2` was created inside it; existing R, Python,
ComfyUI, global PATH and shell configuration were not changed.

## Exact compatible setup

```powershell
wsl -d Ubuntu-24.04 -- /root/miniconda3/bin/conda create -y --prefix /root/.cache/openbio-monocle2 --override-channels -c conda-forge -c bioconda bioconductor-monocle=2.34.0 r-igraph=2.0.3 r-dplyr=1.1.4 r-jsonlite
```

Resolved R 4.4.3, monocle 2.34.0, DDRTree 0.1.6, igraph 2.0.3, dplyr 1.1.4,
Matrix 1.7-6 and ggplot2 4.0.3. See `conda-linux-64-explicit.txt` for every
package URL and `smoke-output/sessionInfo.txt` for the actual loaded packages.

The igraph and dplyr pins are real upstream compatibility constraints, rather
than scientific input restrictions:

- Monocle 2 source calls `nei()` and `graph.dfs(neimode=...)`, which become
  defunct in newer igraph releases. The current Bioconductor source still
  contains those calls.
- The first actual smoke failed inside `estimateDispersions()` because
  dplyr 1.2 made its `group_by_()` call defunct. dplyr 1.1.4 passes.
- No package source, namespace or function was patched or monkeypatched.

Package sources:

- <https://api.anaconda.org/package/bioconda/bioconductor-monocle>
- <https://api.anaconda.org/package/conda-forge/r-igraph>
- <https://bioconductor.org/packages/3.20/bioc/html/monocle.html>
- <https://raw.githubusercontent.com/bioc/monocle/RELEASE_3_22/R/order_cells.R>

## Reproduce the actual pipeline

```powershell
wsl -d Ubuntu-24.04 -- /root/.cache/openbio-monocle2/bin/Rscript --vanilla /mnt/d/learn/openbio-singlecell/.scratch/monocle2-subpopulation/runtime-setup/smoke_monocle2.R /mnt/d/learn/openbio-singlecell/.scratch/monocle2-subpopulation/runtime-setup/smoke-output
```

The deterministic fixture contains 180 cells and 400 genes, simulated with
negative-binomial counts and one progenitor population splitting into two
lineages. The real pipeline runs `newCellDataSet`, `estimateSizeFactors`,
`estimateDispersions`, `detectGenes`, `setOrderingFilter`, `reduceDimension`
with DDRTree, `orderCells`, an explicit `root_state` reroot, and
`differentialGeneTest(~sm.ns(Pseudotime))` on 12 genes.

Verified: three states, a genuine Y-shaped tree, 180/180 finite pseudotimes,
12 differential gene rows, RDS roundtrip, and an actual PNG rendered from
`plot_cell_trajectory`. The reroot uses the state enriched in the fixture's
known progenitors; automatic orientation would not establish biological time.

Separately verified `BEAM(cds[1:12,], branch_point=1,
progenitor_method="duplicate", cores=1)` on the saved branched CDS. All 12
rows returned `status="OK"`; results are in `smoke-output/BEAM.csv` and
the command output is in `beam.log`.

For Python integration tests, existing `/root/miniconda3/bin/python` supplies
Python 3.12.2, pytest 7.4.4, anndata 0.12.10, NumPy 2.2.4 and SciPy 1.15.0.
Use Linux repository paths (`/mnt/d/learn/openbio-singlecell`) when invoking R
from that Python process.

An initial Monocle 3 install was interrupted when the user clarified Monocle 2.
No Monocle 3 production node or runtime dependency was added.
