# Third-party notices

`openbio-singlecell` is based on and designed to run with the following upstream projects. Their copyright and license terms remain in effect. The complete GNU GPL version 3 text is included in `LICENSE`.

## ComfyUI

- Project: https://github.com/comfyanonymous/ComfyUI
- License: GNU General Public License version 3
- The paired source baseline is recorded in `release_manifest.json`.

## ComfyUI_frontend

- Project: https://github.com/Comfy-Org/ComfyUI_frontend
- License: GNU General Public License version 3 only
- The paired source baseline is recorded in `release_manifest.json`.
- A prebuilt paired frontend must include `dist/LICENSE` and the generated `dist/THIRD_PARTY_NOTICES.md`. Generate and verify them with `scripts/generate_frontend_notices.mjs` after every build. The generated notice inventories the frozen production dependency graph and deduplicates identical license texts.

## Scanpy

- Project: https://github.com/scverse/scanpy
- License: BSD 3-Clause License
- Required version: `>=1.12.3,<1.13`, including its `leiden` and `scrublet` extras

## igraph

- Project: https://github.com/igraph/python-igraph
- License: GNU General Public License version 2 or later
- Installed through Scanpy's `leiden` extra and used for validated graph layouts and the weighted Leiden backend

## leidenalg

- Project: https://github.com/vtraag/leidenalg
- License: GNU General Public License version 3 or later
- Installed through Scanpy's `leiden` extra; retained as part of the declared Scanpy clustering dependency set

## fa2-modified

- Project: https://github.com/AminAlam/fa2_modified
- License: GNU General Public License version 3
- Optional exact version: `0.4`
- Optional dependency used only when the Force-Directed Graph node explicitly requests ForceAtlas2; the node fails
  rather than silently substituting another layout when this package is unavailable

## harmonypy

- Project: https://github.com/slowkow/harmonypy
- License: GNU General Public License version 3 or later
- Optional exact dependency: `harmonypy==2.0.0`, source tag `v2.0.0` / commit
  `7335470f68f60a9719a0bcfc1230ff0ed0a8cd34`
- Official wheels are published for Linux and macOS, not native Windows; Windows users need WSL or a
  separately validated supported source build

## scvi-tools

- Project: https://github.com/scverse/scvi-tools
- License: BSD 3-Clause License
- Optional reviewed range: `>=1.5,<1.6`
- Imported only when an scVI analysis executes

## CellTypist

- Project: https://github.com/Teichlab/celltypist
- License: MIT License
- Optional reviewed range: `>=1.7,<2`
- Optional dependency imported only when CellTypist Annotation executes; model files are not bundled and must be
  installed locally before execution

## decoupler

- Project: https://github.com/scverse/decoupler
- License: BSD 3-Clause License
- Optional exact dependency: `decoupler==2.2.0`
- Imported by scoring, enrichment, CollecTRI, TF-activity and pseudobulk nodes; no knowledge resource is bundled

## Pertpy

- Project: https://github.com/scverse/pertpy
- License: MIT License
- Optional exact dependency: `pertpy==1.3.0`
- Used by Augur, drug perturbation, scCODA/tascCODA, Milo, edgeR and PyDESeq2 adapters; the composition, Milo and
  differential-expression extras select only the upstream optional stacks required by those families
- The Milo and edgeR paths also require a separately installed R/Bioconductor environment; pip and Pertpy do not
  install that external runtime

## Monocle 2 / DDRTree and R file exchange

- Monocle project: https://bioconductor.org/packages/monocle/ ; license: Artistic-2.0
- DDRTree project: https://cran.r-project.org/package=DDRTree ; license: Artistic License 2.0
- Integration tested with R 4.4.3, Monocle 2.34.0, DDRTree 0.1.6, igraph 2.3.3 and dplyr 1.2.1;
  unmodified Monocle with igraph 2.0.3 / dplyr 1.1.4 supplies the scientific reference
- File exchange uses R Matrix (GPL >=2 / its LICENCE) and jsonlite (MIT + LICENSE)
- These are optional, separately installed R packages. The repository ships its R driver and dependency API adapter.
  The adapter copies selected installed Monocle functions into a private R environment and translates obsolete calls;
  it does not modify the installed packages or their namespaces. Upstream algorithm source is not bundled.
  Monocle source: https://github.com/bioc/monocle/tree/RELEASE_3_20/R ; its source, copyright and Artistic-2.0 license
  remain with the upstream project. Standalone reproduction includes the same adapter and requires installed Monocle.

## PyDESeq2

- Project: https://github.com/owkin/PyDESeq2
- License: MIT License
- Optional reviewed range: `>=0.5,<0.6`
- Installed by the pseudobulk extra for the explicitly selected PyDESeq2 engine

## LIANA

- Project: https://github.com/saezlab/liana-py
- License: BSD 3-Clause License
- Optional exact dependency: `liana==1.9.0`, with its required `pandas<3` compatibility constraint
- Imported only by LIANA communication analysis and its typed-result plot adapter

## scVelo

- Project: https://github.com/theislab/scvelo
- License: BSD 3-Clause License
- Optional exact dependency: `scvelo==0.3.4`
- Imported only by the staged RNA-velocity processing and analysis nodes

## infercnvpy

- Project: https://github.com/icbi-lab/infercnvpy
- License: BSD 3-Clause License
- Optional exact version: `0.6.1`
- Imported only by the reviewed Infer CNV, CNV PCA, and CNV Group Score nodes; no reference genome or annotation
  resource is bundled

## cNMF

- Project: https://github.com/dylkot/cNMF
- License: MIT License
- Optional exact method-author distribution: `cnmf==1.7.1`
- Used by the managed-workspace cNMF Rank Survey and Consensus Programs adapter
- OmicVerse 2.3.1 is intentionally not installed: it requires `anndata<0.12.0`, conflicting with this release's
  `anndata>=0.13.2,<0.14` core contract

## Cassiopeia

- Project: https://github.com/andrecossa5/Cassiopeia
- License: MIT License
- Optional exact distribution on the verified Linux/WSL release boundary:
  `cassiopeia-mt==2.1.3; platform_system == 'Linux'`, source tag `v2.1.3-mt` / commit
  `09868dd04c073d81dc60611f6f91cdbab17783f1`
- Native Windows and macOS execution are not verified or claimed by this release
- Imported only by the typed lineage-character, tree-reconstruction, clonal-expansion and plasticity chain

## Schist

- Project: https://github.com/dawe/schist
- License: BSD 3-Clause License
- Audited exact version: `0.10.0`, installed from conda-forge or the official source together with graph-tool
- The PyPI distribution named `schist` is unrelated and must not be installed for these nodes

## graph-tool

- Project: https://git.skewed.de/count0/graph-tool
- License: GNU General Public License version 3
- External compiled backend required by Schist; official ordinary installation uses conda-forge on Linux/macOS or
  WSL on Windows

## pySCENIC

- Project: https://github.com/aertslab/pySCENIC
- License: GNU General Public License version 3 or later
- Audited external-workflow version: `0.12.1`, tag commit
  `ce41b61b6570490949bd12b514e9f6de46d19c1f`; it is not installed in the ComfyUI Python environment
- Run from the official `aertslab/pyscenic:0.12.1` image or a dedicated Python 3.10 environment pinned to
  `pyscenic==0.12.1`, `ctxcore==0.2.0`, `arboreto==0.1.6`, `loompy==3.0.8`,
  `numpy>=1.21,<1.24`, and `pandas>=1.3.5,<2`, then import the complete versioned and hash-bound result bundle
- `openbio_singlecell/scenic_binarization.py` contains GPL-3.0-or-later compatibility rewrites of that commit's
  `src/pyscenic/diptest.py` (`_gcm_`, `_lcm_`, `_touch_diffs_`, `dip`, and `diptst`) and a compatibility
  transcription of the HDT branch of `src/pyscenic/binarization.py::derive_threshold`. The only mechanical
  compatibility changes in the latter are replacing removed NumPy aliases (`np.msort` with `np.sort` and
  `np.float` with built-in `float`) while retaining the seed, dip-test, Gaussian-mixture, KDE, and minimization
  semantics. Standalone code generated by the SCENIC Binarization node reproduces these derived routines and
  carries the same source and license notice.
- Upstream `src/pyscenic/diptest.py` credits Johannes Bauer's `tatome/dip_test`, specifically commit
  `a0e3d448a4b266f54ec63a5b3d5be351fbd1db1c`, and Benjamin Doran's `BenjaminDoran/unidip`; those credits are
  retained here.

## ctxcore

- Project: https://github.com/aertslab/ctxcore
- License: GNU General Public License version 3
- Audited external pySCENIC dependency: `0.2.0`; Feather-v2 ranking databases and their licenses are not bundled

## scikit-image

- Project: https://github.com/scikit-image/scikit-image
- License: BSD 3-Clause License
- Installed through Scanpy's `scrublet` extra for automatic Scrublet threshold selection

## AnnData

- Project: https://github.com/scverse/anndata
- License: BSD 3-Clause License
- Required version: `>=0.13.2,<0.14`

## Pillow

- Project: https://github.com/python-pillow/Pillow
- License: Historical Permission Notice and Disclaimer (HPND)
- Required version: `>=11,<13`; used to fully verify PNG report artifacts before atomic preview or export

Python packages installed as transitive dependencies are not vendored by this repository and retain their own notices and license terms. Redistributors must include the applicable upstream license texts and notices with the packaged artifacts. This source-package notice is not a substitute for the generated notice shipped with a prebuilt frontend.
