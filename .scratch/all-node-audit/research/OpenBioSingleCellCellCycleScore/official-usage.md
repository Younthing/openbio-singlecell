# Cell Cycle Score: official usage and scientific practice

## Audited runtime and primary sources

The repository locks Scanpy 1.12.3. Its installed public interface is:

```python
scanpy.tl.score_genes_cell_cycle(
    adata, *, s_genes, g2m_genes, copy=False, **kwargs
)
```

The keyword arguments are forwarded to `scanpy.tl.score_genes`; in 1.12.3 that includes `gene_pool`, `n_bins=25`, `random_state=0`, `use_raw`, and `layer`. `ctrl_size` is deliberately fixed to the smaller of the two phase lists. The function writes `S_score`, `G2M_score`, and `phase` to `obs`.

Primary sources:

- [Scanpy `score_genes_cell_cycle`](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.score_genes_cell_cycle.html)
- [Scanpy cell-cycle tutorial](https://scanpy.readthedocs.io/en/stable/how-to/cell-cycle.html)
- [Scanpy `score_genes`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.tl.score_genes.html)
- [Seurat `CellCycleScoring`](https://satijalab.org/seurat/reference/cellcyclescoring)
- [Seurat cell-cycle vignette and marker provenance](https://satijalab.org/seurat/articles/cell_cycle_vignette.html)
- Satija et al., *Nature Biotechnology* 2015, DOI [10.1038/nbt.3192](https://doi.org/10.1038/nbt.3192)
- Tirosh et al., *Science* 2016, DOI [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)
- Wolf et al., *Genome Biology* 2018, DOI [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

## Official gene-set and expression-state practice

Scanpy's maintained tutorial uses the Regev-lab list of 97 human gene symbols: the first 43 are S-phase genes and the remaining 54 are G2/M genes. It explicitly keeps the full gene set rather than an HVG-only view, because losing roughly 70 genes made the score ineffective in the tutorial. It normalizes total counts, applies `log1p`, and then scores.

The current node's seven S and seven G2/M defaults are only the leading entries of that list. They are not the official complete set and materially change both the score and the matched control-gene sample. This is a release-blocking default error.

The maintained Scanpy workflow recommends full-gene log-normalized/log-transformed values. That recommendation must
be the default and must be prominent in the report, but a storage slot is not proof of biological state and an open
expert tool must not infer or enforce that state through mutable workflow history. Explicit Raw, X, and layer choices
are therefore valid execution choices when their matrix and axes are usable. Count-like, signed/residual-like,
constant, or otherwise unverified states are disclosed as interpretation warnings; selecting Raw is itself the
complete source-selection action and does not require a Raw-binding or history-attestation check.

The bundled Regev/Seurat markers are human gene symbols. Uppercasing mouse symbols is not an orthology mapping. A
scientific wrapper must declare the organism and gene-set identity and use exact feature identifiers without silent
case conversion. Human is the recommended and matched use. If an expert deliberately applies the literal human
program to a dataset declared as mouse/other, the calculation remains defined and should run with a prominent
resource/organism mismatch warning; no orthology or species-validity claim may be made. Version suffixes and
duplicate gene identifiers still need explicit handling.

## Algorithm and randomness

For each phase, `score_genes` subtracts the mean of expression-matched control genes from the mean phase-gene expression. Control genes are sampled within expression bins; the seed and candidate `gene_pool` therefore affect results. Both requested phase sets and the genes actually found must be disclosed. A very small surviving set materially changes the intended program and control size, so it must produce a prominent warning and cautious report language. It is not a computational error while at least one exact gene remains in each phase program and the backend can calculate a result.

Scores and the derived G1/S/G2M label are descriptive per-cell annotations. They are not cell-cycle-time estimates, synchronized-time measurements, or evidence for a `Condition` difference. Regressing these scores is a separate modeling decision and must not be performed in this node.

## Preconditions and reporting

- In-memory AnnData, unique observation and feature identifiers, nonempty cells/features.
- A real numeric, finite, axis-aligned source explicitly selected from Raw, X, or a layer. Full-gene log-normalized
  expression is recommended; count-like, residual-like, unverified, or feature-restricted sources remain expert
  choices and are disclosed rather than blocked.
- Explicit organism/gene identifier namespace and no string-normalization collisions. A nonhuman declaration with
  the literal human resource is an expert override reported as a mismatch, not a runtime error.
- Unique exact identifiers within each S and G2/M list and at least one observed gene in each phase. Cross-phase
  overlap and coverage below the documented 20/half bundled or 10-gene custom practice thresholds are warnings.
- No output-column collision unless overwrite is explicit.
- Report requested/found/missing genes, coverage by phase, control size, source state, seed, score quantiles, phase cell counts/proportions, warnings, references, and runtime versions.

The report must state that the result is a supervised program score, that cells are not biological replicates, and that no `Sample`-level `Condition` inference was performed.
