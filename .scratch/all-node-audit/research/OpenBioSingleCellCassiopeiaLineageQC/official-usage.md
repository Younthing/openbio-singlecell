# Cassiopeia character preparation and lineage QC: official usage

## Audited, installable implementation identity

The upstream repository currently declares `cassiopeia-lineage` 2.1.0 but does not publish that implementation on PyPI: the package named `cassiopeia-lineage` on PyPI stops at the incompatible 2019 release 1.0.4. The selected implementation is **`cassiopeia-mt==2.1.3`**, an MIT-licensed PyPI redistribution of the current YosefLab interface. The audited source is tag `v2.1.3-mt`, commit `09868dd04c073d81dc60611f6f91cdbab17783f1`; its sdist SHA-256 is `5934043af6506a6439576816be49e170bc41dbe3766eaf49d14d05c3577a3e24`. Runtime code must require exactly this distribution and inspect the public function signatures before reading or converting data.

The distribution metadata declares Python `>=3.10,<4`, NumPy `>=1.22,<3`, and pandas `>=2.2.1`, but metadata is not an installation guarantee. The real end-to-end smoke passed on Linux (WSL Ubuntu 24.04, Python 3.12.3, NumPy 2.5.2, pandas 3.0.5). Native Windows/Python 3.13 installation was tested and failed: the transitive `pysam` dependency had no compatible wheel and its source build failed, while `cassiopeia-mt --no-deps` also failed its own sdist build because of a hard-coded/missing Windows build path. Therefore the fail-closed release contract is Linux-only, with Windows users running a Linux/WSL environment; macOS and Linux Python 3.13 remain unverified. An optional dependency must use the exact pin and a Linux platform marker rather than advertising the package on every platform.

Sources:

- [Official Cassiopeia repository](https://github.com/YosefLab/Cassiopeia)
- [Pinned compatibility source tag](https://github.com/andrecossa5/Cassiopeia/tree/v2.1.3-mt)
- [Pinned PyPI distribution](https://pypi.org/project/cassiopeia-mt/2.1.3/)

## Public character-conversion contract

The maintained public calls are:

```python
indel_priors = cassiopeia.pp.compute_empirical_indel_priors(
    allele_table,
    grouping_variables=["intBC", "LineageGroup"],
    cut_sites=["r1", "r2", "r3"],
)

character_matrix, priors, state_to_indel = (
    cassiopeia.pp.convert_alleletable_to_character_matrix(
        one_lineage_allele_table,
        ignore_intbcs=[],
        allele_rep_thresh=0.98,
        missing_data_allele=None,
        missing_data_state=-1,
        mutation_priors=indel_priors,
        cut_sites=["r1", "r2", "r3"],
        collapse_duplicates=True,
    )
)
```

The character matrix is cells × retained target characters. Integer `0` means an observed uncut state, `-1` means missing, and positive integers identify observed mutation alleles separately within each character. The returned `priors` and `state_to_indel` are character-position keyed and must travel with the matrix. The converter's character order and positive-state numbering depend on allele-table encounter order, so reproducible execution requires a canonical row order before calling it.

`collapse_duplicates=True` is safe here only after preflight has rejected conflicting alleles for the same `(lineage, cellBC, intBC, cut_site)` identity; it then collapses identical technical duplicate records rather than silently resolving biological conflicts. The audited backend treats exact `"NONE"` and strings containing `"None"` as uncut, so non-exact reserved-token matches must be rejected before conversion. CSV nulls and an explicitly declared missing allele map to `-1` and may never be conflated with `0`.

Empirical-prior grouping variables define putatively independent mutation occurrences. `intBC` alone is appropriate only for one clonal population; adding the lineage/clone identity is appropriate only when those lineages are independent experiments. Priors are estimated once from the declared population, not separately from each target lineage, and the grouping columns, group count and prior fingerprint must be reported.

Sources:

- [Official allele-table conversion API](https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.pp.convert_alleletable_to_character_matrix.html)
- [Official empirical indel-prior API](https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.pp.compute_empirical_indel_priors.html)
- [Official reconstruction tutorial](https://cassiopeia-lineage.readthedocs.io/en/latest/notebooks/reconstruct.html)
- Jones et al., *Genome Biology* 2020, DOI [10.1186/s13059-020-02000-8](https://doi.org/10.1186/s13059-020-02000-8)

## QC denominators and blocking conditions

The converted matrix, not `n_unique_intBC × 3`, is the authoritative character set. For a cell with `M` retained converted characters:

- missing fraction = `count(state == -1) / M`;
- observed count = `count(state != -1)`;
- uncut fraction among observed = `count(state == 0) / observed_count`;
- mutation fraction among observed = `count(state > 0) / observed_count`.

An all-missing cell has no defined observed-state denominator, is explicitly counted, and has missing fraction 1; it follows the declared missing-fraction filter rather than an additional hidden gate. An all-uncut cell has uncut fraction 1 and no mutation information and likewise follows the declared uncut filter. At lineage level, unique-profile fraction uses retained cells as its denominator, while informative-character fraction uses all post-conversion characters as its denominator and counts columns with at least one positive mutation after cell QC. Every lineage in the file receives one terminal pass/warning/unavailable row; small or low-quality lineages must not disappear, and a computable warning lineage remains available to experts.

Cell identifiers must be unique globally across lineage IDs. Reusing one `cellBC` in two lineages makes the later AnnData join ambiguous and is blocking. Null/blank lineage, cell or intBC IDs, missing cut-site columns, non-string allele values other than true nulls, undeclared conflicts, no backend-usable retained character, non-finite/out-of-domain fraction parameters, oversized input, and dependency/signature mismatch also fail closed. Prior grouping choices and unmet QC recommendations are warnings with full summary disclosure, not execution gates.

Thresholds are declared workflow policies rather than universal biological standards. The report must disclose raw and retained denominators, allele/conflict handling, prior population, exact versions and fingerprints, and must not present a QC pass as proof that a reconstructed topology is true.
