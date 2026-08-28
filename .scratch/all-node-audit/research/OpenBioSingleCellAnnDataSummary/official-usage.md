# AnnData Summary: official usage research

Researched: 2026-08-28

## Official AnnData structure

AnnData is an annotated observation-by-variable matrix. The official interface defines:

- `X`: the primary observations-by-variables matrix;
- `obs` / `var`: one-dimensional axis annotations and their indexes;
- `obsm` / `varm`: multi-dimensional annotations aligned to observations/variables;
- `obsp` / `varp`: pairwise matrices aligned to one axis;
- `layers`: alternate observation-by-variable matrices;
- `raw`: a stored `X`/`var` snapshot;
- `uns`: unstructured annotations;
- `is_view`, `isbacked`, and `filename`: storage/state properties.

- Official `AnnData` API: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html
- Official on-disk schema: https://anndata.readthedocs.io/en/stable/fileformat-prose.html
- Getting-started examples: https://anndata.readthedocs.io/en/stable/tutorials/notebooks/getting-started.html

## Software reference

Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data matrices. *Journal of Open Source Software*. 2024;9(101):4371. https://doi.org/10.21105/joss.04371

## Reporting contract

This node performs a structural inventory, not a biological analysis. It should report enough state to review whether downstream methods have the expected representation without interpreting expression values as counts, normalized abundance, cell quality, cell type, or differential signal.

The structured `summary` should include:

- shape; `X` presence, dtype, storage class, sparse/dense state, and matrix shape;
- observation and variable index names, dtypes, uniqueness, missing/blank counts, and bounded examples;
- bounded inventories of `obs`, `var`, `layers`, `obsm`, `varm`, `obsp`, `varp`, and `uns`, with shapes/dtypes/storage classes where defined;
- `raw` presence and dimensions; view/backed status;
- bounded OpenBio display/source/history provenance when available;
- methods text, results text, key structural results, warnings/limitations, AnnData references, and runtime software versions.

All values must be strict JSON data: no NumPy/Pandas scalars, arrays, `NaN`, or infinity. Large key/name collections are bounded and marked with total/truncation counts. The equivalent `code` output should implement the same read-only inventory and return the same scientific/structural payload; it need not recreate plugin-only result wrappers.

