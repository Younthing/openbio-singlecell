# Normalize Gene Names: official usage research

Researched: 2026-08-28

## Official storage operations versus gene identity

Pandas supports mechanical string replacement through `Series.str.replace`/`Index.str.replace`, and AnnData permits replacing `var_names`. AnnData also provides `var_names_make_unique()`, which keeps the first occurrence and appends numeric suffixes to later duplicates. These are data-structure operations; neither defines a scientifically valid, organism-independent normalization of gene identifiers.

- pandas string replacement: https://pandas.pydata.org/docs/reference/api/pandas.Series.str.replace.html
- AnnData unique variable names: https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.var_names_make_unique.html
- AnnData variable-axis semantics: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html

Official nomenclature sources instead distinguish curated symbols from stable identifiers. HGNC assigns approved human gene symbols and stable HGNC IDs and advises authors to quote identifiers because symbols can change. Ensembl similarly distinguishes stable IDs from names, which can change as annotation knowledge changes.

- HGNC nomenclature guidelines: https://www.genenames.org/about/guidelines/
- HGNC instructions to authors: https://www.genenames.org/useful/instructions-to-authors/
- Ensembl stable IDs: https://mart.ensembl.org/info/genome/stable_ids/index.html

GENCODE gene IDs can carry a numeric version and a pseudoautosomal suffix, for example an identifier shaped like `ENSG....7_PAR_Y`. Version reconciliation is therefore a narrow identifier-matching operation; an arbitrary replacement such as changing every period to a hyphen is not equivalent to removing the terminal version and can corrupt valid identity strings.

- GENCODE GTF format and identifier conventions: https://www.gencodegenes.org/pages/data_format.html

## Method and database references

There is no publication supporting generic regular-expression replacement as “gene-name normalization.” Relevant authoritative references are:

- Seal RL et al. Genenames.org: the HGNC resources in 2023. *Nucleic Acids Research*. 2023;51(D1):D1003-D1009. https://doi.org/10.1093/nar/gkac888
- Frankish A et al. GENCODE 2025. *Nucleic Acids Research*. 2025. https://doi.org/10.1093/nar/gkae1078
- Martin FJ et al. Ensembl 2025. *Nucleic Acids Research*. 2025. https://doi.org/10.1093/nar/gkae1071

## Scientific implications

- Gene symbols, Ensembl IDs, HGNC IDs, transcript IDs, and vendor feature names are different identifier systems; one regex cannot normalize them all.
- `make_unique()` creates technical index labels, not new biological identifiers. Artificial suffixes must never be reported as approved gene symbols.
- Rewriting only current `var_names` after a Raw snapshot exists creates incompatible current and `raw` feature namespaces.
- Ensembl version handling should occur only while matching a declared stable-ID source to a declared annotation resource, with exact matches preferred and version-normalized matches counted.
- Consumers that require gene symbols should select an explicit `.var` annotation or adapt a private working copy; they should not require a global feature-axis rewrite.

## Open expert-boundary re-review (2026-08-28)

Pandas string replacement and AnnData index assignment are defined for empty patterns, empty axes, duplicate results,
and objects carrying an independent `raw` snapshot. Those states can be scientifically confusing but are not storage
errors. During this deprecated node's compatibility period, preserve the expert's exact replacement choice, allow
duplicate output when `make_unique=False`, and leave `raw` unchanged with explicit warnings. Do not infer or enforce
that current and Raw feature names must match.
