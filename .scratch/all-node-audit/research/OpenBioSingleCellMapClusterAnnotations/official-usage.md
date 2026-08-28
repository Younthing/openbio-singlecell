# Map Cluster Annotations: official usage research

Researched: 2026-08-28

## Audited implementation and schema

The current node accepts an `AnnData`, a source observation column (`groupby`, default `leiden`), a JSON object,
an output observation column (`cell_type`), and `keep_unmapped`. It parses the JSON with the default
`json.loads`, coerces every key and value with `str`, maps `obs[groupby].astype("string")`, optionally fills
unmapped results with the original cluster label, constructs `pandas.Categorical` without explicit categories,
and writes the column into a full `adata.copy()`.

The primary output remains in slot 0, but there is no `summary` or equivalent `code` output. There are no direct
tests for this node. `finish_adata` records the normalized mapping and parameters in the generic OpenBio analysis
history, but the node does not record annotation status, source/output category order, effective versus unused
mapping entries, missing/unmapped cells, many-to-one merges, or overwrite decisions.

## Official pandas mapping semantics

`pandas.Series.map` substitutes values through a mapping while preserving the caller's index. With an ordinary
dictionary, values absent from the mapping become missing. `na_action="ignore"` propagates missing input values
without passing them through a callable. For categorical input, pandas can apply a function once per category.

- Official API: https://pandas.pydata.org/docs/reference/api/pandas.Series.map.html
- Categorical mapping API: https://pandas.pydata.org/docs/reference/api/pandas.Categorical.map.html

`Categorical.map` is not itself a complete implementation for this node. Its official contract says a one-to-one
mapping can remain categorical, but a many-to-one mapping returns an `Index`; any unmapped category also leads to
an `Index`. Many-to-one cluster-to-identity maps are legitimate and common, so the node must resolve values first
and explicitly construct the final categorical rather than rely on categorical rename/map return types.

The current `mapped.fillna(group_values)` is particularly important to disclose: it treats every missing result as
"preserve the source label." That includes an observed source cluster absent from the mapping. Source values that
were already missing happen to remain missing only because the fill value is also missing. These are distinct cases
and must be counted separately.

## Official categorical semantics

The pandas categorical guide and API establish the following:

- every value is either one of the declared categories or missing;
- missing is represented in values/codes and is not a category;
- categories must be unique;
- category order is the explicit `categories` order, not lexical value order;
- if categories are not supplied, pandas infers them (sorted if possible, otherwise by appearance);
- `Series.unique()` contains only observed values in first-appearance order, while `cat.categories` can include
  unused categories;
- operations such as `value_counts()` can include unused categorical levels with count zero;
- a categorical is unordered unless `ordered=True` is explicitly supplied.

- Official categorical guide: https://pandas.pydata.org/docs/user_guide/categorical.html
- Official constructor: https://pandas.pydata.org/docs/reference/api/pandas.Categorical.html
- Official category accessor: https://pandas.pydata.org/docs/reference/api/pandas.Series.cat.html

The current call `pd.Categorical(mapped.to_numpy())` delegates the public category order to inference. It drops the
source categorical's declared order and all unused source categories. For string outputs it will commonly sort the
realized labels, but that is an implementation-derived order rather than a scientific node contract. The refactor
must provide explicit, deterministic categories and separately report declared, observed, and unused source levels.

Scanpy's official Leiden interface stores cluster labels in `adata.obs[key_added]` as dtype `category`, so preserving
the declared source category order is the normal path for the default `groupby="leiden"` input.

- Official Leiden API: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html

## Official AnnData copy and storage semantics

`AnnData.obs` is the observation-aligned pandas `DataFrame`. `AnnData.copy()` is documented as a full copy,
optionally on disk. Creating an output copy and assigning only one observation column is the correct nonmutating
node shape; it also means the expression matrix, layers, raw snapshot, embeddings, graphs, variable annotations,
and observation order must remain scientifically unchanged.

- AnnData object and `obs`: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html and
  https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.obs.html
- Full copy: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.copy.html
- AnnData categorical on-disk format: https://anndata.readthedocs.io/en/stable/fileformat-prose.html#categorical-arrays

The AnnData categorical encoding stores the category array, the `ordered` flag, and integer codes; code `-1`
represents missing. Consequently category order, ordered/unordered status, and missingness are durable public data,
not cosmetic implementation details.

For predictable semantics the node should reject backed inputs and ask callers to use `to_memory()` first. This
avoids making the meaning and cost of `copy()` depend on a backing file and keeps generated code equivalent to the
in-memory plugin path.

## JSON mapping validation

RFC 8259 defines a JSON object as an unordered collection and says member names should be unique. It warns that
duplicate names are handled inconsistently across implementations. Category order must therefore never depend on
JSON object member order.

- JSON standard: https://www.rfc-editor.org/info/rfc8259/ (DOI: 10.17487/RFC8259)
- Python decoder: https://docs.python.org/3.12/library/json.html

Python's default decoder retains only one duplicate member and also accepts `NaN`/`Infinity` constants even though
they are outside the JSON standard. Its documented `object_pairs_hook` exposes ordered pairs so duplicate keys can
be rejected, and `parse_constant` can reject non-standard numeric constants.

The current coercion `{str(key): str(value)}` silently turns `null` into `"None"`, booleans into `"True"`/`"False"`,
numbers into text, and arrays/objects into Python repr-like strings. Those are not validated biological labels.
The reviewed mapping contract should instead require a non-empty JSON object with unique names and nonblank string
values, reject non-standard constants and nested/non-string values, and bound the input size/entry count. Duplicate
target values must remain allowed because they explicitly merge multiple source clusters into one annotation.

## Source labels, unknown keys, missing values, and unused categories

Source cluster labels may be strings, numeric categoricals, or other scalar categorical values. Because JSON object
keys are strings, the adapter must render each non-missing source label to one canonical trimmed string and reject
distinct source values that collapse to the same rendered key (for example integer `1` and string `"1"`). Missing
source labels are not the string `"nan"`, `"None"`, or `"<NA>"`; they remain true missing categorical codes and are
counted separately.

The implementation should distinguish four cases:

1. **Observed mapped level**: a source cluster with cells and an explicit mapping; apply it.
2. **Observed unmapped level**: handle through an explicit constrained policy, defaulting to error.
3. **Declared but unused categorical level**: no cell currently has the level; report it but do not create an output
   category solely because it exists.
4. **Mapping key unknown to both observed and declared source levels**: reject it as a likely typo or stale map.

A mapping entry for a declared-but-unused level may be accepted but must be listed as unused and must not make the
operation look effective. At least one observed level must be mapped. Under a preserve-source policy, a mapped
target label that equals an unmapped preserved cluster label must be rejected; otherwise two scientifically distinct
cases silently collapse into one output category.

## Deterministic output categories

The category-order policy must not depend on cell row order or JSON member order:

1. for categorical source data, use the declared source category order after canonicalization;
2. for a noncategorical source, use sorted canonical observed labels;
3. exclude unused source levels from the output traversal;
4. resolve each observed level through the mapping/unmapped policy;
5. deduplicate resolved annotation labels in first occurrence along that source-level traversal;
6. construct `pd.Categorical(values, categories=resolved_categories, ordered=False)` explicitly.

This supports intentional many-to-one cluster merges while producing a stable category array. Cell-identity labels
are nominal here, so `ordered=False` should be fixed. A biologically ordered hierarchy or trajectory is a different
analysis and should not be inferred from cluster IDs or JSON order.

## Scientific interpretation

Cluster IDs are algorithmic partitions, not cell-type truth. Luecken and Theis describe annotation as interpretation
using marker genes and external/reference knowledge, warn that clustering resolution changes the identity granularity,
and recommend retaining manual marker review even when automated annotation is used.

- Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems
  Biology*. 2019;15:e8746. https://doi.org/10.15252/msb.20188746

This agrees with the repository domain: Cluster marker evidence supports review; Provisional annotation must remain
distinct from Curated annotation; formal Condition contrast requires Sample-level inference within a reviewed
population. A deterministic mapping operation cannot verify marker evidence or expert review. Annotation status must
therefore be an explicit caller declaration, default to `provisional`, and be carried in AnnData provenance and the
summary. Selecting `curated` records an assertion, not software validation.

## Stable output and disclosure contract

The primary output should remain a copied `AnnData`, followed by strict JSON `summary` and equivalent Python `code`.
The summary should disclose:

- source/output columns, annotation status, mapping and unmapped/collision policies;
- source dtype, ordered flag, declared category order, observed levels and counts, unused declared levels, and missing
  source cells;
- supplied/effective/unused mapping entries, unknown/rejected entries, identity mappings, and intentional many-to-one
  target merges;
- mapped/unmapped/missing cell counts and proportions;
- output category order, per-category counts, output dtype/ordered flag, and overwrite status;
- that every cell in a source cluster receives the same label and that no expression or marker evidence was examined;
- limitations concerning cluster-resolution dependence, within-cluster heterogeneity, manual-review uncertainty,
  provisional-versus-curated status, and lack of ontology validation;
- the annotation-practice, pandas, AnnData, and JSON references plus dynamic Python, openbio-singlecell, pandas, and
  AnnData versions.

The generated function must parse and validate the resolved JSON under the same strict rules, preserve true missing
values, implement the same deterministic categorical order and collision policy, reject backed input, return a new
AnnData, and leave the caller unchanged. It should not require ComfyUI or plugin-only result wrappers.

## Software references

- McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the 9th Python in Science
  Conference*. 2010:56-61. https://doi.org/10.25080/Majora-92bf1922-00a
- Virshup I et al. anndata: Access and store annotated data matrices. *Journal of Open Source Software*.
  2024;9:4371. https://doi.org/10.21105/joss.04371
