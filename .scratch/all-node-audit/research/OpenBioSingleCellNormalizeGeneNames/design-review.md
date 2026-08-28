# Normalize Gene Names: module design review

> Final current-only disposition (2026-08-29): deleted completely and not registered. Evidence-backed GTF
> annotation remains the supported identity operation; no compatibility node is retained.

## Current module

The node applies an arbitrary literal or regular-expression replacement to every `var_name` and optionally calls `var_names_make_unique()`. Its default changes every period into a hyphen. This is a generic string editor presented as a biological identity operation.

## Correctness and interface problems

- The default transformation is neither Ensembl version removal nor approved gene-symbol normalization.
- Arbitrary regex can alter valid symbols, vendor IDs, antibody features, guide IDs, or multi-species prefixes.
- `make_unique=True` silently creates artificial names and hides collisions rather than resolving gene identity.
- Existing `.raw.var_names` are not changed, creating two incompatible feature namespaces.
- The large interface (`pattern`, `replacement`, regex mode, collision behavior) exposes string-processing mechanics while providing no organism, identifier system, annotation release, or mapping evidence.

## Decision: delete; absorb only narrow ID matching into the GTF module

There is no coherent scientific contract for a general “normalize gene names” node. It should not be merged wholesale into `MapGeneIdsFromGTF`; only the narrow, evidence-backed operation of reconciling terminal Ensembl numeric versions belongs there. The GTF matcher must try exact identifiers first and may then remove only a terminal numeric version while retaining `_PAR_Y`, for example with `re.sub(r"\.\d+(?=_PAR_Y$|$)", "", identifier)`.

Deletion must wait for these migration prerequisites:

1. Inventory workflows and consumers that currently rely on rewritten `var_names`, including annotation tools that expect symbols.
2. Move Ensembl-version reconciliation into `MapGeneIdsFromGTF`, where the input identifier field, GTF file, match mode, unmatched count, and annotation release are reportable.
3. Store gene symbols as an explicit variable annotation such as `.var["gene_symbols"]`; keep stable gene IDs as the principal identity.
4. For packages that insist on symbols in `var_names`, construct a private working copy inside that package adapter, validate duplicate symbols there, and map results back by stable ID.
5. Migrate example and saved workflows. Add tests that no supported workflow uses this node and that all symbol-dependent consumers declare their identifier source.
6. Ensure no migration changes current `.var` without either occurring before Raw snapshot creation or deliberately preserving/synchronizing the Raw snapshot namespace.
7. Remove registration only after the GTF and consumer migrations pass; document arbitrary regex replacement as an unsupported data-cleaning step rather than silently emulating it.

## Compatibility-period contract

If one release of compatibility is required, relabel the node `Deprecated: Replace Feature Names` and describe it as data cleaning, not normalization. Preserve old positional inputs, append `summary` and `code`, reject execution when `.raw` already exists, save the original index in a collision-safe `.var` column, and change the default collision policy to `error` for new workflows. A legacy `make_unique` path may remain only as an explicit compatibility choice.

The report must include changed-name count, unchanged count, collision count, artificial-suffix count, input/output uniqueness, pattern/replacement, identifier-system limitation, a deprecation warning, relevant software documentation, and dynamically collected Python/openbio-singlecell/AnnData/Pandas versions. It must not cite HGNC/GENCODE as if an arbitrary regex implemented those resources.

## Parameter policy

- Target state after migration: no standalone parameters because the node is deleted.
- GTF module: expose identifier source and gene-symbol destination; keep the precise version regex hidden, fixed, tested, and disclosed.
- Symbol-dependent consumer: expose the annotation column only when it affects the scientific input; hide temporary index adaptation.
- Never expose a general `make_unique` switch as resolution of biological identity ambiguity.

## Cohesion and coupling assessment

The current module is a shallow pass-through from regex parameters to pandas. It couples every downstream analysis to a mutated global feature axis and duplicates identifier knowledge that belongs beside the annotation resource. The higher-cohesion design is to delete the node, localize stable-ID reconciliation in the GTF module, and localize package-specific symbol adaptation in each consumer.

## Verification plan

- Registry/workflow tests prove no supported flow depends on the node before removal.
- GTF tests cover exact IDs, versioned IDs, and `_PAR_Y` preservation without generic regex mutation.
- Symbol-dependent consumer tests use `.var["gene_symbols"]` or a private working index and map results back by stable identity.
- If the deprecated compatibility node remains, test literal/regex replacement, collisions, input immutability, Raw snapshot rejection, original-index preservation, strict-JSON summary, and compiling/equivalent code.

## Open expert-boundary decision (2026-08-28)

Revise the compatibility contract: Raw presence, duplicate output, empty axes, and an empty pattern become warnings or
ordinary declared behavior. Hard errors remain for backed mutation, invalid regular-expression syntax, and collision
with the reserved provenance column. Raw is an independent AnnData snapshot; the node records that its namespace was
preserved rather than asserting a binding to current `var_names`. No Raw-binding test is retained.
