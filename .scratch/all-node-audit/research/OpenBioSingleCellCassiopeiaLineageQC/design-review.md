# Cassiopeia character preparation and lineage QC: module design review

## Decision

**Retain and deepen `OpenBioSingleCellCassiopeiaLineageQC` as the sole allele-table-to-character-artifact module.** Its interface is one bounded file plus declared protocol semantics; its implementation owns parsing, conflict checks, empirical priors, public conversion, cell/character QC, canonicalization and tamper evidence. This is a deep module: deleting it would force every solver to duplicate scientifically coupled conversion rules. Tree solving remains downstream because it is a separate inference operation.

The module creates a typed process-local `OPENBIO_CASSIOPEIA_CHARACTERS` artifact. QC thresholds are expert screening policies: a lineage that remains computationally valid is retained with `status="warning"` and disclosed flags. Only a lineage for which the backend produces no usable character payload is `unavailable` downstream.

## Current interface and invariants

Inputs, in order:

1. `allele_table_file`;
2. `first_column_as_index` (file-format adapter only);
3. exact `lineage_column`, `cell_barcode_column`, `integration_barcode_column`;
4. explicit comma-separated `cut_site_columns` (no silent cassette discovery);
5. explicit comma-separated `prior_grouping_columns`;
6. `missing_data_allele` (blank means true CSV nulls only; never the uncut token);
7. `allele_representation_threshold=0.98`;
8. `minimum_cells=2`, `maximum_missing_fraction=0.8`, `maximum_uncut_fraction=0.8`, `minimum_unique_fraction=0.05`, and `minimum_informative_character_fraction=0.2`;
9. advanced `max_file_mib=512` and `max_matrix_gib=2` resource limits.

Outputs are `characters`, `table`, `summary`, and `code`.

The implementation reads and hashes the exact same bounded bytes once, canonicalizes identifiers and row order, normalizes only declared columns to Cassiopeia's public names, rejects global duplicate cell identity and per-character allele conflicts, and invokes only `cassiopeia.pp.compute_empirical_indel_priors` and `cassiopeia.pp.convert_alleletable_to_character_matrix`. It verifies backend non-mutation, integer state semantics (`-1`, `0`, positive mutations), exact axes, observed-state denominators, prior/state-map closure and finite resource estimates. All-missing/all-uncut rows follow the declared fraction filters rather than an extra hidden gate. Missing lineage/intBC prior grouping or cut-site grouping is honored with a prominent expert-review warning. QC pass/warning is separate from computational availability.

The immutable artifact owns defensive matrix/prior/state-map copies and exact ordered cell/character axes. Its metadata binds the pinned distribution/tag/sdist, source bytes, column/cut-site schema, prior population and fingerprint, QC policies, per-lineage content hashes, QC-table hash and one canonical artifact fingerprint. Consumers validate current contents against all hashes before receiving a copy.

`summary` uses the strict report schema (`methods`, report-ready `results`, `key_results`, `references`, `software_versions`, parameters, warnings and limitations) with no NaN/Infinity. `code` is standalone source that accepts a file path and returns `(characters_artifact, qc_table, summary_dict)` with the same validation and numerical results.

## Legacy/current migration matrix

Legacy output `table` becomes current output index 1; output index 0 is the new typed artifact. Legacy `tumor_column` maps by name to `lineage_column`; `minimum_cells_for_summary`, `maximum_uncut_fraction`, `allele_representation_threshold`, and `percent_unique_threshold` have direct value mappings. The following legacy fields are **not** semantically equivalent and make automatic migration fail closed:

- `cut_sites_per_intbc`: a count cannot prove exact cut-site column identities;
- `minimum_intbc_fraction`: depended on a private filtering function removed from the audited module;
- `lineage_size_threshold`: mixed descriptive size with reconstruction eligibility;
- `percent_unsaturated_threshold`: used an incorrect hard-coded cassette denominator;
- absence of missingness, conflict, prior-grouping and resource policies.

Migration may proceed only when a workflow-owned review payload explicitly supplies exact cut sites, prior grouping, missing allele and all new policies. Otherwise the graph is left untouched with a rerun/review error. Migration must be atomic and idempotent and must preserve array/object links, subgraphs, exposures and global IDs. Every serialized input port, including a widget converted to a connected input, must retain the reviewed widget type; matching names alone are insufficient evidence of the current or legacy schema.

## Verification gate

Tests cover CSV/TSV bytes and TOCTOU identity, exact columns, non-three-site cassettes, prior grouping, duplicate cells across lineages, identical duplicate rows versus conflicts, reserved uncut/missing tokens, all-missing/all-uncut cells, denominator fixtures, failed-lineage accounting, deterministic row/state order, backend mutation, artifact tampering and defensive copies, memory limits, exact package/signature gating, strict JSON, compiled standalone code and numerical/structural parity. A real smoke uses `cassiopeia-mt==2.1.3`.
