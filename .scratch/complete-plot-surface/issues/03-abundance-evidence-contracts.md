# Retain Milo and composition-model diagnostics

Type: task
Status: resolved

## Scope

- Portable Milo result with membership, graph, coordinates, axes, fingerprints, and canonical result table.
- Portable scCODA/tascCODA model result with canonical table, hierarchy/model metadata, and validated InferenceData/posterior diagnostics.

## Acceptance criteria

- New typed artifacts have strict codecs, content/axis validation, persistence policy, Worker round-trip tests, and real downstream Plot consumers.
- Existing table outputs remain current primary evidence and are not recomputed by consumers.

## Comments

- 2026-09-04: Claimed because the user explicitly authorized expanding these producer contracts.
- 2026-09-04: Resolved with vertical red-green slices; downstream Plot consumption was handed to ticket 07.

## Answer

- Milo now publishes an immutable portable result containing the unchanged canonical table, exact observation and
  neighborhood axes, binary membership, the neighborhood-overlap graph, selected-representation index-cell
  coordinates, semantic provenance, and current-content fingerprints. Its H5AD/table codec rejects schema, axis,
  graph, coordinate, and payload tampering.
- scCODA and tascCODA now publish an immutable portable result containing the unchanged canonical table, complete
  focal posterior draws, per-draw NUTS energy/step diagnostics, method/model/annotation metadata, and the complete
  tascCODA hierarchy manifest. The codec uses the installed ArviZ/xarray Zarr representation and verifies exact
  groups, variables, dimensions, coordinates, and fingerprints on read.
- Both current producer schemas return `result, table, summary, code`; both artifact kinds are accepted by Persist
  Artifact. Codec, tamper, persistence, generated-code parity, and Worker round-trip tests pass without changing the
  input artifacts.
