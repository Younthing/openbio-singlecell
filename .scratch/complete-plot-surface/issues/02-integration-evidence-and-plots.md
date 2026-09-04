# Retain integration diagnostics and plot them

Type: task
Status: resolved

## Scope

- Harmony complete objective/convergence evidence and Harmony Convergence Plot.
- scVI complete epoch metric evidence in the native model artifact and scVI Training Plot.

## Acceptance criteria

- Diagnostic evidence round-trips with exact identity and tamper checks.
- Plots never retrain/reintegrate and preserve native Worker affinity.
- Runtime, summary, and generated code remain equivalent.

## Comments

- 2026-09-04: Claimed as an upstream-contract prerequisite.
- 2026-09-04: Resolved with focused red-green slices; 97 integration/scVI tests and Ruff pass. The full suite
  remained temporarily blocked during parallel work by an unrelated unimplemented Schist plot test.

## Answer

- Harmony now retains and validates the complete public harmonypy objective contract (`objective_harmony`,
  `objective_kmeans`, and `kmeans_rounds`) with outer-boundary consistency plus source/adjusted basis,
  observation-axis, and Technical-batch fingerprints. Harmony Convergence Plot reads only this evidence.
- scVI now retains all six base native epoch-indexed train/validation metric series plus any exposed global-KL
  train/validation series in a strict JSON sidecar inside the Worker-affine native model
  artifact. The sidecar fingerprints its payload and `model.pt`; runtime and standalone scVI Training Plot code read
  only that verified artifact evidence, including after native save/load.
- Runtime/generated PNG parity, codec round trips, input/global-state immutability, tamper rejection, and static
  plot resource limits are covered by focused tests.
