# Quality-control node audit and refactor

Type: task
Status: resolved

## Scope

- `OpenBioSingleCellCalculateQC`
- `OpenBioSingleCellFilterCells`
- `OpenBioSingleCellFilterGenes`
- `OpenBioSingleCellQCPlots`

## Acceptance criteria

- The two required pre-change documents exist for every node.
- Count-based calculations disclose their expression source and do not silently describe transformed values as counts.
- Combined filtering remains an atomic intersection of the configured criteria.
- Every node returns a report-ready structured `summary` and equivalent-function `code` in addition to its primary result.
- Summaries disclose input/output dimensions, retained/removed observations or variables, QC distributions, configured thresholds, warnings, references, and software versions as applicable.
- Tests cover schemas, JSON-serializability, report contents, code contents, sparse matrices, input immutability, and scientific edge cases.

## Baseline

The repository-local `.venv` cannot collect ComfyUI-dependent tests because it lacks `torch`. The configured ComfyUI environment is valid: `D:\learn\ComfyUI\.venv\Scripts\python.exe -m pytest -q` passed all 138 tests before this batch.

## Comments

- 2026-08-28: Claimed as the first batch because the four modules share raw-count and QC-metric semantics while retaining atomic external interfaces.
- 2026-08-28: Independent standards/spec review found and drove fixes for legacy widget ordering, sparse explicit-zero and duplicate-entry semantics, mixed mitochondrial denominators, unavailable mitochondrial evidence, non-finite plotting data, and generated-code error equivalence.

## Answer

Enhanced all four nodes without merging their atomic responsibilities. Nodes now expose and report `X`/`raw`/layer selection, preserve the caller's sparse storage, and hard-fail only structural or mathematically invalid inputs. Signed/non-count-like expert selections and contradictory/all-removed filtering outcomes remain executable with explicit warnings and result disclosure. QC feature sets with no matching genes are reported as unavailable rather than biological zero; plots use coherent metric sources and omit unavailable/non-finite mitochondrial panels.

Each node now returns its primary result, an `OPENBIO_SINGLE_CELL_SUMMARY` artifact with a strict JSON-compatible `.summary` payload, and compiled equivalent-function source on the `code` port. The shared reporting seam records Methods, Results, key results, resolved parameters, references, limitations/warnings, and runtime software versions. Generated functions reproduce scientific validation and primary results while intentionally omitting plugin-only history bookkeeping.

Validation: `153 passed` for the full suite; Ruff passed; the independent final review reported no remaining P0/P1 findings.
