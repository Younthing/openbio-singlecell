# Preview Result: module design review

## Current module

The node accepts one of the three concrete OpenBio report artifacts and returns a ComfyUI UI
payload. Summary and table values stay in memory. Plot bytes are written to a deterministic
temporary path so ComfyUI can display the image.

## Correctness and cohesion findings

- The plot path is written directly. An exception or process interruption can expose a
  truncated preview and overwrite the previous valid preview.
- `PlotResult` currently proves only the eight-byte PNG signature. The adapter must reject a
  truncated or corrupt payload before publishing it.
- Workflow and node IDs are framework plumbing and should remain hidden; they are not user
  parameters.
- Previewing is one coherent UI responsibility. It should not also export durable files or
  change the supplied report artifact.

## Decision: enhance and keep separate

Keep this node as a terminal preview adapter. Do not merge it with Save PNG: preview files
have deterministic temporary identity and replacement semantics, whereas durable exports
use user-selected names and collision policy. Do not merge it with analysis nodes because
that would couple scientific computation to UI execution.

This is not an analysis/processing node. It intentionally has no graph outputs and does not
add `summary` or `code`; the supplied report already contains the scientific disclosure.

## Target interface

- Visible input: typed `result` (`summary`, `table`, or `plot`).
- Hidden inputs: ComfyUI workflow metadata and node identity.
- Graph outputs: none; UI payload only.
- Hidden fixed behavior: payload normalization, deterministic identity hashing, PNG
  verification, same-directory staging, flush/fsync, atomic replacement, and cleanup.

## Invariants and tests

Reject values outside the three report contracts. For a plot, verify and decode the PNG
before any destination change. Publish only a complete staged file and preserve the prior
preview if validation, staging, or replacement fails. Tests should cover workflow isolation,
stable reruns, real PNG decoding, corrupt/truncated data, cleanup, and failure atomicity.
