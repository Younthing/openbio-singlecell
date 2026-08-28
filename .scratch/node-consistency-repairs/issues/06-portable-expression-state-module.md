# Design one portable expression-state Module

Type: task
Blocked by: 05

## Question

What single state vocabulary and dependency-free resolver can serve package runtime callers and be embedded verbatim
in standalone `code` outputs without confusing a log transform with verified normalization provenance?

## Acceptance criteria

- `logged` versus `logged_unverified` has one documented scientific meaning.
- One pure resolver owns operation-to-state transitions; package wrappers and emitted source reuse that implementation.
- Raw remains an ordinary explicit source selection; do not add provenance-binding or history-tamper gates.
- Runtime/generated contract tests cover supported state transitions without adding Raw-binding or history-tamper edge
  cases.
