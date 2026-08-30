# File artifacts, ArtifactTickets, and one-shot Workers

Implement the approved workflow-wide redesign from fixed point `d553573ffec73a7eb4d1e4eabd1326550ba443d5`: large values cross node boundaries only as immutable process-local ArtifactTickets, every scientific node invocation executes in a fresh one-shot Worker, all File artifacts from one run publish atomically under a shared RunLease, portable File artifacts can become Persisted artifacts, and ComfyUI's workflow signature remains the cache key.

The RunLease must keep the shared run directory available while any cache entry or running consumer retains one of its ArtifactTickets, then own cleanup after the final reference is released. Cache hits reuse tickets without starting a Worker; Classic, None, and LRU are supported, while the default RAM-pressure cache is rejected because ticket size cannot account for retained files.

The registry must switch atomically with no live-object compatibility mode. Native session-only artifacts are limited to the originating runtime session and producer's probed Python interpreter identity and cannot become Persisted artifacts or be restored. Persist Artifact is the only generic durability boundary for portable File artifacts; it returns no workflow value, and there is no generic Restore node.

The accepted codecs, Worker trust boundary, node-by-node copy audit, tests, acceptance criteria, and explicit exclusions remain those in the originating user plan. Preserve algorithm-required scientific copies and Raw snapshot semantics, and leave the pre-existing `AGENTS.md` worktree modification untouched.
