# Monocle 2 dependency API adapter

## Request and scope

The user approved trying a small R adapter for the removed igraph and dplyr APIs used by Monocle 2.
Keep the Monocle 2 algorithms, native scientific parameters, expert choices, node schemas and file protocol.
Adapt the affected function calls in a private R environment. Do not modify installed packages or their global
namespaces, silently retry failed scientific analyses, introduce version hard gates, or add a generic patch framework.
This supersedes the earlier subpopulation spec's exclusion of local dependency adaptation.

## Implementation and public seams

- Add `openbio_singlecell/r/monocle2_compat.R`, loaded by the existing native R driver.
- Cover the real removed calls in dispersion estimation, DDRTree ordering/projection, and BEAM.
- Include the adapter in standalone reproduction code, shipped R resources, cache identity and execution reports.
- Use the existing `run_analysis`, emitted reproduction code and native R driver as the behavior test seams.
- Create a separate modern-dependency R environment; retain the previously verified native environment as reference.
- Compare graph structure, root and parent assignments, State/Pseudotime, dispersions, DE and BEAM results.
  Report numerical differences or remaining limitations instead of treating successful execution as equivalence.
- Update the installation guidance only after real verification. Preserve unrelated worktree edits; do not commit.

## Fixed point and pre-existing work

HEAD: `57eb64e48ec853bb231c36c9c6a769dd25a453c5`.
The complete uncommitted Monocle 2 feature and unrelated cache-pressure changes already existed at turn start.
Their file contents were saved for incremental review at
`C:/Users/admin/AppData/Local/Temp/openbio-monocle2-compat-before-r36x5oks` (`snapshot.json` lists the files).
Review this adapter increment against that snapshot and the complete affected Monocle feature against HEAD;
exclude the unrelated cache policy, ADR and architecture edits from this request.

## Verification

Use vertical red-green behavior tests, then focused Python/R regressions, emitted-code execution in a fresh R
process, actual wheel contents, and independent Standards and Spec review. Record commands and results in
`verification.md`; environment details belong in `runtime/README.md`.
