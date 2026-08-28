# OpenBioSingleCellDrugScores — design review

## Final open-expert boundary

Expression-state history is a reporting input, not an authority gate. Any finite explicitly selected source remains
executable, with nonrecommended scale disclosed. No history-tamper or Raw/current binding contract is introduced.

## Decision

**Enhance and narrow to one-drug descriptive target-expression score.** It consumes the explicit artifact from the
DGIdb resource loader. Do not merge it with drug ORA/GSEA: scoring cells, testing an explicit selected set, and
testing a complete rank are different scientific operations.

## Atomic interface

Inputs, in order:

1. `adata`;
2. typed DGIdb `resource` artifact;
3. required exact `drug` identifier;
4. explicit expression `source` (`X`, named layer, or Raw snapshot), with log-normalized expression recommended and
   Raw/count-like selection disclosed as an expert choice;
5. nonblank `output_key`, default `drug_target_score`;
6. `min_matched_targets`, default 1;
7. explicit overwrite policy.

Outputs:

```text
adata, summary, code
```

Fix Pertpy `method="mean"`, `nested=False`, and one target dictionary. No group, Condition, Sample, Technical batch,
p-value, control-size/bin, random seed, direction, dose, or online resource selector is exposed. The internal adapter
uses a private copy, extracts Pertpy's one-column score, and stores a finite Series in `output.obs[output_key]`;
backend `uns` scratch keys do not leak unless intentionally namespaced and documented.

## Validation/provenance

Validate resource producer/schema/content fingerprint, exact human HGNC namespace, resource release/license/hash,
drug identity, target uniqueness, exact expression feature IDs, finite aligned numeric data, target overlap,
observation
order, and output collision. Report requested/matched/unmatched targets and contributing DGIdb sources. Inputs are
immutable and failure is atomic.

`summary` follows the strict analysis-report schema with report-ready method/results, bounded score extremes and
distribution, target/resource accounting, non-inferential and non-clinical limitations, citations, warnings, and
dynamic OpenBio/Pertpy/AnnData/Scanpy/NumPy/pandas/SciPy versions. `code` defines an equivalent function accepting
the resource DataFrame/metadata and returning `(AnnData, summary_dict)` with exact summary equality.

## Scientific boundary

Per-cell scores characterize one target set only. Sample remains the inference unit for downstream Condition
contrasts and Technical batch a nuisance variable; neither belongs inside the score. A high score is not drug
sensitivity or a recommendation.

## Migration and tests

Legacy workflows add a resource link and drug choice; scoring “all drugs” cannot be silently preserved. The old node
usually fails before creating output, so migration should stop with a clear selection requirement rather than emulate
the broken behavior.

Tests: fresh DGIdb initialization regression; exact Pertpy 1.3 call/storage extraction; one/multiple/no targets;
resource artifact/hash tampering; namespace/release/license validation; normalized and explicit Raw/count-like source
warnings; all-zero observations with finite zero scores; finite score hand fixture; sparse/dense equality; output
collision; no network; immutability/atomic failure; strict JSON;
code compilation; and runtime/generated score+summary equality.

## Cohesion assessment

The module owns one score for one declared drug target set. Resource lifecycle is upstream and statistical
interpretation downstream, while Pertpy storage quirks and validation remain hidden implementation details.
