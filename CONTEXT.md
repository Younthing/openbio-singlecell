# OpenBio Single-Cell Analysis

Shared language for the scientific objects and claims represented by OpenBio single-cell workflows.

## Language

**Sample**:
An independent biological specimen used as the replicate unit for quality assessment and condition-level inference.
_Avoid_: Cell batch, group

**Technical batch**:
A technical source of unwanted variation that integration may model without removing the biological condition of interest.
_Avoid_: Condition, sample

**Condition**:
The biological cohort, disease state, or treatment whose differences are the subject of inference.
_Avoid_: Batch

**Raw snapshot**:
The post-QC, full-gene count state retained before highly variable gene selection so the accepted counts can be recovered.
_Avoid_: Unfiltered input, log-normalized data

**Cluster marker evidence**:
Ranked expression features that help a reviewer annotate a discovered cell cluster; it is exploratory evidence, not confirmatory condition-level inference.
_Avoid_: Differential condition result, cell-type truth

**Provisional annotation**:
An automated reference-model label used to keep an analysis flow complete before expert review; it remains visibly separate from the curated annotation.
_Avoid_: Final cell type, ground truth

**Curated annotation**:
An expert-reviewed cell-population label used for formal population-level interpretation and inference.
_Avoid_: Provisional annotation, classifier prediction

**Condition contrast**:
A replicate-aware comparison between biological conditions within a defined cell population.
_Avoid_: Cell-level cluster marker test
