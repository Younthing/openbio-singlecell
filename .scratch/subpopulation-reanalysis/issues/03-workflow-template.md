# Subpopulation workflow template

Type: task
Status: resolved

Add and verify a complete explicit subpopulation reclustering, evidence, annotation, and save workflow.

## Comments

- Do not introduce a monolithic Subpopulation Analysis node.

## Answer

Added the ordinary-node `Subpopulation Reclustering and Annotation` workflow, its cover and release entries. The
template restores Raw, recreates explicit expression snapshots, recomputes full subpopulation geometry and marker
evidence, requires manual population/mapping choices, saves the subpopulation, and transfers only its reviewed
subtype annotation to the parent.
