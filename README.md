# Chronos

This repository is organized into two main sections:

## 1) WHY CHRONOS AND NOT TREEPL

This tab is the evidence: a 720-run benchmark showing, in the same simulated conditions, that chronos had lower dating error and fewer failures than treePL (conditions tested: 4 clock regimes × 3 extinction levels × 2 heterotachy levels × 30 replicates). It is there to justify method choice with data, not opinion.

- Folder: `WHY_CHRONOS_AND_NOT_TREEPL`
- Start here: `WHY_CHRONOS_AND_NOT_TREEPL/README.md`

## 2) CHRONOS -- CUSTOM DATING TREE PIPELINE

This tab is the workflow: a practical wrapper around `ape::chronos` that automates things ape alone does not give you out-of-the-box in one run, including calibration handling (manual or congruification), model/lambda grid search, robust model selection, threshold sensitivity checks, branching-tempo diagnostics, optional large-tree subset tuning, and clean reproducible outputs.

- Folder: `CHRONOS_CUSTOM_DATING_TREE_PIPELINE`
- Main script: `CHRONOS_CUSTOM_DATING_TREE_PIPELINE/Run_chronos_pipeline.R`
- Pipeline guide: `CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md`
- Example files: `CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES`
