# Chronos

This repository is organized into two main sections:

## 1) [WHY CHRONOS AND NOT TREEPL](1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

This tab is the evidence: a 720-run benchmark showing, in the same simulated conditions, that chronos had lower dating error and fewer failures than treePL (conditions tested: 4 clock regimes × 3 extinction levels × 2 heterotachy levels × 30 replicates). It is there to justify method choice with data.

- Folder: [`1_WHY_CHRONOS_AND_NOT_TREEPL`](1_WHY_CHRONOS_AND_NOT_TREEPL)
- Start here: [`1_WHY_CHRONOS_AND_NOT_TREEPL/README.md`](1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

## 2) [CHRONOS -- CUSTOM DATING TREE PIPELINE](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md)

This tab is the workflow: a practical wrapper around `ape::chronos` that automates things ape alone does not give you out-of-the-box in one run, including calibration handling (manual or congruification), model/lambda grid search, robust model selection, threshold sensitivity checks, branching-tempo diagnostics, optional large-tree subset tuning, and clean reproducible outputs.

- Folder: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE)
- Main script: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/Run_chronos_pipeline.R`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/Run_chronos_pipeline.R)
- Pipeline guide: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md)
- Example files: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES)
