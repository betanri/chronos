# Chronos

This repository is organized into two main sections:

## 1) [WHY CHRONOS AND NOT TREEPL](1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

A 720-run benchmark showing, in the same simulated conditions, that chronos had lower dating error and fewer failures than treePL (conditions tested: 4 clock regimes x 3 extinction levels x 2 heterotachy levels x 30 replicates). It is there to justify method choice with data.

This section now also incorporates the newer comparative metric framework used later in the project:

- `pulse preservation` is the main diversification-rhythm diagnostic
- `gap burden` captures fossil-gap burden or calibration slack when that information exists
- `rate plausibility` captures how extreme the implied branchwise rate shifts are

For the 720 simulation benchmark specifically, the additional interpretation is:

- `pulse preservation`: chronos still outperforms treePL on the representative-tree comparison
- `gap burden`: not applicable here because the benchmark uses only a fixed root calibration
- `rate plausibility`: chronos also outperforms treePL on the representative-tree comparison

- Folder: [`1_WHY_CHRONOS_AND_NOT_TREEPL`](1_WHY_CHRONOS_AND_NOT_TREEPL)
- Start here: [`1_WHY_CHRONOS_AND_NOT_TREEPL/README.md`](1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

## 2) [CHRONOS -- CUSTOM DATING TREE PIPELINE](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md)

This tab is the workflow: a practical wrapper around `ape::chronos` that automates things ape alone does not give you out-of-the-box in one run, including calibration handling (manual or congruification), model/lambda grid search, robust model selection, threshold sensitivity checks, pulse-preservation diagnostics, optional large-tree subset tuning, and clean reproducible outputs.

The main update relative to the older version is a broader evaluation layer. The pipeline now reads candidate trees through four complementary lenses:

- model-fit preference
- pulse preservation
- gap burden
- rate plausibility

- Folder: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE)
- Main script: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/Run_chronos_pipeline.R`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/Run_chronos_pipeline.R)
- Pipeline guide: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md)
- Example files: [`2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES`](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES)
