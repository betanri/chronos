# Chronos

This repository is organized into three main sections:

## 1) [WHY CHRONOS AND NOT TREEPL](1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

A 720-run benchmark showing, in the same simulated conditions, that chronos had lower dating error and fewer failures than treePL (conditions tested: 4 clock regimes x 3 extinction levels x 2 heterotachy levels x 30 replicates). It is there to justify method choice with data.

## 2) [CHRONOS -- CUSTOM DATING TREE PIPELINE](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md)

This tab is the workflow: a practical wrapper around `ape::chronos` that automates things ape alone does not give you out-of-the-box in one run. It handles calibration setup, clock-model fitting, lambda tuning, and the production of dated-tree outputs.

## 3) [POST-FIT EVALUATION METRICS](3_POST_FIT_EVALUATION_METRICS/README.md)

This tab is the comparison layer applied after fitting is finished. It focuses on how the resulting chronograms behave biologically under `pulse preservation`, `gap burden`, and `rate plausibility`, using the Terap example as the worked case.
