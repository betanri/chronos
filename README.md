# Chronos

This repository is organized into two main sections:

## 1) [WHY CHRONOS AND NOT TREEPL](1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

A 720-run benchmark showing, in the same simulated conditions, that chronos had lower dating error and fewer failures than treePL (conditions tested: 4 clock regimes x 3 extinction levels x 2 heterotachy levels x 30 replicates). It is there to justify method choice with data.

## 2) [CHRONOS -- CUSTOM DATING TREE PIPELINE](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md)

This tab is the workflow: a practical wrapper around `ape::chronos` that automates things ape alone does not give you out-of-the-box in one run. It keeps three layers distinct: clock-model fitting, lambda tuning, and post-fit evaluation metrics. Those post-fit metrics are not part of clock fitting or lambda tuning; they are a separate evaluation layer applied to the resulting chronograms.
