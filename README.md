# Chronos

This repository is organized into two main sections inside `Chronos` plus one companion repo:

## 1) [WHY CHRONOS AND NOT TREEPL](1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

A 720-run benchmark showing, in the same simulated conditions, that chronos had lower dating error and fewer failures than treePL (conditions tested: 4 clock regimes x 3 extinction levels x 2 heterotachy levels x 30 replicates). It is there to justify method choice with data.

## 2) [CHRONOS -- CUSTOM DATING TREE PIPELINE](2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/README.md)

This tab is the workflow: a practical wrapper around `ape::chronos` that automates things ape alone does not give you out-of-the-box in one run. It handles calibration setup, clock-model fitting, lambda tuning, and the production of dated-tree outputs. The post-fit evaluation layer now lives in the companion repo listed below.

## 3) [POST-FIT EVALUATION METRICS: PhyloChronoRank (PCR)](https://github.com/betanri/PhyloChronoRank)

This companion repo, `PhyloChronoRank (PCR)`, is the comparison layer applied after fitting is finished. It asks a different question from model fitting: once you already have dated trees, which one behaves most plausibly as a chronogram?

It evaluates trees under `pulse preservation`, `gap burden`, and `rate plausibility`, and shows how those metrics can be read across two empirical datasets. In plain language, it turns visual impressions about branching bursts, calibration fit, and implied rate behavior into explicit side-by-side comparisons.
