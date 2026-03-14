# Chronos

This repository centers on a chronos-vs-treePL benchmark and points to the companion `PhyloChronoRank (PCR)` repo for shared dating-grid generation and post-fit chronogram comparison.

## [WHY CHRONOS AND NOT TREEPL](WHY_CHRONOS_AND_NOT_TREEPL/README.md)

A 720-run benchmark showing, under an exact-root and calibration-sparse simulation design, that chronos had lower dating error and fewer failures than treePL (conditions tested: 4 clock regimes x 3 extinction levels x 2 heterotachy levels x 30 replicates). It is there to justify method choice with data, but it should be read as an idealized benchmark rather than as a full proxy for empirical dating problems with uncertain root ages and multiple internal calibrations.

## [SHARED DATING GRID PIPELINE: PhyloChronoRank (PCR)](https://github.com/betanri/PhyloChronoRank/blob/main/README_Dating_Grid.md)

This guide runs `chronos`, `treePL`, and `RelTime` from one shared calibration source, using either a calibration CSV or a reference backbone time tree for congruification. It is the active dating workflow paired with the chronos benchmark in this repo.

## [POST-FIT EVALUATION METRICS: PhyloChronoRank (PCR)](https://github.com/betanri/PhyloChronoRank)

This companion repo, `PhyloChronoRank (PCR)`, is the comparison layer applied after fitting is finished. It asks a different question from model fitting: once you already have dated trees, which one behaves most plausibly as a chronogram?

It evaluates trees under `pulse preservation`, `gap burden`, and `rate irregularity`, and shows how those metrics can be read across two empirical datasets. In plain language, it turns visual impressions about branching bursts, calibration fit, and implied rate behavior into explicit side-by-side comparisons.
