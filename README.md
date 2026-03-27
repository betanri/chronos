# Chronos

This repository now centers on a **multi-benchmark comparison of `chronos`, `treePL`, and `RelTime`**, paired with a compact `treePL` environment diagnostic and a companion post-fit evaluation layer in `PhyloChronoRank (PCR)`.

## Main benchmark suite

## [BENCHMARK SUITE: chronos, treePL, and RelTime](WHY_CHRONOS_AND_NOT_TREEPL_OR_RELTIME/README.md)

This is the active benchmark write-up. It now reports finished `A-E` results and lays out the broader suite structure with:

- `A-E` main benchmarks
- `RelTime` treated as a full equal to `chronos` and `treePL`
- future `P1/P2` Paradis-style sequence-simulation extensions
- compact `treePL` environment diagnostics (`Tenv`)
- a linked `PCR` post-fit interpretation layer

The benchmark page is now the main place where we summarize:

- method accuracy by benchmark
- `chronos` model-specific behavior
- `RelTime` performance on equal footing with `chronos` and `treePL`
- clock-model recovery
- `treePL` environment sensitivity
- how the empirical `PCR` layer connects back to the simulation benchmark

## Shared dating workflow

## [SHARED DATING GRID PIPELINE: PhyloChronoRank (PCR)](https://github.com/betanri/PhyloChronoRank/blob/main/1_PCR_CUSTOM_DATING_PIPELINE_FROM_PHYLOGRAMS/README.md)

This is the shared fitting workflow used after calibration inputs have been defined. It runs `chronos`, `treePL`, and `RelTime` from one calibration source, either from a calibration table or from a reference time tree by congruification.

## Post-fit evaluation layer

## [POST-FIT EVALUATION METRICS: PhyloChronoRank (PCR)](https://github.com/betanri/PhyloChronoRank/blob/main/2_PCR_POSTFIT_METRICS/README.md)

`PCR` is the post-fit comparison layer. It evaluates dated trees after fitting is already done, using metrics such as:

- pulse preservation
- gap burden
- rate irregularity

In plain terms, the benchmark in this repo asks **which method dates better under controlled simulations**, while `PCR` asks **which fitted chronogram behaves more plausibly once it already exists**.
