# Chronos Empirical Pipeline

This pipeline dates a phylogram with `ape::chronos` and gives you:

- one fit-based selected chronogram
- one chronogram per clock model (`clock`, `correlated`, `relaxed`, `discrete`)
- an explicit comparison between:
  - model-fit preference
  - branching-tempo similarity to the original phylogram

The goal is to help you decide which dated tree is most defensible for your biological question.

## What The Pipeline Does

Given a target phylogram, the script:

1. Reads the tree (single tree, multi-tree, or named-Newick).
2. Builds calibrations either:
   - from a reference timetree (congruification), or
   - from a manual CSV.
3. Fits chronos across clock models and lambda values.
4. Applies robust model selection with threshold sensitivity (`1` and `2`).
5. Writes:
   - threshold-selected trees
   - one tree per model
   - fit tables
   - branching-tempo metric table
   - a plain-language interpretation text file.

## Main Script

- `Run_chronos_pipeline.R`

## Required R Packages

```r
install.packages(c("ape", "geiger", "phytools"))
```

## Inputs

### Option A: Reference-tree calibrations (congruification)

Set:

```r
USE_CONGRUIF <- TRUE
TARGET_TREE_FILE <- "path/to/target_phylogram.tre"
REFERENCE_TIME_TREE <- "path/to/reference_time_tree.tre"
```

### Option B: Manual calibrations (CSV)

Set:

```r
USE_CONGRUIF <- FALSE
TARGET_TREE_FILE <- "path/to/target_phylogram.tre"
MANUAL_CAL_CSV <- "path/to/manual_calibrations.csv"
```

Manual CSV must contain:

- `taxonA`
- `taxonB`
- `age_min`
- `age_max`

## How To Run

```r
setwd("/Users/ricardobetancur/Desktop/Proxy_Misplaced/chronosPL/Terapontoid_Trees")
source("Run_chronos_pipeline.R")
```

## Model-Fit Sensitivity

The script evaluates and reports both:

- `PLOG_CLOCK_SWITCH_THRESH = 1` (default strict)
- `PLOG_CLOCK_SWITCH_THRESH = 2` (stricter)

It does **not** use threshold `0` in the default empirical workflow.

## Branching-Tempo Metric (Biology-Oriented Diagnostic)

This metric compares each chronogram to the original phylogram in terms of branching tempo.

How it works:

1. Match internal nodes by clade identity (same descendant tip set).
2. Normalize node heights to `[0,1]` (root deep, tips shallow).
3. Compute:
   - `tempo_mae_all` (all internal nodes)
   - `tempo_mae_early_q75` (deep/early nodes only, top quartile)
   - `tempo_median_early_q75`

Lower values mean the chronogram preserves the phylogram’s branching-tempo pattern better.

## Outputs (Dedicated Run Folder)

Outputs are written to:

- `<OUT_BASE_DIR>/<OUT_PREFIX>/`

with subfolders:

- `tables/`
- `trees/`
- `logs/`
- `checkpoints/`

Main files:

- `tables/summary_<OUT_PREFIX>.csv`
  - fit-favored model (default threshold row)
  - whether fit and tempo diagnostics agree
- `tables/summary_<OUT_PREFIX>_sensitivity.csv`
  - selected model for thresholds `1` and `2`
- `tables/summary_<OUT_PREFIX>_model_fits.csv`
  - per-model fit stats + branching-tempo metrics
- `tables/interpretation_<OUT_PREFIX>.txt`
  - plain-language conclusion:
    - fit-favored model
    - lowest overall tempo-error model
    - lowest early-tempo-error model
    - agreement/disagreement statement
- `trees/<target_id>_chronos_dated_clockthresh1.tre`
- `trees/<target_id>_chronos_dated_clockthresh2.tre`
- `trees/<target_id>_chronos_dated_modelclock.tre`
- `trees/<target_id>_chronos_dated_modelcorrelated.tre`
- `trees/<target_id>_chronos_dated_modelrelaxed.tre`
- `trees/<target_id>_chronos_dated_modeldiscrete.tre`

## Practical Reading Of Results

- If fit and tempo both favor the same model: strong convergence.
- If fit favors one model but tempo favors another: report both explicitly and choose based on your biological objective (clock-model fit vs preservation of branching-tempo pattern).

## Common Issues

1. Many calibration pairs but few calibration nodes:
   - some taxa are missing in target tree, or multiple pairs collapse to same MRCA.
2. Congruify warning about non-ultrametric reference:
   - the script handles this internally by ultrametricizing the reference.
3. No usable calibrations:
   - check exact label matching between tree and calibrations.
