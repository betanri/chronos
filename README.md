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
3. Optionally subsets large trees before fitting (while preserving calibration signal and tempo extremes).
4. Fits chronos across clock models and lambda values.
5. Applies robust model selection with threshold sensitivity (`1` and `2`).
6. Computes branching-tempo metrics (overall and early-tempo) to compare chronograms against the input phylogram branching pattern.
7. Writes:
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

## Optional Subset Mode For Large Trees

Use this when full-tree fitting is too slow.

```r
USE_SUBSET <- TRUE
SUBSET_N <- 400L
SUBSET_EXTREME_FRAC <- 0.05
SUBSET_SEED <- 1L
```

Subset strategy:

1. Always keep all taxa used in calibration pairs (`taxonA`, `taxonB`).
2. Add root-to-tip extremes from the phylogram (shortest and longest paths).
3. Fill remaining tips using diversified topological spread (ladderized traversal).
4. If needed, fill final slots at random from remaining tips.

If calibration taxa alone exceed `SUBSET_N`, the script stops and asks you to increase `SUBSET_N`.

## How To Run

```r
setwd("path/to/project/Terapontoid_Trees")
source("Run_chronos_pipeline.R")
```

## Model-Fit Sensitivity

The script evaluates and reports both:

- `PLOG_CLOCK_SWITCH_THRESH = 1` (default strict)
- `PLOG_CLOCK_SWITCH_THRESH = 2` (stricter)


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

Output folder structure:

1. `main_files/`
   - compact deliverables for direct review:
     - threshold-selected trees
     - one tree per model
     - fit summary tables
     - branching-tempo metric table
     - interpretation text
     - subset tip list (when subset mode is enabled)
2. `all_files/`
   - full run contents:
     - `tables/`
     - `trees/`
     - `logs/`
     - `checkpoints/`

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
