# Chronos Empirical Pipeline (Terapontoids)

This folder contains an empirical dating workflow using `ape::chronos` with:

1. **Congruification-derived internal calibrations** (from a reference time tree), or  
2. **Manual internal calibrations** (CSV file).

Main script:

- `Run_chronos_pipeline.R`

---

## What This Script Does

1. Loads a **target phylogram** (single tree, multi-tree, or named-Newick file).
2. Builds chronos calibrations from either:
   - `USE_CONGRUIF <- TRUE`: derive calibrations from a reference time tree via `geiger::congruify.phylo`, or
   - `USE_CONGRUIF <- FALSE`: read calibrations from CSV.
3. Fits chronos across model/lambda options.
4. Selects final model using robust selector logic (`c1_n2_t2` defaults):
   - clock -> non-clock switch threshold
   - correlated vs relaxed switch threshold
   - tie handling toward correlated
5. Writes dated trees, model-fit tables, branching-tempo diagnostics, interpretation text, and checkpoint outputs.

---

## Required R Packages

Install once:

```r
install.packages(c("ape", "geiger", "phytools"))
```

Load in session:

```r
library(ape)
library(geiger)
library(phytools)
```

---

## Input Modes

### Mode A: Congruification (reference time tree)

Set in `Run_chronos_pipeline.R`:

```r
USE_CONGRUIF <- TRUE
TARGET_TREE_FILE <- "..."
REFERENCE_TIME_TREE <- "..."
```

Notes:
- If reference is non-ultrametric, script uses `force.ultrametric(..., method="extend")` before congruification.
- If `ROOT_AGE` is `NA`, root age is set from reference crown age.

### Mode B: Manual calibrations

Set:

```r
USE_CONGRUIF <- FALSE
TARGET_TREE_FILE <- "..."
MANUAL_CAL_CSV <- "manual_calibrations.csv"
```

Manual CSV required columns:

- `taxonA`
- `taxonB`
- `age_min`
- `age_max`

Extra columns are allowed and ignored.

---

## Target Tree Selection

If target file contains many trees:

- set `TARGET_TREE_NAME` (preferred), or
- set `TARGET_TREE_INDEX`.

If file is a single tree, these are ignored.

---

## How To Run

From R:

```r
setwd("/Users/ricardobetancur/Desktop/Proxy_Misplaced/chronosPL/Terapontoid_Trees")
source("Run_chronos_pipeline.R")
```

---

## Key Tuning Parameters (current defaults)

- `LAMBDA_GRID <- c(0.1, 1, 10)`
- `CHRONOS_MODELS <- c("clock","correlated","relaxed","discrete")`
- `K_FIT_GRID <- c(2L, 3L, 5L, 10L)`
- `N_RETRIES <- 2L`
- `PLOG_CLOCK_SWITCH_THRESH <- 1`
- `PLOG_NONCLOCK_SWITCH_THRESH <- 2`
- `PLOG_TIE_EPS <- 2`
- `CLOCK_SWITCH_SENSITIVITY <- c(1, 2)`

These correspond to the current robust selector used in your simulation benchmarks (`c1_n2_t2`).

---

## Outputs

Configured by:

- `OUT_BASE_DIR`
- `OUT_PREFIX`

Generated files (in dedicated run folder):

- `<OUT_BASE_DIR>/<OUT_PREFIX>/tables/...`
- `<OUT_BASE_DIR>/<OUT_PREFIX>/trees/...`
- `<OUT_BASE_DIR>/<OUT_PREFIX>/logs/...`
- `<OUT_BASE_DIR>/<OUT_PREFIX>/checkpoints/...`

Key files:

- `tables/summary_<OUT_PREFIX>.csv` : default selected row (threshold 1)
- `tables/summary_<OUT_PREFIX>_sensitivity.csv` : threshold 1 vs 2 side-by-side
- `tables/summary_<OUT_PREFIX>_model_fits.csv` : per-model fit summary plus branching-tempo metrics
- `tables/interpretation_<OUT_PREFIX>.txt` : explicit fit-vs-tempo conclusion
- `tables/results_<OUT_PREFIX>.rds` : full fit object + calibration object
- `tables/<target_id>_calibrations_used.csv` : calibration pairs used
- `trees/<target_id>_chronos_dated_clockthresh1.tre` : dated tree at threshold 1
- `trees/<target_id>_chronos_dated_clockthresh2.tre` : dated tree at threshold 2
- `trees/<target_id>_chronos_dated_modelclock.tre` : one tree for clock model
- `trees/<target_id>_chronos_dated_modelcorrelated.tre` : one tree for correlated model
- `trees/<target_id>_chronos_dated_modelrelaxed.tre` : one tree for relaxed model
- `trees/<target_id>_chronos_dated_modeldiscrete.tre` : one tree for discrete model
- `logs/run_<OUT_PREFIX>.log` : run log
- `checkpoints/checkpoint_<OUT_PREFIX>.rds` : lightweight checkpoint

If `CLEAN_PREVIOUS_PREFIX_OUTPUTS <- TRUE`, prior outputs for the same prefix are moved to:
- `_archive/<timestamp>/`

## Empirical Model Sensitivity (recommended)

For empirical analyses, report both threshold settings:

- `PLOG_CLOCK_SWITCH_THRESH = 1` (default strict)
- `PLOG_CLOCK_SWITCH_THRESH = 2` (stricter)

Do not use threshold `0` in the default empirical workflow.

## Branching-Tempo Diagnostic (new)

In addition to model fit (PHIIC/ploglik), the pipeline computes a branching-tempo diagnostic against the input phylogram for each model tree.

Node matching:
- Internal nodes are matched by clade identity (descendant tip set).

Tempo variables:
- Node heights are normalized to `[0,1]` using max root-to-tip depth.

Reported metrics in `summary_<OUT_PREFIX>_model_fits.csv`:
- `tempo_mae_all` : mean absolute error across all matched internal-node normalized heights.
- `tempo_mae_early_q75` : MAE restricted to early/deep nodes (top quartile of phylogram node heights).
- `tempo_median_early_q75` : median absolute error for the same early/deep node subset.

Interpretation logic:
- Fit-favored model comes from selector thresholds (1 and 2).
- Tempo-favored model is the model with lowest `tempo_mae_all` and/or `tempo_mae_early_q75`.
- `interpretation_<OUT_PREFIX>.txt` reports both and states whether they agree.

---

## Common Issues

1. **`pairs=...` is high but `Final calibration nodes on target tree` is low**
   - Many calibration tip pairs are not present in the target tree, or collapse onto same MRCA.

2. **Congruify warning about non-ultrametric reference**
   - Script handles this automatically by ultrametricizing reference before congruification.

3. **No usable calibration nodes**
   - Check taxon labels in target and calibration source.
   - Ensure names match exactly.

---

## Reproducibility Notes

- Script currently contains a hardcoded `setwd(...)`.  
  If sharing with others, replace with project-relative paths or remove `setwd`.
- Keep tree/csv file names stable for reproducible outputs.
