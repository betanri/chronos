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
5. Writes dated tree + summary + checkpoint outputs.

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

These correspond to the current robust selector used in your simulation benchmarks (`c1_n2_t2`).

---

## Outputs

Configured by:

- `OUT_DIR`
- `OUT_PREFIX`

Generated files:

- `summary_<OUT_PREFIX>.csv` : selected model/lambda + calibration counts
- `results_<OUT_PREFIX>.rds` : full fit object + calibration object
- `checkpoint_<OUT_PREFIX>.rds` : lightweight checkpoint
- `<target_id>_chronos_dated.tre` : dated output tree
- `run_<OUT_PREFIX>.log` : run log

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

