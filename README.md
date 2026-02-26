# chronos pipeline

Simple, reproducible pipeline for divergence-time dating with `ape::chronos`.

It supports two calibration workflows:

1. **Reference-tree workflow (congruification)**  
   Build internal calibrations from a reference time tree.

2. **Manual workflow**  
   Provide internal calibrations directly in a CSV file.

The script is designed to be easy to edit and run on any dataset.

---

## What this does

Given a target phylogram:

- reads the tree (single-tree, multi-tree, or named-Newick)
- builds chronos calibrations (congruify or manual)
- fits chronos across model/lambda combinations
- applies a robust molecular clock model selector
- writes dated tree + summary + checkpoint files

---

## Main script

- `Run_chronos_pipeline.R`



---

## Requirements

R packages:

- `ape`
- `geiger`
- `phytools`

Install once:

```r
install.packages(c("ape", "geiger", "phytools"))
```

---

## Input options

### Option A: Congruification (reference time tree)

Set in script:

- `USE_CONGRUIF <- TRUE`
- `TARGET_TREE_FILE <- "path/to/target_phylogram.tre"`
- `REFERENCE_TIME_TREE <- "path/to/reference_time_tree.tre"`

The script will:
- make reference ultrametric if needed
- run congruification
- convert resulting constraints into chronos calibrations

---

### Option B: Manual calibrations (CSV)

Set in script:

- `USE_CONGRUIF <- FALSE`
- `TARGET_TREE_FILE <- "path/to/target_phylogram.tre"`
- `MANUAL_CAL_CSV <- "path/to/manual_calibrations.csv"`

CSV required columns:

- `taxonA`
- `taxonB`
- `age_min`
- `age_max`

Each row defines one calibration interval on MRCA(taxonA, taxonB).

---

## Run

In R:

```r
source("Run_chronos_pipeline.R")
```

---

## Empirical model-sensitivity protocol (recommended)

For empirical trees, run a separate model-sensitivity protocol and report it.
Run the same tree with:

- `PLOG_CLOCK_SWITCH_THRESH = 1` (current strict)
- `PLOG_CLOCK_SWITCH_THRESH = 0` (less strict)
- `PLOG_CLOCK_SWITCH_THRESH = 2` (more strict)

Keep other settings fixed (`PLOG_NONCLOCK_SWITCH_THRESH`, `PLOG_TIE_EPS`, `LAMBDA_GRID`, `CHRONOS_MODELS`, `K_FIT_GRID`) and compare selected model/lambda and dated-tree stability.

Minimal pattern:

```r
# strict
PLOG_CLOCK_SWITCH_THRESH <- 1
source("Run_chronos_pipeline.R")

# less strict
PLOG_CLOCK_SWITCH_THRESH <- 0
source("Run_chronos_pipeline.R")

# more strict
PLOG_CLOCK_SWITCH_THRESH <- 2
source("Run_chronos_pipeline.R")
```

---

## Outputs

Configured by `OUT_DIR` and `OUT_PREFIX` in the script.

Typical outputs:

- `summary_<prefix>.csv` (selected model/lambda and run summary)
- `results_<prefix>.rds` (full fit/calibration objects)
- `checkpoint_<prefix>.rds` (restart-friendly checkpoint)
- `<target_id>_chronos_dated.tre` (dated output tree)
- `run_<prefix>.log` (run log)

---

## Key settings you can change

- `LAMBDA_GRID`
- `CHRONOS_MODELS`
- `K_FIT_GRID` (for discrete model)
- robust selector thresholds:
  - `PLOG_CLOCK_SWITCH_THRESH`
  - `PLOG_NONCLOCK_SWITCH_THRESH`
  - `PLOG_TIE_EPS`

---

## Notes

- The current default selector corresponds to the robust configuration used in prior benchmark tuning.
- Selection is based on model-fit quantities from the target phylogram (not true ages), so it is suitable for empirical datasets.
