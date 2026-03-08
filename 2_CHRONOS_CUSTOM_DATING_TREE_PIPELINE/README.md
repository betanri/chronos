# Chronos -- Custom Dating Tree Pipeline

## Benchmark Evidence

- [**Why Chronos and not treePL?**](../1_WHY_CHRONOS_AND_NOT_TREEPL/README.md)

This pipeline dates a phylogram with `ape::chronos` and gives you:

- one chronogram per clock model (`clock`, `correlated`, `relaxed`, `discrete`)
- an explicit comparison between:
  - model-fit preference
  - pulse-preservation similarity to the original phylogram
  - gap burden against the calibration information
  - rate plausibility relative to the source phylogram

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
6. Computes pulse-preservation metrics to compare chronograms against the input phylogram branching pattern.
7. Writes:
   - one tree per model
   - fit tables
   - pulse-preservation metric table
   - a plain-language interpretation text file.

## Main Script

- `2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/Run_chronos_pipeline.R`

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

Use this when full-tree model-grid fitting is too slow.

```r
USE_SUBSET <- TRUE
SUBSET_N <- 400L
SUBSET_EXTREME_FRAC <- 0.05
SUBSET_SEED <- 1L
SUBSET_TUNE_ON_SUBSET_ONLY <- TRUE
```

Subset strategy:

1. Always keep all taxa used in calibration pairs (`taxonA`, `taxonB`), so calibration MRCAs are retained.
2. Add root-to-tip extremes from the phylogram (shortest and longest paths).
3. Fill remaining tips using diversified topological spread (ladderized traversal).
4. If needed, fill final slots at random from remaining tips.

If calibration taxa alone exceed `SUBSET_N`, the script stops and asks you to increase `SUBSET_N`.

How this is used in the pipeline:

1. The subset tree is used only to tune model/lambda/K and threshold sensitivity.
2. The selected settings are then applied to the full tree for final dated outputs.
3. This avoids running the full model grid on the full tree while still producing a full-tree chronogram.

## How To Run

```r
setwd("path/to/project/2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE")
source("Run_chronos_pipeline.R")
```

## Model-Fit Sensitivity

The script evaluates and reports both:

- `PLOG_CLOCK_SWITCH_THRESH = 1` (default strict)
- `PLOG_CLOCK_SWITCH_THRESH = 2` (stricter)

## Pulse-Preservation Metric (Biology-Oriented Diagnostic)

Detailed walkthrough:

- [Pulse-Preservation Metric Guide](BRANCHING_TEMPO_METRIC_GUIDE.md)

This metric family asks whether a dated tree preserves clustering of branching bursts and quiet intervals rather than flattening them into an unrealistically even schedule.

How the updated metric family is interpreted:

1. Match internal nodes by clade identity when node-based comparisons are needed.
2. Compare the distribution and clustering of branching events through relative time.
3. Penalize loss of burst structure and reward preservation of diversification pulses.
4. Read these results together with gap burden and rate plausibility rather than as a standalone winner metric.

## Three Complementary Metrics

The updated comparison layer is meant to be read through three complementary metrics.

### 1. Pulse preservation

- asks whether a dated tree preserves branching bursts and quiet intervals from the source phylogram
- this is the main tree-shape diagnostic in the current workflow

### 2. Gap burden

- if calibrations are fossil minima or bounded intervals, this behaves as implied ghost-lineage / fossil-gap burden
- if calibrations are exact ages, this behaves as calibration slack
- lower values mean the dated tree sits closer to the calibration information

### 3. Rate plausibility

- converts phylogram branch lengths and chronogram branch durations into implied branchwise rates
- penalizes extreme rate dispersion, abrupt parent-child jumps, weak local autocorrelation, and too many outlier rates

## Outputs (Dedicated Run Folder)

Outputs are written to:

- `<OUT_BASE_DIR>/`

Output folder structure:

1. `main_files/`
   - compact deliverables for direct review:
     - phylogram used for dating (`*_phylogram_used.tree`)
     - one tree per model
     - model-fit + pulse summary table (`summary_*_model_fits.csv`)
     - pulse-preservation metric table
     - interpretation text
     - subset tip list (when subset mode is enabled)
2. `all_files/`
   - full run contents:
     - `tables/`
     - `trees/`
     - `logs/`
     - `checkpoints/`
   - includes threshold-selected trees (`clockthresh*`, and `fulltree_clockthresh*` when subset-tuning mode is enabled)

Notes:

- `OUT_BASE_DIR` defaults to `CHRONOS_OUT`.
- `OUT_PREFIX` is used for file names (summary, interpretation, logs), not as a subfolder name.

## Practical Reading Of Results

- If fit and pulse preservation favor the same model: strong convergence.
- If fit favors one model but pulse preservation favors another: report both explicitly.
- If gap burden or rate plausibility disagree with both fit and pulse preservation, treat that as biologically meaningful conflict rather than noise.
- Do not choose from fit alone when the dated trees imply very different branching rhythm or very different rate behavior.

## Common Issues

1. Many calibration pairs but few calibration nodes:
   - some taxa are missing in target tree, or multiple pairs collapse to same MRCA.
2. Congruify warning about non-ultrametric reference:
   - the script handles this internally by ultrametricizing the reference.
3. No usable calibrations:
   - check exact label matching between tree and calibrations.
