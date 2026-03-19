# Benchmark suite input schema

These files are the intended drop-in summaries for the rebuilt `A-E` benchmark page.

## `benchmark_suite_method_summary.csv`

One row per benchmark-method entry used for the main benchmark summary plots.

Required columns:

- `benchmark`
- `benchmark_label`
- `completed_reps`
- `method`
- `mae`
- `rank`
- `group`

Suggested `method` values:

- `treePL`
- `RelTime`
- `chronos overall`
- `chronos-clock`
- `chronos-discrete`
- `chronos-correlated`
- `chronos-relaxed`

Suggested `group` values:

- `core`
- `chronos_model`

`RelTime` should be treated as a core method alongside `chronos overall` and `treePL`, not as a supplementary series.

## `benchmark_suite_recovery_summary.csv`

One row per benchmark and true `chronos` model.

Required columns:

- `benchmark`
- `benchmark_label`
- `completed_reps`
- `true_model`
- `recovery_rate`
- `n_total`
- `n_recovered`

Suggested `true_model` values:

- `clock`
- `discrete`
- `correlated`
- `relaxed`

## `benchmark_suite_rank_summary.csv`

One row per benchmark and method family for compact rank visualizations.

Required columns:

- `benchmark`
- `benchmark_label`
- `method_family`
- `mean_rank`

Suggested `method_family` values:

- `chronos`
- `RelTime`
- `treePL`

## `tenv_summary.csv`

Compact `treePL` environment-diagnostic summary.

Required columns:

- `run`
- `label`
- `bundles_complete`
- `mean_treePL_mae`
- `median_treePL_mae`
- `mean_delta_vs_baseline`
- `mean_runtime_sec`

Suggested `run` values:

- `tenv_local`
- `tenvp_local`
- `tenv_oscer`

## `representative_tree_panel_manifest.csv`

One row per small tree panel used in `Fig 4`.

Required columns:

- `benchmark`
- `benchmark_label`
- `row_order`
- `panel_order`
- `panel_label`
- `method`
- `condition`
- `image_path`

Suggested `method` values:

- `reference`
- `chronos`
- `RelTime`
- `treePL`

Intended use:

- one row per benchmark (`A-E`)
- the same four panel columns within each row:
  - `reference`
  - `chronos`
  - `RelTime`
  - `treePL`
- same output figure can later be extended to `P1/P2` if needed

## `pcr_summary_table.csv`

Compact `PCR` summary table shown on the benchmark page.

Required columns:

- `dataset`
- `metric_family`
- `best_method`
- `summary_value`
- `interpretation`

Suggested `metric_family` values:

- `pulse preservation`
- `gap burden`
- `rate irregularity`
- `overall post-fit rank`

## Legacy exact-root input

The historical plotting script [../make_figures_and_summary.R](../make_figures_and_summary.R) still reads the old exact-root benchmark CSV directly.

Minimum expected legacy columns:

- `condition`
- `clock_model`
- `true_chronos_model`
- `mu`
- `heterotachy`
- `replicate`
- `treePL_MAE`
- `chronos_MAE`

Optional legacy column:

- `reltime_MAE`

If `reltime_MAE` is present, the legacy summaries and figures should include `RelTime` on equal footing with `chronos` and `treePL` rather than silently dropping it.
