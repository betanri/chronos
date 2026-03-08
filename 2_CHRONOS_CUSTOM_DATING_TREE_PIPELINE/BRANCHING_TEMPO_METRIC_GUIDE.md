# Post-Fit Evaluation Metrics Guide (Terap Example Output)

This page shows how to read the Terap example after `clock fitting` and `lambda tuning` are finished. The goal here is to compare the dated trees that come out of the pipeline, not to rerun model fitting.

## What this guide adds

The pipeline keeps three layers separate:

1. `clock fitting`: which chronos clock process is favored by the fit statistics
2. `lambda tuning`: which smoothing strength is favored within the chosen model search
3. `post-fit evaluation`: how the resulting dated trees behave biologically after fitting is done

That third layer is the focus here. It uses three complementary metric families:

- `pulse preservation`: does the chronogram keep the branching bursts and quiet intervals seen in the source phylogram?
- `gap burden`: how much slack does the dated tree imply relative to the calibration information?
- `rate plausibility`: do the implied branchwise rate changes stay within a biologically defensible range?

## Example files used

- Input phylogram: `2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES/OUTPUT_DEMO/TERAP_ML_MAIN_phylogram_used.tree`
- Bundled chronos trees:
  - `.../TERAP_ML_MAIN_chronos_dated_modelclock.tre`
  - `.../TERAP_ML_MAIN_chronos_dated_modelcorrelated.tre`
  - `.../TERAP_ML_MAIN_chronos_dated_modelrelaxed.tre`
  - `.../TERAP_ML_MAIN_chronos_dated_modeldiscrete.tre`
- Fit summary table:
  - `2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES/OUTPUT_DEMO/summary_terap_empirical_model_fits.csv`
- Post-fit summary table:
  - `2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES/OUTPUT_DEMO/summary_terap_empirical_postfit_metrics.csv`

The post-fit summary table includes the same four chronos trees plus the `treePL` comparator used in the Terap comparison.

## Terap example: fit layer vs post-fit layer

The fit layer and the post-fit layer are close here, but they are not identical.

- `clock` has the best `PHIIC` in the bundled fit summary
- `discrete` has the best penalized log-likelihood and the best pulse-focused composite from the older output
- in the post-fit comparison, `chronos_discrete` is the best overall tree across the four ranking components used here
- `chronos_clock` is a near-tie second and is the best tree for `rate plausibility`

So this example is not a case where one model dominates every decision layer. It is a case where `clock` and `discrete` stay at the top, but for partly different reasons.

## Ranked post-fit results (lower is better)

In this Terap example, `gap burden` behaves as `point-calibration slack`, not as fossil-minimum ghost-lineage burden, because the comparison uses point calibrations.

| candidate | pulse preservation (default) | pulse preservation (burst) | gap burden | rate plausibility | overall mean rank |
| --- | ---: | ---: | ---: | ---: | ---: |
| `chronos_discrete` | `0.1603` | `0.1462` | `0.0733` | `2.5982` | `1.25` |
| `chronos_clock` | `0.1605` | `0.1464` | `0.0736` | `2.5861` | `1.75` |
| `chronos_correlated` | `0.1722` | `0.1478` | `0.1058` | `3.4566` | `3.25` |
| `treePL` | `0.1729` | `0.1604` | `0.1615` | `3.4631` | `4.25` |
| `chronos_relaxed` | `0.1949` | `0.1684` | `0.1030` | `4.5535` | `4.50` |

Metric-family leaders in this example:

- `pulse preservation (default)`: `chronos_discrete`
- `pulse preservation (burst)`: `chronos_discrete`
- `gap burden`: `chronos_discrete`
- `rate plausibility`: `chronos_clock`

## Figure A: Pulse-layer tree-shape comparison among bundled chronos trees

![Pulse preservation tree panel](figures/branching_tempo_tree_panel_clean_v3.png)

This figure is useful because it shows the pulse layer directly on the bundled `chronos` trees only; `treePL` is not shown in this panel. It helps explain why `discrete` and `clock` sit at the top of the pulse-preservation ranking. But this panel is only for the pulse issue. It does not show the `gap burden` or `rate plausibility` parts of the broader post-fit comparison.

## Figure B: Post-fit comparison across metric families

![Post-fit evaluation metric families](figures/postfit_metric_family_values.png)

How to read:

- each panel is one post-fit metric family
- lower values are better in every panel
- the last panel summarizes the mean rank across the four ranking components

Interpretation for this example:

- `chronos_discrete` is the overall post-fit winner because it leads both pulse summaries and gap burden while staying near-best on rate plausibility
- `chronos_clock` is essentially tied at the top on pulse preservation, nearly tied on gap burden, and is the best tree on rate plausibility
- `chronos_correlated` sits in the middle
- `treePL` beats `chronos_relaxed` on both pulse summaries and on rate plausibility, but it has the highest gap burden in this comparison
- `chronos_relaxed` is worst overall because it combines the weakest pulse preservation with the highest rate irregularity

## Practical decision rule

For the Terap example:

1. If you want one overall post-fit winner, choose `chronos_discrete`.
2. If you want the best implied rate behavior, choose `chronos_clock`.
3. If fit-based selection and post-fit evaluation point to different trees, report both explicitly rather than collapsing them into one claim.
4. In this example, `treePL` is not the leading solution under the post-fit layer.

## Reproducibility

- Figure script:
  - `2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/scripts/make_postfit_metric_guide_figures.R`
- Post-fit summary table:
  - `2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES/OUTPUT_DEMO/summary_terap_empirical_postfit_metrics.csv`
- Fit summary table:
  - `2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE/EXAMPLE_FILES/OUTPUT_DEMO/summary_terap_empirical_model_fits.csv`
