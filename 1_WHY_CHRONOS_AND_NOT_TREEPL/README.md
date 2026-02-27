# Why Chronos and not treePL?

`treePL` is very widely used in divergence-time studies, so the default expectation is often to use it first.
Outside of computationally demanding Bayesian approaches (often not feasible for very large trees), this benchmark tests whether that default is justified under the simulation conditions I care about.

I compared **treePL** and **chronos** under the same topology, calibration, and replicate grid.

- Total tests: **720** (4 true clock regimes x 3 extinction levels x 2 heterotachy levels x 30 replicates)
- Same topology and root calibration for both methods
- treePL tuning: smooth in `{0.1, 1, 10}`
- chronos tuning: model in `{clock, correlated, relaxed, discrete}`, lambda in `{0.1, 1, 10}`, robust selector

## Headline result

- Mean MAE: treePL **1.8113** vs chronos **0.4966**
- Head-to-head (both finite, n=635): chronos better in **558/635 (87.9%)**
- treePL failed to return finite MAE in **85/720** tests

## By true clock model (mean MAE; lower is better)

- strict (n both=180): treePL **1.0888** | chronos **0.0041**
- autocorrelated (n both=136): treePL **4.1675** | chronos **0.9149**
- independent (n both=160): treePL **1.1317** | chronos **0.4363**
- discrete (n both=159): treePL **1.2980** | chronos **0.6309**

## By heterotachy (mean MAE)

- H=0.05 (n both=341): treePL **1.2876** | chronos **0.0885**
- H=0.25 (n both=294): treePL **2.4188** | chronos **0.9046**

## By extinction (mu; mean MAE)

- mu=0.0 (n both=207): treePL **0.7540** | chronos **0.2256**
- mu=0.5 (n both=211): treePL **1.1017** | chronos **0.3711**
- mu=0.8 (n both=217): treePL **3.5100** | chronos **0.8930**

## Chronos model recovery (true simulated model vs selected model)

- Overall exact recovery: **313/720 (43.5%)**
- clock: **179/180 (99.4%)**
- correlated: **134/180 (74.4%)**
- discrete: **0/180 (0.0%)**
- relaxed: **0/180 (0.0%)**

This means chronos is much better on age accuracy in this benchmark, but its selector tends to collapse relaxed/discrete scenarios into clock/correlated solutions under this robust setup.

## Figures

![Overall MAE boxplot](figures/fig1_overall_mae_boxplot.png)

![MAE by true clock](figures/fig2_mae_by_true_clock.png)

![MAE heatmaps by heterotachy and extinction (same scale)](figures/fig3_mae_heatmap_mu_heterotachy.png)

![Heterotachy x extinction interaction plot](figures/fig4_interaction_heterotachy_extinction.png)

## Practical interpretation

- For these conditions, chronos is the better default for dating accuracy and stability.
- treePL can still be informative, but it degrades strongly in harder regimes.
- For empirical analyses, keep model-sensitivity reporting explicit even when chronos is used as the primary method.

## Data files

- `by_clock_summary.csv`
- `by_mu_summary.csv`
- `by_heterotachy_summary.csv`
- `chronos_recovery_summary.csv`
- `chronos_recovery_confusion_matrix.csv`
