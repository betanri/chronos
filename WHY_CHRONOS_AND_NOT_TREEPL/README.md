# Benchmark suite: chronos, treePL, and RelTime

This page is now the **main benchmark-suite summary** for `chronos`, `treePL`, and `RelTime`.

The folder name is historical. The benchmark itself is no longer framed as only “chronos versus treePL”; `RelTime` is part of the main comparison layer throughout this page.

The old exact-root 720-run benchmark was useful, but it is no longer enough on its own. The active benchmark story is now broader:

- `A-E` main benchmarks across different calibration regimes
- `RelTime` included as a full equal, not as an afterthought
- future `P1/P2` extensions for tree-size scaling
- compact `treePL` environment diagnostics (`Tenv`)
- one linked `PCR` section for post-fit interpretation

Final numbers will be plugged in as the current runs finish. The point of this page is to put the **new structure** in place now, so the repo already matches the benchmark we are actually running.

## Benchmark suite

| Benchmark | Purpose | Calibration design | Status |
|---|---|---|---|
| `A` | idealized exact-root benchmark | exact root age only | active |
| `B` | sparse wide-bracket benchmark | root bracket + 5 internal brackets | active |
| `C` | sparse tight-bracket benchmark | root bracket + 5 tight internal brackets | active |
| `D` | sparse minimum-only benchmark | root bracket + 5 internal minimums | active |
| `E` | richer wide-bracket benchmark | root bracket + 10 internal brackets | active |
| `P1` | small-tree extension | same logic on smaller trees | planned |
| `P2` | large-tree extension | same logic on larger trees | planned |

## Reporting layers

This suite keeps four layers separate:

1. `dating accuracy`
   - MAE by method and by benchmark
2. `chronos model behavior`
   - `clock`, `discrete`, `correlated`, and `relaxed` reported separately
3. `treePL environment sensitivity`
   - compact `Tenv` comparison on identical saved bundles
4. `post-fit chronogram plausibility`
   - summarized through the linked `PCR` framework

That separation matters. A method can date well while still recovering clock models poorly. A dated tree can also look plausible or implausible after fitting, which is a different question again. `RelTime` is part of the main dating comparison, not a side method.

## Main figure set

The repo is being reorganized around a compact main figure set:

- `Fig 1`: overall cross-benchmark rank summary across `A-E`
- `Fig 2`: per-benchmark MAE panels with `chronos` split by model, plus `RelTime` and `treePL`
- `Fig 3`: clock-model recovery across benchmarks
- `Fig 4`: representative tree-shape comparison panel with **one row per benchmark (`A-E`)**
- `PCR` figure: compact empirical post-fit summary
- `PCR` table: compact metric/dataset interpretation table

The old single-benchmark figure stack is not the main story anymore, so it is no longer surfaced here.

## What the benchmark is testing

The main benchmark asks:

- when root age is exact, sparse, bracketed, or minimum-only, which method dates best?
- how much does the answer depend on the true clock regime?
- how often does `chronos` recover the correct clock model?
- does `treePL` behave consistently across environments?

The benchmark is therefore no longer just “chronos versus treePL under one favorable setup.” It is now a structured test of **calibration regime**, **clock regime**, **heterotachy**, **extinction**, and later **tree size**.

## Current benchmark outputs to report

The final page will report at least these tables:

- benchmark-level MAE summary across `A-E`
- benchmark-level MAE summary with `RelTime` reported on equal footing
- model-specific `chronos` MAE by benchmark
- clock-model recovery summary by benchmark
- method ranking summary across benchmarks
- compact `PCR` metric table

The expected input layout for those tables and figures is documented in [data/README.md](data/README.md).

## Representative tree panel

`Fig 4` will be a compact multi-row panel:

- one row each for `A`, `B`, `C`, `D`, and `E`
- the same four columns in every row:
  - `reference`
  - `chronos`
  - `RelTime`
  - `treePL`
- intended to show the benchmark-specific shape differences without taking over the whole page

This panel is meant to stay interpretable at a glance. It is not supposed to be a full gallery of every condition.

## `treePL` environment diagnostic

<details>
<summary><strong>Tenv: compact treePL environment diagnostic</strong></summary>

`Tenv` is a small side benchmark used to test whether `treePL` behaves the same way across environments on the **same exact saved input bundles**.

It is not a replacement for the main `A-E` benchmark. It is a targeted diagnostic layer that asks a different question:

- do local and OSCER reruns of the same `treePL` input reproduce one another?
- when they disagree, is the discrepancy mild or catastrophic?

The reporting here should stay compact:

- local vs OSCER summary
- baseline delta summary
- one short interpretation paragraph

</details>

## `PCR` post-fit interpretation

The benchmark page also needs one compact `PCR` section, because the benchmark and `PCR` answer different but connected questions.

The benchmark asks:

- which method dates better under controlled simulation?

`PCR` asks:

- once trees have been fitted, which chronograms behave more plausibly?

This page will therefore include:

- one `PCR` summary figure
- one `PCR` summary table

The point is not to duplicate the whole `PCR` repo here. The point is to make the link explicit, visible, and interpretable from this benchmark page.

## Future extensions

`P1` and `P2` are reserved for the later tree-size extension. They should be integrated into the same page structure once those runs are ready, not split into separate benchmark pages.

## Figure/data scaffold

The new suite-level figure scaffold lives here:

- [make_benchmark_suite_figures.R](make_benchmark_suite_figures.R)
- [data/README.md](data/README.md)

The legacy exact-root plotting script also remains available for the historical 720-run dataset:

- [make_figures_and_summary.R](make_figures_and_summary.R)

That legacy script should now be read as a small benchmark utility for `chronos`, `treePL`, and `RelTime` whenever a `reltime_MAE` column is present in the input CSV. It is no longer framed as a two-method plotting path.

Those files define the expected summary inputs so the final values can be dropped in without another repo redesign.
