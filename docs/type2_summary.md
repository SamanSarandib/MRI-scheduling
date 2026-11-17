## Type 2 patients (`R/11_type2.R`)

All numbers below come straight from the data, with a non-parametric bootstrap for the 95% intervals; a smoothed
bootstrap gives practically the same intervals. Four distributions (exponential, gamma, Weibull, lognormal) were also
fitted and checked with a Kolmogorov–Smirnov (KS) test, histograms and QQ plots.
Results below are for 209 Type 2 scans over 21 working days.

### Scan durations

A Type 2 scan takes 40.7 minutes on average, and 95% of scans finish within about 63 minutes.

| Quantity | Estimate | 95% interval |
| --- | --- | --- |
| Mean duration | 40.7 min | 39.1 to 42.4 |
| 95% of scans finish within | 63.0 min | 58.6 to 67.6 |
| Share of scans over 45 min | 33.5% | 27% to 40% |
| Share of scans over 60 min | 6.2% | 3.3% to 9.6% |

The gamma distribution fits best (KS p = 0.97) and the lognormal fits too (p = 0.63); the Weibull, normal and
exponential are rejected. The gamma gives almost the same numbers as the data (95% quantile 62.3 against 63.0 min).

### Calls

Unlike Type 1, Type 2 calls are not random (not Poisson). They are very regular: about 10 a day, and every day
had between 8 and 11 calls, where random calls would give roughly 5 to 15.

| Quantity | Estimate | 95% interval | If calls were random |
| --- | --- | --- | --- |
| Calls per working day | 9.95 | 9.5 to 10.3 | |
| Variance / mean of daily calls | 0.08 | 0.04 to 0.12 | 1 |
| Mean time between calls | 54 working min | 51 to 57 | |
| CV of time between calls | 0.37 | 0.34 to 0.41 | 1 |

A Weibull fits the time between calls (KS p = 0.21); the exponential and lognormal are rejected (p = 0.001).
For the simulation, either the Weibull or the non-parametric gaps can be used.

### What it means for scheduling

Type 2 needs about 6.75 scan-hours a day (6.4 to 7.2) and no day exceeded 9 hours, so one machine can handle
the workload. The slot length is a trade-off:

| Slot length | Slots per 9-hour day | Scans that overrun | Days with more calls than slots |
| --- | --- | --- | --- |
| 45 min | 12 | 33.5% | 0% |
| 60 min | 9 | 6.2% | 71% (52% to 91%) |

The best slot length is probably between 45 and 60 minutes, which the simulation can test.

### Outputs (`output/`)

- `tab11_type2_durations.csv`, `tab11_type2_arrivals.csv`, `tab11_type2_gaps.csv`: results with 95% intervals
- `tab11_type2_gof.csv`, `tab11_type2_gof_gaps.csv`: KS fit tests
- `fig12`: bootstrap distributions; `fig13` / `fig15`: histograms with fitted models (durations / time between
  calls); `fig14` / `fig16`: QQ plots
- `type2.rds`: everything for the simulation and the Monte Carlo

### Using the results

- **Simulation (Part 2):** `readRDS("output/type2.rds")$inputs` has generators in minutes. Durations: `dur2_gamma`
  (or `dur2`, non-parametric). Time between calls: `gap2` (non-parametric) or `gap2_weibull`. Not an exponential
  as for Type 1. Slot lengths to test: between 45 and 60 minutes.
- **Monte Carlo:** after `source("R/11_type2.R")` you have `fits` / `gap_fits` (fitted parameters), `models`
  (with `r()` to simulate), `model_stats()` (true values under a model) and `dur_stats()` + `boot_table()`
  (the estimation method). Good "true" distributions: gamma and lognormal for durations, Weibull for the gaps.