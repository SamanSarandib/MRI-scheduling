# MRI scheduling

Should the hospital let both patient types use either of its two MRI scanners, and how long
should the appointment slots be? Part I estimates the input models from one month of call
records (August 2025) with bootstrap uncertainty; Part II simulates the old and the new
booking policy day by day.

## How to run

Open `MRI-scheduling.Rproj` (or set the working directory to this folder) and run

```r
source("run_all.R")
```

or `Rscript run_all.R` from a terminal. Every table and figure is written to `output/`,
which is regenerated on each run and not tracked by git. Only base R is used.

## Layout

| File | What it does | Status |
|---|---|---|
| `run_all.R` | sources the scripts below in order | done |
| `data/raw/ScanRecords.csv` | the call records, never edited by hand | done |
| `R/00_load_clean.R` | dates, minutes, working-hour clock, interarrivals, daily aggregates, checks | done |
| `R/01_eda.R` | summary table and figures for the data section | done |
| `R/10_type1.R` | normal/Poisson estimates, parametric and nonparametric bootstrap | done |
| `R/11_type2.R` | plug-in estimates, iid, smoothed and parametric bootstraps, ML fits of four models, KS fit tests | done |
| `R/20_boot_functions.R` | percentile, percentile-t, parametric, smoothed interval functions | done |
| `R/30_mc_study.R` | Monte Carlo coverage study of the intervals in 10 and 11, parallelised | done |
| `R/40_des_engine.R` | `generate_calls()`, `book_call()`, `run_day()`, `simulate()` | to do |
| `R/41_des_experiments.R` | slot grid, policy comparison, input-uncertainty loop | to do |
| `R/50_report_outputs.R` | tables and figures for the report | to do |

## Conventions

- Run everything from the repository root; relative paths only.
- `set.seed()` at the top of every script that draws random numbers.
- Bootstrap loops vectorised where possible (one column per bootstrap sample).
- Durations are in hours in the data and in minutes in tables, the working-hour clock
  ignores time outside 8:00 to 17:00, so a working day is 9 hours.
- Results are saved as `.rds`, figures as PDF.
