# run_all.R: regenerates every number, table and figure. Run from the repository root:
#   Rscript run_all.R
# Scripts are added here as they are written (see README.md for the full layout).

scripts <- c("R/00_load_clean.R",   # data, checks, working-hour clock, daily aggregates
             "R/01_eda.R",          # data section: summary table and figures
             "R/10_type1.R")        # Type 1: estimates, parametric and nonparametric bootstrap

for (f in scripts) {
  cat("\n=====", f, "=====\n")
  source(f)
}
cat("\nDone. Results are in output/.\n")
