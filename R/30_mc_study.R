# 30_mc_study.R
# Can the hospital trust our intervals? Invent a truth, draw fake data from it, let every method
# estimate without seeing the truth, and count how often its 95% interval catches the truth:
# simulate, apply, evaluate, summarise, as in BootstrapCode.R. The methods are the interval
# rules of 10_type1.R and 11_type2.R, built with the same boot_table() calls:
#   Type 2 durations, four invented truths with the Type 2 mean and sd, n = 50 to 500:
#     M1 EDF, percentile          what 11 uses for the sd; its plain interval for q90 and q95
#     M2 EDF, percentile-t        what 11 uses for the mean and the shares (the plan's rule)
#     M3 normal fit               the Type 1 recipe of 10 applied to Type 2 (the naive analyst)
#     M4 gamma fit                as in 11, with Efron's interval for the shares
#     M5 smoothed bootstrap       as in 11; what its main table reports for q90 and q95
#     M6 EDF, Efron percentile    11's earlier rule for the shares (slightly too narrow)
#     M7 Weibull fit              as M4
#     M8 lognormal fit            as M4
#   M4, M7 and M8 under all four truths show what a model-based interval is worth when the model
#   is wrong (tab35).
#   Type 1: every quantity of 10 when its own model is the truth, at its sample sizes.
#   Type 2 gaps between calls: the mean gap, its CV and the calls per day from the gaps, with
#     the non-parametric and the model-based bootstraps of 11, under three gap shapes (tab34).
#   Type 2 days: the mean daily count (resampling days as 11, or assuming Poisson as 10), the
#     dispersion, the scan-hours per day, p_work9 and the shares of days of 11, when calls
#     follow a gamma renewal process or a Poisson process and durations a gamma or lognormal (tab31).
# Everything runs on four workers; the study takes about 10 minutes.
# Output: output/tab30-35_mc_*.csv, output/mc_study.rds, output/fig30_mc_coverage.pdf.

source("R/20_boot_functions.R")
library(parallel)
dir.create("output", showWarnings = FALSE)
set.seed(2026)                           # also keeps this session on the default generator after the cluster
n_sim <- 1000                            # Monte Carlo samples per cell (200 while developing)
B <- 499                                 # bootstrap draws per interval
alpha <- 0.05
n_grid <- c(50, 100, 209, 500)           # sample sizes; 209 is ours
probs <- c(0.90, 0.95); slots <- c(45, 60)                       # Type 2 quantiles and slots, as in 11
q_names <- paste0("q", 100 * probs); p_names <- paste0("p_over_", slots)
targets <- c("mean", "sd", q_names, p_names)
methods <- c(M1 = "EDF, percentile", M2 = "EDF, percentile-t", M3 = "normal fit", M4 = "gamma fit",
             M5 = "smoothed bootstrap", M6 = "EDF, Efron percentile", M7 = "Weibull fit", M8 = "lognormal fit")
rows <- function(L) do.call(rbind, L)

# ---- Type 2 durations -----------------------------------------------------------------------
# Step 1, the truths: four shapes with the Type 2 mean and sd (40.7 and 12.1 minutes). The
# abs() only guards the mixture against a negative draw (probability below 1e-8).
m <- 40.7; s <- 12.1
sl <- sqrt(log(1 + (s / m)^2))                                                                    # lognormal sdlog
kw <- uniroot(function(k) gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1 - (s / m)^2, c(1, 20), tol = 1e-10)$root  # Weibull shape
truths <- list(
  gamma     = list(r = function(n) rgamma(n, (m / s)^2, m / s^2), p = function(x) pgamma(x, (m / s)^2, m / s^2), sd = s),
  lognormal = list(r = function(n) rlnorm(n, log(m) - sl^2 / 2, sl), p = function(x) plnorm(x, log(m) - sl^2 / 2, sl), sd = s),
  weibull   = list(r = function(n) rweibull(n, kw, m / gamma(1 + 1 / kw)), p = function(x) pweibull(x, kw, m / gamma(1 + 1 / kw)),
                   sd = s),
  mixture   = list(r = function(n) abs(ifelse(runif(n) < 0.7, rnorm(n, 34.1, 6), rnorm(n, 56.1, 8))),
                   p = function(x) 0.7 * pnorm(x, 34.1, 6) + 0.3 * pnorm(x, 56.1, 8),
                   sd = sqrt(0.7 * (6^2 + 34.1^2) + 0.3 * (8^2 + 56.1^2) - m^2)))
true_q <- function(tr, p) uniroot(function(x) tr$p(x) - p, c(1, 200), tol = 1e-10)$root
truth_values <- sapply(truths, function(tr) c(mean = m, sd = tr$sd, setNames(sapply(probs, function(p) true_q(tr, p)), q_names),
                                              setNames(1 - tr$p(slots), p_names)))

# Step 2, the methods. Statistics for every column of X (one sample per column), with the
# standard errors ("se.") that make boot_table() choose percentile-t: the EDF as dur_stats()
# in 11, the normal fit as type1_quantities() in 10, and the model fits as model_stats() and
# tab_mod in 11, with the models of 20_boot_functions.R. col_quantiles() gives the same order
# statistics as col_quantile() for both probabilities, from one sort of the whole matrix.
col_quantiles <- function(X, p) { X <- as.matrix(X); matrix(X[order(col(X), X)], nrow(X))[ceiling(nrow(X) * p), , drop = FALSE] }
edf_stats <- function(X) {
  X <- as.matrix(X); P <- rows(lapply(slots, function(a) colMeans(X > a))); sdv <- col_sd(X)
  out <- rbind(colMeans(X), sdv, col_quantiles(X, probs), P, sdv / sqrt(nrow(X)), se_prop(P, nrow(X)))
  rownames(out) <- c(targets, "se.mean", paste0("se.", p_names)); out
}
normal_stats <- function(X) {
  X <- as.matrix(X); mu <- colMeans(X); sdv <- col_sd(X); Z <- rows(lapply(slots, function(a) (a - mu) / sdv))
  out <- rbind(mu, sdv, rows(lapply(probs, function(p) mu + qnorm(p) * sdv)), pnorm(Z, lower.tail = FALSE),
               sdv / sqrt(nrow(X)), dnorm(Z) * sqrt((1 + Z^2 / 2) / nrow(X)))
  rownames(out) <- c(targets, "se.mean", paste0("se.", p_names)); out
}
fit_stats <- function(M, P) {            # the targets implied by fitted parameters, one column per parameter set
  out <- rbind(M$mean(P), M$sd(P), rows(lapply(probs, function(p) M$q(p, P))), rows(lapply(slots, function(a) M$sf(a, P))))
  rownames(out) <- targets; out
}
fit_table <- function(M, x) {            # parametric bootstrap of one model: simulate from the ML fit and refit
  P <- M$fit(x)
  boot_table(fit_stats(M, P)[, 1], fit_stats(M, M$fit(boot_param(function(k) M$r(k, P), length(x), B))), alpha, efron = p_names)
}

# All methods on one sample x: estimate and interval per method and target, NA where a method
# is not used (M2 only where there is a standard error, M6 for the shares). The smoothed
# bootstrap is centred on the sd, quantiles and shares of the smoothed distribution, as in 11.
mc_one <- function(x) {
  n <- length(x); th <- edf_stats(x)[, 1]; est <- th[targets]; star <- edf_stats(boot_iid(x, B))
  h <- bw.nrd0(x); sm <- edf_stats(boot_smooth(x, B, h))
  centre <- c(sd = sqrt(mean((x - mean(x))^2) + h^2), setNames(sapply(probs, function(p) kde_quantile(x, p, h)), q_names),
              setNames(sapply(slots, function(a) mean(pnorm((x - a) / h))), p_names))
  sm[names(centre), ] <- sm[names(centre), ] - (centre - th[names(centre)])
  tabs <- list(M1 = boot_table(est, star, alpha),
               M2 = boot_table(th[c("mean", p_names, "se.mean", paste0("se.", p_names))], star, alpha),
               M3 = boot_table(normal_stats(x)[, 1], normal_stats(boot_param(rnorm, n, B, mean(x), sd(x))), alpha),
               M4 = fit_table(models$gamma, x),
               M5 = boot_table(th, sm, alpha),
               M6 = boot_table(est[p_names], star, alpha, efron = p_names),
               M7 = fit_table(models$weibull, x),
               M8 = fit_table(models$lognormal, x))
  out <- array(NA, c(length(methods), length(targets), 3), list(names(methods), targets, c("estimate", "lower", "upper")))
  for (M in names(tabs)) out[M, rownames(tabs[[M]]), ] <- as.matrix(tabs[[M]][, 1:3])
  out
}

# Steps 3 and 4: bias and RMSE of the estimates, coverage and average length of the intervals,
# and truth_above, the share of intervals whose upper limit lies below the truth (2.5% for an
# exact interval): an upper limit that is too low understates the risk of an overrun.
evaluate <- function(a, theta) {         # a: method x target x (estimate, lower, upper) x n_sim
  avg <- function(A) c(apply(A, 1:2, mean))
  err <- sweep(a[, , "estimate", ], 2, theta)
  data.frame(expand.grid(method = dimnames(a)[[1]], target = dimnames(a)[[2]]), true_value = rep(theta, each = dim(a)[1]),
             bias = avg(err), rmse = sqrt(avg(err^2)),
             coverage = avg(sweep(a[, , "lower", ], 2, theta, "<=") & sweep(a[, , "upper", ], 2, theta, ">=")),
             truth_above = avg(sweep(a[, , "upper", ], 2, theta, "<")), length = avg(a[, , "upper", ] - a[, , "lower", ]))
}
run_truth <- function(name) do.call(rbind, lapply(n_grid, function(n)
  data.frame(truth = name, n = n, evaluate(replicate(n_sim, mc_one(truths[[name]]$r(n))), truth_values[, name]))))
# Runs `one()` n_sim times split over the workers and returns the method x target x 3 x n_sim array.
run_split <- function(cl, one) {
  a <- parLapply(cl, lengths(splitIndices(n_sim, length(cl))), function(r, one) replicate(r, one()), one = one)
  array(unlist(a), c(dim(a[[1]])[1:3], n_sim), dimnames(a[[1]]))
}

# ---- Type 1 ---------------------------------------------------------------------------------
# 10_type1.R assumes normal durations and Poisson calls. With that model as the truth, at the
# Type 1 estimates (N(25.7, 5.84^2) minutes, 16.95 calls a day) and sample sizes (356 scans,
# 21 days, 355 gaps), does every interval of 10 cover 95%? Both of its bootstraps are tested;
# the plain percentile interval for p_work9 shows why 10 uses Efron's interval there.
mu1 <- 25.7; s1 <- 5.84; lam1 <- 16.95; slots1 <- c(25, 30, 35); n1 <- c(D = 356, N = 21, G = 355)
work9 <- function(lam, mu, s, k = 1:100)  # p_work9() of 10, one value per entry of lam, mu and s
  colSums(dpois(k, rep(lam, each = length(k))) * pnorm(540, k %o% mu, sqrt(k) %o% s, lower.tail = FALSE))
t1_stats <- function(D, N, G) {          # type1_quantities() of 10, one column per sample
  D <- as.matrix(D); N <- as.matrix(N); G <- as.matrix(G)
  mu <- colMeans(D); s <- col_sd(D); lam <- colMeans(N); lg <- 9 / colMeans(G); nD <- nrow(D)
  Z <- rows(lapply(slots1, function(a) (a - mu) / s))
  out <- rbind(lam, lg, mu, s, mu + qnorm(0.90) * s, mu + qnorm(0.95) * s, pnorm(Z, lower.tail = FALSE),
               lam * mu / 60, work9(lam, mu, s), sqrt(lam / nrow(N)), lg / sqrt(nrow(G)), s / sqrt(nD),
               dnorm(Z) * sqrt((1 + Z^2 / 2) / nD), sqrt(mu^2 * lam / nrow(N) + lam^2 * s^2 / nD) / 60)
  rownames(out) <- c("lambda", "lambda_gaps", "mu", "sigma", "q90", "q95", paste0("p_over_", slots1), "hours", "p_work9",
                     "se.lambda", "se.lambda_gaps", "se.mu", paste0("se.p_over_", slots1), "se.hours")
  out
}
truth1 <- c(lambda = lam1, lambda_gaps = lam1, mu = mu1, sigma = s1, q90 = mu1 + qnorm(0.90) * s1,
            q95 = mu1 + qnorm(0.95) * s1, setNames(pnorm(slots1, mu1, s1, lower.tail = FALSE), paste0("p_over_", slots1)),
            hours = lam1 * mu1 / 60, p_work9 = work9(lam1, mu1, s1))
mc_type1_one <- function() {
  D <- rnorm(n1[["D"]], mu1, s1); N <- rpois(n1[["N"]], lam1); G <- rexp(n1[["G"]], lam1 / 9)
  th <- t1_stats(D, N, G)[, 1]
  par <- t1_stats(boot_param(rnorm, length(D), B, th[["mu"]], th[["sigma"]]), boot_param(rpois, length(N), B, th[["lambda"]]),
                  boot_param(rexp, length(G), B, th[["lambda_gaps"]] / 9))
  np <- t1_stats(boot_iid(D, B), boot_iid(N, B), boot_iid(G, B))
  out <- array(NA, c(3, length(truth1), 3), list(c("parametric", "nonparametric", "parametric, percentile"), names(truth1),
                                                 c("estimate", "lower", "upper")))
  out[1, , ] <- as.matrix(boot_table(th, par, alpha, efron = "p_work9")[, 1:3])
  out[2, , ] <- as.matrix(boot_table(th, np, alpha, efron = "p_work9")[, 1:3])
  out[3, "p_work9", ] <- c(th[["p_work9"]], ci_percentile(th[["p_work9"]], par["p_work9", ], alpha))
  out
}

# ---- Type 2 gaps between calls ----------------------------------------------------------------
# 11 estimates the mean gap (working minutes), its CV and the calls per day from the gaps
# (lambda_gaps = 540 / mean gap) by resampling the gaps and by a parametric bootstrap of each
# model fitted to them. Truths: three shapes with the Type 2 mean gap and CV (54.0 minutes and
# 0.371) and 208 gaps, as in the data; gaps below 0.3 minutes are raised to 0.3 before a fit,
# as in 11. The exponential is left out: its CV is 1 by construction, against 0.371 in the data.
mg <- 54.0; cvg <- 0.371; ng <- 208
slg <- sqrt(log(1 + cvg^2))
kg <- uniroot(function(k) gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1 - cvg^2, c(1, 20), tol = 1e-10)$root
gap_truths <- list(gamma     = function(n) rgamma(n, 1 / cvg^2, 1 / (cvg^2 * mg)),
                   weibull   = function(n) rweibull(n, kg, mg / gamma(1 + 1 / kg)),
                   lognormal = function(n) rlnorm(n, log(mg) - slg^2 / 2, slg))
gap_targets <- c("mean_gap", "gap_cv", "lambda_gaps")
truth_gaps <- c(mean_gap = mg, gap_cv = cvg, lambda_gaps = 540 / mg)
gap_methods <- c(G1 = "resample gaps (tab11_type2_gaps)", G2 = "resample gaps, percentile-t (tab11_type2_arrivals)",
                 G3 = "gamma fit", G4 = "Weibull fit", G5 = "lognormal fit")
gap_np_stats <- function(X) {            # gap_np_stats() of 11, plus the delta-method se of lambda_gaps (gap_stats() of 11)
  X <- as.matrix(X); mx <- colMeans(X); sx <- col_sd(X)
  rbind(mean_gap = mx, se.mean_gap = sx / sqrt(nrow(X)), gap_cv = sx / mx, lambda_gaps = 540 / mx,
        se.lambda_gaps = 540 / mx * (sx / mx) / sqrt(nrow(X)))
}
gap_fit_table <- function(M, x) {        # tab_gap_mod of 11 for one model
  st <- function(P) rbind(mean_gap = M$mean(P), gap_cv = M$sd(P) / M$mean(P), lambda_gaps = 540 / M$mean(P))
  P <- M$fit(pmax(x, 0.3))
  boot_table(st(P)[, 1], st(M$fit(pmax(boot_param(function(k) M$r(k, P), length(x), B), 0.3))), alpha)
}
mc_gap_one <- function(x) {
  th <- gap_np_stats(x)[, 1]; star <- gap_np_stats(boot_iid(x, B))
  tabs <- list(G1 = boot_table(th[c("mean_gap", "se.mean_gap", "gap_cv", "lambda_gaps")], star, alpha),
               G2 = boot_table(th[c("lambda_gaps", "se.lambda_gaps")], star, alpha),
               G3 = gap_fit_table(models$gamma, x), G4 = gap_fit_table(models$weibull, x),
               G5 = gap_fit_table(models$lognormal, x))
  out <- array(NA, c(length(gap_methods), length(gap_targets), 3), list(names(gap_methods), gap_targets, c("estimate", "lower", "upper")))
  for (M in names(tabs)) out[M, rownames(tabs[[M]]), ] <- as.matrix(tabs[[M]][, 1:3])
  out
}

# ---- Type 2 days: calls, scan-hours and the 9-hour risk ------------------------------------------
# 21 days of calls at 10 a day on the working-hour clock, from a gamma renewal process as regular
# as Type 2 (shape 7) or a Poisson process, or with every day's number of calls drawn from the
# August Type 2 days (8 to 11 calls, more regular than the renewal process), and with gamma or
# lognormal scan durations of the Type 2 mean and sd. The estimators of 11, computed from the 21
# days and their scans:
#   lambda (calls per day):      D1 resample days, percentile; D2 resample days, percentile-t (as 11);
#                                D3 Poisson model, percentile-t (as 10)
#   dispersion, hours:           D1 (as 11 for the dispersion), D2 (as 11 for the hours)
#   p_work9 (counts and durations combined, p_day_over()): D4 Efron percentile, D5 bias-corrected (as 11)
#   share_days_<event>:          D6 exact binomial (as 11)
# The true values are those of a long stretch of the process (stationary): the share of 9-hour
# days and p_work9 estimate the same probability. p_work9 is not computed under the Poisson
# truth, which the Type 2 counts reject (dispersion 0.075), as it takes far longer there.
day_counts <- function(k, days = 21, rate = 10) tabulate(ceiling(cumsum(rgamma(4 * rate * days, k, k * rate / 9)) / 9), days)
obs_counts <- c(`8` = 1, `9` = 5, `10` = 8, `11` = 7) / 21    # the Type 2 days of August: 8 to 11 calls
day_truths <- list(
  `gamma renewal, gamma durations`     = list(k = 7, r = truths$gamma$r, p = truths$gamma$p, p9 = TRUE),
  `gamma renewal, lognormal durations` = list(k = 7, r = truths$lognormal$r, p = truths$lognormal$p, p9 = TRUE),
  `August counts, gamma durations`     = list(k = NA, r = truths$gamma$r, p = truths$gamma$p, p9 = TRUE),
  `Poisson, gamma durations`           = list(k = 1, r = truths$gamma$r, p = truths$gamma$p, p9 = FALSE))
draw_counts <- function(tr, days = 21)   # the August counts are drawn day by day from their own distribution
  if (is.na(tr$k)) sample(as.integer(names(obs_counts)), days, TRUE, obs_counts) else day_counts(tr$k, days)
cap <- 9 * 60 / slots                    # 12 and 9 calls: one machine's slots of 45 and 60 minutes
day_targets <- c("lambda", "dispersion", "hours", "p_work9", paste0("share_days_calls_over_", cap), "share_days_work9")
day_methods <- c(D1 = "resample days, percentile", D2 = "resample days, percentile-t", D3 = "Poisson model, percentile-t",
                 D4 = "resample days and scans, Efron", D5 = "resample days and scans, bias-corrected",
                 D6 = "exact binomial")
# True values. Counts of the renewal processes: a stationary stretch of 200,000 days; of the
# August truth: its own distribution. P(k scans take more than 540 minutes): by FFT on a
# 0.01-minute grid of the true duration distribution, as in p_day_over().
tail_sum <- function(cdf, k, limit = 540, step = 0.01, top = 250) {
  pm <- diff(cdf(c(0, seq(step / 2, top + step / 2, by = step))))   # mass of every grid point 0, step, ..., top
  L <- nextn(max(k) * length(pm)); f <- fft(c(pm, numeric(L - length(pm))))
  cw <- Conj(fft(pmin(pmax(seq_len(L) - 1 - limit / step + 0.5, 0), 1)))
  sapply(k, function(j) Re(sum(f^j * cw)) / L)
}
day_truth_values <- function(tr, days = 2e5) {
  N <- if (is.na(tr$k)) rep(as.integer(names(obs_counts)), round(obs_counts * 21)) else
    tabulate(ceiling(cumsum(rgamma(11 * days, tr$k, tr$k * 10 / 9)) / 9), days)[-1]   # drop the first day: stationary
  pk <- tabulate(N + 1) / length(N); k <- which(pk > 0) - 1
  p9 <- sum(pk[k + 1] * tail_sum(tr$p, k))
  mk <- sum(pk * (seq_along(pk) - 1)); disp <- sum(pk * (seq_along(pk) - 1 - mk)^2) / mk    # every truth has 10 calls a day
  c(lambda = 10, dispersion = disp, hours = 10 * m / 60, p_work9 = p9,
    setNames(sapply(cap, function(c) mean(N > c)), paste0("share_days_calls_over_", cap)), share_days_work9 = p9)
}
mc_days_one <- function(tr) {
  N <- draw_counts(tr); nd <- length(N); D <- tr$r(sum(N))
  H <- vapply(split(D, factor(rep(seq_len(nd), N), levels = seq_len(nd))), sum, 0) / 60     # scan-hours per day
  I <- boot_iid(seq_len(nd), B); Nm <- matrix(N[I], nd); Hm <- matrix(H[I], nd)
  lam <- mean(N); lb <- colMeans(Nm); lp <- colMeans(boot_param(rpois, nd, B, lam))
  out <- array(NA, c(length(day_methods), length(day_targets), 3), list(names(day_methods), day_targets, c("estimate", "lower", "upper")))
  out["D1", "lambda", ] <- c(lam, ci_percentile(lam, lb, alpha))
  out["D2", "lambda", ] <- c(lam, ci_percentile_t(lam, sd(N) / sqrt(nd), lb, col_sd(Nm) / sqrt(nd), alpha))
  out["D3", "lambda", ] <- c(lam, ci_percentile_t(lam, sqrt(lam / nd), lp, sqrt(lp / nd), alpha))
  out["D1", "dispersion", ] <- c(var(N) / lam, ci_percentile(var(N) / lam, col_sd(Nm)^2 / lb, alpha))
  out["D2", "hours", ] <- c(mean(H), ci_percentile_t(mean(H), sd(H) / sqrt(nd), colMeans(Hm), col_sd(Hm) / sqrt(nd), alpha))
  if (tr$p9) {
    p9 <- p_day_over(D, N); p9b <- p_day_over(boot_iid(D, B), Nm)
    out["D4", "p_work9", ] <- c(p9, ci_efron(p9b, alpha))
    out["D5", "p_work9", ] <- c(p9, ci_bc(p9, p9b, alpha))
  }
  events <- c(lapply(cap, function(c) N > c), list(H > 9)); names(events) <- tail(day_targets, 3)
  for (e in names(events)) {
    x <- sum(events[[e]]); out["D6", e, ] <- c(x / nd, binom.test(x, nd, conf.level = 1 - alpha)$conf.int)
  }
  out
}

# ---- Run --------------------------------------------------------------------------------------
cl <- makeCluster(4)                     # one worker per duration truth
clusterSetRNGStream(cl, 2026)            # independent random-number streams
# Each worker gets only what it uses: the functions of 20_boot_functions.R, sourced there, and
# the settings, truths and functions above. Not the whole workspace, which under run_all.R holds
# every result of 10 and 11, so the study runs the same on its own as from run_all.R.
invisible(clusterCall(cl, source, "R/20_boot_functions.R"))
clusterExport(cl, c("n_sim", "B", "alpha", "n_grid", "probs", "slots", "q_names", "p_names", "targets", "methods", "rows",
                    "m", "s", "sl", "kw", "truths", "truth_values", "col_quantiles", "edf_stats", "normal_stats",
                    "fit_stats", "fit_table", "mc_one", "evaluate", "run_truth",
                    "mu1", "s1", "lam1", "slots1", "n1", "work9", "t1_stats", "truth1", "mc_type1_one",
                    "mg", "cvg", "ng", "slg", "kg", "gap_truths", "truth_gaps", "gap_targets", "gap_methods", "gap_np_stats",
                    "gap_fit_table", "mc_gap_one",
                    "day_counts", "obs_counts", "draw_counts", "cap", "day_targets", "day_methods", "mc_days_one"))
tab30 <- do.call(rbind, parLapply(cl, names(truths), run_truth))
clusterSetRNGStream(cl, 2028)            # Type 1 results do not depend on the settings of the Type 2 part
a1 <- run_split(cl, mc_type1_one)
clusterSetRNGStream(cl, 2029)
tab34 <- do.call(rbind, parLapply(cl, names(gap_truths), function(name)
  data.frame(truth = name, n = ng, evaluate(replicate(n_sim, mc_gap_one(gap_truths[[name]](ng))), truth_gaps))))
set.seed(2027)                           # the true values of the day study, drawn here
truth_days <- sapply(day_truths, day_truth_values)
clusterSetRNGStream(cl, 2030)
tab31 <- do.call(rbind, lapply(names(day_truths), function(name) {
  tr <- day_truths[[name]]
  data.frame(truth = name, evaluate(run_split(cl, function() mc_days_one(tr)), truth_days[, name]))
}))
stopCluster(cl)
tab30 <- tab30[!is.na(tab30$coverage), ]
tab30$mc_se <- sqrt(tab30$coverage * (1 - tab30$coverage) / n_sim)   # Monte Carlo standard error of the coverage
tab33 <- evaluate(a1, truth1)
tab33 <- tab33[!is.na(tab33$coverage), ]
tab33$mc_se <- sqrt(tab33$coverage * (1 - tab33$coverage) / n_sim)
tab34 <- tab34[!is.na(tab34$coverage), ]
tab34$mc_se <- sqrt(tab34$coverage * (1 - tab34$coverage) / n_sim)
tab31 <- tab31[!is.na(tab31$coverage), ]
tab31$mc_se <- sqrt(tab31$coverage * (1 - tab31$coverage) / n_sim)

# The report tables: coverage (%) at our sample size, n = 209, for every target and method
# (tab32), and for the quantities that set the slots, q95 and P(scan > 60 min), the model fits
# next to the intervals of 11's main table under every truth (tab35): a model fit is only safe
# when its model is the truth.
tab32 <- with(subset(tab30, n == 209), tapply(round(100 * coverage, 1), list(paste(target, method), factor(truth, names(truths))), sum))
tab35 <- with(subset(tab30, n == 209 & target %in% c("q95", "p_over_60") &
                       (method %in% c("M4", "M7", "M8") | (target == "q95" & method == "M5") | (target == "p_over_60" & method == "M2"))),
              tapply(round(100 * coverage, 1), list(paste(target, methods[as.character(method)]), factor(truth, names(truths))), sum))
print(tab32)
cat("\nModel fits against the main-table intervals of 11, coverage (%) at n = 209:\n"); print(tab35)
cat("\nType 1, coverage when the model of 10_type1.R is the truth:\n")
print(data.frame(tab33[1:2], round(tab33[c("true_value", "coverage", "truth_above")], 3)), row.names = FALSE)
cat("\nType 2 gaps between calls (n = 208):\n")
print(data.frame(tab34[c("truth", "method", "target")], round(tab34[c("true_value", "bias", "coverage", "truth_above")], 3)),
      row.names = FALSE)
cat("\nType 2 days (21 days):\n")
print(data.frame(tab31[c("truth", "method", "target")], round(tab31[c("true_value", "bias", "coverage", "truth_above", "length")], 4)),
      row.names = FALSE)

# ---- Output ---------------------------------------------------------------------------------
write.csv(tab30, "output/tab30_mc_durations.csv", row.names = FALSE)
write.csv(tab31, "output/tab31_mc_days.csv", row.names = FALSE)
write.csv(tab32, "output/tab32_mc_coverage_n209.csv")
write.csv(tab33, "output/tab33_mc_type1.csv", row.names = FALSE)
write.csv(tab34, "output/tab34_mc_gaps.csv", row.names = FALSE)
write.csv(tab35, "output/tab35_mc_models_n209.csv")
saveRDS(list(durations = tab30, days = tab31, coverage_n209 = tab32, type1 = tab33, gaps = tab34, models_n209 = tab35,
             truth_values = truth_values, truth_days = truth_days, truth_gaps = truth_gaps, truth_type1 = truth1,
             n_sim = n_sim, B = B), "output/mc_study.rds")

# Figure 30: coverage against n per method (matplot, as in BootstrapCode.R), one panel per
# target and truth; dashed: 95%, dotted: two Monte Carlo standard errors either side.
pdf("output/fig30_mc_coverage.pdf", width = 10, height = 13)
par(mfrow = c(length(targets), length(truths)), mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.6, 0))
for (tg in targets) for (tn in names(truths)) {
  d <- subset(tab30, target == tg & truth == tn)
  matplot(n_grid, tapply(d$coverage, list(d$n, d$method), sum), type = "b", lty = 1, pch = 16, col = seq_along(methods),
          ylim = c(0, 1), log = "x", xlab = "n", ylab = "Coverage", main = paste(tg, "|", tn), cex.main = 0.9)
  abline(h = 0.95 + c(-2, 0, 2) * sqrt(0.95 * 0.05 / n_sim), lty = c(3, 2, 3))
  if (tg == targets[1] && tn == names(truths)[1])
    legend("bottomleft", paste(names(methods), methods), col = seq_along(methods), lty = 1, pch = 16, bty = "n", cex = 0.6)
}
invisible(dev.off())
