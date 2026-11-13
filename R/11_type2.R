# 11_type2.R
# Type 2 patients. The hospital does not know the distributions of their scan durations or
# arrivals, so the quantities it cares about are estimated by plug-in (no model assumed),
# with uncertainty from bootstraps of the durations:
#   - non-parametric: resample the observed durations (boot_iid)
#   - smoothed:       resample and add N(0, h^2) noise (boot_smooth), centred on the KDE
#   - parametric:     for each candidate model (exponential, gamma, Weibull, lognormal),
#                     simulate from its ML fit and refit every sample (boot_param)
# Arrivals (calls and scan-hours per day, gaps between calls) get the non-parametric bootstrap.
# A parametric-bootstrap KS test checks which of these models (and the normal) fits, both for the
# durations and for the time between calls; the gaps also get a non-parametric and a parametric
# bootstrap per model.
# Names follow 10_type1.R (lambda, lambda_gaps, hours, p_work9, p_over_<slot>).
# Output: output/tab11_type2_{durations,arrivals,gaps,gof,gof_gaps}.csv, output/fig12_type2_bootstrap.pdf,
# fig13-14 (durations) and fig15-16 (time between calls): histograms and QQ plots of the fitted
# models, as fig11 for Type 1; figures are also drawn in RStudio's Plots pane.
# output/type2.rds, and `type2_inputs` (generators in minutes) for the simulation in Part 2.

if (!exists("scans")) source("R/00_load_clean.R")
source("R/20_boot_functions.R")

set.seed(2025)
B     <- 999                 # bootstrap samples: 999 gives exact 2.5% / 97.5% order statistics
alpha <- 0.05                # 95% intervals
probs <- c(0.90, 0.95)       # quantiles: slot length that fits 90% / 95% of the scans
slots <- c(45, 60)           # candidate slot lengths in minutes (dotted lines in fig02)
cap   <- 9 * 60 / slots      # Type 2 slots in one 9-hour machine day: 12 and 9

D  <- scans$dur_min[scans$type == 2]   # scan durations in minutes
n  <- length(D)
G  <- gaps[["2"]]                        # gaps between calls in working hours
N  <- daily$calls2                      # calls per working day
H  <- daily$hours2                      # scan-hours per working day
nd <- nrow(daily)

q_names <- paste0("q", 100 * probs)     # q90, q95
p_names <- paste0("p_over_", slots)     # p_over_45, p_over_60: P(scan longer than the slot)
rows    <- function(L) do.call(rbind, L)

# ---- Scan durations ---------------------------------------------------------------------
# Plug-in statistics for every column of X (one bootstrap sample per column).
dur_stats <- function(X) {
  X <- as.matrix(X); s <- col_sd(X)
  out <- rbind(mean = colMeans(X), se.mean = s / sqrt(nrow(X)), sd = s,
               rows(lapply(probs, function(p) col_quantile(X, p))),
               rows(lapply(slots, function(a) colMeans(X > a))))
  rownames(out) <- c("mean", "se.mean", "sd", q_names, p_names)
  out
}
# Candidate models for the durations. Each has an ML fit for every column of a matrix (one
# parameter column per sample), a random generator, quantile and survival functions, and its
# mean and sd. Exponential and lognormal ML are closed form, the gamma uses gamma_ml from
# 20_boot_functions.R and the Weibull weibull_ml below (kept here, not in the shared file 20).
# The normal is included only in the fit test, as the Type 1 assumption.
# Weibull ML: the shape k solves 1/k + mean(log x) - sum(x^k log x) / sum(x^k) = 0 (Newton from
# the usual CV-based start), then scale = mean(x^k)^(1/k); it matches MASS::fitdistr.
weibull_ml <- function(X) {
  X <- as.matrix(X); L <- log(X); mL <- colMeans(L)
  k <- (col_sd(X) / colMeans(X))^-1.086
  for (i in 1:20) {
    Xk <- X^rep(k, each = nrow(X)); S0 <- colSums(Xk); S1 <- colSums(Xk * L); S2 <- colSums(Xk * L^2)
    k <- k - (1 / k + mL - S1 / S0) / (-1 / k^2 - (S2 / S0 - (S1 / S0)^2))
  }
  rbind(shape = k, scale = colMeans(X^rep(k, each = nrow(X)))^(1 / k))
}
models <- list(
  exponential = list(fit  = function(X) rbind(rate = 1 / colMeans(as.matrix(X))),
                     r    = function(m, p) rexp(m, p[1, ]),
                     d    = function(x, p) dexp(x, p[1, ]),
                     q    = function(pr, p) qexp(pr, p[1, ]),
                     sf   = function(a, p) pexp(a, p[1, ], lower.tail = FALSE),
                     mean = function(p) 1 / p[1, ], sd = function(p) 1 / p[1, ]),
  gamma       = list(fit  = gamma_ml,
                     r    = function(m, p) rgamma(m, p[1, ], p[2, ]),
                     d    = function(x, p) dgamma(x, p[1, ], p[2, ]),
                     q    = function(pr, p) qgamma(pr, p[1, ], p[2, ]),
                     sf   = function(a, p) pgamma(a, p[1, ], p[2, ], lower.tail = FALSE),
                     mean = function(p) p[1, ] / p[2, ], sd = function(p) sqrt(p[1, ]) / p[2, ]),
  weibull     = list(fit  = weibull_ml,
                     r    = function(m, p) rweibull(m, p[1, ], p[2, ]),
                     d    = function(x, p) dweibull(x, p[1, ], p[2, ]),
                     q    = function(pr, p) qweibull(pr, p[1, ], p[2, ]),
                     sf   = function(a, p) pweibull(a, p[1, ], p[2, ], lower.tail = FALSE),
                     mean = function(p) p[2, ] * gamma(1 + 1 / p[1, ]),
                     sd   = function(p) p[2, ] * sqrt(gamma(1 + 2 / p[1, ]) - gamma(1 + 1 / p[1, ])^2)),
  lognormal   = list(fit  = function(X) { L <- log(as.matrix(X)); m <- colMeans(L)
  rbind(meanlog = m, sdlog = sqrt(colMeans(sweep(L, 2, m)^2))) },
  r    = function(m, p) rlnorm(m, p[1, ], p[2, ]),
  d    = function(x, p) dlnorm(x, p[1, ], p[2, ]),
  q    = function(pr, p) qlnorm(pr, p[1, ], p[2, ]),
  sf   = function(a, p) plnorm(a, p[1, ], p[2, ], lower.tail = FALSE),
  mean = function(p) exp(p[1, ] + p[2, ]^2 / 2),
  sd   = function(p) exp(p[1, ] + p[2, ]^2 / 2) * sqrt(exp(p[2, ]^2) - 1)))
normal <- list(fit = function(X) rbind(mean = colMeans(as.matrix(X)), sd = col_sd(X)),
               r   = function(m, p) rnorm(m, p[1, ], p[2, ]),
               d   = function(x, p) dnorm(x, p[1, ], p[2, ]),
               q   = function(pr, p) qnorm(pr, p[1, ], p[2, ]),
               sf  = function(a, p) pnorm(a, p[1, ], p[2, ], lower.tail = FALSE))

# The hospital quantities implied by fitted parameters (one column per parameter set).
model_stats <- function(M, par) {
  out <- rbind(M$mean(par), M$sd(par), rows(lapply(probs, function(pr) M$q(pr, par))),
               rows(lapply(slots, function(a) M$sf(a, par))))
  rownames(out) <- c("mean", "sd", q_names, p_names)
  out
}

th <- dur_stats(D)[, 1]

# 1. Non-parametric bootstrap: percentile-t for the mean, percentile for sd and quantiles,
#    Efron's percentile for the tail probabilities (skewed, bounded by 0 and 1).
star_np <- dur_stats(boot_iid(D, B))
tab_np  <- boot_table(th, star_np, alpha, efron = p_names)

# 2. Smoothed bootstrap. Its samples come from the EDF convolved with N(0, h^2), whose sd,
#    quantiles and tail probabilities differ from the plug-in estimates. Shifting the
#    bootstrap estimates by (smoothed value - estimate) centres them on the right value, which
#    is what ci_percentile's `center` does; the mean needs no shift.
h <- bw.nrd0(D)
star_sm <- dur_stats(boot_smooth(D, B, h))
centre <- c(sd = sqrt(mean((D - mean(D))^2) + h^2),
            setNames(sapply(probs, function(p) kde_quantile(D, p, h)), q_names),
            setNames(sapply(slots, function(a) mean(pnorm((D - a) / h))), p_names))
star_sm[names(centre), ] <- star_sm[names(centre), ] - (centre - th[names(centre)])
tab_sm <- boot_table(th, star_sm, alpha, efron = p_names)

# 3. Parametric bootstraps, one per model: estimates are the model-implied quantities, and
#    each bootstrap sample is simulated from the fit and refitted by ML. If a model is right,
#    its intervals are narrower than the non-parametric ones; if it is wrong (the exponential,
#    with CV 1, clearly is) its estimates of the tail quantities are biased.
fits     <- lapply(models, function(M) M$fit(D))
par_star <- lapply(names(models), function(m) models[[m]]$fit(boot_param(function(k) models[[m]]$r(k, fits[[m]]), n, B)))
names(par_star) <- names(models)
tab_mod  <- lapply(names(models), function(m)
  boot_table(model_stats(models[[m]], fits[[m]])[, 1], model_stats(models[[m]], par_star[[m]]), alpha, efron = p_names))
names(tab_mod) <- names(models)

# Long table: one row per quantity and method (plug-in with the non-parametric and smoothed
# intervals, then each model's estimate with its parametric interval).
long <- function(tab, method) data.frame(quantity = rownames(tab), method = method, tab, row.names = NULL)
tab_dur <- do.call(rbind, c(list(long(tab_np, "non-parametric"), long(tab_sm, "smoothed")),
                            lapply(names(tab_mod), function(m) long(tab_mod[[m]], m))))
tab_dur <- tab_dur[order(match(tab_dur$quantity, rownames(tab_np))), ]
rownames(tab_dur) <- NULL

# ---- Goodness of fit: KS distance to each fitted model ------------------------------------
# The parameters are estimated, so the standard KS p-value is wrong (too large); the null
# distribution comes from a parametric bootstrap that refits the parameters on every sample.
ks <- function(x, cdf) { x <- sort(x); m <- length(x); F <- cdf(x); max((1:m) / m - F, F - (0:(m - 1)) / m) }
ks_model <- function(M, x) {
  stat <- function(y) { p <- M$fit(y); ks(y, function(q) 1 - M$sf(q, p)) }
  obs  <- stat(x)
  X    <- boot_param(function(k) M$r(k, M$fit(x)), length(x), B)
  c(KS = obs, p_value = (1 + sum(apply(X, 2, stat) >= obs)) / (B + 1))
}
all_models <- c(models, normal = list(normal))
all_fits   <- lapply(all_models, function(M) M$fit(D))
tab_gof    <- t(sapply(all_models, ks_model, x = D))

# ---- Arrivals ------------------------------------------------------------------------------
# Days are resampled whole, so a day's calls and scan-hours stay together; gaps are resampled
# on their own (as in 10_type1.R). Names as in Type 1: lambda = calls per day, lambda_gaps =
# calls per day from the gaps (9 / mean gap), hours = scan-hours per day, p_work9 = P(a day's
# Type 2 scans take more than 9 h). Without a model these are plug-in: p_work9 is the share of
# such days. dispersion (variance / mean) is 1 for Poisson counts, gap_cv 1 for exponential
# gaps. p_calls_over_<cap>: the chance a day has more calls than one machine has 45- or 60-minute slots.
day_stats <- function(I) {
  I <- as.matrix(I); Nm <- matrix(N[I], nrow(I)); Hm <- matrix(H[I], nrow(I))
  sN <- col_sd(Nm); sH <- col_sd(Hm)
  out <- rbind(colMeans(Nm), sN^2 / colMeans(Nm), rows(lapply(cap, function(c) colMeans(Nm > c))),
               colMeans(Hm), colMeans(Hm > 9), sN / sqrt(nrow(I)), sH / sqrt(nrow(I)))
  rownames(out) <- c("lambda", "dispersion", paste0("p_calls_over_", cap), "hours", "p_work9", "se.lambda", "se.hours")
  out
}
gap_stats <- function(X) {
  X <- as.matrix(X); m <- colMeans(X); cv <- col_sd(X) / m
  rbind(lambda_gaps = 9 / m, gap_cv = cv, se.lambda_gaps = 9 / m * cv / sqrt(nrow(X)))   # delta method
}
th_arr   <- c(day_stats(seq_len(nd))[, 1], gap_stats(G)[, 1])
star_arr <- rbind(day_stats(boot_iid(seq_len(nd), B)), gap_stats(boot_iid(G, B)))
tab_arr  <- boot_table(th_arr, star_arr, alpha, efron = c(paste0("p_calls_over_", cap), "p_work9"))

# Time between calls: the same four models as for the durations, fitted by ML to the gaps in
# working minutes, with the KS test (figures 15 and 16 show the fit). Exponential gaps mean
# Poisson arrivals, as assumed for Type 1; gamma and Weibull contain the exponential as shape 1,
# so their fitted shape shows how far the gaps are from it (below 1: more bunched calls).
# Calls in the same 0.01 h give a gap of 0, where log-based fits break down; for fitting those
# are set to half the recording unit (0.3 minutes).
Gm        <- pmax(60 * G, 0.3)
gap_fits  <- lapply(models, function(M) M$fit(Gm))
tab_gof_gaps <- t(sapply(models, ks_model, x = Gm))

# Bootstraps of the time between calls, in the layout of the durations table: the mean gap in
# working minutes, its CV (1 for exponential gaps) and the implied calls per working day
# (540 / mean gap, as lambda_gaps). Non-parametric: resample the observed gaps. Parametric, one
# per model: simulate gaps from the ML fit and refit (as Type 1 does with the exponential).
# Under the exponential the CV is 1 by construction, so its interval is [1, 1].
gap_np_stats <- function(X) {
  X <- as.matrix(X); m <- colMeans(X); s <- col_sd(X)
  rbind(mean_gap = m, se.mean_gap = s / sqrt(nrow(X)), gap_cv = s / m, lambda_gaps = 540 / m)
}
gap_model_stats <- function(M, par) {
  m <- M$mean(par)
  rbind(mean_gap = m, gap_cv = M$sd(par) / m, lambda_gaps = 540 / m)
}
th_gap      <- gap_np_stats(60 * G)[, 1]
tab_gap_np  <- boot_table(th_gap, gap_np_stats(boot_iid(60 * G, B)), alpha)
gap_par_star <- lapply(names(models), function(m)
  models[[m]]$fit(pmax(boot_param(function(k) models[[m]]$r(k, gap_fits[[m]]), length(Gm), B), 0.3)))
names(gap_par_star) <- names(models)
tab_gap_mod <- lapply(names(models), function(m)
  boot_table(gap_model_stats(models[[m]], gap_fits[[m]])[, 1], gap_model_stats(models[[m]], gap_par_star[[m]]), alpha))
names(tab_gap_mod) <- names(models)
tab_gaps <- do.call(rbind, c(list(long(tab_gap_np, "non-parametric")),
                             lapply(names(tab_gap_mod), function(m) long(tab_gap_mod[[m]], m))))
tab_gaps <- tab_gaps[order(match(tab_gaps$quantity, rownames(tab_gap_np))), ]
rownames(tab_gaps) <- NULL

# Busy-day range: 5% and 95% of a day's calls, from the data itself and, for comparison, from
# a Poisson(lambda) averaged over the bootstrap lambdas (what Type 1 assumes).
k <- 0:60
F_mix <- sapply(k, function(j) mean(ppois(j, star_arr["lambda", ])))
busy <- rbind(plug_in = c(col_quantile(N, 0.05), col_quantile(N, 0.95)),
              poisson_bootstrap = k[c(which(F_mix >= 0.05)[1], which(F_mix >= 0.95)[1])],
              observed_min_max = range(N))
colnames(busy) <- c("5%", "95%")

# Diagnostics, as for Type 1: dispersion test of the daily counts (chi-square, nd - 1 df),
# lag-1 autocorrelations in call order (iid check), and the KS fit tests from above.
X2 <- sum((N - mean(N))^2) / mean(N)
diagnostics <- c(dispersion_X2 = X2, p_over_dispersed = pchisq(X2, nd - 1, lower.tail = FALSE),
                 p_under_dispersed = pchisq(X2, nd - 1),
                 acf1_durations = acf(D, plot = FALSE)$acf[2], acf1_gaps = acf(G, plot = FALSE)$acf[2],
                 setNames(tab_gof[, "p_value"], paste0("ks_p_", rownames(tab_gof))))

# ---- Output ---------------------------------------------------------------------------------
round_df <- function(d) { num <- sapply(d, is.numeric); d[num] <- lapply(d[num], round, 3); d }
cat(sprintf("Type 2: %d scans, %.1f calls per day; bandwidth h = %.2f min.\n", n, mean(N), h))
cat("\nML fits of the duration models:\n"); for (m in names(fits)) cat(sprintf("  %-12s %s\n", m,
                                                                               paste(sprintf("%s = %.3f", rownames(fits[[m]]), fits[[m]][, 1]), collapse = ", ")))
cat("\nScan durations (minutes; p_over = share of scans longer than the slot):\n"); print(round_df(tab_dur))
cat("\nGoodness of fit of the durations (parametric-bootstrap KS test):\n"); print(round(tab_gof, 3))
cat("\nTime between calls (working minutes): ML fits and KS test:\n")
for (m in names(gap_fits)) cat(sprintf("  %-12s %s\n", m,
                                       paste(sprintf("%s = %.3f", rownames(gap_fits[[m]]), gap_fits[[m]][, 1]), collapse = ", ")))
print(round(tab_gof_gaps, 3))
cat("\nTime between calls: mean gap (working minutes), CV and calls per day, per method:\n"); print(round_df(tab_gaps))
cat("\nArrivals (per working day; lambda_gaps from the gaps between calls):\n"); print(round_df(tab_arr))
cat("\nBusy-day range (calls per day):\n"); print(busy)
cat("\nDiagnostics:\n"); print(round(diagnostics, 3))
write.csv(round_df(tab_dur), "output/tab11_type2_durations.csv", row.names = FALSE)
write.csv(round_df(tab_arr), "output/tab11_type2_arrivals.csv")
write.csv(round(tab_gof, 3), "output/tab11_type2_gof.csv")
write.csv(round(tab_gof_gaps, 3), "output/tab11_type2_gof_gaps.csv")
write.csv(round_df(tab_gaps), "output/tab11_type2_gaps.csv", row.names = FALSE)

# ---- Figures ----------------------------------------------------------------------------------
# Every figure is written to output/ as a PDF and, when the script runs in RStudio, also drawn
# in the Plots pane (use the arrows there to page through them).
show_fig <- function(file, width, height, draw) {
  pdf(file, width = width, height = height); draw(); invisible(dev.off())
  if (interactive()) draw()
}
mcols <- c(exponential = "firebrick", gamma = "darkorange3", weibull = "darkgreen",
           lognormal = "purple", normal = "steelblue4")

# Figure 12: bootstrap distributions of the 95% quantile (slot length) and of P(scan > 60 min)
# under every bootstrap, with the plug-in estimate as a dashed line. The exponential is left out:
# it is far off (see the tables and the fit test) and would squash the other curves.
in_fig <- setdiff(names(models), "exponential")
cols <- c("grey40", "steelblue4", mcols[in_fig])
labs <- c("non-parametric", "smoothed", in_fig)
show_fig("output/fig12_type2_bootstrap.pdf", 8, 4, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  for (q in c("q95", "p_over_60")) {
    dens <- lapply(c(list(star_np[q, ], star_sm[q, ]),
                     lapply(in_fig, function(m) model_stats(models[[m]], par_star[[m]])[q, ])), density)
    plot(NULL, xlim = range(sapply(dens, function(d) range(d$x))), ylim = c(0, max(sapply(dens, function(d) max(d$y)))),
         main = if (q == "q95") "95% quantile of scan duration" else "P(scan longer than 60 min)",
         xlab = if (q == "q95") "Minutes" else "Probability", ylab = "Bootstrap density")
    for (i in seq_along(dens)) lines(dens[[i]], col = cols[i], lwd = 2)
    abline(v = th[q], lty = 2)
    legend("topright", labs, col = cols, lwd = 2, bty = "n", cex = 0.7)
  }
})

# Model fit, the visual check next to the KS test (as fig11_type1_diagnostics.pdf for Type 1).
# Histogram: the data with each fitted density, one panel per model and a last panel with all.
# QQ: observed quantiles against each model's; points on the line mean the model fits, a bend
# at the right end means it gets the long values wrong.
fit_hist <- function(x, mods, fits, gof, binw, xlab) {
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  hb <- hist(x, breaks = seq(0, binw * ceiling(max(x) / binw), binw), plot = FALSE)
  xs <- seq(0.01 * binw, max(hb$breaks), length.out = 400)
  ymax <- min(2 * max(hb$density), 1.15 * max(hb$density, sapply(names(mods), function(m) max(mods[[m]]$d(xs, fits[[m]])))))
  for (m in c(names(mods), "all")) {
    plot(hb, freq = FALSE, ylim = c(0, ymax), col = "grey85", border = "white", xlab = xlab, ylab = "Density",
         main = if (m == "all") "All models" else sprintf("%s (KS p = %.3f)", m, gof[m, "p_value"]))
    for (j in if (m == "all") names(mods) else m) lines(xs, mods[[j]]$d(xs, fits[[j]]), col = mcols[[j]], lwd = 2)
    if (m == "all") legend("topright", names(mods), col = mcols[names(mods)], lwd = 2, bty = "n", cex = 0.8)
  }
}
fit_qq <- function(x, mods, fits, gof, unit) {
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  pp <- ppoints(length(x))
  for (m in names(mods)) {
    qqplot(mods[[m]]$q(pp, fits[[m]]), x, main = sprintf("%s (KS p = %.3f)", m, gof[m, "p_value"]),
           pch = 16, cex = 0.5, col = mcols[[m]], xlab = paste0("Model quantiles (", unit, ")"),
           ylab = paste0("Observed (", unit, ")"))
    abline(0, 1)
  }
}

# Figures 13 and 14: scan durations, five models (normal included as the Type 1 assumption).
show_fig("output/fig13_type2_dur_histograms.pdf", 9, 6, function()
  fit_hist(D, all_models, all_fits, tab_gof, 5, "Scan duration (minutes)"))
show_fig("output/fig14_type2_dur_qq.pdf", 9, 6, function()
  fit_qq(D, all_models, all_fits, tab_gof, "minutes"))

# Figures 15 and 16: time between calls (working minutes), four models.
show_fig("output/fig15_type2_gap_histograms.pdf", 9, 6, function()
  fit_hist(Gm, models, gap_fits, tab_gof_gaps, 15, "Time between calls (working minutes)"))
show_fig("output/fig16_type2_gap_qq.pdf", 9, 6, function()
  fit_qq(Gm, models, gap_fits, tab_gof_gaps, "minutes"))

# Hand-over to Part 2 in the same format as type1_inputs: minutes, gaps on the working-minute
# clock (540 a day). dur2 draws from the smoothed EDF (no model assumed), dur2_<model> from
# each fitted model (e.g. dur2_gamma); gap2 resamples the observed gaps, gap2_<model> draws
# from the models fitted to the gaps. Values are fixed here, as in Type 1.
type2_inputs <- local({
  d <- D; hh <- h; g <- 60 * G; ms <- models; fs <- fits; gf <- gap_fits
  c(list(dur2 = function(n) pmax(sample(d, n, replace = TRUE) + hh * rnorm(n), 1),
         gap2 = function(n) sample(g, n, replace = TRUE)),
    setNames(lapply(names(ms), function(m) { M <- ms[[m]]; p <- fs[[m]]; function(n) M$r(n, p) }),
             paste0("dur2_", names(ms))),
    setNames(lapply(names(ms), function(m) { M <- ms[[m]]; p <- gf[[m]]; function(n) M$r(n, p) }),
             paste0("gap2_", names(ms))))
})

saveRDS(list(estimates = th, durations = tab_dur, arrivals = tab_arr, gof = tab_gof, gof_gaps = tab_gof_gaps,
             gap_fits = lapply(gap_fits, function(p) p[, 1]), gaps = tab_gaps,
             gap_par_star = lapply(gap_par_star, t), busy_day = busy,
             diagnostics = diagnostics, np_star = star_np, sm_star = star_sm, arr_star = star_arr,
             fits = lapply(fits, function(p) p[, 1]), par_star = lapply(par_star, t),   # B parameter draws per model
             inputs = type2_inputs), "output/type2.rds")
