# 30_mc_study.R
# Can the hospital trust our intervals? Invent a truth, draw fake data from it, let every method
# estimate without seeing the truth, and count how often its 95% interval catches the truth:
# simulate, apply, evaluate, summarise, as in BootstrapCode.R. The methods are the interval
# rules of 10_type1.R and 11_type2.R, built with the same boot_table() calls:
#   Type 2 durations, four invented truths with the Type 2 mean and sd, n = 50 to 500:
#     M1 EDF, percentile          what 11 uses for q90 and q95
#     M2 EDF, percentile-t        what 11 uses for the mean; the plan's rule for the shares
#     M3 normal fit               the Type 1 recipe of 10 applied to Type 2 (the naive analyst)
#     M4 gamma fit                as in 11, with Efron's interval for the shares
#     M5 smoothed bootstrap       as in 11, for the quantiles
#     M6 EDF, Efron percentile    what 11 uses for the shares
#   Type 1: every quantity of 10 when its own model is the truth, at its sample sizes.
#   Type 2 arrivals: the mean daily count from 21 days, resampling days (as 11) or assuming
#     Poisson (as 10), under a gamma renewal truth and a Poisson truth.
# The four duration truths run in parallel, one worker each; the Type 1 study is split over the
# same four workers.
# Output: output/tab30-33_mc_*.csv, output/mc_study.rds, output/fig30_mc_coverage.pdf.

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
targets <- c("mean", q_names, p_names)
methods <- c(M1 = "EDF, percentile", M2 = "EDF, percentile-t", M3 = "normal fit", M4 = "gamma fit",
             M5 = "smoothed bootstrap", M6 = "EDF, Efron percentile")
rows <- function(L) do.call(rbind, L)

# ---- Type 2 durations -----------------------------------------------------------------------
# Step 1, the truths: four shapes with the Type 2 mean and sd (40.7 and 12.1 minutes). The
# abs() only guards the mixture against a negative draw (probability below 1e-8).
m <- 40.7; s <- 12.1
sl <- sqrt(log(1 + (s / m)^2))                                                                    # lognormal sdlog
kw <- uniroot(function(k) gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1 - (s / m)^2, c(1, 20), tol = 1e-10)$root  # Weibull shape
truths <- list(
  gamma     = list(r = function(n) rgamma(n, (m / s)^2, m / s^2), p = function(x) pgamma(x, (m / s)^2, m / s^2)),
  lognormal = list(r = function(n) rlnorm(n, log(m) - sl^2 / 2, sl), p = function(x) plnorm(x, log(m) - sl^2 / 2, sl)),
  weibull   = list(r = function(n) rweibull(n, kw, m / gamma(1 + 1 / kw)), p = function(x) pweibull(x, kw, m / gamma(1 + 1 / kw))),
  mixture   = list(r = function(n) abs(ifelse(runif(n) < 0.7, rnorm(n, 34.1, 6), rnorm(n, 56.1, 8))),
                   p = function(x) 0.7 * pnorm(x, 34.1, 6) + 0.3 * pnorm(x, 56.1, 8)))
true_q <- function(tr, p) uniroot(function(x) tr$p(x) - p, c(1, 200), tol = 1e-10)$root
truth_values <- sapply(truths, function(tr) c(mean = m, setNames(sapply(probs, function(p) true_q(tr, p)), q_names),
                                              setNames(1 - tr$p(slots), p_names)))

# Step 2, the methods. Statistics for every column of X (one sample per column), with the
# standard errors ("se.") that make boot_table() choose percentile-t: the EDF as dur_stats()
# in 11 (plus se_prop() for the shares), the normal fit as type1_quantities() in 10, and the
# gamma fit as model_stats() in 11. col_quantiles() gives the same order statistics as
# col_quantile() for both probabilities, from one sort of the whole matrix.
col_quantiles <- function(X, p) { X <- as.matrix(X); matrix(X[order(col(X), X)], nrow(X))[ceiling(nrow(X) * p), , drop = FALSE] }
edf_stats <- function(X) {
  X <- as.matrix(X); P <- rows(lapply(slots, function(a) colMeans(X > a)))
  out <- rbind(colMeans(X), col_quantiles(X, probs), P, col_sd(X) / sqrt(nrow(X)), se_prop(P, nrow(X)))
  rownames(out) <- c(targets, "se.mean", paste0("se.", p_names)); out
}
normal_stats <- function(X) {
  X <- as.matrix(X); mu <- colMeans(X); sdv <- col_sd(X); Z <- rows(lapply(slots, function(a) (a - mu) / sdv))
  out <- rbind(mu, rows(lapply(probs, function(p) mu + qnorm(p) * sdv)), pnorm(Z, lower.tail = FALSE),
               sdv / sqrt(nrow(X)), dnorm(Z) * sqrt((1 + Z^2 / 2) / nrow(X)))
  rownames(out) <- c(targets, "se.mean", paste0("se.", p_names)); out
}
gamma_stats <- function(X) {
  g <- gamma_ml(X)
  out <- rbind(g[1, ] / g[2, ], rows(lapply(probs, function(p) qgamma(p, g[1, ], g[2, ]))),
               rows(lapply(slots, function(a) pgamma(a, g[1, ], g[2, ], lower.tail = FALSE))))
  rownames(out) <- targets; out
}

# All methods on one sample x: estimate and interval per method and target, NA where a method
# is not used (M2 only where there is a standard error, M5 for quantiles, M6 for the shares).
# The smoothed quantiles are centred on those of the smoothed distribution, as in 11.
mc_one <- function(x) {
  n <- length(x); th <- edf_stats(x)[, 1]; est <- th[targets]; star <- edf_stats(boot_iid(x, B))
  h <- bw.nrd0(x); sm <- col_quantiles(boot_smooth(x, B, h), probs)
  sm <- sm - (sapply(probs, function(p) kde_quantile(x, p, h)) - est[q_names]); rownames(sm) <- q_names
  g <- gamma_ml(x)
  tabs <- list(M1 = boot_table(est, star, alpha),
               M2 = boot_table(th[!names(th) %in% q_names], star, alpha),
               M3 = boot_table(normal_stats(x)[, 1], normal_stats(boot_param(rnorm, n, B, mean(x), sd(x))), alpha),
               M4 = boot_table(gamma_stats(x)[, 1], gamma_stats(boot_param(rgamma, n, B, g[1], g[2])), alpha, efron = p_names),
               M5 = boot_table(est[q_names], sm, alpha),
               M6 = boot_table(est[p_names], star, alpha, efron = p_names))
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

# ---- Run --------------------------------------------------------------------------------------
cl <- makeCluster(4)                     # one worker per truth
clusterSetRNGStream(cl, 2026)            # independent random-number streams
clusterExport(cl, ls())                  # the functions and settings above, and those of 20_boot_functions.R
tab30 <- do.call(rbind, parLapply(cl, names(truths), run_truth))
clusterSetRNGStream(cl, 2028)            # Type 1 results do not depend on the settings of the Type 2 part
a1 <- parLapply(cl, lengths(splitIndices(n_sim, 4)), function(r) replicate(r, mc_type1_one()))
stopCluster(cl)
tab30 <- tab30[!is.na(tab30$coverage), ]
tab30$mc_se <- sqrt(tab30$coverage * (1 - tab30$coverage) / n_sim)   # Monte Carlo standard error of the coverage
tab33 <- evaluate(array(unlist(a1), c(dim(a1[[1]])[1:3], n_sim), dimnames(a1[[1]])), truth1)
tab33 <- tab33[!is.na(tab33$coverage), ]
tab33$mc_se <- sqrt(tab33$coverage * (1 - tab33$coverage) / n_sim)

# The report table: coverage (%) at our sample size, n = 209.
tab32 <- with(subset(tab30, n == 209), tapply(round(100 * coverage, 1), list(paste(target, method), factor(truth, names(truths))), sum))
print(tab32)
cat("\nType 1, coverage when the model of 10_type1.R is the truth:\n")
print(data.frame(tab33[1:2], round(tab33[c("true_value", "coverage", "truth_above")], 3)), row.names = FALSE)

# ---- Type 2 arrivals ------------------------------------------------------------------------
# The mean daily count from 21 days, when calls follow a gamma renewal process as regular as
# Type 2 (shape 7) or a Poisson process, both 10 calls a day on the working-hour clock.
day_counts <- function(k, days = 21, rate = 10) tabulate(ceiling(cumsum(rgamma(4 * rate * days, k, k * rate / 9)) / 9), days)
mc_counts_one <- function(Nd) {
  d <- length(Nd); lam <- mean(Nd); Nb <- boot_iid(Nd, B); lb <- colMeans(Nb); lp <- colMeans(boot_param(rpois, d, B, lam))
  cbind(lam, rbind(c(ci_percentile(lam, lb, alpha)), c(ci_percentile_t(lam, sd(Nd) / sqrt(d), lb, col_sd(Nb) / sqrt(d), alpha)),
                   c(ci_percentile_t(lam, sqrt(lam / d), lp, sqrt(lp / d), alpha))))
}
set.seed(2027)
arrival_truths <- c(`gamma renewal, shape 7` = 7, Poisson = 1)
tab31 <- do.call(rbind, lapply(names(arrival_truths), function(name) {
  a <- replicate(n_sim, mc_counts_one(day_counts(arrival_truths[[name]])))      # method x (estimate, lower, upper) x n_sim
  data.frame(truth = name, method = c("resample days, percentile", "resample days, percentile-t (as 11)",
                                      "Poisson model, percentile-t (as 10)"),
             bias = rowMeans(a[, 1, ]) - 10, rmse = sqrt(rowMeans((a[, 1, ] - 10)^2)),
             coverage = rowMeans(a[, 2, ] <= 10 & a[, 3, ] >= 10), length = rowMeans(a[, 3, ] - a[, 2, ]))
}))
tab31$mc_se <- sqrt(tab31$coverage * (1 - tab31$coverage) / n_sim)
print(tab31)

# ---- Output ---------------------------------------------------------------------------------
write.csv(tab30, "output/tab30_mc_durations.csv", row.names = FALSE)
write.csv(tab31, "output/tab31_mc_arrivals.csv", row.names = FALSE)
write.csv(tab32, "output/tab32_mc_coverage_n209.csv")
write.csv(tab33, "output/tab33_mc_type1.csv", row.names = FALSE)
saveRDS(list(durations = tab30, arrivals = tab31, coverage_n209 = tab32, type1 = tab33, truth_values = truth_values,
             truth_type1 = truth1, n_sim = n_sim, B = B), "output/mc_study.rds")

# Figure 30: coverage against n per method (matplot, as in BootstrapCode.R), one panel per
# target and truth; dashed: 95%, dotted: two Monte Carlo standard errors either side.
pdf("output/fig30_mc_coverage.pdf", width = 10, height = 11)
par(mfrow = c(length(targets), length(truths)), mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.6, 0))
for (tg in targets) for (tn in names(truths)) {
  d <- subset(tab30, target == tg & truth == tn)
  matplot(n_grid, tapply(d$coverage, list(d$n, d$method), sum), type = "b", lty = 1, pch = 16, col = seq_along(methods),
          ylim = c(0, 1), log = "x", xlab = "n", ylab = "Coverage", main = paste(tg, "|", tn), cex.main = 0.9)
  abline(h = 0.95 + c(-2, 0, 2) * sqrt(0.95 * 0.05 / n_sim), lty = c(3, 2, 3))
  if (tg == targets[1] && tn == names(truths)[1])
    legend("bottomleft", paste(names(methods), methods), col = seq_along(methods), lty = 1, pch = 16, bty = "n", cex = 0.7)
}
invisible(dev.off())
