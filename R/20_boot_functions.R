# 20_boot_functions.R
# Bootstrap tools shared by 10_type1.R, 11_type2.R and 30_mc_study.R, the candidate models for
# Type 2 (11 and 30) and p_day_over(), the chance that a day's scans take more than 9 hours
# (10, 11 and 30). Functions only: sourcing this file runs nothing. Bootstrap samples are
# matrices with one sample per column, so statistics are computed for all B samples at once.
# Notation as in Boot.R: Q* is the bootstrap quantity and c*_a its (1 - a)-quantile.

# Resampling: iid (Boot.R with all B index draws at once), smoothed (D* = D_J + h eps,
# Silverman's bandwidth) and parametric (any r-function with its parameters).
boot_iid    <- function(x, B) matrix(x[sample.int(length(x), length(x) * B, replace = TRUE)], length(x))
boot_smooth <- function(x, B, h = bw.nrd0(x)) boot_iid(x, B) + h * rnorm(length(x) * B)
boot_param  <- function(rdist, n, B, ...) matrix(rdist(n * B, ...), n)

# Column statistics: sd, the EDF quantile F_n^{-1}(p) = X_(ceiling(np)), and the standard
# error of a share p with half a count added, so a share of 0 still has a positive one.
col_sd       <- function(X) { X <- as.matrix(X); sqrt(colSums(sweep(X, 2, colMeans(X))^2) / (nrow(X) - 1)) }
col_quantile <- function(X, p) apply(as.matrix(X), 2, function(x) { k <- ceiling(length(x) * p); sort(x, partial = k)[k] })
se_prop      <- function(p, n) { pt <- (n * p + 0.5) / (n + 1); sqrt(pt * (1 - pt) / n) }

# p-quantile of the smoothed distribution itself (the EDF convolved with N(0, h^2)): the
# value that the smoothed bootstrap is centred on.
kde_quantile <- function(x, p, h = bw.nrd0(x)) uniroot(function(q) mean(pnorm((q - x) / h)) - p, range(x) + c(-5, 5) * h)$root

# Gamma maximum likelihood for each column (positive data): the shape solves
# log(k) - digamma(k) = log(mean) - mean(log), by Newton from Minka's closed-form start.
gamma_ml <- function(X) {
  X <- as.matrix(X); m <- colMeans(X); s <- log(m) - colMeans(log(X))
  k <- (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
  for (i in 1:6) k <- k - (log(k) - digamma(k) - s) / (1 / k - trigamma(k))
  rbind(shape = k, rate = k / m)
}

# Candidate models for the durations. Each has an ML fit for every column of a matrix (one
# parameter column per sample), a random generator, quantile and survival functions, and its
# mean and sd. Exponential and lognormal ML are closed form, the gamma uses gamma_ml above
# and the Weibull weibull_ml below. Used by 11_type2.R for the durations and the gaps between
# calls, and by 30_mc_study.R. normal_model (no mean and sd functions) is the Type 1
# assumption, which 11 includes only in the fit test.
# Weibull ML: the shape k solves 1/k + mean(log x) - sum(x^k log x) / sum(x^k) = 0 (Newton from
# the usual CV-based start), then scale = mean(x^k)^(1/k); it matches MASS::fitdistr. Newton
# stops once every column's step is below `tol` times its shape; a step that would make the
# shape negative is replaced by halving it, and a column that has not converged after
# `max_iter` steps stops the script with an error rather than giving a wrong fit.
weibull_ml <- function(X, tol = 1e-10, max_iter = 50) {
  X <- as.matrix(X); L <- log(X); mL <- colMeans(L)
  k <- (col_sd(X) / colMeans(X))^-1.086
  for (i in seq_len(max_iter)) {
    Xk <- X^rep(k, each = nrow(X)); S0 <- colSums(Xk); S1 <- colSums(Xk * L); S2 <- colSums(Xk * L^2)
    step <- (1 / k + mL - S1 / S0) / (-1 / k^2 - (S2 / S0 - (S1 / S0)^2))
    k <- ifelse(k - step > 0, k - step, k / 2)
    if (isTRUE(all(abs(step) < tol * k))) break
  }
  ok <- is.finite(k) & is.finite(step) & abs(step) < tol * k
  if (!all(ok)) stop(sprintf("weibull_ml: no convergence after %d Newton steps for %d of %d samples",
                             max_iter, sum(!ok), length(k)))
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
normal_model <- list(fit = function(X) rbind(mean = colMeans(as.matrix(X)), sd = col_sd(X)),
                     r   = function(m, p) rnorm(m, p[1, ], p[2, ]),
                     d   = function(x, p) dnorm(x, p[1, ], p[2, ]),
                     q   = function(pr, p) qnorm(pr, p[1, ], p[2, ]),
                     sf  = function(a, p) pnorm(a, p[1, ], p[2, ], lower.tail = FALSE))

# Equal-tailed cut-offs c*_{alpha/2} and c*_{1 - alpha/2}, one row per statistic. Quantile
# type 6 takes order statistic p (B + 1): exact when alpha / 2 * (B + 1) is an integer (B = 999
# gives 25 and 975, the lecture's rule); with B = 499 it interpolates between two neighbours.
cuts <- function(Q, alpha) t(apply(rbind(Q), 1, quantile, probs = c(1 - alpha / 2, alpha / 2), type = 6, names = FALSE))
lower_upper <- function(M) { colnames(M) <- c("lower", "upper"); M }

# Intervals, one row per statistic: `est` has one entry per row of `star` (the bootstrap
# estimates). Percentile: Q* = theta* - theta, centred on `center`, the value of the
# statistic in the bootstrap world (the estimate itself, except for the smoothed bootstrap).
# Percentile-t: Q* = (theta* - theta) / se*. Efron: the quantiles of theta* themselves.
# Bias-corrected (Efron's BC): the quantiles of theta* at Phi(2 z0 -/+ z_{1 - alpha/2}), where
# z0 = Phi^{-1}(share of theta* below the estimate) corrects for a bootstrap distribution
# that is not centred on the estimate (z0 = 0 gives Efron's interval).
ci_percentile   <- function(est, star, alpha = 0.05, center = est) lower_upper(est - cuts(star - center, alpha))
ci_percentile_t <- function(est, se, star, se_star, alpha = 0.05) {
  Q <- (star - est) / se_star
  Q[is.nan(Q)] <- 0                      # 0/0: a resample without spread that equals the estimate
  lower_upper(est - cuts(Q, alpha) * se)
}
ci_efron        <- function(star, alpha = 0.05) lower_upper(cuts(star, alpha)[, 2:1, drop = FALSE])
ci_bc           <- function(est, star, alpha = 0.05) {
  z0 <- qnorm(mean(star < est))
  lower_upper(rbind(quantile(star, pnorm(2 * z0 + qnorm(c(alpha / 2, 1 - alpha / 2))), type = 6, names = FALSE)))
}

# Table of estimates and intervals. `th` holds the estimates plus "se.<name>" for every
# quantity that gets a percentile-t interval; `star` has the same rows, one column per
# bootstrap sample. Quantities in `efron` (skewed tail probabilities) get Efron's interval,
# those in `bc` the bias-corrected one, all others the percentile interval of the plan.
boot_table <- function(th, star, alpha = 0.05, efron = character(), bc = character()) {
  qs <- names(th)[!startsWith(names(th), "se.")]
  type <- ifelse(paste0("se.", qs) %in% names(th), "percentile-t",
                 ifelse(qs %in% efron, "Efron percentile", ifelse(qs %in% bc, "bias-corrected percentile", "percentile")))
  ci <- t(sapply(seq_along(qs), function(i) {
    q <- qs[i]; se <- paste0("se.", q)
    switch(type[i],
           "percentile-t" = ci_percentile_t(th[[q]], th[[se]], star[q, ], star[se, ], alpha),
           "percentile"   = ci_percentile(th[[q]], star[q, ], alpha),
           "bias-corrected percentile" = ci_bc(th[[q]], star[q, ], alpha),
           ci_efron(star[q, ], alpha))
  }))
  data.frame(estimate = th[qs], lower = ci[, 1], upper = ci[, 2], interval = type, row.names = qs)
}

# P(one day's scans take more than `limit` minutes) when the number of calls that day is drawn
# from the observed daily counts N and every scan from the observed durations D (minutes), all
# independently: p_work9 of 10_type1.R with the observed distributions in place of the Poisson
# and the normal. One value per column of D and N (one sample per column, as above). Exact for
# those distributions up to a grid: the durations are rounded to `step` minutes and the grid
# point at the limit counts half (a continuity correction); against a 0.002-minute grid the
# answer moves by less than 0.1%. The sum of k scans has the k-fold convolution of the duration
# distribution, whose Fourier transform is f^k, and P(sum > limit) is that convolution summed
# with the weights of the grid points (1 above the limit, 0 below), which Parseval's identity
# gives as sum(f^k * Conj(fft(weights))) / L, so no inverse transform is needed.
p_day_over <- function(D, N, limit = 540, step = 0.25) {
  D <- as.matrix(D); N <- as.matrix(N); i <- round(D / step)     # grid index of every scan, 0 = 0 minutes
  L <- nextn(max(N) * max(i) + 1)                                # FFT length: room for the longest day
  cw <- Conj(fft(pmin(pmax(seq_len(L) - 1 - limit / step + 0.5, 0), 1)))   # transform of the weights
  sapply(seq_len(ncol(D)), function(b) {
    f <- fft(tabulate(i[, b] + 1, L) / nrow(D))                  # transform of the duration distribution
    nk <- tabulate(N[, b] + 1); k <- which(nk > 0) - 1           # the numbers of calls that occur
    p <- numeric(length(k)); fk <- f^k[1]
    for (j in seq_along(k)) {
      if (j > 1) fk <- fk * f^(k[j] - k[j - 1])
      p[j] <- Re(sum(fk * cw)) / L                               # P(sum of k[j] scans > limit)
    }
    sum(nk[k + 1] / nrow(N) * p)
  })
}
