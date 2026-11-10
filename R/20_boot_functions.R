# 20_boot_functions.R
# Bootstrap tools shared by 10_type1.R, 11_type2.R and 30_mc_study.R. Functions only:
# sourcing this file runs nothing. Bootstrap samples are matrices with one sample per
# column, so statistics are computed for all B samples at once. Notation as in Boot.R:
# Q* is the bootstrap quantity and c*_a its (1 - a)-quantile.

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

# Equal-tailed cut-offs c*_{alpha/2} and c*_{1 - alpha/2}, one row per statistic. Quantile
# type 6 takes order statistic p (B + 1): exact when alpha / 2 * (B + 1) is an integer (B = 999
# gives 25 and 975, the lecture's rule); with B = 499 it interpolates between two neighbours.
cuts <- function(Q, alpha) t(apply(rbind(Q), 1, quantile, probs = c(1 - alpha / 2, alpha / 2), type = 6, names = FALSE))
lower_upper <- function(M) { colnames(M) <- c("lower", "upper"); M }

# Intervals, one row per statistic: `est` has one entry per row of `star` (the bootstrap
# estimates). Percentile: Q* = theta* - theta, centred on `center`, the value of the
# statistic in the bootstrap world (the estimate itself, except for the smoothed bootstrap).
# Percentile-t: Q* = (theta* - theta) / se*. Efron: the quantiles of theta* themselves.
ci_percentile   <- function(est, star, alpha = 0.05, center = est) lower_upper(est - cuts(star - center, alpha))
ci_percentile_t <- function(est, se, star, se_star, alpha = 0.05) {
  Q <- (star - est) / se_star
  Q[is.nan(Q)] <- 0                      # 0/0: a resample without spread that equals the estimate
  lower_upper(est - cuts(Q, alpha) * se)
}
ci_efron        <- function(star, alpha = 0.05) lower_upper(cuts(star, alpha)[, 2:1, drop = FALSE])

# Table of estimates and intervals. `th` holds the estimates plus "se.<name>" for every
# quantity that gets a percentile-t interval; `star` has the same rows, one column per
# bootstrap sample. Quantities in `efron` (skewed tail probabilities) get Efron's interval,
# all others the percentile interval of the plan.
boot_table <- function(th, star, alpha = 0.05, efron = character()) {
  qs <- names(th)[!startsWith(names(th), "se.")]
  type <- ifelse(paste0("se.", qs) %in% names(th), "percentile-t", ifelse(qs %in% efron, "Efron percentile", "percentile"))
  ci <- t(sapply(seq_along(qs), function(i) {
    q <- qs[i]; se <- paste0("se.", q)
    switch(type[i],
           "percentile-t" = ci_percentile_t(th[[q]], th[[se]], star[q, ], star[se, ], alpha),
           "percentile"   = ci_percentile(th[[q]], star[q, ], alpha),
           ci_efron(star[q, ], alpha))
  }))
  data.frame(estimate = th[qs], lower = ci[, 1], upper = ci[, 2], interval = type, row.names = qs)
}
