# 10_type1.R
# Type 1: normal durations and Poisson calls, as the hospital says. Maximum likelihood
# estimates, plug-in hospital quantities, a parametric bootstrap for the intervals and a
# nonparametric bootstrap as a check, plus the diagnostics for the appendix.

# Output: output/tab10_type1.csv, output/type1.rds, output/fig10-11_*.pdf, and
# `type1_inputs` (generators in minutes) for the simulation in Part 2.

if (!exists("scans")) source("R/00_load_clean.R")
source("R/20_boot_functions.R")
set.seed(2025)
B <- 999                                 # alpha / 2 * (B + 1) = 25 Int
alpha <- 0.05
slots <- c(25, 30, 35)                   # candidate Type 1 slot lengths, minutes

D <- scans$dur_min[scans$type == 1]      # 356 durations, minutes
N <- daily$calls1                        # 21 daily counts
G <- gaps[["1"]]                         # 355 gaps between calls, working hours

# P(a day's scan time > 9 h), exact under the model: given k calls it is N(k mu, k sigma^2).
p_work9 <- function(lam, mu, s, k = 1:100) sum(dpois(k, lam) * pnorm(540, k * mu, s * sqrt(k), lower.tail = FALSE))

# Every quantity from one data set (D, N, G), plus the standard errors ("se") of those
# that get a percentile-t interval: model-based for lambda, delta method for the rest.
type1_quantities <- function(D, N, G) {
  mu <- mean(D); s <- sd(D); lam <- mean(N); lam_g <- 9 / mean(G); z <- (slots - mu) / s
  c(lambda = lam, lambda_gaps = lam_g, mu = mu, sigma = s,
    q90 = mu + qnorm(0.90) * s, q95 = mu + qnorm(0.95) * s,
    setNames(pnorm(z, lower.tail = FALSE), paste0("p_over_", slots)),
    hours = lam * mu / 60, p_work9 = p_work9(lam, mu, s),
    se.lambda = sqrt(lam / length(N)), se.lambda_gaps = lam_g / sqrt(length(G)), se.mu = s / sqrt(length(D)),
    setNames(dnorm(z) * sqrt((1 + z^2 / 2) / length(D)), paste0("se.p_over_", slots)),
    se.hours = sqrt(mu^2 * lam / length(N) + lam^2 * s^2 / length(D)) / 60)
}
th <- type1_quantities(D, N, G)

# B fake Augusts, one per column; then every quantity per fake August.
boot <- function(Ds, Ns, Gs) sapply(1:B, function(b) type1_quantities(Ds[, b], Ns[, b], Gs[, b]))
par_star <- boot(boot_param(rnorm, length(D), B, th[["mu"]], th[["sigma"]]),   # primary: parametric
                 boot_param(rpois, length(N), B, th[["lambda"]]),
                 boot_param(rexp, length(G), B, th[["lambda_gaps"]] / 9))
np_star <- boot(boot_iid(D, B), boot_iid(N, B), boot_iid(G, B))               # check: nonparametric

# Equal-tailed 1 - alpha intervals (20_boot_functions.R): percentile-t when a standard error
# exists, the percentile interval of the plan otherwise. P(workload > 9 h) has a skewed
# bootstrap distribution, where that reflected interval under-covers (about 88%), so it gets
# Efron's percentile interval.
par_tab <- boot_table(th, par_star, alpha, efron = "p_work9")
np_tab <- boot_table(th, np_star, alpha, efron = "p_work9")
tab10 <- data.frame(par_tab["estimate"], par_lower = par_tab$lower, par_upper = par_tab$upper,
                    np_lower = np_tab$lower, np_upper = np_tab$upper, par_tab["interval"])
print(format(tab10, digits = 3))
write.csv(tab10, "output/tab10_type1.csv")

# Busy-day range: 90% prediction interval for a day's calls, Poisson(lambda) averaged over
# the bootstrap lambdas so that the estimation error widens it.
k <- 0:60
F_mix <- sapply(k, function(j) mean(ppois(j, par_star["lambda", ])))
busy <- rbind(plug_in = qpois(c(0.05, 0.95), th[["lambda"]]),
              bootstrap = k[c(which(F_mix >= 0.05)[1], which(F_mix >= 0.95)[1])],
              observed_min_max = range(N))
colnames(busy) <- c("5%", "95%")
print(busy)

# Diagnostics: dispersion test of the daily counts (chi-square, 20 df), lag-1 autocorrelations
# in call order, Shapiro-Wilk, and the observed overrun shares next to the normal ones.
X2 <- sum((N - mean(N))^2) / mean(N)
diagnostics <- c(dispersion_X2 = X2, p_over_dispersed = pchisq(X2, length(N) - 1, lower.tail = FALSE),
                 p_under_dispersed = pchisq(X2, length(N) - 1),
                 acf1_durations = acf(D, plot = FALSE)$acf[2], acf1_gaps = acf(G, plot = FALSE)$acf[2],
                 shapiro_p = shapiro.test(D)$p.value,
                 setNames(sapply(slots, function(s) mean(D > s)), paste0("observed_over_", slots)))
print(round(diagnostics, 3))

# Figure 10: durations with the fitted normal and the candidate slots (the report figure).
pdf("output/fig10_type1_durations.pdf", width = 7, height = 4)
h <- hist(D, breaks = seq(0, 50, 2.5), plot = FALSE)
plot(h, ylim = c(0, 1.1 * max(h$counts, 2.5 * length(D) * dnorm(0, 0, th[["sigma"]]))), col = "grey70",
     border = "white", main = "Type 1 scan durations", xlab = "Minutes", ylab = "Scans per bin")
curve(2.5 * length(D) * dnorm(x, th[["mu"]], th[["sigma"]]), add = TRUE, lwd = 2)
abline(v = slots, lty = 3)
legend("topright", sprintf("%d min: %.1f%%", slots, 100 * th[paste0("p_over_", slots)]), bty = "n",
       title = "Fitted share over the slot", cex = 0.8)
invisible(dev.off())

# Figure 11: diagnostics, normal QQ of durations, exponential QQ of gaps, daily counts against Poisson.
pdf("output/fig11_type1_diagnostics.pdf", width = 9, height = 3.5)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
qqnorm(D, main = "Durations: normal QQ", pch = 16, cex = 0.5); qqline(D)
qqplot(qexp(ppoints(length(G)), th[["lambda_gaps"]] / 9), G, main = "Gaps: exponential QQ", pch = 16, cex = 0.5,
       xlab = "Theoretical quantiles (h)", ylab = "Sample quantiles (h)"); abline(0, 1)
plot(table(factor(N, levels = 5:30)) / length(N), ylab = "Share of days", xlab = "Calls per day",
     main = "Daily counts vs Poisson")
points(5:30, dpois(5:30, th[["lambda"]]), pch = 16, col = "steelblue4")
invisible(dev.off())

# Hand-over to Part 2 in the plan's format: minutes, gaps on the working-minute clock (540 a
# day). Values are fixed here, so a later `th` in the workspace cannot change them.
type1_inputs <- local({
  mu <- th[["mu"]]; s <- th[["sigma"]]; lam <- th[["lambda"]]
  list(dur1 = function(n) rnorm(n, mu, s), gap1 = function(n) rexp(n, lam / 540))
})

saveRDS(list(estimates = th, table = tab10, busy_day = busy, diagnostics = diagnostics,
             par_star = par_star, np_star = np_star, inputs = type1_inputs), "output/type1.rds")
