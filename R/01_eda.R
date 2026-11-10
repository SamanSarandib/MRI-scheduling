# 01_eda.R
# Tables and figures for the data section: calls per day, scan durations and the
# arrival pattern of both patient types. Output: 
# output/tab01_data_summary.csv and output/fig01-03_*.pdf.

if (!exists("scans")) source("R/00_load_clean.R")

skew <- function(x) mean((x - mean(x))^3) / mean((x - mean(x))^2)^1.5

# Table: one column per type
summarise_type <- function(k) {
  D <- scans$dur_min[scans$type == k]; N <- daily[[paste0("calls", k)]]
  H <- daily[[paste0("hours", k)]];    G <- gaps[[k]]
  c(calls = length(D), calls_per_day = mean(N), calls_min = min(N), calls_max = max(N),
    var_to_mean = var(N) / mean(N), gap_mean_h = mean(G), gap_cv = sd(G) / mean(G),
    dur_mean = mean(D), dur_sd = sd(D), dur_skew = skew(D), dur_max = max(D),
    hours_per_day = mean(H), hours_max = max(H))
}
tab01 <- sapply(1:2, summarise_type)
colnames(tab01) <- c("Type 1", "Type 2")
print(round(tab01, 2))
write.csv(round(tab01, 3), "output/tab01_data_summary.csv")

pooled <- daily$hours1 + daily$hours2
cat(sprintf("Both types together: %.1f scan-hours a day on average, %.1f at most, against 18 machine-hours.\n",
            mean(pooled), max(pooled)))

cols <- c("grey40", "steelblue4")

# Figure 1: calls per working day.
pdf("output/fig01_calls_per_day.pdf", width = 7, height = 4)
matplot(daily[c("calls1", "calls2")], type = "o", pch = 16, lty = 1, col = cols, ylim = c(0, 25),
        xaxt = "n", xlab = "Working day in August 2025 (day of the month)", ylab = "Calls per day")
axis(1, at = seq_len(nrow(daily)), labels = format(daily$date, "%d"), cex.axis = 0.8, gap.axis = 0)
abline(h = colMeans(daily[c("calls1", "calls2")]), col = cols, lty = 2)
legend("bottomleft", c("Type 1", "Type 2"), col = cols, lty = 1, pch = 16, bty = "n")
invisible(dev.off())

# Figure 2: durations in 5-minute bins, with fitted curves scaled to counts and candidate slots.
pdf("output/fig02_durations.pdf", width = 7, height = 7)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
for (k in 1:2) {
  D <- scans$dur_min[scans$type == k]; m <- mean(D); s <- sd(D)
  h <- hist(D, breaks = seq(0, 90, 5), plot = FALSE)
  plot(h, ylim = c(0, 1.1 * max(h$counts, 5 * length(D) * dnorm(m, m, s))), col = cols[k], border = "white",
       main = paste0("Type ", k, ": ", length(D), " scans"), xlab = "Scan duration (minutes)", ylab = "Scans per bin")
  curve(5 * length(D) * dnorm(x, m, s), add = TRUE, lwd = 2)
  if (k == 2) curve(5 * length(D) * dgamma(x, shape = m^2 / s^2, rate = m / s^2), add = TRUE, lwd = 2, lty = 2)
  abline(v = if (k == 1) c(30, 35) else c(45, 60), lty = 3)
  legend("topright", c("normal fit", if (k == 2) "gamma fit (moments)"), lty = 1:k, lwd = 2, bty = "n")
}
invisible(dev.off())

# Figure 3: gaps between calls and calls by hour.
pdf("output/fig03_arrivals.pdf", width = 8, height = 6)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (k in 1:2) {
  G <- gaps[[k]]
  h <- hist(G, breaks = seq(0, ceiling(max(G) * 4) / 4, 0.25), plot = FALSE)
  plot(h, freq = FALSE, ylim = c(0, 1.1 * max(h$density, 1 / mean(G))), col = cols[k], border = "white",
       main = sprintf("Type %d: gaps between calls (CV %.2f)", k, sd(G) / mean(G)), xlab = "Working hours")
  curve(dexp(x, 1 / mean(G)), add = TRUE, lwd = 2)
}
for (k in 1:2) {
  barplot(table(factor(floor(scans$Time[scans$type == k]), levels = 8:16)) / nrow(daily), col = cols[k],
          border = "white", main = sprintf("Type %d: calls per hour of the day", k), xlab = "Hour", ylab = "Calls per day")
}
invisible(dev.off())
