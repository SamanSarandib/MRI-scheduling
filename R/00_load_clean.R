# 00_load_clean.R
# Reads the August 2025 call records, checks them, and builds what the rest of the
# project uses: the working-hour clock, the gaps between calls and daily aggregates.

dir.create("output", showWarnings = FALSE)

scans <- read.csv("data/raw/ScanRecords.csv", stringsAsFactors = FALSE)
scans$Date <- as.Date(scans$Date)
scans$type <- ifelse(scans$PatientType == "Type 1", 1L, 2L)
scans$dur_min <- 60 * scans$Duration                    # minutes for reporting, hours kept for computing
scans <- scans[order(scans$Date, scans$Time), ]


# Checks on load: no missing values, calls within 8:00-17:00, positive durations,
# and every weekday between the first and last date present (no weekends, no gaps), 21 in all.
days <- sort(unique(scans$Date))
all_days <- seq(min(days), max(days), by = "day")
weekdays_only <- all_days[format(all_days, "%u") <= "5"]
stopifnot(!anyNA(scans), all(scans$PatientType %in% c("Type 1", "Type 2")),
          all(scans$Time >= 8 & scans$Time <= 17), all(scans$Duration > 0),
          length(days) == length(weekdays_only) && all(days == weekdays_only),
          length(days) == 21)

# Working-hour clock: time outside 8:00-17:00 does not count, so a day is 9 hours long.
scans$day <- match(scans$Date, days)                    # working-day index 1, ..., 21
scans$t <- 9 * (scans$day - 1) + (scans$Time - 8)       # hours since 8:00 on day 1

# Gaps between consecutive calls of the same type, in working hours (also across days).
gaps <- lapply(split(scans$t, scans$type), diff)

# Daily aggregates per type: calls and scan-hours (zero on a day without calls).
daily <- data.frame(date = days, weekday = format(days, "%a"))
for (k in 1:2) {
  s <- scans[scans$type == k, ]
  daily[[paste0("calls", k)]] <- tabulate(s$day, nbins = length(days))
  daily[[paste0("hours", k)]] <- as.vector(tapply(s$Duration, factor(s$day, levels = seq_along(days)),
                                                  sum, default = 0))
}

# Quirk to document: identical (date, time, type) triplets are kept as separate
# patients, since times are only recorded to 0.01 h.
n_dup <- sum(duplicated(scans[c("Date", "Time", "type")]))
cat(sprintf("Loaded %d calls (%d Type 1, %d Type 2) on %d working days; %d repeated (date, time, type) triplets kept.\n",
            nrow(scans), sum(scans$type == 1), sum(scans$type == 2), length(days), n_dup))

saveRDS(list(scans = scans, gaps = gaps, daily = daily), "output/clean.rds")
