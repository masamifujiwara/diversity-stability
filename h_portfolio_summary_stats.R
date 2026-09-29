# ──────────────────────────────────────────────────────────────
#  Script:  h_portfolio_summary_stats.R
#  Purpose: Compute descriptive statistics for the Portfolio effect
#           (mean, median, SE, bootstrap CI, Wilcoxon test vs. 1)
#           This reproduces the statistics block at the end of
#           e3_portfolio.R as a standalone, reportable summary.
#  Author:  Masami Fujiwara
#  Date:    2026-08-21
# ──────────────────────────────────────────────────────────────

rm(list = ls())

# ──────────────────────────────────────────────────────────────
# 0) Dependencies & working directory
# ──────────────────────────────────────────────────────────────
library(tidyverse)

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    !is.null(rstudioapi::getActiveDocumentContext()$path)) {
  
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  
} else if (!is.null(sys.frames()[[1]]$ofile)) {
  
  script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile))
  
} else {
  stop("Cannot determine script directory. Are you running this interactively?")
}

setwd(script_dir)

## This script assumes portfolio_ps_*.csv is saved in a subdirectory named "results",
## matching the convention used in e3_portfolio.R. Adjust results_dir if needed.
results_dir <- file.path(script_dir, "results")

# ──────────────────────────────────────────────────────────────
# 1) Identify and read the most-recent portfolio_ps CSV
# ──────────────────────────────────────────────────────────────
all_csv <- list.files(
  results_dir,
  pattern = "^portfolio_ps_\\d{4}-\\d{2}-\\d{2}\\.csv$",
  full.names = TRUE
)

if (length(all_csv) == 0) {
  stop("No portfolio_ps_*.csv file found in: ", results_dir)
}

files_df <- tibble(
  path = all_csv,
  fname = basename(all_csv),
  date = as.Date(str_extract(basename(all_csv), "\\d{4}-\\d{2}-\\d{2}"))
) %>%
  arrange(desc(date), desc(fname))

latest_path <- files_df$path[1]
cat("Using file:", basename(latest_path), "\n\n")

portfolio <- readr::read_csv(latest_path, show_col_types = FALSE)

# ──────────────────────────────────────────────────────────────
# 2) Descriptive statistics on the raw (untransformed) Portfolio
#    variable, matching the statistics block in e3_portfolio.R
# ──────────────────────────────────────────────────────────────
x <- portfolio$portfolio
x <- x[is.finite(x)]        # drop NA/NaN/Inf
n <- length(x)

med_x   <- median(x)
mean_x  <- mean(x)
se_mean <- sd(x) / sqrt(n)  # standard error of the mean

# ──────────────────────────────────────────────────────────────
# 3) Bootstrap SE and 95% CI for the median (reproducible seed)
# ──────────────────────────────────────────────────────────────
set.seed(2025)
B <- 10000
med_boot <- replicate(B, median(sample(x, n, replace = TRUE)))
se_med   <- sd(med_boot)
ci_med   <- quantile(med_boot, c(0.025, 0.975))

# ──────────────────────────────────────────────────────────────
# 4) Wilcoxon signed-rank test vs. 1 (one-sided: greater than 1)
# ──────────────────────────────────────────────────────────────
wilc <- wilcox.test(
  x, mu = 1, alternative = "greater",
  conf.int = TRUE, conf.level = 0.95, exact = FALSE
)

# ──────────────────────────────────────────────────────────────
# 5) Report results
# ──────────────────────────────────────────────────────────────
summary_out <- list(
  n                 = n,
  mean              = mean_x,
  se_mean           = se_mean,
  median            = med_x,
  median_boot_se    = se_med,
  median_boot_CI95  = ci_med,
  range             = range(x),
  wilcoxon = list(
    statistic = unname(wilc$statistic),
    p_value   = wilc$p.value,
    conf_int  = wilc$conf.int,   # CI for the location shift (x - 1)
    estimate  = wilc$estimate    # pseudo-median of (x - 1)
  )
)

cat("──────────────────────────────────────────────\n")
cat(" Portfolio effect: descriptive summary (n =", n, ")\n")
cat("──────────────────────────────────────────────\n")
cat(sprintf("  Mean            : %.4f (SE = %.4f)\n", mean_x, se_mean))
cat(sprintf("  Median          : %.4f\n", med_x))
cat(sprintf("  Median boot SE  : %.4f\n", se_med))
cat(sprintf("  Median 95%% CI   : (%.4f, %.4f)\n", ci_med[1], ci_med[2]))
cat(sprintf("  Range           : (%.4f, %.4f)\n", summary_out$range[1], summary_out$range[2]))
cat(sprintf("  Wilcoxon W      : %.2f\n", summary_out$wilcoxon$statistic))
cat(sprintf("  Wilcoxon p-value: %.3g\n", summary_out$wilcoxon$p_value))
cat("──────────────────────────────────────────────\n\n")

print(summary_out)

# ──────────────────────────────────────────────────────────────
# 6) Save summary to CSV for easy inclusion in the manuscript
# ──────────────────────────────────────────────────────────────
out_row <- tibble(
  n                = n,
  mean             = mean_x,
  se_mean          = se_mean,
  median           = med_x,
  median_boot_se   = se_med,
  median_ci_lower  = ci_med[1],
  median_ci_upper  = ci_med[2],
  range_min        = summary_out$range[1],
  range_max        = summary_out$range[2],
  wilcoxon_W       = summary_out$wilcoxon$statistic,
  wilcoxon_p       = summary_out$wilcoxon$p_value
)

out_path <- file.path(results_dir, paste0("portfolio_summary_stats_", Sys.Date(), ".csv"))
readr::write_csv(out_row, out_path)
cat("Summary saved to:", out_path, "\n")