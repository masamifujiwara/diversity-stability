# ──────────────────────────────────────────────────────────────
#  Script:  f3_portfolio.R
#  Purpose: Load latest CSVs and analyze the portfolio effect
#           using log10 transform only (time trend + diversity models)
#
#  Author:  Masami Fujiwara with assistance of ChatGPT 5.0 in debugging 
#     Most of the annotations were added by ChatGPT for readability. 
#  Date:    2026-08-24   
# ──────────────────────────────────────────────────────────────

rm(list = ls())
P <- 1  # show figures? 1=YES, 0=NO

# ──────────────────────────────────────────────────────────────
# 0) Dependencies & working directory
# ──────────────────────────────────────────────────────────────
library(tidyverse)
library(mgcv)
library(lme4)
library(lmerTest)
library(nlme)

# Determine script directory so relative paths work anywhere
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
## This script assumes that all of the scv files are saved subdirectory names "results"
results_dir <- file.path(script_dir, "results")


# ──────────────────────────────────────────────────────────────
# 1.  Identify the most‐recent date tag among saved CSVs
# ──────────────────────────────────────────────────────────────
base_names <- c(
  "population_invariability",
  "community_invariability",
  "portfolio_ps",
  "synchrony_ps",
  "L5_alpha_diversity_observed",
  "L5_gamma_diversity_observed",
  "L5_beta_diversity_observed",
  "beta_bray_pairs_monthly",
  "beta_bray_bay_period_season"
)

file_regex <- paste0(
  "^(", paste(base_names, collapse = "|"), ")_\\d{4}-\\d{2}-\\d{2}\\.csv$"
)

all_csv   <- list.files(results_dir, pattern = file_regex, full.names = TRUE)

files_df <- tibble(
  path  = all_csv,
  fname = basename(all_csv),
  base  = str_replace(basename(all_csv), "_\\d{4}-\\d{2}-\\d{2}\\.csv$", ""),
  date  = as.Date(str_extract(basename(all_csv), "\\d{4}-\\d{2}-\\d{2}"))
)

# Keep only bases we care about (in case directory has extras)
files_df <- files_df |> dplyr::filter(base %in% base_names)

# One latest row per base (break ties deterministically by fname)
latest_per_base <- files_df |>
  dplyr::arrange(base, dplyr::desc(date), dplyr::desc(fname)) |>
  dplyr::group_by(base) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()

# Warn if any requested base has no file
missing_bases <- setdiff(base_names, unique(latest_per_base$base))
if (length(missing_bases) > 0) {
  warning("No matching CSV found for: ",
          paste(missing_bases, collapse = ", "))
}


# ──────────────────────────────────────────────────────────────
# 2. Read each latest file into a named list
# ──────────────────────────────────────────────────────────────
latest_data <- latest_per_base |>
  dplyr::mutate(data = purrr::map(path, ~ readr::read_csv(.x, show_col_types = FALSE))) |>
  dplyr::select(base, data) |>
  tibble::deframe()

# ──────────────────────────────────────────────────────────────
# 3.  Attach to global environment (or rename as you like)
# ──────────────────────────────────────────────────────────────

# Note 1:
# Make sure that all of the results saved in the csv files
# are saved in a sub-directory named "results"

# Note 2:
# This will create the following objects:
#   "alpha_div", "gamma_div", "beta_div", "beta_bray",
#   "synchrony",  "portfolio",  "pop_inv",  "com_inv"

list2env(latest_data, envir = .GlobalEnv)

# Optionally, rename your diversity tables for clarity:
alpha_div <- L5_alpha_diversity_observed
gamma_div <- L5_gamma_diversity_observed
beta_div  <- L5_beta_diversity_observed
beta_bray <- beta_bray_bay_period_season
beta_bray_montn <- beta_bray_pairs_monthly
synchrony <- synchrony_ps
portfolio <- portfolio_ps
pop_inv <- population_invariability
com_inv <- community_invariability

# Specify the names of the objects you want to keep
keep_objs <- c(
  "alpha_div", "gamma_div", "beta_div", "beta_bray",
  "synchrony",  "portfolio",  "pop_inv",  "com_inv", "P", "script_dir", "beta_bray_montn"
)

# Remove all other objects from the environment
rm(list = setdiff(ls(), keep_objs))
# ──────────────────────────────────────────────────────────────
# 2) Prepare portfolio data (log10 transform + mid-year)
#    Periods are 3-year bins: P1→1992–1994 (mid≈1993), etc.
#    mid-year = 1991 + 3*P
# ──────────────────────────────────────────────────────────────
port_clean <- portfolio %>%
  mutate(
    period_num = as.numeric(sub("^P", "", period)),
    period_mid = 1991 + 3 * period_num,
    season     = factor(season),
    major_area = factor(major_area),
    port_t     = ifelse(is.finite(portfolio) & portfolio > 0, log10(portfolio), NA_real_)
  ) %>%
  filter(is.finite(port_t), is.finite(period_mid), !is.na(season), !is.na(major_area))

# ──────────────────────────────────────────────────────────────
# 3) Time trend in log10(Portfolio effect)
# ──────────────────────────────────────────────────────────────

# (a) GAM with smooth time + season + bay
gam_port <- gam(
  port_t ~ s(period_mid, k = 6) + season + major_area,
  data   = port_clean, method = "REML"
)
print(summary(gam_port))

# (b) Simple linear trend (centered time), controlling for season + bay
dat_lin <- port_clean %>% mutate(t_c = as.numeric(scale(period_mid, center = TRUE, scale = FALSE)))

m_lin <- gls(
  port_t ~ t_c + season + major_area,
  data = dat_lin, method = "REML"
)
cat("\n[Linear trend on log10(portfolio)]\n")
print(summary(m_lin)$tTable["t_c", , drop = FALSE])  # estimate, SE, t, p

# ── Figures (time) ────────────────────────────────────────────

fig_port_time <- ggplot(port_clean, aes(x = period_mid, y = port_t)) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    method  = "gam", formula = y ~ s(x, k = 6),
    se = TRUE, color = "steelblue", fill = "lightblue", linewidth = 1
  ) +
  labs(
    x = "Mid-period year", y = expression(log[10](Portfolio~effect)),
    title = "Temporal trend in Portfolio Effect (log10 scale)",
    subtitle = "GAM smooth (k = 6) with 95% CI"
  ) +
  theme_minimal()

results_dir <- file.path(script_dir, "results")

ggsave(filename = file.path(results_dir, "S3_fig_PORTFOLIO_time.tif"),
       fig_port_time , width = 16, height = 10, units = "cm",
       dpi = 600, device = "tiff", compression = "lzw")

print(fig_port_time)


# ──────────────────────────────────────────────────────────────
# 4) Diversity → log10(Portfolio effect) mixed models
# ──────────────────────────────────────────────────────────────
df_port_div <- portfolio %>%
  rename(portfolio_effect = portfolio) %>%
  left_join(alpha_div, by = c("major_area", "period", "season")) %>%
  mutate(
    major_area = factor(major_area),
    season     = factor(season),
    port_t     = ifelse(is.finite(portfolio_effect) & portfolio_effect > 0,
                        log10(portfolio_effect), NA_real_)
  ) %>%
  drop_na(port_t, richness, shannon)

# Per-metric LMMs
m_rich_t <- lmer(port_t ~ richness + (1 | major_area) + (1 | season),
                 data = df_port_div)
m_shan_t <- lmer(port_t ~ shannon + (1 | major_area) + (1 | season),
                 data = df_port_div)

# Combined LMM
m_both_t <- lmer(port_t ~ richness + shannon +
                   (1 | major_area) + (1 | season),
                 data = df_port_div)

cat("\n[Richness → log10(Portfolio)]\n");  print(summary(m_rich_t))
cat("\n[Shannon  → log10(Portfolio)]\n");  print(summary(m_shan_t))
cat("\n[Both     → log10(Portfolio)]\n");  print(summary(m_both_t))

# ── Figures (diversity → log10 PE) ────────────────────────────
if (P == 1) {
  fig_rich_t <- ggplot(df_port_div, aes(x = richness, y = port_t)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(
      x = "Richness (q = 0)",
      y = expression(log[10](Portfolio~effect)),
    #  title = "Log10(Portfolio) vs. Richness"
    ) +
    theme_minimal()
  
  fig_shan_t <- ggplot(df_port_div, aes(x = shannon, y = port_t)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(
      x = "Shannon diversity (q = 1)",
      y = expression(log[10](Portfolio~effect)),
     # title = "Log10(Portfolio) vs. Shannon"
    ) +
    theme_minimal()
  
  library(patchwork)
  fig_DS_two_panel <- (fig_rich_t + labs(tag = "a")) / (fig_shan_t + labs(tag = "b"))
  print(fig_DS_two_panel)
}

results_dir <- file.path(script_dir, "results")

ggsave(filename = file.path(results_dir, "fig_2_PORTFOLIO_both.tif"),
   fig_DS_two_panel , width = 8.5, height = 12, units = "cm",
        dpi = 600, device = "tiff", compression = "lzw")

if (requireNamespace("MuMIn", quietly = TRUE)) {
  library(MuMIn)
  r2_rich   <- MuMIn::r.squaredGLMM(m_rich_t)
  r2_shan   <- MuMIn::r.squaredGLMM(m_shan_t)
  r2_both   <- MuMIn::r.squaredGLMM(m_both_t)
  
  print(r2_rich)
  print(r2_shan)
  print(r2_both)
}

# ──────────────────────────────────────────────────────────────
# 5) Minimal reporting helpers 
# ──────────────────────────────────────────────────────────────
log10_mult <- function(beta) 10^beta
cat("\n[Back-transform hints]\n")
cat("  • 1 unit increase in Shannon multiplies PE by ~", round(log10_mult(fixef(m_shan_t)["shannon"]), 2), "x\n", sep = "")
cat("  • +1 species multiplies PE by ~", round(log10_mult(fixef(m_rich_t)["richness"]), 3), "x\n", sep = "")

## Some statistics
# Vector
x <- portfolio$portfolio
x <- x[is.finite(x)]         
n <- length(x)

# Descriptives
med_x  <- median(x)
mean_x <- mean(x)
se_mean <- sd(x) / sqrt(n)   # standard error of the mean

# (Optional) Bootstrap SE and CI for the median
set.seed(2025)
B <- 10000
med_boot <- replicate(B, median(sample(x, n, replace = TRUE)))
se_med   <- sd(med_boot)
ci_med   <- quantile(med_boot, c(0.025, 0.975))

# Wilcoxon signed-rank test vs 1 (one-sided: greater than 1)
wilc <- wilcox.test(x, mu = 1, alternative = "greater",
                    conf.int = TRUE, conf.level = 0.95, exact = FALSE)

# Print a tidy summary
list(
  n = n,
  median = med_x,
  mean = mean_x,
  se_mean = se_mean,
  median_boot_se = se_med,
  median_boot_CI95 = ci_med,
  wilcoxon = list(
    statistic = unname(wilc$statistic),
    p_value = wilc$p.value,
    conf_int = wilc$conf.int,   # CI for the location shift (x - 1)
    estimate = wilc$estimate    # pseudo-median of (x - 1)
  )
)

