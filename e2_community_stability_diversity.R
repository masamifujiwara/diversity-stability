
# ──────────────────────────────────────────────────────────────
#  Script:  e_diversity_stability_analysis.R
#  Purpose: Load the most recent ecological diversity & stability CSVs,
#           perform trend analyses, mixed‐effects modeling,
#           piecewise SEM, and generate visualizations for diversity–stability relationships
#  Author:   Masami Fujiwara
#  Created:  2025-08-07
#
#  Dependencies:
#    • tidyverse (dplyr, ggplot2, tidyr)
#    • readr, stringr, purrr
#    • mgcv, trend
#    • lme4, lmerTest, broom.mixed
#    • MASS, nlme
#    • piecewiseSEM
#    • car
#
#  Inputs:
#    • CSV files in “results/” directory named 
#      `<base_name>_YYYY-MM-DD.csv`
#      where `base_name` ∈ {
#        population_invariability, community_invariability,
#        portfolio_effect, synchrony_ps,
#        L5_alpha_diversity_rarefied, L5_gamma_diversity_rarefied,
#        L5_beta_diversity_rarefied, beta_bray_consecutive_years
#      }
#
#  Outputs (in global environment):
#    • Data frames:
#        – alpha_div, gamma_div, beta_div, beta_bray,
#          synchrony, portfolio, pop_inv, com_inv
#    • Analysis objects:
#        – gamma_4y, gam_gamma, mk.test results
#        – df_CS, m_rr, m_sh
#        – port, gam_port, t.test results
#        – df_port_div, m_rich, m_shan, m_both
#        – df_SEM, sem_out; df_SEM_shan, sem_out_shan
#        – df_beta, cor.test results
#        – aov_IC, TukeyHSD
#        – m_pop, ranef(m_pop)
#    • ggplot figures:
#        – fig_gamma_1, fig_gamma_2, fig1_community, fig2_community,
#          fig_port_time, fig_port_season, fig_port_box_season, fig_port_bay,
#          fig_rich, fig_shan, fig_both, fig_phi_rich, fig_logIC_phi,
#          fig_logIC_rich, fig_phi_shan, fig_logIC_phi_shan, fig_logIC_shan,
#          fig_turnover_stability, fig_turnover_by_season, fig_trends,
#          fig_season, fig_bay, fig_season_summary, fig_interaction,
#          fig_rich_t, fig_shan_t, fig_pop_season, fig_bay_re, fig_spec_re
#
#  Sections:
#    0.  Load dependencies & set working directory
#    1.  Identify most‐recent date tag among saved CSVs
#    2.  Read latest CSVs into a named list & attach to environment
#    3.  Rename objects & clean workspace
#    4.  Regional γ‐diversity trend analysis (GAM & Mann–Kendall test)
#    5.  Community‐level diversity → stability (linear mixed models)
#    6.  Portfolio effect through time (GAM, t‐test)
#    7.  Diversity → synchrony → stability (piecewise SEM)
#    8.  Compositional turnover vs. community stability (Kendall’s τ)
#    9.  Seasonal & bay effects (ANOVA, mixed models)
#   10.  Abundance–variability scaling (Taylor’s law)
#   11.  Population‐level mixed‐effects model
# ──────────────────────────────────────────────────────────────
rm(list = ls())

# Display figures? 1 YES, 0 NO
P <- 1

# ──────────────────────────────────────────────────────────────
# 0.  Load dependencies
# ──────────────────────────────────────────────────────────────
library(tidyverse)
library(patchwork)
library(lme4)
library(lmerTest)

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

# -----------------------------------------------------------------------------
# Data frame: alpha_div
#   • Scope: Local‐scale (α) rarefied diversity
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - rarefied_richness : q = 0 richness, standardized to a common sample size
#       - rarefied_shannon  : q = 1 Shannon diversity on rarefied counts
#   • Use: Controls for unequal sampling when tracking how local species diversity changes
#     across subregions (areas), time periods, and seasons.
#
# Data frame: gamma_div
#   • Scope: Regional‐scale (γ) rarefied diversity within each subregion
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - rarefied_richness : pooled species richness (q = 0)
#       - rarefied_shannon  : pooled Shannon diversity (q = 1)
#   • Use: Assesses how the total species pool within each area shifts over time/season
#     once sampling effort is equalized.
#
# Data frame: beta_div
#   • Scope: Spatial turnover (β) of rarefied diversity among sites
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - beta_richness     : multiplicative β (γ / α) for q = 0
#       - beta_shannon      : multiplicative β for q = 1
#   • Use: Quantifies changes in community heterogeneity within each area over
#     time and season.
#
# Data frame: beta_bray
#   • Scope: Year‐to‐year compositional turnover via Bray–Curtis dissimilarity
#   • Rows indexed by: major_area, season, year2
#   • Key columns:
#       - bray_consecutive  : Bray–Curtis dissimilarity between year2 and year2−1
#   • Use: Tracks temporal stability (or change) in species composition at annual resolution.
#
# Data frame: synchrony
#   • Scope: Community‐wide species synchrony index
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - phi               : Loreau & de Mazancourt’s synchrony index (0–1)
#   • Use: Measures the extent to which species’ abundances fluctuate in concert
#     (high φ) versus compensatory (low φ), driving portfolio effects.
#
# Data frame: portfolio
#   • Scope: Portfolio effect on community stability
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - portfolio_effect  : ratio of community invariability to mean population invariability
#   • Use: Quantifies how diversity and asynchrony combine to buffer community‐level
#     fluctuations relative to individual populations.
#
# Data frame: pop_inv
#   • Scope: Population‐level invariability (stability) of individual species
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - mean_invariability: average of (mean² / variance) across species (and optionally median)
#   • Use: Captures how stable each species’ abundance time series is; provides the
#     baseline for the portfolio calculation.
#
# Data frame: com_inv
#   • Scope: Community‐level invariability of total abundance
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - community_invariability : invariability of the summed abundance time series
#   • Use: Reflects overall community stability, which is compared to pop_inv to yield
#     the portfolio effect.
# -----------------------------------------------------------------------------


# ──────────────────────────────────────────────────────────────
# 2a. Diversity–stability figures (all seasons pooled)
#     Uses your fitted LMMs m_sh (Shannon) and m_rr (Richness)
# ──────────────────────────────────────────────────────────────


df_CS <- com_inv %>%
  left_join(alpha_div, by = c("major_area", "period", "season")) %>%
  dplyr::mutate(log_IC = log10(I_C))

# Richness model
m_rr <- lmer(log_IC ~ richness + (1 | major_area) + (1 | season),
             data = df_CS)
summary(m_rr)

# Shannon model
m_sh <- lmer(log_IC ~ shannon  + (1 | major_area) + (1 | season),
             data = df_CS)
summary(m_sh)


# Pull slope + p from the LMMs for nice subtitles
sh_slope <- unname(lme4::fixef(m_sh)["shannon"])
sh_p     <- summary(m_sh)$coefficients["shannon", "Pr(>|t|)"]

rr_slope <- unname(lme4::fixef(m_rr)["richness"])
rr_p     <- summary(m_rr)$coefficients["richness","Pr(>|t|)"]

fig_DS_shannon_all <- ggplot(df_CS, aes(x = shannon, y = log_IC)) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", na.rm = TRUE) +
  labs(
    x = "Shannon diversity (q = 1)",
    y = expression(log[10](I[C])),
   ) +
  theme_minimal()

fig_DS_richness_all <- ggplot(df_CS, aes(x = richness, y = log_IC)) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", na.rm = TRUE) +
  labs(
    x = "Richness (q = 0)",
    y = expression(log[10](I[C])),
  ) +
  theme_minimal()

# ──────────────────────────────────────────────────────────────
# Combine pooled-season figures into a two-panel figure
# ──────────────────────────────────────────────────────────────

tag_theme <- theme(plot.tag = element_text(face = "bold"))

fig_DS_both <- 
  (
    fig_DS_richness_all + labs(tag = "a") + tag_theme
  ) /
(
  fig_DS_shannon_all + labs(tag = "b") + tag_theme
)  
results_dir <- file.path(script_dir, "results")

ggsave(filename = file.path(results_dir, "fig_DS_both.tif"),
       fig_DS_both, width = 8.5, height = 12, units = "cm",
        dpi = 600, device = "tiff", compression = "lzw")

if (P == 1) print(fig_DS_both)

# ──────────────────────────────────────────────────────────────
# 2b. Season-specific slope tests for Shannon & richness
#     Robust extraction via summary(fit) (lmerTest provides p-values)
# ──────────────────────────────────────────────────────────────

# Helper to fit per-season LMM and extract the slope row safely
.fit_season_lmm <- function(d, season_val, pred_col, pred_label) {
  dd <- d %>%
    dplyr::filter(season == season_val) %>%
    tidyr::drop_na(log_IC, !!rlang::sym(pred_col), major_area)
  
  # not enough data or no variation in predictor
  if (nrow(dd) < 5 || sd(dd[[pred_col]], na.rm = TRUE) == 0) {
    return(tibble(
      season    = season_val,
      metric    = pred_label,
      n         = nrow(dd),
      estimate  = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df        = NA_real_,
      p.value   = NA_real_,
      conf.low  = NA_real_,
      conf.high = NA_real_
    ))
  }
  
  fm  <- as.formula(paste("log_IC ~", pred_col, "+ (1 | major_area)"))
  fit <- lmerTest::lmer(fm, data = dd)
  
  # lmerTest adds df and p-values to coef(summary(fit))
  s  <- summary(fit)
  ct <- as.data.frame(coef(s))  # columns: Estimate, Std. Error, df, t value, Pr(>|t|)
  if (!pred_col %in% rownames(ct)) {
    return(tibble(
      season    = season_val,
      metric    = pred_label,
      n         = nrow(dd),
      estimate  = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df        = NA_real_,
      p.value   = NA_real_,
      conf.low  = NA_real_,
      conf.high = NA_real_
    ))
  }
  
  est <- ct[pred_col, "Estimate"]
  se  <- ct[pred_col, "Std. Error"]
  t   <- ct[pred_col, "t value"]
  dfc <- if ("df" %in% colnames(ct)) ct[pred_col, "df"] else NA_real_
  p   <- if ("Pr(>|t|)" %in% colnames(ct)) ct[pred_col, "Pr(>|t|)"] else NA_real_
  crit <- if (is.finite(dfc)) qt(0.975, df = dfc) else 1.96
  
  tibble(
    season    = season_val,
    metric    = pred_label,
    n         = nrow(dd),
    estimate  = est,
    std.error = se,
    statistic = t,
    df        = dfc,
    p.value   = p,
    conf.low  = est - crit * se,
    conf.high = est + crit * se
  )
}

# Seasons to loop and two metrics
seasons_vec <- sort(unique(df_CS$season))
metrics <- tibble(
  pred_col     = c("shannon", "richness"),
  metric_label = c("Shannon (q = 1)",  "Richness (q = 0)")
)

# Run all 8 tests (4 seasons × 2 metrics), add BH-adjusted q-values
ds_slope_tests <- purrr::map_dfr(
  seasons_vec,
  function(sea) purrr::map_dfr(
    seq_len(nrow(metrics)),
    ~ .fit_season_lmm(df_CS, sea,
                      metrics$pred_col[.x],
                      metrics$metric_label[.x])
  )
) %>%
  dplyr::mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(metric, season)

if (P == 1) {
  message("Season-specific diversity–stability slope tests (LMM):")
  print(ds_slope_tests, n = nrow(ds_slope_tests))
}
