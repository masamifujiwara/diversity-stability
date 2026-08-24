
# ──────────────────────────────────────────────────────────────
#  Script:  e2_diversity_stability_analysis.R
#  Purpose: Load the most recent ecological diversity & stability CSVs,
#           perform trend analyses, mixed‐effects modeling,
#           piecewise SEM, and generate visualizations for diversity–stability relationships
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
#        L5_alpha_diversity, L5_gamma_diversity,
#        L5_beta_diversity, beta_bray_consecutive_years
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
#
#  Author:  Masami Fujiwara with assistance of ChatGPT 5.0 in debugging 
#     Most of the annotations were added by ChatGPT for readability. 
#  Date:    2026-08-24   
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
# 3.  Attach to global environment 
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
  "alpha_div", "gamma_div", "beta_div", "beta_bray","results_dir",
  "synchrony",  "portfolio",  "pop_inv",  "com_inv", "P", "script_dir", "beta_bray_montn"
)

# Remove all other objects from the environment
rm(list = setdiff(ls(), keep_objs))

# -----------------------------------------------------------------------------
# Data frame: alpha_div
#   • Scope: Local‐scale (α) diversity
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - richness : q = 0 richness, standardized to a common sample size
#       - shannon  : q = 1 Shannon diversity 
#   • Use: Controls for unequal sampling when tracking how local species diversity changes
#     across subregions (areas), time periods, and seasons.
#
# Data frame: gamma_div
#   • Scope: Regional‐scale (γ) diversity within each subregion
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - richness : pooled species richness (q = 0)
#       - shannon  : pooled Shannon diversity (q = 1)
#   • Use: Assesses how the total species pool within each area shifts over time/season
#     once sampling effort is equalized.
#
# Data frame: beta_div
#   • Scope: Spatial turnover (β) of diversity among sites
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
# 2a. Diversity–stability (all seasons pooled): fit + display + save
# ──────────────────────────────────────────────────────────────

stamp <- as.character(Sys.Date())

# Map columns robustly
alpha_div_clean <- alpha_div %>%
  select(major_area, period, season, richness, shannon)

com_inv_clean <- com_inv %>%
  dplyr::select(major_area, period, season, I_C)

# Build analysis frame
df_CS <- com_inv_clean %>%
  left_join(alpha_div_clean, by = c("major_area","period","season")) %>%
  mutate(
    season     = as.character(season),   # safer for printing/saving
    log_IC     = log10(I_C)
  ) %>%
  tidyr::drop_na(log_IC, richness, shannon, major_area, season)

# Fit LMMs
m_rr <- lmer(log_IC ~ richness + (1 | major_area) + (1 | season), data = df_CS)
m_sh <- lmer(log_IC ~ shannon  + (1 | major_area) + (1 | season), data = df_CS)

# Extract fixed-effect rows for predictor of interest
.extract_fixed <- function(fit, term_label) {
  s  <- summary(fit)
  ct <- as.data.frame(coef(s))
  out <- tibble::tibble(
    term      = term_label,
    estimate  = ct[term_label, "Estimate"],
    std.error = ct[term_label, "Std. Error"],
    t_value   = ct[term_label, "t value"],
    df        = if ("df" %in% colnames(ct)) ct[term_label, "df"] else NA_real_,
    p_value   = if ("Pr(>|t|)" %in% colnames(ct)) ct[term_label, "Pr(>|t|)"] else NA_real_
  )
  out
}

# Model fit stats, groups, AIC, (optional) R2m/R2c if MuMIn is available
.extract_fit <- function(fit, label) {
  s   <- summary(fit)
  ng  <- lme4::ngrps(fit)
  vc  <- as.data.frame(lme4::VarCorr(fit)) %>% select(grp, sdcor)
  r2m <- r2c <- NA_real_
  if (requireNamespace("MuMIn", quietly = TRUE)) {
    r2    <- MuMIn::r.squaredGLMM(fit)
    r2m   <- as.numeric(r2[1])
    r2c   <- as.numeric(r2[2])
  }
  tibble::tibble(
    model   = label,
    AIC     = AIC(fit),
    n       = nobs(fit),
    groups_major_area = ng[["major_area"]],
    groups_season     = ng[["season"]],
    R2_marginal       = r2m,
    R2_conditional    = r2c,
    resid_sd          = s$sigma
  ) -> fitrow
  
  # also return random effects SDs in a tidy table
  re_tab <- tibble::as_tibble(vc) %>%
    transmute(model = label, effect = grp, sd = sdcor)
  
  list(fit = fitrow, re = re_tab)
}

# Build tidy outputs
fx_rich <- .extract_fixed(m_rr, "richness") %>% mutate(model = "logIC ~ richness")
fx_shan <- .extract_fixed(m_sh, "shannon")  %>% mutate(model = "logIC ~ shannon")

fits_r  <- .extract_fit(m_rr, "logIC ~ richness")
fits_s  <- .extract_fit(m_sh, "logIC ~ shannon")

ds_overall_fixed  <- bind_rows(fx_rich, fx_shan) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  select(model, term, estimate, std.error, t_value, df, p_value, q_value)

ds_overall_fit    <- bind_rows(fits_r$fit, fits_s$fit)
ds_overall_random <- bind_rows(fits_r$re,  fits_s$re)

# DISPLAY (console) when P==1
if (P == 1) {
  cat("\n[Overall LMM — fixed effects]\n"); print(ds_overall_fixed, n = nrow(ds_overall_fixed))
  cat("\n[Overall LMM — fit stats]\n");    print(ds_overall_fit)
  cat("\n[Overall LMM — random effect SDs]\n"); print(ds_overall_random, n = nrow(ds_overall_random))
}

# SAVE
readr::write_csv(ds_overall_fixed,  file.path(results_dir, paste0("ds_overall_fixed_",  stamp, ".csv")))
readr::write_csv(ds_overall_fit,    file.path(results_dir, paste0("ds_overall_fit_",    stamp, ".csv")))
readr::write_csv(ds_overall_random, file.path(results_dir, paste0("ds_overall_random_", stamp, ".csv")))

# ──────────────────────────────────────────────────────────────
# 2b. Season-specific slope tests for Shannon & richness (per-season LMMs)
# ──────────────────────────────────────────────────────────────

.fit_season_lmm <- function(d, season_val, pred_col, pred_label) {
  dd <- d %>%
    dplyr::filter(season == season_val) %>%
    tidyr::drop_na(log_IC, !!rlang::sym(pred_col), major_area)
  
  if (nrow(dd) < 5 || sd(dd[[pred_col]], na.rm = TRUE) == 0) {
    return(tibble::tibble(
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
  
  fm  <- stats::as.formula(paste("log_IC ~", pred_col, "+ (1 | major_area)"))
  fit <- lmerTest::lmer(fm, data = dd)
  s   <- summary(fit)
  ct  <- as.data.frame(coef(s))
  
  est <- ct[pred_col, "Estimate"]
  se  <- ct[pred_col, "Std. Error"]
  t   <- ct[pred_col, "t value"]
  dfc <- if ("df" %in% colnames(ct)) ct[pred_col, "df"] else NA_real_
  p   <- if ("Pr(>|t|)" %in% colnames(ct)) ct[pred_col, "Pr(>|t|)"] else NA_real_
  crit <- if (is.finite(dfc)) qt(0.975, df = dfc) else 1.96
  
  tibble::tibble(
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

seasons_vec <- sort(unique(df_CS$season))
metrics <- tibble::tibble(
  pred_col     = c("shannon", "richness"),
  metric_label = c("Shannon (q = 1)",  "Richness (q = 0)")
)

ds_slope_tests <- purrr::map_dfr(
  seasons_vec,
  function(sea) purrr::map_dfr(
    seq_len(nrow(metrics)),
    ~ .fit_season_lmm(df_CS, sea, metrics$pred_col[.x], metrics$metric_label[.x])
  )
) %>%
  dplyr::mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(metric, season)

# DISPLAY and SAVE
if (P == 1) {
  cat("\n[Season-specific diversity–stability slopes (LMM)]\n")
  print(ds_slope_tests, n = nrow(ds_slope_tests))
}
readr::write_csv(ds_slope_tests, file.path(results_dir, paste0("ds_season_slopes_", stamp, ".csv")))

# ──────────────────────────────────────────────────────────────
# Combine pooled-season figures into a two-panel figure
# ──────────────────────────────────────────────────────────────

# Ensure df_CS exists and has needed columns
stopifnot(exists("df_CS"))
stopifnot(all(c("richness","shannon","log_IC") %in% names(df_CS)))

# Pooled-season scatter + OLS line (visual only; LMMs reported elsewhere)
fig_DS_richness_all <- ggplot(df_CS, aes(x = richness, y = log_IC)) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", na.rm = TRUE) +
  labs(
    x = "Richness (q = 0, observed; no rarefaction)",
    y = expression(log[10](I[C]))
  ) +
  theme_minimal()

fig_DS_shannon_all <- ggplot(df_CS, aes(x = shannon, y = log_IC)) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", na.rm = TRUE) +
  labs(
    x = "Shannon diversity (q = 1, observed; no rarefaction)",
    y = expression(log[10](I[C]))
  ) +
  theme_minimal()


tag_theme <- theme(plot.tag = element_text(face = "bold"))

fig_DS_both <- 
  (
    fig_DS_richness_all + labs(tag = "a") + tag_theme
  ) /
(
  fig_DS_shannon_all + labs(tag = "b") + tag_theme
)  

ggsave(filename = file.path(results_dir, "fig_1_DS_both.tif"),
       fig_DS_both, width = 8.5, height = 12, units = "cm",
        dpi = 600, device = "tiff", compression = "lzw")

if (P == 1) print(fig_DS_both)

# ──────────────────────────────────────────────────────────────
# 3. Diversity vs Community Stability — per-season results
#    Adds slope ± SE, t, p/q, R² (OLS), n; and optional LMM R²m/R²c
# ──────────────────────────────────────────────────────────────

stamp <- as.character(Sys.Date())

# Safety checks
stopifnot(all(c("season","major_area","log_IC","richness","shannon") %in% names(df_CS)))

# ---------- Per-season OLS helpers ----------
.fit_ols <- function(d, y, x) {
  fml <- as.formula(paste(y, "~", x))
  fit <- lm(fml, data = d)
  s   <- summary(fit)
  co  <- s$coefficients
  # slope row is named by predictor x
  slope <- co[x, "Estimate"]
  se    <- co[x, "Std. Error"]
  tval  <- co[x, "t value"]
  pval  <- co[x, "Pr(>|t|)"]
  n     <- nrow(model.frame(fit))
  r2    <- s$r.squared
  r2a   <- s$adj.r.squared
  fstat <- unname(s$fstatistic["value"])
  df1   <- unname(s$fstatistic["numdf"])
  df2   <- unname(s$fstatistic["dendf"])
  tibble::tibble(
    slope = slope, se = se, t = tval, p = pval,
    R2 = r2, R2_adj = r2a, n = n,
    F = fstat, df1 = df1, df2 = df2
  )
}

# Build per-season OLS tables for richness and Shannon
ols_rich <- df_CS %>%
  dplyr::group_by(season) %>%
  dplyr::group_modify(~ .fit_ols(.x, y = "log_IC", x = "richness")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metric = "Richness (q = 0)")

ols_shan <- df_CS %>%
  dplyr::group_by(season) %>%
  dplyr::group_modify(~ .fit_ols(.x, y = "log_IC", x = "shannon")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metric = "Shannon (q = 1)")

# Add BH-adjusted q-values within each metric across the 4 seasons
ds_season_ols <- dplyr::bind_rows(ols_rich, ols_shan) %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(q = p.adjust(p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(metric, season)

# ---------- Optional: per-season LMM with R²m/R²c ----------
# log_IC ~ predictor + (1 | major_area)
.get_lmm_r2 <- function(d, predictor) {
  fml <- as.formula(paste("log_IC ~", predictor, "+ (1 | major_area)"))
  fit <- lmerTest::lmer(fml, data = d)
  r2m <- r2c <- NA_real_
  if (requireNamespace("MuMIn", quietly = TRUE)) {
    r2  <- MuMIn::r.squaredGLMM(fit)
    r2m <- as.numeric(r2[1]); r2c <- as.numeric(r2[2])
  }
  tibble::tibble(R2_marginal = r2m, R2_conditional = r2c)
}

lmm_rich <- df_CS %>%
  dplyr::group_by(season) %>%
  dplyr::group_modify(~ .get_lmm_r2(.x, "richness")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metric = "Richness (q = 0)")

lmm_shan <- df_CS %>%
  dplyr::group_by(season) %>%
  dplyr::group_modify(~ .get_lmm_r2(.x, "shannon")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(metric = "Shannon (q = 1)")

ds_season_lmm_r2 <- dplyr::bind_rows(lmm_rich, lmm_shan) %>%
  dplyr::relocate(metric, season)

# ---------- Display + Save ----------
if (P == 1) {
  cat("\n[Per-season OLS: slope ± SE, t, p/q, R², n]\n")
  print(ds_season_ols, n = nrow(ds_season_ols))
  cat("\n[Per-season LMM R² (marginal/conditional); NA if MuMIn not installed]\n")
  print(ds_season_lmm_r2, n = nrow(ds_season_lmm_r2))
}

readr::write_csv(ds_season_ols,     file.path(results_dir, paste0("ds_season_OLS_", stamp, ".csv")))
readr::write_csv(ds_season_lmm_r2,  file.path(results_dir, paste0("ds_season_LMM_R2_", stamp, ".csv")))

### ──────────────────────────────────────────────────────────────
### Autocorrelation Function (ACF) Check on Model Residuals
### ──────────────────────────────────────────────────────────────
# Extract the residuals from your primary Shannon diversity model
res_sh <- residuals(m_sh)

# Run the ACF test and plot the results
acf(res_sh, main = "ACF of Shannon-Stability Model Residuals")

# You can also check the richness model just to be thorough
res_rr <- residuals(m_rr)
acf(res_rr, main = "ACF of Richness-Stability Model Residuals")

### ──────────────────────────────────────────────────────────────
### Extract Numeric Output for ACF (Autocorrelation)
### ──────────────────────────────────────────────────────────────
# 1. Extract residuals from the main Shannon-Stability model
res_sh <- residuals(m_sh)

# 2. Compute ACF values (plot = FALSE to just get the numbers)
acf_out <- acf(res_sh, plot = FALSE)

# 3. Calculate the 95% confidence bounds
ci_bound <- 1.96 / sqrt(length(res_sh))

# 4. Print the results
cat("\n[ACF Numeric Results]\n")
cat("95% CI Bounds: +/-", round(ci_bound, 4), "\n\n")

# Create a clean table to show Lag, Autocorrelation, and Significance
acf_df <- data.frame(
  Lag = acf_out$lag[,1,1],
  Autocorrelation = round(acf_out$acf[,1,1], 4)
)

# Flag any lags (greater than 0) that exceed the CI bound
acf_df$Significant <- ifelse(acf_df$Lag > 0 & abs(acf_df$Autocorrelation) > ci_bound, "YES", "No")

print(acf_df)
