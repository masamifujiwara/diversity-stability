
# ──────────────────────────────────────────────────────────────
#  Script:  f4_synchrony.R
#  Purpose: 
#    • Load the most recent ecological diversity & stability CSV files
#      from the “results/” directory
#    • Assemble data frames for α, β, γ diversity, synchrony, portfolio, 
#      and stability metrics
#    • Fit piecewise Structural Equation Models (SEMs) linking diversity, 
#      synchrony, and community stability using mixed-effects models
#    • Quantify direct and indirect pathways of richness and Shannon 
#      diversity on stability
#    • Produce diagnostic summaries and ggplot2 figures for interpretation
#
#  Author:   Masami Fujiwara
#  Created:  2025-08-07
#
#  Dependencies (loaded here):
#    • tidyverse   – data wrangling (dplyr, tibble, purrr, stringr, readr) and ggplot2 graphics
#    • nlme        – linear mixed-effects models for SEM component models
#    • piecewiseSEM– piecewise SEM fitting and evaluation (psem, coefs, rsquared, fisherC)
#    • patchwork   – combining multiple ggplot2 figures into panels
#
#  Inputs:
#    • CSV files in “results/” directory, named as
#        <base_name>_YYYY-MM-DD.csv
#      where <base_name> ∈ {
#        population_invariability, community_invariability,
#        portfolio_ps, synchrony_ps,
#        L5_alpha_diversity_observed, L5_gamma_diversity_observed,
#        L5_beta_diversity_observed, beta_bray_pairs_monthly,
#        beta_bray_bay_period_season
#      }
#
#  Outputs (in global environment):
#    • Data frames:
#        – alpha_div, gamma_div, beta_div, beta_bray, synchrony,
#          portfolio, pop_inv, com_inv
#    • SEM objects:
#        – sem_rich (richness-based SEM)
#        – sem_shan (Shannon-based SEM)
#    • SEM summaries:
#        – raw & standardized coefficients
#        – R² values
#        – Fisher’s C tests
#        – indirect effect estimates
#    • ggplot figures:
#        – fig_logphi_rich, fig_logphi_shan, fig_logphi_both
#
#  Notes:
#    • Synchrony (φ) is log-transformed with epsilon adjustment to handle zeros
#    • Richness and Shannon models are run separately to compare direct vs.
#      indirect stabilizing effects
#    • Figures display relationships between diversity and synchrony, and
#      are combined into multi-panel layouts with patchwork
# ──────────────────────────────────────────────────────────────

rm(list = ls())

# Display figures? 1 YES, 0 NO
P <- 1

# ──────────────────────────────────────────────────────────────
# 0.  Load dependencies
# ──────────────────────────────────────────────────────────────
library(tidyverse)
library(nlme)           # lme()
library(piecewiseSEM)   # psem(), coefs(), rsquared(), fisherC()
library(patchwork)

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
#   • Scope: Local‐scale (α)  diversity
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - richness : q = 0 richness, standardized to a common sample size
#       - shannon  : q = 1 Shannon diversity on  counts
#   • Use: Controls for unequal sampling when tracking how local species diversity changes
#     across subregions (areas), time periods, and seasons.
#
# Data frame: gamma_div
#   • Scope: Regional‐scale (γ)  diversity within each subregion
#   • Rows indexed by: major_area, period, season
#   • Key columns:
#       - richness : pooled species richness (q = 0)
#       - shannon  : pooled Shannon diversity (q = 1)
#   • Use: Assesses how the total species pool within each area shifts over time/season
#     once sampling effort is equalized.
#
# Data frame: beta_div
#   • Scope: Spatial turnover (β) of  diversity among sites
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
# Log-transform synchrony (φ), SEMs (richness vs. Shannon), and plots
# ──────────────────────────────────────────────────────────────

## 1) Analysis dataset with log-transformed synchrony
##    (Guard against φ = 0 by adding a tiny epsilon before log)
eps <- 1e-6

df_sem <- com_inv %>%
  left_join(synchrony, by = c("major_area","period","season")) %>%
  left_join(alpha_div, by = c("major_area","period","season")) %>%
  transmute(
    major_area,
    season,
    period,
    #I_C = community_invariability,
    phi,                               # synchrony in [0,1]
    log_phi = log(pmax(phi, eps)),     # log-transform
    log_IC  = log10(I_C),              # stability on log10 scale (as in prior work)
    richness,
    shannon
  ) %>%
  filter(is.finite(log_phi), is.finite(log_IC),
         is.finite(richness), is.finite(shannon))

# ensure grouping factors
df_sem <- df_sem %>%
  mutate(
    major_area = factor(major_area),
    season     = factor(season)
  )

## 2) SEM with RICHNESS
m_sync_rich <- lme(
  log_phi ~ richness,
  random = list(major_area = ~1, season = ~1),
  data   = df_sem,
  na.action = na.omit,
  method = "REML"
)

m_ic_rich <- lme(
  log_IC ~ richness + log_phi,
  random = list(major_area = ~1, season = ~1),
  data   = df_sem,
  na.action = na.omit,
  method = "REML"
)

sem_rich <- psem(m_sync_rich, m_ic_rich, data = df_sem)

# Useful outputs
sem_rich_coefs_raw  <- coefs(sem_rich, standardize = "none")
sem_rich_coefs_std  <- coefs(sem_rich, standardize = "scale")  # standardized (β)
sem_rich_r2         <- rsquared(sem_rich)
sem_rich_fisherC    <- fisherC(sem_rich)

if (P == 1) {
  message("\n[SEM (Richness) – raw coefficients]"); print(sem_rich_coefs_raw)
  message("\n[SEM (Richness) – standardized coefficients]"); print(sem_rich_coefs_std)
  message("\n[SEM (Richness) – R^2]"); print(sem_rich_r2)
  message("\n[SEM (Richness) – Fisher's C test]"); print(sem_rich_fisherC)
}

## 3) SEM with SHANNON
m_sync_shan <- lme(
  log_phi ~ shannon,
  random = list(major_area = ~1, season = ~1),
  data   = df_sem,
  na.action = na.omit,
  method = "REML"
)

m_ic_shan <- lme(
  log_IC ~ shannon + log_phi,
  random = list(major_area = ~1, season = ~1),
  data   = df_sem,
  na.action = na.omit,
  method = "REML"
)

sem_shan <- psem(m_sync_shan, m_ic_shan, data = df_sem)

# Useful outputs
sem_shan_coefs_raw  <- coefs(sem_shan, standardize = "none")
sem_shan_coefs_std  <- coefs(sem_shan, standardize = "scale")
sem_shan_r2         <- rsquared(sem_shan)
sem_shan_fisherC    <- fisherC(sem_shan)

if (P == 1) {
  message("\n[SEM (Shannon) – raw coefficients]"); print(sem_shan_coefs_raw)
  message("\n[SEM (Shannon) – standardized coefficients]"); print(sem_shan_coefs_std)
  message("\n[SEM (Shannon) – R^2]"); print(sem_shan_r2)
  message("\n[SEM (Shannon) – Fisher's C test]"); print(sem_shan_fisherC)
}

## 4) Plots: log(synchrony) vs diversity (separate + combined)

# (a) Richness vs log(φ)
fig_logphi_rich <- ggplot(df_sem, aes(x = richness, y = log_phi)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  labs(
    x = "Richness (q = 0)",
    y = expression(log(phi))
  ) +
  theme_minimal()

# (b) Shannon vs log(φ)
fig_logphi_shan <- ggplot(df_sem, aes(x = shannon, y = log_phi)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  labs(
    x = "Shannon diversity (q = 1)",
    y = expression(log(phi))
  ) +
  theme_minimal()

# Combined two-panel figure
tag_theme <- theme(plot.tag = element_text(face = "bold"))
fig_logphi_both <- (fig_logphi_rich + labs(tag = "a") + tag_theme) /
  (fig_logphi_shan + labs(tag = "b") + tag_theme)

results_dir <- file.path(script_dir, "results")

ggsave(filename = file.path(results_dir, "fig_Synchrony_both.tif"),
  fig_logphi_both, width = 8.5, height = 12, units = "cm",
       dpi = 600, device = "tiff", compression = "lzw")


if (P == 1) {
  print(fig_logphi_rich)
  print(fig_logphi_shan)
  print(fig_logphi_both)
}

## 5) (Optional) Indirect effects via synchrony
##    Unstandardized indirect effect = (diversity → log_phi) * (log_phi → log_IC)
get_indirect <- function(sem_obj, div_term) {
  tb <- coefs(sem_obj, standardize = "none")
  a  <- tb$Estimate[tb$Response == "log_phi" & tb$Predictor == div_term]
  b  <- tb$Estimate[tb$Response == "log_IC"  & tb$Predictor == "log_phi"]
  a * b
}
indirect_rich <- get_indirect(sem_rich, "richness")
indirect_shan <- get_indirect(sem_shan, "shannon")
if (P == 1) {
  message("\n[Indirect effect via log_phi]  Richness: ", round(indirect_rich, 4),
          "   Shannon: ", round(indirect_shan, 4))
}

