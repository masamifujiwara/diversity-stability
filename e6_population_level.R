# ──────────────────────────────────────────────────────────────
#  Script:  e_diversity_stability_analysis.R
#  Section: Population-Level Mixed-Effects Model
#
#  Purpose:
#    • Fit linear mixed-effects models to population-level stability (Iᵢ).
#    • Partition variation in invariability (log10(Iᵢ)) across
#      – Seasons (fixed effects),
#      – Bays (random intercepts by major_area),
#      – Species identity (random intercepts by species_code).
#    • Quantify relative contributions of seasonality, spatial variation,
#      and species identity to population stability.
#    • Visualize seasonal differences, bay-level random effects,
#      and species-level random effects.
#
#  Author:   Masami Fujiwara
#  Created:  2025-08-07
#
#  Dependencies (this section only):
#    • tidyverse   – data wrangling (dplyr) and visualization (ggplot2)
#    • lme4        – fit mixed-effects models (lmer)
#    • lmerTest    – add denominator df and p-values to lmer outputs
#    • broom.mixed – extract and tidy random-effects estimates for plotting
#
#  Inputs:
#    • Data frame `pop_inv` containing:
#        – I_i          : population-level invariability (μ² / σ²)
#        – species_code : unique species identifier
#        – major_area   : bay / spatial unit
#        – season       : categorical (Fall, Spring, Summer, Winter)
#    • Optional: `species` lookup table (species_code ↔ common name)
#
#  Outputs:
#    • Model object:
#        – m_pop (lmer fit: log10(Iᵢ) ~ season + (1|major_area) + (1|species_code))
#    • Random-effects tables (bay and species)
#    • Figures:
#        – fig_pop_season : seasonal boxplots of log10(Iᵢ)
#        – fig_bay_re     : caterpillar plot of bay random intercepts
#        – fig_spec_re    : caterpillar plot of species random intercepts (by code)
#        – fig_spec_re_c  : caterpillar plot of species random intercepts (by common name)
#
#  Ecological significance:
#    • Reveals how stability varies by season, species, and spatial unit.
#    • Identifies species with unusually high or low baseline stability.
#    • Provides evidence on whether seasonal or species-level processes
#      dominate population stability patterns in coastal communities.
# ──────────────────────────────────────────────────────────────

rm(list = ls())

# Display figures? 1 YES, 0 NO
P <- 1

# ──────────────────────────────────────────────────────────────
# 0.  Load dependencies
# ──────────────────────────────────────────────────────────────
library(tidyverse)
library(broom.mixed)
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


###############################################################################
#  8. POPULATION-LEVEL MIXED MODEL
###############################################################################
#   • Season + random effects of bay & species_code

m_pop <- lmer(log10(I_i) ~ season + (1 | major_area) + (1 | species_code),
              data = pop_inv)
summary(m_pop)
ranef(m_pop)$major_area

# Linear mixed model for log10(I_i) ~ season + (1 | major_area) + (1 | species_code)
#
# Random effects:
#   species_code (Intercept):
#     Variance = 0.08759, SD = 0.29596
#       → Species identity accounts for substantial variation in population stability.
#   major_area (Intercept):
#     Variance = 0.000522, SD = 0.02286
#       → Bay-to-bay differences are minimal.
#   Residual:
#     Variance = 0.10651, SD = 0.32636
#
# Fixed effects (reference = Fall):
#   (Intercept):
#     Estimate = -2.056, SE = 0.0196, t = -104.76, p < 2e-16
#       → Baseline log10(I_i) in Fall corresponds to I_i ≈ 10^-2.056 ≈ 0.0088.
#   seasonSpring:
#     Estimate = -0.0356, SE = 0.0058, t = -6.13, p < 1e-9
#       → Spring stability ≈ 10^-0.0356 ≈ 0.92 × Fall (8% lower).
#   seasonSummer:
#     Estimate =  0.0280, SE = 0.0055, t =  5.07, p < 1e-6
#       → Summer stability ≈ 10^0.0280 ≈ 1.07 × Fall (7% higher).
#   seasonWinter:
#     Estimate = -0.1551, SE = 0.0063, t = -24.67, p < 2e-16
#       → Winter stability ≈ 10^-0.1551 ≈ 0.70 × Fall (30% lower).
#
# Random intercepts for major_area (deviation in log10(I_i)):
#   Bay 1: -0.01733 (≈4% lower than average Fall stability)
#   Bay 2:  0.03039 (≈7% higher)
#   Bay 3:  0.00872 (≈2% higher)
#   Bay 4:  0.02096 (≈5% higher)
#   Bay 5:  0.01094 (≈2.5% higher)
#   Bay 6: -0.00116 (≈0% difference)
#   Bay 7: -0.03575 (≈8% lower)
#   Bay 8: -0.01677 (≈4% lower)
#
# Overall interpretation:
#   • Population stability varies significantly by season:
#       – Highest in Summer, intermediate in Fall, lower in Spring, lowest in Winter.
#   • Species identity is the main source of variability; spatial (bay) effects are minor.



# 1. Seasonal differences in population stability
fig_pop_season <- pop_inv %>%
  filter(I_i > 0) %>% 
  ggplot(aes(x = season, y = log10(I_i), fill = season)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    x     = "Season",
    y     = expression(log[10](I[i])),
    title = "Population-level Invariability by Season"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 2. Caterpillar plot of bay random intercepts
bay_re <- tidy(m_pop, effects = "ran_vals", conf.int = TRUE) %>%
  filter(group == "major_area", term == "(Intercept)")

fig_bay_re <- ggplot(bay_re, aes(x = estimate, y = level)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = expression(
      paste(
        "Bay random intercept (", 
        log[10](I[i]), 
        ") ", 
        "\u00B1 95% CI"
      )
    ),
    y     = "Major area",
    title = "Between-Bay Variation in Population Stability"
  )+
  theme_minimal()

# 3. Caterpillar plot of species random intercepts (top 15 most extreme)
spec_re <- tidy(m_pop, effects = "ran_vals", conf.int = TRUE) %>%
  filter(group == "species_code", term == "(Intercept)") %>%
  arrange(estimate) %>%
  { 
    n_tot <- nrow(.)
    if (n_tot <= 30) . else slice(., c(1:15, (n_tot-14):n_tot))
  } %>%
  dplyr::mutate(level = factor(level, levels = rev(level)))

load("data.Rdata")

# 1) join on species code, assuming 'level' holds the code:
spec_re2 <- spec_re %>%
  dplyr::mutate(level = as.character(level)) %>%
  left_join(
    species %>%
      dplyr::mutate(species_code = as.character(species_code)) %>%
      dplyr::select(species_code, com_name),
    by = c("level" = "species_code")
  ) %>%
  dplyr::mutate(com_name = factor(com_name, levels = com_name))

# 3) plot using the common name
fig_spec_re_c <- ggplot(spec_re2, aes(x = estimate, y = com_name)) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = expression(
      paste(
        "Species random intercept (",
        log[10](I[i]),
        ") ± 95% CI"
      )
    ),
    y     = "Common name"  ) +
  theme_minimal(base_size = 10)


fig_spec_re <- ggplot(spec_re, aes(x = estimate, y = level)) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = expression(
      paste(
        "Species random intercept (",
        log[10](I[i]),
        ") ",
        "\u00B1 95% CI"
      )
    ),
    y     = "Species code",
    title = "Top 15 Species with Highest and Lowest Baseline Stability"
  )+
  theme_minimal(base_size = 8)

# To display:
if (P==1){
  print(fig_pop_season)
  print(fig_bay_re)
  print(fig_spec_re)
  print(fig_spec_re_c)
}

results_dir <- file.path(script_dir, "results")

ggsave(filename = file.path(results_dir, "population_variability_600dpi.tif"),
  plot = fig_spec_re_c,
  width = 18, height = 10.0, units = "cm",
  dpi = 600, device = "tiff", compression = "lzw",
  bg = "white"
)
