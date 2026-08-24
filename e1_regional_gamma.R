# ──────────────────────────────────────────────────────────────
#  Script:  f1_regional_gamma.R
#  Purpose: Load most recent diversity & stability CSV outputs,
#           extract γ-diversity (richness, Shannon),
#           and perform temporal trend analyses:
#             • Linear regressions (all seasons combined)
#             • GAM fits (per season, with CI ribbons)
#           Also generates summary tables of regression & GAM results.
#
#  Author:   Masami Fujiwara
#  Created:  2025-09-09
#
#  Dependencies:
#    • tidyverse  – data wrangling & plotting
#    • mgcv       – GAM models
#    • broom      – tidy model summaries
#    • patchwork  – (optional) arranging plots
#
#  Inputs:
#    • Latest CSVs from “results/” directory:
#         population_invariability, community_invariability,
#         portfolio_ps, synchrony_ps,
#         L5_alpha_diversity_observed,
#         L5_gamma_diversity_observed,
#         L5_beta_diversity_observed,
#         beta_bray_pairs_monthly,
#         beta_bray_bay_period_season
#
#  Outputs:
#    • Data frames:
#         – gamma_div (regional γ-diversity, richness & Shannon)
#         – lm_results (linear regression slopes, SE, q-values)
#         – gam_summary_table (parametric & smooth terms with q-values)
#    • Figures:
#         – p_richness (GAM smooths for γ Richness, faceted by season)
#         – p_shannon  (GAM smooths for γ Shannon, faceted by season)
#
#  Sections:
#    0.  Load dependencies & set working directory
#    1.  Identify most recent CSV files by base name
#    2.  Read and attach selected data frames
#    3.  Prepare γ-diversity data (add mid-year values, factor seasons)
#    4.  Linear regression (all seasons combined)
#    5.  GAM fits for richness & Shannon (per season)
#    6.  Extract tidy GAM summaries (parametric & smooth terms)
#    7.  Visualization and optional export
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
library(mgcv)

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
# This was written to keep all old results. Hand made version controlling scheme. 
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
#  1. REGIONAL TREND IN γ-DIVERSITY (Richness & Shannon)
###############################################################################
#   • Input: gamma_div (data frame with: period, season, richness, shannon)
#   • Analyses:
#       - Linear regression: richness & Shannon vs. mid-year (per season + overall)
#       - GAM: richness & Shannon vs. mid-year (per season)
#   • Output:
#       - Two figures: GAM smooths (4 panels each, one per season)
#       - One regression results table (linear models for both metrics)
###############################################################################

# ── Prepare data ──────────────────────────────────────────────────────────────
# Add numeric mid-year (for 3-yr periods: 1992–1994 → 1993, etc.)
period_midyears <- tibble(
  period = factor(paste0("P", 1:11), levels = paste0("P", 1:11)),
  mid_year = seq(1993, 2023, by = 3)
)

gamma_div <- gamma_div %>%
  left_join(period_midyears, by = "period")%>%
  mutate(
    season = factor(season, levels = c("Winter","Spring","Summer","Fall"))  # or levels = unique(season)
  )

# ── Linear regression results (ALL SEASONS COMBINED ONLY) ────────────────────
make_lm_table_overall <- function(df, response) {
  # drop missing values for the target metric
  d <- df %>%
    dplyr::filter(is.finite(.data[[response]]), is.finite(mid_year))
  
  # fit linear model
  fit <- lm(stats::as.formula(paste(response, "~ mid_year")), data = d)
  
  # extract slope, se, p
  slope <- coef(fit)["mid_year"]
  se    <- summary(fit)$coefficients["mid_year", "Std. Error"]
  p_val <- summary(fit)$coefficients["mid_year", "Pr(>|t|)"]
  q_val <- p.adjust(p_val, method = "BH")  # with one test, q == p
  
  tibble::tibble(
    metric   = response,
    n        = nrow(d),
    slope_lm = unname(slope),
    se_lm    = unname(se),
    q_lm     = unname(q_val)
  )
}

# Run for richness and shannon
lm_richness <- make_lm_table_overall(gamma_div, "richness")
lm_shannon  <- make_lm_table_overall(gamma_div, "shannon")

# Combine results
lm_results <- dplyr::bind_rows(lm_richness, lm_shannon)

print(lm_results)


# ── GAMs  ─────────────────────────────────────────────────────────
## Richness
gam_richness <- mgcv::gam(richness ~ s(mid_year, by = season) + season,
                          data = gamma_div, method = "REML")

pred_grid <- gamma_div %>%
  group_by(season) %>%
  reframe(mid_year = seq(min(mid_year), max(mid_year), length.out = 100)) %>%
  ungroup()

pr_rich <- predict(gam_richness, newdata = pred_grid, se.fit = TRUE)
pred_richness <- pred_grid %>%
  mutate(pred = pr_rich$fit, se = pr_rich$se.fit)


## Shannon
gam_shannon <- mgcv::gam(shannon ~ s(mid_year, by = season) + season,
                         data = gamma_div, method = "REML")

pr_shan <- predict(gam_shannon, newdata = pred_grid, se.fit = TRUE)
pred_shannon <- pred_grid %>%
  mutate(pred = pr_shan$fit, se = pr_shan$se.fit)

# Helper to extract GAM summary into a tidy table
# SAFE extractor for mgcv::gam summaries
extract_gam_summary <- function(model, metric) {
  s <- summary(model)
  
  # Parametric terms
  if (!is.null(s$p.table) && nrow(s$p.table) > 0) {
    param <- as.data.frame(s$p.table)
    param$term <- rownames(s$p.table)
    param <- tibble::as_tibble(param) |>
      dplyr::transmute(
        metric    = metric,
        component = "parametric",
        term      = term,
        edf       = NA_real_,
        statistic = `t value`,
        p_val     = `Pr(>|t|)`,
        estimate  = Estimate,
        std.error = `Std. Error`
      )
  } else {
    param <- tibble::tibble(
      metric=character(), component=character(), term=character(),
      edf=double(), statistic=double(), p_val=double(),
      estimate=double(), std.error=double()
    )
  }
  
  # Smooth terms
  if (!is.null(s$s.table) && nrow(s$s.table) > 0) {
    sm <- as.data.frame(s$s.table)
    sm$term <- rownames(s$s.table)
    sm <- tibble::as_tibble(sm) |>
      dplyr::transmute(
        metric    = metric,
        component = "smooth",
        term      = term,
        edf       = edf,
        statistic = F,
        p_val     = `p-value`,
        estimate  = NA_real_,
        std.error = NA_real_
      )
  } else {
    sm <- tibble::tibble(
      metric=character(), component=character(), term=character(),
      edf=double(), statistic=double(), p_val=double(),
      estimate=double(), std.error=double()
    )
  }
  
  dplyr::bind_rows(param, sm) |>
    dplyr::mutate(q_val = p.adjust(p_val, method = "BH")) |>
    dplyr::arrange(metric, component, term)
}


# ---- Use it on your models ----
rich_summ <- extract_gam_summary(gam_richness, "Richness")
shan_summ <- extract_gam_summary(gam_shannon,  "Shannon")

gam_summary_table <- dplyr::bind_rows(rich_summ, shan_summ)

extract_gam_summary <- function(model, metric) {
  s <- summary(model)

  # Parametric terms
  if (!is.null(s$p.table) && nrow(s$p.table) > 0) {
    param <- as.data.frame(s$p.table)
    param$term <- rownames(s$p.table)
    param <- tibble::as_tibble(param) |>
      dplyr::transmute(
        metric    = metric,
        component = "parametric",
        term      = term,
        edf       = NA_real_,
        statistic = `t value`,
        p_val     = `Pr(>|t|)`,
        estimate  = Estimate,
        std.error = `Std. Error`
      )
  } else {
    param <- tibble::tibble(
      metric=character(), component=character(), term=character(),
      edf=double(), statistic=double(), p_val=double(),
      estimate=double(), std.error=double()
    )
  }

  # Smooth terms
  if (!is.null(s$s.table) && nrow(s$s.table) > 0) {
    sm <- as.data.frame(s$s.table)
    sm$term <- rownames(s$s.table)
    sm <- tibble::as_tibble(sm) |>
      dplyr::transmute(
        metric    = metric,
        component = "smooth",
        term      = term,
        edf       = edf,
        statistic = F,
        p_val     = `p-value`,
        estimate  = NA_real_,
        std.error = NA_real_
      )
  } else {
    sm <- tibble::tibble(
      metric=character(), component=character(), term=character(),
      edf=double(), statistic=double(), p_val=double(),
      estimate=double(), std.error=double()
    )
  }

  dplyr::bind_rows(param, sm) |>
    dplyr::mutate(q_val = p.adjust(p_val, method = "BH")) |>
    dplyr::arrange(metric, component, term)
}


# Extract summaries
rich_summ <- extract_gam_summary(gam_richness, "Richness")
shan_summ <- extract_gam_summary(gam_shannon, "Shannon")

# Combine smooth + parametric terms
gam_summary_table <- dplyr::bind_rows(rich_summ, shan_summ)
print(gam_summary_table)

# Optional: select compact columns
gam_summary_compact <- gam_summary_table |>
  dplyr::select(metric, component, term, edf, statistic, p_val, q_val, estimate, std.error)
print(gam_summary_compact)


# # ── FIGURES  ─────────────────────────────────────────────────────────
# 
# p_richness <- ggplot(gamma_div, aes(x = mid_year, y = richness)) +
#   geom_point(alpha = 0.6) +
#   geom_line(data = pred_richness,
#             inherit.aes = FALSE,
#             aes(x = mid_year, y = pred),
#             color = "darkred") +
#   geom_ribbon(data = pred_richness,
#               inherit.aes = FALSE,
#               aes(x = mid_year, ymin = pred - 2*se, ymax = pred + 2*se),
#               fill = "pink", alpha = 0.3) +
#   facet_wrap(~ season) +
#   labs(x = "Mid-year", y = "γ Richness",
#        title = "GAM trends in γ Richness (per season)") +
#   theme_minimal()
# 
# p_shannon <- ggplot(gamma_div, aes(x = mid_year, y = shannon)) +
#   geom_point(alpha = 0.6) +
#   geom_line(data = pred_shannon,
#             inherit.aes = FALSE,
#             aes(x = mid_year, y = pred),
#             color = "darkgreen") +
#   geom_ribbon(data = pred_shannon,
#               inherit.aes = FALSE,
#               aes(x = mid_year, ymin = pred - 2*se, ymax = pred + 2*se),
#               fill = "lightgreen", alpha = 0.3) +
#   facet_wrap(~ season) +
#   labs(x = "Mid-year", y = "γ Shannon (H′)",
#        title = "GAM trends in γ Shannon (per season)") +
#   theme_minimal()
# 
# print(p_richness)
# print(p_shannon)

# ggsave("gamma_richness_GAM.png", p_richness, width = 8, height = 6, dpi = 300)
# ggsave("gamma_shannon_GAM.png", p_shannon, width = 8, height = 6, dpi = 300)
# write_csv(lm_results_out, "gamma_div_linear_regression_results.csv")


library(patchwork)

pt_mm   <- 0.3527778   # 1 point in mm
fig_w   <- 18          # cm (two-column)
fig_h   <- 12          # cm (adjust as you like; width stays 18 cm)

# Season order
season_levels <- c("Fall","Winter","Spring","Summer")
gamma_div     <- gamma_div %>% mutate(season = factor(season, levels = season_levels))
pred_richness <- pred_richness %>% mutate(season = factor(season, levels = season_levels))
pred_shannon  <- pred_shannon  %>% mutate(season = factor(season, levels = season_levels))

# Facet strip labels with tags
lab_rich <- setNames(paste0("(", letters[1:4], ") ", season_levels), season_levels)  # (a)-(d)
lab_shan <- setNames(paste0("(", letters[5:8], ") ", season_levels), season_levels)  # (e)-(h)

# ---- γ Richness (top row) ----
p_rich <- ggplot(gamma_div, aes(x = mid_year, y = richness)) +
  geom_point(color = "black", shape = 16, size = 0.8, alpha = 0.6) +
  geom_ribbon(data = pred_richness,
              inherit.aes = FALSE,
              aes(x = mid_year, ymin = pred - 2*se, ymax = pred + 2*se),
              fill = "grey80", alpha = 0.6) +
  geom_line(data = pred_richness,
            inherit.aes = FALSE,
            aes(x = mid_year, y = pred),
            color = "black", linewidth = pt_mm) +        # true 1-pt line
  facet_wrap(~ season, ncol = 4, labeller = labeller(season = lab_rich)) +
  labs(x = NULL, y = expression(gamma~Richness)) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position  = "none",
    plot.title       = element_blank(),
    strip.background = element_blank(),                  # <- remove box behind “Fall/Winter/…”
    strip.text       = element_text(color = "black", size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
    panel.spacing    = unit(3, "mm")
  )

# ---- γ Shannon (bottom row) ----
p_shan <- ggplot(gamma_div, aes(x = mid_year, y = shannon)) +
  geom_point(color = "black", shape = 16, size = 0.8, alpha = 0.6) +
  geom_ribbon(data = pred_shannon,
              inherit.aes = FALSE,
              aes(x = mid_year, ymin = pred - 2*se, ymax = pred + 2*se),
              fill = "grey80", alpha = 0.6) +
  geom_line(data = pred_shannon,
            inherit.aes = FALSE,
            aes(x = mid_year, y = pred),
            color = "black", linewidth = pt_mm) +        # true 1-pt line
  facet_wrap(~ season, ncol = 4, labeller = labeller(season = lab_shan)) +
  labs(x = "Mid-year", y = expression(gamma~Shannon~(H*minute))) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position  = "none",
    plot.title       = element_blank(),
    strip.background = element_blank(),                  # <- remove box
    strip.text       = element_text(color = "black", size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
    panel.spacing    = unit(3, "mm")
  )

# ---- Combine into one 8-panel figure (4 columns × 2 rows) ----
p_8panel <- p_rich / p_shan + plot_layout(heights = c(1, 1))

# ---- Save: TIFF, 18 cm wide, 600 dpi ----
results_dir <- file.path(script_dir, "results")

ggsave(
  filename = file.path(results_dir, "gamma_richness_shannon_8panels_minimal_18cm_600dpi.tif"),
  plot   = p_8panel,
  width  = fig_w, height = fig_h, units = "cm",
  dpi    = 600, device = "tiff", compression = "lzw"
)
