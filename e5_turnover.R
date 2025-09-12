# ──────────────────────────────────────────────────────────────
# Turnover (Bray–Curtis) vs Community Stability (I_C)
# Using updated strata-level turnover in `beta_bray`
# ──────────────────────────────────────────────────────────────
library(tidyverse)
library(lme4)
library(lmerTest)

# ──────────────────────────────────────────────────────────────
# 0.  Load dependencies
# ──────────────────────────────────────────────────────────────


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
# Turnover (Bray–Curtis; bray_bc) vs Community Stability (I_C)
# Using updated strata-level turnover in `beta_bray`
# ──────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)

# 1) Standardize columns & align turnover to I_C at the same (bay×period×season) stratum
beta_bray_std <- beta_bray %>%
  transmute(
    major_area,
    period,
    season = as.character(season),
    turnover = bray_bc,                 # <- updated column
    n_pairs,
    mean_months_per_pair
  ) %>%
  filter(is.finite(turnover))

# (Optional) QC filter if desired:
# beta_bray_std <- beta_bray_std %>% filter(n_pairs >= 2, mean_months_per_pair >= 3)

# Community invariability table
com_inv_std <- com_inv %>%
  select(major_area, period, season, I_C) %>%
  mutate(season = as.character(season)) %>%
  filter(is.finite(I_C))

# 2) Join turnover with I_C at the same stratum
df_turnover_ic <- com_inv_std %>%
  inner_join(beta_bray_std, by = c("major_area","period","season")) %>%
  filter(I_C > 0) %>%                                  # needed for log10
  mutate(
    season = factor(season, levels = c("Fall","Spring","Summer","Winter")),
    log10_Ic = log10(I_C)
  )

stopifnot(nrow(df_turnover_ic) > 0)

# 3) Mixed-effects model: log10(I_C) ~ turnover + season + (1|major_area)
mod_season <- lmer(log10_Ic ~ turnover + season + (1 | major_area),
                   data = df_turnover_ic)

mod_no_season <- lmer(log10_Ic ~ turnover + (1 | major_area),
                      data = df_turnover_ic)

# Likelihood-ratio test for season effect
lrt <- anova(mod_no_season, mod_season)

cat("\n[Mixed-effects model: log10(I_C) ~ turnover + season + (1|major_area)]\n")
print(summary(mod_season))
cat("\n[Likelihood-ratio test: add(season)]\n")
print(lrt)

# 4) Table-style summary with back-transformed effects
fixef_tab <- as.data.frame(coef(summary(mod_season)))
fixef_tab$term <- rownames(fixef_tab); rownames(fixef_tab) <- NULL

fixef_tab <- fixef_tab %>%
  transmute(
    term,
    estimate_log10 = Estimate,
    se = `Std. Error`,
    df = df,
    t  = `t value`,
    p  = `Pr(>|t|)`,
    effect_mult = 10^estimate_log10,          # multiplicative effect on I_C
    effect_pct  = (effect_mult - 1) * 100     # % change vs reference
  ) %>%
  mutate(
    term = dplyr::recode(
      term,
      `(Intercept)`  = "Intercept (Fall)",
      `seasonSpring` = "Spring",
      `seasonSummer` = "Summer",
      `seasonWinter` = "Winter",
      .default       = term
    )
  )

# Random-effect and residual variances
re_var <- as.data.frame(VarCorr(mod_season))
rand_info <- tibble::tibble(
  component = c("Random intercept (bay)", "Residual"),
  variance  = c(re_var$vcov[re_var$grp == "major_area"][1],
                re_var$vcov[re_var$grp == "Residual"][1]),
  sd        = sqrt(variance)
)

cat("\n[Fixed effects (with back-transformed effects and % change)]\n")
print(fixef_tab)
cat("\n[Random-effect and residual variances]\n")
print(rand_info)

# 5) Figure: population-level lines by season (fixed effects only)
xrng <- range(df_turnover_ic$turnover, na.rm = TRUE)
season_levels <- levels(df_turnover_ic$season)

pred_grid <- tidyr::crossing(
  turnover = seq(xrng[1], xrng[2], length.out = 200),
  season   = factor(season_levels, levels = season_levels)
)

pred_grid$log10_Ic_hat <- predict(mod_season, newdata = pred_grid, re.form = NA)

fig_turnover_stability <- ggplot() +
  geom_point(
    data = df_turnover_ic,
    aes(x = turnover, y = log10_Ic, color = season),
    alpha = 0.35, size = 1.5
  ) +
  geom_line(
    data = pred_grid,
    aes(x = turnover, y = log10_Ic_hat, color = season),
    linewidth = 1.1
  ) +
  labs(
    x = "Bray–Curtis turnover (bray_bc)",
    y = expression(log[10](I[C])),
 #   title = "Stability ~ Turnover across seasons (mixed model)",
    color = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(fig_turnover_stability)


## Figure for publication
season_levels <- c("Fall","Winter","Spring","Summer")
df_turnover_ic <- df_turnover_ic %>% mutate(season = factor(season, levels = season_levels))
pred_grid      <- pred_grid      %>% mutate(season = factor(season, levels = season_levels))

library(grid)

# 1) Tweak text sizes for an 8.5 cm-wide figure
base_pt <- 7  # good for single-column width

fig_turnover_stability_small <- ggplot() +
  geom_point(
    data = df_turnover_ic,
    aes(x = turnover, y = log10_Ic, color = season, shape = season),
    alpha = 0.5, size = 1.4, stroke = 0.4
  ) +
  geom_line(
    data = pred_grid,
    aes(x = turnover, y = log10_Ic_hat, color = season, linetype = season),
    linewidth = 0.3528   # true 1-pt line (in mm)
  ) +
  scale_color_grey(start = 0.15, end = 0.65, name = "Season") +
  scale_shape_manual(values = c(16, 17, 15, 3), name = "Season") +
  scale_linetype_manual(values = c("solid","dashed","dotdash","twodash"), name = "Season") +
  labs(
    x = "Bray–Curtis turnover (bray_bc)",
    y = expression(log[10](I[C]))
  ) +
  theme_minimal(base_size = base_pt) +
  theme(
    legend.position   = "bottom",
    axis.title        = element_text(size = base_pt),
    axis.text         = element_text(size = base_pt * 0.9),
    legend.title      = element_text(size = base_pt * 0.95),
    legend.text       = element_text(size = base_pt * 0.9),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.2, colour = "grey85"),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA)
  )

# 2) Save: TIFF, 8.5 cm width, 600 dpi, white background
results_dir <- file.path(script_dir, "results")

ggsave(filename = file.path(results_dir, "fig_turnover_stability_8p5cm_600dpi.tif"),
  plot = fig_turnover_stability_small,
  width = 8.5, height = 6.0, units = "cm",
  dpi = 600, device = "tiff", compression = "lzw",
  bg = "white"
)
# (Optional) Save outputs
# readr::write_csv(fixef_tab, file.path(results_dir, "turnover_vs_stability_LMM_table.csv"))
# ggsave(file.path(results_dir, "turnover_vs_stability_mixed_model.png"),
#        fig_turnover_stability, width = 7.5, height = 5.5, dpi = 300)

