# ──────────────────────────────────────────────────────────────
#  Section: Abundance–Variability Scaling (Taylor’s law)
#  Purpose:
#    • Examine the relationship between mean abundance (μ) and variance (σ²)
#      across species, seasons, and bays.
#    • Evaluate Taylor’s law: log(σ²) = a + b·log(μ).
#    • Express the equivalent invariability form:
#        Iᵢ = μ² / σ²  ⇒  log(Iᵢ) ~ m·log(μ), where m = 2 – b.
#    • Quantify how population stability (Iᵢ) scales with abundance.
#    • Visualize scaling with log–log plots and fit linear models
#      (overall and stratified by season).
#
#
#  Author:  Masami Fujiwara with assistance of ChatGPT 5.0 in debugging 
#     Most of the annotations were added by ChatGPT for readability. 
#  Date:    2026-08-24   
#
#  Dependencies (this section only):
#    • tidyverse   – data wrangling (dplyr, tibble) and plotting (ggplot2)
#    • patchwork   – combine multiple ggplot panels
#
#  Input:
#    • Data frame `pop_inv` containing population-level invariability metrics:
#        – mu   : mean abundance per species × season × bay
#        – I_i  : invariability (μ² / σ²)
#
#  Outputs:
#    • Model fits:
#        – lm(log10(σ²) ~ log10(μ)) → slope b (Taylor’s law)
#        – lm(log10(Iᵢ) ~ log10(μ)) → slope m (invariability form)
#    • Confidence intervals for slopes b and m
#    • Figures:
#        – Scatterplots with log–log scaling
#        – Smooth OLS fits by season
#        – Multi-panel figure (classic Taylor’s law vs. invariability form)
#
#  Ecological significance:
#    • Tests whether abundant species are disproportionately more stable.
#    • Links variance–mean scaling to stability theory across populations.
# ──────────────────────────────────────────────────────────────

rm(list = ls())

# Display figures? 1 YES, 0 NO
P <- 1

# ──────────────────────────────────────────────────────────────
# 0.  Load dependencies
# ──────────────────────────────────────────────────────────────
library(tidyverse)
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

###############################################################################
#  7. ABUNDANCE–VARIABILITY SCALING (Taylor’s law)
###############################################################################
#   • Using `pop_inv` (population_invariability)

pop_inv <- pop_inv %>% filter(is.finite(mu), is.finite(I_i))

pop_inv %>%
  filter(mu > 0, I_i > 0) %>%     # remove zeros
  ggplot(aes(x = mu, y = I_i)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  facet_wrap(~ season) +
  labs(
    x = "Mean abundance (μ_i)",
    y = "Population invariability (I_i)",
    title = "Scaling of population invariability (zeros removed)"
  ) 

summary(
  lm(log10(I_i) ~ log10(mu), 
     data = pop_inv %>% filter(mu > 0, I_i > 0))
)


pop_inv %>%
  filter(mu > 0, I_i > 0) %>%
  ggplot(aes(x = mu, y = I_i)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(
    method  = "lm",
    formula = y ~ x,
    se      = TRUE,
    color   = "steelblue",
    na.rm = TRUE
  ) +
  facet_wrap(~ season) +
  labs(
    x     = expression(Mean~abundance~(mu[i])),
    y     = expression(Population~invariability~(I[i])),
    title = "Scaling of population invariability per Taylor’s law"
  ) +
  theme_minimal()

if (P==1){
  print(pop_inv )
}

# --- Data prep (zeros removed; finite only) ---
df_taylor <- pop_inv %>%
  filter(is.finite(mu), is.finite(I_i), mu > 0, I_i > 0) %>%
  mutate(
    sigma2     = (mu^2) / I_i,      # implied variance
    log_mu     = log10(mu),
    log_sigma2 = log10(sigma2),
    log_I      = log10(I_i)
  )

### --- Pooled fits (optional: for reporting in caption) ---
### --- Pooled fits (optional: for reporting in caption) ---
fit_sigma <- lm(log_sigma2 ~ log_mu, data = df_taylor)   # classic Taylor: slope = b 
fit_I     <- lm(log_I ~ log_mu, data = df_taylor)        # invariability: slope = m = 2 - b

b    <- unname(coef(fit_sigma)["log_mu"]) 
ci_b <- unname(confint(fit_sigma)["log_mu", ]) 
t_b  <- summary(fit_sigma)$coefficients["log_mu", "t value"]
df_b <- fit_sigma$df.residual

m    <- unname(coef(fit_I)["log_mu"]) 
ci_m <- unname(confint(fit_I)["log_mu", ]) 
t_m  <- summary(fit_I)$coefficients["log_mu", "t value"]
df_m <- fit_I$df.residual

message(sprintf("Taylor slope b = %.3f (95%% CI [%.3f, %.3f], t_%d = %.1f)", b, ci_b[1], ci_b[2], df_b, t_b))
message(sprintf("Invariability slope m = %.3f (95%% CI [%.3f, %.3f], t_%d = %.1f)", m, ci_m[1], ci_m[2], df_m, t_m))

# --- Panel (a): Classic Taylor’s law: log10(σ²) vs log10(μ) ---
p_a <- ggplot(df_taylor, aes(x = mu, y = sigma2, color = season)) +
  geom_point(alpha = 0.25, size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(
    aes(x = mu, y = sigma2),
    method = "lm", formula = y ~ x, se = TRUE, linewidth = 1
  ) +
  labs(
    x = expression(Mean~abundance~(mu)),
    y = expression(Variance~(sigma^2))
    ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_a)

pt_mm <- 0.3527778  # 1 point in mm
pt_mm <- 0.3527778  # 1 pt in mm

p_a_bw <- ggplot(df_taylor, aes(x = mu, y = sigma2)) +
  # points by season in grey scale (tiny)
  geom_point(aes(color = season), alpha = 0.35, size = 0.1) +
  scale_color_grey(start = 0.15, end = 0.65, guide = "none") +
  # log–log axes
  scale_x_log10() +
  scale_y_log10() +
  # ONE overall regression line + grey SE ribbon (no season mapping here)
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    color = "black", fill = "grey70",
    linewidth = pt_mm
  ) +
  labs(
    x = expression(Mean~abundance~(mu)),
    y = expression(Variance~(sigma^2))
  ) +
  theme_minimal(base_size = 7) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey85"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# Save: TIFF, 8.5 cm wide, 600 dpi, white background
results_dir <- file.path(script_dir, "results")

ggsave(filename = file.path(results_dir, "Fig_6_taylor_loglog_8p5cm_600dpi.tif"),
  plot = p_a_bw,
  width = 8.5*1.5, height = 6.0*1.5, units = "cm",
  dpi = 600, device = "tiff", compression = "lzw",
  bg = "white"
)
