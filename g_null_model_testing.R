### ──────────────────────────────────────────────────────────────
### Null Model Test for Statistical Coupling
###
###  Revision note (what changed relative to the previous version):
###   1. The randomization now permutes each species' MONTHLY abundance
###      time series within a stratum, rather than its 180 individual haul
###      counts. Synchrony (Eq. 6) and community invariability (Eq. 3) are
###      both defined on monthly series, so the null must be generated at
###      that same grain. Permuting monthly totals preserves each species'
###      monthly variance EXACTLY, which is the quantity that enters the
###      denominator of Eq. 6, while still destroying between-species
###      covariance. The previous version computed the denominator from
###      haul-level variances while the numerator came from monthly totals,
###      which put the two halves of phi on different scales and produced
###      null phi values outside the [0, 1] bound of the Loreau-de
###      Mazancourt index.
###   2. The observed SEM slopes are now estimated inside this script using
###      the same function that fits the null models, instead of being typed
###      in by hand. Observed and null coefficients are therefore guaranteed
###      to come from identical code.
###   3. Permutation p-values use the (1 + r) / (1 + n) estimator, which
###      cannot return exactly 0 or 1.
###   4. Assertions verify that the permutation preserves species variances
###      and that null phi stays within [0, 1].
###
###
###  Author:  Masami Fujiwara with assistance of ChatGPT 5.0 in debugging 
###     The original code had errors. Code revisions were suggested by
###     Claude (Opus 4.9). Then, the revisions were implemented.
###  Date:    2026-08-24   
### ──────────────────────────────────────────────────────────────

library(tidyverse)
library(nlme)

### ──────────────────────────────────────────────────────────────
### 1. Environment Setup & Data Loading
### ──────────────────────────────────────────────────────────────

script_dir <- NULL
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  try({
    path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(path) && path != "") script_dir <- dirname(path)
  }, silent = TRUE)
}
if (is.null(script_dir)) {
  try({
    if (!is.null(sys.frames()[[1]]$ofile))
      script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile))
  }, silent = TRUE)
}
if (is.null(script_dir)) script_dir <- getwd()

setwd(script_dir)
results_dir <- file.path(script_dir, "results")

# --- A. Load alpha_div from the latest CSV ---
alpha_files <- list.files(results_dir,
                          pattern = "^L5_alpha_diversity_observed_\\d{4}-\\d{2}-\\d{2}\\.csv$", full.names = TRUE)
if (length(alpha_files) == 0)
  stop("L5_alpha_diversity_observed CSV not found in results folder.")
alpha_file_latest <- sort(alpha_files, decreasing = TRUE)[1]
alpha_div <- read_csv(alpha_file_latest, show_col_types = FALSE)
cat("Loaded alpha_div from:", basename(alpha_file_latest), "\n")

# --- B. Reconstruct monthly_sp_valid from data.Rdata ---
load("data.Rdata")
station <- STATION_BS %>% dplyr::select(sample_id, major_area, year, month) %>%
  filter(major_area %in% 1:8)
data <- bagseine %>% dplyr::select(sample_id, major_area, year, month, catch, species_code)

rare_sp <- data %>% count(species_code, wt = (catch > 0)) %>% filter(n <= 5) %>% pull(species_code)
data_trim <- data %>% filter(!species_code %in% rare_sp)

sample_species <- data_trim %>% group_by(sample_id, species_code) %>%
  summarise(abund = sum(catch), .groups = "drop")
species_wide <- sample_species %>%
  pivot_wider(names_from = species_code, values_from = abund,
              values_fill = 0, names_prefix = "sp_")
final_df_wide <- station %>% distinct(sample_id, major_area, year, month) %>%
  left_join(species_wide, by = "sample_id") %>%
  mutate(across(starts_with("sp_"), ~ replace_na(.x, 0)))

start_year <- 1992
end_year   <- 2024
period_lengths <- rep(3, 11)
year_period_lut <- tibble(
  year = seq(start_year, end_year),
  period_num = rep(seq_along(period_lengths), times = period_lengths)
) %>% mutate(period = factor(paste0("P", period_num), levels = paste0("P", 1:11)))

monthly_sp <- final_df_wide %>%
  filter(year >= start_year, year <= end_year) %>%
  left_join(year_period_lut, by = "year") %>%
  pivot_longer(cols = starts_with("sp_"), names_to = "species_code",
               names_pattern = "sp_(.*)",
               names_transform = list(species_code = as.integer),
               values_to = "abund") %>%
  mutate(species_code = as.integer(species_code),
         date = as.Date(sprintf("%d-%02d-15", year, month)),
         season = case_when(month %in% c(12, 1, 2) ~ "Winter",
                            month %in% 3:5         ~ "Spring",
                            month %in% 6:8         ~ "Summer",
                            month %in% 9:11        ~ "Fall"),
         season = factor(season, levels = c("Winter","Spring","Summer","Fall"))) %>%
  filter(!is.na(period))

station_strat <- station %>% filter(year >= start_year, year <= end_year) %>%
  left_join(year_period_lut %>% select(year, period), by = "year") %>%
  mutate(season = case_when(month %in% c(12, 1, 2) ~ "Winter",
                            month %in% 3:5         ~ "Spring",
                            month %in% 6:8         ~ "Summer",
                            month %in% 9:11        ~ "Fall"),
         season = factor(season, levels = c("Winter","Spring","Summer","Fall"))) %>%
  filter(!is.na(period))

valid_strata <- station_strat %>%
  group_by(major_area, period, season) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop") %>%
  filter(n_samples >= 180) %>%                       # >= rather than == (see audit item G-4)
  select(major_area, period, season)

monthly_sp_valid <- monthly_sp %>%
  semi_join(valid_strata, by = c("major_area","period","season"))
cat("Constructed monthly_sp_valid:", nrow(valid_strata), "fully sampled strata.\n\n")

### ──────────────────────────────────────────────────────────────
### 2. Monthly species time series (the grain at which phi is defined)
### ──────────────────────────────────────────────────────────────

sp_monthly <- monthly_sp_valid %>%
  group_by(major_area, period, season, date, species_code) %>%
  summarise(abund_date = sum(abund), .groups = "drop")

strata_key <- sp_monthly %>% distinct(major_area, period, season) %>%
  mutate(sid = row_number())

sp_monthly <- sp_monthly %>%
  left_join(strata_key, by = c("major_area","period","season")) %>%
  group_by(sid) %>% mutate(mo = as.integer(factor(date))) %>% ungroup()

n_months <- n_distinct(sp_monthly$mo)
stopifnot(n_months == 9)

# Species present at all in a stratum. Species absent throughout contribute
# zero to both the numerator and the denominator of Eq. 6 and are inert.
active <- sp_monthly %>% group_by(sid, species_code) %>%
  filter(sum(abund_date) > 0) %>% ungroup()

Xwide <- active %>% select(sid, species_code, mo, abund_date) %>%
  pivot_wider(names_from = mo, values_from = abund_date, values_fill = 0) %>%
  arrange(sid, species_code)
sid_vec <- Xwide$sid
X <- as.matrix(Xwide[, as.character(seq_len(n_months))])
N_rows <- nrow(X)
cat("Species x stratum series:", N_rows, "rows x", n_months, "months\n")

### ──────────────────────────────────────────────────────────────
### 3. Per-stratum constants (unchanged by the permutation)
### ──────────────────────────────────────────────────────────────
# mu_C is a sum of species means and so is invariant to permutation.
# sum_sd_i is a sum of species standard deviations, each of which a
# within-species permutation preserves exactly. Both are therefore held
# fixed across simulations by construction rather than by assumption.

sd_i <- apply(X, 1, sd)
stratum_base <- tibble(sid = sid_vec, sd_i = sd_i) %>%
  group_by(sid) %>% summarise(sum_sd_i = sum(sd_i), .groups = "drop") %>%
  left_join(strata_key, by = "sid") %>%
  left_join(alpha_div %>% select(major_area, period, season, shannon),
            by = c("major_area","period","season")) %>%
  mutate(mu_C = as.numeric(rowsum(rowSums(X), sid_vec)) / n_months,
         major_area = factor(major_area),
         season     = factor(as.character(season))) %>%
  arrange(sid)

### ──────────────────────────────────────────────────────────────
### 4. Helpers: permutation, metric assembly, SEM fitting
### ──────────────────────────────────────────────────────────────

# Independently permute each row (species x stratum series) across months.
grp_cm <- rep(seq_len(N_rows), times = n_months)   # column-major row index
col_cm <- rep(seq_len(n_months), each = N_rows)    # column-major column index
row_rm <- rep(seq_len(N_rows), each = n_months)    # row-major row index

permute_rows <- function() {
  ord <- order(grp_cm, runif(N_rows * n_months))   # random order within each row
  matrix(X[cbind(row_rm, col_cm[ord])],
         nrow = N_rows, ncol = n_months, byrow = TRUE)
}

# Assemble the stratum-level metrics from a species x month matrix.
build_metrics <- function(M) {
  comm <- rowsum(M, sid_vec)                       # strata x months community totals
  stratum_base %>%
    mutate(var_C   = apply(comm, 1, var),
           I_C     = if_else(var_C > 0, mu_C^2 / var_C, NA_real_),
           phi     = if_else(sum_sd_i > 0, var_C / (sum_sd_i^2), NA_real_),
           log_IC  = log10(I_C),
           log_phi = log(pmax(phi, 1e-6)))
}

# Fit the two SEM submodels on standardized variables and return STD betas.
fit_std_sem <- function(d) {
  d <- d %>%
    filter(is.finite(log_IC), is.finite(log_phi), is.finite(shannon)) %>%
    mutate(shannon_std = as.numeric(scale(shannon)),
           log_phi_std = as.numeric(scale(log_phi)),
           log_IC_std  = as.numeric(scale(log_IC)))
  ctl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
  m_sync <- tryCatch(lme(log_phi_std ~ shannon_std,
                         random = list(major_area = ~1, season = ~1),
                         data = d, na.action = na.omit, method = "REML", control = ctl),
                     error = function(e) NULL)
  m_stab <- tryCatch(lme(log_IC_std ~ shannon_std + log_phi_std,
                         random = list(major_area = ~1, season = ~1),
                         data = d, na.action = na.omit, method = "REML", control = ctl),
                     error = function(e) NULL)
  c(synchrony = if (is.null(m_sync)) NA_real_ else unname(fixef(m_sync)["shannon_std"]),
    stability = if (is.null(m_stab)) NA_real_ else unname(fixef(m_stab)["shannon_std"]))
}

### ──────────────────────────────────────────────────────────────
### 5. Observed slopes, estimated with the same code as the null
### ──────────────────────────────────────────────────────────────

observed <- build_metrics(X)
stopifnot(all(observed$phi >= 0 & observed$phi <= 1, na.rm = TRUE))
obs <- fit_std_sem(observed)
cat("\n[Observed standardized slopes, from this script]\n")
cat("  Shannon -> log(phi)   :", round(obs["synchrony"], 4), "\n")
cat("  Shannon -> log10(I_C) :", round(obs["stability"], 4), "\n")
cat("  (these should match the piecewiseSEM output in e4_synchrony.R)\n")

### ──────────────────────────────────────────────────────────────
### 6. Null model
### ──────────────────────────────────────────────────────────────

n_sims <- 1000
set.seed(2025)

# One-off verification that the permutation behaves as intended.
P_chk <- permute_rows()
stopifnot(all(abs(rowSums(P_chk)   - rowSums(X))       < 1e-8))   # totals preserved
stopifnot(all(abs(apply(P_chk,1,var) - apply(X,1,var)) < 1e-8))   # variances preserved
chk <- build_metrics(P_chk)
stopifnot(all(chk$phi >= 0 & chk$phi <= 1, na.rm = TRUE))         # phi stays bounded
cat("\nPermutation checks passed (species totals, variances, and phi bound).\n")

cat("Starting Null Model Simulations...\n")
null_results <- tibble(sim = integer(),
                       slope_stability_std = numeric(),
                       slope_synchrony_std = numeric())
for (i in seq_len(n_sims)) {
  if (i %% 100 == 0) cat("  simulation", i, "of", n_sims, "\n")
  s <- fit_std_sem(build_metrics(permute_rows()))
  null_results <- bind_rows(null_results,
                            tibble(sim = i, slope_stability_std = s["stability"], slope_synchrony_std = s["synchrony"]))
}

### ──────────────────────────────────────────────────────────────
### 7. Compare null slopes to observed
### ──────────────────────────────────────────────────────────────

obs_slope_stability <- unname(obs["stability"])
obs_slope_synchrony <- unname(obs["synchrony"])

# (1 + r) / (1 + n) permutation p-values
p_val_stability <- (1 + sum(null_results$slope_stability_std >= obs_slope_stability, na.rm = TRUE)) /
  (1 + sum(!is.na(null_results$slope_stability_std)))
p_val_synchrony <- (1 + sum(null_results$slope_synchrony_std <= obs_slope_synchrony, na.rm = TRUE)) /
  (1 + sum(!is.na(null_results$slope_synchrony_std)))

ci <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)

cat("\n[Null Model Results]\n")
cat(sprintf("Stability  observed %.4f | null mean %.4f | null 95%% CI [%.4f, %.4f] | p = %.4f\n",
            obs_slope_stability, mean(null_results$slope_stability_std, na.rm = TRUE),
            ci(null_results$slope_stability_std)[1], ci(null_results$slope_stability_std)[2],
            p_val_stability))
cat(sprintf("Synchrony  observed %.4f | null mean %.4f | null 95%% CI [%.4f, %.4f] | p = %.4f\n",
            obs_slope_synchrony, mean(null_results$slope_synchrony_std, na.rm = TRUE),
            ci(null_results$slope_synchrony_std)[1], ci(null_results$slope_synchrony_std)[2],
            p_val_synchrony))
cat("Simulations that converged:",
    sum(!is.na(null_results$slope_stability_std)), "of", n_sims, "\n")

### ──────────────────────────────────────────────────────────────
### 8. Save
### ──────────────────────────────────────────────────────────────
today_tag <- format(Sys.Date(), "%Y-%m-%d")
write_csv(null_results, file.path(results_dir, paste0("null_model_slopes_std_", today_tag, ".csv")))

null_summary <- tibble(
  metric              = c("Stability", "Synchrony"),
  observed_slope_std  = c(obs_slope_stability, obs_slope_synchrony),
  mean_null_slope_std = c(mean(null_results$slope_stability_std, na.rm = TRUE),
                          mean(null_results$slope_synchrony_std, na.rm = TRUE)),
  null_ci_low         = c(ci(null_results$slope_stability_std)[1],
                          ci(null_results$slope_synchrony_std)[1]),
  null_ci_high        = c(ci(null_results$slope_stability_std)[2],
                          ci(null_results$slope_synchrony_std)[2]),
  p_value             = c(p_val_stability, p_val_synchrony),
  n_sims              = n_sims
)
write_csv(null_summary, file.path(results_dir, paste0("null_model_summary_std_", today_tag, ".csv")))
print(null_summary)

cat("\nNull model results and summary saved to the 'results' directory.\n")
