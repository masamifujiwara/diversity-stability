###############################################################################
#  Script purpose
#  --------------
#  Diversity & Stability: Population and Community Invariability
#
#  1. Pre-process Texas bay bag-seine data
#        • drop ultra-rare species (≤ 5 non-zero records)
#        • build a complete sample × species matrix (zeros = true absences)
#        • convert to long format and add period (3-year bins) and season (DJF/MAM/JJA/SON)
#
#  2. Compute invariability metrics within each bay × period × season:
#        • Population invariability (I_i)
#        • Community invariability (I_C)
#        • Synchrony (φ_LM)
#        • Portfolio effect (I_C / mean I_i, abundance-weighted)
#
#  Outputs
#  -------
#   1) population_invariability_<date>.csv   (pop_inv_ps)
#        • Grain: species × bay × period × season
#        • Columns:
#            major_area | period | season | species_code
#            mu (mean abundance), var (variance), I_i (μ²/σ²)
#
#   2) community_invariability_<date>.csv   (comm_inv_ps)
#        • Grain: bay × period × season
#        • Columns:
#            major_area | period | season
#            mu_C (mean total abundance), var_C (variance), I_C (μ²/σ²)
#
#   3) synchrony_ps_<date>.csv   (synchrony_ps)
#        • Grain: bay × period × season
#        • Columns:
#            sum_sd_i (Σ species SDs), var_C (community variance), φ_LM
#            φ_LM = σ_C² / (Σ σ_i)²  ∈ [0,1]
#
#   4) portfolio_ps_<date>.csv   (portfolio_fix2)
#        • Grain: bay × period × season
#        • Columns:
#            major_area | period | season
#            weighted_mean_I_i, I_C, portfolio (= I_C / mean_I_i(weighted))
#        • Interpretation:
#            > 1 → diversity buffers variability (portfolio effect)
#            < 1 → diversity amplifies variability (no buffering)
#
#  Notes
#  -----
#   • Periods: 11 bins (1992–1994 … 2022–2024), each 3 years.
#   • Seasons: Winter (Dec–Feb), Spring (Mar–May), Summer (Jun–Aug), Fall (Sep–Nov).
#   • Population invariability I_i and community invariability I_C are based
#     on **monthly temporal series** of abundances within each stratum.
#   • Synchrony and portfolio effect are computed from the same strata.
#   • Missing months shorten the time series; mean/variance estimates remain
#     unbiased under MCAR missingness.
#
#  Author:  Masami Fujiwara with assistance of ChatGPT 5.0 in debugging 
#     Most of the annotations were added by ChatGPT for readability. 
#  Date:    2026-08-24   
###############################################################################


# ── 0. Packages & Initialization───────────────────────────────────────────────
library(tidyverse)         # dplyr, tidyr, readr, ggplot2 …

rm(list=ls())  # Erase all variables created previously

# retrieves the file path of the script file.
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

# ── Load data and keep only 'station`, `data`, and 'species' ──────────────────
load("data.Rdata")   # provides STATION_BS, bagseine (names from your .Rdata)

station <- STATION_BS %>% 
  dplyr::select(sample_id, major_area, year, month) %>% 
  filter(major_area %in% 1:8) # These are the eight major bays. 

data <- bagseine %>% 
  dplyr::select(sample_id, major_area, year, month, catch, species_code)

rm(list = setdiff(ls(), c("station", "data","species","script_dir")))  # keep only what we need

## 1 ── Remove ultra-rare species ------------------------------------------------
rare_sp <- data %>% 
  count(species_code, wt = (catch > 0)) %>%   # TRUE/FALSE summed → n pos obs
  filter(n <= 5) %>% 
  pull(species_code)

data_trim <- data %>% 
  filter(!species_code %in% rare_sp)

## 2 ── Sample ✕ species matrix with zeros --------------------------------------
sample_species <- data_trim %>% 
  group_by(sample_id, species_code) %>% 
  summarise(abund = sum(catch), .groups = "drop")   # total per haul & species

species_wide <- sample_species %>% 
  pivot_wider(names_from  = species_code,
              values_from = abund,
              values_fill = 0,
              names_prefix = "sp_")                 # sp_101, sp_207, …

final_df_wide <- station %>%                       # every sampling occasion
  distinct(sample_id, major_area, year, month) %>% 
  left_join(species_wide, by = "sample_id") %>% 
  mutate(across(starts_with("sp_"), \(x) replace_na(x, 0)))

# ── Period scheme: 11 periods, all 3-year bins over 1992–2024 ────────────────
start_year <- 1992
end_year   <- 2024
period_lengths <- rep(3, 11)  # P1..P11, each 3 years
stopifnot(sum(period_lengths) == (end_year - start_year + 1))

# Build a year → period lookup
year_period_lut <- tibble(
  year       = seq(start_year, end_year),
  period_num = rep(seq_along(period_lengths), times = period_lengths)
) %>%
  mutate(
    period = factor(paste0("P", period_num), levels = paste0("P", 1:11)),
    period_is_partial = FALSE
  )

## 3 ── Long data with numeric species_code -------------------------------------
monthly_sp <- final_df_wide %>% 
  # keep only years covered by the scheme
  filter(year >= start_year, year <= end_year) %>% 
  # attach period mapping
  left_join(year_period_lut, by = "year") %>% 
  pivot_longer(
    cols            = starts_with("sp_"),
    names_to        = "species_code",
    names_pattern   = "sp_(.*)",
    names_transform = list(species_code = as.integer),
    values_to       = "abund"
  ) %>% 
  mutate(
    species_code = as.integer(species_code),
    date  = as.Date(sprintf("%d-%02d-15", year, month)),
    
    # Seasons for the Texas coast
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",   # DJF
      month %in% 3:5         ~ "Spring",   # MAM
      month %in% 6:8         ~ "Summer",   # JJA
      TRUE                   ~ "Fall"      # SON
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Fall"))
  ) %>% 
  # safety
  filter(!is.na(period))


# ── Population invariability  I_i  ──────────────────────────────────────────
pop_inv_ps <- monthly_sp %>% 
  group_by(major_area, period, season, species_code) %>% 
  summarise(
    mu  = mean(abund),
    var = var(abund),
    I_i = if_else(var > 0, mu^2 / var, NA_real_),
    .groups = "drop"
  )

# ── Community invariability  I_C  ──────────────────────────────────────────
comm_inv_ps <- monthly_sp %>% 
  group_by(major_area, period, season, date) %>%            # monthly totals
  summarise(total_abund = sum(abund), .groups = "drop") %>% 
  group_by(major_area, period, season) %>% 
  summarise(
    mu_C  = mean(total_abund),
    var_C = var(total_abund),
    I_C   = if_else(var_C > 0, mu_C^2 / var_C, NA_real_),
    .groups = "drop"
  )

# ── Synchrony φ_LM (Loreau–de Mazancourt 2008) ───────────────────────────────
# Compute on aligned monthly series (same dates, same species set, same units)

# 1) Monthly totals per species (temporal series; align units for all)
sp_monthly <- monthly_sp %>%
  group_by(major_area, period, season, date, species_code) %>%
  summarise(abund_date = sum(abund), .groups = "drop")

# ---- choose the species set used in synchrony ----
# Option A: all species with finite temporal variance in the stratum
sp_var_temporal <- sp_monthly %>%
  group_by(major_area, period, season, species_code) %>%
  summarise(var_i = var(abund_date, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(var_i))   # keep only species with definable temporal variance

# (Option B: use your prevalence filter or abundance-weighted set;
#  just replace sp_var_temporal with a filtered subset before proceeding.)

# 2) Rebuild the community monthly series from the SAME species set
sp_set <- sp_var_temporal %>%
  dplyr::select(major_area, period, season, species_code)

comm_monthly_kept <- sp_monthly %>%
  inner_join(sp_set, by = c("major_area","period","season","species_code")) %>%
  group_by(major_area, period, season, date) %>%
  summarise(total_abund_kept = sum(abund_date), .groups = "drop")

# 3) Compute per-stratum pieces on the SAME dates
var_parts <- sp_var_temporal %>%
  group_by(major_area, period, season) %>%
  summarise(
    sum_var_i = sum(var_i, na.rm = TRUE),
    sum_sd_i  = sum(sqrt(pmax(var_i,0)), na.rm = TRUE),
    .groups = "drop"
  )

comm_var_kept <- comm_monthly_kept %>%
  group_by(major_area, period, season) %>%
  summarise(var_C_kept = var(total_abund_kept, na.rm = TRUE),
            n_dates = n(), .groups = "drop")

# 4) Synchrony: variance-ratio and bounded 0–1 variant
synchrony_fixed <- var_parts %>%
  inner_join(comm_var_kept, by = c("major_area","period","season")) %>%
  mutate(
    phi_var = if_else(sum_var_i > 0, var_C_kept / sum_var_i, NA_real_),
    phi  = if_else(sum_sd_i  > 0, var_C_kept / (sum_sd_i^2), NA_real_)
  )

# ── Portfolio effect  I_C / mean(I_i) ──────────────────────────────────────

# ==========================================================
# Portfolio: Temporal + Abundance-weighted mean(I_i)
# ==========================================================

pop_inv_ps_temporal <- monthly_sp %>%
  group_by(major_area, period, season, date, species_code) %>%
  summarise(abund_date = sum(abund), .groups = "drop") %>%
  group_by(major_area, period, season, species_code) %>%
  summarise(
    mu  = mean(abund_date, na.rm = TRUE),
    var = var(abund_date,  na.rm = TRUE),
    I_i = if_else(var > 0, mu^2 / var, NA_real_),
    .groups = "drop"
  )

weights_temporal <- monthly_sp %>%
  group_by(major_area, period, season, date, species_code) %>%
  summarise(sum_abund = sum(abund), .groups = "drop") %>%
  group_by(major_area, period, season, species_code) %>%
  summarise(w_raw = mean(sum_abund, na.rm = TRUE), .groups = "drop") %>%
  group_by(major_area, period, season) %>%
  mutate(sum_w = sum(w_raw, na.rm = TRUE),
         w = if_else(sum_w > 0, w_raw / sum_w, 0.0)) %>%
  ungroup() %>%
  select(major_area, period, season, species_code, w)

Ii_temporal_weighted <- pop_inv_ps_temporal %>%
  left_join(weights_temporal, by = c("major_area","period","season","species_code")) %>%
  group_by(major_area, period, season) %>%
  summarise(
    weighted_mean_I_i = {
      idx <- is.finite(I_i) & is.finite(w) & w > 0
      if (!any(idx)) NA_real_ else weighted.mean(I_i[idx], w[idx])
    },
    n_I_i_used = sum(is.finite(I_i) & is.finite(w) & w > 0),
    .groups = "drop"
  )

portfolio_fix2 <- Ii_temporal_weighted %>%
  left_join(comm_inv_ps %>% select(major_area, period, season, I_C),
            by = c("major_area","period","season")) %>%
  mutate(portfolio = if_else(is.finite(weighted_mean_I_i) & weighted_mean_I_i > 0,
                             I_C / weighted_mean_I_i, NA_real_)) %>%
  mutate(version = "Fix 2: Temporal + weighted")


# ── Strata counts: identify fully sampled strata ──────────────────────────────
station_strat <- station %>%
  filter(year >= start_year, year <= end_year) %>%
  left_join(year_period_lut %>% select(year, period), by = "year") %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% 3:5         ~ "Spring",
      month %in% 6:8         ~ "Summer",
      TRUE                   ~ "Fall"
    ),
    season = factor(season, levels = c("Winter","Spring","Summer","Fall")),
    ym = sprintf("%d-%02d", year, month)
  ) %>%
  filter(!is.na(period))

# Count hauls per bay × period × season
strata_counts_bps <- station_strat %>%
  group_by(major_area, period, season) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    n_months  = n_distinct(ym),
    months_observed = paste(sort(unique(ym)), collapse = ", "),
    .groups = "drop"
  )

# Keep only strata with exactly 180 hauls
valid_strata <- strata_counts_bps %>%
  filter(n_samples == 180) %>%
  select(major_area, period, season)

# ── Apply filter to all results ───────────────────────────────────────────────
pop_inv_ps <- pop_inv_ps %>%
  semi_join(valid_strata, by = c("major_area","period","season"))

comm_inv_ps <- comm_inv_ps %>%
  semi_join(valid_strata, by = c("major_area","period","season"))

synchrony_ps <- synchrony_fixed %>%
  semi_join(valid_strata, by = c("major_area","period","season"))

portfolio_ps <- portfolio_fix2 %>%
  semi_join(valid_strata, by = c("major_area","period","season"))

results_dir <- file.path(script_dir, "results")
today_tag <- format(Sys.Date(), "%Y-%m-%d")

write_csv(
  pop_inv_ps,
  file = file.path(results_dir, paste0("population_invariability_", today_tag, ".csv"))
)

write_csv(
  comm_inv_ps,
  file = file.path(results_dir, paste0("community_invariability_", today_tag, ".csv"))
)

write_csv(
  synchrony_ps,
  file = file.path(results_dir, paste0("synchrony_ps_", today_tag, ".csv"))
)

write_csv(
  portfolio_ps,
  file = file.path(results_dir, paste0("portfolio_ps_", today_tag, ".csv"))
)
