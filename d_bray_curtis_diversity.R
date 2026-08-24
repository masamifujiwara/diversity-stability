###############################################################################
#  Script purpose
#  --------------
#  1.  Pre-process Texas bay bag-seine data
#        • drop ultra-rare species (≤ 5 non-zero records)
#        • build a complete sample × species matrix (zeros = true absences)
#        • convert to long format and add DJF–MAM–JJA–SON “season” factors
#
#  2.  Temporal β-diversity (Bray–Curtis) consistent with invariability metrics
#      ------------------------------------------------------------------------
#      • For each bay × period × season × month × year × species:
#          – compute **mean per-haul abundance** (total abundance ÷ hauls).
#      • Within each month, compute Bray–Curtis dissimilarity between
#        consecutive years (year t vs year t+1).
#      • Average Bray–Curtis values across the 3 months of the season
#        to obtain a seasonal turnover value for each year-pair.
#      • Within each 3-year period, average across the two consecutive
#        year-pairs (e.g., 1992–1993 and 1993–1994) to obtain a single
#        stratum-level Bray–Curtis dissimilarity.
#
#  Outputs
#  -------
#    1) beta_bray_pairs_monthly_<YYYY-MM-DD>.csv
#         • One row per bay × period × season × month × consecutive-year pair
#         • Columns: major_area | period | season | month | year1 | year2 | bray_consecutive
#
#    2) beta_bray_bay_period_season_<YYYY-MM-DD>.csv
#         • One row per bay × period × season
#         • Columns: major_area | period | season | bray_bc
#         • bray_bc = mean Bray–Curtis across months and consecutive-year pairs
#
#  Notes
#  -----
#    • bray values ∈ [0,1]: 0 = identical composition, 1 = no shared species.
#    • This approach preserves **within-season structure** (month-matched
#      comparisons) while integrating across years, making it conceptually
#      consistent with invariability measures.
#
#
#  Author:  Masami Fujiwara with assistance of ChatGPT 5.0 in debugging 
#     Most of the annotations were added by ChatGPT for readability. 
#  Date:    2026-08-24   
###############################################################################


###############################################################################
#  0.  Packages and initialisation
###############################################################################
library(tidyverse)
library(vegan)      # Bray–Curtis

rm(list = ls())  # clean workspace

# Resolve script dir
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

# ── Load data ────────────────────────────────────────────────────────────────
load("data.Rdata")   # provides STATION_BS, bagseine

station <- STATION_BS %>%
  dplyr::select(sample_id, major_area, year, month) %>%
  filter(major_area %in% 1:8)

data <- bagseine %>%
  dplyr::select(sample_id, major_area, year, month, catch, species_code)

rm(list = setdiff(ls(), c("station", "data", "species", "script_dir")))  # keep only what we need

## 1 ── Remove ultra-rare species ------------------------------------------------
rare_sp <- data %>%
  count(species_code, wt = (catch > 0)) %>%   # TRUE/FALSE summed → n pos obs
  filter(n <= 5) %>%
  pull(species_code)

data_trim <- data %>%
  filter(!species_code %in% rare_sp)

## 2 ── Sample × species matrix with zeros --------------------------------------
sample_species <- data_trim %>%
  group_by(sample_id, species_code) %>%
  summarise(abund = sum(catch), .groups = "drop")

species_wide <- sample_species %>%
  pivot_wider(names_from  = species_code,
              values_from = abund,
              values_fill = 0,
              names_prefix = "sp_")

final_df_wide <- station %>%
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

## 3 ── Long data with species_code numeric + seasons + PERIOD ------------------
monthly_sp <- final_df_wide %>%
  filter(year >= start_year, year <= end_year) %>%
  left_join(year_period_lut, by = "year") %>%   # <-- add PERIOD to monthly data
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
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",   # DJF
      month %in% 3:5         ~ "Spring",   # MAM
      month %in% 6:8         ~ "Summer",   # JJA
      month %in% 9:11         ~ "Fall"      # SON
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Fall"))
  ) %>%
  filter(!is.na(period))

## 4 ── Identify and KEEP only fully sampled strata (n_samples == 180) ----------
station_strat <- station %>%
  filter(year >= start_year, year <= end_year) %>%
  left_join(year_period_lut, by = "year") %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% 3:5         ~ "Spring",
      month %in% 6:8         ~ "Summer",
      month %in% 9:11          ~ "Fall"
    ),
    season = factor(season, levels = c("Winter","Spring","Summer","Fall"))
  ) %>%
  filter(!is.na(period))

valid_strata <- station_strat %>%
  group_by(major_area, period, season) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop") %>%
  filter(n_samples == 180) %>%
  select(major_area, period, season)

message(
  "Retaining ", nrow(valid_strata), " strata with 180 samples. Dropped ",
  nrow(station_strat %>% distinct(major_area, period, season)) - nrow(valid_strata),
  " strata with < 180."
)

# Apply filter to analysis data
monthly_sp <- monthly_sp %>%
  semi_join(valid_strata, by = c("major_area","period","season"))

###############################################################################
#  Temporal β-diversity (Bray–Curtis) consistent with invariability
#  ----------------------------------------------------------------
#  • For each bay × period × season × month:
#      – build a community vector of mean-per-haul abundance per species
#      – compute Bray–Curtis between consecutive years (t vs t+1)
#  • Within each bay × period × season:
#      – average Bray–Curtis across the 3 months in the season
#      – if 3-year period → also average across the 2 consecutive-year pairs
#  • Outputs:
#      1) beta_bray_pairs_monthly : per bay × period × season × month × (year1, year2)
#      2) beta_bray_stratum       : one row per bay × period × season (mean across months & pairs)
###############################################################################

# ---- Assumes you already have: station, data, monthly_sp with period & season ----
# start_year <- 1992; end_year <- 2024; year_period_lut <- ...
# monthly_sp columns used: major_area, period, season, year, month, sample_id, species_code, abund

# ============ (Optional but recommended) keep only fully sampled strata ============
# Count hauls in station by bay × period × season and keep n_samples == 180
station_strat <- station %>%
  filter(year >= start_year, year <= end_year) %>%
  left_join(year_period_lut %>% select(year, period), by = "year") %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% 3:5         ~ "Spring",
      month %in% 6:8         ~ "Summer",
      month %in% 9:11          ~ "Fall"
    ),
    season = factor(season, levels = c("Winter","Spring","Summer","Fall"))
  ) %>%
  filter(!is.na(period))

valid_strata <- station_strat %>%
  group_by(major_area, period, season) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop") %>%
  filter(n_samples == 180) %>%
  select(major_area, period, season)

# Filter monthly data to fully sampled strata
monthly_sp_valid <- monthly_sp %>%
  semi_join(valid_strata, by = c("major_area","period","season"))

# ===================== Build month-level mean-per-haul communities =================
# For each bay × period × season × month × year × species, compute mean per haul
mps <- monthly_sp_valid %>%
  group_by(major_area, period, season, month, year, species_code) %>%
  summarise(
    total_abund = sum(abund, na.rm = TRUE),
    n_hauls     = n_distinct(sample_id),
    .groups = "drop"
  ) %>%
  mutate(mean_per_haul = if_else(n_hauls > 0, total_abund / n_hauls, 0))

# Helper to compute consecutive-year Bray–Curtis for one (bay, period, season, month)
bray_consecutive_one_month <- function(df) {
  # df has columns: year, species_code, mean_per_haul
  if (n_distinct(df$year) < 2) {
    return(tibble(year1 = integer(), year2 = integer(), bray = double()))
  }
  # Wide community matrix: rows = year, cols = species
  X <- df %>%
    select(year, species_code, mean_per_haul) %>%
    pivot_wider(names_from = species_code, values_from = mean_per_haul, values_fill = 0) %>%
    arrange(year)
  yrs <- X$year
  sp_mat <- as.matrix(X[ , setdiff(names(X), "year"), drop = FALSE])
  rownames(sp_mat) <- as.character(yrs)
  
  # Bray–Curtis on all pairs, then pick consecutive (t, t+1)
  bc <- as.matrix(vegan::vegdist(sp_mat, method = "bray"))
  pairs <- tibble(year1 = yrs, year2 = yrs + 1) %>%
    filter(year2 %in% yrs)
  
  if (nrow(pairs) == 0) {
    return(tibble(year1 = integer(), year2 = integer(), bray = double()))
  }
  
  pairs %>%
    mutate(bray = purrr::map2_dbl(year1, year2, ~ bc[as.character(.x), as.character(.y)]))
}

# --------------------- Month-by-month consecutive-year Bray–Curtis -----------------
beta_bray_pairs_monthly <- mps %>%
  group_by(major_area, period, season, month) %>%
  group_modify(~{
    df <- .x %>% select(year, species_code, mean_per_haul)
    br <- bray_consecutive_one_month(df)
    if (nrow(br) == 0) return(tibble(year1 = integer(), year2 = integer(), bray = double()))
    br
  }) %>%
  ungroup() %>%
  rename(bray_consecutive = bray) %>%
  arrange(major_area, period, season, month, year1, year2)

# ---------------------- Aggregate to the stratum (bay × period × season) ----------
# For each bay × period × season:
#   1) average across the 3 months (Dec/Jan/Feb etc.) for each consecutive pair
#   2) then average across the consecutive pairs inside the 3-year period
beta_bray_stratum <- beta_bray_pairs_monthly %>%
  group_by(major_area, period, season, year1, year2) %>%
  summarise(
    bray_mean_across_months = mean(bray_consecutive, na.rm = TRUE),
    n_months_used = sum(is.finite(bray_consecutive)),
    .groups = "drop"
  ) %>%
  group_by(major_area, period, season) %>%
  summarise(
    bray_bc = mean(bray_mean_across_months, na.rm = TRUE), # final stratum value
    n_pairs = n(),                                         # should be 2 inside a 3-yr period
    mean_months_per_pair = mean(n_months_used),
    .groups = "drop"
  ) %>%
  arrange(major_area, period, season)

# ----------------------------- Export -----------------------------------
results_dir <- file.path(script_dir, "results")
if (!dir.exists(results_dir)) dir.create(results_dir)
today_tag <- format(Sys.Date(), "%Y-%m-%d")

readr::write_csv(
  beta_bray_pairs_monthly,
  file.path(results_dir, paste0("beta_bray_pairs_monthly_", today_tag, ".csv"))
)
readr::write_csv(
  beta_bray_stratum,
  file.path(results_dir, paste0("beta_bray_bay_period_season_", today_tag, ".csv"))
)
