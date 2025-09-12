###############################################################################
#  Observed richness & Shannon diversity (no rarefaction)
#  ───────────────────────────────────────────────────────
#  INPUT  – monthly_sp : species_code, major_area, year, month, sample_id,
#                        abund, date, period, season (zeros included)
#
#  OUTPUT – three tidy tables:
#    1) alpha_div : bay × period × season
#         major_area | period | season |
#         richness | shannon | hill_shannon
#    2) gamma_div : period × season (pooled across retained bays)
#         period | season |
#         richness | shannon | hill_shannon
#    3) beta_div  : multiplicative β
#         period | season |
#         beta_richness | beta_shannon
#       where beta_shannon uses Hill numbers: ^1D = exp(H')
###############################################################################

library(tidyverse)
library(lubridate)
library(vegan)

# Start with a clean workspace
rm(list = ls())

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
  filter(major_area %in% 1:8)

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

# ---- Identify fully sampled strata (n_samples == 180) ------------------------
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
    season = factor(season, levels = c("Winter","Spring","Summer","Fall"))
  ) %>%
  filter(!is.na(period))

strata_counts_bps <- station_strat %>%
  group_by(major_area, period, season) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop")

valid_strata <- strata_counts_bps %>%
  filter(n_samples == 180) %>%
  select(major_area, period, season)

message("Strata retained (n=180): ", nrow(valid_strata),
        " | Strata dropped: ",
        nrow(strata_counts_bps) - nrow(valid_strata))

# ---- Restrict monthly data to fully sampled strata ---------------------------
monthly_sp_valid <- monthly_sp %>%
  semi_join(valid_strata, by = c("major_area","period","season"))

###############################################################################
# 4.  Compute diversity metrics (observed only; no rarefaction)
###############################################################################

# Helper to build per-stratum community vectors (species totals)
build_comm <- function(df) tapply(df$abund, df$species_code, sum)

# -------------------------
# α-diversity (bay × P × S)
# -------------------------
alpha_div <- monthly_sp_valid %>%
  group_by(period, season, major_area) %>%
  summarise(comm = list(build_comm(across(everything()))), .groups = "drop") %>%
  mutate(
    richness     = purrr::map_dbl(comm, vegan::specnumber),
    shannon      = purrr::map_dbl(comm, ~ vegan::diversity(.x, index = "shannon")),
    hill_shannon = exp(shannon)  # ^1D
  ) %>%
  select(period, season, major_area, richness, shannon, hill_shannon)

# ------------------------------
# γ-diversity (pooled across bays)
# ------------------------------
gamma_div <- monthly_sp_valid %>%
  group_by(period, season) %>%
  summarise(comm = list(build_comm(across(everything()))), .groups = "drop") %>%
  mutate(
    richness     = purrr::map_dbl(comm, vegan::specnumber),
    shannon      = purrr::map_dbl(comm, ~ vegan::diversity(.x, index = "shannon")),
    hill_shannon = exp(shannon)  # ^1D
  ) %>%
  select(period, season, richness, shannon, hill_shannon)

# -------------------------------
# β-diversity (multiplicative β)
# -------------------------------
beta_div <- alpha_div %>%
  group_by(period, season) %>%
  summarise(
    mean_alpha_rich = mean(richness,     na.rm = TRUE),
    mean_alpha_hill = mean(hill_shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    gamma_div %>% rename(
      richness_gamma     = richness,
      hill_shannon_gamma = hill_shannon
    ),
    by = c("period","season")
  ) %>%
  transmute(
    period, season,
    beta_richness = richness_gamma     / mean_alpha_rich,
    beta_shannon  = hill_shannon_gamma / mean_alpha_hill  # uses ^1D
  )

###############################################################################
# 5.  Export results (observed-only names)
###############################################################################
results_dir <- file.path(script_dir, "results")
if (!dir.exists(results_dir)) dir.create(results_dir)
today_tag <- format(Sys.Date(), "%Y-%m-%d")

write_csv(
  alpha_div,
  file.path(results_dir, paste0("L5_alpha_diversity_observed_", today_tag, ".csv"))
)
write_csv(
  gamma_div,
  file.path(results_dir, paste0("L5_gamma_diversity_observed_", today_tag, ".csv"))
)
write_csv(
  beta_div,
  file.path(results_dir, paste0("L5_beta_diversity_observed_", today_tag, ".csv"))
)
