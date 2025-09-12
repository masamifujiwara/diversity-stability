# diversity-stability
# Evenness and Taylor’s law scaling shape biodiversity–stability relationships in subtropical estuarine communities

This repository contains all code and data required to reproduce the analyses for the manuscript  
**“Evenness and Taylor’s law scaling shape biodiversity–stability relationships in subtropical estuarine communities”** (in preparation for submission to *Ecology*).

## Repository Contents

The repository includes one data file and a series of R scripts.  
Scripts are named alphabetically to indicate the order in which they should be run.

### Data File

**data.RData** – R workspace file containing three tibble data frames:

* **bagseine** (`493,954 × 6`)  
  - `sample_id` (num): Unique identifier for each sampling event  
  - `major_area` (num): Texas coastal bay/estuary code  
  - `year` (num): Year of sampling  
  - `month` (num): Month of sampling (1–12)  
  - `species_code` (num): Species identifier (matches `species$species_code`)  
  - `catch` (num): Number of individuals captured  

* **species** (`740 × 4`)  
  - `species_code` (num): Species identifier  
  - `sci_name` (chr): Scientific name  
  - `com_name` (chr): Common name  
  - `taxa` (num): Broad taxonomic group code  

* **STATION_BS** (`73,600 × 4`)  
  - `sample_id` (num): Identifier matching `bagseine$sample_id`  
  - `major_area` (num): Texas coastal bay/estuary code  
  - `year` (num): Year of sampling  
  - `month` (num): Month of sampling (1–12)  

These data represent long-term fish and invertebrate community monitoring in Texas estuarine systems.

### File Descriptions

- **data.RData** – All processed data used in the analyses.
- **a_initialize.R** – *Initialize the analysis environment.*  
  Sets the working directory to the script’s location and creates a `results/` subfolder where all outputs will be saved.  
  Key features:
  - Determines the directory in which the script resides  
  - Sets that directory as the working directory  
  - Creates a subdirectory named `results` if it does not exist
- **b_invariability_metrics_PS.R** – Calculates population and community invariability metrics following Wang & Loreau (2016).
- **c_diversity_metrics.R** – Computes species diversity metrics (e.g., rarefied richness, Shannon diversity).
- **d_bray_curtis_diversity.R** – Computes Bray–Curtis dissimilarities among samples.
- **e1_regional_gamma.R** – Estimates regional (γ) diversity.
- **e2_community_stability_diversity.R** – Examines the relationship between diversity and community stability.
- **e3_portfolio.R** – Evaluates portfolio effects of diversity on stability.
- **e4_synchrony.R** – Quantifies species synchrony/asynchrony.
- **e5_turnover.R** – Calculates species turnover rates.
- **e6_population_level.R** – Performs population-level stability analyses.
- **e7_taylor.R** – Fits Taylor’s law relationships between mean abundance and variance.
- **f5_turnover.R** – Additional turnover analyses and figure generation.

## Reproducibility

1. **Clone or download** this repository.
2. Open `a_initialize.R` in R or RStudio and source it.  
   - This will set the working directory to the script location and create a `results/` subfolder.
3. Run the remaining scripts in alphabetical order (as listed above).  
   - Each script saves its outputs into the `results/` directory.

All scripts are written for R ≥ 4.2 and use only publicly available packages from CRAN.

## Citation

If you use these data or scripts, please cite:

> Fujiwara, M. et al. *Evenness and Taylor’s law scaling shape biodiversity–stability relationships in subtropical estuarine communities.* *Ecology*. In preparation.

## License

This project is released under the **MIT License**.  
You are free to use, modify, and distribute the code and data with appropriate attribution.

---

For questions, please contact:  
**Masami Fujiwara**  
Department of Ecology & Conservation Biology  
Texas A&M University  
[add email or ORCID if desired]
