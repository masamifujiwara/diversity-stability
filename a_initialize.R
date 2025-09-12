# ──────────────────────────────────────────────────────────────
#  Script:  a_initialize.R
#  Purpose: Initialize the analysis environment by
#           • Determining the directory in which this script resides
#           • Setting that directory as the working directory
#           • Creating a subdirectory named "results" for all output files
#
#  Usage:
#    • Source this script at the start of any analysis to ensure
#      that file paths are set relative to the script location.
#    • All results, figures, and intermediate data will be saved
#      in the automatically created "results" folder.
#
#  Inputs:  None (the script detects its own location)
#  Outputs: Creates a "results" subdirectory under the script directory
#
#  Author:  Masami Fujiwara
#  Date:    2025-09-12   # ← update with current date
# ──────────────────────────────────────────────────────────────

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

# --- Create "results" subdirectory if it does not exist ---
results_dir <- file.path(script_dir, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}