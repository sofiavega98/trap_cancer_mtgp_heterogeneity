################################################################################
# Sensitivity Analysis Script
# 
# This script runs the sensitivity analysis portion of the lymphoma subsets study.
# It should be run after the main analysis is complete.
#
# The sensitivity analysis includes:
# 1. Extended study period (1983-2003) vs main analysis (1988-2003)
# 2. Earlier treatment year (1990) vs main analysis (1995)
# 3. Matrix Completion (MC) models as an alternative to MTGP
# 4. Additional robustness checks
# 5. Comparison of results across different model specifications
#
# INPUTS:
# - Raw SEER data (same as main analysis)
# - Functions from functions.R
#
# OUTPUTS:
# - Sensitivity analysis results and tables
# - Comparison with main analysis results
#
# CONFIGURATION:
# All paths and parameters are managed through config.R
# This ensures portability across different systems
################################################################################

# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(501)

# Load required libraries
library(tidyverse)
library(gsynth)
library(rstan)

################################################################################
# LOAD CONFIGURATION
# This loads all base paths and parameters
################################################################################

# Load configuration file
source("config.R")

# Create output directories if they don't exist
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

cat("=== LYMPHOMA SUBSETS SENSITIVITY ANALYSIS ===\n")
cat("Analysis date:", DATE, "\n")
cat("Main study period:", min(YEARS), "-", max(YEARS), "\n")
cat("Main treatment year:", TREATED_YEAR, "\n")
cat("Sensitivity analysis: Extended period (1983-2003), earlier treatment (1990)\n")
cat("Base directory:", BASE_DIR, "\n\n")

################################################################################
# STEP 1: Change to SA directory and run complete sensitivity analysis
################################################################################

cat("STEP 1: Running complete sensitivity analysis pipeline...\n")

# Change to SA directory
setwd(SA_DIR)

# Run complete sensitivity analysis pipeline
# Note: Sensitivity analysis uses extended study period (1983-2003) and earlier treatment year (1990)

cat("  - Cleaning sensitivity analysis data (extended period)...\n")
# Clean data for sensitivity analysis (extended period)
source("0_cleaning_data_total.R")      # Overall population
source("0_cleaning_data_age.R")        # Age groups
source("0_cleaning_data_sex.R")        # Sex groups
source("0_cleaning_data_race.R")       # Race groups
source("0_cleaning_data_hisp.R")       # Hispanic status groups
source("0_cleaning_data_type.R")       # Lymphoma type groups

cat("  - Fitting MC models for sensitivity analysis...\n")
# Fit MC models for sensitivity analysis
source("1_model_run_MC_total.R")       # Overall population
source("1_model_run_MC_age.R")         # Age groups
source("1_model_run_MC_sex.R")         # Sex groups
source("1_model_run_MC_race.R")        # Race groups
source("1_model_run_MC_hisp.R")        # Hispanic status groups
source("1_model_run_MC_type.R")        # Lymphoma type groups

cat("  - Fitting MTGP models for sensitivity analysis...\n")
# Fit MTGP models for sensitivity analysis
source("2_model_run_MTGP_total.R")     # Overall population
source("2_model_run_MTGP_age.R")       # Age groups
source("2_model_run_MTGP_sex.R")       # Sex groups
source("2_model_run_MTGP_race.R")      # Race groups
source("2_model_run_MTGP_hisp.R")      # Hispanic status groups
source("2_model_run_MTGP_type.R")      # Lymphoma type groups

cat("  - Exporting sensitivity analysis results...\n")
# Export sensitivity analysis results
source("3_export_results.R")

# Return to main directory
setwd(BASE_DIR)

cat("✓ Sensitivity analysis completed\n\n")

################################################################################
# COMPLETION
################################################################################

cat("=== SENSITIVITY ANALYSIS COMPLETE ===\n")
cat("Sensitivity analysis results have been generated and saved.\n")
cat("Check the SA/ directory for outputs.\n\n")

cat("Key sensitivity analysis components:\n")
cat("- Extended study period: 1983-2003 (vs main analysis: 1988-2003)\n")
cat("- Earlier treatment year: 1990 (vs main analysis: 1995)\n")
cat("- Matrix Completion (MC) models as alternative to MTGP\n")
cat("- Alternative model specifications\n")
cat("- Robustness checks\n")
cat("- Comparison with main MTGP results\n\n")

cat("Sensitivity analysis outputs:\n")
cat("- SA/results/: Contains all sensitivity analysis model objects\n")
cat("- SA/results/tables/: Contains sensitivity analysis tables\n")
cat("- SA/results/figures/: Contains sensitivity analysis figures\n")
cat("- Comparison tables with main analysis results\n\n")

cat("Configuration used:\n")
cat("- Base directory:", BASE_DIR, "\n")
cat("- SA directory:", SA_DIR, "\n")
cat("- Analysis date:", DATE, "\n\n")

cat("For questions about the sensitivity analysis, please refer to the README.md file.\n") 
