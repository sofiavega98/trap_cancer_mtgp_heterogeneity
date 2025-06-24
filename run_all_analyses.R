#!/usr/bin/env Rscript
################################################################################
# Master Script: Lymphoma Subsets Analysis
# 
# This script reproduces all analyses from the paper:
# "MTGP_Lymphoma_Strata_paper.pdf"
#
# Author: Sofia Vega
# Date: 2025
#
# This script runs the complete analysis pipeline:
# 0. Data quality checks and validation
# 1. Data cleaning for all demographic subsets
# 2. Descriptive statistics and visualizations
# 3. GSC (Generalized Synthetic Control) model fitting
# 4. MTGP (Multi-Task Gaussian Process) model fitting  
# 5. Results export and table/figure generation
# 6. Sensitivity analysis with extended study period (1983-2003)
#
# OUTPUTS:
# - Tables: Main results tables for the paper
# - Figures: Treatment effect plots and diagnostics
# - Model fits: Saved Stan and GSC model objects
# - Sensitivity analysis results with extended study period
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
library(ggplot2)
library(readr)
library(gsynth)
library(rstan)
library(kableExtra)
library(cowplot)
library(posterior)

################################################################################
# LOAD CONFIGURATION
# This loads all base paths and parameters
################################################################################

# Load configuration file
source("config.R")

cat("=== LYMPHOMA SUBSETS ANALYSIS ===\n")
cat("Analysis date:", DATE, "\n")
cat("Study period:", min(YEARS), "-", max(YEARS), "\n")
cat("Treatment year:", TREATED_YEAR, "\n")
cat("Base directory:", BASE_DIR, "\n\n")

################################################################################
# STEP 0: DATA QUALITY CHECKS
# This step validates data consistency across all demographic subsets
################################################################################

cat("STEP 0: Running data quality checks and validation...\n")

# Run comprehensive data quality checks
# This validates that all demographic breakdowns are consistent and complete
source("data_checks.R")

cat("✓ Data quality checks completed\n\n")

################################################################################
# STEP 1: DATA CLEANING
# This step prepares the SEER lymphoma data for analysis by demographic subsets
################################################################################

cat("STEP 1: Cleaning data for all demographic subsets...\n")

# Clean data for each demographic subset
# Each script creates a cleaned dataset for a specific demographic group
# All scripts use configuration-based paths for portability
source("0_cleaning_data_total.R")      # Overall population
source("0_cleaning_data_age.R")        # Age groups (0-19, 20-29)
source("0_cleaning_data_sex.R")        # Sex (Male, Female)  
source("0_cleaning_data_race.R")       # Race (White, Black, Other)
source("0_cleaning_data_hisp.R")       # Hispanic status (Hispanic, Non-Hispanic)
source("0_cleaning_data_type.R")       # Lymphoma type (Hodgkin, Non-Hodgkin)

cat("✓ Data cleaning completed\n\n")

################################################################################
# STEP 2: DESCRIPTIVE STATISTICS
# This step creates descriptive statistics and visualizations
################################################################################

cat("STEP 2: Creating descriptive statistics and visualizations...\n")

# Create descriptive statistics and map
# This script generates Table 1 and Figure S2 for the paper
source("0_create_descriptive_statistics.R")

cat("✓ Descriptive statistics completed\n\n")

################################################################################
# STEP 3: GSC MODEL FITTING
# This step fits Generalized Synthetic Control models for each subset
################################################################################

cat("STEP 3: Fitting GSC models...\n")

# Fit GSC models for each demographic subset
# These models provide baseline treatment effect estimates
source("1_model_run_GSC_total.R")      # Overall population
source("1_model_run_GSC_age.R")        # Age groups
source("1_model_run_GSC_sex.R")        # Sex groups
source("1_model_run_GSC_race.R")       # Race groups
source("1_model_run_GSC_hisp.R")       # Hispanic status groups
source("1_model_run_GSC_type.R")       # Lymphoma type groups

cat("✓ GSC model fitting completed\n\n")

################################################################################
# STEP 4: MTGP MODEL FITTING  
# This step fits Multi-Task Gaussian Process models for each subset
################################################################################

cat("STEP 4: Fitting MTGP models...\n")

# Fit MTGP models for each demographic subset
# These models provide the main treatment effect estimates for the paper
source("2_model_run_MTGP_total.R")     # Overall population
source("2_model_run_MTGP_age.R")       # Age groups
source("2_model_run_MTGP_sex.R")       # Sex groups
source("2_model_run_MTGP_race.R")      # Race groups
source("2_model_run_MTGP_hisp.R")      # Hispanic status groups
source("2_model_run_MTGP_type.R")      # Lymphoma type groups

cat("✓ MTGP model fitting completed\n\n")

################################################################################
# STEP 5: RESULTS EXPORT
# This step generates all tables and figures for the paper
################################################################################

cat("STEP 5: Generating results tables and figures...\n")

# Generate all tables and figures
# This script creates the main results presented in the paper
source("3_export_results.R")

cat("✓ Results export completed\n\n")

################################################################################
# STEP 6: SENSITIVITY ANALYSIS (OPTIONAL)
# This step runs additional sensitivity analyses with extended study period
################################################################################

cat("STEP 6: Running sensitivity analyses...\n")

# Change to SA directory
setwd(SA_DIR)

# Run complete sensitivity analysis pipeline
# Note: Sensitivity analysis uses extended study period (1983-2003) and earlier treatment year (1990)

cat("  - Cleaning sensitivity analysis data...\n")
# Clean data for sensitivity analysis (extended period)
source("0_cleaning_data_total.R")      # Overall population
source("0_cleaning_data_age.R")        # Age groups
source("0_cleaning_data_sex.R")        # Sex groups
source("0_cleaning_data_race.R")       # Race groups
source("0_cleaning_data_hisp.R")       # Hispanic status groups
source("0_cleaning_data_type.R")       # Lymphoma type groups

cat("  - Fitting sensitivity analysis models...\n")
# Fit MC models for sensitivity analysis
source("1_model_run_MC_total.R")       # Overall population
source("1_model_run_MC_age.R")         # Age groups
source("1_model_run_MC_sex.R")         # Sex groups
source("1_model_run_MC_race.R")        # Race groups
source("1_model_run_MC_hisp.R")        # Hispanic status groups
source("1_model_run_MC_type.R")        # Lymphoma type groups

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

cat("=== ANALYSIS COMPLETE ===\n")
cat("All results have been generated and saved to the 'results' directory.\n")
cat("Check the following files:\n")
cat("- results/tables/: Contains all LaTeX tables for the paper\n")
cat("- results/figures/: Contains all figures for the paper\n")
cat("- results/: Contains saved model objects\n")
cat("- SA/results/: Contains sensitivity analysis results\n\n")

cat("Key outputs generated:\n")
cat("- Data quality validation tables (case_counts.tex, case_counts_grouped.tex)\n")
cat("- Table 1: Descriptive statistics by demographic group and treatment status\n")
cat("- Figure S2: Map of treated vs control counties\n")
cat("- Table 2: Main MTGP treatment effect estimates\n")
cat("- Table 3: Treatment effect ratios\n") 
cat("- Table 4: Sum of treatment effects across treated units/time periods\n")
cat("- Figures: Treatment effect plots and diagnostics\n")
cat("- Sensitivity analysis: Extended study period (1983-2003) results\n\n")

cat("Configuration used:\n")
cat("- Base directory:", BASE_DIR, "\n")
cat("- Data directory:", DATA_DIR, "\n")
cat("- Results directory:", RESULTS_DIR, "\n")
cat("- Analysis date:", DATE, "\n\n") 
