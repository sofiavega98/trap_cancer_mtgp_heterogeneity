################################################################################
# Configuration File for Lymphoma Subsets Analysis
# 
# This file defines all base paths and configuration parameters
# that are used across all analysis scripts.
#
# USAGE: 
# 1. Source this file at the beginning of each script: source("config.R")
# 2. Use the defined variables instead of hardcoded paths
# 3. Modify BASE_DIR below to match your local directory structure
#
# Author: Sofia Vega
# Date: 2024
################################################################################

# Set the base directory for the project
# This should point to the root of the cloned repository
# Users should modify this path to match their local directory structure
BASE_DIR <- getwd()  # Current working directory (recommended)
# Alternative: Set manually if needed
# BASE_DIR <- "~/path/to/your/Lymphoma_Subsets/github"

# Define subdirectories relative to BASE_DIR
DATA_DIR <- file.path(BASE_DIR, "data")
RESULTS_DIR <- file.path(BASE_DIR, "results")
CODE_DIR <- file.path(BASE_DIR, "code")
SA_DIR <- file.path(BASE_DIR, "SA")

# Define specific subdirectories
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")

# Define analysis parameters
YEARS <- c(1988:2003)
TREATED_YEAR <- 1995
DATE <- format(Sys.Date(), "%m%d%Y")

# Create output directories if they don't exist
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# Print configuration for verification
cat("Configuration loaded:\n")
cat("Base directory:", BASE_DIR, "\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Results directory:", RESULTS_DIR, "\n")
cat("Code directory:", CODE_DIR, "\n")
cat("SA directory:", SA_DIR, "\n")
cat("Study period:", min(YEARS), "-", max(YEARS), "\n")
cat("Treatment year:", TREATED_YEAR, "\n")
cat("Analysis date:", DATE, "\n\n") 
