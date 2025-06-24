################################################################################
source("config.R")

# Dependency Check Script
# 
# This script checks if all required packages are installed and provides
# installation instructions if needed.
#
# Run this script before running the main analysis to ensure all dependencies
# are properly installed.
################################################################################

cat("=== LYMPHOMA SUBSETS ANALYSIS - DEPENDENCY CHECK ===\n\n")

# List of required packages
required_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "ggplot2",      # Plotting (included in tidyverse but listed separately)
  "readr",        # Data reading (included in tidyverse)
  "gsynth",       # Generalized synthetic control
  "rstan",        # Stan interface for Bayesian modeling
  "kableExtra",   # Table formatting
  "cowplot",      # Plot arrangement
  "posterior"     # Posterior analysis
)

# Check which packages are installed
installed_packages <- installed.packages()[,"Package"]
missing_packages <- required_packages[!required_packages %in% installed_packages]

if (length(missing_packages) == 0) {
  cat("✓ All required packages are installed!\n")
  cat("You can now run the main analysis using:\n")
  cat("  source('run_all_analyses.R')\n\n")
} else {
  cat("✗ The following packages are missing:\n")
  for (pkg in missing_packages) {
    cat("  -", pkg, "\n")
  }
  cat("\n")
  
  cat("To install missing packages, run:\n")
  cat("install.packages(c(", paste0("'", missing_packages, "'", collapse = ", "), "))\n\n")
  
  # Special instructions for rstan
  if ("rstan" %in% missing_packages) {
    cat("Note: rstan may require additional setup. If installation fails, try:\n")
    cat("1. Install C++ compiler (Xcode on Mac, Rtools on Windows)\n")
    cat("2. Restart R session\n")
    cat("3. Run: install.packages('rstan')\n\n")
  }
}

# Check R version
r_version <- R.version.string
cat("R version:", r_version, "\n")

# Check if Stan is working (if rstan is installed)
if ("rstan" %in% installed_packages) {
  cat("\nTesting Stan installation...\n")
  tryCatch({
    library(rstan)
    cat("✓ Stan is working correctly\n")
  }, error = function(e) {
    cat("✗ Stan installation may have issues. Error:", e$message, "\n")
    cat("Try restarting R and running: library(rstan)\n")
  })
}

cat("\n=== DEPENDENCY CHECK COMPLETE ===\n") 
