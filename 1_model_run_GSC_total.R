################################################################################
# Model Fitting: GSC (Interactive Fixed Effects) for Total Population
# 
# This script fits Interactive Fixed Effects (GSC) models using the gsynth package
# for the total population analysis. It analyzes treatment effects on lymphoma
# incidence rates for the overall population (all demographic groups combined).
#
# INPUTS:
# - Cleaned data: "data/total_Lymphoma_1988_2003_0-29.RData"
#
# OUTPUTS:
# - GSC model fit: "results/[date]/total_Lymphoma_1988_2003_0-29_GSC.RData"
#
# MODEL SPECIFICATION:
# - Method: Interactive Fixed Effects (GSC) via gsynth package
# - Outcome: Overall lymphoma incidence rates per 100,000
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Fixed effects: Two-way (county and year)
# - Inference: Nonparametric bootstrap (1000 iterations)
# - Cross-validation: Yes (for factor number selection)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Population: All ages 0-29 combined
# - Geographic scope: SEER registry counties
################################################################################

# Load configuration
source("config.R")

# Load required libraries
library(gsynth)
library(tidyverse)

# Define study parameters
years <- YEARS
treated_year <- TREATED_YEAR
date <- DATE

cat("Fitting GSC model for total population...\n")
cat("Study period:", min(years), "-", max(years), "\n")
cat("Treatment year:", treated_year, "\n")

# Load cleaned data
load(file.path(DATA_DIR, paste0("total_Lymphoma_",years[1],"_",years[length(years)],"_0-29.RData")))

# Calculate incidence rates per 100,000 population
data_full_pop$rate <- 100000 * (data_full_pop$CL_CASES/data_full_pop$POP)

# Fit GSC (Interactive Fixed Effects) model
# This is the primary model used in main analysis
fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "ife", 
                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                  se = TRUE, k=2, nboots = 1000, parallel = FALSE, inference = "nonparametric")

# Create results directory if it doesn't exist
dir.create(file.path(RESULTS_DIR, date), showWarnings = FALSE, recursive = TRUE)

# Save GSC model results
save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('total_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))

cat("✓ GSC model for total population completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")
  
  




