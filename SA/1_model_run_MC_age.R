################################################################################
source("config.R")

# Sensitivity Analysis: Model Fitting - MC (Matrix Completion) for Age Groups
# 
# This script fits Matrix Completion (MC) models using the gsynth package for
# age group sensitivity analysis. It analyzes treatment effects on lymphoma
# incidence rates for young (0-19) and old (20-29) age groups using an extended
# study period and earlier treatment year.
#
# INPUTS:
# - Age-stratified data: "~/Documents/Harvard/Lymphoma_Subsets/data/SA/age_Lymphoma_1983_2003_0-29.RData"
#
# OUTPUTS:
# - Young age MC results: "~/Documents/Harvard/Lymphoma_Subsets/results/09102024/SA/young_Lymphoma_1983_2003_0-29_GSC.RData"
# - Old age MC results: "~/Documents/Harvard/Lymphoma_Subsets/results/09102024/SA/old_Lymphoma_1983_2003_0-29_GSC.RData"
#
# MODEL SPECIFICATION:
# - Method: Matrix Completion (MC) via gsynth package
# - Outcome: Age-adjusted lymphoma incidence rates per 100,000
# - Treatment: Binary indicator for CA/CT counties after 1990
# - Fixed effects: Two-way (county and year)
# - Inference: Nonparametric bootstrap (1000 iterations)
# - Cross-validation: Yes (for factor number selection)
#
# SENSITIVITY ANALYSIS PARAMETERS:
# - Study period: 1983-2003 (21 years)
# - Treatment year: 1990 (vs 1995 in main analysis)
# - Age groups: 0-19 years (young) and 20-29 years (old)
# - Purpose: Test robustness under different time assumptions
################################################################################

## Run MC in gsynth package for age group sensitivity analysis

# Load required libraries
library(gsynth)
library(tidyverse)

# Define analysis date for file naming
date <- "09102024"

# Define sensitivity analysis parameters
years <- c(1983:2003)          # Extended study period (21 years)
treated_year <- 1990           # Earlier treatment year for sensitivity analysis

# Load age-stratified sensitivity analysis data
load(file.path(DATA_DIR, "SA/age_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_full <- data_full_hisp_noZero

# Calculate age-adjusted incidence rates per 100,000 population
data_full_hisp_full$rate <- 100000 * (data_full_hisp_full$CL_CASES/data_full_hisp_full$POP)

# Split data by age groups for separate analysis
data_full_hisp_young <- data_full_hisp_full %>% filter(AGE == "Younger Age (0-4)")
data_full_hisp_old <- data_full_hisp_full %>% filter(AGE == "Older Age (5-6)")

# Save age-stratified datasets for model fitting
save(data_full_hisp_young, file = file.path(DATA_DIR, paste0('young_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_hisp_old, file = file.path(DATA_DIR, paste0('old_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit MC models for each age group
for(i in c("Younger Age (0-4)","Older Age (5-6)")){  
  # Load appropriate age group data
  if(i == "Younger Age (0-4)"){
    load(file.path(DATA_DIR, paste0('young_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_young
  }else{
    load(file.path(DATA_DIR, paste0('old_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_old
  }
  
  # Note: Original MC model commented out for reference
  # Run MC Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_pop, estimator = "mc", 
  #                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
  #if(i == "Younger Age (0-4)"){
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('young_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}else{
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('old_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}
  
  # Run MC (Matrix Completion) in Gsynth
  # This is the primary model used in sensitivity analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "mc", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, k=2, nboots = 1000, parallel = FALSE, inference = "nonparametric")
  
  # Save results by age group
  if(i == "Younger Age (0-4)"){
    # Save young age group results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('young_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  }else{
    # Save old age group results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('old_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  } 
}

# Print completion message
cat("MC models for age groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")



