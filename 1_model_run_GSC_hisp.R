################################################################################
source("config.R")

# Model Fitting: GSC (Interactive Fixed Effects) for Hispanic Ethnicity
# 
# This script fits Interactive Fixed Effects (GSC) models using the gsynth package
# for Hispanic ethnicity analysis. It analyzes treatment effects on lymphoma
# incidence rates for Hispanic and Non-Hispanic groups separately.
#
# INPUTS:
# - Hispanic ethnicity data: "data/hispanic_Lymphoma_1988_2003_0-29.RData"
#
# OUTPUTS:
# - Non-Hispanic GSC results: "results/[date]/nonHisp_Lymphoma_1988_2003_0-29_GSC.RData"
# - Hispanic GSC results: "results/[date]/Hisp_Lymphoma_1988_2003_0-29_GSC.RData"
#
# MODEL SPECIFICATION:
# - Method: Interactive Fixed Effects (GSC) via gsynth package
# - Outcome: Hispanic ethnicity-adjusted lymphoma incidence rates per 100,000
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Fixed effects: Two-way (county and year)
# - Inference: Nonparametric bootstrap (1000 iterations)
# - Cross-validation: Yes (for factor number selection)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Ethnicity groups: Non-Hispanic (0) and Hispanic (1)
# - Geographic scope: SEER registry counties
################################################################################

## Run GSC in gsynth package 
## Feature: Hispanic Ethnicity (Hispanic vs Non-Hispanic)

# Load required libraries
library(gsynth)
library(tidyverse)

# Define analysis date for file naming
date <- "09102024"

# Define study parameters
years <- c(1988:2003)          # Study period
treated_year <- 1995           # Year treatment was implemented

# Load Hispanic ethnicity data
load(file.path(DATA_DIR, "hispanic_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_full <- data_full_pop_noZero

# Calculate Hispanic ethnicity-adjusted incidence rates per 100,000 population
data_full_pop_full$rate <- 100000 * (data_full_pop_full$CL_CASES/data_full_pop_full$POP)

# Split data by Hispanic ethnicity for separate analysis
data_full_pop_nonHisp <- data_full_pop_full %>% filter(HISPANIC == 0)
data_full_pop_Hisp <- data_full_pop_full %>% filter(HISPANIC == 1)

# Save Hispanic ethnicity-stratified datasets for model fitting
save(data_full_pop_nonHisp, file = file.path(DATA_DIR, paste0('nonHisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_pop_Hisp, file = file.path(DATA_DIR, paste0('Hisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit GSC models for each Hispanic ethnicity group
for(i in c(0,1)){
  # Load appropriate Hispanic ethnicity group data
  if(i == 0){
    load(file.path(DATA_DIR, paste0('nonHisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_nonHisp
  }else{
    load(file.path(DATA_DIR, paste0('Hisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_Hisp
  }  
  
  # Note: Original MC model commented out for reference
  # Run MC Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_pop, estimator = "mc", 
  #                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
  #if(i == 0){
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('nonHisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}else{
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('Hisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}
  
  # Run GSC (Interactive Fixed Effects) in Gsynth
  # This is the primary model used in main analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "ife", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, k=2, nboots = 1000, parallel = FALSE, inference = "nonparametric")
  
  # Save results by Hispanic ethnicity group
  if(i == 0){
    # Save Non-Hispanic results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('nonHisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  }else{
    # Save Hispanic results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('Hisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  } 
}

# Print completion message
cat("GSC models for Hispanic ethnicity groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")

