################################################################################
source("config.R")

# Model Fitting: GSC (Interactive Fixed Effects) for Race
# 
# This script fits Interactive Fixed Effects (GSC) models using the gsynth package
# for race-based analysis. It analyzes treatment effects on lymphoma incidence
# rates for White, Black, and Other racial groups separately.
#
# INPUTS:
# - Race-stratified data: "data/race_x3_Lymphoma_1988_2003_0-29.RData"
#
# OUTPUTS:
# - White GSC results: "results/[date]/white_Lymphoma_1988_2003_0-29_race_GSC.RData"
# - Black GSC results: "results/[date]/black_Lymphoma_1988_2003_0-29_race_GSC.RData"
# - Other GSC results: "results/[date]/other_Lymphoma_1988_2003_0-29_race_GSC.RData"
#
# MODEL SPECIFICATION:
# - Method: Interactive Fixed Effects (GSC) via gsynth package
# - Outcome: Race-adjusted lymphoma incidence rates per 100,000
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Fixed effects: Two-way (county and year)
# - Inference: Nonparametric bootstrap (1000 iterations)
# - Cross-validation: Yes (for factor number selection)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Race groups: White (1), Black (2), Other (3)
# - Geographic scope: SEER registry counties
################################################################################

## Run GSC in gsynth package 
## Feature: Race (White, Black, Other)

# Load required libraries
library(gsynth)
library(tidyverse)

# Define analysis date for file naming
date <- "09102024"

# Define study parameters
years <- c(1988:2003)          # Study period
treated_year <- 1995           # Year treatment was implemented

# Load race-stratified data
load(file.path(DATA_DIR, "race_x3_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_full <- data_full_pop

# Calculate race-adjusted incidence rates per 100,000 population
data_full_pop_full$rate <- 100000 * (data_full_pop_full$CL_CASES/data_full_pop_full$POP)

# Split data by race for separate analysis
data_full_pop_white <- data_full_pop_full %>% filter(RACE == 1)
data_full_pop_black <- data_full_pop_full %>% filter(RACE == 2)
data_full_pop_other <- data_full_pop_full %>% filter(RACE == 3)

# Save race-stratified datasets for model fitting
save(data_full_pop_white, file = file.path(DATA_DIR, paste0('white_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_pop_black, file = file.path(DATA_DIR, paste0('black_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_pop_other, file = file.path(DATA_DIR, paste0('other_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit GSC models for each race group
for(i in c(1,2,3)){  
  # Load appropriate race group data
  if(i == 1){
    load(file.path(DATA_DIR, paste0('white_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_white
  } else if ( i == 2) {
    load(file.path(DATA_DIR, paste0('black_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_black  
  }else{
    load(file.path(DATA_DIR, paste0('other_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_other
  }
  
  # Note: Original MC model commented out for reference
  # Run MC in Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_hisp, estimator = "mc", 
  #                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
  # Save results
  #if(i == 1){
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('white_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #} else if ( i == 2) {
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('black_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #} else {
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('other_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}
  
  # Run GSC (Interactive Fixed Effects) in Gsynth
  # This is the primary model used in main analysis
  fit_gsc <- gsynth(data = data_full_hisp, Y = "rate", D = "C", estimator = "ife", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, nboots = 1000, parallel = FALSE, inference = "nonparametric")
  
  # Save results by race group
  if(i == 1){
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('white_Lymphoma_',years[1],'_',years[length(years)],'_0-29_race_GSC.RData')))
  } else if ( i == 2) {
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('black_Lymphoma_',years[1],'_',years[length(years)],'_0-29_race_GSC.RData')))
  } else {
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('other_Lymphoma_',years[1],'_',years[length(years)],'_0-29_race_GSC.RData')))
  }
}

# Print completion message
cat("GSC models for race groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")
