################################################################################
source("config.R")

# Model Fitting: GSC (Interactive Fixed Effects) for Age Groups
# 
# This script fits Interactive Fixed Effects (GSC) models using the gsynth package
# for age group analysis. It analyzes treatment effects on lymphoma incidence
# rates for young (0-19) and old (20-29) age groups separately.
#
# INPUTS:
# - Age-stratified data: "data/age_Lymphoma_1988_2003_0-29.RData"
#
# OUTPUTS:
# - Young age GSC results: "results/[date]/young_Lymphoma_1988_2003_0-29_GSC.RData"
# - Old age GSC results: "results/[date]/old_Lymphoma_1988_2003_0-29_GSC.RData"
#
# MODEL SPECIFICATION:
# - Method: Interactive Fixed Effects (GSC) via gsynth package
# - Outcome: Age-adjusted lymphoma incidence rates per 100,000
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Fixed effects: Two-way (county and year)
# - Inference: Nonparametric bootstrap (1000 iterations)
# - Cross-validation: Yes (for factor number selection)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Age groups: 0-19 years (young) and 20-29 years (old)
# - Geographic scope: SEER registry counties
################################################################################

## Run GSC in gsynth package 
## Feature: Age Groups (0-19 vs 20-29 years)

# Load required libraries
library(gsynth)
library(tidyverse)

# Define analysis date for file naming
date <- "09102024"

# Define study parameters
years <- c(1988:2003)          # Study period
treated_year <- 1995           # Year treatment was implemented

# Load age-stratified data
load(file.path(DATA_DIR, "age_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_full <- data_full_pop

# Calculate age-adjusted incidence rates per 100,000 population
data_full_pop_full$rate <- 100000 * (data_full_pop_full$CL_CASES/data_full_pop_full$POP)

# Split data by age groups for separate analysis
data_full_pop_young <- data_full_pop_full %>% filter(AGE == "Younger Age (0-4)")
data_full_pop_old <- data_full_pop_full %>% filter(AGE == "Older Age (5-6)")

# Save age-stratified datasets for model fitting
save(data_full_pop_young, file = file.path(DATA_DIR, paste0('young_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_pop_old, file = file.path(DATA_DIR, paste0('old_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit GSC models for each age group
for(i in c("Younger Age (0-4)","Older Age (5-6)")){  
  # Load appropriate age group data
  if(i == "Younger Age (0-4)"){
    load(file.path(DATA_DIR, paste0('young_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_young
  }else{
    load(file.path(DATA_DIR, paste0('old_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_old
  }
  
  # Note: Original MC model commented out for reference
  # Run MC Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_pop, estimator = "mc", 
   #                 index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
 # if(i == "Younger Age (0-4)"){
  #  # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('young_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
 # }else{
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('old_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}
  
  # Run GSC (Interactive Fixed Effects) in Gsynth
  # This is the primary model used in main analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "ife", 
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
cat("GSC models for age groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")



