################################################################################
source("config.R")

# Model Fitting: GSC (Interactive Fixed Effects) for Sex
# 
# This script fits Interactive Fixed Effects (GSC) models using the gsynth package
# for sex-based analysis. It analyzes treatment effects on lymphoma incidence
# rates for males and females separately.
#
# INPUTS:
# - Sex-stratified data: "data/sex_Lymphoma_1988_2003_0-29.RData"
#
# OUTPUTS:
# - Male GSC results: "results/[date]/male_Lymphoma_1988_2003_0-29_GSC.RData"
# - Female GSC results: "results/[date]/female_Lymphoma_1988_2003_0-29_GSC.RData"
#
# MODEL SPECIFICATION:
# - Method: Interactive Fixed Effects (GSC) via gsynth package
# - Outcome: Sex-adjusted lymphoma incidence rates per 100,000
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Fixed effects: Two-way (county and year)
# - Inference: Nonparametric bootstrap (1000 iterations)
# - Cross-validation: Yes (for factor number selection)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Sex groups: Male (1) and Female (2)
# - Geographic scope: SEER registry counties
################################################################################

## Run GSC in gsynth package 
## Feature: Sex (Male vs Female)

# Load required libraries
library(gsynth)
library(tidyverse)

# Define analysis date for file naming
date <- "09102024"

# Define study parameters
years <- c(1988:2003)          # Study period
treated_year <- 1995           # Year treatment was implemented

# Load sex-stratified data
load(file.path(DATA_DIR, "sex_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_full <- data_full_pop

# Calculate sex-adjusted incidence rates per 100,000 population
data_full_pop_full$rate <- 100000 * (data_full_pop_full$CL_CASES/data_full_pop_full$POP)
  
# Split data by sex for separate analysis
data_full_pop_male <- data_full_pop_full %>% filter(SEX == 1)
data_full_pop_female <- data_full_pop_full %>% filter(SEX == 2)

# Save sex-stratified datasets for model fitting
save(data_full_pop_male, file = file.path(DATA_DIR, paste0('male_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_pop_female, file = file.path(DATA_DIR, paste0('female_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit GSC models for each sex group
for(i in c(1,2)){
  # Load appropriate sex group data
  if(i == 1){
    load(file.path(DATA_DIR, paste0('male_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_male
  }else{
    load(file.path(DATA_DIR, paste0('female_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_female
  }
  
  # Note: Original MC model commented out for reference
  # Run MC Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_hisp, estimator = "mc", 
  #                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
  #if(i == 1){
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('male_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}else{
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('female_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
 
  
  # Run GSC (Interactive Fixed Effects) in Gsynth
  # This is the primary model used in main analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "ife", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, k=2, nboots = 1000, parallel = FALSE, inference = "nonparametric")
  
  # Save results by sex group
  if(i == 1){
    # Save male results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('male_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  }else{
    # Save female results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('female_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  } 
}

# Print completion message
cat("GSC models for sex groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")

