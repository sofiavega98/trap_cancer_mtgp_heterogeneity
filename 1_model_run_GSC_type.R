################################################################################
source("config.R")

# Model Fitting: GSC (Interactive Fixed Effects) for Lymphoma Type
# 
# This script fits Interactive Fixed Effects (GSC) models using the gsynth package
# for lymphoma type analysis. It analyzes treatment effects on lymphoma incidence
# rates for Hodgkin and Non-Hodgkin lymphoma separately.
#
# INPUTS:
# - Lymphoma type data: "~/Documents/Harvard/Lymphoma_Subsets/data/type_Lymphoma_1988_2003_0-29.RData"
#
# OUTPUTS:
# - Hodgkin GSC results: "~/Documents/Harvard/Lymphoma_Subsets/results/09102024/H_Lymphoma_1988_2003_0-29_GSC.RData"
# - Non-Hodgkin GSC results: "~/Documents/Harvard/Lymphoma_Subsets/results/09102024/NH_Lymphoma_1988_2003_0-29_GSC.RData"
#
# MODEL SPECIFICATION:
# - Method: Interactive Fixed Effects (GSC) via gsynth package
# - Outcome: Lymphoma type-adjusted incidence rates per 100,000
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Fixed effects: Two-way (county and year)
# - Inference: Nonparametric bootstrap (1000 iterations)
# - Cross-validation: Yes (for factor number selection)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Lymphoma types: Hodgkin Lymphoma and Non-Hodgkin Lymphoma
# - Geographic scope: SEER registry counties
################################################################################

## Run GSC in gsynth package 
## Feature: Lymphoma Type (Hodgkin vs Non-Hodgkin)

# Load required libraries
library(gsynth)
library(tidyverse)

# Define analysis date for file naming
date <- "09102024"

# Define study parameters
years <- c(1988:2003)          # Study period
treated_year <- 1995           # Year treatment was implemented

# Load lymphoma type data
load(file.path(DATA_DIR, "type_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_full <- data_full_pop

# Calculate lymphoma type-adjusted incidence rates per 100,000 population
data_full_pop_full$rate <- 100000 * (data_full_pop_full$CL_CASES/data_full_pop_full$POP)

# Split data by lymphoma type for separate analysis
data_full_pop_H <- data_full_pop_full %>% filter(TYPE == "Hodgkin Lymphoma")
data_full_pop_NH <- data_full_pop_full %>% filter(TYPE == "Non-Hodgkin Lymphoma")

# Save lymphoma type-stratified datasets for model fitting
save(data_full_pop_H, file = file.path(DATA_DIR, paste0('HL_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_pop_NH, file = file.path(DATA_DIR, paste0('NHL_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit GSC models for each lymphoma type group
for(i in c("HL","NHL")){
  # Load appropriate lymphoma type group data
  if(i == "HL"){
    load(file.path(DATA_DIR, paste0('HL_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_H
  }else{
    load(file.path(DATA_DIR, paste0('NHL_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_pop_NH
  }
  
  # Note: Original MC model commented out for reference
  # Run MC Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_pop, estimator = "mc", 
  #                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
  #if(i == "HL"){
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('HL_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}else{
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('NHL_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}
  
  # Run GSC (Interactive Fixed Effects) in Gsynth
  # This is the primary model used in main analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "ife", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, k=2, nboots = 1000, parallel = FALSE, inference = "nonparametric")
  
  # Save results by lymphoma type group
  if(i == "HL"){
    # Save Hodgkin Lymphoma results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('HL_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  }else{
    # Save Non-Hodgkin Lymphoma results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('NHL_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  } 
}

# Print completion message
cat("GSC models for lymphoma type groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")



