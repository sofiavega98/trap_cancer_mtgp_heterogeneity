## Run MC in gsynth package 
source("config.R")


# Load libraries
library(gsynth)
library(tidyverse)

# Define date
date <- "09102024"

# Define years of interest
years <- c(1983:2003)

# Define treated year
treated_year <- 1990

# Load data
load(file.path(DATA_DIR, "SA/sex_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_full <- data_full_hisp_noZero
# Convert cases into rates for gsynth
data_full_hisp_full$rate <- 100000 * (data_full_hisp_full$CL_CASES/data_full_hisp_full$POP)
  
# Subset to a certain sex (1 = male ; 2 = female)
data_full_hisp_male <- data_full_hisp_full %>% filter(SEX == 1)
data_full_hisp_female <- data_full_hisp_full %>% filter(SEX == 2)

# Save sex-stratified datasets for model fitting
save(data_full_hisp_male, file = file.path(DATA_DIR, paste0('male_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_hisp_female, file = file.path(DATA_DIR, paste0('female_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit MC models for each sex group
for(i in c(1,2)){
  # Load appropriate sex group data
  if(i == 1){
    load(file.path(DATA_DIR, paste0('male_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_male
  }else{
    load(file.path(DATA_DIR, paste0('female_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_female
  }
  
  # Note: Original MC model commented out for reference
  # Run MC Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_pop, estimator = "mc", 
  #                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
  #if(i == 1){
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('male_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}else{
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('female_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}
  
  # Run MC (Matrix Completion) in Gsynth
  # This is the primary model used in sensitivity analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "mc", 
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
cat("MC models for sex groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")

