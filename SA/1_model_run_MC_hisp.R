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
load(file.path(DATA_DIR, "SA/hispanic_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_full <- data_full_hisp_noZero
# Convert cases into rates for gsynth
data_full_hisp_full$rate <- 100000 * (data_full_hisp_full$CL_CASES/data_full_hisp_full$POP)


data_full_hisp_nonHisp <- data_full_hisp_full %>% filter(HISPANIC == 0)
data_full_hisp_Hisp <- data_full_hisp_full %>% filter(HISPANIC == 1)

# Save Hispanic ethnicity-stratified datasets for model fitting
save(data_full_hisp_nonHisp, file = file.path(DATA_DIR, paste0('nonHisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_hisp_Hisp, file = file.path(DATA_DIR, paste0('Hisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit MC models for each Hispanic ethnicity group
for(i in c(0,1)){
  # Load appropriate Hispanic ethnicity group data
  if(i == 0){
    load(file.path(DATA_DIR, paste0('nonHisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_nonHisp
  }else{
    load(file.path(DATA_DIR, paste0('Hisp_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_Hisp
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
  
  # Run MC (Matrix Completion) in Gsynth
  # This is the primary model used in sensitivity analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "mc", 
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
cat("MC models for Hispanic ethnicity groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")

