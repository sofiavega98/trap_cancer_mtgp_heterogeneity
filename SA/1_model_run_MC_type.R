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
load(file.path(DATA_DIR, "SA/type_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_full <- data_full_hisp_noZero
# Convert cases into rates for gsynth
data_full_hisp_full$rate <- 100000 * (data_full_hisp_full$CL_CASES/data_full_hisp_full$POP)


# Subset to a certain sex (1 = male ; 2 = female)
data_full_hisp_H <- data_full_hisp_full %>% filter(TYPE == "Hodgkin Lymphoma")
data_full_hisp_NH <- data_full_hisp_full %>% filter(TYPE == "Non-Hodgkin Lymphoma")

# Save lymphoma type-stratified datasets for model fitting
save(data_full_hisp_H, file = file.path(DATA_DIR, paste0('H_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
save(data_full_hisp_NH, file = file.path(DATA_DIR, paste0('NH_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Fit MC models for each lymphoma type group
for(i in c("H","NH")){
  # Load appropriate lymphoma type group data
  if(i == "H"){
    load(file.path(DATA_DIR, paste0('H_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_H
  }else{
    load(file.path(DATA_DIR, paste0('NH_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
    data_full_pop <- data_full_hisp_NH
  }
  
  # Note: Original MC model commented out for reference
  # Run MC Gysnth (Matrix Completion)
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_pop, estimator = "mc", 
  #                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
  #if(i == "H"){
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('H_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}else{
    # Save results
  #  save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('NH_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))
  #}
  
  # Run MC (Matrix Completion) in Gsynth
  # This is the primary model used in sensitivity analysis
  fit_gsc <- gsynth(data = data_full_pop, Y = "rate", D = "C", estimator = "mc", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, k=2, nboots = 1000, parallel = FALSE, inference = "nonparametric")
  
  # Save results by lymphoma type group
  if(i == "H"){
    # Save Hodgkin Lymphoma results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('H_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  }else{
    # Save Non-Hodgkin Lymphoma results
    save(fit_gsc, file = file.path(RESULTS_DIR, date, paste0('NH_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData')))
  } 
}

# Print completion message
cat("MC models for lymphoma type groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")



