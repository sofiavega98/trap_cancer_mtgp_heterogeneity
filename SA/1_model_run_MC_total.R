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
load(file.path(DATA_DIR, "total_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_full <- data_full_hisp_noZero
# Convert cases into rates for gsynth
data_full_hisp_full$rate <- 100000 * (data_full_hisp_full$CL_CASES/data_full_hisp_full$POP)


  

  
  # Run MC Gysnth
  #fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_hisp, estimator = "mc", 
   #                 index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
  #                  se = TRUE, k=2, nboots = 1000, parallel = FALSE)
  
 # if(i == "Younger Age (0-4)"){
  #  # Save results
  #  save(fit_gsc ,file = paste0('~/Documents/Harvard/Lymphoma_Subsets/results/',date,'/young_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData'))
 # }else{
    # Save results
  #  save(fit_gsc ,file = paste0('~/Documents/Harvard/Lymphoma_Subsets/results/',date,'/old_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData'))
  #}
  
  # Run GSC in Gsynth
  
  fit_gsc <- gsynth(data = data_full_hisp_full, Y = "rate", D = "C", estimator = "ife", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, k=2, nboots = 1000, parallel = FALSE, inference = "nonparametric")
  # Save results

    # Save results
    save(fit_gsc ,file = paste0('~/Documents/Harvard/Lymphoma_Subsets/results/',date,'/SA/total_Lymphoma_',years[1],'_',years[length(years)],'_0-29_GSC.RData'))
  
  




