########################################################
source("config.R")

# Run Lymphoma data subsets on Ben-Michael et al. MTGP #
# Author: Sofia Vega                                   #
# Date: 08/13/24                                       #
########################################################

# Load libraries
library(tidyverse)
library(rstan)

# Define date
date <- "11132024"

# Load data
load(file.path(DATA_DIR, "total_Lymphoma_1983_2003_0-29.RData"))

# Combine into one dataframe
data_full <- data_full_hisp_noZero

# NOTE: I removed rows of all 0s for NH and H, for now, remove those rows in other datasets
# Identify counties with 0 rows
#young_0 <- data_full %>%
#  filter(!(FIPS %in% data_full_hisp_young$FIPS)) %>% distinct(FIPS)
#old_0 <- data_full%>%
#  filter(!(FIPS %in% data_full_hisp_old$FIPS)) %>% distinct(FIPS)
# Remove from other datasets
#data_full <- data_full %>% filter(!(FIPS %in% c(young_0$FIPS,old_0$FIPS)))
# Save for future use
#save(data_full, file=paste0("~/Documents/Harvard/Lymphoma_Subsets/data/data_outcome_full_",date,"_type.RData"))

#example cohort to get things like number treated etc
#ex_cohort <- data_full %>% filter(Cohort == "Older Population (20-29)")

# Create y

# Add a helper column to differentiate cases within each cohort
#df <- data_full %>% select(Cohort, CL_CASES) %>%
#  group_by(Cohort) %>%
#  mutate(CaseNumber = row_number()) %>%
#  ungroup()

# Pivot the data to wide format
#y <- df %>%
#  pivot_wider(names_from = CaseNumber, values_from = CL_CASES) %>% column_to_rownames(var = "Cohort")



# Try x as a vector of times
# Create input data for STAN
N =  length(unique(data_full$YEAR_DX))     #number of observations (years)
D = length(unique(data_full$FIPS))    #number of units (counties)
num_outcomes = 1 # number of outcomes we're modeling
n_k_f = 15     # number of latent functions for f??
n_k_d = 10   # number of latent functions for units??
x = 1983:2003    # univariate covariate (vector length N) should be year based on example use
population = matrix(data_full$POP, nrow = 1)# matrix[num_ooutcomes, N, D]
y =  matrix(data_full$CL_CASES, nrow=1 )       # target variable matrix [num_outcomes, N * D]
num_treated = sum(data_full$C)
control_idx = which(data_full$C == 0) # [N * D - num_treated] vector?? this is in MC format

stan_in <- list(N = N,
                D = D,
                num_outcomes = num_outcomes,
                n_k_f = n_k_f,
                n_k_d = n_k_d,
                x = x,
                population = population,
                y = y,
                num_treated = num_treated,
                control_idx = control_idx)
# Load stan
dir = "."
mod <- stan_model(paste0(dir,"mt_gp_pois_multigroup.stan"),auto_write=F)

# Run STAN Model
fit <- sampling(object = mod, data = stan_in, chains = 1, iter = 1000, 
                warmup=500, pars = c('f_mean'), init_r = .1, seed = 300, save_warmup = FALSE)

dir_out <- RESULTS_DIR

# Save MTGP results
save(fit, file = file.path(RESULTS_DIR, date, "SA/mtgp_fit_total_year.RData"))

# Print completion message
cat("MTGP model for total population (sensitivity analysis) completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date, "SA/mtgp_fit_total_year.RData"), "\n")
  
