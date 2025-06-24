################################################################################
source("config.R")

# Sensitivity Analysis: Model Fitting - MTGP (Multi-Task Gaussian Process) for Age Groups
# 
# This script fits Multi-Task Gaussian Process (MTGP) models using Stan for
# age group sensitivity analysis. It analyzes treatment effects on lymphoma
# incidence rates for young (0-19) and old (20-29) age groups using an extended
# study period and earlier treatment year.
#
# INPUTS:
# - Young age data: "data/SA/young_Lymphoma_1983_2003_0-29.RData"
# - Old age data: "data/SA/old_Lymphoma_1983_2003_0-29.RData"
# - Stan model: "mt_gp_pois_multigroup.stan"
#
# OUTPUTS:
# - MTGP results: "results/[date]/SA/mtgp_fit_age_year.RData"
#
# MODEL SPECIFICATION:
# - Method: Multi-Task Gaussian Process (MTGP) via Stan
# - Outcome: Lymphoma case counts (Poisson distribution)
# - Treatment: Binary indicator for CA/CT counties after 1990
# - Population: Offset for rate calculations
# - Latent functions: 15 for time trends, 10 for unit effects
# - Inference: MCMC sampling (1 chain, 1000 iterations, 500 warmup)
#
# SENSITIVITY ANALYSIS PARAMETERS:
# - Study period: 1983-2003 (21 years)
# - Treatment year: 1990 (vs 1995 in main analysis)
# - Age groups: 0-19 years (young) and 20-29 years (old)
# - Purpose: Test robustness under different time assumptions
################################################################################

########################################################
# Run Lymphoma data subsets on Ben-Michael et al. MTGP #
# Author: Sofia Vega                                   #
# Date: 08/13/24                                       #
# Sensitivity Analysis: Age Groups                     #
########################################################

# Load required libraries
library(tidyverse)
library(rstan)

# Define analysis date for file naming
date <- "11132024"

# Load age-stratified sensitivity analysis data
load(file.path(DATA_DIR, "SA/young_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_young <- data_full_hisp_young %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP")
data_full_hisp_young$Cohort <- "Younger Population (0-19)"

load(file.path(DATA_DIR, "SA/old_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_old <- data_full_hisp_old %>% select("FIPS", "YEAR_DX", "CL_CASES", "C",  "POP")
data_full_hisp_old$Cohort <- "Older Population (20-29)"

# Combine age groups into one dataframe for multi-task modeling
data_full <- rbind(data_full_hisp_young, data_full_hisp_old)

# Use example cohort to get treatment information
ex_cohort <- data_full %>% filter(Cohort == "Older Population (20-29)")

# Prepare outcome matrix (y) for Stan
# Add a helper column to differentiate cases within each cohort
df <- data_full %>% select(Cohort, CL_CASES) %>%
  group_by(Cohort) %>%
  mutate(CaseNumber = row_number()) %>%
  ungroup()

# Pivot the data to wide format for Stan input
y <- df %>%
  pivot_wider(names_from = CaseNumber, values_from = CL_CASES) %>% column_to_rownames(var = "Cohort")

# Prepare population matrix for Stan
# Add a helper column to differentiate population within each cohort
df2 <- data_full %>% select(Cohort, POP) %>%
  group_by(Cohort) %>%
  mutate(CaseNumber = row_number()) %>%
  ungroup()
pop_array <- df2 %>%
  pivot_wider(names_from = CaseNumber, values_from = POP) %>% column_to_rownames(var = "Cohort")

# Create input data for Stan model
N =  length(unique(data_full$YEAR_DX))     # Number of observations (years)
D = length(unique(data_full$FIPS))        # Number of units (counties)
num_outcomes = length(unique(data_full$Cohort)) # Number of outcomes (age groups)
n_k_f = 15     # Number of latent functions for time trends
n_k_d = 10     # Number of latent functions for unit effects
x = 1983:2003  # Time covariate (vector length N)
population = pop_array # Population matrix [num_outcomes, N, D]
y =  y        # Target variable matrix [num_outcomes, N * D]
num_treated = sum(ex_cohort$C)  # Number of treated units
control_idx = which(ex_cohort$C == 0) # Control unit indices

# Prepare Stan input data
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

# Load Stan model
dir = "."
mod <- stan_model(paste0(dir,"mt_gp_pois_multigroup.stan"),auto_write=F)

# Run Stan model for sensitivity analysis
fit <- sampling(object = mod, data = stan_in, chains = 1, iter = 1000, 
                warmup=500, pars = c('f_mean'), init_r = .1, seed = 300, save_warmup = FALSE)

# Save results
save(fit, file = file.path(RESULTS_DIR, date, "SA/mtgp_fit_age_year.RData"))

# Print completion message
cat("Sensitivity analysis MTGP model for age groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date, "SA/mtgp_fit_age_year.RData"), "\n") 
  
