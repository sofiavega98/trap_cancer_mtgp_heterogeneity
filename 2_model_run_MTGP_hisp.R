################################################################################
source("config.R")

# Model Fitting: MTGP (Multi-Task Gaussian Process) for Hispanic Ethnicity
# 
# This script fits Multi-Task Gaussian Process (MTGP) models using Stan for
# Hispanic ethnicity analysis. It analyzes treatment effects on lymphoma
# incidence rates for Hispanic and Non-Hispanic groups simultaneously.
#
# INPUTS:
# - Non-Hispanic data: "~/Documents/Harvard/Lymphoma_Subsets/data/nonHisp_Lymphoma_1988_2003_0-29.RData"
# - Hispanic data: "~/Documents/Harvard/Lymphoma_Subsets/data/Hisp_Lymphoma_1988_2003_0-29.RData"
# - Stan model: "~/Documents/Harvard/Lymphoma_Subsets/code/mt_gp_pois_multigroup.stan"
#
# OUTPUTS:
# - MTGP Hispanic ethnicity results: "~/Documents/Harvard/Lymphoma_Subsets/results/11132024/mtgp_fit_hisp_year.RData"
#
# MODEL SPECIFICATION:
# - Method: Multi-Task Gaussian Process (MTGP) via Stan
# - Outcome: Hispanic ethnicity-stratified lymphoma case counts (Poisson)
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Population: Hispanic ethnicity-specific population denominators
# - Latent functions: 15 for f, 10 for units
# - Inference: MCMC sampling (1000 iterations, 500 warmup)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Ethnicity groups: Non-Hispanic and Hispanic
# - Geographic scope: SEER registry counties (181 unique counties)
################################################################################

########################################################
# Run Lymphoma data subsets on Ben-Michael et al. MTGP #
# Author: Sofia Vega                                   #
# Date: 08/13/24                                       #
########################################################

## Feature: Hispanic Ethnicity (Hispanic vs Non-Hispanic)

# Load required libraries
library(tidyverse)
library(rstan)

# Define analysis date for file naming
date <- "11132024"

# Load Hispanic ethnicity-stratified data
load(file.path(DATA_DIR, "nonHisp_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_nonHisp <- data_full_pop_nonHisp %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP")
data_full_pop_nonHisp$Cohort <- "nonHisp"

load(file.path(DATA_DIR, "Hisp_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_Hisp <- data_full_pop_Hisp %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP")
data_full_pop_Hisp$Cohort <- "Hisp"

# Combine Hispanic ethnicity groups into one dataframe for MTGP analysis
data_full <- rbind(data_full_pop_nonHisp, data_full_pop_Hisp)
# 181 unique counties

# Use example cohort to get treatment information
ex_cohort <- data_full %>% filter(Cohort == "nonHisp")

# Create population matrices for each Hispanic ethnicity group
pop_mat <- list()
pop_mat[[1]] <- data_full %>% filter(Cohort == "nonHisp") %>% pull(POP) %>% matrix(nrow = length(unique(data_full$YEAR_DX)) , ncol = length(unique(data_full$FIPS))  , byrow = F)
pop_mat[[2]] <- data_full %>% filter(Cohort == "Hisp") %>% pull(POP) %>% matrix(nrow = length(unique(data_full$YEAR_DX)) , ncol = length(unique(data_full$FIPS))  , byrow = F)
pop_array <- simplify2array(pop_mat)
pop_array <- aperm(pop_array, c(3, 1, 2)) # Permute dimensions to [2, 16, 197] (cohorts, years, counties)

# Create case count matrix (y) for MTGP
# Add helper column to differentiate cases within each cohort
df <- data_full %>% select(Cohort, CL_CASES) %>%
  group_by(Cohort) %>%
  mutate(CaseNumber = row_number()) %>%
  ungroup()

# Pivot to wide format for Stan input
y <- df %>%
  pivot_wider(names_from = CaseNumber, values_from = CL_CASES) %>% column_to_rownames(var = "Cohort")

# Create population array for Stan input
# Add helper column to differentiate population within each cohort
df2 <- data_full %>% select(Cohort, POP) %>%
  group_by(Cohort) %>%
  mutate(CaseNumber = row_number()) %>%
  ungroup()
pop_array <- df2 %>%
  pivot_wider(names_from = CaseNumber, values_from = POP) %>% column_to_rownames(var = "Cohort")

# Define Stan model parameters
N =  length(unique(data_full$YEAR_DX))     # Number of observations (years)
D = length(unique(data_full$FIPS))    # Number of units (counties)
num_outcomes = length(unique(data_full$Cohort)) # Number of outcomes (Hispanic ethnicity groups)
n_k_f = 15     # Number of latent functions for f
n_k_d = 10   # Number of latent functions for units
x = 1988:2003    # Univariate covariate (years)
population = pop_array # Population matrix [num_outcomes, N, D]
y =  y        # Target variable matrix [num_outcomes, N * D]
num_treated = sum(ex_cohort$C) # Number of treated units
control_idx = which(ex_cohort$C == 0) # Control unit indices

# Create Stan input data list
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

# Run Stan MTGP model
fit <- sampling(object = mod, data = stan_in, iter = 1000, 
                warmup=500, pars = c('f_mean'), init_r = .1, seed = 501, save_warmup = FALSE,
                chains = 1) # Adjust adapt_delta here
#, control = list(adapt_delta = 0.99, max_treedepth = 15)

# Save MTGP results
save(fit, file = file.path(RESULTS_DIR, date, "mtgp_fit_hisp_year.RData"))

# Print completion message
cat("MTGP model for Hispanic ethnicity groups completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date, "mtgp_fit_hisp_year.RData"), "\n")
  
