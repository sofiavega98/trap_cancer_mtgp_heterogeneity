################################################################################
# Model Fitting: MTGP (Multi-Task Gaussian Process) for Total Population
# 
# This script fits Multi-Task Gaussian Process (MTGP) models using Stan for
# the total population analysis. It analyzes treatment effects on lymphoma
# incidence rates for the overall population (all demographic groups combined).
#
# INPUTS:
# - Cleaned data: "data/total_Lymphoma_1988_2003_0-29.RData"
# - Stan model: "mt_gp_pois_multigroup.stan"
#
# OUTPUTS:
# - MTGP model fit: "results/[date]/mtgp_fit_total_year.RData"
#
# MODEL SPECIFICATION:
# - Method: Multi-Task Gaussian Process (MTGP) via Stan
# - Outcome: Overall lymphoma case counts (Poisson)
# - Treatment: Binary indicator for CA/CT counties after 1995
# - Population: Overall population denominators
# - Latent functions: 15 for f, 10 for units
# - Inference: MCMC sampling (1000 iterations, 500 warmup)
#
# STUDY DESIGN:
# - Study period: 1988-2003 (16 years)
# - Treatment year: 1995 (California and Connecticut implement treatment)
# - Population: All ages 0-29 combined
# - Geographic scope: SEER registry counties
################################################################################

# Load configuration
source("config.R")

# Load required libraries
library(tidyverse)
library(rstan)

# Define study parameters
years <- YEARS
treated_year <- TREATED_YEAR
date <- DATE

cat("Fitting MTGP model for total population...\n")
cat("Study period:", min(years), "-", max(years), "\n")
cat("Treatment year:", treated_year, "\n")

# Load cleaned data
load(file.path(DATA_DIR, paste0("total_Lymphoma_",years[1],"_",years[length(years)],"_0-29.RData")))

# Prepare data for Stan MTGP model
# Select required columns and rename for consistency
data_full_pop <- data_full_pop %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP")

# Create case count matrix (y) for MTGP
# Add helper column to differentiate cases
df <- data_full_pop %>% select(CL_CASES) %>%
  mutate(CaseNumber = row_number())

# Pivot to wide format for Stan input
y <- df %>%
  pivot_wider(names_from = CaseNumber, values_from = CL_CASES) %>% 
  as.matrix()

# Create population array for Stan input
# Add helper column to differentiate population
df2 <- data_full_pop %>% select(POP) %>%
  mutate(CaseNumber = row_number())

pop_array <- df2 %>%
  pivot_wider(names_from = CaseNumber, values_from = POP) %>% 
  as.matrix()

# Define Stan model parameters
N = length(unique(data_full_pop$YEAR_DX))     # Number of observations (years)
D = length(unique(data_full_pop$FIPS))       # Number of units (counties)
num_outcomes = 1                              # Number of outcomes (total population)
n_k_f = 15                                   # Number of latent functions for f
n_k_d = 10                                   # Number of latent functions for units
x = years                                    # Univariate covariate (years)
population = pop_array                       # Population matrix [num_outcomes, N, D]
y = y                                        # Target variable matrix [num_outcomes, N * D]
num_treated = sum(data_full_pop$C)           # Number of treated units
control_idx = which(data_full_pop$C == 0)    # Control unit indices

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
mod <- stan_model("mt_gp_pois_multigroup.stan", auto_write = FALSE)

# Run Stan MTGP model
fit <- sampling(object = mod, data = stan_in, iter = 1000, 
                warmup = 500, pars = c('f_mean'), init_r = .1, seed = 501, 
                save_warmup = FALSE, chains = 1)

# Create results directory if it doesn't exist
dir.create(file.path(RESULTS_DIR, date), showWarnings = FALSE, recursive = TRUE)

# Save MTGP results
save(fit, file = file.path(RESULTS_DIR, date, "mtgp_fit_total_year.RData"))

cat("✓ MTGP model for total population completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date), "\n")
  
