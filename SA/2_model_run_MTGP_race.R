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
load(file.path(DATA_DIR, "SA/white_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_white <- data_full_hisp_white %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP")
data_full_hisp_white$Cohort <- "white"

load(file.path(DATA_DIR, "SA/black_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_black <- data_full_hisp_black %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP")
data_full_hisp_black$Cohort <- "black"

load(file.path(DATA_DIR, "SA/other_Lymphoma_1983_2003_0-29.RData"))
data_full_hisp_other <- data_full_hisp_other %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP")
data_full_hisp_other$Cohort <- "other"

# Combine into one dataframe
data_full <- rbind(data_full_hisp_white, data_full_hisp_black, data_full_hisp_other)
#181 unique counties


#example cohort to get things like number treated etc
ex_cohort <- data_full %>% filter(Cohort == "white")


# Create y

# Add a helper column to differentiate cases within each cohort
df <- data_full %>% select(Cohort, CL_CASES) %>%
  group_by(Cohort) %>%
  mutate(CaseNumber = row_number()) %>%
  ungroup()

# Pivot the data to wide format
y <- df %>%
  pivot_wider(names_from = CaseNumber, values_from = CL_CASES) %>% column_to_rownames(var = "Cohort")

# Add a helper column to differentiate cases within each cohort
df2 <- data_full %>% select(Cohort, POP) %>%
  group_by(Cohort) %>%
  mutate(CaseNumber = row_number()) %>%
  ungroup()
pop_array <- df2 %>%
  pivot_wider(names_from = CaseNumber, values_from = POP) %>% column_to_rownames(var = "Cohort")



# Try x as a vector of times
# Create input data for STAN
N =  length(unique(data_full$YEAR_DX))     #number of observations (years)
D = length(unique(data_full$FIPS))    #number of units (counties)
num_outcomes = length(unique(data_full$Cohort)) # number of outcomes we're modeling
n_k_f = 15     # number of latent functions for f??
n_k_d = 10   # number of latent functions for units??
x = 1983:2003    # univariate covariate (vector length N) should be year based on example use
population = pop_array # matrix[num_outcomes, N, D]
y =  y        # target variable matrix [num_outcomes, N * D]
num_treated = sum(ex_cohort$C)
control_idx = which(ex_cohort$C == 0) # [N * D - num_treated] vector?? this is in MC format

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


# Save MTGP results
save(fit, file = file.path(RESULTS_DIR, date, "SA/mtgp_fit_race_year.RData"))

# Print completion message
cat("MTGP model for race groups (sensitivity analysis) completed!\n")
cat("Results saved in:", file.path(RESULTS_DIR, date, "SA/mtgp_fit_race_year.RData"), "\n")
  
