################################################################################
source("config.R")

# Sensitivity Analysis: MTGP Model Runner Script
# 
# This script runs Multi-Task Gaussian Process (MTGP) models for all demographic
# groups in the sensitivity analysis. It executes the MTGP models sequentially
# for age, Hispanic ethnicity, race, sex, total population, and lymphoma type
# using the extended study period (1983-2003) and earlier treatment year (1990).
#
# EXECUTION ORDER:
# 1. Age groups (0-19 vs 20-29 years)
# 2. Hispanic ethnicity (Hispanic vs Non-Hispanic)
# 3. Race (White, Black, Other)
# 4. Sex (Male vs Female)
# 5. Total population (all groups combined)
# 6. Lymphoma type (Hodgkin vs Non-Hodgkin)
#
# SENSITIVITY ANALYSIS PARAMETERS:
# - Study period: 1983-2003 (21 years vs 16 years in main analysis)
# - Treatment year: 1990 (vs 1995 in main analysis)
# - Purpose: Test robustness of treatment effects under different assumptions
#
# OUTPUTS:
# - MTGP results for each demographic group saved in results directory
# - Each model run clears memory to prevent conflicts
################################################################################

# Run MTGP models for all demographic groups in sensitivity analysis

print("Running MTGP models for age groups...")
source("SA/2_model_run_MTGP_age.R")
rm(list = ls())

print("Running MTGP models for Hispanic ethnicity...")
source("SA/2_model_run_MTGP_hisp.R")
rm(list = ls())

print("Running MTGP models for race...")
source("SA/2_model_run_MTGP_race.R")
rm(list = ls())

print("Running MTGP models for sex...")
source("SA/2_model_run_MTGP_sex.R")
rm(list = ls())

print("Running MTGP models for total population...")
source("SA/2_model_run_MTGP_total.R")
rm(list = ls())

print("Running MTGP models for lymphoma type...")
source("SA/2_model_run_MTGP_type.R")
rm(list = ls())

print("Sensitivity analysis MTGP models completed for all demographic groups!")
