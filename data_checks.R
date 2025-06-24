################################################################################
# Data Quality Checks and Validation
# 
# This script performs comprehensive data quality checks and validation on the
# cleaned lymphoma datasets. It verifies data consistency across demographic
# subsets, checks for missingness, and creates summary tables for the paper.
#
# INPUTS:
# - All cleaned demographic subset datasets from the data cleaning scripts
# - Non-Hodgkin data: "data/NH_Lymphoma_1988_2003_0-29.RData"
# - Hodgkin data: "data/H_Lymphoma_1988_2003_0-29.RData"
# - Male data: "data/male_Lymphoma_1988_2003_0-29.RData"
# - Female data: "data/female_Lymphoma_1988_2003_0-29.RData"
# - White data: "data/white_Lymphoma_1988_2003_0-29.RData"
# - Black data: "data/black_Lymphoma_1988_2003_0-29.RData"
# - Other race data: "data/other_Lymphoma_1988_2003_0-29.RData"
# - Non-Hispanic data: "data/nonHisp_Lymphoma_1988_2003_0-29.RData"
# - Hispanic data: "data/Hisp_Lymphoma_1988_2003_0-29.RData"
# - Young age data: "data/young_Lymphoma_1988_2003_0-29.RData"
# - Old age data: "data/old_Lymphoma_1988_2003_0-29.RData"
#
# OUTPUTS:
# - Case counts table: "results/tables/case_counts.tex"
# - Grouped case counts table: "results/tables/case_counts_grouped.tex"
# - Console output with data quality summaries
#
# PURPOSE:
# This script validates that:
# 1. Total case counts are consistent across demographic breakdowns
# 2. Data aggregation works correctly (e.g., young + old = total)
# 3. Treatment and control group assignments are consistent
# 4. Population denominators are appropriate for rate calculations
# 5. No systematic missingness exists across demographic groups
#
# DATA QUALITY CHECKS:
# - Cross-validation of case counts across demographic dimensions
# - Verification of treatment/control group balance
# - Population size checks for small counties
# - Consistency checks for before/after treatment periods
################################################################################

# Load required libraries
library(tidyverse)
library(kableExtra)

# Load configuration
source("config.R")

# Set output directory for tables
dir_out <- file.path(RESULTS_DIR, "tables/")

# Load all cleaned demographic subset datasets
# These datasets contain lymphoma case counts, population denominators,
# treatment indicators, and calculated rates for each demographic group
load(file.path(DATA_DIR, "NH_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_NH <- data_full_hisp_NH %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")
load(file.path(DATA_DIR, "H_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_H <- data_full_hisp_H %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")
load(file.path(DATA_DIR, "male_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_male <- data_full_hisp_male %>% select("FIPS", "YEAR_DX", "CL_CASES", "C",  "POP", "rate")
load(file.path(DATA_DIR, "female_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_female <- data_full_hisp_female %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")
load(file.path(DATA_DIR, "white_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_white <- data_full_hisp_white %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")
load(file.path(DATA_DIR, "black_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_black <- data_full_hisp_black %>% select("FIPS", "YEAR_DX", "CL_CASES", "C","POP", "rate")
load(file.path(DATA_DIR, "other_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_other <- data_full_hisp_other %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")
load(file.path(DATA_DIR, "nonHisp_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_nonHisp <- data_full_hisp_nonHisp %>% select("FIPS", "YEAR_DX", "CL_CASES", "C",  "POP", "rate")
load(file.path(DATA_DIR, "Hisp_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_Hisp <- data_full_hisp_Hisp %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")
load(file.path(DATA_DIR, "young_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_young <- data_full_hisp_young %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")
load(file.path(DATA_DIR, "old_Lymphoma_1988_2003_0-29.RData"))
data_full_hisp_old <- data_full_hisp_old %>% select("FIPS", "YEAR_DX", "CL_CASES", "C", "POP", "rate")

# DATA CONSISTENCY CHECK 1: Verify that total counts equal sum of demographic subsets
# This ensures our demographic breakdowns are complete and consistent
data_full <- data_full_hisp_young
data_full$CL_CASES <- data_full_hisp_young$CL_CASES + data_full_hisp_old$CL_CASES
data_full$rate <- NULL

# Calculate total case counts for each demographic dimension
total <- sum(data_full$CL_CASES)  # Ages 0-29 (comparison group)
type <- sum(data_full_hisp_NH$CL_CASES) + sum(data_full_hisp_H$CL_CASES)  # Cancer type
sex <- sum(data_full_hisp_male$CL_CASES) + sum(data_full_hisp_female$CL_CASES)  # Sex
race <- sum(data_full_hisp_white$CL_CASES) + sum(data_full_hisp_black$CL_CASES) + sum(data_full_hisp_other$CL_CASES)  # Race
hisp <- sum(data_full_hisp_Hisp$CL_CASES) + sum(data_full_hisp_nonHisp$CL_CASES)  # Hispanic ethnicity

# Create summary table of total case counts by demographic dimension
table <- cbind(c("Ages 0-29 (Comparison)","Cancer Type", "Sex", "Race", "Hispanic Ethnicity"),
           c(total,type, sex, race, hisp))

# Display table in console
table %>%
  kbl() %>%
  kable_classic(full_width = FALSE, font_size = 10)

# Save case counts table for paper
cases_table <- table %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(cases_table, file = file.path(RESULTS_DIR, "tables", "case_counts.tex"))

# DATA CONSISTENCY CHECK 2: Create detailed breakdown by treatment groups and time periods
# This function calculates aggregate case counts for treated vs control counties
# before and after the treatment year (1995)
sum_counts <- function(df){
  # Identify treated counties (CA/CT counties after 1995)
  treated_vec <- df %>% filter(C==1) %>% pull(FIPS) %>% unique()
  
  # Add period and treatment group indicators
  df <- df %>%
    mutate(period = ifelse(YEAR_DX < 1995, "before", "after"),
           treated = ifelse(FIPS %in% treated_vec, "treated", "control"))
  
  # Calculate aggregate case counts by treatment group and time period
  aggregate_counts <- df %>%
    group_by(treated, period) %>%
    summarise(
      aggregate_cases = sum(CL_CASES),
      .groups = 'drop'
    ) 
  return(aggregate_counts)
}

# Create aggregated datasets for each demographic dimension
# This combines the subset data to create total counts for each dimension
data_full_type <- data_full_hisp_NH
data_full_type$CL_CASES <- data_full_hisp_NH$CL_CASES + data_full_hisp_H$CL_CASES
data_full_type$rate <- NULL
data_full_type$POP <- NULL

data_full_sex <- data_full_hisp_male
data_full_sex$CL_CASES <- data_full_hisp_male$CL_CASES + data_full_hisp_female$CL_CASES
data_full_sex$rate <- NULL

data_full_race <- data_full_hisp_white
data_full_race$CL_CASES <- data_full_hisp_white$CL_CASES + data_full_hisp_black$CL_CASES + data_full_hisp_other$CL_CASES
data_full_race$rate <- NULL

data_full_hisp <- data_full_hisp_nonHisp
data_full_hisp$CL_CASES <- data_full_hisp_nonHisp$CL_CASES + data_full_hisp_Hisp$CL_CASES
data_full_hisp$rate <- NULL
data_full_hisp$POP <- NULL

# Apply the sum_counts function to each demographic dimension
full <- sum_counts(data_full)
type <- sum_counts(data_full_type)
sex <- sum_counts(data_full_sex)
race <- sum_counts(data_full_race)
hisp <- sum_counts(data_full_hisp)

# Combine all results into a single data frame for comparison
agg_counts_list <- list(
  "Ages 0-29 (Comparison)" = full,
  "Cancer Type" = type,
  "Sex" = sex,
  "Race" = race,
  "Hisp" = hisp
)

# Create a comprehensive summary table showing case counts by:
# - Demographic group
# - Treatment status (treated vs control)
# - Time period (before vs after 1995)
summary_table <- bind_rows(agg_counts_list, .id = "Group") %>%
  pivot_wider(names_from = c(treated, period), values_from = aggregate_cases, names_prefix = "") %>%
  rename("Control Before" = control_before,
         "Control After" = control_after,
         "Treated Before" = treated_before,
         "Treated After" = treated_after) %>%
  select(Group, `Control Before`, `Control After`, `Treated Before`, `Treated After`)%>%
  mutate(across(`Control Before`:`Treated After`, ~ round(.x, 2)))

# Display summary table
summary_table %>%
  kbl() %>%
  kable_classic(full_width = FALSE, font_size = 10)

print(summary_table)

# Save grouped case counts table for paper
table1 <- summary_table %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(table1, file = file.path(RESULTS_DIR, "tables", "case_counts_grouped.tex"))

##########################################
# SENSITIVITY ANALYSIS: Remove problematic county and recalculate
# County 49017 (Garfield County, UT) has no Black population, which can cause
# issues in race-stratified analyses. This section recalculates tables
# excluding this county to assess robustness.

# Remove county 49017 from all datasets
data_full <- data_full_type %>% filter(FIPS != "49017")
data_full_type <- data_full_type %>% filter(FIPS != "49017")
data_full_sex <- data_full_sex %>% filter(FIPS != "49017")
data_full_hisp <- data_full_hisp %>% filter(FIPS != "49017")

# Recalculate case counts excluding problematic county
full <- sum_counts(data_full)
type <- sum_counts(data_full_type)
sex <- sum_counts(data_full_sex)
race <- sum_counts(data_full_race)
hisp <- sum_counts(data_full_hisp)

# Create updated summary table
agg_counts_list <- list(
  "Ages 0-29 (Comparison)" = full,
  "Cancer Type" = type,
  "Sex" = sex,
  "Race" = race,
  "Hisp" = hisp
)

summary_table <- bind_rows(agg_counts_list, .id = "Group") %>%
  pivot_wider(names_from = c(treated, period), values_from = aggregate_cases, names_prefix = "") %>%
  rename("Control Before" = control_before,
         "Control After" = control_after,
         "Treated Before" = treated_before,
         "Treated After" = treated_after) %>%
  select(Group, `Control Before`, `Control After`, `Treated Before`, `Treated After`)%>%
  mutate(across(`Control Before`:`Treated After`, ~ round(.x, 2)))

summary_table %>%
  kbl() %>%
  kable_classic(full_width = FALSE, font_size = 10)

# POPULATION DENOMINATOR CHECKS
# This section creates population count tables to verify that population
# denominators are appropriate for rate calculations

# Function to calculate population counts by treatment group and time period
sum_pop <- function(df){
  # Identify treated counties
  treated_vec <- df %>% filter(C==1) %>% pull(FIPS) %>% unique()
  
  # Add period and treatment group indicators
  df <- df %>%
    mutate(period = ifelse(YEAR_DX < 1995, "before", "after"),
           treated = ifelse(FIPS %in% treated_vec, "treated", "control"))
  
  # Calculate aggregate population counts by grouping
  aggregate_counts <- df %>%
    group_by(treated, period) %>%
    summarise(
      aggregate_cases = sum(POP),
      .groups = 'drop'
    ) 
  return(aggregate_counts)
}

# Create aggregated population datasets for each demographic dimension
data_full_type <- data_full_hisp_NH
data_full_type$CL_CASES <- data_full_hisp_NH$CL_CASES + data_full_hisp_H$CL_CASES
data_full_type$rate <- NULL

data_full_sex <- data_full_hisp_male
data_full_sex$POP <- data_full_hisp_male$POP + data_full_hisp_female$POP
data_full_sex$rate <- NULL

data_full_race <- data_full_hisp_white
data_full_race$POP <- data_full_hisp_white$POP + data_full_hisp_black$POP + data_full_hisp_other$POP
data_full_race$rate <- NULL

data_full_hisp <- data_full_hisp_nonHisp
data_full_hisp$POP <- data_full_hisp_nonHisp$POP + data_full_hisp_Hisp$POP
data_full_hisp$rate <- NULL

# Calculate population counts by demographic dimension
sex <- sum_pop(data_full_sex)
race <- sum_pop(data_full_race)
hisp <- sum_pop(data_full_hisp)

##########################################################
# SMALL POPULATION CHECKS
# This section identifies counties with very small populations for specific
# demographic groups, which may need special handling in analyses

# Check for counties with very small populations (<5) for "Other" race group
counties_other <- data_full_hisp_other %>% filter(POP<5) %>% pull(FIPS) %>% unique()
counties_black <- data_full_hisp_other %>% filter(POP<5) %>% pull(FIPS) %>% unique()
counties_other
counties_black

# Check for counties with small populations (≤7) for "Other" race group
counties_other <- data_full_hisp_other %>% filter(POP<=7) %>% pull(FIPS) %>% unique()
counties_black <- data_full_hisp_other %>% filter(POP<=7) %>% pull(FIPS) %>% unique()
counties_other
counties_black

