################################################################################
# Sensitivity Analysis: Results Export and Tables Generation
# 
# This script exports and formats results from the sensitivity analysis models
# (MC and MTGP) for all demographic groups. It generates tables and figures
# comparing treatment effects across different assumptions and demographic
# subgroups using the extended study period (1983-2003) and earlier treatment
# year (1990).
#
# INPUTS:
# - MTGP model results from 2_model_run_MTGP_*.R scripts
# - MC model results from 1_model_run_MC_*.R scripts
# - Demographic-stratified datasets
# - Functions from functions.R
#
# OUTPUTS:
# - MTGP results table: "results/tables/SA/mtgp_table.tex"
# - R-hat diagnostic table: "results/tables/SA/rhat_table.tex"
# - GSC results table: "results/tables/SA/gsc_table.tex"
# - Treatment effect plots: "results/figures/ATT_plot_*_trt1990.png"
# - Additional tables and figures for sensitivity analysis
#
# SENSITIVITY ANALYSIS PARAMETERS:
# - Study period: 1983-2003 (21 years vs 16 years in main analysis)
# - Treatment year: 1990 (vs 1995 in main analysis)
# - Purpose: Test robustness of treatment effects under different assumptions
#
# TABLES GENERATED:
# 1. MTGP treatment effects by demographic group
# 2. Model convergence diagnostics (R-hat values)
# 3. GSC treatment effects for comparison
# 4. Comparison with main analysis results
#
# CONFIGURATION:
# All paths and parameters are managed through config.R
# This ensures portability across different systems
################################################################################

# Load configuration
source("config.R")

# Load required libraries
library(tidyverse)
library(kableExtra)
library(cowplot)
library(gsynth)

# Load demographic-stratified sensitivity analysis data
# These datasets contain the extended study period (1983-2003) data
load(file.path(DATA_DIR, "total_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "NH_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "H_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "male_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "female_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "white_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "black_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "other_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "nonHisp_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "hisp_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "young_Lymphoma_1983_2003_0-29.RData"))
load(file.path(DATA_DIR, "old_Lymphoma_1983_2003_0-29.RData"))

# Load utility functions from parent directory
source(file.path(BASE_DIR, "functions.R"))

# Set output directories for sensitivity analysis
dir_tables <- file.path(RESULTS_DIR, "tables", "SA")
dir_figures <- file.path(RESULTS_DIR, "figures", "SA")

# Create output directories if they don't exist
dir.create(dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

################################################################################
# UTILITY FUNCTIONS
################################################################################

# Function to get samples from the posterior distribution
# This generates random samples from Poisson distribution based on model means
get_posterior_pois <- function(means) {
  
  # Generate random samples from Poisson distribution (n=1)
  posts <- array(rpois(length(means), lambda = means), dim = dim(means))
  
  # Return posterior samples of same dimension as means
  return(posts)
}

# Function to convert 3D arrays to 2D matrices for analysis
# This reshapes the posterior samples for easier processing
to_2D <- function(Y_pred){
  # First, calculate the new dimensions
  n_rows <- dim(Y_pred)[1]  # Number of MCMC samples
  n_cols <- prod(dim(Y_pred)[-1])  # Total number of county-time combinations
  
  # Reshape the array: county 1 time 1, county 1 time 2, ....
  Y_pred_2d <- matrix(Y_pred, nrow = n_rows, ncol = n_cols)
  
  return(Y_pred_2d)
}

################################################################################
# LOAD MTGP MODEL RESULTS
# Load Stan fit results from sensitivity analysis MTGP models
################################################################################

# Set model output directory using configuration
dir_out <- file.path(RESULTS_DIR, DATE, "SA")

# Load total population results
load(file.path(dir_out, "mtgp_fit_total_year.RData"))
fit_total <- fit
post_samp <- rstan::extract(fit) # posterior samples 
means <- post_samp$f # Extract f and store results
Y_pred <- get_posterior_pois(means) # Get posterior samples

# Load lymphoma type results
load(file.path(dir_out, "mtgp_fit_type_year.RData"))
fit_type <- fit
post_samp <- rstan::extract(fit) # posterior samples 
means_NH <- post_samp$f[,1,,] # Extract f for Non-Hodgkin
means_H <- post_samp$f[,2,,] # Extract f for Hodgkin
Y_pred_NH <- get_posterior_pois(means_NH) # Get posterior samples
Y_pred_H <- get_posterior_pois(means_H) # Get posterior samples

# Load sex results
load(file.path(dir_out, "mtgp_fit_sex_year.RData"))
fit_sex <- fit
post_samp <- rstan::extract(fit) # posterior samples 
means_male <- post_samp$f[,1,,] # Extract f for males
means_female <- post_samp$f[,2,,] # Extract f for females
Y_pred_male <- get_posterior_pois(means_male) # Get posterior samples
Y_pred_female <- get_posterior_pois(means_female) # Get posterior samples

# Load race results
load(file.path(dir_out, "mtgp_fit_race_year.RData"))
fit_race <- fit
post_samp <- rstan::extract(fit_race) # posterior samples 
means_white <- post_samp$f[,1,,] # Extract f for whites
means_black <- post_samp$f[,2,,] # Extract f for blacks
means_other <- post_samp$f[,3,,] # Extract f for other races
Y_pred_white <- get_posterior_pois(means_white) # Get posterior samples
Y_pred_black <- get_posterior_pois(means_black) # Get posterior samples
Y_pred_other <- get_posterior_pois(means_other) # Get posterior samples

# Load Hispanic ethnicity results
load(file.path(dir_out, "mtgp_fit_hisp_year.RData"))
fit_hisp <- fit
post_samp <- rstan::extract(fit) # posterior samples 
means_nonHisp <- post_samp$f[,1,,] # Extract f for non-Hispanics
means_Hisp <- post_samp$f[,2,,] # Extract f for Hispanics
Y_pred_nonHisp <- get_posterior_pois(means_nonHisp) # Get posterior samples
Y_pred_Hisp <- get_posterior_pois(means_Hisp) # Get posterior samples

# Load age group results
load(file.path(dir_out, "mtgp_fit_age_year.RData"))
fit_age <- fit
post_samp <- rstan::extract(fit) # posterior samples 
means_young <- post_samp$f[,1,,] # Extract f for young age group
means_old <- post_samp$f[,2,,] # Extract f for old age group
Y_pred_young <- get_posterior_pois(means_young) # Get posterior samples
Y_pred_old <- get_posterior_pois(means_old) # Get posterior samples

################################################################################
# CONVERT PREDICTIONS TO 2D FORMAT
# Reshape all predictions for easier analysis
################################################################################

# Convert all predictions to 2D format for analysis
Y_pred_2d <- to_2D(Y_pred) 
Y_pred_NH_2d <- to_2D(Y_pred_NH) 
Y_pred_H_2d <- to_2D(Y_pred_H)
Y_pred_male_2d <- to_2D(Y_pred_male)
Y_pred_female_2d <- to_2D(Y_pred_female)
Y_pred_white_2d <- to_2D(Y_pred_white)
Y_pred_black_2d <- to_2D(Y_pred_black)
Y_pred_other_2d <- to_2D(Y_pred_other)
Y_pred_nonHisp_2d <- to_2D(Y_pred_nonHisp)
Y_pred_Hisp_2d <- to_2D(Y_pred_Hisp)
Y_pred_young_2d <- to_2D(Y_pred_young)
Y_pred_old_2d <- to_2D(Y_pred_old)

################################################################################
# TABLE 1: SENSITIVITY ANALYSIS MTGP RESULTS
# Calculate ATT for each demographic group using sensitivity analysis parameters
################################################################################

# Calculate treatment effects for each demographic subgroup
# Note: Using treatment year 1990 for sensitivity analysis (vs 1995 in main analysis)
total_df <- ATT_CACT_overall_table(data_full_hisp_noZero, Y_pred_2d, "1990")
NH_df <- ATT_CACT_overall_table(data_full_hisp_NH, Y_pred_NH_2d, "1990")
H_df <- ATT_CACT_overall_table(data_full_hisp_H, Y_pred_H_2d, "1990")
male_df <- ATT_CACT_overall_table(data_full_hisp_male, Y_pred_male_2d, "1990")
female_df <- ATT_CACT_overall_table(data_full_hisp_female, Y_pred_female_2d, "1990")
white_df <- ATT_CACT_overall_table(data_full_hisp_white, Y_pred_white_2d, "1990")
black_df <- ATT_CACT_overall_table(data_full_hisp_black, Y_pred_black_2d, "1990")
other_df <- ATT_CACT_overall_table(data_full_hisp_other, Y_pred_other_2d, "1990")
nonHisp_df <- ATT_CACT_overall_table(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "1990")
Hisp_df <- ATT_CACT_overall_table(data_full_hisp_Hisp, Y_pred_Hisp_2d, "1990")
young_df <- ATT_CACT_overall_table(data_full_hisp_young, Y_pred_young_2d, "1990")
old_df <- ATT_CACT_overall_table(data_full_hisp_old, Y_pred_old_2d, "1990")

# Combine all demographic group results into a single table
df_mtgp_year <- rbind(total_df, NH_df, H_df, male_df, female_df, white_df, black_df, other_df,
                      nonHisp_df, Hisp_df, young_df, old_df)
df_mtgp_year <- df_mtgp_year[,2:7]
df_mtgp_year <- cbind(c("Overall","non-Hodgkin", "Hodgkin", "Male", "Female",
                        "White", "Black", "Other", "non-Hispanic", "Hispanic", "Ages 0-19",
                        "Ages 20-29"), df_mtgp_year)
colnames(df_mtgp_year) <- c("Subset", "CA ATT", "CA 95% CI", 
                            "CT ATT", "CT 95% CI", "Overall ATT", "Overall 95% CI" )

# Reorder columns for better presentation (Overall first, then CA/CT)
df_mtgp_year <- df_mtgp_year[,c(1,6,7,2,3,4,5)]

# Generate LaTeX table for sensitivity analysis results
mtgp_table <- df_mtgp_year %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(mtgp_table, file = file.path(dir_tables, "mtgp_table.tex"))

################################################################################
# TABLE 2: SENSITIVITY ANALYSIS MTGP R-HAT DIAGNOSTICS
# Calculate R-hat values to assess model convergence
################################################################################

library(rstan)

# Get treated counties for diagnostic analysis
trt_fips <- data_full_hisp_NH %>% filter(C == 1) %>%
  distinct(FIPS) %>% pull(FIPS)

trt_year = 1990  # Sensitivity analysis treatment year

# Get observed data for diagnostic calculations
Y0_obs <- get_Y0_obs(data_full_hisp_NH, trt_fips, trt_year = trt_year)[[1]]
Y1 <- get_Y1_obs1(data_full_hisp_NH) # Observed treated values
pop_mat <- get_pop_mat(data_full_hisp_NH)

# Calculate dimensions for indexing
n_trt <- length(trt_fips)
pop_trt <- as.vector(t(pop_mat[1:n_trt,]))
m_trt <- sum(is.na(Y0_obs)*1)/n_trt
m <- length(unique(data_full_hisp_NH$YEAR_DX))

# Get indices for treated times
ind <- c()
for(i in 0:(n_trt-1)){
  ind <- append(ind, seq(((m-m_trt)+1), m) + (i*m))
}

# Calculate R-hat diagnostics for model convergence
# R-hat values close to 1.0 indicate good convergence
total_summary <- summary(fit_total, pars = c("f_mean"))$summary[c(ind),]
total_mean <- mean(total_summary[,"Rhat"])

type_summary <- summary(fit_type, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_NH$FIPS))*m)),]
type_mean <- mean(type_summary[,"Rhat"])

sex_summary <- summary(fit_sex, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_male$FIPS))*m)),]
sex_mean <- mean(sex_summary[,"Rhat"])

race_summary <- summary(fit_race, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_white$FIPS))*m)),]
race_mean <- mean(race_summary[,"Rhat"])

hisp_summary <- summary(fit_hisp, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_nonHisp$FIPS))*m)),]
hisp_mean <- mean(hisp_summary[,"Rhat"])

age_summary <- summary(fit_age, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_old$FIPS))*m)),]
age_mean <- mean(age_summary[,"Rhat"])

# Combine R-hat results into a table
df_rhat <- round(rbind(total_mean, type_mean, sex_mean, race_mean, hisp_mean, age_mean), 3)
rownames(df_rhat) <- NULL
df_rhat <- cbind(c("Overall","Cancer Type", "Sex", "Race", "Hispanic Ethnicity", "Age of Diagnosis"), df_rhat)
colnames(df_rhat) <- c("Category","Rhat")

# Save R-hat diagnostic table
rhat_table <- df_rhat %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(rhat_table, file = file.path(dir_tables, "rhat_table.tex"))

################################################################################
# FIGURES: SENSITIVITY ANALYSIS TREATMENT EFFECT PLOTS
# Generate treatment effect plots over time for each demographic subgroup
################################################################################

# Create individual treatment effect plots for each demographic group
# These show ATT over time with 95% confidence intervals
NH_plot <- ATT_overall_plot(data_full_hisp_NH, Y_pred_NH_2d, "NH", trt_year=1990, years=c(1983:2003)) + ggtitle("Non-Hodgkin")
H_plot <- ATT_overall_plot(data_full_hisp_H, Y_pred_H_2d, "H", trt_year=1990, years=c(1983:2003)) + ggtitle("Hodgkin")
male_plot <- ATT_overall_plot(data_full_hisp_male, Y_pred_male_2d, "male", trt_year=1990, years=c(1983:2003)) + ggtitle("Male")
female_plot <- ATT_overall_plot(data_full_hisp_female, Y_pred_female_2d, "female", trt_year=1990, years=c(1983:2003)) + ggtitle("Female")
white_plot <- ATT_overall_plot(data_full_hisp_white, Y_pred_white_2d, "white", trt_year=1990, years=c(1983:2003)) + ggtitle("White")
black_plot <- ATT_overall_plot(data_full_hisp_black, Y_pred_black_2d, "black", trt_year=1990, years=c(1983:2003)) + ggtitle("Black")
other_plot <- ATT_overall_plot(data_full_hisp_other, Y_pred_other_2d, "other", trt_year=1990, years=c(1983:2003)) + ggtitle("Other")
nonHisp_plot <- ATT_overall_plot(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "nonHisp", trt_year=1990, years=c(1983:2003)) + ggtitle("Non-Hispanic")
Hisp_plot <- ATT_overall_plot(data_full_hisp_Hisp, Y_pred_Hisp_2d, "Hisp", trt_year=1990, years=c(1983:2003)) + ggtitle("Hispanic")
young_plot <- ATT_overall_plot(data_full_hisp_young, Y_pred_young_2d, "young", trt_year=1990, years=c(1983:2003)) + ggtitle("Ages 0-19")
old_plot <- ATT_overall_plot(data_full_hisp_old, Y_pred_old_2d, "old", trt_year=1990, years=c(1983:2003)) + ggtitle("Ages 20-29")

# Create a comprehensive grid layout combining all demographic plots
top_grid <- plot_grid(NH_plot, H_plot, NULL, male_plot, female_plot, NULL, 
                      nonHisp_plot, Hisp_plot, NULL, young_plot, old_plot, NULL, ncol = 3)
bottom_row <- plot_grid(white_plot, black_plot, other_plot, ncol = 3)

final_plot <- plot_grid(top_grid, bottom_row, ncol = 1, rel_heights = c(4, 1))

# Save the combined treatment effect plot
final_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_overall_trt1990.png"),
         height = 14, width = 10, units = "in")

################################################################################
# INDIVIDUAL TREATMENT EFFECT PLOTS
# Generate separate plots for each demographic subgroup
################################################################################

# Save individual treatment effect plots for each demographic group
NH_plot <- ATT_full_plot(data_full_hisp_NH, Y_pred_NH_2d, "NH", trt_year=1990, years=c(1983:2003))
NH_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_NH_trt1990.png"),
         height = 5, width = 14, units = "in")

H_plot <- ATT_full_plot(data_full_hisp_H, Y_pred_H_2d, "H", trt_year=1990, years=c(1983:2003))
H_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_H_trt1990.png"),
         height = 5, width = 14, units = "in")

male_plot <- ATT_full_plot(data_full_hisp_male, Y_pred_male_2d, "male", trt_year=1990, years=c(1983:2003))
male_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_male_trt1990.png"),
         height = 5, width = 14, units = "in")

female_plot <- ATT_full_plot(data_full_hisp_female, Y_pred_female_2d, "female", trt_year=1990, years=c(1983:2003))
female_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_female_trt1990.png"),
         height = 5, width = 14, units = "in")

white_plot <- ATT_full_plot(data_full_hisp_white, Y_pred_white_2d, "white", trt_year=1990, years=c(1983:2003))
white_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_white_trt1990.png"),
         height = 5, width = 14, units = "in")

black_plot <- ATT_full_plot(data_full_hisp_black, Y_pred_black_2d, "black", trt_year=1990, years=c(1983:2003))
black_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_black_trt1990.png"),
         height = 5, width = 14, units = "in")

other_plot <- ATT_full_plot(data_full_hisp_other, Y_pred_other_2d, "other", trt_year=1990, years=c(1983:2003))
other_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_other_trt1990.png"),
         height = 5, width = 14, units = "in")

nonHisp_plot <- ATT_full_plot(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "nonHisp", trt_year=1990, years=c(1983:2003))
nonHisp_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_nonHisp_trt1990.png"),
         height = 5, width = 14, units = "in")

Hisp_plot <- ATT_full_plot(data_full_hisp_Hisp, Y_pred_Hisp_2d, "Hisp", trt_year=1990, years=c(1983:2003))
Hisp_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_Hisp_trt1990.png"),
         height = 5, width = 14, units = "in")

young_plot <- ATT_full_plot(data_full_hisp_young, Y_pred_young_2d, "young", trt_year=1990, years=c(1983:2003))
young_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_young_trt1990.png"),
         height = 5, width = 14, units = "in")

old_plot <- ATT_full_plot(data_full_hisp_old, Y_pred_old_2d, "old", trt_year=1990, years=c(1983:2003))
old_plot %>%
  ggsave(filename = file.path(dir_figures, "ATT_plot_old_trt1990.png"),
         height = 5, width = 14, units = "in")

################################################################################
# TABLE 3: SENSITIVITY ANALYSIS GSC RESULTS
# Load and process GSC model results for comparison with MTGP
################################################################################

# Load GSC model results from sensitivity analysis
# These provide an alternative estimation approach for comparison
load(file.path(RESULTS_DIR, DATE, "SA", "total_Lymphoma_1983_2003_0-29_GSC.RData"))
total <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "NH_Lymphoma_1983_2003_0-29_GSC.RData"))
NH <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "H_Lymphoma_1983_2003_0-29_GSC.RData"))
H <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "male_Lymphoma_1983_2003_0-29_GSC.RData"))
male <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "female_Lymphoma_1983_2003_0-29_GSC.RData"))
female <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "white_Lymphoma_1983_2003_0-29_race_GSC.RData"))
white <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "black_Lymphoma_1983_2003_0-29_race_GSC.RData"))
black <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "other_Lymphoma_1983_2003_0-29_race_GSC.RData"))
other <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "nonHisp_Lymphoma_1983_2003_0-29_GSC.RData"))
nonHisp <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "Hisp_Lymphoma_1983_2003_0-29_GSC.RData"))
Hisp <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "young_Lymphoma_1983_2003_0-29_GSC.RData"))
young <- fit_gsc

load(file.path(RESULTS_DIR, DATE, "SA", "old_Lymphoma_1983_2003_0-29_GSC.RData"))
old <- fit_gsc

# Calculate ATT for each demographic group using GSC models
total_df <- ATT_CACT_overall_df(data_full_hisp_noZero, total, trt_year=1990)
NH_df <- ATT_CACT_overall_df(data_full_hisp_NH, NH, trt_year=1990)
H_df <- ATT_CACT_overall_df(data_full_hisp_H, H, trt_year=1990)
male_df <- ATT_CACT_overall_df(data_full_hisp_male, male, trt_year=1990)
female_df <- ATT_CACT_overall_df(data_full_hisp_female, female, trt_year=1990)
white_df <- ATT_CACT_overall_df(data_full_hisp_white, white, trt_year=1990)
black_df <- ATT_CACT_overall_df(data_full_hisp_black, black, trt_year=1990)
other_df <- ATT_CACT_overall_df(data_full_hisp_other, other, trt_year=1990)
nonHisp_df <- ATT_CACT_overall_df(data_full_hisp_nonHisp, nonHisp, trt_year=1990)
Hisp_df <- ATT_CACT_overall_df(data_full_hisp_Hisp, Hisp, trt_year=1990)
young_df <- ATT_CACT_overall_df(data_full_hisp_young, young, trt_year=1990)
old_df <- ATT_CACT_overall_df(data_full_hisp_old, old, trt_year=1990)

# Combine GSC results into a single table
df_gsc <- rbind(total_df, NH_df, H_df, male_df, female_df, white_df, black_df, other_df,
                nonHisp_df, Hisp_df, young_df, old_df)
df_gsc <- df_gsc[,2:7]
df_gsc <- cbind(c("Overall","non-Hodgkin", "Hodgkin", "Male", "Female",
                  "White", "Black", "Other", "non-Hispanic", "Hispanic", "Ages 0-19",
                  "Ages 20-29"), df_gsc)
colnames(df_gsc) <- c("Subset", "CA ATT", "CA 95% CI", 
                      "CT ATT", "CT 95% CI", "Overall ATT", "Overall 95% CI" )
# Reorder columns for consistency with MTGP table
df_gsc <- df_gsc[,c(1,6,7,2,3,4,5)]

# Save GSC results table
gsc_table <- df_gsc %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(gsc_table, file = file.path(dir_tables, "gsc_table.tex"))

################################################################################
# COMPLETION
################################################################################

cat("=== SENSITIVITY ANALYSIS RESULTS EXPORT COMPLETE ===\n")
cat("Sensitivity analysis results have been generated and saved.\n")
cat("Check the following directories for outputs:\n")
cat("- Tables:", dir_tables, "\n")
cat("- Figures:", dir_figures, "\n\n")

cat("Key sensitivity analysis outputs generated:\n")
cat("- MTGP treatment effects table: mtgp_table.tex\n")
cat("- Model convergence diagnostics: rhat_table.tex\n")
cat("- GSC comparison results: gsc_table.tex\n")
cat("- Treatment effect plots: ATT_plot_*_trt1990.png\n\n")

cat("Sensitivity analysis parameters used:\n")
cat("- Study period: 1983-2003 (21 years)\n")
cat("- Treatment year: 1990\n")
cat("- Extended period for robustness testing\n\n")


