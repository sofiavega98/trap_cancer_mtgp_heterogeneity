################################################################################
source("config.R")

# Results Export Script: Generate Tables and Figures for Paper
# 
# This script generates all tables and figures presented in the paper:
# "MTGP_Lymphoma_Strata_paper.pdf"
#
# INPUTS:
# - MTGP model fits from 2_model_run_MTGP_*.R scripts
# - GSC model fits from 1_model_run_GSC_*.R scripts  
# - Cleaned data from 0_cleaning_data_*.R scripts
# - Functions from functions.R
#
# OUTPUTS:
# TABLES (saved to results/tables/):
# - mtgp_table.tex: Table 2 - Main MTGP treatment effect estimates
# - mtgp_table_sum_TE.tex: Table 4 - Sum of treatment effects across treated units/time periods
# - mtgp_table_ratio.tex: Table 3 - Treatment effect ratios
# - rhat_table.tex: Table S4 - MCMC convergence diagnostics (R-hat values)
# - gsc_table.tex: Table S5 - Sensitivity analysis GSC results
#
# FIGURES (saved to results/figures/):
# - ATT_plot_overall_trt1995.png: Figure 1 - Treatment effects over time (all subsets)
# - ATT_plot_[subset]_trt1995.png: Figures S4-S7 - Treatment effects by demographic subset
#
# This script processes posterior samples from MTGP models to calculate
# Average Treatment Effects (ATT) and generates publication-ready tables and figures.
################################################################################

# Load required libraries
library(tidyverse)
library(kableExtra)
library(cowplot)
library(gsynth)
library(posterior)

# Load all cleaned datasets for ATT calculation
# These datasets contain case counts, population, and treatment indicators by demographic subset
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/total_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/NH_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/H_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/male_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/female_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/white_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/black_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/other_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/nonHisp_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/hisp_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/young_Lymphoma_1988_2003_0-29.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/data/old_Lymphoma_1988_2003_0-29.RData")

# Load helper functions for ATT calculation and plotting
source("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/code/functions.R")

# Set output directories for tables and figures
dir_tables <- "~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/results/tables/"
dir_figures <- "~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/results/figures/"

# Set random seed for reproducibility
set.seed(501)

# Function to generate posterior samples from Poisson distribution
# This converts posterior means to posterior samples for uncertainty quantification
get_posterior_pois <- function(means) {
  
  # Generate random samples from Poisson distribution (n=1)
  posts <- array(rpois(length(means), lambda = means), dim = dim(means))
  
  # Return posterior samples of same dimension as means
  return(posts)
}

#######################################
## Table 2: Main Analysis MTGP Results ##
#######################################
# This table presents the main treatment effect estimates from the MTGP model
# Load Stan fit from 2_model_run_MTGP_total.R
dir_out <- "~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/Lymphoma_Subsets/results/11132024"

# Load and process MTGP model fits for each demographic subset
# Extract posterior samples and convert to Poisson samples for uncertainty quantification

# Overall population
load(paste0(dir_out,"/mtgp_fit_total_year.RData"))
fit_total <- fit
post_samp <-rstan::extract(fit) # posterior samples 
means <- post_samp$f # Extract f and store results
Y_pred <- get_posterior_pois(means) # Get posterior samples

# Lymphoma type (Hodgkin vs Non-Hodgkin)
load(paste0(dir_out,"/mtgp_fit_type_year.RData"))
fit_type <- fit
post_samp <-rstan::extract(fit) # posterior samples 
means_NH <- post_samp$f[,1,,] # Extract f and store results
means_H <- post_samp$f[,2,,] # Extract f and store results
Y_pred_NH <- get_posterior_pois(means_NH) # Get posterior samples
Y_pred_H <- get_posterior_pois(means_H) # Get posterior samples

# Sex (Male vs Female)
load(paste0(dir_out,"/mtgp_fit_sex_year.RData"))
fit_sex <- fit
post_samp <-rstan::extract(fit) # posterior samples 
means_male <- post_samp$f[,1,,] # Extract f and store results
means_female <- post_samp$f[,2,,] # Extract f and store results
Y_pred_male <- get_posterior_pois(means_male) # Get posterior samples
Y_pred_female <- get_posterior_pois(means_female) # Get posterior samples

# Race (White, Black, Other)
load(paste0(dir_out,"/mtgp_fit_race_year.RData"))
fit_race <- fit
post_samp <-rstan::extract(fit) # posterior samples 
means_white <- post_samp$f_mean[,1,,] # Extract f and store results
means_black <- post_samp$f_mean[,2,,] # Extract f and store results
means_other <- post_samp$f_mean[,3,,] # Extract f and store results
Y_pred_white <- get_posterior_pois(means_white) # Get posterior samples
Y_pred_black <- get_posterior_pois(means_black) # Get posterior samples
Y_pred_other <- get_posterior_pois(means_other) # Get posterior samples

# Hispanic ethnicity (Hispanic vs Non-Hispanic)
load(paste0(dir_out,"/mtgp_fit_hisp_year.RData"))
fit_hisp <- fit
post_samp <-rstan::extract(fit) # posterior samples 
means_nonHisp <- post_samp$f[,1,,] # Extract f and store results
means_Hisp <- post_samp$f[,2,,] # Extract f and store results
Y_pred_nonHisp <- get_posterior_pois(means_nonHisp) # Get posterior samples
Y_pred_Hisp <- get_posterior_pois(means_Hisp) # Get posterior samples

# Age groups (0-19 vs 20-29)
load(paste0(dir_out,"/mtgp_fit_age_year.RData"))
fit_age <- fit
post_samp <-rstan::extract(fit) # posterior samples 
means_young <- post_samp$f[,1,,] # Extract f and store results
means_old <- post_samp$f[,2,,] # Extract f and store results
Y_pred_young <- get_posterior_pois(means_young) # Get posterior samples
Y_pred_old <- get_posterior_pois(means_old) # Get posterior samples

# Convert 3D arrays to 2D matrices for ATT calculation
# This reshapes the posterior samples to match the expected format
to_2D <- function(Y_pred){
  # First, calculate the new dimensions
  n_rows <- dim(Y_pred)[1]  # 500
  n_cols <- prod(dim(Y_pred)[-1])  # 16 * 141 = 2256
  
  # Reshape the array
  # county 1 time 1, county 1 time 2, ....
  Y_pred_2d <- matrix(Y_pred, nrow = n_rows, ncol = n_cols)
  
  # Return
  return(Y_pred_2d)
}

# Apply 2D conversion to all posterior samples
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

# Calculate ATT for each demographic subset using treated units and time periods
# This applies the ATT_CACT_overall_table function to each subset
total_df <- ATT_CACT_overall_table(data_full_hisp_noZero, Y_pred_2d, "1995")
NH_df <- ATT_CACT_overall_table(data_full_hisp_NH, Y_pred_NH_2d, "1995")
H_df <- ATT_CACT_overall_table(data_full_hisp_H, Y_pred_H_2d, "1995")
male_df <- ATT_CACT_overall_table(data_full_hisp_male, Y_pred_male_2d, "1995")
female_df <- ATT_CACT_overall_table(data_full_hisp_female, Y_pred_female_2d, "1995")
white_df <- ATT_CACT_overall_table(data_full_hisp_white, Y_pred_white_2d, "1995")
black_df <- ATT_CACT_overall_table(data_full_hisp_black, Y_pred_black_2d, "1995")
other_df <- ATT_CACT_overall_table(data_full_hisp_other, Y_pred_other_2d, "1995")
nonHisp_df <- ATT_CACT_overall_table(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "1995")
Hisp_df <- ATT_CACT_overall_table(data_full_hisp_Hisp, Y_pred_Hisp_2d, "1995")
young_df <- ATT_CACT_overall_table(data_full_hisp_young, Y_pred_young_2d, "1995")
old_df <- ATT_CACT_overall_table(data_full_hisp_old, Y_pred_old_2d, "1995")

# Combine all demographic subsets into one table
df_mtgp_year <- rbind(total_df, NH_df, H_df, male_df, female_df, white_df, black_df, other_df,
                      nonHisp_df, Hisp_df, young_df, old_df)
df_mtgp_year <- df_mtgp_year[,2:7]
df_mtgp_year <- cbind(c("Overall","non-Hodgkin", "Hodgkin", "Male", "Female",
                        "White", "Black", "Other", "non-Hispanic", "Hispanic", "Ages 0-19",
                        "Ages 20-29"), df_mtgp_year)
colnames(df_mtgp_year) <- c("Subset", "CA ATT", "CA 95% CI", 
                            "CT ATT", "CT 95% CI", "Overall ATT", "Overall 95% CI" )
# Reorder columns to put overall results first
df_mtgp_year <- df_mtgp_year[,c(1,6,7,2,3,4,5)]

# Create and save LaTeX table
mtgp_table <- df_mtgp_year %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(mtgp_table, file = paste0(dir_tables,"mtgp_table.tex"))

cat("✓ Table 2 (Main MTGP Results) generated and saved\n")

##########################################################################
## Table 4: sum of the Y(1)-\hat{Y}(0) across treated units/time periods ##
#########################################################################
# This table shows the cumulative treatment effect across all treated observations
# Apply function to calculate sum of treatment effects
total_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_noZero, Y_pred_2d, "1995")
NH_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_NH, Y_pred_NH_2d, "1995")
H_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_H, Y_pred_H_2d, "1995")
male_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_male, Y_pred_male_2d, "1995")
female_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_female, Y_pred_female_2d, "1995")
white_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_white, Y_pred_white_2d, "1995")
black_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_black, Y_pred_black_2d, "1995")
other_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_other, Y_pred_other_2d, "1995")
nonHisp_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "1995")
Hisp_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_Hisp, Y_pred_Hisp_2d, "1995")
young_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_young, Y_pred_young_2d, "1995")
old_df <- ATT_CACT_overall_table_sum_TE(data_full_hisp_old, Y_pred_old_2d, "1995")

# Combine data
df_mtgp_year <- rbind(total_df, NH_df, H_df, male_df, female_df, white_df, black_df, other_df,
                      nonHisp_df, Hisp_df, young_df, old_df)
df_mtgp_year <- df_mtgp_year[,2:7]
df_mtgp_year <- cbind(c("Overall","non-Hodgkin", "Hodgkin", "Male", "Female",
                        "White", "Black", "Other", "non-Hispanic", "Hispanic", "Ages 0-19",
                        "Ages 20-29"), df_mtgp_year)
colnames(df_mtgp_year) <- c("Subset", "CA ATT", "CA 95% CI", 
                            "CT ATT", "CT 95% CI", "Overall ATT", "Overall 95% CI" )
# Reorder columns
df_mtgp_year <- df_mtgp_year[,c(1,6,7,2,3,4,5)]

# Make table
mtgp_table <- df_mtgp_year %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(mtgp_table, file = paste0(dir_tables,"mtgp_table_sum_TE.tex"))

cat("✓ Table 4 (Sum of Treatment Effects) generated and saved\n")

###################################
## Table 3: Treatment Effect Ratio ##
###################################
# This table shows treatment effect ratios (relative effects)

# Apply function to calculate ATT ratios
total_df <- ATT_CACT_overall_table_ratio(data_full_hisp_noZero, Y_pred_2d, "1995")
NH_df <- ATT_CACT_overall_table_ratio(data_full_hisp_NH, Y_pred_NH_2d, "1995")
H_df <- ATT_CACT_overall_table_ratio(data_full_hisp_H, Y_pred_H_2d, "1995")
male_df <- ATT_CACT_overall_table_ratio(data_full_hisp_male, Y_pred_male_2d, "1995")
female_df <- ATT_CACT_overall_table_ratio(data_full_hisp_female, Y_pred_female_2d, "1995")
white_df <- ATT_CACT_overall_table_ratio(data_full_hisp_white, Y_pred_white_2d, "1995")
black_df <- ATT_CACT_overall_table_ratio(data_full_hisp_black, Y_pred_black_2d, "1995")
other_df <- ATT_CACT_overall_table_ratio(data_full_hisp_other, Y_pred_other_2d, "1995")
nonHisp_df <- ATT_CACT_overall_table_ratio(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "1995")
Hisp_df <- ATT_CACT_overall_table_ratio(data_full_hisp_Hisp, Y_pred_Hisp_2d, "1995")
young_df <- ATT_CACT_overall_table_ratio(data_full_hisp_young, Y_pred_young_2d, "1995")
old_df <- ATT_CACT_overall_table_ratio(data_full_hisp_old, Y_pred_old_2d, "1995")

# Combine data
df_mtgp_year <- rbind(total_df, NH_df, H_df, male_df, female_df, white_df, black_df, other_df,
                      nonHisp_df, Hisp_df, young_df, old_df)
df_mtgp_year <- df_mtgp_year[,2:7]
df_mtgp_year <- cbind(c("Overall","non-Hodgkin", "Hodgkin", "Male", "Female",
                        "White", "Black", "Other", "non-Hispanic", "Hispanic", "Ages 0-19",
                        "Ages 20-29"), df_mtgp_year)
colnames(df_mtgp_year) <- c("Subset", "CA ATT", "CA 95% CI", 
                            "CT ATT", "CT 95% CI", "Overall ATT", "Overall 95% CI" )
# Reorder columns
df_mtgp_year <- df_mtgp_year[,c(1,6,7,2,3,4,5)]

# Make table
mtgp_table <- df_mtgp_year %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(mtgp_table, file = paste0(dir_tables,"mtgp_table_ratio.tex"))

########################################
## Table S4: Main Analysis MTGP R-Hat ##
########################################
library(rstan)

trt_fips <- data_full_hisp_NH %>% filter(C == 1) %>%
  distinct(FIPS) %>% pull(FIPS)

trt_year = 1995

# get Y0_obs
Y0_obs <- get_Y0_obs(data_full_hisp_NH,trt_fips,trt_year = trt_year)[[1]]

# get Y(1) i.e. observed treated values
Y1 <- get_Y1_obs1(data_full_hisp_NH) # 13*7=91 treated values; vector of length 91 (county1 time1, county1 time2,...)

# get population matrix
pop_mat <- get_pop_mat(data_full_hisp_NH)

# get population matrix for treated counties
n_trt <- length(trt_fips)

#get populations for treated counties
pop_trt <- as.vector(t(pop_mat[1:n_trt,]))

# get number of treated time points
m_trt <- sum(is.na(Y0_obs)*1)/n_trt

# get number of time points
m <- length(unique(data_full_hisp_NH$YEAR_DX))

# Get indices for treated times
ind <- c()
for(i in 0:(n_trt-1)){
  ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
}

total_summary <- summary(fit_total, pars = c("f_mean"))$summary[c(ind),]
total_mean <- mean(total_summary[,"Rhat"])

type_summary <- summary(fit_type, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_NH$FIPS))*m)),] # 1 outcome
type_mean <- mean(type_summary[,"Rhat"])

sex_summary <- summary(fit_sex, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_male$FIPS))*m)),] # 2 outcomes
sex_mean <- mean(sex_summary[,"Rhat"])

race_summary <- summary(fit_race, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_white$FIPS))*m), ind+(length(unique(data_full_hisp_white$FIPS))*(2*m))),] # 3 outcomes
race_mean <- mean(race_summary[,"Rhat"])

hisp_summary <- summary(fit_hisp, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_nonHisp$FIPS))*m)),]  # 2 outcomes
hisp_mean <- mean(hisp_summary[,"Rhat"])

age_summary <- summary(fit_age, pars = c("f_mean"))$summary[c(ind,ind+(length(unique(data_full_hisp_old$FIPS))*m)),]  # 2 outcomes
age_mean <- mean(age_summary[,"Rhat"])

# Combine data
df_rhat <- round(rbind(total_mean,type_mean, sex_mean, race_mean, hisp_mean, age_mean),3)
rownames(df_rhat) <- NULL
df_rhat <- cbind(c("Overall","Cancer Type", "Sex", "Race", "Hispanic Ethnicity", "Age of Diagnosis"), df_rhat)
colnames(df_rhat) <- c("Category","Rhat")

# Save table
rhat_table <- df_rhat %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(rhat_table, file = paste0(dir_tables,"rhat_table.tex"))
#rhat(fit)

############################################################
## Figure 1: Main Analysis MTGP Results Over Time Overall ##
############################################################
NH_plot <- ATT_overall_plot(data_full_hisp_NH, Y_pred_NH_2d, "NH", trt_year=1995, years=c(1988:2003)) + ggtitle("Non-Hodgkin")
H_plot <- ATT_overall_plot(data_full_hisp_H, Y_pred_H_2d, "H", trt_year=1995, years=c(1988:2003)) + ggtitle("Hodgkin")
male_plot <- ATT_overall_plot(data_full_hisp_male, Y_pred_male_2d, "male", trt_year=1995, years=c(1988:2003)) + ggtitle("Male")
female_plot <- ATT_overall_plot(data_full_hisp_female, Y_pred_female_2d, "female", trt_year=1995, years=c(1988:2003)) + ggtitle("Female")
white_plot <- ATT_overall_plot(data_full_hisp_white, Y_pred_white_2d, "white", trt_year=1995, years=c(1988:2003)) + ggtitle("White")
black_plot <- ATT_overall_plot(data_full_hisp_black, Y_pred_black_2d, "black", trt_year=1995, years=c(1988:2003)) + ggtitle("Black")
other_plot <- ATT_overall_plot(data_full_hisp_other, Y_pred_other_2d, "other", trt_year=1995, years=c(1988:2003)) + ggtitle("Other")
nonHisp_plot <- ATT_overall_plot(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "nonHisp", trt_year=1995, years=c(1988:2003)) + ggtitle("Non-Hispanic")
Hisp_plot <- ATT_overall_plot(data_full_hisp_Hisp, Y_pred_Hisp_2d, "Hisp", trt_year=1995, years=c(1988:2003)) + ggtitle("Hispanic")
young_plot <- ATT_overall_plot(data_full_hisp_young, Y_pred_young_2d, "young", trt_year=1995, years=c(1988:2003)) + ggtitle("Ages 0-19")
old_plot <- ATT_overall_plot(data_full_hisp_old, Y_pred_old_2d, "old", trt_year=1995, years=c(1988:2003)) + ggtitle("Ages 20-29")

# Create the layout
top_grid <- plot_grid(NH_plot, H_plot,NULL ,male_plot, female_plot, NULL, nonHisp_plot, Hisp_plot, NULL,young_plot, old_plot, NULL, ncol = 3)
bottom_row <- plot_grid(white_plot, black_plot, other_plot, ncol = 3)

final_plot <- plot_grid(top_grid, bottom_row, ncol = 1, rel_heights = c(4, 1))

# Display the final plot
final_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_overall_trt1995.png"),
         height = 12,
         width = 10,
         units = "in")


####################################################################
## Figures S4-S7: Main Analysis MTGP Results Over Time By Feature ##
####################################################################

NH_plot <- ATT_full_plot(data_full_hisp_NH, Y_pred_NH_2d, "NH", trt_year=1995, years=c(1988:2003))
NH_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_NH_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

H_plot <- ATT_full_plot(data_full_hisp_H, Y_pred_H_2d, "H", trt_year=1995, years=c(1988:2003))
H_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_H_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

male_plot <- ATT_full_plot(data_full_hisp_male, Y_pred_male_2d, "male", trt_year=1995, years=c(1988:2003))
male_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_male_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

female_plot <- ATT_full_plot(data_full_hisp_female, Y_pred_female_2d, "female", trt_year=1995, years=c(1988:2003))
female_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_female_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

white_plot <- ATT_full_plot(data_full_hisp_white, Y_pred_white_2d, "white", trt_year=1995, years=c(1988:2003))
white_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_white_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

black_plot <- ATT_full_plot(data_full_hisp_black, Y_pred_black_2d, "black", trt_year=1995, years=c(1988:2003))
black_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_black_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

other_plot <- ATT_full_plot(data_full_hisp_other, Y_pred_other_2d, "other", trt_year=1995, years=c(1988:2003))
other_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_other_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

nonHisp_plot <- ATT_full_plot(data_full_hisp_nonHisp, Y_pred_nonHisp_2d, "nonHisp", trt_year=1995, years=c(1988:2003))
nonHisp_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_nonHisp_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

Hisp_plot <- ATT_full_plot(data_full_hisp_Hisp, Y_pred_Hisp_2d, "Hisp", trt_year=1995, years=c(1988:2003))
Hisp_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_Hisp_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

young_plot <- ATT_full_plot(data_full_hisp_young, Y_pred_young_2d, "young", trt_year=1995, years=c(1988:2003))
young_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_young_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")

old_plot <- ATT_full_plot(data_full_hisp_old, Y_pred_old_2d, "old", trt_year=1995, years=c(1988:2003))
old_plot %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "ATT_plot_old_trt1995.png"),
         height = 5,
         width = 14,
         units = "in")


################################################
## Table S5: Sensitivity Analysis GSC Results ##
################################################
# Load GSC data
load(file.path(RESULTS_DIR, "09102024/total_Lymphoma_1988_2003_0-29_GSC.RData"))
total <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/NH_Lymphoma_1988_2003_0-29_GSC.RData"))
NH <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/H_Lymphoma_1988_2003_0-29_GSC.RData"))
H <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/male_Lymphoma_1988_2003_0-29_GSC.RData"))
male <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/female_Lymphoma_1988_2003_0-29_GSC.RData"))
female <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/white_Lymphoma_1988_2003_0-29_race_GSC.RData"))
white <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/black_Lymphoma_1988_2003_0-29_race_GSC.RData"))
black <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/other_Lymphoma_1988_2003_0-29_race_GSC.RData"))
other <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/nonHisp_Lymphoma_1988_2003_0-29_GSC.RData"))
nonHisp <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/Hisp_Lymphoma_1988_2003_0-29_GSC.RData"))
Hisp <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/young_Lymphoma_1988_2003_0-29_GSC.RData"))
young <- fit_gsc

load(file.path(RESULTS_DIR, "09102024/old_Lymphoma_1988_2003_0-29_GSC.RData"))
old <- fit_gsc

# Get ATT 
total_df <- ATT_CACT_overall_df(data_full_hisp_noZero, total, trt_year=1995)
NH_df <- ATT_CACT_overall_df(data_full_hisp_NH, NH, trt_year=1995)
H_df <- ATT_CACT_overall_df(data_full_hisp_H, H, trt_year=1995)
male_df <- ATT_CACT_overall_df(data_full_hisp_male, male, trt_year=1995)
female_df <- ATT_CACT_overall_df(data_full_hisp_female, female, trt_year=1995)
white_df <- ATT_CACT_overall_df(data_full_hisp_white, white, trt_year=1995)
black_df <- ATT_CACT_overall_df(data_full_hisp_black, black, trt_year=1995)
other_df <- ATT_CACT_overall_df(data_full_hisp_other, other, trt_year=1995)
nonHisp_df <- ATT_CACT_overall_df(data_full_hisp_nonHisp, nonHisp, trt_year=1995)
Hisp_df <- ATT_CACT_overall_df(data_full_hisp_Hisp, Hisp, trt_year=1995)
young_df <- ATT_CACT_overall_df(data_full_hisp_young, young, trt_year=1995)
old_df <- ATT_CACT_overall_df(data_full_hisp_old, old, trt_year=1995)

# Combine data
df_gsc <- rbind(total_df,NH_df, H_df, male_df, female_df, white_df, black_df, other_df,
                nonHisp_df, Hisp_df, young_df, old_df)
df_gsc <- df_gsc[,2:7]
df_gsc <- cbind(c("Overall","non-Hodgkin", "Hodgkin", "Male", "Female",
                  "White", "Black", "Other", "non-Hispanic", "Hispanic", "Ages 0-19",
                  "Ages 20-29"), df_gsc)
colnames(df_gsc) <- c("Subset", "CA ATT", "CA 95% CI", 
                      "CT ATT", "CT 95% CI", "Overall ATT", "Overall 95% CI" )
# Reorder columns
df_gsc <- df_gsc[,c(1,6,7,2,3,4,5)]

# Save table
gsc_table <- df_gsc %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(gsc_table, file = paste0(dir_tables,"gsc_table.tex"))


