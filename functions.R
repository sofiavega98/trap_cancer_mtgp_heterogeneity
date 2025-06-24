################################################################################
# Helper Functions for Lymphoma Subsets Analysis
# 
# This file contains utility functions for:
# 1. GSC (Generalized Synthetic Control) analysis
# 2. ATT (Average Treatment Effect) calculation
# 3. Data processing and manipulation
# 4. Plotting and visualization
#
# These functions are used by the main analysis scripts to:
# - Calculate treatment effects from model outputs
# - Process posterior samples from MTGP models
# - Generate tables and figures for the paper
# - Handle data transformations and aggregations
################################################################################

# Load configuration
source("config.R")

################################################################################
# GSC (Generalized Synthetic Control) Functions
################################################################################

# Function for aggregate ATT for gsynth
# Calculates the average treatment effect from GSC model outputs
#
# INPUTS:
#   Y0_est_rate: Vector of estimated counterfactual rates from GSC model
#   Y1: Vector of observed treated outcomes (case counts)
#   pop_trt: Vector of population denominators for treated counties
#   n_trt: Number of treated counties
#   fit_gsc: GSC model fit object containing confidence intervals
#
# OUTPUTS:
#   List containing:
#     - ATT: Average Treatment Effect (rate difference)
#     - CI: 95% confidence interval from GSC model
#
# USAGE:
#   result <- get_agg_ATT_gsc(Y0_est_rate, Y1, pop_trt, n_trt, fit_gsc)
#   att <- result$ATT
#   ci <- result$CI
get_agg_ATT_gsc <- function(Y0_est_rate, Y1, pop_trt, n_trt, fit_gsc) {
  
  # Convert to rates
  #Y0_est_rate <- as.vector(t(fit_gsc$Y.ct)[,(m_trt-1):m])
  
  # subset population to be treated counties at treated times
  #pop_trt_temp <- matrix(pop_trt,nrow = n_trt, byrow = T)[,(m-m_trt+1):m] #(13 x 7)
  pop_trt_temp <- matrix(pop_trt, nrow = n_trt, byrow = T)
  
  # convert Y1 into rate per 100,000 population
  Y1_rate <- 100000 * Y1/as.vector(t(pop_trt_temp))
  
  # Calculate ATT as mean difference between observed and counterfactual rates
  ATT <- mean(as.vector(Y1_rate) - Y0_est_rate) # this matches gsynth
  
  # Get 95% CI from GSC model (Note: I can only get CIs for all treated years)
  CI <- fit_gsc$est.avg[3:4]
  
  return(list(ATT = ATT, CI = c(CI)))
}

################################################################################
# Data Processing Functions
################################################################################

# Extract observed treated values from dataset
# This function pulls the observed treated values using the treatment indicator variable
#
# INPUTS:
#   data: Dataframe with columns FIPS, CL_CASES, YEAR_DX, C
#         - FIPS: County FIPS codes
#         - CL_CASES: Case counts
#         - YEAR_DX: Diagnosis year
#         - C: Treatment indicator (1 = treated, 0 = control)
#
# OUTPUTS:
#   Y1_obs1: Vector of observed case counts for treated counties
#
# USAGE:
#   treated_cases <- get_Y1_obs1(data_full_hisp)
get_Y1_obs1 <- function(data) {
  #Input:
  ## data: the data set with columns FIPS, CL_CASES, YEAR_DX, C
  ## C is an indicator column indicating treatment
  #Output:
  ## Y1_obs: a vector of treated outcomes
  
  Y1_obs1 <- data[which(data$C==1),]$CL_CASES
  
  return(Y1_obs1)
}

# Create matrix of observed non-treated values
# This obtains a matrix of observed non-treated values with NAs for treated observations
#
# INPUTS:
#   data: Dataframe with columns FIPS, CL_CASES, YEAR_DX
#   trt_fips: Vector of treated county FIPS codes
#   trt_year: Treatment start year (default "1996")
#
# OUTPUTS:
#   List containing:
#     - Y0_obs1: Matrix of observed non-treated values with NAs for treated values
#     - Y0_obs2: Matrix of observed non-treated values with 99999 for treated values
#                (for Stan compatibility)
#
# USAGE:
#   result <- get_Y0_obs(data_full_hisp, trt_fips, "1995")
#   Y0_obs1 <- result[[1]]  # Matrix with NAs
#   Y0_obs2 <- result[[2]]  # Matrix with 99999 placeholders
get_Y0_obs <- function(data, trt_fips, trt_year = "1996") {
  # Input:
  ## data: the data set with columns FIPS, CL_CASES, YEAR_DX
  ## trt_fips: a vector of treated counties
  ## trt_year: the year where treatment starts (default "1996")
  # Output:
  ## Y0_obs1 matrix of observed nontreated values with NAs for trt values
  ## Y0_obs2 matrix of observed nontreated values with 99999 for trt values
  
  #observed non treated values
  Y0_obs1 <- pivot_wider(data, id_cols = FIPS, 
                         names_from = YEAR_DX, values_from = CL_CASES) 
  
  #put NAs where counties are "treated" (Note: the +1 accounts for the FIPS column)
  Y0_obs1[which(Y0_obs1$FIPS %in% trt_fips),
          (which(unique(data$YEAR_DX)==trt_year)+1):(length(unique(data$YEAR_DX))+1)] <-NA
  #remove FIPS column
  Y0_obs1 <- Y0_obs1[,-1]
  
  ## stan doesn't allow missing values so we put a numeric "placeholder" in the missings
  Y0_obs2 <- Y0_obs1
  Y0_obs2[is.na(Y0_obs1)]<-(99999)
  
  return(list(Y0_obs1,Y0_obs2))
}

# Convert population data to matrix format
# This function pivots the population column in the data to be a matrix
# with counties in the rows and years in the columns
#
# INPUTS:
#   data: Dataframe with columns FIPS, YEAR_DX, POP
#         - FIPS: County FIPS codes
#         - YEAR_DX: Diagnosis year
#         - POP: Population denominators
#
# OUTPUTS:
#   pop_mat: Matrix with counties as rows and years as columns
#            Dimensions: [n_counties x n_years]
#
# USAGE:
#   population_matrix <- get_pop_mat(data_full_hisp)
get_pop_mat <- function(data) {
  
  pop_mat <- matrix(data$POP, nrow=length(unique(data$FIPS)),
                    ncol=length(unique(data$YEAR_DX)), byrow=T)
  
  return(pop_mat)
}

# Create missing value indicator matrix
# Creates a matrix indicating which observations are missing (treated)
#
# INPUTS:
#   Y0_obs1: Matrix of observed values with NAs for treated values
#
# OUTPUTS:
#   y_miss: Matrix with 1s for treated values and 0s for non-treated values
#           Same dimensions as input matrix
#
# USAGE:
#   missing_indicator <- get_y_miss(Y0_obs1)
get_y_miss <- function(Y0_obs1) {
  #Input: 
  ## Y0_obs1: a matrix of observed values with NAs for trt values
  #Output: 
  ## y_miss: a matrix with 1s for treated values and 0s for non-treated values
  y_miss <- as.matrix(is.na(Y0_obs1),
                      nrow=nrow(Y0_obs1), ncol=ncol(Y0_obs1), byrow=T)*1
  return(y_miss)
}

################################################################################
# ATT Calculation Functions for GSC Models
################################################################################

# Main function to calculate ATT from GSC model outputs
# Separates California and Connecticut effects and calculates overall ATT
#
# INPUTS:
#   data_full_hisp: Dataframe with demographic data including FIPS, CL_CASES, YEAR_DX, C
#   fit_gsc: GSC model fit object containing counterfactual predictions
#   trt_year: Treatment start year (e.g., "1995")
#
# OUTPUTS:
#   df: Dataframe with columns:
#       - Model: Model name ("Gsynth")
#       - ATT: Treatment effects for CA, CT, and Overall
#       - 95% CI: Confidence intervals (NA for CA/CT, actual CI for Overall)
#
# USAGE:
#   gsc_results <- ATT_CACT_overall_df(data_full_hisp, fit_gsc, "1995")
ATT_CACT_overall_df <- function(data_full_hisp, fit_gsc, trt_year){
  # Edit 0912: seq((m_trt-1),m) to seq(((m-m_trt)+1),m)
  
  ############
  ## Set-Up ##
  ############
  # get treated counties
  trt_fips <- data_full_hisp %>% filter(C == 1) %>%
    distinct(FIPS) %>% pull(FIPS)
  
  # get Y0_obs
  Y0_obs <- get_Y0_obs(data_full_hisp,trt_fips,trt_year = trt_year)[[1]]
  
  # get Y(1) i.e. observed treated values
  Y1 <- get_Y1_obs1(data_full_hisp) # 13*7=91 treated values; vector of length 91 (county1 time1, county1 time2,...)
  
  # get population matrix
  pop_mat <- get_pop_mat(data_full_hisp)
  
  # get population matrix for treated counties
  n_trt <- length(trt_fips)
  
  #get populations for treated counties
  pop_trt <- as.vector(t(pop_mat[1:n_trt,]))
  
  # get number of treated time points
  m_trt <- sum(is.na(Y0_obs)*1)/n_trt
  
  # get number of time points
  m <- length(unique(data_full_hisp$YEAR_DX))
  
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  # Split CA and CT for separate analysis
  
  # Get CA (California counties)
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  # Find which indices are CA
  ind <- c()
  for(i in 0:(n_trt_CA-1)){
    ind <- append(ind,seq((m-m_trt+1),m) + (i*m))
  }
  
  # Subset CA counterfactual predictions
  Mu_trt_gsc_CA <- as.vector(t(fit_gsc$Y.ct)[1:(length(unique(treated1$FIPS))),(m-m_trt+1):m])
  
  # Subset Population for CA
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CA <- pop_trt_CA$POP
  
  # Get CT (Connecticut counties)
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq((m-m_trt+1),m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_gsc_CT <- as.vector(t(fit_gsc$Y.ct)[(n_trt_CA+1):(n_trt_CA + n_trt_CT),(m-m_trt+1):m])
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CT <- pop_trt_CT$POP
  
  # Calculate ATT for California
  agg_ATT_gsc_CA_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CA,Y1_CA,pop_trt_CA,n_trt_CA,fit_gsc)$ATT
  
  # Calculate ATT for Connecticut
  agg_ATT_gsc_CT_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CT,Y1_CT,pop_trt_CT,n_trt_CT,fit_gsc)$ATT
  
  # Helper function to format confidence intervals
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  # Calculate overall ATT
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq((m-m_trt+1),m) + (i*m))
  }
  
  Mu_trt_gsc <- as.vector(t(fit_gsc$Y.ct)[,(m-m_trt+1):m])
  
  agg_ATT_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$ATT
  agg_bounds_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$CI
  
  CI_gsc <- CI_fun(agg_bounds_gsc_sub)
  
  # Make kable table
  Model <- c("Gsynth")
  df <- cbind(Model, 
              round(c(agg_ATT_gsc_CA_sub), 2),
              c("NA"),
              round(c(agg_ATT_gsc_CT_sub), 2),
              c("NA"), 
              round(c(agg_ATT_gsc_sub), 2),
              c(CI_gsc))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI", "ATT", "95% CI")
  
  return(df)
}

################################################################################
# ATT Calculation Functions for MTGP Models
################################################################################

# Calculate aggregate ATT from MTGP posterior samples
# Converts predictions to rates and calculates treatment effects across MCMC samples
#
# INPUTS:
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   Y1: Vector of observed treated outcomes
#   pop_trt: Vector of population denominators
#
# OUTPUTS:
#   List containing:
#     - ATT_vec: Vector of ATT values for each MCMC sample
#     - ATT: Median ATT across MCMC samples
#     - CI: 95% confidence interval (2.5% and 97.5% quantiles)
#
# USAGE:
#   result <- get_agg_ATT(Mu_trt, Y1, pop_trt)
#   att <- result$ATT
#   ci <- result$CI
get_agg_ATT <- function(Mu_trt,Y1,pop_trt) {
  
  # Convert to rates
  Mu_trt_rate <- 100000 * sweep(Mu_trt, 2, pop_trt, "/" )
  Y1_rate <- 100000 * Y1/pop_trt
  
  # Get an ATT for each MCMC iteration
  ATT_vec <- rowMeans(-1*sweep(Mu_trt_rate, 2, Y1_rate, "-")) #ATT= Y1-Y0 = -1*(Y0-Y1)
  
  # Get median of MCMC ATT
  ATT <- median(ATT_vec)
  
  # Get 95% CI
  CI <- quantile(ATT_vec, probs = c(.025,.975),na.rm = TRUE)
  
  return(list(ATT_vec = ATT_vec, ATT = ATT, CI = c(CI)))
}

# Calculate aggregate ATT ratio from MTGP posterior samples
# Computes the ratio of observed to counterfactual total cases
#
# INPUTS:
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   Y1: Vector of observed treated outcomes
#   pop_trt: Vector of population denominators (not used in ratio calculation)
#
# OUTPUTS:
#   List containing:
#     - ATT_vec: Vector of ATT ratio values for each MCMC sample
#     - ATT: Median ATT ratio across MCMC samples
#     - CI: 95% confidence interval (2.5% and 97.5% quantiles)
#
# USAGE:
#   result <- get_agg_ATT_ratio(Mu_trt, Y1, pop_trt)
#   ratio <- result$ATT
#   ci <- result$CI
get_agg_ATT_ratio <- function(Mu_trt,Y1,pop_trt) {
  
  # Convert to rates
  #Mu_trt_rate <- 100000 * sweep(Mu_trt, 2, pop_trt, "/" )
  #Y1_rate <- 100000 * Y1/pop_trt
  
  # Get an ATT for each MCMC iteration
  #ATT_vec <- rowMeans((sweep(Mu_trt_rate, 2, Y1_rate, "/"))^(-1)) #ATT= Y1/Y0 = (Y0/Y1)^(-1)
  
  # Get median of MCMC ATT
  #ATT <- median(ATT_vec)
  
  # Get 95% CI
  #CI <- quantile(ATT_vec, probs = c(.025,.975),na.rm = TRUE)
  
  # Sum Y1 over all county time points
  Y1_sum <- sum(Y1)
  
  Y0_sum_vec <- rowSums(Mu_trt)
  
  # ATT ratio vector
  ATT_vec <- Y1_sum/Y0_sum_vec
  
  # Get median of MCMC ATT
  ATT <- median(ATT_vec)
  
  # Get 95% CI
  CI <- quantile(ATT_vec, probs = c(.025,.975),na.rm = TRUE)
  
  return(list(ATT_vec = ATT_vec, ATT = ATT, CI = c(CI)))
}

# Calculate sum of treatment effects across counties and years
# Computes the total number of additional cases due to treatment
#
# INPUTS:
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   Y1: Vector of observed treated outcomes
#   pop_trt: Vector of population denominators (not used in sum calculation)
#
# OUTPUTS:
#   List containing:
#     - TE_vec: Vector of total treatment effect values for each MCMC sample
#     - ATT: Median total treatment effect across MCMC samples
#     - CI: 95% confidence interval (2.5% and 97.5% quantiles)
#
# USAGE:
#   result <- get_agg_sum_TE(Mu_trt, Y1, pop_trt)
#   total_te <- result$ATT
#   ci <- result$CI
get_agg_sum_TE<- function(Mu_trt,Y1,pop_trt) {
  
  
  # Get TE
  TE_mat <- -1*sweep(Mu_trt, 2, Y1, "-") #Y1 - Y0
    #sum(Y1-Mu_trt)
  
  # Get sum TE for each MCMC
  TE_vec <- rowSums(TE_mat)
  
  # Get median of MCMC ATT
  ATT <- median(TE_vec)
  
  # Get 95% CI
  CI <- quantile(TE_vec, probs = c(.025,.975),na.rm = TRUE)
  
  return(list(TE_vec = TE_vec, ATT = ATT, CI = c(CI)))
}

################################################################################
# Table Generation Functions
################################################################################

# Create a table with the total ATT for CA and CT and Overall
# Generates formatted tables for MTGP treatment effect results
#
# INPUTS:
#   data_full_hisp: Dataframe with demographic data including FIPS, CL_CASES, YEAR_DX, C
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   trt_year: Treatment start year (e.g., "1995")
#
# OUTPUTS:
#   df: Dataframe with columns:
#       - Model: Model name ("MTGP")
#       - ATT: Treatment effects for CA, CT, and Overall
#       - 95% CI: Confidence intervals for each group
#
# USAGE:
#   mtgp_table <- ATT_CACT_overall_table(data_full_hisp, Mu_trt, "1995")
ATT_CACT_overall_table <- function(data_full_hisp, Mu_trt, trt_year){
  # Edit 0912: seq((m_trt-1),m) to seq(((m-m_trt)+1),m)
  library(kableExtra)
  
  
  ############
  ## Set-Up ##
  ############
  # get treated counties
  trt_fips <- data_full_hisp %>% filter(C == 1) %>%
    distinct(FIPS) %>% pull(FIPS)
  
  # get Y0_obs
  Y0_obs <- get_Y0_obs(data_full_hisp,trt_fips,trt_year = trt_year)[[1]]
  
  # get Y(1) i.e. observed treated values
  Y1 <- get_Y1_obs1(data_full_hisp) # 13*7=91 treated values; vector of length 91 (county1 time1, county1 time2,...)
  
  # get population matrix
  pop_mat <- get_pop_mat(data_full_hisp)
  
  # get population matrix for treated counties
  n_trt <- length(trt_fips)
  
  #get populations for treated counties
  pop_trt <- as.vector(t(pop_mat[1:n_trt,]))
  
  # get number of treated time points
  m_trt <- sum(is.na(Y0_obs)*1)/n_trt
  
  # get number of time points
  m <- length(unique(data_full_hisp$YEAR_DX))
  
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  
  # Split CA and CT
  
  # Get CA
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  # Find which indices are CA
  ind_CA <- c()
  for(i in 0:(n_trt_CA-1)){
    ind_CA <- append(ind_CA,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  # Subset CA
  Mu_trt_CA <- Mu_trt[,ind_CA]
 
  # Subset Population for CA
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CA <- pop_trt_CA$POP
  
  # Get CT
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(((m-m_trt)+1),m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_CT <- Mu_trt[,ind_CT]
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CT <- pop_trt_CT$POP
  
  # California
  
  agg_ATT_CA_sub <- get_agg_ATT(Mu_trt_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_CA_sub <-  get_agg_ATT(Mu_trt_CA,Y1_CA,pop_trt_CA)$CI
  
 
  # Connecticut
  agg_ATT_CT_sub <- get_agg_ATT(Mu_trt_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_CT_sub <-  get_agg_ATT(Mu_trt_CT,Y1_CT,pop_trt_CT)$CI
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_CA <- CI_fun(agg_bounds_CA_sub)
  CI_CT <- CI_fun(agg_bounds_CT_sub)
  
  
  # Overall ATT
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  agg_ATT_sub <- get_agg_ATT(Mu_trt[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_sub <-  get_agg_ATT(Mu_trt[,ind],Y1,pop_trt[ind])$CI
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI <- CI_fun(agg_bounds_sub)
  
  # Make kable table

  Model <- c("MTGP")
  df <- cbind(Model, 
              round(c(agg_ATT_CA_sub), 2),
              c(CI_CA),
              round(c(agg_ATT_CT_sub), 2),
              c(CI_CT), 
              round(c(agg_ATT_sub), 2),
              c(CI))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI", "ATT", "95% CI")
  
  
  
  return(df)
}

# Create a table with the total treatment effect ratio for CA and CT and Overall
# Generates formatted tables for MTGP treatment effect ratios
#
# INPUTS:
#   data_full_hisp: Dataframe with demographic data including FIPS, CL_CASES, YEAR_DX, C
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   trt_year: Treatment start year (e.g., "1995")
#
# OUTPUTS:
#   df: Dataframe with columns:
#       - Model: Model name ("MTGP")
#       - ATT: Treatment effect ratios for CA, CT, and Overall
#       - 95% CI: Confidence intervals for each group
#
# USAGE:
#   ratio_table <- ATT_CACT_overall_table_ratio(data_full_hisp, Mu_trt, "1995")
ATT_CACT_overall_table_ratio <- function(data_full_hisp, Mu_trt, trt_year){
  # Edit 0912: seq((m_trt-1),m) to seq(((m-m_trt)+1),m)
  library(kableExtra)
  
  
  ############
  ## Set-Up ##
  ############
  # get treated counties
  trt_fips <- data_full_hisp %>% filter(C == 1) %>%
    distinct(FIPS) %>% pull(FIPS)
  
  # get Y0_obs
  Y0_obs <- get_Y0_obs(data_full_hisp,trt_fips,trt_year = trt_year)[[1]]
  
  # get Y(1) i.e. observed treated values
  Y1 <- get_Y1_obs1(data_full_hisp) # 13*7=91 treated values; vector of length 91 (county1 time1, county1 time2,...)
  
  # get population matrix
  pop_mat <- get_pop_mat(data_full_hisp)
  
  # get population matrix for treated counties
  n_trt <- length(trt_fips)
  
  #get populations for treated counties
  pop_trt <- as.vector(t(pop_mat[1:n_trt,]))
  
  # get number of treated time points
  m_trt <- sum(is.na(Y0_obs)*1)/n_trt
  
  # get number of time points
  m <- length(unique(data_full_hisp$YEAR_DX))
  
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  
  # Split CA and CT
  
  # Get CA
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  # Find which indices are CA
  ind_CA <- c()
  for(i in 0:(n_trt_CA-1)){
    ind_CA <- append(ind_CA,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  # Subset CA
  Mu_trt_CA <- Mu_trt[,ind_CA]
  
  # Subset Population for CA
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CA <- pop_trt_CA$POP
  
  # Get CT
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(((m-m_trt)+1),m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_CT <- Mu_trt[,ind_CT]
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CT <- pop_trt_CT$POP
  
  # California
  
  agg_ATT_CA_sub <- get_agg_ATT_ratio(Mu_trt_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_CA_sub <-  get_agg_ATT_ratio(Mu_trt_CA,Y1_CA,pop_trt_CA)$CI
  
  
  # Connecticut
  agg_ATT_CT_sub <- get_agg_ATT_ratio(Mu_trt_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_CT_sub <-  get_agg_ATT_ratio(Mu_trt_CT,Y1_CT,pop_trt_CT)$CI
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_CA <- CI_fun(agg_bounds_CA_sub)
  CI_CT <- CI_fun(agg_bounds_CT_sub)
  
  
  # Overall ATT
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  agg_ATT_sub <- get_agg_ATT_ratio(Mu_trt[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_sub <-  get_agg_ATT_ratio(Mu_trt[,ind],Y1,pop_trt[ind])$CI
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI <- CI_fun(agg_bounds_sub)
  
  # Make kable table
  
  Model <- c("MTGP")
  df <- cbind(Model, 
              round(c(agg_ATT_CA_sub), 2),
              c(CI_CA),
              round(c(agg_ATT_CT_sub), 2),
              c(CI_CT), 
              round(c(agg_ATT_sub), 2),
              c(CI))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI", "ATT", "95% CI")
  
  
  
  return(df)
}

# Create a table of the sum of the treatment effects
# Generates formatted tables for total treatment effects (case counts)
#
# INPUTS:
#   data_full_hisp: Dataframe with demographic data including FIPS, CL_CASES, YEAR_DX, C
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   trt_year: Treatment start year (e.g., "1995")
#
# OUTPUTS:
#   df: Dataframe with columns:
#       - Model: Model name ("MTGP")
#       - ATT: Total treatment effects for CA, CT, and Overall (case counts)
#       - 95% CI: Confidence intervals for each group
#
# USAGE:
#   sum_table <- ATT_CACT_overall_table_sum_TE(data_full_hisp, Mu_trt, "1995")
ATT_CACT_overall_table_sum_TE <- function(data_full_hisp, Mu_trt, trt_year){
  # Edit 0912: seq((m_trt-1),m) to seq(((m-m_trt)+1),m)
  library(kableExtra)
  
  
  ############
  ## Set-Up ##
  ############
  # get treated counties
  trt_fips <- data_full_hisp %>% filter(C == 1) %>%
    distinct(FIPS) %>% pull(FIPS)
  
  # get Y0_obs
  Y0_obs <- get_Y0_obs(data_full_hisp,trt_fips,trt_year = trt_year)[[1]]
  
  # get Y(1) i.e. observed treated values
  Y1 <- get_Y1_obs1(data_full_hisp) # 13*7=91 treated values; vector of length 91 (county1 time1, county1 time2,...)
  
  # get population matrix
  pop_mat <- get_pop_mat(data_full_hisp)
  
  # get population matrix for treated counties
  n_trt <- length(trt_fips)
  
  #get populations for treated counties
  pop_trt <- as.vector(t(pop_mat[1:n_trt,]))
  
  # get number of treated time points
  m_trt <- sum(is.na(Y0_obs)*1)/n_trt
  
  # get number of time points
  m <- length(unique(data_full_hisp$YEAR_DX))
  
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  
  # Split CA and CT
  
  # Get CA
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  # Find which indices are CA
  ind_CA <- c()
  for(i in 0:(n_trt_CA-1)){
    ind_CA <- append(ind_CA,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  # Subset CA
  Mu_trt_CA <- Mu_trt[,ind_CA]
  
  # Subset Population for CA
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CA <- pop_trt_CA$POP
  
  # Get CT
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(((m-m_trt)+1),m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_CT <- Mu_trt[,ind_CT]
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CT <- pop_trt_CT$POP
  
  # California
  
  agg_ATT_CA_sub <- get_agg_sum_TE(Mu_trt_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_CA_sub <-  get_agg_sum_TE(Mu_trt_CA,Y1_CA,pop_trt_CA)$CI
  
  
  # Connecticut
  agg_ATT_CT_sub <- get_agg_sum_TE(Mu_trt_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_CT_sub <-  get_agg_sum_TE(Mu_trt_CT,Y1_CT,pop_trt_CT)$CI
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1])," , ",round(bounds[2]),")")
    return(out)
  }
  
  CI_CA <- CI_fun(agg_bounds_CA_sub)
  CI_CT <- CI_fun(agg_bounds_CT_sub)
  
  
  # Overall ATT
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq(((m-m_trt)+1),m) + (i*m))
  }
  
  agg_ATT_sub <- get_agg_sum_TE(Mu_trt[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_sub <-  get_agg_sum_TE(Mu_trt[,ind],Y1,pop_trt[ind])$CI
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1])," , ",round(bounds[2]),")")
    return(out)
  }
  
  CI <- CI_fun(agg_bounds_sub)
  
  # Make kable table
  
  Model <- c("MTGP")
  df <- cbind(Model, 
              round(c(agg_ATT_CA_sub)),
              c(CI_CA),
              round(c(agg_ATT_CT_sub)),
              c(CI_CT), 
              round(c(agg_ATT_sub)),
              c(CI))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI", "ATT", "95% CI")
  
  
  
  return(df)
}

################################################################################
# Plotting Functions
################################################################################

# Create full ATT plots for CA, CT, and combined states
# Generates time series plots of treatment effects with confidence intervals
#
# INPUTS:
#   data_full_hisp: Dataframe with demographic data including FIPS, CL_CASES, YEAR_DX, C
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   Model: Model name for plot titles (e.g., "MTGP")
#   trt_year: Treatment start year (e.g., "1995")
#   years: Vector of years for x-axis labels
#
# OUTPUTS:
#   plot: Combined ggplot object with three panels (CA, CT, Overall)
#         Each panel shows ATT over time with 95% confidence intervals
#
# USAGE:
#   full_plot <- ATT_full_plot(data_full_hisp, Mu_trt, "MTGP", "1995", years)
ATT_full_plot <- function(data_full_hisp, Mu_trt, Model, trt_year, years){
  library(matrixStats)
  # Full ATT plot for CA
  
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- treated1$CL_CASES
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  ind <- c()
  for(i in 0:(n_trt_CA-1)){
    ind <- append(ind,seq((1),m) + (i*m))
  }
  
  Mu_trt_CA <- Mu_trt[,ind]
  
  # Subset Population
  pop_trt_CA <- treated1$POP
  
  
  # Calculate ATT for all MCMC simulations
  ATT_CA <- get_ATT_full(Mu_trt_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m)
  
  # Summarize ATT
  med <- rowMedians(ATT_CA, na.rm = TRUE)
  bounds <- apply(ATT_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  ## Plot
  model <- rep(Model,m)
  #year <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003")
  year <- years
  LB <-  c(bounds[1,])
  UB <-  c(bounds[2,])
  ATT <- c(med)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  plot_CA <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = "red"))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year), alpha = 0.1)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=trt_year)+
    #scale_colour_manual(ATT,values="blue")+
    #scale_fill_manual("",values="grey12")
    theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+ 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Estimated ATT for California")
  
  plot_CA
  
  # Full ATT plot for CT
  
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  
  Y1_CT <- treated1$CL_CASES
  
  n_trt_CT <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(1,m) + (j*m))
  }
  
  Mu_trt_CT <- Mu_trt[,ind_CT]
  
  # Subset Population
  pop_trt_CT <- treated1$POP
  
  
  # Calculate ATT for all MCMC simulations
  ATT_CT <- get_ATT_full(Mu_trt_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m)
  
  # Summarize ATT
  med <- rowMedians(ATT_CT, na.rm = TRUE)
  bounds <- apply(ATT_CT, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  ## Plot
  model <- rep(Model,m)
  #year <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003")
  year <- years
  LB <-  c(bounds[1,])
  UB <-  c(bounds[2,])
  ATT <- c(med)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  plot_CT <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = "red"))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year), alpha = 0.1)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=trt_year)+
    #scale_colour_manual(ATT,values="blue")+
    #scale_fill_manual("",values="grey12")
    theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+ 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Estimated ATT for Connecticut")
  
  plot_CT
  
  
  # Full ATT plot for both states
  
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")) | (str_starts(FIPS, "09")))
  
  Y1 <- treated1$CL_CASES
  
  n_trt <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  
  
  Mu_trt_CACT <- Mu_trt[,c(ind,ind_CT)]
  
  # Subset Population
  pop_trt <- treated1$POP
  
  
  # Calculate ATT for all MCMC simulations
  ATT_CACT <- get_ATT_full(Mu_trt_CACT,Y1,pop_trt,n_trt,m_trt,m)
  
  # Summarize ATT
  med <- rowMedians(ATT_CACT, na.rm = TRUE)
  bounds <- apply(ATT_CACT, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  ## Plot
  model <- rep(Model,m)
  #year <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003")
  year <- years
  LB <-  c(bounds[1,])
  UB <-  c(bounds[2,])
  ATT <- c(med)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  plot_CACT <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = "red"))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year), alpha = 0.1)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=trt_year)+
    #scale_colour_manual(ATT,values="blue")+
    #scale_fill_manual("",values="grey12")
    theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+ 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Overall Estimated ATT")
  
  plot_CACT
  
  library(cowplot)
  plot <-plot_grid(plot_CA, plot_CT, plot_CACT,nrow = 1)
  
  return(plot)
}

# Create overall ATT plot for combined states
# Generates a single time series plot of treatment effects for all treated counties
#
# INPUTS:
#   data_full_hisp: Dataframe with demographic data including FIPS, CL_CASES, YEAR_DX, C
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   Model: Model name for plot title (e.g., "MTGP")
#   trt_year: Treatment start year (e.g., "1995")
#   years: Vector of years for x-axis labels
#
# OUTPUTS:
#   plot_CACT: Single ggplot object showing overall ATT over time with 95% CI
#
# USAGE:
#   overall_plot <- ATT_overall_plot(data_full_hisp, Mu_trt, "MTGP", "1995", years)
ATT_overall_plot <- function(data_full_hisp, Mu_trt, Model, trt_year, years){
  library(matrixStats)
  
  # Full ATT plot for both states
  
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")) | (str_starts(FIPS, "09")))
  
  Y1 <- treated1$CL_CASES
  
  n_trt <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  
  
  Mu_trt_CACT <- Mu_trt[,c(ind)]
  
  # Subset Population
  pop_trt <- treated1$POP
  
  
  # Calculate ATT for all MCMC simulations
  ATT_CACT <- get_ATT_full(Mu_trt_CACT,Y1,pop_trt,n_trt,m_trt,m)
  
  # Summarize ATT
  med <- rowMedians(ATT_CACT, na.rm = TRUE)
  bounds <- apply(ATT_CACT, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  ## Plot
  model <- rep(Model,m)
  #year <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003")
  year <- years
  LB <-  c(bounds[1,])
  UB <-  c(bounds[2,])
  ATT <- c(med)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  plot_CACT <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = "red"))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year), alpha = 0.1)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=trt_year)+
    #scale_colour_manual(ATT,values="blue")+
    #scale_fill_manual("",values="grey12")
    theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+ 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Overall Estimated ATT")
  
  plot_CACT
  
  return(plot_CACT)
}

# Calculate full ATT matrix for all time points
# Computes treatment effects for each time point across all MCMC samples
#
# INPUTS:
#   Mu_trt: Matrix of posterior predictions [n_mcmc x n_observations]
#   Y1: Vector of observed treated outcomes
#   pop_trt: Vector of population denominators
#   n_trt: Number of treated counties
#   m_trt: Number of treated time points
#   m: Total number of time points
#
# OUTPUTS:
#   ATT_mat: Matrix of ATT values [n_time_points x n_mcmc_samples]
#            Each row represents a time point, each column an MCMC sample
#
# USAGE:
#   att_matrix <- get_ATT_full(Mu_trt, Y1, pop_trt, n_trt, m_trt, m)
get_ATT_full <- function(Mu_trt,Y1,pop_trt,n_trt,m_trt,m) {
  ATT_mat <- matrix(NA,nrow=m, ncol = nrow(Mu_trt))
  Y1_temp <- matrix(Y1, nrow=n_trt, byrow = T) #put y1 in matrix (13x7)
  pop_trt_temp <- matrix(pop_trt, nrow = n_trt, byrow = T)[,1:m] #(13 x 7)
  
  for(i in 1:nrow(Mu_trt)){ #for all MCMC samples
    
    Mu_trt_temp <- Mu_trt[i,] #(1 x 208) 195 = 13 counties x m time points take one iteration
    Mu_trt_temp_mat <- matrix(Mu_trt_temp, nrow = n_trt, byrow = T) #turn into (13 x 15) matrix
    Y0_est <- Mu_trt_temp_mat #estimated counterfactuals (13x7) matrix (trt counties x trt time)
    
    diff <- Y1_temp - Y0_est # observed Y(1) - estimated Y(0)
    rate <- 100000 * (diff / pop_trt_temp) # get rate (13 x 7)
    ATT_mat[,i] <- colMeans(rate)  #get avg of all counties in each time point  (1 x 7) and store into column
  }
  
  return(ATT_mat)
}

#summary(fit)$summary[ind,]
  
  
  
