/*******************************************************************************
 * Multi-Task Gaussian Process Model for Treatment Effect Estimation
 * 
 * This Stan model implements a Multi-Task Gaussian Process (MTGP) approach
 * for estimating treatment effects on lymphoma incidence across demographic
 * subgroups. The model is adapted from Ben-Michael et al.'s multi-outcome
 * setting to a multi-group setting.
 *
 * MODEL STRUCTURE:
 * - Multi-task learning across demographic subgroups (age, sex, race, etc.)
 * - Gaussian Process priors for time trends and unit effects
 * - Poisson likelihood for count data (lymphoma cases)
 * - Treatment effect estimation via counterfactual prediction
 *
 * KEY FEATURES:
 * - Handles multiple demographic outcomes simultaneously
 * - Accounts for temporal correlation via GP time trends
 * - Models unit-specific effects via GP unit effects
 * - Estimates treatment effects for treated units (CA/CT counties after 1995)
 * - Uses only control units for model fitting (synthetic control approach)
 *
 * PARAMETERS:
 * - lengthscale_f: GP lengthscale for time trends
 * - sigma_f: GP scale parameter for time trends
 * - lengthscale_global: GP lengthscale for global effects
 * - sigma_global: GP scale parameter for global effects
 * - z_f: Latent functions for time trends
 * - k_f: Outcome-specific coefficients for time trends
 * - k_d: Unit-specific coefficients
 * - intercepts: Outcome-specific intercepts
 *
 * INPUTS:
 * - N: Number of time points (years)
 * - D: Number of units (counties)
 * - num_outcomes: Number of demographic subgroups
 * - n_k_f: Number of latent functions for time trends
 * - n_k_d: Number of latent functions for unit effects
 * - x: Time covariate (years)
 * - population: Population denominators for each outcome/unit/time
 * - y: Observed case counts for each outcome/unit/time
 * - num_treated: Number of treated units
 * - control_idx: Indices of control units (for fitting)
 *
 * OUTPUTS:
 * - f_mean: Predicted case counts for all units (including treated)
 * - Treatment effects calculated as difference between observed and predicted
 *******************************************************************************/

// Data block: Define all input variables and their constraints
data {
  int<lower=1> N;                    // Number of time points (years: 1988-2003)
  int<lower=1> D;                    // Number of units (counties)
  int<lower=1> num_outcomes;         // Number of demographic outcomes (subgroups)
  int<lower=1> n_k_f;                // Number of latent functions for time trends
  int<lower=1> n_k_d;                // Number of latent functions for unit effects
  vector[N] x;                       // Time covariate (years)
  int<lower=0> population[num_outcomes, N * D];  // Population denominators
  int<lower=0> y[num_outcomes, N * D];           // Observed case counts
  int num_treated;                   // Number of treated units (CA/CT counties)
  int control_idx[N * D - num_treated];          // Indices of control units
}

// Transformed data block: Preprocess and normalize data
transformed data {
  // Normalize time covariate to improve numerical stability
  real xmean = mean(x);
  real xsd = sd(x);
  real xn[N] = to_array_1d((x - xmean)/xsd);
  
  // Add small jitter to covariance matrices for numerical stability
  vector[N] jitter = rep_vector(1e-9, N);
}

// Parameters block: Define all model parameters
parameters {
  // Gaussian Process parameters for time trends
  real<lower=0> lengthscale_f;       // Lengthscale for time trend GP
  real<lower=0> sigma_f;             // Scale parameter for time trend GP

  // Gaussian Process parameters for global effects
  real<lower=0> lengthscale_global;  // Lengthscale for global effect GP
  real<lower=0> sigma_global;        // Scale parameter for global effect GP
  
  // Global effect latent functions for each outcome
  vector[N] z_global[num_outcomes];

  // Multi-task learning parameters
  matrix[N, n_k_f] z_f;              // Latent functions for time trends
  matrix[n_k_f, n_k_d] k_f[num_outcomes];  // Outcome-specific coefficients for time trends
  matrix[n_k_d, D] k_d;              // Unit-specific coefficients
  row_vector[D] intercepts[num_outcomes];   // Outcome-specific intercepts
}

// Model block: Define the likelihood and priors
model {
  // Construct Gaussian Process covariance matrices
  matrix[N, N] K_f = gp_exp_quad_cov(xn, sigma_f, lengthscale_f);
  matrix[N, N] L_f = cholesky_decompose(add_diag(K_f, jitter));
  matrix[N, N] K_global = gp_exp_quad_cov(xn, sigma_global, lengthscale_global);
  matrix[N, N] L_global = cholesky_decompose(add_diag(K_global, jitter));
  
  // Initialize function values for each outcome
  matrix[N, D] f[num_outcomes];
  vector[N] f_global[num_outcomes];

  // Compute function values for each demographic outcome
  for(i in 1:num_outcomes) {
    // f[i] combines time trends and global effects for outcome i
    f[i] = L_f * z_f * k_f[i] * k_d + rep_matrix(L_global * z_global[i], D);
  }
  
  // PRIORS
  // Standard normal priors for latent functions and coefficients
  to_vector(z_f) ~ std_normal();
  for(i in 1:num_outcomes) {
    to_vector(k_f[i]) ~ std_normal();
    z_global[i] ~ std_normal();
    intercepts[i] ~ std_normal();
  }

  // GP hyperparameter priors
  lengthscale_f ~ inv_gamma(5, 5);      // Prior for time trend lengthscale
  sigma_f ~ std_normal();               // Prior for time trend scale
  lengthscale_global ~ inv_gamma(5, 5); // Prior for global effect lengthscale
  sigma_global ~ std_normal();          // Prior for global effect scale
  to_vector(k_d) ~ std_normal();        // Prior for unit-specific coefficients

  // LIKELIHOOD
  // Only use control units for model fitting (synthetic control approach)
  // This allows us to learn the counterfactual trends for treated units
  for(i in 1:num_outcomes) {
    y[i][control_idx] ~ poisson_log(
      log(to_vector(population[i])[control_idx]) + 
      to_vector(rep_matrix(intercepts[i], N) + f[i])[control_idx]
    );
  }
}

// Generated quantities block: Compute predictions for all units
generated quantities {
  matrix[N, D] f_mean[num_outcomes];  // Predicted case counts for all units
  
  {
    // Reconstruct covariance matrices for prediction
    matrix[N, N] K_f = gp_exp_quad_cov(xn, sigma_f, lengthscale_f);
    matrix[N, N] L_f = cholesky_decompose(add_diag(K_f, jitter));
    matrix[N, N] K_global = gp_exp_quad_cov(xn, sigma_global, lengthscale_global);
    matrix[N, N] L_global = cholesky_decompose(add_diag(K_global, jitter));
    
    // Generate predictions for each demographic outcome
    for(i in 1:num_outcomes) {
      // Reshape population data back to matrix format
      matrix[N, D] population_matrix = to_matrix(population[i], N, D);
      
      // Compute predicted case counts (including treated units)
      // This gives us the counterfactual predictions needed for treatment effect estimation
      f_mean[i] = exp(
        log(population_matrix) + 
        rep_matrix(intercepts[i], N) + 
        L_f * z_f * k_f[i] * k_d + 
        rep_matrix(L_global * z_global[i], D)
      );
    }
  }
}
