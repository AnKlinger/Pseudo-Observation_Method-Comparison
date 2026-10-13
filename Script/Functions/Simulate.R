# ========================================================================================================================================
# Simulate function
# ========================================================================================================================================

Simulate <- function(sample_sizes = 100:1000, 
                     baseline_lambda = 0.05,
                     beta_X1 = -0.3,
                     beta_X4 = 0.45, 
                     beta_treatment = -0.7,
                     treatment_shape_effect = 0,
                     censoring = TRUE,
                     censoring_rate = 0.05,
                     max_censoring_time_unif = 15,
                     p_treat_manual = NULL,
                     seed = NULL,
                     tau = 10,
                     n_boot = 100,
                     n_workers = future::availableCores() - 1) {
  
  # Load required package
  if (!require("future.apply")) install.packages("future.apply", dependencies=TRUE)
  
  # Load external functions
  source("Script/Functions/Wrapper.R")
  source("Script/Functions/Calculate_true_treatment_effect.R")
  
  # Calculate the true treatment effect once
  true_treatment_effect <- Calculate_true_treatment_effect(
    baseline_lambda = baseline_lambda,
    beta_X1 = beta_X1,
    beta_X4 = beta_X4, 
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring = censoring,
    censoring_rate = censoring_rate,
    max_censoring_time_unif = max_censoring_time_unif,
    p_treat_manual = p_treat_manual,
    seed = seed,
    tau = tau
  )
  
  # Plan parallel backend
  future::plan(future::multisession, workers = n_workers)
  
  # Define the wrapper function for future_lapply
  run_wrapper <- function(n_val) {
    result <- Wrapper(
      n = n_val,
      baseline_lambda = baseline_lambda,
      beta_X1 = beta_X1,
      beta_X4 = beta_X4, 
      beta_treatment = beta_treatment,
      treatment_shape_effect = treatment_shape_effect,
      censoring = censoring,
      censoring_rate = censoring_rate,
      max_censoring_time_unif = max_censoring_time_unif,
      p_treat_manual = p_treat_manual,
      seed = seed,
      tau = tau,
      n_boot = n_boot,
      true_treatment_effect = true_treatment_effect
    )
    result$n <- n_val
    return(result)
  }
  
  # Run parallel lapply
  results_list <- future.apply::future_lapply(sample_sizes, run_wrapper)
  
  # Combine all results into a data frame
  results_df <- do.call(rbind, results_list)
  
  # Return results
  return(results_df)
}
