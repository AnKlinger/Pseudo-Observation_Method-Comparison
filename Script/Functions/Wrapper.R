# ========================================================================================================================================
# Wrapper Function
# ========================================================================================================================================

Wrapper <- function(n = 1000, 
                    baseline_lambda = 0.05,
                    beta_X1 = - 0.3,
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
                    true_treatment_effect) {
  
  # Load necessary functions
  source("Script/Functions/Generate_Data.R")
  source("Script/Functions/Performance.R")
  
  # Load necessary libaries
  if (!require("pseudo")) {install.packages("pseudo", dependencies=TRUE); library(pseudo)}
  if (!require("dplyr")) {install.packages("dplyr", dependencies=TRUE); library(dplyr)}
  if (!require("survival")) {install.packages("survival", dependencies=TRUE); library(survival)}
  
  
  # Set seed 
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  
  # Generate the data =====================================================================================================================
  data <- Generate_Data(n = n, 
                        baseline_lambda = baseline_lambda,
                        beta_X1 = beta_X1,
                        beta_X4 = beta_X4, 
                        beta_treatment = beta_treatment,
                        treatment_shape_effect = treatment_shape_effect,
                        censoring = censoring,
                        censoring_rate = censoring_rate,
                        max_censoring_time_unif = max_censoring_time_unif,
                        p_treat_manual = p_treat_manual,
                        seed = seed)
  # Calculate corresponding pseudo-observations:
  data$pseudo <- pseudo::pseudomean(time = data$ObservedTime,
                                    event = data$Status,
                                    tmax = tau)

  
  # Get performance measures ==============================================================================================================
  
  # tryCatch because there can be rare errors in the Performance() function
  result <- tryCatch({
    Performance(data,
                true_treatment_effect,
                n_boot,
                tau)
  }, error = function(e) {
    # Return a data frame with correct structure and NAs
    data.frame(
      model = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"),
      treatment_effect = NA_real_,
      bias = NA_real_,
      MSPE = NA_real_,
      CI_excludes_0 = NA,
      CI_includes_true_value = NA,
      CI_width = NA_real_,
      Prob_ATE_smaller_equal_0 = NA_real_,
      Prob_ATE_smaller_equal_true_treatment_effect = NA_real_,
      true_treatment_effect = true_treatment_effect
    )
  })
  
  return(result)
  
  
}




