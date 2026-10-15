# ========================================================================================================================================
# Calculate_true_treatment_effect
# ========================================================================================================================================


Calculate_true_treatment_effect <- function(baseline_lambda = 0.05,
                                            beta_X1 = - 0.3,
                                            beta_X4 = 0.45, 
                                            beta_treatment = -0.7,
                                            treatment_shape_effect = 0,
                                            censoring = FALSE,
                                            censoring_rate = 0.05,
                                            max_censoring_time_unif = 15,
                                            p_treat_manual = NULL,
                                            seed = NULL,
                                            tau = 10) {
  
  # Load necessary functions and packages
  source("Script/Functions/Generate_Data.R")
  if (!require("survival")) {install.packages("survival", dependencies=TRUE); library(survival)}
  
  set.seed(1)
  data_all_treatment <- Generate_Data(n = 100000, 
                                      baseline_lambda = baseline_lambda,
                                      beta_X1 = beta_X1,
                                      beta_X4 = beta_X4, 
                                      beta_treatment = beta_treatment,
                                      treatment_shape_effect = treatment_shape_effect,
                                      censoring = censoring,
                                      censoring_rate = censoring_rate,
                                      max_censoring_time_unif = max_censoring_time_unif,
                                      p_treat_manual = 1)
  
  set.seed(1)
  data_no_treatment <- Generate_Data(n = 100000, 
                                     baseline_lambda = baseline_lambda,
                                     beta_X1 = beta_X1,
                                     beta_X4 = beta_X4, 
                                     beta_treatment = beta_treatment,
                                     treatment_shape_effect = treatment_shape_effect,
                                     censoring = censoring,
                                     censoring_rate = censoring_rate,
                                     max_censoring_time_unif = max_censoring_time_unif,
                                     p_treat_manual = 0)
  
  set.seed(rnorm(1)) # reset seed the Simulate() function doesn't get the same data every time
  
  # Compute the Kaplan-Meier curves
  km_fit_all_treatment <- survfit(Surv(ObservedTime, Status) ~ 1, data = data_all_treatment)  # One curve for all individuals
  km_fit_no_treatment <- survfit(Surv(ObservedTime, Status) ~ 1, data = data_no_treatment)  # One curve for all individuals
  
  # Extract RMST at time tau
  RMST_all_treatment <- survival:::survmean(km_fit_all_treatment, rmean=tau)[[1]]["rmean"]
  RMST_no_treatment <- survival:::survmean(km_fit_no_treatment, rmean=tau)[[1]]["rmean"] 
  
  return(unname(RMST_all_treatment - RMST_no_treatment)) # remove row name "rmean"
  
}