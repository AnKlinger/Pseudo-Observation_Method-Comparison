# ========================================================================================================================================
# Calculate_Weibull_RMST_Treatment_effect
# ========================================================================================================================================

Calculate_Weibull_RMST_Treatment_effect <- function(data, tau) {
  
  ### This is the Oracle Model!!! ###
  
  # Load necessary libraries
  if (!require("flexsurv")) {install.packages("flexsurv", dependencies = TRUE); library(flexsurv)}
  
  # Fit a Weibull model with shape modeled as function of X1, X4, Treatment
  flex_model <- flexsurvreg(
    Surv(ObservedTime, Status) ~ X1 + X4 + Treatment,
    anc = list(shape = ~ Treatment),
    data = data,
    dist = "weibull"
  )
  
  # Create dataframe to predict on
  mean_profile <- data.frame(
    X1 = 0,
    X4 = 0,
    Treatment = c(0, 1)
  )
  
  # Compute RMST at time tau for each profile
  rmst0 <- summary(flex_model, newdata = mean_profile[1, ], t = tau, type = "rmst")[[1]]$est
  rmst1 <- summary(flex_model, newdata = mean_profile[2, ], t = tau, type = "rmst")[[1]]$est

  return(rmst1 - rmst0)
}







