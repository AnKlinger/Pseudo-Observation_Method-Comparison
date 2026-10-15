# ========================================================================================================================================
# bias_big_dataset Function
# ========================================================================================================================================

bias_big_dataset <- function(beta_treatment,
                             treatment_shape_effect,
                             censoring_rate,
                             max_censoring_time_unif) {
  
  # Load required packages
  if (!require("rms")) {install.packages("rms", dependencies=TRUE); library(rms)}
  if (!require("geepack")) {install.packages("geepack", dependencies=TRUE); library(geepack)}
  if (!require("grf")) {install.packages("grf", dependencies=TRUE); library(grf)}
  if (!require("randomForest")) {install.packages("randomForest", dependencies=TRUE); library(randomForest)}
  
  # Load necessary functions
  source("Script/Functions/Generate_Data.R")
  source("Script/Functions/Calculate_Weibull_RMST_Treatment_effect.R")
  source("Script/Functions/Calculate_true_treatment_effect.R")
  
  # Generate the data
  data <- Generate_Data(n = 1500, 
                        beta_treatment = beta_treatment,
                        treatment_shape_effect = treatment_shape_effect,
                        censoring_rate = censoring_rate,
                        max_censoring_time_unif = max_censoring_time_unif)
  # Calculate corresponding pseudo-observations:
  data$pseudo <- pseudo::pseudomean(time = data$ObservedTime,
                                    event = data$Status,
                                    tmax = 10)
  
  # Calculate true treatment effect
  true_treatment_effect <- Calculate_true_treatment_effect( 
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring_rate = censoring_rate,
    max_censoring_time_unif = max_censoring_time_unif,
    tau = 10
  )
  
  
  
  # ================================================================================================
  ##### Proportional Odds Model #####
  # ================================================================================================
  dd <- datadist(data)
  options(datadist = dd)
  model <- orm(pseudo ~ X1 + X2 + X3 + X4 + Treatment, data = data)
  
  no_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  no_treatment_data$Treatment <- 0
  all_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  all_treatment_data$Treatment <- 1
  
  PO_Treatment_effect <- mean(predict(model, type = "mean", newdata = all_treatment_data) -
                                predict(model, type = "mean", newdata = no_treatment_data))
  
  PO_bias <- PO_Treatment_effect - true_treatment_effect
  
  # ================================================================================================
  ##### Proportional Odds Model with Splines #####
  # ================================================================================================
  dd <- datadist(data)
  options(datadist = dd)
  model <- orm(pseudo ~ rcs(X1, 3)*Treatment + rcs(X2, 3)*Treatment +
                 rcs(X3, 3)*Treatment + rcs(X4, 3)*Treatment + Treatment,
               data = data)
  
  no_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  no_treatment_data$Treatment <- 0
  all_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  all_treatment_data$Treatment <- 1
  
  PO_splines_Treatment_effect <- mean(predict(model, type = "mean", newdata = all_treatment_data) -
                                        predict(model, type = "mean", newdata = no_treatment_data))
  
  PO_splines_bias <- PO_splines_Treatment_effect - true_treatment_effect
  
  # ================================================================================================
  ##### GEE Model #####
  # ================================================================================================
  data$id <- 1:nrow(data)
  model <- geepack::geese(pseudo ~ X1 + X2 + X3 + X4 + Treatment,
                          data = data, id = id, jack = TRUE, family = gaussian,
                          corstr = "independence", scale.fix = FALSE)
  
  GEE_Treatment_effect <- model$beta["Treatment"]
  GEE_bias <- GEE_Treatment_effect - true_treatment_effect
  
  # ================================================================================================
  ##### GEE Model with Splines #####
  # ================================================================================================
  
  data$id <- 1:nrow(data)
  model <- geepack::geese(pseudo ~ rcs(X1, 3)*Treatment + rcs(X2, 3)*Treatment +
                            rcs(X3, 3)*Treatment + rcs(X4, 3)*Treatment + Treatment,
                          data = data, id = id, jack = TRUE, family = gaussian,
                          corstr = "independence", scale.fix = FALSE)
  
  # Want the average treatment effect (ATE), since splines are not linear, I have to get it via
  # predictions
  no_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  no_treatment_data$Treatment <- 0
  all_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  all_treatment_data$Treatment <- 1
  
  # Predict on no_treatment data
  X_mat_no_treatment <- model.matrix(~ rcs(X1, 3)*Treatment + rcs(X2, 3)*Treatment +
                                       rcs(X3, 3)*Treatment + rcs(X4, 3)*Treatment + Treatment,
                                     data = no_treatment_data)
  preds_no_treatment <- as.vector(X_mat_no_treatment %*% model$beta)
  
  # Predict on all_treatment_data
  X_mat_all_treatment <- model.matrix(~ rcs(X1, 3)*Treatment + rcs(X2, 3)*Treatment +
                                        rcs(X3, 3)*Treatment + rcs(X4, 3)*Treatment + Treatment,
                                      data = all_treatment_data)
  preds_all_treatment <- as.vector(X_mat_all_treatment %*% model$beta)
  
  GEE_splines_Treatment_effect <- mean(preds_all_treatment - preds_no_treatment)
  GEE_splines_bias <- GEE_splines_Treatment_effect - true_treatment_effect
  
  # ================================================================================================
  ##### Random Forest Model #####
  # ================================================================================================
  
  # data$id <- 1:nrow(data)
  # 
  # # Fit Causal Random Forest
  # Y = data$pseudo
  # X = data[, c("X1", "X2", "X3", "X4")]
  # W = data$Treatment
  # 
  # tau.forest <- causal_forest(X, Y, W, num.trees = 4000) # 4000 recommended in vignette for CI
  # ate <- average_treatment_effect(tau.forest)
  # # ate [1] is estimate of the average treatment effect
  # # ate [2] is the standard error
  # 
  # RF_Treatment_effect <- ate [1]
  # RF_bias <- RF_Treatment_effect - true_treatment_effect
  
  
  rf <- randomForest::randomForest(pseudo ~ X1 + X2 + X3 + X4 + Treatment, data = data)
  
  no_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  no_treatment_data$Treatment <- 0
  all_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  all_treatment_data$Treatment <- 1
  
  RF_Treatment_effect <- mean(predict(rf, newdata = all_treatment_data) -
                                predict(rf, newdata = no_treatment_data))
  
  RF_bias <- RF_Treatment_effect - true_treatment_effect
  
  # ================================================================================================
  # Causal Random Forest Model
  # ================================================================================================
  
  data$id <- 1:nrow(data)
  
  # Fit Causal Random Forest
  Y = data$pseudo
  X = data[, c("X1", "X2", "X3", "X4")]
  W = data$Treatment
  
  tau.forest <- causal_forest(X, Y, W, num.trees = 4000) # 4000 recommended in vignette for CI
  ate <- average_treatment_effect(tau.forest)
  # ate [1] is estimate of the average treatment effect
  # ate [2] is the standard error
  
  CRF_Treatment_effect <- ate [1]
  CRF_bias <- CRF_Treatment_effect - true_treatment_effect
  
  
  # ================================================================================================
  ##### Weibull Model #####
  # ================================================================================================
  
  WB_Treatment_effect <- Calculate_Weibull_RMST_Treatment_effect(data = data, tau = 10)
  
  WB_bias <- WB_Treatment_effect - true_treatment_effect
  
  
  return(data.frame(true_treatment_effect = true_treatment_effect,
                    PO_bias = PO_bias,
                    PO_splines_bias = PO_splines_bias,
                    GEE_bias = GEE_bias,
                    GEE_splines_bias = GEE_splines_bias,
                    RF_bias = RF_bias,
                    WB_bias = WB_bias))
  
}