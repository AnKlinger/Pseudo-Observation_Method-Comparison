# ========================================================================================================================================
# Performance Function
# ========================================================================================================================================

Performance <- function(data, true_treatment_effect, n_boot, tau) { 
  
  # Load required packages
  if (!require("rms")) {install.packages("rms", dependencies=TRUE); library(rms)}
  if (!require("geepack")) {install.packages("geepack", dependencies=TRUE); library(geepack)}
  if (!require("grf")) {install.packages("grf", dependencies=TRUE); library(grf)}
  if (!require("randomForest")) {install.packages("randomForest", dependencies=TRUE); library(randomForest)}
  
  # Load function for Weibull model
  source("Script/Functions/Calculate_Weibull_RMST_Treatment_effect.R")
  
  # ================================================================================================
  # Helper function for bootstrapping CIs
  # ================================================================================================
  bootstrap_CI <- function(data, model_type) {
    effects <- numeric(n_boot)
    
    for (i in 1:n_boot) {
      boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
      
      # Calculate pseudo observations for the bootstrapped data
      boot_data$pseudo <- pseudo::pseudomean(time = boot_data$ObservedTime,
                                             event = boot_data$Status,
                                             tmax = tau)
      
      if (model_type == "PO") {
        dd <- datadist(boot_data)
        options(datadist = dd)
        model <- orm(pseudo ~ X1 + X2 + X3 + X4 + Treatment, data = boot_data)
        
        no_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
        no_treatment_data$Treatment <- 0
        all_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
        all_treatment_data$Treatment <- 1
        
        pred_diff <- predict(model, type = "mean", newdata = all_treatment_data) -
          predict(model, type = "mean", newdata = no_treatment_data)
        effects[i] <- mean(pred_diff)
        
      } else if (model_type == "PO_splines") {
        dd <- datadist(boot_data)
        options(datadist = dd)
        
        # Wrap everything inside tryCatch to catch model fitting or prediction errors
        effects[i] <- tryCatch({
          # Capture all output, warnings, and messages silently
          invisible(
            capture.output({
              suppressWarnings(
                suppressMessages(
                  model <- orm(pseudo ~ rcs(X1, 3)*Treatment + rcs(X2, 3)*Treatment +
                                 rcs(X3, 3)*Treatment + rcs(X4, 3)*Treatment + Treatment,
                               data = boot_data)
                )
              )
            })
          )
          
          no_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
          no_treatment_data$Treatment <- 0
          all_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
          all_treatment_data$Treatment <- 1
          
          pred_diff <- predict(model, type = "mean", newdata = all_treatment_data) -
            predict(model, type = "mean", newdata = no_treatment_data)
          
          mean(pred_diff)
          
        }, error = function(e) {
          NA
        })
        
        


        
      } else if (model_type == "GEE_splines") {
        boot_data$id <- 1:nrow(boot_data)
        model <- geepack::geese(pseudo ~ rcs(X1, 3)*Treatment + rcs(X2, 3)*Treatment +
                                  rcs(X3, 3)*Treatment + rcs(X4, 3)*Treatment + Treatment,
                                data = boot_data, id = id, jack = TRUE, family = gaussian,
                                corstr = "independence", scale.fix = FALSE)
        
        # Want the average treatment effect (ATE), since splines are not linear, I have to get it via
        # predictions
        no_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
        no_treatment_data$Treatment <- 0
        all_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
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
        
        pred_diff <- preds_all_treatment - preds_no_treatment
        effects[i] <- mean(pred_diff)
      } else if (model_type == "RF") {
        rf <- randomForest::randomForest(pseudo ~ X1 + X2 + X3 + X4 + Treatment, data = boot_data,
                                         mtry = 2)
        
        no_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
        no_treatment_data$Treatment <- 0
        all_treatment_data <- boot_data[, c("X1", "X2", "X3", "X4")]
        all_treatment_data$Treatment <- 1
        
        pred_diff <- predict(rf, newdata = all_treatment_data) - predict(rf, newdata = no_treatment_data)
        effects[i] <- mean(pred_diff)
        
      } else if (model_type == "WB") {
        effects[i] <- Calculate_Weibull_RMST_Treatment_effect(data = boot_data, tau = tau)
      }
    }
    
    CI <- quantile(effects, probs = c(0.025, 0.975), na.rm = TRUE)
    list(CI = CI, effects = effects)
  }
  
  # ================================================================================================
  # Proportional Odds Model
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
  PO_MSPE <- mean((data$pseudo - predict(model, type = "mean", newdata = data))^2)
  PO_CI_excludes_0 <- ifelse(confint(model)["Treatment", 1] >= 0 | confint(model)["Treatment", 2] <= 0, TRUE, FALSE)
  
  # Bootstrap CI
  PO_boot <- bootstrap_CI(data, model_type = "PO")
  PO_CI <- PO_boot$CI
  PO_CI_includes_true_value <- (PO_CI[1] <= true_treatment_effect) & (PO_CI[2] >= true_treatment_effect)
  PO_CI_width <- PO_CI[2] - PO_CI[1]
  
  PO_Prob_ATE_smaller_equal_0 <- sum(PO_boot$effects <= 0) / n_boot
  PO_Prob_ATE_smaller_equal_true_treatment_effect <- sum(PO_boot$effects <= true_treatment_effect) / n_boot
  
  # ================================================================================================
  # Proportional Odds Model with Splines
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
  PO_splines_MSPE <- mean((data$pseudo - predict(model, type = "mean", newdata = data))^2)
  
  # Bootstrap CI
  PO_splines_boot <- bootstrap_CI(data, model_type = "PO_splines")
  PO_splines_CI <- PO_splines_boot$CI
  PO_splines_CI_excludes_0 <- (PO_splines_CI[1] > 0 | PO_splines_CI[2] < 0)
  PO_splines_CI_includes_true_value <- (PO_splines_CI[1] <= true_treatment_effect) &
    (PO_splines_CI[2] >= true_treatment_effect)
  PO_splines_CI_width <- PO_splines_CI[2] - PO_splines_CI[1]
  
  PO_splines_Prob_ATE_smaller_equal_0 <- sum(PO_splines_boot$effects <= 0) / n_boot
  PO_splines_Prob_ATE_smaller_equal_true_treatment_effect <- sum(PO_splines_boot$effects <= true_treatment_effect) / n_boot
  
  # ================================================================================================
  # GEE Model
  # ================================================================================================
  data$id <- 1:nrow(data)
  model <- geepack::geese(pseudo ~ X1 + X2 + X3 + X4 + Treatment,
                          data = data, id = id, jack = TRUE, family = gaussian,
                          corstr = "independence", scale.fix = FALSE)
  
  GEE_Treatment_effect <- model$beta["Treatment"]
  GEE_bias <- GEE_Treatment_effect - true_treatment_effect
  
  X_mat <- model.matrix(~ X1 + X2 + X3 + X4 + Treatment, data = data)
  preds <- as.vector(X_mat %*% model$beta)
  GEE_MSPE <- mean((preds - data$pseudo)^2)
  
  GEE_sd <- summary(model)$mean$san.se[rownames(summary(model)$mean) == "Treatment"]
  GEE_CI_excludes_0 <- ifelse((GEE_Treatment_effect - 1.96 * GEE_sd >= 0) |
                                (GEE_Treatment_effect + 1.96 * GEE_sd <= 0), TRUE, FALSE)
  GEE_CI_includes_true_value <- (GEE_Treatment_effect - 1.96 * GEE_sd <= true_treatment_effect) &
    (GEE_Treatment_effect + 1.96 * GEE_sd >= true_treatment_effect)
  GEE_CI_width <- 1.96 * 2 * GEE_sd
  
  # P(estimated ATE <= 0) can be taken from the normal distribution, since the estimate is
  # normally distributed under assumptions of the GEE model
  GEE_Prob_ATE_smaller_equal_0 <- pnorm(0, mean = GEE_Treatment_effect, sd = GEE_sd)
  # P(estimated ATE <= true ATE) can also be taken from the normal distribution for the same reasons
  GEE_Prob_ATE_smaller_equal_true_treatment_effect <- pnorm(true_treatment_effect,
                                                            mean = GEE_Treatment_effect,
                                                            sd = GEE_sd)
  
  # ================================================================================================
  # GEE Model with Splines
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
  
  # Predict on original data
  X_mat <- model.matrix(~ rcs(X1, 3)*Treatment + rcs(X2, 3)*Treatment +
                          rcs(X3, 3)*Treatment + rcs(X4, 3)*Treatment + Treatment,
                        data = data)
  preds <- as.vector(X_mat %*% model$beta)
  GEE_splines_MSPE <- mean((preds - data$pseudo)^2)
  
  # Bootstrap CI
  GEE_splines_boot <- bootstrap_CI(data, model_type = "GEE_splines")
  GEE_splines_CI <- GEE_splines_boot$CI
  GEE_splines_CI_excludes_0 <- (GEE_splines_CI[1] > 0 | GEE_splines_CI[2] < 0)
  GEE_splines_CI_includes_true_value <- (GEE_splines_CI[1] <= true_treatment_effect) &
    (GEE_splines_CI[2] >= true_treatment_effect)
  GEE_splines_CI_width <- GEE_splines_CI[2] - GEE_splines_CI[1]
  
  GEE_splines_Prob_ATE_smaller_equal_0 <- sum(GEE_splines_boot$effects <= 0) / n_boot
  GEE_splines_Prob_ATE_smaller_equal_true_treatment_effect <- sum(GEE_splines_boot$effects <= true_treatment_effect) / n_boot
  
  
  
  # ================================================================================================
  # Random Forest Model
  # ================================================================================================
  
  rf <- randomForest::randomForest(pseudo ~ X1 + X2 + X3 + X4 + Treatment, data = data,
                                   mtry = 2)

  no_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  no_treatment_data$Treatment <- 0
  all_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
  all_treatment_data$Treatment <- 1

  RF_Treatment_effect <- mean(predict(rf, newdata = all_treatment_data) -
                                predict(rf, newdata = no_treatment_data))

  RF_bias <- RF_Treatment_effect - true_treatment_effect
  RF_MSPE <- mean((data$pseudo - predict(rf, newdata = data))^2)

  # Bootstrap CI
  RF_boot <- bootstrap_CI(data, model_type = "RF")
  RF_CI <- RF_boot$CI
  RF_CI_excludes_0 <- (RF_CI[1] > 0 | RF_CI[2] < 0)
  RF_CI_includes_true_value <- (RF_CI[1] <= true_treatment_effect) & (RF_CI[2] >= true_treatment_effect)
  RF_CI_width <- RF_CI[2] - RF_CI[1]
  
  RF_Prob_ATE_smaller_equal_0 <- sum(RF_boot$effects <= 0) / n_boot
  RF_Prob_ATE_smaller_equal_true_treatment_effect <- sum(RF_boot$effects <= true_treatment_effect) / n_boot
  
  # ================================================================================================
  # Causal Random Forest Model
  # ================================================================================================
  
  data$id <- 1:nrow(data)

  # Fit Causal Random Forest
  Y = data$pseudo
  X = data[, c("X1", "X2", "X3", "X4")]
  W = data$Treatment

  tau.forest <- causal_forest(X, Y, W, num.trees = 4000,
                              mtry = 2) # 4000 recommended in vignette for CI
  ate <- average_treatment_effect(tau.forest)
  # ate [1] is estimate of the average treatment effect
  # ate [2] is the standard error

  CRF_Treatment_effect <- ate [1]
  CRF_bias <- CRF_Treatment_effect - true_treatment_effect
  CRF_MSPE <- NA # causal random forest doesn't predict, also not needed since rf will be better anyway

  CRF_CI <- c(ate [1] - 1.96 * ate [2], ate [1] + 1.96 * ate [2])
  CRF_CI_excludes_0 <- (CRF_CI[1] > 0 | CRF_CI[2] < 0)
  CRF_CI_includes_true_value <- (CRF_CI[1] <= true_treatment_effect) & (CRF_CI[2] >= true_treatment_effect)
  CRF_CI_width <- CRF_CI[2] - CRF_CI[1]
  
  # P(estimated ATE <= 0) can be taken from the normal distribution, since the estimate is
  # normally distributed by Estimation and Inference of Heterogeneous
  # Treatment Effects using Random Forests by Wager 2017 Theorem 1
  CRF_Prob_ATE_smaller_equal_0 <- pnorm(0, mean = ate [1], sd = ate [2])
  # P(estimated ATE <= true ATE) can also be taken from the normal distribution for the same reasons
  CRF_Prob_ATE_smaller_equal_true_treatment_effect <- pnorm(true_treatment_effect,
                                                            mean = ate [1],
                                                            sd = ate [2])
  
  # ================================================================================================
  # Weibull Model
  # ================================================================================================
  
  WB_Treatment_effect <- Calculate_Weibull_RMST_Treatment_effect(data = data, tau = tau)
  
  WB_bias <- WB_Treatment_effect - true_treatment_effect
  WB_MSPE <- NA # the model doesn't predict the pseudo-observations
  
  # Bootstrap CI
  WB_boot <- bootstrap_CI(data, model_type = "WB")
  WB_CI <- WB_boot$CI
  WB_CI_excludes_0 <- (WB_CI[1] > 0 | WB_CI[2] < 0)
  WB_CI_includes_true_value <- (WB_CI[1] <= true_treatment_effect) & (WB_CI[2] >= true_treatment_effect)
  WB_CI_width <- WB_CI[2] - WB_CI[1]
  
  WB_Prob_ATE_smaller_equal_0 <- sum(WB_boot$effects <= 0) / n_boot
  WB_Prob_ATE_smaller_equal_true_treatment_effect <- sum(WB_boot$effects <= true_treatment_effect) / n_boot
  
  # ================================================================================================
  return(data.frame(
    model = rep(c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"), each = length(PO_bias)),
    treatment_effect = c(PO_Treatment_effect, PO_splines_Treatment_effect,
                         GEE_Treatment_effect, GEE_splines_Treatment_effect,
                         RF_Treatment_effect, CRF_Treatment_effect,
                         WB_Treatment_effect),
    bias = c(PO_bias, PO_splines_bias,
             GEE_bias, GEE_splines_bias, 
             RF_bias, CRF_bias,
             WB_bias),
    MSPE = c(PO_MSPE, PO_splines_MSPE,
             GEE_MSPE, GEE_splines_MSPE,
             RF_MSPE, CRF_MSPE,
             WB_MSPE),
    CI_excludes_0 = c(PO_CI_excludes_0, PO_splines_CI_excludes_0,
                      GEE_CI_excludes_0, GEE_splines_CI_excludes_0, 
                      RF_CI_excludes_0, CRF_CI_excludes_0,
                      WB_CI_excludes_0),
    CI_includes_true_value = c(PO_CI_includes_true_value, PO_splines_CI_includes_true_value, 
                               GEE_CI_includes_true_value, GEE_splines_CI_includes_true_value,
                               RF_CI_includes_true_value, CRF_CI_includes_true_value,
                               WB_CI_includes_true_value),
    CI_width = c(PO_CI_width, PO_splines_CI_width, 
                 GEE_CI_width, GEE_splines_CI_width, 
                 RF_CI_width, CRF_CI_width,
                 WB_CI_width),
    Prob_ATE_smaller_equal_0 = c(PO_Prob_ATE_smaller_equal_0, PO_splines_Prob_ATE_smaller_equal_0,
                                 GEE_Prob_ATE_smaller_equal_0, GEE_splines_Prob_ATE_smaller_equal_0,
                                 RF_Prob_ATE_smaller_equal_0, CRF_Prob_ATE_smaller_equal_0,
                                 WB_Prob_ATE_smaller_equal_0),
    Prob_ATE_smaller_equal_true_treatment_effect = c(PO_Prob_ATE_smaller_equal_true_treatment_effect,
                                                     PO_splines_Prob_ATE_smaller_equal_true_treatment_effect,
                                                     GEE_Prob_ATE_smaller_equal_true_treatment_effect,
                                                     GEE_splines_Prob_ATE_smaller_equal_true_treatment_effect,
                                                     RF_Prob_ATE_smaller_equal_true_treatment_effect,
                                                     CRF_Prob_ATE_smaller_equal_true_treatment_effect,
                                                     WB_Prob_ATE_smaller_equal_true_treatment_effect),
    true_treatment_effect = true_treatment_effect
  ))
  
}






# # Causal Forest model ========================================================================================
# Y = data$pseudo
# X = data[, c("X1", "X2", "X3", "X4")]
# W = data$Treatment
# 
# tau.forest <- causal_forest(X, Y, W, num.trees = 2000)
# ate <- average_treatment_effect(tau.forest) # average treatment effect
# 
# # Calculate treatment effect
# CF_Treatment_effect <- ate[1]
# 
# # Calculate bias
# CF_bias <- CF_Treatment_effect - true_treatment_effect
# 
# # Calculate MSE (cannot be done by causal forest model)
# mean( (predict(tau.forest, X) - data$pseudo)^2 )





# # Random Forest
# 
# source("Script/Functions/Generate_Data.R")
# library(grf)
# 
# data <- Generate_Data()
# data$pseudo <- pseudo::pseudomean(time = data$ObservedTime,
#                                   event = data$Status,
#                                   tmax = 10)
# 
# # Code taken and modified from grf vignette at https://grf-labs.github.io/grf/
# 
# # Prepare Data
# Y <- data$pseudo
# X <- data[, c("X1", "X2", "X3", "X4")]
# W <- data$Treatment
# 
# # Add confidence intervals for heterogeneous treatment effects; growing more trees is now recommended.
# tau.forest <- causal_forest(X, Y, W, num.trees = 4000)
# tau.hat <- predict(tau.forest, X, estimate.variance = TRUE)
# mu <- mean(tau.hat$predictions)
# sigma.hat <- sqrt(tau.hat$variance.estimates)
# plot(x = 1:10, y = tau.hat$predictions [seq(1,1000, length.out = 10)],
#      ylim = range(tau.hat$predictions[seq(1,1000, length.out = 10)] +
#                     1.96 * sigma.hat[seq(1,1000, length.out = 10)],
#                   tau.hat$predictions[seq(1,1000, length.out = 10)] -
#                     1.96 * sigma.hat[seq(1,1000, length.out = 10)]),
#      xlab = "x", ylab = "tau", type = "l")
# lines(1:10, tau.hat$predictions[seq(1,1000, length.out = 10)] +
#         1.96 * sigma.hat[seq(1,1000, length.out = 10)], col = 1, lty = 2)
# lines(1:10, tau.hat$predictions[seq(1,1000, length.out = 10)] -
#         1.96 * sigma.hat[seq(1,1000, length.out = 10)], col = 1, lty = 2)
# 
# 
# 
# average_treatment_effect(tau.forest)












# # Random Forest
# source("Script/Functions/Generate_Data.R")
# library(grf)
# 
# # Generate Data
# data <- Generate_Data(n = 5000, beta_treatment = -0.3)
# # add pseudo observations for rmst
# data$pseudo <- pseudo::pseudomean(time = data$ObservedTime,
#                                   event = data$Status,
#                                   tmax = 10)
# # add id to data
# data$id <- 1:nrow(data)
# 
# # Fit Causal Random Forest
# Y = data$ObservedTime
# X = data[, c("X1", "X2", "X3", "X4")]
# W = data$Treatment
# 
# tau.forest <- causal_forest(X, Y, W, num.trees = 2000)
# tau.hat <- predict(tau.forest, estimate.variance = TRUE)
# 
# ate <- average_treatment_effect(tau.forest)
# 
# # ate [1] is estimate of the average treatment effect
# # ate [2] is the standard error
# 
# # Fit GLM as control
# model <- geepack::geese(pseudo ~ X1 + X2 + X3 + X4 + Treatment,
#                         data = data, id = id, jack = TRUE, family=gaussian, 
#                         corstr="independence", scale.fix=FALSE)
# summary(model)
# 
# #MSE
# X_mat <- model.matrix(~ X1 + X2 + X3 + X4 + Treatment, data = data)
# preds <- as.vector(X_mat %*% model$beta)
# 
# mean( ( preds - data$pseudo )^2 )
# 
# 







# rf <- randomForest::randomForest(pseudo ~ X1 + X2 + X3 + X4 + Treatment,
#                            data = data, importance = TRUE)
# predict(rf, newdata = data)
# # Get importance and standard errors
# randomForest::importance(rf, type = 1, scale = TRUE)
# 
# # Calculate treatment effect
# no_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
# no_treatment_data$Treatment <- 0
# 
# all_treatment_data <- data[, c("X1", "X2", "X3", "X4")]
# all_treatment_data$Treatment <- 1
# 
# RF_Treatment_effect <- mean( predict(rf, newdata = all_treatment_data) -
#                                predict(rf, newdata = no_treatment_data) )


