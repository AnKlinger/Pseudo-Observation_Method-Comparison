# ========================================================================================================================================
# Generate_Data Function
# ========================================================================================================================================

Generate_Data <- function(n = 1000, 
                          baseline_lambda = 0.05,
                          beta_X1 = - 0.3,
                          beta_X4 = 0.45, 
                          beta_treatment = -0.7,
                          treatment_shape_effect = 0,
                          censoring = TRUE,
                          censoring_rate = 0.05,
                          max_censoring_time_unif = 15,
                          p_treat_manual = NULL,
                          seed = NULL) { 
  
  # Load required package
  if (!require("simdata")) install.packages("simdata", dependencies=TRUE)
  
  # Set seed for reproducibility
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  ### v1 should have a strong effect on treatment allocation and outcome ###
  ### v2 should have a strong effect on treatment allocation ###
  ### v3 should have no effect on either treatment allocation or outcome, but should have a strong correlation with the rest ###
  ### v4 should have a strong effect on outcome ###
  
  # Step 1: Generate confounders =========================================================================================================
  
  # Create Correlation Matrix for 4 initial variables z
  correlation_matrix = simdata::cor_from_upper(
    4, 
    rbind(c(1,2,0.1), c(1,3,0.5), c(1,4,0.1),
          c(2,3,0.6), c(2,4,-0.2),
          c(3,4,-0.5)))
  
  # define Transformation for initial variables z into v1,v2,v3,v4
  transformation = simdata::function_list(
    X1 = function(z) z[,1], 
    X2 = function(z) z[,2], 
    X3 = function(z) z[,3],
    X4 = function(z) z[,4] 
  )
  
  # Define simulation design with correlation, transformation, and from a multivariate normal distribution
  sim_design = simdata::simdesign_mvtnorm(
    relations = correlation_matrix,
    transform_initial = transformation,
    prefix_final = NULL
  )
  
  # Simulate dataset sim_data with variables v1,v2,v3,v4
  sim_data = simdata::simulate_data(sim_design, n_obs = n)
  
  
  
  # Step 2: Generate treatment assignment ================================================================================================
  if (is.null(p_treat_manual)) {
    # Default behavior: treatment probability depends on confounders v1 and v2
    logit_p <- -0.5 + 0.5 * sim_data$X1 + 0.8 * sim_data$X2  # Logistic function linear predictor
    p_treat <- 1 / (1 + exp(-logit_p))     # Convert to probability
  } else {
    # Manual setting: all individuals have the same probability of treatment
    p_treat <- rep(p_treat_manual, n)
  }
  
  Treatment <- rbinom(n, 1, p_treat)  # Assign treatment based on probability
  
  # # Step 3: Define shape parameter to be different based on X1 (positive effect on survival)
  # # and X4 (negative effect on survival) =================================================================================================
  # p_weibull <- rep(1, n) # base shape is 1, which is equal to the exponential distribution
  # 
  # # Adjust based on X1 (positive effect), X4 (negative effect) and treatment_shape_effect (positive or negative depends on value)
  # p_weibull <- p_weibull - 0.08 * sim_data$X1 + 0.2 * sim_data$X4 + treatment_shape_effect * Treatment
  
  # # Keep p_weibull within reasonable bounds
  # p_weibull <- pmin(pmax(p_weibull, 0.5), 1.5)  # bound between 0.5 and 1.5 to avoid instability
  
  # Step 3: Adjust shape parameter to be different based on treatment
  p_weibull <- rep(1, n) + treatment_shape_effect * Treatment
  
  # Step 4: Generate Weibull survival times ==============================================================================================
  lin_pred <- beta_X1 * sim_data$X1 + beta_X4 * sim_data$X4 + beta_treatment * Treatment
  scale_param <- exp(lin_pred)  # Convert to scale parameter
  
  # Weibull survival times
  SurvivalTime <- ( (-log(runif(n))) / (baseline_lambda * scale_param) )^(1/p_weibull)
  
  # Step 5: Introduce censoring if censoring TRUE, otherwise no censoring ================================================================
  if (censoring == TRUE) {
    uniform_censoring <- runif(n, min = 0, max = max_censoring_time_unif)
    exponential_censoring <- rexp(n, rate = censoring_rate)
    censoring_time <- 0.7 * uniform_censoring + 0.3 * pmin(exponential_censoring, 30)
    ObservedTime <- pmin(SurvivalTime, censoring_time)  
    Status <- ifelse(SurvivalTime <= censoring_time, 1, 0)  # 1 = event, 0 = censored
  } else {
    ObservedTime <- SurvivalTime  # No censoring, observed time is survival time
    Status <- rep(1, n)           # All individuals experience the event
    censoring_time <- NA          # censoring_time is NA if censoring != TRUE
  }
  
  # Step 6: Create a data frame ==========================================================================================================
  survival_data <- data.frame(
    ObservedTime = ObservedTime,
    Status = Status,
    X1 = sim_data$X1,
    X2 = sim_data$X2,
    X3 = sim_data$X3,
    X4 = sim_data$X4,
    Treatment = Treatment,
    SurvivalTime = SurvivalTime,
    censoring_time = censoring_time
  )
  
  return(survival_data)
}


# Example usage:
#sim_data <- simulate_survival_data(baseline_lambda = 0.07)

# Check first few rows
#head(sim_data)




# define Transformation for initial variables z into v1,v2,v3,v4
# transformation = simdata::function_list(
#   v1 = function(z) floor(10 * z[,1] + 55), 
#   v2 = function(z) exp(0.5 * z[,2] + 1.5), 
#   v3 = function(z) exp(0.4 * z[,3] + 3),
#   v4 = function(z) floor(pmax(0, 5 * exp(z[,4]) - 1)) # was 25 * exp(z[,4]) - 5 in vignette
# )
