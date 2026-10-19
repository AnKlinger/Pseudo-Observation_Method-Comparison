# ========================================================================================================================================
##### Main #####
# ========================================================================================================================================


# ========================================================================================================================================
##### Simulation #####
# ========================================================================================================================================


rm(list = ls())

source("Script/Functions/Simulate.R")

n <- 4

# Calculate uniformly distributed sample sizes on the log2 scale
uniform_sequence_on_log_scale <- seq(from = log(100.1, base = 2),
                                     to = log(1050, base = 2),
                                     length.out = n)
# Calculate real sample sizes for simulation
sample_sizes <- round(2^uniform_sequence_on_log_scale)
# Add anchor points for point estimates and spline knots
sample_sizes <- c(sample_sizes,
                  rep(125, max(round(n/10), 1)),
                  rep(250, max(round(n/10), 1)),
                  rep(500, max(round(n/10), 1)),
                  rep(1000, max(round(n/10), 1)))

# Define scenarios with fixed scenario numbers
scenarios <- list(
  "1" = list(beta_treatment = -0.5, treatment_shape_effect = 0),
  "2" = list(beta_treatment = 0, treatment_shape_effect = 0),
  "3" = list(beta_treatment = 0.5, treatment_shape_effect = -0.3),
  "4" = list(beta_treatment = -1, treatment_shape_effect = 0.49)
)

# Optional: Comment out scenarios you don't want to run
# scenarios["1"] <- NULL  # Uncomment this line to skip scenario 1
scenarios["2"] <- NULL  # Uncomment this line to skip scenario 2
scenarios["3"] <- NULL  # Uncomment this line to skip scenario 3
scenarios["4"] <- NULL  # Uncomment this line to skip scenario 4

# Define censoring settings
censoring_settings <- list(
  list(max_censoring_time_unif = 20, censoring_rate = 0.02, label = 20),
  list(max_censoring_time_unif = 20, censoring_rate = 0.06, label = 40),
  list(max_censoring_time_unif = 15, censoring_rate = 0.08, label = 60)
)

start <- Sys.time()

# Loop over scenarios and censoring settings
for (scenario_number in names(scenarios)) {
  scenario <- scenarios[[scenario_number]]
  
  for (censoring in censoring_settings) {
    
    # Track start time
    start_time <- Sys.time()
    
    # Run simulation
    results_df <- Simulate(
      sample_sizes = sample_sizes,
      beta_treatment = scenario$beta_treatment,
      treatment_shape_effect = scenario$treatment_shape_effect,
      censoring_rate = censoring$censoring_rate,
      max_censoring_time_unif = censoring$max_censoring_time_unif
    )
    
    # Track end time and calculate duration
    end_time <- Sys.time()
    simulation_duration <- as.numeric(difftime(end_time, start_time, units = "hours"))
    
    # Add time column to results_df
    results_df$simulation_time_hours <- simulation_duration
    
    # Create filename and save
    # filename <- paste0("Output/Scenario_", scenario_number, 
    #                    "_censoring_", censoring$label, ".RData")
    filename <- paste0("Output/Scenario_", scenario_number, 
                       "_censoring_", censoring$label, "_test", ".RData")
    save(results_df, file = filename)
  }
}

end <- Sys.time()
end-start











# ========================================================================================================================================
##### Visualization of Scenarios with Survival Curves #####
# ========================================================================================================================================

library(survival)
source("Script/Functions/Generate_Data.R")

simulate_and_plot <- function(beta_treatment, treatment_shape_effect, Scenario) {
  
  set.seed(1)
  data_all_treatment <- Generate_Data(
    n = 20000,
    p_treat_manual = 1,
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring = FALSE
  )
  
  set.seed(1)
  data_no_treatment <- Generate_Data(
    n = 20000,
    p_treat_manual = 0,
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring = FALSE
  )
  
  data <- rbind(data_all_treatment, data_no_treatment)
  
  fit_true <- survfit(Surv(ObservedTime, Status) ~ Treatment, data = data)
  
  data_simulated_normally <- Generate_Data(
    n = 20000,
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring = FALSE
  )
  
  fit_descriptive <- survfit(Surv(ObservedTime, Status) ~ Treatment, data = data_simulated_normally)
  
  # Kaplan Meier Plot
  svg(paste0("Output/Plots/Scenario Visualization/Kaplan_Meier_", Scenario, ".svg"),
      width = 5, height = 3.75)
  
  par(mar = c(4, 4, 4, 1) + 0.1) # default is c(5, 4, 4, 2) + 0.1
  
  plot(fit_true, col = c("red", "blue"), lty = c(1,1), lwd = 1.5, xlab = "Time",
       xlim = c(0, 11), ylim = c(0.3, 1))
       #main = paste0("Survival Curves of Treatment (", "Population ", Scenario, ")"))
  
  lines(fit_descriptive, col = c("red", "blue"), lty = c(2,2))
  
  segments(x0 = 10, y0 = 0.3, x1 = 10, y1 = 1, lty = 2, col = "black")
  
  legend("bottomleft",
         legend = c("Control (Descriptive)", "Treatment (Descriptive)",
                    "Control (True)", "Treatment (True)"),
         col = c("red", "blue", "red", "blue"),
         lty = c(2, 2, 1, 1), bty = "n")

  dev.off()
}

# Define parameter combinations
params <- data.frame(
  beta_treatment = c(-0.5, 0, 0.5, -1),
  treatment_shape_effect = c(0, 0, -0.3, 0.49),
  Scenario = c(1, 2, 3, 4)
)

# Loop over each row
for (i in 1:nrow(params)) {
  simulate_and_plot(
    beta_treatment = params$beta_treatment[i],
    treatment_shape_effect = params$treatment_shape_effect[i],
    Scenario = params$Scenario[i]
  )
}










# ========================================================================================================================================
##### Calculate Raw RMST for each population #####
# ========================================================================================================================================

source("Script/Functions/Generate_Data.R")
library(survival)

scenarios <- list(
  "1" = list(beta_treatment = -0.5, treatment_shape_effect = 0),
  "2" = list(beta_treatment = 0, treatment_shape_effect = 0),
  "3" = list(beta_treatment = 0.5, treatment_shape_effect = -0.3),
  "4" = list(beta_treatment = -1, treatment_shape_effect = 0.49)
)
tau = 10

for(i in 1:4) {
  set.seed(1)
  scenario <- scenarios[[i]]
  data <- Generate_Data(n = 100000,
                        beta_treatment = scenario$beta_treatment,
                        treatment_shape_effect = scenario$treatment_shape_effect,
                        censoring = FALSE
  ) 
  
  
  # Compute the Kaplan-Meier curves
  km_fit<- survfit(Surv(ObservedTime, Status) ~ Treatment, data = data)  # One curve for all individuals
  
  # Extract RMST at time tau
  RMST <- survival:::survmean(km_fit, rmean=tau)$matrix[1:2, "rmean"]
  print(paste0("Population ", i, ": ", RMST["Treatment=1"] - RMST["Treatment=0"]))
}









# # ========================================================================================================================================
# ##### Visualize Proportion of Successes in Binary Variables #####
# # ========================================================================================================================================
# 
# 
# rm(list = ls())
# 
# load("Output/Scenario_4_censoring_60.RData")
# 
# library(tidyverse)
# 
# GEE_results_df <- results_df %>% filter(model == "GEE")
# 
# plot(y = GEE_results_df$CI_includes_true_value, x = GEE_results_df$n)
# 
# # Create 20 equal-width bins of 'n'
# GEE_results_df <- GEE_results_df %>%
#   mutate(n_bin = cut(n, breaks = 20))
# 
# # Calculate proportion of CI_includes_true_value == 1 in each bin
# prop_by_bin <- GEE_results_df %>%
#   group_by(n_bin) %>%
#   summarise(
#     count = n(),
#     num_success = sum(CI_includes_true_value == 1),
#     proportion = num_success / count
#   )
# 
# # Plot the proportions
# ggplot(prop_by_bin, aes(x = n_bin, y = proportion)) +
#   geom_col(fill = "steelblue") +
#   labs(x = "n (binned into 20 intervals)", 
#        y = "Proportion where CI includes true value", 
#        title = "Proportion of CI_includes_true_value == 1 by n bin") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))






# ========================================================================================================================================
##### Visualization of Pseudo-Observations #####
# ========================================================================================================================================

library(survival)
library(tidyverse)
library(patchwork)
source("Script/Functions/Generate_Data.R")

set.seed(101)
low_censoring <- Generate_Data(beta_treatment = -0.5,
                               treatment_shape_effect = 0,
                               max_censoring_time_unif = 20,
                               censoring_rate = 0.02)
low_censoring$pseudo <- pseudo::pseudomean(time = low_censoring$ObservedTime,
                                           event = low_censoring$Status,
                                           tmax = 10)

set.seed(101)
mid_censoring <- Generate_Data(beta_treatment = -0.5,
                               treatment_shape_effect = 0,
                               max_censoring_time_unif = 20,
                               censoring_rate = 0.06)
mid_censoring$pseudo <- pseudo::pseudomean(time = mid_censoring$ObservedTime,
                                           event = mid_censoring$Status,
                                           tmax = 10)

set.seed(101)
high_censoring <- Generate_Data(beta_treatment = -0.5,
                                treatment_shape_effect = 0,
                                max_censoring_time_unif = 15,
                                censoring_rate = 0.08)
high_censoring$pseudo <- pseudo::pseudomean(time = high_censoring$ObservedTime,
                                           event = high_censoring$Status,
                                           tmax = 10)

# Create a histogram of the pseudo-observations
p1 <- ggplot(low_censoring, aes(x = pseudo)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7, binwidth = 0.01) +
  labs(title = "RMST(10) Pseudo-Observations\nwith 20% Censoring",
       x = "Value of Pseudo-Observation",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(limits = c(0,11))
p2 <- ggplot(mid_censoring, aes(x = pseudo)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7, binwidth = 0.01) +
  labs(title = "RMST(10) Pseudo-Observations\nwith 40% Censoring",
       x = "Value of Pseudo-Observation",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(limits = c(0,11))
p3 <- ggplot(high_censoring, aes(x = pseudo)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7, binwidth = 0.01) +
  labs(title = "RMST(10) Pseudo-Observations\nwith 60% Censoring",
       x = "Value of Pseudo-Observation",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(limits = c(0,11))

# Save plots as svg
ggsave(filename = "Output/Plots/PO examples/PO_RMST_low_censoring.svg",
       plot = p1,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
ggsave(filename = "Output/Plots/PO examples/PO_RMST_high_censoring.svg",
       plot = p3,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6












# # ========================================================================================================================================
# ##### Bias of Models #####
# # ========================================================================================================================================
# 
# 
# # To know if models are biased, use them on a big dataset
# 
# 
# rm(list = ls())
# 
# source("Script/Functions/bias_big_dataset.R")
# 
# 
# ##### Bias Simulation #####
# 
# set.seed(3932)
# 
# # Empty list to store all results
# results_list <- list()
# 
# scenarios <- list(
#   "1" = list(beta_treatment = -0.5, treatment_shape_effect = 0),
#   "2" = list(beta_treatment = 0, treatment_shape_effect = 0),
#   "3" = list(beta_treatment = 0.5, treatment_shape_effect = -0.3),
#   "4" = list(beta_treatment = -1, treatment_shape_effect = 0.49)
# )
# 
# # Define censoring settings
# censoring_settings <- list(
#   list(max_censoring_time_unif = 20, censoring_rate = 0.02, label = 20),
#   list(max_censoring_time_unif = 20, censoring_rate = 0.06, label = 40),
#   list(max_censoring_time_unif = 15, censoring_rate = 0.08, label = 60)
# )
# 
# # Optional: Comment out scenarios you don't want to run
# # scenarios["1"] <- NULL  # Uncomment this line to skip scenario 1
# # scenarios["2"] <- NULL  # Uncomment this line to skip scenario 2
# # scenarios["3"] <- NULL  # Uncomment this line to skip scenario 3
# # scenarios["4"] <- NULL  # Uncomment this line to skip scenario 4
# 
# library(future)
# library(future.apply)
# 
# # Set parallel plan
# plan(multisession, workers = availableCores() - 1)  # or multicore on Linux/Mac, multisession is cross-platform
# 
# # Create list of combinations to iterate over
# param_grid <- expand.grid(
#   scenario_number = names(scenarios),
#   censoring_idx = seq_along(censoring_settings),
#   replicate = 1:1000
# )
# 
# start <- Sys.time()
# 
# # Use future_lapply to parallelize the loop
# results_list <- future_lapply(1:nrow(param_grid), function(i) {
#   row <- param_grid[i, ]
#   
#   scenario_number <- as.character(row$scenario_number)
#   scenario <- scenarios[[scenario_number]]
#   censoring <- censoring_settings[[row$censoring_idx]]
#   replicate <- row$replicate
#   
#   result_df <- bias_big_dataset(
#     beta_treatment = scenario$beta_treatment,
#     treatment_shape_effect = scenario$treatment_shape_effect,
#     censoring_rate = censoring$censoring_rate,
#     max_censoring_time_unif = censoring$max_censoring_time_unif
#   )
#   
#   # Add identifying columns
#   result_df$Scenario <- paste0("Scenario_", scenario_number, "_censoring_", censoring$label)
#   result_df$Replicate <- replicate
#   
#   return(result_df)
# })
# 
# end <- Sys.time()
# end - start
# 
# 
# # Combine all rows into one final data frame
# df <- do.call(rbind, results_list)
# 
# save(df, file = "Output/model_bias.RData")
# 
# 
# 
# 
# # ========================================================================================================================================
# ##### Visualization of Bias #####
# # ========================================================================================================================================
# 
# rm(list = ls())
# load("Output/model_bias.RData")
# library(tidyverse)
# 
# # Average bias over replicates
# df <- df %>%
#   group_by(Scenario) %>%
#   summarise(
#     Scenario = first(Scenario),
#     true_treatment_effect = first(true_treatment_effect),
#     PO_bias = mean(PO_bias, na.rm = TRUE),
#     PO_splines_bias = mean(PO_splines_bias, na.rm = TRUE),
#     GEE_bias = mean(GEE_bias, na.rm = TRUE),
#     GEE_splines_bias = mean(GEE_splines_bias, na.rm = TRUE),
#     RF_bias = mean(RF_bias, na.rm = TRUE),
#     WB_bias = mean(WB_bias, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # Reshape to long format
# df_long <- df %>%
#   pivot_longer(
#     cols = ends_with("_bias"),
#     names_to = "model",
#     values_to = "bias"
#   )
# 
# # Ensure Scenario is a factor with levels in desired top-to-bottom order
# df_long <- df_long %>%
#   mutate(Scenario = factor(Scenario, levels = rev(unique(Scenario))))
# 
# # Extract true_treatment_effect per scenario (assuming it's constant within scenario)
# treat_effect_labels <- df_long %>%
#   group_by(Scenario) %>%
#   summarise(true_treatment_effect = unique(true_treatment_effect)) %>%
#   arrange(Scenario)
# 
# # Create the plot
# model_bias_plot <- ggplot(df_long, aes(x = bias, y = Scenario, color = model, shape = model)) +
#   geom_point(size = 4, position = position_dodge(width = 0.7)) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   labs(
#     title = "Bias of Models per Scenario",
#     x = "Bias",
#     y = "Scenario"
#   ) +
#   scale_y_discrete(
#     labels = paste0(treat_effect_labels$Scenario, " (True TE: ", round(treat_effect_labels$true_treatment_effect, 3), ")")
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(size = 10),
#     legend.title = element_blank()
#   )
# 
# save(model_bias_plot, file = "Output/model_bias_plot.RData")
  









# ========================================================================================================================================
##### Visualization of Performance Measures #####
# ========================================================================================================================================

# Create all Plots ============================================================================

library(tidyverse)
library(patchwork)

source("Script/Functions/plot_continuous.R")
source("Script/Functions/plot_logistic.R")
source("Script/Functions/create_plots.R")

load("Output/Scenario_1_censoring_20.RData")
all_plots_1_20 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.45,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0, 0.5),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_1_censoring_40.RData")
all_plots_1_40 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.45,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0, 0.5),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_1_censoring_60.RData")
all_plots_1_60 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.5,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0, 0.5),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_2_censoring_20.RData")
all_plots_2_20 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0.2, 0.8),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_2_censoring_40.RData")
all_plots_2_40 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0.2, 0.8),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_2_censoring_60.RData")
all_plots_2_60 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0.2, 0.8),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_3_censoring_20.RData")
all_plots_3_20 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0.15, 0.75),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_3_censoring_40.RData")
all_plots_3_40 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0.15, 0.75),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_3_censoring_60.RData")
all_plots_3_60 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0.15, 0.75),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_4_censoring_20.RData")
all_plots_4_20 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0, 0.7),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_4_censoring_40.RData")
all_plots_4_40 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0, 0.7),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

load("Output/Scenario_4_censoring_60.RData")
all_plots_4_60 <- create_plots(data = results_df,
                               treatment_effect_added_y_limits = 0.35,
                               CI_width_y_limits = c(0, 3.2),
                               Prob_ATE_smaller_equal_0_y_limits = c(0, 0.7),
                               Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1))

# =============================================================================================



# RMSE ========================================================================================
# Loop through the scenarios 1 to 4
for (i in 1:4) {
  # Loop through the censoring levels: 20, 40, 60
  for (j in c(20, 40, 60)) {
    
    # Construct the object name as a string
    obj_name <- paste0("all_plots_", i, "_", j)
    
    # Get object
    obj <- get(obj_name)
    
    p <- obj$RMSE_plot[[1]] + ggtitle(paste("Population", i, "Censoring", j, "%"))
    
    ggsave(filename = paste0("Output/Plots/RMSE/RMSE_", "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
    
  }
}
# =============================================================================================




# Estimated ATE ===============================================================================
# Loop through the scenarios 1 to 4
for (i in 1:4) {
  # Loop through the censoring levels: 20, 40, 60
  for (j in c(20, 40, 60)) {
    
    # Construct the object name as a string
    obj_name <- paste0("all_plots_", i, "_", j)
    
    # Get object
    obj <- get(obj_name)
    
    p <- obj$treatment_effect[[1]] + plot_annotation(title = paste("Population", i, "Censoring", j, "%")) 
    
    ggsave(filename = paste0("Output/Plots/Estimated_ATE/Estimated_ATE_grid_", "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 5) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
    
    p <- obj$treatment_effect[[2]] + ggtitle(paste("Population", i, "Censoring", j, "%"))
    
    ggsave(filename = paste0("Output/Plots/Estimated_ATE/Estimated_ATE_", "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
    
    
  }
}    
# =============================================================================================




# CI Width ====================================================================================
# Loop through the scenarios 1 to 4
for (i in 1:4) {
  # Loop through the censoring levels: 20, 40, 60
  for (j in c(20, 40, 60)) {
    
    # Construct the object name as a string
    obj_name <- paste0("all_plots_", i, "_", j)
    
    # Get object
    obj <- get(obj_name)
    
    p <- obj$CI_width[[1]] + plot_annotation(title = paste("Population", i, "Censoring", j, "%"))
    
    ggsave(filename = paste0("Output/Plots/CI_Width/CI_Width_grid_", "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 5) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
    

    p <- obj$CI_width[[2]] + ggtitle(paste("Population", i, "Censoring", j, "%"))
    
    ggsave(filename = paste0("Output/Plots/CI_Width/CI_Width_", "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

  }
} 
# =============================================================================================




# Prob ATE smaller equal 0 ====================================================================
# Loop through the scenarios 1 to 4
for (i in 1:4) {
  # Loop through the censoring levels: 20, 40, 60
  for (j in c(20, 40, 60)) {
    
    # Construct the object name as a string
    obj_name <- paste0("all_plots_", i, "_", j)
    
    # Get object
    obj <- get(obj_name)

    p <- obj$Prob_ATE_smaller_equal_0[[1]] + plot_annotation(title = paste("Population", i, "Censoring", j, "%"))

    ggsave(filename = paste0("Output/Plots/Prob_smaller_0/Prob_ATE_smaller_equal_0_grid_",
                             "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 5) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
    
    p <- obj$Prob_ATE_smaller_equal_0[[2]] + ggtitle(paste("Population", i, "Censoring", j, "%"))
    
    ggsave(filename = paste0("Output/Plots/Prob_smaller_0/Prob_ATE_smaller_equal_0_",
                             "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

    
  }
} 
# =============================================================================================




# Prob ATE smaller equal true ATE =============================================================
# Loop through the scenarios 1 to 4
for (i in 1:4) {
  # Loop through the censoring levels: 20, 40, 60
  for (j in c(20, 40, 60)) {
    
    # Construct the object name as a string
    obj_name <- paste0("all_plots_", i, "_", j)
    
    # Get object
    obj <- get(obj_name)

    p <- obj$Prob_ATE_smaller_equal_true_treatment_effect[[1]] + plot_annotation(title = paste("Population", i, "Censoring", j, "%"))

    ggsave(filename = paste0("Output/Plots/Prob_smaller_ATE/Prob_ATE_smaller_equal_true_ATE_grid_",
                             "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 5) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
    
    p <- obj$Prob_ATE_smaller_equal_true_treatment_effect[[2]] + ggtitle(paste("Population", i, "Censoring", j, "%"))

    ggsave(filename = paste0("Output/Plots/Prob_smaller_ATE/Prob_ATE_smaller_equal_true_ATE_",
                             "Scenario_", i, "_censoring_", j, ".svg"),
           plot = p,
           units = "in",
           width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
    
  }
} 
# =============================================================================================





# ========================================================================================================================================
##### Visualize Censoring Distribution #####
# ========================================================================================================================================

library(tidyverse)

source("Script/Functions/Generate_Data.R")

set.seed(12334)

n <- 8000

censoring_distribution_data <- rbind(Generate_Data(n = n, max_censoring_time_unif = 20, censoring_rate = 0.02),
                                     Generate_Data(n = n, max_censoring_time_unif = 20, censoring_rate = 0.06),
                                     Generate_Data(n = n, max_censoring_time_unif = 15, censoring_rate = 0.08))

censoring_distribution_data$censoring_group <- rep(c("20%", "40%", "60%"), each = n)

censoring_distribution_data <- censoring_distribution_data %>%
  group_by(censoring_group) %>%
  arrange(censoring_time, .by_group = TRUE) %>%
  mutate(cdf = seq_along(censoring_time) / length(censoring_time))

p <- ggplot(censoring_distribution_data, aes(x = censoring_time, y = cdf, linetype = censoring_group)) +
  geom_step(lwd = 1.2) +
  labs(
    title = "Cumulative Distribution of Censoring Times",
    x = "Time",
    y = "Proportion Censored"
  ) +
  theme_classic() +
  geom_vline(xintercept = 10, linetype = "dotted") +
  scale_linetype_manual(values = c("dotted", "dotdash", "dashed")) +
  guides(linetype = "none")# +
  #scale_linetype_manual(name = "Amount of Censoring", # change legend label
                        #values = c("dotted", "dashed", "solid")) 

ggsave(filename = paste0("Output/Plots/Censoring Visualization/CDF_censoring.svg"),
       plot = p,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6



# ========================================================================================================================================
##### Visualize Variable Distributions #####
# ========================================================================================================================================

library(tidyverse)

source("Script/Functions/Generate_Data.R")

set.seed(12334)

n <- 30000

df <- Generate_Data(n = n) 
df$A <- factor(df$Treatment)
df <- df %>% select(A, X1, X2, X3, X4)

df_long <- gather(df, key = variable, value = Value, X1:X4, factor_key = TRUE)

p <- ggplot(data = df_long, aes(x = A, y = Value, color = A)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_classic() +
  scale_color_manual(values=c("red", "blue"), name = "", labels = c("A = 0", "A = 1")) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ variable) +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(0.5,"cm")
  ) #+
  #scale_color_discrete(name = "", labels = c("A = 0", "A = 1"))

ggsave(filename = paste0("Output/Plots/Variable Distributions/Variable Distributions.svg"),
       plot = p,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6










# ========================================================================================================================================
##### Show confidence bands for RCS smoothing with RF and GEE_splines #####
# ========================================================================================================================================

# load libaries
library(tidyverse)
library(rms)

# load needed functions
source("Script/Functions/plot_continuous.R")

# load data and prepare it
load("Output/Scenario_1_censoring_40.RData")
results_df <- results_df %>%
  filter(model %in% c("RF", "GEE_splines")) %>%
  mutate(error_time_root_n = abs(bias) * sqrt(n))

model_colors = c(
  "GEE" = "#1b9e77",       # Teal
  "GEE_splines" = "#d95f02", # Orange
  "PO" = "#7570b3",        # Purple
  "PO_splines" = "#e7298a", # Pink
  "WB" = "#e6ab02",        # Mustard
  "RF" = "#66a61e",        # Green
  "CRF" = "#A0522D"        # Brown
)


# RMSE
p1 <- ggplot(results_df, aes(x = n, y = error_time_root_n)) + 
  stat_smooth(data = results_df,
              aes(color = model),
              method = "lm",
              se = TRUE,
              formula = y ~ rcs(x, c(125, 250, 500, 1000)),
              size = 1) +
  labs(
    title = "Population 1 Censoring 40 %",
    y = "SMAD",
    x = "Sample Size"
  ) +
  scale_color_manual(
    name = "Model",  # change
    values = model_colors,
    labels = c(
      "PO" = "PO",
      "PO_splines" = expression(PO[complex]), # change
      "GEE" = "GEE",
      "GEE_splines" = expression(GEE[complex]), # change
      "RF" = "RF",
      "CRF" = "CRF",
      "WB" = "WB"
    )
  ) +
  scale_fill_manual(values = model_colors) +
  theme_minimal()


# treatment effect
p2 <- ggplot(results_df, aes(x = log2(n), y = treatment_effect)) +
  stat_smooth(aes(color = model),
    method = "lm",
    se = TRUE,
    formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
    size = 1
  ) +
  labs(
    title = "Population 1 Censoring 40 %",
    y = "MATE",
    x = "Sample Size"
  ) +
  scale_x_continuous(
    breaks = log2(c(125, 250, 500, 1000)),
    labels = c("125", "250", "500", "1000")
  ) +
  scale_color_manual(
    name = "Model", 
    values = model_colors,
    labels = c(
      "PO" = "PO",
      "PO_splines" = expression(PO[complex]), 
      "GEE" = "GEE",
      "GEE_splines" = expression(GEE[complex]), 
      "RF" = "RF",
      "CRF" = "CRF",
      "WB" = "WB"
    )
  ) +
  theme_minimal() +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect)) +
  scale_y_continuous(sec.axis = sec_axis(transform = ~ . - first(results_df$true_treatment_effect),
                                         name = "Bias")
  ) +
  coord_cartesian(ylim = c(results_df$true_treatment_effect[1] - 0.45, results_df$true_treatment_effect[1] + 0.45))


# CI width
p3 <- ggplot(results_df, aes(x = log2(n), y = CI_width)) +
  stat_smooth(aes(color = model),
              method = "lm",
              se = TRUE,
              formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
              size = 1
  ) +
  labs(
    title = "Population 1 Censoring 40 %",
    y = "MCW",
    x = "Sample Size"
  ) +
  scale_x_continuous(
    breaks = log2(c(125, 250, 500, 1000)),
    labels = c("125", "250", "500", "1000")
  ) +
  scale_color_manual(
    name = "Model", 
    values = model_colors,
    labels = c(
      "PO" = "PO",
      "PO_splines" = expression(PO[complex]), 
      "GEE" = "GEE",
      "GEE_splines" = expression(GEE[complex]), 
      "RF" = "RF",
      "CRF" = "CRF",
      "WB" = "WB"
    )
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 3.2))


# Prob ATE smaller equal 0
p4 <- ggplot(results_df, aes(x = log2(n), y = Prob_ATE_smaller_equal_0)) +
  stat_smooth(aes(color = model),
              method = "lm",
              se = TRUE,
              formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
              size = 1
  ) +
  labs(
    title = "Population 1 Censoring 40 %",
    y = "Mean p-value",
    x = "Sample Size"
  ) +
  scale_x_continuous(
    breaks = log2(c(125, 250, 500, 1000)),
    labels = c("125", "250", "500", "1000")
  ) +
  scale_color_manual(
    name = "Model", 
    values = model_colors,
    labels = c(
      "PO" = "PO",
      "PO_splines" = expression(PO[complex]), 
      "GEE" = "GEE",
      "GEE_splines" = expression(GEE[complex]), 
      "RF" = "RF",
      "CRF" = "CRF",
      "WB" = "WB"
    )
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.5))

# Prob ATE smaller equal true value
p5 <- ggplot(results_df, aes(x = log2(n), y = Prob_ATE_smaller_equal_true_treatment_effect)) +
  stat_smooth(aes(color = model),
              method = "lm",
              se = TRUE,
              formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
              size = 1
  ) +
  labs(
    title = "Population 1 Censoring 40 %",
    y = expression("Mean " * hat(P)[B](ATE)),
    x = "Sample Size"
  ) +
  scale_x_continuous(
    breaks = log2(c(125, 250, 500, 1000)),
    labels = c("125", "250", "500", "1000")
  ) +
  scale_color_manual(
    name = "Model", 
    values = model_colors,
    labels = c(
      "PO" = "PO",
      "PO_splines" = expression(PO[complex]), 
      "GEE" = "GEE",
      "GEE_splines" = expression(GEE[complex]), 
      "RF" = "RF",
      "CRF" = "CRF",
      "WB" = "WB"
    )
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1))

# save plots as svg

ggsave(filename = "Output/Plots/CI showcase/SMAD.svg",
       plot = p1,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/MATE.svg",
       plot = p2,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/MCW.svg",
       plot = p3,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/Mean_p_value.svg",
       plot = p4,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/Mean_P.svg",
       plot = p5,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

