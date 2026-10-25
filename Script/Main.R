# ========================================================================================================================================
##### Main #####
# ========================================================================================================================================


# ========================================================================================================================================
##### Simulation #####
# ========================================================================================================================================


rm(list = ls())

source("Script/Functions/Simulate.R")

n <- 1000

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
# scenarios["2"] <- NULL  # Uncomment this line to skip scenario 2
# scenarios["3"] <- NULL  # Uncomment this line to skip scenario 3
# scenarios["4"] <- NULL  # Uncomment this line to skip scenario 4

# Define censoring settings
censoring_settings <- list(
  list(max_censoring_time_unif = 20, censoring_rate = 0.02, label = 20),
  list(max_censoring_time_unif = 20, censoring_rate = 0.06, label = 40),
  list(max_censoring_time_unif = 15, censoring_rate = 0.08, label = 60)
)

start <- Sys.time()

set.seed(327)

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
                       "_censoring_", censoring$label, ".RData")
    save(results_df, file = filename)
  }
}

end <- Sys.time()
end-start











# ========================================================================================================================================
##### Visualization of Scenarios with Survival Curves #####
# ========================================================================================================================================

if (!require("survival")) {install.packages("survival", dependencies=TRUE); library(survival)}
source("Script/Functions/Generate_Data.R")

## Simulate the true and the descriptive survival curves for one population:
simulate_fits <- function(beta_treatment, treatment_shape_effect) {
  
  ## Simulate everyone treated:
  set.seed(1)
  data_all_treatment <- Generate_Data(
    n = 20000,
    p_treat_manual = 1,
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring = FALSE
  )
  
  ## Simulate no one treated:
  set.seed(1)
  data_no_treatment <- Generate_Data(
    n = 20000,
    p_treat_manual = 0,
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring = FALSE
  )
  
  ## Fit the true curves on the combined data:
  data <- rbind(data_all_treatment, data_no_treatment)
  fit_true <- survfit(Surv(ObservedTime, Status) ~ Treatment, data = data)
  
  ## Fit the descriptive curves on normally simulated data:
  data_simulated_normally <- Generate_Data(
    n = 20000,
    beta_treatment = beta_treatment,
    treatment_shape_effect = treatment_shape_effect,
    censoring = FALSE
  )
  fit_descriptive <- survfit(Surv(ObservedTime, Status) ~ Treatment, data = data_simulated_normally)
  
  list(true = fit_true, descriptive = fit_descriptive)
}

## Draw one Kaplan Meier panel:
plot_panel <- function(fits, Scenario) {
  
  ## Plot the true curves:
  plot(fits$true, col = c("red", "blue"), lty = c(1, 1), lwd = 1.5, xlab = "Time",
       xlim = c(0, 11), ylim = c(0.3, 1), main = paste("Population", Scenario))
  
  ## Add the descriptive curves:
  lines(fits$descriptive, col = c("red", "blue"), lty = c(2, 2))
  
  ## Add the vertical line at time 10:
  segments(x0 = 10, y0 = 0.3, x1 = 10, y1 = 1, lty = 2, col = "black")
}

## Define parameter combinations:
params <- data.frame(
  beta_treatment = c(-0.5, 0, 0.5, -1),
  treatment_shape_effect = c(0, 0, -0.3, 0.49),
  Scenario = c(1, 2, 3, 4)
)

## Open the graphics device and set up a 2x2 layout with space for a shared legend at the bottom:
png("Output/Plots/Scenario Visualization/Kaplan_Meier_all_populations.png",
    width = 7, height = 6.5, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(3, 0, 0, 0))

## Create one panel per population:
for (i in 1:nrow(params)) {
  fits <- simulate_fits(
    beta_treatment = params$beta_treatment[i],
    treatment_shape_effect = params$treatment_shape_effect[i]
  )
  plot_panel(fits, params$Scenario[i])
}

## Add one shared legend in the bottom margin:
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom", ncol = 2, bty = "n", xpd = TRUE,
       legend = c("Control (True)", "Treatment (True)",
                  "Control (Descriptive)", "Treatment (Descriptive)"),
       col = c("red", "blue", "red", "blue"),
       lty = c(1, 1, 2, 2))

dev.off()










# # ========================================================================================================================================
# ##### Calculate Raw ATE for each population #####
# # ========================================================================================================================================
# 
# source("Script/Functions/Generate_Data.R")
# if (!require("survival")) {install.packages("survival", dependencies=TRUE); library(survival)}
# 
# scenarios <- list(
#   "1" = list(beta_treatment = -0.5, treatment_shape_effect = 0),
#   "2" = list(beta_treatment = 0, treatment_shape_effect = 0),
#   "3" = list(beta_treatment = 0.5, treatment_shape_effect = -0.3),
#   "4" = list(beta_treatment = -1, treatment_shape_effect = 0.49)
# )
# tau = 10
# 
# for(i in 1:4) {
#   set.seed(1)
#   scenario <- scenarios[[i]]
#   data <- Generate_Data(n = 100000,
#                         beta_treatment = scenario$beta_treatment,
#                         treatment_shape_effect = scenario$treatment_shape_effect,
#                         censoring = FALSE
#   ) 
#   
#   
#   # Compute the Kaplan-Meier curves
#   km_fit<- survfit(Surv(ObservedTime, Status) ~ Treatment, data = data)  # One curve for all individuals
#   
#   # Extract RMST at time tau
#   RMST <- survival:::survmean(km_fit, rmean=tau)$matrix[1:2, "rmean"]
#   print(paste0("Population ", i, ": ", RMST["Treatment=1"] - RMST["Treatment=0"]))
# }










# ========================================================================================================================================
##### Visualization of Pseudo-Observations #####
# ========================================================================================================================================

if (!require("survival")) {install.packages("survival", dependencies=TRUE); library(survival)}
if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}
if (!require("patchwork")) {install.packages("patchwork", dependencies=TRUE); library(patchwork)}
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
ggsave(filename = "Output/Plots/PO examples/PO_RMST_low_censoring.png",
       plot = p1,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6
ggsave(filename = "Output/Plots/PO examples/PO_RMST_high_censoring.png",
       plot = p3,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6













# ========================================================================================================================================
##### Visualization of Performance Measures #####
# ========================================================================================================================================



# Create combined 2x2 plots (one figure per performance measure, 4 populations)

## Load packages:
if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}
if (!require("patchwork")) {install.packages("patchwork", dependencies=TRUE); library(patchwork)}
if (!require("rms")) {install.packages("rms", dependencies=TRUE); library(rms)}

## Load helper functions:
source("Script/Functions/plot_continuous.R")
source("Script/Functions/plot_logistic.R")
source("Script/Functions/create_plots.R")

## Choose the censoring levels to create figures for:
censoring_levels <- c(20, 40, 60)

## Define the y-axis limits of the Mean p-value plot per population (copied from the Main file):
p0_limits <- list(
  "1" = c(0, 0.5),
  "2" = c(0.2, 0.8),
  "3" = c(0.15, 0.75),
  "4" = c(0, 0.7)
)

## Define the y-axis range around the true ATE in the MATE plot (population 1 at 60% censoring is wider):
ate_added_limit <- function(population, censoring) {
  if (population == 1 && censoring == 60) return(0.5)
  if (population == 1) return(0.45)
  return(0.35)
}

## Define the four measures (name in create_plots output, index of the all-methods plot, file name):
measures <- tibble::tribble(
  ~plot_name,                  ~plot_index, ~file,
  "RMSE_plot",                 1,           "SMAD",
  "treatment_effect",          2,           "MATE",
  "CI_width",                  2,           "MCW",
  "Prob_ATE_smaller_equal_0",  2,           "Mean_p_value"
)

## Function that combines the 4 populations of one measure into a 2x2 figure:
combine_populations <- function(all_plots, plot_name, plot_index) {
  
  ## Take the all-methods-in-one plot of each population and add the population title:
  panels <- lapply(1:4, function(i) {
    all_plots[[i]][[plot_name]][[plot_index]] +
      ggtitle(paste("Population", i)) +
      theme(plot.title = element_text(size = 10, face = "bold"))
  })
  
  ## Arrange as 2x2 with one shared legend at the bottom and tags A-D:
  wrap_plots(panels, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
}

## Create the figures for each censoring level:
## Note: create_plots() reads results_df from the global environment, so load() must run at top level:
dir.create("Output/Plots/Combined", recursive = TRUE, showWarnings = FALSE)

for (censoring in censoring_levels) {
  
  ## Create all individual plots for each population:
  all_plots <- list()
  
  for (i in 1:4) {
    
    ## Load the results of population i (this creates the object results_df):
    load(paste0("Output/Scenario_", i, "_censoring_", censoring, ".RData"))
    
    ## Create the plots of population i:
    all_plots[[i]] <- create_plots(
      data = results_df,
      treatment_effect_added_y_limits = ate_added_limit(i, censoring),
      CI_width_y_limits = c(0, 3.2),
      Prob_ATE_smaller_equal_0_y_limits = p0_limits[[as.character(i)]],
      Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1)
    )
  }
  
  ## Create and save one 2x2 figure per measure:
  for (k in seq_len(nrow(measures))) {
    
    p <- combine_populations(all_plots, measures$plot_name[k], measures$plot_index[k])
    
    ggsave(filename = paste0("Output/Plots/Combined/", measures$file[k], "_censoring_", censoring, ".png"),
           plot = p,
           units = "in",
           width = 7, height = 6)
  }
}



# ========================================================================================================================================
##### Visualize Censoring Distribution #####
# ========================================================================================================================================

if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}

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

ggsave(filename = paste0("Output/Plots/Censoring Visualization/CDF_censoring.png"),
       plot = p,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6



# ========================================================================================================================================
##### Visualize Variable Distributions #####
# ========================================================================================================================================

if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}

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

ggsave(filename = paste0("Output/Plots/Variable Distributions/Variable Distributions.png"),
       plot = p,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6










# ========================================================================================================================================
##### Show confidence bands for RCS smoothing with RF and GEE_splines #####
# ========================================================================================================================================

# load libaries
if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}
if (!require("rms")) {install.packages("rms", dependencies=TRUE); library(rms)}

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

ggsave(filename = "Output/Plots/CI showcase/SMAD.png",
       plot = p1,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/MATE.png",
       plot = p2,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/MCW.png",
       plot = p3,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/Mean_p_value.png",
       plot = p4,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6

ggsave(filename = "Output/Plots/CI showcase/Mean_P.png",
       plot = p5,
       units = "in",
       width = 4.7, height = 3) # 4.7 is exactly the width of the latex document, 3 is ~4.7/1.6





# ========================================================================================================================================
##### Literature Review #####
# ========================================================================================================================================

if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}
if (!require("magrittr")) {install.packages("magrittr", dependencies=TRUE); library(magrittr)}
if (!require("readxl")) {install.packages("readxl", dependencies=TRUE); library(readxl)}
literature_review <- readxl::read_excel("Data/Papers that cited Andersen RMST 2004 paper.xlsx")

literature_review <- literature_review %>%
  mutate(row_id = row_number()) %>%
  rename(used_PO = `Used PO in Analysis`,
         used_GLM_on_PO = `Used GLM framework to analyse PO`,
         used_other_on_PO = `Used another method to analyse PO`) %T>%
  {print(nrow(.)) } %>%
  filter(!is.na(used_PO)) %T>% # there were some citations that were not publications, these are them
  {print(nrow(.)) } %>%
  filter(used_GLM_on_PO != 'unclear') %T>% # in one paper it wasn't clear, this is taken out
  {print(nrow(.)) } %>%
  mutate(used_GLM_on_PO = as.numeric(used_GLM_on_PO),
         used_other_on_PO = as.numeric(used_other_on_PO))

# save(literature_review, file = "Output/literature_review_cleaned.RData")
load("Output/literature_review_cleaned.RData")

ggplot(data = literature_review, aes(x = year)) +
  geom_bar()

df <- literature_review %>%
  group_by(year) %>%
  summarise(sum_used_PO = sum(used_PO),
            sum_used_GLM_on_PO = sum(used_GLM_on_PO),
            sum_used_other_on_PO = sum(used_other_on_PO))

ggplot(data = df, aes(y = sum_used_other_on_PO, x = year)) +
  geom_point()

literature_review %>%
  filter(is.na(used_other_on_PO))

# share of other methods in PO papers
literature_review %>%
  filter(used_PO == 1 & year %in% c(2022:2026)) %>%
  summarise(number_other = sum(used_other_on_PO),
            number_GLM = sum(used_GLM_on_PO),
            number_PO = n())

literature_review %>%
  filter(used_PO == 1 & year %in% c(2007:2013)) %>%
  summarise(number_other = sum(used_other_on_PO),
            number_GLM = sum(used_GLM_on_PO),
            number_PO = n())

literature_review %>%
  filter(used_PO == 1 & year %in% c(2014:2017)) %>%
  summarise(number_other = sum(used_other_on_PO),
            number_GLM = sum(used_GLM_on_PO),
            number_PO = n())

literature_review %>%
  filter(used_PO == 1 & year %in% c(2018:2021)) %>%
  summarise(number_other = sum(used_other_on_PO),
            number_GLM = sum(used_GLM_on_PO),
            number_PO = n())

  
