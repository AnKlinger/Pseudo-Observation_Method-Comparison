# ========================================================================================================================================
# create_plots Function
# ========================================================================================================================================

create_plots <- function(data,
                         treatment_effect_added_y_limits = 0.7,
                         CI_width_y_limits = c(0, 3.2),
                         Prob_ATE_smaller_equal_0_y_limits = c(0, 0.5),
                         Prob_ATE_smaller_equal_true_treatment_effect_y_limits = c(0, 1),
                         model_colors = c(
                           "GEE" = "#1b9e77",       # Teal
                           "GEE_splines" = "#d95f02", # Orange
                           "PO" = "#7570b3",        # Purple
                           "PO_splines" = "#e7298a", # Pink
                           "WB" = "#e6ab02",        # Mustard
                           "RF" = "#66a61e",        # Green
                           "CRF" = "#A0522D"        # Brown
                         )) {

if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}
if (!require("patchwork")) {install.packages("patchwork", dependencies=TRUE); library(patchwork)}
  
# Add factor levels to results_df for legend
results_df$model <- factor(results_df$model,
                           levels = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"))

# Create empty list to save results
all_plots <- list()
  
# Create two plots for treatment_effect

p1 <- plot_continuous(results_df %>% filter(model %in% c("PO", "WB")),
                      "treatment_effect",
                      "MATE",
                      highlight_models = c("PO"),
                      show_points = FALSE,
                      y_limits = c(results_df$true_treatment_effect[1] - treatment_effect_added_y_limits,
                                   results_df$true_treatment_effect[1] + treatment_effect_added_y_limits),
                      boxplots = TRUE) +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect),
              linetype = "longdash")

p2 <- plot_continuous(results_df %>% filter(model %in% c("PO_splines", "WB")),
                      "treatment_effect",
                      "MATE",
                      highlight_models = c("PO_splines"),
                      show_points = FALSE,
                      y_limits = c(results_df$true_treatment_effect[1] - treatment_effect_added_y_limits,
                                   results_df$true_treatment_effect[1] + treatment_effect_added_y_limits),
                      boxplots = TRUE) +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect),
              linetype = "longdash")

p3 <- plot_continuous(results_df %>% filter(model %in% c("GEE", "WB")),
                      "treatment_effect",
                      "MATE",
                      highlight_models = c("GEE"),
                      show_points = FALSE,
                      y_limits = c(results_df$true_treatment_effect[1] - treatment_effect_added_y_limits,
                                   results_df$true_treatment_effect[1] + treatment_effect_added_y_limits),
                      boxplots = TRUE) +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect),
              linetype = "longdash")

p4 <- plot_continuous(results_df %>% filter(model %in% c("GEE_splines", "WB")),
                      "treatment_effect",
                      "MATE",
                      highlight_models = c("GEE_splines"),
                      show_points = FALSE,
                      y_limits = c(results_df$true_treatment_effect[1] - treatment_effect_added_y_limits,
                                   results_df$true_treatment_effect[1] + treatment_effect_added_y_limits),
                      boxplots = TRUE) +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect),
              linetype = "longdash")

p5 <- plot_continuous(results_df %>% filter(model %in% c("RF", "WB")),
                      "treatment_effect",
                      "MATE",
                      highlight_models = c("RF"),
                      show_points = FALSE,
                      y_limits = c(results_df$true_treatment_effect[1] - treatment_effect_added_y_limits,
                                   results_df$true_treatment_effect[1] + treatment_effect_added_y_limits),
                      boxplots = TRUE) +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect),
              linetype = "longdash")

p6 <- plot_continuous(results_df %>% filter(model %in% c("CRF", "WB")),
                      "treatment_effect",
                      "MATE",
                      highlight_models = c("CRF"),
                      show_points = FALSE,
                      y_limits = c(results_df$true_treatment_effect[1] - treatment_effect_added_y_limits,
                                   results_df$true_treatment_effect[1] + treatment_effect_added_y_limits),
                      boxplots = TRUE) +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect),
              linetype = "longdash")

# Remove legends and title
p1 <- p1 + theme(legend.position = "none") + ggtitle("")
p2 <- p2 + theme(legend.position = "none") + ggtitle("")
p3 <- p3 + theme(legend.position = "none") + ggtitle("")
p4 <- p4 + theme(legend.position = "none") + ggtitle("")
p5 <- p5 + theme(legend.position = "none") + ggtitle("")
p6 <- p6 + theme(legend.position = "none") + ggtitle("")

# Create a legend by hand
legend_data <- data.frame(
  method = factor(c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"),
                  levels = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB")),
  x = 1:7,
  y = 1
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = method)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = model_colors,
                     labels = c(
                       "PO" = "PO",
                       "PO_splines" = expression(PO[complex]), # change
                       "GEE" = "GEE",
                       "GEE_splines" = expression(GEE[complex]), # change
                       "RF" = "RF",
                       "CRF" = "CRF",
                       "WB" = "WB"
                     )) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

Six_plots_grid <- (p1 + p2 + p3 + p4 + p5 + p6) + 
  plot_layout(ncol = 2, axes = "collect")

Six_plots_grid <- Six_plots_grid / 
  legend_plot + plot_layout(heights = c(10, 0)) # the number 10 is irrelevant, it's the same as long as no 0 is there

Six_methods_in_one <- plot_continuous(results_df,
                                      "treatment_effect",
                                      "MATE",
                                      highlight_models = c(),
                                      show_points = FALSE,
                                      y_limits = c(results_df$true_treatment_effect[1] - treatment_effect_added_y_limits,
                                                   results_df$true_treatment_effect[1] + treatment_effect_added_y_limits),
                                      boxplots = TRUE) +
  geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect),
              linetype = "longdash") +
  scale_y_continuous(sec.axis = sec_axis(transform = ~ . - first(results_df$true_treatment_effect),
                                         name = "Bias")
  )

# rename legend stuff

Six_methods_in_one <- Six_methods_in_one +
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
  )

treatment_effect <- list(Six_plots_grid, Six_methods_in_one)
all_plots[["treatment_effect"]] <- treatment_effect

# Create two plots for for CI_excludes_0

p1 <- plot_logistic(results_df %>% filter(model %in% c("PO", "WB")),
                    "CI_excludes_0",
                    "CI Excludes 0",
                    highlight_models = c("PO"),
                    show_points = FALSE)

p2 <- plot_logistic(results_df %>% filter(model %in% c("PO_splines", "WB")),
                    "CI_excludes_0",
                    "CI Excludes 0",
                    highlight_models = c("PO_splines"),
                    show_points = FALSE)

p3 <- plot_logistic(results_df %>% filter(model %in% c("GEE", "WB")),
                    "CI_excludes_0",
                    "CI Excludes 0",
                    highlight_models = c("GEE"),
                    show_points = FALSE)

p4 <- plot_logistic(results_df %>% filter(model %in% c("GEE_splines", "WB")),
                    "CI_excludes_0",
                    "CI Excludes 0",
                    highlight_models = c("GEE_splines"),
                    show_points = FALSE)

p5 <- plot_logistic(results_df %>% filter(model %in% c("RF", "WB")),
                    "CI_excludes_0",
                    "CI Excludes 0",
                    highlight_models = c("RF"),
                    show_points = FALSE)

p6 <- plot_logistic(results_df %>% filter(model %in% c("CRF" ,"WB")),
                    "CI_excludes_0",
                    "CI Excludes 0",
                    highlight_models = c("CRF"),
                    show_points = FALSE)

# Remove legends and title
p1 <- p1 + theme(legend.position = "none") + ggtitle("")
p2 <- p2 + theme(legend.position = "none") + ggtitle("")
p3 <- p3 + theme(legend.position = "none") + ggtitle("")
p4 <- p4 + theme(legend.position = "none") + ggtitle("")
p5 <- p5 + theme(legend.position = "none") + ggtitle("")
p6 <- p6 + theme(legend.position = "none") + ggtitle("")

# Create a legend by hand
legend_data <- data.frame(
  method = factor(c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"),
                  levels = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB")),
  x = 1:7,
  y = 1
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = method)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = model_colors,
                     labels = c(
                       "PO" = "PO",
                       "PO_splines" = expression(PO[complex]), # change
                       "GEE" = "GEE",
                       "GEE_splines" = expression(GEE[complex]), # change
                       "RF" = "RF",
                       "CRF" = "CRF",
                       "WB" = "WB"
                     )) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

Six_plots_grid <- (p1 + p2 + p3 + p4 + p5 + p6) + 
  plot_layout(ncol = 2, axes = "collect")

Six_plots_grid <- Six_plots_grid / 
  legend_plot + plot_layout(heights = c(10, 0)) # the number 10 is irrelevant, it's the same as long as no 0 is there

Six_methods_in_one <- plot_logistic(results_df,
                                    "CI_excludes_0",
                                    "CI Excludes 0",
                                    highlight_models = c(),
                                    show_points = FALSE)

CI_excludes_0 <- list(Six_plots_grid, Six_methods_in_one)
all_plots[["CI_excludes_0"]] <- CI_excludes_0

# Create two plots for for CI_includes_true_value

p1 <- plot_logistic(results_df %>% filter(model %in% c("PO", "WB")),
                    "CI_includes_true_value",
                    "CI Includes True Value",
                    highlight_models = c("PO"),
                    show_points = FALSE)

p2 <- plot_logistic(results_df %>% filter(model %in% c("PO_splines", "WB")),
                    "CI_includes_true_value",
                    "CI Includes True Value",
                    highlight_models = c("PO_splines"),
                    show_points = FALSE)

p3 <- plot_logistic(results_df %>% filter(model %in% c("GEE", "WB")),
                    "CI_includes_true_value",
                    "CI Includes True Value",
                    highlight_models = c("GEE"),
                    show_points = FALSE)

p4 <- plot_logistic(results_df %>% filter(model %in% c("GEE_splines", "WB")),
                    "CI_includes_true_value",
                    "CI Includes True Value",
                    highlight_models = c("GEE_splines"),
                    show_points = FALSE)

p5 <- plot_logistic(results_df %>% filter(model %in% c("RF", "WB")),
                    "CI_includes_true_value",
                    "CI Includes True Value",
                    highlight_models = c("RF"),
                    show_points = FALSE)

p6 <- plot_logistic(results_df %>% filter(model %in% c("CRF" ,"WB")),
                    "CI_includes_true_value",
                    "CI Includes True Value",
                    highlight_models = c("CRF"),
                    show_points = FALSE)

# Remove legends and title
p1 <- p1 + theme(legend.position = "none") + ggtitle("")
p2 <- p2 + theme(legend.position = "none") + ggtitle("")
p3 <- p3 + theme(legend.position = "none") + ggtitle("")
p4 <- p4 + theme(legend.position = "none") + ggtitle("")
p5 <- p5 + theme(legend.position = "none") + ggtitle("")
p6 <- p6 + theme(legend.position = "none") + ggtitle("")

# Create a legend by hand
legend_data <- data.frame(
  method = factor(c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"),
                  levels = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB")),
  x = 1:7,
  y = 1
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = method)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = model_colors,
                     labels = c(
                       "PO" = "PO",
                       "PO_splines" = expression(PO[complex]), # change
                       "GEE" = "GEE",
                       "GEE_splines" = expression(GEE[complex]), # change
                       "RF" = "RF",
                       "CRF" = "CRF",
                       "WB" = "WB"
                     )) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

Six_plots_grid <- (p1 + p2 + p3 + p4 + p5 + p6) + 
  plot_layout(ncol = 2, axes = "collect")

Six_plots_grid <- Six_plots_grid / 
  legend_plot + plot_layout(heights = c(10, 0)) # the number 10 is irrelevant, it's the same as long as no 0 is there

Six_methods_in_one <- plot_logistic(results_df,
                                    "CI_includes_true_value",
                                    "CI Includes True Value",
                                    highlight_models = c(),
                                    show_points = FALSE)

# rename legend stuff

Six_methods_in_one <- Six_methods_in_one +
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
  )

CI_includes_true_value <- list(Six_plots_grid, Six_methods_in_one)
all_plots[["CI_includes_true_value"]] <- CI_includes_true_value

# Create two plots for CI_width

p1 <- plot_continuous(results_df %>% filter(model %in% c("PO", "WB")),
                      "CI_width",
                      "MCW",
                      highlight_models = c("PO"),
                      show_points = FALSE,
                      y_limits = CI_width_y_limits,
                      boxplots = TRUE)

p2 <- plot_continuous(results_df %>% filter(model %in% c("PO_splines", "WB")),
                      "CI_width",
                      "MCW",
                      highlight_models = c("PO_splines"),
                      show_points = FALSE,
                      y_limits = CI_width_y_limits,
                      boxplots = TRUE)

p3 <- plot_continuous(results_df %>% filter(model %in% c("GEE", "WB")),
                      "CI_width",
                      "MCW",
                      highlight_models = c("GEE"),
                      show_points = FALSE,
                      y_limits = CI_width_y_limits,
                      boxplots = TRUE)

p4 <- plot_continuous(results_df %>% filter(model %in% c("GEE_splines", "WB")),
                      "CI_width",
                      "MCW",
                      highlight_models = c("GEE_splines"),
                      show_points = FALSE,
                      y_limits = CI_width_y_limits,
                      boxplots = TRUE)

p5 <- plot_continuous(results_df %>% filter(model %in% c("RF", "WB")),
                      "CI_width",
                      "MCW",
                      highlight_models = c("RF"),
                      show_points = FALSE,
                      y_limits = CI_width_y_limits,
                      boxplots = TRUE)

p6 <- plot_continuous(results_df %>% filter(model %in% c("CRF" ,"WB")),
                      "CI_width",
                      "MCW",
                      highlight_models = c("CRF"),
                      show_points = FALSE,
                      y_limits = CI_width_y_limits,
                      boxplots = TRUE)

# Remove legends and title
p1 <- p1 + theme(legend.position = "none") + ggtitle("")
p2 <- p2 + theme(legend.position = "none") + ggtitle("")
p3 <- p3 + theme(legend.position = "none") + ggtitle("")
p4 <- p4 + theme(legend.position = "none") + ggtitle("")
p5 <- p5 + theme(legend.position = "none") + ggtitle("")
p6 <- p6 + theme(legend.position = "none") + ggtitle("")

# Create a legend by hand
legend_data <- data.frame(
  method = factor(c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"),
                  levels = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB")),
  x = 1:7,
  y = 1
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = method)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = model_colors,
                     labels = c(
                       "PO" = "PO",
                       "PO_splines" = expression(PO[complex]), # change
                       "GEE" = "GEE",
                       "GEE_splines" = expression(GEE[complex]), # change
                       "RF" = "RF",
                       "CRF" = "CRF",
                       "WB" = "WB"
                     )) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

Six_plots_grid <- (p1 + p2 + p3 + p4 + p5 + p6) + 
  plot_layout(ncol = 2, axes = "collect")

Six_plots_grid <- Six_plots_grid / 
  legend_plot + plot_layout(heights = c(10, 0)) # the number 10 is irrelevant, it's the same as long as no 0 is there

Six_methods_in_one <- plot_continuous(results_df,
                                      "CI_width",
                                      "MCW",
                                      highlight_models = c(),
                                      show_points = FALSE,
                                      y_limits = CI_width_y_limits,
                                      boxplots = TRUE)

# rename legend stuff

Six_methods_in_one <- Six_methods_in_one +
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
  )

CI_width <- list(Six_plots_grid, Six_methods_in_one)
all_plots[["CI_width"]] <- CI_width

# Create two plots for Prob_ATE_smaller_equal_0

p1 <- plot_continuous(results_df %>% filter(model %in% c("PO", "WB")),
                      "Prob_ATE_smaller_equal_0",
                      "Mean p-value",
                      highlight_models = c("PO"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_0_y_limits,
                      boxplots = TRUE)

#p1 <- p1 + labs(y = expression(hat(P)[B](0)))

p2 <- plot_continuous(results_df %>% filter(model %in% c("PO_splines", "WB")),
                      "Prob_ATE_smaller_equal_0",
                      "Mean p-value",
                      highlight_models = c("PO_splines"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_0_y_limits,
                      boxplots = TRUE)

#p2 <- p2 + labs(y = expression(hat(P)[B](0)))

p3 <- plot_continuous(results_df %>% filter(model %in% c("GEE", "WB")),
                      "Prob_ATE_smaller_equal_0",
                      "Mean p-value",
                      highlight_models = c("GEE"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_0_y_limits,
                      boxplots = TRUE)

#p3 <- p3 + labs(y = expression(hat(P)[B](0)))

p4 <- plot_continuous(results_df %>% filter(model %in% c("GEE_splines", "WB")),
                      "Prob_ATE_smaller_equal_0",
                      "Mean p-value",
                      highlight_models = c("GEE_splines"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_0_y_limits,
                      boxplots = TRUE)

#p4 <- p4 + labs(y = expression(hat(P)[B](0)))

p5 <- plot_continuous(results_df %>% filter(model %in% c("RF", "WB")),
                      "Prob_ATE_smaller_equal_0",
                      "Mean p-value",
                      highlight_models = c("RF"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_0_y_limits,
                      boxplots = TRUE)

#p5 <- p5 + labs(y = expression(hat(P)[B](0)))

p6 <- plot_continuous(results_df %>% filter(model %in% c("CRF" ,"WB")),
                      "Prob_ATE_smaller_equal_0",
                      "Mean p-value",
                      highlight_models = c("CRF"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_0_y_limits,
                      boxplots = TRUE)

#p6 <- p6 + labs(y = expression(hat(P)[B](0)))

# Remove legends and title
p1 <- p1 + theme(legend.position = "none") + ggtitle("")
p2 <- p2 + theme(legend.position = "none") + ggtitle("")
p3 <- p3 + theme(legend.position = "none") + ggtitle("")
p4 <- p4 + theme(legend.position = "none") + ggtitle("")
p5 <- p5 + theme(legend.position = "none") + ggtitle("")
p6 <- p6 + theme(legend.position = "none") + ggtitle("")

# Create a legend by hand
legend_data <- data.frame(
  method = factor(c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"),
                  levels = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB")),
  x = 1:7,
  y = 1
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = method)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = model_colors,
                     labels = c(
                       "PO" = "PO",
                       "PO_splines" = expression(PO[complex]), # change
                       "GEE" = "GEE",
                       "GEE_splines" = expression(GEE[complex]), # change
                       "RF" = "RF",
                       "CRF" = "CRF",
                       "WB" = "WB"
                     )) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

Six_plots_grid <- (p1 + p2 + p3 + p4 + p5 + p6) + 
  plot_layout(ncol = 2, axes = "collect")

Six_plots_grid <- Six_plots_grid / 
  legend_plot + plot_layout(heights = c(10, 0)) # the number 10 is irrelevant, it's the same as long as no 0 is there

Six_methods_in_one <- plot_continuous(results_df,
                                      "Prob_ATE_smaller_equal_0",
                                      "Mean p-value",
                                      highlight_models = c(),
                                      show_points = FALSE,
                                      y_limits = Prob_ATE_smaller_equal_0_y_limits,
                                      boxplots = TRUE)

#Six_methods_in_one <- Six_methods_in_one + labs(y = expression(hat(P)[B](0)))

# rename legend stuff

Six_methods_in_one <- Six_methods_in_one +
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
  )

Prob_ATE_smaller_equal_0 <- list(Six_plots_grid, Six_methods_in_one)
all_plots[["Prob_ATE_smaller_equal_0"]] <- Prob_ATE_smaller_equal_0


# Create two plots for Prob_ATE_smaller_equal_true_treatment_effect

p1 <- plot_continuous(results_df %>% filter(model %in% c("PO", "WB")),
                      "Prob_ATE_smaller_equal_true_treatment_effect",
                      "Mean expression(hat(P)[B](ATE))",
                      highlight_models = c("PO"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_true_treatment_effect_y_limits,
                      boxplots = TRUE)

p1 <- p1 + labs(y = expression("Mean " * hat(P)[B](ATE)))

p2 <- plot_continuous(results_df %>% filter(model %in% c("PO_splines", "WB")),
                      "Prob_ATE_smaller_equal_true_treatment_effect",
                      "Mean expression(hat(P)[B](ATE))",
                      highlight_models = c("PO_splines"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_true_treatment_effect_y_limits,
                      boxplots = TRUE)

p2 <- p2 + labs(y = expression("Mean " * hat(P)[B](ATE)))

p3 <- plot_continuous(results_df %>% filter(model %in% c("GEE", "WB")),
                      "Prob_ATE_smaller_equal_true_treatment_effect",
                      "Mean expression(hat(P)[B](ATE))",
                      highlight_models = c("GEE"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_true_treatment_effect_y_limits,
                      boxplots = TRUE)

p3 <- p3 + labs(y = expression("Mean " * hat(P)[B](ATE)))

p4 <- plot_continuous(results_df %>% filter(model %in% c("GEE_splines", "WB")),
                      "Prob_ATE_smaller_equal_true_treatment_effect",
                      "Mean expression(hat(P)[B](ATE))",
                      highlight_models = c("GEE_splines"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_true_treatment_effect_y_limits,
                      boxplots = TRUE)

p4 <- p4 + labs(y = expression("Mean " * hat(P)[B](ATE)))

p5 <- plot_continuous(results_df %>% filter(model %in% c("RF", "WB")),
                      "Prob_ATE_smaller_equal_true_treatment_effect",
                      "Mean expression(hat(P)[B](ATE))",
                      highlight_models = c("RF"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_true_treatment_effect_y_limits,
                      boxplots = TRUE)

p5 <- p5 + labs(y = expression("Mean " * hat(P)[B](ATE)))

p6 <- plot_continuous(results_df %>% filter(model %in% c("CRF" ,"WB")),
                      "Prob_ATE_smaller_equal_true_treatment_effect",
                      "Mean expression(hat(P)[B](ATE))",
                      highlight_models = c("WB"),
                      show_points = FALSE,
                      y_limits = Prob_ATE_smaller_equal_true_treatment_effect_y_limits,
                      boxplots = TRUE)

p6 <- p6 + labs(y = expression("Mean " * hat(P)[B](ATE)))

# Remove legends and title
p1 <- p1 + theme(legend.position = "none") + ggtitle("")
p2 <- p2 + theme(legend.position = "none") + ggtitle("")
p3 <- p3 + theme(legend.position = "none") + ggtitle("")
p4 <- p4 + theme(legend.position = "none") + ggtitle("")
p5 <- p5 + theme(legend.position = "none") + ggtitle("")
p6 <- p6 + theme(legend.position = "none") + ggtitle("")

# Create a legend by hand
legend_data <- data.frame(
  method = factor(c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB"),
                  levels = c("PO", "PO_splines", "GEE", "GEE_splines", "RF", "CRF", "WB")),
  x = 1:7,
  y = 1
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = method)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = model_colors,
                     labels = c(
                       "PO" = "PO",
                       "PO_splines" = expression(PO[complex]), # change
                       "GEE" = "GEE",
                       "GEE_splines" = expression(GEE[complex]), # change
                       "RF" = "RF",
                       "CRF" = "CRF",
                       "WB" = "WB"
                     )) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

Six_plots_grid <- (p1 + p2 + p3 + p4 + p5 + p6) + 
  plot_layout(ncol = 2, axes = "collect")

Six_plots_grid <- Six_plots_grid / 
  legend_plot + plot_layout(heights = c(10, 0)) # the number 10 is irrelevant, it's the same as long as no 0 is there

Six_methods_in_one <- plot_continuous(results_df,
                                      "Prob_ATE_smaller_equal_true_treatment_effect",
                                      "Mean expression(hat(P)[B](ATE))",
                                      highlight_models = c(),
                                      show_points = FALSE,
                                      y_limits = Prob_ATE_smaller_equal_true_treatment_effect_y_limits,
                                      boxplots = TRUE)

Six_methods_in_one <- Six_methods_in_one + labs(y = expression("Mean " * hat(P)[B](ATE)))

# rename legend stuff

Six_methods_in_one <- Six_methods_in_one +
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
  )

Prob_ATE_smaller_equal_true_treatment_effect <- list(Six_plots_grid, Six_methods_in_one)
all_plots[["Prob_ATE_smaller_equal_true_treatment_effect"]] <- Prob_ATE_smaller_equal_true_treatment_effect

# Create plot that shows sqrt(n) * RMSE based on sample size (not log2)

results_df <- results_df %>%
  mutate(error_time_root_n = abs(bias) * sqrt(n))


RMSE_plot <- ggplot(results_df, aes(x = n, y = error_time_root_n)) + 
  stat_smooth(data = results_df,
              aes(color = model),
              method = "lm",
              se = FALSE,
              formula = y ~ rcs(x, c(125, 250, 500, 1000)),
              size = 1) +
  labs(
    title = "Smoothed RMSE Times Square Root of Sample Size",
    y = "SMAD",
    x = "Sample Size"
  ) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  theme_minimal()

# rename legend stuff

RMSE_plot <- RMSE_plot +
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
  )

all_plots[["RMSE_plot"]] <- list(RMSE_plot)  # wrap single plot in a list for consistency








return(all_plots)



}

