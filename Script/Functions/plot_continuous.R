# ========================================================================================================================================
# plot_continuous Function
# ========================================================================================================================================

plot_continuous <- function(data, variable, label,
                            highlight_models = NULL,
                            show_points = TRUE,
                            y_limits = NULL,
                            boxplots = FALSE,
                            model_colors = c(
                              "GEE" = "#1b9e77",       # Teal
                              "GEE_splines" = "#d95f02", # Orange
                              "PO" = "#7570b3",        # Purple
                              "PO_splines" = "#e7298a", # Pink
                              "WB" = "#e6ab02",        # Mustard
                              "RF" = "#66a61e",        # Green
                              "CRF" = "#A0522D"        # Brown
                            )) {

  
  # Load required packages
  if (!require("rms")) {install.packages("rms", dependencies=TRUE); library(rms)}
  if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE); library(tidyverse)}
  
  data$highlight <- ifelse(data$model %in% highlight_models, "highlight", "background")
  
  p <- ggplot(data, aes(x = log2(n), y = .data[[variable]]))
  
  # Boxplots for highlighted models (BEHIND all points/lines)
  if (boxplots) {
    box_data <- subset(data, highlight == "highlight" & n %in% c(125, 250, 500, 1000))
    box_data$log2_n <- log2(box_data$n)
    p <- p + geom_boxplot(
      data = box_data,
      aes(x = log2_n, y = .data[[variable]], group = interaction(log2_n, model), fill = model),
      alpha = 0.6,
      width = 0.12,
      color = "black",  # Use black borders for visibility
      position = position_dodge(width = 0.2),
      outlier.shape = NA  # Optional: hide outliers for cleaner plot
    )
  }
  
  # Background models
  if (show_points) {
    p <- p + geom_point(
      data = subset(data, highlight == "background"),
      aes(color = model),
      alpha = 0.12
    )
  }
  p <- p + stat_smooth(
    data = subset(data, highlight == "background"),
    aes(color = model),
    method = "lm",
    se = FALSE,
    formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
    size = 1
  )
  
  # Highlighted models
  if (show_points) {
    p <- p + geom_point(
      data = subset(data, highlight == "highlight"),
      aes(color = model),
      alpha = 0.2
    )
  }
  p <- p + stat_smooth(
    data = subset(data, highlight == "highlight"),
    aes(color = model, fill = model),
    method = "lm",
    se = TRUE,
    formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
    size = 1
  )
  
  p <- p +
    labs(
      title = paste(label, "by Sample Size"),
      y = label,
      x = "Sample Size"
    ) +
    scale_x_continuous(
      breaks = log2(c(125, 250, 500, 1000)),
      labels = c("125", "250", "500", "1000")
    ) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors) +
    theme_minimal()
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}




### EXAMPLES ###

# plot_continuous(results_df %>% filter(model %in% c("PO", "WB")),
#                 "treatment_effect",
#                 "Estimated Treatment Effect",
#                 highlight_models = c("PO"),
#                 show_points = FALSE,
#                 y_limits = c(results_df$true_treatment_effect[1] - 0.7, results_df$true_treatment_effect[1] + 0.7),
#                 boxplots = TRUE) +
#   geom_abline(slope = 0, intercept = first(results_df$true_treatment_effect))
# 
# plot_continuous(results_df %>% filter(model %in% c("PO", "WB")),
#                 "Prob_ATE_smaller_equal_true_treatment_effect",
#                 "Probability that the ATE is Smaller or Equal to the True ATE",
#                 highlight_models = c("PO"),
#                 show_points = FALSE,
#                 y_limits = c(0, 1),
#                 boxplots = TRUE)