# ========================================================================================================================================
# plot_logistic Function
# ========================================================================================================================================

plot_logistic <- function(data, variable, label,
                          highlight_models = NULL,
                          show_points = TRUE,
                          y_limits = c(0, 1),
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
  
  p <- ggplot(data, aes(x = log2(n), y = as.numeric(.data[[variable]])))
  
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
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE,
    formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
    size = 1.5
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
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE,
    formula = y ~ rcs(x, log2(c(125, 250, 500, 1000))),
    size = 1
  )
  
  p <- p +
    labs(
      title = paste("Probability of", label, "by Sample Size"),
      y = paste("P(", label, ")"),
      x = "Log2-Transformed Sample Size"
    ) +
    scale_x_continuous(
      breaks = log2(c(125, 250, 500, 1000)),
      labels = paste0(round(log2(c(125, 250, 500, 1000)), digits = 1), "(n=", c("125", "250", "500", "1000"), ")")
    ) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors) +
    theme_minimal()
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}


### EXAMPLE ###
# plot_logistic(results_df %>% filter(model %in% c("RF", "WB")),
#               "CI_excludes_0",
#               "CI Excludes 0",
#               highlight_models = c("RF"),
#               show_points = FALSE)