# Load required libraries 
library(dplyr)
library(ggplot2)

# GEBV result and model performance 
analyze_results = function(results, variance_data = NULL) {
  # Bind results into one data frame
  results_df <- bind_rows(results)
  
  # Correlation per trait
  trait_cor = results_df %>%
    group_by(Trait, Model) %>%
    summarise(Correlation = cor(Observed, Predicted, use = "complete.obs"),
              .groups = "drop")
  
  # Compute r per trait (for labels)
  reg_labels = results_df %>%
    group_by(Trait, Model) %>%
    summarise(
      r = cor(Observed, Predicted, use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("r = ", r),
      x = -Inf,  
      y = Inf
    )


 # Plot the correlation 
 plot = ggplot(results_df, aes(x = Predicted, y = Observed, color = Model)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
  facet_grid(Model ~ Trait, scales = "free") +
  geom_text(
    data = reg_labels,
    aes(x = x, y = y, label = label, color = Model),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1.5,
    size = 4, fontface = "italic"
  ) +
  labs(
    x = "Predicted mid-parent trait value",
    y = "Observed trait value",
    title = "Inbred model prediction model accuracies"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

return(list(
  results_df = results_df,
  trait_cor = trait_cor,
  reg_labels = reg_labels,
  plot = plot
))
}

# 5fold CV results and model performance 
analyze_5foldcv_results = function(results_list_CV, model_name) {
  # Combine fold-level correlations and filter for the specified model
  fold_cor_all <- do.call(rbind, lapply(results_list_CV, function(x) {
    df <- x[["fold_cor_df"]]
    df[df$Model == model_name, ]
  }))
  
  # Combine mean correlations and filter for the specified model
  mean_cor_all <- do.call(rbind, lapply(results_list_CV, function(x) {
    df <- x[["mean_cor_df"]]
    df[df$Model == model_name, ]
  }))
  
  # Plot average 5-fold CV correlation per trait
  plot_mean = ggplot(mean_cor_all, aes(x = Trait, y = MeanMidParentCorrelation)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = round(MeanMidParentCorrelation, 2)), vjust = -0.5) +
    ylim(0, 1) +
    labs(
      title = paste("Average 5-Fold CV Correlation per Trait -", model_name),
      x = "Trait",
      y = "Mean Pearson Correlation (GEBV vs Observed)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot fold-level CV correlations per trait
  plot_folds = ggplot(fold_cor_all, aes(x = factor(Fold), y = MidParentCorrelation, fill = Trait)) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_text(aes(label = round(MidParentCorrelation, 2)), 
              position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
    ylim(0, 1) +
    labs(
      title = paste("Fold-Level 5-Fold CV Correlations per Trait -", model_name),
      x = "Fold",
      y = "Pearson Correlation (GEBV vs Observed)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 11),
      legend.position = "bottom"
    )
  
  return(list(
    fold_cor_all = fold_cor_all,
    mean_cor_all = mean_cor_all,
    plot_mean = plot_mean,
    plot_folds = plot_folds
  ))
}
