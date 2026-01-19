# Load required libraries 
library(dplyr)
library(ggplot2)
library(reshape2)

# Function to compute and plot across-environment trait correlations
plot_trait_corr_across_env = function(data, traits, Output_dir = "Output/") {
  
  # Compute across-environment means per genotype
  mean_data = data %>%
    group_by(genotype) %>%
    summarise(across(all_of(traits), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  # Compute correlation matrix
  cor_mat = cor(mean_data[ , traits], use = "pairwise.complete.obs", method = "pearson")
  
  melted = melt(cor_mat, na.rm = TRUE)
  # Retain lower triangle
  melted = melted[as.numeric(melted$Var1) > as.numeric(melted$Var2), ] 
  
  # Plot correlation heatmap
  p = ggplot(melted, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          panel.grid = element_blank()) +
    labs(title = "Pearson Correlations Between Trait Means Across Environments")
  
  # Save plot to PDF
  pdf_path = file.path(Output_dir, "Trait_Correlations_AcrossEnvironemtMeans.pdf")
  pdf(pdf_path, width = 7, height = 6)
  print(p)
  dev.off()

  return(list(cor_matrix = cor_mat, plot = p))
}

# Function to compute and plot environment correlations per trait
plot_env_corr_per_traits = function(data, traits, Output_dir = "Output/") {
  
  # Per trait environment correlation heatmap
  plot_env_corr_one = function(tr) {
    cat("Processing:", tr, "\n")
    
    mat = data %>%
      dplyr::select(genotype, environment, all_of(tr)) %>%
      rename(value = !!sym(tr)) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::group_by(genotype, environment) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = environment, values_from = value)
    
    if (ncol(mat) <= 2) return(NULL)  
    
    cor_mat = cor(mat[ , -1], use = "pairwise.complete.obs", method = "pearson")
    melted = melt(cor_mat, na.rm = TRUE)
    # Retain lower triangle
    melted <- melted[as.numeric(melted$Var1) > as.numeric(melted$Var2), ]  
    
    ggplot(melted, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, limits = c(-1, 1)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank(),
            panel.grid = element_blank()) +
      labs(title = paste("Environment Correlations -", tr))
  }
  
  # Generate all plots
  plots = lapply(traits, plot_env_corr_one)

  # Save all to one PDF
  pdf_path = file.path(Output_dir, "Env_Correlations_per_trait.pdf")
  pdf(pdf_path, width = 7, height = 6)
  for (p in plots) print(p)
  dev.off()
}
 
