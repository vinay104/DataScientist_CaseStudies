# Load required libraries 
library(arrow)
library(parallel)
library(dplyr)
library(sommer)
library(tidyverse)

# Load phenotypes
source("Phenotype_pull_process.R")
source("InbredLine_Genotype_pullprocesss.r")
source("InbredLine_TraitProcessing&ModelPerformance_Function.r")
source("InbredLine_extrractresults_n_plotting.r")

# Results output folder 
Output_dir = "Output/"

# Assign and split cores for parallel processing
n_cores = detectCores()
cl_cv = makeCluster(floor(n_cores / 2))
cl_gebv = makeCluster(floor(n_cores / 2))

# Load packages and objects
clusterEvalQ(cl_cv,   { library(arrow); library(dplyr); library(sommer); library(tidyverse) })
clusterEvalQ(cl_gebv, { library(arrow); library(dplyr); library(sommer); library(tidyverse) })

clusterExport(cl_cv,   c("shared_objects_5fold_CV", "process_trait_worker_CV"))
clusterExport(cl_gebv, c("process_trait_worker", "process_trait", "shared_objects_GEBVs"))

cv_future = parallel:::mcparallel(parLapply(cl_cv, traits, process_trait_worker_CV, shared = shared_objects_5fold_CV))
gebv_future = parallel:::mcparallel(parLapply(cl_gebv, traits, process_trait_worker, shared = shared_objects_GEBVs))

# Collect results
results_list_CV = parallel:::mccollect(cv_future)[[1]]
results_list_GEBV = parallel:::mccollect(gebv_future)[[1]]

stopCluster(cl_cv)
stopCluster(cl_gebv)

# Row bind the results and pivot them into wide-format for GEBVs
results_df_GEBV <- do.call(rbind, lapply(results_list_GEBV, `[[`, "BLUPs"))

results_df_GEBV_GxE = results_df_GEBV %>%
  dplyr::filter(Model == "GxE")

GEBVs_wide_GxE = results_df_GEBV_GxE %>%
  tidyr::pivot_wider(
    id_cols = genotype,
    names_from = Trait,
    values_from = GEBV
  ) %>%
  dplyr::filter(genotype %in% TestIDs)

write.csv(GEBVs_wide_GxE,
          file = file.path(Output_dir, "InbredLine_PredictedValues_GxE.csv"),
          row.names = FALSE)

# Combine variance/heritability
variance_results <- get_variance_results(results_list_GEBV)

# Model Performance (No-fold CV)
models = unique(results_df_GEBV$Model)

corr_plots <- list()

for (m in models) {
  res_df <- results_df_GEBV %>% dplyr::filter(Model == m)
  results <- get_results(
    traits = traits,
    Pheno_final = Pheno_final,
    common_cols = common_cols,
    results_df_GEBV = res_df,
    TrainIDs = TrainIDs
  )
  
  Model_performance <- analyze_results(results, variance_data = variance_results)
  
  corr_plots[[m]] <- Model_performance$plot + ggtitle("InbredLineTraitCorrelations", m)
  
  write.csv(Model_performance$trait_cor, file.path(Output_dir, paste0(m, "_InbredLineTraitCorrelations.csv")), row.names = FALSE)
  write.csv(Model_performance$all_variance, file.path(Output_dir, paste0(m, "_InbredLineVarianceComponents.csv")), row.names = FALSE)
  write.csv(Model_performance$variance_summary, file.path(Output_dir, paste0(m, "_InbredLineVarianceHeritability.csv")), row.names = FALSE)
  write.csv(Model_performance$merged_summary, file.path(Output_dir, paste0(m, "_InbredLineCompleteSummary.csv")), row.names = FALSE)
}

TraitCorrelation_per_model <- marrangeGrob(grobs = corr_plots, nrow = 1, ncol = 1)
ggsave(filename = file.path(Output_dir, "AllModels_InbredLineTraitCorrelations.pdf"), TraitCorrelation_per_model, width = 8, height = 5)

# 5foldCV model performance 
all_results = lapply(models, function(m) {
  CV_ModelPerformance = analyze_5foldcv_results(results_list_CV, model_name = m)
  
  # Save plots into PDF for each model
  pdf(file = file.path(Output_dir,paste0(m, "model_Inbredline5FoldCV_results_", ".pdf")), width = 10, height = 6)
  print(CV_ModelPerformance$plot_mean)
  print(CV_ModelPerformance$plot_folds)
  dev.off()
  
  # Save fold level and mean correlations of 5Fold CV
  write.csv(CV_ModelPerformance$fold_cor_all, file = file.path(Output_dir, paste0(m,"model_Inbredline5FoldCV_fold_results_", ".csv")), row.names = FALSE)
  write.csv(CV_ModelPerformance$mean_cor_all, file = file.path(Output_dir, paste0(m, "model_Inbredline5FoldCV_mean_results_", ".csv")), row.names = FALSE)
  
  
  return(CV_ModelPerformance)
})
