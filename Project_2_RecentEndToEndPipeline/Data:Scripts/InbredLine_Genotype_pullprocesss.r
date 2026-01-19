# Load required libraries 
library(arrow)
library(dplyr)
library(tidyverse)

# Pull the genotypes from the whole folder of partitioned parquet files
ds = open_dataset("s3://xxxxxxxx/", format = "parquet")

result = ds %>%
  collect()

parents = unique(c(Pheno_final$mother_plant,
                    Pheno_final$father_plant))

# Read test genotypes
Test_data = read.csv("xxxxx.csv", header = FALSE)
Test = data.frame(genotype = Test_data[[1]])
TestIDs = unique(as.factor(Test$genotype))

AllIDs   =   union(TrainIDs, TestIDs)

result_subset = result %>%
  filter(samplename %in% AllIDs)

# Inputs: results_subset (rows = samples, cols = markers + metadata), Pheno_final
id_col = "samplename"

# Marker columns = all numeric columns (adjust if you have numeric metadata you want to exclude)
marker_cols = names(result_subset)[vapply(result_subset, is.numeric, logical(1))]
if (!length(marker_cols)) stop("No numeric marker columns detected in results_subset.")

# Deduplicate samples after cleaning names
rs1 = result_subset %>%
  rowwise() %>%
  dplyr::mutate(n_missing = sum(is.na(c_across(all_of(marker_cols))))) %>%
  ungroup() %>%
  group_by(samplename) %>%
  slice_min(order_by = n_missing, n = 1, with_ties = FALSE) %>%  # pick least missing
  ungroup() %>%
  dplyr::select(-n_missing)

# Build genotype matrix: rows = clean sample IDs, cols = markers
G = rs1 %>%
  dplyr::select(samplename, all_of(marker_cols)+1) %>%
  column_to_rownames("samplename") %>%
  as.matrix()
mode(G) = "numeric"
