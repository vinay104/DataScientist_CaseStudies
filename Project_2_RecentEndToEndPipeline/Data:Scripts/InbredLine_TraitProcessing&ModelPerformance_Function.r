# Read train and test genotypes
TrainIDs = unique(as.factor(c(Pheno_final$male, Pheno_final$female))) %>% droplevels(c())
Test_data = read.csv("xxxx.csv", header = FALSE)
Test = data.frame(genotype = Test_data[[1]])
TestIDs = unique(as.factor(Test$genotype))

AllIDs   =   union(TrainIDs, TestIDs)
common_cols = intersect(AllIDs, row.names(G))

Geno_final   = G[common_cols, drop = FALSE, ]
Geno_center = scale(Geno_final, center = TRUE, scale = FALSE)
GRM    =  A.mat(Geno_center)

# Build GRM using training set genotypes 
IDs = unique(c(Pheno_final$male, Pheno_final$female))

TrainIDsCV = intersect(IDs, rownames(G))

G_train = G[rownames(G) %in% TrainIDsCV, , drop = FALSE]
G_train = G_train[match(TrainIDsCV, rownames(G_train)), , drop = FALSE]

Geno_center_train = scale(G_train, center = TRUE, scale = FALSE)
GRM_TrainIDs = A.mat(Geno_center_train)
GRM_TrainIDs = GRM_TrainIDs[TrainIDsCV, TrainIDsCV]

# Trait list
traits = c("DMC", 
           "DOR" ,
           "EMY28PLUS",
           "EMY35PLUS",
           "EMY40PLUS",
           "EMY40_60",
           "EMY60PLUS",
           "MBD",
           "RT",
           "SG",
           "STC",
           "STY",
           "TDMC",
           "TN",
           "TY",
           "TV",
           "UWW") 

### Shared objects for 5-fold CV
shared_objects_5fold_CV = list(
  Pheno_final = Pheno_final,
  GRM_TrainIDs = GRM_TrainIDs,
  TrainIDsCV = TrainIDsCV,
  common_cols = common_cols
)

#### Save shared objects
shared_objects_GEBVs = list(
  Pheno_final = Pheno_final,
  GRM = GRM,
  common_cols = common_cols,
  TestIDs = TestIDs,
  AllIDs = AllIDs
)

######################################
process_trait = function(tr) {
  cat("\n===== Processing trait:", tr, "=====\n")
  
  # Prepare phenotypes
  Pheno_F = Pheno_final %>%
    dplyr::mutate(Trait = .data[[tr]]) %>%
    dplyr::filter(female %in% common_cols, male %in% common_cols) %>%
    dplyr::select(genotype, female, male, Trait, environment, Year, loc) %>%
    arrange(environment)
  
  
  Pheno_F = Pheno_F %>%
    dplyr::mutate(
      env_female = interaction(environment, female, drop = TRUE),
      env_male = interaction(environment, male, drop = TRUE),
      year_female = interaction(Year, female, drop = TRUE),
      year_male = interaction(Year, male, drop = TRUE),
      location_female = interaction(loc, female, drop = TRUE),
      location_male = interaction(loc, male, drop = TRUE),
      units = as.factor(seq_len(n()))
    )
  
  blups_all = list()
  variance_all = list() 
  
  # Fit the model
  fit_model = function(formula, random, rcov) {
    M = tryCatch(
      mmes(
        formula,
        random = random,
        rcov   = rcov,
        data   = Pheno_F,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(M) || (!is.null(M$singular) && M$singular)) {
      M = tryCatch(
        mmes(
          formula,
          random = random,
          rcov   = rcov,
          data   = Pheno_F,
          tolParInv = 1000,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
    }
    return(M)
  }

  models_to_fit = list(
    G_only = list(
      formula = Trait ~ environment,
      random = ~ vsm(ism(overlay(female, male)), Gu = GRM) +
        genotype,
      rcov    = ~ vsm(dsm(environment), ism(units))
    ),
    GxE = list(
      formula = Trait ~ environment,
      random = ~ vsm(ism(overlay(female, male)), Gu = GRM) +
        vsm(ism(overlay(env_female, env_male))) +
        genotype,      
      rcov    = ~ vsm(dsm(environment), ism(units))
    )
  )
  
 # Extract BLUPs
  extract_blups = function(Model) {
    if (is.null(Model)) return(NULL)
    U = Model$u
    blups = data.frame(
      genotype = rownames(U),
      Predicted_PHENO = U[,1]
    )
    blups = blups[blups$genotype %in% AllIDs, ]
    
    b = Model$b
    intercepts = b[1] + b[2:length(b)]
    Intercept = mean(intercepts)
    blups$GEBV = Intercept + blups$Predicted_PHENO
    return(blups)
  }
  
  for (mname in names(models_to_fit)) {
    cat("  Fitting model:", mname, "\n")
    m = models_to_fit[[mname]]
    Model = fit_model(m$formula, m$random, m$rcov)
    if (is.null(Model) || is.null(Model$theta)) next
    
    vc = Model$theta
    n_comp = length(vc)
    if (n_comp == 0) next
    
    # Residual is always last
    residual_var = mean(vc[[n_comp]], na.rm = TRUE)
    genetic_vars = vc[1:(n_comp - 1)]
    genetic_v = as.numeric(unlist(genetic_vars)[which(!is.na(as.numeric(unlist(genetic_vars))))[1]])
    genetic_sum = sum(unlist(genetic_vars), na.rm = TRUE)
    
    # Heritability
    heritability = if (!is.na(residual_var) && (genetic_sum + residual_var) && (genetic_v) > 0) {
      genetic_v / (genetic_sum + residual_var)
    } else {
      NA
    }
  
    # Store variance components
    vc_df = data.frame(
      Trait = tr,
      Model = mname,
      residual_var = residual_var,
      heritability = heritability
    )
    
    for (i in 1:3) {
      vc_df[[paste0("genetic_var", i)]] <- if (i <= length(genetic_vars)) genetic_vars[[i]] else NA
    }
    
    variance_all[[mname]] <- vc_df
    
    # Extract BLUPs
    blups = extract_blups(Model)
    if (!is.null(blups)) {
      blups$Trait = tr
      blups$Model = mname
      blups_all[[mname]] = blups
    }
  }
  
  if (length(variance_all) > 0) {
    variance_all <- lapply(variance_all, function(df) {
      df <- df[, c("Trait", "Model", "genetic_var1", "genetic_var2", "genetic_var3", "residual_var", "heritability")]
      return(df)
    })
    variance_all <- do.call(rbind, variance_all)
  } else {
    variance_all <- NULL
  }
  
  if (length(blups_all) > 0) {
    blups_all = lapply(blups_all, function(df) {
      df = df[, c("genotype", "Predicted_PHENO", "GEBV", "Trait", "Model")]
      return(df)
    })
    blups_all = do.call(rbind, blups_all)
  } else {
    blups_all = NULL
  }
  
  return(list(
    BLUPs = blups_all,
    Variance = variance_all
  ))
  
}

# Parallel execution per trait 
process_trait_worker = function(tr, shared) {
  attach(shared, warn.conflicts = FALSE)
  res = process_trait(tr)
  detach(shared)
  return(res)
}

# Extract the results 
# Combine variance components across traits 
get_variance_results = function(process_results) {
  all_variance = do.call(rbind, lapply(process_results, `[[`, "Variance"))
  
  variance_summary = all_variance %>%
    group_by(Trait, Model) %>%
    summarise(
      across(starts_with("genetic_var"), ~ mean(.x, na.rm = TRUE)),
      residual_var = mean(residual_var, na.rm = TRUE),
      heritability = mean(heritability, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    All_Variance = all_variance,
    Variance_Summary = variance_summary
  )
}


################# Trait worker function for CV
process_trait_worker_CV = function(tr, shared, k = 5) {
  attach(shared, warn.conflicts = FALSE)
  
  parent_ids = TrainIDsCV
  hybrid_ids = unique(Pheno_final$genotype)   
  
  folds = split(sample(hybrid_ids), rep(1:k, length.out = length(hybrid_ids)))
  
  # Fit model with fallback for singularity 
  fit_model = function(formula, random, rcov, Pheno_F) {
    M = tryCatch(
      mmes(formula, 
           random = random, 
           rcov = rcov, 
           data = Pheno_F, 
           verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(M) || (!is.null(M$singular) && M$singular)) {
      M = tryCatch(
        mmes(formula, 
             random = random, 
             rcov = rcov, 
             data = Pheno_F, 
             tolParInv = 1000, 
             verbose = FALSE),
        error = function(e) NULL
      )
    }
    return(M)
  }
  
  # Calculate mid-parent correlation 
  get_midparent_cor = function(model, Pheno_F, fold_ids) {
    if (is.null(model)) return(NA)
    
    U = model$u
    if (is.null(U) || nrow(U) == 0) return(NA)
    
    blups = data.frame(effect = rownames(U), Predicted_PHENO = U[,1])
    blups$genotype = sub(".*:", "", blups$effect)
    
    # Compute intercept (mean of fixed effects)
    compute_intercept = function(model, n) {
      b = model$b
      intercepts = b[1] + b[2:n]
      mean(intercepts)
    }
    
    Intercept = compute_intercept(model, length(model$b))
    blups$GEBV = Intercept + blups$Predicted_PHENO
    
    # Restrict to fold hybrids
    Pheno_fold = Pheno_F %>% filter(genotype %in% fold_ids)
    
    midparent_df = Pheno_fold %>%
      left_join(blups %>% rename(female = genotype, GEBV_female = GEBV), by = "female") %>%
      left_join(blups %>% rename(male   = genotype, GEBV_male   = GEBV), by = "male") %>%
      mutate(MidParentValue = (GEBV_female + GEBV_male)/2)
    
    # Observed trait values
    obs_df = Pheno_final %>%
      group_by(genotype, female, male) %>%
      summarise(Observed = mean(.data[[tr]], na.rm = TRUE), .groups = "drop")
    
    midparent_df = left_join(midparent_df, obs_df, by = c("genotype","female","male"))
    
    if (nrow(midparent_df) == 0) return(NA)
    
    cor_val = suppressWarnings(cor(midparent_df$MidParentValue, midparent_df$Observed, use="complete.obs"))
    if (is.na(cor_val) || length(cor_val) == 0) return(NA)
    return(cor_val)
  }
  
  # CV folds 
  fold_cor_list = lapply(seq_along(folds), function(i) {
    fold_ids = folds[[i]]
    
    # Prepare phenotype data
    Pheno_F = Pheno_final %>%
      mutate(Trait = .data[[tr]]) %>%
      filter(female %in% TrainIDsCV, male %in% TrainIDsCV) %>%
      mutate(
        mask = genotype %in% fold_ids,
        TraitMasked = ifelse(mask, NA, Trait)
      ) %>%
      dplyr::select(genotype, female, male, TraitMasked, environment, Year, loc) %>%
      arrange(environment) %>%
      mutate(
        env_female      = interaction(environment, female, drop = TRUE),
        env_male        = interaction(environment, male, drop = TRUE),
        year_female     = interaction(Year, female, drop = TRUE),
        year_male       = interaction(Year, male, drop = TRUE),
        location_female = interaction(loc, female, drop = TRUE),
        location_male   = interaction(loc, male, drop = TRUE),
        units           = as.factor(seq_len(n())),
        geno_env        = interaction(environment, genotype, drop = TRUE)
      )
    
    if (nrow(Pheno_F) == 0 || length(unique(Pheno_F$female)) < 2 || length(unique(Pheno_F$male)) < 2) {
      return(data.frame(Trait = tr, Fold = i, Model = NA, MidParentCorrelation = NA))
    }
    
    # Fit models 
    Model_gxe = fit_model(
      TraitMasked ~ environment,
      random = ~ vsm(ism(overlay(female, male)), Gu = GRM_TrainIDs) +
        vsm(ism(overlay(env_female, env_male))) +
        genotype,
      rcov   = ~ vsm(dsm(environment), ism(units)),
      Pheno_F = Pheno_F
    )
    
    Model_noGxE = fit_model(
      TraitMasked ~ environment,
      random = ~ vsm(ism(overlay(female, male)), Gu = GRM_TrainIDs) + 
        genotype,
      rcov   = ~ vsm(dsm(environment), ism(units)),
      Pheno_F = Pheno_F
    )
    
    # Compute correlations
    rbind(
      data.frame(Trait = tr, Fold = i, Model = "GxE",
                 MidParentCorrelation = get_midparent_cor(Model_gxe, Pheno_F, fold_ids)),
      data.frame(Trait = tr, Fold = i, Model = "Baseline",
                 MidParentCorrelation = get_midparent_cor(Model_noGxE, Pheno_F, fold_ids))
    )
  })
  
  fold_cor_df = do.call(rbind, fold_cor_list)
  
  mean_cor_df = fold_cor_df %>%
    group_by(Trait, Model) %>%
    summarise(MeanMidParentCorrelation = mean(MidParentCorrelation, na.rm = TRUE), .groups = "drop")
  
  detach(shared)
  
  list(fold_cor_df = fold_cor_df, mean_cor_df = mean_cor_df)
}

# Hybrid performance function
get_results = function(traits, 
                        Pheno_final, 
                        common_cols, 
                        results_df_GEBV, 
                        TrainIDs) {
  
  results = lapply(traits, function(tr) {
    # Subset phenotype data for the trait
    Pheno_F = Pheno_final %>%
      mutate(TraitValue = .data[[tr]]) %>%
      filter(female %in% common_cols, male %in% common_cols) %>%
      dplyr::select(genotype, female, male, TraitValue, environment, Year, loc) %>%
      arrange(environment)
    
    # Average observed values per hybrid
    Observed = Pheno_F %>%
      group_by(genotype, female, male) %>%
      summarise(
        Observed = mean(TraitValue, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Trait = tr)
    
    # Get GEBVs for training IDs
    GEBVs_train <- results_df_GEBV %>%
      dplyr::filter(genotype %in% TrainIDs, Trait == tr) %>%
      dplyr::select(genotype, Model, GEBV)
    
    # Loop over all models
    merged_list <- lapply(unique(GEBVs_train$Model), function(mdl) {
      GEBVs_model <- GEBVs_train %>%
        dplyr::filter(Model == mdl) %>%
        dplyr::select(genotype, GEBV)
    
    # Join parent GEBVs for female and male
      hybrids <- Observed %>%
        dplyr::left_join(GEBVs_model, by = c("female" = "genotype")) %>%
        dplyr::rename(pred_p1 = GEBV) %>%
        dplyr::left_join(GEBVs_model, by = c("male" = "genotype")) %>%
        dplyr::rename(pred_p2 = GEBV) %>%
        dplyr::mutate(
          pred_p1 = as.numeric(pred_p1),
          pred_p2 = as.numeric(pred_p2),
          mid_parent = (pred_p1 + pred_p2)/2,
          Model = mdl
        )
    
      # Predicted values
      Predicted <- hybrids %>%
        dplyr::select(genotype, mid_parent, Model) %>%
        dplyr::rename(Predicted = mid_parent) %>%
        dplyr::mutate(Trait = tr)
      
      # Keep only common genotypes
      common_ids <- intersect(as.factor(Observed$genotype), Predicted$genotype)
      Observed_filtered <- Observed %>% dplyr::filter(genotype %in% common_ids)
      Predicted_filtered <- Predicted %>% dplyr::filter(genotype %in% common_ids)
      
      # Merge observed and predicted
      merged <- Observed_filtered %>%
        dplyr::left_join(Predicted_filtered, by = c("genotype", "Trait"))
      
      return(merged)
    })
    
    # Combine all models for this trait
    do.call(rbind, merged_list)
  })
  
  # Combine all traits
  do.call(rbind, results)
}
