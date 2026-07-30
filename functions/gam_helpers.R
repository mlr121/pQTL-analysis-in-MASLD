################################################################################
# Script: gam_helpers.R
# Description: Functions to compute GAM associations and Kendall's tau.
################################################################################

library(mgcv)

compute_gam_associations <- function(all_pheno, pheno_cols) {
  if (!all(pheno_cols %in% colnames(all_pheno))) {
    missing_cols <- pheno_cols[!pheno_cols %in% colnames(all_pheno)]
    stop("These pheno_cols are not in all_pheno: ", paste(missing_cols, collapse = ", "))
  }
  if (!"Sex" %in% colnames(all_pheno)) stop("all_pheno must contain a 'Sex' column.")
  
  sex_factor <- as.factor(all_pheno$Sex)
  numeric_cols <- colnames(all_pheno)[sapply(all_pheno, is.numeric)]
  protein_cols <- setdiff(numeric_cols, pheno_cols)
  
  if (length(protein_cols) == 0) stop("No protein columns detected.")
  
  results_list <- vector("list", length = length(pheno_cols))
  names(results_list) <- pheno_cols
  
  for (trait in pheno_cols) {
    message("Processing trait: ", trait)
    y <- all_pheno[[trait]]
    
    trait_res <- lapply(protein_cols, function(prot) {
      protein_vals <- all_pheno[[prot]]
      if (sum(!is.na(protein_vals) & !is.na(y)) < 3) return(NULL)
      
      df_model <- data.frame(trait = y, sex = sex_factor, protein = protein_vals)
      
      fit <- tryCatch({
        mgcv::gam(trait ~ s(protein) + sex, data = df_model)
      }, error = function(e) return(NULL))
      
      if (is.null(fit)) return(NULL)
      
      s_summary <- summary(fit)
      pval <- tryCatch({ s_summary$s.table["s(protein)", "p-value"] }, error = function(e) NA_real_)
      dev_expl <- tryCatch({ s_summary$dev.expl }, error = function(e) NA_real_)
      tau_val <- tryCatch({
        cor(df_model$protein, df_model$trait, method = "kendall", use = "complete.obs")
      }, error = function(e) NA_real_)
      
      direction <- if (is.na(tau_val)) NA_character_ else if (tau_val > 0) "1" else if (tau_val < 0) "-1" else "neutral"
      
      data.frame(
        Protein = prot, p_value = pval, deviance_explained = dev_expl,
        tau = tau_val, direction = direction, stringsAsFactors = FALSE
      )
    })
    
    trait_df <- do.call(rbind, trait_res)
    results_list[[trait]] <- if (!is.null(trait_df) && nrow(trait_df) > 0) trait_df else NULL
  }
  return(results_list)
}

compute_ubr1_associations <- function(all_pheno, ubr1_col = "Ubr1") {
  if (!ubr1_col %in% colnames(all_pheno)) stop("ubr1_col not found in all_pheno: ", ubr1_col)
  if (!"Sex" %in% colnames(all_pheno)) stop("all_pheno must contain a 'Sex' column.")
  
  sex_factor <- as.factor(all_pheno$Sex)
  numeric_cols <- colnames(all_pheno)[sapply(all_pheno, is.numeric)]
  protein_cols <- setdiff(numeric_cols, ubr1_col)
  
  if (length(protein_cols) == 0) stop("No protein outcome columns detected.")
  
  ubr1_vals <- all_pheno[[ubr1_col]]
  message("Computing associations for ", ubr1_col)
  
  results <- lapply(protein_cols, function(prot) {
    y <- all_pheno[[prot]]
    if (sum(!is.na(ubr1_vals) & !is.na(y)) < 3) return(NULL)
    
    df_model <- data.frame(outcome = y, sex = sex_factor, ubr1 = ubr1_vals)
    
    fit <- tryCatch({ mgcv::gam(outcome ~ s(ubr1) + sex, data = df_model) }, error = function(e) return(NULL))
    if (is.null(fit)) return(NULL)
    
    s_summary <- summary(fit)
    pval <- tryCatch({ s_summary$s.table["s(ubr1)", "p-value"] }, error = function(e) NA_real_)
    dev_expl <- tryCatch({ s_summary$dev.expl }, error = function(e) NA_real_)
    tau_val <- tryCatch({
      cor(df_model$ubr1, df_model$outcome, method = "kendall", use = "complete.obs")
    }, error = function(e) NA_real_)
    
    direction <- if (is.na(tau_val)) NA_character_ else if (tau_val > 0) "1" else if (tau_val < 0) "-1" else "neutral"
    
    data.frame(
      Protein = prot, p_value = pval, deviance_explained = dev_expl,
      tau = tau_val, direction = direction, stringsAsFactors = FALSE
    )
  })
  return(do.call(rbind, results))
}