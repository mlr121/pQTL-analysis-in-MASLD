compute_gam_nocov <- function(all_pheno, pheno_cols) {
  # Requires mgcv
  if (!"mgcv" %in% .packages(all.available = TRUE)) {
    suppressPackageStartupMessages(require(mgcv))
  } else {
    suppressPackageStartupMessages(library(mgcv))
  }
  
  # Basic checks
  if (!all(pheno_cols %in% colnames(all_pheno))) {
    missing_cols <- pheno_cols[!pheno_cols %in% colnames(all_pheno)]
    stop("These pheno_cols are not in all_pheno: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Identify numeric columns
  numeric_cols <- colnames(all_pheno)[sapply(all_pheno, is.numeric)]
  
  # Proteins = numeric columns that are NOT phenotypes
  protein_cols <- setdiff(numeric_cols, pheno_cols)
  
  if (length(protein_cols) == 0)
    stop("No protein columns detected. Check pheno_cols selection.")
  
  results_list <- vector("list", length = length(pheno_cols))
  names(results_list) <- pheno_cols
  
  # Loop through phenotypes
  for (trait in pheno_cols) {
    message("Processing trait: ", trait)
    y <- all_pheno[[trait]]
    
    trait_res <- lapply(protein_cols, function(prot) {
      protein_vals <- all_pheno[[prot]]
      
      # Require >= 3 paired non-missing observations
      if (sum(!is.na(protein_vals) & !is.na(y)) < 3) return(NULL)
      
      df_model <- data.frame(
        trait = y,
        protein = protein_vals
      )
      
      fit <- tryCatch({
        mgcv::gam(trait ~ s(protein), data = df_model)
      }, error = function(e) return(NULL))
      
      if (is.null(fit)) return(NULL)
      
      s_summary <- summary(fit)
      
      pval <- tryCatch({
        s_summary$s.table["s(protein)", "p-value"]
      }, error = function(e) NA_real_)
      
      dev_expl <- tryCatch({
        s_summary$dev.expl
      }, error = function(e) NA_real_)
      
      tau_val <- tryCatch({
        cor(df_model$protein, df_model$trait,
            method = "kendall", use = "complete.obs")
      }, error = function(e) NA_real_)
      
      direction <- if (is.na(tau_val)) {
        NA_character_
      } else if (tau_val > 0) {
        "1"
      } else if (tau_val < 0) {
        "-1"
      } else {
        "neutral"
      }
      
      data.frame(
        Protein = prot,
        p_value = pval,
        deviance_explained = dev_expl,
        tau = tau_val,
        direction = direction,
        stringsAsFactors = FALSE
      )
    })
    
    trait_df <- do.call(rbind, trait_res)
    results_list[[trait]] <- if (!is.null(trait_df) && nrow(trait_df) > 0) trait_df else NULL
  }
  
  return(results_list)
}
