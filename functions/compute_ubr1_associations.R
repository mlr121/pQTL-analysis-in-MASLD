compute_ubr1_associations <- function(all_pheno, ubr1_col = "Ubr1") {
  
  # Requires mgcv
  if (!"mgcv" %in% .packages(all.available = TRUE)) {
    suppressPackageStartupMessages(require(mgcv))
  } else {
    suppressPackageStartupMessages(library(mgcv))
  }
  
  # Basic checks
  if (!ubr1_col %in% colnames(all_pheno)) {
    stop("ubr1_col not found in all_pheno: ", ubr1_col)
  }
  if (!"Sex" %in% colnames(all_pheno)) {
    stop("all_pheno must contain a 'Sex' column.")
  }
  
  sex_factor <- as.factor(all_pheno$Sex)
  
  # Identify numeric columns
  numeric_cols <- colnames(all_pheno)[sapply(all_pheno, is.numeric)]
  
  # Outcomes = all numeric proteins except Ubr1
  protein_cols <- setdiff(numeric_cols, ubr1_col)
  
  if (length(protein_cols) == 0) {
    stop("No protein outcome columns detected.")
  }
  
  ubr1_vals <- all_pheno[[ubr1_col]]
  
  message("Computing associations for ", ubr1_col)
  
  results <- lapply(protein_cols, function(prot) {
    
    y <- all_pheno[[prot]]
    
    # Require >= 3 paired non-missing observations
    if (sum(!is.na(ubr1_vals) & !is.na(y)) < 3) return(NULL)
    
    df_model <- data.frame(
      outcome = y,
      sex = sex_factor,
      ubr1 = ubr1_vals
    )
    
    fit <- tryCatch({
      mgcv::gam(outcome ~ s(ubr1) + sex, data = df_model)
    }, error = function(e) return(NULL))
    
    if (is.null(fit)) return(NULL)
    
    s_summary <- summary(fit)
    
    pval <- tryCatch({
      s_summary$s.table["s(ubr1)", "p-value"]
    }, error = function(e) NA_real_)
    
    dev_expl <- tryCatch({
      s_summary$dev.expl
    }, error = function(e) NA_real_)
    
    tau_val <- tryCatch({
      cor(df_model$ubr1, df_model$outcome,
          method = "kendall", use = "complete.obs")
    }, error = function(e) NA_real_)
    
    direction <- if (is.na(tau_val)) NA_character_
    else if (tau_val > 0) "1"
    else if (tau_val < 0) "-1"
    else "neutral"
    
    data.frame(
      Protein = prot,
      p_value = pval,
      deviance_explained = dev_expl,
      tau = tau_val,
      direction = direction,
      stringsAsFactors = FALSE
    )
  })
  
  # bind + drop NULLs
  results_df <- do.call(rbind, results)
  
  return(results_df)
}